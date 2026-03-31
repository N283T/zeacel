// Aligned FASTA (AFA) format parser and writer.
//
// AFA is syntactically identical to FASTA but all sequences must have the same
// length, including gap characters ('-'). Gaps are preserved as alignment
// columns when building an Msa.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const fasta = @import("fasta.zig");
const Msa = @import("../msa.zig").Msa;

/// Parse aligned FASTA into an Msa.
/// Unlike regular FASTA, gaps are preserved as alignment columns.
/// Returns error.InvalidInput if sequences have different lengths or if there
/// are fewer than 1 sequence.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    // Use the existing FASTA parser to read all records.  Gap characters are
    // valid in the alphabet (digitize handles '-'), so they are kept as-is.
    const seqs = try fasta.parseAll(allocator, abc, data);
    defer {
        for (seqs) |*s| @constCast(s).deinit();
        allocator.free(seqs);
    }

    if (seqs.len == 0) return error.InvalidInput;

    const alen = seqs[0].dsq.len;
    for (seqs) |s| {
        if (s.dsq.len != alen) return error.InvalidInput;
    }

    // Build name and text_seq slices for Msa.fromText.
    // We pass the already-digitized sequences through textize so fromText can
    // re-digitize them — keeping a consistent construction path.
    const n = seqs.len;

    var names = try allocator.alloc([]const u8, n);
    defer allocator.free(names);

    var text_seqs = try allocator.alloc([]u8, n);
    defer {
        for (text_seqs) |t| allocator.free(t);
        allocator.free(text_seqs);
    }

    for (0..n) |i| {
        names[i] = seqs[i].name;
        text_seqs[i] = try abc.textize(allocator, seqs[i].dsq);
    }

    return Msa.fromText(allocator, abc, names, @ptrCast(text_seqs));
}

/// Write an Msa as aligned FASTA.
/// Each sequence is written with gap characters preserved.
/// line_width controls residues per line (0 = no wrapping).
pub fn write(dest: std.io.AnyWriter, m: Msa, line_width: usize) !void {
    for (0..m.nseq()) |i| {
        try dest.writeByte('>');
        try dest.writeAll(m.names[i]);
        try dest.writeByte('\n');

        const text = try m.abc.textize(m.allocator, m.seqs[i]);
        defer m.allocator.free(text);

        if (line_width == 0 or text.len == 0) {
            try dest.writeAll(text);
            if (text.len > 0) try dest.writeByte('\n');
            continue;
        }

        var col: usize = 0;
        while (col < text.len) {
            const end = @min(col + line_width, text.len);
            try dest.writeAll(text[col..end]);
            try dest.writeByte('\n');
            col = end;
        }
    }
}

// --- Tests ---

const alphabet_mod = @import("../alphabet.zig");

test "parse: aligned FASTA with gaps produces correct alen" {
    const allocator = std.testing.allocator;
    const data =
        \\>seq1
        \\AC-GT
        \\>seq2
        \\A-CGT
        \\
    ;

    var msa = try parse(allocator, &alphabet_mod.dna, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 5), msa.alen);
    try std.testing.expectEqualStrings("seq1", msa.names[0]);

    // Gap code for DNA is 4.
    try std.testing.expectEqual(@as(u8, 4), msa.seqs[0][2]); // '-' at position 2
}

test "parse: unequal lengths return error" {
    const allocator = std.testing.allocator;
    const data =
        \\>seq1
        \\ACGT
        \\>seq2
        \\ACG
        \\
    ;

    try std.testing.expectError(
        error.InvalidInput,
        parse(allocator, &alphabet_mod.dna, data),
    );
}

test "parse: empty input returns error" {
    const allocator = std.testing.allocator;
    try std.testing.expectError(
        error.InvalidInput,
        parse(allocator, &alphabet_mod.dna, ""),
    );
}

test "write: preserves gaps" {
    const allocator = std.testing.allocator;
    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "AC-GT", "A-CGT" };

    var msa = try Msa.fromText(allocator, &alphabet_mod.dna, &names, &seqs);
    defer msa.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa, 60);

    try std.testing.expect(std.mem.indexOf(u8, buf.items, ">seq1") != null);
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "AC-GT") != null);
    try std.testing.expect(std.mem.indexOf(u8, buf.items, "A-CGT") != null);
}

test "round trip" {
    const allocator = std.testing.allocator;
    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACDE-FGHI", "ACD--FGHI" };

    var msa = try Msa.fromText(allocator, &alphabet_mod.amino, &names, &seqs);
    defer msa.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa, 60);

    var msa2 = try parse(allocator, &alphabet_mod.amino, buf.items);
    defer msa2.deinit();

    try std.testing.expectEqual(msa.nseq(), msa2.nseq());
    try std.testing.expectEqual(msa.alen, msa2.alen);
    for (0..msa.nseq()) |i| {
        try std.testing.expectEqualStrings(msa.names[i], msa2.names[i]);
        try std.testing.expectEqualSlices(u8, msa.seqs[i], msa2.seqs[i]);
    }
}
