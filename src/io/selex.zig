// SELEX MSA format parser/writer.
//
// SELEX is a block-based format where each line is: name  sequence
// Spaces in sequence data are treated as gaps (mapped to '-').
// Annotation lines: #=RF, #=CS, #=SS, #=SA can appear in blocks.
// Blocks are separated by blank lines.
//
// Per-sequence annotation semantics (matching Easel esl_msafile_selex.c):
//   #=SS and #=SA follow the sequence they annotate (the most recently seen
//   sequence line). They are stored as per-sequence annotations in msa.ss /
//   msa.sa, indexed by sequence index.
//
//   #=RF is per-alignment (consensus), stored in msa.reference.
//   #=CS is per-alignment (consensus secondary structure), stored in msa.consensus_ss.
//   #=SA with no preceding sequence is an error in strict mode; here it is
//   silently ignored to be lenient with malformed files.
//
// Note on name-length parsing: sequence data is extracted by trimming the
// name token and then treating the remainder as the sequence field. This
// works correctly when names contain no embedded spaces. Files with
// inconsistent name-column widths across blocks may be misparsed if names
// have trailing spaces; such files are non-standard.

const std = @import("std");
const Allocator = std.mem.Allocator;
const msa_mod = @import("../msa.zig");
const Msa = msa_mod.Msa;
const GcEntry = msa_mod.GcEntry;
const Alphabet = @import("../alphabet.zig").Alphabet;

/// Parse a SELEX format alignment.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    var name_list: std.ArrayList([]const u8) = .empty;
    defer {
        for (name_list.items) |n| allocator.free(n);
        name_list.deinit(allocator);
    }
    var seq_bufs: std.ArrayList(std.ArrayList(u8)) = .empty;
    defer {
        for (seq_bufs.items) |*s| s.deinit(allocator);
        seq_bufs.deinit(allocator);
    }

    var name_to_idx = std.StringHashMap(usize).init(allocator);
    defer name_to_idx.deinit();

    // Buffers for per-alignment annotation lines (accumulated across blocks)
    var rf_buf: std.ArrayList(u8) = .empty;
    defer rf_buf.deinit(allocator);
    var cs_buf: std.ArrayList(u8) = .empty;
    defer cs_buf.deinit(allocator);

    // Per-sequence SS and SA buffers, indexed by sequence index.
    // Allocated lazily when the first #=SS or #=SA line is encountered.
    var ss_bufs: ?std.ArrayList(std.ArrayList(u8)) = null;
    defer if (ss_bufs) |*b| {
        for (b.items) |*s| s.deinit(allocator);
        b.deinit(allocator);
    };
    var sa_bufs: ?std.ArrayList(std.ArrayList(u8)) = null;
    defer if (sa_bufs) |*b| {
        for (b.items) |*s| s.deinit(allocator);
        b.deinit(allocator);
    };

    // Sentinel used to append a new empty per-sequence buffer slot.
    const empty_buf: std.ArrayList(u8) = .empty;

    // Index of the most recently parsed sequence line (null = no sequence seen yet).
    var last_seq_idx: ?usize = null;

    var lines = std.mem.splitScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        // Parse annotation lines starting with '#'
        if (trimmed[0] == '#') {
            // Parse known annotation tags: #=RF, #=CS, #=SS, #=SA
            var it = std.mem.tokenizeAny(u8, trimmed, " \t");
            const tag = it.next() orelse continue;
            const annotation = it.rest();
            if (annotation.len == 0) continue;

            if (std.mem.eql(u8, tag, "#=RF")) {
                try rf_buf.appendSlice(allocator, annotation);
            } else if (std.mem.eql(u8, tag, "#=CS")) {
                try cs_buf.appendSlice(allocator, annotation);
            } else if (std.mem.eql(u8, tag, "#=SS")) {
                // #=SS annotates the most recently seen sequence.
                // Silently ignore if no sequence has been seen yet.
                if (last_seq_idx) |si| {
                    // Lazily allocate per-sequence SS buffers.
                    if (ss_bufs == null) {
                        ss_bufs = .empty;
                        // Pre-populate with empty buffers for all sequences seen so far.
                        for (0..name_list.items.len) |_| {
                            try ss_bufs.?.append(allocator, empty_buf);
                        }
                    }
                    // Ensure there is a buffer for si (may have been added after init).
                    while (ss_bufs.?.items.len <= si) {
                        try ss_bufs.?.append(allocator, empty_buf);
                    }
                    try ss_bufs.?.items[si].appendSlice(allocator, annotation);
                }
            } else if (std.mem.eql(u8, tag, "#=SA")) {
                // #=SA annotates the most recently seen sequence.
                // Silently ignore if no sequence has been seen yet.
                if (last_seq_idx) |si| {
                    // Lazily allocate per-sequence SA buffers.
                    if (sa_bufs == null) {
                        sa_bufs = .empty;
                        for (0..name_list.items.len) |_| {
                            try sa_bufs.?.append(allocator, empty_buf);
                        }
                    }
                    while (sa_bufs.?.items.len <= si) {
                        try sa_bufs.?.append(allocator, empty_buf);
                    }
                    try sa_bufs.?.items[si].appendSlice(allocator, annotation);
                }
            }
            // Other comment lines are silently skipped
            continue;
        }

        // Find where name ends — first whitespace
        var name_end: usize = 0;
        while (name_end < trimmed.len and trimmed[name_end] != ' ' and trimmed[name_end] != '\t') : (name_end += 1) {}
        if (name_end == 0) continue;
        const name = trimmed[0..name_end];

        // Rest is sequence (spaces become gaps)
        const seq_part = if (name_end < trimmed.len) trimmed[name_end..] else "";

        // Extract sequence: skip leading whitespace, then map spaces to gaps
        const seq_trimmed = std.mem.trimLeft(u8, seq_part, " \t");

        if (name_to_idx.get(name)) |idx| {
            for (seq_trimmed) |ch| {
                if (ch == '\t') continue;
                const out: u8 = if (ch == ' ') '-' else ch;
                try seq_bufs.items[idx].append(allocator, out);
            }
            last_seq_idx = idx;
        } else {
            const name_copy = try allocator.dupe(u8, name);
            errdefer allocator.free(name_copy);
            const idx = name_list.items.len;
            try name_to_idx.put(name_copy, idx);
            try name_list.append(allocator, name_copy);

            var buf: std.ArrayList(u8) = .empty;
            for (seq_trimmed) |ch| {
                if (ch == '\t') continue;
                const out: u8 = if (ch == ' ') '-' else ch;
                try buf.append(allocator, out);
            }
            try seq_bufs.append(allocator, buf);

            // Extend ss_bufs / sa_bufs if they have already been allocated,
            // so that buffer index always matches sequence index.
            if (ss_bufs) |*b| {
                while (b.items.len <= idx) {
                    try b.append(allocator, empty_buf);
                }
            }
            if (sa_bufs) |*b| {
                while (b.items.len <= idx) {
                    try b.append(allocator, empty_buf);
                }
            }

            last_seq_idx = idx;
        }
    }

    if (name_list.items.len == 0) return error.InvalidFormat;

    const n = name_list.items.len;
    var text_seqs = try allocator.alloc([]const u8, n);
    defer allocator.free(text_seqs);
    for (0..n) |i| {
        text_seqs[i] = seq_bufs.items[i].items;
    }

    var result = try Msa.fromText(allocator, abc, name_list.items, text_seqs);
    errdefer result.deinit();

    for (name_list.items) |nm| allocator.free(nm);
    name_list.items.len = 0;

    // Store parsed per-alignment annotations
    if (rf_buf.items.len > 0) {
        result.reference = try allocator.dupe(u8, rf_buf.items);
    }
    if (cs_buf.items.len > 0) {
        result.consensus_ss = try allocator.dupe(u8, cs_buf.items);
    }

    // Store per-sequence SS annotations
    if (ss_bufs) |*b| {
        for (b.items, 0..) |*ss, i| {
            if (ss.items.len > 0) {
                try result.setSeqSS(i, ss.items);
            }
        }
    }

    // Store per-sequence SA annotations
    if (sa_bufs) |*b| {
        for (b.items, 0..) |*sa, i| {
            if (sa.items.len > 0) {
                try result.setSeqSA(i, sa.items);
            }
        }
    }

    return result;
}

/// Write an MSA in SELEX format (single block).
pub fn write(msa: Msa, dest: std.io.AnyWriter) !void {
    // Find max name length for alignment
    var max_name: usize = 0;
    for (msa.names) |name| {
        if (name.len > max_name) max_name = name.len;
    }

    for (0..msa.nseq()) |i| {
        try dest.writeAll(msa.names[i]);
        for (0..max_name - msa.names[i].len + 2) |_| try dest.writeByte(' ');
        for (msa.seqs[i]) |code| {
            try dest.writeByte(msa.abc.decode(code));
        }
        try dest.writeByte('\n');
    }
}

// --- Tests ---

test "parse: simple SELEX" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\seq1  ACGT
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 4), msa.alen);
}

test "parse: spaces as gaps" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\seq1  AC-GT
        \\seq2  TGCAA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 5), msa.alen);
}

test "parse: skips plain comment lines" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\# comment
        \\seq1  ACGT
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(?[]const u8, null), msa.reference);
}

test "parse: RF annotation stored" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\seq1  ACGT
        \\#=RF  xx.x
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    const rf = msa.reference orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("xx.x", rf);
}

test "parse: SS annotation stored as per-sequence" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    // #=SS follows seq1, so it should be stored as seq1's per-sequence SS.
    const data =
        \\seq1  ACGT
        \\#=SS  <<>>
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    // consensus_ss must NOT be set — #=SS is per-sequence, not consensus
    try std.testing.expectEqual(@as(?[]const u8, null), msa.consensus_ss);
    // seq1 (index 0) should have the SS annotation
    const ss = msa.ss orelse return error.TestUnexpectedResult;
    const seq1_ss = ss[0] orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("<<>>", seq1_ss);
    // seq2 (index 1) has no SS
    try std.testing.expectEqual(@as(?[]const u8, null), ss[1]);
}

test "parse: CS annotation stored as consensus_ss" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\seq1  ACGT
        \\#=CS  <<>>
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    const cs = msa.consensus_ss orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("<<>>", cs);
}

test "parse: SA annotation stored as per-sequence" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    // #=SA follows seq1, annotating it specifically.
    const data =
        \\seq1  ACGT
        \\#=SA  1234
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    // Should NOT appear in gc_markup
    try std.testing.expectEqual(@as(?[]GcEntry, null), msa.gc_markup);
    // seq1 (index 0) should have the SA annotation
    const sa = msa.sa orelse return error.TestUnexpectedResult;
    const seq1_sa = sa[0] orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("1234", seq1_sa);
    // seq2 (index 1) has no SA
    try std.testing.expectEqual(@as(?[]const u8, null), sa[1]);
}

test "parse: multi-sequence each with own SS" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    // Each sequence has its own #=SS line immediately following it.
    const data =
        \\seq1  ACGT
        \\#=SS  <<>>
        \\seq2  TGCA
        \\#=SS  ....
        \\seq3  AAAA
        \\#=SS  ::::
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 3), msa.nseq());
    try std.testing.expectEqual(@as(?[]const u8, null), msa.consensus_ss);

    const ss = msa.ss orelse return error.TestUnexpectedResult;
    const ss0 = ss[0] orelse return error.TestUnexpectedResult;
    const ss1 = ss[1] orelse return error.TestUnexpectedResult;
    const ss2 = ss[2] orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("<<>>", ss0);
    try std.testing.expectEqualStrings("....", ss1);
    try std.testing.expectEqualStrings("::::", ss2);
}

test "parse: multi-block per-sequence SS accumulated" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    // Two blocks; each block contributes SS annotations for seq1 only.
    const data = "seq1  AC\n#=SS  <<\nseq2  TG\n\nseq1  GT\n#=SS  >>\nseq2  CA\n";

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 4), msa.alen);

    const ss = msa.ss orelse return error.TestUnexpectedResult;
    const ss0 = ss[0] orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("<<>>", ss0);
    try std.testing.expectEqual(@as(?[]const u8, null), ss[1]);
}

test "parse: SS before any sequence is ignored" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    // #=SS at the top, before any sequence line, should be silently ignored.
    const data =
        \\#=SS  ????
        \\seq1  ACGT
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(?[]const u8, null), msa.consensus_ss);
    try std.testing.expectEqual(@as(?[]?[]const u8, null), msa.ss);
}

test "parse: CS and SA stored correctly (CS=consensus, SA=per-seq)" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data =
        \\seq1  ACGT
        \\#=CS  HHhh
        \\#=SA  1234
        \\seq2  TGCA
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    // #=CS → consensus_ss
    const cs = msa.consensus_ss orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("HHhh", cs);

    // #=SA → per-sequence SA for seq1 (index 0)
    const sa = msa.sa orelse return error.TestUnexpectedResult;
    const seq1_sa = sa[0] orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("1234", seq1_sa);
    try std.testing.expectEqual(@as(?[]const u8, null), sa[1]);
}

test "parse: multi-block RF annotation accumulated" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const data = "seq1  AC\n#=RF  xx\nseq2  TG\n\nseq1  GT\n#=RF  .x\nseq2  CA\n";

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try std.testing.expectEqual(@as(usize, 2), msa.nseq());
    try std.testing.expectEqual(@as(usize, 4), msa.alen);
    const rf = msa.reference orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("xx.x", rf);
}

test "write: round trip" {
    const allocator = std.testing.allocator;
    const abc = &@import("../alphabet.zig").dna;

    const names = [_][]const u8{ "s1", "s2" };
    const seqs = [_][]const u8{ "ACGT", "TGCA" };
    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var buf: std.ArrayList(u8) = .empty;
    defer buf.deinit(allocator);
    try write(msa, buf.writer(allocator).any());

    var msa2 = try parse(allocator, abc, buf.items);
    defer msa2.deinit();

    try std.testing.expect(msa.compare(msa2));
}
