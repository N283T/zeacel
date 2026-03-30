// Stockholm format parser and writer.
//
// Format specification:
//   - First line: "# STOCKHOLM 1.0"
//   - Last line: "//"
//   - #=GF <tag> <value>  : per-file annotation
//   - #=GS <name> <tag> <value> : per-sequence annotation (skipped)
//   - #=GC <tag> <annotation>   : per-column annotation (skipped)
//   - #=GR <name> <tag> <annotation> : per-residue annotation (skipped)
//   - <name>  <sequence> : aligned sequence data
//   - Interleaved (multi-block) alignments are supported.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const Msa = @import("../msa.zig").Msa;

/// Parse a Stockholm format alignment from a byte buffer.
/// #=GF ID, AC, DE are stored in msa.name, msa.accession, msa.description.
/// All other markup lines are ignored.
pub fn parse(allocator: Allocator, abc: *const Alphabet, data: []const u8) !Msa {
    var lines = std.mem.splitScalar(u8, data, '\n');

    // Verify first non-blank line is the Stockholm header.
    var found_header = false;
    while (lines.next()) |line| {
        const trimmed = std.mem.trimRight(u8, line, "\r");
        if (trimmed.len == 0) continue;
        if (!std.mem.eql(u8, trimmed, "# STOCKHOLM 1.0")) return error.InvalidFormat;
        found_header = true;
        break;
    }
    if (!found_header) return error.InvalidFormat;

    // Metadata fields.
    var msa_name: ?[]const u8 = null;
    var msa_accession: ?[]const u8 = null;
    var msa_description: ?[]const u8 = null;
    errdefer {
        if (msa_name) |n| allocator.free(n);
        if (msa_accession) |a| allocator.free(a);
        if (msa_description) |d| allocator.free(d);
    }

    // Ordered list of names (first-seen order).
    var name_list = std.ArrayList([]const u8){};
    defer {
        for (name_list.items) |n| allocator.free(n);
        name_list.deinit(allocator);
    }

    // name_index maps sequence name slice (owned by name_list) to its index.
    var name_index = std.StringHashMap(usize).init(allocator);
    defer name_index.deinit();

    // seq_parts[i] collects all text fragments for sequence i across blocks.
    var seq_parts = std.ArrayList(std.ArrayList(u8)){};
    defer {
        for (seq_parts.items) |*parts| parts.deinit(allocator);
        seq_parts.deinit(allocator);
    }

    var ended = false;

    while (lines.next()) |raw_line| {
        const line = std.mem.trimRight(u8, raw_line, "\r");

        // End of alignment record.
        if (std.mem.eql(u8, line, "//")) {
            ended = true;
            break;
        }

        // Blank line — block separator.
        if (line.len == 0) continue;

        // Markup lines.
        if (line[0] == '#') {
            // #=GF tag value
            if (std.mem.startsWith(u8, line, "#=GF ") or std.mem.startsWith(u8, line, "#=GF\t")) {
                const rest = trimLeft(line[5..]);
                // rest is "TAG value..."
                const tag_end = indexOfWhitespace(rest) orelse rest.len;
                const tag = rest[0..tag_end];
                const value = if (tag_end < rest.len) trimLeft(rest[tag_end..]) else "";

                if (std.mem.eql(u8, tag, "ID")) {
                    if (msa_name) |old| allocator.free(old);
                    msa_name = try allocator.dupe(u8, value);
                } else if (std.mem.eql(u8, tag, "AC")) {
                    if (msa_accession) |old| allocator.free(old);
                    msa_accession = try allocator.dupe(u8, value);
                } else if (std.mem.eql(u8, tag, "DE")) {
                    if (msa_description) |old| allocator.free(old);
                    msa_description = try allocator.dupe(u8, value);
                }
                // All other #=GF tags are ignored.
            }
            // #=GS, #=GC, #=GR — skip.
            continue;
        }

        // Sequence line: "<name> <sequence>"
        const name_end = indexOfWhitespace(line) orelse return error.InvalidFormat;
        const seq_name = line[0..name_end];
        const seq_text = trimLeft(line[name_end..]);

        if (seq_text.len == 0) return error.InvalidFormat;

        if (name_index.get(seq_name)) |idx| {
            // Append to existing sequence.
            try seq_parts.items[idx].appendSlice(allocator, seq_text);
        } else {
            // New sequence.
            const idx = name_list.items.len;
            const owned_name = try allocator.dupe(u8, seq_name);
            errdefer allocator.free(owned_name);

            try name_list.append(allocator, owned_name);
            errdefer _ = name_list.pop();

            try name_index.put(owned_name, idx);

            var parts = std.ArrayList(u8){};
            errdefer parts.deinit(allocator);

            try parts.appendSlice(allocator, seq_text);
            try seq_parts.append(allocator, parts);
        }
    }

    if (!ended) return error.InvalidFormat;
    if (name_list.items.len == 0) return error.InvalidInput;

    // Build slices for Msa.fromText.
    const n = name_list.items.len;

    const names_slice = try allocator.alloc([]const u8, n);
    defer allocator.free(names_slice);
    for (name_list.items, 0..) |name, i| {
        names_slice[i] = name;
    }

    const text_seqs = try allocator.alloc([]const u8, n);
    defer allocator.free(text_seqs);
    for (seq_parts.items, 0..) |*parts, i| {
        text_seqs[i] = parts.items;
    }

    var msa = try Msa.fromText(allocator, abc, names_slice, text_seqs);
    errdefer msa.deinit();

    // Transfer ownership of metadata strings to msa.
    msa.name = msa_name;
    msa_name = null;
    msa.accession = msa_accession;
    msa_accession = null;
    msa.description = msa_description;
    msa_description = null;

    return msa;
}

/// Write an Msa in Stockholm format.
pub fn write(dest: std.io.AnyWriter, m: Msa) !void {
    try dest.writeAll("# STOCKHOLM 1.0\n");

    // Write per-file annotations if present.
    if (m.name) |name| {
        try dest.print("#=GF ID   {s}\n", .{name});
    }
    if (m.accession) |acc| {
        try dest.print("#=GF AC   {s}\n", .{acc});
    }
    if (m.description) |desc| {
        try dest.print("#=GF DE   {s}\n", .{desc});
    }

    try dest.writeByte('\n');

    // Compute column width for name padding.
    var max_name_len: usize = 0;
    for (m.names) |name| {
        if (name.len > max_name_len) max_name_len = name.len;
    }
    // Pad names to max_name_len + 2 spaces minimum separation.
    const col_width = max_name_len + 2;

    for (0..m.nseq()) |i| {
        const name = m.names[i];
        const text = try m.abc.textize(m.allocator, m.seqs[i]);
        defer m.allocator.free(text);

        try dest.writeAll(name);
        // Pad with spaces.
        const padding = col_width - name.len;
        var p: usize = 0;
        while (p < padding) : (p += 1) {
            try dest.writeByte(' ');
        }
        try dest.writeAll(text);
        try dest.writeByte('\n');
    }

    try dest.writeAll("//\n");
}

// --- Helpers ---

fn trimLeft(s: []const u8) []const u8 {
    var i: usize = 0;
    while (i < s.len and (s[i] == ' ' or s[i] == '\t')) i += 1;
    return s[i..];
}

fn indexOfWhitespace(s: []const u8) ?usize {
    for (s, 0..) |c, i| {
        if (c == ' ' or c == '\t') return i;
    }
    return null;
}

// --- Tests ---

const alphabet_mod = @import("../alphabet.zig");
const testing = std.testing;

test "parse: simple alignment" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  ACGA
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expectEqual(@as(usize, 2), msa.nseq());
    try testing.expectEqual(@as(usize, 4), msa.alen);
    try testing.expectEqualStrings("seq1", msa.names[0]);
    try testing.expectEqualStrings("seq2", msa.names[1]);

    const text0 = try abc.textize(allocator, msa.seqs[0]);
    defer allocator.free(text0);
    try testing.expectEqualStrings("ACGT", text0);

    const text1 = try abc.textize(allocator, msa.seqs[1]);
    defer allocator.free(text1);
    try testing.expectEqualStrings("ACGA", text1);
}

test "parse: with GF metadata" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\#=GF ID   myalignment
        \\#=GF AC   PF00001
        \\#=GF DE   Test alignment
        \\
        \\seq1  ACGT
        \\seq2  ACGA
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expectEqual(@as(usize, 2), msa.nseq());
    try testing.expectEqualStrings("myalignment", msa.name.?);
    try testing.expectEqualStrings("PF00001", msa.accession.?);
    try testing.expectEqualStrings("Test alignment", msa.description.?);
}

test "parse: interleaved multi-block alignment" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  ACGA
        \\
        \\seq1  TTTT
        \\seq2  CCCC
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expectEqual(@as(usize, 2), msa.nseq());
    try testing.expectEqual(@as(usize, 8), msa.alen);

    const text0 = try abc.textize(allocator, msa.seqs[0]);
    defer allocator.free(text0);
    try testing.expectEqualStrings("ACGTTTTT", text0);

    const text1 = try abc.textize(allocator, msa.seqs[1]);
    defer allocator.free(text1);
    try testing.expectEqualStrings("ACGACCCC", text1);
}

test "parse: skips GS, GC, GR markup" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data =
        \\# STOCKHOLM 1.0
        \\
        \\#=GS seq1 WT 0.42
        \\#=GC SS_cons ....
        \\seq1  ACGT
        \\#=GR seq1 SS ....
        \\seq2  TTTT
        \\//
    ;

    var msa = try parse(allocator, abc, data);
    defer msa.deinit();

    try testing.expectEqual(@as(usize, 2), msa.nseq());
    try testing.expectEqual(@as(usize, 4), msa.alen);
}

test "parse: missing header returns error" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data = "seq1  ACGT\n//\n";
    try testing.expectError(error.InvalidFormat, parse(allocator, abc, data));
}

test "parse: missing terminator returns error" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const data = "# STOCKHOLM 1.0\n\nseq1  ACGT\n";
    try testing.expectError(error.InvalidFormat, parse(allocator, abc, data));
}

test "write: basic output" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "seq1", "seq2" };
    const seqs = [_][]const u8{ "ACGT", "TTTT" };

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa);

    const expected =
        \\# STOCKHOLM 1.0
        \\
        \\seq1  ACGT
        \\seq2  TTTT
        \\//
        \\
    ;
    try testing.expectEqualStrings(expected, buf.items);
}

test "write: with metadata" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{"seq1"};
    const seqs = [_][]const u8{"ACGT"};

    var msa = try Msa.fromText(allocator, abc, &names, &seqs);
    defer msa.deinit();
    msa.name = try allocator.dupe(u8, "myaln");
    msa.accession = try allocator.dupe(u8, "PF00001");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), msa);

    try testing.expect(std.mem.indexOf(u8, buf.items, "#=GF ID   myaln") != null);
    try testing.expect(std.mem.indexOf(u8, buf.items, "#=GF AC   PF00001") != null);
}

test "round-trip: write then parse" {
    const allocator = testing.allocator;
    const abc = &alphabet_mod.dna;

    const names = [_][]const u8{ "alpha", "beta" };
    const seqs = [_][]const u8{ "ACGT-A", "ACGA-T" };

    var original = try Msa.fromText(allocator, abc, &names, &seqs);
    defer original.deinit();
    original.name = try allocator.dupe(u8, "testaln");

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try write(buf.writer(allocator).any(), original);

    var restored = try parse(allocator, abc, buf.items);
    defer restored.deinit();

    try testing.expectEqual(original.nseq(), restored.nseq());
    try testing.expectEqual(original.alen, restored.alen);
    try testing.expectEqualStrings(original.names[0], restored.names[0]);
    try testing.expectEqualStrings(original.names[1], restored.names[1]);
    try testing.expectEqualSlices(u8, original.seqs[0], restored.seqs[0]);
    try testing.expectEqualSlices(u8, original.seqs[1], restored.seqs[1]);
    try testing.expectEqualStrings(original.name.?, restored.name.?);
}
