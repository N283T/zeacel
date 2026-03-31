// SSI (Simple Sequence Index) — fast random access to sequences by name.
//
// Binary format:
//   Header:
//     magic:       [4]u8  = "ZSSI"
//     version:     u32    = 1   (little-endian)
//     num_entries: u64          (little-endian)
//   Per entry:
//     name_len:    u32          (little-endian)
//     name:        [name_len]u8
//     offset:      u64          (little-endian) — byte offset of '>' header
//     data_offset: u64          (little-endian) — byte offset of first residue
//     seq_len:     u64          (little-endian) — number of residues (excl. gaps/whitespace)

const std = @import("std");
const Allocator = std.mem.Allocator;

const MAGIC = "ZSSI";
const VERSION: u32 = 1;

pub const SsiEntry = struct {
    name: []const u8,
    offset: u64,
    data_offset: u64,
    seq_len: u64,
};

pub const SsiIndex = struct {
    entries: []SsiEntry,
    name_map: std.StringHashMap(usize),
    allocator: Allocator,

    /// Build an index from a FASTA byte buffer.
    /// Scans the data recording byte offsets of each '>' header.
    pub fn buildFromFasta(allocator: Allocator, data: []const u8) !SsiIndex {
        var entry_list = std.ArrayList(SsiEntry){};
        errdefer {
            for (entry_list.items) |e| allocator.free(e.name);
            entry_list.deinit(allocator);
        }

        var pos: usize = 0;
        while (pos < data.len) {
            // Skip whitespace / blank lines between records.
            while (pos < data.len) {
                const c = data[pos];
                if (c == '\n' or c == '\r' or c == ' ' or c == '\t') {
                    pos += 1;
                } else {
                    break;
                }
            }
            if (pos >= data.len) break;
            if (data[pos] != '>') return error.InvalidFormat;

            const record_offset = pos;

            // Advance past '>'.
            pos += 1;

            // Parse name (up to first whitespace or newline).
            const name_start = pos;
            while (pos < data.len and data[pos] != ' ' and data[pos] != '\t' and
                data[pos] != '\n' and data[pos] != '\r')
            {
                pos += 1;
            }
            const name = try allocator.dupe(u8, data[name_start..pos]);
            errdefer allocator.free(name);

            // Skip rest of header line.
            while (pos < data.len and data[pos] != '\n' and data[pos] != '\r') {
                pos += 1;
            }
            if (pos < data.len and data[pos] == '\r') pos += 1;
            if (pos < data.len and data[pos] == '\n') pos += 1;

            const data_start = pos;

            // Count residues and advance to next '>'.
            var seq_len: u64 = 0;
            while (pos < data.len and data[pos] != '>') {
                const c = data[pos];
                if (c != '\n' and c != '\r' and c != ' ' and c != '\t') {
                    seq_len += 1;
                }
                pos += 1;
            }

            try entry_list.append(allocator, SsiEntry{
                .name = name,
                .offset = record_offset,
                .data_offset = data_start,
                .seq_len = seq_len,
            });
        }

        const entries = try entry_list.toOwnedSlice(allocator);
        errdefer {
            for (entries) |e| allocator.free(e.name);
            allocator.free(entries);
        }

        var name_map = std.StringHashMap(usize).init(allocator);
        errdefer name_map.deinit();

        for (entries, 0..) |e, i| {
            try name_map.put(e.name, i);
        }

        return SsiIndex{
            .entries = entries,
            .name_map = name_map,
            .allocator = allocator,
        };
    }

    /// Look up a sequence by name. Returns the entry or null if not found.
    pub fn lookup(self: SsiIndex, name: []const u8) ?SsiEntry {
        const idx = self.name_map.get(name) orelse return null;
        return self.entries[idx];
    }

    /// Write the index to a binary stream (little-endian).
    pub fn write(self: SsiIndex, dest: std.io.AnyWriter) !void {
        // Magic + version + num_entries
        try dest.writeAll(MAGIC);
        try dest.writeInt(u32, VERSION, .little);
        try dest.writeInt(u64, @intCast(self.entries.len), .little);

        for (self.entries) |e| {
            try dest.writeInt(u32, @intCast(e.name.len), .little);
            try dest.writeAll(e.name);
            try dest.writeInt(u64, e.offset, .little);
            try dest.writeInt(u64, e.data_offset, .little);
            try dest.writeInt(u64, e.seq_len, .little);
        }
    }

    /// Read an index from a binary stream.
    /// Caller owns the returned SsiIndex and must call deinit().
    pub fn read(allocator: Allocator, src: std.io.AnyReader) !SsiIndex {
        var magic: [4]u8 = undefined;
        _ = try src.readAll(&magic);
        if (!std.mem.eql(u8, &magic, MAGIC)) return error.InvalidFormat;

        const version = try src.readInt(u32, .little);
        if (version != VERSION) return error.UnsupportedVersion;

        const num_entries = try src.readInt(u64, .little);

        const entries = try allocator.alloc(SsiEntry, num_entries);
        var entries_done: usize = 0;
        errdefer {
            for (0..entries_done) |i| allocator.free(entries[i].name);
            allocator.free(entries);
        }

        for (0..num_entries) |i| {
            const name_len = try src.readInt(u32, .little);
            const name = try allocator.alloc(u8, name_len);
            errdefer allocator.free(name);
            _ = try src.readAll(name);

            const offset = try src.readInt(u64, .little);
            const data_offset = try src.readInt(u64, .little);
            const seq_len = try src.readInt(u64, .little);

            entries[i] = SsiEntry{
                .name = name,
                .offset = offset,
                .data_offset = data_offset,
                .seq_len = seq_len,
            };
            entries_done += 1;
        }

        var name_map = std.StringHashMap(usize).init(allocator);
        errdefer name_map.deinit();

        for (entries, 0..) |e, i| {
            try name_map.put(e.name, i);
        }

        return SsiIndex{
            .entries = entries,
            .name_map = name_map,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *SsiIndex) void {
        for (self.entries) |e| self.allocator.free(e.name);
        self.allocator.free(self.entries);
        self.name_map.deinit();
    }
};

// --- Tests ---

test "buildFromFasta: correct offsets" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectEqual(@as(usize, 2), idx.entries.len);
    try std.testing.expectEqualStrings("seq1", idx.entries[0].name);
    try std.testing.expectEqualStrings("seq2", idx.entries[1].name);

    // ">seq1" starts at byte 0.
    try std.testing.expectEqual(@as(u64, 0), idx.entries[0].offset);
    // Data starts after ">seq1\n" = 6 bytes.
    try std.testing.expectEqual(@as(u64, 6), idx.entries[0].data_offset);
    // seq_len is 4 (ACGT).
    try std.testing.expectEqual(@as(u64, 4), idx.entries[0].seq_len);

    // ">seq2" starts at byte 11 (">seq1\nACGT\n" = 11 bytes).
    try std.testing.expectEqual(@as(u64, 11), idx.entries[1].offset);
    try std.testing.expectEqual(@as(u64, 4), idx.entries[1].seq_len);
}

test "lookup: finds entry by name" {
    const allocator = std.testing.allocator;
    const data = ">alpha\nACGT\n>beta\nGGGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    const e = idx.lookup("alpha") orelse return error.TestEntryNotFound;
    try std.testing.expectEqualStrings("alpha", e.name);
    try std.testing.expectEqual(@as(u64, 4), e.seq_len);
}

test "lookup: missing name returns null" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectEqual(@as(?SsiEntry, null), idx.lookup("missing"));
}

test "write and read: round trip" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\n>seq2\nGGGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try idx.write(buf.writer(allocator).any());

    var stream = std.io.fixedBufferStream(buf.items);
    var idx2 = try SsiIndex.read(allocator, stream.reader().any());
    defer idx2.deinit();

    try std.testing.expectEqual(idx.entries.len, idx2.entries.len);
    for (0..idx.entries.len) |i| {
        try std.testing.expectEqualStrings(idx.entries[i].name, idx2.entries[i].name);
        try std.testing.expectEqual(idx.entries[i].offset, idx2.entries[i].offset);
        try std.testing.expectEqual(idx.entries[i].data_offset, idx2.entries[i].data_offset);
        try std.testing.expectEqual(idx.entries[i].seq_len, idx2.entries[i].seq_len);
    }

    // lookup works after round trip.
    const e = idx2.lookup("seq1") orelse return error.TestEntryNotFound;
    try std.testing.expectEqual(@as(u64, 0), e.offset);
    try std.testing.expectEqual(@as(u64, 4), e.seq_len);
}

test "buildFromFasta: multi-line sequences" {
    const allocator = std.testing.allocator;
    const data = ">seq1\nACGT\nACGT\n>seq2\nGG\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectEqual(@as(u64, 8), idx.entries[0].seq_len);
    try std.testing.expectEqual(@as(u64, 2), idx.entries[1].seq_len);
}

test "buildFromFasta: sequence with description" {
    const allocator = std.testing.allocator;
    const data = ">seq1 my description\nACGT\n";

    var idx = try SsiIndex.buildFromFasta(allocator, data);
    defer idx.deinit();

    try std.testing.expectEqualStrings("seq1", idx.entries[0].name);
    try std.testing.expectEqual(@as(u64, 4), idx.entries[0].seq_len);
}
