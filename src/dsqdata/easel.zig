// Easel 4-file dsqdata format reader (read-only).
//
// Reads databases created by Easel's esl_dsqdata module, which stores
// sequences across three binary files (.dsqi, .dsqm, .dsqs) plus a
// stub file. All integers are native byte order (little-endian on x86/ARM).

const std = @import("std");
const Allocator = std.mem.Allocator;
const alphabet_mod = @import("../alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;
const Sequence = @import("../sequence.zig").Sequence;
const pack = @import("pack.zig");

// Easel dsqdata magic number (little-endian).
const MAGIC: u32 = 0xc4d3d1b1;
const MAGIC_SWAPPED: u32 = 0xb1d1d3c4;

pub const EaselDsqData = struct {
    allocator: Allocator,
    idx_file: std.fs.File, // .dsqi
    meta_file: std.fs.File, // .dsqm
    seq_file: std.fs.File, // .dsqs

    // From .dsqi header (56 bytes)
    uniquetag: u32,
    alphabet_type: alphabet_mod.AlphabetType,
    max_name_len: u32,
    max_acc_len: u32,
    max_desc_len: u32,
    max_seq_len: u64,
    num_sequences: u64,
    total_residues: u64,

    // Read cursor
    current_seq: u64,

    pub fn open(allocator: Allocator, basename: []const u8) !EaselDsqData {
        // Open .dsqi
        const dsqi_path = try std.fmt.allocPrint(allocator, "{s}.dsqi", .{basename});
        defer allocator.free(dsqi_path);
        const idx_file = try std.fs.cwd().openFile(dsqi_path, .{});
        errdefer idx_file.close();

        // Open .dsqm
        const dsqm_path = try std.fmt.allocPrint(allocator, "{s}.dsqm", .{basename});
        defer allocator.free(dsqm_path);
        const meta_file = try std.fs.cwd().openFile(dsqm_path, .{});
        errdefer meta_file.close();

        // Open .dsqs
        const dsqs_path = try std.fmt.allocPrint(allocator, "{s}.dsqs", .{basename});
        defer allocator.free(dsqs_path);
        const seq_file = try std.fs.cwd().openFile(dsqs_path, .{});
        errdefer seq_file.close();

        // Read .dsqi header (56 bytes) using positional read
        var idx_hdr: [56]u8 = undefined;
        _ = idx_file.preadAll(&idx_hdr, 0) catch return error.InvalidFormat;

        const magic = std.mem.readInt(u32, idx_hdr[0..4], .little);
        if (magic == MAGIC_SWAPPED) return error.ByteswapNotSupported;
        if (magic != MAGIC) return error.InvalidFormat;

        const uniquetag = std.mem.readInt(u32, idx_hdr[4..8], .little);
        const alphatype_raw = std.mem.readInt(u32, idx_hdr[8..12], .little);
        // idx_hdr[12..16] = flags (ignored)
        const max_name_len = std.mem.readInt(u32, idx_hdr[16..20], .little);
        const max_acc_len = std.mem.readInt(u32, idx_hdr[20..24], .little);
        const max_desc_len = std.mem.readInt(u32, idx_hdr[24..28], .little);
        // idx_hdr[28..32] = padding (ignored)
        const max_seq_len = std.mem.readInt(u64, idx_hdr[32..40], .little);
        const nseq = std.mem.readInt(u64, idx_hdr[40..48], .little);
        const nres = std.mem.readInt(u64, idx_hdr[48..56], .little);

        // Read .dsqm header (8 bytes: magic + uniquetag)
        var meta_hdr: [8]u8 = undefined;
        _ = meta_file.preadAll(&meta_hdr, 0) catch return error.InvalidFormat;

        const meta_magic = std.mem.readInt(u32, meta_hdr[0..4], .little);
        if (meta_magic == MAGIC_SWAPPED) return error.ByteswapNotSupported;
        if (meta_magic != MAGIC) return error.InvalidFormat;
        const meta_tag = std.mem.readInt(u32, meta_hdr[4..8], .little);
        if (meta_tag != uniquetag) return error.TagMismatch;

        // Read .dsqs header (8 bytes: magic + uniquetag)
        var seq_hdr: [8]u8 = undefined;
        _ = seq_file.preadAll(&seq_hdr, 0) catch return error.InvalidFormat;

        const seq_magic = std.mem.readInt(u32, seq_hdr[0..4], .little);
        if (seq_magic == MAGIC_SWAPPED) return error.ByteswapNotSupported;
        if (seq_magic != MAGIC) return error.InvalidFormat;
        const seq_tag = std.mem.readInt(u32, seq_hdr[4..8], .little);
        if (seq_tag != uniquetag) return error.TagMismatch;

        // Seek past the headers for sequential reading.
        // preadAll does not advance file position, so we seek manually.
        meta_file.seekTo(8) catch return error.InvalidFormat;
        seq_file.seekTo(8) catch return error.InvalidFormat;

        // Map Easel alphabet type: 1=dna, 2=rna, 3=amino
        const alphabet_type: alphabet_mod.AlphabetType = switch (alphatype_raw) {
            1 => .dna,
            2 => .rna,
            3 => .amino,
            else => return error.InvalidFormat,
        };

        return EaselDsqData{
            .allocator = allocator,
            .idx_file = idx_file,
            .meta_file = meta_file,
            .seq_file = seq_file,
            .uniquetag = uniquetag,
            .alphabet_type = alphabet_type,
            .max_name_len = max_name_len,
            .max_acc_len = max_acc_len,
            .max_desc_len = max_desc_len,
            .max_seq_len = max_seq_len,
            .num_sequences = nseq,
            .total_residues = nres,
            .current_seq = 0,
        };
    }

    pub fn readNext(self: *EaselDsqData, allocator: Allocator, abc: *const Alphabet) !?Sequence {
        if (self.current_seq >= self.num_sequences) return null;

        const meta_reader = self.meta_file.deprecatedReader();
        const seq_reader = self.seq_file.deprecatedReader();

        // Read metadata from .dsqm: name\0, accession\0, description\0, taxonomy_id(i32)
        const name = try readNullTerminated(allocator, meta_reader);
        errdefer allocator.free(name);

        const acc_raw = try readNullTerminated(allocator, meta_reader);
        errdefer allocator.free(acc_raw);

        const desc_raw = try readNullTerminated(allocator, meta_reader);
        errdefer allocator.free(desc_raw);

        const taxonomy_id = try meta_reader.readInt(i32, .little);

        // Read packed packets from .dsqs until EOD bit is set
        var packet_list: std.ArrayList(u32) = .empty;
        defer packet_list.deinit(allocator);

        while (true) {
            const pkt = try seq_reader.readInt(u32, .little);
            try packet_list.append(allocator, pkt);
            if ((pkt & pack.EOD) != 0) break;
        }

        const packets = packet_list.items;

        // Count residues from packets
        const is_amino = (self.alphabet_type == .amino);
        const seq_len = countResiduesFromPackets(packets, is_amino);

        // Unpack
        const dsq = if (is_amino)
            try pack.unpack5(allocator, packets, seq_len)
        else
            try pack.unpack2(allocator, packets, seq_len);
        errdefer allocator.free(dsq);

        // Convert empty strings to null
        const accession: ?[]const u8 = if (acc_raw.len == 0) blk: {
            allocator.free(acc_raw);
            break :blk null;
        } else acc_raw;

        const description: ?[]const u8 = if (desc_raw.len == 0) blk: {
            allocator.free(desc_raw);
            break :blk null;
        } else desc_raw;

        const tax_id: ?i32 = if (taxonomy_id == -1) null else taxonomy_id;

        self.current_seq += 1;

        return Sequence{
            .name = name,
            .accession = accession,
            .description = description,
            .taxonomy_id = tax_id,
            .dsq = dsq,
            .secondary_structure = null,
            .source = null,
            .abc = abc,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *EaselDsqData) void {
        self.idx_file.close();
        self.meta_file.close();
        self.seq_file.close();
    }
};

fn readNullTerminated(allocator: Allocator, reader: anytype) ![]u8 {
    var buf: std.ArrayList(u8) = .empty;
    errdefer buf.deinit(allocator);

    while (true) {
        const byte = try reader.readByte();
        if (byte == 0) break;
        try buf.append(allocator, byte);
    }

    return buf.toOwnedSlice(allocator);
}

fn countResiduesFromPackets(packets: []const u32, is_amino: bool) u64 {
    var count: u64 = 0;

    for (packets) |pkt| {
        const is_eod = (pkt & pack.EOD) != 0;
        const is_5bit = (pkt & pack.FIVEBIT) != 0;

        if (!is_eod) {
            if (is_5bit) {
                count += 6;
            } else {
                count += 15;
            }
        } else {
            if (is_5bit) {
                // Count non-sentinel positions
                for (0..6) |slot| {
                    const shift: u5 = @intCast((5 - slot) * 5);
                    const r: u8 = @intCast((pkt >> shift) & 0x1F);
                    if (r == 31) break;
                    count += 1;
                }
            } else {
                // 2-bit EOD is always 15 residues
                count += 15;
            }
        }
    }

    _ = is_amino;

    return count;
}

// --- Test helper ---

/// Write a minimal Easel database to a directory for testing.
/// Creates basename.dsqi, basename.dsqm, basename.dsqs, and basename stub.
pub fn writeTestEaselDb(
    dir: std.fs.Dir,
    basename: []const u8,
    uniquetag: u32,
    alphatype: u32,
    names: []const []const u8,
    packet_seqs: []const []const u32,
    seq_lens: []const u64,
) !void {
    _ = seq_lens;

    const nseq = names.len;

    // Compute total residues and max lengths
    var total_res: u64 = 0;
    var max_name_len: u32 = 0;
    for (0..nseq) |i| {
        for (packet_seqs[i]) |pkt| {
            const is_eod = (pkt & pack.EOD) != 0;
            const is_5bit = (pkt & pack.FIVEBIT) != 0;
            if (!is_eod) {
                total_res += if (is_5bit) 6 else 15;
            } else {
                if (is_5bit) {
                    for (0..6) |slot| {
                        const shift: u5 = @intCast((5 - slot) * 5);
                        const r: u8 = @intCast((pkt >> shift) & 0x1F);
                        if (r == 31) break;
                        total_res += 1;
                    }
                } else {
                    total_res += 15;
                }
            }
        }
        if (names[i].len > max_name_len) max_name_len = @intCast(names[i].len);
    }

    var path_buf: [256]u8 = undefined;

    // --- Write .dsqi ---
    {
        const path = try std.fmt.bufPrint(&path_buf, "{s}.dsqi", .{basename});
        const file = try dir.createFile(path, .{});
        defer file.close();
        const w = file.deprecatedWriter();

        // 56-byte header
        try w.writeInt(u32, MAGIC, .little);
        try w.writeInt(u32, uniquetag, .little);
        try w.writeInt(u32, alphatype, .little);
        try w.writeInt(u32, 0, .little); // flags
        try w.writeInt(u32, max_name_len, .little);
        try w.writeInt(u32, 0, .little); // max_acclen
        try w.writeInt(u32, 0, .little); // max_desclen
        try w.writeInt(u32, 0, .little); // padding
        try w.writeInt(u64, 0, .little); // max_seqlen
        try w.writeInt(u64, @intCast(nseq), .little);
        try w.writeInt(u64, total_res, .little);

        // Index records: metadata_end(i64) + psq_end(i64)
        var meta_offset: i64 = 8;
        var psq_offset: i64 = 8;

        for (0..nseq) |i| {
            const meta_size: i64 = @intCast(names[i].len + 1 + 1 + 1 + 4);
            meta_offset += meta_size;

            const psq_size: i64 = @intCast(packet_seqs[i].len * 4);
            psq_offset += psq_size;

            try w.writeInt(i64, meta_offset, .little);
            try w.writeInt(i64, psq_offset, .little);
        }
    }

    // --- Write .dsqm ---
    {
        const path = try std.fmt.bufPrint(&path_buf, "{s}.dsqm", .{basename});
        const file = try dir.createFile(path, .{});
        defer file.close();
        const w = file.deprecatedWriter();

        try w.writeInt(u32, MAGIC, .little);
        try w.writeInt(u32, uniquetag, .little);

        for (0..nseq) |i| {
            try w.writeAll(names[i]);
            try w.writeByte(0);
            try w.writeByte(0); // empty accession
            try w.writeByte(0); // empty description
            try w.writeInt(i32, -1, .little); // no taxonomy_id
        }
    }

    // --- Write .dsqs ---
    {
        const path = try std.fmt.bufPrint(&path_buf, "{s}.dsqs", .{basename});
        const file = try dir.createFile(path, .{});
        defer file.close();
        const w = file.deprecatedWriter();

        try w.writeInt(u32, MAGIC, .little);
        try w.writeInt(u32, uniquetag, .little);

        for (0..nseq) |i| {
            for (packet_seqs[i]) |pkt| {
                try w.writeInt(u32, pkt, .little);
            }
        }
    }

    // --- Write stub file ---
    {
        const file = try dir.createFile(basename, .{});
        defer file.close();
    }
}

// --- Tests ---

test "EaselDsqData.open: reads header from 3-file database" {
    const allocator = std.testing.allocator;

    // Create 1-seq amino db with 5 amino acids: A=0, C=1, D=2, E=3, F=4
    const pkt: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 0) << 25) |
        (@as(u32, 1) << 20) |
        (@as(u32, 2) << 15) |
        (@as(u32, 3) << 10) |
        (@as(u32, 4) << 5) |
        (@as(u32, 31) << 0);

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const packets = [_]u32{pkt};
    const names = [_][]const u8{"testseq"};
    const pkt_seqs = [_][]const u32{&packets};
    const lens = [_]u64{5};

    try writeTestEaselDb(tmp.dir, "test", 0x12345678, 3, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test", .{dir_path});
    defer allocator.free(full_path);

    var db = try EaselDsqData.open(allocator, full_path);
    defer db.deinit();

    try std.testing.expectEqual(@as(u32, 0x12345678), db.uniquetag);
    try std.testing.expectEqual(@as(u64, 1), db.num_sequences);
    try std.testing.expectEqual(alphabet_mod.AlphabetType.amino, db.alphabet_type);
}

test "EaselDsqData.readNext: reads amino acid sequence" {
    const allocator = std.testing.allocator;

    const pkt: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 0) << 25) |
        (@as(u32, 1) << 20) |
        (@as(u32, 2) << 15) |
        (@as(u32, 3) << 10) |
        (@as(u32, 4) << 5) |
        (@as(u32, 31) << 0);

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const packets = [_]u32{pkt};
    const names = [_][]const u8{"myprotein"};
    const pkt_seqs = [_][]const u32{&packets};
    const lens = [_]u64{5};

    try writeTestEaselDb(tmp.dir, "test", 0xAABBCCDD, 3, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/test", .{dir_path});
    defer allocator.free(full_path);

    var db = try EaselDsqData.open(allocator, full_path);
    defer db.deinit();

    var seq = (try db.readNext(allocator, &alphabet_mod.amino)).?;
    defer seq.deinit();

    try std.testing.expectEqualStrings("myprotein", seq.name);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3, 4 }, seq.dsq);
    try std.testing.expect(seq.accession == null);
    try std.testing.expect(seq.taxonomy_id == null);

    // No more sequences
    const next = try db.readNext(allocator, &alphabet_mod.amino);
    try std.testing.expect(next == null);
}

test "EaselDsqData.readNext: reads multiple DNA sequences" {
    const allocator = std.testing.allocator;

    // Sequence 1: 4 DNA residues via 5-bit encoding
    const pkt1: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 0) << 25) |
        (@as(u32, 1) << 20) |
        (@as(u32, 2) << 15) |
        (@as(u32, 3) << 10) |
        (@as(u32, 31) << 5) |
        (@as(u32, 31) << 0);

    // Sequence 2: 15 A's (code 0) via 2-bit EOD packet
    const pkt2: u32 = pack.EOD | 0;

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const packets1 = [_]u32{pkt1};
    const packets2 = [_]u32{pkt2};
    const names = [_][]const u8{ "dna_seq1", "dna_seq2" };
    const pkt_seqs = [_][]const u32{ &packets1, &packets2 };
    const lens = [_]u64{ 4, 15 };

    try writeTestEaselDb(tmp.dir, "dnatest", 0x11223344, 1, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/dnatest", .{dir_path});
    defer allocator.free(full_path);

    var db = try EaselDsqData.open(allocator, full_path);
    defer db.deinit();

    // Read first sequence (5-bit)
    var seq1 = (try db.readNext(allocator, &alphabet_mod.dna)).?;
    defer seq1.deinit();

    try std.testing.expectEqualStrings("dna_seq1", seq1.name);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1, 2, 3 }, seq1.dsq);

    // Read second sequence (2-bit, 15 A's)
    var seq2 = (try db.readNext(allocator, &alphabet_mod.dna)).?;
    defer seq2.deinit();

    try std.testing.expectEqualStrings("dna_seq2", seq2.name);
    const expected_15_zeros = [_]u8{0} ** 15;
    try std.testing.expectEqualSlices(u8, &expected_15_zeros, seq2.dsq);

    // No more
    const next = try db.readNext(allocator, &alphabet_mod.dna);
    try std.testing.expect(next == null);
}
