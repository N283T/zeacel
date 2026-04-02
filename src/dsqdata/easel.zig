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

// .dsqi header: 7 × u32 (28 bytes) + 3 × u64 (24 bytes) = 52 bytes. No padding.
const IDX_HEADER_SIZE: u64 = 52;
// .dsqm / .dsqs header: magic(u32) + uniquetag(u32) = 8 bytes.
const DATA_HEADER_SIZE: u64 = 8;
// Each index record: metadata_end(i64) + psq_end(i64) = 16 bytes.
const IDX_RECORD_SIZE: u64 = 16;

pub const EaselDsqData = struct {
    allocator: Allocator,
    idx_file: std.fs.File, // .dsqi
    meta_file: std.fs.File, // .dsqm
    seq_file: std.fs.File, // .dsqs

    // From .dsqi header (52 bytes)
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
        // Open stub file and read uniquetag.
        const uniquetag = try readStubTag(allocator, basename);

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

        // Read .dsqi header (52 bytes): 7×u32 + 3×u64, no padding.
        var idx_hdr: [IDX_HEADER_SIZE]u8 = undefined;
        _ = idx_file.preadAll(&idx_hdr, 0) catch return error.InvalidFormat;

        const magic = std.mem.readInt(u32, idx_hdr[0..4], .little);
        if (magic == MAGIC_SWAPPED) return error.ByteswapNotSupported;
        if (magic != MAGIC) return error.InvalidFormat;

        const idx_tag = std.mem.readInt(u32, idx_hdr[4..8], .little);
        if (idx_tag != uniquetag) return error.TagMismatch;

        const alphatype_raw = std.mem.readInt(u32, idx_hdr[8..12], .little);
        // idx_hdr[12..16] = flags (ignored)
        const max_name_len = std.mem.readInt(u32, idx_hdr[16..20], .little);
        const max_acc_len = std.mem.readInt(u32, idx_hdr[20..24], .little);
        const max_desc_len = std.mem.readInt(u32, idx_hdr[24..28], .little);
        // No padding — u64s start at offset 28.
        const max_seq_len = std.mem.readInt(u64, idx_hdr[28..36], .little);
        const nseq = std.mem.readInt(u64, idx_hdr[36..44], .little);
        const nres = std.mem.readInt(u64, idx_hdr[44..52], .little);

        // Read .dsqm header (8 bytes: magic + uniquetag)
        var meta_hdr: [DATA_HEADER_SIZE]u8 = undefined;
        _ = meta_file.preadAll(&meta_hdr, 0) catch return error.InvalidFormat;

        const meta_magic = std.mem.readInt(u32, meta_hdr[0..4], .little);
        if (meta_magic == MAGIC_SWAPPED) return error.ByteswapNotSupported;
        if (meta_magic != MAGIC) return error.InvalidFormat;
        const meta_tag = std.mem.readInt(u32, meta_hdr[4..8], .little);
        if (meta_tag != uniquetag) return error.TagMismatch;

        // Read .dsqs header (8 bytes: magic + uniquetag)
        var seq_hdr: [DATA_HEADER_SIZE]u8 = undefined;
        _ = seq_file.preadAll(&seq_hdr, 0) catch return error.InvalidFormat;

        const seq_magic = std.mem.readInt(u32, seq_hdr[0..4], .little);
        if (seq_magic == MAGIC_SWAPPED) return error.ByteswapNotSupported;
        if (seq_magic != MAGIC) return error.InvalidFormat;
        const seq_tag = std.mem.readInt(u32, seq_hdr[4..8], .little);
        if (seq_tag != uniquetag) return error.TagMismatch;

        // Seek past the headers for sequential reading.
        // preadAll does not advance file position, so we seek manually.
        meta_file.seekTo(DATA_HEADER_SIZE) catch return error.InvalidFormat;
        seq_file.seekTo(DATA_HEADER_SIZE) catch return error.InvalidFormat;

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

    /// Seek to and read a specific sequence by index (0-based).
    /// Returns null if seq_idx is out of range.
    /// Does NOT update current_seq (random access doesn't affect sequential cursor).
    pub fn readSequence(self: *EaselDsqData, allocator: Allocator, abc: *const Alphabet, seq_idx: u64) !?Sequence {
        if (seq_idx >= self.num_sequences) return null;

        // Index records start immediately after the 52-byte header.
        // Each record is 16 bytes: metadata_end(i64) + psq_end(i64).
        //
        // metadata_end: 0-based cumulative byte end index into the metadata data
        //   section (the region of .dsqm after the 8-byte file header).
        // psq_end: 0-based cumulative packet (u32) end index into the sequence
        //   data section (the region of .dsqs after the 8-byte file header).
        //
        // For seq i:
        //   meta data occupies bytes [meta_start_data .. meta_end_data_inclusive]
        //     where meta_start_data = (i==0) ? 0 : idx[i-1].metadata_end + 1
        //           meta_end_data   = idx[i].metadata_end        (inclusive end)
        //   psq data occupies packets [psq_start_pkt .. psq_end_pkt_inclusive]
        //     where psq_start_pkt = (i==0) ? 0 : idx[i-1].psq_end + 1
        //           psq_end_pkt   = idx[i].psq_end                (inclusive end)
        //
        // File offsets add DATA_HEADER_SIZE (8) to the data-section offsets.

        const rec_offset = IDX_HEADER_SIZE + seq_idx * IDX_RECORD_SIZE;

        var idx_buf: [IDX_RECORD_SIZE]u8 = undefined;
        _ = try self.idx_file.preadAll(&idx_buf, rec_offset);

        const meta_end_data = std.mem.readInt(i64, idx_buf[0..8], .little);
        const psq_end_pkt = std.mem.readInt(i64, idx_buf[8..16], .little);

        var meta_start_data: i64 = 0;
        var psq_start_pkt: i64 = 0;

        if (seq_idx > 0) {
            const prev_rec_offset = IDX_HEADER_SIZE + (seq_idx - 1) * IDX_RECORD_SIZE;
            var prev_buf: [IDX_RECORD_SIZE]u8 = undefined;
            _ = try self.idx_file.preadAll(&prev_buf, prev_rec_offset);
            meta_start_data = std.mem.readInt(i64, prev_buf[0..8], .little) + 1;
            psq_start_pkt = std.mem.readInt(i64, prev_buf[8..16], .little) + 1;
        }

        // Convert data-section offsets to absolute file offsets.
        const meta_file_start: u64 = @intCast(meta_start_data + @as(i64, @intCast(DATA_HEADER_SIZE)));
        const meta_file_end: u64 = @intCast(meta_end_data + @as(i64, @intCast(DATA_HEADER_SIZE)) + 1);
        const meta_len = meta_file_end - meta_file_start;

        const psq_file_start: u64 = @intCast(psq_start_pkt * 4 + @as(i64, @intCast(DATA_HEADER_SIZE)));
        const psq_num_pkts: u64 = @intCast(psq_end_pkt - psq_start_pkt + 1);
        const psq_byte_len: u64 = psq_num_pkts * 4;

        // Read metadata via positional read
        const meta_data = try allocator.alloc(u8, meta_len);
        defer allocator.free(meta_data);
        _ = try self.meta_file.preadAll(meta_data, meta_file_start);

        // Parse metadata: name\0, accession\0, description\0, taxonomy_id(i32)
        var pos: usize = 0;
        const name_end = std.mem.indexOfScalarPos(u8, meta_data, pos, 0) orelse return error.InvalidFormat;
        const name = try allocator.dupe(u8, meta_data[pos..name_end]);
        errdefer allocator.free(name);
        pos = name_end + 1;

        const acc_end = std.mem.indexOfScalarPos(u8, meta_data, pos, 0) orelse return error.InvalidFormat;
        const acc_raw = try allocator.dupe(u8, meta_data[pos..acc_end]);
        errdefer allocator.free(acc_raw);
        pos = acc_end + 1;

        const desc_end = std.mem.indexOfScalarPos(u8, meta_data, pos, 0) orelse return error.InvalidFormat;
        const desc_raw = try allocator.dupe(u8, meta_data[pos..desc_end]);
        errdefer allocator.free(desc_raw);
        pos = desc_end + 1;

        if (pos + 4 > meta_data.len) return error.InvalidFormat;
        const taxonomy_id = std.mem.readInt(i32, meta_data[pos..][0..4], .little);

        // Read packet data via positional read
        const psq_data = try allocator.alloc(u8, psq_byte_len);
        defer allocator.free(psq_data);
        _ = try self.seq_file.preadAll(psq_data, psq_file_start);

        // Parse packets (each is 4 bytes, little-endian)
        const num_packets = psq_byte_len / 4;
        var packets = try allocator.alloc(u32, num_packets);
        defer allocator.free(packets);
        for (0..num_packets) |i| {
            const offset = i * 4;
            packets[i] = std.mem.readInt(u32, psq_data[offset..][0..4], .little);
        }

        // Count residues and unpack
        const is_amino = (self.alphabet_type == .amino);
        const seq_len = countResiduesFromPackets(packets, is_amino);
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

    /// Read a chunk of sequences starting at seq_idx.
    /// Returns up to max_seqs sequences. Caller owns the returned slice and each Sequence.
    pub fn readChunk(self: *EaselDsqData, allocator: Allocator, abc: *const Alphabet, seq_idx: u64, max_seqs: usize) ![]Sequence {
        if (seq_idx >= self.num_sequences) return allocator.alloc(Sequence, 0);

        const available = self.num_sequences - seq_idx;
        const count: usize = @intCast(@min(available, @as(u64, max_seqs)));

        var sequences = try allocator.alloc(Sequence, count);
        var read_count: usize = 0;
        errdefer {
            for (0..read_count) |i| {
                var s = sequences[i];
                s.deinit();
            }
            allocator.free(sequences);
        }

        for (0..count) |i| {
            const seq = try self.readSequence(allocator, abc, seq_idx + @as(u64, @intCast(i))) orelse break;
            sequences[i] = seq;
            read_count += 1;
        }

        if (read_count < count) {
            sequences = try allocator.realloc(sequences, read_count);
        }

        return sequences;
    }

    pub fn deinit(self: *EaselDsqData) void {
        self.idx_file.close();
        self.meta_file.close();
        self.seq_file.close();
    }
};

/// Read the uniquetag from the Easel stub file.
/// The first line has the format: "Easel dsqdata v1 x<uniquetag>\n"
fn readStubTag(allocator: Allocator, basename: []const u8) !u32 {
    const stub_path = try std.fmt.allocPrint(allocator, "{s}", .{basename});
    defer allocator.free(stub_path);

    const stub_file = std.fs.cwd().openFile(stub_path, .{}) catch return error.InvalidFormat;
    defer stub_file.close();

    var buf: [256]u8 = undefined;
    const reader = stub_file.deprecatedReader();
    const line = (reader.readUntilDelimiterOrEof(&buf, '\n') catch return error.InvalidFormat) orelse return error.InvalidFormat;

    // Expected: "Easel dsqdata v1 x<uniquetag>"
    // Split on whitespace: tokens = ["Easel", "dsqdata", "v1", "x<uniquetag>"]
    var it = std.mem.tokenizeAny(u8, line, " \t\r");

    const tok0 = it.next() orelse return error.InvalidFormat;
    if (!std.mem.eql(u8, tok0, "Easel")) return error.InvalidFormat;

    const tok1 = it.next() orelse return error.InvalidFormat;
    if (!std.mem.eql(u8, tok1, "dsqdata")) return error.InvalidFormat;

    const tok2 = it.next() orelse return error.InvalidFormat;
    if (tok2.len < 2 or tok2[0] != 'v') return error.InvalidFormat;

    const tok3 = it.next() orelse return error.InvalidFormat;
    if (tok3.len < 2 or tok3[0] != 'x') return error.InvalidFormat;

    const tag = std.fmt.parseInt(u32, tok3[1..], 10) catch return error.InvalidFormat;
    return tag;
}

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
/// Creates basename (stub), basename.dsqi, basename.dsqm, basename.dsqs.
///
/// Index record semantics (matching Easel's esl_dsqdata.c):
///   metadata_end: 0-based cumulative end byte index into the metadata data
///     section (the bytes of .dsqm after the 8-byte file header).
///     Computed as: sum of all metadata bytes written so far, minus 1.
///   psq_end: 0-based cumulative end packet index into the sequence data
///     section (the u32 packets of .dsqs after the 8-byte file header).
///     Computed as: total packets written so far, minus 1.
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

    // Compute summary statistics for the header.
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

    // --- Write stub file ---
    // Format: "Easel dsqdata v1 x<uniquetag>\n\n..."
    {
        const file = try dir.createFile(basename, .{});
        defer file.close();
        const w = file.deprecatedWriter();
        try w.print("Easel dsqdata v1 x{d}\n", .{uniquetag});
        try w.print("\n", .{});
    }

    // --- Write .dsqi ---
    // Header: 7 × u32 + 3 × u64 = 52 bytes (no padding).
    // Followed by nseq index records: metadata_end(i64) + psq_end(i64).
    {
        const path = try std.fmt.bufPrint(&path_buf, "{s}.dsqi", .{basename});
        const file = try dir.createFile(path, .{});
        defer file.close();
        const w = file.deprecatedWriter();

        // 52-byte header
        try w.writeInt(u32, MAGIC, .little);
        try w.writeInt(u32, uniquetag, .little);
        try w.writeInt(u32, alphatype, .little);
        try w.writeInt(u32, 0, .little); // flags
        try w.writeInt(u32, max_name_len, .little);
        try w.writeInt(u32, 0, .little); // max_acclen
        try w.writeInt(u32, 0, .little); // max_desclen
        // No padding — u64s follow immediately.
        try w.writeInt(u64, 0, .little); // max_seqlen (not critical for tests)
        try w.writeInt(u64, @intCast(nseq), .little);
        try w.writeInt(u64, total_res, .little);

        // Index records.
        // psq_end  = cumulative packet count - 1 (0-based end index).
        // meta_end = cumulative metadata byte count - 1 (0-based end index).
        // Both count from the start of their data section (after the 8-byte header).
        // metadata per sequence: name\0 + \0(acc) + \0(desc) + 4(tax_id) bytes.
        var meta_cumulative: i64 = 0; // bytes written to metadata data section so far
        var psq_cumulative: i64 = 0; // packets written to seq data section so far

        for (0..nseq) |i| {
            const meta_size: i64 = @intCast(names[i].len + 1 + 1 + 1 + 4);
            meta_cumulative += meta_size;

            const psq_count: i64 = @intCast(packet_seqs[i].len);
            psq_cumulative += psq_count;

            // Easel writes metadata_end first, then psq_end.
            try w.writeInt(i64, meta_cumulative - 1, .little);
            try w.writeInt(i64, psq_cumulative - 1, .little);
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

test "EaselDsqData.readSequence: random access by index" {
    const allocator = std.testing.allocator;

    // Create a 3-sequence amino database.
    // Seq 0: A(0), C(1) -> 2 residues
    const pkt0: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 0) << 25) |
        (@as(u32, 1) << 20) |
        (@as(u32, 31) << 15) |
        (@as(u32, 31) << 10) |
        (@as(u32, 31) << 5) |
        (@as(u32, 31) << 0);

    // Seq 1: D(2), E(3), F(4) -> 3 residues
    const pkt1: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 2) << 25) |
        (@as(u32, 3) << 20) |
        (@as(u32, 4) << 15) |
        (@as(u32, 31) << 10) |
        (@as(u32, 31) << 5) |
        (@as(u32, 31) << 0);

    // Seq 2: G(5), H(6), I(7), K(8) -> 4 residues
    const pkt2: u32 = pack.EOD | pack.FIVEBIT |
        (@as(u32, 5) << 25) |
        (@as(u32, 6) << 20) |
        (@as(u32, 7) << 15) |
        (@as(u32, 8) << 10) |
        (@as(u32, 31) << 5) |
        (@as(u32, 31) << 0);

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const packets0 = [_]u32{pkt0};
    const packets1 = [_]u32{pkt1};
    const packets2 = [_]u32{pkt2};
    const names = [_][]const u8{ "seq_alpha", "seq_beta", "seq_gamma" };
    const pkt_seqs = [_][]const u32{ &packets0, &packets1, &packets2 };
    const lens = [_]u64{ 2, 3, 4 };

    try writeTestEaselDb(tmp.dir, "randtest", 0xDEADBEEF, 3, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/randtest", .{dir_path});
    defer allocator.free(full_path);

    var db = try EaselDsqData.open(allocator, full_path);
    defer db.deinit();

    // Read sequence 2 (0-indexed) first — out of order
    {
        var seq = (try db.readSequence(allocator, &alphabet_mod.amino, 2)).?;
        defer seq.deinit();
        try std.testing.expectEqualStrings("seq_gamma", seq.name);
        try std.testing.expectEqualSlices(u8, &[_]u8{ 5, 6, 7, 8 }, seq.dsq);
    }

    // Read sequence 0
    {
        var seq = (try db.readSequence(allocator, &alphabet_mod.amino, 0)).?;
        defer seq.deinit();
        try std.testing.expectEqualStrings("seq_alpha", seq.name);
        try std.testing.expectEqualSlices(u8, &[_]u8{ 0, 1 }, seq.dsq);
    }

    // Read sequence 1
    {
        var seq = (try db.readSequence(allocator, &alphabet_mod.amino, 1)).?;
        defer seq.deinit();
        try std.testing.expectEqualStrings("seq_beta", seq.name);
        try std.testing.expectEqualSlices(u8, &[_]u8{ 2, 3, 4 }, seq.dsq);
    }

    // Out-of-range returns null
    const out = try db.readSequence(allocator, &alphabet_mod.amino, 3);
    try std.testing.expect(out == null);

    // Verify sequential cursor is unaffected
    try std.testing.expectEqual(@as(u64, 0), db.current_seq);
}

test "EaselDsqData.readChunk: read subset of sequences" {
    const allocator = std.testing.allocator;

    // Create a 5-sequence amino database.
    const make_pkt = struct {
        fn f(code: u32) u32 {
            return pack.EOD | pack.FIVEBIT |
                (code << 25) |
                (@as(u32, 31) << 20) |
                (@as(u32, 31) << 15) |
                (@as(u32, 31) << 10) |
                (@as(u32, 31) << 5) |
                (@as(u32, 31) << 0);
        }
    }.f;

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const pkts = [5][1]u32{
        .{make_pkt(0)},
        .{make_pkt(1)},
        .{make_pkt(2)},
        .{make_pkt(3)},
        .{make_pkt(4)},
    };

    const names = [_][]const u8{ "s0", "s1", "s2", "s3", "s4" };
    const pkt_seqs = [_][]const u32{ &pkts[0], &pkts[1], &pkts[2], &pkts[3], &pkts[4] };
    const lens = [_]u64{ 1, 1, 1, 1, 1 };

    try writeTestEaselDb(tmp.dir, "chunktest", 0xCAFEBABE, 3, &names, &pkt_seqs, &lens);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_path = try tmp.dir.realpath(".", &path_buf);
    const full_path = try std.fmt.allocPrint(allocator, "{s}/chunktest", .{dir_path});
    defer allocator.free(full_path);

    var db = try EaselDsqData.open(allocator, full_path);
    defer db.deinit();

    // Read chunk of 2 starting at index 1
    const chunk = try db.readChunk(allocator, &alphabet_mod.amino, 1, 2);
    defer {
        for (chunk) |*s| {
            var seq = s.*;
            seq.deinit();
        }
        allocator.free(chunk);
    }

    try std.testing.expectEqual(@as(usize, 2), chunk.len);
    try std.testing.expectEqualStrings("s1", chunk[0].name);
    try std.testing.expectEqualSlices(u8, &[_]u8{1}, chunk[0].dsq);
    try std.testing.expectEqualStrings("s2", chunk[1].name);
    try std.testing.expectEqualSlices(u8, &[_]u8{2}, chunk[1].dsq);

    // Read chunk that extends past end
    const tail = try db.readChunk(allocator, &alphabet_mod.amino, 3, 10);
    defer {
        for (tail) |*s| {
            var seq = s.*;
            seq.deinit();
        }
        allocator.free(tail);
    }

    try std.testing.expectEqual(@as(usize, 2), tail.len);
    try std.testing.expectEqualStrings("s3", tail[0].name);
    try std.testing.expectEqualStrings("s4", tail[1].name);

    // Read chunk at out-of-range index returns empty slice
    const empty = try db.readChunk(allocator, &alphabet_mod.amino, 100, 5);
    defer allocator.free(empty);
    try std.testing.expectEqual(@as(usize, 0), empty.len);
}
