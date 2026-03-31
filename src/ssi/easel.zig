// EaselIndex — parser for Easel SSI v3.0 binary index files.
//
// The Easel SSI format is big-endian by default but may be byteswapped
// if the index was created on a little-endian machine. The magic value
// determines endianness: 0xd3d3c9b3 (native/big-endian) or 0xb3c9d3d3
// (byteswapped/little-endian).

const std = @import("std");
const Allocator = std.mem.Allocator;

const ssi = @import("../ssi.zig");

pub const EASEL_MAGIC: u32 = 0xd3d3c9b3;
pub const EASEL_MAGIC_SWAPPED: u32 = 0xb3c9d3d3;

pub const FileInfo = struct {
    name: []const u8,
    format: u32,
    flags: u32,
    bpl: u32, // bytes per line
    rpl: u32, // residues per line
};

pub const EaselIndex = struct {
    file: std.fs.File,
    byteswap: bool,
    offsz: u8,
    nfiles: u16,
    nprimary: u64,
    nsecondary: u64,
    plen: u32,
    slen: u32,
    precsize: u32,
    srecsize: u32,
    poffset: u64,
    soffset: u64,
    files: []FileInfo,
    allocator: Allocator,

    // Maximum header size: fixed (54 bytes) + 3 * 8-byte offsets = 78 bytes.
    const max_header_size = 78;

    /// Parse an Easel SSI v3.0 index from an open file handle.
    /// The file handle is kept open for later disk-based lookups.
    /// Caller owns the returned EaselIndex and must call deinit().
    pub fn read(allocator: Allocator, file: std.fs.File) !EaselIndex {
        // Read the maximum possible header into a stack buffer.
        var hdr_buf: [max_header_size]u8 = undefined;
        const n_read = try file.preadAll(&hdr_buf, 0);
        if (n_read < 14) return error.InvalidFormat; // Need at least magic + flags + offsz + nfiles

        var pos: usize = 0;

        // Read magic (always big-endian on disk).
        const raw_magic = readIntFromBuf(u32, hdr_buf[pos..], .big);
        pos += 4;

        const byteswap = if (raw_magic == EASEL_MAGIC)
            false
        else if (raw_magic == EASEL_MAGIC_SWAPPED)
            true
        else
            return error.InvalidFormat;

        const endian: std.builtin.Endian = if (byteswap) .little else .big;

        // flags (skip)
        pos += 4;

        const offsz_raw = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        if (offsz_raw != 4 and offsz_raw != 8) return error.InvalidFormat;
        const offsz: u8 = @intCast(offsz_raw);

        // Verify we read enough header bytes for the variable-size offsets.
        const actual_header_size: usize = 54 + 3 * @as(usize, offsz);
        if (n_read < actual_header_size) return error.InvalidFormat;

        const nfiles = readIntFromBuf(u16, hdr_buf[pos..], endian);
        pos += 2;
        const nprimary = readIntFromBuf(u64, hdr_buf[pos..], endian);
        pos += 8;
        const nsecondary = readIntFromBuf(u64, hdr_buf[pos..], endian);
        pos += 8;
        const flen = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const plen = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const slen = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const frecsize = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const precsize = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;
        const srecsize = readIntFromBuf(u32, hdr_buf[pos..], endian);
        pos += 4;

        const foffset = readOffsetFromBuf(hdr_buf[pos..], offsz, endian);
        pos += offsz;
        const poffset = readOffsetFromBuf(hdr_buf[pos..], offsz, endian);
        pos += offsz;
        const soffset = readOffsetFromBuf(hdr_buf[pos..], offsz, endian);

        // Read file info section from disk.
        const file_section_size = @as(usize, frecsize) * @as(usize, nfiles);
        const file_section_buf = try allocator.alloc(u8, file_section_size);
        defer allocator.free(file_section_buf);

        const file_n_read = try file.preadAll(file_section_buf, foffset);
        if (file_n_read < file_section_size) return error.InvalidFormat;

        // Parse file info entries.
        const files = try allocator.alloc(FileInfo, nfiles);
        var files_done: usize = 0;
        errdefer {
            for (files[0..files_done]) |f| allocator.free(f.name);
            allocator.free(files);
        }

        var fpos: usize = 0;
        for (0..nfiles) |i| {
            const rec_start = fpos;
            const name_bytes = file_section_buf[fpos..][0..flen];
            fpos += flen;

            // Find actual string length (up to first null).
            var name_len: usize = 0;
            for (name_bytes) |c| {
                if (c == 0) break;
                name_len += 1;
            }
            const name = try allocator.dupe(u8, name_bytes[0..name_len]);
            errdefer allocator.free(name);

            const format = readIntFromBuf(u32, file_section_buf[fpos..], endian);
            fpos += 4;
            const file_flags = readIntFromBuf(u32, file_section_buf[fpos..], endian);
            fpos += 4;
            const bpl = readIntFromBuf(u32, file_section_buf[fpos..], endian);
            fpos += 4;
            const rpl = readIntFromBuf(u32, file_section_buf[fpos..], endian);
            fpos += 4;

            // Advance to next record using frecsize, which may be larger than
            // flen + 16 if a future Easel version adds extra fields.
            fpos = rec_start + @as(usize, frecsize);

            files[i] = FileInfo{
                .name = name,
                .format = format,
                .flags = file_flags,
                .bpl = bpl,
                .rpl = rpl,
            };
            files_done += 1;
        }

        return EaselIndex{
            .file = file,
            .byteswap = byteswap,
            .offsz = offsz,
            .nfiles = nfiles,
            .nprimary = nprimary,
            .nsecondary = nsecondary,
            .plen = plen,
            .slen = slen,
            .precsize = precsize,
            .srecsize = srecsize,
            .poffset = poffset,
            .soffset = soffset,
            .files = files,
            .allocator = allocator,
        };
    }

    /// Free all allocated memory and close the file handle.
    pub fn deinit(self: *EaselIndex) void {
        for (self.files) |f| self.allocator.free(f.name);
        self.allocator.free(self.files);
        self.file.close();
    }
};

/// Read an integer of type T from a byte slice at the given endianness.
fn readIntFromBuf(comptime T: type, buf: []const u8, endian: std.builtin.Endian) T {
    const n = @divExact(@typeInfo(T).int.bits, 8);
    return std.mem.readInt(T, buf[0..n], endian);
}

/// Read a 4-byte or 8-byte offset from a byte slice.
fn readOffsetFromBuf(buf: []const u8, offsz: u8, endian: std.builtin.Endian) u64 {
    return switch (offsz) {
        4 => @as(u64, readIntFromBuf(u32, buf, endian)),
        8 => readIntFromBuf(u64, buf, endian),
        else => unreachable,
    };
}

/// Write a valid Easel SSI v3.0 header and file info section for testing.
/// Calculates section offsets automatically based on the provided parameters.
pub fn writeEaselHeader(
    buf: *std.ArrayList(u8),
    allocator: Allocator,
    opts: struct {
        nfiles: u16 = 1,
        nprimary: u64 = 0,
        nsecondary: u64 = 0,
        flen: u32 = 32,
        plen: u32 = 32,
        slen: u32 = 32,
        offsz: u8 = 8,
        byteswap: bool = false,
        file_names: []const []const u8 = &.{"testfile.fa"},
        file_formats: []const u32 = &.{0},
        file_flags: []const u32 = &.{0},
        file_bpls: []const u32 = &.{80},
        file_rpls: []const u32 = &.{60},
    },
) !void {
    const writer = buf.writer(allocator);

    const endian: std.builtin.Endian = if (opts.byteswap) .little else .big;
    const offsz: u8 = opts.offsz;

    // Calculate sizes.
    const frecsize: u32 = opts.flen + 16; // flen + 4*u32
    const precsize: u32 = opts.plen + 2 + 2 * @as(u32, offsz) + 8;
    const srecsize: u32 = opts.slen + opts.plen;

    // Header size: fixed fields + 3 offsets.
    const header_size: u64 = 4 + 4 + 4 + 2 + 8 + 8 + 4 + 4 + 4 + 4 + 4 + 4 + 3 * @as(u64, offsz);
    const foffset: u64 = header_size;
    const poffset: u64 = foffset + @as(u64, frecsize) * @as(u64, opts.nfiles);
    const soffset: u64 = poffset + @as(u64, precsize) * opts.nprimary;

    // Write magic.
    const magic: u32 = if (opts.byteswap) EASEL_MAGIC_SWAPPED else EASEL_MAGIC;
    try writer.writeInt(u32, magic, .big);

    // flags
    try writer.writeInt(u32, 0, endian);
    // offsz
    try writer.writeInt(u32, @as(u32, offsz), endian);
    // nfiles
    try writer.writeInt(u16, opts.nfiles, endian);
    // nprimary
    try writer.writeInt(u64, opts.nprimary, endian);
    // nsecondary
    try writer.writeInt(u64, opts.nsecondary, endian);
    // flen, plen, slen
    try writer.writeInt(u32, opts.flen, endian);
    try writer.writeInt(u32, opts.plen, endian);
    try writer.writeInt(u32, opts.slen, endian);
    // frecsize, precsize, srecsize
    try writer.writeInt(u32, frecsize, endian);
    try writer.writeInt(u32, precsize, endian);
    try writer.writeInt(u32, srecsize, endian);

    // offsets (foffset, poffset, soffset)
    try writeOffset(writer, foffset, offsz, endian);
    try writeOffset(writer, poffset, offsz, endian);
    try writeOffset(writer, soffset, offsz, endian);

    // File info section.
    for (0..opts.nfiles) |i| {
        const name = if (i < opts.file_names.len) opts.file_names[i] else "unknown";
        const fmt = if (i < opts.file_formats.len) opts.file_formats[i] else 0;
        const flg = if (i < opts.file_flags.len) opts.file_flags[i] else 0;
        const bpl = if (i < opts.file_bpls.len) opts.file_bpls[i] else 80;
        const rpl = if (i < opts.file_rpls.len) opts.file_rpls[i] else 60;

        // Write null-padded filename.
        try writer.writeAll(name);
        const pad_len = opts.flen - @as(u32, @intCast(name.len));
        for (0..pad_len) |_| try writer.writeByte(0);

        try writer.writeInt(u32, fmt, endian);
        try writer.writeInt(u32, flg, endian);
        try writer.writeInt(u32, bpl, endian);
        try writer.writeInt(u32, rpl, endian);
    }
}

pub fn writeOffset(writer: anytype, value: u64, offsz: u8, endian: std.builtin.Endian) !void {
    switch (offsz) {
        4 => try writer.writeInt(u32, @intCast(value), endian),
        8 => try writer.writeInt(u64, value, endian),
        else => return error.InvalidFormat,
    }
}

// --- Tests ---

test "EaselIndex.read: parses header and file info" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{
        .nfiles = 1,
        .nprimary = 3,
        .nsecondary = 0,
        .flen = 32,
        .plen = 32,
        .slen = 32,
        .offsz = 8,
        .file_names = &.{"sequences.fasta"},
        .file_bpls = &.{80},
        .file_rpls = &.{60},
    });

    // Write to a temp file since EaselIndex needs a seekable file handle.
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(.{ .sub_path = "test.ssi", .data = buf.items });
    const file = try tmp_dir.dir.openFile("test.ssi", .{});
    // File will be closed by deinit.

    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    try std.testing.expectEqual(false, idx.byteswap);
    try std.testing.expectEqual(@as(u8, 8), idx.offsz);
    try std.testing.expectEqual(@as(u16, 1), idx.nfiles);
    try std.testing.expectEqual(@as(u64, 3), idx.nprimary);
    try std.testing.expectEqual(@as(u64, 0), idx.nsecondary);
    try std.testing.expectEqual(@as(u32, 32), idx.plen);
    try std.testing.expectEqual(@as(u32, 32), idx.slen);

    // File info.
    try std.testing.expectEqual(@as(usize, 1), idx.files.len);
    try std.testing.expectEqualStrings("sequences.fasta", idx.files[0].name);
    try std.testing.expectEqual(@as(u32, 80), idx.files[0].bpl);
    try std.testing.expectEqual(@as(u32, 60), idx.files[0].rpl);
}

test "EaselIndex.read: rejects invalid magic" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // Write "ZSSI" followed by padding — not a valid Easel magic.
    var bad_data: [64]u8 = undefined;
    @memset(&bad_data, 0);
    @memcpy(bad_data[0..4], "ZSSI");

    try tmp_dir.dir.writeFile(.{ .sub_path = "bad.ssi", .data = &bad_data });
    const file = try tmp_dir.dir.openFile("bad.ssi", .{});
    defer file.close();

    try std.testing.expectError(error.InvalidFormat, EaselIndex.read(allocator, file));
}

test "EaselIndex.read: parses byteswapped header" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{
        .nfiles = 1,
        .nprimary = 5,
        .nsecondary = 2,
        .byteswap = true,
        .file_names = &.{"swapped.fa"},
        .file_bpls = &.{100},
        .file_rpls = &.{80},
    });

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(.{ .sub_path = "swap.ssi", .data = buf.items });
    const file = try tmp_dir.dir.openFile("swap.ssi", .{});

    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    try std.testing.expectEqual(true, idx.byteswap);
    try std.testing.expectEqual(@as(u64, 5), idx.nprimary);
    try std.testing.expectEqual(@as(u64, 2), idx.nsecondary);
    try std.testing.expectEqualStrings("swapped.fa", idx.files[0].name);
    try std.testing.expectEqual(@as(u32, 100), idx.files[0].bpl);
    try std.testing.expectEqual(@as(u32, 80), idx.files[0].rpl);
}

test "EaselIndex.read: parses header with 4-byte offsets" {
    const allocator = std.testing.allocator;

    var buf = std.ArrayList(u8){};
    defer buf.deinit(allocator);

    try writeEaselHeader(&buf, allocator, .{
        .nfiles = 2,
        .nprimary = 10,
        .nsecondary = 3,
        .flen = 16,
        .plen = 16,
        .slen = 16,
        .offsz = 4,
        .file_names = &.{ "small1.fa", "small2.fa" },
        .file_formats = &.{ 1, 2 },
        .file_flags = &.{ 0, 0 },
        .file_bpls = &.{ 60, 70 },
        .file_rpls = &.{ 50, 55 },
    });

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(.{ .sub_path = "off4.ssi", .data = buf.items });
    const file = try tmp_dir.dir.openFile("off4.ssi", .{});

    var idx = try EaselIndex.read(allocator, file);
    defer idx.deinit();

    try std.testing.expectEqual(false, idx.byteswap);
    try std.testing.expectEqual(@as(u8, 4), idx.offsz);
    try std.testing.expectEqual(@as(u16, 2), idx.nfiles);
    try std.testing.expectEqual(@as(u64, 10), idx.nprimary);
    try std.testing.expectEqual(@as(u64, 3), idx.nsecondary);
    try std.testing.expectEqual(@as(u32, 16), idx.plen);
    try std.testing.expectEqual(@as(u32, 16), idx.slen);

    // File info.
    try std.testing.expectEqual(@as(usize, 2), idx.files.len);
    try std.testing.expectEqualStrings("small1.fa", idx.files[0].name);
    try std.testing.expectEqual(@as(u32, 1), idx.files[0].format);
    try std.testing.expectEqual(@as(u32, 60), idx.files[0].bpl);
    try std.testing.expectEqual(@as(u32, 50), idx.files[0].rpl);
    try std.testing.expectEqualStrings("small2.fa", idx.files[1].name);
    try std.testing.expectEqual(@as(u32, 2), idx.files[1].format);
    try std.testing.expectEqual(@as(u32, 70), idx.files[1].bpl);
    try std.testing.expectEqual(@as(u32, 55), idx.files[1].rpl);
}
