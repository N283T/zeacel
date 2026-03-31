// Windowed/streaming sequence reader for chromosome-scale sequences.
//
// Reads FASTA sequences in overlapping windows without loading entire
// sequences into memory. Designed for HMMER-style sliding-window search
// over large genomes.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Alphabet = @import("../alphabet.zig").Alphabet;
const Sequence = @import("../sequence.zig").Sequence;

/// Result of reading one window from a sequence.
pub const WindowResult = struct {
    /// The window sequence (digitally encoded). Caller owns this.
    seq: Sequence,
    /// 1-based start coordinate in the full source sequence.
    window_start: i64,
    /// 1-based end coordinate in the full source sequence.
    window_end: i64,
    /// True if this is the last window for the current sequence.
    complete: bool,
};

/// Streaming reader that yields overlapping windows of residues from a
/// FASTA file. Only a single window's worth of data is held in memory
/// at any time, making it suitable for chromosome-scale sequences.
pub const WindowReader = struct {
    abc: *const Alphabet,
    allocator: Allocator,
    file: std.fs.File,
    /// Number of new residues per window.
    window_size: usize,
    /// Number of overlap residues carried from the previous window.
    context_size: usize,

    // Manual read buffer
    buf: [4096]u8 = undefined,
    buf_len: usize = 0,
    buf_pos: usize = 0,

    // Internal state
    /// Name of the current sequence (owned).
    seq_name: ?[]u8 = null,
    /// Description of the current sequence (owned).
    seq_desc: ?[]u8 = null,
    /// Context residues carried over from the previous window (digital codes).
    context_buf: ?[]u8 = null,
    /// Running count of residues read so far in the current source sequence.
    seq_pos: i64 = 0,
    /// True when we have consumed the header for the current sequence.
    header_parsed: bool = false,
    /// True when the current sequence's residues are exhausted (hit '>' or EOF).
    seq_done: bool = false,
    /// A single byte that was read ahead and needs to be returned first.
    pending_byte: ?u8 = null,
    /// True when the underlying file has reached EOF.
    eof: bool = false,

    /// Open a FASTA file for windowed reading.
    ///
    /// `window_size` is the number of *new* residues per window (W).
    /// `context_size` is the number of overlap residues between consecutive
    /// windows (C). The first window contains up to W+C residues; subsequent
    /// windows contain C context residues followed by up to W new residues.
    pub fn init(
        allocator: Allocator,
        abc: *const Alphabet,
        path: []const u8,
        window_size: usize,
        context_size: usize,
    ) !WindowReader {
        const file = try std.fs.cwd().openFile(path, .{});
        errdefer file.close();

        return WindowReader{
            .abc = abc,
            .allocator = allocator,
            .file = file,
            .window_size = window_size,
            .context_size = context_size,
        };
    }

    /// Return the next window of residues, or null if all sequences have
    /// been consumed.
    pub fn nextWindow(self: *WindowReader) !?WindowResult {
        // If the previous sequence is done, advance to the next header.
        if (self.seq_done) {
            self.resetSequenceState();
        }

        // Parse header if needed.
        if (!self.header_parsed) {
            if (!try self.parseHeader()) return null; // EOF, no more sequences
        }

        // Determine how many new residues to collect.
        const target = if (self.seq_pos == 0)
            self.window_size + self.context_size // first window: W+C
        else
            self.window_size; // subsequent windows: W new residues

        // Build the window buffer.
        var dsq_buf = std.ArrayList(u8){};
        defer dsq_buf.deinit(self.allocator);

        // Prepend context from previous window.
        if (self.context_buf) |ctx| {
            try dsq_buf.appendSlice(self.allocator, ctx);
        }

        var new_count: usize = 0;
        var complete = false;

        // Read residues from the file.
        while (new_count < target) {
            const byte = self.nextByte() orelse {
                complete = true;
                break;
            };

            if (byte == '>') {
                // Hit the start of the next sequence.
                self.seq_done = true;
                complete = true;
                break;
            }

            // Skip whitespace.
            if (byte == '\n' or byte == '\r' or byte == ' ' or byte == '\t') continue;

            // Digitize the character.
            const code = self.abc.encode(byte) catch return error.InvalidCharacter;
            try dsq_buf.append(self.allocator, code);
            new_count += 1;
        }

        // If we got no new residues and no context, this sequence produced
        // nothing (empty or already fully consumed). Move to next sequence.
        if (new_count == 0 and dsq_buf.items.len == 0) {
            self.seq_done = true;
            self.resetSequenceState();
            return self.nextWindow();
        }

        // If we filled the window exactly, peek ahead to determine if
        // the sequence is complete.
        if (!complete) {
            complete = self.peekForSequenceEnd();
        }

        const total_residues = dsq_buf.items.len;

        // Compute 1-based coordinates.
        const context_len: i64 = if (self.context_buf) |ctx| @intCast(ctx.len) else 0;
        const window_start = self.seq_pos + 1 - context_len;
        const window_end = self.seq_pos + @as(i64, @intCast(new_count));

        // Advance the running position by the number of NEW residues.
        self.seq_pos += @as(i64, @intCast(new_count));

        // Save context for next window: the last C residues of this window.
        if (self.context_buf) |old| self.allocator.free(old);
        if (self.context_size > 0 and total_residues >= self.context_size) {
            self.context_buf = try self.allocator.dupe(u8, dsq_buf.items[total_residues - self.context_size ..]);
        } else if (self.context_size > 0 and total_residues > 0) {
            self.context_buf = try self.allocator.dupe(u8, dsq_buf.items);
        } else {
            self.context_buf = null;
        }

        if (complete) {
            self.seq_done = true;
        }

        // Build the Sequence.
        const dsq = try self.allocator.dupe(u8, dsq_buf.items);
        errdefer self.allocator.free(dsq);

        const name_copy = try self.allocator.dupe(u8, self.seq_name.?);
        errdefer self.allocator.free(name_copy);

        const src_name = try self.allocator.dupe(u8, self.seq_name.?);
        errdefer self.allocator.free(src_name);

        var desc_copy: ?[]const u8 = null;
        if (self.seq_desc) |d| {
            desc_copy = try self.allocator.dupe(u8, d);
        }
        errdefer if (desc_copy) |d| self.allocator.free(d);

        const full_length: i64 = if (complete) self.seq_pos else 0;

        return WindowResult{
            .seq = Sequence{
                .name = name_copy,
                .accession = null,
                .description = desc_copy,
                .taxonomy_id = null,
                .dsq = dsq,
                .secondary_structure = null,
                .source = Sequence.Source{
                    .name = src_name,
                    .start = if (window_start < 1) 1 else window_start,
                    .end = window_end,
                    .full_length = full_length,
                },
                .abc = self.abc,
                .allocator = self.allocator,
            },
            .window_start = if (window_start < 1) 1 else window_start,
            .window_end = window_end,
            .complete = complete,
        };
    }

    /// Free all owned resources and close the underlying file.
    pub fn deinit(self: *WindowReader) void {
        if (self.seq_name) |n| self.allocator.free(n);
        if (self.seq_desc) |d| self.allocator.free(d);
        if (self.context_buf) |b| self.allocator.free(b);
        self.file.close();
    }

    // --- Internal helpers ---

    /// Read a single byte from the file, using manual buffering.
    /// Returns null on EOF.
    fn nextByte(self: *WindowReader) ?u8 {
        if (self.pending_byte) |b| {
            self.pending_byte = null;
            return b;
        }
        if (self.buf_pos < self.buf_len) {
            const b = self.buf[self.buf_pos];
            self.buf_pos += 1;
            return b;
        }
        if (self.eof) return null;
        // Refill buffer.
        const n = self.file.read(&self.buf) catch return null;
        if (n == 0) {
            self.eof = true;
            return null;
        }
        self.buf_len = n;
        self.buf_pos = 1;
        return self.buf[0];
    }

    /// Peek ahead through whitespace to determine if the current sequence
    /// is complete (next non-whitespace is '>' or EOF). If a residue
    /// character is found, it is saved as pending_byte for the next read.
    fn peekForSequenceEnd(self: *WindowReader) bool {
        while (true) {
            const byte = self.nextByte() orelse return true; // EOF
            if (byte == '\n' or byte == '\r' or byte == ' ' or byte == '\t') continue;
            if (byte == '>') {
                self.seq_done = true;
                return true;
            }
            // It's a residue character; push it back.
            self.pending_byte = byte;
            return false;
        }
    }

    /// Parse the next FASTA header line. Returns false if EOF with no header found.
    fn parseHeader(self: *WindowReader) !bool {
        // If seq_done was set because we hit '>', the '>' has already been
        // consumed by the read loop. Otherwise, skip to the first '>'.
        if (!self.seq_done) {
            while (true) {
                const byte = self.nextByte() orelse return false;
                if (byte == '>') break;
            }
        }

        // Read the header line.
        var header_buf = std.ArrayList(u8){};
        defer header_buf.deinit(self.allocator);

        while (true) {
            const byte = self.nextByte() orelse break;
            if (byte == '\n') break;
            if (byte == '\r') {
                // Consume optional \n after \r.
                const next = self.nextByte();
                if (next) |n| {
                    if (n != '\n') self.pending_byte = n;
                }
                break;
            }
            try header_buf.append(self.allocator, byte);
        }

        if (header_buf.items.len == 0) return false;

        // Split into name and description.
        const header = header_buf.items;
        var name_end: usize = header.len;
        var desc_start: ?usize = null;
        for (header, 0..) |c, i| {
            if (c == ' ' or c == '\t') {
                name_end = i;
                var ds = i + 1;
                while (ds < header.len and (header[ds] == ' ' or header[ds] == '\t')) ds += 1;
                if (ds < header.len) desc_start = ds;
                break;
            }
        }

        if (self.seq_name) |old| self.allocator.free(old);
        self.seq_name = try self.allocator.dupe(u8, header[0..name_end]);

        if (self.seq_desc) |old| self.allocator.free(old);
        if (desc_start) |ds| {
            self.seq_desc = try self.allocator.dupe(u8, header[ds..]);
        } else {
            self.seq_desc = null;
        }

        self.header_parsed = true;
        self.seq_pos = 0;
        self.seq_done = false;
        self.pending_byte = null;
        if (self.context_buf) |old| {
            self.allocator.free(old);
            self.context_buf = null;
        }

        return true;
    }

    /// Reset per-sequence state for advancing to the next sequence.
    /// Preserves seq_done so parseHeader knows whether '>' was already consumed.
    fn resetSequenceState(self: *WindowReader) void {
        self.header_parsed = false;
        // Do NOT reset seq_done here: parseHeader needs to know whether
        // the '>' delimiter was already consumed by the read loop.
        if (self.context_buf) |old| {
            self.allocator.free(old);
            self.context_buf = null;
        }
    }
};

// --- Tests ---

test "WindowReader: single sequence, single window" {
    const allocator = std.testing.allocator;
    const alphabet = @import("../alphabet.zig");

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const file = try tmp_dir.dir.createFile("test.fa", .{});
    try file.writeAll(">seq1 test sequence\nACGTACGT\n");
    file.close();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try tmp_dir.dir.realpath("test.fa", &path_buf);

    var wr = try WindowReader.init(allocator, &alphabet.dna, path, 100, 0);
    defer wr.deinit();

    // Should get a single window with all 8 residues.
    const w1 = (try wr.nextWindow()).?;
    var seq1 = w1.seq;
    defer seq1.deinit();

    try std.testing.expectEqualStrings("seq1", seq1.name);
    try std.testing.expectEqual(@as(usize, 8), seq1.len());
    try std.testing.expectEqual(@as(i64, 1), w1.window_start);
    try std.testing.expectEqual(@as(i64, 8), w1.window_end);
    try std.testing.expect(w1.complete);

    // Description should be parsed.
    try std.testing.expectEqualStrings("test sequence", seq1.description.?);

    // Verify digitized content round-trips.
    const text = try seq1.toText();
    defer allocator.free(text);
    try std.testing.expectEqualStrings("ACGTACGT", text);

    // Source tracking: full_length is known on the last window.
    try std.testing.expectEqual(@as(i64, 8), seq1.source.?.full_length);

    // No more windows.
    try std.testing.expectEqual(@as(?WindowResult, null), try wr.nextWindow());
}

test "WindowReader: overlapping windows" {
    const allocator = std.testing.allocator;
    const alphabet = @import("../alphabet.zig");

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const file = try tmp_dir.dir.createFile("test.fa", .{});
    // 10 residues: ACGTACGTAC
    try file.writeAll(">seq1\nACGTACGTAC\n");
    file.close();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try tmp_dir.dir.realpath("test.fa", &path_buf);

    // W=4, C=2: windows should overlap by 2.
    var wr = try WindowReader.init(allocator, &alphabet.dna, path, 4, 2);
    defer wr.deinit();

    // Window 1: first window gets W+C = 6 residues (ACGTAC), positions 1..6.
    const w1 = (try wr.nextWindow()).?;
    var seq1 = w1.seq;
    defer seq1.deinit();

    try std.testing.expectEqual(@as(usize, 6), seq1.len());
    try std.testing.expectEqual(@as(i64, 1), w1.window_start);
    try std.testing.expectEqual(@as(i64, 6), w1.window_end);
    try std.testing.expect(!w1.complete);

    const text1 = try seq1.toText();
    defer allocator.free(text1);
    try std.testing.expectEqualStrings("ACGTAC", text1);

    // Window 2: C=2 context (AC) + 4 new residues (GTAC) = 6 residues total.
    // Context starts at position 5, new residues at 7..10.
    const w2 = (try wr.nextWindow()).?;
    var seq2 = w2.seq;
    defer seq2.deinit();

    try std.testing.expectEqual(@as(usize, 6), seq2.len());
    try std.testing.expectEqual(@as(i64, 5), w2.window_start);
    try std.testing.expectEqual(@as(i64, 10), w2.window_end);
    try std.testing.expect(w2.complete);

    const text2 = try seq2.toText();
    defer allocator.free(text2);
    try std.testing.expectEqualStrings("ACGTAC", text2);

    // No more windows.
    try std.testing.expectEqual(@as(?WindowResult, null), try wr.nextWindow());
}

test "WindowReader: multiple sequences" {
    const allocator = std.testing.allocator;
    const alphabet = @import("../alphabet.zig");

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const file = try tmp_dir.dir.createFile("test.fa", .{});
    try file.writeAll(">seq1\nACGT\n>seq2\nTGCA\n");
    file.close();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try tmp_dir.dir.realpath("test.fa", &path_buf);

    // Large window to get each sequence in one go.
    var wr = try WindowReader.init(allocator, &alphabet.dna, path, 100, 0);
    defer wr.deinit();

    // Sequence 1.
    const w1 = (try wr.nextWindow()).?;
    var s1 = w1.seq;
    defer s1.deinit();
    try std.testing.expectEqualStrings("seq1", s1.name);
    try std.testing.expect(w1.complete);

    const t1 = try s1.toText();
    defer allocator.free(t1);
    try std.testing.expectEqualStrings("ACGT", t1);

    // Sequence 2.
    const w2 = (try wr.nextWindow()).?;
    var s2 = w2.seq;
    defer s2.deinit();
    try std.testing.expectEqualStrings("seq2", s2.name);
    try std.testing.expect(w2.complete);

    const t2 = try s2.toText();
    defer allocator.free(t2);
    try std.testing.expectEqualStrings("TGCA", t2);

    // Done.
    try std.testing.expectEqual(@as(?WindowResult, null), try wr.nextWindow());
}

test "WindowReader: multi-line FASTA" {
    const allocator = std.testing.allocator;
    const alphabet = @import("../alphabet.zig");

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const file = try tmp_dir.dir.createFile("test.fa", .{});
    try file.writeAll(">seq1\nACGT\nACGT\nAC\n");
    file.close();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try tmp_dir.dir.realpath("test.fa", &path_buf);

    // W=5, C=2.
    var wr = try WindowReader.init(allocator, &alphabet.dna, path, 5, 2);
    defer wr.deinit();

    // Window 1: W+C = 7 residues (ACGTACG), from multi-line input.
    const w1 = (try wr.nextWindow()).?;
    var s1 = w1.seq;
    defer s1.deinit();

    try std.testing.expectEqual(@as(usize, 7), s1.len());
    try std.testing.expectEqual(@as(i64, 1), w1.window_start);
    try std.testing.expectEqual(@as(i64, 7), w1.window_end);
    try std.testing.expect(!w1.complete);

    const t1 = try s1.toText();
    defer allocator.free(t1);
    try std.testing.expectEqualStrings("ACGTACG", t1);

    // Window 2: 2 context (CG) + up to 5 new residues, but only 3 remain (TAC).
    // Total = 5 residues: CGTAC.
    const w2 = (try wr.nextWindow()).?;
    var s2 = w2.seq;
    defer s2.deinit();

    try std.testing.expectEqual(@as(i64, 6), w2.window_start);
    try std.testing.expectEqual(@as(i64, 10), w2.window_end);
    try std.testing.expect(w2.complete);

    const t2 = try s2.toText();
    defer allocator.free(t2);
    try std.testing.expectEqualStrings("CGTAC", t2);
}

test "WindowReader: no context (C=0)" {
    const allocator = std.testing.allocator;
    const alphabet = @import("../alphabet.zig");

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const file = try tmp_dir.dir.createFile("test.fa", .{});
    try file.writeAll(">seq1\nACGTACGT\n");
    file.close();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try tmp_dir.dir.realpath("test.fa", &path_buf);

    // W=3, C=0: non-overlapping windows.
    var wr = try WindowReader.init(allocator, &alphabet.dna, path, 3, 0);
    defer wr.deinit();

    // Window 1: 3 residues (ACG).
    const w1 = (try wr.nextWindow()).?;
    var s1 = w1.seq;
    defer s1.deinit();
    const t1 = try s1.toText();
    defer allocator.free(t1);
    try std.testing.expectEqualStrings("ACG", t1);
    try std.testing.expectEqual(@as(i64, 1), w1.window_start);
    try std.testing.expectEqual(@as(i64, 3), w1.window_end);
    try std.testing.expect(!w1.complete);

    // Window 2: 3 residues (TAC).
    const w2 = (try wr.nextWindow()).?;
    var s2 = w2.seq;
    defer s2.deinit();
    const t2 = try s2.toText();
    defer allocator.free(t2);
    try std.testing.expectEqualStrings("TAC", t2);
    try std.testing.expectEqual(@as(i64, 4), w2.window_start);
    try std.testing.expectEqual(@as(i64, 6), w2.window_end);
    try std.testing.expect(!w2.complete);

    // Window 3: 2 residues (GT) -- partial last window.
    const w3 = (try wr.nextWindow()).?;
    var s3 = w3.seq;
    defer s3.deinit();
    const t3 = try s3.toText();
    defer allocator.free(t3);
    try std.testing.expectEqualStrings("GT", t3);
    try std.testing.expectEqual(@as(i64, 7), w3.window_start);
    try std.testing.expectEqual(@as(i64, 8), w3.window_end);
    try std.testing.expect(w3.complete);

    try std.testing.expectEqual(@as(?WindowResult, null), try wr.nextWindow());
}

test "WindowReader: multiple sequences with windows" {
    const allocator = std.testing.allocator;
    const alphabet = @import("../alphabet.zig");

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const file = try tmp_dir.dir.createFile("test.fa", .{});
    try file.writeAll(">chr1\nACGTACGT\n>chr2\nTGCATGCA\n");
    file.close();

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try tmp_dir.dir.realpath("test.fa", &path_buf);

    // W=5, C=2.
    var wr = try WindowReader.init(allocator, &alphabet.dna, path, 5, 2);
    defer wr.deinit();

    // chr1 window 1: W+C = 7 residues (ACGTACG).
    const w1 = (try wr.nextWindow()).?;
    var s1 = w1.seq;
    defer s1.deinit();
    try std.testing.expectEqualStrings("chr1", s1.name);
    try std.testing.expectEqual(@as(usize, 7), s1.len());
    try std.testing.expect(!w1.complete);

    // chr1 window 2: 2 context + 1 new residue (T) = 3 residues (CGT).
    const w2 = (try wr.nextWindow()).?;
    var s2 = w2.seq;
    defer s2.deinit();
    try std.testing.expectEqualStrings("chr1", s2.name);
    try std.testing.expect(w2.complete);

    const t2 = try s2.toText();
    defer allocator.free(t2);
    try std.testing.expectEqualStrings("CGT", t2);

    // chr2 window 1: W+C = 7 residues (TGCATGC).
    const w3 = (try wr.nextWindow()).?;
    var s3 = w3.seq;
    defer s3.deinit();
    try std.testing.expectEqualStrings("chr2", s3.name);
    try std.testing.expectEqual(@as(usize, 7), s3.len());
    try std.testing.expect(!w3.complete);

    // chr2 window 2: 2 context + 1 new (A) = 3 residues (GCA).
    const w4 = (try wr.nextWindow()).?;
    var s4 = w4.seq;
    defer s4.deinit();
    try std.testing.expectEqualStrings("chr2", s4.name);
    try std.testing.expect(w4.complete);

    const t4 = try s4.toText();
    defer allocator.free(t4);
    try std.testing.expectEqualStrings("GCA", t4);

    try std.testing.expectEqual(@as(?WindowResult, null), try wr.nextWindow());
}
