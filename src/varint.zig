// Variable-length integer encoding/decoding.
//
// Supports Google varint (base-128) encoding used in compressed
// data formats. Each byte uses 7 bits for data and 1 bit (MSB)
// as a continuation flag.

const std = @import("std");

/// Encode a u64 value as a Google varint (base-128, little-endian).
/// Returns the number of bytes written.
pub fn encode(value: u64, buf: []u8) usize {
    var v = value;
    var i: usize = 0;
    while (v >= 0x80) {
        buf[i] = @intCast((v & 0x7F) | 0x80);
        v >>= 7;
        i += 1;
    }
    buf[i] = @intCast(v & 0x7F);
    return i + 1;
}

/// Decode a Google varint from a byte buffer.
/// Returns the decoded value and the number of bytes consumed.
pub fn decode(buf: []const u8) !struct { value: u64, bytes_read: usize } {
    var value: u64 = 0;
    var shift: u6 = 0;
    for (buf, 0..) |byte, i| {
        value |= @as(u64, byte & 0x7F) << shift;
        if (byte & 0x80 == 0) {
            return .{ .value = value, .bytes_read = i + 1 };
        }
        shift +|= 7;
        if (shift > 63) return error.Overflow;
    }
    return error.UnexpectedEndOfInput;
}

/// Return the number of bytes needed to encode a value.
pub fn encodedSize(value: u64) usize {
    if (value == 0) return 1;
    var v = value;
    var n: usize = 0;
    while (v > 0) {
        n += 1;
        v >>= 7;
    }
    return n;
}

// --- Tests ---

test "encode/decode: zero" {
    var buf: [10]u8 = undefined;
    const n = encode(0, &buf);
    try std.testing.expectEqual(@as(usize, 1), n);
    try std.testing.expectEqual(@as(u8, 0), buf[0]);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(@as(u64, 0), result.value);
    try std.testing.expectEqual(@as(usize, 1), result.bytes_read);
}

test "encode/decode: small value (127)" {
    var buf: [10]u8 = undefined;
    const n = encode(127, &buf);
    try std.testing.expectEqual(@as(usize, 1), n);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(@as(u64, 127), result.value);
}

test "encode/decode: 128 (needs 2 bytes)" {
    var buf: [10]u8 = undefined;
    const n = encode(128, &buf);
    try std.testing.expectEqual(@as(usize, 2), n);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(@as(u64, 128), result.value);
}

test "encode/decode: large value" {
    var buf: [10]u8 = undefined;
    const val: u64 = 123456789;
    const n = encode(val, &buf);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(val, result.value);
}

test "encode/decode: max u64" {
    var buf: [10]u8 = undefined;
    const val: u64 = std.math.maxInt(u64);
    const n = encode(val, &buf);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(val, result.value);
}

test "encodedSize" {
    try std.testing.expectEqual(@as(usize, 1), encodedSize(0));
    try std.testing.expectEqual(@as(usize, 1), encodedSize(127));
    try std.testing.expectEqual(@as(usize, 2), encodedSize(128));
    try std.testing.expectEqual(@as(usize, 5), encodedSize(0xFFFFFFFF));
}
