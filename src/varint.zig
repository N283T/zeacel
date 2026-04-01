// Variable-length integer encoding/decoding.
//
// Supports Google varint (base-128) encoding used in compressed
// data formats. Each byte uses 7 bits for data and 1 bit (MSB)
// as a continuation flag.

const std = @import("std");

/// Encode a u64 value as a Google varint (base-128, little-endian).
/// Returns the number of bytes written, or error.BufferTooSmall if
/// the buffer is not large enough (up to 10 bytes needed for u64 max).
pub fn encode(value: u64, buf: []u8) error{BufferTooSmall}!usize {
    const required = encodedSize(value);
    if (buf.len < required) return error.BufferTooSmall;

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

// --- Exponential-Golomb coding ---

/// Encode a non-negative integer using exponential-Golomb-k code.
/// Returns the codeword (right-flushed in the u64) and its bit length.
/// k=0 is standard exponential-Golomb; k>0 are generalized codes.
pub fn expGolombEncode(v: u64, k: u5) !struct { code: u64, bits: u6 } {
    const m: u64 = @as(u64, 1) << k;
    const q: u64 = v / m;
    const r: u64 = v % m;

    // Determine bit width of (q+1)
    const qp1 = q + 1;
    var a: u6 = 0;
    var tmp = qp1;
    while (tmp > 0) {
        a += 1;
        tmp >>= 1;
    }
    // Exp-Golomb-0 part: (a-1) leading zeros + a bits of q+1 = 2*a - 1 bits
    const n_wide: u7 = @as(u7, 2) * @as(u7, a) - 1 + @as(u7, k);
    if (n_wide > 64) return error.Overflow;
    const n: u6 = @intCast(n_wide);

    var code: u64 = qp1;
    code = (code << k) | r;
    return .{ .code = code, .bits = n };
}

/// Decode an exponential-Golomb-k codeword from the high-order (left-flushed) bits of `code`.
/// Returns the decoded value and the number of bits consumed.
pub fn expGolombDecode(code: u64, k: u5) !struct { value: u64, bits: u6 } {
    if (code == 0) return error.UnexpectedEndOfInput;

    // Find position of leftmost 1 bit
    var b: u6 = 63;
    while (b > 0 and (code & (@as(u64, 1) << b)) == 0) {
        b -= 1;
    }
    if ((code & (@as(u64, 1) << b)) == 0) return error.UnexpectedEndOfInput;

    // b is position of MSB (63..0). Number of leading zeros = 63 - b.
    // Codeword length: 2*(64-b) + k - 1
    const n_wide: u7 = @as(u7, 2) * (@as(u7, 64) - @as(u7, b)) + @as(u7, k) - 1;
    if (n_wide > 64) return error.Overflow;
    const n: u6 = @intCast(n_wide);

    const shift_amt: u7 = @as(u7, 64) - @as(u7, n);
    var shifted = code >> @as(u6, @intCast(shift_amt));
    const r = shifted & ((@as(u64, 1) << k) - 1);
    shifted >>= k;
    const v = (shifted - 1) * (@as(u64, 1) << k) + r;

    return .{ .value = v, .bits = n };
}

// --- Golomb-Rice coding ---

/// Encode a non-negative integer using Golomb-Rice-k code.
/// Returns the codeword (right-flushed) and its bit length.
/// The code is: q ones, a zero, then k bits of remainder.
pub fn riceEncode(v: u64, k: u5) !struct { code: u64, bits: u6 } {
    const m: u64 = @as(u64, 1) << k;
    const q: u64 = v / m;
    const r: u64 = v % m;

    const n_wide: u7 = @as(u7, @intCast(q)) + @as(u7, k) + 1;
    if (n_wide > 64) return error.Overflow;
    const n: u6 = @intCast(n_wide);

    // q ones followed by a 0 followed by k bits of r
    var code: u64 = if (q > 0) (@as(u64, 1) << @intCast(q)) - 1 else 0;
    code = (code << (1 + @as(u6, k))) | r;

    return .{ .code = code, .bits = n };
}

/// Decode a Golomb-Rice-k codeword from the high-order (left-flushed) bits of `code`.
/// Returns the decoded value and the number of bits consumed.
pub fn riceDecode(code: u64, k: u5) !struct { value: u64, bits: u6 } {
    // A valid Rice code must have at least one 0 bit
    if (~code == 0) return error.UnexpectedEndOfInput;

    // Count leading 1 bits = quotient
    var q: u6 = 0;
    while (@as(u7, q) <= 64 - (@as(u7, k) + 1)) {
        if ((~code) & (@as(u64, 1) << (@as(u6, 63) - q)) != 0) break;
        q += 1;
    }

    const n_wide: u7 = @as(u7, q) + @as(u7, k) + 1;
    if (n_wide > 64) return error.Overflow;
    const n: u6 = @intCast(n_wide);

    const shift_amt: u7 = @as(u7, 64) - @as(u7, n);
    const shifted = code >> @as(u6, @intCast(shift_amt));
    const r = shifted & ((@as(u64, 1) << k) - 1);
    const v = @as(u64, q) * (@as(u64, 1) << k) + r;

    return .{ .value = v, .bits = n };
}

// --- Tests ---

test "encode/decode: zero" {
    var buf: [10]u8 = undefined;
    const n = try encode(0, &buf);
    try std.testing.expectEqual(@as(usize, 1), n);
    try std.testing.expectEqual(@as(u8, 0), buf[0]);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(@as(u64, 0), result.value);
    try std.testing.expectEqual(@as(usize, 1), result.bytes_read);
}

test "encode/decode: small value (127)" {
    var buf: [10]u8 = undefined;
    const n = try encode(127, &buf);
    try std.testing.expectEqual(@as(usize, 1), n);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(@as(u64, 127), result.value);
}

test "encode/decode: 128 (needs 2 bytes)" {
    var buf: [10]u8 = undefined;
    const n = try encode(128, &buf);
    try std.testing.expectEqual(@as(usize, 2), n);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(@as(u64, 128), result.value);
}

test "encode/decode: large value" {
    var buf: [10]u8 = undefined;
    const val: u64 = 123456789;
    const n = try encode(val, &buf);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(val, result.value);
}

test "encode/decode: max u64" {
    var buf: [10]u8 = undefined;
    const val: u64 = std.math.maxInt(u64);
    const n = try encode(val, &buf);

    const result = try decode(buf[0..n]);
    try std.testing.expectEqual(val, result.value);
}

test "encode: buffer too small returns error" {
    // u64 max needs 10 bytes; a 9-byte buffer should fail
    var buf: [9]u8 = undefined;
    try std.testing.expectError(error.BufferTooSmall, encode(std.math.maxInt(u64), &buf));

    // 128 needs 2 bytes; a 1-byte buffer should fail
    var buf1: [1]u8 = undefined;
    try std.testing.expectError(error.BufferTooSmall, encode(128, &buf1));

    // 0 needs 1 byte; a 0-length buffer should fail
    var buf0: [0]u8 = undefined;
    try std.testing.expectError(error.BufferTooSmall, encode(0, &buf0));

    // 127 fits in 1 byte; a 1-byte buffer should succeed
    var buf_ok: [1]u8 = undefined;
    const n = try encode(127, &buf_ok);
    try std.testing.expectEqual(@as(usize, 1), n);
}

test "encodedSize" {
    try std.testing.expectEqual(@as(usize, 1), encodedSize(0));
    try std.testing.expectEqual(@as(usize, 1), encodedSize(127));
    try std.testing.expectEqual(@as(usize, 2), encodedSize(128));
    try std.testing.expectEqual(@as(usize, 5), encodedSize(0xFFFFFFFF));
}

/// Helper: left-flush a right-flushed codeword of `bits` length to MSB position in u64.
fn leftFlush(code: u64, bits: u6) u64 {
    if (bits == 0) return 0;
    const shift: u7 = @as(u7, 64) - @as(u7, bits);
    return code << @as(u6, @intCast(shift));
}

test "expGolomb: encode/decode k=0 roundtrip" {
    const test_values = [_]u64{ 0, 1, 2, 3, 4, 10, 100, 1000 };
    for (test_values) |v| {
        const enc = try expGolombEncode(v, 0);
        const left_flushed = leftFlush(enc.code, enc.bits);
        const dec = try expGolombDecode(left_flushed, 0);
        try std.testing.expectEqual(v, dec.value);
        try std.testing.expectEqual(enc.bits, dec.bits);
    }
}

test "expGolomb: encode/decode k=3 roundtrip" {
    const test_values = [_]u64{ 0, 1, 7, 8, 15, 16, 255 };
    for (test_values) |v| {
        const enc = try expGolombEncode(v, 3);
        const left_flushed = leftFlush(enc.code, enc.bits);
        const dec = try expGolombDecode(left_flushed, 3);
        try std.testing.expectEqual(v, dec.value);
        try std.testing.expectEqual(enc.bits, dec.bits);
    }
}

test "expGolomb: known code for v=0 k=0 is 1 (1 bit)" {
    const enc = try expGolombEncode(0, 0);
    try std.testing.expectEqual(@as(u64, 1), enc.code);
    try std.testing.expectEqual(@as(u6, 1), enc.bits);
}

test "expGolomb: known code for v=4 k=0 is 00101 (5 bits)" {
    const enc = try expGolombEncode(4, 0);
    try std.testing.expectEqual(@as(u64, 5), enc.code);
    try std.testing.expectEqual(@as(u6, 5), enc.bits);
}

test "rice: encode/decode k=2 roundtrip" {
    const test_values = [_]u64{ 0, 1, 2, 3, 4, 7, 8, 15, 20 };
    for (test_values) |v| {
        const enc = try riceEncode(v, 2);
        const left_flushed = leftFlush(enc.code, enc.bits);
        const dec = try riceDecode(left_flushed, 2);
        try std.testing.expectEqual(v, dec.value);
        try std.testing.expectEqual(enc.bits, dec.bits);
    }
}

test "rice: encode/decode k=0 roundtrip" {
    const test_values = [_]u64{ 0, 1, 2, 3, 10 };
    for (test_values) |v| {
        const enc = try riceEncode(v, 0);
        const left_flushed = leftFlush(enc.code, enc.bits);
        const dec = try riceDecode(left_flushed, 0);
        try std.testing.expectEqual(v, dec.value);
        try std.testing.expectEqual(enc.bits, dec.bits);
    }
}

test "rice: known code for v=0 k=2 is 000 (3 bits)" {
    const enc = try riceEncode(0, 2);
    try std.testing.expectEqual(@as(u64, 0), enc.code);
    try std.testing.expectEqual(@as(u6, 3), enc.bits);
}
