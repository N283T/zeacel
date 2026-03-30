// Biological sequence alphabets (DNA, RNA, Amino Acid).

const std = @import("std");

pub const AlphabetType = enum {
    dna,
    rna,
    amino,
};

pub const INVALID_CODE: u8 = 0xFF;

pub const Alphabet = struct {
    kind: AlphabetType,
    k: u8, // Canonical residue count (4 for DNA/RNA, 20 for amino)
    kp: u8, // Total symbol count
    symbols: []const u8,
    encode_map: [128]u8, // ASCII -> digital code. INVALID_CODE for unmapped.
    complement: ?[]const u8, // Complement table (DNA/RNA only)

    pub fn encode(self: *const Alphabet, char: u8) error{InvalidCharacter}!u8 {
        if (char > 127) return error.InvalidCharacter;
        const code = self.encode_map[char];
        if (code == INVALID_CODE) return error.InvalidCharacter;
        return code;
    }

    pub fn decode(self: *const Alphabet, code: u8) u8 {
        return self.symbols[code];
    }

    pub fn isCanonical(self: *const Alphabet, code: u8) bool {
        return code < self.k;
    }

    pub fn isGap(self: *const Alphabet, code: u8) bool {
        return code == self.k;
    }

    pub fn isDegenerate(self: *const Alphabet, code: u8) bool {
        return code > self.k and code < self.kp - 2;
    }

    pub fn isUnknown(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 3;
    }

    pub fn isNonresidue(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 2;
    }

    pub fn isMissing(self: *const Alphabet, code: u8) bool {
        return code == self.kp - 1;
    }

    pub fn gapCode(self: *const Alphabet) u8 {
        return self.k;
    }

    pub fn unknownCode(self: *const Alphabet) u8 {
        return self.kp - 3;
    }
};

// --- Tests for Task 2 ---

test "Alphabet: struct size constants for DNA" {
    // Manually construct a minimal DNA-like alphabet for testing.
    // dna_symbols = "ACGT-RYMKSWHBVDN*~" has 18 chars: k=4, kp=18
    const symbols = "ACGT-RYMKSWHBVDN*~";
    var encode_map = [_]u8{INVALID_CODE} ** 128;
    for (symbols, 0..) |sym, i| {
        encode_map[sym] = @intCast(i);
    }

    const alpha = Alphabet{
        .kind = .dna,
        .k = 4,
        .kp = 18,
        .symbols = symbols,
        .encode_map = encode_map,
        .complement = null,
    };

    try std.testing.expectEqual(@as(u8, 4), alpha.k);
    try std.testing.expectEqual(@as(u8, 18), alpha.kp);
    try std.testing.expectEqual(@as(u8, 4), alpha.gapCode()); // gap at index k=4
    try std.testing.expectEqual(@as(u8, 15), alpha.unknownCode()); // kp-3 = 18-3 = 15
}

test "Alphabet: classification methods" {
    const symbols = "ACGT-RYMKSWHBVDN*~";
    var encode_map = [_]u8{INVALID_CODE} ** 128;
    for (symbols, 0..) |sym, i| {
        encode_map[sym] = @intCast(i);
    }

    const alpha = Alphabet{
        .kind = .dna,
        .k = 4,
        .kp = 18,
        .symbols = symbols,
        .encode_map = encode_map,
        .complement = null,
    };

    // Canonical: codes 0..3 (A, C, G, T)
    try std.testing.expect(alpha.isCanonical(0)); // A
    try std.testing.expect(alpha.isCanonical(3)); // T
    try std.testing.expect(!alpha.isCanonical(4)); // gap

    // Gap: code == k == 4
    try std.testing.expect(alpha.isGap(4));
    try std.testing.expect(!alpha.isGap(0));
    try std.testing.expect(!alpha.isGap(5));

    // Degenerate: codes > k and < kp-2, i.e., 5..14 (R,Y,M,K,S,W,H,B,V,D)
    try std.testing.expect(alpha.isDegenerate(5)); // R
    try std.testing.expect(alpha.isDegenerate(14)); // D
    try std.testing.expect(!alpha.isDegenerate(4)); // gap
    try std.testing.expect(!alpha.isDegenerate(15)); // N (unknown)
    try std.testing.expect(!alpha.isDegenerate(16)); // * (nonresidue)

    // Unknown: code == kp-3 == 15 (N)
    try std.testing.expect(alpha.isUnknown(15)); // N
    try std.testing.expect(!alpha.isUnknown(14));
    try std.testing.expect(!alpha.isUnknown(16));

    // Nonresidue: code == kp-2 == 16 (*)
    try std.testing.expect(alpha.isNonresidue(16)); // *
    try std.testing.expect(!alpha.isNonresidue(15));
    try std.testing.expect(!alpha.isNonresidue(17));

    // Missing: code == kp-1 == 17 (~)
    try std.testing.expect(alpha.isMissing(17)); // ~
    try std.testing.expect(!alpha.isMissing(16));
    try std.testing.expect(!alpha.isMissing(0));
}
