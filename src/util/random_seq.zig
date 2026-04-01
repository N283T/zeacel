// Random biological sequence generation and in-place shuffling.

const std = @import("std");
const Allocator = std.mem.Allocator;
const alphabet_mod = @import("../alphabet.zig");
const Alphabet = alphabet_mod.Alphabet;
const Sequence = @import("../sequence.zig").Sequence;
const Random = @import("random.zig").Random;

/// Generate a random i.i.d. sequence from a residue frequency distribution.
///
/// freq[0..K-1] gives the probability of each canonical residue.
/// The caller owns the returned Sequence and must call deinit() on it.
pub fn randomIid(
    allocator: Allocator,
    rng: *Random,
    abc: *const Alphabet,
    name: []const u8,
    length: usize,
    freq: []const f64,
) !Sequence {
    const dsq = try allocator.alloc(u8, length);
    errdefer allocator.free(dsq);

    for (0..length) |i| {
        dsq[i] = @intCast(rng.choose(freq));
    }

    const name_copy = try allocator.dupe(u8, name);
    errdefer allocator.free(name_copy);

    return Sequence{
        .name = name_copy,
        .accession = null,
        .description = null,
        .taxonomy_id = null,
        .dsq = dsq,
        .secondary_structure = null,
        .source = null,
        .abc = abc,
        .allocator = allocator,
    };
}

/// Generate a random sequence with uniform residue frequencies over K canonical residues.
///
/// The caller owns the returned Sequence and must call deinit() on it.
pub fn randomUniform(
    allocator: Allocator,
    rng: *Random,
    abc: *const Alphabet,
    name: []const u8,
    length: usize,
) !Sequence {
    const k = abc.k;
    const freq_value = 1.0 / @as(f64, @floatFromInt(k));

    // Build uniform frequency array on the stack for small K, heap for large K.
    // Amino acids have K=20 — always fits on the stack.
    var freq_buf: [32]f64 = undefined;
    if (k > freq_buf.len) return error.AlphabetTooLarge;
    for (0..k) |i| freq_buf[i] = freq_value;

    return randomIid(allocator, rng, abc, name, length, freq_buf[0..k]);
}

/// Shuffle a digital sequence in place using the Fisher-Yates algorithm.
/// Preserves mono-residue composition.
pub fn shuffle(rng: *Random, dsq: []u8) void {
    var i = dsq.len;
    while (i > 1) {
        i -= 1;
        const j = rng.uniformInt(@intCast(i + 1));
        const tmp = dsq[i];
        dsq[i] = dsq[j];
        dsq[j] = tmp;
    }
}

/// Shuffle preserving mono-residue composition.
/// Alias for shuffle — included for API clarity matching Easel naming.
pub const shuffleMono = shuffle;

/// Maximum alphabet size supported by the dinucleotide shuffle.
/// Amino acid alphabets have K=20; this leaves headroom.
const MAX_K: usize = 32;

/// Shuffle preserving dinucleotide composition (Altschul-Erickson algorithm).
///
/// Maintains the count of each dinucleotide pair while randomizing the
/// sequence.  `k` is the alphabet size (number of canonical residues).
///
/// For sequences of length 0 or 1, the sequence is unchanged.
/// For length 2, the sequence is unchanged (only one dinucleotide).
///
/// Algorithm (Altschul & Erickson 1985; Easel esl_rsq_XShuffleDP):
///   1. Count directed-edge frequencies E[i][j] from consecutive pairs.
///   2. For each node with outgoing edges, randomly mark one edge as
///      "last" (reserve it so an Eulerian path is guaranteed).
///   3. Traverse from dsq[0]: at each node, pick a random non-last edge;
///      when only the last edge remains, use it.  Append the destination
///      residue and continue until all edges are consumed.
pub fn shuffleDi(rng: *Random, dsq: []u8, k: u8) void {
    const K: usize = @intCast(k);
    if (K > MAX_K) return; // safety guard
    if (dsq.len <= 2) return; // nothing to shuffle

    // --- Step 1: count edge frequencies E[i][j] ---
    var E: [MAX_K][MAX_K]u16 = .{.{0} ** MAX_K} ** MAX_K;
    for (0..dsq.len - 1) |pos| {
        E[dsq[pos]][dsq[pos + 1]] += 1;
    }

    // --- Step 2: for each node, randomly choose one outgoing edge to be "last" ---
    // iE[i] = the destination j of the reserved last edge from node i.
    // We decrement E[i][j] for that edge so it won't be picked during normal traversal.
    var iE: [MAX_K]u8 = .{0} ** MAX_K;
    for (0..K) |i| {
        // Count total outgoing edges from node i.
        var total: u32 = 0;
        for (0..K) |j| total += E[i][j];
        if (total == 0) continue;

        // Pick a random edge uniformly among all outgoing edges.
        var pick: u32 = rng.uniformInt(total);
        for (0..K) |j| {
            if (E[i][j] == 0) continue;
            if (pick < E[i][j]) {
                iE[i] = @intCast(j);
                E[i][j] -= 1;
                break;
            }
            pick -= E[i][j];
        }
    }

    // --- Step 3: traverse Eulerian path ---
    var pos: usize = 0;
    dsq[pos] = dsq[0]; // first residue stays the same
    var x: usize = dsq[0];

    while (pos < dsq.len - 2) {
        // Count remaining (non-last) outgoing edges from x.
        var total: u32 = 0;
        for (0..K) |j| total += E[x][j];

        if (total == 0) {
            // Only the reserved last edge remains; use it.
            pos += 1;
            dsq[pos] = iE[x];
            x = iE[x];
        } else {
            // Pick a random edge from the remaining ones.
            var pick: u32 = rng.uniformInt(total);
            for (0..K) |j| {
                if (E[x][j] == 0) continue;
                if (pick < E[x][j]) {
                    pos += 1;
                    dsq[pos] = @intCast(j);
                    E[x][j] -= 1;
                    x = j;
                    break;
                }
                pick -= E[x][j];
            }
        }
    }

    // The last residue must match the original last residue.
    // With a correct Eulerian path it is guaranteed by the reserved last edges,
    // but we set it explicitly for clarity.  The final "last edge" from the
    // penultimate node leads to the original last residue.
    dsq[dsq.len - 1] = iE[x];
}

/// Generate a random sequence preserving 0th-order (mono-residue) composition.
///
/// Counts residue frequencies in `dsq` (only codes 0..K-1), normalizes to
/// probabilities, then generates a new sequence by sampling from that distribution.
/// The caller owns the returned slice and must free it with `allocator`.
pub fn markov0(allocator: Allocator, rng: *Random, dsq: []const u8, K: usize) ![]u8 {
    if (dsq.len == 0) return try allocator.alloc(u8, 0);

    // Count residue frequencies.
    var counts: [MAX_K]u64 = .{0} ** MAX_K;
    var total: u64 = 0;
    for (dsq) |code| {
        if (code < K) {
            counts[code] += 1;
            total += 1;
        }
    }

    // Normalize to probabilities.
    var freq: [MAX_K]f64 = .{0.0} ** MAX_K;
    if (total > 0) {
        for (0..K) |i| {
            freq[i] = @as(f64, @floatFromInt(counts[i])) / @as(f64, @floatFromInt(total));
        }
    } else {
        // All residues were out-of-range; use uniform.
        const uniform_p = 1.0 / @as(f64, @floatFromInt(K));
        for (0..K) |i| freq[i] = uniform_p;
    }

    // Generate new sequence.
    const result = try allocator.alloc(u8, dsq.len);
    for (0..dsq.len) |i| {
        result[i] = @intCast(rng.choose(freq[0..K]));
    }
    return result;
}

/// Generate a random sequence preserving 1st-order (di-residue) Markov transitions.
///
/// Counts the transition matrix T[a][b] from consecutive pairs in `dsq`,
/// normalizes each row to transition probabilities, then generates a new
/// sequence: the first residue is sampled from the initial distribution,
/// each subsequent residue is sampled from P[prev][.].
/// The caller owns the returned slice and must free it with `allocator`.
pub fn markov1(allocator: Allocator, rng: *Random, dsq: []const u8, K: usize) ![]u8 {
    if (dsq.len == 0) return try allocator.alloc(u8, 0);

    // Count initial residue frequencies and transition matrix.
    var f_counts: [MAX_K]u64 = .{0} ** MAX_K;
    var f_total: u64 = 0;
    var T: [MAX_K][MAX_K]u64 = .{.{0} ** MAX_K} ** MAX_K;

    for (dsq) |code| {
        if (code < K) {
            f_counts[code] += 1;
            f_total += 1;
        }
    }

    for (0..dsq.len - 1) |pos| {
        const a = dsq[pos];
        const b = dsq[pos + 1];
        if (a < K and b < K) {
            T[a][b] += 1;
        }
    }

    // Normalize initial distribution.
    var f: [MAX_K]f64 = .{0.0} ** MAX_K;
    if (f_total > 0) {
        for (0..K) |i| {
            f[i] = @as(f64, @floatFromInt(f_counts[i])) / @as(f64, @floatFromInt(f_total));
        }
    } else {
        const uniform_p = 1.0 / @as(f64, @floatFromInt(K));
        for (0..K) |i| f[i] = uniform_p;
    }

    // Normalize each row of T to transition probabilities P[a][.].
    var P: [MAX_K][MAX_K]f64 = .{.{0.0} ** MAX_K} ** MAX_K;
    for (0..K) |a| {
        var row_total: u64 = 0;
        for (0..K) |b| row_total += T[a][b];
        if (row_total > 0) {
            for (0..K) |b| {
                P[a][b] = @as(f64, @floatFromInt(T[a][b])) / @as(f64, @floatFromInt(row_total));
            }
        } else {
            // No transitions observed from this residue; fall back to initial distribution.
            for (0..K) |b| P[a][b] = f[b];
        }
    }

    // Generate sequence.
    const result = try allocator.alloc(u8, dsq.len);

    // First residue: sample from initial distribution.
    result[0] = @intCast(rng.choose(f[0..K]));

    // Subsequent residues: sample from transition row of previous residue.
    for (1..dsq.len) |i| {
        const prev: usize = result[i - 1];
        result[i] = @intCast(rng.choose(P[prev][0..K]));
    }

    return result;
}

// --- Tests ---

test "randomIid: length and all codes < K for DNA" {
    const allocator = std.testing.allocator;
    var rng = Random.init(1);
    const freq = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    var seq = try randomIid(allocator, &rng, &alphabet_mod.dna, "test", 100, &freq);
    defer seq.deinit();

    try std.testing.expectEqual(@as(usize, 100), seq.len());
    for (seq.dsq) |code| {
        try std.testing.expect(code < alphabet_mod.dna.k);
    }
}

test "randomIid: zero-length sequence" {
    const allocator = std.testing.allocator;
    var rng = Random.init(2);
    const freq = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    var seq = try randomIid(allocator, &rng, &alphabet_mod.dna, "empty", 0, &freq);
    defer seq.deinit();

    try std.testing.expectEqual(@as(usize, 0), seq.len());
}

test "randomUniform: amino acid — correct length and canonical codes" {
    const allocator = std.testing.allocator;
    var rng = Random.init(3);
    var seq = try randomUniform(allocator, &rng, &alphabet_mod.amino, "prot", 200);
    defer seq.deinit();

    try std.testing.expectEqual(@as(usize, 200), seq.len());
    for (seq.dsq) |code| {
        try std.testing.expect(code < alphabet_mod.amino.k);
    }
}

test "randomUniform: DNA has correct length and canonical codes" {
    const allocator = std.testing.allocator;
    var rng = Random.init(4);
    var seq = try randomUniform(allocator, &rng, &alphabet_mod.dna, "dna_seq", 50);
    defer seq.deinit();

    try std.testing.expectEqual(@as(usize, 50), seq.len());
    for (seq.dsq) |code| {
        try std.testing.expect(code < alphabet_mod.dna.k);
    }
}

test "shuffle: preserves mono-residue composition" {
    var rng = Random.init(5);
    var dsq = [_]u8{ 0, 0, 1, 1, 2, 2, 3, 3 };

    var counts_before = [_]usize{0} ** 4;
    for (dsq) |c| counts_before[c] += 1;

    shuffle(&rng, &dsq);

    var counts_after = [_]usize{0} ** 4;
    for (dsq) |c| counts_after[c] += 1;

    try std.testing.expectEqualSlices(usize, &counts_before, &counts_after);
}

test "shuffle: length preserved" {
    var rng = Random.init(6);
    var dsq = [_]u8{ 0, 1, 2, 3, 0, 1 };
    const original_len = dsq.len;
    shuffle(&rng, &dsq);
    try std.testing.expectEqual(original_len, dsq.len);
}

test "shuffle: empty slice does not panic" {
    var rng = Random.init(7);
    var dsq = [_]u8{};
    shuffle(&rng, &dsq);
}

test "shuffle: single element unchanged" {
    var rng = Random.init(8);
    var dsq = [_]u8{3};
    shuffle(&rng, &dsq);
    try std.testing.expectEqual(@as(u8, 3), dsq[0]);
}

test "reproducibility: same seed produces same random sequence" {
    const allocator = std.testing.allocator;

    var rng1 = Random.init(42);
    var seq1 = try randomUniform(allocator, &rng1, &alphabet_mod.dna, "s1", 50);
    defer seq1.deinit();

    var rng2 = Random.init(42);
    var seq2 = try randomUniform(allocator, &rng2, &alphabet_mod.dna, "s1", 50);
    defer seq2.deinit();

    try std.testing.expectEqualSlices(u8, seq1.dsq, seq2.dsq);
}

test "reproducibility: different seeds produce different sequences" {
    const allocator = std.testing.allocator;

    var rng1 = Random.init(100);
    var seq1 = try randomUniform(allocator, &rng1, &alphabet_mod.dna, "s1", 50);
    defer seq1.deinit();

    var rng2 = Random.init(200);
    var seq2 = try randomUniform(allocator, &rng2, &alphabet_mod.dna, "s2", 50);
    defer seq2.deinit();

    // Extremely unlikely to match by chance with length 50.
    try std.testing.expect(!std.mem.eql(u8, seq1.dsq, seq2.dsq));
}

test "markov0: preserves mono-residue composition" {
    const allocator = std.testing.allocator;
    var rng = Random.init(12345);

    // Create a DNA sequence with known biased composition: heavy on A (0) and T (3).
    var input: [2000]u8 = undefined;
    for (0..2000) |i| {
        input[i] = switch (i % 10) {
            0, 1, 2, 3 => 0, // 40% A
            4, 5, 6 => 1, // 30% C
            7, 8 => 2, // 20% G
            9 => 3, // 10% T
            else => unreachable,
        };
    }

    const result = try markov0(allocator, &rng, &input, 4);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 2000), result.len);

    // Count composition in result.
    var counts = [_]usize{0} ** 4;
    for (result) |code| {
        try std.testing.expect(code < 4);
        counts[code] += 1;
    }

    // Expected: ~800 A, ~600 C, ~400 G, ~200 T. Tolerate +-80 (~4%).
    try std.testing.expect(counts[0] > 700 and counts[0] < 900);
    try std.testing.expect(counts[1] > 500 and counts[1] < 700);
    try std.testing.expect(counts[2] > 300 and counts[2] < 500);
    try std.testing.expect(counts[3] > 120 and counts[3] < 280);
}

test "markov0: empty sequence returns empty" {
    const allocator = std.testing.allocator;
    var rng = Random.init(99);
    const empty: []const u8 = &.{};
    const result = try markov0(allocator, &rng, empty, 4);
    defer allocator.free(result);
    try std.testing.expectEqual(@as(usize, 0), result.len);
}

test "markov0: single residue" {
    const allocator = std.testing.allocator;
    var rng = Random.init(77);
    const input = [_]u8{2};
    const result = try markov0(allocator, &rng, &input, 4);
    defer allocator.free(result);
    try std.testing.expectEqual(@as(usize, 1), result.len);
    // Only residue 2 was observed, so output must be 2.
    try std.testing.expectEqual(@as(u8, 2), result[0]);
}

test "markov1: preserves di-residue transitions" {
    const allocator = std.testing.allocator;
    var rng = Random.init(54321);

    // Create a sequence with a strong pattern: alternating 0-1-0-1-...
    // This creates strong transitions: 0->1 and 1->0.
    var input: [2000]u8 = undefined;
    for (0..2000) |i| {
        input[i] = @intCast(i % 2);
    }

    const result = try markov1(allocator, &rng, &input, 4);
    defer allocator.free(result);

    try std.testing.expectEqual(@as(usize, 2000), result.len);

    // Count transitions in the result.
    var T_result = [_][4]usize{.{0} ** 4} ** 4;
    for (0..result.len - 1) |pos| {
        T_result[result[pos]][result[pos + 1]] += 1;
    }

    // The input has transitions: 0->1 (999), 1->0 (1000).
    // The result should overwhelmingly consist of 0->1 and 1->0 transitions.
    const dominant = T_result[0][1] + T_result[1][0];
    try std.testing.expect(dominant > 1900);
}

test "markov1: empty sequence returns empty" {
    const allocator = std.testing.allocator;
    var rng = Random.init(88);
    const empty: []const u8 = &.{};
    const result = try markov1(allocator, &rng, empty, 4);
    defer allocator.free(result);
    try std.testing.expectEqual(@as(usize, 0), result.len);
}

test "markov1: single residue returns single residue" {
    const allocator = std.testing.allocator;
    var rng = Random.init(66);
    const input = [_]u8{3};
    const result = try markov1(allocator, &rng, &input, 4);
    defer allocator.free(result);
    try std.testing.expectEqual(@as(usize, 1), result.len);
    // Only residue 3 observed, so output must be 3.
    try std.testing.expectEqual(@as(u8, 3), result[0]);
}

test "markov1: skips out-of-range codes" {
    const allocator = std.testing.allocator;
    var rng = Random.init(55);
    // Sequence with some codes >= K (treated as gaps/degenerate).
    const input = [_]u8{ 0, 1, 255, 0, 1, 0, 1, 254, 0, 1 };
    const result = try markov1(allocator, &rng, &input, 4);
    defer allocator.free(result);
    try std.testing.expectEqual(@as(usize, 10), result.len);
    // All output codes must be < K.
    for (result) |code| {
        try std.testing.expect(code < 4);
    }
}
