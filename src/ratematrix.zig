// Evolutionary rate matrix Q and probability matrix P = exp(tQ).
//
// A rate matrix Q is a square matrix where:
//   - Off-diagonal Q[i][j] >= 0 (rate of change from i to j)
//   - Rows sum to zero: Q[i][i] = -sum_{j!=i} Q[i][j]
//
// The transition probability matrix P(t) = exp(tQ) gives the
// probability of changing from residue i to j in time t.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Matrix = @import("matrix.zig").Matrix;

/// Create a rate matrix Q from an exchangeability matrix S and
/// stationary frequencies pi. Q[i][j] = S[i][j] * pi[j] for i != j,
/// and Q[i][i] set so rows sum to zero.
/// S must be symmetric. Returns a new n*n Matrix.
pub fn fromExchangeability(allocator: Allocator, s: Matrix, pi: []const f64) !Matrix {
    const n = s.rows;
    std.debug.assert(s.cols == n);
    std.debug.assert(pi.len == n);

    var q = try Matrix.init(allocator, n, n);
    errdefer q.deinit();

    for (0..n) |i| {
        var row_sum: f64 = 0;
        for (0..n) |j| {
            if (i != j) {
                const rate = s.get(i, j) * pi[j];
                q.set(i, j, rate);
                row_sum += rate;
            }
        }
        q.set(i, i, -row_sum);
    }

    return q;
}

/// Normalize a rate matrix Q so that the expected number of
/// substitutions per unit time is 1: -sum_i pi[i] * Q[i][i] = 1.
pub fn normalize(q: *Matrix, pi: []const f64) void {
    const n = q.rows;
    var rate: f64 = 0;
    for (0..n) |i| {
        rate -= pi[i] * q.get(i, i);
    }
    if (rate > 0) {
        q.scale(1.0 / rate);
    }
}

/// Compute the probability matrix P(t) = exp(tQ).
/// Caller owns the returned Matrix.
pub fn probMatrix(_: Allocator, q: Matrix, t: f64) !Matrix {
    return q.exp(t);
}

/// Validate a rate matrix: rows should sum to ~0.
pub fn isValid(q: Matrix, tol: f64) bool {
    if (q.rows != q.cols) return false;
    const n = q.rows;
    for (0..n) |i| {
        var row_sum: f64 = 0;
        for (0..n) |j| {
            row_sum += q.get(i, j);
        }
        if (@abs(row_sum) > tol) return false;
        // Diagonal should be negative
        if (q.get(i, i) > tol) return false;
        // Off-diagonals should be non-negative
        for (0..n) |j| {
            if (i != j and q.get(i, j) < -tol) return false;
        }
    }
    return true;
}

// --- Tests ---

test "fromExchangeability: 2x2" {
    const allocator = std.testing.allocator;
    // S = [[0, 1], [1, 0]], pi = [0.6, 0.4]
    // Q[0][1] = 1 * 0.4 = 0.4, Q[1][0] = 1 * 0.6 = 0.6
    // Q[0][0] = -0.4, Q[1][1] = -0.6
    var s = try Matrix.init(allocator, 2, 2);
    defer s.deinit();
    s.set(0, 1, 1.0);
    s.set(1, 0, 1.0);

    const pi = [_]f64{ 0.6, 0.4 };
    var q = try fromExchangeability(allocator, s, &pi);
    defer q.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, -0.4), q.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), q.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), q.get(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, -0.6), q.get(1, 1), 1e-10);
    try std.testing.expect(isValid(q, 1e-10));
}

test "normalize: expected rate becomes 1" {
    const allocator = std.testing.allocator;
    var s = try Matrix.init(allocator, 2, 2);
    defer s.deinit();
    s.set(0, 1, 2.0);
    s.set(1, 0, 2.0);

    const pi = [_]f64{ 0.5, 0.5 };
    var q = try fromExchangeability(allocator, s, &pi);
    defer q.deinit();

    normalize(&q, &pi);

    // Check: -sum(pi[i] * Q[i][i]) should be 1.0
    var rate: f64 = 0;
    for (0..2) |i| rate -= pi[i] * q.get(i, i);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rate, 1e-10);
}

test "probMatrix: P(0) = I" {
    const allocator = std.testing.allocator;
    var q = try Matrix.init(allocator, 2, 2);
    defer q.deinit();
    q.set(0, 0, -1.0);
    q.set(0, 1, 1.0);
    q.set(1, 0, 1.0);
    q.set(1, 1, -1.0);

    var p = try probMatrix(allocator, q, 0.0);
    defer p.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, 1.0), p.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), p.get(0, 1), 1e-10);
}

test "probMatrix: rows sum to 1" {
    const allocator = std.testing.allocator;
    var q = try Matrix.init(allocator, 3, 3);
    defer q.deinit();
    q.set(0, 0, -2.0);
    q.set(0, 1, 1.0);
    q.set(0, 2, 1.0);
    q.set(1, 0, 1.0);
    q.set(1, 1, -2.0);
    q.set(1, 2, 1.0);
    q.set(2, 0, 1.0);
    q.set(2, 1, 1.0);
    q.set(2, 2, -2.0);

    var p = try probMatrix(allocator, q, 0.5);
    defer p.deinit();

    for (0..3) |i| {
        var row_sum: f64 = 0;
        for (0..3) |j| row_sum += p.get(i, j);
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), row_sum, 1e-8);
    }
}
