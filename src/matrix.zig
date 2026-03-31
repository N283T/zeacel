// Dense matrix of f64 values with row-major storage.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const Matrix = struct {
    data: []f64,
    rows: usize,
    cols: usize,
    allocator: Allocator,

    /// Allocate a rows×cols matrix initialised to zero.
    pub fn init(allocator: Allocator, rows: usize, cols: usize) !Matrix {
        const data = try allocator.alloc(f64, rows * cols);
        @memset(data, 0.0);
        return Matrix{ .data = data, .rows = rows, .cols = cols, .allocator = allocator };
    }

    /// Allocate an n×n identity matrix.
    pub fn initIdentity(allocator: Allocator, n: usize) !Matrix {
        var m = try Matrix.init(allocator, n, n);
        for (0..n) |i| m.data[i * n + i] = 1.0;
        return m;
    }

    pub fn deinit(self: *Matrix) void {
        self.allocator.free(self.data);
    }

    pub fn get(self: Matrix, i: usize, j: usize) f64 {
        return self.data[i * self.cols + j];
    }

    pub fn set(self: *Matrix, i: usize, j: usize, val: f64) void {
        self.data[i * self.cols + j] = val;
    }

    /// Return a new matrix equal to self + other. Dimensions must match.
    pub fn add(self: Matrix, other: Matrix) !Matrix {
        if (self.rows != other.rows or self.cols != other.cols) return error.DimensionMismatch;
        const result = try Matrix.init(self.allocator, self.rows, self.cols);
        for (result.data, self.data, other.data) |*r, a, b| r.* = a + b;
        return result;
    }

    /// Scale all elements in-place by factor.
    pub fn scale(self: *Matrix, factor: f64) void {
        for (self.data) |*v| v.* *= factor;
    }

    /// Return a new matrix C = A * B.
    pub fn multiply(allocator: Allocator, a: Matrix, b: Matrix) !Matrix {
        if (a.cols != b.rows) return error.DimensionMismatch;
        var c = try Matrix.init(allocator, a.rows, b.cols);
        for (0..a.rows) |i| {
            for (0..b.cols) |j| {
                var sum: f64 = 0.0;
                for (0..a.cols) |k| sum += a.get(i, k) * b.get(k, j);
                c.set(i, j, sum);
            }
        }
        return c;
    }

    /// Return a new matrix equal to self transposed.
    pub fn transpose(self: Matrix) !Matrix {
        var t = try Matrix.init(self.allocator, self.cols, self.rows);
        for (0..self.rows) |i| {
            for (0..self.cols) |j| {
                t.set(j, i, self.get(i, j));
            }
        }
        return t;
    }

    /// Return a deep copy of this matrix (same allocator).
    pub fn clone(self: Matrix) !Matrix {
        const data = try self.allocator.dupe(f64, self.data);
        return Matrix{ .data = data, .rows = self.rows, .cols = self.cols, .allocator = self.allocator };
    }

    /// Return true when self is square and |a[i,j] - a[j,i]| <= tol for all i,j.
    pub fn isSymmetric(self: Matrix, tol: f64) bool {
        if (self.rows != self.cols) return false;
        const n = self.rows;
        for (0..n) |i| {
            for (i + 1..n) |j| {
                if (@abs(self.get(i, j) - self.get(j, i)) > tol) return false;
            }
        }
        return true;
    }
};

// --- Tests ---

test "Matrix.init: zero-initialised" {
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 2, 3);
    defer m.deinit();
    for (m.data) |v| try std.testing.expectEqual(@as(f64, 0.0), v);
    try std.testing.expectEqual(@as(usize, 2), m.rows);
    try std.testing.expectEqual(@as(usize, 3), m.cols);
}

test "Matrix.initIdentity" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 3);
    defer m.deinit();
    for (0..3) |i| {
        for (0..3) |j| {
            const expected: f64 = if (i == j) 1.0 else 0.0;
            try std.testing.expectEqual(expected, m.get(i, j));
        }
    }
}

test "Matrix.get and set" {
    const allocator = std.testing.allocator;
    var m = try Matrix.init(allocator, 2, 2);
    defer m.deinit();
    m.set(0, 1, 3.14);
    try std.testing.expectApproxEqAbs(@as(f64, 3.14), m.get(0, 1), 1e-10);
    try std.testing.expectEqual(@as(f64, 0.0), m.get(1, 0));
}

test "Matrix.add" {
    const allocator = std.testing.allocator;
    var a = try Matrix.init(allocator, 2, 2);
    defer a.deinit();
    var b = try Matrix.init(allocator, 2, 2);
    defer b.deinit();
    a.set(0, 0, 1.0);
    a.set(1, 1, 2.0);
    b.set(0, 0, 3.0);
    b.set(0, 1, 4.0);
    var c = try a.add(b);
    defer c.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), c.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), c.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), c.get(1, 1), 1e-10);
}

test "Matrix.scale in-place" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 2);
    defer m.deinit();
    m.scale(3.0);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), m.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), m.get(0, 1), 1e-10);
}

test "Matrix.multiply: 2x2 * 2x2" {
    const allocator = std.testing.allocator;
    // A = [[1,2],[3,4]]  B = [[5,6],[7,8]]
    // C = [[1*5+2*7, 1*6+2*8],[3*5+4*7, 3*6+4*8]] = [[19,22],[43,50]]
    var a = try Matrix.init(allocator, 2, 2);
    defer a.deinit();
    a.set(0, 0, 1);
    a.set(0, 1, 2);
    a.set(1, 0, 3);
    a.set(1, 1, 4);
    var b = try Matrix.init(allocator, 2, 2);
    defer b.deinit();
    b.set(0, 0, 5);
    b.set(0, 1, 6);
    b.set(1, 0, 7);
    b.set(1, 1, 8);
    var c = try Matrix.multiply(allocator, a, b);
    defer c.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 19.0), c.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 22.0), c.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 43.0), c.get(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 50.0), c.get(1, 1), 1e-10);
}

test "Matrix.multiply: dimension mismatch" {
    const allocator = std.testing.allocator;
    var a = try Matrix.init(allocator, 2, 3);
    defer a.deinit();
    var b = try Matrix.init(allocator, 2, 2);
    defer b.deinit();
    try std.testing.expectError(error.DimensionMismatch, Matrix.multiply(allocator, a, b));
}

test "Matrix.transpose" {
    const allocator = std.testing.allocator;
    // [[1,2,3],[4,5,6]] -> [[1,4],[2,5],[3,6]]
    var m = try Matrix.init(allocator, 2, 3);
    defer m.deinit();
    m.set(0, 0, 1);
    m.set(0, 1, 2);
    m.set(0, 2, 3);
    m.set(1, 0, 4);
    m.set(1, 1, 5);
    m.set(1, 2, 6);
    var t = try m.transpose();
    defer t.deinit();
    try std.testing.expectEqual(@as(usize, 3), t.rows);
    try std.testing.expectEqual(@as(usize, 2), t.cols);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), t.get(0, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), t.get(0, 1), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), t.get(1, 0), 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), t.get(2, 1), 1e-10);
}

test "Matrix.clone is independent" {
    const allocator = std.testing.allocator;
    var m = try Matrix.initIdentity(allocator, 2);
    defer m.deinit();
    var c = try m.clone();
    defer c.deinit();
    c.set(0, 0, 99.0);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), m.get(0, 0), 1e-10);
}

test "Matrix.isSymmetric" {
    const allocator = std.testing.allocator;
    var sym = try Matrix.initIdentity(allocator, 3);
    defer sym.deinit();
    sym.set(0, 1, 2.0);
    sym.set(1, 0, 2.0);
    try std.testing.expect(sym.isSymmetric(1e-10));

    var asym = try Matrix.initIdentity(allocator, 2);
    defer asym.deinit();
    asym.set(0, 1, 1.0);
    // asym[1,0] stays 0, so not symmetric
    try std.testing.expect(!asym.isSymmetric(1e-10));

    // Non-square is not symmetric
    var rect = try Matrix.init(allocator, 2, 3);
    defer rect.deinit();
    try std.testing.expect(!rect.isSymmetric(1e-10));
}
