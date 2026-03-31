// Seeded pseudo-random number generator wrapper.
//
// Wraps std.Random.Pcg to provide a reproducible, seedable RNG with
// utility methods matching Easel-style usage patterns.

const std = @import("std");
const math = std.math;

pub const Random = struct {
    pcg: std.Random.Pcg,

    /// Create a new RNG with a specific seed.
    pub fn init(seed: u64) Random {
        return Random{ .pcg = std.Random.Pcg.init(seed) };
    }

    /// Return the underlying std.Random interface.
    fn rng(self: *Random) std.Random {
        return self.pcg.random();
    }

    /// Uniform random float in [0, 1).
    pub fn uniform(self: *Random) f64 {
        return self.pcg.random().float(f64);
    }

    /// Uniform random integer in [0, n).
    pub fn uniformInt(self: *Random, n: u32) u32 {
        return self.pcg.random().intRangeLessThan(u32, 0, n);
    }

    /// Gaussian (normal) random variate with mean 0, standard deviation 1.
    /// Uses the Box-Muller transform. Consumes two uniform samples per call.
    pub fn gaussian(self: *Random) f64 {
        // Box-Muller: given u1, u2 uniform in (0, 1),
        // z = sqrt(-2 * ln(u1)) * cos(2*pi*u2) is standard normal.
        // Guard against u1 == 0 to avoid ln(0).
        var s: f64 = 0.0;
        while (s == 0.0) {
            s = self.uniform();
        }
        const t = self.uniform();
        const mag = math.sqrt(-2.0 * math.log(f64, math.e, s));
        return mag * math.cos(2.0 * math.pi * t);
    }

    /// Uniform random float in (0, 1), excluding 0.
    /// Important for log-based sampling where ln(0) is undefined.
    pub fn uniformPositive(self: *Random) f64 {
        var x: f64 = 0.0;
        while (x == 0.0) {
            x = self.uniform();
        }
        return x;
    }

    /// Gamma-distributed random deviate with shape parameter `alpha` and scale 1.
    /// Uses Knuth's algorithms: integer gamma for small integer alpha,
    /// Ahrens' method for alpha >= 3, Knuth's exercise 16 for fractional parts.
    pub fn gamma(self: *Random, alpha: f64) f64 {
        std.debug.assert(alpha > 0.0);

        const aint = @floor(alpha);
        if (alpha == aint and alpha < 12.0) {
            return gammaInteger(self, @intFromFloat(alpha));
        } else if (alpha > 3.0) {
            return gammaAhrens(self, alpha);
        } else if (alpha < 1.0) {
            return gammaFraction(self, alpha);
        } else {
            return gammaInteger(self, @intFromFloat(aint)) + gammaFraction(self, alpha - aint);
        }
    }

    /// Gamma deviate for small positive integer alpha (alpha < 12).
    /// Product of alpha uniform positive deviates, then -log.
    fn gammaInteger(self: *Random, a: u32) f64 {
        var u: f64 = 1.0;
        for (0..a) |_| {
            u *= self.uniformPositive();
        }
        return -math.log(f64, math.e, u);
    }

    /// Gamma deviate for alpha >= 3 using Ahrens-Dieter acceptance-rejection.
    fn gammaAhrens(self: *Random, a: f64) f64 {
        while (true) {
            var y: f64 = undefined;
            var x: f64 = undefined;
            while (true) {
                y = math.tan(math.pi * self.uniform());
                x = y * math.sqrt(2.0 * a - 1.0) + a - 1.0;
                if (x > 0.0) break;
            }
            const v = self.uniform();
            const test_val = (1.0 + y * y) * math.exp((a - 1.0) * math.log(f64, math.e, x / (a - 1.0)) - y * math.sqrt(2.0 * a - 1.0));
            if (v <= test_val) return x;
        }
    }

    /// Gamma deviate for fractional alpha in (0, 1) using Knuth 3.4.1 exercise 16.
    fn gammaFraction(self: *Random, a: f64) f64 {
        const p = math.e / (a + math.e);
        while (true) {
            const u = self.uniform();
            const v = self.uniformPositive();
            var x: f64 = undefined;
            var q: f64 = undefined;
            if (u < p) {
                x = math.pow(f64, v, 1.0 / a);
                q = math.exp(-x);
            } else {
                x = 1.0 - math.log(f64, math.e, v);
                q = math.pow(f64, x, a - 1.0);
            }
            const u_accept = self.uniform();
            if (u_accept < q) return x;
        }
    }

    /// Random choice from a probability distribution p[0..n-1].
    /// p must sum to approximately 1.0 and all values must be non-negative.
    /// Returns the index of the chosen residue.
    pub fn choose(self: *Random, p: []const f64) usize {
        const u = self.uniform();
        var cumulative: f64 = 0.0;
        for (p, 0..) |prob, i| {
            cumulative += prob;
            if (u < cumulative) return i;
        }
        // Return last index if floating-point rounding pushes u >= sum(p).
        return p.len - 1;
    }
};

// --- Tests ---

test "init: reproducibility — same seed produces same sequence" {
    var r1 = Random.init(42);
    var r2 = Random.init(42);

    for (0..20) |_| {
        try std.testing.expectEqual(r1.uniform(), r2.uniform());
    }
}

test "init: different seeds produce different sequences" {
    var r1 = Random.init(1);
    var r2 = Random.init(2);

    // Extremely unlikely all 10 values match by chance.
    var all_equal = true;
    for (0..10) |_| {
        if (r1.uniform() != r2.uniform()) {
            all_equal = false;
            break;
        }
    }
    try std.testing.expect(!all_equal);
}

test "uniform: all values in [0, 1)" {
    var rng = Random.init(123);
    for (0..1000) |_| {
        const v = rng.uniform();
        try std.testing.expect(v >= 0.0);
        try std.testing.expect(v < 1.0);
    }
}

test "uniformInt: all values in [0, n)" {
    var rng = Random.init(456);
    const n: u32 = 7;
    for (0..1000) |_| {
        const v = rng.uniformInt(n);
        try std.testing.expect(v < n);
    }
}

test "uniformInt: n=1 always returns 0" {
    var rng = Random.init(789);
    for (0..100) |_| {
        try std.testing.expectEqual(@as(u32, 0), rng.uniformInt(1));
    }
}

test "choose: uniform distribution — all 4 choices appear over 1000 draws" {
    var rng = Random.init(999);
    const freq = [_]f64{ 0.25, 0.25, 0.25, 0.25 };
    var counts = [_]usize{0} ** 4;

    for (0..1000) |_| {
        const idx = rng.choose(&freq);
        counts[idx] += 1;
    }

    for (counts) |c| {
        try std.testing.expect(c > 0);
    }
}

test "choose: deterministic with biased distribution" {
    var rng = Random.init(0);
    // Always pick index 2 — probability 1.0
    const freq = [_]f64{ 0.0, 0.0, 1.0, 0.0 };
    for (0..100) |_| {
        try std.testing.expectEqual(@as(usize, 2), rng.choose(&freq));
    }
}

test "uniformPositive: all values in (0, 1)" {
    var rng = Random.init(42);
    for (0..1000) |_| {
        const v = rng.uniformPositive();
        try std.testing.expect(v > 0.0);
        try std.testing.expect(v < 1.0);
    }
}

test "gamma: integer alpha — mean near alpha over 2000 samples" {
    var rng = Random.init(7777);
    const alpha = 5.0;
    var sum: f64 = 0.0;
    const n: f64 = 2000.0;
    for (0..@intFromFloat(n)) |_| {
        const v = rng.gamma(alpha);
        try std.testing.expect(v > 0.0);
        sum += v;
    }
    const mean = sum / n;
    // Gamma(alpha, 1) has mean = alpha. Tolerance ~0.5 for 2000 samples.
    try std.testing.expect(@abs(mean - alpha) < 0.5);
}

test "gamma: fractional alpha — mean near alpha over 2000 samples" {
    var rng = Random.init(8888);
    const alpha = 0.5;
    var sum: f64 = 0.0;
    const n: f64 = 2000.0;
    for (0..@intFromFloat(n)) |_| {
        const v = rng.gamma(alpha);
        try std.testing.expect(v > 0.0);
        sum += v;
    }
    const mean = sum / n;
    try std.testing.expect(@abs(mean - alpha) < 0.3);
}

test "gamma: large alpha (Ahrens) — mean near alpha" {
    var rng = Random.init(9999);
    const alpha = 10.0;
    var sum: f64 = 0.0;
    const n: f64 = 2000.0;
    for (0..@intFromFloat(n)) |_| {
        const v = rng.gamma(alpha);
        try std.testing.expect(v > 0.0);
        sum += v;
    }
    const mean = sum / n;
    try std.testing.expect(@abs(mean - alpha) < 1.0);
}

test "gaussian: mean near 0 and std near 1 over 1000 samples" {
    var rng = Random.init(314159);
    var sum: f64 = 0.0;
    var sum_sq: f64 = 0.0;
    const n: f64 = 1000.0;

    for (0..@intFromFloat(n)) |_| {
        const v = rng.gaussian();
        sum += v;
        sum_sq += v * v;
    }

    const mean = sum / n;
    const variance = sum_sq / n - mean * mean;
    const std_dev = @sqrt(variance);

    // With 1000 samples the CLT guarantees these tolerances hold with very high
    // probability (mean within 0.15, std_dev within 0.15 of 1.0).
    try std.testing.expect(@abs(mean) < 0.15);
    try std.testing.expect(@abs(std_dev - 1.0) < 0.15);
}
