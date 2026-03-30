/// Gamma distribution and incomplete gamma function.
///
/// Parameters:
///   mu     - location (shift)
///   lambda - rate (lambda > 0)
///   a      - shape (a > 0)
///
/// The incomplete gamma functions are computed using series expansion
/// (for x < a+1) or continued fraction (for x >= a+1), following
/// Numerical Recipes §6.2.
const std = @import("std");
const math = std.math;

/// Maximum iterations for series / continued-fraction expansions.
const MAX_ITER: usize = 200;
/// Convergence threshold.
const EPSILON: f64 = 3e-15;
/// Small float to avoid division by zero in Lentz's method.
const FPMIN: f64 = 1e-300;

/// Regularized lower incomplete gamma P(a, x) = gamma(a,x) / Gamma(a).
/// Returns a value in [0, 1].
/// Uses series expansion for x < a+1, continued fraction otherwise.
pub fn incompleteGamma(a: f64, x: f64) f64 {
    if (x < 0.0 or a <= 0.0) return 0.0;
    if (x == 0.0) return 0.0;
    if (x < a + 1.0) {
        return gammaSeries(a, x);
    } else {
        return 1.0 - gammaContFrac(a, x);
    }
}

/// Series expansion for the regularized lower incomplete gamma.
/// Converges for x < a+1.
fn gammaSeries(a: f64, x: f64) f64 {
    const log_gamma_a = math.lgamma(f64, a);
    const ln_x = @log(x);

    // ap = a, sum = 1/a, del = 1/a
    var ap = a;
    var del: f64 = 1.0 / a;
    var sum = del;

    for (0..MAX_ITER) |_| {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (@abs(del) < @abs(sum) * EPSILON) break;
    }

    return sum * @exp(-x + a * ln_x - log_gamma_a);
}

/// Continued fraction expansion for the regularized upper incomplete gamma.
/// Returns Q(a, x) = 1 - P(a, x). Converges for x >= a+1.
/// Uses Lentz's modified continued fraction method.
fn gammaContFrac(a: f64, x: f64) f64 {
    const log_gamma_a = math.lgamma(f64, a);
    const ln_x = @log(x);

    // Set up for evaluating CF via modified Lentz's method.
    var b = x + 1.0 - a;
    var c: f64 = 1.0 / FPMIN;
    var d = 1.0 / b;
    var h = d;

    var i: f64 = 1.0;
    while (i <= @as(f64, @floatFromInt(MAX_ITER))) : (i += 1.0) {
        const an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (@abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (@abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        const delta = d * c;
        h *= delta;
        if (@abs(delta - 1.0) < EPSILON) break;
    }

    return @exp(-x + a * ln_x - log_gamma_a) * h;
}

/// PDF: (lambda^a / Gamma(a)) * (x-mu)^(a-1) * exp(-lambda*(x-mu))
/// Returns 0 for x <= mu.
pub fn pdf(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    if (x <= mu) return 0.0;
    const t = lambda * (x - mu);
    // log pdf = a*log(lambda) - lgamma(a) + (a-1)*log(x-mu) - lambda*(x-mu)
    //         = a*log(lambda) - lgamma(a) + (a-1)*log(t/lambda) - t
    //         = log(lambda) - lgamma(a) + (a-1)*log(t) - t
    const log_p = @log(lambda) - math.lgamma(f64, a) + (a - 1.0) * @log(t) - t;
    return @exp(log_p);
}

/// CDF via regularized incomplete gamma: P(a, lambda*(x-mu))
pub fn cdf(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    if (x <= mu) return 0.0;
    return incompleteGamma(a, lambda * (x - mu));
}

/// Survival = 1 - CDF
pub fn surv(x: f64, mu: f64, lambda: f64, a: f64) f64 {
    return 1.0 - cdf(x, mu, lambda, a);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "incompleteGamma(1, 1) = 1 - exp(-1)" {
    // P(1, x) = 1 - exp(-x)
    const result = incompleteGamma(1.0, 1.0);
    const expected = 1.0 - @exp(-1.0); // ≈ 0.63212056
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "incompleteGamma(2, 1) ≈ 0.26424" {
    // P(2, 1) = 1 - 2*exp(-1) ≈ 0.26424111...
    const result = incompleteGamma(2.0, 1.0);
    const expected: f64 = 1.0 - 2.0 * @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "incompleteGamma boundary: P(a,0) = 0" {
    try std.testing.expect(math.approxEqAbs(f64, incompleteGamma(2.0, 0.0), 0.0, 1e-14));
}

test "incompleteGamma large x approaches 1" {
    const result = incompleteGamma(2.0, 100.0);
    try std.testing.expect(math.approxEqAbs(f64, result, 1.0, 1e-10));
}

test "incompleteGamma continued fraction path (x >= a+1)" {
    // For a=2, x=5: uses continued fraction
    // Known: P(2,5) = 1 - 6*exp(-5) ≈ 0.95957...
    const result = incompleteGamma(2.0, 5.0);
    const expected = 1.0 - 6.0 * @exp(-5.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "gamma pdf positive for x > mu" {
    // With a=1, the gamma reduces to exponential: pdf = lambda*exp(-lambda*(x-mu))
    const result = pdf(1.0, 0.0, 1.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, @exp(-1.0), 1e-6));
}

test "gamma pdf zero at or below mu" {
    try std.testing.expectEqual(@as(f64, 0.0), pdf(0.0, 0.0, 1.0, 2.0));
    try std.testing.expectEqual(@as(f64, 0.0), pdf(-1.0, 0.0, 1.0, 2.0));
}

test "gamma cdf + surv = 1" {
    const x: f64 = 2.0;
    const c = cdf(x, 0.0, 1.0, 2.0);
    const s = surv(x, 0.0, 1.0, 2.0);
    try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-14));
}
