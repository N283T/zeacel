/// Exponential distribution with location (mu) and rate (lambda).
///
/// Defined for x >= mu.
/// PDF: lambda * exp(-lambda * (x - mu))
const std = @import("std");
const math = std.math;

/// PDF: lambda * exp(-lambda * (x - mu)) for x >= mu, else 0.
pub fn pdf(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return 0.0;
    return lambda * @exp(-lambda * (x - mu));
}

/// CDF: 1 - exp(-lambda * (x - mu)) for x >= mu, else 0.
pub fn cdf(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return 0.0;
    return 1.0 - @exp(-lambda * (x - mu));
}

/// Survival function: exp(-lambda * (x - mu)) for x >= mu, else 1.
pub fn surv(x: f64, mu: f64, lambda: f64) f64 {
    if (x < mu) return 1.0;
    return @exp(-lambda * (x - mu));
}

/// Inverse CDF: mu - (1/lambda) * log(1 - p)
pub fn invcdf(p: f64, mu: f64, lambda: f64) f64 {
    return mu - (1.0 / lambda) * @log(1.0 - p);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "exponential pdf at x=1, mu=0, lambda=1" {
    // pdf = exp(-1) ≈ 0.36787944
    const result = pdf(1.0, 0.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, @exp(-1.0), 1e-6));
}

test "exponential cdf at x=1, mu=0, lambda=1" {
    // cdf = 1 - exp(-1) ≈ 0.63212056
    const result = cdf(1.0, 0.0, 1.0);
    const expected = 1.0 - @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-6));
}

test "exponential pdf below mu is zero" {
    try std.testing.expectEqual(@as(f64, 0.0), pdf(-1.0, 0.0, 1.0));
}

test "exponential cdf below mu is zero" {
    try std.testing.expectEqual(@as(f64, 0.0), cdf(-1.0, 0.0, 1.0));
}

test "exponential surv below mu is one" {
    try std.testing.expectEqual(@as(f64, 1.0), surv(-1.0, 0.0, 1.0));
}

test "exponential cdf + surv = 1" {
    const xs = [_]f64{ 0.0, 0.5, 1.0, 3.0 };
    for (xs) |x| {
        const c = cdf(x, 0.0, 1.0);
        const s = surv(x, 0.0, 1.0);
        try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-14));
    }
}

test "exponential invcdf round-trip" {
    const mu: f64 = 2.0;
    const lambda: f64 = 1.5;
    const xs = [_]f64{ 2.0, 2.5, 3.0, 5.0 };
    for (xs) |x| {
        const p = cdf(x, mu, lambda);
        const x2 = invcdf(p, mu, lambda);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-6));
    }
}
