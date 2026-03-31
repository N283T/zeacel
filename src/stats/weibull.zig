/// Weibull distribution.
///
/// Parameters:
///   mu     - location parameter
///   lambda - scale parameter (must be > 0)
///   tau    - shape parameter (must be > 0)
///
/// Defined for x >= mu.
/// PDF: (tau * lambda) * (lambda*(x-mu))^(tau-1) * exp(-(lambda*(x-mu))^tau)
const std = @import("std");
const math = std.math;

/// PDF: (tau * lambda) * (lambda*(x-mu))^(tau-1) * exp(-(lambda*(x-mu))^tau)  for x >= mu, else 0.
pub fn pdf(x: f64, mu: f64, lambda: f64, tau: f64) f64 {
    if (x < mu) return 0.0;
    const u = lambda * (x - mu);
    if (u == 0.0) {
        // Limit depends on tau: pdf(mu) = tau*lambda when tau=1, 0 when tau>1, inf when tau<1
        if (tau < 1.0) return math.inf(f64);
        if (tau > 1.0) return 0.0;
        return lambda; // tau == 1: pure exponential
    }
    return (tau * lambda) * math.pow(f64, u, tau - 1.0) * @exp(-math.pow(f64, u, tau));
}

/// CDF: 1 - exp(-(lambda*(x-mu))^tau)  for x >= mu, else 0.
pub fn cdf(x: f64, mu: f64, lambda: f64, tau: f64) f64 {
    if (x < mu) return 0.0;
    const u = lambda * (x - mu);
    return 1.0 - @exp(-math.pow(f64, u, tau));
}

/// Survival function: exp(-(lambda*(x-mu))^tau)  for x >= mu, else 1.
pub fn surv(x: f64, mu: f64, lambda: f64, tau: f64) f64 {
    if (x < mu) return 1.0;
    const u = lambda * (x - mu);
    return @exp(-math.pow(f64, u, tau));
}

/// Inverse CDF: mu + (1/lambda) * (-log(1-p))^(1/tau)
pub fn invcdf(p: f64, mu: f64, lambda: f64, tau: f64) f64 {
    return mu + (1.0 / lambda) * math.pow(f64, -@log(1.0 - p), 1.0 / tau);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "weibull cdf at x=mu is 0" {
    try std.testing.expectEqual(@as(f64, 0.0), cdf(0.0, 0.0, 1.0, 1.0));
    try std.testing.expectEqual(@as(f64, 0.0), cdf(2.0, 2.0, 1.5, 2.0));
}

test "weibull cdf below mu is 0" {
    try std.testing.expectEqual(@as(f64, 0.0), cdf(-1.0, 0.0, 1.0, 2.0));
}

test "weibull surv below mu is 1" {
    try std.testing.expectEqual(@as(f64, 1.0), surv(-1.0, 0.0, 1.0, 2.0));
}

test "weibull cdf tau=1 is exponential special case" {
    // When tau=1, Weibull reduces to Exponential(lambda)
    // cdf(1, 0, 1, 1) = 1 - exp(-1) ≈ 0.6321
    const result = cdf(1.0, 0.0, 1.0, 1.0);
    const expected = 1.0 - @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-10));
}

test "weibull pdf at known value tau=1" {
    // pdf(1, 0, 1, 1) = lambda * exp(-1) = exp(-1) ≈ 0.36787944
    const result = pdf(1.0, 0.0, 1.0, 1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, @exp(-1.0), 1e-10));
}

test "weibull pdf at known value tau=2" {
    // pdf(1, 0, 1, 2) = 2*1*(1)^1 * exp(-1) = 2*exp(-1)
    const result = pdf(1.0, 0.0, 1.0, 2.0);
    const expected = 2.0 * @exp(-1.0);
    try std.testing.expect(math.approxEqAbs(f64, result, expected, 1e-10));
}

test "weibull pdf below mu is zero" {
    try std.testing.expectEqual(@as(f64, 0.0), pdf(-1.0, 0.0, 1.0, 2.0));
}

test "weibull cdf + surv = 1" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 5.0 };
    for (xs) |x| {
        const c = cdf(x, 0.0, 1.0, 2.0);
        const s = surv(x, 0.0, 1.0, 2.0);
        try std.testing.expect(math.approxEqAbs(f64, c + s, 1.0, 1e-14));
    }
}

test "weibull invcdf round-trip" {
    const mu: f64 = 1.0;
    const lambda: f64 = 0.5;
    const tau: f64 = 2.0;
    const xs = [_]f64{ 1.5, 2.0, 3.0, 5.0 };
    for (xs) |x| {
        const p = cdf(x, mu, lambda, tau);
        const x2 = invcdf(p, mu, lambda, tau);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-10));
    }
}

test "weibull invcdf tau=1 round-trip" {
    const xs = [_]f64{ 0.5, 1.0, 2.0, 4.0 };
    for (xs) |x| {
        const p = cdf(x, 0.0, 1.0, 1.0);
        const x2 = invcdf(p, 0.0, 1.0, 1.0);
        try std.testing.expect(math.approxEqAbs(f64, x, x2, 1e-10));
    }
}
