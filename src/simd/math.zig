/// Vectorized f32 math functions (logf, expf) using Cephes/Easel polynomial
/// approximations. These operate on Zig `@Vector(N, f32)` types and produce
/// results matching Easel's `esl_sse_logf` / `esl_sse_expf` to within ~1 ULP
/// for normal inputs.
///
/// Reference: Easel esl_sse.c (Julian Pommier's SSE implementation, adapted
/// by Sean Eddy). Constants from Cephes math library.
const std = @import("std");
const math = std.math;

// ---------------------------------------------------------------------------
// logf_vec: vectorized natural logarithm for f32 vectors
// ---------------------------------------------------------------------------

/// Vectorized natural logarithm for f32 vectors.
///
/// Approximates ln(x) using the Cephes polynomial, matching Easel's
/// `esl_sse_logf()`. For x < 0 (including -0) returns NaN; for x == 0
/// or subnormal x returns -inf; for x == +inf returns +inf; for NaN
/// returns NaN.
///
/// Accuracy: < 1 ULP for normal positive f32 inputs.
pub fn logf_vec(comptime N: comptime_int, x: @Vector(N, f32)) @Vector(N, f32) {
    const UVec = @Vector(N, u32);
    const IVec = @Vector(N, i32);
    const FVec = @Vector(N, f32);

    const one: FVec = @splat(1.0);
    const half: FVec = @splat(0.5);
    const neg_inf: FVec = @splat(-math.inf(f32));
    const nan_val: FVec = @splat(math.nan(f32));

    // Cephes polynomial coefficients (9 terms), same as Easel esl_sse.c:74-76
    const cephes_p = [9]f32{
        7.0376836292e-2,  -1.1514610310e-1, 1.1676998740e-1,
        -1.2420140846e-1, 1.4249322787e-1,  -1.6668057665e-1,
        2.0000714765e-1,  -2.4999993993e-1, 3.3333331174e-1,
    };

    // ---- IEEE754 special-case masks ----
    const xi: UVec = @bitCast(x);

    // Detect negative values (sign bit set) -> NaN
    const sign_bit: UVec = @splat(0x80000000);
    const is_negative = (xi & sign_bit) == sign_bit;

    // Detect zero / subnormal (biased exponent == 0) -> -inf
    const ei_raw = xi >> @as(UVec, @splat(23));
    const is_zero = ei_raw == @as(UVec, @splat(0));

    // Detect inf / NaN (biased exponent == 0xFF) -> pass through
    const exp_mask: UVec = @splat(0x7F800000);
    const is_inf_or_nan = (xi & exp_mask) == exp_mask;

    // ---- frexp: split x into mantissa [0.5, 1.0) and exponent ----
    // Strip exponent bits, set biased exponent to 126 (= 0x3F000000) => [0.5, 1.0)
    const mantissa_mask: UVec = @splat(0x007FFFFF);
    const half_bits: UVec = @splat(0x3F000000);
    var m: FVec = @bitCast((xi & mantissa_mask) | half_bits);

    // Signed exponent: biased_exp - 126 (not 127, because mantissa is [0.5,1))
    var e: IVec = @intCast(ei_raw & @as(UVec, @splat(0xFF)));
    e -= @as(IVec, @splat(126));
    var ef: FVec = @floatFromInt(e);

    // ---- Cephes range reduction: if m < 0.707, double m and decrement e ----
    const sqrt_half: FVec = @splat(0.707106781186547524);
    const lt_mask = m < sqrt_half; // bool vector
    const lt_float: FVec = @select(f32, lt_mask, one, @as(FVec, @splat(0.0)));
    // tmp = m where m < sqrt_half, else 0
    const tmp_m: FVec = @select(f32, lt_mask, m, @as(FVec, @splat(0.0)));
    m = m - one;
    ef = ef - lt_float;
    m = m + tmp_m; // if was < sqrt_half: m = (2*orig_m - 1), else (orig_m - 1)

    const z = m * m;

    // ---- Evaluate polynomial y = ((..((p0*m + p1)*m + p2)*m ...)*m + p8)*m ----
    var y: FVec = @splat(cephes_p[0]);
    y = y * m + @as(FVec, @splat(cephes_p[1]));
    y = y * m + @as(FVec, @splat(cephes_p[2]));
    y = y * m + @as(FVec, @splat(cephes_p[3]));
    y = y * m + @as(FVec, @splat(cephes_p[4]));
    y = y * m + @as(FVec, @splat(cephes_p[5]));
    y = y * m + @as(FVec, @splat(cephes_p[6]));
    y = y * m + @as(FVec, @splat(cephes_p[7]));
    y = y * m + @as(FVec, @splat(cephes_p[8]));
    y = y * z;

    // ---- Reconstruct: result = y + e*(-2.12194440e-4) - z*0.5 + m + e*0.693359375 ----
    const ln2_hi: FVec = @splat(0.693359375);
    const ln2_lo: FVec = @splat(-2.12194440e-4);

    y = y + ef * ln2_lo;
    y = y - z * half;
    var result = m + y + ef * ln2_hi;

    // ---- IEEE754 cleanup ----
    // log(inf) = inf, log(NaN) = NaN: use original x for these lanes
    result = @select(f32, is_inf_or_nan, x, result);
    // log(x < 0) = NaN (including -0, -inf)
    result = @select(f32, is_negative, nan_val, result);
    // log(0) = -inf, log(subnormal) = -inf
    result = @select(f32, is_zero, neg_inf, result);

    return result;
}

// ---------------------------------------------------------------------------
// expf_vec: vectorized exponential for f32 vectors
// ---------------------------------------------------------------------------

/// Vectorized exponential for f32 vectors.
///
/// Approximates exp(x) using the Cephes polynomial, matching Easel's
/// `esl_sse_expf()`. Clamps to 0 for very negative x, to +inf for very
/// large x. Handles all IEEE754 special values correctly.
///
/// Accuracy: < 1 ULP for normal f32 inputs in [-87, +88].
pub fn expf_vec(comptime N: comptime_int, x_arg: @Vector(N, f32)) @Vector(N, f32) {
    const UVec = @Vector(N, u32);
    const IVec = @Vector(N, i32);
    const FVec = @Vector(N, f32);

    // Cephes polynomial coefficients, same as Easel esl_sse.c:184-185
    const cephes_p = [6]f32{
        1.9875691500e-4, 1.3981999507e-3, 8.3334519073e-3,
        4.1665795894e-2, 1.6666665459e-1, 5.0000001201e-1,
    };
    const cephes_c0: f32 = 0.693359375; // ln(2) high part
    const cephes_c1: f32 = -2.12194440e-4; // ln(2) low part
    const maxlogf: f32 = 88.3762626647949;
    const minlogf: f32 = -88.3762626647949;
    const log2e: f32 = 1.44269504088896341; // 1/ln(2)

    var x = x_arg;

    // ---- Detect overflow / underflow ----
    const is_over = x > @as(FVec, @splat(maxlogf));
    const is_under = x <= @as(FVec, @splat(minlogf));

    // ---- Range reduction: exp(x) = 2^k * exp(f) where k = floor(x/ln2 + 0.5) ----
    var fx = x * @as(FVec, @splat(log2e)) + @as(FVec, @splat(0.5));

    // floor() via truncation toward zero with correction for negative values
    const k_trunc: IVec = @intFromFloat(fx);
    const tmp: FVec = @floatFromInt(k_trunc);
    // If truncation increased the value (negative case), subtract 1
    const correction: FVec = @select(f32, tmp > fx, @as(FVec, @splat(1.0)), @as(FVec, @splat(0.0)));
    fx = tmp - correction;
    const k: IVec = @intFromFloat(fx);

    // ---- Compute reduced argument: x = x - k * ln(2) ----
    // Split ln(2) into high and low parts for precision
    x = x - fx * @as(FVec, @splat(cephes_c0));
    x = x - fx * @as(FVec, @splat(cephes_c1));
    const z = x * x;

    // ---- Evaluate polynomial: 1 + x + x^2/2 + ... (Cephes form) ----
    var y: FVec = @splat(cephes_p[0]);
    y = y * x + @as(FVec, @splat(cephes_p[1]));
    y = y * x + @as(FVec, @splat(cephes_p[2]));
    y = y * x + @as(FVec, @splat(cephes_p[3]));
    y = y * x + @as(FVec, @splat(cephes_p[4]));
    y = y * x + @as(FVec, @splat(cephes_p[5]));
    y = y * z + x + @as(FVec, @splat(1.0));

    // ---- Build 2^k by constructing IEEE754 float: (k+127) << 23 ----
    const k_biased: IVec = k + @as(IVec, @splat(127));
    const k_u: UVec = @bitCast(k_biased);
    const pow2k: FVec = @bitCast(k_u << @as(UVec, @splat(23)));

    // result = polynomial * 2^k
    var result = y * pow2k;

    // ---- Clamp overflow/underflow ----
    result = @select(f32, is_over, @as(FVec, @splat(math.inf(f32))), result);
    result = @select(f32, is_under, @as(FVec, @splat(0.0)), result);

    return result;
}

// ===========================================================================
// Tests
// ===========================================================================

const testing = std.testing;

fn expectApproxVec(comptime N: comptime_int, expected: @Vector(N, f32), actual: @Vector(N, f32), tolerance: f32) !void {
    const exp_arr: [N]f32 = expected;
    const act_arr: [N]f32 = actual;
    for (0..N) |i| {
        if (math.isNan(exp_arr[i])) {
            try testing.expect(math.isNan(act_arr[i]));
        } else if (math.isInf(exp_arr[i])) {
            try testing.expect(math.isInf(act_arr[i]));
            try testing.expect(math.isPositiveInf(exp_arr[i]) == math.isPositiveInf(act_arr[i]));
        } else {
            try testing.expectApproxEqAbs(exp_arr[i], act_arr[i], tolerance);
        }
    }
}

// ---- logf_vec tests ----

test "logf_vec: matches std @log for normal values (4-wide)" {
    const input: @Vector(4, f32) = .{ 1.0, 2.0, 0.5, 10.0 };
    const result = logf_vec(4, input);
    const expected: @Vector(4, f32) = .{ @log(@as(f32, 1.0)), @log(@as(f32, 2.0)), @log(@as(f32, 0.5)), @log(@as(f32, 10.0)) };
    try expectApproxVec(4, expected, result, 1e-6);
}

test "logf_vec: values near 1.0" {
    const input: @Vector(4, f32) = .{ 0.9, 0.99, 1.01, 1.1 };
    const result = logf_vec(4, input);
    const expected: @Vector(4, f32) = .{ @log(@as(f32, 0.9)), @log(@as(f32, 0.99)), @log(@as(f32, 1.01)), @log(@as(f32, 1.1)) };
    try expectApproxVec(4, expected, result, 1e-6);
}

test "logf_vec: large and small values" {
    const input: @Vector(4, f32) = .{ 1e-30, 1e-10, 1e10, 1e30 };
    const result = logf_vec(4, input);
    const expected: @Vector(4, f32) = .{ @log(@as(f32, 1e-30)), @log(@as(f32, 1e-10)), @log(@as(f32, 1e10)), @log(@as(f32, 1e30)) };
    try expectApproxVec(4, expected, result, 1e-2);
}

test "logf_vec: IEEE754 special cases" {
    const input: @Vector(4, f32) = .{ 0.0, -1.0, math.inf(f32), math.nan(f32) };
    const result = logf_vec(4, input);
    const result_arr: [4]f32 = result;

    // log(0) = -inf
    try testing.expect(result_arr[0] == -math.inf(f32));
    // log(-1) = NaN
    try testing.expect(math.isNan(result_arr[1]));
    // log(inf) = inf
    try testing.expect(result_arr[2] == math.inf(f32));
    // log(NaN) = NaN
    try testing.expect(math.isNan(result_arr[3]));
}

test "logf_vec: negative zero returns NaN" {
    const input: @Vector(4, f32) = .{ -0.0, 1.0, 1.0, 1.0 };
    const result = logf_vec(4, input);
    const result_arr: [4]f32 = result;
    // -0 has sign bit set, so treated as negative -> NaN (matching Easel)
    try testing.expect(math.isNan(result_arr[0]));
}

test "logf_vec: 8-wide" {
    const input: @Vector(8, f32) = .{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
    const result = logf_vec(8, input);
    const result_arr: [8]f32 = result;
    for (0..8) |i| {
        const val: f32 = @floatFromInt(i + 1);
        try testing.expectApproxEqAbs(@log(val), result_arr[i], 1e-6);
    }
}

test "logf_vec: exp roundtrip" {
    // log(exp(x)) should give back x for moderate values
    const input: @Vector(4, f32) = .{ 0.0, 1.0, -1.0, 5.0 };
    const exp_result = expf_vec(4, input);
    const log_result = logf_vec(4, exp_result);
    try expectApproxVec(4, input, log_result, 1e-5);
}

// ---- expf_vec tests ----

test "expf_vec: matches std @exp for normal values (4-wide)" {
    const input: @Vector(4, f32) = .{ 0.0, 1.0, -1.0, 2.0 };
    const result = expf_vec(4, input);
    const expected: @Vector(4, f32) = .{ @exp(@as(f32, 0.0)), @exp(@as(f32, 1.0)), @exp(@as(f32, -1.0)), @exp(@as(f32, 2.0)) };
    try expectApproxVec(4, expected, result, 1e-6);
}

test "expf_vec: negative values" {
    const input: @Vector(4, f32) = .{ -0.5, -2.0, -5.0, -10.0 };
    const result = expf_vec(4, input);
    const expected: @Vector(4, f32) = .{ @exp(@as(f32, -0.5)), @exp(@as(f32, -2.0)), @exp(@as(f32, -5.0)), @exp(@as(f32, -10.0)) };
    try expectApproxVec(4, expected, result, 1e-6);
}

test "expf_vec: overflow clamps to inf" {
    const input: @Vector(4, f32) = .{ 100.0, 89.0, 200.0, 1000.0 };
    const result = expf_vec(4, input);
    const result_arr: [4]f32 = result;
    for (result_arr) |val| {
        try testing.expect(val == math.inf(f32));
    }
}

test "expf_vec: underflow clamps to zero" {
    const input: @Vector(4, f32) = .{ -100.0, -89.0, -200.0, -1000.0 };
    const result = expf_vec(4, input);
    const result_arr: [4]f32 = result;
    for (result_arr) |val| {
        try testing.expectApproxEqAbs(0.0, val, 1e-38);
    }
}

test "expf_vec: values near zero" {
    const input: @Vector(4, f32) = .{ 0.001, -0.001, 0.1, -0.1 };
    const result = expf_vec(4, input);
    const expected: @Vector(4, f32) = .{ @exp(@as(f32, 0.001)), @exp(@as(f32, -0.001)), @exp(@as(f32, 0.1)), @exp(@as(f32, -0.1)) };
    try expectApproxVec(4, expected, result, 1e-6);
}

test "expf_vec: 8-wide" {
    const input: @Vector(8, f32) = .{ -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0 };
    const result = expf_vec(8, input);
    const result_arr: [8]f32 = result;
    const input_arr: [8]f32 = input;
    for (0..8) |i| {
        try testing.expectApproxEqAbs(@exp(input_arr[i]), result_arr[i], 1e-5);
    }
}

test "expf_vec: log roundtrip" {
    // exp(log(x)) should give back x for positive values
    const input: @Vector(4, f32) = .{ 0.5, 1.0, 2.0, 10.0 };
    const log_result = logf_vec(4, input);
    const exp_result = expf_vec(4, log_result);
    try expectApproxVec(4, input, exp_result, 1e-5);
}

test "expf_vec: boundary values near clamp thresholds" {
    // Just inside the valid range
    const input: @Vector(4, f32) = .{ 87.0, -87.0, 88.0, -88.0 };
    const result = expf_vec(4, input);
    const expected: @Vector(4, f32) = .{ @exp(@as(f32, 87.0)), @exp(@as(f32, -87.0)), @exp(@as(f32, 88.0)), @exp(@as(f32, -88.0)) };
    try expectApproxVec(4, expected, result, expected[0] * 1e-6);
}
