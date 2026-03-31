// CPU feature detection for runtime SIMD dispatch.
//
// At comptime, Zig knows the target features. This module provides
// runtime detection for cases where the binary needs to select
// the fastest code path on the actual hardware.

const std = @import("std");
const builtin = @import("builtin");

/// SIMD instruction set support flags.
pub const SimdSupport = struct {
    sse: bool = false,
    sse2: bool = false,
    sse4_1: bool = false,
    avx: bool = false,
    avx2: bool = false,
    avx512f: bool = false,
    neon: bool = false,
};

/// Detect available SIMD features on the current CPU.
/// Uses comptime target info (Zig cross-compilation model means
/// runtime CPUID is rarely needed — the compiler already knows).
pub fn detect() SimdSupport {
    var support = SimdSupport{};

    switch (builtin.cpu.arch) {
        .x86_64, .x86 => {
            const features = builtin.cpu.features;
            support.sse = features.isEnabled(@intFromEnum(std.Target.x86.Feature.sse));
            support.sse2 = features.isEnabled(@intFromEnum(std.Target.x86.Feature.sse2));
            support.sse4_1 = features.isEnabled(@intFromEnum(std.Target.x86.Feature.sse4_1));
            support.avx = features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx));
            support.avx2 = features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx2));
            support.avx512f = features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx512f));
        },
        .aarch64 => {
            support.neon = true; // AArch64 always has NEON
        },
        else => {},
    }

    return support;
}

/// Return the number of available CPU cores.
pub fn cpuCount() usize {
    return std.Thread.getCpuCount() catch 1;
}

/// Return the best SIMD vector width in bytes for this platform.
pub fn bestVectorWidth() usize {
    const simd = detect();
    if (simd.avx512f) return 64;
    if (simd.avx2 or simd.avx) return 32;
    if (simd.sse or simd.neon) return 16;
    return 8; // scalar fallback
}

// --- Tests ---

test "detect: returns valid struct" {
    const simd = detect();
    // On x86_64, at minimum SSE2 should be available
    if (builtin.cpu.arch == .x86_64) {
        try std.testing.expect(simd.sse2);
    }
    // On aarch64, NEON should be available
    if (builtin.cpu.arch == .aarch64) {
        try std.testing.expect(simd.neon);
    }
}

test "cpuCount: at least 1" {
    try std.testing.expect(cpuCount() >= 1);
}

test "bestVectorWidth: at least 8" {
    try std.testing.expect(bestVectorWidth() >= 8);
}
