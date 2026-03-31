const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("zeasel", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const lib = b.addLibrary(.{
        .name = "zeasel",
        .linkage = .static,
        .root_module = mod,
    });

    b.installArtifact(lib);

    const mod_tests = b.addTest(.{
        .root_module = mod,
    });

    const run_mod_tests = b.addRunArtifact(mod_tests);

    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_mod_tests.step);

    // CLI tools
    const tools_step = b.step("tools", "Build CLI tools");

    const tool_defs = [_]struct { name: []const u8, src: []const u8 }{
        .{ .name = "zeasel-seqstat", .src = "src/tools/seqstat.zig" },
        .{ .name = "zeasel-reformat", .src = "src/tools/reformat.zig" },
        .{ .name = "zeasel-seqfetch", .src = "src/tools/seqfetch.zig" },
    };

    inline for (tool_defs) |tool| {
        const exe = b.addExecutable(.{
            .name = tool.name,
            .root_module = b.createModule(.{
                .root_source_file = b.path(tool.src),
                .target = target,
                .optimize = optimize,
                .imports = &.{
                    .{ .name = "zeasel", .module = mod },
                },
            }),
        });
        b.installArtifact(exe);
        tools_step.dependOn(&exe.step);
    }
}
