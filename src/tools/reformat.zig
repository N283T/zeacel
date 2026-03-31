// zeasel-reformat: convert between sequence file formats.
// Equivalent to Easel's esl-reformat.
//
// Usage: zeasel-reformat <output-format> <input-file>
// Where output-format is one of: fasta, genbank, embl

const std = @import("std");
const zeasel = @import("zeasel");

fn parseOutputFormat(name: []const u8) ?zeasel.io.Format {
    if (std.mem.eql(u8, name, "fasta")) return .fasta;
    if (std.mem.eql(u8, name, "genbank")) return .genbank;
    if (std.mem.eql(u8, name, "embl")) return .embl;
    return null;
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 3) {
        std.debug.print("Usage: zeasel-reformat <output-format> <input-file>\n", .{});
        std.debug.print("Output formats: fasta, genbank, embl\n", .{});
        std.process.exit(1);
    }

    const out_fmt_name = args[1];
    const path = args[2];

    const out_format = parseOutputFormat(out_fmt_name) orelse {
        std.debug.print("Error: unknown output format '{s}'\n", .{out_fmt_name});
        std.debug.print("Output formats: fasta, genbank, embl\n", .{});
        std.process.exit(1);
    };

    // Read entire file.
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    const data = try file.readToEndAlloc(allocator, 512 * 1024 * 1024);
    defer allocator.free(data);

    // Detect input format.
    const in_format = zeasel.io.Format.detect(data) orelse {
        std.debug.print("Error: cannot detect input file format\n", .{});
        std.process.exit(1);
    };

    // Guess alphabet from file content.
    const abc_type = zeasel.alphabet.guessType(data) orelse blk: {
        std.debug.print("Warning: could not detect alphabet, defaulting to DNA\n", .{});
        break :blk .dna;
    };
    const abc: *const zeasel.alphabet.Alphabet = switch (abc_type) {
        .dna => &zeasel.alphabet.dna,
        .rna => &zeasel.alphabet.rna,
        .amino => &zeasel.alphabet.amino,
    };

    // Parse all sequences.
    var reader = try zeasel.io.Reader.fromMemory(allocator, abc, data, in_format);
    defer reader.deinit();

    const sequences = try reader.readAll();
    defer {
        for (sequences) |*seq| @constCast(seq).deinit();
        allocator.free(sequences);
    }

    // Write sequences to stdout in the requested format.
    const stdout_file = std.fs.File.stdout();
    // TODO(M1): deprecatedWriter() is deprecated; migrate to stdout_file.writer(buffer)
    // once the codebase adopts the new std.Io.Writer API (requires an explicit buffer).
    const stdout = stdout_file.deprecatedWriter();
    var writer = zeasel.io.Writer.init(stdout.any(), out_format, 60);
    try writer.writeAll(sequences);
}
