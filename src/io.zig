// I/O subsystem: parsers and writers for biological sequence file formats.

pub const fasta = @import("io/fasta.zig");
pub const Reader = @import("io/reader.zig").Reader;
pub const Writer = @import("io/writer.zig").Writer;
pub const Format = @import("io/reader.zig").Format;

const reader_mod = @import("io/reader.zig");
const writer_mod = @import("io/writer.zig");

// Include tests from all I/O sub-modules in `zig build test`.
test {
    _ = fasta;
    _ = reader_mod;
    _ = writer_mod;
}
