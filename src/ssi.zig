// SSI (Simple Sequence Index) — module hub for sequence index formats.
//
// Re-exports ZeaselIndex (native format) and shared types.

const zeasel_ssi = @import("ssi/zeasel.zig");

pub const ZeaselIndex = zeasel_ssi.ZeaselIndex;

/// Temporary alias: will become a tagged union in a later task.
pub const SsiIndex = ZeaselIndex;

pub const SsiEntry = struct {
    name: []const u8,
    offset: u64,
    data_offset: u64,
    seq_len: u64,
    file_id: u16 = 0, // Index into file_names array for multi-file support
};

test {
    _ = zeasel_ssi;
}
