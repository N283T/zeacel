pub const alphabet = @import("alphabet.zig");
pub const sequence = @import("sequence.zig");
pub const io = @import("io.zig");

// Include tests from all sub-modules in `zig build test`.
test {
    _ = alphabet;
    _ = sequence;
    _ = io;
}
