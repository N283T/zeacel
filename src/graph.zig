// Graph algorithms.
//
// Provides maximum bipartite matching using the Hopcroft-Karp algorithm.
// Used by mixture Dirichlet fitting (mixdchlet).

const std = @import("std");
const Allocator = std.mem.Allocator;

/// Maximum bipartite matching using augmenting paths (Hungarian method variant).
/// adj[i] is the list of vertices in the right set that left vertex i can match to.
/// n_left and n_right are the sizes of the two vertex sets.
/// Returns an array match_left[0..n_left-1] where match_left[i] is the matched
/// right vertex, or null if unmatched. Caller owns the returned slice.
pub fn maxBipartiteMatching(
    allocator: Allocator,
    n_left: usize,
    n_right: usize,
    adj: []const []const usize,
) ![]?usize {
    const match_left = try allocator.alloc(?usize, n_left);
    @memset(match_left, null);
    errdefer allocator.free(match_left);

    const match_right = try allocator.alloc(?usize, n_right);
    defer allocator.free(match_right);
    @memset(match_right, null);

    const visited = try allocator.alloc(bool, n_right);
    defer allocator.free(visited);

    // Try to find augmenting paths for each left vertex
    for (0..n_left) |u| {
        @memset(visited, false);
        _ = augment(u, adj, match_left, match_right, visited);
    }

    return match_left;
}

fn augment(
    u: usize,
    adj: []const []const usize,
    match_left: []?usize,
    match_right: []?usize,
    visited: []bool,
) bool {
    for (adj[u]) |v| {
        if (visited[v]) continue;
        visited[v] = true;

        // If v is free or we can reroute its current match
        if (match_right[v] == null or augment(match_right[v].?, adj, match_left, match_right, visited)) {
            match_left[u] = v;
            match_right[v] = u;
            return true;
        }
    }
    return false;
}

/// Count the number of matched pairs in a matching result.
pub fn matchingSize(matching: []const ?usize) usize {
    var count: usize = 0;
    for (matching) |m| {
        if (m != null) count += 1;
    }
    return count;
}

// --- Tests ---

test "maxBipartiteMatching: complete bipartite K3,3" {
    const allocator = std.testing.allocator;
    // K3,3: each left vertex connects to all 3 right vertices
    const adj0 = [_]usize{ 0, 1, 2 };
    const adj1 = [_]usize{ 0, 1, 2 };
    const adj2 = [_]usize{ 0, 1, 2 };
    const adj = [_][]const usize{ &adj0, &adj1, &adj2 };

    const matching = try maxBipartiteMatching(allocator, 3, 3, &adj);
    defer allocator.free(matching);

    // Perfect matching: all 3 should be matched
    try std.testing.expectEqual(@as(usize, 3), matchingSize(matching));
    // Each left vertex has a different right vertex
    try std.testing.expect(matching[0] != matching[1]);
    try std.testing.expect(matching[1] != matching[2]);
    try std.testing.expect(matching[0] != matching[2]);
}

test "maxBipartiteMatching: no edges" {
    const allocator = std.testing.allocator;
    const empty: [0]usize = .{};
    const adj = [_][]const usize{ &empty, &empty };

    const matching = try maxBipartiteMatching(allocator, 2, 2, &adj);
    defer allocator.free(matching);

    try std.testing.expectEqual(@as(usize, 0), matchingSize(matching));
}

test "maxBipartiteMatching: partial matching" {
    const allocator = std.testing.allocator;
    // Left 0 -> Right 0, Left 1 -> Right 0 (conflict), Left 2 -> Right 1
    const adj0 = [_]usize{0};
    const adj1 = [_]usize{0};
    const adj2 = [_]usize{1};
    const adj = [_][]const usize{ &adj0, &adj1, &adj2 };

    const matching = try maxBipartiteMatching(allocator, 3, 2, &adj);
    defer allocator.free(matching);

    // At most 2 matches possible (only 2 right vertices)
    try std.testing.expectEqual(@as(usize, 2), matchingSize(matching));
}
