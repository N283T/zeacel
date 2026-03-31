// Generalized single-linkage clustering.
//
// Provides a generic clustering function that works on any data type
// via a caller-supplied distance or linkage callback.

const std = @import("std");
const Allocator = std.mem.Allocator;

/// Single-linkage clustering using a pairwise linkage predicate.
/// `n` is the number of items. `linked` returns true if items i and j
/// should be in the same cluster (e.g., pairwise identity >= threshold).
/// Returns cluster assignment array: result[i] = cluster ID (0-indexed, contiguous).
/// Caller owns the returned slice.
pub fn singleLinkage(
    allocator: Allocator,
    n: usize,
    linked: *const fn (i: usize, j: usize) bool,
) ![]usize {
    // Union-find
    var parent = try allocator.alloc(usize, n);
    defer allocator.free(parent);
    for (0..n) |i| parent[i] = i;

    for (0..n) |i| {
        for (i + 1..n) |j| {
            if (linked(i, j)) {
                unionSets(parent, i, j);
            }
        }
    }

    return assignClusterIds(allocator, parent, n);
}

/// Single-linkage clustering from a flat n*n distance matrix.
/// Items with distance <= threshold are linked.
/// Returns cluster assignment array. Caller owns the returned slice.
pub fn singleLinkageFromDist(
    allocator: Allocator,
    dist: []const f64,
    n: usize,
    threshold: f64,
) ![]usize {
    var parent = try allocator.alloc(usize, n);
    defer allocator.free(parent);
    for (0..n) |i| parent[i] = i;

    for (0..n) |i| {
        for (i + 1..n) |j| {
            if (dist[i * n + j] <= threshold) {
                unionSets(parent, i, j);
            }
        }
    }

    return assignClusterIds(allocator, parent, n);
}

/// Count the number of distinct clusters in a cluster assignment array.
pub fn countClusters(assignments: []const usize) usize {
    if (assignments.len == 0) return 0;
    var max_id: usize = 0;
    for (assignments) |id| {
        if (id > max_id) max_id = id;
    }
    return max_id + 1;
}

fn find(parent: []usize, i: usize) usize {
    var x = i;
    while (parent[x] != x) x = parent[x];
    return x;
}

fn unionSets(parent: []usize, i: usize, j: usize) void {
    const ri = find(parent, i);
    const rj = find(parent, j);
    if (ri != rj) parent[ri] = rj;
}

fn assignClusterIds(allocator: Allocator, parent: []usize, n: usize) ![]usize {
    const cluster_id = try allocator.alloc(usize, n);
    errdefer allocator.free(cluster_id);

    const root_to_id = try allocator.alloc(usize, n);
    defer allocator.free(root_to_id);
    @memset(root_to_id, std.math.maxInt(usize));

    var next_id: usize = 0;
    for (0..n) |i| {
        const root = find(parent, i);
        if (root_to_id[root] == std.math.maxInt(usize)) {
            root_to_id[root] = next_id;
            next_id += 1;
        }
        cluster_id[i] = root_to_id[root];
    }

    return cluster_id;
}

// --- Tests ---

test "singleLinkageFromDist: 2 clusters" {
    const allocator = std.testing.allocator;
    // 4 items: (0,1) close, (2,3) close, far between groups
    const dist = [_]f64{
        0.0, 0.1, 0.9, 0.9,
        0.1, 0.0, 0.9, 0.9,
        0.9, 0.9, 0.0, 0.2,
        0.9, 0.9, 0.2, 0.0,
    };

    const clusters = try singleLinkageFromDist(allocator, &dist, 4, 0.5);
    defer allocator.free(clusters);

    try std.testing.expectEqual(clusters[0], clusters[1]);
    try std.testing.expectEqual(clusters[2], clusters[3]);
    try std.testing.expect(clusters[0] != clusters[2]);
    try std.testing.expectEqual(@as(usize, 2), countClusters(clusters));
}

test "singleLinkageFromDist: all in one cluster" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.3, 0.5,
        0.3, 0.0, 0.4,
        0.5, 0.4, 0.0,
    };

    const clusters = try singleLinkageFromDist(allocator, &dist, 3, 0.5);
    defer allocator.free(clusters);

    try std.testing.expectEqual(clusters[0], clusters[1]);
    try std.testing.expectEqual(clusters[1], clusters[2]);
    try std.testing.expectEqual(@as(usize, 1), countClusters(clusters));
}

test "singleLinkageFromDist: all singletons" {
    const allocator = std.testing.allocator;
    const dist = [_]f64{
        0.0, 0.9, 0.9,
        0.9, 0.0, 0.9,
        0.9, 0.9, 0.0,
    };

    const clusters = try singleLinkageFromDist(allocator, &dist, 3, 0.1);
    defer allocator.free(clusters);

    try std.testing.expect(clusters[0] != clusters[1]);
    try std.testing.expect(clusters[1] != clusters[2]);
    try std.testing.expectEqual(@as(usize, 3), countClusters(clusters));
}
