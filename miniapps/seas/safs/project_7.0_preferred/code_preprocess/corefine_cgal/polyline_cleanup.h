#pragma once

// polyline_cleanup.h — deterministic shared-polyline edge collapse across
// multiple meshes, preserving pairwise conformality.
//
// Why this exists: PMP::corefine inserts the *raw* exact-predicate
// triangle-triangle intersection polyline into both meshes of a pair.
// Where the intersection passes near a triangle vertex (or where two
// intersection polylines meet at a triple junction), the raw polyline can
// contain segments far below the user-specified `min_edge` floor.  The
// downstream PMP::isotropic_remeshing call must use protect_constraints
// =true to keep the polyline conformal across the pair, so it cannot
// collapse those short edges itself.
//
// Approach: explicitly collapse polyline edges shorter than `min_edge` in
// every mesh that contains them, with bit-identical decisions and
// placements across all containing meshes:
//
//   (a) An edge is identified across meshes by a canonical key built
//       from its two endpoint coordinates (lex-sorted, integer-quantized
//       to coord_tol = 1e-6 m).  Same coords in mesh A and mesh B map to
//       the same key.  Conformality after corefine guarantees that the
//       polyline endpoint coords are bit-identical across A and B
//       (max_diff_m == 0 in the manifest).
//
//   (b) Placement uses the midpoint of the two endpoint coords.  Since
//       the endpoints are bit-identical across meshes, so is the midpoint.
//
//   (c) Atomic per-edge: an edge is collapsed in EVERY mesh containing
//       it, or NONE.  The link-condition test (does_satisfy_link_condition)
//       depends on local mesh topology, which differs between mesh A and
//       B around the polyline.  If the test fails in any one containing
//       mesh, we skip the collapse in all of them — preserving
//       conformality at the cost of leaving that one polyline edge sub-
//       floor.  The caller decides what to do with residual short edges
//       (the SAFS pipeline hard-fails).

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

namespace polyline_cleanup {

struct Stats {
    std::size_t n_iterations          = 0;
    std::size_t n_collapses_done      = 0;
    std::size_t n_residual_short      = 0;
    std::size_t n_precondition_skips  = 0;  // border-ear configs (see below)
    double      residual_min_length   = 0.0;
    double      residual_max_below    = 0.0;  // (min_edge - residual_min_length)
};

// Release-mode safety check for CGAL::Euler::collapse_edge.  Its
// preconditions are compiled out under NDEBUG; the unchecked case that
// SILENTLY CORRUPTS the mesh (measured 2026-06-12: SIGBUS / mega-degree
// umbrella spins on the SAFv4 fault x bottom border polyline) is an
// incident triangle whose third edge is ALSO a border edge (an "ear" at
// a patch border): collapse_edge then calls remove_face on a border
// halfedge.  Mirror the preconditions at lines 1590/1617 of
// CGAL/boost/graph/Euler_operations.h.
template <class Mesh, class EdgeDesc>
inline bool collapse_preconditions_ok(const Mesh& m, EdgeDesc ed) {
    auto pq = halfedge(ed, m);
    auto qp = opposite(pq, m);
    const bool edge_on_border = is_border(pq, m) || is_border(qp, m);
    if (!is_border(pq, m) && CGAL::is_triangle(pq, m)) {
        auto pt = opposite(prev(pq, m), m);
        if (is_border(opposite(pt, m), m)) return false;
    }
    if (!is_border(qp, m) && CGAL::is_triangle(qp, m)) {
        auto qb = opposite(prev(qp, m), m);
        if (is_border(opposite(qb, m), m)) return false;
    }
    // Pinch guard: collapsing an INTERIOR edge whose endpoints BOTH lie on
    // the border merges two border curves through the interior -> a
    // non-manifold vertex.  Surface_mesh tolerates the pinch silently and
    // later halfedge walks corrupt or spin (measured 2026-06-12 on the
    // SAFv4 fault x bottom polylines: SIGBUS in a later collapse_edge /
    // mega-degree umbrella).  Border-edge collapses remain allowed (their
    // midpoint stays on the border curve), as do interior-edge collapses
    // with at most one border endpoint (the May fixture relies on them).
    if (!edge_on_border) {
        if (CGAL::is_border(source(pq, m), m)
            && CGAL::is_border(target(pq, m), m)) return false;
    }
    return true;
}

// Canonical edge key: two endpoint coords, integer-quantized to coord_tol
// and lex-sorted so (a,b) and (b,a) map to the same key.
struct EdgeKey {
    std::array<std::int64_t, 6> data{};
    bool operator<(const EdgeKey& o) const { return data < o.data; }
    bool operator==(const EdgeKey& o) const { return data == o.data; }
    bool operator!=(const EdgeKey& o) const { return data != o.data; }
};

template <class Point>
inline EdgeKey make_key(const Point& a, const Point& b, double coord_tol) {
    auto qz = [coord_tol](double x) -> std::int64_t {
        return static_cast<std::int64_t>(std::llround(x / coord_tol));
    };
    std::array<std::int64_t, 3> p0 = {qz(a.x()), qz(a.y()), qz(a.z())};
    std::array<std::int64_t, 3> p1 = {qz(b.x()), qz(b.y()), qz(b.z())};
    if (p1 < p0) std::swap(p0, p1);
    return EdgeKey{{p0[0], p0[1], p0[2], p1[0], p1[1], p1[2]}};
}

// Collapse polyline edges shorter than `min_edge`.  Returns Stats; the
// caller should treat n_residual_short > 0 as a hard error (the quality
// contract is violated and downstream gmsh meshing will inherit the
// short edges).
template <class Mesh, class ECM>
Stats collapse_short_polyline_edges_pairwise(
    std::vector<Mesh>& meshes,
    std::vector<ECM>&  ecms,
    double min_edge,
    double coord_tol      = 1e-6,
    unsigned max_iterations = 50)
{
    using EdgeDesc = typename boost::graph_traits<Mesh>::edge_descriptor;
    using Point    = typename Mesh::Point;

    Stats result;
    const double min_edge_sq = min_edge * min_edge;

    for (unsigned iter = 0; iter < max_iterations; ++iter) {
        ++result.n_iterations;

        // Build global registry: edge_key -> list<(mesh_idx, edge_descriptor)>.
        // Also record one representative length per key.
        std::map<EdgeKey, std::vector<std::pair<std::size_t, EdgeDesc>>> registry;
        std::map<EdgeKey, double> sq_length;
        for (std::size_t i = 0; i < meshes.size(); ++i) {
            for (auto e : meshes[i].edges()) {
                if (meshes[i].is_removed(e)) continue;
                if (!get(ecms[i], e)) continue;
                auto h = meshes[i].halfedge(e);
                const auto& a = meshes[i].point(meshes[i].source(h));
                const auto& b = meshes[i].point(meshes[i].target(h));
                EdgeKey k = make_key(a, b, coord_tol);
                registry[k].push_back({i, e});
                sq_length[k] = CGAL::to_double(CGAL::squared_distance(a, b));
            }
        }

        // Sub-floor candidates, sorted by (length asc, key) for determinism.
        std::vector<std::pair<double, EdgeKey>> candidates;
        candidates.reserve(sq_length.size());
        for (const auto& [k, L2] : sq_length) {
            if (L2 < min_edge_sq) candidates.push_back({L2, k});
        }
        if (candidates.empty()) break;
        std::sort(candidates.begin(), candidates.end(),
                  [](const auto& l, const auto& r) {
                      if (l.first != r.first) return l.first < r.first;
                      return l.second < r.second;
                  });

        // Track vertices touched by collapses this iteration; skip any
        // candidate whose endpoints touch a touched vertex (their link
        // conditions may be invalidated; defer to next iteration).
        std::map<std::pair<std::size_t, std::int64_t>, bool> touched_v;
        auto vkey = [coord_tol](std::size_t mi, const Point& p) {
            std::int64_t k = std::llround(p.x() / coord_tol) * 73856093LL
                           ^ std::llround(p.y() / coord_tol) * 19349663LL
                           ^ std::llround(p.z() / coord_tol) * 83492791LL;
            return std::pair<std::size_t, std::int64_t>{mi, k};
        };

        bool progress = false;
        for (const auto& [L2, key] : candidates) {
            const auto& edges_in_meshes = registry[key];
            if (edges_in_meshes.empty()) continue;

            // Validate every (mesh, edge) is still live and link-conditional.
            std::vector<std::pair<std::size_t, EdgeDesc>> live;
            Point midpoint{};
            bool found_midpoint = false;
            bool all_ok = true;
            bool touched = false;
            for (const auto& [i, ed] : edges_in_meshes) {
                if (meshes[i].is_removed(ed)) { all_ok = false; break; }
                if (!get(ecms[i], ed))         { all_ok = false; break; }
                auto h = meshes[i].halfedge(ed);
                auto vs = meshes[i].source(h);
                auto vt = meshes[i].target(h);
                const auto& a = meshes[i].point(vs);
                const auto& b = meshes[i].point(vt);
                if (make_key(a, b, coord_tol) != key) { all_ok = false; break; }
                if (touched_v.count(vkey(i, a)) || touched_v.count(vkey(i, b))) {
                    touched = true; break;
                }
                if (!CGAL::Euler::does_satisfy_link_condition(ed, meshes[i])) {
                    all_ok = false; break;
                }
                if (!collapse_preconditions_ok(meshes[i], ed)) {
                    ++result.n_precondition_skips;
                    all_ok = false; break;
                }
                if (!found_midpoint) {
                    midpoint = Point((a.x() + b.x()) / 2.0,
                                     (a.y() + b.y()) / 2.0,
                                     (a.z() + b.z()) / 2.0);
                    found_midpoint = true;
                }
                live.push_back({i, ed});
            }
            if (touched) continue;            // try again next iteration
            if (!all_ok || !found_midpoint) continue;
            if (live.size() != edges_in_meshes.size()) continue;

            // Same coordinate edge appearing TWICE in one mesh (duplicated
            // pinch wedges from a soup-orient load): collapsing the first
            // invalidates the second descriptor mid-loop.  Skip such
            // candidates entirely (measured 2026-06-12: stale-descriptor
            // SIGBUS via Euler_operations.h:1567 on the SAFv4 inputs).
            {
                std::vector<std::size_t> mesh_ids;
                bool dup_same_mesh = false;
                for (const auto& [i, ed] : live) {
                    if (std::find(mesh_ids.begin(), mesh_ids.end(), i)
                        != mesh_ids.end()) { dup_same_mesh = true; break; }
                    mesh_ids.push_back(i);
                }
                if (dup_same_mesh) { ++result.n_precondition_skips; continue; }
            }

            // Collapse in all containing meshes; place at midpoint.
            for (const auto& [i, ed] : live) {
                // Just-in-time revalidation: a collapse earlier in this
                // loop (other mesh) cannot invalidate `ed`, but stay
                // defensive — a removed/invalid descriptor here would
                // corrupt the mesh silently in Release builds.
                if (meshes[i].is_removed(ed)) { continue; }
                auto h = meshes[i].halfedge(ed);
                Point a = meshes[i].point(meshes[i].source(h));
                Point b = meshes[i].point(meshes[i].target(h));
                touched_v[vkey(i, a)] = true;
                touched_v[vkey(i, b)] = true;
                auto v_kept = CGAL::Euler::collapse_edge(ed, meshes[i]);
                meshes[i].point(v_kept) = midpoint;
            }
            ++result.n_collapses_done;
            progress = true;
        }

        for (auto& m : meshes) m.collect_garbage();
        if (!progress) break;
    }

    // Tally residuals.
    double rmin = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        for (auto e : meshes[i].edges()) {
            if (meshes[i].is_removed(e)) continue;
            if (!get(ecms[i], e)) continue;
            auto h = meshes[i].halfedge(e);
            const auto& a = meshes[i].point(meshes[i].source(h));
            const auto& b = meshes[i].point(meshes[i].target(h));
            double L = std::sqrt(CGAL::to_double(CGAL::squared_distance(a, b)));
            if (L < min_edge) {
                ++result.n_residual_short;
                if (L < rmin) rmin = L;
            }
        }
    }
    if (rmin != std::numeric_limits<double>::infinity()) {
        result.residual_min_length = rmin;
        result.residual_max_below  = min_edge - rmin;
    }
    return result;
}

// ---------------------------------------------------------------------------
// Phase 2: cross-polyline cluster snap.
//
// After the polyline-edge cleanup above, any remaining sub-min_edge edges
// are typically *unconstrained* edges spanning a triple junction —
// where two polyline endpoints from DIFFERENT polylines (e.g. the (0,3)
// polyline endpoint and the (0,4) polyline endpoint) lie within min_edge
// of each other.  isotropic_remeshing then meshes the surface between
// them, creating a sub-floor unconstrained edge that quality_repair
// cannot collapse (vertex_is_constrained_map blocks both endpoints).
//
// Resolution: identify clusters of polyline vertices across ALL polylines
// (regardless of which pair they came from) and snap each cluster to the
// cluster centroid in EVERY mesh containing any of the cluster vertices.
//
// Conformality preservation: cluster identity is determined purely from
// vertex coordinates (union-find with `cluster_radius`); the centroid is
// the unweighted mean of all distinct cluster coords; both are
// independent of mesh order, so the snapped position is bit-identical
// across all meshes that touch the junction.

struct ClusterStats {
    std::size_t n_iterations          = 0;
    std::size_t n_clusters_processed  = 0;
    std::size_t n_vertices_snapped    = 0;
    std::size_t n_intra_mesh_collapses = 0;
};

struct CoordKey {
    std::array<std::int64_t, 3> data{};
    bool operator<(const CoordKey& o) const { return data < o.data; }
    bool operator==(const CoordKey& o) const { return data == o.data; }
    bool operator!=(const CoordKey& o) const { return data != o.data; }
};

template <class Point>
inline CoordKey coord_key(const Point& p, double coord_tol) {
    auto qz = [coord_tol](double x) -> std::int64_t {
        return static_cast<std::int64_t>(std::llround(x / coord_tol));
    };
    return CoordKey{{qz(p.x()), qz(p.y()), qz(p.z())}};
}

template <class Point>
inline Point decode(const CoordKey& k, double coord_tol) {
    return Point(k.data[0] * coord_tol,
                 k.data[1] * coord_tol,
                 k.data[2] * coord_tol);
}

template <class Mesh, class ECM>
ClusterStats snap_cross_polyline_clusters_pairwise(
    std::vector<Mesh>& meshes,
    std::vector<ECM>&  ecms,
    double cluster_radius,
    double coord_tol      = 1e-6,
    unsigned max_iterations = 10)
{
    using VertexDesc = typename boost::graph_traits<Mesh>::vertex_descriptor;
    using EdgeDesc   = typename boost::graph_traits<Mesh>::edge_descriptor;
    using Point      = typename Mesh::Point;

    ClusterStats result;
    const double r2 = cluster_radius * cluster_radius;

    auto is_polyline_vertex = [](const Mesh& m, const ECM& ecm, VertexDesc v) {
        for (auto h : CGAL::halfedges_around_target(m.halfedge(v), m)) {
            if (m.is_removed(h)) continue;
            if (get(ecm, m.edge(h))) return true;
        }
        return false;
    };

    for (unsigned iter = 0; iter < max_iterations; ++iter) {
        ++result.n_iterations;

        // Collect all polyline vertices, keyed by coord.
        std::map<CoordKey, std::vector<std::pair<std::size_t, VertexDesc>>> reg;
        for (std::size_t i = 0; i < meshes.size(); ++i) {
            for (auto v : meshes[i].vertices()) {
                if (meshes[i].is_removed(v)) continue;
                if (!is_polyline_vertex(meshes[i], ecms[i], v)) continue;
                reg[coord_key(meshes[i].point(v), coord_tol)]
                    .push_back({i, v});
            }
        }

        // Union-find on distinct coords; unite if within cluster_radius.
        std::vector<CoordKey> keys;
        keys.reserve(reg.size());
        for (const auto& [k, _] : reg) keys.push_back(k);
        std::map<CoordKey, std::size_t> idx;
        for (std::size_t i = 0; i < keys.size(); ++i) idx[keys[i]] = i;
        std::vector<std::size_t> parent(keys.size());
        std::iota(parent.begin(), parent.end(), std::size_t{0});
        std::function<std::size_t(std::size_t)> find =
            [&](std::size_t x) {
                while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
                return x;
            };
        for (std::size_t a = 0; a < keys.size(); ++a) {
            const Point pa = decode<Point>(keys[a], coord_tol);
            for (std::size_t b = a + 1; b < keys.size(); ++b) {
                const Point pb = decode<Point>(keys[b], coord_tol);
                if (CGAL::to_double(CGAL::squared_distance(pa, pb)) <= r2) {
                    std::size_t ra = find(a), rb = find(b);
                    if (ra != rb) parent[ra] = rb;
                }
            }
        }

        // Group keys by cluster root.
        std::map<std::size_t, std::vector<CoordKey>> clusters;
        for (std::size_t i = 0; i < keys.size(); ++i) {
            clusters[find(i)].push_back(keys[i]);
        }

        bool progress = false;
        for (auto& [_, ks] : clusters) {
            if (ks.size() <= 1) continue;          // singleton — already in-cluster

            // Centroid of distinct cluster coords (unweighted, deterministic).
            double cx = 0, cy = 0, cz = 0;
            for (const auto& k : ks) {
                Point p = decode<Point>(k, coord_tol);
                cx += p.x(); cy += p.y(); cz += p.z();
            }
            const double n = static_cast<double>(ks.size());
            const Point centroid(cx / n, cy / n, cz / n);
            ++result.n_clusters_processed;

            // For each mesh, gather all polyline vertices belonging to the
            // cluster (across all coord keys it spans).
            std::map<std::size_t, std::vector<VertexDesc>> per_mesh;
            for (const auto& k : ks) {
                for (auto& [mi, v] : reg[k]) per_mesh[mi].push_back(v);
            }

            for (auto& [mi, verts] : per_mesh) {
                if (verts.size() == 1) {
                    // Single vertex in this mesh: snap to centroid.
                    if (meshes[mi].point(verts[0]) != centroid) {
                        meshes[mi].point(verts[0]) = centroid;
                        ++result.n_vertices_snapped;
                        progress = true;
                    }
                } else {
                    // Multiple vertices in this mesh must merge.  Find a
                    // direct edge between two of them and collapse it.
                    bool merged_any = true;
                    while (merged_any && verts.size() > 1) {
                        merged_any = false;
                        for (std::size_t a = 0; a < verts.size() && !merged_any; ++a) {
                            for (std::size_t b = a + 1; b < verts.size() && !merged_any; ++b) {
                                auto epair = edge(verts[a], verts[b], meshes[mi]);
                                if (!epair.second) continue;
                                EdgeDesc e = epair.first;
                                if (!CGAL::Euler::does_satisfy_link_condition(e, meshes[mi]))
                                    continue;
                                if (!collapse_preconditions_ok(meshes[mi], e))
                                    continue;
                                VertexDesc kept =
                                    CGAL::Euler::collapse_edge(e, meshes[mi]);
                                meshes[mi].point(kept) = centroid;
                                // Remove the merged-out vertex from list.
                                VertexDesc removed = (kept == verts[a]) ? verts[b] : verts[a];
                                verts.erase(std::remove(verts.begin(), verts.end(), removed),
                                            verts.end());
                                // Make sure the kept descriptor is in the list.
                                if (std::find(verts.begin(), verts.end(), kept) == verts.end()) {
                                    verts.push_back(kept);
                                }
                                ++result.n_intra_mesh_collapses;
                                merged_any = true;
                                progress = true;
                            }
                        }
                    }
                    // After loop: ideally verts.size() == 1.  If not,
                    // remaining vertices are not directly edge-connected;
                    // snap each to centroid (they are now at the same coord
                    // and may be merged later by a downstream pass).
                    for (auto v : verts) {
                        if (meshes[mi].point(v) != centroid) {
                            meshes[mi].point(v) = centroid;
                            ++result.n_vertices_snapped;
                            progress = true;
                        }
                    }
                }
            }
        }

        for (auto& m : meshes) m.collect_garbage();
        if (!progress) break;
    }
    return result;
}

}  // namespace polyline_cleanup
