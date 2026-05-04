// SAFS corefine_faults — CGAL corefinement wrapper implementations.

#include "corefine.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <boost/graph/graph_traits.hpp>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace safs::corefine {

namespace {

// 64-bit hash of a triangle's three vertex coordinates, order-
// invariant.  We sort the three Point_3 values lexicographically
// then mix 9 doubles with FNV-style accumulation.
struct TriKey {
    std::array<double, 9> sorted_xyz;

    bool operator==(const TriKey& o) const {
        return sorted_xyz == o.sorted_xyz;
    }
};

// Canonicalise +0.0 vs -0.0 so the bit-cast hash agrees with the
// double-equality semantics IEEE-754 mandates (R-009).
double canon_zero(double v) noexcept {
    return v == 0.0 ? 0.0 : v;
}

struct TriKeyHash {
    std::size_t operator()(const TriKey& k) const noexcept {
        std::uint64_t h = 1469598103934665603ULL; // FNV offset
        for (double v : k.sorted_xyz) {
            v = canon_zero(v);
            std::uint64_t bits;
            static_assert(sizeof(bits) == sizeof(v));
            std::memcpy(&bits, &v, sizeof(bits));
            h ^= bits;
            h *= 1099511628211ULL;
        }
        return static_cast<std::size_t>(h);
    }
};

TriKey make_tri_key(const Kernel::Point_3& p0,
                    const Kernel::Point_3& p1,
                    const Kernel::Point_3& p2) {
    std::array<std::array<double, 3>, 3> v = {{
        {canon_zero(p0.x()), canon_zero(p0.y()), canon_zero(p0.z())},
        {canon_zero(p1.x()), canon_zero(p1.y()), canon_zero(p1.z())},
        {canon_zero(p2.x()), canon_zero(p2.y()), canon_zero(p2.z())}
    }};
    std::sort(v.begin(), v.end());
    TriKey k;
    for (int i = 0; i < 3; ++i) {
        k.sorted_xyz[3*i+0] = v[i][0];
        k.sorted_xyz[3*i+1] = v[i][1];
        k.sorted_xyz[3*i+2] = v[i][2];
    }
    return k;
}

// Walk the three target() vertices of a face's halfedges.
std::array<Kernel::Point_3, 3> face_points(
    const Mesh& m, Mesh::Face_index f) {
    auto h0 = m.halfedge(f);
    auto h1 = m.next(h0);
    auto h2 = m.next(h1);
    return {m.point(m.target(h0)),
            m.point(m.target(h1)),
            m.point(m.target(h2))};
}

} // namespace

std::pair<std::int64_t, std::int64_t>
find_coincident_triangle(const Mesh& A, const Mesh& B) {
    // Map sorted-triangle-key → (face_index, three Point_3) for A;
    // probe with B faces.  Collisions trigger an explicit Point_3
    // comparison so there are no false positives.
    struct Entry {
        std::int64_t face_idx;
        std::array<Kernel::Point_3, 3> sorted_pts;
    };
    std::unordered_map<TriKey, std::vector<Entry>, TriKeyHash> table;
    table.reserve(A.number_of_faces() * 2 + 1);

    // Canonicalise +0/-0 in the lex comparator too (R-009), so the
    // explicit collision-resolution comparison agrees with the hash.
    auto p3_lex = [](const Kernel::Point_3& p, const Kernel::Point_3& q) {
        const double px = canon_zero(p.x()), qx = canon_zero(q.x());
        if (px != qx) return px < qx;
        const double py = canon_zero(p.y()), qy = canon_zero(q.y());
        if (py != qy) return py < qy;
        return canon_zero(p.z()) < canon_zero(q.z());
    };

    for (auto f : A.faces()) {
        auto pts = face_points(A, f);
        auto sorted = pts;
        std::sort(sorted.begin(), sorted.end(), p3_lex);
        auto k = make_tri_key(pts[0], pts[1], pts[2]);
        // R-006: use the stable face index rather than an enumeration
        // counter, so the reported index remains correct under any
        // future remove_face / garbage_collect cycle.
        table[k].push_back(
            Entry{static_cast<std::int64_t>(f.idx()), sorted});
    }

    auto p3_eq_zero_aware = [](const Kernel::Point_3& p,
                                const Kernel::Point_3& q) {
        return canon_zero(p.x()) == canon_zero(q.x())
            && canon_zero(p.y()) == canon_zero(q.y())
            && canon_zero(p.z()) == canon_zero(q.z());
    };

    for (auto f : B.faces()) {
        auto pts = face_points(B, f);
        auto sorted = pts;
        std::sort(sorted.begin(), sorted.end(), p3_lex);
        auto k = make_tri_key(pts[0], pts[1], pts[2]);
        auto it = table.find(k);
        if (it != table.end()) {
            for (const auto& e : it->second) {
                bool same = true;
                for (int i = 0; i < 3; ++i) {
                    if (!p3_eq_zero_aware(e.sorted_pts[i], sorted[i])) {
                        same = false; break;
                    }
                }
                if (same) {
                    return {e.face_idx,
                            static_cast<std::int64_t>(f.idx())};
                }
            }
        }
    }

    return {-1, -1};
}

std::int64_t count_constrained_edges(const Mesh& m,
                                     const EdgeConstrainedMap& ecm) {
    std::int64_t n = 0;
    for (auto e : m.edges()) {
        if (get(ecm, e)) ++n;
    }
    return n;
}

CorefineResult corefine_pair(
    Mesh& A, Mesh& B,
    const EdgeConstrainedMap& ecm_A,
    const EdgeConstrainedMap& ecm_B,
    std::string_view short_a, std::string_view short_b) {

    CorefineResult r;
    r.pre_n_tri_A = static_cast<std::int64_t>(A.number_of_faces());
    r.pre_n_tri_B = static_cast<std::int64_t>(B.number_of_faces());

    try {
        PMP::corefine(
            A, B,
            CGAL::parameters::edge_is_constrained_map(ecm_A)
                             .throw_on_self_intersection(true),
            CGAL::parameters::edge_is_constrained_map(ecm_B));
    } catch (const PMP::Corefinement::Self_intersection_exception&) {
        std::string msg = "[corefine_faults] self-intersection detected "
                          "during corefine of pair (";
        msg += std::string(short_a);
        msg += ", ";
        msg += std::string(short_b);
        msg += "); cannot determine which fault is at fault from CGAL "
               "exception alone — re-run ts_to_stl.py with explicit "
               "self-intersection diagnostics.";
        throw std::runtime_error(msg);
    }

    r.post_n_tri_A = static_cast<std::int64_t>(A.number_of_faces());
    r.post_n_tri_B = static_cast<std::int64_t>(B.number_of_faces());
    r.n_constrained_edges_A = count_constrained_edges(A, ecm_A);
    r.n_constrained_edges_B = count_constrained_edges(B, ecm_B);
    return r;
}

namespace {

// Canonical-order endpoint pair for a constrained edge: the two
// 3-D points sorted lexicographically so the same physical edge
// produces the same key from either side.
using EndpointPair = std::pair<std::array<double, 3>,
                                std::array<double, 3>>;

// Collect canonical-ordered endpoint pairs of every constrained
// edge in M.  Border edges are skipped so the gate works
// uniformly pre- and post-remesh (post-remesh, callers pass
// protect_ecm which holds polyline ∪ boundary; border edges of A
// never match B's by construction, so filtering them isolates the
// polyline portion).
std::set<EndpointPair> collect_constrained_endpoints(
    const Mesh& M, const EdgeConstrainedMap& ecm) {
    std::set<EndpointPair> out;
    for (auto e : M.edges()) {
        if (!get(ecm, e)) continue;
        if (CGAL::is_border(e, M)) continue;
        auto h = M.halfedge(e);
        const auto& p = M.point(M.source(h));
        const auto& q = M.point(M.target(h));
        // R-205: apply canon_zero so +0/-0 endpoint variants between
        // fault A and fault B compare equal.
        std::array<double, 3> a{canon_zero(p.x()), canon_zero(p.y()),
                                 canon_zero(p.z())};
        std::array<double, 3> b{canon_zero(q.x()), canon_zero(q.y()),
                                 canon_zero(q.z())};
        if (b < a) std::swap(a, b);
        out.emplace(std::move(a), std::move(b));
    }
    return out;
}

} // namespace

bool polyline_edge_coincidence_gate(
    const Mesh& A, const EdgeConstrainedMap& ecm_A,
    const Mesh& B, const EdgeConstrainedMap& ecm_B) {
    const auto sa = collect_constrained_endpoints(A, ecm_A);
    const auto sb = collect_constrained_endpoints(B, ecm_B);
    return sa == sb;
}

std::set<EndpointPairDiag> diagnose_constrained_endpoints(
    const Mesh& M, const EdgeConstrainedMap& ecm) {
    return collect_constrained_endpoints(M, ecm);
}

// ---------------------------------------------------------------------------
// Phase 2 helpers
// ---------------------------------------------------------------------------

namespace {

// Mark every border edge of M in `ecm`.  Idempotent w.r.t. existing
// marks (already-marked edges remain marked).
void mark_border_edges_constrained(Mesh& M,
                                   const EdgeConstrainedMap& ecm) {
    for (auto e : M.edges()) {
        if (CGAL::is_border(e, M)) {
            put(ecm, e, true);
        }
    }
}

// True iff any halfedge incident on `v` is a constrained edge.
// R-104: uses CGAL's documented `halfedges_around_target` iterator
// which handles boundary vertices correctly (the manual
// `opposite(next(h))` rotation can run off the open ring on a
// boundary vertex).
bool vertex_is_constrained(const Mesh& M,
                           const EdgeConstrainedMap& ecm,
                           Mesh::Vertex_index v) {
    if (M.halfedge(v) == Mesh::null_halfedge()) return false;
    for (auto h : CGAL::halfedges_around_target(v, M)) {
        if (get(ecm, M.edge(h))) return true;
    }
    return false;
}

// Cross product magnitude (face normal length, unnormalised) for a
// triangle face — used as a sign sentinel for inversion detection.
Kernel::Vector_3 face_normal_unit(const Mesh& M, Mesh::Face_index f) {
    auto h0 = M.halfedge(f);
    auto h1 = M.next(h0);
    auto h2 = M.next(h1);
    const auto& a = M.point(M.target(h0));
    const auto& b = M.point(M.target(h1));
    const auto& c = M.point(M.target(h2));
    auto n = CGAL::cross_product(b - a, c - a);
    const double len = std::sqrt(n.squared_length());
    if (len > 0.0) {
        n = Kernel::Vector_3(n.x() / len, n.y() / len, n.z() / len);
    }
    return n;
}

} // namespace

namespace {

// Canonical-ordered endpoint pair for an edge.
using EndpointPair3 = std::pair<std::array<double, 3>,
                                 std::array<double, 3>>;

// R-205: apply canon_zero (defined further up in this file's
// anonymous namespace) so the gate / collapse path agrees with
// `make_tri_key` and `find_coincident_triangle` on the +0.0 / -0.0
// canonical form.  Without this, a polyline endpoint at z = 0 whose
// signed-zero bit differs between fault A and fault B would falsely
// fail the polyline_edge_coincidence_gate.
EndpointPair3 endpoint_pair_of(const Mesh& M,
                               Mesh::Edge_index e) {
    auto h = M.halfedge(e);
    const auto& p = M.point(M.source(h));
    const auto& q = M.point(M.target(h));
    std::array<double, 3> a{canon_zero(p.x()), canon_zero(p.y()),
                             canon_zero(p.z())};
    std::array<double, 3> b{canon_zero(q.x()), canon_zero(q.y()),
                             canon_zero(q.z())};
    if (b < a) std::swap(a, b);
    return {a, b};
}

// FNV-style hash for an EndpointPair3 (R-204).
struct EndpointPair3Hash {
    std::size_t operator()(const EndpointPair3& p) const noexcept {
        std::uint64_t h = 1469598103934665603ULL;  // FNV offset
        auto mix = [&](double v) {
            std::uint64_t bits;
            static_assert(sizeof(bits) == sizeof(v));
            std::memcpy(&bits, &v, sizeof(bits));
            h ^= bits;
            h *= 1099511628211ULL;
        };
        for (double v : p.first)  mix(v);
        for (double v : p.second) mix(v);
        return static_cast<std::size_t>(h);
    }
};

// Orient the halfedge of `e` so source(h) has the smaller-lex coords
// — i.e., source(h) corresponds to `key.first`.
Mesh::Halfedge_index halfedge_with_smaller_source(
    const Mesh& M, Mesh::Edge_index e, const EndpointPair3& key) {
    auto h = M.halfedge(e);
    const auto& p = M.point(M.source(h));
    std::array<double, 3> a{canon_zero(p.x()), canon_zero(p.y()),
                             canon_zero(p.z())};
    if (a == key.first) return h;
    return M.opposite(h);
}

double edge_length(const Mesh& M, Mesh::Edge_index e) {
    auto h = M.halfedge(e);
    const auto d = M.point(M.target(h)) - M.point(M.source(h));
    return std::sqrt(d.squared_length());
}

} // namespace

int collapse_short_polyline_edges_symmetric(
    Mesh& A, const EdgeConstrainedMap& polyline_ecm_A,
    Mesh& B, const EdgeConstrainedMap& polyline_ecm_B,
    double threshold) {
    if (threshold <= 0.0) return 0;

    // Snapshot candidate edge endpoint-pairs (with lengths) once.
    // After any collapse, edge_descriptors are invalidated, so we
    // re-find by endpoint coords each iteration.  We use the UNION
    // of A's and B's short-polyline edges since corefine guarantees
    // the polyline edges are bit-identical between the two meshes;
    // an edge present on one side should be present on the other.
    struct Cand {
        EndpointPair3 key;
        double length;
    };
    std::vector<Cand> cands;
    auto collect = [&](const Mesh& M, const EdgeConstrainedMap& pe) {
        for (auto e : M.edges()) {
            if (!get(pe, e)) continue;
            const double L = edge_length(M, e);
            if (L < threshold) cands.push_back({endpoint_pair_of(M, e), L});
        }
    };
    collect(A, polyline_ecm_A);
    collect(B, polyline_ecm_B);

    // Deduplicate candidate pairs (an edge that appears on both
    // sides has the same canonical key).
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) {
                  if (a.key != b.key) return a.key < b.key;
                  return a.length < b.length;
              });
    cands.erase(std::unique(cands.begin(), cands.end(),
                            [](const Cand& a, const Cand& b) {
                                return a.key == b.key;
                            }),
                cands.end());

    // Process shortest first so longer edges adjacent to a
    // just-collapsed vertex see the new geometry.
    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) {
                  return a.length < b.length;
              });

    // R-204: O(1) endpoint→edge index per fault, rebuilt after each
    // successful collapse so stale descriptors do not survive.
    using Index = std::unordered_map<EndpointPair3, Mesh::Edge_index,
                                      EndpointPair3Hash>;
    auto build_index = [](const Mesh& M) {
        Index idx;
        idx.reserve(M.number_of_edges() * 2 + 1);
        for (auto e : M.edges()) idx.emplace(endpoint_pair_of(M, e), e);
        return idx;
    };
    auto lookup = [](const Mesh& M, const Index& idx,
                     const EndpointPair3& key) -> Mesh::Edge_index {
        auto it = idx.find(key);
        if (it == idx.end()) return Mesh::null_edge();
        // Verify the indexed descriptor is still alive AND still has
        // the same endpoints — defensive against any stale entry that
        // could remain if a previous collapse path didn't rebuild.
        if (endpoint_pair_of(M, it->second) != key) return Mesh::null_edge();
        return it->second;
    };
    Index idx_A = build_index(A);
    Index idx_B = build_index(B);

    int collapsed = 0;
    for (const auto& c : cands) {
        auto eA = lookup(A, idx_A, c.key);
        auto eB = lookup(B, idx_B, c.key);
        if (eA == Mesh::null_edge() || eB == Mesh::null_edge()) {
            // The edge no longer exists on at least one side (a
            // previous collapse merged it in).  Skip silently — the
            // polyline conformity check post-cleanup will surface any
            // asymmetry.
            continue;
        }
        if (!CGAL::Euler::does_satisfy_link_condition(eA, A)
            || !CGAL::Euler::does_satisfy_link_condition(eB, B)) {
            continue;
        }
        // Skip if either side's matching edge is no longer marked
        // constrained (e.g., subdivided by an earlier collapse + the
        // mark wasn't propagated).
        if (!get(polyline_ecm_A, eA) || !get(polyline_ecm_B, eB)) {
            continue;
        }
        // Verify lengths are still small — a previous collapse may
        // have moved one endpoint and lengthened this edge.
        if (edge_length(A, eA) >= threshold
            && edge_length(B, eB) >= threshold) {
            continue;
        }
        // Orient halfedges so source = the lex-smaller endpoint
        // (deterministically the same on both sides).
        auto hA = halfedge_with_smaller_source(A, eA, c.key);
        auto hB = halfedge_with_smaller_source(B, eB, c.key);
        try {
            CGAL::Euler::collapse_edge(A.edge(hA), A);
            CGAL::Euler::collapse_edge(B.edge(hB), B);
            ++collapsed;
            // R-204: rebuild indices so subsequent lookups don't
            // return stale descriptors.  Total cost over K collapses
            // is O(K · E), still O(E) amortised per collapse vs the
            // previous O(E) per candidate.
            idx_A = build_index(A);
            idx_B = build_index(B);
        } catch (...) {
            // Defensive: any CGAL collapse failure aborts this pair
            // but does not invalidate previously-collapsed pairs.
            continue;
        }
    }
    return collapsed;
}

UniformizeResult uniformize_polyline_edges_symmetric(
    Mesh& A, const EdgeConstrainedMap& polyline_ecm_A,
    Mesh& B, const EdgeConstrainedMap& polyline_ecm_B,
    double target_edge_m,
    double lower_band,
    double upper_band,
    int    max_iter) {
    namespace PMP_local = CGAL::Polygon_mesh_processing;
    namespace params    = CGAL::parameters;

    UniformizeResult r;
    if (target_edge_m <= 0.0 || lower_band <= 0.0 || upper_band <= lower_band) {
        return r;
    }
    const double L_max = target_edge_m * upper_band;
    const double L_min = target_edge_m * lower_band;

    auto split_long_in = [&](Mesh& M, const EdgeConstrainedMap& ecm) {
        std::vector<EdgeDescriptor> long_edges;
        long_edges.reserve(M.number_of_edges());
        for (auto e : M.edges()) {
            if (!get(ecm, e)) continue;
            if (edge_length(M, e) > L_max) long_edges.push_back(e);
        }
        const int n = static_cast<int>(long_edges.size());
        if (n > 0) {
            PMP_local::split_long_edges(
                long_edges, L_max, M,
                params::edge_is_constrained_map(ecm));
        }
        return n;
    };

    for (int it = 0; it < max_iter; ++it) {
        r.n_iter = it + 1;

        // (a) per-fault split — deterministic across A and B because
        // polyline edges have bit-identical endpoints, so midpoint
        // bisection produces bit-identical sub-edges.
        const int s_a = split_long_in(A, polyline_ecm_A);
        const int s_b = split_long_in(B, polyline_ecm_B);
        r.n_splits += s_a + s_b;

        // (b) symmetric collapse of polyline edges still below L_min.
        const int c = collapse_short_polyline_edges_symmetric(
            A, polyline_ecm_A, B, polyline_ecm_B, L_min);
        r.n_collapses += c;

        if (s_a == 0 && s_b == 0 && c == 0) {
            r.converged = true;
            break;
        }
    }
    return r;
}

namespace {

// KD-tree-style 3-D bucketed lookup for finding near-vertex pairs
// without a full O(N×M) scan.  We bucket each vertex into a uniform
// 3-D grid with cell size = `tol`; querying a vertex against the
// 27 neighbouring cells is O(1) on average.
struct GridKey {
    std::int64_t i, j, k;
    bool operator==(const GridKey& o) const {
        return i == o.i && j == o.j && k == o.k;
    }
};
struct GridKeyHash {
    std::size_t operator()(const GridKey& g) const noexcept {
        std::uint64_t h = static_cast<std::uint64_t>(g.i) * 1469598103934665603ULL;
        h ^= static_cast<std::uint64_t>(g.j); h *= 1099511628211ULL;
        h ^= static_cast<std::uint64_t>(g.k); h *= 1099511628211ULL;
        return static_cast<std::size_t>(h);
    }
};

GridKey grid_of(const Kernel::Point_3& p, double cell) {
    return {static_cast<std::int64_t>(std::floor(p.x() / cell)),
            static_cast<std::int64_t>(std::floor(p.y() / cell)),
            static_cast<std::int64_t>(std::floor(p.z() / cell))};
}

} // namespace

int preweld_near_vertices(Mesh& A, Mesh& B, double tol) {
    if (tol <= 0.0) return 0;
    const double tol_sq = tol * tol;

    // Bucket B's vertices.
    std::unordered_map<GridKey, std::vector<Mesh::Vertex_index>, GridKeyHash>
        buckets_B;
    buckets_B.reserve(B.number_of_vertices() * 2);
    for (auto v : B.vertices()) {
        buckets_B[grid_of(B.point(v), tol)].push_back(v);
    }

    // For each vertex of A, find the closest vertex of B within `tol`.
    // If both are still un-welded, snap them to the midpoint.
    int welded = 0;
    std::vector<bool> b_used(B.number_of_vertices(), false);

    for (auto va : A.vertices()) {
        const auto pa = A.point(va);
        const auto k = grid_of(pa, tol);

        Mesh::Vertex_index best_b = Mesh::null_vertex();
        double best_d2 = tol_sq;

        for (std::int64_t di = -1; di <= 1; ++di)
        for (std::int64_t dj = -1; dj <= 1; ++dj)
        for (std::int64_t dk = -1; dk <= 1; ++dk) {
            const auto it = buckets_B.find({k.i + di, k.j + dj, k.k + dk});
            if (it == buckets_B.end()) continue;
            for (auto vb : it->second) {
                if (b_used[static_cast<std::size_t>(vb.idx())]) continue;
                const auto pb = B.point(vb);
                const double dx = pa.x() - pb.x();
                const double dy = pa.y() - pb.y();
                const double dz = pa.z() - pb.z();
                const double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_b = vb;
                }
            }
        }

        if (best_b == Mesh::null_vertex()) continue;
        // Snap both to the midpoint — bit-identical 3-D position
        // preserves cross-fault conformity by construction.
        const auto pb = B.point(best_b);
        const Kernel::Point_3 mid(
            0.5 * (pa.x() + pb.x()),
            0.5 * (pa.y() + pb.y()),
            0.5 * (pa.z() + pb.z()));
        A.point(va)    = mid;
        B.point(best_b) = mid;
        b_used[static_cast<std::size_t>(best_b.idx())] = true;
        ++welded;
    }
    return welded;
}

bool cleanup_non_polyline_slivers(
    Mesh& M,
    const EdgeConstrainedMap& polyline_ecm,
    double sliver_collapse_length,
    double needle_threshold,
    double cap_threshold_deg) {
    namespace PMP_local = CGAL::Polygon_mesh_processing;
    namespace params    = CGAL::parameters;

    // CGAL's PMP::remove_almost_degenerate_faces expects:
    //   cap_threshold       — cosine of the apex angle above which a
    //                         face is a "cap"; range [-1, 0].  At
    //                         170°: cos(170°) ≈ -0.9848.
    //   needle_threshold    — max(edge)/min(edge) ratio above which a
    //                         face is a "needle".  Default 4.
    //   collapse_length_threshold — max length of edge to collapse;
    //                               0 = no upper bound.
    const double cap_cos =
        std::cos(cap_threshold_deg * 3.14159265358979323846 / 180.0);

    // Cross-fault polyline conformity: every vertex incident to a
    // polyline edge must remain in place across the cleanup.  Without
    // a vertex_is_constrained_map, CGAL may collapse a non-polyline
    // sliver edge incident to a polyline endpoint and reposition that
    // endpoint, breaking bit-equality with fault B's polyline.
    using VertexDescriptor =
        boost::graph_traits<Mesh>::vertex_descriptor;
    auto vcm_pm = M.add_property_map<VertexDescriptor, bool>(
        "v:cleanup_polyline_endpoint", false);
    if (!vcm_pm.second) {
        // Already exists from a prior call — defensive cleanup
        // attempts on the same fault are not expected.
        M.remove_property_map(vcm_pm.first);
        vcm_pm = M.add_property_map<VertexDescriptor, bool>(
            "v:cleanup_polyline_endpoint", false);
    }
    auto vcm = vcm_pm.first;
    for (auto e : M.edges()) {
        if (!get(polyline_ecm, e)) continue;
        auto h = M.halfedge(e);
        put(vcm, M.source(h), true);
        put(vcm, M.target(h), true);
    }

    const bool ok = PMP_local::remove_almost_degenerate_faces(
        M,
        params::cap_threshold(cap_cos)
               .needle_threshold(needle_threshold)
               .collapse_length_threshold(sliver_collapse_length)
               .edge_is_constrained_map(polyline_ecm)
               .vertex_is_constrained_map(vcm));

    M.remove_property_map(vcm);
    return ok;
}

RemeshResult remesh_one_fault(
    Mesh& M,
    const EdgeConstrainedMap& polyline_ecm,
    const EdgeConstrainedMap& protect_ecm,
    const RemeshParams& p) {
    namespace PMP_local = CGAL::Polygon_mesh_processing;
    namespace params    = CGAL::parameters;

    if (p.target_edge_m <= 0.0) {
        throw std::runtime_error(
            "[corefine_faults] remesh_one_fault: target_edge_m must "
            "be > 0 (callers requesting --no-remesh should skip this "
            "function entirely)");
    }
    if (p.n_iterations < 0) {
        throw std::runtime_error(
            "[corefine_faults] remesh_one_fault: n_iterations must "
            "be >= 0");
    }
    // Verify the union invariant: every polyline edge must already be
    // marked in protect_ecm before this function runs.
    for (auto e : M.edges()) {
        if (get(polyline_ecm, e) && !get(protect_ecm, e)) {
            throw std::runtime_error(
                "[corefine_faults] remesh_one_fault: caller failed to "
                "seed protect_ecm with polyline_ecm marks (R-105 union "
                "invariant violated)");
        }
    }

    RemeshResult r;
    r.n_tri_pre = static_cast<std::int64_t>(M.number_of_faces());
    r.n_constrained_pre_split = count_constrained_edges(M, protect_ecm);

    // Step 1 — mark border edges in protect_ecm only (R-105: must NOT
    // pollute polyline_ecm with border-edge marks).
    mark_border_edges_constrained(M, protect_ecm);

    // Step 2 — split long constrained edges so the remesher does not
    // lock in slivers along oversized border edges.  Polyline edges
    // are expected to have been uniformized via
    // `uniformize_polyline_edges_symmetric` (R-301) BEFORE this step,
    // so this split pass acts mostly on borders.
    //
    // R-302: split ceiling reverted to `target_edge_m` (was target × 4/3
    // in R-201).  The 4/3 multiplier widened the polyline-edge spread,
    // doubling the worst-case aspect on SBMT-SAF (7.6 → 15.19) when
    // short polyline edges escaped the collapse pass.  Now that R-301
    // drives polyline edges into [target × 0.7, target × 1.3] BEFORE
    // remesh, the remaining splits inside this function only affect
    // border edges; ceiling = target keeps surrounding non-constrained
    // edges close to target, preserving aspect-ratio quality.
    const double split_ceiling = p.target_edge_m;
    {
        std::vector<EdgeDescriptor> long_constrained;
        long_constrained.reserve(M.number_of_edges());
        for (auto e : M.edges()) {
            if (get(protect_ecm, e)) long_constrained.push_back(e);
        }
        PMP_local::split_long_edges(
            long_constrained, split_ceiling, M,
            params::edge_is_constrained_map(protect_ecm));
    }

    // Step 3 — snapshot constrained-edge count before remeshing.
    r.n_constrained_post_split = count_constrained_edges(M, protect_ecm);

    // Step 4 — isotropic remeshing with constraints fully protected.
    if (p.n_iterations > 0) {
        PMP_local::isotropic_remeshing(
            faces(M), p.target_edge_m, M,
            params::edge_is_constrained_map(protect_ecm)
                   .number_of_iterations(p.n_iterations)
                   .protect_constraints(true)
                   .relax_constraints(false));
    }

    // Step 5 — verify constraints survived.
    r.n_constrained_post_remesh = count_constrained_edges(M, protect_ecm);
    if (r.n_constrained_post_remesh != r.n_constrained_post_split) {
        throw std::runtime_error(
            "[corefine_faults] remesh_one_fault: constrained-edge "
            "count changed during isotropic_remeshing despite "
            "protect_constraints(true). Was "
            + std::to_string(r.n_constrained_post_split)
            + " before, now "
            + std::to_string(r.n_constrained_post_remesh)
            + ". This is a CGAL invariant violation; see Phase 2 "
              "edge-case §1.");
    }

    r.n_tri_post = static_cast<std::int64_t>(M.number_of_faces());
    return r;
}

ClearanceReport check_free_surface_clearance(
    const Mesh& M,
    const EdgeConstrainedMap& ecm,
    double clearance_m,
    double tolerance_m) {
    ClearanceReport rep;
    rep.worst_v = Mesh::null_vertex();
    for (auto v : M.vertices()) {
        const double dz = M.point(v).z() + clearance_m;  // > 0 means above the plane
        if (dz > rep.max_dz) {
            rep.max_dz = dz;
            rep.worst_v = v;
            rep.worst_v_constrained = vertex_is_constrained(M, ecm, v);
        }
    }
    rep.any_violation = (rep.max_dz > tolerance_m);
    return rep;
}

bool project_to_free_surface_safe(
    Mesh& M,
    const EdgeConstrainedMap& ecm,
    double clearance_m) {
    // First, refuse if any constrained vertex would be displaced.
    for (auto v : M.vertices()) {
        if (M.point(v).z() > -clearance_m
            && vertex_is_constrained(M, ecm, v)) {
            return false;
        }
    }

    // Inversion-safety check: take a copy, apply the projection
    // there, compare per-face normals.  Only commit if no flips.
    Mesh trial = M;
    for (auto v : trial.vertices()) {
        const auto& q = trial.point(v);
        if (q.z() > -clearance_m) {
            trial.point(v) = Kernel::Point_3(q.x(), q.y(), -clearance_m);
        }
    }

    // Per-face normal-sign comparison.  Use unit normals so length
    // changes don't affect the dot product sign.
    //
    // R-103: distinguish input degeneracy (skip — not the clamp's
    // fault) from output degeneracy (refuse — the clamp would
    // silently emit a zero-area triangle).
    for (auto f_pre : M.faces()) {
        Mesh::Face_index f_post(f_pre);
        const auto n_pre  = face_normal_unit(M, f_pre);
        const auto n_post = face_normal_unit(trial, f_post);
        if (n_pre.squared_length() == 0.0) {
            continue;            // input was degenerate; not the clamp's fault
        }
        if (n_post.squared_length() == 0.0) {
            return false;        // R-103: clamp would create a degenerate face
        }
        if (CGAL::scalar_product(n_pre, n_post) <= 0.0) {
            return false;  // projection inverts at least one triangle
        }
    }

    // Commit.
    M = std::move(trial);
    return true;
}

} // namespace safs::corefine
