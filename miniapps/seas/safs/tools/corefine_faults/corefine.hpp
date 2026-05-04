// SAFS corefine_faults — CGAL corefinement wrapper.
//
// Phase 1 of PLAN_cgal_corefine.md.  See plan §Phase 1 → Detailed
// requirements §4 (Corefine wrapper).
//
// Important — ECM ownership pattern (P-002).  The
// `edge_is_constrained_map`s on faults A/B are NOT created inside
// `corefine_pair`.  They are installed on each fault ONCE by the
// caller before the first corefine, and the same map is passed to
// every `corefine_pair` call that fault participates in.  When a
// previously-marked edge is split during a later corefine, both
// sub-edges inherit the mark (CGAL guarantee), so the mark cascades
// correctly across 3+ faults.

#ifndef SAFS_TOOLS_COREFINE_FAULTS_COREFINE_HPP
#define SAFS_TOOLS_COREFINE_FAULTS_COREFINE_HPP

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <boost/graph/graph_traits.hpp>

namespace safs::corefine {

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh   = CGAL::Surface_mesh<Kernel::Point_3>;
using EdgeDescriptor = boost::graph_traits<Mesh>::edge_descriptor;
using EdgeConstrainedMap = Mesh::Property_map<EdgeDescriptor, bool>;

struct CorefineResult {
    std::int64_t n_constrained_edges_A = 0;
    std::int64_t n_constrained_edges_B = 0;
    std::int64_t pre_n_tri_A  = 0;
    std::int64_t post_n_tri_A = 0;
    std::int64_t pre_n_tri_B  = 0;
    std::int64_t post_n_tri_B = 0;
};

// Coincident-triangle pre-flight (P-011).
// Returns the (a_face, b_face) indices of the first triangle of A
// and triangle of B that share all 3 vertex coordinates (in any
// order), or {-1, -1} if there are none.  Hash collision triggers
// an explicit 3-D coordinate comparison so there are no false
// positives.
std::pair<std::int64_t, std::int64_t>
find_coincident_triangle(const Mesh& A, const Mesh& B);

// Count constrained edges by walking every edge.
std::int64_t count_constrained_edges(const Mesh& m,
                                     const EdgeConstrainedMap& ecm);

// Wrap PMP::corefine for one unordered pair.  Throws
// std::runtime_error on self-intersection (with the offending
// fault's short name embedded in the message).
//
// ECMs are taken by const reference (R-011) to match the plan.
// CGAL's Property_map is a handle with value semantics, so the
// underlying mark storage is mesh-owned and shared regardless of
// how the parameter is passed; the const& form just makes the
// "shared across pairs (P-002)" intent unambiguous.
CorefineResult corefine_pair(
    Mesh& A, Mesh& B,
    const EdgeConstrainedMap& ecm_A,         // SHARED across pairs (P-002)
    const EdgeConstrainedMap& ecm_B,         // SHARED across pairs (P-002)
    std::string_view short_a, std::string_view short_b);

// R-004 polyline-edge-coincidence gate.  For every constrained
// edge in A, verify that an edge with the same two 3-D endpoints
// (canonical-ordered) exists as a constrained edge in B.  Returns
// true iff the constrained-edge endpoint sets agree.
//
// CGAL docs say this is a tautology of corefine; a direct check
// protects against future CGAL regressions and ECM mis-wiring.
bool polyline_edge_coincidence_gate(
    const Mesh& A, const EdgeConstrainedMap& ecm_A,
    const Mesh& B, const EdgeConstrainedMap& ecm_B);

// Diagnostic helper — exposes the same constrained-endpoint set the
// gate compares.  Used by `--verbose` failure logging to print
// asymmetric edges.
using EndpointPairDiag = std::pair<std::array<double, 3>,
                                    std::array<double, 3>>;
std::set<EndpointPairDiag> diagnose_constrained_endpoints(
    const Mesh& M, const EdgeConstrainedMap& ecm);

// ---------------------------------------------------------------------------
// Phase 2 — isotropic remeshing with constrained-edge protection
// ---------------------------------------------------------------------------

// PLAN §Phase 2 §1 — RemeshParams.
struct RemeshParams {
    double target_edge_m  = 1000.0;
    int    n_iterations   = 3;        // PMP default — sufficient per CGAL docs.
    double clearance_m    = 0.0;      // free-surface assertion floor (z ≤ -clearance_m).
    bool   allow_z_clamp  = false;    // P-008: opt-in clamp; default OFF.
};

struct RemeshResult {
    std::int64_t n_constrained_pre_split   = 0;
    std::int64_t n_constrained_post_split  = 0;
    std::int64_t n_constrained_post_remesh = 0;
    std::int64_t n_tri_pre   = 0;
    std::int64_t n_tri_post  = 0;
};

// PLAN §Phase 2 §2–3 + edge-case (boundary marking).
//
// Steps performed, in order:
//   1. Mark every border edge of M in `ecm` (idempotent — already-
//      marked edges stay marked).
//   2. `PMP::split_long_edges` over the filtered range of constrained
//      edges, using `target_edge_m` as the bisection ceiling.
//   3. Snapshot the constrained-edge count.
//   4. `PMP::isotropic_remeshing` with
//      `protect_constraints(true)` and `relax_constraints(false)` so
//      no constrained edge is split / collapsed / flipped, and no
//      constrained vertex is relaxed.
//   5. Verify the constrained-edge count is unchanged from step 3;
//      throw `std::runtime_error` on mismatch (acceptance criterion
//      "n_constrained_post_remesh == n_constrained_post_split").
//
// Throws on any protocol violation; caller maps to exit code 2.
//
// R-102 + R-105: two-ECM split.
//   `polyline_ecm` — read-only set of polyline (intersection) edges
//                    inherited from corefine.  Used to look up which
//                    edges are polyline (vs boundary) for downstream
//                    polyline-conformity verification.  Border-edge
//                    marks must NOT be added here.
//   `protect_ecm`  — mutable union (polyline ∪ boundary) used as
//                    the protect-and-don't-touch set during
//                    `split_long_edges` and `isotropic_remeshing`.
//                    Border edges are added here at step 1.
//
// CGAL's Property_map is a mesh-owned handle, so a const reference
// still permits mutation of the underlying mark storage via
// `put()`; the const-ref form just documents which handle is
// expected to stay clean and which is expected to be mutated.
//
// On entry, the implementation REQUIRES that every edge marked in
// `polyline_ecm` is also marked in `protect_ecm` (the union
// invariant).  Caller is responsible for seeding `protect_ecm`
// before calling.
RemeshResult remesh_one_fault(
    Mesh& M,
    const EdgeConstrainedMap& polyline_ecm,
    const EdgeConstrainedMap& protect_ecm,
    const RemeshParams& p);

// R-101 — symmetric collapse of short polyline edges.
//
// `PMP::corefine` produces sub-target constrained edges where the
// intersection polyline passes within ~1 m of an existing fault
// vertex.  `protect_constraints(true)` then prevents `isotropic_
// remeshing` from collapsing those edges, locking in needle
// triangles with aspect ratios up to 1000:1.
//
// This helper finds every polyline edge with length < threshold on
// either fault, deterministically picks the lex-smaller endpoint as
// the kept vertex, and calls `Euler::collapse_edge` on the matching
// edge in BOTH faults.  The matching is by endpoint coordinates
// (corefine guarantees bit-equal endpoints).  Edges where either
// side fails CGAL's link condition are skipped.
//
// Polyline conformity is preserved by construction because both
// sides keep the same lex-smaller endpoint, whose coordinates are
// bit-identical across the two meshes.  Caller should still verify
// this with `polyline_edge_coincidence_gate` post-call.
//
// Returns the number of edge pairs successfully collapsed.
int collapse_short_polyline_edges_symmetric(
    Mesh& A, const EdgeConstrainedMap& polyline_ecm_A,
    Mesh& B, const EdgeConstrainedMap& polyline_ecm_B,
    double threshold);

// R-301 — uniformize polyline edges to ~target_edge_m, symmetrically
// across both faults.
//
// Iterates two passes until stable (or `max_iter` reached):
//   (a) per-fault split: `PMP::split_long_edges` on polyline-only
//       edges, max_length = target × upper_band.  Deterministic
//       across faults because polyline edges have bit-identical
//       endpoints — both faults' splits compute the same midpoint.
//   (b) symmetric collapse: existing
//       `collapse_short_polyline_edges_symmetric` with threshold =
//       target × lower_band.
//
// After this call, every polyline edge is in approximately
// `(target × lower_band, target × upper_band]`, mean ≈ target.
// `lower_band` defaults to 0.7, `upper_band` to 1.3 — a tight ±30%
// band centred at the target.  A `runtime_error` is NOT thrown if
// `max_iter` is reached without full convergence; the caller should
// check `UniformizeResult::converged` and decide what to do (the
// non-converged case is rare and means a polyline kink is shorter
// than `lower_band × target` AND its collapse fails the link
// condition, leaving an irreducible short edge).
struct UniformizeResult {
    int  n_splits    = 0;
    int  n_collapses = 0;
    int  n_iter      = 0;
    bool converged   = false;
};
UniformizeResult uniformize_polyline_edges_symmetric(
    Mesh& A, const EdgeConstrainedMap& polyline_ecm_A,
    Mesh& B, const EdgeConstrainedMap& polyline_ecm_B,
    double target_edge_m,
    double lower_band = 0.8,
    double upper_band = 1.2,
    int    max_iter   = 16);

// R-303 — pre-corefine vertex weld.
//
// `PMP::corefine` introduces sub-target constrained edges when the
// intersection polyline pierces a parent triangle within `tol` of an
// existing vertex of either fault.  This helper, called BEFORE
// corefine, finds vertex pairs (a ∈ A, b ∈ B) within `tol` of each
// other (Euclidean distance) and snaps both to a common position
// (the midpoint).  Cross-fault conformity is preserved by
// construction because the midpoint is the same on both sides.
//
// Returns the number of vertex pairs welded.  If a vertex of A is
// within `tol` of multiple vertices of B (or vice versa) the
// closest pair is chosen and the others are skipped this round —
// callers may iterate if they want exhaustive welding.
int preweld_near_vertices(Mesh& A, Mesh& B, double tol);

// R-202: per-fault non-polyline sliver cleanup.
//
// `collapse_short_polyline_edges_symmetric` only handles short edges
// that are marked in `polyline_ecm`.  Corefine also produces NON-
// polyline sliver triangles when the polyline pierces a parent
// triangle close to one of its existing vertices: the parent gets
// split into a polyline-bounded child plus a small near-vertex
// child whose edges are NOT polyline-marked.  Those edges survive
// `collapse_short_polyline_edges_symmetric` (which filters by
// polyline mark) and survive `isotropic_remeshing` (the small edge
// MAY be a constrained-protect edge sub-divided by split_long_edges,
// or it may share a polyline endpoint and avoid collapse for that
// reason).
//
// This helper runs `PMP::remove_almost_degenerate_faces` per-fault
// with `polyline_ecm` as the constraint set, so polyline edges are
// preserved (and cross-fault conformity is preserved by extension).
// Non-polyline slivers — the residual class — get collapsed.
//
// Returns true if PMP succeeded in removing every almost-degenerate
// face (false if some were topologically un-removable, per the
// PMP::remove_almost_degenerate_faces contract).
bool cleanup_non_polyline_slivers(
    Mesh& M,
    const EdgeConstrainedMap& polyline_ecm,
    double sliver_collapse_length,
    double needle_threshold = 4.0,
    double cap_threshold_deg = 170.0);

// P-008 free-surface clearance check.
struct ClearanceReport {
    bool  any_violation       = false;
    double max_dz             = 0.0;   // = max(z(v) + clearance_m)
    Mesh::Vertex_index worst_v;
    bool  worst_v_constrained = false;
};

// Report the worst free-surface clearance violation.  A vertex is
// considered "constrained" if any incident edge is marked in `ecm`.
// Tolerance defaults to 1 mm to absorb FP round-off; `max_dz` is
// reported regardless of tolerance so the caller can log the value.
ClearanceReport check_free_surface_clearance(
    const Mesh& M,
    const EdgeConstrainedMap& ecm,
    double clearance_m,
    double tolerance_m = 1e-3);

// P-008 opt-in clamp.  Project every non-constrained vertex with
// `z > -clearance_m` down to exactly `z = -clearance_m`.
//
// Refuses (returns false, mesh unchanged) if:
//   - any constrained vertex is itself in violation (would have to
//     be moved, which would break the bit-equality of polyline
//     endpoints between fault A and fault B); OR
//   - the projection would invert any triangle (face normal flips
//     sign).  Detected by comparing pre- and post-projection face
//     normals on a copy.
//
// On success, the mesh is modified in place and the function
// returns true.
bool project_to_free_surface_safe(
    Mesh& M,
    const EdgeConstrainedMap& ecm,
    double clearance_m);

} // namespace safs::corefine

#endif // SAFS_TOOLS_COREFINE_FAULTS_COREFINE_HPP
