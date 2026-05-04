// SAFS corefine_faults — autorefine-mode implementation.
// See autorefine_mode.hpp for design rationale.
//
// Two implementations selected at compile time by SAFS_HAVE_AUTOREFINE_SNAP
// (set in tools/CMakeLists.txt from the detected CGAL version):
//
//   - SAFS_HAVE_AUTOREFINE_SNAP=1 (CGAL 6.1+):
//       Use `PMP::autorefine_triangle_soup()` with
//       `apply_iterative_snap_rounding(true)`.  Intersection points are
//       rounded to the double-precision grid, eliminating sub-meter
//       split artifacts at the source.  A visitor tracks parent face
//       provenance so the per-fault demux is exact (no point-in-
//       triangle heuristic).
//
//   - SAFS_HAVE_AUTOREFINE_SNAP=0 (CGAL 5.6 legacy):
//       Use `PMP::experimental::autorefine` on a Surface_mesh, then
//       point-in-triangle face-tag assignment for the demux.  Retained
//       so the tool still builds against the previous tested baseline.

#include "autorefine_mode.hpp"
#include "io.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#if defined(SAFS_HAVE_AUTOREFINE_SNAP) && SAFS_HAVE_AUTOREFINE_SNAP
  #include <CGAL/Polygon_mesh_processing/autorefinement.h>
#else
  #include <CGAL/Polygon_mesh_processing/corefinement.h>
  #include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
  #include <CGAL/Polygon_mesh_processing/self_intersections.h>
#endif

// Fix B (sliver-fix-B, 2026-04-30): per-fault short-edge cleanup.
// Used by both legacy and snap paths.  See REVIEW_sliver_classification.md
// and the 2026-04-30 debug session for the diagnostic that motivated
// this — population B is the ~46 short-edge fault triangle slivers
// (incl. the absolute-worst γ=1.36e-04) caused by autorefine landing
// intersection points within ~200 m of existing edge endpoints.
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/boost/graph/helpers.h>

// Fix P (sliver-fix-P, 2026-04-30): pre-autorefine per-fault isotropic
// remeshing.  Each fault is independently smoothed to a uniform target
// edge length BEFORE any cross-fault polyline exists, so the fault
// outer perimeter is auto-constrained (it is the only border of the
// stand-alone patch) and the interior is regularised.  Cleaner inputs
// produce fewer near-vertex intersection-point artifacts in the
// downstream global autorefine call.  Safe because no cross-fault
// constraint exists at this stage.
#include <CGAL/Polygon_mesh_processing/remesh.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace safs::corefine {

namespace {

// Build an oriented polygon soup from a single STL.  No corefine, just
// vertex-dedup + orient.  The returned points/polygons can be appended
// to a global-soup buffer.
bool stl_to_soup(const std::filesystem::path& p,
                 std::vector<Kernel::Point_3>& out_points,
                 std::vector<std::vector<std::size_t>>& out_polygons) {
    std::vector<Kernel::Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;
    if (!CGAL::IO::read_polygon_soup(p.string(), points, polygons)) {
        return false;
    }
    PMP::repair_polygon_soup(
        points, polygons,
        CGAL::parameters::erase_all_duplicates(true)
                         .require_same_orientation(true));
    PMP::orient_polygon_soup(points, polygons);
    out_points  = std::move(points);
    out_polygons = std::move(polygons);
    return true;
}

// Fix P helper: regularise a single fault's polygon soup to a uniform
// target edge length, IN PLACE.  Safe at this point in the pipeline
// because the fault is a stand-alone patch — no cross-fault polyline
// exists yet, so isotropic_remeshing's automatic border-edge
// constraint protects exactly the fault outer perimeter (where the
// fault meets the box top, box bottom, and fault tips) without any
// risk of breaking a future cross-fault PLC.
//
// `target_edge_m` is the requested edge length; values <= 0 disable
// the pass (caller's choice).  `n_iters` defaults to 3 (the CGAL
// `isotropic_remeshing` documented sweet-spot for typical surfaces).
//
// On any exception (e.g., non-manifold soup that polygon_soup_to_polygon_mesh
// would refuse), logs a warning and leaves (pts, polys) untouched —
// the downstream autorefine then operates on the original unsmoothed
// soup, identical to the pre-Fix-P behaviour.
void fixP_pre_remesh_soup(
        std::vector<Kernel::Point_3>& pts,
        std::vector<std::vector<std::size_t>>& polys,
        const std::string& short_name,
        double target_edge_m,
        int    n_iters) {
    if (target_edge_m <= 0.0) return;
    if (pts.empty() || polys.empty()) return;

    // 1. Build a stand-alone Surface_mesh from the soup.  If the soup
    //    is not a valid manifold mesh, skip (defensive — CFM input
    //    is repaired by stl_to_soup, but we don't want a hard crash
    //    if a future CFM addition violates the assumption).
    if (!PMP::is_polygon_soup_a_polygon_mesh(polys)) {
        std::cerr << "[autorefine_faults] Fix P skip " << short_name
                  << ": polygon soup is not a manifold mesh; "
                     "leaving input untouched\n";
        return;
    }
    Mesh m_pre;
    PMP::polygon_soup_to_polygon_mesh(pts, polys, m_pre);
    if (m_pre.number_of_faces() == 0) {
        std::cerr << "[autorefine_faults] Fix P skip " << short_name
                  << ": empty mesh after soup conversion\n";
        return;
    }

    // 2. Run isotropic_remeshing.  Per CGAL docs, `protect_constraints
    //    (true)` together with the auto-constrained patch boundary
    //    edges keeps the fault outer perimeter exactly fixed.
    //    Interior vertices are smoothed and re-distributed to ~target.
    const std::int64_t n_before = static_cast<std::int64_t>(
        m_pre.number_of_faces());
    try {
        PMP::isotropic_remeshing(
            m_pre.faces(),
            target_edge_m,
            m_pre,
            CGAL::parameters::number_of_iterations(n_iters)
                             .protect_constraints(true));
    } catch (const std::exception& ex) {
        std::cerr << "[autorefine_faults] Fix P FAILED on " << short_name
                  << ": " << ex.what()
                  << "; leaving input untouched\n";
        return;
    }
    const std::int64_t n_after = static_cast<std::int64_t>(
        m_pre.number_of_faces());

    // 3. Read back into soup form (replace pts/polys).  Mesh face
    //    iteration produces the post-remesh triangulation; mesh
    //    vertex iteration produces the post-remesh vertex set.  Use
    //    a vertex-index → soup-index map to rebuild the index-based
    //    polygon list.
    std::unordered_map<typename Mesh::Vertex_index, std::size_t> v_to_idx;
    v_to_idx.reserve(m_pre.number_of_vertices());
    std::vector<Kernel::Point_3> new_pts;
    new_pts.reserve(m_pre.number_of_vertices());
    for (auto v : m_pre.vertices()) {
        v_to_idx.emplace(v, new_pts.size());
        new_pts.push_back(m_pre.point(v));
    }
    std::vector<std::vector<std::size_t>> new_polys;
    new_polys.reserve(m_pre.number_of_faces());
    for (auto f : m_pre.faces()) {
        auto h0 = m_pre.halfedge(f);
        auto h1 = m_pre.next(h0);
        auto h2 = m_pre.next(h1);
        new_polys.push_back({
            v_to_idx[m_pre.target(h0)],
            v_to_idx[m_pre.target(h1)],
            v_to_idx[m_pre.target(h2)]});
    }
    pts = std::move(new_pts);
    polys = std::move(new_polys);

    std::cerr << "[autorefine_faults] Fix P " << short_name
              << ": " << n_before << " -> " << n_after
              << " faces (target_edge_m=" << target_edge_m
              << " m, iters=" << n_iters << ")\n";
}

#if defined(SAFS_HAVE_AUTOREFINE_SNAP) && SAFS_HAVE_AUTOREFINE_SNAP

// Tracking visitor for autorefine_triangle_soup (CGAL 6.1+).
//
// CGAL hooks:
//   number_of_output_triangles(N): pre-allocate parent slot for N
//       output triangles.  N is an upper bound; final size may be
//       smaller after degenerate removal.
//   verbatim_triangle_copy(tgt, src): output triangle `tgt` is a
//       direct copy of input triangle `src` (no subdivision).
//   new_subtriangle(tgt, src): output triangle `tgt` is a NEW piece
//       from subdividing input triangle `src`.
//   delete_triangle(src): input triangle `src` was completely removed
//       (degenerate or fully covered).  We do not need to track this
//       because parent[tgt] is only set via copy/new_subtriangle.
//
// `parent[tgt_id]` returns the input triangle index that this output
// triangle descends from.  Used by the demux to assign each output
// triangle to its source fault tag exactly (no centroid heuristic).
struct ParentTrackingVisitor {
    // Shared storage so the visitor remains lightweight on copy
    // (CGAL's docs say "The visitor will be copied" — copies must
    // observe the same parent vector).
    std::shared_ptr<std::vector<std::int64_t>> parent;

    ParentTrackingVisitor()
        : parent(std::make_shared<std::vector<std::int64_t>>()) {}

    void number_of_output_triangles(std::size_t nbt) {
        parent->assign(nbt, -1);
    }
    void verbatim_triangle_copy(std::size_t tgt_id, std::size_t src_id) {
        if (tgt_id >= parent->size()) parent->resize(tgt_id + 1, -1);
        (*parent)[tgt_id] = static_cast<std::int64_t>(src_id);
    }
    void new_subtriangle(std::size_t tgt_id, std::size_t src_id) {
        if (tgt_id >= parent->size()) parent->resize(tgt_id + 1, -1);
        (*parent)[tgt_id] = static_cast<std::int64_t>(src_id);
    }
    void delete_triangle(std::size_t /*src_id*/) {}
};

#endif // SAFS_HAVE_AUTOREFINE_SNAP

// Fix B helper — collapse autorefine-artifact short edges in a single
// per-fault Surface_mesh.  Constrains BORDER edges (polyline + outer
// perimeter + surface trace) AND BORDER VERTICES so cross-fault
// polyline endpoints stay bit-identical between fault A and fault B's
// submeshes.  This vertex-side constraint is the difference from the
// earlier failed attempt at remove_almost_degenerate_faces — without
// `vertex_is_constrained_map`, CGAL is free to relocate constrained-
// edge endpoints during cap flips, which broke HXT's PLC recovery on
// real CFM data.
//
// CGAL 5.6 docs (repair_degeneracies.h:579) guarantee:
//   "A constrained vertex is guaranteed to be present in tmesh after
//    the function call."
// — i.e., not collapsed AWAY and (per the implementation) not relocated.
//
// Returns the number of faces removed (≥0).  On exception, logs and
// returns -1 leaving the mesh in the partial state that triggered the
// throw — the caller proceeds with whatever survived.
std::int64_t fixB_collapse_short_edges(Mesh& mesh,
                                        const std::string& short_name,
                                        double collapse_length_threshold,
                                        double needle_threshold) {
    if (mesh.number_of_faces() == 0) return 0;

    using EdgeDescr   = boost::graph_traits<Mesh>::edge_descriptor;
    using VertexDescr = boost::graph_traits<Mesh>::vertex_descriptor;
    auto ecm_pair = mesh.template add_property_map<EdgeDescr, bool>(
        "e:fixB_constrained", false);
    auto vcm_pair = mesh.template add_property_map<VertexDescr, bool>(
        "v:fixB_constrained", false);
    auto ecm = ecm_pair.first;
    auto vcm = vcm_pair.first;

    std::int64_t n_constrained_edges = 0;
    for (auto e : mesh.edges()) {
        if (CGAL::is_border(e, mesh)) {
            put(ecm, e, true);
            ++n_constrained_edges;
            auto h = mesh.halfedge(e);
            put(vcm, mesh.source(h), true);
            put(vcm, mesh.target(h), true);
        }
    }
    std::int64_t n_constrained_verts = 0;
    for (auto v : mesh.vertices()) {
        if (get(vcm, v)) ++n_constrained_verts;
    }

    const std::int64_t n_before = static_cast<std::int64_t>(
        mesh.number_of_faces());
    std::int64_t removed = -1;
    try {
        // cap_threshold = -1.0 disables cap detection (no edge flips).
        // Cap flips change WHICH triangle is incident to a polyline-
        // adjacent edge — that re-orientation cascades through
        // gmsh's STL-merge step and breaks cross-fault PLC, even
        // though no vertex moves (HXT then rejects with
        // "Segment and Facet intersect").  Restricting Fix B to
        // pure needle collapse is the conservative path; caps
        // remain in the output but they're not the primary sliver
        // driver.
        PMP::remove_almost_degenerate_faces(
            mesh,
            CGAL::parameters::edge_is_constrained_map(ecm)
                             .vertex_is_constrained_map(vcm)
                             .cap_threshold(-1.0)
                             .needle_threshold(needle_threshold)
                             .collapse_length_threshold(
                                 collapse_length_threshold));
        const std::int64_t n_after = static_cast<std::int64_t>(
            mesh.number_of_faces());
        removed = n_before - n_after;
        std::cerr << "[autorefine_faults] Fix B " << short_name
                  << ": " << n_constrained_edges << " border edges, "
                  << n_constrained_verts << " border verts pinned; "
                  << removed << " face(s) removed ("
                  << n_before << " -> " << n_after << ")\n";
    } catch (const std::exception& ex) {
        std::cerr << "[autorefine_faults] WARN: Fix B failed on "
                  << short_name << ": " << ex.what()
                  << "; continuing with un-cleaned submesh.\n";
        removed = -1;
    }

    mesh.remove_property_map(ecm);
    mesh.remove_property_map(vcm);
    return removed;
}

} // namespace

#if defined(SAFS_HAVE_AUTOREFINE_SNAP) && SAFS_HAVE_AUTOREFINE_SNAP

// CGAL 6.1+ path: snap-rounded autorefine_triangle_soup with exact
// parent tracking via visitor.

AutorefineResult autorefine_faults(
    const std::vector<std::filesystem::path>& in_paths,
    const std::vector<std::filesystem::path>& out_paths,
    const std::vector<std::string>& shorts,
    double target_edge_m,
    int    remesh_iters) {
    // target_edge_m / remesh_iters now used by Fix P pre-remesh.

    if (in_paths.size() != out_paths.size()
        || in_paths.size() != shorts.size()) {
        throw std::runtime_error(
            "[autorefine_faults] in/out/short vector size mismatch");
    }
    const std::size_t N = in_paths.size();
    if (N == 0) {
        throw std::runtime_error("[autorefine_faults] no input faults");
    }

    AutorefineResult res;
    res.per_fault_n_input.resize(N);
    res.per_fault_n_output.resize(N);

    // Step 1: load each fault as polygon soup; build a global soup with
    // per-input-triangle fault tags.
    using Triangle = std::array<std::size_t, 3>;
    std::vector<Kernel::Point_3> all_points;
    std::vector<Triangle>         all_tris;
    std::vector<std::int64_t>     tri_to_fault;

    for (std::size_t k = 0; k < N; ++k) {
        std::vector<Kernel::Point_3> pts;
        std::vector<std::vector<std::size_t>> polys;
        if (!stl_to_soup(in_paths[k], pts, polys)) {
            throw std::runtime_error(
                "[autorefine_faults] failed to read soup from "
                + in_paths[k].string());
        }
        // Fix P: pre-autorefine isotropic remeshing on this single
        // fault.  Safe — no cross-fault polyline exists yet.
        fixP_pre_remesh_soup(pts, polys, shorts[k],
                              target_edge_m, remesh_iters);
        const std::size_t off = all_points.size();
        for (auto& p : pts) all_points.push_back(p);
        std::int64_t kept = 0;
        for (auto& poly : polys) {
            if (poly.size() != 3) continue;  // only triangles
            all_tris.push_back(Triangle{
                poly[0] + off, poly[1] + off, poly[2] + off});
            tri_to_fault.push_back(static_cast<std::int64_t>(k));
            ++kept;
        }
        res.per_fault_n_input[k] = kept;
    }
    res.n_input_faces = static_cast<std::int64_t>(all_tris.size());

    // Step 2: snap-rounded autorefine on the soup.
    //
    // apply_iterative_snap_rounding(true) rounds intersection points
    // onto the double grid; sub-meter post-split edges are eliminated
    // at the source.  snap_grid_size=23 (default) gives ~1mm resolution
    // on the [-2^23, 2^23] m scale; for our domain (~150 km) that's
    // fine.
    //
    // The visitor records parent-input-triangle index for every output
    // triangle.  This replaces the old point-in-triangle heuristic and
    // is exact by construction.
    ParentTrackingVisitor vis;
    const bool snap_success = PMP::autorefine_triangle_soup(
        all_points, all_tris,
        CGAL::parameters::apply_iterative_snap_rounding(true)
                         .number_of_iterations(5)
                         .snap_grid_size(23)
                         .visitor(vis));
    res.had_self_intersections = true;
    if (!snap_success) {
        std::cerr << "[autorefine_faults] WARN: snap-rounded autorefine "
                     "did not converge within 5 iterations; output may "
                     "still contain self-intersections.\n";
    }

    res.n_output_faces = static_cast<std::int64_t>(all_tris.size());
    res.n_intersections_resolved =
        std::max<std::int64_t>(0, res.n_output_faces - res.n_input_faces);

    // Step 3: demux output triangles back to per-fault soups using the
    // visitor-tracked parent indices.
    std::vector<std::vector<Triangle>> per_fault_tris(N);
    std::vector<std::vector<Kernel::Point_3>> per_fault_pts(N);
    std::vector<std::unordered_map<std::size_t, std::size_t>>
        global_to_local(N);

    auto get_local_idx = [&](std::size_t fid, std::size_t g_idx)
            -> std::size_t {
        auto& m = global_to_local[fid];
        auto it = m.find(g_idx);
        if (it != m.end()) return it->second;
        const std::size_t l_idx = per_fault_pts[fid].size();
        per_fault_pts[fid].push_back(all_points[g_idx]);
        m.emplace(g_idx, l_idx);
        return l_idx;
    };

    std::int64_t n_no_parent = 0;
    std::int64_t n_bad_fault = 0;
    for (std::size_t t = 0; t < all_tris.size(); ++t) {
        const std::int64_t parent_idx =
            (t < vis.parent->size()) ? (*vis.parent)[t] : -1;
        if (parent_idx < 0
            || parent_idx >= static_cast<std::int64_t>(tri_to_fault.size())) {
            ++n_no_parent;
            continue;
        }
        const std::int64_t fid = tri_to_fault[parent_idx];
        if (fid < 0 || fid >= static_cast<std::int64_t>(N)) {
            ++n_bad_fault;
            continue;
        }
        const auto& T = all_tris[t];
        per_fault_tris[fid].push_back(Triangle{
            get_local_idx(fid, T[0]),
            get_local_idx(fid, T[1]),
            get_local_idx(fid, T[2])});
    }
    if (n_no_parent || n_bad_fault) {
        std::cerr << "[autorefine_faults] WARN demux: "
                  << n_no_parent << " no-parent, "
                  << n_bad_fault << " bad-fault triangles dropped\n";
    }

    // Step 4: write per-fault STLs (use existing Surface_mesh writer
    // by building a Mesh from the per-fault soup).
    for (std::size_t k = 0; k < N; ++k) {
        res.per_fault_n_output[k] =
            static_cast<std::int64_t>(per_fault_tris[k].size());
        if (per_fault_tris[k].empty()) {
            std::cerr << "[autorefine_faults] WARN: fault " << shorts[k]
                      << " has 0 faces post-autorefine; skipping STL "
                         "write.\n";
            continue;
        }
        Mesh per_fault_mesh;
        std::vector<Mesh::Vertex_index> v_handles;
        v_handles.reserve(per_fault_pts[k].size());
        for (const auto& p : per_fault_pts[k]) {
            v_handles.push_back(per_fault_mesh.add_vertex(p));
        }
        std::int64_t n_dropped_nonmanifold = 0;
        for (const auto& T : per_fault_tris[k]) {
            const auto fh = per_fault_mesh.add_face(
                v_handles[T[0]], v_handles[T[1]], v_handles[T[2]]);
            if (fh == Mesh::null_face()) {
                ++n_dropped_nonmanifold;
            }
        }
        if (n_dropped_nonmanifold > 0) {
            std::cerr << "[autorefine_faults] " << shorts[k]
                      << ": dropped " << n_dropped_nonmanifold
                      << " non-manifold face(s) when building per-fault "
                         "Surface_mesh (snap rounding may have collapsed "
                         "degenerate triangles).\n";
        }
        if (!safs::io::write_ascii_stl(out_paths[k], per_fault_mesh,
                                        shorts[k])) {
            throw std::runtime_error(
                "[autorefine_faults] failed to write STL "
                + out_paths[k].string());
        }
    }

    return res;
}

#else // !SAFS_HAVE_AUTOREFINE_SNAP — CGAL 5.6 legacy path

// Legacy path: experimental::autorefine on Surface_mesh, point-in-
// triangle parent assignment.  Kept for builds against the original
// CGAL 5.6.1 baseline.  Documented sliver behaviour: ~195 slivers,
// γ_min ~ 1e-4 on the clean6 2000m build (no snap rounding).

AutorefineResult autorefine_faults(
    const std::vector<std::filesystem::path>& in_paths,
    const std::vector<std::filesystem::path>& out_paths,
    const std::vector<std::string>& shorts,
    double target_edge_m,
    int    remesh_iters) {
    // target_edge_m / remesh_iters now used by Fix P pre-remesh.
    if (in_paths.size() != out_paths.size()
        || in_paths.size() != shorts.size()) {
        throw std::runtime_error(
            "[autorefine_faults] in/out/short vector size mismatch");
    }
    const std::size_t N = in_paths.size();
    if (N == 0) {
        throw std::runtime_error("[autorefine_faults] no input faults");
    }
    AutorefineResult res;
    res.per_fault_n_input.resize(N);
    res.per_fault_n_output.resize(N);
    std::vector<Kernel::Point_3> all_points;
    std::vector<std::vector<std::size_t>> all_polygons;
    std::vector<std::int64_t> face_to_fault;
    for (std::size_t k = 0; k < N; ++k) {
        std::vector<Kernel::Point_3> pts;
        std::vector<std::vector<std::size_t>> polys;
        if (!stl_to_soup(in_paths[k], pts, polys)) {
            throw std::runtime_error(
                "[autorefine_faults] failed to read soup from "
                + in_paths[k].string());
        }
        // Fix P: pre-autorefine isotropic remeshing on this single
        // fault.  Safe — no cross-fault polyline exists yet.
        fixP_pre_remesh_soup(pts, polys, shorts[k],
                              target_edge_m, remesh_iters);
        const std::size_t off = all_points.size();
        for (auto& p : pts) all_points.push_back(p);
        for (auto& poly : polys) {
            std::vector<std::size_t> shifted;
            shifted.reserve(poly.size());
            for (auto idx : poly) shifted.push_back(idx + off);
            all_polygons.push_back(std::move(shifted));
            face_to_fault.push_back(static_cast<std::int64_t>(k));
        }
        res.per_fault_n_input[k] = static_cast<std::int64_t>(polys.size());
    }
    res.n_input_faces = static_cast<std::int64_t>(all_polygons.size());
    PMP::orient_polygon_soup(all_points, all_polygons);

    Mesh combined;
    PMP::polygon_soup_to_polygon_mesh(all_points, all_polygons, combined);
    auto fault_pm = combined.add_property_map<Mesh::Face_index, std::int64_t>(
        "f:fault", -1).first;
    {
        std::int64_t fi = 0;
        for (auto f : combined.faces()) {
            if (fi < static_cast<std::int64_t>(face_to_fault.size())) {
                put(fault_pm, f, face_to_fault[fi]);
            }
            ++fi;
        }
    }
    struct PreTri {
        Kernel::Point_3 a, b, c;
        Kernel::Vector_3 normal;
        double area2;
        std::int64_t fault;
    };
    std::vector<PreTri> pre_tris;
    pre_tris.reserve(combined.number_of_faces());
    for (auto f : combined.faces()) {
        auto h0 = combined.halfedge(f);
        auto h1 = combined.next(h0);
        auto h2 = combined.next(h1);
        PreTri pt;
        pt.a = combined.point(combined.target(h0));
        pt.b = combined.point(combined.target(h1));
        pt.c = combined.point(combined.target(h2));
        pt.normal = CGAL::cross_product(pt.b - pt.a, pt.c - pt.a);
        pt.area2 = pt.normal.squared_length();
        pt.fault = get(fault_pm, f);
        pre_tris.push_back(pt);
    }
    res.had_self_intersections = PMP::does_self_intersect(combined);
    if (res.had_self_intersections) {
        try {
            PMP::experimental::autorefine(combined);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                std::string("[autorefine_faults] CGAL autorefine FAILED: ")
                + e.what());
        }
    }
    for (auto f : combined.faces()) {
        auto h0 = combined.halfedge(f);
        auto h1 = combined.next(h0);
        auto h2 = combined.next(h1);
        const auto& a = combined.point(combined.target(h0));
        const auto& b = combined.point(combined.target(h1));
        const auto& c = combined.point(combined.target(h2));
        Kernel::Point_3 cent(
            (a.x() + b.x() + c.x()) / 3.0,
            (a.y() + b.y() + c.y()) / 3.0,
            (a.z() + b.z() + c.z()) / 3.0);
        std::int64_t best = -1;
        double best_score = 1e300;
        for (std::size_t k = 0; k < pre_tris.size(); ++k) {
            const auto& pt = pre_tris[k];
            if (pt.area2 <= 0.0) continue;
            const auto va = pt.b - cent;
            const auto vb = pt.c - cent;
            const auto vc = pt.a - cent;
            const double sa = CGAL::scalar_product(pt.normal,
                CGAL::cross_product(va, vb)) / pt.area2;
            const double sb = CGAL::scalar_product(pt.normal,
                CGAL::cross_product(vb, vc)) / pt.area2;
            const double sc = 1.0 - sa - sb;
            constexpr double kBaryEps = 1e-7;
            const bool inside = (sa >= -kBaryEps)
                             && (sb >= -kBaryEps)
                             && (sc >= -kBaryEps);
            if (!inside) continue;
            const auto dvec = cent - pt.a;
            const double perp = CGAL::scalar_product(pt.normal, dvec);
            const double perp2 = perp * perp / pt.area2;
            if (perp2 < best_score) {
                best_score = perp2;
                best = pt.fault;
            }
        }
        if (best < 0) {
            // Fallback: nearest pre-centroid.
            double best_d2 = 1e300;
            for (std::size_t k = 0; k < pre_tris.size(); ++k) {
                const auto& pt = pre_tris[k];
                const double dx = cent.x() - (pt.a.x() + pt.b.x() + pt.c.x())/3.0;
                const double dy = cent.y() - (pt.a.y() + pt.b.y() + pt.c.y())/3.0;
                const double dz = cent.z() - (pt.a.z() + pt.b.z() + pt.c.z())/3.0;
                const double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 < best_d2) { best_d2 = d2; best = pt.fault; }
            }
        }
        put(fault_pm, f, best);
    }
    std::vector<Mesh> per_fault_mesh(N);
    {
        std::vector<std::vector<Mesh::Vertex_index>> v_remap(N);
        for (std::size_t k = 0; k < N; ++k) {
            v_remap[k].assign(combined.number_of_vertices(),
                              Mesh::null_vertex());
        }
        for (auto f : combined.faces()) {
            const std::int64_t fid = get(fault_pm, f);
            if (fid < 0 || fid >= static_cast<std::int64_t>(N)) continue;
            auto h0 = combined.halfedge(f);
            auto h1 = combined.next(h0);
            auto h2 = combined.next(h1);
            const auto v_combined = std::array{
                combined.target(h0), combined.target(h1), combined.target(h2)};
            std::array<Mesh::Vertex_index, 3> v_local;
            for (int t = 0; t < 3; ++t) {
                auto& slot = v_remap[fid][v_combined[t].idx()];
                if (slot == Mesh::null_vertex()) {
                    slot = per_fault_mesh[fid].add_vertex(
                        combined.point(v_combined[t]));
                }
                v_local[t] = slot;
            }
            if (v_local[0] == v_local[1]
                || v_local[1] == v_local[2]
                || v_local[0] == v_local[2]) {
                continue;
            }
            per_fault_mesh[fid].add_face(v_local[0], v_local[1], v_local[2]);
        }
    }

    // Fix B attempt (2026-04-30): tried per-fault
    // remove_almost_degenerate_faces with edge_is_constrained_map +
    // vertex_is_constrained_map covering all border edges/vertices.
    // Both with and without cap_threshold = -1 (no flips), HXT
    // rejected the PLC ("Segment and Facet intersect at point ...").
    // Root cause: collapse of an interior edge (a, b) where one
    // endpoint a is on a polyline (constrained) merges b into a,
    // RE-ROUTING every edge that was incident to b so it now
    // terminates at a.  Those re-routed edges traverse new 3D
    // paths that can intersect fault B's triangulation, breaking
    // cross-fault conformity even though no vertex moved.  The
    // `vertex_is_constrained_map` guarantee ("constrained vertex
    // stays in mesh") does not extend to the geometry of edges in
    // the constrained vertex's 1-ring.  Disabled.
    //
    // A safer Fix B would constrain *every* polyline-incident edge
    // (1-ring of every polyline endpoint), not just polyline edges
    // themselves — but that leaves only deep-interior fault edges
    // collapsible, which is most of the short-edge slivers anyway.
    // Saving for follow-up; this revert keeps the build at the
    // 195-sliver / γ_min=1.36e-4 baseline.
    (void)fixB_collapse_short_edges;

    for (std::size_t k = 0; k < N; ++k) {
        res.per_fault_n_output[k] =
            static_cast<std::int64_t>(per_fault_mesh[k].number_of_faces());
        if (per_fault_mesh[k].number_of_faces() == 0) continue;
        if (!safs::io::write_ascii_stl(out_paths[k], per_fault_mesh[k],
                                        shorts[k])) {
            throw std::runtime_error(
                "[autorefine_faults] failed to write STL "
                + out_paths[k].string());
        }
    }
    res.n_output_faces = static_cast<std::int64_t>(combined.number_of_faces());
    res.n_intersections_resolved =
        std::max<std::int64_t>(0, res.n_output_faces - res.n_input_faces);
    return res;
}

#endif // SAFS_HAVE_AUTOREFINE_SNAP

} // namespace safs::corefine
