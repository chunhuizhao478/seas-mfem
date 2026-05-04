// SAFS corefine_faults — I/O helpers (STL ASCII + diagnostic JSON).
//
// Phase 1 of PLAN_cgal_corefine.md.  See plan §Phase 1 → Detailed
// requirements §2 (STL I/O) and §3 (JSON I/O).
//
// STL writer mirrors `conformalize_faults.py:_write_ascii_stl`
// (header `solid SAFS:<short>:conformal`, `%+.17e` triple precision)
// so downstream `_combine_stls` snap-dedup sees no float32 truncation.
//
// JSON schemas are defined in PLAN_cgal_corefine.md §Constraints /
// Interface constraints.  Hand-rolled because the schemas are fixed
// and small — no external JSON library dependency.

#ifndef SAFS_TOOLS_COREFINE_FAULTS_IO_HPP
#define SAFS_TOOLS_COREFINE_FAULTS_IO_HPP

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

namespace safs::io {

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh   = CGAL::Surface_mesh<Kernel::Point_3>;

// Standard CGAL polygon-soup recipe (PLAN §Phase 1 §2):
// read_polygon_soup → repair_polygon_soup → orient_polygon_soup →
// polygon_soup_to_polygon_mesh.  ASCII STL stores three vertex
// slots per triangle; without `repair_polygon_soup` every shared
// vertex appears K times and the soup is non-manifold.
bool read_ascii_stl(const std::filesystem::path& p,
                    Mesh& out_mesh,
                    std::string* out_solid_name = nullptr);

// %+.17e triple precision; framing
//   solid SAFS:<short>:conformal
//   ...
//   endsolid SAFS:<short>:conformal
bool write_ascii_stl(const std::filesystem::path& p,
                     const Mesh& mesh,
                     std::string_view solid_name);

struct PerFaultStats {
    std::string short_name;
    std::int64_t n_vertices   = 0;
    std::int64_t n_triangles  = 0;
    std::filesystem::path stl_path;
};

struct PairStats {
    std::string short_a, short_b;
    std::int64_t n_constrained_edges_A = 0;
    std::int64_t n_constrained_edges_B = 0;
    std::int64_t pre_n_tri_A  = 0;
    std::int64_t post_n_tri_A = 0;
    std::int64_t pre_n_tri_B  = 0;
    std::int64_t post_n_tri_B = 0;
    int  remesh_iters         = 0;   // 0 in Phase 1
    int  n_t_junctions_A      = 0;   // stub in Phase 1
    int  n_t_junctions_B      = 0;
    std::string gate_manifold_A           = "PASS";
    std::string gate_manifold_B           = "PASS";
    std::string gate_polyline_coincidence = "PASS";
    std::string gate_interior_only        = "PASS";
};

// triangle_to_fault.json — presence-only marker downstream.
// Schema-invariant: ranges contiguous + partition [0, n_total).
void write_triangle_to_fault_json(
    const std::filesystem::path& p,
    const std::vector<PerFaultStats>& per_fault);

// intersection_report.json — diagnostic only.  See plan §Constraints
// / Schema delta vs Python output (P-007).
void write_intersection_report_json(
    const std::filesystem::path& p,
    const std::vector<PerFaultStats>& per_fault,
    const std::vector<PairStats>& pairs,
    double snap_m, double clearance_m, double target_edge_length_m);

// Validate the just-written triangle_to_fault.json against the
// in-memory `per_fault` vector (R-001 + R-005).  Re-reads the file
// so a writer corruption is caught.  Returns 0 on success, 3 on
// any schema violation (invariant the implementer wants to expose
// as exit code 3 in main).  Logs failures to stderr.
//
// Checks performed:
//   (a) `n_total_triangles` key present and equal to
//       sum(per_fault[*].n_triangles);
//   (b) Per-fault range arrays form a contiguous partition of
//       [0, n_total_triangles) in the order `per_fault` was given.
int validate_triangle_to_fault_json(
    const std::filesystem::path& p,
    const std::vector<PerFaultStats>& per_fault);

} // namespace safs::io

#endif // SAFS_TOOLS_COREFINE_FAULTS_IO_HPP
