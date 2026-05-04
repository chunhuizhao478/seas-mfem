// SAFS corefine_faults — autorefine-mode driver.
//
// Replaces the pairwise-corefine cascade (which accumulates near-tangent
// subdivision artifacts on faults that participate in 3+ pairs and
// produces post-corefine self-intersection — see REVIEW notes in the
// 6-fault debug session, 2026-04-30) with a single-pass approach:
// concatenate all faults' triangles into one polygon soup, run
// `PMP::autorefine` (CGAL 5.6+), then split back per-fault by face-tag
// range.
//
// This trades the per-pair `polyline_ecm` granularity for ONE-shot
// global refinement.  Cross-fault conformity is preserved by
// construction because autorefine produces shared edges along every
// pair-of-fault intersection in a single combined operation, with no
// pairwise cascade to accumulate artifacts.

#ifndef SAFS_TOOLS_COREFINE_FAULTS_AUTOREFINE_MODE_HPP
#define SAFS_TOOLS_COREFINE_FAULTS_AUTOREFINE_MODE_HPP

#include <filesystem>
#include <string>
#include <vector>

#include "corefine.hpp"

namespace safs::corefine {

struct AutorefineResult {
    std::int64_t n_input_faces  = 0;
    std::int64_t n_output_faces = 0;
    std::int64_t n_intersections_resolved = 0;
    bool         had_self_intersections   = false;
    std::vector<std::int64_t> per_fault_n_input;
    std::vector<std::int64_t> per_fault_n_output;
};

// Run autorefine-mode on N input STLs.
//   in_paths   — N STL files (per-fault input).
//   out_paths  — N STL files (per-fault output, in 1:1 correspondence).
//   shorts     — fault short names (for the JSON sidecar).
// On entry the output directory must already exist.  The function
// writes per-fault output STLs and returns the result statistics.
//
// Throws std::runtime_error on read/write failure or on autorefine
// algorithmic failure.
AutorefineResult autorefine_faults(
    const std::vector<std::filesystem::path>& in_paths,
    const std::vector<std::filesystem::path>& out_paths,
    const std::vector<std::string>& shorts,
    double target_edge_m = 0.0,    // 0 = no per-fault remesh; >0 =
                                    // post-autorefine cleanup +
                                    // isotropic_remesh to target.
    int    remesh_iters  = 3);

} // namespace safs::corefine

#endif // SAFS_TOOLS_COREFINE_FAULTS_AUTOREFINE_MODE_HPP
