// SAFS corefine_faults — driver.
//
// Phase 1 of PLAN_cgal_corefine.md.  See plan §Phase 1 → Detailed
// requirements §1 (CLI), §5 (pairwise orchestration), §6 (gates).
//
// Usage:
//   corefine_faults
//     --in-stl-dir   <p>
//     --out-stl-dir  <p>
//     --include-fault <short>  (1+, order matters)
//     [--snap-m       <f>]   default 1e-3
//     [--clearance-m  <f>]   default 0.0
//     [--target-edge-m <f>]  Phase 2 — Phase 1 errors if non-zero
//     [--remesh-iters  <i>]  Phase 2 — Phase 1 errors if non-zero
//     [--verbose]
//
// Exit codes (per plan §Constraints / Convention constraints):
//   0  success
//   1  user error (bad flags / missing files)
//   2  CGAL runtime error (corefine failed / non-manifold)
//   3  schema-validation error
//
// Numerical kernel: Exact_predicates_inexact_constructions_kernel
// (EPICK).  Per CGAL docs and plan §Numerical constraints, EPICK
// is sufficient for a single corefinement; EPECK is only needed
// for *consecutive* boolean operations.

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include "corefine.hpp"
#include "io.hpp"
#include "autorefine_mode.hpp"

namespace fs = std::filesystem;
namespace PMP = CGAL::Polygon_mesh_processing;

using Mesh   = safs::corefine::Mesh;
using ECMap  = safs::corefine::EdgeConstrainedMap;

namespace {

struct Args {
    fs::path in_stl_dir;
    fs::path out_stl_dir;
    std::vector<std::string> included;
    double snap_m         = 1e-3;
    double clearance_m    = 0.0;
    // Phase 2 defaults — PLAN_cgal_corefine.md §Phase 2 §6.
    double target_edge_m  = 1000.0;
    int    remesh_iters   = 3;
    bool   no_remesh      = false;
    bool   allow_z_clamp  = false;
    bool   allow_polyline_asymmetry = false;   // R-cascade: convert
                                                // post-corefine
                                                // delta-gate FAIL into
                                                // a warning so cascade
                                                // pairs that hit CGAL's
                                                // multi-pair near-tangency
                                                // asymmetry can proceed.
    bool   verbose        = false;
    // Mode selector: "cascade" (default — pairwise corefine, the
    // pre-existing path) or "autorefine" (single-pass; required for
    // 3+ faults to avoid cumulative subdivision artifacts that break
    // the cascade gate).
    std::string mode = "cascade";
};

void usage() {
    std::cerr <<
        "Usage: corefine_faults --in-stl-dir <p> --out-stl-dir <p> \\\n"
        "                       --include-fault <s> [--include-fault <s>]... \\\n"
        "                       [--snap-m <f>]         (default 1e-3)\n"
        "                       [--clearance-m <f>]    (default 0)\n"
        "                       [--target-edge-m <f>]  (default 1000; 0 disables)\n"
        "                       [--remesh-iters <i>]   (default 3)\n"
        "                       [--no-remesh]          (Phase-1-equivalent: skip remeshing)\n"
        "                       [--allow-z-clamp]      (P-008 opt-in: project violators to -clearance_m)\n"
        "                       [--verbose]\n";
}

bool parse_double(const char* s, double& out) {
    try {
        std::size_t pos = 0;
        out = std::stod(s, &pos);
        return pos == std::string_view(s).size();
    } catch (...) { return false; }
}

bool parse_int(const char* s, int& out) {
    try {
        std::size_t pos = 0;
        out = std::stoi(s, &pos);
        return pos == std::string_view(s).size();
    } catch (...) { return false; }
}

int parse_args(int argc, char** argv, Args& out) {
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        auto take_val = [&](const char* flag) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "[corefine_faults] " << flag
                          << " requires a value.\n";
                return nullptr;
            }
            return argv[++i];
        };

        if (a == "--in-stl-dir") {
            const char* v = take_val("--in-stl-dir");
            if (!v) return 1;
            out.in_stl_dir = v;
        } else if (a == "--out-stl-dir") {
            const char* v = take_val("--out-stl-dir");
            if (!v) return 1;
            out.out_stl_dir = v;
        } else if (a == "--include-fault") {
            const char* v = take_val("--include-fault");
            if (!v) return 1;
            out.included.emplace_back(v);
        } else if (a == "--snap-m") {
            const char* v = take_val("--snap-m");
            if (!v || !parse_double(v, out.snap_m)) {
                std::cerr << "[corefine_faults] bad --snap-m\n";
                return 1;
            }
        } else if (a == "--clearance-m") {
            const char* v = take_val("--clearance-m");
            if (!v || !parse_double(v, out.clearance_m)) {
                std::cerr << "[corefine_faults] bad --clearance-m\n";
                return 1;
            }
        } else if (a == "--target-edge-m") {
            const char* v = take_val("--target-edge-m");
            if (!v || !parse_double(v, out.target_edge_m)) {
                std::cerr << "[corefine_faults] bad --target-edge-m\n";
                return 1;
            }
        } else if (a == "--remesh-iters") {
            const char* v = take_val("--remesh-iters");
            if (!v || !parse_int(v, out.remesh_iters)) {
                std::cerr << "[corefine_faults] bad --remesh-iters\n";
                return 1;
            }
        } else if (a == "--no-remesh") {
            out.no_remesh = true;
        } else if (a == "--allow-z-clamp") {
            out.allow_z_clamp = true;
        } else if (a == "--allow-polyline-asymmetry") {
            out.allow_polyline_asymmetry = true;
        } else if (a == "--mode") {
            const char* v = take_val("--mode");
            if (!v) return 1;
            std::string m = v;
            if (m != "cascade" && m != "autorefine") {
                std::cerr << "[corefine_faults] --mode must be 'cascade' "
                             "or 'autorefine'; got '" << m << "'.\n";
                return 1;
            }
            out.mode = m;
        } else if (a == "--verbose") {
            out.verbose = true;
        } else if (a == "-h" || a == "--help") {
            usage();
            return -1;          // sentinel: print and exit 0
        } else {
            std::cerr << "[corefine_faults] unknown flag: " << a << '\n';
            usage();
            return 1;
        }
    }

    if (out.in_stl_dir.empty()) {
        std::cerr << "[corefine_faults] --in-stl-dir is required.\n";
        return 1;
    }
    if (out.out_stl_dir.empty()) {
        std::cerr << "[corefine_faults] --out-stl-dir is required.\n";
        return 1;
    }
    if (out.included.empty()) {
        std::cerr << "[corefine_faults] must specify at least one "
                     "--include-fault.\n";
        return 1;
    }
    if (out.target_edge_m < 0.0) {
        std::cerr << "[corefine_faults] --target-edge-m must be >= 0\n";
        return 1;
    }
    if (out.remesh_iters < 0) {
        std::cerr << "[corefine_faults] --remesh-iters must be >= 0\n";
        return 1;
    }
    if (out.no_remesh) {
        // --no-remesh is a convenience alias that overrides
        // --target-edge-m to 0 (back-compat to Phase 1 behaviour).
        out.target_edge_m = 0.0;
        out.remesh_iters  = 0;
    }
    return 0;
}

// PLAN §Phase 1 §6 — manifold gate.  CGAL exposes an explicit
// non-manifold check; "polyhedral" here means closed + manifold,
// which is too strong for open fault patches.  Use the disjunction
// of `does_self_intersect` (must be false) and `is_valid_polygon_
// mesh` (must be true).
std::string manifold_gate(const Mesh& m) {
    if (!CGAL::is_valid_polygon_mesh(m)) {
        return "FAIL_INVALID_POLYGON_MESH";
    }
    if (PMP::does_self_intersect(m)) {
        return "FAIL_SELF_INTERSECT";
    }
    return "PASS";
}

// PLAN §Phase 1 §6 — interior-crossing-only gate.
// Re-corefine each pair IN PLACE with FRESH empty constrained maps;
// assert that no new triangles are introduced (i.e., the re-corefine
// doesn't change the triangulation).
//
// R-106: drops the deep-copy of the meshes (~36 MB peak on the
// all-8 build) and instead installs a temporary ECM on the live
// meshes.  A correctly-conformal pair is a no-op for corefine —
// face count is unchanged — so the live meshes are not mutated.
// If a previous bug causes corefine to mutate, we surface that as
// a gate failure AND the meshes are at most as fine-grained as
// they would have been; the next pair's corefine sees the same
// intersection topology either way.
//
// Deviation from the plan's claimed invariant: the plan said the
// re-corefine should produce 0 constrained edges.  Empirically
// (CGAL 5.6.1), `PMP::corefine` re-marks the already-existing
// polyline edges in the ECM — the ECM marks every intersection
// edge regardless of whether the edge was newly introduced.  So
// the *geometric* invariant we actually want is that the
// triangulation is unchanged.
bool interior_crossing_only(Mesh& A, Mesh& B,
                            std::string_view short_a,
                            std::string_view short_b) {
    auto pa = A.add_property_map<safs::corefine::EdgeDescriptor, bool>(
                  "e:gate_idem", false);
    auto pb = B.add_property_map<safs::corefine::EdgeDescriptor, bool>(
                  "e:gate_idem", false);
    if (!pa.second || !pb.second) {
        std::cerr << "[corefine_faults] gate: failed to install fresh "
                     "ECM on live meshes for pair ("
                  << short_a << ", " << short_b << ").\n";
        return false;
    }
    const auto a_pre = A.number_of_faces();
    const auto b_pre = B.number_of_faces();
    bool ok = true;
    try {
        safs::corefine::corefine_pair(A, B, pa.first, pb.first,
                                      short_a, short_b);
    } catch (const std::exception& e) {
        std::cerr << "[corefine_faults] gate: " << e.what() << '\n';
        ok = false;
    }
    A.remove_property_map(pa.first);
    B.remove_property_map(pb.first);
    if (!ok) return false;
    return A.number_of_faces() == a_pre
        && B.number_of_faces() == b_pre;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    int prc = parse_args(argc, argv, args);
    if (prc == -1) return EXIT_SUCCESS;
    if (prc != 0) return prc;

    if (args.verbose) {
        std::cerr << "[corefine_faults] in_stl_dir=" << args.in_stl_dir
                  << " out_stl_dir=" << args.out_stl_dir
                  << " snap_m=" << args.snap_m
                  << " clearance_m=" << args.clearance_m << '\n';
        for (const auto& s : args.included) {
            std::cerr << "[corefine_faults]   --include-fault " << s << '\n';
        }
    }

    const std::size_t N = args.included.size();
    std::vector<Mesh> meshes(N);
    std::vector<safs::io::PerFaultStats> per_fault(N);

    // 1. Read every fault.
    for (std::size_t i = 0; i < N; ++i) {
        const fs::path in = args.in_stl_dir / (args.included[i] + ".stl");
        if (!fs::exists(in)) {
            std::cerr << "[corefine_faults] missing input STL: " << in << '\n';
            return 1;
        }
        if (!safs::io::read_ascii_stl(in, meshes[i])) {
            std::cerr << "[corefine_faults] failed to read " << in << '\n';
            return 1;
        }
        if (args.verbose) {
            std::cerr << "[corefine_faults] read " << in
                      << " V=" << meshes[i].number_of_vertices()
                      << " F=" << meshes[i].number_of_faces() << '\n';
        }
        per_fault[i].short_name = args.included[i];
    }

    // 2. Install TWO ECMs per fault BEFORE any corefine runs (P-002 +
    //    R-102 + R-105):
    //      polyline_ecms[i] — read-only set of intersection-polyline
    //                          edges; mutated only by corefine_pair
    //                          and (later) collapse_short_polyline_
    //                          edges_symmetric.  Used by the post-
    //                          remesh polyline-coincidence gate.
    //      protect_ecms[i]  — mutable union (polyline ∪ boundary)
    //                          used as the protect-and-don't-touch
    //                          set during split_long_edges +
    //                          isotropic_remeshing.
    std::vector<ECMap> polyline_ecms(N);
    std::vector<ECMap> protect_ecms(N);
    for (std::size_t i = 0; i < N; ++i) {
        auto pl = meshes[i]
            .add_property_map<safs::corefine::EdgeDescriptor, bool>(
                "e:polyline", false);
        auto pr = meshes[i]
            .add_property_map<safs::corefine::EdgeDescriptor, bool>(
                "e:protect", false);
        if (!pl.second || !pr.second) {
            // R-207: this is a programmer error (we own the property
            // tag string), not a CGAL runtime error.  Exit 1 (user /
            // programmer error) rather than 2 (CGAL runtime).
            std::cerr << "[corefine_faults] internal: ECM property "
                         "map already exists on fault "
                      << args.included[i] << " — duplicate tag "
                         "elsewhere in the codebase.\n";
            return 1;
        }
        polyline_ecms[i] = pl.first;
        protect_ecms[i]  = pr.first;
    }

    // ============================================================
    // Autorefine mode short-circuit
    // ============================================================
    // For 3+ faults the pairwise cascade accumulates near-tangent
    // subdivision artifacts that produce post-corefine self-intersection
    // (debug session 2026-04-30).  In autorefine mode, concatenate all
    // faults into one polygon soup, run CGAL's
    // `autorefine_and_remove_self_intersections` once, then split back
    // by fault tag — no cascade, no accumulation.
    if (args.mode == "autorefine") {
        std::cerr << "[corefine_faults] mode=autorefine: single-pass "
                     "global refinement (replaces pairwise cascade).\n";
        std::vector<fs::path> in_paths(N), out_paths(N);
        std::vector<std::string> shorts = args.included;
        std::error_code ec;
        fs::create_directories(args.out_stl_dir, ec);
        for (std::size_t k = 0; k < N; ++k) {
            in_paths[k]  = args.in_stl_dir  / (shorts[k] + ".stl");
            out_paths[k] = args.out_stl_dir / (shorts[k] + ".stl");
            if (!fs::exists(in_paths[k])) {
                std::cerr << "[corefine_faults] missing input STL: "
                          << in_paths[k] << '\n';
                return 1;
            }
        }
        try {
            auto rr = safs::corefine::autorefine_faults(
                in_paths, out_paths, shorts,
                args.no_remesh ? 0.0 : args.target_edge_m,
                args.no_remesh ? 0   : args.remesh_iters);
            std::cerr << "[corefine_faults] autorefine: "
                      << rr.n_input_faces << " input faces → "
                      << rr.n_output_faces << " output faces ("
                      << rr.n_intersections_resolved << " new from "
                         "intersection refinement); had self-int: "
                      << (rr.had_self_intersections ? "yes" : "no")
                      << ".\n";
            for (std::size_t k = 0; k < N; ++k) {
                std::cerr << "[corefine_faults]   " << shorts[k] << ": "
                          << rr.per_fault_n_input[k] << " → "
                          << rr.per_fault_n_output[k] << " faces\n";
                per_fault[k].n_vertices = 0;  // populated post-write below
                per_fault[k].n_triangles = rr.per_fault_n_output[k];
                per_fault[k].stl_path    = out_paths[k];
            }
        } catch (const std::exception& e) {
            std::cerr << "[corefine_faults] autorefine FAILED: "
                      << e.what() << '\n';
            return 2;
        }

        // Schema sidecars (no per-pair stats in autorefine mode).
        const fs::path tri2fault = args.out_stl_dir / "triangle_to_fault.json";
        const fs::path report    = args.out_stl_dir / "intersection_report.json";
        std::vector<safs::io::PairStats> empty_pairs;
        try {
            // Re-read each STL to populate accurate vertex counts, since
            // autorefine_faults wrote them but we didn't track vertex
            // counts in AutorefineResult.
            for (std::size_t k = 0; k < N; ++k) {
                Mesh m;
                if (safs::io::read_ascii_stl(out_paths[k], m)) {
                    per_fault[k].n_vertices = static_cast<std::int64_t>(
                        m.number_of_vertices());
                }
            }
            safs::io::write_triangle_to_fault_json(tri2fault, per_fault);
            safs::io::write_intersection_report_json(
                report, per_fault, empty_pairs,
                args.snap_m, args.clearance_m, args.target_edge_m);
        } catch (const std::exception& e) {
            std::cerr << "[corefine_faults] schema write failed: "
                      << e.what() << '\n';
            return 3;
        }
        if (int rc = safs::io::validate_triangle_to_fault_json(
                tri2fault, per_fault); rc != 0) {
            return rc;
        }
        std::cerr << "[corefine_faults] autorefine done.\n";
        return EXIT_SUCCESS;
    }

    // 2b. R-303 — pre-corefine vertex weld.  Snap pairs of vertices
    //     (a ∈ fault i, b ∈ fault j) within `weld_tol` of each other
    //     to a common midpoint position, eliminating the near-vertex
    //     pierce that produces sub-target sliver edges in corefine
    //     output.  Cross-fault conformity is preserved by construction
    //     because the midpoint is the same on both sides.
    //
    //     Skipped when `--target-edge-m == 0` (no remeshing planned →
    //     no need to keep input vertex spacing close to target).
    if (args.target_edge_m > 0.0) {
        const double weld_tol = args.target_edge_m / 4.0;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i + 1; j < N; ++j) {
                const int n_welded = safs::corefine::preweld_near_vertices(
                    meshes[i], meshes[j], weld_tol);
                if (args.verbose && n_welded > 0) {
                    std::cerr << "[corefine_faults] pre-weld pair "
                              << args.included[i] << " x "
                              << args.included[j] << ": " << n_welded
                              << " vertex pair(s) welded at tol "
                              << weld_tol << " m (R-303)\n";
                }
            }
        }
    }

    // 3. Coincident-triangle preflight (P-011) on every unordered pair.
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            auto [ai, bj] = safs::corefine::find_coincident_triangle(
                meshes[i], meshes[j]);
            if (ai >= 0) {
                std::cerr << "[corefine_faults] fault "
                          << args.included[i] << " tri " << ai
                          << " and fault " << args.included[j]
                          << " tri " << bj
                          << " share all 3 vertices.  CGAL corefine() "
                             "has undefined behaviour on coincident "
                             "triangles.  Investigate the CFM source — "
                             "this typically indicates a triangulation "
                             "that should not have been split into two "
                             "faults (PLAN_multifault_intersections.md "
                             "§Caveats).\n";
                return 2;
            }
        }
    }

    // 4. Pairwise corefine in ascending (i, j) order.
    std::vector<safs::io::PairStats> pair_stats;
    pair_stats.reserve(N * (N - 1) / 2);

    auto t0 = std::chrono::steady_clock::now();

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            safs::io::PairStats ps;
            ps.short_a = args.included[i];
            ps.short_b = args.included[j];

            // Cascade gate fix: snapshot the polyline-endpoint sets
            // BEFORE corefine_pair so we can compute the DELTA (edges
            // added by this pair only).  Without the snapshot, the
            // gate compares accumulated marks from all prior pairs
            // and falsely fails on multi-pair runs (a fault's ECM
            // has marks from earlier pairs that the current partner
            // cannot share).
            const auto pre_a = safs::corefine::diagnose_constrained_endpoints(
                meshes[i], polyline_ecms[i]);
            const auto pre_b = safs::corefine::diagnose_constrained_endpoints(
                meshes[j], polyline_ecms[j]);

            try {
                auto cr = safs::corefine::corefine_pair(
                    meshes[i], meshes[j],
                    polyline_ecms[i], polyline_ecms[j],
                    ps.short_a, ps.short_b);
                ps.n_constrained_edges_A = cr.n_constrained_edges_A;
                ps.n_constrained_edges_B = cr.n_constrained_edges_B;
                ps.pre_n_tri_A  = cr.pre_n_tri_A;
                ps.post_n_tri_A = cr.post_n_tri_A;
                ps.pre_n_tri_B  = cr.pre_n_tri_B;
                ps.post_n_tri_B = cr.post_n_tri_B;
            } catch (const std::exception& e) {
                std::cerr << "[corefine_faults] pair (" << ps.short_a
                          << ", " << ps.short_b << "): " << e.what() << '\n';
                return 2;
            }

            // Manifold gate (post-corefine).  R-003: any FAIL is fatal.
            ps.gate_manifold_A = manifold_gate(meshes[i]);
            ps.gate_manifold_B = manifold_gate(meshes[j]);
            if (ps.gate_manifold_A != "PASS") {
                std::cerr << "[corefine_faults] manifold gate FAILED on "
                             "fault " << ps.short_a << " after corefine "
                             "with " << ps.short_b << ": "
                          << ps.gate_manifold_A << '\n';
                return 2;
            }
            if (ps.gate_manifold_B != "PASS") {
                std::cerr << "[corefine_faults] manifold gate FAILED on "
                             "fault " << ps.short_b << " after corefine "
                             "with " << ps.short_a << ": "
                          << ps.gate_manifold_B << '\n';
                return 2;
            }

            // Polyline-edge-coincidence gate (R-004): the polyline
            // edges ADDED BY THIS PAIR must be the same set of 3-D
            // endpoint pairs on both A and B.  We compute the delta
            // (post − pre) on each side and compare those, NOT the
            // full accumulated ECM (which contains marks from
            // earlier pairs that the current partner cannot have).
            const auto post_a = safs::corefine::diagnose_constrained_endpoints(
                meshes[i], polyline_ecms[i]);
            const auto post_b = safs::corefine::diagnose_constrained_endpoints(
                meshes[j], polyline_ecms[j]);
            std::set<safs::corefine::EndpointPairDiag> a_new, b_new;
            std::set_difference(post_a.begin(), post_a.end(),
                                pre_a.begin(),  pre_a.end(),
                                std::inserter(a_new, a_new.end()));
            std::set_difference(post_b.begin(), post_b.end(),
                                pre_b.begin(),  pre_b.end(),
                                std::inserter(b_new, b_new.end()));
            ps.gate_polyline_coincidence = (a_new == b_new) ? "PASS" : "FAIL";
            if (ps.gate_polyline_coincidence != "PASS") {
                std::cerr << "[corefine_faults] polyline-edge-"
                             "coincidence gate FAILED on pair ("
                          << ps.short_a << ", " << ps.short_b
                          << "): the polyline edges ADDED BY THIS "
                             "PAIR (delta over pre-corefine ECMs) "
                             "differ between fault A ("
                          << a_new.size() << " new edges) and fault B ("
                          << b_new.size() << " new edges).  This "
                             "happens in multi-pair cascades when CGAL's "
                             "corefine produces slightly asymmetric "
                             "output on a fault that has already been "
                             "corefined with prior partners (different "
                             "vertex layouts on A vs B at near-tangent "
                             "intersections).\n";
                if (!args.allow_polyline_asymmetry) {
                    std::cerr << "[corefine_faults] To proceed despite "
                                 "the asymmetry, pass "
                                 "--allow-polyline-asymmetry; the bulk "
                                 "mesh's HXT will mesh through small "
                                 "asymmetries (a few edges) but the "
                                 "post-mesh validator's check_5 may "
                                 "report a few non-conformal triangles.\n";
                    return 2;
                }
                std::cerr << "[corefine_faults] --allow-polyline-"
                             "asymmetry: continuing despite "
                          << std::abs((long)a_new.size() - (long)b_new.size())
                          << " asymmetric polyline edges.\n";
            }

            // Interior-crossing-only gate (idempotency).  R-003: FAIL is fatal.
            ps.gate_interior_only = interior_crossing_only(
                meshes[i], meshes[j], ps.short_a, ps.short_b)
                ? "PASS" : "FAIL";
            if (ps.gate_interior_only != "PASS") {
                std::cerr << "[corefine_faults] interior-crossing-only "
                             "gate FAILED on pair (" << ps.short_a
                          << ", " << ps.short_b
                          << "): re-corefine introduced new triangles, "
                             "indicating the first corefine did not "
                             "produce a fully-conformal pair.\n";
                return 2;
            }

            if (args.verbose) {
                std::cerr << "[corefine_faults] pair "
                          << ps.short_a << " x " << ps.short_b
                          << ": ce_A=" << ps.n_constrained_edges_A
                          << " ce_B=" << ps.n_constrained_edges_B
                          << " (pre A,B=" << ps.pre_n_tri_A << ','
                          << ps.pre_n_tri_B
                          << " → post=" << ps.post_n_tri_A << ','
                          << ps.post_n_tri_B << ")\n";
            }

            pair_stats.push_back(std::move(ps));
        }
    }

    // 4b. Phase 2 — isotropic remeshing per fault.  Skipped when
    //     `--no-remesh` (or `--target-edge-m 0`) is set.
    const bool do_remesh = (args.target_edge_m > 0.0);
    for (auto& ps : pair_stats) {
        ps.remesh_iters = do_remesh ? args.remesh_iters : 0;
    }
    if (do_remesh) {
        // R-301: iterative uniform polyline resampling.  Drives every
        // polyline edge into [target × 0.7, target × 1.3] by
        // alternating per-fault `PMP::split_long_edges` (deterministic
        // across faults — bit-identical endpoints produce bit-identical
        // midpoint splits) with the symmetric collapse helper.  Replaces
        // the previous "single collapse pass at target/4" which left
        // 60% of polyline edges in a 250–667 m dead zone, producing the
        // user-visible polyline-corridor density grade (R-301 in
        // REVIEW_phase02_code_v3.md).
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i + 1; j < N; ++j) {
                const auto ur = safs::corefine::
                    uniformize_polyline_edges_symmetric(
                        meshes[i], polyline_ecms[i],
                        meshes[j], polyline_ecms[j],
                        args.target_edge_m);
                if (args.verbose) {
                    std::cerr << "[corefine_faults] uniformize pair "
                              << args.included[i] << " x "
                              << args.included[j] << ": "
                              << ur.n_splits << " split(s) + "
                              << ur.n_collapses << " collapse(s) over "
                              << ur.n_iter << " iter (converged="
                              << (ur.converged ? "yes" : "no") << ")\n";
                }
                if (!ur.converged && (ur.n_splits + ur.n_collapses > 0)) {
                    std::cerr << "[corefine_faults] WARN: uniformize did "
                                 "not fully converge on pair ("
                              << args.included[i] << ", "
                              << args.included[j]
                              << ") — some polyline edges remain outside "
                                 "the [0.7, 1.3] × target band.  The most "
                                 "common cause is a polyline kink shorter "
                                 "than 0.7 × target whose collapse fails "
                                 "the link condition.  Mesh quality is "
                                 "degraded but cross-fault conformity is "
                                 "preserved.\n";
                }
            }
        }

        // Seed protect_ecms[i] with polyline_ecms[i] marks (R-105
        // union invariant; remesh_one_fault asserts this).
        for (std::size_t i = 0; i < N; ++i) {
            for (auto e : meshes[i].edges()) {
                if (get(polyline_ecms[i], e)) {
                    put(protect_ecms[i], e, true);
                }
            }
        }

        for (std::size_t i = 0; i < N; ++i) {
            safs::corefine::RemeshParams rp;
            rp.target_edge_m  = args.target_edge_m;
            rp.n_iterations   = args.remesh_iters;
            rp.clearance_m    = args.clearance_m;
            rp.allow_z_clamp  = args.allow_z_clamp;

            try {
                auto rr = safs::corefine::remesh_one_fault(
                    meshes[i], polyline_ecms[i], protect_ecms[i], rp);
                if (args.verbose) {
                    std::cerr << "[corefine_faults] remesh "
                              << args.included[i]
                              << ": F " << rr.n_tri_pre << " → "
                              << rr.n_tri_post
                              << ", protect_ce " << rr.n_constrained_pre_split
                              << " → " << rr.n_constrained_post_split
                              << " → " << rr.n_constrained_post_remesh
                              << '\n';
                }
            } catch (const std::exception& e) {
                std::cerr << "[corefine_faults] remesh on fault "
                          << args.included[i] << ": " << e.what()
                          << '\n';
                return 2;
            }

            // P-008 free-surface clearance enforcement.  Use
            // protect_ecms (polyline ∪ boundary) for "must not move
            // this vertex" semantics.
            auto cr = safs::corefine::check_free_surface_clearance(
                meshes[i], protect_ecms[i], args.clearance_m);
            if (cr.any_violation) {
                if (!args.allow_z_clamp) {
                    const auto& q = meshes[i].point(cr.worst_v);
                    std::cerr << "[corefine_faults] FATAL: vertex ("
                              << q.x() << ", " << q.y() << ", " << q.z()
                              << ") on fault " << args.included[i]
                              << " is above z = -" << args.clearance_m
                              << " m by " << cr.max_dz
                              << " m after remeshing.  This indicates "
                                 "either (a) the input STL was not "
                                 "clamped by ts_to_stl.py to the same "
                                 "clearance, or (b) a CGAL remesher "
                                 "interior-smoothing pulled a "
                                 "non-locked vertex above its 1-ring's "
                                 "neighbours.  Re-run ts_to_stl.py "
                                 "with the matching "
                                 "--free-surface-clearance value, or "
                                 "pass --allow-z-clamp to project "
                                 "violators down to the clearance "
                                 "plane (see Phase 2 §P-008 for the "
                                 "trade-offs).\n";
                    return 2;
                }
                // --allow-z-clamp set: try the safe projection.
                if (cr.worst_v_constrained) {
                    const auto& q = meshes[i].point(cr.worst_v);
                    std::cerr << "[corefine_faults] --allow-z-clamp "
                                 "refused: a CONSTRAINED vertex ("
                              << q.x() << ", " << q.y() << ", "
                              << q.z() << ") on fault "
                              << args.included[i]
                              << " is in violation; clamping it would "
                                 "break bit-equality of polyline "
                                 "endpoints between fault A and "
                                 "fault B.  Re-run ts_to_stl.py with "
                                 "the matching --free-surface-"
                                 "clearance value.\n";
                    return 1;
                }
                if (!safs::corefine::project_to_free_surface_safe(
                        meshes[i], protect_ecms[i], args.clearance_m)) {
                    std::cerr << "[corefine_faults] --allow-z-clamp "
                                 "refused: projection would invert "
                                 "at least one triangle on fault "
                              << args.included[i] << ", or a "
                                 "constrained vertex is in violation. "
                                 " Re-run ts_to_stl.py with the "
                                 "matching --free-surface-clearance "
                                 "value.\n";
                    return 1;
                }
                if (args.verbose) {
                    std::cerr << "[corefine_faults] --allow-z-clamp: "
                                 "projected interior violators on "
                                 "fault " << args.included[i]
                              << " to z = -" << args.clearance_m
                              << " m.\n";
                }
            }
        }

        // R-102: re-run polyline-edge-coincidence on EVERY pair AFTER
        // remesh.  The gate uses protect_ecms (polyline ∪ boundary)
        // and internally filters out border edges, so the comparison
        // reduces to polyline endpoints only.  We use protect_ecms
        // rather than polyline_ecms here because Property_map storage
        // is indexed by edge_index and isotropic_remeshing recycles
        // indices, leaving polyline_ecms with stale marks; protect_ecms
        // is propagated correctly by split_long_edges + isotropic_
        // remeshing's edge_is_constrained_map plumbing.  Any mismatch
        // is fatal.
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i + 1; j < N; ++j) {
                if (!safs::corefine::polyline_edge_coincidence_gate(
                        meshes[i], protect_ecms[i],
                        meshes[j], protect_ecms[j])) {
                    std::cerr << "[corefine_faults] post-remesh "
                                 "polyline-edge-coincidence FAILED "
                                 "on pair (" << args.included[i]
                              << ", " << args.included[j]
                              << "): the polyline edges of the two "
                                 "remeshed faults disagree on at "
                                 "least one 3-D endpoint.  This is a "
                                 "regression vs the post-corefine "
                                 "gate; investigate split_long_edges "
                                 "or sliver-collapse asymmetry.\n";
                    return 2;
                }
            }
        }
    }

    // 5. Write output STLs + populate per-fault stats with final counts.
    std::error_code ec;
    fs::create_directories(args.out_stl_dir, ec);

    for (std::size_t i = 0; i < N; ++i) {
        const fs::path out = args.out_stl_dir / (args.included[i] + ".stl");
        per_fault[i].n_vertices  = static_cast<std::int64_t>(
                                       meshes[i].number_of_vertices());
        per_fault[i].n_triangles = static_cast<std::int64_t>(
                                       meshes[i].number_of_faces());
        per_fault[i].stl_path    = out;

        if (!safs::io::write_ascii_stl(out, meshes[i], args.included[i])) {
            std::cerr << "[corefine_faults] failed to write " << out << '\n';
            return 1;
        }
        if (args.verbose) {
            std::cerr << "[corefine_faults] wrote " << out
                      << " V=" << per_fault[i].n_vertices
                      << " F=" << per_fault[i].n_triangles << '\n';
        }
    }

    // 6. Write JSON sidecars.
    const fs::path tri2fault = args.out_stl_dir / "triangle_to_fault.json";
    const fs::path report    = args.out_stl_dir / "intersection_report.json";

    try {
        safs::io::write_triangle_to_fault_json(tri2fault, per_fault);
        safs::io::write_intersection_report_json(
            report, per_fault, pair_stats,
            args.snap_m, args.clearance_m, args.target_edge_m);
    } catch (const std::exception& e) {
        std::cerr << "[corefine_faults] schema write failed: "
                  << e.what() << '\n';
        return 3;
    }

    // 7. Schema-invariant cross-check (R-005): the in-memory
    //    `per_fault` triangle sum must equal the live
    //    `meshes[i].number_of_faces()` sum.  Catches a stale-snapshot
    //    bug between mesh mutation and per_fault population.
    {
        std::int64_t mesh_total = 0;
        for (const auto& m : meshes) {
            mesh_total += static_cast<std::int64_t>(m.number_of_faces());
        }
        std::int64_t pf_total = 0;
        for (const auto& f : per_fault) pf_total += f.n_triangles;
        if (mesh_total != pf_total) {
            std::cerr << "[corefine_faults] schema cross-check: "
                         "mesh total " << mesh_total
                      << " != per_fault sum " << pf_total << '\n';
            return 3;
        }
    }

    // 8. Schema-invariant check on triangle_to_fault.json (R-001):
    //    re-read the on-disk JSON and verify n_total_triangles +
    //    contiguous range partition.  Replaces the previous
    //    tautological loop that compared the same vector to itself.
    if (int rc = safs::io::validate_triangle_to_fault_json(
            tri2fault, per_fault); rc != 0) {
        return rc;
    }

    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "[corefine_faults] done in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     t1 - t0).count()
              << " ms ("
              << N << " faults, " << pair_stats.size() << " pairs).\n";
    return EXIT_SUCCESS;
}
