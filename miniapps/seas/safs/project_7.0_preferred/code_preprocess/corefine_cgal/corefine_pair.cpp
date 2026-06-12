// corefine_pair — single-pair CGAL corefine + Delaunay remesh + repair.
//
// Usage: corefine_pair IN_A IN_B OUT_A OUT_B
//        [--mesh-edge-size SIZE] (default 1500.0)
//        [--min-edge SIZE]       (default 100.0)
//        [--polyline-spacing S]  (default 0.5*mesh-edge-size, >= min-edge)
//        [--features-angle-bound DEG] (default 60.0)
//        [--max-iterations N]    (default 5)
//        [--verbose]
//
// Implements PLAN_cgal_corefine_multifault.md Phase 1.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/surface_Delaunay_remeshing.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "io_helpers.h"
#include "polyline_resample.h"
#include "quality_repair.h"

using K     = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Mesh  = CGAL::Surface_mesh<Point>;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace pp  = CGAL::parameters;

struct Args {
    std::string in_a, in_b, out_a, out_b;
    double mesh_edge_size = 1500.0;
    double min_edge       = 100.0;
    double polyline_spacing = -1.0;  // sentinel: derive at parse time
    double features_angle_bound = 60.0;
    unsigned max_iterations = 5;
    bool verbose = false;
};

static int usage(const char* argv0, int code) {
    std::cerr << "usage: " << argv0
              << " IN_A IN_B OUT_A OUT_B"
              << " [--mesh-edge-size SIZE]"
              << " [--min-edge SIZE]"
              << " [--polyline-spacing SPACING]"
              << " [--features-angle-bound DEG]"
              << " [--max-iterations N]"
              << " [--verbose]" << std::endl;
    return code;
}

static bool parse_double(const char* s, double& out) {
    char* end = nullptr;
    out = std::strtod(s, &end);
    return end != s;
}

static bool parse_uint(const char* s, unsigned& out) {
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (end == s || v < 0) return false;
    out = static_cast<unsigned>(v);
    return true;
}

static int parse(int argc, char** argv, Args& a) {
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        auto next = [&](double& x) {
            if (++i >= argc) { std::cerr << "missing value for " << s << "\n"; std::exit(2); }
            if (!parse_double(argv[i], x)) {
                std::cerr << "invalid double for " << s << ": " << argv[i] << "\n"; std::exit(2);
            }
        };
        auto next_u = [&](unsigned& x) {
            if (++i >= argc) { std::cerr << "missing value for " << s << "\n"; std::exit(2); }
            if (!parse_uint(argv[i], x)) {
                std::cerr << "invalid uint for " << s << ": " << argv[i] << "\n"; std::exit(2);
            }
        };
        if (s == "-h" || s == "--help") return usage(argv[0], 0);
        else if (s == "--mesh-edge-size")        next(a.mesh_edge_size);
        else if (s == "--min-edge")              next(a.min_edge);
        else if (s == "--polyline-spacing")      next(a.polyline_spacing);
        else if (s == "--features-angle-bound") next(a.features_angle_bound);
        else if (s == "--max-iterations")        next_u(a.max_iterations);
        else if (s == "--verbose")               a.verbose = true;
        else pos.push_back(s);
    }
    if (pos.size() != 4) return usage(argv[0], 2);
    a.in_a = pos[0]; a.in_b = pos[1]; a.out_a = pos[2]; a.out_b = pos[3];
    if (a.polyline_spacing < 0.0) {
        a.polyline_spacing = std::max(0.5 * a.mesh_edge_size, a.min_edge);
    } else if (a.polyline_spacing < a.min_edge) {
        a.polyline_spacing = a.min_edge;
    }
    return -1;
}

// Mark every edge of `m` whose endpoints both match a polyline vertex
// (within `tol`) and which is itself colinear with the polyline locally.
// Approximation: we flag every edge whose source or target is within
// `tol` of any polyline vertex.  This is conservative — extra edges may
// be marked constrained — which only protects MORE structure, never
// less.
static auto build_polyline_ecm(Mesh& m,
                               const std::vector<std::vector<Point>>& polylines,
                               double tol) {
    auto ecm = get(CGAL::dynamic_edge_property_t<bool>(), m);
    std::unordered_set<std::size_t> tagged_v;
    const double tol2 = tol * tol;
    for (auto v : m.vertices()) {
        const auto& p = m.point(v);
        for (const auto& pl : polylines) {
            for (const auto& q : pl) {
                if (CGAL::squared_distance(p, q) <= tol2) {
                    tagged_v.insert(v.idx());
                    goto next_vertex;
                }
            }
        }
        next_vertex:;
    }
    for (auto e : m.edges()) {
        const auto h = m.halfedge(e);
        const std::size_t s_idx = m.source(h).idx();
        const std::size_t t_idx = m.target(h).idx();
        const bool s_in = tagged_v.count(s_idx) > 0;
        const bool t_in = tagged_v.count(t_idx) > 0;
        put(ecm, e, s_in && t_in);
    }
    return ecm;
}

// Count A_out vertices that have a bit-identical (within tol_match)
// match in B_out.  surface_Delaunay_remeshing resamples the polyline
// internally; A_out and B_out share the function-resampled vertices,
// not the input polyline vertices.
static std::pair<std::size_t, double>
count_shared_vertices(const Mesh& m_a, const Mesh& m_b, double tol_match) {
    std::vector<Point> b_pts;
    b_pts.reserve(m_b.number_of_vertices());
    for (auto u : m_b.vertices()) b_pts.push_back(m_b.point(u));
    const double tol2 = tol_match * tol_match;
    std::size_t shared = 0;
    double max_match_diff = 0.0;
    for (auto v : m_a.vertices()) {
        const auto& pa = m_a.point(v);
        double best2 = std::numeric_limits<double>::infinity();
        for (const auto& pb : b_pts) {
            const double d2 = CGAL::squared_distance(pa, pb);
            if (d2 < best2) best2 = d2;
            if (d2 <= tol2) break;
        }
        if (best2 <= tol2) {
            ++shared;
            const double d = std::sqrt(best2);
            if (d > max_match_diff) max_match_diff = d;
        }
    }
    return {shared, max_match_diff};
}

int main(int argc, char** argv) {
    Args a;
    if (int r = parse(argc, argv, a); r >= 0) return r;

    Mesh A, B;
    if (!io_helpers::read_polygon_mesh_any(a.in_a, A)) {
        std::cerr << "error: cannot read " << a.in_a << "\n"; return 6;
    }
    if (!io_helpers::read_polygon_mesh_any(a.in_b, B)) {
        std::cerr << "error: cannot read " << a.in_b << "\n"; return 6;
    }
    if (a.verbose) {
        std::cerr << "loaded A: " << A.number_of_vertices() << " V, "
                  << A.number_of_faces() << " F\n";
        std::cerr << "loaded B: " << B.number_of_vertices() << " V, "
                  << B.number_of_faces() << " F\n";
    }

    // Self-intersection guard.
    if (PMP::does_self_intersect(A)) {
        if (a.verbose) std::cerr << "  A self-intersects: running autorefine\n";
        PMP::autorefine(A);
    }
    if (PMP::does_self_intersect(B)) {
        if (a.verbose) std::cerr << "  B self-intersects: running autorefine\n";
        PMP::autorefine(B);
    }

    // 4. Compute polylines (read-only).
    std::vector<std::vector<Point>> polylines;
    PMP::surface_intersection(A, B, std::back_inserter(polylines));
    if (polylines.empty()) {
        if (a.verbose) std::cerr << "  no intersection between A and B\n";
        if (!io_helpers::write_polygon_mesh_high_precision(a.out_a, A)) return 6;
        if (!io_helpers::write_polygon_mesh_high_precision(a.out_b, B)) return 6;
        return 2;
    }
    if (a.verbose) {
        std::size_t total_pts = 0;
        for (const auto& pl : polylines) total_pts += pl.size();
        std::cerr << "  intersection: " << polylines.size()
                  << " polyline(s), " << total_pts << " vertices\n";
    }

    // 5. Resample polylines.
    auto resampled = polyline_resample::resample_all(
        polylines, a.polyline_spacing, a.min_edge);
    if (resampled.empty()) {
        std::cerr << "error: all polylines collapsed during resample\n";
        return 5;
    }

    // 6. Corefine.
    auto ecmA = get(CGAL::dynamic_edge_property_t<bool>(), A);
    auto ecmB = get(CGAL::dynamic_edge_property_t<bool>(), B);
    PMP::corefine(A, B,
        pp::edge_is_constrained_map(ecmA),
        pp::edge_is_constrained_map(ecmB));

    // 7. Surface Delaunay remeshing on each, with the resampled polyline
    //    as the shared protected polyline_constraints.  We must set
    //    mesh_facet_size/angle/distance too — `mesh_edge_size` only
    //    constrains 1-D feature edges (the polyline); without facet
    //    criteria, the 2-D surface is left at unbounded element size
    //    (CGAL surface_Delaunay_remeshing.h:281, fsize default 0).
    const double facet_size     = a.mesh_edge_size;
    const double facet_angle    = 25.0;            // degrees, Mesh_3 default
    const double facet_distance = a.mesh_edge_size * 0.1;
    Mesh A_out = PMP::surface_Delaunay_remeshing<Mesh>(
        A,
        pp::polyline_constraints(resampled)
          .protect_constraints(true)
          .mesh_edge_size(a.mesh_edge_size)
          .mesh_facet_size(facet_size)
          .mesh_facet_angle(facet_angle)
          .mesh_facet_distance(facet_distance)
          .features_angle_bound(a.features_angle_bound));
    Mesh B_out = PMP::surface_Delaunay_remeshing<Mesh>(
        B,
        pp::polyline_constraints(resampled)
          .protect_constraints(true)
          .mesh_edge_size(a.mesh_edge_size)
          .mesh_facet_size(facet_size)
          .mesh_facet_angle(facet_angle)
          .mesh_facet_distance(facet_distance)
          .features_angle_bound(a.features_angle_bound));

    // 8. Re-derive constraint maps on the remeshed outputs.
    auto ecmA_out = build_polyline_ecm(A_out, resampled, 1e-9);
    auto ecmB_out = build_polyline_ecm(B_out, resampled, 1e-9);

    // 9. Quality repair.  surface_Delaunay_remeshing already protects
    //    polyline endpoints intrinsically (protect_constraints=true), so
    //    an empty vertex constraint map is fine here.  corefine_set.cpp
    //    uses isotropic_remeshing instead and must build a real VCM.
    auto vcmA_out = quality_repair::make_empty_vcm(A_out);
    auto vcmB_out = quality_repair::make_empty_vcm(B_out);
    auto stA = quality_repair::run(A_out, ecmA_out, vcmA_out, a.min_edge);
    auto stB = quality_repair::run(B_out, ecmB_out, vcmB_out, a.min_edge);

    // 10. Conformality check (A_out ↔ B_out, not against input polyline).
    auto [shared_a, drift_a] = count_shared_vertices(A_out, B_out, 1e-6);
    auto [shared_b, drift_b] = count_shared_vertices(B_out, A_out, 1e-6);
    const bool conformal = (shared_a == shared_b)
                        && (shared_a > 0)
                        && (drift_a < 1e-6) && (drift_b < 1e-6);
    if (!conformal) {
        std::cerr << "error: pairwise conformality lost: shared_a="
                  << shared_a << " shared_b=" << shared_b
                  << " drift_a=" << drift_a << " drift_b=" << drift_b << "\n";
        return 3;
    }

    // 11. Hard-floor validation (warn, don't fail — the polyline floor is
    //     enforced by resample but corefine can produce sub-floor edges
    //     that survive remesh; report them in stats).
    if (stA.n_edges_below_floor > 0 || stB.n_edges_below_floor > 0) {
        std::cerr << "warning: " << stA.n_edges_below_floor
                  << "+" << stB.n_edges_below_floor
                  << " edges below min-edge floor (" << a.min_edge << " m)\n";
    }

    // 12. Write outputs.
    if (!io_helpers::write_polygon_mesh_high_precision(a.out_a, A_out)) {
        std::cerr << "error: cannot write " << a.out_a << "\n"; return 6;
    }
    if (!io_helpers::write_polygon_mesh_high_precision(a.out_b, B_out)) {
        std::cerr << "error: cannot write " << a.out_b << "\n"; return 6;
    }

    // 13. Print stats block.
    std::size_t total_resampled = 0;
    for (const auto& pl : resampled) total_resampled += pl.size();
    std::printf(
        "A:  V_in=%zu  F_in=%zu  V_out=%zu  F_out=%zu\n"
        "    edge:  min=%.3f  p1=%.3f  median=%.3f  max=%.3f\n"
        "    q_tri: min=%.4f  p1=%.4f  median=%.4f\n",
        static_cast<std::size_t>(A.number_of_vertices()),
        static_cast<std::size_t>(A.number_of_faces()),
        static_cast<std::size_t>(stA.n_verts),
        static_cast<std::size_t>(stA.n_faces),
        stA.edge_len_min, stA.edge_len_p1, stA.edge_len_med, stA.edge_len_max,
        stA.tri_q_min, stA.tri_q_p1, stA.tri_q_med);
    std::printf(
        "B:  V_in=%zu  F_in=%zu  V_out=%zu  F_out=%zu\n"
        "    edge:  min=%.3f  p1=%.3f  median=%.3f  max=%.3f\n"
        "    q_tri: min=%.4f  p1=%.4f  median=%.4f\n",
        static_cast<std::size_t>(B.number_of_vertices()),
        static_cast<std::size_t>(B.number_of_faces()),
        static_cast<std::size_t>(stB.n_verts),
        static_cast<std::size_t>(stB.n_faces),
        stB.edge_len_min, stB.edge_len_p1, stB.edge_len_med, stB.edge_len_max,
        stB.tri_q_min, stB.tri_q_p1, stB.tri_q_med);
    std::printf("Polyline:  count=%zu  resampled_count=%zu\n"
                "Conformality: %zu shared verts in A↔B (max match diff %.3e m),\n"
                "              %zu shared verts in B↔A (max match diff %.3e m).\n",
                polylines.size(), total_resampled,
                shared_a, drift_a, shared_b, drift_b);
    std::printf("Quality floor passed: %s\n",
                (stA.n_edges_below_floor == 0 && stB.n_edges_below_floor == 0
                 && stA.tri_q_min >= 0.3 && stB.tri_q_min >= 0.3)
                ? "yes" : "no");
    return 0;
}
