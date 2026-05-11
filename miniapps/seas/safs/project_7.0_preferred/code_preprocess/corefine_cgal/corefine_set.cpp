// corefine_set — N-fault driver: pairwise corefine + per-mesh Delaunay remesh.
//
// Usage: corefine_set IN_DIR OUT_DIR
//        [--ext .off]
//        [--mesh-edge-size SIZE]    (default 1500.0)
//        [--min-edge SIZE]          (default 100.0)
//        [--polyline-spacing S]     (default 0.5*mesh-edge-size, >= min-edge)
//        [--features-angle-bound DEG] (default 60.0)
//        [--manifest path.json]     (default OUT_DIR/manifest.json)
//        [--verbose]
//
// Implements PLAN_cgal_corefine_multifault.md Phase 1 (multi-fault).

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <boost/property_map/property_map.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "io_helpers.h"
#include "polyline_resample.h"
#include "quality_repair.h"
#include "intersection_graph.h"
#include "polyline_cleanup.h"

#include <array>
#include <functional>
#include <map>
#include <numeric>
#include <set>

using K     = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Mesh  = CGAL::Surface_mesh<Point>;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace pp  = CGAL::parameters;
namespace fs  = std::filesystem;

struct Args {
    std::string in_dir, out_dir;
    std::string ext = ".off";
    std::string manifest_path;
    double mesh_edge_size = 1500.0;
    double min_edge       = 100.0;
    double polyline_spacing = -1.0;
    double features_angle_bound = 60.0;
    bool verbose = false;
};

static int usage(const char* argv0, int code) {
    std::cerr << "usage: " << argv0 << " IN_DIR OUT_DIR"
              << " [--ext .off] [--mesh-edge-size SIZE] [--min-edge SIZE]"
              << " [--polyline-spacing S] [--features-angle-bound DEG]"
              << " [--manifest PATH] [--verbose]" << std::endl;
    return code;
}

static int parse(int argc, char** argv, Args& a) {
    std::vector<std::string> pos;
    auto need_double = [&](int& i, double& x, const std::string& name) {
        if (++i >= argc) { std::cerr << "missing value for " << name << "\n"; std::exit(2); }
        char* end = nullptr;
        x = std::strtod(argv[i], &end);
        if (end == argv[i]) {
            std::cerr << "invalid double for " << name << ": " << argv[i] << "\n";
            std::exit(2);
        }
    };
    auto need_str = [&](int& i, std::string& x, const std::string& name) {
        if (++i >= argc) { std::cerr << "missing value for " << name << "\n"; std::exit(2); }
        x = argv[i];
    };
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "-h" || s == "--help") return usage(argv[0], 0);
        else if (s == "--ext")              need_str(i, a.ext, s);
        else if (s == "--mesh-edge-size")   need_double(i, a.mesh_edge_size, s);
        else if (s == "--min-edge")         need_double(i, a.min_edge, s);
        else if (s == "--polyline-spacing") need_double(i, a.polyline_spacing, s);
        else if (s == "--features-angle-bound") need_double(i, a.features_angle_bound, s);
        else if (s == "--manifest")         need_str(i, a.manifest_path, s);
        else if (s == "--verbose")          a.verbose = true;
        else pos.push_back(s);
    }
    if (pos.size() != 2) return usage(argv[0], 2);
    a.in_dir = pos[0]; a.out_dir = pos[1];
    if (a.polyline_spacing < 0.0) {
        a.polyline_spacing = std::max(0.5 * a.mesh_edge_size, a.min_edge);
    } else if (a.polyline_spacing < a.min_edge) {
        a.polyline_spacing = a.min_edge;
    }
    if (a.manifest_path.empty()) {
        a.manifest_path = (fs::path(a.out_dir) / "manifest.json").string();
    }
    return -1;
}

// Mark every edge of `m` whose both endpoints match a polyline vertex
// within `tol`.  See note in corefine_pair.cpp (conservative).
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
                    goto next_v;
                }
            }
        }
        next_v:;
    }
    for (auto e : m.edges()) {
        const auto h = m.halfedge(e);
        const bool s_in = tagged_v.count(m.source(h).idx()) > 0;
        const bool t_in = tagged_v.count(m.target(h).idx()) > 0;
        put(ecm, e, s_in && t_in);
    }
    return ecm;
}

// Count A_out vertices that have a bit-identical (within tol_match) match
// in B_out.  surface_Delaunay_remeshing resamples polyline_constraints
// internally, so the input polyline vertices are NOT preserved verbatim;
// what matters for gmsh dedup is that the function-resampled polyline
// vertices are SHARED between A_out and B_out.  This test counts that
// shared subset directly.
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

// Snap polyline vertices within `snap_radius` to their cluster centroid.
// At triple junctions (e.g. mesh A intersected with both B and C), the
// two surface_intersection calls compute the A-side endpoint
// independently and may disagree by ~mm-m due to exact-predicate
// rounding.  Without snapping, surface_Delaunay_remeshing protects each
// near-coincident vertex as a separate feature and creates a fan of
// needle triangles around the cluster (edge_min ~ 0.1 m).  Snapping
// merges the cluster into one common coord.
//
// Implementation: union-find on point pairs with squared distance ≤
// snap_radius^2; replace each clustered vertex by the cluster centroid.
static void snap_polyline_clusters(std::vector<std::vector<Point>>& polylines,
                                    double snap_radius) {
    struct Entry { std::size_t li, vi; };
    std::vector<Entry> entries;
    std::vector<Point> points;
    for (std::size_t li = 0; li < polylines.size(); ++li) {
        for (std::size_t vi = 0; vi < polylines[li].size(); ++vi) {
            entries.push_back({li, vi});
            points.push_back(polylines[li][vi]);
        }
    }
    const std::size_t N = entries.size();
    if (N < 2) return;

    std::vector<std::size_t> parent(N);
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    std::function<std::size_t(std::size_t)> find = [&](std::size_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    const double r2 = snap_radius * snap_radius;
    for (std::size_t a = 0; a < N; ++a) {
        for (std::size_t b = a + 1; b < N; ++b) {
            if (CGAL::squared_distance(points[a], points[b]) <= r2) {
                const std::size_t ra = find(a);
                const std::size_t rb = find(b);
                if (ra != rb) parent[ra] = rb;
            }
        }
    }

    // Per-cluster centroid.
    std::map<std::size_t, std::pair<int, std::array<double, 3>>> clusters;
    for (std::size_t i = 0; i < N; ++i) {
        const std::size_t root = find(i);
        auto& [count, sum] = clusters[root];
        ++count;
        sum[0] += points[i].x();
        sum[1] += points[i].y();
        sum[2] += points[i].z();
    }

    // Replace vertices in clusters of size > 1.
    std::size_t n_snapped = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const std::size_t root = find(i);
        const auto& [count, sum] = clusters[root];
        if (count > 1) {
            const double cx = sum[0] / count;
            const double cy = sum[1] / count;
            const double cz = sum[2] / count;
            polylines[entries[i].li][entries[i].vi] = Point(cx, cy, cz);
            ++n_snapped;
        }
    }
    (void)n_snapped;  // available for diagnostic logging if desired
}

// JSON formatting helpers (no external dep).
static void json_str(std::ostream& os, const std::string& s) {
    os << '"';
    for (char c : s) {
        if (c == '"' || c == '\\') os << '\\' << c;
        else if (c == '\n') os << "\\n";
        else os << c;
    }
    os << '"';
}

static void json_double(std::ostream& os, double x) {
    if (std::isnan(x) || std::isinf(x)) os << "null";
    else os << std::setprecision(15) << x;
}

int main(int argc, char** argv) {
    Args a;
    if (int r = parse(argc, argv, a); r >= 0) return r;

    fs::create_directories(a.out_dir);

    // 1. Load all .off files in IN_DIR (alphabetical).
    std::vector<std::string> basenames;
    std::vector<std::string> input_paths;
    for (const auto& entry : fs::directory_iterator(a.in_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == a.ext) {
            input_paths.push_back(entry.path().string());
        }
    }
    std::sort(input_paths.begin(), input_paths.end());
    if (input_paths.empty()) {
        std::cerr << "no inputs in " << a.in_dir << " with ext " << a.ext << "\n";
        return 1;
    }

    std::vector<Mesh> meshes(input_paths.size());
    std::vector<std::size_t> v_in(input_paths.size(), 0);
    std::vector<std::size_t> f_in(input_paths.size(), 0);
    for (std::size_t i = 0; i < input_paths.size(); ++i) {
        if (!io_helpers::read_polygon_mesh_any(input_paths[i], meshes[i])) {
            std::cerr << "error: cannot read " << input_paths[i] << "\n"; return 6;
        }
        basenames.push_back(fs::path(input_paths[i]).stem().string());
        v_in[i] = meshes[i].number_of_vertices();
        f_in[i] = meshes[i].number_of_faces();
        if (a.verbose) {
            std::cerr << "[" << i << "] " << basenames[i]
                      << "  V=" << v_in[i] << "  F=" << f_in[i] << "\n";
        }
        if (PMP::does_self_intersect(meshes[i])) {
            if (a.verbose) std::cerr << "    self-intersects: autorefine\n";
            PMP::autorefine(meshes[i]);
        }
    }

    // 2-3. Discover and order intersecting pairs.
    auto pairs = intersection_graph::find_intersecting_pairs(meshes);
    intersection_graph::sort_by_overlap_volume(pairs, meshes);
    if (a.verbose) {
        std::cerr << "found " << pairs.size() << " intersecting pair(s)\n";
        for (const auto& p : pairs) {
            std::cerr << "  (" << p.first << ", " << p.second << ")  bbox-vol="
                      << intersection_graph::bbox_overlap_volume(meshes[p.first], meshes[p.second])
                      << "\n";
        }
    }

    // 4-5. Per-pair: surface_intersection -> resample -> corefine.
    //      Accumulate per-mesh polyline sets and per-pair polylines for manifest.
    //
    //      Persistent constraint maps: corefine writes `true` into ecm_i
    //      and ecm_j for every newly-inserted polyline edge.  We need
    //      ONE map per mesh that accumulates constraints across all
    //      corefine calls involving that mesh.
    //
    //      Use Surface_mesh's internal property_map<edge_descriptor, bool>
    //      with default value `false`.  Critical: a property map that
    //      asserts/UB on missing-key reads (like
    //      boost::associative_property_map<std::map>) silently breaks
    //      isotropic_remeshing's "do not collapse edge with two
    //      constrained vertices" check (remesh_impl.h:701-705): the
    //      check reads vertex_is_constrained_map for vertices that may
    //      not be in our explicit list, gets undefined values in
    //      release builds, and skips protection randomly.  Surface_mesh
    //      property maps return the documented default for unknown
    //      keys.
    using EdgeDesc   = boost::graph_traits<Mesh>::edge_descriptor;
    using VertexDesc = boost::graph_traits<Mesh>::vertex_descriptor;
    using ECMProp = Mesh::Property_map<EdgeDesc,   bool>;
    using VCMProp = Mesh::Property_map<VertexDesc, bool>;
    std::vector<ECMProp> ecm_maps(meshes.size());
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        bool ok;
        std::tie(ecm_maps[i], ok) =
            meshes[i].template add_property_map<EdgeDesc, bool>(
                "e:is_corefine_constrained", false);
        (void)ok;
    }

    std::vector<std::vector<std::vector<Point>>> per_mesh_polylines(meshes.size());
    std::vector<std::vector<std::vector<Point>>> per_pair_polylines(pairs.size());

    for (std::size_t k = 0; k < pairs.size(); ++k) {
        const auto [i, j] = pairs[k];
        std::vector<std::vector<Point>> raw_polylines;
        PMP::surface_intersection(meshes[i], meshes[j],
                                  std::back_inserter(raw_polylines));
        if (raw_polylines.empty()) {
            if (a.verbose) std::cerr << "  pair (" << i << "," << j << "): empty\n";
            continue;
        }
        auto resampled = polyline_resample::resample_all(
            raw_polylines, a.polyline_spacing, a.min_edge);
        if (resampled.empty()) {
            if (a.verbose) std::cerr << "  pair (" << i << "," << j << "): resample collapsed\n";
            continue;
        }
        per_pair_polylines[k] = resampled;

        PMP::corefine(meshes[i], meshes[j],
                      pp::edge_is_constrained_map(ecm_maps[i]),
                      pp::edge_is_constrained_map(ecm_maps[j]));

        for (const auto& pl : resampled) {
            per_mesh_polylines[i].push_back(pl);
            per_mesh_polylines[j].push_back(pl);
        }
        if (a.verbose) {
            std::size_t npts = 0;
            for (const auto& pl : resampled) npts += pl.size();
            std::size_t n_constrained_i = 0;
            for (auto e : meshes[i].edges()) if (get(ecm_maps[i], e)) ++n_constrained_i;
            std::size_t n_constrained_j = 0;
            for (auto e : meshes[j].edges()) if (get(ecm_maps[j], e)) ++n_constrained_j;
            std::cerr << "  pair (" << i << "," << j << "): "
                      << resampled.size() << " polyline(s), " << npts << " verts"
                      << " | constrained edges in mesh " << i << ": " << n_constrained_i
                      << ", in mesh " << j << ": " << n_constrained_j << "\n";
        }
    }

    // 5.5. Deterministic shared-polyline cleanup: collapse polyline edges
    //      shorter than min_edge in EVERY mesh that contains them, with
    //      bit-identical decisions and midpoint placement so pairwise
    //      conformality is preserved.  See polyline_cleanup.h for the full
    //      rationale.  Corrects the symptom that PMP::corefine inserts the
    //      raw exact-predicate polyline (which can have 0.1-m segments at
    //      triple junctions or near-vertex hits) and protect_constraints
    //      =true in the downstream isotropic_remeshing prevents those edges
    //      from being collapsed by the remesher.
    polyline_cleanup::Stats cleanup =
        polyline_cleanup::collapse_short_polyline_edges_pairwise(
            meshes, ecm_maps, a.min_edge);
    if (a.verbose) {
        std::cerr << "polyline cleanup: " << cleanup.n_collapses_done
                  << " collapse(s) over " << cleanup.n_iterations
                  << " iter(s)";
        if (cleanup.n_residual_short > 0) {
            std::cerr << "  RESIDUAL " << cleanup.n_residual_short
                      << " short polyline edge(s), shortest = "
                      << cleanup.residual_min_length << " m";
        }
        std::cerr << "\n";
    }
    if (cleanup.n_residual_short > 0) {
        std::cerr << "ERROR: " << cleanup.n_residual_short
                  << " polyline edge(s) below min_edge=" << a.min_edge
                  << " m could not be collapsed (link condition failed in"
                     " at least one mesh of the pair).  Shortest residual: "
                  << cleanup.residual_min_length << " m. Aborting to"
                     " preserve the hard-constraint contract — a sliver-"
                     "tolerant pipeline must use a smaller min_edge.\n";
        return 8;
    }

    // 5.6. Cross-polyline cluster snap.  Phase 1 (cleanup above) handles
    //      sub-floor edges within a single polyline.  At triple junctions,
    //      polyline endpoints from DIFFERENT polylines (e.g. (0,3) and
    //      (0,4)) can sit within min_edge of each other; isotropic_remeshing
    //      then meshes between them and creates an unconstrained sub-floor
    //      edge that downstream quality_repair cannot collapse.  This pass
    //      snaps every cross-polyline cluster of polyline vertices within
    //      min_edge to the cluster centroid in EVERY mesh that touches it,
    //      preserving conformality by construction.
    polyline_cleanup::ClusterStats snap_stats =
        polyline_cleanup::snap_cross_polyline_clusters_pairwise(
            meshes, ecm_maps, a.min_edge);
    if (a.verbose) {
        std::cerr << "cross-polyline snap: "
                  << snap_stats.n_clusters_processed << " cluster(s), "
                  << snap_stats.n_vertices_snapped << " vert snap(s), "
                  << snap_stats.n_intra_mesh_collapses
                  << " intra-mesh collapse(s) over "
                  << snap_stats.n_iterations << " iter(s)\n";
    }

    // 6. Per-mesh remesh via PMP::isotropic_remeshing (NOT
    //    surface_Delaunay_remeshing).  Reasons:
    //      * isotropic_remeshing constrains borders by default, so the
    //        fault perimeter is preserved exactly — eliminates the
    //        zig-zag we observed on the San_Andreas fault output (a
    //        symptom of surface_Delaunay_remeshing's free-border
    //        reconstruction within facet_distance).
    //      * isotropic_remeshing modifies the input in-place via VCG
    //        split/collapse/swap/smooth; doesn't rebuild from scratch
    //        via Mesh_3.  Orders of magnitude faster on these inputs
    //        (the surface_Delaunay_remeshing path was hanging > 15 min
    //        with hundreds of protected polyline vertices).
    //      * The polyline edges from the corefine step are passed via
    //        edge_is_constrained_map (built from per_mesh_polylines via
    //        proximity match in build_polyline_ecm).  We use
    //        protect_constraints=false / collapse_constraints=true so
    //        polyline sub-edges shorter than 4/5 * mesh_edge_size get
    //        collapsed — handles the triple-junction needle clusters
    //        without losing the polyline as a topological feature.
    //    Use the PERSISTENT corefine constraint map (ecm_maps[i]) — NOT
    //    a proximity-rebuilt map.  Reason: PMP::corefine inserts the RAW
    //    intersection polyline (computed by exact predicates) into the
    //    mesh, marking those edges true in ecm_i.  The RESAMPLED
    //    polyline used elsewhere is at uniform arc-length spacing and
    //    DOES NOT coincide with the raw polyline vertices.  An ECM
    //    rebuilt by proximity to the resampled polyline misses most of
    //    the actual corefine-inserted edges, leaving them unprotected
    //    and isotropic_remeshing collapses them away — that was the
    //    pair (3,5) conformality bug.
    //
    //    NB: meshes[] still holds the post-corefine state at this
    //    point.  We must NOT copy meshes[i] to output_meshes[i] before
    //    isotropic_remeshing because the EdgeDesc keys in ecm_maps[i]
    //    refer to meshes[i]'s edges, not output_meshes[i]'s.  Operate
    //    directly in-place on meshes[i] and move into output_meshes
    //    afterward.

    // DEBUG: dump the constrained edges (post-corefine, pre-remesh) of
    // each mesh as a JSON file with (x,y,z) endpoints.  Used to verify
    // that corefine produced identical polyline edges in both meshes
    // of every intersecting pair, and that protect_constraints in
    // isotropic_remeshing preserves them.
    if (a.verbose) {
        std::ofstream dump((fs::path(a.out_dir) / "debug_constrained_edges.json"));
        dump << std::setprecision(15);
        dump << "{\n  \"meshes\": [\n";
        for (std::size_t i = 0; i < meshes.size(); ++i) {
            dump << "    {\n      \"index\": " << i
                 << ",\n      \"basename\": \"" << basenames[i]
                 << "\",\n      \"constrained_edges_pre_remesh\": [\n";
            const Mesh& m = meshes[i];
            std::size_t cnt = 0;
            for (auto e : m.edges()) {
                if (!get(ecm_maps[i], e)) continue;
                const auto h = m.halfedge(e);
                const auto& p = m.point(m.source(h));
                const auto& q = m.point(m.target(h));
                if (cnt++ > 0) dump << ",\n";
                dump << "        [["<<p.x()<<", "<<p.y()<<", "<<p.z()<<"], ["<<q.x()<<", "<<q.y()<<", "<<q.z()<<"]]";
            }
            dump << "\n      ]\n    }" << (i+1 < meshes.size() ? "," : "") << "\n";
        }
        dump << "  ]\n}\n";
    }

    std::vector<Mesh> output_meshes(meshes.size());
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        if (per_mesh_polylines[i].empty()) {
            output_meshes[i] = meshes[i];
            if (a.verbose) std::cerr << "[" << i << "] no intersections; no remesh\n";
            continue;
        }

        // Build vertex_is_constrained_map by marking every endpoint of
        // every constrained edge.  Without this, isotropic_remeshing's
        // edge-collapse step (remesh_impl.h:701-705) treats polyline
        // vertices as unconstrained and may collapse an unconstrained
        // edge incident on them — geometrically MOVING the polyline
        // endpoint.  When this decision differs between mesh A and B,
        // the polyline endpoint diverges and conformality is lost.
        VCMProp vc_map;
        bool ok;
        std::tie(vc_map, ok) =
            meshes[i].template add_property_map<VertexDesc, bool>(
                "v:is_corefine_constrained", false);
        (void)ok;
        std::size_t n_constrained_e = 0, n_constrained_v = 0;
        for (auto e : meshes[i].edges()) {
            if (!get(ecm_maps[i], e)) continue;
            ++n_constrained_e;
            const auto h = meshes[i].halfedge(e);
            put(vc_map, meshes[i].source(h), true);
            put(vc_map, meshes[i].target(h), true);
        }
        for (auto v : meshes[i].vertices()) if (get(vc_map, v)) ++n_constrained_v;

        if (a.verbose) {
            std::cerr << "[" << i << "] isotropic_remeshing target="
                      << a.mesh_edge_size << "m  constrained_edges="
                      << n_constrained_e << "  constrained_verts=" << n_constrained_v
                      << "  polylines=" << per_mesh_polylines[i].size() << "\n";
        }
        PMP::isotropic_remeshing(
            faces(meshes[i]),
            a.mesh_edge_size,
            meshes[i],
            pp::edge_is_constrained_map(ecm_maps[i])
              .vertex_is_constrained_map(vc_map)
              .protect_constraints(true)
              .number_of_iterations(5));
        output_meshes[i] = std::move(meshes[i]);
    }

    // 7. Quality repair on each output: removes any free unconstrained
    //    needle/cap triangles left over from isotropic_remeshing's
    //    protect_constraints=true squeeze.  CRITICAL: pass the per-mesh
    //    vertex_is_constrained_map so polyline endpoints are not moved
    //    by remove_almost_degenerate_faces (collapsing an unconstrained
    //    edge incident on a polyline vertex would geometrically move the
    //    polyline vertex and break pairwise conformality across the pair).
    //
    //    Property maps must be re-acquired on output_meshes[i] by name
    //    because the std::move above invalidated the handles cached in
    //    ecm_maps[i].
    std::vector<quality_repair::Stats> stats(meshes.size());
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        if (per_mesh_polylines[i].empty()) {
            // No intersections: no polyline to protect.  Use empty maps.
            auto vcm_empty = quality_repair::make_empty_vcm(output_meshes[i]);
            auto ecm_empty = get(CGAL::dynamic_edge_property_t<bool>(),
                                 output_meshes[i]);
            for (auto e : output_meshes[i].edges()) put(ecm_empty, e, false);
            stats[i] = quality_repair::run(output_meshes[i], ecm_empty,
                                            vcm_empty, a.min_edge);
            continue;
        }
        auto ecm_out_opt = output_meshes[i].template
            property_map<EdgeDesc, bool>("e:is_corefine_constrained");
        if (!ecm_out_opt.has_value()) {
            std::cerr << "internal error: ECM lost on output_meshes["
                      << i << "]\n";
            return 9;
        }
        auto ecm_out = *ecm_out_opt;
        // Build VCM from constrained-edge endpoints in the output mesh.
        auto vcm_out = quality_repair::make_empty_vcm(output_meshes[i]);
        for (auto e : output_meshes[i].edges()) {
            if (!get(ecm_out, e)) continue;
            auto h = output_meshes[i].halfedge(e);
            put(vcm_out, output_meshes[i].source(h), true);
            put(vcm_out, output_meshes[i].target(h), true);
        }
        stats[i] = quality_repair::run(output_meshes[i], ecm_out, vcm_out,
                                        a.min_edge);
    }

    // 8. Pairwise conformality re-check (A_out ↔ B_out, not against
    //    input polyline — see count_shared_vertices comment).
    struct PairResult {
        std::size_t i, j, n_polyline_verts_in;
        std::size_t shared_in_a, shared_in_b;
        bool conformal;
        double max_diff;
    };
    std::vector<PairResult> pair_results;
    pair_results.reserve(pairs.size());
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        if (per_pair_polylines[k].empty()) continue;
        const auto [i, j] = pairs[k];
        std::size_t n_in = 0;
        for (const auto& pl : per_pair_polylines[k]) n_in += pl.size();
        // tol_match = 1e-6 m absolute (R-009); 1000x tighter than gmsh
        // Geometry.Tolerance = 1e-3.
        auto [shared_a, d_a] = count_shared_vertices(output_meshes[i],
                                                      output_meshes[j],
                                                      /*tol_match=*/1e-6);
        auto [shared_b, d_b] = count_shared_vertices(output_meshes[j],
                                                      output_meshes[i],
                                                      /*tol_match=*/1e-6);
        // Symmetric and at least one shared vertex — that's the contract.
        const bool conformal = (shared_a == shared_b)
                            && (shared_a > 0)
                            && (d_a < 1e-6) && (d_b < 1e-6);
        pair_results.push_back({i, j, n_in, shared_a, shared_b,
                                conformal, std::max(d_a, d_b)});
        if (a.verbose) {
            std::cerr << "  pair (" << i << "," << j << "): "
                      << (conformal ? "CONFORMAL" : "NON-CONFORMAL")
                      << "  shared=" << shared_a << "/" << shared_b
                      << "  polyline_in=" << n_in
                      << "  max_diff=" << std::max(d_a, d_b) << "\n";
        }
    }

    // 8.5. Final hard-floor check on the OUTPUT meshes.  The polyline
    //      cleanup (step 5.5) zeroes out sub-floor *polyline* edges, but
    //      isotropic_remeshing can introduce sub-floor *unconstrained*
    //      edges near triple junctions where two polyline endpoints from
    //      DIFFERENT polylines lie within min_edge of each other (the
    //      vcm protection on both endpoints prevents quality_repair from
    //      collapsing them, and resolving them would require coordinating
    //      a vertex move across ≥3 meshes).  Hard-fail per project policy
    //      so the user is informed of the contract violation.
    std::size_t total_below = 0;
    double      worst       = std::numeric_limits<double>::infinity();
    std::size_t worst_mesh  = 0;
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        if (stats[i].n_edges_below_floor > 0) {
            total_below += stats[i].n_edges_below_floor;
            if (stats[i].edge_len_min < worst) {
                worst      = stats[i].edge_len_min;
                worst_mesh = i;
            }
        }
    }
    if (total_below > 0) {
        std::cerr << "ERROR: " << total_below
                  << " output edge(s) below min_edge=" << a.min_edge
                  << " m after cleanup + isotropic_remeshing + quality_"
                     "repair.  Worst: " << worst << " m in mesh "
                  << worst_mesh << " (" << basenames[worst_mesh] << ").\n"
                  << "  These are unconstrained edges between polyline "
                     "endpoints from DIFFERENT polylines (triple junctions). "
                     "The current shared-polyline cleanup does not snap "
                     "across-polyline vertex clusters. Aborting.\n";
        return 8;
    }

    // 9. Write outputs and the manifest.
    std::vector<std::string> output_paths(meshes.size());
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        output_paths[i] = (fs::path(a.out_dir) /
                           (basenames[i] + "_corefined" + a.ext)).string();
        if (!io_helpers::write_polygon_mesh_high_precision(
                output_paths[i], output_meshes[i])) {
            std::cerr << "error: cannot write " << output_paths[i] << "\n";
            return 6;
        }
    }

    std::ofstream mfest(a.manifest_path);
    if (!mfest) {
        std::cerr << "error: cannot write manifest " << a.manifest_path << "\n";
        return 6;
    }
    mfest << "{\n";
    mfest << "  \"ext\": "; json_str(mfest, a.ext); mfest << ",\n";
    mfest << "  \"mesh_edge_size\": "; json_double(mfest, a.mesh_edge_size); mfest << ",\n";
    mfest << "  \"min_edge\": ";       json_double(mfest, a.min_edge);       mfest << ",\n";
    mfest << "  \"polyline_spacing\": "; json_double(mfest, a.polyline_spacing); mfest << ",\n";
    mfest << "  \"features_angle_bound\": "; json_double(mfest, a.features_angle_bound); mfest << ",\n";
    mfest << "  \"polyline_cleanup\": {";
    mfest << " \"n_iterations\": " << cleanup.n_iterations;
    mfest << ", \"n_collapses_done\": " << cleanup.n_collapses_done;
    mfest << ", \"n_residual_short\": " << cleanup.n_residual_short;
    mfest << ", \"residual_min_length_m\": "; json_double(mfest, cleanup.residual_min_length);
    mfest << " },\n";
    mfest << "  \"cross_polyline_snap\": {";
    mfest << " \"n_iterations\": " << snap_stats.n_iterations;
    mfest << ", \"n_clusters_processed\": " << snap_stats.n_clusters_processed;
    mfest << ", \"n_vertices_snapped\": " << snap_stats.n_vertices_snapped;
    mfest << ", \"n_intra_mesh_collapses\": " << snap_stats.n_intra_mesh_collapses;
    mfest << " },\n";
    mfest << "  \"meshes\": [\n";
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        const auto& s = stats[i];
        // Compute output bbox.
        auto bb = PMP::bbox(output_meshes[i]);
        mfest << "    {\n";
        mfest << "      \"index\": " << i << ",\n";
        mfest << "      \"input\": ";  json_str(mfest, fs::path(input_paths[i]).filename().string());  mfest << ",\n";
        mfest << "      \"output\": "; json_str(mfest, fs::path(output_paths[i]).filename().string()); mfest << ",\n";
        mfest << "      \"basename\": "; json_str(mfest, basenames[i]); mfest << ",\n";
        mfest << "      \"n_verts_in\": " << v_in[i] << ",\n";
        mfest << "      \"n_faces_in\": " << f_in[i] << ",\n";
        mfest << "      \"n_verts\": " << s.n_verts << ",\n";
        mfest << "      \"n_faces\": " << s.n_faces << ",\n";
        mfest << "      \"bbox\": {\n";
        mfest << "        \"x_lo\": "; json_double(mfest, bb.xmin()); mfest << ", \"x_hi\": "; json_double(mfest, bb.xmax()); mfest << ",\n";
        mfest << "        \"y_lo\": "; json_double(mfest, bb.ymin()); mfest << ", \"y_hi\": "; json_double(mfest, bb.ymax()); mfest << ",\n";
        mfest << "        \"z_lo\": "; json_double(mfest, bb.zmin()); mfest << ", \"z_hi\": "; json_double(mfest, bb.zmax()); mfest << "\n";
        mfest << "      },\n";
        mfest << "      \"edge_stats\": {\"min\": "; json_double(mfest, s.edge_len_min);
        mfest << ", \"p1\": ";   json_double(mfest, s.edge_len_p1);
        mfest << ", \"median\": "; json_double(mfest, s.edge_len_med);
        mfest << ", \"max\": ";   json_double(mfest, s.edge_len_max); mfest << "},\n";
        mfest << "      \"tri_q_stats\": {\"min\": "; json_double(mfest, s.tri_q_min);
        mfest << ", \"p1\": ";   json_double(mfest, s.tri_q_p1);
        mfest << ", \"median\": "; json_double(mfest, s.tri_q_med); mfest << "},\n";
        mfest << "      \"n_edges_below_floor\": " << s.n_edges_below_floor << ",\n";
        mfest << "      \"n_tri_q_below_threshold\": " << s.n_tri_q_below_threshold << "\n";
        mfest << "    }" << (i + 1 < meshes.size() ? "," : "") << "\n";
    }
    mfest << "  ],\n";
    mfest << "  \"pairs\": [\n";
    for (std::size_t k = 0; k < pair_results.size(); ++k) {
        const auto& pr = pair_results[k];
        mfest << "    {\n";
        mfest << "      \"i\": " << pr.i << ",\n";
        mfest << "      \"j\": " << pr.j << ",\n";
        mfest << "      \"n_polyline_verts_in\": " << pr.n_polyline_verts_in << ",\n";
        mfest << "      \"shared_a\": " << pr.shared_in_a << ",\n";
        mfest << "      \"shared_b\": " << pr.shared_in_b << ",\n";
        mfest << "      \"conformal\": " << (pr.conformal ? "true" : "false") << ",\n";
        mfest << "      \"max_diff_m\": "; json_double(mfest, pr.max_diff); mfest << ",\n";
        mfest << "      \"polylines\": [\n";
        // Find the matching per-pair polyline list by (i,j).
        std::size_t pair_idx = 0;
        for (std::size_t kk = 0; kk < pairs.size(); ++kk) {
            if (pairs[kk].first == pr.i && pairs[kk].second == pr.j) {
                pair_idx = kk; break;
            }
        }
        const auto& pls = per_pair_polylines[pair_idx];
        for (std::size_t pi = 0; pi < pls.size(); ++pi) {
            mfest << "        [";
            for (std::size_t vi = 0; vi < pls[pi].size(); ++vi) {
                const auto& p = pls[pi][vi];
                mfest << "[";
                json_double(mfest, p.x()); mfest << ", ";
                json_double(mfest, p.y()); mfest << ", ";
                json_double(mfest, p.z()); mfest << "]";
                if (vi + 1 < pls[pi].size()) mfest << ", ";
            }
            mfest << "]" << (pi + 1 < pls.size() ? "," : "") << "\n";
        }
        mfest << "      ]\n";
        mfest << "    }" << (k + 1 < pair_results.size() ? "," : "") << "\n";
    }
    mfest << "  ]\n";
    mfest << "}\n";
    mfest.close();

    if (a.verbose) {
        std::cerr << "wrote " << meshes.size() << " mesh(es) and manifest "
                  << a.manifest_path << "\n";
    }

    // Check for any non-conformal pair.
    bool all_conformal = true;
    for (const auto& pr : pair_results) {
        if (!pr.conformal) { all_conformal = false; break; }
    }
    return all_conformal ? 0 : 7;
}
