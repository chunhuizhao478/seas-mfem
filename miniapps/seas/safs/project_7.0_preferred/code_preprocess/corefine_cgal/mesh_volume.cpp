// mesh_volume — Build a 3D tet mesh of a bounding-box volume with N corefined
// fault surfaces embedded as INTERNAL constraints, using CGAL Mesh_3's
// Polyhedral_complex_mesh_domain_3.  Replaces the gmsh `Surface{} In Volume{}`
// path in PLAN_cgal_corefine_multifault.md §Phase 3 because tetgen
// (Algorithm3D=1) rejects multi-fault PLCs at shared polyline edges.
//
// CGAL's Polyhedral_complex_mesh_domain_3 is designed for exactly this
// scenario: "two polyhedral surfaces of the complex are either disjoint or
// share an intersection that is itself a polyhedral surface in the complex".
// Our corefine pipeline guarantees the precondition.
//
// Usage:
//   mesh_volume MANIFEST OUTPUT_MEDIT_MESH
//     [--lc-near 1500] [--lc-min 100] [--lc-far 15000]
//     [--pad-xy 50000] [--pad-top 100] [--pad-bottom 25000]
//     [--dist-inner 3000] [--dist-outer 40000]
//     [--verbose]
//
// Output: MEDIT .mesh format.  Use CGAL or meshio to convert to gmsh .msh
// or VTK .vtu downstream.
//
// Subdomain layout:
//   - subdomain 0 = outside of bounding box (not meshed)
//   - subdomain k+1 = "rock" sub-region created by fault k (k=0..N-1)
// All sub-regions share the same physical interpretation ("rock"); the
// post-processor relabels them to a single attribute for MFEM consumption.
// Splitting by fault gives Mesh_3 a clear topological story and lets the
// downstream pipeline read per-fault triangle attributes for friction BCs.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Polyhedral_complex_mesh_domain_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/STL.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include "distance_sizing_field.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace fs   = std::filesystem;
namespace PMP  = CGAL::Polygon_mesh_processing;

using K          = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Mesh_polyhedron_3<K>::type;
// Polyhedral_complex_mesh_domain_3<IGT, Polyhedron> — only 2 template args.
using MeshDomain = CGAL::Polyhedral_complex_mesh_domain_3<K, Polyhedron>;

// CGAL_MESH_3_USE_RELAXED_HEURISTICS / Concurrent_tag: keep simple sequential.
// Pass K explicitly because Polyhedral_complex_mesh_domain_3 doesn't expose
// the nested ::R typedef that Kernel_traits<> defaults to.
using Tr   = CGAL::Mesh_triangulation_3<MeshDomain,
                                        K,
                                        CGAL::Sequential_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr,
                                MeshDomain::Corner_index,
                                MeshDomain::Curve_index>;
using Criteria = CGAL::Mesh_criteria_3<Tr>;

namespace pp = CGAL::parameters;

// -----------------------------------------------------------------------------
// CLI args.
struct Args {
    std::string manifest_path;
    std::string output_path;
    double lc_near    = 1500.0;
    double lc_min     =  100.0;
    double lc_far     = 15000.0;
    double pad_xy     = 50000.0;
    double pad_top    =   100.0;
    double pad_bottom = 25000.0;
    double dist_inner =  3000.0;
    double dist_outer = 40000.0;
    double feature_angle_deg = 60.0;  // for detect_features
    bool   verbose    = false;
};

static int usage(const char* argv0, int code) {
    std::cerr << "usage: " << argv0
              << " MANIFEST OUTPUT_MEDIT_MESH"
              << " [--lc-near SIZE] [--lc-min SIZE] [--lc-far SIZE]"
              << " [--pad-xy M] [--pad-top M] [--pad-bottom M]"
              << " [--dist-inner M] [--dist-outer M]"
              << " [--feature-angle DEG] [--verbose]\n";
    return code;
}

static int parse_args(int argc, char** argv, Args& a) {
    auto need_double = [&](int& i, double& x, const std::string& name) {
        if (++i >= argc) { std::cerr << "missing value for " << name << "\n"; std::exit(2); }
        char* end = nullptr;
        x = std::strtod(argv[i], &end);
        if (end == argv[i]) {
            std::cerr << "invalid number for " << name << ": " << argv[i] << "\n";
            std::exit(2);
        }
    };
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "-h" || s == "--help") return usage(argv[0], 0);
        else if (s == "--lc-near")     need_double(i, a.lc_near, s);
        else if (s == "--lc-min")      need_double(i, a.lc_min, s);
        else if (s == "--lc-far")      need_double(i, a.lc_far, s);
        else if (s == "--pad-xy")      need_double(i, a.pad_xy, s);
        else if (s == "--pad-top")     need_double(i, a.pad_top, s);
        else if (s == "--pad-bottom")  need_double(i, a.pad_bottom, s);
        else if (s == "--dist-inner")  need_double(i, a.dist_inner, s);
        else if (s == "--dist-outer")  need_double(i, a.dist_outer, s);
        else if (s == "--feature-angle") need_double(i, a.feature_angle_deg, s);
        else if (s == "--verbose")     a.verbose = true;
        else pos.push_back(s);
    }
    if (pos.size() != 2) return usage(argv[0], 2);
    a.manifest_path = pos[0];
    a.output_path   = pos[1];
    return -1;
}

// -----------------------------------------------------------------------------
// Trivial JSON parser (we only need a flat key->value lookup + meshes[].basename
// + bbox).  Avoids pulling in a JSON dep.
struct Manifest {
    double mesh_edge_size  = 1500.0;
    double min_edge        =  100.0;
    std::vector<std::string> basenames;
    double xlo, xhi, ylo, yhi, zlo, zhi;
};

static bool find_double(const std::string& s, const std::string& key, double& out) {
    auto p = s.find("\"" + key + "\"");
    if (p == std::string::npos) return false;
    p = s.find(':', p);
    if (p == std::string::npos) return false;
    out = std::strtod(s.c_str() + p + 1, nullptr);
    return true;
}

static Manifest load_manifest(const std::string& path,
                              const fs::path& corefined_dir) {
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open manifest: " << path << "\n"; std::exit(1); }
    std::string content((std::istreambuf_iterator<char>(f)),
                        std::istreambuf_iterator<char>());
    Manifest m;
    if (!find_double(content, "mesh_edge_size", m.mesh_edge_size))
        std::cerr << "warning: manifest missing mesh_edge_size; using default 1500\n";
    if (!find_double(content, "min_edge", m.min_edge))
        std::cerr << "warning: manifest missing min_edge; using default 100\n";

    // Walk meshes[].basename and meshes[].bbox.
    bool first_bb = true;
    std::size_t p = 0;
    while (true) {
        auto bp = content.find("\"basename\"", p);
        if (bp == std::string::npos) break;
        bp = content.find('"', bp + 11);                  // skip past "basename":
        if (bp == std::string::npos) break;
        auto bp_end = content.find('"', bp + 1);
        if (bp_end == std::string::npos) break;
        m.basenames.push_back(content.substr(bp + 1, bp_end - bp - 1));

        // Find the next bbox after this basename.
        auto bbox_p = content.find("\"bbox\"", bp_end);
        if (bbox_p == std::string::npos) break;
        double xlo, xhi, ylo, yhi, zlo, zhi;
        std::string sub = content.substr(bbox_p, 400);
        find_double(sub, "x_lo", xlo); find_double(sub, "x_hi", xhi);
        find_double(sub, "y_lo", ylo); find_double(sub, "y_hi", yhi);
        find_double(sub, "z_lo", zlo); find_double(sub, "z_hi", zhi);
        if (first_bb) {
            m.xlo = xlo; m.xhi = xhi; m.ylo = ylo; m.yhi = yhi;
            m.zlo = zlo; m.zhi = zhi; first_bb = false;
        } else {
            m.xlo = std::min(m.xlo, xlo); m.xhi = std::max(m.xhi, xhi);
            m.ylo = std::min(m.ylo, ylo); m.yhi = std::max(m.yhi, yhi);
            m.zlo = std::min(m.zlo, zlo); m.zhi = std::max(m.zhi, zhi);
        }
        p = bbox_p + 1;
    }
    if (m.basenames.empty()) {
        std::cerr << "manifest contained no meshes[]\n";
        std::exit(1);
    }
    (void)corefined_dir;  // path-resolution done by caller
    return m;
}

// -----------------------------------------------------------------------------
// Build the closed bounding-box polyhedron with outward-oriented faces.
// 8 corners, 12 triangular faces.  Used as the OUTER boundary in the
// polyhedral complex.
static void build_box_polyhedron(double xmin, double xmax,
                                 double ymin, double ymax,
                                 double zmin, double zmax,
                                 Polyhedron& out) {
    using P = K::Point_3;
    std::vector<P> pts = {
        P(xmin, ymin, zmin), P(xmax, ymin, zmin),
        P(xmax, ymax, zmin), P(xmin, ymax, zmin),
        P(xmin, ymin, zmax), P(xmax, ymin, zmax),
        P(xmax, ymax, zmax), P(xmin, ymax, zmax),
    };
    // Triangulate each face with outward-pointing normals.
    // Vertex layout: 0..3 bottom (z=zmin), 4..7 top (z=zmax).
    std::vector<std::vector<std::size_t>> tris = {
        // bottom (normal -z)
        {0, 3, 2}, {0, 2, 1},
        // top (normal +z)
        {4, 5, 6}, {4, 6, 7},
        // front (y=ymin, normal -y)
        {0, 1, 5}, {0, 5, 4},
        // back (y=ymax, normal +y)
        {2, 3, 7}, {2, 7, 6},
        // left (x=xmin, normal -x)
        {3, 0, 4}, {3, 4, 7},
        // right (x=xmax, normal +x)
        {1, 2, 6}, {1, 6, 5},
    };
    PMP::orient_polygon_soup(pts, tris);
    PMP::polygon_soup_to_polygon_mesh(pts, tris, out);
    if (!CGAL::is_closed(out)) {
        std::cerr << "internal: box polyhedron is not closed\n"; std::exit(1);
    }
    if (PMP::is_outward_oriented(out) == false) {
        PMP::reverse_face_orientations(out);
    }
}

// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
    Args a;
    if (int r = parse_args(argc, argv, a); r >= 0) return r;

    if (a.lc_min <= 0.0) {
        std::cerr << "error: lc_min must be > 0 (got " << a.lc_min << ")\n";
        return 2;
    }
    if (a.lc_min > a.lc_near) {
        std::cerr << "error: lc_min (" << a.lc_min
                  << ") > lc_near (" << a.lc_near
                  << "); the hard floor cannot exceed the target size\n";
        return 2;
    }
    if (a.lc_near > a.lc_far) {
        std::cerr << "error: lc_near (" << a.lc_near
                  << ") > lc_far (" << a.lc_far << ")\n";
        return 2;
    }

    fs::path manifest_p(a.manifest_path);
    fs::path corefined_dir = manifest_p.parent_path();
    Manifest m = load_manifest(a.manifest_path, corefined_dir);

    if (a.verbose) {
        std::cerr << "manifest: " << a.manifest_path << "\n"
                  << "  faults: " << m.basenames.size() << "\n"
                  << "  mesh_edge_size = " << m.mesh_edge_size << "\n"
                  << "  min_edge       = " << m.min_edge << "\n"
                  << "  union bbox: x=[" << m.xlo << ", " << m.xhi << "]"
                  << "  y=[" << m.ylo << ", " << m.yhi << "]"
                  << "  z=[" << m.zlo << ", " << m.zhi << "]\n";
    }

    // Build the bounding-box polyhedron with paddings.
    const double xmin = m.xlo - a.pad_xy;
    const double xmax = m.xhi + a.pad_xy;
    const double ymin = m.ylo - a.pad_xy;
    const double ymax = m.yhi + a.pad_xy;
    const double zmin = m.zlo - a.pad_bottom;
    const double zmax = m.zhi + a.pad_top;
    Polyhedron box;
    build_box_polyhedron(xmin, xmax, ymin, ymax, zmin, zmax, box);
    if (a.verbose) {
        std::cerr << "box polyhedron V=" << num_vertices(box)
                  << " F=" << num_faces(box) << "\n";
    }

    // Read each fault STL as a Polyhedron.
    std::vector<Polyhedron> fault_polys(m.basenames.size());
    for (std::size_t i = 0; i < m.basenames.size(); ++i) {
        fs::path stl_p = corefined_dir /
            (m.basenames[i] + std::string("_corefined.stl"));
        std::vector<K::Point_3> pts;
        std::vector<std::vector<std::size_t>> polys;
        if (!CGAL::IO::read_STL(stl_p.string(), pts, polys)) {
            std::cerr << "cannot read STL: " << stl_p << "\n"; return 3;
        }
        PMP::orient_polygon_soup(pts, polys);
        PMP::polygon_soup_to_polygon_mesh(pts, polys, fault_polys[i]);
        if (a.verbose) {
            std::cerr << "fault " << i << ": " << m.basenames[i]
                      << "  V=" << num_vertices(fault_polys[i])
                      << "  F=" << num_faces(fault_polys[i]) << "\n";
        }
    }

    // Build the polyhedron list and the (positive,negative) subdomain pair list.
    // Subdomain layout:
    //   - 0 = outside (default; outward face of box maps to (1, 0))
    //   - 1 = rock-side-A
    //   - 2 = rock-side-B
    // Each fault MUST separate two DISTINCT subdomains.  Earlier we used pair
    // (1, 1) for "embedded constraint, both sides rock" — that turned out to
    // create cavities (12.6% of box volume left unmeshed; 540k missing tets;
    // 716k extra interior boundary faces) because Polyhedral_complex_mesh_
    // domain_3 cannot resolve adjacent-cell labels around a (i, i) facet and
    // assigns them to subdomain 0.  Using (1, 2) puts cells on either side of
    // every fault into a distinct subdomain; downstream relabels {1, 2} → rock.
    // The per-fault triangle distinction is preserved via surface_patch_index
    // (written to MEDIT output via show_patches=true below).
    std::vector<Polyhedron> all_polys;
    all_polys.reserve(1 + fault_polys.size());
    all_polys.push_back(std::move(box));
    for (auto& p : fault_polys) all_polys.push_back(std::move(p));

    std::vector<std::pair<MeshDomain::Subdomain_index,
                          MeshDomain::Subdomain_index>> indices_pairs;
    indices_pairs.reserve(all_polys.size());
    // Box: positive side (outward) = subdomain 0 (outside); negative side (inward) = 1 (rock-A).
    indices_pairs.emplace_back(0, 1);
    // Each fault: positive side = 1 (rock-A), negative side = 2 (rock-B).
    for (std::size_t i = 0; i < fault_polys.size(); ++i) {
        indices_pairs.emplace_back(1, 2);
    }

    MeshDomain domain(all_polys.begin(), all_polys.end(),
                      indices_pairs.begin(), indices_pairs.end());
    // detect_features and detect_borders both produce protection balls
    // around feature/border edges that force sub-100m local refinement,
    // creating sliver tets concentrated within 1km of the fault surfaces.
    // Skipping both: rely on patch-id mechanism alone for surface
    // preservation.  Box corners may smooth slightly (acceptable since
    // PAD_XY=50km ≫ LC_NEAR=1500m).

    if (a.verbose) {
        std::cerr << "domain built: " << all_polys.size()
                  << " polyhedra, feature angle = "
                  << a.feature_angle_deg << " deg\n";
    }

    // Spatially-varying sizing field: lc_near near faults, ramping to
    // lc_far far from faults.  Replaces the previous uniform constant size
    // (which made the bulk uniformly fine at ~1500m everywhere).  See
    // distance_sizing_field.h header for the ramp definition.  The field
    // holds pointers to all_polys[1..N] (the faults; index 0 is the box);
    // these polyhedra outlive the field because they are held by `all_polys`
    // for the duration of make_mesh_3.
    std::vector<const Polyhedron*> fault_pointers;
    for (std::size_t i = 1; i < all_polys.size(); ++i) {
        fault_pointers.push_back(&all_polys[i]);
    }
    Distance_sizing_field<K, Polyhedron, MeshDomain::Index>
        sizing(fault_pointers,
               a.lc_min, a.lc_near, a.lc_far,
               a.dist_inner, a.dist_outer);
    if (a.verbose) {
        std::cerr << "sizing field: lc_min=" << a.lc_min
                  << " lc_near=" << a.lc_near
                  << " lc_far=" << a.lc_far
                  << " dist_inner=" << a.dist_inner
                  << " dist_outer=" << a.dist_outer << "\n";
    }

    // Mesh criteria.  facet_size and cell_size are now SPATIAL FIELDS (graded);
    // edge_size stays uniform (it bounds 1D feature edges only); facet_min_
    // size, cell_min_size, edge_min_size are hard lower bounds (from Phase 1
    // of PLAN_min_edge_enforcement.md).
    Criteria criteria(
        pp::edge_size(a.lc_near),
        pp::edge_min_size(a.lc_min),
        pp::facet_angle(25.0),
        pp::facet_size(sizing),                // graded: lc_near→lc_far
        pp::facet_min_size(a.lc_min),
        pp::facet_distance(a.lc_near * 0.1),
        pp::cell_radius_edge_ratio(3.0),
        pp::cell_size(sizing),                 // graded: lc_near→lc_far
        pp::cell_min_size(a.lc_min)
    );

    if (a.verbose) {
        std::cerr << "calling make_mesh_3 (lc_near=" << a.lc_near
                  << " lc_min=" << a.lc_min << ")\n"
                  << "  hard floors: edge_min=" << a.lc_min
                  << " facet_min=" << a.lc_min
                  << " cell_min=" << a.lc_min << "\n"
                  << "  this can take several minutes for large domains\n";
    }
    // make_mesh_3 defaults: pp::no_lloyd, pp::no_odt, pp::perturb, pp::exude.
    // Tried adding pp::lloyd + pp::odt with time_limit=60 each — they
    // catastrophically degraded the mesh (median per-tet min-edge dropped
    // from 1560 m to 0.46 m).  Sticking with defaults.
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    // Write output (MEDIT .mesh format).
    std::ofstream out(a.output_path);
    if (!out) {
        std::cerr << "cannot open output: " << a.output_path << "\n"; return 4;
    }
    // rebind=true:        compact subdomain indices in the MEDIT output
    // show_patches=true:  write surface_patch_index per triangle as the MEDIT ref
    //                     (so downstream can identify per-fault triangles via
    //                      distinct patch indices; one ref = box surface).
    // all_cells=true:     CRITICAL — include subdomain-0 cells too.  With our
    //                     pair (1, 2) for each fault, CGAL Mesh_3 leaves a
    //                     fraction of cells in subdomain 0 due to ambiguous
    //                     labeling around our open fault patches (5.8% missing
    //                     coverage observed with all_cells=false).  Including
    //                     them here recovers full box volume coverage; the
    //                     downstream converter filters out cells whose centroid
    //                     falls outside the box (the Delaunay convex-hull
    //                     buffer).
    CGAL::IO::output_to_medit(out, c3t3,
                              /*rebind=*/true,
                              /*show_patches=*/true,
                              /*all_vertices=*/true,
                              /*all_cells=*/true);
    if (a.verbose) {
        std::cerr << "wrote MEDIT mesh: " << a.output_path << "\n"
                  << "  triangulation: V=" << c3t3.triangulation().number_of_vertices()
                  << "  C=" << c3t3.number_of_cells_in_complex()
                  << "  F=" << c3t3.number_of_facets_in_complex() << "\n";
    }
    return 0;
}
