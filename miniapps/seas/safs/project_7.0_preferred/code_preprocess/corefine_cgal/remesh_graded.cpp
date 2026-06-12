// remesh_graded — graded isotropic remeshing of ONE surface with the rim
// (border) byte-identical pre/post (Phase 3 of
// PLAN_mesh_quality_safv4_remesh_2026-06-12.md).
//
//   remesh_graded IN OUT --trace-points FILE.xyz
//       [--h-near 500] [--h-far 2500] [--d-near 1000] [--d-far 9000]
//       [--iterations 5] [--no-protect-border] [--verbose]
//
// Target edge length h(d) = h_near + (h_far - h_near) *
// clamp((d - d_near)/(d_far - d_near), 0, 1), d = distance to the nearest
// trace point (one "x y z" per line in FILE.xyz; emitted by
// extend_fault_borders.py from the Phase 2 manifest trace chains).
//
// Border edges are constrained (edge_is_constrained_map) and border
// vertices fixed (vertex_is_constrained_map) with protect_constraints =
// true; the tool verifies the border vertex coordinate set is identical
// pre/post and exits 4 if not.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "io_helpers.h"
#include "graded_sizing_field.h"

using K     = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Mesh  = CGAL::Surface_mesh<Point>;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace pp  = CGAL::parameters;

struct Args {
    std::string in_path, out_path, trace_points;
    double h_near = 500.0, h_far = 2500.0, d_near = 1000.0, d_far = 9000.0;
    int iterations = 5;
    bool protect_border = true;
    bool verbose = false;
};

static int usage(const char* argv0, int code) {
    std::cerr << "usage: " << argv0 << " IN OUT --trace-points FILE.xyz"
              << " [--h-near M] [--h-far M] [--d-near M] [--d-far M]"
              << " [--iterations N] [--no-protect-border] [--verbose]\n";
    return code;
}

int main(int argc, char** argv) {
    Args a;
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        auto need = [&](double& x) {
            if (++i >= argc) { std::cerr << "missing value for " << s << "\n"; std::exit(2); }
            x = std::strtod(argv[i], nullptr);
        };
        if (s == "-h" || s == "--help") return usage(argv[0], 0);
        else if (s == "--h-near") need(a.h_near);
        else if (s == "--h-far")  need(a.h_far);
        else if (s == "--d-near") need(a.d_near);
        else if (s == "--d-far")  need(a.d_far);
        else if (s == "--iterations") {
            if (++i >= argc) { std::cerr << "missing value\n"; return 2; }
            a.iterations = std::atoi(argv[i]);
        }
        else if (s == "--trace-points") {
            if (++i >= argc) { std::cerr << "missing value\n"; return 2; }
            a.trace_points = argv[i];
        }
        else if (s == "--no-protect-border") a.protect_border = false;
        else if (s == "--verbose") a.verbose = true;
        else pos.push_back(s);
    }
    if (pos.size() != 2 || a.trace_points.empty()) return usage(argv[0], 2);
    a.in_path = pos[0];
    a.out_path = pos[1];

    Mesh m;
    if (!io_helpers::read_polygon_mesh_any(a.in_path, m)) {
        std::cerr << "cannot read " << a.in_path << "\n"; return 3;
    }

    std::vector<Point> sites;
    {
        std::ifstream f(a.trace_points);
        if (!f) { std::cerr << "cannot read " << a.trace_points << "\n"; return 3; }
        double x, y, z;
        while (f >> x >> y >> z) sites.emplace_back(x, y, z);
    }
    if (sites.empty()) {
        std::cerr << "no trace points in " << a.trace_points << "\n"; return 3;
    }
    if (a.verbose) {
        std::cerr << "in: V=" << m.number_of_vertices()
                  << " F=" << m.number_of_faces()
                  << "  sites=" << sites.size()
                  << "  h=[" << a.h_near << ", " << a.h_far << "] m"
                  << "  d=[" << a.d_near << ", " << a.d_far << "] m\n";
    }

    // Constrain the rim: border edges + border vertices.
    using EdgeDesc   = boost::graph_traits<Mesh>::edge_descriptor;
    using VertexDesc = boost::graph_traits<Mesh>::vertex_descriptor;
    Mesh::Property_map<EdgeDesc, bool> ecm;
    Mesh::Property_map<VertexDesc, bool> vcm;
    std::tie(ecm, std::ignore) =
        m.add_property_map<EdgeDesc, bool>("e:border_constrained", false);
    std::tie(vcm, std::ignore) =
        m.add_property_map<VertexDesc, bool>("v:border_constrained", false);

    std::set<std::tuple<double, double, double>> rim_before;
    std::size_t n_border_e = 0;
    for (auto h : m.halfedges()) {
        if (!m.is_border(h)) continue;
        put(ecm, m.edge(h), true);
        put(vcm, m.target(h), true);
        ++n_border_e;
        const auto& p = m.point(m.target(h));
        rim_before.insert({p.x(), p.y(), p.z()});
    }
    if (a.verbose) {
        std::cerr << "border: " << n_border_e << " halfedges, "
                  << rim_before.size() << " rim vertices\n";
    }

    Graded_polyline_sizing_field<Mesh> field(a.h_near, a.h_far,
                                             a.d_near, a.d_far, sites, m);
    if (a.protect_border) {
        PMP::isotropic_remeshing(faces(m), field, m,
            pp::edge_is_constrained_map(ecm)
              .vertex_is_constrained_map(vcm)
              .protect_constraints(true)
              .number_of_iterations(a.iterations));
    } else {
        PMP::isotropic_remeshing(faces(m), field, m,
            pp::number_of_iterations(a.iterations));
    }

    // Verify the rim vertex coordinate set is unchanged.
    std::set<std::tuple<double, double, double>> rim_after;
    for (auto h : m.halfedges()) {
        if (!m.is_border(h)) continue;
        const auto& p = m.point(m.target(h));
        rim_after.insert({p.x(), p.y(), p.z()});
    }
    if (a.protect_border && rim_after != rim_before) {
        std::cerr << "ERROR: rim vertex set changed ("
                  << rim_before.size() << " -> " << rim_after.size()
                  << " distinct coords). Aborting.\n";
        return 4;
    }

    if (!io_helpers::write_polygon_mesh_high_precision(a.out_path, m)) {
        std::cerr << "cannot write " << a.out_path << "\n"; return 3;
    }
    if (a.verbose) {
        std::cerr << "out: V=" << m.number_of_vertices()
                  << " F=" << m.number_of_faces()
                  << "  rim byte-identical: "
                  << (a.protect_border ? "yes" : "n/a") << "\n";
    }
    std::cout << "remesh_graded: " << a.in_path << " -> " << a.out_path
              << "  V=" << m.number_of_vertices()
              << " F=" << m.number_of_faces() << "\n";
    return 0;
}
