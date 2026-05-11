#pragma once

// io_helpers.h — read/write helpers for the corefine pipeline.
//
// CGAL's default OFF/STL writers use operator<< on std::ostream with the
// default precision (~6 sig digits), which is insufficient for UTM
// coordinates ~10^6 m and silently destroys polyline conformality at every
// ASCII round-trip.  Both write_off_high_precision and
// write_stl_high_precision below force 15 significant digits (~10^-9 m at
// UTM Y) — six orders of magnitude tighter than gmsh's
// Geometry.Tolerance = 1e-3.

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace io_helpers {

// Lowercase the path's extension for dispatch.
inline std::string ext_lower(const std::string& path) {
    auto dot = path.find_last_of('.');
    if (dot == std::string::npos) return "";
    std::string e = path.substr(dot);
    std::transform(e.begin(), e.end(), e.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return e;
}

// Read OFF/STL/PLY/OBJ — auto-detected by extension.  Returns true on success.
template <class Mesh>
bool read_polygon_mesh_any(const std::string& path, Mesh& m) {
    return CGAL::IO::read_polygon_mesh(path, m);
}

// Write OFF at 15 significant digits.  CGAL's default OFF writer uses
// `<<` at 6 digits, which is unusable for UTM coordinates.  We bypass it
// and emit OFF directly.
template <class Mesh>
bool write_off_high_precision(const std::string& path, const Mesh& m) {
    using Point = typename Mesh::Point;
    std::ofstream out(path);
    if (!out) return false;
    out << std::setprecision(15);
    out << "OFF\n";
    out << m.number_of_vertices() << " " << m.number_of_faces() << " 0\n";
    // Map vertex_descriptor -> dense 0..N-1 index (Surface_mesh's
    // vertex_index() returns the descriptor id which is dense for non-removed
    // verts; we still go through num_vertices() for clarity).
    std::vector<typename Mesh::Vertex_index> vlist;
    vlist.reserve(m.number_of_vertices());
    for (auto v : m.vertices()) vlist.push_back(v);
    std::vector<std::size_t> idx_of(m.num_vertices(), 0);
    for (std::size_t i = 0; i < vlist.size(); ++i) {
        idx_of[vlist[i].idx()] = i;
        const Point& p = m.point(vlist[i]);
        out << static_cast<double>(p.x()) << " "
            << static_cast<double>(p.y()) << " "
            << static_cast<double>(p.z()) << "\n";
    }
    for (auto f : m.faces()) {
        std::vector<std::size_t> ids;
        for (auto v : CGAL::vertices_around_face(m.halfedge(f), m)) {
            ids.push_back(idx_of[v.idx()]);
        }
        out << ids.size();
        for (auto i : ids) out << " " << i;
        out << "\n";
    }
    return out.good();
}

// Write ASCII STL at 15 significant digits.  STL has no shared-vertex
// notion; we emit each triangle as its three vertex coordinates.
template <class Mesh>
bool write_stl_high_precision(const std::string& path, const Mesh& m) {
    using Point = typename Mesh::Point;
    std::ofstream out(path);
    if (!out) return false;
    out << std::setprecision(15);
    out << "solid corefined\n";
    for (auto f : m.faces()) {
        std::vector<Point> pts;
        for (auto v : CGAL::vertices_around_face(m.halfedge(f), m)) {
            pts.push_back(m.point(v));
        }
        if (pts.size() != 3) continue;  // skip non-triangle faces (shouldn't occur)
        // Compute outward normal (right-hand rule); STL spec wants it.
        const double ux = pts[1].x() - pts[0].x();
        const double uy = pts[1].y() - pts[0].y();
        const double uz = pts[1].z() - pts[0].z();
        const double vx = pts[2].x() - pts[0].x();
        const double vy = pts[2].y() - pts[0].y();
        const double vz = pts[2].z() - pts[0].z();
        double nx = uy * vz - uz * vy;
        double ny = uz * vx - ux * vz;
        double nz = ux * vy - uy * vx;
        const double mag = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (mag > 0.0) {
            nx /= mag; ny /= mag; nz /= mag;
        }
        out << "facet normal " << nx << " " << ny << " " << nz << "\n";
        out << "  outer loop\n";
        for (const auto& p : pts) {
            out << "    vertex "
                << static_cast<double>(p.x()) << " "
                << static_cast<double>(p.y()) << " "
                << static_cast<double>(p.z()) << "\n";
        }
        out << "  endloop\n";
        out << "endfacet\n";
    }
    out << "endsolid corefined\n";
    return out.good();
}

// Dispatch by extension.
template <class Mesh>
bool write_polygon_mesh_high_precision(const std::string& path, const Mesh& m) {
    auto e = ext_lower(path);
    if (e == ".off") return write_off_high_precision(path, m);
    if (e == ".stl") return write_stl_high_precision(path, m);
    // Fallback: CGAL default writer for less precision-sensitive formats.
    return CGAL::IO::write_polygon_mesh(path, m);
}

}  // namespace io_helpers
