// SAFS corefine_faults — I/O implementations.

#include "io.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/helpers.h>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace safs::io {

bool read_ascii_stl(const std::filesystem::path& p,
                    Mesh& out_mesh,
                    std::string* out_solid_name) {
    out_mesh.clear();

    if (!std::filesystem::exists(p)) {
        std::cerr << "[corefine_faults] STL not found: " << p << '\n';
        return false;
    }

    if (out_solid_name) {
        // Best-effort: parse the first `solid <name>` line for the
        // human-readable name.  We do NOT use CGAL's solid-name reader
        // because the polygon-soup recipe doesn't expose it.
        std::ifstream fh(p);
        std::string token;
        if (fh >> token && token == "solid") {
            std::string rest;
            std::getline(fh, rest);
            // strip leading whitespace
            const auto first = rest.find_first_not_of(" \t\r\n");
            if (first != std::string::npos) {
                rest.erase(0, first);
            }
            // strip trailing whitespace
            const auto last = rest.find_last_not_of(" \t\r\n");
            if (last != std::string::npos) {
                rest.erase(last + 1);
            }
            *out_solid_name = rest;
        } else {
            out_solid_name->clear();
        }
    }

    std::vector<Kernel::Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;
    if (!CGAL::IO::read_polygon_soup(p.string(), points, polygons)) {
        std::cerr << "[corefine_faults] " << p
                  << ": read_polygon_soup failed.\n";
        return false;
    }

    // ASCII STL stores three vertex slots per triangle so the soup is
    // non-manifold by construction; repair_polygon_soup dedups shared
    // vertices and drops degenerate triangles.
    //
    // R-010: require_same_orientation(true) — two triangles with
    // identical vertex sets but opposite winding (e.g., a buggy
    // ts_to_stl.py emitting both sides of a closed surface) must NOT
    // be silently merged into a non-orientable mesh.  Fault patches
    // are open with a single defined normal, so this never fires on
    // well-formed inputs; it is a correctness guard against bad data.
    PMP::repair_polygon_soup(
        points, polygons,
        CGAL::parameters::erase_all_duplicates(true)
                         .require_same_orientation(true));
    PMP::orient_polygon_soup(points, polygons);

    if (!PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
        std::cerr << "[corefine_faults] " << p
                  << ": polygon soup is non-manifold even after "
                  << "repair; ts_to_stl.py output is corrupt.\n";
        return false;
    }

    PMP::polygon_soup_to_polygon_mesh(points, polygons, out_mesh);
    return CGAL::is_valid_polygon_mesh(out_mesh);
}

bool write_ascii_stl(const std::filesystem::path& p,
                     const Mesh& mesh,
                     std::string_view solid_name) {
    if (mesh.number_of_faces() == 0) {
        std::cerr << "[corefine_faults] refusing to write empty STL: "
                  << p << '\n';
        return false;
    }

    std::error_code ec;
    std::filesystem::create_directories(p.parent_path(), ec);

    std::FILE* fh = std::fopen(p.string().c_str(), "w");
    if (!fh) {
        std::cerr << "[corefine_faults] cannot open " << p
                  << " for writing.\n";
        return false;
    }

    std::string name(solid_name);
    std::fprintf(fh, "solid SAFS:%s:conformal\n", name.c_str());

    for (auto f : mesh.faces()) {
        // Walk the three half-edges of the triangle.
        auto h0 = mesh.halfedge(f);
        auto h1 = mesh.next(h0);
        auto h2 = mesh.next(h1);
        const auto& a = mesh.point(mesh.target(h0));
        const auto& b = mesh.point(mesh.target(h1));
        const auto& c = mesh.point(mesh.target(h2));

        // Mesh halfedge convention: target(h0) = first vertex,
        // target(h1) = second, target(h2) = third.  But we want the
        // outgoing convention so vertices wind counter-clockwise.
        // For Surface_mesh, faces are CCW by default and walking
        // next() from any halfedge visits the three targets in CCW
        // order, so {a, b, c} as above is correct.

        const double ax = a.x(), ay = a.y(), az = a.z();
        const double bx = b.x(), by = b.y(), bz = b.z();
        const double cx = c.x(), cy = c.y(), cz = c.z();

        // Robust outward normal from cross(b-a, c-a).
        const double ux = bx - ax, uy = by - ay, uz = bz - az;
        const double vx = cx - ax, vy = cy - ay, vz = cz - az;
        double nx = uy * vz - uz * vy;
        double ny = uz * vx - ux * vz;
        double nz = ux * vy - uy * vx;
        const double nrm = std::sqrt(nx*nx + ny*ny + nz*nz);
        if (nrm > 0.0) {
            nx /= nrm; ny /= nrm; nz /= nrm;
        } // else leave (0,0,0) — same as the Python implementation.

        std::fprintf(fh,
            "  facet normal %+.17e %+.17e %+.17e\n"
            "    outer loop\n"
            "      vertex %+.17e %+.17e %+.17e\n"
            "      vertex %+.17e %+.17e %+.17e\n"
            "      vertex %+.17e %+.17e %+.17e\n"
            "    endloop\n  endfacet\n",
            nx, ny, nz,
            ax, ay, az,
            bx, by, bz,
            cx, cy, cz);
    }

    std::fprintf(fh, "endsolid SAFS:%s:conformal\n", name.c_str());
    std::fclose(fh);
    return true;
}

namespace {

// Format a double with %.17g (round-trippable).  Hand-rolled so we
// don't depend on std::to_chars on older toolchains.
std::string fmt_g(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

// Escape a JSON string.  Only handles the characters that can appear
// in fault short names and filesystem paths on macOS / Linux.
std::string json_escape(std::string_view s) {
    std::string out;
    out.reserve(s.size() + 2);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b";  break;
            case '\f': out += "\\f";  break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x",
                                  static_cast<unsigned>(c));
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    return out;
}

void write_or_die(const std::filesystem::path& p,
                  std::string_view content) {
    std::error_code ec;
    std::filesystem::create_directories(p.parent_path(), ec);
    std::ofstream fh(p);
    if (!fh) {
        throw std::runtime_error("[corefine_faults] cannot open " +
                                 p.string() + " for writing");
    }
    fh << content;
}

} // namespace

void write_triangle_to_fault_json(
    const std::filesystem::path& p,
    const std::vector<PerFaultStats>& per_fault) {

    std::ostringstream os;
    os << "{\n";
    os << "  \"schema_version\": 1,\n";

    os << "  \"faults\": {\n";
    std::int64_t cursor = 0;
    for (std::size_t i = 0; i < per_fault.size(); ++i) {
        const auto& f = per_fault[i];
        const std::int64_t lo = cursor;
        const std::int64_t hi = cursor + f.n_triangles;
        cursor = hi;

        os << "    \"" << json_escape(f.short_name) << "\": {"
           << "\"n_triangles\": " << f.n_triangles
           << ", \"range\": [" << lo << ", " << hi << "]}";
        if (i + 1 < per_fault.size()) os << ',';
        os << '\n';
    }
    os << "  },\n";
    os << "  \"n_total_triangles\": " << cursor << '\n';
    os << "}\n";

    write_or_die(p, os.str());
}

void write_intersection_report_json(
    const std::filesystem::path& p,
    const std::vector<PerFaultStats>& per_fault,
    const std::vector<PairStats>& pairs,
    double snap_m, double clearance_m, double target_edge_length_m) {

    std::ostringstream os;
    os << "{\n";
    os << "  \"schema_version\": 1,\n";
    os << "  \"backend\": \"cgal\",\n";
    os << "  \"target_edge_length_m\": " << fmt_g(target_edge_length_m) << ",\n";
    os << "  \"snap_m\": "                << fmt_g(snap_m)               << ",\n";
    os << "  \"clearance_m\": "           << fmt_g(clearance_m)          << ",\n";

    os << "  \"pairs\": {\n";
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto& pr = pairs[i];
        const std::string key = pr.short_a + "__x__" + pr.short_b;
        os << "    \"" << json_escape(key) << "\": {\n";
        os << "      \"n_constrained_edges_A\": " << pr.n_constrained_edges_A << ",\n";
        os << "      \"n_constrained_edges_B\": " << pr.n_constrained_edges_B << ",\n";
        os << "      \"pre_split_n_tri_A\": "     << pr.pre_n_tri_A           << ",\n";
        os << "      \"post_split_n_tri_A\": "    << pr.post_n_tri_A          << ",\n";
        os << "      \"pre_split_n_tri_B\": "     << pr.pre_n_tri_B           << ",\n";
        os << "      \"post_split_n_tri_B\": "    << pr.post_n_tri_B          << ",\n";
        os << "      \"remesh_iters\": "          << pr.remesh_iters          << ",\n";
        os << "      \"gates\": {\n";
        os << "        \"manifold_A\": \""                << json_escape(pr.gate_manifold_A)           << "\",\n";
        os << "        \"manifold_B\": \""                << json_escape(pr.gate_manifold_B)           << "\",\n";
        os << "        \"polyline_edge_coincidence\": \"" << json_escape(pr.gate_polyline_coincidence) << "\",\n";
        os << "        \"interior_crossing_only\": \""    << json_escape(pr.gate_interior_only)        << "\"\n";
        os << "      }\n";
        os << "    }";
        if (i + 1 < pairs.size()) os << ',';
        os << '\n';
    }
    os << "  },\n";

    os << "  \"per_fault\": {\n";
    for (std::size_t i = 0; i < per_fault.size(); ++i) {
        const auto& f = per_fault[i];
        // R-007: emit just the filename (the JSON lives next to the
        // STLs, so a filename-only path remains correct after the
        // output directory is moved).
        os << "    \"" << json_escape(f.short_name) << "\": {"
           << "\"n_vertices\": "  << f.n_vertices
           << ", \"n_triangles\": " << f.n_triangles
           << ", \"stl_path\": \""
           << json_escape(
                  std::filesystem::path(f.stl_path).filename().string())
           << "\"}";
        if (i + 1 < per_fault.size()) os << ',';
        os << '\n';
    }
    os << "  }\n";
    os << "}\n";

    write_or_die(p, os.str());
}

int validate_triangle_to_fault_json(
    const std::filesystem::path& p,
    const std::vector<PerFaultStats>& per_fault) {
    std::ifstream fh(p);
    if (!fh) {
        std::cerr << "[corefine_faults] schema: cannot re-open " << p
                  << " for validation.\n";
        return 3;
    }
    const std::string body((std::istreambuf_iterator<char>(fh)),
                            std::istreambuf_iterator<char>());

    // (a) Pull n_total_triangles and verify it matches the in-memory sum.
    const auto pos = body.find("\"n_total_triangles\":");
    if (pos == std::string::npos) {
        std::cerr << "[corefine_faults] schema: n_total_triangles "
                     "key missing from " << p << '\n';
        return 3;
    }
    long long on_disk_total_ll = 0;
    if (std::sscanf(body.c_str() + pos,
                    "\"n_total_triangles\": %lld",
                    &on_disk_total_ll) != 1) {
        std::cerr << "[corefine_faults] schema: n_total_triangles "
                     "parse failed in " << p << '\n';
        return 3;
    }
    const std::int64_t on_disk_total = on_disk_total_ll;

    std::int64_t in_mem_sum = 0;
    for (const auto& f : per_fault) in_mem_sum += f.n_triangles;
    if (on_disk_total != in_mem_sum) {
        std::cerr << "[corefine_faults] schema: on-disk total "
                  << on_disk_total << " != in-memory sum "
                  << in_mem_sum << '\n';
        return 3;
    }

    // (b) Per-fault ranges must form a contiguous partition.
    std::int64_t cursor = 0;
    for (const auto& f : per_fault) {
        const std::string fkey = "\"" + f.short_name + "\":";
        const auto fpos = body.find(fkey);
        if (fpos == std::string::npos) {
            std::cerr << "[corefine_faults] schema: fault "
                      << f.short_name << " missing from " << p << '\n';
            return 3;
        }
        long long n_tri_ll = 0, lo_ll = -1, hi_ll = -1;
        if (std::sscanf(body.c_str() + fpos,
                "\"%*[^\"]\": {\"n_triangles\": %lld, \"range\": [%lld, %lld]}",
                &n_tri_ll, &lo_ll, &hi_ll) != 3) {
            std::cerr << "[corefine_faults] schema: fault "
                      << f.short_name << " entry parse failed in "
                      << p << '\n';
            return 3;
        }
        if (n_tri_ll != f.n_triangles) {
            std::cerr << "[corefine_faults] schema: fault "
                      << f.short_name << " n_triangles on disk "
                      << n_tri_ll << " != in-memory " << f.n_triangles
                      << '\n';
            return 3;
        }
        if (lo_ll != cursor || hi_ll != cursor + f.n_triangles) {
            std::cerr << "[corefine_faults] schema: fault "
                      << f.short_name << " range [" << lo_ll << ", "
                      << hi_ll << "] does not match expected ["
                      << cursor << ", " << cursor + f.n_triangles
                      << "] in " << p << '\n';
            return 3;
        }
        cursor = hi_ll;
    }
    if (cursor != on_disk_total) {
        std::cerr << "[corefine_faults] schema: cumulative "
                  << cursor << " != on-disk total "
                  << on_disk_total << " in " << p << '\n';
        return 3;
    }
    return 0;
}

} // namespace safs::io
