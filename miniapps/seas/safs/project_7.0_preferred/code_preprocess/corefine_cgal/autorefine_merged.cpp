// autorefine_merged — Merge 6 corefined fault STLs + a bounding-box surface
// into ONE polygon soup, run CGAL's autorefine_triangle_soup to insert
// Steiner points at any residual intersections (segment-facet crossings,
// near-duplicate triangles, or self-intersections that survived corefine
// + isotropic_remeshing), and output a clean merged STL + a per-triangle
// marker JSON for tetgen_mesh.py to consume.
//
// Why this tool exists:
//   - corefine_set produces per-fault corefined STLs that share polylines
//     conformally between adjacent fault pairs (manifest pairs[*].max_diff
//     ≈ 0).  But isotropic_remeshing in step 6 modifies each mesh
//     independently and can introduce small residual crossings AT
//     LOCATIONS NOT PROTECTED BY THE PAIR-WISE ECM (e.g., where a
//     remeshed vertex from mesh A falls into a triangle interior of
//     mesh B that was not on a shared polyline).  On the SAFS 2000 m
//     fixture this is 1 fully-duplicate triangle (Garnet × SBMT-SAF) +
//     8 segment-facet crossings (3 unique pairs).  tetgen rejects
//     such inputs as "Invalid PLC: self-intersection".
//   - autorefine_triangle_soup is the documented PLC-cleanup primitive
//     for exactly this situation.  It splits intersecting triangles at
//     the intersection edge/point so the output has no interior
//     self-intersections, ready to feed to tetgen.
//   - Per-triangle markers are propagated via a custom visitor so the
//     downstream tool can still tell which output triangle came from
//     which fault (via the integer marker).
//
// Usage:
//   autorefine_merged MANIFEST OUTPUT_STL OUTPUT_MARKERS_JSON
//     [--pad-xy 50000] [--pad-top 100] [--pad-bottom 25000]
//     [--box-marker 100] [--verbose]
//
// Markers (must match tetgen_mesh.py's expectation):
//   - 100 for box surface triangles
//   - 1..N for each fault (alphabetical order of basenames)
// where N = manifest.meshes.size().

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "io_helpers.h"

namespace fs   = std::filesystem;
namespace PMP  = CGAL::Polygon_mesh_processing;
namespace pp   = CGAL::parameters;

using K       = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point   = K::Point_3;

// -----------------------------------------------------------------------------
// Marker-propagating visitor.  CGAL's PMPAutorefinementVisitor concept:
//   number_of_output_triangles(n)     — total output count, called first
//   verbatim_triangle_copy(out, src)  — output triangle is an unmodified
//                                       copy of input triangle src
//   new_subtriangle(out, src)         — output triangle is a Steiner-
//                                       inserted subdivision of input src
//   delete_triangle(src)              — input triangle dropped (degenerate)
struct Marker_propagating_visitor {
    const std::vector<int>* in;   // input markers
    std::vector<int>*       out;  // output markers (sized in callback)
    int n_faults = 0;             // markers 1..n_faults are faults
    std::size_t n_fault_subtriangles = 0;   // fault-marked splits (gate)

    void number_of_output_triangles(std::size_t n) { out->assign(n, -1); }
    void verbatim_triangle_copy(std::size_t tgt, std::size_t src) {
        (*out)[tgt] = (*in)[src];
    }
    void new_subtriangle(std::size_t tgt, std::size_t src) {
        (*out)[tgt] = (*in)[src];
        const int m = (*in)[src];
        if (m >= 1 && m <= n_faults) ++n_fault_subtriangles;
    }
    void delete_triangle(std::size_t /*src*/) {}
};

// -----------------------------------------------------------------------------
// Trivial JSON parsing for the manifest (basename + bbox per mesh).  Mirrors
// mesh_volume.cpp's Manifest loader; kept duplicated so this tool is a
// standalone executable.
struct ManifestEntry { std::string basename; double xlo, xhi, ylo, yhi, zlo, zhi; };

static std::vector<ManifestEntry>
load_manifest_entries(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open manifest: " << path << "\n"; std::exit(1); }
    std::string content((std::istreambuf_iterator<char>(f)),
                        std::istreambuf_iterator<char>());
    std::vector<ManifestEntry> out;
    std::size_t p = 0;
    auto find_double_in = [&](const std::string& sub, const std::string& key,
                              double& v) {
        auto kp = sub.find("\"" + key + "\"");
        if (kp == std::string::npos) return false;
        kp = sub.find(':', kp);
        if (kp == std::string::npos) return false;
        v = std::strtod(sub.c_str() + kp + 1, nullptr);
        return true;
    };
    while (true) {
        auto bp = content.find("\"basename\"", p);
        if (bp == std::string::npos) break;
        bp = content.find('"', bp + 11);
        if (bp == std::string::npos) break;
        auto bp_end = content.find('"', bp + 1);
        if (bp_end == std::string::npos) break;
        ManifestEntry e;
        e.basename = content.substr(bp + 1, bp_end - bp - 1);
        auto bbox_p = content.find("\"bbox\"", bp_end);
        if (bbox_p == std::string::npos) break;
        std::string sub = content.substr(bbox_p, 400);
        find_double_in(sub, "x_lo", e.xlo); find_double_in(sub, "x_hi", e.xhi);
        find_double_in(sub, "y_lo", e.ylo); find_double_in(sub, "y_hi", e.yhi);
        find_double_in(sub, "z_lo", e.zlo); find_double_in(sub, "z_hi", e.zhi);
        out.push_back(e);
        p = bbox_p + 1;
    }
    if (out.empty()) {
        std::cerr << "manifest has no meshes[]\n"; std::exit(1);
    }
    std::sort(out.begin(), out.end(),
              [](const ManifestEntry& a, const ManifestEntry& b) {
                  return a.basename < b.basename;
              });
    return out;
}

// -----------------------------------------------------------------------------
struct Args {
    std::string manifest_path;
    std::string out_stl;
    std::string out_markers_json;
    // Boundary mode (default when any --boundary is given): the box is
    // replaced by user-supplied boundary surfaces (markers 100+i in CLI
    // order).  --box restores the flat-box fixture path.  --fault-prefix
    // filters manifest entries: only basenames starting with the prefix
    // are loaded as faults (the SAFv4 corefine manifest also lists the
    // dem/bottom participants).
    std::vector<std::string> boundary_paths;
    bool box_mode = false;
    std::string fault_prefix;
    // Soup mode: re-autorefine an EXISTING soup STL + per-triangle marker
    // JSON (e.g. the clipped merged soup after geometry edits), skipping
    // assembly entirely.  MANIFEST is still required for fault basenames.
    std::string soup_in, soup_markers;
    // pad_top = 0 by default so the box top coincides with the geological
    // free surface (z=0) where every input fault outcrops.  With pad_top
    // > 0 the fault tops would sit INSIDE the bulk, breaking free-surface
    // BCs and preventing autorefine_triangle_soup from inserting Steiner
    // points along fault outcrops on the box top face.
    double pad_xy = 50000.0, pad_top = 0.0, pad_bottom = 25000.0;
    // box_edge_size = uniform target edge length for box-face triangulation
    // BEFORE autorefine.  Refining the 12 corner triangles into a grid lets
    // autorefine_triangle_soup resolve fault-perimeter / box-top intersection
    // by inserting Steiner points where fault edges cross box-top facets,
    // producing local refinement on the free surface near each outcrop.
    double box_edge_size = 10000.0;
    int    box_marker = 100;
    bool   verbose = false;
};

static int usage(const char* argv0, int code) {
    std::cerr << "usage: " << argv0
              << " MANIFEST OUTPUT_STL OUTPUT_MARKERS_JSON"
              << " [--boundary FILE]... [--box] [--fault-prefix STR]"
              << " [--pad-xy M] [--pad-top M] [--pad-bottom M]"
              << " [--box-edge-size M] [--box-marker N] [--verbose]\n"
              << "  boundary mode (any --boundary given): the box is replaced"
              << " by the listed boundary surfaces (markers 100+i, CLI"
              << " order); faults = manifest entries whose basename starts"
              << " with --fault-prefix (markers 1..N, alphabetical).\n";
    return code;
}

static int parse(int argc, char** argv, Args& a) {
    auto need_double = [&](int& i, double& x, const std::string& n) {
        if (++i >= argc) { std::cerr << "missing value for " << n << "\n"; std::exit(2); }
        x = std::strtod(argv[i], nullptr);
    };
    auto need_int = [&](int& i, int& x, const std::string& n) {
        if (++i >= argc) { std::cerr << "missing value for " << n << "\n"; std::exit(2); }
        x = std::atoi(argv[i]);
    };
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "-h" || s == "--help") return usage(argv[0], 0);
        else if (s == "--pad-xy")     need_double(i, a.pad_xy, s);
        else if (s == "--pad-top")    need_double(i, a.pad_top, s);
        else if (s == "--pad-bottom") need_double(i, a.pad_bottom, s);
        else if (s == "--box-edge-size") need_double(i, a.box_edge_size, s);
        else if (s == "--box-marker") need_int(i, a.box_marker, s);
        else if (s == "--boundary") {
            if (++i >= argc) { std::cerr << "missing value for --boundary\n"; std::exit(2); }
            a.boundary_paths.push_back(argv[i]);
        }
        else if (s == "--box")        a.box_mode = true;
        else if (s == "--fault-prefix") {
            if (++i >= argc) { std::cerr << "missing value for --fault-prefix\n"; std::exit(2); }
            a.fault_prefix = argv[i];
        }
        else if (s == "--soup-in") {
            if (++i >= argc) { std::cerr << "missing value for --soup-in\n"; std::exit(2); }
            a.soup_in = argv[i];
        }
        else if (s == "--soup-markers") {
            if (++i >= argc) { std::cerr << "missing value for --soup-markers\n"; std::exit(2); }
            a.soup_markers = argv[i];
        }
        else if (s == "--verbose")    a.verbose = true;
        else pos.push_back(s);
    }
    if (pos.size() != 3) return usage(argv[0], 2);
    a.manifest_path    = pos[0];
    a.out_stl          = pos[1];
    a.out_markers_json = pos[2];
    return -1;
}

// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
    Args a;
    if (int r = parse(argc, argv, a); r >= 0) return r;

    fs::path manifest_p(a.manifest_path);
    fs::path corefined_dir = manifest_p.parent_path();
    auto entries = load_manifest_entries(a.manifest_path);
    if (a.boundary_paths.empty()) a.box_mode = true;   // backward compat
    if (!a.fault_prefix.empty()) {
        entries.erase(std::remove_if(entries.begin(), entries.end(),
            [&](const ManifestEntry& e) {
                return e.basename.rfind(a.fault_prefix, 0) != 0;
            }), entries.end());
        if (entries.empty()) {
            std::cerr << "no manifest entries match --fault-prefix "
                      << a.fault_prefix << "\n";
            return 2;
        }
    }

    // Compute box bounds (union of fault bboxes + paddings).
    double xlo = entries[0].xlo, xhi = entries[0].xhi;
    double ylo = entries[0].ylo, yhi = entries[0].yhi;
    double zlo = entries[0].zlo, zhi = entries[0].zhi;
    for (const auto& e : entries) {
        xlo = std::min(xlo, e.xlo); xhi = std::max(xhi, e.xhi);
        ylo = std::min(ylo, e.ylo); yhi = std::max(yhi, e.yhi);
        zlo = std::min(zlo, e.zlo); zhi = std::max(zhi, e.zhi);
    }
    const double xmin = xlo - a.pad_xy, xmax = xhi + a.pad_xy;
    const double ymin = ylo - a.pad_xy, ymax = yhi + a.pad_xy;
    const double zmin = zlo - a.pad_bottom, zmax = zhi + a.pad_top;
    if (a.verbose) {
        std::cerr << "manifest: " << a.manifest_path << "\n"
                  << "  faults: " << entries.size() << " (alphabetical)\n";
        for (std::size_t i = 0; i < entries.size(); ++i) {
            std::cerr << "    marker " << (i + 1) << ": " << entries[i].basename << "\n";
        }
        std::cerr << "  box marker: " << a.box_marker << "\n"
                  << "  box bounds: x=[" << xmin << ", " << xmax << "]"
                  << " y=[" << ymin << ", " << ymax << "]"
                  << " z=[" << zmin << ", " << zmax << "]\n";
    }

    // Build polygon soup with vertex dedup at 1e-6 m (matches stitch_combined_stl.py).
    std::vector<Point> points;
    std::vector<std::vector<std::size_t>> triangles;
    std::vector<int> markers;

    auto coord_key = [](const Point& p) {
        return std::tuple<long long, long long, long long>{
            static_cast<long long>(std::llround(p.x() * 1e6)),
            static_cast<long long>(std::llround(p.y() * 1e6)),
            static_cast<long long>(std::llround(p.z() * 1e6))};
    };
    std::map<std::tuple<long long, long long, long long>, std::size_t> point_idx;

    auto intern = [&](const Point& p) {
        auto k = coord_key(p);
        auto it = point_idx.find(k);
        if (it != point_idx.end()) return it->second;
        std::size_t idx = points.size();
        point_idx[k] = idx;
        points.push_back(p);
        return idx;
    };

    // Box: 6 faces, each refined into a uniform (nu × nv) quad grid →
    // 2*nu*nv outward-oriented triangles per face.  Pre-refining the box
    // surface gives autorefine_triangle_soup something to "cut" along the
    // fault perimeter (where fault edges intersect box-top facets at z=0
    // when pad_top=0): autorefine inserts Steiner points at every fault-
    // box intersection, producing local refinement of the box top near
    // each fault outcrop without our needing to know the outcrop shape
    // a priori.
    auto add_face = [&](int fixed_axis, double fixed_val,
                        int u_axis, double u_lo, double u_hi,
                        int v_axis, double v_lo, double v_hi,
                        int normal_sign) {
        const int nu = std::max(1, static_cast<int>(
            std::ceil((u_hi - u_lo) / a.box_edge_size)));
        const int nv = std::max(1, static_cast<int>(
            std::ceil((v_hi - v_lo) / a.box_edge_size)));
        // Local (iu, iv) → global vertex index via intern().
        std::vector<std::vector<std::size_t>> grid(nu + 1,
            std::vector<std::size_t>(nv + 1, 0));
        for (int iu = 0; iu <= nu; ++iu) {
            const double u = u_lo + (u_hi - u_lo) * iu / nu;
            for (int iv = 0; iv <= nv; ++iv) {
                const double v = v_lo + (v_hi - v_lo) * iv / nv;
                double xyz[3] = {0, 0, 0};
                xyz[fixed_axis] = fixed_val;
                xyz[u_axis] = u;
                xyz[v_axis] = v;
                grid[iu][iv] = intern(Point(xyz[0], xyz[1], xyz[2]));
            }
        }
        for (int iu = 0; iu < nu; ++iu) {
            for (int iv = 0; iv < nv; ++iv) {
                const auto a0 = grid[iu  ][iv  ];
                const auto a1 = grid[iu+1][iv  ];
                const auto a2 = grid[iu+1][iv+1];
                const auto a3 = grid[iu  ][iv+1];
                if (normal_sign > 0) {
                    triangles.push_back({a0, a1, a2});
                    triangles.push_back({a0, a2, a3});
                } else {
                    triangles.push_back({a0, a2, a1});
                    triangles.push_back({a0, a3, a2});
                }
                markers.push_back(a.box_marker);
                markers.push_back(a.box_marker);
            }
        }
    };
    std::vector<std::string> boundary_names;
    if (!a.soup_in.empty()) {
        // Soup mode: load an existing soup + markers verbatim.
        if (a.soup_markers.empty()) {
            std::cerr << "--soup-in requires --soup-markers\n"; return 2;
        }
        std::vector<Point> spts;
        std::vector<std::vector<std::size_t>> spolys;
        if (!CGAL::IO::read_polygon_soup(a.soup_in, spts, spolys)) {
            std::cerr << "cannot read soup STL: " << a.soup_in << "\n";
            return 3;
        }
        std::ifstream mf(a.soup_markers);
        std::string mcontent((std::istreambuf_iterator<char>(mf)),
                             std::istreambuf_iterator<char>());
        auto mp = mcontent.find("\"markers\"");
        mp = mcontent.find('[', mp);
        std::vector<int> smarks;
        const char* c = mcontent.c_str() + mp + 1;
        while (*c && *c != ']') {
            while (*c == ' ' || *c == ',' || *c == '\n') ++c;
            if (*c == ']' || !*c) break;
            smarks.push_back(std::atoi(c));
            while (*c && *c != ',' && *c != ']') ++c;
        }
        if (smarks.size() != spolys.size()) {
            std::cerr << "soup markers (" << smarks.size()
                      << ") != soup tris (" << spolys.size() << ")\n";
            return 4;
        }
        // Carry the fault basenames verbatim from the input soup markers
        // (the soup's marker semantics are already fixed; regenerating
        // from the unfiltered manifest would mislabel marker 1 as the
        // first ALPHABETICAL mesh, e.g. dem_graded).
        {
            auto fbp = mcontent.find("\"fault_basenames\"");
            if (fbp != std::string::npos) {
                auto lb = mcontent.find('[', fbp);
                auto rb = mcontent.find(']', lb);
                std::string arr = mcontent.substr(lb + 1, rb - lb - 1);
                std::vector<ManifestEntry> fent;
                std::size_t q0 = 0;
                while ((q0 = arr.find('"', q0)) != std::string::npos) {
                    auto q1 = arr.find('"', q0 + 1);
                    if (q1 == std::string::npos) break;
                    ManifestEntry e;
                    e.basename = arr.substr(q0 + 1, q1 - q0 - 1);
                    fent.push_back(e);
                    q0 = q1 + 1;
                }
                if (!fent.empty()) entries = fent;
            }
        }
        // Carry boundary_names through (the .msh writer needs them).
        auto bnp = mcontent.find("\"boundary_names\"");
        if (bnp != std::string::npos) {
            auto lb = mcontent.find('[', bnp);
            auto rb = mcontent.find(']', lb);
            std::string arr = mcontent.substr(lb + 1, rb - lb - 1);
            std::size_t q0 = 0;
            while ((q0 = arr.find('"', q0)) != std::string::npos) {
                auto q1 = arr.find('"', q0 + 1);
                if (q1 == std::string::npos) break;
                boundary_names.push_back(arr.substr(q0 + 1, q1 - q0 - 1));
                q0 = q1 + 1;
            }
        }
        std::size_t n_degen2 = 0;
        for (std::size_t pi2 = 0; pi2 < spolys.size(); ++pi2) {
            const auto& poly = spolys[pi2];
            if (poly.size() != 3) continue;
            const std::size_t i0 = intern(spts[poly[0]]);
            const std::size_t i1 = intern(spts[poly[1]]);
            const std::size_t i2 = intern(spts[poly[2]]);
            if (i0 == i1 || i1 == i2 || i0 == i2) { ++n_degen2; continue; }
            triangles.push_back({i0, i1, i2});
            markers.push_back(smarks[pi2]);
        }
        if (a.verbose) {
            std::cerr << "  soup mode: " << triangles.size() << " tris ("
                      << n_degen2 << " degenerate dropped)\n";
        }
    } else if (a.box_mode) {
        // axis indices: 0=x, 1=y, 2=z
        add_face(2, zmin, 0, xmin, xmax, 1, ymin, ymax, -1);  // bottom (-z)
        add_face(2, zmax, 0, xmin, xmax, 1, ymin, ymax, +1);  // top    (+z)  ← free surface
        add_face(1, ymin, 0, xmin, xmax, 2, zmin, zmax, -1);  // front  (-y)
        add_face(1, ymax, 0, xmin, xmax, 2, zmin, zmax, +1);  // back   (+y)
        add_face(0, xmin, 1, ymin, ymax, 2, zmin, zmax, -1);  // left   (-x)
        add_face(0, xmax, 1, ymin, ymax, 2, zmin, zmax, +1);  // right  (+x)
    } else {
        // Boundary mode: user-supplied shell surfaces, marker 100+i in
        // CLI order, interned through the same 1e-6 dedup so the welded
        // rims share vertices.
        for (std::size_t bi = 0; bi < a.boundary_paths.size(); ++bi) {
            const std::string& bp = a.boundary_paths[bi];
            std::vector<Point> bpts;
            std::vector<std::vector<std::size_t>> bpolys;
            if (!CGAL::IO::read_polygon_soup(bp, bpts, bpolys)) {
                std::cerr << "cannot read boundary STL: " << bp << "\n";
                return 3;
            }
            const int marker = a.box_marker + static_cast<int>(bi);
            std::size_t n_tri = 0;
            for (const auto& poly : bpolys) {
                if (poly.size() != 3) continue;
                const std::size_t i0 = intern(bpts[poly[0]]);
                const std::size_t i1 = intern(bpts[poly[1]]);
                const std::size_t i2 = intern(bpts[poly[2]]);
                if (i0 == i1 || i1 == i2 || i0 == i2) continue;
                triangles.push_back({i0, i1, i2});
                markers.push_back(marker);
                ++n_tri;
            }
            boundary_names.push_back(fs::path(bp).stem().string());
            if (a.verbose) {
                std::cerr << "  boundary marker " << marker << ": "
                          << fs::path(bp).filename().string()
                          << "  F=" << n_tri << "\n";
            }
        }
    }

    // Faults in alphabetical order (skipped in soup mode: already inside).
    for (std::size_t fi = 0; a.soup_in.empty() && fi < entries.size(); ++fi) {
        fs::path stl_p = corefined_dir / (entries[fi].basename + "_corefined.stl");
        std::vector<Point> fpts;
        std::vector<std::vector<std::size_t>> fpolys;
        if (!CGAL::IO::read_STL(stl_p.string(), fpts, fpolys)) {
            std::cerr << "cannot read STL: " << stl_p << "\n"; return 3;
        }
        std::size_t n_degen = 0;
        for (const auto& poly : fpolys) {
            if (poly.size() != 3) continue;
            const std::size_t i0 = intern(fpts[poly[0]]);
            const std::size_t i1 = intern(fpts[poly[1]]);
            const std::size_t i2 = intern(fpts[poly[2]]);
            // The 1e-6 intern can merge coincident duplicated vertices
            // (pinch corners, 0-length constrained edges): drop the
            // resulting degenerate triangles instead of feeding them to
            // autorefine.
            if (i0 == i1 || i1 == i2 || i0 == i2) { ++n_degen; continue; }
            triangles.push_back({i0, i1, i2});
            markers.push_back(static_cast<int>(fi + 1));
        }
        if (a.verbose && n_degen > 0) {
            std::cerr << "    dropped " << n_degen
                      << " degenerate triangle(s) after intern dedup\n";
        }
        if (a.verbose) {
            std::cerr << "  loaded " << stl_p.filename().string()
                      << "  V=" << fpts.size() << "  F=" << fpolys.size() << "\n";
        }
    }

    if (a.verbose) {
        std::cerr << "merged input: V=" << points.size()
                  << "  F=" << triangles.size()
                  << "  markers=" << markers.size() << "\n";
    }

    // Run autorefine_triangle_soup with marker-propagating visitor.
    std::vector<int> markers_out;
    Marker_propagating_visitor visitor;
    visitor.in = &markers;
    visitor.out = &markers_out;
    visitor.n_faults = static_cast<int>(entries.size());

    if (a.verbose) {
        std::cerr << "running PMP::autorefine_triangle_soup ...\n";
    }
    bool ok = PMP::autorefine_triangle_soup(points, triangles,
                                            pp::visitor(visitor));
    (void)ok;
    if (a.verbose) {
        std::cerr << "autorefine done: V=" << points.size()
                  << "  F=" << triangles.size()
                  << "  markers_out=" << markers_out.size() << "\n";
    }
    if (markers_out.size() != triangles.size()) {
        std::cerr << "internal error: markers_out.size() = " << markers_out.size()
                  << " != triangles.size() = " << triangles.size() << "\n";
        return 4;
    }
    // Phase 4 gate: in boundary mode all fault x fault and fault x DEM
    // intersections were corefined in Phase 3 — autorefine must not have
    // had to split any fault-marked triangle.  (Expected residual work:
    // boundary-side splits and exact-duplicate removal only.)
    if (a.verbose || visitor.n_fault_subtriangles > 0) {
        std::cerr << "fault-marked new_subtriangle count: "
                  << visitor.n_fault_subtriangles << "\n";
    }
    // Verify all markers were assigned.
    std::size_t n_unassigned = 0;
    for (int m : markers_out) if (m < 0) ++n_unassigned;
    if (n_unassigned > 0) {
        std::cerr << "warning: " << n_unassigned
                  << " output triangles have no marker (visitor missed them)\n";
    }

    // Write merged STL.
    {
        fs::path out_p(a.out_stl);
        fs::create_directories(out_p.parent_path());
        std::ofstream f(out_p);
        if (!f) { std::cerr << "cannot write " << a.out_stl << "\n"; return 5; }
        f << std::setprecision(15);
        f << "solid safs_autorefined\n";
        for (const auto& t : triangles) {
            const auto& A = points[t[0]];
            const auto& B = points[t[1]];
            const auto& C = points[t[2]];
            // Outward normal not strictly needed (tetgen ignores it) but
            // emit it for STL viewer compatibility.
            const double ux = B.x() - A.x(), uy = B.y() - A.y(), uz = B.z() - A.z();
            const double vx = C.x() - A.x(), vy = C.y() - A.y(), vz = C.z() - A.z();
            double nx = uy * vz - uz * vy;
            double ny = uz * vx - ux * vz;
            double nz = ux * vy - uy * vx;
            const double nm = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (nm > 0) { nx /= nm; ny /= nm; nz /= nm; }
            f << "facet normal " << nx << " " << ny << " " << nz << "\n";
            f << "  outer loop\n";
            f << "    vertex " << A.x() << " " << A.y() << " " << A.z() << "\n";
            f << "    vertex " << B.x() << " " << B.y() << " " << B.z() << "\n";
            f << "    vertex " << C.x() << " " << C.y() << " " << C.z() << "\n";
            f << "  endloop\n";
            f << "endfacet\n";
        }
        f << "endsolid safs_autorefined\n";
        if (a.verbose) std::cerr << "wrote " << a.out_stl << "\n";
    }

    // Write markers JSON: per-triangle integer marker.
    // Format: {"box_marker": 100,
    //          "fault_basenames": ["...", ...],   // index = marker-1
    //          "n_box_triangles_input": 12,
    //          "n_input_faults": 6,
    //          "n_output_triangles": <N>,
    //          "n_output_vertices":  <V>,
    //          "markers": [m_0, m_1, ...]}
    {
        std::ofstream f(a.out_markers_json);
        if (!f) { std::cerr << "cannot write " << a.out_markers_json << "\n"; return 6; }
        f << "{\n";
        f << "  \"box_marker\": " << a.box_marker << ",\n";
        f << "  \"mode\": \"" << (a.box_mode ? "box" : "boundary") << "\",\n";
        f << "  \"boundary_names\": [";
        for (std::size_t i = 0; i < boundary_names.size(); ++i) {
            if (i > 0) f << ", ";
            f << "\"" << boundary_names[i] << "\"";
        }
        f << "],\n";
        f << "  \"fault_basenames\": [";
        for (std::size_t i = 0; i < entries.size(); ++i) {
            if (i > 0) f << ", ";
            f << "\"" << entries[i].basename << "\"";
        }
        f << "],\n";
        f << "  \"n_input_faults\": " << entries.size() << ",\n";
        f << "  \"n_box_triangles_input\": 12,\n";
        f << "  \"n_fault_new_subtriangles\": "
          << visitor.n_fault_subtriangles << ",\n";
        f << "  \"n_output_triangles\": " << triangles.size() << ",\n";
        f << "  \"n_output_vertices\": "  << points.size()    << ",\n";
        f << "  \"markers\": [";
        for (std::size_t i = 0; i < markers_out.size(); ++i) {
            if (i > 0) f << (i % 30 == 0 ? ",\n    " : ", ");
            f << markers_out[i];
        }
        f << "]\n";
        f << "}\n";
        if (a.verbose) std::cerr << "wrote " << a.out_markers_json << "\n";
    }

    if (!a.box_mode && a.soup_in.empty() && visitor.n_fault_subtriangles > 0) {
        std::cerr << "ERROR: autorefine split " << visitor.n_fault_subtriangles
                  << " fault-marked triangle(s) — a Phase 3 corefine pair "
                     "leaked an unresolved intersection.  Outputs were "
                     "written for inspection; STOP and report.\n";
        return 9;
    }
    // Soup mode is a REPAIR pass: fault splits are the point; report only.
    return 0;
}
