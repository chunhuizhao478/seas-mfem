// check_self_intersect — diagnostic: report any self-intersecting face pairs
// in each *.stl in IN_DIR.  Used to identify whether corefine_set's outputs
// have geometric self-intersections (the gmsh PLC error symptom).
//
// Usage: check_self_intersect IN_DIR [--ext .stl]

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include "io_helpers.h"

#include <filesystem>
#include <iostream>
#include <vector>

namespace fs = std::filesystem;
using K     = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh  = CGAL::Surface_mesh<K::Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " IN_DIR [--ext .stl]\n";
        return 2;
    }
    std::string ext = ".stl";
    for (int i = 2; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "--ext" && i + 1 < argc) { ext = argv[++i]; }
    }
    std::vector<std::string> paths;
    for (const auto& e : fs::directory_iterator(argv[1])) {
        if (e.is_regular_file() && e.path().extension() == ext)
            paths.push_back(e.path().string());
    }
    std::sort(paths.begin(), paths.end());

    int total_si = 0;
    for (const auto& p : paths) {
        Mesh m;
        if (!io_helpers::read_polygon_mesh_any(p, m)) {
            std::cerr << "could not read " << p << "\n";
            return 3;
        }
        std::vector<std::pair<Mesh::Face_index, Mesh::Face_index>> inters;
        PMP::self_intersections(faces(m), m, std::back_inserter(inters));
        total_si += static_cast<int>(inters.size());
        std::cout << fs::path(p).stem().string()
                  << "  V=" << m.number_of_vertices()
                  << "  F=" << m.number_of_faces()
                  << "  self_intersect_pairs=" << inters.size() << "\n";
        // Print the first few intersection pairs and their centroid coords.
        std::size_t shown = 0;
        for (const auto& [fa, fb] : inters) {
            if (shown >= 5) { std::cout << "  ... (+ more)\n"; break; }
            auto centroid = [&](Mesh::Face_index f) {
                double x = 0, y = 0, z = 0;
                for (auto v : CGAL::vertices_around_face(m.halfedge(f), m)) {
                    auto& pt = m.point(v);
                    x += pt.x(); y += pt.y(); z += pt.z();
                }
                return std::array<double, 3>{x / 3, y / 3, z / 3};
            };
            auto ca = centroid(fa);
            auto cb = centroid(fb);
            std::cout << "  pair faces " << fa << " <-> " << fb
                      << "  cA=(" << ca[0] << "," << ca[1] << "," << ca[2] << ")"
                      << "  cB=(" << cb[0] << "," << cb[1] << "," << cb[2] << ")\n";
            ++shown;
        }
    }
    std::cout << "TOTAL self-intersection pairs across all meshes: "
              << total_si << "\n";
    return total_si > 0 ? 1 : 0;
}
