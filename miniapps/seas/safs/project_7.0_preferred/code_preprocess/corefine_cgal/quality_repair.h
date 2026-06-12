#pragma once

// quality_repair.h — degenerate-face removal + statistics for the
// corefine pipeline.  Stops at PMP::remove_almost_degenerate_faces and
// does NOT call PMP::tangential_relaxation: the latter would move
// polyline vertices unless an explicit vertex_is_constrained_map were
// passed, and the marginal smoothing gain does not justify the
// conformality risk (see PLAN_cgal_corefine_multifault.md R-001).

#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Surface_mesh.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace quality_repair {

namespace PMP = CGAL::Polygon_mesh_processing;

struct Stats {
    std::size_t n_verts = 0;
    std::size_t n_faces = 0;
    double      edge_len_min   = 0.0;
    double      edge_len_p1    = 0.0;
    double      edge_len_med   = 0.0;
    double      edge_len_max   = 0.0;
    double      tri_q_min      = 0.0;
    double      tri_q_p1       = 0.0;
    double      tri_q_med      = 0.0;
    std::size_t n_edges_below_floor = 0;
    std::size_t n_tri_q_below_threshold = 0;
};

// q_tri = 4 * sqrt(3) * area / sum_of_squared_edge_lengths
// 1.0 = equilateral, 0.0 = degenerate.  Matches msh_to_vtu.py.
template <class Point>
double tri_quality(const Point& a, const Point& b, const Point& c) {
    const double ab2 = CGAL::squared_distance(a, b);
    const double bc2 = CGAL::squared_distance(b, c);
    const double ca2 = CGAL::squared_distance(c, a);
    const double sum_sq = ab2 + bc2 + ca2;
    if (sum_sq <= 0.0) return 0.0;
    // Area via cross product.
    const double ux = b.x() - a.x();
    const double uy = b.y() - a.y();
    const double uz = b.z() - a.z();
    const double vx = c.x() - a.x();
    const double vy = c.y() - a.y();
    const double vz = c.z() - a.z();
    const double cx = uy * vz - uz * vy;
    const double cy = uz * vx - ux * vz;
    const double cz = ux * vy - uy * vx;
    const double area = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
    return 4.0 * std::sqrt(3.0) * area / sum_sq;
}

template <class Mesh>
Stats compute_stats(const Mesh& m,
                    double min_edge_size,
                    double q_floor) {
    Stats s;
    s.n_verts = m.number_of_vertices();
    s.n_faces = m.number_of_faces();
    if (s.n_faces == 0) return s;

    std::vector<double> edges;
    edges.reserve(s.n_faces * 3);
    std::vector<double> qs;
    qs.reserve(s.n_faces);

    for (auto e : m.edges()) {
        const auto h = m.halfedge(e);
        const auto& a = m.point(m.source(h));
        const auto& b = m.point(m.target(h));
        edges.push_back(std::sqrt(CGAL::squared_distance(a, b)));
    }
    for (auto f : m.faces()) {
        std::vector<typename Mesh::Point> tri_pts;
        for (auto v : CGAL::vertices_around_face(m.halfedge(f), m)) {
            tri_pts.push_back(m.point(v));
        }
        if (tri_pts.size() == 3) {
            qs.push_back(tri_quality(tri_pts[0], tri_pts[1], tri_pts[2]));
        }
    }

    auto percentile = [](std::vector<double>& v, double p) -> double {
        if (v.empty()) return 0.0;
        const std::size_t n = v.size();
        const std::size_t k =
            std::min<std::size_t>(n - 1,
                                  static_cast<std::size_t>(p * (n - 1)));
        std::nth_element(v.begin(), v.begin() + k, v.end());
        return v[k];
    };

    if (!edges.empty()) {
        std::vector<double> e_sorted = edges;
        std::sort(e_sorted.begin(), e_sorted.end());
        s.edge_len_min = e_sorted.front();
        s.edge_len_max = e_sorted.back();
        s.edge_len_p1  = percentile(e_sorted, 0.01);
        s.edge_len_med = percentile(e_sorted, 0.50);
        for (double e : edges) if (e < min_edge_size) ++s.n_edges_below_floor;
    }
    if (!qs.empty()) {
        std::vector<double> q_sorted = qs;
        std::sort(q_sorted.begin(), q_sorted.end());
        s.tri_q_min = q_sorted.front();
        s.tri_q_p1  = percentile(q_sorted, 0.01);
        s.tri_q_med = percentile(q_sorted, 0.50);
        for (double q : qs) if (q < q_floor) ++s.n_tri_q_below_threshold;
    }
    return s;
}

// Iterate PMP::remove_almost_degenerate_faces until no more changes or
// max_iterations reached.  ECM is the constraint map marking the
// resampled polyline edges.  VCM is the vertex-constraint map marking
// every endpoint of any polyline edge — required to prevent
// remove_almost_degenerate_faces from collapsing an unconstrained edge
// incident on a polyline vertex (which would geometrically MOVE the
// polyline endpoint and break pairwise conformality across the pair).
template <class Mesh, class ECM, class VCM>
Stats run(Mesh& m,
          ECM   ecm,
          VCM   vcm,
          double min_edge_size,
          double q_floor                = 0.3,
          double cap_threshold_cos      = std::cos(160.0 * 3.14159265358979323846 / 180.0),
          double needle_threshold_ratio = 4.0,
          unsigned max_iterations       = 5) {
    namespace pp = PMP::parameters;
    for (unsigned iter = 0; iter < max_iterations; ++iter) {
        std::size_t f_before = m.number_of_faces();
        PMP::remove_almost_degenerate_faces(
            faces(m), m,
            pp::edge_is_constrained_map(ecm)
              .vertex_is_constrained_map(vcm)
              .needle_threshold(needle_threshold_ratio)
              .cap_threshold(cap_threshold_cos)
              .collapse_length_threshold(1.5 * min_edge_size));
        m.collect_garbage();
        if (m.number_of_faces() == f_before) break;
    }
    return compute_stats(m, min_edge_size, q_floor);
}

template <class Mesh, class ECM>
bool validate(const Mesh& m, ECM /*ecm*/,
              double min_edge_size, double q_floor) {
    Stats s = compute_stats(m, min_edge_size, q_floor);
    return (s.n_edges_below_floor == 0) && (s.tri_q_min >= q_floor);
}

// Helper: build an all-false dynamic vertex constraint map.  Used by
// callers that want the legacy single-ECM `run` behaviour (no vertex
// pinning).  The map is keyed by Mesh's vertex_descriptor.
template <class Mesh>
auto make_empty_vcm(const Mesh& m) {
    auto vcm = get(CGAL::dynamic_vertex_property_t<bool>(), m);
    for (auto v : m.vertices()) put(vcm, v, false);
    return vcm;
}

}  // namespace quality_repair
