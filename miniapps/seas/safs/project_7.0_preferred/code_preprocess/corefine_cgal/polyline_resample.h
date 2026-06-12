#pragma once

// polyline_resample.h — arc-length-uniform polyline resampler.
//
// Both meshes A and B in a corefine pair must see the SAME polyline
// vertex sequence after Delaunay remeshing.  Pre-resampling the polyline
// at uniform spacing guarantees that both surface_Delaunay_remeshing
// calls receive identical input; with protect_constraints=true the
// vertex sequence is preserved through the remesh, so A_out and B_out
// share polyline vertex coords bit-identically.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cmath>
#include <vector>

namespace polyline_resample {

template <class Point>
inline double dist(const Point& a, const Point& b) {
    const double dx = a.x() - b.x();
    const double dy = a.y() - b.y();
    const double dz = a.z() - b.z();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Resample a single polyline at arc-length-uniform `target_spacing`.
// Vertices closer than `min_edge_along_polyline` are merged.
//
// Closed polyline detection: first == last (within 1e-9 m).
template <class Point>
std::vector<Point>
resample(const std::vector<Point>& polyline,
         double target_spacing,
         double min_edge_along_polyline) {
    if (polyline.size() < 2) return polyline;

    const bool closed =
        (dist(polyline.front(), polyline.back()) < 1e-9);

    // Cumulative arc length.
    std::vector<double> s(polyline.size(), 0.0);
    for (std::size_t i = 1; i < polyline.size(); ++i) {
        s[i] = s[i - 1] + dist(polyline[i - 1], polyline[i]);
    }
    const double L = s.back();
    if (L < min_edge_along_polyline) {
        // Polyline is shorter than the minimum edge: collapse to endpoints.
        std::vector<Point> out;
        out.push_back(polyline.front());
        if (!closed) out.push_back(polyline.back());
        return out;
    }

    // Sample at j * target_spacing for j = 0..floor(L / target_spacing),
    // plus the final endpoint at s = L.
    std::vector<Point> sampled;
    const std::size_t n_steps =
        static_cast<std::size_t>(std::floor(L / target_spacing));
    for (std::size_t j = 0; j <= n_steps; ++j) {
        const double sj = j * target_spacing;
        // Locate segment containing sj.
        auto it = std::upper_bound(s.begin(), s.end(), sj);
        std::size_t k = (it == s.begin()) ? 0
                                          : std::distance(s.begin(), it) - 1;
        if (k >= polyline.size() - 1) k = polyline.size() - 2;
        const double seg_len = s[k + 1] - s[k];
        const double t = (seg_len > 0.0) ? (sj - s[k]) / seg_len : 0.0;
        const auto& a = polyline[k];
        const auto& b = polyline[k + 1];
        sampled.emplace_back(
            a.x() + t * (b.x() - a.x()),
            a.y() + t * (b.y() - a.y()),
            a.z() + t * (b.z() - a.z()));
    }
    // Always include the actual last point unless we'd duplicate it.
    if (!closed) {
        if (sampled.empty() || dist(sampled.back(), polyline.back()) > 1e-9) {
            sampled.push_back(polyline.back());
        }
    } else {
        // For closed polylines, ensure the sequence ends with the first
        // vertex (closing the loop).
        if (sampled.empty() || dist(sampled.back(), polyline.front()) > 1e-9) {
            sampled.push_back(polyline.front());
        }
    }

    // Merge sub-edges < min_edge_along_polyline.  Walk forward; if the
    // distance to the previously-kept vertex is below the floor, drop the
    // current vertex (keep the running anchor).  Always keep the final
    // vertex (its drop would change topology).
    if (sampled.size() < 3) return sampled;
    std::vector<Point> merged;
    merged.reserve(sampled.size());
    merged.push_back(sampled.front());
    for (std::size_t i = 1; i + 1 < sampled.size(); ++i) {
        if (dist(merged.back(), sampled[i]) >= min_edge_along_polyline) {
            merged.push_back(sampled[i]);
        }
    }
    // Make sure the final vertex is at >= min_edge from the previous.
    if (dist(merged.back(), sampled.back()) < min_edge_along_polyline
        && merged.size() > 1) {
        // Pop the previous-kept vertex to avoid a sub-floor closing edge.
        merged.pop_back();
    }
    merged.push_back(sampled.back());
    return merged;
}

template <class Point>
std::vector<std::vector<Point>>
resample_all(const std::vector<std::vector<Point>>& polylines,
             double target_spacing,
             double min_edge_along_polyline) {
    std::vector<std::vector<Point>> out;
    out.reserve(polylines.size());
    for (const auto& pl : polylines) {
        if (pl.size() < 2) continue;  // skip degenerate
        auto r = resample(pl, target_spacing, min_edge_along_polyline);
        if (r.size() >= 2) out.push_back(std::move(r));
    }
    return out;
}

}  // namespace polyline_resample
