#pragma once

// intersection_graph.h — discover intersecting mesh pairs and order them
// by descending bbox-overlap volume with a deterministic tie-breaker.

#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace intersection_graph {

namespace PMP = CGAL::Polygon_mesh_processing;

template <class Mesh>
std::vector<std::pair<std::size_t, std::size_t>>
find_intersecting_pairs(const std::vector<Mesh>& meshes) {
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    PMP::intersecting_meshes(meshes, std::back_inserter(pairs));
    // Normalize ordering: i < j.
    for (auto& p : pairs) {
        if (p.first > p.second) std::swap(p.first, p.second);
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}

template <class Mesh>
double bbox_overlap_volume(const Mesh& a, const Mesh& b) {
    auto ba = PMP::bbox(a);
    auto bb = PMP::bbox(b);
    const double ix = std::max(0.0,
        std::min(ba.xmax(), bb.xmax()) - std::max(ba.xmin(), bb.xmin()));
    const double iy = std::max(0.0,
        std::min(ba.ymax(), bb.ymax()) - std::max(ba.ymin(), bb.ymin()));
    const double iz = std::max(0.0,
        std::min(ba.zmax(), bb.zmax()) - std::max(ba.zmin(), bb.zmin()));
    return ix * iy * iz;
}

// Sort pairs by descending bbox-overlap volume.  Tie-breaker (per
// PLAN_cgal_corefine_multifault.md R-007): when |vol_a - vol_b| <
// 1e-12 * max(|vol_a|, |vol_b|), order by lexicographic ascending
// (i, j).  Makes the corefine sequence reproducible across compilers.
template <class Mesh>
void sort_by_overlap_volume(
    std::vector<std::pair<std::size_t, std::size_t>>& pairs,
    const std::vector<Mesh>& meshes) {
    std::vector<double> vols(pairs.size(), 0.0);
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        vols[k] = bbox_overlap_volume(meshes[pairs[k].first],
                                       meshes[pairs[k].second]);
    }
    // Index sort (so we can reorder both `pairs` and `vols`).
    std::vector<std::size_t> idx(pairs.size());
    for (std::size_t k = 0; k < pairs.size(); ++k) idx[k] = k;
    std::sort(idx.begin(), idx.end(),
        [&](std::size_t a, std::size_t b) {
            const double va = vols[a], vb = vols[b];
            const double maxv = std::max(std::abs(va), std::abs(vb));
            if (std::abs(va - vb) > 1e-12 * std::max(1.0, maxv)) {
                return va > vb;  // descending
            }
            return pairs[a] < pairs[b];  // lexicographic tie-breaker
        });
    std::vector<std::pair<std::size_t, std::size_t>> sorted_pairs(pairs.size());
    for (std::size_t k = 0; k < idx.size(); ++k) sorted_pairs[k] = pairs[idx[k]];
    pairs.swap(sorted_pairs);
}

}  // namespace intersection_graph
