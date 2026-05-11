#pragma once

// distance_sizing_field.h — spatially-varying sizing function for CGAL
// Mesh_3.  Returns a size proportional to the distance from the query
// point to the nearest fault surface:
//
//   d ≤ dist_inner             → lc_near
//   dist_inner ≤ d ≤ dist_outer → linear ramp lc_near → lc_far
//   d ≥ dist_outer             → lc_far
//
// Implements the gmsh-plan-§Phase-3.5a Distance+Threshold semantics in
// CGAL Mesh_3.  Without this field, mesh_volume.cpp produced a uniformly
// fine ~1500 m mesh across the entire 357×247×43 km box (5M tets); the
// graded version is expected to drop the bulk count by 50-100x while
// preserving near-fault resolution.

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

template <class GT, class Polyhedron, class Index_>
class Distance_sizing_field {
public:
    using FT      = typename GT::FT;
    using Point_3 = typename GT::Point_3;
    using Index   = Index_;

private:
    using Primitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using Traits    = CGAL::AABB_traits_3<GT, Primitive>;
    using Tree      = CGAL::AABB_tree<Traits>;

    // AABB_tree is non-copyable but Mesh_3 internally copies sizing fields,
    // so the tree is held via shared_ptr (cheap copy, shared ownership).
    // The tree references faces of the input polyhedra; those polyhedra
    // MUST outlive the field — the caller keeps them in `all_polys` until
    // make_mesh_3 returns.
    std::shared_ptr<Tree> tree_;
    double lc_min_, lc_near_, lc_far_;
    double dist_inner_, dist_outer_;

public:
    /// @param fault_polys range of polyhedra to compute distance against.
    ///                    Pass a vector or initializer-list of `const Polyhedron*`.
    /// @param lc_min       hard lower bound on the returned size
    /// @param lc_near      target size near (within dist_inner of) any fault
    /// @param lc_far       target size far (beyond dist_outer of) any fault
    /// @param dist_inner   distance below which the size is lc_near
    /// @param dist_outer   distance above which the size is lc_far
    template <class PolyhedronPointerRange>
    Distance_sizing_field(const PolyhedronPointerRange& fault_polys,
                          double lc_min, double lc_near, double lc_far,
                          double dist_inner, double dist_outer)
      : tree_(std::make_shared<Tree>()),
        lc_min_(lc_min), lc_near_(lc_near), lc_far_(lc_far),
        dist_inner_(dist_inner), dist_outer_(dist_outer)
    {
        for (const Polyhedron* p : fault_polys) {
            tree_->insert(faces(*p).first, faces(*p).second, *p);
        }
        tree_->build();
        tree_->accelerate_distance_queries();
    }

    FT operator()(const Point_3& p, int /*dim*/, const Index& /*idx*/) const {
        const double d = std::sqrt(CGAL::to_double(tree_->squared_distance(p)));
        double s;
        if (d <= dist_inner_)       s = lc_near_;
        else if (d >= dist_outer_)  s = lc_far_;
        else {
            const double t = (d - dist_inner_) / (dist_outer_ - dist_inner_);
            s = lc_near_ + t * (lc_far_ - lc_near_);
        }
        // Hard lower bound: never request a size smaller than lc_min.
        return FT(std::max(s, lc_min_));
    }
};
