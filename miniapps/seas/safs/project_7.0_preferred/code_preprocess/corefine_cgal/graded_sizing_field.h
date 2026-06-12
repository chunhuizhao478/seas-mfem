#pragma once

// graded_sizing_field.h — distance-graded target edge length for
// PMP::isotropic_remeshing (model of PMPSizingField, CGAL >= 5.6 API,
// mirrors Uniform_sizing_field / Adaptive_sizing_field of CGAL 6.1).
//
//   h(d) = h_near + (h_far - h_near) * clamp((d - d_near)/(d_far - d_near), 0, 1)
//
// where d = Euclidean distance to the nearest "site" (trace polyline
// vertex), queried through a CGAL Kd-tree.  Edges longer than 4/3 h(mid)
// split; shorter than 4/5 h(mid) collapse — the same band semantics as
// the uniform field, evaluated at the edge midpoint.

#include <CGAL/Polygon_mesh_processing/internal/Sizing_field_base.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

template <class PolygonMesh,
          class VPMap = typename boost::property_map<PolygonMesh,
                            CGAL::vertex_point_t>::const_type>
class Graded_polyline_sizing_field
    : public CGAL::Polygon_mesh_processing::internal::
          Sizing_field_base<PolygonMesh, VPMap>
{
    typedef CGAL::Polygon_mesh_processing::internal::
        Sizing_field_base<PolygonMesh, VPMap> Base;

public:
    typedef typename Base::FT                  FT;
    typedef typename Base::Point_3             Point_3;
    typedef typename Base::halfedge_descriptor halfedge_descriptor;
    typedef typename Base::vertex_descriptor   vertex_descriptor;

private:
    typedef typename CGAL::Kernel_traits<Point_3>::Kernel  GK;
    typedef CGAL::Search_traits_3<GK>                      SearchTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits> Neighbor_search;
    typedef typename Neighbor_search::Tree                 Tree;

public:
    Graded_polyline_sizing_field(const FT h_near, const FT h_far,
                                 const FT d_near, const FT d_far,
                                 const std::vector<Point_3>& sites,
                                 const PolygonMesh& pmesh)
        : m_h_near(h_near), m_h_far(h_far),
          m_d_near(d_near), m_d_far(d_far),
          m_vpmap(get(CGAL::vertex_point, pmesh)),
          m_tree(std::make_shared<Tree>(sites.begin(), sites.end()))
    {
        m_tree->build();
    }

    FT target_at(const Point_3& p) const {
        Neighbor_search search(*m_tree, p, 1);
        if (search.begin() == search.end()) return m_h_far;
        const FT d = std::sqrt(CGAL::to_double(search.begin()->second));
        if (d <= m_d_near) return m_h_near;
        if (d >= m_d_far)  return m_h_far;
        const FT t = (d - m_d_near) / (m_d_far - m_d_near);
        return m_h_near + (m_h_far - m_h_near) * t;
    }

    FT at(const vertex_descriptor v, const PolygonMesh& /*pmesh*/) const {
        return target_at(get(m_vpmap, v));
    }

    std::optional<FT> is_too_long(const vertex_descriptor va,
                                  const vertex_descriptor vb,
                                  const PolygonMesh& /*pmesh*/) const {
        const Point_3& pa = get(m_vpmap, va);
        const Point_3& pb = get(m_vpmap, vb);
        const FT sqlen = FT(CGAL::squared_distance(pa, pb));
        const FT h = target_at(CGAL::midpoint(pa, pb));
        const FT sq_long = CGAL::square(FT(4.0 / 3.0) * h);
        if (sqlen > sq_long) return sqlen / sq_long;
        return std::nullopt;
    }

    std::optional<FT> is_too_short(const halfedge_descriptor h,
                                   const PolygonMesh& pmesh) const {
        const Point_3& pa = get(m_vpmap, target(h, pmesh));
        const Point_3& pb = get(m_vpmap, source(h, pmesh));
        const FT sqlen = FT(CGAL::squared_distance(pa, pb));
        const FT ht = target_at(CGAL::midpoint(pa, pb));
        const FT sq_short = CGAL::square(FT(4.0 / 5.0) * ht);
        if (sqlen < sq_short) return sqlen / sq_short;
        return std::nullopt;
    }

    Point_3 split_placement(const halfedge_descriptor h,
                            const PolygonMesh& pmesh) const {
        return CGAL::midpoint(get(m_vpmap, target(h, pmesh)),
                              get(m_vpmap, source(h, pmesh)));
    }

    void register_split_vertex(const vertex_descriptor /*v*/,
                               const PolygonMesh& /*pmesh*/) {}

private:
    const FT m_h_near, m_h_far, m_d_near, m_d_far;
    const VPMap m_vpmap;
    std::shared_ptr<Tree> m_tree;
};
