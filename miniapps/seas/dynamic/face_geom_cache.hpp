// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// face_geom_cache.hpp — non-fault interior face geometry/shape cache
// (document/code_optimization_dev/action_plan_faceflux_cache_2026-06-24.md).
//
// Builds, ONCE at setup, the per-(interior non-fault face, quadrature point)
// geometry that WaveOperator::ComputeADERFaceFluxRHS otherwise recomputes every
// macro-step:
//
//     n_q   = unit outward normal of Elem1 at QP q   (CalcOrtho, normalised)
//     w_q   = ip.weight * |CalcOrtho(J_face)|         (physical quad weight)
//     φ1_q  = Elem1 shape values at QP q              (CalcShape on Loc1)
//     φ2_q  = Elem2 shape values at QP q              (CalcShape on Loc2)
//
// On an affine (straight-sided) tet mesh these are independent of the time step,
// so caching them removes the per-step GetFaceElementTransformations / CalcOrtho
// / CalcShape from the corrector's inner loop while leaving the flux dispatch
// (InteriorFaceFlux_) untouched.
//
// Geometry-only: no NUM_STATE, no material, no flux — a single cache serves both
// WaveOperator and BimaterialWaveOperator (mirrors elem_derivative_cache.hpp).
//
// REVIEW R-002 round-off note: like the deriv-cache, the cached apply is NOT
// bit-identical to the on-the-fly kernel (storing the normalised normal / shape
// re-associates round-off); the equivalence gate is <= 1e-12 relative.

#ifndef MFEM_SEAS_FACE_GEOM_CACHE_HPP
#define MFEM_SEAS_FACE_GEOM_CACHE_HPP

#include "mfem.hpp"

#include <array>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Cached geometry for one interior non-fault face.
struct FaceGeomEntry
{
   int face_index = -1;          ///< mesh face index (key for InteriorFaceFlux_)
   int e1 = -1, e2 = -1;         ///< Elem1No / Elem2No (e2 >= 0: interior)
   int ndof1 = 0, ndof2 = 0;     ///< dofs/element (homogeneous order => equal)
   int nqp = 0;
   std::vector<std::array<real_t, 3>> nor;     ///< [q] UNIT normal of e1 (signed)
   std::vector<real_t>                w;       ///< [q] ip.weight * |CalcOrtho|
   std::vector<mfem::Vector>          shape1;  ///< [q] phi^{e1}(x_q), size ndof1
   std::vector<mfem::Vector>          shape2;  ///< [q] phi^{e2}(x_q), size ndof2
};

/// @brief Build the geometry cache for every interior non-fault face.
///
/// The face-selection predicate and the per-QP geometry computation replicate
/// WaveOperator::ComputeADERFaceFluxRHS (wave_operator.inl:3962-4008, :4036-4043)
/// EXACTLY, so the cached values reproduce the on-the-fly kernel to round-off and
/// the cached set is precisely the set the corrector's fast-path will own.
///
/// Boundary faces (e2 < 0), fault faces (bdr_attr == fault_attr), and cross-rank
/// shared faces are EXCLUDED — they keep the on-the-fly path.
///
/// @param mesh                 Mesh (non-const: GetFaceElementTransformations).
/// @param fes                  L2 FE space (homogeneous order).
/// @param face_bdr_attr        Per-face boundary attribute (0 on interior),
///                             == WaveOperator::face_bdr_attr_.
/// @param fault_attr           bc_.fault_attr (0 => no fault).
/// @param shared_mesh_face_set Cross-rank shared mesh-face indices (ParMesh).
/// @param quad_order           Integration-rule order; pass 2*fe_order to match
///                             the non-fault on-the-fly rule.
/// @param[out] out             face_index -> FaceGeomEntry (cleared first).
inline void BuildNonFaultInteriorFaceGeomCache(
   mfem::Mesh &mesh,
   const mfem::FiniteElementSpace &fes,
   const std::vector<int> &face_bdr_attr,
   int fault_attr,
   const std::set<int> &shared_mesh_face_set,
   int quad_order,
   std::unordered_map<int, FaceGeomEntry> &out)
{
   out.clear();
   MFEM_VERIFY(static_cast<int>(face_bdr_attr.size()) == mesh.GetNumFaces(),
               "BuildNonFaultInteriorFaceGeomCache: face_bdr_attr size "
               << face_bdr_attr.size() << " != mesh.GetNumFaces() "
               << mesh.GetNumFaces() << " (call after face_bdr_attr_ is built).");

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      mfem::FaceElementTransformations *ftr =
         mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      const int bdr_attr = face_bdr_attr[f];

      // --- predicate: identical to ComputeADERFaceFluxRHS:3972-3986 ---
      if (e2 < 0 && shared_mesh_face_set.count(f) > 0) { continue; }
      const bool is_boundary = (e2 < 0) && (bdr_attr > 0);
      if (e2 < 0 && bdr_attr == 0) { continue; }
      if (is_boundary) { continue; }          // boundary -> on-the-fly path
      if (e2 < 0) { continue; }               // defensive: any other e2<0
      const bool face_is_fault =
         (bdr_attr == fault_attr) && (fault_attr > 0);
      if (face_is_fault) { continue; }        // fault -> on-the-fly path
      // qualifies: interior, non-fault, non-shared.

      const mfem::FiniteElement *fe1 = fes.GetFE(e1);
      const mfem::FiniteElement *fe2 = fes.GetFE(e2);
      MFEM_VERIFY(fe1->GetDof() == fe2->GetDof(),
                  "BuildNonFaultInteriorFaceGeomCache: heterogeneous FE order "
                  "(ndof1=" << fe1->GetDof() << " ndof2=" << fe2->GetDof()
                  << ") unsupported; the operator assumes homogeneous order.");

      FaceGeomEntry ent;
      ent.face_index = f;
      ent.e1 = e1;
      ent.e2 = e2;
      ent.ndof1 = fe1->GetDof();
      ent.ndof2 = fe2->GetDof();

      const mfem::IntegrationRule &ir =
         mfem::IntRules.Get(ftr->GetGeometryType(), quad_order);
      ent.nqp = ir.GetNPoints();
      ent.nor.resize(ent.nqp);
      ent.w.resize(ent.nqp);
      ent.shape1.resize(ent.nqp);
      ent.shape2.resize(ent.nqp);

      for (int q = 0; q < ent.nqp; q++)
      {
         const mfem::IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         mfem::Vector nor_vec(3);
         mfem::CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0) { nor_vec /= nor_len; }
         ent.w[q] = ip.weight * nor_len;
         ent.nor[q] = { nor_vec(0), nor_vec(1), nor_vec(2) };

         mfem::IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         ent.shape1[q].SetSize(ent.ndof1);
         fe1->CalcShape(ip1, ent.shape1[q]);

         mfem::IntegrationPoint ip2;
         ftr->Loc2.Transform(ip, ip2);
         ent.shape2[q].SetSize(ent.ndof2);
         fe2->CalcShape(ip2, ent.shape2[q]);
      }

      out.emplace(f, std::move(ent));
   }
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FACE_GEOM_CACHE_HPP
