// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dg_absorbing_bc_integrator.hpp — Nitsche dynamic rupture, plan item N-06.
//
// First-order absorbing (Lysmer-Kuhlemeyer viscous dashpot) boundary for the
// explicit IP elastodynamics.  On an absorbing face with outward unit normal
// n the traction is  t = -rho ( c_p (v.n) n + c_s (v - (v.n) n) ),  which adds
// the velocity-only damping term C v to  M a + C v + K u = F.  For each
// absorbing face F of element e the element block is
//
//    C_e(c n_e + i, d n_e + j) = int_F  rho_e c_p,e (phi_i n_c)(phi_j n_d)
//                                     + rho_e c_s,e phi_i phi_j (delta_cd - n_c n_d)
//
// (size D n_e x D n_e, D = the element dimension 3 or 2, element-local
// byNODES layout: component c of local dof i at c n_e + i, matching
// NitscheMass).  In 2D the tangential part v - (v.n) n is one component
// along rot90(n); the expression is the same.  rho, c_p, c_s are the element's
// values from NitscheMaterialSample (not the operator's fallback scalars).
// C is not part of K.
//
// Scope (Chunhui Zhao, 2026-10-14): tetrahedral meshes only, triangles in 2D
// (2D in-plane design rev 3, phase C); PML is future work.

#ifndef MFEM_SEAS_DG_ABSORBING_BC_INTEGRATOR_HPP
#define MFEM_SEAS_DG_ABSORBING_BC_INTEGRATOR_HPP

#include "mfem.hpp"

#include <cmath>

namespace mfem
{
namespace seas
{

class DGAbsorbingBCIntegrator
{
public:
   /// Add the dashpot block of one boundary face to C (D nd x D nd).
   /// `fe` is the scalar finite element of Tr.Elem1; Tr is the boundary
   /// face transformation (Elem2 absent).  Quadrature: face rule of order
   /// 2 p (exact for affine faces).
   static void AddFaceBlock(const FiniteElement &fe,
                            FaceElementTransformations &Tr,
                            real_t rho, real_t cp, real_t cs,
                            DenseMatrix &C)
   {
      const int nd = fe.GetDof();
      const int dim = fe.GetDim();
      MFEM_VERIFY(C.Height() == dim * nd && C.Width() == dim * nd,
                  "DGAbsorbingBCIntegrator: block must be " << dim << " nd x "
                  << dim << " nd");
      const IntegrationRule &ir =
         IntRules.Get(Tr.GetGeometryType(), 2 * fe.GetOrder());
      Vector shape(nd), nor(dim);
      for (int q = 0; q < ir.GetNPoints(); ++q)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetAllIntPoints(&ip);
         fe.CalcShape(Tr.GetElement1IntPoint(), shape);
         CalcOrtho(Tr.Jacobian(), nor);
         const real_t nl = nor.Norml2();
         nor /= nl;
         const real_t w = ip.weight * nl;        // |J_face| = |n_raw|
         const real_t zp = rho * cp, zs = rho * cs;
         for (int c = 0; c < dim; ++c)
         {
            for (int d = 0; d < dim; ++d)
            {
               const real_t kcd = zp * nor(c) * nor(d)
                                  + zs * ((c == d ? 1.0 : 0.0) - nor(c) * nor(d));
               for (int i = 0; i < nd; ++i)
               {
                  const real_t wi = w * kcd * shape(i);
                  for (int j = 0; j < nd; ++j)
                  {
                     C(c * nd + i, d * nd + j) += wi * shape(j);
                  }
               }
            }
         }
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DG_ABSORBING_BC_INTEGRATOR_HPP
