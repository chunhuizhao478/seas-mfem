// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_stress_projection.hpp — the projection of a global stress field on
// the fault DOFs: ONE function for the dynamic (Nitsche) operator's prestress
// and nucleation and for the fully dynamic cycle's tensor prestress (Chunhui
// Zhao, 2026-10-29: "we should keep the function the same as in dynamic
// phase, to make it consistent and work for 2D and 3D").
//
// ProjectStressToFaultKernel(E, sig, K): on every local fault face of the
// elasticity operator E (GetFaultInteriorFaces, then GetFaultSharedFaces),
// sigma(x_q) n_q at the face quadrature points (rule 2 p + 1), with n_q the
// face's unit normal CalcOrtho / nl_q, the normal ComputeTractionAtQuadPoints
// uses for the kernel's traction, projected on the face DOFs with
// ProjectTractionToFaultDOFs (the same rule, basis, nl_q weights and kernel
// triad as the kernel's traction).  K holds, per local fault DOF k (face
// k / nbf), the kernel components of that traction:
//   K(3 k)     = tau_1   the dip slot
//   K(3 k + 1) = tau_2   the strike slot
//   K(3 k + 2) = sigma_n compression positive (minus the normal component)
// the kernel's own convention for its traction and its prestress.  The
// Nitsche operator turns them into global three-vectors with FromKernel and
// eps_f (NitscheWaveOperator::ProjectStress_, whose body this was until
// 2026-10-29, arithmetic unchanged); the cycle hands the owned DOFs' values
// to the QD kernel as tau0 and sigma0 and, through hook H2, to the dynamic
// regime (cycle/prestress.hpp, BuildTensorPrestress).
//
// sig(x, s): the Cauchy tensor at the global point x (z = 0 on a 2-D mesh),
// row major, tension positive.  On a 2-D mesh (the in-plane model) an
// anti-plane component (sigma_xz, sigma_yz) is refused: it has no unknown to
// act on.
//
// Scope (Chunhui Zhao, 2026-10-29): the faces of
// ElasticityDomainOperator<ParMesh>, 2-D and 3-D.

#ifndef MFEM_SEAS_FAULT_STRESS_PROJECTION_HPP
#define MFEM_SEAS_FAULT_STRESS_PROJECTION_HPP

#include "mfem.hpp"
#include "elasticity_operator.hpp"

#include <functional>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

/// See the file comment.  K: 3 values per local fault DOF of E.
inline void ProjectStressToFaultKernel(ElasticityDomainOperator<ParMesh> &E,
                                       const std::function<void(const real_t[3], real_t[9])> &sig,
                                       Vector &K)
{
   ParMesh &mesh = E.GetMesh();
   const int dim = mesh.Dimension();
   const int order = E.GetOrder();
   const int nbf = E.GetNbfPerFace();
   const int nd = E.GetNumFaultDOFs();
   const Array<int> &fint = E.GetFaultInteriorFaces();
   const Array<int> &fsh = E.GetFaultSharedFaces();
   const int nfaces = fint.Size() + fsh.Size();
   MFEM_VERIFY(nfaces * nbf == nd, "ProjectStressToFaultKernel: " << nfaces << " fault faces x "
               << nbf << " DOFs per face != " << nd << " local fault DOFs");
   K.SetSize(3 * nd);
   K = 0.0;
   if (nd == 0) { return; }
   const FaceQuadrature *fq = E.GetFaceQuadrature();
   const DenseMatrix &e_q = fq->BasisAtQuadPoints();
   const FaultBasis *fb = E.GetFaultBasis();
   Vector nor(dim), xq, Tq, nlq, out;
   real_t x3[3];   // padded quadrature point (2D: z = 0)
   for (int fi = 0; fi < nfaces; fi++)
   {
      FaceElementTransformations *Tr =
         (fi < fint.Size()) ? mesh.GetInteriorFaceTransformations(fint[fi])
                            : mesh.GetSharedFaceTransformations(fsh[fi - fint.Size()]);
      MFEM_VERIFY(Tr != nullptr, "ProjectStressToFaultKernel: no transformation for "
                  "local fault face " << fi);
      const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 2 * order + 1);
      const int nq = ir.GetNPoints();
      MFEM_VERIFY(e_q.Width() == nq && e_q.Height() == nbf,
                  "ProjectStressToFaultKernel: face basis at quadrature points is "
                  << e_q.Height() << " x " << e_q.Width() << ", expected "
                  << nbf << " x " << nq);
      Tq.SetSize(dim * nq);
      nlq.SetSize(nq);
      for (int q = 0; q < nq; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetAllIntPoints(&ip);
         CalcOrtho(Tr->Jacobian(), nor);
         const real_t nl = nor.Norml2();
         nlq(q) = nl;
         const real_t nq3[3] = {nor(0) / nl, nor(1) / nl,
                                dim == 3 ? nor(2) / nl : 0.0};
         Tr->Transform(ip, xq);
         for (int c = 0; c < 3; c++) { x3[c] = (c < xq.Size()) ? xq(c) : 0.0; }
         real_t s[9];
         sig(x3, s);
         if (dim == 2)
         {
            // in-plane: the anti-plane traction sigma_{zx} n_x + sigma_{zy} n_y
            // has no unknown to act on (file comment)
            MFEM_VERIFY(s[2] == 0.0 && s[5] == 0.0 && s[6] == 0.0 && s[7] == 0.0,
                        "ProjectStressToFaultKernel: 2D mesh but the stress at fault point ("
                        << x3[0] << ", " << x3[1] << ") has sigma_xz = " << s[2]
                        << ", sigma_yz = " << s[5] << " (sigma_zx = " << s[6]
                        << ", sigma_zy = " << s[7] << "); the in-plane model "
                        "carries no anti-plane stress");
         }
         for (int a = 0; a < dim; a++)
         {
            Tq(a * nq + q) = s[3 * a] * nq3[0] + s[3 * a + 1] * nq3[1]
                             + s[3 * a + 2] * nq3[2];
         }
      }
      const FaultBasisData &b = fb->GetBasis(fi);
      real_t bv[3][3];
      for (int c = 0; c < 3; c++)
      {
         bv[0][c] = b.normal[c];
         bv[1][c] = b.tangent1[c];
         bv[2][c] = b.tangent2[c];
      }
      DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
         dim, 3, Tq, nlq, ir, nbf, e_q, bv, out,
         b.qp_data.empty() ? nullptr : &b.qp_data);
      for (int kk = 0; kk < nbf; kk++)
      {
         const int k = fi * nbf + kk;
         K(3 * k) = out(nbf + kk);
         K(3 * k + 1) = out(2 * nbf + kk);
         // sigma_n = -(projected normal component), as ComputeTractionImpl.
         K(3 * k + 2) = -out(kk);
      }
   }
}

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_STRESS_PROJECTION_HPP
