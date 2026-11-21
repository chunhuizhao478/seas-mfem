// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// split_force_kernels.hpp — internal force of the split elasticity for the
// explicit scheme, as one element loop (mfem::forall) over E-vectors.
//
// MOOSE: ADStressDivergenceTensors / ADDynamicStressDivergenceTensors residual
// sigma.row(c) . grad(test), sigma = g(d) sigma+ + sigma- from
// SmallDeformationIsotropicElasticity (RACCOON), quadrature order 2 (MOOSE's
// default SECOND for FIRST Lagrange: 3 points on triangles, 4 on tetrahedra).
//
// Linear simplices only (P1 on triangles and tetrahedra): the strain is
// constant per element and sigma is affine in g, so
//
//   sum_q w_q sigma(eps_e, g(d(x_q))) = |T| sigma(eps_e, gbar_e),
//   gbar_e = sum_q w_q g(d(x_q)) / |T|,
//
// exactly; the same holds for psi = g psi+ + psi-.  One split per element.
//
// Data: B(a, j, e) = dN_a/dx_j (constant per element), V(e) = |T|, and the
// reference shape values N(a, q) with the weights w(q) / sum w for gbar.
// Layout of E-vectors: MFEM ElementRestriction, (dof, vdim, element).

#ifndef MFEM_PFF_SPLIT_FORCE_KERNELS_HPP
#define MFEM_PFF_SPLIT_FORCE_KERNELS_HPP

#include "mfem.hpp"
#include "../constitutive/strain_split.hpp"
#include "../phasefield/phase_field_functions.hpp"

namespace mfem
{
namespace pff
{

class SplitForceKernel
{
public:
   /// vfes: H1 order 1 vector space (byNODES); dfes: H1 order 1 scalar space
   /// of the damage on the same mesh.
   SplitForceKernel(ParFiniteElementSpace &vfes, ParFiniteElementSpace &dfes,
                    real_t K, real_t G, Decomposition dec, const PowerDegradation &deg,
                    int quadrature_order = 2)
      : vfes_(vfes), dfes_(dfes), K_(K), G_(G), dec_(dec), deg_(deg)
   {
      ParMesh &mesh = *vfes.GetParMesh();
      dim_ = mesh.Dimension();
      ne_ = mesh.GetNE();
      MFEM_VERIFY(vfes.GetVDim() == dim_ && vfes.GetOrdering() == Ordering::byNODES,
                  "SplitForceKernel: the displacement space must be vector H1, byNODES");
      const Geometry::Type geom = mesh.GetTypicalElementGeometry();
      MFEM_VERIFY(geom == Geometry::TRIANGLE || geom == Geometry::TETRAHEDRON,
                  "SplitForceKernel: only linear simplices (triangles, tetrahedra)");
      MFEM_VERIFY(ne_ == 0 || (vfes.GetFE(0)->GetOrder() == 1 && dfes.GetFE(0)->GetOrder() == 1),
                  "SplitForceKernel: only P1 displacement and damage");
      for (int e = 0; e < ne_; e++)
      {
         MFEM_VERIFY(mesh.GetElementGeometry(e) == geom, "SplitForceKernel: mixed element types");
      }
      nd_ = dim_ + 1;
      const IntegrationRule &ir = IntRules.Get(geom, quadrature_order);
      nq_ = ir.GetNPoints();

      B_.SetSize(nd_ * dim_ * ne_);
      V_.SetSize(ne_);
      Nq_.SetSize(nd_ * nq_);
      Wq_.SetSize(nq_);
      {
         const FiniteElement &fe = *dfes.GetFE(0);
         Vector shape(nd_);
         real_t wsum = 0.0;
         for (int q = 0; q < nq_; q++) { wsum += ir.IntPoint(q).weight; }
         for (int q = 0; q < nq_; q++)
         {
            fe.CalcShape(ir.IntPoint(q), shape);
            for (int a = 0; a < nd_; a++) { Nq_(a + nd_ * q) = shape(a); }
            Wq_(q) = ir.IntPoint(q).weight / wsum;
         }
      }
      DenseMatrix dshape(nd_, dim_);
      IntegrationPoint c;
      c.Init(0);
      Geometries.GetCenter(geom).Get(&c.x, dim_);
      auto Bh = Reshape(B_.HostWrite(), nd_, dim_, ne_);
      auto Vh = V_.HostWrite();
      for (int e = 0; e < ne_; e++)
      {
         ElementTransformation &T = *vfes.GetElementTransformation(e);
         T.SetIntPoint(&c);
         vfes.GetFE(e)->CalcPhysDShape(T, dshape);
         for (int j = 0; j < dim_; j++)
            for (int a = 0; a < nd_; a++) { Bh(a, j, e) = dshape(a, j); }
         Vh[e] = T.Weight() * Geometry::Volume[geom];
      }

      Ru_ = vfes.GetElementRestriction(ElementDofOrdering::NATIVE);
      Rd_ = dfes.GetElementRestriction(ElementDofOrdering::NATIVE);
      ue_.SetSize(Ru_->Height());
      fe_.SetSize(Ru_->Height());
      de_.SetSize(Rd_->Height());
      uL_.SetSize(vfes.GetVSize());
      fL_.SetSize(vfes.GetVSize());
      dL_.SetSize(dfes.GetVSize());
      psi_e_.SetSize(ne_);
      ue_.UseDevice(true); fe_.UseDevice(true); de_.UseDevice(true);
      uL_.UseDevice(true); fL_.UseDevice(true); dL_.UseDevice(true);
      psi_e_.UseDevice(true);
   }

   /// F = P^T R^T f_e(R P U, R P D): the assembled internal force on true dofs.
   /// When psi_total is given, it receives the global int (g psi+ + psi-).
   void Mult(const Vector &U, const Vector &D, Vector &F, real_t *psi_total = nullptr) const
   {
      vfes_.GetProlongationMatrix()->Mult(U, uL_);
      dfes_.GetProlongationMatrix()->Mult(D, dL_);
      Ru_->Mult(uL_, ue_);
      Rd_->Mult(dL_, de_);
      ElementForces(ue_, de_, fe_, psi_e_);
      Ru_->MultTranspose(fe_, fL_);
      F.SetSize(vfes_.GetTrueVSize());
      vfes_.GetProlongationMatrix()->MultTranspose(fL_, F);
      if (psi_total)
      {
         real_t loc = psi_e_.Sum(), glob = 0.0;
         MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM,
                       vfes_.GetComm());
         *psi_total = glob;
      }
   }

   /// Element values of psi_e^+ (undegraded, constant per P1 element): the
   /// CONSTANT MONOMIAL psie_active that MOOSE transfers to the fracture app.
   /// psi_e has one entry per local element.
   void ElementActiveEnergy(const Vector &U, Vector &psi_e) const
   {
      vfes_.GetProlongationMatrix()->Mult(U, uL_);
      Ru_->Mult(uL_, ue_);
      dL_ = 0.0;
      Rd_->Mult(dL_, de_);
      ElementForces(ue_, de_, fe_, psi_e_, &psi_e);
   }

   int NumQuadraturePoints() const { return nq_; }

private:
   void ElementForces(const Vector &ue, const Vector &de, Vector &fe, Vector &psi_e,
                      Vector *psi_active = nullptr) const
   {
      const bool want_active = psi_active != nullptr;
      if (want_active) { psi_active->SetSize(ne_); psi_active->UseDevice(true); }
      real_t *PA = want_active ? psi_active->Write() : nullptr;
      const int ND = nd_, DIM = dim_, NE = ne_;
      const real_t K = K_, G = G_;
      const Decomposition dec = dec_;
      const PowerDegradation deg = deg_;
      const auto B = Reshape(B_.Read(), ND, DIM, NE);
      const auto Vol = V_.Read();
      const auto N = Reshape(Nq_.Read(), ND, nq_);
      const auto Wq = Wq_.Read();
      const int NQ = nq_;
      const auto U = Reshape(ue.Read(), ND, DIM, NE);
      const auto Dv = Reshape(de.Read(), ND, NE);
      auto Fo = Reshape(fe.Write(), ND, DIM, NE);
      auto P = psi_e.Write();
      mfem::forall(NE, [=] MFEM_HOST_DEVICE (int e)
      {
         real_t eps[9];
         for (int i = 0; i < 9; i++) { eps[i] = 0.0; }
         for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
            {
               real_t gij = 0.0, gji = 0.0;
               for (int a = 0; a < ND; a++)
               {
                  gij += U(a, i, e) * B(a, j, e);
                  gji += U(a, j, e) * B(a, i, e);
               }
               eps[i + 3 * j] = 0.5 * (gij + gji);
            }
         real_t gbar = 0.0;
         for (int q = 0; q < NQ; q++)
         {
            real_t dq = 0.0;
            for (int a = 0; a < ND; a++) { dq += N(a, q) * Dv(a, e); }
            gbar += Wq[q] * deg.g(dq);
         }
         SplitResult r;
         ComputeSplit(dec, eps, K, G, gbar, r);
         const real_t v = Vol[e];
         for (int a = 0; a < ND; a++)
            for (int i = 0; i < DIM; i++)
            {
               real_t s = 0.0;
               for (int j = 0; j < DIM; j++) { s += r.sigma[i + 3 * j] * B(a, j, e); }
               Fo(a, i, e) = v * s;
            }
         P[e] = v * r.psi;
         if (want_active) { PA[e] = r.psi_active; }
      });
   }

   ParFiniteElementSpace &vfes_, &dfes_;
   real_t K_, G_;
   Decomposition dec_;
   PowerDegradation deg_;
   int dim_ = 0, ne_ = 0, nd_ = 0, nq_ = 0;
   Vector B_, V_, Nq_, Wq_;
   const Operator *Ru_ = nullptr, *Rd_ = nullptr;
   mutable Vector ue_, fe_, de_, uL_, fL_, dL_, psi_e_;
};

} // namespace pff
} // namespace mfem

#endif
