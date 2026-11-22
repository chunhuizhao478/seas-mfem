// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// explicit_elasticity_operator.hpp — the pieces of one explicit displacement
// step: lumped mass, lumped dashpot blocks, the pulse load vector and the
// split internal force (kernels/split_force_kernels.hpp).
//
// MOOSE objects (farms@cdms pulsepower/spe_hftc/code_{2d,3d}/elasticity_absorbing.i):
//   ADInertialForce            consistent mass rho (Newmark); pff lumps it (row sum).
//   ADDynamicStressDivergenceTensors (alpha = 0)
//                              sigma.grad(test) (1 + zeta/dt) - zeta/dt sigma_old.grad(test),
//                              i.e. F_int + zeta (F_int - F_int_old)/dt with the stored
//                              (degraded) stress of the previous step.
//   [BCs]/[Pressure] (ADPressure, use_displaced_mesh = false)
//                              residual p(t, x) n_c test, n the outward normal, so the
//                              load on the solid is -p n: b = -int Peak(x) n phi, times Shape(t).
//   FarmsNonReflectDashpotBC   int test rho (c_p v_n n_c + c_s (v_c - v_n n_c)) ds; pff
//                              lumps it per node (row sum over each face, normals of
//                              each face kept), a dim x dim block per node.
//
// All vectors here are true-dof vectors. The vector space is byNODES, so true dof
// (c, i) of the displacement is c * nt + i, with i the true dof of the scalar space.

#ifndef MFEM_PFF_EXPLICIT_ELASTICITY_OPERATOR_HPP
#define MFEM_PFF_EXPLICIT_ELASTICITY_OPERATOR_HPP

#include "mfem.hpp"
#include "../kernels/split_force_kernels.hpp"
#include "../loading/pulse_load.hpp"

#include <cmath>
#include <vector>

namespace mfem
{
namespace pff
{

struct ExplicitElasticityParams
{
   real_t K = 0.0, G = 0.0, density = 0.0;
   Decomposition decomposition = Decomposition::Spectral;
   PowerDegradation deg;
   Array<int> pressure_marker;    ///< boundary attributes with the pulse pressure
   Array<int> absorbing_marker;   ///< boundary attributes with Lysmer dashpots
   int quadrature_order = 2;
   real_t Cp() const { return std::sqrt((K + 4.0 / 3.0 * G) / density); }
   real_t Cs() const { return std::sqrt(G / density); }
};

class ExplicitElasticityOperator
{
public:
   ExplicitElasticityOperator(ParFiniteElementSpace &vfes, ParFiniteElementSpace &sfes,
                              const ExplicitElasticityParams &p, const PulseLoad &pulse,
                              const Array<int> &fixed_tdofs)
      : vfes_(vfes), sfes_(sfes), p_(p), pulse_(pulse), fixed_(fixed_tdofs),
        force_(vfes, sfes, p.K, p.G, p.decomposition, p.deg, p.quadrature_order)
   {
      dim_ = vfes.GetParMesh()->Dimension();
      nt_ = sfes.GetTrueVSize();
      MFEM_VERIFY(vfes.GetTrueVSize() == dim_ * nt_, "vector/scalar true dofs do not match");
      AssembleLumpedMass();
      AssembleDashpots();
      AssemblePressure();
   }

   int Dim() const { return dim_; }
   int ScalarTrueSize() const { return nt_; }
   const Vector &LumpedMass() const { return M_; }     ///< per scalar true dof
   const Vector &PressureVector() const { return b_; } ///< -int Peak n phi (true dofs)
   const Array<int> &FixedTrueDofs() const { return fixed_; }
   const SplitForceKernel &Force() const { return force_; }

   /// F_int(U; D) and optionally int psi.
   void InternalForce(const Vector &U, const Vector &D, Vector &F, real_t *psi = nullptr) const
   {
      force_.Mult(U, D, F, psi);
   }

   /// f = C v (lumped dashpot force).
   void DashpotForce(const Vector &v, Vector &f) const
   {
      f.SetSize(v.Size());
      f = 0.0;
      const int D = dim_;
      for (int i = 0; i < nt_; i++)
         for (int r = 0; r < D; r++)
         {
            real_t s = 0.0;
            for (int c = 0; c < D; c++) { s += Cb_(Block(i, r, c)) * v(c * nt_ + i); }
            f(r * nt_ + i) = s;
         }
   }

   /// f += alpha M v with the lumped mass.
   void AddMassForce(real_t alpha, const Vector &v, Vector &f) const
   {
      for (int c = 0; c < dim_; c++)
         for (int i = 0; i < nt_; i++) { f(c * nt_ + i) += alpha * M_(i) * v(c * nt_ + i); }
   }

   /// a = ((1 + h alpha) M + h C)^{-1} rhs per node (dim x dim solve where the
   /// node has a dashpot block), a = 0 on the fixed dofs.
   void SolveAcceleration(const Vector &rhs, real_t h, Vector &a, real_t alpha = 0.0) const
   {
      a.SetSize(rhs.Size());
      const int D = dim_;
      for (int i = 0; i < nt_; i++)
      {
         const real_t m = (1.0 + h * alpha) * M_(i);
         if (!has_dashpot_[i])
         {
            for (int c = 0; c < D; c++) { a(c * nt_ + i) = rhs(c * nt_ + i) / m; }
            continue;
         }
         DenseMatrix A(D);
         Vector r(D), x(D);
         for (int rr = 0; rr < D; rr++)
         {
            r(rr) = rhs(rr * nt_ + i);
            for (int c = 0; c < D; c++) { A(rr, c) = h * Cb_(Block(i, rr, c)) + (rr == c ? m : 0.0); }
         }
         A.Invert();
         A.Mult(r, x);
         for (int c = 0; c < D; c++) { a(c * nt_ + i) = x(c); }
      }
      for (int t : fixed_) { a(t) = 0.0; }
   }

   /// 1/2 v^T M v (lumped).
   real_t KineticEnergy(const Vector &v) const
   {
      real_t s = 0.0;
      for (int c = 0; c < dim_; c++)
         for (int i = 0; i < nt_; i++) { s += 0.5 * M_(i) * v(c * nt_ + i) * v(c * nt_ + i); }
      real_t g = 0.0;
      MPI_Allreduce(&s, &g, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes_.GetComm());
      return g;
   }

   /// max over nodes of |v|.
   real_t MaxNodalSpeed(const Vector &v) const
   {
      real_t m = 0.0;
      for (int i = 0; i < nt_; i++)
      {
         real_t s = 0.0;
         for (int c = 0; c < dim_; c++) { s += v(c * nt_ + i) * v(c * nt_ + i); }
         m = std::max(m, s);
      }
      real_t g = 0.0;
      MPI_Allreduce(&m, &g, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, vfes_.GetComm());
      return std::sqrt(g);
   }

   /// omega_max^2 of M^{-1} K (K the undamaged stiffness, D = 0) by power
   /// iteration; the fixed dofs are held at 0.
   real_t MaxEigenvalue(int iterations) const
   {
      Vector x(vfes_.GetTrueVSize()), y(x.Size()), D0(nt_);
      D0 = 0.0;
      // Start vector from the node coordinates: independent of the partition.
      {
         const int D = dim_;
         VectorFunctionCoefficient f(D, [D](const Vector &X, Vector &y)
         {
            real_t s = 0.0;
            for (int k = 0; k < X.Size(); k++) { s += (1234.5 + 1000.0 * k) * X(k); }
            for (int c = 0; c < D; c++) { y(c) = std::sin(s * (1.0 + 0.37 * c) + 0.3 + c); }
         });
         ParGridFunction xg(&vfes_);
         xg.ProjectCoefficient(f);
         xg.GetTrueDofs(x);
      }
      for (int t : fixed_) { x(t) = 0.0; }
      real_t lam = 0.0;
      for (int k = 0; k < iterations; k++)
      {
         // normalise in the M norm
         real_t xm = 0.0;
         for (int c = 0; c < dim_; c++)
            for (int i = 0; i < nt_; i++) { xm += M_(i) * x(c * nt_ + i) * x(c * nt_ + i); }
         real_t xmg = 0.0;
         MPI_Allreduce(&xm, &xmg, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes_.GetComm());
         x /= std::sqrt(xmg);
         force_.Mult(x, D0, y);
         lam = InnerProduct(vfes_.GetComm(), x, y);   // Rayleigh quotient, x^T M x = 1
         for (int c = 0; c < dim_; c++)
            for (int i = 0; i < nt_; i++) { x(c * nt_ + i) = y(c * nt_ + i) / M_(i); }
         for (int t : fixed_) { x(t) = 0.0; }
      }
      return lam;
   }

   /// Peak pressure on the pressure boundary (max over its nodes of Peak(x)).
   real_t PeakOnPressureBoundary() const { return peak_max_; }

private:
   int Block(int i, int r, int c) const { return (r * dim_ + c) * nt_ + i; }

   void AssembleLumpedMass()
   {
      ParLinearForm m(&sfes_);
      ConstantCoefficient rho(p_.density);
      m.AddDomainIntegrator(new DomainLFIntegrator(rho));
      m.Assemble();
      M_.SetSize(nt_);
      m.ParallelAssemble(M_);
      for (int i = 0; i < nt_; i++)
      {
         MFEM_VERIFY(M_(i) > 0.0, "lumped mass is not positive at true dof " << i);
      }
   }

   void AssembleDashpots()
   {
      ParMesh &mesh = *sfes_.GetParMesh();
      const int D = dim_;
      const real_t rho = p_.density, cp = p_.Cp(), cs = p_.Cs();
      std::vector<Vector> L(D * D);
      for (auto &v : L) { v.SetSize(sfes_.GetVSize()); v = 0.0; }
      Vector shape, nor(D);
      Array<int> dofs;
      for (int be = 0; be < mesh.GetNBE(); be++)
      {
         const int attr = mesh.GetBdrAttribute(be);
         if (p_.absorbing_marker.Size() == 0 || attr > p_.absorbing_marker.Size() ||
             !p_.absorbing_marker[attr - 1]) { continue; }
         const FiniteElement &fe = *sfes_.GetBE(be);
         ElementTransformation &T = *mesh.GetBdrElementTransformation(be);
         sfes_.GetBdrElementDofs(be, dofs);
         shape.SetSize(fe.GetDof());
         const IntegrationRule &ir = IntRules.Get(fe.GetGeomType(), 2);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            T.SetIntPoint(&ip);
            CalcOrtho(T.Jacobian(), nor);
            const real_t ds = nor.Norml2();
            fe.CalcShape(ip, shape);
            for (int r = 0; r < D; r++)
               for (int c = 0; c < D; c++)
               {
                  const real_t nn = nor(r) * nor(c) / (ds * ds);
                  const real_t blk = rho * (cp * nn + cs * ((r == c ? 1.0 : 0.0) - nn));
                  for (int a = 0; a < fe.GetDof(); a++)
                  {
                     L[r * D + c](dofs[a]) += ip.weight * ds * shape(a) * blk;
                  }
               }
         }
      }
      Cb_.SetSize(D * D * nt_);
      Vector tmp(nt_);
      for (int k = 0; k < D * D; k++)
      {
         sfes_.GetProlongationMatrix()->MultTranspose(L[k], tmp);
         for (int i = 0; i < nt_; i++) { Cb_(k * nt_ + i) = tmp(i); }
      }
      has_dashpot_.assign(nt_, false);
      for (int i = 0; i < nt_; i++)
         for (int k = 0; k < D * D; k++)
            if (Cb_(k * nt_ + i) != 0.0) { has_dashpot_[i] = true; }
   }

   void AssemblePressure()
   {
      ParMesh &mesh = *vfes_.GetParMesh();
      const int D = dim_;
      Vector bL(vfes_.GetVSize());
      bL = 0.0;
      Vector shape, nor(D);
      Array<int> vdofs;
      real_t peak = 0.0;
      for (int be = 0; be < mesh.GetNBE(); be++)
      {
         const int attr = mesh.GetBdrAttribute(be);
         if (p_.pressure_marker.Size() == 0 || attr > p_.pressure_marker.Size() ||
             !p_.pressure_marker[attr - 1]) { continue; }
         FaceElementTransformations *ft = mesh.GetBdrFaceTransformations(be);
         MFEM_VERIFY(ft, "pressure boundary element " << be << " has no face transformation");
         const FiniteElement &fe = *vfes_.GetBE(be);
         vfes_.GetBdrElementVDofs(be, vdofs);
         const int nd = fe.GetDof();
         shape.SetSize(nd);
         const IntegrationRule &ir = IntRules.Get(fe.GetGeomType(), 2);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ft->SetAllIntPoints(&ip);
            // outward normal of element 1 (the solid), |nor| = surface Jacobian
            CalcOrtho(ft->Jacobian(), nor);
            Vector x(3);
            x = 0.0;
            Vector xd(x.GetData(), D);
            ft->Transform(ip, xd);
            const real_t P = pulse_.Peak(x.GetData());
            peak = std::max(peak, P);
            fe.CalcShape(ip, shape);
            for (int c = 0; c < D; c++)
               for (int a = 0; a < nd; a++)
               {
                  bL(vdofs[a + c * nd]) -= ip.weight * P * nor(c) * shape(a);
               }
         }
      }
      b_.SetSize(vfes_.GetTrueVSize());
      vfes_.GetProlongationMatrix()->MultTranspose(bL, b_);
      MPI_Allreduce(&peak, &peak_max_, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                    vfes_.GetComm());
   }

   ParFiniteElementSpace &vfes_, &sfes_;
   ExplicitElasticityParams p_;
   const PulseLoad &pulse_;
   Array<int> fixed_;
   SplitForceKernel force_;
   int dim_ = 0, nt_ = 0;
   Vector M_, Cb_, b_;
   std::vector<bool> has_dashpot_;
   real_t peak_max_ = 0.0;
};

} // namespace pff
} // namespace mfem

#endif
