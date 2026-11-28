// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// pore_pressure_operator.hpp — the pore-pressure side of the u-p Biot model of
// cmame_revision hydromech (phase C4): one backward-Euler pressure step and the
// pore-pressure force on the skeleton, with the damage-dependent properties of
// poro_damage.hpp evaluated at the quadrature points from the nodal damage.
//
// MOOSE objects (farms@cdms pulsepower/cmame_revision/2d_hydromech/benchmark/
// elasticity_E1d25.i, src/kernels/twofieldHM):
//   ElkPorousFlowFullySaturatedMassTimeDerivative (HydroMechanical,
//   multiply_by_density = false, use_damaged_biot = true)
//        test [ (p - p_old)/(dt M) + alpha (eps_v - eps_v_old)/dt ]
//        M = PorousFlow_constant_biot_modulus_qp (ElkPorousFlowDamagedBiotModulus),
//        eps_v the trace of the total small strain (PorousFlowVolumetricStrain).
//   PorousFlowFullySaturatedDarcyBase (multiply_by_density = false, no gravity)
//        grad(test) . (k / mu) grad p, k = effective_perm
//        (ElkPorousFlowPermeabilityDamaged), relative permeability 1.
//   ElkPorousFlowEffectiveStressCoupling (use_damaged_biot = true)
//        residual of disp_c: -alpha p d(test)/dx_c, so the internal force of the
//        skeleton gets F_p = -int alpha(d) p div(phi) and M a = f_ext - F_int - F_p.
//   [boundaries] drained: Dirichlet p (the drained variants of the MOOSE inputs
//   put the pulse pressure on the borehole pp).
//
// One step (PLAN.md, single pass): with u^{n+1} from the explicit step,
//   [ int N (1/M)/dt N + int grad N . (k/mu) grad N ] p^{n+1}
//        = int N (1/M)/dt p^n - int N alpha (eps_v^{n+1} - eps_v^n)/dt,
// all properties at the damage passed in. SPD: PCG with BoomerAMG ("amg", the
// PLAN choice) or diagonal scaling ("jacobi"; the cmame steps are dominated by
// the storage term, diffusion number k M dt / (mu h^2) < 0.1).
//
// The quadrature is the rule of the force kernel (order 2: MOOSE's default for
// FIRST order variables).
//
// permeability_model = normal_strain (cmame 3d_hydromech/em_5J, Heider 2021 eqs.
// 46-48 as in NDSmallDeformationIsotropicElasticity::updatePermeabilityForCracking):
// the tensor k = k0 I + d^b (w^2/12)(I - n n) at each point, from grad d, the
// total small strain of u^{n+1} (MOOSE "mechanical_strain" of the current
// iterate; single pass here) and h_c = the element size. MOOSE takes h_c from
// ElementLengthAux method = max, libMesh Elem::hmax(): the largest distance
// between two vertices of the element. Not wired (MFEM_VERIFY in the
// constructor): the strain porosity model.

#ifndef MFEM_PFF_PORE_PRESSURE_OPERATOR_HPP
#define MFEM_PFF_PORE_PRESSURE_OPERATOR_HPP

#include "mfem.hpp"
#include "poro_damage.hpp"

#include <memory>

namespace mfem
{
namespace pff
{

enum class PoroPreconditioner { AMG, Jacobi };

/// Properties at one damage value. mobility is the isotropic value; for the
/// normal_strain model it is the matrix value k0 / mu (zero strain, no normal),
/// and the tensor comes from PorePressureOperator::MobilityAt.
struct PoroPoint
{
   real_t alpha, phi, inv_M, mobility;   ///< mobility = k / mu
};

inline PoroPoint EvalPoroPoint(const PoroDamageParams &P, real_t d)
{
   PoroPoint q;
   q.alpha = BiotCoefficient(P, d);
   q.phi = Porosity(P, d);
   q.inv_M = 1.0 / BiotModulus(P, q.alpha, q.phi);
   real_t k[9];
   const real_t zero[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
   EffectivePermeability(P, d, zero, zero, 0.0, k);
   q.mobility = k[0] / P.viscosity;
   return q;
}

/// Largest (g(d) Hd + alpha(d)^2 M(d)) / Hd over d in [0, 1], Hd = lambda + 2G/dim:
/// a bound of omega_max^2(undrained) / omega_max^2(drained) for the CFL step
/// (eps : C : eps >= Hd tr(eps)^2, and the nodal p of the storage-only step is
/// an L2 projection of -alpha M eps_v). g2 = (1-d)^2 as in the Biot coefficient.
inline real_t UndrainedStiffnessRatio(const PoroDamageParams &P, real_t G, int dim)
{
   const real_t Hd = P.K - 2.0 * G / 3.0 + 2.0 * G / dim;
   real_t r = 1.0;
   for (int i = 0; i <= 1000; i++)
   {
      const real_t d = i / 1000.0;
      const PoroPoint q = EvalPoroPoint(P, d);
      const real_t g = (1.0 - d) * (1.0 - d);
      r = std::max(r, (g * Hd + q.alpha * q.alpha / q.inv_M) / Hd);
   }
   return r;
}

class PorePressureOperator
{
public:
   /// drained_tdofs: true dofs of the scalar space with prescribed p.
   PorePressureOperator(ParFiniteElementSpace &vfes, ParFiniteElementSpace &sfes,
                        const PoroDamageParams &P, const Array<int> &drained_tdofs,
                        const IntegrationRule &ir, PoroPreconditioner pc = PoroPreconditioner::AMG)
      : vfes_(vfes), sfes_(sfes), P_(P), ir_(ir), pc_(pc), d_gf_(&sfes), p_gf_(&sfes),
        u_gf_(&vfes), uold_gf_(&vfes)
   {
      MFEM_VERIFY(P.porosity_model == PorosityModel::Damage,
                  "PorePressureOperator: porosity_model = 'strain' is not wired yet");
      MFEM_VERIFY(P.K > 0.0 && P.K_s > 0.0 && P.K_f > 0.0 && P.viscosity > 0.0,
                  "PorePressureOperator: K, K_s, K_f and viscosity must be positive");
      drained_.Append(drained_tdofs);
      dim_ = vfes.GetParMesh()->Dimension();
      tensor_perm_ = P.perm_model == PermeabilityModel::NormalStrain;
      u_gf_ = 0.0;
      uold_gf_ = 0.0;
      if (tensor_perm_)
      {
         // libMesh Elem::hmax(): largest vertex-to-vertex distance.
         ParMesh &mesh = *vfes.GetParMesh();
         hmax_.SetSize(mesh.GetNE());
         Array<int> v;
         for (int e = 0; e < mesh.GetNE(); e++)
         {
            mesh.GetElementVertices(e, v);
            real_t h2 = 0.0;
            for (int a = 0; a < v.Size(); a++)
               for (int b = a + 1; b < v.Size(); b++)
               {
                  const real_t *xa = mesh.GetVertex(v[a]), *xb = mesh.GetVertex(v[b]);
                  real_t s = 0.0;
                  for (int c = 0; c < dim_; c++) { s += (xa[c] - xb[c]) * (xa[c] - xb[c]); }
                  h2 = std::max(h2, s);
               }
            hmax_(e) = std::sqrt(h2);
         }
      }
   }

   /// Element size h_c of the normal_strain model (empty for the other models).
   const Vector &ElementHmax() const { return hmax_; }

   /// Set the damage and displacement that MobilityAt and Energies read (Step
   /// does this itself with u^{n+1}).
   void SetState(const Vector &D, const Vector &U) const
   {
      d_gf_.SetFromTrueDofs(D);
      u_gf_.SetFromTrueDofs(U);
   }

   /// k / mu (dim x dim) at the point T is set to, from the damage and the
   /// displacement last passed to Step (zero before the first step).
   void MobilityAt(ElementTransformation &T, const IntegrationPoint &ip, DenseMatrix &K) const
   {
      K.SetSize(dim_);
      real_t k[9], gd[3] = {0.0, 0.0, 0.0}, eps[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      real_t h = 0.0;
      if (tensor_perm_)
      {
         Vector g(gd, dim_);
         d_gf_.GetGradient(T, g);
         DenseMatrix gu;
         u_gf_.GetVectorGradient(T, gu);
         for (int j = 0; j < dim_; j++)
            for (int i = 0; i < dim_; i++) { eps[i + 3 * j] = 0.5 * (gu(i, j) + gu(j, i)); }
         h = hmax_(T.ElementNo);
      }
      EffectivePermeability(P_, d_gf_.GetValue(T, ip), gd, eps, h, k);
      for (int j = 0; j < dim_; j++)
         for (int i = 0; i < dim_; i++) { K(i, j) = k[i + 3 * j] / P_.viscosity; }
   }

   const Array<int> &DrainedTrueDofs() const { return drained_; }

   /// F += -int alpha(d) p div(phi) (true dofs of the vector space).
   void AddCouplingForce(const Vector &D, const Vector &Pp, Vector &F) const
   {
      d_gf_.SetFromTrueDofs(D);
      p_gf_.SetFromTrueDofs(Pp);
      ParLinearForm lf(&vfes_);
      lf.AddDomainIntegrator(new CouplingIntegrator(*this));
      lf.Assemble();
      Vector t(vfes_.GetTrueVSize());
      lf.ParallelAssemble(t);
      F += t;
   }

   struct StepInfo
   {
      int iterations = 0;
      real_t fluid_energy = 0.0;      ///< 1/2 int p^2 / M at p^{n+1}
      real_t dissipation_rate = 0.0;  ///< int grad p . (k/mu) grad p at p^{n+1}
      real_t boundary_power = 0.0;    ///< sum over drained dofs of p_i r_i, r = A p - b:
                                      ///< the power the drained boundary puts into the fluid
   };

   /// Backward Euler step: Pp holds p^n on entry (with the drained values of
   /// t^{n+1} already written on the drained dofs) and p^{n+1} on return.
   void Step(const Vector &D, const Vector &U_new, const Vector &U_old, real_t dt, Vector &Pp,
             StepInfo &info)
   {
      MFEM_VERIFY(dt > 0.0, "PorePressureOperator::Step: dt must be positive");
      d_gf_.SetFromTrueDofs(D);
      u_gf_.SetFromTrueDofs(U_new);
      uold_gf_.SetFromTrueDofs(U_old);
      // p^n in the storage term: the incoming vector with the OLD values on the
      // drained dofs would be the exact backward difference there, but those
      // rows are eliminated, and the interior rows couple to the drained
      // columns only through the storage mass matrix; MOOSE uses pp_old.
      p_gf_.SetFromTrueDofs(p_old_valid_ ? p_old_ : Pp);
      dt_ = dt;

      a_.reset(new ParBilinearForm(&sfes_));
      StorageCoefficient storage(*this);
      MobilityCoefficient mobility(*this);
      MobilityTensorCoefficient mobility_tensor(*this);
      auto *mi = new MassIntegrator(storage);
      mi->SetIntRule(&ir_);
      a_->AddDomainIntegrator(mi);
      auto *di = tensor_perm_ ? new DiffusionIntegrator(mobility_tensor)
                 : new DiffusionIntegrator(mobility);
      di->SetIntRule(&ir_);
      a_->AddDomainIntegrator(di);
      a_->Assemble();
      a_->Finalize();

      ParLinearForm b(&sfes_);
      SourceCoefficient source(*this);
      auto *li = new DomainLFIntegrator(source);
      li->SetIntRule(&ir_);
      b.AddDomainIntegrator(li);
      b.Assemble();

      // Unconstrained operator and right-hand side, for the boundary reaction.
      std::unique_ptr<HypreParMatrix> A0(a_->ParallelAssemble());
      Vector B0(sfes_.GetTrueVSize());
      b.ParallelAssemble(B0);

      // Eliminated copy for the solve (drained rows and columns).
      Vector B(B0), X(Pp);
      std::unique_ptr<HypreParMatrix> Ac(new HypreParMatrix(*A0)), Ae;
      if (HasGlobalDrained())
      {
         Ae.reset(Ac->EliminateRowsCols(drained_));
         Ac->EliminateBC(*Ae, drained_, X, B);
      }

      CGSolver cg(sfes_.GetComm());
      cg.SetRelTol(1e-12);
      cg.SetAbsTol(0.0);
      cg.SetMaxIter(2000);
      cg.SetPrintLevel(0);
      std::unique_ptr<Solver> prec;
      if (pc_ == PoroPreconditioner::AMG)
      {
         auto *amg = new HypreBoomerAMG(*Ac);
         amg->SetPrintLevel(0);
         prec.reset(amg);
      }
      else
      {
         prec.reset(new HypreDiagScale(*Ac));
      }
      cg.SetPreconditioner(*prec);
      cg.SetOperator(*Ac);
      cg.Mult(B, X);
      MFEM_VERIFY(cg.GetConverged(), "pore-pressure PCG did not converge ("
                  << cg.GetNumIterations() << " iterations, |r| = " << cg.GetFinalNorm() << ")");
      info.iterations = cg.GetNumIterations();

      // Diagnostics at p^{n+1}.
      Vector r(X.Size());
      A0->Mult(X, r);
      r -= B0;
      real_t loc[1] = {0.0}, glob[1];
      for (int t : drained_) { loc[0] += X(t) * r(t); }
      MPI_Allreduce(loc, glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, sfes_.GetComm());
      info.boundary_power = glob[0];
      Pp = X;
      p_old_ = X;
      p_old_valid_ = true;
      Energies(D, Pp, info.fluid_energy, info.dissipation_rate);
   }

   /// 1/2 int p^2/M and int grad p . (k/mu) grad p at damage D.
   void Energies(const Vector &D, const Vector &Pp, real_t &fluid_energy, real_t &diss_rate) const
   {
      d_gf_.SetFromTrueDofs(D);
      ParGridFunction pg(&sfes_);
      pg.SetFromTrueDofs(Pp);
      ParMesh &mesh = *sfes_.GetParMesh();
      real_t loc[2] = {0.0, 0.0}, glob[2];
      Vector grad(dim_), Kg(dim_);
      DenseMatrix Kt;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         ElementTransformation &T = *mesh.GetElementTransformation(e);
         for (int q = 0; q < ir_.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir_.IntPoint(q);
            T.SetIntPoint(&ip);
            const PoroPoint pt = EvalPoroPoint(P_, d_gf_.GetValue(T, ip));
            const real_t p = pg.GetValue(T, ip);
            pg.GetGradient(T, grad);
            const real_t w = ip.weight * T.Weight();
            loc[0] += w * 0.5 * pt.inv_M * p * p;
            if (tensor_perm_)
            {
               MobilityAt(T, ip, Kt);
               Kt.Mult(grad, Kg);
               loc[1] += w * (Kg * grad);
            }
            else
            {
               loc[1] += w * pt.mobility * (grad * grad);
            }
         }
      }
      MPI_Allreduce(loc, glob, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, sfes_.GetComm());
      fluid_energy = glob[0];
      diss_rate = glob[1];
   }

   /// Forget the stored p^n (call after the state is changed from outside).
   void Reset() { p_old_valid_ = false; }

   real_t PressureAt(ElementTransformation &T, const IntegrationPoint &ip) const
   { return p_gf_.GetValue(T, ip); }
   real_t DamageAt(ElementTransformation &T, const IntegrationPoint &ip) const
   { return d_gf_.GetValue(T, ip); }
   const PoroDamageParams &Params() const { return P_; }

private:
   bool HasGlobalDrained() const
   {
      int loc = drained_.Size(), glob = 0;
      MPI_Allreduce(&loc, &glob, 1, MPI_INT, MPI_SUM, sfes_.GetComm());
      return glob > 0;
   }

   struct StorageCoefficient : public Coefficient
   {
      const PorePressureOperator &o;
      StorageCoefficient(const PorePressureOperator &o_) : o(o_) {}
      real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      { return EvalPoroPoint(o.P_, o.d_gf_.GetValue(T, ip)).inv_M / o.dt_; }
   };
   struct MobilityCoefficient : public Coefficient
   {
      const PorePressureOperator &o;
      MobilityCoefficient(const PorePressureOperator &o_) : o(o_) {}
      real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      { return EvalPoroPoint(o.P_, o.d_gf_.GetValue(T, ip)).mobility; }
   };
   struct MobilityTensorCoefficient : public MatrixCoefficient
   {
      const PorePressureOperator &o;
      MobilityTensorCoefficient(const PorePressureOperator &o_) : MatrixCoefficient(o_.dim_), o(o_) {}
      void Eval(DenseMatrix &K, ElementTransformation &T, const IntegrationPoint &ip) override
      { o.MobilityAt(T, ip, K); }
   };
   struct SourceCoefficient : public Coefficient
   {
      const PorePressureOperator &o;
      SourceCoefficient(const PorePressureOperator &o_) : o(o_) {}
      real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
         const PoroPoint pt = EvalPoroPoint(o.P_, o.d_gf_.GetValue(T, ip));
         const real_t dev = o.u_gf_.GetDivergence(T) - o.uold_gf_.GetDivergence(T);
         return (pt.inv_M * o.p_gf_.GetValue(T, ip) - pt.alpha * dev) / o.dt_;
      }
   };
   /// elvect(a + c nd) = -int alpha(d) p dN_a/dx_c.
   struct CouplingIntegrator : public LinearFormIntegrator
   {
      const PorePressureOperator &o;
      DenseMatrix dshape;
      CouplingIntegrator(const PorePressureOperator &o_) : o(o_) {}
      void AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &T,
                                  Vector &elvect) override
      {
         const int nd = el.GetDof(), dim = el.GetDim();
         dshape.SetSize(nd, dim);
         elvect.SetSize(nd * dim);
         elvect = 0.0;
         for (int q = 0; q < o.ir_.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = o.ir_.IntPoint(q);
            T.SetIntPoint(&ip);
            el.CalcPhysDShape(T, dshape);
            const real_t alpha = BiotCoefficient(o.P_, o.d_gf_.GetValue(T, ip));
            const real_t w = -ip.weight * T.Weight() * alpha * o.p_gf_.GetValue(T, ip);
            for (int c = 0; c < dim; c++)
               for (int a = 0; a < nd; a++) { elvect(a + c * nd) += w * dshape(a, c); }
         }
      }
   };

   ParFiniteElementSpace &vfes_, &sfes_;
   PoroDamageParams P_;
   const IntegrationRule &ir_;
   PoroPreconditioner pc_;
   Array<int> drained_;
   int dim_ = 0;
   bool tensor_perm_ = false;
   Vector hmax_;
   real_t dt_ = 1.0;
   mutable ParGridFunction d_gf_, p_gf_, u_gf_, uold_gf_;
   std::unique_ptr<ParBilinearForm> a_;
   Vector p_old_;
   bool p_old_valid_ = false;
};

} // namespace pff
} // namespace mfem

#endif
