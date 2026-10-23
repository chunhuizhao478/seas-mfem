// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_newmark.hpp — poroelastic option, step P2.5 (plan rev 3).
//
// PoroNewmarkStepper advances x = [u; w], v and a by one step of derivation
// rev 1.2, section 9: Newmark with beta = 0 and a velocity coefficient matrix
// Gamma (9.1), the acceleration solve (9.2), Gamma of (9.5),
//
//   1   x^{n+1} = x^n + dt v^n + dt^2/2 a^n                          (9.1)
//   2   v*      = v^n + dt (I - Gamma) a^n                           (9.1)
//   2b  V^{n+1} from the PoroSource, when one is set                 (6.6)
//   3   the fault advance at t_{n+1}: boxes 3-11 of the elastic operator
//       with u^{n+1} and v*_u (LtsAdvanceFault, the entry ImplicitSolve
//       uses; hook H2 as set on that operator)
//   4   F = PoroNitscheOperator::Force(x^{n+1}, v*, d^{n+1}) at t_{n+1}
//   5   a^{n+1} = PoroMass(Gamma).Apply(F)                           (9.2)
//   6   v^{n+1} = v* + dt Gamma a^{n+1}                              (9.1)
//
// (plan rev 3, section 5.2).  Gamma per element (9.5): the w row is
// [0, gamma_w], with gamma_w = FittedGamma(z), z = dt / tau (8.4), once per
// dt ("fitted"), or 1/2 ("half", test U2-15 only); the u row is [1/2, g_uw],
// g_uw = 0 ("diagonal") or (rho_f / rho)(1/2 - gamma_w) ("momentum", F3).
//
// The u rows are MFEM's NewmarkSolver::Step (linalg/ode.cpp:1271-1304) with
// beta = 0 and gamma = 1/2, which is the elastic driver's
// CentralDifferenceSolver, operation by operation: the same Vector::Add calls
// with the same factors (fac0 dt dt, fac2 dt, fac4 dt) in the same order, and
// MFEM's x.Add(beta dt^2, a) after the solve, which Vector::Add skips for a
// zero factor.  The g_uw terms of Gamma_mom follow element by element where
// g_uw != 0.  So with alpha = 0, rho_f = 0 and "half" the u block is today's
// elastic step bit for bit (test U2-15).  Kelvin-Voigt damping with Dirichlet
// velocity data (the quasi-dynamic loading, kv_vD_) is refused: its velocity
// time shift belongs to ImplicitSolve.
//
// The start.  The driver starts with a^0 = 0 (derivation section 9.4).
// ConsistentStart gives a^0 = M^{-1} F(x^0, v^0, d^0) instead, the first pass
// of NewmarkSolver::Step: the u rows of today's driver step from rest, and
// the start Proposition 9.1 assumes.  ApplyMassInverse is that M^{-1} (no
// drag, no dashpot).
//
// Collective: the constructor (PoroMass builds a ParFiniteElementSpace),
// ConsistentStart and Step (the operator's force; the fault advance).
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#ifndef MFEM_SEAS_PORO_NEWMARK_HPP
#define MFEM_SEAS_PORO_NEWMARK_HPP

#include "mfem.hpp"
#include "poro_mass.hpp"
#include "poro_nitsche_operator.hpp"
#include "poro_source.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace mfem
{
namespace seas
{
namespace poro
{

/// The fitted coefficient of derivation (9.4): gamma_w(z) = 1 / (1 - e^{-z}) -
/// 1 / z, with which one step of the scalar drag model decays by e^{-z}
/// exactly (Proposition 9.1).  For z < 1e-2 the series 1/2 + z/12 - z^3/720 +
/// z^5/30240 avoids the cancellation of the closed form; z = 0 gives 1/2.
inline real_t FittedGamma(real_t z)
{
   MFEM_VERIFY(z >= 0.0, "FittedGamma: z = " << z << " must be >= 0.");
   if (z < 1e-2) { return 0.5 + z / 12.0 - z * z * z / 720.0 + z * z * z * z * z / 30240.0; }
   return 1.0 / (-std::expm1(-z)) - 1.0 / z;
}

/// [multiphysics.poroelastic] drag_gamma.
enum class DragGamma
{
   kFitted,   ///< gamma_w = FittedGamma(dt / tau), "fitted"
   kHalf      ///< gamma_w = 1/2, "half" (test U2-15)
};

#ifdef MFEM_USE_MPI

class PoroNewmarkStepper
{
public:
   /// The stepper of `op` for the step dt, which must be the elastic
   /// operator's (op.Elastic().SetTimeStep(dt)).  `src`, when given, is
   /// advanced in line 2b; it must be the one set on `op` (SetSource).
   PoroNewmarkStepper(PoroNitscheOperator &op, real_t dt, DragGamma drag_gamma,
                      VelocityCoupling coupling, PoroSource *src = nullptr)
      : op_(op), dt_(dt), drag_gamma_(drag_gamma), coupling_(coupling), src_(src)
   {
      NitscheWaveOperator &el = op.Elastic();
      MFEM_VERIFY(dt > 0.0, "PoroNewmarkStepper: dt (" << dt << ") must be > 0.");
      MFEM_VERIFY(el.GetTimeStep() == dt, "PoroNewmarkStepper: dt = " << dt
                  << " but the elastic operator's step is " << el.GetTimeStep()
                  << "; call op.Elastic().SetTimeStep(dt) first.");
      c_ = ElementConstants(op.Sample());
      const int ne = static_cast<int>(c_.size());
      z_.SetSize(ne);
      gamma_w_.SetSize(ne);
      guw_dt_.SetSize(ne);
      cw_pred_.SetSize(ne);
      cw_corr_.SetSize(ne);
      for (int e = 0; e < ne; e++)
      {
         const BiotConstants &k = c_[e];
         // z = dt / tau with 1 / tau = b rho / (rho m - rho_f^2), eq. (8.4);
         // b = 0 gives z = 0
         z_(e) = dt * k.b * k.rho / (k.rho * k.m - k.rho_f * k.rho_f);
         gamma_w_(e) = drag_gamma == DragGamma::kHalf ? 0.5 : FittedGamma(z_(e));
         guw_dt_(e) = GammaUW(coupling, k.rho, k.rho_f, gamma_w_(e)) * dt;
         cw_pred_(e) = (1.0 - gamma_w_(e)) * dt;
         cw_corr_(e) = gamma_w_(e) * dt;
      }
      mass_ = std::make_unique<PoroMass>(op.GetMesh(), el.GetOrder(), c_, dt, gamma_w_,
                                         coupling, el.HasDamping() ? &el.DampingBlocks()
                                         : nullptr);
      n_ = op.NumU();
      F_.SetSize(2 * n_);
   }

   PoroNewmarkStepper(const PoroNewmarkStepper &) = delete;
   PoroNewmarkStepper &operator=(const PoroNewmarkStepper &) = delete;

   /// One step from t to t + dt (lines 1 to 6 of the header); x, v, a hold
   /// [u; w] (size 2 n) and are advanced in place; t += dt.
   void Step(Vector &x, Vector &v, Vector &a, real_t &t, real_t dt)
   {
      MFEM_VERIFY(dt == dt_, "PoroNewmarkStepper::Step: dt = " << dt
                  << " but the stepper was built for " << dt_ << ".");
      MFEM_VERIFY(x.Size() == 2 * n_ && v.Size() == 2 * n_ && a.Size() == 2 * n_,
                  "PoroNewmarkStepper::Step: sizes " << x.Size() << ", " << v.Size()
                  << ", " << a.Size() << " != 2 n = " << 2 * n_);
      NitscheWaveOperator &el = op_.Elastic();
      MFEM_VERIFY(!el.KelvinVoigtHasDirichletVelocity(), "PoroNewmarkStepper: "
                  "Kelvin-Voigt damping with Dirichlet velocity data is not wired into "
                  "the poroelastic step (its time shift belongs to ImplicitSolve).");
      op_.SetTime(t + dt);
      Predict(x, v, a);
      if (src_) { src_->Advance(t + dt); }
      if (el.HasFault())
      {
         Vector xu, vu;
         xu.MakeRef(x, 0, n_);
         vu.MakeRef(v, 0, n_);
         el.LtsAdvanceFault(xu, t + dt, &vu);
      }
      op_.Force(x, &v, el.HasFault() ? &el.Slip() : nullptr, F_);
      Solve(F_, a);
      Correct(x, v, a);
      t += dt;
   }

   /// a = M^{-1} F(x, v, d) at op.GetTime() with the stored slip d: the first
   /// pass of NewmarkSolver::Step (header).
   void ConsistentStart(const Vector &x, const Vector &v, Vector &a)
   {
      NitscheWaveOperator &el = op_.Elastic();
      op_.Force(x, &v, el.HasFault() ? &el.Slip() : nullptr, F_);
      ApplyMassInverse(F_, a);
   }

   /// Lines 1 and 2: x += dt v + dt^2/2 a; v becomes v* = v + dt (I - Gamma) a.
   void Predict(Vector &x, Vector &v, const Vector &a) const
   {
      const real_t dt = dt_;
      const real_t fac0 = 0.5 - kBeta, fac2 = 1.0 - kGammaU;
      x.Add(dt, v);
      x.Add(fac0*dt*dt, a);
      Vector vu, au;
      vu.MakeRef(v, 0, n_);
      au.MakeRef(const_cast<Vector &>(a), 0, n_);
      vu.Add(fac2*dt, au);
      AddFluid_(v, a, cw_pred_, -1.0);
   }

   /// Line 5: a = (M + dt D Gamma)^{-1} F.
   void Solve(const Vector &F, Vector &a) const { mass_->Apply(F, a); }

   /// After the solve: MFEM's x.Add(beta dt^2, a) (a no-op, beta = 0) and line
   /// 6, v = v* + dt Gamma a.
   void Correct(Vector &x, Vector &v, const Vector &a) const
   {
      const real_t dt = dt_;
      const real_t fac3 = kBeta, fac4 = kGammaU;
      x.Add(fac3*dt*dt, a);
      Vector vu, au;
      vu.MakeRef(v, 0, n_);
      au.MakeRef(const_cast<Vector &>(a), 0, n_);
      vu.Add(fac4*dt, au);
      AddFluid_(v, a, cw_corr_, 1.0);
   }

   /// a = M^{-1} F: the mass alone, without drag or dashpot.
   void ApplyMassInverse(const Vector &F, Vector &a)
   {
      if (!mass0_)
      {
         std::vector<BiotConstants> c0 = c_;
         for (BiotConstants &k : c0) { k.b = 0.0; }
         mass0_ = std::make_unique<PoroMass>(op_.GetMesh(), op_.Elastic().GetOrder(), c0,
                                             dt_, gamma_w_, coupling_, nullptr);
      }
      mass0_->Apply(F, a);
   }

   real_t Dt() const { return dt_; }
   DragGamma Drag() const { return drag_gamma_; }
   VelocityCoupling Coupling() const { return coupling_; }
   /// z = dt / tau and gamma_w per local element.
   const Vector &Z() const { return z_; }
   const Vector &GammaW() const { return gamma_w_; }
   /// g_uw dt per local element (zero for "diagonal" and where rho_f = 0).
   const Vector &GuwDt() const { return guw_dt_; }
   const PoroMass &Mass() const { return *mass_; }
   const std::vector<BiotConstants> &Constants() const { return c_; }

private:
   static constexpr real_t kBeta = 0.0;     ///< Newmark beta
   static constexpr real_t kGammaU = 0.5;   ///< the (u, u) entry of Gamma

   /// The terms of Gamma off the u-u entry, element by element:
   ///   v_w += c_e a_w                  (c_e = (1 - gamma_w) dt, or gamma_w dt)
   ///   v_u += s g_uw dt a_w            where g_uw != 0
   /// with s = -1 for (I - Gamma) and +1 for Gamma.
   void AddFluid_(Vector &v, const Vector &a, const Vector &cw, real_t s) const
   {
      const NitscheMass &mass = mass_->SchurMass();
      const int ns = mass.NumScalarDofs();
      for (int e = 0; e < cw.Size(); e++)
      {
         const Array<int> &dofs = mass.ElementScalarDofs(e);
         const real_t c = cw(e), g = s * guw_dt_(e);
         for (int comp = 0; comp < 3; comp++)
         {
            const int off = comp * ns;
            for (int k = 0; k < dofs.Size(); k++)
            {
               const int i = off + dofs[k];
               v(n_ + i) += c * a(n_ + i);
               if (g != 0.0) { v(i) += g * a(n_ + i); }
            }
         }
      }
   }

   PoroNitscheOperator &op_;
   real_t dt_;
   DragGamma drag_gamma_;
   VelocityCoupling coupling_;
   PoroSource *src_ = nullptr;
   std::vector<BiotConstants> c_;
   Vector z_, gamma_w_, guw_dt_, cw_pred_, cw_corr_;
   std::unique_ptr<PoroMass> mass_, mass0_;
   int n_ = 0;
   Vector F_;
};

#endif // MFEM_USE_MPI

} // namespace poro
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PORO_NEWMARK_HPP
