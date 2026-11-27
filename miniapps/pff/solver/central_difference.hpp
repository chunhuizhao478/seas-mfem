// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// central_difference.hpp — one explicit step with lumped mass (half-step
// velocities, variable step; Belytschko, Liu & Moran, box 6.1):
//
//   a^n     = ((1 + h alpha) M + h C_b)^{-1}
//             [ s(t^n) b + b_c - F_int(u^n) - F_KV^n - (C_b + alpha M) v^{n-1/2} ],
//             h = dt_{n-1/2} / 2 (0 at the first step, where v holds v^0),
//             so the dashpot sees v^n = v^{n-1/2} + h a^n;
//   dt_avg = (dt_{n-1/2} + dt_{n+1/2}) / 2 (dt_{1/2} / 2 at the first step);
//   v^{n+1/2} = v^{n-1/2} + dt_avg a^n,   u^{n+1} = u^n + dt_{n+1/2} v^{n+1/2}.
//
// Kelvin-Voigt term (ADDynamicStressDivergenceTensors, alpha = 0):
//   F_KV^n = zeta (F_int^n - F_int^{n-1}) / dt_{n-1/2},
// Mass-proportional damping (quench = "mass", C = alpha M with the lumped M)
// is treated like the dashpot, at v^n: F_M^n = alpha M v^n.
// the backward difference of the stored degraded internal force, which is
// MOOSE's zeta (sigma - sigma_old)/dt . grad(test) lagged by one step. At
// fixed damage it equals zeta K v^{n-1/2}.
//
// Fixed dofs keep a = v = 0 (u = 0 there).

#ifndef MFEM_PFF_CENTRAL_DIFFERENCE_HPP
#define MFEM_PFF_CENTRAL_DIFFERENCE_HPP

#include "mfem.hpp"
#include "../domain/explicit_elasticity_operator.hpp"

namespace mfem
{
namespace pff
{

struct CentralDifferenceState
{
   Vector u, v_half;      ///< u^n, v^{n-1/2}
   Vector F_prev;         ///< F_int^{n-1}
   real_t t = 0.0;
   real_t dt_prev = 0.0;  ///< dt_{n-1/2}; 0 before the first step
   long step = 0;
};

/// Diagnostics of one step, all at t^n (before the update).
struct CentralDifferenceStepInfo
{
   real_t psi = 0.0;      ///< int psi(u^n)
   Vector v_n;            ///< v^n = v^{n-1/2} + h a^n
   Vector F_int;          ///< F_int^n
   Vector f_ext;          ///< s(t^n) b + b_c (pulse plus confinement)
   Vector f_conf;         ///< b_c (constant confinement)
   Vector f_damp;         ///< C_b v^n
   Vector f_kv;           ///< F_KV^n
   Vector f_mass;         ///< alpha M v^n
};

class CentralDifference
{
public:
   CentralDifference(const ExplicitElasticityOperator &op) : op_(op) {}

   /// Advance state from t^n to t^n + dt with damage D (true dofs, held
   /// fixed over the step), pulse factor s = Shape(t^n), damping zeta(t^n),
   /// mass-proportional damping alpha(t^n). f_add (optional) is a further
   /// internal force at t^n, e.g. the pore-pressure force of the hydromech
   /// model (rhs -= f_add; it is not part of F_int, psi or the KV term).
   void Step(CentralDifferenceState &S, const Vector &D, real_t s, real_t zeta,
             real_t dt, CentralDifferenceStepInfo &info, real_t alpha = 0.0,
             const Vector *f_add = nullptr) const
   {
      const int n = S.u.Size();
      op_.InternalForce(S.u, D, info.F_int, &info.psi);
      info.f_kv.SetSize(n);
      info.f_kv = 0.0;
      if (zeta > 0.0 && S.dt_prev > 0.0)
      {
         for (int i = 0; i < n; i++)
         {
            info.f_kv(i) = zeta * (info.F_int(i) - S.F_prev(i)) / S.dt_prev;
         }
      }
      info.f_ext.SetSize(n);
      info.f_ext.Set(s, op_.PressureVector());
      info.f_conf = op_.ConfinementVector();
      info.f_ext += info.f_conf;

      Vector cv;
      op_.DashpotForce(S.v_half, cv);
      if (alpha > 0.0) { op_.AddMassForce(alpha, S.v_half, cv); }
      Vector rhs(n);
      for (int i = 0; i < n; i++)
      {
         rhs(i) = info.f_ext(i) - info.F_int(i) - info.f_kv(i) - cv(i);
      }
      if (f_add) { rhs -= *f_add; }
      const real_t dt_avg = (S.dt_prev > 0.0) ? 0.5 * (S.dt_prev + dt) : 0.5 * dt;
      // v^n = v^{n-1/2} + (t^n - t^{n-1/2}) a^n; at the first step v_half holds v^0.
      const real_t h = 0.5 * S.dt_prev;
      Vector a(n);
      op_.SolveAcceleration(rhs, h, a, alpha);

      info.v_n.SetSize(n);
      add(S.v_half, h, a, info.v_n);
      for (int t : op_.FixedTrueDofs()) { info.v_n(t) = 0.0; }
      op_.DashpotForce(info.v_n, info.f_damp);
      info.f_mass.SetSize(n);
      info.f_mass = 0.0;
      if (alpha > 0.0) { op_.AddMassForce(alpha, info.v_n, info.f_mass); }

      S.v_half.Add(dt_avg, a);
      for (int t : op_.FixedTrueDofs()) { S.v_half(t) = 0.0; }
      S.u.Add(dt, S.v_half);
      S.F_prev = info.F_int;
      S.dt_prev = dt;
      S.t += dt;
      S.step++;
   }

private:
   const ExplicitElasticityOperator &op_;
};

} // namespace pff
} // namespace mfem

#endif
