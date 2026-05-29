// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/rk_time_stepper.hpp — Phase 14 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// An explicit, tableau-driven Runge–Kutta time integrator for the SAFS
// dynamic-rupture driver — a runtime-selectable alternative to the ADER-DG
// sub-step path so that mixed/central flux can run stably at p=1 (central flux
// is non-dissipative ⇒ fault-adjacent modes sit on the imaginary axis, which
// ADER-O2's stability function does not contain but an RK of order ≥ 3 does).
//
// Two schemes ship, both through ONE tableau-driven loop:
//   - classical RK4 (the validated stepping-stone; reproduces the proven
//     git:8461c67 tpv102_driver.cpp coupled-RK4 idiom bit-for-bit on the bulk);
//   - fixed-step Dormand–Prince RK45 (the user-requested 5th-order integrator).
//
// The integrator is a HAND-WRITTEN coupled RK loop over (Q, ψ, slip): ψ/slip
// live in `DOFData`, not the MFEM `Vector`, so a stock `mfem::ODESolver` cannot
// drive the coupled system.  Each stage writes the stage-local ψ into
// `dof_data[m].psi` BEFORE calling `wave.Mult(Q^(i), ·)` so the instantaneous,
// ψ-stateless fault Riemann solve inside `Mult` reads the correct staged value.
//
// Scope (plan §14): RATE-AND-STATE ONLY (`Mult`'s fault solve is the RS
// `Evaluate`; there is no instantaneous LSW solve on the Mult path) and SCALAR
// interior flux only.  ADER is untouched; the RK path is reachable only via
// `[numerics].time_integrator = "rk4"|"rk45"` (CLI `--time-integrator`).

#ifndef MFEM_SEAS_RK_TIME_STEPPER_HPP
#define MFEM_SEAS_RK_TIME_STEPPER_HPP

#include "mfem.hpp"

#include "wave_operator.hpp"                      // WaveOperator<MeshType>
#include "fault_face_flux.hpp"                    // DOFData
#include "nucleation_method.hpp"                  // INucleationMethod
#include "../spatial/code/spatial_friction.hpp"   // spatial::RateState*
#include "../friction/state_evolution.hpp"        // AgingLawPsi (PsiRateEvaluator)
#include "../friction/slip_law_srw_psi.hpp"       // SlipLawSRWPsi (PsiRateEvaluator)

#include <vector>

namespace mfem
{
namespace seas
{

// =====================================================================
// Butcher tableau
// =====================================================================

/// @brief Explicit Runge–Kutta Butcher tableau.  `a` is strictly lower
/// triangular (explicit ⇒ a[i][j] = 0 for j >= i).  `b` is the primary
/// (highest-order) weight row; `bhat` is the embedded lower-order row used
/// only when `has_embedded` (optional adaptive control, plan §14.6).
struct RKTableau
{
   int                              stages = 0;
   std::vector<std::vector<real_t>> a;        ///< stage coefficients (s×s, strict-lower)
   std::vector<real_t>              b;         ///< primary weight row (size s)
   std::vector<real_t>              bhat;      ///< embedded weight row (size s; empty if none)
   std::vector<real_t>              c;         ///< stage abscissae (size s)
   bool                             has_embedded = false;
   const char*                      name = "";
};

/// Classical 4-stage RK4: c=(0,½,½,1), b=(1,2,2,1)/6, standard lower-triangular
/// a.  `has_embedded = false`.
RKTableau MakeRK4Tableau();

/// 7-stage Dormand–Prince DP(4,5) (FSAL): c=(0,1/5,3/10,4/5,8/9,1,1), `b` = the
/// 5th-order row, `bhat` = the 4th-order row.  `has_embedded = true`.
RKTableau MakeDormandPrinceRK45Tableau();

/// MFEM_VERIFY: Σb = 1, (Σbhat = 1 if embedded), c_i = Σ_j a_ij, strict-lower
/// triangularity, and the row/column sizes are mutually consistent with
/// `stages`.  Cheap; safe to call at the top of `AdvanceRKCoupled_Spatial`.
void ValidateTableau(const RKTableau& tab);

// =====================================================================
// Rate-and-state ψ rate (plan §14.2)
// =====================================================================

/// @brief dψ/dt at fault DOF `m`, dispatched on `rs_cfg.state_evolution`:
///   AgingLaw                    -> AgingLawPsi::Rate(V, d.psi, d.Dc) — the
///                                  ψ-space aging rate (state_evolution.hpp:190),
///                                  with global b/V0/f0 from `rs_cfg` (matching
///                                  the production aging substep iterator).
///   SlipLawStrongRateWeakening  -> SlipLawSRWPsi::Rate_SRW(V, d.psi, d.Dc,
///                                  rs.V_w(m), rs.a(m)) with the law in
///                                  production mode (R-001) and per-QP V_w/a
///                                  from the resolved `rs` (same source as
///                                  tpv104_substep_iterator.cpp).
/// `V` is the stage slip-rate magnitude (`dof_data[m].slip_rate` from the
/// stage `Mult`).  NOTE: the aging arm uses the ψ-SPACE rate (not the
/// θ-space `1 − Vθ/Dc`), because `DOFData::psi` is the ψ-space state variable
/// the whole RS infrastructure integrates; the plan's §14.2 citation
/// (state_evolution.hpp:190 == AgingLawPsi::Rate) is authoritative.
real_t PsiRate(const spatial::RateStateBlock&       rs_cfg,
               const DOFData&                        d,
               real_t                                V,
               const spatial::RateStatePerDOFParams& rs,
               int                                   m);

/// @brief Caches the rate-and-state law object(s) so `AdvanceRKCoupled_Spatial`
/// can evaluate dψ/dt per fault DOF per RK stage WITHOUT reconstructing the law
/// on every call (R-002).  The law is a stateless value-holder over the global
/// config scalars, so ONE instance built at the top of the macro-step serves
/// every (stage, DOF) — mirroring the git:8461c67 reference idiom that builds
/// `AgingLawPsi` once outside the time loop.  The semantics are identical to
/// the free `PsiRate` above, which delegates to this evaluator so the dispatch
/// lives in exactly one place.  Holds references to `rs_cfg`/`rs`; do not let
/// it outlive them.
class PsiRateEvaluator
{
public:
   PsiRateEvaluator(const spatial::RateStateBlock&        rs_cfg,
                    const spatial::RateStatePerDOFParams& rs)
      : rs_cfg_(rs_cfg), rs_(rs),
        aging_(rs_cfg.b_default, rs_cfg.V_0_default, rs_cfg.f_0_default),
        srw_(rs_cfg.a_default, rs_cfg.b_default, rs_cfg.V_0_default,
             rs_cfg.f_0_default, rs_cfg.f_w_default, rs_cfg.V_w_default)
   {
      // Production mode forbids the scalar-V_w base virtual (R-001); the SRW
      // arm uses the per-QP Rate_SRW overload.
      srw_.SetProductionMode();
   }

   /// dψ/dt at fault DOF `m`; `V` is the stage slip-rate magnitude
   /// (`dof_data[m].slip_rate` from the stage `Mult`).  See `PsiRate` for the
   /// per-arm formula citations.
   real_t operator()(const DOFData& d, real_t V, int m) const
   {
      switch (rs_cfg_.state_evolution)
      {
         case spatial::StateEvolutionKind::AgingLaw:
            return aging_.Rate(V, d.psi, d.Dc);
         case spatial::StateEvolutionKind::SlipLawStrongRateWeakening:
            MFEM_VERIFY(rs_.V_w.Size() > m && rs_.a.Size() > m,
                        "PsiRateEvaluator(SRW): rs.V_w / rs.a must be sized to "
                        "the fault DOF count; got V_w.Size()=" << rs_.V_w.Size()
                        << ", a.Size()=" << rs_.a.Size() << ", m=" << m);
            return srw_.Rate_SRW(V, d.psi, d.Dc, rs_.V_w(m), rs_.a(m));
      }
      MFEM_ABORT("PsiRateEvaluator: unhandled state_evolution = "
                 << static_cast<int>(rs_cfg_.state_evolution));
      return 0.0;
   }

private:
   const spatial::RateStateBlock&        rs_cfg_;
   const spatial::RateStatePerDOFParams& rs_;
   AgingLawPsi                           aging_;
   SlipLawSRWPsi                         srw_;
};

// =====================================================================
// Coupled RK macro-step (plan §14.1/§14.2/§14.3/§14.5)
// =====================================================================

/// @brief Advance (Q, ψ, slip) by ONE macro-step `dt_step` with the Butcher
/// tableau `tab`, coupling the bulk wave field and the rate-and-state fault
/// state through the SAME weights.
///
/// For each stage i:
///   Q^(i)            = Q + dt·Σ_{j<i} a_ij k_j
///   dof_data[m].psi  = ψ_n[m] + dt·Σ_{j<i} a_ij psi_k[j][m]   (stage-local ψ)
///   nuc->ApplyAbsolute(dof_data, t_step_start + c_i·dt)        (§14.3)
///   wave.Mult(Q^(i), k_i)   — runs the ψ-stateless fault Riemann solve,
///                             writing dof_data[m].{slip_rate,V1,V2}
///   capture sr_k/V1_k/V2_k[i][m] and psi_k[i][m] = PsiRate(...)
/// Final combine (b row):
///   Q_new            = Q + dt·Σ_i b_i k_i
///   dof_data[m].psi  = ψ_n[m] + dt·Σ_i b_i psi_k[i][m]
///   dof_data[m].slip{1,2} += dt·Σ_i b_i V{1,2}_k[i][m]   (dslip/dt = V)
///
/// `dof_data[m].slip_rate_substep_max` is the RK-stage max-|V| reduction
/// (§14.5 analogue of the ADER per-sub-step peak; the driver resets it to 0
/// before the step).  `nuc` may be nullptr (e.g. a frictionless bulk test).
///
/// Templated on `MeshType` so it serves both the serial-`Mesh` unit tests and
/// the production `ParMesh` driver.  Writes `Q_new` + `dof_data`; returns void.
template <typename MeshType>
void AdvanceRKCoupled_Spatial(WaveOperator<MeshType>&               wave,
                              std::vector<DOFData>&                 dof_data,
                              const spatial::RateStateBlock&        rs_cfg,
                              const spatial::RateStatePerDOFParams& rs,
                              const Vector&                         Q,
                              real_t                                dt_step,
                              real_t                                t_step_start,
                              Vector&                               Q_new,
                              const RKTableau&                      tab,
                              INucleationMethod*                    nuc)
{
   MFEM_VERIFY(dt_step > 0.0,
               "AdvanceRKCoupled_Spatial: dt_step must be > 0, got " << dt_step);
   ValidateTableau(tab);

   // Symmetric friction-law guard (Phase 14 RK+LSW, plan §Phase 3 req 4): this
   // stepper integrates ψ via the rate-and-state PsiRate, and the Mult-path
   // fault dispatch keys on wave.GetFaultFrictionLaw().  If the operator is
   // LSW-flagged, Mult would run EvaluateLSW (slip-stateless) while this loop
   // integrates a ψ the fault solve ignores — silent wrong physics.  Fail loud.
   // Additive MFEM_VERIFY (no arithmetic) ⇒ byte-exact RS-RK behaviour preserved.
   MFEM_VERIFY(wave.GetFaultFrictionLaw() == FaultFrictionLaw::RateAndState,
               "AdvanceRKCoupled_Spatial integrates psi, but the WaveOperator's "
               "fault_friction_law_ is not RateAndState; Mult would run the wrong "
               "fault kernel.  Use AdvanceRKCoupledLSW_Spatial for slip_weakening.");

   const int s      = tab.stages;
   const int height = Q.Size();
   const int n      = static_cast<int>(dof_data.size());

   Q_new.SetSize(height);

   // Stage-height bulk RHS samples k_0..k_{s-1} and the working stage state.
   std::vector<Vector> k(s);
   for (int i = 0; i < s; ++i) { k[i].SetSize(height); }
   Vector Q_stage(height);

   // Per-DOF ψ snapshot at step start + per-stage (ψ-rate, V1, V2) samples.
   std::vector<real_t>               psi_n(n);
   std::vector<std::vector<real_t>>  psi_k(s, std::vector<real_t>(n, 0.0));
   std::vector<std::vector<real_t>>  V1_k(s, std::vector<real_t>(n, 0.0));
   std::vector<std::vector<real_t>>  V2_k(s, std::vector<real_t>(n, 0.0));
   for (int m = 0; m < n; ++m) { psi_n[m] = dof_data[m].psi; }

   // R-002: build the ψ-rate law ONCE per macro-step (stateless over the
   // global config scalars) instead of reconstructing it per DOF per stage.
   const PsiRateEvaluator psi_rate(rs_cfg, rs);

   for (int i = 0; i < s; ++i)
   {
      // Q^(i) = Q + dt·Σ_{j<i} a_ij k_j.
      Q_stage = Q;
      for (int j = 0; j < i; ++j)
      {
         const real_t aij = tab.a[i][j];
         if (aij != 0.0) { Q_stage.Add(dt_step * aij, k[j]); }
      }

      // Stage-local ψ = ψ_n + dt·Σ_{j<i} a_ij psi_k[j]; written BEFORE Mult so
      // the fault Riemann solve sees the staged value (plan §14.1 req 2).
      for (int m = 0; m < n; ++m)
      {
         real_t psi_stage = psi_n[m];
         for (int j = 0; j < i; ++j)
         {
            psi_stage += dt_step * tab.a[i][j] * psi_k[j][m];
         }
         dof_data[m].psi = psi_stage;
      }

      // §14.3: absolute nucleation forcing at this stage's time.
      if (nuc) { nuc->ApplyAbsolute(dof_data, t_step_start + tab.c[i] * dt_step); }

      wave.Mult(Q_stage, k[i]);

      // Capture per-stage fault samples + the ψ rate at the stage state.
      for (int m = 0; m < n; ++m)
      {
         const DOFData& d = dof_data[m];
         V1_k[i][m]  = d.V1;
         V2_k[i][m]  = d.V2;
         psi_k[i][m] = psi_rate(d, d.slip_rate, m);
         // §14.5: RK-stage max-|V| reduction (driver reset this to 0 pre-step).
         if (d.slip_rate > dof_data[m].slip_rate_substep_max)
         {
            dof_data[m].slip_rate_substep_max = d.slip_rate;
         }
      }
   }

   // Final combine with the primary (b) row.
   Q_new = Q;
   for (int i = 0; i < s; ++i)
   {
      const real_t bi = tab.b[i];
      if (bi != 0.0) { Q_new.Add(dt_step * bi, k[i]); }
   }

   for (int m = 0; m < n; ++m)
   {
      real_t psi_new   = psi_n[m];
      real_t slip1_inc = 0.0;
      real_t slip2_inc = 0.0;
      for (int i = 0; i < s; ++i)
      {
         const real_t bi = tab.b[i];
         psi_new   += dt_step * bi * psi_k[i][m];
         slip1_inc += dt_step * bi * V1_k[i][m];
         slip2_inc += dt_step * bi * V2_k[i][m];
      }
      dof_data[m].psi    = psi_new;
      dof_data[m].slip1 += slip1_inc;
      dof_data[m].slip2 += slip2_inc;
   }

   // Endpoint re-evaluation of the instantaneous fault observables (R-001).
   // Only an FSAL tableau (final a-row == b row) evaluates its last stage AT
   // the step endpoint; for a non-FSAL tableau (classical RK4) the final stage
   // is the predictor Q + dt·k_{s-2}, NOT Q_new, so dof_data[m].{slip_rate,V1,
   // V2,tau*_corr,sigma_n_corr} would be phase-lagged from Q(t+dt) (the
   // git:8461c67 R-V92-K01 half-step-lag bug, which that reference fixed with
   // exactly this endpoint re-eval).  Re-run the ψ-stateless fault Riemann
   // solve once at the true endpoint (Q_new, psi_new, nuc(t+dt)) so the driver's
   // V_max diagnostic + ParaView read self-consistent end-of-step values.  DP45
   // is FSAL ⇒ its final stage already evaluated here ⇒ skipped.  ψ/slip
   // accumulators above are untouched (Evaluate is ψ-stateless and never writes
   // slip1/slip2); dof_data[m].psi already holds psi_new.
   bool last_stage_is_endpoint = true;
   for (int j = 0; j < s; ++j)
   {
      if (tab.a[s - 1][j] != tab.b[j]) { last_stage_is_endpoint = false; break; }
   }
   if (n > 0 && !last_stage_is_endpoint)
   {
      if (nuc) { nuc->ApplyAbsolute(dof_data, t_step_start + dt_step); }
      Vector k_endpoint(height);
      wave.Mult(Q_new, k_endpoint);
      for (int m = 0; m < n; ++m)
      {
         if (dof_data[m].slip_rate > dof_data[m].slip_rate_substep_max)
         {
            dof_data[m].slip_rate_substep_max = dof_data[m].slip_rate;
         }
      }
   }
}

// =====================================================================
// Coupled RK macro-step for LINEAR SLIP-WEAKENING (plan §Phase 3)
// =====================================================================

/// @brief Advance (Q, slip) by ONE macro-step `dt_step` with the Butcher
/// tableau `tab`, coupling the bulk wave field and the LSW fault slip through
/// the SAME weights.  This is the LSW sibling of `AdvanceRKCoupled_Spatial`.
///
/// LSW has NO ψ state variable; the stepped fault state is the accumulated slip
/// δ (the friction coefficient μ(δ) depends on it).  Slip is therefore BOTH the
/// stage input AND the integrated output, so the stage-local slip is written
/// into `dof_data` BEFORE each `wave.Mult` — the LSW analogue of
/// `AdvanceRKCoupled_Spatial`'s stage-local ψ write.  `Mult`'s instantaneous,
/// slip-stateless `EvaluateLSW` reads that staged δ.
///
/// For each stage i:
///   Q^(i)              = Q + dt·Σ_{j<i} a_ij k_j
///   dof_data[m].slip{1,2} = slip{1,2}_n[m] + dt·Σ_{j<i} a_ij V{1,2}_k[j][m]
///                                                       (stage-local, before Mult)
///   nuc->ApplyAbsolute(dof_data, t_step_start + c_i·dt)        (no-op for TPV205)
///   wave.Mult(Q^(i), k_i)   — runs the slip-stateless EvaluateLSW at slip^(i),
///                             writing dof_data[m].{slip_rate,V1,V2}
///   capture V{1,2}_k[i][m] = dof_data[m].V{1,2}
/// Final combine (b row):
///   Q_new                 = Q + dt·Σ_i b_i k_i
///   dof_data[m].slip{1,2} = slip{1,2}_n[m] + dt·Σ_i b_i V{1,2}_k[i][m]
///
/// COMBINE OVERWRITES, DOES NOT ACCUMULATE: because the stage loop overwrites
/// dof_data[m].slip{1,2} every stage, at combine time those fields hold the LAST
/// stage's staged slip — NOT slip*_n.  The combine therefore writes the ABSOLUTE
/// value from the step-start snapshot `slip*_n` (NOT `+=`, unlike the RS sibling
/// at the analogous site, where slip is unstaged so `+=` is correct).
///
/// `dof_data[m].slip_rate_substep_max` is the RK-stage max-|V| reduction (the
/// driver resets it to 0 before the step).  `nuc` may be nullptr.  Requires
/// `wave.GetFaultFrictionLaw() == FaultFrictionLaw::LSW`.
///
/// Templated on `MeshType` so it serves both the serial-`Mesh` unit tests and
/// the production `ParMesh` driver.  Writes `Q_new` + `dof_data`; returns void.
template <typename MeshType>
void AdvanceRKCoupledLSW_Spatial(WaveOperator<MeshType>& wave,
                                 std::vector<DOFData>&   dof_data,
                                 const Vector&           Q,
                                 real_t                  dt_step,
                                 real_t                  t_step_start,
                                 Vector&                 Q_new,
                                 const RKTableau&        tab,
                                 INucleationMethod*      nuc)
{
   MFEM_VERIFY(dt_step > 0.0,
               "AdvanceRKCoupledLSW_Spatial: dt_step must be > 0, got "
               << dt_step);
   ValidateTableau(tab);

   // Symmetric friction-law guard (plan §Phase 3 req 4): this stepper integrates
   // slip and relies on Mult dispatching to the slip-stateless EvaluateLSW.  If
   // the operator is RateAndState-flagged, Mult would run the RS `Evaluate` on
   // LSW DOFData (reads data.psi/data.b defaults → NaN) while this loop stages
   // slip the RS solve ignores — silent wrong physics.  Fail loud.
   MFEM_VERIFY(wave.GetFaultFrictionLaw() == FaultFrictionLaw::LSW,
               "AdvanceRKCoupledLSW_Spatial integrates slip and dispatches to "
               "EvaluateLSW, but the WaveOperator's fault_friction_law_ is not "
               "LSW.  Use AdvanceRKCoupled_Spatial for rate_state.");

   const int s      = tab.stages;
   const int height = Q.Size();
   const int n      = static_cast<int>(dof_data.size());

   Q_new.SetSize(height);

   // Stage-height bulk RHS samples k_0..k_{s-1} and the working stage state.
   std::vector<Vector> k(s);
   for (int i = 0; i < s; ++i) { k[i].SetSize(height); }
   Vector Q_stage(height);

   // Per-DOF slip snapshot at step start + per-stage (V1, V2) samples.  NO ψ.
   std::vector<real_t>              slip1_n(n), slip2_n(n);
   std::vector<std::vector<real_t>> V1_k(s, std::vector<real_t>(n, 0.0));
   std::vector<std::vector<real_t>> V2_k(s, std::vector<real_t>(n, 0.0));
   for (int m = 0; m < n; ++m)
   {
      slip1_n[m] = dof_data[m].slip1;
      slip2_n[m] = dof_data[m].slip2;
   }

   for (int i = 0; i < s; ++i)
   {
      // Q^(i) = Q + dt·Σ_{j<i} a_ij k_j.
      Q_stage = Q;
      for (int j = 0; j < i; ++j)
      {
         const real_t aij = tab.a[i][j];
         if (aij != 0.0) { Q_stage.Add(dt_step * aij, k[j]); }
      }

      // Stage-local slip = slip_n + dt·Σ_{j<i} a_ij V_k[j]; written BEFORE Mult
      // so the slip-stateless EvaluateLSW sees the staged δ (μ(δ) is evaluated
      // at the stage slip).  This is the LSW analogue of the RS stepper's
      // stage-local ψ write.
      for (int m = 0; m < n; ++m)
      {
         real_t slip1_stage = slip1_n[m];
         real_t slip2_stage = slip2_n[m];
         for (int j = 0; j < i; ++j)
         {
            const real_t aij = tab.a[i][j];
            slip1_stage += dt_step * aij * V1_k[j][m];
            slip2_stage += dt_step * aij * V2_k[j][m];
         }
         dof_data[m].slip1 = slip1_stage;
         dof_data[m].slip2 = slip2_stage;
      }

      // §14.3: absolute nucleation forcing at this stage's time (no-op for
      // TPV205, which has no [nucleation] block).
      if (nuc) { nuc->ApplyAbsolute(dof_data, t_step_start + tab.c[i] * dt_step); }

      wave.Mult(Q_stage, k[i]);

      // Capture per-stage fault slip-rate samples.
      for (int m = 0; m < n; ++m)
      {
         V1_k[i][m] = dof_data[m].V1;
         V2_k[i][m] = dof_data[m].V2;
         // RK-stage max-|V| reduction (driver reset this to 0 pre-step).
         if (dof_data[m].slip_rate > dof_data[m].slip_rate_substep_max)
         {
            dof_data[m].slip_rate_substep_max = dof_data[m].slip_rate;
         }
      }
   }

   // Final combine with the primary (b) row.
   Q_new = Q;
   for (int i = 0; i < s; ++i)
   {
      const real_t bi = tab.b[i];
      if (bi != 0.0) { Q_new.Add(dt_step * bi, k[i]); }
   }

   // Slip combine: ABSOLUTE from the step-start snapshot (NOT `+=`).  The stage
   // loop overwrote dof_data[m].slip{1,2}, so they no longer hold slip*_n; a
   // `+=` here would double-count the last stage's partial staged sum.
   for (int m = 0; m < n; ++m)
   {
      real_t slip1_new = slip1_n[m];
      real_t slip2_new = slip2_n[m];
      for (int i = 0; i < s; ++i)
      {
         const real_t bi = tab.b[i];
         slip1_new += dt_step * bi * V1_k[i][m];
         slip2_new += dt_step * bi * V2_k[i][m];
      }
      dof_data[m].slip1 = slip1_new;
      dof_data[m].slip2 = slip2_new;
   }

   // Endpoint re-evaluation of the instantaneous fault observables (mirrors the
   // RS sibling).  Only an FSAL tableau (final a-row == b row) evaluates its
   // last stage AT the step endpoint; for a non-FSAL tableau (classical RK4) the
   // final stage is the predictor Q + dt·k_{s-2}, NOT Q_new, so dof_data[m].
   // {slip_rate,V1,V2,tau*_corr,sigma_n_corr} would be phase-lagged from
   // Q(t+dt).  Re-run the slip-stateless EvaluateLSW once at the true endpoint
   // (Q_new, slip_new, nuc(t+dt)) so the driver's V_max diagnostic + ParaView
   // read self-consistent end-of-step values.  The slip combine above already
   // set dof_data[m].slip{1,2} = slip*_new, so EvaluateLSW sees δ at t+dt; being
   // slip-stateless it does not disturb the integrated slip.  DP45 is FSAL ⇒ its
   // final stage already evaluated at the endpoint ⇒ skipped.
   bool last_stage_is_endpoint = true;
   for (int j = 0; j < s; ++j)
   {
      if (tab.a[s - 1][j] != tab.b[j]) { last_stage_is_endpoint = false; break; }
   }
   if (n > 0 && !last_stage_is_endpoint)
   {
      if (nuc) { nuc->ApplyAbsolute(dof_data, t_step_start + dt_step); }
      Vector k_endpoint(height);
      wave.Mult(Q_new, k_endpoint);
      for (int m = 0; m < n; ++m)
      {
         if (dof_data[m].slip_rate > dof_data[m].slip_rate_substep_max)
         {
            dof_data[m].slip_rate_substep_max = dof_data[m].slip_rate;
         }
      }
   }
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RK_TIME_STEPPER_HPP
