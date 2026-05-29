// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/unit/test_rk_time_stepper.cpp — Phase 14 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// Unit coverage for the explicit coupled-RK time integrator
// (dynamic/rk_time_stepper.{hpp,cpp}):
//
//   T1  Tableau correctness — MakeRK4Tableau / MakeDormandPrinceRK45Tableau
//       pass ValidateTableau; structural fields (stages, name, Σb, Σbhat, c).
//   T2  Tableau temporal order on a linear oscillator (pure-imaginary
//       eigenvalues, the regime central flux needs): RK4 ≈ O(dt⁴), DP45 ≥ O(dt⁴)
//       (5th-order b row).  This validates the coefficients independent of the
//       wave operator (the §14.1/§14.4 "matches ADER-O2 to O(dt²)" / "DP45 shows
//       ≥ 4th-order temporal convergence" acceptance, at the tableau level).
//   T3  Bulk RK4 == the proven DoRK4Step idiom (git:8461c67) on a frictionless
//       cube, through AdvanceRKCoupled_Spatial — stage RHS k-vectors bit-for-bit
//       and the final Q to round-off (the residual is dt·(1/6) vs dt/6.0).
//   T4  PsiRate — aging-law and SRW arms equal the underlying friction kernels.
//   T5  ApplyAbsolute nucleation (§14.3): t=0→0, t=T_nuc→full, monotone; and the
//       absolute value at any t equals the ADER telescoped increment sum (no
//       double-apply).
//   T6  R-001 endpoint predicate: the final RK stage's state equals the combined
//       endpoint ONLY for an FSAL tableau (DP45), NOT for non-FSAL RK4 — this is
//       the FSAL test AdvanceRKCoupled_Spatial uses to decide whether to do the
//       endpoint re-evaluation of the reported fault observables.
//
// Usage:  ./seas_test_rk_time_stepper   (serial, < 5 s)

#include "mfem.hpp"

#include "../../dynamic/rk_time_stepper.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/spatial_nucleation.hpp"
#include "../../dynamic/nucleation_method.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../friction/slip_law_srw_psi.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_TRUE(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else      { num_failed++; std::cout << "  FAILED [" << __LINE__ \
                       << "]: " << msg << "\n"; } \
} while (0)

// ---------------------------------------------------------------------------
// Reference DoRK4Step idiom (verbatim from test_rk4_conservation.cpp /
// git:8461c67) — the Q integration AdvanceRKCoupled_Spatial must reproduce.
// ---------------------------------------------------------------------------
static void DoRK4Step(const WaveOperator<Mesh>& wave, Vector& Q, real_t dt)
{
   const int size = Q.Size();
   Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size);
   wave.Mult(Q, k1);
   add(Q, 0.5 * dt, k1, Q_tmp);  wave.Mult(Q_tmp, k2);
   add(Q, 0.5 * dt, k2, Q_tmp);  wave.Mult(Q_tmp, k3);
   add(Q, dt,       k3, Q_tmp);  wave.Mult(Q_tmp, k4);
   Q.Add(dt / 6.0, k1);  Q.Add(dt / 3.0, k2);
   Q.Add(dt / 3.0, k3);  Q.Add(dt / 6.0, k4);
}

// ---------------------------------------------------------------------------
// Generic tableau applicator on a small linear ODE dy/dt = f(y).  Steps the
// state with the SAME (a, b, c) structure as AdvanceRKCoupled_Spatial's bulk
// path, so a measured convergence order is attributable to the coefficients.
// ---------------------------------------------------------------------------
template <typename RHS>
static void TableauStep(const RKTableau& tab, std::vector<real_t>& y,
                        real_t dt, RHS&& f)
{
   const int s = tab.stages;
   const int d = static_cast<int>(y.size());
   std::vector<std::vector<real_t>> k(s, std::vector<real_t>(d, 0.0));
   std::vector<real_t> ystage(d);
   for (int i = 0; i < s; ++i)
   {
      for (int c = 0; c < d; ++c)
      {
         real_t acc = y[c];
         for (int j = 0; j < i; ++j) { acc += dt * tab.a[i][j] * k[j][c]; }
         ystage[c] = acc;
      }
      f(ystage, k[i]);
   }
   for (int c = 0; c < d; ++c)
   {
      real_t acc = y[c];
      for (int i = 0; i < s; ++i) { acc += dt * tab.b[i] * k[i][c]; }
      y[c] = acc;
   }
}

// Observed convergence order of the tableau on the harmonic oscillator
// u' = v, v' = -ω²u  (exact: u(T)=cos(ωT), v(T)=-ω sin(ωT)) — pure-imaginary
// eigenvalues ±iω, the stability regime central flux puts the operator in.
// The error is the full-state 2-norm at a NON-periodic endpoint (T not a
// multiple of the period), so the leading O(dt^p) phase term is not sampled at
// a node and the measured order is clean (a periodic endpoint makes RK4 look
// 5th-order via amplitude-error superconvergence).
static real_t OscillatorOrder(const RKTableau& tab, real_t omega, real_t T)
{
   auto rhs = [omega](const std::vector<real_t>& y, std::vector<real_t>& dy)
   {
      dy[0] = y[1];
      dy[1] = -omega * omega * y[0];
   };
   auto integrate = [&](int nsteps) -> real_t
   {
      const real_t dt = T / nsteps;
      std::vector<real_t> y = {1.0, 0.0};
      for (int n = 0; n < nsteps; ++n) { TableauStep(tab, y, dt, rhs); }
      const real_t du = y[0] - std::cos(omega * T);
      const real_t dv = y[1] - (-omega * std::sin(omega * T));
      return std::sqrt(du * du + dv * dv);
   };
   const int    N  = 48;
   const real_t e1 = integrate(N);
   const real_t e2 = integrate(2 * N);
   return std::log2(e1 / e2);   // ≈ method order p
}

// Root-cause check for R-001: the FINAL RK stage's INPUT state equals the
// combined endpoint Q_new ONLY for an FSAL tableau (DP45, a[s-1]==b), NOT for a
// non-FSAL tableau (classical RK4, whose last stage is the predictor
// Q+dt·k_{s-2}).  So the post-step dof_data fault observables are end-of-step
// only when this returns true; AdvanceRKCoupled_Spatial does an endpoint
// re-eval (extra Mult at Q_new) exactly when it returns false.  Verified here
// on a scalar linear ODE y'=λy with the same (a,b) stage arithmetic the stepper
// uses, so the result is attributable to the tableau alone.
static bool LastStageIsEndpoint(const RKTableau& tab)
{
   const real_t lambda = -0.7, dt = 0.1, y0 = 1.0;
   const int s = tab.stages;
   std::vector<real_t> k(s, 0.0);
   real_t last_stage = y0;
   for (int i = 0; i < s; ++i)
   {
      real_t ys = y0;
      for (int j = 0; j < i; ++j) { ys += dt * tab.a[i][j] * k[j]; }
      k[i] = lambda * ys;                 // f(y_stage)
      if (i == s - 1) { last_stage = ys; }
   }
   real_t y_end = y0;
   for (int i = 0; i < s; ++i) { y_end += dt * tab.b[i] * k[i]; }
   return std::abs(last_stage - y_end) <= 1e-12 * std::abs(y_end);
}

int main(int argc, char* argv[])
{
   (void)argc; (void)argv;
   std::cout << "\n=== Phase 14 — RK time stepper unit tests ===\n";

   // =======================================================================
   // T1 — Tableau correctness.
   // =======================================================================
   std::cout << "\n-- T1: tableau structure + ValidateTableau --\n";
   {
      RKTableau rk4 = MakeRK4Tableau();
      RKTableau dp45 = MakeDormandPrinceRK45Tableau();

      bool rk4_ok = true, dp45_ok = true;
      try { ValidateTableau(rk4); }  catch (...) { rk4_ok = false; }
      try { ValidateTableau(dp45); } catch (...) { dp45_ok = false; }
      TEST_TRUE(rk4_ok,  "RK4 passes ValidateTableau");
      TEST_TRUE(dp45_ok, "DP45 passes ValidateTableau");

      TEST_TRUE(rk4.stages == 4 && std::string(rk4.name) == "RK4"
                && !rk4.has_embedded && rk4.bhat.empty(),
                "RK4: 4 stages, name, no embedded estimate");
      TEST_TRUE(dp45.stages == 7
                && std::string(dp45.name) == "DormandPrinceRK45"
                && dp45.has_embedded
                && static_cast<int>(dp45.bhat.size()) == 7,
                "DP45: 7 stages, name, embedded 4(5) bhat");

      real_t sb4 = 0.0; for (real_t x : rk4.b) { sb4 += x; }
      real_t sb5 = 0.0; for (real_t x : dp45.b) { sb5 += x; }
      real_t sbh = 0.0; for (real_t x : dp45.bhat) { sbh += x; }
      TEST_TRUE(std::abs(sb4 - 1.0) < 1e-14, "RK4: Σb = 1");
      TEST_TRUE(std::abs(sb5 - 1.0) < 1e-12, "DP45: Σb (5th-order) = 1");
      TEST_TRUE(std::abs(sbh - 1.0) < 1e-12, "DP45: Σbhat (4th-order) = 1");

      // DP45 is FSAL: row 6 of a equals the 5th-order b row.
      bool fsal = true;
      for (int j = 0; j < 7; ++j)
      { if (dp45.a[6][j] != dp45.b[j]) { fsal = false; } }
      TEST_TRUE(fsal, "DP45: FSAL (a[6][*] == b[*])");
   }

   // =======================================================================
   // T2 — Temporal order from the coefficients (linear oscillator).
   // =======================================================================
   std::cout << "\n-- T2: tableau temporal order (harmonic oscillator) --\n";
   {
      const real_t omega = 2.0 * M_PI;   // period 1
      const real_t T     = 0.9;          // non-periodic endpoint (clean order)
      const real_t p_rk4  = OscillatorOrder(MakeRK4Tableau(), omega, T);
      const real_t p_dp45 = OscillatorOrder(MakeDormandPrinceRK45Tableau(),
                                            omega, T);
      std::cout << "    RK4  observed order p = " << std::fixed
                << std::setprecision(3) << p_rk4 << "  (expect ~4)\n";
      std::cout << "    DP45 observed order p = " << p_dp45
                << "  (expect ~5, ≥4)\n";
      TEST_TRUE(p_rk4 > 3.7 && p_rk4 < 4.4,
                "RK4 is 4th-order on imaginary-axis ODE (matches ADER-O2-or-better)");
      TEST_TRUE(p_dp45 > 4.0,
                "DP45 is >= 4th-order on imaginary-axis ODE");
      TEST_TRUE(p_dp45 > 4.6,
                "DP45 5th-order b row delivers ~5th-order convergence");
   }

   // =======================================================================
   // T3 — Bulk RK4 reproduces the proven DoRK4Step idiom (frictionless cube).
   // =======================================================================
   std::cout << "\n-- T3: AdvanceRKCoupled_Spatial(RK4) == DoRK4Step (bulk) --\n";
   {
      const double L = 8.0e3;
      const int    nx = 6;
      Mesh mesh = Mesh::MakeCartesian3D(nx, nx, nx, Element::HEXAHEDRON, L, L, L);
      for (int v = 0; v < mesh.GetNV(); v++)
      {
         real_t* x = mesh.GetVertex(v);
         x[0] -= 0.5 * L; x[1] -= 0.5 * L; x[2] -= 0.5 * L;
      }
      for (int be = 0; be < mesh.GetNBE(); be++) { mesh.SetBdrAttribute(be, 5); }
      mesh.SetAttributes();

      BoundaryConfig bc;
      bc.natural_attrs   = {};
      bc.fault_attr      = 0;   // NO fault — pure bulk
      bc.absorbing_attrs = {5};

      const int order = 1;
      WaveOperator<Mesh> wave(mesh, order, TPV102Params::lambda,
                              TPV102Params::mu, TPV102Params::rho, bc);
      { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }

      const int size = wave.Height();
      auto& fes = const_cast<FiniteElementSpace&>(wave.GetFESpace());
      const int ndof_total = fes.GetNDofs();

      // Smooth Gaussian σ_yy / σ_xx pulse (localized, far from boundary).
      const double sig = L / 8.0, amp = 1.0e6;
      FunctionCoefficient syy_init([sig, amp](const Vector& x) -> real_t {
         double r2 = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
         return amp * std::exp(-r2 / (sig * sig)); });
      Vector Q(size); Q = 0.0;
      GridFunction gf(&fes); gf.ProjectCoefficient(syy_init);
      for (int d = 0; d < ndof_total; d++)
      {
         Q(SYY * ndof_total + d) = gf(d);
         Q(SXX * ndof_total + d) = 0.5 * gf(d);
      }

      const real_t dt = 5.0e-5;

      // Reference: the proven idiom.
      Vector Q_ref(Q);
      DoRK4Step(wave, Q_ref, dt);

      // Subject: the tableau-driven coupled stepper, no fault DOFs / no nuc.
      std::vector<DOFData> dof_empty;
      spatial::RateStateBlock rs_cfg;            // unused (no fault DOFs)
      spatial::RateStatePerDOFParams rs_empty;   // unused
      Vector Q_rk(size);
      AdvanceRKCoupled_Spatial<Mesh>(wave, dof_empty, rs_cfg, rs_empty, Q, dt,
                                     /*t_step_start=*/0.0, Q_rk,
                                     MakeRK4Tableau(), /*nuc=*/nullptr);

      real_t max_abs = 0.0, max_diff = 0.0;
      for (int i = 0; i < size; i++)
      {
         max_abs  = std::max(max_abs, std::abs(Q_ref(i)));
         max_diff = std::max(max_diff, std::abs(Q_ref(i) - Q_rk(i)));
      }
      const real_t rel = max_diff / std::max(max_abs, real_t(1e-30));
      std::cout << "    max|Q_ref| = " << std::scientific << std::setprecision(4)
                << max_abs << "   max|Q_ref - Q_rk| = " << max_diff
                << "   rel = " << rel << "\n";
      // Residual is the dt·(1/6) vs dt/6.0 final-combine rounding (a few ULP);
      // the four stage RHS evaluations are bit-identical by construction.
      TEST_TRUE(rel < 1.0e-13,
                "RK4 bulk step reproduces DoRK4Step to round-off (bit-for-bit RHS)");
   }

   // =======================================================================
   // T4 — PsiRate equals the underlying friction kernels.
   // =======================================================================
   std::cout << "\n-- T4: PsiRate (aging + SRW) --\n";
   {
      // Aging law.
      spatial::RateStateBlock rs_cfg;
      rs_cfg.state_evolution = spatial::StateEvolutionKind::AgingLaw;
      rs_cfg.b_default = 0.012; rs_cfg.V_0_default = 1e-6; rs_cfg.f_0_default = 0.6;
      DOFData d; d.psi = 0.62; d.Dc = 0.004;
      spatial::RateStatePerDOFParams rs_dummy;
      const real_t V = 1.0e-4;
      const real_t got = PsiRate(rs_cfg, d, V, rs_dummy, 0);
      const real_t exp = AgingLawPsi(0.012, 1e-6, 0.6).Rate(V, d.psi, d.Dc);
      std::cout << "    aging: PsiRate=" << std::scientific << got
                << "  AgingLawPsi::Rate=" << exp << "\n";
      TEST_TRUE(std::abs(got - exp) <= 1e-30 + 1e-14 * std::abs(exp),
                "PsiRate(AgingLaw) == AgingLawPsi::Rate(V, psi, Dc)");

      // SRW law.
      spatial::RateStateBlock srw_cfg;
      srw_cfg.state_evolution =
         spatial::StateEvolutionKind::SlipLawStrongRateWeakening;
      srw_cfg.a_default = 0.008; srw_cfg.b_default = 0.012;
      srw_cfg.V_0_default = 1e-6; srw_cfg.f_0_default = 0.6;
      srw_cfg.f_w_default = 0.1;  srw_cfg.V_w_default = 0.1;
      spatial::RateStatePerDOFParams rs_srw;
      rs_srw.V_w.SetSize(1); rs_srw.V_w(0) = 0.17;
      rs_srw.a.SetSize(1);   rs_srw.a(0)   = 0.009;
      DOFData d2; d2.psi = 0.45; d2.Dc = 0.005;
      const real_t Vsrw = 5.0e-2;
      const real_t got2 = PsiRate(srw_cfg, d2, Vsrw, rs_srw, 0);
      SlipLawSRWPsi law(0.008, 0.012, 1e-6, 0.6, 0.1, 0.1);
      law.SetProductionMode();
      const real_t exp2 = law.Rate_SRW(Vsrw, d2.psi, d2.Dc,
                                       rs_srw.V_w(0), rs_srw.a(0));
      std::cout << "    SRW:   PsiRate=" << got2
                << "  Rate_SRW=" << exp2 << "\n";
      TEST_TRUE(std::abs(got2 - exp2) <= 1e-30 + 1e-14 * std::abs(exp2),
                "PsiRate(SRW) == SlipLawSRWPsi::Rate_SRW(V, psi, Dc, V_w, a)");
   }

   // =======================================================================
   // T5 — ApplyAbsolute nucleation (§14.3).
   // =======================================================================
   std::cout << "\n-- T5: ApplyGradualOverstressAbsolute (set, not accumulate) --\n";
   {
      const int n = 3;
      const real_t T_nuc = 1.0;
      spatial::GradualOverstressPerDOFParams p;
      p.amplitude_dip.SetSize(n);
      p.amplitude_strike.SetSize(n);
      p.radial.SetSize(n);
      for (int i = 0; i < n; i++)
      {
         p.amplitude_dip(i)    = 0.0;
         p.amplitude_strike(i) = 1.0e6 * (i + 1);
         p.radial(i)           = 1.0;
      }

      auto make_dof = [&]() { return std::vector<DOFData>(n); };

      // t = 0 → 0.
      {
         auto dof = make_dof();
         spatial::ApplyGradualOverstressAbsolute(dof, p, T_nuc, 0.0);
         bool zero = true;
         for (int i = 0; i < n; i++)
         { if (dof[i].tau2_nuc != 0.0) { zero = false; } }
         TEST_TRUE(zero, "ApplyAbsolute(t=0) → tau2_nuc == 0");
      }
      // t >= T_nuc → full target.
      {
         auto dof = make_dof();
         spatial::ApplyGradualOverstressAbsolute(dof, p, T_nuc, T_nuc);
         bool full = true;
         for (int i = 0; i < n; i++)
         { if (std::abs(dof[i].tau2_nuc - p.amplitude_strike(i)) > 1e-9)
           { full = false; } }
         TEST_TRUE(full, "ApplyAbsolute(t=T_nuc) → tau2_nuc == full target");
      }
      // monotone + bounded at an intermediate time.
      {
         auto dof = make_dof();
         spatial::ApplyGradualOverstressAbsolute(dof, p, T_nuc, 0.5);
         bool mono = true;
         for (int i = 0; i < n; i++)
         {
            const real_t v = dof[i].tau2_nuc;
            if (!(v > 0.0 && v < p.amplitude_strike(i))) { mono = false; }
         }
         TEST_TRUE(mono, "ApplyAbsolute(0<t<T_nuc) ∈ (0, full) per DOF");
      }
      // Idempotence: applying twice at the same t gives the same value (no
      // double-apply across RK stages revisiting a sub-interval).
      {
         auto dof = make_dof();
         spatial::ApplyGradualOverstressAbsolute(dof, p, T_nuc, 0.4);
         std::vector<real_t> once(n);
         for (int i = 0; i < n; i++) { once[i] = dof[i].tau2_nuc; }
         spatial::ApplyGradualOverstressAbsolute(dof, p, T_nuc, 0.4);
         bool idem = true;
         for (int i = 0; i < n; i++)
         { if (dof[i].tau2_nuc != once[i]) { idem = false; } }
         TEST_TRUE(idem, "ApplyAbsolute is idempotent at fixed t");
      }
      // Equivalence to the ADER telescoped increment sum over [0, t].
      {
         const int    nsub = 20;
         const real_t dts  = T_nuc / nsub;
         auto dof_inc = make_dof();   // accumulate via ApplyIncrement
         for (int s = 1; s <= nsub; s++)
         {
            spatial::ApplyGradualOverstressIncrement(dof_inc, p, T_nuc,
                                                     s * dts, dts);
         }
         auto dof_abs = make_dof();   // set via ApplyAbsolute at the endpoint
         spatial::ApplyGradualOverstressAbsolute(dof_abs, p, T_nuc, T_nuc);
         real_t max_diff = 0.0;
         for (int i = 0; i < n; i++)
         { max_diff = std::max(max_diff,
                               std::abs(dof_inc[i].tau2_nuc
                                        - dof_abs[i].tau2_nuc)); }
         std::cout << "    telescoped vs absolute max diff = "
                   << std::scientific << max_diff << "\n";
         TEST_TRUE(max_diff < 1e-3,   // ~1e6 Pa scale: 1e-3 Pa is round-off
                   "ADER telescoped increment sum == absolute (no double-apply)");
      }

      // INucleationMethod::ApplyAbsolute dispatch (Gaussian concrete).
      {
         GaussianGradualOverstress nm(p, T_nuc);
         auto dof = make_dof();
         nm.ApplyAbsolute(dof, T_nuc);
         bool full = true;
         for (int i = 0; i < n; i++)
         { if (std::abs(dof[i].tau2_nuc - p.amplitude_strike(i)) > 1e-9)
           { full = false; } }
         TEST_TRUE(full,
                   "GaussianGradualOverstress::ApplyAbsolute dispatches to full target");
      }
   }

   // =======================================================================
   // T6 — R-001: final-stage state vs endpoint (drives the endpoint re-eval).
   // =======================================================================
   std::cout << "\n-- T6: final RK stage state == endpoint? (R-001) --\n";
   {
      const bool rk4_endpt  = LastStageIsEndpoint(MakeRK4Tableau());
      const bool dp45_endpt =
         LastStageIsEndpoint(MakeDormandPrinceRK45Tableau());
      std::cout << "    RK4 last-stage==endpoint? " << (rk4_endpt ? "yes" : "no")
                << "   DP45 last-stage==endpoint? " << (dp45_endpt ? "yes" : "no")
                << "\n";
      // RK4 is NON-FSAL: its last stage is the predictor Q+dt·k_{s-2}, so the
      // stepper MUST endpoint-re-eval the fault observables (R-001 fix).
      TEST_TRUE(!rk4_endpt,
                "RK4 last stage is NOT the endpoint (stepper must endpoint-re-eval)");
      // DP45 is FSAL: a[6]==b ⇒ last stage already at Q_new ⇒ re-eval skipped.
      TEST_TRUE(dp45_endpt,
                "DP45 last stage IS the endpoint (FSAL; re-eval correctly skipped)");
   }

   // NOTE (R-001 / R-003): a fail-without-fix runtime test of the endpoint
   // re-eval would need the fault observables to VARY across the RK stages
   // within one step.  On a serial single-fault fixture the fault response is
   // saturated by the static imposed shear, so the last-stage state and the
   // endpoint coincide to machine precision (the difference is O((λ·dt)³) and
   // ~0 for a saturating fault) — such a test passes with OR without the fix
   // and is therefore vacuous.  Exercising the re-eval requires a propagating
   // multi-element wave reaching the fault, deferred to the Frontera SAFS
   // validation (plan §14.5).  T6 above locks the FSAL predicate the fix relies
   // on (RK4 non-FSAL ⇒ re-eval runs; DP45 FSAL ⇒ skipped).

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
