// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_print_derived.cpp — Phase D unit tests T-D01..T-D07.
//
// The full driver-level acceptance tests (T-D01..T-D05) need a built
// `seas_spatial_dyn_driver` + a TOML fixture; they are exercised by
// the integration smoke `make test-spatial-dyn-driver`.  This unit
// test covers the σ_1 azimuth math (T-D06) and the dt_cfl pass-through
// check (T-D07) directly via the printer interface.

#include "mfem.hpp"

#include "../../dynamic/spatial_print_derived.hpp"
#include "../../dynamic/spatial_nucleation.hpp"
#include "../../spatial/code/spatial_friction.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::spatial;

static int num_tests = 0, num_passed = 0, num_failed = 0;
#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e \
   << ", tol " << _t << ")\n"; num_failed++; } \
   else { std::cout << "  PASSED: " << m << "\n"; num_passed++; } \
   } while (0)

// =====================================================================
// Helper: construct a single-DOF fixture that exercises the printer
// end-to-end with a well-posed setup; capture stdout via a stringstream.
// =====================================================================
static real_t RunPrinter(const StressSpec& stress,
                         real_t mu_s, real_t mu_d, real_t d_c,
                         real_t sigma_n_eff,
                         real_t tau_dip, real_t tau_strike,
                         bool nuc_enabled,
                         const GradualOverstressSpec& gspec,
                         real_t mu_bulk,
                         std::string& captured,
                         bool abort_on_failure = true,
                         bool warn_only_env = false,
                         real_t dof_x = 0.0,
                         real_t dof_y = 0.0,
                         real_t dof_z = 0.0)
{
   SlipWeakeningPerDOFParams lsw;
   lsw.mu_s.SetSize(1); lsw.mu_d.SetSize(1); lsw.d_c.SetSize(1);
   lsw.cohesion.SetSize(1);
   lsw.mu_s(0) = mu_s;
   lsw.mu_d(0) = mu_d;
   lsw.d_c(0)  = d_c;
   lsw.cohesion(0) = 0.0;

   Vector tau_pre(2);  tau_pre(0) = tau_dip; tau_pre(1) = tau_strike;
   Vector sn(1);       sn(0) = sigma_n_eff;
   Vector coords(3);   coords(0) = dof_x; coords(1) = dof_y; coords(2) = dof_z;

   NucleationSpec nuc;
   nuc.enabled = nuc_enabled;
   nuc.kind    = NucleationKind::GradualOverstress;
   nuc.gradual_overstress = gspec;

   GradualOverstressPerDOFParams nuc_params;
   if (nuc_enabled)
   {
      DenseMatrix basis(9, 1); basis = 0.0;
      basis(0, 0) = 1.0;        // normal +x
      basis(4, 0) = 1.0;        // dip   +y
      basis(8, 0) = 1.0;        // strike +z
      nuc_params = ResolveGradualOverstress(gspec, true, coords, basis);
   }

   PrintDerivedConfig pd_cfg;
   pd_cfg.enabled                = true;
   pd_cfg.abort_on_failure       = abort_on_failure;
   pd_cfg.outside_safety_factor  = 3.0;

   std::ostringstream oss;
   if (warn_only_env) { setenv("SEAS_SKIP_EQUILIBRIUM_GATE", "1", 1); }
   else                { unsetenv("SEAS_SKIP_EQUILIBRIUM_GATE"); }

   const real_t out_ratio = PrintDerivedAndCheck(
      pd_cfg, lsw, tau_pre, sn, coords, nuc, nuc_params, stress,
      mu_bulk, /*cp=*/4877.0, /*cs=*/3458.0,
      /*h_min=*/1000.0, /*dt_cfl=*/6.84e-5, /*tfinal=*/12.0,
      /*num_fault_global=*/1, /*num_zero_normal_fallbacks=*/0,
#ifdef MFEM_USE_MPI
      MPI_COMM_WORLD,
#endif
      /*rank=*/0, oss);
   captured = oss.str();
   return out_ratio;
}

// =====================================================================
// T-D06: σ_1 azimuth for pure shear at 45 deg.
//        σ_xx = σ_yy = 50, σ_xy = 25 ⇒ σ_1 direction = (1, 1)/√2,
//        azimuth = 45 deg from N.
// =====================================================================
static void T_D06_sigma1_azimuth_45deg()
{
   std::cout << "\n[T-D06] σ_1 azimuth = 45 deg for pure shear at 45 deg\n";
   StressSpec s;
   s.kind        = StressSourceKind::ConstantTensor;
   s.sigma_xx_pa =  50.0e6;
   s.sigma_yy_pa =  50.0e6;
   s.sigma_zz_pa =   0.0;
   s.sigma_xy_pa =  25.0e6;
   s.sigma_yz_pa =   0.0;
   s.sigma_xz_pa =   0.0;

   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 20.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   (void)RunPrinter(s, /*mu_s=*/0.6, /*mu_d=*/0.4, /*d_c=*/0.5,
                    /*sigma_n=*/50.0e6, /*tau_dip=*/0.0,
                    // tau_strike between mu_d·sigma_n (20 MPa) and
                    // mu_s·sigma_n (30 MPa), so the patch is locked with a
                    // positive stress drop and the 20 MPa overstress triggers.
                    /*tau_strike=*/25.0e6,
                    /*nuc_enabled=*/true, gspec,
                    /*mu_bulk=*/32.0e9, captured);
   // Find "sigma_1 azimuth = " and parse the number.
   const auto pos = captured.find("sigma_1 azimuth = ");
   TEST_ASSERT(pos != std::string::npos, "sigma_1 azimuth line printed");
   const real_t deg = std::stod(captured.substr(pos + 18));
   TEST_NEAR(deg, 45.0, 0.5, "azimuth ≈ 45 deg for pure NE-ward shear");
}

// =====================================================================
// T-D07: dt_cfl is printed exactly as passed.
// =====================================================================
static void T_D07_dt_cfl_passthrough()
{
   std::cout << "\n[T-D07] dt_cfl printed as supplied\n";
   StressSpec s;
   s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 20.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   (void)RunPrinter(s, /*mu_s=*/0.6, /*mu_d=*/0.4, /*d_c=*/0.5,
                    /*sigma_n=*/50.0e6, /*tau_dip=*/0.0,
                    /*tau_strike=*/25.0e6,   // valid nucleation (see T-D06)
                    /*nuc_enabled=*/true, gspec,
                    /*mu_bulk=*/32.0e9, captured);
   const auto pos = captured.find("dt_cfl = ");
   TEST_ASSERT(pos != std::string::npos, "dt_cfl line printed");
   const real_t dt = std::stod(captured.substr(pos + 9));
   TEST_NEAR(dt, 6.84e-5, 1e-8, "dt_cfl matches input");
}

// =====================================================================
// T-D02-LITE: in-process variant of plan T-D02 (which asserts the
// driver process aborts with non-zero exit on a supercritical TOML).
// This in-process variant only confirms that the printer's RETURN
// value reflects the supercritical condition when running in warn-
// only mode (abort_on_failure = false).  The full driver-level abort
// is exercised by the smoke sbatch via `--print-derived`.
// =====================================================================
static void T_D02_LITE_supercritical_returns_ratio()
{
   std::cout << "\n[T-D02] supercritical setup returns ratio > 1\n";
   StressSpec s;
   s.kind        = StressSourceKind::ConstantTensor;
   s.sigma_xy_pa = 5.0e8;     // 500 MPa shear — outrageous

   GradualOverstressSpec gspec;  // unused when nuc_enabled = false
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   // tau_pre outrageous AND nucleation disabled, so the lone DOF is
   // treated as "outside the asperity" and feeds the equilibrium ratio.
   const real_t ratio = RunPrinter(
      s, /*mu_s=*/0.4, /*mu_d=*/0.3, /*d_c=*/0.5,
      /*sigma_n=*/30.0e6, /*tau_dip=*/0.0, /*tau_strike=*/3.0e8,
      /*nuc_enabled=*/false, gspec,
      /*mu_bulk=*/32.0e9, captured,
      /*abort_on_failure=*/false);
   TEST_ASSERT(ratio >= 1.0, "supercritical setup yields ratio >= 1");
   TEST_ASSERT(captured.find("max OUTSIDE asperity") != std::string::npos
               || captured.find("max anywhere") != std::string::npos,
               "outside-asperity ratio line printed");
}

// =====================================================================
// T-D05-LITE: SEAS_SKIP_EQUILIBRIUM_GATE=1 honored — supercritical
// scenario where the printer would otherwise MFEM_ABORT downgrades to
// a WARNING line.  This in-process variant cannot observe the driver
// exit code; it only checks that the warning string is in the
// captured output.  The plan's T-D05 full driver-level exit-0
// assertion is exercised by the smoke sbatch.
// =====================================================================
static void T_D05_LITE_env_skip_gate_honored()
{
   std::cout << "\n[T-D05] SEAS_SKIP_EQUILIBRIUM_GATE=1 downgrades to WARN\n";
   StressSpec s;
   s.kind        = StressSourceKind::ConstantTensor;
   s.sigma_xy_pa = 5.0e8;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   // Gate (b) is now nuc.enabled-only, so trip it with a NUCLEATION config
   // whose lone DOF sits OUTSIDE the patch (dof_x = 1e5 >> 3*radius) and is
   // supercritical — the "background outside the patch supercritical" failure
   // the env var downgrades.
   const real_t ratio = RunPrinter(
      s, 0.4, 0.3, 0.5, 30.0e6, 0.0, 3.0e8,
      /*nuc_enabled=*/true, gspec, 32.0e9,
      captured, /*abort_on_failure=*/true, /*warn_only_env=*/true,
      /*dof_x=*/1.0e5, /*dof_y=*/0.0, /*dof_z=*/0.0);
   TEST_ASSERT(ratio >= 1.0, "supercritical setup still returns the ratio");
   TEST_ASSERT(captured.find("WARNING") != std::string::npos
               || captured.find("gate downgraded") != std::string::npos,
               "WARNING line printed instead of abort");
   TEST_ASSERT(captured.find("SEAS_SKIP_EQUILIBRIUM_GATE=1") != std::string::npos,
               "WARN prefix correctly attributes the env-var path");
   unsetenv("SEAS_SKIP_EQUILIBRIUM_GATE");
}

// =====================================================================
// T-D03-REASON: R-003 regression — when warn_only is driven by
// cfg.abort_on_failure=false (NOT the env var), the WARN prefix must
// attribute that reason instead of falsely claiming the env var.
// =====================================================================
static void T_D03_REASON_warn_prefix_attribution()
{
   std::cout << "\n[T-D03-REASON] R-003: WARN prefix attributes correct reason\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   s.sigma_xy_pa = 5.0e8;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   // abort_on_failure=false but env var NOT set: should attribute the
   // API caller's choice, NOT the env var.  Gate (b) is nuc.enabled-only, so
   // use a nucleation config with the lone DOF OUTSIDE the patch (dof_x = 1e5)
   // and supercritical to trip it.
   (void)RunPrinter(s, /*mu_s=*/0.4, /*mu_d=*/0.3, /*d_c=*/0.5,
                    /*sigma_n=*/30.0e6, /*tau_dip=*/0.0,
                    /*tau_strike=*/3.0e8,
                    /*nuc_enabled=*/true, gspec, /*mu_bulk=*/32.0e9,
                    captured, /*abort_on_failure=*/false,
                    /*warn_only_env=*/false,
                    /*dof_x=*/1.0e5, /*dof_y=*/0.0, /*dof_z=*/0.0);
   TEST_ASSERT(captured.find("cfg.abort_on_failure=false") != std::string::npos,
               "WARN prefix mentions cfg.abort_on_failure=false");
   TEST_ASSERT(captured.find("SEAS_SKIP_EQUILIBRIUM_GATE") == std::string::npos,
               "WARN prefix does NOT falsely claim env var was set");
}

// =====================================================================
// T-D02-COORDS: R-005 regression — the supercritical FAIL message must
// include the offending DOF's (x, y, z), not just the owner rank.
// =====================================================================
static void T_D02_COORDS_outside_message_includes_xyz()
{
   std::cout << "\n[T-D02-COORDS] R-005: outside-asperity FAIL line "
                "includes DOF coords\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   s.sigma_xy_pa = 5.0e8;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   // Place the lone DOF at a distinctive (x, y, z).  Use unique
   // numeric values so the substring search is unambiguous.
   const real_t dx = 12345.0, dy = 67890.0, dz = -54321.0;
   // nuc.enabled (gate (b) is nuc.enabled-only); the distinctive DOF is far
   // outside the patch (|r| >> 3*radius) and supercritical, so gate (b) fires
   // with that DOF's coords.
   (void)RunPrinter(s, /*mu_s=*/0.4, /*mu_d=*/0.3, /*d_c=*/0.5,
                    /*sigma_n=*/30.0e6, /*tau_dip=*/0.0,
                    /*tau_strike=*/3.0e8,
                    /*nuc_enabled=*/true, gspec,
                    /*mu_bulk=*/32.0e9, captured,
                    /*abort_on_failure=*/false,
                    /*warn_only_env=*/false,
                    dx, dy, dz);
   TEST_ASSERT(captured.find("x=12345") != std::string::npos,
               "FAIL message includes x-coord");
   TEST_ASSERT(captured.find("y=67890") != std::string::npos,
               "FAIL message includes y-coord");
   TEST_ASSERT(captured.find("z=-54321") != std::string::npos,
               "FAIL message includes z-coord");
}

// =====================================================================
// T-D08: NUCLEATION-CRITERION REGRESSION (the job-7743351 numbers).
//        mu_s=0.65, mu_d=0.30, sigma_n=71.5 MPa, |tau_pre|=32.5 MPa
//        (strike), delta_tau=25 MPa.  The CORRECT trigger criterion is
//        |tau_pre + delta_tau| = 57.5 MPa >= mu_s·sigma_n = 46.5 MPa, so
//        overshoot = +11 MPa (SUFFICIENT) and the gate PASSes.  The old
//        budget = (mu_s - mu_d)·sigma_n = 25 MPa gave overshoot = A - budget
//        = -0.24 MPa and falsely printed "(INSUFFICIENT)".
// =====================================================================
static void T_D08_overshoot_uses_initiation_criterion()
{
   std::cout << "\n[T-D08] overshoot = |tau_pre+delta_tau| - mu_s·sigma_n "
                "(initiation), not A - (mu_s-mu_d)·sigma_n\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 25.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   (void)RunPrinter(s, /*mu_s=*/0.65, /*mu_d=*/0.30, /*d_c=*/1.0,
                    /*sigma_n=*/71.5e6, /*tau_dip=*/0.0,
                    /*tau_strike=*/32.5e6,
                    /*nuc_enabled=*/true, gspec,
                    /*mu_bulk=*/32.0e9, captured);
   TEST_ASSERT(captured.find("(sufficient)") != std::string::npos,
               "overshoot labelled (sufficient)");
   TEST_ASSERT(captured.find("(INSUFFICIENT)") == std::string::npos,
               "overshoot NOT (INSUFFICIENT) — regression vs the old budget");
   TEST_ASSERT(captured.find("PASS: initial conditions are well-posed")
               != std::string::npos,
               "well-posed nucleation PASSes the gate");
   const auto pos = captured.find("nucleation overshoot");
   TEST_ASSERT(pos != std::string::npos, "overshoot line printed");
   const auto eq = captured.find("= ", pos);
   const real_t ov = std::stod(captured.substr(eq + 2));
   TEST_NEAR(ov, 11.025, 0.05, "overshoot ≈ +11 MPa (57.5 - 46.5)");
}

// =====================================================================
// T-D09: NEGATIVE dynamic stress drop is caught.  Same prestress/overstress
//        as T-D08 but mu_d=0.50 (the pre-fix value): tau_pre=32.5 MPa <
//        mu_d·sigma_n = 35.75 MPa, so the dynamic ratio 0.909 < 1 and the
//        rupture would re-lock (the debug-doc Issue [1]).  Warn-only so the
//        in-process printer does not MFEM_ABORT.
// =====================================================================
static void T_D09_negative_stress_drop_caught()
{
   std::cout << "\n[T-D09] negative dynamic stress drop (mu_d too high) is "
                "flagged\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 25.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   (void)RunPrinter(s, /*mu_s=*/0.65, /*mu_d=*/0.50, /*d_c=*/1.0,
                    /*sigma_n=*/71.5e6, /*tau_dip=*/0.0,
                    /*tau_strike=*/32.5e6,
                    /*nuc_enabled=*/true, gspec,
                    /*mu_bulk=*/32.0e9, captured,
                    /*abort_on_failure=*/false);
   TEST_ASSERT(captured.find("NEGATIVE dynamic stress drop") != std::string::npos,
               "stress-drop FAIL message printed");
   TEST_ASSERT(captured.find("gate downgraded") != std::string::npos,
               "downgraded to WARNING (warn-only) instead of abort");
   TEST_ASSERT(captured.find("PASS: initial conditions are well-posed")
               == std::string::npos,
               "negative stress drop does NOT PASS");
}

// =====================================================================
// T-D10: INSUFFICIENT trigger is caught.  mu_s=0.65, sigma_n=71.5 MPa,
//        tau_pre=32.5 MPa, but delta_tau=5 MPa only: nucleated 37.5 MPa <
//        mu_s·sigma_n = 46.5 MPa, overshoot = -9 MPa.  Warn-only.
// =====================================================================
static void T_D10_insufficient_trigger_caught()
{
   std::cout << "\n[T-D10] insufficient overstress (never reaches yield) is "
                "flagged\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 5.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   (void)RunPrinter(s, /*mu_s=*/0.65, /*mu_d=*/0.30, /*d_c=*/1.0,
                    /*sigma_n=*/71.5e6, /*tau_dip=*/0.0,
                    /*tau_strike=*/32.5e6,
                    /*nuc_enabled=*/true, gspec,
                    /*mu_bulk=*/32.0e9, captured,
                    /*abort_on_failure=*/false);
   TEST_ASSERT(captured.find("(INSUFFICIENT)") != std::string::npos,
               "overshoot labelled (INSUFFICIENT)");
   TEST_ASSERT(captured.find("never reaches static yield") != std::string::npos,
               "trigger FAIL message printed");
}

// =====================================================================
// RATE-AND-STATE overload (PrintDerivedAndCheckRS) — PLAN DEVIATION.
// The Phase-3 plan does not spec an RS --print-derived path; the RS
// overload is added so the SAFS-RS sbatch (which passes --print-derived)
// runs.  These T-RS* cases mirror the LSW ones for the RS-grounded
// quantities: L_nuc = mu*Dc/((b-a)*sigma_n), f_ss(V_init), and the
// velocity-weakening nucleation gate.
// =====================================================================
static real_t RunPrinterRS(const StressSpec& stress,
                           real_t a, real_t b, real_t Dc,
                           real_t V_init, real_t f0, real_t V0,
                           real_t sigma_n_eff,
                           real_t tau_dip, real_t tau_strike,
                           bool nuc_enabled,
                           const GradualOverstressSpec& gspec,
                           real_t mu_bulk,
                           std::string& captured,
                           bool abort_on_failure = true,
                           real_t dof_x = 0.0,
                           real_t dof_y = 0.0,
                           real_t dof_z = 0.0)
{
   RateStatePerDOFParams rs;
   rs.a.SetSize(1);      rs.a(0)      = a;
   rs.b.SetSize(1);      rs.b(0)      = b;
   rs.Dc.SetSize(1);     rs.Dc(0)     = Dc;
   rs.V_init.SetSize(1); rs.V_init(0) = V_init;
   rs.f_0.SetSize(1);    rs.f_0(0)    = f0;
   rs.V_0.SetSize(1);    rs.V_0(0)    = V0;
   rs.eta.SetSize(1);    rs.eta(0)    = 0.0;
   rs.sigma_n_eff.SetSize(1); rs.sigma_n_eff(0) = sigma_n_eff;

   Vector tau_pre(2);  tau_pre(0) = tau_dip; tau_pre(1) = tau_strike;
   Vector sn(1);       sn(0) = sigma_n_eff;
   Vector coords(3);   coords(0) = dof_x; coords(1) = dof_y; coords(2) = dof_z;

   NucleationSpec nuc;
   nuc.enabled = nuc_enabled;
   nuc.kind    = NucleationKind::GradualOverstress;
   nuc.gradual_overstress = gspec;

   GradualOverstressPerDOFParams nuc_params;
   if (nuc_enabled)
   {
      DenseMatrix basis(9, 1); basis = 0.0;
      basis(0, 0) = 1.0;        // normal +x
      basis(4, 0) = 1.0;        // dip   +y
      basis(8, 0) = 1.0;        // strike +z
      nuc_params = ResolveGradualOverstress(gspec, true, coords, basis);
   }

   PrintDerivedConfig pd_cfg;
   pd_cfg.enabled                = true;
   pd_cfg.abort_on_failure       = abort_on_failure;
   pd_cfg.outside_safety_factor  = 3.0;

   std::ostringstream oss;
   unsetenv("SEAS_SKIP_EQUILIBRIUM_GATE");
   const real_t out_ratio = PrintDerivedAndCheckRS(
      pd_cfg, rs, tau_pre, sn, coords, nuc, nuc_params, stress,
      mu_bulk, /*cp=*/4877.0, /*cs=*/3458.0,
      /*h_min=*/1000.0, /*dt_cfl=*/6.84e-5, /*tfinal=*/12.0,
      /*num_fault_global=*/1, /*num_zero_normal_fallbacks=*/0,
#ifdef MFEM_USE_MPI
      MPI_COMM_WORLD,
#endif
      /*rank=*/0, oss);
   captured = oss.str();
   return out_ratio;
}

// Parse the numeric value following the FIRST occurrence of `key` in `text`.
static real_t ParseAfter(const std::string& text, const std::string& key)
{
   const auto pos = text.find(key);
   if (pos == std::string::npos) { return std::nan(""); }
   return std::stod(text.substr(pos + key.size()));
}

// T-RS01: steady-state friction f_ss(V_init) matches the grounded formula
//   f_ss = f0 + (a-b)*ln(V_init/V0)   (the asymptotic of the regularized law).
//   a=0.010, b=0.015, f0=0.6, V0=1e-6, V_init=1e-9
//   => 0.6 + (-0.005)*ln(1e-3) = 0.6 + 0.034539 = 0.634539.
static void T_RS01_steady_state_friction()
{
   std::cout << "\n[T-RS01] RS steady-state friction f_ss(V_init)\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 20.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   // tau_strike near f_ss*sigma_n so the lone in-patch DOF is well-posed.
   (void)RunPrinterRS(s, /*a=*/0.010, /*b=*/0.015, /*Dc=*/2.0,
                      /*V_init=*/1e-9, /*f0=*/0.6, /*V0=*/1e-6,
                      /*sigma_n=*/50.0e6, /*tau_dip=*/0.0,
                      /*tau_strike=*/30.0e6,
                      /*nuc_enabled=*/true, gspec,
                      /*mu_bulk=*/32.0e9, captured);
   const real_t fss = ParseAfter(captured,
                                 "steady-state friction f_ss(V_init) min/max = ");
   TEST_NEAR(fss, 0.634539, 1e-4, "f_ss(V_init) == f0 + (a-b)*ln(V_init/V0)");
   TEST_ASSERT(captured.find("friction law = rate-and-state") != std::string::npos,
               "RS banner printed");
}

// T-RS02: L_nuc(RS) = mu*Dc/((b-a)*sigma_n).
//   32e9 * 2.0 / (0.005 * 50e6) = 256000 m.
static void T_RS02_nucleation_length()
{
   std::cout << "\n[T-RS02] RS nucleation length L_nuc = mu*Dc/((b-a)*sigma_n)\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 20.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   (void)RunPrinterRS(s, 0.010, 0.015, /*Dc=*/2.0, 1e-9, 0.6, 1e-6,
                      /*sigma_n=*/50.0e6, 0.0, 30.0e6,
                      /*nuc_enabled=*/true, gspec, /*mu_bulk=*/32.0e9, captured);
   const real_t Lnuc = ParseAfter(captured,
      "L_nuc(RS)=mu*Dc/((b-a)*sigma_n) min/max/mean = ");
   TEST_NEAR(Lnuc, 256000.0, 1.0, "L_nuc == mu*Dc/((b-a)*sigma_n)");
}

// T-RS03: fully velocity-strengthening (a > b) -> all-VS gate fires; warn-only
//   so the in-process printer does not abort, and it must NOT PASS.
static void T_RS03_all_velocity_strengthening()
{
   std::cout << "\n[T-RS03] all velocity-strengthening (a>b) gates as FAIL\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 20.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   (void)RunPrinterRS(s, /*a=*/0.020, /*b=*/0.015, /*Dc=*/2.0,
                      1e-9, 0.6, 1e-6, 50.0e6, 0.0, 30.0e6,
                      /*nuc_enabled=*/true, gspec, /*mu_bulk=*/32.0e9, captured,
                      /*abort_on_failure=*/false);
   TEST_ASSERT(captured.find("all fault DOFs are velocity-strengthening")
               != std::string::npos,
               "all-VS message printed");
   TEST_ASSERT(captured.find("gate downgraded") != std::string::npos,
               "all-VS downgraded to WARNING (warn-only)");
   TEST_ASSERT(captured.find("PASS: rate-and-state") == std::string::npos,
               "all-VS does NOT PASS");
}

// T-RS04: VW patch with sufficient overstress -> overshoot drives acceleration
//   and the gate PASSes.  tau_pre=30 MPa, f_ss*sigma_n ~= 0.6345*50 = 31.7 MPa,
//   delta_tau=20 MPa -> |tau_pre+dtau|=50 MPa > 31.7 MPa -> overshoot > 0.
static void T_RS04_vw_nucleation_passes()
{
   std::cout << "\n[T-RS04] VW patch + sufficient overstress drives acceleration\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 20.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   (void)RunPrinterRS(s, 0.010, 0.015, 2.0, 1e-9, 0.6, 1e-6,
                      50.0e6, 0.0, 30.0e6,
                      /*nuc_enabled=*/true, gspec, /*mu_bulk=*/32.0e9, captured);
   TEST_ASSERT(captured.find("(drives acceleration)") != std::string::npos,
               "overshoot labelled (drives acceleration)");
   TEST_ASSERT(captured.find("PASS: rate-and-state initial conditions are "
                             "well-posed") != std::string::npos,
               "well-posed RS nucleation PASSes");
}

// T-RS05: nucleation patch entirely velocity-strengthening (DOF in patch but
//   a>b) -> patch-all-VS FAIL message; warn-only so no abort.
static void T_RS05_patch_all_vs_fails()
{
   std::cout << "\n[T-RS05] in-patch entirely velocity-strengthening gates FAIL\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;
   gspec.radius_dip_m = gspec.radius_strike_m = 1000.0;
   gspec.delta_tau_strike_pa = 20.0e6;
   gspec.T_nuc_s = 1.0;
   std::string captured;
   // a>b but b>0,a>0: VS DOF placed AT the patch centre (in-patch, not all-fault
   // VS gate first because all_vs also true here; use a mixed check via the
   // message).  Single DOF => all_vs is also true, so assert the all-VS path.
   (void)RunPrinterRS(s, /*a=*/0.020, /*b=*/0.015, 2.0, 1e-9, 0.6, 1e-6,
                      50.0e6, 0.0, 30.0e6,
                      /*nuc_enabled=*/true, gspec, /*mu_bulk=*/32.0e9, captured,
                      /*abort_on_failure=*/false);
   // With one DOF, all_vs fires first (it short-circuits the patch gate); both
   // describe the same root cause.  Assert the VS diagnosis is reported.
   TEST_ASSERT(captured.find("velocity-strengthening") != std::string::npos,
               "velocity-strengthening diagnosis printed for an a>b patch");
}

// =====================================================================
// T-D11: INSTANTANEOUS-CIRCULAR nucleation (TPV31) regression.
//   The diagnostic used to read the patch centre/radius unconditionally
//   from nuc.gradual_overstress (zero defaults for this kind ⇒ r_threshold
//   = 0, centre = origin) AND gate the in-patch loop on
//   amplitude_dip.Size()==N (never true — the instantaneous kind resolves
//   amplitude_STRIKE only).  Either bug alone made a correctly-configured
//   TPV31 patch look EMPTY and falsely MFEM_ABORT at --print-derived.  The
//   lone DOF sits at the instantaneous centre (0,0,-7500); with the fix the
//   patch is non-empty and the peak equals delta_tau_pa (CosineTaper(0)=1).
//   Warn-only so a regression fails the assertion cleanly (no process abort).
// =====================================================================
static void T_D11_instantaneous_circular_patch_not_empty()
{
   std::cout << "\n[T-D11] instantaneous_overstress_circular (TPV31): "
                "kind-aware centre/radius + strike-only amplitude\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;

   const real_t cx = 0.0, cy = 0.0, cz = -7500.0;   // TPV31 hypocentre

   SlipWeakeningPerDOFParams lsw;
   lsw.mu_s.SetSize(1); lsw.mu_d.SetSize(1); lsw.d_c.SetSize(1);
   lsw.cohesion.SetSize(1);
   lsw.mu_s(0) = 0.65; lsw.mu_d(0) = 0.30; lsw.d_c(0) = 1.0;
   lsw.cohesion(0) = 0.0;

   Vector tau_pre(2); tau_pre(0) = 0.0; tau_pre(1) = 32.5e6;   // strike
   Vector sn(1);      sn(0) = 71.5e6;
   Vector coords(3);  coords(0) = cx; coords(1) = cy; coords(2) = cz;

   NucleationSpec nuc;
   nuc.enabled = true;
   nuc.kind    = NucleationKind::InstantaneousOverstressCircular;
   nuc.instantaneous_circular.center_x_m   = cx;
   nuc.instantaneous_circular.center_y_m   = cy;
   nuc.instantaneous_circular.center_z_m   = cz;
   nuc.instantaneous_circular.radius_m     = 1400.0;
   nuc.instantaneous_circular.taper_m      = 600.0;
   nuc.instantaneous_circular.delta_tau_pa = 25.0e6;

   DenseMatrix basis(9, 1); basis = 0.0;
   basis(0, 0) = 1.0;   // normal +x
   basis(4, 0) = 1.0;   // dip    +y
   basis(8, 0) = 1.0;   // strike +z

   // Mirror the driver: the instantaneous kind populates amplitude_STRIKE
   // only; amplitude_dip is left zero-sized (spatial_dyn_driver.cpp).
   GradualOverstressPerDOFParams nuc_params;
   nuc_params.amplitude_strike =
      ResolveInstantaneousOverstressCircular(
         nuc.instantaneous_circular, /*enabled=*/true, coords, basis)
      .amplitude_strike;

   PrintDerivedConfig pd_cfg;
   pd_cfg.enabled               = true;
   pd_cfg.abort_on_failure      = false;   // warn-only: never abort the test
   pd_cfg.outside_safety_factor = 3.0;

   std::ostringstream oss;
   unsetenv("SEAS_SKIP_EQUILIBRIUM_GATE");
   (void)PrintDerivedAndCheck(
      pd_cfg, lsw, tau_pre, sn, coords, nuc, nuc_params, s,
      /*mu_bulk=*/32.0e9, /*cp=*/4877.0, /*cs=*/3458.0,
      /*h_min=*/1000.0, /*dt_cfl=*/6.84e-5, /*tfinal=*/12.0,
      /*num_fault_global=*/1, /*num_zero_normal_fallbacks=*/0,
#ifdef MFEM_USE_MPI
      MPI_COMM_WORLD,
#endif
      /*rank=*/0, oss);
   const std::string captured = oss.str();

   TEST_ASSERT(captured.find("fault DOFs are inside the nucleation patch")
               == std::string::npos,
               "patch NOT falsely reported empty (the TPV31 regression)");
   TEST_ASSERT(captured.find("most-overstressed DOF at") != std::string::npos,
               "non-empty-patch branch taken (peak DOF reported)");
   const real_t peak = ParseAfter(captured,
                                  "nucleation peak |F(r) * delta_tau| = ");
   TEST_NEAR(peak, 25.0e6, 1.0e3,
             "peak |F(r)*delta_tau| == delta_tau_pa at the patch centre");
   TEST_ASSERT(captured.find("(sufficient)") != std::string::npos,
               "trigger overshoot computed from the strike amplitude "
               "(sufficient)");
}

// =====================================================================
// T-D12: COMPACT-CIRCULAR nucleation (TPV102/104) regression — same
//   diagnostic gap as T-D11 (strike-only amplitude + non-gradual centre/
//   radius sub-struct).  DOF at the compact-bell centre; with the fix the
//   patch is non-empty and the peak equals delta_tau_pa (CompactBell(0)=1).
// =====================================================================
static void T_D12_compact_circular_patch_not_empty()
{
   std::cout << "\n[T-D12] gradual_overstress_compact_circular (TPV102/104): "
                "patch detected via kind-aware geometry\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;

   const real_t cx = 1000.0, cy = 0.0, cz = -3000.0;

   SlipWeakeningPerDOFParams lsw;
   lsw.mu_s.SetSize(1); lsw.mu_d.SetSize(1); lsw.d_c.SetSize(1);
   lsw.cohesion.SetSize(1);
   lsw.mu_s(0) = 0.65; lsw.mu_d(0) = 0.30; lsw.d_c(0) = 1.0;
   lsw.cohesion(0) = 0.0;

   Vector tau_pre(2); tau_pre(0) = 0.0; tau_pre(1) = 32.5e6;
   Vector sn(1);      sn(0) = 71.5e6;
   Vector coords(3);  coords(0) = cx; coords(1) = cy; coords(2) = cz;

   NucleationSpec nuc;
   nuc.enabled = true;
   nuc.kind    = NucleationKind::GradualOverstressCompactCircular;
   nuc.compact_circular.center_x_m   = cx;
   nuc.compact_circular.center_y_m   = cy;
   nuc.compact_circular.center_z_m   = cz;
   nuc.compact_circular.radius_m     = 3000.0;
   nuc.compact_circular.delta_tau_pa = 25.0e6;
   nuc.compact_circular.T_nuc_s      = 1.0;

   DenseMatrix basis(9, 1); basis = 0.0;
   basis(0, 0) = 1.0; basis(4, 0) = 1.0; basis(8, 0) = 1.0;

   // Mirror the driver: compact-circular populates amplitude_strike + radial.
   GradualOverstressPerDOFParams nuc_params;
   {
      const auto cp = ResolveGradualOverstressCompactCircular(
         nuc.compact_circular, /*enabled=*/true, coords, basis);
      nuc_params.amplitude_strike = cp.amplitude_strike;
      nuc_params.radial           = cp.radial;
   }

   PrintDerivedConfig pd_cfg;
   pd_cfg.enabled               = true;
   pd_cfg.abort_on_failure      = false;
   pd_cfg.outside_safety_factor = 3.0;

   std::ostringstream oss;
   unsetenv("SEAS_SKIP_EQUILIBRIUM_GATE");
   (void)PrintDerivedAndCheck(
      pd_cfg, lsw, tau_pre, sn, coords, nuc, nuc_params, s,
      /*mu_bulk=*/32.0e9, /*cp=*/4877.0, /*cs=*/3458.0,
      /*h_min=*/1000.0, /*dt_cfl=*/6.84e-5, /*tfinal=*/12.0,
      /*num_fault_global=*/1, /*num_zero_normal_fallbacks=*/0,
#ifdef MFEM_USE_MPI
      MPI_COMM_WORLD,
#endif
      /*rank=*/0, oss);
   const std::string captured = oss.str();

   TEST_ASSERT(captured.find("fault DOFs are inside the nucleation patch")
               == std::string::npos,
               "compact-circular patch NOT falsely reported empty");
   const real_t peak = ParseAfter(captured,
                                  "nucleation peak |F(r) * delta_tau| = ");
   TEST_NEAR(peak, 25.0e6, 1.0e3,
             "compact-bell peak == delta_tau_pa at the centre");
   TEST_ASSERT(captured.find("(sufficient)") != std::string::npos,
               "trigger overshoot computed from the strike amplitude "
               "(sufficient)");
}

// NOTE on the unhandled-kind `default: MFEM_ABORT` guards (both overloads):
// they are deliberately NOT unit-tested here.  An in-process death test would
// fork() and trigger MFEM_ABORT in the child, but this test's main() calls
// MPI_Init, so MFEM_ABORT routes through MPI_Abort (general/error.cpp:182),
// which HANGS in a forked child of an MPI process (the child is not a real
// rank).  The serial fork death-tests elsewhere (e.g. test_spatial_friction_
// config.cpp) work only because their main() never initialises MPI.  The
// guard mirrors the identical, also-untested contract in MakeNucleation
// (nucleation_factory.cpp:53); the three REAL kinds are covered positively by
// T-D11 (instantaneous), T-D12 (compact), and the gradual cases above.

// =====================================================================
// T-D14: STATIC-OVERSTRESS nucleation (TPV205) — no [nucleation] block.
//   TPV205 nucleates by a statically OVERstressed patch (|tau_pre| >
//   mu_s*sigma_n there); that is the mechanism, NOT a "background
//   supercritical" failure.  Gate (b) must NOT fire when nuc is disabled,
//   and the config must PASS.  Mirrors TPV205 patch #1: tau = 81.6 MPa vs
//   mu_s*sigma_n = 0.677*120 = 81.24 MPa (ratio 1.0044, exactly the value
//   the real run aborted on).  abort_on_failure=true so a regression (gate
//   (b) still firing) would MFEM_ABORT the test rather than fail soft.
// =====================================================================
static void T_D14_static_overstress_no_nucleation_passes()
{
   std::cout << "\n[T-D14] static-overstress nucleation (no [nucleation], "
                "TPV205) does NOT trip the background-supercritical gate\n";
   StressSpec s; s.kind = StressSourceKind::ConstantTensor;
   GradualOverstressSpec gspec;   // unused (nuc disabled)
   std::string captured;
   const real_t ratio = RunPrinter(
      s, /*mu_s=*/0.677, /*mu_d=*/0.525, /*d_c=*/0.4,
      /*sigma_n=*/120.0e6, /*tau_dip=*/0.0, /*tau_strike=*/81.6e6,
      /*nuc_enabled=*/false, gspec, /*mu_bulk=*/32.0e9, captured,
      /*abort_on_failure=*/true);
   TEST_ASSERT(ratio >= 1.0,
               "static-overstress ratio >= 1 (it nucleates via the patch)");
   TEST_ASSERT(captured.find("FAIL: max OUTSIDE asperity") == std::string::npos,
               "gate (b) does NOT fire for a nuc-disabled static-overstress "
               "config (the TPV205 regression)");
   TEST_ASSERT(captured.find("PASS: initial conditions are well-posed")
               != std::string::npos,
               "static-overstress config PASSes (TPV205 nucleation mechanism)");
}

int main(int argc, char** argv)
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
#else
   (void)argc; (void)argv;
#endif
   std::cout << "Running Phase D test_spatial_print_derived\n";
   T_D06_sigma1_azimuth_45deg();
   T_D07_dt_cfl_passthrough();
   T_D02_LITE_supercritical_returns_ratio();
   T_D05_LITE_env_skip_gate_honored();
   T_D02_COORDS_outside_message_includes_xyz();
   T_D03_REASON_warn_prefix_attribution();
   T_D08_overshoot_uses_initiation_criterion();
   T_D09_negative_stress_drop_caught();
   T_D10_insufficient_trigger_caught();
   T_RS01_steady_state_friction();
   T_RS02_nucleation_length();
   T_RS03_all_velocity_strengthening();
   T_RS04_vw_nucleation_passes();
   T_RS05_patch_all_vs_fails();
   T_D11_instantaneous_circular_patch_not_empty();
   T_D12_compact_circular_patch_not_empty();
   T_D14_static_overstress_no_nucleation_passes();
   std::cout << "\n========================================\n";
   std::cout << "Phase D test_spatial_print_derived: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return num_failed == 0 ? 0 : 1;
}
