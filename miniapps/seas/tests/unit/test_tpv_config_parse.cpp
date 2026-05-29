// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_tpv_config_parse.cpp — Phase 8 (PLAN_tpv_regression_via_spatial_
// dyn_driver_2026-05-24.md): the three TPV TOMLs authored in Phase 8
// (tpv{205,102,104}/configs/*.toml) parse through LoadSpatialFrictionConfig
// with no abort, and every Phase-8 field lands in SpatialFrictionConfig as
// authored.  This is the deterministic, run-independent half of the plan's
// "_review" tests + the AC "each config parses".  The end-to-end driver
// dry-run / 8-rank smoke / gold regression are deferred to a run session
// (and additionally blocked by the hardcoded driver boundary attrs — see
// the TOML headers' DRIVER-BOUNDARY GAP note).
//
// Framework: the project's custom TEST_ASSERT / TEST_NEAR harness (no gtest),
// matching tests/unit/test_compute_safs_params.cpp et al.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <sys/stat.h>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e << ", tol " << _t \
   << ", line " << __LINE__ << ")\n"; num_failed++; } else { \
   std::cout << "  PASSED: " << m << "\n"; num_passed++; } } while (0)

namespace
{
bool FileExists(const std::string& p) { struct stat s; return stat(p.c_str(), &s) == 0; }

// Locate a Phase-8 TPV config relative to a few plausible CWDs (the test
// binary normally runs from miniapps/seas).  Returns "" if not found.
std::string FindConfig(const std::string& rel)
{
   const std::vector<std::string> prefixes = { "", "../", "../../" };
   for (const auto& pre : prefixes)
   {
      const std::string p = pre + rel;
      if (FileExists(p)) { return p; }
   }
   return "";
}
} // anon

// ---------------------------------------------------------------------
// TPV205 (LSW)
// ---------------------------------------------------------------------
static void T_TPV205(const std::string& path)
{
   std::cout << "\n[TPV205] parse " << path << "\n";
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(path);

   TEST_ASSERT(cfg.schema_version == 1, "schema_version == 1");
   TEST_ASSERT(cfg.law == spatial::FrictionLawKind::SlipWeakening,
               "law == slip_weakening");

   // Stress: fault-local prestress, right-lateral / compression POSITIVE.
   TEST_ASSERT(cfg.stress.kind == spatial::StressSourceKind::FaultLocalPrestress,
               "stress.kind == fault_local_prestress");
   TEST_NEAR(cfg.stress.tau_strike_pa, 70.0e6, 1.0, "tau_strike_pa == +70 MPa (background)");
   TEST_NEAR(cfg.stress.tau_dip_pa,    0.0,    1.0, "tau_dip_pa == 0");
   TEST_NEAR(cfg.stress.sigma_n_pa,    120.0e6,1.0, "sigma_n_pa == 120 MPa (compression +)");

   // Three SCEC stress patches (nucleation/left/right).
   TEST_ASSERT(cfg.stress.fault_local_patches.size() == 3u, "3 stress patches");
   if (cfg.stress.fault_local_patches.size() == 3u)
   {
      TEST_NEAR(cfg.stress.fault_local_patches[0].tau_strike_pa, 81.6e6, 1.0,
                "patch[0] (nucleation) == +81.6 MPa");
      TEST_NEAR(cfg.stress.fault_local_patches[1].tau_strike_pa, 78.0e6, 1.0,
                "patch[1] (left) == +78 MPa");
      TEST_NEAR(cfg.stress.fault_local_patches[2].tau_strike_pa, 62.0e6, 1.0,
                "patch[2] (right) == +62 MPa");
   }

   // Numerics: scalar Riemann, DG CFL safety, one-shot LSW (Phase 8 req 1a).
   TEST_ASSERT(cfg.numerics.interior_flux == spatial::InteriorFlux::Scalar,
               "interior_flux == scalar");
   TEST_ASSERT(cfg.numerics.cfl_safety == spatial::CflSafety::Dg,
               "cfl_safety == dg");
   // REVIEW R-003: the spatial driver always sub-steps (one-shot is not
   // implemented + the O2-substep gold), so the TPV205 TOML uses "substep".
   TEST_ASSERT(cfg.numerics.fault_iterator == spatial::FaultIteratorKind::Substep,
               "fault_iterator == substep (one-shot unsupported; matches O2 gold)");
   TEST_ASSERT(cfg.numerics.mixed_flux == "adjacent", "mixed_flux == adjacent");

   // LSW block + barrier rules.
   TEST_ASSERT(cfg.slip_weakening.has_value(), "[friction.slip_weakening] present");
   TEST_ASSERT(!cfg.rate_state.has_value(), "no [friction.rate_state]");
   if (cfg.slip_weakening.has_value())
   {
      const auto& sw = *cfg.slip_weakening;
      TEST_NEAR(sw.mu_s_default, 0.677, 1e-12, "mu_s_default == 0.677");
      TEST_NEAR(sw.mu_d_default, 0.525, 1e-12, "mu_d_default == 0.525");
      TEST_NEAR(sw.d_c_default,  0.40,  1e-12, "d_c_default == 0.40 m");
      int n_barrier = 0;
      for (const auto& r : sw.spatial)
      { if (r.kind == spatial::SpatialRule::Kind::Barrier) { ++n_barrier; } }
      TEST_ASSERT(n_barrier == 3, "3 barrier spatial rules");
   }

   // TPV205 has no [nucleation] block (static stress-patch nucleation).
   TEST_ASSERT(!cfg.nucleation.enabled, "nucleation disabled (static patches)");

   // TPV205 mesh boundary convention: fault=103, free=101, absorbing=105
   // (differs from TPV102/104's 1/3/5; matches native tpv205_driver).
   TEST_ASSERT(cfg.boundary.fault_attr == 103,
               "boundary.fault_attr == 103 (TPV205 mesh)");
}

// ---------------------------------------------------------------------
// TPV102 (RS aging)
// ---------------------------------------------------------------------
static void T_TPV102(const std::string& path)
{
   std::cout << "\n[TPV102] parse " << path << "\n";
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(path);

   TEST_ASSERT(cfg.law == spatial::FrictionLawKind::RateState, "law == rate_state");
   TEST_ASSERT(cfg.stress.kind == spatial::StressSourceKind::FaultLocalPrestress,
               "stress.kind == fault_local_prestress");
   TEST_NEAR(cfg.stress.tau_strike_pa, 75.0e6, 1.0, "tau_strike_pa == +75 MPa (tau_ini)");
   TEST_NEAR(cfg.stress.sigma_n_pa,    120.0e6,1.0, "sigma_n_pa == 120 MPa");

   TEST_ASSERT(cfg.numerics.fault_iterator == spatial::FaultIteratorKind::Substep,
               "fault_iterator == substep (RS)");

   TEST_ASSERT(cfg.rate_state.has_value(), "[friction.rate_state] present");
   TEST_ASSERT(!cfg.slip_weakening.has_value(), "no [friction.slip_weakening]");
   if (cfg.rate_state.has_value())
   {
      const auto& rs = *cfg.rate_state;
      // Inside-out: VW value (a_vw=0.008) is the default; VS border (0.016)
      // via box rules.  a_default < b_default satisfies the parser guard.
      TEST_NEAR(rs.a_default,  0.008, 1e-12, "a_default == 0.008 (a_vw, VW core)");
      TEST_NEAR(rs.b_default,  0.012, 1e-12, "b_default == 0.012");
      TEST_NEAR(rs.Dc_default, 0.02,  1e-12, "Dc_default == 0.02 m");
      TEST_NEAR(rs.V_0_default, 1.0e-6, 1e-18, "V_0_default == 1e-6");
      TEST_ASSERT(rs.state_evolution == spatial::StateEvolutionKind::AgingLaw,
                  "state_evolution == aging_law");
      // VS border box rules set a = a_vs = 0.016 (> b: velocity-strengthening).
      int n_vs = 0;
      for (const auto& r : rs.spatial)
      {
         if (r.kind == spatial::SpatialRule::Kind::Box && !std::isnan(r.a))
         {
            ++n_vs;
            TEST_NEAR(r.a, 0.016, 1e-12, "VS-border box a == 0.016 (a_vs)");
         }
      }
      TEST_ASSERT(n_vs == 3, "3 VS-border box rules present");
   }

   TEST_ASSERT(cfg.nucleation.enabled, "nucleation enabled");
   TEST_ASSERT(cfg.nucleation.kind ==
               spatial::NucleationKind::GradualOverstressCompactCircular,
               "nucleation == gradual_overstress_compact_circular");
   TEST_NEAR(cfg.nucleation.compact_circular.delta_tau_pa, 25.0e6, 1.0,
             "compact-circular delta_tau == 25 MPa");
   TEST_NEAR(cfg.nucleation.compact_circular.radius_m, 3000.0, 1e-6,
             "compact-circular radius == 3 km");
}

// ---------------------------------------------------------------------
// TPV104 (RS slip-law strong rate weakening)
// ---------------------------------------------------------------------
static void T_TPV104(const std::string& path)
{
   std::cout << "\n[TPV104] parse " << path << "\n";
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(path);

   TEST_ASSERT(cfg.law == spatial::FrictionLawKind::RateState, "law == rate_state");
   TEST_NEAR(cfg.stress.tau_strike_pa, 40.0e6, 1.0, "tau_strike_pa == +40 MPa (tau_ini)");

   TEST_ASSERT(cfg.rate_state.has_value(), "[friction.rate_state] present");
   if (cfg.rate_state.has_value())
   {
      const auto& rs = *cfg.rate_state;
      TEST_ASSERT(rs.state_evolution ==
                  spatial::StateEvolutionKind::SlipLawStrongRateWeakening,
                  "state_evolution == slip_law_strong_rate_weakening");
      // Inside-out: VW core (a_in=0.01, V_w_in=0.1) default; VS border via box.
      TEST_NEAR(rs.a_default,   0.01,  1e-12, "a_default == 0.01 (a_in, VW core)");
      TEST_NEAR(rs.b_default,   0.014, 1e-12, "b_default == 0.014");
      TEST_NEAR(rs.Dc_default,  0.4,   1e-12, "Dc_default == 0.4 m (L)");
      TEST_NEAR(rs.f_w_default, 0.2,   1e-12, "f_w_default == 0.2 (muW)");
      TEST_NEAR(rs.V_w_default, 0.1,   1e-12, "V_w_default == 0.1 (V_w_in, VW core)");
      // VS border box rules: a = a_out = 0.02, V_w = V_w_out = 1.0.
      int n_vs = 0;
      for (const auto& r : rs.spatial)
      {
         if (r.kind == spatial::SpatialRule::Kind::Box && !std::isnan(r.a))
         {
            ++n_vs;
            TEST_NEAR(r.a,   0.02, 1e-12, "VS-border box a == 0.02 (a_out)");
            TEST_NEAR(r.V_w, 1.0,  1e-12, "VS-border box V_w == 1.0 (V_w_out)");
         }
      }
      TEST_ASSERT(n_vs == 3, "3 VS-border box rules present");
   }

   TEST_ASSERT(cfg.nucleation.kind ==
               spatial::NucleationKind::GradualOverstressCompactCircular,
               "nucleation == gradual_overstress_compact_circular");
   TEST_NEAR(cfg.nucleation.compact_circular.delta_tau_pa, 45.0e6, 1.0,
             "compact-circular delta_tau == 45 MPa (TPV104)");
}

// REVIEW R-001/R-002/R-003: the [numerics] method selectors are now consumed
// by the driver via these pure decision helpers (spatial_friction.hpp).  These
// tests pin the discriminating behavior the driver relies on — they would have
// failed before the fix (the helpers/guards did not exist and the selectors
// were dead).  Mesh-free / file-free, so they always run.
static void T_numerics_dispatch_helpers()
{
   std::cout << "\n[NUM-DISPATCH] cfl_safety / interior_flux / fault_iterator helpers\n";
   spatial::SpatialFrictionConfig cfg;
   cfg.mesh.order = 1;
   cfg.numerics.cfl = 0.25;

   // R-002: CflSafetyFactor — Dg applies 1/(3*(2p+1)); Raw is 1.0.
   cfg.numerics.cfl_safety = spatial::CflSafety::Dg;
   TEST_NEAR(spatial::CflSafetyFactor(cfg), 1.0 / 9.0, 1e-15,
             "CflSafetyFactor(Dg, p=1) == 1/(3*3) = 1/9");
   cfg.numerics.cfl_safety = spatial::CflSafety::Raw;
   TEST_NEAR(spatial::CflSafetyFactor(cfg), 1.0, 1e-15,
             "CflSafetyFactor(Raw) == 1.0 (no DG factor)");
   cfg.mesh.order = 2;
   cfg.numerics.cfl_safety = spatial::CflSafety::Dg;
   TEST_NEAR(spatial::CflSafetyFactor(cfg), 1.0 / 15.0, 1e-15,
             "CflSafetyFactor(Dg, p=2) == 1/(3*5) = 1/15");

   // (Phase 9: the InteriorFluxSupported Phase-8 stopgap is removed — the
   // driver now branches on interior_flux to build the scalar or matrix
   // WaveOperator ctor, so matrix is supported.)

   // R-003: FaultIteratorSupported — only Substep is implemented.
   cfg.numerics.fault_iterator = spatial::FaultIteratorKind::Substep;
   TEST_ASSERT(spatial::FaultIteratorSupported(cfg),
               "fault_iterator=substep supported");
   cfg.numerics.fault_iterator = spatial::FaultIteratorKind::OneShot;
   TEST_ASSERT(!spatial::FaultIteratorSupported(cfg),
               "fault_iterator=one-shot rejected");

   // R-002 SAFS-safe default: a default-constructed config is Dg, so a config
   // that omits cfl_safety keeps the driver's always-DG-factored behavior.
   spatial::SpatialFrictionConfig def;
   TEST_ASSERT(def.numerics.cfl_safety == spatial::CflSafety::Dg,
               "cfl_safety struct default == Dg (R-002: no SAFS regression)");
   // R-001 SAFS-safe default: a config that omits fault_iterator must default
   // to the SUPPORTED mode (Substep), else the driver's FaultIteratorSupported
   // guard aborts every SAFS run.
   TEST_ASSERT(def.numerics.fault_iterator == spatial::FaultIteratorKind::Substep,
               "fault_iterator struct default == Substep (R-001: no SAFS abort)");
   TEST_ASSERT(spatial::FaultIteratorSupported(def),
               "default-constructed config passes FaultIteratorSupported (R-001)");
}

int main(int, char**)
{
   std::cout << "Running Phase 8 TPV config-parse tests\n";

   // Numerics-selector helper tests run unconditionally (no files needed).
   T_numerics_dispatch_helpers();

   const std::string p205 = FindConfig("tpv205/configs/tpv205_spatial.toml");
   const std::string p102 = FindConfig("tpv102/configs/tpv102_spatial.toml");
   const std::string p104 = FindConfig("tpv104/configs/tpv104_spatial.toml");

   if (p205.empty() || p102.empty() || p104.empty())
   {
      std::cout << "(TPV config TOMLs not found relative to CWD; "
                   "file-based parse skipped — helper tests still ran)\n";
      std::cout << "\nPhase 8 config-parse: " << num_passed << " / " << num_tests
                << " passed, " << num_failed << " failed (helpers only)\n";
      return (num_failed == 0) ? 0 : 1;
   }

   T_TPV205(p205);
   T_TPV102(p102);
   T_TPV104(p104);

   std::cout << "\n========================================\n";
   std::cout << "Phase 8 config-parse: " << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
