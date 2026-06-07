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
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
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

std::string ReadFile(const std::string& p)
{
   std::ifstream in(p, std::ios::binary);
   std::ostringstream ss; ss << in.rdbuf();
   return ss.str();
}
void WriteFile(const std::string& p, const std::string& s)
{
   std::ofstream out(p, std::ios::binary | std::ios::trunc);
   out << s;
}
// Replace the FIRST occurrence of `from` with `to`; returns false if `from`
// was not present (so a test can fail loud rather than silently no-op).
bool ReplaceFirst(std::string& s, const std::string& from, const std::string& to)
{
   const std::size_t pos = s.find(from);
   if (pos == std::string::npos) { return false; }
   s.replace(pos, from.size(), to);
   return true;
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
      // via the boxcar_taper rule below.  a_default < b_default satisfies the
      // parser guard.
      TEST_NEAR(rs.a_default,  0.008, 1e-12, "a_default == 0.008 (a_vw, VW core)");
      TEST_NEAR(rs.b_default,  0.012, 1e-12, "b_default == 0.012");
      TEST_NEAR(rs.Dc_default, 0.02,  1e-12, "Dc_default == 0.02 m");
      TEST_NEAR(rs.V_0_default, 1.0e-6, 1e-18, "V_0_default == 1e-6");
      TEST_ASSERT(rs.state_evolution == spatial::StateEvolutionKind::AgingLaw,
                  "state_evolution == aging_law");
      // Phase 8 (smooth SCEC taper): the velocity-strengthening border is the
      // single kind="boxcar_taper" rate-state rule (a ramps a_inner=a_vw ->
      // a_outer=a_vs over the 3 km tanh margin), NOT the removed hard Box rules.
      // Assert exactly one BoxcarTaper rule carrying the VS-border (a_outer)
      // and VW-core (a_inner) endpoints.
      int n_taper = 0;
      for (const auto& r : rs.spatial)
      {
         if (r.kind == spatial::SpatialRule::Kind::BoxcarTaper &&
             !std::isnan(r.a_outer))
         {
            ++n_taper;
            TEST_NEAR(r.a_outer, 0.016, 1e-12,
                      "boxcar_taper a_outer == 0.016 (a_vs, VS border)");
            TEST_NEAR(r.a_inner, 0.008, 1e-12,
                      "boxcar_taper a_inner == 0.008 (a_vw, VW core)");
         }
      }
      TEST_ASSERT(n_taper == 1, "1 VS-border boxcar_taper rule present");
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
      // Inside-out: VW core (a_in=0.01, V_w_in=0.1) default; VS border via the
      // boxcar_taper rule below.
      TEST_NEAR(rs.a_default,   0.01,  1e-12, "a_default == 0.01 (a_in, VW core)");
      TEST_NEAR(rs.b_default,   0.014, 1e-12, "b_default == 0.014");
      TEST_NEAR(rs.Dc_default,  0.4,   1e-12, "Dc_default == 0.4 m (L)");
      TEST_NEAR(rs.f_w_default, 0.2,   1e-12, "f_w_default == 0.2 (muW)");
      TEST_NEAR(rs.V_w_default, 0.1,   1e-12, "V_w_default == 0.1 (V_w_in, VW core)");
      // Phase 8 (smooth SCEC taper): the velocity-strengthening border is the
      // single kind="boxcar_taper" rate-state rule (a ramps a_in -> a_out and
      // V_w ramps V_w_in -> V_w_out over the 3 km tanh margin), NOT the removed
      // hard Box rules.  Assert exactly one BoxcarTaper rule carrying the a /
      // V_w endpoints.
      int n_taper = 0;
      for (const auto& r : rs.spatial)
      {
         if (r.kind == spatial::SpatialRule::Kind::BoxcarTaper &&
             !std::isnan(r.a_outer))
         {
            ++n_taper;
            TEST_NEAR(r.a_outer,   0.02, 1e-12,
                      "boxcar_taper a_outer == 0.02 (a_out, VS border)");
            TEST_NEAR(r.a_inner,   0.01, 1e-12,
                      "boxcar_taper a_inner == 0.01 (a_in, VW core)");
            TEST_NEAR(r.V_w_outer, 1.0,  1e-12,
                      "boxcar_taper V_w_outer == 1.0 (V_w_out)");
            TEST_NEAR(r.V_w_inner, 0.1,  1e-12,
                      "boxcar_taper V_w_inner == 0.1 (V_w_in)");
         }
      }
      TEST_ASSERT(n_taper == 1, "1 VS-border boxcar_taper rule present");
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
   // (Unified bi-material plan, Part A) the central-flux contrast guard defaults to
   // DISABLED (< 0): a config that omits mixed_flux_contrast_tol is byte-exact.
   TEST_ASSERT(def.numerics.mixed_flux_contrast_tol < 0.0,
               "mixed_flux_contrast_tol struct default < 0 (Part A: guard disabled "
               "by default -> byte-exact)");
   // R-001 SAFS-safe default: a config that omits fault_iterator must default
   // to the SUPPORTED mode (Substep), else the driver's FaultIteratorSupported
   // guard aborts every SAFS run.
   TEST_ASSERT(def.numerics.fault_iterator == spatial::FaultIteratorKind::Substep,
               "fault_iterator struct default == Substep (R-001: no SAFS abort)");
   TEST_ASSERT(spatial::FaultIteratorSupported(def),
               "default-constructed config passes FaultIteratorSupported (R-001)");
}

// ---------------------------------------------------------------------
// (Phase 5, BUG-21/22) TPV31 (matrix + depth_profile_1d) parser regression.
// Test 5.3c: the additive seam_continuous parse + the G1 relaxation must not
// perturb the existing TPV31 config fields, and seam_continuous defaults false.
// ---------------------------------------------------------------------
static void T_TPV31_Regression(const std::string& path)
{
   std::cout << "\n[TPV31 regression] parse " << path << "\n";
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(path);

   // Key fields the Phase-5 parser edits sit adjacent to — unchanged.
   TEST_ASSERT(cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix,
               "tpv31 interior_flux == matrix (unchanged)");
   TEST_ASSERT(cfg.numerics.mixed_flux == "none",
               "tpv31 mixed_flux == none (unchanged; G1 relaxation did not flip it)");
   TEST_ASSERT(cfg.material.kind == spatial::MaterialKind::DepthProfile1D,
               "tpv31 material.kind == depth_profile_1d (unchanged)");
   TEST_ASSERT(cfg.material.depth_axis == 'z',
               "tpv31 material.depth_axis == 'z' (unchanged)");
   TEST_ASSERT(!cfg.material.profile_layers.empty(),
               "tpv31 material profile layers parsed (unchanged)");
   // The additive field: absent in the TOML ⇒ defaults false.
   TEST_ASSERT(cfg.material.seam_continuous == false,
               "tpv31 material.seam_continuous defaults false (key omitted)");
}

// ---------------------------------------------------------------------
// Test 5.3a + 5.3b: derive a matrix + adjacent + seam_continuous=true config
// from the real TPV31 TOML (string-mutated copy; the parser stores referenced
// paths as strings without resolving them, so an in-place temp copy parses
// identically).  Asserts G1 no longer aborts on matrix+adjacent and that
// seam_continuous round-trips into MaterialSpec.
// ---------------------------------------------------------------------
static void T_TPV31_MatrixAdjacentSeam(const std::string& path)
{
   std::cout << "\n[TPV31 matrix+adjacent+seam] derive from " << path << "\n";
   std::string toml = ReadFile(path);
   const bool r1 = ReplaceFirst(toml, "mixed_flux     = \"none\"",
                                       "mixed_flux     = \"adjacent\"");
   // Newline-anchored so we hit the standalone `[material]` config line, NOT
   // the `In canonical: depth_axis = "z" ...` comment earlier in the file.
   const bool r2 = ReplaceFirst(toml, "\ndepth_axis = \"z\"\n",
                                       "\ndepth_axis = \"z\"\nseam_continuous = true\n");
   TEST_ASSERT(r1, "5.3a: located + flipped mixed_flux none->adjacent in TPV31 TOML");
   TEST_ASSERT(r2, "5.3b: located the depth_axis config line to inject seam_continuous=true");

   // (P5-3) Write to TMPDIR, NOT the source tree: LoadSpatialFrictionConfig
   // MFEM_ABORTs (does not throw) on a parse failure, so a post-parse remove
   // would leak the file into the repo on failure.  TMPDIR keeps any leak out
   // of the tracked tree.  The parser stores referenced paths as strings
   // without resolving them, so the temp location is irrelevant to parsing.
   const char* td = std::getenv("TMPDIR");
   const std::string tmpdir = (td && *td) ? std::string(td) : std::string("/tmp");
   const std::string tmp = tmpdir + "/seas_tpv31_matrixadj_TESTGEN.toml";
   WriteFile(tmp, toml);
   // G1 (spatial_friction.cpp:1162) would have aborted here pre-Phase-5.
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(tmp);
   std::remove(tmp.c_str());

   TEST_ASSERT(cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix,
               "5.3a: interior_flux == matrix");
   TEST_ASSERT(cfg.numerics.mixed_flux == "adjacent",
               "5.3a: matrix + mixed_flux=adjacent parses WITHOUT G1 abort (BUG-21)");
   TEST_ASSERT(cfg.material.seam_continuous == true,
               "5.3b: [material].seam_continuous=true round-trips into MaterialSpec (BUG-22)");
}

// ---------------------------------------------------------------------
// (Phase 7) The shipped tpv31_rk_mixedflux.toml ARTIFACT parses past G1 and the
// G2 decision is correct: matrix + adjacent + seam_continuous, REJECTED under
// the TOML's default ADER, ALLOWED under the job's --time-integrator rk45.
// ---------------------------------------------------------------------
static void T_TPV31_RkMixedFlux(const std::string& path)
{
   std::cout << "\n[TPV31 rk_mixedflux artifact] parse " << path << "\n";
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(path);
   TEST_ASSERT(cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix,
               "tpv31_rk_mixedflux interior_flux == matrix (preserved)");
   TEST_ASSERT(cfg.numerics.mixed_flux == "adjacent",
               "tpv31_rk_mixedflux mixed_flux == adjacent (parses past G1)");
   TEST_ASSERT(cfg.material.seam_continuous == true,
               "tpv31_rk_mixedflux [material].seam_continuous == true");
   TEST_ASSERT(cfg.material.kind == spatial::MaterialKind::DepthProfile1D,
               "tpv31_rk_mixedflux material.kind == depth_profile_1d (unchanged)");
   // G2: under the TOML default (ADER) the driver would REJECT; the job's CLI
   // --time-integrator rk45 (applied before G2) flips it to ALLOWED.
   TEST_ASSERT(spatial::MatrixMixedFluxUnderAder(cfg),
               "tpv31_rk_mixedflux: G2 REJECTS under ADER (the TOML default)");
   cfg.numerics.time_integrator = spatial::TimeIntegratorKind::RK45;
   TEST_ASSERT(!spatial::MatrixMixedFluxUnderAder(cfg),
               "tpv31_rk_mixedflux: G2 ALLOWS under --time-integrator rk45 (the job CLI)");
}

// ---------------------------------------------------------------------
// (Phase 5, Tests 5.1/5.2) Driver G2 decision — table-test the PURE predicate
// spatial::MatrixMixedFluxUnderAder that the driver's G2 MFEM_VERIFY uses
// verbatim.  Covers the exact guard decision (matrix+mixed+ADER -> reject;
// matrix+mixed+RK -> allow; matrix+none -> allow; scalar -> allow) WITHOUT
// building the driver/operator/mesh.  No files needed.
// ---------------------------------------------------------------------
static void T_G2_MatrixMixedFluxUnderAder()
{
   std::cout << "\n[Phase 5 G2] MatrixMixedFluxUnderAder truth table\n";
   using spatial::InteriorFlux;
   using spatial::TimeIntegratorKind;
   auto make = [](InteriorFlux ifx, const std::string& mf, TimeIntegratorKind ti)
   {
      spatial::SpatialFrictionConfig c;
      c.numerics.interior_flux   = ifx;
      c.numerics.mixed_flux      = mf;
      c.numerics.time_integrator = ti;
      return c;
   };
   // Test 5.2: matrix + mixed + ADER -> G2 rejects (predicate true).
   TEST_ASSERT(spatial::MatrixMixedFluxUnderAder(
                  make(InteriorFlux::Matrix, "adjacent", TimeIntegratorKind::ADER)),
               "5.2: matrix + adjacent + ADER -> G2 rejects (predicate true)");
   TEST_ASSERT(spatial::MatrixMixedFluxUnderAder(
                  make(InteriorFlux::Matrix, "all_continuous", TimeIntegratorKind::ADER)),
               "5.2: matrix + all_continuous + ADER -> G2 rejects");
   // Test 5.1: matrix + mixed + RK -> G2 allows (predicate false).
   TEST_ASSERT(!spatial::MatrixMixedFluxUnderAder(
                  make(InteriorFlux::Matrix, "adjacent", TimeIntegratorKind::RK4)),
               "5.1: matrix + adjacent + RK4 -> G2 allows (predicate false)");
   TEST_ASSERT(!spatial::MatrixMixedFluxUnderAder(
                  make(InteriorFlux::Matrix, "adjacent", TimeIntegratorKind::RK45)),
               "5.1: matrix + adjacent + RK45 -> G2 allows");
   // matrix + none -> allowed under either integrator (bi-material Godunov upwind).
   TEST_ASSERT(!spatial::MatrixMixedFluxUnderAder(
                  make(InteriorFlux::Matrix, "none", TimeIntegratorKind::ADER)),
               "matrix + none + ADER -> G2 allows (the existing TPV31 path)");
   TEST_ASSERT(!spatial::MatrixMixedFluxUnderAder(
                  make(InteriorFlux::Matrix, "none", TimeIntegratorKind::RK4)),
               "matrix + none + RK4 -> G2 allows");
   // scalar -> G2 never constrains (scalar mixed+ADER caught downstream by ComputeMaxDt).
   TEST_ASSERT(!spatial::MatrixMixedFluxUnderAder(
                  make(InteriorFlux::Scalar, "adjacent", TimeIntegratorKind::ADER)),
               "scalar + adjacent + ADER -> G2 allows (scalar not constrained by G2)");
}

// ---------------------------------------------------------------------
// Part C / C2 (TPV6/TPV7): the bi-material-fault config artifact parses —
// halfspace_across_fault material + matrix flux + the TPV205-style static
// stress-patch nucleation (NO [nucleation] block) + LSW + barriers.
// ---------------------------------------------------------------------
static void T_TPV6(const std::string& path, const char* label,
                   double vp_far_expected)
{
   std::cout << "\n[" << label << "] parse " << path << "\n";
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(path);

   // Bi-material across the fault + per-side Riemann (matrix).
   TEST_ASSERT(cfg.material.kind == spatial::MaterialKind::HalfspaceAcrossFault,
               "material.kind == halfspace_across_fault");
   TEST_NEAR(cfg.material.halfspace.vp_near, 6000.0, 1.0,
             "halfspace vp_near == 6000 (fast/near side)");
   TEST_NEAR(cfg.material.halfspace.vp_far, vp_far_expected, 1.0,
             "halfspace vp_far (per problem)");
   TEST_NEAR(cfg.material.halfspace.normal[1], -1.0, 1e-12,
             "halfspace n_y == -1 (matches ref_normal)");
   TEST_ASSERT(cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix,
               "interior_flux == matrix (REQUIRED for per-side fault Riemann)");

   // Nucleation = TPV205-style STATIC stress patch, NOT a time-varying driver.
   TEST_ASSERT(!cfg.nucleation.enabled,
               "no [nucleation] block (static stress-patch nucleation, like TPV205)");
   TEST_ASSERT(cfg.stress.kind == spatial::StressSourceKind::FaultLocalPrestress,
               "stress.kind == fault_local_prestress");
   TEST_NEAR(cfg.stress.tau_strike_pa, 70.0e6, 1.0, "background tau_strike == 70 MPa");
   TEST_NEAR(cfg.stress.sigma_n_pa, 120.0e6, 1.0, "sigma_n == 120 MPa");
   TEST_ASSERT(cfg.stress.fault_local_patches.size() == 1u,
               "exactly 1 stress patch (the nucleation square)");
   if (cfg.stress.fault_local_patches.size() == 1u)
   {
      const auto& pch = cfg.stress.fault_local_patches[0];
      TEST_NEAR(pch.tau_strike_pa, 81.6e6, 1.0,
                "nucleation patch tau_strike == 81.6 MPa (> yield 0.677*120 = 81.24)");
      TEST_NEAR(pch.center_x_m, 0.0, 1e-6, "patch center_x == 0 (config frame)");
      TEST_NEAR(pch.center_z_m, -7500.0, 1e-6, "patch center_z == -7500");
      TEST_NEAR(pch.half_x_m, 1500.0, 1e-6, "patch half_x == 1500 (3000 m square)");
      TEST_NEAR(pch.half_z_m, 1500.0, 1e-6, "patch half_z == 1500");
   }

   // LSW friction + barriers.
   TEST_ASSERT(cfg.law == spatial::FrictionLawKind::SlipWeakening,
               "law == slip_weakening");
   TEST_ASSERT(cfg.slip_weakening.has_value(),
               "[friction.slip_weakening] block present");
   if (cfg.slip_weakening.has_value())
   {
      TEST_NEAR(cfg.slip_weakening->mu_s_default, 0.677, 1e-9, "mu_s == 0.677");
      TEST_NEAR(cfg.slip_weakening->mu_d_default, 0.525, 1e-9, "mu_d == 0.525");
      TEST_NEAR(cfg.slip_weakening->d_c_default, 0.40, 1e-9, "d_c == 0.40");
      TEST_ASSERT(cfg.slip_weakening->spatial.size() == 3u,
                  "3 barrier rules (deep + 2 strike; free surface on top, no barrier)");
   }
}

int main(int, char**)
{
   std::cout << "Running Phase 8 TPV config-parse tests\n";

   // Numerics-selector helper tests run unconditionally (no files needed).
   T_numerics_dispatch_helpers();
   T_G2_MatrixMixedFluxUnderAder();   // Phase 5 Tests 5.1/5.2 (no files needed)

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

   // (Phase 5, Test 5.3) TPV31 (matrix + depth_profile_1d): regression + the
   // new matrix+adjacent+seam_continuous parse.
   const std::string p31    = FindConfig("tpv31/configs/tpv31.toml");
   const std::string p31_p2 = FindConfig("tpv31/configs/tpv31_p2.toml");
   const std::string p31_p3 = FindConfig("tpv31/configs/tpv31_p3.toml");
   const std::string p31_rkmf = FindConfig("tpv31/configs/tpv31_rk_mixedflux.toml");
   const std::string p31_p2_rkmf = FindConfig("tpv31/configs/tpv31_p2_rk_mixedflux.toml");
   if (!p31.empty())
   {
      T_TPV31_Regression(p31);          // 5.3c
      T_TPV31_MatrixAdjacentSeam(p31);  // 5.3a + 5.3b
      if (!p31_rkmf.empty()) { T_TPV31_RkMixedFlux(p31_rkmf); }  // Phase 7 artifact (p1)
      // Phase 7 p2 artifact: same matrix+adjacent+seam_continuous + order=2.
      if (!p31_p2_rkmf.empty())
      {
         std::cout << "\n[TPV31 p2 rk_mixedflux artifact] parse " << p31_p2_rkmf << "\n";
         spatial::SpatialFrictionConfig cfg =
            spatial::LoadSpatialFrictionConfig(p31_p2_rkmf);
         TEST_ASSERT(cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix,
                     "tpv31_p2_rk_mixedflux interior_flux == matrix");
         TEST_ASSERT(cfg.numerics.mixed_flux == "adjacent",
                     "tpv31_p2_rk_mixedflux mixed_flux == adjacent (parses past G1)");
         TEST_ASSERT(cfg.material.seam_continuous == true,
                     "tpv31_p2_rk_mixedflux [material].seam_continuous == true");
         TEST_ASSERT(cfg.mesh.order == 2,
                     "tpv31_p2_rk_mixedflux [mesh].order == 2 (p2)");
      }
      // p2/p3 share the parser; assert only the additive field's default
      // (their other fields are not part of the previously-asserted corpus).
      for (const std::string& pv : {p31_p2, p31_p3})
      {
         if (pv.empty()) { continue; }
         std::cout << "\n[TPV31 p-variant] parse " << pv << "\n";
         spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(pv);
         TEST_ASSERT(cfg.material.seam_continuous == false,
                     "tpv31 p-variant seam_continuous defaults false (key omitted)");
         TEST_ASSERT(cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix,
                     "tpv31 p-variant interior_flux == matrix (unchanged)");
      }
   }
   else
   {
      std::cout << "(tpv31.toml not found relative to CWD; TPV31 parse tests "
                   "skipped)\n";
   }

   // Part C / C2: TPV6 / TPV7 bi-material-fault config artifacts.
   const std::string p6 = FindConfig("tpv6/configs/tpv6.toml");
   const std::string p7 = FindConfig("tpv7/configs/tpv7.toml");
   if (!p6.empty()) { T_TPV6(p6, "TPV6 config",  3750.0); }  // far = slow (high contrast)
   else { std::cout << "(tpv6/configs/tpv6.toml not found; TPV6 parse skipped)\n"; }
   if (!p7.empty()) { T_TPV6(p7, "TPV7 config",  5000.0); }  // far = low-contrast
   else { std::cout << "(tpv7/configs/tpv7.toml not found; TPV7 parse skipped)\n"; }

   std::cout << "\n========================================\n";
   std::cout << "Phase 8 config-parse: " << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
