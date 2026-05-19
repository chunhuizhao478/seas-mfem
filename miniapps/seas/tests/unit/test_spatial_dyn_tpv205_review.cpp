// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_dyn_tpv205_review.cpp — REVIEW R-001..R-008 unit-test
// coverage for the (`seas_spatial_dyn_driver`,
// `tpv205/configs/tpv205.toml`) parity contract.  See
// `miniapps/seas/debug_document/general_driver_debug_document/REVIEW.md`
// for the full set of findings each test guards against.
//
// Tests load the on-disk TPV205 TOML (relative path
// `tpv205/configs/tpv205.toml`) so any regression that re-flips the
// patch z-sign / barrier z-sign / ref_normal / hypocenter z, etc. fails
// this test BEFORE it reaches a production run.

#include "mfem.hpp"

#include "../../config/tpv205_params.hpp"
#include "../../spatial/code/spatial_friction.hpp"
#include "../../spatial/code/spatial_stress.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

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

// ---------------------------------------------------------------------
// Locate the TPV205 TOML on disk relative to a few candidate working
// directories (Makefile-default vs ctest CWD vs invoker-from-anywhere).
// ---------------------------------------------------------------------
static std::string find_tpv205_toml()
{
   const char* candidates[] = {
      "miniapps/seas/tpv205/configs/tpv205.toml",   // from MFEM root
      "tpv205/configs/tpv205.toml",                 // from miniapps/seas
      "../tpv205/configs/tpv205.toml",              // from miniapps/seas/tests/unit
      "../../tpv205/configs/tpv205.toml",
   };
   for (const char* p : candidates)
   {
      std::ifstream f(p);
      if (f.good()) { return std::string(p); }
   }
   return std::string();
}

// ---------------------------------------------------------------------
// R-001  TPV205 TOML stress patches use MESH-z coords; the central
//        nucleation patch lights up the hypocenter DOF at (0,0,-7500).
// ---------------------------------------------------------------------
static void R001_patches_match_mesh_z()
{
   std::cout << "\n[R-001] TPV205 patches placed at mesh z = -7500 "
                "(not spec depth +7500)\n";
   const std::string toml_path = find_tpv205_toml();
   TEST_ASSERT(!toml_path.empty(),
               "tpv205.toml found on disk");
   if (toml_path.empty()) { return; }

   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(toml_path);

   ConstantTensorWithPatchesStressSource src(
      cfg.stress.sigma_xx_pa, cfg.stress.sigma_yy_pa, cfg.stress.sigma_zz_pa,
      cfg.stress.sigma_xy_pa, cfg.stress.sigma_yz_pa, cfg.stress.sigma_xz_pa,
      cfg.stress.patches);

   // Hypocenter DOF (mesh frame: z negative below the free surface).
   const auto S_hypo = src.Evaluate(0.0, 0.0, -7500.0);
   TEST_NEAR(S_hypo(0, 1), TPV205Params::tau_nuc, 1.0,
             "central patch fires at mesh (0,0,-7500): sigma_xy = 81.6 MPa");

   // Left release patch (spec §9).
   const auto S_left = src.Evaluate(-7500.0, 0.0, -7500.0);
   TEST_NEAR(S_left(0, 1), TPV205Params::tau_left, 1.0,
             "left patch fires at mesh (-7500,0,-7500): sigma_xy = 78 MPa");

   // Right release patch (spec §8).
   const auto S_right = src.Evaluate(+7500.0, 0.0, -7500.0);
   TEST_NEAR(S_right(0, 1), TPV205Params::tau_right, 1.0,
             "right patch fires at mesh (+7500,0,-7500): sigma_xy = 62 MPa");

   // Background everywhere else.
   const auto S_bg = src.Evaluate(0.0, 0.0, -100.0);
   TEST_NEAR(S_bg(0, 1), TPV205Params::tau_back, 1.0,
             "background at shallow mesh z: sigma_xy = 70 MPa");

   // R-001 regression guard: the OLD spec-positive coord MUST miss.
   const auto S_old = src.Evaluate(0.0, 0.0, +7500.0);
   TEST_NEAR(S_old(0, 1), TPV205Params::tau_back, 1.0,
             "spec-positive z = +7500 misses patches (mesh-coord contract)");

   // Verify hypocenter strike pre-stress exceeds yield (rupture must
   // initiate at t = 0+).
   const real_t yield_at_hypo =
      TPV205Params::mu_s * TPV205Params::sigma_n;
   TEST_ASSERT(S_hypo(0, 1) > yield_at_hypo,
               "hypocenter sigma_xy > mu_s * sigma_n (rupture initiates)");
}

// ---------------------------------------------------------------------
// R-002  TPV205 bottom-barrier rule uses z_max_m = -15000 and fires for
//        deep mesh DOFs.
// ---------------------------------------------------------------------
static void R002_bottom_barrier_fires_deep()
{
   std::cout << "\n[R-002] TPV205 bottom-barrier locks DOFs at mesh z < -15 km\n";
   const std::string toml_path = find_tpv205_toml();
   if (toml_path.empty())
   {
      std::cerr << "  SKIP: tpv205.toml not found\n"; return;
   }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(toml_path);
   TEST_ASSERT(cfg.slip_weakening.has_value(),
               "tpv205.toml has a [friction.slip_weakening] block");
   if (!cfg.slip_weakening.has_value()) { return; }

   // Synthesize three DOFs:
   //  (a) inside the rupture area  → mu_s = mu_s_default = 0.677
   //  (b) deep (mesh z = -16 km)   → bottom barrier → mu_s = sentinel
   //  (c) along-strike outside (|x| > 15 km) → side barrier → sentinel
   Vector coords(3 * 3);
   coords(0) = 0.0;       coords(1) = 0.0;  coords(2) = -5000.0;   // (a)
   coords(3) = 0.0;       coords(4) = 0.0;  coords(5) = -16000.0;  // (b)
   coords(6) = -16000.0;  coords(7) = 0.0;  coords(8) = -5000.0;   // (c)
   Array<int> attr(3);
   attr[0] = cfg.boundary.fault_attr;
   attr[1] = cfg.boundary.fault_attr;
   attr[2] = cfg.boundary.fault_attr;

   SpatialFrictionResolver resolver;
   const SlipWeakeningPerDOFParams p =
      resolver.ResolveSlipWeakening(*cfg.slip_weakening, coords, attr);

   TEST_NEAR(p.mu_s(0), 0.677, 1e-9,
             "(a) inside rupture area: mu_s = default 0.677");
   TEST_ASSERT(p.mu_s(1) > 5000.0,
               "(b) deep DOF z = -16 km: bottom-barrier sentinel mu_s > 5000");
   TEST_ASSERT(p.mu_s(2) > 5000.0,
               "(c) along-strike DOF |x|>15 km: side-barrier sentinel mu_s > 5000");

   // Regression guard: the OLD `z_min_m = 15000` would not have matched
   // the deep DOF at mesh z = -16000; explicitly confirm the fix flipped
   // the sign.  (Encoded by the (b) assertion above; this comment makes
   // the intent grep-able.)
}

// ---------------------------------------------------------------------
// R-003  TPV205 TOML ref_normal / up inherit the canonical
//        (0,-1,0)/(0,0,1) — NOT the bogus (0,+1,0) SAFS default.
// ---------------------------------------------------------------------
static void R003_canonical_ref_normal()
{
   std::cout << "\n[R-003] TPV205 ref_normal matches the canonical "
                "Tandem convention\n";
   const std::string toml_path = find_tpv205_toml();
   if (toml_path.empty())
   {
      std::cerr << "  SKIP: tpv205.toml not found\n"; return;
   }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(toml_path);
   TEST_NEAR(cfg.fault_geometry.ref_normal[0], 0.0, 1e-12,
             "ref_normal[0] = 0");
   TEST_NEAR(cfg.fault_geometry.ref_normal[1], -1.0, 1e-12,
             "ref_normal[1] = -1 (Tandem convention; matches wave op)");
   TEST_NEAR(cfg.fault_geometry.ref_normal[2], 0.0, 1e-12,
             "ref_normal[2] = 0");
   TEST_NEAR(cfg.fault_geometry.up[0], 0.0, 1e-12, "up[0] = 0");
   TEST_NEAR(cfg.fault_geometry.up[1], 0.0, 1e-12, "up[1] = 0");
   TEST_NEAR(cfg.fault_geometry.up[2], 1.0, 1e-12, "up[2] = 1");
}

// ---------------------------------------------------------------------
// R-004/R-005/R-006  TPV205 TOML opts in to native-driver dispatch.
//
// Updated 2026-05-19 per review of branch feature/heterogeneous_riemann_solver
// (debug_document/general_driver_debug_document/
//  spatial_dyn_heterogeneous_riemann_review_2026-05-19.md):
// R-006 used to assert `interior_flux = "scalar"` (byte-parity opt-out),
// but the new branch review's R-001 flagged that as defeating the
// branch purpose — the canonical TOML must default to "bimaterial"
// (heterogeneous Riemann) so the new WaveOperator ctor is exercised.
// Updated below.
// ---------------------------------------------------------------------
static void R004_R005_R006_dispatch_opt_ins()
{
   std::cout << "\n[R-004/005/006] TPV205 numerics opt in to native dispatch\n";
   const std::string toml_path = find_tpv205_toml();
   if (toml_path.empty())
   {
      std::cerr << "  SKIP: tpv205.toml not found\n"; return;
   }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(toml_path);
   TEST_ASSERT(cfg.numerics.cfl_safety == "dg",
               "R-004: cfl_safety = dg (matches native cfl/9 scaling)");
   TEST_ASSERT(cfg.numerics.fault_iterator == "one-shot",
               "R-005: fault_iterator = one-shot (native TPV205 default)");
   TEST_ASSERT(cfg.numerics.interior_flux == "bimaterial",
               "R-006 (post-2026-05-19): interior_flux defaults to "
               "\"bimaterial\" — heterogeneous Riemann path is the "
               "branch's intended dispatch; \"scalar\" is reserved for "
               "byte-parity regressions only.");
}

// ---------------------------------------------------------------------
// R-007  [numerics].mixed_flux remains "none" (REQUIRED for bi-material
//         runs; the parser enforces the valid set).
// ---------------------------------------------------------------------
static void R007_mixed_flux_none()
{
   std::cout << "\n[R-007] TPV205 mixed_flux = none\n";
   const std::string toml_path = find_tpv205_toml();
   if (toml_path.empty())
   {
      std::cerr << "  SKIP: tpv205.toml not found\n"; return;
   }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(toml_path);
   TEST_ASSERT(cfg.numerics.mixed_flux == "none",
               "TPV205 mixed_flux = none");
}

// ---------------------------------------------------------------------
// R-008  Hypocenter z lives in mesh frame (z = -7500, not +7500).  The
//         parser aborts if a positive z is set with up_z > 0.
// ---------------------------------------------------------------------
static void R008_hypocenter_mesh_z()
{
   std::cout << "\n[R-008] TPV205 hypocenter z = -7500 (mesh frame)\n";
   const std::string toml_path = find_tpv205_toml();
   if (toml_path.empty())
   {
      std::cerr << "  SKIP: tpv205.toml not found\n"; return;
   }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(toml_path);
   TEST_NEAR(cfg.hypocenter.x,  0.0,    1e-12, "hypo x = 0");
   TEST_NEAR(cfg.hypocenter.y,  0.0,    1e-12, "hypo y = 0");
   TEST_NEAR(cfg.hypocenter.z, -7500.0, 1e-9,
             "hypo z = -7500 (mesh frame; was +7500 pre-R-008)");
}

// NOTE: a "parser aborts on positive hypocenter z" death test is NOT
// included here because MFEM is built without MFEM_USE_EXCEPTIONS, so
// MFEM_VERIFY calls abort() rather than throwing — try/catch can't
// observe the abort.  The R-008 guard's behaviour is instead asserted
// in two ways:
//   1. The OK case (R008_hypocenter_mesh_z above) loads the real
//      tpv205.toml and confirms z = -7500 parses without aborting.
//   2. A manual smoke test: editing tpv205.toml to z = +7500 and
//      running `seas_test_spatial_dyn_tpv205_review` should terminate
//      with the R-008 abort message.  See
//      `debug_document/general_driver_debug_document/REVIEW.md`
//      R-008 for the diagnostic the abort emits.

int main(int, char**)
{
   std::cout << "Running test_spatial_dyn_tpv205_review "
                "(REVIEW.md R-001..R-008 coverage)\n";
   R001_patches_match_mesh_z();
   R002_bottom_barrier_fires_deep();
   R003_canonical_ref_normal();
   R004_R005_R006_dispatch_opt_ins();
   R007_mixed_flux_none();
   R008_hypocenter_mesh_z();
   std::cout << "\n========================================\n";
   std::cout << "test_spatial_dyn_tpv205_review: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
