// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_review_2026_05_19_fixes.cpp — regression tests for the code
// review findings dated 2026-05-19
// (debug_document/general_driver_debug_document/
//  spatial_dyn_heterogeneous_riemann_review_2026-05-19.md).
//
// Each test covers one finding from the review and would have failed
// before the corresponding fix was applied.  The tests deliberately
// avoid spinning up MPI or loading meshes — they exercise the parser /
// resolver / DepthProfile1DMaterial::eval_at_xyz directly.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"
#include "../../dynamic/heterogeneous_material.hpp"

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

static std::string find_toml(const char *relpath_under_seas)
{
   const std::string base[] = {
      "miniapps/seas/",
      "",
      "../",
      "../../",
   };
   for (const auto &b : base)
   {
      const std::string full = b + std::string(relpath_under_seas);
      std::ifstream f(full);
      if (f.good()) { return full; }
   }
   return std::string();
}

// ---------------------------------------------------------------------
// R-001:  tpv205.toml and tpv31.toml both default `interior_flux` to
//         "bimaterial" (heterogeneous Riemann solver) — neither TOML
//         re-introduces the `interior_flux = "scalar"` opt-out that
//         defeats the branch's purpose.
// ---------------------------------------------------------------------
static void T_R001_interior_flux_defaults_bimaterial()
{
   std::cout << "\n[T-R001] tpv205/tpv31 TOMLs default to bimaterial "
             << "interior flux (heterogeneous Riemann path)\n";
   for (const char *rel : { "tpv205/configs/tpv205.toml",
                             "tpv31/configs/tpv31.toml" })
   {
      const std::string path = find_toml(rel);
      TEST_ASSERT(!path.empty(), std::string(rel) + " on disk");
      if (path.empty()) { continue; }
      SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
      TEST_ASSERT(cfg.numerics.interior_flux == "bimaterial",
                  std::string(rel) + " uses interior_flux=bimaterial");
   }
}

// ---------------------------------------------------------------------
// R-004:  all four canonical TPV TOMLs flip the paraview_enabled
//         master gate to true so the per-collection mode strings
//         actually drive output.
// ---------------------------------------------------------------------
static void T_R004_paraview_enabled_set_in_all_tomls()
{
   std::cout << "\n[T-R004] All four canonical TOMLs set paraview_enabled = true\n";
   for (const char *rel : { "tpv102/configs/tpv102.toml",
                             "tpv104/configs/tpv104.toml",
                             "tpv205/configs/tpv205.toml",
                             "tpv31/configs/tpv31.toml" })
   {
      const std::string path = find_toml(rel);
      TEST_ASSERT(!path.empty(), std::string(rel) + " on disk");
      if (path.empty()) { continue; }
      SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
      TEST_ASSERT(cfg.output.paraview_enabled,
                  std::string(rel) + " sets paraview_enabled=true");
   }
}

// ---------------------------------------------------------------------
// R-003 / R-008 support:  DepthProfile1DMaterial exposes a populated
// `eval_at_xyz` callback so the depth-proportional stress source can
// look up μ at any (x, y, z) without an ElementTransformation.
// ---------------------------------------------------------------------
static void T_R003_depth_profile_eval_at_xyz_populated()
{
   std::cout << "\n[T-R003] DepthProfile1DMaterial::eval_at_xyz returns "
             << "spec-correct (λ, μ, ρ) along axis='z'\n";
   std::vector<DepthProfileLayer> layers;
   // Spec p. 4 row 1 (constant 0–2400 m).
   layers.push_back({0.0,    2400.0,  4050.0, 2250.0, 2580.0, "constant"});
   // Spec p. 4 row 2 (linear 2400–5000 m to (5200, 3050, 2620)).
   layers.push_back({2400.0, 5000.0,  4450.0, 2550.0, 2600.0, "linear"});
   // Linear endpoint (5000− values).
   layers.push_back({5000.0, 5000.0,  5200.0, 3050.0, 2620.0, "constant"});
   // 5000–10 000 m constant.
   layers.push_back({5000.0, 10000.0, 5750.0, 3450.0, 2720.0, "constant"});
   // ≥10 km constant.
   layers.push_back({10000.0, 30000.0, 6500.0, 3800.0, 3000.0, "constant"});

   auto wrapper = MakeDepthProfile1DMaterial(layers, /*axis=*/'z');
   TEST_ASSERT(wrapper != nullptr, "MakeDepthProfile1DMaterial returns wrapper");
   if (!wrapper) { return; }
   TEST_ASSERT(static_cast<bool>(wrapper->eval_at_xyz),
               "wrapper->eval_at_xyz is populated");
   if (!wrapper->eval_at_xyz) { return; }

   // Probe at surface (z = 0, depth = 0): layer 0 values.
   real_t lam, mu, rho;
   wrapper->eval_at_xyz(0.0, 0.0, 0.0, lam, mu, rho);
   const real_t expected_mu_top = 2580.0 * 2250.0 * 2250.0;     // ρ·c_s²
   TEST_NEAR(mu, expected_mu_top, 1.0e-3 * expected_mu_top,
             "μ at surface = ρ·c_s² (layer 0)");
   TEST_NEAR(rho, 2580.0, 1e-9, "ρ at surface = 2580 kg/m³");

   // Probe at depth 7500 m (z = -7500, hypocenter depth): layer 3
   // values (constant 5000-10000 m).
   wrapper->eval_at_xyz(0.0, 0.0, -7500.0, lam, mu, rho);
   const real_t expected_mu_hypo = 2720.0 * 3450.0 * 3450.0;
   TEST_NEAR(mu, expected_mu_hypo, 1.0e-3 * expected_mu_hypo,
             "μ at depth 7.5 km (layer 3)");
   TEST_NEAR(rho, 2720.0, 1e-9, "ρ at depth 7.5 km = 2720");

   // Probe above the free surface (z = +100): clamps to layer 0 (depth=0).
   wrapper->eval_at_xyz(0.0, 0.0, 100.0, lam, mu, rho);
   TEST_NEAR(mu, expected_mu_top, 1.0e-3 * expected_mu_top,
             "μ above surface clamps to layer 0");
}

// ---------------------------------------------------------------------
// R-009:  tpv205.toml barrier rules use epsilon-nudged bounds so the
//         resolver does NOT mark DOFs at exactly x=±15000 / z=−15000
//         as barrier.  We probe ResolveSlipWeakening directly.
// ---------------------------------------------------------------------
static void T_R009_barrier_boundary_inside_rupture_area()
{
   std::cout << "\n[T-R009] DOFs at x = ±15 km / z = -15 km stay inside "
             << "the TPV205 rupture area (μ_s = 0.677, NOT the 1.0e6 barrier)\n";
   const std::string path = find_toml("tpv205/configs/tpv205.toml");
   TEST_ASSERT(!path.empty(), "tpv205.toml on disk");
   if (path.empty()) { return; }
   SpatialFrictionConfig cfg = LoadSpatialFrictionConfig(path);
   TEST_ASSERT(cfg.slip_weakening.has_value(),
               "has [friction.slip_weakening]");
   if (!cfg.slip_weakening.has_value()) { return; }

   // Probe DOFs:
   //   (x, y, z) = (-15000, 0, -7500)  expect mu_s = 0.677  (left boundary)
   //   (x, y, z) = (+15000, 0, -7500)  expect mu_s = 0.677  (right boundary)
   //   (x, y, z) = (0,      0, -15000) expect mu_s = 0.677  (bottom boundary)
   //   (x, y, z) = (-15001, 0, -7500)  expect mu_s = 1.0e6  (just outside)
   const int N = 4;
   Vector coords(3 * N);
   coords(0) = -15000.0; coords(1) = 0.0; coords(2)  = -7500.0;
   coords(3) =  15000.0; coords(4) = 0.0; coords(5)  = -7500.0;
   coords(6) =      0.0; coords(7) = 0.0; coords(8)  = -15000.0;
   coords(9) = -15001.0; coords(10)= 0.0; coords(11) = -7500.0;
   Array<int> attrs(N);
   for (int i = 0; i < N; ++i) { attrs[i] = cfg.boundary.fault_attr; }

   SpatialFrictionResolver resolver;
   auto p = resolver.ResolveSlipWeakening(*cfg.slip_weakening,
                                          coords, attrs);

   TEST_NEAR(p.mu_s(0), 0.677, 1e-9, "DOF at x = -15000 is INSIDE rupture area");
   TEST_NEAR(p.mu_s(1), 0.677, 1e-9, "DOF at x = +15000 is INSIDE rupture area");
   TEST_NEAR(p.mu_s(2), 0.677, 1e-9, "DOF at z = -15000 is INSIDE rupture area");
   TEST_ASSERT(p.mu_s(3) > 1.0e5,    "DOF at x = -15001 is barrier (μ_s > 1e5)");
}

int main()
{
   std::cout << "test_review_2026_05_19_fixes — regression for the "
             << "spatial_dyn_driver heterogeneous-Riemann review.\n";

   T_R001_interior_flux_defaults_bimaterial();
   T_R004_paraview_enabled_set_in_all_tomls();
   T_R003_depth_profile_eval_at_xyz_populated();
   T_R009_barrier_boundary_inside_rupture_area();

   std::cout << "\n========================================\n"
             << "Tests:  total=" << num_tests
             << " passed=" << num_passed
             << " failed=" << num_failed << "\n";
   return (num_failed == 0) ? 0 : 1;
}
