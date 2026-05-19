// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_friction_config_phaser.cpp — Phase R.3 of
// PLAN_phase_R_exact_bimaterial_riemann_rev3.md.
//
// R.3.T-1 + R.3.T-2 acceptance tests for the new TOML sections added
// in Phase R.3 step 1:
//   * [problem]
//   * [boundary]
//   * [fault_geometry]
//   * [hypocenter]
//   * [material] + [[material_profile.layer]]
//   * [stress] kind = "depth_proportional" + [stress.depth_proportional]
//
// Without SEAS_USE_TOML these tests are stubs (the parser aborts
// at runtime with a "rebuild with toml11" message; the test framework
// detects the missing dependency and prints SKIPPED).

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace mfem;
using namespace mfem::seas::spatial;

static int g_num_tests = 0, g_num_passed = 0, g_num_failed = 0;

#define TEST_ASSERT(c, m) do { g_num_tests++; if (!(c)) {                  \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n";           \
   g_num_failed++; } else { std::cout << "  PASSED: " << m << "\n";        \
   g_num_passed++; } } while (0)

namespace
{
std::string CommonHeader()
{
   return
      "[meta]\n"
      "schema_version = 1\n"
      "law = \"slip_weakening\"\n"
      "[material_constant_fallback]\n"
      "lambda = 32e9\nmu = 32e9\nrho = 2670.0\n"
      "[pore_pressure]\n"
      "P_p_pa = 0.0\nP_p_grad_pa_per_m = 0.0\nmin_sigma_n_pa = 0.0\n"
      "[mesh]\npath = \"x.msh\"\norder = 1\n"
      "[numerics]\nader_order = 2\nmixed_flux = \"none\"\ncfl = 0.5\n"
      "use_pml = false\n"
      "[time]\ntfinal = \"12s\"\nt_initial = \"0s\"\n"
      "dt_initial = \"auto\"\ndt_max = \"0.01s\"\n"
      "[output]\noutput_dir = \"out\"\n"
      "restart_prefix = \"cp\"\n"
      "max_snapshots = 100\ncheckpoint_every_steps = 100\n"
      "[friction.slip_weakening]\n"
      "mu_s_default = 0.677\nmu_d_default = 0.525\n"
      "d_c_default = 0.4\ncohesion_default = 0.0\n";
}
}  // namespace

#ifdef SEAS_USE_TOML
static void T_problem_defaults()
{
   std::cout << "\n[T-problem] default-tag round-trip\n";
   std::string toml = CommonHeader() +
      "[stress]\nkind = \"constant_tensor\"\n"
      "sigma_xx_pa = 0\nsigma_yy_pa = 0\nsigma_zz_pa = 0\n"
      "sigma_xy_pa = 0\nsigma_yz_pa = 0\nsigma_xz_pa = 0\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.problem.tag == "safs",
               "absent [problem] block defaults tag = 'safs'");
   TEST_ASSERT(cfg.boundary.fault_attr == 101,
               "absent [boundary] block defaults fault_attr = 101");
   TEST_ASSERT(cfg.fault_geometry.kind == "bp5_safs",
               "absent [fault_geometry] defaults kind = 'bp5_safs'");
   TEST_ASSERT(cfg.material.kind == MaterialKind::Constant,
               "absent [material] defaults kind = Constant");
}

static void T_problem_tag()
{
   std::cout << "\n[T-tag] [problem].tag round-trip\n";
   std::string toml = CommonHeader() +
      "[stress]\nkind = \"constant_tensor\"\n"
      "sigma_xx_pa = 0\nsigma_yy_pa = 0\nsigma_zz_pa = 0\n"
      "sigma_xy_pa = 0\nsigma_yz_pa = 0\nsigma_xz_pa = 0\n"
      "[problem]\ntag = \"tpv205\"\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.problem.tag == "tpv205",
               "tag = 'tpv205' round-trips");
}

static void T_boundary_round_trip()
{
   std::cout << "\n[T-boundary] explicit ids round-trip\n";
   std::string toml = CommonHeader() +
      "[stress]\nkind = \"constant_tensor\"\n"
      "sigma_xx_pa = 0\nsigma_yy_pa = 0\nsigma_zz_pa = 0\n"
      "sigma_xy_pa = 0\nsigma_yz_pa = 0\nsigma_xz_pa = 0\n"
      "[boundary]\n"
      "fault_attr = 103\n"
      "natural_attrs = [101]\n"
      "absorbing_attrs = [105]\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.boundary.fault_attr == 103,
               "fault_attr = 103 round-trips");
   TEST_ASSERT(cfg.boundary.natural_attrs.size() == 1
               && cfg.boundary.natural_attrs[0] == 101,
               "natural_attrs = [101] round-trips");
   TEST_ASSERT(cfg.boundary.absorbing_attrs.size() == 1
               && cfg.boundary.absorbing_attrs[0] == 105,
               "absorbing_attrs = [105] round-trips");
}

static void T_fault_geometry_round_trip()
{
   std::cout << "\n[T-fault_geometry] orientation round-trip\n";
   std::string toml = CommonHeader() +
      "[stress]\nkind = \"constant_tensor\"\n"
      "sigma_xx_pa = 0\nsigma_yy_pa = 0\nsigma_zz_pa = 0\n"
      "sigma_xy_pa = 0\nsigma_yz_pa = 0\nsigma_xz_pa = 0\n"
      "[fault_geometry]\n"
      "ref_normal = [0.0, 0.0, 1.0]\n"
      "up         = [0.0, 1.0, 0.0]\n"
      "kind       = \"tpv31_lsw\"\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(std::abs(cfg.fault_geometry.ref_normal[2] - 1.0) < 1e-15,
               "ref_normal = (0,0,1) round-trips");
   TEST_ASSERT(std::abs(cfg.fault_geometry.up[1] - 1.0) < 1e-15,
               "up = (0,1,0) round-trips");
   TEST_ASSERT(cfg.fault_geometry.kind == "tpv31_lsw",
               "kind = 'tpv31_lsw' round-trips");
}

static void T_hypocenter_round_trip()
{
   std::cout << "\n[T-hypocenter] round-trip\n";
   std::string toml = CommonHeader() +
      "[stress]\nkind = \"constant_tensor\"\n"
      "sigma_xx_pa = 0\nsigma_yy_pa = 0\nsigma_zz_pa = 0\n"
      "sigma_xy_pa = 0\nsigma_yz_pa = 0\nsigma_xz_pa = 0\n"
      "[hypocenter]\n"
      "x = 0.0\ny = 7500.0\nz = 0.0\n"
      "nucleation_radius_m = 1400.0\n"
      "nucleation_taper_m  = 600.0\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(std::abs(cfg.hypocenter.y - 7500.0) < 1e-9,
               "hypocenter.y = 7500 round-trips");
   TEST_ASSERT(std::abs(cfg.hypocenter.nucleation_radius_m - 1400.0) < 1e-9,
               "nucleation_radius_m = 1400 round-trips");
   TEST_ASSERT(std::abs(cfg.hypocenter.nucleation_taper_m - 600.0) < 1e-9,
               "nucleation_taper_m = 600 round-trips");
}

static void T_material_depth_profile()
{
   std::cout << "\n[T-material] depth_profile_1d round-trip\n";
   std::string toml = CommonHeader() +
      "[stress]\nkind = \"constant_tensor\"\n"
      "sigma_xx_pa = 0\nsigma_yy_pa = 0\nsigma_zz_pa = 0\n"
      "sigma_xy_pa = 0\nsigma_yz_pa = 0\nsigma_xz_pa = 0\n"
      "[material]\nkind = \"depth_profile_1d\"\ndepth_axis = \"y\"\n"
      "[[material_profile.layer]]\n"
      "depth_top_m = 0.0\ndepth_bot_m = 2400.0\n"
      "vp_ms = 4050.0\nvs_ms = 2250.0\nrho_kgm3 = 2580.0\n"
      "interp = \"constant\"\n"
      "[[material_profile.layer]]\n"
      "depth_top_m = 2400.0\ndepth_bot_m = 5000.0\n"
      "vp_ms = 5500.0\nvs_ms = 3000.0\nrho_kgm3 = 2670.0\n"
      "interp = \"constant\"\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.material.kind == MaterialKind::DepthProfile1D,
               "material.kind = 'depth_profile_1d' round-trips");
   TEST_ASSERT(cfg.material.depth_axis == 'y',
               "depth_axis = 'y' round-trips");
   TEST_ASSERT(cfg.material.profile_layers.size() == 2,
               "2 profile layers round-trip");
   TEST_ASSERT(std::abs(cfg.material.profile_layers[0].vp_ms - 4050.0) < 1e-9,
               "layer[0].vp_ms = 4050 round-trips");
   TEST_ASSERT(std::abs(cfg.material.profile_layers[1].rho_kgm3 - 2670.0) < 1e-9,
               "layer[1].rho_kgm3 = 2670 round-trips");
}

static void T_stress_depth_proportional()
{
   std::cout << "\n[T-stress.depth_proportional] round-trip\n";
   std::string toml = CommonHeader() +
      "[stress]\nkind = \"depth_proportional\"\n"
      "[stress.depth_proportional]\n"
      "sigma_xx_per_mu = -60.0\nsigma_yy_per_mu = 0.0\n"
      "sigma_zz_per_mu = -60.0\nsigma_xy_per_mu = 0.0\n"
      "sigma_yz_per_mu = 0.0\nsigma_xz_per_mu = 30.0\n"
      "mu_ref_pa = 32.03812032e9\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.stress.kind == StressSourceKind::DepthProportionalToShearModulus,
               "stress.kind = 'depth_proportional' round-trips");
   // MPa → Pa conversion is done inside the parser.
   TEST_ASSERT(std::abs(cfg.stress.depth_proportional.sigma_xx_per_mu
                        - (-60.0e6)) < 1e-3,
               "sigma_xx_per_mu = -60 MPa converted to -60e6 Pa");
   TEST_ASSERT(std::abs(cfg.stress.depth_proportional.sigma_xz_per_mu
                        - 30.0e6) < 1e-3,
               "sigma_xz_per_mu = +30 MPa converted to +30e6 Pa");
   TEST_ASSERT(std::abs(cfg.stress.depth_proportional.mu_ref_pa
                        - 32.03812032e9) < 1.0,
               "mu_ref_pa round-trips");
}
#endif

int main()
{
   std::cout << "========================================\n";
   std::cout << "test_spatial_friction_config_phaser (R.3.T-1/T-2)\n";
   std::cout << "========================================\n";

#ifndef SEAS_USE_TOML
   std::cout << "SEAS_USE_TOML not defined — skipping (toml11 not in this build).\n";
   return 0;
#else
   T_problem_defaults();
   T_problem_tag();
   T_boundary_round_trip();
   T_fault_geometry_round_trip();
   T_hypocenter_round_trip();
   T_material_depth_profile();
   T_stress_depth_proportional();

   std::cout << "\n========================================\n";
   std::cout << "Results: " << g_num_passed << "/" << g_num_tests
             << " passed, " << g_num_failed << " failed\n";
   std::cout << "========================================\n";
   return g_num_failed ? 1 : 0;
#endif
}
