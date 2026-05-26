// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_constant_tensor_sign.cpp — Phase 3 (D3.1, R-007) golden sign test.
//
// Plan: PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24, §9.3 req 8.
//
// For each shipped SAFS constant_tensor config, build its stress block via
// the DRIVER's factory construction path (the D3.1 `-sigma_xy_pa`
// negation), project onto a canonical-frame fault DOF using the safs
// no-flip rule, and assert the on-fault components equal the pre-change
// GOLDEN literals captured from the current configs.  Because the test
// compares against FIXED golden numbers (not the config's own value), a
// half-applied D3.1 change — factory negation present but a config NOT
// flipped, or vice-versa — produces tau_strike of the wrong sign and FAILS
// loudly on that config (R-007 atomicity guard).
//
// Canonical y=0 vertical strike-slip fault frame (D3.1):
//   n  = ( 0, -1,  0)  (fault normal)
//   t1 = ( 0,  0, -1)  (dip)
//   t2 = (+1,  0,  0)  (strike)
// No-flip rule:  T = S·n ;  sigma_n = n·T ;  tau_dip = t1·T ;  tau_strike = t2·T.

#include "mfem.hpp"

#include "../../spatial/code/spatial_stress.hpp"     // ConstantTensorStressSource
#include "../../spatial/code/spatial_friction.hpp"   // LoadSpatialFrictionConfig

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

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
   << " (got " << _v << ", expected " << _e \
   << ", tol " << _t << ")\n"; num_failed++; } \
   else { std::cout << "  PASSED: " << m << "\n"; num_passed++; } \
   } while (0)

namespace
{
// Golden on-fault literals captured from the current (post-flip) SAFS
// configs — the values every flipped config + the driver negation must
// reproduce.  sigma_n = +sigma_yy, tau_dip = +sigma_yz, tau_strike =
// +|sigma_xy| (right-lateral positive).
constexpr real_t kGoldenSigmaN    = 1.1143380852144492e+08;
constexpr real_t kGoldenTauStrike = 9.888543819998317e+06;
constexpr real_t kGoldenTauDip    = 0.0;
constexpr real_t kAbsTol          = 1.0;   // 1 Pa (≪ the 2·|sigma_xy| flip delta)

const char* kConfigDir =
   "safs/project_7.0_alternative/config/";

struct NamedConfig { const char* label; std::string path; };

// Project the stored Cauchy tensor S onto the canonical fault frame.
void project_on_fault(const DenseMatrix& S,
                      real_t& sigma_n, real_t& tau_dip, real_t& tau_strike)
{
   const real_t n[3]  = { 0.0, -1.0,  0.0 };
   const real_t t1[3] = { 0.0,  0.0, -1.0 };   // dip
   const real_t t2[3] = { 1.0,  0.0,  0.0 };   // strike
   real_t T[3];
   for (int r = 0; r < 3; ++r)
   {
      T[r] = S(r, 0) * n[0] + S(r, 1) * n[1] + S(r, 2) * n[2];
   }
   sigma_n    = n[0]  * T[0] + n[1]  * T[1] + n[2]  * T[2];
   tau_dip    = t1[0] * T[0] + t1[1] * T[1] + t1[2] * T[2];
   tau_strike = t2[0] * T[0] + t2[1] * T[1] + t2[2] * T[2];
}

void check_config(const NamedConfig& nc)
{
   std::cout << "\n[D3.1/" << nc.label << "] " << nc.path << "\n";
   spatial::SpatialFrictionConfig cfg =
      spatial::LoadSpatialFrictionConfig(nc.path);
   TEST_ASSERT(cfg.stress.kind == spatial::StressSourceKind::ConstantTensor,
               "config uses constant_tensor stress");

   // Driver factory construction path: D3.1 negates sigma_xy_pa.
   spatial::ConstantTensorStressSource src(
      cfg.stress.sigma_xx_pa, cfg.stress.sigma_yy_pa, cfg.stress.sigma_zz_pa,
      -cfg.stress.sigma_xy_pa,
      cfg.stress.sigma_yz_pa, cfg.stress.sigma_xz_pa);

   real_t sigma_n, tau_dip, tau_strike;
   project_on_fault(src.Tensor(), sigma_n, tau_dip, tau_strike);

   TEST_NEAR(tau_strike, kGoldenTauStrike, kAbsTol,
             "tau_strike == +golden (right-lateral positive; flip applied)");
   TEST_NEAR(sigma_n, kGoldenSigmaN, kAbsTol,
             "sigma_n == +sigma_yy (compression positive)");
   TEST_NEAR(tau_dip, kGoldenTauDip, kAbsTol,
             "tau_dip == +sigma_yz");
   // The flipped config stores sigma_xy_pa right-lateral-POSITIVE.
   TEST_ASSERT(cfg.stress.sigma_xy_pa > 0.0,
               "config sigma_xy_pa is right-lateral-positive (flip applied)");
}
}  // namespace

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase 3 test_constant_tensor_sign\n";

   // Precondition: the canonical frame this test assumes.
   const real_t dip[3]    = { 0.0, 0.0, -1.0 };
   const real_t strike[3] = { 1.0, 0.0,  0.0 };
   TEST_ASSERT(dip[0] == 0.0 && dip[1] == 0.0 && dip[2] == -1.0,
               "precondition dip == (0,0,-1)");
   TEST_ASSERT(strike[0] == 1.0 && strike[1] == 0.0 && strike[2] == 0.0,
               "precondition strike == (+1,0,0)");

   const std::string dir = kConfigDir;
   std::vector<NamedConfig> configs = {
      { "base",
        dir + "spatial_friction_slip_weakening_safs_projected_stress.toml" },
      { "Dc2",
        dir + "spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml" },
      { "Dc8",
        dir + "spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc8.toml" },
      { "Dc10",
        dir + "spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc10.toml" },
      { "rate_state",
        dir + "spatial_friction_rate_state_safs_projected_stress_resolution_Dc2.toml" },
   };
   for (const auto& nc : configs) { check_config(nc); }

   std::cout << "\n========================================\n";
   std::cout << "Phase 3 test_constant_tensor_sign: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}
