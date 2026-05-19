// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_stress_with_patches.cpp — runtime-side coverage for
// `StressSourceKind::ConstantTensorWithPatches` + the matching
// `ConstantTensorWithPatchesStressSource` introduced to give the SAFS
// spatial driver a spec-exact pre-stress encoding for SCEC TPV205.
//
// Also covers `SquareOverstressSpec` / `ResolveSquareOverstress` /
// `ApplySquareOverstressIncrement` — the sibling rectangular-shape
// nucleation kind (smoothStep ramp + indicator function).

#include "mfem.hpp"

#include "../../config/tpv205_params.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/spatial_nucleation.hpp"
#include "../../spatial/code/spatial_friction.hpp"
#include "../../spatial/code/spatial_stress.hpp"

#include <cmath>
#include <iostream>
#include <limits>

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
// P-1  ConstantTensorWithPatchesStressSource degenerates to the
//      background tensor when the patch list is empty.
// ---------------------------------------------------------------------
static void P_1_empty_patches_degenerate()
{
   std::cout << "\n[P-1] Empty patches reproduce ConstantTensorStressSource\n";
   const real_t sxx = 0.0, syy = 120e6, szz = 0.0;
   const real_t sxy = 70e6, syz = 0.0,   sxz = 0.0;
   ConstantTensorWithPatchesStressSource src(sxx, syy, szz, sxy, syz, sxz,
                                             /*patches=*/{});
   const auto S = src.Evaluate(0, 0, 0);
   TEST_NEAR(S(0, 0), sxx, 0.0, "sxx");
   TEST_NEAR(S(1, 1), syy, 0.0, "syy");
   TEST_NEAR(S(2, 2), szz, 0.0, "szz");
   TEST_NEAR(S(0, 1), sxy, 0.0, "sxy");
   TEST_NEAR(S(1, 0), sxy, 0.0, "syx (symmetric)");
   TEST_NEAR(S(1, 2), syz, 0.0, "syz");
   TEST_NEAR(S(0, 2), sxz, 0.0, "sxz");
   const auto S_far = src.Evaluate(1e6, -1e6, 1e6);
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
         TEST_NEAR(S_far(i, j), S(i, j), 0.0, "far-field == origin");
      }
}

// ---------------------------------------------------------------------
// P-2  Inside a TPV205-style nucleation patch (centred at the
//      hypocenter, 3 km square), the along-strike shear (sigma_xy)
//      is overridden but the normal stress (sigma_yy) is inherited.
// ---------------------------------------------------------------------
static void P_2_tpv205_central_patch()
{
   std::cout << "\n[P-2] TPV205 central nucleation patch override\n";
   StressPatch nuc;
   nuc.center_x_m = 0.0;
   nuc.center_y_m = 0.0;
   nuc.center_z_m = 7500.0;
   nuc.half_x_m   = 1500.0;
   nuc.half_z_m   = 1500.0;
   // half_y_m left at +inf — y-bound unconstrained (TPV205 patches are
   // 2-D rectangles on the y=0 fault plane).
   nuc.sigma_xy_pa = TPV205Params::tau_nuc;   // 81.6 MPa
   // All other sigma_*_pa stay NaN → inherit background.

   ConstantTensorWithPatchesStressSource src(0.0, 120e6, 0.0,
                                              70e6, 0.0, 0.0,
                                              {nuc});

   // (a) DOF at the patch centre → sigma_xy = 81.6e6, sigma_yy = 120e6.
   {
      const auto S = src.Evaluate(0.0, 0.0, 7500.0);
      TEST_NEAR(S(0, 1), TPV205Params::tau_nuc, 0.0,
                "centre: sigma_xy overridden");
      TEST_NEAR(S(1, 1), 120e6, 0.0,
                "centre: sigma_yy inherited from background");
   }
   // (b) DOF at the patch edge (|dx| = half) → still inside.
   {
      const auto S = src.Evaluate(1500.0, 0.0, 7500.0);
      TEST_NEAR(S(0, 1), TPV205Params::tau_nuc, 0.0,
                "edge x=+half: still inside, override applied");
   }
   // (c) DOF just outside (|dx| > half) → background only.
   {
      const auto S = src.Evaluate(1500.0 + 1e-3, 0.0, 7500.0);
      TEST_NEAR(S(0, 1), 70e6, 0.0,
                "just outside x=+half: background sigma_xy");
   }
   // (d) DOF at z outside the patch but x inside → background.
   {
      const auto S = src.Evaluate(0.0, 0.0, 7500.0 + 1500.0 + 1.0);
      TEST_NEAR(S(0, 1), 70e6, 0.0,
                "outside z range: background sigma_xy");
   }
}

// ---------------------------------------------------------------------
// P-3  All three TPV205 patches active; the 62/78/81.6 MPa amplitudes
//      land at their respective centres; the background lands at the
//      far-field DOF.
// ---------------------------------------------------------------------
static void P_3_tpv205_three_patches_match_spec()
{
   std::cout << "\n[P-3] TPV205 3-patch encoding matches "
                "TPV205Params::ComputeTau2_0\n";
   struct Entry { real_t cx; real_t amp; };
   const Entry entries[] = {
      { 0.0,     TPV205Params::tau_nuc   },
      { -7500.0, TPV205Params::tau_left  },
      {  7500.0, TPV205Params::tau_right },
   };
   std::vector<StressPatch> patches;
   for (const Entry& e : entries)
   {
      StressPatch p;
      p.center_x_m  = e.cx;
      p.center_y_m  = 0.0;
      p.center_z_m  = 7500.0;
      p.half_x_m    = TPV205Params::patch_half;
      p.half_z_m    = TPV205Params::patch_half;
      p.sigma_xy_pa = e.amp;
      patches.push_back(p);
   }
   ConstantTensorWithPatchesStressSource src(0.0, 120e6, 0.0,
                                              TPV205Params::tau_back, 0.0, 0.0,
                                              patches);
   struct Sample { real_t x, z; real_t expect_sxy; const char* name; };
   const Sample samples[] = {
      { 0.0,     7500.0, TPV205Params::tau_nuc,   "nucleation centre" },
      { -7500.0, 7500.0, TPV205Params::tau_left,  "left centre"       },
      {  7500.0, 7500.0, TPV205Params::tau_right, "right centre"      },
      { 0.0,     0.0,    TPV205Params::tau_back,  "free-surface bg"   },
      { 12000.0, 7500.0, TPV205Params::tau_back,  "beyond patches"    },
      { -7500.0, 5800.0, TPV205Params::tau_back,  "near left, off z"  },
   };
   for (const auto& s : samples)
   {
      const auto S = src.Evaluate(s.x, 0.0, s.z);
      TEST_NEAR(S(0, 1), s.expect_sxy, 0.0, std::string("sigma_xy at ") + s.name);
   }
}

// ---------------------------------------------------------------------
// P-4  Last-match-wins on overlapping patches with the same component
//      override.
// ---------------------------------------------------------------------
static void P_4_last_match_wins()
{
   std::cout << "\n[P-4] Overlapping patches: last in list wins per component\n";
   StressPatch a;
   a.center_x_m = 0.0;  a.center_y_m = 0.0;  a.center_z_m = 0.0;
   a.half_x_m = 100.0;  a.half_z_m = 100.0;
   a.sigma_xy_pa = 10e6;
   StressPatch b = a;
   b.sigma_xy_pa = 20e6;
   ConstantTensorWithPatchesStressSource src(0,0,0, 0,0,0, {a, b});
   const auto S = src.Evaluate(0, 0, 0);
   TEST_NEAR(S(0, 1), 20e6, 0.0, "second patch's sigma_xy wins");
}

// ---------------------------------------------------------------------
// P-5  Symmetry: writing sigma_xy override also lands at S(1, 0).
// ---------------------------------------------------------------------
static void P_5_patch_override_is_symmetric()
{
   std::cout << "\n[P-5] Override preserves Cauchy symmetry\n";
   StressPatch p;
   p.center_x_m = 0; p.center_y_m = 0; p.center_z_m = 0;
   p.half_x_m = 1; p.half_z_m = 1;
   p.sigma_xy_pa = 42e6;
   p.sigma_yz_pa = 17e6;
   p.sigma_xz_pa = 5e6;
   ConstantTensorWithPatchesStressSource src(0,0,0, 0,0,0, {p});
   const auto S = src.Evaluate(0, 0, 0);
   TEST_NEAR(S(0, 1), S(1, 0), 0.0, "sigma_xy == sigma_yx");
   TEST_NEAR(S(1, 2), S(2, 1), 0.0, "sigma_yz == sigma_zy");
   TEST_NEAR(S(0, 2), S(2, 0), 0.0, "sigma_xz == sigma_zx");
   TEST_NEAR(S(0, 1), 42e6, 0.0, "sigma_xy value");
}

// ---------------------------------------------------------------------
// N-1  SquareOverstress resolver: indicator function across patches.
//      Inside a patch ⇒ amplitude == delta_tau; outside ⇒ 0.
// ---------------------------------------------------------------------
static void N_1_square_resolver_indicator()
{
   std::cout << "\n[N-1] ResolveSquareOverstress indicator behaviour\n";
   SquareOverstressSpec spec;
   spec.T_nuc_s = 0.5;
   {
      SquareOverstressPatch p;
      p.center_x_m = 0; p.center_y_m = 0; p.center_z_m = 7500;
      p.half_x_m = 1500; p.half_z_m = 1500;
      p.delta_tau_strike_pa = 11.6e6;
      spec.patches.push_back(p);
   }
   // 3 DOFs: inside, on edge, outside.
   Vector coords(9);
   coords(0) =    0.0; coords(1) = 0.0; coords(2) = 7500.0;   // inside
   coords(3) = 1500.0; coords(4) = 0.0; coords(5) = 7500.0;   // edge → inside
   coords(6) = 2000.0; coords(7) = 0.0; coords(8) = 7500.0;   // outside
   const auto p = ResolveSquareOverstress(spec, /*enabled=*/true, coords);
   TEST_ASSERT(p.amplitude_strike.Size() == 3, "amplitude_strike sized");
   TEST_NEAR(p.amplitude_strike(0), 11.6e6, 0.0, "centre amplitude");
   TEST_NEAR(p.amplitude_strike(1), 11.6e6, 0.0, "edge amplitude (inclusive)");
   TEST_NEAR(p.amplitude_strike(2),     0.0, 0.0, "outside amplitude");
   TEST_NEAR(p.radial(0), 1.0, 0.0, "centre radial = 1");
   TEST_NEAR(p.radial(2), 0.0, 0.0, "outside radial = 0");
}

// ---------------------------------------------------------------------
// N-2  SquareOverstress applicator telescopes to full delta_tau over
//      multiple sub-steps spanning [0, T_nuc_s].
// ---------------------------------------------------------------------
static void N_2_square_applicator_telescopes()
{
   std::cout << "\n[N-2] ApplySquareOverstressIncrement telescopes "
                "to full delta_tau\n";
   SquareOverstressPerDOFParams params;
   params.amplitude_dip.SetSize(1);    params.amplitude_dip    = 3.0;
   params.amplitude_strike.SetSize(1); params.amplitude_strike = 5.0;
   params.radial.SetSize(1);           params.radial           = 1.0;
   std::vector<DOFData> dof(1);
   const real_t T = 0.5;
   const int N = 5000;
   const real_t dt = T / N;
   for (int i = 0; i < N; ++i)
   {
      const real_t t_end = (i + 1) * dt;
      ApplySquareOverstressIncrement(dof, params, T, t_end, dt);
   }
   // After traversing [0, T], the accumulator should equal the full
   // delta_tau (smoothStep telescopes exactly).
   TEST_NEAR(dof[0].tau1_nuc, 3.0, 1e-9, "tau1_nuc telescoped to 3.0");
   TEST_NEAR(dof[0].tau2_nuc, 5.0, 1e-9, "tau2_nuc telescoped to 5.0");
}

// ---------------------------------------------------------------------
// N-3  SquareOverstress applicator early-returns past T_nuc_s.
// ---------------------------------------------------------------------
static void N_3_square_applicator_post_T_noop()
{
   std::cout << "\n[N-3] Applicator past T_nuc_s is a no-op\n";
   SquareOverstressPerDOFParams params;
   params.amplitude_dip.SetSize(1);    params.amplitude_dip    = 1.0;
   params.amplitude_strike.SetSize(1); params.amplitude_strike = 1.0;
   std::vector<DOFData> dof(1);
   dof[0].tau1_nuc = 99.0;
   dof[0].tau2_nuc = 99.0;
   // t_substep_end = T + 2 dt → applicator early-returns.
   ApplySquareOverstressIncrement(dof, params, /*T=*/0.5,
                                  /*t_end=*/1.0, /*dt=*/0.01);
   TEST_NEAR(dof[0].tau1_nuc, 99.0, 0.0, "tau1_nuc unchanged past T");
   TEST_NEAR(dof[0].tau2_nuc, 99.0, 0.0, "tau2_nuc unchanged past T");
}

int main(int, char**)
{
   std::cout << "Running test_spatial_stress_with_patches "
                "(stress patches + square nucleation)\n";
   P_1_empty_patches_degenerate();
   P_2_tpv205_central_patch();
   P_3_tpv205_three_patches_match_spec();
   P_4_last_match_wins();
   P_5_patch_override_is_symmetric();
   N_1_square_resolver_indicator();
   N_2_square_applicator_telescopes();
   N_3_square_applicator_post_T_noop();
   std::cout << "\n========================================\n";
   std::cout << "test_spatial_stress_with_patches: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
