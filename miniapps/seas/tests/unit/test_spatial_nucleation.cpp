// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_nucleation.cpp — Phase N of
// safs/project_7.0_alternative/document/05_18_2026/PLAN_first_safs_run.md.
//
// 12 unit tests T-N01..T-N12 covering:
//   - SmoothStep boundary values + telescope identity
//   - SmoothStepIncrement dt finiteness assertion
//   - GaussianFactorFaceLocal centre value + half-amplitude radius
//     + anisotropy
//   - ResolveGradualOverstress disabled / 5-DOF fixture
//   - ApplyGradualOverstressIncrement telescope + post-T_nuc no-op
//   - Tpv205SubStepIterator callback overload invocation contract.

#include "mfem.hpp"

#include "../../dynamic/spatial_nucleation.hpp"
#include "../../dynamic/tpv205_substep_iterator.hpp"
#include "../../dynamic/fault_face_flux.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <utility>
#include <vector>

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
// T-N01: SmoothStep boundary values
// =====================================================================
static void T_N01_smoothstep_boundary()
{
   std::cout << "\n[T-N01] SmoothStep boundary values\n";
   TEST_NEAR(SmoothStep(0.0,   1.0), 0.0, 0.0, "SmoothStep(0, 1) = 0");
   TEST_NEAR(SmoothStep(1.0,   1.0), 1.0, 0.0, "SmoothStep(1, 1) = 1");
   TEST_NEAR(SmoothStep(2.0,   1.0), 1.0, 0.0, "SmoothStep(2, 1) = 1");
   TEST_NEAR(SmoothStep(-0.5,  1.0), 0.0, 0.0, "SmoothStep(-0.5, 1) = 0");
   // Monotonicity check at an interior point.
   TEST_ASSERT(SmoothStep(0.25, 1.0) < SmoothStep(0.5, 1.0),
               "SmoothStep monotone increasing on (0, t0)");
   TEST_ASSERT(SmoothStep(0.5,  1.0) < SmoothStep(0.75, 1.0),
               "SmoothStep monotone increasing on (0, t0)");
}

// =====================================================================
// T-N02: SmoothStepIncrement uniform telescope
// =====================================================================
static void T_N02_smoothstep_telescope_uniform()
{
   std::cout << "\n[T-N02] SmoothStepIncrement uniform partition telescopes\n";
   const real_t t0 = 1.0;
   const int    n  = 100;
   const real_t dt = t0 / n;
   real_t sum = 0.0;
   for (int k = 1; k <= n; ++k)
   {
      sum += SmoothStepIncrement(k * dt, dt, t0);
   }
   TEST_NEAR(sum, 1.0, 1e-12,
             "Σ_{k=1..100} SmoothStepIncrement(k·dt, dt, 1) = 1");
}

// =====================================================================
// T-N03: SmoothStepIncrement non-uniform telescope
// =====================================================================
static void T_N03_smoothstep_telescope_nonuniform()
{
   std::cout << "\n[T-N03] SmoothStepIncrement non-uniform partition telescopes\n";
   const real_t t0 = 1.0;
   // Random positive partition summing to 1.
   std::mt19937 rng(12345);
   std::uniform_real_distribution<real_t> u(0.01, 0.05);
   std::vector<real_t> dts;
   real_t total = 0.0;
   while (total < 1.0)
   {
      real_t d = u(rng);
      if (total + d > 1.0) { d = 1.0 - total; }
      dts.push_back(d);
      total += d;
   }
   real_t sum = 0.0;
   real_t t   = 0.0;
   for (real_t d : dts)
   {
      t   += d;
      sum += SmoothStepIncrement(t, d, t0);
   }
   TEST_NEAR(sum, 1.0, 1e-12,
             "non-uniform partition telescopes to 1");
}

// =====================================================================
// T-N04: SmoothStepIncrement dt finiteness MFEM_ASSERT
// =====================================================================
static bool RunInChild_(const std::function<void()>& f)
{
   pid_t pid = fork();
   if (pid == 0)
   {
      // Silence the abort message.
      std::fclose(stderr);
      f();
      std::_Exit(0);
   }
   int status = 0;
   waitpid(pid, &status, 0);
   return WIFSIGNALED(status) || (WIFEXITED(status) && WEXITSTATUS(status) != 0);
}

static void T_N04_smoothstep_increment_dt_assert()
{
   std::cout << "\n[T-N04] SmoothStepIncrement dt finite/positive assertion\n";
#ifdef MFEM_DEBUG
   const bool zero_aborted = RunInChild_([](){
      (void)SmoothStepIncrement(0.5, 0.0, 1.0);
   });
   TEST_ASSERT(zero_aborted, "dt = 0 aborts");
   const bool nan_aborted = RunInChild_([](){
      (void)SmoothStepIncrement(0.5,
                                std::numeric_limits<real_t>::quiet_NaN(),
                                1.0);
   });
   TEST_ASSERT(nan_aborted, "dt = NaN aborts");
#else
   std::cout << "  SKIP (MFEM_DEBUG not defined; MFEM_ASSERT compiled out)\n";
#endif
}

// =====================================================================
// T-N05: GaussianFactor at centre = 1
// =====================================================================
static void T_N05_gaussian_centre_value()
{
   std::cout << "\n[T-N05] GaussianFactorFaceLocal(centre) = 1\n";
   const real_t dof_xyz[3] = { 10.0, 0.0, 0.0 };
   const real_t center[3]  = { 10.0, 0.0, 0.0 };
   const real_t dip[3]     = { 0.0, 1.0, 0.0 };
   const real_t strike[3]  = { 0.0, 0.0, 1.0 };
   const real_t F = GaussianFactorFaceLocal(dof_xyz, dip, strike, center,
                                            1000.0, 1000.0);
   TEST_NEAR(F, 1.0, 1e-15, "F at centre = 1");
}

// =====================================================================
// T-N06: GaussianFactor half-amplitude radius
// =====================================================================
static void T_N06_gaussian_half_amplitude()
{
   std::cout << "\n[T-N06] GaussianFactorFaceLocal at r = radius·sqrt(ln 2) = 0.5\n";
   const real_t radius = 1000.0;
   const real_t offset = radius * std::sqrt(std::log(2.0));   // F = exp(-ln2) = 0.5
   // Place DOF offset along the dip direction.
   const real_t dip[3]     = { 0.0, 1.0, 0.0 };
   const real_t strike[3]  = { 0.0, 0.0, 1.0 };
   const real_t center[3]  = { 0.0, 0.0, 0.0 };
   const real_t dof_xyz[3] = { 0.0, offset, 0.0 };
   const real_t F = GaussianFactorFaceLocal(dof_xyz, dip, strike, center,
                                            radius, radius);
   TEST_NEAR(F, 0.5, 1e-12, "F at half-amplitude radius along dip = 0.5");
}

// =====================================================================
// T-N07: anisotropy — radius_dip ≠ radius_strike
// =====================================================================
static void T_N07_gaussian_anisotropy()
{
   std::cout << "\n[T-N07] GaussianFactorFaceLocal anisotropy "
                "(radius_dip = 2 · radius_strike)\n";
   const real_t r_dip    = 2000.0;
   const real_t r_strike = 1000.0;
   const real_t dip[3]     = { 0.0, 1.0, 0.0 };
   const real_t strike[3]  = { 0.0, 0.0, 1.0 };
   const real_t center[3]  = { 0.0, 0.0, 0.0 };
   // Same offset distance, once along dip, once along strike.
   const real_t off = 1000.0;
   const real_t dof_dip   [3] = { 0.0, off,  0.0 };
   const real_t dof_strike[3] = { 0.0, 0.0,  off };
   const real_t F_dip    = GaussianFactorFaceLocal(dof_dip,    dip, strike,
                                                   center, r_dip, r_strike);
   const real_t F_strike = GaussianFactorFaceLocal(dof_strike, dip, strike,
                                                   center, r_dip, r_strike);
   // F_dip = exp(-(1000/2000)^2) = exp(-0.25)
   // F_strike = exp(-(1000/1000)^2) = exp(-1)
   TEST_NEAR(F_dip,    std::exp(-0.25), 1e-12, "F along dip uses radius_dip");
   TEST_NEAR(F_strike, std::exp(-1.0),  1e-12, "F along strike uses radius_strike");
   TEST_ASSERT(F_dip > F_strike,
               "F_dip > F_strike for larger dip radius at same offset");
}

// =====================================================================
// T-N08: ResolveGradualOverstress disabled returns empty
// =====================================================================
static void T_N08_resolve_disabled_empty()
{
   std::cout << "\n[T-N08] ResolveGradualOverstress(enabled=false) returns empty\n";
   GradualOverstressSpec spec;
   spec.radius_dip_m = spec.radius_strike_m = 1000.0;
   spec.T_nuc_s = 1.0;
   const int N = 5;
   Vector dofs(3 * N); dofs = 0.0;
   DenseMatrix basis(9, N); basis = 0.0;
   const auto p = ResolveGradualOverstress(spec, /*enabled=*/false,
                                           dofs, basis);
   TEST_ASSERT(p.amplitude_dip.Size()    == 0, "amplitude_dip empty");
   TEST_ASSERT(p.amplitude_strike.Size() == 0, "amplitude_strike empty");
   TEST_ASSERT(p.radial.Size()           == 0, "radial empty");
}

// =====================================================================
// T-N09: ResolveGradualOverstress 5-DOF fixture (centre at DOF 2)
//        Validates the (9, N) column-major basis layout (R-001 fix).
// =====================================================================
static void T_N09_resolve_5dof_fixture()
{
   std::cout << "\n[T-N09] ResolveGradualOverstress 5-DOF symmetric fixture\n";
   const int N = 5;
   // 5 DOFs spaced along +y direction, centred on DOF 2.
   Vector dofs(3 * N);
   for (int i = 0; i < N; ++i)
   {
      dofs(3 * i + 0) = 100.0;
      dofs(3 * i + 1) = static_cast<real_t>(i - 2) * 200.0;  // -400..+400
      dofs(3 * i + 2) = -50.0;
   }
   // All DOFs share the same basis: normal=+x, dip=+y, strike=+z.
   DenseMatrix basis(9, N);
   basis = 0.0;
   for (int i = 0; i < N; ++i)
   {
      basis(0, i) = 1.0;   // normal x
      basis(4, i) = 1.0;   // dip   y
      basis(8, i) = 1.0;   // strike z
   }
   GradualOverstressSpec spec;
   spec.center_x_m = 100.0;
   spec.center_y_m = 0.0;
   spec.center_z_m = -50.0;
   spec.radius_dip_m = spec.radius_strike_m = 500.0;
   spec.delta_tau_dip_pa    = 1.0e6;
   spec.delta_tau_strike_pa = 2.0e7;
   spec.T_nuc_s = 1.0;

   const auto p = ResolveGradualOverstress(spec, /*enabled=*/true,
                                           dofs, basis);
   TEST_ASSERT(p.amplitude_dip.Size()    == N, "amplitude_dip sized N");
   TEST_ASSERT(p.amplitude_strike.Size() == N, "amplitude_strike sized N");
   TEST_ASSERT(p.radial.Size()           == N, "radial sized N");

   // Centre DOF: F = 1 exactly.
   TEST_NEAR(p.radial(2),           1.0,                1e-15,
             "F at centre DOF 2 = 1");
   TEST_NEAR(p.amplitude_strike(2), spec.delta_tau_strike_pa, 1e-6,
             "amplitude_strike at centre = delta_tau_strike_pa");
   TEST_NEAR(p.amplitude_dip(2),    spec.delta_tau_dip_pa,    1e-6,
             "amplitude_dip at centre = delta_tau_dip_pa");
   // Symmetry around the centre DOF.
   TEST_NEAR(p.radial(0), p.radial(4), 1e-15, "radial symmetry 0 vs 4");
   TEST_NEAR(p.radial(1), p.radial(3), 1e-15, "radial symmetry 1 vs 3");
   // Monotonic decrease from centre outward.
   TEST_ASSERT(p.radial(2) > p.radial(1), "F decreases off centre");
   TEST_ASSERT(p.radial(1) > p.radial(0), "F decreases farther off");
}

// =====================================================================
// T-N10: ApplyGradualOverstressIncrement telescopes per DOF over [0, T_nuc]
// =====================================================================
static void T_N10_accumulator_telescope()
{
   std::cout << "\n[T-N10] ApplyGradualOverstressIncrement telescopes to full target\n";
   const int N = 3;
   GradualOverstressPerDOFParams params;
   params.amplitude_dip.SetSize(N);
   params.amplitude_strike.SetSize(N);
   params.radial.SetSize(N);
   for (int i = 0; i < N; ++i)
   {
      params.amplitude_dip(i)    = 1.0e6 * (i + 1);    // 1e6, 2e6, 3e6
      params.amplitude_strike(i) = 2.0e6 * (i + 1);    // 2e6, 4e6, 6e6
      params.radial(i)           = 1.0;
   }
   std::vector<DOFData> dof_data(N);
   const real_t T_nuc = 1.0;
   const int    nsub  = 100;
   const real_t dt    = T_nuc / nsub;
   real_t t = 0.0;
   for (int k = 1; k <= nsub; ++k)
   {
      t += dt;
      ApplyGradualOverstressIncrement(dof_data, params, T_nuc, t, dt);
   }
   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(dof_data[i].tau1_nuc, params.amplitude_dip(i),    1e-3,
                "tau1_nuc telescopes to amplitude_dip");
      TEST_NEAR(dof_data[i].tau2_nuc, params.amplitude_strike(i), 1e-3,
                "tau2_nuc telescopes to amplitude_strike");
      TEST_NEAR(dof_data[i].sigma_n_nuc, 0.0, 0.0,
                "sigma_n_nuc untouched (mechanism does not perturb σ_n)");
   }
}

// =====================================================================
// T-N11: post-T_nuc no-op (accumulator does not move)
// =====================================================================
static void T_N11_accumulator_post_tnuc_noop()
{
   std::cout << "\n[T-N11] ApplyGradualOverstressIncrement post-T_nuc no-op\n";
   const int N = 1;
   GradualOverstressPerDOFParams params;
   params.amplitude_dip.SetSize(N);
   params.amplitude_strike.SetSize(N);
   params.radial.SetSize(N);
   params.amplitude_dip(0)    = 3.0e6;
   params.amplitude_strike(0) = 4.0e6;
   params.radial(0)           = 1.0;

   std::vector<DOFData> dof_data(N);
   // Telescope to full target.
   const real_t T_nuc = 1.0;
   const int    nsub  = 50;
   const real_t dt    = T_nuc / nsub;
   real_t t = 0.0;
   for (int k = 1; k <= nsub; ++k)
   {
      t += dt;
      ApplyGradualOverstressIncrement(dof_data, params, T_nuc, t, dt);
   }
   const real_t tau1_after = dof_data[0].tau1_nuc;
   const real_t tau2_after = dof_data[0].tau2_nuc;
   // Now advance well past T_nuc.
   for (int k = 0; k < 20; ++k)
   {
      t += dt;
      ApplyGradualOverstressIncrement(dof_data, params, T_nuc, t, dt);
   }
   TEST_NEAR(dof_data[0].tau1_nuc, tau1_after, 0.0,
             "tau1_nuc unchanged after T_nuc");
   TEST_NEAR(dof_data[0].tau2_nuc, tau2_after, 0.0,
             "tau2_nuc unchanged after T_nuc");
}

// =====================================================================
// T-N12: Tpv205SubStepIterator callback overload — invocation contract
//        (the original no-callback overload also routes through this
//        with a no-op callback, R-N-005 — when the dof_data is empty
//        the iterator validation is still exercised end-to-end).
// =====================================================================
static void T_N12_iterator_callback_invocation()
{
   std::cout << "\n[T-N12] Tpv205SubStepIterator callback overload contract\n";
   // O = 4 sub-steps.  Use an empty fault (n = 0); the per-QP loop is
   // skipped but the callback should still fire once per sub-step.
   FaultFaceFlux flux(2670.0, 6000.0, 3500.0);
   Tpv205SubStepIterator iter(flux);
   const int O = 4;
   const real_t dt_macro = 1.0;
   std::vector<real_t> deltaT(O, dt_macro / O);
   std::vector<real_t> weights(O, 1.0 / O);
   iter.SetSubSteps(deltaT, weights);

   std::vector<DOFData>     dof_data;          // empty fault
   std::vector<Vector>      fault_coords;
   std::vector<std::vector<real_t>> Qp(O), Qm(O);
   real_t I_imp_plus  = 0.0;   // single dummy slot — n_words = 0
   real_t I_imp_minus = 0.0;

   std::vector<std::pair<real_t, real_t>> calls;
   auto cb = [&calls](real_t t_end, real_t dt_sub)
   {
      calls.emplace_back(t_end, dt_sub);
   };

   const real_t t_macro_start = 7.0;
   iter.AdvanceWithSubStepStates(dof_data, fault_coords, Qp, Qm,
                                 dt_macro, t_macro_start,
                                 &I_imp_plus, &I_imp_minus, cb);
   TEST_ASSERT(static_cast<int>(calls.size()) == O,
               "callback fired exactly O = 4 times");
   for (int k = 0; k < O; ++k)
   {
      const real_t expected_t_end =
         t_macro_start + (k + 1) * (dt_macro / O);
      TEST_NEAR(calls[k].first,  expected_t_end, 1e-12,
                "callback t_substep_end matches monotone schedule");
      TEST_NEAR(calls[k].second, dt_macro / O,   1e-12,
                "callback dt_substep matches sub-step size");
   }
}

// =====================================================================
// T-N13: GradualOverstressCompactCircular — compact-support bell at
// r ≥ R returns 0; F(0) = 1; symmetric in (dip, strike).
// =====================================================================
static void T_N13_compact_circular_compact_support()
{
   std::cout << "\n[T-N13] ResolveGradualOverstressCompactCircular compact "
             << "support + F(0) = 1\n";
   const int N = 5;
   // 5 DOFs along the strike axis at distances 0, 0.5R, R, 1.5R, 2R
   // from the centre.  R = 1000 m.
   const real_t R = 1000.0;
   Vector dofs(3 * N);
   const real_t offsets[N] = { 0.0, 500.0, 1000.0, 1500.0, 2000.0 };
   for (int i = 0; i < N; ++i)
   {
      dofs(3 * i + 0) = 0.0;
      dofs(3 * i + 1) = 0.0;
      dofs(3 * i + 2) = offsets[i];   // along strike (basis row 6..8 below)
   }
   DenseMatrix basis(9, N);
   basis = 0.0;
   for (int i = 0; i < N; ++i)
   {
      basis(0, i) = 1.0;   // normal  = +x
      basis(4, i) = 1.0;   // dip     = +y
      basis(8, i) = 1.0;   // strike  = +z
   }
   GradualOverstressCompactCircularSpec spec;
   spec.center_x_m = 0.0;
   spec.center_y_m = 0.0;
   spec.center_z_m = 0.0;
   spec.radius_m   = R;
   spec.delta_tau_strike_pa = 25.0e6;
   spec.T_nuc_s    = 1.0;

   const auto p = ResolveGradualOverstressCompactCircular(
      spec, /*enabled=*/true, dofs, basis);
   TEST_ASSERT(p.amplitude_strike.Size() == N, "amplitude_strike sized N");
   TEST_NEAR(p.radial(0), 1.0,             1e-15, "F(r=0) = 1");
   TEST_ASSERT(p.radial(1) > 0.0 && p.radial(1) < 1.0,
               "F(r=0.5R) in (0, 1)");
   TEST_NEAR(p.radial(2), 0.0,             0.0,   "F(r=R)  = 0 (compact)");
   TEST_NEAR(p.radial(3), 0.0,             0.0,   "F(r>R)  = 0");
   TEST_NEAR(p.radial(4), 0.0,             0.0,   "F(r=2R) = 0");
   TEST_NEAR(p.amplitude_strike(0),
             spec.delta_tau_strike_pa, 1e-6,
             "amplitude_strike at centre = Δτ_strike_pa");
   TEST_NEAR(p.amplitude_strike(2), 0.0, 0.0,
             "amplitude_strike at r=R = 0");
}

// =====================================================================
// T-N14: ApplyGradualOverstressCompactCircularIncrement telescopes
// to the full target over [0, T_nuc] (smoothStep ramp).
// =====================================================================
static void T_N14_compact_circular_telescope()
{
   std::cout << "\n[T-N14] ApplyGradualOverstressCompactCircularIncrement "
             << "telescopes to full target\n";
   const int N = 2;
   GradualOverstressCompactCircularPerDOFParams params;
   params.amplitude_dip.SetSize(N);
   params.amplitude_strike.SetSize(N);
   params.radial.SetSize(N);
   for (int i = 0; i < N; ++i)
   {
      params.amplitude_dip(i)    = 0.0;
      params.amplitude_strike(i) = 25.0e6 * (i + 1);   // 25e6, 50e6
      params.radial(i)           = 1.0;
   }
   std::vector<DOFData> dof_data(N);
   const real_t T_nuc = 1.0;
   const int    nsub  = 200;
   const real_t dt    = T_nuc / nsub;
   real_t t = 0.0;
   for (int k = 1; k <= nsub; ++k)
   {
      t += dt;
      ApplyGradualOverstressCompactCircularIncrement(dof_data, params,
                                                     T_nuc, t, dt);
   }
   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(dof_data[i].tau2_nuc, params.amplitude_strike(i), 1e-3,
                "tau2_nuc telescopes to amplitude_strike");
      TEST_NEAR(dof_data[i].tau1_nuc, 0.0, 0.0,
                "tau1_nuc untouched (Δτ_dip = 0)");
   }
}

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase N test_spatial_nucleation\n";
   T_N01_smoothstep_boundary();
   T_N02_smoothstep_telescope_uniform();
   T_N03_smoothstep_telescope_nonuniform();
   T_N04_smoothstep_increment_dt_assert();
   T_N05_gaussian_centre_value();
   T_N06_gaussian_half_amplitude();
   T_N07_gaussian_anisotropy();
   T_N08_resolve_disabled_empty();
   T_N09_resolve_5dof_fixture();
   T_N10_accumulator_telescope();
   T_N11_accumulator_post_tnuc_noop();
   T_N12_iterator_callback_invocation();
   T_N13_compact_circular_compact_support();
   T_N14_compact_circular_telescope();

   std::cout << "\n========================================\n";
   std::cout << "Phase N test_spatial_nucleation: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}
