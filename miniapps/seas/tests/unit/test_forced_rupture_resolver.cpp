// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_forced_rupture_resolver.cpp — Phase 2 of
// PLAN_TPV26_27_spatial_dyn_driver_2026-06-28.md §5.
//
// Golden tests for spatial::ResolveForcedRupture / ForcedRuptureTime, the
// SCEC TPV26/27 forced-rupture time T(r) (spec Part 5; PLAN §1.4):
//
//   T(r) = r/(vr*Vs) + 0.081*rcrit/(vr*Vs) * (1/(1 - (r/rcrit)^2) - 1)   r <  rcrit
//        = 1e9  ("never forced" sentinel)                                r >= rcrit
//
// ADAPTED (not revived) from the archived `#if 0` block in
// tests/unit/test_spatial_friction_resolver.cpp:904-1203.  The T(r)
// expression and the r = {0,1000,2000,3000,4500} golden layout are reused
// verbatim; the CALL SHAPE is new:
//   * the removed `SpatialFrictionResolver::ResolveForcedRupture(nuc, dofs,
//     elem, mat, mesh)` is now the free function
//     `spatial::ResolveForcedRupture(spec, enabled, dof_coords_3d)`;
//   * the deleted `NucleationSpec::hypocenter_x_m / r_crit_m / t0_decay_s`
//     fields are now `spatial::ForcedRuptureSpec`;
//   * Vs is carried EXPLICITLY by the spec (no material lookup), so the
//     archived F-3 "centroid vs corner material Vs" regression (R-002) no
//     longer applies — there is no material evaluation to get wrong.
//
// Standalone: no MPI, no mesh — the resolver is pure arithmetic on coords.

#include "mfem.hpp"

#include "../../dynamic/spatial_nucleation.hpp"

#include "test_macros.hpp"

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas::spatial;

namespace
{

constexpr real_t kNeverForced = 1.0e9;

/// Runs `f` in a forked child and reports whether it aborted (MFEM_VERIFY).
/// Mirrors the `RunInChild_` idiom in tests/unit/test_friction_iterator_factory.cpp.
bool AbortsInChild(const std::function<void()>& f)
{
   std::fflush(stdout);
   const pid_t pid = fork();
   if (pid == 0)
   {
      std::fclose(stderr);   // silence the abort message
      f();
      std::_Exit(0);
   }
   int status = 0;
   waitpid(pid, &status, 0);
   return WIFSIGNALED(status) || (WIFEXITED(status) && WEXITSTATUS(status) != 0);
}

// Spec defaults (PLAN §1.4): rcrit = 4 km, Vs = 3464 m/s, vr = 0.7, t0 = 0.5 s.
ForcedRuptureSpec SpecDefaults()
{
   ForcedRuptureSpec s;
   s.hypocenter_x_m = 0.0;
   s.hypocenter_y_m = 0.0;
   s.hypocenter_z_m = 0.0;
   s.rcrit_m        = 4000.0;
   s.vs             = 3464.0;
   s.vr_factor      = 0.7;
   s.t0_s           = 0.5;
   return s;
}

// The analytic T(r) — written out independently of the implementation
// (verbatim from the archived test, :941-944).
real_t T_analytic(real_t r, real_t rcrit, real_t Vs, real_t vr)
{
   const real_t denom = 1.0 - (r / rcrit) * (r / rcrit);
   return r / (vr * Vs) + 0.081 * rcrit / (vr * Vs) * (1.0 / denom - 1.0);
}

}  // namespace

// r = {0, 1000, 2000, 3000} match the analytic T(r); r = 4500 > rcrit is the
// 1e9 sentinel; t0 is uniform; T(0) = 0; T is strictly increasing on [0,rcrit).
void Test_T_of_r()
{
   const ForcedRuptureSpec spec = SpecDefaults();
   const int N = 5;
   const real_t rs[N] = { 0.0, 1000.0, 2000.0, 3000.0, 4500.0 };

   Vector coords(3 * N);
   for (int i = 0; i < N; ++i)
   {
      coords(3 * i + 0) = rs[i];   // along-strike offset from hypocenter
      coords(3 * i + 1) = 0.0;
      coords(3 * i + 2) = 0.0;
   }

   const ForcedRupturePerDOFParams p =
      ResolveForcedRupture(spec, /*enabled=*/true, coords);

   TEST_ASSERT(p.T_forced_s.Size() == N && p.t0_decay_s.Size() == N,
               "ResolveForcedRupture sizes == num fault DOFs");

   for (int i = 0; i < 4; ++i)
   {
      const real_t T_expect =
         T_analytic(rs[i], spec.rcrit_m, spec.vs, spec.vr_factor);
      TEST_NEAR(p.T_forced_s(i), T_expect, 1e-6,
                "T(r=" + std::to_string(int(rs[i])) + ") matches analytic");
   }

   TEST_NEAR(p.T_forced_s(0), 0.0, 1e-12, "T(0) = 0 (forced from t=0)");
   TEST_NEAR(p.T_forced_s(4), kNeverForced, 0.0,
             "T(r=4500 >= rcrit) = 1e9 sentinel");

   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(p.t0_decay_s(i), 0.5, 0.0, "t0_decay = 0.5 s at every DOF");
   }

   TEST_ASSERT(p.T_forced_s(0) < p.T_forced_s(1)
               && p.T_forced_s(1) < p.T_forced_s(2)
               && p.T_forced_s(2) < p.T_forced_s(3),
               "T(r) strictly increasing on [0, rcrit)");
}

// The resolver must honour the sibling-resolver contract: disabled ⇒
// zero-sized Vectors (the driver then keeps the never-forced sentinel).
void Test_Disabled_Returns_Empty()
{
   Vector coords(3 * 3);  coords = 0.0;
   const ForcedRupturePerDOFParams p =
      ResolveForcedRupture(SpecDefaults(), /*enabled=*/false, coords);
   TEST_ASSERT(p.T_forced_s.Size() == 0 && p.t0_decay_s.Size() == 0,
               "disabled ⇒ zero-sized per-DOF Vectors");
}

// r is measured in the (x, z) fault plane about the hypocenter: a DOF at the
// hypocenter gets T = 0; a DOF at in-plane distance >= rcrit is never forced.
void Test_Hypocenter_Offset_And_Plane()
{
   ForcedRuptureSpec spec = SpecDefaults();
   spec.hypocenter_x_m = -5000.0;
   spec.hypocenter_z_m = -10000.0;

   const int N = 2;
   Vector coords(3 * N);
   // DOF 0: exactly at the hypocenter -> r = 0.
   coords(0) = -5000.0; coords(1) = 0.0; coords(2) = -10000.0;
   // DOF 1: hypocenter + (3000 in x, 4000 in z) -> r = 5000 >= rcrit.
   coords(3) = -2000.0; coords(4) = 0.0; coords(5) = -6000.0;

   const ForcedRupturePerDOFParams p =
      ResolveForcedRupture(spec, /*enabled=*/true, coords);

   TEST_NEAR(p.T_forced_s(0), 0.0, 1e-12, "T = 0 exactly at the hypocenter");
   TEST_NEAR(p.T_forced_s(1), kNeverForced, 0.0,
             "in-plane r = 5000 >= rcrit -> never-forced sentinel");
}

// R-004: `r` ignores y, so a curved / fault-normal-offset fault would have its
// entire forced-rupture front silently mis-placed.  The resolver must REJECT
// such a fault rather than accept it.  (Before the fix this case was silently
// treated as r = 0, i.e. "forced at t = 0".)
void Test_R004_NonPlanar_Fault_Rejected()
{
   ForcedRuptureSpec spec = SpecDefaults();   // hypocenter_y_m = 0

   Vector planar(3);
   planar(0) = 0.0; planar(1) = 0.0; planar(2) = 0.0;
   TEST_ASSERT(!AbortsInChild([&] {
                  (void) ResolveForcedRupture(spec, true, planar); }),
               "planar y=0 fault DOF is accepted");

   Vector off_plane(3);
   off_plane(0) = 0.0; off_plane(1) = 3000.0; off_plane(2) = 0.0;
   TEST_ASSERT(AbortsInChild([&] {
                  (void) ResolveForcedRupture(spec, true, off_plane); }),
               "fault DOF 3 km off the y-plane is REJECTED (curved fault guard)");
}

// r -> rcrit^- makes T blow up; the resolver clamps to the 1e9 sentinel so no
// non-finite forced-rupture time ever reaches DOFData.
void Test_Asymptote_Clamped_Finite()
{
   const ForcedRuptureSpec spec = SpecDefaults();
   const real_t r_near = spec.rcrit_m * (1.0 - 1e-14);
   const real_t T = ForcedRuptureTime(r_near, spec.rcrit_m, spec.vs,
                                      spec.vr_factor);
   TEST_ASSERT(std::isfinite(T), "T(r -> rcrit^-) is finite (not inf/NaN)");
   TEST_ASSERT(T <= kNeverForced, "T(r -> rcrit^-) clamped to <= 1e9 sentinel");
   TEST_NEAR(ForcedRuptureTime(spec.rcrit_m, spec.rcrit_m, spec.vs,
                               spec.vr_factor),
             kNeverForced, 0.0, "T(r == rcrit) = 1e9 sentinel");
}

int main(int /*argc*/, char * /*argv*/[])
{
   std::cout << "=== test_forced_rupture_resolver (Phase 2) ===\n";
   Test_T_of_r();
   Test_Disabled_Returns_Empty();
   Test_Hypocenter_Offset_And_Plane();
   Test_R004_NonPlanar_Fault_Rejected();
   Test_Asymptote_Clamped_Finite();
   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
