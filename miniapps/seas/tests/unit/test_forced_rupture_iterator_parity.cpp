// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_forced_rupture_iterator_parity.cpp — Phase 2 of
// PLAN_TPV26_27_spatial_dyn_driver_2026-06-28.md §5.
//
// Guards the "round-6" fix in LinearSlipWeakeningIterator::StepOneQP_
// (dynamic/friction_substep_iterator.cpp).  That fix swaps the per-QP
// friction coefficient, when forced rupture is active, from
//
//     LSWFrictionCoefficient_TPV205(delta, mu_s, mu_d, d_c)          [plain]
// to
//     spatial::LSWFrictionCoefficient_ForcedRupture(delta, mu_s, mu_d, d_c,
//                            t_sub_end, d.T_forced_rupture, d.t0_decay_forced)
//
// so INTERIOR fault QPs see the same time-aware mu(delta, t) that the
// seam/inline LSW_ForcedRupture flux path already applies.
//
// These tests pin the mu contract the iterator now depends on:
//   (a) BYTE-EXACT reduction to the plain TPV205 formula at the
//       T_forced >= 1e8 "never forced" sentinel — this is what keeps every
//       existing LSW regression (TPV205 native + spatial) unchanged, and what
//       makes r >= rcrit DOFs behave as plain LSW inside a forced-rupture run;
//   (b) the f_2(t) forced front ramps mu from mu_s to mu_d over
//       [T_forced, T_forced + t0];
//   (c) mu = mu_s + (mu_d - mu_s) * MAX(f_1(delta), f_2(t)) — slip-weakening
//       and forced rupture compose by max, neither is lost;
//   (d) the t0 = 0 and barrier edge cases do not divide by zero / weaken.
//
// Standalone, header-only: both coefficients are inline free functions.

#include "mfem.hpp"

#include "../../dynamic/tpv205_friction.hpp"        // LSWFrictionCoefficient_TPV205
#include "../../spatial/code/spatial_friction.hpp"  // ..._ForcedRupture

#include "test_macros.hpp"

#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{
constexpr real_t kNeverForced = 1.0e9;
constexpr real_t kMuS = 0.18;   // TPV26/27 spec Part 4
constexpr real_t kMuD = 0.12;
constexpr real_t kDc  = 0.30;
}  // namespace

// (a) At the never-forced sentinel the time-aware helper must be BIT-EXACT
// with the plain TPV205 formula for every slip magnitude and every time.
// This is the byte-exact contract the round-6 gate relies on.
void Test_ForcedRupture_Reduces_To_LSW()
{
   const real_t deltas[] = { 0.0, 0.1 * kDc, kDc, 2.0 * kDc };
   const real_t times[]  = { 0.0, 1.0, 5.0, 100.0 };

   for (real_t delta : deltas)
   {
      const real_t plain =
         LSWFrictionCoefficient_TPV205(delta, kMuS, kMuD, kDc);
      for (real_t t : times)
      {
         const real_t forced = spatial::LSWFrictionCoefficient_ForcedRupture(
            delta, kMuS, kMuD, kDc, t, /*T_forced=*/kNeverForced,
            /*t0_decay=*/0.0);
         TEST_NEAR(forced, plain, 0.0,
                   "T_forced=1e9 is BIT-EXACT with plain LSW (delta="
                   + std::to_string(delta) + ", t=" + std::to_string(t) + ")");
      }
   }
}

// (b) The forced front: mu holds at mu_s until T_forced, ramps linearly to
// mu_d over [T_forced, T_forced + t0], then stays at mu_d.  A never-forced
// DOF at the same time stays at plain-LSW mu.
void Test_ForcedRupture_Front()
{
   const real_t T  = 1.0;
   const real_t t0 = 0.5;
   auto mu_at = [&](real_t t, real_t T_forced)
   {
      return spatial::LSWFrictionCoefficient_ForcedRupture(
         /*delta=*/0.0, kMuS, kMuD, kDc, t, T_forced, t0);
   };

   TEST_NEAR(mu_at(0.5, T), kMuS, 1e-15, "t < T_forced: mu = mu_s (unweakened)");
   TEST_NEAR(mu_at(1.0, T), kMuS, 1e-15, "t == T_forced: f_2 = 0, mu = mu_s");
   TEST_NEAR(mu_at(1.25, T), 0.15, 1e-15,
             "mid-ramp (f_2 = 0.5): mu = 0.15");
   TEST_NEAR(mu_at(1.5, T), kMuD, 1e-15,
             "t == T_forced + t0: f_2 = 1, mu = mu_d");
   TEST_NEAR(mu_at(2.0, T), kMuD, 1e-15, "t > T_forced + t0: mu stays mu_d");

   // A DOF outside rcrit (T = 1e9) never weakens by time.
   TEST_NEAR(mu_at(2.0, kNeverForced), kMuS, 0.0,
             "never-forced DOF stays at plain-LSW mu at the same time");
}

// (c) mu uses max(f_1, f_2): slip-weakening still works before the forced
// front arrives, and the forced front still works with zero slip.
void Test_Max_Of_f1_f2()
{
   // delta >= d_c  =>  f_1 = 1 dominates even at t << T_forced.
   TEST_NEAR(spatial::LSWFrictionCoefficient_ForcedRupture(
                kDc, kMuS, kMuD, kDc, /*t=*/0.0, /*T=*/1.0, /*t0=*/0.5),
             kMuD, 1e-15, "f_1 = 1 (delta = d_c) dominates before the front");

   // f_1 = 0.5 and f_2 = 0.5 -> max = 0.5 -> mu = 0.15.
   TEST_NEAR(spatial::LSWFrictionCoefficient_ForcedRupture(
                0.5 * kDc, kMuS, kMuD, kDc, /*t=*/1.25, /*T=*/1.0, /*t0=*/0.5),
             0.15, 1e-15, "f_1 == f_2 == 0.5 -> mu = 0.15");

   // f_1 = 0.5 but f_2 = 1 -> max = 1 -> mu = mu_d.
   TEST_NEAR(spatial::LSWFrictionCoefficient_ForcedRupture(
                0.5 * kDc, kMuS, kMuD, kDc, /*t=*/1.5, /*T=*/1.0, /*t0=*/0.5),
             kMuD, 1e-15, "f_2 = 1 overrides partial slip-weakening");
}

// (d) Edge cases: t0 = 0 is a step (no division by zero); a barrier DOF
// (mu_s at the 1e6 sentinel) short-circuits and never weakens.
void Test_Edge_Cases()
{
   // t0 = 0: mu_s strictly before T_forced, mu_d at/after it.
   TEST_NEAR(spatial::LSWFrictionCoefficient_ForcedRupture(
                0.0, kMuS, kMuD, kDc, /*t=*/0.999, /*T=*/1.0, /*t0=*/0.0),
             kMuS, 1e-15, "t0 = 0: mu = mu_s just before T_forced");
   TEST_NEAR(spatial::LSWFrictionCoefficient_ForcedRupture(
                0.0, kMuS, kMuD, kDc, /*t=*/1.0, /*T=*/1.0, /*t0=*/0.0),
             kMuD, 1e-15, "t0 = 0: step to mu_d at T_forced (no div-by-zero)");

   // Barrier sentinel (mu_s = 1e6): returns mu_s regardless of delta/time.
   const real_t barrier = 1.0e6;
   TEST_NEAR(spatial::LSWFrictionCoefficient_ForcedRupture(
                10.0 * kDc, barrier, kMuD, kDc, /*t=*/1e6, /*T=*/0.0,
                /*t0=*/0.5),
             barrier, 0.0, "barrier DOF short-circuits: mu = mu_s, never forced");
}

// Unify-plan Phase 2+4 (retires the former R-002 bounded-lag test).  Seam ==
// interior mu is STRUCTURAL since Phase 5: on the ADER substep path both face
// classes consume the friction iterator's per-substep buffer through ONE call
// site (ComputeFaultQPImposedStatesCanonical_ -> StepOneQP_), so both sides'
// mu comes from the SAME LSWFrictionCoefficient_ForcedRupture call at the
// SAME t_sub_end — that property is enforced by construction and guarded at
// integration level (test_shared_fault_substep_parity_np2, tet2x2 np=2==np=1
// harness).  What CAN regress at unit level is the composition itself, so pin
// it bit-exactly (REVIEW_phase4_5_unify R-101): the closed-form reference
// below mirrors spatial_friction.hpp's spec composition op-for-op (barrier
// aside — covered by Test_Edge_Cases).  Any conscious change to the mu
// formula must update this pin.  Includes the step ramp (t0 = 0), which is
// now safe at any rank count.
void Test_OneClock_Mu_Composition_Pinned()
{
   const real_t T = 1.0;
   const real_t deltas[] = { 0.0, 0.1 * kDc, 0.5 * kDc, kDc };
   const real_t times[]  = { 0.5, 1.0, 1.25, 1.5, 2.0 };
   const real_t t0s[]    = { 0.5, 0.0 };   // ramp and step ramp

   for (real_t t0 : t0s)
   {
      for (real_t delta : deltas)
      {
         for (real_t t : times)
         {
            // Independent closed-form reference: f_1 slip-weakening factor,
            // f_2 forced front, mu = mu_s + (mu_d - mu_s) * max(f_1, f_2).
            real_t f1;
            if (delta <= 0.0)      { f1 = 0.0; }
            else if (delta >= kDc) { f1 = 1.0; }
            else                   { f1 = delta / kDc; }
            real_t f2;
            if (t < T)             { f2 = 0.0; }
            else if (t < T + t0)   { f2 = (t - T) / t0; }
            else                   { f2 = 1.0; }
            const real_t factor = (f1 > f2) ? f1 : f2;
            const real_t mu_ref = kMuS + (kMuD - kMuS) * factor;

            const real_t mu_impl =
               spatial::LSWFrictionCoefficient_ForcedRupture(
                  delta, kMuS, kMuD, kDc, t, T, t0);
            TEST_NEAR(mu_impl, mu_ref, 0.0,
                      "one-clock mu contract: impl == closed-form spec "
                      "reference, bit-exact (the SINGLE composition both "
                      "seam and interior consume via the unified dispatch)");
         }
      }
   }
}

int main(int /*argc*/, char * /*argv*/[])
{
   std::cout << "=== test_forced_rupture_iterator_parity (Phase 2) ===\n";
   Test_ForcedRupture_Reduces_To_LSW();
   Test_ForcedRupture_Front();
   Test_Max_Of_f1_f2();
   Test_Edge_Cases();
   Test_OneClock_Mu_Composition_Pinned();
   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
