// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.0.0 Phase 1 §3.1c — Godunov identity under normal reversal.
//
// Checks whether the GodunovFlux::Interior is consistent under the
// reflection (n, L, R) → (-n, R, L):
//
//   F(+n, L, R) = −F(−n, R, L)    ⇔    F_plus + F_minus = 0
//
// Post-v9.0.0 Pelties-9 fix, NEITHER the interior-fault nor the
// shared-fault branch relies on this identity any longer — both now use
// the canonical normal `can_n` with per-side `A·T·Q` assembly (plan
// §14.2 and §14.3).  This test is retained as a unit invariant of
// GodunovFlux::Interior: a violation would indicate a drift in
// BuildFrame / BuildRotation under `n → −n`, which would silently
// break `flux_.Interior(n, Q, Q)` identity (§3.1e T-E) that the
// production per-side fix depends on.
//
// H-V9-I fires if this test FAILS beyond ULP.  The suspected
// mechanism (REVIEW R-1004) is `BuildFrame`'s Gram-Schmidt choice:
// `t1 = up × nor` — under `nor → −nor` we get `t1 → −t1`, then
// `t2 = nor × t1` gives `t2 → t2` (sign unchanged).  The resulting
// rotation matrices for ±n are not related by a simple sign flip,
// so the split-flux identity may break.
//
// Test fixture:
//   - Random but physically sane Q_A, Q_B (O(MPa) stress, O(m/s) velocity).
//   - Three non-trivial normals (axis-aligned and oblique).
//   - Component-wise assertion `|F+_c + F-_c| ≤ 1e-12 * (|F+_c| + |F-_c|)`.

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../config/tpv102_params.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static void TestGodunovIdentity(const real_t *nor,
                                const real_t *Q_A, const real_t *Q_B,
                                const std::string &label)
{
   std::cout << "  normal = (" << nor[0] << ", " << nor[1] << ", " << nor[2]
             << ")  [" << label << "]\n";

   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);

   real_t nor_plus[3]  = {nor[0],  nor[1],  nor[2]};
   real_t nor_minus[3] = {-nor[0], -nor[1], -nor[2]};

   real_t F_plus[NUM_STATE], F_minus[NUM_STATE];
   flux.Interior(nor_plus,  Q_A, Q_B, F_plus);
   flux.Interior(nor_minus, Q_B, Q_A, F_minus);

   // Component-wise conservation check.  Scale tolerance by the actual
   // magnitude of each flux component so near-zero components don't
   // trigger false fails, but also don't mask real drift on large ones.
   bool identity_ok = true;
   real_t worst_rel = 0.0;
   int worst_idx = -1;
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t scale = std::abs(F_plus[c]) + std::abs(F_minus[c]);
      real_t sum   = F_plus[c] + F_minus[c];
      if (scale > 0)
      {
         real_t rel = std::abs(sum) / scale;
         if (rel > worst_rel) { worst_rel = rel; worst_idx = c; }
         if (rel > 1e-12) { identity_ok = false; }
      }
      else
      {
         if (std::abs(sum) > 1e-12) { identity_ok = false; }
      }
   }

   std::cout << "    F_plus  = [";
   for (int c = 0; c < NUM_STATE; c++)
   {
      std::cout << F_plus[c] << (c < NUM_STATE-1 ? ", " : "");
   }
   std::cout << "]\n";
   std::cout << "    F_minus = [";
   for (int c = 0; c < NUM_STATE; c++)
   {
      std::cout << F_minus[c] << (c < NUM_STATE-1 ? ", " : "");
   }
   std::cout << "]\n";
   std::cout << "    F+ + F-  = [";
   for (int c = 0; c < NUM_STATE; c++)
   {
      std::cout << (F_plus[c] + F_minus[c])
                << (c < NUM_STATE-1 ? ", " : "");
   }
   std::cout << "]\n";
   std::cout << "    worst relative drift: " << worst_rel
             << " at component " << worst_idx << "\n";

   TEST_ASSERT(identity_ok,
               "Godunov identity F(+n,L,R) + F(-n,R,L) = 0 to ULP (" + label + ")");
}

static void FillRandomSaneState(real_t *Q, unsigned seed)
{
   std::srand(seed);
   auto r01 = []() { return (double)std::rand() / RAND_MAX; };
   // Physically sane magnitudes for elastic waves near a rupture front.
   Q[SXX] = (r01() - 0.5) * 2e6;   // O(MPa) pressure perturbation
   Q[SYY] = (r01() - 0.5) * 2e6;
   Q[SZZ] = (r01() - 0.5) * 2e6;
   Q[SXY] = (r01() - 0.5) * 4e6;   // O(MPa) shear
   Q[SYZ] = (r01() - 0.5) * 4e6;
   Q[SXZ] = (r01() - 0.5) * 4e6;
   Q[VX]  = (r01() - 0.5) * 2.0;   // O(m/s)
   Q[VY]  = (r01() - 0.5) * 2.0;
   Q[VZ]  = (r01() - 0.5) * 2.0;
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 v9.0.0 Phase 1 §3.1c Godunov identity test\n";
   std::cout << "========================================\n\n";

   real_t Q_A[NUM_STATE], Q_B[NUM_STATE];
   FillRandomSaneState(Q_A, 1);
   FillRandomSaneState(Q_B, 2);

   std::cout << "Test 3.1c: F(+n,L,R) + F(-n,R,L) = 0 (Godunov conservation)\n";

   // Case 1: axis-aligned normal (TPV102 vertical fault).
   real_t nor_y[3] = {0.0, -1.0, 0.0};
   TestGodunovIdentity(nor_y, Q_A, Q_B, "y-axis (TPV102 fault)");

   // Case 2: axis-aligned x.
   real_t nor_x[3] = {1.0, 0.0, 0.0};
   TestGodunovIdentity(nor_x, Q_A, Q_B, "x-axis");

   // Case 3: axis-aligned z (triggers the 'up'-fallback branch in BuildFrame).
   real_t nor_z[3] = {0.0, 0.0, 1.0};
   TestGodunovIdentity(nor_z, Q_A, Q_B, "z-axis (up-fallback branch)");

   // Case 4: generic oblique normal.
   real_t s = 1.0 / std::sqrt(3.0);
   real_t nor_obl[3] = {s, s, s};
   TestGodunovIdentity(nor_obl, Q_A, Q_B, "oblique (1,1,1)/sqrt(3)");

   // Case 5: tilted in xy-plane (no z component — independent of up-fallback).
   real_t nor_xy[3] = {0.6, 0.8, 0.0};
   TestGodunovIdentity(nor_xy, Q_A, Q_B, "tilted (0.6, 0.8, 0)");

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";

   return (num_failed > 0) ? 1 : 0;
}
