// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 R-I06-round3 R-006: direct unit tests for the total-Q
// free-surface BC variants (FreeSurfaceTotal, FreeSurfaceGodunovTotal).
//
// Previously these variants were only exercised transitively through
// test_free_surface_godunov_driver.cpp (which inherits the dispatch
// logic from wave_operator.inl and buries the failure mode in a large
// simulation).  This test exercises the flux kernels DIRECTLY on
// synthetic inputs, making it easy to localize regressions.
//
// Gates:
//   1. FreeSurfaceTotal(Q_bg=0) == FreeSurface                    (~10 ULP)
//   2. FreeSurfaceGodunovTotal(Q_bg=0) == FreeSurfaceGodunov       (~10 ULP)
//   3. FreeSurfaceTotal(Q_self=Q_bg, tilted normal) = Interior flux
//      (equilibrium preserved on a tilted free surface under total-Q)
//   4. FreeSurfaceGodunovTotal(Q_self=Q_bg, tilted normal) = Interior flux
//
// Gates 3+4 are the load-bearing R-I06-005 check: under total-Q, a
// uniform pre-stress field must NOT radiate from tilted free-surface
// faces.  On TPV102's horizontal z=0 surface the gamma-mirror already
// produces zero flux on uniform pre-stress (sigma_pre . n_z = 0), so
// the tilted-normal case is where the total-aware variant actually
// deviates from the base variant.

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

static constexpr double kEps = 2.220446049250313e-16;

static real_t MaxNorm(const real_t v[NUM_STATE])
{
   real_t m = 0;
   for (int c = 0; c < NUM_STATE; c++) { m = std::max(m, std::abs(v[c])); }
   return m;
}

static real_t MaxDiff(const real_t a[NUM_STATE], const real_t b[NUM_STATE])
{
   real_t m = 0;
   for (int c = 0; c < NUM_STATE; c++) { m = std::max(m, std::abs(a[c] - b[c])); }
   return m;
}

int main()
{
   std::cout << "\n=== R-I06-round3 R-006: FreeSurfaceTotal / "
             << "FreeSurfaceGodunovTotal direct tests ===\n";

   const real_t lambda = 32.04e9;
   const real_t mu     = 32.04e9;
   const real_t rho    = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   // Axis-aligned normal first — used for the Q_bg=0 equivalence gates.
   real_t nor_axis[3] = {0.0, 0.0, 1.0};

   // ------------------------------------------------------------------
   // Gate 1: FreeSurfaceTotal(Q_bg=0) == FreeSurface
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate 1: Q_bg = 0 ⇒ FreeSurfaceTotal ≡ FreeSurface --\n";

      // A representative Q_self with non-zero stresses and velocities.
      real_t Q_self[NUM_STATE] = {1.0e7, 2.0e7, -1.5e7, 3.0e6, -2.5e6,
                                  4.0e6, 1.0, -2.0, 0.5};
      real_t Q_bg_zero[NUM_STATE] = {0};

      real_t F_base[NUM_STATE], F_total[NUM_STATE];
      flux.FreeSurface     (nor_axis, Q_self, F_base);
      flux.FreeSurfaceTotal(nor_axis, Q_self, Q_bg_zero, F_total);

      const real_t err   = MaxDiff(F_base, F_total);
      const real_t scale = std::max(MaxNorm(F_base), MaxNorm(F_total));
      // Algebraically the two should be BIT-EQUAL: with Q_bg=0, Q_pert =
      // Q_self, and the rotated ghost = 0 + gamma*Q_self_rot = gamma*Q_self_rot,
      // which is exactly what FreeSurface constructs.  Allow a small
      // budget for mat-vec reordering.
      TEST_LE(err, 16 * kEps * scale,
              "FreeSurfaceTotal(Q_bg=0) matches FreeSurface to 16 ULP");
   }

   // ------------------------------------------------------------------
   // Gate 2: FreeSurfaceGodunovTotal(Q_bg=0) == FreeSurfaceGodunov
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate 2: Q_bg = 0 ⇒ FreeSurfaceGodunovTotal ≡ "
                << "FreeSurfaceGodunov --\n";

      real_t Q_self[NUM_STATE] = {5.0e7, -1.0e7, 2.0e6, 1.0e7, 0.5e6,
                                  -2.0e7, 0.3, 1.5, -0.7};
      real_t Q_bg_zero[NUM_STATE] = {0};

      real_t F_base[NUM_STATE], F_total[NUM_STATE];
      flux.FreeSurfaceGodunov     (nor_axis, Q_self, F_base);
      flux.FreeSurfaceGodunovTotal(nor_axis, Q_self, Q_bg_zero, F_total);

      const real_t err   = MaxDiff(F_base, F_total);
      const real_t scale = std::max(MaxNorm(F_base), MaxNorm(F_total));
      TEST_LE(err, 16 * kEps * scale,
              "FreeSurfaceGodunovTotal(Q_bg=0) matches FreeSurfaceGodunov "
              "to 16 ULP");
   }

   // ------------------------------------------------------------------
   // Gate 3: Q_self = Q_bg on a TILTED normal — FreeSurfaceTotal flux
   //         must equal Interior(Q_bg, Q_bg) = A_n·Q_bg.  This is the
   //         core R-I06-005 property: uniform pre-stress under total-Q
   //         does not radiate from tilted free surfaces.
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate 3: Q_self = Q_bg on tilted normal ⇒ "
                << "FreeSurfaceTotal = Interior (no radiation) --\n";

      // Tilted normal (1,1,0)/sqrt(2) — pre-stress σ_yy, σ_xy have
      // non-zero n-projection here (σ_pre · n ≠ 0), so the base
      // FreeSurface would radiate; FreeSurfaceTotal must NOT.
      real_t nor_tilt[3] = {1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0};

      // Bulk pre-stress state (SYY=+σ_n0, SXY=-τ_ini as in TPV102 init)
      real_t Q_bg[NUM_STATE] = {0};
      Q_bg[SYY] =  1.20e8;
      Q_bg[SXY] = -7.50e7;
      real_t Q_self[NUM_STATE];
      std::memcpy(Q_self, Q_bg, NUM_STATE * sizeof(real_t));

      real_t F_total[NUM_STATE], F_interior[NUM_STATE];
      flux.FreeSurfaceTotal(nor_tilt, Q_self, Q_bg, F_total);
      flux.Interior        (nor_tilt, Q_bg,   Q_bg, F_interior);

      const real_t err   = MaxDiff(F_total, F_interior);
      const real_t scale = std::max(MaxNorm(F_total), MaxNorm(F_interior));
      std::cout << "  F_interior max|component| = " << std::scientific
                << std::setprecision(3) << MaxNorm(F_interior) << " Pa\n";
      std::cout << "  F_total    max|component| = "
                << MaxNorm(F_total) << " Pa\n";
      std::cout << "  max|diff|                 = " << err << " Pa\n";
      TEST_LE(err, 32 * kEps * scale,
              "R-I06-005: FreeSurfaceTotal = Interior on uniform Q=Q_bg "
              "(tilted normal; no spurious radiation)");
   }

   // ------------------------------------------------------------------
   // Gate 4: same for FreeSurfaceGodunovTotal.
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate 4: Q_self = Q_bg on tilted normal ⇒ "
                << "FreeSurfaceGodunovTotal = Interior --\n";

      real_t nor_tilt[3] = {1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0};
      real_t Q_bg[NUM_STATE] = {0};
      Q_bg[SYY] =  1.20e8;
      Q_bg[SXY] = -7.50e7;
      real_t Q_self[NUM_STATE];
      std::memcpy(Q_self, Q_bg, NUM_STATE * sizeof(real_t));

      real_t F_total[NUM_STATE], F_interior[NUM_STATE];
      flux.FreeSurfaceGodunovTotal(nor_tilt, Q_self, Q_bg, F_total);
      flux.Interior               (nor_tilt, Q_bg,   Q_bg, F_interior);

      const real_t err   = MaxDiff(F_total, F_interior);
      const real_t scale = std::max(MaxNorm(F_total), MaxNorm(F_interior));
      std::cout << "  max|diff|  = " << std::scientific
                << std::setprecision(3) << err << " Pa\n";
      TEST_LE(err, 32 * kEps * scale,
              "R-I06-005: FreeSurfaceGodunovTotal = Interior on uniform "
              "Q=Q_bg (tilted normal)");
   }

   // ------------------------------------------------------------------
   // Gate 5: Q_self = Q_bg on a 45-deg YZ tilted normal — orthogonal to
   //         the (1,1,0)/sqrt(2) axis above, checks frame-invariance.
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate 5: second tilted normal ⇒ FreeSurfaceTotal "
                << "= Interior --\n";

      real_t nor_tilt[3] = {0.0, 1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0)};
      real_t Q_bg[NUM_STATE] = {0};
      Q_bg[SYY] =  1.20e8;
      Q_bg[SXY] = -7.50e7;
      real_t Q_self[NUM_STATE];
      std::memcpy(Q_self, Q_bg, NUM_STATE * sizeof(real_t));

      real_t F_total[NUM_STATE], F_interior[NUM_STATE];
      flux.FreeSurfaceTotal(nor_tilt, Q_self, Q_bg, F_total);
      flux.Interior        (nor_tilt, Q_bg,   Q_bg, F_interior);

      const real_t err   = MaxDiff(F_total, F_interior);
      const real_t scale = std::max(MaxNorm(F_total), MaxNorm(F_interior));
      TEST_LE(err, 32 * kEps * scale,
              "R-I06-005: FreeSurfaceTotal on (0,1,1)/sqrt(2) normal also "
              "preserves uniform Q=Q_bg");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
