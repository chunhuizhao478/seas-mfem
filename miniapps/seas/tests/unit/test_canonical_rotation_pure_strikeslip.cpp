// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.2 (regression gate): canonical rotation of a
// pure strike-slip global state Q_g onto the BP5 canonical face frame.
//
// Builds a global Q_g with a pure mode-II strike-slip Riemann state
// (only Q_g[VX]=V_ini and Q_g[SXY]=tau_ini nonzero), applies the
// BP5-canonical `Tinv` from GodunovFlux::BuildRotationInverse, and
// verifies that:
//
//   (a) Q_c[VX] (normal-direction velocity) stays 0 — no normal flux
//       from a pure strike-slip jump.
//   (b) Q_c[VY] (dip-direction velocity, t1) stays 0.
//   (c) Q_c[VZ] (strike-direction velocity, t2) picks up V_ini.
//   (d) Q_c[SXX] (normal-normal stress, Q_c[SXX]=σ_{nn}) stays 0.
//   (e) Q_c[SXY] (n·t1 shear) stays 0.
//   (f) Q_c[SXZ] (n·t2 shear, = dip·strike shear) picks up τ_ini.
//
// The BP5 canonical Q is a signed permutation (entries ∈ {−1, 0, +1}),
// so every intermediate product in Tinv is integer arithmetic and the
// assertions are **bit-exact** (tolerance = 0.0).  REVIEW.md R-002
// proved this formally; this test encodes the proof as a regression
// gate.  It subsumes the retracted rev-1 H-V92-O hypothesis.
//
// Note on QIndex convention (wave_state.hpp): Voigt order is
//   SXX=0, SYY=1, SZZ=2, SXY=3, SYZ=4, SXZ=5, VX=6, VY=7, VZ=8.
// In the face-local rotated frame Q_c, the rotation maps:
//   - Global axis (x, y, z) → face-local (n, t1, t2).
// For the BP5 canonical frame (n=-ŷ, t1=-ẑ, t2=x̂), the rotation Q has
// rows (n, t1, t2) = ((0,-1,0), (0,0,-1), (1,0,0)).

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_EQ_EXACT(val, exp, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp); \
   if (v_ == e_) { num_passed++; std::cout << "  PASSED: " << msg \
                   << " (got " << v_ << ")\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " \
          << msg << " (got " << std::setprecision(17) << v_ \
          << ", expected " << e_ << ")\n"; } \
} while (0)

int main()
{
   std::cout << "\n=== TPV102 v9.2.0 §4.2 Regression: "
             << "Canonical Rotation on Pure Strike-Slip ===\n";

   // TPV102 pre-stress scales (documentation only — doesn't need to match
   // exact simulation values; bit-exactness is independent of magnitude).
   const real_t V_ini   = 1.0e-12;   // [m/s]
   const real_t tau_ini = 75.0e6;    // [Pa]

   // -------- BP5 canonical frame (signed permutation) -----------------------
   real_t n [3] = {0.0, -1.0,  0.0};
   real_t t1[3] = {0.0,  0.0, -1.0};
   real_t t2[3] = {1.0,  0.0,  0.0};

   DenseMatrix Tinv;
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);

   // -------- Pure strike-slip global state Q_g ------------------------------
   // Only Q_g[VX] and Q_g[SXY] nonzero.  SXY is σ_xy in global.
   // This represents the Pelties-9 pre-stress state at the start of
   // TPV102 (strike-slip load along x, fault plane at y=0).
   Vector Q_g(NUM_STATE);
   Q_g = 0.0;
   Q_g(VX)  = V_ini;
   Q_g(SXY) = tau_ini;

   // -------- Apply Tinv -----------------------------------------------------
   Vector Q_c(NUM_STATE);
   Tinv.Mult(Q_g, Q_c);

   // -------- Dump for transparency ------------------------------------------
   const char *labels[NUM_STATE] = {
      "SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ", "VX", "VY", "VZ"
   };
   std::cout << "  Input  Q_g: pure strike-slip (VX=" << V_ini
             << ", SXY=" << tau_ini << ")\n";
   std::cout << "  Output Q_c (face-local, BP5 canonical frame):\n";
   for (int i = 0; i < NUM_STATE; i++)
   {
      std::cout << "    Q_c[" << std::setw(3) << std::left << labels[i]
                << "] = " << std::setw(14) << std::setprecision(6)
                << std::scientific << Q_c(i) << "\n";
   }

   // -------- Assertions: bit-exact expectations -----------------------------
   // Decomposition (BP5 canonical frame):
   //   Q_c[VX] = n · v_global    = (0,-1,0) · (V_ini,0,0) = 0
   //   Q_c[VY] = t1 · v_global   = (0,0,-1) · (V_ini,0,0) = 0
   //   Q_c[VZ] = t2 · v_global   = (1,0,0)  · (V_ini,0,0) = V_ini
   //   Q_c[SXX]= σ_{nn}          = n^T · σ · n
   //     σ_global = [[0, tau, 0], [tau, 0, 0], [0, 0, 0]]
   //     σ · n = σ · (0,-1,0)^T = (-tau, 0, 0)
   //     n^T σ n = (0,-1,0) · (-tau, 0, 0) = 0
   //   Q_c[SYY]= σ_{t1 t1}       = t1^T σ t1
   //     σ · t1 = σ · (0,0,-1)^T = (0, 0, 0)
   //     t1^T σ t1 = 0
   //   Q_c[SZZ]= σ_{t2 t2}       = t2^T σ t2
   //     σ · t2 = σ · (1,0,0)^T = (0, tau, 0)
   //     t2^T σ t2 = (1,0,0) · (0,tau,0) = 0
   //   Q_c[SXY]= σ_{n t1}        = n^T σ t1
   //     σ · t1 = (0,0,0), so n^T σ t1 = 0
   //   Q_c[SYZ]= σ_{t1 t2}       = t1^T σ t2
   //     σ · t2 = (0, tau, 0), t1^T σ t2 = (0,0,-1)·(0,tau,0) = 0
   //   Q_c[SXZ]= σ_{n t2}        = n^T σ t2
   //     σ · t2 = (0, tau, 0), n^T σ t2 = (0,-1,0)·(0,tau,0) = -tau
   TEST_EQ_EXACT(Q_c(VX),   0.0,              "Q_c[VX] = 0     (no normal vel)");
   TEST_EQ_EXACT(Q_c(VY),   0.0,              "Q_c[VY] = 0     (no dip vel)");
   TEST_EQ_EXACT(Q_c(VZ),   V_ini,            "Q_c[VZ] = V_ini (strike vel)");

   TEST_EQ_EXACT(Q_c(SXX),  0.0,              "Q_c[SXX] = 0    (no normal stress)");
   TEST_EQ_EXACT(Q_c(SYY),  0.0,              "Q_c[SYY] = 0    (no dip-dip stress)");
   TEST_EQ_EXACT(Q_c(SZZ),  0.0,              "Q_c[SZZ] = 0    (no strike-strike)");
   TEST_EQ_EXACT(Q_c(SXY),  0.0,              "Q_c[SXY] = 0    (no n·dip shear)");
   TEST_EQ_EXACT(Q_c(SYZ),  0.0,              "Q_c[SYZ] = 0    (no dip·strike shear)");
   TEST_EQ_EXACT(Q_c(SXZ), -tau_ini,          "Q_c[SXZ] = -tau (n·strike shear)");

   // -------- Diagnostic: no leakage into normal / dip channels --------------
   // H-V92-O rev-1 hypothesized "Tinv_can leakage" of ULP-scale values into
   // Q_c[SXX] or Q_c[SXY].  With signed-permutation Q the arithmetic is
   // exact: any nonzero leakage is a bug.  Above assertions already cover
   // this; restated here for documentation clarity.
   double leakage_normal = std::abs(Q_c(SXX)) + std::abs(Q_c(VX));
   double leakage_dip    = std::abs(Q_c(SXY)) + std::abs(Q_c(VY));
   TEST_EQ_EXACT(leakage_normal, 0.0, "No leakage into normal channel");
   TEST_EQ_EXACT(leakage_dip,    0.0, "No leakage into dip channel");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";

   return (num_failed == 0) ? 0 : 1;
}
