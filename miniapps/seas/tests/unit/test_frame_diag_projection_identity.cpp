// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// PLAN_frame_orthonormality_diag_2026-05-23 §Phase 1, AC #3.
//
// Guards the arithmetic identity behind the [FRAME] diagnostic added to
// WaveOperator::EvaluateBulkAtFaultQPsCanonical (wave_operator.inl).  That
// trace reports
//      dv_n = can_n · (vel_minus_global - vel_plus_global)
// and claims it equals the iterator's normal-velocity jump
//      sn_vjump / eta_p = (Q~minus - Q~plus)[VX]
// where Q~ = Tinv_can · Q_global is the canonical-frame state the substep
// iterator consumes (VX = normal component).  This test proves the claim:
//
//   (Tinv_can · Q_global)[VX] == can_n · vel(Q_global)            (identity I)
//
// for the BP5 canonical frame (exact) and a non-axis-aligned frame (<= a few
// ULP), so the [FRAME] dv_n is a faithful cross-check of [SLIP] sn_vjump.  If
// BuildRotationInverse's velocity row is ever regressed (wrong row, stress/
// velocity coupling), this fails first.
//
// Also sanity-checks the orthonormality metric the trace prints
// (n.t1, n.t2, t1.t2 -> 0 ; |.|-1 -> 0) on an orthonormal frame.

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

static constexpr double kEps = 2.220446049250313e-16;

// ---------------------------------------------------------------------------
// Reproduce the [FRAME] trace formula exactly: dv_n = can_n . (vel_m - vel_p),
// with vel(Q) = (Q[VX], Q[VY], Q[VZ]) the GLOBAL velocity sub-vector.
// ---------------------------------------------------------------------------
static real_t FrameTraceDvN(const real_t can_n[3],
                            const real_t *Q_plus_g, const real_t *Q_minus_g)
{
   const real_t dvg[3] = { Q_minus_g[VX] - Q_plus_g[VX],
                           Q_minus_g[VY] - Q_plus_g[VY],
                           Q_minus_g[VZ] - Q_plus_g[VZ] };
   return can_n[0]*dvg[0] + can_n[1]*dvg[1] + can_n[2]*dvg[2];
}

// ---------------------------------------------------------------------------
// Reproduce the iterator's canonical normal-velocity jump:
//   (Tinv_can . Q_minus_g)[VX] - (Tinv_can . Q_plus_g)[VX].
// This is sn_vjump / eta_p (the [SLIP] decomposition basis).
// ---------------------------------------------------------------------------
static real_t IteratorCanonicalVnJump(const real_t can_n[3],
                                      const real_t can_t1[3],
                                      const real_t can_t2[3],
                                      const real_t *Q_plus_g,
                                      const real_t *Q_minus_g)
{
   DenseMatrix Tinv;
   GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv);
   real_t vn_p = 0.0, vn_m = 0.0;
   for (int k = 0; k < NUM_STATE; k++)
   {
      vn_p += Tinv(VX, k) * Q_plus_g[k];
      vn_m += Tinv(VX, k) * Q_minus_g[k];
   }
   return vn_m - vn_p;
}

static real_t Dot3(const real_t a[3], const real_t b[3])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

// ---------------------------------------------------------------------------
// Identity (I) for one frame + one (Q_plus, Q_minus) global pair.
// ---------------------------------------------------------------------------
static void CheckIdentity(const std::string &label,
                          const real_t can_n[3], const real_t can_t1[3],
                          const real_t can_t2[3],
                          const real_t *Q_plus_g, const real_t *Q_minus_g,
                          double tol)
{
   const real_t dv_n     = FrameTraceDvN(can_n, Q_plus_g, Q_minus_g);
   const real_t vn_jump  = IteratorCanonicalVnJump(can_n, can_t1, can_t2,
                                                   Q_plus_g, Q_minus_g);
   const double scale = std::max<double>(1.0, std::abs(vn_jump));
   const double rel = std::abs(dv_n - vn_jump) / scale;
   std::cout << "    [" << label << "] dv_n=" << dv_n
             << " vn_jump=" << vn_jump << " rel=" << rel << "\n";
   TEST_ASSERT(rel <= tol,
               label + ": [FRAME] dv_n == iterator (Q~m-Q~p)[VX] (rel "
               + std::to_string(rel) + " <= " + std::to_string(tol) + ")");
}

int main()
{
   std::cout << "=== test_frame_diag_projection_identity ===\n";

   // A representative global 9-state on each side: distinct stresses AND
   // velocities so a spurious stress->velocity-row coupling in Tinv would
   // perturb the VX row and break identity (I).
   real_t Qp[NUM_STATE] = { 1.0, -2.0, 3.0, 0.5, -0.7, 0.9, 4.0, -5.0, 6.0 };
   real_t Qm[NUM_STATE] = { -1.5, 2.5, -3.5, 0.2, 0.3, -0.4, -2.0, 7.0, -8.0 };

   // ---- Frame 1: BP5 canonical (signed permutation; exact) ----------------
   {
      const real_t n[3]  = { 0.0, -1.0, 0.0 };
      const real_t t1[3] = { 0.0,  0.0, -1.0 };
      const real_t t2[3] = { 1.0,  0.0,  0.0 };
      // orthonormality sanity (the metric the trace prints)
      TEST_ASSERT(std::abs(Dot3(n,t1)) <= 4*kEps &&
                  std::abs(Dot3(n,t2)) <= 4*kEps &&
                  std::abs(Dot3(t1,t2)) <= 4*kEps,
                  "BP5: frame orthogonal (input sanity — hand-built frame, "
                  "not production code; the real guard is CheckIdentity)");
      TEST_ASSERT(std::abs(std::sqrt(Dot3(n,n))-1.0)  <= 4*kEps &&
                  std::abs(std::sqrt(Dot3(t1,t1))-1.0) <= 4*kEps &&
                  std::abs(std::sqrt(Dot3(t2,t2))-1.0) <= 4*kEps,
                  "BP5: frame unit-norm (input sanity — hand-built frame)");
      CheckIdentity("BP5-canonical", n, t1, t2, Qp, Qm, /*tol*/ 16*kEps);
   }

   // ---- Frame 2: a non-axis-aligned orthonormal frame ----------------------
   // n along (1,2,2)/3; t1 = normalize(up x n) analog; t2 = n x t1.
   {
      real_t n[3] = { 1.0/3.0, 2.0/3.0, 2.0/3.0 };
      // pick an arbitrary vector not parallel to n, Gram-Schmidt a tangent
      real_t a[3] = { 0.0, 0.0, 1.0 };
      const real_t an = Dot3(a, n);
      real_t t1[3] = { a[0]-an*n[0], a[1]-an*n[1], a[2]-an*n[2] };
      const real_t t1l = std::sqrt(Dot3(t1,t1));
      t1[0]/=t1l; t1[1]/=t1l; t1[2]/=t1l;
      // t2 = n x t1
      real_t t2[3] = { n[1]*t1[2]-n[2]*t1[1],
                       n[2]*t1[0]-n[0]*t1[2],
                       n[0]*t1[1]-n[1]*t1[0] };
      TEST_ASSERT(std::abs(Dot3(n,t1)) <= 8*kEps &&
                  std::abs(Dot3(n,t2)) <= 8*kEps &&
                  std::abs(Dot3(t1,t2)) <= 8*kEps,
                  "tilted: frame orthogonal (input sanity — hand-built frame)");
      CheckIdentity("tilted-frame", n, t1, t2, Qp, Qm, /*tol*/ 64*kEps);
   }

   // ---- Frame 3: sign_flipped frame (negate all three) ---------------------
   // The diagnostic uses can_* AFTER the sign_flipped negation; identity (I)
   // must hold for the negated frame too (dv_n flips sign with can_n, and so
   // does the iterator jump, so they still agree).
   {
      const real_t n[3]  = { 0.0,  1.0, 0.0 };
      const real_t t1[3] = { 0.0,  0.0, 1.0 };
      const real_t t2[3] = { -1.0, 0.0, 0.0 };
      CheckIdentity("sign-flipped", n, t1, t2, Qp, Qm, /*tol*/ 16*kEps);
   }

   std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
             << " passed, " << num_failed << " failed ===\n";
   return (num_failed == 0) ? 0 : 1;
}
