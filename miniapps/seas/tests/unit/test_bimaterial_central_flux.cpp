// Phase 1 (PLAN_mixed_flux_hetero_riemann.md): BimaterialFlux::
// BuildPerFaceCentralMatricesGlobal primitive correctness tests.
//
// Verifies the bi-material CENTRAL-flux builder produces ½·A_self / ½·A_nbr
// (global frame) such that, via the existing ApplyPerFaceFlux,
//   F* = ½·A_self·Q_self + ½·A_nbr·Q_nbr.
//
//   Test 1.1 — homogeneous limit equals GodunovFlux::Central to <= 1e-11
//              relative (>=10 materials x >=50 (Q_self,Q_nbr,nor)).
//   Test 1.2 — consistency: Q_self==Q_nbr==Q, one material => F* == A·Q,
//              compared to a GetAx-based global face-normal Jacobian × Q.
//   Test 1.3 — heterogeneous averaging (BUG-3 primitive anchor): centralSelf
//              uses A_self and centralNbr uses the NEIGHBOUR's A_nbr (NOT
//              A_self for both); plus a discriminating-power check that the
//              two sides' contributions genuinely differ across the contrast.
//
// Note (PLAN BUG-7): the single-valuedness invariant F_h_e1==F_h_e2 is an
// OPERATOR-level property tested in Phase 3 (Test 3.2), NOT here — so this
// file deliberately contains no "swap (self,nbr) and states" primitive test
// (that would only check commutativity of addition and would pass with the
// BUG-3 dispatch defect present).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/godunov_flux_bimaterial.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

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

#define TEST_GE(v, thr, msg) do { \
   num_tests++; \
   double vv = (v), tt = (thr); \
   if (vv >= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", thr " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected >= " << tt << ")\n"; } \
} while (0)

namespace
{
struct Mat { real_t lam, mu, rho; };

// out = M · x  (M is NUM_STATE × NUM_STATE).
void MatVec(const DenseMatrix &M, const real_t *x, real_t *out)
{
   for (int i = 0; i < NUM_STATE; ++i)
   {
      real_t s = 0.0;
      for (int j = 0; j < NUM_STATE; ++j) { s += M(i, j) * x[j]; }
      out[i] = s;
   }
}

real_t OutScale(const real_t *a, const real_t *b)
{
   real_t s = 0.0;
   for (int i = 0; i < NUM_STATE; ++i)
   {
      s = std::max(s, std::abs(a[i]));
      s = std::max(s, std::abs(b[i]));
   }
   return s;
}

// Worst per-component error relative to the output magnitude (scale-invariant;
// matches the relative-tolerance convention of test_godunov_central_flux).
real_t WorstRel(const real_t *a, const real_t *b)
{
   const real_t scale = OutScale(a, b);
   if (scale == 0.0) { return 0.0; }
   real_t m = 0.0;
   for (int i = 0; i < NUM_STATE; ++i)
   {
      m = std::max(m, std::abs(a[i] - b[i]) / scale);
   }
   return m;
}

void NormalizeVec(real_t *n)
{
   const real_t len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   if (len > 0.0) { n[0] /= len; n[1] /= len; n[2] /= len; }
}

// (IMPL-6) Portable uniform real in [lo, hi] from mt19937.  std::uniform_real_
// distribution is NOT specified to produce the same sequence across stdlib
// implementations (libstdc++ / libc++ / MSVC), so we map mt19937's 32-bit
// output directly — fully portable AND deterministic.  (This is a tolerance
// test asserting `worst <= tol`, not a golden-value test, so only
// reproducibility — not statistical optimality — is required of the map.)
real_t URand(std::mt19937 &rng, real_t lo, real_t hi)
{
   const real_t u =
      static_cast<real_t>(rng() - std::mt19937::min())
      / static_cast<real_t>(std::mt19937::max() - std::mt19937::min());
   return lo + (hi - lo) * u;
}

// Random physical isotropic material: rho, cs, cp = (cp/cs)·cs, derive
// (lambda, mu).  Ranges are solid-earth SI (mu ~ 4e9..7e10 >> mu_eps).
Mat RandomMaterialParams(std::mt19937 &rng)
{
   Mat m;
   m.rho = URand(rng, 2000.0, 3500.0);
   const real_t cs = URand(rng, 1500.0, 4500.0);
   const real_t cp = URand(rng, 1.55, 2.05) * cs;   // cp/cs ratio
   m.mu  = m.rho * cs * cs;
   m.lam = m.rho * cp * cp - 2.0 * m.mu;
   return m;
}

// Random state, physically-scaled per channel (normal stresses ~1e7, shear
// stresses ~1e6, velocities ~1) so all channels are exercised.
void RandomState(std::mt19937 &rng, real_t *Q)
{
   for (int c = 0; c < NUM_STATE; ++c)
   {
      const real_t scale = (c < SXY) ? 1.0e7 : (c < VX) ? 1.0e6 : 1.0;
      Q[c] = scale * URand(rng, -1.0, 1.0);
   }
}

void RandomUnitNormal(std::mt19937 &rng, real_t *n)
{
   do {
      n[0] = URand(rng, -1.0, 1.0);
      n[1] = URand(rng, -1.0, 1.0);
      n[2] = URand(rng, -1.0, 1.0);
   } while (n[0]*n[0] + n[1]*n[1] + n[2]*n[2] < 1.0e-4);
   NormalizeVec(n);
}

// Global-frame face-normal Jacobian from GetAx(): A_g = T · GetAx() · T^{-1}.
// Independent reference for Test 1.2 (the analytic physical-flux Jacobian).
void BuildGlobalJacobianFromGetAx(const GodunovFlux &flux, const real_t *nor,
                                  DenseMatrix &A_g)
{
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(nor, t1, t2);
   DenseMatrix T(NUM_STATE, NUM_STATE), Tinv(NUM_STATE, NUM_STATE),
               tmp(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotation(nor, t1, t2, T);
   GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);
   mfem::Mult(T, flux.GetAx(), tmp);
   A_g.SetSize(NUM_STATE, NUM_STATE);
   mfem::Mult(tmp, Tinv, A_g);
}

} // anonymous namespace

int main()
{
   std::cout << "\n=== Phase 1: BimaterialFlux::BuildPerFaceCentralMatricesGlobal ===\n";

   std::mt19937 rng(20260604u);     // fixed seed -> deterministic
   const int kNumMaterials = 12;    // >= 10 (PLAN Test 1.1)
   const int kNumStates    = 60;    // >= 50 (PLAN Test 1.1)

   // -----------------------------------------------------------------
   // Test 1.1: homogeneous limit equals GodunovFlux::Central (<= 1e-11 rel).
   // Both sides use (Ax_plus_ + Ax_minus_), so the only difference from
   // Central is FP operation order (matrix-matrix-then-matvec here vs
   // rotate-apply-rotate in Central) -> ~1e-13, comfortably under 1e-11.
   // -----------------------------------------------------------------
   std::cout << "\n-- Test 1.1: homogeneous limit == GodunovFlux::Central --\n";
   real_t worst_11 = 0.0;
   for (int m = 0; m < kNumMaterials; ++m)
   {
      const Mat p = RandomMaterialParams(rng);
      GodunovFlux g(p.lam, p.mu, p.rho);
      for (int s = 0; s < kNumStates; ++s)
      {
         real_t Qs[NUM_STATE], Qn[NUM_STATE], nor[3];
         RandomState(rng, Qs);
         RandomState(rng, Qn);
         RandomUnitNormal(rng, nor);

         DenseMatrix cS, cN;
         BimaterialFlux::BuildPerFaceCentralMatricesGlobal(nor, g, g, cS, cN);

         real_t Fc[NUM_STATE], Fg[NUM_STATE];
         BimaterialFlux::ApplyPerFaceFlux(cS, cN, Qs, Qn, Fc);
         g.Central(nor, Qs, Qn, Fg);

         worst_11 = std::max(worst_11, WorstRel(Fc, Fg));
      }
   }
   TEST_LE(worst_11, 1.0e-11,
           "Test 1.1: ApplyPerFaceFlux(half-A,half-A,Qs,Qn) == "
           "GodunovFlux::Central in the homogeneous limit "
           "(12 materials x 60 states)");

   // -----------------------------------------------------------------
   // Test 1.2: consistency Q_self==Q_nbr==Q => F* == A·Q (physical flux).
   // With one material, F* = half-A·Q + half-A·Q = A·Q.  Reference A_g is
   // the GetAx-based global Jacobian.
   // -----------------------------------------------------------------
   std::cout << "\n-- Test 1.2: consistency F* == A·Q (physical flux) --\n";
   real_t worst_12 = 0.0;
   for (int m = 0; m < kNumMaterials; ++m)
   {
      const Mat p = RandomMaterialParams(rng);
      GodunovFlux g(p.lam, p.mu, p.rho);
      for (int s = 0; s < kNumStates; ++s)
      {
         real_t Q[NUM_STATE], nor[3];
         RandomState(rng, Q);
         RandomUnitNormal(rng, nor);

         DenseMatrix cS, cN;
         BimaterialFlux::BuildPerFaceCentralMatricesGlobal(nor, g, g, cS, cN);
         real_t Fc[NUM_STATE];
         BimaterialFlux::ApplyPerFaceFlux(cS, cN, Q, Q, Fc);  // = A·Q (homog)

         DenseMatrix A_g;
         BuildGlobalJacobianFromGetAx(g, nor, A_g);
         real_t Fref[NUM_STATE];
         MatVec(A_g, Q, Fref);

         worst_12 = std::max(worst_12, WorstRel(Fc, Fref));
      }
   }
   TEST_LE(worst_12, 1.0e-11,
           "Test 1.2: F*(Q,Q) == A·Q (GetAx-based global Jacobian); one "
           "material, <= 1e-11 relative");

   // -----------------------------------------------------------------
   // Test 1.3 (BUG-3 primitive anchor): heterogeneous averaging.  The
   // builder must use A_self for the self term and the NEIGHBOUR's A_nbr
   // for the neighbour term, NOT A_self for both.  Reference half-A_e1 /
   // half-A_e2 come from the SAME builder in its homogeneous call, so the
   // comparison isolates the (self,nbr) WIRING from the matrix math.
   //
   // The neighbour material is forced to contrast with the self material
   // (cs2 = factor·cs1, factor in [1.3, 2.0]) so half-A_e1·Q and
   // half-A_e2·Q genuinely differ — otherwise the "uses A_nbr" check could
   // pass even if centralNbr wrongly used A_self.
   // -----------------------------------------------------------------
   std::cout << "\n-- Test 1.3: heterogeneous averaging (uses NEIGHBOUR A_nbr) --\n";
   real_t worst_13_self = 0.0, worst_13_nbr = 0.0, worst_13_anchor = 0.0,
          min_contrast = 1.0e30;
   const real_t zero[NUM_STATE] = {0};
   for (int m = 0; m < kNumMaterials; ++m)
   {
      const Mat p1 = RandomMaterialParams(rng);
      GodunovFlux g1(p1.lam, p1.mu, p1.rho);

      // g2: forced material contrast vs g1.
      const real_t f   = URand(rng, 1.3, 2.0);
      const real_t cs1 = std::sqrt(p1.mu / p1.rho);
      const real_t cs2 = f * cs1;
      const real_t cp2 = 1.8 * cs2;
      const real_t rho2 = p1.rho;
      const real_t mu2  = rho2 * cs2 * cs2;
      const real_t lam2 = rho2 * cp2 * cp2 - 2.0 * mu2;
      GodunovFlux g2(lam2, mu2, rho2);

      for (int s = 0; s < kNumStates; ++s)
      {
         real_t Q[NUM_STATE], nor[3];
         RandomState(rng, Q);
         RandomUnitNormal(rng, nor);

         DenseMatrix cS12, cN12;   // heterogeneous build (self=g1, nbr=g2)
         BimaterialFlux::BuildPerFaceCentralMatricesGlobal(nor, g1, g2, cS12, cN12);

         // Reference half-A_e1 / half-A_e2 from the same builder (homog call).
         DenseMatrix halfA1, halfA2, dummy;
         BimaterialFlux::BuildPerFaceCentralMatricesGlobal(nor, g1, g1, halfA1, dummy);
         BimaterialFlux::BuildPerFaceCentralMatricesGlobal(nor, g2, g2, halfA2, dummy);

         // Self term: ApplyPerFaceFlux(cS12,cN12,Q,0) = cS12·Q ; expect half-A_e1·Q.
         real_t F_self[NUM_STATE], ref_self[NUM_STATE];
         BimaterialFlux::ApplyPerFaceFlux(cS12, cN12, Q, zero, F_self);
         MatVec(halfA1, Q, ref_self);
         worst_13_self = std::max(worst_13_self, WorstRel(F_self, ref_self));

         // Neighbour term: ApplyPerFaceFlux(cS12,cN12,0,Q) = cN12·Q ; expect
         // half-A_e2·Q (the NEIGHBOUR Jacobian, not half-A_e1·Q).
         real_t F_nbr[NUM_STATE], ref_nbr[NUM_STATE];
         BimaterialFlux::ApplyPerFaceFlux(cS12, cN12, zero, Q, F_nbr);
         MatVec(halfA2, Q, ref_nbr);
         worst_13_nbr = std::max(worst_13_nbr, WorstRel(F_nbr, ref_nbr));

         // Discriminating power: half-A_e1·Q and half-A_e2·Q must differ.
         real_t contrast = WorstRel(ref_self, ref_nbr);
         min_contrast = std::min(min_contrast, contrast);

         // (IMPL-2) INDEPENDENT het anchor: build A_g1, A_g2 from GetAx (a
         // different code path than BuildPerFaceCentralMatricesGlobal) and
         // assert ApplyPerFaceFlux(cS12,cN12,Q,Q) == 0.5*(A_g1 + A_g2)*Q.  This
         // anchors the heterogeneous OUTPUT against GetAx, NOT against the same
         // builder, so a rotation/transpose bug common to all builder calls
         // (which the self/nbr refs above would cancel) is caught here.
         DenseMatrix A_g1, A_g2;
         BuildGlobalJacobianFromGetAx(g1, nor, A_g1);
         BuildGlobalJacobianFromGetAx(g2, nor, A_g2);
         real_t F_het[NUM_STATE], a1[NUM_STATE], a2[NUM_STATE], ref_het[NUM_STATE];
         BimaterialFlux::ApplyPerFaceFlux(cS12, cN12, Q, Q, F_het);
         MatVec(A_g1, Q, a1);
         MatVec(A_g2, Q, a2);
         for (int c = 0; c < NUM_STATE; ++c) { ref_het[c] = 0.5 * (a1[c] + a2[c]); }
         worst_13_anchor = std::max(worst_13_anchor, WorstRel(F_het, ref_het));
      }
   }
   TEST_LE(worst_13_self, 1.0e-11,
           "Test 1.3: centralSelf·Q == half-A_self·Q (self term uses A_self)");
   TEST_LE(worst_13_nbr, 1.0e-11,
           "Test 1.3: centralNbr·Q == half-A_nbr·Q (neighbour term uses the "
           "NEIGHBOUR Jacobian, NOT A_self -- BUG-3 anchor)");
   TEST_LE(worst_13_anchor, 1.0e-11,
           "Test 1.3 (IMPL-2): het F*(Q,Q) == 0.5*(A_g1+A_g2)*Q vs an INDEPENDENT "
           "GetAx-based Jacobian (not the same builder)");
   TEST_GE(min_contrast, 1.0e-2,
           "Test 1.3 discriminates: half-A_self·Q and half-A_nbr·Q differ "
           "materially across the contrast (so the A_nbr check is non-trivial)");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
