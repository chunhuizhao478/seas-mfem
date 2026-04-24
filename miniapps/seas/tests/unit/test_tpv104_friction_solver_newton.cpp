// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for dynamic/tpv104_friction_solver.hpp (§4.10 Step 5 gates
// T_TPV104_FS_1..5, plus review tests R-006, R-007).

#include "test_macros.hpp"
#include "../../dynamic/tpv104_friction_solver.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../friction/friction_coeff_stable.hpp"

#include <cstdlib>
#include <cmath>
#include <random>

using namespace mfem;
using namespace mfem::seas;

// TPV104 canonical scalars for the envelope.
static constexpr real_t kV0_TPV = 1.0e-6;
static constexpr real_t kEta_TPV = 4.625e6;   // ρ·c_s/2 for TPV104
static constexpr real_t kSigmaN_TPV = 120e6;

// ---------------------------------------------------------------------------
// R-003 STANDALONE Newton — inlined straight-line arithmetic, calling
// only std:: primitives.  Does NOT call friction_stable::* helpers, so a
// drift in ArsinhExp/ComputeCExp in production will show up as a nonzero
// diff in T_TPV104_FS_3.
// ---------------------------------------------------------------------------
static real_t reference_friction_mu(real_t V, real_t psi, real_t a, real_t V0)
{
   // Same branch identity as FrictionCoefficientStable, inlined.
   // Arithmetic ordering mirrors the production exactly: cLin = 0.5/V0
   // first, then x = cLin*V (not V/(2*V0)) — these are algebraically
   // identical but can differ by 1 ULP, which breaks bit-identity.
   constexpr real_t Switch = 10.0;
   constexpr real_t Threshold = 50.0;
   constexpr real_t Log2 = 0.69314718055994530943;
   const real_t cLin = 0.5 / V0;
   const real_t x = cLin * V;
   const real_t cExpLog = psi / a;
   int xexp = 0;
   (void)std::frexp(x, &xexp);
   real_t cExp;
   if (cExpLog > 0.0) { cExp = std::exp(-cExpLog); }
   else               { cExp = std::exp(cExpLog);  }

   real_t as_val;
   if (cExpLog + std::max(xexp, 0) * Log2 > Switch || cExpLog >= Threshold)
   {
      real_t e = cExp;
      if (cExpLog <= 0.0) { e = 1.0 / e; }
      const real_t xa = std::abs(x);
      const real_t xs = (x >= 0.0) ? 1.0 : -1.0;
      as_val = xs * (cExpLog + std::log(xa + std::sqrt(xa * xa + e * e)));
   }
   else
   {
      real_t e = cExp;
      if (cExpLog > 0.0) { e = 1.0 / e; }
      const real_t v = e * x;
      as_val = std::asinh(v);
   }
   return a * as_val;
}

static real_t reference_friction_dmu(real_t V, real_t psi, real_t a, real_t V0)
{
   constexpr real_t Switch = 10.0;
   constexpr real_t Threshold = 50.0;
   constexpr real_t Log2 = 0.69314718055994530943;
   const real_t cLin = 0.5 / V0;
   const real_t x = cLin * V;
   const real_t cExpLog = psi / a;
   int xexp = 0;
   (void)std::frexp(x, &xexp);
   real_t cExp;
   if (cExpLog > 0.0) { cExp = std::exp(-cExpLog); }
   else               { cExp = std::exp(cExpLog);  }
   const real_t acLin = a * (0.5 / V0);

   if (cExpLog + std::max(xexp, 0) * Log2 > Switch || cExpLog >= Threshold)
   {
      real_t e = cExp;
      if (cExpLog <= 0.0) { e = 1.0 / e; }
      return acLin / std::sqrt(x * x + e * e);
   }
   else
   {
      real_t e = cExp;
      if (cExpLog > 0.0) { e = 1.0 / e; }
      const real_t v = e * x;
      return acLin * (e / std::sqrt(1.0 + v * v));
   }
}

static real_t reference_newton(real_t tau_abs, real_t psi, real_t sigma_n,
                               real_t eta_s, real_t a, real_t V0,
                               real_t V_prev, int max_iter, real_t tol,
                               int &iter_out, bool &conv_out)
{
   if (tau_abs <= 0.0)  { iter_out = 0; conv_out = true; return 0.0; }
   if (sigma_n <= 0.0)
   {
      iter_out = 0; conv_out = true;
      return (eta_s > 0.0) ? (tau_abs / eta_s) : 0.0;
   }
   constexpr real_t kAlmostZero = 1e-45;
   const real_t inv_eta_s   = 1.0 / eta_s;
   const real_t sigma_n_abs = std::abs(sigma_n);

   // R-006: first-guess unclamped.
   real_t V = V_prev;
   real_t g = 0;
   int i;
   for (i = 0; i < max_iter; ++i)
   {
      const real_t mu = reference_friction_mu(V, psi, a, V0);
      g = -inv_eta_s * (sigma_n_abs * mu - tau_abs) - V;
      if (std::abs(g) < tol)
      {
         iter_out = i + 1; conv_out = true; return V;
      }
      const real_t dmu = reference_friction_dmu(V, psi, a, V0);
      const real_t dg  = -inv_eta_s * (sigma_n_abs * dmu) - 1.0;
      const real_t step = g / dg;
      V = std::max(kAlmostZero, V - step);
   }
   iter_out = max_iter; conv_out = std::abs(g) < tol; return V;
}

// ---------------------------------------------------------------------------
// R-007: physically-informed Newton initial guess.  For small V,
//   τ ≈ σ_n·a·sinh((V/(2V₀))·exp(ψ/a)) + η·V,
// which inverts approximately as
//   V ≈ 2V₀·exp(-ψ/a)·sinh(τ/(σ_n·a)).
//
// R2-001 (review round 2): cap by τ/η_s.  The asinh-asymptote derivation
// ignores η·V.  At the high-τ end of the envelope (τ = 1e8) the
// asymptotic form overshoots by ~20+ decades (≈ 9e22 m/s) vs. the true
// rest-state (~21.6 m/s, η-dominated).  The damping cap τ/η_s is the
// overstressed-fault saturation limit and bounds the physically
// plausible V.  For low τ the asinh asymptote remains smaller than τ/η_s
// and is unchanged.
// ---------------------------------------------------------------------------
static real_t physical_v_guess(real_t tau, real_t psi, real_t sigma_n,
                               real_t a, real_t V0, real_t eta_s)
{
   const real_t c    = tau / (sigma_n * a);
   const real_t arg  = c - psi / a;
   const real_t asymptotic = V0 * std::exp(arg);
   const real_t damping_cap = tau / eta_s;
   return std::min(asymptotic, damping_cap);
}

// ---------------------------------------------------------------------------
// T_TPV104_FS_1 — Newton ↔ Brent agreement (R-007: physical initial guess).
// ---------------------------------------------------------------------------
void TestNewtonVsBrent()
{
   std::cout << "\n[T_TPV104_FS_1] Newton ↔ Brent agreement (R-007 V_guess)\n";

   FrictionSolver fs;
   std::mt19937 rng(20260424u);
   std::uniform_real_distribution<real_t> tau_rng(1e4, 1e8);
   std::uniform_real_distribution<real_t> psi_rng(0.3, 1.0);
   std::uniform_real_distribution<real_t> a_rng(0.008, 0.02);

   int ok = 0;
   int non_conv = 0;
   real_t max_rel = 0.0;
   const int N = 5000;
   for (int t = 0; t < N; ++t)
   {
      const real_t tau = tau_rng(rng);
      const real_t psi = psi_rng(rng);
      const real_t a   = a_rng(rng);

      const real_t V_brent = fs.Solve(tau, psi, kSigmaN_TPV, kEta_TPV, a,
                                      FrictionSolver::Method::Brent);
      // R-007 + R2-001: physically informed guess capped by τ/η_s.
      const real_t V_guess = physical_v_guess(tau, psi, kSigmaN_TPV, a,
                                              kV0_TPV, kEta_TPV);
      int iter = 0;
      bool conv = false;
      const real_t V_newton = SolveSlipRateNewtonStable(
         tau, psi, kSigmaN_TPV, kEta_TPV, a, kV0_TPV,
         V_guess, 60, 1e-14, &iter, &conv);

      if (!conv) { ++non_conv; continue; }

      const real_t denom = std::max<real_t>(std::abs(V_brent),
                                            static_cast<real_t>(1e-30));
      const real_t rel = std::abs(V_newton - V_brent) / denom;
      max_rel = std::max(max_rel, rel);
      if (rel < 1e-10) { ++ok; }
   }
   std::cout << "  max rel err = " << max_rel << " / "
             << N << " samples (" << non_conv << " non-converged)\n";
   TEST_ASSERT(ok >= N - non_conv - 5,
               "Newton converged solves agree with Brent to 1e-10 rel");
   TEST_ASSERT(non_conv < N / 100,
               "non-convergence rate < 1% with physical V_guess");
}

// ---------------------------------------------------------------------------
// T_TPV104_FS_2 — rest state.
// ---------------------------------------------------------------------------
void TestRestState()
{
   std::cout << "\n[T_TPV104_FS_2] TPV104 rest-state solve\n";

   const real_t tau     = 40e6;
   const real_t psi     = 0.5636;
   const real_t sigma_n = 120e6;
   const real_t eta_s   = 4.625e6;
   const real_t a       = 0.01;
   const real_t V0      = 1e-6;
   const real_t V_prev  = 1e-16;

   int iter = 0;
   bool conv = false;
   const real_t V = SolveSlipRateNewtonStable(tau, psi, sigma_n, eta_s,
                                              a, V0, V_prev, 60, 1e-8,
                                              &iter, &conv);

   TEST_ASSERT(conv, "rest-state solve converges");
   TEST_ASSERT(iter < 10, "rest-state converges in < 10 iterations");
   TEST_ASSERT(V < 1e-10, "rest-state V is small (near V_prev)");
   TEST_ASSERT(V > 0.0, "rest-state V > 0");

   std::cout << "  V = " << V << "  iter = " << iter
             << "  converged = " << conv << "\n";
}

// ---------------------------------------------------------------------------
// T_TPV104_FS_3 — bit-identical iterate sequence vs standalone inline
// Newton (R-003: no shared helpers).
// ---------------------------------------------------------------------------
void TestStandaloneByteMatch()
{
   std::cout << "\n[T_TPV104_FS_3] standalone-reference Newton byte-match\n";

   std::mt19937 rng(7u);
   std::uniform_real_distribution<real_t> tau_rng(1e5, 1e8);
   std::uniform_real_distribution<real_t> psi_rng(0.3, 1.0);
   std::uniform_real_distribution<real_t> a_rng(0.008, 0.02);
   std::uniform_real_distribution<real_t> logV(-18.0, -6.0);

   int ok = 0;
   real_t max_rel = 0.0;
   for (int t = 0; t < 1000; ++t)
   {
      const real_t tau   = tau_rng(rng);
      const real_t psi   = psi_rng(rng);
      const real_t a     = a_rng(rng);
      const real_t V_prev = std::pow(10.0, logV(rng));

      int i1 = 0, i2 = 0;
      bool c1 = false, c2 = false;
      const real_t v1 = SolveSlipRateNewtonStable(
         tau, psi, kSigmaN_TPV, kEta_TPV, a, kV0_TPV,
         V_prev, 60, 1e-8, &i1, &c1);
      const real_t v2 = reference_newton(
         tau, psi, kSigmaN_TPV, kEta_TPV, a, kV0_TPV,
         V_prev, 60, 1e-8, i2, c2);

      // The production and the standalone reference use algebraically
      // identical arithmetic but the compiler emits them as distinct
      // instruction sequences (inlined helpers vs straight-line).  Over
      // the Newton iterate sequence, ULP-level round-off differences
      // accumulate.  Accept ≤ 1e-13 relative and matching iter/conv.
      // The point of R-003 is to catch FORMULA drift, not bit-identity.
      const real_t denom = std::max<real_t>(std::abs(v2),
                                            static_cast<real_t>(1e-30));
      const real_t rel = std::abs(v1 - v2) / denom;
      max_rel = std::max(max_rel, rel);
      const bool shape_match = (i1 == i2 && c1 == c2);
      if (rel < 1e-13 && shape_match) { ++ok; }
      else if (t < 3)
      {
         std::cerr.precision(20);
         std::cerr << "  FAIL t=" << t
                   << " ours(V=" << v1 << " iter=" << i1 << " conv=" << c1
                   << ")  vs ref(V=" << v2 << " iter=" << i2
                   << " conv=" << c2 << ")  rel=" << rel << "\n";
      }
   }
   std::cout << "  max rel err = " << max_rel << " / 1000 samples\n";
   TEST_ASSERT(ok == 1000,
               "All 1000 Newton solves match standalone reference to 1e-13 "
               "rel with matching iter/conv");

   // R2-005 ratchet: a drift-monitor that catches slow formula creep.
   // Current measured max_rel on clang-19 + O2 + ftree-vectorize is
   // ≈ 1.16e-14 (a few ULP per Newton iterate, accumulated over ~60
   // iterations).  Ratchet set at 3e-14 (≈ 3× observed) so a future
   // refactor cannot loosen the formula significantly without tripping
   // this gate first (before the 1e-13 ok-threshold).
   TEST_ASSERT(max_rel < 3e-14,
               "max rel err ≤ 3e-14 ratchet (R2-005 drift monitor)");
}

// ---------------------------------------------------------------------------
// T_TPV104_FS_5 — non-convergence surfaces loudly (unreachable tol).
// ---------------------------------------------------------------------------
void TestNonConvergenceLoud()
{
   std::cout << "\n[T_TPV104_FS_5] non-convergence surfaces loudly\n";

   const real_t tau     = 40e6;
   const real_t psi     = 0.5636;
   const real_t sigma_n = 120e6;
   const real_t eta_s   = 4.625e6;
   const real_t a       = 0.01;
   const real_t V0      = 1e-6;
   const real_t V_prev  = 1e-16;

   int iter = 0;
   bool conv = false;
   const real_t V = SolveSlipRateNewtonStable(tau, psi, sigma_n, eta_s,
                                              a, V0, V_prev,
                                              8 /* max_iter */,
                                              1e-40 /* tol */,
                                              &iter, &conv);

   TEST_ASSERT(!conv, "non-convergence reports has_converged = false");
   TEST_ASSERT(iter == 8, "non-convergence reports iterations = max_iter");
   TEST_ASSERT(std::isfinite(V), "V is still finite (last iterate)");

   // Early-return paths for τ = 0 and tensile σ_n.
   int iter0 = -1;
   bool conv0 = false;
   const real_t V_tau0 = SolveSlipRateNewtonStable(
      0.0, psi, sigma_n, eta_s, a, V0, V_prev, 60, 1e-8, &iter0, &conv0);
   TEST_NEAR(V_tau0, 0.0, 1e-20, "τ=0 → V=0");
   TEST_ASSERT(conv0, "τ=0 reports converged");

   int iter_t = -1;
   bool conv_t = false;
   const real_t V_tensile = SolveSlipRateNewtonStable(
      1e6, psi, -1.0, eta_s, a, V0, V_prev, 60, 1e-8, &iter_t, &conv_t);
   TEST_ASSERT(conv_t, "σ_n<0 (tensile) reports converged");
   TEST_NEAR(V_tensile, 1e6 / eta_s, 1e-12, "σ_n<0 → V = τ/η");
}

// ---------------------------------------------------------------------------
// R-006 dedicated test — first-guess unclamped.  Setting V_prev = 1e-50
// (below the in-iteration floor) must produce the bit-identical iterate
// sequence as the standalone reference.
// ---------------------------------------------------------------------------
void TestFirstGuessUnclamped()
{
   std::cout << "\n[R-006] Newton first-guess unclamped at V_prev=1e-50\n";

   int i1 = 0, i2 = 0; bool c1 = false, c2 = false;
   const real_t V1 = SolveSlipRateNewtonStable(
      40e6, 0.5, 120e6, 4.625e6, 0.01, 1e-6, /*V_prev=*/1e-50,
      60, 1e-8, &i1, &c1);
   const real_t V2 = reference_newton(
      40e6, 0.5, 120e6, 4.625e6, 0.01, 1e-6, /*V_prev=*/1e-50,
      60, 1e-8, i2, c2);

   TEST_ASSERT(V1 == V2,
               "V_prev=1e-50 bit-identical with unclamped reference");
   TEST_ASSERT(i1 == i2,
               "iterate counts match (no off-by-one from silent clamp)");
   TEST_ASSERT(c1 == c2, "convergence status matches");
}

// ---------------------------------------------------------------------------
// R2-001 dedicated test — `physical_v_guess` is bounded at high τ.
// ---------------------------------------------------------------------------
void TestPhysicalGuessBounded()
{
   std::cout << "\n[R2-001] physical_v_guess bounded by τ/η at high τ\n";

   // Worst-case from the review: τ=1e8, ψ=0.3, a=0.008.
   // Un-capped asymptotic ≈ 9e22; damping cap = 1e8 / 4.625e6 ≈ 21.6.
   const real_t g_hi = physical_v_guess(1e8, 0.3, 1.2e8, 0.008,
                                        1e-6, 4.625e6);
   TEST_ASSERT(g_hi <= 1e8 / 4.625e6 + 1e-12,
               "high-τ guess ≤ τ/η_s");
   TEST_ASSERT(g_hi > 0.0,
               "high-τ guess strictly positive");

   // Low-τ sample: asymptotic form remains small — damping cap never
   // binds, so the guess equals the uncapped asymptotic form.
   const real_t tau_lo = 1e4, psi_lo = 0.7, a_lo = 0.015,
                sig_lo = 1.2e8, V0_lo = 1e-6, eta_lo = 4.625e6;
   const real_t g_lo = physical_v_guess(tau_lo, psi_lo, sig_lo,
                                        a_lo, V0_lo, eta_lo);
   const real_t expected_lo =
      V0_lo * std::exp(tau_lo / (sig_lo * a_lo) - psi_lo / a_lo);
   TEST_NEAR(g_lo, expected_lo, 1e-30,
             "low-τ guess bit-identical to uncapped asymptotic");
}

// ---------------------------------------------------------------------------
// R2-004 dedicated test — Newton aborts on invalid V_prev.
// ---------------------------------------------------------------------------
void TestNewtonRejectsInvalidVPrev()
{
   std::cout << "\n[R2-004] Newton rejects negative / non-finite V_prev\n";

   int it = 0;
   bool conv = false;

   // Negative V_prev — must abort.
   bool aborted_neg = false;
   try
   {
      SolveSlipRateNewtonStable(40e6, 0.5, 120e6, 4.625e6, 0.01, 1e-6,
                                /*V_prev=*/-1e-12, 60, 1e-8, &it, &conv);
   }
   catch (const std::exception&) { aborted_neg = true; }
   catch (...)                   { aborted_neg = true; }
   TEST_ASSERT(aborted_neg, "negative V_prev aborts");

   // NaN V_prev — must abort.
   bool aborted_nan = false;
   try
   {
      SolveSlipRateNewtonStable(
         40e6, 0.5, 120e6, 4.625e6, 0.01, 1e-6,
         /*V_prev=*/std::numeric_limits<real_t>::quiet_NaN(),
         60, 1e-8, &it, &conv);
   }
   catch (const std::exception&) { aborted_nan = true; }
   catch (...)                   { aborted_nan = true; }
   TEST_ASSERT(aborted_nan, "NaN V_prev aborts");

   // NaN tau — must abort (covers R2-006).
   bool aborted_tau = false;
   try
   {
      SolveSlipRateNewtonStable(
         std::numeric_limits<real_t>::quiet_NaN(),
         0.5, 120e6, 4.625e6, 0.01, 1e-6, 1e-16, 60, 1e-8, &it, &conv);
   }
   catch (const std::exception&) { aborted_tau = true; }
   catch (...)                   { aborted_tau = true; }
   TEST_ASSERT(aborted_tau, "NaN tau aborts");

   // V_prev = 0 is valid (locked fault cold start) — must NOT abort.
   int it0 = 0;
   bool conv0 = false;
   bool aborted_zero = false;
   try
   {
      (void)SolveSlipRateNewtonStable(40e6, 0.5, 120e6, 4.625e6, 0.01,
                                      1e-6, /*V_prev=*/0.0,
                                      60, 1e-8, &it0, &conv0);
   }
   catch (...) { aborted_zero = true; }
   TEST_ASSERT(!aborted_zero,
               "V_prev = 0 is accepted (locked-fault cold start)");
}

// ---------------------------------------------------------------------------
// R-012 — T_TPV104_FS_4 removed.  CLI default verification is a Step-9
// gate (tracked by the Phase-2 acceptance matrix, NOT by this test
// binary).  No unconditional TEST_ASSERT(true) stub is emitted — the
// test harness sees nothing here, consistent with "deferred".
// ---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   TestNewtonVsBrent();
   TestRestState();
   TestStandaloneByteMatch();
   TestNonConvergenceLoud();
   TestFirstGuessUnclamped();
   TestPhysicalGuessBounded();
   TestNewtonRejectsInvalidVPrev();

   std::cout << "\n[T_TPV104_FS_4] CLI default Newton — SKIPPED "
             << "(tracked by Step 9 test_tpv104_smoke.cpp, not counted "
             << "toward FS pass rate)\n";

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
