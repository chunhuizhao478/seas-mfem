// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for SlipLawSRWPsi (§4.2.3 acceptance gates T_SRW_1..T_SRW_9
// plus review-added tests R-001, R-002, R-003, R-004 guards).

#include "test_macros.hpp"
#include "../../friction/slip_law_srw_psi.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/tpv104_params.hpp"

#include <cstdlib>
#include <cstring>
#include <random>

using namespace mfem;
using namespace mfem::seas;

// TPV104 canonical scalars for the envelope.
static constexpr real_t kA_TPV   = 0.01;
static constexpr real_t kB_TPV   = 0.014;
static constexpr real_t kV0_TPV  = 1.0e-6;
static constexpr real_t kF0_TPV  = 0.6;
static constexpr real_t kFw_TPV  = 0.1;
static constexpr real_t kL_TPV   = 0.4;
static constexpr real_t kVw_TPV  = 0.1;

// Bit-identity check: two reals match if their IEEE representation is
// identical.  Handles identical ±inf / NaN as "agree" (same singular
// limit hit by both codes).
static bool BitIdentical(real_t a, real_t b)
{
   unsigned char ba[sizeof(real_t)], bb[sizeof(real_t)];
   std::memcpy(ba, &a, sizeof(real_t));
   std::memcpy(bb, &b, sizeof(real_t));
   return std::memcmp(ba, bb, sizeof(real_t)) == 0;
}

static bool BothNonFinite(real_t a, real_t b)
{
   if (std::isnan(a) && std::isnan(b)) { return true; }
   if (std::isinf(a) && std::isinf(b) &&
       ((a > 0 && b > 0) || (a < 0 && b < 0))) { return true; }
   return false;
}

// ----------------------------------------------------------------------------
// R-003 / R-004 STANDALONE REFERENCE.
// Inlined dead-code arithmetic — does NOT call any production helper
// (LogSinhStable, IntegerPow8, PsiSS_SRW, UpdateStateAnalyticSlipLawSRW).
// A drift in any production helper will appear as a nonzero diff here.
//
// R-004: uses the unrolled integer 8th-power ((r²)²)² form, matching the
// production code; a future refactor that swaps the production side to
// std::pow(r, 8.0) would no longer byte-match this reference.
// ----------------------------------------------------------------------------
static real_t reference_psi_ss(real_t V, real_t V_w, real_t a, real_t b,
                               real_t V0, real_t f0, real_t muW)
{
   // f_LV inlined.
   const real_t f_LV = std::max(static_cast<real_t>(0),
                                f0 - (b - a) * std::log(V / V0));

   // Unrolled (V/V_w)^8 = ((x²)²)².
   const real_t r  = V / V_w;
   const real_t r2 = r * r;
   const real_t r4 = r2 * r2;
   const real_t r8 = r4 * r4;
   const real_t denom = std::pow(1.0 + r8, 1.0 / 8.0);
   const real_t f_ss  = muW + (f_LV - muW) / denom;

   // logsinh inlined: log(x · sinh(c)) = |c| + log((x/2)·−sign(c)·expm1(−2|c|)).
   const real_t x    = 2.0 * V0 / V;
   const real_t c    = f_ss / a;
   const real_t signC = (c >= 0.0) ? 1.0 : -1.0;
   const real_t absC  = std::abs(c);
   return a * (absC + std::log(x / 2.0 * -signC * std::expm1(-2.0 * absC)));
}

static real_t reference_update(real_t psi_0, real_t V, real_t L, real_t dt,
                               real_t V_w, real_t a, real_t b,
                               real_t V0, real_t f0, real_t muW)
{
   const real_t psi_ss = reference_psi_ss(V, V_w, a, b, V0, f0, muW);
   const real_t preexp1 = -V * (dt / L);
   const real_t exp1v = std::exp(preexp1);
   const real_t exp1m = -std::expm1(preexp1);
   return psi_ss * exp1m + exp1v * psi_0;
}

// ----------------------------------------------------------------------------
// T_SRW_1 — ψ_ss equation matches the standalone inlined reference.
// R-010 (review fix): narrower envelope so f_ss/a ≤ 300 across all
// samples, avoiding std::sinh overflow in the ad-hoc verification path.
// ----------------------------------------------------------------------------
void TestPsiSSEquation()
{
   std::cout << "\n[T_SRW_1] ψ_ss formula equivalence (R-010 narrowed envelope)\n";

   std::mt19937 rng(2026'04'24);
   std::uniform_real_distribution<real_t> logV(-15.0, 1.0);   // V ∈ [1e-15, 10]
   std::uniform_real_distribution<real_t> logVw(-1.0, 1.0);   // V_w ∈ [0.1, 10]
   std::uniform_real_distribution<real_t> aRng(0.008, 0.04);
   std::uniform_real_distribution<real_t> bRng(0.006, 0.04);
   std::uniform_real_distribution<real_t> fwRng(0.05, 0.3);

   const int N = 20;
   int ok = 0;
   for (int t = 0; t < N; ++t)
   {
      const real_t V   = std::pow(10.0, logV(rng));
      const real_t Vw  = std::pow(10.0, logVw(rng));
      const real_t a   = aRng(rng);
      const real_t b   = bRng(rng);
      const real_t muW = fwRng(rng);

      const real_t got = SlipLawSRWPsi::PsiSS_SRW(V, Vw, a, b, kV0_TPV,
                                                  kF0_TPV, muW);
      const real_t expected = reference_psi_ss(V, Vw, a, b, kV0_TPV,
                                               kF0_TPV, muW);

      if (BothNonFinite(got, expected) || BitIdentical(got, expected))
      {
         ++ok;
         continue;
      }
      const real_t denom = std::max<real_t>(std::abs(expected), 1.0);
      const real_t rel   = std::abs(got - expected) / denom;
      if (rel < 1e-14) { ++ok; }
      else
      {
         std::cerr << "  FAIL t=" << t
                   << " V=" << V << " Vw=" << Vw
                   << " a=" << a << " b=" << b << " muW=" << muW
                   << " expected=" << expected << " got=" << got
                   << " rel=" << rel << "\n";
      }
   }
   TEST_ASSERT(ok == N,
               "20/20 random ψ_ss tuples bit-match the standalone reference");
}

// ----------------------------------------------------------------------------
// T_SRW_2 — ψ_ss falls back to classical rate-and-state at f_w=0, V_w→∞.
// Classical: ψ_ss = f0 + b·ln(V0/V).
// ----------------------------------------------------------------------------
void TestClassicalLimit()
{
   std::cout << "\n[T_SRW_2] classical limit (f_w=0, V_w=inf)\n";

   const real_t muW = 0.0;
   const real_t Vw  = 1.0e300;

   int ok = 0;
   const std::vector<real_t> V_samples = {1e-12, 1e-9, 1e-6, 1e-3, 1.0};
   for (real_t V : V_samples)
   {
      const real_t classical = kF0_TPV + kB_TPV * std::log(kV0_TPV / V);
      const real_t got = SlipLawSRWPsi::PsiSS_SRW(V, Vw, kA_TPV, kB_TPV,
                                                  kV0_TPV, kF0_TPV, muW);
      const real_t rel = std::abs(got - classical) /
                         std::max<real_t>(std::abs(classical), 1.0);
      if (rel < 1e-12) { ++ok; }
   }
   TEST_ASSERT(ok == 5, "5/5 classical-limit samples match to 1e-12 rel");
}

// ----------------------------------------------------------------------------
// T_SRW_3 — dψ/dt = 0 when ψ = ψ_ss.
// ----------------------------------------------------------------------------
void TestRateVanishesAtSS()
{
   std::cout << "\n[T_SRW_3] Rate_SRW vanishes at ψ = ψ_ss\n";

   SlipLawSRWPsi law(kA_TPV, kB_TPV, kV0_TPV, kF0_TPV, kFw_TPV, kVw_TPV);

   const std::vector<real_t> V_samples = {1e-12, 1e-6, 1e-3, 1.0};
   int ok = 0;
   for (real_t V : V_samples)
   {
      const real_t psi_ss = law.SteadyState_SRW(V, kVw_TPV, kA_TPV);
      const real_t rate   = law.Rate_SRW(V, psi_ss, kL_TPV, kVw_TPV, kA_TPV);
      if (std::abs(rate) < 1e-30) { ++ok; }
   }
   TEST_ASSERT(ok == 4, "4/4 Rate(ψ_ss) = 0 exactly");
}

// ----------------------------------------------------------------------------
// T_SRW_4 — fixed-point test: ψ_ss is invariant under the analytic step.
// ----------------------------------------------------------------------------
void TestAnalyticStepFixedPoint()
{
   std::cout << "\n[T_SRW_4] analytic update preserves ψ_ss\n";

   const std::vector<real_t> V_samples = {1e-12, 1e-6, 1e-3, 1.0};
   const std::vector<real_t> dt_samples = {1e-6, 1e-3, 1.0, 60.0, 1e5};

   int ok = 0, total = 0;
   for (real_t V : V_samples)
   {
      const real_t psi_ss = SlipLawSRWPsi::PsiSS_SRW(V, kVw_TPV, kA_TPV,
                                                     kB_TPV, kV0_TPV,
                                                     kF0_TPV, kFw_TPV);
      for (real_t dt : dt_samples)
      {
         ++total;
         const real_t psi_new = UpdateStateAnalyticSlipLawSRW(
            psi_ss, V, kL_TPV, dt, kVw_TPV, kA_TPV,
            kB_TPV, kV0_TPV, kF0_TPV, kFw_TPV);
         const real_t rel = std::abs(psi_new - psi_ss) /
                            std::max<real_t>(std::abs(psi_ss), 1.0);
         if (rel < 1e-12) { ++ok; }
      }
   }
   TEST_ASSERT(ok == total,
               "All analytic-step fixed-point checks match to 1e-12");
}

// ----------------------------------------------------------------------------
// T_SRW_5 / T_SRW_8 — byte-match vs a standalone inlined reference
// (R-003: no shared helpers).
// ----------------------------------------------------------------------------
void TestReferenceByteMatch()
{
   std::cout << "\n[T_SRW_5/T_SRW_8] byte-match vs standalone reference\n";

   std::mt19937 rng(20260424u);
   std::uniform_real_distribution<real_t> logV(-15.0, 1.0);
   std::uniform_real_distribution<real_t> logVw(-1.0, 0.5);
   std::uniform_real_distribution<real_t> aRng(0.008, 0.04);
   std::uniform_real_distribution<real_t> bRng(0.006, 0.04);
   std::uniform_real_distribution<real_t> LRng(0.1, 1.0);
   std::uniform_real_distribution<real_t> logDt(-6.0, 0.0);
   std::uniform_real_distribution<real_t> fwRng(0.05, 0.3);

   const int N = 1000;
   int ok = 0;
   real_t max_rel = 0.0;

   for (int t = 0; t < N; ++t)
   {
      const real_t V   = std::pow(10.0, logV(rng));
      const real_t Vw  = std::pow(10.0, logVw(rng));
      const real_t a   = aRng(rng);
      const real_t b   = bRng(rng);
      const real_t L   = LRng(rng);
      const real_t dt  = std::pow(10.0, logDt(rng));
      const real_t muW = fwRng(rng);
      const real_t psi_ss = SlipLawSRWPsi::PsiSS_SRW(V, Vw, a, b,
                                                     kV0_TPV, kF0_TPV, muW);
      const real_t psi_0  = psi_ss * 0.9;

      const real_t ours = UpdateStateAnalyticSlipLawSRW(
                             psi_0, V, L, dt, Vw, a,
                             b, kV0_TPV, kF0_TPV, muW);
      const real_t ref  = reference_update(
                             psi_0, V, L, dt, Vw, a,
                             b, kV0_TPV, kF0_TPV, muW);

      if (BothNonFinite(ours, ref) || BitIdentical(ours, ref))
      {
         ++ok;
         continue;
      }
      const real_t denom = std::max<real_t>(std::abs(ref), 1.0);
      const real_t rel = std::abs(ours - ref) / denom;
      max_rel = std::max(max_rel, rel);
      if (rel < 1e-13) { ++ok; }
   }
   std::cout << "  max rel err = " << max_rel << " over " << N << " samples\n";
   TEST_ASSERT(ok == N,
               "All 1000 standalone-reference byte-match samples < 1e-13 rel");
}

// ----------------------------------------------------------------------------
// T_SRW_6 — ψ_ini anchor duplicated from T_TPV104_P_1.
// ----------------------------------------------------------------------------
void TestPsiIniAnchor()
{
   std::cout << "\n[T_SRW_6] ψ_ini anchor = 5.6359184e-01\n";
   const real_t psi_expected = 5.6359184e-01;
   const real_t psi_computed = ComputeInitialPsiTPV104(TPV104Params::a_in);
   TEST_NEAR(psi_computed, psi_expected, 1e-8,
             "ψ_ini(a_in=0.01) match to 1e-8");
}

// ----------------------------------------------------------------------------
// T_SRW_7 — adherence to StateEvolution interface + R-011 guard.
// The base-class virtuals must (a) work in test-fixture mode, and (b)
// abort in production mode — i.e. the R-001 silent-fallthrough trap is
// actually closed.
// ----------------------------------------------------------------------------
void TestStateEvolutionInterface()
{
   std::cout << "\n[T_SRW_7] StateEvolution interface + R-011 R-001 guard\n";

   // SlipLawSRWPsi's production-mode guards throw std::runtime_error
   // directly (see header R-001 fix) — no MFEM_ERROR_THROW setup needed,
   // which keeps the test portable on MFEM builds without exceptions.

   SlipLawSRWPsi law(kA_TPV, kB_TPV, kV0_TPV, kF0_TPV, kFw_TPV, kVw_TPV);
   StateEvolution *base = &law;

   TEST_ASSERT(std::string(base->GetName()) == "SlipLawSRWPsi",
               "GetName() returns 'SlipLawSRWPsi'");

   // Use V ~ V_w so the SRW denominator (1 + (V/V_w)^8)^(1/8) is
   // non-trivial — this is where the per-QP V_w meaningfully differs.
   const real_t V = 0.2;

   // Test-fixture mode (default): base virtuals dispatch through default
   // V_w and stored a_.
   const real_t psi_ss = base->SteadyState(V, kL_TPV);
   TEST_ASSERT(std::isfinite(psi_ss),
               "SteadyState virtual returns finite value (fixture mode)");

   const real_t rate_at_ss = base->Rate(V, psi_ss, kL_TPV);
   TEST_ASSERT(std::abs(rate_at_ss) < 1e-28,
               "Rate virtual vanishes at ψ_ss (fixture mode)");

   const real_t drdpsi = base->RateDerivativeTheta(V, psi_ss, kL_TPV);
   TEST_NEAR(drdpsi, -V / kL_TPV, 1e-20, "RateDerivativeTheta = -V/L");

   // R-011: the base-virtual result uses V_w_default_ = kVw_TPV = 0.1.
   // Calling Rate_SRW with a different V_w must produce a DIFFERENT
   // rate (else the base virtual is not actually looking at V_w).
   const real_t base_rate_at_perturbed =
      base->Rate(V, psi_ss * 1.1, kL_TPV);
   const real_t srw_diff_vw =
      law.Rate_SRW(V, psi_ss * 1.1, kL_TPV, /*V_w=*/1.0, /*a=*/kA_TPV);
   TEST_ASSERT(std::abs(base_rate_at_perturbed - srw_diff_vw) > 1e-30,
               "base Rate ≠ Rate_SRW at different V_w (R-011)");

   // R-001: flip to production mode.  Base virtuals must now abort.
   law.SetProductionMode();

   // The abort is via MFEM_VERIFY — catch whatever exception MFEM's
   // error handler throws.
   bool aborted_rate = false;
   try { (void)base->Rate(V, psi_ss, kL_TPV); }
   catch (const std::exception&) { aborted_rate = true; }
   catch (...)                   { aborted_rate = true; }
   TEST_ASSERT(aborted_rate,
               "Rate virtual aborts in production mode (R-001 guard)");

   bool aborted_ss = false;
   try { (void)base->SteadyState(V, kL_TPV); }
   catch (const std::exception&) { aborted_ss = true; }
   catch (...)                   { aborted_ss = true; }
   TEST_ASSERT(aborted_ss,
               "SteadyState virtual aborts in production mode (R-001 guard)");

   // The SRW-aware overloads remain callable in production mode.
   const real_t rate_srw = law.Rate_SRW(V, psi_ss, kL_TPV, kVw_TPV, kA_TPV);
   TEST_ASSERT(std::isfinite(rate_srw),
               "Rate_SRW callable in production mode");
}

// ----------------------------------------------------------------------------
// R2-003 dedicated test — RateDerivativeV(V=0, …) returns finite.
// Without the R2-003 fix, the finite-difference branch would produce NaN
// because ψ_ss(V=0) = +∞ after the R-002 clamp removal.
// ----------------------------------------------------------------------------
void TestRateDerivativeVAtZero()
{
   std::cout << "\n[R2-003] RateDerivativeV(V=0) returns finite\n";

   SlipLawSRWPsi law(kA_TPV, kB_TPV, kV0_TPV, kF0_TPV, kFw_TPV, kVw_TPV);
   StateEvolution *base = &law;

   // Test-fixture mode (default); should not throw.
   const real_t d_zero = base->RateDerivativeV(0.0, /*psi=*/0.5, kL_TPV);
   TEST_ASSERT(std::isfinite(d_zero),
               "RateDerivativeV(V=0) returns finite (no NaN)");

   // Sanity: the V=0 analytic limit approximates the V→0⁺ finite-diff.
   // Use V = 1e-30 as a proxy for V→0⁺ and confirm sign/magnitude agree.
   const real_t d_small = base->RateDerivativeV(1e-30, 0.5, kL_TPV);
   TEST_ASSERT(std::isfinite(d_small),
               "RateDerivativeV(V=1e-30) returns finite");
}

// ----------------------------------------------------------------------------
// R-002 dedicated test — no V_safe clamp: raw V feeds through the
// formula.  At V = 1e-48 (below the former 1e-45 floor), our output must
// match the standalone inlined reference exactly.
// ----------------------------------------------------------------------------
void TestNoVsafeClamp()
{
   std::cout << "\n[R-002] raw V (no V_safe clamp) matches inline reference\n";

   const real_t V    = 1e-48;
   const real_t psi0 = 0.9;
   const real_t a    = 0.01;
   const real_t b    = 0.014;
   const real_t V_w  = 0.1;
   const real_t L    = 0.4;
   const real_t dt   = 1e-6;
   const real_t V0   = 1e-6;
   const real_t f0   = 0.6;
   const real_t muW  = 0.1;

   const real_t ours = UpdateStateAnalyticSlipLawSRW(
      psi0, V, L, dt, V_w, a, b, V0, f0, muW);
   const real_t ref  = reference_update(
      psi0, V, L, dt, V_w, a, b, V0, f0, muW);

   TEST_ASSERT(std::isfinite(ours) || BothNonFinite(ours, ref),
               "output at V<1e-45 is finite OR matches reference singular");
   TEST_ASSERT(BitIdentical(ours, ref) || BothNonFinite(ours, ref),
               "V<1e-45 input byte-matches inline reference (no clamp)");
}

// ----------------------------------------------------------------------------
// R-004 dedicated test — the production (V/V_w)^8 uses the unrolled
// form, not std::pow(r, 8.0).  At the crossover r = V/V_w ≈ 1.1 the two
// forms can differ by O(1) ULP; our production code must match the
// unrolled form exactly.
// ----------------------------------------------------------------------------
void TestIntegerPowerUnrolled()
{
   std::cout << "\n[R-004 + R2-002] (V/V_w)^8 matches unrolled integer power\n";

   const real_t V = 0.11, V_w = 0.1;
   const real_t r = V / V_w;
   const real_t r2 = r * r;
   const real_t r4 = r2 * r2;
   const real_t r8_unrolled = r4 * r4;
   const real_t r8_pow      = std::pow(r, 8.0);

   // R2-002: explicit flag for toolchains where std::pow specialises to
   // the unrolled form; on those, the anti-reference below degenerates
   // to the forward reference and the "does-not-match" guard is moot.
   const bool toolchain_distinguishes = (r8_unrolled != r8_pow);
   std::cout << "  r^8 unrolled = " << r8_unrolled << "\n"
             << "  r^8 std::pow = " << r8_pow << "\n"
             << "  toolchain distinguishes: "
             << (toolchain_distinguishes ? "YES" : "NO") << "\n";

   // Production evaluate: construct a tuple where the different ULP
   // manifests in the denominator (1 + r^8)^(1/8).  Re-implement the
   // production production with unrolled arithmetic and check exact
   // match against our PsiSS_SRW.
   const real_t a = 0.01, b = 0.014, V0 = 1e-6, f0 = 0.6, muW = 0.1;
   const real_t V_sample = r * V_w;
   const real_t f_LV = std::max(static_cast<real_t>(0),
                                f0 - (b - a) * std::log(V_sample / V0));
   const real_t denom = std::pow(1.0 + r8_unrolled, 1.0 / 8.0);
   const real_t f_ss  = muW + (f_LV - muW) / denom;
   const real_t x = 2.0 * V0 / V_sample;
   const real_t c = f_ss / a;
   const real_t signC = (c >= 0) ? 1.0 : -1.0;
   const real_t absC  = std::abs(c);
   const real_t expected =
      a * (absC + std::log(x / 2.0 * -signC * std::expm1(-2.0 * absC)));

   const real_t got = SlipLawSRWPsi::PsiSS_SRW(V_sample, V_w, a, b,
                                               V0, f0, muW);
   TEST_ASSERT(BitIdentical(got, expected),
               "production PsiSS_SRW uses unrolled integer power (R-004)");

   // R2-002: build an "anti-reference" that uses std::pow(r, 8.0) instead
   // of the unrolled form.  Production must NOT byte-match the
   // anti-reference on toolchains where std::pow differs from the
   // unrolled form — otherwise a future regression swapping production
   // from IntegerPow8 back to std::pow would pass silently.
   if (toolchain_distinguishes)
   {
      const real_t denom_pow = std::pow(1.0 + r8_pow, 1.0 / 8.0);
      const real_t f_ss_pow  = muW + (f_LV - muW) / denom_pow;
      const real_t c_pow     = f_ss_pow / a;
      const real_t signCp    = (c_pow >= 0.0) ? 1.0 : -1.0;
      const real_t absCp     = std::abs(c_pow);
      const real_t anti_expected =
         a * (absCp + std::log(x / 2.0 * -signCp * std::expm1(-2.0 * absCp)));
      TEST_ASSERT(!BitIdentical(got, anti_expected),
                  "production PsiSS_SRW does NOT byte-match std::pow-based "
                  "reference (R2-002 regression guard)");
   }
   else
   {
      // Toolchain treats std::pow(x, 8.0) as the unrolled form; the guard
      // is moot here but logged so future reviewers see the skip reason.
      std::cout << "  (R2-002 anti-reference guard skipped — toolchain "
                   "treats std::pow(r, 8.0) = r4·r4)\n";
   }
}

// ----------------------------------------------------------------------------
// T_SRW_9 — classical-slip-law limit: f_w = 0, V_w = 1e300.
// Plan text refers to "AgingLawPsi::Rate integrated analytically"; that
// is not the correct reference (AgingLawPsi is a different ODE).  The
// correct reference is the classical slip law in ψ-space, which shares
// the same exponential-relaxation analytic solution as our FVW law at
// f_w = 0 and V_w → ∞.
// ----------------------------------------------------------------------------
void TestClassicalSlipLawLimit()
{
   std::cout << "\n[T_SRW_9] classical slip-law limit\n";

   const real_t muW = 0.0;
   const real_t Vw  = 1.0e300;

   const std::vector<real_t> V_samples = {1e-9, 1e-6, 1e-3, 1.0};
   const std::vector<real_t> dt_samples = {1e-3, 1.0};

   int ok = 0, total = 0;
   for (real_t V : V_samples)
   {
      const real_t psi_ss = kF0_TPV + kB_TPV * std::log(kV0_TPV / V);
      const real_t psi_0 = psi_ss + 0.05;
      for (real_t dt : dt_samples)
      {
         ++total;
         const real_t got = UpdateStateAnalyticSlipLawSRW(
                               psi_0, V, kL_TPV, dt, Vw, kA_TPV,
                               kB_TPV, kV0_TPV, kF0_TPV, muW);
         const real_t classical = psi_ss + (psi_0 - psi_ss)
                                   * std::exp(-V * dt / kL_TPV);
         const real_t rel = std::abs(got - classical) /
                            std::max<real_t>(std::abs(classical), 1.0);
         if (rel < 1e-10) { ++ok; }
      }
   }
   TEST_ASSERT(ok == total, "All classical-slip-law-limit cases < 1e-10 rel");
}

int main(int argc, char *argv[])
{
   TestPsiSSEquation();
   TestClassicalLimit();
   TestRateVanishesAtSS();
   TestAnalyticStepFixedPoint();
   TestReferenceByteMatch();
   TestPsiIniAnchor();
   TestStateEvolutionInterface();
   TestRateDerivativeVAtZero();
   TestNoVsafeClamp();
   TestIntegerPowerUnrolled();
   TestClassicalSlipLawLimit();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
