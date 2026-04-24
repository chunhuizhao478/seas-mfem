// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for friction_coeff_stable.hpp (§4.10 Step 4 gates
// T_TPV104_FC_1..4 + review-added R-005 test).

#include "test_macros.hpp"
#include "../../friction/friction_coeff_stable.hpp"
#include "../../friction/dieterich_ruina.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <string>

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::friction_stable;

// ---------------------------------------------------------------------------
// R-003 STANDALONE REFERENCE: fully inlined straight-line arithmetic for
// asinh(x · exp(c)) using the same branch identity as the production
// header, but calling only std:: primitives — no dispatch into
// friction_stable::* functions.  A regression in ArsinhExp / ComputeCExp /
// the switch thresholds in the production code will manifest as a
// nonzero diff here.
// ---------------------------------------------------------------------------
static real_t reference_arsinhexp(real_t x, real_t cExpLog)
{
   constexpr real_t Switch    = 10.0;
   constexpr real_t Threshold = 50.0;
   constexpr real_t Log2      = 0.69314718055994530943;
   int xexp = 0;
   (void)std::frexp(x, &xexp);

   // cExp = exp(-|cExpLog|) via the same sign branch as production.
   real_t cExp;
   if (cExpLog > 0.0) { cExp = std::exp(-cExpLog); }
   else               { cExp = std::exp(cExpLog);  }

   if (cExpLog + std::max(xexp, 0) * Log2 > Switch || cExpLog >= Threshold)
   {
      // Stable branch: asinh(x·e^c) = sign(x)·(c + log(|x| + sqrt(x² + e^{-2c}))).
      real_t e = cExp;
      if (cExpLog <= 0.0) { e = 1.0 / e; }
      const real_t xa = std::abs(x);
      const real_t xs = (x >= 0.0) ? 1.0 : -1.0;
      return xs * (cExpLog + std::log(xa + std::sqrt(xa * xa + e * e)));
   }
   else
   {
      // Normal branch: v = e^c · x; return asinh(v).
      real_t e = cExp;
      if (cExpLog > 0.0) { e = 1.0 / e; }
      const real_t v = e * x;
      return std::asinh(v);
   }
}

// ----------------------------------------------------------------------------
// T_TPV104_FC_1 — standalone-reference byte-match on 10 000 samples.
// ----------------------------------------------------------------------------
void TestReferenceByteMatch()
{
   std::cout << "\n[T_TPV104_FC_1] friction coefficient ↔ standalone reference\n";

   std::mt19937 rng(20260424u);
   std::uniform_real_distribution<real_t> logV(-18.0, 1.0);
   std::uniform_real_distribution<real_t> psiRng(0.1, 2.0);
   std::uniform_real_distribution<real_t> aRng(0.005, 0.04);

   const real_t V0 = 1.0e-6;

   int ok = 0;
   real_t max_rel = 0.0;
   for (int t = 0; t < 10000; ++t)
   {
      const real_t V   = std::pow(10.0, logV(rng));
      const real_t psi = psiRng(rng);
      const real_t a   = aRng(rng);

      const real_t ours = FrictionCoefficientStable(V, psi, a, V0);
      const real_t ref  = a * reference_arsinhexp(V / (2.0 * V0), psi / a);

      if (std::isnan(ours) && std::isnan(ref))   { ++ok; continue; }
      if (std::isinf(ours) && std::isinf(ref) &&
          ((ours > 0 && ref > 0) || (ours < 0 && ref < 0))) { ++ok; continue; }

      const real_t denom = std::max<real_t>(std::abs(ref), 1.0);
      const real_t rel = std::abs(ours - ref) / denom;
      max_rel = std::max(max_rel, rel);
      if (rel < 1e-15) { ++ok; }
   }
   std::cout << "  max rel err = " << max_rel << " / 10000 samples\n";
   TEST_ASSERT(ok == 10000,
               "All 10000 (V, ψ, a) samples byte-match to 1e-15 rel");
}

// ----------------------------------------------------------------------------
// T_TPV104_FC_2 — equivalence with MFEM's existing FrictionCoefficientPsi
// in the well-conditioned envelope ψ/a ∈ [28, 80].
// ----------------------------------------------------------------------------
void TestMFEMEquivalence()
{
   std::cout << "\n[T_TPV104_FC_2] equivalence with MFEM FrictionCoefficientPsi\n";

   DieterichRuinaFriction dr(
      DieterichRuinaFriction::Constants{1e-6, 0.6, 0.012, 0.02});

   std::mt19937 rng(42u);
   std::uniform_real_distribution<real_t> logV(-18.0, 1.0);
   std::uniform_real_distribution<real_t> ratioRng(28.0, 80.0);
   std::uniform_real_distribution<real_t> aRng(0.005, 0.04);

   int ok = 0;
   real_t max_rel = 0.0;
   for (int t = 0; t < 500; ++t)
   {
      const real_t V   = std::pow(10.0, logV(rng));
      const real_t a   = aRng(rng);
      const real_t psi = ratioRng(rng) * a;

      const real_t stable  = FrictionCoefficientStable(V, psi, a, 1e-6);
      const real_t mfem_qd = dr.FrictionCoefficientPsi(V, psi, a);

      const real_t denom = std::max<real_t>(std::abs(mfem_qd), 1.0);
      const real_t rel = std::abs(stable - mfem_qd) / denom;
      max_rel = std::max(max_rel, rel);
      if (rel < 1e-12) { ++ok; }
   }
   std::cout << "  max rel err = " << max_rel << " / 500 samples\n";
   TEST_ASSERT(ok == 500,
               "stable μ equals MFEM FrictionCoefficientPsi to 1e-12 "
               "inside ψ/a ∈ [28, 80]");
}

// ----------------------------------------------------------------------------
// T_TPV104_FC_3 — derivative matches forward FD.
// ----------------------------------------------------------------------------
void TestDerivativeFD()
{
   std::cout << "\n[T_TPV104_FC_3] derivative byte-match vs forward FD\n";

   std::mt19937 rng(7u);
   std::uniform_real_distribution<real_t> logV(-12.0, 0.0);
   std::uniform_real_distribution<real_t> psiRng(0.3, 1.0);
   std::uniform_real_distribution<real_t> aRng(0.008, 0.02);

   const real_t V0 = 1.0e-6;

   int ok = 0;
   real_t max_rel = 0.0;
   for (int t = 0; t < 10000; ++t)
   {
      const real_t V   = std::pow(10.0, logV(rng));
      const real_t psi = psiRng(rng);
      const real_t a   = aRng(rng);

      const real_t h = V * 1e-6;
      const real_t mu_plus = FrictionCoefficientStable(V + h, psi, a, V0);
      const real_t mu      = FrictionCoefficientStable(V,     psi, a, V0);
      const real_t fd      = (mu_plus - mu) / h;

      const real_t analytic = FrictionCoefficientStableDerivV(V, psi, a, V0);

      const real_t denom = std::max<real_t>(std::abs(analytic),
                                            static_cast<real_t>(1e-30));
      const real_t rel = std::abs(analytic - fd) / denom;
      max_rel = std::max(max_rel, rel);
      if (rel < 1e-5) { ++ok; }
   }
   std::cout << "  max rel err = " << max_rel << " / 10000 samples\n";
   TEST_ASSERT(ok >= 9990,
               "≥ 99.9% of analytic derivatives match FD to 1e-5");
}

// ----------------------------------------------------------------------------
// T_TPV104_FC_4 — static grep: no literal `700` constant in the header.
// ----------------------------------------------------------------------------
void TestNoHardCoded700()
{
   std::cout << "\n[T_TPV104_FC_4] no hard-coded 700 in friction_coeff_stable.hpp\n";

   std::ifstream in("friction/friction_coeff_stable.hpp");
   if (!in) { in.open("../friction/friction_coeff_stable.hpp"); }
   if (!in) { in.open("../../friction/friction_coeff_stable.hpp"); }
   TEST_ASSERT(in.is_open(),
               "can locate friction/friction_coeff_stable.hpp for grep");

   std::ostringstream oss;
   oss << in.rdbuf();
   const std::string body = oss.str();

   const size_t code_start = body.find("namespace friction_stable");
   std::string code = body.substr(code_start != std::string::npos
                                  ? code_start : 0);
   size_t pos = 0;
   int hits = 0;
   while ((pos = code.find("700", pos)) != std::string::npos)
   {
      const size_t line_start = code.rfind('\n', pos);
      const size_t line_begin = (line_start == std::string::npos)
                                ? 0 : line_start + 1;
      const std::string line = code.substr(line_begin,
                                           code.find('\n', pos) - line_begin);
      const size_t cs = line.find("//");
      const bool in_comment = cs != std::string::npos
                              && cs <= (pos - line_begin);
      if (!in_comment) { ++hits; }
      pos += 3;
   }
   TEST_ASSERT(hits == 0,
               "0 literal-700 constants outside comments");
}

// ----------------------------------------------------------------------------
// R-005 dedicated test — no V≤0 early-return.  At V = 0, asinh branch
// evaluates  (V/2V0)·exp(ψ/a) = 0  →  asinh(0) = 0  →  μ = 0 algebraically.
// At V < 0, asinh returns a negative value, matching the canonical FVW
// reference (which has no V≤0 guard).
// ----------------------------------------------------------------------------
void TestNoVNonPositiveGuard()
{
   std::cout << "\n[R-005] no V≤0 early-return in FrictionCoefficientStable\n";

   const real_t psi = 0.5, a = 0.01, V0 = 1e-6;

   // V = 0 → μ = 0 by algebraic identity (asinh(0) = 0), NOT by guard.
   const real_t mu_zero = FrictionCoefficientStable(0.0, psi, a, V0);
   TEST_NEAR(mu_zero, 0.0, 1e-15,
             "μ(V=0) = 0 via asinh branch, not early-return");

   // V < 0 → μ < 0 (matches reference: no guard).
   const real_t mu_neg = FrictionCoefficientStable(-1e-9, psi, a, V0);
   TEST_ASSERT(mu_neg < 0.0,
               "μ(V<0) < 0 (no early-return coerces to 0)");

   // μ odd in V: μ(-V) = -μ(V) for the stable branch.
   const real_t mu_pos = FrictionCoefficientStable(+1e-9, psi, a, V0);
   TEST_NEAR(mu_neg + mu_pos, 0.0, 1e-14,
             "μ is an odd function of V (no V≤0 guard)");
}

int main(int argc, char *argv[])
{
   TestReferenceByteMatch();
   TestMFEMEquivalence();
   TestDerivativeFD();
   TestNoHardCoded700();
   TestNoVNonPositiveGuard();
   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
