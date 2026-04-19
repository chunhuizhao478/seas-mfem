// Phase 4a setup unit tests for TPV102 benchmark parameters.

#include "mfem.hpp"
#include "../../config/tpv102_params.hpp"
#include <iostream>
#include <cmath>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) < t_) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ << ", diff " << std::abs(v_-e_) << ")\n"; } \
} while(0)

// ===== Test 1: Boxcar function (SCEC Eq. 5, tanh transition) =====
void TestBoxcarFunction()
{
   std::cout << "Test 1: TestBoxcarFunction\n";

   real_t W = 15e3, w = 3e3;

   // Center: B(0, W, w) = 1
   TEST_NEAR(Boxcar(0.0, W, w), 1.0, 1e-12, "B(0, W, w) = 1");

   // Well inside: B(10km, 15km, 3km) = 1 (10 < W = 15)
   TEST_NEAR(Boxcar(10e3, W, w), 1.0, 1e-12, "B(10km) = 1 (inside)");

   // At edge of flat zone: B(15km, 15km, 3km) = 1 (SCEC: |x| <= W → 1)
   TEST_NEAR(Boxcar(15e3, W, w), 1.0, 1e-12, "B(15km) = 1 (at W boundary)");

   // In transition: B(16.5km, 15km, 3km) should be between 0 and 1
   real_t B_mid = Boxcar(16.5e3, W, w);
   TEST_ASSERT(B_mid > 0.0 && B_mid < 1.0,
               "B(16.5km) in transition (" + std::to_string(B_mid) + ")");

   // At outer edge: B(18km, 15km, 3km) = 0 (exactly at W+w)
   TEST_NEAR(Boxcar(18e3, W, w), 0.0, 1e-12, "B(18km) = 0 (outer edge)");

   // Well outside: B(25km, 15km, 3km) = 0
   TEST_NEAR(Boxcar(25e3, W, w), 0.0, 1e-12, "B(25km) = 0 (outside)");

   // Negative x: B(-5km, 15km, 3km) = 1 (symmetric)
   TEST_NEAR(Boxcar(-5e3, W, w), 1.0, 1e-12, "B(-5km) = 1 (symmetric)");
}

// ===== Test 2: Spatially varying a parameter =====
void TestSpatiallyVaryingA()
{
   std::cout << "Test 2: TestSpatiallyVaryingA\n";

   // At hypocenter (0, 7.5 km): center of VW zone → a = a_vw = 0.008
   real_t a_hypo = ComputeA(0.0, 7.5e3);
   TEST_NEAR(a_hypo, TPV102Params::a_vw, 1e-12,
             "a(0, 7.5km) = a_vw = 0.008");

   // SCEC: VW covers entire fault depth (0-15 km) and along-strike (-15 to 15 km)
   // At fault bottom (0, 15km): still VW (dip Boxcar centered at 7.5, half-width 7.5)
   real_t a_bottom = ComputeA(0.0, 15e3);
   TEST_NEAR(a_bottom, TPV102Params::a_vw, 1e-12,
             "a(0, 15km) = a_vw (SCEC: VW covers entire fault depth)");

   // At fault edge along-strike (14km, 7.5km): still VW (strike Boxcar W=15km)
   real_t a_edge = ComputeA(14e3, 7.5e3);
   TEST_NEAR(a_edge, TPV102Params::a_vw, 1e-12,
             "a(14km, 7.5km) = a_vw (SCEC: VW covers entire fault along-strike)");

   // Deep VS zone (0, 20km): outside VW → a = a_vs = 0.016
   real_t a_deep = ComputeA(0.0, 20e3);
   TEST_NEAR(a_deep, TPV102Params::a_vs, 1e-6,
             "a(0, 20km) = a_vs = 0.016");

   // Far along-strike (25km, 7.5km): outside VW along-strike
   real_t a_far = ComputeA(25e3, 7.5e3);
   TEST_NEAR(a_far, TPV102Params::a_vs, 1e-6,
             "a(25km, 7.5km) = a_vs = 0.016");

   // Surface (z=0): a = a_vw (VW extends to surface per SCEC spec)
   real_t a_surface = ComputeA(0.0, 0.0);
   TEST_NEAR(a_surface, TPV102Params::a_vw, 1e-12,
             "a(0, surface) = a_vw = 0.008 (VW at surface, per SCEC spec)");

   // Transition outside fault: a(0, 17km) should be between a_vw and a_vs
   real_t a_trans = ComputeA(0.0, 17e3);
   TEST_ASSERT(a_trans > TPV102Params::a_vw && a_trans < TPV102Params::a_vs,
               "a(0, 17km) in transition outside fault (" + std::to_string(a_trans) + ")");
}

// ===== Test 3: Initial state equilibrium =====
void TestInitialStateEquilibrium()
{
   std::cout << "Test 3: TestInitialStateEquilibrium\n";

   real_t a = TPV102Params::a_vw;
   real_t psi = ComputeInitialPsi(a);

   // Verify: tau_ini = sigma_n * a * asinh(V_ini/(2*V0) * exp(psi/a))
   real_t C = std::exp(psi / a) / (2.0 * TPV102Params::V0);
   real_t f = a * std::asinh(TPV102Params::V_ini * C);
   real_t tau_check = TPV102Params::sigma_n * f;

   real_t rel = std::abs(tau_check - TPV102Params::tau_ini) / TPV102Params::tau_ini;
   TEST_ASSERT(rel < 1e-10,
               "Initial psi satisfies equilibrium (rel " + std::to_string(rel) +
               ", psi = " + std::to_string(psi) + ")");

   // psi should be positive and finite
   TEST_ASSERT(std::isfinite(psi) && psi > 0,
               "Initial psi is finite and positive (" + std::to_string(psi) + ")");
}

// ===== Test 4: Nucleation perturbation =====
void TestNucleationPerturbation()
{
   std::cout << "Test 4: TestNucleationPerturbation\n";

   // At hypocenter, t = 0: dtau = 0 (temporal factor is 0 at t=0)
   real_t dtau_t0 = NucleationPerturbation(0.0, 7.5e3, 0.0);
   TEST_NEAR(dtau_t0, 0.0, 1e-6, "Nucleation: dtau = 0 at t=0");

   // At hypocenter, t > T: dtau = dtau0 * F(0) = dtau0 * exp(0) = dtau0 * 1 = dtau0
   real_t dtau_full = NucleationPerturbation(0.0, 7.5e3, 2.0);
   // F(0) = exp(0/(0-R²)) = exp(0) = 1
   TEST_NEAR(dtau_full, TPV102Params::nuc_dtau, 1e-6,
             "Nucleation: dtau = dtau0 at r=0, t>T (" + std::to_string(dtau_full) + ")");

   // Outside nucleation zone: dtau = 0
   real_t dtau_outside = NucleationPerturbation(0.0, 7.5e3 + 5e3, 2.0);
   TEST_NEAR(dtau_outside, 0.0, 1e-6, "Nucleation: dtau = 0 for r > R");

   // At t = T/2: intermediate temporal value
   real_t dtau_half = NucleationPerturbation(0.0, 7.5e3, 0.5);
   TEST_ASSERT(dtau_half > 0 && dtau_half < TPV102Params::nuc_dtau,
               "Nucleation: 0 < dtau < dtau0 at t=T/2 (" + std::to_string(dtau_half) + ")");
}

// ===== Test 5: Material parameter consistency =====
void TestMaterialConsistency()
{
   std::cout << "Test 5: TestMaterialConsistency\n";

   // Verify wave speeds from Lame parameters
   real_t cp_check = std::sqrt((TPV102Params::lambda + 2*TPV102Params::mu) / TPV102Params::rho);
   real_t cs_check = std::sqrt(TPV102Params::mu / TPV102Params::rho);

   TEST_NEAR(cp_check, TPV102Params::cp, 1.0,
             "cp consistent: " + std::to_string(cp_check));
   TEST_NEAR(cs_check, TPV102Params::cs, 1.0,
             "cs consistent: " + std::to_string(cs_check));

   // Verify Poisson's ratio ≈ 0.25
   real_t nu = TPV102Params::lambda / (2.0 * (TPV102Params::lambda + TPV102Params::mu));
   TEST_NEAR(nu, 0.25, 0.01, "Poisson's ratio ≈ 0.25 (" + std::to_string(nu) + ")");

   // Verify eta_s
   TEST_NEAR(TPV102Params::eta_s, TPV102Params::Zs / 2.0, 1e-6,
             "eta_s = Zs/2 (homogeneous)");
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 Setup Unit Tests (Phase 4a)\n";
   std::cout << "========================================\n\n";

   TestBoxcarFunction();
   TestSpatiallyVaryingA();
   TestInitialStateEquilibrium();
   TestNucleationPerturbation();
   TestMaterialConsistency();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";

   return (num_failed > 0) ? 1 : 0;
}
