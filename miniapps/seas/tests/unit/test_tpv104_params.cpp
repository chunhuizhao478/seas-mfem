// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for TPV104 parameter header (§4.1 of tpv104 debug plan).
//
// Tests T_TPV104_P_1..P_5 from §4.10 Step 1.

#include "test_macros.hpp"
#include "../../config/tpv104_params.hpp"

#include <cstdlib>

using namespace mfem;
using namespace mfem::seas;

// ----------------------------------------------------------------------------
// T_TPV104_P_1 — ψ_ini anchor match to the TPV104 reference benchmark
// trace col-9 row-1.  Expected value: 5.6359184e-01 to 1e-8 absolute.
// Anchor is given in plan §4.1.
// ----------------------------------------------------------------------------
void TestPsiIniAnchor()
{
   std::cout << "\n[T_TPV104_P_1] ψ_ini anchor\n";

   const real_t psi_expected = 5.6359184e-01;
   const real_t psi_computed = ComputeInitialPsiTPV104(TPV104Params::a_in);

   // Plan §4.1 specifies tolerance 1e-8 absolute.  The exact value has
   // 7 significant figures (5.6359184), so this is a very tight match.
   TEST_NEAR(psi_computed, psi_expected, 1e-8,
             "ψ_ini(a_in=0.01) = 5.6359184e-01 anchor");

   // R-009 (review 2026-04-24): verify the stable logsinh form does not
   // overflow at small a where the previous std::sinh(tau/(sigma_n*a))
   // form would return +inf.  At a = 0.001, arg = tau/(sigma_n*a)
   // = 4e7/(1.2e8*0.001) = 333.3 — sinh(333) overflows double; the
   // stable form must remain finite and non-NaN.
   const real_t psi_small_a = ComputeInitialPsiTPV104(0.001);
   TEST_ASSERT(std::isfinite(psi_small_a),
               "ComputeInitialPsiTPV104(a=0.001) finite (R-009 stability)");
   // Analytic expectation for large c:
   //   ψ = a · log(x · sinh(c)) ≈ a · (c + log(x/2))
   //   c = 333.333..., x = 2·V0/V_ini = 2e10, log(x/2) = log(1e10) = 23.026
   //   ψ ≈ 0.001 · (333.333 + 23.026) = 0.35636
   TEST_NEAR(psi_small_a, 0.356359184, 1e-6,
             "ψ_ini(a=0.001) = 0.356359... (large-c asymptotic)");
}

// ----------------------------------------------------------------------------
// T_TPV104_P_2 — derived-scalar static asserts (a_out, V_w_out, lambda, eta_s).
// These are compile-time `static_assert`s in the header; this test exists to
// anchor the invariants at runtime so a maintainer sees a visible pass line.
// ----------------------------------------------------------------------------
void TestDerivedScalars()
{
   std::cout << "\n[T_TPV104_P_2] derived scalars\n";

   TEST_NEAR(TPV104Params::a_out,
             TPV104Params::a_in + TPV104Params::da, 1e-15,
             "a_out == a_in + da");
   TEST_NEAR(TPV104Params::V_w_out,
             TPV104Params::V_w_in + TPV104Params::dV_w, 1e-15,
             "V_w_out == V_w_in + dV_w");
   TEST_NEAR(TPV104Params::lambda,
             TPV104Params::rho * TPV104Params::cp * TPV104Params::cp
             - 2.0 * TPV104Params::mu,
             1e-4,  // scale O(1e10)
             "lambda == rho*cp^2 - 2*mu");
   TEST_NEAR(TPV104Params::eta_s,
             TPV104Params::rho * TPV104Params::cs / 2.0, 1e-6,
             "eta_s == rho*cs/2");

   // Numerical anchors: lambda and mu for Poisson ratio = 0.25 should
   // give lambda == mu.  With cp = 6000, cs = 3464, ν ≈ 0.249 — very
   // close to but not exactly 0.25, so lambda ≈ mu but not equal.
   const real_t lambda_expected_min = 3.1e10;
   const real_t lambda_expected_max = 3.3e10;
   TEST_ASSERT(TPV104Params::lambda > lambda_expected_min &&
               TPV104Params::lambda < lambda_expected_max,
               "lambda in expected range [3.1e10, 3.3e10]");
}

// ----------------------------------------------------------------------------
// T_TPV104_P_3 — ComputeA_TPV104 spatial distribution.
//   (0, 7.5e3)      -> a_in       (hypocenter, centre of VW core)
//   (20e3, 20e3)    -> a_out      (far outside)
//   |x| = 13.5 km, z = 7.5 km — inside transition
// ----------------------------------------------------------------------------
void TestComputeA()
{
   std::cout << "\n[T_TPV104_P_3] ComputeA_TPV104 spatial distribution\n";

   // Centre of VW zone: a = a_in.
   real_t a_centre = ComputeA_TPV104(0.0, TPV104Params::hypo_down_dip);
   TEST_NEAR(a_centre, TPV104Params::a_in, 1e-15,
             "a at hypocenter = a_in");

   // Far outside: a = a_out.  (20 km, 20 km) lies outside both strike
   // and dip boxcars (W=15 km half-width in both directions, plus w=3 km
   // transition — so |x|>18 km is strictly outside strike boxcar;
   // |z-7.5|>10.5 km is strictly outside dip boxcar).
   real_t a_far = ComputeA_TPV104(20.0e3, 20.0e3);
   TEST_NEAR(a_far, TPV104Params::a_out, 1e-15,
             "a at (20 km, 20 km) = a_out");

   // Strike-side transition: the boxcar in strike has half-width L_s=15 km
   // and transition w_s=3 km, so a point at x=13.5 km, z=7.5 km is still
   // inside the strike boxcar (|13.5| ≤ 15).  The dip boxcar at
   // z=7.5 km (= hypo) is centre of VW (B_dip=1).  So this point is in
   // the VW core and a = a_in.  Use a point outside VW on the strike axis
   // instead: (16.5e3, 7.5e3) is in the strike-transition zone, inside
   // dip VW.  We expect a strictly between a_in and a_out.
   real_t a_transition = ComputeA_TPV104(16.5e3, TPV104Params::hypo_down_dip);
   TEST_ASSERT(a_transition > TPV104Params::a_in &&
               a_transition < TPV104Params::a_out,
               "a at (16.5 km, 7.5 km) ∈ (a_in, a_out) strictly");

   // Continuity: ComputeA is C∞ across the boundary — sample on both
   // sides and confirm the transition is monotone.
   real_t a_edge_in   = ComputeA_TPV104(14.99e3, TPV104Params::hypo_down_dip);
   real_t a_edge_out  = ComputeA_TPV104(15.01e3, TPV104Params::hypo_down_dip);
   TEST_ASSERT(a_edge_in <= a_edge_out,
               "a increases monotonically across the VW edge");
}

// ----------------------------------------------------------------------------
// T_TPV104_P_4 — ComputeVw_TPV104 returns V_w_in in VW core, V_w_out outside.
// ----------------------------------------------------------------------------
void TestComputeVw()
{
   std::cout << "\n[T_TPV104_P_4] ComputeVw_TPV104 spatial distribution\n";

   real_t vw_core = ComputeVw_TPV104(0.0, TPV104Params::hypo_down_dip);
   TEST_NEAR(vw_core, TPV104Params::V_w_in, 1e-15,
             "V_w at hypocenter (VW core) = 0.1 m/s");

   real_t vw_far = ComputeVw_TPV104(20.0e3, 20.0e3);
   TEST_NEAR(vw_far, TPV104Params::V_w_out, 1e-15,
             "V_w at (20 km, 20 km) (strengthening) = 1.0 m/s");

   // Intermediate: inside strike boxcar but in dip transition.
   real_t vw_trans = ComputeVw_TPV104(0.0, 17.0e3);
   TEST_ASSERT(vw_trans > TPV104Params::V_w_in &&
               vw_trans < TPV104Params::V_w_out,
               "V_w in dip transition ∈ (V_w_in, V_w_out) strictly");
}

// ----------------------------------------------------------------------------
// T_TPV104_P_5 — station list completeness.
// The nine SCEC TPV104 stations are (x2, x3) in km ∈
//   {0, ±9, ±12} × {3, 7.5, 12}
// but only the nine pairs listed in the SCEC benchmark-trace naming.
// ----------------------------------------------------------------------------
void TestStations()
{
   std::cout << "\n[T_TPV104_P_5] station list completeness\n";

   constexpr int N = 9;
   TEST_ASSERT(sizeof(kStationsTPV104)/sizeof(kStationsTPV104[0]) == N,
               "kStationsTPV104 has 9 entries");

   // Each label must be non-null and start with "x2_".
   for (int i = 0; i < N; ++i)
   {
      TEST_ASSERT(kStationsTPV104[i].label != nullptr,
                  "station label not null");
      std::string s(kStationsTPV104[i].label);
      TEST_ASSERT(s.rfind("x2_", 0) == 0,
                  "station label starts with 'x2_'");
   }

   // Spot-check a few coordinates.  Hypocenter station is at (0, 7.5 km).
   bool found_hypo = false;
   for (int i = 0; i < N; ++i)
   {
      if (kStationsTPV104[i].x2 == 0.0 &&
          kStationsTPV104[i].x3 == TPV104Params::hypo_down_dip)
      {
         found_hypo = true;
         break;
      }
   }
   TEST_ASSERT(found_hypo,
               "hypocenter station (0, 7.5 km) is in the list");

   // Ensure (12 km, 3 km) is present — reference station for Phase 3 probes.
   bool found_12_3 = false;
   for (int i = 0; i < N; ++i)
   {
      if (kStationsTPV104[i].x2 == 12.0e3 &&
          kStationsTPV104[i].x3 == 3.0e3)
      {
         found_12_3 = true;
         break;
      }
   }
   TEST_ASSERT(found_12_3,
               "station (12 km, 3 km) is in the list");
}

// ----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   TestPsiIniAnchor();
   TestDerivedScalars();
   TestComputeA();
   TestComputeVw();
   TestStations();

   TEST_PRINT_RESULTS();

   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
