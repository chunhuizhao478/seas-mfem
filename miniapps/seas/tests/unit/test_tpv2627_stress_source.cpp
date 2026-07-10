// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_tpv2627_stress_source.cpp — Phase 1 of
// PLAN_TPV26_27_spatial_dyn_driver_2026-06-28.md §5.
//
// Golden tests for spatial::Tpv2627DepthStressSource (the SCEC TPV26/27
// depth-dependent initial stress, spec Part 3 / PLAN §1.2).  The source's
// Evaluate() returns a compression-POSITIVE Cauchy tensor in the code
// frame; these tests hand-project it onto the canonical fault basis
//   n  = (0, -1, 0)   (fault normal, ref_normal)
//   t1 = (0,  0, -1)  (down-dip)
//   t2 = (1,  0,  0)  (along-strike)
// exactly as FaultGeometry::ComputeParams does (fault_geometry_safs_
// templated.inl:83-127), then subtract the hydrostatic pore pressure
// P_p = P_p_grad * depth with P_p_grad = water_density*g = 9800 Pa/m to
// recover the effective normal stress the spec's "n-stress" reports.
//
// Standalone (no MPI, no mesh): Evaluate is a pure analytic function.
//
// NOTE on tolerances: the spec quotes sigma_n_eff = 175.64 MPa and
// tau_strike = 27.66 MPa to 4 significant figures (0.01 MPa precision).
// The exact analytic values from b11/b33/b13 are 175.6409 MPa and
// 27.6633 MPa, so the stress asserts use tol = 1e4 Pa (0.01 MPa), which
// matches the spec's quoted precision.  (The plan §5 sketch wrote 1e3;
// that is tighter than its own rounded golden value for tau_strike —
// diff 3.3e3 Pa — so 1e4 is used.  See the completion report.)

#include "mfem.hpp"

#include "../../spatial/code/spatial_stress.hpp"
#include "../../spatial/code/spatial_friction.hpp"

#include "test_macros.hpp"

#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas::spatial;

namespace
{

// Canonical fault-basis projection of a code-frame Cauchy tensor S,
// mirroring FaultGeometry::ComputeParams (no external sigma_xy flip).
struct Proj
{
   real_t sigma_n_total;   // compression-positive normal
   real_t tau_dip;         // t1 component
   real_t tau_strike;      // t2 component (right-lateral positive)
};

Proj project(const mfem::DenseMatrix& S)
{
   const real_t n[3]  = { 0.0, -1.0, 0.0 };
   const real_t t1[3] = { 0.0,  0.0, -1.0 };   // dip
   const real_t t2[3] = { 1.0,  0.0,  0.0 };   // strike
   real_t Sn[3];
   for (int r = 0; r < 3; ++r)
   {
      Sn[r] = S(r, 0) * n[0] + S(r, 1) * n[1] + S(r, 2) * n[2];
   }
   Proj p;
   p.sigma_n_total = n[0]  * Sn[0] + n[1]  * Sn[1] + n[2]  * Sn[2];
   p.tau_dip       = t1[0] * Sn[0] + t1[1] * Sn[1] + t1[2] * Sn[2];
   p.tau_strike    = t2[0] * Sn[0] + t2[1] * Sn[1] + t2[2] * Sn[2];
   return p;
}

// Source with the spec-default parameters (what a bare
// [stress] kind="tpv2627_depth" config produces).
Tpv2627DepthStressSource SpecSource()
{
   return Tpv2627DepthStressSource(/*rho*/ 2670.0, /*g*/ 9.8,
                                   /*water_density*/ 1000.0,
                                   /*b11*/ 0.926793, /*b33*/ 1.073206,
                                   /*b13*/ -0.169029,
                                   /*omega_top_m*/ 15000.0,
                                   /*omega_bot_m*/ 20000.0);
}

// Effective pore pressure ComputeParams would subtract (P_p_grad = 9800).
real_t Pf(real_t depth) { return 9800.0 * depth; }

}  // namespace

// Golden check at the hypocenter depth (10 km, Omega = 1): the on-fault
// shear/effective-normal ratio must equal the spec's 0.1575, strictly
// between mu_d = 0.12 and mu_s = 0.18.
void Test_Hypocenter_Ratio()
{
   auto src = SpecSource();
   const real_t depth = 10000.0;
   const mfem::DenseMatrix S = src.Evaluate(-5000.0, 0.0, -depth);
   const Proj p = project(S);
   const real_t sigma_n_eff = p.sigma_n_total - Pf(depth);

   TEST_NEAR(sigma_n_eff, 175.64e6, 1e4,
             "hypocenter sigma_n_eff = 175.64 MPa (compression-positive)");
   TEST_NEAR(p.tau_strike, 27.66e6, 1e4,
             "hypocenter tau_strike = 27.66 MPa (right-lateral positive)");
   const real_t ratio = p.tau_strike / sigma_n_eff;
   TEST_NEAR(ratio, 0.1575, 1e-3, "hypocenter tau/sigma_n_eff = 0.1575");
   TEST_ASSERT(0.12 < ratio && ratio < 0.18,
               "hypocenter ratio strictly between mu_d=0.12 and mu_s=0.18");
}

// The right-lateral sign must be baked into Evaluate (requirement 3):
// after the plain no-flip projection, tau_strike > 0.
void Test_RightLateral_Sign()
{
   auto src = SpecSource();
   const Proj p = project(src.Evaluate(-5000.0, 0.0, -10000.0));
   TEST_ASSERT(p.tau_strike > 0.0,
               "tau_strike right-lateral POSITIVE after no-flip projection");
}

// Omega taper: depths 15000 / 17500 / 20000 m -> Omega = 1 / 0.5 / 0.
// At Omega = 0 the tensor is isotropic (sigma11 = sigma33 = sigma22) with
// zero on-fault shear.
void Test_Omega_Taper()
{
   auto src = SpecSource();

   // Omega = 1 boundary (depth = omega_top_m): full shear.
   const Proj p15 = project(src.Evaluate(0.0, 0.0, -15000.0));
   TEST_NEAR(p15.tau_strike, 4.1495e7, 1e4,
             "Omega=1 tau_strike @ 15 km = 41.495 MPa");

   // Omega = 0.5 (midpoint): half-weighted shear.
   const Proj p175 = project(src.Evaluate(0.0, 0.0, -17500.0));
   TEST_NEAR(p175.tau_strike, 2.4205e7, 1e4,
             "Omega=0.5 tau_strike @ 17.5 km = 24.205 MPa");

   // Omega = 0 (depth = omega_bot_m): isotropic, no shear.
   const mfem::DenseMatrix S20 = src.Evaluate(0.0, 0.0, -20000.0);
   TEST_NEAR(S20(0, 1), 0.0, 1e-6, "Omega=0 on-fault shear -> 0");
   TEST_NEAR(S20(0, 0), S20(2, 2), 1.0,
             "Omega=0 sigma11 -> sigma22 (S_xx == S_zz)");
   TEST_NEAR(S20(1, 1), S20(2, 2), 1.0,
             "Omega=0 sigma33 -> sigma22 (S_yy == S_zz)");
}

// Surface-breaking top node (z = 0): depth = 0 => Pf = 0, sigma22 = 0, so
// the whole tensor is zero (finite, no NaN); sigma_n_total = 0.
void Test_Surface_ZeroPf()
{
   auto src = SpecSource();
   const mfem::DenseMatrix S0 = src.Evaluate(0.0, 0.0, 0.0);
   const Proj p0 = project(S0);
   TEST_NEAR(p0.sigma_n_total, 0.0, 1e-6, "surface (z=0) sigma_n_total = 0");
   TEST_ASSERT(std::isfinite(S0(0, 0)) && std::isfinite(S0(1, 1))
               && std::isfinite(S0(0, 1)),
               "surface tensor finite (no NaN)");
}

// R-001: the depth profile builds sigma11/sigma33/sigma13 from the EFFECTIVE
// vertical stress (sigma22 + Pf).  ComputeParams subtracts a SEPARATE pore
// pressure P_p = P_p_grad*depth.  If P_p_grad != water_density*g, the on-fault
// ratio silently leaves (mu_d, mu_s) and the fault can never slip.  This test
// documents WHY the parser now enforces the equality.
void Test_R001_Pf_Grad_Must_Match_WaterDensity_g()
{
   auto src = SpecSource();
   const real_t depth = 10000.0;
   const Proj p = project(src.Evaluate(-5000.0, 0.0, -depth));

   const real_t ratio_ok  = p.tau_strike / (p.sigma_n_total - 9800.0 * depth);
   const real_t ratio_bad = p.tau_strike / (p.sigma_n_total - 0.0 * depth);

   TEST_ASSERT(0.12 < ratio_ok && ratio_ok < 0.18,
               "P_p_grad = water_density*g -> ratio inside (mu_d, mu_s)");
   TEST_ASSERT(ratio_bad < 0.12,
               "P_p_grad = 0 -> ratio BELOW mu_d: the fault could never slip "
               "(this is what the [stress] tpv2627_depth parser guard rejects)");
}

// R-007: water_density = 0 is the sanctioned "no pore pressure" setting and
// must still produce a right-lateral-positive shear; negatives are rejected by
// the ctor (MFEM_VERIFY), which would otherwise flip the background sense.
void Test_R007_ZeroWaterDensity_Is_TotalStress()
{
   Tpv2627DepthStressSource src(2670.0, 9.8, /*water_density=*/0.0,
                                0.926793, 1.073206, -0.169029, 15000.0, 20000.0);
   const Proj p = project(src.Evaluate(0.0, 0.0, -10000.0));
   TEST_ASSERT(p.tau_strike > 0.0,
               "water_density = 0 still yields right-lateral-positive shear");
   TEST_ASSERT(p.sigma_n_total > 0.0,
               "water_density = 0 yields a compressive (total) normal stress");
}

// Pure strike-slip: sigma12 = sigma23 = 0 => no down-dip prestress.
void Test_Dip_Shear_Zero()
{
   auto src = SpecSource();
   const Proj p = project(src.Evaluate(-5000.0, 0.0, -10000.0));
   TEST_NEAR(p.tau_dip, 0.0, 1e-3, "dip prestress ~ 0 (pure strike-slip)");
}

int main(int /*argc*/, char * /*argv*/[])
{
   std::cout << "=== test_tpv2627_stress_source (Phase 1) ===\n";
   Test_Hypocenter_Ratio();
   Test_RightLateral_Sign();
   Test_Omega_Taper();
   Test_Surface_ZeroPf();
   Test_R001_Pf_Grad_Must_Match_WaterDensity_g();
   Test_R007_ZeroWaterDensity_Is_TotalStress();
   Test_Dip_Shear_Zero();
   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
