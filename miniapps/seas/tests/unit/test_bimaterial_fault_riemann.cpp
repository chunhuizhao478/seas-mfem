// Part B / Phase B0 (GATE): verify the per-side bi-material FAULT Riemann.
//
// The fault-flux math (FaultFaceFlux::ComputeTrialTraction / BuildImposedState)
// already consumes per-side impedances Zp_plus/Zp_minus, Zs_plus/Zs_minus.  This
// test proves it is the CORRECT bi-material welded Riemann for UNEQUAL impedance
// (Zp_plus != Zp_minus), in the LOCKED limit (V=0), so Part B can safely relax the
// homogeneity guards.  If this fails, B2 is BLOCKED.
//
// Frame (fault-local): SXX = sigma_n, SXY = tau1 (dip), SXZ = tau2 (strike);
// VX = v_n, VY = v_t1, VZ = v_t2.  Convention: side "minus" = side 1 (Z1),
// side "plus" = side 2 (Z2); n points to side 2 (plus); compression POSITIVE.
//
// Welded-interface star (derived; matches ComputeTrialTraction algebraically):
//   sigma* = (Z2*sigma_L + Z1*sigma_R + Z1*Z2*(v_L - v_R)) / (Z1 + Z2)
//   v*     = (Z1*v_L + Z2*v_R + (sigma_L - sigma_R)) / (Z1 + Z2)   [continuous]
// For a right-going incident P-wave (sigma_i = +Z1*v_i) from side 1, side 2 at
// rest, this gives R_sigma=(Z2-Z1)/(Z1+Z2) and T_v=2*Z1/(Z1+Z2).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/fault_face_flux.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_NEAR(v, ref, rel, msg) do {                                       \
   num_tests++;                                                                \
   const double vv = (v), rr = (ref);                                          \
   const double denom = std::max(1.0, std::abs(rr));                           \
   const double e = std::abs(vv - rr) / denom;                                 \
   if (e <= (rel)) { num_passed++;                                             \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific         \
                << std::setprecision(6) << vv << ", ref " << rr                \
                << ", rel " << e << ")\n"; }                                   \
   else { num_failed++;                                                        \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "  (got " << vv \
                << ", ref " << rr << ", rel " << e << " > " << (rel) << ")\n"; }\
} while (0)

namespace
{
// Side impedances Zp = rho*cp, Zs = rho*cs.
struct Mat { real_t rho, cp, cs; };
real_t Zp(const Mat &m) { return m.rho * m.cp; }
real_t Zs(const Mat &m) { return m.rho * m.cs; }

// Fill a DOFData's per-side impedances (side1=minus, side2=plus).
DOFData MakeBimatDOF(const Mat &side1, const Mat &side2)
{
   DOFData d;
   d.Zp_minus = Zp(side1); d.Zp_plus = Zp(side2);
   d.Zs_minus = Zs(side1); d.Zs_plus = Zs(side2);
   d.eta_p = d.Zp_plus * d.Zp_minus / (d.Zp_plus + d.Zp_minus);
   d.eta_s = d.Zs_plus * d.Zs_minus / (d.Zs_plus + d.Zs_minus);
   d.sigma_n0 = 0.0; d.tau1_0 = 0.0; d.tau2_0 = 0.0;
   d.sigma_n_nuc = 0.0; d.tau1_nuc = 0.0; d.tau2_nuc = 0.0;
   return d;
}

real_t StarTraction(real_t Z1, real_t Z2, real_t s_L, real_t s_R,
                    real_t v_L, real_t v_R)
{
   return (Z2 * s_L + Z1 * s_R + Z1 * Z2 * (v_L - v_R)) / (Z1 + Z2);
}
} // anonymous namespace

int main()
{
   std::cout << "\n=== Part B / B0: per-side bi-material FAULT Riemann ===\n";

   // TPV6 contrast: far side (1) vp/vs/rho = 3750/2165/2225 ; near side (2) = 6000/3464/2670.
   const Mat side1 = {2225.0, 3750.0, 2165.0};
   const Mat side2 = {2670.0, 6000.0, 3464.0};
   const real_t Z1p = Zp(side1), Z2p = Zp(side2);
   const real_t Z1s = Zs(side1), Z2s = Zs(side2);
   FaultFaceFlux ff(side2.rho, side2.cp, side2.cs);   // ctor material unused by the
                                                      // per-side trial/imposed path

   // -----------------------------------------------------------------------
   std::cout << "\n-- Test 1: trial traction == analytic welded star (general state) --\n";
   {
      DOFData d = MakeBimatDOF(side1, side2);
      real_t Qm[NUM_STATE] = {0}, Qp[NUM_STATE] = {0};
      // Arbitrary non-trivial bulk states on each side.
      Qm[SXX] = 7.0e6; Qm[SXY] = -2.0e6; Qm[SXZ] = 1.5e6;
      Qm[VX]  = 0.8;   Qm[VY]  = -0.3;   Qm[VZ]  = 0.5;
      Qp[SXX] = -3.0e6; Qp[SXY] = 1.0e6; Qp[SXZ] = -2.5e6;
      Qp[VX]  = -0.4;   Qp[VY]  = 0.6;   Qp[VZ]  = -0.2;

      real_t sn_tr, t1_tr, t2_tr;
      FaultFaceFlux::ComputeTrialTraction(d, Qp, Qm, sn_tr, t1_tr, t2_tr);

      const real_t sn_star = StarTraction(Z1p, Z2p, Qm[SXX], Qp[SXX], Qm[VX], Qp[VX]);
      const real_t t1_star = StarTraction(Z1s, Z2s, Qm[SXY], Qp[SXY], Qm[VY], Qp[VY]);
      const real_t t2_star = StarTraction(Z1s, Z2s, Qm[SXZ], Qp[SXZ], Qm[VZ], Qp[VZ]);
      TEST_NEAR(sn_tr, sn_star, 1.0e-12, "sigma_n_trial == welded normal star (P, per-side Zp)");
      TEST_NEAR(t1_tr, t1_star, 1.0e-12, "tau1_trial == welded shear star   (S, per-side Zs)");
      TEST_NEAR(t2_tr, t2_star, 1.0e-12, "tau2_trial == welded shear star   (S, per-side Zs)");
   }

   // -----------------------------------------------------------------------
   std::cout << "\n-- Test 2: incident right-going P-wave => R_sigma, T_v (locked) --\n";
   real_t Tv_check = 0.0;
   {
      DOFData d = MakeBimatDOF(side1, side2);
      const real_t v_i = 1.0;
      real_t Qm[NUM_STATE] = {0}, Qp[NUM_STATE] = {0};
      Qm[VX]  = v_i;          // incident velocity (side 1)
      Qm[SXX] = Z1p * v_i;    // right-going invariant: sigma_i = +Z1 v_i
      // Qp = 0 (side 2 at rest)

      real_t sn_tr, t1_tr, t2_tr;
      FaultFaceFlux::ComputeTrialTraction(d, Qp, Qm, sn_tr, t1_tr, t2_tr);

      // Interface (transmitted) traction; reflection coeff.
      const real_t sigma_star = 2.0 * Z1p * Z2p / (Z1p + Z2p) * v_i;
      TEST_NEAR(sn_tr, sigma_star, 1.0e-12,
                "interface normal traction == 2*Z1*Z2/(Z1+Z2)*v_i (transmitted)");
      const real_t R_sigma = (sn_tr - Z1p * v_i) / (Z1p * v_i);
      TEST_NEAR(R_sigma, (Z2p - Z1p) / (Z1p + Z2p), 1.0e-12,
                "stress reflection R_sigma == (Z2-Z1)/(Z1+Z2)");

      // Locked imposed state: V=0 in the NORMAL channel always; tau*_corr=0 here.
      EvalStageState s;
      s.sigma_n_corr = sn_tr;   // sigma_n0=0 => corr == trial (no friction in normal)
      s.tau1_corr = 0.0; s.tau2_corr = 0.0;
      real_t Qip[NUM_STATE], Qim[NUM_STATE];
      ff.BuildImposedState(d, s, Qp, Qm, Qip, Qim);

      const real_t Tv = 2.0 * Z1p / (Z1p + Z2p);   // velocity transmission
      Tv_check = Tv;
      TEST_NEAR(Qip[VX], Tv * v_i, 1.0e-9,
                "imposed +side normal velocity == T_v*v_i = 2*Z1/(Z1+Z2)*v_i");
      TEST_NEAR(Qim[VX], Qip[VX], 1.0e-9,
                "imposed normal velocity is CONTINUOUS across the locked fault "
                "(welded: v_imp_plus == v_imp_minus)");
      // Imposed traction is single-valued (continuous) on both sides.
      TEST_NEAR(Qip[SXX], Qim[SXX], 1.0e-12,
                "imposed normal traction is single-valued (Qip[SXX]==Qim[SXX])");
   }

   // -----------------------------------------------------------------------
   std::cout << "\n-- Test 3: homogeneous reduction (Z1==Z2) is the standard Riemann --\n";
   {
      DOFData d = MakeBimatDOF(side2, side2);   // both sides identical
      real_t Qm[NUM_STATE] = {0}, Qp[NUM_STATE] = {0};
      Qm[SXX] = 5.0e6; Qm[VX] = 0.7;
      Qp[SXX] = -2.0e6; Qp[VX] = -0.3;
      real_t sn_tr, t1_tr, t2_tr;
      FaultFaceFlux::ComputeTrialTraction(d, Qp, Qm, sn_tr, t1_tr, t2_tr);
      const real_t Z = Z2p;
      const real_t homog = 0.5 * (Qm[SXX] + Qp[SXX]) + 0.5 * Z * (Qm[VX] - Qp[VX]);
      TEST_NEAR(sn_tr, homog, 1.0e-12,
                "Z1==Z2: sigma_n_trial == (sigma_L+sigma_R)/2 + Z(v_L-v_R)/2 "
                "(homogeneous Riemann; byte-exact reduction)");
      EvalStageState s; s.sigma_n_corr = sn_tr;
      real_t Qip[NUM_STATE], Qim[NUM_STATE];
      ff.BuildImposedState(d, s, Qp, Qm, Qip, Qim);
      TEST_NEAR(Qim[VX], Qip[VX], 1.0e-12,
                "Z1==Z2: imposed normal velocity continuous (welded)");
   }

   // -----------------------------------------------------------------------
   // Test 4 (B2): with SetPerSideFluxApplied(true), an Evaluate variant runs on a
   // BIMATERIAL fault face WITHOUT aborting (the homogeneity guard is relaxed).  If
   // the relaxation were broken the guard's MFEM_ABORT would kill this process, so
   // reaching the post-call asserts IS the test.  (The complementary "guard still
   // aborts when the flag is NOT set" is covered by seas_test_fault_face_flux_
   // bimaterial_guard.)
   std::cout << "\n-- Test 4 (B2): per-side flag relaxes the bimaterial guard --\n";
   {
      FaultFaceFlux ff_relaxed(side2.rho, side2.cp, side2.cs);
      ff_relaxed.SetPerSideFluxApplied(true);
      TEST_NEAR(ff_relaxed.GetPerSideFluxApplied() ? 1.0 : 0.0, 1.0, 0.0,
                "SetPerSideFluxApplied(true) is reflected by the getter");

      DOFData d = MakeBimatDOF(side1, side2);       // Zp_plus != Zp_minus
      d.sigma_n0 = 120.0e6;                          // compressive background
      d.tau2_0   = 70.0e6;                           // strike pre-stress
      d.lsw_mu_s = 0.677; d.lsw_mu_d = 0.525; d.lsw_d_c = 0.40;
      d.slip1 = 0.0; d.slip2 = 0.0;

      real_t Qm[NUM_STATE] = {0}, Qp[NUM_STATE] = {0};
      Qm[VZ] = 0.2; Qm[SXZ] = 1.0e6;                 // some strike-channel bulk state
      real_t Qip[NUM_STATE], Qim[NUM_STATE];
      ff_relaxed.EvaluateLSW(d, Qp, Qm, Qip, Qim);   // MUST NOT abort (guard relaxed)

      bool finite = true;
      for (int c = 0; c < NUM_STATE; ++c)
      { finite = finite && std::isfinite(Qip[c]) && std::isfinite(Qim[c]); }
      TEST_NEAR(finite ? 1.0 : 0.0, 1.0, 0.0,
                "EvaluateLSW on a bimaterial fault (flag set) ran without abort and "
                "produced finite imposed states");
   }

   (void)Tv_check;
   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
