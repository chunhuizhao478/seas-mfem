// D2 + D3 diagnostic — TPV104 σ_n invariance + mirror-symmetry.
//
// Motivation: the t12 dev run (results_t12_dev_job7677022) shows σ_n
// perturbations of 1-4 MPa at all 9 stations during/after rupture, with
// SIGN-FLIPPED magnitudes between mirror-symmetric station pairs (e.g.
// x2_-12_x3_3 = +3.7 MPa vs x2_+12_x3_3 = -4.2 MPa).  SeisSol shows σ_n
// flat at 120 MPa at every station.  This test set narrows the source
// by isolating ComputeTrialTraction in two regimes:
//
//   D2 (static-leak):  with bulk Q = 0 on both sides and TPV104-init
//                      DOFData (sigma_n0=120 MPa, tau2_0=40 MPa, ...),
//                      the trial traction MUST be exactly (0, 0, 0).
//                      Any non-zero output is a static frame-rotation
//                      bias / pre-stress contamination.  Such a bias
//                      is RULED OUT by the t12 data (σ_n exact at 120
//                      MPa pre-rupture); this test confirms the
//                      arithmetic agrees with that empirical finding.
//
//   D3 (mirror-symmetry):  for any homogeneous bulk-Q jump configuration
//                          on TPV104-init DOFData, the result must be
//                          deterministic — running with the same
//                          (Q_plus, Q_minus) inputs must produce the
//                          same outputs.  More tellingly: under a sign
//                          flip of the velocity-jump and a sign flip of
//                          the shear-stress jump (i.e. mirroring the
//                          off-fault wave field across x2=0), σ_n_trial
//                          must be unchanged because σ_n is even under
//                          x2-reflection while τ_2 is odd.
//
// The bug we suspect (per the σ_n analysis) is in the per-QP
// `qpd.sign_flipped` value flowing through `wave_operator.inl:2117-2131`,
// which determines which mesh element is the "+" side at each fault QP.
// That bug is at the wave-operator level, not the FaultFaceFlux level —
// so D2 and D3 alone CANNOT trip it.  But D2 + D3 establish the local
// "no leak" baseline at the FaultFaceFlux level: if D2/D3 PASS, the
// bug is upstream (rotation / element-side assignment).  If D2/D3
// FAIL, there's an additional bug downstream as well.
//
// Built locally; no MPI, no mesh, no time loop.  Per CLAUDE.md
// "feedback_no_local_reproducer.md", this is a unit test on a single
// QP fixture, not a TPV104 simulation.

#include "mfem.hpp"
#include "../../config/tpv104_params.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) <= t_) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg \
                << " (got " << v_ << ", expected " << e_ << ")\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << " (got " << v_ << ", expected " << e_ \
                << ", diff " << std::abs(v_ - e_) \
                << ", tol " << t_ << ")\n"; \
   } \
} while(0)

static const real_t RHO = TPV104Params::rho;
static const real_t CP  = TPV104Params::cp;
static const real_t CS  = TPV104Params::cs;

/// DOFData populated as `InitializeFaultDOFs_TPV104` would set it for the
/// hypocenter QP at TPV104 init: σ_n0 = 120 MPa, τ2_0 = 40 MPa (strike
/// pre-stress), τ1_0 = 0 (no dip pre-stress), nucleation channels = 0.
static DOFData MakeTPV104InitDOF()
{
   DOFData d;
   d.Zp_plus  = TPV104Params::Zp;  d.Zp_minus = TPV104Params::Zp;
   d.Zs_plus  = TPV104Params::Zs;  d.Zs_minus = TPV104Params::Zs;
   d.eta_p    = TPV104Params::eta_p;
   d.eta_s    = TPV104Params::eta_s;
   d.sigma_n0 = TPV104Params::sigma_n;
   d.tau1_0   = 0.0;
   d.tau2_0   = TPV104Params::tau_ini;
   d.sigma_n_nuc = 0.0;
   d.tau1_nuc    = 0.0;
   d.tau2_nuc    = 0.0;
   d.a   = TPV104Params::a_in;
   d.Dc  = TPV104Params::L;
   d.psi = ComputeInitialPsiTPV104(TPV104Params::a_in);
   d.slip_rate = TPV104Params::V_ini;
   d.V1 = 0.0;  d.V2 = TPV104Params::V_ini;
   d.slip1 = 0.0; d.slip2 = 0.0;
   d.tau1_corr = 0.0; d.tau2_corr = TPV104Params::tau_ini;
   d.sigma_n_corr = TPV104Params::sigma_n;
   return d;
}

// =====================================================================
// D2 — static invariance: bulk Q = 0 → trial = 0
// =====================================================================
//
// In MFEM's fluctuation-Q dispatch (the production R7-001(b) path,
// drivers/tpv104_driver.cpp:597-598 sets Q ≡ 0 at init), bulk Q
// represents the wave fluctuation only.  At t = 0 the fluctuation is
// zero everywhere.  ComputeTrialTraction reads bulk Q; with Q ≡ 0 it
// MUST return (0, 0, 0) — pre-stress lives in DOFData.sigma_n0/
// tau1_0/tau2_0 and is summed in by `CompleteFromTrial`, NOT here.
//
// If this test FAILS, there's a static frame-rotation bias contaminating
// σ_n_trial from somewhere — either (a) ComputeTrialTraction is reading
// pre-stress fields it shouldn't, or (b) the impedance values have been
// corrupted in init.  Either way, NOT what the t12 data showed (σ_n
// was exact at 120 MPa pre-rupture), so this test should PASS, providing
// the "no static leak" baseline against which D3 is interpreted.
void TestD2_StaticInvariance()
{
   std::cout << "Test D2: static invariance — bulk Q = 0 → trial = 0\n";

   DOFData d = MakeTPV104InitDOF();

   real_t Q_plus[NUM_STATE]  = {};
   real_t Q_minus[NUM_STATE] = {};

   real_t sn_trial = std::nan("0"), t1_trial = std::nan("0"),
          t2_trial = std::nan("0");
   FaultFaceFlux::ComputeTrialTraction(d, Q_plus, Q_minus,
                                       sn_trial, t1_trial, t2_trial);

   TEST_NEAR(sn_trial, 0.0, 1e-12,
             "D2.a sigma_n_trial == 0 with Q ≡ 0");
   TEST_NEAR(t1_trial, 0.0, 1e-12,
             "D2.b tau1_trial    == 0 with Q ≡ 0");
   TEST_NEAR(t2_trial, 0.0, 1e-12,
             "D2.c tau2_trial    == 0 with Q ≡ 0");
}

// =====================================================================
// D3a — mirror-symmetry of σ_n under tangential-velocity sign flip
// =====================================================================
//
// Take a homogeneous bulk-Q configuration with a non-zero tangential
// velocity jump (mode-II strike-slip pulse signature):
//   Q_plus = (Vstrike, 0, 0,  ..., 0)    in fault-local indexing
//   Q_minus= (-Vstrike, 0, 0, ..., 0)    in fault-local indexing
// ComputeTrialTraction reads Q[VX] = fault-normal velocity (==0 here)
// and Q[VY] = fault-tangent1 velocity, Q[VZ] = fault-tangent2 velocity
// (per fault_face_flux.cpp:53-65).  In MFEM convention t1=dip, t2=strike,
// so a strike-direction velocity jump puts (Vstrike, -Vstrike) into
// Q[VZ] (== t2-velocity).
//
// The trial is:
//   sigma_n_trial = etaP * (Q_minus[VX] - Q_plus[VX] + ...) = 0
//   tau1_trial    = etaS * (Q_minus[VY] - Q_plus[VY] + ...) = 0
//   tau2_trial    = etaS * (Q_minus[VZ] - Q_plus[VZ] + ...) = etaS * (-2 Vstrike)
// (with bulk SXX=SXY=SXZ=0; only the velocity jump enters)
//
// Now flip the sign of Vstrike (mirror across x2=0): Q'_plus = -Q_plus,
// Q'_minus = -Q_minus on the t2-velocity component.  By symmetry σ_n
// is even (unchanged), τ_2 is odd (sign-flips), τ_1 is unchanged
// (also zero anyway).
void TestD3a_MirrorSymmetryStrikeJump()
{
   std::cout << "Test D3a: mirror-symmetry of strike-tangential velocity jump\n";

   DOFData d = MakeTPV104InitDOF();
   const real_t Vstrike = 5.0;  // m/s, peak-rupture-scale magnitude

   // Forward configuration: Q_plus tangent2-velocity = +Vstrike
   real_t Q_plus_fwd[NUM_STATE]  = {};
   real_t Q_minus_fwd[NUM_STATE] = {};
   Q_plus_fwd [VZ] = +Vstrike;
   Q_minus_fwd[VZ] = -Vstrike;

   real_t sn_fwd, t1_fwd, t2_fwd;
   FaultFaceFlux::ComputeTrialTraction(d, Q_plus_fwd, Q_minus_fwd,
                                       sn_fwd, t1_fwd, t2_fwd);

   // Mirror configuration: flip the sign of the t2-velocity on both sides
   real_t Q_plus_mir[NUM_STATE]  = {};
   real_t Q_minus_mir[NUM_STATE] = {};
   Q_plus_mir [VZ] = -Vstrike;
   Q_minus_mir[VZ] = +Vstrike;

   real_t sn_mir, t1_mir, t2_mir;
   FaultFaceFlux::ComputeTrialTraction(d, Q_plus_mir, Q_minus_mir,
                                       sn_mir, t1_mir, t2_mir);

   // σ_n is even under x2 mirroring → unchanged
   TEST_NEAR(sn_fwd, sn_mir, 1e-12,
             "D3a.a sigma_n_trial mirror-invariant under strike-velocity flip");
   // τ_1 is identically zero in both → trivially equal
   TEST_NEAR(t1_fwd, t1_mir, 1e-12,
             "D3a.b tau1_trial unchanged under strike-velocity flip");
   // τ_2 is odd under x2 mirroring → sign-flipped
   TEST_NEAR(t2_fwd, -t2_mir, 1e-12,
             "D3a.c tau2_trial sign-flips under strike-velocity flip");

   // Sanity: σ_n_trial == 0 in this purely-tangential configuration
   TEST_NEAR(sn_fwd, 0.0, 1e-12,
             "D3a.d sigma_n_trial == 0 in pure-strike velocity-jump config");
}

// =====================================================================
// D3b — mirror-symmetry of σ_n under tangential-stress sign flip
// =====================================================================
//
// The other half of the strike-slip rupture-front signature is a
// SXY/SXZ stress jump.  Under x2 mirror, σ_n (== SYY in fault-local) is
// even, σ_xz (== shear in t2) is odd.  ComputeTrialTraction with
// mirrored Q[SXZ] should sign-flip τ_2 and leave σ_n unchanged.
void TestD3b_MirrorSymmetryStrikeStressJump()
{
   std::cout << "Test D3b: mirror-symmetry of strike-shear stress jump\n";

   DOFData d = MakeTPV104InitDOF();
   const real_t TauStrike = 10e6;  // 10 MPa shear-stress jump magnitude

   // Forward
   real_t Q_plus_fwd[NUM_STATE]  = {};
   real_t Q_minus_fwd[NUM_STATE] = {};
   Q_plus_fwd [SXZ] = +TauStrike;
   Q_minus_fwd[SXZ] = -TauStrike;

   real_t sn_fwd, t1_fwd, t2_fwd;
   FaultFaceFlux::ComputeTrialTraction(d, Q_plus_fwd, Q_minus_fwd,
                                       sn_fwd, t1_fwd, t2_fwd);

   // Mirror
   real_t Q_plus_mir[NUM_STATE]  = {};
   real_t Q_minus_mir[NUM_STATE] = {};
   Q_plus_mir [SXZ] = -TauStrike;
   Q_minus_mir[SXZ] = +TauStrike;

   real_t sn_mir, t1_mir, t2_mir;
   FaultFaceFlux::ComputeTrialTraction(d, Q_plus_mir, Q_minus_mir,
                                       sn_mir, t1_mir, t2_mir);

   TEST_NEAR(sn_fwd, sn_mir, 1e-12,
             "D3b.a sigma_n_trial mirror-invariant under strike-stress flip");
   TEST_NEAR(t1_fwd, t1_mir, 1e-12,
             "D3b.b tau1_trial unchanged under strike-stress flip");
   TEST_NEAR(t2_fwd, -t2_mir, 1e-12,
             "D3b.c tau2_trial sign-flips under strike-stress flip");
   TEST_NEAR(sn_fwd, 0.0, 1e-12,
             "D3b.d sigma_n_trial == 0 in pure-strike stress-jump config");
}

// =====================================================================
// D3c — coupled velocity + stress strike-slip pulse
// =====================================================================
//
// A real rupture-front transit at a fault QP couples both the velocity
// and stress jumps coherently (via the wave operator's Riemann
// imposed-state construction).  Test that σ_n_trial remains exactly
// zero for the combined rupture-front signature.  Any non-zero result
// here would mean ComputeTrialTraction has a strike→normal coupling
// in its arithmetic — which would be a static bug visible at every
// step of the simulation, NOT something we'd attribute to the per-QP
// sign-flip mismatch in wave_operator.inl.
void TestD3c_PureStrikeRuptureFront()
{
   std::cout << "Test D3c: pure-strike rupture-front signature, σ_n_trial == 0\n";

   DOFData d = MakeTPV104InitDOF();
   const real_t Vstrike   = 5.0;
   const real_t TauStrike = 10e6;

   real_t Q_plus[NUM_STATE]  = {};
   real_t Q_minus[NUM_STATE] = {};

   // Strike velocity (t2 = z) and strike stress (SXZ).  Both anti-symmetric
   // across the fault as a wave radiating in opposite directions on the
   // two sides.
   Q_plus [VZ]  = +Vstrike;
   Q_minus[VZ]  = -Vstrike;
   Q_plus [SXZ] = +TauStrike;
   Q_minus[SXZ] = -TauStrike;

   real_t sn_trial, t1_trial, t2_trial;
   FaultFaceFlux::ComputeTrialTraction(d, Q_plus, Q_minus,
                                       sn_trial, t1_trial, t2_trial);

   TEST_NEAR(sn_trial, 0.0, 1e-12,
             "D3c.a sigma_n_trial == 0 with pure-strike (V + tau) jump");
   TEST_NEAR(t1_trial, 0.0, 1e-12,
             "D3c.b tau1_trial == 0 with pure-strike (V + tau) jump");
   // τ_2 picks up both contributions; the magnitude is the formula
   // (etaS*(Q_minus[VZ] - Q_plus[VZ] + Q_plus[SXZ]/Zs + Q_minus[SXZ]/Zs))
   //   = etaS*(-2*Vstrike + 0)             // shear jump cancels for symmetric Q
   //   = (Zs/2) * (-2*Vstrike) = -Zs*Vstrike
   const real_t expected_t2 = -TPV104Params::Zs * Vstrike;
   TEST_NEAR(t2_trial, expected_t2, 1.0,
             "D3c.c tau2_trial = -Zs*Vstrike for symmetric stress jump");
}

// =====================================================================
// D3d — fault-NORMAL velocity jump should be the ONLY way to get σ_n != 0
// =====================================================================
//
// A non-zero σ_n_trial requires a fault-normal velocity jump (Q[VX]
// in fault-local) or a fault-normal stress jump (Q[SXX] in fault-local).
// This test verifies that — and inversely confirms that the σ_n channel
// is NOT cross-coupled to the tangential channels.
void TestD3d_NormalChannelIsolation()
{
   std::cout << "Test D3d: σ_n_trial responds to fault-normal jumps only\n";

   DOFData d = MakeTPV104InitDOF();

   // Pure fault-normal velocity jump → σ_n_trial != 0; τ_1, τ_2 = 0
   {
      real_t Q_plus[NUM_STATE]  = {};
      real_t Q_minus[NUM_STATE] = {};
      Q_plus [VX] = +0.01;
      Q_minus[VX] = -0.01;

      real_t sn, t1, t2;
      FaultFaceFlux::ComputeTrialTraction(d, Q_plus, Q_minus, sn, t1, t2);

      // sigma_n_trial = etaP * (Q_minus[VX] - Q_plus[VX]) = (Zp/2) * (-0.02)
      const real_t expected_sn = -TPV104Params::Zp * 0.01;
      TEST_NEAR(sn, expected_sn, 1.0,
                "D3d.a sigma_n responds linearly to fault-normal V jump");
      TEST_NEAR(t1, 0.0, 1e-12, "D3d.b tau1 == 0 for pure-normal V jump");
      TEST_NEAR(t2, 0.0, 1e-12, "D3d.c tau2 == 0 for pure-normal V jump");
   }

   // Pure fault-normal stress jump → σ_n_trial != 0; τ_1, τ_2 = 0
   {
      real_t Q_plus[NUM_STATE]  = {};
      real_t Q_minus[NUM_STATE] = {};
      Q_plus [SXX] = +1e6;
      Q_minus[SXX] = +1e6;  // both sides feel the same compression bump

      real_t sn, t1, t2;
      FaultFaceFlux::ComputeTrialTraction(d, Q_plus, Q_minus, sn, t1, t2);

      // For homogeneous (Zp_plus = Zp_minus = Zp), eta_p = Zp/2:
      // sigma_n_trial = (Zp/2) * (1e6/Zp + 1e6/Zp) = (Zp/2) * (2e6/Zp) = 1e6
      TEST_NEAR(sn, 1e6, 1.0,
                "D3d.d sigma_n trial = stress when both sides agree");
      TEST_NEAR(t1, 0.0, 1e-12, "D3d.e tau1 == 0 for pure-normal stress");
      TEST_NEAR(t2, 0.0, 1e-12, "D3d.f tau2 == 0 for pure-normal stress");
   }
}

int main()
{
   std::cout << "=== TPV104 σ_n invariance + mirror-symmetry tests ===\n";
   std::cout << "D2: static invariance at TPV104 init (bulk Q = 0)\n";
   std::cout << "D3: mirror-symmetry of trial under x2 → -x2\n\n";

   TestD2_StaticInvariance();
   TestD3a_MirrorSymmetryStrikeJump();
   TestD3b_MirrorSymmetryStrikeStressJump();
   TestD3c_PureStrikeRuptureFront();
   TestD3d_NormalChannelIsolation();

   std::cout << "\n=== Summary: " << num_passed << " passed, "
             << num_failed << " failed of " << num_tests << " ===\n";
   return num_failed == 0 ? 0 : 1;
}
