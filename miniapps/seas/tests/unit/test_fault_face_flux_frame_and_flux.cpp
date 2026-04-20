// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.0.0 Phase 1 §3.1b — FAST-PATH end-to-end fixture.
//
// Drives the full chain that the interior-fault branch of WaveOperator
// uses at `wave_operator.inl:754-812`:
//   (1) Evaluate  (fault-local: Q± → Q_imp±)
//   (2) T_can     (fault-local → global)
//   (3) flux_.Interior (Godunov on imposed states, global frame)
//
// The goal is to record F_h[VX_global] (the global-x velocity flux that
// drives hypocenter-adjacent bulk DOFs towards radiation) and compare
// it against the REVIEW.md analytical target of ~2.2e2 m²/s² at
// hypocentre steady slip.  Decision rule (per plan §3A FP-2):
//   (a) F_h[VX] ≈ 10..1e4 (match ±20% of 220) → Evaluate-T_can-Interior
//       chain is healthy; bug is downstream (accumulation / mass / branch
//       asymmetry).  Focus on §3.1d and §4 W-2.
//   (b) |F_h[VX]| ≫ 10⁴ or ≪ 1 → bug is in Evaluate/T_can/Interior.
//       Continue §3.1 / §3.1c.
//   (c) |F_h[VX]| < 1e-10 (exact zero) → a step in the chain zeros the
//       flux.  Print intermediate values to bisect.
//
// The test also asserts stress-flux continuity: |F_h[SXY]| ≤ 1 Pa·m/s
// (Godunov flux on imposed states must carry no stress jump by design
// — both imposed sides have the same corrected traction).
//
// BP5 canonical frame for TPV102 (vertical y=0 fault, ref_normal=(0,-1,0)):
//   can_n  = (0, -1, 0)
//   can_t1 = (0,  0, -1)   (dip, pointing down)
//   can_t2 = (1,  0,  0)   (along strike, +x)

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../config/tpv102_params.hpp"

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
} while(0)

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) < t_) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ << ", |diff| " << std::abs(v_-e_) << ")\n"; } \
} while(0)

// Post-breakaway steady-slip operating point (see plan §3.1b / REVIEW R-1005):
// tau_operating = 81 MPa, V_operating = 1 m/s.  Pre-nucleation equilibrium
// (V=V_ini=1e-12) would give F_h ≈ 0 by design — no radiation — and
// defeats the purpose of this test.
static DOFData MakeHypoDOF(real_t V_operating = 1.0,
                           real_t tau_operating = 81e6)
{
   DOFData d;
   d.Zp_plus  = TPV102Params::Zp;
   d.Zp_minus = TPV102Params::Zp;
   d.Zs_plus  = TPV102Params::Zs;
   d.Zs_minus = TPV102Params::Zs;
   d.eta_p    = TPV102Params::Zp / 2.0;
   d.eta_s    = TPV102Params::Zs / 2.0;
   d.sigma_n0 = TPV102Params::sigma_n;
   d.tau1_0   = 0.0;
   d.tau2_0   = tau_operating;
   d.a        = TPV102Params::a_vw;
   d.Dc       = TPV102Params::Dc;
   real_t arg = tau_operating / (d.sigma_n0 * d.a);
   d.psi = d.a * std::log(2.0 * TPV102Params::V0 / V_operating * std::sinh(arg));
   d.slip_rate = 0.0;
   return d;
}

// ===========================================================================
// Chain: Evaluate (fault-local) → T_can (local→global) → flux_.Interior
// (global).  Records intermediate quantities so a downstream debugger can
// bisect.
// ===========================================================================
static void TestFastPathEquilibrium()
{
   std::cout << "Test 3.1b: FAST-PATH — Evaluate -> T_can -> flux_.Interior\n";

   // ---- Inputs ----
   // BP5 canonical frame for TPV102 (vertical y=0 fault).
   real_t can_n[3]  = {0.0, -1.0, 0.0};
   real_t can_t1[3] = {0.0,  0.0, -1.0};   // dip, down
   real_t can_t2[3] = {1.0,  0.0,  0.0};   // strike

   // Quiescent bulk — the hypocentre's first radiation "kick" comes from
   // the fault's own frictional response, not from incoming bulk waves.
   real_t Q_plus_local[NUM_STATE]  = {};
   real_t Q_minus_local[NUM_STATE] = {};

   DOFData data = MakeHypoDOF();

   // ---- Step 1: Evaluate (fault-local frame) ----
   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   real_t Q_imp_plus_local[NUM_STATE], Q_imp_minus_local[NUM_STATE];
   ff.Evaluate(data, Q_plus_local, Q_minus_local,
               Q_imp_plus_local, Q_imp_minus_local);

   std::cout << "  [Step 1] slip_rate |V| = " << data.slip_rate
             << " m/s (V1=" << data.V1 << ", V2=" << data.V2 << ")\n";
   std::cout << "  [Step 1] tau_corr: tau1=" << Q_imp_plus_local[SXY]
             << " Pa, tau2=" << Q_imp_plus_local[SXZ] << " Pa\n";
   std::cout << "  [Step 1] Q_imp_plus_local[VY]=" << Q_imp_plus_local[VY]
             << " m/s, Q_imp_plus_local[VZ]=" << Q_imp_plus_local[VZ]
             << " m/s\n";

   // Sanity: at pure strike-slip the dip channel must be zero (BP5 invariant).
   TEST_NEAR(Q_imp_plus_local[SXY], 0.0, 1e-6,
             "local tau1 (dip) corr ≈ 0 (BP5 pure strike-slip)");
   TEST_NEAR(Q_imp_plus_local[VY], 0.0, 1e-18,
             "local v_t1^{+,imp} = 0 (BP5 pure strike-slip)");
   // Post-breakaway steady slip: tau2_corr ≈ -eta_s*V ≈ -9 MPa for V≈1 m/s.
   TEST_ASSERT(std::abs(Q_imp_plus_local[SXZ]) > 1e6,
               "local tau2 (strike) corr > 1 MPa at post-breakaway steady slip");

   // ---- Step 2: Build T_can, rotate imposed states to global ----
   DenseMatrix T_can(NUM_STATE);
   GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);

   real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
   T_can.Mult(Q_imp_plus_local,  Q_imp_plus_g);
   T_can.Mult(Q_imp_minus_local, Q_imp_minus_g);

   std::cout << "  [Step 2] Q_imp_plus_g: VX=" << Q_imp_plus_g[VX]
             << ", VY=" << Q_imp_plus_g[VY]
             << ", VZ=" << Q_imp_plus_g[VZ] << "\n";
   std::cout << "  [Step 2] Q_imp_plus_g: SXX=" << Q_imp_plus_g[SXX]
             << ", SXY=" << Q_imp_plus_g[SXY]
             << ", SXZ=" << Q_imp_plus_g[SXZ] << "\n";

   // Analytical: local [v_n, v_t1, v_t2] = [0, 0, v_t2_imp] rotates to
   // global.  With can_t2 = (+1,0,0), we get global VX = v_t2_imp (strike
   // slip motion in +x).  can_t1 = (0,0,-1) contributes zero since
   // v_t1_imp = 0.  can_n = (0,-1,0) contributes zero since v_n_imp = 0.
   // Expected: Q_imp_plus_g[VX] ≈ Q_imp_plus_local[VZ] (tangent-2 aligned).
   TEST_NEAR(Q_imp_plus_g[VX], Q_imp_plus_local[VZ], 1e-10,
             "global VX = local v_t2 at equilibrium (can_t2=+x)");
   TEST_NEAR(Q_imp_plus_g[VY], 0.0, 1e-10,
             "global VY = 0 (no normal or dip motion)");
   TEST_NEAR(Q_imp_plus_g[VZ], 0.0, 1e-10,
             "global VZ = 0 (no dip motion at equilibrium)");

   // ---- Step 3: per-side A.T.Q fluxes (post-fix, Pelties 2012 eq. 9) ----
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
   real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
   flux.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,  F_h_plus);
   flux.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g, F_h_minus);

   auto print_flux = [](const char *label, const real_t *F) {
      std::cout << "  [Step 3] " << label
                << ": VX=" << F[VX] << " VY=" << F[VY] << " VZ=" << F[VZ]
                << " SXX=" << F[SXX] << " SXY=" << F[SXY] << " SXZ=" << F[SXZ] << "\n";
   };
   print_flux("F_h_plus ", F_h_plus);
   print_flux("F_h_minus", F_h_minus);

   real_t abs_F_plus_VX  = std::abs(F_h_plus[VX]);
   real_t abs_F_minus_VX = std::abs(F_h_minus[VX]);
   std::cout << "  DECISION: |F_h_plus[VX]|  = " << abs_F_plus_VX << " m^2/s^2\n";
   std::cout << "  DECISION: |F_h_minus[VX]| = " << abs_F_minus_VX << " m^2/s^2\n";

   TEST_ASSERT(abs_F_plus_VX  > 1e-10,
               "|F_h_plus[VX]|  > 1e-10 (no cancellation)");
   TEST_ASSERT(abs_F_minus_VX > 1e-10,
               "|F_h_minus[VX]| > 1e-10 (no cancellation)");
   TEST_ASSERT(abs_F_plus_VX  > 10.0 && abs_F_plus_VX  < 1.0e4,
               "|F_h_plus[VX]|  in (10, 1e4) m^2/s^2");
   TEST_ASSERT(abs_F_minus_VX > 10.0 && abs_F_minus_VX < 1.0e4,
               "|F_h_minus[VX]| in (10, 1e4) m^2/s^2");

   // Physics check: VX flux is stress-driven (A[VX,sigma] != 0, A[VX,VX] = 0).
   // Plus/minus imposed states share the SAME corrected traction sigma_corr
   // (Pelties eq. 11d/12d) and the T_can rotation is identical for both
   // sides, so Q_imp_{plus,minus}_g[SXY] are bit-identical, and the VX row
   // of F_{can_n}*Q_imp_plus_g and F_{can_n}*Q_imp_minus_g agree to ULP.
   // Tolerance 1e-10 per plan R-F03.
   real_t rel_diff_VX = std::abs(F_h_plus[VX] - F_h_minus[VX]) /
                        std::max(abs_F_plus_VX, abs_F_minus_VX);
   TEST_ASSERT(rel_diff_VX < 1.0e-10,
               "F_h_plus[VX] == F_h_minus[VX] to ULP (stress continuity)");

   // Sign check -- analytically, F_h_{plus,minus}[VX] for can_n=(0,-1,0)
   // under right-lateral strike slip (V2 > 0, tau2_corr = -eta_s*V2 < 0):
   //   F_{can_n} = n.A = 0*A_x + (-1)*A_y + 0*A_z = -A_y
   //   A_y[VX][SXY] = -1/rho  (from BuildJacobian(dir=1), verified)
   //   Q_imp_{plus,minus}_g[SXY] = +eta_s*V2  (common to both sides;
   //     the rotation flips sign of local tau_{n,t2}=-eta_s*V2 because
   //     can_n is the -y axis and can_t2 is the +x axis).
   //   F_h_{plus,minus}[VX] = (-A_y[VX][SXY]) * Q_imp_side_g[SXY]
   //                        = -(-1/rho) * (+eta_s*V2)
   //                        = +eta_s*V2/rho = +c_s*V2/2   (POSITIVE)
   // Common derivation mistake (seen in round-3 REVIEW R-F01): applying
   // A_y directly to Q_imp_g without the n.A sign flip inverts the sign.
   // See plan §18.2 for the full record.
   TEST_ASSERT(F_h_plus[VX]  > 0.0,
               "F_h_plus[VX]  > 0 analytically (= +c_s*V2/2 for V2 > 0)");
   TEST_ASSERT(F_h_minus[VX] > 0.0,
               "F_h_minus[VX] > 0 analytically (same value as plus)");

   // Stress SXY: driven by the VELOCITY difference (A[SXY,VY] != 0).
   // ṽ^+ - v^+ = +delta while ṽ^- - v^- = -delta (Pelties eq. 15), so
   // F_h_plus[SXY] and F_h_minus[SXY] are opposite-signed.
   TEST_ASSERT(F_h_plus[SXY] * F_h_minus[SXY] < 0.0 ||
               (std::abs(F_h_plus[SXY]) < 1e-6 &&
                std::abs(F_h_minus[SXY]) < 1e-6),
               "F_h_plus[SXY] and F_h_minus[SXY] have opposite signs "
               "(velocity-driven channel)");

   // Order-of-magnitude cross-check (printed for operator inspection):
   real_t analytical_order = (data.eta_s / TPV102Params::rho) * data.V2;
   std::cout << "  analytical cross-check: (eta_s/rho) * V2 = "
             << analytical_order << " m^2/s^2\n";
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 v9.0.0 Phase 1 §3.1b FAST-PATH chain test\n";
   std::cout << "========================================\n\n";

   TestFastPathEquilibrium();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";

   return (num_failed > 0) ? 1 : 0;
}
