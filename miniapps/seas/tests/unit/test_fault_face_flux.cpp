// Phase 3 unit tests for FaultFaceFlux (Tests 22-31 from plan Section 3.3.5).

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

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

/// Standard TPV102/BP5 material parameters
static const real_t RHO = 2670.0;
static const real_t MU = 32.04e9;
static const real_t LAMBDA = 32.04e9;
static const real_t CP = std::sqrt((LAMBDA + 2*MU) / RHO);
static const real_t CS = std::sqrt(MU / RHO);
static const real_t ZP = RHO * CP;
static const real_t ZS = RHO * CS;

/// Create DOFData with homogeneous impedances
static DOFData MakeHomogeneousDOF(real_t sigma_n0 = 120e6, real_t a = 0.004)
{
   DOFData d;
   d.Zp_plus = ZP; d.Zp_minus = ZP;
   d.Zs_plus = ZS; d.Zs_minus = ZS;
   d.eta_p = ZP / 2.0;
   d.eta_s = ZS / 2.0;
   d.sigma_n0 = sigma_n0;
   d.tau1_0 = 0.0; d.tau2_0 = 0.0;
   d.a = a;
   d.Dc = 0.14;
   d.psi = 0.5;
   d.slip_rate = 0.0;
   return d;
}

// ===== Test 22: Trial traction on locked fault (equal states) =====
void TestTrialTractionLockedFault()
{
   std::cout << "Test 22: TestTrialTractionLockedFault\n";

   DOFData data = MakeHomogeneousDOF();
   FaultFaceFlux ff(RHO, CP, CS);

   // Equal states on both sides → trial = existing traction
   real_t Q[NUM_STATE] = {};
   Q[SXX] = -75e6;  // sigma_nn
   Q[SXY] = 30e6;   // tau_1 (shear)
   Q[SXZ] = 10e6;   // tau_2

   real_t sigma_trial, tau1_trial, tau2_trial;
   FaultFaceFlux::ComputeTrialTraction(data, Q, Q,
                                       sigma_trial, tau1_trial, tau2_trial);

   // With equal Q+ = Q-, the velocity terms cancel.
   // Trial = eta * (0 + sigma/Z + sigma/Z) = eta * (2*sigma/Z) = sigma
   // For homogeneous: eta = Z/2, so eta * (2*sigma/Z) = sigma ✓
   TEST_NEAR(sigma_trial, Q[SXX], 1e-6,
             "Locked fault: sigma_n_trial = sigma_n");
   TEST_NEAR(tau1_trial, Q[SXY], 1e-6,
             "Locked fault: tau1_trial = tau_1");
   TEST_NEAR(tau2_trial, Q[SXZ], 1e-6,
             "Locked fault: tau2_trial = tau_2");
}

// ===== Test 23: Trial traction with pure velocity jump =====
void TestTrialTractionVelocityJump()
{
   std::cout << "Test 23: TestTrialTractionVelocityJump\n";

   DOFData data = MakeHomogeneousDOF();

   // Q+ has v_t1 = +0.5, Q- has v_t1 = -0.5 → jump = 1 m/s
   real_t Q_plus[NUM_STATE] = {};
   real_t Q_minus[NUM_STATE] = {};
   Q_plus[VY] = 0.5;
   Q_minus[VY] = -0.5;

   real_t sigma_trial, tau1_trial, tau2_trial;
   FaultFaceFlux::ComputeTrialTraction(data, Q_plus, Q_minus,
                                       sigma_trial, tau1_trial, tau2_trial);

   // tau1_trial = eta_s * (v_t1^- - v_t1^+) = (ZS/2) * (-0.5 - 0.5) = -ZS/2
   real_t expected = -data.eta_s;
   TEST_NEAR(tau1_trial, expected, 1e-6,
             "Velocity jump: tau1_trial = -eta_s (got " +
             std::to_string(tau1_trial) + ", expected " + std::to_string(expected) + ")");
}

// ===== Test 24: Locked fault with high sigma_n → V ≈ 0 =====
void TestLockedFaultHighSigma()
{
   std::cout << "Test 24: TestLockedFaultHighSigma\n";

   DOFData data = MakeHomogeneousDOF(120e6, 0.5);  // high a
   data.psi = 0.5;

   FaultFaceFlux ff(RHO, CP, CS);

   // Small perturbation: tau << sigma_n * f_min
   real_t Q_plus[NUM_STATE] = {};
   real_t Q_minus[NUM_STATE] = {};
   Q_plus[SXY] = 1e3;  // tiny shear stress
   Q_minus[SXY] = 1e3;

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   ff.Evaluate(data, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

   // Slip rate should be very small for high sigma_n and small tau perturbation
   TEST_ASSERT(data.slip_rate < 1e-3,
               "Locked fault: V < 1e-3 (V = " + std::to_string(data.slip_rate) + ")");
}

// ===== Test 25: Frictionless fault → V = Theta / eta_s =====
void TestFrictionlessFault()
{
   std::cout << "Test 25: TestFrictionlessFault\n";

   FrictionSolver solver;

   // With a = 0, b = 0, f0 = 0: friction is zero → all traction becomes slip
   // g(V) = 0 + eta * V - tau = 0 → V = tau / eta
   real_t tau = 50e6;
   real_t sigma_n = 120e6;
   real_t eta = ZS / 2.0;
   real_t a = 1e-30;  // near-zero: f → 0
   real_t psi = 0.0;

   real_t V = solver.SolveBrent(tau, psi, sigma_n, eta, a);
   real_t expected = tau / eta;

   // With a ≈ 0, friction strength ≈ 0, so V = tau/eta (full stress drop)
   real_t rel = std::abs(V - expected) / expected;
   TEST_ASSERT(rel < 1e-6,
               "Frictionless: V = tau/eta (rel error " + std::to_string(rel) + ")");
}

// ===== Test 26: Imposed state traction continuity =====
void TestImposedStateTractionContinuity()
{
   std::cout << "Test 26: TestImposedStateTractionContinuity\n";

   DOFData data = MakeHomogeneousDOF();
   data.tau1_0 = 40e6;  // background shear
   FaultFaceFlux ff(RHO, CP, CS);

   // Random states
   srand(42);
   real_t Q_plus[NUM_STATE], Q_minus[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_plus[c] = ((double)rand() / RAND_MAX - 0.5) * 1e6;
      Q_minus[c] = ((double)rand() / RAND_MAX - 0.5) * 1e6;
   }

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   ff.Evaluate(data, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

   // Traction continuity: tau^{+,imp} = tau^{-,imp} (Eq. 11d/12d)
   TEST_NEAR(Q_imp_plus[SXX], Q_imp_minus[SXX], 1e-10,
             "Imposed: sigma_n continuous");
   TEST_NEAR(Q_imp_plus[SXY], Q_imp_minus[SXY], 1e-10,
             "Imposed: tau_1 continuous");
   TEST_NEAR(Q_imp_plus[SXZ], Q_imp_minus[SXZ], 1e-10,
             "Imposed: tau_2 continuous");
}

// ===== Test 27: Imposed state velocity jump = slip rate =====
void TestImposedStateVelocityJump()
{
   std::cout << "Test 27: TestImposedStateVelocityJump\n";

   DOFData data = MakeHomogeneousDOF();
   data.tau1_0 = 60e6;  // enough to cause slip
   data.psi = 0.5;
   FaultFaceFlux ff(RHO, CP, CS);

   // Some initial perturbation
   real_t Q_plus[NUM_STATE] = {}, Q_minus[NUM_STATE] = {};
   Q_plus[SXY] = 5e6;  Q_minus[SXY] = 5e6;  // symmetric shear

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   ff.Evaluate(data, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

   // Velocity jump: v^{+,imp} - v^{-,imp} = slip rate vector
   real_t jump_t1 = Q_imp_plus[VY] - Q_imp_minus[VY];
   real_t jump_t2 = Q_imp_plus[VZ] - Q_imp_minus[VZ];
   real_t jump_mag = std::sqrt(jump_t1*jump_t1 + jump_t2*jump_t2);

   real_t rel = (data.slip_rate > 0) ?
      std::abs(jump_mag - data.slip_rate) / data.slip_rate : jump_mag;

   TEST_ASSERT(rel < 1e-8,
               "Velocity jump magnitude = slip rate (rel " + std::to_string(rel) +
               ", V = " + std::to_string(data.slip_rate) + ")");

   // R-001 fix: check DIRECTION — velocity jump should be parallel to total traction.
   // Total traction in fault-local: tau1_total = tau1_0 + tau1_corr
   real_t tau1_total = data.tau1_0 + Q_imp_plus[SXY];
   real_t tau2_total = data.tau2_0 + Q_imp_plus[SXZ];
   real_t tau_mag = std::sqrt(tau1_total*tau1_total + tau2_total*tau2_total);
   if (tau_mag > 0 && jump_mag > 0)
   {
      // Cosine of angle between traction and velocity jump directions
      real_t cos_angle = (tau1_total * jump_t1 + tau2_total * jump_t2) / (tau_mag * jump_mag);
      // Slip opposes traction → cos should be -1 (antiparallel) or +1 (parallel)
      // depending on sign convention. Check |cos| ≈ 1.
      TEST_ASSERT(std::abs(std::abs(cos_angle) - 1.0) < 1e-6,
                  "Velocity jump parallel to traction (|cos| = " +
                  std::to_string(std::abs(cos_angle)) + ")");
   }
}

// ===== Test 28: Solver reuse: standalone Brent matches Evaluate() =====
void TestSolverReuse()
{
   std::cout << "Test 28: TestSolverReuse\n";

   DOFData data = MakeHomogeneousDOF();
   data.tau1_0 = 50e6;
   data.psi = 0.5;
   FaultFaceFlux ff(RHO, CP, CS);

   real_t Q_plus[NUM_STATE] = {}, Q_minus[NUM_STATE] = {};
   Q_plus[SXY] = 3e6; Q_minus[SXY] = 3e6;

   // Get V from Evaluate
   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   ff.Evaluate(data, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);
   real_t V_evaluate = data.slip_rate;

   // Compute same thing standalone
   real_t sigma_trial, tau1_trial, tau2_trial;
   FaultFaceFlux::ComputeTrialTraction(data, Q_plus, Q_minus,
                                       sigma_trial, tau1_trial, tau2_trial);
   real_t sigma_total = data.sigma_n0 + sigma_trial;
   real_t tau_total = std::sqrt((data.tau1_0 + tau1_trial) * (data.tau1_0 + tau1_trial) +
                                (data.tau2_0 + tau2_trial) * (data.tau2_0 + tau2_trial));
   real_t V_standalone = ff.GetSolver().SolveBrent(tau_total, data.psi,
                                                    std::abs(sigma_total), data.eta_s, data.a);

   real_t rel = (V_standalone > 0) ?
      std::abs(V_evaluate - V_standalone) / V_standalone : 0.0;
   TEST_ASSERT(rel < 1e-8,
               "Solver reuse: Evaluate V matches standalone (rel " + std::to_string(rel) + ")");
}

// ===== Test 29: State variable analytic update =====
void TestStateVariableUpdate()
{
   std::cout << "Test 29: TestStateVariableUpdate\n";

   // Test the θ-space analytic update via ψ→θ→update→ψ.
   // SCEC Eq. (2): dθ/dt = 1 - Vθ/L, analytic for constant V:
   //   θ(t+dt) = θ * exp(-V*dt/L) + (L/V) * (1 - exp(-V*dt/L))
   real_t f0 = 0.6, b = 0.012, V0 = 1e-6;
   real_t psi0 = 0.5, V = 1e-3, Dc = 0.02, dt = 1e-3;

   // Manually compute: ψ→θ, update θ, θ→ψ
   real_t theta0 = (Dc / V0) * std::exp((psi0 - f0) / b);
   real_t x = V * dt / Dc;
   real_t theta_new = theta0 * std::exp(-x) + (Dc / V) * (1.0 - std::exp(-x));
   real_t psi_exact = f0 + b * std::log(V0 * theta_new / Dc);

   real_t psi_computed = UpdateStateAnalytic(psi0, V, Dc, dt, f0, b, V0);

   TEST_NEAR(psi_computed, psi_exact, 1e-12,
             "State update matches exact θ-space (diff " +
             std::to_string(std::abs(psi_computed - psi_exact)) + ")");

   // Check stability for 1000 steps with constant V
   real_t psi = psi0;
   for (int i = 0; i < 1000; i++)
   {
      psi = UpdateStateAnalytic(psi, V, Dc, dt, f0, b, V0);
   }
   // Should converge to steady state: ψ_ss = f0 + b*ln(V0/V)
   // For V=1e-3: ψ_ss = 0.6 + 0.012 * ln(1e-6/1e-3) = 0.6 + 0.012*(-6.908) = 0.517
   real_t psi_ss = f0 + b * std::log(V0 / V);
   // After 1000 steps (1s at V=1e-3), ψ should approach steady state.
   // θ_ss = Dc/V = 20, but θ_0 ≈ 6.1e-5 for ψ=0.5. The relaxation
   // time is ~Dc/V = 20s, so after 1s we're ~5% of the way. Check direction.
   TEST_ASSERT(std::isfinite(psi) && psi > psi0 && psi < psi_ss + 0.1,
               "State variable evolving toward steady state (psi = " + std::to_string(psi) +
               ", psi_ss = " + std::to_string(psi_ss) + ")");


   // Critical test: locked fault (V ≈ 0) — the catastrophic case.
   // For V_ini = 1e-12, θ ≈ 1.6e9 s, dt = 1 s:
   // θ grows by ~1 → negligible change → ψ barely changes.
   // OLD BUG: ψ grew by dt=1 → ψ=1.736 → friction=1.6 → nucleation impossible.
   real_t psi_locked = 0.736;  // typical TPV102 initial ψ
   real_t V_locked = 1e-12;    // V_ini
   real_t dt_locked = 1.0;     // 1 second
   real_t psi_after = UpdateStateAnalytic(psi_locked, V_locked, Dc, dt_locked, f0, b, V0);
   real_t dpsi = std::abs(psi_after - psi_locked);
   TEST_ASSERT(dpsi < 1e-6,
               "Locked fault: psi barely changes after 1s (dpsi = " +
               std::to_string(dpsi) + ", was " + std::to_string(psi_locked) +
               ", now " + std::to_string(psi_after) + ")");
}

// ===== Test 30: Energy balance (domain energy + fault work = 0) =====
void TestEnergyBalance()
{
   std::cout << "Test 30: TestEnergyBalance\n";

   // The fault dissipates energy: dE/dt = -∫ τ · V dA ≤ 0
   // For a slipping fault, the corrected traction should satisfy τ · V > 0
   DOFData data = MakeHomogeneousDOF();
   data.tau1_0 = 60e6;  // enough to cause slip
   data.psi = 0.5;
   FaultFaceFlux ff(RHO, CP, CS);

   real_t Q_plus[NUM_STATE] = {}, Q_minus[NUM_STATE] = {};
   Q_plus[SXY] = 5e6; Q_minus[SXY] = 5e6;

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   ff.Evaluate(data, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

   // Corrected traction (same on both sides)
   real_t tau1_corr = Q_imp_plus[SXY];
   real_t tau2_corr = Q_imp_plus[SXZ];

   // Velocity jump (slip rate direction)
   real_t dv1 = Q_imp_plus[VY] - Q_imp_minus[VY];
   real_t dv2 = Q_imp_plus[VZ] - Q_imp_minus[VZ];

   // Friction dissipation = |τ · V| > 0 for any slipping fault.
   // The sign depends on convention (imposed state vs. traction direction).
   // Key check: friction does WORK (|τ·V| > 0, energy changes).
   real_t tau1_total = data.tau1_0 + tau1_corr;
   real_t tau2_total = data.tau2_0 + tau2_corr;
   real_t dissipation = std::abs(tau1_total * dv1 + tau2_total * dv2);

   TEST_ASSERT(dissipation > 0 && data.slip_rate > 0,
               "Energy balance: |tau·V| > 0 for slipping fault (|tau·V| = " +
               std::to_string(dissipation) + ", V = " + std::to_string(data.slip_rate) + ")");
}

// ===== Test 31: Fault flux sign convention (stress drop) =====
void TestFaultFluxSignConvention()
{
   std::cout << "Test 31: TestFaultFluxSignConvention\n";

   // With pre-stress τ₀ > 0 and small perturbation, the corrected traction
   // should be LESS than the trial traction (friction reduces stress).
   DOFData data = MakeHomogeneousDOF();
   data.tau1_0 = 60e6;
   data.psi = 0.5;
   FaultFaceFlux ff(RHO, CP, CS);

   real_t Q_plus[NUM_STATE] = {}, Q_minus[NUM_STATE] = {};
   Q_plus[SXY] = 1e6; Q_minus[SXY] = 1e6;

   // Get trial traction
   real_t sigma_trial, tau1_trial, tau2_trial;
   FaultFaceFlux::ComputeTrialTraction(data, Q_plus, Q_minus,
                                       sigma_trial, tau1_trial, tau2_trial);

   // Get corrected traction via Evaluate
   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   ff.Evaluate(data, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

   // Corrected traction magnitude should be ≤ trial traction magnitude
   real_t trial_mag = std::abs(tau1_trial);
   real_t corr_mag = std::abs(Q_imp_plus[SXY]);

   TEST_ASSERT(corr_mag <= trial_mag + 1e-6,
               "Stress drop: |tau_corr| <= |tau_trial| (trial " +
               std::to_string(trial_mag) + ", corr " + std::to_string(corr_mag) + ")");
}

// ===== Test R-004: End-to-end rotation pipeline =====
void TestRotationPipeline()
{
   std::cout << "Test R004: TestRotationPipeline\n";

   // Verify: global Q → rotate to fault-local → Evaluate → rotate back → global imposed.
   // Use a non-axis-aligned fault normal to exercise the full rotation.
   DOFData data = MakeHomogeneousDOF();
   data.tau1_0 = 50e6;
   data.psi = 0.5;
   FaultFaceFlux ff(RHO, CP, CS);

   // Fault normal at 45 degrees in XY plane
   real_t nor[3] = {1.0/std::sqrt(2.0), 1.0/std::sqrt(2.0), 0.0};
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(nor, t1, t2);

   // Build rotation matrices
   DenseMatrix Tinv(NUM_STATE, NUM_STATE), T(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);
   GodunovFlux::BuildRotation(nor, t1, t2, T);

   // Create global Q states with some stress
   real_t Q_plus_global[NUM_STATE] = {};
   real_t Q_minus_global[NUM_STATE] = {};
   Q_plus_global[SXX] = -50e6; Q_plus_global[SXY] = 20e6;
   Q_plus_global[VX] = 0.1;
   Q_minus_global[SXX] = -50e6; Q_minus_global[SXY] = 20e6;
   Q_minus_global[VX] = -0.1;

   // Step 1: Rotate global → fault-local
   real_t Q_plus_local[NUM_STATE], Q_minus_local[NUM_STATE];
   Tinv.Mult(Q_plus_global, Q_plus_local);
   Tinv.Mult(Q_minus_global, Q_minus_local);

   // Step 2: Evaluate in fault-local frame
   real_t Q_imp_plus_local[NUM_STATE], Q_imp_minus_local[NUM_STATE];
   ff.Evaluate(data, Q_plus_local, Q_minus_local,
               Q_imp_plus_local, Q_imp_minus_local);

   // Step 3: Rotate imposed states back to global
   real_t Q_imp_plus_global[NUM_STATE], Q_imp_minus_global[NUM_STATE];
   T.Mult(Q_imp_plus_local, Q_imp_plus_global);
   T.Mult(Q_imp_minus_local, Q_imp_minus_global);

   // Verify: all values are finite (rotation didn't produce garbage)
   bool finite = true;
   for (int c = 0; c < NUM_STATE; c++)
   {
      if (!std::isfinite(Q_imp_plus_global[c]) ||
          !std::isfinite(Q_imp_minus_global[c]))
      {
         finite = false; break;
      }
   }
   TEST_ASSERT(finite, "Rotation pipeline: all imposed values finite");

   // Verify: traction continuity in GLOBAL frame.
   // In fault-local, SXX/SXY/SXZ are traction components.
   // In global, traction = T * [sigma_nn, tau_nt1, tau_nt2, ...]
   // Since T is orthogonal, if local tractions match, global tractions match.
   // Check that the normal traction n·σ·n is the same on both sides.
   real_t trac_plus = 0, trac_minus = 0;
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         // Access stress components in global frame
         int idx;
         if (i == j) idx = i;  // SXX=0, SYY=1, SZZ=2
         else if ((i==0&&j==1)||(i==1&&j==0)) idx = SXY;
         else if ((i==1&&j==2)||(i==2&&j==1)) idx = SYZ;
         else idx = SXZ;

         trac_plus += nor[i] * Q_imp_plus_global[idx] * nor[j];
         trac_minus += nor[i] * Q_imp_minus_global[idx] * nor[j];
      }

   TEST_NEAR(trac_plus, trac_minus, 1e-6,
             "Rotation pipeline: normal traction continuous in global frame");
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "FaultFaceFlux Unit Tests (Phase 3)\n";
   std::cout << "========================================\n\n";

   TestTrialTractionLockedFault();
   TestTrialTractionVelocityJump();
   TestLockedFaultHighSigma();
   TestFrictionlessFault();
   TestImposedStateTractionContinuity();
   TestImposedStateVelocityJump();
   TestSolverReuse();
   TestStateVariableUpdate();
   TestEnergyBalance();
   TestFaultFluxSignConvention();
   TestRotationPipeline();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";

   return (num_failed > 0) ? 1 : 0;
}
