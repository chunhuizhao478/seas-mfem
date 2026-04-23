// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "fault_face_flux.hpp"
#include "seas_diag_rank.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>

namespace mfem
{
namespace seas
{

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
FaultFaceFlux::FaultFaceFlux(real_t rho, real_t cp, real_t cs)
   : rho_(rho), cp_(cp), cs_(cs)
{
   Zp_ = rho * cp;
   Zs_ = rho * cs;
}

// ---------------------------------------------------------------------------
// Trial traction: Eq. (7a-c)
// In fault-local coordinates: x = normal, y = tangent1, z = tangent2.
// So SXX = sigma_nn, SXY = tau_1 (nt1), SXZ = tau_2 (nt2).
// VX = v_n, VY = v_t1, VZ = v_t2.
// ---------------------------------------------------------------------------
void FaultFaceFlux::ComputeTrialTraction(const DOFData &data,
                                         const real_t *Q_plus,
                                         const real_t *Q_minus,
                                         real_t &sigma_n_trial,
                                         real_t &tau1_trial,
                                         real_t &tau2_trial)
{
   real_t invZp_plus = 1.0 / data.Zp_plus;
   real_t invZp_minus = 1.0 / data.Zp_minus;
   real_t invZs_plus = 1.0 / data.Zs_plus;
   real_t invZs_minus = 1.0 / data.Zs_minus;

   // Eq. (7a): sigma_n^trial = eta_p * (v_n^- - v_n^+ + sigma_n^+/Zp^+ + sigma_n^-/Zp^-)
   sigma_n_trial = data.eta_p * (Q_minus[VX] - Q_plus[VX]
                                 + Q_plus[SXX] * invZp_plus
                                 + Q_minus[SXX] * invZp_minus);

   // Eq. (7b): tau_1^trial = eta_s * (v_t1^- - v_t1^+ + tau_1^+/Zs^+ + tau_1^-/Zs^-)
   tau1_trial = data.eta_s * (Q_minus[VY] - Q_plus[VY]
                              + Q_plus[SXY] * invZs_plus
                              + Q_minus[SXY] * invZs_minus);

   // Eq. (7c): tau_2^trial = eta_s * (v_t2^- - v_t2^+ + tau_2^+/Zs^+ + tau_2^-/Zs^-)
   tau2_trial = data.eta_s * (Q_minus[VZ] - Q_plus[VZ]
                              + Q_plus[SXZ] * invZs_plus
                              + Q_minus[SXZ] * invZs_minus);
}

// ---------------------------------------------------------------------------
// Full Evaluate pipeline: Eq. (7) → (8) → (9) → (10) → (11)-(12)
// ---------------------------------------------------------------------------
void FaultFaceFlux::Evaluate(DOFData &data,
                             const real_t *Q_plus, const real_t *Q_minus,
                             real_t *Q_imp_plus, real_t *Q_imp_minus,
                             FrictionSolver::Method method) const
{
   // v9.2.0 F01+F02 invariant (REVIEW R-V92-H07): the driver's coupled
   // RK4 on (Q, psi) in `drivers/tpv102_driver.cpp` reads `data.psi` at
   // entry to this call as the stage-local psi and relies on it being
   // UNCHANGED on return.  A future refactor that adds a psi write here
   // would silently invalidate the RK4 stage arithmetic (the driver
   // would double-integrate psi).  Keep this guard whenever psi is not
   // a formal input to the Riemann solver.
#ifndef NDEBUG
   const real_t psi_at_entry = data.psi;
#endif

   // v9.0.0 Pelties-9 per-side flux (see
   // debug_document/tpv102_debug_document/tpv102_debug_v9.0.0_seissol_flux_comparison.md
   // §10.1 and §18 R-F08).  The call site in wave_operator.inl applies a
   // SINGLE GodunovFlux instance to both sides of the fault, which is
   // correct only when A_plus == A_minus (homogeneous material).  Guard
   // against silent use on a bimaterial face.
   // R-008: allow ~1-ULP drift — reject only when the two sides are
   // observably different materials, not when they're the same material
   // reached via two different expression chains (e.g. a future driver
   // computing Zp_plus, Zp_minus from separate rho*cp expressions).
   auto homog_ok = [](real_t a, real_t b)
   {
      return std::abs(a - b) <=
             1e-12 * std::max(std::abs(a), std::abs(b));
   };
   MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) &&
               homog_ok(data.Zs_plus, data.Zs_minus),
               "Bimaterial fault face detected (Zp_plus=" << data.Zp_plus
               << " Zp_minus=" << data.Zp_minus
               << " Zs_plus=" << data.Zs_plus
               << " Zs_minus=" << data.Zs_minus
               << ").  v9.0.0 Pelties-9 per-side flux assumes "
               "homogeneous material.  Extend GodunovFlux to per-side A "
               "before running this configuration.");

   // Step 1: Trial traction (Eq. 7)
   real_t sigma_n_trial, tau1_trial, tau2_trial;
   ComputeTrialTraction(data, Q_plus, Q_minus,
                        sigma_n_trial, tau1_trial, tau2_trial);

   // Total traction = pre-stress + trial (Eq. 8 setup)
   real_t sigma_n_total = data.sigma_n0 + sigma_n_trial;
   real_t tau1_total = data.tau1_0 + tau1_trial;
   real_t tau2_total = data.tau2_0 + tau2_trial;

#ifdef SEAS_DIAG_FAULT_FLUX
   // C-1 EVAL: v9.0.0 §0.5 checkpoint — printf only, no MPI.  Prints only
   // on DOFs the driver flagged as diagnostic (typically the hypocenter QP).
   if (data.diag_print)
   {
      std::fprintf(stderr,
         "[C-1 EVAL] rank=%d  tau1_trial=%+.3e Pa  tau2_trial=%+.3e Pa  "
         "|Q_plus[VY]|=%.3e  |Q_plus[SXY]|=%.3e  |Q_plus[VZ]|=%.3e  "
         "|Q_plus[SXZ]|=%.3e  psi=%.3e\n",
         g_seas_my_rank,
         tau1_trial, tau2_trial,
         std::abs(Q_plus[VY]),  std::abs(Q_plus[SXY]),
         std::abs(Q_plus[VZ]),  std::abs(Q_plus[SXZ]),
         data.psi);
   }
#endif


   // Traction magnitude Θ = sqrt(tau1_total² + tau2_total²)
   real_t Theta = std::sqrt(tau1_total * tau1_total + tau2_total * tau2_total);

   // Step 2: Solve friction equation for |V̂| (Eq. 8)
   real_t V_abs = 0.0;
   if (Theta > 0.0)
   {
      V_abs = solver_.Solve(Theta, data.psi, std::abs(sigma_n_total),
                            data.eta_s, data.a, method);
   }

   // Step 3: Decompose slip rate into components (Eq. 9)
   real_t V1 = 0.0, V2 = 0.0;
   real_t tau1_corr = tau1_trial, tau2_corr = tau2_trial;

   if (Theta > 0.0 && V_abs > 0.0)
   {
      // Friction strength
      real_t C = std::exp(data.psi / data.a) / (2.0 * FrictionSolver::V0);
      real_t f_V = data.a * std::asinh(V_abs * C);
      real_t strength = std::abs(sigma_n_total) * f_V;

      // Slip rate decomposition (Eq. 9)
      V1 = V_abs * (tau1_total) / (strength + data.eta_s * V_abs);
      V2 = V_abs * (tau2_total) / (strength + data.eta_s * V_abs);

      // Step 4: Corrected traction (Eq. 10)
      tau1_corr = tau1_trial - data.eta_s * V1;
      tau2_corr = tau2_trial - data.eta_s * V2;
   }

   real_t sigma_n_corr = sigma_n_trial;  // normal traction unchanged by slip

   // Step 5: Construct imposed states (Eq. 11-12)
   // Initialize with current state
   std::memcpy(Q_imp_minus, Q_minus, NUM_STATE * sizeof(real_t));
   std::memcpy(Q_imp_plus, Q_plus, NUM_STATE * sizeof(real_t));

   // Minus side (Eq. 11a-d): v^{-,imp} = v^- - (1/Z)(sigma_corr - sigma^-)
   real_t invZp_m = 1.0 / data.Zp_minus;
   real_t invZs_m = 1.0 / data.Zs_minus;

   Q_imp_minus[VX]  = Q_minus[VX]  - invZp_m * (sigma_n_corr - Q_minus[SXX]);
   Q_imp_minus[VY]  = Q_minus[VY]  - invZs_m * (tau1_corr - Q_minus[SXY]);
   Q_imp_minus[VZ]  = Q_minus[VZ]  - invZs_m * (tau2_corr - Q_minus[SXZ]);

   // Plus side (Eq. 12a-d): v^{+,imp} = v^+ + (1/Z)(sigma_corr - sigma^+)
   real_t invZp_p = 1.0 / data.Zp_plus;
   real_t invZs_p = 1.0 / data.Zs_plus;

   Q_imp_plus[VX]  = Q_plus[VX]  + invZp_p * (sigma_n_corr - Q_plus[SXX]);
   Q_imp_plus[VY]  = Q_plus[VY]  + invZs_p * (tau1_corr - Q_plus[SXY]);
   Q_imp_plus[VZ]  = Q_plus[VZ]  + invZs_p * (tau2_corr - Q_plus[SXZ]);

   // Both sides: imposed stress = corrected traction (Eq. 11d/12d)
   Q_imp_minus[SXX] = sigma_n_corr;
   Q_imp_minus[SXY] = tau1_corr;
   Q_imp_minus[SXZ] = tau2_corr;

   Q_imp_plus[SXX]  = sigma_n_corr;
   Q_imp_plus[SXY]  = tau1_corr;
   Q_imp_plus[SXZ]  = tau2_corr;

   // Non-normal stresses are unchanged (Eq. 11d/12d: sigma_yy, sigma_zz, sigma_yz)
   // Already copied from Q_plus/Q_minus above.

   // Update fault state
   data.slip_rate = V_abs;
   data.V1 = V1;
   data.V2 = V2;
   data.tau1_corr = data.tau1_0 + tau1_corr;  // total corrected traction
   data.tau2_corr = data.tau2_0 + tau2_corr;
   data.sigma_n_corr = data.sigma_n0 + sigma_n_corr;

#ifndef NDEBUG
   // R-V92-H07: psi must be pristine for the driver's coupled RK4 to be
   // correct.  Bit-exact equality is the right check — no arithmetic on
   // psi has happened between entry and exit.
   MFEM_ASSERT(data.psi == psi_at_entry,
               "FaultFaceFlux::Evaluate mutated data.psi "
               "(before = " << psi_at_entry << ", after = " << data.psi
               << ").  The driver's coupled-RK4-on-psi integrator assumes "
               "this function is psi-pure.  See plan Step 5 + REVIEW "
               "R-V92-H07.");
#endif
}

// ---------------------------------------------------------------------------
// Total-stress pipeline: Evaluate operating on Q that already carries the
// pre-stress tensor in every bulk DOF (v9.3.0 Phase 3, I-06 part A).
//
// Differences from Evaluate:
//   - Inputs Q_plus / Q_minus are TOTAL (pre-stress + fluctuation).
//     ComputeTrialTraction on TOTAL Q therefore produces a TOTAL trial
//     traction directly — no additional +sigma_n0/+tau_i_0 shift needed.
//   - DOFData.sigma_n0 / tau1_0 / tau2_0 are assumed ZEROED by the driver
//     under the Phase 4 migration so pre-stress is not double-counted.
//   - DOFData.sigma_n_corr / tau1_corr / tau2_corr are stored as TOTAL
//     directly (no sigma_n0 add).  This matches SeisSol's
//     RateAndState.h:222-240 convention.
//
// The imposed-state construction (Step 4) is byte-identical to the
// fluctuation path: the differences (sigma_n_corr - Q_{plus,minus}[SXX])
// are invariant under a constant pre-stress shift of both sides.
//
// R-003 guards: #ifndef NDEBUG psi-invariance entry/exit MFEM_ASSERT +
// SEAS_DIAG_FAULT_FLUX diag print block after ComputeTrialTraction.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateTotal(DOFData &data,
                                  const real_t *Q_plus, const real_t *Q_minus,
                                  real_t *Q_imp_plus, real_t *Q_imp_minus,
                                  FrictionSolver::Method method) const
{
#ifndef NDEBUG
   // R-003 (review-incorporation plan): the driver's coupled-RK4-on-psi
   // integrator requires this function to be psi-pure.  Enforce the same
   // bit-exact guard as in Evaluate.
   const real_t psi_at_entry = data.psi;

   // R-I06-007 (review-round-2): EvaluateTotal's correctness contract is
   // that bulk Q carries the pre-stress, so DOFData's pre-stress fields
   // MUST be zero — the fluctuation-mode seeding done by
   // InitializeFaultDOFs must be followed by ZeroDOFDataPreStressTotal
   // before any EvaluateTotal call.  A future refactor that removes the
   // zeroing (or a "unification" that folds pre-stress back into
   // ComputeTrialTraction) would silently double-count the pre-stress
   // in the trial traction.  Assert the contract at runtime.
   MFEM_ASSERT(data.sigma_n0 == 0.0 &&
               data.tau1_0   == 0.0 &&
               data.tau2_0   == 0.0,
               "EvaluateTotal contract violated: DOFData pre-stress "
               "fields must be zero (sigma_n0=" << data.sigma_n0 <<
               ", tau1_0=" << data.tau1_0 <<
               ", tau2_0=" << data.tau2_0 <<
               ").  Total-stress drivers must call "
               "ZeroDOFDataPreStressTotal after InitializeFaultDOFs.");
#endif

   // Homogeneous check — identical to Evaluate.  v9.0.0 Pelties-9 per-side
   // flux assumes A_plus == A_minus on a fault face.
   auto homog_ok = [](real_t a, real_t b)
   {
      return std::abs(a - b) <=
             1e-12 * std::max(std::abs(a), std::abs(b));
   };
   MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) &&
               homog_ok(data.Zs_plus, data.Zs_minus),
               "Bimaterial fault face detected (Zp_plus=" << data.Zp_plus
               << " Zp_minus=" << data.Zp_minus
               << " Zs_plus=" << data.Zs_plus
               << " Zs_minus=" << data.Zs_minus
               << ").  EvaluateTotal assumes homogeneous material.  "
               "Extend GodunovFlux / FaultFaceFlux to per-side A before "
               "running this configuration.");

   // Step 1: Trial traction from bulk Q.  Under total-Q, bulk Q
   //         carries the static background prestress, so this trial is
   //         "background + wave" — but does NOT yet include the
   //         persistent nucleation channel.
   real_t sigma_n_trial, tau1_trial, tau2_trial;
   ComputeTrialTraction(data, Q_plus, Q_minus,
                        sigma_n_trial, tau1_trial, tau2_trial);

   // Step 1b: friction inputs = trial + persistent-nucleation channel.
   // Mirrors SeisSol's RateAndState.h:222-240 pattern, where
   // `totalTraction* = initialStressInFaultCS* + faultStresses.traction*`
   // is a temporary used ONLY for the friction solve and the
   // slip-rate divisor — never fed back into `tractionResults`
   // (which is the analog of our `tau*_corr` consumed by the Riemann
   // imposed state below).  Mutating `tau*_trial` here would leak
   // `tau*_nuc` into the Riemann velocity-jump
   // `(2/Zs)·(tau_corr - Q_bulk[SXY])`, producing a ~24× spurious
   // shear pulse into bulk Q every step at full nucleation
   // (tau2_nuc=25 MPa, V_abs≈0.22 m/s, Zs≈9.25 MPa·s/m).  Keep the
   // two scales separate.
   const real_t sigma_n_fric = sigma_n_trial + data.sigma_n_nuc;
   const real_t tau1_fric    = tau1_trial    + data.tau1_nuc;
   const real_t tau2_fric    = tau2_trial    + data.tau2_nuc;

#ifdef SEAS_DIAG_FAULT_FLUX
   // R-003: mirror Evaluate's C-1 checkpoint for total mode.  Same
   // diag_print gate, same fprintf semantics (no MPI).  Prints only
   // on DOFs the driver flagged as diagnostic.  Print the friction-
   // input total (what the solver actually sees); a separate trial
   // value can be reconstructed as tau*_fric - data.tau*_nuc if needed.
   if (data.diag_print)
   {
      std::fprintf(stderr,
         "[C-1 EVAL-TOTAL] rank=%d  tau1_fric=%+.3e Pa  "
         "tau2_fric=%+.3e Pa  sigma_n_fric=%+.3e Pa  psi=%.3e  "
         "a=%.3e  eta_s=%.3e  (tau2_nuc=%+.3e, sigma_n_nuc=%+.3e)\n",
         g_seas_my_rank, tau1_fric, tau2_fric, sigma_n_fric, data.psi,
         data.a, data.eta_s, data.tau2_nuc, data.sigma_n_nuc);
   }
#endif

   // Step 2: Solve friction equation for |V| on the FRICTION-INPUT
   //         total (trial + nuc).  Theta is the magnitude of the
   //         total tangential traction the friction law must balance.
   const real_t Theta = std::sqrt(tau1_fric * tau1_fric
                                + tau2_fric * tau2_fric);
   real_t V_abs = 0.0;
   if (Theta > 0.0)
   {
      V_abs = solver_.Solve(Theta, data.psi, std::abs(sigma_n_fric),
                            data.eta_s, data.a, method);
   }

   // Step 3: Slip-rate decomposition uses the TOTAL traction in the
   //         numerator (SeisSol numerator is `totalTraction*`).
   //         Riemann-side corrected traction tau*_corr stays on the
   //         TRIAL scale (no nuc).
   real_t V1 = 0.0, V2 = 0.0;
   real_t tau1_corr = tau1_trial, tau2_corr = tau2_trial;
   if (Theta > 0.0 && V_abs > 0.0)
   {
      real_t C = std::exp(data.psi / data.a) / (2.0 * FrictionSolver::V0);
      real_t f_V = data.a * std::asinh(V_abs * C);
      real_t strength = std::abs(sigma_n_fric) * f_V;
      V1 = V_abs * tau1_fric / (strength + data.eta_s * V_abs);
      V2 = V_abs * tau2_fric / (strength + data.eta_s * V_abs);
      tau1_corr = tau1_trial - data.eta_s * V1;
      tau2_corr = tau2_trial - data.eta_s * V2;
   }

   // Step 4: Imposed states.  `*_corr` here is TRIAL-scale (matches
   //         SeisSol `tractionResults.traction*`), so the Riemann
   //         velocity jump `(2/Zs)·(tau_corr - Q_bulk[SXY])` carries
   //         only the friction reaction `eta_s·V` — the persistent
   //         nucleation channel does NOT radiate through bulk Q.
   const real_t sigma_n_corr = sigma_n_trial;
   std::memcpy(Q_imp_minus, Q_minus, NUM_STATE * sizeof(real_t));
   std::memcpy(Q_imp_plus,  Q_plus,  NUM_STATE * sizeof(real_t));

   const real_t invZp_m = 1.0 / data.Zp_minus;
   const real_t invZs_m = 1.0 / data.Zs_minus;
   const real_t invZp_p = 1.0 / data.Zp_plus;
   const real_t invZs_p = 1.0 / data.Zs_plus;

   Q_imp_minus[VX] = Q_minus[VX] - invZp_m * (sigma_n_corr - Q_minus[SXX]);
   Q_imp_minus[VY] = Q_minus[VY] - invZs_m * (tau1_corr    - Q_minus[SXY]);
   Q_imp_minus[VZ] = Q_minus[VZ] - invZs_m * (tau2_corr    - Q_minus[SXZ]);
   Q_imp_plus[VX]  = Q_plus[VX]  + invZp_p * (sigma_n_corr - Q_plus[SXX]);
   Q_imp_plus[VY]  = Q_plus[VY]  + invZs_p * (tau1_corr    - Q_plus[SXY]);
   Q_imp_plus[VZ]  = Q_plus[VZ]  + invZs_p * (tau2_corr    - Q_plus[SXZ]);

   Q_imp_minus[SXX] = sigma_n_corr;  Q_imp_plus[SXX] = sigma_n_corr;
   Q_imp_minus[SXY] = tau1_corr;     Q_imp_plus[SXY] = tau1_corr;
   Q_imp_minus[SXZ] = tau2_corr;     Q_imp_plus[SXZ] = tau2_corr;

   // Step 5: Update DOFData — store TOTAL (trial + nuc) for station
   //         output and post-step diagnostics.  Mirrors the
   //         fluctuation-path `data.tau*_corr = data.tau*_0 + tau*_corr`
   //         pattern: post-call DOFData carries the full physical
   //         traction the fault is bearing, while the Riemann path
   //         above operated on the trial scale.
   data.slip_rate    = V_abs;
   data.V1           = V1;
   data.V2           = V2;
   data.tau1_corr    = tau1_corr   + data.tau1_nuc;
   data.tau2_corr    = tau2_corr   + data.tau2_nuc;
   data.sigma_n_corr = sigma_n_corr + data.sigma_n_nuc;

#ifndef NDEBUG
   MFEM_ASSERT(data.psi == psi_at_entry,
               "FaultFaceFlux::EvaluateTotal mutated data.psi "
               "(before = " << psi_at_entry << ", after = " << data.psi
               << ").  R-V92-H07: the driver's coupled-RK4-on-psi "
               "integrator assumes this function is psi-pure.");
#endif
}

// ---------------------------------------------------------------------------
// FACE-AVERAGED EvaluateTotal (option 1, 2026-04-22 pepper fix).
//
// Calls EvaluateTotal on the FACE-AVERAGED tau_trial / sigma_n_trial
// using a representative DOFData snapshot (averaged psi, averaged
// tau*_nuc) for the friction solve, then writes the friction outputs
// uniformly to every per-QP DOFData entry of the face.  The wave
// operator can then deposit the face-uniform Q_imp via shape1·F_h at
// every QP — all per-QP F_h values are identical.
//
// Effect: eliminates per-QP rhs deposition variation that compounds
// through the DG×nonlinear-friction amplification chain (per-DOF
// non-uniformity → per-QP friction outputs differ → per-QP rhs
// variation → next stage sees more variation, etc.).  Diagnosed in
// test_adjacent_triangle_fault_uniformity; see commits e686867 +
// 81ae35a for the full diagnostic chain.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateTotalFaceAveraged(
   DOFData *dof_data, int nqp_per_face,
   const real_t *Q_plus_avg, const real_t *Q_minus_avg,
   real_t *Q_imp_plus, real_t *Q_imp_minus,
   FrictionSolver::Method method) const
{
   MFEM_VERIFY(dof_data != nullptr,
               "EvaluateTotalFaceAveraged: dof_data must not be null");
   MFEM_VERIFY(nqp_per_face > 0,
               "EvaluateTotalFaceAveraged: nqp_per_face must be > 0, got "
               << nqp_per_face);

   // Build a face-representative DOFData by averaging the per-QP fields
   // that the friction solver reads (psi, tau*_nuc, sigma_n_nuc).
   // Material parameters (Zp, Zs, eta_p, eta_s, a, Dc) are uniform per
   // face by physics; assert the first QP matches the rest as a sanity
   // gate.  V1, V2, slip_rate, slip1, slip2, tau*_corr, sigma_n_corr
   // are OUTPUTS — overwritten by EvaluateTotal.
   DOFData face_data = dof_data[0];
   real_t psi_avg = 0.0, tau1_nuc_avg = 0.0, tau2_nuc_avg = 0.0;
   real_t sigma_n_nuc_avg = 0.0;
   for (int q = 0; q < nqp_per_face; q++)
   {
      psi_avg          += dof_data[q].psi;
      tau1_nuc_avg     += dof_data[q].tau1_nuc;
      tau2_nuc_avg     += dof_data[q].tau2_nuc;
      sigma_n_nuc_avg  += dof_data[q].sigma_n_nuc;
   }
   const real_t inv_nqp = 1.0 / static_cast<real_t>(nqp_per_face);
   face_data.psi         = psi_avg          * inv_nqp;
   face_data.tau1_nuc    = tau1_nuc_avg     * inv_nqp;
   face_data.tau2_nuc    = tau2_nuc_avg     * inv_nqp;
   face_data.sigma_n_nuc = sigma_n_nuc_avg  * inv_nqp;
   face_data.tau1_0      = 0.0;  // total-Q contract: must be zero
   face_data.tau2_0      = 0.0;
   face_data.sigma_n0    = 0.0;

   // Run the standard EvaluateTotal on face-averaged inputs.  The
   // outputs (V1, V2, slip_rate, tau*_corr, sigma_n_corr) are written
   // into face_data; Q_imp_plus / Q_imp_minus are the face-uniform
   // imposed states.
   EvaluateTotal(face_data, Q_plus_avg, Q_minus_avg,
                 Q_imp_plus, Q_imp_minus, method);

   // Distribute the friction outputs UNIFORMLY across every per-QP
   // DOFData entry of this face.  psi is intentionally NOT
   // overwritten: the caller's coupled-RK4-on-psi integrator owns
   // psi evolution per QP, and the EvaluateTotal contract guarantees
   // psi is unchanged inside.  We restore each QP's per-QP psi
   // (already preserved by EvaluateTotal contract on face_data; we
   // just leave per-QP dof_data[q].psi untouched).
   for (int q = 0; q < nqp_per_face; q++)
   {
      dof_data[q].slip_rate    = face_data.slip_rate;
      dof_data[q].V1           = face_data.V1;
      dof_data[q].V2           = face_data.V2;
      dof_data[q].tau1_corr    = face_data.tau1_corr;
      dof_data[q].tau2_corr    = face_data.tau2_corr;
      dof_data[q].sigma_n_corr = face_data.sigma_n_corr;
      // psi: untouched (per-QP integrated by driver's coupled RK4).
      // {tau*_0, tau*_nuc, sigma_n_*} are INPUTS to friction; not
      // overwritten here.
   }
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 5: time-integrated Riemann solve (fluctuation variant).
// ---------------------------------------------------------------------------
// I± = ∫_0^{dt} Q±(τ) dτ  ⇒  Q̄± = I±/dt.  Call the standard Evaluate on
// the time-averaged state (which now represents the mean Q over [t_n,
// t_n+dt]); the friction solver returns tractions / slip-rate that are
// themselves time-averaged in the same sense.  Finally rescale the
// imposed Q-states back to time-integrated form: I_imp± = dt · Q_imp±.
//
// Proof of O(dt²) consistency with a midpoint-RK4 call: let
// Q(τ) = Q_m + O(dt),  Q_m = Q(t_n + dt/2).  Then Q̄ = Q_m + O(dt²) by
// Simpson's rule on smooth Q, so feeding Q̄ into the nonlinear friction
// equation gives the same root as feeding Q_m up to O(dt²) (Brent's
// root of a Lipschitz-smooth residual commutes with O(dt²) input
// perturbations).  Imposed states scale linearly in the inputs modulo
// the friction-law nonlinearity, again O(dt²).  Plan §Phase 5 §4.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateADER(DOFData &data,
                                 const real_t *I_plus, const real_t *I_minus,
                                 real_t dt,
                                 real_t *I_imp_plus, real_t *I_imp_minus,
                                 FrictionSolver::Method method) const
{
   MFEM_VERIFY(dt > 0.0,
               "FaultFaceFlux::EvaluateADER: dt must be > 0, got " << dt);

   real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
   const real_t inv_dt = 1.0 / dt;
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_avg_plus[c]  = I_plus[c]  * inv_dt;
      Q_avg_minus[c] = I_minus[c] * inv_dt;
   }

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   Evaluate(data, Q_avg_plus, Q_avg_minus, Q_imp_plus, Q_imp_minus, method);

   for (int c = 0; c < NUM_STATE; c++)
   {
      I_imp_plus[c]  = Q_imp_plus[c]  * dt;
      I_imp_minus[c] = Q_imp_minus[c] * dt;
   }
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 5 + v9.3.0 §Phase 7: time-integrated Riemann solve
// (total-stress variant).  Same 1/dt ↔ dt scaling wrapping pattern as
// EvaluateADER, but dispatches to EvaluateTotal so TPV102's post-I-06
// total-Q path has a direct ADER entry point without the workaround the
// v9.3.0 plan describes.
// ---------------------------------------------------------------------------
void FaultFaceFlux::EvaluateADERTotal(DOFData &data,
                                      const real_t *I_plus_tot,
                                      const real_t *I_minus_tot,
                                      real_t dt,
                                      real_t *I_imp_plus_tot,
                                      real_t *I_imp_minus_tot,
                                      FrictionSolver::Method method) const
{
   MFEM_VERIFY(dt > 0.0,
               "FaultFaceFlux::EvaluateADERTotal: dt must be > 0, got " << dt);

   real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
   const real_t inv_dt = 1.0 / dt;
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_avg_plus[c]  = I_plus_tot[c]  * inv_dt;
      Q_avg_minus[c] = I_minus_tot[c] * inv_dt;
   }

   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
   EvaluateTotal(data, Q_avg_plus, Q_avg_minus, Q_imp_plus, Q_imp_minus,
                 method);

   for (int c = 0; c < NUM_STATE; c++)
   {
      I_imp_plus_tot[c]  = Q_imp_plus[c]  * dt;
      I_imp_minus_tot[c] = Q_imp_minus[c] * dt;
   }
}

} // namespace seas
} // namespace mfem
