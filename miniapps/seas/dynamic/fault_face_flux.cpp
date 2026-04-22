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

} // namespace seas
} // namespace mfem
