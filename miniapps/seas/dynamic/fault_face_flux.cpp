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
#include <cmath>
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
   // Step 1: Trial traction (Eq. 7)
   real_t sigma_n_trial, tau1_trial, tau2_trial;
   ComputeTrialTraction(data, Q_plus, Q_minus,
                        sigma_n_trial, tau1_trial, tau2_trial);

   // Total traction = pre-stress + trial (Eq. 8 setup)
   real_t sigma_n_total = data.sigma_n0 + sigma_n_trial;
   real_t tau1_total = data.tau1_0 + tau1_trial;
   real_t tau2_total = data.tau2_0 + tau2_trial;

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
}

} // namespace seas
} // namespace mfem
