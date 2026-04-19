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

#ifndef MFEM_SEAS_FAULT_FACE_FLUX_HPP
#define MFEM_SEAS_FAULT_FACE_FLUX_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "friction_solver.hpp"

namespace mfem
{
namespace seas
{

/// Per-DOF fault data: impedances, initial stress, state variable, slip rate.
struct DOFData
{
   real_t Zp_plus = 0, Zp_minus = 0;    ///< P-impedance on ± sides
   real_t Zs_plus = 0, Zs_minus = 0;    ///< S-impedance on ± sides
   real_t eta_p = 0, eta_s = 0;          ///< Harmonic mean impedances
   real_t sigma_n0 = 0;                  ///< Background normal stress (>0 compression)
   real_t tau1_0 = 0, tau2_0 = 0;        ///< Background shear pre-stress
   real_t a = 0.004;                      ///< Direct effect parameter
   real_t Dc = 0.14;                      ///< Critical slip distance [m]
   real_t psi = 0;                        ///< State variable (logarithmic)
   real_t slip_rate = 0;                  ///< Current slip rate |V| [m/s]
   real_t V1 = 0, V2 = 0;               ///< Slip rate components [m/s] (from Eq. 9)
   real_t slip1 = 0, slip2 = 0;          ///< Accumulated slip components

   // Corrected traction from Riemann solver (populated by FaultFaceFlux::Evaluate)
   real_t tau1_corr = 0, tau2_corr = 0;  ///< Corrected tangential traction [Pa]
   real_t sigma_n_corr = 0;              ///< Corrected normal traction [Pa]
};

/// @brief Fault-face Riemann solver for dynamic rupture.
///
/// Implements the trial-and-correction approach (plan Section 2.4):
/// 1. Compute trial traction from Godunov state (Eq. 7)
/// 2. Solve friction equation for slip rate (Eq. 8-9)
/// 3. Compute corrected traction (Eq. 10)
/// 4. Construct imposed states (Eq. 11-12)
///
/// Reuses FaultBasis for coordinate transforms and FrictionSolver for Eq. 8.
/// Reference: SeisSol FrictionSolverCommon.h, de la Puente et al. (2009).
class FaultFaceFlux
{
public:
   /// @brief Construct with material parameters.
   ///
   /// @param[in] rho  Density [kg/m³].
   /// @param[in] cp  P-wave speed [m/s].
   /// @param[in] cs  S-wave speed [m/s].
   FaultFaceFlux(real_t rho, real_t cp, real_t cs);

   /// @brief Compute trial traction from Q± states in fault-local coordinates.
   ///
   /// Implements Eq. (7a-c) from the plan:
   ///   σ_n^trial = η_p * (v_n⁻ - v_n⁺ + σ_n⁺/Z_p⁺ + σ_n⁻/Z_p⁻)
   ///   τ_1^trial = η_s * (v_t1⁻ - v_t1⁺ + τ_1⁺/Z_s⁺ + τ_1⁻/Z_s⁻)
   ///   τ_2^trial = η_s * (v_t2⁻ - v_t2⁺ + τ_2⁺/Z_s⁺ + τ_2⁻/Z_s⁻)
   ///
   /// @param[in] data  Per-DOF impedance data.
   /// @param[in] Q_plus  Rotated state on + side (9 components, fault-local).
   /// @param[in] Q_minus  Rotated state on − side (9 components, fault-local).
   /// @param[out] sigma_n_trial  Trial normal stress.
   /// @param[out] tau1_trial  Trial tangential traction (component 1).
   /// @param[out] tau2_trial  Trial tangential traction (component 2).
   static void ComputeTrialTraction(const DOFData &data,
                                    const real_t *Q_plus,
                                    const real_t *Q_minus,
                                    real_t &sigma_n_trial,
                                    real_t &tau1_trial,
                                    real_t &tau2_trial);

   /// @brief Full fault-face Riemann solver pipeline.
   ///
   /// Given Q± in fault-local coordinates, computes the imposed states
   /// Q^{±,imp} that satisfy friction and the characteristic compatibility
   /// relations. Returns the Godunov flux for Elem1 and Elem2.
   ///
   /// Pipeline: Eq. (7) → (8) → (9) → (10) → (11)-(12).
   ///
   /// @param[in,out] data  Per-DOF data (psi, slip_rate updated).
   /// @param[in] Q_plus  State on + side (fault-local, 9 components).
   /// @param[in] Q_minus  State on − side (fault-local, 9 components).
   /// @param[out] Q_imp_plus  Imposed state on + side (9 components).
   /// @param[out] Q_imp_minus  Imposed state on − side (9 components).
   /// @param[in] method  Friction solver method.
   void Evaluate(DOFData &data,
                 const real_t *Q_plus, const real_t *Q_minus,
                 real_t *Q_imp_plus, real_t *Q_imp_minus,
                 FrictionSolver::Method method = FrictionSolver::Method::Brent) const;

   /// Access the friction solver.
   const FrictionSolver &GetSolver() const { return solver_; }

private:
   real_t rho_, cp_, cs_;
   real_t Zp_, Zs_;  ///< Impedances (homogeneous)
   FrictionSolver solver_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_FACE_FLUX_HPP
