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

#ifdef SEAS_DIAG_FAULT_FLUX
   // v9.0.0 §0.5 DIAG gate for the C-1 / C-2 / C-3 bisection checkpoints.
   // Set true on exactly one hypocenter DOF (and optionally one off-hypo
   // witness) by the driver at init; all other DOFs keep diag_print=false
   // and the fprintf blocks are skipped.  Entire field is compiled out in
   // production builds so struct layout matches pre-change byte-for-byte.
   bool diag_print = false;
#endif
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

   /// Total-stress variant of Evaluate (I-06 Phase 3).  Expects Q_plus,
   /// Q_minus to carry TOTAL stresses (pre-stress + fluctuation); outputs
   /// Q_imp_± also carry TOTAL stresses.  On a symmetric-pre-stress fault
   /// with homogeneous material, mathematically equivalent to
   /// Evaluate(data, Q_fluc_±, ...) + pre-stress shift (plan §3).
   ///
   /// DOFData.sigma_n_corr, tau1_corr, tau2_corr are stored as TOTAL.
   /// Under the v9.3.0 migration, DOFData.sigma_n0, tau1_0, tau2_0 are
   /// ZEROED by the driver (Phase 4) so that EvaluateTotal does not
   /// double-count pre-stress through ComputeTrialTraction on TOTAL Q.
   ///
   /// R-003 guards preserved: psi-invariance is asserted at entry/exit
   /// (the driver's coupled RK4-on-psi integrator relies on psi-purity),
   /// and the SEAS_DIAG_FAULT_FLUX diagnostic print block mirrors the
   /// fluctuation path so C-1 / C-2 / C-3 bisection still functions.
   ///
   /// @param[in,out] data  Per-DOF state.  sigma_n_corr, tau1_corr,
   ///                      tau2_corr stored as TOTAL values.
   /// @param[in]  Q_plus  Total-stress state on + side (fault-local).
   /// @param[in]  Q_minus Total-stress state on - side (fault-local).
   /// @param[out] Q_imp_plus  Total-stress imposed state on + side.
   /// @param[out] Q_imp_minus Total-stress imposed state on - side.
   /// @param[in]  method Friction solver choice.
   void EvaluateTotal(DOFData &data,
                      const real_t *Q_plus, const real_t *Q_minus,
                      real_t *Q_imp_plus, real_t *Q_imp_minus,
                      FrictionSolver::Method method = FrictionSolver::Method::Brent) const;

   /// ADER Phase 5: time-integrated friction solve (fluctuation-Q variant).
   ///
   /// Inputs I± = ∫_0^{dt} Q±(τ) dτ in fault-local coordinates.  Converts
   /// to time-averaged Q̄± = I±/dt, calls `Evaluate` once on the averaged
   /// state (the ADER one-shot replacement for the 4 per-stage RK4 calls),
   /// and rescales the imposed outputs back to time-integrated form:
   ///   I_imp± = dt · Q_imp±.
   ///
   /// After the call, `data.{slip_rate,V1,V2,tau1_corr,tau2_corr,sigma_n_corr}`
   /// carry the TIME-AVERAGED values over [t_n, t_n+dt], as specified by
   /// plan §Phase 5 §4.  `data.psi` is NOT updated — the caller must call
   /// `UpdateStateAnalytic` with the returned time-averaged slip rate.
   ///
   /// In the dt→0 limit, `EvaluateADER(data, I±, dt)` ≈
   /// `Evaluate(data, Q±(t_n + dt/2), dt) · dt` to O(dt²)
   /// (averaging-then-solving vs solving-at-midpoint for the nonlinear
   /// friction law).
   ///
   /// @param[in,out] data  Per-DOF state.  Slip-rate & traction fields
   ///                      updated to time-averaged values.
   /// @param[in]  I_plus   Time-integrated + side state (9 components,
   ///                      fault-local).
   /// @param[in]  I_minus  Time-integrated − side state.
   /// @param[in]  dt       Time step.  Must be > 0.
   /// @param[out] I_imp_plus   Time-integrated imposed + state.
   /// @param[out] I_imp_minus  Time-integrated imposed − state.
   /// @param[in]  method   Friction solver choice.
   void EvaluateADER(DOFData &data,
                     const real_t *I_plus, const real_t *I_minus,
                     real_t dt,
                     real_t *I_imp_plus, real_t *I_imp_minus,
                     FrictionSolver::Method method = FrictionSolver::Method::Brent) const;

   /// ADER Phase 5 + v9.3.0 §Phase 7 follow-up: time-integrated friction
   /// solve with TOTAL-stress inputs.  Mirrors `EvaluateADER` but wraps
   /// `EvaluateTotal` internally so the ADER corrector can drive the
   /// post-I-06 TPV102 dispatch path without the explicit 1/dt-then-dt
   /// workaround mentioned in the v9.3.0 plan.
   ///
   /// @param[in,out] data  Per-DOF state; slip/traction fields updated
   ///                      to time-averaged TOTAL values.
   /// @param[in]  I_plus_tot  Total-stress time-integrated + state.
   /// @param[in]  I_minus_tot Total-stress time-integrated − state.
   /// @param[in]  dt          Time step (> 0).
   /// @param[out] I_imp_plus_tot   Total-stress imposed + state.
   /// @param[out] I_imp_minus_tot  Total-stress imposed − state.
   /// @param[in]  method      Friction solver choice.
   void EvaluateADERTotal(DOFData &data,
                          const real_t *I_plus_tot, const real_t *I_minus_tot,
                          real_t dt,
                          real_t *I_imp_plus_tot, real_t *I_imp_minus_tot,
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
