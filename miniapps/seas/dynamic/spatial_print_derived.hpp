// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/spatial_print_derived.hpp — Phase D of
// safs/project_7.0_alternative/document/05_18_2026/PLAN_first_safs_run.md.
//
// Pre-flight printer for `seas_spatial_dyn_driver`.  Replaces the 3-line
// stub in the driver with a hard gate that prints per-region
// histograms + L_nuc + max-ratio statistics on rank 0, and ABORTS on
// initial-equilibrium violation when `cfg.abort_on_failure == true`.
//
// A user that runs `--dry-run --print-derived` on Frontera login gets a
// hard fail before submitting if their stress / friction / nucleation
// setup would spontaneously rupture or fail to nucleate.

#ifndef MFEM_SEAS_SPATIAL_PRINT_DERIVED_HPP
#define MFEM_SEAS_SPATIAL_PRINT_DERIVED_HPP

#include "mfem.hpp"

#include "spatial_nucleation.hpp"
#include "../spatial/code/spatial_friction.hpp"

#include <iostream>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace mfem
{
namespace seas
{
namespace spatial
{

struct PrintDerivedConfig
{
   bool   enabled                = false;  ///< mirrors --print-derived
   bool   abort_on_failure       = true;   ///< false ⇒ warn-only
   real_t outside_safety_factor  = 3.0;    ///< r > F · max(radius_*) ⇒ outside
};

/// @brief Pre-flight pass — prints per-region histograms + L_nuc +
/// max-ratio statistics on rank 0; aborts (MFEM_ABORT) on initial-
/// equilibrium violation when `cfg.abort_on_failure == true`.
///
/// @returns the maximum `|τ_pre| / (μ_s · σ_n_eff)` ratio observed at
/// any DOF OUTSIDE the nucleation Gaussian's support
/// (`r > outside_safety_factor · max(radius_dip, radius_strike)`).
/// When `cfg.enabled == false`, no I/O; returns 0.0 and is a no-op.
///
/// Honors the env var `SEAS_SKIP_EQUILIBRIUM_GATE=1`: any failure that
/// would otherwise MFEM_ABORT is downgraded to a rank-0 WARNING and
/// the function returns the offending ratio.
real_t PrintDerivedAndCheck(
   const PrintDerivedConfig&                cfg,
   const SlipWeakeningPerDOFParams&         lsw,
   const Vector&                            tau_pre_per_dof,
   const Vector&                            sigma_n_eff_per_dof,
   const Vector&                            dof_coords_3d,
   const NucleationSpec&                    nuc,
   const GradualOverstressPerDOFParams&     nuc_params,
   const StressSpec&                        stress,
   real_t                                   mu_bulk,
   real_t                                   cp,
   real_t                                   cs,
   real_t                                   h_min_global,
   real_t                                   dt_cfl,
   real_t                                   tfinal,
   int                                      num_fault_global,
   int                                      num_zero_normal_fallbacks
#ifdef MFEM_USE_MPI
   , MPI_Comm comm
#endif
   , int                                    rank,
   std::ostream&                            out = std::cout);

/// @brief Rate-and-state (aging-law) pre-flight pass — the RS sibling of the
/// LSW @ref PrintDerivedAndCheck above.
///
/// **PLAN DEVIATION (documented):** the Phase-3 plan
/// (PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24) does NOT specify a
/// rate-and-state `--print-derived` path; the LSW overload aborts on an RS run
/// because it dereferences `lsw.mu_s/mu_d/d_c` (absent for RS).  This overload
/// is added so the Phase-3 SAFS-RS sbatch (which passes `--print-derived`) runs
/// to completion.  All RS formulas are grounded in this repo's own machinery
/// (no Tandem reference, per the SAFS-dynamic convention):
///   * RS nucleation length  L_nuc = μ · Dc / ((b − a) · σ_n_eff)  — the direct
///     rate-and-state analog of the LSW `μ · d_c / ((μ_s − μ_d) · σ_n)` formula
///     (raw, no Day/Andrews prefactor), evaluated only on velocity-WEAKENING
///     (b > a) DOFs.
///   * Steady-state friction  f_ss(V) = a · asinh[(V / 2V_0) · exp(ψ_ss / a)]
///     with ψ_ss = f_0 + b · ln(V_0 / V)  — the repo's regularized Dieterich-
///     Ruina coefficient (friction/dieterich_ruina.hpp) at the seeded
///     steady-state ψ (friction/state_evolution.hpp:196), evaluated at V_init.
///
/// @returns the maximum `|τ_pre| / (f_ss(V_init) · σ_n_eff)` ratio observed at
/// any DOF OUTSIDE the nucleation patch (the RS analog of the LSW
/// `|τ_pre| / (μ_s · σ_n)` return); 0.0 when `cfg.enabled == false`.
///
/// Gates (abort iff `cfg.abort_on_failure` and not `SEAS_SKIP_EQUILIBRIUM_GATE=1`):
///   (a) every DOF is velocity-STRENGTHENING (b ≤ a) — nothing can nucleate;
///   (b) nucleation enabled but ZERO fault DOFs inside the patch;
///   (c) nucleation enabled but the in-patch region is entirely
///       velocity-strengthening (b ≤ a) — the patch cannot host an instability.
/// Resolution (L_nuc/h_min) and overstress-drive checks are WARN-only.
real_t PrintDerivedAndCheckRS(
   const PrintDerivedConfig&                cfg,
   const RateStatePerDOFParams&             rs,
   const Vector&                            tau_pre_per_dof,
   const Vector&                            sigma_n_eff_per_dof,
   const Vector&                            dof_coords_3d,
   const NucleationSpec&                    nuc,
   const GradualOverstressPerDOFParams&     nuc_params,
   const StressSpec&                        stress,
   real_t                                   mu_bulk,
   real_t                                   cp,
   real_t                                   cs,
   real_t                                   h_min_global,
   real_t                                   dt_cfl,
   real_t                                   tfinal,
   int                                      num_fault_global,
   int                                      num_zero_normal_fallbacks
#ifdef MFEM_USE_MPI
   , MPI_Comm comm
#endif
   , int                                    rank,
   std::ostream&                            out = std::cout);

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_SPATIAL_PRINT_DERIVED_HPP
