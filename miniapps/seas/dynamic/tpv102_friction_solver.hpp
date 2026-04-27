// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV102 friction solver — Newton-Raphson on the regularized
// rate-and-state friction balance.
//
// The TPV102 friction COEFFICIENT
//   τ = a · σ · arcsinh[ V/(2V0) · exp(ψ/a) ]
// has the same functional form as TPV104's stable-asinh μ (`f_w` and
// `V_w` enter only through the slip-law-with-SRW STATE evolution, not
// through the coefficient).  Therefore the per-QP Newton solver from
// `tpv104_friction_solver.hpp` is bit-identically reusable for TPV102 —
// the only TPV-specific knob is `(a, V0)`, which the caller already
// supplies as scalar arguments.
//
// This header is a thin alias so the TPV102 driver's include surface
// matches TPV104's symmetrically.  No new code, no behavioural drift.

#ifndef MFEM_SEAS_TPV102_FRICTION_SOLVER_HPP
#define MFEM_SEAS_TPV102_FRICTION_SOLVER_HPP

#include "tpv104_friction_solver.hpp"

namespace mfem
{
namespace seas
{

inline constexpr real_t kTpv102AlmostZero = kTpv104AlmostZero;

/// Newton-Raphson slip-rate solver for TPV102 (aging-law rate-and-state).
/// Same μ + Newton machinery as TPV104; aliased here for symmetry.
inline real_t SolveSlipRateNewtonStable_TPV102(
   real_t tau_abs, real_t psi, real_t sigma_n, real_t eta_s,
   real_t a, real_t V0,
   real_t V_prev,
   int max_iter = 60,
   real_t tol   = 1e-8,
   int *iterations = nullptr,
   bool *has_converged = nullptr)
{
   return SolveSlipRateNewtonStable(tau_abs, psi, sigma_n, eta_s,
                                    a, V0, V_prev, max_iter, tol,
                                    iterations, has_converged);
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV102_FRICTION_SOLVER_HPP
