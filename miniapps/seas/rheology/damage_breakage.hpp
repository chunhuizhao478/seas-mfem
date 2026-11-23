// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// damage_breakage.hpp — continuum damage-breakage, phase R8 of the rheology
// plan (rev 1, Appendix D).  NOT IMPLEMENTED: ParseRheology refuses
// [rheology.damage_breakage] enabled = true until the R8 plan revision,
// which takes the equations from Chunhui Zhao's model and its references.
//
// Outline: the stress is a nonlinear function of the elastic strain with
// moduli that depend on a damage and a breakage variable, both with their
// own evolution laws.  In the IRheologyModel interface those variables are
// the state, Delta sigma = sigma_DB - C_0 : eps, and StiffnessBound must
// stay <= 1 (the model may only soften relative to C_0) or the coupler
// stops the run.  Combines only with attenuation (rheology_stack.hpp).
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_DAMAGE_BREAKAGE_HPP
#define MFEM_SEAS_DAMAGE_BREAKAGE_HPP

namespace mfem
{
namespace seas
{

/// Plan phase that implements this model (named in the parser's refusal).
constexpr const char *kDamageBreakagePhase = "R8";

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DAMAGE_BREAKAGE_HPP
