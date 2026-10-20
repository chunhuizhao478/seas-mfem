// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// viscoelastic.hpp — viscoelastic attenuation, phase R7 of the rheology plan
// (rev 1, Appendix D).  NOT IMPLEMENTED: ParseRheology refuses
// [rheology.viscoelastic] enabled = true until the R7 plan revision.
//
// Outline: a generalized Maxwell body with L mechanisms,
//    sigma = C_u : eps - sum_l zeta_l,   d zeta_l/dt + omega_l zeta_l = omega_l Y_l : eps,
// C_u the unrelaxed moduli (they become C_0 and K), omega_l and Y_l fitted
// to Q_p, Q_s over a frequency band (Emmerich and Korn 1987; Kaser et al.
// 2007).  In the IRheologyModel interface: 6 L state variables per node,
// Delta sigma = -sum_l zeta_l, an exact exponential update over dt, and
// IsLinearInStrain() = true.  First in the stack, before plasticity or
// damage-breakage (rheology_stack.hpp).
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_VISCOELASTIC_HPP
#define MFEM_SEAS_VISCOELASTIC_HPP

namespace mfem
{
namespace seas
{

/// Plan phase that implements this model (named in the parser's refusal).
constexpr const char *kViscoelasticPhase = "R7";

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_VISCOELASTIC_HPP
