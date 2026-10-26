// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_normal_stress_modifier.hpp — 2D in-plane design rev 3, phase D
// (section 5.9): the interface of the Nitsche operator's normal-stress hook
// (NitscheWaveOperator::SetNormalStressModifier).
//
// A modifier owns a state per owned fault DOF that lowers the normal stress
// the friction law sees.  In box 9 of the fault step the closure at owned
// DOF j sees
//
//    sigma_law = sigma_eff - NormalStressReduction(j)
//
// (sigma_eff = the Prakash-Clifton sigma* when that layer is on), with the
// reduction from the end of the previous step; after the closures of every
// owned DOF the operator calls Advance once with the step's dt, the shear
// traction |tau_j| the law balanced (|X_tan + gamma V|, the station traction)
// and the slip rate |V_j^{n+1}|, so the modifier lags one step, like the
// Prakash-Clifton and Kelvin-Voigt layers.  Thermal pressurization
// (fault/thermal_pressurization_modifier.hpp, SCEC TPV105) is the first one:
// the reduction is the pore-pressure change on the fault, p(z = 0) - p_ini.
//
// Scope (Chunhui Zhao, 2026-10-26): owned fault DOFs in the operator's order;
// the modifier never changes the traction the bulk sees, only the strength.

#ifndef MFEM_SEAS_FAULT_NORMAL_STRESS_MODIFIER_HPP
#define MFEM_SEAS_FAULT_NORMAL_STRESS_MODIFIER_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

class IFaultNormalStressModifier
{
public:
   virtual ~IFaultNormalStressModifier() = default;

   /// Owned fault DOFs it covers (the operator checks this against its own).
   virtual int NumPoints() const = 0;

   /// Amount taken off sigma_n_eff (compression positive, Pa) at owned DOF j
   /// before the closure.
   virtual real_t NormalStressReduction(int j) const = 0;

   /// One fault step of length dt (s): tau[j] the shear traction the law
   /// balanced (Pa, >= 0) and V[j] = |V^{n+1}| (m/s) per owned DOF.
   virtual void Advance(real_t dt, const real_t *tau, const real_t *V) = 0;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_NORMAL_STRESS_MODIFIER_HPP
