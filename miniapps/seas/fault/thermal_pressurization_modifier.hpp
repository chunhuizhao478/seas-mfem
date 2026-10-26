// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// thermal_pressurization_modifier.hpp — 2D in-plane design rev 3, phase D
// (section 5.9): thermal pressurization (fault/thermal_pressurization.hpp)
// behind the Nitsche operator's normal-stress hook
// (fault/fault_normal_stress_modifier.hpp).  The law sees
// sigma_eff - (p(z = 0) - p_ini); the module is advanced with the step's
// tau V after the closures.  Point j of the module is owned fault DOF j of
// the operator.
//
// Scope (Chunhui Zhao, 2026-10-26).

#ifndef MFEM_SEAS_THERMAL_PRESSURIZATION_MODIFIER_HPP
#define MFEM_SEAS_THERMAL_PRESSURIZATION_MODIFIER_HPP

#include "fault_normal_stress_modifier.hpp"
#include "thermal_pressurization.hpp"

namespace mfem
{
namespace seas
{

class ThermalPressurizationModifier : public IFaultNormalStressModifier
{
public:
   /// `tp` must outlive this object.
   explicit ThermalPressurizationModifier(ThermalPressurization &tp) : tp_(tp) {}

   int NumPoints() const override { return tp_.NumPoints(); }
   real_t NormalStressReduction(int j) const override { return tp_.PorePressureChange(j); }
   void Advance(real_t dt, const real_t *tau, const real_t *V) override
   {
      tp_.Advance(dt, tau, V);
   }

   ThermalPressurization &Module() { return tp_; }
   const ThermalPressurization &Module() const { return tp_; }

private:
   ThermalPressurization &tp_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_THERMAL_PRESSURIZATION_MODIFIER_HPP
