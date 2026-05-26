// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/fault_state_channel.hpp — single source of truth for the ParaView
// fault "state" channel value (R-005 / R-024).
//
// The spatial_dyn_driver ParaView snapshot writer and the R-005 unit test
// MUST agree on what the per-DOF "state" scalar is:
//   * LSW : the slip-weakening friction coefficient
//           LSWFrictionCoefficient_TPV205(|slip|, mu_s, mu_d, d_c).
//   * RS  : the rate-and-state state variable psi (InitializeFaultDOFs_Spatial_RS
//           never sets lsw_*, so the LSW formula would emit 0/NaN for an RS run).
// Extracting it here (R-024) lets the test exercise the real code path instead
// of mirroring the driver's inline branch.

#ifndef MFEM_SEAS_FAULT_STATE_CHANNEL_HPP
#define MFEM_SEAS_FAULT_STATE_CHANNEL_HPP

#include "mfem.hpp"

#include "fault_face_flux.hpp"   // DOFData
#include "tpv205_friction.hpp"   // LSWFrictionCoefficient_TPV205

#include <cmath>

namespace mfem
{
namespace seas
{

/// Per-DOF ParaView "state" channel value, selected by friction law.
/// @param is_lsw  true → LSW friction coefficient; false → RS state var psi.
/// @param d       the fault DOF record (reads slip1/slip2 + lsw_* for LSW, psi for RS).
inline real_t FaultStateChannelValue(bool is_lsw, const DOFData &d)
{
   if (is_lsw)
   {
      const real_t delta_norm = std::sqrt(d.slip1 * d.slip1 + d.slip2 * d.slip2);
      return LSWFrictionCoefficient_TPV205(delta_norm, d.lsw_mu_s,
                                           d.lsw_mu_d, d.lsw_d_c);
   }
   return d.psi;   // RS state variable (matches what tau*_corr was solved with)
}

}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_FAULT_STATE_CHANNEL_HPP
