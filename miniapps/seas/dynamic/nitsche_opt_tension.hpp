// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_opt_tension.hpp — option O2 of the complex-geometries plan
// (rev 2, Section 5.3; Chunhui Zhao 2026-10-20): a tension floor, for TPV12
// (and TPV13/15/22 should their normal stress turn tensile).
//
// Deck (Nitsche driver only; absent, or mode = "abort", = every bit as before):
//
//   [fault_tension]
//   mode = "zero_strength_floor"     # or "abort" (today's behaviour)
//
// With "zero_strength_floor" the friction law sees max(sigma, 0), sigma the
// compression-positive normal stress after Prakash-Clifton, so the shear
// strength in tension is the cohesion, as the TPV12, TPV15 and TPV22
// descriptions define it.  sigma* (Prakash-Clifton) and the recorded
// sigma_n keep the true value; NaN still aborts.  The fault does not open:
// the slip stays tangential and the Nitsche penalty alone holds the normal
// jump (about sigma / sigma_F).  NitscheWaveOperator::SetTensionFloor does
// the work; the run log counts the floored DOF-steps.
//
// Refused with a non-zero [pore_pressure] min_sigma_n_pa (the floor replaces
// the zero range check only) and with rate and state (slip weakening only).

#ifndef MFEM_SEAS_NITSCHE_OPT_TENSION_HPP
#define MFEM_SEAS_NITSCHE_OPT_TENSION_HPP

#include "mfem.hpp"
#include "nitsche_opt_common.hpp"

#include <string>

namespace mfem
{
namespace seas
{
namespace nitsche_opt
{

struct TensionConfig
{
   static constexpr const char *kTable = "fault_tension";
   bool floor = false;   ///< mode = "zero_strength_floor"
};

inline TensionConfig ParseTension(const toml::value &t)
{
   using T = TensionConfig;
   OptRefuseUnknownKeys(t, T::kTable, {"mode"});
   const std::string m = OptStr(t, T::kTable, "mode");
   MFEM_VERIFY(m == "zero_strength_floor" || m == "abort", "[" << T::kTable << "] mode = '"
               << m << "': expected \"zero_strength_floor\" or \"abort\"");
   T c;
   c.floor = (m == "zero_strength_floor");
   return c;
}

} // namespace nitsche_opt
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_OPT_TENSION_HPP
