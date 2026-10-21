// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_opt_slip_path.hpp — option O3 of the complex-geometries plan
// (rev 2, Section 5.4; Chunhui Zhao 2026-10-20): slip weakening from the
// slip path, for TPV12, TPV15 and TPV22 ("total slip path distance").
//
// Deck (Nitsche driver only; absent, or measure = "slip_norm", = every bit
// as before):
//
//   [slip_weakening_distance]
//   measure = "slip_path"            # or "slip_norm" (today: |d^n|)
//
// With "slip_path" mu is computed from psi^n = int_0^{t_n} |V| dt, the slip
// path at the start of the step that the operator already integrates
// (trapezoid rule), instead of |d^n|.  The two agree for slip in a fixed
// direction up to that quadrature (the TPV205 gold ends with path
// 5.5524675 m against h-slip 5.5524646 m) and differ when the slip turns or
// reverses.  NitscheWaveOperator::SetSlipPathWeakening does the work, under
// global and local time stepping alike.  Slip weakening only: for
// rate-and-state psi is the state, not a slip path.

#ifndef MFEM_SEAS_NITSCHE_OPT_SLIP_PATH_HPP
#define MFEM_SEAS_NITSCHE_OPT_SLIP_PATH_HPP

#include "mfem.hpp"
#include "nitsche_opt_common.hpp"

#include <string>

namespace mfem
{
namespace seas
{
namespace nitsche_opt
{

struct SlipPathConfig
{
   static constexpr const char *kTable = "slip_weakening_distance";
   bool path = false;   ///< measure = "slip_path"
};

inline SlipPathConfig ParseSlipPath(const toml::value &t)
{
   using T = SlipPathConfig;
   OptRefuseUnknownKeys(t, T::kTable, {"measure"});
   const std::string m = OptStr(t, T::kTable, "measure");
   MFEM_VERIFY(m == "slip_path" || m == "slip_norm", "[" << T::kTable << "] measure = '"
               << m << "': expected \"slip_path\" or \"slip_norm\"");
   T c;
   c.path = (m == "slip_path");
   return c;
}

} // namespace nitsche_opt
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_OPT_SLIP_PATH_HPP
