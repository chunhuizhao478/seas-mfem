// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_BP1_PARAMS_HPP
#define MFEM_SEAS_BP1_PARAMS_HPP

#include "bp2_params.hpp"
#include <vector>

namespace mfem
{
namespace seas
{

/// Create a BP2Params struct with BP1-specific parameter overrides.
///
/// BP1 differs from BP2 only in:
/// - Dc = 0.008 m (vs 0.004 m for BP2)
/// - t_final = 3000 years (vs 1200 years for BP2)
///
/// All other parameters (mu, rho, cs, sigma_n, a0, amax, b, H, h, Wf,
/// f0, V0, Vp) are identical between BP1 and BP2.
///
/// Reference: https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf
inline BP2Params MakeBP1Params()
{
   BP2Params params;
   params.Dc = 0.008;  // BP1: 8 mm (vs BP2: 4 mm)
   params.t_final = 3000.0 * BP2Params::seconds_per_year;  // 3000 years
   return params;
}

/// Get the 15 BP1 probe depths in meters.
///
/// BP1 uses 2.5 km spacing: 0, 2.5, 5.0, 7.5, ..., 35.0 km
/// (15 probes total, vs BP2's 12 irregular probes)
///
/// Convention: z=0 at surface, z<0 at depth.
inline std::vector<double> GetBP1ProbeDepths()
{
   return {
      0.0, -2500.0, -5000.0, -7500.0, -10000.0,
      -12500.0, -15000.0, -17500.0, -20000.0,
      -22500.0, -25000.0, -27500.0, -30000.0,
      -32500.0, -35000.0
   };
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP1_PARAMS_HPP
