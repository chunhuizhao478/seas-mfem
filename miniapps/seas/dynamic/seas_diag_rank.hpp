// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// SEAS_DIAG_FAULT_FLUX rank cache: used by the C-1 / C-2 / C-3 / C-4
// bisection checkpoints (TPV102 v9.0.0 debug plan §0.5) to tag per-rank
// DIAG stderr without threading `rank` through FaultFaceFlux / GodunovFlux
// APIs.
//
// Definition lives in drivers/tpv102_driver.cpp; the value is written
// once at MPI init and is read-only thereafter.  Entire file is compiled
// out when SEAS_DIAG_FAULT_FLUX is not defined so production builds are
// byte-identical to pre-change builds.

#ifndef MFEM_SEAS_DIAG_RANK_HPP
#define MFEM_SEAS_DIAG_RANK_HPP

#ifdef SEAS_DIAG_FAULT_FLUX
namespace mfem
{
namespace seas
{
extern int g_seas_my_rank;
} // namespace seas
} // namespace mfem
#endif // SEAS_DIAG_FAULT_FLUX

#endif // MFEM_SEAS_DIAG_RANK_HPP
