// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_HDF5_ERROR_FILTER_HPP
#define MFEM_SEAS_HDF5_ERROR_FILTER_HPP

// Suppress HDF5's auto-print of internal error stacks.  On the Frontera
// build, two HDF5 ABIs are loaded in the same process: HDF5 1.10
// (libhdf5.so.200, pulled in by PETSc 3.15 via libpetsc.so) and HDF5
// 1.14 (libhdf5.so.310, built from source for VTKHDF support).  Each
// has its own static ID table; cross-instance closes look up unknown
// IDs in the other library's table and trigger noisy but harmless
// "can't locate ID (already closed?)" diagnostic stacks on stderr.
//
// Real HDF5 failures (disk full, permission denied, missing dataset)
// are surfaced by seas's MFEM-wrapper code paths
// (ParaViewHDFDataCollection, FaultHDFState) via MFEM_VERIFY +
// exceptions, so suppressing the auto-print of the error stack does
// NOT hide actionable errors.  See
// debug_document/paraview_output_debug_document/hdf5_diag_noise_2026-05-17.md
// for the full root-cause analysis.

#ifdef MFEM_USE_HDF5
#include <hdf5.h>
#endif

namespace mfem
{
namespace seas
{

/// @brief Install the seas HDF5 error-output suppressor.
///
/// Call ONCE at driver startup after MPI_Init in every binary that
/// opens HDF5 files (BP5 / TPV104 drivers).  On builds without
/// MFEM_USE_HDF5 this is a compile-time no-op so the same call site
/// is safe in both HDF5 and non-HDF5 configurations.
inline void InstallHdf5ErrorFilter()
{
#ifdef MFEM_USE_HDF5
   // H5E_DEFAULT is the default error stack; passing (NULL, NULL)
   // disables HDF5's auto-print of error stacks to stderr.  Errors
   // are still recorded in the stack and accessible programmatically
   // via H5Eget_current_stack / H5Ewalk2 if a future code path wants
   // to inspect them.
   H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
#endif
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_HDF5_ERROR_FILTER_HPP
