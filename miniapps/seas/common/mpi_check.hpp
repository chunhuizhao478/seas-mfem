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

#ifndef MFEM_SEAS_MPI_CHECK_HPP
#define MFEM_SEAS_MPI_CHECK_HPP

#include "mfem.hpp"

/// @brief Assert that an MPI call returns MPI_SUCCESS (debug builds only).
///
/// In debug builds (MFEM_DEBUG defined), evaluates the MPI call and aborts
/// with a diagnostic message if the return code is not MPI_SUCCESS.
/// In release builds, evaluates the call but discards the return value.
///
/// Usage:
///   MFEM_SEAS_MPI_CHECK(MPI_Isend(buf, n, MPI_DOUBLE, dest, tag, comm, &req));
#ifdef SEAS_USE_MPI
#  ifdef MFEM_DEBUG
#    define MFEM_SEAS_MPI_CHECK(call)                                         \
      do {                                                                     \
         int _mpi_err_ = (call);                                               \
         MFEM_VERIFY(_mpi_err_ == MPI_SUCCESS,                                \
                     "MPI call failed with error code " << _mpi_err_           \
                     << " in " << __FILE__ << ":" << __LINE__);                \
      } while (0)
#  else
#    define MFEM_SEAS_MPI_CHECK(call) (void)(call)
#  endif
#else
   // Serial builds: no-op
#  define MFEM_SEAS_MPI_CHECK(call) ((void)0)
#endif

#endif // MFEM_SEAS_MPI_CHECK_HPP
