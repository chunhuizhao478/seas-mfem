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

#ifndef MFEM_SEAS_MPI_TAGS_HPP
#define MFEM_SEAS_MPI_TAGS_HPP

namespace mfem
{
namespace seas
{

/// Named MPI tag constants for all point-to-point communication.
/// Centralised here so tag collisions are easy to spot.
namespace MPITag
{
   /// Fault data scatter: owned → ghost DOF exchange in ExpandOwnedToLocalFault
   constexpr int kFaultScatter = 27183;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_MPI_TAGS_HPP
