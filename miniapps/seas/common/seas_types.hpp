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

#ifndef MFEM_SEAS_TYPES_HPP
#define MFEM_SEAS_TYPES_HPP

#include "mfem.hpp"
#include <type_traits>

namespace mfem
{
namespace seas
{

/// @brief Type aliases for serial/parallel MFEM class switching
///
/// When SEAS_USE_MPI is defined, parallel types are used.
/// Otherwise, serial types are used. This allows the same code
/// to compile for both serial and parallel execution.

#ifdef SEAS_USE_MPI

// Parallel type aliases
using SEASMesh = ParMesh;
using SEASFiniteElementSpace = ParFiniteElementSpace;
using SEASBilinearForm = ParBilinearForm;
using SEASLinearForm = ParLinearForm;
using SEASGridFunction = ParGridFunction;


/// Get MPI communicator
inline MPI_Comm GetSEASComm() { return MPI_COMM_WORLD; }

/// Get MPI rank
inline int GetSEASRank()
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   return rank;
}

/// Get MPI size
inline int GetSEASSize()
{
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   return size;
}

/// Check if this is the root process
inline bool IsSEASRoot() { return GetSEASRank() == 0; }

#else

// Serial type aliases
using SEASMesh = Mesh;
using SEASFiniteElementSpace = FiniteElementSpace;
using SEASBilinearForm = BilinearForm;
using SEASLinearForm = LinearForm;
using SEASGridFunction = GridFunction;


// Serial stubs for MPI functions
inline int GetSEASRank() { return 0; }
inline int GetSEASSize() { return 1; }
inline bool IsSEASRoot() { return true; }

#endif

/// @brief Type traits for compile-time mesh type detection
template <typename MeshType>
struct IsParallelMesh : std::false_type {};

#ifdef MFEM_USE_MPI
template <>
struct IsParallelMesh<ParMesh> : std::true_type {};
#endif

/// @brief Conditional type selection based on mesh type
template <typename MeshType, typename SerialType, typename ParallelType>
using MeshConditional = typename std::conditional<
   IsParallelMesh<MeshType>::value,
   ParallelType,
   SerialType>::type;

/// @brief Get the appropriate FE space type for a given mesh type
template <typename MeshType>
using FESpaceForMesh = MeshConditional<MeshType,
                                       FiniteElementSpace,
#ifdef MFEM_USE_MPI
                                       ParFiniteElementSpace>;
#else
                                       FiniteElementSpace>;
#endif

/// @brief Get the appropriate bilinear form type for a given mesh type
template <typename MeshType>
using BilinearFormForMesh = MeshConditional<MeshType,
                                            BilinearForm,
#ifdef MFEM_USE_MPI
                                            ParBilinearForm>;
#else
                                            BilinearForm>;
#endif

/// @brief Get the appropriate linear form type for a given mesh type
template <typename MeshType>
using LinearFormForMesh = MeshConditional<MeshType,
                                          LinearForm,
#ifdef MFEM_USE_MPI
                                          ParLinearForm>;
#else
                                          LinearForm>;
#endif

/// @brief Get the appropriate grid function type for a given mesh type
template <typename MeshType>
using GridFunctionForMesh = MeshConditional<MeshType,
                                            GridFunction,
#ifdef MFEM_USE_MPI
                                            ParGridFunction>;
#else
                                            GridFunction>;
#endif

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TYPES_HPP
