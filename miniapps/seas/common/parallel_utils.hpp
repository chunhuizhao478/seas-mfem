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

#ifndef MFEM_SEAS_PARALLEL_UTILS_HPP
#define MFEM_SEAS_PARALLEL_UTILS_HPP

#include "mfem.hpp"
#include "mpi_check.hpp"
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef SEAS_USE_MPI

/// @brief Gather local vectors from all ranks to root (rank 0)
///
/// Each rank provides a local vector. On root, the output contains
/// the concatenation of all local vectors in rank order.
/// On non-root ranks, the output vector is empty.
///
/// Safety: MPI_Gather and MPI_Gatherv are collective — all ranks must call
/// this function or deadlock results.
inline void GatherVectorToRoot(const std::vector<real_t> &local,
                               std::vector<real_t> &global,
                               MPI_Comm comm = MPI_COMM_WORLD)
{
   int rank, size;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &size);

   int local_count = static_cast<int>(local.size());
   MFEM_ASSERT(local_count >= 0, "GatherVectorToRoot: negative local count");

   // Gather counts from all ranks to root (collective)
   std::vector<int> counts(size);
   MFEM_SEAS_MPI_CHECK(
      MPI_Gather(&local_count, 1, MPI_INT,
                 counts.data(), 1, MPI_INT, 0, comm));

   // Compute displacements on root
   std::vector<int> displs(size, 0);
   if (rank == 0)
   {
      for (int i = 1; i < size; i++)
      {
         displs[i] = displs[i-1] + counts[i-1];
      }
      int total = displs[size-1] + counts[size-1];
      MFEM_ASSERT(total >= 0, "GatherVectorToRoot: overflow in total count");
      global.resize(total);
   }

   MFEM_SEAS_MPI_CHECK(
      MPI_Gatherv(local.data(), local_count, MPI_DOUBLE,
                  rank == 0 ? global.data() : nullptr,
                  counts.data(), displs.data(), MPI_DOUBLE, 0, comm));
}

/// @brief Scatter a global vector from root to all ranks
///
/// Root provides the global vector and the per-rank counts.
/// Each rank receives its slice of the global vector.
///
/// Safety: MPI_Bcast and MPI_Scatterv are collective — all ranks must call
/// this function or deadlock results.
inline void ScatterVectorFromRoot(const std::vector<real_t> &global,
                                  const std::vector<int> &counts,
                                  std::vector<real_t> &local,
                                  MPI_Comm comm = MPI_COMM_WORLD)
{
   int rank, size;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &size);

   MFEM_ASSERT(rank != 0 || static_cast<int>(counts.size()) == size,
               "ScatterVectorFromRoot: counts.size() must equal comm size on root");

   // Broadcast counts to all ranks so each knows its local size (collective)
   std::vector<int> all_counts(size);
   if (rank == 0)
   {
      all_counts = counts;
   }
   MFEM_SEAS_MPI_CHECK(
      MPI_Bcast(all_counts.data(), size, MPI_INT, 0, comm));

   MFEM_ASSERT(all_counts[rank] >= 0,
               "ScatterVectorFromRoot: negative count for rank " << rank);
   local.resize(all_counts[rank]);

   // Compute displacements on root
   std::vector<int> displs(size, 0);
   for (int i = 1; i < size; i++)
   {
      displs[i] = displs[i-1] + all_counts[i-1];
   }

   MFEM_SEAS_MPI_CHECK(
      MPI_Scatterv(rank == 0 ? global.data() : nullptr,
                   all_counts.data(), displs.data(), MPI_DOUBLE,
                   local.data(), all_counts[rank], MPI_DOUBLE, 0, comm));
}

/// @brief Broadcast a value from root to all ranks
template <typename T>
inline void BroadcastFromRoot(T &value, MPI_Comm comm = MPI_COMM_WORLD);

template <>
inline void BroadcastFromRoot<int>(int &value, MPI_Comm comm)
{
   MPI_Bcast(&value, 1, MPI_INT, 0, comm);
}

template <>
inline void BroadcastFromRoot<real_t>(real_t &value, MPI_Comm comm)
{
   MPI_Bcast(&value, 1, MPI_DOUBLE, 0, comm);
}

/// @brief Debug helper: check if all ranks agree on a value
///
/// Returns true if min == max across all ranks (i.e., all have the same value).
/// Safety: Uses two MPI_Allreduce calls — collective, all ranks must participate.
template <typename T>
inline bool AllRanksAgree(T local_val, MPI_Comm comm = MPI_COMM_WORLD);

template <>
inline bool AllRanksAgree<int>(int local_val, MPI_Comm comm)
{
   int global_min, global_max;
   MPI_Allreduce(&local_val, &global_min, 1, MPI_INT, MPI_MIN, comm);
   MPI_Allreduce(&local_val, &global_max, 1, MPI_INT, MPI_MAX, comm);
   return global_min == global_max;
}

template <>
inline bool AllRanksAgree<real_t>(real_t local_val, MPI_Comm comm)
{
   real_t global_min, global_max;
   MPI_Allreduce(&local_val, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm);
   MPI_Allreduce(&local_val, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
   return global_min == global_max;
}

#else // Serial stubs

inline void GatherVectorToRoot(const std::vector<real_t> &local,
                               std::vector<real_t> &global)
{
   global = local;
}

inline void ScatterVectorFromRoot(const std::vector<real_t> &global,
                                  const std::vector<int> &counts,
                                  std::vector<real_t> &local)
{
   local = global;
}

template <typename T>
inline void BroadcastFromRoot(T &value)
{
   // No-op in serial
}

template <typename T>
inline bool AllRanksAgree(T local_val)
{
   return true;
}

#endif // SEAS_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PARALLEL_UTILS_HPP
