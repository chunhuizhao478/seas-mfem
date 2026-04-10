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

#ifndef MFEM_SEAS_MPI_CONTEXT_HPP
#define MFEM_SEAS_MPI_CONTEXT_HPP

#include "mfem.hpp"
#include "mpi_check.hpp"
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief RAII wrapper for MPI initialization with stored communicator
///
/// In serial mode (when SEAS_USE_MPI is not defined), this class
/// provides stub implementations that do nothing, allowing the same
/// code to compile and run correctly in serial mode.
///
/// All methods use the stored communicator (comm_), not MPI_COMM_WORLD.
class MPIContext
{
public:
#ifdef SEAS_USE_MPI
   /// Construct with MPI_COMM_WORLD (default).
   MPIContext(int *argc, char ***argv)
   {
      if (!Mpi::IsInitialized())
      {
         Mpi::Init(*argc, *argv);
      }
      Hypre::Init();
      comm_ = MPI_COMM_WORLD;
      MPI_Comm_rank(comm_, &rank_);
      MPI_Comm_size(comm_, &size_);
   }

   /// Construct with a specific communicator (duplicated internally).
   MPIContext(MPI_Comm comm)
   {
      MPI_Comm_dup(comm, &comm_);
      owns_comm_ = true;
      MPI_Comm_rank(comm_, &rank_);
      MPI_Comm_size(comm_, &size_);
   }

   ~MPIContext()
   {
      if (owns_comm_)
      {
         MPI_Comm_free(&comm_);
      }
   }

   /// Get stored MPI communicator
   MPI_Comm GetComm() const { return comm_; }
#else
   /// Serial mode constructor (no-op)
   MPIContext(int *argc, char ***argv) : rank_(0), size_(1) {}

   /// Serial mode destructor (no-op)
   ~MPIContext() = default;
#endif

   // Non-copyable
   MPIContext(const MPIContext &) = delete;
   MPIContext &operator=(const MPIContext &) = delete;

   /// Get this process's rank
   int Rank() const { return rank_; }

   /// Get total number of processes
   int Size() const { return size_; }

   /// Check if this is the root process (rank 0)
   bool IsRoot() const { return rank_ == 0; }

   // =========================================================================
   // Reductions
   // =========================================================================

   real_t GlobalMax(real_t local_val) const
   {
#ifdef SEAS_USE_MPI
      real_t global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MAX, comm_);
      return global_val;
#else
      return local_val;
#endif
   }

   real_t GlobalMin(real_t local_val) const
   {
#ifdef SEAS_USE_MPI
      real_t global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MIN, comm_);
      return global_val;
#else
      return local_val;
#endif
   }

   real_t GlobalSum(real_t local_val) const
   {
#ifdef SEAS_USE_MPI
      real_t global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, comm_);
      return global_val;
#else
      return local_val;
#endif
   }

   int GlobalMinInt(int local_val) const
   {
#ifdef SEAS_USE_MPI
      int global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_INT, MPI_MIN, comm_);
      return global_val;
#else
      return local_val;
#endif
   }

   int GlobalMaxInt(int local_val) const
   {
#ifdef SEAS_USE_MPI
      int global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_INT, MPI_MAX, comm_);
      return global_val;
#else
      return local_val;
#endif
   }

   int GlobalSumInt(int local_val) const
   {
#ifdef SEAS_USE_MPI
      int global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_INT, MPI_SUM, comm_);
      return global_val;
#else
      return local_val;
#endif
   }

   // =========================================================================
   // Synchronization
   // =========================================================================

   void Barrier() const
   {
#ifdef SEAS_USE_MPI
      MPI_Barrier(comm_);
#endif
   }

   // =========================================================================
   // Broadcast
   // =========================================================================

   void Bcast(int &value) const
   {
#ifdef SEAS_USE_MPI
      MPI_Bcast(&value, 1, MPI_INT, 0, comm_);
#endif
   }

   void Bcast(real_t &value) const
   {
#ifdef SEAS_USE_MPI
      MPI_Bcast(&value, 1, MPI_DOUBLE, 0, comm_);
#endif
   }

   void Bcast(Vector &vec) const
   {
#ifdef SEAS_USE_MPI
      int n = vec.Size();
      MPI_Bcast(&n, 1, MPI_INT, 0, comm_);
      if (!IsRoot())
      {
         vec.SetSize(n);
      }
      MPI_Bcast(vec.GetData(), n, MPI_DOUBLE, 0, comm_);
#endif
   }

   // =========================================================================
   // Gather / Allgather
   // =========================================================================

   /// @brief Gather local vectors from all ranks to root (rank 0).
   ///
   /// Each rank provides a local MFEM Vector. On root, the output contains
   /// the concatenation of all local vectors in rank order.
   /// On non-root ranks, the output vector is resized to 0.
   void GatherToRoot(const Vector &local, Vector &global) const
   {
#ifdef SEAS_USE_MPI
      int local_count = local.Size();

      // Gather counts from all ranks to root
      std::vector<int> counts(size_);
      MFEM_SEAS_MPI_CHECK(
         MPI_Gather(&local_count, 1, MPI_INT,
                    counts.data(), 1, MPI_INT, 0, comm_));

      // Compute displacements on root
      std::vector<int> displs(size_, 0);
      if (IsRoot())
      {
         for (int i = 1; i < size_; i++)
         {
            displs[i] = displs[i - 1] + counts[i - 1];
         }
         int total = displs[size_ - 1] + counts[size_ - 1];
         MFEM_ASSERT(total >= 0, "GatherToRoot: overflow in total count");
         global.SetSize(total);
      }
      else
      {
         global.SetSize(0);
      }

      MFEM_SEAS_MPI_CHECK(
         MPI_Gatherv(local.GetData(), local_count, MPI_DOUBLE,
                     IsRoot() ? global.GetData() : nullptr,
                     counts.data(), displs.data(), MPI_DOUBLE, 0, comm_));
#else
      global = local;
#endif
   }

   /// @brief Allgather: concatenate local vectors from all ranks.
   ///
   /// Each rank provides a local MFEM Vector. On return, every rank
   /// holds the concatenation of all local vectors in rank order.
   void AllgatherVec(const Vector &local, Vector &global) const
   {
#ifdef SEAS_USE_MPI
      int local_count = local.Size();

      // Allgather counts
      std::vector<int> counts(size_);
      MFEM_SEAS_MPI_CHECK(
         MPI_Allgather(&local_count, 1, MPI_INT,
                       counts.data(), 1, MPI_INT, comm_));

      // Compute displacements
      std::vector<int> displs(size_, 0);
      for (int i = 1; i < size_; i++)
      {
         displs[i] = displs[i - 1] + counts[i - 1];
      }
      int total = displs[size_ - 1] + counts[size_ - 1];
      MFEM_ASSERT(total >= 0, "AllgatherVec: overflow in total count");
      global.SetSize(total);

      MFEM_SEAS_MPI_CHECK(
         MPI_Allgatherv(local.GetData(), local_count, MPI_DOUBLE,
                        global.GetData(), counts.data(), displs.data(),
                        MPI_DOUBLE, comm_));
#else
      global = local;
#endif
   }

private:
   int rank_;
   int size_;
#ifdef SEAS_USE_MPI
   MPI_Comm comm_;
   bool owns_comm_ = false;
#endif
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_MPI_CONTEXT_HPP
