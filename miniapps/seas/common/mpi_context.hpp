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

namespace mfem
{
namespace seas
{

/// @brief RAII wrapper for MPI initialization
///
/// In serial mode (when SEAS_USE_MPI is not defined), this class
/// provides stub implementations that do nothing, allowing the same
/// code to compile and run correctly in serial mode.
class MPIContext
{
public:
#ifdef SEAS_USE_MPI
   /// Construct and initialize MPI
   MPIContext(int *argc, char ***argv)
   {
      MPI_Init(argc, argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
      MPI_Comm_size(MPI_COMM_WORLD, &size_);
   }

   /// Finalize MPI
   ~MPIContext() { MPI_Finalize(); }

   /// Get MPI communicator
   MPI_Comm GetComm() const { return MPI_COMM_WORLD; }
#else
   /// Serial mode constructor (no-op)
   MPIContext(int *argc, char ***argv) : rank_(0), size_(1) {}

   /// Serial mode destructor (no-op)
   ~MPIContext() = default;
#endif

   /// Get this process's rank
   int Rank() const { return rank_; }

   /// Get total number of processes
   int Size() const { return size_; }

   /// Check if this is the root process (rank 0)
   bool IsRoot() const { return rank_ == 0; }

   /// @brief Global reduction for maximum value
   ///
   /// In serial mode, returns the input value unchanged.
   /// In parallel mode, performs MPI_Allreduce with MPI_MAX.
   real_t GlobalMax(real_t local_val) const
   {
#ifdef SEAS_USE_MPI
      real_t global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
      return global_val;
#else
      return local_val;
#endif
   }

   /// @brief Global reduction for minimum value
   ///
   /// In serial mode, returns the input value unchanged.
   /// In parallel mode, performs MPI_Allreduce with MPI_MIN.
   real_t GlobalMin(real_t local_val) const
   {
#ifdef SEAS_USE_MPI
      real_t global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MIN,
                    MPI_COMM_WORLD);
      return global_val;
#else
      return local_val;
#endif
   }

   /// @brief Global reduction for sum
   ///
   /// In serial mode, returns the input value unchanged.
   /// In parallel mode, performs MPI_Allreduce with MPI_SUM.
   real_t GlobalSum(real_t local_val) const
   {
#ifdef SEAS_USE_MPI
      real_t global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      return global_val;
#else
      return local_val;
#endif
   }

   /// @brief Global reduction for integer sum
   ///
   /// In serial mode, returns the input value unchanged.
   /// In parallel mode, performs MPI_Allreduce with MPI_SUM.
   int GlobalSumInt(int local_val) const
   {
#ifdef SEAS_USE_MPI
      int global_val;
      MPI_Allreduce(&local_val, &global_val, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      return global_val;
#else
      return local_val;
#endif
   }

   /// @brief Barrier synchronization
   ///
   /// In serial mode, this is a no-op.
   /// In parallel mode, performs MPI_Barrier.
   void Barrier() const
   {
#ifdef SEAS_USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
   }

   /// @brief Broadcast a value from root to all processes
   ///
   /// In serial mode, this is a no-op.
   /// In parallel mode, performs MPI_Bcast from rank 0.
   void Bcast(real_t &value) const
   {
#ifdef SEAS_USE_MPI
      MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
   }

   /// @brief Broadcast a vector from root to all processes
   ///
   /// In serial mode, this is a no-op.
   /// In parallel mode, performs MPI_Bcast from rank 0.
   void Bcast(Vector &vec) const
   {
#ifdef SEAS_USE_MPI
      int n = vec.Size();
      MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (!IsRoot())
      {
         vec.SetSize(n);
      }
      MPI_Bcast(vec.GetData(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
   }

private:
   int rank_;
   int size_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_MPI_CONTEXT_HPP
