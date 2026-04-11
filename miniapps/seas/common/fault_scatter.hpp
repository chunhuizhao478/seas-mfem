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

#ifndef MFEM_SEAS_FAULT_SCATTER_HPP
#define MFEM_SEAS_FAULT_SCATTER_HPP

#include "mfem.hpp"
#include "mpi_check.hpp"
#include "mpi_tags.hpp"
#include <vector>

namespace mfem
{
namespace seas
{

/// Communication block for one neighbor rank in the fault scatter.
/// Each block lists which owned faces to send and which local (ghost)
/// faces to receive from that neighbor.
struct SharedFaultCommBlock
{
   int neighbor_rank = -1;
   std::vector<int> send_owned_faces;
   std::vector<int> recv_local_faces;
};

/// @brief Encapsulates MPI point-to-point communication for fault DOF scatter.
///
/// Replaces the raw MPI_Irecv/MPI_Isend/MPI_Waitall pattern in
/// ExpandOwnedToLocalFault with a reusable, testable abstraction.
///
/// Workflow:
///   1. Construct with comm blocks and communicator
///   2. BeginScatter(owned_data, comps_per_dof) — posts non-blocking sends/recvs
///   3. WaitScatter() — completes communication, returns received buffers
///
/// The caller is responsible for unpacking received data (applying permutations
/// etc.), since that logic depends on the domain operator's DOF numbering.
class FaultScatter
{
public:
   // Non-copyable, non-movable (owns MPI_Request handles, stores reference)
   FaultScatter(const FaultScatter &) = delete;
   FaultScatter &operator=(const FaultScatter &) = delete;

#ifdef SEAS_USE_MPI
   FaultScatter(const std::vector<SharedFaultCommBlock> &blocks,
                MPI_Comm comm)
      : blocks_(blocks), comm_(comm),
        send_bufs_(blocks.size()), recv_bufs_(blocks.size()) {}

   /// @brief Post non-blocking receives and sends for all comm blocks.
   ///
   /// @param owned_data     Owned fault data vector (contiguous, interleaved)
   /// @param comps_per_dof  Number of components per DOF
   /// @param nbf_per_face   Number of basis functions per face
   ///
   /// The send buffer is packed: for each send face j, for each basis function
   /// kk in [0, nbf_per_face), for each component c in [0, comps_per_dof):
   ///   send_buf[j * block_size + kk * comps_per_dof + c]
   ///     = owned_data[comps_per_dof * (send_face * nbf_per_face + kk) + c]
   void BeginScatter(const Vector &owned_data,
                     int comps_per_dof,
                     int nbf_per_face)
   {
      requests_.clear();
      requests_.reserve(2 * blocks_.size());

      const int block_size = comps_per_dof * nbf_per_face;

      for (int bi = 0; bi < static_cast<int>(blocks_.size()); bi++)
      {
         const auto &block = blocks_[bi];

         // Post receives
         if (!block.recv_local_faces.empty())
         {
            recv_bufs_[bi].SetSize(
               block_size * static_cast<int>(block.recv_local_faces.size()));
            MPI_Request req;
            MFEM_SEAS_MPI_CHECK(
               MPI_Irecv(recv_bufs_[bi].GetData(), recv_bufs_[bi].Size(),
                         MPI_DOUBLE, block.neighbor_rank,
                         MPITag::kFaultScatter, comm_, &req));
            requests_.push_back(req);
         }

         // Pack and send
         if (!block.send_owned_faces.empty())
         {
            send_bufs_[bi].SetSize(
               block_size * static_cast<int>(block.send_owned_faces.size()));
            for (int j = 0;
                 j < static_cast<int>(block.send_owned_faces.size()); j++)
            {
               const int owned_face = block.send_owned_faces[j];
               for (int kk = 0; kk < nbf_per_face; kk++)
               {
                  const int owned_dof = owned_face * nbf_per_face + kk;
                  for (int c = 0; c < comps_per_dof; c++)
                  {
                     send_bufs_[bi](j * block_size + kk * comps_per_dof + c) =
                        owned_data(comps_per_dof * owned_dof + c);
                  }
               }
            }

            MPI_Request req;
            MFEM_SEAS_MPI_CHECK(
               MPI_Isend(send_bufs_[bi].GetData(), send_bufs_[bi].Size(),
                         MPI_DOUBLE, block.neighbor_rank,
                         MPITag::kFaultScatter, comm_, &req));
            requests_.push_back(req);
         }
      }
   }

   /// @brief Non-blocking test: returns true if all requests are complete.
   bool TestScatter()
   {
      if (requests_.empty()) { return true; }
      int flag = 0;
      MFEM_SEAS_MPI_CHECK(
         MPI_Testall(static_cast<int>(requests_.size()), requests_.data(),
                     &flag, MPI_STATUSES_IGNORE));
      return flag != 0;
   }

   /// @brief Block until all sends/receives complete.
   ///
   /// After this returns, recv_bufs_ contain the received data.
   /// Access via GetRecvBuffer(block_index).
   void WaitScatter()
   {
      if (!requests_.empty())
      {
         MFEM_SEAS_MPI_CHECK(
            MPI_Waitall(static_cast<int>(requests_.size()),
                        requests_.data(), MPI_STATUSES_IGNORE));
      }
   }

   /// Access the receive buffer for a given comm block after WaitScatter().
   const Vector &GetRecvBuffer(int block_idx) const
   {
      return recv_bufs_[block_idx];
   }

   /// Number of comm blocks.
   int NumBlocks() const { return static_cast<int>(blocks_.size()); }

   /// Access comm block metadata.
   const SharedFaultCommBlock &GetBlock(int idx) const { return blocks_[idx]; }

#else
   // Serial stub: no communication needed.
   FaultScatter(const std::vector<SharedFaultCommBlock> &blocks) {}

   void BeginScatter(const Vector &, int, int) {}
   bool TestScatter() { return true; }
   void WaitScatter() {}
   const Vector &GetRecvBuffer(int block_idx) const
   {
      static Vector empty;
      return empty;
   }
   int NumBlocks() const { return 0; }
   const SharedFaultCommBlock &GetBlock(int) const
   {
      MFEM_ABORT("No blocks in serial mode");
      static SharedFaultCommBlock dummy;
      return dummy;
   }
#endif

private:
#ifdef SEAS_USE_MPI
   const std::vector<SharedFaultCommBlock> &blocks_;
   MPI_Comm comm_;
   std::vector<Vector> send_bufs_;
   std::vector<Vector> recv_bufs_;
   std::vector<MPI_Request> requests_;
#endif
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_SCATTER_HPP
