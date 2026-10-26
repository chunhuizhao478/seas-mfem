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

// qd_node_slip.hpp — QD speed item B (qd-b-design-rev1): one slip vector
// per node for the Green's-function stage (--qd-green-gather node).
//
// QDGreenOperator::Apply needs all 2·N_f slips on every rank.  An
// MPI_Allgatherv over all ranks (--qd-green-gather rank) delivers one full
// copy to each rank, so 144 copies land on every Vista node at each stage.
// Here each node holds one copy in an MPI shared-memory window instead:
//   1. every rank writes its owned slip into the window at 2·own_lo,
//   2. node barrier,
//   3. with more than one node, node rank 0 of each node swaps the node's
//      part with the other nodes' rank 0 (MPI_Allgatherv over those ranks),
//   4. node barrier; every rank reads the window.
// The window then holds exactly what the all-gather would deliver, and no
// arithmetic is done, so the product is bitwise the same.
//
// Two windows alternate between calls.  A rank may leave a stage early
// (the traction return does not synchronize a node), so it can write the
// next stage's slip while another rank still reads this stage's.  Calls k
// and k + 2 share a window, and every rank passed call k + 1's first
// barrier only after it stopped reading call k's.  Hence the pointer that
// Gather returns stays valid until the caller's next Gather.
//
// Nodes come from MPI_Comm_split_type(MPI_COMM_TYPE_SHARED), or, as in
// QDHMatrix::SetupNodes, from SEAS_QD_HMAT_NODE_SIZE=n (each n consecutive
// ranks are one node), so a single machine runs the multi-node exchange.

#ifndef MFEM_SEAS_QD_NODE_SLIP_HPP
#define MFEM_SEAS_QD_NODE_SLIP_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

class QDNodeSlip
{
public:
   /// @a own_lo[r], @a own_n[r]: fault DOFs [own_lo, own_lo + own_n) that
   /// rank r owns, global indices in [0, @a nf); their slip is global
   /// [2·own_lo, 2·(own_lo + own_n)).  Collective.
   QDNodeSlip(MPI_Comm comm, const std::vector<long long> &own_lo,
              const std::vector<long long> &own_n, long long nf)
      : comm_(comm), own_lo_(own_lo), own_n_(own_n), nf_(nf)
   {
      MPI_Comm_size(comm_, &nranks_);
      MPI_Comm_rank(comm_, &rank_);
      MFEM_VERIFY(static_cast<int>(own_lo_.size()) == nranks_ &&
                  static_cast<int>(own_n_.size()) == nranks_,
                  "QD node slip: owned ranges of " << own_lo_.size()
                  << " ranks, expected " << nranks_);
      MFEM_VERIFY(2 * nf_ <= INT_MAX, "QD node slip: 2·N_f = " << 2 * nf_
                  << " exceeds the int counts of MPI_Allgatherv");
      SetupNodes();
      PlanExchange();
      const MPI_Aint bytes = node_rank_ == 0 ?
                             static_cast<MPI_Aint>(2 * nf_ * sizeof(real_t)) :
                             0;
      for (int b = 0; b < 2; b++)
      {
         real_t *mine = nullptr;
         MPI_Win_allocate_shared(bytes, sizeof(real_t), MPI_INFO_NULL,
                                 node_comm_, &mine, &win_[b]);
         MPI_Aint size = 0;
         int disp = 0;
         MPI_Win_shared_query(win_[b], 0, &size, &disp, &base_[b]);
         MFEM_VERIFY(size == static_cast<MPI_Aint>(2 * nf_ * sizeof(real_t)),
                     "QD node slip: shared window of " << size
                     << " B, expected " << 2 * nf_ * sizeof(real_t));
         MPI_Win_lock_all(MPI_MODE_NOCHECK, win_[b]);
         if (node_rank_ == 0) { std::fill(base_[b], base_[b] + 2 * nf_, 0.0); }
         MPI_Win_sync(win_[b]);
      }
      MPI_Barrier(node_comm_);
      for (int b = 0; b < 2; b++) { MPI_Win_sync(win_[b]); }
   }

   ~QDNodeSlip()
   {
      int done = 0;
      MPI_Finalized(&done);
      if (done) { return; }
      for (int b = 0; b < 2; b++)
      {
         if (win_[b] != MPI_WIN_NULL)
         {
            MPI_Win_unlock_all(win_[b]);
            MPI_Win_free(&win_[b]);
         }
      }
      if (lead_comm_ != MPI_COMM_NULL) { MPI_Comm_free(&lead_comm_); }
      if (node_comm_ != MPI_COMM_NULL) { MPI_Comm_free(&node_comm_); }
   }
   QDNodeSlip(const QDNodeSlip &) = delete;
   QDNodeSlip &operator=(const QDNodeSlip &) = delete;

   /// All 2·N_f slips in global order from this rank's @a slip_owned
   /// (2·own_n[rank] values).  The pointer stays valid until this rank's
   /// next Gather.  Collective; every rank calls it equally often.
   const real_t *Gather(const real_t *slip_owned)
   {
      const int b = static_cast<int>(calls_ & 1);
      calls_++;
      real_t *w = base_[b];
      if (own_n_[rank_] > 0)
      {
         std::copy(slip_owned, slip_owned + 2 * own_n_[rank_],
                   w + 2 * own_lo_[rank_]);
      }
      MPI_Win_sync(win_[b]);
      MPI_Barrier(node_comm_);
      MPI_Win_sync(win_[b]);
      if (nnodes_ > 1)
      {
         if (node_rank_ == 0) { Exchange(w); }
         MPI_Win_sync(win_[b]);
         MPI_Barrier(node_comm_);
         MPI_Win_sync(win_[b]);
      }
      return w;
   }

   int NumNodes() const { return nnodes_; }
   bool Contiguous() const { return contiguous_; }

   /// One line for the run log (print on rank 0).  Collective.
   std::string Summary() const
   {
      int nmax = 0, ncur = 0;
      MPI_Comm_size(node_comm_, &ncur);
      MPI_Allreduce(&ncur, &nmax, 1, MPI_INT, MPI_MAX, comm_);
      std::ostringstream o;
      o << "[QD-GREEN] slip gather node: " << nnodes_ << " node(s), up to "
        << nmax << " ranks per node, 2 windows of " << std::fixed
        << std::setprecision(3) << 2.0 * nf_ * sizeof(real_t) / 1.0e6
        << " MB per node, exchange "
        << (nnodes_ == 1 ? "none (one node)" :
            contiguous_ ? "contiguous" : "packed") << "\n";
      return o.str();
   }

private:
   void SetupNodes()
   {
      const char *fake = std::getenv("SEAS_QD_HMAT_NODE_SIZE");
      if (fake && std::atoi(fake) > 0)
      {
         MPI_Comm_split(comm_, rank_ / std::atoi(fake), rank_, &node_comm_);
      }
      else
      {
         MPI_Comm_split_type(comm_, MPI_COMM_TYPE_SHARED, rank_,
                             MPI_INFO_NULL, &node_comm_);
      }
      MPI_Comm_rank(node_comm_, &node_rank_);
      MPI_Comm_split(comm_, node_rank_ == 0 ? 0 : MPI_UNDEFINED, rank_,
                     &lead_comm_);
      // Both orders are by world rank, so node n is rank n of lead_comm_
      // node_of_[r]: the node of rank r, numbered by its rank 0's order in
      // lead_comm_ (world order)
      int leader = rank_;
      MPI_Bcast(&leader, 1, MPI_INT, 0, node_comm_);
      std::vector<int> leaders(nranks_);
      MPI_Allgather(&leader, 1, MPI_INT, leaders.data(), 1, MPI_INT, comm_);
      std::vector<int> ids(leaders);
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      nnodes_ = static_cast<int>(ids.size());
      node_of_.resize(nranks_);
      node_ranks_.assign(nnodes_, {});
      for (int r = 0; r < nranks_; r++)
      {
         node_of_[r] = static_cast<int>(std::lower_bound(ids.begin(),
                                                         ids.end(),
                                                         leaders[r]) -
                                        ids.begin());
         node_ranks_[node_of_[r]].push_back(r);
      }
      if (node_rank_ == 0)
      {
         int lead_rank = -1;
         MPI_Comm_rank(lead_comm_, &lead_rank);
         MFEM_VERIFY(lead_rank == node_of_[rank_], "QD node slip: node "
                     << node_of_[rank_] << " has rank " << lead_rank
                     << " among the nodes' rank 0");
      }
   }

   /// Contiguous when each node's ranks own one unbroken range in rank
   /// order (ranks that own nothing aside): then node n's slip is global
   /// [2·node_lo, 2·(node_lo + node_n)) and the leaders swap it in place.
   /// Otherwise they pack and unpack.
   void PlanExchange()
   {
      contiguous_ = true;
      node_lo_.assign(nnodes_, 0);
      node_n_.assign(nnodes_, 0);
      for (int n = 0; n < nnodes_; n++)
      {
         bool first = true;
         for (int r : node_ranks_[n])
         {
            if (own_n_[r] == 0) { continue; }
            if (first) { node_lo_[n] = own_lo_[r]; first = false; }
            if (own_lo_[r] != node_lo_[n] + node_n_[n]) { contiguous_ = false; }
            node_n_[n] += own_n_[r];
         }
      }
      cnt_.assign(nnodes_, 0);
      dsp_.assign(nnodes_, 0);
      int total = 0;
      for (int n = 0; n < nnodes_; n++)
      {
         cnt_[n] = static_cast<int>(2 * node_n_[n]);
         dsp_[n] = contiguous_ ? static_cast<int>(2 * node_lo_[n]) : total;
         total += cnt_[n];
      }
      MFEM_VERIFY(total == 2 * nf_, "QD node slip: the ranks own "
                  << total / 2 << " fault DOFs, expected " << nf_);
      if (!contiguous_ && node_rank_ == 0)
      {
         send_.resize(cnt_[node_of_[rank_]]);
         recv_.resize(2 * nf_);
      }
   }

   /// Node rank 0: swap the node's part of @a w with the other nodes'.
   void Exchange(real_t *w)
   {
      const MPI_Datatype t = MPITypeMap<real_t>::mpi_type;
      if (contiguous_)
      {
         MPI_Allgatherv(MPI_IN_PLACE, 0, t, w, cnt_.data(), dsp_.data(), t,
                        lead_comm_);
         return;
      }
      size_t k = 0;
      for (int r : node_ranks_[node_of_[rank_]])
      {
         if (own_n_[r] == 0) { continue; }
         const real_t *s = w + 2 * own_lo_[r];
         std::copy(s, s + 2 * own_n_[r], send_.data() + k);
         k += 2 * own_n_[r];
      }
      MPI_Allgatherv(send_.data(), static_cast<int>(send_.size()), t,
                     recv_.data(), cnt_.data(), dsp_.data(), t, lead_comm_);
      for (int n = 0; n < nnodes_; n++)
      {
         if (n == node_of_[rank_]) { continue; }
         k = static_cast<size_t>(dsp_[n]);
         for (int r : node_ranks_[n])
         {
            if (own_n_[r] == 0) { continue; }
            std::copy(recv_.data() + k, recv_.data() + k + 2 * own_n_[r],
                      w + 2 * own_lo_[r]);
            k += 2 * own_n_[r];
         }
      }
   }

   MPI_Comm comm_;
   int nranks_ = 1, rank_ = 0;
   std::vector<long long> own_lo_, own_n_;
   long long nf_;
   MPI_Comm node_comm_ = MPI_COMM_NULL;   ///< the ranks of this node
   MPI_Comm lead_comm_ = MPI_COMM_NULL;   ///< node rank 0 of every node
   int node_rank_ = 0, nnodes_ = 1;
   std::vector<int> node_of_;
   std::vector<std::vector<int>> node_ranks_;
   std::vector<long long> node_lo_, node_n_;
   bool contiguous_ = true;
   std::vector<int> cnt_, dsp_;           ///< leaders' Allgatherv, per node
   std::vector<real_t> send_, recv_;      ///< packed exchange only
   MPI_Win win_[2] = {MPI_WIN_NULL, MPI_WIN_NULL};
   real_t *base_[2] = {nullptr, nullptr};
   long long calls_ = 0;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_SEAS_QD_NODE_SLIP_HPP
