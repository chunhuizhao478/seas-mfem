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

// test_qd_node_slip.cpp — QD speed item B (qd-b-design-rev1): the node slip
// copy of qd_node_slip.hpp holds exactly what MPI_Allgatherv delivers, and
// stays valid until the caller's next Gather.  No mesh.
//
//   mpirun -np N ./seas_test_qd_node_slip [--nf NF] [--stages S]
//
//   For nodes from the machine and fake nodes of 1, 2 and 4 ranks
//   (SEAS_QD_HMAT_NODE_SIZE), and three ownerships of NF fault DOFs (rank
//   order with uneven sizes and some empty ranks, reversed rank order,
//   shuffled):
//   1. the node count and the contiguity flag are the expected ones,
//   2. for S stages of hash-valued slips, the window equals the slip in
//      global order byte for byte (memcmp), read after a random sleep of
//      up to 300 us, so ranks that leave a stage early have already
//      written the next one,
//   3. at the first and the last stage it equals an MPI_Allgatherv
//      byte for byte.
//
// Exit code 0 when every check passes.

#include "mfem.hpp"
#include "../../domain/qd_node_slip.hpp"

#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_fail = 0;

void Check(bool ok, const std::string &what, int rank)
{
   if (!ok) { g_fail = 1; }
   if (rank == 0)
   {
      std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
   }
}

/// Slip value i of stage @a stage: any double, from a hash.
real_t Value(long long stage, long long i)
{
   uint64_t z = static_cast<uint64_t>(stage) * 0x9E3779B97F4A7C15ull +
                static_cast<uint64_t>(i) + 0x632BE59BD9B4E019ull;
   z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
   z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
   z ^= z >> 31;
   return (static_cast<real_t>(z >> 11) / 9007199254740992.0 - 0.5) *
          std::ldexp(1.0, static_cast<int>(z & 31) - 16);
}

/// Ownership of @a nf fault DOFs: uneven sizes (every third rank empty),
/// ranges laid out in rank order (kind 0), reversed (1) or shuffled (2).
void Ownership(int kind, long long nf, int P, std::vector<long long> &lo,
               std::vector<long long> &n)
{
   std::vector<long long> w(P);
   long long wsum = 0;
   for (int r = 0; r < P; r++)
   {
      w[r] = (r % 3 == 2) ? 0 : 1 + (7 * r) % 5;
      wsum += w[r];
   }
   n.assign(P, 0);
   long long acc = 0;
   for (int r = 0; r < P; r++)
   {
      n[r] = nf * (acc + w[r]) / wsum - nf * acc / wsum;
      acc += w[r];
   }
   std::vector<int> order(P);
   for (int r = 0; r < P; r++) { order[r] = r; }
   if (kind == 1) { std::reverse(order.begin(), order.end()); }
   if (kind == 2)
   {
      std::mt19937 gen(12345);
      for (int i = P - 1; i > 0; i--)
      {
         const int j = static_cast<int>(gen() % static_cast<unsigned>(i + 1));
         std::swap(order[i], order[j]);
      }
   }
   lo.assign(P, 0);
   long long at = 0;
   for (int r : order) { lo[r] = at; at += n[r]; }
}

}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank(), P = Mpi::WorldSize();
   long long nf = 1000;
   int nstages = 200;
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--nf" && i + 1 < argc) { nf = std::atoll(argv[++i]); }
      else if (a == "--stages" && i + 1 < argc)
      {
         nstages = std::atoi(argv[++i]);
      }
   }
   if (rank == 0)
   {
      std::cout << "test_qd_node_slip: " << P << " ranks, N_f " << nf
                << ", " << nstages << " stages" << std::endl;
   }
   // Machine nodes, for the expected node of each rank without the hook
   int machine_leader = rank;
   {
      MPI_Comm shm;
      MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
                          &shm);
      MPI_Allreduce(&rank, &machine_leader, 1, MPI_INT, MPI_MIN, shm);
      MPI_Comm_free(&shm);
   }
   std::vector<int> machine_leaders(P);
   MPI_Allgather(&machine_leader, 1, MPI_INT, machine_leaders.data(), 1,
                 MPI_INT, comm);

   const char *sizes[] = {"", "1", "2", "4"};
   const char *kinds[] = {"rank order", "reversed", "shuffled"};
   for (const char *ns : sizes)
   {
      if (ns[0]) { setenv("SEAS_QD_HMAT_NODE_SIZE", ns, 1); }
      else { unsetenv("SEAS_QD_HMAT_NODE_SIZE"); }
      // Expected node of each rank: the smallest rank on it
      std::vector<int> leader(P);
      for (int r = 0; r < P; r++)
      {
         leader[r] = ns[0] ? (r / std::atoi(ns)) * std::atoi(ns)
                     : machine_leaders[r];
      }
      std::vector<int> ids(leader);
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      const int nnodes = static_cast<int>(ids.size());
      const std::string nodes = ns[0] ? std::string("nodes of ") + ns +
                                std::string(" rank(s)") :
                                std::string("machine nodes");
      for (int kind = 0; kind < 3; kind++)
      {
         std::vector<long long> lo, n;
         Ownership(kind, nf, P, lo, n);
         // Expected contiguity: on every node, the non-empty ranks in rank
         // order own increasing ranges that leave no gap
         bool contiguous = true;
         for (int id : ids)
         {
            long long first = -1, end = -1, sum = 0;
            for (int r = 0; r < P; r++)
            {
               if (leader[r] != id || n[r] == 0) { continue; }
               if (first < 0) { first = lo[r]; }
               else if (lo[r] < end) { contiguous = false; }
               end = std::max(end, lo[r] + n[r]);
               sum += n[r];
            }
            if (first >= 0 && end - first != sum) { contiguous = false; }
         }
         const std::string what = nodes + ", " + kinds[kind];
         QDNodeSlip ns_op(comm, lo, n, nf);
         const std::string sum = ns_op.Summary();
         if (rank == 0) { std::cout << "  " << what << ": " << sum; }
         Check(ns_op.NumNodes() == nnodes, what + ": " +
               std::to_string(nnodes) + " node(s)", rank);
         Check(ns_op.Contiguous() == contiguous, what + ": exchange " +
               (contiguous ? "contiguous" : "packed"), rank);

         std::vector<real_t> mine(2 * n[rank]), want(2 * nf), ref(2 * nf);
         std::vector<int> cnt(P), dsp(P);
         for (int r = 0; r < P; r++)
         {
            cnt[r] = static_cast<int>(2 * n[r]);
            dsp[r] = static_cast<int>(2 * lo[r]);
         }
         std::mt19937 gen(1000u * static_cast<unsigned>(rank) + 17u * kind +
                          static_cast<unsigned>(std::atoi(ns)));
         int bad_stage = -1, bad_ref = 0;
         for (int s = 0; s < nstages; s++)
         {
            for (long long i = 0; i < 2 * n[rank]; i++)
            {
               mine[i] = Value(s, 2 * lo[rank] + i);
            }
            for (long long i = 0; i < 2 * nf; i++) { want[i] = Value(s, i); }
            const real_t *w = ns_op.Gather(mine.data());
            if (gen() % 2) { usleep(gen() % 300); }
            if (bad_stage < 0 &&
                std::memcmp(w, want.data(), 2 * nf * sizeof(real_t)) != 0)
            {
               bad_stage = s;
            }
            if (s == 0 || s == nstages - 1)
            {
               MPI_Allgatherv(mine.data(), cnt[rank],
                              MPITypeMap<real_t>::mpi_type, ref.data(),
                              cnt.data(), dsp.data(),
                              MPITypeMap<real_t>::mpi_type, comm);
               if (std::memcmp(w, ref.data(), 2 * nf * sizeof(real_t)) != 0)
               {
                  bad_ref = 1;
               }
            }
         }
         int bad_all = 0, ref_all = 0;
         int bad = bad_stage >= 0 ? 1 : 0;
         MPI_Allreduce(&bad, &bad_all, 1, MPI_INT, MPI_MAX, comm);
         MPI_Allreduce(&bad_ref, &ref_all, 1, MPI_INT, MPI_MAX, comm);
         Check(!bad_all, what + ": " + std::to_string(nstages) +
               " stages byte for byte, read after a random sleep", rank);
         Check(!ref_all, what + ": equals MPI_Allgatherv at the first and "
               "last stage", rank);
      }
   }
   unsetenv("SEAS_QD_HMAT_NODE_SIZE");
   int fail_all = 0;
   MPI_Allreduce(&g_fail, &fail_all, 1, MPI_INT, MPI_MAX, comm);
   if (rank == 0)
   {
      std::cout << (fail_all ? "FAILED" : "ALL PASSED") << std::endl;
   }
   return fail_all ? 1 : 0;
}
