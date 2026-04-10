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

// Unit tests for FaultScatter parallel communication
// Run with: mpirun -np 2 ./seas_test_fault_scatter
//       or: mpirun -np 4 ./seas_test_fault_scatter

#include "mfem.hpp"
#include "../../common/fault_scatter.hpp"
#include "../../common/mpi_context.hpp"
#include <iostream>
#include <cmath>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0;
static int num_passed = 0;

#define TEST_CHECK(name, condition) \
   do { \
      num_tests++; \
      if (condition) { num_passed++; } \
      else if (ctx.IsRoot()) { \
         std::cerr << "FAIL: " << name << std::endl; \
      } \
   } while(0)

int main(int argc, char *argv[])
{
   MPIContext ctx(&argc, &argv);

   int rank = ctx.Rank();
   int size = ctx.Size();

   if (ctx.IsRoot())
   {
      std::cout << "Running FaultScatter tests with " << size
                << " processes" << std::endl;
   }

#ifdef SEAS_USE_MPI
   // ========================================================================
   // Test 1: Empty blocks — no communication needed
   // ========================================================================
   {
      std::vector<SharedFaultCommBlock> empty_blocks;
      FaultScatter scatter(empty_blocks, MPI_COMM_WORLD);
      Vector owned(0);
      scatter.BeginScatter(owned, 1, 1);
      TEST_CHECK("empty blocks TestScatter", scatter.TestScatter());
      scatter.WaitScatter();
      TEST_CHECK("empty blocks NumBlocks", scatter.NumBlocks() == 0);
   }

   // ========================================================================
   // Test 2: Two-rank exchange — rank 0 sends to rank 1 and vice versa
   // ========================================================================
   if (size >= 2)
   {
      // Each rank owns 2 faces, 1 basis function per face, 1 component.
      // Rank 0 sends face 0 to rank 1; rank 1 sends face 0 to rank 0.
      const int nbf_per_face = 1;
      const int comps = 1;
      const int num_owned_faces = 2;

      int partner = (rank == 0) ? 1 : (rank == 1) ? 0 : -1;

      std::vector<SharedFaultCommBlock> blocks;
      if (rank == 0 || rank == 1)
      {
         SharedFaultCommBlock block;
         block.neighbor_rank = partner;
         block.send_owned_faces = {0};      // send face 0
         block.recv_local_faces = {0};      // receive into local position 0
         blocks.push_back(block);
      }

      // Owned data: face 0 = rank*10, face 1 = rank*10+1
      Vector owned(num_owned_faces * nbf_per_face * comps);
      owned(0) = rank * 10.0;
      owned(1) = rank * 10.0 + 1.0;

      if (rank == 0 || rank == 1)
      {
         FaultScatter scatter(blocks, MPI_COMM_WORLD);
         scatter.BeginScatter(owned, comps, nbf_per_face);
         scatter.WaitScatter();

         // Rank 0 should receive rank 1's face 0 data (= 10.0)
         // Rank 1 should receive rank 0's face 0 data (= 0.0)
         const Vector &recv = scatter.GetRecvBuffer(0);
         TEST_CHECK("two-rank recv size", recv.Size() == 1);

         double expected = partner * 10.0;
         TEST_CHECK("two-rank recv value",
                    std::abs(recv(0) - expected) < 1e-14);
      }
      else
      {
         // Ranks 2+ don't participate but test still runs
         TEST_CHECK("two-rank non-participant", true);
         TEST_CHECK("two-rank non-participant 2", true);
      }
   }

   // ========================================================================
   // Test 3: Multi-component exchange (3 components per DOF)
   // ========================================================================
   if (size >= 2)
   {
      const int nbf_per_face = 1;
      const int comps = 3;
      const int num_owned_faces = 1;

      int partner = (rank == 0) ? 1 : (rank == 1) ? 0 : -1;

      std::vector<SharedFaultCommBlock> blocks;
      if (rank == 0 || rank == 1)
      {
         SharedFaultCommBlock block;
         block.neighbor_rank = partner;
         block.send_owned_faces = {0};
         block.recv_local_faces = {0};
         blocks.push_back(block);
      }

      Vector owned(num_owned_faces * nbf_per_face * comps);
      for (int c = 0; c < comps; c++)
      {
         owned(c) = rank * 100.0 + c;
      }

      if (rank == 0 || rank == 1)
      {
         FaultScatter scatter(blocks, MPI_COMM_WORLD);
         scatter.BeginScatter(owned, comps, nbf_per_face);
         scatter.WaitScatter();

         const Vector &recv = scatter.GetRecvBuffer(0);
         TEST_CHECK("multi-comp recv size", recv.Size() == comps);

         bool ok = true;
         for (int c = 0; c < comps; c++)
         {
            double expected = partner * 100.0 + c;
            if (std::abs(recv(c) - expected) > 1e-14) { ok = false; }
         }
         TEST_CHECK("multi-comp recv values", ok);
      }
      else
      {
         TEST_CHECK("multi-comp non-participant", true);
         TEST_CHECK("multi-comp non-participant 2", true);
      }
   }

   // ========================================================================
   // Test 4: Asymmetric block sizes — rank 0 sends 2 faces, rank 1 sends 1
   // ========================================================================
   if (size >= 2)
   {
      const int nbf_per_face = 1;
      const int comps = 1;

      std::vector<SharedFaultCommBlock> blocks;
      if (rank == 0)
      {
         // Rank 0: sends faces 0,1 to rank 1; receives 1 face
         SharedFaultCommBlock block;
         block.neighbor_rank = 1;
         block.send_owned_faces = {0, 1};
         block.recv_local_faces = {0};
         blocks.push_back(block);
      }
      else if (rank == 1)
      {
         // Rank 1: sends face 0 to rank 0; receives 2 faces
         SharedFaultCommBlock block;
         block.neighbor_rank = 0;
         block.send_owned_faces = {0};
         block.recv_local_faces = {0, 1};
         blocks.push_back(block);
      }

      // Rank 0 owns 3 faces, rank 1 owns 2 faces
      int num_owned = (rank == 0) ? 3 : 2;
      Vector owned(num_owned * nbf_per_face * comps);
      for (int i = 0; i < owned.Size(); i++)
      {
         owned(i) = rank * 100.0 + i;
      }

      if (rank == 0 || rank == 1)
      {
         FaultScatter scatter(blocks, MPI_COMM_WORLD);
         scatter.BeginScatter(owned, comps, nbf_per_face);
         scatter.WaitScatter();

         const Vector &recv = scatter.GetRecvBuffer(0);
         if (rank == 0)
         {
            // Receives 1 face from rank 1: value = 100.0
            TEST_CHECK("asymmetric r0 recv size", recv.Size() == 1);
            TEST_CHECK("asymmetric r0 recv value",
                       std::abs(recv(0) - 100.0) < 1e-14);
         }
         else
         {
            // Receives 2 faces from rank 0: values = 0.0, 1.0
            TEST_CHECK("asymmetric r1 recv size", recv.Size() == 2);
            TEST_CHECK("asymmetric r1 recv values",
                       std::abs(recv(0) - 0.0) < 1e-14 &&
                       std::abs(recv(1) - 1.0) < 1e-14);
         }
      }
      else
      {
         TEST_CHECK("asymmetric non-participant", true);
         TEST_CHECK("asymmetric non-participant 2", true);
      }
   }

   // ========================================================================
   // Test 5: Multi-rank ring exchange (all ranks participate)
   // ========================================================================
   if (size >= 2)
   {
      const int nbf_per_face = 1;
      const int comps = 1;
      const int num_owned = 1;

      // Each rank sends to (rank+1)%size and receives from (rank-1+size)%size
      int send_to = (rank + 1) % size;
      int recv_from = (rank - 1 + size) % size;

      std::vector<SharedFaultCommBlock> blocks;
      SharedFaultCommBlock block;
      block.neighbor_rank = send_to;
      block.send_owned_faces = {0};
      block.recv_local_faces = {};  // No recv from send_to
      blocks.push_back(block);

      SharedFaultCommBlock block2;
      block2.neighbor_rank = recv_from;
      block2.send_owned_faces = {};  // No send to recv_from
      block2.recv_local_faces = {0};
      blocks.push_back(block2);

      Vector owned(num_owned * nbf_per_face * comps);
      owned(0) = rank * 1000.0;

      FaultScatter scatter(blocks, MPI_COMM_WORLD);
      scatter.BeginScatter(owned, comps, nbf_per_face);
      scatter.WaitScatter();

      // block 1 (recv_from) should have data from rank (rank-1+size)%size
      const Vector &recv = scatter.GetRecvBuffer(1);
      TEST_CHECK("ring recv size", recv.Size() == 1);
      double expected = recv_from * 1000.0;
      TEST_CHECK("ring recv value", std::abs(recv(0) - expected) < 1e-14);
   }

   // ========================================================================
   // Test 6: Multiple basis functions per face (nbf_per_face = 3)
   // ========================================================================
   if (size >= 2)
   {
      const int nbf_per_face = 3;
      const int comps = 2;

      int partner = (rank == 0) ? 1 : (rank == 1) ? 0 : -1;

      std::vector<SharedFaultCommBlock> blocks;
      if (rank == 0 || rank == 1)
      {
         SharedFaultCommBlock block;
         block.neighbor_rank = partner;
         block.send_owned_faces = {0};
         block.recv_local_faces = {0};
         blocks.push_back(block);
      }

      // 1 face * 3 dofs/face * 2 comps = 6 values
      Vector owned(1 * nbf_per_face * comps);
      for (int kk = 0; kk < nbf_per_face; kk++)
      {
         for (int c = 0; c < comps; c++)
         {
            owned(kk * comps + c) = rank * 1000.0 + kk * 10.0 + c;
         }
      }

      if (rank == 0 || rank == 1)
      {
         FaultScatter scatter(blocks, MPI_COMM_WORLD);
         scatter.BeginScatter(owned, comps, nbf_per_face);
         scatter.WaitScatter();

         const Vector &recv = scatter.GetRecvBuffer(0);
         TEST_CHECK("multi-bf recv size", recv.Size() == nbf_per_face * comps);

         bool ok = true;
         for (int kk = 0; kk < nbf_per_face; kk++)
         {
            for (int c = 0; c < comps; c++)
            {
               double expected = partner * 1000.0 + kk * 10.0 + c;
               if (std::abs(recv(kk * comps + c) - expected) > 1e-14)
               {
                  ok = false;
               }
            }
         }
         TEST_CHECK("multi-bf recv values", ok);
      }
      else
      {
         TEST_CHECK("multi-bf non-participant", true);
         TEST_CHECK("multi-bf non-participant 2", true);
      }
   }

#else
   // Serial mode: just verify construction works
   {
      std::vector<SharedFaultCommBlock> empty_blocks;
      FaultScatter scatter(empty_blocks);
      Vector dummy;
      scatter.BeginScatter(dummy, 1, 1);
      TEST_CHECK("serial TestScatter", scatter.TestScatter());
      scatter.WaitScatter();
      TEST_CHECK("serial NumBlocks", scatter.NumBlocks() == 0);
   }
#endif

   // Print results from root
   if (ctx.IsRoot())
   {
      std::cout << "FaultScatter tests: " << num_passed << " / " << num_tests
                << " passed" << std::endl;
   }

   return (num_passed == num_tests) ? 0 : 1;
}
