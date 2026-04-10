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

// Unit tests for MPIContext parallel functionality
// Run with: mpirun -np 2 ./seas_test_mpi_context
//       or: mpirun -np 4 ./seas_test_mpi_context

#include "mfem.hpp"
#include "../../common/mpi_context.hpp"
#include "../../common/mpi_check.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

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
      std::cout << "Running MPIContext tests with " << size
                << " processes" << std::endl;
   }

   // Test 1: Initialization
   TEST_CHECK("rank >= 0", rank >= 0);
   TEST_CHECK("rank < size", rank < size);
   TEST_CHECK("size > 0", size > 0);

   // Test 2: GlobalMax — each rank sends its rank, result should be size-1
   {
      real_t result = ctx.GlobalMax(static_cast<real_t>(rank));
      TEST_CHECK("GlobalMax", std::abs(result - (size - 1)) < 1e-14);
   }

   // Test 3: GlobalMin — each rank sends its rank, result should be 0
   {
      real_t result = ctx.GlobalMin(static_cast<real_t>(rank));
      TEST_CHECK("GlobalMin", std::abs(result) < 1e-14);
   }

   // Test 4: GlobalSum — each rank sends 1.0, result should be size
   {
      real_t result = ctx.GlobalSum(1.0);
      TEST_CHECK("GlobalSum", std::abs(result - size) < 1e-14);
   }

   // Test 5: GlobalSumInt — each rank sends 1, result should be size
   {
      int result = ctx.GlobalSumInt(1);
      TEST_CHECK("GlobalSumInt", result == size);
   }

   // Test 6: Barrier — all ranks reach barrier without deadlock
   {
      ctx.Barrier();
      TEST_CHECK("Barrier", true);
   }

   // Test 7: Bcast scalar — root sets value, all ranks receive it
   {
      real_t value = ctx.IsRoot() ? 42.0 : 0.0;
      ctx.Bcast(value);
      TEST_CHECK("Bcast scalar", std::abs(value - 42.0) < 1e-14);
   }

   // Test 8: Bcast vector — root sets vector, all ranks receive it
   {
      Vector vec;
      if (ctx.IsRoot())
      {
         vec.SetSize(3);
         vec(0) = 1.0;
         vec(1) = 2.0;
         vec(2) = 3.0;
      }
      ctx.Bcast(vec);
      TEST_CHECK("Bcast vector size", vec.Size() == 3);
      TEST_CHECK("Bcast vector values",
                 std::abs(vec(0) - 1.0) < 1e-14 &&
                 std::abs(vec(1) - 2.0) < 1e-14 &&
                 std::abs(vec(2) - 3.0) < 1e-14);
   }

   // Test 9: Serial comparison — single-rank results should match serial stubs
   if (size == 1)
   {
      real_t val = 7.5;
      TEST_CHECK("serial GlobalMax identity", std::abs(ctx.GlobalMax(val) - val) < 1e-14);
      TEST_CHECK("serial GlobalMin identity", std::abs(ctx.GlobalMin(val) - val) < 1e-14);
      TEST_CHECK("serial GlobalSum identity", std::abs(ctx.GlobalSum(val) - val) < 1e-14);
      TEST_CHECK("serial GlobalSumInt identity", ctx.GlobalSumInt(3) == 3);
   }

   // Test 10: Bcast int — root sets value, all ranks receive it
   {
      int ival = ctx.IsRoot() ? 99 : 0;
      ctx.Bcast(ival);
      TEST_CHECK("Bcast int", ival == 99);
   }

   // Test 11: GetComm — stored communicator works
   {
#ifdef SEAS_USE_MPI
      int comm_rank = -1;
      MPI_Comm_rank(ctx.GetComm(), &comm_rank);
      TEST_CHECK("GetComm rank matches", comm_rank == rank);
#else
      TEST_CHECK("GetComm serial stub", true);
#endif
   }

   // Test 12: GatherToRoot — each rank sends its rank, root gets [0,1,...,size-1]
   {
      Vector local(1);
      local(0) = static_cast<real_t>(rank);
      Vector global;
      ctx.GatherToRoot(local, global);
      if (ctx.IsRoot())
      {
         bool ok = (global.Size() == size);
         for (int i = 0; i < size && ok; i++)
         {
            ok = (std::abs(global(i) - i) < 1e-14);
         }
         TEST_CHECK("GatherToRoot values", ok);
      }
      else
      {
         TEST_CHECK("GatherToRoot non-root empty", global.Size() == 0);
      }
   }

   // Test 13: GatherToRoot with varying sizes
   {
      Vector local(rank + 1);  // rank 0 sends 1 element, rank 1 sends 2, etc.
      for (int i = 0; i < local.Size(); i++)
      {
         local(i) = rank * 100.0 + i;
      }
      Vector global;
      ctx.GatherToRoot(local, global);
      if (ctx.IsRoot())
      {
         int expected_total = size * (size + 1) / 2;
         TEST_CHECK("GatherToRoot varying sizes", global.Size() == expected_total);
      }
      else
      {
         TEST_CHECK("GatherToRoot varying non-root", global.Size() == 0);
      }
   }

   // Test 14: AllgatherVec — each rank sends 1 element, all get [0,1,...,size-1]
   {
      Vector local(1);
      local(0) = static_cast<real_t>(rank);
      Vector global;
      ctx.AllgatherVec(local, global);
      bool ok = (global.Size() == size);
      for (int i = 0; i < size && ok; i++)
      {
         ok = (std::abs(global(i) - i) < 1e-14);
      }
      TEST_CHECK("AllgatherVec values", ok);
   }

   // Test 15: AllgatherVec with varying sizes
   {
      Vector local(rank + 1);
      for (int i = 0; i < local.Size(); i++)
      {
         local(i) = rank * 100.0 + i;
      }
      Vector global;
      ctx.AllgatherVec(local, global);
      int expected_total = size * (size + 1) / 2;
      TEST_CHECK("AllgatherVec varying sizes total", global.Size() == expected_total);
      // Verify first element from each rank
      int offset = 0;
      bool ok = true;
      for (int r = 0; r < size; r++)
      {
         if (std::abs(global(offset) - r * 100.0) > 1e-14) { ok = false; }
         offset += r + 1;
      }
      TEST_CHECK("AllgatherVec varying values", ok);
   }

   // Test 16: Stored communicator constructor (dup'd comm)
#ifdef SEAS_USE_MPI
   {
      MPIContext ctx2(MPI_COMM_WORLD);
      TEST_CHECK("dup'd comm rank matches", ctx2.Rank() == rank);
      TEST_CHECK("dup'd comm size matches", ctx2.Size() == size);
      real_t result = ctx2.GlobalSum(1.0);
      TEST_CHECK("dup'd comm GlobalSum", std::abs(result - size) < 1e-14);
   }
#endif

   // Test 17: MFEM_SEAS_MPI_CHECK compiles and works with valid call
   {
#ifdef SEAS_USE_MPI
      MFEM_SEAS_MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
#endif
      TEST_CHECK("MPI_CHECK compiles with valid call", true);
   }

   // Test 18: GatherToRoot with empty input
   {
      Vector empty_local(0);
      Vector global;
      ctx.GatherToRoot(empty_local, global);
      if (ctx.IsRoot())
      {
         TEST_CHECK("GatherToRoot empty", global.Size() == 0);
      }
      else
      {
         TEST_CHECK("GatherToRoot empty non-root", global.Size() == 0);
      }
   }

   // Print results from root
   if (ctx.IsRoot())
   {
      std::cout << "MPIContext tests: " << num_passed << " / " << num_tests
                << " passed" << std::endl;
   }

   return (num_passed == num_tests) ? 0 : 1;
}
