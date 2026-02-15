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

// Unit tests for parallel utility functions
// Run with: mpirun -np 2 ./seas_test_parallel_utils
//       or: mpirun -np 4 ./seas_test_parallel_utils

#include "mfem.hpp"
#include "../../common/mpi_context.hpp"
#include "../../common/parallel_utils.hpp"
#include <iostream>
#include <cmath>
#include <numeric>
#include <vector>

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
      std::cout << "Running parallel_utils tests with " << size
                << " processes" << std::endl;
   }

   // Test 1: GatherVectorToRoot
   // Each rank creates a local vector [rank*10, rank*10+1, ..., rank*10+(rank)]
   // so rank 0 has 1 element, rank 1 has 2 elements, etc.
   {
      int local_size = rank + 1;
      std::vector<real_t> local(local_size);
      for (int i = 0; i < local_size; i++)
      {
         local[i] = rank * 10.0 + i;
      }

      std::vector<real_t> global;
      GatherVectorToRoot(local, global, ctx.GetComm());

      if (ctx.IsRoot())
      {
         // Total size = 1 + 2 + ... + size = size*(size+1)/2
         int expected_size = size * (size + 1) / 2;
         TEST_CHECK("GatherVectorToRoot size",
                    static_cast<int>(global.size()) == expected_size);

         // Check first element from each rank
         bool values_ok = true;
         int offset = 0;
         for (int r = 0; r < size; r++)
         {
            if (std::abs(global[offset] - r * 10.0) > 1e-14)
            {
               values_ok = false;
            }
            offset += r + 1;
         }
         TEST_CHECK("GatherVectorToRoot values", values_ok);
      }
      else
      {
         // Non-root: just count these tests as passed (root validates)
         num_tests += 2;
         num_passed += 2;
      }
   }

   // Test 2: ScatterVectorFromRoot
   // Root creates global vector [0, 1, 2, ..., 2*size-1]
   // Each rank gets 2 elements
   {
      std::vector<real_t> global;
      std::vector<int> counts(size, 2);

      if (ctx.IsRoot())
      {
         global.resize(2 * size);
         for (int i = 0; i < 2 * size; i++)
         {
            global[i] = static_cast<real_t>(i);
         }
      }

      std::vector<real_t> local;
      ScatterVectorFromRoot(global, counts, local, ctx.GetComm());

      TEST_CHECK("ScatterVectorFromRoot size",
                 static_cast<int>(local.size()) == 2);
      TEST_CHECK("ScatterVectorFromRoot values",
                 std::abs(local[0] - 2.0 * rank) < 1e-14 &&
                 std::abs(local[1] - (2.0 * rank + 1)) < 1e-14);
   }

   // Test 3: BroadcastFromRoot
   {
      int int_val = ctx.IsRoot() ? 99 : 0;
      BroadcastFromRoot(int_val, ctx.GetComm());
      TEST_CHECK("BroadcastFromRoot int", int_val == 99);

      real_t real_val = ctx.IsRoot() ? 3.14 : 0.0;
      BroadcastFromRoot(real_val, ctx.GetComm());
      TEST_CHECK("BroadcastFromRoot real_t", std::abs(real_val - 3.14) < 1e-14);
   }

   // Test 4: AllRanksAgree
   {
      // All ranks have same value => should agree
      bool agree = AllRanksAgree(42, ctx.GetComm());
      TEST_CHECK("AllRanksAgree same value", agree);

      // Each rank has different value => should NOT agree (if size > 1)
      bool disagree = AllRanksAgree(rank, ctx.GetComm());
      if (size > 1)
      {
         TEST_CHECK("AllRanksAgree different values", !disagree);
      }
      else
      {
         TEST_CHECK("AllRanksAgree single rank", disagree);
      }
   }

   // Test 5: Serial-parallel comparison
   // Compute sum of [1..N] both via parallel reduction and directly
   {
      int N = 100;
      // Distribute range across ranks
      int per_rank = N / size;
      int start = rank * per_rank + 1;
      int end = (rank == size - 1) ? N : (rank + 1) * per_rank;

      real_t local_sum = 0.0;
      for (int i = start; i <= end; i++)
      {
         local_sum += static_cast<real_t>(i);
      }

      real_t parallel_sum = ctx.GlobalSum(local_sum);
      real_t expected = N * (N + 1) / 2.0;
      TEST_CHECK("serial-parallel sum comparison",
                 std::abs(parallel_sum - expected) < 1e-10);
   }

   // Print results from root
   if (ctx.IsRoot())
   {
      std::cout << "parallel_utils tests: " << num_passed << " / "
                << num_tests << " passed" << std::endl;
   }

   return (num_passed == num_tests) ? 0 : 1;
}
