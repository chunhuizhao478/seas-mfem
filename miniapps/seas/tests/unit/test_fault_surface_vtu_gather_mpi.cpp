// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 1 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
//
// MPI gather correctness for `vtu::GatherFaultPackToRoot`.  Run with
// `mpirun -np 4 seas_test_fault_surface_vtu_gather_mpi`.
//
// Each rank constructs a small `LocalFaultPack` carrying:
//   - 2 triangles (so the pack has > 1 cell per rank);
//   - 2 scalar fields whose values encode the rank id and the local
//     triangle id (so a missing or duplicated entry is immediately
//     visible after gather).
//
// On rank 0, after the gather we assert:
//   1. Total cell count == sum over ranks (nranks * 2).
//   2. Every (rank, local_idx) pair appears exactly once in the
//      gathered field arrays.
//   3. Vertex ordering matches the triangle ordering: triangle k owns
//      vertices 3k, 3k+1, 3k+2.
//   4. Round-tripping through `WriteFaultPackVTU` -> `ParseFaultVTUCellData`
//      produces the same field arrays bit-exactly (modulo
//      base64 encode/decode).

#include "mfem.hpp"
#include "../../io/fault_vtu_binary.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <sys/stat.h>
#include <unistd.h>     // ::rmdir
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   int rank = 0, nranks = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   if (rank == 0)
   {
      std::cout << "=== test_fault_surface_vtu_gather_mpi ("
                << nranks << " ranks) ===\n";
   }

   // 1. Per-rank LocalFaultPack: 2 triangles, 2 fields.
   const int n_local_cells = 2;
   vtu::LocalFaultPack pack;
   pack.vertices.reserve(3 * n_local_cells);
   pack.triangles.reserve(n_local_cells);
   for (int c = 0; c < n_local_cells; ++c)
   {
      // Triangles in a unit square at z = rank, offset in x by c.
      const double z = static_cast<double>(rank);
      const double x0 = static_cast<double>(c);
      pack.vertices.push_back({x0,        0.0, z});
      pack.vertices.push_back({x0 + 1.0,  0.0, z});
      pack.vertices.push_back({x0 + 0.5,  1.0, z});
      pack.triangles.push_back({3*c, 3*c + 1, 3*c + 2});
   }

   // Field 1: rank id (constant across all of this rank's cells).
   // Field 2: encoded id = rank * 1000 + local_idx.
   pack.field_names = {"rank_id", "global_id"};
   pack.field_arrays.resize(2);
   pack.field_arrays[0].resize(n_local_cells, static_cast<double>(rank));
   pack.field_arrays[1].resize(n_local_cells);
   for (int c = 0; c < n_local_cells; ++c)
   {
      pack.field_arrays[1][c] = static_cast<double>(rank * 1000 + c);
   }

   // 2. Gather collectively.
   vtu::GatheredFaultPack gathered =
      vtu::GatherFaultPackToRoot(pack, rank, nranks, MPI_COMM_WORLD);

   if (rank == 0)
   {
      // (a) Total cell count.
      const int expected_total = nranks * n_local_cells;
      TEST_ASSERT(static_cast<int>(gathered.triangles.size()) == expected_total,
                  "gathered cell count == nranks * 2");
      TEST_ASSERT(static_cast<int>(gathered.vertices.size()) == 3 * expected_total,
                  "gathered vertex count == 3 * cell count");

      // (b) Triangle indexing: each triangle k owns 3k, 3k+1, 3k+2.
      bool tri_layout_ok = true;
      for (std::size_t k = 0; k < gathered.triangles.size(); ++k)
      {
         const auto &t = gathered.triangles[k];
         if (t[0] != static_cast<int>(3 * k)
             || t[1] != static_cast<int>(3 * k + 1)
             || t[2] != static_cast<int>(3 * k + 2))
         {
            tri_layout_ok = false;
            break;
         }
      }
      TEST_ASSERT(tri_layout_ok, "triangle indices are 3k, 3k+1, 3k+2");

      // (c) Every (rank, local_idx) appears exactly once.
      TEST_ASSERT(gathered.field_arrays.size() == 2u,
                  "gathered has 2 fields");
      std::set<int> seen;
      for (std::size_t k = 0; k < gathered.field_arrays[1].size(); ++k)
      {
         seen.insert(static_cast<int>(gathered.field_arrays[1][k]));
      }
      TEST_ASSERT(static_cast<int>(seen.size()) == expected_total,
                  "no duplicate global_id values");
      bool full_set = true;
      for (int r = 0; r < nranks; ++r)
      {
         for (int c = 0; c < n_local_cells; ++c)
         {
            if (seen.count(r * 1000 + c) == 0) { full_set = false; }
         }
      }
      TEST_ASSERT(full_set, "every (rank, local_idx) appears exactly once");

      // (d) Round-trip through WriteFaultPackVTU + ParseFaultVTUCellData.
      const std::string out_dir = "/tmp/test_fault_vtu_gather_mpi";
      ::mkdir(out_dir.c_str(), 0755);
      const std::string vtu_path = out_dir + "/gathered.vtu";
      vtu::WriteFaultPackVTU(vtu_path, gathered, VTKFormat::BINARY, /*lvl=*/0);
      auto roundtrip = vtu::ParseFaultVTUCellData(vtu_path, "global_id");
      TEST_ASSERT(static_cast<int>(roundtrip.size()) == expected_total,
                  "round-trip parsed array size matches");
      bool roundtrip_ok = (roundtrip.size() == gathered.field_arrays[1].size());
      for (std::size_t i = 0; i < roundtrip.size() && roundtrip_ok; ++i)
      {
         if (std::abs(roundtrip[i] - gathered.field_arrays[1][i]) > 1e-12)
         {
            roundtrip_ok = false;
         }
      }
      TEST_ASSERT(roundtrip_ok, "round-trip values bit-exact (within 1e-12)");

      ::remove(vtu_path.c_str());
      ::rmdir(out_dir.c_str());
   }

   // -----------------------------------------------------------------
   // R-001 regression: off-fault rank can pack 0 cells with an
   // empty `_k4`-bearing pack.  Every rank packs the same 17-name
   // field list (the producer-side OR-reduce on `has_k4` makes this
   // uniform); off-fault ranks contribute 0 cells via empty vectors.
   // Pre-fix the gather aborted on the count-mismatch MFEM_VERIFY
   // because the empty rank packed 12 fields and rank 0 packed 17.
   // -----------------------------------------------------------------
   {
      const std::vector<std::string> base_names = {
         "slip_dip","slip_strike","slip_rate_dip","slip_rate_strike",
         "traction_dip","traction_strike","state_variable","normal_stress",
         "param_a","param_Dc","fault_x2","fault_x3"};
      const std::vector<std::string> k4_names = {
         "slip_rate_dip_k4","slip_rate_strike_k4","traction_dip_k4",
         "traction_strike_k4","normal_stress_k4"};

      vtu::LocalFaultPack pack;
      pack.field_names = base_names;
      pack.field_names.insert(pack.field_names.end(),
                              k4_names.begin(), k4_names.end());

      if (rank == 0)
      {
         // Rank 0 has 1 fault triangle with 17 populated fields.
         pack.vertices  = {{0,0,0},{1,0,0},{0,1,0}};
         pack.triangles = {{0,1,2}};
         pack.field_arrays.assign(17, std::vector<double>{42.0});
      }
      else
      {
         // Off-fault ranks: 0 cells, 17 empty fields (matching shape
         // post-R-001 fix).
         pack.field_arrays.assign(17, std::vector<double>{});
      }
      vtu::GatheredFaultPack g =
         vtu::GatherFaultPackToRoot(pack, rank, nranks, MPI_COMM_WORLD);
      if (rank == 0)
      {
         TEST_ASSERT(g.field_names.size() == 17,
                     "R-001: gather completes with 17 fields when only "
                     "rank 0 has fault data");
         TEST_ASSERT(g.triangles.size() == 1u,
                     "R-001: gather receives 1 triangle from rank 0 + 0 "
                     "from off-fault ranks");
      }
   }

   // -----------------------------------------------------------------
   // R-006 regression: validate the post-fix data shape on the
   // contract-violating driver scenario.  Pre-fix (REVIEW.md round 2)
   // a rank with local fault faces but EMPTY `_k4` Vectors entered the
   // per-face loop (gated by global `has_k4`) and OOB-read the empty
   // Vectors.  Post-fix the loop is gated by `has_k4_local` so the
   // accumulators stay at 0.0 on that rank; the lambda's push at
   // paraview_output.hpp:788 then writes 0.0 into c_*_k4; the pack-add
   // (still gated by global `has_k4`) still adds the 5 fields, with
   // those zero values.
   //
   // This block synthesizes the post-loop pack state and verifies the
   // gather concatenates it correctly: rank 0 contributes a populated
   // k4 value, ranks 1..n-1 contribute 0.0.  The end-to-end driver
   // path is exercised by the dynamic-rupture sbatches; this test
   // documents the data invariant the fix enforces.
   // -----------------------------------------------------------------
   {
      const std::vector<std::string> base_names = {
         "slip_dip","slip_strike","slip_rate_dip","slip_rate_strike",
         "traction_dip","traction_strike","state_variable","normal_stress",
         "param_a","param_Dc","fault_x2","fault_x3"};
      const std::vector<std::string> k4_names = {
         "slip_rate_dip_k4","slip_rate_strike_k4","traction_dip_k4",
         "traction_strike_k4","normal_stress_k4"};

      vtu::LocalFaultPack pack;
      pack.field_names = base_names;
      pack.field_names.insert(pack.field_names.end(),
                              k4_names.begin(), k4_names.end());

      // Every rank packs 1 triangle (so the gather concatenates n_ranks
      // cells) and 17 field arrays.  The k4 arrays differ:
      //   rank 0 — populated (post-loop with has_k4_local=true)
      //   rank >0 — zero      (post-loop with has_k4_local=false; the
      //                        accumulators were untouched by the loop
      //                        and the lambda pushed 0.0)
      pack.vertices  = {{static_cast<double>(rank), 0.0, 0.0},
                        {static_cast<double>(rank) + 1.0, 0.0, 0.0},
                        {static_cast<double>(rank) + 0.5, 1.0, 0.0}};
      pack.triangles = {{0, 1, 2}};
      pack.field_arrays.assign(17, std::vector<double>{0.0});
      // Mark rank 0's k4 fields with a recognisable sentinel.
      if (rank == 0)
      {
         for (int k = 12; k < 17; ++k)
         {
            pack.field_arrays[k][0] = 7.0;
         }
      }

      vtu::GatheredFaultPack g =
         vtu::GatherFaultPackToRoot(pack, rank, nranks, MPI_COMM_WORLD);

      if (rank == 0)
      {
         TEST_ASSERT(g.field_names.size() == 17,
                     "R-006: gather has 17 fields after empty-k4-on-other"
                     "-ranks scenario");
         TEST_ASSERT(g.field_arrays.size() == 17u,
                     "R-006: 17 field arrays");
         TEST_ASSERT(g.field_arrays[12].size() == static_cast<size_t>(nranks),
                     "R-006: k4 field has one entry per rank");
         TEST_ASSERT(std::abs(g.field_arrays[12][0] - 7.0) < 1e-12,
                     "R-006: rank-0 contribution is 7.0 (populated k4)");
         bool tail_zero = true;
         for (int r = 1; r < nranks; ++r)
         {
            if (std::abs(g.field_arrays[12][r]) > 1e-12) { tail_zero = false; }
         }
         TEST_ASSERT(tail_zero,
                     "R-006: ranks 1..n-1 contribute 0.0 (post-fix loop "
                     "skipped OOB reads)");
      }
   }

   // -----------------------------------------------------------------
   // R-002 regression note: a death-test for field-name reordering
   // would require fork()ing or expecting MFEM_ABORT, which the
   // existing harness does not support.  The hash check is
   // exercised by the previous gathers (matching names ⇒ matching
   // hashes ⇒ no abort).  Coverage is best added when the death-test
   // harness lands.
   // -----------------------------------------------------------------

   if (rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << " / " << num_tests
                << " passed; " << num_failed << " failed ===\n";
   }

   const int local_failed = (rank == 0) ? num_failed : 0;
   int total_failed = 0;
   MPI_Allreduce(&local_failed, &total_failed, 1, MPI_INT, MPI_MAX,
                 MPI_COMM_WORLD);
   MPI_Finalize();
   return total_failed > 0 ? 1 : 0;
#else
   (void)argc; (void)argv;
   std::cout << "test_fault_surface_vtu_gather_mpi requires MFEM_USE_MPI\n";
   return 0;
#endif
}
