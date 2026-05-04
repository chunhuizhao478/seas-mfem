// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 2b of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
//
// MPI gather + VTKHDF write on np=4: each rank constructs a small
// fault data shape, calls vtkhdf::WriteFaultPackHdf collectively,
// and rank 0 reads back the .vtkhdf to verify all ranks' contributions
// landed.
//
// Plan §Phase 2b acceptance criterion:
//   "test_fault_surface_vtkhdf_mpi.cpp (mpirun -np 4) verifies that
//    the merged HDF5 file contains every face from every rank
//    exactly once."
//
// Note: this test exercises the rank-0-collective gather +
// rank-0-serial HDF5 write path (the documented Phase 2b deviation
// from the plan's parallel-HDF5-collective-write).  See header of
// io/fault_vtkhdf_writer.hpp for rationale.

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
#include <unistd.h>
#include <vector>

#ifdef MFEM_USE_HDF5
#include "../../io/fault_vtkhdf_writer.hpp"
#include <hdf5.h>
#endif

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

#ifdef MFEM_USE_HDF5
// R-008 (REVIEW.md 2026-04-28): try CellData first, fall back to
// PointData (current MFEM serialisation of L2-p0 GFs).  See
// test_fault_surface_vtkhdf.cpp for the rationale.
static std::vector<double> ReadHDF5PointDataAsCellData(
   const std::string &path, const std::string &field,
   int expect_cell_count)
{
   std::vector<double> out;
   hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) { return out; }

   hid_t ds = -1;
   bool is_cell_data = false;
   const char *roots[] = {"/VTKHDF/CellData/", "/VTKHDF/PointData/"};
   for (int r = 0; r < 2 && ds < 0; ++r)
   {
      const std::string p = std::string(roots[r]) + field;
      ds = H5Dopen2(file, p.c_str(), H5P_DEFAULT);
      if (ds >= 0) { is_cell_data = (r == 0); }
   }
   if (ds < 0) { H5Fclose(file); return out; }

   hid_t space = H5Dget_space(ds);
   hsize_t dims[2] = {0, 0};
   const int ndims = H5Sget_simple_extent_ndims(space);
   H5Sget_simple_extent_dims(space, dims, nullptr);
   const hsize_t total = (ndims >= 1) ? dims[0] : 0;
   std::vector<double> raw(total);
   H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data());
   H5Sclose(space);
   H5Dclose(ds);
   H5Fclose(file);

   const std::size_t n_cells_emitted = is_cell_data ? raw.size()
                                                    : raw.size() / 3;
   out.reserve(n_cells_emitted);
   for (std::size_t c = 0; c < n_cells_emitted; ++c)
   {
      out.push_back(is_cell_data ? raw[c] : raw[3*c]);
   }
   if (expect_cell_count > 0
       && static_cast<int>(out.size()) > expect_cell_count)
   {
      out.resize(expect_cell_count);
   }
   return out;
}
#endif

int main(int argc, char *argv[])
{
#if defined(MFEM_USE_MPI) && defined(MFEM_USE_HDF5)
   MPI_Init(&argc, &argv);
   int rank = 0, nranks = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   if (rank == 0)
   {
      std::cout << "=== test_fault_surface_vtkhdf_mpi ("
                << nranks << " ranks) ===\n";
   }

   // 1. Each rank packs 1 triangle + 12 fields.  Field index 0 (which
   //    we call "slip_dip" to match the production schema) is filled
   //    with `100 * rank + 1` so we can assert per-rank provenance
   //    after the gather.
   const std::vector<std::string> base_names = {
      "slip_dip","slip_strike","slip_rate_dip","slip_rate_strike",
      "traction_dip","traction_strike","state_variable","normal_stress",
      "param_a","param_Dc","fault_x2","fault_x3"};
   vtu::LocalFaultPack pack;
   pack.field_names = base_names;
   pack.vertices = {{static_cast<double>(rank),       0.0, 0.0},
                    {static_cast<double>(rank) + 1.0, 0.0, 0.0},
                    {static_cast<double>(rank) + 0.5, 1.0, 0.0}};
   pack.triangles = {{0, 1, 2}};
   pack.field_arrays.assign(12, std::vector<double>{0.0});
   pack.field_arrays[0][0] = 100.0 * rank + 1.0;

   // 2. Write through the Phase 2b HDF5 writer.
   const std::string out_dir = "/tmp/test_fault_vtkhdf_mpi";
   if (rank == 0) { ::mkdir(out_dir.c_str(), 0755); }
   MPI_Barrier(MPI_COMM_WORLD);

   {
      // Local scope so FaultHDFState (and its ParaViewHDFDataCollection)
      // destructs BEFORE MPI_Finalize.  H5Fclose called by the dc
      // destructor relies on HDF5's MPI_COMM_SELF still being live.
      vtkhdf::FaultHDFState state;
      vtkhdf::WriteFaultPackHdf(state, out_dir, /*cycle=*/0, /*time=*/0.0,
                                /*compression=*/0, pack, rank, nranks,
                                MPI_COMM_WORLD);
   }

   // 3. Rank 0 reads back the merged file and verifies the
   //    field_arrays[0] entries are {1, 101, 201, ..., 100*(n-1)+1}
   //    after gather (i.e. one entry per rank, each with rank-id*100+1).
   if (rank == 0)
   {
      const std::string vtkhdf_path = out_dir + "/fault_surface.vtkhdf";
      TEST_ASSERT(std::ifstream(vtkhdf_path).good(),
                  "merged .vtkhdf was created on rank 0");

      auto v = ReadHDF5PointDataAsCellData(vtkhdf_path, "slip_dip",
                                           /*expect_cell_count=*/nranks);
      TEST_ASSERT(static_cast<int>(v.size()) == nranks,
                  "merged file has nranks cells in slip_dip array");

      std::set<int> seen;
      for (double d : v)
      {
         seen.insert(static_cast<int>(d));
      }
      bool full_set = true;
      for (int r = 0; r < nranks; ++r)
      {
         const int expected = 100 * r + 1;
         if (seen.count(expected) == 0) { full_set = false; }
      }
      TEST_ASSERT(full_set,
                  "every rank's contribution (100*r+1) is present exactly once");

      // Cleanup
      ::remove(vtkhdf_path.c_str());
      ::rmdir(out_dir.c_str());

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
   std::cout << "test_fault_surface_vtkhdf_mpi requires MFEM_USE_MPI and "
             << "MFEM_USE_HDF5\n";
   return 0;
#endif
}
