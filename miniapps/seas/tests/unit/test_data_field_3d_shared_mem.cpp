// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests for the MPI-3 shared-memory sidecar mode of DataField3D /
// StressField3D (PLAN_sidecar_mpi_shared_memory_2026-07-17.md §6).
//
// Run with 2+ ranks on ONE node (the Makefile target runs np=2 and np=4):
//   mpirun -np 2 ./seas_test_data_field_3d_shared_mem
//
// What is verified:
//   T-SM-1  Shared-mode DataField3D (node comm) evaluates BIT-IDENTICAL
//           to the classic per-rank path (MPI_COMM_NULL 3-arg ctor) at
//           every voxel corner AND at interior points, on every rank.
//   T-SM-2  The payload is genuinely the same across node-local ranks:
//           an elementwise MIN/MAX Allreduce over the corner samples
//           agrees bit-for-bit (min == max at every sample).
//   T-SM-3  OOBPolicy::Clamp composes with shared mode (edge-hold value
//           matches the per-rank Clamp reader exactly).
//   T-SM-4  A single-rank node comm (MPI_COMM_SELF) works — the
//           degenerate "1 rank per node" window.
//   T-SM-5  StressField3D's comm-forwarding ctor: shared-mode tensor
//           Evaluate bit-identical to the per-rank StressField3D.
//
// Fixtures are synthetic schema-v1 HDF5 files written by WORLD rank 0
// (same writer idiom as test_data_field_3d.cpp); the path is broadcast
// so every rank opens the same file.

#include "mfem.hpp"

#include "../../io/data_field_3d.hpp"
#include "../../io/stress_field_3d.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

// ----------------------------------------------------------------------
// Tiny test framework (matches test_data_field_3d.cpp idiom; counts are
// per-rank and reduced in main).
// ----------------------------------------------------------------------

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED (rank " << g_world_rank << "): " << message \
                   << " (line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         num_passed++; \
      } \
   } while (0)

static int g_world_rank = 0;

// ----------------------------------------------------------------------
// HDF5 fixture helpers (same schema-v1 writer as test_data_field_3d.cpp)
// ----------------------------------------------------------------------

namespace
{

void write_string_attr(hid_t loc, const char* name, const std::string& v)
{
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t atype = H5Tcopy(H5T_C_S1);
   H5Tset_size(atype, v.size() + 1);
   H5Tset_strpad(atype, H5T_STR_NULLTERM);
   hid_t aid = H5Acreate2(loc, name, atype, aspace, H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(aid, atype, v.c_str());
   H5Aclose(aid);
   H5Tclose(atype);
   H5Sclose(aspace);
}

void write_double_attr(hid_t loc, const char* name, double v)
{
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t aid = H5Acreate2(loc, name, H5T_NATIVE_DOUBLE, aspace,
                          H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(aid, H5T_NATIVE_DOUBLE, &v);
   H5Aclose(aid);
   H5Sclose(aspace);
}

void write_axis(hid_t group, const char* name, const std::vector<double>& v)
{
   hsize_t dim = v.size();
   hid_t sid = H5Screate_simple(1, &dim, nullptr);
   hid_t did = H5Dcreate2(group, name, H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
   H5Dclose(did);
   H5Sclose(sid);
}

void write_field(hid_t flds, const std::string& name,
                 const std::vector<double>& flat_xyz,
                 hsize_t nx, hsize_t ny, hsize_t nz,
                 double min_v, double max_v)
{
   hsize_t shape[3] = { nx, ny, nz };
   hid_t sid = H5Screate_simple(3, shape, nullptr);
   hid_t did = H5Dcreate2(flds, name.c_str(), H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            flat_xyz.data());
   write_string_attr(did, "units", "unit");
   write_double_attr(did, "min_value", min_v);
   write_double_attr(did, "max_value", max_v);
   H5Dclose(did);
   H5Sclose(sid);
}

// Deterministic, transposition-sensitive cell value: distinct at every
// (i, j, k) and different per field via `salt`.  Corner Evaluate()
// returns these values EXACTLY (trilinear is interpolating).
double cell_value(size_t i, size_t j, size_t k, double salt)
{
   return 0.5 + salt
          + 1.0    * static_cast<double>(i)
          + 10.0   * static_cast<double>(j)
          + 100.0  * static_cast<double>(k)
          + 0.001  * static_cast<double>((3 * i + 5 * j + 7 * k) % 11);
}

const std::vector<double>& axis_x()
{
   static const std::vector<double> v = {0.0, 1.0, 2.5, 4.0, 5.0};
   return v;
}
const std::vector<double>& axis_y()
{
   static const std::vector<double> v = {-2.0, -1.0, 0.5, 2.0};
   return v;
}
const std::vector<double>& axis_z()
{
   static const std::vector<double> v = {-3.0, -1.5, 0.0};
   return v;
}

std::vector<double> field_values(double salt)
{
   const size_t nx = axis_x().size();
   const size_t ny = axis_y().size();
   const size_t nz = axis_z().size();
   std::vector<double> data(nx * ny * nz, 0.0);
   for (size_t i = 0; i < nx; ++i)
      for (size_t j = 0; j < ny; ++j)
         for (size_t k = 0; k < nz; ++k)
         {
            data[(i * ny + j) * nz + k] = cell_value(i, j, k, salt);
         }
   return data;
}

/// Written by WORLD rank 0 only.  One single-field file ("F") and one
/// six-component stress file sharing the same grid.
void write_fixtures(const std::string& field_path,
                    const std::string& stress_path)
{
   const auto& x = axis_x();
   const auto& y = axis_y();
   const auto& z = axis_z();

   auto open_schema_file = [&](const std::string& path)
   {
      hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                             H5P_DEFAULT);
      write_string_attr(file, "schema_version", "data_projection_v1");
      write_string_attr(file, "crs",            "EPSG:32611");
      write_string_attr(file, "units",          "m");
      write_string_attr(file, "z_positive",     "elevation");
      write_string_attr(file, "created_at",     "1970-01-01T00:00:00Z");
      hid_t grid = H5Gcreate2(file, "/grid", H5P_DEFAULT, H5P_DEFAULT,
                              H5P_DEFAULT);
      write_axis(grid, "x", x);
      write_axis(grid, "y", y);
      write_axis(grid, "z", z);
      H5Gclose(grid);
      return file;
   };

   // Values span ~[0, 320]; give generous bounds.
   const double lo = -1.0e3, hi = 1.0e3;

   {
      hid_t file = open_schema_file(field_path);
      hid_t flds = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT,
                              H5P_DEFAULT);
      write_field(flds, "F", field_values(/*salt=*/0.0),
                  x.size(), y.size(), z.size(), lo, hi);
      H5Gclose(flds);
      H5Fclose(file);
   }
   {
      hid_t file = open_schema_file(stress_path);
      hid_t flds = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT,
                              H5P_DEFAULT);
      const char* comps[6] = { "sigma_xx", "sigma_yy", "sigma_zz",
                               "sigma_xy", "sigma_yz", "sigma_xz" };
      for (int c = 0; c < 6; ++c)
      {
         write_field(flds, comps[c], field_values(/*salt=*/1000.0 * (c + 1)),
                     x.size(), y.size(), z.size(), lo, 1.0e4);
      }
      H5Gclose(flds);
      H5Fclose(file);
   }
}

/// Broadcast a path chosen by rank 0 so every rank opens the SAME file
/// (make-tmp-by-pid would differ per rank).
std::string bcast_path(const std::string& rank0_path, MPI_Comm comm)
{
   int len = static_cast<int>(rank0_path.size());
   MPI_Bcast(&len, 1, MPI_INT, 0, comm);
   std::vector<char> buf(len + 1, '\0');
   if (g_world_rank == 0)
   {
      std::memcpy(buf.data(), rank0_path.c_str(), len);
   }
   MPI_Bcast(buf.data(), len, MPI_CHAR, 0, comm);
   return std::string(buf.data(), len);
}

std::string make_tmp_path(const char* stem)
{
   const char* tmpdir = std::getenv("TMPDIR");
   if (!tmpdir || !*tmpdir) { tmpdir = "/tmp"; }
   return std::string(tmpdir) + "/seas_test_df3d_sharedmem_"
          + std::string(stem) + "_" + std::to_string(::getpid()) + ".h5";
}

/// REVIEW.md R-002: fixture whose "F" carries ONE NaN cell (bounds
/// otherwise valid).  Loading it in shared mode must abort the WHOLE
/// job (node-rank-0 validation fires before the publish barrier and
/// MFEM_ABORT reaches MPI_Abort), never hang the peers at the barrier.
void write_nan_fixture(const std::string& path)
{
   const auto& x = axis_x();
   const auto& y = axis_y();
   const auto& z = axis_z();
   hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                          H5P_DEFAULT);
   write_string_attr(file, "schema_version", "data_projection_v1");
   write_string_attr(file, "crs",            "EPSG:32611");
   write_string_attr(file, "units",          "m");
   write_string_attr(file, "z_positive",     "elevation");
   write_string_attr(file, "created_at",     "1970-01-01T00:00:00Z");
   hid_t grid = H5Gcreate2(file, "/grid", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
   write_axis(grid, "x", x);
   write_axis(grid, "y", y);
   write_axis(grid, "z", z);
   H5Gclose(grid);
   hid_t flds = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
   std::vector<double> data = field_values(/*salt=*/0.0);
   data[data.size() / 2] = std::nan("");
   write_field(flds, "F", data, x.size(), y.size(), z.size(),
               -1.0e3, 1.0e3);
   H5Gclose(flds);
   H5Fclose(file);
}

/// All voxel-corner sample points plus a fixed set of interior points.
std::vector<std::array<real_t, 3>> sample_points()
{
   std::vector<std::array<real_t, 3>> pts;
   for (double xv : axis_x())
      for (double yv : axis_y())
         for (double zv : axis_z())
         {
            pts.push_back({ static_cast<real_t>(xv),
                            static_cast<real_t>(yv),
                            static_cast<real_t>(zv) });
         }
   // Interior (strictly inside every axis).
   pts.push_back({ real_t(0.3),  real_t(-1.7), real_t(-2.2) });
   pts.push_back({ real_t(1.9),  real_t(0.1),  real_t(-0.4) });
   pts.push_back({ real_t(3.1),  real_t(1.2),  real_t(-1.1) });
   pts.push_back({ real_t(4.99), real_t(1.99), real_t(-0.01) });
   return pts;
}

} // anonymous namespace


int main(int argc, char* argv[])
{
   MPI_Init(&argc, &argv);
   int world_size = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &g_world_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   // Node-local communicator — on the single-machine test runner this is
   // simply MPI_COMM_WORLD's ranks; key = world rank keeps the order.
   MPI_Comm node_comm = MPI_COMM_NULL;
   MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, g_world_rank,
                       MPI_INFO_NULL, &node_comm);
   // ------------------------------------------------------------------
   // REVIEW.md R-002 self-test mode (driven by the Makefile with a
   // NEGATED, time-bounded invocation): load a NaN fixture in shared
   // mode.  EXPECTED OUTCOME = the whole mpirun ABORTS (nonzero exit);
   // reaching the `return 0` below means the validation guard is broken.
   // NOTE: the fixture file intentionally leaks in TMPDIR (a few KB) —
   // there is no post-abort hook to unlink it.
   // ------------------------------------------------------------------
   if (argc > 1 && std::string(argv[1]) == "--nan-fixture-abort-child")
   {
      std::string p = make_tmp_path("nanchild");
      if (g_world_rank == 0) { write_nan_fixture(p); }
      p = bcast_path(p, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      {
         DataField3D bad(p, "F", OOBPolicy::Abort, node_comm);   // must abort
         std::cerr << "R-002 FAILURE: NaN shared-mode fixture did NOT abort "
                   << "(Evaluate(0,-1,-1.5) = "
                   << bad.Evaluate(real_t(0.0), real_t(-1.0), real_t(-1.5))
                   << ")\n";
      }
      // Reaching here == the guard is BROKEN.  Exit CLEANLY with 0 so the
      // Makefile's negated (`!`) invocation reports the failure loudly
      // (MPI_Abort with a zero code has implementation-defined mpirun
      // status and could mask it).
      MPI_Comm_free(&node_comm);
      MPI_Finalize();
      return 0;
   }

   int node_size = 0;
   MPI_Comm_size(node_comm, &node_size);
   if (g_world_rank == 0 && node_size < 2 && world_size >= 2)
   {
      std::cerr << "WARNING: np=" << world_size << " but node comm size is "
                << node_size << " — ranks do not share a node; the "
                << "cross-rank sharing assertions degenerate.\n";
   }

   // Fixtures: written once by world rank 0, path broadcast to all.
   std::string field_path  = make_tmp_path("field");
   std::string stress_path = make_tmp_path("stress");
   if (g_world_rank == 0) { write_fixtures(field_path, stress_path); }
   field_path  = bcast_path(field_path,  MPI_COMM_WORLD);
   stress_path = bcast_path(stress_path, MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);

   const auto pts = sample_points();

   // ------------------------------------------------------------------
   // T-SM-1: shared-mode reader == per-rank reader, bit-for-bit.
   // ------------------------------------------------------------------
   {
      DataField3D ref(field_path, "F", OOBPolicy::Abort);
      DataField3D shared(field_path, "F", OOBPolicy::Abort, node_comm);

      TEST_ASSERT(shared.NumX() == ref.NumX() &&
                  shared.NumY() == ref.NumY() &&
                  shared.NumZ() == ref.NumZ(),
                  "T-SM-1 axis sizes match");
      bool all_equal = true;
      for (const auto& p : pts)
      {
         const real_t a = ref.Evaluate(p[0], p[1], p[2]);
         const real_t b = shared.Evaluate(p[0], p[1], p[2]);
         if (!(a == b)) { all_equal = false; }
      }
      TEST_ASSERT(all_equal,
                  "T-SM-1 shared-mode Evaluate bit-identical to per-rank "
                  "path at all corner + interior samples");

      // ---------------------------------------------------------------
      // T-SM-2: the payload agrees bit-for-bit ACROSS node-local ranks:
      // elementwise MIN and MAX of every corner sample must coincide.
      // (Schema-v1 forbids NaN, so min==max <=> identical bits.)
      // ---------------------------------------------------------------
      std::vector<double> mine;
      mine.reserve(pts.size());
      for (const auto& p : pts)
      {
         mine.push_back(static_cast<double>(
                           shared.Evaluate(p[0], p[1], p[2])));
      }
      std::vector<double> vmin(mine), vmax(mine);
      MPI_Allreduce(MPI_IN_PLACE, vmin.data(), static_cast<int>(vmin.size()),
                    MPI_DOUBLE, MPI_MIN, node_comm);
      MPI_Allreduce(MPI_IN_PLACE, vmax.data(), static_cast<int>(vmax.size()),
                    MPI_DOUBLE, MPI_MAX, node_comm);
      bool cross_rank_equal = true;
      for (size_t i = 0; i < vmin.size(); ++i)
      {
         if (vmin[i] != vmax[i]) { cross_rank_equal = false; }
      }
      TEST_ASSERT(cross_rank_equal,
                  "T-SM-2 shared payload identical across node-local ranks");
   }

   // ------------------------------------------------------------------
   // T-SM-3: OOBPolicy::Clamp composes with shared mode.
   // ------------------------------------------------------------------
   {
      DataField3D ref(field_path, "F", OOBPolicy::Clamp);
      DataField3D shared(field_path, "F", OOBPolicy::Clamp, node_comm);
      // Outside the hull on every axis -> edge-hold value.
      const real_t a = ref.Evaluate(real_t(99.0), real_t(-99.0), real_t(9.0));
      const real_t b = shared.Evaluate(real_t(99.0), real_t(-99.0), real_t(9.0));
      TEST_ASSERT(a == b, "T-SM-3 Clamp edge-hold bit-identical");
   }

   // ------------------------------------------------------------------
   // T-SM-4: degenerate single-rank node comm (MPI_COMM_SELF).
   // ------------------------------------------------------------------
   {
      DataField3D ref(field_path, "F", OOBPolicy::Abort);
      DataField3D self_shared(field_path, "F", OOBPolicy::Abort,
                              MPI_COMM_SELF);
      bool all_equal = true;
      for (const auto& p : pts)
      {
         if (!(ref.Evaluate(p[0], p[1], p[2]) ==
               self_shared.Evaluate(p[0], p[1], p[2])))
         {
            all_equal = false;
         }
      }
      TEST_ASSERT(all_equal, "T-SM-4 MPI_COMM_SELF window == per-rank path");
   }

   // ------------------------------------------------------------------
   // T-SM-5: StressField3D comm forwarding (six shared windows).
   // ------------------------------------------------------------------
   {
      StressField3D ref(stress_path);
      StressField3D shared(stress_path, OOBPolicy::Abort, node_comm);
      bool all_equal = true;
      for (const auto& p : pts)
      {
         const DenseMatrix A = ref.Evaluate(p[0], p[1], p[2]);
         const DenseMatrix B = shared.Evaluate(p[0], p[1], p[2]);
         for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
            {
               if (!(A(r, c) == B(r, c))) { all_equal = false; }
            }
      }
      TEST_ASSERT(all_equal,
                  "T-SM-5 StressField3D shared-mode tensor bit-identical");
   }

   // Cleanup (windows freed by the dtors above, BEFORE MPI_Finalize).
   MPI_Barrier(MPI_COMM_WORLD);
   if (g_world_rank == 0)
   {
      ::unlink(field_path.c_str());
      ::unlink(stress_path.c_str());
   }
   MPI_Comm_free(&node_comm);

   int failed_global = 0;
   MPI_Allreduce(&num_failed, &failed_global, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
   int tests_global = 0, passed_global = 0;
   MPI_Allreduce(&num_tests, &tests_global, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
   MPI_Allreduce(&num_passed, &passed_global, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
   if (g_world_rank == 0)
   {
      std::cout << "\ntest_data_field_3d_shared_mem (np=" << world_size
                << "): " << passed_global << "/" << tests_global
                << " rank-local assertions passed, " << failed_global
                << " failed.\n";
   }

   MPI_Finalize();
   return failed_global == 0 ? 0 : 1;
}
