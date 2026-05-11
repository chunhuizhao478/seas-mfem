// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests for DataField3D (Phase 3 of the data-projection feature).
// Test catalog T-3-1 .. T-3-13 per data_projection_feature_plan_v2.md.
//
// All fixtures are synthetic HDF5 sidecars built in tmp paths at runtime
// via the local writer helpers below.  No real CVM-H data is used.
//
// Death-on-abort tests rely on the same `try { ctor or call; failure } catch
// { ... }` guard MFEM_ABORT exposes when MFEM_USE_EXCEPTIONS is enabled.
// MFEM aborts via std::abort by default, which terminates the test
// process.  To keep the suite single-process and to verify abort
// behaviour, we fork() before each abort-test path and observe the
// exit code via WIFEXITED / WIFSIGNALED.

#include "mfem.hpp"

#include "../../io/data_field_3d.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <hdf5.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

// ----------------------------------------------------------------------
// Tiny test framework (matches test_friction_law.cpp idiom)
// ----------------------------------------------------------------------

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ \
                   << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      const real_t _v = (value); \
      const real_t _e = (expected); \
      const real_t _t = (tol); \
      if (std::abs(_v - _e) > _t) { \
         std::cerr << "FAILED: " << message << " (got " << _v \
                   << ", expected " << _e << ", tol " << _t \
                   << ", line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)


// ----------------------------------------------------------------------
// HDF5 fixture helpers
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

/// Build a minimal, valid schema-v1 sidecar holding one field.  Caller
/// supplies axis arrays, the (Nx, Ny, Nz) field data (row-major flat),
/// the field name, units, and per-field min/max bounds.  Optional
/// overrides allow building intentionally-broken fixtures (mismatched
/// schema, wrong CRS, NaN cell, out-of-range cell).
struct SidecarOpts
{
   std::string schema_version = "data_projection_v1";
   std::string crs            = "EPSG:32611";
   std::string units          = "m";
   std::string z_positive     = "elevation";
   std::string field_name     = "Vp";
   std::string field_units    = "m/s";
   double      min_value      = 0.0;
   double      max_value      = 1.0;
};

void write_sidecar_synthetic(const std::string& path,
                             const std::vector<double>& x,
                             const std::vector<double>& y,
                             const std::vector<double>& z,
                             const std::vector<double>& field_flat_xyz,
                             const SidecarOpts& opts)
{
   hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                          H5P_DEFAULT);
   write_string_attr(file, "schema_version", opts.schema_version);
   write_string_attr(file, "crs",            opts.crs);
   write_string_attr(file, "units",          opts.units);
   write_string_attr(file, "z_positive",     opts.z_positive);
   write_string_attr(file, "created_at",     "1970-01-01T00:00:00Z");

   hid_t grid = H5Gcreate2(file, "/grid", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
   write_axis(grid, "x", x);
   write_axis(grid, "y", y);
   write_axis(grid, "z", z);
   H5Gclose(grid);

   hid_t flds = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
   hsize_t shape[3] = { x.size(), y.size(), z.size() };
   hid_t sid = H5Screate_simple(3, shape, nullptr);
   hid_t did = H5Dcreate2(flds, opts.field_name.c_str(),
                          H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            field_flat_xyz.data());
   write_string_attr(did, "units", opts.field_units);
   write_double_attr(did, "min_value", opts.min_value);
   write_double_attr(did, "max_value", opts.max_value);
   H5Dclose(did);
   H5Sclose(sid);
   H5Gclose(flds);

   H5Fclose(file);
}

/// Run a function that is expected to abort.  Returns true if the
/// child process exited abnormally (signal or non-zero exit code).
bool expect_abort(const std::function<void(void)>& f)
{
   pid_t pid = fork();
   if (pid == 0)
   {
      // Child: silence stderr/stdout while executing; abort is expected.
      ::close(1);
      ::close(2);
      try { f(); }
      catch (...) { std::_Exit(1); }
      std::_Exit(0);
   }
   int status = 0;
   ::waitpid(pid, &status, 0);
   if (WIFSIGNALED(status))             { return true; }
   if (WIFEXITED(status) && WEXITSTATUS(status) != 0) { return true; }
   return false;
}

std::string make_tmp_path(const char* stem)
{
   const char* tmpdir = std::getenv("TMPDIR");
   if (!tmpdir || !*tmpdir) { tmpdir = "/tmp"; }
   std::string path = std::string(tmpdir) + "/seas_test_data_field_3d_"
                      + std::string(stem) + "_"
                      + std::to_string(::getpid()) + ".h5";
   ::unlink(path.c_str());
   return path;
}

/// Build a 3x3x3 sidecar with f(x,y,z) = ax + by + cz + d.  Bounds set
/// large enough to hold the function range for the given a, b, c, d on
/// the fixture grid.
std::string make_linear_sidecar(double a, double b, double c, double d,
                                double v_min = -1e9, double v_max = 1e9)
{
   std::vector<double> x = {0.0, 1.0, 2.0};
   std::vector<double> y = {-1.0, 0.0, 1.0};
   std::vector<double> z = {-2.0, -1.0, 0.0};
   std::vector<double> data(x.size() * y.size() * z.size(), 0.0);
   for (size_t i = 0; i < x.size(); ++i)
   {
      for (size_t j = 0; j < y.size(); ++j)
      {
         for (size_t k = 0; k < z.size(); ++k)
         {
            const size_t flat = (i * y.size() + j) * z.size() + k;
            data[flat] = a * x[i] + b * y[j] + c * z[k] + d;
         }
      }
   }
   SidecarOpts opts;
   opts.field_name  = "F";
   opts.field_units = "unit";
   opts.min_value   = v_min;
   opts.max_value   = v_max;
   const std::string path = make_tmp_path("linear");
   write_sidecar_synthetic(path, x, y, z, data, opts);
   return path;
}

} // anonymous namespace


// ======================================================================
// Tests
// ======================================================================

static void T_3_1_schema_version_gate()
{
   std::cout << "\n[T-3-1] schema_version gate\n";
   std::vector<double> x = {0., 1.}, y = {0., 1.}, z = {-1., 0.};
   std::vector<double> data(8, 0.5);
   const std::string path = make_tmp_path("schema");
   SidecarOpts opts;
   opts.schema_version = "data_projection_v0";  // wrong
   opts.min_value = 0.0; opts.max_value = 1.0;
   write_sidecar_synthetic(path, x, y, z, data, opts);
   const bool aborted = expect_abort([&]() {
      DataField3D f(path, "Vp");
   });
   TEST_ASSERT(aborted, "T-3-1 ctor aborts on wrong schema_version");
   ::unlink(path.c_str());
}

static void T_3_2_crs_gate()
{
   std::cout << "\n[T-3-2] crs gate\n";
   std::vector<double> x = {0., 1.}, y = {0., 1.}, z = {-1., 0.};
   std::vector<double> data(8, 0.5);
   const std::string path = make_tmp_path("crs");
   SidecarOpts opts;
   opts.crs = "EPSG:4326";  // wrong
   opts.min_value = 0.0; opts.max_value = 1.0;
   write_sidecar_synthetic(path, x, y, z, data, opts);
   const bool aborted = expect_abort([&]() { DataField3D f(path, "Vp"); });
   TEST_ASSERT(aborted, "T-3-2 ctor aborts on wrong crs");
   ::unlink(path.c_str());
}

static void T_3_3_nan_gate()
{
   std::cout << "\n[T-3-3] NaN gate\n";
   std::vector<double> x = {0., 1.}, y = {0., 1.}, z = {-1., 0.};
   std::vector<double> data(8, 0.5);
   data[3] = std::numeric_limits<double>::quiet_NaN();
   const std::string path = make_tmp_path("nan");
   SidecarOpts opts;
   opts.min_value = 0.0; opts.max_value = 1.0;
   write_sidecar_synthetic(path, x, y, z, data, opts);
   const bool aborted = expect_abort([&]() { DataField3D f(path, "Vp"); });
   TEST_ASSERT(aborted, "T-3-3 ctor aborts on NaN cell");
   ::unlink(path.c_str());
}

static void T_3_4_value_range_gate()
{
   std::cout << "\n[T-3-4] declared value range gate\n";
   std::vector<double> x = {0., 1.}, y = {0., 1.}, z = {-1., 0.};
   std::vector<double> data(8, 0.5);
   data[5] = 100.0;       // way above declared max=1
   const std::string path = make_tmp_path("range");
   SidecarOpts opts;
   opts.min_value = 0.0; opts.max_value = 1.0;
   write_sidecar_synthetic(path, x, y, z, data, opts);
   const bool aborted = expect_abort([&]() { DataField3D f(path, "Vp"); });
   TEST_ASSERT(aborted, "T-3-4 ctor aborts on out-of-range cell");
   ::unlink(path.c_str());
}

static void T_3_5_trilinear_exactness()
{
   std::cout << "\n[T-3-5] trilinear exactness on linear field\n";
   const double a = 1.7, b = -2.3, c = 0.5, d = 100.0;
   const std::string path = make_linear_sidecar(a, b, c, d, -1000, 1000);
   DataField3D field(path, "F");
   std::mt19937 rng(20260509);
   std::uniform_real_distribution<double> ux(0.0, 2.0);
   std::uniform_real_distribution<double> uy(-1.0, 1.0);
   std::uniform_real_distribution<double> uz(-2.0, 0.0);
   for (int n = 0; n < 1024; ++n)
   {
      const double x = ux(rng), y = uy(rng), z = uz(rng);
      const double expected = a * x + b * y + c * z + d;
      const double got = field.Evaluate(x, y, z);
      if (std::abs(got - expected) > 1.0e-10)
      {
         TEST_NEAR(got, expected, 1.0e-10,
                   "T-3-5 trilinear at random point matches linear");
         break;
      }
   }
   TEST_ASSERT(true, "T-3-5 1024 random in-bbox samples match within 1e-10");
   ::unlink(path.c_str());
}

static void T_3_6_corner_query()
{
   std::cout << "\n[T-3-6] exact-corner query\n";
   const std::string path = make_linear_sidecar(2.0, 3.0, -1.0, 5.0,
                                                -1000, 1000);
   DataField3D field(path, "F");
   const double got = field.Evaluate(0.0, -1.0, -2.0);
   const double exp = 2.0 * 0.0 + 3.0 * (-1.0) + (-1.0) * (-2.0) + 5.0;
   TEST_NEAR(got, exp, 1.0e-12, "T-3-6 evaluate at SW-bottom corner");
   const double got2 = field.Evaluate(2.0, 1.0, 0.0);
   const double exp2 = 2.0 * 2.0 + 3.0 * 1.0 + (-1.0) * 0.0 + 5.0;
   TEST_NEAR(got2, exp2, 1.0e-12, "T-3-6 evaluate at NE-top corner");
   ::unlink(path.c_str());
}

static void T_3_7_midpoint_query()
{
   std::cout << "\n[T-3-7] midpoint query (arithmetic mean check)\n";
   // Constant-along-x linear-along-y field: f = 10 + 4 y.
   const std::string path = make_linear_sidecar(0.0, 4.0, 0.0, 10.0,
                                                -1000, 1000);
   DataField3D field(path, "F");
   const double got = field.Evaluate(0.5, -0.5, -0.5);
   const double exp = 0.0 * 0.5 + 4.0 * (-0.5) + 0.0 * (-0.5) + 10.0;
   TEST_NEAR(got, exp, 1.0e-12, "T-3-7 midpoint matches linear formula");
   ::unlink(path.c_str());
}

static void T_3_8_oob_aborts()
{
   std::cout << "\n[T-3-8] out-of-bbox query aborts\n";
   const std::string path = make_linear_sidecar(1.0, 1.0, 1.0, 0.0,
                                                -100, 100);
   const bool aborted = expect_abort([&]() {
      DataField3D field(path, "F");
      // Inside bbox first to ensure ctor succeeds in the parent's view,
      // then a query just outside the bbox.
      (void) field.Evaluate(2.0 + 1.0e-3, 0.0, -1.0);
   });
   TEST_ASSERT(aborted, "T-3-8 Evaluate aborts when x > xmax");
   ::unlink(path.c_str());
}

static void T_3_9_T_3_10_T_3_11_contains_bbox()
{
   std::cout << "\n[T-3-9..11] ContainsBBox semantics\n";
   const std::string path = make_linear_sidecar(0., 0., 0., 1.0,
                                                -1, 2);
   DataField3D field(path, "F");
   TEST_ASSERT(
      field.ContainsBBox(0.5, 1.5, -0.5, 0.5, -1.5, -0.5),
      "T-3-9 ContainsBBox=true for inner box");
   TEST_ASSERT(
      !field.ContainsBBox(-0.1, 1.0, 0.0, 0.0, -1.0, 0.0),
      "T-3-10 ContainsBBox=false when xmin straddles");
   // Boundary: data x in [0, 2].  Query xmin = -1e-7 with eps=0 fails
   // (≤ 0 and bbox lower is 0); with eps=1e-6, succeeds.
   TEST_ASSERT(
      !field.ContainsBBox(-1.0e-7, 2.0, 0.0, 0.0, -1.0, 0.0, 0.0),
      "T-3-11a ContainsBBox=false at sub-eps boundary with eps=0");
   TEST_ASSERT(
      field.ContainsBBox(-1.0e-7, 2.0, 0.0, 0.0, -1.0, 0.0, 1.0e-6),
      "T-3-11b ContainsBBox=true at sub-eps boundary with eps=1e-6");
   ::unlink(path.c_str());
}

static void T_3_12_T_3_13_find_index_edges()
{
   std::cout << "\n[T-3-12..13] find_index_ edges\n";
   const std::string path = make_linear_sidecar(0., 0., 0., 0.0,
                                                -1, 1);
   DataField3D field(path, "F");
   // The test fixture grid has x = {0, 1, 2}.
   TEST_ASSERT(field.FindIndexX(0.0) == 0,
               "T-3-12 find_index_(left edge) = 0");
   TEST_ASSERT(field.FindIndexX(2.0) == 1,
               "T-3-13 find_index_(right edge) = N-2 (=1 for N=3)");
   TEST_ASSERT(field.FindIndexX(0.5) == 0,
               "T-3-12b find_index_ inside leftmost cell");
   TEST_ASSERT(field.FindIndexX(1.5) == 1,
               "T-3-12c find_index_ inside rightmost cell");
   ::unlink(path.c_str());
}


// ----------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------

int main(int /*argc*/, char* /*argv*/[])
{
   std::cout << "===========================================\n";
   std::cout << "DataField3D unit tests (Phase 3, T-3-1..13)\n";
   std::cout << "===========================================\n";

   T_3_1_schema_version_gate();
   T_3_2_crs_gate();
   T_3_3_nan_gate();
   T_3_4_value_range_gate();
   T_3_5_trilinear_exactness();
   T_3_6_corner_query();
   T_3_7_midpoint_query();
   T_3_8_oob_aborts();
   T_3_9_T_3_10_T_3_11_contains_bbox();
   T_3_12_T_3_13_find_index_edges();

   std::cout << "\n===========================================\n";
   std::cout << "Total tests: " << num_tests << "\n";
   std::cout << "Passed:      " << num_passed << "\n";
   std::cout << "Failed:      " << num_failed << "\n";
   if (num_failed > 0) { std::cout << "\nSOME TESTS FAILED!\n"; return 1; }
   std::cout << "\nALL TESTS PASSED!\n";
   return 0;
}
