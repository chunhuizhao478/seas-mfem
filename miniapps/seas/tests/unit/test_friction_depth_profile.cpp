// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_friction_depth_profile.cpp — Phase 11b unit tests for the depth-varying
// rate-and-state a(z)/b(z) input.
//
// Plan: PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24, Phase 11b.
//
//   D1  PiecewiseLinear1D: exact at samples, linear at midpoints, flat-clamp
//       outside the range.
//   D2  PiecewiseLinear1D::Validate rejects <2 samples / non-monotonic.
//   D3  Two-CSV loader: columns are (value, depth_km); km->m scaling; the
//       value column is read FIRST.
//   D4  b(depth) = a(depth) - (a-b)(depth) on DIFFERENT depth grids.
//   D5  Loader rejects <2 rows / duplicate depth / non-positive a / 3-field row;
//       skips comments + blank lines; sorts unsorted rows.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas::spatial;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e << ", tol " << _t << ")\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

namespace
{

// Run `body` in a forked child; return true iff it aborted/exited non-zero.
// (MFEM_VERIFY/MFEM_ABORT -> abort() in the serial test.)
bool RunInChild(const std::function<void()>& body)
{
   ::fflush(stdout);
   ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0) { return false; }
   if (pid == 0)
   {
      ::freopen("/dev/null", "w", stderr);
      try { body(); }
      catch (...) { ::_exit(1); }
      ::_exit(0);
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}

// Write `content` to `path` (overwrite).  Returns path for convenience.
std::string write_file(const std::string& path, const std::string& content)
{
   std::ofstream ofs(path, std::ios::trunc);
   ofs << content;
   ofs.close();
   return path;
}

std::string tmp_path(const char* name)
{
   return std::string("/tmp/seas_test_depth_profile_") + name;
}

} // namespace

// =====================================================================
// D1: PiecewiseLinear1D interpolation + flat clamp.
// =====================================================================
static void D1_interp_and_clamp()
{
   std::cout << "\n[D1] PiecewiseLinear1D: samples, midpoints, flat clamp\n";
   PiecewiseLinear1D pl;
   pl.x = { 0.0, 10.0, 50.0 };
   pl.y = { 1.0,  3.0,  9.0 };
   pl.Validate();

   TEST_NEAR(pl(0.0),  1.0, 1e-15, "exact at first sample");
   TEST_NEAR(pl(10.0), 3.0, 1e-15, "exact at interior sample");
   TEST_NEAR(pl(50.0), 9.0, 1e-15, "exact at last sample");
   // Midpoints: linear.
   TEST_NEAR(pl(5.0),  2.0, 1e-15, "linear midpoint [0,10] -> 2.0");
   TEST_NEAR(pl(30.0), 6.0, 1e-15, "linear midpoint [10,50] -> 6.0");
   // Flat clamp outside.
   TEST_NEAR(pl(-100.0), 1.0, 1e-15, "flat clamp below first -> y.front()");
   TEST_NEAR(pl(1e6),    9.0, 1e-15, "flat clamp above last  -> y.back()");
}

// =====================================================================
// D2: Validate rejects degenerate inputs.
// =====================================================================
static void D2_validate_rejects()
{
   std::cout << "\n[D2] PiecewiseLinear1D::Validate rejects degenerate inputs\n";

   const bool one_sample = RunInChild([]() {
      PiecewiseLinear1D pl; pl.x = { 0.0 }; pl.y = { 1.0 }; pl.Validate();
   });
   TEST_ASSERT(one_sample, "Validate aborts on < 2 samples");

   const bool non_monotone = RunInChild([]() {
      PiecewiseLinear1D pl; pl.x = { 0.0, 5.0, 5.0 }; pl.y = { 1.0, 2.0, 3.0 };
      pl.Validate();
   });
   TEST_ASSERT(non_monotone, "Validate aborts on non-strictly-increasing x");

   const bool size_mismatch = RunInChild([]() {
      PiecewiseLinear1D pl; pl.x = { 0.0, 1.0 }; pl.y = { 1.0 }; pl.Validate();
   });
   TEST_ASSERT(size_mismatch, "Validate aborts on x/y size mismatch");
}

// =====================================================================
// D3: Two-CSV loader — (value, depth_km) columns, km->m scaling, value first.
// =====================================================================
static void D3_loader_columns_and_units()
{
   std::cout << "\n[D3] loader: (value, depth_km), km->m scaling, value first\n";
   FrictionDepthProfileSpec spec;
   // value FIRST, depth_km SECOND (matches the committed param_*.csv).
   spec.param_a_csv = write_file(tmp_path("a_d3.csv"),
                                 "0.010, 0\n0.030, 10\n0.150, 50\n");
   spec.param_a_minus_b_csv = write_file(tmp_path("amb_d3.csv"),
                                         "-0.009, 0\n0.130, 60\n");
   spec.depth_to_m = 1000.0;   // km
   FrictionDepthProfile1D prof = LoadFrictionDepthProfileCSVs(spec);

   // value-first: a at the surface is 0.010 (NOT 0, the depth field).
   TEST_NEAR(prof.a(0.0), 0.010, 1e-15, "a(0) = 0.010 (value column read first)");
   // km scaling: the 10 km knot lands at depth 10 000 m, value 0.030.
   TEST_NEAR(prof.a(10000.0), 0.030, 1e-15, "a(10 km) = 0.030 (km->m scaling)");
   // linear between 0 and 10 km: a(5 km) = 0.020.
   TEST_NEAR(prof.a(5000.0), 0.020, 1e-15, "a(5 km) linear -> 0.020");
   // flat clamp above the last a-knot (50 km).
   TEST_NEAR(prof.a(70000.0), 0.150, 1e-15, "a(70 km) flat-clamped -> 0.150");
}

// =====================================================================
// D4: b = a - (a-b) computed on DIFFERENT depth grids.
// =====================================================================
static void D4_b_from_a_minus_amb_diff_grids()
{
   std::cout << "\n[D4] b(z) = a(z) - (a-b)(z) on different depth grids\n";
   FrictionDepthProfileSpec spec;
   // a grid: 0, 10, 50 km ; amb grid: 0, 13, 60 km  (deliberately different).
   spec.param_a_csv = write_file(tmp_path("a_d4.csv"),
                                 "0.010, 0\n0.030, 10\n0.150, 50\n");
   spec.param_a_minus_b_csv = write_file(tmp_path("amb_d4.csv"),
                                         "-0.009, 0\n-0.009, 13\n0.130, 60\n");
   spec.depth_to_m = 1000.0;
   FrictionDepthProfile1D prof = LoadFrictionDepthProfileCSVs(spec);

   // Surface: a=0.010, a-b=-0.009 -> b = 0.019 (VW, a<b).
   TEST_NEAR(prof.b(0.0), 0.019, 1e-15, "b(0) = a - (a-b) = 0.019 (VW)");
   // 5 km: a = 0.020 (a-grid lerp), a-b = -0.009 (amb constant to 13 km) -> 0.029.
   TEST_NEAR(prof.b(5000.0), 0.029, 1e-15, "b(5 km) on differing grids = 0.029");
   // a-b crosses 0 between 13 and 60 km; below that, a-b>0 so a>b (VS).
   TEST_ASSERT(prof.b(60000.0) < prof.a(60000.0),
               "deep: a-b > 0 => b < a (velocity-strengthening)");
}

// =====================================================================
// D5: loader rejections + tolerant parsing.
// =====================================================================
static void D5_loader_rejections_and_tolerant_parse()
{
   std::cout << "\n[D5] loader rejects bad input; tolerates comments/blank/unsorted\n";

   // < 2 data rows.
   const bool too_few = RunInChild([]() {
      FrictionDepthProfileSpec spec;
      spec.param_a_csv = write_file(tmp_path("a_few.csv"), "0.010, 0\n");
      spec.param_a_minus_b_csv = write_file(tmp_path("amb_few.csv"),
                                            "-0.009, 0\n0.130, 60\n");
      (void) LoadFrictionDepthProfileCSVs(spec);
   });
   TEST_ASSERT(too_few, "loader aborts on < 2 data rows");

   // Duplicate depth.
   const bool dup = RunInChild([]() {
      FrictionDepthProfileSpec spec;
      spec.param_a_csv = write_file(tmp_path("a_dup.csv"),
                                    "0.010, 5\n0.030, 5\n");
      spec.param_a_minus_b_csv = write_file(tmp_path("amb_dup.csv"),
                                            "-0.009, 0\n0.130, 60\n");
      (void) LoadFrictionDepthProfileCSVs(spec);
   });
   TEST_ASSERT(dup, "loader aborts on duplicate depth");

   // Non-positive a value.
   const bool nonpos_a = RunInChild([]() {
      FrictionDepthProfileSpec spec;
      spec.param_a_csv = write_file(tmp_path("a_neg.csv"),
                                    "0.010, 0\n-0.001, 10\n");
      spec.param_a_minus_b_csv = write_file(tmp_path("amb_neg.csv"),
                                            "-0.009, 0\n0.130, 60\n");
      (void) LoadFrictionDepthProfileCSVs(spec);
   });
   TEST_ASSERT(nonpos_a, "loader aborts on non-positive a value");

   // 3-field row.
   const bool three_fields = RunInChild([]() {
      FrictionDepthProfileSpec spec;
      spec.param_a_csv = write_file(tmp_path("a_3f.csv"),
                                    "0.010, 0, 99\n0.030, 10\n");
      spec.param_a_minus_b_csv = write_file(tmp_path("amb_3f.csv"),
                                            "-0.009, 0\n0.130, 60\n");
      (void) LoadFrictionDepthProfileCSVs(spec);
   });
   TEST_ASSERT(three_fields, "loader aborts on a 3-field row");

   // Tolerant parse: comments, blank lines, unsorted rows.
   FrictionDepthProfileSpec spec;
   spec.param_a_csv = write_file(tmp_path("a_tol.csv"),
      "# a(z) profile\n\n0.150, 50\n0.010, 0\n  0.030 , 10  \n");  // unsorted, comment, blank
   spec.param_a_minus_b_csv = write_file(tmp_path("amb_tol.csv"),
      "-0.009, 0\n0.130, 60\n");
   spec.depth_to_m = 1000.0;
   FrictionDepthProfile1D prof = LoadFrictionDepthProfileCSVs(spec);
   TEST_NEAR(prof.a(0.0),     0.010, 1e-15, "tolerant: a(0) after sort = 0.010");
   TEST_NEAR(prof.a(10000.0), 0.030, 1e-15, "tolerant: a(10 km) after sort = 0.030");
}

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase 11b test_friction_depth_profile\n";
   D1_interp_and_clamp();
   D2_validate_rejects();
   D3_loader_columns_and_units();
   D4_b_from_a_minus_amb_diff_grids();
   D5_loader_rejections_and_tolerant_parse();

   std::cout << "\n========================================\n";
   std::cout << "Phase 11b test_friction_depth_profile: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}
