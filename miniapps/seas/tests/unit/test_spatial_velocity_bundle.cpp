// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_velocity_bundle.cpp — Phase 2 of
// spatial_dynamic_rupture_plan.md (rev-3).
//
// Plan §Phase 2 §Acceptance: 4 tests — each CVM model resolves the
// right path; override works; missing file aborts cleanly with the
// resolved path printed.

#include "mfem.hpp"

#include "../../spatial/code/spatial_velocity.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
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

// V-1  CVMH model resolves the correct path (no mesh_tag suffix)
static void V_1_cvmh_path()
{
   std::cout << "\n[V-1] CVMH model resolves to .../cvmh/velocity_safs.h5\n";
   VelocitySpec s;
   s.model = VelocityModel::CVMH;
   s.dataset_root = "/some/root";
   const std::string p = ResolveSpatialVelocitySidecarPath(s);
   TEST_ASSERT(p == "/some/root/velocity/results/cvmh/velocity_safs.h5",
               "cvmh path resolution");
   TEST_ASSERT(p.find("1000m") == std::string::npos,
               "no mesh_tag suffix in CVMH filename");
}

// V-2  CVMS_4_26_M01 + MultiscaleStatewise
static void V_2_other_models()
{
   std::cout << "\n[V-2] CVMS_4_26_M01 + MultiscaleStatewise paths\n";
   VelocitySpec s;
   s.dataset_root = "/r";

   s.model = VelocityModel::CVMS_4_26_M01;
   TEST_ASSERT(ResolveSpatialVelocitySidecarPath(s)
               == "/r/velocity/results/cvm_s4.26.m01/velocity_safs.h5",
               "cvm_s4.26.m01 path resolution");

   s.model = VelocityModel::MultiscaleStatewise;
   TEST_ASSERT(ResolveSpatialVelocitySidecarPath(s)
               == "/r/velocity/results/multiscale_statewise_cvm/velocity_safs.h5",
               "multiscale_statewise path resolution");
}

// V-3  override_path bypasses model resolution
static void V_3_override_path()
{
   std::cout << "\n[V-3] override_path bypasses model resolution\n";
   VelocitySpec s;
   s.model = VelocityModel::CVMH;
   s.dataset_root = "/this/should/be/ignored";
   s.override_path = "/elsewhere/explicit.h5";
   TEST_ASSERT(ResolveSpatialVelocitySidecarPath(s)
               == "/elsewhere/explicit.h5",
               "override path used verbatim");
}

namespace
{
// Forked-child abort tester (MFEM_USE_EXCEPTIONS=NO in this build).
bool RunInChild(const std::function<void()>& body)
{
   ::fflush(stdout); ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0)
   {
      std::cerr << "fork() failed: " << std::strerror(errno) << "\n";
      return false;
   }
   if (pid == 0)
   {
      ::freopen("/dev/null", "w", stderr);
      try { body(); } catch (...) { ::_exit(1); }
      ::_exit(0);
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}
}  // namespace

// V-4  missing file aborts (validator fires)
static void V_4_missing_file_aborts()
{
   std::cout << "\n[V-4] missing sidecar file aborts\n";
   const bool aborted = RunInChild([]()
   {
      VelocitySpec s;
      s.model = VelocityModel::CVMH;
      s.dataset_root = "/definitely/does/not/exist/xyzzy";
      mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
                          1, 1, 1, mfem::Element::HEXAHEDRON,
                          1.0, 1.0, 1.0);
      (void)LoadSpatialVelocityBundle(s, mesh);
   });
   TEST_ASSERT(aborted, "missing file must abort");
}

int main(int argc, char** argv)
{
   (void)argc; (void)argv;
   std::cout << "Running Phase 2 test_spatial_velocity_bundle\n";
   V_1_cvmh_path();
   V_2_other_models();
   V_3_override_path();
   V_4_missing_file_aborts();

   std::cout << "\n========================================\n";
   std::cout << "Phase 2 test_spatial_velocity_bundle: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
