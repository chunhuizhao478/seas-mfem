// Phase 4 (PLAN_mixed_flux_hetero_riemann.md): mixed-flux CFL de-rating on the
// matrix path uses the SAME shared factor as the scalar path, and central flux
// is forbidden under ADER on BOTH operators.
//
//   Test 4.1 (BUG-2) — for mode in {None, Adjacent, AllContinuous} x
//                       cfl_rk_aware_ in {false, true}, WaveOperator and
//                       BimaterialWaveOperator return the SAME MixedFluxCflFactor_()
//                       (None->1, Adjacent->0.6/0.9, AllContinuous->0.7/0.4).
//                       Confirms the single source of truth — no residual
//                       divergent 0.9/0.4 switch on the matrix path.
//   Test 4.2          — central+ADER (mixed_flux != None && !cfl_rk_aware_)
//                       ComputeMaxDt aborts on BOTH operators; with RK it does
//                       NOT abort (finite dt); mixed_flux=none under ADER does
//                       NOT abort (byte-exact contract).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../domain/boundary_config.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do {                                          \
   num_tests++;                                                              \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; }     \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: "       \
                                  << msg << "\n"; }                          \
} while (0)

namespace
{
// Expose the protected MixedFluxCflFactor_ for direct probing (Test 4.1).
template <typename MeshType>
class TestableWave : public WaveOperator<MeshType>
{
public:
   using WaveOperator<MeshType>::WaveOperator;
   using WaveOperator<MeshType>::MixedFluxCflFactor_;
};
template <typename MeshType>
class TestableBimat : public BimaterialWaveOperator<MeshType>
{
public:
   using BimaterialWaveOperator<MeshType>::BimaterialWaveOperator;
   using BimaterialWaveOperator<MeshType>::MixedFluxCflFactor_;  // inherited from base
};

// Row of `nx` unit hexes along x; interior face at x == fault_x tagged FAULT
// (attr 3), external faces FREE (attr 1).  (Same fixture as the dispatch test:
// Adjacent needs a fault; this gives one + fault-adjacent interior faces.)
Mesh BuildHexRowFaultMesh(int nx, int fault_x)
{
   Mesh mesh(3, (nx + 1) * 4, nx, 0);
   auto vid = [](int ix, int iy, int iz) { return ix * 4 + iy * 2 + iz; };
   for (int ix = 0; ix <= nx; ++ix)
      for (int iy = 0; iy < 2; ++iy)
         for (int iz = 0; iz < 2; ++iz)
         {
            real_t v[3] = { static_cast<real_t>(ix), static_cast<real_t>(iy),
                            static_cast<real_t>(iz) };
            mesh.AddVertex(v);
         }
   for (int e = 0; e < nx; ++e)
   {
      const int v[8] = {
         vid(e, 0, 0), vid(e + 1, 0, 0), vid(e + 1, 1, 0), vid(e, 1, 0),
         vid(e, 0, 1), vid(e + 1, 0, 1), vid(e + 1, 1, 1), vid(e, 1, 1)
      };
      mesh.AddHex(v, 1);
   }
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cx += mesh.GetVertex(fv[i])[0]; }
      cx /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cx - fault_x) < 1e-9)
         {
            mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3);
         }
         continue;
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Run `body` in a forked child; return true iff the child aborted (non-zero
// exit or signal).  Serial test (no MPI), so MFEM_ABORT -> abort() and fork is
// safe.  (Replicated from test_wave_operator.cpp.)
bool RunAbortsInChild(const std::function<void()> &body)
{
   std::fflush(stdout);
   std::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0) { return false; }
   if (pid == 0)
   {
      std::freopen("/dev/null", "w", stderr);
      body();          // expected to MFEM_ABORT
      ::_exit(0);      // reached only if no abort
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}

} // anonymous namespace

int main()
{
   std::cout << "\n=== Phase 4: mixed-flux CFL de-rating (shared factor + central+ADER guard) ===\n";

   Mesh mesh = BuildHexRowFaultMesh(/*nx=*/5, /*fault_x=*/3);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   const real_t lam = 32.0e9, mu = 32.0e9, rho = 2670.0;   // homogeneous

   TestableWave<Mesh>  w(mesh, order, lam, mu, rho, bc);
   TestableBimat<Mesh> b(mesh, order,
                         MaterialField::MakeConstant(lam, mu, rho), bc);

   // -----------------------------------------------------------------------
   // Test 4.1: MixedFluxCflFactor_() is the single source of truth.
   // -----------------------------------------------------------------------
   std::cout << "\n-- Test 4.1: MixedFluxCflFactor_ matches scalar (all modes x rk) --\n";
   struct Combo { MixedFluxMode mode; bool rk; real_t expected; const char *name; };
   const Combo combos[] = {
      {MixedFluxMode::None,          false, 1.0, "None,ADER"},
      {MixedFluxMode::None,          true,  1.0, "None,RK"},
      {MixedFluxMode::Adjacent,      false, 0.9, "Adjacent,ADER"},
      {MixedFluxMode::Adjacent,      true,  0.6, "Adjacent,RK"},
      {MixedFluxMode::AllContinuous, false, 0.4, "AllContinuous,ADER"},
      {MixedFluxMode::AllContinuous, true,  0.7, "AllContinuous,RK"},
   };
   for (const Combo &c : combos)
   {
      // Setting the mode is independent of cfl_rk_aware_ (the central+ADER
      // ComputeMaxDt guard is NOT hit here — we call only the pure factor).
      w.SetMixedFluxMode(c.mode); w.SetCflRkAware(c.rk);
      b.SetMixedFluxMode(c.mode); b.SetCflRkAware(c.rk);
      const real_t fw = w.MixedFluxCflFactor_();
      const real_t fb = b.MixedFluxCflFactor_();
      TEST_ASSERT(fw == c.expected,
                  std::string("Test 4.1: scalar MixedFluxCflFactor_ == expected for ")
                  + c.name);
      TEST_ASSERT(fb == fw,
                  std::string("Test 4.1: bimaterial == scalar (single source of truth) for ")
                  + c.name);
   }

   // -----------------------------------------------------------------------
   // Test 4.2: central+ADER aborts on BOTH operators; RK does not; none does not.
   // -----------------------------------------------------------------------
   std::cout << "\n-- Test 4.2: central+ADER aborts; RK / none do not --\n";

   b.SetMixedFluxMode(MixedFluxMode::Adjacent); b.SetCflRkAware(false);
   TEST_ASSERT(RunAbortsInChild([&]() { b.ComputeMaxDt(0.5); }),
               "Test 4.2: bimaterial ComputeMaxDt aborts on central+ADER "
               "(mixed flux requires RK)");

   w.SetMixedFluxMode(MixedFluxMode::Adjacent); w.SetCflRkAware(false);
   TEST_ASSERT(RunAbortsInChild([&]() { w.ComputeMaxDt(0.5); }),
               "Test 4.2: scalar ComputeMaxDt aborts on central+ADER "
               "(guard on BOTH operators)");

   // Control: Adjacent + RK -> no abort, finite positive dt.
   b.SetMixedFluxMode(MixedFluxMode::Adjacent); b.SetCflRkAware(true);
   TEST_ASSERT(!RunAbortsInChild([&]() { b.ComputeMaxDt(0.5); }),
               "Test 4.2: bimaterial central+RK does NOT abort");
   const real_t dt_rk = b.ComputeMaxDt(0.5);
   TEST_ASSERT(std::isfinite(dt_rk) && dt_rk > 0.0,
               "Test 4.2: bimaterial central+RK returns a finite positive dt");

   // Control: mixed_flux=none under ADER -> no abort (byte-exact contract).
   b.SetMixedFluxMode(MixedFluxMode::None); b.SetCflRkAware(false);
   const real_t dt_none = b.ComputeMaxDt(0.5);
   TEST_ASSERT(std::isfinite(dt_none) && dt_none > 0.0,
               "Test 4.2: mixed_flux=none under ADER does NOT abort");

   std::cout << "\n==============================================\n";
   std::cout << "Summary: " << num_passed << " / " << num_tests
             << " passed (" << num_failed << " failed)\n";
   std::cout << "==============================================\n";
   return (num_failed > 0) ? 1 : 0;
}
