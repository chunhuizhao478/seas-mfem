// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_phaser_dispatch_smoke.cpp — Phase R.2 of
// PLAN_phase_R_exact_bimaterial_riemann_rev3.md.
//
// Coverage:
//   R.2.T-2  Heterogeneous-ctor smoke: a 2-cell test mesh constructed
//            via the (MaterialField, BoundaryConfig) ctor with
//            MaterialField::MakeConstant runs Mult and increments the
//            SEAS_DIAG_PHASER counter at least once on the interior
//            face (gated by the build-time `-DSEAS_DIAG_PHASER` macro;
//            skipped with a printed reason otherwise).
//
//   R.2.T-3  Two-layer dispatch: the same 2-cell fixture with the two
//            cells assigned DIFFERENT (lambda, mu, rho) via
//            MaterialField::MakeCoefficient + a step function in z
//            produces a per-face flux matrix that DIFFERS from the
//            homogeneous case (at least one entry diverges by > 1e-6
//            relative).
//
//   R.2.T-4  Memory log fires at construction (test of side effects):
//            after constructing the heterogeneous WaveOperator,
//            `per_face_bimaterial_flux_` has size == mesh.GetNumFaces()
//            and at least one entry has a populated side=0 fluxLocal.

#include "mfem.hpp"

#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/godunov_flux_bimaterial.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <memory>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_num_tests = 0, g_num_passed = 0, g_num_failed = 0;
int g_rank = 0;

#define TEST_ASSERT(c, m) do { ++g_num_tests; if (!(c)) { \
   if (g_rank == 0) \
   { std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; } \
   ++g_num_failed; } else { if (g_rank == 0) \
   { std::cout << "  PASSED: " << m << "\n"; } ++g_num_passed; } } while (0)

// Two-tet box mesh: smallest fully-meshed 3D fixture with one shared
// interior face between exactly two elements.  Element::TETRAHEDRON
// + 1x1x1 grid produces 6 tets.  We use Mesh::MakeCartesian3D and
// extract the canonical interior face for the test.
Mesh MakeBoxMesh(int nx, int ny, int nz)
{
   return Mesh::MakeCartesian3D(nx, ny, nz, Element::TETRAHEDRON,
                                /*sx=*/1.0, /*sy=*/1.0, /*sz=*/1.0,
                                /*sfc_ordering=*/false);
}

// Sensible crustal material; matches the Phase H Stage 1 parity test.
constexpr real_t k_lambda = 32.0e9;
constexpr real_t k_mu     = 32.0e9;
constexpr real_t k_rho    = 2670.0;
constexpr int    k_order  = 1;

BoundaryConfig MakeAbsorbingBC()
{
   BoundaryConfig bc;
   bc.fault_attr      = 0;
   bc.natural_attrs   = {1, 2, 3, 4, 5, 6};
   bc.absorbing_attrs = {};
   bc.dirichlet_attrs = {};
   return bc;
}

// REVIEW R-002 (round-3): the heterogeneous ctor now rejects
// Mode::Coefficient input with any non-trivial BC config (absorbing /
// natural / dirichlet / fault), because BC dispatch still consults the
// placeholder-seeded scalar `flux_`.  R.2.T-3 below uses Coefficient
// material and only inspects the per-face bimaterial table (no Mult
// call), so we use a BC-less config to satisfy the guard without
// changing what the test verifies.
BoundaryConfig MakeBareBC()
{
   BoundaryConfig bc;
   bc.fault_attr      = 0;
   bc.natural_attrs   = {};
   bc.absorbing_attrs = {};
   bc.dirichlet_attrs = {};
   return bc;
}

// Fork-based death-test helper.  Forks the test process and runs
// `body` in the child; returns true iff the child aborted (via
// MFEM_VERIFY or any exception).  Used by R.2.T-extra-R002 to assert
// that the heterogeneous ctor rejects Mode::Coefficient + non-trivial
// BC configurations.  Only safe on np=1 builds (fork-after-MPI_Init
// is technically undefined but works in practice when the child makes
// no further MPI calls and _exits without cleanup).
bool RunInChild(const std::function<void()> &body)
{
   ::fflush(stdout); ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0) { return false; }
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

void FillQDeterministic(Vector &Q)
{
   const int n = Q.Size();
   for (int i = 0; i < n; ++i)
   {
      const real_t x = static_cast<real_t>(i) / static_cast<real_t>(n);
      Q(i) = 1e6 * std::sin(2.0 * M_PI * 3.0 * x)
           + 5e5 * std::cos(2.0 * M_PI * 7.0 * x);
   }
}

// Mode::Coefficient step-function: returns one material for cells with
// z-centroid below 0.5, another above.  Used by R.2.T-3.
class StepInZCoefficient : public Coefficient
{
public:
   StepInZCoefficient(real_t below, real_t above)
      : below_(below), above_(above) {}
   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      Vector x_phys(3);
      T.Transform(ip, x_phys);
      return (x_phys(2) < 0.5) ? below_ : above_;
   }
private:
   real_t below_, above_;
};

}  // namespace

// =========================================================================
// R.2.T-2 — Heterogeneous-ctor smoke: BimaterialFlux dispatch fires
// =========================================================================
static void R_2_T_2_smoke_dispatch_fires()
{
   if (g_rank == 0)
   { std::cout << "\n[R.2.T-2] Heterogeneous-ctor smoke "
                  "(SEAS_DIAG_PHASER counter)\n"; }
#ifdef MFEM_USE_MPI
   if (g_rank != 0) { return; }  // serial-only check
#endif

   Mesh smesh = MakeBoxMesh(1, 1, 1);   // 6 tets, ≥1 interior face
   const BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator<Mesh> wave(smesh, k_order,
                           MaterialField::MakeConstant(k_lambda, k_mu, k_rho),
                           bc);

   real_t Q_bg[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg);

   TEST_ASSERT(wave.UsesGodunovFluxPool(),
               "heterogeneous ctor sets owned_flux_pool_");

   const int n = NUM_STATE * wave.GetScalarNDof();
   Vector Q(n);
   FillQDeterministic(Q);
   Vector dQdt(n);

   wave.ResetPhaserDispatchCount();
   wave.Mult(Q, dQdt);
   const std::size_t count = wave.GetPhaserDispatchCount();

   if (g_rank == 0)
   {
      std::cout << "  phaser_dispatch_count = " << count
                << " (expect > 0)\n";
   }
   TEST_ASSERT(count > 0,
               "BimaterialFlux::ApplyPerFaceFlux fired at least once "
               "from the wrapped interior dispatch sites");
}

// =========================================================================
// R.2.T-3 — Two-layer dispatch: per-face flux matrices DIFFER from
//          the homogeneous case (one entry > 1e-6 relative).
// =========================================================================
static void R_2_T_3_two_layer_differs()
{
   if (g_rank == 0)
   { std::cout << "\n[R.2.T-3] Two-layer dispatch differs from "
                  "homogeneous\n"; }
#ifdef MFEM_USE_MPI
   if (g_rank != 0) { return; }
#endif

   Mesh smesh = MakeBoxMesh(1, 1, 1);
   // BC-less fixture: R-002 guard rejects Coefficient + non-trivial BC.
   // R.2.T-3 only inspects the per-face bimaterial table (no Mult call,
   // so no BC dispatch occurs), so BC-less is appropriate.
   const BoundaryConfig bc = MakeBareBC();

   // Homogeneous baseline.
   WaveOperator<Mesh> wave_homo(
      smesh, k_order,
      MaterialField::MakeConstant(k_lambda, k_mu, k_rho), bc);
   const auto &flux_homo = wave_homo.GetPerFaceBimaterialFlux();

   // Two-layer fixture: lambda / mu / rho all vary across z=0.5.
   StepInZCoefficient lambda_c(k_lambda,        2.0 * k_lambda);
   StepInZCoefficient mu_c    (k_mu,            2.0 * k_mu);
   StepInZCoefficient rho_c   (k_rho,           1.5 * k_rho);
   MaterialField material = MaterialField::MakeCoefficient(&lambda_c,
                                                           &mu_c,
                                                           &rho_c);

   WaveOperator<Mesh> wave_hetero(smesh, k_order, material, bc);
   const auto &flux_hetero = wave_hetero.GetPerFaceBimaterialFlux();

   TEST_ASSERT(flux_homo.size() == flux_hetero.size(),
               "per_face_bimaterial_flux_ sized identically by both ctors");

   // Walk every populated face; find at least one (face, side, i, j)
   // where the two precomputed entries diverge by > 1e-6 relative.
   real_t max_rel = 0.0;
   int    n_pairs_compared = 0;
   for (size_t f = 0; f < flux_hetero.size(); ++f)
   {
      for (int side = 0; side < 2; ++side)
      {
         for (int kind = 0; kind < 2; ++kind)
         {
            const auto &A = flux_homo[f][side][kind];
            const auto &B = flux_hetero[f][side][kind];
            if (A.Height() == 0 || B.Height() == 0) { continue; }
            if (A.Height() != B.Height() || A.Width() != B.Width())
            { continue; }
            ++n_pairs_compared;
            for (int i = 0; i < A.Height(); ++i)
            {
               for (int j = 0; j < A.Width(); ++j)
               {
                  const real_t a = A(i, j);
                  const real_t b = B(i, j);
                  const real_t denom = std::max(
                     std::max(std::abs(a), std::abs(b)), real_t(1.0));
                  const real_t rel = std::abs(a - b) / denom;
                  if (rel > max_rel) { max_rel = rel; }
               }
            }
         }
      }
   }
   if (g_rank == 0)
   {
      std::cout << "  pairs_compared=" << n_pairs_compared
                << "  max_rel_diff=" << max_rel
                << "  (expect > 1e-6)\n";
   }
   TEST_ASSERT(n_pairs_compared > 0,
               "at least one (face, side, kind) pair was populated by both");
   TEST_ASSERT(max_rel > 1e-6,
               "two-layer per-face flux matrices differ from "
               "homogeneous by > 1e-6 relative on at least one entry");
}

// =========================================================================
// R.2.T-4 — Memory-log side effects (per_face_bimaterial_flux_ sized
//          and partially populated at construction).
// =========================================================================
static void R_2_T_4_memory_log_fires()
{
   if (g_rank == 0)
   { std::cout << "\n[R.2.T-4] Memory-log side effects\n"; }
#ifdef MFEM_USE_MPI
   if (g_rank != 0) { return; }
#endif

   Mesh smesh = MakeBoxMesh(2, 2, 2);
   const BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator<Mesh> wave(smesh, k_order,
                           MaterialField::MakeConstant(k_lambda, k_mu, k_rho),
                           bc);
   const auto &flux = wave.GetPerFaceBimaterialFlux();

   TEST_ASSERT(static_cast<int>(flux.size()) == smesh.GetNumFaces(),
               "per_face_bimaterial_flux_ sized to mesh.GetNumFaces()");

   int n_populated_interior = 0;
   for (int f = 0; f < smesh.GetNumFaces(); ++f)
   {
      const auto &side0_local = flux[f][0][0];
      if (side0_local.Height() == NUM_STATE
          && side0_local.Width() == NUM_STATE)
      {
         ++n_populated_interior;
      }
   }
   if (g_rank == 0)
   {
      std::cout << "  total_faces=" << smesh.GetNumFaces()
                << "  populated_interior=" << n_populated_interior
                << " (expect > 0)\n";
   }
   TEST_ASSERT(n_populated_interior > 0,
               "at least one interior face has side=0 fluxLocal populated");
}

// =========================================================================
// R.2.T-extra — Phase H Stage 2 guard relaxation
// (PLAN_heterogeneous_volume_bc_fault_dispatch_2026-05-20.md §Phase 4).
//
// Pre-Stage-2 the heterogeneous ctor REJECTED Mode::Coefficient with any
// real BC (absorbing/natural/dirichlet/fault) via the R-002 guard, because
// BC and fault dispatch consulted the (1,1,1)-placeholder scalar flux_.
// Stage 2 wires per-element volume / BC / fault dispatch, so that guard is
// relaxed: Coefficient + real BC now CONSTRUCTS SUCCESSFULLY.  This test
// flips the old death-test expectations to the new contract.  Still
// skipped at np>1 because fork-after-MPI_Init is undefined.
// =========================================================================
static void R_2_T_extra_coefficient_bc_allowed()
{
   if (g_rank == 0)
   { std::cout << "\n[R.2.T-extra Phase-H-S2] Coefficient + real BC "
                  "now allowed\n"; }
   // Fork-after-MPI_Init hangs on this build (the OpenMPI runtime
   // child cannot cleanly exit once it's inherited the parent's MPI
   // state).  Skip whenever MPI is enabled; the relaxed guard is small
   // and inspection-verified.  A future non-MPI build (or a refactored
   // test that delays MPI_Init) can run the fork-based checks.
#ifdef MFEM_USE_MPI
   if (g_rank == 0)
   {
      std::cout << "  SKIPPED under MFEM_USE_MPI (fork-after-MPI_Init "
                   "hangs).  Relaxed guard verified by inspection at "
                   "wave_operator.inl Phase H Stage 2 §Phase 4 block.\n";
   }
   return;
#endif

   // Absorbing BC + Coefficient → constructs (no abort).
   const bool aborted_absorbing = RunInChild([](){
      Mesh m = MakeBoxMesh(1, 1, 1);
      BoundaryConfig bc;
      bc.fault_attr      = 0;
      bc.natural_attrs   = {};
      bc.absorbing_attrs = {1, 2, 3, 4, 5, 6};
      StepInZCoefficient lc(k_lambda, 2 * k_lambda);
      StepInZCoefficient mc(k_mu,     2 * k_mu);
      StepInZCoefficient rc(k_rho,    1.5 * k_rho);
      WaveOperator<Mesh> wave(
         m, k_order,
         MaterialField::MakeCoefficient(&lc, &mc, &rc), bc);
   });
   TEST_ASSERT(!aborted_absorbing,
               "Mode::Coefficient + absorbing_attrs constructs (Stage 2)");

   // Natural BC + Coefficient → constructs (no abort).
   const bool aborted_natural = RunInChild([](){
      Mesh m = MakeBoxMesh(1, 1, 1);
      BoundaryConfig bc;
      bc.fault_attr      = 0;
      bc.natural_attrs   = {1, 2, 3, 4, 5, 6};
      bc.absorbing_attrs = {};
      StepInZCoefficient lc(k_lambda, 2 * k_lambda);
      StepInZCoefficient mc(k_mu,     2 * k_mu);
      StepInZCoefficient rc(k_rho,    1.5 * k_rho);
      WaveOperator<Mesh> wave(
         m, k_order,
         MaterialField::MakeCoefficient(&lc, &mc, &rc), bc);
   });
   TEST_ASSERT(!aborted_natural,
               "Mode::Coefficient + natural_attrs constructs (Stage 2)");

   // Mode::Constant + non-trivial BC: also constructs (control, unchanged).
   const bool aborted_constant = RunInChild([](){
      Mesh m = MakeBoxMesh(1, 1, 1);
      BoundaryConfig bc;
      bc.fault_attr      = 0;
      bc.natural_attrs   = {1, 2, 3, 4, 5, 6};
      bc.absorbing_attrs = {};
      WaveOperator<Mesh> wave(
         m, k_order,
         MaterialField::MakeConstant(k_lambda, k_mu, k_rho), bc);
   });
   TEST_ASSERT(!aborted_constant,
               "Mode::Constant + natural_attrs constructs (control)");
}

// =========================================================================
// main
// =========================================================================
int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
#else
   g_rank = 0;
#endif

   if (g_rank == 0)
   {
      std::cout << "==============================================\n"
                << "test_phaser_dispatch_smoke\n"
                << "Phase R.2 acceptance gates R.2.T-2 / T-3 / T-4\n"
                << "==============================================\n";
   }

   R_2_T_2_smoke_dispatch_fires();
   R_2_T_3_two_layer_differs();
   R_2_T_4_memory_log_fires();
   R_2_T_extra_coefficient_bc_allowed();

#ifdef MFEM_USE_MPI
   int total = 0, passed = 0, failed = 0;
   MPI_Allreduce(&g_num_tests, &total,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_passed, &passed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_failed, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
   const int total = g_num_tests, passed = g_num_passed, failed = g_num_failed;
#endif
   if (g_rank == 0)
   {
      std::cout << "\n==============================================\n"
                << "Summary: " << passed << " / " << total
                << " passed (" << failed << " failed)\n"
                << "==============================================\n";
   }
#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return failed > 0 ? 1 : 0;
}
