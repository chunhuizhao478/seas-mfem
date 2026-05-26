// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_shared_fault_rs_psi_consistency_np2.cpp — R-020(b): a psi-EXPLICIT
// np=2 cross-rank consistency test for the shared-fault rate-and-state path.
//
// ============================================================================
// Why this test exists
// ============================================================================
// Plan R-004 / review R-020 mandate an np=2 check that the rate-and-state state
// variable psi is single-valued across a rank seam for a SHARED fault QP.  The
// fix-round claim that test_shared_fault_reconcile_cross_rank.cpp already covers
// this was overstated: that test asserts an AGGREGATE worst_rel (over 9 DOFData
// fields) and contains zero literal `psi` references.  This test is
// psi-EXPLICIT: it reads and compares DOFData::psi directly across the two ranks.
//
// Fixture: the planar 2-tet fault (y=0) partitioned with the −side on rank 0 and
// the +side on rank 1, so the fault face is SHARED (the R-1601 inline-EvaluateADER
// path).  Rate-and-state (TPV102 aging) is set up identically on both ranks, then
// d.psi is PERTURBED on rank 0 ONLY so the two ranks' shared-QP psi differ.  One
// AdvanceADER macro step is taken: its ComputeADERSharedFaceFluxRHS reconcile
// broadcasts payload[5]=psi across the seam (wave_operator.inl:5304/5327).  The
// test then gathers the shared-QP psi from both ranks and asserts they are
// BIT-IDENTICAL post-reconcile.
//
// Oracle:
//   - BEFORE the step: rank-0 psi (perturbed) != rank-1 psi (unperturbed) — the
//     perturbation took, so the assertion is non-vacuous.
//   - AFTER the step:  rank-0 psi == rank-1 psi (exact) — the reconcile makes
//     the shared-fault RS state single-valued across the seam.
//
// Scope: this is part (a) (cross-rank psi bit-identity).  Part (b) (the
// interior-vs-shared psi temporal-drift bound) needs the multi-step time loop
// plus an np=1 interior twin and remains the deferred close-out — the driver's
// rank-0 "Treat parallel RS results as PRELIMINARY" warning (R-020) gates it.
//
// Usage:  mpirun -np 2 ./seas_test_shared_fault_rs_psi_consistency_np2
//         Requires MFEM built with MPI.  No -DSEAS_TEST_INTERNAL needed: the
//         perturbation writes DOFData::psi directly (a public field) and relies
//         on the PRODUCTION reconcile (not the test-only injection hook).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
#define TEST_TRUE(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace
{
constexpr real_t kL     = 1000.0;
constexpr real_t kDt    = 1.0e-4;
constexpr int    kOrder = 1;

// Planar 2-tet fault (y=0); identical to the reconcile / rake-sweep fixtures.
Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {kL, 0.0, 0.0}, {0.0, 0.0, kL},
      {0.0,  kL, 0.0},   // V3 on +y side
      {0.0, -kL, 0.0},   // V4 on -y side
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);   // tet 0: − side
   mesh.AddTet(0, 1, 2, 3, 1);   // tet 1: + side
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-8)
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }  // fault attr 3
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);         // natural attr 1
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

#ifdef MFEM_USE_MPI
ParMesh MakePartitioned2Tet()
{
   Mesh serial_mesh = BuildTwoTetFaultMesh();
   std::vector<int> part(serial_mesh.GetNE(), 0);
   for (int e = 0; e < serial_mesh.GetNE(); e++)
   {
      real_t cy = 0; Array<int> ev; serial_mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++) { cy += serial_mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      part[e] = (cy < 0.0) ? 0 : 1;   // −side→rank0, +side→rank1
   }
   return ParMesh(MPI_COMM_WORLD, serial_mesh, part.data());
}

// Collect the physical coords of every fault QP (interior + shared) in the
// order InitializeFaultDOFs expects.  Mirrors the reconcile-test helper.
int CollectFaultQPCoords(WaveOperator<ParMesh> &wave, ParMesh &mesh, int order,
                         std::vector<Vector> &coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault FTR null");
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   if (nqp == 0 && shr_faces.Size() > 0)
   {
      auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[0]);
      MFEM_VERIFY(ftr, "shared fault FTR null");
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         coords.push_back(phys);
      }
   }
   for (int i = 0; i < shr_faces.Size(); i++)
   {
      auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         coords.push_back(phys);
      }
   }
   return (int_faces.Size() + shr_faces.Size()) * nqp;
}

// Sorted vector of the local fault QPs' psi (all QPs are shared in this fixture).
std::vector<double> SortedLocalPsi(const std::vector<DOFData> &dof)
{
   std::vector<double> v;
   v.reserve(dof.size());
   for (const auto &d : dof) { v.push_back(static_cast<double>(d.psi)); }
   std::sort(v.begin(), v.end());
   return v;
}

// Send rank-1's sorted psi to rank 0; on rank 0 return whether the two ranks'
// sorted psi multisets are BIT-IDENTICAL (same size + exact ==).  Returns true
// on rank 1 (no-op) so only rank 0's verdict matters.
bool RanksPsiBitIdentical(const std::vector<double> &local_sorted, int rank)
{
   if (rank == 1)
   {
      int n = static_cast<int>(local_sorted.size());
      MPI_Send(&n, 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
      if (n > 0)
      { MPI_Send(local_sorted.data(), n, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD); }
      return true;
   }
   // rank 0
   int n_peer = 0;
   MPI_Recv(&n_peer, 1, MPI_INT, 1, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   std::vector<double> peer(n_peer, 0.0);
   if (n_peer > 0)
   {
      MPI_Recv(peer.data(), n_peer, MPI_DOUBLE, 1, 101, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
   }
   if (static_cast<int>(local_sorted.size()) != n_peer) { return false; }
   for (int i = 0; i < n_peer; i++)
   {
      if (local_sorted[i] != peer[i]) { return false; }   // EXACT bit-identity
   }
   return true;
}
#endif  // MFEM_USE_MPI

}  // namespace

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (rank == 0)
   {
      std::cout << "\n=== shared-fault RS psi cross-rank consistency (R-020b) ==="
                << "\n  np=" << nprocs << ", planar 2-tet fault, "
                << "perturb rank-0 psi then one ADER-2 step\n";
   }
   if (nprocs != 2)
   {
      if (rank == 0) { std::cout << "  SKIPPED: requires np == 2\n"; }
      MPI_Finalize();
      return 77;
   }

   constexpr real_t perturb = 0.05;   // psi offset injected on rank 0 only

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};

   ParMesh pmesh = MakePartitioned2Tet();

   const real_t lambda = TPV102Params::lambda;
   const real_t mu     = TPV102Params::mu;
   const real_t rho    = TPV102Params::rho;
   const real_t cp     = TPV102Params::cp;
   const real_t cs     = TPV102Params::cs;

   WaveOperator<ParMesh> wave(pmesh, kOrder, lambda, mu, rho, bc);
   wave.SetFaultFrictionLaw(FaultFrictionLaw::RateAndState);

   FaultFaceFlux ff(rho, cp, cs);
   wave.SetFaultFlux(&ff);

   std::vector<Vector> coords;
   const int n_fault = CollectFaultQPCoords(wave, pmesh, kOrder, coords);
   int nqp = 0;
   {
      const Array<int> &ifc = wave.GetFaultInteriorFaces();
      const Array<int> &sfc = wave.GetFaultSharedFaces();
      const int nf = ifc.Size() + sfc.Size();
      if (nf > 0 && n_fault > 0) { nqp = n_fault / nf; }
   }

   std::vector<DOFData> dof_data;
   Vector Q(wave.Height()), Q_new(wave.Height());
   Q = 0.0;
   const int ndof_total = wave.GetFESpace().GetNDofs();

   // TPV102 (rate-and-state) — total-Q (prestress in bulk Q).
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, coords);
      ZeroDOFDataPreStressTotal(dof_data, n_fault);
      for (int i = 0; i < n_fault; i++)
      { dof_data[i].tau2_nuc = TPV102Params::nuc_dtau; }
   }
   InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                        TPV102Params::tau_ini);
   wave.SetFaultDOFData(&dof_data, nqp);

   // Perturb psi on rank 0 ONLY so the shared QP psi differs across the seam.
   if (rank == 0)
   {
      for (auto &d : dof_data) { d.psi += perturb; }
   }

   // (1) Pre-step sanity: the two ranks' psi multisets DIFFER (perturbation took).
   {
      const std::vector<double> pre = SortedLocalPsi(dof_data);
      const bool identical = RanksPsiBitIdentical(pre, rank);
      if (rank == 0)
      {
         TEST_TRUE(!identical,
                   "PRE-step: rank-0 (perturbed) vs rank-1 psi DIFFER "
                   "(perturbation took; assertion is non-vacuous)");
      }
   }

   // (2) One ADER-2 macro step — the shared-face reconcile broadcasts psi.
   wave.AdvanceADER(Q, kDt, /*ader_order=*/2, Q_new);
   Q.Swap(Q_new);

   // (3) Post-step: the shared-QP psi must be BIT-IDENTICAL across the ranks.
   {
      const std::vector<double> post = SortedLocalPsi(dof_data);
      const bool identical = RanksPsiBitIdentical(post, rank);
      if (rank == 0)
      {
         TEST_TRUE(identical,
                   "POST-step: shared-fault RS psi is BIT-IDENTICAL across "
                   "ranks after the reconcile (R-020 part a)");
      }
   }

   // Cross-check with the production verifier (psi is field index 5 in its
   // payload, wave_operator.inl:6102): worst_rel == 0 confirms ALL fields,
   // including psi, are cross-rank consistent.
   {
      double worst_rel = 0.0; int worst_field = -1;
      wave.VerifySharedFaultDOFDataConsistency(1.0e-10, &worst_rel, &worst_field,
                                               /*abort_on_fail=*/false);
      double g_wr = worst_rel;
      MPI_Reduce(&worst_rel, &g_wr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (rank == 0)
      {
         TEST_TRUE(g_wr == 0.0,
                   "POST-step: VerifySharedFaultDOFDataConsistency worst_rel==0 "
                   "(all fields incl. psi consistent)");
      }
   }

   int exit_code = 0;
   if (rank == 0)
   {
      std::cout << "\n========================================\n"
                << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n"
                << "========================================\n";
      exit_code = (num_failed == 0) ? 0 : 1;
   }
   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
#else
   (void)argc; (void)argv;
   std::cout << "SKIPPED: requires MFEM with MPI\n";
   return 77;
#endif
}
