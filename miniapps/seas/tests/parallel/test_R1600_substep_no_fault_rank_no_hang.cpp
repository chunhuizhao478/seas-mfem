// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV104 R-1600 regression — `EvaluateBulkAtFaultQPsCanonical` must NOT
// deadlock on a partition where some ranks have shared NON-fault faces
// but zero fault QPs of their own.
//
// History.  R-1003 §1 colocated the per-substep ghost exchange
// `q_gf.ExchangeFaceNbrData()` (NUM_STATE = 9 pairwise MPI calls) inside
// `WaveOperator<ParMesh>::EvaluateBulkAtFaultQPsCanonical`, but gated it
// on FAULT-shared face counts (`fault_shared_faces_.Size() > 0` and
// `n_total_qps > 0`).  At np>1 on a typical TPV104 partition every
// interior rank with shared non-fault faces — but no fault faces of its
// own — silently SKIPPED the collective.  Its fault-adjacent peers
// posted MPI_Irecv and spun in MPI_Wait at 100% CPU forever (the
// production 13.5-min hang at np=10 on the symmirror 1000m mesh).
//
// Fix (R-1600): gate the COLLECTIVE on `pmesh.GetNSharedFaces() > 0`
// (the safe predicate already used by `ComputeADERSharedFaceFluxRHS`)
// and keep the per-fault-face LOOP gated on
// `fault_shared_faces_.Size() > 0`.
//
// Fixture.  Reuses the R-002 1x4x1 cartesian skewed partition: ranks 2
// and 3 hold y-layers [0, L/4] and [3L/4, L] respectively (no fault
// faces), while ranks 0 and 1 hold the inner two y-layers and share the
// fault at y=L/2.  Crucially, rank 2 shares the y=L/4 internal face
// (NON-fault) with rank 0, and rank 3 shares the y=3L/4 internal face
// (NON-fault) with rank 1.  Pre-fix: rank 2 / rank 3 silently skip
// `EvaluateBulkAtFaultQPsCanonical` on `n_total_qps == 0`, leaving
// rank 0 / rank 1's pairwise `ExchangeFaceNbrData` posts unmatched
// → mpirun hangs forever.
//
// Pre-fix: this test hangs (CI must wrap with a wall-clock timeout).
// Post-fix: returns in milliseconds.  We measure elapsed time across
// all ranks via MPI_Allreduce(MPI_MAX) and FAIL if any rank exceeded a
// 5-second budget — defensive net for a future partial-skip regression
// that might slow but not freeze the call.
//
// Run: mpirun -np 4 ./seas_test_R1600_substep_no_fault_rank_no_hang

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

#include <sys/time.h>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

namespace
{

constexpr real_t kL      = 1000.0;
constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;

// Wall-clock seconds (process-relative).  MFEM ships its own StopWatch but
// using gettimeofday keeps this test free of any MFEM-version dependency
// for timer semantics.
double WallSeconds()
{
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return tv.tv_sec + 1e-6 * tv.tv_usec;
}

// 1x4x1 tet-split cartesian mesh with the y=L/2 internal plane tagged
// fault_attr=3 and exterior faces tagged natural_attr=1.  Identical
// fixture to the R-002 regression test — chosen because the y-layer
// partition below produces the exact "no-fault rank with shared non-
// fault face" topology that triggers R-1600.
Mesh BuildLinearFaultMesh()
{
   Mesh mesh = Mesh::MakeCartesian3D(1, 4, 1, Element::TETRAHEDRON,
                                     kL, kL, kL, /*sfc_ordering=*/false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * kL) < 1e-8)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }

   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Skewed partition: rank 2 → y < L/4, rank 0 → y < L/2, rank 1 → y < 3L/4,
// rank 3 → otherwise.  Fault is at y=L/2, shared between rank 0 and rank 1.
// Ranks 2 and 3 hold no fault faces but DO share interior non-fault faces
// (y=L/4 with rank 0; y=3L/4 with rank 1).
std::vector<int> BuildSkewedPartition(const Mesh &mesh)
{
   const int ne = mesh.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; e++)
   {
      Array<int> ev; mesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); v++)
      {
         cy += mesh.GetVertex(ev[v])[1];
      }
      cy /= ev.Size();
      const real_t q1 = 0.25 * kL;
      const real_t q2 = 0.50 * kL;
      const real_t q3 = 0.75 * kL;
      if      (cy < q1) { part[e] = 2; }
      else if (cy < q2) { part[e] = 0; }
      else if (cy < q3) { part[e] = 1; }
      else              { part[e] = 3; }
   }
   return part;
}

} // anonymous

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank, nprocs;
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);
   g_seas_my_rank = rank;

   if (nprocs != 4)
   {
      if (rank == 0)
      {
         std::cerr << "[R-1600] TEST REQUIRES np=4, got " << nprocs
                   << " — skipping.\n";
      }
      MPI_Finalize();
      return 0;
   }

   Mesh serial = BuildLinearFaultMesh();
   std::vector<int> part = BuildSkewedPartition(serial);
   ParMesh pmesh(comm, serial, part.data());

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order, kLambda, kMu, kRho, bc);

   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   const int n_fault_int    = wave.GetFaultInteriorFaces().Size();
   const int n_fault_shared = wave.GetFaultSharedFaces().Size();
   const int n_shared_total = pmesh.GetNSharedFaces();

   for (int r = 0; r < nprocs; r++)
   {
      if (r == rank)
      {
         std::cout << "  rank " << r << ": fault_interior=" << n_fault_int
                   << ", fault_shared="   << n_fault_shared
                   << ", total_shared="   << n_shared_total << "\n";
      }
      MPI_Barrier(comm);
   }

   // Fixture invariants — failure here means the partition didn't produce
   // the intended topology and the deadlock test is vacuous.
   int fixture_ok = 1;
   if (rank == 2 || rank == 3)
   {
      // No fault QPs on this rank.
      if (n_fault_int != 0 || n_fault_shared != 0)
      {
         std::cerr << "[R-1600] FIXTURE ERROR rank " << rank
                   << ": expected zero fault faces, got int="
                   << n_fault_int << ", shared=" << n_fault_shared << "\n";
         fixture_ok = 0;
      }
      // BUT we SHARE non-fault interior faces with the fault-adjacent
      // peer — this is the topology that triggers R-1600 pre-fix.
      if (n_shared_total <= 0)
      {
         std::cerr << "[R-1600] FIXTURE ERROR rank " << rank
                   << ": expected n_shared_total > 0 for the deadlock "
                   << "topology, got " << n_shared_total << "\n";
         fixture_ok = 0;
      }
   }
   if (rank == 0 || rank == 1)
   {
      if (n_fault_shared <= 0)
      {
         std::cerr << "[R-1600] FIXTURE ERROR rank " << rank
                   << ": expected n_fault_shared > 0, got "
                   << n_fault_shared << "\n";
         fixture_ok = 0;
      }
   }
   int fixture_global = 0;
   MPI_Allreduce(&fixture_ok, &fixture_global, 1, MPI_INT, MPI_MIN, comm);
   if (!fixture_global)
   {
      if (rank == 0)
      {
         std::cerr << "[R-1600] FIXTURE INVARIANTS VIOLATED — aborting\n";
      }
      MPI_Finalize();
      return 1;
   }

   // Set up minimum scaffolding for EvaluateBulkAtFaultQPsCanonical:
   // SetFaultDOFData populates fault_face_dof_offset_ and
   // shared_fault_dof_offset_, both of which the function reads.
   int nqp_per_face = 0;
   if (n_fault_int > 0)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(
         wave.GetFaultInteriorFaces()[0]);
      MFEM_VERIFY(ftr, "interior fault face FTR null");
      nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2 * order).GetNPoints();
   }
   if (nqp_per_face == 0 && n_fault_shared > 0)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(
         wave.GetFaultSharedFaces()[0]);
      MFEM_VERIFY(ftr, "shared fault face FTR null");
      nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2 * order).GetNPoints();
   }
   // Ranks 2/3 have neither: nqp_per_face stays 0 here.  Take the global
   // max so SetFaultDOFData receives a consistent value across ranks.
   {
      int local_nqp = nqp_per_face, max_nqp = 0;
      MPI_Allreduce(&local_nqp, &max_nqp, 1, MPI_INT, MPI_MAX, comm);
      nqp_per_face = max_nqp;
   }

   const int n_total_fault_qps =
      (n_fault_int + n_fault_shared) * nqp_per_face;
   std::vector<DOFData> dof_data(n_total_fault_qps);   // zero-initialised
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // Build a non-trivial Q so the per-substep ghost exchange has actual
   // data to transmit.  Zero would also reproduce the deadlock pre-fix
   // (the exchange runs unconditionally in the parallel block), but a
   // non-zero pattern gives a sharper post-fix correctness check should
   // a future iteration tighten the test.
   const int ndof_total = wave.GetScalarNDof();
   Vector Q(NUM_STATE * ndof_total);
   for (int i = 0; i < Q.Size(); i++)
   {
      Q(i) = 1e-3 * static_cast<real_t>((i + rank) % 17);
   }

   std::vector<real_t> Q_plus_flat, Q_minus_flat;

   // Synchronise once so the timer below excludes setup costs that vary
   // across ranks.
   MPI_Barrier(comm);
   const double t0 = WallSeconds();

   // The CALL UNDER TEST.  Pre-R-1600: ranks 2/3 short-circuit on
   // n_total_qps == 0, ranks 0/1 hang in the for-c loop's
   // ExchangeFaceNbrData (mpirun never returns).  Post-R-1600: every
   // rank with `pmesh.GetNSharedFaces() > 0` posts the 9 collectives in
   // matched lockstep; the call returns in milliseconds.
   wave.EvaluateBulkAtFaultQPsCanonical(Q, Q_plus_flat, Q_minus_flat);

   const double t1 = WallSeconds();
   const double elapsed = t1 - t0;

   double max_elapsed = 0.0;
   MPI_Allreduce(&elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, comm);

   // 5-second budget.  In practice post-fix is microseconds; the budget
   // is a defensive net for a future partial-skip regression that slows
   // but doesn't freeze the call.  A frozen call never reaches Allreduce
   // so the test process hangs — the CI runner's outer wall-clock kill
   // is the load-bearing detector for true freezes.
   const double kBudgetSeconds = 5.0;
   int local_pass = (elapsed < kBudgetSeconds) ? 1 : 0;
   int global_pass = 0;
   MPI_Allreduce(&local_pass, &global_pass, 1, MPI_INT, MPI_MIN, comm);

   // Output sizing contract from the function (R-1602 docstring): every
   // rank must observe `Q_*_flat.size() == NUM_STATE * GetNumTotalFaultQPs()`,
   // even on ranks 2/3 where that's zero.
   const size_t expected_words =
      static_cast<size_t>(NUM_STATE) *
      static_cast<size_t>(wave.GetNumTotalFaultQPs());
   int sizing_ok = (Q_plus_flat.size()  == expected_words &&
                    Q_minus_flat.size() == expected_words) ? 1 : 0;
   int global_sizing_ok = 0;
   MPI_Allreduce(&sizing_ok, &global_sizing_ok, 1, MPI_INT, MPI_MIN, comm);

   // R-1601 finite-output sanity: the batched vdim=NUM_STATE ghost
   // exchange replaces 9 sequential per-component exchanges; if the
   // byNODES layout assumption breaks (e.g., MFEM internally uses
   // byVDIM for FaceNbrData), the per-fault-face loop would read
   // mis-strided neighbour data and produce NaN/Inf values on the
   // shared-fault slice.  This catches that.  (Empty buffers on no-
   // fault ranks pass trivially.)
   int finite_ok = 1;
   for (real_t v : Q_plus_flat)  { if (!std::isfinite(v)) { finite_ok = 0; break; } }
   for (real_t v : Q_minus_flat) { if (!std::isfinite(v)) { finite_ok = 0; break; } }
   int global_finite_ok = 0;
   MPI_Allreduce(&finite_ok, &global_finite_ok, 1, MPI_INT, MPI_MIN, comm);

   if (rank == 0)
   {
      std::cout << "[R-1600] elapsed (max across ranks) = " << std::fixed
                << std::setprecision(6) << max_elapsed << " s; "
                << "budget = " << kBudgetSeconds << " s\n";
      if (!global_pass)
      {
         std::cerr << "[R-1600] FAIL — at least one rank exceeded the "
                   << kBudgetSeconds << "s wall-clock budget for "
                   << "EvaluateBulkAtFaultQPsCanonical (likely a "
                   << "re-introduced partial-skip on the per-substep "
                   << "collective).\n";
      }
      if (!global_sizing_ok)
      {
         std::cerr << "[R-1600] FAIL — output buffer sizes inconsistent "
                   << "with NUM_STATE * GetNumTotalFaultQPs() on at "
                   << "least one rank.\n";
      }
      if (!global_finite_ok)
      {
         std::cerr << "[R-1600/R-1601] FAIL — output values contain "
                   << "non-finite entries on at least one rank "
                   << "(likely a byNODES/byVDIM layout mismatch in the "
                   << "R-1601 batched ghost exchange).\n";
      }
      if (global_pass && global_sizing_ok && global_finite_ok)
      {
         std::cout << "[R-1600] PASS — np=4 skewed partition with "
                   << "no-fault ranks completes "
                   << "EvaluateBulkAtFaultQPsCanonical without hanging "
                   << "(R-1600 collective contract honoured) and with "
                   << "finite outputs (R-1601 batched exchange OK).\n";
      }
   }

   MPI_Finalize();
   return (global_pass && global_sizing_ok && global_finite_ok) ? 0 : 1;
}
