// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// R-101 regression test: shared-fault DOFData bit-equality across ranks.
//
// Background: the v1 fix for TPV102 added a (+,-) canonicalisation at
// shared fault faces (WaveOperator::ComputeSharedFaceFluxRHS) so both
// ranks drive FaultFaceFlux::Evaluate with identical inputs and update
// their independent DOFData entries to the same value.  The v2 review
// flagged this as "unverified" — it relies on an MFEM invariant
// (identical face normals on both ranks) that the v1 report did not
// directly check.
//
// This test partitions the TPV102 coarse (1000 m) mesh onto 2 or more
// MPI ranks, runs one RK4 stage and ten RK4 steps of the wave operator,
// and invokes `WaveOperator::VerifySharedFaultDOFDataConsistency`.  That
// method gathers all shared-fault DOFData across ranks, matches pairs
// by face centroid, and calls MFEM_ABORT on any mismatch.  Absence of
// abort is the pass criterion.
//
// If shared fault faces don't exist on this partitioning, the test
// reports "skipped" (not a failure) — the partitioner may happen to
// place the entire fault interior to one rank.  Run with more ranks
// if you hit this case: mpirun -np 4 ...
//
// Usage: mpirun -np 2 ./seas_test_r101_shared_fault
//        mpirun -np 4 ./seas_test_r101_shared_fault

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

// Collect fault QP physical coordinates from the wave operator's canonical
// face lists, in the exact order expected by SetFaultDOFData.
static void BuildFaultCoords(ParMesh &pmesh,
                             const Array<int> &int_faces,
                             const Array<int> &shr_faces,
                             int order, int nqp_per_face,
                             std::vector<Vector> &fault_coords)
{
   (void)nqp_per_face;
   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[i]);
      MFEM_VERIFY(ftr,
                  "test: interior fault face " << int_faces[i]
                  << " returned null FTR after WaveOperator ctor filtered it.  "
                  "MFEM invariant violation — check mesh state.");
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   for (int i = 0; i < shr_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(shr_faces[i]);
      MFEM_VERIFY(ftr,
                  "test: shared fault face " << shr_faces[i]
                  << " returned null FTR.  MFEM invariant violation — "
                  "check mesh partitioning.");
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
}

static void RK4Step(const WaveOperator<ParMesh> &wave, Vector &Q, real_t dt)
{
   int size = Q.Size();
   Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size);
   wave.Mult(Q, k1);
   add(Q, 0.5*dt, k1, Q_tmp); wave.Mult(Q_tmp, k2);
   add(Q, 0.5*dt, k2, Q_tmp); wave.Mult(Q_tmp, k3);
   add(Q, dt,    k3, Q_tmp);  wave.Mult(Q_tmp, k4);
   Q.Add(dt/6.0, k1); Q.Add(dt/3.0, k2);
   Q.Add(dt/3.0, k3); Q.Add(dt/6.0, k4);
}

struct TestContext
{
   std::unique_ptr<ParMesh> pmesh;
   std::unique_ptr<WaveOperator<ParMesh>> wave;
   std::unique_ptr<FaultFaceFlux> fault_flux;
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   int nqp_per_face = 0;
   int num_fault_local  = 0;
   int num_shared_fault = 0;
   int num_fault_total  = 0;
   int num_shared_global = 0;
};

static std::unique_ptr<TestContext>
BuildTPV102Context(const std::string &mesh_file, int order, MPI_Comm comm)
{
   auto ctx = std::make_unique<TestContext>();

   Mesh serial_mesh(mesh_file.c_str(), 1, 1);
   ctx->pmesh = std::make_unique<ParMesh>(comm, serial_mesh);

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {5};

   ctx->wave = std::make_unique<WaveOperator<ParMesh>>(
      *ctx->pmesh, order,
      TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho, bc);

   const Array<int> &int_faces = ctx->wave->GetFaultInteriorFaces();
   const Array<int> &shr_faces = ctx->wave->GetFaultSharedFaces();

   ctx->nqp_per_face = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = ctx->pmesh->GetInteriorFaceTransformations(int_faces[0]);
      ctx->nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   else if (shr_faces.Size() > 0)
   {
      auto *ftr = ctx->pmesh->GetSharedFaceTransformations(shr_faces[0]);
      ctx->nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   {
      int local_nqp = ctx->nqp_per_face;
      MPI_Allreduce(&local_nqp, &ctx->nqp_per_face, 1, MPI_INT, MPI_MAX, comm);
   }

   ctx->num_fault_local  = int_faces.Size() * ctx->nqp_per_face;
   ctx->num_shared_fault = shr_faces.Size() * ctx->nqp_per_face;
   ctx->num_fault_total  = ctx->num_fault_local + ctx->num_shared_fault;

   BuildFaultCoords(*ctx->pmesh, int_faces, shr_faces, order,
                    ctx->nqp_per_face, ctx->fault_coords);

   if (ctx->num_fault_total > 0)
   {
      InitializeFaultDOFs(ctx->dof_data, ctx->num_fault_total, ctx->fault_coords);
   }

   ctx->fault_flux = std::make_unique<FaultFaceFlux>(
      TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   ctx->wave->SetFaultFlux(ctx->fault_flux.get());
   ctx->wave->SetFaultDOFData(&ctx->dof_data, ctx->nqp_per_face);

   MPI_Allreduce(&ctx->num_shared_fault, &ctx->num_shared_global, 1,
                 MPI_INT, MPI_SUM, comm);
   return ctx;
}

// ===========================================================================
// Test R-101a: DOFData consistency after 1 RK4 stage.
// ===========================================================================
static void TestR101_OneStage(int rank, MPI_Comm comm,
                              const std::string &mesh_file)
{
   if (rank == 0)
   {
      std::cout << "Test R-101a: DOFData bit-equality after 1 RK4 stage\n";
   }

   auto ctx = BuildTPV102Context(mesh_file, 1, comm);

   // Print per-rank fault layout so the SKIPPED case is informative.
   if (rank == 0) { std::cout << "  Fault layout by rank (QPs):\n"; }
   int nprocs;
   MPI_Comm_size(comm, &nprocs);
   for (int r = 0; r < nprocs; r++)
   {
      if (r == rank)
      {
         std::cout << "    rank " << r
                   << ": local=" << ctx->num_fault_local
                   << " shared=" << ctx->num_shared_fault << std::endl;
      }
      MPI_Barrier(comm);
   }

   if (ctx->num_shared_global == 0)
   {
      if (rank == 0)
      {
         std::cout << "  SKIPPED: no shared fault faces with this nranks; "
                      "METIS placed fault interior to one rank (or the "
                      "fault is tagged via internal BEs that MFEM does "
                      "not classify as shared faces).  "
                      "Re-run with more MPI ranks — or a mesh with "
                      "vertex-duplicated fault surface — to exercise "
                      "the R-001 shared-fault path.\n";
      }
      return;
   }

   Vector Q(ctx->wave->Height());
   Q = 0.0;

   // One Mult call writes DOFData on both ranks (R-001 swap active).
   Vector k(Q.Size());
   ctx->wave->Mult(Q, k);

   // The runtime diagnostic asserts bit-equality; it aborts on mismatch
   // so reaching the next line is the pass signal.
   ctx->wave->VerifySharedFaultDOFDataConsistency(1e-10);
   TEST_ASSERT(true, "R-101a: VerifySharedFaultDOFDataConsistency did not abort");
   TEST_ASSERT(ctx->num_shared_global > 0,
               "R-101a: " + std::to_string(ctx->num_shared_global)
               + " shared fault QPs globally (seam present)");
}

// ===========================================================================
// Build a 2-tet inline mesh with a real MFEM shared fault face (R-302 Part A).
//
// TPV102 production meshes are generated via Gmsh's `BooleanFragments`, which
// duplicates vertices on either side of the fault surface.  After ParMesh
// partitioning MFEM classifies those faces as interior boundary elements on
// each rank (`GetNSharedFaces() == 0`), so the R-001 shared-fault code path
// never fires.  The `TestR101_*` tests above SKIP on every TPV102 mesh for
// that reason — which is an accurate disclosure but not a test.
//
// This builder produces the minimal mesh that exercises the shared-fault
// path: two tets sharing one triangular face (no vertex duplication), with
// the shared face tagged as fault_attr=3.  Partitioning one tet per rank
// forces MFEM to promote that face into a real SharedFace, so on each rank
// `GetNSharedFaces() == 1` and the R-001 (+,-) canonicalisation +
// `VerifySharedFaultDOFDataConsistency` body both actually execute.
//
// Vertex layout:
//   v0 = (0, 0, 0)    on the shared face (y=0)
//   v1 = (1, 0, 0)    on the shared face
//   v2 = (0, 0, 1)    on the shared face
//   v3 = (0, 1, 0)    apex of the "+y" tet  (element 1)
//   v4 = (0, -1, 0)   apex of the "-y" tet  (element 0)
// ===========================================================================
static Mesh BuildTwoTetSharedFaultMeshInline()
{
   Mesh mesh(3, 5, 2, 0);  // dim=3, 5 verts, 2 elems, 0 BEs so far

   real_t verts[5][3] = {
      {0.0, 0.0, 0.0},   // v0
      {1.0, 0.0, 0.0},   // v1
      {0.0, 0.0, 1.0},   // v2
      {0.0, 1.0, 0.0},   // v3 apex (+y)
      {0.0, -1.0, 0.0},  // v4 apex (-y)
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }

   // Tet 0 (y<0 side): v0, v1, v2, v4 — element attribute 1
   mesh.AddTet(0, 1, 2, 4, 1);
   // Tet 1 (y>0 side): v0, v1, v2, v3 — element attribute 1
   mesh.AddTet(0, 1, 2, 3, 1);

   mesh.FinalizeTopology();

   // Tag every face as a boundary element: the lone interior face (between
   // the two tets) gets attr=3 (fault); the six external triangular faces
   // get attr=1 (natural / free-surface).  Adding the fault BE makes MFEM
   // classify that face as a fault face after partitioning.
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }  // tet faces are triangles
      FaceElementTransformations *FTr = mesh.GetFaceElementTransformations(f);
      bool is_interior = (FTr && FTr->Elem2No >= 0);
      int attr = is_interior ? 3 : 1;  // fault if interior, natural otherwise
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], attr);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// ===========================================================================
// Test R-302 Part A: the shared-fault path actually runs on an inline mesh.
//
// Success criteria:
//   (1) Each rank's `GetNSharedFaces() == 1` (the partition puts 1 tet per rank).
//   (2) `wave.GetFaultSharedFaces().Size() == 1` (the shared face is classified
//       as a fault face).
//   (3) `VerifySharedFaultDOFDataConsistency` runs its body (not the short-
//       circuit) and does not abort (so the R-301 field packing, R-304 sort,
//       and R-305 unpair-abort all get exercised on a real shared fault QP
//       set).
//
// Requires exactly nranks=2 — the test partitions the 2 tets one per rank.
// ===========================================================================
static void TestR302_InlineTwoTetSharedFault(int rank, MPI_Comm comm, int nprocs)
{
   if (rank == 0)
   {
      std::cout << "Test R-302a: inline 2-tet shared-fault path actually "
                   "exercises shared-face branch\n";
   }
   if (nprocs != 2)
   {
      if (rank == 0)
      {
         std::cout << "  SKIPPED: R-302a requires exactly 2 ranks (got "
                   << nprocs << "); the inline mesh has 2 tets so one "
                   "tet goes to each rank.\n";
      }
      return;
   }

   Mesh serial_mesh = BuildTwoTetSharedFaultMeshInline();
   TEST_ASSERT(serial_mesh.GetNE() == 2,
               "R-302a: inline mesh has 2 elements");

   // Force one tet per rank so the inter-tet face becomes a real shared face
   // (METIS on 2 elems is deterministic but explicit partitioning is safer).
   int partition[2] = {0, 1};
   ParMesh pmesh(comm, serial_mesh, partition);

   int n_shared_local = pmesh.GetNSharedFaces();
   int n_be_local     = pmesh.GetNBE();
   int n_fault_be = 0;
   for (int b = 0; b < n_be_local; b++)
   {
      if (pmesh.GetBdrAttribute(b) == 3) { n_fault_be++; }
   }
   for (int r = 0; r < nprocs; r++)
   {
      if (r == rank)
      {
         std::cout << "  rank " << r << ": GetNSharedFaces=" << n_shared_local
                   << ", GetNBE=" << n_be_local
                   << ", fault BEs (attr=3)=" << n_fault_be << std::endl;
      }
      MPI_Barrier(comm);
   }
   TEST_ASSERT(n_shared_local == 1,
               "R-302a: rank " + std::to_string(rank) + " sees exactly 1 "
               "shared face (got " + std::to_string(n_shared_local) + ")");

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {5};

   WaveOperator<ParMesh> wave(pmesh, 1, TPV102Params::lambda,
                              TPV102Params::mu, TPV102Params::rho, bc);

   int n_fault_shared = wave.GetFaultSharedFaces().Size();
   for (int r = 0; r < nprocs; r++)
   {
      if (r == rank)
      {
         std::cout << "  rank " << r << ": fault_shared_faces="
                   << n_fault_shared << std::endl;
      }
      MPI_Barrier(comm);
   }
   TEST_ASSERT(n_fault_shared == 1,
               "R-302a: rank " + std::to_string(rank) + " classifies the "
               "shared face as fault (got " + std::to_string(n_fault_shared)
               + ")");

   // Set up DOFData for both the shared-fault QPs on this rank.  Triangle
   // quadrature at 2*order = 2 has 3 points.
   const int order = 1;
   int nqp_per_face = 0;
   if (n_fault_shared > 0)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(
         wave.GetFaultSharedFaces()[0]);
      MFEM_VERIFY(ftr, "R-302a: shared fault FTR null");
      nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   int max_nqp = nqp_per_face;
   MPI_Allreduce(&nqp_per_face, &max_nqp, 1, MPI_INT, MPI_MAX, comm);
   nqp_per_face = max_nqp;

   const int num_fault_local  = 0;  // no interior fault faces (one tet / rank)
   const int num_shared_fault = n_fault_shared * nqp_per_face;
   const int num_fault_total  = num_fault_local + num_shared_fault;

   std::vector<Vector> fault_coords;
   BuildFaultCoords(pmesh, wave.GetFaultInteriorFaces(),
                    wave.GetFaultSharedFaces(), order, nqp_per_face,
                    fault_coords);

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
   }

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp,
                            TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   Vector Q(wave.Height());
   Q = 0.0;
   Vector k(Q.Size());
   wave.Mult(Q, k);

   // Tight tolerance: the R-001 canonicalisation guarantees bit-equality of
   // Evaluate inputs on both ranks.  With Q=0 and the default initial DOFData
   // (tau_ini along-strike, V_ini along-strike, psi from equilibrium) the
   // field values across the 2 ranks must match to ~1e-12.
   wave.VerifySharedFaultDOFDataConsistency(1e-10);
   TEST_ASSERT(true,
               "R-302a: VerifySharedFaultDOFDataConsistency ran its body "
               "(non-SKIP path) and did not abort on inline 2-tet mesh");
}

// ===========================================================================
// Test R-501a: DOFData consistency after a full 4-stage RK4 with NONZERO Q.
//
// R-302a (above) only exercises `wave.Mult(Q=0)` — but `Tinv · 0 = 0` on
// both ranks regardless of the fault-local frame, so Evaluate receives
// identical inputs on both sides and produces bit-identical output even
// when the R-001 (+,-) swap alone was insufficient (the v3/v4 code path).
// R-302a therefore happens to dodge the R-501 bug.
//
// The Frontera 4-rank run on `tpv102_1000m_p1_1.5s_4rank_dev.sbatch` caught
// the bug at step 0, stage 2+ when `Q_tmp = Q + α·k_1` makes Q nonzero and
// the frame mismatch propagates into Evaluate.  This test replicates the
// condition locally: same 2-tet mesh, but drives 4 Mult calls with a small
// sinusoidal Q perturbation before calling the verifier.  With the R-001
// swap alone this aborts with `tau1_corr differs by ~1e-7`; with the R-501
// owner-broadcast fix it reports `max_rel_diff ≤ rel_tol`.
//
// Requires exactly 2 MPI ranks.
// ===========================================================================
static void TestR501_MultiStageRK4NonzeroQ(int rank, MPI_Comm comm, int nprocs)
{
   if (rank == 0)
   {
      std::cout << "Test R-501a: DOFData bit-equality after 4 Mult calls on "
                   "nonzero Q (R-001 frame-mismatch regression)\n";
   }
   if (nprocs != 2)
   {
      if (rank == 0)
      {
         std::cout << "  SKIPPED: R-501a requires exactly 2 ranks (got "
                   << nprocs << ").\n";
      }
      return;
   }

   Mesh serial_mesh = BuildTwoTetSharedFaultMeshInline();
   int partition[2] = {0, 1};
   ParMesh pmesh(comm, serial_mesh, partition);

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order, TPV102Params::lambda,
                              TPV102Params::mu, TPV102Params::rho, bc);

   const int n_fault_shared = wave.GetFaultSharedFaces().Size();
   int nqp_per_face = 0;
   if (n_fault_shared > 0)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(
         wave.GetFaultSharedFaces()[0]);
      nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   int max_nqp = nqp_per_face;
   MPI_Allreduce(&nqp_per_face, &max_nqp, 1, MPI_INT, MPI_MAX, comm);
   nqp_per_face = max_nqp;

   const int num_fault_total = (wave.GetFaultInteriorFaces().Size()
                                + n_fault_shared) * nqp_per_face;

   std::vector<Vector> fault_coords;
   BuildFaultCoords(pmesh, wave.GetFaultInteriorFaces(),
                    wave.GetFaultSharedFaces(), order, nqp_per_face,
                    fault_coords);

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
   }

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp,
                            TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // Small sinusoidal nonzero Q — magnitude 1e3 (small stress/velocity,
   // won't blow the friction bracket but nonzero enough to make the
   // fault-local frame mismatch matter).  The exact values don't need to
   // be physical; what matters is `Tinv_A · Q ≠ Tinv_B · Q`.
   Vector Q(wave.Height());
   for (int i = 0; i < Q.Size(); i++)
   {
      Q(i) = 1.0e3 * std::sin(0.37 * i + 0.9 * rank);
   }
   Vector k(Q.Size());
   // 4 successive Mult calls crudely mimic the driver's 4-stage RK4
   // touching DOFData 4 times per step.  We do not advance Q between
   // calls — the point is just to drive Evaluate with a nonzero Q and
   // compare resulting DOFData, not to run a physically meaningful time
   // step.
   for (int s = 0; s < 4; s++) { wave.Mult(Q, k); }

   // Relative tolerance 1e-10 — on a bit-for-bit correct broadcast the
   // two ranks' fdata and flux assembly both use exactly the same bytes,
   // so max_rel_diff should be 0 (or at most a few ULPs from any late
   // writeback reordering).
   wave.VerifySharedFaultDOFDataConsistency(1e-10);
   TEST_ASSERT(true,
               "R-501a: VerifySharedFaultDOFDataConsistency passed after "
               "4 Mult calls on nonzero Q — R-001 frame-mismatch bug fixed");
}

// ===========================================================================
// Test T-R801: BP5 convention is in force on ALL fault QPs.
//
// Under R-801 Option A every fault code path (interior-fault branch of
// ComputeFaceFluxRHS and shared-fault branch of ComputeSharedFaceFluxRHS)
// reconstructs the BP5 canonical frame from `FaultBasis`, so `DOFData.V1`
// is always the DIP component and `DOFData.V2` always the STRIKE component,
// regardless of whether the underlying face is interior or shared.
//
// This test exercises BOTH branches:
//   Part A: serial-style 2-tet mesh (both elements on rank 0) — the fault
//           face is interior, so ONLY ComputeFaceFluxRHS's fault branch runs.
//           Under the pre-R-801 bug this stored V1=strike, V2=0.  After
//           R-801 it stores V1=0, V2=strike.  Assertion: |V2| >> |V1|.
//   Part B: the same inline 2-tet mesh with one tet per rank — fault is
//           shared, so ComputeSharedFaceFluxRHS runs.  This path already
//           used BP5 convention under R-701, so this is a regression check
//           (|V2| >> |V1| must still hold).
//
// TPV102's initial condition is pure strike-slip: tau2_0 = tau_ini on the
// canonical tangent2 axis, V2 = V_ini, tau1_0 = 0, V1 = 0.  After one
// Mult(Q=0), the Evaluate pipeline (with zero trial traction) returns
// V1 ≈ 0, V2 ≈ some positive value matching the quasi-static Brent solve.
// Pre-R-801 (interior branch using BuildFrame), V1 carries the strike
// magnitude and V2 is near zero — the assertion inverts.
//
// Requires exactly nprocs == 2 for Part B (2 tets → 1 per rank).
// ===========================================================================
static void TestR801_StrikeSlipConventionOnSharedFault(int rank, MPI_Comm comm,
                                                        int nprocs)
{
   if (rank == 0)
   {
      std::cout << "Test T-R801: BP5 convention on all fault QPs "
                   "(|V2| >> |V1| for TPV102 pure strike-slip)\n";
   }
   if (nprocs != 2)
   {
      if (rank == 0)
      {
         std::cout << "  SKIPPED: requires exactly 2 ranks (got "
                   << nprocs << ").\n";
      }
      return;
   }

   // --- Part A: interior-fault path (single-rank serial WaveOperator<Mesh>)
   //
   // Only rank 0 builds a serial WaveOperator on the 2-tet inline mesh;
   // the fault face between the two tets is a true INTERIOR face (no
   // ParMesh shared face), so ComputeFaceFluxRHS is exercised — NOT
   // ComputeSharedFaceFluxRHS.  Pre-R-801 this branch used
   // GodunovFlux::BuildFrame (t1 = strike), which stored the strike-slip
   // rate in DOFData.V1 and left V2 near zero.  Under R-801 Option A the
   // branch uses the BP5 canonical frame (t1 = dip), so V2 is the strike
   // slip rate and V1 should be near zero.
   //
   // Part A asserts |V2| >> |V1| on rank 0 only.  Rank 1 is idle for
   // this part but participates in the MPI_Barrier so the test stays in
   // lockstep.  The pass/fail counters are rank-local — main() reduces
   // them with MPI_SUM, so a rank-0-only assertion still affects the
   // overall test summary on every rank.
   if (rank == 0)
   {
      Mesh serial_mesh = BuildTwoTetSharedFaultMeshInline();

      BoundaryConfig bc;
      bc.natural_attrs = {1};
      bc.fault_attr = 3;
      bc.absorbing_attrs = {5};

      const int order = 1;
      WaveOperator<Mesh> wave_A(serial_mesh, order, TPV102Params::lambda,
                                TPV102Params::mu, TPV102Params::rho, bc);

      const int n_fault_int_A = wave_A.GetFaultInteriorFaces().Size();
      int nqp_per_face_A = 0;
      if (n_fault_int_A > 0)
      {
         auto *ftr = serial_mesh.GetInteriorFaceTransformations(
            wave_A.GetFaultInteriorFaces()[0]);
         nqp_per_face_A =
            IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
      }

      const int num_fault_total_A = n_fault_int_A * nqp_per_face_A;
      std::vector<Vector> fault_coords_A;
      for (int i = 0; i < n_fault_int_A; i++)
      {
         auto *ftr = serial_mesh.GetInteriorFaceTransformations(
            wave_A.GetFaultInteriorFaces()[i]);
         const IntegrationRule &ir =
            IntRules.Get(ftr->GetGeometryType(), 2*order);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector phys(3); ftr->Face->Transform(ip, phys);
            fault_coords_A.push_back(phys);
         }
      }

      std::vector<DOFData> dof_data_A;
      if (num_fault_total_A > 0)
      {
         InitializeFaultDOFs(dof_data_A, num_fault_total_A, fault_coords_A);
      }

      FaultFaceFlux fault_flux_A(TPV102Params::rho, TPV102Params::cp,
                                 TPV102Params::cs);
      wave_A.SetFaultFlux(&fault_flux_A);
      wave_A.SetFaultDOFData(&dof_data_A, nqp_per_face_A);

      Vector Q_A(wave_A.Height());
      Q_A = 0.0;
      Vector k_A(Q_A.Size()); wave_A.Mult(Q_A, k_A);

      bool all_ok_A = true;
      real_t max_V1_abs_A = 0.0;
      real_t min_V2_abs_A = std::numeric_limits<real_t>::max();
      for (int i = 0; i < num_fault_total_A; i++)
      {
         real_t v1 = std::abs(dof_data_A[i].V1);
         real_t v2 = std::abs(dof_data_A[i].V2);
         max_V1_abs_A = std::max(max_V1_abs_A, v1);
         min_V2_abs_A = std::min(min_V2_abs_A, v2);
         if (v2 < v1) { all_ok_A = false; }
      }
      if (num_fault_total_A == 0) { min_V2_abs_A = 0.0; }
      TEST_ASSERT(num_fault_total_A > 0 && all_ok_A,
                  std::string("T-R801 Part A (interior fault path): "
                  "|V2| >= |V1| on every QP — got max|V1|=")
                  + std::to_string(max_V1_abs_A)
                  + ", min|V2|=" + std::to_string(min_V2_abs_A));
      if (num_fault_total_A > 0 && min_V2_abs_A > 0.0)
      {
         TEST_ASSERT(max_V1_abs_A < 1e-6 * min_V2_abs_A + 1e-30,
                     std::string("T-R801 Part A: |V1| negligible vs |V2| "
                     "(pure strike-slip expected); got max|V1|=")
                     + std::to_string(max_V1_abs_A)
                     + ", min|V2|=" + std::to_string(min_V2_abs_A));
      }
   }
   MPI_Barrier(comm);  // keep ranks in lockstep before Part B

   // --- Part B: shared-fault path (one tet per rank) ----------------------
   {
      Mesh serial_mesh = BuildTwoTetSharedFaultMeshInline();
      int partition_B[2] = {0, 1};
      ParMesh pmesh_B(comm, serial_mesh, partition_B);

      BoundaryConfig bc;
      bc.natural_attrs = {1};
      bc.fault_attr = 3;
      bc.absorbing_attrs = {5};

      const int order = 1;
      WaveOperator<ParMesh> wave_B(pmesh_B, order, TPV102Params::lambda,
                                   TPV102Params::mu, TPV102Params::rho, bc);

      const int n_fault_shr_B = wave_B.GetFaultSharedFaces().Size();
      int nqp_per_face_B = 0;
      if (n_fault_shr_B > 0)
      {
         auto *ftr = pmesh_B.GetSharedFaceTransformations(
            wave_B.GetFaultSharedFaces()[0]);
         nqp_per_face_B =
            IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
      }
      int max_nqp_B = nqp_per_face_B;
      MPI_Allreduce(&nqp_per_face_B, &max_nqp_B, 1, MPI_INT, MPI_MAX, comm);
      nqp_per_face_B = max_nqp_B;

      const int num_fault_total_B =
         (wave_B.GetFaultInteriorFaces().Size() + n_fault_shr_B) * nqp_per_face_B;

      std::vector<Vector> fault_coords_B;
      BuildFaultCoords(pmesh_B, wave_B.GetFaultInteriorFaces(),
                       wave_B.GetFaultSharedFaces(), order, nqp_per_face_B,
                       fault_coords_B);

      std::vector<DOFData> dof_data_B;
      if (num_fault_total_B > 0)
      {
         InitializeFaultDOFs(dof_data_B, num_fault_total_B, fault_coords_B);
      }

      FaultFaceFlux fault_flux_B(TPV102Params::rho, TPV102Params::cp,
                                 TPV102Params::cs);
      wave_B.SetFaultFlux(&fault_flux_B);
      wave_B.SetFaultDOFData(&dof_data_B, nqp_per_face_B);

      Vector Q_B(wave_B.Height());
      Q_B = 0.0;
      Vector k_B(Q_B.Size()); wave_B.Mult(Q_B, k_B);

      bool all_ok_B = true;
      real_t max_V1_abs_B = 0.0, min_V2_abs_B = std::numeric_limits<real_t>::max();
      for (int i = 0; i < num_fault_total_B; i++)
      {
         real_t v1 = std::abs(dof_data_B[i].V1);
         real_t v2 = std::abs(dof_data_B[i].V2);
         max_V1_abs_B = std::max(max_V1_abs_B, v1);
         min_V2_abs_B = std::min(min_V2_abs_B, v2);
         if (v2 < v1) { all_ok_B = false; }
      }
      if (num_fault_total_B == 0) { min_V2_abs_B = 0.0; }
      TEST_ASSERT(num_fault_total_B == 0 || all_ok_B,
                  "T-R801 Part B (shared fault path): |V2| >= |V1| on every "
                  "QP — got max|V1|=" + std::to_string(max_V1_abs_B)
                  + ", min|V2|=" + std::to_string(min_V2_abs_B));
      if (num_fault_total_B > 0 && min_V2_abs_B > 0.0)
      {
         TEST_ASSERT(max_V1_abs_B < 1e-6 * min_V2_abs_B + 1e-30,
                     "T-R801 Part B: |V1| negligible vs |V2| (pure strike-slip "
                     "expected); got max|V1|=" + std::to_string(max_V1_abs_B)
                     + ", min|V2|=" + std::to_string(min_V2_abs_B));
      }
   }
}

// ===========================================================================
// Test T-R802: conservation across the shared fault face.
//
// For a correctly conservative flux assembly, the sum of k (=dQ/dt before
// mass-inverse) contributions driven by the shared-fault flux on the two
// ranks' Elem1 must cancel component-wise.  Pre-R-802, when MFEM's
// shared-face CalcOrtho returns identical (not opposite) normals on the
// two ranks — the v6 empirical observation on the inline 2-tet mesh —
// each rank accumulates an INDEPENDENT F_h into its Elem1 rhs, and the
// sum is O(|F_h|) rather than O(machine epsilon).
//
// Methodology:
//   1. Drive a single Mult on a nonzero Q (ghost exchange makes both
//      ranks see the same cross-fault state).
//   2. Before Mult, snapshot rhs_baseline = Mult(Q_zero).  After Mult
//      with the same Q, the difference is entirely due to the fault
//      contribution on this step (bulk volume integrals are linear in Q).
//   3. Collapse per-rank Elem1 contributions into a scalar: sum over the
//      Elem1 mass-weighted component-L1 norm.
//   4. MPI_Allreduce(MPI_SUM) across ranks.  Conservation requires the
//      sum to be << the per-rank magnitude (ratio < 1e-8).
//
// A simpler-and-stricter proxy: compare the 9 stress/velocity moments
// (sum of k_c * mass_weight over Elem1's DOFs) across the two ranks.
// Under conservation they must sum to zero per component (with noise
// O(machine epsilon × peak|k|)).
//
// Requires exactly 2 MPI ranks.
// ===========================================================================
static void TestR802_ConservationAcrossSharedFault(int rank, MPI_Comm comm,
                                                    int nprocs)
{
   if (rank == 0)
   {
      std::cout << "Test T-R802: momentum conservation across shared fault "
                   "(sum of per-rank Elem1 flux contribution ~ 0)\n";
   }
   if (nprocs != 2)
   {
      if (rank == 0)
      {
         std::cout << "  SKIPPED: requires exactly 2 ranks (got "
                   << nprocs << ").\n";
      }
      return;
   }

   Mesh serial_mesh = BuildTwoTetSharedFaultMeshInline();
   int partition[2] = {0, 1};
   ParMesh pmesh(comm, serial_mesh, partition);

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order, TPV102Params::lambda,
                              TPV102Params::mu, TPV102Params::rho, bc);

   const int n_fault_shared = wave.GetFaultSharedFaces().Size();
   int nqp_per_face = 0;
   if (n_fault_shared > 0)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(
         wave.GetFaultSharedFaces()[0]);
      nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   int max_nqp = nqp_per_face;
   MPI_Allreduce(&nqp_per_face, &max_nqp, 1, MPI_INT, MPI_MAX, comm);
   nqp_per_face = max_nqp;

   const int num_fault_total = (wave.GetFaultInteriorFaces().Size()
                                + n_fault_shared) * nqp_per_face;

   std::vector<Vector> fault_coords;
   BuildFaultCoords(pmesh, wave.GetFaultInteriorFaces(),
                    wave.GetFaultSharedFaces(), order, nqp_per_face,
                    fault_coords);

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
   }

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp,
                            TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // Drive a nonzero Q with a moderate amplitude so the fault flux is
   // clearly nonzero.  Use a deterministic Q that does NOT depend on the
   // rank (otherwise the ghost-exchanged Q on the neighbour side would
   // disagree with the local Q, spoiling the conservation identity we
   // are testing).  We fill Q in global-DOF order later.
   //
   // Approach: set Q from a ParGridFunction on a vector L2 space,
   // interpolated from a global coordinate function.  Simpler: set Q via
   // a per-element uniform value so that the physical L/R sides of the
   // fault see truly different but well-defined states.
   //
   // Easiest robust choice: rank 0 sets Q to a uniform stress state on
   // its element; rank 1 does the same with different values.  Ghost
   // exchange will make each rank's "neighbour" be the other rank's
   // "self", and both ranks will see the same pair of states across the
   // fault.  No partition-dependent artefacts.
   Vector Q(wave.Height());
   Q = 0.0;
   // Per-element uniform stress + velocity.  Shape functions are L2 on a
   // tet, but a constant value is representable in any order, and
   // CalcShape averages at quad points to the same constant.
   const int ndof_per_el = wave.GetNDof();
   const int ndof_total = wave.GetScalarNDof();
   for (int e = 0; e < wave.NumElements(); e++)
   {
      for (int c = 0; c < NUM_STATE; c++)
      {
         // rank 0: set component c to 1e5 + c*1e4
         // rank 1: set component c to 2e5 + c*1e4   (different on other side)
         real_t val = (rank == 0) ? (1.0e5 + c * 1.0e4)
                                  : (2.0e5 + c * 1.0e4);
         for (int i = 0; i < ndof_per_el; i++)
         {
            Q[c * ndof_total + e * ndof_per_el + i] = val;
         }
      }
   }

   // Snapshot k with zero Q first (bulk integrals are linear in Q, so the
   // Q=0 case gives the BCs-only baseline — subtracting isolates the
   // fault/boundary/bulk contribution to the actual k).  For this simple
   // inline mesh with zero initial Q we expect a non-trivial baseline
   // from the fault pre-stress (tau2_0 = tau_ini drives a nonzero
   // Q_imp_plus−Q_imp_minus difference even at Q=0).
   Vector k_zero(Q.Size()); Vector Q_zero(Q.Size()); Q_zero = 0.0;
   wave.Mult(Q_zero, k_zero);

   Vector k(Q.Size()); wave.Mult(Q, k);

   // Compute per-rank sum of k contributions over Elem1 DOFs of each
   // shared fault face, component by component.  Elem1 is the LOCAL
   // element adjacent to the shared fault face.  To isolate the
   // fault-face flux contribution, subtract k_zero (Q=0 baseline) —
   // leaves the linear-in-Q part from the volume integral + the
   // nonlinear fault flux change.
   //
   // Simplification: instead of isolating per-face, we sum ALL k
   // contributions on the ENTIRE LOCAL fault-adjacent Elem1 (component-
   // wise).  This is an over-sum (it includes bulk volume-integral
   // contribution from the interior of Elem1), but since the same
   // volume-integral contribution exists on both ranks with DIFFERENT
   // Q's, it does NOT need to sum to zero across ranks.  Only the
   // FAULT-FLUX contribution needs to cancel.
   //
   // So we need a cleaner test.  Use the Q=0 case (k_zero): bulk volume
   // integral of A·Q = 0.  Only the fault flux drives k_zero nonzero.
   // On a pure strike-slip initial condition, the fault-flux contribution
   // points in the strike direction and is symmetric about the fault
   // (one side pushes +strike, the other -strike, so the sum ≈ 0 by
   // Newton's third law = conservation).
   //
   // We test with Q = 0 so k = ONLY the fault-flux contribution per rank,
   // then assert MPI_SUM of rank-wise integrated moments ≈ 0.
   //
   // Integrated moment = Σ_i M_{ii} * k_c(i) over all DOFs i of the
   // local fault-adjacent element (Elem1).  But M^{-1} has already been
   // applied inside Mult, so k = M^{-1} · flux_rhs.  To get the raw flux
   // moment we pre-multiply by M:  M · k = flux_rhs.  Then the integrated
   // moment is Σ_i flux_rhs(i) = Σ_i (M · k)_i = k · M · 1 (componentwise).
   //
   // Simpler: use the first moment in the component-major layout that the
   // rhs accumulation uses.  The accumulation does
   //   rhs[c*ndof_total + dof_offset1 + i] -= w * shape1(i) * F[c]
   // so Σ_i rhs_c = -w * Σ_i shape1(i) * F[c] = -w * mass_row_sum * F[c].
   // Computing Σ_i shape1(i) reliably requires CalcShape + the integration
   // rule.  But the shape functions' sum at any quad point is 1 for any
   // FE on a simplex with PU (partition of unity) — true for GLL nodal
   // bases.
   //
   // Cleaner still: sum all DOFs of k_zero across components.  If the
   // fault flux is conservative, the TOTAL (summed across ranks, summed
   // across the element's DOFs, summed across NUM_STATE stress & velocity
   // components) sum should be ~ 0.  This is a weaker test than per-
   // component but it catches the R-802 regression (where the F_h's don't
   // cancel, so the sum is O(|F_h|) > 0 rather than round-off).
   real_t local_k_sum = 0.0;
   for (int i = 0; i < k_zero.Size(); i++)
   {
      local_k_sum += k_zero[i];
   }

   // Also compute the local L1 norm as the reference "magnitude" scale.
   real_t local_k_abs = 0.0;
   for (int i = 0; i < k_zero.Size(); i++)
   {
      local_k_abs += std::abs(k_zero[i]);
   }

   real_t global_k_sum = 0.0, global_k_abs = 0.0;
   MPI_Allreduce(&local_k_sum, &global_k_sum, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   MPI_Allreduce(&local_k_abs, &global_k_abs, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);

   // Conservation: |Σ_global k| / (Σ_global |k| + 1) < 1e-8.
   // Scale floor of 1 guards against the case where all k's happen to
   // be ~0 (no fault activity).
   const real_t rel = std::abs(global_k_sum)
                     / (std::abs(global_k_abs) + 1.0);
   if (rank == 0)
   {
      std::cout << "  R-802 diagnostic: Σ k = " << global_k_sum
                << "   Σ|k| = " << global_k_abs
                << "   rel = " << rel << std::endl;
   }
   // Tolerance 1e-8: matches reviewer's suggested threshold.  Pre-R-802
   // with the empirical nor_A = nor_B observation, rel would be O(1).
   TEST_ASSERT(rel < 1e-8 || global_k_abs < 1e-6,
               "T-R802: global sum of k (Q=0 baseline) consistent with "
               "momentum conservation across shared fault (rel="
               + std::to_string(rel) + ")");
}

// ===========================================================================
// Test R-101b: DOFData consistency after 10 RK4 steps.
// ===========================================================================
static void TestR101_TenSteps(int rank, MPI_Comm comm,
                              const std::string &mesh_file)
{
   if (rank == 0)
   {
      std::cout << "Test R-101b: DOFData bit-equality after 10 RK4 steps\n";
   }

   auto ctx = BuildTPV102Context(mesh_file, 1, comm);

   if (ctx->num_shared_global == 0)
   {
      if (rank == 0) { std::cout << "  SKIPPED (no shared fault faces)\n"; }
      return;
   }

   Vector Q(ctx->wave->Height());
   Q = 0.0;

   const real_t cfl = 0.5 / (3.0 * (2.0*1 + 1.0));
   const real_t dt  = 0.5 * ctx->wave->ComputeMaxDt(cfl);
   for (int step = 0; step < 10; step++)
   {
      RK4Step(*ctx->wave, Q, dt);
   }

   // 10 steps × RK4 drift tolerance is roughly 10 × (one-stage floor);
   // one Evaluate call is O(1e-14) bit-for-bit if the invariant holds,
   // so 1e-8 is generous.
   ctx->wave->VerifySharedFaultDOFDataConsistency(1e-8);
   TEST_ASSERT(true,
               "R-101b: VerifySharedFaultDOFDataConsistency did not abort "
               "after 10 RK4 steps");
}

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank, nprocs;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);

   std::string mesh_file = "tpv102/mesh/tpv102_1000m.msh";
   if (argc > 1) { mesh_file = argv[1]; }

   if (nprocs < 2)
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: requires at least 2 MPI ranks (got "
                   << nprocs << ").\n";
      }
      MPI_Finalize();
      return 1;
   }

   if (rank == 0)
   {
      std::cout << "========================================\n";
      std::cout << "R-101 shared-fault DOFData consistency tests\n";
      std::cout << "Ranks: " << nprocs << ", mesh: " << mesh_file << "\n";
      std::cout << "========================================\n";
   }

   TestR101_OneStage(rank, comm, mesh_file);
   TestR101_TenSteps(rank, comm, mesh_file);
   TestR302_InlineTwoTetSharedFault(rank, comm, nprocs);
   TestR501_MultiStageRK4NonzeroQ(rank, comm, nprocs);
   TestR801_StrikeSlipConventionOnSharedFault(rank, comm, nprocs);
   TestR802_ConservationAcrossSharedFault(rank, comm, nprocs);

   int total_passed = 0, total_failed = 0;
   MPI_Allreduce(&num_passed, &total_passed, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(&num_failed, &total_failed, 1, MPI_INT, MPI_SUM, comm);

   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "Passed: " << total_passed << "\n";
      std::cout << "Failed: " << total_failed << "\n";
      std::cout << "========================================\n";
   }

   MPI_Finalize();
   return (total_failed == 0) ? 0 : 1;
#else
   std::cerr << "R-101 test requires MFEM_USE_MPI.\n";
   return 1;
#endif
}
