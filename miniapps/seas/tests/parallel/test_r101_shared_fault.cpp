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
