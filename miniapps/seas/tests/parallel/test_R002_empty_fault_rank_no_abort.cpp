// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.5.0 — REVIEW R-002 regression gate (parallel, np=4).
//
// Purpose: assert that WaveOperator::UsePrecomputedFaceFluxes(true) does
// NOT abort on MPI ranks whose local partition contains no fault-
// adjacent elements.  Pre-fix wave_operator.inl:504-510 ran an
// MFEM_VERIFY that fired whenever fault_interior_faces_ AND
// fault_shared_faces_ were both empty while bc_.fault_attr > 0 — a
// legitimate state for ranks that hold a fault-free slice of the
// domain under sparse-fault partitioning.
//
// Fixture: a 1x4x1 tet-split cartesian mesh (24 tets) with the fault
// attribute at the y=L/2 internal plane.  Explicit partitioning places
// tets in y-layer [0, L/4] on rank 2 and y-layer [3L/4, L] on rank 3,
// so both rank 2 and rank 3 hold NO fault-adjacent faces (neither
// interior nor shared).  Rank 0 owns the y-layer [L/4, L/2] and rank 1
// owns [L/2, 3L/4], so they jointly own the shared fault faces.
//
// Acceptance:
//   * On every rank, UsePrecomputedFaceFluxes(true) returns without
//     aborting.
//   * On ranks 2 and 3 the local fault arrays are in fact empty (the
//     fixture is verified so the test is not vacuous).
//   * The flag is set to true on every rank afterwards.
//
// Run: mpirun -np 4 ./seas_test_R002_empty_fault_rank_no_abort

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
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

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

namespace
{

constexpr real_t kL      = 1000.0;
constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;

// Build a 1x4x1 tet-split cartesian mesh with the plane y=L/2 tagged
// fault_attr=3 and all other exterior faces tagged natural_attr=1.
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
         // Interior face: tag with fault_attr if its centroid lies
         // on the y=L/2 internal plane.
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * kL) < 1e-8)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);
         }
         continue;
      }
      // True domain-boundary face: natural_attr.
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }

   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Assign each element to a rank by its y-centroid.  Rank 2 gets the
// bottom layer (y in [0, L/4]), rank 0 gets [L/4, L/2], rank 1 gets
// [L/2, 3L/4], rank 3 gets [3L/4, L].  Under this partition the
// fault at y=L/2 becomes a shared boundary between ranks 0 and 1,
// and ranks 2 / 3 hold no fault-adjacent faces.
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

   int ret = 0;

   if (nprocs != 4)
   {
      if (rank == 0)
      {
         std::cerr << "[R-002] TEST REQUIRES np=4, got " << nprocs
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

   WaveOperator<ParMesh> wave(pmesh, /*order=*/1, kLambda, kMu, kRho, bc);

   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   const int n_fault_int    = wave.GetFaultInteriorFaces().Size();
   const int n_fault_shared = wave.GetFaultSharedFaces().Size();

   for (int r = 0; r < nprocs; r++)
   {
      if (r == rank)
      {
         std::cout << "  rank " << r << ": fault_interior="
                   << n_fault_int << ", fault_shared="
                   << n_fault_shared << "\n";
      }
      MPI_Barrier(comm);
   }

   // Non-vacuity: ranks 2 and 3 must actually have empty local fault
   // arrays (otherwise the test would pass even without the fix).
   if (rank == 2 || rank == 3)
   {
      if (n_fault_int != 0 || n_fault_shared != 0)
      {
         std::cerr << "[R-002] FIXTURE ERROR on rank " << rank
                   << ": expected empty local fault arrays, got "
                   << "fault_interior=" << n_fault_int
                   << ", fault_shared=" << n_fault_shared << "\n";
         ret = 1;
      }
   }

   // The core gate: UsePrecomputedFaceFluxes(true) must complete on
   // EVERY rank, including the fault-free ones.  Pre-fix, ranks 2/3
   // would abort inside MFEM_VERIFY before reaching this assertion.
   try
   {
      wave.UsePrecomputedFaceFluxes(true);
   }
   catch (const std::exception &e)
   {
      std::cerr << "[R-002] rank " << rank << " threw: " << e.what() << "\n";
      ret = 1;
   }

   const bool flag_on = wave.UsingPrecomputedFaceFluxes();
   if (!flag_on)
   {
      std::cerr << "[R-002] rank " << rank
                << ": flag not set after UsePrecomputedFaceFluxes(true)\n";
      ret = 1;
   }

   // Reduce across ranks: any nonzero ret → global failure.
   int global_ret = 0;
   MPI_Allreduce(&ret, &global_ret, 1, MPI_INT, MPI_MAX, comm);

   if (rank == 0)
   {
      std::cout << "[R-002] " << (global_ret == 0 ? "PASS" : "FAIL")
                << " — UsePrecomputedFaceFluxes(true) survives empty-"
                << "fault ranks under np=4 skewed partition\n";
   }

   MPI_Finalize();
   return global_ret;
}
