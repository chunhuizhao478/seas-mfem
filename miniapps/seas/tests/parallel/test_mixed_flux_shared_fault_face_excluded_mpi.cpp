// Round-11 R-1100 — invariant `central_flux_face_set_ ∩ fault_faces == ∅`
// must hold on EVERY rank, including the non-BE-owner rank of a shared
// fault face (parallel, np=2).
//
// Bug background:
//   `BuildCentralFluxFaceSet_::is_fault_face` originally consulted only
//   `face_bdr_attr_` to classify a face as fault.  In a ParMesh with the
//   fault on the rank seam, MFEM keeps the BE on EXACTLY one of the two
//   ranks (the "BE-owner"); the other rank's `face_bdr_attr_` entry for
//   the shared fault face is 0.  The non-BE-owner's `is_fault_face`
//   returns false on the shared fault face → the local Adjacent walk
//   wrongly inserts the shared fault face into `central_flux_face_set_`.
//   The R-001 MPI key exchange then propagates the wrong key to the
//   BE-owner, so both ranks end up with the shared fault face in the
//   central set — breaking the invariant.
//
//   The bug is currently inert because the shared-face dispatch
//   (`ComputeSharedFaceFluxRHS` / `ComputeADERSharedFaceFluxRHS`)
//   classifies fault status via the merged `shared_face_bdr_attr_`
//   table and the fault branch wins.  But the data structure invariant
//   is broken — a regression hazard.
//
// Fixture: 1x2x1 hex/tet mesh with the fault at y=L/2 (between the
// two hexes).  Manual partition: hex 0 (y < L/2) -> rank 0; hex 1
// (y > L/2) -> rank 1.  The fault face IS the rank seam, so
// `fault_shared_faces_` is populated on BOTH ranks and `face_bdr_attr_`
// for the shared fault face is populated on the BE-owner only.
//
// Acceptance:
//   * Both ranks have at least one fault_shared_face (sanity).
//   * For each fault_shared_face on this rank, its mesh-face index is
//     NOT in central_flux_face_set_, in BOTH Adjacent and AllContinuous
//     modes.
//
// PRE-fix behavior: non-BE-owner inserts the shared fault face → FAIL.
// POST-fix behavior: rank-symmetric `is_fault_face` excludes it → PASS.
//
// Run: mpirun -np 2 ./seas_test_mixed_flux_shared_fault_face_excluded_mpi

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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
constexpr real_t kCp     = 6000.0;
constexpr real_t kCs     = 3464.0;
constexpr real_t kMu     = kRho * kCs * kCs;
constexpr real_t kLambda = kRho * kCp * kCp - 2.0 * kMu;

// 1x2x1 tet-split mesh (2 hexes, 12 tets).  Fault at y=L/2 internal
// plane; all external faces natural_attr=1.
Mesh BuildSeamFaultMesh()
{
   Mesh mesh = Mesh::MakeCartesian3D(1, 2, 1, Element::TETRAHEDRON,
                                     kL, kL, kL, /*sfc_ordering=*/false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
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

// Manual partition: lower-y hex (y<L/2) -> rank 0; upper-y hex (y>L/2)
// -> rank 1.  Fault at y=L/2 IS the rank seam.
std::vector<int> BuildYSplitPartition(const Mesh &mesh)
{
   const int ne = mesh.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; e++)
   {
      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      part[e] = (cy < 0.5 * kL) ? 0 : 1;
   }
   return part;
}

} // anonymous

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 0;
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);
   g_seas_my_rank = rank;

   if (nprocs != 2)
   {
      if (rank == 0)
      {
         std::cerr << "[R-1100/R-1500] requires np=2, got " << nprocs
                   << " — TEST FAILED.  This test is meaningless at "
                   "np!=2; CI must invoke via `mpirun -np 2`.  Returning "
                   "1 (FAIL) instead of 0 (PASS) so CI misconfigurations "
                   "surface (R-1500, mirrors the R-1412 fix in "
                   "dispatch_adjacent_mpi).\n";
      }
      MPI_Finalize();
      return 1;
   }

   if (rank == 0)
   {
      std::cout
         << "\n=== Round-11 R-1100: shared fault face must NOT enter "
         << "central_flux_face_set_ on either rank (np=2) ===\n";
   }

   Mesh serial = BuildSeamFaultMesh();
   std::vector<int> part = BuildYSplitPartition(serial);
   ParMesh pmesh(comm, serial, part.data());

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   WaveOperator<ParMesh> wave(pmesh, /*order=*/1, kLambda, kMu, kRho, bc);
   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   const int n_local_fault_int = wave.GetFaultInteriorFaces().Size();
   const int n_local_fault_shr = wave.GetFaultSharedFaces().Size();
   if (rank == 0)
   {
      std::cout << "  rank=0  fault_interior_faces=" << n_local_fault_int
                << "  fault_shared_faces=" << n_local_fault_shr << "\n";
   }
   MPI_Barrier(comm);
   if (rank == 1)
   {
      std::cout << "  rank=1  fault_interior_faces=" << n_local_fault_int
                << "  fault_shared_faces=" << n_local_fault_shr << "\n";
   }
   MPI_Barrier(comm);

   // Fixture-gate: BOTH ranks must have a non-empty fault_shared_faces_
   // (the fault IS the rank seam).  No interior fault faces (each rank
   // owns only one hex, no fault interior to it).
   const int my_fixture_ok =
      (n_local_fault_int == 0) && (n_local_fault_shr > 0);
   int all_fixture_ok = 0;
   MPI_Allreduce(&my_fixture_ok, &all_fixture_ok, 1, MPI_INT, MPI_LAND, comm);
   if (rank == 0)
   {
      if (all_fixture_ok)
      {
         std::cout << "  [fixture-gate] both ranks have shared fault "
                   << "faces and no interior fault faces OK\n";
      }
      else
      {
         std::cerr << "  [fixture-gate] FAILED — partition did not "
                   << "produce shared fault faces; test is vacuous.\n";
      }
   }
   if (!all_fixture_ok)
   {
      MPI_Finalize();
      return 1;
   }

   int total_failures = 0;

   for (int mode_idx = 0; mode_idx < 2; mode_idx++)
   {
      const MixedFluxMode mode = (mode_idx == 0)
         ? MixedFluxMode::Adjacent
         : MixedFluxMode::AllContinuous;
      const char *mode_name = (mode == MixedFluxMode::Adjacent)
         ? "Adjacent"
         : "AllContinuous";

      wave.SetMixedFluxMode(mode);

      const auto &S = wave.GetCentralFluxFaceSet();
      const Array<int> &fsf = wave.GetFaultSharedFaces();

      // For each fault_shared_face on this rank, get its mesh face
      // index; assert it's NOT in central_flux_face_set_.
      int my_violations = 0;
      for (int sf_i = 0; sf_i < fsf.Size(); sf_i++)
      {
         const int sf = fsf[sf_i];
         const int f  = pmesh.GetSharedFace(sf);
         if (S.count(f) > 0)
         {
            my_violations++;
            std::fprintf(stderr,
               "  [R-1100 VIOLATION] rank=%d  mode=%s  shared fault "
               "face f=%d (sf=%d) WRONGLY in central_flux_face_set_\n",
               rank, mode_name, f, sf);
         }
      }
      int total_violations = 0;
      MPI_Allreduce(&my_violations, &total_violations, 1, MPI_INT, MPI_SUM,
                    comm);
      if (rank == 0)
      {
         std::cout << "  [" << mode_name
                   << "] |central_set ∩ fault_shared_faces| (across all "
                   << "ranks) = " << total_violations;
         if (total_violations == 0)
         {
            std::cout
               << "  PASSED: no shared fault face is in central_flux_face_set_ "
               << "(R-1100 invariant holds on every rank)\n";
         }
         else
         {
            std::cout
               << "\n  FAILED: " << total_violations
               << " shared fault face(s) wrongly classified as central; "
               << "is_fault_face is rank-asymmetric.\n";
         }
      }
      if (total_violations > 0) { total_failures++; }
   }

   if (rank == 0)
   {
      std::cout << "\n========================================\n"
                << "  Results: " << (total_failures == 0 ? 2 : (2 - total_failures))
                << " passed, " << total_failures
                << " failed out of 2 tests\n"
                << "========================================\n";
   }

   MPI_Finalize();
   return (total_failures == 0) ? 0 : 1;
}
