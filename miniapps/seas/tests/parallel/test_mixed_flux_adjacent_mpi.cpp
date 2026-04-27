// Round-11 R-001 / R-002 — MPI cross-rank consistency of Adjacent mixed-
// flux face set (parallel, np=2).
//
// Purpose: verify that `WaveOperator::BuildCentralFluxFaceSet_` in
// `MixedFluxMode::Adjacent` produces a CROSS-RANK CONSISTENT
// `central_flux_face_set_`.  For every shared non-fault face, both
// ranks must agree on membership; otherwise rank A dispatches `Central`
// and rank B dispatches `Interior` on the same physical face, breaking
// conservation across the rank seam by `0.5·|A_n|·jump` per face per
// macro step.
//
// Fixture: 1x2x2 tet-split cartesian mesh (4 hexes, 24 tets).  Fault is
// the y=L/2 internal plane (between cj=0 and cj=1 hexes).  Manual
// partition: lower-z hexes (z<L/2) -> rank 1; upper-z hexes (z>L/2) ->
// rank 0.  The fault is interior to BOTH ranks' z-slabs.  The rank seam
// at z=L/2 carries non-fault shared faces.
//
// In the 6-tet split of each fault-adjacent hex, the tet `pat[1] = {0,
// 3, 2, 7}` has BOTH a fault face (+y) AND a face on the rank seam
// (-z for upper-z hex, +z for lower-z hex).  This tet IS in E_fault_adj
// on its rank.  Walking its faces, the local Adjacent algorithm reaches
// the rank-seam face and (correctly) inserts it.
//
// The PEER tet across the rank seam is a DIFFERENT tet of the lower-z
// (or upper-z) hex — one whose own boundary faces do NOT include the
// fault.  That peer tet is NOT in E_fault_adj on its rank, so the local
// walk on that rank misses the shared face.  Rank A inserts; rank B
// does not.  R-001 trigger.
//
// PRE-fix behavior (without R-001):
//   - Rank 0 walks E_fault_adj = {tets-of-hex0, tets-of-hex1}; hits the
//     y=L/2 shared face via hex-1's faces; inserts.
//   - Rank 1's E_fault_adj is empty (no fault on rank 1).  Set is empty.
//   - Asymmetric => test FAILS.
//
// POST-fix behavior (R-001 Allgatherv exchange applied):
//   - Rank 0 contributes the y=L/2 face's global vertex key.
//   - Rank 1 sees the key in `global_keys`, walks its shared faces,
//     finds the matching one, inserts.
//   - Symmetric => test PASSES.
//
// Acceptance:
//   * For every shared face, the set of reported memberships across
//     ranks contains exactly one value (either {0} or {1}, never {0,1}).
//
// Run: mpirun -np 2 ./seas_test_mixed_flux_adjacent_mpi

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
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

// 1x2x2 tet-split mesh (4 hexes, 24 tets).  Fault on y=L/2 internal
// plane; all external faces natural_attr=1.
Mesh BuildOrthogonalFaultMesh()
{
   Mesh mesh = Mesh::MakeCartesian3D(1, 2, 2, Element::TETRAHEDRON,
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

// Partition: lower-z hexes (z < L/2) -> rank 1; upper-z hexes (z > L/2)
// -> rank 0.  Fault at y=L/2 is INTERIOR to BOTH ranks (each rank holds
// 2 hexes, 1 on each side of the fault, sharing the fault face).  The
// rank seam at z=L/2 carries non-fault shared faces, and pat[1] of
// each fault-adjacent hex has BOTH a fault face (+y) AND a face on
// the rank seam (the -z or +z face).
std::vector<int> BuildSkewedPartition(const Mesh &mesh)
{
   const int ne = mesh.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; e++)
   {
      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t cz = 0.0;
      for (int v = 0; v < ev.Size(); v++) { cz += mesh.GetVertex(ev[v])[2]; }
      cz /= ev.Size();
      part[e] = (cz < 0.5 * kL) ? 1 : 0;
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
         std::cerr << "[R-001/R-002/R-1500 mixed-flux MPI test] requires "
                   "np=2, got " << nprocs << " — TEST FAILED.  This test "
                   "is meaningless at np!=2; CI must invoke via "
                   "`mpirun -np 2`.  Returning 1 (FAIL) instead of 0 "
                   "(PASS) so CI misconfigurations surface (R-1500, "
                   "mirrors the R-1412 fix in dispatch_adjacent_mpi).\n";
      }
      MPI_Finalize();
      return 1;
   }

   if (rank == 0)
   {
      std::cout
         << "\n=== Round-11 R-001/R-002: mixed-flux Adjacent MPI consistency"
         << " (np=2) ===\n";
   }

   Mesh serial = BuildOrthogonalFaultMesh();
   std::vector<int> part = BuildSkewedPartition(serial);
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
      std::cout
         << "  rank=0  fault_interior_faces.Size() = " << n_local_fault_int
         << "  fault_shared_faces.Size() = " << n_local_fault_shr << "\n";
   }
   if (rank == 1)
   {
      std::cout
         << "  rank=1  fault_interior_faces.Size() = " << n_local_fault_int
         << "  fault_shared_faces.Size() = " << n_local_fault_shr << "\n";
   }

   // Sanity: BOTH ranks should have fault as interior (each rank holds
   // 2 hexes, 1 on each side of the fault); neither has shared-fault
   // faces (fault perpendicular to rank seam).
   const int my_fixture_ok =
      (n_local_fault_int > 0) && (n_local_fault_shr == 0);
   int all_fixture_ok = 0;
   MPI_Allreduce(&my_fixture_ok, &all_fixture_ok, 1, MPI_INT, MPI_LAND, comm);
   if (rank == 0)
   {
      if (all_fixture_ok)
      {
         std::cout
            << "  [fixture-gate] both ranks own fault interior with no "
            << "shared-fault faces (fault perpendicular to rank seam) OK\n";
      }
      else
      {
         std::cerr
            << "  [fixture-gate] FAILED: expected both ranks to have "
            << "fault_interior_faces > 0 and no shared fault faces — "
            << "actual values logged above.  Test cannot proceed.\n";
      }
   }
   if (!all_fixture_ok)
   {
      MPI_Finalize();
      return 1;
   }

   // Engage Adjacent mode (the dispatch path R-001 fixed).
   wave.SetMixedFluxMode(MixedFluxMode::Adjacent);

   const std::unordered_set<int> &S = wave.GetCentralFluxFaceSet();
   if (rank == 0)
   {
      std::cout << "  rank=0  |central_flux_face_set_| = " << S.size() << "\n";
   }
   if (rank == 1)
   {
      std::cout << "  rank=1  |central_flux_face_set_| = " << S.size() << "\n";
   }

   // Cross-rank consistency check.
   //
   // Build per-shared-face (key, member) tuples on this rank.  The KEY
   // is the sorted global-vertex-id 4-tuple (same convention as the
   // production exchange in BuildCentralFluxFaceSet_).  The MEMBER is
   // 1 if this rank's central_flux_face_set_ contains the shared face,
   // 0 otherwise.  Allgatherv to every rank, then group by key and
   // verify all reports of any given key have the same MEMBER value.
   Array<HYPRE_BigInt> gvi;
   pmesh.GetGlobalVertexIndices(gvi);
   auto make_global_key = [&](const Array<int> &verts)
   {
      std::array<HYPRE_BigInt, 4> key = {0, 0, 0, 0};
      for (int v = 0; v < std::min(verts.Size(), 4); v++)
      {
         key[v] = gvi[verts[v]];
      }
      std::sort(key.begin(), key.begin() + verts.Size());
      return key;
   };

   const int n_shared = pmesh.GetNSharedFaces();
   std::vector<HYPRE_BigInt> local_flat;  // 5 ints per tuple: 4 key + 1 flag
   local_flat.reserve(static_cast<size_t>(n_shared) * 5);
   for (int sf = 0; sf < n_shared; sf++)
   {
      const int f = pmesh.GetSharedFace(sf);
      Array<int> verts;
      pmesh.GetFaceVertices(f, verts);
      auto key = make_global_key(verts);
      const HYPRE_BigInt member = (S.count(f) > 0) ? 1 : 0;
      for (int i = 0; i < 4; i++) { local_flat.push_back(key[i]); }
      local_flat.push_back(member);
   }
   const int my_size = static_cast<int>(local_flat.size());
   std::vector<int> sizes(nprocs), displs(nprocs);
   MPI_Allgather(&my_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);
   int total = 0;
   for (int r = 0; r < nprocs; r++) { displs[r] = total; total += sizes[r]; }
   std::vector<HYPRE_BigInt> all_flat(total);
   MPI_Allgatherv(local_flat.data(), my_size,
                  MPITypeMap<HYPRE_BigInt>::mpi_type,
                  all_flat.data(), sizes.data(), displs.data(),
                  MPITypeMap<HYPRE_BigInt>::mpi_type, comm);

   // R-1209: group by key into a VECTOR of reports (not a set), so we
   // verify the per-shared-face report count is exactly 2 (one from
   // each rank that owns the seam) AND all reports agree.  A `std::set`
   // collapsed both checks into "distinct values > 1", which would
   // silently pass if a future MFEM refactor changed shared-face
   // exchange to a one-sided broadcast.
   std::map<std::array<HYPRE_BigInt, 4>, std::vector<int>> per_key;
   for (int i = 0; i < total; i += 5)
   {
      std::array<HYPRE_BigInt, 4> key = {
         all_flat[i], all_flat[i+1], all_flat[i+2], all_flat[i+3]
      };
      const int member = static_cast<int>(all_flat[i+4]);
      per_key[key].push_back(member);
   }

   int num_inconsistent = 0;
   int num_shared_keys  = 0;
   int num_wrong_count  = 0;
   for (const auto &kv : per_key)
   {
      // Each shared face MUST be reported by exactly 2 ranks (the two
      // sharing it); anything else is a fixture or MFEM-semantics bug.
      if (kv.second.size() != 2)
      {
         num_wrong_count++;
         std::fprintf(stderr,
            "[R-1209] key reported %zu times (expected 2)\n",
            kv.second.size());
         continue;
      }
      if (kv.second[0] != kv.second[1]) { num_inconsistent++; }
      num_shared_keys++;
   }

   if (rank == 0)
   {
      std::cout << "  [consistency-gate] shared-face keys reported globally: "
                << num_shared_keys
                << "  inconsistent (rank-disagree on member): "
                << num_inconsistent
                << "  wrong-count (not exactly 2 reports): "
                << num_wrong_count << "\n";
   }

   int local_pass = (num_inconsistent == 0 && num_wrong_count == 0) ? 1 : 0;
   int all_pass = 0;
   MPI_Allreduce(&local_pass, &all_pass, 1, MPI_INT, MPI_LAND, comm);

   if (rank == 0)
   {
      if (all_pass)
      {
         std::cout
            << "  PASSED: every shared non-fault face has the same "
            << "central_flux_face_set_ membership on both ranks (R-001 "
            << "post-fix consistency holds)\n";
      }
      else
      {
         std::cout
            << "  FAILED: at least one shared non-fault face has "
            << "different central_flux_face_set_ membership across "
            << "ranks (R-001 bug is present)\n";
      }
      std::cout << "\n========================================\n"
                << "  Results: " << (all_pass ? 1 : 0) << " passed, "
                << (all_pass ? 0 : 1) << " failed out of 1 test\n"
                << "========================================\n";
   }

   MPI_Finalize();
   return all_pass ? 0 : 1;
}
