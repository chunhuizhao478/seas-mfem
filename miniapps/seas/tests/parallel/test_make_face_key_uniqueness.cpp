// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug unit-test plan (2026-04-22) — Phase 4 P4-T1:
// MakeFaceKey uniqueness probe on shared fault faces.
//
// ============================================================================
// Motivation
// ============================================================================
// Debug plan §16.1 identifies `MakeFaceKey` as the **integer-exact**
// cross-rank pairing primitive that BP5 uses but TPV102 does not yet
// route into its runtime path.  Under a correct ParMesh partitioning,
// every shared fault face has exactly two sharing ranks, and the
// FaceVertexKey produced by `dynamic::MakeFaceKey` on each side must
// be bit-identical (sorted triple of global vertex IDs).
//
// This test:
//   1. Builds a 4-km Cartesian tet fixture (same as Phase D).
//   2. np=4.  Each rank enumerates its shared fault faces, computes
//      MakeFaceKey on each.
//   3. Allgathers all keys + source ranks.
//   4. Asserts EVERY observed key appears exactly TWICE across the
//      global pool (once per sharing rank).  A key appearing once
//      or three+ times indicates a ParMesh global-vertex-index mismatch.
//
// ============================================================================
// Usage:  mpirun -np 4 ./seas_test_make_face_key_uniqueness

#include "mfem.hpp"
#include "../../dynamic/shared_fault_key.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static Mesh BuildFixtureMesh()
{
   const int nx = 2, ny = 2, nz = 2;
   const double Lx = 2.0e3, Ly = 2.0e3, Lz = 4.0e3;
   Mesh mesh = Mesh::MakeCartesian3D(2*nx, 2*ny, nz,
                                     Element::TETRAHEDRON,
                                     2.0*Lx, 2.0*Ly, Lz);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *x = mesh.GetVertex(v);
      x[0] -= Lx; x[1] -= Ly; x[2] -= Lz;
   }
   const real_t tol = 1e-6;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      auto *Tr = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetIntPoint(&ip);
      Vector c(3); Tr->Transform(ip, c);
      if (std::abs(c(2)) < tol) { mesh.SetBdrAttribute(be, 1); }
      else                      { mesh.SetBdrAttribute(be, 5); }
   }
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(f);
      if (!ftr) { continue; }
      const IntegrationPoint &ip = Geometries.GetCenter(ftr->GetGeometryType());
      ftr->Face->SetIntPoint(&ip);
      Vector c(3); ftr->Face->Transform(ip, c);
      if (std::abs(c(1)) > tol) { continue; }
      Array<int> verts; mesh.GetFaceVertices(f, verts);
      if (verts.Size() == 3)
      { mesh.AddBdrTriangle(verts[0], verts[1], verts[2], 3); }
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (rank == 0)
   {
      std::cout << "\n=== TPV102 Pepper-Bug Phase 4 P4-T1: "
                << "MakeFaceKey Uniqueness ===\n  MPI ranks: "
                << nprocs << "\n";
   }

   if (nprocs < 2)
   {
      if (rank == 0)
      { std::cout << "  SKIPPED: requires np >= 2\n"; }
      MPI_Finalize();
      return 77;
   }

   Mesh serial_mesh = BuildFixtureMesh();
   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);

   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   const int n_shr = shr_faces.Size();

   int total_shared = 0;
   MPI_Allreduce(&n_shr, &total_shared, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::cout << "  total shared fault faces across ranks: "
                << total_shared << "\n";
   }
   if (total_shared == 0)
   {
      if (rank == 0) { std::cout << "  SKIPPED: no shared fault faces\n"; }
      MPI_Finalize();
      return 77;
   }

   // Compute MakeFaceKey on this rank's shared fault faces.
   Array<HYPRE_BigInt> gvi;
   pmesh.GetGlobalVertexIndices(gvi);

   std::vector<HYPRE_BigInt> local_keys_flat; // 3 per face
   std::vector<int>          local_src_rank;  // 1 per face
   for (int i = 0; i < n_shr; i++)
   {
      const int local_face = pmesh.GetSharedFace(shr_faces[i]);
      dynamic::FaceVertexKey k =
         dynamic::MakeFaceKey(local_face, gvi, pmesh);
      local_keys_flat.push_back(k.v[0]);
      local_keys_flat.push_back(k.v[1]);
      local_keys_flat.push_back(k.v[2]);
      local_src_rank.push_back(rank);
   }

   int my_size = (int)local_src_rank.size();
   std::vector<int> sizes(nprocs), key_byte_sizes(nprocs), key_byte_displs(nprocs);
   std::vector<int> src_byte_sizes(nprocs), src_byte_displs(nprocs);
   MPI_Allgather(&my_size, 1, MPI_INT, sizes.data(), 1, MPI_INT,
                 MPI_COMM_WORLD);
   int total_faces = 0;
   for (int r = 0; r < nprocs; r++)
   {
      key_byte_sizes[r]  = sizes[r] * 3 * (int)sizeof(HYPRE_BigInt);
      key_byte_displs[r] = total_faces * 3 * (int)sizeof(HYPRE_BigInt);
      src_byte_sizes[r]  = sizes[r] * (int)sizeof(int);
      src_byte_displs[r] = total_faces * (int)sizeof(int);
      total_faces       += sizes[r];
   }

   std::vector<HYPRE_BigInt> all_keys(total_faces * 3);
   std::vector<int>          all_src(total_faces);
   MPI_Allgatherv(local_keys_flat.data(),
                  my_size * 3 * (int)sizeof(HYPRE_BigInt),
                  MPI_BYTE,
                  all_keys.data(),
                  key_byte_sizes.data(), key_byte_displs.data(),
                  MPI_BYTE, MPI_COMM_WORLD);
   MPI_Allgatherv(local_src_rank.data(),
                  my_size * (int)sizeof(int), MPI_BYTE,
                  all_src.data(),
                  src_byte_sizes.data(), src_byte_displs.data(),
                  MPI_BYTE, MPI_COMM_WORLD);

   // Histogram by key.  Correct behaviour: every key appears exactly twice.
   std::map<dynamic::FaceVertexKey, std::vector<int>> by_key;
   for (int f = 0; f < total_faces; f++)
   {
      dynamic::FaceVertexKey k;
      k.v[0] = all_keys[3*f + 0];
      k.v[1] = all_keys[3*f + 1];
      k.v[2] = all_keys[3*f + 2];
      by_key[k].push_back(all_src[f]);
   }

   int n_bad_count = 0;
   int n_self_pair = 0;
   int worst_count = 0;
   for (const auto &p : by_key)
   {
      worst_count = std::max(worst_count, (int)p.second.size());
      if (p.second.size() != 2) { n_bad_count++; }
      // Also ensure the two ranks are DISTINCT (can't pair a rank to itself)
      if (p.second.size() == 2 && p.second[0] == p.second[1]) { n_self_pair++; }
   }

   int exit_code = 0;
   if (rank == 0)
   {
      std::cout << "  unique keys observed : " << by_key.size() << "\n";
      std::cout << "  worst count per key  : " << worst_count << "\n";
      std::cout << "  keys with count != 2 : " << n_bad_count << "\n";
      std::cout << "  keys paired to self  : " << n_self_pair << "\n";

      const bool pass_count  = (n_bad_count == 0);
      const bool pass_distinct = (n_self_pair == 0);
      const bool pass_worst  = (worst_count == 2);
      if (pass_count && pass_distinct && pass_worst)
      {
         std::cout << "  PASSED: every MakeFaceKey appears exactly twice, "
                   << "on distinct ranks\n";
      }
      else
      {
         std::cout << "  FAILED: MakeFaceKey cross-rank pairing invariant\n";
         exit_code = 1;
      }
      std::cout << "========================================\n";
   }

   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
}
