// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// shared_fault_fixtures.hpp — small meshes whose fault face is INTERIOR at np=1
// and SHARED at np=2, plus the deterministic partitioner that makes it so.
//
// Extracted (review R-006) from tests/parallel/test_bimaterial_seam_fault_np2.cpp
// so that test_fault_frame_bit_census.cpp and
// test_shared_fault_substep_parity_np2.cpp can share one definition instead of
// carrying three copies.  Header-only; no MPI calls at namespace scope.

#ifndef MFEM_SEAS_TESTS_SHARED_FAULT_FIXTURES_HPP
#define MFEM_SEAS_TESTS_SHARED_FAULT_FIXTURES_HPP

#include "mfem.hpp"

#include <vector>

namespace mfem
{
namespace seas
{
namespace test_fixtures
{

/// Two tets sharing the triangular fault face v0v1v2 on the y = 0 plane.
/// Interior faces are tagged attribute 3 (FAULT); exterior faces attribute 1.
/// At np=1 the fault face is INTERIOR; under `PartitionByYSign` at np=2 it is
/// SHARED.  That contrast is the whole point of the fixture.
inline Mesh BuildTwoTetSharedFault()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0,  0.0, 0.0},   // v0 (fault, y=0)
      {1.0,  0.0, 0.0},   // v1 (fault)
      {0.0,  0.0, 1.0},   // v2 (fault)
      {0.0,  1.0, 0.0},   // v3 apex (+y)
      {0.0, -1.0, 0.0},   // v4 apex (-y)
   };
   for (int v = 0; v < 5; ++v) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);   // elem 0: y < 0
   mesh.AddTet(0, 1, 2, 3, 1);   // elem 1: y > 0
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      const int attr = (ftr && ftr->Elem2No >= 0) ? 3 : 1;   // interior => fault
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], attr);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

/// 4 unit hexes: x columns [0,1] and [1,2]; y halves [-1,0] and [0,1]; z [0,1].
/// The y = 0 plane carries TWO interior faces (one per x column), both tagged
/// FAULT (attr 3).  Extracted from test_bimaterial_seam_fault_np2.cpp.
/// Element order: e0=(ix0,iy0) e1=(ix0,iy1) e2=(ix1,iy0) e3=(ix1,iy1).
inline Mesh BuildFourHexTwoFault()
{
   const real_t xs[3] = {0.0, 1.0, 2.0};
   const real_t ys[3] = {-1.0, 0.0, 1.0};
   const real_t zs[2] = {0.0, 1.0};
   Mesh mesh(3, 3 * 3 * 2, 4, 0);
   auto vid = [](int ix, int iy, int iz) { return (ix * 3 + iy) * 2 + iz; };
   for (int ix = 0; ix < 3; ++ix)
      for (int iy = 0; iy < 3; ++iy)
         for (int iz = 0; iz < 2; ++iz)
         {
            real_t v[3] = { xs[ix], ys[iy], zs[iz] };
            mesh.AddVertex(v);
         }
   for (int ix = 0; ix < 2; ++ix)
      for (int iy = 0; iy < 2; ++iy)
      {
         const int v[8] = {
            vid(ix,   iy,   0), vid(ix+1, iy,   0), vid(ix+1, iy+1, 0), vid(ix, iy+1, 0),
            vid(ix,   iy,   1), vid(ix+1, iy,   1), vid(ix+1, iy+1, 1), vid(ix, iy+1, 1)
         };
         mesh.AddHex(v, 1);
      }
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      real_t cy = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cy += mesh.GetVertex(fv[i])[1]; }
      cy /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cy) < 1e-9)   // y = 0 interior face => FAULT
         { mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3); }
         continue;
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);    // FREE (external)
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

/// TWO disconnected 2-tet diamonds (pair A at x in [0,1], pair B translated by
/// +2 in x), each carrying one TRIANGULAR y = 0 fault face (attr 3).  Element
/// order: eA0 (y<0), eA1 (y>0), eB0 (y<0), eB1 (y>0).  With
/// `PartitionOneInteriorOneShared` at np == 2, fault A is INTERIOR to rank 0
/// while fault B is SHARED — rank 0 then owns interior AND shared fault QPs,
/// so the shared substep-buffer index base (GetNumLocalFaultQPs()) is nonzero:
/// the configuration that distinguishes ABSOLUTE from locally-REBASED buffer
/// indexing (R-1304).  Both faults are congruent with identical seeds, so at
/// np = 1 they must evolve identically — a built-in consistency check.
/// (Triangular faces are required: the shared-fault reconcile key,
/// dynamic/shared_fault_key.hpp MakeFaceKey, keys on a sorted vertex TRIPLE
/// and does not support quad faults at np > 1.)
inline Mesh BuildTwoTetPairsTwoFaults()
{
   Mesh mesh(3, 10, 4, 0);
   const real_t shift = 2.0;
   const real_t base[5][3] = {
      {0.0,  0.0, 0.0},   // v0 (fault, y=0)
      {1.0,  0.0, 0.0},   // v1 (fault)
      {0.0,  0.0, 1.0},   // v2 (fault)
      {0.0,  1.0, 0.0},   // v3 apex (+y)
      {0.0, -1.0, 0.0},   // v4 apex (-y)
   };
   for (int p = 0; p < 2; ++p)
   {
      for (int v = 0; v < 5; ++v)
      {
         real_t xyz[3] = { base[v][0] + p * shift, base[v][1], base[v][2] };
         mesh.AddVertex(xyz);
      }
   }
   mesh.AddTet(0, 1, 2, 4, 1);       // eA0: pair A, y < 0
   mesh.AddTet(0, 1, 2, 3, 1);       // eA1: pair A, y > 0
   mesh.AddTet(5, 6, 7, 9, 1);       // eB0: pair B, y < 0
   mesh.AddTet(5, 6, 7, 8, 1);       // eB1: pair B, y > 0
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      const int attr = (ftr && ftr->Elem2No >= 0) ? 3 : 1;   // interior => fault
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], attr);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

#ifdef MFEM_USE_MPI
/// Split by element-centroid y-sign so every y = 0 fault face straddles the
/// seam.  At np == 1 every element MUST map to rank 0: a rank-1 entry makes
/// ParMesh post an MPI_Isend to a nonexistent rank (MPI_ERR_RANK).
inline ParMesh PartitionByYSign(Mesh &&smesh_in)
{
   Mesh smesh(std::move(smesh_in));
   const int ne     = smesh.GetNE();
   const int nprocs = Mpi::WorldSize();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; ++e)
   {
      if (nprocs < 2) { part[e] = 0; continue; }
      Array<int> ev; smesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); ++v) { cy += smesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      part[e] = (cy < 0.0) ? 0 : 1;
   }
   ParMesh pmesh(MPI_COMM_WORLD, smesh, part.data());
   smesh.Clear();
   return pmesh;
}

/// MIXED partition for BuildFourHexTwoFault at np == 2: only e3 = (ix1, iy1)
/// goes to rank 1, so
///   * x-column 0 fault face (e0|e1)  ->  INTERIOR on rank 0,
///   * x-column 1 fault face (e2|e3)  ->  SHARED across the seam.
/// Rank 0 then owns interior AND shared fault QPs simultaneously — the
/// configuration that distinguishes ABSOLUTE from locally-REBASED indexing
/// into the substep I_imp buffer (R-1304), which the 2-tet fixture cannot
/// (its shared offset base is always 0).  At np == 1 everything is rank 0.
inline ParMesh PartitionMixedInteriorShared(Mesh &&smesh_in)
{
   Mesh smesh(std::move(smesh_in));
   const int ne     = smesh.GetNE();
   const int nprocs = Mpi::WorldSize();
   std::vector<int> part(ne, 0);
   if (nprocs >= 2)
   {
      MFEM_VERIFY(ne == 4, "PartitionMixedInteriorShared expects the 4-hex "
                  "fixture; got ne=" << ne);
      part[3] = 1;   // e3 = (ix=1, iy=1)
   }
   ParMesh pmesh(MPI_COMM_WORLD, smesh, part.data());
   smesh.Clear();
   return pmesh;
}

/// Partition for BuildTwoTetPairsTwoFaults at np == 2: only eB1 goes to
/// rank 1, so fault A (eA0|eA1) is INTERIOR on rank 0 and fault B (eB0|eB1)
/// is SHARED across the seam.  At np == 1 everything is rank 0.
inline ParMesh PartitionOneInteriorOneShared(Mesh &&smesh_in)
{
   Mesh smesh(std::move(smesh_in));
   const int ne     = smesh.GetNE();
   const int nprocs = Mpi::WorldSize();
   std::vector<int> part(ne, 0);
   if (nprocs >= 2)
   {
      MFEM_VERIFY(ne == 4, "PartitionOneInteriorOneShared expects the 4-tet "
                  "two-pair fixture; got ne=" << ne);
      part[3] = 1;   // eB1 (pair B, y > 0)
   }
   ParMesh pmesh(MPI_COMM_WORLD, smesh, part.data());
   smesh.Clear();
   return pmesh;
}
#endif  // MFEM_USE_MPI

}  // namespace test_fixtures
}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_TESTS_SHARED_FAULT_FIXTURES_HPP
