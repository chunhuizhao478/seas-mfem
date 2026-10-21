// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_fixtures.hpp — what the poroelastic phase tests share (plan rev 3,
// section 6.1): the rank-aware check macro and comparison helpers, and the
// meshes.
//
//   F-B  TPV102Box(nx, ny, nz): the TPV102 inline box (24 tets per hex, fault
//        attribute 3 on Y = 0), as in test_nitsche_hook_spy.cpp.
//   F-C  MarkSidesFC(mesh): F-B with element attribute 1 for x < 0 (medium A)
//        and 2 for x > 0 (medium B).  With nx even, x = 0 is a plane of faces,
//        so the material jump lies on interior faces.
//
// Moved from test_poro_p0.cpp in step P1.1 (unchanged), so every phase test
// uses one copy.

#ifndef MFEM_SEAS_PORO_FIXTURES_HPP
#define MFEM_SEAS_PORO_FIXTURES_HPP

#include "mfem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

using namespace mfem;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;

static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static long LSum(long v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG, MPI_SUM, c); return v; }

static long BitDiff(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return -1; }
   long n = 0;
   for (int i = 0; i < a.Size(); i++)
   {
      if (std::memcmp(&a(i), &b(i), sizeof(real_t)) != 0) { n++; }
   }
   return n;
}

/// max_i |a_i - b_i| / max_i |b_i| over all ranks.
static real_t RelDiff(const Vector &a, const Vector &b, MPI_Comm c)
{
   real_t num = 0.0, den = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      num = std::max(num, std::abs(a(i) - b(i)));
      den = std::max(den, std::abs(b(i)));
   }
   num = RMax(num, c);
   den = RMax(den, c);
   return den > 0.0 ? num / den : num;
}

// ---------------------------------------------------------------- fixture
// F-B of plan rev 3, section 6.1: the TPV102 inline box of
// test_nitsche_hook_spy.cpp (24 tets per hex, fault on Y = 0).

static Mesh *TPV102Box(int nx, int ny, int nz)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2 * kL, 2 * kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   const real_t tol = kL / nz * 0.1;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < tol) ? 1 : 5);
   }
   const real_t ytol = 1e-6 * kL;
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         if (std::abs(mesh->GetVertex(fv[i])[1]) > ytol) { on_plane = false; }
      }
      if (on_plane) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}

/// F-C: element attribute 1 where the centroid has x < 0, 2 where x > 0.
static void MarkSidesFC(Mesh &mesh)
{
   Vector x(3);
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      mesh.GetElementTransformation(e)->Transform(
         Geometries.GetCenter(mesh.GetElementBaseGeometry(e)), x);
      mesh.SetAttribute(e, x(0) < 0.0 ? 1 : 2);
   }
   mesh.SetAttributes();
}

#endif // MFEM_SEAS_PORO_FIXTURES_HPP
