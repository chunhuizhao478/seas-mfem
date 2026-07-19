// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_fault_reorder.cpp — LTS Phase 3, Step 0 (P-006) construction gate.
//
// The wave operator reorders its interior fault-face list into cluster-
// contiguous order at construction when a per-element cluster-id vector is
// passed (⇔ lts != "off").  This test proves, on a tiny Cartesian fault
// fixture, that the reorder is a pure PERMUTATION of the per-face QP blocks:
//
//   (1) the reordered face list is a permutation of the canonical (lts-off)
//       list — same multiset of mesh faces;
//   (2) the reordered list is cluster-CONTIGUOUS (the per-face cluster sequence
//       is non-decreasing), and within a cluster is sorted by mesh-face id;
//   (3) the canonical-permutation table inverts the reorder exactly: the
//       per-QP fault_coords of the face at new slot i equal the canonical
//       operator's fault_coords of the face at canonical position perm[i]
//       (geometry preserved, per-face QP blocks intact);
//   (4) lts-off leaves the list in canonical order with an EMPTY perm.
//
// Everything the driver keys off the fault ordering (DOFData, nucleation,
// stations, ParaView, checkpoints) derives from GetFaultInteriorFaces(), so a
// per-face-block permutation of the geometry is exactly the invariance the
// still-GTS canary asserts end-to-end.

#include "mfem.hpp"

#include "../../dynamic/wave_operator.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int g_checks = 0, g_fails = 0;
#define CHECK(cond, msg)                                                     \
   do { ++g_checks; if (!(cond)) { ++g_fails;                               \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg); } } while (0)

namespace
{
constexpr real_t kL = 1000.0;
constexpr int    kOrder = 2;
constexpr real_t kLambda = 3.2044e10, kMu = 3.2038e10, kRho = 2670.0;

// 2×2×2 Cartesian tet mesh; fault attribute 3 on the interior triangles at
// y = L/2, free-surface attribute 1 elsewhere.
Mesh BuildCartesianFaultMesh(int n)
{
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON,
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
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * kL) < 1e-8)
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Per-QP physical coords of every interior fault face, walked in the operator's
// GetFaultInteriorFaces() order (mirrors the driver's fault_coords build).
std::vector<Vector> FaultCoords(const WaveOperator<Mesh> &wave, Mesh &mesh)
{
   std::vector<Vector> out;
   const Array<int> &faces = wave.GetFaultInteriorFaces();
   for (int i = 0; i < faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2 * kOrder);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         out.push_back(phys);
      }
   }
   return out;
}

// x-centroid of an interior fault face.
real_t FaceCx(Mesh &mesh, int f)
{
   Array<int> fv; mesh.GetFaceVertices(f, fv);
   real_t cx = 0;
   for (int v = 0; v < fv.Size(); v++) { cx += mesh.GetVertex(fv[v])[0]; }
   return cx / fv.Size();
}
}  // namespace

int main()
{
   const int n = 4;                    // 4×4×4 → 16 hex-faces → 32 fault tris
   Mesh mesh = BuildCartesianFaultMesh(n);
   BoundaryConfig bc; bc.fault_attr = 3;

   // -------- canonical (lts-off) operator --------
   WaveOperator<Mesh> wave_off(mesh, kOrder, kLambda, kMu, kRho, bc);
   const Array<int> faces_off = wave_off.GetFaultInteriorFaces();
   const int nfi = faces_off.Size();
   CHECK(nfi > 0, "fixture produced no interior fault faces");
   CHECK(!wave_off.FaultFacesReordered(),
         "lts-off operator reports a reorder");
   CHECK(wave_off.GetFaultFaceCanonicalPerm().empty(),
         "lts-off canonical perm is non-empty");

   // -------- per-element cluster ids: 3-way x-split, both fault-face
   //          elements forced same-cluster (D-1). --------
   std::vector<int> cluster(mesh.GetNE(), 0);
   for (int i = 0; i < nfi; i++)
   {
      const int f = faces_off[i];
      auto *ftr = mesh.GetInteriorFaceTransformations(f);
      const real_t cx = FaceCx(mesh, f);
      int c = static_cast<int>(std::floor(3.0 * cx / kL));
      c = std::max(0, std::min(2, c));
      cluster[ftr->Elem1No] = c;
      cluster[ftr->Elem2No] = c;
   }

   // -------- reordered (lts-on) operator --------
   WaveOperator<Mesh> wave_on(mesh, kOrder, kLambda, kMu, kRho, bc, &cluster);
   const Array<int> faces_on = wave_on.GetFaultInteriorFaces();
   const std::vector<int> &perm = wave_on.GetFaultFaceCanonicalPerm();

   CHECK(wave_on.FaultFacesReordered(), "lts-on operator did not reorder");
   CHECK(faces_on.Size() == nfi, "reorder changed the fault-face count");
   CHECK(static_cast<int>(perm.size()) == nfi,
         "canonical perm length != fault-face count");

   // (1) permutation of the same face multiset.
   {
      Array<int> a(faces_off), b(faces_on);
      a.Sort(); b.Sort();
      bool same = (a.Size() == b.Size());
      for (int i = 0; same && i < a.Size(); i++) { same = (a[i] == b[i]); }
      CHECK(same, "reordered list is not a permutation of the canonical list");
   }

   // Face -> cluster lookup (via the canonical operator's element data).
   auto face_cluster = [&](int f)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(f);
      return cluster[ftr->Elem1No];
   };

   // (2) cluster-contiguity + within-cluster mesh-face-id sort.
   {
      bool ok = true;
      for (int i = 1; i < nfi; i++)
      {
         const int c0 = face_cluster(faces_on[i - 1]);
         const int c1 = face_cluster(faces_on[i]);
         if (c1 < c0) { ok = false; break; }
         if (c1 == c0 && faces_on[i] < faces_on[i - 1]) { ok = false; break; }
      }
      CHECK(ok, "reordered list is not (cluster, face-id)-contiguous");
   }

   // Non-triviality: >1 cluster present AND the perm actually moved a face
   // (else the test would be vacuously satisfied by an identity reorder).
   {
      bool multi = false, moved = false;
      for (int i = 0; i < nfi; i++)
      {
         if (face_cluster(faces_on[i]) != face_cluster(faces_on[0])) { multi = true; }
         if (perm[i] != i) { moved = true; }
      }
      CHECK(multi, "fixture spans only one cluster (test is weak)");
      CHECK(moved, "reorder was the identity (test is weak)");
   }

   // (3) canonical-perm inversion + geometry invariance (per-face QP blocks).
   {
      // perm[i] = canonical position of the face now at new slot i.
      bool face_ok = true;
      for (int i = 0; i < nfi; i++)
      {
         if (faces_on[i] != faces_off[perm[i]]) { face_ok = false; break; }
      }
      CHECK(face_ok, "canonical perm does not map new slot -> canonical face");

      const std::vector<Vector> co_off = FaultCoords(wave_off, mesh);
      const std::vector<Vector> co_on  = FaultCoords(wave_on, mesh);
      const int nqp = static_cast<int>(co_off.size()) / nfi;
      CHECK(nqp > 0 && static_cast<int>(co_off.size()) == nfi * nqp,
            "fault_coords size not a multiple of face count");
      CHECK(co_on.size() == co_off.size(), "fault_coords size changed");

      real_t md = 0.0;
      for (int i = 0; i < nfi; i++)
      {
         const int cpos = perm[i];
         for (int q = 0; q < nqp; q++)
         {
            const Vector &a = co_on[i * nqp + q];
            const Vector &b = co_off[cpos * nqp + q];
            for (int d = 0; d < 3; d++)
            { md = std::max(md, std::abs(a(d) - b(d))); }
         }
      }
      CHECK(md == 0.0,
            "per-face QP-block geometry not preserved under the reorder");
   }

   std::printf("test_lts_fault_reorder: %d checks, %d failures\n",
               g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}
