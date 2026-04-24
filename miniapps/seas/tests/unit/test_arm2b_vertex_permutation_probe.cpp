// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 3 STOP — Arm 2b: vertex-permutation verification probe.
//
// Hypothesis: on the hand-built D4-equivariant fixture, MFEM's
// `FinalizeMesh` / `SetAttributes` reorders tet vertex indices so
// orbit-pair elements no longer have mirror-related vertex lists, even
// though their nodal-DOF coordinates still do.  That silent permutation
// is what breaks `CalcOrtho`-based face normals across orbit.
//
// Probe: for each orbit pair (e_a, e_b) under y = L/2 mirror, examine
//   (1) the vertex coordinates stored for each element
//   (2) whether the per-element face normals (from CalcOrtho at QP=0)
//       on orbit-matched local sides are y-mirror images of each other
//   (3) the face-Jacobian determinant sign on matched faces
//
// If (2) reveals non-mirror-related normals, the hypothesis is confirmed:
// MFEM's internal storage doesn't preserve our D4 structure, and any
// flux dispatch reading normals via CalcOrtho will see inconsistent
// orientations across orbit.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define VERDICT(name, ok, extra) do {                                     \
   num_tests++;                                                            \
   if (ok) { num_passed++;                                                 \
      std::cout << "PROBE " << (name) << ": PASS " << extra << "\n"; }     \
   else { num_failed++;                                                    \
      std::cout << "PROBE " << (name) << ": FAIL " << extra << "\n"; }     \
} while (0)

namespace
{

constexpr real_t kL = 1000.0;

// Same D4 builder as Arm 2.
Mesh BuildD4EquivariantMesh()
{
   const int n = 2;
   const real_t h = kL / n;
   Mesh mesh(3, (n+1)*(n+1)*(n+1), 0, 0, 3);
   auto vid = [&](int ix, int iy, int iz) {
      return ix + (n+1) * (iy + (n+1) * iz);
   };
   for (int iz = 0; iz <= n; iz++)
      for (int iy = 0; iy <= n; iy++)
         for (int ix = 0; ix <= n; ix++)
            mesh.AddVertex(ix*h, iy*h, iz*h);

   auto add_hex = [&](int ix, int iy, int iz, int orient) {
      const int v000=vid(ix,iy,iz),v100=vid(ix+1,iy,iz),v010=vid(ix,iy+1,iz),v110=vid(ix+1,iy+1,iz);
      const int v001=vid(ix,iy,iz+1),v101=vid(ix+1,iy,iz+1),v011=vid(ix,iy+1,iz+1),v111=vid(ix+1,iy+1,iz+1);
      auto v = [&](int k) -> int {
         if (orient > 0) {
            switch(k){case 0:return v000;case 1:return v100;case 2:return v010;case 3:return v110;
                      case 4:return v001;case 5:return v101;case 6:return v011;case 7:return v111;}
         } else {
            switch(k){case 0:return v010;case 1:return v110;case 2:return v000;case 3:return v100;
                      case 4:return v011;case 5:return v111;case 6:return v001;case 7:return v101;}
         }
         return -1;
      };
      int T[6][4] = {
         {v(0),v(1),v(3),v(7)},{v(0),v(2),v(3),v(7)},{v(0),v(2),v(6),v(7)},
         {v(0),v(4),v(6),v(7)},{v(0),v(4),v(5),v(7)},{v(0),v(1),v(5),v(7)}
      };
      for (int t = 0; t < 6; t++) mesh.AddTet(T[t], 1);
   };

   for (int iz = 0; iz < n; iz++)
      for (int ix = 0; ix < n; ix++) {
         add_hex(ix, 0, iz, +1);
         add_hex(ix, 1, iz, -1);
      }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

} // anonymous

int main()
{
   std::cout << "=== Arm 2b: vertex-permutation verification on D4 fixture ===\n";

   Mesh mesh = BuildD4EquivariantMesh();
   std::cout << "  NE=" << mesh.GetNE() << "  NV=" << mesh.GetNV()
             << "  NF=" << mesh.GetNumFaces() << "\n";

   // Dump element vertex indices + coords for every element.
   const real_t scale = 1e6;
   auto q = [&](real_t v) { return (long long)std::round(v*scale); };

   // Build orbit buckets by centroid (cx, min(cy,L-cy), cz).
   struct OrbitKey { long long x, ymir, z; };
   std::map<std::tuple<long long,long long,long long>, std::vector<int>> buckets;
   for (int e = 0; e < mesh.GetNE(); e++) {
      Array<int> ev; mesh.GetElementVertices(e, ev);
      real_t cx=0,cy=0,cz=0;
      for (int i = 0; i < ev.Size(); i++) {
         const real_t *v = mesh.GetVertex(ev[i]);
         cx+=v[0]; cy+=v[1]; cz+=v[2];
      }
      cx/=ev.Size(); cy/=ev.Size(); cz/=ev.Size();
      real_t cy_mir = std::min(cy, kL-cy);
      buckets[std::make_tuple(q(cx),q(cy_mir),q(cz))].push_back(e);
   }

   std::cout << "\n-- Orbit-pair vertex-list analysis --\n";
   int total_pairs = 0;
   int mirror_preserving = 0;
   int mirror_breaking = 0;
   int examined_detail = 0;

   for (const auto &kv : buckets) {
      if (kv.second.size() < 2) continue;
      total_pairs++;
      int ea = kv.second[0];
      int eb = kv.second[1];

      Array<int> va, vb;
      mesh.GetElementVertices(ea, va);
      mesh.GetElementVertices(eb, vb);

      // Compute the vertex COORDINATES (after any finalization permutation).
      std::array<std::array<real_t,3>,4> Pa, Pb;
      for (int i = 0; i < 4; i++) {
         const real_t *pa = mesh.GetVertex(va[i]);
         const real_t *pb = mesh.GetVertex(vb[i]);
         for (int j = 0; j < 3; j++) { Pa[i][j] = pa[j]; Pb[i][j] = pb[j]; }
      }

      // Check: is Pb[i] equal to y-mirror(Pa[i]) for every i?
      bool ordered_mirror = true;
      for (int i = 0; i < 4; i++) {
         const real_t ax = Pa[i][0], ay = Pa[i][1], az = Pa[i][2];
         const real_t bx = Pb[i][0], by = Pb[i][1], bz = Pb[i][2];
         const real_t amy = kL - ay;  // y-mirror of ax,ay,az
         if (std::abs(bx - ax) > 1e-8 || std::abs(by - amy) > 1e-8 ||
             std::abs(bz - az) > 1e-8) {
            ordered_mirror = false; break;
         }
      }

      // Weaker: is {Pb[.]} the same SET as {y-mirror(Pa[.])} but possibly in
      // a different order?  Match pairwise by coordinate equality.
      bool set_mirror = true;
      std::array<bool,4> used = {false,false,false,false};
      for (int i = 0; i < 4; i++) {
         const real_t ax = Pa[i][0], amy = kL - Pa[i][1], az = Pa[i][2];
         int found = -1;
         for (int j = 0; j < 4; j++) {
            if (used[j]) continue;
            if (std::abs(Pb[j][0]-ax)<1e-8 && std::abs(Pb[j][1]-amy)<1e-8 &&
                std::abs(Pb[j][2]-az)<1e-8) { found = j; used[j] = true; break; }
         }
         if (found < 0) { set_mirror = false; break; }
      }

      if (ordered_mirror)      mirror_preserving++;
      else if (set_mirror)     mirror_breaking++;   // set match but ordering permuted
      else                     {} // not even set-matched

      if (examined_detail < 3) {
         std::cout << "\n  orbit key (" << std::get<0>(kv.first) << ","
                   << std::get<1>(kv.first) << "," << std::get<2>(kv.first)
                   << ")  ea=" << ea << "  eb=" << eb << "\n";
         std::cout << "    ea verts: ";
         for (int i = 0; i < 4; i++) {
            std::cout << "[" << std::fixed << std::setprecision(1)
                      << Pa[i][0] << "," << Pa[i][1] << "," << Pa[i][2] << "] ";
         }
         std::cout << "\n    eb verts: ";
         for (int i = 0; i < 4; i++) {
            std::cout << "[" << std::fixed << std::setprecision(1)
                      << Pb[i][0] << "," << Pb[i][1] << "," << Pb[i][2] << "] ";
         }
         std::cout << "\n    ordered-mirror: " << (ordered_mirror ? "YES" : "NO")
                   << "  set-mirror: " << (set_mirror ? "YES" : "NO") << "\n";
         examined_detail++;
      }
   }

   std::cout << "\n  total orbit pairs:       " << total_pairs << "\n";
   std::cout << "  ordered-mirror preserving: " << mirror_preserving << "\n";
   std::cout << "  set-mirror but reordered:  " << mirror_breaking << "\n";
   std::cout << "  neither:                   "
             << (total_pairs - mirror_preserving - mirror_breaking) << "\n";

   VERDICT("ORDERED_VERTEX_MIRROR",
           mirror_preserving == total_pairs,
           "(preserving=" + std::to_string(mirror_preserving) + "/" +
           std::to_string(total_pairs) + ")");
   VERDICT("SET_VERTEX_MIRROR",
           (mirror_preserving + mirror_breaking) == total_pairs,
           "(set-match=" + std::to_string(mirror_preserving + mirror_breaking)
           + "/" + std::to_string(total_pairs) + ")");

   // --- Face-normal orientation check -------------------------------------
   // For orbit-pair (ea, eb) examine CalcOrtho at QP=0 on each of the 4
   // local sides.  We match sides by their face centroid (under y-mirror);
   // on matched sides, the CalcOrtho normals should be y-mirror images.
   std::cout << "\n-- Face-normal orientation check --\n";

   auto face_centroid = [&](int face_idx) {
      Array<int> fv; mesh.GetFaceVertices(face_idx, fv);
      std::array<real_t,3> c{0,0,0};
      for (int i = 0; i < fv.Size(); i++) {
         const real_t *v = mesh.GetVertex(fv[i]);
         c[0]+=v[0]; c[1]+=v[1]; c[2]+=v[2];
      }
      if (fv.Size() > 0) { c[0]/=fv.Size(); c[1]/=fv.Size(); c[2]/=fv.Size(); }
      return c;
   };

   int pair_count = 0;
   int pairs_with_matching_normals = 0;
   int total_matched_sides = 0;
   int sides_with_mirror_normal = 0;
   int max_shown = 2;

   for (const auto &kv : buckets) {
      if (kv.second.size() < 2) continue;
      pair_count++;
      int ea = kv.second[0];
      int eb = kv.second[1];

      Array<int> fa, ob, fb, obb;
      mesh.GetElementFaces(ea, fa, ob);
      mesh.GetElementFaces(eb, fb, obb);

      int matched = 0, mirrored = 0;
      bool show = (pair_count <= max_shown);

      // For each side of ea, find the side of eb whose face centroid is
      // the y-mirror of ea's side centroid.
      for (int sa = 0; sa < fa.Size(); sa++) {
         auto ca = face_centroid(fa[sa]);
         int sb_match = -1;
         for (int sb = 0; sb < fb.Size(); sb++) {
            auto cb = face_centroid(fb[sb]);
            if (std::abs(cb[0]-ca[0])<1e-8 &&
                std::abs(cb[1]-(kL-ca[1]))<1e-8 &&
                std::abs(cb[2]-ca[2])<1e-8) { sb_match = sb; break; }
         }
         if (sb_match < 0) continue;
         matched++;

         // CalcOrtho at IP(0) on each side.
         auto *ftra = const_cast<Mesh&>(mesh).GetFaceElementTransformations(fa[sa]);
         auto *ftrb = const_cast<Mesh&>(mesh).GetFaceElementTransformations(fb[sb_match]);
         if (!ftra || !ftrb) continue;

         const IntegrationRule &ir = IntRules.Get(ftra->GetGeometryType(), 2);
         const IntegrationPoint &ip = ir.IntPoint(0);
         ftra->SetAllIntPoints(&ip);
         ftrb->SetAllIntPoints(&ip);

         Vector na(3), nb(3);
         CalcOrtho(ftra->Face->Jacobian(), na);
         CalcOrtho(ftrb->Face->Jacobian(), nb);

         // Expect nb ≈ y-mirror(na) = (na[0], -na[1], na[2])
         real_t dx = nb(0) - na(0);
         real_t dy = nb(1) - (-na(1));
         real_t dz = nb(2) - na(2);
         real_t err = std::sqrt(dx*dx + dy*dy + dz*dz);

         bool is_mirror = err < 1e-6 * std::max({na.Norml2(), nb.Norml2(), real_t(1.0)});
         if (is_mirror) mirrored++;

         if (show && matched <= 2) {
            std::cout << "  ea=" << ea << " sa=" << sa
                      << " n=(" << std::fixed << std::setprecision(3)
                      << na(0) << "," << na(1) << "," << na(2) << ")"
                      << "   eb=" << eb << " sb=" << sb_match
                      << " n=(" << nb(0) << "," << nb(1) << "," << nb(2) << ")"
                      << "  mirror?=" << (is_mirror ? "YES" : "NO")
                      << "  err=" << std::scientific << std::setprecision(2) << err << "\n";
         }
      }

      total_matched_sides   += matched;
      sides_with_mirror_normal += mirrored;
      if (matched > 0 && matched == mirrored) pairs_with_matching_normals++;
   }

   std::cout << "\n  pairs examined:              " << pair_count << "\n";
   std::cout << "  matched-side count (total):  " << total_matched_sides << "\n";
   std::cout << "  sides with mirror normal:    " << sides_with_mirror_normal << "\n";
   std::cout << "  pairs with ALL sides mirrored: " << pairs_with_matching_normals << "\n";

   VERDICT("FACE_NORMAL_ORBIT_COVARIANT",
           sides_with_mirror_normal == total_matched_sides,
           "(" + std::to_string(sides_with_mirror_normal) + "/" +
           std::to_string(total_matched_sides) + ")");

   std::cout << "\n=== Summary ===\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";

   std::cout << "\nInterpretation:\n";
   std::cout << "  SET_VERTEX_MIRROR PASS, ORDERED_VERTEX_MIRROR FAIL\n"
             << "    ==> MFEM permuted vertex lists (vertex SET is mirror\n"
             << "        but ordering was reshuffled by Finalize).  The\n"
             << "        reshuffling is the hypothesis's smoking gun.\n";
   std::cout << "  FACE_NORMAL_ORBIT_COVARIANT FAIL\n"
             << "    ==> CalcOrtho produces non-mirror normals on\n"
             << "        orbit-matched faces.  Flux dispatches reading\n"
             << "        normals via CalcOrtho WILL see asymmetric\n"
             << "        orientations across orbit — direct evidence of\n"
             << "        the mechanism behind Gate 14/14'.\n";
   std::cout << "  Both PASS ==> hypothesis wrong; bug is elsewhere.\n";

   return (num_failed == 0) ? 0 : 1;
}
