// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 2D (2026-04-23 post-round-2 findings §C): shared D4-equivariant
// tet mesh fixture for Arm 1 probes, the arm3d state-sampling test, the
// Phase 1 fixture-equivariance test, and the first-step audit.
//
// Previously BuildD4Mesh was duplicated in three test files (with one
// slight drift risk per the round-2 review's unreviewed-areas note).
// This header hosts the single authoritative definition.
//
// Invariant (verified by test_d4_fixture_vertex_equivariance):
// orbit-paired tets have slot-wise y-mirror-equivariant vertex lists.
// The `add_fault` flag tags the y=L/2 interior plane with attr 3 (same
// convention as BuildM0FaultMesh).

#ifndef MFEM_SEAS_D4_TET_MESH_HPP
#define MFEM_SEAS_D4_TET_MESH_HPP

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

namespace mfem
{
namespace seas
{

/// Build a 2x2x2 D4-equivariant tet mesh with edge length `L`.
/// - Vertex positions: {0, L/2, L}^3 (27 vertices).
/// - 8 hexes, split into 6 tets each (48 tets total).
/// - Orient +1 for (iy=0) hexes, orient -1 for (iy=1) hexes — chosen so
///   that orbit-paired tets (centroid y-mirror reflections) have
///   slot-wise y-mirror-equivariant vertex lists.
/// - `add_fault == true`: every triangle at y = L/2 is added as a
///   boundary triangle with attr 3.  Every exterior triangle gets
///   attr 1.
/// - `add_fault == false`: only exterior triangles are added (attr 1);
///   the y = L/2 internal plane is left as an interior face with
///   bdr_attr = 0.
inline Mesh BuildD4Mesh(bool add_fault, real_t L = 1000.0)
{
   const int n = 2;
   const real_t h = L / n;
   Mesh mesh(3, (n+1)*(n+1)*(n+1), 0, 0, 3);
   auto vid = [&](int ix, int iy, int iz)
   { return ix + (n+1)*(iy + (n+1)*iz); };
   for (int iz = 0; iz <= n; iz++)
      for (int iy = 0; iy <= n; iy++)
         for (int ix = 0; ix <= n; ix++)
            mesh.AddVertex(ix*h, iy*h, iz*h);
   auto add_hex = [&](int ix, int iy, int iz, int orient)
   {
      const int v000=vid(ix,iy,iz),  v100=vid(ix+1,iy,iz);
      const int v010=vid(ix,iy+1,iz),v110=vid(ix+1,iy+1,iz);
      const int v001=vid(ix,iy,iz+1),v101=vid(ix+1,iy,iz+1);
      const int v011=vid(ix,iy+1,iz+1),v111=vid(ix+1,iy+1,iz+1);
      auto v = [&](int k) -> int
      {
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
   // R4-R002 Path A (2026-04-23): call CheckElementOrientation(true)
   // explicitly.  The default Mesh::Finalize(refine=false,
   // fix_orientation=true) parameter does NOT actually fix tet
   // orientations for mesh objects constructed via AddTet (verified
   // empirically — round-3 R-003 Jacobian audit saw 24/48 negatives
   // AFTER Finalize()).  An explicit CheckElementOrientation(true)
   // call forces MFEM to swap two vertices in every negatively-
   // oriented tet.  This inherently conflicts with slot-wise y-mirror
   // equivariance (y-reflection is orientation-reversing in 3D, so
   // no mesh with uniform positive orientation AND slot-wise y-mirror
   // equivariance exists).  AssertD4FixtureValid below reports which
   // invariants are preserved.
   mesh.CheckElementOrientation(/*fix_it=*/true);
   mesh.Finalize();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= fv.Size();
      if (ftr && ftr->Elem2No >= 0)
      {
         if (add_fault && std::abs(cy - 0.5 * L) < 1e-8)
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

/// Consolidated fixture-quality gate for D4-using tests (round-4 R-001 + R-007).
///
/// Single call that verifies ALL invariants Arm 1 / Arm 2 / arm3d / audit
/// probes depend on.  Each invariant is a separate check with its own
/// PASS/FAIL report; the overall return is 0 iff every check passes.
///
/// Invariants currently checked:
///   1. Jacobian positivity: every element's det(J) > 0.  Negative dets
///      cause DG volume rhs to sign-flip on orientation-reversed tets
///      (found in round-3 R-003).
///   2. Slot-wise y-mirror equivariance: for every orbit-paired (e_a, e_b)
///      pair, the four vertices of e_b are the y-mirrors of e_a's four
///      vertices in matching local-vertex slots.  Required for the
///      G_ORBIT_DOF_MAP basis sub-probe to test basis covariance.
///   3. Set-wise y-mirror equivariance: each vertex of e_a y-mirrors to
///      SOME vertex of e_b (weaker than slot-wise; required for vertex
///      coincidence at all).
///
/// Usage: call from main() as the first line of every D4-using test.
/// Abort-on-fail behavior: prints diagnostic to stderr, returns nonzero
/// exit code via `std::exit(2)` if `abort_on_fail == true` (default).
/// Tests may pass `abort_on_fail = false` to collect the verdict
/// without aborting (for diagnostic tests that themselves audit the
/// fixture).
///
/// Extensibility: add new invariants as additional `check_*` blocks.
/// Every D4-using test picks them up automatically.
inline int AssertD4FixtureValid(Mesh &mesh, real_t L,
                                bool abort_on_fail = true,
                                std::ostream &out = std::cout)
{
   const char *HDR = "[AssertD4FixtureValid]";
   int n_failures = 0;

   // ---- Check 1: Jacobian positivity ----
   {
      int n_pos = 0, n_neg = 0, n_zero = 0;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         ElementTransformation *Tr = mesh.GetElementTransformation(e);
         IntegrationPoint ip;  ip.Init(0);
         ip.x = 0.25; ip.y = 0.25; ip.z = 0.25; ip.weight = 1.0;
         Tr->SetIntPoint(&ip);
         const real_t detJ = Tr->Jacobian().Det();
         if (detJ > 0)      { n_pos++; }
         else if (detJ < 0) { n_neg++; }
         else               { n_zero++; }
      }
      const bool pass = (n_neg == 0 && n_zero == 0);
      out << "  " << HDR << " Jacobian positivity: "
          << (pass ? "PASS" : "FAIL")
          << " (pos=" << n_pos << " neg=" << n_neg
          << " zero=" << n_zero << ")\n";
      if (!pass) { n_failures++; }
   }

   // ---- Checks 2/3: y-mirror equivariance (slot-wise + set-wise) ----
   {
      const real_t s_scale = 1e6;
      auto q = [&](real_t v)
      { return static_cast<long long>(std::round(v*s_scale)); };

      // Orbit bucketing by (cx, min(cy,L-cy), cz).
      struct Centroid { int e; real_t cx, cy, cz; };
      std::map<std::tuple<long long,long long,long long>, std::vector<Centroid>> bk;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         Array<int> ev; mesh.GetElementVertices(e, ev);
         real_t cx = 0, cy = 0, cz = 0;
         for (int i = 0; i < ev.Size(); i++)
         {
            const real_t *v = mesh.GetVertex(ev[i]);
            cx += v[0]; cy += v[1]; cz += v[2];
         }
         cx /= ev.Size(); cy /= ev.Size(); cz /= ev.Size();
         bk[std::make_tuple(q(cx), q(std::min(cy, L-cy)), q(cz))]
            .push_back({e, cx, cy, cz});
      }

      int n_pairs_checked = 0, n_slot_mismatches = 0, n_setwise_mismatches = 0;
      for (const auto &kv : bk)
      {
         const auto &orbit = kv.second;
         if (orbit.size() < 2) { continue; }
         // Split into lower/upper halves and try to pair them.
         std::vector<int> lower, upper;
         for (const auto &c : orbit)
         {
            if (c.cy < 0.5 * L) { lower.push_back(c.e); }
            else                { upper.push_back(c.e); }
         }
         if (lower.empty() || upper.empty()) { continue; }

         // For each lower-upper pair: verify slot-wise + set-wise.
         for (int ea : lower)
         {
            for (int eb : upper)
            {
               Array<int> va, vb;
               mesh.GetElementVertices(ea, va);
               mesh.GetElementVertices(eb, vb);
               if (va.Size() != 4 || vb.Size() != 4) { continue; }
               n_pairs_checked++;

               // Slot-wise check: va[i] y-mirrors to vb[i].
               bool slotwise_ok = true;
               for (int i = 0; i < 4; i++)
               {
                  const real_t *xa = mesh.GetVertex(va[i]);
                  const real_t *xb = mesh.GetVertex(vb[i]);
                  if (std::abs(xa[0] - xb[0]) > 1e-10
                   || std::abs(xa[1] - (L - xb[1])) > 1e-10
                   || std::abs(xa[2] - xb[2]) > 1e-10)
                  {
                     slotwise_ok = false; break;
                  }
               }
               if (!slotwise_ok) { n_slot_mismatches++; }

               // Set-wise: every xa y-mirrors to SOME xb (weaker).
               bool setwise_ok = true;
               for (int i = 0; i < 4; i++)
               {
                  const real_t *xa = mesh.GetVertex(va[i]);
                  bool found = false;
                  for (int j = 0; j < 4; j++)
                  {
                     const real_t *xb = mesh.GetVertex(vb[j]);
                     if (std::abs(xa[0] - xb[0]) < 1e-10
                      && std::abs(xa[1] - (L - xb[1])) < 1e-10
                      && std::abs(xa[2] - xb[2]) < 1e-10)
                     { found = true; break; }
                  }
                  if (!found) { setwise_ok = false; break; }
               }
               if (!setwise_ok) { n_setwise_mismatches++; }
            }
         }
      }

      const bool slot_pass = (n_slot_mismatches == 0) && (n_pairs_checked > 0);
      const bool set_pass  = (n_setwise_mismatches == 0) && (n_pairs_checked > 0);
      out << "  " << HDR << " slot-wise y-mirror equivariance: "
          << (slot_pass ? "PASS" : "FAIL")
          << " (" << n_slot_mismatches << " mismatches / "
          << n_pairs_checked << " pairs)\n";
      out << "  " << HDR << " set-wise y-mirror equivariance:  "
          << (set_pass ? "PASS" : "FAIL")
          << " (" << n_setwise_mismatches << " mismatches / "
          << n_pairs_checked << " pairs)\n";

      if (n_pairs_checked == 0)
      {
         out << "  " << HDR << " WARNING: no orbit-paired tets found — "
                "fixture may not be a D4 mesh\n";
         n_failures++;
      }
      // Hard gates: Jacobian (above) and set-wise equivariance.
      // Set-wise is the minimum for orbit-paired DOF position matching.
      if (!set_pass) { n_failures++; }

      // Slot-wise is REPORTED but does NOT contribute to n_failures.
      // After round-4 R-002 (2026-04-23): in 3D, y-reflection is
      // orientation-reversing, so a D4 mesh cannot have BOTH uniform
      // positive orientation AND slot-wise y-mirror equivariance
      // simultaneously.  CheckElementOrientation(true) applied to a
      // slot-equivariant fixture breaks slot-wise by swapping two
      // vertices per orient-reversed tet.  Set-wise and Jacobian
      // positivity are the invariants every D4-consuming probe needs;
      // slot-wise is a DIFFERENT invariant the basis sub-probe
      // (G_ORBIT_DOF_MAP) needs.  Tests requiring slot-wise must
      // check it explicitly; they cannot rely on this gate alone.
      if (!slot_pass)
      {
         out << "  " << HDR << " NOTE: slot-wise equivariance is FAIL,\n"
             << "      but this is expected under Path A (orientation-\n"
             << "      uniform fixture).  Tests that require slot-wise\n"
             << "      (e.g. G_ORBIT_DOF_MAP basis sub-probe) must use\n"
             << "      a slot-equivariant fixture variant or skip.\n";
      }
   }

   out << "  " << HDR << " overall: "
       << (n_failures == 0 ? "PASS" : "FAIL")
       << " (" << n_failures << " invariant(s) failed)\n";

   if (n_failures > 0 && abort_on_fail)
   {
      std::cerr << "  " << HDR << " aborting — fixture invariants not "
                   "satisfied; re-run on a fixed BuildD4Mesh before "
                   "interpreting results.\n";
      std::exit(2);
   }
   return n_failures;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_D4_TET_MESH_HPP
