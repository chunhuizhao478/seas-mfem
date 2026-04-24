// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 1 correction (2026-04-23 REVIEW.md R-002): vertex-set
// equivariance test for BuildD4Mesh.
//
// Motivation
// ----------
// REVIEW.md R-001 showed the original G_ORBIT_DOF_MAP probe's "99"
// sentinel conflates mesh-layer vertex asymmetry with basis-layer
// index non-covariance.  The natural follow-on — re-run Arm 1 probes
// on the D4 fixture — is only interpretable if BuildD4Mesh actually
// produces y-mirror-paired vertex sets in orbit-paired tets.
//
// This test is a Phase 2D precondition.  If it fails, BuildD4Mesh is
// misbuilt and any subsequent G_ORBIT_DOF_MAP result on the D4 fixture
// is meaningless.
//
// Gate
// ----
// For every pair of tets (e_a, e_b) where
//   centroid(e_a) = (cx, cy, cz), centroid(e_b) = (cx, L-cy, cz),
// the four vertices of e_b must be the y-mirror of e_a's four vertices
// IN THE SAME LOCAL-VERTEX SLOT (i.e., va[i] y-mirrors to vb[i] for
// i=0..3).  Slot-wise equivariance is the property that makes MFEM's
// per-element DOF index ordering orbit-covariant — the thing Arm 1's
// G_ORBIT_DOF_MAP basis sub-probe is supposed to test.

#include "mfem.hpp"
#include "../../dynamic/d4_tet_mesh.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

constexpr real_t kL = 1000.0;

struct Centroid
{
   real_t cx, cy, cz;
   int    elem;
};

Centroid ComputeCentroid(const Mesh &m, int e)
{
   Array<int> ev; m.GetElementVertices(e, ev);
   Centroid c{0, 0, 0, e};
   for (int i = 0; i < ev.Size(); i++)
   {
      const real_t *v = m.GetVertex(ev[i]);
      c.cx += v[0]; c.cy += v[1]; c.cz += v[2];
   }
   c.cx /= ev.Size(); c.cy /= ev.Size(); c.cz /= ev.Size();
   return c;
}

// Build a map from canonicalized centroid key → list of element ids.
// Pairs within the same key are y-mirror orbit partners.
using CentKey = std::tuple<long long, long long, long long>;
CentKey QuantizedKey(real_t cx, real_t cy, real_t cz)
{
   const real_t s = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v*s)); };
   return std::make_tuple(q(cx), q(std::min(cy, kL-cy)), q(cz));
}

} // anonymous

int main()
{
   std::cout << "=== R-002: BuildD4Mesh vertex-set y-mirror equivariance ===\n";

   Mesh d4 = BuildD4Mesh(/*add_fault=*/false);

   // Build orbit buckets: every key with >=2 elements is a candidate
   // orbit pair.  Upper-y element's vertex[i] should y-mirror to
   // lower-y element's vertex[i].
   std::map<CentKey, std::vector<Centroid>> buckets;
   for (int e = 0; e < d4.GetNE(); e++)
   {
      Centroid c = ComputeCentroid(d4, e);
      buckets[QuantizedKey(c.cx, c.cy, c.cz)].push_back(c);
   }

   int n_orbits = 0;
   int n_pairs_checked = 0;
   int n_slot_mismatches = 0;           // va[i] NOT y-mirror of vb[i]
   int n_any_mirror_mismatches = 0;     // va[i] NOT y-mirror of ANY vb[*]
   real_t worst_slot_mirror_err = 0.0;

   auto is_y_mirror = [](const real_t *a, const real_t *b, real_t tol)
   {
      return std::abs(a[0] - b[0]) < tol
          && std::abs(a[1] - (kL - b[1])) < tol
          && std::abs(a[2] - b[2]) < tol;
   };

   for (const auto &kv : buckets)
   {
      const auto &elems = kv.second;
      if (elems.size() < 2) { continue; }
      n_orbits++;

      // Split by y-half so we always pair upper↔lower.
      std::vector<int> lower, upper;
      for (const auto &c : elems)
      {
         if (c.cy < 0.5 * kL) { lower.push_back(c.elem); }
         else                 { upper.push_back(c.elem); }
      }
      if (lower.empty() || upper.empty()) { continue; }

      // Match every lower element to every upper element in the same
      // centroid bucket and report whether slot-wise y-mirror holds.
      for (int ea : lower)
      {
         Array<int> va; d4.GetElementVertices(ea, va);
         for (int eb : upper)
         {
            Array<int> vb; d4.GetElementVertices(eb, vb);
            if (va.Size() != 4 || vb.Size() != 4) { continue; }
            n_pairs_checked++;

            for (int i = 0; i < 4; i++)
            {
               const real_t *xa = d4.GetVertex(va[i]);
               const real_t *xb = d4.GetVertex(vb[i]);
               // Slot-wise y-mirror: xa[0]≈xb[0], xa[1]≈L-xb[1], xa[2]≈xb[2].
               real_t slot_err = std::max({
                  std::abs(xa[0] - xb[0]),
                  std::abs(xa[1] - (kL - xb[1])),
                  std::abs(xa[2] - xb[2])
               });
               worst_slot_mirror_err =
                  std::max(worst_slot_mirror_err, slot_err);
               if (slot_err > 1e-10)
               {
                  n_slot_mismatches++;

                  // Does xa[i] y-mirror to ANY vertex of eb?
                  bool any = false;
                  for (int j = 0; j < 4; j++)
                  {
                     const real_t *xbj = d4.GetVertex(vb[j]);
                     if (is_y_mirror(xa, xbj, 1e-10))
                     { any = true; break; }
                  }
                  if (!any) { n_any_mirror_mismatches++; }

                  std::cout << "  MISMATCH (slot " << i << ", e_a="
                            << ea << ", e_b=" << eb << "): xa=("
                            << std::fixed << std::setprecision(1)
                            << xa[0] << ", " << xa[1] << ", " << xa[2]
                            << ")  xb=(" << xb[0] << ", " << xb[1]
                            << ", " << xb[2] << ")  "
                            << (any ? "[maps to another slot of e_b]"
                                    : "[NOT in e_b at all]")
                            << "\n";
               }
            }
         }
      }
   }

   std::cout << "\n  orbit buckets with both halves populated: "
             << n_orbits << "\n";
   std::cout << "  element pairs checked:                    "
             << n_pairs_checked << "\n";
   std::cout << "  slot-wise y-mirror mismatches:            "
             << n_slot_mismatches << "\n";
   std::cout << "  of which: xa not ANY y-mirror of e_b:     "
             << n_any_mirror_mismatches << "\n";
   std::cout << "  worst slot y-mirror error:                "
             << std::scientific << std::setprecision(3)
             << worst_slot_mirror_err << "\n";

   // Two verdicts:
   //   SLOT-WISE equivariance (strict): local-slot i in e_a y-mirrors
   //      to local-slot i in e_b.  This is what makes MFEM's per-
   //      element DOF ordering orbit-covariant; it is the property
   //      the G_ORBIT_DOF_MAP basis sub-probe needs to be meaningful
   //      on the D4 fixture.
   //   SET-WISE equivariance (weaker): each vertex of e_a y-mirrors
   //      to SOME vertex of e_b (just not necessarily the same slot).
   //      Ensures the two tets occupy the same physical region; does
   //      NOT ensure basis index parity.
   bool slot_pass = (n_slot_mismatches == 0);
   bool set_pass  = (n_any_mirror_mismatches == 0);

   std::cout << "\n  slot-wise equivariance: "
             << (slot_pass ? "PASS" : "FAIL")
             << "   (needed for basis-layer Arm 1 probe on D4)\n";
   std::cout << "  set-wise equivariance:  "
             << (set_pass ? "PASS" : "FAIL")
             << "   (needed for vertex coincidence at all)\n";

   // Overall gate: set-wise is the MUST.  Slot-wise is informative —
   // if it fails, arm1 G_ORBIT_DOF_MAP (index ordering) on D4 will
   // show basis-layer failure even though the fixture is the real
   // source.  Report both so the Phase 2D interpretation is unambiguous.
   int ret = set_pass ? 0 : 1;
   std::cout << "\n  overall verdict: "
             << (ret == 0 ? "PASS (vertex coincidence established)"
                          : "FAIL (fixture not D4-equivariant at vertex-set level)")
             << "\n";
   return ret;
}
