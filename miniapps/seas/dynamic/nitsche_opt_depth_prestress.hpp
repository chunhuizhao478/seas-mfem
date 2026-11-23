// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_opt_depth_prestress.hpp — option O1 of the complex-geometries plan
// (rev 2, Section 5.2; Chunhui Zhao 2026-10-20): a prestress tensor that
// varies with depth, for TPV12 and TPV22.
//
// Deck (Nitsche driver only; absent = every bit as before):
//
//   [prestress_depth_table]
//   stress_measure = "effective"         # or "total"; required
//   depth_m        = [0.0, 11951.15, 30000.0]   # nodes, strictly increasing
//   sigma_xx_pa    = [...]               # all six global components, one value
//   sigma_yy_pa    = [...]               # per node, tension positive (spec I1)
//   sigma_zz_pa    = [...]
//   sigma_xy_pa    = [...]
//   sigma_yz_pa    = [...]
//   sigma_xz_pa    = [...]
//   jump_depths_m  = [11951.15]          # optional; each an interior node
//   sigma_xx_pa_below = [...]            # with jumps: per jump, the value just
//   ...                                  # below it (all six components)
//
// Each component is piecewise linear in depth = -(x . up) (up from
// [fault_geometry], the free surface through the origin).  At a jump node
// the plain list holds the value from above and the _below list the value
// from below; a point exactly at a jump depth takes the value below.
// Depths outside [depth_m first, last] abort (no extrapolation; a 1e-9
// relative round-off margin is clamped to the end values).
//
// The tensor is a global Cauchy tensor, so the function has the signature of
// NitscheFaultSetup::prestress and the operator projects it face by face as
// today; [[stress_region]] boxes still replace it inside each box.  The
// projection is an L2 projection per face, so a face that straddles a jump
// would see an overshoot, not the spec's value: the jump depth must be a
// mesh line on the fault, and the driver refuses a mesh with a fault face
// across a jump (CountFaultFacesAcrossDepths, on the serial mesh).
//
// stress_measure = "effective" declares sigma - P_f I (TPV12's Method 1).
// It is exact only while the prestress enters as a fault traction and the
// bulk is linear elastic; the driver refuses it with any bulk rheology.
// "total" is the total stress (no pore pressure; the driver requires P_p = 0
// either way).

#ifndef MFEM_SEAS_NITSCHE_OPT_DEPTH_PRESTRESS_HPP
#define MFEM_SEAS_NITSCHE_OPT_DEPTH_PRESTRESS_HPP

#include "mfem.hpp"
#include "nitsche_opt_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace nitsche_opt
{

struct DepthPrestressTable
{
   enum class Measure { Effective, Total };
   /// Component order of the deck keys (as [[stress_region]]).
   static constexpr const char *kKeys[6] = {"sigma_xx_pa", "sigma_yy_pa", "sigma_zz_pa",
                                            "sigma_xy_pa", "sigma_yz_pa", "sigma_xz_pa"};
   static constexpr const char *kTable = "prestress_depth_table";

   Measure measure = Measure::Total;
   std::vector<real_t> depth;                ///< nodes (m), strictly increasing
   std::vector<real_t> jump;                 ///< jump depths, interior nodes
   std::array<std::vector<real_t>, 6> node;  ///< per component, value per node (from above)
   std::array<std::vector<real_t>, 6> below; ///< per component, value below each jump
   real_t up[3] = {0.0, 0.0, 1.0};           ///< unit vector; depth = -(x . up)

   /// Empty when the table is consistent, else the reason.
   std::string Check() const
   {
      std::ostringstream e;
      const size_t n = depth.size();
      if (n < 2) { e << "depth_m needs at least two nodes"; return e.str(); }
      for (size_t k = 0; k < n; k++)
      {
         if (!std::isfinite(depth[k])) { e << "depth_m[" << k << "] is not finite"; return e.str(); }
         if (k > 0 && !(depth[k] > depth[k - 1]))
         {
            e << "depth_m must be strictly increasing (depth_m[" << k - 1 << "] = "
              << depth[k - 1] << ", depth_m[" << k << "] = " << depth[k]
              << "); list a discontinuity in jump_depths_m instead";
            return e.str();
         }
      }
      for (int c = 0; c < 6; c++)
      {
         if (node[c].size() != n)
         {
            e << kKeys[c] << " has " << node[c].size() << " values; depth_m has " << n;
            return e.str();
         }
         if (below[c].size() != jump.size())
         {
            e << kKeys[c] << "_below has " << below[c].size() << " values; jump_depths_m has "
              << jump.size();
            return e.str();
         }
      }
      for (size_t i = 0; i < jump.size(); i++)
      {
         if (i > 0 && !(jump[i] > jump[i - 1]))
         {
            e << "jump_depths_m must be strictly increasing";
            return e.str();
         }
         if (NodeOf(jump[i]) < 0)
         {
            e << "jump_depths_m[" << i << "] = " << jump[i] << " is not an interior node of "
              "depth_m (it must equal one of depth_m[1] .. depth_m[" << n - 2 << "] exactly)";
            return e.str();
         }
      }
      const real_t u = std::sqrt(up[0] * up[0] + up[1] * up[1] + up[2] * up[2]);
      if (!(std::abs(u - 1.0) <= 1e-12)) { e << "up is not a unit vector (|up| = " << u << ")"; }
      return e.str();
   }

   /// Index of the interior node at depth d (exact match), or -1.
   int NodeOf(real_t d) const
   {
      for (size_t k = 1; k + 1 < depth.size(); k++)
      {
         if (depth[k] == d) { return static_cast<int>(k); }
      }
      return -1;
   }

   real_t Depth(const real_t x[3]) const
   {
      return -(x[0] * up[0] + x[1] * up[1] + x[2] * up[2]);
   }

   /// sigma(x), row-major 3x3, tension positive.
   void Eval(const real_t x[3], real_t s[9]) const
   {
      const size_t n = depth.size();
      // round-off tolerance at the two ends (e.g. a point on the free surface)
      const real_t tol = 1e-9 * (depth[n - 1] - depth[0]);
      real_t d = Depth(x);
      MFEM_VERIFY(d >= depth[0] - tol && d <= depth[n - 1] + tol, "[" << kTable
                  << "]: the point (" << x[0] << ", " << x[1] << ", " << x[2]
                  << ") at depth " << d << " m lies outside depth_m [" << depth[0] << ", "
                  << depth[n - 1] << "] (no extrapolation)");
      d = std::min(std::max(d, depth[0]), depth[n - 1]);
      // interval [k, k+1] with depth[k] <= d < depth[k+1] (k = n-2 at the last node)
      size_t k = static_cast<size_t>(std::upper_bound(depth.begin(), depth.end(), d)
                                     - depth.begin());
      k = (k == 0) ? 0 : k - 1;
      if (k > n - 2) { k = n - 2; }
      const real_t f = (d - depth[k]) / (depth[k + 1] - depth[k]);
      int jk = -1;   // jump at node k: the interval starts from the value below
      for (size_t i = 0; i < jump.size(); i++) { if (jump[i] == depth[k]) { jk = static_cast<int>(i); } }
      real_t v[6];
      for (int c = 0; c < 6; c++)
      {
         const real_t a = (jk >= 0) ? below[c][jk] : node[c][k];
         v[c] = a + (node[c][k + 1] - a) * f;
      }
      const real_t full[9] = {v[0], v[3], v[5], v[3], v[1], v[4], v[5], v[4], v[2]};
      for (int i = 0; i < 9; i++) { s[i] = full[i]; }
   }
};

/// Read [prestress_depth_table] from `t`; `up` from [fault_geometry].
/// Aborts on any unknown key, missing key or inconsistency.
inline DepthPrestressTable ParseDepthPrestress(const toml::value &t, const real_t up[3])
{
   using T = DepthPrestressTable;
   OptRefuseUnknownKeys(t, T::kTable,
                        {"stress_measure", "depth_m", "jump_depths_m",
                         "sigma_xx_pa", "sigma_yy_pa", "sigma_zz_pa",
                         "sigma_xy_pa", "sigma_yz_pa", "sigma_xz_pa",
                         "sigma_xx_pa_below", "sigma_yy_pa_below", "sigma_zz_pa_below",
                         "sigma_xy_pa_below", "sigma_yz_pa_below", "sigma_xz_pa_below"});
   T tab;
   const std::string m = OptStr(t, T::kTable, "stress_measure");
   MFEM_VERIFY(m == "effective" || m == "total", "[" << T::kTable << "] stress_measure = '"
               << m << "': expected \"effective\" or \"total\"");
   tab.measure = (m == "effective") ? T::Measure::Effective : T::Measure::Total;
   tab.depth = OptRealArray(t, T::kTable, "depth_m");
   const bool has_jumps = t.contains("jump_depths_m");
   if (has_jumps) { tab.jump = OptRealArray(t, T::kTable, "jump_depths_m"); }
   for (int c = 0; c < 6; c++)
   {
      tab.node[c] = OptRealArray(t, T::kTable, T::kKeys[c]);
      const std::string kb = std::string(T::kKeys[c]) + "_below";
      if (!tab.jump.empty()) { tab.below[c] = OptRealArray(t, T::kTable, kb.c_str()); }
      else
      {
         MFEM_VERIFY(!t.contains(kb), "[" << T::kTable << "] " << kb
                     << " without jump_depths_m");
      }
   }
   const real_t u = std::sqrt(up[0] * up[0] + up[1] * up[1] + up[2] * up[2]);
   MFEM_VERIFY(u > 0.0, "[" << T::kTable << "]: [fault_geometry] up is zero");
   for (int i = 0; i < 3; i++) { tab.up[i] = up[i] / u; }
   const std::string err = tab.Check();
   MFEM_VERIFY(err.empty(), "[" << T::kTable << "]: " << err);
   return tab;
}

/// Fault faces (boundary elements with attribute fault_attr) whose vertices
/// lie on both sides of a depth in `depths`, more than `tol` m away.  Run on
/// the serial mesh, so the count does not depend on the partition.
inline long long CountFaultFacesAcrossDepths(const Mesh &mesh, int fault_attr,
                                             const std::vector<real_t> &depths,
                                             const real_t up[3], real_t tol,
                                             std::string *first = nullptr)
{
   long long count = 0;
   Array<int> v;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      if (mesh.GetBdrAttribute(be) != fault_attr) { continue; }
      mesh.GetBdrElementVertices(be, v);
      real_t lo = std::numeric_limits<real_t>::max(), hi = -lo;
      for (int k = 0; k < v.Size(); k++)
      {
         const real_t *x = mesh.GetVertex(v[k]);
         const real_t d = -(x[0] * up[0] + x[1] * up[1] + x[2] * up[2]);
         lo = std::min(lo, d);
         hi = std::max(hi, d);
      }
      for (real_t dj : depths)
      {
         if (lo < dj - tol && hi > dj + tol)
         {
            if (count == 0 && first)
            {
               std::ostringstream s;
               s << "fault boundary element " << be << " spans depths [" << lo << ", " << hi
                 << "] m across the jump at " << dj << " m";
               *first = s.str();
            }
            count++;
            break;
         }
      }
   }
   return count;
}

} // namespace nitsche_opt
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_OPT_DEPTH_PRESTRESS_HPP
