// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_chains.hpp — the fault lines of a 2-D mesh as chains, and the
// velocity-strengthening tip-zone field on them (Chunhui Zhao, 2026-10-29:
// "let's write this in the document, and prepare to add this function";
// UCERF 2-D report rev 15, section 3.1).
//
// FaultChains gathers the fault faces of every rank (the boundary elements of
// the fault attribute, segments in 2-D), joins them into chains at their
// shared vertices, and places a point of the fault on its chain: the chain,
// the arc length s from the chain's first end, and the distance d along the
// chain to the nearer free tip.
//   - Vertices are matched by their coordinates, which ParMesh copies bit for
//     bit from the serial mesh, so every rank and every rank count builds the
//     same chains (as FaultOrientation's vertex keys, fault_orientation.hpp).
//   - A vertex on three or more fault faces (a junction, a branch of the
//     fault line) and a closed chain are refused.
//   - A chain end is a free tip unless it lies on the outer boundary (a face
//     with one element) or on a boundary element of another attribute (a
//     creep or prescribed-slip surface): a fault that reaches the free
//     surface (BP3) has no tip there.  A chain without a free tip is refused.
//   - Each chain runs from its end with the smaller (x, y) (lexicographic) and
//     s is summed in that order, so the numbers do not depend on the
//     partition.
//   - A point goes to the nearest fault segment (a sparse grid of cells as
//     large as the longest segment) and must lie on it to 1e-6 of the
//     segment's length; a point off every fault segment is refused.  The
//     cycle's meshes have straight faces (LoadScaledMesh refuses a nodal grid
//     function), so the fault DOFs lie on them.
//
// TipZoneField, deck [cycle.fields.a] kind = "tip_zone":
//   L      the length of the point's chain
//   W_eff  = min(width, max_fraction L)   (width when max_fraction is not set)
//   value  = tip                                          d <= W_eff - taper
//          = tip + (inside - tip) (d - W_eff + taper) / taper   in between
//          = inside                                       d >= W_eff
//   (taper = 0: tip for d <= W_eff, else inside).  W_eff < taper on a chain is
//   refused (CheckTipZone): the tip itself would not get the tip value.
//
// 2-D meshes only: in 3-D the fault is a surface and its edge a curve.

#ifndef MFEM_SEAS_CYCLE_FAULT_CHAINS_HPP
#define MFEM_SEAS_CYCLE_FAULT_CHAINS_HPP

#include "mfem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// The tip-zone field of [cycle.fields.<name>] kind = "tip_zone".
struct TipZoneSpec
{
   real_t width = 0.0;          ///< W [m], measured along the chain from the tip
   real_t taper = 0.0;          ///< the linear ramp over [W - taper, W] [m]
   real_t inside = 0.0;         ///< value from W on (away from the tips)
   real_t tip = 0.0;            ///< value within W - taper of a free tip
   real_t max_fraction = 0.0;   ///< > 0: W_eff = min(W, max_fraction L); 0 = no cap
   bool set = false;            ///< the deck gave this field

   /// Refuses width <= 0, a taper outside [0, width], max_fraction < 0.
   void Validate(const std::string &name) const
   {
      MFEM_VERIFY(width > 0.0, "[cycle.fields." << name << "] width must be > 0");
      MFEM_VERIFY(taper >= 0.0 && taper <= width,
                  "[cycle.fields." << name << "] taper must be in [0, width]");
      MFEM_VERIFY(max_fraction >= 0.0,
                  "[cycle.fields." << name << "] max_fraction must be > 0");
   }

   /// W_eff on a chain of length L.
   real_t Width(real_t L) const
   {
      return (max_fraction > 0.0) ? std::min(width, max_fraction * L) : width;
   }

   /// The value at distance d from the nearer free tip on a chain of length L.
   real_t Value(real_t d, real_t L) const
   {
      const real_t W = Width(L);
      if (taper == 0.0) { return (d <= W) ? tip : inside; }
      const real_t r = (d - (W - taper)) / taper;
      if (r <= 0.0) { return tip; }
      if (r >= 1.0) { return inside; }
      return tip + r * (inside - tip);
   }
};

#ifdef MFEM_USE_MPI

class FaultChains
{
public:
   using XY = std::array<real_t, 2>;

   /// Collective on pmesh's communicator.
   FaultChains(const ParMesh &pmesh, int fault_attr)
   {
      MPI_Comm comm = pmesh.GetComm();
      MFEM_VERIFY(pmesh.Dimension() == 2 && pmesh.SpaceDimension() == 2,
                  "FaultChains: needs a 2-D mesh (fault lines with tips); the mesh is "
                  << pmesh.Dimension() << "-D");
      // ---- this rank's fault segments: x0, y0, x1, y1
      std::vector<double> loc;
      for (int be = 0; be < pmesh.GetNBE(); be++)
      {
         if (pmesh.GetBdrAttribute(be) != fault_attr) { continue; }
         MFEM_VERIFY(pmesh.GetBdrElementGeometry(be) == Geometry::SEGMENT,
                     "FaultChains: fault boundary element " << be << " is not a segment");
         Array<int> v;
         pmesh.GetBdrElementVertices(be, v);
         for (int k = 0; k < 2; k++)
         {
            const real_t *c = pmesh.GetVertex(v[k]);
            loc.push_back(c[0]);
            loc.push_back(c[1]);
         }
      }
      int np = 1;
      MPI_Comm_size(comm, &np);
      int nl = static_cast<int>(loc.size());
      std::vector<int> cnt(np), dsp(np, 0);
      MPI_Allgather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
      for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      std::vector<double> all(std::max(1, dsp[np - 1] + cnt[np - 1]));
      MPI_Allgatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(),
                     MPI_DOUBLE, comm);
      const int ntot = dsp[np - 1] + cnt[np - 1];

      // ---- unique segments (a shared face may come from two ranks), vertices
      std::set<std::pair<XY, XY>> uniq;
      for (int o = 0; o + 4 <= ntot; o += 4)
      {
         XY a = {{all[o], all[o + 1]}}, b = {{all[o + 2], all[o + 3]}};
         MFEM_VERIFY(a != b, "FaultChains: a fault segment of zero length at ("
                     << a[0] << ", " << a[1] << ")");
         if (b < a) { std::swap(a, b); }
         uniq.insert({a, b});
      }
      MFEM_VERIFY(!uniq.empty(), "FaultChains: the mesh has no fault faces (attribute "
                  << fault_attr << ")");
      std::map<XY, std::vector<int>> vseg;   // vertex -> its segments (lexicographic)
      std::vector<std::pair<XY, XY>> seg(uniq.begin(), uniq.end());
      for (int k = 0; k < static_cast<int>(seg.size()); k++)
      {
         vseg[seg[k].first].push_back(k);
         vseg[seg[k].second].push_back(k);
      }
      int junctions = 0;
      XY first_junction = {{0.0, 0.0}};
      for (const auto &v : vseg)
      {
         if (v.second.size() >= 3)
         {
            if (junctions == 0) { first_junction = v.first; }
            junctions++;
         }
      }
      MFEM_VERIFY(junctions == 0, "FaultChains: " << junctions << " fault vertices carry "
                  "three or more fault faces (a junction, first at (" << first_junction[0]
                  << ", " << first_junction[1] << ")); the tip distance needs separate "
                  "chains");

      // ---- chains, walked from their lexicographically smaller end
      std::vector<char> used(seg.size(), 0);
      for (const auto &v : vseg)
      {
         if (v.second.size() != 1 || used[v.second[0]]) { continue; }
         Chain ch;
         ch.first = static_cast<int>(pts_.size());
         XY p = v.first;
         int k = v.second[0];
         pts_.push_back(p);
         s_.push_back(0.0);
         while (true)
         {
            used[k] = 1;
            const XY q = (seg[k].first == p) ? seg[k].second : seg[k].first;
            const real_t dx = q[0] - p[0], dy = q[1] - p[1];
            s_.push_back(s_.back() + std::sqrt(dx * dx + dy * dy));
            pts_.push_back(q);
            chain_of_seg_.push_back(static_cast<int>(chains_.size()));
            p = q;
            const std::vector<int> &next = vseg[p];
            if (next.size() == 1) { break; }   // the other end
            k = (next[0] == k) ? next[1] : next[0];
         }
         ch.last = static_cast<int>(pts_.size()) - 1;
         ch.length = s_[ch.last];
         chains_.push_back(ch);
      }
      int loose = 0;
      for (char u : used) { loose += u ? 0 : 1; }
      MFEM_VERIFY(loose == 0, "FaultChains: " << loose << " fault segments form closed "
                  "chains (no tip)");
      // pts_ holds each chain's vertices in order; segment i of the walk runs
      // pts_[i + c] -> pts_[i + c + 1] with c the number of earlier chains.

      // ---- free tips: an end on the outer boundary or on another boundary
      //      attribute is not one
      std::set<XY> other;
      for (int f = 0; f < pmesh.GetNumFaces(); f++)
      {
         if (!pmesh.GetFaceInformation(f).IsBoundary()) { continue; }
         Array<int> v;
         pmesh.GetFaceVertices(f, v);
         for (int k = 0; k < v.Size(); k++)
         {
            other.insert({{pmesh.GetVertex(v[k])[0], pmesh.GetVertex(v[k])[1]}});
         }
      }
      for (int be = 0; be < pmesh.GetNBE(); be++)
      {
         if (pmesh.GetBdrAttribute(be) == fault_attr) { continue; }
         Array<int> v;
         pmesh.GetBdrElementVertices(be, v);
         for (int k = 0; k < v.Size(); k++)
         {
            other.insert({{pmesh.GetVertex(v[k])[0], pmesh.GetVertex(v[k])[1]}});
         }
      }
      const int nc = static_cast<int>(chains_.size());
      std::vector<int> bound(2 * nc, 0), gbound(2 * nc, 0);
      for (int c = 0; c < nc; c++)
      {
         bound[2 * c] = other.count(pts_[chains_[c].first]) ? 1 : 0;
         bound[2 * c + 1] = other.count(pts_[chains_[c].last]) ? 1 : 0;
      }
      MPI_Allreduce(bound.data(), gbound.data(), 2 * nc, MPI_INT, MPI_MAX, comm);
      for (int c = 0; c < nc; c++)
      {
         Chain &ch = chains_[c];
         ch.free_first = (gbound[2 * c] == 0);
         ch.free_last = (gbound[2 * c + 1] == 0);
         const XY &a = pts_[ch.first], &b = pts_[ch.last];
         MFEM_VERIFY(ch.free_first || ch.free_last,
                     "FaultChains: chain " << c << " from (" << a[0] << ", " << a[1] << ") to ("
                     << b[0] << ", " << b[1] << ") has no free tip (both ends lie on the outer "
                     "boundary or on another boundary attribute)");
      }

      // ---- the sparse grid of segments
      for (const XY &p : pts_)
      {
         lo_[0] = std::min(lo_[0], p[0]);
         lo_[1] = std::min(lo_[1], p[1]);
      }
      h_ = 0.0;
      ForEachSegment_([&](int, int i0) { h_ = std::max(h_, s_[i0 + 1] - s_[i0]); });
      ForEachSegment_([&](int g, int i0)
      {
         const XY &a = pts_[i0], &b = pts_[i0 + 1];
         const long long x0 = Cell_(std::min(a[0], b[0]), 0), x1 = Cell_(std::max(a[0], b[0]), 0);
         const long long y0 = Cell_(std::min(a[1], b[1]), 1), y1 = Cell_(std::max(a[1], b[1]), 1);
         for (long long ix = x0; ix <= x1; ix++)
            for (long long iy = y0; iy <= y1; iy++) { grid_[{ix, iy}].push_back(g); }
      });
   }

   int NumChains() const { return static_cast<int>(chains_.size()); }
   real_t Length(int c) const { return chains_[c].length; }
   /// The chain's ends in walk order (the first has the smaller (x, y)).
   const XY &FirstEnd(int c) const { return pts_[chains_[c].first]; }
   const XY &LastEnd(int c) const { return pts_[chains_[c].last]; }
   bool FreeFirst(int c) const { return chains_[c].free_first; }
   bool FreeLast(int c) const { return chains_[c].free_last; }
   int NumSegments() const { return static_cast<int>(chain_of_seg_.size()); }

   /// The distance along a chain of length L at arc length s to its nearer
   /// free tip.
   real_t TipDistance(int c, real_t s) const
   {
      const Chain &ch = chains_[c];
      const real_t d0 = ch.free_first ? s : std::numeric_limits<real_t>::max();
      const real_t d1 = ch.free_last ? ch.length - s : std::numeric_limits<real_t>::max();
      return std::min(d0, d1);
   }

   /// Places the point x (z = 0) on its fault segment: the chain, the arc
   /// length s from the chain's first end, the distance d to the nearer free
   /// tip.  Refuses a point off every fault segment.
   void Locate(const real_t x[3], int &chain, real_t &s, real_t &d) const
   {
      MFEM_VERIFY(x[2] == 0.0, "FaultChains::Locate: a 2-D point has z = 0, got " << x[2]);
      const long long cx = Cell_(x[0], 0), cy = Cell_(x[1], 1);
      int best = -1, best_i0 = -1;
      real_t best_dist = std::numeric_limits<real_t>::max(), best_t = 0.0;
      for (long long ix = cx - 1; ix <= cx + 1; ix++)
         for (long long iy = cy - 1; iy <= cy + 1; iy++)
         {
            auto it = grid_.find({ix, iy});
            if (it == grid_.end()) { continue; }
            for (int g : it->second)
            {
               const int i0 = g + chain_of_seg_[g];
               const XY &a = pts_[i0], &b = pts_[i0 + 1];
               const real_t ex = b[0] - a[0], ey = b[1] - a[1];
               const real_t px = x[0] - a[0], py = x[1] - a[1];
               real_t t = (px * ex + py * ey) / (ex * ex + ey * ey);
               t = std::min(std::max(t, real_t(0.0)), real_t(1.0));
               const real_t qx = px - t * ex, qy = py - t * ey;
               const real_t dist = std::sqrt(qx * qx + qy * qy);
               if (dist < best_dist || (dist == best_dist && g < best))
               {
                  best = g;
                  best_i0 = i0;
                  best_dist = dist;
                  best_t = t;
               }
            }
         }
      const real_t len = (best < 0) ? 0.0 : s_[best_i0 + 1] - s_[best_i0];
      MFEM_VERIFY(best >= 0 && best_dist <= 1e-6 * len,
                  "FaultChains::Locate: the point (" << x[0] << ", " << x[1] << ") lies on no "
                  "fault segment (" << ((best < 0) ? std::string("none within a cell") :
                  "nearest " + std::to_string(best_dist) + " m away") << ")");
      chain = chain_of_seg_[best];
      // the segment's ends give s exactly: t = 1 is s_[i0 + 1]
      s = (best_t == 1.0) ? s_[best_i0 + 1] : s_[best_i0] + best_t * len;
      d = TipDistance(chain, s);
   }

   /// One line per chain: length, ends, free tips (and W_eff, the tip-zone
   /// length with spec).
   std::string Report(const TipZoneSpec *spec = nullptr) const
   {
      std::ostringstream o;
      int free = 0;
      for (const Chain &ch : chains_) { free += (ch.free_first ? 1 : 0) + (ch.free_last ? 1 : 0); }
      o << NumChains() << " chains, " << NumSegments() << " segments, " << free
        << " free tips\n";
      for (int c = 0; c < NumChains(); c++)
      {
         const Chain &ch = chains_[c];
         const XY &a = pts_[ch.first], &b = pts_[ch.last];
         o << std::setprecision(6) << "    chain " << c << ": L = " << ch.length / 1e3
           << " km, (" << a[0] / 1e3 << ", " << a[1] / 1e3 << ") km" << (ch.free_first ? " tip" : " bound")
           << " to (" << b[0] / 1e3 << ", " << b[1] / 1e3 << ") km" << (ch.free_last ? " tip" : " bound");
         if (spec)
         {
            const real_t W = spec->Width(ch.length);
            const int nt = (ch.free_first ? 1 : 0) + (ch.free_last ? 1 : 0);
            o << "; W_eff = " << W / 1e3 << " km, a = tip value on "
              << std::min(ch.length, nt * std::max(W - spec->taper, real_t(0.0))) / 1e3
              << " km, weakening (d >= W_eff) on "
              << std::max(ch.length - nt * W, real_t(0.0)) / 1e3 << " km";
         }
         o << "\n";
      }
      return o.str();
   }

private:
   struct Chain
   {
      int first = 0, last = 0;   ///< indices into pts_ / s_
      real_t length = 0.0;
      bool free_first = true, free_last = true;
   };

   std::vector<XY> pts_;            ///< every chain's vertices, chain after chain
   std::vector<real_t> s_;          ///< arc length of pts_ from its chain's first end
   std::vector<int> chain_of_seg_;  ///< walk segment -> chain
   std::vector<Chain> chains_;
   std::map<std::pair<long long, long long>, std::vector<int>> grid_;
   XY lo_ = {{std::numeric_limits<real_t>::max(), std::numeric_limits<real_t>::max()}};
   real_t h_ = 1.0;

   long long Cell_(real_t v, int axis) const
   {
      return static_cast<long long>(std::floor((v - lo_[axis]) / h_));
   }

   /// f(g, i0) for every walk segment g, running pts_[i0] -> pts_[i0 + 1].
   template <typename F> void ForEachSegment_(F f) const
   {
      for (int g = 0; g < static_cast<int>(chain_of_seg_.size()); g++)
      {
         f(g, g + chain_of_seg_[g]);
      }
   }
};

/// The field value of a tip zone at a fault point x.
inline real_t TipZoneField(const TipZoneSpec &f, const FaultChains &chains, const real_t x[3])
{
   int c = -1;
   real_t s = 0.0, d = 0.0;
   chains.Locate(x, c, s, d);
   return f.Value(d, chains.Length(c));
}

/// Refuses a chain on which W_eff < taper (the tip would not get the tip value).
inline void CheckTipZone(const TipZoneSpec &f, const FaultChains &chains, const std::string &name)
{
   for (int c = 0; c < chains.NumChains(); c++)
   {
      const real_t W = f.Width(chains.Length(c));
      MFEM_VERIFY(W >= f.taper, "[cycle.fields." << name << "] kind = \"tip_zone\": on chain "
                  << c << " (L = " << chains.Length(c) << " m) max_fraction L = " << W
                  << " m is shorter than the taper (" << f.taper << " m)");
   }
}

#endif // MFEM_USE_MPI

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FAULT_CHAINS_HPP
