// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_surface_zones.hpp — the velocity-strengthening tip and depth zones of
// [cycle.fields.a] kind = "tip_zone" on 3-D fault surfaces (UCERF 3-D plan
// rev5, P1.2; Chunhui Zhao, 2026-11-25: "start P1").  The 2-D tip zone
// (fault_chains.hpp) is unchanged; this file is its surface version.
//
// The surface.  FaultSurfaceZones gathers the fault triangles of every rank
// (the boundary elements of the fault attribute), de-duplicated by their
// vertex coordinates (ParMesh copies them bit for bit from the serial mesh,
// so every rank and rank count builds the same surface), and splits them
// into patches (triangles joined through shared edges).  An edge with three
// or more fault triangles (a junction) is refused.
//
// Lateral tip edges.  A boundary edge of the surface (one fault triangle) is
// free unless both its vertices lie on the outer boundary or on a boundary
// element of another attribute (the free surface z = 0, a creep or
// prescribed-slip surface): the 2-D rule for chain ends, applied to an edge.
// A free edge is lateral when it descends at least 30 degrees from the
// horizontal (|dz| >= L / 2): the ends of the fault along strike.  A free
// edge flatter than that (the bottom of a fault with no creep surface below,
// loading option B) is not a tip; the depth zone below handles the bottom.
//
// The tip distance d of a point x is the straight-line (3-D) distance to the
// nearest lateral edge segment of x's own patch (decision card 11-25,
// "Straight-line to the tip edges"; exact on planar faults; on a fault
// curved with radius R it is short of the along-surface distance by about
// d^3 / (24 R^2), 33 m at d = 20 km and R = 70 km).  A patch with no lateral
// edge has d = infinity.  x is placed on its patch through the nearest fault
// triangle (a sparse grid of cells as large as the longest edge) and must lie
// on it to 1e-6 of that size; a point off every triangle is refused.
//
// The value.  With depth z_d = -z (z = 0 at the free surface):
//    w_tip = 1 for d <= W - taper, linear to 0 at d = W        (W = width)
//    w_top = 1 for z_d <= top_depth, linear to 0 at top_depth + top_taper
//    w_bot = 0 for z_d <= D(x), linear to 1 at D(x) + bottom_taper
//    a     = inside + max(w_tip, w_top, w_bot) (tip - inside)
// (a taper of 0 is a step: 1 at and inside the zone).  The depth zones are
// optional: top_depth and lower_depth / lower_depth_file each turn one on.
// D(x) is the lower seismogenic depth: lower_depth (constant), or from
// lower_depth_file, a CSV with the header "chain,x_m,y_m,lower_depth_m" and
// one row per trace vertex, the rows of a chain in order along it.  Each
// patch takes the chain whose trace is nearest (mean horizontal distance) to
// the patch's shallowest vertices; D(x) is the lower depth interpolated
// linearly along that trace at the horizontal projection of x (clamped at the
// trace ends).  The 2 km tapers between sections are rows of the file (the
// script that writes it puts a vertex 1 km either side of each boundary).
//
// 3-D meshes only (FaultChains is the 2-D version).

#ifndef MFEM_SEAS_CYCLE_FAULT_SURFACE_ZONES_HPP
#define MFEM_SEAS_CYCLE_FAULT_SURFACE_ZONES_HPP

#include "mfem.hpp"
#include "fault_chains.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
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

#ifdef MFEM_USE_MPI

class FaultSurfaceZones
{
public:
   using P3 = std::array<real_t, 3>;

   /// Collective on pmesh's communicator.
   FaultSurfaceZones(const ParMesh &pmesh, int fault_attr, const TipZoneSpec &spec)
      : spec_(spec)
   {
      MPI_Comm comm = pmesh.GetComm();
      MFEM_VERIFY(pmesh.Dimension() == 3 && pmesh.SpaceDimension() == 3,
                  "FaultSurfaceZones: needs a 3-D mesh; the mesh is " << pmesh.Dimension()
                  << "-D");
      // ---- every rank's fault triangles (9 coordinates each)
      std::vector<double> loc;
      for (int be = 0; be < pmesh.GetNBE(); be++)
      {
         if (pmesh.GetBdrAttribute(be) != fault_attr) { continue; }
         MFEM_VERIFY(pmesh.GetBdrElementGeometry(be) == Geometry::TRIANGLE,
                     "FaultSurfaceZones: fault boundary element " << be << " is not a triangle");
         Array<int> v;
         pmesh.GetBdrElementVertices(be, v);
         for (int k = 0; k < 3; k++)
            for (int d = 0; d < 3; d++) { loc.push_back(pmesh.GetVertex(v[k])[d]); }
      }
      std::vector<double> all = Gather_(loc, comm);
      std::set<std::array<P3, 3>> uniq;
      for (size_t o = 0; o + 9 <= all.size(); o += 9)
      {
         std::array<P3, 3> t;
         for (int k = 0; k < 3; k++) { t[k] = {all[o + 3 * k], all[o + 3 * k + 1], all[o + 3 * k + 2]}; }
         std::sort(t.begin(), t.end());
         uniq.insert(t);
      }
      MFEM_VERIFY(!uniq.empty(), "FaultSurfaceZones: the mesh has no fault faces (attribute "
                  << fault_attr << ")");
      std::map<P3, int> vid;
      for (const auto &t : uniq)
         for (const P3 &p : t) { vid.emplace(p, 0); }
      int nv = 0;
      for (auto &kv : vid) { kv.second = nv++; xv_.push_back(kv.first); }
      for (const auto &t : uniq) { tri_.push_back({vid.at(t[0]), vid.at(t[1]), vid.at(t[2])}); }
      const int nt = static_cast<int>(tri_.size());

      // ---- edges, junctions, patches
      std::map<std::pair<int, int>, std::vector<int>> edges;
      for (int t = 0; t < nt; t++)
         for (int a = 0; a < 3; a++)
         {
            int i = tri_[t][a], j = tri_[t][(a + 1) % 3];
            if (j < i) { std::swap(i, j); }
            edges[{i, j}].push_back(t);
         }
      std::vector<int> comp(nt);
      for (int t = 0; t < nt; t++) { comp[t] = t; }
      auto root = [&](int t) { while (comp[t] != t) { t = comp[t] = comp[comp[t]]; } return t; };
      int junctions = 0;
      for (const auto &e : edges)
      {
         if (e.second.size() >= 3) { junctions++; continue; }
         if (e.second.size() == 2)
         {
            const int a = root(e.second[0]), b = root(e.second[1]);
            if (a != b) { comp[std::max(a, b)] = std::min(a, b); }
         }
      }
      MFEM_VERIFY(junctions == 0, "FaultSurfaceZones: " << junctions << " fault edges carry "
                  "three or more fault faces (a junction); the tip distance needs separate "
                  "surfaces");
      std::map<int, int> pid;
      patch_.resize(nt);
      for (int t = 0; t < nt; t++)
      {
         const int r = root(t);
         auto it = pid.find(r);
         if (it == pid.end()) { it = pid.emplace(r, static_cast<int>(pid.size())).first; }
         patch_[t] = it->second;
      }
      np_ = static_cast<int>(pid.size());

      // ---- vertices on the outer boundary or another attribute (any rank)
      std::vector<int> other(nv, 0), gother(nv, 0);
      auto mark = [&](const Array<int> &v)
      {
         for (int k = 0; k < v.Size(); k++)
         {
            const real_t *c = pmesh.GetVertex(v[k]);
            auto it = vid.find({c[0], c[1], c[2]});
            if (it != vid.end()) { other[it->second] = 1; }
         }
      };
      // A face between ranks reads IsBoundary() until ExchangeFaceNbrData
      // (measured on T-plane, 4 ranks: 877 of rank 0's 1502), so skip those.
      std::vector<char> shared(pmesh.GetNumFaces(), 0);
      for (int s = 0; s < pmesh.GetNSharedFaces(); s++) { shared[pmesh.GetSharedFace(s)] = 1; }
      for (int f = 0; f < pmesh.GetNumFaces(); f++)
      {
         if (shared[f] || !pmesh.GetFaceInformation(f).IsBoundary()) { continue; }
         Array<int> v;
         pmesh.GetFaceVertices(f, v);
         mark(v);
      }
      for (int be = 0; be < pmesh.GetNBE(); be++)
      {
         if (pmesh.GetBdrAttribute(be) == fault_attr) { continue; }
         Array<int> v;
         pmesh.GetBdrElementVertices(be, v);
         mark(v);
      }
      MPI_Allreduce(other.data(), gother.data(), nv, MPI_INT, MPI_MAX, comm);

      // ---- lateral free edges, per patch
      lat_.assign(np_, {});
      for (const auto &e : edges)
      {
         if (e.second.size() != 1) { continue; }
         const int i = e.first.first, j = e.first.second;
         if (gother[i] && gother[j]) { bound_edges_++; continue; }
         const P3 &a = xv_[i], &b = xv_[j];
         const real_t L = Dist_(a, b);
         if (std::abs(b[2] - a[2]) >= 0.5 * L) { lat_[patch_[e.second[0]]].push_back({a, b}); }
         else { flat_free_edges_++; }
      }

      // ---- the sparse grid of triangles
      h_ = 0.0;
      for (const auto &e : edges) { h_ = std::max(h_, Dist_(xv_[e.first.first], xv_[e.first.second])); }
      for (const P3 &p : xv_)
         for (int d = 0; d < 3; d++) { lo_[d] = std::min(lo_[d], p[d]); }
      for (int t = 0; t < nt; t++)
      {
         std::array<long long, 3> c0, c1;
         for (int d = 0; d < 3; d++)
         {
            real_t mn = xv_[tri_[t][0]][d], mx = mn;
            for (int a = 1; a < 3; a++) { mn = std::min(mn, xv_[tri_[t][a]][d]); mx = std::max(mx, xv_[tri_[t][a]][d]); }
            c0[d] = Cell_(mn, d);
            c1[d] = Cell_(mx, d);
         }
         for (long long i = c0[0]; i <= c1[0]; i++)
            for (long long j = c0[1]; j <= c1[1]; j++)
               for (long long k = c0[2]; k <= c1[2]; k++) { grid_[{i, j, k}].push_back(t); }
      }

      // ---- the lower depth
      if (!spec_.lower_depth_file.empty()) { ReadLowerDepth_(); }
   }

   int NumPatches() const { return np_; }
   int NumTriangles() const { return static_cast<int>(tri_.size()); }
   int NumLateralEdges(int p) const { return static_cast<int>(lat_[p].size()); }

   /// The patch of fault point x and its straight-line distance to the
   /// nearest lateral edge of that patch (infinity without one).
   void Locate(const real_t x[3], int &patch, real_t &d) const
   {
      const int t = NearestTriangle_(x);
      patch = patch_[t];
      d = std::numeric_limits<real_t>::infinity();
      const P3 p = {x[0], x[1], x[2]};
      for (const auto &s : lat_[patch]) { d = std::min(d, SegDist_(p, s[0], s[1])); }
   }

   /// D(x), the lower seismogenic depth at fault point x (positive, m);
   /// -1 when the deck gives none.
   real_t LowerDepth(const real_t x[3], int patch) const
   {
      if (spec_.lower_depth > 0.0) { return spec_.lower_depth; }
      if (chains_.empty()) { return -1.0; }
      return OnChain_(chain_of_patch_[patch], x);
   }

   /// The a value at fault point x (file comment).
   real_t Value(const real_t x[3]) const
   {
      int p = 0;
      real_t d = 0.0;
      Locate(x, p, d);
      real_t w = Ramp_(spec_.width - spec_.taper, spec_.width, d);
      const real_t zd = -x[2];
      if (spec_.top_depth >= 0.0)
      {
         w = std::max(w, Ramp_(spec_.top_depth, spec_.top_depth + spec_.top_taper, zd));
      }
      const real_t D = LowerDepth(x, p);
      if (D > 0.0) { w = std::max(w, 1.0 - Ramp_(D, D + spec_.bottom_taper, zd)); }
      return spec_.inside + w * (spec_.tip - spec_.inside);
   }

   /// One line per patch: triangles, lateral edges and their extent, the chain.
   std::string Report() const
   {
      std::ostringstream o;
      int nl = 0;
      for (const auto &l : lat_) { nl += static_cast<int>(l.size()); }
      o << np_ << " fault surface(s), " << tri_.size() << " triangles, " << nl
        << " lateral tip edges, " << flat_free_edges_ << " flat free edges (not tips), "
        << bound_edges_ << " edges on the outer boundary or another attribute\n";
      std::vector<int> ntp(np_, 0);
      for (int p : patch_) { ntp[p]++; }
      for (int p = 0; p < np_; p++)
      {
         o << std::setprecision(6) << "    surface " << p << ": " << ntp[p] << " triangles, "
           << lat_[p].size() << " lateral edges";
         if (!chains_.empty()) { o << ", lower depth from chain \"" << chain_names_[chain_of_patch_[p]] << "\""; }
         o << "\n";
      }
      o << "    a = " << spec_.inside << " + max(w_tip, w_top, w_bot) (" << spec_.tip << " - "
        << spec_.inside << "); tip zone W = " << spec_.width << " m, taper " << spec_.taper << " m";
      if (spec_.top_depth >= 0.0)
      {
         o << "; VS above " << spec_.top_depth << " m, taper " << spec_.top_taper << " m";
      }
      if (spec_.lower_depth > 0.0) { o << "; VS below " << spec_.lower_depth << " m"; }
      else if (!chains_.empty()) { o << "; VS below the lower depth of " << spec_.lower_depth_file; }
      if (spec_.lower_depth > 0.0 || !chains_.empty()) { o << ", taper " << spec_.bottom_taper << " m"; }
      o << "\n";
      return o.str();
   }

private:
   static std::vector<double> Gather_(const std::vector<double> &loc, MPI_Comm comm)
   {
      int np = 1;
      MPI_Comm_size(comm, &np);
      int nl = static_cast<int>(loc.size());
      std::vector<int> cnt(np), dsp(np, 0);
      MPI_Allgather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
      for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      std::vector<double> all(std::max(1, dsp[np - 1] + cnt[np - 1]));
      MPI_Allgatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(),
                     MPI_DOUBLE, comm);
      all.resize(dsp[np - 1] + cnt[np - 1]);
      return all;
   }

   static real_t Dist_(const P3 &a, const P3 &b)
   {
      return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
                       + (a[2] - b[2]) * (a[2] - b[2]));
   }

   static real_t SegDist_(const P3 &p, const P3 &a, const P3 &b)
   {
      real_t e[3], q[3], ee = 0.0, qe = 0.0;
      for (int d = 0; d < 3; d++)
      {
         e[d] = b[d] - a[d];
         q[d] = p[d] - a[d];
         ee += e[d] * e[d];
         qe += q[d] * e[d];
      }
      const real_t t = std::min(std::max(qe / ee, real_t(0.0)), real_t(1.0));
      real_t s = 0.0;
      for (int d = 0; d < 3; d++) { s += (q[d] - t * e[d]) * (q[d] - t * e[d]); }
      return std::sqrt(s);
   }

   /// 1 for v <= lo, linear to 0 at hi (a step at lo when hi == lo).
   static real_t Ramp_(real_t lo, real_t hi, real_t v)
   {
      if (v <= lo) { return 1.0; }
      if (v >= hi) { return 0.0; }
      return (hi - v) / (hi - lo);
   }

   long long Cell_(real_t v, int d) const
   {
      return static_cast<long long>(std::floor((v - lo_[d]) / h_));
   }

   /// The distance from p to triangle t (Ericson, closest point on a triangle).
   real_t TriDist_(const P3 &p, int t) const
   {
      const P3 &a = xv_[tri_[t][0]], &b = xv_[tri_[t][1]], &c = xv_[tri_[t][2]];
      auto sub = [](const P3 &u, const P3 &v) { return P3{u[0] - v[0], u[1] - v[1], u[2] - v[2]}; };
      auto dot = [](const P3 &u, const P3 &v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
      const P3 ab = sub(b, a), ac = sub(c, a), ap = sub(p, a);
      const real_t d1 = dot(ab, ap), d2 = dot(ac, ap);
      if (d1 <= 0 && d2 <= 0) { return Dist_(p, a); }
      const P3 bp = sub(p, b);
      const real_t d3 = dot(ab, bp), d4 = dot(ac, bp);
      if (d3 >= 0 && d4 <= d3) { return Dist_(p, b); }
      const real_t vc = d1 * d4 - d3 * d2;
      if (vc <= 0 && d1 >= 0 && d3 <= 0) { return SegDist_(p, a, b); }
      const P3 cp = sub(p, c);
      const real_t d5 = dot(ab, cp), d6 = dot(ac, cp);
      if (d6 >= 0 && d5 <= d6) { return Dist_(p, c); }
      const real_t vb = d5 * d2 - d1 * d6;
      if (vb <= 0 && d2 >= 0 && d6 <= 0) { return SegDist_(p, a, c); }
      const real_t va = d3 * d6 - d5 * d4;
      if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) { return SegDist_(p, b, c); }
      const real_t den = 1.0 / (va + vb + vc), v = vb * den, w = vc * den;
      const P3 q = {a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w,
                    a[2] + ab[2] * v + ac[2] * w};
      return Dist_(p, q);
   }

   int NearestTriangle_(const real_t x[3]) const
   {
      const P3 p = {x[0], x[1], x[2]};
      const long long c[3] = {Cell_(x[0], 0), Cell_(x[1], 1), Cell_(x[2], 2)};
      int best = -1;
      real_t bd = std::numeric_limits<real_t>::max();
      for (long long i = c[0] - 1; i <= c[0] + 1; i++)
         for (long long j = c[1] - 1; j <= c[1] + 1; j++)
            for (long long k = c[2] - 1; k <= c[2] + 1; k++)
            {
               auto it = grid_.find({i, j, k});
               if (it == grid_.end()) { continue; }
               for (int t : it->second)
               {
                  const real_t d = TriDist_(p, t);
                  if (d < bd || (d == bd && t < best)) { bd = d; best = t; }
               }
            }
      MFEM_VERIFY(best >= 0 && bd <= 1e-6 * h_, "FaultSurfaceZones: the point (" << x[0]
                  << ", " << x[1] << ", " << x[2] << ") lies on no fault triangle ("
                  << ((best < 0) ? std::string("none within a cell")
                      : "nearest " + std::to_string(bd) + " m away") << ")");
      return best;
   }

   // lower_depth_file: chain,x_m,y_m,lower_depth_m
   void ReadLowerDepth_()
   {
      std::ifstream f(spec_.lower_depth_file);
      MFEM_VERIFY(f.good(), "[cycle.fields.a] lower_depth_file: cannot read \""
                  << spec_.lower_depth_file << "\"");
      std::string line;
      std::getline(f, line);
      MFEM_VERIFY(line.rfind("chain,x_m,y_m,lower_depth_m", 0) == 0,
                  "[cycle.fields.a] lower_depth_file " << spec_.lower_depth_file
                  << ": the header must be \"chain,x_m,y_m,lower_depth_m\", got \"" << line << "\"");
      std::map<std::string, int> cid;
      int row = 1;
      while (std::getline(f, line))
      {
         row++;
         if (line.empty()) { continue; }
         std::stringstream ss(line);
         std::string name, a, b, c;
         MFEM_VERIFY(std::getline(ss, name, ',') && std::getline(ss, a, ',')
                     && std::getline(ss, b, ',') && std::getline(ss, c),
                     "[cycle.fields.a] lower_depth_file row " << row << ": 4 columns");
         auto it = cid.find(name);
         if (it == cid.end())
         {
            it = cid.emplace(name, static_cast<int>(chains_.size())).first;
            chains_.push_back({});
            chain_names_.push_back(name);
         }
         const real_t D = std::stod(c);
         MFEM_VERIFY(D > 0.0, "[cycle.fields.a] lower_depth_file row " << row
                     << ": lower_depth_m must be > 0");
         chains_[it->second].push_back({std::stod(a), std::stod(b), D});
      }
      MFEM_VERIFY(!chains_.empty(), "[cycle.fields.a] lower_depth_file "
                  << spec_.lower_depth_file << " has no rows");
      // each patch: the chain nearest its shallowest vertices
      std::vector<real_t> top(np_, -std::numeric_limits<real_t>::max());
      for (size_t t = 0; t < tri_.size(); t++)
         for (int a = 0; a < 3; a++) { top[patch_[t]] = std::max(top[patch_[t]], xv_[tri_[t][a]][2]); }
      std::vector<std::vector<real_t>> sum(np_, std::vector<real_t>(chains_.size(), 0.0));
      std::vector<std::set<int>> seen(np_);
      for (size_t t = 0; t < tri_.size(); t++)
         for (int a = 0; a < 3; a++)
         {
            const int v = tri_[t][a], p = patch_[t];
            if (xv_[v][2] < top[p] - 1e-6 * h_ || !seen[p].insert(v).second) { continue; }
            for (size_t c = 0; c < chains_.size(); c++)
            {
               real_t s = 0.0, dist = 0.0;
               Project_(static_cast<int>(c), xv_[v].data(), s, dist);
               sum[p][c] += dist;
            }
         }
      chain_of_patch_.assign(np_, 0);
      for (int p = 0; p < np_; p++)
      {
         chain_of_patch_[p] = static_cast<int>(std::min_element(sum[p].begin(), sum[p].end())
                                               - sum[p].begin());
      }
   }

   /// Horizontal projection of x on chain c: the lower depth there (return
   /// value of OnChain_), its horizontal distance (dist).
   real_t Project_(int c, const real_t x[3], real_t &D, real_t &dist) const
   {
      const auto &pts = chains_[c];
      dist = std::numeric_limits<real_t>::max();
      D = pts[0][2];
      for (size_t i = 0; i + 1 < pts.size(); i++)
      {
         const real_t ex = pts[i + 1][0] - pts[i][0], ey = pts[i + 1][1] - pts[i][1];
         const real_t px = x[0] - pts[i][0], py = x[1] - pts[i][1];
         const real_t ee = ex * ex + ey * ey;
         const real_t t = (ee > 0.0) ? std::min(std::max((px * ex + py * ey) / ee, real_t(0.0)),
                                                real_t(1.0)) : 0.0;
         const real_t qx = px - t * ex, qy = py - t * ey, d = std::sqrt(qx * qx + qy * qy);
         if (d < dist) { dist = d; D = pts[i][2] + t * (pts[i + 1][2] - pts[i][2]); }
      }
      if (pts.size() == 1)
      {
         dist = std::hypot(x[0] - pts[0][0], x[1] - pts[0][1]);
      }
      return D;
   }

   real_t OnChain_(int c, const real_t x[3]) const
   {
      real_t D = 0.0, dist = 0.0;
      Project_(c, x, D, dist);
      return D;
   }

   TipZoneSpec spec_;
   std::vector<P3> xv_;                          ///< fault vertices (sorted)
   std::vector<std::array<int, 3>> tri_;         ///< fault triangles
   std::vector<int> patch_;                      ///< triangle -> patch
   int np_ = 0;
   std::vector<std::vector<std::array<P3, 2>>> lat_;   ///< lateral edges per patch
   int bound_edges_ = 0, flat_free_edges_ = 0;
   std::map<std::array<long long, 3>, std::vector<int>> grid_;
   P3 lo_ = {{std::numeric_limits<real_t>::max(), std::numeric_limits<real_t>::max(),
              std::numeric_limits<real_t>::max()}};
   real_t h_ = 1.0;
   std::vector<std::vector<P3>> chains_;         ///< lower-depth traces: x, y, D
   std::vector<std::string> chain_names_;
   std::vector<int> chain_of_patch_;
};

#endif // MFEM_USE_MPI

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FAULT_SURFACE_ZONES_HPP
