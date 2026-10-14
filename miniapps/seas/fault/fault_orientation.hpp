// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_orientation.hpp — Nitsche dynamic rupture, plan item N-08: the
// oriented fault normal field and side labels of the fault coordinate
// specification rev 8.3 (sections 3 and 5.1).
//
//   J1  junction edges: edges carrying m >= 3 fault faces (m = 1: free tip).
//   O1  face adjacency across non-junction shared edges (m = 2).
//   O2  per patch, seed = face with the smallest sorted vertex key.
//   O3  propagation by the shared edge: canonical cycle C_f (vertices in
//       increasing key order for a triangle), tau_f(e) = +1 if C_f runs e
//       from the smaller to the larger key, s_g = -s_f tau_f(e) tau_g(e).
//   O4/O5  one pin per patch (this header: "reference", J6-gated; "point").
//   O6  an inconsistent loop aborts (non-orientable).
//   J3  parity BFS over H_S (face-incident elements + fans of interior
//       edges; flip = fault face, keep = other shared face), per patch and
//       over all faces; reported.
//   P1  every decision is computed from gathered, partition-invariant data
//       by the same operations on every rank.
//   P2  Orient(fi, n_raw) is sign-equivariant: the same answer for +-n_raw.
//   S3  side test by element centroids, both incident elements, opposite
//       signs and margins above d_F + floor asserted; S3(b) planarity.
//   S4  never Elem1/Elem2 ordering.
//
// Vertex keys (implementation decision, Chunhui Zhao, 2026-10-14).  The spec
// (rev 8.3, frozen) keys seeds (O2), canonical cycles (O3), the F6 fallback
// edge and the internal patch key (J5b) by global vertex ids and lists them
// as partition-invariant input (P1).  MFEM's ParMesh::GetGlobalVertexIndices
// (mesh/pmesh.cpp:6629-6643) returns the order-1 H1 global true-DOF number,
// numbered by owning rank, so it changes with the partition (probe: the
// vertex (1,1,1) of a 24-tet box is 17 / 30 / 28 at 1 / 2 / 4 ranks);
// pmesh.hpp:416-418 says these indices "have no meaning or significance for
// ParMesh", and the serial-mesh numbering map built in the ParMesh
// constructor (vert_global_local, pmesh.cpp:177) is discarded.  Spec P1 also
// lists vertex coordinates as invariant input, and they are copied bit for
// bit from the serial mesh, so the key of a fault vertex here is its rank in
// the lexicographic (x, y, z) order of the gathered fault vertices: the same
// on every rank and every rank count.  MFEM ids are used only to gather and
// de-duplicate faces within one run.  This is an implementation note, not a
// spec change.
//
// Supported geometry: planar affine triangular fault faces (every TPV102 and
// BP5 tet mesh).  Not implemented here: the plus_region pin, F8 transport,
// curved-face representation, junction policies beyond the census.
//
// Scope (Chunhui Zhao, 2026-10-14): tetrahedral meshes only.

#ifndef MFEM_SEAS_FAULT_ORIENTATION_HPP
#define MFEM_SEAS_FAULT_ORIENTATION_HPP

#include "mfem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

/// The small configuration N-08 needs (the plan's FaultGeometryConfig does
/// not exist in the repository).
struct FaultOrientationConfig
{
   enum class Pin { Reference, Point };
   Pin pin = Pin::Reference;
   real_t ref_normal[3] = {0.0, -1.0, 0.0};   ///< TPV102 / BP5 reference pin
   real_t up[3] = {0.0, 0.0, 1.0};
   real_t pin_point[3] = {0.0, 0.0, 0.0};     ///< Pin::Point: n points away
   real_t c_min = 0.1;                ///< J6: |c| below this refuses the run
   real_t side_floor_frac = 0.25;     ///< S3(c): floor = frac * (h_perp / 4)
   real_t planarity_rel_tol = 1e-12;  ///< S3(b): d_F <= tol * face diameter
   std::set<int> fault_attrs = {3};   ///< boundary attrs of fault faces (J5a)
};

struct FaultOrientationReport
{
   int faces = 0;              ///< global fault faces
   int components = 0;         ///< patches
   int junction_edges = 0;     ///< edges with m >= 3
   real_t flipped_vs_raw = 0;  ///< faces with n_raw . n < 0 (shared faces count 1/2)
   real_t min_pin_margin = 1;  ///< min over patches of |c| (J6)
   real_t max_face_departure = 0;   ///< max d_F (S3(b))
   real_t min_side_margin = 1e300;  ///< min |(c_E - c_F).n| / (h_perp / 4)
   bool parity_per_patch_ok = true;  ///< J3 with S = each patch
   bool parity_global_ok = true;     ///< J3 with S = all fault faces
};

class FaultOrientation
{
public:
   using Key3 = std::array<long long, 3>;

   FaultOrientation(ParMesh &pm, const Array<int> &interior_faces,
                    const Array<int> &shared_faces,
                    const FaultOrientationConfig &cfg = {})
      : pm_(pm), cfg_(cfg)
   {
      comm_ = pm.GetComm();
      pm.ExchangeFaceNbrData();
      n_int_ = interior_faces.Size();
      const int nloc = n_int_ + shared_faces.Size();
      local_face_.resize(nloc);
      for (int i = 0; i < n_int_; i++) { local_face_[i] = interior_faces[i]; }
      for (int i = 0; i < shared_faces.Size(); i++)
      {
         local_face_[n_int_ + i] = pm.GetSharedFace(shared_faces[i]);
         shared_index_.push_back(shared_faces[i]);
      }
      Gather_();
      BuildKeys_();
      Census_();
      Propagate_();
      Pin_();
      Parity_();
      Local_();
   }

   // ---- per local fault face fi (interior faces first, then shared) ----
   int NumLocalFaces() const { return static_cast<int>(local_face_.size()); }
   const real_t *Normal(int fi) const { return gface_[lglob_[fi]].n; }
   int Patch(int fi) const { return gface_[lglob_[fi]].patch; }
   int PersistentId(int fi) const { return gface_[lglob_[fi]].attr; }
   /// Local element on the + side (S1, S3), or -1 if it is on another rank.
   int PlusElement(int fi) const { return plus_elem_[fi]; }
   real_t SideMargin(int fi) const { return side_margin_[fi]; }
   const Key3 &FaceKey(int fi) const { return gface_[lglob_[fi]].key; }
   /// The face's vertices (keys + coordinates) for the F6 fallback.
   void FaceVertexData(int fi, long long key[3], real_t x[3][3]) const
   {
      const GFace &g = gface_[lglob_[fi]];
      for (int a = 0; a < 3; a++)
      {
         key[a] = g.key[a];
         for (int d = 0; d < 3; d++) { x[a][d] = g.x[a][d]; }
      }
   }
   /// P2: the oriented unit normal from any raw normal of face fi; returns
   /// the same vector for n_raw and -n_raw.
   void Orient(int fi, const real_t n_raw[3], real_t out[3]) const
   {
      const real_t *n = Normal(fi);
      const real_t c = n_raw[0] * n[0] + n_raw[1] * n[1] + n_raw[2] * n[2];
      const real_t nl = std::sqrt(n_raw[0] * n_raw[0] + n_raw[1] * n_raw[1]
                                  + n_raw[2] * n_raw[2]);
      const real_t s = (c >= 0.0) ? 1.0 : -1.0;
      for (int d = 0; d < 3; d++) { out[d] = s * n_raw[d] / nl; }
   }

   const FaultOrientationReport &GetReport() const { return rep_; }

private:
   struct GFace
   {
      Key3 gid;                 // sorted MFEM global ids (gather / dedupe only)
      Key3 key;                 // sorted partition-invariant vertex keys
      real_t x[3][3];           // coordinates, in key order
      int attr = -1;            // persistent id (J5a)
      real_t nC[3];             // right-hand normal of the canonical cycle
      real_t area = 0;
      int patch = -1;
      int s = 0;                // O3 sign
      real_t n[3];              // oriented normal
   };

   // ---------------------------------------------------------------- gather
   void Gather_()
   {
      Array<HYPRE_BigInt> gv;
      pm_.GetGlobalVertexIndices(gv);
      std::map<int, int> face_attr;
      for (int be = 0; be < pm_.GetNBE(); be++)
      {
         const int a = pm_.GetBdrAttribute(be);
         if (cfg_.fault_attrs.count(a))
         {
            face_attr[pm_.GetBdrElementFaceIndex(be)] = a;
         }
      }
      // record: 3 gids, 9 coords, attr
      const int R = 13;
      std::vector<double> loc;
      for (int fi = 0; fi < NumLocalFaces(); fi++)
      {
         const int f = local_face_[fi];
         MFEM_VERIFY(pm_.GetFaceGeometry(f) == Geometry::TRIANGLE,
                     "FaultOrientation: fault face " << f << " is a "
                     << Geometry::Name[pm_.GetFaceGeometry(f)]
                     << "; triangular fault faces only.");
         Array<int> fv;
         pm_.GetFaceVertices(f, fv);
         for (int a = 0; a < 3; a++) { loc.push_back(static_cast<double>(gv[fv[a]])); }
         for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++) { loc.push_back(pm_.GetVertex(fv[a])[d]); }
         auto it = face_attr.find(f);
         loc.push_back(it == face_attr.end() ? -1.0 : static_cast<double>(it->second));
      }
      int nl = static_cast<int>(loc.size()), np = 0;
      MPI_Comm_size(comm_, &np);
      std::vector<int> cnt(np), dsp(np, 0);
      MPI_Allgather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm_);
      for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      std::vector<double> all(dsp[np - 1] + cnt[np - 1]);
      MPI_Allgatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(),
                     dsp.data(), MPI_DOUBLE, comm_);
      // de-duplicate (shared faces appear on two ranks)
      std::map<Key3, int> seen;
      for (size_t o = 0; o + R <= all.size(); o += R)
      {
         long long g[3] = {std::llround(all[o]), std::llround(all[o + 1]),
                           std::llround(all[o + 2])};
         int ord[3] = {0, 1, 2};
         std::sort(ord, ord + 3, [&](int a, int b) { return g[a] < g[b]; });
         Key3 k = {g[ord[0]], g[ord[1]], g[ord[2]]};
         const int attr = static_cast<int>(std::llround(all[o + 12]));
         auto it = seen.find(k);
         if (it != seen.end())
         {
            gface_[it->second].attr = std::max(gface_[it->second].attr, attr);
            continue;
         }
         GFace F;
         F.gid = k;
         for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++) { F.x[a][d] = all[o + 3 + 3 * ord[a] + d]; }
         F.attr = attr;
         seen[k] = static_cast<int>(gface_.size());
         gface_.push_back(F);
      }
      // local face -> gathered face (by MFEM gids, consistent within the run)
      lglob_.resize(NumLocalFaces());
      for (int fi = 0; fi < NumLocalFaces(); fi++)
      {
         Array<int> fv;
         pm_.GetFaceVertices(local_face_[fi], fv);
         long long g[3] = {gv[fv[0]], gv[fv[1]], gv[fv[2]]};
         std::sort(g, g + 3);
         lglob_[fi] = seen.at(Key3{g[0], g[1], g[2]});
      }
   }

   // ------------------------------------------------- partition-invariant keys
   void BuildKeys_()
   {
      std::map<long long, std::array<real_t, 3>> vx;
      for (const GFace &F : gface_)
         for (int a = 0; a < 3; a++) { vx[F.gid[a]] = {F.x[a][0], F.x[a][1], F.x[a][2]}; }
      // key = rank in lexicographic (x, y, z) order (see the header note)
      std::map<long long, long long> key;
      std::vector<std::pair<std::array<real_t, 3>, long long>> v;
      for (const auto &kv : vx) { v.push_back({kv.second, kv.first}); }
      std::sort(v.begin(), v.end());
      for (size_t i = 0; i < v.size(); i++) { key[v[i].second] = static_cast<long long>(i); }
      for (GFace &F : gface_)
      {
         long long k[3] = {key[F.gid[0]], key[F.gid[1]], key[F.gid[2]]};
         real_t xs[3][3];
         int ord[3] = {0, 1, 2};
         std::sort(ord, ord + 3, [&](int a, int b) { return k[a] < k[b]; });
         for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++) { xs[a][d] = F.x[ord[a]][d]; }
         F.key = {k[ord[0]], k[ord[1]], k[ord[2]]};
         for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++) { F.x[a][d] = xs[a][d]; }
         // canonical cycle (k0 < k1 < k2): right-hand normal, area
         real_t e1[3], e2[3], c[3];
         for (int d = 0; d < 3; d++) { e1[d] = F.x[1][d] - F.x[0][d]; e2[d] = F.x[2][d] - F.x[0][d]; }
         c[0] = e1[1] * e2[2] - e1[2] * e2[1];
         c[1] = e1[2] * e2[0] - e1[0] * e2[2];
         c[2] = e1[0] * e2[1] - e1[1] * e2[0];
         const real_t cl = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
         MFEM_VERIFY(cl > 0.0, "FaultOrientation: degenerate fault face");
         for (int d = 0; d < 3; d++) { F.nC[d] = c[d] / cl; }
         F.area = 0.5 * cl;
      }
      // deterministic global face order: by sorted key triple
      std::vector<int> perm(gface_.size());
      for (size_t i = 0; i < perm.size(); i++) { perm[i] = static_cast<int>(i); }
      std::sort(perm.begin(), perm.end(),
                [&](int a, int b) { return gface_[a].key < gface_[b].key; });
      std::vector<GFace> sorted(gface_.size());
      std::vector<int> inv(gface_.size());
      for (size_t i = 0; i < perm.size(); i++) { sorted[i] = gface_[perm[i]]; inv[perm[i]] = static_cast<int>(i); }
      gface_.swap(sorted);
      for (int &g : lglob_) { g = inv[g]; }
      rep_.faces = static_cast<int>(gface_.size());
   }

   // canonical-cycle edge direction tau_f(e) for edge (a, b), a < b, of face F
   static int Tau_(const GFace &F, long long a, long long b)
   {
      // cycle k0 -> k1 -> k2 -> k0
      if (a == F.key[0] && b == F.key[2]) { return -1; }
      return +1;
   }

   // --------------------------------------------------------- J1 census
   void Census_()
   {
      for (int g = 0; g < static_cast<int>(gface_.size()); g++)
      {
         const Key3 &k = gface_[g].key;
         edges_[{k[0], k[1]}].push_back(g);
         edges_[{k[1], k[2]}].push_back(g);
         edges_[{k[0], k[2]}].push_back(g);
      }
      for (const auto &e : edges_) { if (e.second.size() >= 3) { rep_.junction_edges++; } }
   }

   // ------------------------------------------------ O1-O3, O5, O6
   void Propagate_()
   {
      const int nf = static_cast<int>(gface_.size());
      std::vector<std::vector<std::pair<int, std::pair<long long, long long>>>> adj(nf);
      for (const auto &e : edges_)
      {
         if (e.second.size() != 2) { continue; }   // tips and junction edges cut
         const int a = e.second[0], b = e.second[1];
         adj[a].push_back({b, e.first});
         adj[b].push_back({a, e.first});
      }
      int patch = 0;
      for (int seed = 0; seed < nf; seed++)       // faces sorted by key: O2
      {
         if (gface_[seed].patch >= 0) { continue; }
         std::vector<int> stack = {seed};
         gface_[seed].patch = patch;
         gface_[seed].s = +1;
         while (!stack.empty())
         {
            const int f = stack.back();
            stack.pop_back();
            for (const auto &nb : adj[f])
            {
               const int g = nb.first;
               const long long ea = nb.second.first, eb = nb.second.second;
               const int sg = -gface_[f].s * Tau_(gface_[f], ea, eb) * Tau_(gface_[g], ea, eb);
               if (gface_[g].patch < 0)
               {
                  gface_[g].patch = patch;
                  gface_[g].s = sg;
                  stack.push_back(g);
               }
               else
               {
                  MFEM_VERIFY(gface_[g].s == sg,
                     "FaultOrientation (O6): patch " << patch << " is not "
                     "orientable; inconsistent loop at the edge between "
                     "faces with keys (" << gface_[f].key[0] << ","
                     << gface_[f].key[1] << "," << gface_[f].key[2] << ") and ("
                     << gface_[g].key[0] << "," << gface_[g].key[1] << ","
                     << gface_[g].key[2] << ").");
               }
            }
         }
         patch++;
      }
      rep_.components = patch;
   }

   // ----------------------------------------------------------- O4 pin, J6
   void Pin_()
   {
      std::vector<real_t> num(rep_.components, 0.0), den(rep_.components, 0.0);
      std::vector<int> seed(rep_.components, -1);
      for (int g = 0; g < static_cast<int>(gface_.size()); g++)   // fixed order
      {
         const GFace &F = gface_[g];
         if (seed[F.patch] < 0) { seed[F.patch] = g; }
         const real_t dn = F.s * (F.nC[0] * cfg_.ref_normal[0] + F.nC[1] * cfg_.ref_normal[1]
                                  + F.nC[2] * cfg_.ref_normal[2]);
         num[F.patch] += F.area * dn;
         den[F.patch] += F.area;
      }
      std::vector<int> pin(rep_.components, 1);
      rep_.min_pin_margin = 1.0;
      for (int p = 0; p < rep_.components; p++)
      {
         if (cfg_.pin == FaultOrientationConfig::Pin::Reference)
         {
            const real_t c = num[p] / den[p];
            MFEM_VERIFY(std::abs(c) >= cfg_.c_min,
               "FaultOrientation (J6): the reference pin of patch " << p
               << " is ill-conditioned: c = " << c << ", |c| < c_min = "
               << cfg_.c_min << ". Refusing the run.");
            pin[p] = (c < 0.0) ? -1 : 1;
            rep_.min_pin_margin = std::min(rep_.min_pin_margin, std::abs(c));
         }
         else   // Point: n points away from the given point (seed face)
         {
            const GFace &F = gface_[seed[p]];
            real_t cf[3], v = 0.0;
            for (int d = 0; d < 3; d++)
            {
               cf[d] = (F.x[0][d] + F.x[1][d] + F.x[2][d]) / 3.0;
               v += (cf[d] - cfg_.pin_point[d]) * F.s * F.nC[d];
            }
            MFEM_VERIFY(v != 0.0, "FaultOrientation: point pin undecided for patch " << p);
            pin[p] = (v < 0.0) ? -1 : 1;
         }
      }
      for (GFace &F : gface_)
      {
         for (int d = 0; d < 3; d++) { F.n[d] = pin[F.patch] * F.s * F.nC[d]; }
      }
   }

   // ------------------------------------------------------------- J3 parity
   void Parity_()
   {
      Array<HYPRE_BigInt> gv;
      pm_.GetGlobalVertexIndices(gv);
      // fault edges (by MFEM gid pairs, consistent within the run) with m >= 2
      std::map<std::pair<long long, long long>, int> edge_m;
      std::set<Key3> fault_tri;
      for (const GFace &F : gface_)
      {
         fault_tri.insert(F.gid);
         edge_m[{F.gid[0], F.gid[1]}]++;
         edge_m[{F.gid[1], F.gid[2]}]++;
         edge_m[{F.gid[0], F.gid[2]}]++;
      }
      std::set<long long> fverts;
      for (const GFace &F : gface_) { for (int a = 0; a < 3; a++) { fverts.insert(F.gid[a]); } }
      // local elements touching the fault (a fault face or an interior fault edge)
      std::vector<double> loc;
      for (int e = 0; e < pm_.GetNE(); e++)
      {
         Array<int> ev;
         pm_.GetElementVertices(e, ev);
         std::vector<long long> g;
         int nf = 0;
         for (int a = 0; a < ev.Size(); a++) { g.push_back(gv[ev[a]]); if (fverts.count(gv[ev[a]])) { nf++; } }
         if (nf < 2) { continue; }
         std::sort(g.begin(), g.end());
         bool touch = false;
         for (int a = 0; a < 4 && !touch; a++)
            for (int b = a + 1; b < 4 && !touch; b++)
            {
               auto it = edge_m.find({g[a], g[b]});
               if (it != edge_m.end() && it->second >= 2) { touch = true; }
            }
         for (int a = 0; a < 4 && !touch; a++)
         {
            Key3 t;
            int c = 0;
            for (int b = 0; b < 4; b++) { if (b != a) { t[c++] = g[b]; } }
            if (fault_tri.count(t)) { touch = true; }
         }
         if (touch) { for (int a = 0; a < 4; a++) { loc.push_back(static_cast<double>(g[a])); } }
      }
      int nl = static_cast<int>(loc.size()), np = 0;
      MPI_Comm_size(comm_, &np);
      std::vector<int> cnt(np), dsp(np, 0);
      MPI_Allgather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm_);
      for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      std::vector<double> all(dsp[np - 1] + cnt[np - 1]);
      MPI_Allgatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(),
                     dsp.data(), MPI_DOUBLE, comm_);
      std::set<std::array<long long, 4>> elems;
      for (size_t o = 0; o + 4 <= all.size(); o += 4)
      {
         elems.insert({std::llround(all[o]), std::llround(all[o + 1]),
                       std::llround(all[o + 2]), std::llround(all[o + 3])});
      }
      std::vector<std::array<long long, 4>> E(elems.begin(), elems.end());
      // per patch, then global (patch = -1)
      std::map<Key3, int> patch_of;
      for (const GFace &F : gface_) { patch_of[F.gid] = F.patch; }
      for (int p = -1; p < rep_.components; p++)
      {
         const bool ok = ParityOK_(E, patch_of, p);
         if (p < 0) { rep_.parity_global_ok = ok; }
         else if (!ok) { rep_.parity_per_patch_ok = false; }
      }
   }

   // J3 parity BFS for S = patch p (p < 0: all faces).
   static bool ParityOK_(const std::vector<std::array<long long, 4>> &E,
                         const std::map<Key3, int> &patch_of, int p)
   {
      auto inS = [&](const Key3 &t) {
         auto it = patch_of.find(t);
         return it != patch_of.end() && (p < 0 || it->second == p);
      };
      // interior edges of S: in >= 2 faces of S
      std::map<std::pair<long long, long long>, int> em;
      for (const auto &kv : patch_of)
      {
         if (!inS(kv.first)) { continue; }
         const Key3 &k = kv.first;
         em[{k[0], k[1]}]++; em[{k[1], k[2]}]++; em[{k[0], k[2]}]++;
      }
      // N_S
      std::vector<int> N;
      for (int i = 0; i < static_cast<int>(E.size()); i++)
      {
         const auto &g = E[i];
         bool in = false;
         for (int a = 0; a < 4 && !in; a++)
         {
            Key3 t; int c = 0;
            for (int b = 0; b < 4; b++) { if (b != a) { t[c++] = g[b]; } }
            if (inS(t)) { in = true; }
         }
         for (int a = 0; a < 4 && !in; a++)
            for (int b = a + 1; b < 4 && !in; b++)
            {
               auto it = em.find({g[a], g[b]});
               if (it != em.end() && it->second >= 2) { in = true; }
            }
         if (in) { N.push_back(i); }
      }
      // H_S edges via shared triangles
      std::map<Key3, std::vector<int>> tri;
      for (int idx = 0; idx < static_cast<int>(N.size()); idx++)
      {
         const auto &g = E[N[idx]];
         for (int a = 0; a < 4; a++)
         {
            Key3 t; int c = 0;
            for (int b = 0; b < 4; b++) { if (b != a) { t[c++] = g[b]; } }
            tri[t].push_back(idx);
         }
      }
      std::vector<std::vector<std::pair<int, int>>> adj(N.size());
      for (const auto &kv : tri)
      {
         if (kv.second.size() != 2) { continue; }
         const int flip = inS(kv.first) ? 1 : 0;
         adj[kv.second[0]].push_back({kv.second[1], flip});
         adj[kv.second[1]].push_back({kv.second[0], flip});
      }
      std::vector<int> col(N.size(), -1);
      for (size_t s = 0; s < N.size(); s++)
      {
         if (col[s] >= 0) { continue; }
         col[s] = 0;
         std::vector<int> st = {static_cast<int>(s)};
         while (!st.empty())
         {
            const int u = st.back(); st.pop_back();
            for (const auto &e : adj[u])
            {
               const int want = col[u] ^ e.second;
               if (col[e.first] < 0) { col[e.first] = want; st.push_back(e.first); }
               else if (col[e.first] != want) { return false; }
            }
         }
      }
      return true;
   }

   // ------------------------------------------- per local face: S3, report
   void Local_()
   {
      const int nloc = NumLocalFaces();
      plus_elem_.assign(nloc, -1);
      side_margin_.assign(nloc, 0.0);
      real_t flipped = 0.0, max_dF = 0.0, min_side = 1e300;
      for (int fi = 0; fi < nloc; fi++)
      {
         const bool shared = (fi >= n_int_);
         const int f = local_face_[fi];
         FaceElementTransformations *Tr = shared
            ? pm_.GetSharedFaceTransformations(shared_index_[fi - n_int_])
            : pm_.GetInteriorFaceTransformations(f);
         MFEM_VERIFY(Tr, "FaultOrientation: no face transformation for face " << f);
         const GFace &G = gface_[lglob_[fi]];
         // face centre image and the oriented normal there
         const IntegrationPoint &ipc = Geometries.GetCenter(Tr->GetGeometryType());
         Tr->SetAllIntPoints(&ipc);
         Vector cF(3), nr(3);
         Tr->Transform(ipc, cF);
         CalcOrtho(Tr->Jacobian(), nr);
         real_t nc[3];
         Orient(fi, nr.GetData(), nc);
         const real_t s_c = (nr * Vector(const_cast<real_t *>(G.n), 3) >= 0.0) ? 1.0 : -1.0;
         if (nr * Vector(const_cast<real_t *>(G.n), 3) < 0.0) { flipped += shared ? 0.5 : 1.0; }
         // S3(b): planarity and normal agreement over quadrature + vertices
         real_t dF = 0.0, diam = 0.0;
         for (int a = 0; a < 3; a++)
            for (int b = a + 1; b < 3; b++)
            {
               real_t l = 0;
               for (int d = 0; d < 3; d++) { l += (G.x[a][d] - G.x[b][d]) * (G.x[a][d] - G.x[b][d]); }
               diam = std::max(diam, std::sqrt(l));
            }
         const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 6);
         for (int q = 0; q < ir.GetNPoints() + 3; q++)
         {
            IntegrationPoint ip;
            if (q < ir.GetNPoints()) { ip = ir.IntPoint(q); }
            else { const int a = q - ir.GetNPoints(); ip.Set2(a == 1 ? 1.0 : 0.0, a == 2 ? 1.0 : 0.0); }
            Tr->SetAllIntPoints(&ip);
            Vector xq(3), nq(3);
            Tr->Transform(ip, xq);
            real_t dep = 0.0;
            for (int d = 0; d < 3; d++) { dep += (xq(d) - cF(d)) * nc[d]; }
            dF = std::max(dF, std::abs(dep));
            CalcOrtho(Tr->Jacobian(), nq);
            // n(x_q): the mapped normal with the centre's orientation sign
            // s_c (the face map's normal is continuous on the face).
            real_t agree = 0.0;
            for (int d = 0; d < 3; d++) { agree += s_c * nq(d) / nq.Norml2() * nc[d]; }
            MFEM_VERIFY(agree > 0.0,
               "FaultOrientation (S3(b)): face " << f << " (keys " << G.key[0] << ","
               << G.key[1] << "," << G.key[2] << "): mapped normal disagrees with "
               "the face-centre normal at a quadrature point (n(x_q).n_c = " << agree << ")");
         }
         max_dF = std::max(max_dF, dF);
         MFEM_VERIFY(dF <= cfg_.planarity_rel_tol * diam,
            "FaultOrientation (S3(b)): face " << f << " (keys " << G.key[0] << ","
            << G.key[1] << "," << G.key[2] << ") is not planar: d_F = " << dF
            << " m > " << cfg_.planarity_rel_tol << " * diameter " << diam
            << " m. Refusing to start (planar affine fault faces only).");
         // S3(a),(c): side test by centroids along the oriented normal.
         auto margin = [&](ElementTransformation *T, bool local, int e, real_t &hq) {
            Array<int> ev;
            real_t ce[3] = {0, 0, 0};
            if (local)
            {
               pm_.GetElementVertices(e, ev);
               for (int a = 0; a < ev.Size(); a++)
                  for (int d = 0; d < 3; d++) { ce[d] += pm_.GetVertex(ev[a])[d] / ev.Size(); }
            }
            else
            {
               const IntegrationPoint &ipe = Geometries.GetCenter(T->GetGeometryType());
               Vector x(3);
               T->Transform(ipe, x);
               for (int d = 0; d < 3; d++) { ce[d] = x(d); }
            }
            // tet: h_perp = 3 V / A_face, expected centroid margin h_perp / 4
            const IntegrationPoint &ipe = Geometries.GetCenter(T->GetGeometryType());
            T->SetIntPoint(&ipe);
            const real_t V = T->Weight() / 6.0;
            hq = 3.0 * V / G.area / 4.0;
            real_t m = 0.0;
            for (int d = 0; d < 3; d++) { m += (ce[d] - cF(d)) * nc[d]; }
            return m;
         };
         real_t h1, h2;
         const real_t m1 = margin(Tr->Elem1, true, Tr->Elem1No, h1);
         const real_t floor1 = cfg_.side_floor_frac * h1;
         MFEM_VERIFY(std::abs(m1) > dF + floor1,
            "FaultOrientation (S3(c)): face " << f << ": element " << Tr->Elem1No
            << " centroid margin " << std::abs(m1) << " <= d_F + floor = " << dF + floor1);
         real_t smin = std::abs(m1) / h1;
         int plus = (m1 < 0.0) ? Tr->Elem1No : -1;
         if (!shared)
         {
            const real_t m2 = margin(Tr->Elem2, true, Tr->Elem2No, h2);
            const real_t floor2 = cfg_.side_floor_frac * h2;
            MFEM_VERIFY(std::abs(m2) > dF + floor2,
               "FaultOrientation (S3(c)): face " << f << ": element " << Tr->Elem2No
               << " centroid margin " << std::abs(m2) << " <= d_F + floor = " << dF + floor2);
            MFEM_VERIFY((m1 < 0.0) != (m2 < 0.0),
               "FaultOrientation (S3(c)): face " << f << ": both incident elements on "
               "the same side (margins " << m1 << ", " << m2 << "); inverted or "
               "degenerate element.");
            smin = std::min(smin, std::abs(m2) / h2);
            if (m2 < 0.0) { plus = Tr->Elem2No; }
         }
         plus_elem_[fi] = plus;
         side_margin_[fi] = smin;
         min_side = std::min(min_side, smin);
      }
      MPI_Allreduce(MPI_IN_PLACE, &flipped, 1, MPI_DOUBLE, MPI_SUM, comm_);
      MPI_Allreduce(MPI_IN_PLACE, &max_dF, 1, MPI_DOUBLE, MPI_MAX, comm_);
      MPI_Allreduce(MPI_IN_PLACE, &min_side, 1, MPI_DOUBLE, MPI_MIN, comm_);
      rep_.flipped_vs_raw = flipped;
      rep_.max_face_departure = max_dF;
      rep_.min_side_margin = min_side;
   }

   ParMesh &pm_;
   FaultOrientationConfig cfg_;
   MPI_Comm comm_;
   int n_int_ = 0;
   std::vector<int> local_face_;
   std::vector<int> shared_index_;
   std::vector<GFace> gface_;
   std::vector<int> lglob_;
   std::map<std::pair<long long, long long>, std::vector<int>> edges_;
   std::vector<int> plus_elem_;
   std::vector<real_t> side_margin_;
   FaultOrientationReport rep_;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_ORIENTATION_HPP
