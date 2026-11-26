// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/ucerf3d/test_a_field_3d.cpp — UCERF 3-D plan rev5, P1.2 (Chunhui
// Zhao, 2026-11-25: "start P1"; decision card 11-25 "Straight-line to the tip
// edges"): the tip and depth a-field of [cycle.fields.a] kind = "tip_zone" on
// 3-D fault surfaces (cycle/fault_surface_zones.hpp).
//
//   U-A1  on T-plane and T-dip (ucerf3d/mesh; the 40 x 20 km fault, vertical
//         and dipping 58 deg, creep surface below, attr 6): one surface, its
//         lateral tip edges the two ends x = +-20 km, the free-surface and
//         creep-seam edges not tips; at every vertex, edge midpoint and
//         centroid of the local fault triangles, a = the analytic value with
//         d = 20 km - |x| and z_d = -z (W 8 km, taper 2 km; VS above 2 km,
//         taper 2 km; VS below 15 km, taper 2 km) to 1e-12
//   U-A1b the same on T-plane with lower_depth_file: one chain y = 0, D
//         10 km for x < -1 km, 14 km for x > 1 km, linear in between
//   U-A1c loading option B: a box with a buried fault and nothing below it:
//         its flat bottom edge is free but not a tip, its two vertical ends are
//         tips; a = the analytic value
//   U-A1d partition invariance: rank 0 prints "AF3-DIGEST" (a hash of every
//         (x, a) sample, sorted); the make target checks 1, 2 and 4 ranks
//
//   seas_test_a_field_3d              (1, 2 or 4 ranks; AF3_DUMP=DIR writes
//                                     each case's samples x y z a to DIR)
//   seas_test_a_field_3d --refuse N   one refusal, must exit non-zero:
//     1 a T-shaped fault ("three or more fault faces"),
//     2 a point off the fault ("lies on no fault triangle")
//
//   make -f tests/ucerf3d/ucerf3d.mk test-a-field-3d
//
// Scope (Chunhui Zhao, 2026-11-25): tetrahedral meshes, straight faces.

#include "mfem.hpp"
#include "../../cycle/fault_surface_zones.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kIn = 0.004, kTip = 0.08;

static cycle::TipZoneSpec Spec()
{
   cycle::TipZoneSpec s;
   s.width = 8e3;
   s.taper = 2e3;
   s.inside = kIn;
   s.tip = kTip;
   s.top_depth = 2e3;
   s.top_taper = 2e3;
   s.lower_depth = 15e3;
   s.bottom_taper = 2e3;
   s.set = true;
   s.Validate("a");
   return s;
}

static real_t Ramp(real_t lo, real_t hi, real_t v)
{
   if (v <= lo) { return 1.0; }
   if (v >= hi) { return 0.0; }
   return (hi - v) / (hi - lo);
}

/// The analytic value: d to the nearer end x = +-xe, depth -z, lower depth D.
static real_t Exact(const cycle::TipZoneSpec &s, const real_t x[3], real_t xe, real_t D)
{
   const real_t d = xe - std::abs(x[0]), zd = -x[2];
   real_t w = Ramp(s.width - s.taper, s.width, d);
   w = std::max(w, Ramp(s.top_depth, s.top_depth + s.top_taper, zd));
   w = std::max(w, 1.0 - Ramp(D, D + s.bottom_taper, zd));
   return s.inside + w * (s.tip - s.inside);
}

/// Sample points of this rank's fault triangles: vertices, edge midpoints,
/// centroids.
static std::vector<std::array<real_t, 3>> Samples(const ParMesh &pm, int attr)
{
   std::vector<std::array<real_t, 3>> s;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) != attr) { continue; }
      Array<int> v;
      pm.GetBdrElementVertices(be, v);
      const real_t *p[3] = {pm.GetVertex(v[0]), pm.GetVertex(v[1]), pm.GetVertex(v[2])};
      for (int a = 0; a < 3; a++)
      {
         s.push_back({p[a][0], p[a][1], p[a][2]});
         const real_t *q = p[(a + 1) % 3];
         s.push_back({0.5 * (p[a][0] + q[0]), 0.5 * (p[a][1] + q[1]), 0.5 * (p[a][2] + q[2])});
      }
      s.push_back({(p[0][0] + p[1][0] + p[2][0]) / 3, (p[0][1] + p[1][1] + p[2][1]) / 3,
                   (p[0][2] + p[1][2] + p[2][2]) / 3});
   }
   return s;
}

/// FNV-1a over every rank's (x, a), sorted (rank 0; "" elsewhere).
static std::string Digest(std::vector<std::array<real_t, 4>> rows, MPI_Comm comm,
                          const std::string &tag)
{
   std::vector<double> loc;
   for (const auto &r : rows) { loc.insert(loc.end(), r.begin(), r.end()); }
   int np = 1, nl = static_cast<int>(loc.size());
   MPI_Comm_size(comm, &np);
   std::vector<int> cnt(np), dsp(np, 0);
   MPI_Gather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm);
   std::vector<double> all;
   if (g_rank == 0)
   {
      for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      all.resize(dsp[np - 1] + cnt[np - 1]);
   }
   MPI_Gatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(), MPI_DOUBLE, 0,
               comm);
   if (g_rank != 0) { return ""; }
   std::vector<std::array<real_t, 4>> g;
   for (size_t o = 0; o + 4 <= all.size(); o += 4) { g.push_back({all[o], all[o + 1], all[o + 2], all[o + 3]}); }
   std::sort(g.begin(), g.end());
   g.erase(std::unique(g.begin(), g.end()), g.end());
   if (const char *dir = std::getenv("AF3_DUMP"))   // x y z a per sample, for the review figures
   {
      std::string f = tag;
      std::replace(f.begin(), f.end(), ' ', '_');
      std::ofstream o(std::string(dir) + "/" + f + ".txt");
      o.precision(17);
      for (const auto &r : g) { o << r[0] << " " << r[1] << " " << r[2] << " " << r[3] << "\n"; }
   }
   uint64_t h = 1469598103934665603ULL;
   for (const auto &r : g)
      for (real_t v : r)
      {
         unsigned char b[sizeof(double)];
         const double dv = v;
         std::memcpy(b, &dv, sizeof(double));
         for (unsigned char c : b) { h ^= c; h *= 1099511628211ULL; }
      }
   char s[64];
   std::snprintf(s, sizeof(s), "%016llx points %zu", static_cast<unsigned long long>(h), g.size());
   return s;
}

static void CheckMesh(const std::string &tag, ParMesh &pm, const cycle::TipZoneSpec &spec,
                      real_t xe, int patches, int lat_want, int flat_free,
                      real_t (*D)(const real_t *), MPI_Comm comm)
{
   cycle::FaultSurfaceZones Z(pm, 3, spec);
   if (g_rank == 0) { std::cout << "  " << Z.Report(); }
   int lat = 0;
   for (int p = 0; p < Z.NumPatches(); p++) { lat += Z.NumLateralEdges(p); }
   const std::string rep = Z.Report();
   const bool flat_ok = rep.find(std::to_string(flat_free) + " flat free edges") != std::string::npos;
   TEST_ASSERT(Z.NumPatches() == patches && lat == lat_want && flat_ok,
               tag << ": " << Z.NumPatches() << " surface(s), " << lat << " lateral tip edges, "
               << flat_free << " flat free edges (not tips)");
   real_t worst = 0.0, worst_d = 0.0;
   int n = 0, nt = 0, nr = 0, ni = 0;
   std::vector<std::array<real_t, 4>> rows;
   for (const auto &x : Samples(pm, 3))
   {
      int p = 0;
      real_t d = 0.0;
      Z.Locate(x.data(), p, d);
      worst_d = std::max(worst_d, std::abs(d - (xe - std::abs(x[0]))));
      const real_t a = Z.Value(x.data()), e = Exact(spec, x.data(), xe, D(x.data()));
      worst = std::max(worst, std::abs(a - e));
      n++;
      if (a == kTip) { nt++; } else if (a == kIn) { ni++; } else { nr++; }
      rows.push_back({x[0], x[1], x[2], a});
   }
   MPI_Allreduce(MPI_IN_PLACE, &worst, 1, MPI_DOUBLE, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, &worst_d, 1, MPI_DOUBLE, MPI_MAX, comm);
   int c[4] = {n, nt, nr, ni};
   MPI_Allreduce(MPI_IN_PLACE, c, 4, MPI_INT, MPI_SUM, comm);
   TEST_ASSERT(worst_d <= 1e-9 * xe, tag << ": tip distance d = " << xe / 1e3
               << " km - |x| at every sample (max error " << worst_d << " m)");
   TEST_ASSERT(worst <= 1e-12 && c[1] > 0 && c[2] > 0 && c[3] > 0,
               tag << ": a = the analytic value at " << c[0] << " samples (tip / ramp / inside "
               << c[1] << " / " << c[2] << " / " << c[3] << "), max |error| " << worst);
   const std::string dg = Digest(rows, comm, tag);
   if (g_rank == 0) { std::cout << "AF3-DIGEST " << tag << " " << dg << "\n"; }
}

static real_t D15(const real_t *) { return 15e3; }
static real_t DFile(const real_t *x)
{
   if (x[0] <= -1e3) { return 10e3; }
   if (x[0] >= 1e3) { return 14e3; }
   return 10e3 + (x[0] + 1e3) / 2e3 * 4e3;
}

// A Kuhn box (6 tets per cell, the face j = const split along (i, k)-(i+1,
// k+1)), cells of h = 1 km: x in [-10, 10] km, y in [-5, 5] km, z in [-16, 0]
// km; the fault y = 0, |x| <= 6 km, z >= -12 km (buried bottom, nothing below:
// loading option B).  Outer faces attr 5 (the top z = 0 too).  branch = true
// adds the half plane x = 0, 0 <= y <= 2 km, z >= -12 km (refusal 1).
static Mesh *BoxB(bool branch = false)
{
   const int ni = 20, nj = 10, nk = 16;
   const real_t h = 1e3;
   Mesh *m = new Mesh(3, (ni + 1) * (nj + 1) * (nk + 1), 6 * ni * nj * nk, 0, 3);
   auto v = [&](int i, int j, int k) { return i + (ni + 1) * (j + (nj + 1) * k); };
   for (int k = 0; k <= nk; k++)
      for (int j = 0; j <= nj; j++)
         for (int i = 0; i <= ni; i++) { m->AddVertex(-10e3 + i * h, -5e3 + j * h, -16e3 + k * h); }
   const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
   for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
         for (int i = 0; i < ni; i++)
            for (int p = 0; p < 6; p++)
            {
               int c[3] = {i, j, k}, t[4];
               t[0] = v(c[0], c[1], c[2]);
               for (int s = 0; s < 3; s++) { c[perm[p][s]]++; t[s + 1] = v(c[0], c[1], c[2]); }
               m->AddTet(t[0], t[1], t[2], t[3], 1);
            }
   // the quad of the grid face normal to `axis` at cell index (i, j, k), split
   // along its (low, low)-(high, high) diagonal as the Kuhn tets split it
   auto quad = [&](int axis, int i, int j, int k, int attr)
   {
      const int p = (axis == 0) ? 1 : 0, q = (axis == 2) ? 1 : 2;
      int c0[3] = {i, j, k}, c1[3] = {i, j, k}, c2[3] = {i, j, k}, c3[3] = {i, j, k};
      c1[p]++; c2[q]++; c3[p]++; c3[q]++;
      m->AddBdrTriangle(v(c0[0], c0[1], c0[2]), v(c1[0], c1[1], c1[2]), v(c3[0], c3[1], c3[2]), attr);
      m->AddBdrTriangle(v(c0[0], c0[1], c0[2]), v(c2[0], c2[1], c2[2]), v(c3[0], c3[1], c3[2]), attr);
   };
   const int n[3] = {ni, nj, nk};
   for (int axis = 0; axis < 3; axis++)
   {
      const int p = (axis == 0) ? 1 : 0, q = (axis == 2) ? 1 : 2;
      for (int a = 0; a <= n[axis]; a += n[axis])
         for (int b = 0; b < n[p]; b++)
            for (int c = 0; c < n[q]; c++)
            {
               int idx[3];
               idx[axis] = a; idx[p] = b; idx[q] = c;
               quad(axis, idx[0], idx[1], idx[2], 5);
            }
   }
   for (int i = 4; i < 16; i++)
      for (int k = 4; k < 16; k++) { quad(1, i, 5, k, 3); }
   if (branch)
      for (int j = 5; j < 7; j++)
         for (int k = 4; k < 16; k++) { quad(0, 10, j, k, 3); }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

/// The box's spec: the fault is 12 x 12 km, so a narrower tip zone (W 3 km,
/// taper 1 km) and VS below 9 km leave an inside region.
static cycle::TipZoneSpec BoxSpec()
{
   cycle::TipZoneSpec s = Spec();
   s.width = 3e3;
   s.taper = 1e3;
   s.lower_depth = 9e3;
   s.Validate("a");
   return s;
}
static real_t D9(const real_t *) { return 9e3; }

static int Refuse(int which, MPI_Comm comm)
{
   const cycle::TipZoneSpec spec = Spec();
   if (which == 1)
   {
      std::unique_ptr<Mesh> s(BoxB(true));
      ParMesh pm(comm, *s);
      cycle::FaultSurfaceZones Z(pm, 3, spec);   // aborts here
   }
   else if (which == 2)
   {
      std::unique_ptr<Mesh> s(BoxB());
      ParMesh pm(comm, *s);
      cycle::FaultSurfaceZones Z(pm, 3, spec);
      const real_t x[3] = {0.0, 3e3, -5e3};
      (void)Z.Value(x);   // aborts here
   }
   if (g_rank == 0) { std::cout << "refusal " << which << " did NOT fire\n"; }
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   g_rank = Mpi::WorldRank();
   if (argc == 3 && std::strcmp(argv[1], "--refuse") == 0)
   {
      return Refuse(std::atoi(argv[2]), comm) == 0 ? 0 : 3;
   }
   if (g_rank == 0)
   {
      std::cout << "test_a_field_3d (UCERF 3-D plan rev5, P1.2) on " << Mpi::WorldSize()
                << " rank(s)\n";
   }
   const cycle::TipZoneSpec spec = Spec();
   for (const std::string name : {"T-plane", "T-dip"})
   {
      if (g_rank == 0) { std::cout << "\n=== U-A1 " << name << " ===\n"; }
      Mesh sm(("ucerf3d/mesh/" + name + ".msh").c_str(), 1, 1);
      ParMesh pm(comm, sm);
      sm.Clear();
      CheckMesh("U-A1 " + name, pm, spec, 20e3, 1, (name == "T-plane") ? 28 : 32, 0, D15, comm);
   }
   {
      if (g_rank == 0) { std::cout << "\n=== U-A1b T-plane, lower_depth_file ===\n"; }
      const char *tmp = std::getenv("TMPDIR");
      const std::string csv = std::string(tmp ? tmp : "/tmp") + "/seas_test_a_field_3d_lower_depth.csv";
      if (g_rank == 0)
      {
         std::ofstream f(csv);
         f << "chain,x_m,y_m,lower_depth_m\n" << "c1,-30000,0,10000\nc1,-1000,0,10000\n"
           << "c1,1000,0,14000\nc1,30000,0,14000\n";
      }
      MPI_Barrier(comm);
      cycle::TipZoneSpec s2 = spec;
      s2.lower_depth = -1.0;
      s2.lower_depth_file = csv;
      s2.Validate("a");
      Mesh sm("ucerf3d/mesh/T-plane.msh", 1, 1);
      ParMesh pm(comm, sm);
      sm.Clear();
      CheckMesh("U-A1b T-plane file", pm, s2, 20e3, 1, 28, 0, DFile, comm);
   }
   {
      if (g_rank == 0) { std::cout << "\n=== U-A1c buried fault, nothing below (option B) ===\n"; }
      std::unique_ptr<Mesh> sm(BoxB());
      ParMesh pm(comm, *sm);
      CheckMesh("U-A1c box", pm, BoxSpec(), 6e3, 1, 24, 12, D9, comm);
   }
   if (g_rank == 0)
   {
      std::cout << "\n" << (num_tests - num_failed) << "/" << num_tests << " tests passed\n";
   }
   return num_failed == 0 ? 0 : 1;
}
