// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_vertex_normals.cpp — vertex-normal plan rev1, phase
// P1 (Chunhui Zhao, 2026-11-01: "send this document to 2d thread and start
// implementation and test"; /mnt/project-files/2d-inplane/
// nonplanar-literature/vertex-normal-plan-rev1.pdf): the averaged vertex
// normals and the frame along a face of FaultFrameField
// (fault/fault_frame_field.hpp, [fault_geometry] vertex_normals = "averaged"),
// the plan's tests T2 (geometry) and T3 (orientation).
//
// Mesh: an annulus sector around the centre (0, R), R = 10 km, angle phi in
// [-24, 24] deg in 24 cells of 2 deg, radius r in [R - 3 km, R + 3 km] in 6
// cells; a vertex at (r sin phi, R - r cos phi), so the grid line r = R is a
// polygon inscribed in the circle, through (0, 0).  The fault (attr 3) is that
// line for |phi| <= 12 deg: 12 straight faces, 11 interior vertices with a
// turn of 2 deg each, 2 tips.  Outer sides attr 5.  Triangles: each cell split
// along its diagonal.
//
//   T2a  report: 11 smooth vertices (largest turn 2 deg), 0 kinks, 2 tips
//   T2b  every smooth vertex normal is the circle's radial unit vector
//        (sin phi, -cos phi) to 1e-13; a tip keeps its face normal bit for bit
//   T2c  FrameAt: at a vertex n = the vertex normal; at every point n, d, s
//        orthonormal, d x s = n (1e-14); at the chord's midpoint n = the
//        normalized mean of the two vertex normals (1e-14), and the slot rule
//        (up = e_z: d = -+e_z; up = e_x: s = -+e_z): the out-of-plane vector
//        has exactly zero in-plane components, n and the in-plane tangent
//        exactly zero z components
//   T2d  threshold 1 deg: 0 smooth, 11 kinks kept (largest turn 2 deg); every
//        VertexNormal and FrameAt is the face's own, bit for bit
//   T2e  the option off ("face"): no counts, FrameAt = Frame(fi) bit for bit
//   T2f  a straight fault (y = 0) with the option on: every vertex smooth with
//        turn 0, every VertexNormal and FrameAt the face's own bit for bit
//   T2g  partition invariance: rank 0 prints "VN-DIGEST" (every vertex normal,
//        sorted by position, in hex) for METIS's partition and for one that
//        cuts along the fault (elements with r > R on one rank, r < R on the
//        other, so every fault face is a rank-boundary face; on 4 ranks also
//        split at phi = 0); the make target checks that every line of the
//        1-, 2- and 4-rank runs is the same
//   T3   orientation: with ref_normal (0, -1), (0, 1) and (1, -1)/sqrt 2, and
//        up e_z and e_x: n_v . n_A > 0 and n_v . n_B > 0 at every smooth
//        vertex, n_v = +-radial with the sign of the faces, and the frames
//        follow the field's oriented normal
//
//   seas_test_vertex_normals              (1, 2 or 4 ranks)
//   seas_test_vertex_normals --refuse N   one refusal, must exit non-zero:
//     1 a Y-shaped fault (a vertex with three fault faces: "junction"),
//     2 a 3-D mesh ("2-D meshes only"), 3 vertex_normal_max_turn_deg = 0
//     ("must be > 0")
//
//   make -f tests/2d_inplane/2d_inplane.mk test-vertex-normals
//
// Scope (Chunhui Zhao, 2026-11-01): 2-D triangle meshes, p = 1 geometry.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"
#include "frame_box.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
namespace ff = mfem::seas::fault_frame;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kLambda = 32.04e9, kMu = 32.04e9;
static const real_t kR = 10e3;            // circle radius [m]
static const real_t kDphi = 2.0;          // cell angle [deg]
static const int kNphi = 24, kNr = 6;     // cells in phi and r
static const real_t kFaultPhi = 12.0;     // the fault: |phi| <= 12 deg

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static bool AllTrue(bool b, MPI_Comm c) { return ISum(b ? 0 : 1, c) == 0; }

// ------------------------------------------------------------------ meshes
// The annulus sector of the file comment; straight = true: the same grid with
// x = s (arc length R phi), y = r - R, so the fault is the line y = 0.
static Mesh *ArcMesh(bool straight = false)
{
   const int nvx = kNphi + 1, nvy = kNr + 1;
   Mesh *m = new Mesh(2, nvx * nvy, 2 * kNphi * kNr, 0, 2);
   const real_t rad = M_PI / 180.0;
   for (int j = 0; j < nvy; j++)
      for (int i = 0; i < nvx; i++)
      {
         const real_t phi = (-0.5 * kNphi + i) * kDphi * rad;
         const real_t r = kR - 3e3 + 6e3 * j / kNr;
         if (straight) { m->AddVertex(kR * phi, r - kR); }
         else { m->AddVertex(r * std::sin(phi), kR - r * std::cos(phi)); }
      }
   auto v = [&](int i, int j) { return i + nvx * j; };
   for (int j = 0; j < kNr; j++)
      for (int i = 0; i < kNphi; i++)
      {
         m->AddTriangle(v(i, j), v(i + 1, j), v(i + 1, j + 1), 1);
         m->AddTriangle(v(i, j), v(i + 1, j + 1), v(i, j + 1), 1);
      }
   for (int i = 0; i < kNphi; i++)
   {
      m->AddBdrSegment(v(i, 0), v(i + 1, 0), 5);
      m->AddBdrSegment(v(i, kNr), v(i + 1, kNr), 5);
      const real_t pm = (-0.5 * kNphi + i + 0.5) * kDphi;
      if (std::abs(pm) < kFaultPhi) { m->AddBdrSegment(v(i, kNr / 2), v(i + 1, kNr / 2), 3); }
   }
   for (int j = 0; j < kNr; j++)
   {
      m->AddBdrSegment(v(0, j), v(0, j + 1), 5);
      m->AddBdrSegment(v(kNphi, j), v(kNphi, j + 1), 5);
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

// A Y-shaped fault (test_2d_fault_layout.cpp's YMesh): y = 0 for |x| < L/2
// and the branch y = x for 0 < x < L/2, all attr 3; three faces at (0, 0).
static Mesh *YMesh()
{
   const int n = 8;
   const real_t L = 10e3;
   Mesh *m = new Mesh(2, (n + 1) * (n + 1), 2 * n * n, 0, 2);
   for (int j = 0; j <= n; j++)
      for (int i = 0; i <= n; i++) { m->AddVertex(-L + 2.0 * L * i / n, -L + 2.0 * L * j / n); }
   auto v = [&](int i, int j) { return i + (n + 1) * j; };
   for (int j = 0; j < n; j++)
      for (int i = 0; i < n; i++)
      {
         m->AddTriangle(v(i, j), v(i + 1, j), v(i + 1, j + 1), 1);
         m->AddTriangle(v(i, j), v(i + 1, j + 1), v(i, j + 1), 1);
      }
   for (int k = 0; k < n; k++)
   {
      m->AddBdrSegment(v(k, 0), v(k + 1, 0), 5);
      m->AddBdrSegment(v(k, n), v(k + 1, n), 5);
      m->AddBdrSegment(v(0, k), v(0, k + 1), 5);
      m->AddBdrSegment(v(n, k), v(n, k + 1), 5);
   }
   const int c = n / 2;
   for (int i = c - n / 4; i < c + n / 4; i++) { m->AddBdrSegment(v(i, c), v(i + 1, c), 3); }
   for (int k = 0; k < n / 4; k++) { m->AddBdrSegment(v(c + k, c + k), v(c + k + 1, c + k + 1), 3); }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

// ---------------------------------------------------------------- operator
static FaultFrameSource Source(const real_t ref[3], const real_t up[3], bool averaged,
                               real_t max_turn = 5.0)
{
   auto c = std::make_shared<FaultFrameConfig>();
   for (int d = 0; d < 3; d++)
   {
      c->orientation.ref_normal[d] = ref[d];
      c->orientation.up[d] = up[d];
      c->up[d] = up[d];
   }
   c->orientation.fault_attrs = {3};
   c->averaged_vertex_normals = averaged;
   c->vertex_normal_max_turn_deg = max_turn;
   FaultFrameSource s;
   s.config = c;
   return s;
}

using Op = ElasticityDomainOperator<ParMesh>;

static std::unique_ptr<Op> MakeOp(ParMesh &pm, FaultFrameSource src,
                                  std::set<int> dirichlet = {5})
{
   LinearElastic *model = new LinearElastic(kLambda, kMu);   // leaked: test lifetime
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = dirichlet;
   const int dim = pm.Dimension();
   bc.default_dirichlet_func = [dim](const Vector &, real_t, Vector &u) { u.SetSize(dim); u = 0.0; };
   DomainConfig dc;
   dc.fault_frames = src;
   return std::make_unique<Op>(pm, 1, *model, 0.0, 1e3, 1e3, bc, DGMethod::IP,
                               SolverType::MUMPS, dc);
}

static real_t Dist3(const real_t a[3], const real_t b[3])
{
   return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
                    + (a[2] - b[2]) * (a[2] - b[2]));
}
static bool Same3(const real_t a[3], const real_t b[3])
{
   return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}
static bool SameFrame(const ff::Frame &a, const ff::Frame &b)
{
   return Same3(a.n, b.n) && Same3(a.d, b.d) && Same3(a.s, b.s) && a.fallback == b.fallback;
}
// largest deviation of (n, d, s) from an orthonormal triad with d x s = n
static real_t FrameError(const ff::Frame &f)
{
   real_t e = 0.0, c[3];
   const real_t *v[3] = {f.n, f.d, f.s};
   for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++)
      {
         e = std::max(e, std::abs(ff::Dot(v[a], v[b]) - (a == b ? 1.0 : 0.0)));
      }
   ff::Cross(f.d, f.s, c);
   return std::max(e, Dist3(c, f.n));
}

// the circle's unit radial vector at the point x on r = R (away from the
// centre (0, R))
static void Radial(const real_t x[3], real_t out[3])
{
   const real_t dx = x[0], dy = x[1] - kR, l = std::sqrt(dx * dx + dy * dy);
   out[0] = dx / l;
   out[1] = dy / l;
   out[2] = 0.0;
}

// ------------------------------------------------------------------- T2, T3
static void TestArc(MPI_Comm comm, const real_t ref[3], const real_t up[3], const char *tag,
                    bool digest, bool cut = false)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== T2/T3 inscribed polygon, " << tag << (cut ? ", cut along the fault" : "")
                << " ===\n";
   }
   std::unique_ptr<Mesh> s(ArcMesh());
   int np = 1;
   MPI_Comm_size(comm, &np);
   std::vector<int> part(s->GetNE(), 0);
   if (cut && np > 1)
   {
      for (int e = 0; e < s->GetNE(); e++)
      {
         Array<int> ev;
         s->GetElementVertices(e, ev);
         real_t cx = 0.0, cy = 0.0;
         for (int a = 0; a < ev.Size(); a++)
         {
            cx += s->GetVertex(ev[a])[0] / ev.Size();
            cy += s->GetVertex(ev[a])[1] / ev.Size();
         }
         const real_t r = std::sqrt(cx * cx + (cy - kR) * (cy - kR));
         part[e] = ((r > kR) ? 1 : 0) + ((np >= 4 && cx > 0.0) ? 2 : 0);
      }
   }
   ParMesh pm(comm, *s, cut ? part.data() : nullptr);
   if (cut)
   {
      int sh = ISum(pm.GetNSharedFaces(), comm);
      if (g_rank == 0) { std::cout << "  cut partition: " << sh << " shared faces in all\n"; }
   }
   auto op = MakeOp(pm, Source(ref, up, true));
   const FaultFrameField &F = *op->FrameField();
   const FaultFrameReport &R = F.Report();
   TEST_ASSERT(F.AveragedVertexNormals(), tag << ": the option is on");
   TEST_ASSERT(R.smooth_vertices == 11 && R.kink_vertices == 0 && R.tip_vertices == 2
               && std::abs(R.max_turn_smooth_deg - kDphi) < 1e-9,
               tag << ": T2a 11 smooth vertices (largest turn " << R.max_turn_smooth_deg
               << " deg), 0 kinks, 2 tips (" << R.smooth_vertices << ", " << R.kink_vertices
               << ", " << R.tip_vertices << ")");
   if (g_rank == 0) { F.Print(std::cout); }

   real_t e_rad = 0.0, e_frame = 0.0, e_vert = 0.0, e_mid = 0.0;
   real_t min_side = 1.0;
   bool tips = true, slots = true, kinds = true, sign_ok = true;
   const int inplane = F.InPlaneSlot();
   for (int fi = 0; fi < F.NumFaultFaces(); fi++)
   {
      const real_t *nf = F.Normal(fi);
      // the faces' orientation relative to the outward radial direction
      real_t mid[3], rm[3];
      for (int d = 0; d < 3; d++) { mid[d] = 0.5 * (F.VertexX(fi, 0)[d] + F.VertexX(fi, 1)[d]); }
      Radial(mid, rm);
      const real_t sgn = (ff::Dot(nf, rm) > 0.0) ? 1.0 : -1.0;
      for (int a = 0; a < 2; a++)
      {
         const real_t *x = F.VertexX(fi, a), *nv = F.VertexNormal(fi, a);
         const int k = F.VertexKind(fi, a);
         if (k == 0)
         {
            if (!Same3(nv, nf)) { tips = false; }
            continue;
         }
         if (k != 1 || std::abs(F.VertexTurnDeg(fi, a) - kDphi) > 1e-9) { kinds = false; }
         real_t r[3];
         Radial(x, r);
         for (int d = 0; d < 3; d++) { r[d] *= sgn; }
         e_rad = std::max(e_rad, Dist3(nv, r));
         // T3: the vertex normal points to the side of the face normal
         min_side = std::min(min_side, ff::Dot(nv, nf));
         ff::Frame fr;
         F.FrameAt(fi, x, fr);
         e_vert = std::max(e_vert, Dist3(fr.n, nv));
         e_frame = std::max(e_frame, FrameError(fr));
      }
      ff::Frame fm;
      F.FrameAt(fi, mid, fm);
      real_t want[3], wl = 0.0;
      for (int d = 0; d < 3; d++)
      {
         want[d] = 0.5 * (F.VertexNormal(fi, 0)[d] + F.VertexNormal(fi, 1)[d]);
         wl += want[d] * want[d];
      }
      for (int d = 0; d < 3; d++) { want[d] /= std::sqrt(wl); }
      e_mid = std::max(e_mid, Dist3(fm.n, want));
      e_frame = std::max(e_frame, FrameError(fm));
      if (ff::Dot(fm.n, nf) <= 0.0) { sign_ok = false; }
      // the slot rule: the out-of-plane vector has exactly zero in-plane
      // components (unit to rounding), the in-plane tangent and n exactly
      // zero z components, so the two slots never mix
      const real_t *oop = (inplane == 1) ? fm.d : fm.s, *inp = (inplane == 1) ? fm.s : fm.d;
      if (!(oop[0] == 0.0 && oop[1] == 0.0 && std::abs(std::abs(oop[2]) - 1.0) < 1e-15
            && inp[2] == 0.0 && fm.n[2] == 0.0)) { slots = false; }
   }
   e_rad = RMax(e_rad, comm);
   e_frame = RMax(e_frame, comm);
   e_vert = RMax(e_vert, comm);
   e_mid = RMax(e_mid, comm);
   min_side = -RMax(-min_side, comm);
   TEST_ASSERT(e_rad < 1e-13, tag << ": T2b every smooth vertex normal is the radial unit "
               "vector, with the faces' sign (largest error " << e_rad << ")");
   TEST_ASSERT(AllTrue(tips, comm), tag << ": T2b the two tips keep their face normal bit for bit");
   TEST_ASSERT(AllTrue(kinds, comm), tag << ": T2b every interior vertex smooth, turn 2 deg");
   TEST_ASSERT(e_vert < 1e-14 && e_mid < 1e-14 && e_frame < 1e-14,
               tag << ": T2c FrameAt: n = n_v at the vertices (" << e_vert << "), the "
               "normalized mean at the midpoints (" << e_mid << "), orthonormal with d x s = n ("
               << e_frame << ")");
   TEST_ASSERT(AllTrue(slots, comm), tag << ": T2c the out-of-plane slot is exactly +-e_z "
               "(slot " << inplane << ")");
   TEST_ASSERT(min_side > std::cos(kDphi * M_PI / 180.0) - 1e-12 && AllTrue(sign_ok, comm),
               tag << ": T3 n_v . n_face = cos(turn / 2) > 0 on both faces (smallest "
               << min_side << "), the midpoint frame on the face normal's side");

   if (digest)
   {
      // T2g: every vertex normal, sorted by position, in hex
      std::vector<double> loc;
      for (int fi = 0; fi < F.NumFaultFaces(); fi++)
         for (int a = 0; a < 2; a++)
         {
            for (int d = 0; d < 2; d++) { loc.push_back(F.VertexX(fi, a)[d]); }
            for (int d = 0; d < 2; d++) { loc.push_back(F.VertexNormal(fi, a)[d]); }
         }
      int nl = static_cast<int>(loc.size()), np = 0;
      MPI_Comm_size(comm, &np);
      std::vector<int> cnt(np), dsp(np, 0);
      MPI_Gather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm);
      if (g_rank == 0) { for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; } }
      std::vector<double> all(g_rank == 0 ? dsp[np - 1] + cnt[np - 1] : 0);
      MPI_Gatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(), MPI_DOUBLE, 0,
                  comm);
      if (g_rank == 0)
      {
         std::set<std::array<double, 4>> rows;
         for (size_t o = 0; o + 4 <= all.size(); o += 4)
         {
            rows.insert({all[o], all[o + 1], all[o + 2], all[o + 3]});
         }
         std::cout << "VN-DIGEST " << tag;
         for (const auto &r : rows) { std::printf(" %a,%a", r[2], r[3]); }
         std::cout << "\n";
      }
   }
}

static void TestThreshold(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== T2d threshold 1 deg: every vertex a kink ===\n"; }
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   std::unique_ptr<Mesh> s(ArcMesh());
   ParMesh pm(comm, *s);
   auto op = MakeOp(pm, Source(ref, up, true, 1.0));
   const FaultFrameField &F = *op->FrameField();
   const FaultFrameReport &R = F.Report();
   TEST_ASSERT(R.smooth_vertices == 0 && R.kink_vertices == 11 && R.tip_vertices == 2
               && std::abs(R.max_turn_kink_deg - kDphi) < 1e-9,
               "T2d 0 smooth, 11 kinks kept (largest turn " << R.max_turn_kink_deg << " deg), 2 tips");
   bool same = true;
   for (int fi = 0; fi < F.NumFaultFaces(); fi++)
   {
      for (int a = 0; a < 2; a++)
      {
         if (!Same3(F.VertexNormal(fi, a), F.Normal(fi))) { same = false; }
         ff::Frame fr;
         F.FrameAt(fi, F.VertexX(fi, a), fr);
         if (!SameFrame(fr, F.Frame(fi))) { same = false; }
      }
   }
   TEST_ASSERT(AllTrue(same, comm), "T2d every VertexNormal and FrameAt is the face's own, bit for bit");
}

static void TestOff(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== T2e the option off (\"face\") ===\n"; }
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   std::unique_ptr<Mesh> s(ArcMesh());
   ParMesh pm(comm, *s);
   auto op = MakeOp(pm, Source(ref, up, false));
   const FaultFrameField &F = *op->FrameField();
   const FaultFrameReport &R = F.Report();
   TEST_ASSERT(!F.AveragedVertexNormals() && R.smooth_vertices == 0 && R.kink_vertices == 0
               && R.tip_vertices == 0, "T2e off: no vertex counts");
   bool same = true;
   for (int fi = 0; fi < F.NumFaultFaces(); fi++)
   {
      real_t mid[3];
      for (int d = 0; d < 3; d++) { mid[d] = 0.5 * (F.VertexX(fi, 0)[d] + F.VertexX(fi, 1)[d]); }
      ff::Frame fr;
      F.FrameAt(fi, mid, fr);
      if (!SameFrame(fr, F.Frame(fi))) { same = false; }
      for (int a = 0; a < 2; a++) { if (!Same3(F.VertexNormal(fi, a), F.Normal(fi))) { same = false; } }
   }
   TEST_ASSERT(AllTrue(same, comm), "T2e off: FrameAt = Frame(fi) and VertexNormal = Normal(fi), bit for bit");
}

static void TestStraight(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== T2f a straight fault with the option on ===\n"; }
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   std::unique_ptr<Mesh> s(ArcMesh(true));
   ParMesh pm(comm, *s);
   auto op = MakeOp(pm, Source(ref, up, true));
   const FaultFrameField &F = *op->FrameField();
   const FaultFrameReport &R = F.Report();
   TEST_ASSERT(R.smooth_vertices == 11 && R.kink_vertices == 0 && R.tip_vertices == 2
               && R.max_turn_smooth_deg == 0.0,
               "T2f 11 smooth vertices with turn 0, 2 tips (" << R.smooth_vertices << ", "
               << R.max_turn_smooth_deg << ")");
   bool same = true;
   for (int fi = 0; fi < F.NumFaultFaces(); fi++)
   {
      real_t mid[3];
      for (int d = 0; d < 3; d++) { mid[d] = 0.5 * (F.VertexX(fi, 0)[d] + F.VertexX(fi, 1)[d]); }
      for (const real_t *x : {F.VertexX(fi, 0), F.VertexX(fi, 1), static_cast<const real_t *>(mid)})
      {
         ff::Frame fr;
         F.FrameAt(fi, x, fr);
         if (!SameFrame(fr, F.Frame(fi))) { same = false; }
      }
      for (int a = 0; a < 2; a++) { if (!Same3(F.VertexNormal(fi, a), F.Normal(fi))) { same = false; } }
   }
   TEST_ASSERT(AllTrue(same, comm), "T2f every VertexNormal and FrameAt is the face's own, bit for bit");
}

// ---------------------------------------------------------------- refusals
static int Refuse(int which, MPI_Comm comm)
{
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   if (which == 1)
   {
      std::unique_ptr<Mesh> s(YMesh());
      ParMesh pm(comm, *s);
      auto op = MakeOp(pm, Source(ref, up, true));
      (void)op;
   }
   else if (which == 2)
   {
      std::unique_ptr<Mesh> s(frame_box::FrameBox(0, frame_box::Rot()));
      ParMesh pm(comm, *s);
      auto op = MakeOp(pm, Source(ref, up, true), {5, 6});   // 6: the band, as in test_fault_frame_field
      (void)op;
   }
   else if (which == 3)
   {
      std::unique_ptr<Mesh> s(ArcMesh());
      ParMesh pm(comm, *s);
      auto op = MakeOp(pm, Source(ref, up, true, 0.0));
      (void)op;
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
   const int np = Mpi::WorldSize();
   if (argc == 3 && std::strcmp(argv[1], "--refuse") == 0)
   {
      // exit 0 only when the refusal did NOT fire; the make target expects
      // a non-zero exit with the refusal's message
      return Refuse(std::atoi(argv[2]), comm) == 0 ? 0 : 3;
   }
   if (g_rank == 0)
   {
      std::cout << "test_vertex_normals (vertex-normal plan rev1, phase P1) on " << np
                << " rank(s)\n";
   }
   const real_t s2 = 1.0 / std::sqrt(2.0);
   const real_t ez[3] = {0.0, 0.0, 1.0}, ex[3] = {1.0, 0.0, 0.0};
   const real_t rdown[3] = {0.0, -1.0, 0.0}, rup[3] = {0.0, 1.0, 0.0}, rtilt[3] = {s2, -s2, 0.0};
   TestArc(comm, rdown, ez, "ref (0,-1) up e_z", true);
   TestArc(comm, rdown, ez, "ref (0,-1) up e_z", true, true);
   TestArc(comm, rup, ez, "ref (0,1) up e_z", false);
   TestArc(comm, rtilt, ez, "ref (1,-1)/sqrt2 up e_z", false);
   TestArc(comm, rdown, ex, "ref (0,-1) up e_x", false);
   TestThreshold(comm);
   TestOff(comm);
   TestStraight(comm);
   if (g_rank == 0)
   {
      std::cout << "\n" << num_tests - num_failed << " / " << num_tests << " tests passed\n";
   }
   return num_failed == 0 ? 0 : 1;
}
