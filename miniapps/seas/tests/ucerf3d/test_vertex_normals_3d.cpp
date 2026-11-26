// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/ucerf3d/test_vertex_normals_3d.cpp — UCERF 3-D plan rev5, P1.1
// (Chunhui Zhao, 2026-11-25: "start P1"): the averaged vertex normals of
// FaultFrameField on 3-D fault surfaces ([fault_geometry] vertex_normals =
// "averaged"; fault/fault_frame_field.hpp, (1') and (2')).
//
//   U-N1  exact on the planar test meshes of plan section 10.1
//         (ucerf3d/mesh/T-plane.msh, T-dip.msh; fault attr 3, creep attr 6
//         as the interior Dirichlet faces): no kinks, every vertex normal and
//         every FrameAt normal is the face normal (1e-14) and the analytic
//         plane normal (1e-12); T-bend.msh: the vertices on the 20 deg bend
//         line are kinks (turn 20 deg) and keep the face normals, every
//         other vertex keeps its plane's normal (1e-14)
//   U-N2  h^2 on a cylinder (radius R = 10 km, axis e_z; a graded
//         structured grid in (phi, z), n = 16, 32, 64, 128 cells (n = 8 turns
//         5.9 deg at a vertex, above the 5 deg threshold), each cell split
//         into 6 Kuhn tets; the
//         fault is the patch r = R, |phi| <= 7.5 deg, |z| <= L/4): on 4
//         levels the largest angle between the exact radial normal and
//         (a) n_v at the smooth vertices, (b) FrameAt's normal at points of
//         the faces whose vertices are all smooth; both must fall as h^2
//         (rate >= 1.8 between the two finest levels), and face mode's error
//         (the face normal at the same points) is printed: it falls as h
//   U-N1c on the cylinder: frames orthonormal with d x s = n (1e-14), n_v . n_f > 0
//         over every fan, the two faces of an interior edge give the same
//         FrameAt bits at its midpoint, and the operator's per-quadrature-
//         point and per-DOF triads are eps_f FrameAt bit for bit (p = 1, 2)
//   U-N1d partition invariance: rank 0 prints "VN3-DIGEST" (a hash of the
//         FaultVertexRows table, which is sorted by key) for the cylinder and
//         T-dip; the make target checks that 1, 2 and 4 ranks print the same
//
//   seas_test_vertex_normals_3d              (1, 2 or 4 ranks)
//   seas_test_vertex_normals_3d --refuse N   one refusal, must exit non-zero:
//     1 a T-shaped fault (an edge with three fault faces: "junction"),
//     2 vertex_normal_max_turn_deg = 0 ("must be > 0")
//
//   make -f tests/ucerf3d/ucerf3d.mk test-vertex-normals-3d
//
// Scope (Chunhui Zhao, 2026-11-25): tetrahedral meshes, p = 1 geometry.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"
#include "../../cycle/cycle_frames.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
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
static const real_t kR = 10e3;                  // cylinder radius [m]
static const real_t kPhi = 15.0 * M_PI / 180;   // mesh half angle
static const real_t kL = 2.0 * kPhi * kR;       // mesh height [m]
static const real_t kDr = 1e3;                  // radial cell size [m]

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static bool AllTrue(bool b, MPI_Comm c) { return ISum(b ? 0 : 1, c) == 0; }

static real_t Dist3(const real_t a[3], const real_t b[3])
{
   return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
                    + (a[2] - b[2]) * (a[2] - b[2]));
}
static real_t TriadError(const real_t n[3], const real_t d[3], const real_t s[3])
{
   real_t e = 0.0, c[3];
   const real_t *v[3] = {n, d, s};
   for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++)
      {
         e = std::max(e, std::abs(ff::Dot(v[a], v[b]) - (a == b ? 1.0 : 0.0)));
      }
   ff::Cross(d, s, c);
   return std::max(e, Dist3(c, n));
}
/// The angle [rad] between the lines of two unit vectors (sign-free).
static real_t LineAngle(const real_t a[3], const real_t b[3])
{
   real_t c[3];
   ff::Cross(a, b, c);
   return std::atan2(std::sqrt(ff::Dot(c, c)), std::abs(ff::Dot(a, b)));
}

// ------------------------------------------------------------------ meshes
// A structured grid of ni x nj x nk cells, vertex (i, j, k) at X(i, j, k),
// each cell split into the 6 Kuhn tets (paths from corner 000 to 111), so
// that the faces j = const are split along the diagonal (i, k)-(i+1, k+1)
// and the faces i = const along (j, k)-(j+1, k+1).  fault(i, j, k, axis)
// says whether the cell face of the vertex (i, j, k) normal to axis (0: i,
// 1: j) is a fault face (attr 3).  Outer faces: attr 5.
template <class XF, class FF>
static Mesh *KuhnMesh(int ni, int nj, int nk, XF X, FF fault)
{
   const int nv = (ni + 1) * (nj + 1) * (nk + 1);
   Mesh *m = new Mesh(3, nv, 6 * ni * nj * nk, 0, 3);
   auto v = [&](int i, int j, int k) { return i + (ni + 1) * (j + (nj + 1) * k); };
   for (int k = 0; k <= nk; k++)
      for (int j = 0; j <= nj; j++)
         for (int i = 0; i <= ni; i++)
         {
            real_t x[3];
            X(i, j, k, x);
            m->AddVertex(x[0], x[1], x[2]);
         }
   const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
   for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
         for (int i = 0; i < ni; i++)
            for (int p = 0; p < 6; p++)
            {
               int c[3] = {i, j, k}, t[4];
               t[0] = v(c[0], c[1], c[2]);
               for (int s = 0; s < 3; s++)
               {
                  c[perm[p][s]]++;
                  t[s + 1] = v(c[0], c[1], c[2]);
               }
               m->AddTet(t[0], t[1], t[2], t[3], 1);
            }
   // the cell face of vertex (i, j, k) normal to axis, split like the tets:
   // (V, V + e_p, V + e_p + e_q) and (V, V + e_q, V + e_p + e_q), p < q
   auto quad = [&](int axis, int i, int j, int k, int attr)
   {
      const int p = (axis == 0) ? 1 : 0, q = (axis == 2) ? 1 : 2;
      int c0[3] = {i, j, k}, c1[3] = {i, j, k}, c2[3] = {i, j, k}, c3[3] = {i, j, k};
      c1[p]++; c2[q]++; c3[p]++; c3[q]++;
      const int v0 = v(c0[0], c0[1], c0[2]), v1 = v(c1[0], c1[1], c1[2]),
                v2 = v(c2[0], c2[1], c2[2]), v3 = v(c3[0], c3[1], c3[2]);
      m->AddBdrTriangle(v0, v1, v3, attr);
      m->AddBdrTriangle(v0, v2, v3, attr);
   };
   const int n[3] = {ni, nj, nk};
   for (int axis = 0; axis < 3; axis++)
   {
      const int p = (axis == 0) ? 1 : 0, q = (axis == 2) ? 1 : 2;
      for (int a = 0; a <= n[axis]; a++)
         for (int b = 0; b < n[p]; b++)
            for (int c = 0; c < n[q]; c++)
            {
               int idx[3];
               idx[axis] = a; idx[p] = b; idx[q] = c;
               if (a == 0 || a == n[axis]) { quad(axis, idx[0], idx[1], idx[2], 5); }
               else if (axis < 2 && fault(idx[0], idx[1], idx[2], axis))
               {
                  quad(axis, idx[0], idx[1], idx[2], 3);
               }
            }
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

// The graded parameter of grid line i of n: s in [-1, 1], smooth, monotone.
static real_t Graded(int i, int n)
{
   const real_t s = -1.0 + 2.0 * i / n;
   return s + 0.2 * std::sin(M_PI * s);
}

// U-N2's cylinder at level n: n cells in phi and z, 2 in r (r = R - Dr, R,
// R + Dr); the fault patch |phi| <= Phi / 2, |z| <= L / 4 on r = R.
static Mesh *CylMesh(int n)
{
   auto X = [n](int i, int j, int k, real_t x[3])
   {
      const real_t phi = kPhi * Graded(i, n), r = kR + (j - 1) * kDr;
      x[0] = r * std::sin(phi);
      x[1] = kR - r * std::cos(phi);
      x[2] = 0.5 * kL * Graded(k, n);
   };
   auto F = [n](int i, int j, int k, int axis)
   {
      return axis == 1 && j == 1 && i >= n / 4 && i < 3 * n / 4 && k >= n / 4 && k < 3 * n / 4;
   };
   return KuhnMesh(n, 2, n, X, F);
}

// Refusal 1: a T-shaped fault in a 8^3 box: the plane j = 4 for 2 <= i < 6
// and the half plane i = 4 for j >= 4 (two cells), 2 <= k < 6: the edge
// (4, 4, k) has three fault faces.
static Mesh *TMesh()
{
   const real_t h = 1e3;
   auto X = [h](int i, int j, int k, real_t x[3]) { x[0] = i * h; x[1] = j * h; x[2] = k * h; };
   auto F = [](int i, int j, int k, int axis)
   {
      if (k < 2 || k >= 6) { return false; }
      if (axis == 1) { return j == 4 && i >= 2 && i < 6; }
      return i == 4 && j >= 4 && j < 6;
   };
   return KuhnMesh(8, 8, 8, X, F);
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

static std::unique_ptr<Op> MakeOp(ParMesh &pm, FaultFrameSource src, int order = 1,
                                  std::set<int> dirichlet = {5})
{
   LinearElastic *model = new LinearElastic(kLambda, kMu);   // leaked: test lifetime
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = dirichlet;
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   DomainConfig dc;
   dc.fault_frames = src;
   return std::make_unique<Op>(pm, order, *model, 0.0, 1e3, 1e3, bc, DGMethod::IP,
                               SolverType::MUMPS, dc);
}

static FaceElementTransformations *FaceTr(ParMesh &pm, const Op &op, int fi)
{
   const int n_int = op.GetFaultInteriorFaces().Size();
   return (fi < n_int) ? pm.GetInteriorFaceTransformations(op.GetFaultInteriorFaces()[fi])
                       : pm.GetSharedFaceTransformations(op.GetFaultSharedFaces()[fi - n_int]);
}

/// Sample points of a triangle (barycentric): the vertices, the edge
/// midpoints, the centroid and a degree-4 rule's points.
static std::vector<std::array<real_t, 3>> SampleBary()
{
   std::vector<std::array<real_t, 3>> b = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0},
                                           {0, 0.5, 0.5}, {0.5, 0, 0.5}, {1. / 3, 1. / 3, 1. / 3}};
   const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 4);
   for (int q = 0; q < ir.GetNPoints(); q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      b.push_back({1.0 - ip.x - ip.y, ip.x, ip.y});
   }
   return b;
}
static void BaryPoint(const FaultFrameField &F, int fi, const std::array<real_t, 3> &l, real_t x[3])
{
   for (int d = 0; d < 3; d++)
   {
      x[d] = l[0] * F.VertexX(fi, 0)[d] + l[1] * F.VertexX(fi, 1)[d] + l[2] * F.VertexX(fi, 2)[d];
   }
}

/// FNV-1a over the FaultVertexRows table (rank 0; "" elsewhere).
static std::string Digest(const FaultFrameField &F, MPI_Comm comm)
{
   Vector rows;
   cycle::FaultVertexRows(F, comm, rows);
   if (g_rank != 0) { return ""; }
   uint64_t h = 1469598103934665603ULL;
   for (int i = 0; i < rows.Size(); i++)
   {
      unsigned char b[sizeof(double)];
      const double v = rows(i);
      std::memcpy(b, &v, sizeof(double));
      for (unsigned char c : b) { h ^= c; h *= 1099511628211ULL; }
   }
   char s[64];
   std::snprintf(s, sizeof(s), "%016llx rows %d", static_cast<unsigned long long>(h),
                 rows.Size() / 9);
   return s;
}

// ---------------------------------------------------------------- U-N1
static void TestPlane(const std::string &name, const real_t n_exact[3], MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== U-N1 " << name << " (planar) ===\n"; }
   Mesh smesh(("ucerf3d/mesh/" + name + ".msh").c_str(), 1, 1);
   ParMesh pm(comm, smesh);
   smesh.Clear();
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   auto op = MakeOp(pm, Source(ref, up, true), 1, {5, 6});
   const FaultFrameField &F = *op->FrameField();
   const FaultFrameReport &r = F.Report();
   const int nf = F.NumFaultFaces();
   real_t e_vn = 0.0, e_fr = 0.0, e_an = 0.0;
   const auto bary = SampleBary();
   int varies = 0;
   for (int fi = 0; fi < nf; fi++)
   {
      if (F.FrameVaries(fi)) { varies++; }
      const real_t *nf_ = F.Normal(fi);
      e_an = std::max(e_an, LineAngle(nf_, n_exact));
      for (int a = 0; a < 3; a++) { e_vn = std::max(e_vn, Dist3(F.VertexNormal(fi, a), nf_)); }
      const ff::Frame &f0 = F.Frame(fi);
      for (const auto &l : bary)
      {
         real_t x[3];
         BaryPoint(F, fi, l, x);
         ff::Frame fr;
         F.FrameAt(fi, x, fr);
         e_fr = std::max({e_fr, Dist3(fr.n, f0.n), Dist3(fr.d, f0.d), Dist3(fr.s, f0.s)});
      }
   }
   e_vn = RMax(e_vn, comm); e_fr = RMax(e_fr, comm); e_an = RMax(e_an, comm);
   varies = ISum(varies, comm);
   if (g_rank == 0)
   {
      F.Print(std::cout);
      std::printf("  %s: faces %d, vary %d; max |n_v - n_f| %.3e, max FrameAt - Frame %.3e, "
                  "max angle(n_f, exact) %.3e rad\n", name.c_str(), r.fault.faces, varies, e_vn,
                  e_fr, e_an);
   }
   TEST_ASSERT(r.kink_vertices == 0 && r.max_turn_smooth_deg < 1e-9,
               "U-N1 " << name << ": no kinks, largest turn " << r.max_turn_smooth_deg << " deg");
   TEST_ASSERT(r.smooth_vertices > 0 && r.edge_vertices > 0,
               "U-N1 " << name << ": " << r.smooth_vertices << " smooth, " << r.edge_vertices
               << " edge, " << r.tip_vertices << " tip vertices");
   TEST_ASSERT(e_vn <= 1e-14, "U-N1 " << name << ": every vertex normal is the face normal ("
               << e_vn << " <= 1e-14)");
   TEST_ASSERT(e_fr <= 1e-14, "U-N1 " << name << ": FrameAt (n, d, s) = Frame(fi) at 19 points "
               "per face (" << e_fr << " <= 1e-14)");
   TEST_ASSERT(e_an <= 1e-12, "U-N1 " << name << ": the face normals are the analytic plane "
               "normal (" << e_an << " rad <= 1e-12)");
   const std::string dg = Digest(F, comm);
   if (g_rank == 0) { std::cout << "VN3-DIGEST " << name << " " << dg << "\n"; }
}

static void TestBend(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== U-N1 T-bend (two planes, 20 deg bend) ===\n"; }
   Mesh smesh("ucerf3d/mesh/T-bend.msh", 1, 1);
   ParMesh pm(comm, smesh);
   smesh.Clear();
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   auto op = MakeOp(pm, Source(ref, up, true), 1, {5, 6});
   const FaultFrameField &F = *op->FrameField();
   const FaultFrameReport &r = F.Report();
   real_t e_vn = 0.0, off_line = 0.0;
   for (int fi = 0; fi < F.NumFaultFaces(); fi++)
      for (int a = 0; a < 3; a++)
      {
         const real_t *x = F.VertexX(fi, a);
         if (F.VertexKind(fi, a) == 2) { off_line = std::max(off_line, std::hypot(x[0], x[1])); }
         e_vn = std::max(e_vn, Dist3(F.VertexNormal(fi, a), F.Normal(fi)));
      }
   e_vn = RMax(e_vn, comm); off_line = RMax(off_line, comm);
   if (g_rank == 0) { F.Print(std::cout); }
   TEST_ASSERT(r.kink_vertices > 0 && std::abs(r.max_turn_kink_deg - 20.0) < 1e-6,
               "U-N1 T-bend: " << r.kink_vertices << " kinks, turn " << r.max_turn_kink_deg
               << " deg (20)");
   TEST_ASSERT(off_line < 1e-6, "U-N1 T-bend: every kink lies on the bend line x = y = 0 ("
               << off_line << " m)");
   TEST_ASSERT(r.max_turn_smooth_deg < 1e-9, "U-N1 T-bend: every other vertex has turn 0 ("
               << r.max_turn_smooth_deg << " deg)");
   TEST_ASSERT(e_vn <= 1e-14, "U-N1 T-bend: every vertex keeps its plane's normal (" << e_vn
               << " <= 1e-14)");
}

// ---------------------------------------------------------------- U-N2
struct CylErr { real_t h = 0, ev = 0, ef = 0, eface = 0, eedge = 0; int smooth = 0, edge = 0; };

static void Radial(const real_t x[3], real_t n[3])
{
   const real_t a = x[0], b = x[1] - kR, l = std::hypot(a, b);
   n[0] = a / l; n[1] = b / l; n[2] = 0.0;
}

static CylErr Cylinder(int n, MPI_Comm comm)
{
   std::unique_ptr<Mesh> sm(CylMesh(n));
   ParMesh pm(comm, *sm);
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   auto op = MakeOp(pm, Source(ref, up, true));
   const FaultFrameField &F = *op->FrameField();
   CylErr e;
   e.h = kL / n;
   e.smooth = F.Report().smooth_vertices;
   e.edge = F.Report().edge_vertices;
   const auto bary = SampleBary();
   for (int fi = 0; fi < F.NumFaultFaces(); fi++)
   {
      bool all_smooth = true;
      for (int a = 0; a < 3; a++)
      {
         real_t rn[3];
         Radial(F.VertexX(fi, a), rn);
         const int k = F.VertexKind(fi, a);
         if (k == 1) { e.ev = std::max(e.ev, LineAngle(F.VertexNormal(fi, a), rn)); }
         if (k == 3) { e.eedge = std::max(e.eedge, LineAngle(F.VertexNormal(fi, a), rn)); }
         all_smooth = all_smooth && (k == 1);
      }
      if (!all_smooth) { continue; }
      for (const auto &l : bary)
      {
         real_t x[3], rn[3];
         BaryPoint(F, fi, l, x);
         Radial(x, rn);
         ff::Frame fr;
         F.FrameAt(fi, x, fr);
         e.ef = std::max(e.ef, LineAngle(fr.n, rn));
         e.eface = std::max(e.eface, LineAngle(F.Normal(fi), rn));
      }
   }
   e.ev = RMax(e.ev, comm); e.ef = RMax(e.ef, comm);
   e.eface = RMax(e.eface, comm); e.eedge = RMax(e.eedge, comm);
   return e;
}

static void TestCylinderRates(MPI_Comm comm)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== U-N2 cylinder R = 10 km: angle to the radial normal [rad] ===\n"
                << "      n       h(m)  smooth  edge     n_v(smooth)  rate   FrameAt n  rate"
                   "   face mode  rate   n_v(edge)  rate\n";
   }
   const int levels[4] = {16, 32, 64, 128};
   std::vector<CylErr> E;
   for (int n : levels) { E.push_back(Cylinder(n, comm)); }
   auto rate = [](real_t a, real_t b) { return std::log2(a / b); };
   for (size_t i = 0; i < E.size(); i++)
   {
      if (g_rank != 0) { continue; }
      const CylErr &c = E[i];
      if (i == 0)
      {
         std::printf("  %5d %10.2f %7d %5d   %11.3e     -   %9.3e     -   %9.3e     -   %9.3e     -\n",
                     levels[i], c.h, c.smooth, c.edge, c.ev, c.ef, c.eface, c.eedge);
      }
      else
      {
         const CylErr &p = E[i - 1];
         std::printf("  %5d %10.2f %7d %5d   %11.3e  %4.2f   %9.3e  %4.2f   %9.3e  %4.2f   %9.3e  %4.2f\n",
                     levels[i], c.h, c.smooth, c.edge, c.ev, rate(p.ev, c.ev), c.ef,
                     rate(p.ef, c.ef), c.eface, rate(p.eface, c.eface), c.eedge,
                     rate(p.eedge, c.eedge));
      }
   }
   const CylErr &a = E[2], &b = E[3];
   TEST_ASSERT(rate(a.ev, b.ev) >= 1.8, "U-N2 n_v at the smooth vertices falls as h^2 (rate "
               << rate(a.ev, b.ev) << " >= 1.8)");
   TEST_ASSERT(rate(a.ef, b.ef) >= 1.8, "U-N2 FrameAt's normal inside the faces falls as h^2 "
               "(rate " << rate(a.ef, b.ef) << " >= 1.8)");
   TEST_ASSERT(b.ef < b.eface, "U-N2 averaged error below face mode's on the finest level ("
               << b.ef << " < " << b.eface << ")");
}

// ---------------------------------------------------------------- U-N1c
static void TestCylinderFrames(int order, MPI_Comm comm)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== U-N1c cylinder n = 16, p = " << order << ": frames and kernels ===\n";
   }
   std::unique_ptr<Mesh> sm(CylMesh(16));
   ParMesh pm(comm, *sm);
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   auto op = MakeOp(pm, Source(ref, up, true), order);
   const FaultFrameField &F = *op->FrameField();
   const FaultBasis &fb = *op->GetFaultBasis();
   const int nf = op->GetNumFaultFaces(), nbf = op->GetNbfPerFace();
   const IntegrationRule &nir = op->GetFaceQuadrature()->GetNodalRule();
   DenseMatrix B;
   op->GetFaultDOFBasis(B);
   bool qp_same = true, dof_same = true, fan_ok = true, edge_same = true;
   real_t orth = 0.0;
   int varies = 0;
   std::map<std::array<long long, 3>, std::array<real_t, 3>> mid;   // edge midpoint -> n
   Vector x;
   for (int fi = 0; fi < nf; fi++)
   {
      if (F.FrameVaries(fi)) { varies++; }
      for (int a = 0; a < 3; a++)
      {
         if (F.VertexKind(fi, a) == 1 || F.VertexKind(fi, a) == 3)
         {
            fan_ok = fan_ok && ff::Dot(F.VertexNormal(fi, a), F.Normal(fi)) > 0.0;
         }
      }
      for (int a = 0; a < 3; a++)
      {
         const real_t *p = F.VertexX(fi, a), *q = F.VertexX(fi, (a + 1) % 3);
         real_t m[3] = {0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2])};
         ff::Frame fr;
         F.FrameAt(fi, m, fr);
         const std::array<long long, 3> k{std::llround(m[0] * 1e3), std::llround(m[1] * 1e3),
                                          std::llround(m[2] * 1e3)};
         auto it = mid.find(k);
         if (it == mid.end()) { mid[k] = {fr.n[0], fr.n[1], fr.n[2]}; }
         else
         {
            edge_same = edge_same && it->second[0] == fr.n[0] && it->second[1] == fr.n[1]
                        && it->second[2] == fr.n[2];
         }
      }
      const FaultBasisData &b = fb.GetBasis(fi);
      const real_t eps = b.sign_flipped ? -1.0 : 1.0;
      FaceElementTransformations *Tr = FaceTr(pm, *op, fi);
      const IntegrationRule *fir = &IntRules.Get(Geometry::TRIANGLE, 2 * order + 1);
      if (F.FrameVaries(fi))
      {
         for (int q = 0; q < static_cast<int>(b.qp_data.size()); q++)
         {
            Tr->Face->Transform(fir->IntPoint(q), x);
            const real_t x3[3] = {x(0), x(1), x(2)};
            ff::Frame fr;
            F.FrameAt(fi, x3, fr);
            const FaultBasisQPData &qd = b.qp_data[q];
            for (int d = 0; d < 3; d++)
            {
               qp_same = qp_same && qd.normal[d] == eps * fr.n[d]
                         && qd.tangent1[d] == eps * fr.d[d] && qd.tangent2[d] == eps * fr.s[d];
            }
            orth = std::max(orth, TriadError(fr.n, fr.d, fr.s));
         }
      }
      for (int k = 0; k < nbf; k++)
      {
         Tr->Face->Transform(nir.IntPoint(k), x);
         const real_t x3[3] = {x(0), x(1), x(2)};
         ff::Frame fr;
         F.FrameAt(fi, x3, fr);
         orth = std::max(orth, TriadError(fr.n, fr.d, fr.s));
         const int col = fi * nbf + k;
         for (int d = 0; d < 3; d++)
         {
            dof_same = dof_same && B(d, col) == eps * fr.n[d] && B(3 + d, col) == eps * fr.d[d]
                       && B(6 + d, col) == eps * fr.s[d];
         }
      }
   }
   varies = ISum(varies, comm);
   orth = RMax(orth, comm);
   TEST_ASSERT(varies > 0, "U-N1c " << varies << " fault face copies carry varying frames");
   TEST_ASSERT(orth < 1e-14, "U-N1c every triad orthonormal with d x s = n (" << orth << ")");
   TEST_ASSERT(AllTrue(fan_ok, comm), "U-N1c n_v . n_f > 0 for every face of every fan");
   TEST_ASSERT(AllTrue(edge_same, comm), "U-N1c the faces of an interior edge give the same "
               "FrameAt bits at its midpoint (faces on one rank)");
   TEST_ASSERT(AllTrue(qp_same, comm), "U-N1c qp_data = eps_f FrameAt at the quadrature points, "
               "bit for bit");
   TEST_ASSERT(AllTrue(dof_same, comm), "U-N1c GetFaultDOFBasis = eps_f FrameAt at the nodal "
               "points, bit for bit");
   if (order == 1)
   {
      const std::string dg = Digest(F, comm);
      if (g_rank == 0) { std::cout << "VN3-DIGEST cylinder16 " << dg << "\n"; }
   }
}

// ---------------------------------------------------------------- refusals
static int Refuse(int which, MPI_Comm comm)
{
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   if (which == 1)
   {
      std::unique_ptr<Mesh> s(TMesh());
      ParMesh pm(comm, *s);
      // both sheets (normals e_y and e_x) pinned well by one reference normal
      const real_t ref_t[3] = {M_SQRT1_2, -M_SQRT1_2, 0.0};
      auto op = MakeOp(pm, Source(ref_t, up, true));
      (void)op;
   }
   else if (which == 2)
   {
      std::unique_ptr<Mesh> s(CylMesh(8));
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
   if (argc == 3 && std::strcmp(argv[1], "--refuse") == 0)
   {
      return Refuse(std::atoi(argv[2]), comm) == 0 ? 0 : 3;
   }
   if (g_rank == 0)
   {
      std::cout << "test_vertex_normals_3d (UCERF 3-D plan rev5, P1.1) on " << Mpi::WorldSize()
                << " rank(s)\n";
   }
   const real_t n_plane[3] = {0.0, 1.0, 0.0};
   const real_t dip = 58.0 * M_PI / 180.0;
   const real_t n_dip[3] = {0.0, std::sin(dip), std::cos(dip)};
   TestPlane("T-plane", n_plane, comm);
   TestPlane("T-dip", n_dip, comm);
   TestBend(comm);
   TestCylinderRates(comm);
   TestCylinderFrames(1, comm);
   TestCylinderFrames(2, comm);
   if (g_rank == 0)
   {
      std::cout << "\n" << (num_tests - num_failed) << "/" << num_tests << " tests passed\n";
   }
   return num_failed == 0 ? 0 : 1;
}
