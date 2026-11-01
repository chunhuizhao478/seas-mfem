// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_vertex_normal_kernels.cpp — vertex-normal plan rev1,
// phase P2 (Chunhui Zhao, 2026-11-01: "send this document to 2d thread and
// start implementation and test"; decision 2026-11-01: "Per-QP frame from
// eq. 2"; /mnt/project-files/2d-inplane/nonplanar-literature/
// vertex-normal-plan-rev1.pdf): the kernels of ElasticityDomainOperator with
// averaged vertex normals.  On a face whose frame varies, SetupVaryingFrames_
// (domain/elasticity_operator_setup.inl) puts the kernel triad
// eps_f FrameAt(x) at every quadrature point (FaultBasis qp_data) and at every
// face DOF (dof_triad, returned by GetFaultDOFBasis); the slip embedding and
// the traction projection read qp_data, unchanged.
//
// Mesh: test_vertex_normals.cpp's annulus sector (the same function): the
// fault is a polygon of 12 faces inscribed in the circle of radius
// R = 10 km around (0, R), a turn of 2 deg at each of its 11 interior
// vertices; straight = true maps it to y = 0.  Outer sides Dirichlet (u = 0).
//
//   F    frames (arc, averaged): every fault face varies; qp_data[q] =
//        eps_f FrameAt(x_q) bit for bit, nl and sign_flipped as the face
//        mode's; GetFaultDOFBasis = dof_triad = eps_f FrameAt(x_k); at every
//        interior vertex the oriented normals of all its DOFs (both faces,
//        every copy) are the same bit for bit and radial to 1e-13 (the
//        plan's T2); every triad orthonormal (1e-15)
//   T1   regression: a straight fault with the option on, and the arc with
//        threshold 1 deg (every vertex a kink), equal the option off bit for
//        bit (qp_data, GetFaultDOFBasis, ComputeTraction of a random u and
//        slip); the arc with the option off has no per-DOF triads
//   T4   adjoint: with u = 0 the traction is sigma_F times the projection of
//        the embedded slip, so S(a, b) = sum_faces a^T M tau(b) (M the face
//        mass matrix, nl_q-weighted) must be symmetric: |S(a, b) - S(b, a)| <=
//        1e-12 |S|, for the arc averaged and face mode; the slip gives no
//        normal traction at u = 0 (relative 1e-12); the per-face operators
//        built after the frames (ComputeTractionFast, AddSlipLoadFast) equal
//        ComputeTraction and AssembleSlipOnlyRHS (relative 1e-11)
//   S    static sawtooth (a small form of the plan's T5): a tapered slip
//        delta(s) = cos^2(pi s / 2 L) m along the fault (L its half length),
//        one static solve; the vertex jump in sigma_n (the spread of the
//        normal traction among the DOFs at one interior vertex) for the
//        straight fault, the arc in face mode and the arc averaged.  The
//        straight fault's jump at the same arc length is the discretization's
//        own (about 3.7 MPa here, h about 350 m); the arc's excess over it is
//        the sawtooth.  Face mode: the excess must exceed the straight jump;
//        averaged: below a tenth of face mode's excess
//
//   seas_test_vertex_normal_kernels       (1, 2 or 4 ranks)
//   make -f tests/2d_inplane/2d_inplane.mk test-vertex-normal-kernels
//
// Scope (Chunhui Zhao, 2026-11-01): 2-D triangle meshes, p = 1, IP.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <random>
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
static real_t RSum(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static bool AllTrue(bool b, MPI_Comm c) { return ISum(b ? 0 : 1, c) == 0; }

// ------------------------------------------------------------------ meshes
// test_vertex_normals.cpp's ArcMesh (the same function).
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

// ---------------------------------------------------------------- operator
static FaultFrameSource Source(bool averaged, real_t max_turn = 5.0)
{
   auto c = std::make_shared<FaultFrameConfig>();
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
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

static std::unique_ptr<Op> MakeOp(ParMesh &pm, FaultFrameSource src)
{
   LinearElastic *model = new LinearElastic(kLambda, kMu);   // leaked: test lifetime
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(2); u = 0.0; };
   DomainConfig dc;
   dc.fault_frames = src;
   return std::make_unique<Op>(pm, 1, *model, 0.0, 1e3, 1e3, bc, DGMethod::IP,
                               SolverType::MUMPS, dc);
}

static ParFiniteElementSpace *Space(Op &op)
{
   return const_cast<ParFiniteElementSpace *>(&op.GetFESpace());
}

static FaceElementTransformations *FaceTr(ParMesh &pm, const Op &op, int fi)
{
   const int n_int = op.GetFaultInteriorFaces().Size();
   return (fi < n_int) ? pm.GetInteriorFaceTransformations(op.GetFaultInteriorFaces()[fi])
                       : pm.GetSharedFaceTransformations(op.GetFaultSharedFaces()[fi - n_int]);
}

static bool Same3(const real_t a[3], const real_t b[3])
{
   return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}
static real_t Dist3(const real_t a[3], const real_t b[3])
{
   return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
                    + (a[2] - b[2]) * (a[2] - b[2]));
}
// largest deviation of (n, d, s) from an orthonormal triad with d x s = n
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

static void Random(Vector &v, unsigned seed, real_t scale)
{
   std::mt19937_64 rng(seed + 1000u * static_cast<unsigned>(g_rank));
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   for (int i = 0; i < v.Size(); i++) { v(i) = scale * U(rng); }
}

static real_t RelDiff(const Vector &a, const Vector &b, MPI_Comm c)
{
   real_t num = 0.0, den = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      num = std::max(num, std::abs(a(i) - b(i)));
      den = std::max(den, std::abs(b(i)));
   }
   num = RMax(num, c);
   den = RMax(den, c);
   return (den > 0.0) ? num / den : num;
}
static bool SameVec(const Vector &a, const Vector &b, MPI_Comm c)
{
   bool ok = (a.Size() == b.Size());
   for (int i = 0; ok && i < a.Size(); i++) { ok = (a(i) == b(i)); }
   return AllTrue(ok, c);
}

// every local fault DOF's (x, y) and values, gathered on every rank
static std::vector<std::vector<real_t>> GatherRows(const std::vector<std::vector<real_t>> &rows,
                                                   int width, MPI_Comm c)
{
   std::vector<double> loc;
   for (const auto &r : rows) { loc.insert(loc.end(), r.begin(), r.end()); }
   int nl = static_cast<int>(loc.size()), np = 0;
   MPI_Comm_size(c, &np);
   std::vector<int> cnt(np), dsp(np, 0);
   MPI_Allgather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, c);
   for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
   std::vector<double> all(dsp[np - 1] + cnt[np - 1]);
   MPI_Allgatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(), MPI_DOUBLE, c);
   std::vector<std::vector<real_t>> out;
   for (size_t o = 0; o + width <= all.size(); o += width)
   {
      out.emplace_back(all.begin() + o, all.begin() + o + width);
   }
   return out;
}

// rows grouped by position (x, y) rounded to 1 mm
static std::map<std::pair<long long, long long>, std::vector<std::vector<real_t>>>
ByPosition(const std::vector<std::vector<real_t>> &rows)
{
   std::map<std::pair<long long, long long>, std::vector<std::vector<real_t>>> g;
   for (const auto &r : rows)
   {
      g[{std::llround(r[0] * 1e3), std::llround(r[1] * 1e3)}].push_back(r);
   }
   return g;
}

// ------------------------------------------------------------------ F
static void TestFrames(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== F: the frames of the arc, averaged ===\n"; }
   std::unique_ptr<Mesh> sm(ArcMesh());
   ParMesh pm(comm, *sm);
   auto on = MakeOp(pm, Source(true));
   auto off = MakeOp(pm, Source(false));
   const FaultFrameField &field = *on->FrameField();
   const FaultBasis &fb = *on->GetFaultBasis(), &fb0 = *off->GetFaultBasis();
   const int nf = on->GetNumFaultFaces(), nbf = on->GetNbfPerFace();
   const IntegrationRule &ir = IntRules.Get(Geometry::SEGMENT, 3);
   int varies = 0;
   bool qp_same = true, keep_same = true, dof_same = true;
   real_t orth = 0.0;
   DenseMatrix B;
   on->GetFaultDOFBasis(B);
   const IntegrationRule &nir = on->GetFaceQuadrature()->GetNodalRule();
   std::vector<std::vector<real_t>> rows;   // x, y, oriented normal (3), tip flag
   Vector x;
   for (int fi = 0; fi < nf; fi++)
   {
      const FaultBasisData &b = fb.GetBasis(fi), &b0 = fb0.GetBasis(fi);
      if (field.FrameVaries(fi)) { varies++; }
      if (static_cast<int>(b.dof_triad.size()) != nbf) { dof_same = false; continue; }
      const real_t eps = b.sign_flipped ? -1.0 : 1.0;
      FaceElementTransformations *Tr = FaceTr(pm, *on, fi);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         Tr->Face->Transform(ir.IntPoint(q), x);
         const real_t x3[3] = {x(0), x(1), 0.0};
         ff::Frame fr;
         field.FrameAt(fi, x3, fr);
         const FaultBasisQPData &qd = b.qp_data[q], &q0 = b0.qp_data[q];
         for (int d = 0; d < 3; d++)
         {
            qp_same = qp_same && qd.normal[d] == eps * fr.n[d] && qd.tangent1[d] == eps * fr.d[d]
                      && qd.tangent2[d] == eps * fr.s[d];
         }
         keep_same = keep_same && qd.nl == q0.nl && qd.sign_flipped == q0.sign_flipped;
         orth = std::max(orth, TriadError(fr.n, fr.d, fr.s));
      }
      for (int k = 0; k < nbf; k++)
      {
         Tr->Face->Transform(nir.IntPoint(k), x);
         const real_t x3[3] = {x(0), x(1), 0.0};
         ff::Frame fr;
         field.FrameAt(fi, x3, fr);
         const int col = fi * nbf + k;
         for (int d = 0; d < 3; d++)
         {
            dof_same = dof_same && B(d, col) == eps * fr.n[d] && B(3 + d, col) == eps * fr.d[d]
                       && B(6 + d, col) == eps * fr.s[d]
                       && B(d, col) == b.dof_triad[k].normal[d];
         }
         rows.push_back({x(0), x(1), eps * B(0, col), eps * B(1, col), eps * B(2, col)});
      }
   }
   const bool all_vary = (varies == nf);
   varies = ISum(varies, comm);
   TEST_ASSERT(AllTrue(all_vary, comm) && varies >= 12, "F every fault face varies and carries "
               "per-DOF triads (" << varies << " face copies)");
   TEST_ASSERT(AllTrue(dof_same, comm), "F GetFaultDOFBasis = dof_triad = eps_f FrameAt at the "
               "nodal points, bit for bit");
   TEST_ASSERT(AllTrue(qp_same, comm), "F qp_data = eps_f FrameAt at the quadrature points, bit "
               "for bit");
   TEST_ASSERT(AllTrue(keep_same, comm), "F qp_data nl and sign_flipped as in face mode");
   orth = RMax(orth, comm);
   TEST_ASSERT(orth < 1e-15, "F every quadrature-point triad orthonormal with d x s = n ("
               << orth << ")");
   // vertices: every DOF there (both faces, every copy) has the same
   // oriented normal; interior vertices radial
   int nvert = 0, ntip = 0;
   bool same = true;
   real_t rad_err = 0.0;
   for (const auto &g : ByPosition(GatherRows(rows, 5, comm)))
   {
      const auto &v = g.second;
      for (const auto &r : v)
      {
         same = same && r[2] == v[0][2] && r[3] == v[0][3] && r[4] == v[0][4];
      }
      const real_t px = v[0][0], py = v[0][1];
      const real_t phi = std::atan2(px, kR - py) * 180.0 / M_PI;
      if (std::abs(std::abs(phi) - kFaultPhi) < 1e-6) { ntip++; continue; }
      nvert++;
      const real_t l = std::sqrt(px * px + (py - kR) * (py - kR));
      const real_t rd[3] = {px / l, (py - kR) / l, 0.0};
      const real_t sgn = (ff::Dot(rd, &v[0][2]) < 0.0) ? -1.0 : 1.0;
      const real_t rs[3] = {sgn * rd[0], sgn * rd[1], 0.0};
      rad_err = std::max(rad_err, Dist3(rs, &v[0][2]));
   }
   TEST_ASSERT(nvert == 11 && ntip == 2, "F 11 interior vertices and 2 tips (" << nvert << ", "
               << ntip << "), " << varies << " face copies vary");
   TEST_ASSERT(same, "F at every vertex all DOFs (both faces, every copy) have the same oriented "
               "normal, bit for bit");
   TEST_ASSERT(rad_err < 1e-13, "F the interior vertex normals are radial to 1e-13, the plan's "
               "T2 (" << rad_err << ")");
}

// ------------------------------------------------------------------ T1
static void CompareOps(MPI_Comm comm, ParMesh &pm, Op &a, Op &b, const std::string &label)
{
   const FaultBasis &fa = *a.GetFaultBasis(), &fb = *b.GetFaultBasis();
   const int nf = a.GetNumFaultFaces();
   bool qp = (nf == b.GetNumFaultFaces()), none = true;
   for (int fi = 0; qp && fi < nf; fi++)
   {
      const FaultBasisData &x = fa.GetBasis(fi), &y = fb.GetBasis(fi);
      none = none && x.dof_triad.empty() && y.dof_triad.empty();
      qp = qp && x.qp_data.size() == y.qp_data.size();
      for (size_t q = 0; qp && q < x.qp_data.size(); q++)
      {
         const FaultBasisQPData &p = x.qp_data[q], &r = y.qp_data[q];
         qp = Same3(p.normal, r.normal) && Same3(p.tangent1, r.tangent1)
              && Same3(p.tangent2, r.tangent2) && p.nl == r.nl && p.sign_flipped == r.sign_flipped;
      }
   }
   TEST_ASSERT(AllTrue(none, comm), "T1 " << label << ": no face carries per-DOF triads");
   TEST_ASSERT(AllTrue(qp, comm), "T1 " << label << ": qp_data identical to the option off, bit "
               "for bit");
   DenseMatrix Ba, Bb;
   a.GetFaultDOFBasis(Ba);
   b.GetFaultDOFBasis(Bb);
   Vector va(Ba.Data(), Ba.Height() * Ba.Width()), vb(Bb.Data(), Bb.Height() * Bb.Width());
   TEST_ASSERT(SameVec(va, vb, comm), "T1 " << label << ": GetFaultDOFBasis identical, bit for "
               "bit");
   const int nd = a.GetNumFaultDOFs();
   ParGridFunction u(Space(a));
   Vector slip(2 * nd);
   Random(u, 11, 1e-3);
   Random(slip, 12, 1.0);
   Vector ta, tna, tb, tnb;
   a.ComputeTraction(u, slip, ta, &tna);
   b.ComputeTraction(u, slip, tb, &tnb);
   TEST_ASSERT(SameVec(ta, tb, comm) && SameVec(tna, tnb, comm), "T1 " << label
               << ": ComputeTraction of a random u and slip identical, bit for bit");
   (void)pm;
}

static void TestRegression(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== T1: unchanged where the frame does not vary ===\n"; }
   {
      std::unique_ptr<Mesh> sm(ArcMesh(true));
      ParMesh pm(comm, *sm);
      auto on = MakeOp(pm, Source(true));
      auto off = MakeOp(pm, Source(false));
      CompareOps(comm, pm, *on, *off, "straight fault, option on");
   }
   {
      std::unique_ptr<Mesh> sm(ArcMesh());
      ParMesh pm(comm, *sm);
      auto on = MakeOp(pm, Source(true, 1.0));
      auto off = MakeOp(pm, Source(false));
      CompareOps(comm, pm, *on, *off, "arc, threshold 1 deg (all kinks)");
      bool none = true;
      for (int fi = 0; fi < off->GetNumFaultFaces(); fi++)
      {
         none = none && off->GetFaultBasis()->GetBasis(fi).dof_triad.empty();
      }
      TEST_ASSERT(AllTrue(none, comm), "T1 arc, option off: no per-DOF triads (face mode "
                  "untouched)");
   }
}

// ------------------------------------------------------------------ T4
// S(a, b) = sum over local fault faces of a^T M tau(b), the tangential slots
static real_t WorkPair(const Op &op, const Vector &a, const Vector &tau, MPI_Comm comm)
{
   const int nf = op.GetNumFaultFaces(), nbf = op.GetNbfPerFace();
   const DenseMatrix &e = op.GetFaceQuadrature()->BasisAtQuadPoints();
   const IntegrationRule &ir = IntRules.Get(Geometry::SEGMENT, 3);
   real_t s = 0.0;
   for (int fi = 0; fi < nf; fi++)
   {
      const FaultBasisData &b = op.GetFaultBasis()->GetBasis(fi);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const real_t wn = ir.IntPoint(q).weight * b.qp_data[q].nl;
         for (int c = 0; c < 2; c++)
         {
            real_t aq = 0.0, tq = 0.0;
            for (int k = 0; k < nbf; k++)
            {
               aq += e(k, q) * a(2 * (fi * nbf + k) + c);
               tq += e(k, q) * tau(2 * (fi * nbf + k) + c);
            }
            s += wn * aq * tq;
         }
      }
   }
   return RSum(s, comm);
}

static void TestAdjoint(MPI_Comm comm, bool averaged)
{
   const std::string label = averaged ? "arc averaged" : "arc face mode";
   if (g_rank == 0) { std::cout << "\n=== T4: adjoint, " << label << " ===\n"; }
   std::unique_ptr<Mesh> sm(ArcMesh());
   ParMesh pm(comm, *sm);
   auto op = MakeOp(pm, Source(averaged));
   const int nd = op->GetNumFaultDOFs();
   ParGridFunction zero(Space(*op));
   zero = 0.0;
   Vector a(2 * nd), b(2 * nd), ta, tb, na, nbv;
   Random(a, 21, 1.0);
   Random(b, 22, 1.0);
   op->ComputeTraction(zero, a, ta, &na);
   op->ComputeTraction(zero, b, tb, &nbv);
   const real_t sab = WorkPair(*op, a, tb, comm), sba = WorkPair(*op, b, ta, comm);
   const real_t rel = std::abs(sab - sba) / std::max(std::abs(sab), std::abs(sba));
   TEST_ASSERT(rel < 1e-12, "T4 " << label << ": a^T M tau(b) = b^T M tau(a) (" << sab << ", "
               << sba << ", relative " << rel << ")");
   real_t nmax = std::max(na.Normlinf(), nbv.Normlinf()), tmax = std::max(ta.Normlinf(),
                                                                          tb.Normlinf());
   nmax = RMax(nmax, comm);
   tmax = RMax(tmax, comm);
   TEST_ASSERT(nmax <= 1e-12 * tmax, "T4 " << label << ": no normal traction from the slip at "
               "u = 0 (" << nmax << " vs shear " << tmax << ")");
   // the per-face operators, built after the frames
   op->BuildFaultFaceOps();
   ParGridFunction u(Space(*op));
   Random(u, 23, 1e-3);
   Vector tr, tn, trf, tnf;
   op->ComputeTraction(u, a, tr, &tn);
   op->ComputeTractionFast(u, a, trf, tnf);
   const real_t e1 = RelDiff(trf, tr, comm), e2 = RelDiff(tnf, tn, comm);
   TEST_ASSERT(e1 < 1e-11 && e2 < 1e-11, "T4 " << label << ": ComputeTractionFast = "
               "ComputeTraction (shear " << e1 << ", normal " << e2 << ")");
   Vector g1(Space(*op)->GetTrueVSize()), g2(g1.Size());
   op->AssembleSlipOnlyRHS(g1, a);
   g2 = 0.0;
   op->AddSlipLoadFast(a, g2);
   const real_t e3 = RelDiff(g2, g1, comm);
   TEST_ASSERT(e3 < 1e-11, "T4 " << label << ": AddSlipLoadFast = AssembleSlipOnlyRHS ("
               << e3 << ")");
}

// ------------------------------------------------------------------ S
// the spread of sigma_n among the DOFs of each interior vertex [Pa] (largest
// minus smallest, both faces, every copy), keyed by the vertex's arc length s
// rounded to 1 m, after one static solve with the tapered slip
static std::map<long long, real_t> StaticJumps(MPI_Comm comm, bool straight, bool averaged)
{
   std::unique_ptr<Mesh> sm(ArcMesh(straight));
   ParMesh pm(comm, *sm);
   auto op = MakeOp(pm, Source(averaged));
   const int nd = op->GetNumFaultDOFs();
   const int slot = op->FrameField()->InPlaneSlot();
   Vector xyz;
   op->GetFaultDOFCoords3D(xyz);
   const real_t half = kR * kFaultPhi * M_PI / 180.0;
   auto arc_s = [&](real_t px, real_t py) { return straight ? px : kR * std::atan2(px, kR - py); };
   Vector slip(2 * nd);
   slip = 0.0;
   for (int k = 0; k < nd; k++)
   {
      const real_t c = std::cos(0.5 * M_PI * arc_s(xyz(3 * k), xyz(3 * k + 1)) / half);
      slip(2 * k + slot) = c * c;
   }
   ParGridFunction u(Space(*op));
   op->Solve(0.0, slip, u);
   Vector tr, tn;
   op->ComputeTraction(u, slip, tr, &tn);
   std::vector<std::vector<real_t>> rows;
   for (int k = 0; k < nd; k++) { rows.push_back({xyz(3 * k), xyz(3 * k + 1), tn(k)}); }
   std::map<long long, real_t> jumps;
   for (const auto &g : ByPosition(GatherRows(rows, 3, comm)))
   {
      const auto &v = g.second;
      const real_t s = arc_s(v[0][0], v[0][1]);
      if (std::abs(std::abs(s) - half) < 1.0) { continue; }   // a tip
      real_t lo = v[0][2], hi = v[0][2];
      for (const auto &r : v) { lo = std::min(lo, r[2]); hi = std::max(hi, r[2]); }
      jumps[std::llround(s)] = hi - lo;
   }
   return jumps;
}

// The straight fault is the same grid, so its jump at the same s is the
// discretization's own (slip gradient, p = 1, h about 350 m); the excess
// over it is the part the vertex turn adds.
static void TestSawtooth(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== S: static sawtooth, tapered slip (1 m peak) ===\n"; }
   const auto js = StaticJumps(comm, true, false);
   const auto jf = StaticJumps(comm, false, false);
   const auto ja = StaticJumps(comm, false, true);
   bool keys = (js.size() == 11 && jf.size() == 11 && ja.size() == 11);
   real_t ms = 0.0, mf = 0.0, ma = 0.0, xf = 0.0, xa = 0.0;
   if (g_rank == 0) { std::printf("  s [km]   straight   arc face   arc averaged   (MPa)\n"); }
   for (const auto &kv : js)
   {
      if (!jf.count(kv.first) || !ja.count(kv.first)) { keys = false; continue; }
      const real_t a = kv.second, f = jf.at(kv.first), v = ja.at(kv.first);
      ms = std::max(ms, a);
      mf = std::max(mf, f);
      ma = std::max(ma, v);
      xf = std::max(xf, std::abs(f - a));
      xa = std::max(xa, std::abs(v - a));
      if (g_rank == 0)
      {
         std::printf("  %6.3f   %8.4f   %8.4f   %8.4f\n", kv.first / 1e3, a / 1e6, f / 1e6,
                     v / 1e6);
      }
   }
   if (g_rank == 0)
   {
      std::printf("  largest jump: straight %.4f, arc face %.4f, arc averaged %.4f MPa; largest "
                  "excess over the straight fault: face %.4f, averaged %.4f MPa\n", ms / 1e6,
                  mf / 1e6, ma / 1e6, xf / 1e6, xa / 1e6);
   }
   TEST_ASSERT(keys, "S the same 11 interior vertices (by arc length) on the three meshes");
   TEST_ASSERT(xf > ms, "S the arc in face mode adds a sawtooth larger than the straight fault's "
               "own jump (excess " << xf / 1e6 << " MPa, straight " << ms / 1e6 << " MPa)");
   TEST_ASSERT(xa < 0.1 * xf, "S averaged vertex normals remove at least nine tenths of it "
               "(excess " << xa / 1e6 << " MPa vs " << xf / 1e6 << " MPa in face mode)");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "test_vertex_normal_kernels (vertex-normal plan rev1, phase P2) on "
                << Mpi::WorldSize() << " rank(s)\n";
   }
   TestFrames(comm);
   TestRegression(comm);
   TestAdjoint(comm, true);
   TestAdjoint(comm, false);
   TestSawtooth(comm);
   if (g_rank == 0)
   {
      std::cout << "\n" << num_tests - num_failed << " / " << num_tests << " tests passed\n";
   }
   return num_failed == 0 ? 0 : 1;
}
