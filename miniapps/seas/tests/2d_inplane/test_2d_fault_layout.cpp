// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_2d_fault_layout.cpp — 2D in-plane design rev 3,
// phase B: the dimension switch on 2D triangle meshes built here (part 1)
// and on the committed 2D meshes (part 2).
//
//   L1  LoadScaledMesh reads a 2D triangle mesh: Dimension 2, the element and
//       boundary counts and tags of the file, every vertex times the scale.
//   C1  census on the BP3 layout in miniature (a dipping fault from the free
//       surface to W_f, its creep continuation to the bottom), dips 30, 60,
//       90: one patch, 0 junction vertices, 2 tips; the creep segment one
//       patch; 1 seam vertex; n = (sin, cos, 0) with ref_normal (1, 0); the
//       S3 plus element on the side n points away from; the report the same
//       at every rank count (P1).
//   C2  a Y-shaped fault: 1 junction vertex, 3 tips, 3 patches.
//   C3  the 2D slice of the phase A box: the fault census, the oriented
//       normals, the frames and the plus side equal the 3D field's.
//   S1  the slot rule (derivation 10.2): up = [0, 1] (vertical section) gives
//       s = -+e_z and d = the down-dip tangent; up = [0, 0, 1] (horizontal
//       slice) gives d = -+e_z and s in the plane; the out-of-plane slot is an
//       exact zero in the field's frames, in GetFaultDOFBasis and in
//       FaultBasis::ProjectTraction / EmbedSlip.
//   S2  the F6 fallback in 2D (a fault segment parallel to an in-plane up):
//       s along the segment, d = +-e_z, the out-of-plane slot still exact.
//   D1  the fault DOF layout at p = 1, 2, 3: p + 1 nodes per face, every fault
//       face owned once, the owned -> local scatter puts f(x) on every local
//       DOF (the ghost copies of shared faces included), and the table of
//       (position, n, d, s) equals the one of a 1-rank run.
//   M1  (part 2) the committed meshes bp3/mesh/bp3_d{30,60,90}_25m.msh,
//       tpv205/mesh/tpv205_2d_100m.msh, tpv105/mesh/tpv105_2d_100m.msh:
//       Dimension 2, the triangle count and the tag census of the file
//       (segments, exterior, interior per tag), the fault's length and
//       bounding box; the census: one patch, 0 junction vertices, 2 tips
//       (BP3: the creep segment one patch, 1 seam vertex), the oriented
//       normal, the plus side and the slot of the deck's up (BP3 [0, 1]: dip;
//       TPV [0, 0, 1]: strike), the report the same as on 1 rank.
//   D2  (part 2) D1 on the BP3 mesh (dip 60, 25 m), p = 1, 2, 3, with
//       METIS's partition and with one that cuts along the fault (every
//       fault segment a rank-boundary face).
//   Q1  (phase F) the QD stage operators and the dense Green's function in
//       2D: on the BP3 layout in miniature (dip 60, p = 1, 2) with a plate
//       loading linear in t, the traction and normal traction of a test
//       slip at t = 100 yr from the face loops (Solve, ComputeTraction), the
//       stage operators (SolveStage) and G (GreenStage) agree to 1e-10 of
//       their size, with and without a strike slot in the slip (which
//       embeds to no displacement jump in 2D).
//
//   seas_test_2d_fault_layout             (1, 2 or 4 ranks)
//   seas_test_2d_fault_layout --refuse N  one refusal, must exit non-zero:
//     1 a quadrilateral mesh in LoadScaledMesh, 2 a 2D fault without the
//     frame field, 3 up off the plane, 4 ref_normal off the plane, 5 a G
//     cache file in 2D (the cache and the H-matrix key fault DOFs by 3D
//     position; the stage operators and the dense G open in phase F), 6 up
//     parallel to n in the QD basis, 7 BR2 in 2D.
//
// Scope (Chunhui Zhao, 2026-10-25).

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "frame_box.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
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

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static bool AllTrue(bool b, MPI_Comm c) { return ISum(b ? 0 : 1, c) == 0; }

// ------------------------------------------------------------------ meshes
// The BP3 layout in miniature: a grid (s, y), s in [-Ls, Ls], y in [-D, 0],
// mapped by x = s - y cot(dip), so the grid line s = 0 is the fault
// x = -y cot(dip) (y up, the free surface y = 0).  On s = 0: attr 3 for
// depth -y < Wf, attr 5 below (the creep segment); top attr 1, the other
// sides attr 5.  Triangles: each cell split along its (s, y) diagonal.
static const real_t kLs = 20e3, kD = 20e3, kWf = 12e3;
static const int kNs = 8, kNy = 10;

static Mesh *ShearBox(real_t dip_deg, int quad = 0)
{
   const real_t c = (dip_deg == 90.0) ? 0.0 : 1.0 / std::tan(dip_deg * M_PI / 180.0);
   const int nvx = kNs + 1, nvy = kNy + 1;
   const int ne = quad ? kNs * kNy : 2 * kNs * kNy;
   Mesh *m = new Mesh(2, nvx * nvy, ne, 0, 2);
   for (int j = 0; j < nvy; j++)
      for (int i = 0; i < nvx; i++)
      {
         const real_t s = -kLs + 2.0 * kLs * i / kNs, y = -kD + kD * j / kNy;
         m->AddVertex(s - y * c, y);
      }
   auto v = [&](int i, int j) { return i + nvx * j; };
   for (int j = 0; j < kNy; j++)
      for (int i = 0; i < kNs; i++)
      {
         if (quad) { m->AddQuad(v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1), 1); continue; }
         m->AddTriangle(v(i, j), v(i + 1, j), v(i + 1, j + 1), 1);
         m->AddTriangle(v(i, j), v(i + 1, j + 1), v(i, j + 1), 1);
      }
   for (int i = 0; i < kNs; i++)
   {
      m->AddBdrSegment(v(i, 0), v(i + 1, 0), 5);
      m->AddBdrSegment(v(i, kNy), v(i + 1, kNy), 1);
   }
   for (int j = 0; j < kNy; j++)
   {
      m->AddBdrSegment(v(0, j), v(0, j + 1), 5);
      m->AddBdrSegment(v(kNs, j), v(kNs, j + 1), 5);
      const real_t ymid = -kD + kD * (j + 0.5) / kNy;
      m->AddBdrSegment(v(kNs / 2, j), v(kNs / 2, j + 1), (-ymid < kWf) ? 3 : 5);
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

// A Y-shaped fault in the square [-L, L]^2 (n x n cells, split along the
// diagonal (i, j)-(i+1, j+1)): the line y = 0 for |x| < L/2 and the branch
// y = x for 0 < x < L/2, all attr 3; outer sides attr 5.  The branch point
// (0, 0) carries three fault edges.
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

// The 2D slice of frame_box kind 0 (the coplanar box, unrotated): the same
// 4 x 4 grid in (x, y), the fault on y = 0 for |x| < L/2 (attr 3), the band
// on the rest of y = 0 (attr 6), the side y = L attr 1, the others attr 5.
static Mesh *BoxSlice()
{
   const int n = 4;
   const real_t L = frame_box::kL;
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
      m->AddBdrSegment(v(k, n), v(k + 1, n), 1);
      m->AddBdrSegment(v(0, k), v(0, k + 1), 5);
      m->AddBdrSegment(v(n, k), v(n, k + 1), 5);
      const real_t xm = -L + 2.0 * L * (k + 0.5) / n;
      m->AddBdrSegment(v(k, n / 2), v(k + 1, n / 2), (std::abs(xm) < 0.5 * L) ? 3 : 6);
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

// ---------------------------------------------------------------- operator
static std::shared_ptr<const FaultFrameConfig> FrameConfig(const real_t ref[3], const real_t up[3],
                                                           int fault_attr = 3)
{
   auto c = std::make_shared<FaultFrameConfig>();
   for (int d = 0; d < 3; d++)
   {
      c->orientation.ref_normal[d] = ref[d];
      c->orientation.up[d] = up[d];
      c->up[d] = up[d];
   }
   c->orientation.fault_attrs = {fault_attr};
   return c;
}

static BoundaryConfig MakeBC(int dim, std::set<int> dirichlet, int fault_attr = 3,
                             std::set<int> natural = {1}, std::set<int> absorbing = {})
{
   BoundaryConfig bc;
   bc.fault_attr = fault_attr;
   bc.natural_attrs = natural;
   bc.absorbing_attrs = absorbing;
   bc.dirichlet_attrs = dirichlet;
   bc.default_dirichlet_func = [dim](const Vector &, real_t, Vector &u) { u.SetSize(dim); u = 0.0; };
   return bc;
}

using Op = ElasticityDomainOperator<ParMesh>;

static std::unique_ptr<Op> MakeOp(ParMesh &pm, int order, FaultFrameSource src,
                                  std::set<int> dirichlet = {5},
                                  DGMethod method = DGMethod::IP, int fault_attr = 3,
                                  std::set<int> natural = {1}, std::set<int> absorbing = {})
{
   LinearElastic *model = new LinearElastic(kLambda, kMu);   // leaked: test lifetime
   DomainConfig dc;
   dc.fault_frames = src;
   return std::make_unique<Op>(pm, order, *model, 0.0, kWf, 2.0 * kLs,
                               MakeBC(pm.Dimension(), dirichlet, fault_attr, natural, absorbing), method,
                               SolverType::MUMPS, dc);
}

static FaultFrameSource Source(const real_t ref[3], const real_t up[3], int fault_attr = 3)
{
   FaultFrameSource s;
   s.config = FrameConfig(ref, up, fault_attr);
   return s;
}

static real_t Dist3(const real_t a[3], const real_t b[3])
{
   return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
                    + (a[2] - b[2]) * (a[2] - b[2]));
}

// centroid of local element e (padded to 3)
static void Centroid(ParMesh &pm, int e, real_t c[3])
{
   Array<int> ev;
   pm.GetElementVertices(e, ev);
   c[0] = c[1] = c[2] = 0.0;
   for (int a = 0; a < ev.Size(); a++)
      for (int d = 0; d < pm.SpaceDimension(); d++) { c[d] += pm.GetVertex(ev[a])[d] / ev.Size(); }
}

// face centre of local fault face fi of the operator (padded to 3)
static void FaceCentre(ParMesh &pm, Op &op, int fi, real_t c[3])
{
   const Array<int> &fint = op.GetFaultInteriorFaces();
   const Array<int> &fsh = op.GetFaultSharedFaces();
   const int lf = (fi < fint.Size()) ? fint[fi] : pm.GetSharedFace(fsh[fi - fint.Size()]);
   Array<int> fv;
   pm.GetFaceVertices(lf, fv);
   c[0] = c[1] = c[2] = 0.0;
   for (int a = 0; a < fv.Size(); a++)
      for (int d = 0; d < pm.SpaceDimension(); d++) { c[d] += pm.GetVertex(fv[a])[d] / fv.Size(); }
}

// --------------------------------------------------------------------- L1
static void TestLoader(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== L1 LoadScaledMesh on a 2D triangle mesh ===\n"; }
   const std::string path = "test_2d_fault_layout_L1.mesh";
   std::unique_ptr<Mesh> src(ShearBox(60.0));
   if (g_rank == 0)
   {
      std::ofstream ofs(path);
      ofs.precision(17);
      src->Print(ofs);
   }
   MPI_Barrier(comm);
   auto m1 = cycle::LoadScaledMesh(path, 1.0);
   auto mk = cycle::LoadScaledMesh(path, 1000.0);
   TEST_ASSERT(mk->Dimension() == 2 && mk->SpaceDimension() == 2,
               "L1 Dimension 2, SpaceDimension 2");
   TEST_ASSERT(mk->GetNE() == src->GetNE() && mk->GetNBE() == src->GetNBE(),
               "L1 " << mk->GetNE() << " triangles, " << mk->GetNBE() << " boundary segments");
   bool scaled = true;
   for (int i = 0; i < mk->GetNV(); i++)
      for (int d = 0; d < 2; d++)
      {
         if (mk->GetVertex(i)[d] != 1000.0 * m1->GetVertex(i)[d]) { scaled = false; }
      }
   TEST_ASSERT(scaled, "L1 every vertex = 1000 x the file's (bit for bit)");
   const auto census = cycle::CountBoundary(*mk);
   const int n1 = census.count(1) ? census.at(1).faces : 0;
   const int n3 = census.count(3) ? census.at(3).faces : 0;
   const int n5 = census.count(5) ? census.at(5).faces : 0;
   const int i3 = census.count(3) ? census.at(3).interior : 0;
   TEST_ASSERT(n1 == kNs && n3 == 6 && n5 == kNs + 2 * kNy + 4 && i3 == 6
               && census.at(3).simplex_faces == 6,
               "L1 tag census: 1 (free surface) " << n1 << ", 3 (fault, interior) " << n3
               << ", 5 (walls + creep) " << n5);
   TEST_ASSERT(cycle::FaultDOFsPerFace(2, 2) == 3 && cycle::FaultDOFsPerFace(2, 3) == 6,
               "L1 FaultDOFsPerFace: p + 1 in 2D, (p + 1)(p + 2) / 2 in 3D");
   MPI_Barrier(comm);
   if (g_rank == 0) { std::remove(path.c_str()); }
}

// --------------------------------------------------------- C1, S1 (BP3 layout)
struct ReportKey
{
   int faces, comps, junc, tips, cfaces, ccomps, seams;
   real_t pin, side;
};

static ReportKey KeyOf(const FaultFrameReport &R)
{
   return {R.fault.faces, R.fault.components, R.fault.junction_edges, R.fault.tip_edges,
           R.creep.faces, R.creep.components, R.seam_edges, R.fault.min_pin_margin,
           R.fault.min_side_margin};
}

static void TestBP3Census(Mesh &serial, real_t dip_deg, MPI_Comm comm, const ReportKey *ref1)
{
   ParMesh pm(comm, serial);
   const real_t ref[3] = {1, 0, 0}, upB[3] = {0, 1, 0};
   auto op = MakeOp(pm, 1, Source(ref, upB));
   auto F = op->FrameField();
   TEST_ASSERT(F != nullptr && F->Dimension() == 2, "C1 dip " << dip_deg
               << ": the operator built a 2D FaultFrameField");
   if (!F) { return; }
   const FaultFrameReport &R = F->Report();
   TEST_ASSERT(R.fault.faces == 6 && R.fault.components == 1 && R.fault.junction_edges == 0
               && R.fault.tip_edges == 2,
               "C1 dip " << dip_deg << ": fault 6 segments, 1 patch, 0 junction vertices, 2 tips ("
               << R.fault.faces << ", " << R.fault.components << ", " << R.fault.junction_edges
               << ", " << R.fault.tip_edges << ")");
   TEST_ASSERT(R.creep.faces == 4 && R.creep.components == 1 && R.seam_edges == 1,
               "C1 dip " << dip_deg << ": creep 4 segments, 1 patch, 1 seam vertex ("
               << R.creep.faces << ", " << R.creep.components << ", " << R.seam_edges << ")");
   TEST_ASSERT(R.fault.parity_global_ok && R.fault.parity_per_patch_ok, "C1 dip " << dip_deg
               << ": J3 parity per patch and global");
   const real_t ps = std::sin(dip_deg * M_PI / 180.0), pc = std::cos(dip_deg * M_PI / 180.0);
   const real_t n_exp[3] = {ps, pc, 0.0}, d_exp[3] = {pc, -ps, 0.0};
   real_t en = 0, ed = 0, ec = 0;
   bool side_ok = true, slots = true;
   for (int fi = 0; fi < F->NumFaultFaces(); fi++)
   {
      en = std::max(en, Dist3(F->Normal(fi), n_exp));
      const ff::Frame &fr = F->Frame(fi);
      ed = std::max(ed, Dist3(fr.d, d_exp));
      // case B: s = -e_z (in-plane entries exactly 0), d in the plane
      if (!(fr.s[0] == 0.0 && fr.s[1] == 0.0 && fr.s[2] < 0.0 && fr.d[2] == 0.0
            && fr.n[2] == 0.0 && !fr.fallback)) { slots = false; }
      // S3: the plus element lies on the side n points away from
      const int pe = F->FaultCensus().PlusElement(fi);
      if (pe >= 0)
      {
         real_t ce[3], cf[3];
         Centroid(pm, pe, ce);
         FaceCentre(pm, *op, fi, cf);
         real_t m = 0.0;
         for (int d = 0; d < 3; d++) { m += (ce[d] - cf[d]) * F->Normal(fi)[d]; }
         if (!(m < 0.0)) { side_ok = false; }
      }
   }
   const FaultOrientation *C = F->CreepCensus();
   for (int ci = 0; C && ci < C->NumLocalFaces(); ci++) { ec = std::max(ec, Dist3(C->Normal(ci), n_exp)); }
   en = RMax(en, comm); ed = RMax(ed, comm); ec = RMax(ec, comm);
   TEST_ASSERT(en <= 1e-15 && ec <= 1e-15, "C1 dip " << dip_deg
               << ": every oriented normal = (sin, cos, 0) (fault " << en << ", creep " << ec << ")");
   TEST_ASSERT(AllTrue(side_ok, comm), "C1 dip " << dip_deg
               << ": the plus element is on the side n points away from (the footwall)");
   TEST_ASSERT(AllTrue(slots, comm) && ed <= 1e-15, "S1 dip " << dip_deg
               << ": up = [0, 1]: s = -e_z exactly in the plane, d = down-dip (" << ed << ")");
   if (ref1)
   {
      const ReportKey k = KeyOf(R);
      TEST_ASSERT(k.faces == ref1->faces && k.comps == ref1->comps && k.junc == ref1->junc
                  && k.tips == ref1->tips && k.cfaces == ref1->cfaces && k.ccomps == ref1->ccomps
                  && k.seams == ref1->seams && k.pin == ref1->pin
                  && std::abs(k.side - ref1->side) <= 1e-12 * ref1->side,
                  "C1 dip " << dip_deg << ": the report equals the 1-rank report (P1; pin margin "
                  << k.pin << " bit for bit, side margin " << k.side << ")");
   }
}

// one 1-rank report per dip, built on rank 0 alone (MPI_COMM_SELF)
static ReportKey Report1(Mesh &serial)
{
   ParMesh pm(MPI_COMM_SELF, serial);
   const real_t ref[3] = {1, 0, 0}, upB[3] = {0, 1, 0};
   auto op = MakeOp(pm, 1, Source(ref, upB));
   return KeyOf(op->FrameField()->Report());
}

// -------------------------------------------------------------- S1 slots
static void TestSlots(Mesh &serial, MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== S1 the out-of-plane slot is an exact zero ===\n"; }
   ParMesh pm(comm, serial);
   const real_t ref[3] = {1, 0, 0};
   for (int c = 0; c < 2; c++)
   {
      const real_t upB[3] = {0, 1, 0}, upA[3] = {0, 0, 1};
      auto op = MakeOp(pm, 2, Source(ref, c == 0 ? upB : upA));
      const char *lab = (c == 0) ? "up = [0, 1] (slip in the dip slot)"
                                 : "up = [0, 0, 1] (slip in the strike slot)";
      DenseMatrix B;
      op->GetFaultDOFBasis(B);
      bool zero = true, unit = true;
      for (int k = 0; k < B.Width(); k++)
      {
         if (B(2, k) != 0.0) { zero = false; }
         if (c == 0 && !(B(6, k) == 0.0 && B(7, k) == 0.0 && B(5, k) == 0.0)) { zero = false; }
         if (c == 1 && !(B(3, k) == 0.0 && B(4, k) == 0.0 && B(8, k) == 0.0)) { zero = false; }
         const real_t z = (c == 0) ? B(8, k) : B(5, k);
         if (std::abs(std::abs(z) - 1.0) > 1e-15) { unit = false; }
      }
      const int ndof = ISum(B.Width(), comm);   // collectives outside TEST_ASSERT's message
      TEST_ASSERT(AllTrue(zero && unit, comm), "S1 " << lab << ": GetFaultDOFBasis has n_z = 0 and "
                  "the +-e_z tangent's in-plane entries (and the other tangent's z) exactly 0 at "
                  << ndof << " DOFs");
      const FaultBasis *fb = op->GetFaultBasis();
      bool proj = true, embed = true;
      for (int fi = 0; fi < fb->NumFaces(); fi++)
      {
         const real_t T[2] = {1.234e6 * (fi + 1), -7.89e5};
         real_t tau[2];
         fb->ProjectTraction(fi, T, tau);
         if (tau[(c == 0) ? 1 : 0] != 0.0) { proj = false; }
         const real_t s1[2] = {0.37, 0.0}, s2[2] = {0.37, 5.0};
         const real_t s3[2] = {0.0, 0.37}, s4[2] = {5.0, 0.37};
         real_t a[2], b[2];
         fb->EmbedSlip(fi, c == 0 ? s1 : s3, a);
         fb->EmbedSlip(fi, c == 0 ? s2 : s4, b);
         if (a[0] != b[0] || a[1] != b[1]) { embed = false; }
      }
      TEST_ASSERT(AllTrue(proj, comm), "S1 " << lab << ": ProjectTraction of (Tx, Ty) gives the "
                  "out-of-plane slot exactly 0");
      TEST_ASSERT(AllTrue(embed, comm), "S1 " << lab << ": EmbedSlip ignores the out-of-plane "
                  "slot bit for bit");
      // the kernel triad = sigma_f (n, d, s) at every DOF
      auto F = op->FrameField();
      const int nbf = op->GetNbfPerFace();
      real_t et = 0.0;
      for (int fi = 0; fi < F->NumFaultFaces(); fi++)
      {
         const ff::Frame &fr = F->Frame(fi);
         const real_t sg = F->Sigma(fi);
         for (int k = 0; k < nbf; k++)
            for (int d = 0; d < 3; d++)
            {
               const int col = fi * nbf + k;
               et = std::max({et, std::abs(B(d, col) - sg * fr.n[d]),
                              std::abs(B(3 + d, col) - sg * fr.d[d]),
                              std::abs(B(6 + d, col) - sg * fr.s[d])});
            }
      }
      et = RMax(et, comm);
      TEST_ASSERT(et <= 1e-15, "S1 " << lab << ": kernel triad = sigma_f (n, d, s) (" << et << ")");
   }
}

// ------------------------------------------------------------------- C2
static void TestY(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== C2 a Y-shaped fault ===\n"; }
   std::unique_ptr<Mesh> serial(YMesh());
   ParMesh pm(comm, *serial);
   const real_t ref[3] = {0, 1, 0}, upA[3] = {0, 0, 1};
   auto op = MakeOp(pm, 1, Source(ref, upA));
   const FaultFrameReport &R = op->FrameField()->Report();
   TEST_ASSERT(R.fault.faces == 6 && R.fault.junction_edges == 1 && R.fault.tip_edges == 3
               && R.fault.components == 3,
               "C2 6 segments, 1 junction vertex, 3 tips, 3 patches (" << R.fault.faces << ", "
               << R.fault.junction_edges << ", " << R.fault.tip_edges << ", "
               << R.fault.components << ")");
   TEST_ASSERT(R.creep.faces == 0, "C2 no prescribed-slip faces");
}

// ------------------------------------------------------------------- S2
// The F6 fallback in 2D: the Y fault with up = [0, 1] (in the plane), so the
// main fault's normal is parallel to up.  The frame field (the dynamic
// regime's) takes F6 there: s = the segment's direction, d = s x n = +-e_z.
// (The QD basis refuses the same deck: refusal 6.)  Built directly on a
// 1-rank mesh on rank 0, since the operator would refuse first.
static void TestFallback2D()
{
   if (g_rank != 0) { return; }
   std::cout << "\n=== S2 the F6 fallback in 2D ===\n";
   std::unique_ptr<Mesh> serial(YMesh());
   ParMesh pm(MPI_COMM_SELF, *serial);
   Array<int> fint, none;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) != 3) { continue; }
      const int f = pm.GetBdrElementFaceIndex(be);
      if (pm.FaceIsInterior(f)) { fint.Append(f); }
   }
   const real_t ref[3] = {0, 1, 0}, upB[3] = {0, 1, 0};
   FaultFrameField F(pm, fint, none, none, none, {}, *FrameConfig(ref, upB));
   TEST_ASSERT(F.Report().frame_fallbacks == 4, "S2 the 4 segments on y = 0 (n parallel to up) "
               "take F6 (" << F.Report().frame_fallbacks << ")");
   bool slots = true, f1 = true;
   for (int fi = 0; fi < F.NumFaultFaces(); fi++)
   {
      const ff::Frame &fr = F.Frame(fi);
      const bool on_y0 = (fr.n[0] == 0.0);
      if (on_y0 != fr.fallback) { f1 = false; }
      // one tangent exactly in the plane, the other exactly +-e_z
      const bool d_ez = (fr.d[0] == 0.0 && fr.d[1] == 0.0 && std::abs(std::abs(fr.d[2]) - 1.0) <= 1e-15);
      const bool s_ez = (fr.s[0] == 0.0 && fr.s[1] == 0.0 && std::abs(std::abs(fr.s[2]) - 1.0) <= 1e-15);
      if (fr.fallback && !(d_ez && fr.s[2] == 0.0)) { slots = false; }
      if (!fr.fallback && !(s_ez && fr.d[2] == 0.0)) { slots = false; }
   }
   TEST_ASSERT(f1, "S2 F6 on exactly the y = 0 segments; the branch y = x keeps F1");
   TEST_ASSERT(slots, "S2 the out-of-plane slot is exact: F6 gives d = +-e_z and s in the plane, "
               "F1 (up in the plane) s = +-e_z and d in the plane");
}

// ------------------------------------------------------------------- C3
static void TestSlice(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== C3 the 2D slice of the phase A box ===\n"; }
   const real_t ref[3] = {0, -1, 0}, upA[3] = {0, 0, 1};
   frame_box::Rot I;
   std::unique_ptr<Mesh> s3(frame_box::FrameBox(0, I));
   ParMesh pm3(comm, *s3);
   auto op3 = MakeOp(pm3, 1, Source(ref, upA), {5, 6});
   std::unique_ptr<Mesh> s2(BoxSlice());
   ParMesh pm2(comm, *s2);
   auto op2 = MakeOp(pm2, 1, Source(ref, upA), {5, 6});
   auto F3 = op3->FrameField(), F2 = op2->FrameField();
   const FaultFrameReport &R3 = F3->Report(), &R2 = F2->Report();
   TEST_ASSERT(R2.fault.components == R3.fault.components && R2.fault.junction_edges == 0
               && R3.fault.junction_edges == 0,
               "C3 fault: " << R2.fault.components << " patch in 2D, " << R3.fault.components
               << " in 3D, no junctions");
   TEST_ASSERT(R2.creep.components == 2 && R3.creep.components == 1 && R2.seam_edges == 2
               && R3.seam_edges == 4,
               "C3 band: 2 pieces in the slice (x < -L/2, x > L/2) and 2 seam vertices, 1 patch "
               "and 4 seam edges in 3D (" << R2.creep.components << ", " << R2.seam_edges << "; "
               << R3.creep.components << ", " << R3.seam_edges << ")");
   // every fault face of both: the same n, d, s (bit for bit) and the plus
   // element on y > 0, the side n = (0, -1, 0) points away from
   auto check = [&](ParMesh &pm, Op &op, const FaultFrameField &F, real_t fr_ref[9],
                    bool &same, bool &side) {
      for (int fi = 0; fi < F.NumFaultFaces(); fi++)
      {
         const ff::Frame &fr = F.Frame(fi);
         for (int d = 0; d < 3; d++)
         {
            if (fr.n[d] != fr_ref[d] || fr.d[d] != fr_ref[3 + d] || fr.s[d] != fr_ref[6 + d])
            { same = false; }
         }
         const int pe = F.FaultCensus().PlusElement(fi);
         if (pe >= 0)
         {
            real_t ce[3];
            Centroid(pm, pe, ce);
            if (!(ce[1] > 0.0)) { side = false; }
         }
         (void)op;
      }
   };
   real_t fr_ref[9] = {0, -1, 0, 0, 0, -1, 1, 0, 0};   // n, d = -e_z, s = e_x
   bool same3 = true, side3 = true, same2 = true, side2 = true;
   check(pm3, *op3, *F3, fr_ref, same3, side3);
   check(pm2, *op2, *F2, fr_ref, same2, side2);
   TEST_ASSERT(AllTrue(same3 && same2, comm), "C3 n = (0, -1, 0), d = (0, 0, -1), s = (1, 0, 0) "
               "on every fault face, 3D and 2D, bit for bit");
   TEST_ASSERT(AllTrue(side3 && side2, comm), "C3 the plus element on y > 0 in 3D and in 2D");
}

// ------------------------------------------------------------------- D1
struct Row { real_t x[2]; real_t f[9]; };

static std::vector<Row> DofTable(Op &op)
{
   // owned DOFs: position and the oriented (n, d, s) (sigma_f taken out)
   Vector X, Xo;
   op.GetFaultDOFCoords3D(X);
   DenseMatrix B;
   op.GetFaultDOFBasis(B);
   Vector Bv(9 * B.Width());
   for (int k = 0; k < B.Width(); k++)
      for (int r = 0; r < 9; r++) { Bv(9 * k + r) = B(r, k); }
   Vector Bo;
   op.RestrictToOwnedFault(X, Xo, 3);
   op.RestrictToOwnedFault(Bv, Bo, 9);
   auto F = op.FrameField();
   std::vector<Row> rows(op.GetNumOwnedFaultDOFs());
   for (int k = 0; k < op.GetNumOwnedFaultDOFs(); k++)
   {
      Row &r = rows[k];
      r.x[0] = Xo(3 * k);
      r.x[1] = Xo(3 * k + 1);
      // sign of the stored normal against the field's: the face is the same
      // for every DOF of an owned face, so compare with any fault normal of
      // the patch (one straight fault here)
      const real_t *n = F->Normal(0);
      const real_t s = (Bo(9 * k) * n[0] + Bo(9 * k + 1) * n[1] >= 0.0) ? 1.0 : -1.0;
      for (int c = 0; c < 9; c++) { r.f[c] = s * Bo(9 * k + c); }
   }
   return rows;
}

static std::vector<Row> Gather(const std::vector<Row> &loc, MPI_Comm comm)
{
   const int R = 11;
   std::vector<double> b;
   for (const Row &r : loc)
   {
      b.push_back(r.x[0]); b.push_back(r.x[1]);
      for (int c = 0; c < 9; c++) { b.push_back(r.f[c]); }
   }
   int nl = static_cast<int>(b.size()), np = 0;
   MPI_Comm_size(comm, &np);
   std::vector<int> cnt(np), dsp(np, 0);
   MPI_Gather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm);
   if (g_rank == 0) { for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; } }
   std::vector<double> all(g_rank == 0 ? dsp[np - 1] + cnt[np - 1] : 0);
   MPI_Gatherv(b.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(), MPI_DOUBLE, 0, comm);
   std::vector<Row> out;
   for (size_t o = 0; o + R <= all.size(); o += R)
   {
      Row r;
      r.x[0] = all[o]; r.x[1] = all[o + 1];
      for (int c = 0; c < 9; c++) { r.f[c] = all[o + 2 + c]; }
      out.push_back(r);
   }
   // order by position rounded to 1 mm, then by the frame (the coincident
   // end nodes of two adjacent segments carry the two faces' frames)
   std::sort(out.begin(), out.end(), [](const Row &a, const Row &b) {
      for (int d = 0; d < 2; d++)
      {
         const long long ad = std::llround(a.x[d] * 1e3), bd = std::llround(b.x[d] * 1e3);
         if (ad != bd) { return ad < bd; }
      }
      return std::lexicographical_compare(a.f, a.f + 9, b.f, b.f + 9);
   });
   return out;
}

// `part` (optional): the element partition; `tag` names the test (D1, D2)
static void TestLayout(Mesh &serial, MPI_Comm comm, int np, const char *tag = "D1",
                       const std::vector<int> *part = nullptr)
{
   if (g_rank == 0) { std::cout << "\n=== " << tag << " the 2D fault DOF layout ===\n"; }
   const real_t ref[3] = {1, 0, 0}, upB[3] = {0, 1, 0};
   for (int p = 1; p <= 3; p++)
   {
      ParMesh pm(comm, serial, part ? part->data() : nullptr);
      auto op = MakeOp(pm, p, Source(ref, upB));
      const int nf = op->FrameField()->Report().fault.faces;
      TEST_ASSERT(op->GetNbfPerFace() == p + 1, tag << " p = " << p << ": " << op->GetNbfPerFace()
                  << " nodes per fault segment");
      const int nof = ISum(op->GetNumOwnedFaultFaces(), comm);
      const int nsh = ISum(op->GetFaultSharedFaces().Size(), comm);
      const int nod = ISum(op->GetNumOwnedFaultDOFs(), comm);
      TEST_ASSERT(nof == nf && nod == nf * (p + 1), tag << " p = " << p << ": each of the " << nf
                  << " fault segments owned once (" << nod << " owned DOFs; " << nsh
                  << " rank-boundary copies)");
      // owned -> local scatter of f(x)
      Vector X, Xo;
      op->GetFaultDOFCoords3D(X);
      op->RestrictToOwnedFault(X, Xo, 3);
      auto f = [](real_t x, real_t y) { return 1.0 + x / 7.0e3 - 2.0 * y / 9.0e3 + (x * y) / 5.0e7; };
      Vector fo(op->GetNumOwnedFaultDOFs()), fl;
      for (int k = 0; k < fo.Size(); k++) { fo(k) = f(Xo(3 * k), Xo(3 * k + 1)); }
      op->ExpandOwnedToLocalFault(fo, fl, 1);
      real_t es = 0.0;
      for (int k = 0; k < op->GetNumFaultDOFs(); k++)
      {
         es = std::max(es, std::abs(fl(k) - f(X(3 * k), X(3 * k + 1))));
      }
      es = RMax(es, comm);
      TEST_ASSERT(es <= 1e-12, tag << " p = " << p << ": the owned -> local scatter puts f(x) on every "
                  "local DOF, shared segments included (" << es << ")");
      // the table against a 1-rank run
      std::vector<Row> all = Gather(DofTable(*op), comm);
      bool same = true;
      real_t ex = 0.0;
      if (g_rank == 0)
      {
         ParMesh pm1(MPI_COMM_SELF, serial);
         auto op1 = MakeOp(pm1, p, Source(ref, upB));
         std::vector<Row> one = Gather(DofTable(*op1), MPI_COMM_SELF);
         same = (one.size() == all.size());
         for (size_t i = 0; same && i < one.size(); i++)
         {
            ex = std::max({ex, std::abs(one[i].x[0] - all[i].x[0]), std::abs(one[i].x[1] - all[i].x[1])});
            for (int c = 0; c < 9; c++) { if (one[i].f[c] != all[i].f[c]) { same = false; } }
         }
      }
      MPI_Bcast(&ex, 1, MPI_DOUBLE, 0, comm);
      TEST_ASSERT(AllTrue(same, comm) && ex <= 1e-6, tag << " p = " << p << ": the table (position, n, d, s) "
                  "at " << np << " rank(s) equals the 1-rank one (frames bit for bit, positions "
                  << ex << " m)");
   }
}

// ------------------------------------------------------------- M1 (files)
// The committed 2D meshes (Gmsh 4.15.2, msh 2.2, from the .geo files next to
// them).  The expected counts are those of the files, counted outside MFEM
// (the .msh read in Python, each segment's adjacent triangles counted).
struct TagCount { int attr, faces, exterior, interior; };
struct MeshFile
{
   const char *path;
   real_t scale;
   int ne;
   std::vector<TagCount> tags;
   int fault_attr;
   std::set<int> natural, dirichlet, absorbing;
   real_t length;                      // the fault's length [m]
   real_t lo[2], hi[2];                // the fault's bounding box [m]
   real_t ref[3], up[3], n[3];         // the deck's ref_normal and up, the oriented normal
   int creep_faces;                    // BP3: the creep segment (interior tag-5 segments)
};

static std::vector<MeshFile> MeshFiles()
{
   const real_t c30 = std::cos(M_PI / 6.0), c60 = 0.5;
   // BP3: fault from the trace (0, 0) to 40 km down-dip, (40 cos, -40 sin)
   // km; ref_normal (1, 0), up (0, 1): n = (sin, cos), into the hanging wall
   return {
      {"bp3/mesh/bp3_d30_25m.msh", 1000.0, 127843,
       {{1, 242, 242, 0}, {3, 1600, 0, 1600}, {5, 279, 41, 238}}, 3, {1}, {5}, {}, 40e3,
       {0.0, -20e3}, {40e3 * c30, 0.0}, {1, 0, 0}, {0, 1, 0}, {0.5, c30, 0}, 238},
      {"bp3/mesh/bp3_d60_25m.msh", 1000.0, 212659,
       {{1, 210, 210, 0}, {3, 1600, 0, 1600}, {5, 255, 25, 230}}, 3, {1}, {5}, {}, 40e3,
       {0.0, -40e3 * c30}, {40e3 * c60, 0.0}, {1, 0, 0}, {0, 1, 0}, {c30, c60, 0}, 230},
      {"bp3/mesh/bp3_d90_25m.msh", 1000.0, 244958,
       {{1, 210, 210, 0}, {3, 1600, 0, 1600}, {5, 248, 20, 228}}, 3, {1}, {5}, {}, 40e3,
       {0.0, -40e3}, {0.0, 0.0}, {1, 0, 0}, {0, 1, 0}, {1, 0, 0}, 228},
      // TPV205-2D and TPV105-2D: ref_normal (0, 1), up (0, 0, 1) (design rev 3,
      // section 6.1): n = (0, 1)
      {"tpv205/mesh/tpv205_2d_100m.msh", 1.0, 40526,
       {{103, 300, 0, 300}, {105, 136, 136, 0}}, 103, {}, {}, {105}, 30e3,
       {-15e3, 0.0}, {15e3, 0.0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 0}, 0},
      {"tpv105/mesh/tpv105_2d_100m.msh", 1.0, 56652,
       {{3, 604, 0, 604}, {5, 152, 152, 0}}, 3, {}, {}, {5}, 140e3,
       {-70e3, 0.0}, {70e3, 0.0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 0}, 0}};
}

static void TestMeshFile(const MeshFile &m, MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== M1 " << m.path << " ===\n"; }
   std::unique_ptr<Mesh> mesh = cycle::LoadScaledMesh(m.path, m.scale);
   TEST_ASSERT(mesh->Dimension() == 2 && mesh->SpaceDimension() == 2 && mesh->GetNE() == m.ne,
               "M1 Dimension 2, " << mesh->GetNE() << " triangles (" << m.ne << ")");
   const auto census = cycle::CountBoundary(*mesh);
   bool tags = census.size() == m.tags.size();
   for (const TagCount &t : m.tags)
   {
      const auto it = census.find(t.attr);
      const bool ok = it != census.end() && it->second.faces == t.faces
                      && it->second.exterior == t.exterior && it->second.interior == t.interior
                      && it->second.simplex_faces == t.faces;
      tags = tags && ok;
      if (g_rank == 0 && it != census.end())
      {
         std::cout << "    tag " << t.attr << ": " << it->second.faces << " segments ("
                   << it->second.exterior << " exterior, " << it->second.interior << " interior)\n";
      }
   }
   TEST_ASSERT(tags, "M1 tag census as the file (every tag, segments, exterior, interior)");
   // the fault: its length and bounding box
   real_t len = 0.0;
   Array<int> v;
   for (int be = 0; be < mesh->GetNBE(); be++)
   {
      if (mesh->GetBdrAttribute(be) != m.fault_attr) { continue; }
      mesh->GetBdrElementVertices(be, v);
      const real_t *a = mesh->GetVertex(v[0]), *b = mesh->GetVertex(v[1]);
      len += std::hypot(a[0] - b[0], a[1] - b[1]);
   }
   const cycle::BdrAttrCount &fc = census.at(m.fault_attr);
   const real_t box = std::max({std::abs(fc.lo[0] - m.lo[0]), std::abs(fc.lo[1] - m.lo[1]),
                                std::abs(fc.hi[0] - m.hi[0]), std::abs(fc.hi[1] - m.hi[1])});
   TEST_ASSERT(std::abs(len - m.length) <= 1e-9 * m.length && box <= 1e-9,
               "M1 fault (tag " << m.fault_attr << "): length " << len << " m, bounding box ("
               << fc.lo[0] << ", " << fc.lo[1] << ") - (" << fc.hi[0] << ", " << fc.hi[1]
               << ") m (off by " << box << " m)");

   // the census at this rank count and on 1 rank (P1)
   ParMesh pm(comm, *mesh);
   auto op = MakeOp(pm, 1, Source(m.ref, m.up, m.fault_attr), m.dirichlet, DGMethod::IP,
                    m.fault_attr, m.natural, m.absorbing);
   auto F = op->FrameField();
   const FaultFrameReport &R = F->Report();
   const int nf = census.at(m.fault_attr).faces;
   TEST_ASSERT(R.fault.faces == nf && R.fault.components == 1 && R.fault.junction_edges == 0
               && R.fault.tip_edges == 2 && R.fault.parity_global_ok && R.fault.parity_per_patch_ok,
               "M1 census: " << R.fault.faces << " fault segments, " << R.fault.components
               << " patch, " << R.fault.junction_edges << " junction vertices, " << R.fault.tip_edges
               << " tips, J3 parity");
   TEST_ASSERT(R.creep.faces == m.creep_faces && (m.creep_faces == 0 ||
               (R.creep.components == 1 && R.seam_edges == 1)),
               "M1 creep: " << R.creep.faces << " segments, " << R.creep.components << " patch, "
               << R.seam_edges << " seam vertex");
   real_t en = 0.0;
   bool side_ok = true, slots = true;
   const bool in_plane_up = (m.up[2] == 0.0);
   for (int fi = 0; fi < F->NumFaultFaces(); fi++)
   {
      en = std::max(en, Dist3(F->Normal(fi), m.n));
      const ff::Frame &fr = F->Frame(fi);
      // up in the plane: s = -+e_z (slip in the dip slot); up = e_z: d = -+e_z
      // (slip in the strike slot); the other tangent exactly in the plane
      const real_t *oz = in_plane_up ? fr.s : fr.d, *ip = in_plane_up ? fr.d : fr.s;
      if (!(oz[0] == 0.0 && oz[1] == 0.0 && std::abs(oz[2]) == 1.0 && ip[2] == 0.0
            && fr.n[2] == 0.0 && !fr.fallback)) { slots = false; }
      const int pe = F->FaultCensus().PlusElement(fi);
      if (pe >= 0)
      {
         real_t ce[3], cf[3];
         Centroid(pm, pe, ce);
         FaceCentre(pm, *op, fi, cf);
         real_t mg = 0.0;
         for (int d = 0; d < 3; d++) { mg += (ce[d] - cf[d]) * F->Normal(fi)[d]; }
         if (!(mg < 0.0)) { side_ok = false; }
      }
   }
   en = RMax(en, comm);
   // the file's vertices carry 16 digits and Gmsh's rounding along the
   // line: a few 1e-12 m 40 km out, so a few 1e-13 on a 25 m segment's
   // normal (6.4e-13 seen at dip 60)
   TEST_ASSERT(en <= 4e-12, "M1 every oriented normal = (" << m.n[0] << ", " << m.n[1]
               << ", 0) (" << en << ")");
   TEST_ASSERT(AllTrue(side_ok, comm), "M1 the plus element is on the side n points away from");
   TEST_ASSERT(AllTrue(slots, comm), "M1 up = (" << m.up[0] << ", " << m.up[1] << ", " << m.up[2]
               << "): the slip in the " << (in_plane_up ? "dip" : "strike")
               << " slot, the out-of-plane slot exact");
   ReportKey k1{};
   if (g_rank == 0)
   {
      ParMesh pm1(MPI_COMM_SELF, *mesh);
      auto op1 = MakeOp(pm1, 1, Source(m.ref, m.up, m.fault_attr), m.dirichlet, DGMethod::IP,
                        m.fault_attr, m.natural, m.absorbing);
      k1 = KeyOf(op1->FrameField()->Report());
   }
   MPI_Bcast(&k1, sizeof(ReportKey), MPI_BYTE, 0, comm);
   const ReportKey k = KeyOf(R);
   TEST_ASSERT(k.faces == k1.faces && k.comps == k1.comps && k.junc == k1.junc && k.tips == k1.tips
               && k.cfaces == k1.cfaces && k.ccomps == k1.ccomps && k.seams == k1.seams
               && k.pin == k1.pin && std::abs(k.side - k1.side) <= 1e-12 * k1.side,
               "M1 the report equals the 1-rank report (P1; pin margin " << k.pin
               << ", side margin " << k.side << ")");
}

// ------------------------------------------------------------------- Q1
// The plate loading of cycle/loading.hpp on the miniature: the plane is the
// fault's (normal n = (sin 60, cos 60), through the origin), s = up-dip;
// u = +-(Vp t / 2) s off the plane, Vp t s across the creep segment.
static void PlateLoading2D(const Vector &x, real_t t, Vector &u)
{
   const real_t n[2] = {std::sqrt(3.0) / 2.0, 0.5}, sd[2] = {-0.5, std::sqrt(3.0) / 2.0};
   const real_t Vp = 1e-9, d = x(0) * n[0] + x(1) * n[1];
   const real_t Vh = (d > 1.0) ? 0.5 * Vp * t : (d < -1.0) ? -0.5 * Vp * t : Vp * t;
   u.SetSize(2);
   u(0) = Vh * sd[0];
   u(1) = Vh * sd[1];
}

static void TestQDOperators2D(MPI_Comm comm)
{
   const real_t ref[3] = {1, 0, 0}, upB[3] = {0, 1, 0};
   const real_t t = 100.0 * 365.25 * 86400.0;
   for (int p : {1, 2})
   {
      if (g_rank == 0) { std::cout << "\n=== Q1 stage operators and G in 2D, p = " << p << " ===\n"; }
      std::unique_ptr<Mesh> s(ShearBox(60.0));
      ParMesh pm(comm, *s);
      LinearElastic *model = new LinearElastic(kLambda, kMu);   // leaked: test lifetime
      DomainConfig dc;
      dc.fault_frames = Source(ref, upB);
      BoundaryConfig bc = MakeBC(2, {5});
      bc.default_dirichlet_func = PlateLoading2D;
      Op op(pm, p, *model, 0.0, kWf, 2.0 * kLs, bc, DGMethod::IP, SolverType::MUMPS, dc);
      op.BuildQDStageOperators(true);
      op.BuildQDGreen(QDGreenConfig{});   // dense; its own check against SolveStage (1e-8)
      const int nown = op.GetNumOwnedFaultDOFs(), nc = op.NumSlipComponents();
      for (bool strike : {false, true})
      {
         // dip slip of order 1 m on every DOF (the face loops skip a face
         // whose slip is below 1e-15 m everywhere)
         Vector slip(nc * nown);
         for (int i = 0; i < nown; i++)
         {
            slip(nc * i) = 1.0 + 0.5 * std::sin(0.37 * i + 0.1 + g_rank);
            slip(nc * i + 1) = strike ? 0.3 * std::cos(0.11 * i + g_rank) : 0.0;
         }
         Vector local, lt, lnt, T1, Tn1, T2, Tn2, T3, Tn3;
         ParGridFunction u1(&op.GetFESpace()), u2(&op.GetFESpace());
         op.ExpandOwnedToLocalFault(slip, local, nc);
         op.Solve(t, local, u1);
         op.ComputeTraction(u1, local, lt, &lnt);
         op.RestrictToOwnedFault(lt, T1, nc);
         op.RestrictToOwnedFault(lnt, Tn1);
         op.SolveStage(t, slip, u2, T2, &Tn2);
         op.GreenStage(t, slip, T3, &Tn3);
         auto rel = [&](const Vector &a, const Vector &b)
         {
            real_t d = 0.0, m = 0.0;
            for (int i = 0; i < a.Size(); i++)
            {
               d = std::max(d, std::abs(a(i) - b(i)));
               m = std::max(m, std::abs(b(i)));
            }
            d = RMax(d, comm);
            m = RMax(m, comm);
            return m > 0.0 ? d / m : d;
         };
         real_t ts = 0.0;   // the strike traction (exactly zero in 2D)
         for (int i = 0; i < nown; i++)
         {
            ts = std::max({ts, std::abs(T1(nc * i + 1)), std::abs(T2(nc * i + 1)),
                           std::abs(T3(nc * i + 1))});
         }
         ts = RMax(ts, comm);
         const std::string tag = std::string("Q1 p = ") + std::to_string(p)
                                 + (strike ? ", strike slot set" : ", dip slip only");
         const real_t r_ms = rel(T2, T1), r_ns = rel(Tn2, Tn1);
         const real_t r_gs = rel(T3, T1), r_gn = rel(Tn3, Tn1);
         TEST_ASSERT(r_ms <= 1e-10 && r_ns <= 1e-10, tag << ": stage operators = face loops "
                     "(shear " << r_ms << ", normal " << r_ns << ")");
         TEST_ASSERT(r_gs <= 1e-10 && r_gn <= 1e-10, tag << ": G = face loops (shear "
                     << r_gs << ", normal " << r_gn << ")");
         TEST_ASSERT(ts == 0.0, tag << ": strike traction exactly 0 on every path");
      }
   }
}

// --------------------------------------------------------------- refusals
static int Refuse(int n, MPI_Comm comm)
{
   const real_t ref[3] = {1, 0, 0}, upB[3] = {0, 1, 0};
   switch (n)
   {
      case 1:
      {
         std::unique_ptr<Mesh> q(ShearBox(60.0, 1));
         const std::string path = "test_2d_fault_layout_quad.mesh";   // the make target removes it
         if (g_rank == 0) { std::ofstream ofs(path); q->Print(ofs); }
         MPI_Barrier(comm);
         cycle::LoadScaledMesh(path, 1.0);
         break;
      }
      case 2:
      {
         std::unique_ptr<Mesh> s(ShearBox(60.0));
         ParMesh pm(comm, *s);
         MakeOp(pm, 1, FaultFrameSource{});
         break;
      }
      case 3:
      case 4:
      {
         std::unique_ptr<Mesh> s(ShearBox(60.0));
         ParMesh pm(comm, *s);
         const real_t up_bad[3] = {0, 0.6, 0.8}, ref_bad[3] = {0.8, 0, 0.6};
         MakeOp(pm, 1, (n == 3) ? Source(ref, up_bad) : Source(ref_bad, upB));
         break;
      }
      case 5:
      {
         std::unique_ptr<Mesh> s(ShearBox(60.0));
         ParMesh pm(comm, *s);
         auto op = MakeOp(pm, 1, Source(ref, upB));
         op->BuildQDStageOperators(false);
         QDGreenConfig gc;
         gc.cache_path = "test_2d_fault_layout_g.bin";   // refused before it is written
         op->BuildQDGreen(gc);
         break;
      }
      case 6:
      {
         std::unique_ptr<Mesh> y(YMesh());
         ParMesh pm(comm, *y);
         const real_t refY[3] = {0, 1, 0};
         MakeOp(pm, 1, Source(refY, upB));
         break;
      }
      case 7:
      {
         std::unique_ptr<Mesh> s(ShearBox(60.0));
         ParMesh pm(comm, *s);
         MakeOp(pm, 1, Source(ref, upB), {5}, DGMethod::BR2);
         break;
      }
      default: return 2;
   }
   if (g_rank == 0) { std::cout << "refusal " << n << ": NOT refused\n"; }
   return 0;   // reaching here means the refusal did not fire
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
      std::cout << "test_2d_fault_layout (2D in-plane design rev 3, phase B) on " << np
                << " rank(s)\n";
   }
   TestLoader(comm);
   for (real_t dip : {30.0, 60.0, 90.0})
   {
      if (g_rank == 0) { std::cout << "\n=== C1 the BP3 layout, dip " << dip << " ===\n"; }
      std::unique_ptr<Mesh> s(ShearBox(dip));
      ReportKey k1{};
      if (g_rank == 0) { k1 = Report1(*s); }
      MPI_Bcast(&k1, sizeof(ReportKey), MPI_BYTE, 0, comm);
      TestBP3Census(*s, dip, comm, &k1);
   }
   {
      std::unique_ptr<Mesh> s(ShearBox(60.0));
      TestSlots(*s, comm);
   }
   TestY(comm);
   TestFallback2D();
   TestSlice(comm);
   {
      std::unique_ptr<Mesh> s(ShearBox(60.0));
      TestLayout(*s, comm, np);
   }
   TestQDOperators2D(comm);
   for (const MeshFile &m : MeshFiles()) { TestMeshFile(m, comm); }
   {
      // D2: the layout on the BP3 mesh (dip 60, 25 m), with METIS's
      // partition (which keeps the fault inside the parts) and with one that
      // cuts along the fault: the hanging wall and the footwall on different
      // ranks, each side split at 20 km down-dip on 4 ranks, so that every
      // fault segment is a rank-boundary face
      std::unique_ptr<Mesh> s = cycle::LoadScaledMesh("bp3/mesh/bp3_d60_25m.msh", 1000.0);
      if (g_rank == 0) { std::cout << "\n(D2: bp3/mesh/bp3_d60_25m.msh)"; }
      TestLayout(*s, comm, np, "D2 (METIS)");
      const real_t n60[2] = {std::sqrt(3.0) / 2.0, 0.5};
      std::vector<int> cut(s->GetNE());
      for (int e = 0; e < s->GetNE(); e++)
      {
         Array<int> ev;
         s->GetElementVertices(e, ev);
         real_t c[2] = {0.0, 0.0};
         for (int a = 0; a < ev.Size(); a++)
            for (int d = 0; d < 2; d++) { c[d] += s->GetVertex(ev[a])[d] / ev.Size(); }
         const int side = (c[0] * n60[0] + c[1] * n60[1] > 0.0) ? 0 : 1;
         const int deep = (c[1] < -20e3 * 0.5 * std::sqrt(3.0)) ? 1 : 0;
         cut[e] = (np == 1) ? 0 : (np == 2) ? side : std::min(np - 1, 2 * side + deep);
      }
      TestLayout(*s, comm, np, "D2 (fault cut)", &cut);
   }
   if (g_rank == 0)
   {
      std::cout << "\ntest_2d_fault_layout: " << (num_tests - num_failed) << " / " << num_tests
                << " passed\n";
   }
   return num_failed == 0 ? 0 : 1;
}
