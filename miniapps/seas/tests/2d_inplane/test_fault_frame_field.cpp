// tests/2d_inplane/test_fault_frame_field.cpp — 2D in-plane design rev 3,
// phase A (A1 and the first half of A6): the one fault orientation system
// (fault/fault_frame_field.hpp) and the QD operator that reads it.
//
// Mesh: a box [-L, L] x [-L, L] x [-L, 0] of 4 x 4 x 2 hexes split into 24
// tets each (768 tets), whose plane y = 0 carries
// the fault (attr 3, |x| < L/2 and z > -L/2: 8 triangles) and a coplanar
// prescribed-slip band around it (attr 6, interior Dirichlet: 24
// triangles); walls and bottom Dirichlet (attr 5), top natural (attr 1).
// The box is used as built and rotated rigidly:
//   flat     no rotation, ref_normal (0, -1, 0), up (0, 0, 1)
//   dip60    rotated 30 degrees about x: the fault dips 60 degrees;
//            ref_normal R (0, -1, 0), up (0, 0, 1) (true vertical)
//   general  rotated 40 degrees about (1, 2, 3); ref_normal and up rotated
//
//   F1  field census: 8 fault and 24 prescribed-slip faces, one patch each,
//       4 seam edges, no junction edges, no F6 frames; every oriented normal
//       = R (0, -1, 0) to 1e-14; frames orthonormal, d x s = n, and on dip60
//       s = (1, 0, 0), d = (0, sin 30, -cos 30) to 1e-14
//   F2  sigma_f = sign(n_raw . n) on every fault face; the QD kernel triad
//       of every fault DOF (GetFaultDOFBasis) = sigma_f (n, d, s) to 1e-13,
//       and the per-face and per-quadrature-point FaultBasis the same
//   F3  flat box: the operator in field mode and in reference mode (no
//       source, the constants) give the same FaultBasis (per face and per
//       quadrature point), the same kernel triad and the same creep census
//       sign bit for bit; and the same traction (bit for bit on 1 rank, 1e-12
//       relative on more: MUMPS)
//   F4  rotation invariance (A6, first half): the fault traction (dip,
//       strike, normal) of a strike slip with the prescribed-slip band
//       loaded by R (0.3, 0, 0.1) is the same on dip60 and general as on
//       flat to 1e-9 relative
//   F5  a second operator handed the field (DomainConfig::fault_frames.field)
//       reads it: the same kernel triad bit for bit, FrameField() is the
//       same object
//
// Modes: (default) F1-F5;  --seam good|bad: an L-shaped surface (fault on
// y = 0 for -L/2 < x < 0, prescribed slip on x = 0 for 0 < y < L/2) with
// ref_normal (1, -1, 0)/sqrt 2 (consistent: passes) or (-1, -1, 0)/sqrt 2
// (the two censuses orient the seam inconsistently: the field refuses the
// run, exit status != 0, message "seam rule").
//
//   make -f tests/2d_inplane/2d_inplane.mk test-fault-frame-field
//
// Scope (Chunhui Zhao, 2026-10-24).

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
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

static const real_t kL = 10.0e3;
static const real_t kLambda = 32.04e9, kMu = 32.04e9;

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static bool AllTrue(bool b, MPI_Comm c) { int v = b ? 0 : 1; return ISum(v, c) == 0; }

// ---------------------------------------------------------------- rotations
struct Rot
{
   real_t R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
   void Apply(const real_t a[3], real_t b[3]) const
   {
      for (int i = 0; i < 3; i++) { b[i] = R[i][0] * a[0] + R[i][1] * a[1] + R[i][2] * a[2]; }
   }
   // Rodrigues: angle (degrees) about axis.
   static Rot Axis(real_t ax, real_t ay, real_t az, real_t deg)
   {
      const real_t l = std::sqrt(ax * ax + ay * ay + az * az);
      const real_t k[3] = {ax / l, ay / l, az / l};
      const real_t t = deg * M_PI / 180.0, c = std::cos(t), s = std::sin(t);
      Rot r;
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++)
         {
            r.R[i][j] = (i == j ? c : 0.0) + (1 - c) * k[i] * k[j];
         }
      r.R[0][1] -= s * k[2]; r.R[0][2] += s * k[1];
      r.R[1][0] += s * k[2]; r.R[1][2] -= s * k[0];
      r.R[2][0] -= s * k[1]; r.R[2][1] += s * k[0];
      return r;
   }
};

// ------------------------------------------------------------------- meshes
// kind 0: the coplanar box (fault + band on y = 0); kind 1: the L-shaped
// seam surface.  Vertices rotated by R after the faces are tagged.
static Mesh *FrameBox(int kind, const Rot &rot)
{
   const int nx = 4, ny = 4, nz = (kind == 0) ? 2 : 1;
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(nx, ny, nz, 2 * kL, 2 * kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   const real_t tol = 1e-6 * kL;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < tol) ? 1 : 5);
   }
   struct Tri { int v[3]; int attr; };
   std::vector<Tri> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      real_t c[3] = {0, 0, 0};
      bool y0 = true, x0 = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         const real_t *x = mesh->GetVertex(fv[i]);
         for (int d = 0; d < 3; d++) { c[d] += x[d] / fv.Size(); }
         if (std::abs(x[1]) > tol) { y0 = false; }
         if (std::abs(x[0]) > tol) { x0 = false; }
      }
      int attr = 0;
      if (kind == 0 && y0)
      {
         attr = (std::abs(c[0]) < 0.5 * kL && c[2] > -0.5 * kL) ? 3 : 6;
      }
      else if (kind == 1)
      {
         if (y0 && c[0] > -0.5 * kL && c[0] < 0.0) { attr = 3; }
         if (x0 && c[1] > 0.0 && c[1] < 0.5 * kL) { attr = 6; }
      }
      if (attr) { tris.push_back({{fv[0], fv[1], fv[2]}, attr}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t.v[0], t.v[1], t.v[2], t.attr); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      const real_t x[3] = {c[0], c[1], c[2]};
      rot.Apply(x, c);
   }
   return mesh;
}

static BoundaryConfig MakeBC(const Rot &rot)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5, 6};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   // the prescribed-slip band: a constant jump, rotated with the box
   const real_t j0[3] = {0.3, 0.0, 0.1};
   real_t j[3];
   rot.Apply(j0, j);
   const std::array<real_t, 3> ja = {{j[0], j[1], j[2]}};
   bc.dirichlet_funcs[6] = [ja](const Vector &, real_t, Vector &u)
   { u.SetSize(3); for (int d = 0; d < 3; d++) { u(d) = ja[d]; } };
   return bc;
}

static std::shared_ptr<const FaultFrameConfig> FrameConfig(const real_t ref[3], const real_t up[3])
{
   auto c = std::make_shared<FaultFrameConfig>();
   for (int d = 0; d < 3; d++)
   {
      c->orientation.ref_normal[d] = ref[d];
      c->orientation.up[d] = up[d];
      c->up[d] = up[d];
   }
   c->orientation.fault_attrs = {3};
   return c;
}

using Op = ElasticityDomainOperator<ParMesh>;

static std::unique_ptr<Op> MakeOp(ParMesh &pm, const Rot &rot, FaultFrameSource src)
{
   LinearElastic *model = new LinearElastic(kLambda, kMu);   // leaked: test lifetime
   DomainConfig dc;
   dc.fault_frames = src;
   return std::make_unique<Op>(pm, 1, *model, 0.0, kL, 2.0 * kL, MakeBC(rot), DGMethod::IP,
                               SolverType::MUMPS, dc);
}

static real_t Dist3(const real_t a[3], const real_t b[3])
{
   return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
                    + (a[2] - b[2]) * (a[2] - b[2]));
}

// ---------------------------------------------------------------- F1, F2
static void TestCensusAndTriad(ParMesh &pm, const Rot &rot, const real_t up[3], bool dip60,
                               const char *label)
{
   MPI_Comm comm = pm.GetComm();
   if (g_rank == 0) { std::cout << "\n=== F1/F2 on " << label << " ===\n"; }
   const real_t ref0[3] = {0, -1, 0};
   real_t ref[3];
   rot.Apply(ref0, ref);
   FaultFrameSource src;
   src.config = FrameConfig(ref, up);
   auto op = MakeOp(pm, rot, src);
   auto F = op->FrameField();
   TEST_ASSERT(F != nullptr, "the operator built a FaultFrameField");
   if (!F) { return; }
   const FaultFrameReport &R = F->Report();
   TEST_ASSERT(R.fault.faces == 8 && R.creep.faces == 24,
               "F1 8 fault faces and 24 prescribed-slip faces (" << R.fault.faces << ", "
               << R.creep.faces << ")");
   TEST_ASSERT(R.fault.components == 1 && R.creep.components == 1,
               "F1 one patch each (" << R.fault.components << ", " << R.creep.components << ")");
   TEST_ASSERT(R.seam_edges == 4, "F1 4 seam edges (" << R.seam_edges << ")");
   TEST_ASSERT(R.fault.junction_edges == 0 && R.frame_fallbacks == 0,
               "F1 no junction edges, no F6 frames");

   real_t en = 0.0, eorth = 0.0, edip = 0.0, ecreep = 0.0;
   const real_t s_exp[3] = {1, 0, 0};
   const real_t d_exp[3] = {0, std::sin(M_PI / 6), -std::cos(M_PI / 6)};
   for (int fi = 0; fi < F->NumFaultFaces(); fi++)
   {
      const ff::Frame &fr = F->Frame(fi);
      en = std::max(en, Dist3(F->Normal(fi), ref));
      real_t dxs[3];
      ff::Cross(fr.d, fr.s, dxs);
      eorth = std::max({eorth, Dist3(dxs, fr.n), std::abs(ff::Dot(fr.d, fr.s)),
                        std::abs(ff::Dot(fr.d, fr.d) - 1), std::abs(ff::Dot(fr.s, fr.s) - 1)});
      if (dip60) { edip = std::max({edip, Dist3(fr.s, s_exp), Dist3(fr.d, d_exp)}); }
   }
   const FaultOrientation *C = F->CreepCensus();
   if (C)
   {
      for (int ci = 0; ci < C->NumLocalFaces(); ci++) { ecreep = std::max(ecreep, Dist3(C->Normal(ci), ref)); }
   }
   en = RMax(en, comm); eorth = RMax(eorth, comm); edip = RMax(edip, comm);
   ecreep = RMax(ecreep, comm);
   TEST_ASSERT(en <= 1e-14 && ecreep <= 1e-14, "F1 every oriented normal = R (0,-1,0) (fault "
               << en << ", prescribed slip " << ecreep << ")");
   TEST_ASSERT(eorth <= 1e-14, "F1 frames orthonormal, d x s = n (" << eorth << ")");
   if (dip60)
   {
      TEST_ASSERT(edip <= 1e-14, "F1 dip60: s = (1,0,0), d = (0, sin 30, -cos 30) (" << edip << ")");
   }

   // F2: sigma_f and the kernel triad
   const Array<int> &fint = op->GetFaultInteriorFaces();
   const Array<int> &fsh = op->GetFaultSharedFaces();
   bool sig_ok = true;
   for (int fi = 0; fi < F->NumFaultFaces(); fi++)
   {
      const bool sh = fi >= fint.Size();
      FaceElementTransformations *Tr = sh ? pm.GetSharedFaceTransformations(fsh[fi - fint.Size()])
                                          : pm.GetInteriorFaceTransformations(fint[fi]);
      const IntegrationPoint &ipc = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetAllIntPoints(&ipc);
      Vector nr(3);
      CalcOrtho(Tr->Jacobian(), nr);
      const real_t s = (ff::Dot(nr.GetData(), F->Normal(fi)) > 0) ? 1.0 : -1.0;
      if (s != F->Sigma(fi)) { sig_ok = false; }
   }
   TEST_ASSERT(AllTrue(sig_ok, comm), "F2 sigma_f = sign(n_raw . n) on every fault face");

   DenseMatrix Bk;
   op->GetFaultDOFBasis(Bk);
   const int nbf = op->GetNbfPerFace();
   real_t etri = 0.0, eface = 0.0, eqp = 0.0;
   const FaultBasis *fb = op->GetFaultBasis();
   for (int k = 0; k < op->GetNumFaultDOFs(); k++)
   {
      const int fi = k / nbf;
      const ff::Frame &fr = F->Frame(fi);
      const real_t sg = F->Sigma(fi);
      for (int d = 0; d < 3; d++)
      {
         etri = std::max({etri, std::abs(Bk(d, k) - sg * fr.n[d]), std::abs(Bk(3 + d, k) - sg * fr.d[d]),
                          std::abs(Bk(6 + d, k) - sg * fr.s[d])});
      }
   }
   for (int fi = 0; fi < F->NumFaultFaces(); fi++)
   {
      const FaultBasisData &b = fb->GetBasis(fi);
      const ff::Frame &fr = F->Frame(fi);
      const real_t sg = F->Sigma(fi);
      for (int d = 0; d < 3; d++)
      {
         eface = std::max({eface, std::abs(b.normal[d] - sg * fr.n[d]),
                           std::abs(b.tangent1[d] - sg * fr.d[d]), std::abs(b.tangent2[d] - sg * fr.s[d])});
      }
      for (const auto &q : b.qp_data)
      {
         for (int d = 0; d < 3; d++)
         {
            eqp = std::max({eqp, std::abs(q.normal[d] - sg * fr.n[d]),
                            std::abs(q.tangent1[d] - sg * fr.d[d]), std::abs(q.tangent2[d] - sg * fr.s[d])});
         }
      }
   }
   etri = RMax(etri, comm); eface = RMax(eface, comm); eqp = RMax(eqp, comm);
   TEST_ASSERT(etri <= 1e-13, "F2 kernel triad = sigma_f (n, d, s) on every fault DOF (" << etri << ")");
   TEST_ASSERT(eface <= 1e-13 && eqp <= 1e-13, "F2 FaultBasis per face and per QP = sigma_f (n, d, s) ("
               << eface << ", " << eqp << ")");

   // F5: a second operator handed the field
   FaultFrameSource given;
   given.field = F;
   auto op2 = MakeOp(pm, rot, given);
   DenseMatrix Bk2;
   op2->GetFaultDOFBasis(Bk2);
   const bool same = Bk.Height() == Bk2.Height() && Bk.Width() == Bk2.Width()
                     && (Bk.Height() * Bk.Width() == 0
                         || std::memcmp(Bk.Data(), Bk2.Data(), sizeof(real_t) * Bk.Height() * Bk.Width()) == 0);
   TEST_ASSERT(AllTrue(same && op2->FrameField() == F, comm),
               "F5 an operator handed the field reads it (same triad bit for bit, same object)");
}

// ------------------------------------------------------------------- F3, F4
struct TracResult { Vector t, tn, x; };

static TracResult SolveTraction(Op &op)
{
   const int nd = op.GetNumFaultDOFs();
   Vector slip(2 * nd);
   for (int k = 0; k < nd; k++) { slip(2 * k) = 0.0; slip(2 * k + 1) = 1.0; }   // strike slip
   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);
   TracResult r;
   op.ComputeTraction(u, slip, r.t, &r.tn);
   op.GetFaultDOFCoords3D(r.x);
   return r;
}

static bool SameBits(const Vector &a, const Vector &b)
{
   return a.Size() == b.Size()
          && (a.Size() == 0 || std::memcmp(a.GetData(), b.GetData(), sizeof(real_t) * a.Size()) == 0);
}

static void TestFlatParity(ParMesh &pm, int nranks)
{
   MPI_Comm comm = pm.GetComm();
   if (g_rank == 0) { std::cout << "\n=== F3 flat box: field mode against the constants ===\n"; }
   Rot I;
   const real_t ref[3] = {0, -1, 0}, up[3] = {0, 0, 1};
   auto opr = MakeOp(pm, I, FaultFrameSource{});
   FaultFrameSource src;
   src.config = FrameConfig(ref, up);
   auto opf = MakeOp(pm, I, src);
   TEST_ASSERT(opr->FrameField() == nullptr && opf->FrameField() != nullptr,
               "F3 no source: reference mode; config: field mode");

   bool basis_same = true;
   const FaultBasis *a = opr->GetFaultBasis(), *b = opf->GetFaultBasis();
   basis_same = (a->NumFaces() == b->NumFaces());
   for (int fi = 0; basis_same && fi < a->NumFaces(); fi++)
   {
      const FaultBasisData &x = a->GetBasis(fi), &y = b->GetBasis(fi);
      basis_same = std::memcmp(x.normal, y.normal, sizeof x.normal) == 0
                   && std::memcmp(x.tangent1, y.tangent1, sizeof x.tangent1) == 0
                   && std::memcmp(x.tangent2, y.tangent2, sizeof x.tangent2) == 0
                   && x.sign_flipped == y.sign_flipped && x.qp_data.size() == y.qp_data.size();
      for (size_t q = 0; basis_same && q < x.qp_data.size(); q++)
      {
         const auto &p = x.qp_data[q], &r = y.qp_data[q];
         basis_same = std::memcmp(p.normal, r.normal, sizeof p.normal) == 0
                      && std::memcmp(p.tangent1, r.tangent1, sizeof p.tangent1) == 0
                      && std::memcmp(p.tangent2, r.tangent2, sizeof p.tangent2) == 0
                      && p.nl == r.nl && p.sign_flipped == r.sign_flipped;
      }
   }
   TEST_ASSERT(AllTrue(basis_same, comm), "F3 FaultBasis per face and per QP bit for bit");
   DenseMatrix Ba, Bb;
   opr->GetFaultDOFBasis(Ba);
   opf->GetFaultDOFBasis(Bb);
   const bool tri_same = Ba.Height() * Ba.Width() == Bb.Height() * Bb.Width()
      && (Ba.Height() * Ba.Width() == 0
          || std::memcmp(Ba.Data(), Bb.Data(), sizeof(real_t) * Ba.Height() * Ba.Width()) == 0);
   TEST_ASSERT(AllTrue(tri_same, comm), "F3 kernel triad bit for bit");

   // the prescribed-slip sign: sign(n_raw . (0,-1,0)) against sign(n_raw . n_creep)
   auto F = opf->FrameField();
   bool csign = true;
   int ncreep = 0;
   for (int lf = 0; lf < pm.GetNumFaces(); lf++)
   {
      if (!F->IsCreepFace(lf)) { continue; }
      ncreep++;
      FaceElementTransformations *Tr = pm.GetFaceInformation(lf).IsShared()
         ? pm.GetSharedFaceTransformationsByLocalIndex(lf) : pm.GetInteriorFaceTransformations(lf);
      const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 3);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         Tr->SetAllIntPoints(&ir.IntPoint(q));
         Vector nr(3);
         CalcOrtho(Tr->Jacobian(), nr);
         const real_t s_ref = (-nr(1) < 0.0) ? -1.0 : 1.0;   // (0*nr0 + -1*nr1 + 0*nr2)
         const real_t *c = F->CreepNormal(lf);
         const real_t dot = nr(0) * c[0] + nr(1) * c[1] + nr(2) * c[2];
         const real_t s_fld = (dot < 0.0) ? -1.0 : 1.0;
         if (s_ref != s_fld) { csign = false; }
      }
   }
   TEST_ASSERT(AllTrue(csign, comm) && ISum(ncreep, comm) > 0,
               "F3 prescribed-slip sign: the creep census gives the constants' sign at every QP");

   TracResult ra = SolveTraction(*opr), rb = SolveTraction(*opf);
   if (nranks == 1)
   {
      TEST_ASSERT(SameBits(ra.t, rb.t) && SameBits(ra.tn, rb.tn), "F3 traction bit for bit (1 rank)");
   }
   else
   {
      real_t e = 0.0, s = 0.0;
      for (int i = 0; i < ra.t.Size(); i++) { e = std::max(e, std::abs(ra.t(i) - rb.t(i))); s = std::max(s, std::abs(ra.t(i))); }
      for (int i = 0; i < ra.tn.Size(); i++) { e = std::max(e, std::abs(ra.tn(i) - rb.tn(i))); s = std::max(s, std::abs(ra.tn(i))); }
      e = RMax(e, comm); s = RMax(s, comm);
      TEST_ASSERT(e <= 1e-12 * s, "F3 traction to 1e-12 relative (" << e / s << ", MUMPS on "
                  << nranks << " ranks)");
   }
}

static void TestRotation(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== F4 rotation invariance of the fault traction ===\n"; }
   Rot I;
   const real_t ref0[3] = {0, -1, 0}, up0[3] = {0, 0, 1};
   std::unique_ptr<Mesh> m0(FrameBox(0, I));
   ParMesh p0(comm, *m0);
   FaultFrameSource s0;
   s0.config = FrameConfig(ref0, up0);
   auto op0 = MakeOp(p0, I, s0);
   TracResult r0 = SolveTraction(*op0);
   real_t scale = 0.0;
   for (int i = 0; i < r0.t.Size(); i++) { scale = std::max(scale, std::abs(r0.t(i))); }
   for (int i = 0; i < r0.tn.Size(); i++) { scale = std::max(scale, std::abs(r0.tn(i))); }
   scale = RMax(scale, comm);

   struct Case { const char *name; Rot R; bool rotate_up; };
   const Case cases[2] = {{"dip60", Rot::Axis(1, 0, 0, 30.0), false},
                          {"general", Rot::Axis(1, 2, 3, 40.0), true}};
   for (const Case &c : cases)
   {
      std::unique_ptr<Mesh> m(FrameBox(0, c.R));
      ParMesh pm(comm, *m);
      real_t ref[3], up[3];
      c.R.Apply(ref0, ref);
      if (c.rotate_up) { c.R.Apply(up0, up); } else { for (int d = 0; d < 3; d++) { up[d] = up0[d]; } }
      FaultFrameSource s;
      s.config = FrameConfig(ref, up);
      auto op = MakeOp(pm, c.R, s);
      TracResult r = SolveTraction(*op);
      bool layout = (r.t.Size() == r0.t.Size() && r.tn.Size() == r0.tn.Size());
      real_t ex = 0.0, et = 0.0;
      for (int k = 0; layout && k < r0.tn.Size(); k++)
      {
         real_t xr[3];
         c.R.Apply(&r0.x(3 * k), xr);
         ex = std::max(ex, Dist3(xr, &r.x(3 * k)));
      }
      for (int i = 0; layout && i < r0.t.Size(); i++) { et = std::max(et, std::abs(r.t(i) - r0.t(i))); }
      for (int i = 0; layout && i < r0.tn.Size(); i++) { et = std::max(et, std::abs(r.tn(i) - r0.tn(i))); }
      ex = RMax(ex, comm); et = RMax(et, comm);
      TEST_ASSERT(AllTrue(layout, comm) && ex <= 1e-9 * kL,
                  "F4 " << c.name << ": the same fault DOFs, at R x (" << ex << " m)");
      TEST_ASSERT(et <= 1e-9 * scale, "F4 " << c.name << ": traction (dip, strike, normal) = flat's ("
                  << std::scientific << std::setprecision(2) << et / scale << " relative, scale "
                  << scale << " Pa)" << std::defaultfloat);
   }
}

// --------------------------------------------------------------------- seam
static int RunSeam(MPI_Comm comm, bool good)
{
   Rot I;
   const real_t s = 1.0 / std::sqrt(2.0);
   const real_t ref[3] = {good ? s : -s, -s, 0.0}, up[3] = {0, 0, 1};
   std::unique_ptr<Mesh> m(FrameBox(1, I));
   ParMesh pm(comm, *m);
   FaultFrameSource src;
   src.config = FrameConfig(ref, up);
   if (g_rank == 0)
   {
      std::cout << "\n=== seam " << (good ? "good" : "bad") << ": ref_normal ("
                << ref[0] << ", " << ref[1] << ", 0) ===\n";
   }
   auto op = MakeOp(pm, I, src);   // the bad case refuses here
   auto F = op->FrameField();
   TEST_ASSERT(F && F->Report().seam_edges > 0 && F->Report().fault.faces > 0
               && F->Report().creep.faces > 0,
               "seam good: the field accepts the L-shaped surface (" << (F ? F->Report().seam_edges : -1)
               << " seam edges)");
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   g_rank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();
   std::string seam;
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "--seam") == 0 && i + 1 < argc) { seam = argv[++i]; }
   }
   if (g_rank == 0)
   {
      std::cout << "test_fault_frame_field (2D in-plane design rev 3, phase A) on " << nranks
                << " rank(s)\n";
   }
   if (!seam.empty())
   {
      MFEM_VERIFY(seam == "good" || seam == "bad", "--seam good|bad");
      RunSeam(comm, seam == "good");
   }
   else
   {
      {
         Rot I;
         std::unique_ptr<Mesh> m(FrameBox(0, I));
         ParMesh pm(comm, *m);
         const real_t up[3] = {0, 0, 1};
         TestCensusAndTriad(pm, I, up, false, "flat");
         TestFlatParity(pm, nranks);
      }
      {
         const Rot R = Rot::Axis(1, 0, 0, 30.0);
         std::unique_ptr<Mesh> m(FrameBox(0, R));
         ParMesh pm(comm, *m);
         const real_t up[3] = {0, 0, 1};
         TestCensusAndTriad(pm, R, up, true, "dip60");
      }
      {
         const Rot R = Rot::Axis(1, 2, 3, 40.0);
         std::unique_ptr<Mesh> m(FrameBox(0, R));
         ParMesh pm(comm, *m);
         const real_t up0[3] = {0, 0, 1};
         real_t up[3];
         R.Apply(up0, up);
         TestCensusAndTriad(pm, R, up, false, "general");
      }
      TestRotation(comm);
   }
   if (g_rank == 0)
   {
      std::cout << "\ntest_fault_frame_field: " << (num_tests - num_failed) << " / " << num_tests
                << " passed" << (num_failed ? " -- FAILED" : "") << "\n";
   }
   return num_failed ? 1 : 0;
}
