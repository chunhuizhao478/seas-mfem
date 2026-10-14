// Nitsche dynamic rupture, phase PF: the fault coordinate layer (fault
// coordinate spec rev 8.3; plan N-08, N-09, the Part IV dictionary).
//
//   T-14(a) eps_f = -1 count == the kernel's flipped count; one component,
//           zero junction edges.
//   T-14(b) FromKernel(ToKernel(delta)) == delta on random tangential
//           vectors (1e-15).
//   T-14(c) (T.d, T.s, sigma_n) through the dictionary == the kernel's own
//           traction output bit for bit on every fault DOF (spec T1).
//   T-14(d) d x s == n exactly on every face; TPV102 values; 0 fallbacks.
//   P2      Orient(n_raw) == Orient(-n_raw) exactly (spec T3, last part).
//   T11(a)  tangent-basis invariance at the constitutive level (one
//           evaluation, through the P5 RateStateClosure).
//   T12     isotropy guard (static part).
//
// Modes: (default) inline 24-tet box; --gmsh <file.msh> T-14 on a Gmsh mesh;
//        --dump <prefix> per-face data for the T-16 cross-rank comparison.
//
// Scope (Chunhui Zhao, 2026-10-14).

#include "mfem.hpp"
#include "../../fault/fault_orientation.hpp"
#include "../../fault/fault_frame.hpp"
#include "../../friction/rate_state_closure.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../config/tpv102_params.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
namespace ff = mfem::seas::fault_frame;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }

// 24-tet TPV102 box with the fault on Y=0 (attr 3), free surface (attr 1),
// Dirichlet elsewhere (attr 5).
static Mesh *FaultBox(int nx, int ny, int nz)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2*kL, 2*kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   const real_t tol = kL / nz * 0.1;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < tol) ? 1 : 5);
   }
   const real_t ytol = 1e-6 * kL;
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         if (std::abs(mesh->GetVertex(fv[i])[1]) > ytol) { on_plane = false; }
      }
      if (on_plane) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}

static BoundaryConfig MakeBC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   return bc;
}

// Frame of local fault face fi from the oriented normal (N-09).
static bool FaceFrame(const FaultOrientation &fo, int fi, const real_t up[3],
                      ff::Frame &fr)
{
   ff::FaceVertices fv;
   fv.nv = 3;
   fo.FaceVertexData(fi, fv.key, fv.x);
   return ff::ComputeFrame(fo.Normal(fi), up, ff::DegeneracyThreshold(1.0), fv, fr);
}

// ------------------------------------------------------------- T-14 (a)-(d)
// exact = true: the exactly planar inline box; bit-for-bit gates (spec T1:
// bit-for-bit holds only when the frame changes by pure sign flips).
// exact = false: a Gmsh mesh whose fault is planar only to round-off; T1 is
// gated at round-off (relative <= 1e-14, triad <= 1e-15) and the bitwise
// counts are printed as information (decision, Chunhui Zhao, 2026-10-14).
static void TestT14(ParMesh &pm, int p, const char *label, bool exact)
{
   if (g_rank == 0) { std::cout << "\n=== T-14 on " << label << ", IP p = " << p << " ===\n"; }
   MPI_Comm comm = pm.GetComm();
   LinearElastic model(TPV102Params::lambda, TPV102Params::mu);
   ElasticityDomainOperator<ParMesh> op(pm, p, model, 0.0, kL, 2.0 * kL, MakeBC(),
                                        DGMethod::IP, SolverType::CG_AMG);
   FaultOrientationConfig cfg;   // reference pin (0,-1,0), up (0,0,1)
   FaultOrientation fo(pm, op.GetFaultInteriorFaces(), op.GetFaultSharedFaces(), cfg);
   const FaultOrientationReport &R = fo.GetReport();
   const FaultBasis *fb = op.GetFaultBasis();
   const real_t up[3] = {0, 0, 1};
   const real_t ref[3] = {0, -1, 0};

   // (a) eps_f = sign(n . n_k) and the kernel's flipped count (G0.5).
   int n_eps_neg = 0, n_kernel_flip = 0, n_fallback = 0, n_frame_bad = 0, n_tpv = 0;
   const int nloc = fo.NumLocalFaces();
   const int n_int = op.GetFaultInteriorFaces().Size();
   std::vector<real_t> eps(nloc);
   std::vector<ff::Frame> frames(nloc);
   real_t worst_cross = 0.0, worst_tri = 0.0;
   for (int fi = 0; fi < nloc; fi++)
   {
      const FaultBasisData &b = fb->GetBasis(fi);
      eps[fi] = ff::EpsilonF(fo.Normal(fi), b.normal);
      if (fi < n_int)
      {
         if (eps[fi] < 0) { n_eps_neg++; }
         if (ff::Dot(b.normal, ref) < 0.0) { n_kernel_flip++; }
      }
      // (d) frame
      ff::Frame &fr = frames[fi];
      if (!FaceFrame(fo, fi, up, fr)) { n_frame_bad++; continue; }
      if (fr.fallback) { n_fallback++; }
      real_t c[3];
      ff::Cross(fr.d, fr.s, c);
      for (int d = 0; d < 3; d++) { worst_cross = std::max(worst_cross, std::abs(c[d] - fr.n[d])); }
      const real_t tri[9] = {0, -1, 0, 0, 0, -1, 1, 0, 0};
      for (int d = 0; d < 3; d++)
      {
         worst_tri = std::max({worst_tri, std::abs(fr.n[d] - tri[d]),
                               std::abs(fr.d[d] - tri[3 + d]), std::abs(fr.s[d] - tri[6 + d])});
      }
      const bool tpv = fr.n[0] == 0 && fr.n[1] == -1 && fr.n[2] == 0
                       && fr.d[0] == 0 && fr.d[1] == 0 && fr.d[2] == -1
                       && fr.s[0] == 1 && fr.s[1] == 0 && fr.s[2] == 0;
      if (tpv) { n_tpv++; }
   }
   const int g_neg = ISum(n_eps_neg, comm), g_kflip = ISum(n_kernel_flip, comm);
   const int g_fb = ISum(n_fallback, comm), g_bad = ISum(n_frame_bad, comm);
   const int g_tpv = ISum(n_tpv, comm), g_nloc_int = ISum(n_int, comm);
   worst_cross = RMax(worst_cross, comm);
   worst_tri = RMax(worst_tri, comm);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6)
                << "  report: faces " << R.faces << ", components " << R.components
                << ", junction edges " << R.junction_edges << ", flipped_vs_raw "
                << R.flipped_vs_raw << ", min pin margin " << R.min_pin_margin
                << ", max d_F " << R.max_face_departure << " m, min side margin "
                << R.min_side_margin << " (x h_perp/4), parity per patch "
                << R.parity_per_patch_ok << ", global " << R.parity_global_ok << "\n"
                << "  (a) eps_f = -1 on " << g_neg << " interior faces; kernel stored normal "
                   "anti-aligned with (0,-1,0) on " << g_kflip << " (of " << g_nloc_int
                << " interior faces counted locally)\n"
                << "  (d) frames: max |d x s - n| = " << worst_cross << ", F6 fallbacks = "
                << g_fb << ", degenerate = " << g_bad << ", faces equal to the TPV102 triad "
                   "(0,-1,0)/(0,0,-1)/(1,0,0) exactly: " << g_tpv
                << "; max componentwise distance from that triad: " << worst_tri << "\n";
   }
   TEST_ASSERT(g_neg == g_kflip, "T-14(a) eps_f = -1 count == kernel flipped count");
   TEST_ASSERT(R.components == 1 && R.junction_edges == 0,
               "T-14(a) one component, zero junction edges");
   TEST_ASSERT(worst_cross == 0.0 && g_bad == 0, "T-14(d) d x s == n exactly on every face");
   TEST_ASSERT(g_fb == 0, "F7: zero F6 fallback faces");
   if (exact)
   {
      TEST_ASSERT(g_tpv == ISum(nloc, comm),
                  "N-09 TPV102 frame n=(0,-1,0), d=(0,0,-1), s=(1,0,0) exactly on every face");
   }
   else
   {
      TEST_ASSERT(worst_tri <= 1e-15,
                  "N-09 TPV102 frame to round-off (componentwise <= 1e-15) on every face");
      // off-plane magnitude of the fault-face vertices (information)
      real_t ymax = 0.0, yres = 0.0;
      const real_t tilt = 6.123233995736766e-17;   // Cos(Pi/2) in double
      for (int fi = 0; fi < nloc; fi++)
      {
         long long key[3]; real_t x[3][3];
         fo.FaceVertexData(fi, key, x);
         for (int a = 0; a < 3; a++)
         {
            ymax = std::max(ymax, std::abs(x[a][1]));
            yres = std::max(yres, std::abs(x[a][1] + x[a][2] * tilt));
         }
      }
      ymax = RMax(ymax, comm);
      yres = RMax(yres, comm);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  (info) fault-face vertices: max |y| = " << ymax
                   << " m, max |y + z * cos(pi/2)_double| = " << yres << " m\n";
      }
   }

   // (b) FromKernel(ToKernel(delta)) == delta on random tangential vectors.
   std::mt19937 rng(7 + g_rank);
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   real_t worst_rt = 0.0;
   for (int fi = 0; fi < nloc; fi++)
   {
      const FaultBasisData &b = fb->GetBasis(fi);
      for (int r = 0; r < 4; r++)
      {
         const real_t a1 = U(rng), a2 = U(rng);
         real_t delta[3], Tback[3], d1k, d2k;
         for (int d = 0; d < 3; d++) { delta[d] = a1 * frames[fi].d[d] + a2 * frames[fi].s[d]; }
         ff::ToKernel(delta, eps[fi], b.tangent1, b.tangent2, d1k, d2k);
         ff::FromKernel(d1k, d2k, 0.0, eps[fi], b.normal, b.tangent1, b.tangent2, Tback);
         for (int d = 0; d < 3; d++) { worst_rt = std::max(worst_rt, std::abs(Tback[d] - delta[d])); }
      }
   }
   worst_rt = RMax(worst_rt, comm);
   if (g_rank == 0) { std::cout << "  (b) max |FromKernel(ToKernel(delta)) - delta| = " << worst_rt << "\n"; }
   TEST_ASSERT(worst_rt <= 1e-15, "T-14(b) dictionary round trip == identity (1e-15)");

   // (c) kernel traction vs the dictionary path, bit for bit, every DOF.
   ParGridFunction u(&op.GetFESpace());
   u.Randomize(31 + g_rank);
   u -= 0.5;
   const int nd = op.GetNumFaultDOFs(), nbf = op.GetNbfPerFace();
   Vector slip(2 * nd);
   slip.Randomize(41 + g_rank);
   slip -= 0.5;
   Vector T, sn;
   op.ComputeTraction(u, slip, T, &sn);
   DenseMatrix B;
   op.GetFaultDOFBasis(B);
   int nbad = 0;
   real_t worst_c = 0.0, scale_c = 0.0;
   for (int k = 0; k < nd; k++)
   {
      const int fi = k / nbf;
      real_t nk[3], t1k[3], t2k[3], Tg[3];
      for (int d = 0; d < 3; d++) { nk[d] = B(d, k); t1k[d] = B(3 + d, k); t2k[d] = B(6 + d, k); }
      const real_t e = ff::EpsilonF(fo.Normal(fi), nk);
      ff::FromKernel(T(2 * k), T(2 * k + 1), sn(k), e, nk, t1k, t2k, Tg);
      const real_t Td = ff::Dot(Tg, frames[fi].d), Ts = ff::Dot(Tg, frames[fi].s);
      const real_t snd = -ff::Dot(Tg, frames[fi].n);
      if (!(Td == T(2 * k) && Ts == T(2 * k + 1) && snd == sn(k))) { nbad++; }
      worst_c = std::max({worst_c, std::abs(Td - T(2 * k)), std::abs(Ts - T(2 * k + 1)),
                          std::abs(snd - sn(k))});
      scale_c = std::max({scale_c, std::abs(T(2 * k)), std::abs(T(2 * k + 1)), std::abs(sn(k))});
   }
   const int g_nbad = ISum(nbad, comm), g_nd = ISum(nd, comm);
   worst_c = RMax(worst_c, comm);
   scale_c = RMax(scale_c, comm);
   if (g_rank == 0)
   {
      std::cout << "  (c) dictionary (T.d, T.s, sigma_n) vs kernel output: " << g_nbad
                << " of " << g_nd << " fault DOFs differ (bitwise); max |difference| = "
                << worst_c << " (traction scale " << scale_c << ", relative "
                << worst_c / scale_c << ")\n";
   }
   if (exact)
   {
      TEST_ASSERT(g_nbad == 0, "T-14(c) dictionary == kernel traction bit for bit (spec T1)");
   }
   else
   {
      TEST_ASSERT(worst_c <= 1e-14 * scale_c,
                  "T-14(c) dictionary == kernel traction to round-off (relative <= 1e-14; spec T1)");
      if (g_rank == 0)
      {
         std::cout << "  (info) bitwise mismatches: " << g_nbad << " of " << g_nd
                   << " DOFs; max componentwise triad distance " << worst_tri << "\n";
      }
   }

   // P2: Orient(n_raw) == Orient(-n_raw) exactly, on every face.
   int p2bad = 0;
   for (int fi = 0; fi < nloc; fi++)
   {
      const FaultBasisData &b = fb->GetBasis(fi);
      real_t nr[3] = {b.normal[0] * 3.7, b.normal[1] * 3.7, b.normal[2] * 3.7};
      real_t nm[3] = {-nr[0], -nr[1], -nr[2]};
      real_t o1[3], o2[3];
      fo.Orient(fi, nr, o1);
      fo.Orient(fi, nm, o2);
      if (std::memcmp(o1, o2, sizeof(o1)) != 0) { p2bad++; }
   }
   TEST_ASSERT(ISum(p2bad, comm) == 0, "P2 Orient(n_raw) == Orient(-n_raw) exactly");
}

// ------------------------------------------------- T11(a) and T12 (static)
static void TestT11T12(ParMesh &pm)
{
   if (g_rank == 0) { std::cout << "\n=== T11(a) tangent-basis invariance (constitutive level), T12 ===\n"; }
   MPI_Comm comm = pm.GetComm();
   LinearElastic model(TPV102Params::lambda, TPV102Params::mu);
   ElasticityDomainOperator<ParMesh> op(pm, 2, model, 0.0, kL, 2.0 * kL, MakeBC(),
                                        DGMethod::IP, SolverType::CG_AMG);
   FaultOrientation fo(pm, op.GetFaultInteriorFaces(), op.GetFaultSharedFaces());
   const real_t up[3] = {0, 0, 1};
   DieterichRuinaFriction fr({TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                              TPV102Params::Dc});
   RateStateClosure cl(fr, RateStateClosure::StateLaw::Aging);
   const real_t gamma = 7.67e6, sigma_n = 120e6;
   std::mt19937 rng(11 + g_rank);
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   real_t worst_inv = 0.0, worst_rot = 0.0;
   const real_t angles[3] = {0.3, M_PI / 2.0, -1.0};   // -1: per-face random
   for (int ia = 0; ia < 3; ia++)
   {
      for (int fi = 0; fi < fo.NumLocalFaces(); fi++)
      {
         ff::Frame F;
         FaceFrame(fo, fi, up, F);
         real_t th = angles[ia];
         if (th < -0.5)   // per-face angle seeded by the face's vertex key
         {
            const auto &k = fo.FaceKey(fi);
            std::mt19937_64 fr64(static_cast<unsigned long long>(k[0] * 1000003 + k[1] * 7919 + k[2]));
            th = std::uniform_real_distribution<real_t>(0.0, 2 * M_PI)(fr64);
         }
         const real_t c = std::cos(th), s = std::sin(th);
         real_t dQ[3], sQ[3];
         for (int d = 0; d < 3; d++) { dQ[d] = c * F.d[d] + s * F.s[d]; sQ[d] = -s * F.d[d] + c * F.s[d]; }
         // a global traction state (shear 70-80 MPa, random direction) + normal part
         real_t T[3];
         const real_t a1 = 75e6 * U(rng), a2 = 75e6 * U(rng);
         for (int d = 0; d < 3; d++) { T[d] = a1 * F.d[d] + a2 * F.s[d] - 1e6 * U(rng) * F.n[d]; }
         FaultPointState st; st.psi = 0.7; st.a = 0.008; st.Dc = 0.02;
         auto eval = [&](const real_t d_[3], const real_t s_[3], real_t X[2], real_t V[3],
                         real_t &snE, real_t &Vabs, real_t &D) {
            ff::ToLocal(T, sigma_n, F.n, d_, s_, X[0], X[1], snE);
            real_t Vl[2];
            Vabs = cl.SolveSlipRate(X, gamma, snE, st, Vl);
            ff::Embed(Vl[0], Vl[1], d_, s_, V);
            D = -(X[0] * Vl[0] + X[1] * Vl[1]);
         };
         real_t X0[2], V0[3], sn0, Va0, D0, X1[2], V1[3], sn1, Va1, D1;
         eval(F.d, F.s, X0, V0, sn0, Va0, D0);
         eval(dQ, sQ, X1, V1, sn1, Va1, D1);
         const real_t sc = std::max(Va0, 1e-300);
         for (int d = 0; d < 3; d++) { worst_inv = std::max(worst_inv, std::abs(V1[d] - V0[d]) / sc); }
         worst_inv = std::max(worst_inv, std::abs(sn1 - sn0) / std::abs(sn0));
         worst_inv = std::max(worst_inv, std::abs(Va1 - Va0) / sc);
         worst_inv = std::max(worst_inv, std::abs(D1 - D0) / std::max(std::abs(D0), 1e-300));
         // (tau1, tau2) rotated by Q^T
         const real_t Xr0 = c * X0[0] + s * X0[1], Xr1 = -s * X0[0] + c * X0[1];
         const real_t xs = std::hypot(X0[0], X0[1]);
         worst_rot = std::max({worst_rot, std::abs(X1[0] - Xr0) / xs, std::abs(X1[1] - Xr1) / xs});
      }
   }
   worst_inv = RMax(worst_inv, comm);
   worst_rot = RMax(worst_rot, comm);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6)
                << "  Q: global 0.3 rad, 90 deg, per-face random (seeded by vertex key)\n"
                << "  max relative change of V (global), sigma_n, |V|, D: " << worst_inv
                << "\n  max relative error of (tau1, tau2) vs Q^T-rotated: " << worst_rot << "\n";
   }
   TEST_ASSERT(worst_inv < 1e-13, "T11(a) V, sigma_n, |V|, D unchanged under (d, s) -> (d, s) Q");
   TEST_ASSERT(worst_rot < 1e-13, "T11(a) (tau1, tau2) rotated by Q^T");
   // T12 (static): the per-point friction parameters are scalars only.
   static_assert(sizeof(FaultPointState) == 3 * sizeof(real_t),
                 "FaultPointState must stay three scalars (psi, a, Dc): T12");
   TEST_ASSERT(true, "T12 static: FaultPointState is {psi, a, Dc}, no component index");
}

// --------------------------------------------------------------- T-16 dump
static void Dump(ParMesh &pm, const std::string &prefix)
{
   LinearElastic model(TPV102Params::lambda, TPV102Params::mu);
   ElasticityDomainOperator<ParMesh> op(pm, 1, model, 0.0, kL, 2.0 * kL, MakeBC(),
                                        DGMethod::IP, SolverType::CG_AMG);
   FaultOrientation fo(pm, op.GetFaultInteriorFaces(), op.GetFaultSharedFaces());
   const real_t up[3] = {0, 0, 1};
   std::ofstream os(prefix + "_r" + std::to_string(g_rank) + ".txt");
   os << std::setprecision(17);
   for (int fi = 0; fi < fo.NumLocalFaces(); fi++)
   {
      long long key[3]; real_t x[3][3];
      fo.FaceVertexData(fi, key, x);
      ff::Frame F;
      FaceFrame(fo, fi, up, F);
      real_t c[3];
      for (int d = 0; d < 3; d++) { c[d] = (x[0][d] + x[1][d] + x[2][d]) / 3.0; }
      // plus side as the centroid of the plus element (partition invariant)
      std::string plus = "remote";
      if (fo.PlusElement(fi) >= 0)
      {
         Array<int> ev;
         pm.GetElementVertices(fo.PlusElement(fi), ev);
         real_t ce[3] = {0, 0, 0};
         for (int a = 0; a < ev.Size(); a++)
            for (int d = 0; d < 3; d++) { ce[d] += pm.GetVertex(ev[a])[d] / ev.Size(); }
         plus = std::to_string(std::llround(ce[0] * 1000)) + "," + std::to_string(std::llround(ce[1] * 1000))
                + "," + std::to_string(std::llround(ce[2] * 1000));
      }
      os << c[0] << " " << c[1] << " " << c[2] << " "
         << F.n[0] << " " << F.n[1] << " " << F.n[2] << " "
         << F.d[0] << " " << F.d[1] << " " << F.d[2] << " "
         << F.s[0] << " " << F.s[1] << " " << F.s[2] << " "
         << fo.Patch(fi) << " " << fo.PersistentId(fi) << " " << plus << " "
         << key[0] << "," << key[1] << "," << key[2] << "\n";
   }
   if (g_rank == 0)
   {
      const auto &R = fo.GetReport();
      std::cout << "  dump: faces " << R.faces << ", components " << R.components
                << ", flipped_vs_raw " << R.flipped_vs_raw << "\n";
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const std::string mode = (argc > 1) ? argv[1] : "";
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche PF (fault coordinate layer) on " << Mpi::WorldSize()
                << " rank(s) " << mode << " ===\n";
   }
   if (mode == "--gmsh" && argc > 2)
   {
      Mesh m(argv[2], 1, 1);
      int non_tet = 0;
      for (int e = 0; e < m.GetNE(); e++)
      {
         if (m.GetElementBaseGeometry(e) != Geometry::TETRAHEDRON) { non_tet++; }
      }
      if (g_rank == 0)
      {
         std::cout << "  mesh " << argv[2] << ": " << m.GetNE() << " elements, "
                   << non_tet << " non-tetrahedral, bdr attrs";
         for (int a : m.bdr_attributes) { std::cout << " " << a; }
         std::cout << "\n";
      }
      TEST_ASSERT(non_tet == 0, "Gmsh mesh is all-tetrahedral");
      ParMesh pm(MPI_COMM_WORLD, m);
      TestT14(pm, 1, "tpv102_1000m.msh", false);
   }
   else if (mode == "--dump" && argc > 2)
   {
      std::unique_ptr<Mesh> m(FaultBox(4, 4, 2));
      ParMesh pm(MPI_COMM_WORLD, *m);
      Dump(pm, argv[2]);
   }
   else
   {
      std::unique_ptr<Mesh> m(FaultBox(4, 4, 2));
      ParMesh pm(MPI_COMM_WORLD, *m);
      TestT14(pm, 2, "inline 24-tet 4x4x2 box", true);
      TestT11T12(pm);
   }
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                << " passed, " << num_failed << " failed ===\n";
   }
   return (num_failed == 0) ? 0 : 1;
}
