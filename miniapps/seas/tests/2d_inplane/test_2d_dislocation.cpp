// tests/2d_inplane/test_2d_dislocation.cpp — 2D in-plane design rev 3,
// phase E (section 7, row E): the plane-strain static dislocation check, the
// 2D instance of phase A's Okada test (test_okada_dislocation.cpp).  A
// uniform dip-slip patch on a dipping fault in a half-space, imposed through
// the QD operator's slip embedding (ElasticityDomainOperator in field mode,
// fault/fault_frame_field.hpp) on a 2D triangle mesh, against the analytic
// plane-strain solution: the surface displacement and the fault tractions
// (dip shear and normal), for dips 30, 60 and 90 degrees, thrust and normal
// sense, p = 1 and 2, at two resolutions.
//
// The analytic solution is Okada (1992) (okada92.hpp, the transcription of
// his DC3D the phase A test checks) for a rectangle of length 2 L = 2000 km
// along strike, evaluated in its middle plane x_o = 0: plane strain to the
// precision the tests below need (P0 checks that L = 1000 km and 10000 km
// agree and that the strike components vanish).  Okada's frame (x_o along
// strike, y_o horizontal, z_o up) maps to the mesh frame (X, Y, Z) by the
// cyclic permutation X = y_o, Y = z_o, Z = x_o (a proper rotation): the
// mesh lies in the (X, Y) plane, Y up, the medium Y <= 0, the free surface
// Y = 0 (BP3's layout), and the fault's strike is out of the plane.
//
// Problem: lambda = mu = 32.04 GPa (alpha = 2/3); the patch -1 <= eta <= 1
// km along dip (W = 2 km) about depth 2 km, slip 1 m: Okada's dip slip U2 =
// +1 (thrust, the hanging wall up) or -1 (normal), b = U2 (cos, sin) in
// (X, Y).  Mesh: a grid in (x', z), x' in [-1, 1] km, z from 0.5 km below
// the patch to 0, each cell split into 4 triangles about its centre, sheared
// X = x' + (z + 2 km) cot(dip) so that the grid line x' = 0 is the fault
// line; n cells along the patch's dip and cells of the same height
// elsewhere; the patch's edges are the fault (attr 3), the top Y = 0 is free
// (attr 1), the other sides carry Okada's displacement as Dirichlet data
// (attr 5): the difference to Okada is the discretisation's alone.  up =
// (0, 1) (BP3's), ref_normal = the normal into the foot wall, nu_FW = (sin,
// -cos) in (X, Y), unless said otherwise.  With up in the plane the frame
// puts the slip in the dip slot: s = up x n is out of the plane (design rev
// 3, rule R2, derivation 10.2), so the strike slot of the slip and of the
// traction is an exact zero.
//
// Conventions (phase A's): the operator's slip is u+ - u- = S_dip d with the
// plus side the one n points away from, so S = b . d when n points into the
// foot wall and -b . d otherwise; its traction (ComputeTraction) is the
// traction on the plus side, sigma n, in (d, s), and the normal traction
// -n . sigma n (compression positive).
//
//   P0  the plane-strain limit, no mesh (dips 30, 60, 90): the strike
//       components of u and of sigma . nu vanish in x_o = 0 (1e-9 of |b|
//       and of mu |b| / W); L = 1000 and 10000 km agree (5e-6 of the same)
//       up to a rigid offset of the box (vertical, ~ W sin(2 dip) / L);
//       u(hanging wall) - u(foot wall) = b on the patch and 0 off it (1e-3
//       m); sigma . e_Y = 0 on Y = 0 (1e-10 of mu |b| / W)
//   P1  dips 30, 60, 90, thrust, ref_normal = nu_FW: for p = 1 and 2, n =
//       4 and 8: e_u, the relative L2 error of the surface displacement (top
//       edge, order 2p + 2); e_T, the relative RMS error of (tau_dip, normal
//       traction) over the DOFs of the fault faces at least 0.5 km from the
//       patch edges, and the same for the traction's stress part {sigma} n
//       alone.  Pass: e_u < 0.25 and the stress part < 0.5 (a reversed sign
//       gives about 2) at every n; from n = 4 to 8 e_u and the stress part
//       fall and e_T falls by more than 2.  The strike slot of the slip and
//       of the traction is exactly 0 at every DOF.
//   P2  normal sense (U2 = -1), n = 4, p = 1, each dip: the tractions are
//       the thrust case's negated DOF by DOF (1e-9 of the largest)
//   P3  ref_normal = -nu_FW (n flipped, so S = -b . d), dip 60: the physical
//       driving traction tau_dip d reverses and the normal traction stays,
//       DOF by DOF (1e-9 of the largest)
//   P4  rotation in the plane: mesh, Okada's field, ref_normal and up
//       rotated by 30 degrees about Z, dip 60: (tau_dip, tau_strike, normal
//       traction) DOF by DOF the same as unrotated (1e-8 of the largest)
//
// Modes: (default) P0-P4;  --quick: n = 4 only (no fall checked);
// --single DIP P N: one P1 case, printed (no check).
//
//   make -f tests/2d_inplane/2d_inplane.mk test-2d-dislocation
//
// Scope (Chunhui Zhao, 2026-10-26).

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"
#include "frame_box.hpp"
#include "okada92.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using frame_box::Rot;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kLambda = 32.04e9, kMu = 32.04e9;
static const real_t kAlpha = (kLambda + kMu) / (kLambda + 2.0 * kMu);
// patch (km, Okada's DC3D arguments); kL: half the length along strike
static const real_t kDepth = 2.0, kAW1 = -1.0, kAW2 = 1.0, kL = 1000.0;
static const real_t kW = (kAW2 - kAW1) * 1e3;      ///< patch width [m]
static const real_t kSigScale = kMu * 1.0 / kW;    ///< mu |b| / W, |b| = 1 m
// box (km): x' in [-kMarginX, kMarginX], from kMargin below the patch to z = 0
static const real_t kMargin = 0.5, kMarginX = 1.0;
static const real_t kEdgeGap = 0.5;   ///< traction faces at least this far from the patch edges [km]

static real_t RSum(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }

/// Okada's sin and cos of the dip (DCCON0: cos set to 0 below 1e-6).
static void SinCos(real_t dip, real_t &sd, real_t &cd)
{
   sd = std::sin(dip * M_PI / 180.0);
   cd = std::cos(dip * M_PI / 180.0);
   if (std::abs(cd) < 1e-6) { cd = 0.0; sd = (sd > 0.0) ? 1.0 : -1.0; }
}

/// Okada's frame index of each mesh axis: X = y_o, Y = z_o, Z = x_o.
static const int kPerm[3] = {1, 2, 0};

// ------------------------------------------------------------ Okada's field
/// The slip vector b (hanging wall minus foot wall) in (X, Y, Z).
static void SlipVector(real_t U2, real_t dip, real_t b[3])
{
   real_t sd, cd;
   SinCos(dip, sd, cd);
   b[0] = U2 * cd; b[1] = U2 * sd; b[2] = 0.0;
}

/// The normal into the foot wall in (X, Y, Z): Okada's (0, sin, -cos).
static void FootWallNormal(real_t dip, real_t nu[3])
{
   real_t sd, cd;
   SinCos(dip, sd, cd);
   nu[0] = sd; nu[1] = -cd; nu[2] = 0.0;
}

/// Okada's displacement [m] and gradient [1] at X [m] (mesh frame, rotated
/// by R: X = R X_m), for a patch of half-length L [km] along strike, in its
/// middle plane x_o = Z_m (0 for a mesh point): u = R P u_o, G = R P G_o
/// P^T R^T, P the permutation kPerm.
static int OkadaAt(real_t U2, real_t dip, const Rot &rot, const real_t X[3],
                   real_t u[3], real_t G[3][3], real_t L = kL)
{
   real_t Xm[3];
   for (int i = 0; i < 3; i++) { Xm[i] = rot.R[0][i] * X[0] + rot.R[1][i] * X[1] + rot.R[2][i] * X[2]; }
   const real_t xo = Xm[2], yo = Xm[0], zo = std::min(Xm[1], 0.0);
   const okada92::Result r = okada92::DC3D(kAlpha, xo * 1e-3, yo * 1e-3, zo * 1e-3, kDepth,
                                           dip, -L, L, kAW1, kAW2, 0.0, U2, 0.0);
   real_t um[3], Gm[3][3];
   for (int i = 0; i < 3; i++)
   {
      um[i] = r.u[kPerm[i]];
      for (int j = 0; j < 3; j++) { Gm[i][j] = r.du[kPerm[j]][kPerm[i]] * 1e-3; }   // per km -> per m
   }
   for (int i = 0; i < 3; i++)
   {
      u[i] = rot.R[i][0] * um[0] + rot.R[i][1] * um[1] + rot.R[i][2] * um[2];
      for (int j = 0; j < 3; j++)
      {
         real_t s = 0.0;
         for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++) { s += rot.R[i][k] * Gm[k][l] * rot.R[j][l]; }
         G[i][j] = s;
      }
   }
   return r.iret;
}

static void Stress(const real_t G[3][3], real_t S[3][3])
{
   const real_t tr = G[0][0] + G[1][1] + G[2][2];
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         S[i][j] = kMu * (G[i][j] + G[j][i]) + (i == j ? kLambda * tr : 0.0);
      }
}

// ---------------------------------------------------------------- P0
static void TestPlaneStrainLimit()
{
   if (g_rank == 0) { std::cout << "\n=== P0: Okada's plane-strain limit (no mesh) ===\n"; }
   const real_t eps = 1e-5;   // km off the plane
   for (real_t dip : {30.0, 60.0, 90.0})
   {
      real_t sd, cd, b[3], nu[3];
      SinCos(dip, sd, cd);
      SlipVector(1.0, dip, b);
      FootWallNormal(dip, nu);
      real_t strike_u = 0.0, strike_t = 0.0, dL_u = 0.0, dL_s = 0.0, jump_on = 0.0, jump_off = 0.0;
      for (int k = 0; k <= 16; k++)
      {
         const real_t eta = -1.6 + 3.2 * k / 16.0 + 0.0071;
         const bool on = (eta > kAW1 && eta < kAW2);
         const bool near_edge = std::abs(eta - kAW1) < 0.1 || std::abs(eta - kAW2) < 0.1;
         const real_t P[3] = {eta * cd, -kDepth + eta * sd, 0.0};
         if (P[1] > -0.05 || near_edge) { continue; }
         real_t Xh[3], Xf[3];
         for (int c = 0; c < 3; c++) { Xh[c] = (P[c] - eps * nu[c]) * 1e3; Xf[c] = (P[c] + eps * nu[c]) * 1e3; }
         real_t uh[3], uf[3], Gh[3][3], Gf[3][3];
         OkadaAt(1.0, dip, Rot(), Xh, uh, Gh);
         OkadaAt(1.0, dip, Rot(), Xf, uf, Gf);
         real_t e = 0.0;
         for (int c = 0; c < 3; c++) { e = std::max(e, std::abs(uh[c] - uf[c] - (on ? b[c] : 0.0))); }
         (on ? jump_on : jump_off) = std::max(on ? jump_on : jump_off, e);
      }
      // the strike components and the L dependence on a grid off the plane
      std::vector<std::array<real_t, 3>> du;
      for (int i = 0; i <= 10; i++)
         for (int j = 0; j <= 10; j++)
         {
            const real_t X[3] = {(-5.0 + i + 0.031) * 1e3, (-0.2 - 0.5 * j - 0.017) * 1e3, 0.0};
            real_t u1[3], u2[3], G1[3][3], G2[3][3], S1[3][3], S2[3][3];
            OkadaAt(1.0, dip, Rot(), X, u1, G1, 1000.0);
            OkadaAt(1.0, dip, Rot(), X, u2, G2, 10000.0);
            Stress(G1, S1); Stress(G2, S2);
            strike_u = std::max(strike_u, std::abs(u1[2]));
            for (int r = 0; r < 2; r++) { strike_t = std::max(strike_t, std::abs(S1[2][r]) / kSigScale); }
            du.push_back({u1[0] - u2[0], u1[1] - u2[1], u1[2] - u2[2]});
            for (int r = 0; r < 3; r++)
               for (int c = 0; c < 3; c++) { dL_s = std::max(dL_s, std::abs(S1[r][c] - S2[r][c]) / kSigScale); }
         }
      // the finite length moves the whole box by a rigid offset (vertical,
      // about W sin(2 dip) / L: 1.24e-4 m at L = 1000 km against 10000 km
      // for dips 30 and 60, 0 at 90) that strains nothing; what is left
      // after it is the plane-strain error
      real_t mean[3] = {0, 0, 0}, off = 0.0;
      for (const auto &d : du) { for (int c = 0; c < 3; c++) { mean[c] += d[c] / du.size(); } }
      for (const auto &d : du) { for (int c = 0; c < 3; c++) { dL_u = std::max(dL_u, std::abs(d[c] - mean[c])); } }
      for (int c = 0; c < 3; c++) { off = std::max(off, std::abs(mean[c])); }
      real_t free = 0.0;
      for (int i = 0; i <= 24; i++)
      {
         const real_t X[3] = {(-6.0 + 0.5 * i + 0.031) * 1e3, 0.0, 0.0};
         real_t u[3], G[3][3], S[3][3];
         OkadaAt(1.0, dip, Rot(), X, u, G);
         Stress(G, S);
         for (int r = 0; r < 3; r++) { free = std::max(free, std::abs(S[r][1]) / kSigScale); }
      }
      std::ostringstream tag;
      tag << "P0 dip " << dip << ": ";
      TEST_ASSERT(strike_u <= 1e-9 && strike_t <= 1e-9, tag.str() << "strike components vanish in x_o = 0 (u "
                  << strike_u << " m, sigma " << strike_t << " of mu |b| / W)");
      TEST_ASSERT(dL_u <= 5e-6 && dL_s <= 5e-6, tag.str() << "L = 1000 km and 10000 km agree up "
                  "to a rigid offset (" << off << " m): u " << dL_u << " m, sigma " << dL_s
                  << " of mu |b| / W");
      TEST_ASSERT(jump_on <= 1e-3 && jump_off <= 1e-3, tag.str() << "u(hanging wall) - u(foot wall) = b "
                  "on the patch (" << jump_on << " m), 0 off it (" << jump_off << " m)");
      TEST_ASSERT(free <= 1e-10, tag.str() << "sigma . e_Y = 0 on Y = 0 (" << free << " of mu |b| / W)");
   }
}

// ------------------------------------------------------------- the FE side
/// The sheared grid (see the file comment), rotated by R after tagging.
static Mesh *DislocationBox(int n, real_t dip, const Rot &rot)
{
   real_t sd, cd;
   SinCos(dip, sd, cd);
   const real_t hz = (kAW2 - kAW1) * sd / n;
   const real_t zt = -kDepth + kAW2 * sd, zb = -kDepth + kAW1 * sd;
   auto Nodes = [](std::vector<std::pair<real_t, real_t>> segs, real_t h, std::vector<int> fixed)
   {
      std::vector<real_t> v;
      for (size_t s = 0; s < segs.size(); s++)
      {
         const real_t a = segs[s].first, b = segs[s].second;
         const int m = fixed[s] > 0 ? fixed[s] : std::max(1, (int)std::ceil((b - a) / h - 1e-9));
         for (int i = (s == 0 ? 0 : 1); i <= m; i++) { v.push_back(a + (b - a) * i / m); }
      }
      return v;
   };
   const std::vector<real_t> xs = Nodes({{-kMarginX, 0.0}, {0.0, kMarginX}}, hz, {0, 0});
   const std::vector<real_t> zs = Nodes({{zb - kMargin, zb}, {zb, zt}, {zt, 0.0}}, hz, {0, n, 0});
   const int nx = (int)xs.size() - 1, nz = (int)zs.size() - 1;
   int i0 = -1, kb = -1, kt = -1;
   for (int i = 0; i <= nx; i++) { if (xs[i] == 0.0) { i0 = i; } }
   for (int k = 0; k <= nz; k++) { if (zs[k] == zb) { kb = k; } if (zs[k] == zt) { kt = k; } }
   MFEM_VERIFY(i0 > 0 && kb >= 0 && kt > kb, "DislocationBox: patch grid lines");

   const real_t cot = (sd != 0.0) ? cd / sd : 0.0;
   const int ngv = (nx + 1) * (nz + 1);
   Mesh *m = new Mesh(2, ngv + nx * nz, 4 * nx * nz, 2 * (nx + nz) + (kt - kb), 2);
   auto gv = [&](int i, int k) { return i + (nx + 1) * k; };
   auto cv = [&](int i, int k) { return ngv + i + nx * k; };
   auto put = [&](real_t xp, real_t z)
   {
      const real_t X[3] = {(xp + (z + kDepth) * cot) * 1e3, z * 1e3, 0.0};
      real_t Y[3];
      rot.Apply(X, Y);
      m->AddVertex(Y[0], Y[1]);
   };
   for (int k = 0; k <= nz; k++)
      for (int i = 0; i <= nx; i++) { put(xs[i], zs[k]); }
   for (int k = 0; k < nz; k++)
      for (int i = 0; i < nx; i++) { put(0.5 * (xs[i] + xs[i + 1]), 0.5 * (zs[k] + zs[k + 1])); }
   for (int k = 0; k < nz; k++)
      for (int i = 0; i < nx; i++)
      {
         const int c = cv(i, k);
         m->AddTriangle(gv(i, k), gv(i + 1, k), c, 1);
         m->AddTriangle(gv(i + 1, k), gv(i + 1, k + 1), c, 1);
         m->AddTriangle(gv(i + 1, k + 1), gv(i, k + 1), c, 1);
         m->AddTriangle(gv(i, k + 1), gv(i, k), c, 1);
      }
   for (int i = 0; i < nx; i++)
   {
      m->AddBdrSegment(gv(i, 0), gv(i + 1, 0), 5);
      m->AddBdrSegment(gv(i + 1, nz), gv(i, nz), 1);
   }
   for (int k = 0; k < nz; k++)
   {
      m->AddBdrSegment(gv(nx, k), gv(nx, k + 1), 5);
      m->AddBdrSegment(gv(0, k + 1), gv(0, k), 5);
   }
   for (int k = kb; k < kt; k++) { m->AddBdrSegment(gv(i0, k), gv(i0, k + 1), 3); }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

struct Case
{
   real_t dip = 60.0;
   real_t U2 = 1.0;              ///< Okada's dip slip: +1 thrust, -1 normal
   int p = 1, n = 4;
   real_t ref_sign = 1.0;        ///< ref_normal = ref_sign nu_FW
   Rot rot;
};

struct Result
{
   real_t e_u = 0.0, e_T = 0.0;
   real_t e_Ts = 0.0;   ///< the same for the stress part {sigma} n alone (diagnostic)
   real_t jump = 0.0;   ///< largest |[u] - S| at the region's DOFs [m] (diagnostic)
   real_t strike = 0.0; ///< largest |slip or traction in the strike slot| over all DOFs
   int n_dofs = 0, n_trac_dofs = 0, n_fault_faces = 0;
   std::vector<real_t> trac;    ///< owned fault DOFs: tau_dip, tau_strike, normal traction
   std::vector<real_t> tphys;   ///< owned fault DOFs: tau_dip d + tau_strike s (3), normal traction
};

using Op = ElasticityDomainOperator<ParMesh>;

static Result Run(MPI_Comm comm, const Case &cs)
{
   std::unique_ptr<Mesh> m(DislocationBox(cs.n, cs.dip, cs.rot));
   ParMesh pm(comm, *m);
   m.reset();
   real_t nu0[3], nu[3], up0[3] = {0, 1, 0}, up[3], ref[3];
   FootWallNormal(cs.dip, nu0);
   cs.rot.Apply(nu0, nu);
   cs.rot.Apply(up0, up);
   for (int d = 0; d < 3; d++) { ref[d] = cs.ref_sign * nu[d]; }

   const real_t U2 = cs.U2, dip = cs.dip;
   const Rot rot = cs.rot;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [U2, dip, rot](const Vector &x, real_t, Vector &u)
   {
      u.SetSize(2);
      const real_t X[3] = {x(0), x(1), 0.0};
      real_t uu[3], G[3][3];
      OkadaAt(U2, dip, rot, X, uu, G);
      u(0) = uu[0]; u(1) = uu[1];
   };
   auto fc = std::make_shared<FaultFrameConfig>();
   for (int d = 0; d < 3; d++)
   {
      fc->orientation.ref_normal[d] = ref[d];
      fc->orientation.up[d] = up[d];
      fc->up[d] = up[d];
   }
   fc->orientation.fault_attrs = {3};
   DomainConfig dc;
   dc.fault_frames.config = fc;
   LinearElastic model(kLambda, kMu);
   Op op(pm, cs.p, model, 0.0, 1e4, 1e4, bc, DGMethod::IP, SolverType::MUMPS, dc);
   auto F = op.FrameField();
   MFEM_VERIFY(F, "dislocation: no frame field");

   // slip: u+ - u- = J, J = b when n points into the foot wall, else -b
   real_t b0[3], b[3];
   SlipVector(U2, dip, b0);
   rot.Apply(b0, b);
   const int nd = op.GetNumFaultDOFs();
   const int nf = F->NumFaultFaces();
   MFEM_VERIFY(nf == 0 || nd % nf == 0, "dislocation: DOFs per face");
   const int nbf = nf ? nd / nf : 1;
   Vector slip(2 * nd);
   real_t strike = 0.0;
   for (int k = 0; k < nd; k++)
   {
      const fault_frame::Frame &fr = F->Frame(k / nbf);
      const real_t side = (fr.n[0] * nu[0] + fr.n[1] * nu[1] + fr.n[2] * nu[2] > 0.0) ? 1.0 : -1.0;
      slip(2 * k) = side * (b[0] * fr.d[0] + b[1] * fr.d[1] + b[2] * fr.d[2]);
      slip(2 * k + 1) = side * (b[0] * fr.s[0] + b[1] * fr.s[1] + b[2] * fr.s[2]);
      strike = std::max(strike, std::abs(slip(2 * k + 1)));
   }
   ParGridFunction u(&op.GetFESpace());
   op.Solve(0.0, slip, u);
   Vector trac, sn, xyz, trac_s, trac_c, jres, sn_s, sn_c, sn_d;
   op.ComputeTraction(u, slip, trac, &sn);   // the QD regime's traction
   Vector trac_d;
   op.ComputeTractionDiagnostics(u, slip, trac_d, trac_s, trac_c, jres, &sn_d, &sn_s, &sn_c);
   op.GetFaultDOFCoords3D(xyz);

   Result res;
   res.n_dofs = (int)RSum((real_t)op.GetFESpace().GetTrueVSize(), comm);
   res.n_fault_faces = (int)RSum((real_t)op.GetNumOwnedFaultFaces(), comm);

   // fault tractions against Okada's on the faces whose DOFs all lie at
   // least kEdgeGap from the patch edges
   real_t sd, cd;
   SinCos(dip, sd, cd);
   auto Eta = [&](int k)
   {
      const real_t X[3] = {xyz(3 * k), xyz(3 * k + 1), xyz(3 * k + 2)};
      real_t Xm[3];
      for (int i = 0; i < 3; i++) { Xm[i] = (rot.R[0][i] * X[0] + rot.R[1][i] * X[1] + rot.R[2][i] * X[2]) * 1e-3; }
      return Xm[0] * cd + (Xm[1] + kDepth) * sd;
   };
   Vector err2(nd), errs2(nd), ref2(nd), tl(3 * nd), tp(4 * nd), cnt(nd), jmax(nd);
   err2 = 0.0; errs2 = 0.0; ref2 = 0.0; cnt = 0.0; jmax = 0.0;
   for (int f = 0; f < (nbf ? nd / nbf : 0); f++)
   {
      bool in = true;
      for (int kk = 0; kk < nbf; kk++)
      {
         const real_t eta = Eta(f * nbf + kk);
         if (eta < kAW1 + kEdgeGap - 1e-9 || eta > kAW2 - kEdgeGap + 1e-9) { in = false; }
      }
      for (int kk = 0; kk < nbf; kk++)
      {
         const int k = f * nbf + kk;
         const fault_frame::Frame &fr = F->Frame(f);
         const real_t th[3] = {trac(2 * k), trac(2 * k + 1), sn(k)};
         const real_t ts[3] = {trac_s(2 * k), trac_s(2 * k + 1), sn_s(k)};
         strike = std::max(strike, std::abs(th[1]));
         for (int c = 0; c < 3; c++) { tl(3 * k + c) = th[c]; tp(4 * k + c) = th[0] * fr.d[c] + th[1] * fr.s[c]; }
         tp(4 * k + 3) = th[2];
         if (!in) { continue; }
         const real_t X[3] = {xyz(3 * k), xyz(3 * k + 1), xyz(3 * k + 2)};
         real_t uu[3], G[3][3], S[3][3], t[3];
         OkadaAt(U2, dip, rot, X, uu, G);
         Stress(G, S);
         for (int r = 0; r < 3; r++) { t[r] = S[r][0] * fr.n[0] + S[r][1] * fr.n[1] + S[r][2] * fr.n[2]; }
         // the traction on the plus side, sigma n . (d, s) (n points from the
         // plus into the minus side); normal: -n . sigma n (compression > 0)
         const real_t tok[3] = {t[0] * fr.d[0] + t[1] * fr.d[1] + t[2] * fr.d[2],
                                t[0] * fr.s[0] + t[1] * fr.s[1] + t[2] * fr.s[2],
                                -(t[0] * fr.n[0] + t[1] * fr.n[1] + t[2] * fr.n[2])};
         for (int c = 0; c < 3; c++)
         {
            err2(k) += (th[c] - tok[c]) * (th[c] - tok[c]);
            errs2(k) += (ts[c] - tok[c]) * (ts[c] - tok[c]);
            ref2(k) += tok[c] * tok[c];
         }
         jmax(k) = std::max(std::abs(jres(2 * k)), std::abs(jres(2 * k + 1)));
         cnt(k) = 1.0;
      }
   }
   Vector err2_o, errs2_o, ref2_o, tl_o, tp_o, cnt_o, jmax_o;
   op.RestrictToOwnedFault(err2, err2_o, 1);
   op.RestrictToOwnedFault(errs2, errs2_o, 1);
   op.RestrictToOwnedFault(ref2, ref2_o, 1);
   op.RestrictToOwnedFault(cnt, cnt_o, 1);
   op.RestrictToOwnedFault(jmax, jmax_o, 1);
   op.RestrictToOwnedFault(tl, tl_o, 3);
   op.RestrictToOwnedFault(tp, tp_o, 4);
   const real_t ref_sum = RSum(ref2_o.Sum(), comm);
   res.e_Ts = (ref_sum > 0.0) ? std::sqrt(RSum(errs2_o.Sum(), comm) / ref_sum) : 0.0;
   res.jump = RMax(jmax_o.Size() ? jmax_o.Max() : 0.0, comm);
   res.e_T = (ref_sum > 0.0) ? std::sqrt(RSum(err2_o.Sum(), comm) / ref_sum) : 1e300;
   res.n_trac_dofs = (int)RSum(cnt_o.Sum(), comm);
   res.strike = RMax(strike, comm);
   res.trac.assign(tl_o.GetData(), tl_o.GetData() + tl_o.Size());
   res.tphys.assign(tp_o.GetData(), tp_o.GetData() + tp_o.Size());

   // surface displacement against Okada's, L2 over the top edge
   real_t e2 = 0.0, r2 = 0.0;
   Vector val(2), xq(2);
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) != 1) { continue; }
      FaceElementTransformations *FTr = pm.GetBdrFaceTransformations(be);
      if (!FTr) { continue; }
      const IntegrationRule &ir = IntRules.Get(FTr->GetGeometryType(), 2 * cs.p + 2);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         FTr->SetAllIntPoints(&ip);
         FTr->Transform(ip, xq);
         u.GetVectorValue(*FTr->Elem1, FTr->GetElement1IntPoint(), val);
         const real_t X[3] = {xq(0), xq(1), 0.0};
         real_t uu[3], G[3][3];
         OkadaAt(U2, dip, rot, X, uu, G);
         const real_t w = ip.weight * FTr->Weight();
         for (int c = 0; c < 2; c++) { e2 += w * (val(c) - uu[c]) * (val(c) - uu[c]); r2 += w * uu[c] * uu[c]; }
      }
   }
   res.e_u = std::sqrt(RSum(e2, comm) / RSum(r2, comm));
   return res;
}

static void Print(const Case &cs, const Result &r)
{
   if (g_rank != 0) { return; }
   std::cout << "    dip " << std::setw(2) << cs.dip << "  " << (cs.U2 > 0 ? "thrust" : "normal")
             << "  p " << cs.p << "  n " << cs.n << "  DOFs " << std::setw(6) << r.n_dofs
             << "  fault faces " << std::setw(3) << r.n_fault_faces << "  traction DOFs " << std::setw(3)
             << r.n_trac_dofs << "  e_u " << std::scientific << std::setprecision(3) << r.e_u
             << "  e_T " << r.e_T << "  (stress part " << r.e_Ts << ", |[u] - S| <= " << r.jump
             << " m, strike slot " << r.strike << ")" << std::defaultfloat << "\n";
}

static real_t MaxAbs(const std::vector<real_t> &v, MPI_Comm comm)
{
   real_t m = 0.0;
   for (real_t x : v) { m = std::max(m, std::abs(x)); }
   return RMax(m, comm);
}

// ---------------------------------------------------------------- P1, P2
static void TestConvergence(MPI_Comm comm, bool quick)
{
   if (g_rank == 0) { std::cout << "\n=== P1, P2: convergence to the plane-strain solution; the normal sense ===\n"; }
   std::vector<int> ns = quick ? std::vector<int>{4} : std::vector<int>{4, 8};
   for (real_t dip : {30.0, 60.0, 90.0})
   {
      for (int p : {1, 2})
      {
         std::vector<Result> rs;
         for (int n : ns)
         {
            Case cs;
            cs.dip = dip;
            cs.p = p;
            cs.n = n;
            rs.push_back(Run(comm, cs));
            Print(cs, rs.back());
         }
         std::ostringstream tag;
         tag << "P1 dip " << dip << ", p = " << p << ": ";
         for (size_t i = 0; i < rs.size(); i++)
         {
            TEST_ASSERT(rs[i].e_u < 0.25 && rs[i].e_Ts < 0.5,
                        tag.str() << "n = " << ns[i] << ": e_u " << rs[i].e_u << " (< 0.25), stress "
                        "part " << rs[i].e_Ts << " (< 0.5; a reversed sign gives about 2)");
            TEST_ASSERT(rs[i].strike == 0.0, tag.str() << "n = " << ns[i] << ": the strike slot of the "
                        "slip and of the traction is exactly 0 (" << rs[i].strike << ")");
         }
         if (rs.size() == 2)
         {
            TEST_ASSERT(rs[1].e_u < rs[0].e_u && rs[1].e_Ts < rs[0].e_Ts && rs[1].e_T < 0.5 * rs[0].e_T,
                        tag.str() << "n = 4 -> 8: e_u / " << rs[0].e_u / rs[1].e_u << ", stress part / "
                        << rs[0].e_Ts / rs[1].e_Ts << ", e_T / " << rs[0].e_T / rs[1].e_T << " (> 2)");
         }
      }
      // P2: the normal sense negates the tractions
      Case th, nm;
      th.dip = nm.dip = dip;
      nm.U2 = -1.0;
      const Result r0 = Run(comm, th), r1 = Run(comm, nm);
      Print(nm, r1);
      real_t e = 0.0;
      for (size_t k = 0; k < r0.trac.size() && k < r1.trac.size(); k++) { e = std::max(e, std::abs(r1.trac[k] + r0.trac[k])); }
      e = RMax(e, comm) / MaxAbs(r0.trac, comm);
      TEST_ASSERT(r1.trac.size() == r0.trac.size() && e <= 1e-9 && r1.e_u < 0.25,
                  "P2 dip " << dip << ": normal sense: tractions negated (" << e << " of the largest), e_u "
                  << r1.e_u);
   }
}

// ---------------------------------------------------------------- P3, P4
static void TestSignAndRotation(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== P3, P4: the normal reversed; the problem rotated in the plane ===\n"; }
   Case base;
   const Result r0 = Run(comm, base);
   Print(base, r0);

   Case flip = base;
   flip.ref_sign = -1.0;
   const Result r1 = Run(comm, flip);
   real_t e = 0.0;
   for (size_t k = 0; k < r0.tphys.size() && k < r1.tphys.size(); k += 4)
   {
      for (int c = 0; c < 3; c++) { e = std::max(e, std::abs(r1.tphys[k + c] + r0.tphys[k + c])); }
      e = std::max(e, std::abs(r1.tphys[k + 3] - r0.tphys[k + 3]));
   }
   e = RMax(e, comm) / MaxAbs(r0.tphys, comm);
   TEST_ASSERT(r1.tphys.size() == r0.tphys.size() && e <= 1e-9,
               "P3 ref_normal -nu_FW: driving traction reversed, normal traction the same (" << e
               << " of the largest)");

   Case rotc = base;
   rotc.rot = Rot::Axis(0.0, 0.0, 1.0, 30.0);
   const Result r2 = Run(comm, rotc);
   real_t er = 0.0;
   for (size_t k = 0; k < r0.trac.size() && k < r2.trac.size(); k++) { er = std::max(er, std::abs(r2.trac[k] - r0.trac[k])); }
   er = RMax(er, comm) / MaxAbs(r0.trac, comm);
   TEST_ASSERT(r2.trac.size() == r0.trac.size() && er <= 1e-8,
               "P4 rotated by 30 degrees about Z: (tau_dip, tau_strike, normal traction) the same ("
               << er << " of the largest; e_u " << r2.e_u << " vs " << r0.e_u << ")");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   g_rank = Mpi::WorldRank();
   bool quick = false;
   int single_p = 0, single_n = 0;
   real_t single_dip = 0.0;
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "--quick") == 0) { quick = true; }
      if (std::strcmp(argv[i], "--single") == 0 && i + 3 < argc)
      {
         single_dip = std::atof(argv[++i]);
         single_p = std::atoi(argv[++i]);
         single_n = std::atoi(argv[++i]);
      }
   }
   if (single_p > 0)   // one P1 case, printed only
   {
      Case cs;
      cs.dip = single_dip;
      cs.p = single_p;
      cs.n = single_n;
      Print(cs, Run(comm, cs));
      return 0;
   }
   if (g_rank == 0)
   {
      std::cout << "test_2d_dislocation (2D in-plane design rev 3, phase E) on "
                << Mpi::WorldSize() << " rank(s)" << (quick ? ", --quick" : "") << "\n";
   }
   TestPlaneStrainLimit();
   TestConvergence(comm, quick);
   TestSignAndRotation(comm);
   if (g_rank == 0)
   {
      std::cout << "\ntest_2d_dislocation: " << (num_tests - num_failed) << " / " << num_tests
                << " passed" << (num_failed ? " -- FAILED" : "") << "\n";
   }
   return num_failed ? 1 : 0;
}
