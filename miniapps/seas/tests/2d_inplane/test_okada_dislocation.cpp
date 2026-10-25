// tests/2d_inplane/test_okada_dislocation.cpp — 2D in-plane design rev 3,
// phase A (section 7, row A): the static dislocation check against Okada
// (1992).  A uniform-slip rectangular patch on a fault dipping 60 degrees in
// a half-space, imposed through the QD operator's slip embedding
// (ElasticityDomainOperator in field mode, fault/fault_frame_field.hpp),
// against Okada's solution (okada92.hpp, a transcription of his DC3D):
// the surface displacement and the fault tractions, for dip slip and strike
// slip, p = 1 and 2, at two resolutions.
//
// Problem (Okada's frame: x along strike, z up, z <= 0 in the medium):
// lambda = mu = 32.04 GPa (alpha = 2/3); the square patch -1 <= xi <= 1 km
// along strike, -1 <= eta <= 1 km along dip (W = 2 km) about (0, 0, -2 km),
// so its top edge is 1.13 km deep; slip 1 m, the motion of the hanging wall
// relative to the foot wall (Okada's convention): dip slip = (0, cos, sin)
// (thrust), strike slip = (1, 0, 0).  Mesh: hexes split into 6 tets each,
// the box 0.5 km beyond the patch along strike and below it, -1 <= y' <= 1
// km and up to z = 0, sheared y = y' + (z + 2 km) cot(dip) so that the
// grid plane y' = 0 is the fault plane; n x n cells on the patch and about
// the same size outside; the patch's faces are the fault (attr 3), the top
// z = 0 is free (attr 1), the other sides carry Okada's displacement as
// Dirichlet data (attr 5): no truncation error, the difference to Okada is
// the discretisation's alone.  Solver MUMPS BLR at the operator's tolerance
// 1e-12 (full-rank MUMPS runs out of memory at n = 8, p = 2 in 15 GB; BLR
// takes about 7 GB on 1 rank).
//
// Conventions checked.  The field orients the fault by ref_normal (n .
// ref > 0); the operator's slip is u+ - u- = S_dip d + S_strike s with the
// plus side the one n points away from, so S = b (Okada's slip) when n
// points into the foot wall (nu_FW = (0, sin, -cos)) and -b otherwise.  Its
// traction (ComputeTraction) is the traction on the plus side, sigma n,
// in (d, s): tau_dip = d . sigma n, tau_strike = s . sigma n (Tandem's: the
// QD friction law puts the slip rate anti-parallel to tau_pre + tau,
// friction/device_rate_state.hpp RSSolveSlipRateVectorPsi), and the normal
// traction -n . sigma n (compression positive, sigma_n,eff = SnPre +
// normal traction, rate_state_fault.hpp).
//
//   O0  the transcription's physics, no mesh (dip 60 and 90 degrees, dip and
//       strike slip): u(hanging wall) - u(foot wall) = b across the patch
//       and 0 across the plane outside it (1e-3 |b|, points 1e-5 km off the
//       plane); sigma . e_z = 0 on z = 0 (1e-10 of mu |b| / W); div sigma =
//       0 by central differences (1e-4 of the largest stress gradient);
//       sigma . nu the same on both sides of the patch and on it (1e-3 of
//       mu |b| / W).  (In the cloud container the transcription also
//       equals Okada's Fortran DC3D, REAL*8, bit for bit at 1e5 points:
//       review page 2.)
//   O1  dip 60, ref_normal = nu_FW: for dip and strike slip, p = 1 and 2,
//       n = 4 and 8: e_u, the relative L2 error of the surface displacement
//       (top face, order 2p + 2); e_T, the relative RMS error of (tau_dip,
//       tau_strike, normal traction) over the DOFs of the fault faces that
//       lie at least 0.5 km inside the patch edges, and the same for the
//       traction's stress part {sigma} n alone.  Pass: e_u < 0.25 and the
//       stress part < 0.5 (a reversed sign gives about 2) at every n; from
//       n = 4 to 8 e_u and the stress part fall and e_T falls by more than 2.
//       The full traction (stress part + the penalty times [u] - S, what the
//       QD regime uses) converges from far above 1: the edges of a uniform
//       slip patch are singular and pollute the penalty term
//   O2  ref_normal = -nu_FW (n flipped, so S = -b): the physical driving
//       traction tau_dip d + tau_strike s reverses and the normal traction
//       stays, DOF by DOF, to 1e-9 of their largest value (n = 2, p = 1)
//   O3  rotation: mesh, Okada's field, ref_normal and up rotated by R (40
//       degrees about (1, 2, 3)): (tau_dip, tau_strike, normal traction) DOF
//       by DOF the same as unrotated to 1e-8 of the largest (n = 2, p = 1)
//   O4  BP5's frame: a vertical fault (dip 90) with ref_normal (0, -1, 0) and
//       up (0, 0, 1), strike slip, n = 4 and 8, p = 1: the same checks as O1
//       (the conventions above are BP5's)
//
// Modes: (default) O0-O4;  --quick: n = 4 only (no fall checked);
// --single dip|strike P N: one O1 case, printed (no check).
//
//   make -f tests/2d_inplane/2d_inplane.mk test-okada-dislocation
//
// Scope (Chunhui Zhao, 2026-10-25).

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
// patch (km, Okada's DC3D arguments)
static const real_t kDepth = 2.0, kAL1 = -1.0, kAL2 = 1.0, kAW1 = -1.0, kAW2 = 1.0;
static const real_t kW = (kAW2 - kAW1) * 1e3;      ///< patch width [m]
static const real_t kSigScale = kMu * 1.0 / kW;    ///< mu |b| / W, |b| = 1 m
// box (km): the patch's rectangle and kMargin around it along strike and
// below, y' in [-kMarginY, kMarginY], up to z = 0
static const real_t kMargin = 0.5, kMarginY = 1.0;
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

// ------------------------------------------------------------ Okada's field
struct Slip { real_t U1 = 0.0, U2 = 0.0; const char *name = ""; };

/// Okada's slip vector b (hanging wall minus foot wall) in Okada's frame.
static void SlipVector(const Slip &sl, real_t dip, real_t b[3])
{
   real_t sd, cd;
   SinCos(dip, sd, cd);
   b[0] = sl.U1;
   b[1] = sl.U2 * cd;
   b[2] = sl.U2 * sd;
}

/// The normal into the foot wall, (0, sin, -cos), in Okada's frame.
static void FootWallNormal(real_t dip, real_t nu[3])
{
   real_t sd, cd;
   SinCos(dip, sd, cd);
   nu[0] = 0.0; nu[1] = sd; nu[2] = -cd;
}

/// Okada's displacement [m] and gradient [1] at X [m] of the frame rotated
/// by R (X = R X_o): u = R u_o, G = R G_o R^T (G_ij = du_i / dx_j).
static int OkadaAt(const Slip &sl, real_t dip, const Rot &rot, const real_t X[3],
                   real_t u[3], real_t G[3][3])
{
   real_t Xo[3];
   for (int i = 0; i < 3; i++) { Xo[i] = rot.R[0][i] * X[0] + rot.R[1][i] * X[1] + rot.R[2][i] * X[2]; }
   const real_t zo = std::min(Xo[2], 0.0);
   const okada92::Result r = okada92::DC3D(kAlpha, Xo[0] * 1e-3, Xo[1] * 1e-3, zo * 1e-3, kDepth,
                                           dip, kAL1, kAL2, kAW1, kAW2, sl.U1, sl.U2, 0.0);
   real_t Go[3][3];
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) { Go[i][j] = r.du[j][i] * 1e-3; }   // per km -> per m
   for (int i = 0; i < 3; i++)
   {
      u[i] = rot.R[i][0] * r.u[0] + rot.R[i][1] * r.u[1] + rot.R[i][2] * r.u[2];
      for (int j = 0; j < 3; j++)
      {
         real_t s = 0.0;
         for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++) { s += rot.R[i][k] * Go[k][l] * rot.R[j][l]; }
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

static void StressAtKm(const Slip &sl, real_t dip, real_t x, real_t y, real_t z, real_t S[3][3])
{
   const real_t X[3] = {x * 1e3, y * 1e3, z * 1e3};
   real_t u[3], G[3][3];
   OkadaAt(sl, dip, Rot(), X, u, G);
   Stress(G, S);
}

// ---------------------------------------------------------------- O0
static void TestOkadaPhysics()
{
   if (g_rank == 0) { std::cout << "\n=== O0: the transcription's physics (no mesh) ===\n"; }
   const real_t eps = 1e-5;   // km off the plane
   for (real_t dip : {60.0, 90.0})
   {
      for (Slip sl : {Slip{0.0, 1.0, "dip slip"}, Slip{1.0, 0.0, "strike slip"}})
      {
         real_t sd, cd, b[3], nu[3];
         SinCos(dip, sd, cd);
         SlipVector(sl, dip, b);
         FootWallNormal(dip, nu);
         // jump across the plane, on and off the patch; traction continuity
         real_t jump_on = 0.0, jump_off = 0.0, trac = 0.0;
         for (int i = 0; i <= 8; i++)
         {
            for (int k = 0; k <= 8; k++)
            {
               const real_t xi = -3.0 + 6.0 * i / 8.0 + 0.0123, eta = -1.6 + 3.2 * k / 8.0 + 0.0071;
               const real_t P[3] = {xi, eta * cd, -kDepth + eta * sd};
               real_t Xh[3], Xf[3];
               for (int c = 0; c < 3; c++) { Xh[c] = (P[c] - eps * nu[c]) * 1e3; Xf[c] = (P[c] + eps * nu[c]) * 1e3; }
               real_t uh[3], uf[3], Gh[3][3], Gf[3][3], G0[3][3], u0[3];
               OkadaAt(sl, dip, Rot(), Xh, uh, Gh);
               OkadaAt(sl, dip, Rot(), Xf, uf, Gf);
               const real_t X0[3] = {P[0] * 1e3, P[1] * 1e3, P[2] * 1e3};
               OkadaAt(sl, dip, Rot(), X0, u0, G0);
               const bool on = (xi > kAL1 && xi < kAL2 && eta > kAW1 && eta < kAW2);
               const bool near_edge = std::abs(xi - kAL1) < 0.1 || std::abs(xi - kAL2) < 0.1 ||
                                      std::abs(eta - kAW1) < 0.1 || std::abs(eta - kAW2) < 0.1;
               if (near_edge) { continue; }
               real_t e = 0.0;
               for (int c = 0; c < 3; c++) { e = std::max(e, std::abs(uh[c] - uf[c] - (on ? b[c] : 0.0))); }
               (on ? jump_on : jump_off) = std::max(on ? jump_on : jump_off, e);
               if (on)
               {
                  real_t Sh[3][3], Sf[3][3], S0[3][3];
                  Stress(Gh, Sh); Stress(Gf, Sf); Stress(G0, S0);
                  for (int r = 0; r < 3; r++)
                  {
                     real_t th = 0, tf = 0, t0 = 0;
                     for (int c = 0; c < 3; c++) { th += Sh[r][c] * nu[c]; tf += Sf[r][c] * nu[c]; t0 += S0[r][c] * nu[c]; }
                     trac = std::max(trac, std::max(std::abs(th - tf), std::abs(t0 - tf)) / kSigScale);
                  }
               }
            }
         }
         // free surface
         real_t free = 0.0;
         for (int i = 0; i <= 12; i++)
            for (int j = 0; j <= 12; j++)
            {
               real_t S[3][3];
               StressAtKm(sl, dip, -6.0 + i + 0.031, -6.0 + j + 0.017, 0.0, S);
               for (int r = 0; r < 3; r++) { free = std::max(free, std::abs(S[r][2]) / kSigScale); }
            }
         // equilibrium: central differences of sigma, h = 1e-3 km, points >= 0.5 km from the patch
         real_t div = 0.0, grad = 0.0;
         const real_t h = 1e-3;
         for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
               for (int k = 0; k < 4; k++)
               {
                  const real_t x = -4.0 + 2.0 * i + 0.13, y = -4.0 + 2.0 * j + 0.29, z = -0.6 - 1.5 * k;
                  const real_t q = y * sd - (kDepth + z) * cd;
                  if (std::abs(q) < 0.5) { continue; }
                  real_t dS[3][3][3];   // dS[a][r][c] = d sigma_rc / d x_a
                  for (int a = 0; a < 3; a++)
                  {
                     real_t Sp[3][3], Sm[3][3];
                     const real_t dx[3] = {a == 0 ? h : 0.0, a == 1 ? h : 0.0, a == 2 ? h : 0.0};
                     StressAtKm(sl, dip, x + dx[0], y + dx[1], z + dx[2], Sp);
                     StressAtKm(sl, dip, x - dx[0], y - dx[1], z - dx[2], Sm);
                     for (int r = 0; r < 3; r++)
                        for (int c = 0; c < 3; c++)
                        {
                           dS[a][r][c] = (Sp[r][c] - Sm[r][c]) / (2.0 * h);
                           grad = std::max(grad, std::abs(dS[a][r][c]));
                        }
                  }
                  for (int r = 0; r < 3; r++)
                  {
                     div = std::max(div, std::abs(dS[0][r][0] + dS[1][r][1] + dS[2][r][2]));
                  }
               }
         std::ostringstream tag;
         tag << "O0 dip " << dip << ", " << sl.name << ": ";
         TEST_ASSERT(jump_on <= 1e-3, tag.str() << "u(hanging wall) - u(foot wall) = b on the patch ("
                     << jump_on << " m)");
         TEST_ASSERT(jump_off <= 1e-3, tag.str() << "no jump across the plane off the patch ("
                     << jump_off << " m)");
         TEST_ASSERT(free <= 1e-10, tag.str() << "sigma . e_z = 0 on z = 0 (" << free
                     << " of mu |b| / W)");
         TEST_ASSERT(div <= 1e-4 * grad, tag.str() << "div sigma = 0 (" << div / grad
                     << " of the largest stress gradient)");
         TEST_ASSERT(trac <= 1e-3, tag.str() << "sigma . nu continuous across the patch (" << trac
                     << " of mu |b| / W)");
      }
   }
}

// ------------------------------------------------------------- the FE side
/// The sheared box (see the file comment), rotated by R after tagging.
static Mesh *OkadaBox(int n, real_t dip, const Rot &rot)
{
   real_t sd, cd;
   SinCos(dip, sd, cd);
   const real_t hx = (kAL2 - kAL1) / n, hz = (kAW2 - kAW1) * sd / n;
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
   const std::vector<real_t> xs = Nodes({{kAL1 - kMargin, kAL1}, {kAL1, kAL2}, {kAL2, kAL2 + kMargin}},
                                        hx, {0, n, 0});
   const std::vector<real_t> ys = Nodes({{-kMarginY, 0.0}, {0.0, kMarginY}}, hx, {0, 0});
   const std::vector<real_t> zs = Nodes({{zb - kMargin, zb}, {zb, zt}, {zt, 0.0}}, hz, {0, n, 0});
   const int nx = (int)xs.size() - 1, ny = (int)ys.size() - 1, nz = (int)zs.size() - 1;
   int ia = -1, ib = -1, j0 = -1, kb = -1, kt = -1;
   for (int i = 0; i <= nx; i++) { if (xs[i] == kAL1) { ia = i; } if (xs[i] == kAL2) { ib = i; } }
   for (int j = 0; j <= ny; j++) { if (ys[j] == 0.0) { j0 = j; } }
   for (int k = 0; k <= nz; k++) { if (zs[k] == zb) { kb = k; } if (zs[k] == zt) { kt = k; } }
   MFEM_VERIFY(ia >= 0 && ib > ia && j0 >= 0 && kb >= 0 && kt > kb, "OkadaBox: patch grid lines");

   Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(nx, ny, nz, Element::TETRAHEDRON, nx, ny, nz));
   auto Grid = [&](int v, int g[3])
   {
      const real_t *c = mesh->GetVertex(v);
      for (int d = 0; d < 3; d++) { g[d] = (int)std::lround(c[d]); }
   };
   for (int be = 0; be < mesh->GetNBE(); be++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(be, verts);
      bool top = true;
      for (int i = 0; i < verts.Size(); i++) { int g[3]; Grid(verts[i], g); if (g[2] != nz) { top = false; } }
      mesh->SetBdrAttribute(be, top ? 1 : 5);
   }
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool in = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         int g[3]; Grid(fv[i], g);
         if (g[1] != j0 || g[0] < ia || g[0] > ib || g[2] < kb || g[2] > kt) { in = false; }
      }
      if (in) { tris.push_back({{fv[0], fv[1], fv[2]}}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   const real_t cot = (sd != 0.0) ? cd / sd : 0.0;
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      int g[3]; Grid(v, g);
      const real_t x = xs[g[0]], z = zs[g[2]], y = ys[g[1]] + (z + kDepth) * cot;
      const real_t X[3] = {x * 1e3, y * 1e3, z * 1e3};
      rot.Apply(X, mesh->GetVertex(v));
   }
   return mesh;
}

struct Case
{
   real_t dip = 60.0;
   Slip slip;
   int p = 1, n = 2;
   real_t ref_sign = 1.0;            ///< ref_normal = ref_sign nu_FW (unless bp5_frame)
   bool bp5_frame = false;           ///< ref_normal (0, -1, 0)
   Rot rot;
};

struct Result
{
   real_t e_u = 0.0, e_T = 0.0;
   real_t e_Ts = 0.0;   ///< the same for the stress part {sigma} n alone (diagnostic)
   real_t jump = 0.0;   ///< largest |[u] - S| at the region's DOFs [m] (diagnostic)
   int n_dofs = 0, n_trac_dofs = 0, n_fault_faces = 0;
   std::vector<real_t> trac;    ///< owned fault DOFs: tau_dip, tau_strike, normal traction
   std::vector<real_t> tphys;   ///< owned fault DOFs: tau_dip d + tau_strike s (3), normal traction
};

using Op = ElasticityDomainOperator<ParMesh>;

static Result Run(MPI_Comm comm, const Case &cs)
{
   std::unique_ptr<Mesh> m(OkadaBox(cs.n, cs.dip, cs.rot));
   ParMesh pm(comm, *m);
   m.reset();
   real_t nu0[3], nu[3], up0[3] = {0, 0, 1}, up[3], ref[3];
   FootWallNormal(cs.dip, nu0);
   cs.rot.Apply(nu0, nu);
   cs.rot.Apply(up0, up);
   if (cs.bp5_frame) { const real_t r0[3] = {0, -1, 0}; cs.rot.Apply(r0, ref); }
   else { for (int d = 0; d < 3; d++) { ref[d] = cs.ref_sign * nu[d]; } }

   const Slip sl = cs.slip;
   const real_t dip = cs.dip;
   const Rot rot = cs.rot;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [sl, dip, rot](const Vector &x, real_t, Vector &u)
   {
      u.SetSize(3);
      const real_t X[3] = {x(0), x(1), x(2)};
      real_t uu[3], G[3][3];
      OkadaAt(sl, dip, rot, X, uu, G);
      for (int d = 0; d < 3; d++) { u(d) = uu[d]; }
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
   // MUMPS with block low-rank compression at 1e-12, the cycle QD regime's solver
   // (full rank needs about 13 GB at n = 8, p = 2)
   Op op(pm, cs.p, model, 0.0, 1e4, 1e4, bc, DGMethod::IP, SolverType::MUMPS_BLR, dc);
   auto F = op.FrameField();
   MFEM_VERIFY(F, "Okada: no frame field");

   // slip: u+ - u- = J, J = b when n points into the foot wall, else -b
   real_t b0[3], b[3];
   SlipVector(sl, dip, b0);
   rot.Apply(b0, b);
   const int nd = op.GetNumFaultDOFs();
   const int nf = F->NumFaultFaces();
   MFEM_VERIFY(nf == 0 || nd % nf == 0, "Okada: DOFs per face");
   const int nbf = nf ? nd / nf : 1;
   Vector slip(2 * nd);
   for (int k = 0; k < nd; k++)
   {
      const fault_frame::Frame &fr = F->Frame(k / nbf);
      const real_t side = (fr.n[0] * nu[0] + fr.n[1] * nu[1] + fr.n[2] * nu[2] > 0.0) ? 1.0 : -1.0;
      slip(2 * k) = side * (b[0] * fr.d[0] + b[1] * fr.d[1] + b[2] * fr.d[2]);
      slip(2 * k + 1) = side * (b[0] * fr.s[0] + b[1] * fr.s[1] + b[2] * fr.s[2]);
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
   auto XiEta = [&](int k, real_t &xi, real_t &eta)
   {
      const real_t X[3] = {xyz(3 * k), xyz(3 * k + 1), xyz(3 * k + 2)};
      real_t Xo[3];
      for (int i = 0; i < 3; i++) { Xo[i] = (rot.R[0][i] * X[0] + rot.R[1][i] * X[1] + rot.R[2][i] * X[2]) * 1e-3; }
      xi = Xo[0];
      eta = Xo[1] * cd + (Xo[2] + kDepth) * sd;
   };
   Vector err2(nd), errs2(nd), ref2(nd), tl(3 * nd), tp(4 * nd), cnt(nd), jmax(nd);
   err2 = 0.0; errs2 = 0.0; ref2 = 0.0; cnt = 0.0; jmax = 0.0;
   for (int f = 0; f < (nbf ? nd / nbf : 0); f++)
   {
      bool in = true;
      for (int kk = 0; kk < nbf; kk++)
      {
         real_t xi, eta;
         XiEta(f * nbf + kk, xi, eta);
         if (xi < kAL1 + kEdgeGap - 1e-9 || xi > kAL2 - kEdgeGap + 1e-9 ||
             eta < kAW1 + kEdgeGap - 1e-9 || eta > kAW2 - kEdgeGap + 1e-9) { in = false; }
      }
      for (int kk = 0; kk < nbf; kk++)
      {
         const int k = f * nbf + kk;
         const fault_frame::Frame &fr = F->Frame(f);
         const real_t th[3] = {trac(2 * k), trac(2 * k + 1), sn(k)};
         const real_t ts[3] = {trac_s(2 * k), trac_s(2 * k + 1), sn_s(k)};
         for (int c = 0; c < 3; c++) { tl(3 * k + c) = th[c]; tp(4 * k + c) = th[0] * fr.d[c] + th[1] * fr.s[c]; }
         tp(4 * k + 3) = th[2];
         if (!in) { continue; }
         const real_t X[3] = {xyz(3 * k), xyz(3 * k + 1), xyz(3 * k + 2)};
         real_t uu[3], G[3][3], S[3][3], t[3];
         OkadaAt(sl, dip, rot, X, uu, G);
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
   res.trac.assign(tl_o.GetData(), tl_o.GetData() + tl_o.Size());
   res.tphys.assign(tp_o.GetData(), tp_o.GetData() + tp_o.Size());

   // surface displacement against Okada's, L2 over the top face
   real_t e2 = 0.0, r2 = 0.0;
   Vector val(3), xq(3);
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
         const real_t X[3] = {xq(0), xq(1), xq(2)};
         real_t uu[3], G[3][3];
         OkadaAt(sl, dip, rot, X, uu, G);
         const real_t w = ip.weight * FTr->Weight();
         for (int c = 0; c < 3; c++) { e2 += w * (val(c) - uu[c]) * (val(c) - uu[c]); r2 += w * uu[c] * uu[c]; }
      }
   }
   res.e_u = std::sqrt(RSum(e2, comm) / RSum(r2, comm));
   return res;
}

static void Print(const Case &cs, const Result &r)
{
   if (g_rank != 0) { return; }
   std::cout << "    dip " << std::setw(2) << cs.dip << "  " << std::left << std::setw(11) << cs.slip.name
             << std::right << "  p " << cs.p << "  n " << cs.n << "  DOFs " << std::setw(7) << r.n_dofs
             << "  fault faces " << std::setw(4) << r.n_fault_faces << "  traction DOFs " << std::setw(4)
             << r.n_trac_dofs << "  e_u " << std::scientific << std::setprecision(3) << r.e_u
             << "  e_T " << r.e_T << "  (stress part " << r.e_Ts << ", |[u] - S| <= " << r.jump
             << " m)" << std::defaultfloat << "\n";
}

static real_t MaxAbs(const std::vector<real_t> &v, MPI_Comm comm)
{
   real_t m = 0.0;
   for (real_t x : v) { m = std::max(m, std::abs(x)); }
   return RMax(m, comm);
}

// ---------------------------------------------------------------- O1, O4
static void TestConvergence(MPI_Comm comm, bool quick, bool bp5_frame)
{
   if (g_rank == 0)
   {
      std::cout << (bp5_frame ? "\n=== O4: BP5's frame: vertical fault, ref_normal (0, -1, 0) ===\n"
                              : "\n=== O1: dip 60, ref_normal = nu_FW: convergence to Okada ===\n");
   }
   std::vector<Slip> slips = bp5_frame ? std::vector<Slip>{Slip{1.0, 0.0, "strike slip"}}
                                       : std::vector<Slip>{Slip{0.0, 1.0, "dip slip"}, Slip{1.0, 0.0, "strike slip"}};
   std::vector<int> ps = bp5_frame ? std::vector<int>{1} : std::vector<int>{1, 2};
   std::vector<int> ns = quick ? std::vector<int>{4} : std::vector<int>{4, 8};
   for (const Slip &sl : slips)
   {
      for (int p : ps)
      {
         std::vector<Result> rs;
         for (int n : ns)
         {
            Case cs;
            cs.dip = bp5_frame ? 90.0 : 60.0;
            cs.slip = sl;
            cs.p = p;
            cs.n = n;
            cs.bp5_frame = bp5_frame;
            rs.push_back(Run(comm, cs));
            Print(cs, rs.back());
         }
         std::ostringstream tag;
         tag << (bp5_frame ? "O4 " : "O1 ") << sl.name << ", p = " << p << ": ";
         for (size_t i = 0; i < rs.size(); i++)
         {
            TEST_ASSERT(rs[i].e_u < 0.25 && rs[i].e_Ts < 0.5,
                        tag.str() << "n = " << ns[i] << ": e_u " << rs[i].e_u << " (< 0.25), stress "
                        "part " << rs[i].e_Ts << " (< 0.5; a reversed sign gives about 2)");
         }
         if (rs.size() == 2)
         {
            TEST_ASSERT(rs[1].e_u < rs[0].e_u && rs[1].e_Ts < rs[0].e_Ts && rs[1].e_T < 0.5 * rs[0].e_T,
                        tag.str() << "n = 4 -> 8: e_u / " << rs[0].e_u / rs[1].e_u << ", stress part / "
                        << rs[0].e_Ts / rs[1].e_Ts << ", e_T / " << rs[0].e_T / rs[1].e_T << " (> 2)");
         }
      }
   }
}

// ---------------------------------------------------------------- O2, O3
static void TestSignAndRotation(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== O2, O3: the normal reversed; the problem rotated ===\n"; }
   for (Slip sl : {Slip{0.0, 1.0, "dip slip"}, Slip{1.0, 0.0, "strike slip"}})
   {
      Case base;
      base.slip = sl;
      const Result r0 = Run(comm, base);
      Print(base, r0);

      Case flip = base;
      flip.ref_sign = -1.0;
      const Result r1 = Run(comm, flip);
      real_t e = 0.0;
      for (size_t k = 0; k < r0.tphys.size(); k += 4)
      {
         for (int c = 0; c < 3; c++) { e = std::max(e, std::abs(r1.tphys[k + c] + r0.tphys[k + c])); }
         e = std::max(e, std::abs(r1.tphys[k + 3] - r0.tphys[k + 3]));
      }
      e = RMax(e, comm) / MaxAbs(r0.tphys, comm);
      TEST_ASSERT(r1.tphys.size() == r0.tphys.size() && e <= 1e-9,
                  "O2 " << sl.name << ": ref_normal -nu_FW: driving traction reversed, normal "
                  "traction the same (" << e << " of the largest)");

      Case rotc = base;
      rotc.rot = Rot::Axis(1.0, 2.0, 3.0, 40.0);
      const Result r2 = Run(comm, rotc);
      real_t er = 0.0;
      for (size_t k = 0; k < r0.trac.size(); k++) { er = std::max(er, std::abs(r2.trac[k] - r0.trac[k])); }
      er = RMax(er, comm) / MaxAbs(r0.trac, comm);
      TEST_ASSERT(r2.trac.size() == r0.trac.size() && er <= 1e-8,
                  "O3 " << sl.name << ": rotated: (tau_dip, tau_strike, normal traction) the same ("
                  << er << " of the largest; e_u " << r2.e_u << " vs " << r0.e_u << ")");
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   g_rank = Mpi::WorldRank();
   bool quick = false;
   int single_p = 0, single_n = 0;
   const char *single_slip = nullptr;
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "--quick") == 0) { quick = true; }
      if (std::strcmp(argv[i], "--single") == 0 && i + 3 < argc)
      {
         single_slip = argv[++i];
         single_p = std::atoi(argv[++i]);
         single_n = std::atoi(argv[++i]);
      }
   }
   if (single_slip)   // one O1 case, printed only
   {
      Case cs;
      cs.slip = (std::strcmp(single_slip, "strike") == 0) ? Slip{1.0, 0.0, "strike slip"}
                                                          : Slip{0.0, 1.0, "dip slip"};
      cs.p = single_p;
      cs.n = single_n;
      Print(cs, Run(comm, cs));
      return 0;
   }
   if (g_rank == 0)
   {
      std::cout << "test_okada_dislocation (2D in-plane design rev 3, phase A) on "
                << Mpi::WorldSize() << " rank(s)" << (quick ? ", --quick" : "") << "\n";
   }
   TestOkadaPhysics();
   TestConvergence(comm, quick, false);
   TestSignAndRotation(comm);
   TestConvergence(comm, quick, true);
   if (g_rank == 0)
   {
      std::cout << "\ntest_okada_dislocation: " << (num_tests - num_failed) << " / " << num_tests
                << " passed" << (num_failed ? " -- FAILED" : "") << "\n";
   }
   return num_failed ? 1 : 0;
}
