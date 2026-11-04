// tests/pml/test_pml_band.cpp - PML plan rev2, P1: the PmlBand unit test and
// the band integrators.
//
//   A. walls from a tilted mesh: attribute 5 holds two opposite walls, the
//      normals and planes match the rotation (2-D at azimuth 123 deg, 3-D
//      tilted cube); corner (2-D) and edge/corner (3-D) detection
//   B. xi_a at random points against the distance to each wall
//   C. C0 : G of ApplyC0 against the tensor (2.22) built entry by entry;
//      major symmetry; C0 = C off the band; the axis form S/(s_j s_l) C_ijkl
//   D. kappa (3.14) on face bands, 2-D corners, 3-D edges and corners
//   E. element and face matrices outside the bands are the plain ones bit for
//      bit (T0 at the integrator level); band matrices are symmetric
//   F. patch test (Proposition 1): u = E x~(x), x~ = x + sum_a D_a(xi_a)/alpha n_a
//      is an exact static solution of div(C0 grad u) = 0.  With p = 3 it lies
//      in the DG space, so K0 u - b(g = u) = 0 to round-off (2-D and 3-D,
//      tilted boxes, every face Dirichlet).  The same u is not a solution of
//      the plain operator (the check has teeth).
//
// Serial, no MPI ranks needed.  Exit code 0 = all pass.

#include "mfem.hpp"
#include "miniapps/seas/domain/pml_band.hpp"
#include "miniapps/seas/integrator/pml_elasticity_integrators.hpp"
#include "miniapps/seas/integrator/dg_elasticity_ip_combined_integrator.hpp"

#include <cmath>
#include <cstdio>
#include <random>

using namespace mfem;
using namespace mfem::seas;

static int g_fail = 0;
static void Check(bool ok, const char *what, double val = 0.0)
{
   std::printf("  [%s] %s  (%.3e)\n", ok ? "PASS" : "FAIL", what, val);
   if (!ok) { g_fail++; }
}

// rotation of the reference box: 2-D by angle th; 3-D by R = Rz(th) Rx(ph)
struct Rot
{
   int D;
   double R[3][3];
   Rot(int dim, double th, double ph) : D(dim)
   {
      double c = std::cos(th), s = std::sin(th), cp = std::cos(ph), sp = std::sin(ph);
      double Rz[3][3] = {{c, -s, 0}, {s, c, 0}, {0, 0, 1}};
      double Rx[3][3] = {{1, 0, 0}, {0, cp, -sp}, {0, sp, cp}};
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++)
         {
            R[i][j] = 0.0;
            for (int k = 0; k < 3; k++) { R[i][j] += Rz[i][k] * (dim == 3 ? Rx[k][j] : (k == j)); }
         }
   }
   void Apply(const double *a, double *b) const
   {
      for (int i = 0; i < D; i++)
      {
         b[i] = 0.0;
         for (int j = 0; j < D; j++) { b[i] += R[i][j] * a[j]; }
      }
   }
};

// box [0,W]^D (2-D: W x H), shifted by -W/2, rotated by rot, shifted by x0
static Mesh MakeBox(int D, int n, double W, double H, const Rot &rot, const double *x0)
{
   Mesh m = (D == 2) ? Mesh::MakeCartesian2D(n, n * H / W, Element::TRIANGLE, true, W, H)
                     : Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON, W, W, W);
   for (int v = 0; v < m.GetNV(); v++)
   {
      double *x = m.GetVertex(v);
      double a[3] = {x[0] - 0.5 * W, x[1] - 0.5 * (D == 2 ? H : W), D == 3 ? x[2] - 0.5 * W : 0.0};
      double b[3];
      rot.Apply(a, b);
      for (int i = 0; i < D; i++) { x[i] = b[i] + x0[i]; }
   }
   return m;
}

static PmlBandParams Params(double L)
{
   PmlBandParams p;
   p.thickness = L;
   p.power = 2;
   p.R0 = 1e-3;
   p.c_ref = 6000.0;
   p.f_alpha = 0.3 * p.c_ref / L;   // s0 at the wall ~ 4.7 (moderate, for the tests)
   return p;
}

// C_ijkl of (1.3)
static double Ciso(int i, int j, int k, int l, double lam, double mu)
{
   return lam * (i == j) * (k == l) + mu * ((i == k) * (j == l) + (i == l) * (j == k));
}

// ---------------------------------------------------------------- A, B, D
static void TestGeometry2D()
{
   std::printf("A/B/D 2-D: tilted box at azimuth 123 deg\n");
   const double th = 123.0 * M_PI / 180.0, W = 100e3, H = 60e3, L = 20e3;
   const double x0[3] = {3e5, -2e5, 0};
   Rot rot(2, th, 0.0);
   Mesh m = MakeBox(2, 10, W, H, rot, x0);
   // MakeCartesian2D: 1 bottom, 2 right, 3 top, 4 left.  Bottom and top -> 5.
   for (int be = 0; be < m.GetNBE(); be++)
   {
      if (m.GetBdrAttribute(be) == 1 || m.GetBdrAttribute(be) == 3) { m.SetBdrAttribute(be, 5); }
   }
   m.SetAttributes();
   PmlBand band(m, {5, 2}, Params(L));
   const auto &w = band.Walls();
   Check(w.size() == 3, "three walls (attribute 5 split into two, attribute 2)", w.size());
   // expected normals: bottom -e2, top +e2, right +e1, rotated
   double e[3][2] = {{0, -1}, {0, 1}, {1, 0}};
   int attr_exp[3] = {5, 5, 2};
   double off[3] = {0.5 * H, 0.5 * H, 0.5 * W};
   double worst = 0.0, worst_h = 0.0;
   int found = 0;
   for (int k = 0; k < 3; k++)
   {
      double n[2];
      rot.Apply(e[k], n);
      for (const PmlWall &ww : w)
      {
         if (ww.attr != attr_exp[k] || n[0] * ww.n[0] + n[1] * ww.n[1] < 0.99) { continue; }
         found++;
         worst = std::max({worst, std::abs(ww.n[0] - n[0]), std::abs(ww.n[1] - n[1])});
         const double h = off[k] + n[0] * x0[0] + n[1] * x0[1];
         worst_h = std::max(worst_h, std::abs(ww.h - h));
      }
   }
   Check(found == 3, "each expected wall found once", found);
   Check(worst < 1e-12, "wall normals match the rotation", worst);
   Check(worst_h < 1e-6, "wall planes match (m)", worst_h);

   // B: xi at random points; D: kappa
   std::mt19937 gen(7);
   std::uniform_real_distribution<double> U(-0.5, 0.5);
   double err_xi = 0.0, err_k = 0.0;
   int ncorner = 0, nface = 0, nin = 0;
   const double s0w = band.S0Wall();
   for (int t = 0; t < 4000; t++)
   {
      double a[3] = {W * U(gen), H * U(gen), 0}, b[3], x[2];
      rot.Apply(a, b);
      x[0] = b[0] + x0[0];
      x[1] = b[1] + x0[1];
      PmlPoint p = band.At(x);
      // reference: distances in the box frame
      double dist[3] = {0.5 * H + a[1], 0.5 * H - a[1], 0.5 * W - a[0]};
      int nact = 0;
      for (int k = 0; k < 3; k++)
      {
         if (L - dist[k] > 0) { nact++; }
      }
      if (nact != p.nact) { err_xi = 1.0; continue; }
      double sx[3] = {1, 1, 1};
      for (int k = 0; k < p.nact; k++)
      {
         // match the wall by normal
         const PmlWall &ww = w[p.wall[k]];
         for (int j = 0; j < 3; j++)
         {
            double n[2];
            rot.Apply(e[j], n);
            if (ww.attr == attr_exp[j] && n[0] * ww.n[0] + n[1] * ww.n[1] > 0.99)
            {
               err_xi = std::max(err_xi, std::abs(p.xi[k] - (L - dist[j])) / L);
            }
         }
         sx[k] = 1.0 + band.D0() * (p.xi[k] / L) * (p.xi[k] / L) / band.Alpha();
      }
      double kap = 1.0;
      if (p.nact == 1) { kap = sx[0]; nface++; }
      if (p.nact == 2) { kap = std::max(sx[0], sx[1]); ncorner++; }
      if (p.nact == 0) { nin++; }
      err_k = std::max(err_k, std::abs(p.kappa - kap) / kap);
   }
   Check(err_xi < 1e-9, "xi_a = L - dist(x, wall a) at 4000 random points (rel.)", err_xi);
   Check(nface > 0 && ncorner > 0 && nin > 0, "points in face bands, corners and inside",
         ncorner);
   Check(err_k < 1e-12, "kappa: s0 on a face band, max(s0_1, s0_2) in a 2-D corner", err_k);
   const double d0 = 3 * 6000.0 * std::log(1e3) / (2 * L);
   Check(std::abs(band.D0() - d0) < 1e-12 * d0, "d0 = (N+1) c ln(1/R0) / (2L)", band.D0());
   Check(std::abs(band.DeltaL() - 6000.0 * std::log(1e3) / (2 * band.Alpha()))
         < 1e-9 * band.DeltaL(), "Delta L = c ln(1/R0) / (2 alpha) (4.4)", band.DeltaL());
   (void)s0w;
}

static void TestGeometry3D()
{
   std::printf("A/D 3-D: tilted cube, every face a wall\n");
   const double W = 10.0, L = 3.0;
   const double x0[3] = {1, 2, 3};
   Rot rot(3, 0.4, 0.7);
   Mesh m = MakeBox(3, 4, W, W, rot, x0);
   std::vector<int> attrs = {1, 2, 3, 4, 5, 6};
   PmlBand band(m, attrs, Params(L));
   Check(band.Walls().size() == 6, "six walls", band.Walls().size());
   // a point near a corner: in three bands; near an edge: two; face: one
   double a_c[3] = {4.5, 4.5, 4.5}, a_e[3] = {4.5, 4.5, 0.0}, a_f[3] = {4.5, 0.0, 0.0};
   double *as[3] = {a_c, a_e, a_f};
   int nexp[3] = {3, 2, 1};
   for (int t = 0; t < 3; t++)
   {
      double b[3];
      rot.Apply(as[t], b);
      for (int i = 0; i < 3; i++) { b[i] += x0[i]; }
      PmlPoint p = band.At(b);
      double s = 1.0 + band.D0() * std::pow((L - 0.5) / L, 2) / band.Alpha();
      double kap = (t == 0) ? s * s : (t == 1 ? s * s : s);
      Check(p.nact == nexp[t], t == 0 ? "3-D corner: 3 bands" : (t == 1 ? "3-D edge: 2 bands"
            : "3-D face: 1 band"), p.nact);
      Check(std::abs(p.kappa - kap) < 1e-12 * kap, t == 0 ? "kappa corner = S0/s_min"
            : (t == 1 ? "kappa edge = s1 s2" : "kappa face = s0"), p.kappa);
      Check(std::abs(p.S0 - std::pow(s, nexp[t])) < 1e-12 * p.S0, "S0 = prod s0_a", p.S0);
   }
   // perpendicular-wall refusal is exercised by construction (no throw here)
}

// ---------------------------------------------------------------- C
static void TestTensor(int D)
{
   std::printf("C %d-D: C0 : G against (2.22)\n", D);
   const double W = 10.0, L = 3.0;
   Rot rot(D, 0.9, 0.3);
   std::vector<PmlWall> walls;
   for (int a = 0; a < D; a++)
      for (int sgn = -1; sgn <= 1; sgn += 2)
      {
         double e[3] = {0, 0, 0}, n[3] = {0, 0, 0};
         e[a] = sgn;
         rot.Apply(e, n);
         PmlWall w;
         for (int i = 0; i < D; i++) { w.n[i] = n[i]; }
         w.h = 0.5 * W;
         walls.push_back(w);
      }
   PmlBand band(D, walls, Params(L));
   const double lam = 3.2e10, mu = 2.7e10;
   std::mt19937 gen(11);
   std::uniform_real_distribution<double> U(-0.5, 0.5);
   double err = 0.0, asym = 0.0, err_off = 0.0, minor = 0.0, err_axis = 0.0;
   for (int t = 0; t < 500; t++)
   {
      double a[3] = {W * U(gen), W * U(gen), D == 3 ? W * U(gen) : 0.0}, x[3] = {0, 0, 0};
      rot.Apply(a, x);
      PmlPoint p = band.At(x);
      // C0 entry by entry (2.22): S Linv_jq Linv_lr C_iqkr
      double C0[3][3][3][3];
      for (int i = 0; i < D; i++)
         for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
               for (int l = 0; l < D; l++)
               {
                  double v = 0.0;
                  for (int q = 0; q < D; q++)
                     for (int r = 0; r < D; r++)
                     { v += p.Linv[j][q] * p.Linv[l][r] * Ciso(i, q, k, r, lam, mu); }
                  C0[i][j][k][l] = p.S0 * v;
               }
      for (int i = 0; i < D; i++)
         for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
               for (int l = 0; l < D; l++)
               {
                  asym = std::max(asym, std::abs(C0[i][j][k][l] - C0[k][l][i][j]) / mu);
                  if (p.nact == 0)
                  {
                     err_off = std::max(err_off, std::abs(C0[i][j][k][l] - Ciso(i, j, k, l, lam, mu)) / mu);
                  }
                  if (p.nact > 0) { minor = std::max(minor, std::abs(C0[i][j][k][l] - C0[j][i][k][l]) / mu); }
               }
      double G[3][3] = {{0}}, sig[3][3];
      for (int i = 0; i < D; i++)
         for (int j = 0; j < D; j++) { G[i][j] = U(gen); }
      band.ApplyC0(p, lam, mu, G, sig);
      for (int i = 0; i < D; i++)
         for (int j = 0; j < D; j++)
         {
            double v = 0.0;
            for (int k = 0; k < D; k++)
               for (int l = 0; l < D; l++) { v += C0[i][j][k][l] * G[k][l]; }
            err = std::max(err, std::abs(sig[i][j] - v) / (mu * p.S0 * 4));
         }
   }
   // axis-aligned walls: C0_ijkl = S0 / (s_j s_l) C_ijkl
   {
      std::vector<PmlWall> aw;
      for (int a = 0; a < D; a++)
         for (int sgn = -1; sgn <= 1; sgn += 2)
         {
            PmlWall w;
            w.n[a] = sgn;
            w.h = 0.5 * W;
            aw.push_back(w);
         }
      PmlBand ab(D, aw, Params(L));
      for (int t = 0; t < 200; t++)
      {
         double x[3] = {W * U(gen), W * U(gen), D == 3 ? W * U(gen) : 0.0};
         PmlPoint p = ab.At(x);
         double s[3] = {1, 1, 1};
         for (int k = 0; k < p.nact; k++)
         {
            const PmlWall &w = aw[p.wall[k]];
            for (int i = 0; i < D; i++)
            {
               if (std::abs(w.n[i]) > 0.5) { s[i] = p.s0[k]; }
            }
         }
         for (int k = 0; k < D; k++)
            for (int l = 0; l < D; l++)
            {
               double G[3][3] = {{0}}, sig[3][3];
               G[k][l] = 1.0;
               ab.ApplyC0(p, lam, mu, G, sig);
               for (int i = 0; i < D; i++)
                  for (int j = 0; j < D; j++)
                  {
                     double v = p.S0 / (s[j] * s[l]) * Ciso(i, j, k, l, lam, mu);
                     err_axis = std::max(err_axis, std::abs(sig[i][j] - v) / (mu * p.S0));
                  }
            }
      }
   }
   Check(err < 1e-14, "ApplyC0 (three stages) = C0 : G of (2.22)", err);
   Check(asym < 1e-12, "major symmetry C0_ijkl = C0_klij", asym);
   Check(err_off < 1e-15, "C0 = C off the band", err_off);
   Check(minor > 1e-3, "minor symmetry lost in the band (expected)", minor);
   Check(err_axis < 1e-14, "axis walls: C0_ijkl = S0/(s_j s_l) C_ijkl", err_axis);
}

// ---------------------------------------------------------------- E, F
struct XTildeCoef : public VectorCoefficient
{
   const PmlBand &band;
   double E[3][3];
   XTildeCoef(const PmlBand &b, int D) : VectorCoefficient(D), band(b)
   {
      // a full, non-symmetric constant matrix
      double v[9] = {1.0e-3, -0.4e-3, 0.7e-3, 0.3e-3, -0.8e-3, 0.2e-3, 0.5e-3, 0.6e-3, -0.9e-3};
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++) { E[i][j] = v[3 * i + j]; }
   }
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      const int D = vdim;
      Vector x(D);
      T.Transform(ip, x);
      double xt[3] = {0, 0, 0};
      for (int i = 0; i < D; i++) { xt[i] = x(i); }
      PmlPoint p = band.At(x);
      const double L = band.Params().thickness;
      const int N = band.Params().power;
      for (int k = 0; k < p.nact; k++)
      {
         // D_a(xi) / alpha, D_a(xi) = d0 L / (N + 1) (xi / L)^(N + 1)
         const double Da = band.D0() * L / (N + 1) * std::pow(p.xi[k] / L, N + 1);
         for (int i = 0; i < D; i++) { xt[i] += Da / band.Alpha() * band.Walls()[p.wall[k]].n[i]; }
      }
      V.SetSize(D);
      for (int i = 0; i < D; i++)
      {
         V(i) = 0.0;
         for (int j = 0; j < D; j++) { V(i) += E[i][j] * xt[j]; }
      }
   }
};

static double PatchResidual(Mesh &m, const PmlBand *band, int p, bool use_band_ops,
                            int &n_band_el, double &unorm)
{
   const int D = m.Dimension();
   ConstantCoefficient lam(3.2e10), mu(2.7e10);
   L2_FECollection fec(p, D);
   FiniteElementSpace fes(&m, &fec, D);
   Array<int> all(m.bdr_attributes.Max());
   all = 1;
   BilinearForm a(&fes);
   const double eps = -1.0, pf = 1.0;
   if (use_band_ops) { a.AddDomainIntegrator(new PmlBandElasticityIntegrator(lam, mu, *band)); }
   else { a.AddDomainIntegrator(new ElasticityIntegrator(lam, mu)); }
   auto *fi = new DGElasticityIPCombinedIntegrator(lam, mu, D, eps, pf);
   auto *bi = new DGElasticityIPCombinedIntegrator(lam, mu, D, eps, pf);
   if (use_band_ops) { fi->SetPmlBand(band); bi->SetPmlBand(band); }
   a.AddInteriorFaceIntegrator(fi);
   a.AddBdrFaceIntegrator(bi, all);
   a.Assemble(0);
   a.Finalize();

   XTildeCoef uc(*band, D);
   GridFunction u(&fes);
   u.ProjectCoefficient(uc);
   unorm = u.Normlinf();

   // Dirichlet RHS with g = u on every boundary face
   Vector b(fes.GetVSize());
   b = 0.0;
   DGElasticityIPCombinedIntegrator rhs_int(lam, mu, D, eps, pf);
   if (use_band_ops) { rhs_int.SetPmlBand(band); }
   Array<int> vdofs;
   for (int be = 0; be < m.GetNBE(); be++)
   {
      FaceElementTransformations *FTr = m.GetBdrFaceTransformations(be);
      const FiniteElement *fe = fes.GetFE(FTr->Elem1No);
      const int order = 2 * fe->GetOrder() + 1;
      const IntegrationRule &ir = IntRules.Get(FTr->GetGeometryType(), order);
      const int nq = ir.GetNPoints();
      Vector gq(D * nq), gv(D);
      for (int q = 0; q < nq; q++)
      {
         FTr->SetAllIntPoints(&ir.IntPoint(q));
         uc.Eval(gv, *FTr->Elem1, FTr->GetElement1IntPoint());
         for (int c = 0; c < D; c++) { gq(c * nq + q) = gv(c); }
      }
      Vector el;
      rhs_int.AssembleBoundaryFaceRHS(*fe, *FTr, gq, el);
      fes.GetElementVDofs(FTr->Elem1No, vdofs);
      b.AddElementVector(vdofs, el);
   }
   Vector r(fes.GetVSize());
   a.Mult(u, r);
   const double bnorm = b.Normlinf();
   r -= b;
   n_band_el = 0;
   for (int e = 0; e < m.GetNE(); e++) { n_band_el += PmlBandElement(m, e, *band); }
   // scale: the size of one term of K u
   Vector au(fes.GetVSize());
   a.Mult(u, au);
   return r.Normlinf() / std::max(au.Normlinf(), bnorm);
}

static void TestPatch(int D)
{
   std::printf("F %d-D: patch test u = E x~(x), p = 3, tilted box, all walls Dirichlet\n", D);
   // element edges on the bands' inner edges, so u is a cubic on every element
   const double W = 10.0, L = (D == 2) ? 3.0 : 4.0;
   const double x0[3] = {0.3, -0.2, 0.1};
   Rot rot(D, 0.6, 0.35);
   Mesh m = MakeBox(D, D == 2 ? 10 : 5, W, W, rot, x0);
   std::vector<int> attrs;
   for (int i = 1; i <= m.bdr_attributes.Max(); i++) { attrs.push_back(i); }
   PmlBand band(m, attrs, Params(L));
   int nb = 0;
   double un = 0.0;
   const double res = PatchResidual(m, &band, 3, true, nb, un);
   std::printf("    %d of %d elements in a band, |u|_inf = %.3e\n", nb, m.GetNE(), un);
   Check(res < 1e-11, "band operator: |K0 u - b| / |K0 u| (exact solution)", res);
   int nb2 = 0;
   const double res_plain = PatchResidual(m, &band, 3, false, nb2, un);
   Check(res_plain > 1e-4, "plain operator does not reproduce it (the test has teeth)", res_plain);
}

static void TestBitIdentity(int D)
{
   std::printf("E %d-D: matrices outside the band, symmetry inside\n", D);
   const double W = 10.0, L = 1.5;
   const double x0[3] = {0, 0, 0};
   Rot rot(D, 0.25, 0.1);
   Mesh m = MakeBox(D, D == 2 ? 8 : 4, W, W, rot, x0);
   std::vector<int> attrs;
   for (int i = 1; i <= m.bdr_attributes.Max(); i++) { attrs.push_back(i); }
   PmlBand band(m, attrs, Params(L));
   ConstantCoefficient lam(3.2e10), mu(2.7e10);
   L2_FECollection fec(2, D);
   FiniteElementSpace fes(&m, &fec, D);
   PmlBandElasticityIntegrator bv(lam, mu, band);
   ElasticityIntegrator pv(lam, mu);
   double dmax = 0.0, asym = 0.0;
   int nout = 0, nin = 0;
   for (int e = 0; e < m.GetNE(); e++)
   {
      DenseMatrix A, B;
      ElementTransformation *T = m.GetElementTransformation(e);
      bv.AssembleElementMatrix(*fes.GetFE(e), *T, A);
      if (!PmlBandElement(m, e, band))
      {
         pv.AssembleElementMatrix(*fes.GetFE(e), *T, B);
         B -= A;
         dmax = std::max(dmax, B.MaxMaxNorm());
         nout++;
      }
      else
      {
         DenseMatrix At(A);
         At.Transpose();
         At -= A;
         asym = std::max(asym, At.MaxMaxNorm() / A.MaxMaxNorm());
         nin++;
      }
   }
   Check(nout > 0 && dmax == 0.0, "volume: element matrices outside the band identical", dmax);
   Check(nin > 0 && asym < 1e-13, "volume: band element matrices symmetric", asym);

   DGElasticityIPCombinedIntegrator fp(lam, mu, D, -1.0, 1.0), fb(lam, mu, D, -1.0, 1.0);
   fb.SetPmlBand(&band);
   double fdmax = 0.0, fasym = 0.0;
   int fout = 0, fin = 0;
   for (int f = 0; f < m.GetNumFaces(); f++)
   {
      FaceElementTransformations *FTr = m.GetInteriorFaceTransformations(f);
      if (!FTr) { continue; }
      DenseMatrix A, B;
      const FiniteElement &e1 = *fes.GetFE(FTr->Elem1No), &e2 = *fes.GetFE(FTr->Elem2No);
      fb.AssembleFaceMatrix(e1, e2, *FTr, A);
      fp.AssembleFaceMatrix(e1, e2, *FTr, B);
      DenseMatrix Dm(B);
      Dm -= A;
      if (Dm.MaxMaxNorm() == 0.0) { fout++; }
      else
      {
         fin++;
         DenseMatrix At(A);
         At.Transpose();
         At -= A;
         fasym = std::max(fasym, At.MaxMaxNorm());
      }
      fdmax = std::max(fdmax, 0.0);
   }
   Check(fout > 0 && fin > 0, "faces: some identical (outside), some changed (band)", fin);
   Check(fasym == 0.0, "faces: band face matrices symmetric", fasym);
}

// ---------------------------------------------------------------- T1
// Proposition 1 on the mesh: K0 on a box with bands against plain
// elasticity on the same mesh with its band vertices moved to x~(x).  Both
// get the same wall data, g(x~) of a smooth field.  In the inner domain the
// two solutions must agree up to the discretization error, which falls about
// 4x per halving of h (p = 1).
static void XTilde(const PmlBand &band, const double *x, double *xt)
{
   const int D = band.Dim();
   for (int i = 0; i < D; i++) { xt[i] = x[i]; }
   PmlPoint p = band.At(x);
   const double L = band.Params().thickness;
   const int N = band.Params().power;
   for (int k = 0; k < p.nact; k++)
   {
      const double Da = band.D0() * L / (N + 1) * std::pow(p.xi[k] / L, N + 1);
      for (int i = 0; i < D; i++) { xt[i] += Da / band.Alpha() * band.Walls()[p.wall[k]].n[i]; }
   }
}

struct WallData : public VectorCoefficient
{
   const PmlBand *band;   // non-null: evaluate at x~(x) (the band mesh)
   WallData(int D, const PmlBand *b) : VectorCoefficient(D), band(b) {}
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      const int D = vdim;
      Vector x(D);
      T.Transform(ip, x);
      double y[3] = {0, 0, 0};
      if (band) { XTilde(*band, x.GetData(), y); }
      else { for (int i = 0; i < D; i++) { y[i] = x(i); } }
      V.SetSize(D);
      const double k = 2.0 * M_PI / 40.0;   // smooth on the stretched box
      V(0) = 1e-3 * std::sin(k * y[0] + 0.3) * std::cos(0.7 * k * y[1]);
      V(1) = 1e-3 * std::cos(0.6 * k * y[0]) * std::sin(k * y[1] + 0.5);
   }
};

static void SolveWalls(Mesh &m, const PmlBand *band, bool band_ops, GridFunction &u)
{
   const int D = m.Dimension();
   ConstantCoefficient lam(3.2e10), mu(2.7e10);
   FiniteElementSpace &fes = *u.FESpace();
   Array<int> all(m.bdr_attributes.Max());
   all = 1;
   BilinearForm a(&fes);
   if (band_ops) { a.AddDomainIntegrator(new PmlBandElasticityIntegrator(lam, mu, *band)); }
   else { a.AddDomainIntegrator(new ElasticityIntegrator(lam, mu)); }
   auto *fi = new DGElasticityIPCombinedIntegrator(lam, mu, D, -1.0, 1.0);
   auto *bi = new DGElasticityIPCombinedIntegrator(lam, mu, D, -1.0, 1.0);
   if (band_ops) { fi->SetPmlBand(band); bi->SetPmlBand(band); }
   a.AddInteriorFaceIntegrator(fi);
   a.AddBdrFaceIntegrator(bi, all);
   a.Assemble(0);
   a.Finalize();
   WallData g(D, band_ops ? band : nullptr);
   Vector b(fes.GetVSize());
   b = 0.0;
   DGElasticityIPCombinedIntegrator rhs_int(lam, mu, D, -1.0, 1.0);
   if (band_ops) { rhs_int.SetPmlBand(band); }
   Array<int> vdofs;
   for (int be = 0; be < m.GetNBE(); be++)
   {
      FaceElementTransformations *FTr = m.GetBdrFaceTransformations(be);
      const FiniteElement *fe = fes.GetFE(FTr->Elem1No);
      const IntegrationRule &ir = IntRules.Get(FTr->GetGeometryType(), 2 * fe->GetOrder() + 1);
      const int nq = ir.GetNPoints();
      Vector gq(D * nq), gv(D);
      for (int q = 0; q < nq; q++)
      {
         FTr->SetAllIntPoints(&ir.IntPoint(q));
         g.Eval(gv, *FTr->Elem1, FTr->GetElement1IntPoint());
         for (int c = 0; c < D; c++) { gq(c * nq + q) = gv(c); }
      }
      Vector el;
      rhs_int.AssembleBoundaryFaceRHS(*fe, *FTr, gq, el);
      fes.GetElementVDofs(FTr->Elem1No, vdofs);
      b.AddElementVector(vdofs, el);
   }
   u = 0.0;
   GSSmoother M(a.SpMat());
   PCG(a.SpMat(), M, b, u, 0, 20000, 1e-26, 0.0);
}

static void TestMappedMesh()
{
   std::printf("T1 2-D: K0 against plain elasticity on the mapped mesh, p = 1\n");
   const double W = 10.0, L = 3.0;
   const double x0[3] = {0, 0, 0};
   Rot rot(2, 0.45, 0.0);
   double prev = 0.0;
   std::vector<double> errs;
   for (int n : {10, 20, 40})
   {
      Mesh mA = MakeBox(2, n, W, W, rot, x0);
      std::vector<int> attrs = {1, 2, 3, 4};
      PmlBand band(mA, attrs, Params(L));
      Mesh mB(mA);
      for (int v = 0; v < mB.GetNV(); v++)
      {
         double xt[3];
         XTilde(band, mA.GetVertex(v), xt);
         double *x = mB.GetVertex(v);
         x[0] = xt[0];
         x[1] = xt[1];
      }
      L2_FECollection fec(1, 2);
      FiniteElementSpace fA(&mA, &fec, 2), fB(&mB, &fec, 2);
      GridFunction uA(&fA), uB(&fB);
      SolveWalls(mA, &band, true, uA);
      SolveWalls(mB, &band, false, uB);
      double err = 0.0, un = 0.0;
      Array<int> vd;
      for (int e = 0; e < mA.GetNE(); e++)
      {
         if (PmlBandElement(mA, e, band)) { continue; }
         fA.GetElementVDofs(e, vd);
         for (int k = 0; k < vd.Size(); k++)
         {
            err = std::max(err, std::abs(uA(vd[k]) - uB(vd[k])));
            un = std::max(un, std::abs(uB(vd[k])));
         }
      }
      err /= un;
      std::printf("    n = %d: max inner-domain difference / max |u| = %.3e%s\n", n, err,
                  prev > 0 ? "" : "");
      if (prev > 0) { std::printf("      ratio %.2f\n", prev / err); }
      errs.push_back(err);
      prev = err;
   }
   const double r1 = errs[0] / errs[1], r2 = errs[1] / errs[2];
   Check(r1 > 3.0 && r2 > 3.0, "difference falls about 4x per halving of h (both ratios > 3)",
         std::min(r1, r2));
}

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   Mpi::Init(argc, argv);
#endif
   TestGeometry2D();
   TestGeometry3D();
   TestTensor(2);
   TestTensor(3);
   TestBitIdentity(2);
   TestBitIdentity(3);
   TestPatch(2);
   TestPatch(3);
   TestMappedMesh();
   std::printf("%s: %d failure(s)\n", g_fail ? "FAILED" : "ALL PASS", g_fail);
   return g_fail ? 1 : 0;
}
