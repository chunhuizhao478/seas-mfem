// Poroelastic plan rev 3, phase P2.  Step P2.1: DGPoroDivIntegrator
// (dg_poro_div_integrator.hpp), the divergence form B_c of derivation rev 1.2,
// eq. (6.3), with its fluid penalty (6.5), boundary rows (section 6.6) and
// source loads (6.6).  F-A (two tets sharing a face), medium A = tpv5-poro
// (element A), medium B = the V-F3 contrast set (element B), p = 1-3, the
// driver's GaussLobatto DG basis.
//
//   U2-1  Volume block (c = alpha M): symmetric; equal to a separate
//         quadrature loop written here (CalcPhysDShape, a 2p + 2 rule); zero
//         on rigid motions and divergence-free fields; b^T Z a =
//         c trA trB |K| for linear a = A x, b = B x.
//   U2-2  Face block (c = M, kappa_w = 1, so with the penalty): zero for
//         pairs of continuous fields; unchanged when the two elements swap
//         sides; for discontinuous linear fields equal to the closed-form
//         face integral; Penalty() equal to (6.5) by hand from |F| and |K|.
//   Extras, for the parts of P2.1 whose plan tests come later (U2-7, U2-9):
//   X2-1  Boundary rows on the six outer faces against closed forms:
//         kTestNormal and kTrialNormal (c = alpha M), kBothNormal with the
//         one-sided penalty (c = M); the kTrialNormal matrix is the
//         kTestNormal matrix transposed, bit for bit; the interior face
//         matrix and the kBothNormal matrices are symmetric bit for bit.
//   X2-2  Source load (6.6) through LinearForm: for a linear s and linear b
//         the closed form; with s = div a for a continuous a, the load
//         equals Z a (Z = volume + interior face blocks, no penalty).
//   X2-3  Boundary source load -int c s (b.n) against the closed form.
//   X2-4  ParBilinearForm (with KeepNbrBlock) and ParLinearForm on this rank
//         count against the serial assembly, row by row: on 2 ranks the
//         shared face is assembled once from each side.
//   X2-5  (Makefile) --abort interior-only | penalty-one-sided abort with
//         their messages.
//
// Tolerances (plan section 6.3): R(1e-14) for one evaluation and R(1e-12)
// for an assembled vector, relative to the largest term summed; B = bit for
// bit.  The forms b^T Z a of the extras X2-1 to X2-3 are measured against
// the sum of the |terms| b_i Z_ij a_j instead, chosen by Chunhui Zhao on
// 2026-10-21 in the review of step P2.1: on the first run X2-1 p=3 failed at
// 1.01e-14 against the largest term, with an error of 1.6e-16 times that sum
// (below one machine epsilon: round-off in the matrix entries of a 3600-term
// sum; summing in long double gave 1.00e-14).  Closed forms use the exact rules for linear functions on a triangle:
// int f = |F| (f1 + f2 + f3) / 3, int f g = |F| (sum f_v g_v + sum f sum g) / 12.
//
// Step P2.2: PoroNitscheOperator (poro_nitsche_operator.hpp), the bulk (u, w)
// operator of (8.1) and its force (9.2), on the TPV102 box: F-B (tpv5-poro
// everywhere) and F-C (tpv5-poro for x < 0, TPV205-equivalent for x > 0),
// p = 1-3.  The outer faces are drained walls (free top, Dirichlet sides and
// bottom), mixed walls (also Dirichlet + no-flow and absorbing + no-flow) or
// clamped (Dirichlet everywhere), see Walls.
//   U2-3  K symmetric for random x, y; K x on this rank count equal to K x on
//         one rank (MPI_COMM_SELF copy, DOFs matched by position).
//   U2-4  w = 0: the u rows of Force equal the Force of a NitscheWaveOperator
//         built from (lambda_u, mu, rho), bit for bit, with and without slip,
//         velocity and Kelvin-Voigt damping (absorbing faces present).
//   U2-5  Every eigenvalue of K at half the chosen kappa_w (kKappaW) by
//         LAPACK on the 2 x 2 x 1 box (1 rank; F-B and F-C drained, F-B
//         mixed): none below -1e-10 lambda_max, and the count within 1e-10
//         lambda_max of zero equal to the dimension of the divergence-free
//         BDM_p fields with w.n = 0 on no-flow faces (u = 0), the null space
//         of K.  Dense instead of the plan's LOBPCG, which cannot count that
//         null space, chosen by Chunhui Zhao on 2026-10-21 in the review of
//         step P2.2.
//   U2-6  PT-1 on F-C, clamped: the force of the interpolated exact fields
//         is zero; the revision 2.2 jump (alpha_A - alpha_B) u(0) against
//         poro_test_refs.hpp.
//   U2-7  (a) PT-2 (uniform V) likewise; (b) F_init of a quadratic p_init and
//         (c) F_V of a linear V on F-B with mixed walls against their
//         integrated-by-parts forms.
//   U2-8  Consistency on F-B, drained: for random polynomial (u, w) of degree
//         p with u = g on the Dirichlet faces, the force against the strong
//         form (-div sigma, -grad(M zeta)) and its boundary terms (sigma n on
//         the free surface, M zeta q.n on every drained face), integrated
//         against each test function with a 2p + 4 rule.
//   X2-6  (Makefile) --abort dirichlet-interior | no-flow-fault | kappa-zero
//         abort with their messages.
// U2-6 to U2-8 are measured against the largest of max |loads|, max |K x|
// (and max |strong form| in U2-8) over the rows of one field.
//
// Step P2.3: PoroSource (poro_source.hpp), the fluid source of (4.3) and its
// loads (6.6).
//   U2-9  V(t) = rate F(t) for step, ramp and smooth_step against
//         Gauss-Legendre integrals of the plan's f at six times; the
//         trapezoid V of the smooth step at T/2 converges at second order
//         (dt = T/8 to T/64, ratios in [3.5, 4.5]) and is exact for step and
//         ramp on a grid through T; a well inside one element (F-B, and F-C
//         on both sides of x = 0, p = 1-3) loads only that element's rows,
//         b^T L = M (alpha div u + div w)(x_w) for polynomial (u, w) of
//         degree p, and Force adds rate F(t) L bit for bit.  R(1e-15); the
//         well's b^T L against the sum of the |terms| (Chunhui Zhao,
//         2026-10-21, see TestU29).  (Makefile) --abort well-face |
//         well-vertex abort.
//   X2-7  (extra) a gaussian source through SetSource gives the F_V of
//         rate g(x) F(t) set directly (R(1e-12), mixed walls, p = 1, 2);
//         InitialPressure at quadrature points; (Makefile) --abort
//         well-outside | source-time (Force with the source behind the
//         operator's time) abort.
//
// Step P2.4: PoroMass (poro_mass.hpp), the inverse of M + dt D Gamma of (8.2),
// (8.3) and (9.5), on F-C (2 x 2 x 1), p = 1-3, dt from PoroMaxDt at CFL 0.5,
// gamma_w(z) = 1 / (1 - e^-z) - 1 / z with the drag b set so that z = dt / tau
// on every element.  Without absorbing faces (drained walls) and with them
// (mixed walls: x = +L absorbing).
//   U2-10 For z = 1e-2, 1, 2600, 2e5, Gamma_mom and Gamma_diag: PoroMass
//         applied to (M + dt D Gamma) x gives x, R(1e-12) per field;
//         M + dt D Gamma assembled independently (MFEM's
//         VectorMassIntegrator with a 2p + 2 rule for rho, rho_f and
//         m + dt gamma_w b; the elastic operator's dashpot C on the u rows).
//         Without absorbing faces the two Gamma give the same a, bit for bit
//         (Proposition 9.3).  (Makefile) --abort mass-gamma (gamma_w = 0.4)
//         aborts.
//   U2-11 rho_f = 0, z = 2600: the u rows of PoroMass equal NitscheMass(rho)
//         without absorbing faces and the elastic operator's damped mass
//         NitscheMass(rho, C, dt/2) with them, bit for bit, for both Gamma.
//
// Usage: mpirun -np 1|2|4 ./seas_test_poro_p2 [--only U2-3|...|U2-11|X2-7]
//        [--abort interior-only|penalty-one-sided|dirichlet-interior|
//                 no-flow-fault|kappa-zero|well-face|well-vertex|
//                 well-outside|source-time|mass-gamma]

#include "mfem.hpp"
#include "../dg_poro_div_integrator.hpp"
#include "../poro_cfl.hpp"
#include "../poro_mass.hpp"
#include "../poro_nitsche_operator.hpp"
#include "../poro_source.hpp"
#include "../../../friction/rate_state_closure.hpp"
#include "../../../config/tpv102_params.hpp"
#include "poro_fixtures.hpp"
#include "poro_test_refs.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

using mfem::seas::DGPoroDivIntegrator;
using mfem::seas::DGPoroDivSourceIntegrator;

namespace
{

const real_t kR1 = 1e-14;   // R, one evaluation
const real_t kRV = 1e-12;   // R, an assembled vector

// The two media of F-A, by element attribute - 1.
const real_t kAlphaM[2] = {poro_refs::tpv5_poro::alpha * poro_refs::tpv5_poro::M,
                           poro_refs::vf3_contrast::alpha * poro_refs::vf3_contrast::M};
const real_t kM[2] = {poro_refs::tpv5_poro::M, poro_refs::vf3_contrast::M};

std::string Fmt(real_t v)
{
   char b[32];
   std::snprintf(b, sizeof b, "%.2e", v);
   return b;
}

typedef std::function<void(const Vector &, Vector &)> VecFn;

// ------------------------------------------------------------- fields
// Fields in xi = x / kFAL, so every term is of order one.

/// a(x) = A xi + a0; div a = tr A / kFAL.
struct Lin
{
   real_t A[3][3];
   real_t a0[3];
   void Eval(const Vector &x, Vector &v) const
   {
      v.SetSize(3);
      for (int i = 0; i < 3; i++)
      {
         v(i) = a0[i];
         for (int j = 0; j < 3; j++) { v(i) += A[i][j] * x(j) / kFAL; }
      }
   }
   real_t Div() const { return (A[0][0] + A[1][1] + A[2][2]) / kFAL; }
   real_t Normal(const Vector &x, const Vector &n) const
   {
      Vector v; Eval(x, v);
      return v * n;
   }
   VecFn Fn() const { return [this](const Vector &x, Vector &v) { Eval(x, v); }; }
};

// a on element A / B and b on element A / B (attribute - 1).
const Lin kLinA[2] = {{{{0.31, -0.12, 0.05}, {0.08, 0.27, -0.19}, {-0.04, 0.15, 0.22}},
                       {0.11, -0.07, 0.03}},
                      {{{-0.18, 0.09, 0.14}, {0.21, 0.12, 0.06}, {0.03, -0.25, 0.33}},
                       {-0.05, 0.09, 0.12}}};
const Lin kLinB[2] = {{{{0.17, 0.04, -0.11}, {-0.09, 0.35, 0.07}, {0.12, 0.02, -0.06}},
                       {0.02, 0.08, -0.13}},
                      {{{0.26, -0.14, 0.03}, {0.05, -0.08, 0.19}, {-0.16, 0.07, 0.29}},
                       {0.07, -0.04, 0.06}}};

/// A vector polynomial of total degree deg with coefficients in [-1, 1] from
/// a fixed linear congruential sequence (the same on every machine).
struct Poly
{
   std::vector<std::array<int, 3>> mono;
   std::vector<real_t> co[3];

   Poly(int deg, uint32_t seed)
   {
      for (int i = 0; i <= deg; i++)
         for (int j = 0; i + j <= deg; j++)
            for (int k = 0; i + j + k <= deg; k++) { mono.push_back({i, j, k}); }
      uint32_t st = seed;
      for (int c = 0; c < 3; c++)
      {
         for (size_t m = 0; m < mono.size(); m++)
         {
            st = 1664525u * st + 1013904223u;
            co[c].push_back(2.0 * ((st >> 8) / 16777216.0) - 1.0);
         }
      }
   }
   static real_t Pow(real_t x, int e) { real_t r = 1.0; for (int i = 0; i < e; i++) { r *= x; } return r; }
   void Eval(const Vector &x, Vector &v) const
   {
      v.SetSize(3);
      v = 0.0;
      for (size_t m = 0; m < mono.size(); m++)
      {
         const real_t t = Pow(x(0) / kFAL, mono[m][0]) * Pow(x(1) / kFAL, mono[m][1])
                          * Pow(x(2) / kFAL, mono[m][2]);
         for (int c = 0; c < 3; c++) { v(c) += co[c][m] * t; }
      }
   }
   real_t Div(const Vector &x) const
   {
      real_t d = 0.0;
      for (size_t m = 0; m < mono.size(); m++)
      {
         for (int c = 0; c < 3; c++)
         {
            const int e = mono[m][c];
            if (e == 0) { continue; }
            real_t t = e / kFAL;
            for (int i = 0; i < 3; i++) { t *= Pow(x(i) / kFAL, i == c ? e - 1 : mono[m][i]); }
            d += co[c][m] * t;
         }
      }
      return d;
   }
   VecFn Fn() const { return [this](const Vector &x, Vector &v) { Eval(x, v); }; }
};

/// s(x) for the source loads: linear, of order one.
real_t SFn(const Vector &x)
{
   return 0.7 + 0.3 * x(0) / kFAL - 0.2 * x(1) / kFAL + 0.5 * x(2) / kFAL;
}

// ------------------------------------------------------------- geometry

Vector Vtx(const Mesh &m, int v)
{
   Vector x(3);
   for (int i = 0; i < 3; i++) { x(i) = m.GetVertex(v)[i]; }
   return x;
}

/// |K| from the vertex coordinates; also the centroid.
real_t TetVolume(const Mesh &m, int e, Vector &cen)
{
   Array<int> v;
   m.GetElementVertices(e, v);
   Vector x0 = Vtx(m, v[0]);
   DenseMatrix J(3);
   cen = x0;
   for (int k = 1; k < 4; k++)
   {
      Vector xk = Vtx(m, v[k]);
      cen += xk;
      for (int i = 0; i < 3; i++) { J(i, k - 1) = xk(i) - x0(i); }
   }
   cen *= 0.25;
   return std::abs(J.Det()) / 6.0;
}

/// A triangle: vertices, area |F| and the unit normal pointing away from
/// the point `inside`.
struct Tri
{
   Vector x[3];
   real_t area;
   Vector n;
};

Tri MakeTri(const Mesh &m, const Array<int> &fv, const Vector &inside)
{
   Tri t;
   for (int i = 0; i < 3; i++) { t.x[i] = Vtx(m, fv[i]); }
   Vector d1(t.x[1]), d2(t.x[2]), c(3);
   d1 -= t.x[0];
   d2 -= t.x[0];
   c(0) = d1(1) * d2(2) - d1(2) * d2(1);
   c(1) = d1(2) * d2(0) - d1(0) * d2(2);
   c(2) = d1(0) * d2(1) - d1(1) * d2(0);
   const real_t l = c.Norml2();
   t.area = 0.5 * l;
   t.n = c;
   t.n /= l;
   Vector out(t.x[0]);
   out -= inside;
   if (out * t.n < 0.0) { t.n.Neg(); }
   return t;
}

real_t Int1(const Tri &t, const real_t f[3]) { return t.area * (f[0] + f[1] + f[2]) / 3.0; }

real_t Int2(const Tri &t, const real_t f[3], const real_t g[3])
{
   const real_t ff = f[0] * g[0] + f[1] * g[1] + f[2] * g[2];
   return t.area * (ff + (f[0] + f[1] + f[2]) * (g[0] + g[1] + g[2])) / 12.0;
}

// ------------------------------------------------------------- algebra

/// The element vector of f at the element's nodes (MFEM's vector layout:
/// component c of node k at c * ndof + k).
Vector ProjectElem(const FiniteElement &fe, ElementTransformation &T, const VecFn &f)
{
   const IntegrationRule &nodes = fe.GetNodes();
   const int nd = fe.GetDof();
   Vector v(3 * nd), x(3), val(3);
   for (int k = 0; k < nd; k++)
   {
      T.Transform(nodes.IntPoint(k), x);
      f(x, val);
      for (int c = 0; c < 3; c++) { v(c * nd + k) = val(c); }
   }
   return v;
}

/// f on every element, chosen by the element's attribute, as a vector of the
/// space's local vector DOFs.
Vector ProjectGlobal(FiniteElementSpace &fes, const std::function<VecFn(int)> &by_attr)
{
   Mesh &m = *fes.GetMesh();
   Vector g(fes.GetVSize());
   Array<int> vd;
   IsoparametricTransformation T;
   for (int e = 0; e < m.GetNE(); e++)
   {
      m.GetElementTransformation(e, &T);
      const Vector v = ProjectElem(*fes.GetFE(e), T, by_attr(m.GetAttribute(e)));
      fes.GetElementVDofs(e, vd);
      g.SetSubVector(vd, v);
   }
   return g;
}

/// b^T Z a; big = the largest term |b_i Z_ij a_j|, and *sumabs (if given)
/// the sum of all |b_i Z_ij a_j|.
real_t Form(const DenseMatrix &Z, const Vector &a, const Vector &b, real_t &big,
            real_t *sumabs = nullptr)
{
   real_t s = 0.0, sa = 0.0;
   big = 0.0;
   for (int i = 0; i < Z.Height(); i++)
   {
      for (int j = 0; j < Z.Width(); j++)
      {
         const real_t t = b(i) * Z(i, j) * a(j);
         s += t;
         sa += std::abs(t);
         big = std::max(big, std::abs(t));
      }
   }
   if (sumabs) { *sumabs = sa; }
   return s;
}

/// Bit-for-bit symmetry: the number of entries with Z(i, j) != Z(j, i).
long AsymBits(const DenseMatrix &Z)
{
   long n = 0;
   for (int i = 0; i < Z.Height(); i++)
      for (int j = 0; j < Z.Width(); j++)
      {
         if (std::memcmp(&Z(i, j), &Z(j, i), sizeof(real_t)) != 0) { n++; }
      }
   return n;
}

/// max_i |(Z a)_i| / max_ij |Z_ij a_j|.
real_t NullR(const DenseMatrix &Z, const Vector &a)
{
   real_t num = 0.0, big = 0.0;
   for (int i = 0; i < Z.Height(); i++)
   {
      real_t s = 0.0;
      for (int j = 0; j < Z.Width(); j++)
      {
         s += Z(i, j) * a(j);
         big = std::max(big, std::abs(Z(i, j) * a(j)));
      }
      num = std::max(num, std::abs(s));
   }
   return num / big;
}

real_t MaxAbs(const DenseMatrix &Z) { return Z.MaxMaxNorm(); }

int InteriorFace(const Mesh &m)
{
   for (int f = 0; f < m.GetNumFaces(); f++)
   {
      if (m.FaceIsInterior(f)) { return f; }
   }
   return -1;
}

/// Local vector [f1 on side 1; f2 on side 2] of an interior face.
Vector Stack(FiniteElementSpace &fes, int e1, int e2, const VecFn &f1, const VecFn &f2)
{
   Mesh &m = *fes.GetMesh();
   IsoparametricTransformation T1, T2;
   m.GetElementTransformation(e1, &T1);
   m.GetElementTransformation(e2, &T2);
   const Vector v1 = ProjectElem(*fes.GetFE(e1), T1, f1);
   const Vector v2 = ProjectElem(*fes.GetFE(e2), T2, f2);
   Vector v(v1.Size() + v2.Size());
   for (int i = 0; i < v1.Size(); i++) { v(i) = v1(i); }
   for (int i = 0; i < v2.Size(); i++) { v(v1.Size() + i) = v2(i); }
   return v;
}

/// The test's own volume loop: CalcPhysDShape and a 2p + 2 rule.
DenseMatrix RefVolume(const FiniteElement &fe, ElementTransformation &T, real_t c)
{
   const int nd = fe.GetDof();
   DenseMatrix Z(3 * nd), g(nd, 3);
   Z = 0.0;
   const IntegrationRule &ir = IntRules.Get(Geometry::TETRAHEDRON, 2 * fe.GetOrder() + 2);
   for (int q = 0; q < ir.GetNPoints(); q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      T.SetIntPoint(&ip);
      fe.CalcPhysDShape(T, g);
      const real_t w = ip.weight * T.Weight() * c;
      for (int p = 0; p < 3; p++)
         for (int k = 0; k < nd; k++)
            for (int u = 0; u < 3; u++)
               for (int l = 0; l < nd; l++)
               {
                  Z(p * nd + k, u * nd + l) += w * g(k, p) * g(l, u);
               }
   }
   return Z;
}

/// sigma_w of (6.5) by hand on F-A: (D+1) = 4, p(p+D-1)/D = p(p+2)/3 and
/// D |J_F| / |J_K| = 3 (2 |F|) / (6 |K|) = |F| / |K|.
real_t SigmaInterior(real_t kappa, int p, real_t F, real_t K1, real_t M1,
                     real_t K2, real_t M2)
{
   return kappa / 4.0 * 4.0 * p * (p + 2.0) / 3.0 * F * (M1 / K1 + M2 / K2);
}

real_t SigmaBoundary(real_t kappa, int p, real_t F, real_t K, real_t M)
{
   return kappa * 4.0 * p * (p + 2.0) / 3.0 * F * M / K;
}

// ------------------------------------------------------------- U2-1

void TestU21()
{
   if (g_rank == 0) { std::cout << "--- U2-1 volume block of B_c on F-A (c = alpha M)\n"; }
   Mesh *mesh = TwoTets();
   for (int p = 1; p <= 3; p++)
   {
      DG_FECollection fec(p, 3, BasisType::GaussLobatto);
      FiniteElementSpace fes(mesh, &fec, 3, Ordering::byNODES);

      // Rigid motions and divergence-free fields of degree <= p.
      std::vector<VecFn> null_fields;
      for (int d = 0; d < 3; d++)
      {
         null_fields.push_back([d](const Vector &, Vector &v) { v.SetSize(3); v = 0.0; v(d) = 1.0; });
         null_fields.push_back([d](const Vector &x, Vector &v)
         {
            Vector w(3); w = 0.0; w(d) = 1.0;     // w x x
            v.SetSize(3);
            v(0) = (w(1) * x(2) - w(2) * x(1)) / kFAL;
            v(1) = (w(2) * x(0) - w(0) * x(2)) / kFAL;
            v(2) = (w(0) * x(1) - w(1) * x(0)) / kFAL;
         });
      }
      auto xi = [](const Vector &x, int i) { return x(i) / kFAL; };
      null_fields.push_back([&](const Vector &x, Vector &v)
      { v.SetSize(3); v(0) = xi(x, 1); v(1) = xi(x, 2); v(2) = xi(x, 0); });
      if (p >= 2)
      {
         null_fields.push_back([&](const Vector &x, Vector &v)
         { v.SetSize(3); v(0) = xi(x, 0) * xi(x, 0); v(1) = -2.0 * xi(x, 0) * xi(x, 1); v(2) = 0.0; });
         null_fields.push_back([&](const Vector &x, Vector &v)
         { v.SetSize(3); v(0) = xi(x, 1) * xi(x, 2); v(1) = xi(x, 0) * xi(x, 2); v(2) = xi(x, 0) * xi(x, 1); });
      }
      if (p >= 3)
      {
         null_fields.push_back([&](const Vector &x, Vector &v)
         {
            v.SetSize(3);
            v(0) = xi(x, 0) * xi(x, 0) * xi(x, 1);
            v(1) = -xi(x, 0) * xi(x, 1) * xi(x, 1);
            v(2) = 0.0;
         });
         null_fields.push_back([&](const Vector &x, Vector &v)
         {
            v.SetSize(3);
            v(0) = std::pow(xi(x, 1), 3); v(1) = std::pow(xi(x, 2), 3); v(2) = std::pow(xi(x, 0), 3);
         });
      }

      real_t r_sym = 0.0, r_ref = 0.0, r_null = 0.0, r_lin = 0.0;
      for (int e = 0; e < mesh->GetNE(); e++)
      {
         const int att = mesh->GetAttribute(e);
         const real_t cval = kAlphaM[att - 1];
         ConstantCoefficient c(cval);
         DGPoroDivIntegrator B(c);
         const FiniteElement &fe = *fes.GetFE(e);
         IsoparametricTransformation T;
         mesh->GetElementTransformation(e, &T);
         DenseMatrix Z;
         B.AssembleElementMatrix(fe, T, Z);

         real_t d = 0.0;
         for (int i = 0; i < Z.Height(); i++)
            for (int j = 0; j < Z.Width(); j++) { d = std::max(d, std::abs(Z(i, j) - Z(j, i))); }
         r_sym = std::max(r_sym, d / MaxAbs(Z));

         DenseMatrix Zr = RefVolume(fe, T, cval);
         DenseMatrix D(Z);
         D -= Zr;
         r_ref = std::max(r_ref, MaxAbs(D) / MaxAbs(Zr));

         for (const VecFn &f : null_fields)
         {
            r_null = std::max(r_null, NullR(Z, ProjectElem(fe, T, f)));
         }

         Lin a = kLinA[att - 1], b = kLinB[att - 1];
         for (int i = 0; i < 3; i++) { a.a0[i] = 0.0; b.a0[i] = 0.0; }   // a = A x, b = B x
         Vector cen;
         const real_t K = TetVolume(*mesh, e, cen);
         real_t big;
         const real_t val = Form(Z, ProjectElem(fe, T, a.Fn()), ProjectElem(fe, T, b.Fn()), big);
         const real_t ref = cval * a.Div() * b.Div() * K;
         r_lin = std::max(r_lin, std::abs(val - ref) / std::max(big, std::abs(ref)));
      }
      const std::string P = "U2-1 p=" + std::to_string(p) + " ";
      TEST_ASSERT(r_sym <= kR1, P + "volume block symmetric on both elements (max |Z - Z^T| / max |Z| = "
                  + Fmt(r_sym) + ")");
      TEST_ASSERT(r_ref <= kR1, P + "volume block equal to the test's own quadrature loop (R = "
                  + Fmt(r_ref) + ")");
      const size_t ndf = null_fields.size() - 6;
      TEST_ASSERT(r_null <= kR1, P + "Z a = 0 for 3 translations, 3 rotations and "
                  + std::to_string(ndf) + " divergence-free field" + (ndf == 1 ? "" : "s")
                  + " (R = " + Fmt(r_null) + ")");
      TEST_ASSERT(r_lin <= kR1, P + "b^T Z a = c trA trB |K| for a = A x, b = B x (R = "
                  + Fmt(r_lin) + ")");
   }
   delete mesh;
}

// ------------------------------------------------------------- U2-2

void TestU22()
{
   if (g_rank == 0) { std::cout << "--- U2-2 face block of B_c on F-A (c = M, kappa_w = 1)\n"; }
   Mesh *m1 = TwoTets(false), *m2 = TwoTets(true);
   Vector cvals(2);
   cvals(0) = kM[0];
   cvals(1) = kM[1];
   PWConstCoefficient cM(cvals);
   const real_t kappa = 1.0;
   DGPoroDivIntegrator B(cM, kappa);

   // Element volumes and centroids by attribute (A = 0, B = 1) on mesh 1.
   real_t K[2];
   Vector cen[2];
   for (int e = 0; e < 2; e++) { K[m1->GetAttribute(e) - 1] = TetVolume(*m1, e, cen[m1->GetAttribute(e) - 1]); }

   const int f1 = InteriorFace(*m1), f2 = InteriorFace(*m2);
   Array<int> fv;
   m1->GetFaceVertices(f1, fv);
   const Tri tri = MakeTri(*m1, fv, cen[0]);     // n out of element A

   for (int p = 1; p <= 3; p++)
   {
      DG_FECollection fec(p, 3, BasisType::GaussLobatto);
      FiniteElementSpace fes1(m1, &fec, 3, Ordering::byNODES);
      FiniteElementSpace fes2(m2, &fec, 3, Ordering::byNODES);
      const std::string P = "U2-2 p=" + std::to_string(p) + " ";

      // Face matrices on both meshes.
      FaceElementTransformations *T1 = m1->GetInteriorFaceTransformations(f1);
      const int a1 = T1->Elem1No, b1 = T1->Elem2No;
      DenseMatrix Z1;
      B.AssembleFaceMatrix(*fes1.GetFE(a1), *fes1.GetFE(b1), *T1, Z1);
      real_t nl_q, detJ1, detJ2;
      {
         const IntegrationPoint &ip = IntRules.Get(T1->GetGeometryType(), 1).IntPoint(0);
         T1->SetAllIntPoints(&ip);
         Vector nor(3);
         CalcOrtho(T1->Jacobian(), nor);
         nl_q = nor.Norml2();
         detJ1 = T1->Elem1->Weight();
         detJ2 = T1->Elem2->Weight();
      }
      const int att1 = m1->GetAttribute(a1), att2 = m1->GetAttribute(b1);

      FaceElementTransformations *T2 = m2->GetInteriorFaceTransformations(f2);
      const int a2 = T2->Elem1No, b2 = T2->Elem2No;
      DenseMatrix Z2;
      B.AssembleFaceMatrix(*fes2.GetFE(a2), *fes2.GetFE(b2), *T2, Z2);
      const bool swapped = (m2->GetAttribute(a2) != att1);

      // (a) continuous fields: every face term has a jump factor.
      real_t r_cont = 0.0;
      for (int k = 0; k < 3; k++)
      {
         const Poly a(p, 11u + 7u * k), b(p, 101u + 13u * k);
         real_t big;
         const real_t val = Form(Z1, Stack(fes1, a1, b1, a.Fn(), a.Fn()),
                                 Stack(fes1, a1, b1, b.Fn(), b.Fn()), big);
         r_cont = std::max(r_cont, std::abs(val) / big);
      }
      TEST_ASSERT(r_cont <= kR1, P + "b^T Z_F a = 0 for 3 pairs of continuous degree-p fields (R = "
                  + Fmt(r_cont) + ")");

      // (b) swap: the two face matrices, put in mesh 1's global numbering,
      // match entry by entry.  Each element keeps its vertex list, so the
      // same local DOF is the same node (checked).
      Array<int> vd, vd2;
      std::vector<int> perm(fes1.GetVSize(), -1);
      real_t node_gap = 0.0;
      for (int e = 0; e < 2; e++)
      {
         const int att = m1->GetAttribute(e);
         int e2 = (m2->GetAttribute(0) == att) ? 0 : 1;
         fes1.GetElementVDofs(e, vd);
         fes2.GetElementVDofs(e2, vd2);
         for (int i = 0; i < vd.Size(); i++) { perm[vd[i]] = vd2[i]; }
         IsoparametricTransformation Ta, Tb;
         m1->GetElementTransformation(e, &Ta);
         m2->GetElementTransformation(e2, &Tb);
         const IntegrationRule &nodes = fes1.GetFE(e)->GetNodes();
         for (int k = 0; k < nodes.GetNPoints(); k++)
         {
            Vector xa(3), xb(3);
            Ta.Transform(nodes.IntPoint(k), xa);
            Tb.Transform(nodes.IntPoint(k), xb);
            xa -= xb;
            node_gap = std::max(node_gap, xa.Normlinf());
         }
      }
      const int n = fes1.GetVSize();
      DenseMatrix G1(n), G2(n);
      G1 = 0.0;
      G2 = 0.0;
      {
         Array<int> v1, v2, t;
         fes1.GetElementVDofs(a1, v1); fes1.GetElementVDofs(b1, t); v1.Append(t);
         fes2.GetElementVDofs(a2, v2); fes2.GetElementVDofs(b2, t); v2.Append(t);
         G1.AddSubMatrix(v1, Z1);
         G2.AddSubMatrix(v2, Z2);
      }
      real_t d = 0.0;
      for (int i = 0; i < n; i++)
         for (int j = 0; j < n; j++) { d = std::max(d, std::abs(G1(i, j) - G2(perm[i], perm[j]))); }
      const real_t r_swap = d / MaxAbs(G1);
      TEST_ASSERT(swapped && node_gap == 0.0 && r_swap <= kR1,
                  P + "face matrix unchanged when the elements swap sides (sides swapped: "
                  + (swapped ? "yes" : "NO") + ", node gap " + Fmt(node_gap) + " m, R = "
                  + Fmt(r_swap) + ")");

      // (c) discontinuous linear fields against the closed form.
      const VecFn aA = kLinA[0].Fn(), aB = kLinA[1].Fn(), bA = kLinB[0].Fn(), bB = kLinB[1].Fn();
      const Vector av = Stack(fes1, a1, b1, att1 == 1 ? aA : aB, att2 == 1 ? aA : aB);
      const Vector bv = Stack(fes1, a1, b1, att1 == 1 ? bA : bB, att2 == 1 ? bA : bB);
      real_t big;
      const real_t val = Form(Z1, av, bv, big);
      real_t f[3], g[3];
      for (int v = 0; v < 3; v++)
      {
         f[v] = kLinA[0].Normal(tri.x[v], tri.n) - kLinA[1].Normal(tri.x[v], tri.n);   // [a.n]
         g[v] = kLinB[0].Normal(tri.x[v], tri.n) - kLinB[1].Normal(tri.x[v], tri.n);   // [b.n]
      }
      const real_t sw = SigmaInterior(kappa, p, tri.area, K[0], kM[0], K[1], kM[1]);
      const real_t avg_div_a = 0.5 * (kM[0] * kLinA[0].Div() + kM[1] * kLinA[1].Div());
      const real_t avg_div_b = 0.5 * (kM[0] * kLinB[0].Div() + kM[1] * kLinB[1].Div());
      const real_t ref = -avg_div_a * Int1(tri, g) - avg_div_b * Int1(tri, f) + sw * Int2(tri, f, g);
      const real_t r_cf = std::abs(val - ref) / std::max(big, std::abs(ref));
      TEST_ASSERT(r_cf <= kR1, P + "b^T Z_F a for discontinuous linear a, b equal to the closed-form "
                  "face integral (R = " + Fmt(r_cf) + ")");

      // (d) Penalty() against (6.5) by hand, interior and boundary.
      const real_t pen = B.Penalty(3, p, nl_q, detJ1, kM[att1 - 1], detJ2, kM[att2 - 1], true);
      const real_t r_pen = std::abs(pen - sw) / sw;
      FaceElementTransformations *Tb = m1->GetBdrFaceTransformations(0);
      const int eb = Tb->Elem1No;
      Vector cb;
      const real_t Kb = TetVolume(*m1, eb, cb);
      Array<int> bv_;
      m1->GetBdrElementVertices(0, bv_);
      const Tri tb = MakeTri(*m1, bv_, cb);
      real_t nlb, detb;
      {
         const IntegrationPoint &ip = IntRules.Get(Tb->GetGeometryType(), 1).IntPoint(0);
         Tb->SetAllIntPoints(&ip);
         Vector nor(3);
         CalcOrtho(Tb->Jacobian(), nor);
         nlb = nor.Norml2();
         detb = Tb->Elem1->Weight();
      }
      const real_t Mb = kM[m1->GetAttribute(eb) - 1];
      const real_t penb = B.Penalty(3, p, nlb, detb, Mb, 0.0, 0.0, false);
      const real_t swb = SigmaBoundary(kappa, p, tb.area, Kb, Mb);
      const real_t r_penb = std::abs(penb - swb) / swb;
      TEST_ASSERT(r_pen <= kR1 && r_penb <= kR1,
                  P + "sigma_w equal to (6.5) by hand from |F| and |K|: interior "
                  + Fmt(pen) + " Pa/m (R = " + Fmt(r_pen) + "), boundary " + Fmt(penb)
                  + " Pa/m (R = " + Fmt(r_penb) + ")");
   }
   delete m1;
   delete m2;
}

// ------------------------------------------------------------- X2-1

void TestX21()
{
   if (g_rank == 0) { std::cout << "--- X2-1 boundary rows (derivation section 6.6) on F-A\n"; }
   Mesh *mesh = TwoTets();
   Vector va(2), vm(2);
   va(0) = kAlphaM[0]; va(1) = kAlphaM[1];
   vm(0) = kM[0]; vm(1) = kM[1];
   PWConstCoefficient cA(va), cM(vm);
   const real_t kappa = 1.0;
   DGPoroDivIntegrator Btn(cA, 0.0, DGPoroDivIntegrator::kTestNormal);
   DGPoroDivIntegrator Btr(cA, 0.0, DGPoroDivIntegrator::kTrialNormal);
   DGPoroDivIntegrator Bbo(cM, kappa, DGPoroDivIntegrator::kBothNormal);
   for (int p = 1; p <= 3; p++)
   {
      DG_FECollection fec(p, 3, BasisType::GaussLobatto);
      FiniteElementSpace fes(mesh, &fec, 3, Ordering::byNODES);
      real_t r_tn = 0.0, r_tr = 0.0, r_bo = 0.0;
      long nbits = 0, nsym = 0;
      for (int be = 0; be < mesh->GetNBE(); be++)
      {
         FaceElementTransformations *T = mesh->GetBdrFaceTransformations(be);
         const int e = T->Elem1No;
         const int att = mesh->GetAttribute(e);
         const FiniteElement &fe = *fes.GetFE(e);
         DenseMatrix Ztn, Ztr, Zbo;
         Btn.AssembleFaceMatrix(fe, fe, *T, Ztn);
         Btr.AssembleFaceMatrix(fe, fe, *T, Ztr);
         Bbo.AssembleFaceMatrix(fe, fe, *T, Zbo);
         nsym += AsymBits(Zbo);
         for (int i = 0; i < Ztn.Height(); i++)
            for (int j = 0; j < Ztn.Width(); j++)
            {
               if (std::memcmp(&Ztr(i, j), &Ztn(j, i), sizeof(real_t)) != 0) { nbits++; }
            }

         IsoparametricTransformation Te;
         mesh->GetElementTransformation(e, &Te);
         const Lin &a = kLinA[att - 1], &b = kLinB[att - 1];
         const Vector av = ProjectElem(fe, Te, a.Fn()), bv = ProjectElem(fe, Te, b.Fn());
         Vector cen;
         const real_t K = TetVolume(*mesh, e, cen);
         Array<int> fv;
         mesh->GetBdrElementVertices(be, fv);
         const Tri t = MakeTri(*mesh, fv, cen);
         real_t f[3], g[3];
         for (int v = 0; v < 3; v++) { f[v] = a.Normal(t.x[v], t.n); g[v] = b.Normal(t.x[v], t.n); }

         const real_t tn = -kAlphaM[att - 1] * a.Div() * Int1(t, g);
         const real_t tr = -kAlphaM[att - 1] * b.Div() * Int1(t, f);
         const real_t bo = -kM[att - 1] * a.Div() * Int1(t, g) - kM[att - 1] * b.Div() * Int1(t, f)
                           + SigmaBoundary(kappa, p, t.area, K, kM[att - 1]) * Int2(t, f, g);
         real_t big, sa, val;
         val = Form(Ztn, av, bv, big, &sa);
         r_tn = std::max(r_tn, std::abs(val - tn) / std::max(sa, std::abs(tn)));
         val = Form(Ztr, av, bv, big, &sa);
         r_tr = std::max(r_tr, std::abs(val - tr) / std::max(sa, std::abs(tr)));
         val = Form(Zbo, av, bv, big, &sa);
         r_bo = std::max(r_bo, std::abs(val - bo) / std::max(sa, std::abs(bo)));
      }
      {
         const int f = InteriorFace(*mesh);
         FaceElementTransformations *T = mesh->GetInteriorFaceTransformations(f);
         DGPoroDivIntegrator Bin(cM, kappa);
         DenseMatrix Zin;
         Bin.AssembleFaceMatrix(*fes.GetFE(T->Elem1No), *fes.GetFE(T->Elem2No), *T, Zin);
         nsym += AsymBits(Zin);
      }
      const std::string P = "X2-1 p=" + std::to_string(p) + " (extra) ";
      TEST_ASSERT(r_tn <= kR1, P + "kTestNormal: b^T Z a = -int c (div a)(b.n) on the 6 outer faces (R = "
                  + Fmt(r_tn) + " of the sum of |terms|)");
      TEST_ASSERT(r_tr <= kR1, P + "kTrialNormal: b^T Z a = -int c (div b)(a.n) (R = " + Fmt(r_tr) + " of the sum of |terms|)");
      TEST_ASSERT(r_bo <= kR1, P + "kBothNormal with the one-sided penalty against its closed form (R = "
                  + Fmt(r_bo) + " of the sum of |terms|)");
      TEST_ASSERT(nbits == 0, P + "kTrialNormal matrix = kTestNormal matrix transposed, bit for bit ("
                  + std::to_string(nbits) + " entries differ)");
      TEST_ASSERT(nsym == 0, P + "interior face matrix (kappa_w = 1) and the 6 kBothNormal matrices "
                  "symmetric bit for bit (" + std::to_string(nsym) + " entries differ)");
   }
   delete mesh;
}

// ------------------------------------------------------------- X2-2, X2-3

void TestX22()
{
   if (g_rank == 0) { std::cout << "--- X2-2, X2-3 source loads of (6.6) on F-A\n"; }
   Mesh *mesh = TwoTets();
   Vector va(2), vm(2);
   va(0) = kAlphaM[0]; va(1) = kAlphaM[1];
   vm(0) = kM[0]; vm(1) = kM[1];
   PWConstCoefficient cA(va), cM(vm);
   FunctionCoefficient s(SFn);
   DGPoroDivIntegrator BA(cA), BM(cM);
   Array<int> all_bdr(mesh->bdr_attributes.Max());
   all_bdr = 1;

   real_t K[2];
   Vector cen[2];
   for (int e = 0; e < 2; e++) { K[mesh->GetAttribute(e) - 1] = TetVolume(*mesh, e, cen[mesh->GetAttribute(e) - 1]); }
   Array<int> fv;
   mesh->GetFaceVertices(InteriorFace(*mesh), fv);
   const Tri tri = MakeTri(*mesh, fv, cen[0]);

   for (int p = 1; p <= 3; p++)
   {
      DG_FECollection fec(p, 3, BasisType::GaussLobatto);
      FiniteElementSpace fes(mesh, &fec, 3, Ordering::byNODES);
      const std::string P = "X2-2 p=" + std::to_string(p) + " (extra) ";
      const Vector bv = ProjectGlobal(fes, [](int att) { return kLinB[att - 1].Fn(); });

      // (a) linear s, discontinuous linear b: closed form.
      LinearForm Fe(&fes), Ff(&fes);
      Fe.AddDomainIntegrator(new DGPoroDivSourceIntegrator(BA, s));
      Ff.AddInteriorFaceIntegrator(new DGPoroDivSourceIntegrator(BA, s));
      Fe.Assemble();
      Ff.Assemble();
      real_t val = 0.0, big = 0.0;      // big = the sum of the |terms|
      for (int i = 0; i < bv.Size(); i++)
      {
         val += bv(i) * Fe[i] + bv(i) * Ff[i];
         big += std::abs(bv(i) * Fe[i]) + std::abs(bv(i) * Ff[i]);
      }
      real_t sv[3], g[3];
      for (int v = 0; v < 3; v++)
      {
         sv[v] = SFn(tri.x[v]);
         g[v] = kLinB[0].Normal(tri.x[v], tri.n) - kLinB[1].Normal(tri.x[v], tri.n);
      }
      const real_t ref = kAlphaM[0] * kLinB[0].Div() * K[0] * SFn(cen[0])
                         + kAlphaM[1] * kLinB[1].Div() * K[1] * SFn(cen[1])
                         - 0.5 * (kAlphaM[0] + kAlphaM[1]) * Int2(tri, sv, g);
      const real_t r_a = std::abs(val - ref) / std::max(big, std::abs(ref));
      TEST_ASSERT(r_a <= kR1, P + "b . F for linear s and discontinuous linear b equal to "
                  "sum_K c trB int_K s - int_F {c s}[b.n] (R = " + Fmt(r_a) + " of the sum of |terms|)");

      // (b) s = div a for a continuous a: F = Z a.
      const Poly a(p, 7u);
      FunctionCoefficient sdiv([&a](const Vector &x) { return a.Div(x); });
      LinearForm F(&fes);
      F.AddDomainIntegrator(new DGPoroDivSourceIntegrator(BA, sdiv));
      F.AddInteriorFaceIntegrator(new DGPoroDivSourceIntegrator(BA, sdiv));
      F.Assemble();
      BilinearForm Z(&fes);
      Z.AddDomainIntegrator(new DGPoroDivIntegrator(cA));
      Z.AddInteriorFaceIntegrator(new DGPoroDivIntegrator(cA));
      Z.Assemble(0);
      Z.Finalize(0);
      const Vector av = ProjectGlobal(fes, [&a](int) { return a.Fn(); });
      DenseMatrix Zd;
      Z.SpMat().ToDenseMatrix(Zd);
      Vector Za(av.Size());
      Zd.Mult(av, Za);
      real_t num = 0.0, bigz = 0.0;
      for (int i = 0; i < Zd.Height(); i++)
      {
         num = std::max(num, std::abs(F[i] - Za(i)));
         for (int j = 0; j < Zd.Width(); j++) { bigz = std::max(bigz, std::abs(Zd(i, j) * av(j))); }
      }
      TEST_ASSERT(num / bigz <= kRV, P + "with s = div a for a continuous degree-p a, the load equals Z a "
                  "(R = " + Fmt(num / bigz) + ", assembled vector)");

      // X2-3: the boundary load -int c s (b.n) on the six outer faces.
      LinearForm Fb(&fes);
      Fb.AddBdrFaceIntegrator(new DGPoroDivSourceIntegrator(BM, s), all_bdr);
      Fb.Assemble();
      real_t valb = 0.0, bigb = 0.0, refb = 0.0;
      for (int i = 0; i < bv.Size(); i++)
      {
         valb += bv(i) * Fb[i];
         bigb += std::abs(bv(i) * Fb[i]);
      }
      for (int be = 0; be < mesh->GetNBE(); be++)
      {
         FaceElementTransformations *T = mesh->GetBdrFaceTransformations(be);
         const int e = T->Elem1No, att = mesh->GetAttribute(e);
         Vector ce;
         TetVolume(*mesh, e, ce);
         Array<int> bfv;
         mesh->GetBdrElementVertices(be, bfv);
         const Tri t = MakeTri(*mesh, bfv, ce);
         real_t s3[3], g3[3];
         for (int v = 0; v < 3; v++) { s3[v] = SFn(t.x[v]); g3[v] = kLinB[att - 1].Normal(t.x[v], t.n); }
         refb += -kM[att - 1] * Int2(t, s3, g3);
      }
      const real_t r_b = std::abs(valb - refb) / std::max(bigb, std::abs(refb));
      TEST_ASSERT(r_b <= kR1, "X2-3 p=" + std::to_string(p) + " (extra) boundary load b . F = "
                  "-sum_F int c s (b.n) on the 6 outer faces (R = " + Fmt(r_b) + " of the sum of |terms|)");
   }
   delete mesh;
}

// ------------------------------------------------------------- X2-4

/// For each local vector DOF of pfes, the serial vector DOF at the same node
/// and component (elements matched by centroid, nodes by position).
std::vector<int> MatchDofs(ParFiniteElementSpace &pfes, FiniteElementSpace &fes, real_t &gap)
{
   ParMesh &pm = *pfes.GetParMesh();
   Mesh &m = *fes.GetMesh();
   std::vector<int> map(pfes.GetVSize(), -1);
   gap = 0.0;
   Array<int> pvd, svd;
   for (int e = 0; e < pm.GetNE(); e++)
   {
      IsoparametricTransformation Tp;
      pm.GetElementTransformation(e, &Tp);
      Vector cp(3);
      Tp.Transform(Geometries.GetCenter(Geometry::TETRAHEDRON), cp);
      int E = -1;
      for (int k = 0; k < m.GetNE(); k++)
      {
         IsoparametricTransformation Ts;
         m.GetElementTransformation(k, &Ts);
         Vector cs(3);
         Ts.Transform(Geometries.GetCenter(Geometry::TETRAHEDRON), cs);
         cs -= cp;
         if (cs.Normlinf() < 1e-9 * kFAL) { E = k; }
      }
      if (E < 0) { gap = 1e30; return map; }
      IsoparametricTransformation Ts;
      m.GetElementTransformation(E, &Ts);
      const FiniteElement &fe = *pfes.GetFE(e);
      const IntegrationRule &nodes = fe.GetNodes();
      const int nd = fe.GetDof();
      pfes.GetElementVDofs(e, pvd);
      fes.GetElementVDofs(E, svd);
      for (int k = 0; k < nd; k++)
      {
         Vector xp(3);
         Tp.Transform(nodes.IntPoint(k), xp);
         int K = -1;
         real_t best = 1e30;
         for (int l = 0; l < nd; l++)
         {
            Vector xs(3);
            Ts.Transform(nodes.IntPoint(l), xs);
            xs -= xp;
            if (xs.Normlinf() < best) { best = xs.Normlinf(); K = l; }
         }
         gap = std::max(gap, best);
         for (int c = 0; c < 3; c++) { map[pvd[c * nd + k]] = svd[c * nd + K]; }
      }
   }
   return map;
}

void TestX24()
{
   const int np = Mpi::WorldSize();
   if (g_rank == 0) { std::cout << "--- X2-4 parallel assembly on " << np << " rank(s) against the serial one\n"; }
   if (np > 2)
   {
      // F-A has two elements; make test-poro-p2 runs X2-4 on 1 and 2 ranks.
      if (g_rank == 0) { std::cout << "  X2-4 runs on 1 or 2 ranks only (F-A has two elements)\n"; }
      return;
   }
   Mesh *mesh = TwoTets();
   int part[2] = {0, np > 1 ? 1 : 0};
   ParMesh pm(MPI_COMM_WORLD, *mesh, part);
   Vector vm(2);
   vm(0) = kM[0]; vm(1) = kM[1];
   PWConstCoefficient cM(vm);
   FunctionCoefficient s(SFn);
   const real_t kappa = 1.0;
   Array<int> all_bdr(mesh->bdr_attributes.Max());
   all_bdr = 1;
   DGPoroDivIntegrator Bs(cM, kappa);
   for (int p = 1; p <= 3; p++)
   {
      DG_FECollection fec(p, 3, BasisType::GaussLobatto);
      FiniteElementSpace fes(mesh, &fec, 3, Ordering::byNODES);
      ParFiniteElementSpace pfes(&pm, &fec, 3, Ordering::byNODES);
      const std::string P = "X2-4 p=" + std::to_string(p) + " (extra) ";

      BilinearForm K(&fes);
      K.AddDomainIntegrator(new DGPoroDivIntegrator(cM, kappa));
      K.AddInteriorFaceIntegrator(new DGPoroDivIntegrator(cM, kappa));
      K.AddBdrFaceIntegrator(new DGPoroDivIntegrator(cM, kappa, DGPoroDivIntegrator::kBothNormal), all_bdr);
      K.Assemble(0);
      K.Finalize(0);

      ParBilinearForm Kp(&pfes);
      Kp.KeepNbrBlock(true);
      Kp.AddDomainIntegrator(new DGPoroDivIntegrator(cM, kappa));
      Kp.AddInteriorFaceIntegrator(new DGPoroDivIntegrator(cM, kappa));
      Kp.AddBdrFaceIntegrator(new DGPoroDivIntegrator(cM, kappa, DGPoroDivIntegrator::kBothNormal), all_bdr);
      Kp.Assemble(0);
      Kp.Finalize(0);
      HypreParMatrix *Kh = Kp.ParallelAssemble();

      LinearForm F(&fes);
      F.AddDomainIntegrator(new DGPoroDivSourceIntegrator(Bs, s));
      F.AddInteriorFaceIntegrator(new DGPoroDivSourceIntegrator(Bs, s));
      F.AddBdrFaceIntegrator(new DGPoroDivSourceIntegrator(Bs, s), all_bdr);
      F.Assemble();
      ParLinearForm Fp(&pfes);
      Fp.AddDomainIntegrator(new DGPoroDivSourceIntegrator(Bs, s));
      Fp.AddInteriorFaceIntegrator(new DGPoroDivSourceIntegrator(Bs, s));
      Fp.AddBdrFaceIntegrator(new DGPoroDivSourceIntegrator(Bs, s), all_bdr);
      Fp.Assemble();

      const Poly xa(p, 3u), xb(p, 5u);
      auto by_attr = [&](int att) { return att == 1 ? xa.Fn() : xb.Fn(); };
      const Vector xs = ProjectGlobal(fes, by_attr);
      const Vector xp = ProjectGlobal(pfes, by_attr);
      Vector ys(xs.Size()), yp(xp.Size());
      K.Mult(xs, ys);
      Kh->Mult(xp, yp);

      DenseMatrix Kd;
      K.SpMat().ToDenseMatrix(Kd);
      real_t big = 0.0;
      for (int i = 0; i < Kd.Height(); i++)
         for (int j = 0; j < Kd.Width(); j++) { big = std::max(big, std::abs(Kd(i, j) * xs(j))); }

      real_t gap;
      const std::vector<int> map = MatchDofs(pfes, fes, gap);
      real_t dk = 0.0, df = 0.0;
      for (int i = 0; i < yp.Size(); i++)
      {
         dk = std::max(dk, std::abs(yp(i) - ys(map[i])));
         df = std::max(df, std::abs(Fp[i] - F[map[i]]));
      }
      dk = RMax(dk, MPI_COMM_WORLD);
      df = RMax(df, MPI_COMM_WORLD);
      gap = RMax(gap, MPI_COMM_WORLD);
      const real_t r_k = dk / big, r_f = df / F.Normlinf();
      TEST_ASSERT(gap < 1e-9 * kFAL && r_k <= kRV,
                  P + "ParBilinearForm K x (interior, shared and no-flow faces, kappa_w = 1) equals the "
                  "serial K x on every local row (R = " + Fmt(r_k) + ", assembled vector)");
      TEST_ASSERT(gap < 1e-9 * kFAL && r_f <= kRV,
                  P + "ParLinearForm source load equals the serial load on every local row (R = "
                  + Fmt(r_f) + ", assembled vector)");
      delete Kh;
   }
   delete mesh;
}

// ------------------------------------------------------------- X2-5

int RunAbort(const std::string &mode)
{
   Mesh *mesh = TwoTets();
   Vector vm(2);
   vm(0) = kM[0]; vm(1) = kM[1];
   PWConstCoefficient cM(vm);
   if (g_rank == 0) { std::cout << "X2-5 expecting an abort: " << mode << "\n"; }
   if (mode == "interior-only")
   {
      DG_FECollection fec(1, 3, BasisType::GaussLobatto);
      FiniteElementSpace fes(mesh, &fec, 3, Ordering::byNODES);
      Array<int> all_bdr(mesh->bdr_attributes.Max());
      all_bdr = 1;
      BilinearForm K(&fes);
      K.AddBdrFaceIntegrator(new DGPoroDivIntegrator(cM), all_bdr);
      K.Assemble(0);
   }
   else if (mode == "penalty-one-sided")
   {
      DGPoroDivIntegrator B(cM, 1.0, DGPoroDivIntegrator::kTestNormal);
   }
   if (g_rank == 0) { std::cout << "FAILED: X2-5 " << mode << " did not abort\n"; }
   delete mesh;
   return 1;
}

// ================================================================== P2.2
// Step P2.2: PoroNitscheOperator (poro_nitsche_operator.hpp) on F-B
// (tpv5-poro everywhere) and F-C (tpv5-poro for x < 0, the TPV205-equivalent
// set for x > 0), the TPV102 box with its fault on Y = 0.

using mfem::seas::BoundaryConfig;
using mfem::seas::DieterichRuinaFriction;
using mfem::seas::DirichletFunc;
using mfem::seas::FaultPointState;
using mfem::seas::IFaultClosure;
using mfem::seas::MaterialField;
using mfem::seas::NitscheFaultSetup;
using mfem::seas::NitscheWaveOperator;
using mfem::seas::RateStateClosure;
using mfem::seas::TPV102Params;
using mfem::seas::poro::BiotConstants;
using mfem::seas::poro::BiotFluidInputs;
using mfem::seas::poro::BulkDensity;
using mfem::seas::poro::DeriveBiotConstants;
using mfem::seas::poro::PoroBoundaryConfig;
using mfem::seas::poro::PoroMaterialSample;
using mfem::seas::poro::PoroNitscheOperator;

struct Medium
{
   BiotFluidInputs f;
   real_t lam, mu;
};

#define PORO_MEDIUM(NS) Medium{BiotFluidInputs{NS::phi, NS::kappa, NS::eta, NS::T, \
   NS::rho_f, NS::rho_s, NS::Kf, NS::Ks}, NS::lam, NS::mu}

const Medium kMedA = PORO_MEDIUM(poro_refs::tpv5_poro);           // F-B; F-C x < 0
const Medium kMedB = PORO_MEDIUM(poro_refs::tpv205_equivalent);   // F-C x > 0

/// The outer faces (the fault on Y = 0 in every case):
///   kDrained  z = 0 free surface (attribute 1), every other face Dirichlet (5);
///   kMixed    z = 0 free; x = -L Dirichlet and no-flow (6); x = +L absorbing
///             and no-flow (7); y = +-L and z = -L Dirichlet (5);
///   kClamped  every outer face Dirichlet (1 and 5).
/// All of them are drained except the no-flow faces.
enum class Walls { kDrained, kMixed, kClamped };

Mesh *BoxMesh(int nx, int ny, int nz, Walls walls)
{
   Mesh *m = TPV102Box(nx, ny, nz);
   MarkSidesFC(*m);
   if (walls == Walls::kMixed)
   {
      for (int b = 0; b < m->GetNBE(); b++)
      {
         if (m->GetBdrAttribute(b) != 5) { continue; }
         Array<int> v;
         m->GetBdrElementVertices(b, v);
         real_t xc = 0.0;
         for (int i = 0; i < v.Size(); i++) { xc += m->GetVertex(v[i])[0]; }
         xc /= v.Size();
         if (std::abs(xc + kL) < 1e-6 * kL) { m->SetBdrAttribute(b, 6); }
         else if (std::abs(xc - kL) < 1e-6 * kL) { m->SetBdrAttribute(b, 7); }
      }
      m->SetAttributes();
   }
   return m;
}

const IFaultClosure &Closure()
{
   static DieterichRuinaFriction fr{{TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                                     TPV102Params::Dc}};
   static RateStateClosure cl(fr, RateStateClosure::StateLaw::Aging);
   return cl;
}

/// The TPV102 fault of test_poro_p0.cpp (sigma_yy = -120 MPa, its broad
/// gentle nucleation), without a pore pressure.
NitscheFaultSetup FaultSetup()
{
   NitscheFaultSetup fs;
   fs.closure = &Closure();
   fs.prestress = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = TPV102Params::tau_ini;
      s[4] = -TPV102Params::sigma_n;
   };
   fs.nucleation_spatial = [](const real_t x[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      const real_t dx = x[0] - TPV102Params::hypo_along_strike;
      const real_t dz = -x[2] - TPV102Params::hypo_down_dip;
      const real_t r2 = (dx * dx + dz * dz) / (20.0e3 * 20.0e3);
      s[1] = s[3] = 12.0e6 * std::exp(-r2);
   };
   fs.nucleation_temporal = [](real_t t) { return mfem::seas::NucleationTemporal(t); };
   fs.parameters = [](const real_t x[3], FaultPointState &st) {
      st.a = mfem::seas::ComputeA(x[0], -x[2]);
      st.Dc = TPV102Params::Dc;
   };
   fs.sigma_min = 0.0;
   return fs;
}

/// d = 0, V = V_ini along strike, psi in equilibrium (as test_poro_p0.cpp).
void SeedState(NitscheWaveOperator &op)
{
   const int no = op.NumOwnedFaultDOFs();
   Vector d(3 * no), V(3 * no), psi(no), S;
   d = 0.0;
   op.BackgroundTraction(0.0, S);
   for (int j = 0; j < no; j++)
   {
      const real_t *f = &op.FaultFrames()(9 * j);
      real_t X[2], sig;
      mfem::seas::fault_frame::ToLocal(&S(3 * j), 0.0, f, f + 3, f + 6, X[0], X[1], sig);
      FaultPointState st = op.FaultParams(j);
      Closure().SeedEquilibrium(std::hypot(X[0], X[1]), TPV102Params::V_ini, sig, st);
      psi(j) = st.psi;
      for (int c = 0; c < 3; c++) { V(3 * j + c) = TPV102Params::V_ini * f[6 + c]; }
   }
   op.SetFaultState(d, V, psi);
}

/// Everything one operator needs, kept alive together.
struct PoroBox
{
   std::unique_ptr<ParMesh> pm;
   std::unique_ptr<PWConstCoefficient> lam, mu, rho;
   MaterialField mf;
   std::unique_ptr<PoroMaterialSample> sample;
   NitscheFaultSetup fs;
   BoundaryConfig bc;
   PoroBoundaryConfig pbc;
   std::unique_ptr<PoroNitscheOperator> op;
};

typedef std::function<void(BoundaryConfig &, PoroBoundaryConfig &, real_t &)> Tweak;

/// The operator on `serial` (from BoxMesh) over `comm`; two_media = F-C.
std::unique_ptr<PoroBox> MakeBox(Mesh &serial, MPI_Comm comm, int p, bool two_media,
                                 Walls walls, real_t kappa_w, const DirichletFunc &g,
                                 const Tweak &tweak = nullptr)
{
   auto B = std::make_unique<PoroBox>();
   B->pm = std::make_unique<ParMesh>(comm, serial);
   const Medium &mb = two_media ? kMedB : kMedA;
   Vector lv(2), mv(2), rv(2);
   lv(0) = kMedA.lam; lv(1) = mb.lam;
   mv(0) = kMedA.mu;  mv(1) = mb.mu;
   rv(0) = BulkDensity(kMedA.f); rv(1) = BulkDensity(mb.f);
   B->lam = std::make_unique<PWConstCoefficient>(lv);
   B->mu = std::make_unique<PWConstCoefficient>(mv);
   B->rho = std::make_unique<PWConstCoefficient>(rv);
   B->mf.mode = MaterialField::Mode::Coefficient;
   B->mf.lambda_coef = B->lam.get();
   B->mf.mu_coef = B->mu.get();
   B->mf.rho_coef = B->rho.get();
   B->sample = std::make_unique<PoroMaterialSample>(
                  *B->pm, B->mf, [two_media](int, const Vector &x)
   { return (two_media && x(0) > 0.0) ? kMedB.f : kMedA.f; });
   B->fs = FaultSetup();
   B->bc.fault_attr = 3;
   switch (walls)
   {
      case Walls::kDrained:
         B->bc.natural_attrs = {1};
         B->bc.dirichlet_attrs = {5};
         break;
      case Walls::kMixed:
         B->bc.natural_attrs = {1};
         B->bc.dirichlet_attrs = {5, 6};
         B->bc.absorbing_attrs = {7};
         B->pbc.no_flow_attrs = {6, 7};
         break;
      case Walls::kClamped:
         B->bc.dirichlet_attrs = {1, 5};
         break;
   }
   B->bc.default_dirichlet_func = g;
   if (tweak) { tweak(B->bc, B->pbc, kappa_w); }
   B->op = std::make_unique<PoroNitscheOperator>(*B->pm, p, *B->sample, B->bc, kappa_w,
                                                 B->pbc, mfem::seas::DomainConfig{}, &B->fs);
   return B;
}

DirichletFunc ZeroG()
{
   return [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
}

/// A time-dependent far field (as test_poro_p0.cpp), so Force exercises the
/// Dirichlet loads.
DirichletFunc FarField()
{
   return [](const Vector &x, real_t t, Vector &u) {
      u.SetSize(3);
      u(0) = 1e-6 * t * x(1) / kL;
      u(1) = -3e-7 * t * x(2) / kL;
      u(2) = 2e-7 * t;
   };
}

typedef std::function<void(const Vector &cen, const Vector &x, Vector &v)> NodeFn;

/// f(centroid of e, x) at every node of every local element, as a vector of
/// the space's vector DOFs.
Vector ProjectNodes(ParFiniteElementSpace &fes, const NodeFn &f)
{
   ParMesh &m = *fes.GetParMesh();
   Vector g(fes.GetVSize()), cen(3), x(3), val(3);
   Array<int> vd;
   IsoparametricTransformation T;
   for (int e = 0; e < m.GetNE(); e++)
   {
      m.GetElementTransformation(e, &T);
      T.Transform(Geometries.GetCenter(Geometry::TETRAHEDRON), cen);
      const FiniteElement &fe = *fes.GetFE(e);
      const IntegrationRule &nodes = fe.GetNodes();
      const int nd = fe.GetDof();
      fes.GetElementVDofs(e, vd);
      for (int k = 0; k < nd; k++)
      {
         T.Transform(nodes.IntPoint(k), x);
         f(cen, x, val);
         for (int c = 0; c < 3; c++) { g(vd[c * nd + k]) = val(c); }
      }
   }
   return g;
}

Vector Join(const Vector &a, const Vector &b)
{
   Vector x(a.Size() + b.Size());
   for (int i = 0; i < a.Size(); i++) { x(i) = a(i); }
   for (int i = 0; i < b.Size(); i++) { x(a.Size() + i) = b(i); }
   return x;
}

ParFiniteElementSpace &Space(PoroBox &B) { return B.op->Elastic().Elasticity().GetFESpace(); }

/// max_i |v_i| over rows [lo, hi) and all ranks.
real_t MaxAbs(const Vector &v, int lo, int hi, MPI_Comm c)
{
   real_t m = 0.0;
   for (int i = lo; i < hi; i++) { m = std::max(m, std::abs(v(i))); }
   return RMax(m, c);
}

/// max_ij |A_ij x_j| for a one-rank HypreParMatrix (no off-diagonal part).
real_t MaxTerm(const HypreParMatrix &A, const Vector &x)
{
   SparseMatrix D;
   A.GetDiag(D);
   const int *I = D.GetI(), *J = D.GetJ();
   const real_t *a = D.GetData();
   real_t m = 0.0;
   for (int i = 0; i < D.Height(); i++)
   {
      for (int k = I[i]; k < I[i + 1]; k++) { m = std::max(m, std::abs(a[k] * x(J[k]))); }
   }
   return m;
}

/// The Rayleigh quotient after `iters` power iterations on K: a lower bound
/// of |K|_2 (K is symmetric).
real_t LambdaMax(const PoroNitscheOperator &op, MPI_Comm comm, int iters)
{
   Vector x(op.Height()), y;
   x.Randomize(977 + 13 * g_rank);
   x /= std::sqrt(InnerProduct(comm, x, x));
   real_t lam = 0.0;
   for (int it = 0; it < iters; it++)
   {
      op.ApplyStiffness(x, y);
      lam = InnerProduct(comm, x, y);
      x = y;
      x /= std::sqrt(InnerProduct(comm, x, x));
   }
   return lam;
}

/// For each local vector DOF of a, the vector DOF of b at the same node and
/// component (elements matched by centroid, nodes by position); both spaces
/// on partitions of the same serial mesh.
std::vector<int> MatchBoxDofs(ParFiniteElementSpace &a, ParFiniteElementSpace &b, real_t &gap)
{
   ParMesh &ma = *a.GetParMesh(), &mb = *b.GetParMesh();
   std::vector<Vector> cb(mb.GetNE());
   IsoparametricTransformation Ta, Tb;
   for (int k = 0; k < mb.GetNE(); k++)
   {
      mb.GetElementTransformation(k, &Tb);
      cb[k].SetSize(3);
      Tb.Transform(Geometries.GetCenter(Geometry::TETRAHEDRON), cb[k]);
   }
   std::vector<int> map(a.GetVSize(), -1);
   gap = 0.0;
   Array<int> va, vb;
   Vector ca(3), xa(3), xb(3), d(3);
   for (int e = 0; e < ma.GetNE(); e++)
   {
      ma.GetElementTransformation(e, &Ta);
      Ta.Transform(Geometries.GetCenter(Geometry::TETRAHEDRON), ca);
      int E = -1;
      for (int k = 0; k < mb.GetNE() && E < 0; k++)
      {
         d = cb[k];
         d -= ca;
         if (d.Normlinf() < 1e-9 * kL) { E = k; }
      }
      if (E < 0) { gap = 1e30; return map; }
      mb.GetElementTransformation(E, &Tb);
      const IntegrationRule &nodes = a.GetFE(e)->GetNodes();
      const int nd = a.GetFE(e)->GetDof();
      a.GetElementVDofs(e, va);
      b.GetElementVDofs(E, vb);
      for (int k = 0; k < nd; k++)
      {
         Ta.Transform(nodes.IntPoint(k), xa);
         int K = -1;
         real_t best = 1e30;
         for (int l = 0; l < nd; l++)
         {
            Tb.Transform(nodes.IntPoint(l), xb);
            xb -= xa;
            if (xb.Normlinf() < best) { best = xb.Normlinf(); K = l; }
         }
         gap = std::max(gap, best);
         for (int c = 0; c < 3; c++) { map[va[c * nd + k]] = vb[c * nd + K]; }
      }
   }
   return map;
}

std::string Box(bool two_media) { return two_media ? "F-C" : "F-B"; }

// ------------------------------------------------------------- U2-3

void TestU23()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-3: K symmetric and the same on every rank count (mixed walls: "
                   "Dirichlet, Dirichlet + no-flow, absorbing + no-flow, free)\n";
   }
   const int np = Mpi::WorldSize();
   for (int two = 0; two <= 1; two++)
   {
      for (int p = 1; p <= 3; p++)
      {
         const std::string P = "U2-3 " + Box(two) + " p = " + std::to_string(p) + ": ";
         std::unique_ptr<Mesh> m(BoxMesh(4, 4, 2, Walls::kMixed));
         auto box = MakeBox(*m, MPI_COMM_WORLD, p, two, Walls::kMixed, 1.0, ZeroG());
         const PoroNitscheOperator &op = *box->op;
         const int N = op.Height();
         Vector x(N), y(N), Kx, Ky;
         x.Randomize(101 + 7 * g_rank);
         y.Randomize(203 + 7 * g_rank);
         op.ApplyStiffness(x, Kx);
         op.ApplyStiffness(y, Ky);
         const real_t yKx = InnerProduct(MPI_COMM_WORLD, y, Kx);
         const real_t xKy = InnerProduct(MPI_COMM_WORLD, x, Ky);
         const real_t nx = std::sqrt(InnerProduct(MPI_COMM_WORLD, x, x));
         const real_t ny = std::sqrt(InnerProduct(MPI_COMM_WORLD, y, y));
         const real_t lmax = LambdaMax(op, MPI_COMM_WORLD, 40);
         const real_t rs = std::abs(yKx - xKy) / (nx * lmax * ny);
         TEST_ASSERT(rs <= kR1, P + "|y^T K x - x^T K y| / (|x| |K| |y|) = " + Fmt(rs)
                     + " for random x, y (|K| = " + Fmt(lmax) + ", 40 power iterations)");

         // The same K x on one rank: the operator on MPI_COMM_SELF, x a
         // field of the node positions with a jump between elements.
         auto self = MakeBox(*m, MPI_COMM_SELF, p, two, Walls::kMixed, 1.0, ZeroG());
         NodeFn fu = [](const Vector &cen, const Vector &x, Vector &v) {
            for (int c = 0; c < 3; c++)
            {
               v(c) = std::sin(1.3 * x(0) / kL + 0.7 * x(1) / kL - 0.4 * x(2) / kL + c)
                      + 0.2 * std::cos(5.1 * cen(0) / kL - 3.3 * cen(2) / kL + 2 * c);
            }
         };
         NodeFn fw = [](const Vector &cen, const Vector &x, Vector &v) {
            for (int c = 0; c < 3; c++)
            {
               v(c) = std::cos(0.9 * x(0) / kL - 1.1 * x(1) / kL + 0.6 * x(2) / kL - c)
                      + 0.3 * std::sin(4.7 * cen(1) / kL + 2.9 * cen(0) / kL + c);
            }
         };
         ParFiniteElementSpace &fp = Space(*box), &fs = Space(*self);
         const Vector xp = Join(ProjectNodes(fp, fu), ProjectNodes(fp, fw));
         const Vector xs = Join(ProjectNodes(fs, fu), ProjectNodes(fs, fw));
         Vector kp, ks;
         op.ApplyStiffness(xp, kp);
         self->op->ApplyStiffness(xs, ks);
         const int n = op.NumU(), ns = self->op->NumU();
         Vector xsu(const_cast<Vector &>(xs), 0, ns), xsw(const_cast<Vector &>(xs), ns, ns);
         PoroNitscheOperator &so = *self->op;
         const real_t big = std::max({MaxTerm(so.Elastic().Elasticity().ParStiffness(), xsu),
                                      MaxTerm(so.Zuw(), xsw), MaxTerm(so.Zwu(), xsu),
                                      MaxTerm(so.Zww(), xsw)});
         real_t gap;
         const std::vector<int> map = MatchBoxDofs(fp, fs, gap);
         real_t diff = 0.0;
         for (int i = 0; i < n; i++)
         {
            diff = std::max(diff, std::abs(kp(i) - ks(map[i])));
            diff = std::max(diff, std::abs(kp(n + i) - ks(ns + map[i])));
         }
         diff = RMax(diff, MPI_COMM_WORLD);
         gap = RMax(gap, MPI_COMM_WORLD);
         const real_t rr = diff / big;
         TEST_ASSERT(gap < 1e-9 * kL && rr <= kRV,
                     P + "K x on " + std::to_string(np) + " rank(s) equals K x on one rank "
                     "(R = " + Fmt(rr) + " of the largest term, assembled vector)");
      }
   }
}

// ------------------------------------------------------------- U2-4

void TestU24()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-4: locked fluid, the u rows against the elastic operator built "
                   "from (lambda_u, mu, rho) (F-B, mixed walls, slip, Dirichlet, damping)\n";
   }
   const real_t rho = BulkDensity(kMedA.f);
   const BiotConstants c = DeriveBiotConstants(kMedA.lam, kMedA.mu, rho, kMedA.f);
   for (int p = 1; p <= 3; p++)
   {
      for (int kv = 0; kv <= 1; kv++)
      {
         const std::string P = "U2-4 p = " + std::to_string(p) + (kv ? ", Kelvin-Voigt" : "")
                               + ": ";
         std::unique_ptr<Mesh> m(BoxMesh(4, 4, 2, Walls::kMixed));
         auto box = MakeBox(*m, MPI_COMM_WORLD, p, false, Walls::kMixed, 1.0, FarField());
         PoroNitscheOperator &op = *box->op;
         NitscheFaultSetup fs = FaultSetup();
         NitscheWaveOperator ref(*box->pm, p, MaterialField::MakeConstant(c.lambda_u, kMedA.mu, rho),
                                 box->bc, {}, &fs);
         if (kv)
         {
            ref.SetKelvinVoigt(0.1, NitscheWaveOperator::KVScale::Node, -1.0);
            op.Elastic().SetKelvinVoigt(0.1, NitscheWaveOperator::KVScale::Node, -1.0);
         }
         SeedState(ref);
         SeedState(op.Elastic());
         ref.SetTime(0.3);
         op.SetTime(0.3);
         const int n = ref.Height(), no = ref.NumOwnedFaultDOFs();
         MFEM_VERIFY(n == op.NumU(), "U2-4: the two operators differ in size");
         Vector xu(n), vu(n), vw(n), w0(n), d(3 * no);
         xu.Randomize(11 + g_rank);
         vu.Randomize(23 + g_rank);
         vw.Randomize(29 + g_rank);
         d.Randomize(37 + g_rank);
         xu *= 1e-3;
         vu *= 1e-2;
         vw *= 1e-2;
         w0 = 0.0;
         const Vector x = Join(xu, w0), v = Join(vu, vw);
         Vector Fref, F;
         ref.Force(xu, &vu, &d, Fref);
         op.Force(x, &v, &d, F);
         Vector Fu(F, 0, n);
         const long b1 = LSum(BitDiff(Fu, Fref), MPI_COMM_WORLD);
         // Kelvin-Voigt needs the velocity (NitscheWaveOperator refuses null).
         ref.Force(xu, kv ? &vu : nullptr, nullptr, Fref);
         op.Force(x, kv ? &v : nullptr, nullptr, F);
         Vector Fu2(F, 0, n);
         const long b2 = LSum(BitDiff(Fu2, Fref), MPI_COMM_WORLD);
         // HasDamping is per rank: some ranks hold no absorbing face.
         const long damped = LSum(ref.HasDamping() ? 1 : 0, MPI_COMM_WORLD);
         const long damped_op = LSum(op.Elastic().HasDamping() ? 1 : 0, MPI_COMM_WORLD);
         TEST_ASSERT(damped > 0 && damped == damped_op,
                     P + "both operators have absorbing faces, on the same "
                     + std::to_string(damped) + " rank(s)");
         TEST_ASSERT(b1 == 0, P + "with slip d, velocity and w = 0, the u rows of Force equal "
                     "the elastic Force bit for bit (" + std::to_string(b1) + " differ)");
         TEST_ASSERT(b2 == 0, P + "without slip" + (kv ? "" : " and velocity")
                     + ", likewise (" + std::to_string(b2) + " differ)");
      }
   }
}

// ------------------------------------------------------------- U2-5

/// Dense copy of a one-rank HypreParMatrix into K at (r0, c0).
void AddDense(const HypreParMatrix &A, DenseMatrix &K, int r0, int c0)
{
   SparseMatrix D;
   A.GetDiag(D);
   const int *I = D.GetI(), *J = D.GetJ();
   const real_t *a = D.GetData();
   for (int i = 0; i < D.Height(); i++)
   {
      for (int k = I[i]; k < I[i + 1]; k++) { K(r0 + i, c0 + J[k]) += a[k]; }
   }
}

/// The fluid penalty kappa_w of (6.5) that U2-5 fixes: half of it must leave K
/// positive semidefinite (plan section 6.3, a margin of 2).  At 1.0, the
/// derivation's starting value, half was not enough (first run: lambda_min /
/// lambda_max = -3.2e-3 to -4.2e-5); K turns semidefinite at 0.953-0.968 for
/// p = 1, 0.603-0.610 for p = 2 and 0.544-0.546 for p = 3 (bisection on this
/// box: F-B and F-C, drained and mixed walls; p = 3 on F-B mixed only).
/// kappa_w = 2 chosen by Chunhui Zhao on 2026-10-21 in the review of step P2.2.
const real_t kKappaW = 2.0;

void TestU25()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-5: coercivity at half the fluid penalty (kappa_w = " << 0.5 * kKappaW
                << " of " << kKappaW << "): every eigenvalue of K by LAPACK on the 2 x 2 x 1 "
                   "box (drained walls; mixed walls on F-B)\n";
   }
   if (Mpi::WorldSize() != 1)
   {
      if (g_rank == 0) { std::cout << "  U2-5 runs on 1 rank only (dense eigenvalues)\n"; }
      return;
   }
#ifndef MFEM_USE_LAPACK
   TEST_ASSERT(false, "U2-5 needs MFEM built with LAPACK (DenseMatrix::Eigenvalues)");
   return;
#else
   for (int c = 0; c < 3; c++)
   {
      // F-B and F-C with drained walls, F-B with mixed walls (the no-flow
      // faces carry the one-sided penalty).
      const bool two = (c == 1);
      const Walls walls = (c == 2) ? Walls::kMixed : Walls::kDrained;
      for (int p = 1; p <= 3; p++)
      {
         const std::string P = "U2-5 " + Box(two) + (c == 2 ? " mixed" : " drained") + " p = "
                               + std::to_string(p) + ": ";
         std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, walls));
         auto box = MakeBox(*m, MPI_COMM_WORLD, p, two, walls, 0.5 * kKappaW, ZeroG());
         PoroNitscheOperator &op = *box->op;
         const int n = op.NumU();
         const double t0 = MPI_Wtime();
         DenseMatrix K(2 * n);
         K = 0.0;
         AddDense(op.Elastic().Elasticity().ParStiffness(), K, 0, 0);
         AddDense(op.Zuw(), K, 0, n);
         AddDense(op.Zwu(), K, n, 0);
         AddDense(op.Zww(), K, n, n);
         Vector ev;
         K.Eigenvalues(ev);
         const real_t lmax = ev.Max(), lmin = ev.Min();
         long zeros = 0;
         real_t first_pos = lmax;
         for (int i = 0; i < ev.Size(); i++)
         {
            if (std::abs(ev(i)) <= 1e-10 * lmax) { zeros++; }
            else if (ev(i) > 0.0) { first_pos = std::min(first_pos, ev(i)); }
         }
         // The zero modes: (0, w) with w in BDM_p (continuous normal
         // component, w.n = 0 on the no-flow faces) and div w = 0 in each
         // element; div maps that space onto the discontinuous P_{p-1}
         // (drained faces leave the total flux free), so their number is
         // its dimension minus dim P_{p-1}.
         long Fnf = 0;
         for (int b = 0; b < box->pm->GetNBE(); b++)
         {
            if (box->pbc.no_flow_attrs.count(box->pm->GetBdrAttribute(b))) { Fnf++; }
         }
         const long T = box->pm->GetNE(), F = box->pm->GetNumFaces() - Fnf;
         const long N0 = F * (p + 1) * (p + 2) / 2 + T * (p + 1) * (p + 2) * (p - 1) / 2
                         - T * p * (p + 1) * (p + 2) / 6;
         if (g_rank == 0)
         {
            std::cout << "  " << P << 2 * n << " unknowns, " << T << " tets, " << F
                      << " faces without no-flow: lambda_min = " << Fmt(lmin) << ", lambda_max = " << Fmt(lmax)
                      << ", smallest positive above the zeros " << Fmt(first_pos) << " ("
                      << Fmt(first_pos / lmax) << " lambda_max); " << Fmt(MPI_Wtime() - t0)
                      << " s\n";
         }
         TEST_ASSERT(lmin >= -1e-10 * lmax, P + "no eigenvalue below -1e-10 lambda_max "
                     "(lambda_min / lambda_max = " + Fmt(lmin / lmax) + ")");
         TEST_ASSERT(zeros == N0, P + std::to_string(zeros) + " eigenvalues within 1e-10 "
                     "lambda_max of zero, the dimension of the divergence-free BDM_p fields ("
                     + std::to_string(N0) + ")");
      }
   }
#endif
}

// ------------------------------------------------------------- U2-6, U2-7

/// Residual of the force F against the parts it sums (loads L and K x), in
/// rows [lo, hi): max |F| / max(max |L|, max |K x|).
real_t ResidualR(const Vector &F, const Vector &L, const Vector &Kx, int lo, int hi)
{
   const real_t s = std::max(MaxAbs(L, lo, hi, MPI_COMM_WORLD), MaxAbs(Kx, lo, hi, MPI_COMM_WORLD));
   return MaxAbs(F, lo, hi, MPI_COMM_WORLD) / s;
}

void TestU26U27a()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-6, U2-7(a): patch tests PT-1 and PT-2 on F-C (exact u on every "
                   "outer face, drained), residual of the interpolated exact solution\n";
   }
   {
      // The revision 2.2 jump [[(alpha u + w).n]] = (alpha_A - alpha_B) u(0).
      const real_t rho = BulkDensity(kMedA.f);
      const real_t aA = DeriveBiotConstants(kMedA.lam, kMedA.mu, rho, kMedA.f).alpha;
      const real_t aB = DeriveBiotConstants(kMedB.lam, kMedB.mu, BulkDensity(kMedB.f),
                                            kMedB.f).alpha;
      const real_t j = (aA - aB) * poro_refs::patch_PT1::u_at_interface;
      const real_t r = std::abs(j - poro_refs::patch_PT1::old_form_jump)
                       / poro_refs::patch_PT1::old_form_jump;
      TEST_ASSERT(r <= kR1, "U2-6 the revision 2.2 form would leave [[(alpha u + w).n]] = "
                  + Fmt(j) + " m on x = 0 (poro_test_refs.hpp " +
                  Fmt(poro_refs::patch_PT1::old_form_jump) + ", R = " + Fmt(r) + ")");
   }
   for (int which = 0; which <= 1; which++)
   {
      const real_t duA = which ? poro_refs::patch_PT2::A_du_dx : poro_refs::patch_PT1::A_du_dx;
      const real_t duB = which ? poro_refs::patch_PT2::B_du_dx : poro_refs::patch_PT1::B_du_dx;
      const real_t dwA = which ? poro_refs::patch_PT2::A_dw_dx : poro_refs::patch_PT1::A_dw_dx;
      const real_t dwB = which ? poro_refs::patch_PT2::B_dw_dx : poro_refs::patch_PT1::B_dw_dx;
      const real_t u0 = poro_refs::patch_PT1::u_at_interface;
      const real_t V = which ? poro_refs::patch_PT2::V : 0.0;
      DirichletFunc g = [=](const Vector &x, real_t, Vector &u) {
         u.SetSize(3);
         u = 0.0;
         u(0) = u0 + (x(0) < 0.0 ? duA : duB) * x(0);
      };
      NodeFn fu = [=](const Vector &cen, const Vector &x, Vector &v) {
         v = 0.0;
         v(0) = u0 + (cen(0) < 0.0 ? duA : duB) * x(0);
      };
      NodeFn fw = [=](const Vector &cen, const Vector &x, Vector &v) {
         v = 0.0;
         v(0) = (cen(0) < 0.0 ? dwA : dwB) * x(0);
      };
      for (int p = 1; p <= 3; p++)
      {
         const std::string P = std::string(which ? "U2-7(a) PT-2" : "U2-6 PT-1") + " p = "
                               + std::to_string(p) + ": ";
         std::unique_ptr<Mesh> m(BoxMesh(4, 4, 2, Walls::kClamped));
         auto box = MakeBox(*m, MPI_COMM_WORLD, p, true, Walls::kClamped, 1.0, g);
         PoroNitscheOperator &op = *box->op;
         ConstantCoefficient Vc(V);
         if (which) { op.SetInjectedVolume(&Vc); }
         op.SetTime(0.0);
         ParFiniteElementSpace &fes = Space(*box);
         const Vector x = Join(ProjectNodes(fes, fu), ProjectNodes(fes, fw));
         Vector zero(x.Size()), F, L, Kx;
         zero = 0.0;
         op.Force(x, nullptr, nullptr, F);
         op.Force(zero, nullptr, nullptr, L);
         op.ApplyStiffness(x, Kx);
         const int n = op.NumU();
         const real_t ru = ResidualR(F, L, Kx, 0, n), rw = ResidualR(F, L, Kx, n, 2 * n);
         TEST_ASSERT(ru <= kRV, P + "u rows: the force of the exact solution is zero (R = "
                     + Fmt(ru) + " of max(|loads|, |K x|))");
         TEST_ASSERT(rw <= kRV, P + "w rows: likewise (R = " + Fmt(rw) + ")");
      }
   }
}

/// int over the local elements of f(x) . phi_k (component c) and, on the
/// boundary faces with an attribute in `faces`, of h(x) phi_k n_c, into the
/// rows [off, off + n) of E (the rule has order 2p + 4).
void AssembleExpected(ParFiniteElementSpace &fes, int off,
                      const std::function<void(const Vector &, Vector &)> &f,
                      const std::set<int> &faces,
                      const std::function<void(const Vector &, const Vector &, Vector &)> &h,
                      Vector &E)
{
   ParMesh &m = *fes.GetParMesh();
   Array<int> vd;
   Vector shape, x(3), val(3), nor(3);
   for (int e = 0; e < m.GetNE(); e++)
   {
      const FiniteElement &fe = *fes.GetFE(e);
      const int nd = fe.GetDof();
      shape.SetSize(nd);
      ElementTransformation *T = m.GetElementTransformation(e);
      const IntegrationRule &ir = IntRules.Get(Geometry::TETRAHEDRON, 2 * fe.GetOrder() + 4);
      fes.GetElementVDofs(e, vd);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         T->SetIntPoint(&ip);
         T->Transform(ip, x);
         fe.CalcShape(ip, shape);
         f(x, val);
         const real_t w = ip.weight * T->Weight();
         for (int c = 0; c < 3; c++)
            for (int k = 0; k < nd; k++) { E(off + vd[c * nd + k]) += w * shape(k) * val(c); }
      }
   }
   for (int b = 0; b < m.GetNBE(); b++)
   {
      if (!faces.count(m.GetBdrAttribute(b))) { continue; }
      FaceElementTransformations *FT = m.GetBdrFaceTransformations(b);
      if (!FT) { continue; }
      const FiniteElement &fe = *fes.GetFE(FT->Elem1No);
      const int nd = fe.GetDof();
      shape.SetSize(nd);
      fes.GetElementVDofs(FT->Elem1No, vd);
      const IntegrationRule &ir = IntRules.Get(FT->GetGeometryType(), 2 * fe.GetOrder() + 4);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         FT->SetAllIntPoints(&ip);
         const IntegrationPoint &eip = FT->GetElement1IntPoint();
         FT->Elem1->Transform(eip, x);
         fe.CalcShape(eip, shape);
         CalcOrtho(FT->Jacobian(), nor);
         const real_t nl = nor.Norml2();
         Vector nhat(nor);
         nhat /= nl;
         h(x, nhat, val);
         const real_t w = ip.weight * nl;
         for (int c = 0; c < 3; c++)
            for (int k = 0; k < nd; k++) { E(off + vd[c * nd + k]) += w * shape(k) * val(c); }
      }
   }
}

void TestU27bc()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-7(b, c): the source loads on F-B with mixed walls against their "
                   "integrated-by-parts forms\n";
   }
   const real_t rho = BulkDensity(kMedA.f);
   const BiotConstants c = DeriveBiotConstants(kMedA.lam, kMedA.mu, rho, kMedA.f);
   const real_t aM = c.alpha * c.M, M = c.M;
   // A quadratic starting pressure and a linear injected volume.
   const real_t P0 = 1e7, V0 = 1e-4;
   auto pinit = [=](const Vector &x) {
      const real_t a = x(0) / kL, b = x(1) / kL, z = x(2) / kL;
      return P0 * (1.0 + 0.4 * a - 0.3 * b + 0.2 * z + 0.5 * a * b - 0.25 * z * z + 0.3 * a * a);
   };
   auto gradp = [=](const Vector &x, Vector &g) {
      const real_t a = x(0) / kL, b = x(1) / kL, z = x(2) / kL;
      g.SetSize(3);
      g(0) = P0 * (0.4 + 0.5 * b + 0.6 * a) / kL;
      g(1) = P0 * (-0.3 + 0.5 * a) / kL;
      g(2) = P0 * (0.2 - 0.5 * z) / kL;
   };
   auto vol = [=](const Vector &x) {
      return V0 * (1.0 + 0.3 * x(0) / kL - 0.2 * x(1) / kL + 0.1 * x(2) / kL);
   };
   const real_t gradV[3] = {0.3 * V0 / kL, -0.2 * V0 / kL, 0.1 * V0 / kL};
   for (int p = 1; p <= 3; p++)
   {
      const std::string P = "U2-7 p = " + std::to_string(p) + ": ";
      std::unique_ptr<Mesh> m(BoxMesh(4, 4, 2, Walls::kMixed));
      auto box = MakeBox(*m, MPI_COMM_WORLD, p, false, Walls::kMixed, 1.0, ZeroG());
      PoroNitscheOperator &op = *box->op;
      ParFiniteElementSpace &fes = Space(*box);
      const int n = op.NumU();
      Vector zero(2 * n), F, E(2 * n);
      zero = 0.0;

      // (b) F_init: -int grad p_init . q + int_drained p_init (q.n); the
      // no-flow faces 6, 7 carry none.
      FunctionCoefficient pc(pinit);
      op.SetInitialPressure(pc);
      op.SetTime(0.0);
      op.Force(zero, nullptr, nullptr, F);
      E = 0.0;
      AssembleExpected(fes, n, [&](const Vector &x, Vector &v) { gradp(x, v); v.Neg(); },
                       {1, 5}, [&](const Vector &x, const Vector &nh, Vector &v)
      { v = nh; v *= pinit(x); }, E);
      real_t rb = 0.0, sb = 0.0;
      for (int i = n; i < 2 * n; i++) { rb = std::max(rb, std::abs(F(i) - E(i))); }
      sb = std::max(MaxAbs(F, n, 2 * n, MPI_COMM_WORLD), MaxAbs(E, n, 2 * n, MPI_COMM_WORLD));
      rb = RMax(rb, MPI_COMM_WORLD) / sb;
      TEST_ASSERT(rb <= kRV, P + "F_init of a quadratic p_init equals -int grad p_init . q + "
                  "int p_init q.n on the drained faces (R = " + Fmt(rb) + ")");

      // (c) F_V of a linear V: solid rows -int alpha M grad V . v + int alpha M V
      // v.n on the faces that are not Dirichlet (1, 7); fluid rows
      // -int M grad V . q + int M V q.n on the drained faces (1, 5).
      auto box2 = MakeBox(*m, MPI_COMM_WORLD, p, false, Walls::kMixed, 1.0, ZeroG());
      FunctionCoefficient Vc(vol);
      box2->op->SetInjectedVolume(&Vc);
      box2->op->SetTime(0.0);
      box2->op->Force(zero, nullptr, nullptr, F);
      ParFiniteElementSpace &fes2 = Space(*box2);
      E = 0.0;
      AssembleExpected(fes2, 0, [&](const Vector &, Vector &v)
      { for (int k = 0; k < 3; k++) { v(k) = -aM * gradV[k]; } },
      {1, 7}, [&](const Vector &x, const Vector &nh, Vector &v) { v = nh; v *= aM * vol(x); }, E);
      AssembleExpected(fes2, n, [&](const Vector &, Vector &v)
      { for (int k = 0; k < 3; k++) { v(k) = -M * gradV[k]; } },
      {1, 5}, [&](const Vector &x, const Vector &nh, Vector &v) { v = nh; v *= M * vol(x); }, E);
      for (int half = 0; half <= 1; half++)
      {
         const int lo = half * n, hi = lo + n;
         real_t r = 0.0;
         for (int i = lo; i < hi; i++) { r = std::max(r, std::abs(F(i) - E(i))); }
         const real_t s = std::max(MaxAbs(F, lo, hi, MPI_COMM_WORLD), MaxAbs(E, lo, hi, MPI_COMM_WORLD));
         r = RMax(r, MPI_COMM_WORLD) / s;
         TEST_ASSERT(r <= kRV, P + (half ? "fluid" : "solid") + " rows of F_V for a linear V "
                     "equal their integrated-by-parts form (R = " + Fmt(r) + ")");
      }
   }
}

// ------------------------------------------------------------- U2-8

/// A scalar polynomial in (x, y, z) / kL.
struct SPoly
{
   std::map<std::array<int, 3>, real_t> t;
   real_t Eval(const Vector &x) const
   {
      real_t s = 0.0;
      for (const auto &kv : t)
      {
         s += kv.second * Poly::Pow(x(0) / kL, kv.first[0]) * Poly::Pow(x(1) / kL, kv.first[1])
              * Poly::Pow(x(2) / kL, kv.first[2]);
      }
      return s;
   }
   /// d/dx_d in physical units.
   SPoly D(int d) const
   {
      SPoly r;
      for (const auto &kv : t)
      {
         if (kv.first[d] == 0) { continue; }
         std::array<int, 3> e = kv.first;
         e[d]--;
         r.t[e] += kv.second * kv.first[d] / kL;
      }
      return r;
   }
   void Add(const SPoly &o, real_t s)
   {
      for (const auto &kv : o.t) { t[kv.first] += s * kv.second; }
   }
};

/// Three random polynomials of total degree deg (coefficients in [-1, 1]).
std::array<SPoly, 3> RandomField(int deg, uint32_t seed)
{
   std::array<SPoly, 3> f;
   uint32_t st = seed;
   for (int c = 0; c < 3; c++)
      for (int i = 0; i <= deg; i++)
         for (int j = 0; i + j <= deg; j++)
            for (int k = 0; i + j + k <= deg; k++)
            {
               st = 1664525u * st + 1013904223u;
               f[c].t[{i, j, k}] = 2.0 * ((st >> 8) / 16777216.0) - 1.0;
            }
   return f;
}

SPoly Div(const std::array<SPoly, 3> &f)
{
   SPoly d = f[0].D(0);
   d.Add(f[1].D(1), 1.0);
   d.Add(f[2].D(2), 1.0);
   return d;
}

void TestU28()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-8: consistency on F-B (drained walls): for polynomial (u, w) of "
                   "degree p, the force against the strong form integrated against each "
                   "test function\n";
   }
   const real_t rho = BulkDensity(kMedA.f);
   const BiotConstants c = DeriveBiotConstants(kMedA.lam, kMedA.mu, rho, kMedA.f);
   const real_t lu = c.lambda_u, mu = c.mu, aM = c.alpha * c.M, M = c.M, al = c.alpha;
   for (int p = 1; p <= 3; p++)
   {
      const std::string P = "U2-8 p = " + std::to_string(p) + ": ";
      const std::array<SPoly, 3> u = RandomField(p, 4321u + p), w = RandomField(p, 8765u + p);
      // r_u = -div sigma = -(mu lap u + (lambda_u + mu) grad div u + alpha M grad div w),
      // r_w = -grad(M zeta) = -M (alpha grad div u + grad div w), and on the
      // faces sigma n and M zeta.
      const SPoly du = Div(u), dw = Div(w);
      std::array<SPoly, 3> ru, rw;
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++) { ru[i].Add(u[i].D(j).D(j), -mu); }
         ru[i].Add(du.D(i), -(lu + mu));
         ru[i].Add(dw.D(i), -aM);
         rw[i].Add(du.D(i), -M * al);
         rw[i].Add(dw.D(i), -M);
      }
      std::array<std::array<SPoly, 3>, 3> grad;
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++) { grad[i][j] = u[i].D(j); }
      auto sigma_n = [&](const Vector &x, const Vector &nh, Vector &v) {
         const real_t dv = du.Eval(x), dww = dw.Eval(x);
         for (int i = 0; i < 3; i++)
         {
            v(i) = (lu * dv + aM * dww) * nh(i);
            for (int j = 0; j < 3; j++) { v(i) += mu * (grad[i][j].Eval(x) + grad[j][i].Eval(x)) * nh(j); }
         }
      };
      auto Mzeta_n = [&](const Vector &x, const Vector &nh, Vector &v) {
         v = nh;
         v *= M * (al * du.Eval(x) + dw.Eval(x));
      };
      DirichletFunc g = [&](const Vector &x, real_t, Vector &v) {
         v.SetSize(3);
         for (int i = 0; i < 3; i++) { v(i) = u[i].Eval(x); }
      };
      std::unique_ptr<Mesh> m(BoxMesh(4, 4, 2, Walls::kDrained));
      auto box = MakeBox(*m, MPI_COMM_WORLD, p, false, Walls::kDrained, 1.0, g);
      PoroNitscheOperator &op = *box->op;
      op.SetTime(0.0);
      ParFiniteElementSpace &fes = Space(*box);
      const int n = op.NumU();
      NodeFn fu = [&](const Vector &, const Vector &x, Vector &v) {
         for (int i = 0; i < 3; i++) { v(i) = u[i].Eval(x); }
      };
      NodeFn fw = [&](const Vector &, const Vector &x, Vector &v) {
         for (int i = 0; i < 3; i++) { v(i) = w[i].Eval(x); }
      };
      const Vector x = Join(ProjectNodes(fes, fu), ProjectNodes(fes, fw));
      Vector zero(2 * n), F, L, Kx, E(2 * n);
      zero = 0.0;
      op.Force(x, nullptr, nullptr, F);
      op.Force(zero, nullptr, nullptr, L);
      op.ApplyStiffness(x, Kx);
      // -F = int r . (v, q) + int_free sigma n . v + int_all M zeta q.n
      E = 0.0;
      AssembleExpected(fes, 0, [&](const Vector &y, Vector &v)
      { for (int i = 0; i < 3; i++) { v(i) = -ru[i].Eval(y); } },
      {1}, [&](const Vector &y, const Vector &nh, Vector &v) { sigma_n(y, nh, v); v.Neg(); }, E);
      AssembleExpected(fes, n, [&](const Vector &y, Vector &v)
      { for (int i = 0; i < 3; i++) { v(i) = -rw[i].Eval(y); } },
      {1, 5}, [&](const Vector &y, const Vector &nh, Vector &v) { Mzeta_n(y, nh, v); v.Neg(); }, E);
      for (int half = 0; half <= 1; half++)
      {
         const int lo = half * n, hi = lo + n;
         real_t r = 0.0;
         for (int i = lo; i < hi; i++) { r = std::max(r, std::abs(F(i) - E(i))); }
         const real_t s = std::max({MaxAbs(L, lo, hi, MPI_COMM_WORLD),
                                    MaxAbs(Kx, lo, hi, MPI_COMM_WORLD),
                                    MaxAbs(E, lo, hi, MPI_COMM_WORLD)});
         r = RMax(r, MPI_COMM_WORLD) / s;
         TEST_ASSERT(r <= kRV, P + (half ? "w" : "u") + " rows of the force equal the strong "
                     "form with its boundary terms (R = " + Fmt(r) + " of max(|loads|, |K x|, "
                     "|strong form|))");
      }
   }
}

// ------------------------------------------------------------- X2-6 (aborts)

int RunAbortP22(const std::string &mode)
{
   if (g_rank == 0) { std::cout << "X2-6 expecting an abort: " << mode << "\n"; }
   std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, Walls::kDrained));
   Tweak tweak;
   if (mode == "dirichlet-interior")
   {
      // attribute 8 on the interior faces of the plane x = 0, made Dirichlet
      std::vector<std::array<int, 3>> tris;
      for (int f = 0; f < m->GetNumFaces(); f++)
      {
         if (m->GetFaceInformation(f).IsBoundary()) { continue; }
         Array<int> fv;
         m->GetFaceVertices(f, fv);
         bool on = true;
         for (int i = 0; i < fv.Size(); i++)
         {
            if (std::abs(m->GetVertex(fv[i])[0]) > 1e-6 * kL) { on = false; }
         }
         if (on) { tris.push_back({fv[0], fv[1], fv[2]}); }
      }
      for (const auto &t : tris) { m->AddBdrTriangle(t[0], t[1], t[2], 8); }
      m->FinalizeTopology();
      m->Finalize();
      m->SetAttributes();
      tweak = [](BoundaryConfig &bc, PoroBoundaryConfig &, real_t &) { bc.dirichlet_attrs = {5, 8}; };
   }
   else if (mode == "no-flow-fault")
   {
      tweak = [](BoundaryConfig &, PoroBoundaryConfig &pbc, real_t &) { pbc.no_flow_attrs = {3}; };
   }
   else if (mode == "kappa-zero")
   {
      tweak = [](BoundaryConfig &, PoroBoundaryConfig &, real_t &k) { k = 0.0; };
   }
   else
   {
      if (g_rank == 0) { std::cout << "FAILED: unknown abort mode " << mode << "\n"; }
      return 1;
   }
   auto box = MakeBox(*m, MPI_COMM_WORLD, 1, false, Walls::kDrained, 1.0, ZeroG(), tweak);
   if (g_rank == 0) { std::cout << "FAILED: X2-6 " << mode << " did not abort\n"; }
   return 1;
}

// ------------------------------------------------------------- U2-9, X2-7

using mfem::seas::poro::InitialPressure;
using mfem::seas::poro::PoroInitialPressure;
using mfem::seas::poro::PoroSource;
using mfem::seas::poro::PoroSourceConfig;
using mfem::seas::poro::SourceIntegration;

const real_t kRS = 1e-15;     // R of U2-9 (plan section 6.4)
const real_t kRise = 0.125;   // s, a power of two, so the trapezoid grids are exact

/// The plan's time functions, written here independently of poro_source.hpp.
real_t RefRate(const std::string &kind, real_t t)
{
   if (t < 0.0) { return 0.0; }
   if (kind == "step" || t >= kRise) { return 1.0; }
   const real_t s = t / kRise;
   if (kind == "ramp") { return s; }
   return 10.0 * s * s * s - 15.0 * s * s * s * s + 6.0 * s * s * s * s * s;
}

/// int_0^t RefRate by Gauss-Legendre on [0, min(t, T)] and [T, t], exact for
/// these polynomials.
real_t RefIntegral(const std::string &kind, real_t t)
{
   const IntegrationRule &ir = IntRules.Get(Geometry::SEGMENT, 11);
   auto gauss = [&](real_t a, real_t b) {
      real_t q = 0.0;
      for (int i = 0; i < ir.GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir.IntPoint(i);
         q += ip.weight * RefRate(kind, a + (b - a) * ip.x);
      }
      return (b - a) * q;
   };
   if (t <= 0.0) { return 0.0; }
   if (kind == "step") { return gauss(0.0, t); }
   return t <= kRise ? gauss(0.0, t) : gauss(0.0, kRise) + gauss(kRise, t);
}

PoroSourceConfig SourceCfg(const std::string &shape, const std::string &kind, real_t rate,
                           const real_t x0[3])
{
   PoroSourceConfig c;
   c.on = true;
   c.shape = shape;
   for (int d = 0; d < 3; d++) { c.origin_m[d] = x0[d]; }
   c.width_per_m2 = shape == "gaussian" ? 1.0 / (0.4 * kL * 0.4 * kL) : 0.0;
   c.rate = rate;
   c.time_function = kind;
   c.rise_time_s = kind == "step" ? 0.0 : kRise;
   return c;
}

void TestU29()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-9: PoroSource: V(t) for each time function, the trapezoid "
                   "rule's order, and a well's point load\n";
   }
   const real_t origin[3] = {0.0, 0.0, -0.5 * kL};
   const char *kinds[3] = {"step", "ramp", "smooth_step"};

   // V(t) exact: the closed form against Gauss-Legendre of the plan's f.
   const real_t times[6] = {0.0, 0.013, 0.05, kRise, 0.1371, 0.25};
   for (const char *kind : kinds)
   {
      PoroSource src(SourceCfg("gaussian", kind, 1.0, origin));
      real_t rv = 0.0, rf = 0.0;
      bool zero = true;
      for (real_t t : times)
      {
         src.Advance(t);
         const real_t ref = RefIntegral(kind, t);
         if (t == 0.0) { zero = src.VolumeFactor() == 0.0; }
         else { rv = std::max(rv, std::abs(src.VolumeFactor() - ref) / std::abs(ref)); }
         rf = std::max(rf, std::abs(src.TimeFunction().Rate(t) - RefRate(kind, t)));
      }
      TEST_ASSERT(zero && rv <= kRS && rf <= kRS,
                  std::string("U2-9 ") + kind + ": V(t) = rate F(t) equals the integral of "
                  "the plan's f at 6 times, 0 at t = 0 (R = " + Fmt(rv) + "; f within "
                  + Fmt(rf) + ")");
   }

   // The trapezoid rule: second order where F is not piecewise linear (the
   // smooth step during its rise), exact on the grid otherwise.
   {
      real_t err[4];
      const real_t tstar = 0.5 * kRise;
      for (int k = 0; k < 4; k++)
      {
         const int nsteps = 4 << k;   // dt = T/8 ... T/64
         const real_t dt = tstar / nsteps;
         PoroSource src(SourceCfg("gaussian", "smooth_step", 1.0, origin),
                        SourceIntegration::kTrapezoid);
         for (int i = 1; i <= nsteps; i++) { src.Advance(i * dt); }
         err[k] = std::abs(src.VolumeFactor() - RefIntegral("smooth_step", tstar));
      }
      bool ok = true;
      std::string rs;
      for (int k = 0; k < 3; k++)
      {
         const real_t r = err[k] / err[k + 1];
         ok = ok && r >= 3.5 && r <= 4.5;
         rs += (k ? ", " : "") + Fmt(r);
      }
      TEST_ASSERT(ok, "U2-9 smooth_step: the trapezoid V at T/2 converges at second order "
                  "(dt = T/8 to T/64: errors " + Fmt(err[0]) + " to " + Fmt(err[3])
                  + ", ratios " + rs + " in [3.5, 4.5])");
      for (const char *kind : {"step", "ramp"})
      {
         PoroSource src(SourceCfg("gaussian", kind, 1.0, origin),
                        SourceIntegration::kTrapezoid);
         const real_t dt = kRise / 8.0;
         for (int i = 1; i <= 16; i++) { src.Advance(i * dt); }
         const real_t ref = RefIntegral(kind, 2.0 * kRise);
         const real_t r = std::abs(src.VolumeFactor() - ref) / ref;
         TEST_ASSERT(r <= kRS, std::string("U2-9 ") + kind + ": the trapezoid V is exact on "
                     "a grid through T (f is linear between grid times; R = " + Fmt(r) + ")");
      }
   }

   // A well: the load is M V_w zeta(v, q)(x_w) on its element's rows only.
   const real_t wells[2][3] = {{-0.613 * kL, 0.287 * kL, -0.391 * kL},
                               {0.442 * kL, -0.533 * kL, -0.718 * kL}};
   for (int fc = 0; fc <= 1; fc++)
   {
      for (int p = 1; p <= 3; p++)
      {
         std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, Walls::kDrained));
         auto box = MakeBox(*m, MPI_COMM_WORLD, p, fc == 1, Walls::kDrained, kKappaW, ZeroG());
         PoroNitscheOperator &op = *box->op;
         ParFiniteElementSpace &fes = Space(*box);
         const int n = op.NumU();
         for (int wi = 0; wi <= fc; wi++)
         {
            const real_t *xw = wells[wi];
            const std::string P = "U2-9 " + Box(fc == 1) + " p = " + std::to_string(p)
                                  + ", well at x " + (xw[0] < 0 ? "< 0" : "> 0") + ": ";
            const Medium &md = xw[0] > 0.0 && fc == 1 ? kMedB : kMedA;
            const BiotConstants c = DeriveBiotConstants(md.lam, md.mu, BulkDensity(md.f), md.f);
            PoroSource well(SourceCfg("well", "smooth_step", 3.0e-3, xw),
                            SourceIntegration::kExact, &fes, box->sample.get());
            const Vector &L = well.WellLoad();
            const int e = well.WellElement();
            const long holders = LSum(e >= 0 ? 1 : 0, MPI_COMM_WORLD);
            // nonzeros outside the element's rows
            std::vector<char> mine(2 * n, 0);
            if (e >= 0)
            {
               Array<int> vd;
               fes.GetElementVDofs(e, vd);
               for (int i = 0; i < vd.Size(); i++) { mine[vd[i]] = mine[n + vd[i]] = 1; }
            }
            long outside = 0;
            for (int i = 0; i < 2 * n; i++) { if (!mine[i] && L(i) != 0.0) { outside++; } }
            outside = LSum(outside, MPI_COMM_WORLD);
            TEST_ASSERT(holders == 1 && outside == 0 && L.Size() == 2 * n,
                        P + "one element on one rank holds it, and every other entry of "
                        "the load is zero (no face terms; " + std::to_string(outside)
                        + " nonzero outside)");

            // b^T L for polynomial (u, w) of degree p = M (alpha div u + div w)(x_w).
            const std::array<SPoly, 3> u = RandomField(p, 2468u + p), w = RandomField(p, 1357u + p);
            NodeFn fu = [&](const Vector &, const Vector &x, Vector &v) {
               for (int i = 0; i < 3; i++) { v(i) = u[i].Eval(x); }
            };
            NodeFn fw = [&](const Vector &, const Vector &x, Vector &v) {
               for (int i = 0; i < 3; i++) { v(i) = w[i].Eval(x); }
            };
            // Measured against the sum of the |terms| x_i L_i, chosen by Chunhui
            // Zhao on 2026-10-21 in the review of step P2.3: on the first run
            // p = 3 was 4.0e-15 to 4.8e-15 of the largest term, from the round-off
            // of MFEM's p = 3 basis gradients at x_w (their sum, zero in exact
            // arithmetic, is 21 eps of the largest one), and 3.6e-16 to 4.6e-16
            // of the sum.
            const Vector x = Join(ProjectNodes(fes, fu), ProjectNodes(fes, fw));
            real_t sums[2] = {0.0, 0.0};   // b^T L, sum |b_i L_i|
            for (int i = 0; i < 2 * n; i++)
            {
               sums[0] += x(i) * L(i);
               sums[1] += std::abs(x(i) * L(i));
            }
            MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            Vector X(const_cast<real_t *>(xw), 3);
            const real_t exact = c.M * (c.alpha * Div(u).Eval(X) + Div(w).Eval(X));
            const real_t r = std::abs(sums[0] - exact) / std::max(sums[1], std::abs(exact));
            TEST_ASSERT(r <= kRS, P + "for polynomial (u, w) of degree p the load times "
                        "their DOFs is M (alpha div u + div w)(x_w) (R = " + Fmt(r)
                        + " of the sum of |terms|)");

            // Through the operator: Force adds rate F(t) times that load.
            const real_t t1 = 0.0731;
            well.Advance(t1);
            op.SetTime(t1);
            Vector F1, F0;
            op.SetSource(&well);
            op.Force(x, nullptr, nullptr, F1);
            op.SetSource(nullptr);
            op.Force(x, nullptr, nullptr, F0);
            well.AddWellLoad(F0);
            long diff = 0;
            for (int i = 0; i < 2 * n; i++) { if (F1(i) != F0(i)) { diff++; } }
            diff = LSum(diff, MPI_COMM_WORLD);
            TEST_ASSERT(diff == 0, P + "Force with the well equals Force without it plus "
                        "rate F(t) times the load, bit for bit");
         }
      }
   }
}

/// X2-7: the gaussian source and p_init reach the operator.
void TestX27()
{
   if (g_rank == 0)
   {
      std::cout << "--- X2-7 (extra): a gaussian PoroSource through SetSource, and "
                   "InitialPressure\n";
   }
   const real_t x0[3] = {-0.21 * kL, 0.13 * kL, -0.47 * kL};
   for (int p = 1; p <= 2; p++)
   {
      const std::string P = "X2-7 F-B p = " + std::to_string(p) + ": ";
      std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, Walls::kMixed));
      auto box = MakeBox(*m, MPI_COMM_WORLD, p, false, Walls::kMixed, kKappaW, ZeroG());
      PoroNitscheOperator &op = *box->op;
      const int n = op.NumU();
      const PoroSourceConfig cfg = SourceCfg("gaussian", "ramp", 2.0e-6, x0);
      PoroSource src(cfg);
      const real_t t1 = 0.0917;
      src.Advance(t1);
      op.SetTime(t1);
      FunctionCoefficient Vref([&](const Vector &x) {
         real_t r2 = 0.0;
         for (int d = 0; d < 3; d++) { r2 += (x(d) - x0[d]) * (x(d) - x0[d]); }
         return cfg.rate * RefIntegral("ramp", t1) * std::exp(-cfg.width_per_m2 * r2);
      });
      Vector zero(2 * n), F1, F0;
      zero = 0.0;
      op.SetSource(&src);
      op.Force(zero, nullptr, nullptr, F1);
      op.SetSource(nullptr);
      op.SetInjectedVolume(&Vref);
      op.Force(zero, nullptr, nullptr, F0);
      op.SetInjectedVolume(nullptr);
      for (int half = 0; half <= 1; half++)
      {
         const int lo = half * n, hi = lo + n;
         real_t r = 0.0;
         for (int i = lo; i < hi; i++) { r = std::max(r, std::abs(F1(i) - F0(i))); }
         r = RMax(r, MPI_COMM_WORLD) / MaxAbs(F0, lo, hi, MPI_COMM_WORLD);
         TEST_ASSERT(r <= kRV, P + (half ? "fluid" : "solid") + " rows: F_V from the "
                     "source equals F_V of rate g(x) F(t) given directly (R = " + Fmt(r) + ")");
      }
   }
   // p_init at points.
   PoroInitialPressure ip;
   ip.on = true;
   ip.origin_m[0] = 0.0; ip.origin_m[1] = -1000.0; ip.origin_m[2] = -7500.0;
   ip.magnitude_pa = 5.0e7;
   ip.width_per_m2 = 1.0e-7;
   std::unique_ptr<Coefficient> pc = InitialPressure(ip);
   std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, Walls::kDrained));
   IsoparametricTransformation T;
   real_t r = 0.0;
   const IntegrationRule &ir = IntRules.Get(Geometry::TETRAHEDRON, 4);
   for (int e = 0; e < m->GetNE(); e += 7)
   {
      m->GetElementTransformation(e, &T);
      for (int k = 0; k < ir.GetNPoints(); k++)
      {
         Vector x;
         T.SetIntPoint(&ir.IntPoint(k));
         T.Transform(ir.IntPoint(k), x);
         real_t r2 = 0.0;
         for (int d = 0; d < 3; d++) { r2 += (x(d) - ip.origin_m[d]) * (x(d) - ip.origin_m[d]); }
         const real_t ref = ip.magnitude_pa * std::exp(-ip.width_per_m2 * r2);
         if (ref > 0.0) { r = std::max(r, std::abs(pc->Eval(T, ir.IntPoint(k)) - ref) / ref); }
      }
   }
   TEST_ASSERT(r <= kRS, "X2-7 InitialPressure is magnitude exp(-width |x - x0|^2) at the "
               "quadrature points (R = " + Fmt(r) + ")");
}

/// U2-9 and X2-7 aborts.
int RunAbortP23(const std::string &mode)
{
   if (g_rank == 0) { std::cout << "expecting an abort: " << mode << "\n"; }
   std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, Walls::kDrained));
   real_t xw[3] = {2.0 * kL, 0.0, -0.5 * kL};   // well-outside
   if (mode == "well-face")
   {
      // the centroid of the first interior face that is not on the fault plane
      for (int f = 0; f < m->GetNumFaces(); f++)
      {
         if (m->GetFaceInformation(f).IsBoundary()) { continue; }
         Array<int> fv;
         m->GetFaceVertices(f, fv);
         real_t c[3] = {0.0, 0.0, 0.0};
         for (int i = 0; i < 3; i++)
            for (int d = 0; d < 3; d++) { c[d] += m->GetVertex(fv[i])[d] / 3.0; }
         if (std::abs(c[1]) < 1e-6 * kL) { continue; }
         for (int d = 0; d < 3; d++) { xw[d] = c[d]; }
         break;
      }
   }
   else if (mode == "well-vertex")
   {
      // the first vertex strictly inside the box
      for (int v = 0; v < m->GetNV(); v++)
      {
         const real_t *c = m->GetVertex(v);
         if (std::abs(c[0]) < 0.99 * kL && std::abs(c[1]) < 0.99 * kL && c[2] < -0.01 * kL
             && c[2] > -0.99 * kL && std::abs(c[1]) > 1e-6 * kL)
         {
            for (int d = 0; d < 3; d++) { xw[d] = c[d]; }
            break;
         }
      }
   }
   else if (mode != "well-outside" && mode != "source-time")
   {
      if (g_rank == 0) { std::cout << "FAILED: unknown abort mode " << mode << "\n"; }
      return 1;
   }
   auto box = MakeBox(*m, MPI_COMM_WORLD, 1, false, Walls::kDrained, kKappaW, ZeroG());
   if (mode == "source-time")
   {
      PoroSource src(SourceCfg("gaussian", "step", 1.0e-6, xw));
      box->op->SetSource(&src);
      box->op->SetTime(0.25);
      Vector x(2 * box->op->NumU()), F;
      x = 0.0;
      box->op->Force(x, nullptr, nullptr, F);
   }
   else
   {
      if (g_rank == 0)
      {
         std::cout << "  well at (" << xw[0] << ", " << xw[1] << ", " << xw[2] << ")\n";
      }
      PoroSource well(SourceCfg("well", "step", 1.0e-3, xw), SourceIntegration::kExact,
                      &Space(*box), box->sample.get());
   }
   if (g_rank == 0) { std::cout << "FAILED: " << mode << " did not abort\n"; }
   return 1;
}

// ------------------------------------------------------------- step P2.4

using mfem::seas::poro::ElementConstants;
using mfem::seas::poro::GammaUW;
using mfem::seas::poro::PoroMass;
using mfem::seas::poro::PoroMaxDt;
using mfem::seas::poro::VelocityCoupling;

/// z = dt / tau of U2-10 (plan section 6.4).
const real_t kZ[4] = {1e-2, 1.0, 2600.0, 2e5};
/// CFL number for dt in U2-10 and U2-11.
const real_t kCfl = 0.5;

/// gamma_w(z) = 1 / (1 - e^-z) - 1 / z, derivation (9.4): an input of U2-10
/// and U2-11 (FittedGamma arrives in step P2.5).
real_t GammaOfZ(real_t z) { return 1.0 / (-std::expm1(-z)) - 1.0 / z; }

/// The sample's constants with the drag b set so that z = dt / tau = dt b rho
/// / (rho m - rho_f^2) (eq. (8.4)) on every element.
std::vector<BiotConstants> ConstantsAtZ(const PoroMaterialSample &s, real_t dt, real_t z)
{
   std::vector<BiotConstants> c = ElementConstants(s);
   for (BiotConstants &k : c) { k.b = z * (k.rho * k.m - k.rho_f * k.rho_f) / (dt * k.rho); }
   return c;
}

/// The vector DG mass with one value per element, by MFEM's
/// VectorMassIntegrator and a 2p + 2 rule (the NitscheMass loop uses 2p).
std::unique_ptr<HypreParMatrix> ElementMass(ParFiniteElementSpace &fes,
                                            const std::vector<real_t> &val)
{
   ParMesh &m = *fes.GetParMesh();
   L2_FECollection fec0(0, 3);
   ParFiniteElementSpace fes0(&m, &fec0);
   ParGridFunction gf(&fes0);
   Array<int> d;
   for (int e = 0; e < m.GetNE(); e++)
   {
      fes0.GetElementDofs(e, d);
      gf(d[0]) = val[e];
   }
   GridFunctionCoefficient q(&gf);
   const int p = fes.GetMaxElementOrder();
   VectorMassIntegrator *vmi = new VectorMassIntegrator(q);
   vmi->SetIntRule(&IntRules.Get(Geometry::TETRAHEDRON, 2 * p + 2));
   ParBilinearForm a(&fes);
   a.AddDomainIntegrator(vmi);
   a.Assemble(0);
   a.Finalize(0);
   return std::unique_ptr<HypreParMatrix>(a.ParallelAssemble());
}

/// F = (M + dt D Gamma) x of (8.2), (9.5), assembled independently of
/// PoroMass: the three element masses above and, on the u rows, the elastic
/// operator's dashpot dt C (1/2 x_u + g_uw x_w) (ApplyDamping).
Vector ApplyMassDrag(PoroNitscheOperator &op, const std::vector<BiotConstants> &c,
                     real_t dt, const Vector &gw, VelocityCoupling coupling,
                     const Vector &x)
{
   ParFiniteElementSpace &fes = op.Elastic().Elasticity().GetFESpace();
   ParMesh &m = *fes.GetParMesh();
   const int n = op.NumU(), ne = m.GetNE();
   std::vector<real_t> rho(ne), rhof(ne), mp(ne), guw(ne);
   for (int e = 0; e < ne; e++)
   {
      rho[e] = c[e].rho;
      rhof[e] = c[e].rho_f;
      mp[e] = c[e].m + dt * gw(e) * c[e].b;
      guw[e] = GammaUW(coupling, c[e].rho, c[e].rho_f, gw(e));
   }
   auto Mr = ElementMass(fes, rho), Mf = ElementMass(fes, rhof), Mm = ElementMass(fes, mp);
   Vector xu(x.GetData(), n), xw(x.GetData() + n, n), F(2 * n), t(n), v(n), Cv;
   Vector Fu(F.GetData(), n), Fw(F.GetData() + n, n);
   Mr->Mult(xu, Fu);
   Mf->Mult(xw, t);
   Fu += t;
   Mf->Mult(xu, Fw);
   Mm->Mult(xw, t);
   Fw += t;
   Array<int> vd;
   for (int e = 0; e < ne; e++)
   {
      fes.GetElementVDofs(e, vd);
      for (int i : vd) { v(i) = 0.5 * xu(i) + guw[e] * xw(i); }
   }
   op.Elastic().ApplyDamping(v, Cv);
   Fu.Add(dt, Cv);
   return F;
}

/// A field that changes from element to element, the same on every rank
/// count: smooth in x plus a jump through the centroid.
Vector MassField(ParFiniteElementSpace &fes, real_t phase)
{
   return ProjectNodes(fes, [phase](const Vector &cen, const Vector &x, Vector &v) {
      for (int c = 0; c < 3; c++)
      {
         v(c) = std::sin(2.1 * x(0) / kL + 1.3 * x(1) / kL - 0.7 * x(2) / kL + c + phase)
                + 0.4 * std::cos(3.7 * cen(0) / kL - 2.3 * cen(1) / kL + 1.1 * cen(2) / kL
                                 + 0.5 * c + 2.0 * phase);
      }
   });
}

/// U2-10: PoroMass against (M + dt D Gamma), as its inverse.
void TestU210()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-10: PoroMass::Apply inverts an independently assembled "
                   "M + dt D Gamma (z = 1e-2, 1, 2600, 2e5; F-C)\n";
   }
   for (Walls walls : {Walls::kDrained, Walls::kMixed})
   {
      const bool absorbing = walls == Walls::kMixed;
      std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, walls));
      for (int p = 1; p <= 3; p++)
      {
         auto box = MakeBox(*m, MPI_COMM_WORLD, p, true, walls, kKappaW, ZeroG());
         PoroNitscheOperator &op = *box->op;
         ParFiniteElementSpace &fes = Space(*box);
         const int n = op.NumU(), ne = box->pm->GetNE();
         const real_t dt = PoroMaxDt(*box->sample, *box->pm, p, kCfl);
         const std::map<int, DenseMatrix> &blocks = op.Elastic().DampingBlocks();
         const long ndamped = LSum(static_cast<long>(blocks.size()), MPI_COMM_WORLD);
         const Vector x = Join(MassField(fes, 0.0), MassField(fes, 1.0));
         const std::string P = std::string("U2-10 F-C p = ") + std::to_string(p)
                               + (absorbing ? ", absorbing faces (" + std::to_string(ndamped)
                                  + " damped elements)" : ", no absorbing faces");
         std::map<VelocityCoupling, std::vector<Vector>> sol;
         for (VelocityCoupling cp : {VelocityCoupling::kMomentum, VelocityCoupling::kDiagonal})
         {
            bool ok = absorbing ? ndamped > 0 : ndamped == 0;
            std::string rs;
            for (real_t z : kZ)
            {
               const std::vector<BiotConstants> c = ConstantsAtZ(*box->sample, dt, z);
               Vector gw(ne);
               gw = GammaOfZ(z);
               PoroMass mass(*box->pm, p, c, dt, gw, cp, &blocks);
               const Vector F = ApplyMassDrag(op, c, dt, gw, cp, x);
               Vector a;
               mass.Apply(F, a);
               sol[cp].push_back(a);
               real_t r[2];
               for (int half = 0; half <= 1; half++)
               {
                  real_t num = 0.0;
                  for (int i = half * n; i < (half + 1) * n; i++)
                  {
                     num = std::max(num, std::abs(a(i) - x(i)));
                  }
                  r[half] = RMax(num, MPI_COMM_WORLD) / MaxAbs(x, half * n, (half + 1) * n,
                                                               MPI_COMM_WORLD);
                  ok = ok && r[half] <= kRV;
               }
               rs += (rs.empty() ? "" : "; ") + Fmt(z) + ": " + Fmt(r[0]) + " / " + Fmt(r[1]);
            }
            TEST_ASSERT(ok, P + std::string(cp == VelocityCoupling::kMomentum ? ", Gamma_mom"
                                            : ", Gamma_diag")
                        + ": a = PoroMass(M + dt D Gamma) x equals x, R of the u / w rows "
                        "at z = " + rs);
         }
         if (!absorbing)
         {
            long diff = 0;
            for (size_t k = 0; k < sol[VelocityCoupling::kMomentum].size(); k++)
            {
               diff += BitDiff(sol[VelocityCoupling::kMomentum][k],
                               sol[VelocityCoupling::kDiagonal][k]);
            }
            diff = LSum(diff, MPI_COMM_WORLD);
            TEST_ASSERT(diff == 0, P + ": Gamma_mom and Gamma_diag give the same a at every "
                        "z, bit for bit (Proposition 9.3; " + std::to_string(diff)
                        + " entries differ)");
         }
      }
   }
}

/// U2-11: with rho_f = 0 the u rows of PoroMass are NitscheMass's.
void TestU211()
{
   if (g_rank == 0)
   {
      std::cout << "--- U2-11: with rho_f = 0 the u rows of PoroMass equal NitscheMass "
                   "(rho, and rho M + dt/2 C where C is present), bit for bit\n";
   }
   const real_t z = 2600.0;
   for (Walls walls : {Walls::kDrained, Walls::kMixed})
   {
      const bool absorbing = walls == Walls::kMixed;
      std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, walls));
      for (int p = 1; p <= 3; p++)
      {
         auto box = MakeBox(*m, MPI_COMM_WORLD, p, true, walls, kKappaW, ZeroG());
         PoroNitscheOperator &op = *box->op;
         ParFiniteElementSpace &fes = Space(*box);
         const int n = op.NumU(), ne = box->pm->GetNE();
         const real_t dt = PoroMaxDt(*box->sample, *box->pm, p, kCfl);
         std::vector<BiotConstants> c = ConstantsAtZ(*box->sample, dt, z);
         for (BiotConstants &k : c) { k.rho_f = 0.0; }
         Vector gw(ne);
         gw = GammaOfZ(z);
         const Vector F = Join(MassField(fes, 2.0), MassField(fes, 3.0));
         const Vector Fu(const_cast<real_t *>(F.GetData()), n);
         Vector ref;
         const mfem::seas::NitscheMass *nm = &op.Elastic().Mass();
         if (absorbing)
         {
            op.Elastic().SetTimeStep(dt);
            nm = op.Elastic().DampedMass();
         }
         nm->Apply(Fu, ref);
         const long ndamped = LSum(static_cast<long>(op.Elastic().DampingBlocks().size()),
                                   MPI_COMM_WORLD);
         for (VelocityCoupling cp : {VelocityCoupling::kMomentum, VelocityCoupling::kDiagonal})
         {
            PoroMass mass(*box->pm, p, c, dt, gw, cp, &op.Elastic().DampingBlocks());
            Vector a;
            mass.Apply(F, a);
            Vector au(a.GetData(), n);
            const long diff = LSum(BitDiff(au, ref), MPI_COMM_WORLD);
            TEST_ASSERT(diff == 0 && (absorbing ? ndamped > 0 : ndamped == 0),
                        std::string("U2-11 F-C p = ") + std::to_string(p)
                        + (absorbing ? ", absorbing faces (" + std::to_string(ndamped)
                           + " damped elements), against NitscheMass(rho, C, dt/2)"
                           : ", no absorbing faces, against NitscheMass(rho)")
                        + (cp == VelocityCoupling::kMomentum ? ", Gamma_mom" : ", Gamma_diag")
                        + ": u rows identical (" + std::to_string(diff) + " differ)");
         }
      }
   }
}

/// U2-10 abort: gamma_w outside [1/2, 1].
int RunAbortP24(const std::string &mode)
{
   if (g_rank == 0) { std::cout << "expecting an abort: " << mode << "\n"; }
   if (mode != "mass-gamma")
   {
      if (g_rank == 0) { std::cout << "FAILED: unknown abort mode " << mode << "\n"; }
      return 1;
   }
   std::unique_ptr<Mesh> m(BoxMesh(2, 2, 1, Walls::kDrained));
   auto box = MakeBox(*m, MPI_COMM_WORLD, 1, false, Walls::kDrained, kKappaW, ZeroG());
   Vector gw(box->pm->GetNE());
   gw = 0.4;
   PoroMass mass(*box->pm, 1, ElementConstants(*box->sample), 1e-3, gw,
                 VelocityCoupling::kMomentum);
   if (g_rank == 0) { std::cout << "FAILED: " << mode << " did not abort\n"; }
   return 1;
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   std::string only;
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--abort" && i + 1 < argc)
      {
         const std::string mode = argv[i + 1];
         if (mode == "interior-only" || mode == "penalty-one-sided") { return RunAbort(mode); }
         if (mode.rfind("well-", 0) == 0 || mode == "source-time") { return RunAbortP23(mode); }
         if (mode.rfind("mass-", 0) == 0) { return RunAbortP24(mode); }
         return RunAbortP22(mode);
      }
      if (a == "--only" && i + 1 < argc) { only = argv[++i]; }
   }
   if (g_rank == 0)
   {
      std::cout << "=== Poroelastic P2 (plan rev 3), steps P2.1 to P2.4, on "
                << Mpi::WorldSize() << " rank(s)" << (only.empty() ? "" : ", only " + only)
                << " ===\n";
   }
   auto run = [&only](const char *id) { return only.empty() || only == id; };
   if (run("U2-1")) { TestU21(); }
   if (run("U2-2")) { TestU22(); }
   if (run("X2-1")) { TestX21(); }
   if (run("X2-2")) { TestX22(); }
   if (run("X2-4")) { TestX24(); }
   if (run("U2-3")) { TestU23(); }
   if (run("U2-4")) { TestU24(); }
   if (run("U2-5")) { TestU25(); }
   if (run("U2-6") || run("U2-7")) { TestU26U27a(); }
   if (run("U2-7")) { TestU27bc(); }
   if (run("U2-8")) { TestU28(); }
   if (run("U2-9")) { TestU29(); }
   if (run("X2-7")) { TestX27(); }
   if (run("U2-10")) { TestU210(); }
   if (run("U2-11")) { TestU211(); }
   int nf = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return nf == 0 ? 0 : 1;
}
