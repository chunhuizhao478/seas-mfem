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
// Usage: mpirun -np 1|2 ./seas_test_poro_p2 [--abort interior-only|penalty-one-sided]

#include "mfem.hpp"
#include "../dg_poro_div_integrator.hpp"
#include "poro_fixtures.hpp"
#include "poro_test_refs.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iostream>
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
      TEST_ASSERT(false, "X2-4 runs on 1 or 2 ranks (F-A has two elements)");
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

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   for (int i = 1; i < argc; i++)
   {
      if (std::string(argv[i]) == "--abort" && i + 1 < argc) { return RunAbort(argv[i + 1]); }
   }
   if (g_rank == 0)
   {
      std::cout << "=== Poroelastic P2 (plan rev 3), step P2.1, on " << Mpi::WorldSize()
                << " rank(s) ===\n";
   }
   TestU21();
   TestU22();
   TestX21();
   TestX22();
   TestX24();
   int nf = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return nf == 0 ? 0 : 1;
}
