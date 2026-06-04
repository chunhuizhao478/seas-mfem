// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 2 of the fault-dealiasing plan: unit test for the fault-face RESAMPLE
// operator R (dynamic/fault_resample.hpp).  Pure reference-element linear
// algebra — no mesh, no simulation.  Verifies the plan §Phase-2 acceptance:
//   (i)   R reproduces any degree-N field EXACTLY (incl. the top mode N),
//   (ii)  R removes degree->N content (changes a degree-(N+1) field),
//   (iii) R is idempotent (R²=R) with rank (N+1)(N+2)/2,
//   (iv)  R == I when #QP == #DOF (minimal rule) ⇒ a no-op without over-int,
//   plus W-self-adjointness (R is the W-orthogonal projector, matching
//   SeisSol's resample matrices).
//
// The R-002 subtlety (prior plan review) is exercised directly: on a TRIANGLE
// the minimal degree-2N rule is unisolvent only at N<=2; at N=3 it is
// over-determined (#QP > #DOF) so R != I even at the minimal rule.

#include "mfem.hpp"
#include "../../dynamic/fault_resample.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_passed = 0, num_failed = 0, num_total = 0;

#define CHECK(cond, msg) do { \
   num_total++; \
   if (cond) { num_passed++; std::cout << "    PASS: " << msg << "\n"; } \
   else { num_failed++; std::cout << "    FAIL [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace {

real_t MaxAbs(const DenseMatrix &A)
{
   real_t m = 0.0;
   for (int i = 0; i < A.Height(); i++)
      for (int j = 0; j < A.Width(); j++) { m = std::max(m, std::abs(A(i, j))); }
   return m;
}

// g[q] = x_q^a * y_q^b sampled at the rule's points (a degree-(a+b) monomial).
Vector MonomialField(const IntegrationRule &ir, int a, int b)
{
   Vector g(ir.GetNPoints());
   for (int q = 0; q < ir.GetNPoints(); q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      g(q) = std::pow(ip.x, a) * std::pow(ip.y, b);
   }
   return g;
}

// An arbitrary degree-N field: Σ_{a+b<=N} c_ab x^a y^b with deterministic c_ab.
Vector RandomDegreeNField(const IntegrationRule &ir, int N)
{
   Vector g(ir.GetNPoints()); g = 0.0;
   int k = 0;
   for (int a = 0; a <= N; a++)
      for (int b = 0; a + b <= N; b++)
      {
         Vector m = MonomialField(ir, a, b);
         g.Add(std::sin(0.7 * (k + 1) + 1.3), m);   // nonzero, well-mixed coeffs
         k++;
      }
   return g;
}

void RunCase(int order, int quad_deg)
{
   const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, quad_deg);
   const int nq   = ir.GetNPoints();
   const int ndof = (order + 1) * (order + 2) / 2;
   const bool expect_identity = (nq == ndof);   // R=I iff the rule is unisolvent

   DenseMatrix R;
   BuildFaultResampleMatrix(Geometry::TRIANGLE, order, ir, R);

   std::cout << "  -- order=" << order << "  quad_deg=" << quad_deg
             << "  (#QP=" << nq << ", #DOF=" << ndof << ")"
             << (expect_identity ? "  [#QP==#DOF]" : "  [over-integrated]")
             << "\n";

   CHECK(R.Height() == nq && R.Width() == nq, "R is nq x nq");

   // (iii) idempotent
   DenseMatrix RR(nq, nq); Mult(R, R, RR); RR -= R;
   CHECK(MaxAbs(RR) < 1e-10, "idempotent (R^2 = R)");

   // (iii) rank = trace = (N+1)(N+2)/2
   real_t tr = 0.0; for (int i = 0; i < nq; i++) { tr += R(i, i); }
   CHECK(std::abs(tr - ndof) < 1e-8, "rank = trace = (N+1)(N+2)/2");

   // W-self-adjoint: w_k R(k,j) is symmetric
   real_t wsa = 0.0;
   for (int k = 0; k < nq; k++)
      for (int j = 0; j < nq; j++)
      {
         const real_t a = ir.IntPoint(k).weight * R(k, j);
         const real_t b = ir.IntPoint(j).weight * R(j, k);
         wsa = std::max(wsa, std::abs(a - b));
      }
   CHECK(wsa < 1e-12, "W-self-adjoint (W R = (W R)^T)");

   // (i) reproduces an arbitrary degree-N field EXACTLY (top mode kept)
   Vector g = RandomDegreeNField(ir, order), Rg(nq);
   R.Mult(g, Rg); Rg -= g;
   CHECK(Rg.Normlinf() < 1e-10, "reproduces any degree-N field exactly");

   // (R-005) Conservation: R is the W-orthogonal projector onto a space that
   // contains constants, so (a) a uniform field is reproduced exactly and
   // (b) the W-weighted mean of ANY field is preserved — the secular resample
   // injects no net drift into the slip/state increment (plan §4.2).
   {
      Vector c(nq); c = 0.37;
      Vector Rc(nq); R.Mult(c, Rc); Rc -= c;
      CHECK(Rc.Normlinf() < 1e-10, "reproduces a uniform field exactly (constant kept)");

      Vector v(nq);
      for (int q = 0; q < nq; q++) { v(q) = std::sin(0.9 * q + 0.2) + 0.4 * q; }
      Vector Rv(nq); R.Mult(v, Rv);
      real_t mv = 0.0, mRv = 0.0;
      for (int q = 0; q < nq; q++)
      {
         mv  += ir.IntPoint(q).weight * v(q);
         mRv += ir.IntPoint(q).weight * Rv(q);
      }
      CHECK(std::abs(mv - mRv) < 1e-10 * (std::abs(mv) + 1.0),
            "R preserves the W-weighted mean (conservation; no net drift)");
   }

   // (iv) / (ii)
   DenseMatrix I(nq, nq); I = 0.0; for (int i = 0; i < nq; i++) { I(i, i) = 1.0; }
   DenseMatrix RmI(R); RmI -= I;
   if (expect_identity)
   {
      CHECK(MaxAbs(RmI) < 1e-10, "R == I at the minimal rule (no-op w/o over-int)");
   }
   else
   {
      CHECK(MaxAbs(RmI) > 1e-3, "R != I when over-integrated (#QP > #DOF)");
      // (ii) R annihilates content OUTSIDE the degree-N space.  Split a generic
      // field v = R v + (v - R v): the remainder is the removed > N part, and R
      // must send it to zero (R projects onto the degree-N space; rank < #QP).
      Vector v(nq);
      for (int i = 0; i < nq; i++)
      { v(i) = std::sin(1.1 * i + 0.4) + 0.5 * std::cos(0.3 * i); }
      Vector Rv(nq); R.Mult(v, Rv);
      Vector rem(v); rem -= Rv;
      CHECK(rem.Normlinf() > 1e-6,
            "R removes a nonzero > N component from a generic field");
      Vector Rrem(nq); R.Mult(rem, Rrem);
      CHECK(Rrem.Normlinf() < 1e-10,
            "R annihilates the removed > N content (R(v - R v) = 0)");
   }
}

} // namespace

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
   std::cout << "\n=== Phase 2: fault-face resample operator R ===\n";

   RunCase(1, 2);   // p1 minimal: 3 QP == 3 DOF  -> R = I
   RunCase(1, 6);   // p1 over-int (--fault-overint 2): rank 3
   RunCase(2, 4);   // p2 minimal: 6 QP == 6 DOF  -> R = I
   RunCase(2, 8);   // p2 over-int: rank 6
   RunCase(3, 6);   // p3 minimal deg-2N rule is OVER-determined (R-002): R != I
   RunCase(3, 8);   // p3 over-int: rank 10

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, " << num_failed
             << " failed out of " << num_total << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
