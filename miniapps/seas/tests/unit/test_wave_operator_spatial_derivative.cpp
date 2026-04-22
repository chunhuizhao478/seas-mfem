// ADER I-05 Phase 1: unit test for WaveOperator::ApplySpatialDerivative.
//
// Gates polynomial exactness of the element-local L2-projected spatial
// derivative: for a monomial Q_c(x) = x^n / n the L2-projection of its
// ∂_dir is n * x^{n-1} δ_{dir, coord}, so the recovered derivative must
// match that at every DG DOF (interpolation).  We run n ∈ {1, 2, 3} at
// order = max(2, p) on a uniform hex mesh, checking dQ/dx_j at every DOF
// for each coordinate direction.

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

// Fill Q_c[dof] with f(x_dof) for a single state component, leaving the
// other components zero.  Uses the FE space's nodal positions at the L2
// Lagrange DOFs.
static void FillQFromMonomial(const FiniteElementSpace &fes,
                              int component, int coord_dir, int power,
                              Vector &Q)
{
   const int NS = NUM_STATE;
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NS * ndof_total);
   Q = 0.0;
   Mesh *mesh = fes.GetMesh();
   const int ne = fes.GetNE();
   for (int e = 0; e < ne; e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const IntegrationRule &nodes = fe->GetNodes();
      const int ndof = fe->GetDof();
      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      for (int j = 0; j < ndof; j++)
      {
         const IntegrationPoint &ip = nodes.IntPoint(j);
         Vector xj(3);
         Tr->Transform(ip, xj);
         real_t val = std::pow(xj(coord_dir), power);
         Q[component * ndof_total + edofs[j]] = val;
      }
   }
}

// Evaluate the DG field Q_c at every DOF by reading Q directly at the
// Lagrange nodes (matches how FillQFromMonomial wrote them).
// Returns max |dQ_c_dg[dof] - expected(x_dof)| where expected is provided.
static real_t MaxDOFError(const FiniteElementSpace &fes,
                          int component,
                          const Vector &dQ,
                          int coord_dir, int power)
{
   const int ndof_total = fes.GetNDofs();
   const int ne = fes.GetNE();
   real_t max_err = 0.0;
   for (int e = 0; e < ne; e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const IntegrationRule &nodes = fe->GetNodes();
      const int ndof = fe->GetDof();
      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      for (int j = 0; j < ndof; j++)
      {
         const IntegrationPoint &ip = nodes.IntPoint(j);
         Vector xj(3);
         Tr->Transform(ip, xj);
         // Expected derivative of x_d^n with respect to x_coord_dir
         // equals n * x_coord_dir^{n-1} * δ_{d, coord_dir}.  The Q we
         // filled has component_data[j] = x_coord_dir(j)^power, so its
         // spatial derivative w.r.t. direction `der_dir` is:
         //    power * x_coord_dir^{power-1}    if der_dir == coord_dir
         //    0                                otherwise.
         real_t expected = 0.0;
         // This function is called with (der_dir = coord_dir, power)
         // meaning we computed ∂_{coord_dir} (x_coord_dir^power).
         if (power == 0)
         {
            expected = 0.0;
         }
         else
         {
            expected = static_cast<real_t>(power)
                     * std::pow(xj(coord_dir), power - 1);
         }
         real_t got = dQ[component * ndof_total + edofs[j]];
         real_t err = std::abs(got - expected);
         max_err = std::max(max_err, err);
      }
   }
   return max_err;
}

// Check that ∂_{der_dir} of (x_coord_dir^power) is zero everywhere
// when der_dir != coord_dir.  Returns max |dQ|.
static real_t MaxAbs(const FiniteElementSpace &fes,
                     int component,
                     const Vector &dQ)
{
   const int ndof_total = fes.GetNDofs();
   real_t m = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      m = std::max(m, std::abs(dQ[component * ndof_total + i]));
   }
   return m;
}

int main()
{
   std::cout << "\n=== ADER I-05 Phase 1: ApplySpatialDerivative "
             << "polynomial-exactness test ===\n";

   // DG polynomial order.  Plan acceptance requires exactness for n ∈
   // {1, 2, 3}, which needs a DG space of polynomial degree ≥ 3.
   const int order = 3;
   Mesh mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON,
                                     1.0, 1.0, 1.0);

   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator<Mesh> wave(mesh, order, lambda, mu, rho, bc);

   const FiniteElementSpace &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const real_t kScale = 1.0;  // unit cube — |x| ≤ 1
   // Polynomial-exactness tolerance.  The derivative is computed via L2
   // projection, which is exact (up to FP roundoff) on monomials of
   // degree ≤ order.  Tolerance scales with the number of cross-terms
   // in the sum-over-QPs × number-of-DOFs mat-vec.
   const real_t kEps = 2.220446049250313e-16;
   const real_t tol = 2048 * kEps * kScale;

   // Plan acceptance: n ∈ {1, 2, 3}.  Exactness holds iff the monomial x^n
   // is in the DG polynomial space of degree `order`, i.e., n ≤ order.  At
   // order = 3 all three cases are exact.
   for (int n = 1; n <= 3; n++)
   {
      std::cout << "\n-- n = " << n << " (monomial x^n) --\n";

      // Place the monomial in component SXX (c=0), coordinate direction i.
      for (int i = 0; i < 3; i++)
      {
         Vector Q;
         FillQFromMonomial(fes, /*component=*/SXX,
                           /*coord_dir=*/i, /*power=*/n, Q);

         // ∂_i (x_i^n) = n x_i^{n-1} — reproduced exactly on DG DOFs when
         // n ≤ order (L2 projection of a polynomial of degree n-1 onto a
         // degree-order space is exact).
         Vector dQ;
         wave.ApplySpatialDerivative(i, Q, dQ);
         real_t err_parallel = MaxDOFError(fes, SXX, dQ, i, n);
         std::ostringstream msg;
         msg << "∂_" << i << " x_" << i << "^" << n << " at DG DOFs";
         TEST_LE(err_parallel, tol, msg.str());

         // ∂_j (x_i^n) = 0 for j != i — orthogonality check.
         for (int j = 0; j < 3; j++)
         {
            if (j == i) { continue; }
            Vector dQj;
            wave.ApplySpatialDerivative(j, Q, dQj);
            real_t m = MaxAbs(fes, SXX, dQj);
            std::ostringstream msg2;
            msg2 << "∂_" << j << " x_" << i << "^" << n << " = 0 (j != i)";
            TEST_LE(m, tol, msg2.str());
         }
      }
   }

   // Sanity: input-size error is MFEM_VERIFY'd.  Test with a
   // well-formed zero input yields zero output.
   {
      Vector Q(NUM_STATE * ndof_total);
      Q = 0.0;
      Vector dQ;
      wave.ApplySpatialDerivative(0, Q, dQ);
      real_t m = MaxAbs(fes, SXX, dQ);
      TEST_LE(m, tol, "∂_x(0) = 0");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
