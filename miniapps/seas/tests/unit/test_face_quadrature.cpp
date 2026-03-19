// Unit tests for FaceQuadrature class
// Tests L2 projection infrastructure for multi-DOF fault discretization.
//
// Mathematical reference (reference triangle vertices (0,0), (1,0), (0,1)):
//   - Reference area = 1/2
//   - Order 0: nbf=1, φ₀=1, M_ref=1/2, M_ref_inv=2
//   - Order 1: nbf=3, φ_k = {1-x-y, x, y}, M_ref = (1/24)[[2,1,1],[1,2,1],[1,1,2]]
//   - Partition of unity: Σ_k φ_k(q) = 1 for all q at all orders

#include "mfem.hpp"
#include "../../fault/face_quadrature.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

using namespace mfem;
using namespace mfem::seas;

// Test infrastructure
static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(cond, msg) \
   do { \
      num_tests++; \
      if (!(cond)) { \
         std::cerr << "  FAIL: " << (msg) << std::endl; \
         num_failed++; \
      } else { \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(val, expected, tol, msg) \
   do { \
      num_tests++; \
      double _v = (val), _e = (expected), _t = (tol); \
      if (std::fabs(_v - _e) > _t) { \
         std::cerr << "  FAIL: " << (msg) << " got=" << _v \
                   << " expected=" << _e << " diff=" << std::fabs(_v - _e) \
                   << " tol=" << _t << std::endl; \
         num_failed++; \
      } else { \
         num_passed++; \
      } \
   } while (0)

// ============================================================================
// Test 1: Construction and dimensions
// ============================================================================
void TestConstruction()
{
   std::cout << "TestConstruction..." << std::endl;

   // Order 0, vol_order 1: nbf=1
   {
      FaceQuadrature fq(0, 1);
      TEST_ASSERT(fq.NumBasisFunctions() == 1,
                  "order=0: nbf should be 1, got " +
                  std::to_string(fq.NumBasisFunctions()));
      TEST_ASSERT(fq.NumQuadPoints() >= 1,
                  "order=0: nq should be >= 1");
   }

   // Order 1, vol_order 1: nbf=3
   {
      FaceQuadrature fq(1, 1);
      TEST_ASSERT(fq.NumBasisFunctions() == 3,
                  "order=1: nbf should be 3, got " +
                  std::to_string(fq.NumBasisFunctions()));
   }

   // Order 2, vol_order 2: nbf=6
   {
      FaceQuadrature fq(2, 2);
      TEST_ASSERT(fq.NumBasisFunctions() == 6,
                  "order=2: nbf should be 6, got " +
                  std::to_string(fq.NumBasisFunctions()));
   }

   // Order 3, vol_order 3: nbf=10
   {
      FaceQuadrature fq(3, 3);
      TEST_ASSERT(fq.NumBasisFunctions() == 10,
                  "order=3: nbf should be 10, got " +
                  std::to_string(fq.NumBasisFunctions()));
   }

   // BasisAtQuadPoints dimensions
   {
      FaceQuadrature fq(2, 2);
      const DenseMatrix &e_q = fq.BasisAtQuadPoints();
      TEST_ASSERT(e_q.Height() == fq.NumBasisFunctions(),
                  "e_q height should be nbf");
      TEST_ASSERT(e_q.Width() == fq.NumQuadPoints(),
                  "e_q width should be nq");
   }

   // RefMassInverse dimensions
   {
      FaceQuadrature fq(2, 2);
      const DenseMatrix &Minv = fq.RefMassInverse();
      TEST_ASSERT(Minv.Height() == fq.NumBasisFunctions(),
                  "Minv height should be nbf");
      TEST_ASSERT(Minv.Width() == fq.NumBasisFunctions(),
                  "Minv width should be nbf");
   }

   std::cout << "  Done." << std::endl;
}

// ============================================================================
// Test 2: Partition of unity at all orders
// ============================================================================
void TestPartitionOfUnity()
{
   std::cout << "TestPartitionOfUnity..." << std::endl;

   int orders[] = {0, 1, 2, 3};
   int vol_orders[] = {1, 1, 2, 3};

   for (int idx = 0; idx < 4; idx++)
   {
      int p = orders[idx];
      int vp = vol_orders[idx];
      FaceQuadrature fq(p, vp);
      int nbf = fq.NumBasisFunctions();
      int nq = fq.NumQuadPoints();
      const DenseMatrix &e_q = fq.BasisAtQuadPoints();

      for (int q = 0; q < nq; q++)
      {
         real_t sum = 0.0;
         for (int k = 0; k < nbf; k++)
         {
            sum += e_q(k, q);
         }
         std::ostringstream msg;
         msg << "order=" << p << " q=" << q
             << ": sum(phi_k)=" << sum;
         TEST_NEAR(sum, 1.0, 1e-14, msg.str());
      }
   }

   std::cout << "  Done." << std::endl;
}

// ============================================================================
// Test 3: Order 0 is constant with correct M_ref_inv
// ============================================================================
void TestOrder0IsConstant()
{
   std::cout << "TestOrder0IsConstant..." << std::endl;

   FaceQuadrature fq(0, 1);
   int nq = fq.NumQuadPoints();
   const DenseMatrix &e_q = fq.BasisAtQuadPoints();
   const DenseMatrix &Minv = fq.RefMassInverse();

   // Basis should be 1.0 at all quad points
   for (int q = 0; q < nq; q++)
   {
      std::ostringstream msg;
      msg << "order=0: e_q(0," << q << ")";
      TEST_NEAR(e_q(0, q), 1.0, 1e-15, msg.str());
   }

   // M_ref_inv should be 1/Σw = 1/(1/2) = 2.0 for reference triangle
   // (reference triangle area = 1/2, and Σ w_q = 1/2 for MFEM's IntRules)
   real_t sum_w = 0.0;
   const IntegrationRule &ir = fq.GetQuadRule();
   for (int q = 0; q < nq; q++)
   {
      sum_w += ir.IntPoint(q).weight;
   }
   TEST_NEAR(sum_w, 0.5, 1e-14, "Sum of quad weights should be 0.5 (ref tri area)");
   TEST_NEAR(Minv(0, 0), 1.0 / sum_w, 1e-14,
             "M_ref_inv(0,0) should be 1/sum_w");
   TEST_NEAR(Minv(0, 0), 2.0, 1e-14,
             "M_ref_inv(0,0) should be 2.0 for reference triangle");

   // Nodal rule for order 0: centroid (1/3, 1/3)
   const IntegrationRule &nodal = fq.GetNodalRule();
   TEST_ASSERT(nodal.GetNPoints() == 1,
               "order=0 nodal rule should have 1 point");
   TEST_NEAR(nodal.IntPoint(0).x, 1.0/3.0, 1e-14,
             "order=0 nodal x should be 1/3");
   TEST_NEAR(nodal.IntPoint(0).y, 1.0/3.0, 1e-14,
             "order=0 nodal y should be 1/3");

   std::cout << "  Done." << std::endl;
}

// ============================================================================
// Test 4: Mass matrix is SPD and M × M_inv = I
// ============================================================================
void TestMassMatrixSPD()
{
   std::cout << "TestMassMatrixSPD..." << std::endl;

   int orders[] = {1, 2, 3};
   int vol_orders[] = {1, 2, 3};

   for (int idx = 0; idx < 3; idx++)
   {
      int p = orders[idx];
      int vp = vol_orders[idx];
      FaceQuadrature fq(p, vp);
      int nbf = fq.NumBasisFunctions();
      const DenseMatrix &Minv = fq.RefMassInverse();
      const DenseMatrix &e_q = fq.BasisAtQuadPoints();
      const IntegrationRule &ir = fq.GetQuadRule();
      int nq = fq.NumQuadPoints();

      // Reconstruct M_ref
      DenseMatrix M_ref(nbf);
      M_ref = 0.0;
      for (int i = 0; i < nbf; i++)
      {
         for (int j = 0; j < nbf; j++)
         {
            real_t val = 0.0;
            for (int q = 0; q < nq; q++)
            {
               val += ir.IntPoint(q).weight * e_q(i, q) * e_q(j, q);
            }
            M_ref(i, j) = val;
         }
      }

      // Check symmetry
      for (int i = 0; i < nbf; i++)
      {
         for (int j = i+1; j < nbf; j++)
         {
            std::ostringstream msg;
            msg << "p=" << p << " M_ref symmetry (" << i << "," << j << ")";
            TEST_NEAR(M_ref(i, j), M_ref(j, i), 1e-14, msg.str());
         }
      }

      // Check M × M_inv = I
      DenseMatrix I_check(nbf);
      Mult(M_ref, Minv, I_check);
      for (int i = 0; i < nbf; i++)
      {
         for (int j = 0; j < nbf; j++)
         {
            real_t expected = (i == j) ? 1.0 : 0.0;
            std::ostringstream msg;
            msg << "p=" << p << " M*Minv(" << i << "," << j << ")";
            TEST_NEAR(I_check(i, j), expected, 1e-12, msg.str());
         }
      }

      // Check M_inv symmetry
      for (int i = 0; i < nbf; i++)
      {
         for (int j = i+1; j < nbf; j++)
         {
            std::ostringstream msg;
            msg << "p=" << p << " M_inv symmetry (" << i << "," << j << ")";
            TEST_NEAR(Minv(i, j), Minv(j, i), 1e-12, msg.str());
         }
      }

      // Check positive diagonal (necessary for SPD)
      for (int i = 0; i < nbf; i++)
      {
         std::ostringstream msg;
         msg << "p=" << p << " M_ref(" << i << "," << i << ") > 0";
         TEST_ASSERT(M_ref(i, i) > 0.0, msg.str());
      }
   }

   std::cout << "  Done." << std::endl;
}

// ============================================================================
// Test 5: Order 1 mass matrix values (hand-calculated)
// ============================================================================
void TestOrder1MassMatrixValues()
{
   std::cout << "TestOrder1MassMatrixValues..." << std::endl;

   // For linear Lagrange on reference triangle with vertices (0,0),(1,0),(0,1):
   //   φ_0 = 1-x-y,  φ_1 = x,  φ_2 = y
   //   ∫ φ_i φ_j dA:
   //     diagonal: ∫ φ_i² dA = 1/12
   //     off-diag: ∫ φ_i φ_j dA = 1/24
   //   M_ref = (1/24) [[2,1,1],[1,2,1],[1,1,2]]

   FaceQuadrature fq(1, 1);
   int nbf = fq.NumBasisFunctions();
   const DenseMatrix &Minv = fq.RefMassInverse();
   const DenseMatrix &e_q = fq.BasisAtQuadPoints();
   const IntegrationRule &ir = fq.GetQuadRule();
   int nq = fq.NumQuadPoints();

   TEST_ASSERT(nbf == 3, "order=1 should have 3 DOFs");

   // Reconstruct M_ref
   DenseMatrix M_ref(nbf);
   M_ref = 0.0;
   for (int i = 0; i < nbf; i++)
   {
      for (int j = 0; j < nbf; j++)
      {
         real_t val = 0.0;
         for (int q = 0; q < nq; q++)
         {
            val += ir.IntPoint(q).weight * e_q(i, q) * e_q(j, q);
         }
         M_ref(i, j) = val;
      }
   }

   // Note: MFEM's H1_TriangleElement with GaussLobatto may order nodes
   // differently from (0,0),(1,0),(0,1). The mass matrix entries should
   // still satisfy: diagonal = 1/12, off-diagonal = 1/24.
   // This is because the mass matrix of ANY set of 3 linear Lagrange
   // functions on the reference triangle gives this result.
   for (int i = 0; i < 3; i++)
   {
      std::ostringstream msg_d;
      msg_d << "M_ref(" << i << "," << i << ") = 1/12";
      TEST_NEAR(M_ref(i, i), 1.0/12.0, 1e-14, msg_d.str());

      for (int j = i+1; j < 3; j++)
      {
         std::ostringstream msg_o;
         msg_o << "M_ref(" << i << "," << j << ") = 1/24";
         TEST_NEAR(M_ref(i, j), 1.0/24.0, 1e-14, msg_o.str());
      }
   }

   // Verify M_ref_inv hand-calculated values:
   // M_ref = (1/24) [[2,1,1],[1,2,1],[1,1,2]]
   // det([[2,1,1],[1,2,1],[1,1,2]]) = 4
   // adj = [[3,-1,-1],[-1,3,-1],[-1,-1,3]]
   // M_ref_inv = 24 * (1/4) * adj = 6 * [[3,-1,-1],[-1,3,-1],[-1,-1,3]]
   //           = [[18,-6,-6],[-6,18,-6],[-6,-6,18]]
   for (int i = 0; i < 3; i++)
   {
      std::ostringstream msg_d;
      msg_d << "M_inv(" << i << "," << i << ") = 18";
      TEST_NEAR(Minv(i, i), 18.0, 1e-10, msg_d.str());

      for (int j = i+1; j < 3; j++)
      {
         std::ostringstream msg_o;
         msg_o << "M_inv(" << i << "," << j << ") = -6";
         TEST_NEAR(Minv(i, j), -6.0, 1e-10, msg_o.str());
      }
   }

   std::cout << "  Done." << std::endl;
}

// ============================================================================
// Test 6: Roundtrip (interpolate then project recovers nodal values)
// ============================================================================
void TestRoundtrip()
{
   std::cout << "TestRoundtrip..." << std::endl;

   int orders[] = {0, 1, 2, 3};
   int vol_orders[] = {1, 1, 2, 3};

   for (int idx = 0; idx < 4; idx++)
   {
      int p = orders[idx];
      int vp = vol_orders[idx];
      FaceQuadrature fq(p, vp);
      int nbf = fq.NumBasisFunctions();
      int nq = fq.NumQuadPoints();
      int ncomp = 2;  // dip + strike

      // --- Constant roundtrip (should work at all orders) ---
      {
         Vector nodal(ncomp * nbf);
         for (int c = 0; c < ncomp; c++)
         {
            real_t constant_val = (c == 0) ? 3.14 : -2.71;
            for (int k = 0; k < nbf; k++)
            {
               nodal(c * nbf + k) = constant_val;
            }
         }

         Vector quad(ncomp * nq);
         fq.InterpolateToQuadPoints(ncomp, nodal, quad);

         // Verify quad values are constant
         for (int c = 0; c < ncomp; c++)
         {
            real_t expected = (c == 0) ? 3.14 : -2.71;
            for (int q = 0; q < nq; q++)
            {
               std::ostringstream msg;
               msg << "p=" << p << " const interp c=" << c << " q=" << q;
               TEST_NEAR(quad(c * nq + q), expected, 1e-13, msg.str());
            }
         }

         // Project back
         Vector recovered(ncomp * nbf);
         fq.GalerkinProject(ncomp, quad, recovered);

         for (int c = 0; c < ncomp; c++)
         {
            real_t expected = (c == 0) ? 3.14 : -2.71;
            for (int k = 0; k < nbf; k++)
            {
               std::ostringstream msg;
               msg << "p=" << p << " const roundtrip c=" << c << " k=" << k;
               TEST_NEAR(recovered(c * nbf + k), expected, 1e-13, msg.str());
            }
         }
      }

      // --- Linear roundtrip (should work at order >= 1) ---
      if (p >= 1)
      {
         // Set nodal values to a linear function evaluated at nodes
         const IntegrationRule &nodes = fq.GetNodalRule();
         Vector nodal(1 * nbf);  // single component
         for (int k = 0; k < nbf; k++)
         {
            real_t x = nodes.IntPoint(k).x;
            real_t y = nodes.IntPoint(k).y;
            nodal(k) = 2.0 * x + 3.0 * y + 1.0;  // linear function
         }

         Vector quad(1 * nq);
         fq.InterpolateToQuadPoints(1, nodal, quad);

         Vector recovered(1 * nbf);
         fq.GalerkinProject(1, quad, recovered);

         for (int k = 0; k < nbf; k++)
         {
            std::ostringstream msg;
            msg << "p=" << p << " linear roundtrip k=" << k;
            TEST_NEAR(recovered(k), nodal(k), 1e-12, msg.str());
         }
      }

      // --- Quadratic roundtrip (should work at order >= 2) ---
      if (p >= 2)
      {
         const IntegrationRule &nodes = fq.GetNodalRule();
         Vector nodal(1 * nbf);
         for (int k = 0; k < nbf; k++)
         {
            real_t x = nodes.IntPoint(k).x;
            real_t y = nodes.IntPoint(k).y;
            nodal(k) = x*x + x*y + y*y + x + y + 1.0;  // quadratic
         }

         Vector quad(1 * nq);
         fq.InterpolateToQuadPoints(1, nodal, quad);

         Vector recovered(1 * nbf);
         fq.GalerkinProject(1, quad, recovered);

         for (int k = 0; k < nbf; k++)
         {
            std::ostringstream msg;
            msg << "p=" << p << " quadratic roundtrip k=" << k;
            TEST_NEAR(recovered(k), nodal(k), 1e-11, msg.str());
         }
      }
   }

   std::cout << "  Done." << std::endl;
}

// ============================================================================
// Test 7: Order 0 GalerkinProject gives weighted average
// ============================================================================
void TestOrder0ProjectIsAverage()
{
   std::cout << "TestOrder0ProjectIsAverage..." << std::endl;

   FaceQuadrature fq(0, 1);
   int nq = fq.NumQuadPoints();
   const IntegrationRule &ir = fq.GetQuadRule();

   // Set quad values to known function
   Vector quad_vals(nq);
   for (int q = 0; q < nq; q++)
   {
      real_t x = ir.IntPoint(q).x;
      real_t y = ir.IntPoint(q).y;
      quad_vals(q) = 1.0 + 2.0*x + 3.0*y;  // some varying function
   }

   // Project to single DOF
   Vector nodal(1);
   fq.GalerkinProject(1, quad_vals, nodal);

   // Expected: weighted average = (Σ w_q * f_q) / (Σ w_q)
   real_t sum_wf = 0.0, sum_w = 0.0;
   for (int q = 0; q < nq; q++)
   {
      sum_wf += ir.IntPoint(q).weight * quad_vals(q);
      sum_w += ir.IntPoint(q).weight;
   }
   real_t expected = sum_wf / sum_w;

   TEST_NEAR(nodal(0), expected, 1e-14,
             "order=0 project should equal weighted average");

   // Also verify with multi-component
   int ncomp = 3;
   Vector quad_multi(ncomp * nq);
   for (int c = 0; c < ncomp; c++)
   {
      for (int q = 0; q < nq; q++)
      {
         quad_multi(c * nq + q) = (c + 1.0) * quad_vals(q);
      }
   }

   Vector nodal_multi(ncomp * 1);
   fq.GalerkinProject(ncomp, quad_multi, nodal_multi);

   for (int c = 0; c < ncomp; c++)
   {
      std::ostringstream msg;
      msg << "multi-comp c=" << c << " average";
      TEST_NEAR(nodal_multi(c), (c + 1.0) * expected, 1e-14, msg.str());
   }

   std::cout << "  Done." << std::endl;
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char *argv[])
{
   std::cout << "=== FaceQuadrature Unit Tests ===" << std::endl;

   TestConstruction();
   TestPartitionOfUnity();
   TestOrder0IsConstant();
   TestMassMatrixSPD();
   TestOrder1MassMatrixValues();
   TestRoundtrip();
   TestOrder0ProjectIsAverage();

   std::cout << "\n=== Results: " << num_passed << "/" << num_tests
             << " passed";
   if (num_failed > 0)
   {
      std::cout << " (" << num_failed << " FAILED)";
   }
   std::cout << " ===" << std::endl;

   return (num_failed > 0) ? 1 : 0;
}
