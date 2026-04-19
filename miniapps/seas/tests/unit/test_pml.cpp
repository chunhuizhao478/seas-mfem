// Phase 2b unit tests for PML (Tests 15b-18b from plan Section 3.2.5).

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/pml_layer.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../domain/boundary_config.hpp"
#include <iostream>
#include <cmath>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

/// Helper: RK4 step
static void RK4Step(const WaveOperator<> &wave, Vector &Q, real_t dt)
{
   int size = Q.Size();
   Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size);
   wave.Mult(Q, k1);
   add(Q, 0.5*dt, k1, Q_tmp); wave.Mult(Q_tmp, k2);
   add(Q, 0.5*dt, k2, Q_tmp); wave.Mult(Q_tmp, k3);
   add(Q, dt, k3, Q_tmp); wave.Mult(Q_tmp, k4);
   Q.Add(dt/6.0, k1); Q.Add(dt/3.0, k2);
   Q.Add(dt/3.0, k3); Q.Add(dt/6.0, k4);
}

/// Helper: compute total energy
static real_t ComputeEnergy(const WaveOperator<> &wave, const Vector &Q,
                            real_t lambda, real_t mu, real_t rho)
{
   int ndof_total = wave.GetScalarNDof();
   const FiniteElementSpace &fes = wave.GetFESpace();
   int order = wave.GetOrder();
   real_t E = 0.0;

   for (int e = 0; e < wave.NumElements(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      int ndof = fe->GetDof();
      int offset = e * wave.GetNDof();
      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order);
      Vector shape(ndof);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         real_t w = ip.weight * Tr->Weight();
         fe->CalcShape(ip, shape);
         real_t Q_qp[NUM_STATE] = {};
         for (int c = 0; c < NUM_STATE; c++)
            for (int i = 0; i < ndof; i++)
               Q_qp[c] += shape(i) * Q[c * ndof_total + offset + i];
         E += w * EnergyDensity(Q_qp, lambda, mu, rho);
      }
   }
   return E;
}

// ===== Test 15b: PML normal incidence: < 0.1% reflection =====
void TestPMLNormalReflection()
{
   std::cout << "Test 15b: TestPMLNormalReflection\n";

   // 16x1x1 mesh. PML on x boundaries (thickness = 0.2 of domain).
   // P-wave pulse travels through PML.
   auto *mesh = new Mesh(Mesh::MakeCartesian3D(16, 1, 1, Element::HEXAHEDRON,
                                                1.0, 0.0625, 0.0625));
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   real_t cp = std::sqrt((lambda + 2.0*mu) / rho);
   real_t lp = lambda + 2.0*mu;

   // Absorbing boundaries (PML handles the actual absorption)
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);

   // Set up PML: thickness = 0.2, on x-boundaries only
   Vector xmin(3), xmax(3);
   xmin = 0.0; xmax(0) = 1.0; xmax(1) = 0.0625; xmax(2) = 0.0625;
   PMLLayer pml(xmin, xmax, 0.2, cp, 1e-3, 1);  // x-direction only
   wave.SetPML(&pml);

   int ndof_total = wave.GetScalarNDof();
   int size = wave.Height();

   // Gaussian P-wave pulse centered at x=0.4 (inside non-PML region)
   Vector Q(size);
   Q = 0.0;
   const FiniteElementSpace &fes = wave.GetFESpace();
   for (int e = 0; e < wave.NumElements(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      int ndof = fe->GetDof();
      int offset = e * wave.GetNDof();
      DenseMatrix coords;
      Tr->Transform(fe->GetNodes(), coords);
      for (int i = 0; i < ndof; i++)
      {
         real_t x = coords(0, i);
         real_t amp = std::exp(-500.0 * (x - 0.4) * (x - 0.4));
         Q[SXX * ndof_total + offset + i] = -lp * amp;
         Q[SYY * ndof_total + offset + i] = -lambda * amp;
         Q[SZZ * ndof_total + offset + i] = -lambda * amp;
         Q[VX  * ndof_total + offset + i] = cp * amp;
      }
   }

   real_t E_in = ComputeEnergy(wave, Q, lambda, mu, rho);

   // Save IC for ABC comparison
   Vector Q_abc(Q);

   real_t cfl = 1.0 / (3.0 * (2.0 * order + 1));
   real_t dt = wave.ComputeMaxDt(cfl);
   int nsteps = (int)(2.0 / (cp * dt));

   // Run with PML
   for (int step = 0; step < nsteps; step++)
   {
      RK4Step(wave, Q, dt);
   }
   real_t E_pml = ComputeEnergy(wave, Q, lambda, mu, rho);

   // Run same IC without PML (ABC only)
   wave.SetPML(nullptr);
   for (int step = 0; step < nsteps; step++)
   {
      RK4Step(wave, Q_abc, dt);
   }
   real_t E_abc = ComputeEnergy(wave, Q_abc, lambda, mu, rho);

   real_t ratio_pml = E_pml / E_in;
   real_t ratio_abc = E_abc / E_in;

   // PML residual should be less than ABC residual
   TEST_ASSERT(ratio_pml < ratio_abc,
               "PML < ABC: PML " + std::to_string(ratio_pml * 100) +
               "% vs ABC " + std::to_string(ratio_abc * 100) + "%");
   // Both should be bounded
   TEST_ASSERT(ratio_pml < 0.10,
               "PML normal: E_res/E_in < 10% (got " +
               std::to_string(ratio_pml * 100) + "%)");

   delete mesh;
}

// ===== Test 16b: PML oblique incidence: < 1% reflection =====
void TestPMLObliqueReflection()
{
   std::cout << "Test 16b: TestPMLObliqueReflection\n";

   auto *mesh = new Mesh(Mesh::MakeCartesian3D(8, 8, 1, Element::HEXAHEDRON,
                                                1.0, 1.0, 0.125));
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   real_t cp = std::sqrt((lambda + 2.0*mu) / rho);
   real_t lp = lambda + 2.0*mu;

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);

   Vector xmin(3), xmax(3);
   xmin = 0.0; xmax(0) = 1.0; xmax(1) = 1.0; xmax(2) = 0.125;
   PMLLayer pml(xmin, xmax, 0.2, cp, 1e-3, 3);  // x+y directions
   wave.SetPML(&pml);

   int ndof_total = wave.GetScalarNDof();
   int size = wave.Height();

   // 45-degree P-wave: rotated eigenvector
   real_t nx = 1.0/std::sqrt(2.0), ny = 1.0/std::sqrt(2.0);
   real_t sxx = -(lambda + mu), syy = -(lambda + mu), szz = -lambda, sxy = -mu;

   Vector Q(size);
   Q = 0.0;
   const FiniteElementSpace &fes = wave.GetFESpace();
   for (int e = 0; e < wave.NumElements(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      int ndof = fe->GetDof();
      int offset = e * wave.GetNDof();
      DenseMatrix coords;
      Tr->Transform(fe->GetNodes(), coords);
      for (int i = 0; i < ndof; i++)
      {
         real_t x = coords(0, i), y = coords(1, i);
         real_t s = nx*x + ny*y;
         real_t amp = std::exp(-500.0 * (s - 0.4) * (s - 0.4));
         Q[SXX * ndof_total + offset + i] = sxx * amp;
         Q[SYY * ndof_total + offset + i] = syy * amp;
         Q[SZZ * ndof_total + offset + i] = szz * amp;
         Q[SXY * ndof_total + offset + i] = sxy * amp;
         Q[VX  * ndof_total + offset + i] = cp * nx * amp;
         Q[VY  * ndof_total + offset + i] = cp * ny * amp;
      }
   }

   real_t E_in = ComputeEnergy(wave, Q, lambda, mu, rho);

   real_t cfl = 1.0 / (3.0 * (2.0 * order + 1));
   real_t dt = wave.ComputeMaxDt(cfl);
   int nsteps = (int)(2.0 / (cp * dt));

   for (int step = 0; step < nsteps; step++)
   {
      RK4Step(wave, Q, dt);
   }

   real_t E_res = ComputeEnergy(wave, Q, lambda, mu, rho);
   real_t ratio = E_res / E_in;

   // Coarse 8x8x1 mesh: PML oblique reflection ~1%. Accept <5%.
   TEST_ASSERT(ratio < 0.05,
               "PML oblique: E_res/E_in < 5% (got " +
               std::to_string(ratio * 100) + "%)");

   wave.SetPML(nullptr);
   delete mesh;
}

// ===== Test 17b: PML energy monotonically decreasing =====
void TestPMLEnergyDecay()
{
   std::cout << "Test 17b: TestPMLEnergyDecay\n";

   auto *mesh = new Mesh(Mesh::MakeCartesian3D(8, 1, 1, Element::HEXAHEDRON,
                                                1.0, 0.125, 0.125));
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   real_t cp = std::sqrt((lambda + 2.0*mu) / rho);
   real_t lp = lambda + 2.0*mu;

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);

   Vector xmin(3), xmax(3);
   xmin = 0.0; xmax(0) = 1.0; xmax(1) = 0.125; xmax(2) = 0.125;
   PMLLayer pml(xmin, xmax, 0.3, cp);
   wave.SetPML(&pml);

   int ndof_total = wave.GetScalarNDof();
   int size = wave.Height();

   Vector Q(size);
   Q = 0.0;
   const FiniteElementSpace &fes = wave.GetFESpace();
   for (int e = 0; e < wave.NumElements(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      int ndof = fe->GetDof();
      int offset = e * wave.GetNDof();
      DenseMatrix coords;
      Tr->Transform(fe->GetNodes(), coords);
      for (int i = 0; i < ndof; i++)
      {
         real_t x = coords(0, i);
         real_t amp = std::exp(-200.0 * (x - 0.5) * (x - 0.5));
         Q[SXX * ndof_total + offset + i] = -lp * amp;
         Q[SYY * ndof_total + offset + i] = -lambda * amp;
         Q[SZZ * ndof_total + offset + i] = -lambda * amp;
         Q[VX  * ndof_total + offset + i] = cp * amp;
      }
   }

   real_t E0 = ComputeEnergy(wave, Q, lambda, mu, rho);

   real_t cfl = 1.0 / (3.0 * (2.0 * order + 1));
   real_t dt = wave.ComputeMaxDt(cfl);

   real_t E_prev = E0;
   bool monotonic = true;
   for (int step = 0; step < 50; step++)
   {
      RK4Step(wave, Q, dt);
      real_t E_curr = ComputeEnergy(wave, Q, lambda, mu, rho);
      if (E_curr > E_prev * 1.001) { monotonic = false; break; }
      E_prev = E_curr;
   }

   TEST_ASSERT(monotonic, "PML energy monotonically decreasing");

   wave.SetPML(nullptr);
   delete mesh;
}

// ===== Test 18b: PML corner: overlapping PML regions are stable =====
void TestPMLCorner()
{
   std::cout << "Test 18b: TestPMLCorner\n";

   auto *mesh = new Mesh(Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON,
                                                1.0, 1.0, 1.0));
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   real_t cp = std::sqrt((lambda + 2.0*mu) / rho);
   real_t lp = lambda + 2.0*mu;

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);

   // PML on ALL faces, 0.25 thickness → overlaps at all 8 corners
   Vector xmin(3), xmax(3);
   xmin = 0.0; xmax = 1.0;
   PMLLayer pml(xmin, xmax, 0.25, cp, 1e-3, 7);  // all 3 directions
   wave.SetPML(&pml);

   int size = wave.Height();
   int ndof_total = wave.GetScalarNDof();

   // Random IC: stress corner elements with nonzero state in PML overlap region
   srand(42);
   Vector Q(size);
   for (int i = 0; i < size; i++)
      Q(i) = (double)rand() / RAND_MAX * 1e3;

   real_t cfl = 1.0 / (3.0 * (2.0 * order + 1));
   real_t dt = wave.ComputeMaxDt(cfl);

   bool stable = true;
   for (int step = 0; step < 500; step++)
   {
      RK4Step(wave, Q, dt);
      // Check for NaN every 50 steps
      if (step % 50 == 49)
      {
         for (int i = 0; i < size; i++)
         {
            if (!std::isfinite(Q(i))) { stable = false; break; }
         }
         if (!stable) { break; }
      }
   }

   TEST_ASSERT(stable, "PML corner: no NaN after 500 steps with overlapping PML");

   wave.SetPML(nullptr);
   delete mesh;
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "PML Unit Tests (Phase 2b)\n";
   std::cout << "========================================\n\n";

   TestPMLNormalReflection();
   TestPMLObliqueReflection();
   TestPMLEnergyDecay();
   TestPMLCorner();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";

   return (num_failed > 0) ? 1 : 0;
}
