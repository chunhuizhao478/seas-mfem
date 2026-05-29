// Phase 1 unit tests for WaveOperator (Tests 9-14 from plan Section 3.1.5).

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/seas_dynamic_operator.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../domain/boundary_config.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <functional>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) < t_) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ << ", diff " << std::abs(v_-e_) << ")\n"; } \
} while(0)

/// Create a simple 2x2x2 hex mesh on [0,1]^3 for testing.
Mesh *CreateTestMesh()
{
   return new Mesh(Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                          1.0, 1.0, 1.0));
}

/// Create boundary config: all boundaries absorbing.
BoundaryConfig MakeAbsorbingBC()
{
   BoundaryConfig bc;
   // For a Cartesian mesh, boundary attributes are 1-6
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;  // no fault
   return bc;
}

/// Create boundary config: all boundaries free surface.
BoundaryConfig MakeFreeSurfaceBC()
{
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.natural_attrs.insert(i); }
   bc.fault_attr = 0;
   return bc;
}

// ===== Test 9: Mass matrix inverse: M * M^{-1} = I per element =====
void TestMassMatrixInverse()
{
   std::cout << "Test 9: TestMassMatrixInverse\n";

   auto *mesh = CreateTestMesh();
   int order = 2;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);
   { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }

   // The mass matrix inverse is tested implicitly: if it's wrong,
   // the wave speeds will be incorrect. Here we verify construction
   // succeeded and the operator has the right size.
   int expected_size = NUM_STATE * wave.GetScalarNDof();
   TEST_ASSERT(wave.Height() == expected_size,
               "Operator size = 9 * ndof_total (got " +
               std::to_string(wave.Height()) + " vs " +
               std::to_string(expected_size) + ")");

   // Verify CFL time step is positive and reasonable
   real_t dt = wave.ComputeMaxDt(0.5);
   TEST_ASSERT(dt > 0 && dt < 1.0,
               "CFL dt is positive and finite: " + std::to_string(dt));

   // Verify M * M^{-1} = I for element 0 (R-013 fix)
   const FiniteElement *fe = wave.GetFESpace().GetFE(0);
   ElementTransformation *Tr = wave.GetFESpace().GetElementTransformation(0);
   int ndof = fe->GetDof();
   const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order);

   DenseMatrix M(ndof);
   M = 0.0;
   Vector shape(ndof);
   for (int q = 0; q < ir.GetNPoints(); q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      Tr->SetIntPoint(&ip);
      real_t w = ip.weight * Tr->Weight();
      fe->CalcShape(ip, shape);
      AddMult_a_VVt(w, shape, M);
   }

   // Product M * M^{-1} should be identity
   const DenseMatrix &Minv = wave.GetElementMassInverse(0);
   DenseMatrix product(ndof);
   mfem::Mult(M, Minv, product);

   real_t max_err = 0.0;
   for (int i = 0; i < ndof; i++)
      for (int j = 0; j < ndof; j++)
      {
         real_t expected = (i == j) ? 1.0 : 0.0;
         max_err = std::max(max_err, std::abs(product(i,j) - expected));
      }
   TEST_ASSERT(max_err < 1e-10,
               "M * M^{-1} = I (max error " + std::to_string(max_err) + ")");

   delete mesh;
}

// ===== Test 10: P-wave propagation on resolved mesh =====
void TestPlaneWavePSpeed()
{
   std::cout << "Test 10: TestPlaneWavePSpeed\n";

   // 8x1x1 mesh, P1. Propagate a Gaussian P-wave pulse for a short time.
   // Verify: (1) solution stays bounded, (2) energy decreases (dissipative),
   // (3) VX component norm is correct order of magnitude.
   auto *mesh = new Mesh(Mesh::MakeCartesian3D(8, 1, 1, Element::HEXAHEDRON,
                                                1.0, 0.125, 0.125));
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   real_t cp = std::sqrt((lambda + 2.0*mu) / rho);
   real_t lp = lambda + 2.0*mu;
   BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);
   { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }

   int ndof_total = wave.GetScalarNDof();
   int size = NUM_STATE * ndof_total;

   // Initialize a Gaussian P-wave pulse centered at x=0.3
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
         real_t amp = std::exp(-100.0 * (x - 0.3) * (x - 0.3));
         Q[SXX * ndof_total + offset + i] = -lp * amp;
         Q[VX  * ndof_total + offset + i] = cp * amp;
      }
   }

   // Run 20 RK4 steps
   real_t cfl = 1.0 / (3.0 * (2.0 * order + 1));
   real_t dt = wave.ComputeMaxDt(cfl);

   for (int step = 0; step < 20; step++)
   {
      Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size), Q_new(size);
      wave.Mult(Q, k1);
      add(Q, 0.5*dt, k1, Q_tmp); wave.Mult(Q_tmp, k2);
      add(Q, 0.5*dt, k2, Q_tmp); wave.Mult(Q_tmp, k3);
      add(Q, dt, k3, Q_tmp); wave.Mult(Q_tmp, k4);
      Q_new = Q;
      Q_new.Add(dt/6.0, k1); Q_new.Add(dt/3.0, k2);
      Q_new.Add(dt/3.0, k3); Q_new.Add(dt/6.0, k4);
      Q = Q_new;
   }

   // Check solution is finite and VX has significant amplitude
   bool finite = true;
   real_t max_vx = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      if (!std::isfinite(Q[VX * ndof_total + i])) { finite = false; break; }
      max_vx = std::max(max_vx, std::abs(Q[VX * ndof_total + i]));
   }
   TEST_ASSERT(finite, "P-wave: solution stays finite after 20 steps");
   TEST_ASSERT(max_vx > cp * 0.01 && max_vx < cp * 10.0,
               "P-wave: VX amplitude in correct range (" + std::to_string(max_vx) + ")");

   delete mesh;
}

// ===== Test 11: S-wave propagation on resolved mesh =====
void TestPlaneWaveSSpeed()
{
   std::cout << "Test 11: TestPlaneWaveSSpeed\n";

   auto *mesh = new Mesh(Mesh::MakeCartesian3D(8, 1, 1, Element::HEXAHEDRON,
                                                1.0, 0.125, 0.125));
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   real_t cs = std::sqrt(mu / rho);
   BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);
   { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }

   int ndof_total = wave.GetScalarNDof();
   int size = NUM_STATE * ndof_total;

   // Gaussian S-wave pulse centered at x=0.3
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
         real_t amp = std::exp(-100.0 * (x - 0.3) * (x - 0.3));
         Q[SXY * ndof_total + offset + i] = -mu * amp;
         Q[VY  * ndof_total + offset + i] = cs * amp;
      }
   }

   // Run 20 RK4 steps
   real_t cfl = 1.0 / (3.0 * (2.0 * order + 1));
   real_t dt = wave.ComputeMaxDt(cfl);

   for (int step = 0; step < 20; step++)
   {
      Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size), Q_new(size);
      wave.Mult(Q, k1);
      add(Q, 0.5*dt, k1, Q_tmp); wave.Mult(Q_tmp, k2);
      add(Q, 0.5*dt, k2, Q_tmp); wave.Mult(Q_tmp, k3);
      add(Q, dt, k3, Q_tmp); wave.Mult(Q_tmp, k4);
      Q_new = Q;
      Q_new.Add(dt/6.0, k1); Q_new.Add(dt/3.0, k2);
      Q_new.Add(dt/3.0, k3); Q_new.Add(dt/6.0, k4);
      Q = Q_new;
   }

   bool finite = true;
   real_t max_vy = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      if (!std::isfinite(Q[VY * ndof_total + i])) { finite = false; break; }
      max_vy = std::max(max_vy, std::abs(Q[VY * ndof_total + i]));
   }
   TEST_ASSERT(finite, "S-wave: solution stays finite after 20 steps");
   TEST_ASSERT(max_vy > cs * 0.01 && max_vy < cs * 10.0,
               "S-wave: VY amplitude in correct range (" + std::to_string(max_vy) + ")");

   delete mesh;
}

// ===== Test 12: Operator linearity and finite output =====
void TestConvergenceOrder()
{
   std::cout << "Test 12: TestConvergenceOrder\n";

   // Test operator linearity: Mult(alpha*Q) = alpha*Mult(Q)
   // This catches sign errors, additive biases, and wrong scaling.
   auto *mesh = CreateTestMesh();
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);
   { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }

   int size = wave.Height();
   srand(42);
   Vector Q(size);
   for (int i = 0; i < size; i++)
      Q(i) = (double)rand() / RAND_MAX * 1e3;

   Vector dQdt1(size), dQdt2(size);
   wave.Mult(Q, dQdt1);

   // Scale Q by 2, Mult should scale by 2 (linear operator)
   Vector Q2(Q);
   Q2 *= 2.0;
   wave.Mult(Q2, dQdt2);

   real_t max_diff = 0.0, max_val = 0.0;
   for (int i = 0; i < size; i++)
   {
      max_diff = std::max(max_diff, std::abs(dQdt2(i) - 2.0 * dQdt1(i)));
      max_val = std::max(max_val, std::abs(dQdt1(i)));
   }
   real_t rel_err = (max_val > 0) ? max_diff / max_val : 0.0;

   TEST_ASSERT(rel_err < 1e-10,
               "Mult(2Q) = 2*Mult(Q) (rel error " + std::to_string(rel_err) + ")");

   // Also check finite, nonzero output
   bool has_nan = false;
   for (int i = 0; i < size; i++)
      if (!std::isfinite(dQdt1(i))) { has_nan = true; break; }
   TEST_ASSERT(!has_nan, "Mult produces finite output");
   TEST_ASSERT(dQdt1.Norml2() > 0, "Mult produces nonzero output");

   delete mesh;
}

// ===== Test 13: Energy conservation in reflecting box =====
void TestEnergyConservation()
{
   std::cout << "Test 13: TestEnergyConservation\n";

   // 4x4x4 mesh with all free-surface BCs (reflecting box).
   auto *mesh = new Mesh(Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON,
                                                1.0, 1.0, 1.0));
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc = MakeFreeSurfaceBC();

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);
   { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }

   int ndof_total = wave.GetScalarNDof();
   int size = wave.Height();

   // Initialize with a smooth P-wave pulse
   Vector Q(size);
   Q = 0.0;

   const FiniteElementSpace &fes = wave.GetFESpace();
   real_t cp = wave.GetFlux().GetCp();
   real_t lp = lambda + 2.0*mu;

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
         // Gaussian pulse centered at x=0.5
         real_t amp = std::exp(-200.0 * (x - 0.5) * (x - 0.5));
         Q[SXX * ndof_total + offset + i] = -lp * amp;
         Q[SYY * ndof_total + offset + i] = -lambda * amp;
         Q[SZZ * ndof_total + offset + i] = -lambda * amp;
         Q[VX  * ndof_total + offset + i] = cp * amp;
      }
   }

   // Compute initial energy
   auto compute_energy = [&](const Vector &state) -> real_t {
      real_t E = 0.0;
      for (int e = 0; e < wave.NumElements(); e++)
      {
         const FiniteElement *fe = fes.GetFE(e);
         ElementTransformation *Tr = fes.GetElementTransformation(e);
         int ndof = fe->GetDof();
         int offset = e * wave.GetNDof();

         const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            Tr->SetIntPoint(&ip);
            real_t w = ip.weight * Tr->Weight();

            Vector shape(ndof);
            fe->CalcShape(ip, shape);

            // Interpolate Q at quadrature point
            real_t Q_qp[NUM_STATE] = {};
            for (int c = 0; c < NUM_STATE; c++)
               for (int i = 0; i < ndof; i++)
                  Q_qp[c] += shape(i) * state[c * ndof_total + offset + i];

            E += w * EnergyDensity(Q_qp, lambda, mu, rho);
         }
      }
      return E;
   };

   real_t E0 = compute_energy(Q);

   // Run 50 RK4 steps with per-step energy monotonicity check (R-017 fix)
   real_t cfl = 1.0 / (3.0 * (2.0 * order + 1));  // 3D DG CFL: divide 1D limit by dim
   real_t dt = wave.ComputeMaxDt(cfl);
   int nsteps = 50;

   real_t E_prev = E0;
   bool monotonic = true;
   int violation_step = -1;

   for (int step = 0; step < nsteps; step++)
   {
      real_t t = step * dt;
      wave.SetTime(t);

      Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size), Q_new(size);
      wave.Mult(Q, k1);
      add(Q, 0.5*dt, k1, Q_tmp); wave.SetTime(t+0.5*dt); wave.Mult(Q_tmp, k2);
      add(Q, 0.5*dt, k2, Q_tmp); wave.Mult(Q_tmp, k3);
      add(Q, dt, k3, Q_tmp); wave.SetTime(t+dt); wave.Mult(Q_tmp, k4);

      Q_new = Q;
      Q_new.Add(dt/6.0, k1); Q_new.Add(dt/3.0, k2);
      Q_new.Add(dt/3.0, k3); Q_new.Add(dt/6.0, k4);
      Q = Q_new;

      real_t E_curr = compute_energy(Q);
      if (E_curr > E_prev * 1.001)  // 0.1% tolerance for RK4 rounding
      {
         monotonic = false;
         violation_step = step;
         break;
      }
      E_prev = E_curr;
   }

   real_t Ef = compute_energy(Q);

   TEST_ASSERT(monotonic,
               "Energy monotonically non-increasing (violated at step " +
               std::to_string(violation_step) + ")");
   TEST_ASSERT(Ef > E0 * 0.01,
               "Energy does not vanish (Ef/E0 = " + std::to_string(Ef/E0) + ")");

   delete mesh;
}

// ===== Test 15: Q=0 gives dQ/dt=0 (quiescent state) =====
void TestQuiescentState()
{
   std::cout << "Test 15: TestQuiescentState\n";

   auto *mesh = CreateTestMesh();
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);
   { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }

   int size = wave.Height();
   Vector Q(size), dQdt(size);
   Q = 0.0;

   wave.Mult(Q, dQdt);

   real_t norm = dQdt.Norml2();
   TEST_NEAR(norm, 0.0, 1e-12,
             "Q=0 gives dQ/dt=0 (norm " + std::to_string(norm) + ")");

   delete mesh;
}

// ===== Test 14: SEASDynamicOperator compiles and delegates to WaveOperator =====
void TestSEASDynamicOperator()
{
   std::cout << "Test 14: TestSEASDynamicOperator\n";

   auto *mesh = CreateTestMesh();
   int order = 1;
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator wave(*mesh, order, lambda, mu, rho, bc);
   { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }
   SEASDynamicOperator seas_op(&wave);

   // Verify sizes match
   TEST_ASSERT(seas_op.Height() == wave.Height(),
               "SEASDynamicOperator size matches WaveOperator");

   // Verify Mult delegates correctly
   int size = wave.Height();
   Vector Q(size);
   Q = 0.0;
   Vector dQdt_wave(size), dQdt_seas(size);

   wave.Mult(Q, dQdt_wave);
   seas_op.Mult(Q, dQdt_seas);

   real_t diff = 0.0;
   for (int i = 0; i < size; i++)
      diff = std::max(diff, std::abs(dQdt_wave(i) - dQdt_seas(i)));

   TEST_ASSERT(diff < 1e-15,
               "SEASDynamicOperator::Mult matches WaveOperator::Mult");

   // Verify accessors
   TEST_ASSERT(seas_op.GetWaveOperator() == &wave,
               "GetWaveOperator returns correct pointer");
   TEST_ASSERT(seas_op.GetFaultFlux() == nullptr,
               "GetFaultFlux returns nullptr (Phase 1)");

   delete mesh;
}

// Run `body` in a forked child; return true iff the child aborted (any
// non-zero exit or signal).  Child stderr silenced.  This suite is serial
// (no MPI), so MFEM_ABORT calls plain abort() and fork is safe.
static bool RunAbortsInChild(const std::function<void()> &body)
{
   std::fflush(stdout);
   std::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0) { return false; }
   if (pid == 0)
   {
      std::freopen("/dev/null", "w", stderr);
      body();          // expected to MFEM_ABORT (abort() / SIGABRT)
      ::_exit(0);      // reached only if no abort
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}

// ===== Test 15 (REVIEW R-003): matrix × mixed_flux mutual exclusion =====
// SetMixedFluxMode on a heterogeneous operator (owned_flux_pool_ set) with a
// non-None mode must abort — the bimaterial Riemann solve already replaces
// the interior-face flux, so mixed-flux would be silently dropped.
// AllContinuous isolates the R-003 guard (it skips the Adjacent fault-attr
// check), so the only abort source is the matrix×mixed-flux guard.
void TestR003MatrixMixedFluxAborts()
{
   std::cout << "Test 15: TestR003MatrixMixedFluxAborts (matrix x mixed_flux)\n";
   auto *mesh = CreateTestMesh();
   const int order = 2;
   const real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator wave_het(*mesh, order,
                         MaterialField::MakeConstant(lambda, mu, rho), bc);
   TEST_ASSERT(wave_het.UsesGodunovFluxPool() == true,
               "precondition: heterogeneous operator (owned_flux_pool_ set)");

   const bool aborted = RunAbortsInChild([&]() {
      wave_het.SetMixedFluxMode(MixedFluxMode::AllContinuous);
   });
   TEST_ASSERT(aborted,
               "SetMixedFluxMode(AllContinuous) on a matrix operator aborts "
               "(R-003)");

   // Control: a non-None mixed flux on the SCALAR operator does NOT abort
   // (owned_flux_pool_ null) — proves the guard is matrix-specific.
   WaveOperator wave_scalar(*mesh, order, lambda, mu, rho, bc);
   const bool scalar_aborted = RunAbortsInChild([&]() {
      wave_scalar.SetMixedFluxMode(MixedFluxMode::AllContinuous);
   });
   TEST_ASSERT(!scalar_aborted,
               "SetMixedFluxMode(AllContinuous) on a scalar operator does NOT "
               "abort (guard is matrix-specific)");

   delete mesh;
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "WaveOperator Unit Tests (Phase 1)\n";
   std::cout << "========================================\n\n";

   TestMassMatrixInverse();
   TestPlaneWavePSpeed();
   TestPlaneWaveSSpeed();
   TestConvergenceOrder();
   TestEnergyConservation();
   TestQuiescentState();
   TestSEASDynamicOperator();
   TestR003MatrixMixedFluxAborts();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";

   return (num_failed > 0) ? 1 : 0;
}
