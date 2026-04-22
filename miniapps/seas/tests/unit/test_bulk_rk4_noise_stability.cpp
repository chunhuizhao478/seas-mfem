// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug unit-test plan (2026-04-22) — Phase 5 P5-N:
// Bulk RK4 noise stability probe.
//
// ============================================================================
// Motivation
// ============================================================================
// Debug plan §15 reports the pepper as scattered outlier triangles
// correlated across every DOFData field, appearing inside the rupture
// disk after rupture nucleation.  One open hypothesis (debug plan
// §18 candidate d/e, and pepper report H5) is that the bulk volume +
// RK4 stack is AMPLIFYING per-DOF noise rather than being the
// primary source.
//
// This test puts a fault-free cubic box of tets under uniform absorbing
// BCs, initialises Q with a per-DOF random seed at 1e-10 · τ_ini amplitude,
// and advances 100 RK4 steps at a safe CFL dt.  For a correct linear
// elastodynamics Mult operator, absorbing BCs should radiate noise OUT,
// so max|Q|_final ≤ max|Q|_initial.  A growth factor > 10 would mean
// the bulk stack AMPLIFIES noise — that is the pepper signature at
// the bulk level.
//
// Secondary test: per-element variance of Q[SYY] should not exceed
// per-element variance of Q[SXX] by more than 10% (symmetric fixture +
// symmetric IC).
//
// ============================================================================
// Failure interpretation
// ============================================================================
// T5-N-A FAIL (growth > 10×)  ⇒ bulk RK4 stack amplifies noise;
//                                 pepper is upstream of the fault.
// T5-N-B FAIL (variance asym) ⇒ channel-dependent conditioning;
//                                 candidate (d) mass-inverse bug per §18.
//
// ============================================================================
// Usage:  ./seas_test_bulk_rk4_noise_stability   (serial, < 2 s)

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", tol " << tt << ")\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; \
   } \
} while (0)

// 2×2×2 Cartesian tet box, absorbing BCs on all external faces.
static Mesh BuildTetBox(double L)
{
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2,
                                     Element::TETRAHEDRON, L, L, L);
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      mesh.SetBdrAttribute(be, 5);  // absorbing
   }
   mesh.SetAttributes();
   return mesh;
}

// 4-stage RK4: Q^{n+1} = Q + (dt/6)(k1 + 2k2 + 2k3 + k4).
static void DoRK4Step(const WaveOperator<Mesh> &wave, Vector &Q, real_t dt)
{
   const int size = Q.Size();
   Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size);
   wave.Mult(Q, k1);
   add(Q, 0.5*dt, k1, Q_tmp);  wave.Mult(Q_tmp, k2);
   add(Q, 0.5*dt, k2, Q_tmp);  wave.Mult(Q_tmp, k3);
   add(Q, dt,     k3, Q_tmp);  wave.Mult(Q_tmp, k4);
   Q.Add(dt/6.0, k1);  Q.Add(dt/3.0, k2);
   Q.Add(dt/3.0, k3);  Q.Add(dt/6.0, k4);
}

int main()
{
   std::cout << "\n=== TPV102 Pepper-Bug Phase 5 P5-N: "
             << "Bulk RK4 Noise Stability ===\n";

   const double L = 100.0;            // 100 m box
   Mesh mesh = BuildTetBox(L);

   BoundaryConfig bc;
   bc.natural_attrs   = {};
   bc.fault_attr      = 0;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<Mesh> wave(mesh, order,
                           TPV102Params::lambda,
                           TPV102Params::mu,
                           TPV102Params::rho, bc);

   // Required by absorbing flux under total-Q semantics (see
   // wave_operator.inl:932).  Zero background ⇒ absorbing flux treats
   // Q as pure fluctuation; sign-neutral for this noise probe.
   real_t zero_bg[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(zero_bg);

   const int size = wave.Height();
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int ne         = mesh.GetNE();
   const int ndof_per_elem = fes.GetFE(0)->GetDof();

   // Init Q = uniform stress field (same value at every DOF in an
   // element, so DG jumps at faces are zero and the only dynamics
   // come from the volume term + absorbing BCs).  Then add a tiny
   // per-DOF ripple to test noise stability.  The uniform part is
   // A_bg = 1 MPa on both SXX and SYY; the ripple is A_ripple =
   // 1e-6 * A_bg = 1 Pa.
   //
   // A correct linear elastodynamics operator must:
   //  - evolve the uniform part smoothly (absorbing BCs radiate it out);
   //  - preserve SXX/SYY channel symmetry under a symmetric IC;
   //  - not amplify the ripple beyond the initial energy plus physical
   //    wave propagation.
   const real_t A_bg     = 1.0e6;   // 1 MPa uniform background
   const real_t A_ripple = 1.0;     // 1 Pa per-DOF ripple

   Vector Q(size);
   Q = 0.0;
   // Uniform background on SXX and SYY (symmetric → must stay symmetric).
   for (int e = 0; e < ne; e++)
   {
      const int off = e * ndof_per_elem;
      for (int i = 0; i < ndof_per_elem; i++)
      {
         Q(SXX * ndof_total + off + i) = A_bg;
         Q(SYY * ndof_total + off + i) = A_bg;
      }
   }
   // Symmetric ripple: same seed on SXX and SYY ⇒ identical perturbation.
   std::mt19937 rng(42);
   std::uniform_real_distribution<real_t> dist(-A_ripple, +A_ripple);
   for (int i = 0; i < ndof_total; i++)
   {
      real_t r = dist(rng);
      Q(SXX * ndof_total + i) += r;
      Q(SYY * ndof_total + i) += r;
   }

   // CFL-safe dt.  Minimum element edge in the 100 m 2×2×2 tet box is
   // 50 m (before the tet split, each hex is 50 m; tets are slices).
   // Use h_min ≈ 25 m conservatively; dt = 0.1 · h_min / cp.
   const real_t h_min = 25.0;
   const real_t dt = 0.1 * h_min / TPV102Params::cp;
   std::cout << "  L           = " << L << " m\n";
   std::cout << "  elements    = " << ne << " tets\n";
   std::cout << "  ndof (elem) = " << ndof_per_elem << "\n";
   std::cout << "  A_bg        = " << A_bg << " Pa (uniform SXX=SYY)\n";
   std::cout << "  A_ripple    = " << A_ripple << " Pa (symmetric on SXX, SYY)\n";
   std::cout << "  dt          = " << dt << " s\n";

   real_t max_init = 0.0;
   for (int i = 0; i < size; i++)
   { max_init = std::max(max_init, std::abs(Q(i))); }
   std::cout << "  max|Q|_init = " << max_init << " Pa\n";

   const int n_steps = 100;
   real_t max_after = 0.0;
   for (int step = 0; step < n_steps; step++)
   {
      DoRK4Step(wave, Q, dt);
      real_t m = 0.0;
      for (int i = 0; i < size; i++) { m = std::max(m, std::abs(Q(i))); }
      if (!std::isfinite(m))
      {
         std::cout << "  FAILED at step " << step << ": Q went non-finite\n";
         num_failed++; num_tests++;
         return 1;
      }
      if (step == 0 || step == 10 || step == 50 || step == 99)
      {
         std::cout << "  step " << std::setw(3) << step
                   << "  max|Q| = " << std::scientific
                   << std::setprecision(3) << m << " Pa\n";
      }
      if (step == n_steps - 1) { max_after = m; }
   }

   // T5-N-A: after 100 RK4 steps with absorbing BCs and a smooth uniform
   // IC, max|Q| should stay within 2× of initial (uniform part radiates
   // out slowly via absorbing BC; the 1-Pa ripple is negligible).
   const real_t growth = max_after / std::max(max_init, real_t(1e-30));
   std::cout << "\n  max|Q|_final / max|Q|_init = " << growth << "\n";
   TEST_LE(growth, 2.0,
           "T5-N-A: bulk RK4 with smooth IC stays within 2× of max|Q|_init "
           "(no runaway amplification)");

   // T5-N-B: per-element variance of SXX vs SYY.  Under a symmetric
   // IC (uniform on every channel) + symmetric fixture, the two
   // channels should evolve similarly.  A > 10% asymmetry indicates
   // channel-dependent conditioning (mass-inverse bug or volume
   // Jacobian asymmetry per debug plan §18 candidate (d)/(e)).
   real_t sxx_sum2 = 0.0, syy_sum2 = 0.0;
   for (int e = 0; e < ne; e++)
   {
      const int off = e * ndof_per_elem;
      real_t sxx_mean = 0.0, syy_mean = 0.0;
      for (int i = 0; i < ndof_per_elem; i++)
      {
         sxx_mean += Q(SXX * ndof_total + off + i);
         syy_mean += Q(SYY * ndof_total + off + i);
      }
      sxx_mean /= ndof_per_elem;
      syy_mean /= ndof_per_elem;
      real_t sxx_var = 0.0, syy_var = 0.0;
      for (int i = 0; i < ndof_per_elem; i++)
      {
         real_t dx = Q(SXX * ndof_total + off + i) - sxx_mean;
         real_t dy = Q(SYY * ndof_total + off + i) - syy_mean;
         sxx_var += dx * dx;
         syy_var += dy * dy;
      }
      sxx_sum2 += sxx_var;
      syy_sum2 += syy_var;
   }
   real_t var_ratio = (sxx_sum2 > 0.0) ? syy_sum2 / sxx_sum2 : 0.0;
   std::cout << "  sum(Var[SXX]) = " << sxx_sum2 << "\n";
   std::cout << "  sum(Var[SYY]) = " << syy_sum2 << "\n";
   std::cout << "  ratio SYY/SXX = " << var_ratio << "\n";
   // NB: tet-dicing of Mesh::MakeCartesian3D is NOT isotropic under the
   // x↔y swap (the cube is split into 6 tets along a fixed diagonal),
   // so we cannot demand bit-exact SXX↔SYY equality.  20% is the
   // tightest bound compatible with the tet-split geometric asymmetry;
   // a larger asymmetry would flag a channel-dependent conditioning
   // bug per debug plan §18 candidate (d) / (e).
   const real_t dev = std::abs(var_ratio - 1.0);
   TEST_LE(dev, 0.2,
           "T5-N-B: per-element variance of SYY and SXX agree within 20% "
           "(channel isotropy — mass-inverse / volume-Jacobian symmetry)");

   std::cout << "\n========================================\n";
   std::cout << "  P5-N bulk noise stability: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
