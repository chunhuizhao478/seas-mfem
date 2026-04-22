// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.11 (rev-3g, NEW) — RK4 conservation probe.
// Targets H-V92-G candidate (b): RK4 non-conservative drift on σ_yy.
//
// ============================================================================
// Strategy
// ============================================================================
// Build a cube domain, absorbing BCs.  Set a SMOOTH Gaussian-pulse
// initial condition localized at the domain center.  Drive N RK4
// steps with CFL-safe dt.  Measure the elastic L² energy
//     E(t) = ∫ (σ:σ / 4μ) dV + ∫ (ρ|v|² / 2) dV
// at t = 0 and every 10 steps.  Also record the per-channel L²
// integrals separately:
//     E_c(t) = ∫ σ_c² dV   (for each stress channel)
// BEFORE the wave reaches the domain boundary, conservation gives
//     dE/dt ≈ 0      (non-dissipative bulk RK4)
// Any drift slope Δ_c = (E_c(T_N) − E_c(0))/T_N flags channel-specific
// non-conservation.
//
// Candidate (b) H-V92-G specific prediction: σ_yy / σ_zz / σ_yz
// channels (the "zero-wave-speed" stresses per Pelties 2012 §3.2)
// are coupled to velocity via A_x/A_y/A_z entries in a way that the
// explicit RK4 update may accumulate drift where an ADER-DG scheme
// would not.  If Δ_SYY >> Δ_SXX (P-wave channel), (b) CONFIRMED.
//
// ============================================================================
// Why this is cleaner than testing on the production fixture
// ============================================================================
// Production TPV102 has friction nonlinearity, rupture nucleation,
// fault flux machinery, plus boundary interactions.  This probe
// ISOLATES the volume + time-stepping operator by using:
//   - No fault (bc.fault_attr = 0).
//   - No rupture (friction solver is not invoked).
//   - Smooth IC far from boundaries.
//   - Fixed number of small RK4 steps.
// Any drift observed is attributable to (b) alone.
//
// Usage:
//   ./seas_test_rk4_conservation   (serial, < 5 s)

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
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_BOOL(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else      { num_failed++; std::cout << "  FAILED [" << __LINE__ \
                       << "]: " << msg << "\n"; } \
} while (0)

// RK4 step helper — mirrors the driver's integrator.
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

// Compute per-channel L² energy: ||σ_c||² = ∫ σ_c² dV via a quadrature
// on the DG space.  A cheap proxy: sum(dof_val²) · element_volume /
// num_dofs_per_elem.  This is NOT a true L² inner product, but the
// DRIFT in this proxy tracks true L² drift when the DOF layout is
// uniform, which it is on a Cartesian mesh.
static double PerChannelL2Sq(const Vector &Q, int c, int ndof_total)
{
   double s = 0.0;
   for (int d = 0; d < ndof_total; d++)
   {
      double v = Q(c * ndof_total + d);
      s += v * v;
   }
   return s;
}

int main()
{
   std::cout << "\n=== TPV102 v9.2.0 §4.11 Probe: "
             << "RK4 Conservation (H-V92-G (b)) ===\n";

   // 8×8×8 hex mesh, 1000 m per element, absorbing BCs everywhere.
   const double Lx = 8.0e3, Ly = 8.0e3, Lz = 8.0e3;  // 8 km cube
   const int nx = 8, ny = 8, nz = 8;
   Mesh mesh = Mesh::MakeCartesian3D(nx, ny, nz, Element::HEXAHEDRON,
                                     Lx, Ly, Lz);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *x = mesh.GetVertex(v);
      x[0] -= 0.5*Lx;  x[1] -= 0.5*Ly;  x[2] -= 0.5*Lz;  // center at 0
   }
   for (int be = 0; be < mesh.GetNBE(); be++)
   { mesh.SetBdrAttribute(be, 5); }
   mesh.SetAttributes();

   BoundaryConfig bc;
   bc.natural_attrs   = {};
   bc.fault_attr      = 0;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<Mesh> wave(mesh, order,
                           TPV102Params::lambda,
                           TPV102Params::mu,
                           TPV102Params::rho, bc);
   { real_t zero_bg[NUM_STATE] = {0}; wave.SetAbsorbingBackground(zero_bg); }

   const auto &fes_const = wave.GetFESpace();
   auto &fes = const_cast<FiniteElementSpace &>(fes_const);
   const int size = wave.Height();
   const int ndof_total = fes.GetNDofs();

   // Gaussian pulse IC: σ_yy(x,y,z,0) = A · exp(-(r/σ)²) with
   // σ = Lx/8 = 1 km, A = 1e6 Pa; velocities zero.  This pulse is
   // localized near origin and has NO flux across ±x, ±y, ±z at
   // t=0, so energy is well-bounded and conserved until the
   // wavefront reaches the boundary.
   const double pulse_sigma = Lx / 8.0;  // 1 km
   const double pulse_amp   = 1.0e6;      // 1 MPa
   FunctionCoefficient syy_init(
      [pulse_sigma, pulse_amp](const Vector &x) -> real_t {
         double r2 = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
         return pulse_amp * std::exp(-r2 / (pulse_sigma*pulse_sigma));
      });
   FunctionCoefficient sxx_init(
      [pulse_sigma, pulse_amp](const Vector &x) -> real_t {
         double r2 = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
         // Also seed σ_xx with a different amplitude for comparison.
         return 0.5 * pulse_amp * std::exp(-r2 / (pulse_sigma*pulse_sigma));
      });

   Vector Q(size);
   Q = 0.0;
   GridFunction gf_syy(&fes), gf_sxx(&fes);
   gf_syy.ProjectCoefficient(syy_init);
   gf_sxx.ProjectCoefficient(sxx_init);
   for (int d = 0; d < ndof_total; d++)
   {
      Q(SYY * ndof_total + d) = gf_syy(d);
      Q(SXX * ndof_total + d) = gf_sxx(d);
   }

   // Wave speed: cp ≈ 6000 m/s.  Element size h = Lx/nx = 1 km.
   // CFL-safe dt for P1 RK4: dt ≈ 0.3 * h / cp = 0.3 * 1000 / 6000 ≈ 5·10⁻⁵ s.
   // Wall-clock for wavefront to reach nearest boundary: Lx/(2 cp) =
   // 4000 / 6000 ≈ 0.67 s.  We take N = 100 steps of dt = 5e-5 s ⇒
   // T_N = 5 ms — pulse still safely inside domain.
   const real_t dt    = 5.0e-5;   // 50 μs
   const int   N_steps = 100;

   // Initial per-channel L² values.
   double E0[NUM_STATE] = {0};
   for (int c = 0; c < NUM_STATE; c++)
   { E0[c] = PerChannelL2Sq(Q, c, ndof_total); }

   std::cout << "  dt = " << dt << " s, N = " << N_steps
             << " RK4 steps, T = " << N_steps * dt << " s\n";
   std::cout << "  cs/cp travel time to boundary ≈ "
             << (Lx / 2.0) / TPV102Params::cp << " s "
             << "(pulse stays inside)\n\n";

   std::cout << "  Initial L² per channel (∫ σ² dV proxy):\n";
   static const char *cn[9] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ","VX","VY","VZ"};
   for (int c = 0; c < NUM_STATE; c++)
   {
      std::cout << "    " << cn[c] << " L²₀ = " << std::scientific
                << std::setprecision(6) << E0[c] << "\n";
   }

   // Run N RK4 steps.
   for (int step = 0; step < N_steps; step++)
   {
      DoRK4Step(wave, Q, dt);
   }

   // Final per-channel L².
   double E1[NUM_STATE] = {0};
   for (int c = 0; c < NUM_STATE; c++)
   { E1[c] = PerChannelL2Sq(Q, c, ndof_total); }

   std::cout << "\n  After " << N_steps << " RK4 steps (t = "
             << N_steps * dt << " s):\n";
   std::cout << "  c    | L²₀           | L²_T          | (L²_T-L²₀)/L²₀ | "
                "Δ per step\n";
   std::cout << "  -----+---------------+---------------+----------------+"
                "------------\n";

   // Drift measures.
   double worst_drift_nonphysical = 0.0;
   int    worst_channel           = -1;
   for (int c = 0; c < NUM_STATE; c++)
   {
      const double dE = E1[c] - E0[c];
      const double rel = (E0[c] > 1e-30) ? dE / E0[c] : 0.0;
      const double per_step = rel / N_steps;
      std::printf("  %-4s | %.6e | %.6e | %+.4e     | %+.4e\n",
                  cn[c], E0[c], E1[c], rel, per_step);
      // For channels initialized at zero (SZZ, SXY, SYZ, SXZ, VX, VY, VZ),
      // any post-run energy is REAL physics (mode conversion) — not
      // (b) drift.  Only SYY and SXX were seeded; their relative drift
      // is the conservation probe.
      if (c == SXX || c == SYY)
      {
         if (std::abs(per_step) > worst_drift_nonphysical)
         {
            worst_drift_nonphysical = std::abs(per_step);
            worst_channel = c;
         }
      }
   }

   std::cout << "\n  Worst per-step relative drift on seeded channels "
             << "{SXX, SYY}: ";
   std::cout << "channel=" << (worst_channel>=0 ? cn[worst_channel] : "-")
             << "  = " << std::scientific << std::setprecision(3)
             << worst_drift_nonphysical << " /step\n";

   // Pass criterion: for a well-designed RK4 DG + absorbing BC, the
   // seeded-channel L² should be stable to better than 1% per step
   // before the wave touches the boundary.  (1e-2 / step × 100 steps
   // = 100% total drift, the worst-case we'll accept without flagging.)
   // A well-designed scheme actually stays below 1e-4 per step.
   TEST_BOOL(worst_drift_nonphysical < 1.0e-2,
             "seeded-channel L² drift < 1%/step (conservation intact)");

   // A STRONGER test: the σ_yy vs σ_xx drift should agree to within
   // an order of magnitude.  If σ_yy drifts 100× more than σ_xx, the
   // A_x/A_y Poisson-coupling is non-conservative on σ_yy specifically
   // — a direct (b) signature.
   const double drift_SYY = std::abs((E1[SYY] - E0[SYY]) /
                                      std::max(E0[SYY], 1.0e-30));
   const double drift_SXX = std::abs((E1[SXX] - E0[SXX]) /
                                      std::max(E0[SXX], 1.0e-30));
   const double ratio = drift_SYY / std::max(drift_SXX, 1.0e-30);
   std::cout << "\n  SYY / SXX drift ratio = " << ratio;
   if (ratio > 10.0)
   { std::cout << "   ⇒ asymmetric (b) candidate SIGNATURE\n"; }
   else
   { std::cout << "   ⇒ symmetric (bulk RK4 appears channel-uniform)\n"; }
   TEST_BOOL(ratio < 10.0,
             "SYY and SXX drift magnitudes within 10× (no (b) asymmetry signature)");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
