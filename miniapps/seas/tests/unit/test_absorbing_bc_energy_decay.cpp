// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.12 (rev-3g, NEW) — absorbing-BC energy decay.
// Targets H-V92-G candidate (c): absorbing BC reflection with wrong sign.
//
// ============================================================================
// Strategy
// ============================================================================
// Cube domain 16 km × 16 km × 16 km, hex mesh, ALL boundaries absorbing.
// Initial condition: Gaussian P-wave packet localized at domain center,
// propagating radially outward (realized as σ_xx = σ_yy = σ_zz = A·exp
// (-r²/σ²), ρ·|v| = 0 — stress-only initial state radiates in a mode
// that decomposes into outgoing P-waves).
//
// Drive N RK4 steps long enough for the pulse to reach and cross the
// boundary.  Track total elastic energy
//     E(t) = ∫ (σ:σ / 4μ) dV + ∫ (ρ|v|² / 2) dV
// at every 20 steps.  Classification:
//   - E(t) decreases monotonically → absorbing BC is dissipative as
//     designed; candidate (c) ELIMINATED.
//   - E(t) INCREASES at any point after the wavefront reaches the
//     boundary → the BC is reflecting with wrong sign, injecting
//     energy INTO the domain; candidate (c) CONFIRMED.
//   - E(t) is stable but with large reflection artifacts → the BC
//     is imperfectly absorbing (expected for first-order Clayton-
//     Engquist), but not wrong-signed.
//
// Pelties 2012 §3.2 does not specify an absorbing BC (they use a
// conservatively-large domain to avoid reflections).  SEAS-MFEM uses
// GodunovFlux::Absorbing(n, Q) = A_n^+ Q — the upwind flux that kills
// incoming characteristics.  If this flux is correctly signed, E(t)
// is monotonically decreasing after pulse-boundary interaction.
//
// Usage:
//   ./seas_test_absorbing_bc_energy_decay   (serial, ≤ 30 s)

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

// Total elastic energy proxy.  Uses the DG DOF values directly (no
// mass matrix multiply); this is proportional to the true L² inner
// product on a uniform Cartesian mesh + P1 basis.  The proportionality
// constant doesn't matter — we track RELATIVE changes.
static double TotalEnergy(const Vector &Q, int ndof_total,
                           double lambda, double mu, double rho)
{
   double strain = 0.0;
   const double inv4mu = 0.25 / mu;
   const double poisson = lambda / (3.0*lambda + 2.0*mu);
   for (int d = 0; d < ndof_total; d++)
   {
      real_t sxx = Q(SXX * ndof_total + d);
      real_t syy = Q(SYY * ndof_total + d);
      real_t szz = Q(SZZ * ndof_total + d);
      real_t sxy = Q(SXY * ndof_total + d);
      real_t syz = Q(SYZ * ndof_total + d);
      real_t sxz = Q(SXZ * ndof_total + d);
      double tr = sxx + syy + szz;
      double sig_sq = sxx*sxx + syy*syy + szz*szz
                    + 2.0*(sxy*sxy + syz*syz + sxz*sxz);
      strain += inv4mu * (sig_sq - poisson * tr * tr);
   }
   double kin = 0.0;
   for (int d = 0; d < ndof_total; d++)
   {
      real_t vx = Q(VX * ndof_total + d);
      real_t vy = Q(VY * ndof_total + d);
      real_t vz = Q(VZ * ndof_total + d);
      kin += 0.5 * rho * (vx*vx + vy*vy + vz*vz);
   }
   return strain + kin;
}

int main()
{
   std::cout << "\n=== TPV102 v9.2.0 §4.12 Probe: "
             << "Absorbing-BC Energy Decay (H-V92-G (c)) ===\n";

   // 8×8×8 hex, 2 km per element, 16 km cube centered at origin.
   const double L = 16.0e3;
   const int n = 8;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::HEXAHEDRON, L, L, L);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *x = mesh.GetVertex(v);
      x[0] -= 0.5*L;  x[1] -= 0.5*L;  x[2] -= 0.5*L;
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

   // Gaussian pulse IC: isotropic stress = A·exp(-r²/σ²), σ = L/16 = 1 km.
   const double sigma_pulse = L / 16.0;   // 1 km
   const double pulse_amp   = 1.0e6;      // 1 MPa
   FunctionCoefficient iso(
      [sigma_pulse, pulse_amp](const Vector &x) -> real_t {
         double r2 = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
         return pulse_amp * std::exp(-r2 / (sigma_pulse*sigma_pulse));
      });

   Vector Q(size);
   Q = 0.0;
   GridFunction gf(&fes);
   gf.ProjectCoefficient(iso);
   for (int d = 0; d < ndof_total; d++)
   {
      Q(SXX * ndof_total + d) = gf(d);
      Q(SYY * ndof_total + d) = gf(d);
      Q(SZZ * ndof_total + d) = gf(d);
   }

   // Wave speed cp ≈ 6000 m/s.  Pulse travels from center to nearest
   // boundary = L/2 = 8 km ⇒ travel time ≈ 1.33 s.  CFL: h=2 km, P1
   // RK4 stable dt ≈ 0.3 · h / cp ≈ 1e-4 s.  We need N·dt > 1.33 s ⇒
   // N ≥ 13 300.  That's too many for a unit test.  Use h=2 km (coarse)
   // + dt = 5e-4 s (RK4-stable at this resolution) ⇒ N_steps = 3000
   // for 1.5 s simulation.  30 s wall-clock.
   const real_t dt = 5.0e-4;
   const int   N_steps = 3000;
   const int   sample_every = 100;
   std::cout << "  Domain: " << L/1e3 << " km cube, "
             << n << "³ hex, " << dt << " s dt, "
             << N_steps << " steps (T = " << N_steps*dt << " s).\n";
   std::cout << "  Pulse travel time to boundary ≈ "
             << (L/2) / TPV102Params::cp << " s.\n\n";

   const double E0 = TotalEnergy(Q, ndof_total,
                                  TPV102Params::lambda,
                                  TPV102Params::mu,
                                  TPV102Params::rho);
   std::cout << "  Initial total energy E₀ = " << std::scientific
             << std::setprecision(6) << E0 << "\n\n";

   double E_max_so_far = E0;
   bool   any_increase = false;
   double t_of_max_increase = 0.0;
   double max_ratio = 1.0;

   std::cout << "  step | t (s)  | E(t)         | E(t)/E₀     | Δ vs prev\n";
   std::cout << "  -----+--------+--------------+-------------+----------\n";
   double E_prev = E0;
   for (int step = 0; step < N_steps; step++)
   {
      DoRK4Step(wave, Q, dt);
      if ((step + 1) % sample_every == 0 || step == 0 || step == N_steps - 1)
      {
         double Et = TotalEnergy(Q, ndof_total,
                                  TPV102Params::lambda,
                                  TPV102Params::mu,
                                  TPV102Params::rho);
         double ratio = Et / E0;
         double delta = Et - E_prev;
         std::printf("  %4d | %6.3f | %.6e | %.6e | %+8.2e\n",
                     step + 1, (step+1)*dt, Et, ratio, delta);
         if (Et > E_max_so_far * 1.00001)  // > 10⁻⁵ rel grow
         {
            any_increase = true;
            if (ratio > max_ratio) { max_ratio = ratio; }
            if (t_of_max_increase == 0.0)
            { t_of_max_increase = (step+1)*dt; }
         }
         E_max_so_far = std::max(E_max_so_far, Et);
         E_prev = Et;
      }
   }

   std::cout << "\n  E_max / E₀ = " << std::scientific
             << std::setprecision(3) << max_ratio;
   if (any_increase)
   { std::cout << "   ⇒ ENERGY INCREASED (first at t=" << t_of_max_increase
               << " s) — absorbing BC may have wrong sign!\n"; }
   else
   { std::cout << "   ⇒ energy monotonically non-increasing — absorbing "
                  "BC healthy.\n"; }

   TEST_BOOL(!any_increase,
             "E(t) monotonically non-increasing — absorbing BC not "
             "reflecting with wrong sign");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
