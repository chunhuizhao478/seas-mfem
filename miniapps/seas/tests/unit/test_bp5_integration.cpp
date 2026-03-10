// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Phase 4 BP5 integration tests: SEASQuasiDynamicOperator wiring.
// Verifies full BP5 coupling stack end-to-end:
//   SetInitialCondition, Mult(), RK4, DormandPrinceRK45, stress balance.

#include "mfem.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp5_params.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <cmath>
#include <memory>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Helper: Create 3D hex mesh (same as test_bp5_fault_operator.cpp)
// =============================================================================
static Mesh Create3DMesh(int nx, int ny, int nz,
                          real_t Lx, real_t Ly, real_t Lz)
{
   Mesh mesh = Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                                      Element::HEXAHEDRON,
                                      2.0 * Lx, 2.0 * Ly, Lz);

   Vector shift(3);
   shift(0) = -Lx;
   shift(1) = -Ly;
   shift(2) = 0.0;

   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      for (int d = 0; d < 3; d++)
      {
         v[d] += shift(d);
      }
   }

   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      ElementTransformation *T = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(T->GetGeometryType());
      T->SetIntPoint(&ip);
      Vector center(3);
      T->Transform(ip, center);

      real_t tol = 1e-6;
      if (std::abs(center(0) - (-Lx)) < tol)      { mesh.SetBdrAttribute(be, 1); }
      else if (std::abs(center(0) - Lx) < tol)     { mesh.SetBdrAttribute(be, 2); }
      else if (std::abs(center(1) - Ly) < tol)     { mesh.SetBdrAttribute(be, 3); }
      else if (std::abs(center(1) - (-Ly)) < tol)  { mesh.SetBdrAttribute(be, 4); }
      else if (std::abs(center(2) - 0.0) < tol)    { mesh.SetBdrAttribute(be, 5); }
      else if (std::abs(center(2) - Lz) < tol)     { mesh.SetBdrAttribute(be, 6); }
   }

   mesh.SetAttributes();
   return mesh;
}

// =============================================================================
// Shared test fixture: mesh + domain + fault + SEAS operator
// =============================================================================
struct BP5IntegrationFixture
{
   BP5Params params;
   std::unique_ptr<Mesh> mesh;
   std::unique_ptr<ElasticityDomainOperator<Mesh>> domain_op;
   std::unique_ptr<FaultGeometry<Mesh>> fault_geom;
   std::unique_ptr<DieterichRuinaFriction> friction;
   std::unique_ptr<AgingLawPsi> evolution;
   std::unique_ptr<RateStateFaultOperator<Mesh, 2>> fault_op;
   std::unique_ptr<BP5SEASOp> seas_op;
   int nf = 0;

   bool Setup()
   {
      real_t Lx = 50e3, Ly = 60e3, Lz = 40e3;
      mesh = std::make_unique<Mesh>(Create3DMesh(1, 1, 1, Lx, Ly, Lz));

      domain_op = std::make_unique<ElasticityDomainOperator<Mesh>>(
         *mesh, 1, params.lambda(), params.mu(),
         params.Vp, params.Wf, params.lf, DGMethod::IP);

      nf = domain_op->GetNumFaultDOFs();
      if (nf == 0) { return false; }

      fault_geom = std::make_unique<FaultGeometry<Mesh>>(*domain_op, params);

      DieterichRuinaFriction::Constants fc;
      fc.V0 = params.V0;
      fc.f0 = params.f0;
      fc.b = params.b;
      fc.Dc = params.L0;
      friction = std::make_unique<DieterichRuinaFriction>(fc);

      evolution = std::make_unique<AgingLawPsi>(params.b, params.V0, params.f0);

      fault_op = std::make_unique<RateStateFaultOperator<Mesh, 2>>(
         fault_geom.get(), friction.get(), evolution.get(), params);

      seas_op = std::make_unique<BP5SEASOp>(domain_op.get(), fault_op.get());

      return true;
   }
};

// =============================================================================
// Test 1: Full Stack Construction
// =============================================================================
void TestBP5FullStackConstruction()
{
   std::cout << "\n--- Test: BP5 Full Stack Construction ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;

   // Template instantiation succeeded — verify sizes
   TEST_ASSERT(fix.fault_op->StateSize() == 3 * N,
               "StateSize = 3*N (2 slip + 1 psi per node)");
   TEST_ASSERT(fix.fault_op->SlipSize() == 2 * N,
               "SlipSize = 2*N");
   TEST_ASSERT(fix.seas_op->Height() == 3 * N,
               "SEAS operator height = StateSize = 3*N");

   // Pointer consistency
   TEST_ASSERT(fix.seas_op->GetDomain() == fix.domain_op.get(),
               "GetDomain returns correct pointer");
   TEST_ASSERT(fix.seas_op->GetFault() == fix.fault_op.get(),
               "GetFault returns correct pointer");

   std::cout << "  N = " << N << ", StateSize = " << fix.fault_op->StateSize()
             << ", SlipSize = " << fix.fault_op->SlipSize() << "\n";
}

// =============================================================================
// Test 2: SetInitialCondition
// =============================================================================
void TestBP5SetInitialCondition()
{
   std::cout << "\n--- Test: BP5 SetInitialCondition ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());

   // 4-phase init: PreInit, domain solve, Init, verify
   // Internal MFEM_VERIFY checks eq_error < 1e-6 and V within 50% of V_ref.
   // If it returns without crashing, those checks passed.
   fix.seas_op->SetInitialCondition(state);
   TEST_ASSERT(true, "SetInitialCondition completed without error");

   // V_max should be positive
   real_t V_max = fix.seas_op->GetMaxSlipRate();
   TEST_ASSERT(V_max > 0.0, "V_max > 0 after init");

   // Initial slip should be zero
   Vector slip;
   fix.fault_op->GetSlip(state, slip);
   TEST_NEAR(slip.Norml2(), 0.0, 1e-15, "Initial slip is zero");

   // Traction should be finite
   const Vector &traction = fix.seas_op->GetTraction();
   bool traction_finite = true;
   for (int i = 0; i < traction.Size(); i++)
   {
      if (!std::isfinite(traction(i)))
      {
         traction_finite = false;
         break;
      }
   }
   TEST_ASSERT(traction_finite, "Initial traction is finite");

   std::cout << "  V_max = " << V_max << " m/s, traction size = "
             << traction.Size() << "\n";
}

// =============================================================================
// Test 3: Mult Consistency
// =============================================================================
void TestBP5MultConsistency()
{
   std::cout << "\n--- Test: BP5 Mult Consistency ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.seas_op->SetInitialCondition(state);

   // Call Mult
   Vector rate(fix.fault_op->StateSize());
   fix.seas_op->Mult(state, rate);

   // Rate vector should be finite and correct size
   TEST_ASSERT(rate.Size() == 3 * N, "Rate size = 3*N");

   bool all_finite = true;
   for (int i = 0; i < rate.Size(); i++)
   {
      if (!std::isfinite(rate(i)))
      {
         all_finite = false;
         break;
      }
   }
   TEST_ASSERT(all_finite, "Rate vector is finite");

   // V > 0 for rate-state zone nodes
   bool all_V_positive = true;
   for (int i = 0; i < N; i++)
   {
      real_t V0 = rate(i * 3 + 0);
      real_t V1 = rate(i * 3 + 1);
      real_t V_abs = std::sqrt(V0 * V0 + V1 * V1);
      if (V_abs <= 0.0)
      {
         all_V_positive = false;
         break;
      }
   }
   TEST_ASSERT(all_V_positive, "|V| > 0 in rate vector");

   // Deterministic: two calls give the same result
   Vector rate2(fix.fault_op->StateSize());
   fix.seas_op->Mult(state, rate2);

   real_t diff = 0.0;
   for (int i = 0; i < rate.Size(); i++)
   {
      diff = std::max(diff, std::abs(rate(i) - rate2(i)));
   }
   TEST_NEAR(diff, 0.0, 1e-15, "Mult is deterministic");

   std::cout << "  Rate Norml2 = " << rate.Norml2() << "\n";
}

// =============================================================================
// Test 4: State Round-Trip
// =============================================================================
void TestBP5StateRoundTrip()
{
   std::cout << "\n--- Test: BP5 State Round-Trip ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.seas_op->SetInitialCondition(state);

   // GetSlip should be zero (2*N)
   Vector slip;
   fix.fault_op->GetSlip(state, slip);
   TEST_ASSERT(slip.Size() == 2 * N, "GetSlip size = 2*N");
   TEST_NEAR(slip.Norml2(), 0.0, 1e-15, "GetSlip = zero after init");

   // GetTheta should be positive (N)
   Vector theta;
   fix.fault_op->GetTheta(state, theta);
   TEST_ASSERT(theta.Size() == N, "GetTheta size = N");
   TEST_ASSERT(theta.Min() > 0.0, "GetTheta all positive");

   // GetTraction should have correct size (2*N)
   const Vector &traction = fix.seas_op->GetTraction();
   TEST_ASSERT(traction.Size() == 2 * N, "GetTraction size = 2*N");

   // GetDisplacement should be valid
   const GridFunction &u = fix.seas_op->GetDisplacement();
   TEST_ASSERT(u.Size() > 0, "GetDisplacement has non-zero size");

   bool u_finite = true;
   for (int i = 0; i < u.Size(); i++)
   {
      if (!std::isfinite(u(i)))
      {
         u_finite = false;
         break;
      }
   }
   TEST_ASSERT(u_finite, "Displacement is finite");

   std::cout << "  Theta range: [" << theta.Min() << ", " << theta.Max()
             << "], |u| = " << u.Norml2() << "\n";
}

// =============================================================================
// Test 5: Short RK4 Run (10 steps, dt=1e3)
// =============================================================================
void TestBP5ShortRK4Run()
{
   std::cout << "\n--- Test: BP5 Short RK4 Run ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.seas_op->SetInitialCondition(state);

   RK4Solver ode_solver;
   ode_solver.Init(*fix.seas_op);

   real_t t = 0.0;
   real_t dt = 1e3;
   int nsteps = 10;

   for (int step = 0; step < nsteps; step++)
   {
      ode_solver.Step(state, t, dt);
   }

   // State should be finite
   bool state_finite = true;
   for (int i = 0; i < state.Size(); i++)
   {
      if (!std::isfinite(state(i)))
      {
         state_finite = false;
         break;
      }
   }
   TEST_ASSERT(state_finite, "State is finite after 10 RK4 steps");

   // V_max should be in (0, 1) m/s (interseismic)
   // Need to call Mult to refresh internal slip rates
   Vector rate(fix.fault_op->StateSize());
   fix.seas_op->Mult(state, rate);
   real_t V_max = fix.seas_op->GetMaxSlipRate();
   TEST_ASSERT(V_max > 0.0 && V_max < 1.0,
               "V_max in (0, 1) m/s after RK4 steps");

   // Slip should have accumulated (small but nonzero)
   Vector slip;
   fix.fault_op->GetSlip(state, slip);
   TEST_ASSERT(slip.Norml2() > 0.0, "Slip accumulated after time stepping");

   // Theta should remain positive
   Vector theta;
   fix.fault_op->GetTheta(state, theta);
   TEST_ASSERT(theta.Min() > 0.0, "Theta remains positive after RK4 steps");

   std::cout << "  t = " << t << " s, V_max = " << V_max
             << " m/s, |slip| = " << slip.Norml2() << "\n";
}

// =============================================================================
// Test 6: Short RK45 Run (DormandPrinceRK45, 5 accepted steps)
// =============================================================================
void TestBP5ShortRK45Run()
{
   std::cout << "\n--- Test: BP5 Short RK45 Run ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   Vector state(fix.fault_op->StateSize());
   fix.seas_op->SetInitialCondition(state);

   DormandPrinceRK45 rk45;
   rk45.SetAbsTol(1e-7);
   rk45.SetDt(1e3);
   rk45.SetDtMax(1e6);
   rk45.SetStatePerNode(3);  // BP5: 3 components per node (slip_dip, slip_strike, psi)
   rk45.Init(*fix.seas_op);

   real_t t = 0.0;
   int accepted = 0;
   int attempts = 0;
   int max_attempts = 50;

   while (accepted < 5 && attempts < max_attempts)
   {
      real_t dt;
      if (rk45.Step(*fix.seas_op, state, t, dt))
      {
         accepted++;
      }
      attempts++;
   }

   TEST_ASSERT(accepted == 5, "RK45 accepted 5 steps within 50 attempts");

   // State should be finite
   bool state_finite = true;
   for (int i = 0; i < state.Size(); i++)
   {
      if (!std::isfinite(state(i)))
      {
         state_finite = false;
         break;
      }
   }
   TEST_ASSERT(state_finite, "State is finite after RK45 steps");

   // V_max should be reasonable
   Vector rate(fix.fault_op->StateSize());
   fix.seas_op->Mult(state, rate);
   real_t V_max = fix.seas_op->GetMaxSlipRate();
   TEST_ASSERT(V_max > 0.0 && V_max < 1.0,
               "V_max reasonable after RK45 steps");

   // dt should be positive
   TEST_ASSERT(rk45.GetDt() > 0.0, "dt > 0 after RK45 steps");

   std::cout << "  t = " << t << " s, accepted = " << accepted
             << ", attempts = " << attempts
             << ", V_max = " << V_max
             << ", dt = " << rk45.GetDt() << "\n";
}

// =============================================================================
// Test 7: Stress Balance During Time Step
// =============================================================================
void TestBP5StressBalanceDuringTimeStep()
{
   std::cout << "\n--- Test: BP5 Stress Balance During Time Step ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   Vector state(fix.fault_op->StateSize());
   fix.seas_op->SetInitialCondition(state);

   RK4Solver ode_solver;
   ode_solver.Init(*fix.seas_op);

   real_t t = 0.0;
   real_t dt = 1e3;
   real_t max_eq_error = 0.0;

   for (int step = 0; step < 3; step++)
   {
      ode_solver.Step(state, t, dt);

      // Call Mult to refresh internal traction
      Vector rate(fix.fault_op->StateSize());
      fix.seas_op->Mult(state, rate);

      // Check stress equilibrium
      real_t eq_error = fix.fault_op->VerifyStressEquilibrium(
         fix.seas_op->GetTraction(), state);
      max_eq_error = std::max(max_eq_error, eq_error);
   }

   TEST_ASSERT(max_eq_error < 1e-4,
               "Stress balance maintained during stepping (error < 1e-4)");
   std::cout << "  Max stress balance error: " << max_eq_error << "\n";
}

// =============================================================================
// Test 8: Traction at t=0 with Zero Slip
// =============================================================================
void TestBP5ZeroSlipTraction()
{
   std::cout << "\n--- Test: BP5 Traction at t=0 (Zero Slip) ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;

   // At t=0 with zero slip and zero Dirichlet loading, the domain solve
   // should produce zero displacement and zero traction.
   Vector slip(2 * N);
   slip = 0.0;

   GridFunction u_gf(&fix.domain_op->GetFESpace());
   u_gf = 0.0;

   fix.domain_op->Solve(0.0, slip, u_gf);

   // Displacement should be zero (no loading at t=0)
   TEST_NEAR(u_gf.Norml2(), 0.0, 1e-12,
             "Displacement is zero at t=0 with zero slip");

   // Compute traction from zero displacement
   Vector traction(2 * N);
   fix.domain_op->ComputeTraction(u_gf, slip, traction);

   // Traction should be zero (or very small)
   real_t max_trac = 0.0;
   for (int i = 0; i < traction.Size(); i++)
   {
      max_trac = std::max(max_trac, std::abs(traction(i)));
   }
   TEST_NEAR(max_trac, 0.0, 1e-6,
             "Traction is near-zero at t=0 with zero slip");

   std::cout << "  Max |traction| = " << max_trac << " Pa\n";
}

// =============================================================================
// Test 9: Pre-stress + Traction Yields Correct Initial Slip Rate
// =============================================================================
void TestBP5InitialSlipRateFromPreStress()
{
   std::cout << "\n--- Test: BP5 Pre-stress → Initial Slip Rate ---\n";

   BP5IntegrationFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.seas_op->SetInitialCondition(state);

   // After initialization, the traction from the domain solve (near-zero at t=0)
   // combined with tau_pre should produce initial slip rates matching V_init/V_nuc.
   const Vector &traction = fix.seas_op->GetTraction();

   // Get fault geometry for coordinate-dependent checks
   const auto &fault_geom = *fix.fault_geom;
   const Vector &coords_x2 = fault_geom.GetCoordsX2();
   const Vector &coords_x3 = fault_geom.GetCoordsX3();
   const Vector &tau_pre = fault_geom.GetTauPre();
   const Vector &V_init_vals = fault_geom.GetVInit();

   real_t max_V_err = 0.0;
   int n_checked = 0;

   for (int i = 0; i < N; i++)
   {
      // Compute total stress: tau_pre + elastic traction
      real_t tau_total[2] = {
         tau_pre(2*i) + traction(2*i),
         tau_pre(2*i+1) + traction(2*i+1)
      };
      real_t tau_abs = std::sqrt(tau_total[0]*tau_total[0] +
                                  tau_total[1]*tau_total[1]);

      real_t Vi_abs = std::sqrt(V_init_vals(2*i) * V_init_vals(2*i) +
                                 V_init_vals(2*i+1) * V_init_vals(2*i+1));

      // Verify tau_abs is in physically reasonable range (5-25 MPa)
      if (tau_abs < 1e3 || tau_abs > 1e8) { continue; }

      real_t a = fault_geom.GetAValues()(i);
      real_t eta = fault_geom.GetEtaValues()(i);

      // Use the friction solver to get V from tau_abs
      real_t psi_ss = fix.params.f0 + fix.params.b *
         std::log(fix.params.V0 / fix.params.Vp);

      real_t V_solved = fix.friction->SolveSlipRatePsi(
         tau_abs, psi_ss, fix.params.sigma_n, eta, a);

      // V_solved should match Vi_abs within tolerance
      real_t rel_err = std::abs(V_solved - Vi_abs) /
         std::max(Vi_abs, 1e-20);
      max_V_err = std::max(max_V_err, rel_err);
      n_checked++;
   }

   TEST_ASSERT(n_checked > 0, "Checked at least one fault DOF");
   TEST_ASSERT(max_V_err < 1e-4,
               "Initial slip rate matches V_init from pre-stress (rel err < 1e-4)");

   std::cout << "  Checked " << n_checked << " DOFs, max rel error = "
             << max_V_err << "\n";
}

// =============================================================================
// Main
// =============================================================================
int main()
{
   std::cout << "========================================\n";
   std::cout << "  BP5 Integration Tests (Phase 4)\n";
   std::cout << "========================================\n";

   TestBP5FullStackConstruction();
   TestBP5SetInitialCondition();
   TestBP5MultConsistency();
   TestBP5StateRoundTrip();
   TestBP5ShortRK4Run();
   TestBP5ShortRK45Run();
   TestBP5StressBalanceDuringTimeStep();
   TestBP5ZeroSlipTraction();
   TestBP5InitialSlipRateFromPreStress();

   TEST_PRINT_RESULTS();

   return num_failed;
}
