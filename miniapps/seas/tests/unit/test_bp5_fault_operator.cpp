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

// Unit tests for Phase 3: BP5 (SlipComponents=2) path of RateStateFaultOperator
// Covers: construction, PreInit, Init, ComputeRHS, state roundtrip,
// stress equilibrium, size queries, SetSlipRate, below-fault handling.

#include "mfem.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"
#include "../../config/bp5_params.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Helper: Create 3D hex mesh (same as test_elasticity_operator.cpp)
// =============================================================================
static Mesh Create3DMesh(int nx, int ny, int nz,
                          real_t Lx, real_t Ly, real_t Lz)
{
   // Match Tandem coordinates: [-Lx, Lx] x [-Ly, Ly] x [-Lz, 0]
   Mesh mesh = Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                                      Element::HEXAHEDRON,
                                      2.0 * Lx, 2.0 * Ly, Lz);

   Vector shift(3);
   shift(0) = -Lx;
   shift(1) = -Ly;
    shift(2) = -Lz;

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
      // attr 3 is the BP5 fault attribute (set by
      // ElasticityDomainOperator's BoundaryConfig).  The +Ly outer
      // face is a 1-sided boundary, NOT a fault — tagging it with 3
      // would collide with the strict fault-validator (REVIEW R-001).
      // Use attr 7 (unused by the BC machinery) for this fixture's
      // +Ly face instead.
      else if (std::abs(center(1) - Ly) < tol)     { mesh.SetBdrAttribute(be, 7); }
      else if (std::abs(center(1) - (-Ly)) < tol)  { mesh.SetBdrAttribute(be, 4); }
      else if (std::abs(center(2) - 0.0) < tol)    { mesh.SetBdrAttribute(be, 5); }
      else if (std::abs(center(2) + Lz) < tol)     { mesh.SetBdrAttribute(be, 6); }
   }

   mesh.SetAttributes();
   return mesh;
}

// =============================================================================
// Shared test fixture: builds mesh + operator + geometry + fault operator
// =============================================================================
struct BP5Fixture
{
   BP5Params params;
   std::unique_ptr<Mesh> mesh;
   std::unique_ptr<ElasticityDomainOperator<Mesh>> domain_op;
   std::unique_ptr<FaultGeometry<Mesh>> fault_geom;
   std::unique_ptr<DieterichRuinaFriction> friction;
   std::unique_ptr<AgingLawPsi> evolution;
   std::unique_ptr<RateStateFaultOperator<Mesh, 2>> fault_op;

   int nf = 0;

   bool Setup()
   {
      // Domain sized to contain fault zone
      real_t Lx = 50e3, Ly = 60e3, Lz = 40e3;
      mesh = std::make_unique<Mesh>(Create3DMesh(1, 1, 1, Lx, Ly, Lz));

      domain_op = std::make_unique<ElasticityDomainOperator<Mesh>>(
         *mesh, 1, params.lambda(), params.mu(),
         params.Vp, params.Wf, params.lf, DGMethod::IP);

      nf = domain_op->GetNumFaultDOFs();
      if (nf == 0) { return false; }

      fault_geom = std::make_unique<FaultGeometry<Mesh>>(*domain_op, params);

      // Friction constants (using BP5 defaults)
      DieterichRuinaFriction::Constants fc;
      fc.V0 = params.V0;
      fc.f0 = params.f0;
      fc.b = params.b;
      fc.Dc = params.L0;  // Default; per-DOF Dc handled by fault operator
      friction = std::make_unique<DieterichRuinaFriction>(fc);

      // Aging law in psi-space (BP5 always uses psi)
      evolution = std::make_unique<AgingLawPsi>(params.b, params.V0, params.f0);

      fault_op = std::make_unique<RateStateFaultOperator<Mesh, 2>>(
         fault_geom.get(), friction.get(), evolution.get(), params);

      return true;
   }
};

static real_t TandemPsiInit(real_t tau_abs, real_t V_abs_init, real_t sigma_n,
                            real_t eta, real_t a, real_t V0)
{
   const real_t arg = (tau_abs - eta * V_abs_init) / (a * sigma_n);
   const real_t s = std::sinh(arg);
   return a * std::log((2.0 * V0 / V_abs_init) * s);
}

static void TandemSlipRateVectorPsi(const real_t tau_vec[2], real_t psi,
                                    real_t sigma_n, real_t eta, real_t a,
                                    real_t V0, DieterichRuinaFriction &friction,
                                    real_t V_vec[2])
{
   const real_t tau_abs = std::sqrt(tau_vec[0] * tau_vec[0] +
                                    tau_vec[1] * tau_vec[1]);
   if (tau_abs <= 0.0)
   {
      V_vec[0] = 0.0;
      V_vec[1] = 0.0;
      return;
   }

   const real_t V_abs = friction.SolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a);
   // v55 D8: V_vec anti-parallel to tau (Tandem convention)
   V_vec[0] = -(V_abs / tau_abs) * tau_vec[0];
   V_vec[1] = -(V_abs / tau_abs) * tau_vec[1];
   (void)V0;
}

// =============================================================================
// Test 1: BP5 Construction — size queries
// =============================================================================
void TestBP5Construction()
{
   std::cout << "\n--- Test: BP5 RateStateFaultOperator Construction ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   TEST_ASSERT(fix.fault_op->NumNodes() == N, "NumNodes matches fault DOF count");
   TEST_ASSERT(fix.fault_op->StateSize() == 3 * N,
               "StateSize = 3*N (2 slip + 1 psi per node)");
   TEST_ASSERT(fix.fault_op->SlipSize() == 2 * N,
               "SlipSize = 2*N (2 components per node)");
   TEST_ASSERT(fix.fault_op->TractionSize() == 2 * N,
               "TractionSize = 2*N (2 components per node)");
   TEST_ASSERT(fix.fault_op->UsePsi(), "BP5 uses psi-space integration");

   std::cout << "  N = " << N << ", StateSize = " << fix.fault_op->StateSize()
             << ", SlipSize = " << fix.fault_op->SlipSize() << "\n";
}

// =============================================================================
// Test 2: BP5 PreInit — slip=0, psi=steady state from per-DOF V/Dc
// =============================================================================
void TestBP5PreInit()
{
   std::cout << "\n--- Test: BP5 PreInit ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.fault_op->PreInit(state);

   // All slips should be zero
   Vector slip;
   fix.fault_op->GetSlip(state, slip);
   TEST_NEAR(slip.Norml2(), 0.0, 1e-15, "PreInit: all slips are zero");

   // Psi values should be positive and reasonable
   // psi_ss = f0 + b * ln(V0 / V_init_magnitude)
   for (int i = 0; i < N; i++)
   {
      real_t psi = state(i * 3 + 2);  // PsiIndex = 2
      TEST_ASSERT(psi > 0.0, "PreInit: psi > 0");
      if (psi <= 0.0) { break; }
   }

   // Theta from psi should be positive
   Vector theta;
   fix.fault_op->GetTheta(state, theta);
   TEST_ASSERT(theta.Min() > 0.0, "PreInit: all theta values > 0");

   std::cout << "  Theta range: [" << theta.Min() << ", " << theta.Max() << "]\n";
}

// =============================================================================
// Test 3: BP5 Init — stress equilibrium and V_max ~ V_nuc
// =============================================================================
void TestBP5Init()
{
   std::cout << "\n--- Test: BP5 Init ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.fault_op->PreInit(state);

   // Zero traction (before domain solve)
   Vector traction(fix.fault_op->TractionSize());
   traction = 0.0;

   real_t V_max = fix.fault_op->Init(traction, state);

   // V_max should be positive and reasonable
   TEST_ASSERT(V_max > 0.0, "Init: V_max > 0");

   // V_max should be close to V_nuc if nucleation zone DOFs exist,
   // or close to V_init otherwise
   real_t V_ref = fix.fault_op->GetReferenceVInit();
   real_t V_rel_err = std::abs(V_max - V_ref) / std::max(V_ref, 1e-30);
   TEST_ASSERT(V_rel_err < 0.5,
               "Init: V_max within 50% of reference V_init");

   // Psi values should be updated by Init
   for (int i = 0; i < N; i++)
   {
      real_t psi = state(i * 3 + 2);
      TEST_ASSERT(std::isfinite(psi), "Init: psi is finite");
      if (!std::isfinite(psi)) { break; }
   }

   std::cout << "  V_max = " << V_max << " m/s, V_ref = " << V_ref << "\n";
}

// =============================================================================
// Test 4: BP5 ComputeRHS — V direction and dpsi/dt sign
// =============================================================================
void TestBP5ComputeRHS()
{
   std::cout << "\n--- Test: BP5 ComputeRHS ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.fault_op->PreInit(state);

   Vector traction(fix.fault_op->TractionSize());
   traction = 0.0;

   fix.fault_op->Init(traction, state);

   // Compute RHS
   Vector rate(fix.fault_op->StateSize());
   real_t V_max = fix.fault_op->ComputeRHS(traction, state, rate);

   TEST_ASSERT(V_max > 0.0, "ComputeRHS: V_max > 0");

   // Check rate vector has correct size
   TEST_ASSERT(rate.Size() == 3 * N, "ComputeRHS: rate size = 3*N");

   // For each node, check velocity magnitude is positive
   bool all_V_positive = true;
   bool all_dpsi_finite = true;
   for (int i = 0; i < N; i++)
   {
      real_t V0 = rate(i * 3 + 0);  // dip rate
      real_t V1 = rate(i * 3 + 1);  // strike rate
      real_t V_abs = std::sqrt(V0 * V0 + V1 * V1);
      real_t dpsi_dt = rate(i * 3 + 2);

      if (V_abs <= 0.0) { all_V_positive = false; }
      if (!std::isfinite(dpsi_dt)) { all_dpsi_finite = false; }
   }

   TEST_ASSERT(all_V_positive, "ComputeRHS: all |V| > 0");
   TEST_ASSERT(all_dpsi_finite, "ComputeRHS: all dpsi/dt are finite");

   // v55 D8: Slip rate direction is anti-parallel to tau (Tandem convention)
   // V_vec = -(V_abs / tau_abs) * tau_vec, so V · tau_pre < 0
   const Vector &tau_pre = fix.fault_geom->GetTauPre();
   bool direction_ok = true;
   for (int i = 0; i < N; i++)
   {
      real_t V0 = rate(i * 3 + 0);
      real_t V1 = rate(i * 3 + 1);
      real_t tp0 = tau_pre(2 * i);
      real_t tp1 = tau_pre(2 * i + 1);

      // tau_total = tau_pre + traction (traction = 0 here)
      real_t dot = V0 * tp0 + V1 * tp1;
      // V should be anti-parallel to tau_pre → dot < 0 (Tandem convention)
      if (dot > 1e-20)
      {
         direction_ok = false;
         break;
      }
   }
   TEST_ASSERT(direction_ok,
               "ComputeRHS: V direction is anti-parallel to tau_pre (Tandem)");

   std::cout << "  V_max from RHS: " << V_max << " m/s\n";
}

// =============================================================================
// Test 4b: BP5 node update matches explicit Tandem-style formulas
// =============================================================================
void TestBP5NodeUpdateMatchesTandemSource()
{
   std::cout << "\n--- Test: BP5 Node Update vs Tandem Source Formulas ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   const int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.fault_op->PreInit(state);

   Vector traction(fix.fault_op->TractionSize());
   traction = 0.0;

   fix.fault_op->Init(traction, state);

   Vector rate(fix.fault_op->StateSize());
   fix.fault_op->ComputeRHS(traction, state, rate);

   const Vector &a_values = fix.fault_geom->GetAValues();
   const Vector &eta_values = fix.fault_geom->GetEtaValues();
   const Vector &dc_values = fix.fault_geom->GetDcValues();
   const Vector &tau_pre = fix.fault_geom->GetTauPre();
   const Vector &V_init = fix.fault_geom->GetVInit();

   int node = -1;
   for (int i = 0; i < N; i++)
   {
      real_t V0i = V_init(2 * i);
      real_t V1i = V_init(2 * i + 1);
      real_t Vabs = std::sqrt(V0i * V0i + V1i * V1i);
      if (Vabs > fix.params.V_init * 10.0)
      {
         node = i;
         break;
      }
   }
   if (node < 0) { node = 0; }

   const real_t tau_vec[2] = {tau_pre(2 * node) + traction(2 * node),
                              tau_pre(2 * node + 1) + traction(2 * node + 1)};
   const real_t tau_abs = std::sqrt(tau_vec[0] * tau_vec[0] +
                                    tau_vec[1] * tau_vec[1]);
   const real_t V_abs_init = std::sqrt(V_init(2 * node) * V_init(2 * node) +
                                       V_init(2 * node + 1) * V_init(2 * node + 1));
   const real_t a = a_values(node);
   const real_t eta = eta_values(node);
   const real_t Dc = dc_values(node);
   const real_t sigma_n = fix.params.sigma_n;

   const real_t psi_expected = TandemPsiInit(
      tau_abs, V_abs_init, sigma_n, eta, a, fix.params.V0);
   const real_t psi_actual = state(node * 3 + 2);

   real_t V_expected[2];
   TandemSlipRateVectorPsi(tau_vec, psi_actual, sigma_n, eta, a,
                           fix.params.V0, *fix.friction, V_expected);
   const real_t dpsi_expected =
      fix.params.b * fix.params.V0 / Dc
      * (std::exp((fix.params.f0 - psi_actual) / fix.params.b)
         - std::sqrt(V_expected[0] * V_expected[0]
                     + V_expected[1] * V_expected[1]) / fix.params.V0);

   const real_t V_actual[2] = {rate(node * 3 + 0), rate(node * 3 + 1)};
   const real_t dpsi_actual = rate(node * 3 + 2);

   std::cout << "  node=" << node
             << " psi_expected=" << psi_expected
             << " psi_actual=" << psi_actual << "\n";
   std::cout << "  V_expected=(" << V_expected[0] << "," << V_expected[1] << ")"
             << " V_actual=(" << V_actual[0] << "," << V_actual[1] << ")\n";
   std::cout << "  dpsi_expected=" << dpsi_expected
             << " dpsi_actual=" << dpsi_actual << "\n";

   TEST_NEAR(psi_actual, psi_expected, 1e-12,
             "Init psi matches Tandem-style psi_init");
   TEST_NEAR(V_actual[0], V_expected[0], 1e-15,
             "ComputeRHS dip rate matches Tandem-style slip_rate");
   TEST_NEAR(V_actual[1], V_expected[1], 1e-15,
             "ComputeRHS strike rate matches Tandem-style slip_rate");
   TEST_NEAR(dpsi_actual, dpsi_expected, 1e-12,
             "ComputeRHS dpsi/dt matches Tandem-style aging law");
}

void TestBP5ComputeRHSMatchesTandemSourceAllNodes()
{
   std::cout << "\n--- Test: BP5 ComputeRHS vs Tandem Source (all nodes, nonzero traction) ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   const int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   fix.fault_op->PreInit(state);

   Vector traction0(fix.fault_op->TractionSize());
   traction0 = 0.0;
   fix.fault_op->Init(traction0, state);

   Vector slip_bc(fix.fault_op->SlipSize());
   slip_bc = 0.0;
   for (int i = 0; i < N; i++)
   {
      slip_bc(2 * i) = 0.02 * ((i % 3) - 1);
      slip_bc(2 * i + 1) = -0.5 - 0.05 * i;
   }
   fix.fault_op->SetSlip(slip_bc, state);

   GridFunction u(&fix.domain_op->GetFESpace());
   u = 0.0;
   fix.domain_op->Solve(0.0, slip_bc, u);

   Vector traction(fix.fault_op->TractionSize());
   Vector normal_traction(N);
   fix.domain_op->ComputeTraction(u, slip_bc, traction, &normal_traction);

   Vector rate(fix.fault_op->StateSize());
   fix.fault_op->ComputeRHS(traction, state, rate, &normal_traction);

   const Vector &a_values = fix.fault_geom->GetAValues();
   const Vector &eta_values = fix.fault_geom->GetEtaValues();
   const Vector &dc_values = fix.fault_geom->GetDcValues();
   const Vector &depths = fix.fault_geom->GetDepths();
   const Vector &tau_pre = fix.fault_geom->GetTauPre();

   real_t max_V_diff = 0.0;
   real_t max_dpsi_diff = 0.0;
   int max_V_node = -1;
   int max_dpsi_node = -1;

   for (int i = 0; i < N; i++)
   {
      real_t V_expected[2];
      real_t dpsi_expected = 0.0;

      // All DOFs (including below-Wf) handled by friction solver naturally.
      {
         const real_t psi = state(i * 3 + 2);
         const real_t tau_vec[2] = {
            tau_pre(2 * i) + traction(2 * i),
            tau_pre(2 * i + 1) + traction(2 * i + 1)
         };
         const real_t a = a_values(i);
         const real_t eta = eta_values(i);
         const real_t Dc = dc_values(i);

         // v55 P1: sigma_n_eff = SnPre + normal_traction (Tandem convention)
         real_t sigma_n_eff = fix.params.sigma_n + normal_traction(i);
         sigma_n_eff = std::max(sigma_n_eff, 0.1 * fix.params.sigma_n);

         TandemSlipRateVectorPsi(tau_vec, psi, sigma_n_eff, eta, a,
                                 fix.params.V0, *fix.friction, V_expected);

         const real_t V_abs = std::sqrt(V_expected[0] * V_expected[0] +
                                        V_expected[1] * V_expected[1]);
         dpsi_expected = fix.params.b * fix.params.V0 / Dc
            * (std::exp((fix.params.f0 - psi) / fix.params.b)
               - V_abs / fix.params.V0);
      }

      const real_t V0_actual = rate(i * 3 + 0);
      const real_t V1_actual = rate(i * 3 + 1);
      const real_t dpsi_actual = rate(i * 3 + 2);

      max_V_diff = std::max(max_V_diff, std::abs(V0_actual - V_expected[0]));
      if (std::abs(V0_actual - V_expected[0]) == max_V_diff) { max_V_node = i; }
      max_V_diff = std::max(max_V_diff, std::abs(V1_actual - V_expected[1]));
      if (std::abs(V1_actual - V_expected[1]) == max_V_diff) { max_V_node = i; }
      if (std::abs(dpsi_actual - dpsi_expected) > max_dpsi_diff)
      {
         max_dpsi_diff = std::abs(dpsi_actual - dpsi_expected);
         max_dpsi_node = i;
      }
   }

   std::cout << "  max |dV|    = " << max_V_diff
             << " at node " << max_V_node << "\n";
   std::cout << "  max |dpsi|  = " << max_dpsi_diff
             << " at node " << max_dpsi_node << "\n";

   TEST_ASSERT(std::isfinite(max_V_diff), "All-node V mismatch is finite");
   TEST_ASSERT(std::isfinite(max_dpsi_diff), "All-node dpsi mismatch is finite");
   TEST_ASSERT(max_V_diff < 1e-13,
               "ComputeRHS slip-rate vector matches Tandem-style formulas at all nodes");
   TEST_ASSERT(max_dpsi_diff < 1e-12,
               "ComputeRHS dpsi matches Tandem-style formulas at all nodes");
}

void TestBP5RejectedStepPurity()
{
   std::cout << "\n--- Test: BP5/IP Rejected-Step Purity ---\n";

   BP5Fixture fix_reject;
   BP5Fixture fix_fresh;
   if (!fix_reject.Setup() || !fix_fresh.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   BP5SEASOp seas_reject(fix_reject.domain_op.get(), fix_reject.fault_op.get());
   BP5SEASOp seas_fresh(fix_fresh.domain_op.get(), fix_fresh.fault_op.get());
   seas_reject.SetElasticSigmaN(true);
   seas_fresh.SetElasticSigmaN(true);

   Vector state0(fix_reject.fault_op->StateSize());
   seas_reject.SetInitialCondition(state0);

   // Use a heterogeneous slip perturbation so the step is nontrivial; the
   // exact steady initial state can accept even extremely large dt on this
   // coarse fixture and does not exercise rejection behavior.
   Vector slip_seed(fix_reject.fault_op->SlipSize());
   slip_seed = 0.0;
   for (int i = 0; i < fix_reject.nf; i++)
   {
      slip_seed(2 * i) = 0.02 * ((i % 3) - 1);
      slip_seed(2 * i + 1) = -0.5 - 0.05 * i;
   }
   fix_reject.fault_op->SetSlip(slip_seed, state0);

   Vector state_reject = state0;
   Vector state_fresh = state0;

   DormandPrinceRK45 rk_reject;
   rk_reject.SetAbsTol(1e-7);
   rk_reject.SetRelTol(1e-50);
   rk_reject.SetDt(1.0e6);
   rk_reject.SetDtMax(1.0e6);
   rk_reject.SetStatePerNode(3);
   rk_reject.SetVGuard(1.05);
   rk_reject.Init(seas_reject);

   real_t t_reject = 0.0;
   real_t accepted_dt = -1.0;
   int attempts = 0;
   while (attempts < 40)
   {
      real_t dt_try = 0.0;
      if (rk_reject.Step(seas_reject, state_reject, t_reject, dt_try))
      {
         accepted_dt = dt_try;
         break;
      }
      attempts++;
   }

   TEST_ASSERT(accepted_dt > 0.0,
               "Rejected-path solver eventually accepts a BP5/IP step");
   TEST_ASSERT(rk_reject.GetTotalRejections() > 0,
               "Rejected-path solver incurred at least one rejection");

   DormandPrinceRK45 rk_fresh;
   rk_fresh.SetAbsTol(1e-7);
   rk_fresh.SetRelTol(1e-50);
   rk_fresh.SetDt(accepted_dt);
   rk_fresh.SetDtMax(1.0e6);
   rk_fresh.SetStatePerNode(3);
   rk_fresh.SetVGuard(1.05);
   rk_fresh.Init(seas_fresh);

   real_t t_fresh = 0.0;
   real_t dt_fresh = 0.0;
   bool accepted_fresh = rk_fresh.Step(seas_fresh, state_fresh, t_fresh, dt_fresh);
   TEST_ASSERT(accepted_fresh,
               "Fresh-path solver accepts the same reduced BP5/IP step");

   Vector state_diff(state_reject.Size());
   state_diff = state_reject;
   state_diff -= state_fresh;
   real_t rel_state_diff = state_diff.Norml2() /
                           std::max(state_reject.Norml2(), 1e-30);

   Vector trac_diff(seas_reject.GetTraction().Size());
   trac_diff = seas_reject.GetTraction();
   trac_diff -= seas_fresh.GetTraction();
   real_t rel_trac_diff = trac_diff.Norml2() /
                          std::max(seas_reject.GetTraction().Norml2(), 1e-30);

   real_t vmax_reject = seas_reject.GetMaxSlipRate();
   real_t vmax_fresh = seas_fresh.GetMaxSlipRate();
   real_t rel_vmax_diff = std::abs(vmax_reject - vmax_fresh) /
                          std::max(std::abs(vmax_reject), 1e-30);

   std::cout << "  rejected attempts = " << rk_reject.GetTotalRejections()
             << ", accepted_dt = " << accepted_dt << "\n";
   std::cout << "  rel_state_diff = " << rel_state_diff
             << ", rel_trac_diff = " << rel_trac_diff
             << ", rel_vmax_diff = " << rel_vmax_diff << "\n";

   TEST_ASSERT(std::abs(t_reject - t_fresh) < 1e-12,
               "Accepted time increment matches after rejection vs fresh path");
   TEST_ASSERT(rel_state_diff < 1e-11,
               "Rejected BP5/IP attempt does not change next accepted state");
   TEST_ASSERT(rel_trac_diff < 1e-11,
               "Rejected BP5/IP attempt does not change next accepted traction");
   TEST_ASSERT(rel_vmax_diff < 1e-11,
               "Rejected BP5/IP attempt does not change next accepted V_max");
}


// =============================================================================
// Test 5: BP5 State Access Roundtrip — GetSlip/SetSlip, GetTheta/SetTheta
// =============================================================================
void TestBP5StateRoundtrip()
{
   std::cout << "\n--- Test: BP5 State Access Roundtrip ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector state(fix.fault_op->StateSize());
   state = 0.0;

   // Set slip values (2 components per node)
   Vector slip_in(2 * N);
   for (int i = 0; i < N; i++)
   {
      slip_in(2 * i) = 0.001 * (i + 1);       // dip slip
      slip_in(2 * i + 1) = 0.002 * (i + 1);   // strike slip
   }
   fix.fault_op->SetSlip(slip_in, state);

   // Set theta values (1 per node)
   Vector theta_in(N);
   for (int i = 0; i < N; i++)
   {
      theta_in(i) = 1000.0 + 100.0 * i;
   }
   fix.fault_op->SetTheta(theta_in, state);

   // Get back and compare
   Vector slip_out, theta_out;
   fix.fault_op->GetSlip(state, slip_out);

   TEST_ASSERT(slip_out.Size() == 2 * N, "GetSlip: correct size (2*N)");

   real_t slip_diff = 0.0;
   for (int i = 0; i < 2 * N; i++)
   {
      slip_diff = std::max(slip_diff, std::abs(slip_in(i) - slip_out(i)));
   }
   TEST_NEAR(slip_diff, 0.0, 1e-15, "SetSlip/GetSlip roundtrip exact");

   // Theta roundtrip: SetTheta stores raw values in state, but GetTheta
   // converts psi→theta via PsiToTheta since use_psi_=true.
   // So we need to test the raw state values directly.
   for (int i = 0; i < N; i++)
   {
      real_t psi_stored = state(i * 3 + 2);
      TEST_NEAR(psi_stored, theta_in(i), 1e-15,
                "SetTheta stores raw value at PsiIndex");
      if (std::abs(psi_stored - theta_in(i)) > 1e-15) { break; }
   }
}

// =============================================================================
// Test 6: BP5 Stress Equilibrium — VerifyStressEquilibrium < 1e-8
// =============================================================================
void TestBP5StressEquilibrium()
{
   std::cout << "\n--- Test: BP5 Stress Equilibrium ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   Vector state(fix.fault_op->StateSize());
   fix.fault_op->PreInit(state);

   Vector traction(fix.fault_op->TractionSize());
   traction = 0.0;

   fix.fault_op->Init(traction, state);

   real_t max_error = fix.fault_op->VerifyStressEquilibrium(traction, state);
   TEST_ASSERT(max_error < 1e-8,
               "Stress equilibrium satisfied after Init (error < 1e-8)");
   std::cout << "  Max stress balance error: " << max_error << "\n";
}

// =============================================================================
// Test 7: BP5 VerifyInitialSlipRate — checks the §2.1 fix
// =============================================================================
void TestBP5VerifyInitialSlipRate()
{
   std::cout << "\n--- Test: BP5 VerifyInitialSlipRate ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   Vector state(fix.fault_op->StateSize());
   fix.fault_op->PreInit(state);

   Vector traction(fix.fault_op->TractionSize());
   traction = 0.0;

   real_t V_max = fix.fault_op->Init(traction, state);

   // GetReferenceVInit should return max |V_init| (accounts for nucleation)
   real_t V_ref = fix.fault_op->GetReferenceVInit();

   // V_ref should be >= bp5_params.V_init (may be V_nuc in nucleation zone)
   TEST_ASSERT(V_ref >= fix.params.V_init,
               "GetReferenceVInit >= V_init (accounts for nucleation)");

   // VerifyInitialSlipRate should NOT crash (this was the §2.1 bug)
   // It should pass since V_max ~ V_ref
   fix.fault_op->VerifyInitialSlipRate(V_max);
   TEST_ASSERT(true, "VerifyInitialSlipRate does not crash for BP5");

   std::cout << "  V_max = " << V_max << ", V_ref = " << V_ref << "\n";
}

// =============================================================================
// Test 8: BP5 SetSlipRate — V_max from vector magnitudes
// =============================================================================
void TestBP5SetSlipRate()
{
   std::cout << "\n--- Test: BP5 SetSlipRate ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   int N = fix.nf;
   Vector V(2 * N);
   for (int i = 0; i < N; i++)
   {
      V(2 * i) = 0.001 * (i + 1);       // dip component
      V(2 * i + 1) = 0.002 * (i + 1);   // strike component
   }

   fix.fault_op->SetSlipRate(V);

   // V_max should be the maximum vector magnitude
   real_t expected_V_max = 0.0;
   for (int i = 0; i < N; i++)
   {
      real_t V_abs = std::sqrt(V(2*i)*V(2*i) + V(2*i+1)*V(2*i+1));
      expected_V_max = std::max(expected_V_max, V_abs);
   }

   TEST_NEAR(fix.fault_op->GetMaxSlipRate(), expected_V_max, 1e-15,
             "SetSlipRate: V_max computed from vector magnitudes");

   // Verify stored slip rate matches
   const Vector &V_stored = fix.fault_op->GetSlipRate();
   TEST_ASSERT(V_stored.Size() == 2 * N, "GetSlipRate: correct size");

   real_t V_diff = 0.0;
   for (int i = 0; i < 2 * N; i++)
   {
      V_diff = std::max(V_diff, std::abs(V(i) - V_stored(i)));
   }
   TEST_NEAR(V_diff, 0.0, 1e-15, "SetSlipRate/GetSlipRate roundtrip exact");
}

// =============================================================================
// Test 9: BP5 PrintState — exercises §4.2 fix
// =============================================================================
void TestBP5PrintState()
{
   std::cout << "\n--- Test: BP5 PrintState ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   Vector state(fix.fault_op->StateSize());
   fix.fault_op->PreInit(state);

   Vector traction(fix.fault_op->TractionSize());
   traction = 0.0;
   fix.fault_op->Init(traction, state);

   std::ostringstream oss;
   fix.fault_op->PrintState(state, oss);
   std::string output = oss.str();

   // Should contain |Slip| and |V| ranges (from §4.2 fix)
   TEST_ASSERT(output.find("|Slip| range") != std::string::npos,
               "PrintState outputs |Slip| range for BP5");
   TEST_ASSERT(output.find("|V| range") != std::string::npos,
               "PrintState outputs |V| range for BP5");
   TEST_ASSERT(output.find("V_max") != std::string::npos,
               "PrintState outputs V_max");

   std::cout << "  PrintState output:\n" << output;
}

// =============================================================================
// Test 10: BP5 GetBP5Params accessor
// =============================================================================
void TestBP5ParamsAccessor()
{
   std::cout << "\n--- Test: BP5 GetBP5Params Accessor ---\n";

   BP5Fixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   const BP5Params &p = fix.fault_op->GetBP5Params();
   TEST_NEAR(p.sigma_n, fix.params.sigma_n, 1e-10,
             "GetBP5Params: sigma_n matches");
   TEST_NEAR(p.b, fix.params.b, 1e-15,
             "GetBP5Params: b matches");
   TEST_NEAR(p.V0, fix.params.V0, 1e-20,
             "GetBP5Params: V0 matches");

   real_t sigma_n = fix.fault_op->GetSigmaN();
   TEST_NEAR(sigma_n, fix.params.sigma_n, 1e-10,
             "GetSigmaN matches BP5 sigma_n");
}

// =============================================================================
// Test 11: BP5 Below-Fault Handling — friction solver handles all DOFs
// naturally (no hardcoded below-Wf branch), matching Tandem's approach.
// Below-fault DOFs are velocity-strengthening and should converge to
// steady-state plate rate via the friction solver.
// =============================================================================
void TestBP5BelowFault()
{
   std::cout << "\n--- Test: BP5 Below-Fault Handling (natural, no hardcoded branch) ---\n";

   BP5Params params;

   // Strategy: construct domain operator with large Wf to detect all fault
   // faces, but pass a smaller Wf to the fault operator via BP5Params so
   // that some DOFs are "below fault zone" (depth > Wf_small + 1.0).
   //
   // Mesh: 50km x 60km x 10km, with 2 elements in z → centroids at ~2.5km, ~7.5km.
   // Domain operator Wf = 10km → detects all faces.
   // BP5Params Wf = 5km → DOFs at ~7.5km are "below fault" (7.5 > 5 + 1 = 6).
   real_t Lx = 50e3, Ly = 60e3, Lz = 10e3;
   real_t domain_Wf = Lz;  // Detect all faces
   Mesh mesh = Create3DMesh(1, 1, 2, Lx, Ly, Lz);

   ElasticityDomainOperator<Mesh> domain_op(
      mesh, 1, params.lambda(), params.mu(),
      params.Vp, domain_Wf, params.lf, DGMethod::IP);

   int nf = domain_op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Build FaultGeometry with the normal params (it uses coords from domain_op)
   FaultGeometry<Mesh> fault_geom(domain_op, params);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.L0;
   DieterichRuinaFriction friction(fc);
   AgingLawPsi evolution(params.b, params.V0, params.f0);

   // Use smaller Wf in BP5Params for the fault operator
   BP5Params params_small_wf = params;
   params_small_wf.Wf = 5.0e3;  // 5 km — smaller than some DOF depths

   RateStateFaultOperator<Mesh, 2> fault_op(
      &fault_geom, &friction, &evolution, params_small_wf);

   int N = fault_op.NumNodes();
   Vector state(fault_op.StateSize());
   fault_op.PreInit(state);

   Vector traction(fault_op.TractionSize());
   traction = 0.0;
   fault_op.Init(traction, state);

   Vector rate(fault_op.StateSize());
   fault_op.ComputeRHS(traction, state, rate);

   // All DOFs (including below-Wf) are handled by the friction solver.
   // At initial equilibrium with zero elastic traction, all DOFs should
   // have V ≈ V_init. For below-fault DOFs: V_init = (0, Vp).
   // The friction solver: V_vec = -(tau_vec/|tau|)*V_abs
   //   tau_pre = (0, -tau0_scalar) → V_vec = (0, +Vp)
   const Vector &depths = fault_geom.GetDepths();
   real_t Wf_small = params_small_wf.Wf;
   int below_count = 0;
   int above_count = 0;
   for (int i = 0; i < N; i++)
   {
      if (depths(i) > Wf_small + 1.0)
      {
         below_count++;
         // Friction solver should give V ≈ (0, +Vp) at steady state
         real_t V_dip = rate(i * 3 + 0);
         real_t V_strike = rate(i * 3 + 1);
         real_t V_abs = std::sqrt(V_dip*V_dip + V_strike*V_strike);
         // V_abs should be close to Vp (within friction solver tolerance)
         TEST_NEAR(V_abs, params_small_wf.Vp, params_small_wf.Vp * 1e-6,
                   "Below-fault: |V| ≈ Vp (friction solver)");
         // Strike component should be positive (right-lateral)
         TEST_ASSERT(V_strike > 0.0,
                     "Below-fault: strike rate > 0 (right-lateral)");
         // Dip component should be near zero
         TEST_NEAR(V_dip, 0.0, params_small_wf.Vp * 1e-6,
                   "Below-fault: dip rate ≈ 0");
      }
      else
      {
         above_count++;
      }
   }

   TEST_ASSERT(below_count > 0,
               "At least one DOF is below fault zone (depth > Wf)");
   TEST_ASSERT(above_count > 0,
               "At least one DOF is within fault zone (depth <= Wf)");
   std::cout << "  " << below_count << " below, " << above_count
             << " within fault zone (Wf=" << Wf_small/1e3 << "km)\n";
}

// =============================================================================
// Main
// =============================================================================
int main()
{
   std::cout << "========================================\n";
   std::cout << "  BP5 Fault Operator Tests (Phase 3)\n";
   std::cout << "========================================\n";

   TestBP5Construction();
   TestBP5PreInit();
   TestBP5Init();
   TestBP5ComputeRHS();
   TestBP5NodeUpdateMatchesTandemSource();
   TestBP5ComputeRHSMatchesTandemSourceAllNodes();
   TestBP5RejectedStepPurity();
   TestBP5StateRoundtrip();
   TestBP5StressEquilibrium();
   TestBP5VerifyInitialSlipRate();
   TestBP5SetSlipRate();
   TestBP5PrintState();
   TestBP5ParamsAccessor();
   TestBP5BelowFault();

   TEST_PRINT_RESULTS();

   // REVIEW R-002 (round-4): a "0 tests ran" outcome (BP5Fixture::Setup
   // returning false because GetNumFaultDOFs()==0) silently passes
   // exit=0 — misleading CI signal.  Treat it as a failure so the
   // user sees the fixture short-circuit.
   if (num_tests == 0)
   {
      std::cerr << "*** FAIL: 0 tests ran (BP5Fixture::Setup likely "
                << "returned false because GetNumFaultDOFs()==0).  "
                << "Fix the fixture to produce fault DOFs or remove "
                << "the test target.\n";
      return 1;
   }
   return num_failed;
}
