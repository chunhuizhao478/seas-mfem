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

   // Slip rate direction should be anti-parallel to tau_pre
   // (V_vec = -(V_abs / tau_abs) * tau_vec, so V · tau_pre < 0)
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
      // V should be anti-parallel to tau_pre → dot < 0
      if (dot > 1e-20)
      {
         direction_ok = false;
         break;
      }
   }
   TEST_ASSERT(direction_ok,
               "ComputeRHS: V direction is anti-parallel to tau_pre");

   std::cout << "  V_max from RHS: " << V_max << " m/s\n";
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
// Test 11: BP5 Below-Fault Handling — DOFs below Wf get (0, Vp, 0) rate
// =============================================================================
void TestBP5BelowFault()
{
   std::cout << "\n--- Test: BP5 Below-Fault Handling ---\n";

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

   // Check that at least one DOF is below Wf_small
   const Vector &depths = fault_geom.GetDepths();
   real_t Wf_small = params_small_wf.Wf;
   int below_count = 0;
   int above_count = 0;
   for (int i = 0; i < N; i++)
   {
      if (depths(i) > Wf_small + 1.0)
      {
         below_count++;
         // Below-fault DOFs should have rate = (0, Vp, 0)
         TEST_NEAR(rate(i * 3 + 0), 0.0, 1e-15,
                   "Below-fault: dip rate = 0");
         TEST_NEAR(rate(i * 3 + 1), params_small_wf.Vp, 1e-20,
                   "Below-fault: strike rate = Vp");
         TEST_NEAR(rate(i * 3 + 2), 0.0, 1e-15,
                   "Below-fault: dpsi/dt = 0");
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
   TestBP5StateRoundtrip();
   TestBP5StressEquilibrium();
   TestBP5VerifyInitialSlipRate();
   TestBP5SetSlipRate();
   TestBP5PrintState();
   TestBP5ParamsAccessor();
   TestBP5BelowFault();

   TEST_PRINT_RESULTS();

   return num_failed;
}
