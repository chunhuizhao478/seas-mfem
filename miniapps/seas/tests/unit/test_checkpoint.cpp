// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

/// @file test_checkpoint.cpp
/// @brief Unit tests for checkpoint/restart functionality
///
/// Tests:
/// 1. Checkpoint file round-trip: write -> read -> verify all fields match
/// 2. Simulation restart consistency: continuous run == checkpoint + restart

#include "mfem.hpp"
#include "../../io/checkpoint.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp2_params.hpp"

#include <iostream>
#include <cmath>
#include <cstdio>
#include <string>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Test infrastructure
// =============================================================================

static int num_tests_passed = 0;
static int num_tests_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      if (!(condition)) { \
         std::cerr << "FAIL: " << (message) << "\n"; \
         std::cerr << "  at " << __FILE__ << ":" << __LINE__ << "\n"; \
         num_tests_failed++; \
      } else { \
         std::cout << "  PASS: " << (message) << "\n"; \
         num_tests_passed++; \
      } \
   } while (0)

#define TEST_ASSERT_NEAR(a, b, tol, message) \
   do { \
      double _a = (a), _b = (b), _tol = (tol); \
      if (std::abs(_a - _b) > _tol) { \
         std::cerr << "FAIL: " << (message) \
                   << " (got " << _a << " vs " << _b \
                   << ", diff=" << std::abs(_a - _b) << ")\n"; \
         num_tests_failed++; \
      } else { \
         std::cout << "  PASS: " << (message) << "\n"; \
         num_tests_passed++; \
      } \
   } while (0)

// =============================================================================
// Helper: compare two vectors exactly
// =============================================================================
bool VectorsEqual(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return false; }
   for (int i = 0; i < a.Size(); i++)
   {
      if (a(i) != b(i)) { return false; }
   }
   return true;
}

double VectorMaxDiff(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return 1e30; }
   double max_diff = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      max_diff = std::max(max_diff, std::abs(a(i) - b(i)));
   }
   return max_diff;
}

// =============================================================================
// Test 1: Checkpoint file round-trip
// =============================================================================
void test_checkpoint_roundtrip()
{
   std::cout << "\n=== Test: Checkpoint File Round-Trip ===\n";

   std::string prefix = "test_ckpt_roundtrip";

   // Create test data with known values
   real_t t_write = 1.23456789012345e8;
   real_t dt_write = 4.567e3;
   int step_write = 42;
   int num_eq_write = 2;
   bool in_eq_write = true;

   int state_size = 20;
   Vector state_write(state_size);
   for (int i = 0; i < state_size; i++)
   {
      state_write(i) = 1.0 + 0.1 * i + 1e-15 * i * i;
   }

   int disp_size = 50;
   Vector disp_write(disp_size);
   for (int i = 0; i < disp_size; i++)
   {
      disp_write(i) = -2.5 + 0.05 * i;
   }

   int traction_size = 10;
   Vector traction_write(traction_size);
   for (int i = 0; i < traction_size; i++)
   {
      traction_write(i) = 25e6 + 1e5 * i;
   }

   int slip_rate_size = 10;
   Vector slip_rate_write(slip_rate_size);
   for (int i = 0; i < slip_rate_size; i++)
   {
      slip_rate_write(i) = 1e-9 + 1e-10 * i;
   }

   bool fsal_write = true;
   Vector k0_write(state_size);
   for (int i = 0; i < state_size; i++)
   {
      k0_write(i) = -1e-8 + 2e-9 * i;
   }

   // Write checkpoint
   WriteCheckpoint(prefix, t_write, dt_write,
                   step_write, num_eq_write, in_eq_write,
                   state_write, disp_write, traction_write, slip_rate_write,
                   fsal_write, k0_write, nullptr);

   // Read checkpoint
   real_t t_read, dt_read;
   int step_read, num_eq_read;
   bool in_eq_read, fsal_read;
   Vector state_read, disp_read, traction_read, slip_rate_read, k0_read;

   bool ok = ReadCheckpoint(prefix, t_read, dt_read,
                            step_read, num_eq_read, in_eq_read,
                            state_read, disp_read, traction_read,
                            slip_rate_read, fsal_read, k0_read, nullptr);

   TEST_ASSERT(ok, "ReadCheckpoint returned true");
   TEST_ASSERT_NEAR(t_read, t_write, 0.0, "Time round-trip exact");
   TEST_ASSERT_NEAR(dt_read, dt_write, 0.0, "dt round-trip exact");
   TEST_ASSERT(step_read == step_write, "Step round-trip");
   TEST_ASSERT(num_eq_read == num_eq_write, "num_eq round-trip");
   TEST_ASSERT(in_eq_read == in_eq_write, "in_eq round-trip");
   TEST_ASSERT(fsal_read == fsal_write, "FSAL flag round-trip");

   TEST_ASSERT(VectorsEqual(state_read, state_write),
               "State vector round-trip exact");
   TEST_ASSERT(VectorsEqual(disp_read, disp_write),
               "Displacement round-trip exact");
   TEST_ASSERT(VectorsEqual(traction_read, traction_write),
               "Traction round-trip exact");
   TEST_ASSERT(VectorsEqual(slip_rate_read, slip_rate_write),
               "Slip rate round-trip exact");
   TEST_ASSERT(VectorsEqual(k0_read, k0_write),
               "K0 round-trip exact");

   // Clean up
   std::remove(CheckpointFilename(prefix, 0).c_str());
}

// =============================================================================
// Test 2: Simulation restart produces identical results
// =============================================================================
void test_restart_consistency()
{
   std::cout << "\n=== Test: Simulation Restart Consistency ===\n";

   // Small mesh for fast testing
   BP2Params params;
   real_t Lx = 2.0e3;    // 2 km (tiny domain)
   real_t Lz = 4.0e3;    // 4 km
   int mesh_nx = 2;
   int mesh_nz = 5;
   real_t grading_x = 1.0;
   real_t grading_z = 1.0;

   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = Lx;
   mesh_params.Lz = Lz;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = mesh_nx;
   mesh_params.nz = mesh_nz;
   mesh_params.grading_x = grading_x;
   mesh_params.grading_z = grading_z;

   // --- Run A: continuous 20 steps ---
   Vector state_A;
   real_t t_A;
   int step_A;
   {
      auto mesh = BP2MeshGenerator::CreateGraded(mesh_params);
      int order = 1;
      AntiplaneDomainOperator<Mesh> domain(
         *mesh, order, params.mu(), params.Vp, params.Wf, DGMethod::BR2);

      FaultGeometry<Mesh> fault_geom(domain, params);

      DieterichRuinaFriction::Constants fc;
      fc.V0 = params.V0; fc.f0 = params.f0;
      fc.b = params.b;   fc.Dc = params.Dc;
      DieterichRuinaFriction friction(fc);
      AgingLaw aging;

      RateStateFaultOperator<Mesh> fault_op(
         &fault_geom, &friction, &aging, params);

      SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault_op);

      Vector state(fault_op.StateSize());
      seas_op.SetInitialCondition(state);

      DormandPrinceRK45 solver;
      solver.SetAbsTol(1e-7);
      solver.SetRelTol(1e-7);
      solver.SetDtMin(1e-6);
      solver.SetDtMax(1e8);
      solver.SetDt(1e3);
      solver.Init(seas_op);

      real_t t = 0.0;
      int step = 0;
      int target_steps = 20;

      while (step < target_steps)
      {
         real_t dt;
         bool accepted = solver.Step(seas_op, state, t, dt);
         if (accepted) { step++; }
      }

      state_A = state;
      t_A = t;
      step_A = step;
   }

   std::cout << "  Run A (continuous): t=" << t_A << " s, step=" << step_A
             << "\n";

   // --- Run B: 10 steps, checkpoint, restart, 10 more steps ---
   Vector state_B;
   real_t t_B;
   int step_B;
   {
      // Phase 1: run 10 steps and checkpoint
      std::string ckpt_prefix = "test_ckpt_restart";
      real_t ckpt_t, ckpt_dt;
      int ckpt_step;
      Vector ckpt_state, ckpt_disp, ckpt_traction, ckpt_slip_rate, ckpt_k0;
      bool ckpt_fsal;
      {
         auto mesh = BP2MeshGenerator::CreateGraded(mesh_params);
         int order = 1;
         AntiplaneDomainOperator<Mesh> domain(
            *mesh, order, params.mu(), params.Vp, params.Wf, DGMethod::BR2);

         FaultGeometry<Mesh> fault_geom(domain, params);

         DieterichRuinaFriction::Constants fc;
         fc.V0 = params.V0; fc.f0 = params.f0;
         fc.b = params.b;   fc.Dc = params.Dc;
         DieterichRuinaFriction friction(fc);
         AgingLaw aging;

         RateStateFaultOperator<Mesh> fault_op(
            &fault_geom, &friction, &aging, params);

         SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault_op);

         Vector state(fault_op.StateSize());
         seas_op.SetInitialCondition(state);

         DormandPrinceRK45 solver;
         solver.SetAbsTol(1e-7);
         solver.SetRelTol(1e-7);
         solver.SetDtMin(1e-6);
         solver.SetDtMax(1e8);
         solver.SetDt(1e3);
         solver.Init(seas_op);

         real_t t = 0.0;
         int step = 0;
         int checkpoint_step = 10;

         while (step < checkpoint_step)
         {
            real_t dt;
            bool accepted = solver.Step(seas_op, state, t, dt);
            if (accepted) { step++; }
         }

         // Save checkpoint
         Vector u_vec;
         u_vec = seas_op.GetDisplacement();
         WriteCheckpoint(ckpt_prefix, t, solver.GetDt(),
                         step, 0, false,
                         state, u_vec,
                         seas_op.GetTraction(), fault_op.GetSlipRate(),
                         solver.IsInitialized(), solver.GetK0(),
                         nullptr);

         ckpt_t = t;
         ckpt_dt = solver.GetDt();
         ckpt_step = step;
         ckpt_state = state;
         ckpt_disp = u_vec;
         ckpt_traction = seas_op.GetTraction();
         ckpt_slip_rate = fault_op.GetSlipRate();
         ckpt_fsal = solver.IsInitialized();
         ckpt_k0 = solver.GetK0();
      }

      std::cout << "  Checkpoint at: t=" << ckpt_t << " s, step=" << ckpt_step
                << "\n";

      // Phase 2: restart from checkpoint and run 10 more steps
      {
         auto mesh = BP2MeshGenerator::CreateGraded(mesh_params);
         int order = 1;
         AntiplaneDomainOperator<Mesh> domain(
            *mesh, order, params.mu(), params.Vp, params.Wf, DGMethod::BR2);

         FaultGeometry<Mesh> fault_geom(domain, params);

         DieterichRuinaFriction::Constants fc;
         fc.V0 = params.V0; fc.f0 = params.f0;
         fc.b = params.b;   fc.Dc = params.Dc;
         DieterichRuinaFriction friction(fc);
         AgingLaw aging;

         RateStateFaultOperator<Mesh> fault_op(
            &fault_geom, &friction, &aging, params);

         SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault_op);

         // Load checkpoint
         real_t t, restart_dt;
         int step, num_eq;
         bool in_eq, fsal_init;
         Vector state, disp, traction, slip_rate, k0;

         bool ok = ReadCheckpoint(ckpt_prefix, t, restart_dt,
                                  step, num_eq, in_eq,
                                  state, disp, traction, slip_rate,
                                  fsal_init, k0, nullptr);
         TEST_ASSERT(ok, "Checkpoint loaded successfully");

         // Verify loaded data matches what was written
         TEST_ASSERT_NEAR(t, ckpt_t, 0.0, "Loaded time matches");
         TEST_ASSERT_NEAR(restart_dt, ckpt_dt, 0.0, "Loaded dt matches");
         TEST_ASSERT(step == ckpt_step, "Loaded step matches");
         TEST_ASSERT(VectorsEqual(state, ckpt_state),
                     "Loaded state matches");

         // Restore simulation state
         fault_op.InitPreStress();  // Must set tau0_ before ComputeRHS
         seas_op.SetDisplacement(disp);
         fault_op.SetSlipRate(slip_rate);

         DormandPrinceRK45 solver;
         solver.SetAbsTol(1e-7);
         solver.SetRelTol(1e-7);
         solver.SetDtMin(1e-6);
         solver.SetDtMax(1e8);
         solver.SetDt(restart_dt);
         solver.Init(seas_op);

         if (fsal_init && k0.Size() > 0)
         {
            solver.RestoreFSAL(k0);
         }

         // Run 10 more steps
         int target_steps = 20;
         while (step < target_steps)
         {
            real_t dt;
            bool accepted = solver.Step(seas_op, state, t, dt);
            if (accepted) { step++; }
         }

         state_B = state;
         t_B = t;
         step_B = step;
      }

      // Clean up checkpoint files
      std::remove(CheckpointFilename(ckpt_prefix, 0).c_str());
   }

   std::cout << "  Run B (restart):    t=" << t_B << " s, step=" << step_B
             << "\n";

   // Compare results
   TEST_ASSERT(step_A == step_B, "Step counts match");
   TEST_ASSERT_NEAR(t_A, t_B, 0.0, "Final times match exactly");

   double max_diff = VectorMaxDiff(state_A, state_B);
   std::cout << "  State max diff: " << max_diff << "\n";
   TEST_ASSERT(VectorsEqual(state_A, state_B),
               "Final states match exactly (bitwise)");
}

// =============================================================================
// Test 3: Missing checkpoint file returns false
// =============================================================================
void test_missing_checkpoint()
{
   std::cout << "\n=== Test: Missing Checkpoint Returns False ===\n";

   real_t t, dt;
   int step, num_eq;
   bool in_eq, fsal;
   Vector state, disp, traction, slip_rate, k0;

   bool ok = ReadCheckpoint("nonexistent_checkpoint_xyz", t, dt,
                            step, num_eq, in_eq,
                            state, disp, traction, slip_rate,
                            fsal, k0, nullptr);
   TEST_ASSERT(!ok, "ReadCheckpoint returns false for missing file");
}

// =============================================================================
// Main
// =============================================================================
int main(int argc, char *argv[])
{
   std::cout << "Checkpoint/Restart Unit Tests\n";
   std::cout << "=============================\n";

   test_checkpoint_roundtrip();
   test_missing_checkpoint();
   test_restart_consistency();

   std::cout << "\n=============================\n";
   std::cout << "Tests passed: " << num_tests_passed << "\n";
   std::cout << "Tests failed: " << num_tests_failed << "\n";

   return num_tests_failed > 0 ? 1 : 0;
}
