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

// BP5 Parallel Smoke Test
//
// Tests the BP5 simulation pipeline without full ODE integration:
// 1. CreateBP5InlineMesh() with coarse resolution
// 2. ElasticityDomainOperator<ParMesh> with BR2
// 3. FaultGeometry<ParMesh> for BP5 parameters
// 4. RateStateFaultOperator<ParMesh, 2> with DR friction
// 5. SEAS operator: single Mult() call (= Solve + ComputeTraction + ComputeRHS)
// 6. Verify: no NaN/Inf, traction bounded, displacement bounded
//
// Usage: mpirun -np N ./seas_test_bp5_parallel_smoke

#include "mfem.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp5_params.hpp"
#include "../../common/mpi_context.hpp"

#include <iostream>
#include <cmath>
#include <memory>

using namespace mfem;
using namespace mfem::seas;

// ============================================================================
// Inline mesh creation (same as bp5_verification_full.cpp)
// ============================================================================

std::unique_ptr<Mesh> CreateBP5InlineMesh(
   int nx, int ny, int nz,
   real_t Lx, real_t Ly, real_t Lz)
{
   auto mesh = std::make_unique<Mesh>(
      Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                            Element::HEXAHEDRON,
                            2.0 * Lx, 2.0 * Ly, Lz));

   for (int i = 0; i < mesh->GetNV(); i++)
   {
      real_t *v = mesh->GetVertex(i);
      v[0] -= Lx;
      v[1] -= Ly;
   }

   const real_t tol = 1e-6 * std::max({Lx, Ly, Lz});

   for (int i = 0; i < mesh->GetNBE(); i++)
   {
      Array<int> vertices;
      mesh->GetBdrElementVertices(i, vertices);

      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int j = 0; j < vertices.Size(); j++)
      {
         const real_t *v = mesh->GetVertex(vertices[j]);
         cx += v[0]; cy += v[1]; cz += v[2];
      }
      cx /= vertices.Size();
      cy /= vertices.Size();
      cz /= vertices.Size();

      int attr;
      if (std::abs(cx - (-Lx)) < tol)      { attr = 1; }
      else if (std::abs(cx - Lx) < tol)     { attr = 2; }
      else if (std::abs(cy - Ly) < tol)     { attr = 3; }
      else if (std::abs(cy - (-Ly)) < tol)  { attr = 4; }
      else if (std::abs(cz) < tol)          { attr = 5; }
      else if (std::abs(cz - Lz) < tol)     { attr = 6; }
      else                                   { attr = 1; }

      mesh->SetBdrAttribute(i, attr);
   }

   mesh->SetAttributes();
   return mesh;
}

// ============================================================================
// Test macros
// ============================================================================

static int num_passed = 0;
static int num_failed = 0;

#define TEST_CHECK(ctx, name, cond) \
   do { \
      bool ok = (cond); \
      int local_ok = ok ? 1 : 0; \
      int global_ok = ctx.GlobalMinInt(local_ok); \
      if (ctx.IsRoot()) { \
         if (global_ok) { \
            std::cout << "  PASS: " << (name) << "\n"; \
            num_passed++; \
         } else { \
            std::cout << "  FAIL: " << (name) << "\n"; \
            num_failed++; \
         } \
      } \
   } while (0)

// ============================================================================
// Main test
// ============================================================================

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);

   if (mpi.IsRoot())
   {
      std::cout << "BP5 Parallel Smoke Test (NP=" << mpi.Size() << ")\n";
      std::cout << std::string(50, '=') << "\n";
   }

   // ========================================================================
   // Create coarse mesh
   // ========================================================================
   BP5Params params;

   // Use small mesh matching parallel_elasticity tests (which work reliably)
   // with Wf/lf sized to capture fault DOFs at these dimensions.
   int nx = 2, ny = 1, nz = 1;
   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;

   auto serial_mesh = CreateBP5InlineMesh(nx, ny, nz, Lx, Ly, Lz);

   if (mpi.IsRoot())
   {
      std::cout << "  Serial mesh: " << serial_mesh->GetNE()
                << " elements\n" << std::flush;
   }

   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   // GetGlobalNE() may be collective — call on all ranks
   long long global_ne = pmesh.GetGlobalNE();
   if (mpi.IsRoot())
   {
      std::cout << "  ParMesh: " << global_ne
                << " global elements\n";
   }

   // ========================================================================
   // Domain operator
   // ========================================================================
   int order = 1;
   // Use mesh-matching Wf/lf to ensure fault DOFs are detected on small mesh
   real_t Wf = Lz, lf_domain = 2.0 * Ly;
   ElasticityDomainOperator<ParMesh> domain(
      pmesh, order, params.lambda(), params.mu(),
      params.Vp, Wf, lf_domain, DGMethod::BR2);

   int local_fault_dofs = domain.GetNumFaultDOFs();
   int global_fault_dofs = mpi.GlobalSumInt(local_fault_dofs);

   if (mpi.IsRoot())
   {
      std::cout << "  Global fault DOFs: " << global_fault_dofs << "\n";
   }

   TEST_CHECK(mpi, "Fault DOFs detected", global_fault_dofs > 0);

   // ========================================================================
   // Fault components
   // ========================================================================
   FaultGeometry<ParMesh> fault_geom(domain, params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.L0;
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(params.b, params.V0, params.f0);

   RateStateFaultOperator<ParMesh, 2> fault_op(
      &fault_geom, &friction, &aging, params, &mpi);

   // ========================================================================
   // SEAS operator + initial condition
   // ========================================================================
   PBP5SEASOp seas_op(&domain, &fault_op, &mpi);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   real_t V_init = seas_op.GetMaxSlipRate();

   if (mpi.IsRoot())
   {
      std::cout << "  V_init = " << V_init << " m/s\n";
      std::cout << "  StateSize = " << fault_op.StateSize() << "\n";
   }

   TEST_CHECK(mpi, "Initial V_max finite",
              std::isfinite(V_init) && V_init >= 0.0);

   // ========================================================================
   // Check initial traction is bounded
   // ========================================================================
   {
      const Vector &trac = seas_op.GetTraction();
      bool trac_ok = true;
      real_t max_tau = 0.0;
      for (int i = 0; i < trac.Size(); i++)
      {
         if (!std::isfinite(trac(i)))
         {
            trac_ok = false;
            break;
         }
         max_tau = std::max(max_tau, std::abs(trac(i)));
      }
      real_t global_max_tau = mpi.GlobalMax(max_tau);

      TEST_CHECK(mpi, "Initial traction finite", trac_ok);
      TEST_CHECK(mpi, "Initial traction bounded (< 1 GPa)",
                 global_max_tau < 1e9);

      if (mpi.IsRoot())
      {
         std::cout << "  Max initial |tau| = "
                   << global_max_tau / 1e6 << " MPa\n";
      }
   }

   // ========================================================================
   // Single SEAS Mult() call — exercises Solve + ComputeTraction + ComputeRHS
   // This is the critical path that blows up on Frontera.
   // ========================================================================
   if (mpi.IsRoot())
   {
      std::cout << "\nSingle Mult() call:\n";
   }

   seas_op.SetTime(0.0);
   Vector rate(fault_op.StateSize());
   rate = 0.0;
   seas_op.Mult(state, rate);

   // Check rate vector for NaN/Inf
   {
      bool rate_ok = true;
      for (int i = 0; i < rate.Size(); i++)
      {
         if (!std::isfinite(rate(i))) { rate_ok = false; break; }
      }
      int local_ok = rate_ok ? 1 : 0;
      int global_ok = mpi.GlobalMinInt(local_ok);
      TEST_CHECK(mpi, "Rate vector finite after Mult()", global_ok == 1);
   }

   // Check traction after Mult()
   {
      const Vector &trac = seas_op.GetTraction();
      real_t max_tau = 0.0;
      bool trac_ok = true;
      for (int i = 0; i < trac.Size(); i++)
      {
         if (!std::isfinite(trac(i))) { trac_ok = false; break; }
         max_tau = std::max(max_tau, std::abs(trac(i)));
      }
      real_t global_max_tau = mpi.GlobalMax(max_tau);

      TEST_CHECK(mpi, "Post-Mult traction finite", trac_ok);
      TEST_CHECK(mpi, "Post-Mult traction bounded (< 1 GPa)",
                 global_max_tau < 1e9);

      if (mpi.IsRoot())
      {
         std::cout << "  Max |tau| after Mult = "
                   << global_max_tau / 1e6 << " MPa\n";
      }
   }

   // Check displacement after Mult()
   {
      const auto &u = seas_op.GetDisplacement();
      real_t u_max = 0.0;
      for (int i = 0; i < u.Size(); i++)
      {
         u_max = std::max(u_max, std::abs(u(i)));
      }
      real_t global_u_max = mpi.GlobalMax(u_max);

      TEST_CHECK(mpi, "Displacement bounded (< 1e6 m)",
                 global_u_max < 1e6);

      if (mpi.IsRoot())
      {
         std::cout << "  ||u||_inf = " << global_u_max << " m\n";
      }
   }

   // Check V_max after Mult()
   {
      real_t V_max = seas_op.GetMaxSlipRate();
      TEST_CHECK(mpi, "V_max finite after Mult()",
                 std::isfinite(V_max) && V_max >= 0.0);
      TEST_CHECK(mpi, "V_max bounded (< 1e4 m/s)", V_max < 1e4);

      if (mpi.IsRoot())
      {
         std::cout << "  V_max = " << V_max << " m/s\n";
      }
   }

   // ========================================================================
   // Explicit Euler step: state += dt * rate, then call Mult() again
   // Tests that a second evaluation doesn't blow up.
   // ========================================================================
   if (mpi.IsRoot())
   {
      std::cout << "\nExplicit Euler step (dt=0.1s):\n";
   }

   real_t dt_euler = 0.1;
   Vector state2(state);
   state2.Add(dt_euler, rate);

   seas_op.SetTime(dt_euler);
   Vector rate2(fault_op.StateSize());
   rate2 = 0.0;
   seas_op.Mult(state2, rate2);

   // Check after second Mult()
   {
      bool rate_ok = true;
      for (int i = 0; i < rate2.Size(); i++)
      {
         if (!std::isfinite(rate2(i))) { rate_ok = false; break; }
      }
      int local_ok = rate_ok ? 1 : 0;
      int global_ok = mpi.GlobalMinInt(local_ok);
      TEST_CHECK(mpi, "Rate2 finite after Euler step", global_ok == 1);

      const Vector &trac = seas_op.GetTraction();
      real_t max_tau = 0.0;
      for (int i = 0; i < trac.Size(); i++)
      {
         max_tau = std::max(max_tau, std::abs(trac(i)));
      }
      real_t global_max_tau = mpi.GlobalMax(max_tau);

      TEST_CHECK(mpi, "Post-Euler traction bounded (< 1 GPa)",
                 global_max_tau < 1e9);

      real_t V_max = seas_op.GetMaxSlipRate();
      TEST_CHECK(mpi, "Post-Euler V_max bounded (< 1e4 m/s)",
                 std::isfinite(V_max) && V_max < 1e4);

      if (mpi.IsRoot())
      {
         std::cout << "  Max |tau| = " << global_max_tau / 1e6
                   << " MPa, V_max = " << V_max << " m/s\n";
      }
   }

   // ========================================================================
   // Summary
   // ========================================================================
   if (mpi.IsRoot())
   {
      std::cout << "\n" << std::string(50, '=') << "\n";
      std::cout << "Results: " << num_passed << " passed, "
                << num_failed << " failed\n";
   }

   return num_failed;
}
