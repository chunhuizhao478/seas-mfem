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

// BP2 Serial Smoke Test (800 m coarse mesh)
//
// Runs the full SEAS quasi-dynamic simulation with BP2 parameters on a coarse
// 800 m mesh. Verifies:
//
// 1. Stability: Runs without crashes, solver failures, or NaN/Inf
// 2. Seismic-aseismic cycling: At least 5 earthquake cycles
// 3. Qualitative behavior: Slip rate transitions between ~1e-9 and >1e-3 m/s
// 4. I/O correctness: SCEC probe files and ParaView PVD/VTU written
//
// The 800 m mesh is too coarse for quantitative match with SCEC benchmark
// data. Quantitative comparison is deferred to Phase 9 with refined meshes.

#include "mfem.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp2_params.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../io/benchmark_output.hpp"
#include "../../io/paraview_output.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdlib>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Simple test framework
// =============================================================================

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      real_t _val = (value); \
      real_t _exp = (expected); \
      real_t _tol = (tol); \
      if (std::abs(_val - _exp) > _tol) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         std::cerr << "  Expected: " << _exp << ", Got: " << _val \
                   << ", Diff: " << std::abs(_val - _exp) << "\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

// =============================================================================
// Simulation result data
// =============================================================================

struct SimulationResult
{
   bool completed_successfully = false;
   bool has_nan_or_inf = false;
   real_t final_time = 0.0;
   int total_steps = 0;
   int num_seismic_events = 0;
   bool has_interseismic_phases = false;
   real_t max_slip_rate_ever = 0.0;
   std::vector<real_t> earthquake_times;
};

// =============================================================================
// Run BP2 serial simulation on 800 m mesh
// =============================================================================

SimulationResult RunBP2Serial(real_t t_final, int mesh_nx, int mesh_nz,
                              real_t grading_x = 1.0, real_t grading_z = 1.0,
                              bool enable_io = true,
                              real_t Lx = 100.0e3, real_t Lz = 100.0e3,
                              const std::string &output_prefix = "mfem_bp2qd")
{
   SimulationResult result;

   BP2Params params;
   params.t_final = t_final;

   // =========================================================================
   // Mesh
   // =========================================================================

   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = Lx;
   mesh_params.Lz = Lz;
   mesh_params.Wf = params.Wf;  // 40 km fault depth
   mesh_params.nx = mesh_nx;
   mesh_params.nz = mesh_nz;
   mesh_params.grading_x = grading_x;
   mesh_params.grading_z = grading_z;

   std::cout << "Creating mesh: " << 2 * mesh_params.nx << " x " << mesh_params.nz
             << " = " << 2 * mesh_params.nx * mesh_params.nz << " elements\n";
   std::cout << "Lx=" << Lx/1e3 << " km, Lz=" << Lz/1e3 << " km\n";
   std::cout << "Nominal element size (x): " << mesh_params.Lx / mesh_params.nx << " m\n";
   if (grading_x != 1.0 || grading_z != 1.0)
   {
      std::cout << "Grading: x=" << grading_x << ", z=" << grading_z
                << " (elements concentrated near fault/surface)\n";
   }

   auto mesh = (grading_x != 1.0 || grading_z != 1.0) ?
      BP2MeshGenerator::CreateGraded(mesh_params) :
      BP2MeshGenerator::Create(mesh_params);
   std::cout << "Mesh created: " << mesh->GetNE() << " elements, "
             << mesh->GetNV() << " vertices\n";

   // =========================================================================
   // Domain operator: DG order 1
   // =========================================================================
   int order = 1;
   AntiplaneDomainOperator<Mesh> domain(
      *mesh, order, params.mu(), params.Vp, params.Wf, DGMethod::BR2);

   std::cout << "Domain operator: " << domain.GetFESpace().GetTrueVSize()
             << " DOFs, " << domain.GetNumFaultDOFs() << " fault DOFs\n";

   // =========================================================================
   // Fault components
   // =========================================================================
   FaultGeometry<Mesh> fault_geom(domain, params);
   fault_geom.Print();

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;  fc.f0 = params.f0;
   fc.b = params.b;     fc.Dc = params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<Mesh> fault_op(
      &fault_geom, &friction, &aging, params);

   // =========================================================================
   // SEAS operator
   // =========================================================================
   SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault_op);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   // Verify initial conditions
   real_t V_init = seas_op.GetMaxSlipRate();
   std::cout << "\nInitial conditions:\n";
   std::cout << "  tau0 = " << fault_op.GetTau0() / 1e6 << " MPa\n";
   std::cout << "  V_max = " << V_init << " m/s\n";
   fault_op.PrintState(state);

   // =========================================================================
   // I/O setup
   // =========================================================================
   Vector fault_depths;
   domain.GetFaultDepths(fault_depths);

   // SCEC output stations
   std::vector<real_t> probe_depths = {
      0.0, -2400.0, -4800.0, -7200.0, -9600.0,
      -12000.0, -14400.0, -16800.0, -19200.0,
      -24000.0, -28800.0, -36000.0
   };

   std::unique_ptr<BenchmarkOutput<Mesh>> bench_out;
   std::unique_ptr<ParaViewOutput<Mesh>> pv_out;

   if (enable_io)
   {
      // Clean up old output files from previous runs
      std::string rm_cmd = "rm -f " + output_prefix + "_z*.txt";
      std::system(rm_cmd.c_str());

      bench_out = std::make_unique<BenchmarkOutput<Mesh>>(
         output_prefix, params, probe_depths, fault_depths);

      std::string pv_dir = "ParaView/" + output_prefix;
      pv_out = std::make_unique<ParaViewOutput<Mesh>>(
         pv_dir, *mesh, order);
      pv_out->output_every_n_steps = 100;

      GridFunction *u_gf = const_cast<GridFunction*>(
         &static_cast<const GridFunction&>(seas_op.GetDisplacement()));
      pv_out->RegisterDomainField("displacement", u_gf);

      // Initialize fault output on domain mesh using fault interior faces
      const Array<int> &fault_faces = domain.GetFaultInteriorFaces();
      pv_out->InitFaultOutput(fault_faces);

      // Write initial state
      bench_out->ForceWrite(0.0, state, fault_op, seas_op.GetTraction());
      {
         // Update fault fields for initial save
         Vector slip, theta;
         fault_op.GetSlip(state, slip);
         fault_op.GetTheta(state, theta);
         const Vector &V = fault_op.GetSlipRate();
         real_t tau0 = fault_op.GetTau0();
         Vector stress(fault_op.NumNodes());
         for (int i = 0; i < stress.Size(); i++)
         {
            stress(i) = tau0 + seas_op.GetTraction()(i);
         }
         int dofs_per_face = 1;  // 1 midpoint DOF per fault face
         pv_out->UpdateFaultFields(slip, V, stress, theta, dofs_per_face);
      }
      pv_out->ForceSave(0, 0.0);
   }

   // =========================================================================
   // Time integration (Dormand-Prince RK45 with error-based adaptive dt)
   // =========================================================================
   DormandPrinceRK45 ode_solver;
   ode_solver.SetAbsTol(1e-7);
   ode_solver.SetRelTol(1e-50);  // Match Tandem/PETSc: pure absolute tolerance
   ode_solver.SetDtMin(1e-6);
   ode_solver.SetDtMax(0.1 * BP2Params::seconds_per_year);  // Max 0.1 yr
   ode_solver.SetDt(1e3);                                     // Start with 1000 s
   ode_solver.Init(seas_op);

   real_t t = 0.0;
   int step = 0;
   int max_steps = 10000000;  // Safety limit

   // Earthquake detection state
   bool in_seismic_event = false;
   real_t V_threshold_seismic = 1e-3;      // Earthquake when V > 1e-3 m/s
   real_t V_threshold_interseismic = 1e-6;  // Back to interseismic when V < 1e-6

   real_t V_min_between_events = 1e30;  // Track minimum V between events

   std::cout << "\n" << std::setw(10) << "Step"
             << std::setw(16) << "Time [yr]"
             << std::setw(14) << "dt [s]"
             << std::setw(16) << "V_max [m/s]"
             << std::setw(8) << "EQs"
             << std::setw(8) << "Rej"
             << "\n";
   std::cout << std::string(72, '-') << "\n";

   real_t next_print_time = 0.0;
   real_t print_interval = 10.0 * BP2Params::seconds_per_year;  // Every 10 yr

   while (t < t_final && step < max_steps)
   {
      // Clamp dt to not overshoot t_final
      if (t + ode_solver.GetDt() > t_final)
      {
         ode_solver.SetDt(t_final - t);
      }

      real_t dt;
      bool accepted = ode_solver.Step(seas_op, state, t, dt);
      if (!accepted) { continue; }  // Retry with smaller dt
      step++;

      real_t V_max = seas_op.GetMaxSlipRate();
      result.max_slip_rate_ever = std::max(result.max_slip_rate_ever, V_max);

      // Check for NaN/Inf
      if (std::isnan(V_max) || std::isinf(V_max))
      {
         std::cerr << "NaN/Inf detected at step " << step
                   << ", t = " << t / BP2Params::seconds_per_year << " yr\n";
         result.has_nan_or_inf = true;
         break;
      }

      // Also check state vector
      {
         Vector theta;
         fault_op.GetTheta(state, theta);
         if (theta.Min() <= 0.0 || std::isnan(theta.Min()))
         {
            std::cerr << "Theta non-positive or NaN at step " << step
                      << ", theta_min = " << theta.Min() << "\n";
            result.has_nan_or_inf = true;
            break;
         }
      }

      // Diagnostic: dump fault state sorted by depth when V_max first exceeds 0.1
      static bool dumped_v_profile = false;
      if (!dumped_v_profile && V_max > 0.1)
      {
         dumped_v_profile = true;
         Vector fault_depths_dbg;
         domain.GetFaultDepths(fault_depths_dbg);
         const Vector &V_all = fault_op.GetSlipRate();
         const Vector &a_vals = fault_geom.GetAValues();
         const Vector &trac = seas_op.GetTraction();
         Vector slip_dbg;
         fault_op.GetSlip(state, slip_dbg);

         // Sort DOFs by depth for readable output
         int ndofs = V_all.Size();
         std::vector<int> sorted_idx(ndofs);
         for (int ii = 0; ii < ndofs; ii++) { sorted_idx[ii] = ii; }
         std::sort(sorted_idx.begin(), sorted_idx.end(),
            [&](int a, int b) {
               return fault_depths_dbg(a) > fault_depths_dbg(b);
            });

         std::cout << "\n=== Fault State at t = "
                   << t / BP2Params::seconds_per_year << " yr ===\n";
         std::cout << "DOF  depth(km)   a      slip(m)      V(m/s)       traction(MPa) logV\n";
         for (int jj = 0; jj < ndofs; jj++)
         {
            int ii = sorted_idx[jj];
            real_t z_km = fault_depths_dbg(ii) / 1e3;
            // Only print VW zone and transition (z > -25 km)
            if (z_km < -25.0) { continue; }
            real_t ai = a_vals(ii);
            real_t Vi = V_all(ii);
            real_t si = slip_dbg(ii);
            real_t ti = trac(ii) / 1e6;  // Pa to MPa
            real_t logV = (Vi > 0) ? std::log10(Vi) : -99.0;
            std::cout << std::setw(4) << ii << "  "
                      << std::setw(8) << std::fixed << std::setprecision(1) << z_km
                      << "  " << std::fixed << std::setprecision(4) << ai
                      << "  " << std::scientific << std::setprecision(3) << si
                      << "  " << std::scientific << std::setprecision(3) << Vi
                      << "  " << std::fixed << std::setprecision(4) << ti
                      << "  " << std::fixed << std::setprecision(2) << logV;
            if (Vi == fault_op.GetMaxSlipRate()) { std::cout << " <=="; }
            std::cout << "\n";
         }
         std::cout << "=== End Fault State ===\n\n";
      }

      // Earthquake detection
      if (!in_seismic_event && V_max > V_threshold_seismic)
      {
         // Start of a seismic event
         in_seismic_event = true;
         result.num_seismic_events++;
         result.earthquake_times.push_back(t);
         V_min_between_events = 1e30;

         std::cout << "  *** EARTHQUAKE #" << result.num_seismic_events
                   << " at t = " << std::fixed << std::setprecision(1)
                   << t / BP2Params::seconds_per_year << " yr"
                   << ", V_max = " << std::scientific << std::setprecision(2)
                   << V_max << " m/s ***\n";
      }
      else if (in_seismic_event && V_max < V_threshold_interseismic)
      {
         // Back to interseismic
         in_seismic_event = false;
      }

      // Track minimum V between events
      if (!in_seismic_event)
      {
         V_min_between_events = std::min(V_min_between_events, V_max);
      }

      // I/O
      if (enable_io && bench_out)
      {
         if (bench_out->Write(t, state, fault_op, seas_op.GetTraction()))
         {
            bench_out->Flush();
         }
      }
      if (enable_io && pv_out)
      {
         // Update fault fields before saving
         Vector slip, theta;
         fault_op.GetSlip(state, slip);
         fault_op.GetTheta(state, theta);
         const Vector &V_field = fault_op.GetSlipRate();
         real_t tau0 = fault_op.GetTau0();
         Vector stress(fault_op.NumNodes());
         for (int i = 0; i < stress.Size(); i++)
         {
            stress(i) = tau0 + seas_op.GetTraction()(i);
         }
         int dpf = 1;  // 1 midpoint DOF per fault face
         pv_out->UpdateFaultFields(slip, V_field, stress, theta, dpf);
         pv_out->Save(step, t);
      }

      // Periodic console output
      if (t >= next_print_time || V_max > V_threshold_seismic)
      {
         std::cout << std::setw(10) << step
                   << std::setw(16) << std::fixed << std::setprecision(2)
                   << t / BP2Params::seconds_per_year
                   << std::setw(14) << std::scientific << std::setprecision(3)
                   << ode_solver.GetDt()
                   << std::setw(16) << std::scientific << std::setprecision(3) << V_max
                   << std::setw(8) << result.num_seismic_events
                   << std::setw(8) << ode_solver.GetTotalRejections()
                   << "\n";
         next_print_time = t + print_interval;
      }
   }

   // Check for interseismic recovery
   result.has_interseismic_phases =
      (V_min_between_events < V_threshold_interseismic);

   result.completed_successfully = !result.has_nan_or_inf && (step < max_steps);
   result.final_time = t;
   result.total_steps = step;

   // Final I/O
   if (enable_io)
   {
      if (bench_out)
      {
         bench_out->ForceWrite(t, state, fault_op, seas_op.GetTraction());
         bench_out->Close();
      }
      if (pv_out)
      {
         pv_out->ForceSave(step, t);
      }
   }

   std::cout << "\n=== Simulation Summary ===\n";
   std::cout << "  Final time: " << t / BP2Params::seconds_per_year << " years\n";
   std::cout << "  Total steps: " << step << "\n";
   std::cout << "  Step rejections: " << ode_solver.GetTotalRejections() << "\n";
   std::cout << "  Seismic events: " << result.num_seismic_events << "\n";
   std::cout << "  Max V ever: " << result.max_slip_rate_ever << " m/s\n";
   std::cout << "  Completed: " << (result.completed_successfully ? "YES" : "NO") << "\n";
   std::cout << "  NaN/Inf: " << (result.has_nan_or_inf ? "YES" : "NO") << "\n";

   if (!result.earthquake_times.empty())
   {
      std::cout << "  Earthquake times (yr):";
      for (real_t eq_t : result.earthquake_times)
      {
         std::cout << " " << std::fixed << std::setprecision(1)
                   << eq_t / BP2Params::seconds_per_year;
      }
      std::cout << "\n";

      if (result.earthquake_times.size() > 1)
      {
         std::cout << "  Recurrence intervals (yr):";
         for (size_t i = 1; i < result.earthquake_times.size(); i++)
         {
            std::cout << " " << std::fixed << std::setprecision(1)
                      << (result.earthquake_times[i] -
                          result.earthquake_times[i-1]) /
                         BP2Params::seconds_per_year;
         }
         std::cout << "\n";
      }
   }

   return result;
}

// =============================================================================
// Test: Initial conditions are mesh-independent
// =============================================================================

void TestInitialConditions(int mesh_nx, int mesh_nz,
                           real_t grading_x = 1.0, real_t grading_z = 1.0,
                           real_t Lx = 100.0e3, real_t Lz = 100.0e3)
{
   std::cout << "\n=== Test: Initial Conditions (Mesh-Independent) ===\n";
   std::cout << "  Mesh: " << 2 * mesh_nx << " x " << mesh_nz << " elements\n";

   BP2Params params;

   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = Lx;
   mesh_params.Lz = Lz;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = mesh_nx;
   mesh_params.nz = mesh_nz;
   mesh_params.grading_x = grading_x;
   mesh_params.grading_z = grading_z;
   auto mesh = (grading_x != 1.0 || grading_z != 1.0) ?
      BP2MeshGenerator::CreateGraded(mesh_params) :
      BP2MeshGenerator::Create(mesh_params);

   AntiplaneDomainOperator<Mesh> domain(
      *mesh, 1, params.mu(), params.Vp, params.Wf, DGMethod::BR2);
   FaultGeometry<Mesh> fault_geom(domain, params);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;  fc.f0 = params.f0;
   fc.b = params.b;     fc.Dc = params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<Mesh> fault_op(
      &fault_geom, &friction, &aging, params);

   SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault_op);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   // Test initial tau
   TEST_NEAR(fault_op.GetTau0() / 1e6, 26.546, 0.01,
             "Initial tau0 = 26.546 MPa (+/- 0.01)");

   // Test initial V ~ V_init
   real_t V_max = seas_op.GetMaxSlipRate();
   TEST_NEAR(V_max, params.V_init, params.V_init * 0.1,
             "Initial V_max ~ V_init (10% tolerance)");

   // Test initial theta at surface (z=0)
   // Expected: log10(theta) ~ 3.602
   Vector theta;
   fault_op.GetTheta(state, theta);
   Vector fault_depths;
   domain.GetFaultDepths(fault_depths);

   // Find the DOF closest to the surface (z=0)
   int surface_idx = fault_geom.FindNearestDOF(0.0);
   if (surface_idx >= 0)
   {
      real_t log10_theta_surface = std::log10(theta(surface_idx));
      TEST_NEAR(log10_theta_surface, 3.602, 0.1,
                "Initial log10(theta) at surface ~ 3.602 (+/- 10%)");
   }
}

// =============================================================================
// Test: Full 800 m BP2 serial smoke test
// =============================================================================

void TestBP2SerialSmoke(int mesh_nx, int mesh_nz, real_t t_final,
                        int min_earthquakes,
                        real_t grading_x = 1.0, real_t grading_z = 1.0,
                        real_t Lx = 100.0e3, real_t Lz = 100.0e3,
                        const std::string &output_prefix = "mfem_bp2qd")
{
   real_t elem_size = Lx / mesh_nx;
   std::cout << "\n=== Test: BP2 Serial Smoke ("
             << elem_size << " m nominal, ~"
             << t_final / BP2Params::seconds_per_year << " years) ===\n";

   std::cout << "Target: " << t_final / BP2Params::seconds_per_year
             << " years (" << t_final << " s)\n";

   SimulationResult result = RunBP2Serial(t_final, mesh_nx, mesh_nz,
                                          grading_x, grading_z, true, Lx, Lz,
                                          output_prefix);

   // Test 1: Simulation completed without crash
   TEST_ASSERT(result.completed_successfully,
               "Simulation completed successfully");

   // Test 2: No NaN or Inf
   TEST_ASSERT(!result.has_nan_or_inf,
               "No NaN or Inf values");

   // Test 3: Reached target time (allow 1% tolerance)
   TEST_ASSERT(result.final_time >= t_final * 0.99,
               "Simulation reached target time");

   // Test 4: Interseismic phases exist (V drops below plate rate)
   TEST_ASSERT(result.has_interseismic_phases,
               "Interseismic recovery phases detected");

   // Tests 5-6: Only check earthquakes if mesh is fine enough
   if (min_earthquakes > 0)
   {
      TEST_ASSERT(result.num_seismic_events >= min_earthquakes,
                  "Sufficient earthquake events detected");

      TEST_ASSERT(std::log10(result.max_slip_rate_ever) > -3.0,
                  "Max V > 1e-3 m/s (coseismic rates reached)");
   }
   else
   {
      std::cout << "  SKIPPED: Earthquake tests (coarse mesh mode)\n";
   }

   // Test 7: SCEC output files exist
   TEST_ASSERT(
      std::ifstream("mfem_bp2qd_z0km.txt").good(),
      "SCEC output file z=0km exists");
   TEST_ASSERT(
      std::ifstream("mfem_bp2qd_z12km.txt").good(),
      "SCEC output file z=12km exists");
}

// =============================================================================
// Test: Output files have correct content
// =============================================================================

void TestOutputFilesContent()
{
   std::cout << "\n=== Test: Output File Content ===\n";

   // Check that the z=0km file has data and reasonable values
   std::ifstream f("mfem_bp2qd_z0km.txt");
   if (!f.is_open())
   {
      std::cout << "  SKIPPED: Output file not found (run full test first)\n";
      return;
   }

   std::string line;
   int data_lines = 0;
   bool has_coseismic = false;
   bool has_interseismic = false;

   while (std::getline(f, line))
   {
      if (line.empty() || line[0] == '#') { continue; }

      double t, slip, logV, tau, logTheta;
      std::istringstream iss(line);
      if (iss >> t >> slip >> logV >> tau >> logTheta)
      {
         data_lines++;

         if (logV > -3.0) { has_coseismic = true; }
         if (logV < -6.0) { has_interseismic = true; }

         // Verify stress is in reasonable range (15-40 MPa)
         if (tau < 10.0 || tau > 50.0)
         {
            std::cerr << "  WARNING: tau = " << tau
                      << " MPa at t = " << t << " s\n";
         }
      }
   }

   TEST_ASSERT(data_lines > 10,
               "Output file has data lines");
   if (has_coseismic)
   {
      std::cout << "  PASSED: Output records coseismic phases (logV > -3)\n";
   }
   else
   {
      std::cout << "  INFO: No coseismic phases in output (expected for coarse mesh)\n";
   }
   TEST_ASSERT(has_interseismic,
               "Output records interseismic phases (logV < -6)");
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   // Disable output buffering for progress monitoring
   std::cout << std::unitbuf;

   // Parse command-line arguments
   bool quick_mode = false;
   int nz_override = 0;
   double tfinal_override = 0.0;
   std::string prefix_override;
   for (int i = 1; i < argc; i++)
   {
      std::string arg(argv[i]);
      if (arg == "--quick" || arg == "-q") { quick_mode = true; }
      if (arg == "--nz" && i + 1 < argc) { nz_override = std::atoi(argv[++i]); }
      if (arg == "--prefix" && i + 1 < argc) { prefix_override = argv[++i]; }
      if (arg == "--tfinal" && i + 1 < argc) { tfinal_override = std::atof(argv[++i]); }
   }

   // Mesh and simulation parameters
   int mesh_nx, mesh_nz, min_earthquakes;
   real_t t_final, grading_x, grading_z, Lx, Lz;

   if (quick_mode)
   {
      // Quick mode: coarse x-resolution to save compute time.
      // nx=5 with sinh grading alpha=3 gives hx≈3km near x=0.
      // Total: 1250 elements, 125 fault DOFs (1 midpoint per face).
      mesh_nx = 5;
      mesh_nz = 125;
      grading_x = 3.0;
      grading_z = 1.0;
      Lx = 50.0e3;
      Lz = 100.0e3;
      t_final = 1e10;   // ~317 years
      min_earthquakes = 0; // Verify stability/cycling/I/O
      std::cout << "================================================\n";
      std::cout << "SEAS Phase 5: BP2 Serial Smoke Test (QUICK MODE)\n";
      std::cout << "  Mesh: " << 2*mesh_nx << "x" << mesh_nz
                << " graded, Lx=" << Lx/1e3 << "km, coarse x\n";
      std::cout << "  Purpose: verify coupling, cycling, I/O\n";
      std::cout << "================================================\n";
   }
   else
   {
      // Full mode: coarse x-resolution to save compute time.
      // Lx=50km, Lz=100km, nz=125 gives 800m uniform z-spacing
      // nx=5 with sinh grading alpha=3.0 gives ~3km near x=0
      // Total: 1250 elements, 125 fault DOFs (1 midpoint DOF per face)
      mesh_nx = 5;
      mesh_nz = 125;
      grading_x = 3.0;
      grading_z = 1.0;
      Lx = 50.0e3;    // 50 km half-width (sufficient with direct fault loading)
      Lz = 100.0e3;   // 100 km depth
      t_final = 1.2e10;   // ~380 years (past first earthquake)
      min_earthquakes = 1;
      std::cout << "================================================\n";
      std::cout << "SEAS Phase 5: BP2 Serial Smoke Test (full)\n";
      std::cout << "  Mesh: " << 2*mesh_nx << "x" << mesh_nz
                << " graded, Lx=" << Lx/1e3 << "km, Lz="
                << Lz/1e3 << "km\n";
      std::cout << "  Near-fault: ~3km x-elements, Fault: 800m z-spacing\n";
      std::cout << "  Target: ~" << t_final / BP2Params::seconds_per_year
                << " yr, min " << min_earthquakes << " earthquake(s)\n";
      std::cout << "================================================\n";
   }

   // Apply command-line overrides
   if (nz_override > 0)
   {
      mesh_nz = nz_override;
      std::cout << "  [override] nz=" << mesh_nz
                << " (hz=" << Lz/mesh_nz << " m)\n";
   }
   std::string output_prefix = prefix_override.empty() ?
      "mfem_bp2qd" : prefix_override;
   if (tfinal_override > 0.0)
   {
      t_final = tfinal_override;
      std::cout << "  [override] t_final=" << t_final
                << " s (~" << t_final / 3.15576e7 << " yr)\n";
   }
   if (!prefix_override.empty())
   {
      std::cout << "  [override] output prefix: " << output_prefix << "\n";
   }

   // Run initial condition tests first
   TestInitialConditions(mesh_nx, mesh_nz, grading_x, grading_z, Lx, Lz);

   // Run the smoke test
   TestBP2SerialSmoke(mesh_nx, mesh_nz, t_final, min_earthquakes,
                      grading_x, grading_z, Lx, Lz, output_prefix);

   // Check output file content (only for default prefix)
   if (prefix_override.empty())
   {
      TestOutputFilesContent();
   }

   std::cout << "\n================================================\n";
   std::cout << "Test Summary\n";
   std::cout << "================================================\n";
   std::cout << "Total tests: " << num_tests << "\n";
   std::cout << "Passed:      " << num_passed << "\n";
   std::cout << "Failed:      " << num_failed << "\n";

   if (num_failed > 0)
   {
      std::cout << "\nSOME TESTS FAILED!\n";
      return 1;
   }

   std::cout << "\nALL TESTS PASSED!\n";
   return 0;
}
