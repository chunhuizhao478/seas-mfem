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

// Unit tests for Phase 5: I/O classes
//
// Tests for:
// - ProbeOutput: file writing
// - ProbeInterpolator: interpolation accuracy
// - BenchmarkOutput: SCEC format output
// - ParaViewOutput: PVD/VTU output

#include "mfem.hpp"
#include "../../io/probe_output.hpp"
#include "../../io/benchmark_output.hpp"
#include "../../io/paraview_output.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../solver/time_stepper.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp2_params.hpp"
#include "../../domain/bp2_mesh.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <cstdio>
#include <sstream>

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
// Helper: check file exists and is non-empty
// =============================================================================

bool FileExistsAndNonEmpty(const std::string &path)
{
   std::ifstream f(path);
   if (!f.is_open()) { return false; }
   f.seekg(0, std::ios::end);
   return f.tellg() > 0;
}

// =============================================================================
// Test 1: ProbeOutput writes a file correctly
// =============================================================================

void TestProbeOutput()
{
   std::cout << "\n=== Test: ProbeOutput ===\n";

   std::string filename = "test_probe_output.txt";

   // Write some data
   {
      std::vector<std::string> columns = {"time(s)", "value1", "value2"};
      ProbeOutput probe(filename, columns);

      TEST_ASSERT(probe.Good(), "ProbeOutput file is open");

      probe.WriteStep({0.0, 1.0, 2.0});
      probe.WriteStep({1.0, 3.0, 4.0});
      probe.WriteStep({2.0, 5.0, 6.0});
      probe.Close();
   }

   // Verify the file
   TEST_ASSERT(FileExistsAndNonEmpty(filename), "Output file exists and is non-empty");

   // Read and verify contents
   {
      std::ifstream f(filename);
      std::string line;

      // First line should be header
      std::getline(f, line);
      TEST_ASSERT(line.find("Columns:") != std::string::npos,
                  "Header contains column names");

      // Second line should be first data row
      std::getline(f, line);
      TEST_ASSERT(line.find("0.000000000000E+00") != std::string::npos,
                  "First data row starts with t=0");

      // Count total data lines
      int count = 1; // Already read one
      while (std::getline(f, line))
      {
         if (!line.empty()) { count++; }
      }
      TEST_ASSERT(count == 3, "File has 3 data rows");
   }

   // Cleanup
   std::remove(filename.c_str());
}

// =============================================================================
// Test 2: ProbeInterpolator accuracy
// =============================================================================

void TestProbeInterpolator()
{
   std::cout << "\n=== Test: ProbeInterpolator ===\n";

   // Create a simple fault depth array: z = -40, -30, -20, -10, 0 km
   Vector fault_depths(5);
   fault_depths(0) = -40000.0;
   fault_depths(1) = -30000.0;
   fault_depths(2) = -20000.0;
   fault_depths(3) = -10000.0;
   fault_depths(4) = 0.0;

   // Create probe depths (some match, some require interpolation)
   std::vector<real_t> probe_depths = {0.0, -15000.0, -40000.0, -25000.0};

   ProbeInterpolator interp(fault_depths, probe_depths);

   TEST_ASSERT(interp.NumProbes() == 4, "4 probes created");

   // Create a linear field: f(z) = z / 1000 (depth in km)
   Vector field(5);
   for (int i = 0; i < 5; i++)
   {
      field(i) = fault_depths(i) / 1000.0;
   }

   // Interpolate
   Vector probe_values;
   interp.Interpolate(field, probe_values);

   // Check values
   TEST_NEAR(probe_values(0), 0.0, 1e-10,
             "Probe at z=0 km: exact match");
   TEST_NEAR(probe_values(1), -15.0, 1e-10,
             "Probe at z=-15 km: linear interpolation");
   TEST_NEAR(probe_values(2), -40.0, 1e-10,
             "Probe at z=-40 km: exact match");
   TEST_NEAR(probe_values(3), -25.0, 1e-10,
             "Probe at z=-25 km: linear interpolation");

   // Test InterpolateOne
   real_t val = interp.InterpolateOne(field, 1);
   TEST_NEAR(val, -15.0, 1e-10, "InterpolateOne at z=-15 km");
}

// =============================================================================
// Test 3: BenchmarkOutput SCEC filename convention
// =============================================================================

void TestBenchmarkOutputFilenames()
{
   std::cout << "\n=== Test: BenchmarkOutput SCEC Filenames ===\n";

   // Test the adaptive output interval computation
   TEST_NEAR(BenchmarkOutput<Mesh>::OutputInterval(1e-10),
             0.01 * BP2Params::seconds_per_year, 1.0,
             "Interseismic output interval ~ 0.01 yr");
   TEST_NEAR(BenchmarkOutput<Mesh>::OutputInterval(1e-4),
             0.1, 1e-10,
             "Nucleation output interval = 0.1 s");
   TEST_NEAR(BenchmarkOutput<Mesh>::OutputInterval(1e-2),
             0.001, 1e-10,
             "Coseismic output interval = 0.001 s");
}

// =============================================================================
// Test 4: BenchmarkOutput integration test
// =============================================================================

void TestBenchmarkOutputIntegration()
{
   std::cout << "\n=== Test: BenchmarkOutput Integration ===\n";

   BP2Params params;

   // Create a small mesh and run components
   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = 50.0e3;
   mesh_params.Lz = 50.0e3;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = 4;
   mesh_params.nz = 8;
   auto mesh = BP2MeshGenerator::Create(mesh_params);

   AntiplaneDomainOperator<Mesh> domain(
      *mesh, 1, params.mu(), params.Vp, params.Wf);

   FaultGeometry<Mesh> fault_geom(domain, params);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;  fc.f0 = params.f0;
   fc.b = params.b;     fc.Dc = params.Dc;
   DieterichRuinaFriction friction(fc);
   AgingLaw aging;

   RateStateFaultOperator<Mesh> fault_op(
      &fault_geom, &friction, &aging, params);

   SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault_op);

   // Initialize
   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   // Set up benchmark output with a couple of probe depths
   Vector fault_depths;
   domain.GetFaultDepths(fault_depths);

   std::vector<real_t> probe_depths = {0.0, -12000.0};
   std::string prefix = "test_mfem_bp2qd";
   BenchmarkOutput<Mesh> bench_out(prefix, params, probe_depths, fault_depths);

   // Write initial state
   bench_out.ForceWrite(0.0, state, fault_op, seas_op.GetTraction());
   bench_out.Flush();

   // Verify files were created
   TEST_ASSERT(FileExistsAndNonEmpty("test_mfem_bp2qd_z0km.txt"),
               "SCEC file for z=0km created");
   TEST_ASSERT(FileExistsAndNonEmpty("test_mfem_bp2qd_z12km.txt"),
               "SCEC file for z=12km created");

   // Verify content of z=0km file
   {
      std::ifstream f("test_mfem_bp2qd_z0km.txt");
      std::string line;
      std::getline(f, line); // description header (probe location)
      std::getline(f, line); // column header
      std::getline(f, line); // first data line

      // Should have 5 columns: time, slip, log10(V), tau(MPa), log10(theta)
      std::istringstream iss(line);
      double t, slip, logV, tau, logTheta;
      iss >> t >> slip >> logV >> tau >> logTheta;

      TEST_NEAR(t, 0.0, 1e-10, "Initial time = 0");
      TEST_NEAR(slip, 0.0, 1e-10, "Initial slip = 0");
      TEST_ASSERT(logV < -7.0, "Initial log10(V) < -7 (near V_init=1e-9)");
      TEST_ASSERT(tau > 20.0 && tau < 35.0,
                  "Initial tau in [20, 35] MPa range");
      TEST_ASSERT(logTheta > 2.0 && logTheta < 15.0,
                  "Initial log10(theta) in reasonable range");
   }

   bench_out.Close();

   // Cleanup
   std::remove("test_mfem_bp2qd_z0km.txt");
   std::remove("test_mfem_bp2qd_z12km.txt");
}

// =============================================================================
// Test 5: ParaViewOutput creates PVD/VTU files
// =============================================================================

void TestParaViewOutput()
{
   std::cout << "\n=== Test: ParaViewOutput ===\n";

   // Create a simple mesh
   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = 10.0e3;
   mesh_params.Lz = 10.0e3;
   mesh_params.nx = 4;
   mesh_params.nz = 4;
   auto mesh = BP2MeshGenerator::Create(mesh_params);

   int order = 1;
   DG_FECollection fec(order, 2, BasisType::GaussLobatto);
   FiniteElementSpace fes(mesh.get(), &fec);

   GridFunction u(&fes);
   u = 0.0;

   // Create ParaView output (single collection, no _domain suffix)
   ParaViewOutput<Mesh> pv("test_pv_bp2", *mesh, order);
   pv.RegisterDomainField("displacement", &u);

   // Save at t=0
   pv.ForceSave(0, 0.0);

   // Set some displacement and save again
   u = 1.0;
   pv.ForceSave(1, 100.0);

   // The PVD file should exist (MFEM nests: prefix/prefix/prefix.pvd)
   TEST_ASSERT(FileExistsAndNonEmpty("test_pv_bp2.pvd") ||
               FileExistsAndNonEmpty("test_pv_bp2/test_pv_bp2.pvd") ||
               FileExistsAndNonEmpty("test_pv_bp2/test_pv_bp2/test_pv_bp2.pvd"),
               "ParaView PVD file created");

   std::cout << "  ParaView output test completed.\n";
}

// =============================================================================
// Test 5b: ParaViewOutput combined domain + fault fields
// =============================================================================

void TestParaViewCombinedOutput()
{
   std::cout << "\n=== Test: ParaView Combined Domain+Fault Output ===\n";

   BP2Params params;

   // Create a small mesh with a fault at x=0
   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = 10.0e3;
   mesh_params.Lz = 10.0e3;
   mesh_params.Wf = 10.0e3;  // Fault covers entire depth
   mesh_params.nx = 4;
   mesh_params.nz = 4;
   auto mesh = BP2MeshGenerator::Create(mesh_params);

   int order = 1;

   // Create domain operator to get fault faces
   AntiplaneDomainOperator<Mesh> domain(
      *mesh, order, params.mu(), params.Vp, mesh_params.Wf);

   const Array<int> &fault_faces = domain.GetFaultInteriorFaces();
   int num_fault_dofs = domain.GetNumFaultDOFs();
   // Per-face DOFs: for order 1, each face has 2 independent DOFs (not shared).
   // Total DOFs = 2 * num_faces.
   int dofs_per_face = 2;  // Per-face DOFs

   TEST_ASSERT(fault_faces.Size() > 0, "Fault interior faces found");
   TEST_ASSERT(num_fault_dofs == 2 * fault_faces.Size(),
               "Fault DOFs = 2 * num_faces (per-face DOFs)");

   // Create displacement GridFunction
   DG_FECollection fec(order, 2, BasisType::GaussLobatto);
   FiniteElementSpace fes(mesh.get(), &fec);
   GridFunction u(&fes);
   u = 0.0;

   // Create combined ParaView output
   ParaViewOutput<Mesh> pv("test_pv_combined", *mesh, order);
   pv.RegisterDomainField("displacement", &u);
   pv.InitFaultOutput(fault_faces);

   TEST_ASSERT(pv.HasFaultOutput(), "Fault output initialized");
   TEST_ASSERT(pv.GetNumFaultFaces() == fault_faces.Size(),
               "Fault face count matches");

   // Create mock fault field data
   Vector slip(num_fault_dofs), slip_rate(num_fault_dofs);
   Vector stress(num_fault_dofs), state(num_fault_dofs);

   slip = 0.5;           // 0.5 m slip
   slip_rate = 1e-6;     // 1 um/s
   stress = 26.0e6;      // 26 MPa
   state = 4000.0;       // 4000 s

   // Update and save
   pv.UpdateFaultFields(slip, slip_rate, stress, state, dofs_per_face);
   pv.ForceSave(0, 0.0);

   // Change fields and save again
   slip = 1.0;
   slip_rate = 1e-3;
   pv.UpdateFaultFields(slip, slip_rate, stress, state, dofs_per_face);
   pv.ForceSave(1, 100.0);

   // Verify single PVD created (not separate _domain and _fault)
   TEST_ASSERT(FileExistsAndNonEmpty("test_pv_combined.pvd") ||
               FileExistsAndNonEmpty("test_pv_combined/test_pv_combined.pvd") ||
               FileExistsAndNonEmpty("test_pv_combined/test_pv_combined/test_pv_combined.pvd"),
               "Combined PVD file created");

   // Verify VTU has all 5 fields
   // Find a VTU file
   std::string vtu_path;
   for (const auto &candidate : {
      std::string("test_pv_combined/test_pv_combined/Cycle000000/proc000000.vtu"),
      std::string("test_pv_combined/Cycle000000/proc000000.vtu")
   })
   {
      std::ifstream f(candidate);
      if (f.good()) { vtu_path = candidate; break; }
   }

   if (!vtu_path.empty())
   {
      // Read VTU and check for field names
      std::ifstream vtu(vtu_path);
      std::string content((std::istreambuf_iterator<char>(vtu)),
                          std::istreambuf_iterator<char>());

      TEST_ASSERT(content.find("displacement") != std::string::npos,
                  "VTU contains displacement field");
      TEST_ASSERT(content.find("slip_rate") != std::string::npos,
                  "VTU contains slip_rate field");
      TEST_ASSERT(content.find("shear_stress") != std::string::npos,
                  "VTU contains shear_stress field");
      TEST_ASSERT(content.find("state_variable") != std::string::npos,
                  "VTU contains state_variable field");
      TEST_ASSERT(content.find("\"slip\"") != std::string::npos,
                  "VTU contains slip field");
   }
   else
   {
      std::cout << "  SKIPPED: VTU file not found for field verification\n";
   }

   std::cout << "  Combined ParaView output test completed.\n";
}

// =============================================================================
// Test 6: Adaptive ParaView output interval
// =============================================================================

void TestParaViewOutputInterval()
{
   std::cout << "\n=== Test: ParaView Adaptive Output Interval ===\n";

   TEST_NEAR(ParaViewOutput<Mesh>::OutputInterval(1e-10),
             1.0 * BP2Params::seconds_per_year, 1.0,
             "PV interseismic interval ~ 1 yr");
   TEST_NEAR(ParaViewOutput<Mesh>::OutputInterval(1e-4),
             1.0, 1e-10,
             "PV nucleation interval = 1.0 s");
   TEST_NEAR(ParaViewOutput<Mesh>::OutputInterval(1e-2),
             0.01, 1e-10,
             "PV coseismic interval = 0.01 s");
}

// =============================================================================
// Test 7: ForceWrite flushes to disk immediately
// =============================================================================

void TestForceWriteFlushes()
{
   std::cout << "\n=== Test: ForceWrite Flushes to Disk ===\n";

   BP2Params params;

   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = 50.0e3;
   mesh_params.Lz = 50.0e3;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = 4;
   mesh_params.nz = 8;
   auto mesh = BP2MeshGenerator::Create(mesh_params);

   AntiplaneDomainOperator<Mesh> domain(
      *mesh, 1, params.mu(), params.Vp, params.Wf);

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

   Vector fault_depths;
   domain.GetFaultDepths(fault_depths);

   std::vector<real_t> probe_depths = {0.0, -12000.0};
   std::string prefix = "test_forcewrite_flush";
   BenchmarkOutput<Mesh> bench_out(prefix, params, probe_depths, fault_depths);

   // ForceWrite should flush automatically — no explicit Flush() needed
   bench_out.ForceWrite(0.0, state, fault_op, seas_op.GetTraction());

   // Check files are non-empty immediately (data on disk, not just buffered)
   TEST_ASSERT(FileExistsAndNonEmpty("test_forcewrite_flush_z0km.txt"),
               "ForceWrite flushes z=0km file to disk");
   TEST_ASSERT(FileExistsAndNonEmpty("test_forcewrite_flush_z12km.txt"),
               "ForceWrite flushes z=12km file to disk");

   // Verify data line is present (not just headers)
   {
      std::ifstream f("test_forcewrite_flush_z0km.txt");
      int line_count = 0;
      std::string line;
      while (std::getline(f, line))
      {
         if (!line.empty()) { line_count++; }
      }
      // Should have: description header + column header + 1 data line = 3 lines
      TEST_ASSERT(line_count >= 3,
                  "ForceWrite wrote header + data (got " +
                  std::to_string(line_count) + " lines)");
   }

   bench_out.Close();

   std::remove("test_forcewrite_flush_z0km.txt");
   std::remove("test_forcewrite_flush_z12km.txt");
}

// =============================================================================
// Test 8: Write respects adaptive output interval
// =============================================================================

void TestWriteAdaptiveInterval()
{
   std::cout << "\n=== Test: Write Adaptive Output Interval ===\n";

   BP2Params params;

   BP2MeshGenerator::Parameters mesh_params;
   mesh_params.Lx = 50.0e3;
   mesh_params.Lz = 50.0e3;
   mesh_params.Wf = params.Wf;
   mesh_params.nx = 4;
   mesh_params.nz = 8;
   auto mesh = BP2MeshGenerator::Create(mesh_params);

   AntiplaneDomainOperator<Mesh> domain(
      *mesh, 1, params.mu(), params.Vp, params.Wf);

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

   Vector fault_depths;
   domain.GetFaultDepths(fault_depths);

   std::vector<real_t> probe_depths = {0.0};
   std::string prefix = "test_adaptive_interval";
   BenchmarkOutput<Mesh> bench_out(prefix, params, probe_depths, fault_depths);

   // ForceWrite at t=0
   bench_out.ForceWrite(0.0, state, fault_op, seas_op.GetTraction());

   // Interseismic: V_max = 1e-9, interval = 0.01 yr ~ 315576 s
   // Write at t=1000s should be skipped (too soon)
   bool wrote1 = bench_out.Write(1000.0, state, fault_op, seas_op.GetTraction());
   TEST_ASSERT(!wrote1, "Write skipped at t=1000s (interseismic, need 0.01yr)");

   // Write at t=0.01yr should succeed
   real_t t_01yr = 0.01 * BP2Params::seconds_per_year;
   bool wrote2 = bench_out.Write(t_01yr, state, fault_op, seas_op.GetTraction());
   TEST_ASSERT(wrote2, "Write accepted at t=0.01yr");

   // Immediately after should be skipped again
   bool wrote3 = bench_out.Write(t_01yr + 100.0, state, fault_op,
                                  seas_op.GetTraction());
   TEST_ASSERT(!wrote3, "Write skipped shortly after previous write");

   // Count data lines in the file
   bench_out.Close();
   {
      std::ifstream f("test_adaptive_interval_z0km.txt");
      int data_lines = 0;
      std::string line;
      while (std::getline(f, line))
      {
         if (!line.empty() && line[0] != '#') { data_lines++; }
      }
      // Should have exactly 2 data lines: t=0 (ForceWrite) + t=0.01yr (Write)
      TEST_ASSERT(data_lines == 2,
                  "File has 2 data lines (got " +
                  std::to_string(data_lines) + ")");
   }

   std::remove("test_adaptive_interval_z0km.txt");
}

// =============================================================================
// Test 9: ProbeOutput flush guarantees on-disk data
// =============================================================================

void TestProbeOutputFlush()
{
   std::cout << "\n=== Test: ProbeOutput Flush ===\n";

   std::string filename = "test_probe_flush.txt";

   {
      std::vector<std::string> columns = {"time(s)", "value"};
      ProbeOutput probe(filename, columns);

      probe.WriteStep({0.0, 1.0});

      // Without flush, file may be empty on disk due to buffering.
      // After flush, data must be on disk.
      probe.Flush();

      // Verify file is non-empty while still open
      {
         std::ifstream f(filename);
         f.seekg(0, std::ios::end);
         TEST_ASSERT(f.tellg() > 0,
                     "ProbeOutput: file non-empty after Flush()");
      }

      // Write more data without flushing
      probe.WriteStep({1.0, 2.0});
      probe.WriteStep({2.0, 3.0});
      probe.Flush();

      // Verify all 3 data lines are on disk
      {
         std::ifstream f(filename);
         int data_lines = 0;
         std::string line;
         while (std::getline(f, line))
         {
            if (!line.empty() && line[0] != '#') { data_lines++; }
         }
         TEST_ASSERT(data_lines == 3,
                     "ProbeOutput: 3 data lines after second Flush()");
      }

      probe.Close();
   }

   std::remove(filename.c_str());
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char *argv[])
{
   std::cout << "================================================\n";
   std::cout << "SEAS Phase 5 Tests: I/O and Output\n";
   std::cout << "================================================\n";

   TestProbeOutput();
   TestProbeInterpolator();
   TestBenchmarkOutputFilenames();
   TestBenchmarkOutputIntegration();
   TestParaViewOutput();
   TestParaViewCombinedOutput();
   TestParaViewOutputInterval();
   TestForceWriteFlushes();
   TestWriteAdaptiveInterval();
   TestProbeOutputFlush();

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
