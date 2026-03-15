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

// Unit tests for BP5 output infrastructure:
//   - Probe2DInterpolator (2D nearest-DOF matching)
//   - BP5BenchmarkOutput (8-column SCEC format output)

#include "mfem.hpp"
#include "test_macros.hpp"
#include "../../io/bp5_benchmark_output.hpp"
#include "../../config/bp5_params.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>

using namespace mfem;
using namespace mfem::seas;

// ============================================================================
// Test: Probe2DInterpolator nearest-DOF matching
// ============================================================================

void TestProbe2DInterpolator_NearestDOF()
{
   std::cout << "\n=== Test: Probe2DInterpolator_NearestDOF ===\n";

   // Create a synthetic 10x8 = 80 DOF grid
   // x2: [-40e3, -30e3, ..., 40e3, 50e3] (10 points, 10 km spacing)
   // x3: [0, 5e3, 10e3, ..., 35e3] (8 points, 5 km spacing)
   const int nx = 10, nz = 8;
   const int n_dofs = nx * nz;
   Vector fault_x2(n_dofs), fault_x3(n_dofs);

   for (int j = 0; j < nz; j++)
   {
      for (int i = 0; i < nx; i++)
      {
         int idx = j * nx + i;
         fault_x2(idx) = (-40.0 + i * 10.0) * 1e3;
         fault_x3(idx) = j * 5.0 * 1e3;
      }
   }

   // Query a few stations
   std::vector<Probe2DInterpolator::Station> stations = {
      {"s1", 0.0, 10e3},     // Should match (0, 10e3)
      {"s2", -24e3, 10e3},   // Should match (-20e3, 10e3)
      {"s3", 50e3, 35e3},    // Should match (50e3, 35e3) - corner DOF
   };

   Probe2DInterpolator interp(fault_x2, fault_x3, stations);

   TEST_ASSERT(interp.NumStations() == 3, "NumStations() == 3");

   // Station s1 at (0, 10e3): nearest DOF at (0, 10e3) = index 2*10+4 = 24
   int dof0 = interp.GetNearestDOF(0);
   TEST_NEAR(fault_x2(dof0), 0.0, 1.0, "s1 matched x2 ≈ 0");
   TEST_NEAR(fault_x3(dof0), 10e3, 1.0, "s1 matched x3 ≈ 10km");
   TEST_NEAR(interp.GetMatchDistance(0), 0.0, 1.0, "s1 exact match distance ≈ 0");

   // Station s2 at (-24e3, 10e3): nearest DOF at (-20e3, 10e3)
   int dof1 = interp.GetNearestDOF(1);
   TEST_NEAR(fault_x2(dof1), -20e3, 1.0, "s2 matched x2 ≈ -20km");
   TEST_NEAR(fault_x3(dof1), 10e3, 1.0, "s2 matched x3 ≈ 10km");
   TEST_NEAR(interp.GetMatchDistance(1), 4e3, 1.0, "s2 match distance ≈ 4km");

   // Station s3 at (50e3, 35e3): exact match at corner
   int dof2 = interp.GetNearestDOF(2);
   TEST_NEAR(fault_x2(dof2), 50e3, 1.0, "s3 matched x2 ≈ 50km");
   TEST_NEAR(fault_x3(dof2), 35e3, 1.0, "s3 matched x3 ≈ 35km");
   TEST_NEAR(interp.GetMatchDistance(2), 0.0, 1.0, "s3 exact match distance ≈ 0");
}

// ============================================================================
// Test: Probe2DInterpolator with empty fault DOFs
// ============================================================================

void TestProbe2DInterpolator_NoFaultDOFs()
{
   std::cout << "\n=== Test: Probe2DInterpolator_NoFaultDOFs ===\n";

   Vector empty_x2(0), empty_x3(0);
   std::vector<Probe2DInterpolator::Station> stations = {
      {"s1", 0.0, 10e3},
   };

   Probe2DInterpolator interp(empty_x2, empty_x3, stations);
   TEST_ASSERT(interp.GetNearestDOF(0) == -1,
               "Empty fault DOFs → nearest DOF = -1");
   TEST_ASSERT(interp.GetMatchDistance(0) > 1e30,
               "Empty fault DOFs → infinite match distance");
}

// ============================================================================
// Test: BP5BenchmarkOutput file creation and headers
// ============================================================================

void TestBP5BenchmarkOutput_FileCreation()
{
   std::cout << "\n=== Test: BP5BenchmarkOutput_FileCreation ===\n";

   // Create a small synthetic fault with 4 DOFs
   Vector x2(4), x3(4);
   x2(0) = 0.0;   x3(0) = 5e3;
   x2(1) = 0.0;   x3(1) = 10e3;
   x2(2) = 0.0;   x3(2) = 15e3;
   x2(3) = 16e3;  x3(3) = 10e3;

   std::vector<Probe2DInterpolator::Station> stations = {
      {"fltst_strk+00dp+05", 0.0, 5e3},
      {"fltst_strk+00dp+10", 0.0, 10e3},
      {"fltst_strk+16dp+10", 16e3, 10e3},
   };

   BP5Params params;
   std::string prefix = "test_bp5out";

   {
      BP5BenchmarkOutput<Mesh> out(prefix, params, stations, x2, x3);

      TEST_ASSERT(out.NumProbes() == 3, "NumProbes() == 3");
   }

   // Verify files exist and have correct headers
   for (const auto &st : stations)
   {
      std::string filename = prefix + "_" + st.name + ".txt";
      std::ifstream file(filename);
      TEST_ASSERT(file.is_open(), "File created: " + filename);

      if (file.is_open())
      {
         std::string line1, line2;
         std::getline(file, line1);
         std::getline(file, line2);

         TEST_ASSERT(line1.find("BP5-QD") != std::string::npos,
                     "Header contains BP5-QD: " + st.name);
         TEST_ASSERT(line2.find("slip_strike") != std::string::npos,
                     "Header contains slip_strike: " + st.name);
         TEST_ASSERT(line2.find("tau_dip") != std::string::npos,
                     "Header contains tau_dip: " + st.name);
         file.close();
      }

      // Clean up
      std::remove(filename.c_str());
   }
}

// ============================================================================
// Test: BP5BenchmarkOutput component swap (dip/strike → strike/dip)
// ============================================================================

void TestBP5BenchmarkOutput_ComponentSwap()
{
   std::cout << "\n=== Test: BP5BenchmarkOutput_ComponentSwap ===\n";

   // 1 DOF, 1 station at exact match
   Vector x2(1), x3(1);
   x2(0) = 0.0; x3(0) = 10e3;

   std::vector<Probe2DInterpolator::Station> stations = {
      {"test_station", 0.0, 10e3},
   };

   BP5Params params;
   std::string prefix = "test_bp5swap";

   // Set known tau_pre for stress computation
   Vector tau_pre_dip(1), tau_pre_strike(1);
   tau_pre_dip(0) = 0.0;
   tau_pre_strike(0) = 0.0;

   {
      BP5BenchmarkOutput<Mesh> out(prefix, params, stations, x2, x3);
      out.SetTauPre(tau_pre_dip, tau_pre_strike);

      // Write one row of known data using WriteFromGlobalData
      // slip: dip=1.0, strike=2.0
      Vector slip_dip(1), slip_strike(1), theta(1);
      Vector V_dip(1), V_strike(1);
      Vector trac_dip(1), trac_strike(1);

      // Internal convention: negative = right-lateral motion.
      // WriteFromGlobalData negates for SCEC output (positive = right-lateral).
      slip_dip(0) = -1.0;
      slip_strike(0) = -2.0;
      theta(0) = 100.0;      // 100 seconds
      V_dip(0) = 1e-6;
      V_strike(0) = 1e-5;
      trac_dip(0) = 0.0;
      trac_strike(0) = 0.0;

      out.WriteFromGlobalData(1.0, slip_dip, slip_strike, theta,
                              V_dip, V_strike, trac_dip, trac_strike);
      out.Flush();
   }

   // Read output and verify column order: strike first, dip second
   std::string filename = prefix + "_test_station.txt";
   std::ifstream file(filename);
   TEST_ASSERT(file.is_open(), "Output file opened for swap test");

   if (file.is_open())
   {
      std::string line;
      // Skip 2 header lines
      std::getline(file, line);
      std::getline(file, line);
      // Read data line
      std::getline(file, line);

      std::istringstream iss(line);
      double t, col2, col3;
      iss >> t >> col2 >> col3;

      // col2 = slip_strike = 2.0, col3 = slip_dip = 1.0 (SWAPPED)
      TEST_NEAR(col2, 2.0, 1e-6, "Column 2 = slip_strike = 2.0 (swapped)");
      TEST_NEAR(col3, 1.0, 1e-6, "Column 3 = slip_dip = 1.0 (swapped)");
      file.close();
   }

   std::remove(filename.c_str());
}

// ============================================================================
// Test: BP5BenchmarkOutput stress computation
// ============================================================================

void TestBP5BenchmarkOutput_StressComputation()
{
   std::cout << "\n=== Test: BP5BenchmarkOutput_StressComputation ===\n";

   Vector x2(1), x3(1);
   x2(0) = 0.0; x3(0) = 10e3;

   std::vector<Probe2DInterpolator::Station> stations = {
      {"test_stress", 0.0, 10e3},
   };

   BP5Params params;
   std::string prefix = "test_bp5stress";

   // Set known tau_pre: dip=-5e6 Pa, strike=-10e6 Pa
   // Internal convention: negative = right-lateral shear.
   // WriteFromGlobalData negates for SCEC output (positive = right-lateral).
   Vector tau_pre_dip(1), tau_pre_strike(1);
   tau_pre_dip(0) = -5e6;
   tau_pre_strike(0) = -10e6;

   {
      BP5BenchmarkOutput<Mesh> out(prefix, params, stations, x2, x3);
      out.SetTauPre(tau_pre_dip, tau_pre_strike);

      Vector slip_dip(1), slip_strike(1), theta(1);
      Vector V_dip(1), V_strike(1);
      Vector trac_dip(1), trac_strike(1);

      slip_dip(0) = 0.0; slip_strike(0) = 0.0;
      theta(0) = 100.0;
      V_dip(0) = 1e-9; V_strike(0) = 1e-9;
      // Elastic traction: dip=-1e6, strike=-2e6
      // Internal convention: negative = right-lateral shear.
      trac_dip(0) = -1e6;
      trac_strike(0) = -2e6;

      out.WriteFromGlobalData(1.0, slip_dip, slip_strike, theta,
                              V_dip, V_strike, trac_dip, trac_strike);
      out.Flush();
   }

   std::string filename = prefix + "_test_stress.txt";
   std::ifstream file(filename);
   TEST_ASSERT(file.is_open(), "Stress output file opened");

   if (file.is_open())
   {
      std::string line;
      std::getline(file, line); // header 1
      std::getline(file, line); // header 2
      std::getline(file, line); // data

      std::istringstream iss(line);
      double t, s1, s2, v1, v2, tau_s, tau_d, state;
      iss >> t >> s1 >> s2 >> v1 >> v2 >> tau_s >> tau_d >> state;

      // tau_strike = (tau_pre_strike + trac_strike) / 1e6 = (10e6 + 2e6) / 1e6 = 12
      // tau_dip = (tau_pre_dip + trac_dip) / 1e6 = (5e6 + 1e6) / 1e6 = 6
      TEST_NEAR(tau_s, 12.0, 1e-6,
                "tau_strike = (10 + 2) MPa = 12 MPa");
      TEST_NEAR(tau_d, 6.0, 1e-6,
                "tau_dip = (5 + 1) MPa = 6 MPa");
      file.close();
   }

   std::remove(filename.c_str());
}

// ============================================================================
// Test: BP5BenchmarkOutput adaptive output interval
// ============================================================================

void TestBP5BenchmarkOutput_AdaptiveOutput()
{
   std::cout << "\n=== Test: BP5BenchmarkOutput_AdaptiveOutput ===\n";

   // Coseismic: V > 1e-3 → 0.1 s (SCEC spec)
   real_t dt1 = BP5BenchmarkOutput<Mesh>::OutputInterval(1.0);
   TEST_NEAR(dt1, 0.1, 1e-10, "V=1.0 → dt=0.1s (coseismic)");

   real_t dt1b = BP5BenchmarkOutput<Mesh>::OutputInterval(5e-3);
   TEST_NEAR(dt1b, 0.1, 1e-10, "V=5e-3 → dt=0.1s (coseismic)");

   // Nucleation: 1e-6 < V < 1e-3 → 0.1 s
   real_t dt2 = BP5BenchmarkOutput<Mesh>::OutputInterval(1e-4);
   TEST_NEAR(dt2, 0.1, 1e-10, "V=1e-4 → dt=0.1s (nucleation)");

   // Interseismic: V < 1e-6 → 0.1 yr (SCEC spec)
   real_t dt3 = BP5BenchmarkOutput<Mesh>::OutputInterval(1e-9);
   real_t expected = 0.1 * BP5Params::seconds_per_year;
   TEST_NEAR(dt3, expected, 1.0, "V=1e-9 → dt=0.1yr (interseismic)");
}

// ============================================================================
// Test: BP5BenchmarkOutput default stations
// ============================================================================

void TestBP5BenchmarkOutput_DefaultStations()
{
   std::cout << "\n=== Test: BP5BenchmarkOutput_DefaultStations ===\n";

   auto stations = BP5BenchmarkOutput<Mesh>::DefaultStations();

   TEST_ASSERT(static_cast<int>(stations.size()) == 10,
               "DefaultStations() returns 10 stations");

   // Verify a few known stations
   // fltst_strk+00dp+10 at (0, 10 km)
   bool found_00_10 = false;
   for (const auto &st : stations)
   {
      if (st.name == "fltst_strk+00dp+10")
      {
         TEST_NEAR(st.x2, 0.0, 1.0, "strk+00dp+10: x2 = 0");
         TEST_NEAR(st.x3, 10e3, 1.0, "strk+00dp+10: x3 = 10 km");
         found_00_10 = true;
      }
   }
   TEST_ASSERT(found_00_10, "Station fltst_strk+00dp+10 found");

   // fltst_strk-24dp+10 at (-24 km, 10 km)
   bool found_m24_10 = false;
   for (const auto &st : stations)
   {
      if (st.name == "fltst_strk-24dp+10")
      {
         TEST_NEAR(st.x2, -24e3, 1.0, "strk-24dp+10: x2 = -24 km");
         TEST_NEAR(st.x3, 10e3, 1.0, "strk-24dp+10: x3 = 10 km");
         found_m24_10 = true;
      }
   }
   TEST_ASSERT(found_m24_10, "Station fltst_strk-24dp+10 found");

   // fltst_strk+00dp+00 at (0, 0) — free surface
   bool found_00_00 = false;
   for (const auto &st : stations)
   {
      if (st.name == "fltst_strk+00dp+00")
      {
         TEST_NEAR(st.x2, 0.0, 1.0, "strk+00dp+00: x2 = 0");
         TEST_NEAR(st.x3, 0.0, 1.0, "strk+00dp+00: x3 = 0");
         found_00_00 = true;
      }
   }
   TEST_ASSERT(found_00_00, "Station fltst_strk+00dp+00 found");
}

// ============================================================================
// Main
// ============================================================================

int main()
{
   std::cout << "BP5 Output Infrastructure Unit Tests\n";
   std::cout << "=====================================\n";

   TestProbe2DInterpolator_NearestDOF();
   TestProbe2DInterpolator_NoFaultDOFs();
   TestBP5BenchmarkOutput_FileCreation();
   TestBP5BenchmarkOutput_ComponentSwap();
   TestBP5BenchmarkOutput_StressComputation();
   TestBP5BenchmarkOutput_AdaptiveOutput();
   TestBP5BenchmarkOutput_DefaultStations();

   TEST_PRINT_RESULTS();
   return num_failed;
}
