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
//   - Probe2DInterpolator (nearest-DOF fallback + exact face interpolation)
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
// Test: Probe2DInterpolator exact p=1 face interpolation
// ============================================================================

void TestProbe2DInterpolator_ExactFaceInterpolationP1()
{
   std::cout << "\n=== Test: Probe2DInterpolator_ExactFaceInterpolationP1 ===\n";

   // Single p=1 triangular fault face with vertices
   // (0,0), (1000,0), (0,1000) in (x2, x3) coordinates.
   Vector x2(3), x3(3);
   x2(0) = 0.0;    x3(0) = 0.0;
   x2(1) = 1e3;    x3(1) = 0.0;
   x2(2) = 0.0;    x3(2) = 1e3;

   std::vector<Probe2DInterpolator::Station> stations = {
      {"inside", 250.0, 250.0},
   };

   Probe2DInterpolator interp(x2, x3, stations, 3);

   TEST_ASSERT(interp.HasExactMatch(0),
               "p=1 triangular face supports exact station interpolation");
   TEST_NEAR(interp.GetMatchDistance(0), 0.0, 1e-12,
             "Exact face interpolation reports zero match distance");

   Vector field(3);
   field(0) = 10.0;
   field(1) = 20.0;
   field(2) = 40.0;

   // Barycentric weights at (250,250) are [0.5, 0.25, 0.25].
   real_t expected = 0.5 * 10.0 + 0.25 * 20.0 + 0.25 * 40.0;
   TEST_NEAR(interp.EvaluateScalar(field, 0), expected, 1e-12,
             "Exact face interpolation reproduces barycentric value");
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
// Test: Multi-DOF improves station matching (Phase 5)
//
// With more DOFs per face (p>=2), station matching should find closer DOFs.
// Simulate: sparse grid (face centers only) vs dense grid (6 DOFs per face).
// ============================================================================

void TestProbe2DInterpolator_MultiDOFCloserMatch()
{
   std::cout << "\n=== Test: Probe2DInterpolator_MultiDOFCloserMatch ===\n";

   // Case 1: 4 faces, 1 DOF per face (face centers) — coarse grid
   // Face centers at (0, 5e3), (0, 15e3), (20e3, 5e3), (20e3, 15e3)
   const int nfaces = 4;
   Vector x2_coarse(nfaces), x3_coarse(nfaces);
   x2_coarse(0) = 0.0;    x3_coarse(0) = 5e3;
   x2_coarse(1) = 0.0;    x3_coarse(1) = 15e3;
   x2_coarse(2) = 20e3;   x3_coarse(2) = 5e3;
   x2_coarse(3) = 20e3;   x3_coarse(3) = 15e3;

   // Case 2: 4 faces, 6 DOFs per face (p=2) — dense grid
   // Each face gets 6 DOFs spread within it
   const int nbf = 6;
   const int ndofs = nfaces * nbf;
   Vector x2_fine(ndofs), x3_fine(ndofs);
   // Face 0 center (0, 5e3): spread DOFs around center
   real_t offsets_x2[] = {-2e3, 0.0, 2e3, -1e3, 1e3, 0.0};
   real_t offsets_x3[] = {-2e3, -2e3, -2e3, 0.0, 0.0, 2e3};
   real_t face_cx2[] = {0.0, 0.0, 20e3, 20e3};
   real_t face_cx3[] = {5e3, 15e3, 5e3, 15e3};
   for (int f = 0; f < nfaces; f++)
   {
      for (int k = 0; k < nbf; k++)
      {
         int idx = f * nbf + k;
         x2_fine(idx) = face_cx2[f] + offsets_x2[k];
         x3_fine(idx) = face_cx3[f] + offsets_x3[k];
      }
   }

   // Station at (1e3, 6e3) — slightly offset from face 0 center
   std::vector<Probe2DInterpolator::Station> stations = {
      {"test_nearby", 1e3, 6e3},
   };

   Probe2DInterpolator interp_coarse(x2_coarse, x3_coarse, stations);
   Probe2DInterpolator interp_fine(x2_fine, x3_fine, stations);

   real_t dist_coarse = interp_coarse.GetMatchDistance(0);
   real_t dist_fine = interp_fine.GetMatchDistance(0);

   std::cout << "  Coarse (nbf=1) match distance: " << dist_coarse << " m\n";
   std::cout << "  Fine   (nbf=6) match distance: " << dist_fine << " m\n";

   // Fine grid should match closer (or equal)
   TEST_ASSERT(dist_fine <= dist_coarse + 1.0,
               "Multi-DOF provides equal or closer station match");

   // Both should find a valid DOF
   TEST_ASSERT(interp_coarse.GetNearestDOF(0) >= 0,
               "Coarse: valid nearest DOF");
   TEST_ASSERT(interp_fine.GetNearestDOF(0) >= 0,
               "Fine: valid nearest DOF");
}

// ============================================================================
// Test: BP5BenchmarkOutput works with multi-DOF expanded vectors (Phase 5)
//
// Verifies that WriteFromGlobalData correctly indexes multi-DOF data
// (6 DOFs per face at p=2). The output pipeline uses per-DOF indexing
// from Probe2DInterpolator, so it should work generically.
// ============================================================================

void TestBP5BenchmarkOutput_MultiDOFWrite()
{
   std::cout << "\n=== Test: BP5BenchmarkOutput_MultiDOFWrite ===\n";

   // 2 faces × 6 DOFs = 12 total DOFs
   const int nfaces = 2;
   const int nbf = 6;
   const int ndofs = nfaces * nbf;

   // DOF coordinates (spread across two faces)
   Vector x2(ndofs), x3(ndofs);
   for (int f = 0; f < nfaces; f++)
   {
      for (int k = 0; k < nbf; k++)
      {
         int idx = f * nbf + k;
         x2(idx) = f * 10e3 + k * 1e3;  // Spread in x2
         x3(idx) = 10e3 + k * 0.5e3;    // Spread in x3
      }
   }

   // Station near DOF index 3 (face 0, DOF 3)
   std::vector<Probe2DInterpolator::Station> stations = {
      {"test_mdof", x2(3), x3(3)},
   };

   BP5Params params;
   std::string prefix = "test_bp5_mdof";

   Vector tau_pre_dip(ndofs), tau_pre_strike(ndofs);
   tau_pre_dip = 0.0;
   tau_pre_strike = 0.0;
   // Give DOF 3 a known tau_pre
   tau_pre_dip(3) = -5e6;
   tau_pre_strike(3) = -13.27e6;

   {
      BP5BenchmarkOutput<Mesh> out(prefix, params, stations, x2, x3);
      out.SetTauPre(tau_pre_dip, tau_pre_strike);

      // Create per-DOF data vectors
      Vector slip_dip(ndofs), slip_strike(ndofs), theta(ndofs);
      Vector V_dip(ndofs), V_strike(ndofs);
      Vector trac_dip(ndofs), trac_strike(ndofs);

      // Set DOF 3 with known values, others with different values
      for (int i = 0; i < ndofs; i++)
      {
         slip_dip(i) = -0.01 * (i + 1);
         slip_strike(i) = -0.02 * (i + 1);
         theta(i) = 100.0 + i;
         V_dip(i) = 1e-9 * (i + 1);
         V_strike(i) = 2e-9 * (i + 1);
         trac_dip(i) = -1e5 * (i + 1);
         trac_strike(i) = -2e5 * (i + 1);
      }

      out.WriteFromGlobalData(1.0, slip_dip, slip_strike, theta,
                              V_dip, V_strike, trac_dip, trac_strike);
      out.Flush();
   }

   // Read output and verify it picked up DOF 3's values
   std::string filename = prefix + "_test_mdof.txt";
   std::ifstream file(filename);
   TEST_ASSERT(file.is_open(), "Multi-DOF output file created");

   if (file.is_open())
   {
      std::string line;
      std::getline(file, line); // header 1
      std::getline(file, line); // header 2
      std::getline(file, line); // data

      std::istringstream iss(line);
      double t, s_strike, s_dip, v_strike, v_dip, tau_s, tau_d, state_val;
      iss >> t >> s_strike >> s_dip >> v_strike >> v_dip
          >> tau_s >> tau_d >> state_val;

      // DOF 3 (i=3): slip_strike negated = 0.02*4 = 0.08
      TEST_NEAR(s_strike, 0.02 * 4, 1e-6,
                "Multi-DOF: slip_strike from DOF 3");
      // DOF 3 (i=3): slip_dip negated = 0.01*4 = 0.04
      TEST_NEAR(s_dip, 0.01 * 4, 1e-6,
                "Multi-DOF: slip_dip from DOF 3");
      // DOF 3: V_strike = abs(2e-9*4) = 8e-9, output = log10(8e-9) ≈ -8.097
      TEST_NEAR(v_strike, std::log10(2e-9 * 4), 0.01,
                "Multi-DOF: log10(V_strike) from DOF 3");
      // DOF 3: theta = 103, output = log10(103) ≈ 2.013
      TEST_NEAR(state_val, std::log10(103.0), 0.01,
                "Multi-DOF: log10(theta) from DOF 3");
      // tau_strike = (tau_pre_strike + trac_strike) / 1e6
      // = (13.27e6 + 2e5*4) / 1e6 = (13.27e6 + 8e5) / 1e6 = 14.07
      TEST_NEAR(tau_s, 14.07, 0.01,
                "Multi-DOF: tau_strike = tau_pre + elastic");

      file.close();
   }

   std::remove(filename.c_str());
}

// ============================================================================
// Test: BP5BenchmarkOutput exact p=1 face interpolation at off-node station
// ============================================================================

void TestBP5BenchmarkOutput_ExactFaceInterpolationP1()
{
   std::cout << "\n=== Test: BP5BenchmarkOutput_ExactFaceInterpolationP1 ===\n";

   Vector x2(3), x3(3);
   x2(0) = 0.0;    x3(0) = 0.0;
   x2(1) = 1e3;    x3(1) = 0.0;
   x2(2) = 0.0;    x3(2) = 1e3;

   std::vector<Probe2DInterpolator::Station> stations = {
      {"test_exact_p1", 250.0, 250.0},
   };

   BP5Params params;
   std::string prefix = "test_bp5_exact_p1";

   Vector tau_pre_dip(3), tau_pre_strike(3);
   tau_pre_dip = 0.0;
   tau_pre_strike = 0.0;

   {
      BP5BenchmarkOutput<Mesh> out(prefix, params, stations, x2, x3, 3);
      out.SetTauPre(tau_pre_dip, tau_pre_strike);

      Vector slip_dip(3), slip_strike(3), theta(3);
      Vector V_dip(3), V_strike(3);
      Vector trac_dip(3), trac_strike(3);

      // Internal values are negated before SCEC output.
      slip_dip(0) = -1.0;   slip_strike(0) = -10.0;
      slip_dip(1) = -2.0;   slip_strike(1) = -20.0;
      slip_dip(2) = -4.0;   slip_strike(2) = -40.0;

      V_dip(0) = 1e-6;      V_strike(0) = 1e-5;
      V_dip(1) = 2e-6;      V_strike(1) = 2e-5;
      V_dip(2) = 4e-6;      V_strike(2) = 4e-5;

      trac_dip = 0.0;
      trac_strike = 0.0;

      theta(0) = 10.0;
      theta(1) = 20.0;
      theta(2) = 40.0;

      out.WriteFromGlobalData(1.0, slip_dip, slip_strike, theta,
                              V_dip, V_strike, trac_dip, trac_strike);
      out.Flush();
   }

   std::string filename = prefix + "_test_exact_p1.txt";
   std::ifstream file(filename);
   TEST_ASSERT(file.is_open(), "Exact p=1 BP5 output file created");

   if (file.is_open())
   {
      std::string line;
      std::getline(file, line);
      std::getline(file, line);
      std::getline(file, line);

      std::istringstream iss(line);
      double t, slip_s, slip_d, logV_s, logV_d, tau_s, tau_d, log_theta;
      iss >> t >> slip_s >> slip_d >> logV_s >> logV_d >> tau_s >> tau_d >> log_theta;

      // Barycentric weights = [0.5, 0.25, 0.25] at (250,250).
      TEST_NEAR(slip_s, 20.0, 1e-12,
                "strike slip interpolated exactly on p=1 face");
      TEST_NEAR(slip_d, 2.0, 1e-12,
                "dip slip interpolated exactly on p=1 face");
      TEST_NEAR(logV_s, std::log10(2e-5), 1e-12,
                "strike slip-rate interpolated exactly on p=1 face");
      TEST_NEAR(logV_d, std::log10(2e-6), 1e-12,
                "dip slip-rate interpolated exactly on p=1 face");
      TEST_NEAR(log_theta, std::log10(20.0), 1e-12,
                "state variable interpolated exactly on p=1 face");

      file.close();
   }

   std::remove(filename.c_str());
}

// ============================================================================
// Test: PrintDiagnostics produces valid output
// ============================================================================

void TestProbe2DInterpolator_PrintDiagnostics()
{
   std::cout << "\n=== Test: Probe2DInterpolator_PrintDiagnostics ===\n";

   // Single p=1 triangle
   Vector x2(3), x3(3);
   x2(0) = 0.0;    x3(0) = 0.0;
   x2(1) = 1e3;    x3(1) = 0.0;
   x2(2) = 0.0;    x3(2) = 1e3;

   std::vector<Probe2DInterpolator::Station> stations = {
      {"inside", 250.0, 250.0},
      {"outside", 2e3, 2e3},
   };

   Probe2DInterpolator interp(x2, x3, stations, 3);

   // Station inside should have exact match
   TEST_ASSERT(interp.HasExactMatch(0),
               "Station inside triangle has exact match");
   TEST_ASSERT(!interp.HasExactMatch(1),
               "Station outside triangle falls back to nearest-DOF");

   // PrintDiagnostics should not crash
   std::ostringstream oss;
   interp.PrintDiagnostics(stations, x2, x3, oss);
   std::string diag = oss.str();
   TEST_ASSERT(diag.find("exact_match = true") != std::string::npos,
               "Diagnostics show exact match for inside station");
   TEST_ASSERT(diag.find("exact_match = false") != std::string::npos,
               "Diagnostics show no exact match for outside station");
   TEST_ASSERT(diag.find("DEGENERATE") == std::string::npos,
               "No degenerate triangle detected");

   std::cout << diag;
}

// ============================================================================
// Test: CountDuplicateFaces detects shared-face duplication
// ============================================================================

void TestProbe2DInterpolator_CountDuplicateFaces()
{
   std::cout << "\n=== Test: Probe2DInterpolator_CountDuplicateFaces ===\n";

   // Simulate 3 faces gathered from 2 ranks:
   //   Rank 0: face A (unique) + face B (shared)
   //   Rank 1: face B (shared, same coords) + face C (unique)
   // After gather: [faceA, faceB_rank0, faceB_rank1, faceC] = 4 face blocks
   // But faceB appears twice.

   const int nbf = 3;
   const int nfaces = 4;
   Vector x2(nfaces * nbf), x3(nfaces * nbf);

   // Face A: (0, 0), (1000, 0), (0, 1000) — unique
   x2(0) = 0.0;    x3(0) = 0.0;
   x2(1) = 1e3;    x3(1) = 0.0;
   x2(2) = 0.0;    x3(2) = 1e3;

   // Face B (rank 0 copy): (2000, 0), (3000, 0), (2000, 1000)
   x2(3) = 2e3;    x3(3) = 0.0;
   x2(4) = 3e3;    x3(4) = 0.0;
   x2(5) = 2e3;    x3(5) = 1e3;

   // Face B (rank 1 copy): same coordinates
   x2(6) = 2e3;    x3(6) = 0.0;
   x2(7) = 3e3;    x3(7) = 0.0;
   x2(8) = 2e3;    x3(8) = 1e3;

   // Face C: (4000, 0), (5000, 0), (4000, 1000) — unique
   x2(9)  = 4e3;   x3(9)  = 0.0;
   x2(10) = 5e3;   x3(10) = 0.0;
   x2(11) = 4e3;   x3(11) = 1e3;

   std::ostringstream oss;
   int dups = Probe2DInterpolator::CountDuplicateFaces(x2, x3, nbf, 1.0, &oss);
   TEST_ASSERT(dups == 1,
               "Detected exactly 1 duplicate face pair");

   std::cout << oss.str();

   // No duplicates when all faces are unique
   x2(6) = 6e3;  // Move rank 1's "face B" to a different location
   int dups2 = Probe2DInterpolator::CountDuplicateFaces(x2, x3, nbf, 1.0);
   TEST_ASSERT(dups2 == 0,
               "No duplicates when all faces are unique");
}

// ============================================================================
// Test: Exact interpolation with duplicated shared faces
// ============================================================================

void TestProbe2DInterpolator_ExactWithDuplicates()
{
   std::cout << "\n=== Test: Probe2DInterpolator_ExactWithDuplicates ===\n";

   // Simulate what happens when a shared face is duplicated in the global
   // array. Station is inside the duplicated face. Both copies should
   // produce the same interpolated value if field data is consistent.

   const int nbf = 3;
   // 3 faces: faceA, faceB_rank0, faceB_rank1 (duplicate)
   Vector x2(9), x3(9);

   // Face A
   x2(0) = 0.0;    x3(0) = 0.0;
   x2(1) = 1e3;    x3(1) = 0.0;
   x2(2) = 0.0;    x3(2) = 1e3;

   // Face B (rank 0)
   x2(3) = 2e3;    x3(3) = 0.0;
   x2(4) = 3e3;    x3(4) = 0.0;
   x2(5) = 2e3;    x3(5) = 1e3;

   // Face B (rank 1) — same coords
   x2(6) = 2e3;    x3(6) = 0.0;
   x2(7) = 3e3;    x3(7) = 0.0;
   x2(8) = 2e3;    x3(8) = 1e3;

   // Station inside face B
   std::vector<Probe2DInterpolator::Station> stations = {
      {"in_faceB", 2500.0, 250.0},
   };

   Probe2DInterpolator interp(x2, x3, stations, nbf);

   TEST_ASSERT(interp.HasExactMatch(0),
               "Station inside duplicated face has exact match");

   // The interpolator should match the FIRST occurrence (rank 0's copy)
   TEST_ASSERT(interp.GetFaceStart(0) == 3,
               "Matched face_start is rank 0's copy (index 3)");

   // Evaluate with consistent field data across both copies
   Vector field(9);
   field(0) = 100.0;  field(1) = 200.0;  field(2) = 300.0;  // face A
   field(3) = 10.0;   field(4) = 20.0;   field(5) = 40.0;   // face B rank 0
   field(6) = 10.0;   field(7) = 20.0;   field(8) = 40.0;   // face B rank 1 (same)

   real_t val = interp.EvaluateScalar(field, 0);

   // Barycentric at (2500, 250) in triangle (2000,0)-(3000,0)-(2000,1000):
   // r = (500*1000 - 250*0) / (1000*1000 - 0*0) = 0.5
   // s = (0*250 - (-1000)*500) / (1000*1000) = ... let me compute properly
   // Using ComputeReferenceIP logic:
   // ax = 3000-2000 = 1000, az = 0-0 = 0
   // bx = 2000-2000 = 0, bz = 1000-0 = 1000
   // px = 2500-2000 = 500, pz = 250-0 = 250
   // det = 1000*1000 - 0*0 = 1e6
   // r = (500*1000 - 250*0) / 1e6 = 0.5
   // s = (1000*250 - 0*500) / 1e6 = 0.25
   // l0 = 1 - 0.5 - 0.25 = 0.25, l1 = 0.5, l2 = 0.25
   // For p=1 H1_TriangleElement with GaussLobatto: shape = [l0, l1, l2]
   // val = 0.25*10 + 0.5*20 + 0.25*40 = 2.5 + 10 + 10 = 22.5
   TEST_NEAR(val, 22.5, 1e-12,
             "Exact interpolation correct with duplicated face (consistent data)");

   // Now test with INCONSISTENT field data across copies (the dangerous case)
   field(6) = 15.0;   field(7) = 25.0;   field(8) = 45.0;  // rank 1 differs

   real_t val_consistent = interp.EvaluateScalar(field, 0);
   // Should still use rank 0's copy (face_start = 3)
   TEST_NEAR(val_consistent, 22.5, 1e-12,
             "Uses rank 0's copy even when rank 1 differs (first-match wins)");
}

// ============================================================================
// Test: Exact interpolation at face boundary and vertex
// ============================================================================

void TestProbe2DInterpolator_ExactBoundaryAndVertex()
{
   std::cout << "\n=== Test: Probe2DInterpolator_ExactBoundaryAndVertex ===\n";

   Vector x2(3), x3(3);
   x2(0) = 0.0;    x3(0) = 0.0;
   x2(1) = 1e3;    x3(1) = 0.0;
   x2(2) = 0.0;    x3(2) = 1e3;

   Vector field(3);
   field(0) = 10.0;  field(1) = 20.0;  field(2) = 30.0;

   // Station exactly at vertex 0
   {
      std::vector<Probe2DInterpolator::Station> stations = {
         {"at_v0", 0.0, 0.0},
      };
      Probe2DInterpolator interp(x2, x3, stations, 3);
      TEST_ASSERT(interp.HasExactMatch(0),
                  "Station at vertex has exact match");
      real_t val = interp.EvaluateScalar(field, 0);
      TEST_NEAR(val, 10.0, 1e-12,
                "Exact interpolation at vertex 0 gives vertex value");
   }

   // Station at midpoint of edge (v0-v1)
   {
      std::vector<Probe2DInterpolator::Station> stations = {
         {"mid_edge", 500.0, 0.0},
      };
      Probe2DInterpolator interp(x2, x3, stations, 3);
      TEST_ASSERT(interp.HasExactMatch(0),
                  "Station at edge midpoint has exact match");
      real_t val = interp.EvaluateScalar(field, 0);
      // At (500, 0): r=0.5, s=0, l0=0.5, l1=0.5, l2=0
      TEST_NEAR(val, 15.0, 1e-12,
                "Exact interpolation at edge midpoint correct");
   }
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
   TestProbe2DInterpolator_ExactFaceInterpolationP1();
   TestBP5BenchmarkOutput_FileCreation();
   TestBP5BenchmarkOutput_ComponentSwap();
   TestBP5BenchmarkOutput_StressComputation();
   TestBP5BenchmarkOutput_AdaptiveOutput();
   TestBP5BenchmarkOutput_DefaultStations();

   // v45 Phase 5: Multi-DOF output tests
   TestProbe2DInterpolator_MultiDOFCloserMatch();
   TestBP5BenchmarkOutput_MultiDOFWrite();
   TestBP5BenchmarkOutput_ExactFaceInterpolationP1();

   // v56: Diagnostics and duplicate detection tests
   TestProbe2DInterpolator_PrintDiagnostics();
   TestProbe2DInterpolator_CountDuplicateFaces();
   TestProbe2DInterpolator_ExactWithDuplicates();
   TestProbe2DInterpolator_ExactBoundaryAndVertex();

   TEST_PRINT_RESULTS();
   return num_failed;
}
