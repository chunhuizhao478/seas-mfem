// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_tpv31_station_writer.cpp — coverage for the TPV31 on-fault
// station-trace writer (dynamic/tpv31_stations.hpp).
//
//   1. DefaultStations_TPV31 returns the spec's 30-station 3×10 grid
//      with the correct SCEC `faultst###dp###` names and canonical
//      (along-strike, down-dip) coordinates.
//   2. FindNearestDOF_TPV31 selects the closest fault DOF using |z| as
//      the down-dip depth.
//   3. WriteStep emits the 9-column TPV205-mirror layout in the right
//      order (h = strike = comp 2, v = dip = comp 1) with the LSW
//      mu_eff(δ) as the trailing column.

#include "mfem.hpp"

#include "../../dynamic/tpv31_stations.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e \
   << ", tol " << _t << ")\n"; num_failed++; } \
   else { std::cout << "  PASSED: " << m << "\n"; num_passed++; } \
   } while (0)

// ---------------------------------------------------------------------
// S-1  Station list: 30 stations, correct names + canonical coords.
// ---------------------------------------------------------------------
static void S1_station_list()
{
   std::cout << "\n[S-1] DefaultStations_TPV31: 30-station 3x10 grid\n";
   const std::vector<TPV31Station> st = DefaultStations_TPV31();
   TEST_ASSERT(st.size() == 30, "30 on-fault stations (spec Part 5)");
   if (st.size() != 30) { return; }

   // First station: strike 0, depth 0.
   TEST_ASSERT(st[0].name == "faultst000dp000", "st[0] name");
   TEST_NEAR(st[0].along_strike, 0.0, 1e-9, "st[0] along_strike");
   TEST_NEAR(st[0].down_dip,     0.0, 1e-9, "st[0] down_dip");

   // The hypocenter station: strike 0, depth 7.5 km (8th depth → idx 7).
   TEST_ASSERT(st[7].name == "faultst000dp075", "st[7] name (hypocenter)");
   TEST_NEAR(st[7].along_strike, 0.0,    1e-9, "st[7] along_strike");
   TEST_NEAR(st[7].down_dip,     7500.0, 1e-9, "st[7] down_dip = 7500");

   // First station of the 6 km along-strike block (idx 10).
   TEST_ASSERT(st[10].name == "faultst060dp000", "st[10] name (6 km strike)");
   TEST_NEAR(st[10].along_strike, 6000.0, 1e-9, "st[10] along_strike = 6000");

   // faultst060dp024 — the spec's worked example, (6000, 2400).  idx 14.
   TEST_ASSERT(st[14].name == "faultst060dp024", "st[14] name (spec example)");
   TEST_NEAR(st[14].along_strike, 6000.0, 1e-9, "st[14] along_strike = 6000");
   TEST_NEAR(st[14].down_dip,     2400.0, 1e-9, "st[14] down_dip = 2400");

   // Last station: strike 12 km, depth 12 km.
   TEST_ASSERT(st[29].name == "faultst120dp120", "st[29] name (deepest/far)");
   TEST_NEAR(st[29].along_strike, 12000.0, 1e-9, "st[29] along_strike = 12000");
   TEST_NEAR(st[29].down_dip,     12000.0, 1e-9, "st[29] down_dip = 12000");

   // dp002 / dp005 / dp010 zero-padding spot checks.
   TEST_ASSERT(st[1].name == "faultst000dp002", "0.2 km → dp002");
   TEST_ASSERT(st[2].name == "faultst000dp005", "0.5 km → dp005");
   TEST_ASSERT(st[3].name == "faultst000dp010", "1.0 km → dp010");
}

// ---------------------------------------------------------------------
// S-2  Nearest-DOF selection uses |z| as down-dip depth.
// ---------------------------------------------------------------------
static void S2_nearest_dof()
{
   std::cout << "\n[S-2] FindNearestDOF_TPV31 picks closest (x, |z|)\n";
   // Three candidate DOFs on the y=0 fault at z<0 (depth = |z|).
   std::vector<Vector> coords(3);
   for (auto &c : coords) { c.SetSize(3); }
   coords[0](0) = 0.0;    coords[0](1) = 0.0; coords[0](2) =     0.0; // (0, 0)
   coords[1](0) = 0.0;    coords[1](1) = 0.0; coords[1](2) = -7500.0; // (0, 7500)
   coords[2](0) = 6000.0; coords[2](1) = 0.0; coords[2](2) = -2400.0; // (6000, 2400)

   const std::vector<TPV31Station> st = DefaultStations_TPV31();
   // st[7] = faultst000dp075 → (0, 7500) → DOF 1.
   TEST_ASSERT(FindNearestDOF_TPV31(st[7], coords, 3) == 1,
               "hypocenter station maps to DOF at z=-7500");
   // st[14] = faultst060dp024 → (6000, 2400) → DOF 2.
   TEST_ASSERT(FindNearestDOF_TPV31(st[14], coords, 3) == 2,
               "faultst060dp024 maps to DOF at (6000, -2400)");
   // st[0] = faultst000dp000 → (0, 0) → DOF 0.
   TEST_ASSERT(FindNearestDOF_TPV31(st[0], coords, 3) == 0,
               "surface station maps to DOF at z=0");
   TEST_ASSERT(FindNearestDOF_TPV31(st[0], coords, 0) == -1,
               "ndof=0 returns -1");
}

// ---------------------------------------------------------------------
// S-3  WriteStep column order + mu_eff.
// ---------------------------------------------------------------------
static void S3_writestep_columns()
{
   std::cout << "\n[S-3] WriteStep emits 9 columns in TPV205-mirror order\n";

   // One station at (0, 7500); one DOF carrying distinct sentinel values
   // so we can detect column-order swaps.
   std::vector<TPV31Station> stations(1);
   stations[0] = { 0.0, 7500.0, "faultst000dp075" };

   std::vector<Vector> coords(1);
   coords[0].SetSize(3);
   coords[0](0) = 0.0; coords[0](1) = 0.0; coords[0](2) = -7500.0;

   std::vector<DOFData> dof(1);
   DOFData &d = dof[0];
   d.slip2     = 2.0;        // h-slip
   d.V2        = 3.0;        // h-slip-rate
   d.tau2_corr = 4.0e7;      // h-shear-stress (Pa)
   d.slip1     = 0.0;        // v-slip
   d.V1        = 0.0;        // v-slip-rate
   d.tau1_corr = 0.0;        // v-shear-stress
   d.sigma_n_corr = 6.0e7;   // n-stress (Pa, compression-positive)
   d.lsw_mu_s  = 0.580;
   d.lsw_mu_d  = 0.450;
   d.lsw_d_c   = 0.18;
   // δ = sqrt(0² + 2²) = 2.0 > d_c → mu_eff = mu_d = 0.450.

   const std::string dir = ".";
   const std::string prefix = "test_tpv31_sw_tmp";
   const std::string fname = dir + "/" + prefix
                             + "_station_faultst000dp075.dat";
   std::remove(fname.c_str());

   {
      TPV31StationWriter w;
      w.Open(dir, prefix, stations, coords, /*ndof=*/1);
      w.WriteStep(0.5, dof);
      w.Flush();
      w.Close();
   }

   std::ifstream in(fname);
   TEST_ASSERT(in.good(), "station file created");
   if (!in.good()) { return; }

   std::string line, last_data;
   while (std::getline(in, line))
   {
      if (!line.empty() && line[0] != '#') { last_data = line; }
   }
   in.close();

   std::istringstream iss(last_data);
   real_t t, h_slip, h_rate, h_shear, v_slip, v_rate, v_shear, n_stress, mu_eff;
   iss >> t >> h_slip >> h_rate >> h_shear >> v_slip >> v_rate
       >> v_shear >> n_stress >> mu_eff;

   TEST_NEAR(t,        0.5,   1e-9, "col1 t");
   TEST_NEAR(h_slip,   2.0,   1e-9, "col2 h-slip = slip2");
   TEST_NEAR(h_rate,   3.0,   1e-9, "col3 h-slip-rate = V2");
   TEST_NEAR(h_shear,  4.0e7, 1.0,  "col4 h-shear-stress = tau2_corr (Pa)");
   TEST_NEAR(v_slip,   0.0,   1e-9, "col5 v-slip = slip1");
   TEST_NEAR(v_rate,   0.0,   1e-9, "col6 v-slip-rate = V1");
   TEST_NEAR(v_shear,  0.0,   1e-9, "col7 v-shear-stress = tau1_corr");
   TEST_NEAR(n_stress, 6.0e7, 1.0,  "col8 n-stress = sigma_n_corr (compression+)");
   TEST_NEAR(mu_eff,   0.450, 1e-9, "col9 mu_eff = mu_d (δ > d_c)");

   std::remove(fname.c_str());
}

// ---------------------------------------------------------------------
// S-4  LSW mu_eff formula at the partial-weakening regime.
// ---------------------------------------------------------------------
static void S4_mu_eff_formula()
{
   std::cout << "\n[S-4] LSWFrictionCoefficient_TPV31 ramp\n";
   const real_t mu_s = 0.580, mu_d = 0.450, d_c = 0.18;
   TEST_NEAR(LSWFrictionCoefficient_TPV31(0.0,    mu_s, mu_d, d_c),
             mu_s, 1e-12, "δ=0 → mu_s");
   TEST_NEAR(LSWFrictionCoefficient_TPV31(d_c,    mu_s, mu_d, d_c),
             mu_d, 1e-12, "δ=d_c → mu_d");
   TEST_NEAR(LSWFrictionCoefficient_TPV31(2*d_c,  mu_s, mu_d, d_c),
             mu_d, 1e-12, "δ>d_c → mu_d");
   TEST_NEAR(LSWFrictionCoefficient_TPV31(0.5*d_c, mu_s, mu_d, d_c),
             mu_s - 0.5 * (mu_s - mu_d), 1e-12, "δ=d_c/2 → midpoint");
}

int main(int, char **)
{
   std::cout << "Running test_tpv31_station_writer\n";
   S1_station_list();
   S2_nearest_dof();
   S3_writestep_columns();
   S4_mu_eff_formula();
   std::cout << "\n========================================\n";
   std::cout << "test_tpv31_station_writer: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
