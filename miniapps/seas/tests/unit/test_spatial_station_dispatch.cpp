// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_station_dispatch.cpp — coverage for the per-problem on-fault
// station writers that seas_spatial_dyn_driver selects by [problem].tag
// (tpv205 / tpv104 / tpv102).  The tpv31 writer is covered separately by
// test_tpv31_station_writer.cpp; this test exercises the three writers that
// were extracted into the lean dynamic/tpv{205,104,102}_stations.hpp headers
// and newly wired into the spatial driver.
//
// For each problem it checks:
//   (a) DefaultStations_* returns the SCEC station count, and
//   (b) the writer emits the benchmark trace columns in the correct order
//       (SCEC layout, BP5 frame: h = strike = component 2, v = dip =
//       component 1), with the problem-specific trailing channel
//       (tpv205: LSW mu_eff; tpv104: psi; tpv102: log10_theta).
//
// Distinct sentinel field values are used so a column-order swap is caught.

#include "mfem.hpp"

#include "../../dynamic/tpv205_stations.hpp"
#include "../../dynamic/tpv104_stations.hpp"
#include "../../dynamic/tpv102_stations.hpp"

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

// Read the last non-comment line of a station file; return false if absent.
static bool LastDataLine(const std::string &fname, std::string &out)
{
   std::ifstream in(fname);
   if (!in.good()) { return false; }
   std::string line;
   bool have = false;
   while (std::getline(in, line))
   {
      if (!line.empty() && line[0] != '#') { out = line; have = true; }
   }
   return have;
}

// One DOF at (0, 7500); single station there; distinct sentinel fields.
static void MakeSingleStationDOF(std::vector<Vector> &coords,
                                 std::vector<DOFData> &dof)
{
   coords.resize(1);
   coords[0].SetSize(3);
   coords[0](0) = 0.0; coords[0](1) = 0.0; coords[0](2) = -7500.0;

   dof.resize(1);
   DOFData &d = dof[0];
   d.slip2     = 2.0;        // h-slip   (strike)
   d.V2        = 3.0;        // h-slip-rate
   d.tau2_corr = 4.0e7;      // h-shear-stress (Pa)
   d.slip1     = 0.5;        // v-slip   (dip)  — nonzero to catch a swap
   d.V1        = 0.7;        // v-slip-rate
   d.tau1_corr = 1.1e7;      // v-shear-stress
   d.sigma_n_corr = 6.0e7;   // n-stress (Pa, compression-positive)
}

// ---------------------------------------------------------------------
// D-1  TPV205 (LSW): 16 stations, SCEC TPV5 columns, trailing mu_eff(δ).
// ---------------------------------------------------------------------
static void D1_tpv205()
{
   std::cout << "\n[D-1] TPV205 station writer (LSW)\n";
   const std::vector<TPV205Station> st = DefaultStations_TPV205();
   TEST_ASSERT(st.size() == 16, "DefaultStations_TPV205 returns 16 stations");

   std::vector<TPV205Station> one(1);
   one[0] = { 0.0, 7500.0, "x2_0_x3_7.5" };

   std::vector<Vector> coords;
   std::vector<DOFData> dof;
   MakeSingleStationDOF(coords, dof);
   dof[0].lsw_mu_s = 0.677;   // SCEC TPV205
   dof[0].lsw_mu_d = 0.525;
   dof[0].lsw_d_c  = 0.40;
   // δ = sqrt(2² + 0.5²) = 2.06 > d_c → mu_eff = mu_d = 0.525.

   const std::string prefix = "test_dispatch_tpv205_tmp";
   const std::string fname = "./" + prefix + "_station_x2_0_x3_7.5.dat";
   std::remove(fname.c_str());
   {
      TPV205StationWriter w;
      w.Open(".", prefix, one, coords, /*ndof=*/1);
      w.WriteStep(0.5, dof);
      w.Flush();
      w.Close();
   }

   std::string data;
   TEST_ASSERT(LastDataLine(fname, data), "tpv205 station file created");
   if (!data.empty())
   {
      std::istringstream iss(data);
      real_t t, hs, hr, hsh, vs, vr, vsh, ns, mu;
      iss >> t >> hs >> hr >> hsh >> vs >> vr >> vsh >> ns >> mu;
      TEST_NEAR(t,   0.5,   1e-9, "tpv205 col1 t");
      TEST_NEAR(hs,  2.0,   1e-9, "tpv205 col2 h-slip = slip2");
      TEST_NEAR(hr,  3.0,   1e-9, "tpv205 col3 h-slip-rate = V2");
      TEST_NEAR(hsh, 4.0e7, 1.0,  "tpv205 col4 h-shear = tau2_corr");
      TEST_NEAR(vs,  0.5,   1e-9, "tpv205 col5 v-slip = slip1");
      TEST_NEAR(vr,  0.7,   1e-9, "tpv205 col6 v-slip-rate = V1");
      TEST_NEAR(vsh, 1.1e7, 1.0,  "tpv205 col7 v-shear = tau1_corr");
      TEST_NEAR(ns,  6.0e7, 1.0,  "tpv205 col8 n-stress = sigma_n_corr");
      TEST_NEAR(mu,  0.525, 1e-9, "tpv205 col9 mu_eff = mu_d (δ > d_c)");
   }
   std::remove(fname.c_str());
}

// ---------------------------------------------------------------------
// D-2  TPV104 (rate-state): 9 stations, SCEC columns, trailing psi.
// ---------------------------------------------------------------------
static void D2_tpv104()
{
   std::cout << "\n[D-2] TPV104 station writer (rate-state)\n";
   const std::vector<TPV104Station> st = DefaultStations_TPV104();
   TEST_ASSERT(st.size() == 9, "DefaultStations_TPV104 returns 9 stations");

   std::vector<TPV104Station> one(1);
   one[0] = { 0.0, 7500.0, "x2_0_x3_7.5" };

   std::vector<Vector> coords;
   std::vector<DOFData> dof;
   MakeSingleStationDOF(coords, dof);
   dof[0].psi = 0.55;   // trailing channel for TPV104

   const std::string prefix = "test_dispatch_tpv104_tmp";
   const std::string fname = "./" + prefix + "_station_x2_0_x3_7.5.dat";
   std::remove(fname.c_str());
   {
      TPV104StationWriter w;
      w.Open(".", prefix, one, coords, /*ndof=*/1);
      w.WriteStep(0.5, dof);
      w.Flush();
      w.Close();
   }

   std::string data;
   TEST_ASSERT(LastDataLine(fname, data), "tpv104 station file created");
   if (!data.empty())
   {
      std::istringstream iss(data);
      real_t t, hs, hr, hsh, vs, vr, vsh, ns, psi;
      iss >> t >> hs >> hr >> hsh >> vs >> vr >> vsh >> ns >> psi;
      TEST_NEAR(t,   0.5,   1e-9, "tpv104 col1 t");
      TEST_NEAR(hs,  2.0,   1e-9, "tpv104 col2 h-slip = slip2");
      TEST_NEAR(hr,  3.0,   1e-9, "tpv104 col3 h-slip-rate = V2");
      TEST_NEAR(hsh, 4.0e7, 1.0,  "tpv104 col4 h-shear = tau2_corr");
      TEST_NEAR(vs,  0.5,   1e-9, "tpv104 col5 v-slip = slip1");
      TEST_NEAR(vr,  0.7,   1e-9, "tpv104 col6 v-slip-rate = V1");
      TEST_NEAR(vsh, 1.1e7, 1.0,  "tpv104 col7 v-shear = tau1_corr");
      TEST_NEAR(ns,  6.0e7, 1.0,  "tpv104 col8 n-stress = sigma_n_corr");
      TEST_NEAR(psi, 0.55,  1e-9, "tpv104 col9 = psi (rate-state)");
   }
   std::remove(fname.c_str());
}

// ---------------------------------------------------------------------
// D-3  TPV102 (rate-state): 9 stations, generic DefaultStations(),
//      columns time slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta.
// ---------------------------------------------------------------------
static void D3_tpv102()
{
   std::cout << "\n[D-3] TPV102 station writer (rate-state, log10_theta)\n";
   const std::vector<TPV102Station> st = DefaultStations();
   TEST_ASSERT(st.size() == 9, "DefaultStations (tpv102) returns 9 stations");

   std::vector<TPV102Station> one(1);
   one[0] = { 0.0, 7500.0, "flt_0_7.5" };

   std::vector<Vector> coords;
   std::vector<DOFData> dof;
   MakeSingleStationDOF(coords, dof);
   dof[0].psi = 0.735;
   dof[0].Dc  = TPV102Params::Dc;
   // Expected trailing column reproduces the writer's psi->theta conversion.
   const real_t theta = (dof[0].Dc / TPV102Params::V0)
                        * std::exp((dof[0].psi - TPV102Params::f0)
                                   / TPV102Params::b);
   const real_t expect_log10_theta = std::log10(std::max(theta, 1e-300));

   const std::string prefix = "test_dispatch_tpv102_tmp";
   const std::string fname = "./" + prefix + "_station_flt_0_7.5.dat";
   std::remove(fname.c_str());
   {
      TPV102StationWriter w;
      w.Open(".", prefix, one, coords, /*ndof=*/1);
      w.WriteStep(0.5, dof);
      w.Flush();
      w.Close();
   }

   std::string data;
   TEST_ASSERT(LastDataLine(fname, data), "tpv102 station file created");
   if (!data.empty())
   {
      // TPV102 layout: time slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta.
      std::istringstream iss(data);
      real_t t, s1, s2, v1, v2, tau1, tau2, ns, lt;
      iss >> t >> s1 >> s2 >> v1 >> v2 >> tau1 >> tau2 >> ns >> lt;
      TEST_NEAR(t,    0.5,   1e-9, "tpv102 col1 t");
      TEST_NEAR(s1,   0.5,   1e-9, "tpv102 col2 slip1");
      TEST_NEAR(s2,   2.0,   1e-9, "tpv102 col3 slip2");
      TEST_NEAR(v1,   0.7,   1e-9, "tpv102 col4 V1");
      TEST_NEAR(v2,   3.0,   1e-9, "tpv102 col5 V2");
      TEST_NEAR(tau1, 1.1e7, 1.0,  "tpv102 col6 tau1_corr");
      TEST_NEAR(tau2, 4.0e7, 1.0,  "tpv102 col7 tau2_corr");
      TEST_NEAR(ns,   6.0e7, 1.0,  "tpv102 col8 sigma_n_corr");
      TEST_NEAR(lt, expect_log10_theta, 1e-6, "tpv102 col9 log10_theta");
   }
   std::remove(fname.c_str());
}

int main(int, char **)
{
   std::cout << "Running test_spatial_station_dispatch\n";
   D1_tpv205();
   D2_tpv104();
   D3_tpv102();
   std::cout << "\n========================================\n";
   std::cout << "test_spatial_station_dispatch: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
