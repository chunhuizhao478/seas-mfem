// Part C / C2: TPV6 per-side on-fault station writer (dynamic/tpv6_stations.hpp).
//
// Verifies, WITHOUT a mesh/driver (the writer is dof_data-only + header-only):
//   - DefaultStations_TPV6 = the 5 drdg3d stations.
//   - nearside trace = the STRONG side (larger Zp), farside = WEAK side, robust to
//     which DOFData side (plus/minus) is strong.
//   - component mapping h=strike(v_imp[2]) / v=dip(v_imp[1]) / n=normal(v_imp[0]);
//     stresses in MPa (compression-positive sigma_n).
//   - displacement = trapezoidal time integral of the per-side velocity.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/tpv6_stations.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
#define TEST_NEAR(v, ref, tol, msg) do {                                       \
   num_tests++; const double e = std::abs((double)(v) - (double)(ref));        \
   if (e <= (tol)) { num_passed++; std::cout << "  PASSED: " << msg            \
        << " (" << (double)(v) << ")\n"; }                                     \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg  \
        << " got " << (double)(v) << " want " << (double)(ref) << "\n"; }      \
} while (0)
#define TEST_TRUE(c, msg) do { num_tests++; if (c) { num_passed++; \
   std::cout << "  PASSED: " << msg << "\n"; } else { num_failed++; \
   std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } } while (0)

namespace
{
// Read the LAST data row (non-#) of a station file into 10 doubles.
bool LastRow(const std::string &path, double col[10])
{
   std::ifstream f(path);
   if (!f.is_open()) { return false; }
   std::string line, last;
   while (std::getline(f, line))
   {
      if (!line.empty() && line[0] != '#') { last = line; }
   }
   if (last.empty()) { return false; }
   std::istringstream is(last);
   for (int i = 0; i < 10; ++i) { if (!(is >> col[i])) { return false; } }
   return true;
}
} // namespace

int main()
{
   std::cout << "\n=== Part C: TPV6 per-side station writer ===\n";
   const char *td = std::getenv("TMPDIR");
   const std::string dir = (td && td[0]) ? std::string(td) : std::string("/tmp");

   // -- DefaultStations --
   auto stations = DefaultStations_TPV6();
   TEST_TRUE(stations.size() == 5u, "DefaultStations_TPV6 = 5 stations");
   TEST_TRUE(stations[0].name == "x2_0_x3_0" && stations[0].along_strike == 0.0
             && stations[0].down_dip == 0.0, "station[0] = (0,0) x2_0_x3_0");

   // Synthetic fault DOFs co-located with the stations (so dof[s] <-> station[s]).
   std::vector<Vector> coords(stations.size());
   for (size_t s = 0; s < stations.size(); ++s)
   {
      coords[s].SetSize(3);
      coords[s](0) = stations[s].along_strike;
      coords[s](1) = 0.0;
      coords[s](2) = -stations[s].down_dip;   // depth -> z (z<0 below)
   }

   // DOFData: minus side = STRONG (larger Zp); plus = weak.
   std::vector<DOFData> dof(stations.size());
   for (auto &d : dof)
   {
      d.Zp_minus = 1.602e7; d.Zp_plus = 8.34e6;     // minus is near/strong
      d.v_imp_minus[0] = 0.0; d.v_imp_minus[1] = 0.1; d.v_imp_minus[2] = 0.5;  // n,dip,strike
      d.v_imp_plus[0]  = 0.0; d.v_imp_plus[1] = -0.2; d.v_imp_plus[2]  = 2.0;
      d.tau2_corr = 70.0e6; d.tau1_corr = 0.0; d.sigma_n_corr = 120.0e6;
   }

   {
      TPV6StationWriter w;
      w.Open(dir, "utest_tpv6", stations, coords, (int)dof.size());
      w.WriteStep(0.0, dof);
      w.WriteStep(0.1, dof);   // constant v -> disp = v * 0.1
      w.Close();
   }

   double nr[10], fr[10];
   const std::string base = dir + "/utest_tpv6_";
   TEST_TRUE(LastRow(base + "nearside_x2_0_x3_0.dat", nr), "nearside file written");
   TEST_TRUE(LastRow(base + "farside_x2_0_x3_0.dat",  fr), "farside file written");

   // cols: 0=t 1=h-disp 2=h-vel 3=h-stress 4=v-disp 5=v-vel 6=v-stress 7=n-disp 8=n-vel 9=n-stress
   TEST_NEAR(nr[0], 0.1, 1e-12, "t == 0.1 (last row)");
   TEST_NEAR(nr[2], 0.5, 1e-12, "nearside h-vel == strong-side strike vel (0.5)");
   TEST_NEAR(fr[2], 2.0, 1e-12, "farside  h-vel == weak-side strike vel (2.0)");
   TEST_NEAR(nr[3], 70.0, 1e-9, "h-stress == tau2_corr in MPa (70)");
   TEST_NEAR(nr[9], 120.0, 1e-9, "n-stress == sigma_n_corr in MPa, compression+ (120)");
   TEST_NEAR(nr[1], 0.05, 1e-12, "nearside h-disp == trapezoid(0.5)*0.1 = 0.05");
   TEST_NEAR(fr[1], 0.20, 1e-12, "farside  h-disp == 2.0*0.1 = 0.20");
   TEST_NEAR(nr[8], 0.0, 1e-12, "nearside n-vel == normal vel (0)");
   TEST_NEAR(nr[5], 0.1, 1e-12, "nearside v-vel == dip vel (0.1)");

   // -- Robustness: FLIP so plus is strong -> nearside must follow the strong side --
   {
      std::vector<DOFData> dof2 = dof;
      for (auto &d : dof2)
      {
         d.Zp_minus = 8.34e6; d.Zp_plus = 1.602e7;   // now PLUS is near/strong
         // strong (plus) strike vel = 0.5, weak (minus) = 2.0 (swap the roles)
         d.v_imp_plus[2]  = 0.5;
         d.v_imp_minus[2] = 2.0;
      }
      TPV6StationWriter w;
      w.Open(dir, "utest_tpv6b", stations, coords, (int)dof2.size());
      w.WriteStep(0.0, dof2);
      w.WriteStep(0.1, dof2);
      w.Close();
      double nr2[10], fr2[10];
      LastRow(dir + "/utest_tpv6b_nearside_x2_0_x3_0.dat", nr2);
      LastRow(dir + "/utest_tpv6b_farside_x2_0_x3_0.dat",  fr2);
      TEST_NEAR(nr2[2], 0.5, 1e-12, "FLIP: nearside still = STRONG side (0.5)");
      TEST_NEAR(fr2[2], 2.0, 1e-12, "FLIP: farside still = WEAK side (2.0)");
   }

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, " << num_failed
             << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
