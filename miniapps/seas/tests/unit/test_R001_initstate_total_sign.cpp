// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 Phase 4 (I-06 part B): R-001 regression — verify the
// signs baked into InitializeStateTotal.
//
// Under MFEM's compression-positive convention (miniapps/seas/CLAUDE.md
// §"Normal stress") and BP5's canonical fault-local frame
//    can_n = (0,-1,0),  can_t1_dip = (0,0,-1),  can_t2_strike = (+1,0,0)
// the TPV102 pre-stress tensor in GLOBAL Cartesian coordinates has:
//    Q[SXX] = 0
//    Q[SYY] = +sigma_n0         (compressive, R-001 fix)
//    Q[SZZ] = 0
//    Q[SXY] = -tau_ini          (R-001 fix: strike shear picks up a sign
//                                flip from the canonical-local rotation
//                                since can_t2 = +x, can_n = -y)
//    Q[SYZ] = 0
//    Q[SXZ] = 0
//    VX = VY = VZ = 0
// at every bulk DOF.
//
// This test was flagged as the primary R-001 regression gate in the
// review-incorporation plan: pre-R-001 the plan wrote
//    Q[SYY] = -sigma_n0,   Q[SXY] = -tau_ini
// with the first sign in conflict with MFEM's compression-positive
// convention.  A bug that flips Q[SYY]'s sign would pass every ULP-based
// equivalence test (R-001 affects the physical state, not the flux
// kernel's arithmetic) — this test catches it directly.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

static void CheckQComponent(const Vector &Q, int ndof_total, int c,
                            real_t expected, const std::string &name)
{
   real_t max_dev = 0.0;
   int    worst_i = -1;
   for (int i = 0; i < ndof_total; i++)
   {
      real_t dev = std::abs(Q(c * ndof_total + i) - expected);
      if (dev > max_dev) { max_dev = dev; worst_i = i; }
   }
   std::cout << "  " << name << ": max|Q[" << c << "] - " << expected
             << "| = " << std::scientific << std::setprecision(3)
             << max_dev << " (worst DOF " << worst_i << ")\n";
   TEST_ASSERT(max_dev == 0.0,
               "R-001: " + name + " exactly equals " + std::to_string(expected));
}

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 Phase 4 (I-06): "
             << "R-001 InitializeStateTotal sign regression ===\n";

   // A handful of realistic pre-stress values per the TPV102 params.
   struct Case
   {
      real_t sigma_n0;
      real_t tau_ini;
      const char *label;
   };
   const Case cases[] = {
      {120.0e6, 75.0e6, "TPV102 default (sigma_n=120 MPa, tau=75 MPa)"},
      {  1.0,    1.0,   "unit dimensionless"},
      {  0.0,   50.0e6, "zero normal stress"},
      { 50.0e6,  0.0,   "zero shear pre-stress"},
   };

   const int ndof_total = 27;    // arbitrary, >1

   for (const Case &c : cases)
   {
      std::cout << "\n-- " << c.label << " --\n";
      Vector Q;
      InitializeStateTotal(Q, ndof_total, c.sigma_n0, c.tau_ini);

      TEST_ASSERT(Q.Size() == NUM_STATE * ndof_total,
                  "Q has size NUM_STATE * ndof_total");

      CheckQComponent(Q, ndof_total, SXX, 0.0,          "Q[SXX]");
      CheckQComponent(Q, ndof_total, SYY, +c.sigma_n0,  "Q[SYY] = +sigma_n0 (R-001)");
      CheckQComponent(Q, ndof_total, SZZ, 0.0,          "Q[SZZ]");
      CheckQComponent(Q, ndof_total, SXY, -c.tau_ini,   "Q[SXY] = -tau_ini (R-001)");
      CheckQComponent(Q, ndof_total, SYZ, 0.0,          "Q[SYZ]");
      CheckQComponent(Q, ndof_total, SXZ, 0.0,          "Q[SXZ]");
      CheckQComponent(Q, ndof_total, VX,  0.0,          "Q[VX]");
      CheckQComponent(Q, ndof_total, VY,  0.0,          "Q[VY]");
      CheckQComponent(Q, ndof_total, VZ,  0.0,          "Q[VZ]");
   }

   // ZeroDOFDataPreStressTotal sanity.
   {
      std::cout << "\n-- ZeroDOFDataPreStressTotal --\n";
      std::vector<DOFData> d(5);
      for (int i = 0; i < 5; i++)
      {
         d[i].sigma_n0 = 120e6;
         d[i].tau1_0   = 1.23e6;
         d[i].tau2_0   = 75e6;
         // sigma_n_corr / tau_corr seeded as the "physical pre-stress" per
         // the InitializeFaultDOFs convention — MUST be left alone by
         // ZeroDOFDataPreStressTotal.
         d[i].sigma_n_corr = 120e6;
         d[i].tau1_corr    = 1.23e6;
         d[i].tau2_corr    = 75e6;
      }
      ZeroDOFDataPreStressTotal(d, 5);
      for (int i = 0; i < 5; i++)
      {
         TEST_ASSERT(d[i].sigma_n0 == 0.0,
                     "DOFData[" + std::to_string(i) + "].sigma_n0 zeroed");
         TEST_ASSERT(d[i].tau1_0   == 0.0,
                     "DOFData[" + std::to_string(i) + "].tau1_0 zeroed");
         TEST_ASSERT(d[i].tau2_0   == 0.0,
                     "DOFData[" + std::to_string(i) + "].tau2_0 zeroed");
         TEST_ASSERT(d[i].sigma_n_corr == 120e6,
                     "DOFData[" + std::to_string(i) +
                     "].sigma_n_corr preserved (OUTPUT field)");
         TEST_ASSERT(d[i].tau2_corr == 75e6,
                     "DOFData[" + std::to_string(i) +
                     "].tau2_corr preserved (OUTPUT field)");
      }
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
