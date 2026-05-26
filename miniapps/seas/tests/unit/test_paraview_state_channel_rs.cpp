// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_paraview_state_channel_rs.cpp — Phase 3 (R-005) unit test for the
// ParaView fault "state" channel branch in spatial_dyn_driver.cpp.
//
// Plan: PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24, §9.3 req 10.
//
// The snapshot writer chooses the per-DOF "state" field by law:
//   is_lsw  -> LSWFrictionCoefficient_TPV205(delta_norm, lsw_mu_s/mu_d/d_c)
//   else    -> d.psi   (the RS state variable)
// InitializeFaultDOFs_Spatial_RS never sets the lsw_* fields (they stay 0),
// so the old unconditional LSW formula returned a uniform 0 for an RS run
// — the very channel meant to show the RS state was garbage.  This test
// mirrors the driver branch and asserts:
//   (a) RS  -> pv_state(i) == dof_data[i].psi for all i, finite at step 0;
//   (b) the old LSW formula on RS DOFData (lsw_* = 0) is NOT psi (0 here),
//       documenting why the branch is needed;
//   (c) LSW -> pv_state == LSWFrictionCoefficient_TPV205(...) byte-identical.

#include "mfem.hpp"

#include "../../dynamic/fault_face_flux.hpp"        // DOFData
#include "../../dynamic/tpv205_friction.hpp"        // LSWFrictionCoefficient_TPV205
#include "../../dynamic/fault_state_channel.hpp"    // R-024: the REAL state-channel selection

#include <cmath>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

namespace
{
// R-024: exercise the REAL driver code path.  The driver's ParaView writer and
// this test both call mfem::seas::FaultStateChannelValue (the shared selection
// in dynamic/fault_state_channel.hpp), so the test now guards the actual code
// instead of a local mirror that could silently drift.
real_t state_channel_value(bool is_lsw, const DOFData& d)
{
   return mfem::seas::FaultStateChannelValue(is_lsw, d);
}
}  // namespace

// =====================================================================
// (a) RS: pv_state == psi for all i, finite at step 0 (delta == 0).
// =====================================================================
static void rs_channel_writes_psi()
{
   std::cout << "\n[R-005a] RS state channel writes d.psi\n";
   const int N = 3;
   std::vector<DOFData> dof(N, DOFData{});
   // RS init shape: lsw_* stay 0 (InitializeFaultDOFs_Spatial_RS), psi set
   // by SeedEquilibriumPsi_RS; slip = 0 at step 0 -> delta_norm = 0.
   for (int i = 0; i < N; ++i)
   {
      dof[i].lsw_mu_s = dof[i].lsw_mu_d = dof[i].lsw_d_c = 0.0;
      dof[i].slip1 = dof[i].slip2 = 0.0;
      dof[i].psi = 0.45 + 0.05 * i;   // distinct finite seeded states
   }

   bool all_match = true, all_finite = true;
   for (int i = 0; i < N; ++i)
   {
      const real_t s = state_channel_value(/*is_lsw=*/false, dof[i]);
      all_match  = all_match  && (s == dof[i].psi);
      all_finite = all_finite && std::isfinite(s);
   }
   TEST_ASSERT(all_match, "RS pv_state(i) == dof_data[i].psi for all i");
   TEST_ASSERT(all_finite, "RS pv_state finite at step 0 (delta == 0)");
}

// =====================================================================
// (b) The old LSW formula on RS DOFData is garbage (not psi).
// =====================================================================
static void old_lsw_formula_on_rs_is_garbage()
{
   std::cout << "\n[R-005b] old LSW formula on RS DOFData != psi (why branch)\n";
   DOFData d{};
   d.lsw_mu_s = d.lsw_mu_d = d.lsw_d_c = 0.0;   // RS init: lsw_* unset
   d.slip1 = d.slip2 = 0.0;
   d.psi = 0.5;

   const real_t old_val = mfem::seas::LSWFrictionCoefficient_TPV205(
                             0.0, d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
   const real_t new_val = state_channel_value(/*is_lsw=*/false, d);
   TEST_ASSERT(old_val != d.psi,
               "old LSW formula on RS DOFData does NOT yield psi (garbage=0)");
   TEST_ASSERT(new_val == d.psi,
               "new RS branch yields psi (the fix)");
}

// =====================================================================
// (c) LSW branch byte-identical to LSWFrictionCoefficient_TPV205.
// =====================================================================
static void lsw_channel_byte_identical()
{
   std::cout << "\n[R-005c] LSW state channel byte-identical to pre-change\n";
   DOFData d{};
   d.lsw_mu_s = 0.6;
   d.lsw_mu_d = 0.3;
   d.lsw_d_c  = 0.4;
   d.slip1 = 0.1; d.slip2 = 0.0;   // delta_norm = 0.1 (mid-weakening)
   d.psi = 0.123;                  // must be IGNORED on the LSW path

   const real_t delta_norm = std::sqrt(d.slip1 * d.slip1 + d.slip2 * d.slip2);
   const real_t expected = mfem::seas::LSWFrictionCoefficient_TPV205(
                              delta_norm, d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
   const real_t got = state_channel_value(/*is_lsw=*/true, d);
   TEST_ASSERT(got == expected,
               "LSW pv_state == LSWFrictionCoefficient_TPV205 (byte-identical)");
   TEST_ASSERT(got != d.psi, "LSW path ignores psi");

   // Finite at step 0 (delta == 0) -> mu_s.
   DOFData d0{};
   d0.lsw_mu_s = 0.6; d0.lsw_mu_d = 0.3; d0.lsw_d_c = 0.4;
   d0.slip1 = d0.slip2 = 0.0;
   const real_t s0 = state_channel_value(/*is_lsw=*/true, d0);
   TEST_ASSERT(std::isfinite(s0) && s0 == 0.6,
               "LSW pv_state finite at step 0 (== mu_s)");
}

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase 3 test_paraview_state_channel_rs\n";
   rs_channel_writes_psi();
   old_lsw_formula_on_rs_is_garbage();
   lsw_channel_byte_identical();

   std::cout << "\n========================================\n";
   std::cout << "Phase 3 test_paraview_state_channel_rs: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}
