// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 Phase 5 (I-06 part C): R-009 regression — nucleation
// injection must update BOTH the + and - side element DOFs at each fault
// QP.
//
// Rationale (review-incorporation plan, R-009): the trial traction in
// EvaluateTotal averages Q+[SXY]/Zs + Q-[SXY]/Zs across the two sides
// (Pelties eq. 7c).  If nucleation updates only ONE side, the perturbation
// on the trial traction is HALVED, which silently reduces the effective
// nucleation amplitude by a factor of 2 — the rupture takes roughly
// twice as long to initiate.
//
// This test covers three scenarios:
//   1. Both sides local (interior fault face): both DOFs updated.
//   2. Only + side local (shared fault face on rank A): only + DOF
//      updated; - is -1 ("peer rank"); unaffected local DOFs stay at
//      initial value.
//   3. Only - side local (shared fault face on rank B): same for -.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) < t_) { num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << " (got " << v_ << ", expected " << e_ << ")\n"; } \
} while (0)

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 Phase 5 (I-06): "
             << "R-009 nucleation both-sides regression ===\n";

   // Fixture: 2 "elements" (representing + and - sides of a fault face),
   // 2 DOFs per element.
   const int ndof_per_el = 2;
   const int ne          = 2;
   const int ndof_total  = ne * ndof_per_el;
   const int plus_elem   = 1;
   const int minus_elem  = 0;
   const int plus_dof    = 0;
   const int minus_dof   = 1;

   std::vector<Vector> fault_coords(1, Vector(3));
   fault_coords[0][0] = TPV102Params::hypo_along_strike;
   fault_coords[0][1] = 0.0;
   fault_coords[0][2] = -TPV102Params::hypo_down_dip;

   const real_t t_saturated = TPV102Params::nuc_T + 0.5;
   const real_t expected_after = -(TPV102Params::tau_ini + TPV102Params::nuc_dtau);

   // Scenario 1: both sides local (interior fault face).
   {
      std::cout << "\n-- Scenario 1: both sides local --\n";
      FaultQPNodalMap map;
      map.qp_to_elem_plus    = {plus_elem};
      map.qp_to_dof_plus     = {plus_dof};
      map.qp_to_elem_minus   = {minus_elem};
      map.qp_to_dof_minus    = {minus_dof};
      map.qp_shape_max_plus  = {1.0};  // interpolatory fixture
      map.qp_shape_max_minus = {1.0};
      TEST_ASSERT(map.Consistent(), "map consistent");

      Vector Q;
      InitializeStateTotal(Q, ndof_total,
                           TPV102Params::sigma_n, TPV102Params::tau_ini);
      FaultQPNucleationState st;
      st.Reset(1);
      ApplyNucleationTotal(Q, map, st, fault_coords, ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_saturated);

      const int off_p = plus_elem  * ndof_per_el + plus_dof;
      const int off_m = minus_elem * ndof_per_el + minus_dof;
      TEST_NEAR(Q[SXY * ndof_total + off_p], expected_after,
                1e-6, "R-009: + side DOF updated");
      TEST_NEAR(Q[SXY * ndof_total + off_m], expected_after,
                1e-6, "R-009: - side DOF updated");
   }

   // Scenario 2: only + side local (e.g., shared fault face on rank A).
   // map.qp_to_elem_minus = -1 signals "peer rank handles the - side".
   {
      std::cout << "\n-- Scenario 2: only + side local --\n";
      FaultQPNodalMap map;
      map.qp_to_elem_plus    = {plus_elem};
      map.qp_to_dof_plus     = {plus_dof};
      map.qp_to_elem_minus   = {-1};
      map.qp_to_dof_minus    = {-1};
      map.qp_shape_max_plus  = {1.0};
      map.qp_shape_max_minus = {0.0};  // non-local; shape_max irrelevant
      TEST_ASSERT(map.Consistent(), "map consistent");

      Vector Q;
      InitializeStateTotal(Q, ndof_total,
                           TPV102Params::sigma_n, TPV102Params::tau_ini);
      FaultQPNucleationState st;
      st.Reset(1);
      ApplyNucleationTotal(Q, map, st, fault_coords, ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_saturated);

      const int off_p = plus_elem  * ndof_per_el + plus_dof;
      TEST_NEAR(Q[SXY * ndof_total + off_p], expected_after,
                1e-6, "R-009: + side DOF updated");

      // Every other DOF should remain at the initial -tau_ini (peer rank
      // updates its own local DOF for the - side; on THIS rank we don't
      // touch anything else).
      for (int e = 0; e < ne; e++)
      {
         for (int d = 0; d < ndof_per_el; d++)
         {
            const int off = e * ndof_per_el + d;
            if (off == off_p) { continue; }
            TEST_NEAR(Q[SXY * ndof_total + off], -TPV102Params::tau_ini,
                      1e-6,
                      "R-009 (+ only local): non-plus DOF unchanged "
                      "(elem=" + std::to_string(e) +
                      " dof=" + std::to_string(d) + ")");
         }
      }
   }

   // Scenario 3: only - side local.
   {
      std::cout << "\n-- Scenario 3: only - side local --\n";
      FaultQPNodalMap map;
      map.qp_to_elem_plus    = {-1};
      map.qp_to_dof_plus     = {-1};
      map.qp_to_elem_minus   = {minus_elem};
      map.qp_to_dof_minus    = {minus_dof};
      map.qp_shape_max_plus  = {0.0};  // non-local
      map.qp_shape_max_minus = {1.0};
      TEST_ASSERT(map.Consistent(), "map consistent");

      Vector Q;
      InitializeStateTotal(Q, ndof_total,
                           TPV102Params::sigma_n, TPV102Params::tau_ini);
      FaultQPNucleationState st;
      st.Reset(1);
      ApplyNucleationTotal(Q, map, st, fault_coords, ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_saturated);

      const int off_m = minus_elem * ndof_per_el + minus_dof;
      TEST_NEAR(Q[SXY * ndof_total + off_m], expected_after,
                1e-6, "R-009: - side DOF updated");

      for (int e = 0; e < ne; e++)
      {
         for (int d = 0; d < ndof_per_el; d++)
         {
            const int off = e * ndof_per_el + d;
            if (off == off_m) { continue; }
            TEST_NEAR(Q[SXY * ndof_total + off], -TPV102Params::tau_ini,
                      1e-6,
                      "R-009 (- only local): non-minus DOF unchanged "
                      "(elem=" + std::to_string(e) +
                      " dof=" + std::to_string(d) + ")");
         }
      }
   }

   // Scenario 4: map inconsistency rejected.
   {
      std::cout << "\n-- Scenario 4: inconsistency guard --\n";
      FaultQPNodalMap bad;
      bad.qp_to_elem_plus    = {plus_elem};
      bad.qp_to_dof_plus     = {-1};           // mismatch: elem valid, dof not
      bad.qp_to_elem_minus   = {minus_elem};
      bad.qp_to_dof_minus    = {minus_dof};
      bad.qp_shape_max_plus  = {1.0};
      bad.qp_shape_max_minus = {1.0};
      TEST_ASSERT(!bad.Consistent(),
                  "R-009: invalid dof with valid elem rejected");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
