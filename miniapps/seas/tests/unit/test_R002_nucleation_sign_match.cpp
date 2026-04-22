// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 Phase 5 (I-06 part C): R-002 regression — the
// nucleation injection into bulk Q must use the correct sign.
//
// Under BP5 canonical fault-local frame (n = (0,-1,0), t2_strike = (+x)),
// the Voigt rotation to global coordinates maps
//     sigma_xy_global = -tau_nt2_local
// Adding dtau to tau_nt2_local therefore requires SUBTRACTING dtau from
// Q[SXY] in global coordinates.  In ApplyNucleationTotal this takes the
// form:
//     Q[SXY, fault_qp_DOF] = -(tau_ini + dtau(x, z, t))   (overwrite)
//
// A bug that uses "+=" instead of the correct "= -(tau_ini + dtau)"
// would drive the fault in the WRONG direction — tau_nt2_local would
// DECREASE, not increase, nucleation would never trigger rupture.
//
// This test injects a known dtau with ApplyNucleationTotal on a tiny
// synthetic fixture and verifies the resulting Q[SXY] matches the
// predicted negative-of-sum.

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
                << " (got " << v_ << ", expected " << e_ \
                << ", diff " << std::abs(v_-e_) << ")\n"; } \
} while (0)

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 Phase 5 (I-06): "
             << "R-002 nucleation sign regression ===\n";

   // Tiny synthetic fixture: 1 fault QP at the hypocenter location, both
   // + and - sides local.  dof_data order: (SXX..SXZ, VX..VZ) per DOF,
   // ndof_per_el = 2 (one for each side's nodal DOF).  Total 2 "elements"
   // x 2 DOFs per element = 4 DOFs per component.
   const int ndof_per_el = 2;
   const int ne          = 2;
   const int ndof_total  = ne * ndof_per_el;
   const int plus_elem   = 1;     // arbitrary choice
   const int minus_elem  = 0;
   const int plus_dof    = 0;     // nodal DOF closest to fault QP on + side
   const int minus_dof   = 1;     // nodal DOF closest to fault QP on - side

   // One fault QP AT the hypocenter so NucleationPerturbation gives a
   // non-zero dtau once t > 0.
   std::vector<Vector> fault_coords(1, Vector(3));
   fault_coords[0][0] = TPV102Params::hypo_along_strike;  // along-strike
   fault_coords[0][1] = 0.0;                               // on fault plane
   fault_coords[0][2] = -TPV102Params::hypo_down_dip;      // depth < 0

   // Interpolatory fixture: shape_max = 1.0 on both sides so the
   // add-delta + 1/shape_max scaling collapses to the original
   // overwrite magnitude (the add-delta test is R-I06-003 below).
   FaultQPNodalMap map;
   map.qp_to_elem_plus     = {plus_elem};
   map.qp_to_dof_plus      = {plus_dof};
   map.qp_to_elem_minus    = {minus_elem};
   map.qp_to_dof_minus     = {minus_dof};
   map.qp_shape_max_plus   = {1.0};
   map.qp_shape_max_minus  = {1.0};
   TEST_ASSERT(map.Consistent(),
               "FaultQPNodalMap::Consistent() on the synthetic fixture");

   // Per-QP add-delta state (R-I06-004).  Driver owns; reset to n_qp.
   FaultQPNucleationState nuc_state;
   nuc_state.Reset(1);

   // Initialize Q with the TPV102 pre-stress (Q[SXY] = -tau_ini everywhere).
   Vector Q;
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);

   // t = 0: dtau_new = 0, delta = 0 - dtau_applied(0) = 0.
   //        ApplyNucleationTotal must be a no-op; Q stays at InitializeStateTotal.
   {
      ApplyNucleationTotal(Q, map, nuc_state, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, 0.0);
      const int off_p = plus_elem  * ndof_per_el + plus_dof;
      const int off_m = minus_elem * ndof_per_el + minus_dof;
      TEST_NEAR(Q[SXY * ndof_total + off_p], -TPV102Params::tau_ini,
                1e-12 * TPV102Params::tau_ini,
                "R-002 t=0: Q[SXY, + side QP] = -tau_ini");
      TEST_NEAR(Q[SXY * ndof_total + off_m], -TPV102Params::tau_ini,
                1e-12 * TPV102Params::tau_ini,
                "R-002 t=0: Q[SXY, - side QP] = -tau_ini");
      TEST_NEAR(nuc_state.dtau_applied_plus[0], 0.0, 1e-12,
                "R-I06-004: nuc_state tracks dtau_applied (t=0 → 0)");
   }

   // t = nuc_T+: dtau_new saturates to nuc_dtau, delta = nuc_dtau - 0 =
   //            nuc_dtau.  Q[SXY] should add -nuc_dtau/shape_max =
   //            -nuc_dtau (shape_max = 1 here).  Total Q[SXY] =
   //            -tau_ini - nuc_dtau.
   {
      const real_t t = TPV102Params::nuc_T + 0.5;  // past the ramp
      const real_t expected = -(TPV102Params::tau_ini + TPV102Params::nuc_dtau);
      ApplyNucleationTotal(Q, map, nuc_state, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t);
      const int off_p = plus_elem  * ndof_per_el + plus_dof;
      const int off_m = minus_elem * ndof_per_el + minus_dof;
      std::cout << "  at t=" << t << " s, expected Q[SXY] = "
                << std::scientific << std::setprecision(3) << expected << " Pa\n";
      std::cout << "    + side: Q[SXY] = " << Q[SXY * ndof_total + off_p] << "\n";
      std::cout << "    - side: Q[SXY] = " << Q[SXY * ndof_total + off_m] << "\n";
      TEST_NEAR(Q[SXY * ndof_total + off_p], expected,
                1e-12 * std::abs(expected),
                "R-002 t=nuc_T+: Q[SXY, + side QP] = -(tau_ini + nuc_dtau)");
      TEST_NEAR(Q[SXY * ndof_total + off_m], expected,
                1e-12 * std::abs(expected),
                "R-002 t=nuc_T+: Q[SXY, - side QP] = -(tau_ini + nuc_dtau)");
      TEST_NEAR(nuc_state.dtau_applied_plus[0], TPV102Params::nuc_dtau,
                1e-12 * TPV102Params::nuc_dtau,
                "R-I06-004: nuc_state tracks dtau_applied (t>nuc_T → nuc_dtau)");
   }

   // Idempotency check (R-I06-004): calling ApplyNucleationTotal AGAIN
   // with the same t should be a no-op (delta = 0).  If the old
   // overwrite semantic leaked back in, we'd see Q[SXY] reset to
   // -(tau_ini + nuc_dtau) regardless of previous wave perturbations.
   {
      const real_t t = TPV102Params::nuc_T + 0.5;
      const int off_p = plus_elem  * ndof_per_el + plus_dof;
      // Simulate a "wave perturbation" landing on the nearest DOF.
      const real_t wave_perturb = 1.5e5;
      Q[SXY * ndof_total + off_p] += wave_perturb;
      const real_t before = Q[SXY * ndof_total + off_p];
      ApplyNucleationTotal(Q, map, nuc_state, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t);
      const real_t after = Q[SXY * ndof_total + off_p];
      TEST_NEAR(after, before, 1e-9 * std::abs(before),
                "R-I06-004: ApplyNucleationTotal preserves wave "
                "perturbations at the fault-QP DOF (add-delta, not "
                "overwrite)");
   }

   // Bug check: confirm the result is NEGATIVE.  A wrong-sign implementation
   // (Q[SXY] = +(tau_ini + nuc_dtau)) would drive rupture in the wrong
   // direction; even a magnitude-only test would miss this.
   {
      const int off_p = plus_elem  * ndof_per_el + plus_dof;
      TEST_ASSERT(Q[SXY * ndof_total + off_p] < 0.0,
                  "R-002: Q[SXY] is negative (sign flip from -tau_nt2_local)");
   }

   // Far-field DOFs (not in the map) must NOT be touched.
   {
      // element 0 (minus_elem) has DOFs {0, 1}.  DOF 0 is NOT in the map
      // (map uses minus_dof = 1).  It should still carry InitializeStateTotal's
      // -tau_ini value.
      const int off_far = 0 * ndof_per_el + 0;  // not in map
      TEST_NEAR(Q[SXY * ndof_total + off_far], -TPV102Params::tau_ini,
                1e-12 * TPV102Params::tau_ini,
                "R-002: non-fault-QP DOF untouched by ApplyNucleationTotal");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
