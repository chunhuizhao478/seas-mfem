// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 R-I06-003 regression gate: ApplyNucleationTotal must
// deliver the requested dtau as the SHAPE-WEIGHTED sum at the fault QP,
// not just as a raw write to one DOF.  Under the original overwrite
// implementation the QP-effective amplitude was shape_max · dtau, i.e.
// under-amplified by (1 - shape_max) at every non-interpolatory face QP.
//
// Shape-max correction (R-I06-003): ApplyNucleationTotal now adds
// -delta / shape_max to the nearest-nodal DOF.  The shape-weighted QP
// sum then increments by -delta exactly (subject to the semantic that
// other DOFs carry their current values — under add-delta those are
// preserved between calls, so successive dtau increments stack
// linearly).
//
// This test constructs a synthetic fixture with a known shape_max < 1
// and checks:
//   1. The DOF value after ApplyNucleationTotal reflects the
//      shape_max-scaled delta (direct arithmetic check).
//   2. The shape-weighted QP sum equals the requested -(tau_ini + dtau)
//      — the physical QP amplitude the friction solver would see.
//   3. Iterative calls at successive dtau values produce stacking
//      deltas that each correctly update the QP-effective amplitude.
//
// On the ORIGINAL (overwrite, no shape_max) implementation this test
// would fail case 2 by a factor of ~shape_max on the QP-effective sum.

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
                << ", |diff| " << std::abs(v_-e_) << ")\n"; } \
} while (0)

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 R-I06-003 (I-06): "
             << "Nucleation shape-max amplitude scaling ===\n";

   // Fixture: 1 element with 3 DOFs.  DOF 0 is "nearest-nodal" to a
   // synthetic face QP where shape_max = 0.55 (simulating a non-
   // interpolatory Dunavant QP typical of order=2 tet face quadrature).
   // DOFs 1, 2 are "other" DOFs; the shape-weighted sum at the QP is
   //    SXY(QP) = shape(0)·Q[SXY,0] + shape(1)·Q[SXY,1] + shape(2)·Q[SXY,2]
   // with shape(0) = 0.55, shape(1) = shape(2) = 0.225  (sums to 1.0).
   const int ndof_per_el = 3;
   const int ne          = 1;
   const int ndof_total  = ne * ndof_per_el;
   const int plus_elem   = 0;
   const int plus_dof    = 0;     // the nearest-nodal DOF
   const real_t shape_max  = 0.55;
   const real_t shape_rest = (1.0 - shape_max) * 0.5;   // 0.225 each

   std::vector<Vector> fault_coords(1, Vector(3));
   fault_coords[0][0] = TPV102Params::hypo_along_strike;
   fault_coords[0][1] = 0.0;
   fault_coords[0][2] = -TPV102Params::hypo_down_dip;

   FaultQPNodalMap map;
   map.qp_to_elem_plus    = {plus_elem};
   map.qp_to_dof_plus     = {plus_dof};
   map.qp_to_elem_minus   = {-1};    // synthetic single-sided fixture
   map.qp_to_dof_minus    = {-1};
   map.qp_shape_max_plus  = {shape_max};
   map.qp_shape_max_minus = {0.0};

   // Compute QP shape-weighted sum helper.
   auto qp_sxy = [&](const Vector &Q) -> real_t {
      const int base = plus_elem * ndof_per_el;
      return shape_max * Q[SXY * ndof_total + base + 0]
           + shape_rest * Q[SXY * ndof_total + base + 1]
           + shape_rest * Q[SXY * ndof_total + base + 2];
   };

   // Initialize: Q[SXY, every DOF] = -tau_ini.
   Vector Q;
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);

   FaultQPNucleationState st;
   st.Reset(1);

   // Check 1: before any nucleation, QP-weighted sum = -tau_ini.
   TEST_NEAR(qp_sxy(Q), -TPV102Params::tau_ini,
             1e-6 * TPV102Params::tau_ini,
             "R-I06-003 t=0: shape-weighted QP sum = -tau_ini");

   // Check 2: apply nucleation at t = nuc_T + 0.5 (dtau saturates to
   // nuc_dtau since we are AT the hypocenter).
   const real_t t_saturated = TPV102Params::nuc_T + 0.5;
   const real_t dtau        = TPV102Params::nuc_dtau;
   ApplyNucleationTotal(Q, map, st, fault_coords,
                        ndof_per_el, ndof_total,
                        TPV102Params::tau_ini, t_saturated);

   // Direct DOF check: the nearest-nodal DOF should have been
   //    Q[SXY, 0] += -dtau / shape_max
   // from its initial value -tau_ini.
   const real_t expected_dof_p = -TPV102Params::tau_ini - dtau / shape_max;
   TEST_NEAR(Q[SXY * ndof_total + 0], expected_dof_p,
             1e-6 * std::abs(expected_dof_p),
             "R-I06-003: Q[SXY, nearest DOF] += -dtau / shape_max");

   // QP-EFFECTIVE check: the shape-weighted sum at the QP should equal
   // -(tau_ini + dtau) — this is the value the friction solver sees.
   //   = shape_max · (-tau_ini - dtau/shape_max) + 2·shape_rest · (-tau_ini)
   //   = -shape_max · tau_ini - dtau - (1 - shape_max) · tau_ini
   //   = -tau_ini - dtau     ✓
   TEST_NEAR(qp_sxy(Q), -(TPV102Params::tau_ini + dtau),
             1e-6 * std::abs(TPV102Params::tau_ini + dtau),
             "R-I06-003: shape-weighted QP sum = -(tau_ini + dtau) — "
             "the amplitude the friction solver actually sees");

   // Check 3: iterative re-application (add-delta).  Apply again at
   // t_saturated (delta = 0 since dtau is already applied); QP should
   // be unchanged.
   const real_t qp_before = qp_sxy(Q);
   ApplyNucleationTotal(Q, map, st, fault_coords,
                        ndof_per_el, ndof_total,
                        TPV102Params::tau_ini, t_saturated);
   TEST_NEAR(qp_sxy(Q), qp_before, 1e-9 * std::abs(qp_before),
             "R-I06-003 + R-I06-004: idempotent at same t (delta = 0)");

   // Check 4: partial ramp before full saturation.  At t inside the
   // nucleation ramp, dtau < nuc_dtau; the QP-effective sum must still
   // match -(tau_ini + dtau(t)).  Reset Q and state for this.
   {
      Vector Q2;
      InitializeStateTotal(Q2, ndof_total,
                           TPV102Params::sigma_n, TPV102Params::tau_ini);
      FaultQPNucleationState st2;
      st2.Reset(1);
      const real_t t_mid = TPV102Params::nuc_T * 0.5;
      const real_t dtau_mid = NucleationPerturbation(
         fault_coords[0][0], std::abs(fault_coords[0][2]), t_mid);
      std::cout << "  mid-ramp t = " << t_mid << " s, dtau_mid = "
                << std::scientific << std::setprecision(3) << dtau_mid
                << " Pa (nuc_dtau=" << TPV102Params::nuc_dtau << " Pa)\n";
      ApplyNucleationTotal(Q2, map, st2, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_mid);
      TEST_NEAR(qp_sxy(Q2), -(TPV102Params::tau_ini + dtau_mid),
                1e-6 * std::abs(TPV102Params::tau_ini + dtau_mid),
                "R-I06-003 mid-ramp: shape-weighted QP amplitude tracks "
                "partial dtau value");
   }

   // Check 5: wave preservation (R-I06-004 companion).  Add a known
   // wave perturbation to a NON-nearest DOF and verify the QP sum
   // shifts by shape_rest · wave_perturb — ApplyNucleationTotal does
   // NOT clobber it.
   {
      Vector Q3;
      InitializeStateTotal(Q3, ndof_total,
                           TPV102Params::sigma_n, TPV102Params::tau_ini);
      FaultQPNucleationState st3;
      st3.Reset(1);
      // Saturate nucleation.
      ApplyNucleationTotal(Q3, map, st3, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_saturated);
      // Simulate wave landing on DOF 1.
      const real_t wave_dof1 = 1.0e5;  // 100 kPa
      Q3[SXY * ndof_total + 1] += wave_dof1;
      const real_t qp_with_wave = qp_sxy(Q3);
      const real_t expected_qp = -(TPV102Params::tau_ini + dtau)
                               + shape_rest * wave_dof1;
      TEST_NEAR(qp_with_wave, expected_qp,
                1e-6 * std::abs(expected_qp),
                "R-I06-003 + R-I06-004: wave on non-nearest DOF "
                "shifts QP sum by shape_rest · wave (preserved across "
                "ApplyNucleationTotal)");
      // Re-apply at same t: still delta = 0, and the wave must stay.
      ApplyNucleationTotal(Q3, map, st3, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_saturated);
      TEST_NEAR(qp_sxy(Q3), qp_with_wave,
                1e-9 * std::abs(qp_with_wave),
                "R-I06-004: add-delta preserves wave perturbation "
                "across idempotent ApplyNucleationTotal");
   }

   // Check 6: per-stage semantics (R-001 round 3).  `update_state=false`
   // applies the delta to Q but does NOT advance state.dtau_applied,
   // enabling per-RK4-stage nucleation injection on transient Q_stage
   // vectors without corrupting the step-base state.
   {
      std::cout << "\n-- Check 6: per-stage (update_state=false) --\n";
      Vector Q6;
      InitializeStateTotal(Q6, ndof_total,
                           TPV102Params::sigma_n, TPV102Params::tau_ini);
      FaultQPNucleationState st6;
      st6.Reset(1);

      // "Stage" call with update_state=false: Q6 is modified, state is not.
      const real_t t_stage = TPV102Params::nuc_T * 0.4;  // partial ramp
      const real_t dtau_stage = NucleationPerturbation(
         fault_coords[0][0], std::abs(fault_coords[0][2]), t_stage);
      ApplyNucleationTotal(Q6, map, st6, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_stage,
                           /*update_state=*/false);
      TEST_NEAR(qp_sxy(Q6), -(TPV102Params::tau_ini + dtau_stage),
                1e-6 * std::abs(TPV102Params::tau_ini + dtau_stage),
                "R-001: update_state=false applies delta to Q correctly");
      TEST_NEAR(st6.dtau_applied_plus[0], 0.0, 1e-12,
                "R-001: update_state=false leaves state.dtau_applied at 0");

      // Calling again with update_state=false at a LATER t: the delta
      // is (dtau_new - 0) again (state didn't advance), so Q6's nearest
      // DOF gets the FULL new target applied on top of the previous
      // stage-delta.  This is INTENDED: per-stage calls are transient;
      // the driver manages state explicitly at post-step.
      const real_t t_stage2 = TPV102Params::nuc_T * 0.8;
      const real_t dtau_stage2 = NucleationPerturbation(
         fault_coords[0][0], std::abs(fault_coords[0][2]), t_stage2);
      ApplyNucleationTotal(Q6, map, st6, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_stage2,
                           /*update_state=*/false);
      // Second call adds -(dtau_stage2 - 0)/shape_max to Q6[nearest].
      // Q6 already has -(dtau_stage)/shape_max added from the first
      // call, so the nearest DOF is -tau_ini - (dtau_stage + dtau_stage2)/shape_max.
      // The shape-weighted QP sum is therefore -tau_ini - dtau_stage - dtau_stage2.
      TEST_NEAR(qp_sxy(Q6),
                -(TPV102Params::tau_ini + dtau_stage + dtau_stage2),
                1e-6 * (TPV102Params::tau_ini + dtau_stage + dtau_stage2),
                "R-001: successive update_state=false calls compose "
                "(driver keeps state consistent at post-step)");

      // Final update_state=true at t_stage2 finalizes the state.
      // Expected behavior: Q6 gets the incremental delta to reach
      // dtau_stage2 relative to state.dtau_applied=0 → -(dtau_stage2)/shape_max added.
      ApplyNucleationTotal(Q6, map, st6, fault_coords,
                           ndof_per_el, ndof_total,
                           TPV102Params::tau_ini, t_stage2,
                           /*update_state=*/true);
      TEST_NEAR(st6.dtau_applied_plus[0], dtau_stage2,
                1e-12 * dtau_stage2,
                "R-001: update_state=true advances state.dtau_applied");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
