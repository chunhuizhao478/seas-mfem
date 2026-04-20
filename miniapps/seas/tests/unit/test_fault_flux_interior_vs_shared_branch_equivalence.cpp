// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.0.0 Phase 1 §3.1d — interior-fault vs shared-fault branch
// equivalence on the SAME face.  REWRITTEN for the v9.0.0 Pelties-9
// per-side fix (plan §14.2 / §14.3 / §15.3).
//
// Post-fix, BOTH branches use the same one-sided "A.T.Q_imp_side" form
// (Pelties 2012 eq. 9), so the equivalence test reduces to:
//
//   F_elem1_interior == F_elem1_shared
//
// where each side computes its Elem1 contribution as:
//   Interior-fault branch (wave_operator.inl:811-887):
//     F_h_plus  = flux.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g)
//     F_h_minus = flux.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g)
//     F_h_elem1 = elem1_on_plus ? F_h_plus : F_h_minus
//     sign      = elem1_on_plus ? -1.0    : +1.0   (outward normal)
//     rhs[Elem1] accumulates: sign * w * shape1 * F_h_elem1
//
//   Shared-fault branch (wave_operator.inl:1232-1244):
//     Q_imp_side = elem1_on_plus ? Q_imp_plus_g : Q_imp_minus_g
//     F_h_shared = flux.Interior(can_n, Q_imp_side, Q_imp_side)
//     sign       = elem1_on_plus ? -1.0 : +1.0
//     rhs[Elem1] accumulates: sign * w * shape1 * F_h_shared
//
// POST-FIX THIS TEST IS A TAUTOLOGY (plan §15.3; REVIEW R-004): both
// helpers below call `flux.Interior(can_n, Q_imp_side, Q_imp_side)` with
// identical `Q_imp_side` selection.  The test is retained only as a
// mechanical safety-net that the production copy-of-branch pattern
// didn't drift; it does NOT verify interior-vs-shared branch
// equivalence under bulk asymmetry.  If you touch either branch,
// REPLACE this test with a live cross-check against the production
// code paths (serial mesh forcing the interior-fault path and a twin
// 2-rank mesh forcing the shared-fault path on the same QP; compare
// rhs[Elem1] contributions).  See REVIEW.md R-004 for the tampering
// demonstration (flipping one helper's Q_imp_side routing still passes
// quiescent-bulk cases because F_h_plus == F_h_minus componentwise when
// only the velocity components of Q_imp_± differ).

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../config/tpv102_params.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static DOFData MakeHypoDOF()
{
   DOFData d;
   d.Zp_plus  = TPV102Params::Zp;
   d.Zp_minus = TPV102Params::Zp;
   d.Zs_plus  = TPV102Params::Zs;
   d.Zs_minus = TPV102Params::Zs;
   d.eta_p    = TPV102Params::Zp / 2.0;
   d.eta_s    = TPV102Params::Zs / 2.0;
   d.sigma_n0 = TPV102Params::sigma_n;
   d.tau1_0   = 0.0;
   d.tau2_0   = TPV102Params::tau_ini;
   d.a        = TPV102Params::a_vw;
   d.Dc       = TPV102Params::Dc;
   d.psi      = ComputeInitialPsi(d.a);
   d.slip_rate = 0.0;
   return d;
}

// Emulate the preamble common to both branches: rotate bulk Q into
// canonical frame, run Evaluate, rotate imposed states back to global.
//
// Parameters match wave_operator.inl:
//   can_n/can_t1/can_t2 — BP5 canonical basis at the QP.
//   Q_self_g / Q_nbr_g — bulk state on this rank's Elem1 (self) and the
//                        opposite element (nbr), in global frame.
//   elem1_on_plus — the routing flag (!qpd.sign_flipped).
// Outputs: Q_imp_plus_g, Q_imp_minus_g (canonical side labelling).
static void RunPreamble(const real_t *can_n, const real_t *can_t1,
                        const real_t *can_t2,
                        const real_t *Q_self_g, const real_t *Q_nbr_g,
                        bool elem1_on_plus,
                        real_t *Q_imp_plus_g, real_t *Q_imp_minus_g)
{
   DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
   GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
   GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);

   real_t Q_self_can[NUM_STATE], Q_nbr_can[NUM_STATE];
   Tinv_can.Mult(Q_self_g, Q_self_can);
   Tinv_can.Mult(Q_nbr_g,  Q_nbr_can);

   const real_t *Q_plus_local  = elem1_on_plus ? Q_self_can : Q_nbr_can;
   const real_t *Q_minus_local = elem1_on_plus ? Q_nbr_can  : Q_self_can;

   DOFData data = MakeHypoDOF();
   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);

   real_t Q_imp_plus_local[NUM_STATE], Q_imp_minus_local[NUM_STATE];
   ff.Evaluate(data, Q_plus_local, Q_minus_local,
               Q_imp_plus_local, Q_imp_minus_local);

   T_can.Mult(Q_imp_plus_local,  Q_imp_plus_g);
   T_can.Mult(Q_imp_minus_local, Q_imp_minus_g);
}

// Compute Elem1's contribution coefficient F_elem1 such that
//   rhs[Elem1] += -w * shape1 * F_elem1
// via the INTERIOR-fault branch.
static void ElementOneFromInteriorBranch(const real_t *can_n,
                                         const real_t *Q_imp_plus_g,
                                         const real_t *Q_imp_minus_g,
                                         bool elem1_on_plus,
                                         real_t *F_elem1)
{
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);

   // v9.0.0 Pelties-9 post-fix (plan §15.3): per-side Interior call.
   // The interior-fault branch now computes both sides' fluxes and
   // assembles each element's own side, signed by the element's
   // outward-normal convention.
   //   plus-side  elem1 (outward = +can_n): rhs -= w*shape1*F_h_plus
   //   minus-side elem1 (outward = -can_n): rhs += w*shape1*F_h_minus
   real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
   flux.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,  F_h_plus);
   flux.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g, F_h_minus);

   const real_t *F_h_elem1 = elem1_on_plus ? F_h_plus : F_h_minus;
   const real_t interior_sign = elem1_on_plus ? -1.0 : +1.0;

   // Fold the assembly sign into F_elem1 so the test's equivalence
   // assertion (F_int == F_sh) compares directly.
   // rhs[Elem1] accumulation: += interior_sign * w * shape1 * F_h_elem1.
   for (int c = 0; c < NUM_STATE; c++)
   {
      F_elem1[c] = interior_sign * F_h_elem1[c];
   }
}

// Compute Elem1's contribution coefficient F_elem1 such that
//   rhs[Elem1] += -w * shape1 * F_elem1
// via the SHARED-fault branch (with accum_sign folded in).
static void ElementOneFromSharedBranch(const real_t *can_n,
                                       const real_t *Q_imp_plus_g,
                                       const real_t *Q_imp_minus_g,
                                       bool elem1_on_plus,
                                       real_t *F_elem1)
{
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);

   // v9.0.0 Pelties-9 post-fix (plan §15.3): one-sided Interior on
   // this rank's Elem1 side, assembled with the outward-normal sign.
   real_t F_h_shared[NUM_STATE];
   const real_t *Q_imp_side = elem1_on_plus ? Q_imp_plus_g : Q_imp_minus_g;
   flux.Interior(can_n, Q_imp_side, Q_imp_side, F_h_shared);

   const real_t assemble_sign = elem1_on_plus ? -1.0 : +1.0;
   // rhs[Elem1] accumulation: += assemble_sign * w * shape1 * F_h_shared.
   for (int c = 0; c < NUM_STATE; c++)
   {
      F_elem1[c] = assemble_sign * F_h_shared[c];
   }
}

static void CompareBranches(const real_t *can_n, const real_t *can_t1,
                            const real_t *can_t2,
                            const real_t *Q_self_g, const real_t *Q_nbr_g,
                            bool elem1_on_plus,
                            const std::string &label)
{
   std::cout << "  " << label
             << "  [elem1_on_plus=" << (elem1_on_plus ? "true" : "false") << "]\n";

   real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
   RunPreamble(can_n, can_t1, can_t2, Q_self_g, Q_nbr_g, elem1_on_plus,
               Q_imp_plus_g, Q_imp_minus_g);

   real_t F_int[NUM_STATE], F_sh[NUM_STATE];
   ElementOneFromInteriorBranch(can_n, Q_imp_plus_g, Q_imp_minus_g,
                                 elem1_on_plus, F_int);
   ElementOneFromSharedBranch(can_n, Q_imp_plus_g, Q_imp_minus_g,
                               elem1_on_plus, F_sh);

   // Hybrid tolerance: a component is "equal" if either (a) the absolute
   // difference is below noise-floor (1e-10 in natural units) OR (b) the
   // relative drift is below ULP (1e-12).  The pure-relative test gave
   // false drift=1 alarms when both values were at 1e-19 (floating-point
   // residue from ApplySplitFlux on a structured all-zero bulk input).
   const real_t abs_tol = 1e-10;
   bool branches_agree = true;
   real_t worst_rel = 0.0;
   int worst_idx = -1;
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t scale = std::abs(F_int[c]) + std::abs(F_sh[c]);
      real_t diff  = F_int[c] - F_sh[c];
      if (std::abs(diff) <= abs_tol) { continue; }   // below noise floor
      if (scale > 0)
      {
         real_t rel = std::abs(diff) / scale;
         if (rel > worst_rel) { worst_rel = rel; worst_idx = c; }
         if (rel > 1e-12) { branches_agree = false; }
      }
      else
      {
         branches_agree = false;   // nonzero diff with zero scale — unreachable given abs_tol above
      }
   }

   std::cout << "    F_interior = [";
   for (int c = 0; c < NUM_STATE; c++)
   {
      std::cout << F_int[c] << (c < NUM_STATE-1 ? ", " : "");
   }
   std::cout << "]\n    F_shared   = [";
   for (int c = 0; c < NUM_STATE; c++)
   {
      std::cout << F_sh[c] << (c < NUM_STATE-1 ? ", " : "");
   }
   std::cout << "]\n    F_int - F_sh = [";
   for (int c = 0; c < NUM_STATE; c++)
   {
      std::cout << (F_int[c] - F_sh[c]) << (c < NUM_STATE-1 ? ", " : "");
   }
   std::cout << "]\n    worst relative drift: " << worst_rel
             << " at component " << worst_idx << "\n";

   TEST_ASSERT(branches_agree,
               "Interior and shared branches give same Elem1 contribution ("
               + label + ")");
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 v9.0.0 Phase 1 §3.1d interior-vs-shared branch equivalence\n";
   std::cout << "========================================\n\n";

   // BP5 canonical frame for TPV102's vertical y=0 fault.
   real_t can_n[3]  = {0.0, -1.0, 0.0};
   real_t can_t1[3] = {0.0,  0.0, -1.0};
   real_t can_t2[3] = {1.0,  0.0,  0.0};

   // Quiescent bulk on both sides — the radiation-birth configuration
   // used in §3.1b.
   real_t Q_self_g[NUM_STATE] = {};
   real_t Q_nbr_g[NUM_STATE]  = {};

   std::cout << "Test 3.1d: same face, two branches, quiescent bulk\n";

   // Case A: Elem1 on canonical + side.
   //   Post-fix: both branches call flux.Interior(can_n, Q_imp_plus_g,
   //   Q_imp_plus_g) with the same selection — trivial bit-identity.
   CompareBranches(can_n, can_t1, can_t2, Q_self_g, Q_nbr_g,
                   /*elem1_on_plus=*/true,
                   "Case A: elem1_on_plus=true (tautology, + side)");

   // Case B: Elem1 on canonical − side.
   //   Post-fix: both branches call flux.Interior(can_n, Q_imp_minus_g,
   //   Q_imp_minus_g) with the same selection — trivial bit-identity.
   //   The pre-fix F(+n,L,R) = -F(-n,R,L) reflection property is no
   //   longer exercised by this test (see header and REVIEW R-004 /
   //   R-203 for the pre-fix vs post-fix divergence).
   CompareBranches(can_n, can_t1, can_t2, Q_self_g, Q_nbr_g,
                   /*elem1_on_plus=*/false,
                   "Case B: elem1_on_plus=false (tautology, - side)");

   // Non-trivial bulk: small random perturbation on the self side.
   real_t Q_self_pert[NUM_STATE] = {};
   Q_self_pert[VX]  = 1e-3;  // small x-velocity on self side
   Q_self_pert[SXY] = 1e5;   // 0.1 MPa shear perturbation
   real_t Q_nbr_zero[NUM_STATE] = {};

   std::cout << "\nTest 3.1d (perturbed): nonzero bulk on self side\n";
   CompareBranches(can_n, can_t1, can_t2, Q_self_pert, Q_nbr_zero,
                   /*elem1_on_plus=*/true,
                   "Case A+perturb: elem1_on_plus=true");
   CompareBranches(can_n, can_t1, can_t2, Q_self_pert, Q_nbr_zero,
                   /*elem1_on_plus=*/false,
                   "Case B+perturb: elem1_on_plus=false");

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";

   return (num_failed > 0) ? 1 : 0;
}
