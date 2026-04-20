// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.0.0 Phase 1 §3.1f — per-side DG assembly for the Pelties-9
// fix.
//
// Reproduces the §14.2 per-side DG assembly with a minimal in-test rhs
// scalar buffer and asserts:
//   (a) both Elem1 and Elem2 receive non-zero contributions (no
//       welded-face cancellation),
//   (b) the plus-side and minus-side rhs are NOT +-related (which was
//       the pre-fix welded-Riemann assumption),
//   (c) the signs match the outward-normal convention: for right-
//       lateral strike-slip (V2 > 0, F_h_plus[VX] > 0 analytically):
//          plus-side  elem1: rhs -= w*shape*F_h_plus  => rhs[VX] < 0
//          minus-side elem1: rhs += w*shape*F_h_minus => rhs[VX] > 0
//
// See plan §15.5 and §18.2 for the full sign derivation and the
// round-3 REVIEW R-F01 history.

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

#define T_OK(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

// Post-breakaway hypocentre fixture matching §3.1b.
static DOFData MakeHypoDOF(real_t V_operating = 1.0,
                           real_t tau_operating = 81e6)
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
   d.tau2_0   = tau_operating;
   d.a        = TPV102Params::a_vw;
   d.Dc       = TPV102Params::Dc;
   real_t arg = tau_operating / (d.sigma_n0 * d.a);
   d.psi = d.a * std::log(2.0 * TPV102Params::V0 / V_operating *
                          std::sinh(arg));
   d.slip_rate = 0.0;
   return d;
}

static void RunCase(bool elem1_on_plus)
{
   const char *tag = elem1_on_plus ? "elem1_on_plus=true"
                                   : "elem1_on_plus=false";
   std::cout << "Case: " << tag << "\n";

   // BP5 canonical frame for TPV102 (vertical y=0 fault).
   real_t can_n[3]  = {0.0, -1.0, 0.0};
   real_t can_t1[3] = {0.0,  0.0, -1.0};  // dip, down
   real_t can_t2[3] = {1.0,  0.0,  0.0};  // strike
   DenseMatrix T_can(NUM_STATE);
   GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);

   DOFData d = MakeHypoDOF();

   // Quiescent bulk on both sides — the radiation-birth configuration.
   real_t Q_plus_local[NUM_STATE]  = {};
   real_t Q_minus_local[NUM_STATE] = {};
   real_t Q_imp_plus_local[NUM_STATE], Q_imp_minus_local[NUM_STATE];
   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp,
                    TPV102Params::cs);
   ff.Evaluate(d, Q_plus_local, Q_minus_local,
               Q_imp_plus_local, Q_imp_minus_local);

   // Rotate imposed states to the global frame (matches wave_operator.inl
   // §14.2 preamble).
   real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
   T_can.Mult(Q_imp_plus_local,  Q_imp_plus_g);
   T_can.Mult(Q_imp_minus_local, Q_imp_minus_g);

   // Per-side Pelties-9 fluxes.
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu,
                    TPV102Params::rho);
   real_t F_plus[NUM_STATE], F_minus[NUM_STATE];
   flux.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,  F_plus);
   flux.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g, F_minus);

   std::cout << "  F_plus[VX]  = " << F_plus[VX]
             << "  F_minus[VX]  = " << F_minus[VX] << "\n";
   std::cout << "  F_plus[SXY] = " << F_plus[SXY]
             << "  F_minus[SXY] = " << F_minus[SXY] << "\n";

   // Replica of the §14.2 assembly (one DOF per element, w=1, shape=1).
   real_t rhs_elem1[NUM_STATE] = {};
   real_t rhs_elem2[NUM_STATE] = {};
   const real_t w = 1.0, shape1 = 1.0, shape2 = 1.0;
   if (elem1_on_plus)
   {
      for (int c = 0; c < NUM_STATE; c++)
      {
         rhs_elem1[c] -= w * shape1 * F_plus[c];
         rhs_elem2[c] += w * shape2 * F_minus[c];
      }
   }
   else
   {
      for (int c = 0; c < NUM_STATE; c++)
      {
         rhs_elem1[c] += w * shape1 * F_minus[c];
         rhs_elem2[c] -= w * shape2 * F_plus[c];
      }
   }

   std::cout << "  rhs_elem1[VX]=" << rhs_elem1[VX]
             << "  rhs_elem2[VX]=" << rhs_elem2[VX] << "\n";

   // (a) non-zero contributions.
   T_OK(std::abs(rhs_elem1[VX]) > 10.0,
        std::string("|rhs_Elem1[VX]| > 10  ") + tag);
   T_OK(std::abs(rhs_elem2[VX]) > 10.0,
        std::string("|rhs_Elem2[VX]| > 10  ") + tag);

   // (b) NOT +-related.  The welded-face assembly would give
   // rhs_elem2 == -rhs_elem1 exactly.  At a fault, SXY is the
   // velocity-driven channel where Q_imp_plus and Q_imp_minus differ
   // most, so their contributions do NOT cancel.
   real_t sum_SXY = rhs_elem1[SXY] + rhs_elem2[SXY];
   T_OK(std::abs(sum_SXY) > 1.0,
        std::string("rhs_Elem1[SXY] + rhs_Elem2[SXY] != 0 (not welded)  ")
        + tag);

   // (c) outward-normal sign convention.
   // SIGN DERIVATION -- important: see plan §18 for history; round-3
   // REVIEW R-F01 proposed inverting these asserts, incorrectly, by
   // applying A_y instead of A_{can_n}.
   //
   // For right-lateral strike-slip TPV102 (V2 > 0, tau2_corr=-eta_s*V2):
   //   can_n = (0, -1, 0)  =>  A_{can_n} = n.A = -A_y
   //   A_y[VX][SXY] = -1/rho   (verified: BuildJacobian(dir=1))
   //   Q_imp_{plus,minus}_g[SXY] = +eta_s*V2 (common; local
   //     tau_{n,t2}=-eta_s*V2 flips sign under the can_n=-y,
   //     can_t2=+x rotation).
   //   F_h_{plus,minus}[VX] = A_{can_n}[VX][SXY] * Q_imp_side_g[SXY]
   //                        = (-A_y[VX][SXY]) * Q_imp_side_g[SXY]
   //                        = -(-1/rho) * (+eta_s*V2)
   //                        = +eta_s*V2/rho = +c_s*V2/2  (POSITIVE)
   //
   // DG rhs assembly (plus side, outward = +can_n):
   //     rhs -= w*shape*F_h_plus  => rhs_plus[VX] < 0.
   // DG rhs assembly (minus side, outward = -can_n):
   //     rhs += w*shape*F_h_minus => rhs_minus[VX] > 0.
   //
   // Common sign-error: computing F_h = A_y*Q directly, omitting the
   // n.A projection.  Do not fall into that trap.
   //
   // R-005 fixture-invariance guard: the sign assertions below assume
   // V2 > 0 (right-lateral strike-slip, tau2_0 > 0).  If a future edit
   // flips tau_operating sign, the derivation still holds but all signs
   // invert; fail loudly on V2 <= 0 so the resulting test failure
   // points at the fixture rather than the sign convention.
   T_OK(d.V2 > 0.0,
        std::string("fixture requires V2 > 0 for the sign-assertion "
                    "block below (tau_operating > 0); see plan §18.2 "
                    "and REVIEW R-005  ") + tag);

   bool plus_side_rhs_sign = elem1_on_plus
                             ? (rhs_elem1[VX] < 0.0)
                             : (rhs_elem2[VX] < 0.0);
   bool minus_side_rhs_sign = elem1_on_plus
                              ? (rhs_elem2[VX] > 0.0)
                              : (rhs_elem1[VX] > 0.0);
   T_OK(plus_side_rhs_sign && minus_side_rhs_sign,
        std::string("rhs_plus[VX]<0 AND rhs_minus[VX]>0 (F_h>0, "
                    "outward-normal)  ") + tag);
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 v9.0.0 Phase 1 §3.1f — per-side assembly\n";
   std::cout << "========================================\n\n";

   RunCase(/*elem1_on_plus=*/true);
   std::cout << "\n";
   RunCase(/*elem1_on_plus=*/false);

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";
   return (num_failed > 0) ? 1 : 0;
}
