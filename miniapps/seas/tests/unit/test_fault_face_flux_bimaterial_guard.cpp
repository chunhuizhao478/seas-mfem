// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.0.0 REVIEW R-007 — bimaterial guard regression test.
//
// Verifies that FaultFaceFlux::Evaluate triggers MFEM_VERIFY when
// Zp_plus != Zp_minus (or Zs_plus != Zs_minus) beyond the 1e-12
// relative tolerance introduced by R-008 — i.e., on a genuinely
// bimaterial fault face, not on two numerically-identical copies
// produced by different expression chains.
//
// This is a death test with an INVERTED exit-code convention (R-201):
// the parent shell target treats any non-zero rc as PASS (abort ⇒ rc≠0
// via SIGABRT=134), and the only way for the child to signal FAIL is
// to return rc=0 on the fall-through path.  Do NOT change the
// fall-through return code to non-zero: the shell cannot then
// distinguish "guard fired" (abort) from "guard gone" (clean return),
// which is exactly the regression this test exists to catch.

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../config/tpv102_params.hpp"

#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

int main()
{
   DOFData d;
   // Bimaterial mismatch far exceeding the 1e-12 tolerance — a 1%
   // Zp imbalance should always trip the guard.
   d.Zp_plus  = TPV102Params::Zp * 1.01;
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
   real_t arg = d.tau2_0 / (d.sigma_n0 * d.a);
   d.psi      = d.a * std::log(2.0 * TPV102Params::V0 / 1.0 *
                               std::sinh(arg));
   d.slip_rate = 0.0;

   real_t Q_plus[NUM_STATE] = {};
   real_t Q_minus[NUM_STATE] = {};
   // R-205: zero-init the output buffers so any future fall-through
   // diagnostic (e.g. a printf of Q_imp_plus[VX]) reads deterministic
   // zeros, not uninitialized stack memory.
   real_t Q_imp_plus[NUM_STATE] = {};
   real_t Q_imp_minus[NUM_STATE] = {};
   FaultFaceFlux ff(TPV102Params::rho,
                    TPV102Params::cp,
                    TPV102Params::cs);

   // Expect MFEM_VERIFY to abort here.  If the call returns at all,
   // the bimaterial guard has been silently removed or weakened —
   // that is the regression this test catches.
   ff.Evaluate(d, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

   // REGRESSION PATH: if execution reaches this line, the guard is
   // gone.  Return 0 so the shell's `[ rc -ne 0 ]` flips to FAIL.  A
   // non-zero exit from abort() is the PASS signal; a clean zero exit
   // here is the UNIQUE marker that Evaluate returned normally.
   std::cerr << "REGRESSION: Evaluate returned despite bimaterial "
             << "DOFData; the R-F08 / R-008 guard is no longer firing.\n";
   return 0;
}
