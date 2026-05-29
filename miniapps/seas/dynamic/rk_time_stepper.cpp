// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/rk_time_stepper.cpp — Phase 14 implementation.
//
// Non-template pieces of the RK integrator: the Butcher tableau factories
// (classical RK4 + Dormand–Prince DP(4,5)), the tableau validator, and the
// rate-and-state ψ rate.  The templated `AdvanceRKCoupled_Spatial` lives in
// the header (instantiated per-TU for serial `Mesh` and parallel `ParMesh`).

#include "rk_time_stepper.hpp"

#include "../friction/state_evolution.hpp"     // AgingLawPsi
#include "../friction/slip_law_srw_psi.hpp"    // SlipLawSRWPsi

#include <cmath>

namespace mfem
{
namespace seas
{

// =====================================================================
// Tableau factories
// =====================================================================

RKTableau MakeRK4Tableau()
{
   RKTableau t;
   t.name         = "RK4";
   t.stages       = 4;
   t.has_embedded = false;
   t.c    = {0.0, 0.5, 0.5, 1.0};
   t.b    = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
   t.bhat = {};                              // no embedded estimate
   // Strictly-lower-triangular stage matrix (a[i][j] = 0 for j >= i).
   t.a = {
      {0.0, 0.0, 0.0, 0.0},
      {0.5, 0.0, 0.0, 0.0},
      {0.0, 0.5, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0}
   };
   return t;
}

RKTableau MakeDormandPrinceRK45Tableau()
{
   RKTableau t;
   t.name         = "DormandPrinceRK45";
   t.stages       = 7;
   t.has_embedded = true;
   t.c = {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0};

   // 5th-order weight row (primary).  Row 6 of `a` equals `b` ⇒ FSAL.
   t.b = {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0,
          -2187.0 / 6784.0, 11.0 / 84.0, 0.0};

   // 4th-order embedded weight row.
   t.bhat = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0,
             -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};

   // Dormand–Prince DP(4,5) coefficients (verbatim; strict-lower-triangular).
   t.a = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
      {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0,
       -212.0 / 729.0, 0.0, 0.0, 0.0},
      {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0,
       -5103.0 / 18656.0, 0.0, 0.0},
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0,
       -2187.0 / 6784.0, 11.0 / 84.0, 0.0}
   };
   return t;
}

// =====================================================================
// Validator
// =====================================================================

void ValidateTableau(const RKTableau& tab)
{
   const int s = tab.stages;
   MFEM_VERIFY(s >= 1, "ValidateTableau: stages must be >= 1, got " << s);
   MFEM_VERIFY(static_cast<int>(tab.a.size()) == s,
               "ValidateTableau[" << tab.name << "]: a.size() (" << tab.a.size()
               << ") != stages (" << s << ")");
   MFEM_VERIFY(static_cast<int>(tab.b.size()) == s,
               "ValidateTableau[" << tab.name << "]: b.size() (" << tab.b.size()
               << ") != stages (" << s << ")");
   MFEM_VERIFY(static_cast<int>(tab.c.size()) == s,
               "ValidateTableau[" << tab.name << "]: c.size() (" << tab.c.size()
               << ") != stages (" << s << ")");
   MFEM_VERIFY(!tab.has_embedded
               || static_cast<int>(tab.bhat.size()) == s,
               "ValidateTableau[" << tab.name << "]: has_embedded but bhat.size() ("
               << tab.bhat.size() << ") != stages (" << s << ")");

   // Round-off tolerance for the rational-coefficient identities (each
   // fraction carries ~1e-16; a length-7 weighted sum accumulates ~1e-14).
   const real_t tol = 1.0e-10;

   // Σ b = 1.
   real_t sum_b = 0.0;
   for (int i = 0; i < s; ++i) { sum_b += tab.b[i]; }
   MFEM_VERIFY(std::abs(sum_b - 1.0) < tol,
               "ValidateTableau[" << tab.name << "]: Σb = " << sum_b
               << " (must be 1)");

   if (tab.has_embedded)
   {
      real_t sum_bhat = 0.0;
      for (int i = 0; i < s; ++i) { sum_bhat += tab.bhat[i]; }
      MFEM_VERIFY(std::abs(sum_bhat - 1.0) < tol,
                  "ValidateTableau[" << tab.name << "]: Σbhat = " << sum_bhat
                  << " (must be 1)");
   }

   for (int i = 0; i < s; ++i)
   {
      MFEM_VERIFY(static_cast<int>(tab.a[i].size()) == s,
                  "ValidateTableau[" << tab.name << "]: a[" << i << "].size() ("
                  << tab.a[i].size() << ") != stages (" << s << ")");
      // Strict-lower-triangular: a[i][j] = 0 for j >= i.
      for (int j = i; j < s; ++j)
      {
         MFEM_VERIFY(tab.a[i][j] == 0.0,
                     "ValidateTableau[" << tab.name << "]: a[" << i << "][" << j
                     << "] = " << tab.a[i][j] << " (must be 0 for explicit RK; "
                     "strict-lower-triangular)");
      }
      // c_i = Σ_j a_ij.
      real_t row_sum = 0.0;
      for (int j = 0; j < i; ++j) { row_sum += tab.a[i][j]; }
      MFEM_VERIFY(std::abs(row_sum - tab.c[i]) < tol,
                  "ValidateTableau[" << tab.name << "]: c[" << i << "] = "
                  << tab.c[i] << " != Σ_j a[" << i << "][j] = " << row_sum);
   }
}

// =====================================================================
// Rate-and-state ψ rate (plan §14.2)
// =====================================================================

real_t PsiRate(const spatial::RateStateBlock&       rs_cfg,
               const DOFData&                        d,
               real_t                                V,
               const spatial::RateStatePerDOFParams& rs,
               int                                   m)
{
   switch (rs_cfg.state_evolution)
   {
      case spatial::StateEvolutionKind::AgingLaw:
      {
         // ψ-space aging rate dψ/dt = (b·V0/Dc)·(exp((f0−ψ)/b) − V/V0).
         // Global b/V0/f0 from the config block (the production aging substep
         // iterator builds AgingLawPsi the same way:
         // friction_iterator_factory.cpp:90); per-QP ψ/Dc from the DOFData.
         const AgingLawPsi law(rs_cfg.b_default, rs_cfg.V_0_default,
                               rs_cfg.f_0_default);
         return law.Rate(V, d.psi, d.Dc);
      }
      case spatial::StateEvolutionKind::SlipLawStrongRateWeakening:
      {
         // SRW slip-law (TPV104).  Global scalars (b, V0, f0, muW=f_w_default)
         // from the config; per-QP V_w(m)/a(m) from the resolved rs — the same
         // side-channel source as the production iterator
         // (friction_iterator_factory.cpp:73-81).  Production mode forbids the
         // scalar-V_w base virtual (R-001); Rate_SRW takes per-QP V_w + a.
         MFEM_VERIFY(rs.V_w.Size() > m && rs.a.Size() > m,
                     "PsiRate(SRW): rs.V_w / rs.a must be sized to the fault "
                     "DOF count; got V_w.Size()=" << rs.V_w.Size()
                     << ", a.Size()=" << rs.a.Size() << ", m=" << m);
         SlipLawSRWPsi law(rs_cfg.a_default, rs_cfg.b_default,
                           rs_cfg.V_0_default, rs_cfg.f_0_default,
                           rs_cfg.f_w_default, rs_cfg.V_w_default);
         law.SetProductionMode();
         return law.Rate_SRW(V, d.psi, d.Dc, rs.V_w(m), rs.a(m));
      }
   }
   MFEM_ABORT("PsiRate: unhandled state_evolution = "
              << static_cast<int>(rs_cfg.state_evolution));
   return 0.0;
}

} // namespace seas
} // namespace mfem
