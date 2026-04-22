// ADER I-05 Phase 5: equivalence tests for FaultFaceFlux::EvaluateADER
// and EvaluateADERTotal (the time-integrated friction-solve wrappers).
//
// Gates (plan §Phase 5 acceptance):
//   A. Locked fault: EvaluateADER on a pure pre-stress + zero-velocity I±
//      gives `data.sigma_n_corr ≈ sigma_n0` within the friction-solver
//      tolerance (no rupture → tractions match initial conditions).
//   B. dt → 0 consistency: |EvaluateADER(I=dt·Q̄, dt) − dt·Evaluate(Q̄)|
//      scales as O(dt²) over three halvings of dt when the state is
//      smooth.
//   C. EvaluateADERTotal at zero fluctuation: inputs that are pure
//      pre-stress + zero velocity yield `data.sigma_n_corr ≈ sigma_n0`
//      and `V1 = V2 = 0`.
//   D. dt → 0 consistency for EvaluateADERTotal: analogous to Gate B
//      against the total-stress Evaluate.

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

namespace
{
constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;
const real_t kCp  = std::sqrt((kLambda + 2.0 * kMu) / kRho);
const real_t kCs  = std::sqrt(kMu / kRho);
const real_t kZp  = kRho * kCp;
const real_t kZs  = kRho * kCs;

DOFData MakeHomogeneousDOF(real_t sigma_n0 = 120e6, real_t a = 0.004)
{
   DOFData d;
   d.Zp_plus = kZp; d.Zp_minus = kZp;
   d.Zs_plus = kZs; d.Zs_minus = kZs;
   d.eta_p = kZp / 2.0;
   d.eta_s = kZs / 2.0;
   d.sigma_n0 = sigma_n0;
   d.tau1_0 = 0.0;
   d.tau2_0 = 0.0;
   d.a = a;
   d.Dc = 0.4;
   d.psi = 0.5;
   d.slip_rate = 0.0;
   return d;
}

real_t MaxAbs(const real_t v[NUM_STATE])
{
   real_t m = 0;
   for (int c = 0; c < NUM_STATE; c++) { m = std::max(m, std::abs(v[c])); }
   return m;
}

real_t MaxAbsDiff(const real_t a[NUM_STATE], const real_t b[NUM_STATE])
{
   real_t m = 0;
   for (int c = 0; c < NUM_STATE; c++) { m = std::max(m, std::abs(a[c] - b[c])); }
   return m;
}
} // anonymous

int main()
{
   std::cout << "\n=== ADER I-05 Phase 5: EvaluateADER / EvaluateADERTotal "
             << "equivalence tests ===\n";

   FaultFaceFlux ff(kRho, kCp, kCs);

   // ------------------------------------------------------------------
   // Gate A: locked-fault EvaluateADER preserves pre-stress.
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate A: locked-fault EvaluateADER --\n";
      DOFData data = MakeHomogeneousDOF();
      data.tau1_0 = 0.0;
      data.tau2_0 = 20e6;  // strike component under the threshold

      // Fluctuation Q± represents the PERTURBATION away from pre-stress.
      // Locked = zero fluctuation and zero velocity on both sides.
      real_t Q_plus[NUM_STATE]  = {0};
      real_t Q_minus[NUM_STATE] = {0};

      const real_t dt = 1.0e-4;
      real_t I_plus[NUM_STATE], I_minus[NUM_STATE];
      real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++)
      {
         I_plus[c]  = dt * Q_plus[c];
         I_minus[c] = dt * Q_minus[c];
      }
      ff.EvaluateADER(data, I_plus, I_minus, dt, I_imp_plus, I_imp_minus);

      // On a locked fault, tau_corr must equal tau_0 and sigma_n_corr
      // must equal sigma_n0 (no slip → no relief of pre-stress).
      TEST_LE(std::abs(data.sigma_n_corr - data.sigma_n0), 1.0,
              "Gate A: sigma_n_corr = sigma_n0 (locked fault, within 1 Pa)");
      TEST_LE(std::abs(data.tau2_corr  - data.tau2_0),  1.0,
              "Gate A: tau2_corr = tau2_0 (locked fault, within 1 Pa)");
      TEST_LE(std::abs(data.slip_rate), 1e-12,
              "Gate A: |slip_rate| < 1e-12 (locked fault)");
   }

   // ------------------------------------------------------------------
   // Gate B: dt → 0 consistency — EvaluateADER on I = dt·Q̄ matches
   //         dt · (Evaluate's Q_imp on Q̄) to O(dt²) (with the same psi).
   //         We use a smooth non-trivial Q̄ and compare across 3 dt values.
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate B: EvaluateADER dt → 0 consistency --\n";
      real_t Q_avg_plus[NUM_STATE]  = {0};
      real_t Q_avg_minus[NUM_STATE] = {0};
      Q_avg_plus[SXX]  = -5e6; Q_avg_plus[SXY]  = 2e6;
      Q_avg_minus[SXX] = -6e6; Q_avg_minus[SXY] = 1.5e6;
      Q_avg_plus[VY]   =  0.05;  // small slip rate perturbation
      Q_avg_minus[VY]  = -0.05;

      auto make_data = [](){
         DOFData d = MakeHomogeneousDOF(60e6);
         d.tau2_0 = 20e6; d.psi = 0.5;
         return d;
      };

      const real_t dt_list[3] = { 1.0e-3, 5.0e-4, 2.5e-4 };
      real_t diff_list[3] = { 0.0, 0.0, 0.0 };

      for (int idt = 0; idt < 3; idt++)
      {
         const real_t dt = dt_list[idt];

         // Reference: run Evaluate on Q̄ directly → Q_imp_ref.  Save psi
         // before so both paths see the same input state.
         DOFData data_ref = make_data();
         real_t Q_imp_plus_ref[NUM_STATE], Q_imp_minus_ref[NUM_STATE];
         ff.Evaluate(data_ref, Q_avg_plus, Q_avg_minus,
                     Q_imp_plus_ref, Q_imp_minus_ref);

         // ADER: run EvaluateADER on I = dt · Q̄ → I_imp; divide by dt.
         DOFData data_ader = make_data();
         real_t I_plus[NUM_STATE], I_minus[NUM_STATE];
         real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_plus[c]  = dt * Q_avg_plus[c];
            I_minus[c] = dt * Q_avg_minus[c];
         }
         ff.EvaluateADER(data_ader, I_plus, I_minus, dt,
                         I_imp_plus, I_imp_minus);
         real_t Q_imp_plus_ader[NUM_STATE], Q_imp_minus_ader[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_imp_plus_ader[c]  = I_imp_plus[c]  / dt;
            Q_imp_minus_ader[c] = I_imp_minus[c] / dt;
         }

         // EvaluateADER on dt·Q̄ should give identically `dt·Evaluate(Q̄)`
         // (pure linear rescaling; the friction solve is evaluated on Q̄
         // exactly in both paths).  Assert within solver tolerance.
         const real_t diff_p = MaxAbsDiff(Q_imp_plus_ader, Q_imp_plus_ref);
         const real_t diff_m = MaxAbsDiff(Q_imp_minus_ader, Q_imp_minus_ref);
         diff_list[idt] = std::max(diff_p, diff_m);
         const real_t scale = std::max(MaxAbs(Q_imp_plus_ref),
                                       MaxAbs(Q_imp_minus_ref));
         std::cout << "  dt = " << std::scientific << std::setprecision(3)
                   << dt << "  |diff|/|Q_imp| = "
                   << diff_list[idt] / std::max(scale, real_t(1.0)) << "\n";
         TEST_LE(diff_list[idt], 1e-6 * std::max(scale, real_t(1.0)),
                 "Gate B: EvaluateADER on I=dt·Q̄ matches Evaluate on Q̄ (solver tolerance)");
      }
   }

   // ------------------------------------------------------------------
   // Gate C: EvaluateADERTotal on pure pre-stress + zero velocity.
   //         Under the total-Q convention, Q_±_tot encodes TOTAL stress
   //         (pre-stress + fluctuation).  Zero-fluctuation input = pure
   //         pre-stress bulk state; no rupture → sigma_n_corr = sigma_n0
   //         (the TOTAL stored in data.sigma_n0 under Phase 4) and
   //         V1 = V2 = 0.
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate C: EvaluateADERTotal locked-fault --\n";
      DOFData data = MakeHomogeneousDOF(120e6);
      // Phase 4 convention: under total-Q, DOFData pre-stress fields are
      // zeroed by the driver so EvaluateTotal doesn't double-count.  Q_tot
      // carries the pre-stress directly.
      data.sigma_n0 = 0.0;
      data.tau1_0   = 0.0;
      data.tau2_0   = 0.0;

      real_t Q_plus_tot[NUM_STATE]  = {0};
      real_t Q_minus_tot[NUM_STATE] = {0};
      // Total-stress pre-stress on both sides: sigma_nn = sigma_n0,
      // tau_nt2 = tau_ini in fault-local canonical frame.
      Q_plus_tot[SXX]  = 120e6;  Q_plus_tot[SXZ]  = 20e6;  // tau_nt2 = SXZ in fault-local
      Q_minus_tot[SXX] = 120e6;  Q_minus_tot[SXZ] = 20e6;

      const real_t dt = 1.0e-4;
      real_t I_plus[NUM_STATE], I_minus[NUM_STATE];
      real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++)
      {
         I_plus[c]  = dt * Q_plus_tot[c];
         I_minus[c] = dt * Q_minus_tot[c];
      }
      ff.EvaluateADERTotal(data, I_plus, I_minus, dt,
                           I_imp_plus, I_imp_minus);

      TEST_LE(std::abs(data.sigma_n_corr - 120e6), 1.0,
              "Gate C: sigma_n_corr = 120 MPa (locked fault, total-Q)");
      TEST_LE(std::abs(data.slip_rate), 1e-12,
              "Gate C: |slip_rate| < 1e-12 (locked fault)");
   }

   // ------------------------------------------------------------------
   // Gate D: dt → 0 consistency for EvaluateADERTotal — analogous to
   //         Gate B but using EvaluateTotal as the reference.
   // ------------------------------------------------------------------
   {
      std::cout << "\n-- Gate D: EvaluateADERTotal dt → 0 consistency --\n";
      real_t Q_avg_plus[NUM_STATE]  = {0};
      real_t Q_avg_minus[NUM_STATE] = {0};
      // Small fluctuation on top of a pre-stress bulk state.
      Q_avg_plus[SXX]  = 120e6 - 1e6;  Q_avg_plus[SXZ]  = 20e6 + 2e5;
      Q_avg_minus[SXX] = 120e6 + 1e6;  Q_avg_minus[SXZ] = 20e6 - 2e5;
      Q_avg_plus[VZ]   =  0.02;
      Q_avg_minus[VZ]  = -0.02;

      auto make_data = [](){
         DOFData d = MakeHomogeneousDOF(120e6);
         d.sigma_n0 = 0.0; d.tau1_0 = 0.0; d.tau2_0 = 0.0;
         d.psi = 0.5;
         return d;
      };

      const real_t dt_list[3] = { 1.0e-3, 5.0e-4, 2.5e-4 };
      for (int idt = 0; idt < 3; idt++)
      {
         const real_t dt = dt_list[idt];

         DOFData data_ref = make_data();
         real_t Q_imp_plus_ref[NUM_STATE], Q_imp_minus_ref[NUM_STATE];
         ff.EvaluateTotal(data_ref, Q_avg_plus, Q_avg_minus,
                          Q_imp_plus_ref, Q_imp_minus_ref);

         DOFData data_ader = make_data();
         real_t I_plus[NUM_STATE], I_minus[NUM_STATE];
         real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_plus[c]  = dt * Q_avg_plus[c];
            I_minus[c] = dt * Q_avg_minus[c];
         }
         ff.EvaluateADERTotal(data_ader, I_plus, I_minus, dt,
                              I_imp_plus, I_imp_minus);
         real_t Q_imp_plus_ader[NUM_STATE], Q_imp_minus_ader[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_imp_plus_ader[c]  = I_imp_plus[c]  / dt;
            Q_imp_minus_ader[c] = I_imp_minus[c] / dt;
         }

         const real_t diff_p = MaxAbsDiff(Q_imp_plus_ader, Q_imp_plus_ref);
         const real_t diff_m = MaxAbsDiff(Q_imp_minus_ader, Q_imp_minus_ref);
         const real_t diff = std::max(diff_p, diff_m);
         const real_t scale = std::max(MaxAbs(Q_imp_plus_ref),
                                       MaxAbs(Q_imp_minus_ref));
         std::cout << "  dt = " << std::scientific << std::setprecision(3)
                   << dt << "  |diff|/|Q_imp| = "
                   << diff / std::max(scale, real_t(1.0)) << "\n";
         TEST_LE(diff, 1e-6 * std::max(scale, real_t(1.0)),
                 "Gate D: EvaluateADERTotal on I=dt·Q̄ matches EvaluateTotal on Q̄");
      }
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
