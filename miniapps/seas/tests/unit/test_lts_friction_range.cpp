// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_friction_range.cpp — Appendix B.9 of
// PLAN_clustered_lts_ader_2026-07-18.md (LTS Phase 3, A.7).
//
// The LTS fault sweep advances each cluster's cluster-contiguous fault-QP range
// at that cluster's own dt via the range overload
//   FrictionSubStepIterator::Advance(qp_begin, qp_end, ...).
// This test pins the three B.9 properties on a rate-state (aging) fixture:
//
//   (1) per-QP equivalence: advancing QP i as part of a range gives the SAME
//       (psi, slip1, slip2, V, I_imp) as advancing that QP standalone (a
//       single-QP whole-vector Advance) with its own dt/quadrature/Q-trace —
//       to 1e-15 (in fact bit-for-bit: the friction solve is per-QP local);
//   (2) I_imp outside the range is left byte-untouched (a sentinel survives),
//       and a second disjoint range does not perturb the first range's output;
//   (3) the Σ deltaT == dt_step guard fires when dt_step ≠ Σ deltaT.
//
// Two ranges are advanced with DIFFERENT dt (dt_A ≠ dt_B), exactly the
// multi-cluster case the driver drives per sync tick.

#include "mfem.hpp"

#include "../../dynamic/fault_face_flux.hpp"       // DOFData, FaultFaceFlux, NUM_STATE
#include "../../dynamic/friction_solver.hpp"       // FrictionSolver::Method
#include "../../dynamic/friction_substep_iterator.hpp"
#include "../../dynamic/wave_state.hpp"            // QIndex (VZ)
#include "../../friction/state_evolution.hpp"      // AgingLawPsi
#include "../../friction/slip_law_srw_psi.hpp"     // SlipLawSRWPsi (SRW sub-case)

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int g_checks = 0, g_fails = 0;
#define CHECK(cond, msg)                                                     \
   do { ++g_checks; if (!(cond)) { ++g_fails;                               \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg); } } while (0)

namespace
{
constexpr real_t kRho = 2670.0, kCp = 6000.0, kCs = 3464.0;
const real_t kZp = kRho * kCp, kZs = kRho * kCs;

// A rate-state (aging) QP fixture; `pert` differentiates the QPs.
DOFData MakeAgingDOF(int pert)
{
   DOFData d;
   d.Zp_plus = kZp; d.Zp_minus = kZp;
   d.Zs_plus = kZs; d.Zs_minus = kZs;
   d.eta_p = 0.5 * kZp; d.eta_s = 0.5 * kZs;
   d.sigma_n0 = 50.0e6;
   d.tau1_0   = 0.0;
   d.tau2_0   = 29.2e6 + pert * 1.0e6;
   d.a   = 0.010;
   d.b   = 0.015;
   d.Dc  = 2.0;
   d.psi = 0.60 + pert * 0.01;
   return d;
}

// Non-uniform quadrature (Σδt = dt, Σw = 1), so a [0]-vs-[o] mis-index would
// be observable.
void MakeQuadrature(int O, real_t dt, std::vector<real_t> &deltaT,
                    std::vector<real_t> &w)
{
   deltaT.assign(O, 0.0); w.assign(O, 0.0);
   real_t ds = 0.0, ws = 0.0;
   for (int o = 0; o < O; ++o)
   { deltaT[o] = 1.0 + 0.5 * o; w[o] = 1.0 + 0.25 * o; ds += deltaT[o]; ws += w[o]; }
   for (int o = 0; o < O; ++o) { deltaT[o] = deltaT[o] / ds * dt; w[o] = w[o] / ws; }
}

// Per-(o, i) varying pointwise-Q field (a VZ velocity jump driving slip).
void MakeQField(int n, int O, std::vector<std::vector<real_t>> &Qp,
                std::vector<std::vector<real_t>> &Qm)
{
   Qp.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
   Qm.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
   for (int o = 0; o < O; ++o)
   {
      for (int i = 0; i < n; ++i)
      {
         const size_t b = static_cast<size_t>(i) * NUM_STATE;
         Qp[o][b + VZ] = -(1.0 + 0.3 * o + 0.05 * i);
         Qm[o][b + VZ] = +(1.0 + 0.3 * o + 0.05 * i);
      }
   }
}

// Max |diff| over the mutated rate-state DOFData fields.
real_t DofDiff(const DOFData &a, const DOFData &b)
{
   real_t m = 0.0;
   m = std::max(m, std::abs(a.psi        - b.psi));
   m = std::max(m, std::abs(a.slip1      - b.slip1));
   m = std::max(m, std::abs(a.slip2      - b.slip2));
   m = std::max(m, std::abs(a.V1         - b.V1));
   m = std::max(m, std::abs(a.V2         - b.V2));
   m = std::max(m, std::abs(a.slip_rate  - b.slip_rate));
   m = std::max(m, std::abs(a.tau1_corr  - b.tau1_corr));
   m = std::max(m, std::abs(a.tau2_corr  - b.tau2_corr));
   m = std::max(m, std::abs(a.sigma_n_corr - b.sigma_n_corr));
   return m;
}

// Extract QP i's NUM_STATE trace from a per-substep global field into a
// single-QP field (for the standalone reference).
std::vector<std::vector<real_t>> SliceQP(
   const std::vector<std::vector<real_t>> &Q, int i, int O)
{
   std::vector<std::vector<real_t>> out(O, std::vector<real_t>(NUM_STATE, 0.0));
   for (int o = 0; o < O; ++o)
   {
      for (int c = 0; c < NUM_STATE; ++c)
      { out[o][c] = Q[o][static_cast<size_t>(i) * NUM_STATE + c]; }
   }
   return out;
}

const auto noop_nuc = [](real_t, real_t) {};
}  // namespace

int main()
{
   const int N = 5;                 // 5 global fault QPs
   const int O = 3;                 // ADER order (sub-steps)
   const std::size_t k = 3;         // range split: [0,3) at dt_A, [3,5) at dt_B
   const real_t dt_A = 1.0e-4, dt_B = 2.0e-4;   // DIFFERENT per range
   const real_t t0 = 3.7;           // non-zero step-start time

   FaultFaceFlux flux(kRho, kCp, kCs);
   const real_t b = 0.015, V0 = 1.0e-6, f0 = 0.6;
   const auto method = FrictionSolver::Method::Brent;

   std::vector<Vector> coords(N, Vector(3));
   for (int i = 0; i < N; ++i) { coords[i] = 0.0; coords[i](2) = -7500.0; }

   std::vector<std::vector<real_t>> Qp, Qm;
   MakeQField(N, O, Qp, Qm);

   std::vector<real_t> deltaT_A, w_A, deltaT_B, w_B;
   MakeQuadrature(O, dt_A, deltaT_A, w_A);
   MakeQuadrature(O, dt_B, deltaT_B, w_B);

   const size_t flat = static_cast<size_t>(NUM_STATE) * N;
   const real_t SENT = -123456.0;   // finite sentinel for the untouched check

   // ---------------------------------------------------------------------
   // RANGED path: one iterator, one global buffer, two disjoint ranges at
   // different dt.  Range A = [0,k) at dt_A; Range B = [k,N) at dt_B.
   // ---------------------------------------------------------------------
   std::vector<DOFData> dof(N);
   for (int i = 0; i < N; ++i) { dof[i] = MakeAgingDOF(i); }
   std::vector<real_t> Ip(flat, SENT), Im(flat, SENT);

   RateStateAgingIterator it(flux, AgingLawPsi(b, V0, f0), method);

   it.SetSubSteps(deltaT_A, w_A);
   it.Advance(0, k, dof, coords, Qp, Qm, dt_A, t0, Ip.data(), Im.data(),
              noop_nuc);

   // (2a) After range A, the [k,N) slice of I_imp must still hold the sentinel.
   {
      bool untouched = true;
      for (size_t j = k * NUM_STATE; j < flat; ++j)
      { untouched = untouched && (Ip[j] == SENT) && (Im[j] == SENT); }
      CHECK(untouched, "range A wrote I_imp outside [0,k)");
   }
   // Snapshot range A's output slice to prove range B does not perturb it.
   std::vector<real_t> Ip_A(Ip.begin(), Ip.begin() + k * NUM_STATE);
   std::vector<real_t> Im_A(Im.begin(), Im.begin() + k * NUM_STATE);
   std::vector<DOFData> dof_A(dof.begin(), dof.begin() + k);

   it.SetSubSteps(deltaT_B, w_B);
   it.Advance(k, N, dof, coords, Qp, Qm, dt_B, t0, Ip.data(), Im.data(),
              noop_nuc);

   // (2b) Range B must not perturb range A's DOFData or I_imp slice.
   {
      bool a_stable = true;
      for (size_t j = 0; j < k * NUM_STATE; ++j)
      { a_stable = a_stable && (Ip[j] == Ip_A[j]) && (Im[j] == Im_A[j]); }
      CHECK(a_stable, "range B perturbed range A's I_imp slice");
      real_t md = 0.0;
      for (size_t i = 0; i < k; ++i) { md = std::max(md, DofDiff(dof[i], dof_A[i])); }
      CHECK(md == 0.0, "range B perturbed range A's DOFData");
   }

   // ---------------------------------------------------------------------
   // (1) REFERENCE: each QP advanced standalone with its own dt/quadrature.
   // ---------------------------------------------------------------------
   real_t maxdiff_dof = 0.0, maxdiff_I = 0.0;
   for (int i = 0; i < N; ++i)
   {
      const bool in_A = (static_cast<std::size_t>(i) < k);
      const real_t dt_i = in_A ? dt_A : dt_B;
      const std::vector<real_t> &dT = in_A ? deltaT_A : deltaT_B;
      const std::vector<real_t> &wt = in_A ? w_A : w_B;

      std::vector<DOFData> d1(1, MakeAgingDOF(i));
      std::vector<Vector> c1(1, coords[i]);
      std::vector<std::vector<real_t>> qp1 = SliceQP(Qp, i, O);
      std::vector<std::vector<real_t>> qm1 = SliceQP(Qm, i, O);
      std::vector<real_t> ip1(NUM_STATE, 0.0), im1(NUM_STATE, 0.0);

      RateStateAgingIterator ref(flux, AgingLawPsi(b, V0, f0), method);
      ref.SetSubSteps(dT, wt);
      ref.Advance(d1, c1, qp1, qm1, dt_i, t0, ip1.data(), im1.data(), noop_nuc);

      maxdiff_dof = std::max(maxdiff_dof, DofDiff(dof[i], d1[0]));
      for (int c = 0; c < NUM_STATE; ++c)
      {
         maxdiff_I = std::max(maxdiff_I,
                              std::abs(Ip[static_cast<size_t>(i) * NUM_STATE + c] - ip1[c]));
         maxdiff_I = std::max(maxdiff_I,
                              std::abs(Im[static_cast<size_t>(i) * NUM_STATE + c] - im1[c]));
      }
   }
   // R-008: the per-QP friction solve is identical-order arithmetic, so the
   // ranged and standalone results are BIT-for-bit equal (not merely ~1e-15).
   CHECK(maxdiff_dof == 0.0, "ranged per-QP DOFData != standalone reference (bit-exact)");
   CHECK(maxdiff_I   == 0.0, "ranged per-QP I_imp != standalone reference (bit-exact)");

   // Non-triviality: the fixture must actually solve (psi moved or slip built).
   CHECK(dof[0].psi != MakeAgingDOF(0).psi || std::abs(dof[0].slip2) > 0.0,
         "fixture is vacuous (no friction/state solve happened)");

   // ---------------------------------------------------------------------
   // R-001: SRW sub-case.  The one range-body operation that depends on the
   // GLOBAL fault-QP index (not `iu - qp_begin`) is the SRW ψ-update's
   // V_w side-channel `(*Vw)(i)`.  A per-QP-VARYING V_w makes a local-vs-global
   // mis-index observable: range B's first QP would read V_w[0] instead of
   // V_w[k].  The aging path above cannot catch it (aging ignores the index).
   // ---------------------------------------------------------------------
   {
      const real_t a0 = 0.008, bb = 0.012, V0s = 1.0e-6, f0s = 0.6;
      const real_t muW = 0.1, Vwd = 0.1;
      const auto srw_method = FrictionSolver::Method::NewtonRaphsonStable;

      mfem::Vector Vw(N);
      for (int i = 0; i < N; ++i) { Vw(i) = 0.1 + 0.02 * i; }   // varies per QP

      auto make_srw_dof = [&](int i)
      {
         DOFData d;
         d.Zp_plus = kZp; d.Zp_minus = kZp;
         d.Zs_plus = kZs; d.Zs_minus = kZs;
         d.eta_p = 0.5 * kZp; d.eta_s = 0.5 * kZs;
         d.sigma_n0 = 50.0e6; d.tau1_0 = 0.0;
         d.tau2_0 = 29.38e6 + i * 1.0e6;
         d.a = 0.008 + i * 1.0e-4; d.b = 0.012; d.Dc = 0.40;
         d.psi = 0.564 + i * 0.01;
         return d;
      };

      // Ranged: one SRW iterator, GLOBAL V_w, two ranges at different dt.
      std::vector<DOFData> sdof(N);
      for (int i = 0; i < N; ++i) { sdof[i] = make_srw_dof(i); }
      std::vector<real_t> sIp(flat, 0.0), sIm(flat, 0.0);
      SlipLawSRWPsi slaw(a0, bb, V0s, f0s, muW, Vwd); slaw.SetProductionMode();
      RateStateSlipLawSrwIterator sit(flux, std::move(slaw), srw_method, &Vw);
      sit.SetSubSteps(deltaT_A, w_A);
      sit.Advance(0, k, sdof, coords, Qp, Qm, dt_A, t0, sIp.data(), sIm.data(),
                  noop_nuc);
      sit.SetSubSteps(deltaT_B, w_B);
      sit.Advance(k, N, sdof, coords, Qp, Qm, dt_B, t0, sIp.data(), sIm.data(),
                  noop_nuc);

      // Reference: per-QP standalone with a SINGLE-element V_w = {Vw(i)} — so a
      // global-index bug (reading V_w[local]) diverges from this reference.
      real_t sdiff = 0.0;
      for (int i = 0; i < N; ++i)
      {
         const bool inA = (static_cast<std::size_t>(i) < k);
         const real_t dt_i = inA ? dt_A : dt_B;
         const std::vector<real_t> &dT = inA ? deltaT_A : deltaT_B;
         const std::vector<real_t> &wt = inA ? w_A : w_B;
         std::vector<DOFData> d1(1, make_srw_dof(i));
         std::vector<Vector> c1(1, coords[i]);
         std::vector<std::vector<real_t>> qp1 = SliceQP(Qp, i, O);
         std::vector<std::vector<real_t>> qm1 = SliceQP(Qm, i, O);
         std::vector<real_t> ip1(NUM_STATE, 0.0), im1(NUM_STATE, 0.0);
         mfem::Vector Vw1(1); Vw1(0) = Vw(i);
         SlipLawSRWPsi slaw1(a0, bb, V0s, f0s, muW, Vwd); slaw1.SetProductionMode();
         RateStateSlipLawSrwIterator ref1(flux, std::move(slaw1), srw_method, &Vw1);
         ref1.SetSubSteps(dT, wt);
         ref1.Advance(d1, c1, qp1, qm1, dt_i, t0, ip1.data(), im1.data(), noop_nuc);
         sdiff = std::max(sdiff, DofDiff(sdof[i], d1[0]));
      }
      CHECK(sdiff == 0.0,
            "SRW ranged per-QP != standalone reference (global V_w index)");
      CHECK(sdof[0].psi != make_srw_dof(0).psi || std::abs(sdof[0].slip2) > 0.0,
            "SRW fixture is vacuous");
   }

   // ---------------------------------------------------------------------
   // (3) Σ deltaT == dt_step guard fires on mismatch.
   // ---------------------------------------------------------------------
   {
      std::vector<DOFData> d = dof;
      std::vector<real_t> ip(flat, 0.0), im(flat, 0.0);
      it.SetSubSteps(deltaT_A, w_A);          // Σ deltaT_A == dt_A
      bool threw = false;
      try
      {
         // dt_step deliberately WRONG (dt_B ≠ Σ deltaT_A).
         it.Advance(0, k, d, coords, Qp, Qm, dt_B, t0, ip.data(), im.data(),
                    noop_nuc);
      }
      catch (const std::exception &) { threw = true; }
      CHECK(threw, "Σ deltaT == dt_step guard did not fire on mismatch");
   }

   // ---------------------------------------------------------------------
   // (3b) Range-bounds guard: qp_end > n throws.
   // ---------------------------------------------------------------------
   {
      std::vector<DOFData> d = dof;
      std::vector<real_t> ip(flat, 0.0), im(flat, 0.0);
      it.SetSubSteps(deltaT_A, w_A);
      bool threw = false;
      try
      {
         it.Advance(0, N + 1, d, coords, Qp, Qm, dt_A, t0, ip.data(), im.data(),
                    noop_nuc);
      }
      catch (const std::exception &) { threw = true; }
      CHECK(threw, "range-bounds guard (qp_end > n) did not fire");
   }

   std::printf("test_lts_friction_range: %d checks, %d failures\n",
               g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}
