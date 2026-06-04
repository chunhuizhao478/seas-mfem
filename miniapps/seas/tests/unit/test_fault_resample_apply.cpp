// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Phase 3 unit test: the rate-state Δψ RESAMPLE apply inside
// RateStateSubStepIterator::Advance (dynamic/friction_substep_iterator.cpp),
// using the Phase-2 projector (dynamic/fault_resample.hpp).  Verifies:
//   T1: with R = I (unisolvent rule, no over-integration) the resample matches
//       no-resample to round-off — i.e. it is inert (the driver additionally
//       DISABLES it there for a true byte-exact no-op).
//   T2: with R != I (over-integrated) the resample applies EXACTLY
//       psi = psi_before + R·Δψ (Δψ = the net macro-step state increment),
//       and it actually changes psi (non-vacuous).
//   T3: the resample acts PER FACE (face 0 is unaffected by face 1's QPs).
//
// Tiny fixture (no mesh); links the iterator + flux + friction solver only.

#include "mfem.hpp"

#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_substep_iterator.hpp"
#include "../../dynamic/fault_resample.hpp"
#include "../../friction/state_evolution.hpp"   // AgingLawPsi
#include "../../friction/slip_law_srw_psi.hpp"  // SlipLawSRWPsi (TPV104 SRW path)

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int n_pass = 0, n_fail = 0, n_tot = 0;
#define CHECK(cond, msg) do { \
   n_tot++; \
   if (cond) { n_pass++; std::cout << "  PASS: " << msg << "\n"; } \
   else { n_fail++; std::cout << "  FAIL [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace {

constexpr real_t kRho = 2670.0, kCp = 6000.0, kCs = 3464.0;
const real_t kZp = kRho * kCp, kZs = kRho * kCs;
const real_t kEtaP = 0.5 * kZp, kEtaS = 0.5 * kZs;

// Aging (TPV102-style) DOFData; `pert` perturbs the per-QP state (lifted from
// test_friction_substep_iterator_parity.cpp's MakeAgingDOF).
DOFData MakeAgingDOF(real_t pert)
{
   DOFData d;
   d.Zp_plus = kZp; d.Zp_minus = kZp; d.Zs_plus = kZs; d.Zs_minus = kZs;
   d.eta_p = kEtaP; d.eta_s = kEtaS;
   d.sigma_n0 = 50.0e6; d.tau1_0 = 0.0; d.tau2_0 = 29.2e6 + pert * 1.0e6;
   d.a = 0.010; d.b = 0.015; d.Dc = 2.0;
   d.psi = 0.60 + pert * 0.01;
   return d;
}

void MakeQuadrature(int O, real_t dt_macro,
                    std::vector<real_t> &dT, std::vector<real_t> &w)
{
   dT.assign(O, 0.0); w.assign(O, 0.0);
   real_t ds = 0.0, ws = 0.0;
   for (int o = 0; o < O; ++o)
   { dT[o] = 1.0 + 0.5 * o; w[o] = 1.0 + 0.25 * o; ds += dT[o]; ws += w[o]; }
   for (int o = 0; o < O; ++o) { dT[o] = dT[o] / ds * dt_macro; w[o] = w[o] / ws; }
}

void MakeQField(int n, int O, std::vector<std::vector<real_t>> &Qp,
                std::vector<std::vector<real_t>> &Qm)
{
   Qp.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
   Qm.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
   for (int o = 0; o < O; ++o)
      for (int i = 0; i < n; ++i)
      {
         const real_t s = 1.0e-3 * (1.0 + 0.1 * i) * (1.0 + 0.3 * o);
         Qp[o][static_cast<size_t>(i) * NUM_STATE + VZ] = +s;
         Qm[o][static_cast<size_t>(i) * NUM_STATE + VZ] = -s;
      }
}

// Run ONE aging macro-step Advance over `perts.size()` QPs; return the
// resulting psi per QP.  `resample`/`R`/`nbf` configure the Δψ resample.
std::vector<real_t> RunAging(const std::vector<real_t> &perts,
                             const DenseMatrix *R, int nbf, bool resample)
{
   const int n = static_cast<int>(perts.size());
   FaultFaceFlux flux(kRho, kCp, kCs);
   std::vector<DOFData> dof(n);
   for (int i = 0; i < n; ++i) { dof[i] = MakeAgingDOF(perts[i]); }
   std::vector<Vector> coords(n, Vector(3));
   for (int i = 0; i < n; ++i) { coords[i] = 0.0; coords[i](2) = -7500.0; }

   const int O = 3; const real_t dt_macro = 1.0e-4, t0 = 3.7;
   std::vector<real_t> dT, w; MakeQuadrature(O, dt_macro, dT, w);
   std::vector<std::vector<real_t>> Qp, Qm; MakeQField(n, O, Qp, Qm);
   std::vector<real_t> Ip(static_cast<size_t>(NUM_STATE) * n, 0.0);
   std::vector<real_t> Im(static_cast<size_t>(NUM_STATE) * n, 0.0);

   RateStateAgingIterator it(flux, AgingLawPsi(0.015, 1.0e-6, 0.6),
                             FrictionSolver::Method::Brent);
   it.SetSubSteps(dT, w);
   it.SetFaultResample(resample ? R : nullptr, nbf, resample);
   it.Advance(dof, coords, Qp, Qm, dt_macro, t0, Ip.data(), Im.data(),
              [](real_t, real_t) {});

   std::vector<real_t> psi(n);
   for (int i = 0; i < n; ++i) { psi[i] = dof[i].psi; }
   return psi;
}

real_t MaxAbsDiff(const std::vector<real_t> &a, const std::vector<real_t> &b)
{
   real_t m = (a.size() == b.size()) ? 0.0 : 1e300;
   for (size_t i = 0; i < std::min(a.size(), b.size()); ++i)
   { m = std::max(m, std::abs(a[i] - b[i])); }
   return m;
}

std::vector<real_t> InitialPsi(const std::vector<real_t> &perts)
{
   std::vector<real_t> p(perts.size());
   for (size_t i = 0; i < perts.size(); ++i) { p[i] = MakeAgingDOF(perts[i]).psi; }
   return p;
}

// Distinct, non-linear per-QP perturbations ⇒ Δψ carries > N content.
std::vector<real_t> PertsOf(int n, int seed)
{
   std::vector<real_t> p(n);
   for (int i = 0; i < n; ++i)
   { p[i] = std::sin(1.7 * i + 0.5 * seed) + 0.3 * std::cos(0.6 * i + seed); }
   return p;
}

// ---- LSW (TPV205-style) fixture for the Phase-3 R-001 slip resample ----
// Super-critical, PURE strike-slip (tau1_0 = 0 ⇒ V1 = slip1 = 0): with a zero
// Q field the only driver is tau2_0 > mu_s·sigma_n, so the fault slips in the
// +strike direction only.  Then over one macro step starting from slip = 0 the
// path-length increment equals the displacement: dslip_mag[q] = slip2_off[q],
// so the resample collapses to the EXACT relation  slip2_on = R · slip2_off
// (per face), which T5 checks.  Lifted from
// test_friction_substep_iterator_parity.cpp:MakeLswDOF.
DOFData MakeLswDOF(real_t pert)
{
   DOFData d;
   d.Zp_plus = kZp; d.Zp_minus = kZp; d.Zs_plus = kZs; d.Zs_minus = kZs;
   d.eta_p = kEtaP; d.eta_s = kEtaS;
   d.sigma_n0 = 120.0e6; d.tau1_0 = 0.0;
   // 90 MPa > mu_s*sigma_n = 0.677*120e6 = 81.2 MPa ⇒ ruptures; per-QP `pert`
   // varies tau2_0 so slip2_off carries > N content for R to remove.
   d.tau2_0   = 90.0e6 + pert * 1.0e6;
   d.lsw_mu_s = 0.677; d.lsw_mu_d = 0.525; d.lsw_d_c = 0.40;
   d.slip1 = 0.0; d.slip2 = 0.0;
   d.slip_rate_substep_max = 0.0;
   d.sigma_n_substep_min   = std::numeric_limits<real_t>::max();
   return d;
}

// Run ONE LSW macro-step Advance over `perts.size()` QPs with a ZERO Q field
// (pure strike-slip); return slip2 (the strike slip) per QP.  Asserts slip1
// stayed 0 (pure mode-II) so the caller's slip2-only relation is valid.
std::vector<real_t> RunLSW(const std::vector<real_t> &perts,
                           const DenseMatrix *R, int nbf, bool resample,
                           real_t *max_slip1_out = nullptr)
{
   const int n = static_cast<int>(perts.size());
   FaultFaceFlux flux(kRho, kCp, kCs);
   std::vector<DOFData> dof(n);
   for (int i = 0; i < n; ++i) { dof[i] = MakeLswDOF(perts[i]); }
   std::vector<Vector> coords(n, Vector(3));
   for (int i = 0; i < n; ++i) { coords[i] = 0.0; coords[i](2) = -7500.0; }

   const int O = 3; const real_t dt_macro = 1.0e-4, t0 = 3.7;
   std::vector<real_t> dT, w; MakeQuadrature(O, dt_macro, dT, w);
   // Zero Q field ⇒ no normal/dip drive; tau2_0 alone ruptures (strike only).
   std::vector<std::vector<real_t>> Qp, Qm;
   Qp.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
   Qm.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
   std::vector<real_t> Ip(static_cast<size_t>(NUM_STATE) * n, 0.0);
   std::vector<real_t> Im(static_cast<size_t>(NUM_STATE) * n, 0.0);

   LinearSlipWeakeningIterator it(flux);
   it.SetSubSteps(dT, w);
   it.SetFaultResample(resample ? R : nullptr, nbf, resample);
   it.Advance(dof, coords, Qp, Qm, dt_macro, t0, Ip.data(), Im.data(),
              [](real_t, real_t) {});

   real_t m1 = 0.0;
   std::vector<real_t> slip2(n);
   for (int i = 0; i < n; ++i)
   { slip2[i] = dof[i].slip2; m1 = std::max(m1, std::abs(dof[i].slip1)); }
   if (max_slip1_out) { *max_slip1_out = m1; }
   return slip2;
}

// ---- SRW (TPV104-style) fixture for the Δψ resample (R-005) ----
// The slip-law-SRW policy goes through the SAME templated
// RateStateSubStepIterator::Advance as aging, so this confirms the resample
// covers BOTH rate-state policies.  Fixture lifted from
// test_friction_substep_iterator_parity.cpp:MakeSrwDOF.
DOFData MakeSrwDOF(real_t pert)
{
   DOFData d;
   d.Zp_plus = kZp; d.Zp_minus = kZp; d.Zs_plus = kZs; d.Zs_minus = kZs;
   d.eta_p = kEtaP; d.eta_s = kEtaS;
   d.sigma_n0 = 50.0e6; d.tau1_0 = 0.0; d.tau2_0 = 29.38e6 + pert * 1.0e6;
   d.a = 0.008 + pert * 0.0001; d.b = 0.012; d.Dc = 0.40;
   d.psi = 0.564 + pert * 0.01;
   return d;
}

std::vector<real_t> RunSRW(const std::vector<real_t> &perts,
                           const DenseMatrix *R, int nbf, bool resample)
{
   const int n = static_cast<int>(perts.size());
   FaultFaceFlux flux(kRho, kCp, kCs);
   std::vector<DOFData> dof(n);
   for (int i = 0; i < n; ++i) { dof[i] = MakeSrwDOF(perts[i]); }
   std::vector<Vector> coords(n, Vector(3));
   for (int i = 0; i < n; ++i) { coords[i] = 0.0; coords[i](2) = -7500.0; }

   const int O = 3; const real_t dt_macro = 1.0e-4, t0 = 3.7;
   std::vector<real_t> dT, w; MakeQuadrature(O, dt_macro, dT, w);
   std::vector<std::vector<real_t>> Qp, Qm; MakeQField(n, O, Qp, Qm);
   std::vector<real_t> Ip(static_cast<size_t>(NUM_STATE) * n, 0.0);
   std::vector<real_t> Im(static_cast<size_t>(NUM_STATE) * n, 0.0);

   mfem::Vector V_w(n); V_w = 0.1;                         // SRW V_w side-channel
   SlipLawSRWPsi law(0.008, 0.012, 1.0e-6, 0.6, 0.1, 0.1);
   law.SetProductionMode();
   RateStateSlipLawSrwIterator it(flux, law,
                                  FrictionSolver::Method::NewtonRaphsonStable, &V_w);
   it.SetSubSteps(dT, w);
   it.SetFaultResample(resample ? R : nullptr, nbf, resample);
   it.Advance(dof, coords, Qp, Qm, dt_macro, t0, Ip.data(), Im.data(),
              [](real_t, real_t) {});

   std::vector<real_t> psi(n);
   for (int i = 0; i < n; ++i) { psi[i] = dof[i].psi; }
   return psi;
}

} // namespace

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
   std::cout << "\n=== Phase 3: rate-state Δψ resample apply ===\n";

   // ---- T1: R = I (p1 minimal degree-2 rule, #QP=#DOF=3) ⇒ inert ----
   {
      const int order = 1;
      const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2 * order);
      const int nbf = ir.GetNPoints();
      DenseMatrix R; BuildFaultResampleMatrix(Geometry::TRIANGLE, order, ir, R);
      real_t offI = 0.0;
      for (int i = 0; i < nbf; ++i)
         for (int j = 0; j < nbf; ++j)
         { offI = std::max(offI, std::abs(R(i, j) - (i == j ? 1.0 : 0.0))); }
      CHECK(offI < 1e-12, "minimal-rule R is the identity (#QP==#DOF)");

      std::vector<real_t> perts = PertsOf(2 * nbf, 1);     // 2 faces
      auto psi_off = RunAging(perts, &R, nbf, false);
      auto psi_on  = RunAging(perts, &R, nbf, true);
      CHECK(MaxAbsDiff(psi_off, psi_on) < 1e-13,
            "resample with R=I matches no-resample to round-off (inert)");
   }

   // ---- T2: R != I (p1 over-int degree 6, #QP=12) ⇒ psi = psi_before + R·Δψ ----
   {
      const int order = 1;
      const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2 * (order + 2));
      const int nbf = ir.GetNPoints();
      DenseMatrix R; BuildFaultResampleMatrix(Geometry::TRIANGLE, order, ir, R);
      CHECK(nbf > (order + 1) * (order + 2) / 2, "over-integrated rule (#QP > #DOF)");

      std::vector<real_t> perts = PertsOf(nbf, 2);          // 1 face
      auto psi_before = InitialPsi(perts);
      auto psi_off = RunAging(perts, &R, nbf, false);        // psi_before + Δψ
      auto psi_on  = RunAging(perts, &R, nbf, true);         // psi_before + R·Δψ

      std::vector<real_t> dpsi(nbf), proj(nbf), expected(nbf);
      for (int q = 0; q < nbf; ++q) { dpsi[q] = psi_off[q] - psi_before[q]; }
      ApplyFaultResample(R, dpsi.data(), proj.data());
      for (int q = 0; q < nbf; ++q) { expected[q] = psi_before[q] + proj[q]; }

      CHECK(MaxAbsDiff(psi_on, expected) < 1e-13,
            "resample applies psi = psi_before + R·Δψ exactly");
      // Non-vacuous: R removed a > N fraction of the increment (R·Δψ != Δψ).
      // Δψ over one macro step is tiny, so measure RELATIVE to ||Δψ||.
      real_t dpsi_inf = 0.0, rem_inf = 0.0;
      for (int q = 0; q < nbf; ++q)
      {
         dpsi_inf = std::max(dpsi_inf, std::abs(dpsi[q]));
         rem_inf  = std::max(rem_inf, std::abs(proj[q] - dpsi[q]));
      }
      CHECK(rem_inf > 1e-3 * dpsi_inf,
            "resample removes a > N fraction of Δψ (R·Δψ != Δψ; non-vacuous)");
   }

   // ---- T3: per-face independence (2 faces, R != I) ----
   {
      const int order = 1;
      const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2 * (order + 2));
      const int nbf = ir.GetNPoints();
      DenseMatrix R; BuildFaultResampleMatrix(Geometry::TRIANGLE, order, ir, R);

      std::vector<real_t> perts2 = PertsOf(2 * nbf, 3);            // 2 faces
      std::vector<real_t> perts0(perts2.begin(), perts2.begin() + nbf);  // face 0
      auto psi_2face = RunAging(perts2, &R, nbf, true);
      auto psi_1face = RunAging(perts0, &R, nbf, true);
      std::vector<real_t> face0(psi_2face.begin(), psi_2face.begin() + nbf);
      CHECK(MaxAbsDiff(face0, psi_1face) < 1e-13,
            "resample is per-face (face 0 unaffected by face 1's QPs)");
   }

   std::cout << "\n=== Phase 3 (R-001): LSW slip-magnitude resample apply ===\n";

   // ---- T4: LSW, R = I (p1 minimal degree-2 rule, #QP=#DOF=3) ⇒ inert ----
   {
      const int order = 1;
      const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2 * order);
      const int nbf = ir.GetNPoints();
      DenseMatrix R; BuildFaultResampleMatrix(Geometry::TRIANGLE, order, ir, R);

      std::vector<real_t> perts = PertsOf(2 * nbf, 5);     // 2 faces
      real_t m1_off = 0.0, m1_on = 0.0;
      auto s2_off = RunLSW(perts, &R, nbf, false, &m1_off);
      auto s2_on  = RunLSW(perts, &R, nbf, true,  &m1_on);
      CHECK(m1_off < 1e-14 && m1_on < 1e-14, "LSW fixture is pure strike-slip (slip1==0)");
      CHECK(MaxAbsDiff(s2_off, s2_on) < 1e-13,
            "LSW resample with R=I matches no-resample to round-off (inert)");
   }

   // ---- T5: LSW, R != I (p1 over-int degree 6, #QP=12) ----
   //   slip2_on == R · slip2_off (per face), EXACTLY (pure strike-slip from
   //   slip=0 ⇒ path length == displacement ⇒ proj == R·slip2_off), and
   //   non-vacuous (R removed > N content of slip2_off).
   {
      const int order = 1;
      const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2 * (order + 2));
      const int nbf = ir.GetNPoints();
      DenseMatrix R; BuildFaultResampleMatrix(Geometry::TRIANGLE, order, ir, R);
      CHECK(nbf > (order + 1) * (order + 2) / 2, "over-integrated rule (#QP > #DOF)");

      std::vector<real_t> perts = PertsOf(nbf, 6);          // 1 face
      real_t m1_off = 0.0, m1_on = 0.0;
      auto s2_off = RunLSW(perts, &R, nbf, false, &m1_off);
      auto s2_on  = RunLSW(perts, &R, nbf, true,  &m1_on);
      CHECK(m1_off < 1e-14 && m1_on < 1e-14, "LSW fixture is pure strike-slip (slip1==0)");

      std::vector<real_t> expected(nbf);
      ApplyFaultResample(R, s2_off.data(), expected.data());   // R · slip2_off

      real_t off_inf = 0.0, rel_err = 0.0, vac = 0.0;
      for (int q = 0; q < nbf; ++q)
      {
         off_inf = std::max(off_inf, std::abs(s2_off[q]));
         rel_err = std::max(rel_err, std::abs(s2_on[q] - expected[q]));
         vac     = std::max(vac, std::abs(s2_on[q] - s2_off[q]));
      }
      CHECK(rel_err < 1e-9 * off_inf,
            "LSW resample yields slip2 = R·slip2_off per face (rescaled to dealiased mag)");
      CHECK(vac > 1e-3 * off_inf,
            "LSW resample removes a > N fraction of the slip magnitude (non-vacuous)");
   }

   std::cout << "\n=== Phase 3 (R-005): SRW (slip-law) Δψ resample apply ===\n";

   // ---- T6: SRW Δψ resample — same templated Advance as aging ----
   //   psi = psi_before + R·Δψ EXACTLY, on the slip-law-SRW (TPV104) path.
   {
      const int order = 1;
      const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2 * (order + 2));
      const int nbf = ir.GetNPoints();
      DenseMatrix R; BuildFaultResampleMatrix(Geometry::TRIANGLE, order, ir, R);

      std::vector<real_t> perts = PertsOf(nbf, 7);          // 1 face
      std::vector<real_t> psi_before(nbf);
      for (int q = 0; q < nbf; ++q) { psi_before[q] = MakeSrwDOF(perts[q]).psi; }
      auto psi_off = RunSRW(perts, &R, nbf, false);          // psi_before + Δψ
      auto psi_on  = RunSRW(perts, &R, nbf, true);           // psi_before + R·Δψ

      std::vector<real_t> dpsi(nbf), proj(nbf), expected(nbf);
      for (int q = 0; q < nbf; ++q) { dpsi[q] = psi_off[q] - psi_before[q]; }
      ApplyFaultResample(R, dpsi.data(), proj.data());
      for (int q = 0; q < nbf; ++q) { expected[q] = psi_before[q] + proj[q]; }
      CHECK(MaxAbsDiff(psi_on, expected) < 1e-13,
            "SRW resample applies psi = psi_before + R·Δψ exactly (slip-law path)");

      real_t dpsi_inf = 0.0, rem_inf = 0.0;
      for (int q = 0; q < nbf; ++q)
      {
         dpsi_inf = std::max(dpsi_inf, std::abs(dpsi[q]));
         rem_inf  = std::max(rem_inf, std::abs(proj[q] - dpsi[q]));
      }
      CHECK(rem_inf > 1e-3 * dpsi_inf,
            "SRW resample removes a > N fraction of Δψ (non-vacuous)");
   }

   std::cout << "\n========================================\n"
             << "  Results: " << n_pass << " passed, " << n_fail
             << " failed out of " << n_tot << " tests\n"
             << "========================================\n";
   return (n_fail == 0) ? 0 : 1;
}
