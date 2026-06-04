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

   std::cout << "\n========================================\n"
             << "  Results: " << n_pass << " passed, " << n_fail
             << " failed out of " << n_tot << " tests\n"
             << "========================================\n";
   return (n_fail == 0) ? 0 : 1;
}
