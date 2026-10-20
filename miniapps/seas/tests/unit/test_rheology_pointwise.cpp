// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Rheology plan rev 1, unit tests U1-U7 (no mesh): the Drucker-Prager
// pointwise model of rheology/drucker_prager.hpp and the stack around it.
//
//   U1  invariants on hand-made tensors; Y clamped at 0 in tension
//   U2  elastic step: state unchanged, Delta sigma = -C_0 : eps_p (bitwise)
//   U3  return map against a transcription of the SCEC sample Fortran
//       (TPV26/27 description v13, p. 16-17), 10^4 random states
//   U4  Tv = 0: the returned stress lies on the yield surface, sm unchanged
//   U5  relaxation at fixed strain: exponential law; N steps of dt/N equal
//       one step of dt (Appendix B)
//   U6  tr eps_p = 0; sigma_0' + C_0 : (eps - eps_p) is the returned stress
//   U7  frame independence under a rotation
//
// Run: ./seas_test_rheology_pointwise

#include "mfem.hpp"
#include "../../rheology/rheology_model.hpp"
#include "../../rheology/drucker_prager.hpp"
#include "../../rheology/rheology_stack.hpp"
#include "test_macros.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <random>

using namespace mfem;
using namespace mfem::seas;

namespace
{

// TPV27 material and plasticity (description v13, Parts 2 and 6).
const real_t kRho = 2670.0, kVs = 3464.0, kVp = 6000.0;
const real_t kMu = kRho * kVs * kVs;
const real_t kLambda = kRho * (kVp * kVp - 2.0 * kVs * kVs);
const real_t kCohesion = 1.36e6, kNu = 0.1934, kTv = 0.03;

DruckerPragerParams Tpv27Params(real_t Tv = kTv)
{
   DruckerPragerParams p;
   p.cohesion = kCohesion;
   p.bulk_friction = kNu;
   p.relaxation_time = Tv;
   return p;
}

RheologyPointContext Ctx(const real_t s0[6])
{
   RheologyPointContext c;
   for (int i = 0; i < 6; i++) { c.sigma0_eff[i] = s0[i]; }
   c.lambda0 = kLambda;
   c.mu0 = kMu;
   return c;
}

real_t MaxAbs6(const real_t a[6])
{
   real_t m = 0.0;
   for (int i = 0; i < 6; i++) { m = std::max(m, std::abs(a[i])); }
   return m;
}

real_t MaxDiff6(const real_t a[6], const real_t b[6])
{
   real_t m = 0.0;
   for (int i = 0; i < 6; i++) { m = std::max(m, std::abs(a[i] - b[i])); }
   return m;
}

/// sigma^{n+1} = sigma_0' + C_0 : eps + dsig
void TotalStress(const RheologyPointContext &c, const real_t eps[6], const real_t dsig[6],
                 real_t sig[6])
{
   sym6::ApplyIsotropic(c.lambda0, c.mu0, eps, sig);
   for (int i = 0; i < 6; i++) { sig[i] += c.sigma0_eff[i] + dsig[i]; }
}

/// Random symmetric tensor with entries uniform in [-a, a].
void RandSym(std::mt19937_64 &g, real_t a, real_t t[6])
{
   std::uniform_real_distribution<real_t> u(-a, a);
   for (int i = 0; i < 6; i++) { t[i] = u(g); }
}

/// Random deviatoric tensor (trace removed).
void RandDev(std::mt19937_64 &g, real_t a, real_t t[6])
{
   RandSym(g, a, t);
   const real_t m = sym6::Trace(t) / 3.0;
   for (int i = 0; i < 3; i++) { t[i] -= m; }
}

/// The SCEC sample Fortran (description v13, p. 16-17), transcribed line by
/// line: TOTAL stress s, strain increment e (tensor shear), fluid pressure.
void ScecSampleReturnMap(real_t s[6], const real_t e[6], real_t deltat, real_t fluidpresh,
                         real_t mu, real_t lambda, real_t cohes, real_t blkfric, real_t tvisc)
{
   real_t &sxx = s[0], &syy = s[1], &szz = s[2], &sxy = s[3], &syz = s[4], &sxz = s[5];
   const real_t exx = e[0], eyy = e[1], ezz = e[2], exy = e[3], eyz = e[4], exz = e[5];
   const real_t angfric = std::atan(blkfric);
   // ---[ Calculate trial stress tensor ]---
   const real_t etrace = exx + eyy + ezz;
   sxx = sxx + lambda * etrace + 2.0 * mu * exx;
   syy = syy + lambda * etrace + 2.0 * mu * eyy;
   szz = szz + lambda * etrace + 2.0 * mu * ezz;
   sxy = sxy + 2.0 * mu * exy;
   syz = syz + 2.0 * mu * eyz;
   sxz = sxz + 2.0 * mu * exz;
   // ---[ Calculate trial yield stress ]---
   const real_t sm = (sxx + syy + szz) / 3.0;
   const real_t sdxx = sxx - sm, sdyy = syy - sm, sdzz = szz - sm;
   const real_t sdxy = sxy, sdyz = syz, sdxz = sxz;
   const real_t secinv = 0.5 * (sdxx * sdxx + sdyy * sdyy + sdzz * sdzz)
                         + sdxy * sdxy + sdyz * sdyz + sdxz * sdxz;
   const real_t tau = std::sqrt(secinv);
   real_t taulim = cohes * std::cos(angfric) - (sm + fluidpresh) * std::sin(angfric);
   taulim = std::max(real_t(0.0), taulim);
   // ---[ Yielding in shear ]---
   if (tau > taulim)
   {
      const real_t decay = std::exp(-deltat / tvisc);
      const real_t yldfac = decay + (1.0 - decay) * taulim / tau;
      sxx = sdxx * yldfac + sm;
      syy = sdyy * yldfac + sm;
      szz = sdzz * yldfac + sm;
      sxy = sdxy * yldfac;
      syz = sdyz * yldfac;
      sxz = sdxz * yldfac;
   }
}

/// Rotation matrix (row-major) from a random unit quaternion.
void RandRotation(std::mt19937_64 &g, real_t R[9])
{
   std::normal_distribution<real_t> n(0.0, 1.0);
   real_t q[4];
   real_t nq = 0.0;
   for (int i = 0; i < 4; i++) { q[i] = n(g); nq += q[i] * q[i]; }
   nq = std::sqrt(nq);
   for (int i = 0; i < 4; i++) { q[i] /= nq; }
   const real_t w = q[0], x = q[1], y = q[2], z = q[3];
   const real_t M[9] = {1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
                        2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
                        2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)
                       };
   for (int i = 0; i < 9; i++) { R[i] = M[i]; }
}

/// a' = R a R^T on the 6-vector layout.
void Rotate6(const real_t R[9], const real_t a[6], real_t out[6])
{
   real_t A[9], T[9], B[9];
   sym6::ToRowMajor(a, A);
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         T[3 * i + j] = 0.0;
         for (int k = 0; k < 3; k++) { T[3 * i + j] += R[3 * i + k] * A[3 * k + j]; }
      }
   }
   for (int i = 0; i < 3; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         B[3 * i + j] = 0.0;
         for (int k = 0; k < 3; k++) { B[3 * i + j] += T[3 * i + k] * R[3 * j + k]; }
      }
   }
   sym6::FromRowMajor(B, out);
}

} // namespace

// ---------------------------------------------------------------- U1
void TestU1Invariants()
{
   std::cout << "\n=== U1: invariants and yield stress ===\n";
   DruckerPrager dp(Tpv27Params());
   // sm = -60 MPa, s = (-40, 0, 40, 10, 0, 0) MPa, J2 = 1700 MPa^2
   const real_t sig[6] = {-100e6, -60e6, -20e6, 10e6, 0.0, 0.0};
   real_t sm, s[6], tau;
   DruckerPrager::Invariants(sig, sm, s, tau);
   const real_t s_ref[6] = {-40e6, 0.0, 40e6, 10e6, 0.0, 0.0};
   TEST_ASSERT(std::abs(sm - (-60e6)) <= 1e-15 * 60e6, "U1 mean stress -60 MPa");
   TEST_ASSERT(MaxDiff6(s, s_ref) <= 1e-15 * 40e6, "U1 deviator");
   TEST_ASSERT(std::abs(tau - std::sqrt(1700.0) * 1e6) <= 1e-15 * tau, "U1 tau = sqrt(1700) MPa");

   // Y by the other trigonometric route: cos(atan nu) = 1/sqrt(1+nu^2)
   const real_t cphi = 1.0 / std::sqrt(1.0 + kNu * kNu), sphi = kNu * cphi;
   const real_t Y_ref = kCohesion * cphi + 60e6 * sphi;
   TEST_ASSERT(std::abs(dp.YieldStress(sm) - Y_ref) <= 1e-15 * Y_ref, "U1 Y at sm = -60 MPa");
   // tension: c cos(phi) - sm sin(phi) < 0 for sm = +10 MPa, so Y = 0 exactly
   TEST_ASSERT(dp.YieldStress(10e6) == 0.0, "U1 Y clamped at 0 in tension");
   // apex of the cone: sm = c cot(phi) gives Y = 0 up to round-off
   TEST_ASSERT(dp.YieldStress(kCohesion / kNu) <= 1e-15 * kCohesion, "U1 Y = 0 at the apex");
}

// ---------------------------------------------------------------- U2
void TestU2ElasticStep()
{
   std::cout << "\n=== U2: elastic step ===\n";
   DruckerPrager dp(Tpv27Params());
   // deep, nearly isotropic compression: far below yield
   const real_t s0[6] = {-200e6, -200e6, -200e6, 1e6, 0.0, 0.0};
   RheologyPointContext c = Ctx(s0);
   real_t q[DruckerPrager::kNumState] = {1e-5, -2e-5, 1e-5, 3e-6, -1e-6, 2e-6, 4e-5};
   real_t q0[DruckerPrager::kNumState];
   std::memcpy(q0, q, sizeof(q));
   const real_t eps[6] = {1e-5, 2e-5, -3e-5, 1e-6, 2e-6, -1e-6};
   real_t dsig[6] = {};
   const bool inel = dp.Update(c, eps, 1e-3, q, dsig);
   real_t c_ep[6];
   sym6::ApplyIsotropic(kLambda, kMu, q0, c_ep);
   bool bitwise = true;
   for (int i = 0; i < 6; i++) { bitwise = bitwise && dsig[i] == -c_ep[i]; }
   TEST_ASSERT(!inel, "U2 reported elastic");
   TEST_ASSERT(std::memcmp(q, q0, sizeof(q)) == 0, "U2 state unchanged (bitwise)");
   TEST_ASSERT(bitwise, "U2 dsig = -C0:eps_p (bitwise)");

   // the stack with one model gives the same bits
   RheologyStack st;
   st.Add(RheologyKind::DruckerPrager, std::unique_ptr<IRheologyModel>(new DruckerPrager(
                                                                           Tpv27Params())));
   real_t q2[DruckerPrager::kNumState];
   std::memcpy(q2, q0, sizeof(q2));
   real_t dsig2[6];
   st.Update(c, eps, 1e-3, q2, dsig2);
   TEST_ASSERT(std::memcmp(dsig2, dsig, sizeof(dsig)) == 0, "U2 stack = model (bitwise)");
}

// ---------------------------------------------------------------- U3
void TestU3ScecSample()
{
   std::cout << "\n=== U3: return map vs the SCEC sample Fortran ===\n";
   DruckerPrager dp(Tpv27Params());
   std::mt19937_64 g(20261020);
   std::uniform_real_distribution<real_t> ud(0.0, 20000.0), ut(1e-4, 0.1);
   const int n = 10000;
   int n_el = 0, n_pl = 0;
   real_t worst = 0.0;
   for (int k = 0; k < n; k++)
   {
      const real_t depth = ud(g);
      const real_t pf = 1000.0 * 9.8 * depth;
      // effective initial stress: lithostatic-like mean plus a deviator
      real_t s0[6], dev[6];
      RandDev(g, 25e6, dev);
      const real_t sm0 = -(2670.0 - 1000.0) * 9.8 * depth;
      for (int i = 0; i < 6; i++) { s0[i] = dev[i] + (i < 3 ? sm0 : 0.0); }
      RheologyPointContext c = Ctx(s0);
      real_t eps_n[6], ep_n[6], de[6];
      RandSym(g, 4e-4, eps_n);
      RandDev(g, 2e-4, ep_n);
      RandSym(g, 4e-4, de);
      const real_t dt = ut(g);

      // ours: state eps_p^n, total strain eps^{n+1}
      real_t q[DruckerPrager::kNumState];
      for (int i = 0; i < 6; i++) { q[i] = ep_n[i]; }
      q[DruckerPrager::kEta] = 0.0;
      real_t eps1[6], dsig[6] = {}, sig_ours[6];
      for (int i = 0; i < 6; i++) { eps1[i] = eps_n[i] + de[i]; }
      const bool inel = dp.Update(c, eps1, dt, q, dsig);
      TotalStress(c, eps1, dsig, sig_ours);
      (inel ? n_pl : n_el)++;

      // Fortran: total stress at t_n, strain increment, fluid pressure
      real_t ee[6], s_n[6];
      for (int i = 0; i < 6; i++) { ee[i] = eps_n[i] - ep_n[i]; }
      sym6::ApplyIsotropic(kLambda, kMu, ee, s_n);
      for (int i = 0; i < 6; i++) { s_n[i] += s0[i] - (i < 3 ? pf : 0.0); }
      ScecSampleReturnMap(s_n, de, dt, pf, kMu, kLambda, kCohesion, kNu, kTv);
      for (int i = 0; i < 3; i++) { s_n[i] += pf; }   // back to effective

      worst = std::max(worst, MaxDiff6(sig_ours, s_n) / std::max(MaxAbs6(s_n), real_t(1.0)));
   }
   std::cout << "  elastic " << n_el << ", viscoplastic " << n_pl << ", worst relative "
             << worst << "\n";
   TEST_ASSERT(n_el > n / 10 && n_pl > n / 10, "U3 both branches sampled (> 10% each)");
   TEST_ASSERT(worst <= 1e-13, "U3 agrees with the sample code to 1e-13 relative");
}

// ---------------------------------------------------------------- U4
void TestU4RateIndependent()
{
   std::cout << "\n=== U4: Tv = 0 returns onto the yield surface ===\n";
   DruckerPrager dp(Tpv27Params(0.0));
   std::mt19937_64 g(4);
   real_t worst_tau = 0.0, worst_sm = 0.0;
   int n_pl = 0;
   for (int k = 0; k < 1000; k++)
   {
      real_t s0[6];
      RandDev(g, 40e6, s0);
      for (int i = 0; i < 3; i++) { s0[i] += -80e6; }
      RheologyPointContext c = Ctx(s0);
      real_t eps[6], q[DruckerPrager::kNumState] = {}, dsig[6] = {};
      RandSym(g, 5e-4, eps);
      // trial stress (state zero)
      real_t tr[6], sm_tr, s_tr[6], tau_tr;
      const real_t zero[6] = {};
      TotalStress(c, eps, zero, tr);
      DruckerPrager::Invariants(tr, sm_tr, s_tr, tau_tr);
      if (!dp.Update(c, eps, 1e-3, q, dsig)) { continue; }
      n_pl++;
      real_t sig[6], sm, s[6], tau;
      TotalStress(c, eps, dsig, sig);
      DruckerPrager::Invariants(sig, sm, s, tau);
      const real_t Y = dp.YieldStress(sm);
      worst_tau = std::max(worst_tau, std::abs(tau - Y) / std::max(Y, kCohesion));
      worst_sm = std::max(worst_sm, std::abs(sm - sm_tr) / std::abs(sm_tr));
   }
   std::cout << "  " << n_pl << " yielding samples\n";
   TEST_ASSERT(n_pl > 100, "U4 enough yielding samples");
   TEST_ASSERT(worst_tau <= 1e-13, "U4 sqrt(J2) = Y after the return");
   TEST_ASSERT(worst_sm <= 1e-13, "U4 mean stress unchanged");
}

// ---------------------------------------------------------------- U5
void TestU5Relaxation()
{
   std::cout << "\n=== U5: relaxation at fixed strain ===\n";
   DruckerPrager dp(Tpv27Params());
   const real_t s0[6] = {-80e6, -120e6, -100e6, 30e6, 0.0, 5e6};
   RheologyPointContext c = Ctx(s0);
   const real_t eps[6] = {1e-4, -5e-5, 0.0, 2e-4, 0.0, 0.0};
   real_t tr[6], sm, s[6], tau0;
   const real_t zero[6] = {};
   TotalStress(c, eps, zero, tr);
   DruckerPrager::Invariants(tr, sm, s, tau0);
   const real_t Y = dp.YieldStress(sm);
   TEST_ASSERT(tau0 > 1.5 * Y, "U5 starts well above yield");

   const real_t dt = 0.02;
   const int N = 10;
   real_t qa[DruckerPrager::kNumState] = {}, qb[DruckerPrager::kNumState] = {};
   real_t da[6] = {}, db[6] = {};
   dp.Update(c, eps, dt, qa, da);
   for (int k = 0; k < N; k++)
   {
      for (int i = 0; i < 6; i++) { db[i] = 0.0; }
      dp.Update(c, eps, dt / N, qb, db);
   }
   TEST_ASSERT(MaxDiff6(da, db) <= 1e-13 * MaxAbs6(da), "U5 N steps of dt/N = one step of dt");
   TEST_ASSERT(std::abs(qa[DruckerPrager::kEta] - qb[DruckerPrager::kEta])
               <= 1e-13 * qa[DruckerPrager::kEta], "U5 eta the same either way");

   real_t sig[6], sm1, s1[6], tau1;
   TotalStress(c, eps, da, sig);
   DruckerPrager::Invariants(sig, sm1, s1, tau1);
   const real_t tau_exact = Y + (tau0 - Y) * std::exp(-dt / kTv);
   TEST_ASSERT(std::abs(tau1 - tau_exact) <= 1e-13 * tau_exact,
               "U5 sqrt(J2) = Y + (sqrt(J2_0) - Y) exp(-t/Tv)");
}

// ---------------------------------------------------------------- U6
void TestU6Consistency()
{
   std::cout << "\n=== U6: deviatoric plastic strain, stress identity ===\n";
   DruckerPrager dp(Tpv27Params());
   std::mt19937_64 g(6);
   real_t worst_tr = 0.0, worst_sig = 0.0;
   for (int k = 0; k < 1000; k++)
   {
      real_t s0[6];
      RandDev(g, 40e6, s0);
      for (int i = 0; i < 3; i++) { s0[i] += -60e6; }
      RheologyPointContext c = Ctx(s0);
      real_t q[DruckerPrager::kNumState] = {}, eps[6], dsig[6] = {};
      for (int step = 0; step < 5; step++)
      {
         real_t de[6];
         RandSym(g, 2e-4, de);
         for (int i = 0; i < 6; i++) { eps[i] = (step ? eps[i] : 0.0) + de[i]; }
         for (int i = 0; i < 6; i++) { dsig[i] = 0.0; }
         dp.Update(c, eps, 5e-3, q, dsig);
      }
      worst_tr = std::max(worst_tr, std::abs(sym6::Trace(q)) / std::max(MaxAbs6(q), real_t(1e-30)));
      real_t sig[6], ee[6], sig_ref[6];
      TotalStress(c, eps, dsig, sig);
      for (int i = 0; i < 6; i++) { ee[i] = eps[i] - q[i]; }
      sym6::ApplyIsotropic(kLambda, kMu, ee, sig_ref);
      for (int i = 0; i < 6; i++) { sig_ref[i] += s0[i]; }
      worst_sig = std::max(worst_sig, MaxDiff6(sig, sig_ref) / MaxAbs6(sig_ref));
   }
   TEST_ASSERT(worst_tr <= 1e-13, "U6 tr eps_p = 0");
   TEST_ASSERT(worst_sig <= 1e-13, "U6 sigma = sigma_0' + C0:(eps - eps_p)");
}

// ---------------------------------------------------------------- U7
void TestU7FrameIndependence()
{
   std::cout << "\n=== U7: frame independence ===\n";
   DruckerPrager dp(Tpv27Params());
   std::mt19937_64 g(7);
   real_t worst = 0.0, worst_eta = 0.0;
   int n_pl = 0;
   for (int k = 0; k < 1000; k++)
   {
      real_t R[9], s0[6], eps[6], ep[6];
      RandRotation(g, R);
      RandDev(g, 40e6, s0);
      for (int i = 0; i < 3; i++) { s0[i] += -70e6; }
      RandSym(g, 5e-4, eps);
      RandDev(g, 1e-4, ep);
      real_t s0r[6], epsr[6], epr[6];
      Rotate6(R, s0, s0r);
      Rotate6(R, eps, epsr);
      Rotate6(R, ep, epr);
      RheologyPointContext c = Ctx(s0), cr = Ctx(s0r);
      real_t q[DruckerPrager::kNumState], qr[DruckerPrager::kNumState];
      for (int i = 0; i < 6; i++) { q[i] = ep[i]; qr[i] = epr[i]; }
      q[DruckerPrager::kEta] = qr[DruckerPrager::kEta] = 0.0;
      real_t d[6] = {}, dr[6] = {}, d_rot[6];
      const bool a = dp.Update(c, eps, 2e-3, q, d);
      const bool b = dp.Update(cr, epsr, 2e-3, qr, dr);
      if (a != b) { continue; }   // a sample on the yield surface to round-off
      n_pl += a;
      Rotate6(R, d, d_rot);
      worst = std::max(worst, MaxDiff6(d_rot, dr) / std::max(MaxAbs6(dr), real_t(1.0)));
      worst_eta = std::max(worst_eta, std::abs(q[DruckerPrager::kEta] - qr[DruckerPrager::kEta])
                           / std::max(q[DruckerPrager::kEta], real_t(1e-30)));
   }
   std::cout << "  " << n_pl << " yielding samples\n";
   TEST_ASSERT(n_pl > 100, "U7 enough yielding samples");
   TEST_ASSERT(worst <= 1e-13, "U7 dsig rotates with the frame");
   TEST_ASSERT(worst_eta <= 1e-13, "U7 eta frame independent");
}

int main()
{
   std::cout << "Rheology pointwise tests (plan rev 1, U1-U7)\n";
   TestU1Invariants();
   TestU2ElasticStep();
   TestU3ScecSample();
   TestU4RateIndependent();
   TestU5Relaxation();
   TestU6Consistency();
   TestU7FrameIndependence();
   std::cout << "\n" << num_passed << " / " << num_tests << " passed, " << num_failed
             << " failed\n";
   return num_failed == 0 ? 0 : 1;
}
