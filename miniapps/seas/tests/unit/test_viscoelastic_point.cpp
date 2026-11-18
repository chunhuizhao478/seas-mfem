// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_viscoelastic_point.cpp — attenuation plan rev6, phase V1: the fit
// (FitGMB) and the pointwise model (ViscoelasticModel), no mesh, no MPI.
//
//   U1  the fit equals an independent numpy least-squares fit (1e-12) and
//       the plan's numbers (2.20b) at UCERF's Q_S = 0.05 V_S = 173.2
//   U2  a strain ramp eps = r t is integrated exactly by (5.3), any dt
//   U3  step relaxation: a held strain relaxes as sum_l C_l eps0 exp(-omega_l t) (2.2a)
//   U4  harmonic loading: Q^{-1} of the stepped model equals (2.16a), and is
//       within 6 % of 1/Q over the band
//   U5  energy: the work in equals the stored energy (2.24) plus the
//       dissipation (2.26), D >= 0
//   U6  refusals of FitGMB
//   U7  moduli = "reference" equals SeisSol's fitAttenuation formulas
//       (Attenuation.h:76-109 at 2cf185ad) and returns the relaxed moduli
//   U8  the stack: viscoelastic state layout, SetEntryStrain, StiffnessBound
//
// Run: ./seas_test_viscoelastic_point          (all, exit 0)
//      ./seas_test_viscoelastic_point --refuse 1   (must abort: moduli "reference")

#include "mfem.hpp"
#include "../../rheology/viscoelastic.hpp"
#include "../../rheology/rheology_stack.hpp"
#include "test_macros.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

const real_t kRho = 2670.0, kMu = 32038120320.0, kLam = 32038120320.0;   // UCERF decks
const real_t kVs = 3464.0;   // sqrt(kMu / kRho)

GMBBand Band(int L) { GMBBand b; b.L = L; b.f_c = 0.1; b.R = 100.0; return b; }

real_t RelDiff(real_t a, real_t b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }

RheologyPointContext Ctx()
{
   RheologyPointContext c;
   c.lambda0 = kLam;
   c.mu0 = kMu;
   c.rho0 = kRho;
   return c;
}

// numpy.linalg.lstsq of the plan's fit (script in the review notes),
// Q_S = 0.05 * 3464 = 173.2, Q_P = 259.8, f_c 0.1 Hz, R 100
void TestU1()
{
   std::cout << "U1: the fit\n";
   struct Ref { int L; std::vector<real_t> w, YP, YS; real_t lamU, muU; };
   const Ref refs[] = {
      {3, {0.06283185307179587, 0.6283185307179586, 6.283185307179587},
       {0.006539073126882695, 0.005355929377432948, 0.006630968138194183},
       {0.009728949522687012, 0.007996243146883547, 0.00993474782394782},
       32029579868.39039, 32949501572.79459},
      {1, {0.6283185307179586}, {0.007668711656441717}, {0.01148105625717566},
       32036682516.274414, 32410223924.181183},
      {4, {0.06283185307179587, 0.29163962761324647, 1.3536712389686336, 6.283185307179587},
       {0.005625383520258073, 0.003393021603072143, 0.003384025421498611, 0.005720127995909973},
       {0.008364519286647415, 0.005069554003079949, 0.005049404207832025, 0.008576721408420244},
       32029966737.988106, 32929190772.106144}};
   for (const Ref &r : refs)
   {
      const real_t QS = 0.05 * kVs, QP = 1.5 * QS;
      const GMBFit f = FitGMB(QP, QS, Band(r.L), kLam, kMu);
      TEST_ASSERT(f.error.empty(), "U1 L = " << r.L << " fit accepted (" << f.error << ")");
      real_t e = 0.0;
      for (int l = 0; l < r.L; l++)
      {
         e = std::max(e, RelDiff(f.omega[l], r.w[l]));
         e = std::max(e, RelDiff(f.YP[l], r.YP[l]));
         e = std::max(e, RelDiff(f.YS[l], r.YS[l]));
      }
      e = std::max(e, RelDiff(f.lam_U, r.lamU));
      e = std::max(e, RelDiff(f.mu_U, r.muU));
      TEST_ASSERT(e <= 1e-12, "U1 L = " << r.L << " omega, Y^P, Y^S, lam_U, mu_U = numpy (rel "
                  << e << ")");
      // (2.12b) and (2.13) consistency
      real_t c = 0.0;
      for (int l = 0; l < r.L; l++)
      {
         c = std::max(c, RelDiff(f.lam_l[l] + 2.0 * f.mu_l[l], f.YP[l] * (f.lam_U + 2.0 * f.mu_U)));
         c = std::max(c, RelDiff(f.mu_l[l], f.YS[l] * f.mu_U));
         c = std::max(c, RelDiff(f.Ylam[l] * f.lam_U, f.lam_l[l]));
      }
      TEST_ASSERT(c <= 1e-13, "U1 L = " << r.L << " (2.12b), (2.13) hold (rel " << c << ")");
   }
   // the plan's (2.20b), printed to 4 digits
   const GMBFit f = FitGMB(259.8, 173.2, Band(3), kLam, kMu);
   TEST_NEAR(f.YS[0], 0.009729, 5e-7, "U1 (2.20b) Y^S_1 = 0.009729");
   TEST_NEAR(f.YS[1], 0.007996, 5e-7, "U1 (2.20b) Y^S_2 = 0.007996");
   TEST_NEAR(f.YS[2], 0.009935, 5e-7, "U1 (2.20b) Y^S_3 = 0.009935");
   // Q from the relation (2.17a) through the model
   ViscoelasticModel m(Band(3), GMBQuality{});
   const GMBFit g = m.Fit(kLam, kMu, kRho);
   TEST_ASSERT(g.error.empty() && RelDiff(g.YS[0], f.YS[0]) < 1e-12,
               "U1 the model's Q_S = 0.05 V_S gives the same fit");
}

// Exact integration of a strain ramp: xi_l(t) = r / omega_l (1 - exp(-omega_l t))
void TestU2()
{
   std::cout << "U2: strain ramp, exact for any dt\n";
   ViscoelasticModel m(Band(3), GMBQuality{});
   const RheologyPointContext ctx = Ctx();
   const real_t rate[6] = {1e-6, -2e-6, 0.5e-6, 3e-7, -1e-7, 2e-7};
   for (real_t dt : {0.001, 0.37, 5.0})
   {
      std::vector<real_t> q(m.NumStateVars());
      m.InitState(ctx, q.data());
      const real_t zero[6] = {0, 0, 0, 0, 0, 0};
      m.SetEntryStrain(ctx, zero, q.data());
      const int n = static_cast<int>(std::ceil(20.0 / dt));
      real_t dsig[6];
      for (int k = 1; k <= n; k++)
      {
         real_t eps[6];
         for (int a = 0; a < 6; a++) { eps[a] = rate[a] * k * dt; }
         for (int a = 0; a < 6; a++) { dsig[a] = 0.0; }
         m.Update(ctx, eps, dt, q.data(), dsig);
      }
      const real_t t = n * dt;
      real_t err = 0.0, scale = 0.0;
      for (int l = 0; l < 3; l++)
      {
         const real_t w = m.Omega()[l];
         for (int a = 0; a < 6; a++)
         {
            const real_t ex = rate[a] / w * (-std::expm1(-w * t));
            err = std::max(err, std::fabs(q[6 * l + a] - ex));
            scale = std::max(scale, std::fabs(ex));
         }
      }
      TEST_ASSERT(err <= 1e-12 * scale, "U2 dt = " << dt << " s: xi_l exact (rel " << err / scale
                  << ")");
   }
}

// Step relaxation (2.2a): eps0 at t = 0, then held.
void TestU3()
{
   std::cout << "U3: step relaxation\n";
   ViscoelasticModel m(Band(3), GMBQuality{});
   const RheologyPointContext ctx = Ctx();
   const GMBFit f = m.Fit(kLam, kMu, kRho);
   std::vector<real_t> q(m.NumStateVars());
   m.InitState(ctx, q.data());
   const real_t zero[6] = {0, 0, 0, 0, 0, 0};
   m.SetEntryStrain(ctx, zero, q.data());
   const real_t e0 = 1e-5;
   const real_t eps[6] = {0, 0, 0, e0, 0, 0};   // pure shear
   // the jump: one step of 1e-12 s puts the whole increment into xi (phi -> 1)
   real_t dsig[6] = {0, 0, 0, 0, 0, 0};
   m.Update(ctx, eps, 1e-12, q.data(), dsig);
   real_t sum_mu = 0.0;
   for (int l = 0; l < 3; l++) { sum_mu += f.mu_l[l]; }
   TEST_ASSERT(RelDiff(dsig[3], 2.0 * sum_mu * e0) < 1e-9,
               "U3 instant: sigma_xy = 2 mu_U eps0 (Delta = 2 sum mu_l eps0)");
   real_t worst = 0.0;
   real_t t = 1e-12;
   const real_t dt = 0.05;
   for (int k = 0; k < 2000; k++)
   {
      for (int a = 0; a < 6; a++) { dsig[a] = 0.0; }
      m.Update(ctx, eps, dt, q.data(), dsig);
      t += dt;
      real_t ex = 0.0;
      for (int l = 0; l < 3; l++) { ex += 2.0 * f.mu_l[l] * e0 * std::exp(-f.omega[l] * t); }
      worst = std::max(worst, std::fabs(dsig[3] - ex) / (2.0 * sum_mu * e0));
   }
   TEST_ASSERT(worst <= 1e-9, "U3 held: Delta sigma_xy = sum 2 mu_l eps0 exp(-omega_l t) (rel "
               << worst << ")");
   TEST_ASSERT(std::fabs(dsig[3]) < 1e-3 * 2.0 * sum_mu * e0, "U3 relaxes to mu_R after 100 s");
}

// Harmonic loading: steady state, project sigma on sin / cos over whole periods.
void TestU4()
{
   std::cout << "U4: Q^{-1} under harmonic loading\n";
   ViscoelasticModel m(Band(3), GMBQuality{});
   const RheologyPointContext ctx = Ctx();
   const GMBFit f = m.Fit(kLam, kMu, kRho);
   const real_t QS = 0.05 * kVs;
   real_t worst_model = 0.0, worst_band = 0.0;
   for (real_t freq : {0.01, 0.03, 0.0866, 0.1, 0.3, 1.0})
   {
      const real_t w = 2.0 * M_PI * freq, T = 1.0 / freq;
      const int nper = 4000;
      const real_t dt = T / nper;
      std::vector<real_t> q(m.NumStateVars());
      m.InitState(ctx, q.data());
      const real_t zero[6] = {0, 0, 0, 0, 0, 0};
      m.SetEntryStrain(ctx, zero, q.data());
      // run 40 slowest relaxation times so the transient is gone, then 4 periods
      const int nskip = static_cast<int>(std::ceil(40.0 / f.omega[0] / T)) * nper;
      real_t a_s = 0.0, a_c = 0.0;
      const real_t e0 = 1e-6;
      for (int k = 1; k <= nskip + 4 * nper; k++)
      {
         const real_t t = k * dt;
         const real_t eps[6] = {0, 0, 0, e0 * std::sin(w * t), 0, 0};
         real_t dsig[6] = {0, 0, 0, 0, 0, 0};
         m.Update(ctx, eps, dt, q.data(), dsig);
         if (k > nskip)
         {
            const real_t s = 2.0 * kMu * eps[3] + dsig[3];   // total sigma_xy
            a_s += s * std::sin(w * t) * dt;
            a_c += s * std::cos(w * t) * dt;
         }
      }
      // sigma = 2 e0 (Re mu sin + Im mu cos)
      const real_t ReM = a_s / (2.0 * T) / e0, ImM = a_c / (2.0 * T) / e0;
      const real_t qinv = ImM / ReM;
      const real_t qinv_th = gmb::InverseQ(f.omega, f.YS, w);
      worst_model = std::max(worst_model, std::fabs(qinv - qinv_th) / qinv_th);
      worst_band = std::max(worst_band, std::fabs(qinv_th * QS - 1.0));
      std::cout << "    f = " << freq << " Hz: Q^-1 stepped " << qinv << ", (2.16a) " << qinv_th
                << ", 1/Q_S " << 1.0 / QS << "\n";
   }
   TEST_ASSERT(worst_model <= 2e-3, "U4 stepped Q^{-1} = (2.16a) (rel " << worst_model << ")");
   TEST_ASSERT(worst_band <= 0.06, "U4 Q^{-1} within 6 % of 1/Q_S over 0.01-1 Hz ("
               << worst_band << ")");
}

// Energy: work in = W(eps, xi) + int D dt, (2.24)-(2.26)
void TestU5()
{
   std::cout << "U5: energy\n";
   ViscoelasticModel m(Band(3), GMBQuality{});
   const RheologyPointContext ctx = Ctx();
   std::vector<real_t> q(m.NumStateVars());
   m.InitState(ctx, q.data());
   const real_t zero[6] = {0, 0, 0, 0, 0, 0};
   m.SetEntryStrain(ctx, zero, q.data());
   const real_t dt = 1e-3;
   real_t work = 0.0, diss = 0.0, eps_old[6] = {0, 0, 0, 0, 0, 0};
   real_t sig_old[6] = {0, 0, 0, 0, 0, 0};
   bool d_nonneg = true;
   real_t D_old = 0.0;
   for (int k = 1; k <= 30000; k++)
   {
      const real_t t = k * dt;
      const real_t eps[6] = {1e-6 * std::sin(2.0 * M_PI * 0.2 * t),
                             -4e-7 * std::sin(2.0 * M_PI * 0.07 * t), 2e-7 * t / 30.0,
                             5e-7 * std::sin(2.0 * M_PI * 1.3 * t), 0.0,
                             3e-7 * (1.0 - std::cos(2.0 * M_PI * 0.5 * t))};
      real_t dsig[6] = {0, 0, 0, 0, 0, 0};
      m.Update(ctx, eps, dt, q.data(), dsig);
      real_t sig[6], de[6];
      sym6::ApplyIsotropic(kLam, kMu, eps, sig);
      for (int a = 0; a < 6; a++) { sig[a] += dsig[a]; de[a] = eps[a] - eps_old[a]; }
      real_t sm[6];
      for (int a = 0; a < 6; a++) { sm[a] = 0.5 * (sig[a] + sig_old[a]); }
      work += sym6::DoubleDot(sm, de);   // trapezoid
      real_t out[2];
      m.Diagnostics(ctx, q.data(), out);
      d_nonneg = d_nonneg && out[1] >= 0.0;
      diss += 0.5 * (out[1] + D_old) * dt;
      D_old = out[1];
      for (int a = 0; a < 6; a++) { eps_old[a] = eps[a]; sig_old[a] = sig[a]; }
   }
   // W = 1/2 eps : C_R : eps + 1/2 sum xi_l : C_l : xi_l
   real_t sR[6];
   sym6::ApplyIsotropic(kLam, kMu, eps_old, sR);
   real_t W = 0.5 * sym6::DoubleDot(eps_old, sR);
   for (int l = 0; l < 3; l++)
   {
      real_t s[6];
      sym6::ApplyIsotropic(q[6 * 3 + 6 + l], q[7 * 3 + 6 + l], &q[6 * l], s);
      W += 0.5 * sym6::DoubleDot(&q[6 * l], s);
   }
   const real_t bal = std::fabs(work - W - diss) / diss;
   TEST_ASSERT(d_nonneg, "U5 D >= 0 at every step");
   TEST_ASSERT(diss > 0.0 && bal < 1e-4, "U5 work = W + int D dt (rel to int D: " << bal << ")");
}

void TestU6()
{
   std::cout << "U6: refusals\n";
   auto has = [](const GMBFit &f, const char *w) { return f.error.find(w) != std::string::npos; };
   TEST_ASSERT(has(FitGMB(100, 50, Band(0), kLam, kMu), "mechanisms"), "U6 L = 0");
   GMBBand nb = Band(3);
   nb.R = 1.0;
   TEST_ASSERT(has(FitGMB(100, 50, nb, kLam, kMu), "R = 1"), "U6 R = 1 with 3 mechanisms");
   TEST_ASSERT(has(FitGMB(1.0, 50, Band(3), kLam, kMu), "need finite Q > 1"), "U6 Q_P = 1");
   TEST_ASSERT(has(FitGMB(100, 50, Band(3), kLam, 0.0), "lambda ="), "U6 mu = 0");
   // sum Y >= 1 is a guard only: in this band sum Y^S -> 1 from below as
   // Q_S -> 1 (0.956 at Q = 1.5, 0.9992 at Q = 1.01, numpy), and Q <= 1 is
   // refused first.  So the strongest allowed Q still fits below 1.
   {
      const GMBFit f = FitGMB(1.01, 1.01, Band(3), kLam, kMu);
      real_t s = 0.0;
      for (real_t y : f.YS) { s += y; }
      TEST_ASSERT(f.error.empty() && s < 1.0 && std::fabs(s - 0.9991917687041328) < 1e-9,
                  "U6 Q = 1.01: accepted, sum Y^S = 0.99919 < 1 (numpy)");
   }
   // Q_P much larger than Q_S with lambda small: kappa_l < 0
   TEST_ASSERT(has(FitGMB(2000.0, 5.0, Band(3), 0.0, kMu), "kappa"),
               "U6 Q_P >> Q_S, lambda 0: negative bulk part");
   TEST_ASSERT(FitGMB(259.8, 173.2, Band(3), kLam, kMu).error.empty(), "U6 UCERF values accepted");
}

// SeisSol fitAttenuation (2cf185ad, Attenuation.h:76-109), written in its own form
void TestU7()
{
   std::cout << "U7: moduli = \"reference\" (SeisSol)\n";
   const GMBBand b = Band(3);
   const real_t QP = 259.8, QS = 173.2;
   const GMBFit f = FitGMB(QP, QS, b, kLam, kMu, GMBModuli::Reference);
   TEST_ASSERT(f.error.empty(), "U7 reference fit accepted");
   const real_t w0 = 2.0 * M_PI * b.f_c;
   auto seissol = [&](const std::vector<real_t> &Y, real_t Mref)
   {
      real_t psi1 = 1.0, psi2 = 0.0;
      for (int l = 0; l < b.L; l++)
      {
         const real_t r = w0 / f.omega[l];
         psi1 -= Y[l] / (1.0 + r * r);
         psi2 += Y[l] * r / (1.0 + r * r);
      }
      const real_t R = std::sqrt(psi1 * psi1 + psi2 * psi2);
      return Mref * (R + psi1) / (2.0 * R * R);
   };
   const real_t PU = seissol(f.YP, kLam + 2.0 * kMu), muU = seissol(f.YS, kMu);
   TEST_ASSERT(RelDiff(f.mu_U, muU) < 1e-14 && RelDiff(f.lam_U, PU - 2.0 * muU) < 1e-13,
               "U7 unrelaxed moduli = SeisSol's");
   // (2.23) is the choice Re sqrt(rho / M(omega_c)) = sqrt(rho / M_ref): the
   // phase velocity at f_c is the deck's.  With M = M_U r e^{i phi}:
   // M_ref = 2 M_U r^2 / (r + psi_1).
   real_t p1, p2;
   gmb::Psi(f.omega, f.YS, w0, p1, p2);
   const real_t r = std::sqrt(p1 * p1 + p2 * p2);
   TEST_ASSERT(RelDiff(2.0 * f.mu_U * r * r / (r + p1), kMu) < 1e-14,
               "U7 the S phase velocity at f_c is the deck's");
   real_t sS = 0.0;
   for (real_t y : f.YS) { sS += y; }
   TEST_ASSERT(RelDiff(f.mu_R, f.mu_U * (1.0 - sS)) < 1e-14 && f.mu_R < kMu,
               "U7 returns the relaxed moduli, below the deck's (" << f.mu_R / kMu - 1.0 << ")");
}

void TestU8()
{
   std::cout << "U8: the model in a stack\n";
   std::unique_ptr<RheologyStack> s(new RheologyStack);
   s->Add(RheologyKind::Viscoelastic,
          std::unique_ptr<IRheologyModel>(new ViscoelasticModel(Band(3), GMBQuality{})));
   TEST_ASSERT(s->NumStateVars() == 8 * 3 + 6 && s->IsLinearInStrain(),
               "U8 state 8 L + 6, linear in the strain");
   const RheologyPointContext ctx = Ctx();
   std::vector<real_t> q(s->NumStateVars());
   s->InitState(ctx, q.data());
   const real_t e0[6] = {1e-5, 2e-5, 0, 3e-6, 0, 0};
   s->SetEntryStrain(ctx, e0, q.data());
   real_t dsig[6];
   s->Update(ctx, e0, 0.01, q.data(), dsig);
   real_t m = 0.0;
   for (int a = 0; a < 6; a++) { m = std::max(m, std::fabs(dsig[a])); }
   TEST_ASSERT(m == 0.0, "U8 entry strain held: Delta sigma = 0 exactly (rest round trip)");
   const GMBFit f = FitGMB(259.8, 173.2, Band(3), kLam, kMu);
   real_t sP = 0.0;
   for (real_t y : f.YP) { sP += y; }
   TEST_ASSERT(RelDiff(s->StiffnessBound(ctx, q.data()), 1.0 / (1.0 - sP)) < 1e-13,
               "U8 StiffnessBound = 1 / (1 - sum Y^P) = " << 1.0 / (1.0 - sP));
}

} // namespace

int main(int argc, char *argv[])
{
   if (argc == 3 && std::strcmp(argv[1], "--refuse") == 0)
   {
      if (std::atoi(argv[2]) == 1)
      {
         ViscoelasticModel m(Band(3), GMBQuality{}, GMBModuli::Reference);
         (void)m;
      }
      return 0;
   }
   std::cout << "Viscoelastic attenuation pointwise tests (plan rev6, V1: U1-U8)\n";
   TestU1();
   TestU2();
   TestU3();
   TestU4();
   TestU5();
   TestU6();
   TestU7();
   TestU8();
   std::cout << "\n" << num_passed << " / " << num_tests << " passed, " << num_failed
             << " failed\n";
   return num_failed == 0 ? 0 : 1;
}
