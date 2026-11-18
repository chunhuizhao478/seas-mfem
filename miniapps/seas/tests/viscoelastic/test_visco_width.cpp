// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Viscoelastic attenuation plan rev6, gate V4: the 2.5D width memory
// (rheology/width_memory.hpp; section 5.6, Result 6, (2.42b), (5.4)-(5.5))
// with the strain memory of V2/V3, on a 2-D box of side 12 km with the
// 2.5D term beta = b mu_R (ElasticityDomainOperator::SetFiniteWidth).
// Material TPV102; band f_c = 0.5 Hz, R = 100, L = 3; Q_S = 20, Q_P = 30;
// b such that the cutoff f_c' = sqrt(b mu_R / rho) / 2 pi = 0.25 Hz.
//
//   T-W1 the uniform mode (k = 0), every face natural, u_pl = 0: the run
//        follows the scalar GMB oscillator rho u'' = -b mu_R u - sum_l b
//        Y^S_l mu_U eta_l, eta_l' = u' - omega_l eta_l (fine RK4) to 1e-3
//        of the amplitude over 8 cutoff periods (dt a quarter of the stable
//        step, so the scheme's phase error stays below that); the decay rate
//        between the peaks of cycles 3 and 8 is Im w of the mode's root of
//        rho w^2 = b mu(w) to 1 %, and 2 Im w / Re w = Q_S^-1(Re w) of (2.16a)
//        to O(1/Q) (mu at the complex root against mu at Re w).  Without the width memory the same mode does not decay
//        (strain only, (2.40) at k = 0: Q_eff = inf).
//   T-W2 damped plane waves, every face Dirichlet with the exact wave:
//        strain memory only, rho w^2 = mu(w) k^2 + mu_R b (2.39), whose
//        loss is (2.40); and both memories, rho w^2 = M(w) k^2 + mu(w) b
//        (Result 7), S and P.  Each run reproduces its own exact wave to
//        within 3 x the elastic 2.5D discretization error and below 0.1 of
//        the gap to the elastic wave.
//   T-W3 energy with the 2.5D terms: a bump in a closed box (u = 0 on every
//        face), both memories: E_h = 1/2 v.M v + 1/2 u.(K + M_beta) u -
//        dt^2/8 r.M^-1 r (r = (K + M_beta) u) + the strain memory energy +
//        1/2 sum_l eta_l.(b Y^S_l mu_U M_1) eta_l never rises (1e-10 E_0),
//        and its drop equals int (D_strain + D_eta) dt to 5 %.
//   T-W4 off: WidthMemory refuses b = 0 (--refuse 1) and the operator
//        refuses it without the 2.5D term (--refuse 2).
//   T-W5 the plate motion: u = u_pl(t) = (t_off + t) U_pl1 (uniform), with
//        t_off = 1000 s, is a steady state: after 200 steps u = u_pl to
//        1e-12 and eta = 0 exactly.
//
// Usage: mpirun -np 1|2 ./seas_test_visco_width [--refuse 1|2]
//
// Scope (Chunhui Zhao, 2026-11-18).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../rheology/rheology_coupler.hpp"
#include "../../rheology/rheology_config.hpp"
#include "../../rheology/width_memory.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using cplx = std::complex<real_t>;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kBox = 12.0e3;
static const real_t kF0 = 0.5, kFc = 0.25;
static const real_t kQs = 20.0, kQp = 30.0;
static const real_t kRho = TPV102Params::rho, kMu = TPV102Params::mu, kLam = TPV102Params::lambda;
static const real_t kB = std::pow(2.0 * M_PI * kFc, 2) * kRho / kMu;   // b [1/m^2]

static real_t RSum(real_t v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return v; }
static real_t RMax(real_t v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return v; }

static Mesh *Box2(int n, int attr)
{
   Mesh *m = new Mesh(Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, kBox, kBox));
   for (int b = 0; b < m->GetNBE(); b++) { m->SetBdrAttribute(b, attr); }
   m->SetAttributes();
   return m;
}

static MaterialField Material() { return MaterialField::MakeConstant(kLam, kMu, kRho); }

static RheologyConfig VEConfig()
{
   RheologyConfig cfg;
   cfg.viscoelastic = true;
   cfg.gmb_band.L = 3;
   cfg.gmb_band.f_c = kF0;
   cfg.gmb_band.R = 100.0;
   cfg.gmb_q.from_vs = false;
   cfg.gmb_q.Q_s = kQs;
   cfg.gmb_q.Q_p = kQp;
   return cfg;
}

static const GMBFit &TheFit()
{
   static const GMBFit f = FitGMB(kQp, kQs, VEConfig().gmb_band, kLam, kMu, GMBModuli::Relaxed);
   return f;
}

/// mu(w) and P(w) of the fit (relaxed + memory).
static cplx MuW(real_t w, bool memory)
{
   const GMBFit &f = TheFit();
   cplx m = f.mu_R;
   if (memory) { for (int l = 0; l < 3; l++) { m += f.mu_l[l] * cplx(0.0, w) / cplx(f.omega[l], w); } }
   return m;
}
static cplx PW(real_t w, bool memory)
{
   const GMBFit &f = TheFit();
   cplx m = f.lam_R + 2.0 * f.mu_R;
   if (memory)
   {
      for (int l = 0; l < 3; l++)
      {
         m += (f.lam_l[l] + 2.0 * f.mu_l[l]) * cplx(0.0, w) / cplx(f.omega[l], w);
      }
   }
   return m;
}

using VecFn = std::function<void(const Vector &, real_t, Vector &)>;

struct Setup
{
   std::unique_ptr<NitscheWaveOperator> op;
   std::unique_ptr<RheologyCoupler> rc;
   std::unique_ptr<WidthMemory> wm;
   std::unique_ptr<ViscoelasticModel> model;
   std::shared_ptr<real_t> t_off;
};

/// The operator with the 2.5D term (U_pl1 uniform with components upl[2]),
/// the strain memory when strain, the width memory when width.
static Setup Make(ParMesh &pm, int p, bool natural, const VecFn &uD, bool strain, bool width,
                  const real_t upl[2] = nullptr, real_t t_off = 0.0)
{
   Setup S;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   if (natural) { bc.natural_attrs = {1}; }
   else
   {
      bc.dirichlet_attrs = {5};
      bc.default_dirichlet_func = [uD](const Vector &x, real_t t, Vector &u) { uD(x, t, u); };
   }
   S.op = std::make_unique<NitscheWaveOperator>(pm, p, Material(), bc);
   auto &fes = S.op->Elasticity().GetFESpace();
   Vector U_pl1(fes.GetVSize());
   U_pl1 = 0.0;
   if (upl)
   {
      const int nd = fes.GetNDofs();
      for (int j = 0; j < nd; j++) { U_pl1(j) = upl[0]; U_pl1(nd + j) = upl[1]; }
   }
   S.t_off = std::make_shared<real_t>(t_off);
   S.op->SetFiniteWidth(kB, U_pl1, S.t_off);
   if (strain)
   {
      auto zero = [](const real_t *, real_t s[9]) { for (int i = 0; i < 9; i++) { s[i] = 0.0; } };
      S.rc = std::make_unique<RheologyCoupler>(fes, S.op->Sample(), S.op->Boundary().dirichlet_attrs,
                                               BuildRheologyStack(VEConfig()), zero);
      S.op->SetRheology(S.rc.get());
   }
   if (width)
   {
      const RheologyConfig cfg = VEConfig();
      S.model = std::make_unique<ViscoelasticModel>(cfg.gmb_band, cfg.gmb_q, cfg.gmb_moduli);
      S.wm = std::make_unique<WidthMemory>(fes, S.op->Sample(), *S.model, kB, U_pl1, S.t_off);
      S.op->SetWidthMemory(S.wm.get());
   }
   return S;
}

static real_t Step(Setup &S, ParMesh &pm, int p)
{
   real_t dt = ComputeMaxDt(S.op->Sample(), pm, p, 0.5 * NitscheStableConstant(2, p));
   // the strain memory's cut (plan 5.9); the width memory adds sum Y^S b mu_U
   // to the stiffest mode's b mu_R, far below the strain part here
   const GMBFit &f = TheFit();
   real_t sY = 0.0;
   for (int l = 0; l < 3; l++) { sY += f.YP[l]; }
   return dt * std::sqrt(1.0 - sY);
}

// ------------------------------------------------------------------ T-W1
/// The scalar oscillator with its memories, RK4 with h; u(t) at the times ts.
static std::vector<real_t> OscillatorRef(real_t u0, bool width, const std::vector<real_t> &ts)
{
   const GMBFit &f = TheFit();
   const real_t muR = f.mu_R, muU = f.mu_U;
   real_t c[3];
   for (int l = 0; l < 3; l++) { c[l] = width ? kB * f.YS[l] * muU : 0.0; }
   auto rhs = [&](const real_t *y, real_t *d) {
      d[0] = y[1];
      real_t F = -kB * muR * y[0];
      for (int l = 0; l < 3; l++) { F -= c[l] * y[2 + l]; }
      d[1] = F / kRho;
      for (int l = 0; l < 3; l++) { d[2 + l] = y[1] - f.omega[l] * y[2 + l]; }
   };
   real_t y[5] = {u0, 0.0, 0.0, 0.0, 0.0};
   const real_t h = 1e-4;
   std::vector<real_t> out;
   real_t t = 0.0;
   for (real_t tt : ts)
   {
      while (t < tt - 1e-12)
      {
         const real_t hh = std::min(h, tt - t);
         real_t k1[5], k2[5], k3[5], k4[5], z[5];
         rhs(y, k1);
         for (int i = 0; i < 5; i++) { z[i] = y[i] + 0.5 * hh * k1[i]; }
         rhs(z, k2);
         for (int i = 0; i < 5; i++) { z[i] = y[i] + 0.5 * hh * k2[i]; }
         rhs(z, k3);
         for (int i = 0; i < 5; i++) { z[i] = y[i] + hh * k3[i]; }
         rhs(z, k4);
         for (int i = 0; i < 5; i++) { y[i] += hh / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]); }
         t += hh;
      }
      out.push_back(y[0]);
   }
   return out;
}

static void TestTW1()
{
   if (g_rank == 0) { std::cout << "\n=== T-W1 the uniform mode decays at Q_S(omega_c) ===\n"; }
   const int p = 2;
   const real_t u0 = 1.0e-3;
   for (int width = 1; width >= 0; width--)
   {
      std::unique_ptr<Mesh> sm(Box2(3, 1));
      ParMesh pm(MPI_COMM_WORLD, *sm);
      sm.reset();
      Setup S = Make(pm, p, true, nullptr, true, width == 1);
      auto &fes = S.op->Elasticity().GetFESpace();
      const int n = fes.GetVSize(), nd = fes.GetNDofs();
      Vector U(n), V(n);
      U = 0.0;
      V = 0.0;
      for (int j = 0; j < nd; j++) { U(j) = u0; }
      S.rc->Reset(U);
      if (S.wm) { S.wm->Reset(U, 0.0); }
      const real_t Tc = 1.0 / kFc, T = 8.0 * Tc;
      real_t dt = 0.25 * Step(S, pm, p);
      const int nsteps = static_cast<int>(std::ceil(T / dt));
      dt = T / nsteps;
      S.op->SetTime(0.0);
      S.op->SetTimeStep(dt);
      std::vector<real_t> ts, uh;
      real_t t = 0.0;
      CentralDifferenceSolver ode;
      ode.Init(*S.op);
      real_t spread = 0.0;
      for (int s = 0; s < nsteps; s++)
      {
         real_t h = dt;
         ode.Step(U, V, t, h);
         ts.push_back(t);
         uh.push_back(U(0));
         for (int j = 0; j < nd; j++)
         {
            spread = std::max(spread, std::max(std::abs(U(j) - U(0)), std::abs(U(nd + j))));
         }
      }
      spread = RMax(spread);
      const std::vector<real_t> ur = OscillatorRef(u0, width == 1, ts);
      real_t err = 0.0;
      for (std::size_t i = 0; i < ur.size(); i++) { err = std::max(err, std::abs(uh[i] - ur[i])); }
      err = RMax(err) / u0;
      // decay rate between the last peak before 3 Tc and the last before 8 Tc
      // (parabolic peak through the three samples around the discrete one)
      auto peak = [&](real_t t1, real_t &tp) {
         std::size_t ip = 1;
         for (std::size_t i = 1; i + 1 < ts.size() && ts[i] < t1; i++)
         {
            if (uh[i] > uh[i - 1] && uh[i] >= uh[i + 1]) { ip = i; }
         }
         const real_t a = uh[ip - 1], b = uh[ip], c = uh[ip + 1], den = a - 2.0 * b + c;
         const real_t s = (den != 0.0) ? 0.5 * (a - c) / den : 0.0;
         tp = ts[ip] + s * dt;
         return b - 0.25 * (a - c) * s;
      };
      real_t t3, t8;
      const real_t A3 = peak(3.0 * Tc, t3), A8 = peak(8.0 * Tc, t8);
      const real_t rate = std::log(A3 / A8) / (t8 - t3);
      // the mode's root of rho w^2 = b mu(w) (width) or b mu_R (strain only), Newton
      cplx w = 2.0 * M_PI * kFc;
      for (int it = 0; it < 50; it++)
      {
         auto g = [&](cplx z) {
            const GMBFit &f = TheFit();
            cplx m = f.mu_R;
            if (width) { for (int l = 0; l < 3; l++) { m += f.mu_l[l] * cplx(0.0, 1.0) * z / (f.omega[l] + cplx(0.0, 1.0) * z); } }
            return kRho * z * z - kB * m;
         };
         const cplx h = 1e-7 * std::abs(w);
         w -= g(w) / ((g(w + h) - g(w - h)) / (2.0 * h));
      }
      // u ~ Re[exp(i w t)]: decay rate Im w for the root with Re w > 0
      const real_t rate_root = std::imag(w);
      const cplx m = MuW(std::real(w), true);
      const real_t Qinv_fit = std::imag(m) / std::real(m), Qinv_root = 2.0 * std::imag(w) / std::real(w);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(5) << "  " << (width ? "both memories" : "strain only")
                   << ": " << nsteps << " steps of " << dt << " s; max |u_h - u_ref| / u0 " << err
                   << ", spread over the DOFs " << spread / u0 << "; decay rate from the peaks "
                   << rate << " 1/s, -Im w of the root " << rate_root << " 1/s (Re w "
                   << std::real(w) << " rad/s, omega_c " << 2.0 * M_PI * kFc << "); 2 Im w / Re w "
                   << Qinv_root << ", Q_S^-1(Re w) " << Qinv_fit << "\n";
      }
      TEST_ASSERT(err < 1e-3 && spread < 1e-12 * u0, "T-W1 " << (width ? "both" : "strain only")
                  << ": the run is the GMB oscillator (1e-3), uniform to round-off");
      if (width)
      {
         TEST_ASSERT(std::abs(rate / rate_root - 1.0) < 0.01,
                     "T-W1 the decay rate is that of the mode's root (1 %)");
         // mu at the complex root and at the real Re w differ by O(Im w / Re w)
         // = O(1/Q), so the two losses agree to that order only
         TEST_ASSERT(std::abs(Qinv_root / Qinv_fit - 1.0) < 1.0 / kQs,
                     "T-W1 the root's loss 2 Im w / Re w is Q_S^-1 at its frequency (to 1/Q)");
      }
      else
      {
         TEST_ASSERT(std::abs(rate) < 1e-3 * std::abs(Qinv_fit) * std::real(w),
                     "T-W1 strain only: the uniform mode does not decay ((2.40), k = 0)");
      }
   }
}

// ------------------------------------------------------------------ T-W2
struct Wave25
{
   real_t n[2], d[2], w, A;
   cplx kappa, mem[3];
   void U(const Vector &x, real_t t, Vector &u) const
   {
      const cplx e = A * std::exp(cplx(0.0, 1.0) * (w * t - kappa * (n[0] * x(0) + n[1] * x(1))));
      u.SetSize(2);
      for (int i = 0; i < 2; i++) { u(i) = std::real(d[i] * e); }
   }
   void V(const Vector &x, real_t t, Vector &v) const
   {
      const cplx e = cplx(0.0, w) * A * std::exp(cplx(0.0, 1.0) * (w * t - kappa * (n[0] * x(0) + n[1] * x(1))));
      v.SetSize(2);
      for (int i = 0; i < 2; i++) { v(i) = std::real(d[i] * e); }
   }
   cplx Phase(const real_t *x) const { return std::exp(cplx(0.0, -1.0) * kappa * (n[0] * x[0] + n[1] * x[1])); }
};

/// strain / width: which memories the wave has (both false: elastic 2.5D).
static Wave25 MakeWave25(bool pwave, bool strain, bool width)
{
   Wave25 W;
   const real_t nn = std::hypot(1.0, 0.7);
   W.n[0] = 1.0 / nn; W.n[1] = 0.7 / nn;
   if (pwave) { W.d[0] = W.n[0]; W.d[1] = W.n[1]; }
   else { W.d[0] = -W.n[1]; W.d[1] = W.n[0]; }
   W.w = 2.0 * M_PI * kF0;
   const GMBFit &f = TheFit();
   for (int l = 0; l < 3; l++) { W.mem[l] = cplx(0.0, W.w) / cplx(f.omega[l], W.w); }
   const cplx Mk = pwave ? PW(W.w, strain) : MuW(W.w, strain);
   const cplx Mb = MuW(W.w, width);
   W.kappa = std::sqrt((kRho * W.w * W.w - Mb * kB) / Mk);
   if (std::real(W.kappa) < 0.0) { W.kappa = -W.kappa; }   // Re > 0; then Im <= 0 (decay)
   W.A = 1.0e-5 / std::abs(W.kappa);
   return W;
}

static real_t RunWave25(int n, int p, const Wave25 &W, bool strain, bool width)
{
   std::unique_ptr<Mesh> sm(Box2(n, 5));
   ParMesh pm(MPI_COMM_WORLD, *sm);
   sm.reset();
   VecFn uD = [W](const Vector &x, real_t t, Vector &u) { W.U(x, t, u); };
   Setup S = Make(pm, p, false, uD, strain, width);
   auto &fes = S.op->Elasticity().GetFESpace();
   ParGridFunction gu(&fes), gv(&fes);
   VectorFunctionCoefficient cu(2, [W](const Vector &x, real_t t, Vector &u) { W.U(x, t, u); });
   VectorFunctionCoefficient cv(2, [W](const Vector &x, real_t t, Vector &v) { W.V(x, t, v); });
   cu.SetTime(0.0);
   cv.SetTime(0.0);
   gu.ProjectCoefficient(cu);
   gv.ProjectCoefficient(cv);
   Vector U(gu), V(gv);
   const GMBFit &f = TheFit();
   if (S.rc)
   {
      S.rc->Reset(U);
      RheologyField &F = S.rc->Field();
      for (int e = 0; e < F.NumElements(); e++)
      {
         for (int i = 0; i < F.NodesPerElement(); i++)
         {
            const cplx g = cplx(0.0, -1.0) * W.kappa * W.A * W.Phase(F.Context(e, i).x);
            cplx eps[6] = {g * W.d[0] * W.n[0], g * W.d[1] * W.n[1], 0.0,
                           0.5 * g * (W.d[0] * W.n[1] + W.d[1] * W.n[0]), 0.0, 0.0
                          };
            real_t *q = F.State(e, i), *ds = F.Dsig(e, i);
            for (int c = 0; c < 6; c++) { ds[c] = 0.0; }
            for (int l = 0; l < 3; l++)
            {
               real_t xi[6], s6[6];
               for (int c = 0; c < 6; c++) { xi[c] = std::real(W.mem[l] * eps[c]); q[6 * l + c] = xi[c]; }
               sym6::ApplyIsotropic(f.lam_l[l], f.mu_l[l], xi, s6);
               for (int c = 0; c < 6; c++) { ds[c] += s6[c]; }
            }
         }
      }
      S.rc->Refresh();
   }
   if (S.wm)
   {
      S.wm->Reset(U, 0.0);
      // eta_l = Re[mem_l u] at every DOF (u_pl = 0): project the complex wave
      for (int l = 0; l < 3; l++)
      {
         const cplx ml = W.mem[l];
         VectorFunctionCoefficient ce(2, [W, ml](const Vector &x, Vector &u) {
            const cplx e = ml * W.A * W.Phase(x.GetData());
            u.SetSize(2);
            for (int i = 0; i < 2; i++) { u(i) = std::real(W.d[i] * e); }
         });
         ParGridFunction ge(&fes);
         ge.ProjectCoefficient(ce);
         std::memcpy(S.wm->Eta(l), ge.GetData(), sizeof(real_t) * ge.Size());
      }
   }
   const real_t T = 1.0 / kF0;
   real_t dt = Step(S, pm, p);
   const int nsteps = static_cast<int>(std::ceil(T / dt));
   dt = T / nsteps;
   S.op->SetTime(0.0);
   S.op->SetTimeStep(dt);
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*S.op);
   for (int s = 0; s < nsteps; s++) { real_t h = dt; ode.Step(U, V, t, h); }
   gu = U;
   VectorFunctionCoefficient ce(2, [W](const Vector &x, real_t tt, Vector &u) { W.U(x, tt, u); });
   ce.SetTime(t);
   ParGridFunction zero(&fes);
   zero = 0.0;
   return gu.ComputeL2Error(ce) / zero.ComputeL2Error(ce);
}

static real_t Gap25(int n, int p, const Wave25 &a, const Wave25 &b)
{
   std::unique_ptr<Mesh> sm(Box2(n, 5));
   ParMesh pm(MPI_COMM_WORLD, *sm);
   H1_FECollection fec(p + 2, 2);
   ParFiniteElementSpace fes(&pm, &fec, 2);
   ParGridFunction g(&fes), zero(&fes);
   VectorFunctionCoefficient ca(2, [a](const Vector &x, real_t t, Vector &u) { a.U(x, t, u); });
   VectorFunctionCoefficient cb(2, [b](const Vector &x, real_t t, Vector &u) { b.U(x, t, u); });
   ca.SetTime(1.0 / kF0);
   cb.SetTime(1.0 / kF0);
   g.ProjectCoefficient(ca);
   zero = 0.0;
   return g.ComputeL2Error(cb) / zero.ComputeL2Error(cb);
}

static void TestTW2()
{
   if (g_rank == 0) { std::cout << "\n=== T-W2 damped 2.5D plane waves ===\n"; }
   const int p = 3, n = 6;
   struct Case { const char *name; bool pwave, strain, width; };
   const Case cases[3] = {{"S, strain memory only (2.39)", false, true, false},
                          {"S, both memories", false, true, true},
                          {"P, both memories", true, true, true}};
   for (const Case &c : cases)
   {
      const Wave25 Wv = MakeWave25(c.pwave, c.strain, c.width);
      const Wave25 We = MakeWave25(c.pwave, false, false);
      const real_t ev = RunWave25(n, p, Wv, c.strain, c.width);
      const real_t ee = RunWave25(n, p, We, false, false);
      const real_t gap = Gap25(n, p, We, Wv);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(5) << "  " << c.name << ": kappa " << Wv.kappa
                   << " 1/m (elastic " << std::real(We.kappa) << "); run vs its wave " << ev
                   << ", elastic run vs elastic wave " << ee << ", gap " << gap << "\n";
      }
      TEST_ASSERT(ev < 0.1 * gap && ev < 3.0 * ee, "T-W2 " << c.name << " reproduces its damped wave");
   }
}

// ------------------------------------------------------------------ T-W3
static void TestTW3()
{
   if (g_rank == 0) { std::cout << "\n=== T-W3 energy with the 2.5D terms ===\n"; }
   const int p = 2, nsteps = 600;
   std::vector<real_t> E[2], D[2];
   real_t dt = 0.0;
   for (int visco = 0; visco <= 1; visco++)
   {
      std::unique_ptr<Mesh> sm(Box2(4, 5));
      ParMesh pm(MPI_COMM_WORLD, *sm);
      sm.reset();
      VecFn zero = [](const Vector &, real_t, Vector &u) { u.SetSize(2); u = 0.0; };
      Setup S = Make(pm, p, false, zero, visco == 1, visco == 1);
      auto &fes = S.op->Elasticity().GetFESpace();
      ParGridFunction gu(&fes);
      VectorFunctionCoefficient cu(2, [](const Vector &x, Vector &u) {
         const real_t r2 = std::pow(x(0) - 0.45 * kBox, 2) + std::pow(x(1) - 0.55 * kBox, 2);
         const real_t g = 1e-2 * std::exp(-r2 / std::pow(0.15 * kBox, 2));
         u.SetSize(2);
         u(0) = g; u(1) = -0.6 * g;
      });
      gu.ProjectCoefficient(cu);
      Vector U(gu), V(U.Size());
      V = 0.0;
      if (S.rc) { S.rc->Reset(U); }
      if (S.wm) { S.wm->Reset(U, 0.0); }
      dt = Step(S, pm, p);
      S.op->SetTime(0.0);
      S.op->SetTimeStep(dt);
      ParBilinearForm mf(&fes);
      ConstantCoefficient rho(kRho);
      mf.AddDomainIntegrator(new VectorMassIntegrator(rho));
      mf.Assemble();
      mf.Finalize();
      Vector Ku(U.Size()), MiKu(U.Size()), Mv(U.Size());
      const GMBFit &f = TheFit();
      auto energy = [&](real_t &En, real_t &Dn) {
         S.op->Elasticity().ApplyStiffness(U, Ku);
         S.op->Mass().Apply(Ku, MiKu);
         mf.SpMat().Mult(V, Mv);
         real_t e = 0.5 * (V * Mv) + 0.5 * (U * Ku) - 0.125 * dt * dt * (Ku * MiKu);
         real_t d = 0.0;
         if (S.rc)
         {
            const RheologyField &F = S.rc->Field();
            MassIntegrator mi;
            DenseMatrix Me;
            for (int el = 0; el < F.NumElements(); el++)
            {
               mi.AssembleElementMatrix(F.NodeElement(), *pm.GetElementTransformation(el), Me);
               for (int l = 0; l < 3; l++)
               {
                  real_t s = 0.0;
                  for (int i = 0; i < F.NodesPerElement(); i++)
                  {
                     real_t ci[6];
                     sym6::ApplyIsotropic(f.lam_l[l], f.mu_l[l], F.State(el, i) + 6 * l, ci);
                     for (int j = 0; j < F.NodesPerElement(); j++)
                     {
                        s += Me(i, j) * sym6::DoubleDot(ci, F.State(el, j) + 6 * l);
                     }
                  }
                  e += 0.5 * s;
                  d += f.omega[l] * s;
               }
            }
         }
         if (S.wm)
         {
            real_t ew, dw;
            S.wm->Energy(ew, dw);
            e += ew;
            d += dw;
         }
         En = RSum(e);
         Dn = RSum(d);
      };
      CentralDifferenceSolver ode;
      ode.Init(*S.op);
      real_t t = 0.0, En, Dn;
      energy(En, Dn);
      E[visco].push_back(En);
      D[visco].push_back(Dn);
      for (int s = 0; s < nsteps; s++)
      {
         real_t h = dt;
         ode.Step(U, V, t, h);
         energy(En, Dn);
         E[visco].push_back(En);
         D[visco].push_back(Dn);
      }
   }
   real_t el_dev = 0.0, rise = 0.0, intD = 0.0;
   for (int s = 1; s <= nsteps; s++)
   {
      el_dev = std::max(el_dev, std::abs(E[0][s] - E[0][0]));
      rise = std::max(rise, E[1][s] - E[1][s - 1]);
      intD += 0.5 * dt * (D[1][s] + D[1][s - 1]);
   }
   const real_t drop = E[1][0] - E[1][nsteps];
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  p = " << p << ", " << nsteps << " steps of " << dt
                << " s: elastic 2.5D |E - E0| / E0 " << el_dev / E[0][0] << "; with both memories: "
                << "largest rise / E0 " << rise / E[1][0] << ", drop / E0 " << drop / E[1][0]
                << ", drop / int D " << drop / intD << "\n";
   }
   TEST_ASSERT(el_dev < 1e-10 * E[0][0], "T-W3 the elastic 2.5D run keeps E_h (1e-10)");
   TEST_ASSERT(rise < 1e-10 * E[1][0], "T-W3 E_h with both memories never rises");
   TEST_ASSERT(drop > 1e-3 * E[1][0] && std::abs(drop / intD - 1.0) < 0.05,
               "T-W3 the drop is int (D_strain + D_eta) dt to 5 %");
}

// ------------------------------------------------------------------ T-W5
static void TestTW5()
{
   if (g_rank == 0) { std::cout << "\n=== T-W5 the plate motion is a steady state ===\n"; }
   const int p = 2, nsteps = 200;
   std::unique_ptr<Mesh> sm(Box2(3, 1));
   ParMesh pm(MPI_COMM_WORLD, *sm);
   sm.reset();
   const real_t upl[2] = {1.0e-9, -0.4e-9};   // [m] at physical time 1 s
   const real_t toff = 1000.0;
   Setup S = Make(pm, p, true, nullptr, true, true, upl, toff);
   auto &fes = S.op->Elasticity().GetFESpace();
   const int n = fes.GetVSize(), nd = fes.GetNDofs();
   Vector U(n), V(n);
   for (int j = 0; j < nd; j++)
   {
      U(j) = toff * upl[0]; U(nd + j) = toff * upl[1];
      V(j) = upl[0]; V(nd + j) = upl[1];
   }
   S.rc->Reset(U);
   S.wm->Reset(U, 0.0);
   const real_t dt = Step(S, pm, p);
   S.op->SetTime(0.0);
   S.op->SetTimeStep(dt);
   real_t t = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(*S.op);
   for (int s = 0; s < nsteps; s++) { real_t h = dt; ode.Step(U, V, t, h); }
   real_t err = 0.0, eta = 0.0;
   for (int j = 0; j < nd; j++)
   {
      err = std::max(err, std::abs(U(j) - (toff + t) * upl[0]));
      err = std::max(err, std::abs(U(nd + j) - (toff + t) * upl[1]));
   }
   for (int l = 0; l < 3; l++) { for (int k = 0; k < n; k++) { eta = std::max(eta, std::abs(S.wm->Eta(l)[k])); } }
   err = RMax(err) / ((toff + t) * std::abs(upl[0]));
   eta = RMax(eta);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(5) << "  t_off " << toff << " s, " << nsteps << " steps: max |u - u_pl| / |u_pl| "
                << err << ", max |eta| " << eta << " m\n";
   }
   TEST_ASSERT(err < 1e-12, "T-W5 u = u_pl(t) to 1e-12");
   TEST_ASSERT(eta < 1e-12 * toff * std::abs(upl[0]), "T-W5 eta = 0 (to round-off of w = u - u_pl)");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   int refuse = 0;
   for (int a = 1; a + 1 < argc; a++)
   {
      if (std::string(argv[a]) == "--refuse") { refuse = std::atoi(argv[a + 1]); }
   }
   if (refuse > 0)
   {
      std::unique_ptr<Mesh> sm(Box2(2, 1));
      ParMesh pm(MPI_COMM_WORLD, *sm);
      BoundaryConfig bc;
      bc.fault_attr = 3;
      bc.natural_attrs = {1};
      NitscheWaveOperator op(pm, 1, Material(), bc);
      auto &fes = op.Elasticity().GetFESpace();
      Vector Upl(fes.GetVSize());
      Upl = 0.0;
      const RheologyConfig cfg = VEConfig();
      ViscoelasticModel model(cfg.gmb_band, cfg.gmb_q, cfg.gmb_moduli);
      if (refuse == 1) { WidthMemory wm(fes, op.Sample(), model, 0.0, Upl, nullptr); }
      if (refuse == 2)
      {
         WidthMemory wm(fes, op.Sample(), model, kB, Upl, nullptr);
         op.SetWidthMemory(&wm);
      }
      return 0;
   }
   if (g_rank == 0)
   {
      std::cout << "=== 2.5D width memory (attenuation plan rev6, V4) on " << Mpi::WorldSize()
                << " rank(s); b = " << kB << " 1/m^2 (cutoff " << kFc << " Hz) ===\n";
   }
   std::string only = (argc > 2 && std::string(argv[1]) == "--only") ? argv[2] : "";
   auto run = [&](const char *nm) { return only.empty() || only.find(nm) != std::string::npos; };
   if (run("TW5")) { TestTW5(); }
   if (run("TW1")) { TestTW1(); }
   if (run("TW3")) { TestTW3(); }
   if (run("TW2")) { TestTW2(); }
   if (g_rank == 0)
   {
      std::cout << "\n" << num_passed << " passed, " << num_failed << " failed, " << num_tests
                << " total\n";
   }
   return num_failed == 0 ? 0 : 1;
}
