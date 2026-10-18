// Nitsche dynamic rupture: the rate-and-state closure with the slip law and
// strong rate weakening (friction/rate_state_srw_closure.hpp; benchmark
// plan rev 1.3, section 4.4, phase B5).
//
//   S-01  the closure is an adapter: SolveSlipRate, StateRate and
//         SeedEquilibrium equal the direct calls (SolveSlipRateVectorPsi
//         with gamma as eta, SlipLawSRWPsi::Rate_SRW, InitialStatePsi with
//         eta = 0) bit for bit at 10^4 random points; CheckParameters.
//   S-02  initial state: psi_ini = 0.563591842632738 at a = 0.01 (TPV104
//         eq. 11, tau = 40 MPa, sigma = 120 MPa, V_ini = 1e-16 m/s) to
//         1e-15, and equal to ComputeInitialPsiTPV104 at a = 0.01, 0.02.
//   S-03  Heun against the exact exponential on a spring-slider driven by
//         the TPV104 overstress (45 MPa, smooth over 1 s), stepped as the
//         operator steps a fault node (predictor psi* with |V^n|, closure
//         with gamma = eta + k dt/2, corrector): Heun and the frozen-V
//         exponential (UpdateStateAnalyticSlipLawSRW with |V^n| in the
//         predictor and |V^{n+1}| in the update) against a fine Heun
//         reference.  Heun is expected at order 2 (asserted: ratios in
//         [3.5, 4.5] for slip and psi), the exponential at order 1
//         (reported).
//   S-04  low-speed robustness: at V_ini = 1e-16 m/s and tau = 40 MPa with
//         gamma from 1e5 to 1e10 Pa s/m the closure returns a finite V
//         equal to V_ini to 1e-9 relative.
//
// Pointwise tests, no mesh.  Usage: ./seas_test_nitsche_srw
//
// Scope (Chunhui Zhao, 2026-10-16).

#include "mfem.hpp"
#include "../../friction/rate_state_srw_closure.hpp"
#include "../../config/tpv104_params.hpp"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

using P = TPV104Params;

static bool SameBits(real_t a, real_t b) { return std::memcmp(&a, &b, sizeof(real_t)) == 0; }

static DieterichRuinaFriction MakeFriction()
{
   return DieterichRuinaFriction({P::V0, P::f0, P::b, P::L});
}

// ------------------------------------------------------------------ S-01
static void TestAdapter()
{
   std::cout << "\n=== S-01 adapter ===\n";
   DieterichRuinaFriction fr = MakeFriction();
   RateStateSRWClosure cl(fr, P::f_w);
   SlipLawSRWPsi ref(0.01, P::b, P::V0, P::f0, P::f_w, 0.1);
   std::mt19937_64 rng(104);
   std::uniform_real_distribution<real_t> U(0.0, 1.0);
   int bad_v = 0, bad_r = 0, bad_s = 0, n = 0;
   for (int k = 0; k < 10000; k++)
   {
      FaultPointState s;
      s.a = 0.005 + 0.02 * U(rng);
      s.Dc = 0.1 + 0.6 * U(rng);
      s.V_w = 0.05 + 1.0 * U(rng);
      s.psi = 0.3 + 0.6 * U(rng);
      const real_t sigma = 50e6 + 100e6 * U(rng);
      const real_t gamma = std::pow(10.0, 5.0 + 4.0 * U(rng));
      const real_t X[2] = {120e6 * (2 * U(rng) - 1), 120e6 * (2 * U(rng) - 1)};
      real_t V[2], Vr[2];
      const real_t Va = cl.SolveSlipRate(X, gamma, sigma, s, V);
      fr.SolveSlipRateVectorPsi(X, s.psi, sigma, gamma, s.a, Vr);
      const real_t Var = std::sqrt(Vr[0] * Vr[0] + Vr[1] * Vr[1]);
      if (!(SameBits(V[0], Vr[0]) && SameBits(V[1], Vr[1]) && SameBits(Va, Var))) { bad_v++; }
      const real_t Vabs = std::pow(10.0, -16.0 + 17.0 * U(rng));
      if (!SameBits(cl.StateRate(Vabs, s), ref.Rate_SRW(Vabs, s.psi, s.Dc, s.V_w, s.a))) { bad_r++; }
      FaultPointState q = s;
      const real_t tau0 = 20e6 + 60e6 * U(rng), Vi = std::pow(10.0, -16.0 + 10.0 * U(rng));
      cl.SeedEquilibrium(tau0, Vi, sigma, q);
      if (!SameBits(q.psi, fr.InitialStatePsi(tau0, Vi, sigma, 0.0, s.a))) { bad_s++; }
      n++;
   }
   std::cout << "  " << n << " points: SolveSlipRate mismatches " << bad_v << ", StateRate "
             << bad_r << ", SeedEquilibrium " << bad_s << "\n";
   TEST_ASSERT(bad_v == 0, "S-01 SolveSlipRate == SolveSlipRateVectorPsi (gamma as eta), bitwise");
   TEST_ASSERT(bad_r == 0, "S-01 StateRate == SlipLawSRWPsi::Rate_SRW(V, psi, Dc, V_w, a), bitwise");
   TEST_ASSERT(bad_s == 0, "S-01 SeedEquilibrium == InitialStatePsi(eta = 0), bitwise");
   TEST_ASSERT(cl.StateLaw().InProductionMode() && cl.WeakeningFriction() == P::f_w,
               "S-01 the SRW law is in production mode with f_w from the constructor");
   FaultPointState g; g.a = 0.01; g.Dc = 0.4; g.V_w = 0.1;
   FaultPointState b1 = g, b2 = g, b3 = g;
   b1.a = 0.0; b2.Dc = 0.0; b3.V_w = 0.0;
   TEST_ASSERT(cl.CheckParameters(g) && !cl.CheckParameters(b1) && !cl.CheckParameters(b2)
               && !cl.CheckParameters(b3), "S-01 CheckParameters: a > 0, Dc > 0, V_w > 0");
}

// ------------------------------------------------------------------ S-02
static void TestInitialState()
{
   std::cout << "\n=== S-02 initial state ===\n";
   DieterichRuinaFriction fr = MakeFriction();
   RateStateSRWClosure cl(fr, P::f_w);
   FaultPointState s;
   s.a = 0.01; s.Dc = P::L; s.V_w = P::V_w_in;
   cl.SeedEquilibrium(P::tau_ini, P::V_ini, P::sigma_n, s);
   FaultPointState s2 = s;
   s2.a = 0.02;
   cl.SeedEquilibrium(P::tau_ini, P::V_ini, P::sigma_n, s2);
   const real_t e1 = std::abs(s.psi - 0.563591842632738);
   const real_t c1 = std::abs(s.psi - ComputeInitialPsiTPV104(0.01));
   const real_t c2 = std::abs(s2.psi - ComputeInitialPsiTPV104(0.02));
   std::cout << std::setprecision(17) << "  psi_ini(a = 0.01) = " << s.psi << " (|diff| to "
             << "0.563591842632738: " << e1 << "; to ComputeInitialPsiTPV104: " << c1
             << "), psi_ini(a = 0.02) = " << s2.psi << " (" << c2 << ")\n";
   TEST_ASSERT(e1 <= 1e-15, "S-02 psi_ini(a = 0.01) = 0.563591842632738 to 1e-15");
   TEST_ASSERT(c1 <= 1e-15 && c2 <= 1e-15, "S-02 equals ComputeInitialPsiTPV104 at a = 0.01, 0.02");
}

// ------------------------------------------------------------------ S-03
struct SliderOut
{
   std::vector<real_t> d, psi, V;   ///< at the dt0 output times
   real_t Vmax = 0.0, ratio_max = 0.0;
};

/// Spring-slider stepped as the operator steps one fault node (boxes 3-11):
/// X = tau_ini + dtau G(t) - k d (k: unloading stiffness), eta radiation
/// damping; scheme 0 = Heun (the operator), 1 = frozen-V exponential.
static SliderOut Slider(int scheme, real_t dt0, int refine, real_t T)
{
   DieterichRuinaFriction fr = MakeFriction();
   RateStateSRWClosure cl(fr, P::f_w);
   const real_t k = 10e6, eta = P::eta_s, dt = dt0 / refine;
   FaultPointState s;
   s.a = P::a_in; s.Dc = P::L; s.V_w = P::V_w_in;
   cl.SeedEquilibrium(P::tau_ini, P::V_ini, P::sigma_n, s);
   real_t d = 0.0, psi = s.psi, V = P::V_ini, t = 0.0;
   const long n = std::lround(T / dt0) * refine;
   SliderOut out;
   for (long i = 0; i < n; i++)
   {
      const real_t ds = d + 0.5 * dt * V;
      FaultPointState st = s;
      st.psi = psi;
      const real_t r0 = cl.StateRate(V, st);
      const real_t ps = (scheme == 0) ? psi + dt * r0
                        : UpdateStateAnalyticSlipLawSRW(psi, V, s.Dc, dt, s.V_w, s.a, P::b,
                                                        P::V0, P::f0, P::f_w);
      t += dt;
      const real_t X[2] = {P::tau_ini + P::nuc_dtau * NucleationTemporal_TPV104(t) - k * ds, 0.0};
      st.psi = ps;
      real_t Vt[2];
      const real_t Vn = cl.SolveSlipRate(X, eta + 0.5 * k * dt, P::sigma_n, st, Vt);
      if (scheme == 0)
      {
         psi = psi + 0.5 * dt * (r0 + cl.StateRate(Vn, st));
      }
      else
      {
         psi = UpdateStateAnalyticSlipLawSRW(psi, Vn, s.Dc, dt, s.V_w, s.a, P::b, P::V0, P::f0,
                                             P::f_w);
      }
      d = ds + 0.5 * dt * Vn;
      V = Vn;
      out.Vmax = std::max(out.Vmax, V);
      out.ratio_max = std::max(out.ratio_max, V * dt / s.Dc);
      if ((i + 1) % refine == 0) { out.d.push_back(d); out.psi.push_back(psi); out.V.push_back(V); }
   }
   return out;
}

static real_t MaxDiff(const std::vector<real_t> &a, const std::vector<real_t> &b)
{
   real_t e = 0.0;
   for (size_t i = 0; i < std::min(a.size(), b.size()); i++) { e = std::max(e, std::abs(a[i] - b[i])); }
   return e;
}

static void TestHeunVsExponential()
{
   std::cout << "\n=== S-03 Heun against the frozen-V exponential (spring-slider) ===\n";
   const real_t dt0 = 0.01, T = 3.0;
   const SliderOut ref = Slider(0, dt0, 256, T);
   std::cout << std::setprecision(6) << "  reference (Heun, dt0/256): max |V| = " << ref.Vmax
             << " m/s, final slip " << ref.d.back() << " m, final psi " << ref.psi.back()
             << ", max |V| dt / L = " << ref.ratio_max << "\n";
   const char *names[2] = {"Heun (operator)", "frozen-V exponential"};
   real_t rH[2][2] = {{0, 0}, {0, 0}};
   for (int sch = 0; sch < 2; sch++)
   {
      real_t ed[3], ep[3], ev[3];
      int i = 0;
      for (int r : {2, 4, 8})
      {
         const SliderOut o = Slider(sch, dt0, r, T);
         ed[i] = MaxDiff(o.d, ref.d);
         ep[i] = MaxDiff(o.psi, ref.psi);
         ev[i] = MaxDiff(o.V, ref.V);
         i++;
      }
      std::cout << std::setprecision(4) << "  " << names[sch] << " at dt0/2, /4, /8:\n"
                << "    slip  err " << ed[0] << ", " << ed[1] << ", " << ed[2] << "  ratios "
                << ed[0] / ed[1] << ", " << ed[1] / ed[2] << "  (order " << std::log2(ed[0] / ed[1])
                << ", " << std::log2(ed[1] / ed[2]) << ")\n"
                << "    psi   err " << ep[0] << ", " << ep[1] << ", " << ep[2] << "  ratios "
                << ep[0] / ep[1] << ", " << ep[1] / ep[2] << "  (order " << std::log2(ep[0] / ep[1])
                << ", " << std::log2(ep[1] / ep[2]) << ")\n"
                << "    |V|   err " << ev[0] << ", " << ev[1] << ", " << ev[2] << "  ratios "
                << ev[0] / ev[1] << ", " << ev[1] / ev[2] << "\n";
      if (sch == 0)
      {
         rH[0][0] = ed[0] / ed[1]; rH[0][1] = ed[1] / ed[2];
         rH[1][0] = ep[0] / ep[1]; rH[1][1] = ep[1] / ep[2];
      }
   }
   auto in = [](real_t r) { return r >= 3.5 && r <= 4.5; };
   TEST_ASSERT(in(rH[0][0]) && in(rH[0][1]) && in(rH[1][0]) && in(rH[1][1]),
               "S-03 Heun (the operator's update) is second order in slip and psi");
}

// ------------------------------------------------------------------ S-04
static void TestLowSpeed()
{
   std::cout << "\n=== S-04 low-speed robustness ===\n";
   DieterichRuinaFriction fr = MakeFriction();
   RateStateSRWClosure cl(fr, P::f_w);
   real_t worst = 0.0;
   int bad = 0, n = 0;
   for (real_t a : {P::a_in, P::a_out})
   {
      FaultPointState s;
      s.a = a; s.Dc = P::L; s.V_w = P::V_w_in;
      cl.SeedEquilibrium(P::tau_ini, P::V_ini, P::sigma_n, s);
      for (int e = 0; e <= 10; e++)
      {
         const real_t gamma = std::pow(10.0, 5.0 + 0.5 * e);
         const real_t X[2] = {P::tau_ini, 0.0};
         real_t V[2];
         const real_t Va = cl.SolveSlipRate(X, gamma, P::sigma_n, s, V);
         if (!std::isfinite(Va) || !(Va > 0.0)) { bad++; continue; }
         worst = std::max(worst, std::abs(Va / P::V_ini - 1.0));
         n++;
      }
   }
   std::cout << std::setprecision(4) << "  " << n << " solves (a = 0.01, 0.02; gamma 1e5..1e10): "
             << bad << " non-finite or zero; max |V / V_ini - 1| = " << worst << "\n";
   TEST_ASSERT(bad == 0 && worst < 1e-9, "S-04 V = V_ini (1e-16 m/s) to 1e-9 at tau = 40 MPa");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   std::cout << "=== Rate-and-state strong-rate-weakening closure (S-01 to S-04) ===\n";
   TestAdapter();
   TestInitialState();
   TestHeunVsExponential();
   TestLowSpeed();
   std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, "
             << num_failed << " failed ===\n";
   return num_failed == 0 ? 0 : 1;
}
