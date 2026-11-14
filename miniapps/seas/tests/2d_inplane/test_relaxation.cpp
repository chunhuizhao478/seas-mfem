// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_relaxation.cpp — the opt-in sigma_n' relaxation of
// the friction law ([cycle.relaxation] t_relax_yr; Chunhui Zhao, 2026-11-02:
// "implement the relaxation, test and verify it, and apply to UCERF 2D plane
// strain, with 10MPa cap"; cycle/normal_stress_relaxation.hpp, design rev 0).
//
//   K1  the QD kernels (friction/device_rate_state.hpp, on the host):
//       RSFaultRHS and RSFaultSlipRate with the reduction R give, bit for
//       bit, what they give without it at sigma_n = (SN + TN) - R; with the
//       floor F as well, at max((SN + TN) - R, F) (the floor acts after the
//       reduction); a node that the reduction takes below F sets
//       kRSSigmaNFloored; a null reduction changes no bit.
//   K2  [cycle.relaxation] opening = "relaxed" (RSFaultInputs::red_in_true)
//       with the tapered floor: the kernels with the reduction R give, bit
//       for bit, what they give without it at sigma_n = (SN + TN) - R with
//       the same tapered floor (the taper reads the relaxed stress); with
//       the flag off the result differs (the taper reads the elastic one);
//       without the taper the flag changes no bit.
//   R1  NormalStressRelaxation against the exact solution for a linear
//       dsig = a t: r = a (t - t_r (1 - exp(-t / t_r))), stepped with
//       accepted steps h; the error is O(h): below a h / 2 + 1e-12 |r| and
//       halving h halves it (ratio 1.8 - 2.2).
//   R2  constant dsig: one step of any length gives the exact
//       r = dsig (1 - exp(-h / t_r)) to 1e-14 (exact for frozen dsig).
//   R3  the stage contract: Fill(t_n) = r_n bit for bit (a stage at the
//       step start, a restart's re-evaluated first stage); Accept(t) =
//       Fill(t) of the step bit for bit; after Restart(t) every Fill is r;
//       SetState(t_n, r, d) then Fill(t_n) = r; after Accept the reduction
//       is r_{n+1}.
//   R4  FrozenNormalStressReduction: NumPoints, the reduction is r(j) and
//       follows the vector it holds; Advance changes nothing.
//   R5  ReductionAt(t, i) = Fill(t)'s Reduction()(i), bit for bit, with and
//       without a frozen mask (the option I stage solves read it per point).
//
// Refusals (--refuse n, each must abort with its message):
//   1  t_relax = 0
//   2  SetDelta with the wrong size
// The deck keys and the driver are checked through seas_fd_cycle_driver
// (test-relaxation in 2d_inplane.mk).
//
// Usage: mpirun -np 1 ./seas_test_relaxation [--refuse n]

#include "mfem.hpp"
#include "../../friction/device_rate_state.hpp"
#include "../../cycle/normal_stress_relaxation.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using mfem::seas::cycle::NormalStressRelaxation;
using mfem::seas::cycle::FrozenNormalStressReduction;

static int num_tests = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kF = 10.0e6;   // the floor [Pa]
static const real_t kYr = 365.25 * 24.0 * 3600.0;

// ---------------------------------------------------------------- K1

/// Per-node inputs of the BP5 fault kernels (state [slip, slip, psi]).
struct KernelInputs
{
   Vector T, X, TN, TP, SN, A, ETA, DC, RED;
   void Resize(int n)
   {
      T.SetSize(2 * n); X.SetSize(3 * n); TN.SetSize(n); TP.SetSize(2 * n);
      SN.SetSize(n); A.SetSize(n); ETA.SetSize(n); DC.SetSize(n); RED.SetSize(n);
   }
   RSFaultInputs In(real_t floor, bool normal, bool red) const
   {
      RSFaultInputs in;
      in.traction = &T;
      in.state = &X;
      in.normal = normal ? &TN : nullptr;
      in.tau_pre = &TP;
      in.sigma_n_dof = &SN;
      in.sigma_n = 0.0;
      in.a = &A;
      in.eta = &ETA;
      in.sigma_floor = floor;
      in.sigma_red = red ? &RED : nullptr;
      return in;
   }
};

static bool SameBits(const Vector &a, const Vector &b, int stride, int n)
{
   for (int i = 0; i < n * stride; i++)
   {
      const real_t x = a(i), y = b(i);
      if (!(x == y || (std::isnan(x) && std::isnan(y)))) { return false; }
   }
   return true;
}

static void TestKernels()
{
   std::cout << "\nK1: the QD fault kernels with a normal-stress reduction\n";
   // SN + TN and the reduction: no change, small, to just above F, below F,
   // into tension, a negative reduction (sigma_law above the elastic one)
   const std::vector<real_t> sig = {80.0e6, 52.0e6, 30.0e6, 45.0e6, 20.0e6, 67.7e6, 15.0e6, 60.0e6};
   const std::vector<real_t> red = {0.0, 3.75e6, 19.0e6, 40.0e6, 25.0e6, -5.0e6, 4.0e6, 50.1e6};
   const int n = static_cast<int>(sig.size());
   KernelInputs k;
   k.Resize(n);
   for (int i = 0; i < n; i++)
   {
      k.TN(i) = -20.0e6 + 3.3e6 * i;        // elastic normal traction
      k.SN(i) = sig[i] - k.TN(i);            // pre-stress, so SN + TN ~ sig
      k.RED(i) = red[i];
      k.TP(2 * i) = 0.1e6 * i;
      k.TP(2 * i + 1) = -0.55 * sig[i] - 0.2e6 * i;
      k.T(2 * i) = 1.0e3 * i;
      k.T(2 * i + 1) = -3.0e3;
      k.X(3 * i) = 0.1 * i;
      k.X(3 * i + 1) = 0.2;
      k.X(3 * i + 2) = 0.55 + 0.02 * i;      // psi
      k.A(i) = (i % 2) ? 0.04 : 0.004;
      k.ETA(i) = 4.62e6;
      k.DC(i) = 0.47;
   }
   // the same nodes without a reduction or normal traction: SN = (SN + TN) - R,
   // in the kernel's order of operations; m2 with the floor applied by hand
   KernelInputs m = k, m2 = k;
   int below = 0;
   for (int i = 0; i < n; i++)
   {
      real_t s = k.SN(i) + k.TN(i);
      s -= k.RED(i);
      m.SN(i) = s;
      m2.SN(i) = std::max(s, kF);
      if (s < kF) { below++; }
   }

   RSDeviceParams prm;
   prm.V0 = 1.0e-6;
   prm.aging_b = 0.03;
   prm.aging_V0 = 1.0e-6;
   prm.aging_f0 = 0.6;
   Array<RSReduceValue> ws;
   Vector Ra(3 * n), Va(2 * n), Rb(3 * n), Vb(2 * n);
   const RSReduceValue a = RSFaultRHS(n, false, prm, k.In(0.0, true, true), k.DC, Ra, Va, ws);
   const RSReduceValue b = RSFaultRHS(n, false, prm, m.In(0.0, false, false), m.DC, Rb, Vb, ws);
   TEST_ASSERT(SameBits(Ra, Rb, 3, n) && SameBits(Va, Vb, 2, n) && a.first == b.first &&
               a.second == b.second,
               "K1 RSFaultRHS: reduction R = the law at (SN + TN) - R, bit for bit (rate, "
               "slip rate, max V, status)");

   Vector Rc(3 * n), Vc(2 * n), Rd(3 * n), Vd(2 * n);
   const RSReduceValue c = RSFaultRHS(n, false, prm, k.In(kF, true, true), k.DC, Rc, Vc, ws);
   const RSReduceValue d = RSFaultRHS(n, false, prm, m2.In(0.0, false, false), m2.DC, Rd, Vd, ws);
   std::cout << "  " << below << " of " << n << " nodes below the floor after the reduction; "
             << "status " << c.second << "\n";
   TEST_ASSERT(below > 0 && SameBits(Rc, Rd, 3, n) && SameBits(Vc, Vd, 2, n) &&
               c.first == d.first,
               "K1 RSFaultRHS: reduction then floor = the law at max((SN + TN) - R, F), bit "
               "for bit");
   TEST_ASSERT((c.second & kRSSigmaNFloored) && !(c.second & kRSSigmaNNonPositive),
               "K1 RSFaultRHS: a node the reduction takes below F sets kRSSigmaNFloored");

   Vector Re(3 * n), Ve(2 * n), Rf(3 * n), Vf(2 * n);
   KernelInputs z = k;
   z.RED = 0.0;
   const RSReduceValue e = RSFaultRHS(n, false, prm, k.In(kF, true, false), k.DC, Re, Ve, ws);
   const RSReduceValue f = RSFaultRHS(n, false, prm, z.In(kF, true, true), z.DC, Rf, Vf, ws);
   TEST_ASSERT(SameBits(Re, Rf, 3, n) && SameBits(Ve, Vf, 2, n) && e.first == f.first &&
               e.second == f.second,
               "K1 RSFaultRHS: a zero reduction = no reduction, bit for bit (r = 0 at t = 0)");

   Vector Sa(2 * n), Sb(2 * n), Sc(2 * n), Sd(2 * n);
   const int sa = RSFaultSlipRate(n, false, prm.V0, k.In(0.0, true, true), Sa, ws);
   const int sb = RSFaultSlipRate(n, false, prm.V0, m.In(0.0, false, false), Sb, ws);
   const int sc = RSFaultSlipRate(n, false, prm.V0, k.In(kF, true, true), Sc, ws);
   const int sd = RSFaultSlipRate(n, false, prm.V0, m2.In(0.0, false, false), Sd, ws);
   TEST_ASSERT(SameBits(Sa, Sb, 2, n) && SameBits(Sc, Sd, 2, n) && SameBits(Sa, Va, 2, n) &&
               SameBits(Sc, Vc, 2, n) && sa == sb && (sc & kRSSigmaNFloored) && sd == 0,
               "K1 RSFaultSlipRate: the same, with and without the floor, and = RSFaultRHS's "
               "slip rate");

   // K2: red_in_true with the tapered floor
   {
      // elastic stress and reduction: above the floor, relaxed below the
      // floor, elastic below the floor with a negative reduction (the
      // relaxed one above it), relaxed into tension, elastic in tension
      // with the relaxed one in compression
      const std::vector<real_t> sg = {80.0e6, 30.0e6, 4.0e6, 20.0e6, -3.0e6, 6.0e6, 15.0e6, 2.0e6};
      const std::vector<real_t> rd = {0.0, 25.0e6, -30.0e6, 25.0e6, -8.0e6, -1.0e6, 4.0e6, 1.5e6};
      KernelInputs q = k;
      for (int i = 0; i < n; i++)
      {
         q.SN(i) = sg[i] - q.TN(i);
         q.RED(i) = rd[i];
      }
      KernelInputs qm = q;
      for (int i = 0; i < n; i++)
      {
         real_t sr = q.SN(i) + q.TN(i);
         sr -= q.RED(i);
         qm.SN(i) = sr;
      }
      RSFaultInputs on = q.In(kF, true, true), off = q.In(kF, true, true),
                    ref = qm.In(kF, false, false), nt_on = q.In(kF, true, true),
                    nt_off = q.In(kF, true, true);
      on.floor_taper = off.floor_taper = ref.floor_taper = true;
      on.red_in_true = nt_on.red_in_true = true;
      Vector R1(3 * n), V1(2 * n), R2(3 * n), V2(2 * n), R3(3 * n), V3(2 * n);
      const RSReduceValue x1 = RSFaultRHS(n, false, prm, on, q.DC, R1, V1, ws);
      const RSReduceValue x2 = RSFaultRHS(n, false, prm, ref, qm.DC, R2, V2, ws);
      const RSReduceValue x3 = RSFaultRHS(n, false, prm, off, q.DC, R3, V3, ws);
      TEST_ASSERT(SameBits(R1, R2, 3, n) && SameBits(V1, V2, 2, n) && x1.first == x2.first &&
                  x1.second == x2.second,
                  "K2 RSFaultRHS: red_in_true with the tapered floor = the tapered law at "
                  "(SN + TN) - R, bit for bit (rate, slip rate, max V, status)");
      TEST_ASSERT(!SameBits(V1, V3, 2, n),
                  "K2 RSFaultRHS: with the flag off the taper reads the elastic stress (a "
                  "different slip rate)");
      Vector S1(2 * n), S2(2 * n);
      const int y1 = RSFaultSlipRate(n, false, prm.V0, on, S1, ws);
      const int y2 = RSFaultSlipRate(n, false, prm.V0, ref, S2, ws);
      TEST_ASSERT(SameBits(S1, S2, 2, n) && SameBits(S1, V1, 2, n) && y1 == y2,
                  "K2 RSFaultSlipRate: the same, and = RSFaultRHS's slip rate");
      Vector R4(3 * n), V4(2 * n), R5(3 * n), V5(2 * n);
      const RSReduceValue x4 = RSFaultRHS(n, false, prm, nt_on, q.DC, R4, V4, ws);
      const RSReduceValue x5 = RSFaultRHS(n, false, prm, nt_off, q.DC, R5, V5, ws);
      TEST_ASSERT(SameBits(R4, R5, 3, n) && SameBits(V4, V5, 2, n) && x4.first == x5.first &&
                  x4.second == x5.second,
                  "K2 RSFaultRHS: without the taper the flag changes no bit");
   }
}

// ---------------------------------------------------------------- R1-R4

/// r of a linear dsig = a t (r(0) = 0): a (t - t_r (1 - exp(-t / t_r))).
static real_t ExactLinear(real_t a, real_t t, real_t tr)
{
   return a * (t + tr * std::expm1(-t / tr));
}

/// Steps 0 -> T with accepted steps h, dsig = a t; returns the max error.
static real_t StepLinear(real_t a, real_t tr, real_t T, real_t h, real_t &r_end)
{
   NormalStressRelaxation rx(1, tr);
   Vector d(1);
   d = 0.0;
   rx.SetDelta(d);
   real_t err = 0.0, t = 0.0;
   const int n = static_cast<int>(std::lround(T / h));
   for (int k = 1; k <= n; k++)
   {
      // a few stages inside the step, as RK does (they do not change r)
      for (real_t c : {0.2, 0.3, 0.8, 8.0 / 9.0}) { rx.Fill(t + c * h); }
      t = k * h;
      d(0) = a * t;
      rx.Accept(t, d);
      err = std::max(err, std::abs(rx.R()(0) - ExactLinear(a, t, tr)));
   }
   r_end = rx.R()(0);
   return err;
}

static void TestRelaxation()
{
   std::cout << "\nR1-R4: NormalStressRelaxation\n";
   const real_t tr = 5000.0 * kYr;              // t_r
   const real_t a = -0.75e3 / kYr;              // UCERF's load drift at B [Pa/s]
   const real_t T = 20000.0 * kYr;
   real_t r1 = 0.0, r2 = 0.0;
   const real_t h1 = 10.0 * kYr, h2 = 5.0 * kYr;
   const real_t e1 = StepLinear(a, tr, T, h1, r1);
   const real_t e2 = StepLinear(a, tr, T, h2, r2);
   const real_t ex = ExactLinear(a, T, tr);
   std::cout << std::scientific << std::setprecision(4)
             << "  linear dsig, t_r 5 kyr, 20 kyr: exact r " << ex / 1e6 << " MPa; max error "
             << e1 / 1e6 << " MPa at h 10 yr, " << e2 / 1e6 << " MPa at h 5 yr (bound a h / 2 = "
             << std::abs(a) * h1 / 2 / 1e6 << " MPa)\n";
   TEST_ASSERT(e1 <= std::abs(a) * h1 / 2 + 1e-12 * std::abs(ex) &&
               e2 <= std::abs(a) * h2 / 2 + 1e-12 * std::abs(ex),
               "R1: error below a h / 2 (dsig frozen at the step start)");
   TEST_ASSERT(e1 / e2 > 1.8 && e1 / e2 < 2.2, "R1: first order in h (ratio " << e1 / e2 << ")");
   // the steady offset: dsig - r -> a t_r
   TEST_ASSERT(std::abs((a * T - r2) - a * tr * (-std::expm1(-T / tr))) < 1e-3 * std::abs(a * tr),
               "R1: dsig - r -> a t_r (1 - exp(-t / t_r)), the steady load offset");

   // R2: constant dsig, one step of any length
   {
      bool ok = true;
      for (real_t h : {1.0, 1.0e7, 0.1 * tr, tr, 30.0 * tr})
      {
         NormalStressRelaxation rx(2, tr);
         Vector d(2);
         d(0) = -40.0e6;
         d(1) = 12.5e6;
         rx.SetDelta(d);
         rx.Accept(h, d);
         for (int i = 0; i < 2; i++)
         {
            const real_t want = -d(i) * std::expm1(-h / tr);
            if (!(std::abs(rx.R()(i) - want) <= 1e-14 * std::abs(d(i)))) { ok = false; }
         }
      }
      TEST_ASSERT(ok, "R2: constant dsig: r = dsig (1 - exp(-h / t_r)) to 1e-14 for h from 1 s "
                  "to 30 t_r");
   }

   // R3: the stage contract
   {
      NormalStressRelaxation rx(3, tr);
      Vector d(3), r(3);
      d(0) = -40.0e6; d(1) = 3.0e6; d(2) = 0.25e6;
      rx.SetDelta(d);
      real_t t = 0.0;
      for (int k = 1; k <= 7; k++)
      {
         const real_t tn = t;
         t = tn + 3.7e8 * k;
         rx.Fill(tn + 0.31 * (t - tn));
         rx.Fill(t);                       // the FSAL stage
         Vector fsal = rx.Reduction();
         d(0) -= 1.0e5; d(1) += 2.0e4 * k; d(2) *= -1.1;
         rx.Accept(t, d);
         TEST_ASSERT(SameBits(fsal, rx.R(), 1, 3) && SameBits(fsal, rx.Reduction(), 1, 3),
                     "R3 step " << k << ": Accept(t) = the FSAL stage's Fill(t), bit for bit");
         rx.Fill(t);
         TEST_ASSERT(SameBits(rx.Reduction(), rx.R(), 1, 3),
                     "R3 step " << k << ": Fill(t_n) = r_n, bit for bit");
      }
      r = rx.R();
      const real_t tn = rx.Time();
      NormalStressRelaxation ry(3, tr);
      ry.SetState(tn, r, rx.Delta());
      ry.Fill(tn);
      TEST_ASSERT(SameBits(ry.Reduction(), r, 1, 3) && ry.Time() == tn,
                  "R3: SetState(t_n, r, dsig) then Fill(t_n) = r, bit for bit (restart)");
      ry.Fill(tn + 1.0e9);
      rx.Fill(tn + 1.0e9);
      TEST_ASSERT(SameBits(ry.Reduction(), rx.Reduction(), 1, 3),
                  "R3: the restarted object's later stages = the original's, bit for bit");
      rx.Restart(tn + 5.0e5);
      bool same = true;
      for (real_t tt : {tn + 5.0e5, tn + 6.0e5, tn + 1.0e11})
      {
         rx.Fill(tt);
         same = same && SameBits(rx.Reduction(), r, 1, 3);
      }
      TEST_ASSERT(same && rx.Time() == tn + 5.0e5,
                  "R3: after Restart(t) (the switch back) every Fill is r, until SetDelta");
   }

   // R4: the frozen modifier
   {
      Vector r(4);
      r(0) = 1.0e6; r(1) = -2.0e6; r(2) = 0.0; r(3) = 40.0e6;
      FrozenNormalStressReduction mod(r);
      const real_t tau[4] = {1, 2, 3, 4}, V[4] = {1, 1, 1, 1};
      mod.Advance(0.01, tau, V);
      bool ok = (mod.NumPoints() == 4);
      for (int j = 0; j < 4; j++) { ok = ok && mod.NormalStressReduction(j) == r(j); }
      r(1) = 7.0e6;
      ok = ok && mod.NormalStressReduction(1) == 7.0e6;
      TEST_ASSERT(ok, "R4: FrozenNormalStressReduction gives r(j), follows its vector, "
                  "Advance changes nothing");
   }

   // R5: ReductionAt = Fill, per point
   {
      NormalStressRelaxation rx(4, tr);
      Vector d(4), mask(4);
      d(0) = -40.0e6; d(1) = 3.0e6; d(2) = 0.25e6; d(3) = -7.0e6;
      mask = 0.0;
      mask(2) = 1.0;
      rx.SetDelta(d);
      rx.Accept(2.0e9, d);
      d(0) -= 3.0e6; d(3) += 1.0e6;
      rx.Accept(5.0e9, d);
      bool ok = true;
      for (int pass = 0; pass < 2; pass++)
      {
         rx.SetFrozen(pass == 0 ? nullptr : &mask);
         for (real_t tt : {5.0e9, 5.0e9 + 1.0, 6.3e9, 5.0e9 + 30.0 * tr})
         {
            rx.Fill(tt);
            for (int i = 0; i < 4; i++)
            {
               ok = ok && rx.ReductionAt(tt, i) == rx.Reduction()(i);
            }
         }
      }
      TEST_ASSERT(ok, "R5: ReductionAt(t, i) = Fill(t)'s Reduction()(i), bit for bit, with "
                  "and without a frozen mask");
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   int refuse = 0;
   OptionsParser args(argc, argv);
   args.AddOption(&refuse, "-r", "--refuse", "Run refusal n (it must abort).");
   args.Parse();
   if (!args.Good()) { args.PrintUsage(std::cout); return 1; }
   if (refuse == 1) { NormalStressRelaxation rx(1, 0.0); return 0; }
   if (refuse == 2)
   {
      NormalStressRelaxation rx(2, 1.0);
      Vector d(3);
      d = 0.0;
      rx.SetDelta(d);
      return 0;
   }
   std::cout << "=== test_relaxation ===\n";
   TestKernels();
   TestRelaxation();
   std::cout << "\n=== " << num_tests - num_failed << " / " << num_tests << " passed ===\n";
   return num_failed == 0 ? 0 : 1;
}
