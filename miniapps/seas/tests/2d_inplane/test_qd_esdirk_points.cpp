// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_qd_esdirk_points.cpp — option I, phase I1
// (cycle/qd_implicit.hpp): the fault system of the implicit QD step on
// small dense G's, PETSc arkimex run fully implicit (plan "Option I:
// implementation and test plan", tests T0 and T1).
//
//   T0a  the Newton matrix matches central differences of F (J x), on a
//        6-point chain with open, closed and closed-with-gap points, two
//        shifts, the whole-law taper (sigma_law = max(sigma, 0)) and the
//        hard floor (max(sigma, 10 MPa)).
//   T0b  the preconditioner inverts it (J P b = b); in PETSc's slope mode
//        (shift 1e16) it is the identity on slip and psi and 0 on the gaps.
//   T0d  consistent start (plan 2.10): from g = 0 with three points in
//        tension, g >= 0, sigma >= 0 and g sigma = 0 at every point, and the
//        gap columns of exactly the points that opened.
//   T0c  a 6-point adaptive run with slip events and opening: the same
//        steps and final state bit for bit on all ranks and on one rank.
//   T1a  one point at zero strength, delta' = -k delta / eta: each fixed
//        step multiplies delta by the stability function R(-h k / eta) of
//        the method's implicit tableau (ARK3 and ARK4), h = 1e-2 .. 1e7 s.
//   T1b  stable spring slider (k = 2 k_c), fixed steps halved 5 times
//        against an RK4 reference with tiny steps: order >= 2.8 (ARK3) and
//        >= 3.6 (ARK4).
//   T1c  one point pulled into tension and back: at every accepted state
//        sigma >= 0, g >= 0 and g sigma = 0; psi frozen while open; one gap
//        column fetched; closed with g = 0 at the end.
//   T1d  unstable spring slider (k = 0.5 k_c) through one event, adaptive
//        (atol 1e-7): the event time within 1% of an RK 5dp reference with
//        atol 1e-11.
//
//   T1e  [cycle.relaxation] opening = "relaxed" (Problem::normal_reduction
//        with cycle/normal_stress_relaxation.hpp, accepted after each step
//        as QDPhase does): one point unclamped at a constant rate rho.
//        (a) rho t_r = sigma0 / 2: the elastic stress goes to -2 sigma0, the
//        relaxed one follows sigma0 - rho t_r (1 - exp(-t / t_r)) and the
//        point never opens.  (b) rho t_r = 2 sigma0 until T/4, then no
//        change: it opens at t* = -t_r ln(1 - sigma0 / (rho t_r)), the gap
//        grows at (rho - sigma0 / t_r) / gamma, shrinks at sigma0 / (t_r
//        gamma) after T/4 and closes; g >= 0, s >= 0, g s = 0 at every
//        accepted state (s the relaxed stress); s is back above 0.9 sigma0
//        at the end.
//
// Run on 1, 2 and 4 ranks (make -f tests/2d_inplane/2d_inplane.mk
// test-qd-esdirk-points).

#include "mfem.hpp"
#include "../../cycle/qd_implicit.hpp"
#include "../../cycle/normal_stress_relaxation.hpp"

#include <cmath>
#include <cstdio>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using mfem::seas::cycle::QDImplicitFault;
using mfem::seas::cycle::QDImplicitSettings;

namespace
{

int rank = 0, nranks = 1, fails = 0;
double newton_tol = QDImplicitSettings().newton_tol;   ///< -qd_newton_tol (cost study)

void Check(bool ok, const std::string &what)
{
   int all = ok ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   if (rank == 0) { std::cout << (all ? "  PASSED: " : "  FAILED: ") << what << std::endl; }
   if (!all) { fails++; }
}

void Say(const std::string &s) { if (rank == 0) { std::cout << s << std::endl; } }

std::string Fmt(double x, int prec = 6)
{
   std::ostringstream o;
   o << std::setprecision(prec) << x;
   return o.str();
}

// ---- physics (BP-like rate and state, aging law) ----------------------------

struct Phys
{
   double a = 0.01, b = 0.015, Dc = 0.01, V0 = 1e-6, f0 = 0.6;
   double eta = 4.6e6;      // ~ mu / (2 cs), Pa s/m
   double sigma = 50e6;     // Pa
   double Vpl = 1e-9;       // m/s
   double kc() const { return sigma * (b - a) / Dc; }
};
const Phys ph;

double PsiSS(double V) { return ph.f0 - ph.b * std::log(V / ph.V0); }
double Fric(double V, double psi) { return seas::RSFrictionCoefficientPsi(V, psi, ph.a, ph.V0); }

/// Signed slip rate of the kernel: tau + sigma f(V, psi) + eta V = 0.
double SlipRate(double tau, double psi, double sigma)
{
   int st = 0;
   const double Va = seas::RSSolveSlipRatePsi(std::abs(tau), psi, sigma, ph.eta, ph.a, ph.V0, st);
   return (tau > 0.0) ? -Va : ((tau < 0.0) ? Va : 0.0);
}

double AgingRate(double V, double psi)
{
   const seas::RSDeviceParams prm{ph.V0, ph.b, ph.V0, ph.f0};
   return seas::RSAgingRatePsi(std::abs(V), psi, ph.Dc, prm);
}

// ---- problems ----------------------------------------------------------------

void Partition(MPI_Comm comm, int N, int &lo, int &n)
{
   int r = 0, p = 1;
   MPI_Comm_rank(comm, &r);
   MPI_Comm_size(comm, &p);
   lo = static_cast<int>(static_cast<long>(N) * r / p);
   n = static_cast<int>(static_cast<long>(N) * (r + 1) / p) - lo;
}

void FillOwned(QDImplicitFault::Problem &P, const std::function<double(int)> &tau_pre,
               const std::function<double(int)> &sigma0)
{
   P.tau_pre.resize(P.n);
   P.sigma0.resize(P.n);
   for (int o = 0; o < P.n; o++)
   {
      P.tau_pre[o] = tau_pre(P.lo + o);
      P.sigma0[o] = sigma0(P.lo + o);
   }
   P.a.assign(P.n, ph.a);
   P.eta.assign(P.n, ph.eta);
   P.Dc.assign(P.n, ph.Dc);
   P.b = ph.b;
   P.f0 = ph.f0;
   P.V0 = ph.V0;
}

/// One point: tau = tau_pre + k delta + gt g + t TDt, sigma = sigma0 + gam g.
QDImplicitFault::Problem OnePoint(MPI_Comm comm, double k, double tau_pre, double sigma0,
                                  double TDt, double gt = 0.0, double gam = 0.0)
{
   QDImplicitFault::Problem P;
   P.N = 1;
   P.Gt = {k};
   P.Gs = {0.0};
   P.TDt = {TDt};
   P.TDs = {0.0};
   Partition(comm, 1, P.lo, P.n);
   FillOwned(P, [=](int) { return tau_pre; }, [=](int) { return sigma0; });
   if (gam > 0.0)
   {
      P.gap_column = [=](int, std::vector<real_t> &t, std::vector<real_t> &s)
      {
         t = {gt};
         s = {gam};
      };
   }
   return P;
}

/// A chain of N points: Gt = k (I - 0.3 L), Gs = 0.05 k (right - left
/// neighbour), gap column l: Ggs = gam (e_l - 0.2 neighbours), Ggt = 0.1 k
/// e_l + 0.02 k neighbours (L: the neighbour matrix).
QDImplicitFault::Problem Chain(MPI_Comm comm, int N, double k, double gam,
                               const std::function<double(int)> &tau_pre,
                               const std::function<double(int)> &sigma0, double TDt_scale)
{
   QDImplicitFault::Problem P;
   P.N = N;
   P.Gt.assign(size_t(N) * N, 0.0);
   P.Gs.assign(size_t(N) * N, 0.0);
   for (int j = 0; j < N; j++)
   {
      P.Gt[size_t(j) * N + j] = k;
      if (j > 0)
      {
         P.Gt[size_t(j) * N + j - 1] = -0.3 * k;
         P.Gs[size_t(j) * N + j - 1] = -0.05 * k;
      }
      if (j + 1 < N)
      {
         P.Gt[size_t(j) * N + j + 1] = -0.3 * k;
         P.Gs[size_t(j) * N + j + 1] = 0.05 * k;
      }
   }
   P.TDt.assign(N, 0.0);
   P.TDs.assign(N, 0.0);
   for (int j = 0; j < N; j++)
   {
      double rs = 0.0;
      for (int l = 0; l < N; l++) { rs += P.Gt[size_t(j) * N + l]; }
      P.TDt[j] = -rs * TDt_scale;
   }
   Partition(comm, N, P.lo, P.n);
   FillOwned(P, tau_pre, sigma0);
   P.gap_column = [=](int l, std::vector<real_t> &t, std::vector<real_t> &s)
   {
      t.assign(N, 0.0);
      s.assign(N, 0.0);
      t[l] = 0.1 * k;
      s[l] = gam;
      for (int m : {l - 1, l + 1})
      {
         if (m < 0 || m >= N) { continue; }
         t[m] = 0.02 * k;
         s[m] = -0.2 * gam;
      }
   };
   return P;
}

std::vector<real_t> GatherAll(MPI_Comm comm, const real_t *mine, int nmine)
{
   int p = 1;
   MPI_Comm_size(comm, &p);
   std::vector<int> c(p), d(p);
   MPI_Allgather(&nmine, 1, MPI_INT, c.data(), 1, MPI_INT, comm);
   int tot = 0;
   for (int r = 0; r < p; r++) { d[r] = tot; tot += c[r]; }
   std::vector<real_t> all(tot);
   MPI_Allgatherv(mine, nmine, MPITypeMap<real_t>::mpi_type, all.data(), c.data(), d.data(),
                  MPITypeMap<real_t>::mpi_type, comm);
   return all;
}

// ---- the TS driver -------------------------------------------------------------

struct RunSpec
{
   double t0 = 0.0, h0 = 1.0, tmax = 1.0;
   PetscInt max_steps = 1000000;
   bool adaptive = true;
   /// After every accepted step, on every rank: the TS, the time and the
   /// gathered global state (3 N).
   std::function<void(TS, QDImplicitFault &, double, const std::vector<real_t> &)> post;
};

struct RunOut
{
   TSConvergedReason reason = TS_CONVERGED_ITERATING;
   PetscErrorCode err = PETSC_SUCCESS;
   PetscInt steps = 0, rejections = 0, snes_fails = 0, snes_its = 0;
   double t_end = 0.0;
   std::vector<real_t> final_all;
   std::vector<double> times;
   QDImplicitFault::Counts counts;
};

struct Ctx
{
   QDImplicitFault *f = nullptr;
   const RunSpec *spec = nullptr;
   std::vector<double> times;
};

PetscErrorCode PostStep(TS ts)
{
   PetscFunctionBeginUser;
   void *c = nullptr;
   PetscCall(TSGetApplicationContext(ts, &c));
   Ctx *ctx = static_cast<Ctx *>(c);
   PetscReal t;
   Vec U;
   PetscCall(TSGetTime(ts, &t));
   PetscCall(TSGetSolution(ts, &U));
   const PetscScalar *u;
   PetscInt n;
   PetscCall(VecGetLocalSize(U, &n));
   PetscCall(VecGetArrayRead(U, &u));
   const std::vector<real_t> all = GatherAll(ctx->f->Comm(), reinterpret_cast<const real_t *>(u),
                                             static_cast<int>(n));
   PetscCall(VecRestoreArrayRead(U, &u));
   ctx->times.push_back(static_cast<double>(t));
   if (ctx->spec->post) { ctx->spec->post(ts, *ctx->f, static_cast<double>(t), all); }
   PetscFunctionReturn(PETSC_SUCCESS);
}

#define TSCHK(call) MFEM_VERIFY((call) == PETSC_SUCCESS, #call)

RunOut Run(QDImplicitFault &f, const std::vector<real_t> &U0own, const RunSpec &spec)
{
   RunOut r;
   MPI_Comm comm = f.Comm();
   const int nl = QDImplicitFault::kPer * f.NumOwned();
   MFEM_VERIFY(U0own.size() == size_t(nl), "Run: U0 has " << U0own.size() << " entries");
   Vec U;
   TS ts;
   TSCHK(VecCreateMPI(comm, nl, PETSC_DETERMINE, &U));
   PetscScalar *u;
   TSCHK(VecGetArray(U, &u));
   for (int i = 0; i < nl; i++) { u[i] = U0own[i]; }
   TSCHK(VecRestoreArray(U, &u));
   TSCHK(TSCreate(comm, &ts));
   TSCHK(f.SetupTS(ts));
   TSCHK(f.SetTolerances(ts));
   TSCHK(TSSetTime(ts, spec.t0));
   TSCHK(TSSetTimeStep(ts, spec.h0));
   TSCHK(TSSetMaxTime(ts, spec.tmax));
   TSCHK(TSSetMaxSteps(ts, spec.max_steps));
   TSCHK(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
   TSCHK(TSSetMaxStepRejections(ts, -1));
   TSAdapt adapt;
   TSCHK(TSGetAdapt(ts, &adapt));
   TSCHK(TSAdaptSetType(adapt, spec.adaptive ? TSADAPTBASIC : TSADAPTNONE));
   TSCHK(TSSetFromOptions(ts));
   Ctx ctx;
   ctx.f = &f;
   ctx.spec = &spec;
   TSCHK(TSSetApplicationContext(ts, &ctx));
   TSCHK(TSSetPostStep(ts, PostStep));
   TSCHK(TSSetSolution(ts, U));
   TSCHK(TSSetUp(ts));
   f.ResetCounts();
   PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
   r.err = TSSolve(ts, U);
   PetscPopErrorHandler();
   PetscReal t;
   TSGetConvergedReason(ts, &r.reason);
   TSGetStepNumber(ts, &r.steps);
   TSGetStepRejections(ts, &r.rejections);
   TSGetSNESFailures(ts, &r.snes_fails);
   TSGetSNESIterations(ts, &r.snes_its);
   TSGetTime(ts, &t);
   r.t_end = static_cast<double>(t);
   const PetscScalar *uu;
   VecGetArrayRead(U, &uu);
   r.final_all = GatherAll(comm, reinterpret_cast<const real_t *>(uu), nl);
   VecRestoreArrayRead(U, &uu);
   r.times = ctx.times;
   r.counts = f.GetCounts();
   TSDestroy(&ts);
   VecDestroy(&U);
   return r;
}

std::string Summary(const RunOut &r)
{
   std::ostringstream o;
   o << "reason " << TSConvergedReasons[r.reason] << ", t " << std::setprecision(8) << r.t_end
     << ", steps " << r.steps << ", rejections " << r.rejections << ", SNES failures "
     << r.snes_fails << ", Newton its " << r.snes_its << ", factors " << r.counts.factors
     << ", Jacobians " << r.counts.jacobians << ", residuals " << r.counts.residuals
     << ", clipped its " << r.counts.clipped_its << ", Newton fails " << r.counts.newton_fails
     << ", branch changes " << r.counts.branch_changes << ", gap columns " << r.counts.columns
     << ", most open " << r.counts.max_open << ", plain updates " << r.counts.plain_applies;
   return o.str();
}

/// Owned part of a global state.
std::vector<real_t> Owned(const QDImplicitFault &f, const std::vector<real_t> &all)
{
   const int k = QDImplicitFault::kPer;
   return std::vector<real_t>(all.begin() + k * f.FirstOwned(),
                              all.begin() + k * (f.FirstOwned() + f.NumOwned()));
}

// ---- tableaus (PETSc arkimex.c, implicit part; stiffly accurate: b = last row)

struct Tableau
{
   std::string name;
   int s;
   std::vector<double> A;   // s x s row-major
};

Tableau ARK3()
{
   const double g = 1767732205903.0 / 4055673282236.0;
   return {"3", 4,
   {
      0, 0, 0, 0,
      g, g, 0, 0,
      2746238789719.0 / 10658868560708.0, -640167445237.0 / 6845629431997.0, g, 0,
      1471266399579.0 / 7840856788654.0, -4482444167858.0 / 7529755066697.0,
      11266239266428.0 / 11593286722821.0, g
   }};
}

Tableau ARK4()
{
   return {"4", 6,
   {
      0, 0, 0, 0, 0, 0,
      1. / 4, 1. / 4, 0, 0, 0, 0,
      8611. / 62500, -1743. / 31250, 1. / 4, 0, 0, 0,
      5012029. / 34652500, -654441. / 2922500, 174375. / 388108, 1. / 4, 0, 0,
      15267082809. / 155376265600, -71443401. / 120774400, 730878875. / 902184768,
      2285395. / 8070912, 1. / 4, 0,
      82889. / 524892, 0, 15625. / 83664, 69875. / 102672, -2260. / 8211, 1. / 4
   }};
}

/// R(z) of a stiffly accurate tableau: the last stage of Y = 1 + z A Y.
double Stability(const Tableau &T, double z)
{
   std::vector<double> Y(T.s);
   for (int i = 0; i < T.s; i++)
   {
      double num = 1.0;
      for (int j = 0; j < i; j++) { num += z * T.A[size_t(i) * T.s + j] * Y[j]; }
      Y[i] = num / (1.0 - z * T.A[size_t(i) * T.s + i]);
   }
   return Y[T.s - 1];
}

bool AllFinite(const std::vector<real_t> &v)
{
   for (real_t x : v) { if (!std::isfinite(x)) { return false; } }
   return true;
}

// ---- T0a / T0b -------------------------------------------------------------------

/// The 6-point state of T0a: point 1 open (sigma < gam g), point 4 closed
/// with g > 0, the rest closed with g = 0; every slip rate clearly nonzero.
struct T0State
{
   QDImplicitFault::Problem P;
   std::vector<real_t> U, Ud;   // global
};

T0State MakeT0(MPI_Comm comm, double floor, bool taper)
{
   const int N = 6;
   const double k = 5e7, gam = 1e8;
   auto s0 = [](int j) { return (j == 1) ? -3e5 : 20e6; };
   auto tp = [&](int j) { return -0.62 * std::abs(s0(j)) - 1e4 * j; };
   T0State S{Chain(comm, N, k, gam, tp, s0, 1e-9), {}, {}};
   S.P.sigma_floor = floor;
   S.P.floor_taper = taper;
   S.U.resize(3 * N);
   S.Ud.resize(3 * N);
   for (int j = 0; j < N; j++)
   {
      S.U[3 * j] = 1e-3 * (j + 1);
      S.U[3 * j + 1] = 0.6 + 0.01 * j;
      S.U[3 * j + 2] = (j == 1) ? 2e-3 : ((j == 4) ? 1e-3 : 0.0);
      S.Ud[3 * j] = 1e-6 * std::sin(0.4 * j + 0.1);
      S.Ud[3 * j + 1] = 1e-4 * std::cos(0.7 * j);
      S.Ud[3 * j + 2] = 1e-5 * std::sin(1.1 * j + 0.3);
   }
   return S;
}

void TestT0ab(MPI_Comm comm)
{
   Say("T0a/T0b: Newton matrix and preconditioner, 6-point chain");
   struct Mode { double floor; bool taper; const char *name; };
   for (const Mode m : {Mode{0.0, false, "no floor"}, Mode{10e6, true, "floor 10 MPa, taper"},
                        Mode{10e6, false, "hard floor 10 MPa"}})
   {
      for (const double shift : {2.3e-3, 10.0})
      for (const bool local : {false, true})
      {
         // local: Newton's own residual G (each point's stage solve) and
         // its matrix; else F and dF/dU + shift dF/dUdot
         T0State S = MakeT0(comm, m.floor, m.taper);
         QDImplicitSettings set;
         set.newton_local = local;
         QDImplicitFault f(comm, S.P, set);
         f.EnsureColumn(1);
         f.EnsureColumn(4);
         const int N = f.NumGlobal(), n = f.NumOwned(), lo = f.FirstOwned();
         auto own = [&](const std::vector<real_t> &all)
         {
            return std::vector<real_t>(all.begin() + 3 * lo, all.begin() + 3 * (lo + n));
         };
         const std::vector<real_t> U = own(S.U), Ud = own(S.Ud);
         // branch of each point
         std::vector<real_t> V, tau, sig;
         std::vector<char> open;
         f.Evaluate(0.0, U.data(), V, tau, sig, open);
         std::vector<real_t> openr(open.begin(), open.end());
         const std::vector<real_t> open_all = GatherAll(comm, openr.data(), n);
         const std::vector<real_t> V_all = GatherAll(comm, V.data(), n);
         const bool branches = open_all[1] == 1.0 && open_all[4] == 0.0 && open_all[0] == 0.0;
         // direction: x_g only where a column is cached (points 1 and 4)
         std::vector<real_t> xall(3 * N);
         for (int j = 0; j < N; j++)
         {
            xall[3 * j] = 1e-4 * std::sin(1.3 * j + 0.7);
            xall[3 * j + 1] = 1e-3 * std::cos(0.9 * j + 0.2);
            xall[3 * j + 2] = (j == 1 || j == 4) ? 1e-4 * (1.0 + 0.5 * j) : 0.0;
         }
         const std::vector<real_t> x = own(xall);
         f.Jacobian(0.0, U.data(), Ud.data(), shift);
         std::vector<real_t> Jx(3 * n), Fp(3 * n), Fm(3 * n), Up(3 * n), Um(3 * n), Udp(3 * n),
             Udm(3 * n);
         if (local) { f.ApplyLocalJacobian(x.data(), Jx.data()); }
         else { f.ApplyJacobian(x.data(), Jx.data()); }
         const double eps = 1e-4;
         for (int i = 0; i < 3 * n; i++)
         {
            Up[i] = U[i] + eps * x[i];
            Um[i] = U[i] - eps * x[i];
            Udp[i] = Ud[i] + shift * eps * x[i];
            Udm[i] = Ud[i] - shift * eps * x[i];
         }
         if (local)
         {
            f.LocalResidual(0.0, Up.data(), Udp.data(), shift, Fp.data());
            f.LocalResidual(0.0, Um.data(), Udm.data(), shift, Fm.data());
         }
         else
         {
            f.Residual(0.0, Up.data(), Udp.data(), Fp.data());
            f.Residual(0.0, Um.data(), Udm.data(), Fm.data());
         }
         double num = 0.0, den = 0.0;
         for (int i = 0; i < 3 * n; i++)
         {
            const double fd = (Fp[i] - Fm[i]) / (2.0 * eps);
            num = std::max(num, std::abs(fd - Jx[i]));
            den = std::max(den, std::abs(Jx[i]));
         }
         double nd[2] = {num, den};
         MPI_Allreduce(MPI_IN_PLACE, nd, 2, MPI_DOUBLE, MPI_MAX, comm);
         const double rel = nd[0] / nd[1];
         // PC: J (P b) = b, b_g only at points with a column
         std::vector<real_t> ball(3 * N);
         for (int j = 0; j < N; j++)
         {
            ball[3 * j] = 1e-6 * std::cos(0.3 * j + 1.0);
            ball[3 * j + 1] = 1e-5 * std::sin(0.8 * j);
            ball[3 * j + 2] = (j == 1 || j == 4) ? 3e3 * (j - 2.0) : 0.0;
         }
         const std::vector<real_t> b = own(ball);
         std::vector<real_t> y(3 * n), Jy(3 * n);
         if (local)
         {
            f.ApplyLocalInverse(b.data(), y.data());
            f.ApplyLocalJacobian(y.data(), Jy.data());
         }
         else
         {
            f.ApplyInverse(b.data(), y.data());
            f.ApplyJacobian(y.data(), Jy.data());
         }
         double pn = 0.0, pd = 0.0;
         for (int i = 0; i < 3 * n; i++)
         {
            pn = std::max(pn, std::abs(Jy[i] - b[i]));
            pd = std::max(pd, std::abs(b[i]));
         }
         double pnd[2] = {pn, pd};
         MPI_Allreduce(MPI_IN_PLACE, pnd, 2, MPI_DOUBLE, MPI_MAX, comm);
         const double prel = pnd[0] / pnd[1];
         std::ostringstream w;
         w << m.name << ", shift " << shift << (local ? ", G (point stage solves)" : ", F");
         Say("    " + w.str() + ": V = " + Fmt(V_all[0], 3) + " .. " + Fmt(V_all[5], 3) +
             " m/s, |J x - FD| / |J x| = " + Fmt(rel, 3) + ", |J P b - b| / |b| = " +
             Fmt(prel, 3));
         Check(branches, "T0a " + w.str() + ": point 1 open, 4 closed with g > 0, 0 closed");
         Check(rel < 1e-6, "T0a " + w.str() + ": Newton matrix = central differences to 1e-6");
         Check(prel < 1e-10, "T0b " + w.str() + ": preconditioner inverts it to 1e-10");
      }
   }
   // slope mode
   T0State S = MakeT0(comm, 0.0, false);
   QDImplicitFault f(comm, S.P, QDImplicitSettings());
   f.EnsureColumn(1);
   f.EnsureColumn(4);
   const int n = f.NumOwned(), lo = f.FirstOwned();
   std::vector<real_t> U(S.U.begin() + 3 * lo, S.U.begin() + 3 * (lo + n));
   std::vector<real_t> Ud(S.Ud.begin() + 3 * lo, S.Ud.begin() + 3 * (lo + n));
   f.Jacobian(0.0, U.data(), Ud.data(), 1.0 / std::numeric_limits<double>::epsilon());
   std::vector<real_t> b(3 * n), y(3 * n);
   for (int i = 0; i < 3 * n; i++) { b[i] = 1.0 + i; }
   f.ApplyInverse(b.data(), y.data());
   bool ok = true;
   for (int o = 0; o < n; o++)
   {
      ok = ok && y[3 * o] == b[3 * o] && y[3 * o + 1] == b[3 * o + 1] && y[3 * o + 2] == 0.0;
   }
   Check(ok && f.GetCounts().factors == 0,
         "T0b slope mode (shift 1/eps): identity on slip and psi, 0 on the gaps, no factor");
}

// ---- T0d ---------------------------------------------------------------------------

void TestT0d(MPI_Comm comm)
{
   Say("T0d: consistent start, 6-point chain, points 1, 2 and 4 in tension");
   const int N = 6;
   const double k = 5e7, gam = 1e8;
   auto s0 = [](int j) { return (j == 1) ? -4e5 : ((j == 2) ? -1e5 : ((j == 4) ? -2e5 : 5e6)); };
   auto tp = [](int) { return -1e6; };
   QDImplicitFault f(comm, Chain(comm, N, k, gam, tp, s0, 0.0), QDImplicitSettings());
   const int n = f.NumOwned();
   std::vector<real_t> U(3 * n, 0.0);
   for (int o = 0; o < n; o++) { U[3 * o] = 1e-3; U[3 * o + 1] = 0.6; }
   const long opened = f.ConsistentStart(0.0, U.data());
   std::vector<real_t> V, tau, sig;
   std::vector<char> open;
   f.Evaluate(0.0, U.data(), V, tau, sig, open);
   std::vector<real_t> g(n);
   for (int o = 0; o < n; o++) { g[o] = U[3 * o + 2]; }
   const std::vector<real_t> g_all = GatherAll(comm, g.data(), n);
   const std::vector<real_t> s_all = GatherAll(comm, sig.data(), n);
   double min_sig = 1e30, min_g = 1e30, max_gs = 0.0;
   std::string line;
   for (int j = 0; j < N; j++)
   {
      min_sig = std::min(min_sig, static_cast<double>(s_all[j]));
      min_g = std::min(min_g, static_cast<double>(g_all[j]));
      max_gs = std::max(max_gs, std::abs(static_cast<double>(g_all[j] * s_all[j])));
      line += " (" + Fmt(g_all[j], 4) + " m, " + Fmt(s_all[j], 4) + " Pa)";
   }
   bool cols = true;
   for (int j = 0; j < N; j++) { cols = cols && (f.HasColumn(j) == (g_all[j] > 0.0)); }
   Say("    open points " + std::to_string(opened) + "; (g, sigma):" + line);
   Check(opened >= 1 && min_g >= 0.0 && min_sig >= -1e-6 && max_gs <= 1e-6 * 5e6 * 1e-2,
         "T0d g >= 0, sigma >= 0, g sigma = 0 at every point");
   Check(cols, "T0d a gap column for exactly the points with g > 0");
}

// ---- T1a ---------------------------------------------------------------------------

void TestT1a(MPI_Comm comm)
{
   Say("T1a: one point at zero strength, delta' = -k delta / eta, fixed steps");
   const double k = 1e6;
   for (const Tableau &T : {ARK3(), ARK4()})
   {
      double worst = 0.0;
      bool finite = true;
      for (const double h : {1e-2, 1.0, 1e2, 1e4, 1e5, 1e7})
      {
         QDImplicitSettings set;
         set.esdirk_type = T.name;
         QDImplicitFault f(comm, OnePoint(comm, k, 0.0, 0.0, 0.0), set);
         std::vector<real_t> U0;
         if (f.NumOwned() == 1) { U0 = {1.0, PsiSS(1e-9), 0.0}; }
         std::vector<double> delta;
         RunSpec spec;
         spec.h0 = h;
         spec.tmax = 1e30;
         spec.max_steps = 3;
         spec.adaptive = false;
         spec.post = [&](TS, QDImplicitFault &, double, const std::vector<real_t> &all)
         {
            delta.push_back(all[0]);
         };
         const RunOut r = Run(f, U0, spec);
         const double z = -h * k / ph.eta;
         const double R = Stability(T, z);
         double prev = 1.0, err = 0.0;
         for (double d : delta)
         {
            err = std::max(err, std::abs(d / prev - R) / std::max(1.0, std::abs(R)));
            prev = d;
         }
         finite = finite && AllFinite(r.final_all) && delta.size() == 3;
         worst = std::max(worst, err);
         Say("    ARK" + T.name + " h " + Fmt(h, 3) + " s: z = " + Fmt(z, 4) + ", R(z) = " +
             Fmt(R, 10) + ", step ratios " + Fmt(delta.size() > 0 ? delta[0] : NAN, 10) +
             ", max |ratio - R| " + Fmt(err, 3) + ", Newton its " + std::to_string(r.snes_its));
      }
      Check(finite, "T1a ARK" + T.name + ": 3 steps at every h, finite");
      Check(worst < 1e-12, "T1a ARK" + T.name + ": every step ratio = R(-h k / eta) to 1e-12");
   }
}

// ---- T1b ---------------------------------------------------------------------------

/// The spring slider (delta, psi) by classical RK4, n steps to T.
void RK4Slider(double k, double tau_pre, double TDt, double psi0, double T, long n,
               double &delta, double &psi)
{
   auto rhs = [&](double t, double d, double p, double &dd, double &dp)
   {
      const double V = SlipRate(tau_pre + k * d + TDt * t, p, ph.sigma);
      dd = V;
      dp = AgingRate(V, p);
   };
   const double h = T / n;
   double d = 0.0, p = psi0;
   for (long i = 0; i < n; i++)
   {
      const double t = i * h;
      double k1d, k1p, k2d, k2p, k3d, k3p, k4d, k4p;
      rhs(t, d, p, k1d, k1p);
      rhs(t + 0.5 * h, d + 0.5 * h * k1d, p + 0.5 * h * k1p, k2d, k2p);
      rhs(t + 0.5 * h, d + 0.5 * h * k2d, p + 0.5 * h * k2p, k3d, k3p);
      rhs(t + h, d + h * k3d, p + h * k3p, k4d, k4p);
      d += h / 6.0 * (k1d + 2 * k2d + 2 * k3d + k4d);
      p += h / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p);
   }
   delta = d;
   psi = p;
}

void TestT1b(MPI_Comm comm)
{
   Say("T1b: stable spring slider (k = 2 k_c), fixed steps halved 5 times vs RK4");
   const double k = 2.0 * ph.kc();
   const double psi0 = PsiSS(ph.Vpl);
   const double Vstart = 10.0 * ph.Vpl;
   const double tau_pre = -(ph.sigma * Fric(Vstart, psi0) + ph.eta * Vstart);
   const double TDt = -k * ph.Vpl;
   const double T = 1e7;
   double dref, pref, dref2, pref2;
   RK4Slider(k, tau_pre, TDt, psi0, T, 40960, dref, pref);
   RK4Slider(k, tau_pre, TDt, psi0, T, 81920, dref2, pref2);
   Say("    RK4 reference: delta(T) = " + Fmt(dref, 15) + " m, psi(T) = " + Fmt(pref, 15) +
       " (halved step changes them by " + Fmt(std::abs(dref - dref2), 2) + ", " +
       Fmt(std::abs(pref - pref2), 2) + ")");
   for (const Tableau &Tb : {ARK3(), ARK4()})
   {
      std::vector<double> err;
      for (int m = 0; m <= 5; m++)
      {
         const long nsteps = 10L << m;
         QDImplicitSettings set;
         set.esdirk_type = Tb.name;
         set.atol = 1e-15;
         set.rtol = 1e-12;
         set.newton_tol = 1e-2;
         QDImplicitFault f(comm, OnePoint(comm, k, tau_pre, ph.sigma, TDt), set);
         std::vector<real_t> U0;
         if (f.NumOwned() == 1) { U0 = {0.0, psi0, 0.0}; }
         RunSpec spec;
         spec.h0 = T / nsteps;
         spec.tmax = T;
         spec.adaptive = false;
         const RunOut r = Run(f, U0, spec);
         const double e = std::max(std::abs(r.final_all[0] - dref),
                                   std::abs(r.final_all[1] - pref) * ph.Dc / ph.b);
         err.push_back(e);
         Say("    ARK" + Tb.name + " " + std::to_string(nsteps) + " steps: t " +
             Fmt(r.t_end, 10) + ", delta " + Fmt(r.final_all[0], 12) + ", psi " +
             Fmt(r.final_all[1], 12) + ", error " + Fmt(e, 3) + (m ? ", order " +
             Fmt(std::log2(err[m - 1] / e), 3) : std::string()) + ", Newton its " +
             std::to_string(r.snes_its));
      }
      double pmin = 1e9;
      for (int m = 3; m <= 5; m++) { pmin = std::min(pmin, std::log2(err[m - 1] / err[m])); }
      const double want = (Tb.name == "3") ? 2.8 : 3.6;
      Check(pmin >= want, "T1b ARK" + Tb.name + ": observed order over the 3 finest pairs " +
            Fmt(pmin, 3) + " >= " + Fmt(want, 2));
   }
}

// ---- T1c ---------------------------------------------------------------------------

void TestT1c(MPI_Comm comm)
{
   Say("T1c: one point pulled into tension and back (sigma0 10 MPa - 20 MPa sin(pi t/T))");
   const double k = 5e7, gam = 2.0 * k, T = 1e6;
   const double tau_pre = -0.5e6, sigma0 = 10e6;
   QDImplicitFault::Problem P = OnePoint(comm, k, tau_pre, sigma0, 0.0, 0.1 * k, gam);
   P.extra_load = [=](real_t t, int, real_t &dt, real_t &ds)
   {
      dt = 0.0;
      ds = -20e6 * std::sin(M_PI * t / T);
   };
   QDImplicitSettings set;
   set.newton_tol = newton_tol;
   QDImplicitFault f(comm, P, set);
   std::vector<real_t> U0;
   if (f.NumOwned() == 1) { U0 = {0.0, 0.6, 0.0}; }
   double min_sig = 1e30, max_gsig = 0.0, min_g = 1e30, max_g = 0.0, max_dpsi_open = 0.0;
   double max_gerr = 0.0;
   long open_steps = 0;
   double prev_psi = 0.0, prev_g = 0.0, t_open = -1.0, t_close = -1.0;
   RunSpec spec;
   spec.h0 = 1e2;
   spec.tmax = T;
   spec.post = [&](TS, QDImplicitFault &ff, double t, const std::vector<real_t> &all)
   {
      std::vector<real_t> V, tau, sig;
      std::vector<char> open;
      const std::vector<real_t> own = Owned(ff, all);
      ff.Evaluate(t, own.data(), V, tau, sig, open);
      const std::vector<real_t> s_all = GatherAll(ff.Comm(), sig.data(), ff.NumOwned());
      const double s = s_all[0], g = all[2], psi = all[1];
      min_sig = std::min(min_sig, s);
      min_g = std::min(min_g, g);
      max_g = std::max(max_g, g);
      max_gsig = std::max(max_gsig, std::abs(g * s));
      // one point, Gs = 0: sigma = sigma0 + load + gam g, so g is known
      const double g_exact = std::max(0.0, -(sigma0 - 20e6 * std::sin(M_PI * t / T))) / gam;
      max_gerr = std::max(max_gerr, std::abs(g - g_exact));
      if (g > 0.0 && t_open < 0.0) { t_open = t; }
      if (t_open >= 0.0 && g == 0.0 && t_close < 0.0) { t_close = t; }
      if (g > 0.0 && prev_g > 0.0)
      {
         open_steps++;
         max_dpsi_open = std::max(max_dpsi_open, std::abs(psi - prev_psi));
      }
      prev_psi = psi;
      prev_g = g;
   };
   const RunOut r = Run(f, U0, spec);
   Say("    " + Summary(r));
   Say("    opened at t = " + Fmt(t_open, 6) + " s (sigma0 + load = 0 at T/6 = " +
       Fmt(T / 6, 6) + "), closed at " + Fmt(t_close, 6) + " s (5T/6 = " + Fmt(5 * T / 6, 6) +
       "), max g " + Fmt(max_g, 6) + " m (" + Fmt(10e6 / gam, 6) + " expected), min sigma " +
       Fmt(min_sig, 3) + " Pa, max |g sigma| " + Fmt(max_gsig, 3) + " N/m, min g " +
       Fmt(min_g, 3) + ", steps open " + std::to_string(open_steps) + ", max |dpsi| open " +
       Fmt(max_dpsi_open, 3) + ", max |g - g_exact| " + Fmt(max_gerr, 3) + " m");
   Check(r.reason == TS_CONVERGED_TIME && AllFinite(r.final_all), "T1c reaches T, finite");
   Check(min_sig >= -1e-9 * sigma0, "T1c sigma >= -1e-9 sigma0 at every accepted state");
   Check(min_g >= -1e-15, "T1c g >= 0 at every accepted state");
   Check(max_gsig <= 1e-9 * sigma0 * max_g, "T1c |g sigma| <= 1e-9 sigma0 max g");
   Check(open_steps > 0 && max_dpsi_open <= 1e-12, "T1c psi frozen while open");
   Check(max_gerr <= 1e-9 * 10e6 / gam && max_g > 0.9 * 10e6 / gam,
         "T1c g = -(sigma0 + load) / gamma at every accepted state, to 1e-9 of its peak");
   Check(r.counts.columns == 1, "T1c one gap column fetched");
   Check(r.final_all[2] == 0.0 && t_close > 0.0, "T1c closed again, g = 0 at T");
}

// ---- T1d ---------------------------------------------------------------------------

struct RefCtx
{
   double k, tau_pre, TDt, vth, t_event = -1.0, vmax = 0.0;
};

PetscErrorCode RefRHS(TS, PetscReal t, Vec U, Vec F, void *c)
{
   PetscFunctionBeginUser;
   RefCtx *ctx = static_cast<RefCtx *>(c);
   const PetscScalar *u;
   PetscScalar *f;
   PetscCall(VecGetArrayRead(U, &u));
   PetscCall(VecGetArray(F, &f));
   const double V = SlipRate(ctx->tau_pre + ctx->k * u[0] + ctx->TDt * t, u[1], ph.sigma);
   f[0] = V;
   f[1] = AgingRate(V, u[1]);
   PetscCall(VecRestoreArray(F, &f));
   PetscCall(VecRestoreArrayRead(U, &u));
   PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RefPost(TS ts)
{
   PetscFunctionBeginUser;
   void *c;
   PetscCall(TSGetApplicationContext(ts, &c));
   RefCtx *ctx = static_cast<RefCtx *>(c);
   PetscReal t;
   Vec U;
   PetscCall(TSGetTime(ts, &t));
   PetscCall(TSGetSolution(ts, &U));
   const PetscScalar *u;
   PetscCall(VecGetArrayRead(U, &u));
   const double V = SlipRate(ctx->tau_pre + ctx->k * u[0] + ctx->TDt * t, u[1], ph.sigma);
   PetscCall(VecRestoreArrayRead(U, &u));
   ctx->vmax = std::max(ctx->vmax, std::abs(V));
   if (ctx->t_event < 0.0 && std::abs(V) >= ctx->vth)
   {
      ctx->t_event = t;
      PetscCall(TSSetConvergedReason(ts, TS_CONVERGED_USER));
   }
   PetscFunctionReturn(PETSC_SUCCESS);
}

void TestT1d(MPI_Comm comm)
{
   Say("T1d: unstable spring slider (k = 0.5 k_c) through one event, adaptive");
   const double k = 0.5 * ph.kc();
   const double psi0 = PsiSS(ph.Vpl);
   const double Vstart = 2.0 * ph.Vpl;
   const double tau_pre = -(ph.sigma * Fric(Vstart, psi0) + ph.eta * Vstart);
   const double TDt = -k * ph.Vpl, vth = 1e-3, tmax = 1e11;
   // reference: RK 5dp on (delta, psi), atol 1e-11, every rank alone
   RefCtx rc{k, tau_pre, TDt, vth};
   {
      Vec U;
      TS ts;
      TSCHK(VecCreateSeq(PETSC_COMM_SELF, 2, &U));
      PetscScalar *u;
      TSCHK(VecGetArray(U, &u));
      u[0] = 0.0;
      u[1] = psi0;
      TSCHK(VecRestoreArray(U, &u));
      TSCHK(TSCreate(PETSC_COMM_SELF, &ts));
      TSCHK(TSSetType(ts, TSRK));
      TSCHK(TSRKSetType(ts, TSRK5DP));
      TSCHK(TSSetRHSFunction(ts, nullptr, RefRHS, &rc));
      TSCHK(TSSetTimeStep(ts, 1e3));
      TSCHK(TSSetMaxTime(ts, tmax));
      TSCHK(TSSetMaxSteps(ts, 10000000));
      TSCHK(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
      TSCHK(TSSetTolerances(ts, 1e-11, nullptr, 0.0, nullptr));
      TSCHK(TSSetMaxStepRejections(ts, -1));
      TSCHK(TSSetFromOptions(ts));
      TSCHK(TSSetApplicationContext(ts, &rc));
      TSCHK(TSSetPostStep(ts, RefPost));
      TSCHK(TSSetSolution(ts, U));
      TSCHK(TSSolve(ts, U));
      PetscInt steps;
      TSGetStepNumber(ts, &steps);
      Say("    RK 5dp reference (atol 1e-11): |V| >= 1e-3 m/s first at t = " +
          Fmt(rc.t_event, 12) + " s, " + std::to_string(steps) + " steps");
      TSDestroy(&ts);
      VecDestroy(&U);
   }
   for (const char *type : {"3", "4"})
   {
      QDImplicitSettings set;
      set.esdirk_type = type;
      set.newton_tol = newton_tol;
      QDImplicitFault f(comm, OnePoint(comm, k, tau_pre, ph.sigma, TDt), set);
      std::vector<real_t> U0;
      if (f.NumOwned() == 1) { U0 = {0.0, psi0, 0.0}; }
      double t_event = -1.0;
      RunSpec spec;
      spec.h0 = 1e3;
      spec.tmax = tmax;
      spec.post = [&](TS ts, QDImplicitFault &ff, double t, const std::vector<real_t> &all)
      {
         std::vector<real_t> V, tau, sig;
         std::vector<char> open;
         const std::vector<real_t> own = Owned(ff, all);
         ff.Evaluate(t, own.data(), V, tau, sig, open);
         const std::vector<real_t> V_all = GatherAll(ff.Comm(), V.data(), ff.NumOwned());
         if (t_event < 0.0 && std::abs(V_all[0]) >= vth)
         {
            t_event = t;
            TSSetConvergedReason(ts, TS_CONVERGED_USER);
         }
      };
      const RunOut r = Run(f, U0, spec);
      const double rel = std::abs(t_event - rc.t_event) / rc.t_event;
      Say("    ARK" + std::string(type) + " (atol 1e-7): event at t = " + Fmt(t_event, 12) +
          " s, relative difference " + Fmt(rel, 3) + "; " + Summary(r));
      Check(r.reason == TS_CONVERGED_USER && t_event > 0.0,
            "T1d ARK" + std::string(type) + " reaches the event");
      Check(rel < 0.01, "T1d ARK" + std::string(type) + " event time within 1% of RK 5dp");
   }
}

// ---- T0c ---------------------------------------------------------------------------

RunOut ChainRun(MPI_Comm comm, std::vector<real_t> *open_counts)
{
   const int N = 6;
   const double k = 0.5 * ph.kc(), gam = 2.0 * k, T = 3e8;
   auto s0 = [](int) { return ph.sigma; };
   auto tp = [&](int j)
   {
      const double V = ph.Vpl * (1.0 + 0.5 * j);
      return -(ph.sigma * Fric(V, PsiSS(ph.Vpl)) + ph.eta * V);
   };
   QDImplicitFault::Problem P = Chain(comm, N, k, gam, tp, s0, ph.Vpl);
   // points 0 and 1 are pulled into tension between T/3 and 2T/3
   P.extra_load = [=](real_t t, int j, real_t &dt, real_t &ds)
   {
      dt = 0.0;
      ds = (j < 2) ? -1.2 * ph.sigma * std::max(0.0, std::sin(M_PI * (3.0 * t / T - 1.0)))
           : 0.0;
   };
   QDImplicitSettings set;
   set.newton_tol = newton_tol;
   QDImplicitFault f(comm, P, set);
   std::vector<real_t> U0(3 * f.NumOwned());
   for (int o = 0; o < f.NumOwned(); o++) { U0[3 * o + 1] = PsiSS(ph.Vpl); }
   long opened = 0;
   RunSpec spec;
   spec.h0 = 1e3;
   spec.tmax = T;
   spec.max_steps = 20000;
   spec.post = [&](TS, QDImplicitFault &, double, const std::vector<real_t> &all)
   {
      for (int j = 0; j < N; j++) { if (all[3 * j + 2] > 0.0) { opened++; } }
   };
   RunOut r = Run(f, U0, spec);
   if (open_counts) { open_counts->push_back(static_cast<real_t>(opened)); }
   return r;
}

void TestT0c()
{
   Say("T0c: 6-point chain with events and opening, all ranks vs one rank");
   std::vector<real_t> oc;
   const RunOut rw = ChainRun(MPI_COMM_WORLD, &oc);
   Say("    " + std::to_string(nranks) + " rank(s): " + Summary(rw) + ", point-steps open " +
       Fmt(oc[0], 8));
   // one rank: rank 0 alone, then broadcast
   int steps = 0, nfinal = 0;
   std::vector<real_t> fin;
   double tend = 0.0;
   if (rank == 0)
   {
      const RunOut rs = ChainRun(MPI_COMM_SELF, nullptr);
      std::cout << "    1 rank: " << Summary(rs) << std::endl;
      steps = static_cast<int>(rs.steps);
      fin = rs.final_all;
      nfinal = static_cast<int>(fin.size());
      tend = rs.t_end;
   }
   MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&nfinal, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&tend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   fin.resize(nfinal);
   MPI_Bcast(fin.data(), nfinal, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
   double maxdiff = 0.0;
   for (int i = 0; i < nfinal && i < static_cast<int>(rw.final_all.size()); i++)
   {
      maxdiff = std::max(maxdiff, std::abs(static_cast<double>(fin[i] - rw.final_all[i])));
   }
   Say("    max |state difference| " + Fmt(maxdiff, 3) + ", steps " + std::to_string(rw.steps) +
       " vs " + std::to_string(steps));
   Check(rw.reason == TS_CONVERGED_TIME && AllFinite(rw.final_all) && oc[0] > 0,
         "T0c the chain reaches T with points opened on the way, finite");
   Check(rw.steps == steps && rw.final_all == fin && rw.t_end == tend,
         "T0c same steps and final state bit for bit on " + std::to_string(nranks) +
         " rank(s) and on 1 rank");
}

} // namespace

// ---- T1e ---------------------------------------------------------------------------

void TestT1e(MPI_Comm comm)
{
   Say("T1e: [cycle.relaxation] opening = \"relaxed\": one point unclamped at a constant rate");
   const double k = 5e7, gam = 2.0 * k, T = 4e6, t_r = 2e5;
   const double tau_pre = -0.5e6, sigma0 = 10e6;
   for (int c = 0; c < 2; c++)
   {
      const double rho = (c == 0 ? 0.5 : 2.0) * sigma0 / t_r;      // Pa/s
      const double T1 = (c == 0) ? T : 0.25 * T;                   // the load stops at T1
      auto load = [=](double t) { return -rho * std::min(t, T1); };
      QDImplicitFault::Problem P = OnePoint(comm, k, tau_pre, sigma0, 0.0, 0.1 * k, gam);
      P.extra_load = [=](real_t t, int, real_t &dt, real_t &ds)
      {
         dt = 0.0;
         ds = load(t);
      };
      seas::cycle::NormalStressRelaxation relax(P.n, t_r);
      P.normal_reduction = [&relax](real_t t, int o) { return relax.ReductionAt(t, o); };
      QDImplicitSettings set;
      set.newton_tol = newton_tol;
      QDImplicitFault f(comm, P, set);
      std::vector<real_t> U0;
      if (f.NumOwned() == 1) { U0 = {0.0, 0.6, 0.0}; }
      double min_s = 1e30, min_g = 1e30, max_g = 0.0, max_gs = 0.0, max_serr = 0.0;
      double min_el = 1e30, t_open = -1.0, t_close = -1.0, g_T1 = -1.0, s_end = 0.0;
      double prev_g = 0.0;
      RunSpec spec;
      spec.h0 = T / 4000.0;
      spec.tmax = T;
      spec.adaptive = false;
      spec.post = [&](TS, QDImplicitFault &ff, double t, const std::vector<real_t> &all)
      {
         std::vector<real_t> V, tau, sig;
         std::vector<char> open;
         const std::vector<real_t> own = Owned(ff, all);
         ff.Evaluate(t, own.data(), V, tau, sig, open);   // sig: the relaxed stress
         // accept the step as QDPhase does: r at t, dsig the elastic change
         Vector dsig(ff.NumOwned());
         for (int o = 0; o < ff.NumOwned(); o++)
         {
            dsig(o) = sig[o] + relax.ReductionAt(t, o) - sigma0;
         }
         std::vector<real_t> el(ff.NumOwned());
         for (int o = 0; o < ff.NumOwned(); o++) { el[o] = sigma0 + dsig(o); }
         relax.Accept(t, dsig);
         const double sr = GatherAll(ff.Comm(), sig.data(), ff.NumOwned())[0];
         const double se = GatherAll(ff.Comm(), el.data(), ff.NumOwned())[0];
         const double g = all[2];
         min_s = std::min(min_s, sr);
         min_el = std::min(min_el, se);
         min_g = std::min(min_g, g);
         max_g = std::max(max_g, g);
         max_gs = std::max(max_gs, std::abs(g * sr));
         if (c == 0)
         {
            const double exact = sigma0 - rho * t_r * (1.0 - std::exp(-t / t_r));
            max_serr = std::max(max_serr, std::abs(sr - exact));
         }
         if (g > 0.0 && t_open < 0.0) { t_open = t; }
         if (t_open >= 0.0 && prev_g > 0.0 && g == 0.0 && t_close < 0.0) { t_close = t; }
         if (g_T1 < 0.0 && t >= T1 * (1.0 - 1e-12)) { g_T1 = g; }
         prev_g = g;
         s_end = sr;
      };
      const RunOut r = Run(f, U0, spec);
      Say("    case " + std::string(c == 0 ? "(a)" : "(b)") + ": " + Summary(r));
      Say("    relaxed stress min " + Fmt(min_s, 6) + " Pa, at the end " + Fmt(s_end, 6) +
          " Pa; elastic stress min " + Fmt(min_el, 6) + " Pa; max g " + Fmt(max_g, 6) +
          " m, min g " + Fmt(min_g, 3) + ", max |g s| " + Fmt(max_gs, 3) + " N/m");
      const std::string w = std::string("T1e ") + (c == 0 ? "(a)" : "(b)");
      Check(r.reason == TS_CONVERGED_TIME && AllFinite(r.final_all), w + " reaches T, finite");
      if (c == 0)
      {
         Say("    max |s - exact| " + Fmt(max_serr, 4) + " Pa (rho t_r = " + Fmt(rho * t_r, 4) +
             " Pa)");
         Check(max_g == 0.0 && r.counts.columns == 0, w + " never opens (no gap column)");
         Check(min_el < -1.9 * sigma0, w + " the elastic stress goes below -1.9 sigma0");
         Check(max_serr <= 0.01 * rho * t_r,
               w + " relaxed stress = sigma0 - rho t_r (1 - exp(-t / t_r)) to 1% of rho t_r");
      }
      else
      {
         const double ts = -t_r * std::log(1.0 - sigma0 / (rho * t_r));
         const double gx = (rho - sigma0 / t_r) * (T1 - ts) / gam;
         const double tc = T1 + gx * gam * t_r / sigma0;
         Say("    opened at " + Fmt(t_open, 6) + " s (t* = " + Fmt(ts, 6) + "), g(T1) " +
             Fmt(g_T1, 6) + " m (" + Fmt(gx, 6) + " expected), closed at " + Fmt(t_close, 6) +
             " s (" + Fmt(tc, 6) + " expected)");
         Check(min_s >= -1e-9 * sigma0 && min_g >= -1e-15 && max_gs <= 1e-9 * sigma0 * max_g,
               w + " s >= 0, g >= 0, g s = 0 at every accepted state");
         Check(std::abs(t_open - ts) <= 3.0 * spec.h0, w + " opens at t* within 3 steps");
         Check(std::abs(g_T1 - gx) <= 0.02 * gx, w + " gap at T1 within 2% of the expected one");
         Check(t_close > 0.0 && std::abs(t_close - tc) <= 0.02 * (tc - T1) + 3.0 * spec.h0,
               w + " closes again at the expected time (2%)");
         Check(r.final_all[2] == 0.0 && s_end > 0.9 * sigma0,
               w + " closed at T with the relaxed stress back above 0.9 sigma0");
      }
   }
}

int main(int argc, char **argv)
{
   Mpi::Init(argc, argv);
   MFEMInitializePetsc(&argc, &argv, nullptr, nullptr);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   Say("test_qd_esdirk_points on " + std::to_string(nranks) + " rank(s), PETSc " +
       std::to_string(PETSC_VERSION_MAJOR) + "." + std::to_string(PETSC_VERSION_MINOR) + "." +
       std::to_string(PETSC_VERSION_SUBMINOR));
   // -ts_adapt_wnormtype infinity, as ucerf2d/petsc_ts_rk45_ucerf.cfg
   MFEM_VERIFY(PetscOptionsSetValue(nullptr, "-ts_adapt_wnormtype", "infinity")
               == PETSC_SUCCESS, "options");
   PetscReal ntol = newton_tol;
   PetscBool has = PETSC_FALSE;
   MFEM_VERIFY(PetscOptionsGetReal(nullptr, nullptr, "-qd_newton_tol", &ntol, &has)
               == PETSC_SUCCESS, "-qd_newton_tol");
   newton_tol = static_cast<double>(ntol);
   Say("newton_tol " + Fmt(newton_tol, 3) + " in T1c, T1d, T0c" + (has ? " (-qd_newton_tol)" : ""));

   TestT0ab(MPI_COMM_WORLD);
   TestT0d(MPI_COMM_WORLD);
   TestT1a(MPI_COMM_WORLD);
   TestT1b(MPI_COMM_WORLD);
   TestT1c(MPI_COMM_WORLD);
   TestT1d(MPI_COMM_WORLD);
   TestT1e(MPI_COMM_WORLD);
   TestT0c();

   Say(fails ? "test_qd_esdirk_points: " + std::to_string(fails) + " FAILED" :
       "test_qd_esdirk_points: all PASSED");
   MFEMFinalizePetsc();
   Mpi::Finalize();
   return fails ? 1 : 0;
}
