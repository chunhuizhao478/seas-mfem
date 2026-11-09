// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// qd_implicit.hpp — option I of the QD step at sigma_n' -> 0 (Chunhui Zhao,
// 2026-11-09: "let's proceed with implementation, take all decisions p1-p6
// by recommended items"; plan /mnt/project-files/fully-dynamic-fault-opening/
// qd-implicit-friction/, Claude Doc "Option I: implementation and test
// plan", rev 14; note qd-implicit-friction-rev2.pdf).  Phase I1: the fault
// system, its Newton matrix and the PETSc set-up, on any dense G; the QD
// phase hands it the planar Green's function in I2.
//
// Every fault point is implicit in every stage of a stiffly accurate ESDIRK
// (PETSc arkimex, fully implicit; P2), contact is decided inside each stage
// solve (P1), an open point's slip row is the radiation row tau + eta V = 0
// (P3), the Newton matrix is factored densely on every rank (P4), on a 2-D
// planar dense G (P5), with the relaxation off (P6: r = 0).
//
// State per point j (global index, owned by one rank), 3 TS entries:
//
//    U_j = (delta_j, psi_j, g_j)      in-plane slip [m] (the kernel's
//                                     in-plane slot), state variable,
//                                     opening [m] (>= 0 at a solution)
//
// Tractions (kernel frame, as the DP5 kernel forms them; sigma > 0 is
// compression):
//
//    tau_j   = tau_pre_j + sum_l Gt_jl delta_l + sum_{l cached} Ggt_jl g_l + t TDt_j
//    sigma_j = sigma0_j  + sum_l Gs_jl delta_l + sum_{l cached} Ggs_jl g_l + t TDs_j
//
// Gt, Gs: the in-plane shear and normal rows of the planar G (in-plane slip
// columns); Ggt, Ggs: the gap columns (traction per unit opening), fetched
// once per point the first time it is on the open branch (D-B6) and kept.
//
// Residual F(t, U, Udot) = 0, the form PETSc's fully implicit arkimex solves
// (plan Section 2, with V eliminated point by point):
//
//    F_delta,j = delta'_j - V*_j(tau_j, psi_j, sigma_law,j)
//    F_psi,j   = psi'_j   - (open_j ? 0 : (b V0 / Dc_j)(exp((f0 - psi_j)/b) - |V*_j| / V0))
//    F_g,j     = min(sigma_j, gamma_j g_j)
//
// V*_j is the friction solve of the DP5 kernel (RSSolveSlipRatePsi: tau +
// sigma_law f(V, psi) + eta V = 0, V anti-parallel to tau), so the row
// delta' = V* is the plan's closed row (2.3) with V solved out, and at
// sigma_law = 0 it is V* = -tau / eta: the radiation row of an open point
// (P3, 2.4) with no branch in the slip row.  sigma_law = Lambda(sigma,
// sigma) of the kernel at r = 0 (RSLawSigmaN): max(sigma, 0) with the
// whole-law taper or without a floor, max(sigma, floor) with the hard floor.
// open_j: sigma_j < gamma_j g_j, the branch of the min (2.5, 2.9); psi is
// frozen while open (D-B3).  gamma_j = Ggs_jj, the point's own normal stress
// per unit opening ([cycle.opening] contact_scale = "auto"); it only enters
// where g_j != 0, which needs the cached column anyway.
//
// Why V is solved out instead of kept in the row (plan 2.3): the slip and
// psi rows are then Udot - H(t, U), an identity mass, and PETSc's explicit
// first stage gets the exact slope on every version we run (3.19 in the
// cloud: YdotI_1 = -F(t, U_n, 0); 3.20 and later: the slope solve F(t_n,
// U_n, Udot) = 0 at each restart, which our preconditioner answers with the
// identity on those rows and 0 on the algebraic gap rows).  The Newton
// matrix is the same as the plan's (2.7)-(2.8) up to the row scale D_j.
//
// Newton matrix (PETSc's J = dF/dU + shift dF/dUdot, shift = 1/(h gamma_d)).
// Per point, at the iterate (V = V*_j):
//
//    D_j  = sigma_law f_V + eta              (slope of the law in V)
//    fs_j = f_j dsigma_law/dsigma            (strength change with sigma)
//    p_j  = sigma_law f_psi / D_j            (dF_delta/dpsi_j)
//    c_j  = (b / Dc_j) sgn(V_j), 0 when open (psi row's pull from |V|)
//    A_j  = shift + (V0 / Dc_j) exp((f0 - psi_j)/b), shift when open
//
//    J x |_delta,j = shift x_delta,j + w_j + p_j x_psi,j
//    J x |_psi,j   = -c_j w_j + (A_j - c_j p_j) x_psi,j
//    J x |_g,j     = closed: gamma_j x_g,j;  open: [Gs x_delta + Ggs x_g]_j
//    w_j = ( [(Gt + fs Gs) x_delta]_j + [(Ggt + fs Ggs) x_g]_j ) / D_j
//
// Eliminating psi per point and the closed gap rows (x_g = b_g / gamma)
// leaves the plan's (2.8) on (delta, g_O), times D_j on the slip rows:
//
//    [ diag(Dt) + Gt + fs Gs     Ggt_:O + fs Ggs_:O ] [y_delta]   [r_delta]
//    [ Gs_O:                     Ggs_OO             ] [y_gO   ] = [r_gO   ]
//    Dt_j = D_j shift (1 - p_j c_j / A_j)
//    r_delta,j = D_j (b_delta,j - p_j (b_psi,j + c_j b_delta,j) / A_j) - closed-g terms
//    y_psi,j   = (b_psi,j + c_j b_delta,j - c_j shift y_delta,j) / A_j
//
// size N + m (m points on the open branch), assembled alike on every rank
// from the replicated G and the gathered per-point numbers, and factored by
// LU with partial pivoting (mfem::LUFactors, LAPACK where MFEM has it), so
// every rank solves the same system bit for bit.  The factor is kept while
// the shift and the branch set stay the same and Newton converges fast
// ([cycle.qd] jacobian_reuse = "step", P4), so the stages of a step share
// it (ESDIRK: one gamma_d).  The SNES Mat is a MATSHELL that is never
// multiplied (KSP preonly + PCSHELL).
//
// Newton in psi with the aging exponential kept (a line-search pre-check).
// Why: the aging rate is (b V0 / Dc) exp((f0 - psi)/b), so with a long step
// (shift = 1/(h gamma_d) small) the stage residual in psi is concave with a
// slope of about shift far above the root, and the full Newton step d_N
// lands far below it, where exp overflows (I1 test T1a: Inf after one
// iteration at h >= 100 s); from below, Newton gains less than b per
// iteration.  The linear model of the psi row at the iterate is
//    F_psi + (other terms) + (s' + kappa) d = 0,  s' = shift - c p,
//    kappa = (V0 / Dc) exp((f0 - psi)/b)  (closed points; 0 when open),
// and Newton's d_N zeroes it.  The pre-check keeps the other terms and the
// shift part linear and replaces kappa d by the exact change of the aging
// term, kappa b (1 - exp(-d / b)), so the psi change d solves
//    s' d + kappa b (1 - exp(-d / b)) = (s' + kappa) d_N
// (Brent's method; one root, between d_N and 0 for d_N < 0 and between d_N
// and d_N (s' + kappa) / s' for d_N > 0, when s' > 0).  For the aging law
// alone (p = 0) this is the exact psi row; with s' -> 0 it is d = -b ln(1
// - d_N / b), Newton in u = exp((f0 - psi)/b), where the aging rate is
// linear.  Where s' <= 0 (fast weakening at fixed traction) the change is
// log-damped beyond s = newton_psi_scale b instead: sign(d_N) s (1 +
// ln(|d_N| / s)).  Every change is capped at newton_psi_cap b (30 b).
// Changes below 1e-3 b, and so Newton's last iterations, are untouched.
//
// Newton stops when the update, in the TS error weights (atol + rtol |x|,
// max over entries), is below newton_tol, or when the error left after it,
// estimated from the contraction theta = w_k / w_(k-1) as theta / (1 -
// theta) w_k (Hairer and Wanner IV.8), is below newton_tol; and the update
// was not damped and the iterate's branch set is the one the factor was
// built with.  A solve that needs more than newton_max_it iterations fails,
// so the TS rejects the step and shrinks it (-ts_max_snes_failures -1; plan
// C').  A non-finite residual fails the solve the same way.  A contraction
// theta > newton_theta_refactor (0.2) on undamped updates marks the factor
// slow, and so does a damped update: the next Jacobian refactors.  (T1a: a
// stale factor made a clipped Newton cycle between +cap and -cap; one kept
// from the start of a long step, with psi 20 b away from the stage
// solution, contracted at theta = 0.42, so 0.5 let Newton run out of
// iterations.)
//
// Consistent start (2.10): ConsistentStart(t, U) solves the gap rows for g
// at fixed delta and psi (the open-set LCP of StageContact_, by active set).

#ifndef MFEM_SEAS_CYCLE_QD_IMPLICIT_HPP
#define MFEM_SEAS_CYCLE_QD_IMPLICIT_HPP

#include "mfem.hpp"
#include "../friction/device_rate_state.hpp"

#if defined(MFEM_USE_PETSC) && defined(MFEM_USE_MPI)
#include "petsc.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// [cycle.qd] keys of option I (plan Section 5).
struct QDImplicitSettings
{
   std::string esdirk_type = "3";     ///< PETSc arkimex type, run fully implicit
   real_t newton_tol = 1e-4;          ///< update / (atol + rtol |x|), max norm
   real_t newton_psi_scale = 1.0;     ///< psi updates beyond it (in b) are log-damped
   real_t newton_psi_cap = 30.0;      ///< |psi update| <= cap * b per iteration
   real_t newton_theta_refactor = 0.2;   ///< contraction above it: refactor
   int newton_max_it = 20;            ///< more: a failed solve, the step is rejected
   std::string jacobian_reuse = "step";   ///< "step", "stage" or "iteration"
   real_t atol = 1e-7;                ///< the TS's absolute tolerance
   real_t rtol = 0.0;                 ///< the TS's relative tolerance
};

/// The fault system of option I: data, residual, Newton matrix, and the
/// PETSc glue.  Every rank holds the whole dense G (replicated) and owns a
/// contiguous range of points, in rank order.
class QDImplicitFault
{
public:
   static constexpr int kPer = 3;   ///< TS entries per point: delta, psi, g

   /// The problem.  Global arrays are replicated on every rank; owned
   /// arrays hold this rank's points [lo, lo + n).
   struct Problem
   {
      int N = 0;                       ///< points, global
      std::vector<real_t> Gt, Gs;      ///< N x N row-major (global)
      std::vector<real_t> TDt, TDs;    ///< loading rates, N (global)
      int lo = 0, n = 0;               ///< owned range
      std::vector<real_t> tau_pre;     ///< in-plane shear pre-stress, n
      std::vector<real_t> sigma0;      ///< normal pre-stress, n (compression > 0)
      std::vector<real_t> a, eta, Dc;  ///< per owned point, n
      real_t b = 0.0, f0 = 0.0, V0 = 0.0;   ///< aging law and friction scalars
      real_t sigma_floor = 0.0;        ///< [cycle].sigma_n_floor_pa
      bool floor_taper = false;        ///< [cycle.opening].floor_taper
      /// Global point l: its gap column, Ggt(:, l) and Ggs(:, l), both of
      /// length N (global).  Collective: every rank calls it with the same l.
      std::function<void(int l, std::vector<real_t> &t, std::vector<real_t> &s)> gap_column;
      /// Optional extra loading at time t on global point j (tests only).
      std::function<void(real_t t, int j, real_t &dtau, real_t &dsig)> extra_load;
   };

   /// What happened since the last ResetCounts (this rank; the decisions
   /// are collective, so the counts agree on every rank).
   struct Counts
   {
      long residuals = 0, jacobians = 0, factors = 0, slope_solves = 0;
      long branch_changes = 0;   ///< Jacobian evaluations whose branch set moved
      long columns = 0;          ///< gap columns fetched
      long newton_its = 0, newton_solves = 0, newton_fails = 0;
      long clipped_its = 0;      ///< Newton iterations whose psi change moved > 10%
      long max_open = 0;         ///< most points on the open branch at a Jacobian
   };

   QDImplicitFault(MPI_Comm comm, Problem p, const QDImplicitSettings &s)
      : comm_(comm), P_(std::move(p)), S_(s)
   {
      MPI_Comm_rank(comm_, &rank_);
      MPI_Comm_size(comm_, &nranks_);
      const int N = P_.N, n = P_.n;
      MFEM_VERIFY(N > 0, "QDImplicitFault: no points");
      MFEM_VERIFY(P_.Gt.size() == size_t(N) * N && P_.Gs.size() == size_t(N) * N,
                  "QDImplicitFault: G is not N x N");
      MFEM_VERIFY(P_.TDt.size() == size_t(N) && P_.TDs.size() == size_t(N),
                  "QDImplicitFault: T_D is not of length N");
      for (const std::vector<real_t> *v : {&P_.tau_pre, &P_.sigma0, &P_.a, &P_.eta, &P_.Dc})
      {
         MFEM_VERIFY(v->size() == size_t(n), "QDImplicitFault: an owned array is not of "
                     "length n = " << n);
      }
      MFEM_VERIFY(P_.b > 0.0 && P_.V0 > 0.0, "QDImplicitFault: b and V0 must be > 0");
      MFEM_VERIFY(S_.jacobian_reuse == "step" || S_.jacobian_reuse == "stage" ||
                  S_.jacobian_reuse == "iteration", "QDImplicitFault: jacobian_reuse \""
                  << S_.jacobian_reuse << "\" (step, stage or iteration)");
      MFEM_VERIFY(S_.newton_tol > 0.0 && S_.newton_max_it >= 1 && S_.newton_psi_scale > 0.0 &&
                  S_.newton_psi_cap >= S_.newton_psi_scale,
                  "QDImplicitFault: newton_tol > 0, newton_max_it >= 1, newton_psi_scale > 0, "
                  "newton_psi_cap >= newton_psi_scale");
      // the owned ranges must tile [0, N) in rank order (the Green's
      // function's global index follows the partition, rank by rank)
      int mine[2] = {P_.lo, n};
      std::vector<int> all(2 * nranks_);
      MPI_Allgather(mine, 2, MPI_INT, all.data(), 2, MPI_INT, comm_);
      cnt_.resize(nranks_);
      dsp_.resize(nranks_);
      int next = 0;
      for (int r = 0; r < nranks_; r++)
      {
         MFEM_VERIFY(all[2 * r] == next, "QDImplicitFault: rank " << r << " owns from "
                     << all[2 * r] << ", expected " << next << " (contiguous, in rank order)");
         dsp_[r] = all[2 * r];
         cnt_[r] = all[2 * r + 1];
         next += cnt_[r];
      }
      MFEM_VERIFY(next == N, "QDImplicitFault: the ranks own " << next << " of " << N
                  << " points");
      prm_.V0 = P_.V0;
      prm_.aging_b = P_.b;
      prm_.aging_V0 = P_.V0;
      prm_.aging_f0 = P_.f0;
      slot_.assign(N, -1);
      delta_all_.resize(N);
      g_all_.resize(N);
      res_branch_.assign(n, 0);
   }

   QDImplicitFault(const QDImplicitFault &) = delete;
   QDImplicitFault &operator=(const QDImplicitFault &) = delete;
   ~QDImplicitFault()
   {
      if (J_) { MatDestroy(&J_); }
   }

   MPI_Comm Comm() const { return comm_; }
   int NumGlobal() const { return P_.N; }
   int NumOwned() const { return P_.n; }
   int FirstOwned() const { return P_.lo; }
   const Problem &GetProblem() const { return P_; }
   const QDImplicitSettings &Settings() const { return S_; }
   const Counts &GetCounts() const { return counts_; }
   void ResetCounts() { counts_ = Counts(); }
   /// The column of global point l is cached.
   bool HasColumn(int l) const { return slot_[l] >= 0; }
   /// gamma_l = Ggs_ll (needs the column).
   real_t Gamma(int l) const
   {
      MFEM_VERIFY(slot_[l] >= 0, "QDImplicitFault: no gap column for point " << l);
      return cols_s_[slot_[l]][l];
   }

   /// Fetch and keep the gap column of global point l if it is not cached
   /// (e.g. the points with g != 0 of a restart state).  Collective: every
   /// rank calls it with the same l, in the same order.
   void EnsureColumn(int l)
   {
      MFEM_VERIFY(0 <= l && l < P_.N, "QDImplicitFault::EnsureColumn: point " << l);
      if (slot_[l] < 0) { FetchColumn_(l); }
   }

   /// Slip rate, tractions and branch at the owned points of state U at
   /// time t (diagnostics and tests).  Collective.
   void Evaluate(real_t t, const real_t *U, std::vector<real_t> &V, std::vector<real_t> &tau,
                 std::vector<real_t> &sig, std::vector<char> &open)
   {
      Gather_(U);
      Tractions(t, tau, sig);
      V.resize(P_.n);
      open.resize(P_.n);
      for (int o = 0; o < P_.n; o++)
      {
         const PointEval e = EvalPoint(o, tau[o], sig[o], U[kPer * o + 1], U[kPer * o + 2]);
         V[o] = e.V;
         open[o] = e.open ? 1 : 0;
      }
   }

   /// The kernel's law stress at r = 0 and its derivative in sigma.
   real_t LawSigma(real_t sigma, real_t &ds) const
   {
      if (P_.sigma_floor > 0.0 && !P_.floor_taper)
      {
         ds = (sigma > P_.sigma_floor) ? 1.0 : 0.0;
         return std::max(sigma, P_.sigma_floor);
      }
      // the whole-law taper at r = 0, or no floor: max(sigma, 0) (the
      // kernel's solve takes sigma <= 0 as zero strength)
      ds = (sigma > 0.0) ? 1.0 : 0.0;
      return std::max(sigma, real_t(0.0));
   }

   /// Signed friction coefficient f(V, psi) = a asinh(V/(2 V0) e^(psi/a))
   /// (odd in V; the kernel's RSFrictionCoefficientPsi for V > 0) and its
   /// derivatives in V and psi.
   void Friction(real_t V, real_t psi, real_t a, real_t &f, real_t &fV, real_t &fpsi) const
   {
      const real_t av = std::abs(V), sg = (V > 0.0) ? 1.0 : ((V < 0.0) ? -1.0 : 0.0);
      const real_t pa = psi / a;
      if (pa > 700.0)
      {
         // the kernel's large-argument form: a ln(V/V0) + psi, floored at 0
         const real_t fu = (av > 0.0) ? std::max(a * std::log(av / P_.V0) + psi, real_t(0.0))
                                      : 0.0;
         f = sg * fu;
         fV = (av > 0.0 && fu > 0.0) ? a / av : a * std::exp(real_t(700.0)) / (2.0 * P_.V0);
         fpsi = (fu > 0.0) ? sg : 0.0;
         return;
      }
      const real_t q = std::exp(pa) / (2.0 * P_.V0);
      const real_t x = av * q;
      const real_t h = std::hypot(real_t(1.0), x);
      f = sg * RSFrictionCoefficientPsi(av, psi, a, P_.V0);
      fV = a * q / h;
      fpsi = sg * (x / h);
   }

   /// Per-point evaluation at (tau, sigma, psi), owned point o.
   struct PointEval
   {
      real_t V = 0.0, sig_law = 0.0, ds = 0.0, f = 0.0, fV = 0.0, fpsi = 0.0;
      real_t dpsi = 0.0;   ///< aging rate (closed) or 0 (open)
      bool open = false;
      int status = 0;
   };
   PointEval EvalPoint(int o, real_t tau, real_t sigma, real_t psi, real_t g) const
   {
      PointEval e;
      const int j = P_.lo + o;
      const real_t gam = (slot_[j] >= 0) ? cols_s_[slot_[j]][j] : 1.0;
      MFEM_VERIFY(g == 0.0 || slot_[j] >= 0, "QDImplicitFault: point " << j << " has g = "
                  << g << " m but no gap column");
      e.open = sigma < gam * g;
      e.sig_law = LawSigma(sigma, e.ds);
      const real_t Vabs = RSSolveSlipRatePsi(std::abs(tau), psi, e.sig_law, P_.eta[o], P_.a[o],
                                             P_.V0, e.status);
      e.V = (tau > 0.0) ? -Vabs : ((tau < 0.0) ? Vabs : 0.0);
      Friction(e.V, psi, P_.a[o], e.f, e.fV, e.fpsi);
      e.dpsi = e.open ? 0.0 : RSAgingRatePsi(std::abs(e.V), psi, P_.Dc[o], prm_);
      return e;
   }

   /// tau and sigma at the owned points from the gathered delta and g.
   void Tractions(real_t t, std::vector<real_t> &tau, std::vector<real_t> &sig) const
   {
      const int N = P_.N, n = P_.n;
      tau.assign(n, 0.0);
      sig.assign(n, 0.0);
      // the gap columns with g != 0
      std::vector<int> gl;
      for (int l = 0; l < N; l++)
      {
         if (g_all_[l] != 0.0)
         {
            MFEM_VERIFY(slot_[l] >= 0, "QDImplicitFault: point " << l << " has g = "
                        << g_all_[l] << " m but no gap column");
            gl.push_back(l);
         }
      }
      for (int o = 0; o < n; o++)
      {
         const int j = P_.lo + o;
         const real_t *gt = P_.Gt.data() + size_t(j) * N;
         const real_t *gs = P_.Gs.data() + size_t(j) * N;
         real_t st = 0.0, ss = 0.0;
         for (int l = 0; l < N; l++)
         {
            st += gt[l] * delta_all_[l];
            ss += gs[l] * delta_all_[l];
         }
         for (int l : gl)
         {
            st += cols_t_[slot_[l]][j] * g_all_[l];
            ss += cols_s_[slot_[l]][j] * g_all_[l];
         }
         real_t et = 0.0, es = 0.0;
         if (P_.extra_load) { P_.extra_load(t, j, et, es); }
         tau[o] = P_.tau_pre[o] + st + t * P_.TDt[j] + et;
         sig[o] = P_.sigma0[o] + ss + t * P_.TDs[j] + es;
      }
   }

   /// F(t, U, Udot) on the owned entries (3 n each).  Collective.
   void Residual(real_t t, const real_t *U, const real_t *Ud, real_t *F)
   {
      counts_.residuals++;
      Gather_(U);
      std::vector<real_t> tau, sig;
      Tractions(t, tau, sig);
      for (int o = 0; o < P_.n; o++)
      {
         const PointEval e = EvalPoint(o, tau[o], sig[o], U[kPer * o + 1], U[kPer * o + 2]);
         const real_t gam = GammaOr1_(P_.lo + o);
         F[kPer * o] = Ud[kPer * o] - e.V;
         F[kPer * o + 1] = Ud[kPer * o + 1] - e.dpsi;
         F[kPer * o + 2] = e.open ? sig[o] : gam * U[kPer * o + 2];
         res_branch_[o] = e.open ? 1 : 0;
      }
   }

   /// The Newton matrix at (t, U, Udot, shift): per-point numbers, the gap
   /// columns of new open points, and a new LU factor when the reuse rule
   /// asks for one (header).  shift >= 1e14 is PETSc's slope solve (3.20
   /// and later): the preconditioner is then the identity on the slip and
   /// psi rows and zero on the gap rows.  Collective.
   void Jacobian(real_t t, const real_t *U, const real_t *Ud, real_t shift)
   {
      (void) Ud;
      counts_.jacobians++;
      const int N = P_.N, n = P_.n;
      if (shift >= 1e14)
      {
         slope_mode_ = true;
         counts_.slope_solves++;
         return;
      }
      slope_mode_ = false;
      Gather_(U);
      std::vector<real_t> tau, sig;
      Tractions(t, tau, sig);
      // per point: D, Dt, fs, p, c, A, open
      constexpr int kRec = 7;
      std::vector<real_t> mine(size_t(kRec) * n), all(size_t(kRec) * N);
      for (int o = 0; o < n; o++)
      {
         const real_t psi = U[kPer * o + 1];
         const PointEval e = EvalPoint(o, tau[o], sig[o], psi, U[kPer * o + 2]);
         const real_t D = e.sig_law * e.fV + P_.eta[o];
         const real_t fs = e.f * e.ds;
         const real_t p = e.sig_law * e.fpsi / D;
         const real_t sgn = (e.V > 0.0) ? 1.0 : ((e.V < 0.0) ? -1.0 : 0.0);
         const real_t c = e.open ? 0.0 : (P_.b / P_.Dc[o]) * sgn;
         const real_t A = shift + (e.open ? 0.0
                                   : (P_.V0 / P_.Dc[o]) * std::exp((P_.f0 - psi) / P_.b));
         real_t *r = mine.data() + size_t(kRec) * o;
         r[0] = D;
         r[1] = D * shift * (1.0 - p * c / A);
         r[2] = fs;
         r[3] = p;
         r[4] = c;
         r[5] = A;
         r[6] = e.open ? 1.0 : 0.0;
      }
      GatherRec_(mine, all, kRec);
      ev_.shift = shift;
      ev_.D.resize(N); ev_.Dt.resize(N); ev_.fs.resize(N); ev_.p.resize(N);
      ev_.c.resize(N); ev_.A.resize(N); ev_.open.assign(N, 0);
      bool finite = true;
      for (int j = 0; j < N; j++)
      {
         const real_t *r = all.data() + size_t(kRec) * j;
         ev_.D[j] = r[0]; ev_.Dt[j] = r[1]; ev_.fs[j] = r[2]; ev_.p[j] = r[3];
         ev_.c[j] = r[4]; ev_.A[j] = r[5]; ev_.open[j] = (r[6] != 0.0) ? 1 : 0;
         for (int k = 0; k < 6; k++) { finite = finite && std::isfinite(r[k]); }
      }
      // the gap columns of points new on the open branch (same list on
      // every rank: from the gathered flags), in global order
      for (int j = 0; j < N; j++)
      {
         if (ev_.open[j] && slot_[j] < 0) { FetchColumn_(j); }
      }
      ev_.O.clear();
      for (int j = 0; j < N; j++) { if (ev_.open[j]) { ev_.O.push_back(j); } }
      counts_.max_open = std::max<long>(counts_.max_open, long(ev_.O.size()));
      ev_.gamma.assign(N, 1.0);
      for (int j = 0; j < N; j++) { if (slot_[j] >= 0) { ev_.gamma[j] = cols_s_[slot_[j]][j]; } }
      ev_.valid = true;
      const bool branch_moved = have_factor_ && (fa_.open != ev_.open);
      if (branch_moved) { counts_.branch_changes++; }
      bool refactor = !have_factor_ || fa_.shift != shift || fa_.open != ev_.open ||
                      S_.jacobian_reuse == "iteration" || slow_ ||
                      (S_.jacobian_reuse == "stage" && fa_.solve != solve_id_);
      if (!finite)
      {
         // a non-finite Newton matrix: no factor; the solve fails at the
         // residual (non-finite) or at the next convergence test
         bad_jacobian_ = true;
         return;
      }
      bad_jacobian_ = false;
      if (refactor) { Factor_(); }
   }

   /// y = J^-1 b with the current factor (slope mode: identity on the slip
   /// and psi rows, 0 on the gap rows).  Owned entries, 3 n.  Collective.
   void ApplyInverse(const real_t *b, real_t *y) const
   {
      const int N = P_.N, n = P_.n;
      if (slope_mode_)
      {
         for (int o = 0; o < n; o++)
         {
            y[kPer * o] = b[kPer * o];
            y[kPer * o + 1] = b[kPer * o + 1];
            y[kPer * o + 2] = 0.0;
         }
         return;
      }
      MFEM_VERIFY(have_factor_, "QDImplicitFault::ApplyInverse: no factor");
      std::vector<real_t> bl(b, b + size_t(kPer) * n), ball(size_t(kPer) * N);
      GatherRec_(bl, ball, kPer);
      const int m = static_cast<int>(fa_.O.size());
      // closed gap rows: y_g = b_g / gamma
      std::vector<real_t> yg(N, 0.0);
      std::vector<int> gc;   // closed points with y_g != 0
      for (int l = 0; l < N; l++)
      {
         if (fa_.open[l]) { continue; }
         const real_t bg = ball[size_t(kPer) * l + 2];
         if (bg != 0.0)
         {
            MFEM_VERIFY(slot_[l] >= 0, "QDImplicitFault: closed point " << l << " has a gap "
                        "residual " << bg << " but no gap column");
            yg[l] = bg / cols_s_[slot_[l]][l];
            gc.push_back(l);
         }
      }
      Vector z(N + m);
      for (int j = 0; j < N; j++)
      {
         const real_t bd = ball[size_t(kPer) * j], bp = ball[size_t(kPer) * j + 1];
         real_t r = fa_.D[j] * (bd - fa_.p[j] * (bp + fa_.c[j] * bd) / fa_.A[j]);
         for (int l : gc)
         {
            r -= (cols_t_[slot_[l]][j] + fa_.fs[j] * cols_s_[slot_[l]][j]) * yg[l];
         }
         z(j) = r;
      }
      for (int k = 0; k < m; k++)
      {
         const int j = fa_.O[k];
         real_t r = ball[size_t(kPer) * j + 2];
         for (int l : gc) { r -= cols_s_[slot_[l]][j] * yg[l]; }
         z(N + k) = r;
      }
      lu_.Solve(N + m, 1, z.GetData());
      for (int k = 0; k < m; k++) { yg[fa_.O[k]] = z(N + k); }
      for (int o = 0; o < n; o++)
      {
         const int j = P_.lo + o;
         const real_t bd = b[kPer * o], bp = b[kPer * o + 1];
         const real_t yd = z(j);
         y[kPer * o] = yd;
         y[kPer * o + 1] = (bp + fa_.c[j] * bd - fa_.c[j] * fa_.shift * yd) / fa_.A[j];
         y[kPer * o + 2] = yg[j];
      }
   }

   /// y = J x with the numbers of the last Jacobian evaluation (exact, no
   /// factor; tests).  Owned entries.  Collective.
   void ApplyJacobian(const real_t *x, real_t *y) const
   {
      MFEM_VERIFY(ev_.valid, "QDImplicitFault::ApplyJacobian: no Jacobian evaluation");
      const int N = P_.N, n = P_.n;
      std::vector<real_t> xl(x, x + size_t(kPer) * n), xall(size_t(kPer) * N);
      GatherRec_(xl, xall, kPer);
      for (int l = 0; l < N; l++)
      {
         MFEM_VERIFY(xall[size_t(kPer) * l + 2] == 0.0 || slot_[l] >= 0,
                     "QDImplicitFault::ApplyJacobian: x_g != 0 at point " << l
                     << " without a gap column");
      }
      for (int o = 0; o < n; o++)
      {
         const int j = P_.lo + o;
         const real_t *gt = P_.Gt.data() + size_t(j) * N;
         const real_t *gs = P_.Gs.data() + size_t(j) * N;
         real_t st = 0.0, ss = 0.0, gtt = 0.0, gss = 0.0;
         for (int l = 0; l < N; l++)
         {
            const real_t xd = xall[size_t(kPer) * l], xg = xall[size_t(kPer) * l + 2];
            st += gt[l] * xd;
            ss += gs[l] * xd;
            if (xg != 0.0)
            {
               gtt += cols_t_[slot_[l]][j] * xg;
               gss += cols_s_[slot_[l]][j] * xg;
            }
         }
         const real_t w = (st + ev_.fs[j] * ss + gtt + ev_.fs[j] * gss) / ev_.D[j];
         const real_t xd = x[kPer * o], xp = x[kPer * o + 1], xg = x[kPer * o + 2];
         y[kPer * o] = ev_.shift * xd + w + ev_.p[j] * xp;
         y[kPer * o + 1] = -ev_.c[j] * w + (ev_.A[j] - ev_.c[j] * ev_.p[j]) * xp;
         y[kPer * o + 2] = ev_.open[j] ? (ss + gss) : ev_.gamma[j] * xg;
      }
   }

   /// Consistent start (plan 2.10): g from the gap rows at fixed delta and
   /// psi, min(sigma, gamma g) = 0 at every point, by active set (open set
   /// O: sigma_O = 0 solved with Ggs_OO; elsewhere g = 0).  Returns the
   /// number of points with g > 0.  Collective.
   long ConsistentStart(real_t t, real_t *U, int max_passes = 50)
   {
      const int N = P_.N, n = P_.n;
      std::vector<char> open(N, 0), prev(N, 2);
      for (int pass = 0;; pass++)
      {
         MFEM_VERIFY(pass < max_passes, "QDImplicitFault::ConsistentStart: the open set did "
                     "not settle in " << max_passes << " passes");
         Gather_(U);
         std::vector<real_t> tau, sig;
         Tractions(t, tau, sig);
         std::vector<real_t> mine(n), sall(N);
         for (int o = 0; o < n; o++) { mine[o] = sig[o]; }
         GatherRec_(mine, sall, 1);
         for (int j = 0; j < N; j++)
         {
            const real_t gam = (slot_[j] >= 0) ? cols_s_[slot_[j]][j] : 1.0;
            // the active-set rule: in O, stays while g > 0 (sigma = 0 there
            // to round-off); out of O (g = 0), enters when sigma < 0
            open[j] = (sall[j] < gam * g_all_[j]) ? 1 : 0;
         }
         if (open == prev) { break; }
         for (int j = 0; j < N; j++) { if (open[j] && slot_[j] < 0) { FetchColumn_(j); } }
         std::vector<int> O;
         for (int j = 0; j < N; j++) { if (open[j]) { O.push_back(j); } }
         const int m = static_cast<int>(O.size());
         std::vector<real_t> gnew(N, 0.0);
         if (m > 0)
         {
            // sigma without the gap terms, at O: sall - sum_l Ggs(:, l) g_l
            DenseMatrix M(m);
            Vector r(m);
            for (int k = 0; k < m; k++)
            {
               const int j = O[k];
               real_t s = sall[j];
               for (int l = 0; l < N; l++)
               {
                  if (g_all_[l] != 0.0) { s -= cols_s_[slot_[l]][j] * g_all_[l]; }
               }
               r(k) = -s;
               for (int kk = 0; kk < m; kk++) { M(k, kk) = cols_s_[slot_[O[kk]]][j]; }
            }
            std::vector<int> piv(m);
            LUFactors lu(M.GetData(), piv.data());
            MFEM_VERIFY(lu.Factor(m), "QDImplicitFault::ConsistentStart: Ggs_OO is singular");
            lu.Solve(m, 1, r.GetData());
            for (int k = 0; k < m; k++) { gnew[O[k]] = r(k); }
         }
         for (int o = 0; o < n; o++) { U[kPer * o + 2] = gnew[P_.lo + o]; }
         prev = open;
      }
      long cnt = 0;
      for (int j = 0; j < N; j++) { if (g_all_[j] > 0.0) { cnt++; } }
      return cnt;
   }

   // ---- PETSc glue -------------------------------------------------------

   /// Set up @a ts (created by the caller, on comm) for option I: arkimex
   /// fully implicit of type esdirk_type, the IFunction and IJacobian of
   /// this object, SNES newtonls with a basic line search and this object's
   /// convergence test, KSP preonly with a PCSHELL (ApplyInverse), and
   /// -ts_max_snes_failures -1.  The TS adapt and tolerances are the
   /// caller's.
   PetscErrorCode SetupTS(TS ts)
   {
      PetscFunctionBeginUser;
      if (!J_)
      {
         PetscCall(MatCreateShell(comm_, kPer * P_.n, kPer * P_.n, PETSC_DETERMINE,
                                  PETSC_DETERMINE, this, &J_));
         PetscCall(MatSetUp(J_));
      }
      PetscCall(TSSetType(ts, TSARKIMEX));
      PetscCall(TSARKIMEXSetType(ts, S_.esdirk_type.c_str()));
      PetscCall(TSARKIMEXSetFullyImplicit(ts, PETSC_TRUE));
      PetscCall(TSSetIFunction(ts, nullptr, &QDImplicitFault::IFunction_, this));
      PetscCall(TSSetIJacobian(ts, J_, J_, &QDImplicitFault::IJacobian_, this));
      PetscCall(TSSetMaxSNESFailures(ts, -1));
      SNES snes;
      PetscCall(TSGetSNES(ts, &snes));
      PetscCall(SNESSetType(snes, SNESNEWTONLS));
      SNESLineSearch ls;
      PetscCall(SNESGetLineSearch(snes, &ls));
      PetscCall(SNESLineSearchSetType(ls, SNESLINESEARCHBASIC));
      PetscCall(SNESLineSearchSetPreCheck(ls, &QDImplicitFault::PreCheck_, this));
#if PETSC_VERSION_GE(3, 22, 0)
      PetscCall(SNESSetTolerances(snes, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT,
                                  S_.newton_max_it, PETSC_CURRENT));
#else
      PetscCall(SNESSetTolerances(snes, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT,
                                  S_.newton_max_it, PETSC_DEFAULT));
#endif
      PetscCall(SNESSetConvergenceTest(snes, &QDImplicitFault::Converged_, this, nullptr));
      KSP ksp;
      PetscCall(SNESGetKSP(snes, &ksp));
      PetscCall(KSPSetType(ksp, KSPPREONLY));
      PC pc;
      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PCSetType(pc, PCSHELL));
      PetscCall(PCShellSetContext(pc, this));
      PetscCall(PCShellSetApply(pc, &QDImplicitFault::PCApply_));
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// The TS error tolerances: atol and rtol on the slip and psi entries;
   /// the gap entries are left out of the local error test (vatol = huge,
   /// vrtol = 0), the way IDA's suppressalg treats algebraic unknowns.  Why:
   /// F_g has no g', so a stage solve fixes g from the gap row whatever the
   /// stage derivative is, and the derivative PETSc keeps for g at the
   /// explicit first stage is not g' (3.19: -F_g(t_n, U_n, 0) = 0; 3.20 and
   /// later: 0 at each restart); in the embedded error that slope error is
   /// O(h g') and would cap the step on an opening point.  g is set by delta
   /// through the gap rows at every stage, so its accuracy follows delta's.
   /// The norm is the caller's (-ts_adapt_wnormtype; the UCERF decks use
   /// infinity).  Collective.
   PetscErrorCode SetTolerances(TS ts)
   {
      PetscFunctionBeginUser;
      Vec va, vr;
      PetscCall(VecCreateMPI(comm_, kPer * P_.n, PETSC_DETERMINE, &va));
      PetscCall(VecDuplicate(va, &vr));
      PetscScalar *a, *r;
      PetscCall(VecGetArray(va, &a));
      PetscCall(VecGetArray(vr, &r));
      for (int o = 0; o < P_.n; o++)
      {
         a[kPer * o] = a[kPer * o + 1] = S_.atol;
         r[kPer * o] = r[kPer * o + 1] = S_.rtol;
         a[kPer * o + 2] = kSuppressed;
         r[kPer * o + 2] = 0.0;
      }
      PetscCall(VecRestoreArray(vr, &r));
      PetscCall(VecRestoreArray(va, &a));
      PetscCall(TSSetTolerances(ts, S_.atol, va, S_.rtol, vr));
      PetscCall(VecDestroy(&va));
      PetscCall(VecDestroy(&vr));
      PetscFunctionReturn(PETSC_SUCCESS);
   }
   static constexpr real_t kSuppressed = 1e250;   ///< vatol of the gap entries

   /// Drop the factor (the next Jacobian refactors); e.g. after the state
   /// was changed outside the TS.
   void DropFactor() { have_factor_ = false; slow_ = false; }

private:
   struct Snapshot
   {
      bool valid = false;
      real_t shift = 0.0;
      long solve = -1;
      std::vector<real_t> D, Dt, fs, p, c, A, gamma;
      std::vector<char> open;
      std::vector<int> O;
   };

   real_t GammaOr1_(int l) const { return (slot_[l] >= 0) ? cols_s_[slot_[l]][l] : 1.0; }

   /// delta_all_ and g_all_ from the owned U.
   void Gather_(const real_t *U)
   {
      const int n = P_.n;
      std::vector<real_t> d(n), g(n);
      for (int o = 0; o < n; o++) { d[o] = U[kPer * o]; g[o] = U[kPer * o + 2]; }
      MPI_Allgatherv(d.data(), n, MPITypeMap<real_t>::mpi_type, delta_all_.data(), cnt_.data(),
                     dsp_.data(), MPITypeMap<real_t>::mpi_type, comm_);
      MPI_Allgatherv(g.data(), n, MPITypeMap<real_t>::mpi_type, g_all_.data(), cnt_.data(),
                     dsp_.data(), MPITypeMap<real_t>::mpi_type, comm_);
   }

   /// k values per point, owned -> global.
   void GatherRec_(const std::vector<real_t> &mine, std::vector<real_t> &all, int k) const
   {
      std::vector<int> c(nranks_), d(nranks_);
      for (int r = 0; r < nranks_; r++) { c[r] = k * cnt_[r]; d[r] = k * dsp_[r]; }
      all.resize(size_t(k) * P_.N);
      MPI_Allgatherv(mine.data(), k * P_.n, MPITypeMap<real_t>::mpi_type, all.data(), c.data(),
                     d.data(), MPITypeMap<real_t>::mpi_type, comm_);
   }

   void FetchColumn_(int l)
   {
      MFEM_VERIFY(P_.gap_column, "QDImplicitFault: point " << l << " is on the open branch, "
                  "but there is no gap_column callback");
      std::vector<real_t> ct, cs;
      P_.gap_column(l, ct, cs);
      MFEM_VERIFY(ct.size() == size_t(P_.N) && cs.size() == size_t(P_.N),
                  "QDImplicitFault: gap column of point " << l << " has " << ct.size()
                  << " / " << cs.size() << " entries for " << P_.N << " points");
      MFEM_VERIFY(cs[l] > 0.0, "QDImplicitFault: gap column of point " << l << ": its own "
                  "normal stress per unit opening is " << cs[l] << " Pa/m (must be > 0)");
      slot_[l] = static_cast<int>(cols_t_.size());
      cols_t_.push_back(std::move(ct));
      cols_s_.push_back(std::move(cs));
      counts_.columns++;
   }

   /// Assemble the reduced matrix from ev_ and factor it.
   void Factor_()
   {
      const int N = P_.N;
      const int m = static_cast<int>(ev_.O.size());
      const int M = N + m;
      lu_mat_.SetSize(M);
      lu_piv_.resize(M);
      real_t *a = lu_mat_.GetData();   // column-major
      for (int j = 0; j < N; j++)
      {
         const real_t *gt = P_.Gt.data() + size_t(j) * N;
         const real_t *gs = P_.Gs.data() + size_t(j) * N;
         const real_t fs = ev_.fs[j];
         for (int l = 0; l < N; l++) { a[size_t(l) * M + j] = gt[l] + fs * gs[l]; }
         a[size_t(j) * M + j] += ev_.Dt[j];
         for (int k = 0; k < m; k++)
         {
            const int l = ev_.O[k];
            a[size_t(N + k) * M + j] = cols_t_[slot_[l]][j] + fs * cols_s_[slot_[l]][j];
         }
      }
      for (int k = 0; k < m; k++)
      {
         const int j = ev_.O[k];
         const real_t *gs = P_.Gs.data() + size_t(j) * N;
         for (int l = 0; l < N; l++) { a[size_t(l) * M + N + k] = gs[l]; }
         for (int kk = 0; kk < m; kk++)
         {
            a[size_t(N + kk) * M + N + k] = cols_s_[slot_[ev_.O[kk]]][j];
         }
      }
      lu_ = LUFactors(a, lu_piv_.data());
      const bool ok = lu_.Factor(M);
      fa_ = ev_;
      fa_.solve = solve_id_;
      have_factor_ = ok;
      slow_ = false;
      counts_.factors++;
      factor_singular_ = !ok;
   }

   static PetscErrorCode IFunction_(TS, PetscReal t, Vec U, Vec Ud, Vec F, void *ctx)
   {
      PetscFunctionBeginUser;
      QDImplicitFault *self = static_cast<QDImplicitFault *>(ctx);
      const PetscScalar *u, *ud;
      PetscScalar *f;
      PetscCall(VecGetArrayRead(U, &u));
      PetscCall(VecGetArrayRead(Ud, &ud));
      PetscCall(VecGetArray(F, &f));
      self->Residual(static_cast<real_t>(t), reinterpret_cast<const real_t *>(u),
                     reinterpret_cast<const real_t *>(ud), reinterpret_cast<real_t *>(f));
      PetscCall(VecRestoreArray(F, &f));
      PetscCall(VecRestoreArrayRead(Ud, &ud));
      PetscCall(VecRestoreArrayRead(U, &u));
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   static PetscErrorCode IJacobian_(TS, PetscReal t, Vec U, Vec Ud, PetscReal shift, Mat A,
                                    Mat B, void *ctx)
   {
      PetscFunctionBeginUser;
      QDImplicitFault *self = static_cast<QDImplicitFault *>(ctx);
      const PetscScalar *u, *ud;
      PetscCall(VecGetArrayRead(U, &u));
      PetscCall(VecGetArrayRead(Ud, &ud));
      self->Jacobian(static_cast<real_t>(t), reinterpret_cast<const real_t *>(u),
                     reinterpret_cast<const real_t *>(ud), static_cast<real_t>(shift));
      PetscCall(VecRestoreArrayRead(Ud, &ud));
      PetscCall(VecRestoreArrayRead(U, &u));
      PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
      PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
      if (B != A)
      {
         PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
         PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
      }
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   static PetscErrorCode PCApply_(PC pc, Vec b, Vec y)
   {
      PetscFunctionBeginUser;
      void *ctx = nullptr;
      PetscCall(PCShellGetContext(pc, &ctx));
      const QDImplicitFault *self = static_cast<const QDImplicitFault *>(ctx);
      if (!self->slope_mode_ && (self->bad_jacobian_ || !self->have_factor_ ||
                                 self->factor_singular_))
      {
         // no usable factor: a NaN update fails the solve, so the TS
         // rejects the step
         PetscCall(VecSet(y, std::numeric_limits<PetscReal>::quiet_NaN()));
         PetscFunctionReturn(PETSC_SUCCESS);
      }
      const PetscScalar *bb;
      PetscScalar *yy;
      PetscCall(VecGetArrayRead(b, &bb));
      PetscCall(VecGetArray(y, &yy));
      self->ApplyInverse(reinterpret_cast<const real_t *>(bb), reinterpret_cast<real_t *>(yy));
      PetscCall(VecRestoreArray(y, &yy));
      PetscCall(VecRestoreArrayRead(b, &bb));
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// The psi change of one point (header): d_N Newton's, s' = shift - c
   /// p, kappa the aging slope, sc the log-damping scale, cap the cap.
   static real_t PsiChange(real_t dN, real_t sp, real_t kappa, real_t b, real_t sc, real_t cap)
   {
      if (!(std::abs(dN) > 1e-3 * b)) { return dN; }
      if (kappa <= 0.0) { return std::max(-cap, std::min(dN, cap)); }   // open: linear row
      if (!(sp > 0.0))
      {
         const real_t a = std::abs(dN);
         const real_t d = (a > sc) ? std::min(sc * (1.0 + std::log(a / sc)), cap) : a;
         return (dN > 0.0) ? d : -d;
      }
      auto phi = [&](real_t d) { return sp * d + kappa * b * (1.0 - std::exp(-d / b))
                                        - (sp + kappa) * dN; };
      real_t lo, hi;
      if (dN < 0.0)
      {
         lo = std::max(dN, -cap);
         hi = 0.0;
         if (phi(lo) >= 0.0) { return lo; }   // the root is below the cap (or at it)
      }
      else
      {
         lo = dN;
         hi = std::min(dN * (sp + kappa) / sp, cap);
         if (!(lo < hi)) { return std::min(dN, cap); }
         if (phi(hi) <= 0.0) { return hi; }
      }
      int st = 0;
      const real_t d = RSZeroIn(lo, hi, phi, st);
      return (st == 0 && std::isfinite(d)) ? d : std::max(-cap, std::min(dN, cap));
   }

   /// The psi change of the Newton update, point by point (header): X_new
   /// = X - Y, so the change is -Y_psi.  Collective.
   static PetscErrorCode PreCheck_(SNESLineSearch, Vec X, Vec Y, PetscBool *changed, void *ctx)
   {
      PetscFunctionBeginUser;
      (void) X;
      QDImplicitFault *self = static_cast<QDImplicitFault *>(ctx);
      *changed = PETSC_FALSE;
      self->clipped_ = false;
      if (self->slope_mode_) { PetscFunctionReturn(PETSC_SUCCESS); }
      if (self->bad_jacobian_ || !self->ev_.valid) { PetscFunctionReturn(PETSC_SUCCESS); }
      const real_t b = self->P_.b;
      const real_t sc = self->S_.newton_psi_scale * b;
      const real_t cap = self->S_.newton_psi_cap * b;
      const Snapshot &ev = self->ev_;   // the Jacobian numbers at this iterate
      PetscScalar *y;
      PetscCall(VecGetArray(Y, &y));
      int clipped = 0;
      for (int o = 0; o < self->P_.n; o++)
      {
         const int j = self->P_.lo + o;
         PetscScalar &yp = y[kPer * o + 1];
         const real_t dN = -static_cast<real_t>(PetscRealPart(yp));
         const real_t kappa = ev.A[j] - ev.shift;
         const real_t sp = ev.shift - ev.c[j] * ev.p[j];
         const real_t d = PsiChange(dN, sp, kappa, b, sc, cap);
         if (d != dN)
         {
            yp = -d;
            if (std::abs(d - dN) > 0.1 * std::abs(dN)) { clipped = 1; }
         }
      }
      PetscCall(VecRestoreArray(Y, &y));
      PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &clipped, 1, MPI_INT, MPI_MAX, self->comm_));
      if (clipped)
      {
         *changed = PETSC_TRUE;
         self->clipped_ = true;
         self->counts_.clipped_its++;
      }
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// The Newton test (header).  Collective.
   static PetscErrorCode Converged_(SNES snes, PetscInt it, PetscReal xnorm, PetscReal gnorm,
                                    PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
   {
      PetscFunctionBeginUser;
      QDImplicitFault *self = static_cast<QDImplicitFault *>(ctx);
      *reason = SNES_CONVERGED_ITERATING;
      if (it == 0)
      {
         self->solve_id_++;
         self->wprev_ = -1.0;
         self->clipped_ = self->prev_clipped_ = false;
         self->counts_.newton_solves++;
         if (!std::isfinite(static_cast<double>(fnorm)))
         {
            *reason = SNES_DIVERGED_FNORM_NAN;
            self->counts_.newton_fails++;
         }
         else if (fnorm == 0.0) { *reason = SNES_CONVERGED_FNORM_ABS; }
         PetscFunctionReturn(PETSC_SUCCESS);
      }
      self->counts_.newton_its++;
      if (!std::isfinite(static_cast<double>(fnorm)))
      {
         *reason = SNES_DIVERGED_FNORM_NAN;
         self->counts_.newton_fails++;
         PetscFunctionReturn(PETSC_SUCCESS);
      }
      if (self->slope_mode_)
      {
         PetscCall(SNESConvergedDefault(snes, it, xnorm, gnorm, fnorm, reason, nullptr));
         PetscFunctionReturn(PETSC_SUCCESS);
      }
      // the update in the TS error weights, max over all entries and ranks
      Vec dx, x;
      PetscCall(SNESGetSolutionUpdate(snes, &dx));
      PetscCall(SNESGetSolution(snes, &x));
      const PetscScalar *d, *xx;
      PetscInt nl = 0;
      PetscCall(VecGetLocalSize(x, &nl));
      PetscCall(VecGetArrayRead(dx, &d));
      PetscCall(VecGetArrayRead(x, &xx));
      double w = 0.0;
      for (PetscInt i = 0; i < nl; i++)
      {
         const double di = std::abs(static_cast<double>(PetscRealPart(d[i])));
         const double xi = std::abs(static_cast<double>(PetscRealPart(xx[i])));
         const double wi = di / (self->S_.atol + self->S_.rtol * xi);
         w = std::isfinite(wi) ? std::max(w, wi) : std::numeric_limits<double>::infinity();
      }
      PetscCall(VecRestoreArrayRead(x, &xx));
      PetscCall(VecRestoreArrayRead(dx, &d));
      // the iterate's branch set (last residual) against the factor's
      int moved = 0;
      if (self->have_factor_)
      {
         for (int o = 0; o < self->P_.n; o++)
         {
            if (self->res_branch_[o] != self->fa_.open[self->P_.lo + o]) { moved = 1; break; }
         }
      }
      double in[2] = {w, static_cast<double>(moved)}, out[2];
      MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(snes));
      PetscCallMPI(MPI_Allreduce(in, out, 2, MPI_DOUBLE, MPI_MAX, comm));
      w = out[0];
      moved = (out[1] != 0.0) ? 1 : 0;
      self->last_w_ = w;
      const bool clipped = self->clipped_;
      const bool rate_ok = !clipped && !self->prev_clipped_ && self->wprev_ > 0.0;
      const double theta = rate_ok ? w / self->wprev_ : -1.0;
      const double tol = self->S_.newton_tol;
      if (!std::isfinite(w))
      {
         *reason = SNES_DIVERGED_FNORM_NAN;
         self->counts_.newton_fails++;
      }
      else if (!moved && !clipped &&
               (w <= tol || (theta > 0.0 && theta < 1.0 && theta / (1.0 - theta) * w <= tol)))
      {
         *reason = SNES_CONVERGED_SNORM_RELATIVE;
      }
      else if (it >= self->S_.newton_max_it)
      {
         *reason = SNES_DIVERGED_MAX_IT;
         self->counts_.newton_fails++;
         self->slow_ = true;   // refactor at the next try
      }
      else if (clipped || theta > self->S_.newton_theta_refactor)
      {
         // a clipped update: the linearization is far off (a factor kept
         // from a distant iterate can make Newton cycle between +cap and
         // -cap, T1a at h = 1 s), so take a fresh one
         self->slow_ = true;   // a fresh factor at the next iteration
      }
      self->wprev_ = w;
      self->prev_clipped_ = clipped;
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   MPI_Comm comm_;
   int rank_ = 0, nranks_ = 1;
   Problem P_;
   QDImplicitSettings S_;
   RSDeviceParams prm_{};
   std::vector<int> cnt_, dsp_;
   std::vector<real_t> delta_all_, g_all_;
   std::vector<int> slot_;                       ///< global point -> column slot
   std::vector<std::vector<real_t>> cols_t_, cols_s_;
   std::vector<char> res_branch_;                ///< owned: open at the last residual
   Snapshot ev_, fa_;                            ///< last evaluation, factor's
   DenseMatrix lu_mat_;
   std::vector<int> lu_piv_;
   LUFactors lu_;
   bool have_factor_ = false, factor_singular_ = false, bad_jacobian_ = false;
   bool slope_mode_ = false, slow_ = false, clipped_ = false, prev_clipped_ = false;
   long solve_id_ = 0;
   double wprev_ = -1.0, last_w_ = 0.0;
   Counts counts_;
   Mat J_ = nullptr;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_USE_PETSC && MFEM_USE_MPI
#endif // MFEM_SEAS_CYCLE_QD_IMPLICIT_HPP
