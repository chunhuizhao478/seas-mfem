// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// qd_implicit.hpp — option I of the QD step at sigma_n' -> 0 (Chunhui Zhao,
// 2026-11-09: "let's proceed with implementation, take all decisions p1-p6
// by recommended items"; plan /mnt/project-files/fully-dynamic-fault-opening/
// qd-implicit-friction/, Claude Doc "Option I: implementation and test
// plan", rev 14; note qd-implicit-friction-rev2.pdf).  Phase I1: the fault
// system, its Newton matrix and the PETSc set-up, on any dense G.  Phase I2:
// PlanarProblem (end of file) builds it from the 2-D planar dense Green's
// function (QDGreenOperator::GatherPlanar) and the solve path's unit
// opening responses (ElasticityDomainOperator::QDUnitColumn); the QD phase
// calls it in I3.
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
// (P3, 2.4) with no branch in the slip row.  sigma_law = Lambda(s, s) of
// the kernel at r = 0 (RSLawSigmaN): max(s, 0) with the whole-law taper or
// without a floor, max(s, floor) with the hard floor, at
//
//    s_j = sigma_j - gamma_j g_j       the normal stress without the point's
//                                      own gap term (the trial normal stress
//                                      of semismooth contact Newton)
//
// open_j: s_j < 0, i.e. sigma_j < gamma_j g_j, the branch of the min (2.5,
// 2.9); psi is frozen while open (D-B3).  gamma_j = Ggs_jj, the point's own
// normal stress per unit opening ([cycle.opening] contact_scale = "auto");
// it only enters where g_j != 0, which needs the cached column anyway.
//
// Why s and not sigma in the law (phase I2, test G6): at a solution they
// give the same strength (closed: g = 0, s = sigma; open: sigma = 0 and
// s < 0, both zero strength, or the hard floor), so the solutions are the
// same.  They differ at Newton iterates of an open point, where sigma is
// not yet 0: with max(sigma, 0) an iterate with sigma > 0 switches the
// whole friction strength sigma f on at tau ~ 0, V* drops from |tau| / eta
// to nothing within sigma ~ |tau| / f (0.01 Pa), and the linearization at
// sigma = 0 answers a 0.1 MPa change of sigma with metres of slip (the
// opening run on the real G: 30 factors and more than one failed solve per
// step while a gap closed).  With s the open point's slip row is tau +
// eta V = 0 at every iterate, linear, and the residual stays continuous at
// the branch (s = 0 on both sides).
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
//    fs_j = f_j dsigma_law/ds                (strength change with s_j)
//    p_j  = sigma_law f_psi / D_j            (dF_delta/dpsi_j)
//    c_j  = (b / Dc_j) sgn(V_j), 0 when open (psi row's pull from |V|)
//    A_j  = shift + (V0 / Dc_j) exp((f0 - psi_j)/b), shift when open
//
//    J x |_delta,j = shift x_delta,j + w_j + p_j x_psi,j
//    J x |_psi,j   = -c_j w_j + (A_j - c_j p_j) x_psi,j
//    J x |_g,j     = closed: gamma_j x_g,j;  open: [Gs x_delta + Ggs x_g]_j
//    w_j = ( [(Gt + fs Gs) x_delta]_j + [(Ggt + fs Ggs~) x_g]_j ) / D_j
//    Ggs~ = Ggs without its diagonal (ds_j/dg_j = 0)
//
// Eliminating psi per point and the closed gap rows (x_g = b_g / gamma)
// leaves the plan's (2.8) on (delta, g_O), times D_j on the slip rows:
//
//    [ diag(Dt) + Gt + fs Gs     Ggt_:O + fs Ggs~_:O] [y_delta]   [r_delta]
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
// [cycle.qd] dense_lu = "shared" (opt-in; "serial" is the default and the
// path above, bit for bit): the same reduced matrix is factored once by all
// ranks together (cycle/par_dense_lu.hpp, ScaLAPACK), each rank assembling
// only the blocks it stores, and every solve goes through that factor.
// Why: with "serial" every rank spends the same 2/3 (N + m)^3 flops on one
// core (UCERF, N = 3,290, 144 ranks: about 5 factors of 0.8 s in a 4.8 s
// step).  The solution agrees with the serial one to round-off, not bit for
// bit, and depends on the rank count in the last digits.
//
// Newton on each point's own stage equation ([cycle.qd] newton_local, on by
// default; phase I2, test G6).  Why: V* is exponential in tau / (a
// sigma_law) (f ~ a ln V + psi), so Newton on F, which linearizes V* in the
// slips, gains one unit of that exponent per iteration from above and
// throws metres of slip from below.  Where a sliding point's sigma falls
// through its last MPa (tau / (a sigma) ~ 65, a 7% change of sigma in a
// step), that took 60 to 80 factors and a rejected step per step.  The cure
// is the one of every implicit rate-and-state solver: solve the point's
// friction equation for V with the point's own elastic response inside it.
// In a stage, delta_j = (delta'_j - const) / shift, so with the other
// points held, the point's row is the scalar equation
//
//    tau_j - keff_j Vit_j + (eta_j + keff_j) V + sigma_law,j f(V, psi_j) = 0,
//    keff_j = max(Gt_jj, 0) / shift,   Vit_j = delta'_j at the iterate
//
// (tau_j - keff_j Vit_j does not depend on delta_j).  Its root V_loc,j is
// the kernel's bracketed solve with eta + keff for eta.  Newton then works
// on
//
//    G_delta,j = delta'_j - V_loc,j        G_psi,j = psi'_j - aging(|V_loc,j|, psi_j)
//    G_g,j     = F_g,j
//
// G = 0 exactly where F = 0 (V_loc = Vit = V* there), so the TS sees the
// same stage solutions, and F itself (the IFunction, the first-stage slope)
// is untouched.  For one point G is solved in one iteration at any step
// size; with many, Newton only has the coupling between points left, which
// is linear in the slips.  Its matrix is the one below with the law
// evaluated at V_loc, behind the row map q_j = D_j / (D_j + keff_j):
//
//    (J_G x)_delta = q (J x)_delta     (J_G x)_psi = (J x)_psi + c (1 - q)(J x)_delta
//
// so the same reduced matrix and factor serve, and the update is
// ApplyInverse of (G_delta / q, G_psi - c (keff / D) G_delta, G_g).  PETSc
// hands the preconditioner F; the slip and psi rows of G are kept from the
// Jacobian evaluation at the same iterate (checked: the same U as the last
// residual, and no residual since; otherwise the update is the plain one
// and Counts::plain_applies counts it).  newton_local = false is Newton on
// F as in phase I1.
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
// (p -> p D / (D + keff) with newton_local: V_loc's slope in psi),
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
#include "../domain/qd_green_operator.hpp"
#include "par_dense_lu.hpp"

#if defined(MFEM_USE_PETSC) && defined(MFEM_USE_MPI)
#include "petsc.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
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
   bool newton_local = true;          ///< Newton on each point's own stage equation
                                      ///< (file comment); false: on F as it stands
   int newton_max_it = 20;            ///< more: a failed solve, the step is rejected
   std::string jacobian_reuse = "step";   ///< "step", "stage" or "iteration"
   std::string dense_lu = "serial";   ///< "serial" (every rank factors) or "shared"
                                      ///< (one ScaLAPACK factor, cycle/par_dense_lu.hpp)
   std::string stage_guess = "previous";   ///< where a stage's Newton starts: "previous"
                                      ///< (the stage before) or "extrapolate" (the last
                                      ///< step's dense output, SetupTS)
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
      long plain_applies = 0;    ///< newton_local: updates that fell back to F (not expected)
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
      MFEM_VERIFY(S_.dense_lu == "serial" || S_.dense_lu == "shared",
                  "QDImplicitFault: dense_lu \"" << S_.dense_lu << "\" (serial or shared)");
      MFEM_VERIFY(S_.stage_guess == "previous" || S_.stage_guess == "extrapolate",
                  "QDImplicitFault: stage_guess \"" << S_.stage_guess
                  << "\" (previous or extrapolate)");
      MFEM_VERIFY(S_.dense_lu == "serial" || ParDenseLU::Available(),
                  "QDImplicitFault: dense_lu = \"shared\" needs ScaLAPACK (a build with "
                  "MFEM_USE_MUMPS, double precision)");
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
   /// With keff > 0 the slip rate is the point's own stage solve (file
   /// comment, "Newton on each point's own stage equation"): tau - keff
   /// Vit + (eta + keff) V + sigma_law f(V, psi) = 0; keff = 0 is the law.
   PointEval EvalPoint(int o, real_t tau, real_t sigma, real_t psi, real_t g,
                       real_t keff = 0.0, real_t Vit = 0.0) const
   {
      PointEval e;
      const int j = P_.lo + o;
      const real_t gam = (slot_[j] >= 0) ? cols_s_[slot_[j]][j] : 1.0;
      // a non-finite g (a Newton iterate gone bad) evaluates to non-finite
      // values and fails the solve; only a finite g needs its column
      MFEM_VERIFY(g == 0.0 || !std::isfinite(g) || slot_[j] >= 0, "QDImplicitFault: point "
                  << j << " has g = " << g << " m but no gap column");
      // the normal stress without the point's own gap term (file comment)
      const real_t st = sigma - gam * g;
      e.open = st < 0.0;
      e.sig_law = LawSigma(st, e.ds);
      const real_t tau0 = (keff > 0.0) ? tau - keff * Vit : tau;
      const real_t eta0 = (keff > 0.0) ? P_.eta[o] + keff : P_.eta[o];
      const real_t Vabs = RSSolveSlipRatePsi(std::abs(tau0), psi, e.sig_law, eta0, P_.a[o],
                                             P_.V0, e.status);
      e.V = (tau0 > 0.0) ? -Vabs : ((tau0 < 0.0) ? Vabs : 0.0);
      Friction(e.V, psi, P_.a[o], e.f, e.fV, e.fpsi);
      e.dpsi = e.open ? 0.0 : RSAgingRatePsi(std::abs(e.V), psi, P_.Dc[o], prm_);
      if (!std::isfinite(tau) || !std::isfinite(sigma) || !std::isfinite(psi))
      {
         // a non-finite iterate: the sign tests above would give V = 0 and
         // a finite residual; keep it non-finite so the solve fails
         e.V = e.dpsi = e.f = std::numeric_limits<real_t>::quiet_NaN();
      }
      return e;
   }

   /// tau and sigma at the owned points from the gathered delta and g.
   void Tractions(real_t t, std::vector<real_t> &tau, std::vector<real_t> &sig) const
   {
      const int N = P_.N, n = P_.n;
      tau.assign(n, 0.0);
      sig.assign(n, 0.0);
      // A non-finite slip or gap anywhere (a Newton iterate gone bad, e.g.
      // after a too-long trial step): non-finite tractions everywhere, so
      // the residual is non-finite, the Newton test reports
      // SNES_DIVERGED_FNORM_NAN and TS rejects the step and shrinks dt.
      // Before this, a non-finite g at a point without a gap column stopped
      // the run in the check below (UCERF T4b, Vista 1133642).  The same
      // result on every rank: delta_all_ and g_all_ are gathered.
      for (int l = 0; l < N; l++)
      {
         if (!std::isfinite(delta_all_[l]) || !std::isfinite(g_all_[l]))
         {
            tau.assign(n, std::numeric_limits<real_t>::quiet_NaN());
            sig.assign(n, std::numeric_limits<real_t>::quiet_NaN());
            return;
         }
      }
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
      res_U_.assign(U, U + size_t(kPer) * P_.n);
      res_serial_++;
   }

   /// newton_local: the residual Newton works on, G(t, U, Udot; shift).  The
   /// slip row is delta' - V_loc and the psi row psi' - aging(|V_loc|, psi),
   /// V_loc the point's own stage solve (EvalPoint with keff = max(Gt_jj,
   /// 0) / shift and Vit = delta'); the gap row is F's.  G = 0 exactly where
   /// F = 0.  Owned entries (3 n).  Collective.
   void LocalResidual(real_t t, const real_t *U, const real_t *Ud, real_t shift, real_t *G)
   {
      Gather_(U);
      std::vector<real_t> tau, sig;
      Tractions(t, tau, sig);
      for (int o = 0; o < P_.n; o++)
      {
         const PointEval e = EvalPoint(o, tau[o], sig[o], U[kPer * o + 1], U[kPer * o + 2],
                                       Keff_(o, shift), Ud[kPer * o]);
         G[kPer * o] = Ud[kPer * o] - e.V;
         G[kPer * o + 1] = Ud[kPer * o + 1] - e.dpsi;
         G[kPer * o + 2] = e.open ? sig[o] : GammaOr1_(P_.lo + o) * U[kPer * o + 2];
      }
   }

   /// The Newton matrix at (t, U, Udot, shift): per-point numbers, the gap
   /// columns of new open points, and a new LU factor when the reuse rule
   /// asks for one (header).  shift >= 1e14 is PETSc's slope solve (3.20
   /// and later): the preconditioner is then the identity on the slip and
   /// psi rows and zero on the gap rows.  Collective.
   void Jacobian(real_t t, const real_t *U, const real_t *Ud, real_t shift)
   {
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
      // per point: D, Dt, fs, p, c, A, open, keff.  newton_local: the law
      // at V_loc, the point's own stage solve; else at V* (the law's)
      constexpr int kRec = 8;
      std::vector<real_t> mine(size_t(kRec) * n), all(size_t(kRec) * N);
      loc_d_.resize(n);
      loc_p_.resize(n);
      for (int o = 0; o < n; o++)
      {
         const real_t psi = U[kPer * o + 1];
         const real_t keff = S_.newton_local ? Keff_(o, shift) : 0.0;
         const PointEval e = EvalPoint(o, tau[o], sig[o], psi, U[kPer * o + 2], keff,
                                       Ud[kPer * o]);
         loc_d_[o] = Ud[kPer * o] - e.V;
         loc_p_[o] = Ud[kPer * o + 1] - e.dpsi;
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
         r[7] = keff;
      }
      // the update is built from loc_d_, loc_p_ only if they belong to the
      // residual PETSc hands the preconditioner: the last one, at this U
      loc_ok_ = S_.newton_local && res_U_.size() == size_t(kPer) * n &&
                std::equal(res_U_.begin(), res_U_.end(), U);
      loc_serial_ = res_serial_;
      GatherRec_(mine, all, kRec);
      ev_.shift = shift;
      ev_.local = S_.newton_local;
      ev_.D.resize(N); ev_.Dt.resize(N); ev_.fs.resize(N); ev_.p.resize(N);
      ev_.c.resize(N); ev_.A.resize(N); ev_.keff.resize(N); ev_.open.assign(N, 0);
      bool finite = true;
      for (int j = 0; j < N; j++)
      {
         const real_t *r = all.data() + size_t(kRec) * j;
         ev_.D[j] = r[0]; ev_.Dt[j] = r[1]; ev_.fs[j] = r[2]; ev_.p[j] = r[3];
         ev_.c[j] = r[4]; ev_.A[j] = r[5]; ev_.open[j] = (r[6] != 0.0) ? 1 : 0;
         ev_.keff[j] = r[7];
         for (int k = 0; k < 6; k++) { finite = finite && std::isfinite(r[k]); }
         finite = finite && std::isfinite(r[7]);
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
            r -= (cols_t_[slot_[l]][j] +
                  ((l == j) ? 0.0 : fa_.fs[j] * cols_s_[slot_[l]][j])) * yg[l];
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
      if (plu_)
      {
         Vector zz(N + m);
         plu_->Solve(z.GetData(), zz.GetData());
         z = zz;
      }
      else { lu_.Solve(N + m, 1, z.GetData()); }
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
   /// newton_local: y = (dG/dU + shift dG/dUdot) x at the last Jacobian
   /// evaluation, G of LocalResidual.  With J x of ApplyJacobian (the same
   /// formulas, the law at V_loc) and q_j = D_j / (D_j + keff_j):
   ///    slip: q (J x)_delta      psi: (J x)_psi + c (1 - q)(J x)_delta
   /// and the gap rows as they are.
   void ApplyLocalJacobian(const real_t *x, real_t *y) const
   {
      ApplyJacobian(x, y);
      for (int o = 0; o < P_.n; o++)
      {
         const int j = P_.lo + o;
         const real_t q = ev_.D[j] / (ev_.D[j] + ev_.keff[j]), jd = y[kPer * o];
         y[kPer * o] = q * jd;
         y[kPer * o + 1] += ev_.c[j] * (1.0 - q) * jd;
      }
   }

   /// newton_local: y = (dG/dU + shift dG/dUdot)^-1 b with the factor's
   /// numbers: ApplyInverse of (b_delta / q, b_psi - c (keff / D) b_delta,
   /// b_g), the inverse of the row map of ApplyLocalJacobian.
   void ApplyLocalInverse(const real_t *b, real_t *y) const
   {
      std::vector<real_t> bh(b, b + size_t(kPer) * P_.n);
      if (!slope_mode_)
      {
         MFEM_VERIFY(have_factor_, "QDImplicitFault::ApplyLocalInverse: no factor");
         for (int o = 0; o < P_.n; o++)
         {
            const int j = P_.lo + o;
            const real_t r = fa_.keff[j] / fa_.D[j], bd = b[kPer * o];
            bh[kPer * o] = (1.0 + r) * bd;
            bh[kPer * o + 1] = b[kPer * o + 1] - fa_.c[j] * r * bd;
         }
      }
      ApplyInverse(bh.data(), y);
   }

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
         const real_t xd = x[kPer * o], xp = x[kPer * o + 1], xg = x[kPer * o + 2];
         // ds_j = d sigma_j - gamma_j dg_j: Ggs without its diagonal
         const real_t gss_t = (xg != 0.0) ? gss - cols_s_[slot_[j]][j] * xg : gss;
         const real_t w = (st + ev_.fs[j] * ss + gtt + ev_.fs[j] * gss_t) / ev_.D[j];
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
   /// caller's.  stage_guess = "extrapolate" switches on PETSc's
   /// -ts_arkimex_initial_guess_extrapolate (StageGuessExtrapolate_).
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
      if (S_.stage_guess == "extrapolate") { PetscCall(StageGuessExtrapolate_(ts)); }
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// stage_guess = "extrapolate": each implicit stage's Newton starts from
   /// the dense output of the last accepted step, evaluated at the stage
   /// time, instead of from the stage before.  PETSc does it (arkimex.c,
   /// TSExtrapolate_ARKIMEX) and falls back to the stage before on the first
   /// step, after TSRestartStep and on the retry of a step whose stage
   /// solve failed.  Only the start of Newton moves: the stage equations,
   /// the convergence test and the error control are the same.
   ///
   /// PETSc has no call for it, only the option, and the option is read by
   /// TSSetFromOptions while the type is arkimex.  The caller's
   /// TSSetFromOptions ran before SetupTS set the type, and running it again
   /// would read the options file's -ts_type again.  So the option is set
   /// and read under a prefix of its own that no other option carries: the
   /// second TSSetFromOptions finds this one option and leaves every other
   /// setting at its current value; the TS's prefix is put back after.
   PetscErrorCode StageGuessExtrapolate_(TS ts)
   {
      PetscFunctionBeginUser;
      const char *old = nullptr;
      PetscCall(TSGetOptionsPrefix(ts, &old));
      const std::string keep = old ? old : "";
      PetscCall(PetscOptionsSetValue(nullptr,
                                     "-qdimp_guess_ts_arkimex_initial_guess_extrapolate",
                                     "true"));
      PetscCall(TSSetOptionsPrefix(ts, "qdimp_guess_"));
      PetscCall(TSSetFromOptions(ts));
      PetscCall(TSSetOptionsPrefix(ts, keep.empty() ? nullptr : keep.c_str()));
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
      bool local = false;               ///< evaluated at V_loc (newton_local)
      std::vector<real_t> D, Dt, fs, p, c, A, gamma;
      std::vector<real_t> keff;         ///< max(Gt_jj, 0) / shift, or 0
      std::vector<char> open;
      std::vector<int> O;
   };

   /// Owned point o: its own shear stiffness over the shift, the weight of
   /// its own slip in its stage equation (0 where Gt_jj <= 0).
   real_t Keff_(int o, real_t shift) const
   {
      const int j = P_.lo + o;
      return std::max(P_.Gt[size_t(j) * P_.N + j], real_t(0.0)) / shift;
   }

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
      if (S_.dense_lu == "shared")
      {
         // the entries of the loops below, one at a time: row i, column l
         auto entry = [this, N](int i, int l) -> real_t
         {
            const int j = (i < N) ? i : ev_.O[i - N];   // the row's point
            if (i < N)
            {
               const real_t fs = ev_.fs[j];
               if (l < N)
               {
                  real_t v = P_.Gt[size_t(j) * N + l] + fs * P_.Gs[size_t(j) * N + l];
                  if (l == j) { v += ev_.Dt[j]; }
                  return v;
               }
               const int lo = ev_.O[l - N];
               return cols_t_[slot_[lo]][j] + ((lo == j) ? 0.0 : fs * cols_s_[slot_[lo]][j]);
            }
            if (l < N) { return P_.Gs[size_t(j) * N + l]; }
            return cols_s_[slot_[ev_.O[l - N]]][j];
         };
         if (!plu_) { plu_ = std::make_unique<ParDenseLU>(comm_); }
         const bool ok = plu_->Factor(M, entry);
         fa_ = ev_;
         fa_.solve = solve_id_;
         have_factor_ = ok;
         slow_ = false;
         counts_.factors++;
         factor_singular_ = !ok;
         return;
      }
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
            // Ggs~: no own term (ds_j/dg_j = 0; fs_l is 0 on the open branch anyway)
            a[size_t(N + k) * M + j] = cols_t_[slot_[l]][j] +
                                       ((l == j) ? 0.0 : fs * cols_s_[slot_[l]][j]);
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
      QDImplicitFault *self = static_cast<QDImplicitFault *>(ctx);
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
      const real_t *br = reinterpret_cast<const real_t *>(bb);
      real_t *yr = reinterpret_cast<real_t *>(yy);
      if (self->slope_mode_ || !self->S_.newton_local) { self->ApplyInverse(br, yr); }
      else
      {
         // Newton on G: b is F at the iterate of the last Jacobian
         // evaluation (SNES newtonls: function, Jacobian, solve), so its
         // slip and psi rows are replaced by G's from that evaluation.  If
         // the order ever differs (another residual in between, or another
         // U), the update falls back to F with the same factor and is
         // counted; every rank decides alike.
         int ok = (self->loc_ok_ && self->loc_serial_ == self->res_serial_ && self->fa_.local)
                  ? 1 : 0;
         PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, self->comm_));
         if (ok)
         {
            const int n = self->P_.n;
            std::vector<real_t> g(br, br + size_t(kPer) * n);
            for (int o = 0; o < n; o++)
            {
               g[kPer * o] = self->loc_d_[o];
               g[kPer * o + 1] = self->loc_p_[o];
            }
            self->ApplyLocalInverse(g.data(), yr);
         }
         else
         {
            self->counts_.plain_applies++;
            self->ApplyInverse(br, yr);
         }
      }
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
         // dV/dpsi at the point: -p (the law), -p D / (D + keff) (its stage solve)
         const real_t pl = (ev.keff[j] > 0.0) ? ev.p[j] * ev.D[j] / (ev.D[j] + ev.keff[j])
                           : ev.p[j];
         const real_t sp = ev.shift - ev.c[j] * pl;
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
   std::vector<real_t> res_U_;                   ///< owned: the last residual's U
   long res_serial_ = 0, loc_serial_ = -1;       ///< residuals so far; at the last Jacobian
   std::vector<real_t> loc_d_, loc_p_;           ///< owned: G's slip and psi rows there
   bool loc_ok_ = false;                         ///< they belong to the last residual
   Snapshot ev_, fa_;                            ///< last evaluation, factor's
   DenseMatrix lu_mat_;
   std::vector<int> lu_piv_;
   LUFactors lu_;
   std::unique_ptr<ParDenseLU> plu_;   ///< dense_lu = "shared"
   bool have_factor_ = false, factor_singular_ = false, bad_jacobian_ = false;
   bool slope_mode_ = false, slow_ = false, clipped_ = false, prev_clipped_ = false;
   long solve_id_ = 0;
   double wprev_ = -1.0, last_w_ = 0.0;
   Counts counts_;
   Mat J_ = nullptr;
};

/// Option I, phase I2 (P5): the owned fault data of the planar problem, in
/// the QD kernel's layouts.
struct PlanarFaultInputs
{
   Vector tau_pre;        ///< 2 per owned DOF (the kernel's TauPre; in-plane at the slot)
   Vector sigma0;         ///< 1 per owned DOF, compression > 0
   Vector a, eta, Dc;     ///< 1 per owned DOF
   real_t b = 0.0, f0 = 0.0, V0 = 0.0;
   real_t sigma_floor = 0.0;    ///< [cycle].sigma_n_floor_pa
   bool floor_taper = false;    ///< [cycle.opening].floor_taper
};

/// The pure traction response of one owned stage-slip vector e (3 entries
/// per owned DOF; entry 3 o + 2 the opening): T (2 per owned DOF), T_n (1
/// per owned DOF).  ElasticityDomainOperator::QDUnitColumn.  Collective.
using UnitResponse = std::function<void(const Vector &e, Vector &T, Vector &Tn)>;

/// Option I, phase I2: the problem of QDImplicitFault on the 2-D planar
/// dense Green's function @a G (P5: refuses the H-matrix, 3-D and a G
/// without normal rows).  Gt, Gs, TDt, TDs: G's in-plane shear and normal
/// rows (in-plane slip columns) and its loading rates, gathered on every
/// rank (QDGreenOperator::GatherPlanar; the global point order is G's, rank
/// by rank).  tau_pre: the in-plane slot of the kernel's pre-stress (an
/// out-of-plane pre-stress is refused: option I has one shear component).
/// gap_column(l): @a unit_g of a unit opening at global point l (zero
/// elsewhere, so no loading enters), its in-plane shear and normal parts
/// gathered on every rank; it refuses a nonzero out-of-plane shear.
/// Collective.
inline QDImplicitFault::Problem PlanarProblem(MPI_Comm comm, const QDGreenOperator &G,
                                              const PlanarFaultInputs &in,
                                              UnitResponse unit_g)
{
   MFEM_VERIFY(G.Planar() && G.NComp() == 3 && !G.Compressed(), "option I (P5): the 2-D "
               "planar dense Green's function with its normal rows only");
   int rank = 0, nranks = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nranks);
   const int n = G.NumOwned(), slot = G.PlanarSlot();
   long long lo = 0, mine = n;
   MPI_Exscan(&mine, &lo, 1, MPI_LONG_LONG, MPI_SUM, comm);
   if (rank == 0) { lo = 0; }
   MFEM_VERIFY(lo == G.FirstOwned(), "option I: the Green's function's owned range starts "
               "at " << G.FirstOwned() << ", not at the rank-order offset " << lo);
   MFEM_VERIFY(in.tau_pre.Size() == 2 * n && in.sigma0.Size() == n && in.a.Size() == n &&
               in.eta.Size() == n && in.Dc.Size() == n, "option I: PlanarFaultInputs "
               "sizes (tau_pre 2 n, the rest n; n = " << n << ")");
   real_t out = 0.0;
   for (int o = 0; o < n; o++) { out = std::max(out, std::abs(in.tau_pre(2 * o + 1 - slot))); }
   MPI_Allreduce(MPI_IN_PLACE, &out, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MFEM_VERIFY(out == 0.0, "option I: an out-of-plane shear pre-stress (largest " << out
               << " Pa); the implicit fault has the in-plane shear only");
   MFEM_VERIFY(unit_g, "option I: no unit opening response");

   QDImplicitFault::Problem P;
   P.N = static_cast<int>(G.NumFaultDofs());
   G.GatherPlanar(P.Gt, P.Gs, P.TDt, P.TDs);
   P.lo = static_cast<int>(lo);
   P.n = n;
   P.tau_pre.resize(n);
   P.sigma0.resize(n);
   P.a.resize(n);
   P.eta.resize(n);
   P.Dc.resize(n);
   for (int o = 0; o < n; o++)
   {
      P.tau_pre[o] = in.tau_pre(2 * o + slot);
      P.sigma0[o] = in.sigma0(o);
      P.a[o] = in.a(o);
      P.eta[o] = in.eta(o);
      P.Dc[o] = in.Dc(o);
   }
   P.b = in.b;
   P.f0 = in.f0;
   P.V0 = in.V0;
   P.sigma_floor = in.sigma_floor;
   P.floor_taper = in.floor_taper;

   std::vector<int> cnt(nranks), dsp(nranks);
   MPI_Allgather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
   for (int r = 0, d = 0; r < nranks; r++) { dsp[r] = d; d += cnt[r]; }
   const int N = P.N, lo_i = P.lo;
   P.gap_column = [comm, N, n, lo_i, slot, cnt, dsp, unit_g](int l, std::vector<real_t> &t,
                                                             std::vector<real_t> &s)
   {
      Vector e(3 * n), T, Tn;
      e = 0.0;
      if (lo_i <= l && l < lo_i + n) { e(3 * (l - lo_i) + 2) = 1.0; }
      unit_g(e, T, Tn);
      MFEM_VERIFY(T.Size() == 2 * n && Tn.Size() == n, "option I: the unit opening "
                  "response has " << T.Size() << " shear and " << Tn.Size()
                  << " normal entries (n = " << n << ")");
      std::vector<real_t> mt(n), ms(n);
      real_t big[2] = {0.0, 0.0};   // out-of-plane, in-plane
      for (int o = 0; o < n; o++)
      {
         mt[o] = T(2 * o + slot);
         ms[o] = Tn(o);
         big[0] = std::max(big[0], std::abs(T(2 * o + 1 - slot)));
         big[1] = std::max(big[1], std::max(std::abs(mt[o]), std::abs(ms[o])));
      }
      MPI_Allreduce(MPI_IN_PLACE, big, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      MFEM_VERIFY(big[0] <= 1e-12 * big[1], "option I: the opening at point " << l
                  << " drives an out-of-plane shear of " << big[0] << " Pa/m (in-plane and "
                  "normal up to " << big[1] << " Pa/m)");
      t.resize(N);
      s.resize(N);
      const MPI_Datatype R = MPITypeMap<real_t>::mpi_type;
      MPI_Allgatherv(mt.data(), n, R, t.data(), cnt.data(), dsp.data(), R, comm);
      MPI_Allgatherv(ms.data(), n, R, s.data(), cnt.data(), dsp.data(), R, comm);
   };
   return P;
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_USE_PETSC && MFEM_USE_MPI
#endif // MFEM_SEAS_CYCLE_QD_IMPLICIT_HPP
