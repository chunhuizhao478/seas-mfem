// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// ts_rock4.hpp — QD speed S3: ROCK4 (A. Abdulle, SIAM J. Sci. Comput. 23(6),
// 2041-2054, 2002) as a PETSc TS type, "-ts_type rock4" (opt-in; the QD
// default stays -ts_type rk -ts_rk_type 5dp).
//
// ROCK4 is an explicit Runge-Kutta method of order 4 whose stability region
// grows with the square of its stage count: with s stages it is stable for
// h·rho(J) <= about 0.35 s².  Each step takes as many stages as h·rho needs,
// so on a stiff problem whose spectrum lies near the negative real axis (the
// QD fault between earthquakes: S1 found real negative eigenvalues and an
// RK5DP step on the stability limit) the step is no longer set by rho.  A
// step of s stages costs s right-hand sides.
//
// One step from (t, y) with step h (rock4.f, and OrdinaryDiffEq.jl's ROCK4):
//   1. rho: the spectral radius of J = df/dy, by a power iteration on f
//      around y (each iteration one right-hand side), every
//      -ts_rock4_eig_interval steps and after a rejected step.
//   2. Stage count: m = floor(sqrt((3 + h·rho)/0.353)) + 1, at least 5
//      (-ts_rock4_min_stages) and at most -ts_rock4_max_stages (124); if h·rho
//      needs more, h is cut to 0.8·(0.353·max² - 3)/rho.  The method used is
//      the first of rock4_tables.hpp's 50 whose Chebyshev degree
//      deg >= m - 4: it has deg + 4 stages.  Above 24 stages the methods
//      come in steps (26, 28, ..., 124), so the cap is rounded down to a
//      method's stage count and no step takes more stages than asked.
//      Why 124 and not the tables' 152: two of the 50 methods, with 133 and
//      152 stages (degree 129 and 148), are not stable on their whole
//      interval; |R(z)| reaches 1.0059 and 1.0077 near z = -8.4 (found
//      10-25 by stepping y' = z y with them and in exact rational arithmetic
//      from the tables; test_ts_rock4 keeps both as its negative control).
//      The other 48 keep |R(z)| < 1 on their whole interval (at most
//      0.9984 away from z = 0); capping at 124 stages (h·rho up to about
//      5,400 per step) leaves both out with one rule.
//   3. Chebyshev part (deg stages, three-term recurrence, mu and kappa from
//      kRock4Recf): u_0 = y, u_1 = y + h mu_1 f(t, y),
//      u_j = h mu_j f(t_{j-1}, u_{j-1}) + (1 + kappa_j) u_{j-1} - kappa_j
//      u_{j-2}, with the stage times t_j by the same recurrence.
//   4. Finishing procedure: a 4-stage explicit RK from u_deg (kRock4Fpa,
//      kRock4Fpb) gives y_{n+1} (order 4); its embedded weights
//      (kRock4Fpbe, the fifth on f(t + h, y_{n+1})) give an order-3
//      solution, and their difference is the error estimate.
//   5. f(t + h, y_{n+1}) is evaluated in every attempt (it is part of the
//      error estimate), so the last right-hand side of an accepted step is at
//      the accepted state (first same as last): the next step starts from it,
//      and whatever the right-hand side caches (the BP5 monitor's V_max and
//      traction) belongs to the accepted state, as with RK 5dp.
//   6. Step control by the TSAdapt in use (basic by default), with the order
//      4 and error estimate of 4.
//
// Right-hand sides per accepted step: deg + 4, plus the power iteration's
// share.  Options:
//   -ts_rock4_max_stages <124>   at most this many stages (5..124)
//   -ts_rock4_min_stages <5>     at least this many stages (5..max)
//   -ts_rock4_eig_interval <25>  steps between spectral radius estimates
//   -ts_rock4_eig_safety <1.2>   factor on the converged estimate
//
// Uses PETSc's private TS header (the TS type interface), as every TS
// implementation does; written against PETSc 3.24.

#ifndef MFEM_SEAS_TS_ROCK4_HPP
#define MFEM_SEAS_TS_ROCK4_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_PETSC

#include "rock4_tables.hpp"

#include <petsc/private/tsimpl.h>

#include <algorithm>
#include <cmath>

namespace mfem
{
namespace seas
{

/// Counts of a rock4 TS since it was created (TSRock4GetStats).
struct Rock4Stats
{
   long long steps = 0;          ///< accepted steps
   long long rejected = 0;       ///< rejected attempts
   long long f_stage = 0;        ///< right-hand sides in the steps
   long long f_eig = 0;          ///< right-hand sides in power iterations
   long long eig_estimates = 0;  ///< power iterations run
   long long eig_unconverged = 0;///< of those, stopped at 50 iterations
   long long stages_sum = 0;     ///< stages of the accepted steps
   int stages_max = 0;           ///< most stages in one accepted step
   long long h_cut = 0;          ///< attempts whose h the stage cap cut
   double eig_last = 0.0;        ///< last spectral radius estimate
};

/// Most stages any step may take: the 124-stage method (degree 120) is the
/// last one used, as the 133- and 152-stage methods after it are not stable
/// on their whole interval (see the top of this file).
constexpr int kRock4MaxStages = 124;

namespace rock4_detail
{

struct TSRock4
{
   Vec F0 = nullptr;      ///< f(t_n, y_n)
   Vec F1 = nullptr;      ///< f(t_n + h, y_{n+1}) of the attempt
   Vec Fprev = nullptr;   ///< f(t_{n-1}, y_{n-1}) after a step (interpolation)
   Vec Yprev = nullptr;   ///< y_{n-1} after a step (rollback, interpolation)
   Vec U[3] = {nullptr, nullptr, nullptr};  ///< Chebyshev recurrence
   Vec Y2 = nullptr, Y3 = nullptr, Y4 = nullptr;  ///< finishing stages
   Vec Ynew = nullptr;    ///< y_{n+1} of the attempt
   Vec E = nullptr;       ///< order-3 minus order-4 solution
   Vec K = nullptr;       ///< a stage's right-hand side
   Vec Z = nullptr;       ///< power iteration vector
   Vec Zdir = nullptr;    ///< last power iteration direction
   PetscBool have_f0 = PETSC_FALSE;
   PetscBool have_zdir = PETSC_FALSE;
   PetscBool have_eig = PETSC_FALSE;
   PetscBool last_rejected = PETSC_FALSE;
   PetscReal eig = 0.0;
   PetscInt eig_age = 0;
   PetscInt eig_interval = 25;
   PetscReal eig_safety = 1.2;
   PetscInt max_stages = 124;
   PetscInt min_stages = 5;
   PetscInt last_stages = 0;
   TSStepStatus status = TS_STEP_COMPLETE;
   Rock4Stats stats;
};

inline PetscErrorCode RHS(TS ts, TSRock4 *rk, PetscReal t, Vec X, Vec F,
                          bool eig)
{
   PetscFunctionBegin;
   PetscCall(TSComputeRHSFunction(ts, t, X, F));
   if (eig) { rk->stats.f_eig++; }
   else { rk->stats.f_stage++; }
   PetscFunctionReturn(PETSC_SUCCESS);
}

/// Spectral radius of df/dy at (t, y = ts->vec_sol) by the power iteration
/// of rock4.f's rho (as OrdinaryDiffEq.jl's maxeig! for ROCK4): difference
/// quotients of f around y, at most 50 iterations, converged when two
/// estimates agree to 1 % of the larger; then times eig_safety.
inline PetscErrorCode EstimateEig(TS ts, TSRock4 *rk, PetscReal t)
{
   PetscFunctionBegin;
   Vec y = ts->vec_sol;
   PetscReal ynorm, znorm;
   PetscCall(VecNorm(y, NORM_2, &ynorm));
   PetscCall(VecCopy(rk->have_zdir ? rk->Zdir : rk->F0, rk->Z));
   PetscCall(VecNorm(rk->Z, NORM_2, &znorm));
   const PetscReal sqrtu = std::sqrt(PETSC_MACHINE_EPSILON);
   PetscReal dz;
   if (ynorm > 0.0 && znorm > 0.0)
   {
      dz = ynorm * sqrtu;
      PetscCall(VecAYPX(rk->Z, dz / znorm, y));    // Z = y + (dz/|z|) z
   }
   else if (ynorm > 0.0)
   {
      dz = ynorm * sqrtu;
      PetscCall(VecAXPBY(rk->Z, 1.0 + sqrtu, 0.0, y));
   }
   else if (znorm > 0.0)
   {
      dz = sqrtu;
      PetscCall(VecScale(rk->Z, dz / znorm));
   }
   else
   {
      dz = sqrtu;
      PetscInt n;
      PetscCall(VecGetSize(rk->Z, &n));
      PetscCall(VecSet(rk->Z, dz / std::sqrt(static_cast<PetscReal>(n))));
   }
   PetscReal est = 0.0;
   PetscBool converged = PETSC_FALSE;
   const PetscReal dt_max = ts->adapt ? ts->adapt->dt_max : PETSC_MAX_REAL;
   for (int it = 1; it <= 50; it++)
   {
      PetscCall(RHS(ts, rk, t, rk->Z, rk->K, true));
      PetscCall(VecAXPY(rk->K, -1.0, rk->F0));     // K = f(Z) - f(y)
      PetscReal delta;
      PetscCall(VecNorm(rk->K, NORM_2, &delta));
      const PetscReal prev = est;
      est = delta / dz;
      if (it >= 2 && std::abs(prev - est) <
          0.01 * std::max(est, static_cast<PetscReal>(1.0 / dt_max)))
      {
         converged = PETSC_TRUE;
         PetscCall(VecWAXPY(rk->Zdir, -1.0, y, rk->Z));
         rk->have_zdir = PETSC_TRUE;
         break;
      }
      if (delta > 0.0)
      {
         PetscCall(VecWAXPY(rk->Z, dz / delta, rk->K, y));  // Z = y + K dz/|K|
      }
      else
      {
         // f did not change: perturb another direction
         PetscCall(VecAXPY(rk->Z, -1.0, y));
         PetscCall(VecScale(rk->Z, -1.0));
         PetscCall(VecAXPY(rk->Z, 1.0, y));
      }
   }
   rk->stats.eig_estimates++;
   if (!converged)
   {
      rk->stats.eig_unconverged++;
      PetscCall(VecWAXPY(rk->Zdir, -1.0, y, rk->Z));
      rk->have_zdir = PETSC_TRUE;
   }
   rk->eig = rk->eig_safety * est;
   rk->have_eig = PETSC_TRUE;
   rk->eig_age = 0;
   rk->stats.eig_last = rk->eig;
   PetscFunctionReturn(PETSC_SUCCESS);
}

/// Largest stage count of a method (degree + 4) that is <= @a stages.
inline int RoundDownToMethod(int stages)
{
   int best = kRock4Ms[0] + 4;
   for (int i = 0; i < kRock4Methods; i++)
   {
      if (kRock4Ms[i] + 4 <= stages) { best = kRock4Ms[i] + 4; }
   }
   return best;
}

/// Index of the first method of rock4_tables.hpp with degree >= @a m, and
/// the offset of its recurrence coefficients in kRock4Recf.
inline void ChooseMethod(int m, int &index, int &start)
{
   start = 0;
   for (int i = 0; i < kRock4Methods; i++)
   {
      if (kRock4Ms[i] >= m) { index = i; return; }
      start += 2 * kRock4Ms[i] - 1;
   }
   index = kRock4Methods - 1;
   start -= 2 * kRock4Ms[index] - 1;
}

inline PetscErrorCode TSStep_Rock4(TS ts)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   TSAdapt adapt;
   PetscBool accept = PETSC_TRUE;
   PetscInt rejections = 0;
   PetscReal next_h = ts->time_step;

   if (ts->steprollback || ts->steprestart) { rk->have_f0 = PETSC_FALSE; }
   if (!rk->have_f0)
   {
      PetscCall(RHS(ts, rk, ts->ptime, ts->vec_sol, rk->F0, false));
      rk->have_f0 = PETSC_TRUE;
   }
   rk->status = TS_STEP_INCOMPLETE;
   while (!ts->reason && rk->status != TS_STEP_COMPLETE)
   {
      const PetscReal t = ts->ptime;
      if (!rk->have_eig || rk->eig_age >= rk->eig_interval ||
          rk->last_rejected)
      {
         PetscCall(EstimateEig(ts, rk, t));
      }
      // Stage count for h, or h for the largest stage count
      PetscReal h = ts->time_step;
      PetscReal mreal = std::floor(std::sqrt((3.0 + h * rk->eig) / 0.353))
                        + 1.0;
      if (mreal > rk->max_stages)
      {
         h = 0.8 * (0.353 * rk->max_stages * rk->max_stages - 3.0) / rk->eig;
         ts->time_step = h;
         mreal = rk->max_stages;
         rk->stats.h_cut++;
      }
      const int m = static_cast<int>(std::max<PetscReal>(mreal,
                                                         rk->min_stages));
      int mi, start;
      ChooseMethod(std::max(m, 5) - 4, mi, start);
      const int deg = kRock4Ms[mi];
      rk->last_stages = deg + 4;
      const double *recf = kRock4Recf + start;
      Vec y = ts->vec_sol;
      PetscBool stageok = PETSC_TRUE;

      // Chebyshev part: uold = y, ucur = y + h mu_1 f(t, y)
      Vec uold = rk->U[0], ucur = rk->U[1], unext = rk->U[2];
      PetscReal t_old = t, t_cur = t + h * recf[0];
      PetscCall(VecCopy(y, uold));
      PetscCall(VecWAXPY(ucur, h * recf[0], rk->F0, y));
      for (int j = 2; j <= deg && stageok; j++)
      {
         const PetscReal mu = recf[2 * (j - 2) + 1];
         const PetscReal kappa = recf[2 * (j - 2) + 2];
         PetscCall(TSPreStage(ts, t_cur));
         PetscCall(RHS(ts, rk, t_cur, ucur, rk->K, false));
         // unext = h mu K + (1 + kappa) ucur - kappa uold
         PetscCall(VecAXPBYPCZ(unext, h * mu, 1.0 + kappa, 0.0, rk->K,
                               ucur));
         PetscCall(VecAXPY(unext, -kappa, uold));
         const PetscReal t_next = h * mu + (1.0 + kappa) * t_cur -
                                  kappa * t_old;
         t_old = t_cur;
         t_cur = t_next;
         Vec tmp = uold;
         uold = ucur;
         ucur = unext;
         unext = tmp;
         PetscCall(TSPostStage(ts, t_cur, j - 1, &ucur));
         PetscCall(TSGetAdapt(ts, &adapt));
         PetscCall(TSAdaptCheckStage(adapt, ts, t_cur, ucur, &stageok));
      }
      if (stageok)
      {
         // Finishing procedure from acc = ucur at t_cur
         Vec acc = ucur;
         const double *a = kRock4Fpa[mi];
         const double *b = kRock4Fpb[mi];
         const double *be = kRock4Fpbe[mi];
         const PetscReal c2 = h * a[0];
         const PetscReal c3 = h * (a[1] + a[2]);
         const PetscReal c4 = h * (a[3] + a[4] + a[5]);
         PetscCall(TSPreStage(ts, t_cur));
         PetscCall(RHS(ts, rk, t_cur, acc, rk->K, false));
         PetscCall(VecWAXPY(rk->Y2, h * a[0], rk->K, acc));
         PetscCall(VecWAXPY(rk->Y3, h * a[1], rk->K, acc));
         PetscCall(VecWAXPY(rk->Y4, h * a[3], rk->K, acc));
         PetscCall(VecWAXPY(rk->Ynew, h * b[0], rk->K, acc));
         PetscCall(VecAXPBY(rk->E, h * (be[0] - b[0]), 0.0, rk->K));
         PetscCall(TSPreStage(ts, t_cur + c2));
         PetscCall(RHS(ts, rk, t_cur + c2, rk->Y2, rk->K, false));
         PetscCall(VecAXPY(rk->Y3, h * a[2], rk->K));
         PetscCall(VecAXPY(rk->Y4, h * a[4], rk->K));
         PetscCall(VecAXPY(rk->Ynew, h * b[1], rk->K));
         PetscCall(VecAXPY(rk->E, h * (be[1] - b[1]), rk->K));
         PetscCall(TSPreStage(ts, t_cur + c3));
         PetscCall(RHS(ts, rk, t_cur + c3, rk->Y3, rk->K, false));
         PetscCall(VecAXPY(rk->Y4, h * a[5], rk->K));
         PetscCall(VecAXPY(rk->Ynew, h * b[2], rk->K));
         PetscCall(VecAXPY(rk->E, h * (be[2] - b[2]), rk->K));
         PetscCall(TSPreStage(ts, t_cur + c4));
         PetscCall(RHS(ts, rk, t_cur + c4, rk->Y4, rk->K, false));
         PetscCall(VecAXPY(rk->Ynew, h * b[3], rk->K));
         PetscCall(VecAXPY(rk->E, h * (be[3] - b[3]), rk->K));
         PetscCall(TSGetAdapt(ts, &adapt));
         PetscCall(TSAdaptCheckStage(adapt, ts, t + h, rk->Ynew, &stageok));
      }
      if (stageok)
      {
         // Last: f at the new solution (error estimate, next step's F0)
         PetscCall(TSPreStage(ts, t + h));
         PetscCall(RHS(ts, rk, t + h, rk->Ynew, rk->F1, false));
         PetscCall(VecAXPY(rk->E, h * kRock4Fpbe[mi][4], rk->F1));

         rk->status = TS_STEP_PENDING;
         PetscCall(TSGetAdapt(ts, &adapt));
         PetscCall(TSAdaptCandidatesClear(adapt));
         PetscCall(TSAdaptCandidateAdd(adapt, "rock4", 4, 1, 1.0,
                                       static_cast<PetscReal>(deg + 4),
                                       PETSC_TRUE));
         PetscCall(TSAdaptChoose(adapt, ts, ts->time_step, nullptr, &next_h,
                                 &accept));
      }
      else { accept = PETSC_FALSE; }
      rk->eig_age++;
      if (!accept)
      {
         rk->status = TS_STEP_INCOMPLETE;
         rk->stats.rejected++;
         rk->last_rejected = PETSC_TRUE;
         if (stageok) { ts->time_step = next_h; }
         ts->reject++;
         if (!ts->reason && ++rejections > ts->max_reject &&
             ts->max_reject >= 0)
         {
            ts->reason = TS_DIVERGED_STEP_REJECTED;
            PetscCall(PetscInfo(ts, "Step=%" PetscInt_FMT ", step rejections "
                                "%" PetscInt_FMT " greater than current TS "
                                "allowed, stopping solve\n", ts->steps,
                                rejections));
         }
         continue;
      }
      // Accept: y_{n-1}, f_{n-1} kept for rollback and interpolation
      PetscCall(VecCopy(y, rk->Yprev));
      PetscCall(VecCopy(rk->Ynew, y));
      Vec f = rk->Fprev;
      rk->Fprev = rk->F0;
      rk->F0 = rk->F1;
      rk->F1 = f;
      rk->status = TS_STEP_COMPLETE;
      rk->last_rejected = PETSC_FALSE;
      rk->stats.steps++;
      rk->stats.stages_sum += deg + 4;
      rk->stats.stages_max = std::max(rk->stats.stages_max, deg + 4);
      ts->ptime += ts->time_step;
      ts->time_step = next_h;
   }
   PetscFunctionReturn(PETSC_SUCCESS);
}

/// The order-4 error estimate of the attempt (TSAdapt basic and dsp).
inline PetscErrorCode TSEvaluateWLTE_Rock4(TS ts, NormType wnormtype,
                                           PetscInt *order, PetscReal *wlte)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   PetscReal wltea, wlter;
   PetscCall(VecWAXPY(rk->K, 1.0, rk->Ynew, rk->E));   // order-3 solution
   PetscCall(TSErrorWeightedNorm(ts, rk->Ynew, rk->K, wnormtype, wlte,
                                 &wltea, &wlter));
   if (order) { *order = 4; }
   PetscFunctionReturn(PETSC_SUCCESS);
}

/// The attempt's order-4 solution (the order 3 one with @a order 3).
inline PetscErrorCode TSEvaluateStep_Rock4(TS ts, PetscInt order, Vec X,
                                           PetscBool *done)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   if (rk->status == TS_STEP_COMPLETE && order == 4)
   {
      PetscCall(VecCopy(ts->vec_sol, X));
   }
   else if (order == 4) { PetscCall(VecCopy(rk->Ynew, X)); }
   else if (order == 3 && rk->status != TS_STEP_COMPLETE)
   {
      PetscCall(VecWAXPY(X, 1.0, rk->Ynew, rk->E));
   }
   else
   {
      PetscCheck(done, PetscObjectComm((PetscObject)ts), PETSC_ERR_SUP,
                 "rock4 cannot evaluate the step at order %" PetscInt_FMT,
                 order);
      *done = PETSC_FALSE;
      PetscFunctionReturn(PETSC_SUCCESS);
   }
   if (done) { *done = PETSC_TRUE; }
   PetscFunctionReturn(PETSC_SUCCESS);
}

inline PetscErrorCode TSRollBack_Rock4(TS ts)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   PetscCall(VecCopy(rk->Yprev, ts->vec_sol));
   Vec f = rk->F0;
   rk->F0 = rk->Fprev;
   rk->Fprev = f;
   rk->have_f0 = PETSC_TRUE;
   PetscFunctionReturn(PETSC_SUCCESS);
}

/// Cubic Hermite interpolation on the last step (y and f at both ends).
inline PetscErrorCode TSInterpolate_Rock4(TS ts, PetscReal itime, Vec X)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   const PetscReal t0 = ts->ptime_prev, h = ts->ptime - ts->ptime_prev;
   PetscCheck(h > 0.0, PetscObjectComm((PetscObject)ts), PETSC_ERR_ORDER,
              "rock4: no completed step to interpolate on");
   const PetscReal s = (itime - t0) / h;
   const PetscReal h00 = (1 + 2 * s) * (1 - s) * (1 - s);
   const PetscReal h10 = s * (1 - s) * (1 - s);
   const PetscReal h01 = s * s * (3 - 2 * s);
   const PetscReal h11 = s * s * (s - 1);
   PetscCall(VecAXPBYPCZ(X, h00, h01, 0.0, rk->Yprev, ts->vec_sol));
   PetscCall(VecAXPBYPCZ(X, h * h10, h * h11, 1.0, rk->Fprev, rk->F0));
   PetscFunctionReturn(PETSC_SUCCESS);
}

inline PetscErrorCode TSSetUp_Rock4(TS ts)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   Vec *all[] = {&rk->F0, &rk->F1, &rk->Fprev, &rk->Yprev, &rk->U[0],
                 &rk->U[1], &rk->U[2], &rk->Y2, &rk->Y3, &rk->Y4, &rk->Ynew,
                 &rk->E, &rk->K, &rk->Z, &rk->Zdir
                };
   for (Vec *v : all)
   {
      if (!*v) { PetscCall(VecDuplicate(ts->vec_sol, v)); }
   }
   PetscCall(TSGetAdapt(ts, &ts->adapt));
   PetscCall(TSAdaptCandidatesClear(ts->adapt));
   PetscCall(TSAdaptCandidateAdd(ts->adapt, "rock4", 4, 1, 1.0, 5.0,
                                 PETSC_TRUE));
   rk->have_f0 = PETSC_FALSE;
   rk->have_eig = PETSC_FALSE;
   rk->have_zdir = PETSC_FALSE;
   PetscFunctionReturn(PETSC_SUCCESS);
}

inline PetscErrorCode TSReset_Rock4(TS ts)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   Vec *all[] = {&rk->F0, &rk->F1, &rk->Fprev, &rk->Yprev, &rk->U[0],
                 &rk->U[1], &rk->U[2], &rk->Y2, &rk->Y3, &rk->Y4, &rk->Ynew,
                 &rk->E, &rk->K, &rk->Z, &rk->Zdir
                };
   for (Vec *v : all) { PetscCall(VecDestroy(v)); }
   rk->have_f0 = PETSC_FALSE;
   rk->have_eig = PETSC_FALSE;
   rk->have_zdir = PETSC_FALSE;
   PetscFunctionReturn(PETSC_SUCCESS);
}

inline PetscErrorCode TSDestroy_Rock4(TS ts)
{
   PetscFunctionBegin;
   PetscCall(TSReset_Rock4(ts));
   delete static_cast<TSRock4 *>(ts->data);
   ts->data = nullptr;
   PetscFunctionReturn(PETSC_SUCCESS);
}

/// The options argument of a TS setfromoptions: PetscOptionItems from
/// PETSc 3.23 on, PetscOptionItems * before.
template <typename F> struct OptionsArg;
template <typename R, typename A, typename B> struct OptionsArg<R (*)(A, B)>
{
   using type = B;
};
using OptionItemsArg = OptionsArg<decltype(_TSOps::setfromoptions)>::type;

inline PetscErrorCode TSSetFromOptions_Rock4(TS ts,
                                             OptionItemsArg PetscOptionsObject)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   PetscOptionsHeadBegin(PetscOptionsObject, "ROCK4 options");
   PetscCall(PetscOptionsInt("-ts_rock4_max_stages", "Most stages per step "
                             "(5..124)", "", rk->max_stages, &rk->max_stages,
                             nullptr));
   PetscCall(PetscOptionsInt("-ts_rock4_min_stages", "Fewest stages per step "
                             "(5..max)", "", rk->min_stages, &rk->min_stages,
                             nullptr));
   PetscCall(PetscOptionsInt("-ts_rock4_eig_interval", "Steps between "
                             "spectral radius estimates", "",
                             rk->eig_interval, &rk->eig_interval, nullptr));
   PetscCall(PetscOptionsReal("-ts_rock4_eig_safety", "Factor on the "
                              "spectral radius estimate", "", rk->eig_safety,
                              &rk->eig_safety, nullptr));
   PetscOptionsHeadEnd();
   PetscCheck(rk->max_stages >= 5 && rk->max_stages <= kRock4MaxStages,
              PetscObjectComm((PetscObject)ts), PETSC_ERR_ARG_OUTOFRANGE,
              "-ts_rock4_max_stages %" PetscInt_FMT ": need 5..%d",
              rk->max_stages, kRock4MaxStages);
   rk->max_stages = RoundDownToMethod(static_cast<int>(rk->max_stages));
   PetscCheck(rk->min_stages >= 5 && rk->min_stages <= rk->max_stages,
              PetscObjectComm((PetscObject)ts), PETSC_ERR_ARG_OUTOFRANGE,
              "-ts_rock4_min_stages %" PetscInt_FMT ": need 5..%" PetscInt_FMT,
              rk->min_stages, rk->max_stages);
   PetscCheck(rk->eig_interval >= 1, PetscObjectComm((PetscObject)ts),
              PETSC_ERR_ARG_OUTOFRANGE, "-ts_rock4_eig_interval %"
              PetscInt_FMT ": need >= 1", rk->eig_interval);
   PetscCheck(rk->eig_safety >= 1.0, PetscObjectComm((PetscObject)ts),
              PETSC_ERR_ARG_OUTOFRANGE, "-ts_rock4_eig_safety %g: need >= 1",
              (double)rk->eig_safety);
   PetscFunctionReturn(PETSC_SUCCESS);
}

inline PetscErrorCode TSView_Rock4(TS ts, PetscViewer viewer)
{
   PetscFunctionBegin;
   TSRock4 *rk = static_cast<TSRock4 *>(ts->data);
   PetscBool ascii;
   PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII,
                                    &ascii));
   if (ascii)
   {
      PetscCall(PetscViewerASCIIPrintf(viewer, "  ROCK4: stages %"
                                       PetscInt_FMT "..%" PetscInt_FMT
                                       ", spectral radius every %"
                                       PetscInt_FMT " steps, safety %g\n",
                                       rk->min_stages, rk->max_stages,
                                       rk->eig_interval,
                                       (double)rk->eig_safety));
   }
   PetscFunctionReturn(PETSC_SUCCESS);
}

inline PetscErrorCode TSCreate_Rock4(TS ts)
{
   PetscFunctionBegin;
   ts->ops->reset = TSReset_Rock4;
   ts->ops->destroy = TSDestroy_Rock4;
   ts->ops->view = TSView_Rock4;
   ts->ops->setup = TSSetUp_Rock4;
   ts->ops->step = TSStep_Rock4;
   ts->ops->interpolate = TSInterpolate_Rock4;
   ts->ops->evaluatewlte = TSEvaluateWLTE_Rock4;
   ts->ops->evaluatestep = TSEvaluateStep_Rock4;
   ts->ops->rollback = TSRollBack_Rock4;
   ts->ops->setfromoptions = TSSetFromOptions_Rock4;
   ts->default_adapt_type = TSADAPTBASIC;   // it has an error estimate
   ts->data = new TSRock4();
   PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace rock4_detail

/// Make "-ts_type rock4" known to PETSc.  Call after PetscInitialize and
/// before TSSetFromOptions; later calls do nothing.
inline void TSRock4Register()
{
   static bool done = false;
   if (done) { return; }
   PetscErrorCode ierr = TSRegister("rock4", rock4_detail::TSCreate_Rock4);
   MFEM_VERIFY(ierr == PETSC_SUCCESS, "TSRegister(rock4) failed");
   done = true;
}

/// Whether @a ts is of type rock4.
inline bool TSIsRock4(TS ts)
{
   PetscBool is = PETSC_FALSE;
   PetscErrorCode ierr = PetscObjectTypeCompare((PetscObject)ts, "rock4",
                                                &is);
   MFEM_VERIFY(ierr == PETSC_SUCCESS, "PetscObjectTypeCompare failed");
   return is == PETSC_TRUE;
}

/// Counts of a rock4 TS (TSIsRock4(ts) must hold).
inline Rock4Stats TSRock4GetStats(TS ts)
{
   MFEM_VERIFY(TSIsRock4(ts), "TSRock4GetStats: the TS is not rock4");
   return static_cast<rock4_detail::TSRock4 *>(ts->data)->stats;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_PETSC

#endif // MFEM_SEAS_TS_ROCK4_HPP
