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

// test_ts_rock4.cpp — QD speed S3: the PETSc TS type rock4 of
// solver/ts_rock4.hpp on the linear test problem
//
//   y_i' = lam_i (y_i - cos t) - sin t,   y_i(0) = 2,
//   y_i(t) = cos t + exp(lam_i t)          (exact)
//
// with N entries spread over the ranks, on a nonlinear one and on y' = z y.
// No mesh.
//
//   mpirun -np N ./seas_test_ts_rock4
//
//   1. Order: lam = -20 for all i, fixed steps (adapt none) h = 0.2, 0.1,
//      0.05, 0.025 to t = 1 with 20 stages in every step; the error falls
//      by at least 2^3.6 per halving (order 4 expected).
//   2. Stiff, adaptive: lam_i log-spaced in [-1e4, -1], atol 1e-8, to t = 10:
//      error below 1e-6, the spectral radius estimate within [1, 1.5]·1e4,
//      more than 5 stages in the largest step, and fewer right-hand sides
//      than rk 5dp at the same tolerance (both printed).
//   3. First same as last: after every accepted step the last right-hand
//      side was evaluated at the accepted time and state, bitwise.
//   4. Stage cap: with -ts_rock4_max_stages 12 the step size is cut to the
//      cap (counted) and the error stays below 1e-6.
//   5. Options outside their ranges stop with an error (min stages 4, max
//      stages 125); a cap between two methods (31) is rounded down to one
//      (30 stages) and never exceeded (stiff problem, fixed h = 1, cut).
//   6. Every method used (stage counts 5..124, 47 methods), forced with
//      -ts_rock4_min_stages = -ts_rock4_max_stages, on the nonlinear,
//      time-dependent problem (one copy per rank)
//        y1' = -y2 + y1 (1 - y1² - y2²) + 0.3 sin t
//        y2' =  y1 + y2 (1 - y1² - y2²)
//        y3' = -5 (y3 - cos t) - sin t + y1 y2,   y(0) = (0.8, 0.1, 2),
//      fixed steps h = 1/8 and 1/16 to t = 1, against y(1) from RK4 with
//      2·10^5 steps (agrees with 10^5 steps to 2e-14): the error falls by at
//      least 2^3.8 per halving.  y(1) at h = 1/16 is printed in full for
//      the cross-check against an independent transcription of
//      OrdinaryDiffEq.jl's ROCK4 step (QD S3 review).
//   7. Stability of every method used: one step h = 1 of y' = z y, y(0) = 1,
//      for z on [-0.98·(0.353 s² - 3), 0) (spacing 0.05 to -30, then 0.25;
//      the h·rho range for which the method is chosen), forced stage count,
//      spectral radius safety 1: one step, no cut, |y(1)| <= 1 + 1e-12, and
//      y(1) equal within 1e-12 to R(z) evaluated from the tables here.
//      The same tables check on the two methods the cap leaves out because
//      they are unstable (133 and 152 stages) finds |R| > 1.005 (negative
//      control); the 142-stage method between them is stable (printed).
//
// Exit code 0 when every check passes.

#include "mfem.hpp"

#ifndef MFEM_USE_PETSC
#include <iostream>
int main()
{
   std::cout << "test_ts_rock4: SKIPPED (MFEM built without PETSc)\n";
   return 0;
}
#else

#include "../../solver/ts_rock4.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_fail = 0;

void Check(bool ok, const std::string &what, int rank)
{
   if (!ok) { g_fail = 1; }
   if (rank == 0)
   {
      std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
   }
}

/// y(1) of the nonlinear problem of check 6 (RK4, 2·10^5 steps).
const double kNlRef[3] = {0.5239712127195414, 0.8747935203735767,
                          0.6390730531240915};

/// R(z) of method @a mi: one step h = 1 of y' = z y from y = 1, written from
/// the tables independently of TSStep_Rock4 (scalar, Abdulle's recurrence
/// and finishing procedure).
double TableR(int mi, double z)
{
   int start = 0;
   for (int i = 0; i < mi; i++) { start += 2 * kRock4Ms[i] - 1; }
   const double *r = kRock4Recf + start;
   double um2 = 1.0, um1 = 1.0 + r[0] * z * 1.0;
   double u = um1;
   for (int j = 2; j <= kRock4Ms[mi]; j++)
   {
      const double mu = r[2 * j - 3], ka = r[2 * j - 2];
      u = mu * z * um1 + (1.0 + ka) * um1 - ka * um2;
      um2 = um1;
      um1 = u;
   }
   const double *a = kRock4Fpa[mi], *b = kRock4Fpb[mi];
   const double k1 = z * u;
   const double k2 = z * (u + a[0] * k1);
   const double k3 = z * (u + a[1] * k1 + a[2] * k2);
   const double k4 = z * (u + a[3] * k1 + a[4] * k2 + a[5] * k3);
   return u + b[0] * k1 + b[1] * k2 + b[2] * k3 + b[3] * k4;
}

/// z values of check 7 for a method of @a s stages.
std::vector<double> StabilityGrid(int s)
{
   std::vector<double> z;
   const double zmin = -0.98 * (0.353 * s * s - 3.0);
   for (int k = 1; -0.05 * k >= std::max(zmin, -30.0); k++)
   {
      z.push_back(-0.05 * k);
   }
   for (double x = -30.25; x >= zmin; x -= 0.25) { z.push_back(x); }
   return z;
}

struct Problem
{
   int kind = 0;              // 0: lam (y - cos t) - sin t; 1: nonlinear;
                              // 2: lam y
   std::vector<double> lam;   // this rank's entries
   long long nf = 0;          // right-hand sides
   PetscReal last_t = -1.0;
   Vec last_x = nullptr;      // copy of the last right-hand side's state
};

PetscErrorCode RHS(TS, PetscReal t, Vec X, Vec F, void *ctx)
{
   PetscFunctionBeginUser;
   Problem *p = static_cast<Problem *>(ctx);
   const PetscScalar *x;
   PetscScalar *f;
   PetscCall(VecGetArrayRead(X, &x));
   PetscCall(VecGetArray(F, &f));
   if (p->kind == 1)
   {
      for (size_t b = 0; b + 2 < p->lam.size(); b += 3)
      {
         const double y1 = x[b], y2 = x[b + 1], y3 = x[b + 2];
         const double r2 = y1 * y1 + y2 * y2;
         f[b] = -y2 + y1 * (1.0 - r2) + 0.3 * std::sin(t);
         f[b + 1] = y1 + y2 * (1.0 - r2);
         f[b + 2] = -5.0 * (y3 - std::cos(t)) - std::sin(t) + y1 * y2;
      }
   }
   else
   {
      for (size_t i = 0; i < p->lam.size(); i++)
      {
         f[i] = p->kind == 2 ? p->lam[i] * x[i]
                : p->lam[i] * (x[i] - std::cos(t)) - std::sin(t);
      }
   }
   PetscCall(VecRestoreArray(F, &f));
   PetscCall(VecRestoreArrayRead(X, &x));
   p->nf++;
   p->last_t = t;
   if (p->last_x) { PetscCall(VecCopy(X, p->last_x)); }
   PetscFunctionReturn(PETSC_SUCCESS);
}

struct FsalCheck
{
   Problem *p;
   long long steps = 0, bad = 0;
};

PetscErrorCode Monitor(TS, PetscInt step, PetscReal t, Vec u, void *ctx)
{
   PetscFunctionBeginUser;
   FsalCheck *c = static_cast<FsalCheck *>(ctx);
   if (step == 0) { PetscFunctionReturn(PETSC_SUCCESS); }
   PetscBool same = PETSC_FALSE;
   PetscCall(VecEqual(u, c->p->last_x, &same));
   c->steps++;
   if (!same || c->p->last_t != t) { c->bad++; }
   PetscFunctionReturn(PETSC_SUCCESS);
}

struct Result
{
   double err = 0.0;
   long long nf = 0;
   Rock4Stats st;
   bool ok = true;
   long long fsal_steps = 0, fsal_bad = 0;
   std::vector<double> y;     // this rank's final state
};

/// Integrate to @a T; @a opts are set in the options database first and
/// cleared after.
Result Run(MPI_Comm comm, const std::vector<double> &lam, double T,
           double dt0, const std::string &opts, bool fsal, int kind = 0)
{
   Result r;
   PetscOptionsClear(nullptr);
   PetscOptionsInsertString(nullptr, opts.c_str());
   Problem p;
   p.kind = kind;
   p.lam = lam;
   const PetscInt n = static_cast<PetscInt>(lam.size());
   Vec y;
   VecCreateMPI(comm, n, PETSC_DETERMINE, &y);
   VecSet(y, kind == 2 ? 1.0 : 2.0);
   if (kind == 1)
   {
      PetscScalar *x;
      VecGetArray(y, &x);
      for (PetscInt b = 0; b + 2 < n; b += 3)
      {
         x[b] = 0.8; x[b + 1] = 0.1; x[b + 2] = 2.0;
      }
      VecRestoreArray(y, &x);
   }
   if (fsal) { VecDuplicate(y, &p.last_x); }
   TS ts;
   TSCreate(comm, &ts);
   TSSetProblemType(ts, TS_NONLINEAR);
   TSSetRHSFunction(ts, nullptr, RHS, &p);
   TSSetTime(ts, 0.0);
   TSSetTimeStep(ts, dt0);
   TSSetMaxTime(ts, T);
   TSSetMaxSteps(ts, 10000000);
   TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
   TSSetTolerances(ts, 1e-8, nullptr, 0.0, nullptr);
   FsalCheck fc{&p};
   if (fsal) { TSMonitorSet(ts, Monitor, &fc, nullptr); }
   PetscErrorCode ierr = TSSetFromOptions(ts);
   if (ierr == PETSC_SUCCESS) { ierr = TSSolve(ts, y); }
   if (ierr != PETSC_SUCCESS)
   {
      r.ok = false;
   }
   else
   {
      PetscReal tf;
      TSGetTime(ts, &tf);
      const PetscScalar *x;
      VecGetArrayRead(y, &x);
      double e = 0.0;
      r.y.assign(x, x + n);
      for (PetscInt i = 0; i < n; i++)
      {
         const double ex = kind == 1 ? kNlRef[i % 3]
                           : kind == 2 ? std::exp(lam[i] * tf)
                           : std::cos(tf) + std::exp(lam[i] * tf);
         e = std::max(e, std::abs(x[i] - ex));
      }
      VecRestoreArrayRead(y, &x);
      MPI_Allreduce(&e, &r.err, 1, MPI_DOUBLE, MPI_MAX, comm);
      if (TSIsRock4(ts)) { r.st = TSRock4GetStats(ts); }
      r.nf = p.nf;
      r.fsal_steps = fc.steps;
      r.fsal_bad = fc.bad;
   }
   TSDestroy(&ts);
   VecDestroy(&y);
   if (p.last_x) { VecDestroy(&p.last_x); }
   PetscOptionsClear(nullptr);
   return r;
}

}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   MFEMInitializePetsc(nullptr, nullptr, nullptr, nullptr);
   PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
   TSRock4Register();
   MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank(), P = Mpi::WorldSize();
   const long long N = 2000;
   const long long lo = N * rank / P, hi = N * (rank + 1) / P;
   if (rank == 0)
   {
      std::cout << "test_ts_rock4: " << P << " ranks, N " << N << std::endl;
   }

   // 1. Order
   {
      std::vector<double> lam(hi - lo, -20.0);
      double prev = 0.0, worst = 1e300;
      for (int k = 0; k < 4; k++)
      {
         const double h = 0.2 / (1 << k);
         Result r = Run(comm, lam, 1.0, h, "-ts_type rock4 -ts_adapt_type "
                        "none -ts_rock4_min_stages 20 -ts_rock4_max_stages 20",
                        false);
         if (rank == 0)
         {
            std::cout << "  order: h " << h << ", error " << r.err
                      << ", stages " << r.st.stages_max << ", steps "
                      << r.st.steps << std::endl;
         }
         if (k > 0) { worst = std::min(worst, std::log2(prev / r.err)); }
         prev = r.err;
         if (!r.ok) { worst = -1.0; }
      }
      Check(worst >= 3.6, "order: error falls by 2^" + std::to_string(worst) +
            " or more per halving (need 2^3.6)", rank);
   }

   // 2. + 3. Stiff, adaptive, with the first-same-as-last check
   std::vector<double> lam(hi - lo);
   for (long long i = lo; i < hi; i++)
   {
      lam[i - lo] = -std::pow(10.0, 4.0 * static_cast<double>(i) / (N - 1));
   }
   {
      Result r4 = Run(comm, lam, 10.0, 1e-6, "-ts_type rock4", true);
      Result r5 = Run(comm, lam, 10.0, 1e-6, "-ts_type rk -ts_rk_type 5dp",
                      false);
      if (rank == 0)
      {
         std::cout << "  stiff: rock4 error " << r4.err << ", "
                   << r4.st.steps << " steps (" << r4.st.rejected
                   << " rejected), stages " << static_cast<double>(
                      r4.st.stages_sum) / std::max(1LL, r4.st.steps)
                   << " mean / " << r4.st.stages_max << " max, f "
                   << r4.nf << " (" << r4.st.f_eig << " in "
                   << r4.st.eig_estimates << " spectral radius estimates, "
                   << r4.st.eig_unconverged << " unconverged), rho "
                   << r4.st.eig_last << "; rk 5dp error " << r5.err
                   << ", f " << r5.nf << std::endl;
      }
      Check(r4.ok && r4.err < 1e-6, "stiff: rock4 error below 1e-6", rank);
      Check(r4.st.eig_last >= 1e4 && r4.st.eig_last <= 1.5e4,
            "stiff: spectral radius estimate within [1, 1.5]e4", rank);
      Check(r4.st.stages_max > 5, "stiff: more than 5 stages in the largest "
            "step", rank);
      Check(r5.ok && r4.nf < r5.nf, "stiff: fewer right-hand sides than rk "
            "5dp", rank);
      Check(r4.fsal_steps == r4.st.steps && r4.fsal_bad == 0,
            "first same as last: the last right-hand side of each of " +
            std::to_string(r4.fsal_steps) + " accepted steps is at the "
            "accepted time and state, bitwise", rank);
   }

   // 4. Stage cap
   {
      Result r = Run(comm, lam, 10.0, 1.0, "-ts_type rock4 "
                     "-ts_rock4_max_stages 12", false);
      if (rank == 0)
      {
         std::cout << "  cap: error " << r.err << ", stages max "
                   << r.st.stages_max << ", h cut " << r.st.h_cut
                   << " times" << std::endl;
      }
      Check(r.ok && r.err < 1e-6 && r.st.stages_max <= 12 && r.st.h_cut > 0,
            "cap: at most 12 stages, h cut to the cap, error below 1e-6",
            rank);
   }

   // 5. Options outside their ranges; a cap between two methods
   {
      Result r = Run(comm, lam, 1.0, 1e-3, "-ts_type rock4 "
                     "-ts_rock4_min_stages 4", false);
      Check(!r.ok, "options: -ts_rock4_min_stages 4 stops with an error",
            rank);
      r = Run(comm, lam, 1.0, 1e-3, "-ts_type rock4 -ts_rock4_max_stages 125",
              false);
      Check(!r.ok, "options: -ts_rock4_max_stages 125 stops with an error",
            rank);
      r = Run(comm, lam, 1.0, 1.0, "-ts_type rock4 -ts_adapt_type none "
              "-ts_rock4_max_stages 31", false);
      if (rank == 0)
      {
         std::cout << "  cap 31: stages max " << r.st.stages_max
                   << ", h cut " << r.st.h_cut << " times" << std::endl;
      }
      Check(r.ok && r.st.stages_max == 30 && r.st.h_cut > 0,
            "options: -ts_rock4_max_stages 31 is rounded down to the "
            "30-stage method and never exceeded", rank);
   }

   // 6. Every method used, nonlinear and time-dependent, fixed steps
   {
      std::vector<double> lam3(3, 0.0);
      double worst = 1e300;
      int worst_s = 0, nmeth = 0;
      for (int mi = 0; mi < kRock4Methods; mi++)
      {
         const int s = kRock4Ms[mi] + 4;
         if (s > kRock4MaxStages) { break; }
         nmeth++;
         const std::string o = "-ts_type rock4 -ts_adapt_type none "
                               "-ts_rock4_min_stages " + std::to_string(s) +
                               " -ts_rock4_max_stages " + std::to_string(s);
         Result r1 = Run(comm, lam3, 1.0, 0.125, o, false, 1);
         Result r2 = Run(comm, lam3, 1.0, 0.0625, o, false, 1);
         const bool ok = r1.ok && r2.ok && r1.st.stages_max == s &&
                         r2.st.stages_max == s && r1.st.steps == 8 &&
                         r2.st.steps == 16;
         const double order = ok ? std::log2(r1.err / r2.err) : -1.0;
         if (order < worst) { worst = order; worst_s = s; }
         if (rank == 0)
         {
            char line[256];
            std::snprintf(line, sizeof line, "  nonlinear %3d stages: error "
                          "%.3e %.3e order %.2f  y(1) %.17g %.17g %.17g", s,
                          r1.err, r2.err, order, r2.y[0], r2.y[1], r2.y[2]);
            std::cout << line << std::endl;
         }
      }
      Check(nmeth == 47 && worst >= 3.8, "nonlinear: all " +
            std::to_string(nmeth) + " methods of 5 to 124 stages, error "
            "falls by 2^" + std::to_string(worst) + " or more per halving ("
            + std::to_string(worst_s) + " stages; need 2^3.8)", rank);
   }

   // 7. Stability of every method used; the two left out fail it
   {
      double worst_r = 0.0, worst_d = 0.0;
      int worst_s = 0, bad_runs = 0;
      for (int mi = 0; mi < kRock4Methods; mi++)
      {
         const int s = kRock4Ms[mi] + 4;
         const std::vector<double> z = StabilityGrid(s);
         double tmax = 0.0;
         for (double zz : z) { tmax = std::max(tmax, std::abs(TableR(mi, zz))); }
         if (s > kRock4MaxStages)
         {
            if (rank == 0)
            {
               std::cout << "  stability " << s << " stages (left out): "
                         << "tables max|R| " << tmax << std::endl;
            }
            if (s == 133 || s == 152)
            {
               Check(tmax > 1.005, "stability: the " + std::to_string(s) +
                     "-stage method the cap leaves out has |R| > 1.005 "
                     "(negative control)", rank);
            }
            continue;
         }
         const long long nz = static_cast<long long>(z.size());
         const long long zlo = nz * rank / P, zhi = nz * (rank + 1) / P;
         std::vector<double> lamz(z.begin() + zlo, z.begin() + zhi);
         Result r = Run(comm, lamz, 1.0, 1.0, "-ts_type rock4 "
                        "-ts_adapt_type none -ts_rock4_eig_safety 1 "
                        "-ts_rock4_min_stages " + std::to_string(s) +
                        " -ts_rock4_max_stages " + std::to_string(s),
                        false, 2);
         double rmax = 0.0, dmax = 0.0;
         for (size_t i = 0; r.ok && i < lamz.size(); i++)
         {
            rmax = std::max(rmax, std::abs(r.y[i]));
            dmax = std::max(dmax, std::abs(r.y[i] - TableR(mi, lamz[i])));
         }
         double loc[2] = {rmax, dmax}, glob[2];
         MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_MAX, comm);
         if (!r.ok || r.st.steps != 1 || r.st.h_cut != 0 ||
             r.st.stages_max != s) { bad_runs++; }
         if (glob[0] > worst_r) { worst_r = glob[0]; worst_s = s; }
         worst_d = std::max(worst_d, glob[1]);
         worst_r = std::max(worst_r, tmax);
      }
      if (rank == 0)
      {
         std::cout << "  stability: 47 methods, max|y(1)| " << worst_r
                   << " (" << worst_s << " stages), largest difference to "
                   << "the tables' R " << worst_d << std::endl;
      }
      Check(bad_runs == 0 && worst_r <= 1.0 + 1e-12 && worst_d <= 1e-12,
            "stability: every method of 5 to 124 stages takes its one step "
            "uncut with |R(z)| <= 1 on its whole interval, and the step "
            "equals the tables' R(z)", rank);
   }

   int fail_all = 0;
   MPI_Allreduce(&g_fail, &fail_all, 1, MPI_INT, MPI_MAX, comm);
   if (rank == 0)
   {
      std::cout << (fail_all ? "FAILED" : "ALL PASSED") << std::endl;
   }
   MFEMFinalizePetsc();
   return fail_all ? 1 : 0;
}

#endif // MFEM_USE_PETSC
