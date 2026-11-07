// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_exact_retry.cpp — [cycle.qd] exact_retry
// (cycle/exact_retry.hpp): a PETSc RK 5dp toy with the same TS setup as
// QDPhase (infinity error norm, atol 1e-7, rtol 1e-50, a stage-state check
// as D-B5's StageStateFinite_), two state entries per rank:
//
//   u0' = -u0 (decays), u1' = 1e-3 (drifts).
//
// A test hook in the toy's own right-hand side puts a wild but finite rate
// into stage 6 of step 1 (u0' = 1e6, the error estimate rejects the step)
// and, optionally, an Inf into stage 7 (u1' = Inf, the overflow of Vista
// 1128426 at 12,222 yr):
//
//   R1  clean, hooks on vs off: every printed step and the final state bit
//       for bit, no retry (a step without a rejection is untouched).
//   R2  PETSc's own rollback with the Inf: the retry starts from NaN, the
//       stage check rejects every retry, and TSSolve ends with
//       DIVERGED_STEP_REJECTED (the crash pattern of Vista 1127480).
//   R3  hooks on with the Inf: one restore, the rollback had left NaN, the
//       run reaches t = 1 with a finite state.
//   R4  hooks on, the same rejected step without the Inf: bit for bit the
//       same steps and final state as R3 (after the restore both retries
//       start from U_n exactly, with the same dt and the same FSAL rate).
//
// Run on 1 and 2 ranks (make -f tests/2d_inplane/2d_inplane.mk
// test-exact-retry).

#include "mfem.hpp"
#include "../../cycle/exact_retry.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

using namespace mfem;
using mfem::seas::cycle::ExactRetry;

namespace
{

struct Toy
{
   int evals = 0;
   int big_at = -1;   ///< RHS call (counted on every rank alike) with u0' = 1e6
   int inf_at = -1;   ///< RHS call with u1' = Inf (rank 0 only)
   bool hooks = false;
   ExactRetry retry;
   long nonfinite_retries = 0;
   std::vector<double> rows;   ///< t, u0, u1 of rank 0 after every step
   int rank = 0;
};

PetscErrorCode RHS(TS, PetscReal, Vec U, Vec F, void *ctx)
{
   PetscFunctionBeginUser;
   Toy *toy = static_cast<Toy *>(ctx);
   toy->evals++;
   const PetscScalar *u;
   PetscScalar *f;
   PetscInt n;
   PetscCall(VecGetLocalSize(U, &n));
   PetscCall(VecGetArrayRead(U, &u));
   PetscCall(VecGetArray(F, &f));
   for (PetscInt i = 0; i + 1 < n; i += 2)
   {
      f[i] = -u[i];
      f[i + 1] = 1.0e-3;
   }
   if (toy->evals == toy->big_at) { f[0] = 1.0e6; }
   if (toy->evals == toy->inf_at && toy->rank == 0) { f[1] = INFINITY; }
   PetscCall(VecRestoreArrayRead(U, &u));
   PetscCall(VecRestoreArray(F, &f));
   PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode StageFinite(TS ts, PetscReal, Vec Y, PetscBool *accept)
{
   PetscFunctionBeginUser;
   const PetscScalar *a;
   PetscInt n;
   PetscCall(VecGetLocalSize(Y, &n));
   PetscCall(VecGetArrayRead(Y, &a));
   int ok = 1;
   for (PetscInt i = 0; i < n; i++)
   {
      if (!std::isfinite(static_cast<double>(PetscRealPart(a[i])))) { ok = 0; }
   }
   PetscCall(VecRestoreArrayRead(Y, &a));
   MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(ts));
   PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, comm));
   *accept = ok ? PETSC_TRUE : PETSC_FALSE;
   PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PreStep(TS ts)
{
   PetscFunctionBeginUser;
   void *ctx;
   PetscCall(TSGetApplicationContext(ts, &ctx));
   PetscCall(static_cast<Toy *>(ctx)->retry.Save(ts));
   PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PreStage(TS ts, PetscReal)
{
   PetscFunctionBeginUser;
   void *ctx;
   PetscCall(TSGetApplicationContext(ts, &ctx));
   Toy *toy = static_cast<Toy *>(ctx);
   bool restored = false;
   long bad = 0;
   PetscCall(toy->retry.RestoreIfRetry(ts, restored, &bad));
   if (restored)
   {
      long all = 0;
      PetscCallMPI(MPI_Allreduce(&bad, &all, 1, MPI_LONG, MPI_SUM,
                                 PetscObjectComm(reinterpret_cast<PetscObject>(ts))));
      if (all > 0) { toy->nonfinite_retries++; }
   }
   PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PostStep(TS ts)
{
   PetscFunctionBeginUser;
   void *ctx;
   PetscCall(TSGetApplicationContext(ts, &ctx));
   Toy *toy = static_cast<Toy *>(ctx);
   PetscReal t;
   Vec U;
   const PetscScalar *u;
   PetscCall(TSGetTime(ts, &t));
   PetscCall(TSGetSolution(ts, &U));
   PetscCall(VecGetArrayRead(U, &u));
   toy->rows.push_back(static_cast<double>(t));
   toy->rows.push_back(static_cast<double>(PetscRealPart(u[0])));
   toy->rows.push_back(static_cast<double>(PetscRealPart(u[1])));
   PetscCall(VecRestoreArrayRead(U, &u));
   PetscFunctionReturn(PETSC_SUCCESS);
}

struct Result
{
   TSConvergedReason reason = TS_CONVERGED_ITERATING;
   PetscInt steps = 0, rejections = 0;
   std::vector<double> rows, final_state;
   long restores = 0, nonfinite_retries = 0;
};

Result Run(int big_at, int inf_at, bool hooks)
{
   Result r;
   Toy toy;
   toy.big_at = big_at;
   toy.inf_at = inf_at;
   toy.hooks = hooks;
   MPI_Comm_rank(MPI_COMM_WORLD, &toy.rank);
   Vec U;
   TS ts;
   MFEM_VERIFY(VecCreateMPI(MPI_COMM_WORLD, 2, PETSC_DETERMINE, &U) == PETSC_SUCCESS, "Vec");
   MFEM_VERIFY(VecSet(U, 1.0) == PETSC_SUCCESS, "VecSet");
   MFEM_VERIFY(TSCreate(MPI_COMM_WORLD, &ts) == PETSC_SUCCESS, "TSCreate");
   MFEM_VERIFY(TSSetType(ts, TSRK) == PETSC_SUCCESS && TSRKSetType(ts, TSRK5DP) == PETSC_SUCCESS,
               "TSRK 5dp");
   MFEM_VERIFY(TSSetRHSFunction(ts, nullptr, RHS, &toy) == PETSC_SUCCESS, "RHS");
   MFEM_VERIFY(TSSetSolution(ts, U) == PETSC_SUCCESS, "solution");
   MFEM_VERIFY(TSSetTimeStep(ts, 0.1) == PETSC_SUCCESS, "dt");
   MFEM_VERIFY(TSSetMaxTime(ts, 1.0) == PETSC_SUCCESS, "max time");
   MFEM_VERIFY(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP) == PETSC_SUCCESS, "final");
   MFEM_VERIFY(TSSetTolerances(ts, 1e-7, nullptr, 1e-50, nullptr) == PETSC_SUCCESS, "tol");
   TSAdapt adapt;
   MFEM_VERIFY(TSGetAdapt(ts, &adapt) == PETSC_SUCCESS, "adapt");
   MFEM_VERIFY(TSAdaptSetType(adapt, TSADAPTBASIC) == PETSC_SUCCESS, "basic");
   // -ts_adapt_wnormtype infinity, as ucerf2d/petsc_ts_rk45_ucerf.cfg
   MFEM_VERIFY(PetscOptionsSetValue(nullptr, "-ts_adapt_wnormtype", "infinity")
               == PETSC_SUCCESS, "options");
   MFEM_VERIFY(TSSetFromOptions(ts) == PETSC_SUCCESS, "TSSetFromOptions");
   MFEM_VERIFY(TSSetFunctionDomainError(ts, StageFinite) == PETSC_SUCCESS, "domain");
   MFEM_VERIFY(TSSetApplicationContext(ts, &toy) == PETSC_SUCCESS, "ctx");
   MFEM_VERIFY(TSSetPostStep(ts, PostStep) == PETSC_SUCCESS, "post");
   if (hooks)
   {
      MFEM_VERIFY(TSSetPreStep(ts, PreStep) == PETSC_SUCCESS, "pre step");
      MFEM_VERIFY(TSSetPreStage(ts, PreStage) == PETSC_SUCCESS, "pre stage");
   }
   MFEM_VERIFY(TSSetUp(ts) == PETSC_SUCCESS, "setup");
   PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
   (void) TSSolve(ts, U);   // DIVERGED_STEP_REJECTED returns an error code
   PetscPopErrorHandler();
   TSGetConvergedReason(ts, &r.reason);
   TSGetStepNumber(ts, &r.steps);
   TSGetStepRejections(ts, &r.rejections);
   const PetscScalar *u;
   PetscInt n;
   VecGetLocalSize(U, &n);
   VecGetArrayRead(U, &u);
   for (PetscInt i = 0; i < n; i++) { r.final_state.push_back(PetscRealPart(u[i])); }
   VecRestoreArrayRead(U, &u);
   r.rows = toy.rows;
   r.restores = toy.retry.Restores();
   r.nonfinite_retries = toy.nonfinite_retries;
   TSDestroy(&ts);
   VecDestroy(&U);
   return r;
}

int fails = 0;
int rank = 0;

void Check(bool ok, const char *what)
{
   int all = ok ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   if (rank == 0) { std::cout << (all ? "  PASSED: " : "  FAILED: ") << what << std::endl; }
   if (!all) { fails++; }
}

bool AllFinite(const std::vector<double> &v)
{
   for (double x : v) { if (!std::isfinite(x)) { return false; } }
   return true;
}

bool SameBits(const Result &a, const Result &b)
{
   return a.reason == b.reason && a.steps == b.steps && a.rejections == b.rejections &&
          a.rows == b.rows && a.final_state == b.final_state;
}

} // namespace

int main(int argc, char **argv)
{
   Mpi::Init(argc, argv);
   MFEMInitializePetsc(&argc, &argv, nullptr, nullptr);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int size = 1;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   if (rank == 0) { std::cout << "test_exact_retry on " << size << " rank(s)" << std::endl; }

   // The first step evaluates 7 stages (calls 1-7), every later one 6
   // (FSAL): step 1 is calls 8-13, its stage 6 is call 12, stage 7 call 13.
   const int big = 12, inf = 13;

   const Result clean_off = Run(-1, -1, false);
   const Result clean_on = Run(-1, -1, true);
   Check(clean_off.reason == TS_CONVERGED_TIME && clean_off.rejections == 0,
         "R1 the clean toy reaches t = 1 without a rejection");
   Check(SameBits(clean_off, clean_on) && clean_on.restores == 0,
         "R1 hooks on = off bit for bit (every step, the final state), no restore");

   const Result petsc_inf = Run(big, inf, false);
   if (rank == 0)
   {
      std::cout << "    R2: reason " << TSConvergedReasons[petsc_inf.reason] << ", steps "
                << petsc_inf.steps << ", rejections " << petsc_inf.rejections
                << ", final u1 (rank 0) " << petsc_inf.final_state[1] << std::endl;
   }
   Check(petsc_inf.reason == TS_DIVERGED_STEP_REJECTED,
         "R2 PETSc's rollback with the Inf: DIVERGED_STEP_REJECTED (the crash pattern)");
   Check(rank != 0 || !std::isfinite(petsc_inf.final_state[1]),
         "R2 the state the retries started from holds NaN (rank 0, u1)");

   const Result retry_inf = Run(big, inf, true);
   if (rank == 0)
   {
      std::cout << "    R3: reason " << TSConvergedReasons[retry_inf.reason] << ", steps "
                << retry_inf.steps << ", rejections " << retry_inf.rejections
                << ", restores " << retry_inf.restores << " (non-finite rollbacks "
                << retry_inf.nonfinite_retries << ")" << std::endl;
   }
   Check(retry_inf.reason == TS_CONVERGED_TIME && AllFinite(retry_inf.final_state) &&
         AllFinite(retry_inf.rows),
         "R3 exact retry with the Inf: reaches t = 1, every state finite");
   Check(retry_inf.rejections == 1 && retry_inf.restores == 1 &&
         retry_inf.nonfinite_retries == 1,
         "R3 one rejection, one restore, and its rollback had left a non-finite entry");

   const Result retry_big = Run(big, -1, true);
   Check(SameBits(retry_inf, retry_big),
         "R4 exact retry: the same rejected step with or without the Inf, bit for bit");
   const Result petsc_big = Run(big, -1, false);
   double dmax = 0.0;
   for (std::size_t i = 0; i < petsc_big.final_state.size() &&
        i < retry_big.final_state.size(); i++)
   {
      dmax = std::max(dmax, std::fabs(petsc_big.final_state[i] - retry_big.final_state[i]));
   }
   if (rank == 0)
   {
      std::cout << "    note: without the Inf, PETSc's rollback vs the copy: final state "
                << "differs by " << dmax << " (the round-off of subtracting the rejected "
                << "update)" << std::endl;
   }
   Check(petsc_big.reason == TS_CONVERGED_TIME,
         "R4 without the Inf, PETSc's own rollback also reaches t = 1");

   if (rank == 0)
   {
      std::cout << (fails ? "test_exact_retry: FAILED" : "test_exact_retry: ALL PASS")
                << std::endl;
   }
   MFEMFinalizePetsc();
   return fails ? 1 : 0;
}
