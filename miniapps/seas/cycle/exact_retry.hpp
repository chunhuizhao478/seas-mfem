// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// exact_retry.hpp — [cycle.qd] exact_retry: undo a rejected QD step by
// copying back the state the step started from, instead of PETSc's own
// rollback (fault opening version B, the UCERF quench run's QD crash of
// 2026-11-07).
//
// Why.  PETSc's RK undoes a rejected step by subtracting it again
// (TSRollBack_RK, PETSc 3.19 and 3.24.6):
//
//    U_n  =  U_{n+1} - h sum_j b_j K_j       (all s stages, b_s = 0 for 5dp)
//
// The update U_{n+1} = U_n + h sum_j b_j K_j multiplies every stage rate by
// its weight, the zero weight too, so a last-stage rate K_s = Inf puts
// 0 * Inf = NaN into U_{n+1}; the rollback multiplies it by 0 again, and the
// NaN stays in U_n.  The step's own error estimate ignores NaN entries
// (VecErrorWeightedNorms, infinity norm: PetscMax(nrm, NaN) = nrm), so it
// still rejects the step on its other entries, but every retry starts from
// the poisoned state, fails the stage check (D-B5) at any dt, and the solve
// ends DIVERGED_STEP_REJECTED.  Vista 1128426 (d9dcd2b5): step 82624 at
// 12,222 yr, the sjf1 point at s = 26 km, a trial stage with sigma_n'
// -2,610 MPa, the law's sigma_n 0, |V| = tau / eta = 87,644 m/s, psi -61 and
// dpsi/dt = Inf; error 1.87e17, rejected; then ten retries from psi = NaN
// down to dt = 3.3 s.
//
// What.  TSPreStep copies the state the step starts from; TSPreStage, at
// the first stage of a retry (PETSc's rejection count rose since the
// copy), copies it back over PETSc's rollback before the stage reads it.
// Then the retry starts from U_n exactly (also free of the round-off of
// subtracting a huge rejected update), while its first rate stays the
// previous step's last one (FSAL), which no rejection touches.  A step
// without a rejection is untouched: no copy back, every bit as before.
//
// The PETSc hooks live in QDPhase (the TS context); this class holds the
// copy and the counts.  Unit test: tests/2d_inplane/test_exact_retry.cpp.

#ifndef MFEM_SEAS_CYCLE_EXACT_RETRY_HPP
#define MFEM_SEAS_CYCLE_EXACT_RETRY_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_PETSC
#include "petsc.h"

#include <cmath>

namespace mfem
{
namespace seas
{
namespace cycle
{

class ExactRetry
{
public:
   ExactRetry() = default;
   ExactRetry(const ExactRetry &) = delete;
   ExactRetry &operator=(const ExactRetry &) = delete;
   ~ExactRetry() { if (saved_) { VecDestroy(&saved_); } }

   /// TSPreStep: copy the state this step starts from.
   PetscErrorCode Save(TS ts)
   {
      PetscFunctionBeginUser;
      Vec U;
      PetscCall(TSGetSolution(ts, &U));
      if (saved_)
      {
         PetscInt n_saved = 0, n_u = 0;
         PetscCall(VecGetLocalSize(saved_, &n_saved));
         PetscCall(VecGetLocalSize(U, &n_u));
         if (n_saved != n_u) { PetscCall(VecDestroy(&saved_)); }
      }
      if (!saved_) { PetscCall(VecDuplicate(U, &saved_)); }
      PetscCall(VecCopy(U, saved_));
      PetscCall(TSGetStepRejections(ts, &rej_seen_));
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// TSPreStage: at the first stage of a retry, copy the saved state back
   /// over PETSc's rollback.  @a restored says whether it did; @a bad (if
   /// not null) gets the number of non-finite entries the rollback had left
   /// on this rank.  Collective only through PETSc's own calls.
   PetscErrorCode RestoreIfRetry(TS ts, bool &restored, long *bad = nullptr)
   {
      PetscFunctionBeginUser;
      restored = false;
      if (bad) { *bad = 0; }
      PetscInt rej = 0;
      PetscCall(TSGetStepRejections(ts, &rej));
      if (!saved_ || rej <= rej_seen_) { PetscFunctionReturn(PETSC_SUCCESS); }
      Vec U;
      PetscCall(TSGetSolution(ts, &U));
      if (bad)
      {
         PetscInt n = 0;
         const PetscScalar *a = nullptr;
         PetscCall(VecGetLocalSize(U, &n));
         PetscCall(VecGetArrayRead(U, &a));
         for (PetscInt k = 0; k < n; k++)
         {
            if (!std::isfinite(static_cast<double>(PetscRealPart(a[k])))) { (*bad)++; }
         }
         PetscCall(VecRestoreArrayRead(U, &a));
      }
      PetscCall(VecCopy(saved_, U));
      rej_seen_ = rej;
      restores_++;
      restored = true;
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// Retries that started from the copy (this rank's count; every rank
   /// sees the same rejections, so the counts agree).
   long Restores() const { return restores_; }

private:
   Vec saved_ = nullptr;
   PetscInt rej_seen_ = 0;
   long restores_ = 0;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_USE_PETSC
#endif // MFEM_SEAS_CYCLE_EXACT_RETRY_HPP
