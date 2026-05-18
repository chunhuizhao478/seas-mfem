// R-003 (REVIEW.md 2026-05-16) regression test.
//
// `bp5_verification_full.cpp:1541` previously did:
//
//     Vector state(fault_op.StateSize());
//
// On a 400-rank Frontera BP5 run, many ranks own zero fault DOFs (the
// strike-slip fault is a thin 2D slab inside a 3D box).  For those ranks,
// `fault_op.StateSize()` is 0 and `Vector(0)` does NOT allocate
// (`vector.hpp:574-582`: `if (s > 0) data.New(s);`).  So
// `state.GetMemory().h_ptr == NULL` and `state.GetMemory().Empty() == true`.
//
// Later, `petsc_ode->Run(state, ...)` calls
// `PetscParVector::PlaceMemory(state.GetMemory(), true)`.  That builds
// `pdata.MakeAliasForSync(mem, 0, 0, ...)` which sets
// `pdata.h_ptr = mem.h_ptr + 0 = NULL`.  At the end of `TSSolve`,
// `PetscParVector::ResetMemory()` (linalg/petsc.cpp:899) checks
// `MFEM_VERIFY(!pdata.Empty(),...)` — which aborts on the NULL pointer.
// User-visible symptom: "PetscParVector::ResetMemory: Vector data is
// empty" on every zero-fault-DOF rank simultaneously.
//
// Driver-side fix: pad to ≥1 element, then shrink back.  SetSize only
// re-allocates when `new_size > capacity` (`vector.hpp:584-602`), so
// the 1-element allocation from `SetSize(1)` survives the immediate
// `SetSize(0)` shrink.  Resulting `state` has size=0 but capacity=1
// and h_ptr != NULL → `GetMemory().Empty() == false`.
//
// This test exercises that invariant without spinning up a full
// PETSc + SEAS operator stack.  The MFEM-side ResetMemory fix (allow
// zero-length aliases) is a separate upstream patch; this is the
// driver-side hot-fix.

#include "mfem.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace mfem;

static int num_tests = 0, num_passed = 0, num_failed = 0;
#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
   std::cout << "=== test_bp5_petsc_ts_zero_fault_rank (R-003) ===\n";

   // -- Pre-condition: the naive constructor pattern produces an empty
   //    Memory when StateSize() is 0.  This is what triggered the
   //    PetscParVector::ResetMemory abort.
   {
      const int actual = 0;     // simulates rank with no owned fault DOFs
      Vector naive(actual);
      TEST_ASSERT(naive.Size() == 0,
                  "R-003 pre-condition: naive Vector(0) has size 0");
      TEST_ASSERT(naive.GetMemory().Empty(),
                  "R-003 pre-condition: naive Vector(0) has empty "
                  "Memory (h_ptr == NULL) — this is what tripped "
                  "PetscParVector::ResetMemory");
   }

   // -- The fix: pad to >= 1, then shrink back.  Post-fix, the Memory
   //    must NOT be empty even when the logical size is 0.
   {
      const int actual = 0;
      Vector state;
      const int padded = std::max(actual, 1);
      state.SetSize(padded);
      state.SetSize(actual);    // shrink back; allocation preserved

      TEST_ASSERT(state.Size() == actual,
                  "R-003: padded-shrink yields the actual size "
                  "expected by SetInitialCondition's verify");
      TEST_ASSERT(!state.GetMemory().Empty(),
                  "R-003: padded-shrink Memory is NOT empty (h_ptr "
                  "is allocated from SetSize(1) and survives the "
                  "SetSize(0) shrink because SetSize only re-allocates "
                  "when new_size > capacity)");
      TEST_ASSERT(state.GetMemory().Capacity() >= 1,
                  "R-003: padded-shrink Memory retains capacity >= 1 "
                  "(the SetSize(1) allocation is preserved by the "
                  "SetSize(0) shrink)");
   }

   // -- Non-zero actual: the trick must be transparent.  When
   //    `actual > 0`, the padded-shrink reduces to a single SetSize
   //    call and behaves identically to the naive Vector(actual).
   {
      const int actual = 7;     // simulates rank with owned fault DOFs
      Vector state;
      const int padded = std::max(actual, 1);
      state.SetSize(padded);
      state.SetSize(actual);

      TEST_ASSERT(state.Size() == actual,
                  "R-003: non-zero case keeps the requested size");
      TEST_ASSERT(!state.GetMemory().Empty(),
                  "R-003: non-zero case has a non-empty Memory (this "
                  "was always true for the naive constructor too; we "
                  "just verify the workaround doesn't regress it)");
   }

   // -- R-204 (REVIEW.md round 3): static-grep check that the BP5
   //    driver still applies the padded-shrink workaround.  The
   //    sub-tests above pin the MFEM-side SetSize invariant the fix
   //    depends on, but say nothing about whether the driver itself
   //    still uses the trick.  Without this check, a refactor that
   //    reverts the driver to `Vector state(fault_op.StateSize())`
   //    would crash on Frontera but still pass this unit test.
   //
   //    The check reads the driver source and asserts the canonical
   //    workaround pattern is present.  It is conservative: it just
   //    requires the literal `padded = std::max(actual, 1)` and
   //    `state.SetSize(actual)` substrings.  A future stylistic
   //    refactor (e.g., renaming `padded` → `padded_size`) would
   //    require updating the search strings here too; that's
   //    intentional — any non-trivial change to the workaround
   //    deserves a second look at this gate.
   //
   //    Path resolution: the seas Makefile invokes the test binary
   //    from `miniapps/seas/`, so the relative path below is
   //    consistent with `make test-bp5-petsc-ts-zero-fault-rank`.
   //    If the file is not found (test invoked from a different
   //    CWD), the assertion is skipped with an INFO message rather
   //    than failing — keeps the unit test runnable in arbitrary
   //    contexts.
   {
      const std::string driver_path =
         "tests/verification/bp5_verification_full.cpp";
      std::ifstream driver(driver_path);
      if (driver.is_open())
      {
         std::stringstream buf;
         buf << driver.rdbuf();
         const std::string src = buf.str();
         const bool has_padded =
            src.find("padded = std::max(actual, 1)") != std::string::npos;
         const bool has_shrink =
            src.find("state.SetSize(actual)") != std::string::npos;
         TEST_ASSERT(has_padded && has_shrink,
                     "R-204: BP5 driver must still contain the R-003 "
                     "padded-shrink workaround at " << driver_path
                     << "; a refactor that removed it would crash on "
                     "Frontera on the next zero-fault-DOF rank.  "
                     "(has_padded=" << has_padded
                     << " has_shrink=" << has_shrink << ")");
      }
      else
      {
         std::cout << "  INFO: R-204 driver-source check skipped — "
                      "could not open " << driver_path
                   << " (test run from unexpected CWD?  Expected "
                      "to run from miniapps/seas/ — `make "
                      "test-bp5-petsc-ts-zero-fault-rank`)\n";
      }
   }

   // -- The full PetscParVector::PlaceMemory + ResetMemory roundtrip
   //    on a real Vec is MPI-only and requires PETSc; it is covered by
   //    a separate MPI integration test (test_bp5_petsc_ts_zero_fault_rank_mpi.cpp,
   //    not yet authored — pending Frontera reproduction of the failing
   //    sbatch with the patched binary).  This serial unit test
   //    guarantees the invariant that the workaround depends on,
   //    which is sufficient to catch regression of the SetSize-shrink
   //    semantics in MFEM.

   std::cout << "\n=== Summary: " << num_passed << " / " << num_tests
             << " passed; " << num_failed << " failed ===\n";
   return num_failed > 0 ? 1 : 0;
}
