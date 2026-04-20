# Code Review: TPV102 v9.0.0 Pelties-9 per-side flux fix — 2026-04-20 (round 2, final pass)

## Review Scope

- Plan: `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v9.0.0_seissol_flux_comparison.md`
- Files re-reviewed (fresh adversarial pass after round-1 fixes were applied):
  - `miniapps/seas/dynamic/wave_operator.inl` (interior-fault §14.2; shared-fault §14.3; including round-1 R-002/R-003 guards)
  - `miniapps/seas/dynamic/fault_face_flux.cpp` (R-007/R-008 bimaterial guard with 1e-12 relative tolerance)
  - `miniapps/seas/dynamic/godunov_flux.cpp` + `.hpp`
  - `miniapps/seas/fault/fault_basis.hpp` (sign-flipping convention verification)
  - `miniapps/seas/tests/unit/test_fault_face_flux_bimaterial_guard.cpp` (NEW for R-007)
  - `miniapps/seas/tests/unit/test_fault_face_flux_frame_and_flux.cpp` (§3.1b)
  - `miniapps/seas/tests/unit/test_fault_flux_interior_vs_shared_branch_equivalence.cpp` (§3.1d)
  - `miniapps/seas/tests/unit/test_fault_face_flux_per_side_assembly.cpp` (§3.1f; round-1 R-005 V2>0 guard applied)
  - `miniapps/seas/tests/unit/test_godunov_identity_normal_reversal.cpp` (§3.1c; round-1 R-006 comment updated)
  - `miniapps/seas/tests/unit/test_godunov_interior_equal_sides_identity.cpp` (§3.1e)
  - `miniapps/seas/Makefile` (test-target registration; round-1 R-001 §3.1 tests added)
- Round-1 findings confirmed applied: R-001 (test: target), R-002 (outer-conditional guard at `wave_operator.inl:710`), R-003 (shared-fault routing guard at `:1179`), R-005 (V2>0 fixture guard), R-006 (stale comment in §3.1c), R-007 (new bimaterial-guard test), R-008 (1e-12 relative tolerance in guard).  R-004 only partially applied (header comment updated but inline Case B label at `:269-270` still references the obsolete "Godunov identity (§3.1c)" premise — see R-203 below).
- Domain context: `miniapps/seas/CLAUDE.md`, `seas-mfem/CLAUDE.md`, plan §14–§18 self-audit, FaultBasis sign convention (`fault_basis.hpp:25-45, :342-451`).

## Findings (round 2 — NEW, not duplicates of round 1)

### [R-201] CRITICAL [test_fault_face_flux_bimaterial_guard.cpp:60-62 + Makefile:1003-1009] — Death-test exit-code logic is broken; test spuriously passes when the guard is removed

**Category:** BUG

**Description:**
The R-007 bimaterial-guard test is a "death test": the child process is expected to abort inside `MFEM_VERIFY`, and the shell target accepts **any non-zero exit code** as PASS.  But the child's *fall-through* path (when `Evaluate` returns WITHOUT aborting — i.e., when the guard has been silently removed by a future edit) also exits with non-zero code because it does `return 1`.  The shell target `[ $rc -ne 0 ]` cannot distinguish "child aborted inside the guard (rc=134 on SIGABRT)" from "child fell through and returned 1 because the guard was gone", and in both cases declares `PASS`.

Consequence: the test's stated purpose — catching a regression of the bimaterial guard — is not fulfilled.  The test exists to detect someone removing `MFEM_VERIFY(homog_ok(...))` from `fault_face_flux.cpp:90-98`.  If that removal happens, the next `make test` run will emit "PASS: bimaterial guard aborted as expected (rc=1)" despite the guard being gone.

The `"FAIL: Evaluate returned ..."` message the child writes to stderr (line 60-61) is redirected to `/dev/null` by the shell target (`2>/dev/null` on Makefile line 1003), so it cannot be used as a signal either.

**Trigger:**
Any future edit that removes, weakens, or comments out the `MFEM_VERIFY(homog_ok(data.Zp_plus, data.Zp_minus) && homog_ok(data.Zs_plus, data.Zs_minus), ...)` guard at `fault_face_flux.cpp:90-98`.  Run `make test-fault-face-flux-bimaterial-guard`.

**Actual behavior:**
```
$ make test-fault-face-flux-bimaterial-guard
PASS: bimaterial guard aborted as expected (rc=1)
```
— emitted both when the guard IS firing (rc=134) and when the guard is GONE and the child returned 1 (rc=1).

**Expected behavior:**
PASS only when the child aborted (rc ≠ 0 AND rc ≠ 1 from return path), or equivalently: the child returns a distinct exit code for "guard gone" that the shell treats as FAIL.

**Suggested fix:**
Swap the child's fall-through return code so `rc==0` (not `rc==1`) signals "guard did not abort":
```diff
diff --git a/miniapps/seas/tests/unit/test_fault_face_flux_bimaterial_guard.cpp b/miniapps/seas/tests/unit/test_fault_face_flux_bimaterial_guard.cpp
@@
    // Expect MFEM_VERIFY to abort here.  If the call returns at all,
    // the bimaterial guard has been silently removed or weakened —
    // that is the regression this test catches.
    ff.Evaluate(d, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

-   std::cerr << "FAIL: Evaluate returned despite bimaterial DOFData; "
-             << "the R-F08 / R-008 guard is no longer firing.\n";
-   return 1;
+   // REGRESSION PATH: if we reach this line, the guard is gone.  Return
+   // 0 so the shell's `[ $rc -ne 0 ]` check flips to FAIL — a non-zero
+   // exit from abort() is the PASS signal; a zero exit here is the
+   // UNIQUE marker that Evaluate returned normally.
+   std::cerr << "REGRESSION: Evaluate returned despite bimaterial "
+             << "DOFData; the R-F08 / R-008 guard is no longer firing.\n";
+   return 0;
 }
```
The shell target on `Makefile:1003-1009` is already correct for this semantics (no change needed there).

**Test case:**
```bash
# Meta-test: prove the fix correctly distinguishes the two cases.
# Run from repo root.
set -e

# Case 1: guard present (current state).  Expect test to PASS.
cd miniapps/seas
make test-fault-face-flux-bimaterial-guard 2>&1 | tee /tmp/guard_case1.log
grep -q "^PASS: " /tmp/guard_case1.log \
   || { echo "FATAL: guard-present case not PASSing"; exit 1; }

# Case 2: simulate guard removal and confirm test FAILs.
cp dynamic/fault_face_flux.cpp /tmp/fault_face_flux.cpp.bak
sed -i.bak 's/MFEM_VERIFY(homog_ok/MFEM_VERIFY(true || homog_ok/' \
   dynamic/fault_face_flux.cpp
rm -f dynamic/fault_face_flux.o
make -j8 seas_test_fault_face_flux_bimaterial_guard 2>&1 >/dev/null
rc=0
make test-fault-face-flux-bimaterial-guard 2>&1 | tee /tmp/guard_case2.log || rc=$?
# Restore original source
cp /tmp/fault_face_flux.cpp.bak dynamic/fault_face_flux.cpp
rm -f dynamic/fault_face_flux.o

if grep -q "^FAIL: " /tmp/guard_case2.log; then
   echo "META-PASS: guard-removed case FAILs as expected"
else
   echo "META-FAIL: guard-removed case should FAIL but was accepted"
   exit 1
fi
```

---

### [R-202] CRITICAL [wave_operator.inl:1194] — Shared-fault branch uses `MFEM_ASSERT` instead of `MFEM_VERIFY` for a bounds check that precedes an out-of-bounds access

**Category:** BUG / ASSUMPTION

**Description:**
In the shared-fault branch, the code reads `bd.qp_data[q]` after a `MFEM_ASSERT` bounds check:

```cpp
const FaultBasisData &bd = fault_basis_->GetBasis(basis_idx);
MFEM_ASSERT(q < static_cast<int>(bd.qp_data.size()),
            "FaultBasis::qp_data not populated for shared "
            "fault face — ComputeQPBasisShared missed "
            "this face");
const FaultBasisQPData &qpd = bd.qp_data[q];
```

`MFEM_ASSERT` is compiled out in release builds (it is only active under `MFEM_DEBUG`), whereas `MFEM_VERIFY` is always active.  In a release build, if `bd.qp_data.size() <= q` (e.g., the ctor's `ComputeQPBasisShared` failed to populate this face for some geometric reason), the assertion is a no-op and the next line (`bd.qp_data[q]`) performs an out-of-bounds `operator[]` on `std::vector` — **undefined behavior**.  Every subsequent field read (`qpd.normal`, `qpd.tangent1`, `qpd.tangent2`, `qpd.sign_flipped`) reads garbage, and `can_n / can_t1 / can_t2` become nonsense, the rotation matrices become non-orthogonal, and the per-side flux is silently wrong.

The interior-fault branch at `wave_operator.inl:740-744` has the correct analogue — an *explicit runtime* `if (q < static_cast<int>(bd.qp_data.size())) { qpd_ptr = &bd.qp_data[q]; }` — which then flows through the outer `if (dof_idx >= 0 && ... && qpd_ptr != nullptr)` to the `MFEM_ABORT` fallback.  The shared-fault branch has no such safety net in release.

This asymmetry was not flagged in round 1.  It is post-v9.0.0 load-bearing because the shared-fault branch now uses `qpd.normal / qpd.tangent1 / qpd.tangent2` to reconstruct `can_n / can_t1 / can_t2` (lines 1204-1213) — a wrong `can_n` corrupts the per-side flux AND the cross-rank consistency guarantee.

**Trigger:**
A release build (`MFEM_DEBUG` not defined) where `FaultBasis::ComputeQPBasisShared` left one shared fault face's `qp_data` under-populated (e.g., a geometric edge case, a stale mesh-regeneration, or a ctor-ordering bug).

**Actual behavior:**
In release: `bd.qp_data[q]` returns garbage; rotation matrices are computed from nonsense `normal / tangent1 / tangent2`; per-side flux on that QP is wrong; cross-rank consistency breaks.  No error is emitted.

In debug: `MFEM_ASSERT` fires loudly; the bug surfaces.

**Expected behavior:**
Same symmetric hard-abort behaviour in release as in debug — and the same behaviour as the interior-fault branch.

**Suggested fix:**
Promote the assertion to `MFEM_VERIFY`:
```diff
diff --git a/miniapps/seas/dynamic/wave_operator.inl b/miniapps/seas/dynamic/wave_operator.inl
@@
                   const FaultBasisData &bd = fault_basis_->GetBasis(basis_idx);
-                  MFEM_ASSERT(q < static_cast<int>(bd.qp_data.size()),
-                              "FaultBasis::qp_data not populated for shared "
-                              "fault face — ComputeQPBasisShared missed "
-                              "this face");
+                  MFEM_VERIFY(q < static_cast<int>(bd.qp_data.size()),
+                              "FaultBasis::qp_data not populated for shared "
+                              "fault face — ComputeQPBasisShared missed "
+                              "this face (basis_idx=" << basis_idx
+                              << ", q=" << q
+                              << ", qp_data.size()=" << bd.qp_data.size()
+                              << "). v9.0.0 Pelties-9 per-side flux reads "
+                              "qpd.normal/tangent1/tangent2 from this "
+                              "entry; a wrong can_n corrupts the flux.");
                   const FaultBasisQPData &qpd = bd.qp_data[q];
```

**Test case:**
```cpp
// tests/unit/test_shared_fault_qp_data_guard.cpp
//
// Verifies that the shared-fault QP-data bounds check is active in
// release builds (MFEM_VERIFY, not MFEM_ASSERT).  Without the fix,
// the test compiles but `bd.qp_data[q]` would be UB.

// Friend-class hook (add to WaveOperator private section under
// #ifdef SEAS_TEST_HOOKS):
//    std::vector<FaultBasisData> *TestGetFaultBasisRaw() {
//       return const_cast<FaultBasis*>(fault_basis_)->TestGetBasisRaw();
//    }
// And in FaultBasis:
//    #ifdef SEAS_TEST_HOOKS
//    std::vector<FaultBasisData> *TestGetBasisRaw() { return &basis_; }
//    #endif

int main()
{
   // ... build 2-rank parallel ParMesh with a shared fault face at y=0 ...
   WaveOperator op(...);
   op.SetFaultFlux(&ff);
   op.SetFaultDOFData(&dof_data, nqp_per_face);

   // Truncate qp_data on one shared fault face.
   auto *basis_vec = op.TestGetFaultBasisRaw();
   int n_int = op.FaultInteriorFacesSize();
   (*basis_vec)[n_int].qp_data.clear();   // force qp_data.size() = 0

   Vector Q(op.GetStateSize()), rhs(op.GetStateSize());
   Q = 0.0; rhs = 0.0;

   // Expect MFEM_VERIFY to abort in RELEASE (MFEM_DEBUG not defined)
   // as well as DEBUG.  Without the fix, release builds would silently
   // proceed with UB.  Run as death test; parent treats non-zero rc
   // as PASS.
   op.Apply(Q, rhs);
   std::cerr << "FAIL: shared-fault qp_data guard did not abort\n";
   return 0;   // 0 = guard did NOT fire = FAIL (per R-201 convention)
}
```
Shell wrapper (death-test with exit-code discrimination):
```bash
# test-shared-fault-qp-data-guard
@./seas_test_shared_fault_qp_data_guard 2>/dev/null; rc=$$?; \
    if [ $$rc -eq 0 ]; then \
        echo "FAIL: qp_data guard did not fire (MFEM_ASSERT regressed to no-op in release)"; \
        exit 1; \
    else \
        echo "PASS: qp_data guard fired (rc=$$rc)"; \
    fi
```

---

### [R-203] MODERATE [test_fault_flux_interior_vs_shared_branch_equivalence.cpp:267-273] — Stale Case B label still references the obsolete "Godunov identity" premise

**Category:** QUALITY / DEVIATION

**Description:**
Round-1 R-004 flagged the entire test as a tautology and recommended updating its documentation.  The **header comment** was updated (line 26-38 now states "POST-FIX THIS TEST IS A TAUTOLOGY") — but the *inline* Case A/B labels at line 261-273 still describe the pre-fix algorithmic dependency:
```cpp
   // Case A: Elem1 on canonical + side (rank-local nor = +can_n).
   //   Interior branch's arguments match shared branch trivially — MATCH.
   CompareBranches(can_n, can_t1, can_t2, Q_self_g, Q_nbr_g,
                   /*elem1_on_plus=*/true,
                   "Case A: elem1_on_plus=true (trivial match)");

   // Case B: Elem1 on canonical - side (rank-local nor = -can_n).
   //   Interior branch invokes flux.Interior(-can_n, Q_-g, Q_+g).
   //   Shared branch invokes flux.Interior(+can_n, Q_+g, Q_-g) with accum_sign=-1.
   //   Equality ⇔ Godunov identity (§3.1c).
   CompareBranches(can_n, can_t1, can_t2, Q_self_g, Q_nbr_g,
                   /*elem1_on_plus=*/false,
                   "Case B: elem1_on_plus=false (Godunov identity)");
```
Post-fix, the interior branch does NOT invoke `flux.Interior(-can_n, Q_-g, Q_+g)` — it uses `+can_n` with per-side `(Q, Q)` calls (line 135-136).  The shared branch similarly uses `+can_n` with a single per-side call (line 165).  No "Godunov identity" is invoked.  The label `"Case B: elem1_on_plus=false (Godunov identity)"` is outright wrong — the test's console output still advertises a dependency that was removed.  A reader skimming a CI log with "Case B (Godunov identity)" will mis-attribute the test's semantics.

**Trigger:**
Read the test source or its console output.

**Actual behavior:**
Labels advertise a behaviour the test no longer exercises.

**Expected behavior:**
Labels reflect the actual post-fix behaviour (trivial tautology at quiescent bulk; non-trivial only with `+perturb` cases).

**Suggested fix:**
```diff
diff --git a/miniapps/seas/tests/unit/test_fault_flux_interior_vs_shared_branch_equivalence.cpp b/miniapps/seas/tests/unit/test_fault_flux_interior_vs_shared_branch_equivalence.cpp
@@
    // Case A: Elem1 on canonical + side (rank-local nor = +can_n).
-   //   Interior branch's arguments match shared branch trivially — MATCH.
+   //   Post-fix: both branches call flux.Interior(can_n, Q_imp_plus_g,
+   //   Q_imp_plus_g) with the same selection — trivial bit-identity.
    CompareBranches(can_n, can_t1, can_t2, Q_self_g, Q_nbr_g,
                    /*elem1_on_plus=*/true,
-                   "Case A: elem1_on_plus=true (trivial match)");
+                   "Case A: elem1_on_plus=true (tautology, + side)");

    // Case B: Elem1 on canonical - side (rank-local nor = -can_n).
-   //   Interior branch invokes flux.Interior(-can_n, Q_-g, Q_+g).
-   //   Shared branch invokes flux.Interior(+can_n, Q_+g, Q_-g) with accum_sign=-1.
-   //   Equality ⇔ Godunov identity (§3.1c).
+   //   Post-fix: both branches call flux.Interior(can_n, Q_imp_minus_g,
+   //   Q_imp_minus_g) with the same selection — trivial bit-identity.
+   //   NO "Godunov identity" is exercised any longer; see header comment
+   //   and REVIEW R-004 / R-203 for the pre-fix vs post-fix divergence.
    CompareBranches(can_n, can_t1, can_t2, Q_self_g, Q_nbr_g,
                    /*elem1_on_plus=*/false,
-                   "Case B: elem1_on_plus=false (Godunov identity)");
+                   "Case B: elem1_on_plus=false (tautology, − side)");
```

**Test case:**
Comment-only edit — a compile-time grep suffices:
```bash
! grep -qE 'Godunov identity|flux\.Interior\(-can_n' \
    miniapps/seas/tests/unit/test_fault_flux_interior_vs_shared_branch_equivalence.cpp \
    || { echo "FAIL: stale Godunov-identity comment still present"; exit 1; }
echo "PASS: test labels reflect post-fix behaviour"
```

---

### [R-204] MODERATE [wave_operator.inl:710-717] — `MFEM_VERIFY` for fault bookkeeping is re-checked inside the QP loop; silent perf penalty, noisy on failure

**Category:** QUALITY / POSSIBLE BUG

**Description:**
The round-1 R-002 guard is placed at `wave_operator.inl:710-717`, inside the `for (int q = 0; q < nqp; q++)` loop at `:595` and inside the face loop at `:562`.  For a fault simulation with tens of thousands of fault faces × ~16 QPs × 6 RK stages × millions of time steps, the `MFEM_VERIFY` is evaluated billions of times — once per (face, QP, stage, step) combination.  The check is cheap (predicted branch on a non-null pointer), so the happy-path cost is negligible on modern CPUs; but if the condition ever fails, the diagnostic fires *per QP*, producing thousands of identical error messages before abort (noise that buries the underlying ctor-population bug).

More importantly, this placement makes the guard *conditional on reaching a fault face* — which is the right scope for catching "fault_flux_ null when a fault face appears" but means the guard is NOT hoisted to the ctor (where it would catch the bug at setup time, before any time-stepping).  A follow-up caller adding a new fault-flux entry point elsewhere would not inherit this guard.

The correct design is a one-time check at the boundary between setup and first face-loop entry — e.g., hoisted to `ComputeInteriorFaceFluxRHS`'s preamble (once per Mult call) or to the ctor (once per WaveOperator construction).

**Trigger:**
Any fault simulation — the guard fires per-QP in the happy case (performance), per-face-per-QP on failure (noise).

**Actual behavior:**
Billions of successful `MFEM_VERIFY` evaluations in the happy case.  On failure, thousands of duplicate abort messages before the process dies.

**Expected behavior:**
Single one-time check per `Mult` call (or per ctor).

**Suggested fix:**
Hoist the check to the top of `ComputeInteriorFaceFluxRHS`, before the face loop.  Keep the per-QP inner check only for debug builds via `MFEM_ASSERT`:
```diff
diff --git a/miniapps/seas/dynamic/wave_operator.inl b/miniapps/seas/dynamic/wave_operator.inl
@@
 template <typename MeshType>
 void WaveOperator<MeshType>::ComputeInteriorFaceFluxRHS(
    const Vector &Q, Vector &rhs) const
 {
+   // R-002 (hoisted from inner QP loop in round 2 / R-204): check fault
+   // bookkeeping ONCE per call, not per (face, QP).  `bc_.fault_attr > 0`
+   // is the condition under which the interior-face loop can enter the
+   // fault branch; if so, fault_flux_ and fault_dof_data_ must be set.
+   if (bc_.fault_attr > 0)
+   {
+      MFEM_VERIFY(fault_flux_ && fault_dof_data_,
+                  "WaveOperator::ComputeInteriorFaceFluxRHS: "
+                  "bc_.fault_attr=" << bc_.fault_attr
+                  << " > 0 (fault mesh configured) but "
+                  "fault_flux_=" << (void*)fault_flux_
+                  << ", fault_dof_data_=" << (void*)fault_dof_data_
+                  << "; ctor did not populate fault bookkeeping. "
+                  "v9.0.0 Pelties-9 per-side flux requires both; "
+                  "welded-flux fallback is no longer physical.");
+   }
    ...
    for (int f = 0; f < mesh_.GetNumFaces(); f++)
    {
       ...
       for (int q = 0; q < nqp; q++)
       {
          ...
             if (is_fault)
             {
-               // R-002: a fault face must have fully populated bookkeeping.
-               // Post-v9.0.0 the welded-flux fallback (outer else below) is
-               // no longer physical — per-side flux is required for fault
-               // radiation, so silently falling back would under-radiate
-               // exactly the way R-F02 guarded against for the inner else.
-               MFEM_VERIFY(fault_flux_ && fault_dof_data_,
-                           "interior fault face f=" << f
-                           << " reached ComputeInteriorFaceFluxRHS but "
-                           "fault_flux_=" << (void*)fault_flux_
-                           << ", fault_dof_data_=" << (void*)fault_dof_data_
-                           << "; ctor did not populate fault bookkeeping. "
-                           "v9.0.0 Pelties-9 per-side flux requires both; "
-                           "welded-flux fallback is no longer physical.");
+               MFEM_ASSERT(fault_flux_ && fault_dof_data_,
+                           "R-204: bookkeeping hoisted check should have "
+                           "fired at the top of ComputeInteriorFaceFluxRHS. "
+                           "If we reach here with null fault_flux_, the "
+                           "hoisted check was removed or bypassed.");
                // Fault face: dispatch to FaultFaceFlux with tracked DOFData.
```
Parallel change in `ComputeSharedFaceFluxRHS` for symmetry (hoist the analogous check for the shared-fault path).

**Test case:**
```cpp
// Assertion on perf behaviour: this fix is a hot-path optimization, not
// a correctness fix.  The CORRECTNESS regression test for R-002 is
// already covered by a prior "test_fault_face_null_bookkeeping_aborts"
// pattern (see round 1).  For R-204, verify the hoist does NOT change
// abort semantics:

// Death test: bc_.fault_attr > 0 but no SetFaultFlux call → MFEM_VERIFY
// should fire at the top of the first Mult call, before any face-loop
// iteration.  Observable: process aborts immediately, not after N
// partial QP writes to rhs.

int main()
{
   WaveOperator op(fes, bc_with_fault_attr, /*order=*/1);
   // Deliberately skip SetFaultFlux.

   Vector Q(op.GetStateSize()), rhs(op.GetStateSize());
   Q = 0.0; rhs = 42.0;   // sentinel value

   // Expect abort BEFORE rhs is touched; the hoisted check runs before
   // the face loop.  A per-QP check would have touched rhs on the first
   // non-fault face; the hoisted one should not.
   op.Apply(Q, rhs);

   // If we reach here, guard didn't fire — FAIL.
   std::cerr << "REGRESSION: hoisted R-002 guard did not fire\n";
   return 0;
}
```

---

### [R-205] LOW [test_fault_face_flux_bimaterial_guard.cpp:50] — Uninitialized `Q_imp_plus / Q_imp_minus` inspected after a guard-bypass regression

**Category:** QUALITY / EDGE_CASE

**Description:**
In the bimaterial-guard test, `Q_imp_plus[NUM_STATE]` and `Q_imp_minus[NUM_STATE]` are declared uninitialized at `:50`:
```cpp
real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
```
On the intended PASS path (guard fires), `Evaluate` aborts before writing to them — no read, no UB.  On the fall-through regression path (guard gone), `Evaluate` writes to them with (possibly junky) bimaterial output, and the test prints to stderr and returns.  No access to the uninitialized data occurs either way, so this is not a correctness bug *in the test* — but combined with R-201 (the fall-through doesn't FAIL), it paints a misleading picture where the test appears to run successfully on uninitialized buffers.  If a future edit adds a diagnostic print of `Q_imp_plus[VX]` in the fall-through path, it would print uninitialized data and cause a compiler warning or, with aggressive optimization, undefined behaviour.

**Trigger:**
A future edit adds a `std::cerr << Q_imp_plus[VX]` in the fall-through error message to help debug what Evaluate "computed".

**Actual behavior:**
`-Wuninitialized` / `-fsanitize=memory` would flag the access.

**Expected behavior:**
Either zero-initialize the buffers so any future debug print reads deterministic zeros, or explicitly mark them uninitialized with `[[maybe_unused]]`.

**Suggested fix:**
```diff
diff --git a/miniapps/seas/tests/unit/test_fault_face_flux_bimaterial_guard.cpp b/miniapps/seas/tests/unit/test_fault_face_flux_bimaterial_guard.cpp
@@
    real_t Q_plus[NUM_STATE] = {};
    real_t Q_minus[NUM_STATE] = {};
-   real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
+   real_t Q_imp_plus[NUM_STATE] = {};
+   real_t Q_imp_minus[NUM_STATE] = {};
    FaultFaceFlux ff(TPV102Params::rho,
```

**Test case:** n/a — defensive hygiene; no runtime behaviour change in the current code path.

---

## Summary

- Critical issues: **2** (R-201: bimaterial-guard death test is broken; R-202: shared-fault MFEM_ASSERT is a no-op in release, allowing OOB access)
- Moderate issues: **2** (R-203: stale Case-B label in §3.1d test; R-204: per-QP `MFEM_VERIFY` should hoist to ctor/call-time)
- Low issues: **1** (R-205: uninitialized `Q_imp_*` in bimaterial-guard test)
- Plan compliance: **PARTIAL-PLUS** — every source edit the plan and round-1 review prescribed is now present; however, R-004's tautology-comment fix is incomplete (R-203), and R-007's new test has an exit-code bug that defeats the test's purpose (R-201).  The shared-fault branch has a release-only UB risk the round-1 review missed (R-202).
- Verdict: **FAIL — must fix R-201 and R-202 before proceeding.**  R-201 leaves the bimaterial-guard regression gate disarmed.  R-202 admits release-only undefined behaviour on a code path the v9.0.0 fix made load-bearing (shared-fault per-side flux depends on `qpd.normal/tangent1/tangent2` being valid).  Both are simple mechanical fixes.  R-203–R-205 are cleanup and can be batched into the same commit or follow up.

## Unreviewed Areas

- Ctor population of `fault_basis_->basis_[].qp_data` for shared fault faces: R-202 assumes this *could* be incomplete.  I did not audit `ComputeQPBasisShared` to verify that every shared fault face always receives `ir.GetNPoints()` entries.  If ComputeQPBasisShared is provably total, R-202's severity drops from CRITICAL to MODERATE (defense-in-depth only).  A quick source audit of `fault_basis.hpp`'s `ComputeQPBasisShared` implementation before or after applying R-202's fix is recommended.
- Performance measurements of R-204's per-QP guard in a TPV102 production run: I did not profile; the "billions of checks" estimate is back-of-envelope.  The hoist is worth doing for correctness-of-design reasons (single check per call, consistent with ctor-invariant style) even if the perf delta is negligible.
- The existence of analogous `MFEM_ASSERT` bounds checks elsewhere in the v9.0.0-touched code: I searched `wave_operator.inl` but not every file.  A focused `grep -n "MFEM_ASSERT" miniapps/seas/dynamic/*.inl miniapps/seas/dynamic/*.cpp | grep -v "MFEM_ASSERT(true"` would flag any other spots where a runtime bounds check was downgraded to a debug-only assertion on a post-v9.0.0 hot path.
- Whether the round-1 R-004 recommendation's DELETE option (removing the §3.1d test entirely rather than merely annotating it as tautological) is preferred by the project owner — R-203's minimal comment fix is the low-risk route; a full replacement with a non-tautological branch-equivalence test (serial + parallel twin mesh) remains deferred.
