# Code Review: Phase 2 + Phase 3 — unified shared-fault substep dispatch — 2026-07-10

## Review Scope
- **Plan:** `document/PLAN_unify_interior_shared_fault_substep_2026-07-09.md` (Phase 2 + 3, per its
  "PHASE 2 + PHASE 3 IMPLEMENTED" banner) as amended by
  `debug_document/tpv104_debug_document/R1601_root_cause_2026-07-10.md`.
- **Files reviewed:** `dynamic/wave_operator.inl` (unified gate, v_imp port, deleted `[R1601-DIAG]`
  block, reconcile assertion), `dynamic/wave_operator.hpp` (History note, tol member + setter),
  `tests/parallel/test_shared_fault_substep_parity_np2.cpp`,
  `tests/parallel/diag_r1601_shared_substep_experiment.cpp`, `Makefile`.
- **Domain context:** CLAUDE.md (fix-revert protocol, sign conventions), the root-cause doc, the
  implementer's completion report (verification matrix 6/6), `ExchangeAndPairSharedFaultQPs`
  implementation (`wave_operator.inl:6278+`).

## Hunt-list verdicts (from the review brief) — verified against code

| # | Question | Verdict |
|---|---|---|
| 1 | `VerifyForcedRuptureTimeReady` only in the else-branch — correct? | **YES.** The interior branch has the identical structure (guard lives inside the inline `LSW_ForcedRupture` arm, `:4338-4348`; the consume path bypasses it there too). The guard protects `GetTime()`, which the buffer path never reads. Symmetric by construction. |
| 1b | DOFData writes the inline solve did that the iterator does not? | **None found.** `EvaluateADER_LSW`/`EvaluateADER` write `V/tau_corr/sigma_n_corr` (+`StoreImposedVelocity_`); the iterator's `WriteBackState` writes all of these (incl. `psi` on RS), and the v_imp port covers the last gap. Empirically pinned by np=10 O∈{2,3,4} runs with the Phase-3 assertion armed: zero aborts over 522 steps × 3. |
| 2 | `payload.data() + local_qp*NPAY` indexing vs exchange record order | **ALIGNED.** Payload is built by iterating `fault_qp_buf` in push order; `fault_qp_buf` is filled by the Pass-1 loop over ascending `sf`; the exchange builder (`:6324+`) iterates the same `fault_shared_faces_` ascending with the same skip conditions (`!ftr`, missing `dof_offset` — both of which would equally have prevented the Pass-1 push). The pre-existing callback already indexes `fault_qp_buf[local_qp]` on the same assumption. See R-004 for a cheap belt-and-braces assert. |
| 2b | Assertion firing on one rank only — deadlock? | **No deadlock.** Both ranks compare the same unordered pair (a,b), so the predicate is symmetric — both abort or neither. Even if only one aborted, `MFEM_ABORT` → `MPI_Abort` kills the job; no hang. |
| 2c | tol member non-mutable, read from const method | **Fine.** Const methods may read non-mutable members; only the (non-const) setter writes it. |
| 3 | Stale refs to deleted env vars / discard behavior | **Two found** (R-002, R-003). The diag harness's mention of the former env var is intentional history and reads correctly. |
| 4 | `[MACRO]` comment stale? | **Yes, doubly** (R-003). |
| 5 | Negative control for the reconcile assertion; one-shot shared coverage | **Assertion never demonstrated to fire** (R-001). One-shot shared path IS exercised (the oracle's `run(false)` at np=2 goes through the inline dispatch every invocation) but only as a baseline; acceptable. |
| 6 | v_imp port | **Correct fields/guard** (`v_imp_plus/minus[0..2]` ← `I_imp[VX/VY/VZ]/dt`, LSW-family-only guard mirrors interior). `1/dt` without a `dt>0` check is the pre-existing interior pattern (`EvaluateADER*` verifies `dt>0` on the inline path; the buffer path relies on the caller) — R-005, LOW. |

---

## Findings

### [R-001] MODERATE [wave_operator.inl:ComputeADERSharedFaceFluxRHS reconcile] — The new assertion has never been shown to fire; as shipped it is unproven guard code

**Category:** ASSUMPTION (missing negative control)

**Description:**
Phase 3's whole point is that cross-rank divergence "can never hide again". The assertion passed
0-abort across the entire verification matrix — which proves it doesn't *false*-trigger, but
nothing proves it *true*-triggers. An inverted comparison, a wrong tolerance scaling, or an
accidentally dead code path would produce exactly the same all-green result. Guard code that has
never fired is untested code.

**Trigger:** any genuine cross-rank disagreement — none exists in the healthy build, so a test must
manufacture one.

**Expected behavior:** a test that corrupts one rank's shared-QP state and observes the abort.

**Suggested fix:** the codebase already ships the tool: `TamperSharedFaultElem1OnPlus(idx)`
(`wave_operator.hpp:849`, `SEAS_TEST_INTERNAL`-gated). Flipping the routing bit on ONE rank swaps
that rank's ± slots, so its payload (I_imp at minimum) disagrees with the peer ⇒ the assertion
must abort. Add to the parity oracle (compiled with `-DSEAS_TEST_INTERNAL` — needs the obj rule
flag) or a small dedicated test:

```cpp
// Negative control (R-001): a manufactured cross-rank inconsistency must be
// CAUGHT by the Phase-3 reconcile assertion, not silently boss-overwritten.
if (g_nprocs == 2 && std::getenv("SEAS_TEST_RECONCILE_NEGATIVE") != nullptr)
{
   if (g_rank == 1) { wave.TamperSharedFaultElem1OnPlus(0); }
   // run one AdvanceADER with the substep buffer installed, inside a forked
   // child (AbortsInChild pattern from test_forced_rupture_resolver.cpp) or
   // simply as a separate env-gated Makefile invocation asserting exit != 0.
}
```
Simplest robust wiring: a second env-gated invocation in the `test-shared-fault-substep-parity-np2`
run target that expects non-zero exit:
```make
	@SEAS_TEST_RECONCILE_NEGATIVE=1 mpirun -np 2 ./seas_test_shared_fault_substep_parity_np2; \
	 if [ $$? -eq 0 ]; then echo "NEGATIVE CONTROL FAILED: reconcile assertion did not fire"; exit 1; \
	 else echo "  negative control OK: reconcile assertion fires on manufactured mismatch"; fi
```

**Test case:** as above — tamper on rank 1, expect abort; untampered run stays green.

---

### [R-002] LOW [Makefile:5871] — Stale comment: the parity target still "characterises the surviving R-1601 defect"

**Category:** QUALITY (stale doc that now asserts the opposite of the code's behavior)

**Description:** the target's comment says the test characterises the *surviving defect*; since the
Phase-2 flip the test asserts the defect is GONE (both rank counts consume). A reader grepping for
open defects will be misled.

**Suggested fix:**
```diff
-# Phase 0(b) ORACLE (review R-001): does the ADER corrector consume the
-# sub-step buffer on interior vs shared fault QPs?  Characterises the surviving
-# R-1601 defect.  -DSEAS_TEST_INTERNAL not needed (uses only public accessors).
+# Phase 0(b)/Phase 2 ORACLE: the ADER corrector must consume the sub-step
+# buffer on BOTH interior and shared fault QPs (unified dispatch; the R-1601
+# discard fallback was retired — see R1601_root_cause_2026-07-10.md).
```

---

### [R-003] LOW [wave_operator.inl:[MACRO] diagnostic block] — Two stale claims in the env-gated diagnostic's comment

**Category:** QUALITY

**Description:** the `SEAS_DIAG_MACRO` block's comment (immediately after the unified gate) still
says (a) the printed traction decomposition is "exactly what EvaluateADER_LSW just consumed" —
false whenever the buffer was consumed (the inline solve did not run; `fdata.sigma_n_corr` is the
iterator's write-back), and (b) "this shared +/- routing uses elem1_on_plus while the iterator uses
sign_flipped" — a frame claim the Phase-0 census showed is not a meaningful producer/consumer
mismatch (both use `sign_flipped` for the frame; `elem1_on_plus` is the routing bit on both sides).
Diagnostic-only, but it planted the wrong root cause once already.

**Suggested fix:**
```diff
-                  // traction decomposition on the TIME-INTEGRATED canonical +/-
-                  // state (I_{plus,minus}_local / dt — exactly what EvaluateADER_
-                  // LSW just consumed via ComputeTrialTraction) and the WRITTEN
-                  // output fdata.sigma_n_corr / slip_rate.
+                  // traction decomposition on the TIME-INTEGRATED canonical +/-
+                  // state (I_{plus,minus}_local / dt — what the inline one-shot
+                  // path would consume; on the unified substep path the flux
+                  // came from the iterator's buffer instead) and the WRITTEN
+                  // output fdata.sigma_n_corr / slip_rate (iterator write-back).
```
and delete the "NB: this shared +/- routing uses elem1_on_plus while the iterator uses
sign_flipped" sentence (superseded by the census; cite
`tests/parallel/test_fault_frame_bit_census.cpp` instead).

---

### [R-004] LOW [wave_operator.inl:reconcile] — Belt-and-braces: assert the payload/callback index alignment the new code inherits

**Category:** ASSUMPTION

**Description:** `mine = payload.data() + local_qp*NPAY` inherits the pre-existing assumption that
the exchange builder's record order equals `fault_qp_buf` order. Verified true today (both ascend
`sf` with identical skip conditions), but the two loops live 600 lines apart and nothing enforces
it. A future `continue` added to one loop but not the other would silently compare wrong pairs —
and the assertion would then fire with a *misleading* message (or worse, pass wrongly).

**Suggested fix:** one cheap structural check inside the callback:
```diff
    const double *mine = payload.data() +
                         static_cast<size_t>(local_qp) * NPAY;
+   MFEM_ASSERT(static_cast<size_t>((local_qp + 1) * NPAY) <= payload.size(),
+               "reconcile: local_qp " << local_qp << " out of payload range");
```
plus a comment at the payload-build loop naming the ordering contract with the exchange builder.

---

### [R-005] LOW [wave_operator.inl:v_imp recovery, shared branch] — `1/dt` without a `dt > 0` verify

**Category:** EDGE_CASE

**Description:** the inline path's `EvaluateADER*` entry points verify `dt > 0`; the consume path
divides by `dt` (v_imp recovery) with no check. Identical to the pre-existing interior consume
branch, so this is inherited, not introduced — but the unified change doubled the surface.
`AdvanceADER` callers always pass the macro `dt > 0` in practice.

**Suggested fix:** either add `MFEM_ASSERT(dt > 0.0, ...)` at the top of
`ComputeADERSharedFaceFluxRHS` (and its interior sibling, same line of reasoning), or explicitly
document the caller contract at both v_imp sites. Assert preferred:
```diff
 void WaveOperator<MeshType>::ComputeADERSharedFaceFluxRHS(...)
 {
+   MFEM_ASSERT(dt > 0.0,
+               "ComputeADERSharedFaceFluxRHS: dt must be > 0; got " << dt);
```

---

## Summary
- Critical issues: **0**
- Moderate issues: **1** (R-001 — unproven assertion)
- Low issues: **4** (R-002, R-003, R-004, R-005)
- **Plan compliance: FULL.** Phase 2 and Phase 3 are implemented as specified (including the
  cleanup mandate: env block, discards, and stale comments retired; grep confirms zero live
  references to `SEAS_DIAG_SHARED_SUBSTEP*`). The verification matrix is genuinely 6/6, and the
  np=10/O=2 exact match to np=1 (12.7387) exceeds the plan's acceptance bar.
- **Verdict: PASS WITH FIXES.** R-001 should land before this is called done — an assertion that
  has never fired is the same category of latent risk (untested guard) that let R-1601's
  mis-attribution survive for months. R-002..R-005 are cheap hygiene.

## Fix Status (applied 2026-07-10, post-review)

| ID | Status | Notes |
|---|---|---|
| R-001 | **FIXED** + verified | Negative control added (`SEAS_TEST_RECONCILE_NEGATIVE`, tamper on rank 1, `-DSEAS_TEST_INTERNAL` + `FAULT_FACE_FLUX_TESTINTERNAL_OBJ` link); wired into the `make` run target expecting non-zero exit. **Deviation from the review's sketch, with cause:** the tamper is invisible on the sentinel-buffer leg (payloads identical by construction) AND on a uniform field (swapping two identical vectors) — the control runs the INLINE leg with a rank-asymmetric velocity field. First attempt with uniform Q reached the fallthrough (assertion did NOT fire), proving the control itself needed the asymmetry — exactly the class of dead-guard bug R-001 exists to catch. Final state: assertion fires with `payload field 4 = V2`; fallthrough never reached; positive runs unaffected (3/3, 6/6). |
| R-002 | **FIXED** | Makefile comment rewritten (unified contract, cites the root-cause doc). |
| R-003 | **FIXED** | Both stale `[MACRO]` claims corrected; frame/routing note now cites the census. |
| R-004 | **FIXED** | Payload-range `MFEM_ASSERT` in the callback + ORDERING CONTRACT comment at the payload build loop. |
| R-005 | **FIXED** (deviation) | `dt > 0` verify added to `ComputeADERSharedFaceFluxRHS` only — `ComputeADERFaceFluxRHS` **already had it** (`:3918`); the review's "add to both" was half-stale. |

**Post-fix verification:** `make test-shared-fault-substep-parity-np2` = positive 3/3 + 6/6 +
negative control fires; census PASS both rank counts; TPV26 smoke `V_max = 0.185543` (np=1
byte-exact); symmirror np=10 O=2 clean, `V_max = 12.7387` (== np=1), 0 reconcile aborts.

## Unreviewed Areas
- The np≥4 attractor (station 4.8123 vs 5.0022) — explicitly out of scope; pre-existing, tracked
  in the plan banner and memory as the Phase-5 problem.
- Bimaterial operator (`BimaterialWaveOperator`) inherits the changed shared branch via the same
  `.inl`; not separately exercised beyond the existing bimaterial test suite (not run here).
- Full `make test` (~230 MFEM-linking targets) not run; the directly-affected set plus the np=10
  production-scale sweep was.
