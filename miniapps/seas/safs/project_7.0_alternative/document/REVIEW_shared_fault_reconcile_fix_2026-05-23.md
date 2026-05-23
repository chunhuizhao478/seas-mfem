# Code Review: shared-fault cross-rank reconcile (Phases 0–3) — 2026-05-23

## Review Scope
- **Plan:** `miniapps/seas/safs/project_7.0_alternative/document/PLAN_shared_fault_reconcile_fix_2026-05-23.md`
- **Commits reviewed (chronological):**
  - `6a3af3e` Phase 0 — `SEAS_DIAG_XRANK` `%.17e` + RAWDOF trace; oblique non-fatal sbatch
  - `e2a2554` Phase 0 verdict (source-1 interpolation) + Phase 1 oracle (injection hook, test, Makefile)
  - `a26f6ac` Phase 2 — `ExchangeAndPairSharedFaultQPs` + two-pass reconcile in `ComputeADERSharedFaceFluxRHS`
  - `7761baf` Phase 3 — `s_seas_test_disable_reconcile` switch + negative test leg + abort-wording update
- **Files reviewed:** `dynamic/wave_operator.{hpp,inl}`, `dynamic/fault_face_flux.{hpp,cpp}`,
  `dynamic/tpv205_substep_iterator.cpp`, `drivers/spatial_dyn_driver.cpp`,
  `tests/unit/test_shared_fault_reconcile_cross_rank.cpp`, `Makefile`,
  `jobs/safs/spatial_dyn_resDc2_XRANKdiag_8N_400r_dev_2hr_safs.sbatch`
- **Domain context:** `CLAUDE.md` (project + miniapp), the plan above, the verify rationale in
  `VerifySharedFaultDOFDataConsistency`.

## Headline finding — WHY the reconcile passes the local test but has no effect on the SAFS production run

The reconcile **does run** in production (`AdvanceADER` → `ComputeADERSharedFaceFluxRHS`
is called unconditionally for parallel meshes, `wave_operator.inl:5255`; the driver binary
*is* rebuilt against `wave_operator.inl` — the Makefile dependency `$(SPATIAL_DYN_DRIVER_OBJ):
… $(SEAS_HEADERS)` covers it, and `SEAS_HEADERS ⊃ DYNAMIC_HEADERS ∋ dynamic/wave_operator.inl`).
So this is **not** a stale-binary problem.

The gap is that **the reconcile, its runtime guard, and the Phase-1 test oracle all operate on
the 8-field set `{tau1_corr, tau2_corr, sigma_n_corr, V1, V2, psi, slip1, slip2}`, while every
quantity the SAFS production run is *judged by* — `V_max`, the rupture-area proxy, the
`[DIAG-ONSET]` dump, and the `fault.vtkhdf` slip-rate "speckle" — is read from
`DOFData::slip_rate`, which is NOT in the payload and NOT checked by the verify** (R-001, R-002).
Compounding this, **the local test never runs the production code path at all** — it calls
`AdvanceADER` directly and bypasses the substep iterator that is production's authoritative
fault-state evolver, so its GREEN result says nothing about the production substep path (R-003).

The result: in production the reconcile silently makes `V1/V2/slip1/slip2` cross-rank-identical,
but leaves `slip_rate` per-rank-divergent — so the diagnostics the user reads are unchanged
("no effect"), while the local oracle (which measures the reconciled fields) is GREEN.

---

## Findings

### [R-001] CRITICAL [wave_operator.inl:ComputeADERSharedFaceFluxRHS reconcile callback] — `slip_rate` is omitted from the reconcile payload; it is the field every production diagnostic reads

**Category:** BUG

**Description:**
`DOFData::slip_rate` (the scalar `|V|`, `fault_face_flux.hpp:50`) is written per-rank by
`WriteBackState` (`data.slip_rate = s.V_abs;`, `fault_face_flux.cpp:307`) on *both* the
substep iterator path (`StepOneQP_` → `WriteBackState`, `tpv205_substep_iterator.cpp:133`) and
the inline Pass-1 solve in `ComputeADERSharedFaceFluxRHS`
(`EvaluateADER_LSW_ForcedRupture` → `WriteBackState`). The reconcile payload (`wave_operator.inl`
~`5113-5141`) and the overwrite callback (~`5148-5163`) carry only the 8 fields
`{tau1_corr, tau2_corr, sigma_n_corr, V1, V2, psi, slip1, slip2}` plus `I_imp_{plus,minus}_can`.
`slip_rate` is never copied.

After the reconcile, the non-boss rank holds the boss's `V1/V2` but **its own** `slip_rate`
— so `slip_rate ≠ sqrt(V1² + V2²)` on the non-boss, and at the LSW kink (one rank slips,
the other locks) the two ranks' `slip_rate` differ by O(1).

Every SAFS production observable reads exactly this un-reconciled field:
- `V_max` blow-up monitor: `V_max_local = max(V_max_local, dof_data[i].slip_rate)` (`spatial_dyn_driver.cpp:1891`)
- rupture-area proxy: `if (d.slip_rate > 0.5) { n_rup_local++; }` (`spatial_dyn_driver.cpp:1916`)
- `[DIAG-ONSET]` dump prints `d.slip_rate` (`spatial_dyn_driver.cpp:1944`)
- the ParaView/`fault.vtkhdf` smoke verifier thresholds the `slip_rate` field
  (`spatial/code/scripts/verify_spatial_dyn_smoke_safs.py:185-196`)

So the run can report a `V_max`/`n_rupturing` blow-up driven by a non-boss QP whose **physics
(`V1/V2` → `I_imp` → bulk) was already reconciled and bounded** — the diagnostic lies, and the
user sees "no effect."

**Trigger:**
np ≥ 2, LSW / LSW_ForcedRupture, a shared fault QP parked at the slip-onset threshold during the
gradual_overstress window — i.e. exactly the SAFS Dc2/zerodip runs. The two ranks split at the
kink; the reconcile fixes `V1/V2` but not `slip_rate`.

**Actual behavior:**
`dof_data[shared].slip_rate` retains the non-boss's per-rank `V_abs`; `V_max`, `n_rupturing`,
`[DIAG-ONSET]`, and the `fault.vtkhdf` speckle remain cross-rank-divergent post-reconcile.

**Expected behavior:**
`slip_rate` must equal the boss's value (`= sqrt(V1² + V2²)`) on both ranks after the reconcile,
so the diagnostics reflect the reconciled physics.

**Suggested fix:** add `slip_rate` to the payload and the callback (preferred — bit-identical to
boss), and add it to the verify so a future regression is caught (see R-002).

```diff
@@ payload build (wave_operator.inl ~5126) @@
             payload.push_back(d.slip1);
             payload.push_back(d.slip2);
+            payload.push_back(d.slip_rate);
             for (int c = 0; c < NUM_STATE; c++) { payload.push_back(qa.I_imp_plus_can[c]); }
             for (int c = 0; c < NUM_STATE; c++) { payload.push_back(qa.I_imp_minus_can[c]); }
@@ NPAY (wave_operator.inl ~5111) @@
-         constexpr int NPAY = 8 + 2 * NUM_STATE;  // 8 DOFData fields + I_imp +/-
+         constexpr int NPAY = 9 + 2 * NUM_STATE;  // 9 DOFData fields (incl slip_rate) + I_imp +/-
@@ overwrite callback (wave_operator.inl ~5160) @@
               d.slip1        = peer[6];
               d.slip2        = peer[7];
-              for (int c = 0; c < NUM_STATE; c++) { qa.I_imp_plus_can[c]  = peer[8 + c]; }
-              for (int c = 0; c < NUM_STATE; c++) { qa.I_imp_minus_can[c] = peer[8 + NUM_STATE + c]; }
+              d.slip_rate    = peer[8];
+              for (int c = 0; c < NUM_STATE; c++) { qa.I_imp_plus_can[c]  = peer[9 + c]; }
+              for (int c = 0; c < NUM_STATE; c++) { qa.I_imp_minus_can[c] = peer[9 + NUM_STATE + c]; }
```

(Equivalent minimal alternative — recompute instead of copy — if you prefer not to grow the
payload: in the callback after `d.V2 = peer[4];` add
`d.slip_rate = std::sqrt(d.V1*d.V1 + d.V2*d.V2);`. This is exact since
`SolveLSW_TPV205`/the rate-state solve define `V_abs = sqrt(V1²+V2²)`. Copying via the payload
is preferred because it is bit-identical to the boss and survives any future change to that
identity.)

**Test case (C++, np=2 — extends the Phase-1 oracle):**
```cpp
// In RunLeg(...), after each step, ALSO gather slip_rate across ranks and assert bit-identity.
// Reuse the verifier's pairing but compare DOFData::slip_rate. Concretely, add a second metric:
//   double wr_sliprate = MaxCrossRankRelDiff(dof_data, &DOFData::slip_rate);  // helper to add
//   res.max_wr_sliprate = std::max(res.max_wr_sliprate, wr_sliprate);
// Assert (post-fix): TEST_TRUE(g_wr_sliprate == 0.0,
//   "cross-rank slip_rate bit-identical after reconcile — LSW");
// RED before this fix (LSW leg: one rank slips, slip_rate=2.16, other locks slip_rate=0 -> wr=1.0),
// GREEN after.
```

---

### [R-002] CRITICAL [wave_operator.inl:VerifySharedFaultDOFDataConsistency + test oracle] — the guard and oracle check a different field set than production reads, so they cannot detect R-001

**Category:** DEVIATION (oracle blind spot)

**Description:**
`VerifySharedFaultDOFDataConsistency` packs exactly `NUM_FIELDS = 8`
(`FIELD_NAMES = {tau1_corr, tau2_corr, sigma_n_corr, V1, V2, psi, slip1, slip2}`,
`wave_operator.inl:5914-5919`). The Phase-1 test asserts `worst_rel == 0` from this verify
(`test_shared_fault_reconcile_cross_rank.cpp:386`). Neither checks `slip_rate`. Therefore the
oracle is **structurally incapable** of catching R-001: the field that drives the production
verdict is invisible to both the runtime guard and the test. This is the direct reason "passes
local test" ≠ "fixes production": the test measures the reconciled fields, the production run is
judged by an un-reconciled one.

The plan's own failure-mode note (PLAN §"Reading the result": "`worst_rel ≤ 1e-13` but the
speckle SURVIVES ⇒ a co-present cause (under-resolution / SSO)") will **mis-attribute** this bug:
the true cause is the un-reconciled `slip_rate`, not under-resolution.

**Trigger:** any post-fix run — the guard reports `max_rel_diff == 0` while `slip_rate` diverges.

**Actual behavior:** verify (and test) GREEN while a mutable, output-bound field is divergent.

**Expected behavior:** the verify's field set must be a **superset** of every mutable per-step
`DOFData` field that is (a) carried across steps, (b) visualized, or (c) used by a driver
diagnostic — minimally it must include `slip_rate`.

**Suggested fix:** add `slip_rate` to the verify (`REC` 16→17, `NUM_FIELDS` 8→9), so the guard
trips if any field the run reads ever diverges.

```diff
@@ wave_operator.inl ~5908 @@
-      constexpr int REC = 16;
+      constexpr int REC = 17;
@@ ~5913 @@
-      constexpr int NUM_FIELDS = 8;
-      constexpr int RANK_OFFSET = 15;
+      constexpr int NUM_FIELDS = 9;
+      constexpr int RANK_OFFSET = 16;
@@ ~5916 @@
       static const char *FIELD_NAMES[NUM_FIELDS] = {
          "tau1_corr", "tau2_corr", "sigma_n_corr",
-         "V1", "V2", "psi", "slip1", "slip2"
+         "V1", "V2", "psi", "slip1", "slip2", "slip_rate"
       };
@@ pack loop ~5969 @@
          local_data.push_back(d.slip1);
          local_data.push_back(d.slip2);
+         local_data.push_back(d.slip_rate);
          local_data.push_back(static_cast<double>(my_rank_));
```
(Keep the reconcile payload and the verify field set in lockstep — both must list `slip_rate`.)

**Test case:** the Phase-1 negative leg already proves the guard trips on a real desync; after
adding `slip_rate` to the verify, the R-001 test above (slip_rate divergence) must also trip the
guard when the reconcile is disabled, and report 0 when enabled.

---

### [R-003] CRITICAL [tests/unit/test_shared_fault_reconcile_cross_rank.cpp] — the local oracle bypasses the substep iterator that performs production's authoritative fault-state evolution

**Category:** DEVIATION (test does not exercise the production path)

**Description:**
Production advances the fault via
`spatial_dyn_driver.cpp:1880 AdvanceADERWithSubStep_Spatial`
→ `iterator.AdvanceWithSubStepStates(...)` (`spatial_dyn_driver.cpp:439`,
`tpv205_substep_iterator.cpp:262`) → `StepOneQP_` **then** `wave.AdvanceADER(...)`
(`spatial_dyn_driver.cpp:460`). The substep iterator is the **authoritative** evolver of
`slip1/slip2` (LSW weakening state): `EvaluateADER_LSW(_ForcedRupture)` deliberately does **not**
accumulate slip — "slip evolution is owned exclusively by the iterator"
(`fault_face_flux.cpp:815-824`).

The test calls `wave.AdvanceADER(...)` **directly** (`test_…:309`) and never constructs or runs a
`Tpv205SubStepIterator`. Consequences:
1. In the test, `slip1/slip2` are **never advanced** (only the iterator advances them, and it is
   absent) — they remain 0 on both ranks for all steps. So the reconcile's copy of `slip1/slip2`
   is a `0→0` no-op and is **not actually tested**, even though the plan's Phase-1 AC claims to
   validate cross-rank slip consistency.
2. The test exercises only the `ComputeADERSharedFaceFluxRHS` inline path, which in production is
   a **secondary** macro-`dt` re-solve (the R-1601 fallback, `wave_operator.inl:4755-4812`)
   layered on top of the iterator's substep solve. The interaction iterator-writes-then-inline-
   re-solves-then-reconciles (and the slip_rate inconsistency it produces, R-001) is untested.

This is why a GREEN local result does not imply a fixed production run.

**Trigger:** running the production driver vs. running the unit test — different call graphs.

**Actual behavior:** test validates a path (`AdvanceADER`-direct, slip frozen at 0) that does not
represent production.

**Expected behavior:** at least one oracle must drive the iterator path (`AdvanceWithSubStepStates`
→ `AdvanceADER`) so that (a) `slip1/slip2` actually evolve and their reconcile is exercised, and
(b) `slip_rate` consistency post-reconcile is checked (R-001/R-002).

**Suggested fix:** add a third leg to the test that mirrors `AdvanceADERWithSubStep_Spatial`:
construct a `Tpv205SubStepIterator`, call `ComputeADERSubStepStates` /
`EvaluateBulkAtFaultQPsCanonical` / `AdvanceWithSubStepStates` / `SetSubStepFaultImposedStates`
before `AdvanceADER` (the driver wrapper at `spatial_dyn_driver.cpp:365-461` is the template),
inject the seed, and assert cross-rank bit-identity of **all** fields including `slip1/slip2` and
`slip_rate`. RED before R-001/R-002 fixes (slip_rate diverges), GREEN after.

```cpp
// Sketch — substep leg:
//   Tpv205SubStepIterator it(...); it.SetSubSteps(deltaT, weights);
//   for step: build Q_per_node via wave.ComputeADERSubStepStates(...);
//             wave.EvaluateBulkAtFaultQPsCanonical(...);
//             it.AdvanceWithSubStepStates(dof_data, coords, Qp, Qm, dt, t,
//                                         Iplus.data(), Iminus.data(), nuc_cb);
//             wave.SetSubStepFaultImposedStates(Iplus.data(), Iminus.data(), nqp_total);
//             wave.AdvanceADER(Q, dt, 2, Q_new);
//             wave.VerifySharedFaultDOFDataConsistency(1e-10, &wr, &wf, false);
//             assert wr == 0 AND cross-rank slip_rate identical;
```

---

### [R-004] MODERATE [Makefile] — Phase-1 reconcile test left standalone, not wired into the parallel CI aggregate

**Category:** DEVIATION

**Description:**
Plan §"Testing Strategy" item 1 (and §Q2) is explicit: the injection test "**MUST be added to the
parallel CI target** (`test-parallel`-class) … **NOT only as a standalone target like the tilted
test**." After Phase 2/3 the oracle is GREEN, so the standalone-only justification (avoid breaking
CI with a RED oracle, `e2a2554` commit msg) no longer applies. The Makefile adds only the
standalone target `test-shared-fault-reconcile-cross-rank` (`Makefile:3459`); it is absent from
`test-parallel` (`Makefile:3774`) and `test` (`Makefile:2960`). The standing local guard the plan
requires therefore does not guard CI.

**Trigger:** `make test-parallel` — the reconcile test does not run; a future reconcile regression
passes CI.

**Suggested fix:**
```diff
@@ Makefile:3774 @@
-test-parallel: test-mpi-context test-parallel-utils test-parallel-domain test-mms-parallel test-parallel-fault
+test-parallel: test-mpi-context test-parallel-utils test-parallel-domain test-mms-parallel test-parallel-fault test-shared-fault-reconcile-cross-rank
```
(`test-shared-fault-reconcile-cross-rank` returns 77/skip when not np=2 or built without
`-DSEAS_TEST_INTERNAL`, so it is safe to aggregate.)

**Test case:** `make test-parallel` invokes `seas_test_shared_fault_reconcile_cross_rank` and
fails the aggregate if it returns non-zero.

---

### [R-005] MODERATE [POSSIBLE] [wave_operator.inl:ComputeADERSharedFaceFluxRHS Pass 2 + ExchangeAndPairSharedFaultQPs] — boss's *canonical-frame* `I_imp` is rotated by the non-boss's *own* `T_can`; correctness relies on can_* being bit-identical across ranks, which the code does not enforce

**Category:** BUG (POSSIBLE — needs the curvilinear run to confirm magnitude/sign)

**Description:**
The reconcile copies the boss's `I_imp_{plus,minus}_can` (boss's canonical frame). Pass 2 then
rebuilds `T_can` per QP from **this rank's own** buffered `qa.can_n/can_t1/can_t2`
(`wave_operator.inl` ~`5172`), which were taken from this rank's `qpd` (`bd.qp_data[q]`) with
**this rank's** `sign_flipped` (`wave_operator.inl:4719-4728`). The non-boss thus rotates the
boss's canonical vector with the non-boss's frame. B.3's stated goal — "Copying the boss's
`Q_imp` makes the flux that BOTH ranks assemble bit-identical" — holds only if both ranks'
canonical frames are bit-identical.

They are not guaranteed to be: `VerifySharedFaultDOFDataConsistency`'s own comment documents that
"MFEM's shared-face orientation can flip between the two ranks owning the face (observed in DIAG
output: rank 5 had `sign_flipped=1`, rank 6 `sign_flipped=0`)" (`wave_operator.inl:5996-5999`).
If the two ranks' `sign_flipped` differ, `can_*` differ by sign; even when they agree, the raw
`CalcOrtho` normal can differ at ~1e-14. Either way Pass 2 re-introduces a cross-rank difference
in the assembled global flux — exactly what the reconcile was meant to remove — re-seeding the
bulk every step. The plan's Risk Assessment listed "the canonical-frame reconstruction … must
stay bit-identical across ranks for the Pass-2 `T_can` rotation" as a known risk; the
implementation neither enforces nor asserts it.

Critically, this distinguishes the **planar 2-tet test fixture** (both ranks compute the same
axis-aligned `can_*` → frames coincide → GREEN) from the **curvilinear SAFS fault** (frames can
differ → residual survives), which is consistent with "passes local test, no effect on
production."

**Trigger:** a shared fault QP where the two owning ranks pick different `sign_flipped` (or
~1e-14-different raw normals) — common on the curvilinear SAFS fault, absent on the planar test.

**Actual behavior:** Pass-2 global fluxes are not bit-identical across ranks; bulk re-seeds.

**Expected behavior:** both ranks assemble from a single canonical frame — either reconcile the
frame too, or have the non-boss rotate with the boss's frame.

**Suggested fix (option A — reconcile the frame):** add `can_n/can_t1/can_t2` (9 doubles) to the
payload; in the callback, when adopting the peer (boss), overwrite `qa.can_n/can_t1/can_t2` with
the boss's; Pass 2 then builds `T_can` and calls `flux_.Interior(qa.can_n, …)` from the boss's
frame on both ranks → bit-identical global flux.

**Suggested fix (option B — assert the assumption, cheaper diagnostic first):** before Pass 2,
add a debug `MFEM_VERIFY` (or a `SEAS_DIAG`-gated print) that the paired peer's `can_*` match this
rank's within a tight tol; this confirms/refutes R-005 on the next curvilinear run without
changing physics. **Recommend running option B first** to decide whether option A is needed.

**Test case:** the R-003 substep leg run on a **tilted** (non-axis-aligned) 2-tet fixture (reuse
`test_rupture_tilted_fault_serial_vs_parallel.cpp`'s mesh) where the two ranks pick different
`sign_flipped`; assert cross-rank bit-identity of the assembled shared-face RHS contribution.
RED if R-005 is real, GREEN after option A.

---

### [R-006] LOW [wave_operator.inl:ExchangeAndPairSharedFaultQPs] — collective short-circuit guard is weaker than the caller's `fault_active`, risking a payload-size abort

**Category:** EDGE_CASE

**Description:**
The helper gates on `any_shared_local = (fault_dof_data_ && fault_shared_faces_.Size() > 0)`
(`wave_operator.inl` ~`5650`). The caller buffers QPs only when `fault_active = (fault_flux_ &&
fault_dof_data_ && nfs > 0 && nbf_per_face_ > 0 && fault_basis_)` (`wave_operator.inl:4570-4572`).
If `fault_dof_data_` and `fault_shared_faces_` are set but `fault_flux_`/`nbf_per_face_`/
`fault_basis_` are not, the caller produces an **empty** payload while the helper builds
`n_local_qp > 0` records, hitting the `MFEM_VERIFY((n_local_qp+1)*npay <= local_payload.size())`
abort. Aborts (loud) rather than corrupts, and the precondition is unusual in production, hence
LOW — but the two guards should match.

**Suggested fix:** make the helper's short-circuit mirror `fault_active`:
```diff
-      int any_shared_local = (fault_dof_data_ &&
-                              fault_shared_faces_.Size() > 0) ? 1 : 0;
+      int any_shared_local = (fault_flux_ && fault_dof_data_ && fault_basis_ &&
+                              nbf_per_face_ > 0 &&
+                              fault_shared_faces_.Size() > 0) ? 1 : 0;
```

**Test case:** construct a `WaveOperator` with `fault_dof_data_` + shared fault faces set but
`SetFaultFlux` not called; call `AdvanceADER` at np=2; assert no abort (helper short-circuits).

---

### [R-007] LOW [test_shared_fault_reconcile_cross_rank.cpp + Makefile comments] — "1-ULP" naming no longer matches the 1e7 Pa LSW injection

**Category:** QUALITY

**Description:**
The file header (`test_…:18-22`) and the Makefile comments (`Makefile:626-627, 2730`) describe a
"1-ULP injection hook," but the implementation injects an absolute `1.0e7` Pa for LSW
(`test_…:297`). The inline block at `:284-296` correctly explains the deviation, but the header
and Makefile still mislead a reader. Cosmetic; update the stale comments to "absolute-Pa
cross-rank seed (1e7 Pa LSW / 1e-6 Pa rate-state)."

**Suggested fix:** edit the header docstring and the two Makefile comments to drop "1-ULP."

---

## Cross-checks that PASSED (no issue found)

- **Caller↔helper QP ordering alignment:** both iterate fault shared faces in ascending `sf` ×
  ascending `q`; `fault_shared_faces_` is appended ascending (`wave_operator.inl:304-311`),
  `sf_to_basis_idx[sf] ≥ 0` for every member when `fault_active` (`:4577-4584`), and a missing
  `shared_fault_dof_offset_` entry hits `MFEM_ABORT` in the caller (`:5040`) rather than silently
  desyncing. `local_qp = local_e - my_entry_base` maps correctly to `fault_qp_buf`. No silent
  misalignment.
- **Payload pack/unpack offsets** (8 DOFData + `NUM_STATE`×2) are symmetric between build and
  callback. `REC`/`PAYLOAD_OFFSET`/`RANK_OFFSET` arithmetic is correct.
- **Boss rule** `peer_rank < my_rank_` ⇒ overwrite is unique and deterministic; `rank_a==rank_b`
  pairs are excluded, so `peer_rank == my_rank_` cannot occur.
- **Phase-3 disable switch:** `#ifdef SEAS_TEST_INTERNAL if(!disable)` wraps only the broadcast;
  Pass 2 assembly is unconditional; compiles out in production. The switch is effective in the
  test because `wave_operator.cpp` is empty and the template is instantiated per-TU (the test TU
  is built `-DSEAS_TEST_INTERNAL`).
- **Collective safety:** `MPI_Allreduce(any_shared)` short-circuit + `MPI_Allgatherv`; ranks with
  no shared fault faces contribute empty buffers. Serial returns 0.
- **`tau*_nuc` not reconciled** is acceptable: they are *inputs* (per-rank gradual_overstress
  accumulators) folded into `tau*_corr` (which IS reconciled) and `I_imp` (reconciled); they are
  not in the verify set. (Worth a one-line note in the plan, but not a bug.)

## Summary
- Critical issues: 3 (R-001 slip_rate omission; R-002 oracle blind to slip_rate; R-003 test
  bypasses the substep iterator)
- Moderate issues: 2 (R-004 not in `test-parallel`; R-005 POSSIBLE Pass-2 per-rank `T_can` frame
  mismatch on curvilinear faults)
- Low issues: 2 (R-006 guard asymmetry; R-007 stale "1-ULP" naming)
- Plan compliance: **PARTIAL** — Phase 2 reconcile + Phase 1/3 tests are implemented and the
  reconcile runs in production, but (a) the payload/verify omit `slip_rate` (the production
  observable), (b) the oracle never exercises the production substep-iterator path, (c) the
  required `test-parallel` wiring is missing, and (d) the plan's own "can_* bit-identical across
  ranks" prerequisite for Pass-2 is not enforced.
- **Verdict: FAIL — must fix R-001, R-002, R-003 before the reconcile can affect the SAFS run.**
  R-001+R-002 are the direct cause of "passes local test, no effect on production." R-005 should
  be confirmed/refuted with a cheap assert (option B) on the next curvilinear run, since it may be
  a *second* reason the bulk still diverges even after R-001/R-002 are fixed.

## Unreviewed Areas
- `FaultBasis::AppendSharedFaces` / `ComputeQPBasisShared` — not opened; R-005 hinges on whether
  these produce bit-identical `can_*` across the two owning ranks (the plan asserts they should
  post-Part-A; the code path was not traced here). Recommend the option-B assert to settle it.
- The Phase-0 `SEAS_DIAG_XRANK` / RAWDOF trace (`6a3af3e`) was reviewed only for env-gating and
  byte-exactness-when-off (confirmed); the trace's numerical content was not independently
  recomputed.
- Frontera-only behavior (np=400 deadlock, actual `worst_rel`/`V_max` post-fix) cannot be checked
  here; the analysis is from the source + the local-vs-production call-graph divergence.
