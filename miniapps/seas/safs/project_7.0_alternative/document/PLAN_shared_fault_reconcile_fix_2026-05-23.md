# Implementation Plan: shared-fault cross-rank DOFData consistency (SAFS)

## TL;DR — there are TWO separate problems; this plan is about the second one

| | Problem A — the blow-up | Problem B — the cross-rank inconsistency |
|---|---|---|
| **What** | The fault +/- side was mislabeled on the tilted (curvilinear) fault, so the nucleation pushed with the wrong sign → V_max → 10³ blow-up. | The same fault quad-point is solved **twice** (once per MPI rank) from inputs that differ at the 16th digit; at the rupture onset that hair's-width difference makes the two ranks disagree (slip vs lock). |
| **Scope** | Geometry (which side is "+"). Independent of friction law / dip / strike. | MPI / floating-point. Triggered the same way for any channel; only **runs away** when there is dip slip. |
| **Status** | **FIXED** — commit `800b3281` (see Part A). Validated by the `zerodip` job. | **OPEN** — this plan. |
| **Test** | `spatial_dyn_zerodip_8N_400r_dev_2hr_safs.sbatch` (no dip slip, guard non-fatal). | New local injection test + the Dc2 `XRANK` Frontera run. |

Documents behind this plan: diagnosis `spatial_dynamic_rupture_speckle_blowup_2026-05-22.md`,
design review `spatial_dynamic_rupture_reconcile_review_2026-05-23.md`.

---

## Part A — the blow-up fix (DONE, commit `800b3281`) — documented for the record

> **This part is already landed and is NOT implemented by this plan.** It is
> documented here because the user asked to keep it separate and because it is
> the reason the `zerodip` run is now stable.

### What was wrong (plain language)
On each fault face the code must decide which of the two touching elements is on
the "+" side and which is on the "−" side (this sets the sign of the nucleation
push and the flux). It decided by measuring which element sits "below" a FIXED
reference direction `ref_normal = (0,−1,0)` — i.e. it projected the
element-to-face offset onto that hardcoded direction. That works for a flat
fault whose normal IS (0,±1,0). The meshed SAFS fault is **curvilinear**: where
the fault swings to a roughly N–S strike, its normal becomes nearly
perpendicular to `(0,−1,0)`, so the projection measures mostly the *sideways*
part of the offset instead of the real "above/below." Past a tilt of ≈71.6° the
sideways part wins and the side label **flips** — both ranks of a shared face
then call themselves "−", the "exactly one + side" rule breaks, the nucleation
is applied with the wrong sign, and the rupture blows up.

### What we changed (commit `800b3281`)
Project the element-to-face offset onto the **actual face normal** (canonicalized
to a rank-independent direction), not the hardcoded `ref_normal`:
- `fault/fault_basis.hpp` — new shared static
  `FaultBasis::NormalNeedsFlipToCanonical(n_raw, ref_normal, dim)` (orient the
  normal LINE to `n·ref_normal>0`, with a largest-component fallback when the
  fault normal is ⟂ `ref_normal`); `ComputeOrientedFrame` routes its sign
  decision through it.
- `dynamic/wave_operator.inl` — both side-determination blocks (interior
  `:459-465`, shared `:519-525`) compute the canonical face normal via
  `CalcOrtho` and project the centroid offset onto it instead of `ref_normal`.
- For a flat fault the canonical normal **equals** `ref_normal`, so TPV/BP5 are
  byte-exact; the change only bites where the fault tilts.

### How Part A is validated (per the user's directive)
Job script: `jobs/safs/spatial_dyn_zerodip_8N_400r_dev_2hr_safs.sbatch`.
- `SEAS_ZERO_DIP_PRESTRESS=1` (default) — zeroes the dip pre-stress at the source
  so there is **no dip slip**; the fault is loaded in strike only.
- `SEAS_R101_NONFATAL=1` — the cross-rank guard logs `worst_rel` each check
  instead of aborting, so **Problem B does not force-abort the Problem-A run**
  (this is the "do not force guard the MPI DOFData inconsistency" requirement;
  the mechanism already exists, commit `fc28454`, and the call site is correct:
  `spatial_dyn_driver.cpp:1994` passes `tol, &worst_rel, &worst_field,
  abort_on_fail=!nonfatal`).
- Expected: the run advances with no V_max blow-up (Part A fixed). The R-101
  `worst_rel` stays bounded because with no dip slip there is nothing for the
  cross-rank inconsistency to run away into (see "why dip, not strike" below).

**Acceptance for Part A (already met / re-confirm):** the zerodip job runs past
the old blow-up window with `V_max` physical (no 10³) and the non-fatal R-101
log shows `worst_rel` bounded (not → O(1)).

---

## Part B — the cross-rank DOFData inconsistency (THIS PLAN)

### B.1 The problem in plain language

A fault quadrature point that lies on the boundary between two MPI ranks is a
**shared** point: rank A owns the element on one side, rank B owns the element
on the other. To advance the fault friction there, each rank needs the stress on
**both** sides. Each rank already has its own side exactly, and gets the other
side as a **ghost copy** from its neighbor. So both ranks have "both sides" and
both run the friction solve for that point. **The same physics is computed
twice, once on each rank.** (This is what R-701 set up: it deleted an older
"one rank computes, tells the other" broadcast, betting that both ranks would
get identical answers anyway.)

That bet is wrong by a hair. The two ranks do **not** feed the friction solve
bit-identical numbers. The reason is pure floating-point bookkeeping, not a
logic bug:

- The shared point is one physical location, but rank A reaches it through its
  element's local coordinates and rank B reaches the *same* point through the
  *neighbor* element's local coordinates. The interpolation weights
  (`shape × DOFs`) are mathematically equal but rounded differently, so the bulk
  stress the two ranks plug in differs at about the **16th significant digit**
  (~1e-14). (Confirmed on Frontera: the two ranks' inputs match to 11 printed
  digits and differ below that.)

Normally a 1e-14 difference in → a 1e-14 difference out: invisible. The trouble
is the **rupture onset is a yes/no switch**. The slip rate is
`V = max(0, (|τ| − strength)/η_s)`: below the strength the point is **locked**
(`V=0`); above it, it **slips**. Right at the threshold `|τ| = strength`, a
hair's-width difference in `|τ|` flips the answer between "locked" and "slips."
The SAFS nucleation ramps the stress up **slowly**, so the point sits **right on
the knife's edge for several sub-steps** — and during that window the 1e-14
difference makes **rank A slip while rank B stays locked** (or vice-versa). From
that instant the two ranks are computing two different earthquakes at that point;
the slip-weakening feedback then drives them apart, and the run is corrupted.

A runtime guard (the R-101 check) compares the two ranks each step and aborts
when they disagree — that is what stops the SAFS run at t≈0.455 s.

**Why benchmarks never caught it:** TPV102/104 use rate-and-state friction,
whose slip rate is a **smooth** function (no yes/no switch) — a 1e-14 input stays
a 1e-14 output forever, so the ranks never visibly disagree. TPV205 uses the
same slip-weakening switch as SAFS, but its rupture front sweeps past a point in
one step (no slow knife's-edge dwell). Only SAFS combines the slip-weakening
switch **and** a slow nucleation that parks a shared point on the threshold.

### B.2 Why the DIP direction blows up but the STRIKE direction is fine

This is the part that confused us, so here it is carefully.

The **trigger** above is **blind to direction**: the yes/no slip switch is on the
*total* shear magnitude `|τ|`, so the slip-vs-lock split corrupts the **whole**
slip vector — both the strike component and the dip component become inconsistent
across ranks at the same instant. (In the Frontera trace, when one rank slipped
it had both `V2` (strike) and `V1` (dip) nonzero while the other rank had both
zero.)

What differs is the **consequence** — whether that inconsistency stays small or
runs away — and that **is** direction-dependent, for a physical reason:

- **Strike slip does not change the fault's normal stress.** Sliding along strike
  is motion *in the fault plane, along the surface*; it doesn't push the two
  sides together or pull them apart. So a strike inconsistency is just a small
  error on a large, well-loaded quantity — it radiates away and stays bounded.
- **Dip slip on the ~32°-dipping SAFS fault DOES change the normal stress**
  (slip ↔ normal-stress coupling: dip motion on a dipping fault has a component
  that compresses/decompresses the fault). And the fault's **normal-stress
  channel is undamped** (unlike the shear channels, which radiation damping
  `η_s·V` holds in check). So a dip inconsistency feeds the *un-damped* normal
  channel: more dip slip → changed σ_n → changed strength → more slip → … a
  positive feedback loop with no brake. That is the runaway that becomes the
  blow-up.

So: **strike inconsistency = small error on a stable channel (bounded); dip
inconsistency = small error injected into an undamped feedback loop (runs away).**

This is exactly what the **zerodip experiment proves empirically**: zero the dip
pre-stress (no dip slip) and the *same code with the same 1e-14 split* runs
without blowing up — the strike-only inconsistency stays bounded (the R-101
non-fatal log shows `worst_rel` small). Restore the dip load and it runs away.
(The zerodip sbatch header states this hypothesis directly: "that dip-slip drive
(via slip↔normal-stress coupling) is what makes the rupture blow up.")

**Important:** the fix below removes the inconsistency **at its source for every
direction**, so we do NOT rely on strike being self-stabilizing. After the fix,
strike and dip are both cross-rank-identical.

### B.3 How the proposed fix works in plain language

Stop computing the shared point twice and hoping the two answers match. Instead,
**pick one rank to be the boss for that point; the boss computes the friction
solve; it tells the other rank the answer; both ranks use the boss's answer.**
Now there is only ONE answer, so the two ranks cannot disagree — the yes/no slip
switch is decided once, by the boss, and copied. The 1e-14 input difference is
irrelevant because the non-boss never uses its own solve.

Concretely (still plain):
1. Both ranks do their friction solve as today (cheap; we don't restructure the
   solver).
2. They **exchange** the results for every shared fault point (a small MPI
   message — only the points on rank boundaries).
3. For each point, the **non-boss overwrites its own answer with the boss's** —
   both the friction state (slip rate, corrected traction, …) **and** the
   "imposed stress" that gets pushed back into the bulk wave solve. (We must copy
   the imposed stress too, not just the stored numbers — otherwise the bulk would
   re-grow a fresh 1e-14 difference next step and we'd be back where we started.)
4. THEN both ranks assemble their side of the flux from that single shared
   answer.

"Boss" = the rank with the **smaller MPI rank number** of the two sharing the
point. It's a fixed, unique, geometry-independent choice (so it can't be confused
by the curvilinear fault the way the +/- side label was).

**Why this is the right kind of fix (general, not a patch):** it does the same
thing for every friction law and every slip direction — it just makes the shared
point single-valued. It does not tune a threshold, special-case the dip channel,
or lean on strike being lucky. (An older "owner-broadcast" doing exactly this
existed and was deleted by R-701; we are restoring it in the path R-701 left
uncovered, now that we've proven R-701's "both ranks agree anyway" assumption is
false.)

**One honest consequence (needs your OK):** making the non-boss adopt the boss's
value changes TPV102/104 by ~1e-14 (they currently differ across ranks by that
much; after the fix they're bit-identical across ranks). So we trade
"bit-identical to the old binary" for "physically identical (≤1e-13) AND now
bit-identical across ranks" — a stronger correctness property, but it does
relax the strict byte-exact rule. Flagged again here.

---

## Constraints
- **Method- and direction-invariant:** the reconcile runs for LSW,
  LSW_ForcedRupture, AND rate-state; no per-law / per-direction / threshold gate.
- **Regression contract (relaxed, see B.3):** TPV102/104/205 + BP5
  physically-exact (`worst_rel ≤ 1e-13` vs pre-fix) AND cross-rank bit-identical
  (R-101 `max_rel_diff == 0`). Overrides CLAUDE.md byte-exact; recorded; needs
  user sign-off.
- **MPI-collective-safe** (R-1600 deadlock class): every rank reaches the
  exchange; ranks with no shared fault faces contribute empty buffers and do not
  hang.
- **Boss = lower MPI rank** (unique, geometry-independent; `MFEM_VERIFY` exactly
  one boss per shared point).
- **Files Requiring Extreme Care:** `wave_operator.inl/.hpp`, `fault_face_flux.cpp`.
- Do not touch the friction solvers; do not add hardcoded constants.

## Phase 1: reproducible local test (injection-based, np=2, both laws)

### Goal
A deterministic np=2 unit test that reproduces the cross-rank split by injecting
a 1-ULP input difference at the slip-onset threshold and asserts the two ranks'
`DOFData` are bit-identical — RED before Phase 2, GREEN after — for BOTH
rate-state and LSW (proving the fix is method-invariant).

### Files to Create
- `tests/unit/test_shared_fault_reconcile_cross_rank.cpp`.

### Files to Modify
- `Makefile` — 5 entries mirroring the tilted-test wiring (SRC ~323, OBJ ~623,
  link target `seas_test_shared_fault_reconcile_cross_rank` with the tilted
  test's object set, compile rule ~2719, run target
  `test-shared-fault-reconcile-cross-rank` at np=2). Not in `make test`.
- `dynamic/fault_face_flux.hpp` — behind the existing `#ifdef SEAS_TEST_INTERNAL`:
  `static real_t s_seas_test_qplus_xz_perturb_ulp = 0.0;` consumed in
  `ComputeTrialTraction` to add `s · ULP · Q_plus[SXZ]` ONLY when set (compiles
  out in production). (Or, if it reproduces, inject by perturbing one rank's
  bulk Q DOF in the test — preferred; decide at implementation.)

### Detailed Requirements
1. Fixture: reuse `BuildTwoTetFaultMesh()` (planar y=0, 2 tets) + the
   `cy<0→rank0/cy≥0→rank1` partition from
   `test_rupture_multistep_serial_vs_parallel.cpp`. One shared fault face, 3 QPs.
2. Drive `|τ|` across `μ_s·σ_n` during the run (derive load from `TPV205Params` /
   `TPV102Params`); assert at least one step has `0 < V_abs < V_small` on one
   rank (the onset is reached), else FAIL "did not cross slip-onset threshold."
3. On rank 0 only, perturb its `Q_plus[SXZ]` by exactly 1 ULP (`std::nextafter`)
   vs rank 1's matching input — the deterministic stand-in for the 1e-14
   interpolation seed.
4. Run `kNSteps` ADER-2 steps through the threshold.
5. Each step gather both ranks' shared-QP `DOFData` via
   `VerifySharedFaultDOFDataConsistency(tol,&wr,&wf,/*abort=*/false)` (added in
   `a9bd4d2`); assert `wr == 0`.
6. Loop over `{RateState, LSW}` in one executable; the assertion is per-law.

### Acceptance Criteria
- [ ] Builds + runs at np=2.
- [ ] Before Phase 2: RED — LSW `wr→O(1)`; rate-state `wr~1e-14`>0 (shows the
      defect is method-invariant, just sub-threshold for rate-state).
- [ ] After Phase 2: GREEN — `wr == 0` for BOTH laws, every step.
- [ ] Production build (`SEAS_TEST_INTERNAL` off): hook compiles out.

### Dependencies
- Depends on: nothing. Required by: Phase 2 (oracle), Phase 4.

## Phase 2: the reconcile (the fix)

### Goal
After the friction solve on a shared fault face, both ranks hold bit-identical
`DOFData` and assemble from the boss's imposed stress — for ALL friction laws —
so the slip/lock decision is made once and copied.

### Files to Modify
- `dynamic/wave_operator.inl` `ComputeADERSharedFaceFluxRHS` (~4671-4797): make
  the shared-fault QP handling **two-pass** with a reconcile between them.
- `dynamic/wave_operator.hpp`: declare the helper + assembly-buffer struct.
- Factor the R-101 record gather/pair (`:5470+`) into a shared helper used by
  both the verify and the reconcile (one matcher).

### Detailed Requirements
1. Helper (parallel-only):
   ```cpp
   int WaveOperator<MeshType>::ExchangeAndPairSharedFaultQPs(
        int npay, const std::vector<double>& local_payload,
        const std::function<void(int local_qp, const double* peer_payload,
                                 int peer_rank)>& cb) const;
   ```
   Gathers one record/local-shared-QP (face-vertex-key + qp_idx + rank + payload)
   via the existing `MPI_Allgatherv`, pairs by (face-key, qp_idx), invokes `cb`
   for each paired local QP. The verify uses `npay=8` + a compare callback; the
   reconcile uses `npay = 8 + 2·NUM_STATE` (8 DOFData fields + `I_imp_plus_can` +
   `I_imp_minus_can`) + an overwrite callback.
2. Per-QP assembly buffer:
   ```cpp
   struct SharedFaultQPAssembly {
      int dof_offset1, ndof, dof_idx; bool elem1_on_plus; real_t w;
      real_t can_n[3], can_t1[3], can_t2[3];
      real_t shape1[MAX_NDOF];
      real_t I_imp_plus_can[NUM_STATE], I_imp_minus_can[NUM_STATE];
   };
   ```
3. Pass 1 — compute (no assembly): run the existing friction dispatch (LSW /
   LSW_ForcedRupture / rate-state) → fill `fault_dof_data_[dof_idx]` and
   `I_imp_*`; store the buffer; do NOT touch `rhs`.
4. Reconcile: payload = 8 DOFData fields + `I_imp_plus_can` + `I_imp_minus_can`.
   Call the helper; in the callback, if `peer_rank < my_rank_` (peer is boss),
   overwrite this rank's `fault_dof_data_[dof_idx]` AND the buffer's
   `I_imp_*_can` with the peer's payload. `MFEM_VERIFY` every local shared QP was
   paired.
5. Pass 2 — assemble from the (possibly overwritten) `I_imp_*_can`, rotated to
   global via `T_can` (rebuilt from buffered `can_*`, bit-identical on both
   ranks), exactly as the current inline assembly.
6. No friction-law gate (method-invariant). Rate-state `else` branch included.
7. Collective safety: the helper uses the same `MPI_Allreduce(any_shared)`
   short-circuit + `MPI_Allgatherv` as the verify; empty contribution from ranks
   with no shared fault QPs.

### Edge Cases to Handle
- Unpaired shared QP → `MFEM_VERIFY` fail (same policy as the verify).
- Same-`elem1_on_plus`-on-both-ranks edge (θ≈90°) → boss-by-rank still unique.
- Serial / np=1 → no shared faces → two-pass degenerates to the current path
  (byte-exact).

### Acceptance Criteria
- [ ] Phase-1 test GREEN for both laws (`wr == 0` every step).
- [ ] TPV102/104/205 + BP5: `worst_rel ≤ 1e-13` vs pre-fix AND R-101
      `max_rel_diff == 0`.
- [ ] Local fault/TPV suite green (fault-basis trio, shared-fault role/dof-data,
      interior-flux ×3, godunov identity, tpv102 locked/absorbing/pepper/
      ader-smoke, multistep + tilted serial-vs-parallel, sign-flipped truth
      table).
- [ ] np=2,4,8 incl. a partition with some rank lacking shared fault faces — no
      deadlock.

### Dependencies
- Depends on: Phase 1. Required by: Phase 4.

## Phase 3: correctness guard (runtime, method-invariant)

### Goal
Keep `VerifySharedFaultDOFDataConsistency` as the permanent runtime sentinel so
any future regression aborts loudly; post-fix it must report `max_rel_diff == 0`.

### Files to Modify
- `dynamic/wave_operator.inl`: refactor the verify onto the shared
  `ExchangeAndPairSharedFaultQPs` helper (one matcher); update the stale R-501
  wording in the abort message to point at the Phase-2 reconcile.
- Keep `SEAS_R101_NONFATAL` (commit `fc28454`) as the diagnostic escape hatch;
  keep the existing gated cadence; add `SEAS_VERIFY_XRANK_EVERY=N` (default
  unchanged) for debugging.

### Detailed Requirements
1. The guard stays method-invariant (checks all shared QPs regardless of law).
2. Document that post-fix the expected `max_rel_diff` is `0`; any nonzero is a
   real regression, not roundoff.
3. Negative test (Phase 1, with the reconcile disabled via a test switch): the
   guard trips — proving it still catches a real desync.

### Acceptance Criteria
- [ ] Post-fix Dc2 `XRANK` Frontera run: guard reports `max_rel_diff == 0`
      through and past t≈0.455 s; no abort; run advances past t=1.0 s.
- [ ] Negative test trips the guard.

### Dependencies
- Depends on: Phase 2. Required by: Phase 4.

## Phase 4: full regression + Frontera re-validation

### Goal
End-to-end: the Dc2 mesh advances past t=1.0 s (no abort, no blow-up); fix is
config-agnostic.

### Detailed Requirements
1. Local regression battery (Phase-2 AC) green; document unchanged pre-existing
   failures (adjacent-triangle pepper-bug, r101 missing-precondition, macOS MUMPS
   Bus error).
2. Re-run the `SEAS_DIAG_XRANK` Dc2 sbatch: at the onset QP both ranks now show
   identical V1/V2/τ*_corr through t=0.478 s; R-101 passes; reaches tfinal.
3. Re-run the physical `D_c=1.0` baseline config (not just the inflated Dc2) to
   confirm the fix is not tuned to Dc2.
4. Re-confirm Part A: the `zerodip` job still runs clean (the reconcile must not
   regress the strike-only path).

### Acceptance Criteria
- [ ] Local regression green (modulo documented pre-existing).
- [ ] Frontera Dc2 AND `D_c=1.0`: no R-101 abort, no blow-up; `V_max` peaks then
      decreases.
- [ ] zerodip job unchanged (still clean).

### Dependencies
- Depends on: Phase 2, Phase 3.

## Testing Strategy
- **Local oracle:** Phase-1 injection test (np=2, both laws), RED→GREEN — proves
  the reconcile absorbs the seed for any friction law.
- **Guard:** Phase-3 verify, post-fix `max_rel_diff==0`; negative test trips it.
- **Regression:** TPV102/104/205 + BP5 — capture pre-fix, assert `≤1e-13` AND
  R-101 `max_rel_diff==0`; stash-rebuild-baseline to confirm the ~1e-14 TPV change
  is the reconcile (non-boss→boss), not a physics change.
- **Part-A non-regression:** the zerodip job stays clean.
- **MPI robustness:** np=2/4/8, incl. a no-shared-fault-on-some-ranks partition.

## Risk Assessment
- **Deadlock (R-1600 class):** the new exchange must be reached by every rank;
  detect via np=4/8 smoke hangs.
- **Per-sub-step exchange cost:** correctness-first; if hot on 400r, move the
  helper to point-to-point over the face-nbr topology (same payload). NO
  kink-proximity gate (re-introduces problem-specificity).
- **Two-pass refactor indexing slip** (dof_offset1/shape1/w): the interior path
  is untouched; tilted + multistep + adjacent-triangle tests cover shared-face
  assembly; assert Pass-2 reproduces the pre-refactor RHS when the reconcile is a
  no-op.
- **Byte-exact relaxation:** the ~1e-14 TPV change overrides a CLAUDE.md
  non-negotiable; recorded; needs sign-off.
- **Tricky existing code:** `ComputeADERSharedFaceFluxRHS` (4671-4797, the
  R-1601 inline fallback + assembly), `VerifySharedFaultDOFDataConsistency`
  (5364+, the matcher being factored), the canonical-frame reconstruction
  (4695-4709, must stay bit-identical across ranks for the Pass-2 `T_can`
  rotation).
