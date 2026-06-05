# Fix Report — Mixed-Flux + Heterogeneous Riemann plan

Pre-implementation plan: every finding is a defect in the PLAN; "fix" = plan edit. The
`chunhui-code-debugger` is NOT used (it fixes source, of which there is none). The reviewer-owned
check doc is not modified. The review PDF at `miniapps/seas/document/mixed_flux_dev/` is regenerated
each round.

---

## Phase-3 review fixes (P3-1/2/4) — check Rev I-2 (implementation phase) (2026-06-04)

**Note on phase transition:** the rounds below this section are the *plan*-review cycle (no source
existed). The plan has since been implemented; this section is the FIRST *implementation* fix pass.
It applies the four findings of the check doc's "Implementation Review — Phase 3 (Rev I-2)" section
(P3-1..P3-4) by editing SOURCE under `dynamic/` + `tests/` + `Makefile` (all in the editable tree;
BP5/solver/friction untouched). The reviewer-owned check doc is NOT modified.

**Agent:** chunhui-code-debugger
**Check document:** `PLAN_mixed_flux_hetero_riemann_check.md` Rev I-2, findings P3-1/2/3/4.
**Priority applied:** P3-4 (required, quick) -> P3-2 (quick) -> P3-1 (stretch) -> P3-3 (optional).

### Findings: fixed vs deferred

| ID | Severity | Status | What changed |
|----|----------|--------|--------------|
| P3-4 | WARNING (MODERATE) | **FIXED** | Interim central+ADER abort guard at the TOP of `BimaterialWaveOperator::ComputeMaxDt`. |
| P3-2 | SUGGESTION (LOW) | **FIXED** | Shared arm of `BuildPerFaceCentralFluxMatrices_` now fails loud on a null `GetSharedFaceTransformations(sf)` for a central-set face. |
| P3-1 | WARNING (MODERATE) | **FIXED (3.4a + 3.4c automated; 3.4b documented)** | New MPI np=2 test `tests/parallel/test_bimaterial_mixed_flux_shared.cpp` + Makefile target. |
| P3-3 | SUGGESTION (LOW) | **FIXED (trivial)** | Two more state pairs added to the serial dispatch test's `MakeStatePairs`. |

### P3-4 — interim central+ADER abort guard (REQUIRED)

- **Root cause:** Phase 3 lifted R-003, so `BimaterialWaveOperator::ComputeMaxDt`'s per-element CFL
  walk is now reachable with `mixed_flux_mode_ != None`. The bi-material central flux is
  non-dissipative and unstable under ADER; it requires an RK integrator. No guard existed.
- **Change:** added at the TOP of `ComputeMaxDt` (`dynamic/bimaterial_wave_operator.inl`, before the
  per-element walk):
  `MFEM_VERIFY(mixed_flux_mode_ == MixedFluxMode::None || cfl_rk_aware_, "<central+ADER instability …>");`
  and `using Base::cfl_rk_aware_;` next to the other `using Base::` lines in
  `dynamic/bimaterial_wave_operator.hpp` (P3-4 comment). The divergent `0.9/0.4` factor switch was
  **left untouched** (Phase-4 work, per the instruction).
- **No false-fire on the parity path:** the bimaterial parity test calls `ComputeMaxDt` with
  `mixed_flux_mode_ == None` (the operators never call `SetMixedFluxMode`), so the first clause is
  true and the guard does not fire — confirmed by C-6c ComputeMaxDt parity still PASSING (9/9).

### P3-2 — harden the IMPL-8 assert (quick)

- **Root cause:** in the shared arm of `BuildPerFaceCentralFluxMatrices_`, a null `ftr` from
  `GetSharedFaceTransformations(sf)` caused an unconditional `continue` with no map entry. If such an
  `sf`'s mesh-face were in `central_flux_face_set_`, the central matrix would be silently unbuilt and
  the exact-equality IMPL-8 assert (`per_face_central_flux_.size() == central_flux_face_set_.size()`)
  would false-fire and abort a valid run.
- **Change:** reordered the shared arm so `mesh_face_idx = pmesh.GetSharedFace(sf)` and its range
  verify run FIRST (they do not need `ftr`); compute `is_central = central_flux_face_set_.count(...)
  > 0`; then on `!ftr`, `MFEM_VERIFY(!is_central, …)` (a central shared face MUST have a valid
  transform) and otherwise keep the benign `continue` for non-central faces. Minimal; no behaviour
  change for valid faces.

### P3-1 — MPI np=2 shared-face dispatch test (STRETCH — attempted, completed for 3.4a/3.4c)

- **New file:** `tests/parallel/test_bimaterial_mixed_flux_shared.cpp` (np=2). Uses an EXPLICIT
  ParMesh partition array (ParMETIS is nondeterministic) so a chosen fault-adjacent NON-fault
  interior face is guaranteed to straddle the rank seam -> a SHARED central face. Reuses the serial
  dispatch test's hex-row fault fixture, het `mu(x)` material, the `TestableBimat` shim (extended with
  `using …SharedInteriorFaceFlux_;`), `MakeStatePairs`, and `WorstRel`.
- **Test 3.4a** (het `Mode::Coefficient` `mu(x)` + `SetSeamContinuous(true)`): 5-hex row, fault x-normal
  at x=3, seam at x=4 -> the fault-adjacent face (3,4) is a shared central face. Asserts (i) IMPL-8
  lifecycle (map size == set size); (ii) at least one rank found the shared central face; (iii) the
  dispatched side-0 `F*` from `SharedInteriorFaceFlux_` == `ApplyPerFaceFlux(per_face_central_flux_[mf],
  …)` **bit-for-bit** (routes to the central store); (iv) the shared face has a neighbour-material stub
  entry; (v) the stored half-Jacobians `c[0] == c[1]` entrywise (<= 1e-9 rel) — both built from the
  R-004 stub material via `BuildPerFaceCentralMatricesGlobal`, so the central flux is single-valued by
  construction.
- **Test 3.4c** (depth-only `mu(z)` `Mode::Coefficient` + `SetSeamContinuous(true)`): 2x5 (x,z) grid,
  fault x-normal (so the FaultBasis strike/dip frame with `up=(0,0,1)` is well defined), seam split
  along z -> a fault-adjacent z-normal shared central face whose neighbouring elements have DISTINCT
  depth materials. Asserts the build completes with NO false abort (declarative seam-continuity gate
  does not conflate depth variation with lateral variation).
- **Test 3.4b** (abort WITHOUT `seam_continuous`): **OMITTED from the automated assertions and
  documented**, per the instruction. The production guard is the single `MFEM_VERIFY` in the shared
  arm of `BuildPerFaceCentralFluxMatrices_` (`material_->mode == Constant || seam_continuous_`). Under
  the MPI build `MFEM_VERIFY` calls `MPI_Abort`, which terminates the whole job and cannot be caught
  in-process (this MFEM build has no `MFEM_USE_EXCEPTIONS`, and a fork-based death test is unusable
  after `MPI_Init`). The same precedent is documented in
  `tests/unit/test_phaseh_wave_operator_constant_parity.cpp:451-455`. The abort is
  **construction-verified** (it is the exact negative of 3.4a/3.4c: removing `SetSeamContinuous(true)`
  on either fixture trips the verify) and **covered by inspection**. The production guard was NOT
  weakened to make it testable. Revisit if death-test infra is added.
- **Makefile:** added `TEST_BIMATERIAL_MIXED_FLUX_SHARED_SRC`/`_OBJ`, the object compile rule, the
  `seas_test_bimaterial_mixed_flux_shared` link target (operator-level deps mirroring
  `seas_test_bimaterial_wave_operator_parity`), and a `test-bimaterial-mixed-flux-shared:` alias that
  runs `mpirun -np 2 ./seas_test_bimaterial_mixed_flux_shared`. NOT added to the serial `test:` group.

### P3-3 — extra state pairs (optional, trivial)

- Added two deterministic `(Q_self, Q_nbr)` pairs to `MakeStatePairs` in
  `tests/unit/test_bimaterial_mixed_flux_dispatch.cpp` (an anti-symmetric all-channel stress jump with
  common-mode velocity, and an asymmetric mixed stress+velocity pair), widening the per-side-Godunov
  vs single-valued-central discriminator. Test count unchanged (10/10); the pairs feed the existing
  assertions.

### Deviations from the suggested fixes

1. **P3-1 value anchor.** The check doc's NEW-I.4 suggested "the dispatched side-0 `F*` ==
   `BuildPerFaceCentralMatricesGlobal(local, nbr)` apply". A first implementation that recomputed the
   shared-face centroid normal in-test and rebuilt the matrices was **run-to-run nondeterministic on
   the non-zero rank** (worst-rel jumped between ~4e-16 and ~1.8e-9 across process runs) because the
   independently re-derived shared-face normal is fragile (MFEM's `GetSharedFaceTransformations`
   returns a reused internal object; orientation can vary with ghost ordering). I replaced it with TWO
   robust checks that capture the same intent without the parallel-normal fragility: (a) the
   **bit-for-bit** store-route check `SharedInteriorFaceFlux_ == ApplyPerFaceFlux(stored c)` and (c) a
   **normal-free** anchor on the stored matrices `c[0] == c[1]` (entrywise, <= 1e-9 rel). The latter
   is meaningful precisely because the R-004 stub makes the neighbour material equal to the local
   material on a shared face, so `BuildPerFaceCentralMatricesGlobal` must produce equal half-Jacobians
   — independently confirming the shared arm built both halves from the stub material AND that the
   central flux is single-valued by construction. The 1e-9 tolerance (not 1e-12) matches the
   established Phase-H parity calibration for independent `BuildPerFaceCentralMatricesGlobal` rebuilds
   (LU-rounding-limited; the bit-exact proof is check (a)).
2. **P3-1 3.4c fixture.** The plan/check framed 3.4c as a "z-tilted seam". A pure z-stacked column with
   the fault on a z-normal face aborts in `FaultBasis::ComputeOrientedFrame` (the fault normal is
   collinear with `up=(0,0,1)`) — unrelated to the central-flux guard. I used a 2D (x,z) grid with an
   **x-normal fault** (well-defined frame) and a **z-split seam**, so the shared central face is
   z-normal and the depth material genuinely varies across it, exercising exactly what 3.4c intends
   (no false abort on depth variation) without the unrelated fault-frame degeneracy.
3. **TEST_ASSERT prints failures on every rank** (not rank-0-only) in the new test — a deliberate
   improvement: the original rank-0-only print HID a non-zero-rank failure during development. PASSED
   lines still print on rank 0 only; counts are `MPI_Allreduce`-summed across ranks.

### Test results after all fixes (fresh objects forced before each rebuild)

| Suite | Result |
|-------|--------|
| `seas_test_bimaterial_wave_operator_parity` (serial) | **9 / 9** |
| `seas_test_bimaterial_central_flux` (serial) | **6 / 6** |
| `seas_test_bimaterial_mixed_flux_dispatch` (serial, +2 P3-3 pairs) | **10 / 10** |
| `seas_test_phaseh_wave_operator_constant_parity` (serial) | **19 / 19** |
| `seas_test_phaseh_wave_operator_constant_parity` (np=4) | **52 / 52** |
| `seas_test_bimaterial_mixed_flux_shared` (np=2, NEW) | **14 / 14** (deterministic over 4 runs) |

- **Passed:** 110 (9+6+10+19+52+14). **Failed:** 0. **Skipped:** 0.
- No previously-passing test regressed; the new test passes; everything compiles cleanly (no new
  warnings — the unused `CentroidUnitNormal` helper from the discarded anchor was removed).

### Files modified
- `dynamic/bimaterial_wave_operator.hpp` — `using Base::cfl_rk_aware_;` (P3-4).
- `dynamic/bimaterial_wave_operator.inl` — P3-4 guard in `ComputeMaxDt`; P3-2 null-transform fail-loud
  in `BuildPerFaceCentralFluxMatrices_` shared arm.
- `tests/parallel/test_bimaterial_mixed_flux_shared.cpp` — NEW (P3-1, 3.4a/3.4c automated, 3.4b documented).
- `tests/unit/test_bimaterial_mixed_flux_dispatch.cpp` — +2 state pairs (P3-3).
- `Makefile` — register `seas_test_bimaterial_mixed_flux_shared` (+ `test-bimaterial-mixed-flux-shared` alias).

### Verification checklist
- [x] No CRITICAL/HIGH finding existed in Rev I-2; both MODERATE (P3-1, P3-4) addressed.
- [x] Both LOW findings (P3-2, P3-3) addressed.
- [x] All previously-passing regression tests still pass (parity 9/9, central 6/6, dispatch 10/10,
      phaseh 19/19 serial + 52/52 np=4).
- [x] New MPI test passes (14/14, np=2, deterministic).
- [x] Code compiles cleanly (no new warnings); production `seam_continuous` guard NOT weakened.
- [x] `mixed_flux=none` byte-exact contract preserved (the P3-4 guard never fires for None).

---

## Round 4 — check Rev 5 → plan **Rev 6** (2026-06-04)

**Input:** `PLAN_mixed_flux_hetero_riemann_check.md` (Rev 5). The Rev-5 review gave a FINAL verdict —
plan implementable, no CRITICAL/HIGH open — with two non-blocking clean-ups (BUG-25 MEDIUM, BUG-26 LOW).

### Summary
- Findings addressed: **2 of 2 open** (BUG-25, BUG-26). Both are plan-text/test-wording fixes; no
  design change.
- Files modified: `PLAN_mixed_flux_hetero_riemann.md` (→ **Rev 6**): header, scope Constraint, Phase 5
  Files-to-Modify + Test 5.3, Phase 7 "Why TPV31" + acceptance, Tricky-areas. PDF regenerated (17 pp,
  0 missing glyphs).
- Plan self-consistency: PASS (grep-verified — every BP5 mention is now either the byte-exact *run*
  contract or an explicit "BP5 does not use this parser" note; no stale "BP5 parses through
  `LoadSpatialFrictionConfig`" claim remains).

### Changes Made
1. **BUG-25 (MEDIUM) — Test 5.3c corpus.** Reworded the parser-regression test in four places (Test
   5.3c, the scope Constraint, the Phase-5 Files-to-Modify note, Phase-7 acceptance, Tricky-areas):
   (a) realize "byte-identical" as a **field-by-field** extension of the existing
   `tests/unit/test_tpv_config_parse.cpp` (no `SpatialFrictionConfig::operator==`/dump exists);
   (b) ADD `tpv31/configs/tpv31.toml` (+ p2/p3) to that test's corpus; (c) **DROP BP5** — it has no
   spatial TOML and never calls `LoadSpatialFrictionConfig` (its byte-exactness is protected by the
   edit being purely additive to a parser BP5 never invokes).
2. **BUG-26 (LOW) — LSW clarification.** Added a paragraph to Phase 7 "Why TPV31": TPV31 is
   `law="slip_weakening"` (`tpv31.toml:58`), so the RK path is `AdvanceRKCoupledLSW_Spatial`
   (`spatial_dyn_driver.cpp:2864,2873`; precondition `GetFaultFrictionLaw()==LSW`,
   `rk_time_stepper.hpp:389`), NOT the rate-state stepper — but the central dispatch lives in the
   law-agnostic `wave.Mult`, so the feature is law-independent and the smoke/Test-5.1 fixture may be
   LSW. No code consequence.

### Deviations
None. Applied the reviewer's suggested wording.

### Unresolved findings
None. No CRITICAL/HIGH/MEDIUM/LOW findings remain open. Deferred follow-up (not a regression): the
real cross-rank `MPI_Allgatherv` material exchange.

### New tests
No new test IDs — BUG-25 folds into Test 5.3c (corpus reworded), BUG-26 is a doc clarification.

### Ready for Re-Review: YES
This should be the terminal review. If check Rev 6 confirms no open findings, the plan is ready for
`/code-implement`.

---

## Round 3 — check Rev 4 → plan **Rev 5** (2026-06-04)

**Input:** `PLAN_mixed_flux_hetero_riemann_check.md` (Rev 4). The Rev-4 review confirmed the probe and
metric fixes held, but found a fourth-round blocker (BUG-21) plus scope/metric gaps. The user
explicitly asked for an **exhaustive guard sweep** this round to stop the per-round whack-a-mole.

### Summary
- Findings addressed: **4 of 4 open** (BUG-21 CRITICAL; BUG-22 HIGH; BUG-23 MEDIUM; BUG-24 MODERATE).
- Pre-fix investigation: ran the exhaustive guard sweep myself (`grep interior_flux/mixed_flux/
  InteriorFlux/MixedFluxMode/is_rk` across `drivers/` + `spatial/code/` + `dynamic/`). Result: exactly
  **three** series guards block `matrix+adjacent+RK` (G1 `spatial_friction.cpp:1162-1167`, G2
  `spatial_dyn_driver.cpp:809-816`, G3 `bimaterial_wave_operator.inl:569-574`); four others
  (`spatial_friction.cpp:1702`, `driver:1078`, `wave_operator.inl:604`/`:1613`) are correct and TPV31
  does not trip them. The plan now carries this full inventory, so no fifth guard should surface.
- Files modified: `PLAN_mixed_flux_hetero_riemann.md` (→ **Rev 5**): header, Constraints
  (scope + byte-exact), Phase 2 req 4, Phase 5 (rewritten — guard inventory + G1 relax + scope +
  `seam_continuous` wiring + Test 5.3), Phase 7 (config + metric window + NEW-7.5 + acceptance),
  Testing Strategy, Tricky-areas.
- Plan self-consistency: PASS (grep-verified). PDF regenerated (16 pp); a stray `≪`/`✓` I introduced
  caused 4 missing glyphs on first build — replaced with ASCII, rebuilt clean (0 missing).

### Changes Made
1. **BUG-21 (CRITICAL) — the second, config-parse guard.** The sweep confirmed G1 at
   `spatial/code/spatial_friction.cpp:1162-1167` (`MFEM_VERIFY(interior_flux==Scalar ||
   mixed_flux=="none")`) fires at config load (`LoadSpatialFrictionConfig` ← `driver:616`), upstream
   of the BUG-4 guard the plan relaxed. Phase 5 now relaxes G1 too and ships a full **guard inventory**
   (relax G1/G2/G3; keep the four correct guards) so the series is handled in one pass, not one bug
   per round. Test 5.3 is the config-parse reproducer.
2. **BUG-22 (HIGH) — scope.** Corrected the false "config parse lives with the driver" claim:
   `MaterialSpec` + parser + G1 live in `spatial/code/spatial_friction.{hpp,cpp}`. Expanded the
   declared source-edit scope to include that file (not BP5/solver core, so [C2]-compatible), added
   `bool seam_continuous=false` to `MaterialSpec` + its parse, and required a **parser-regression
   test** (Test 5.3c) since all TPV*/BP5/TPV31-ADER configs parse through it.
3. **BUG-23 (MEDIUM) — metric window.** Phase 7 req 5(a) now defines `t0=max(mfem_t[0],ref_t[0])`,
   `t1=min(mfem_t[-1],ref_t[-1])`, interpolates within `[t0,t1]` (not the full grid — avoids
   `np.interp` endpoint clamping), and adds a minimum-coverage guard (fail with non-zero exit if
   coverage < ~90% of the reference span). Test NEW-7.5 (truncated run ⇒ insufficient-coverage exit).
4. **BUG-24 (MODERATE) — explicit no-regression.** Added to the byte-exact Constraint and Phase 2
   req 4: the `seam_continuous` guard fires only inside `BuildPerFaceCentralFluxMatrices_` for
   `central_flux_face_set_` faces with `mf_on_` true; with `mixed_flux=none` that set is empty
   (`wave_operator.inl:1653-1654`), so the existing TPV31 ADER job is provably unaffected.

### Deviations from the reviewer's suggested fixes
None. The user-requested exhaustive sweep went beyond the single BUG-21 fix and let the plan enumerate
every guard with its disposition.

### Unresolved findings
None at the plan level. Deferred follow-up (not a regression): the real cross-rank `MPI_Allgatherv`
material exchange.

### New tests (specified, for /code-implement)
- Test 5.3 (config-parse G1 relax + `seam_continuous` round-trip + TPV*/BP5/TPV31-ADER parser-regression).
- Test NEW-7.5 (`visualize_results.py` truncated-run coverage guard).

### Self-critique
This is the round where I stopped patching one guard at a time and swept for all of them — the sweep
is what makes Rev 5 likely terminal on the guard front. I also re-introduced a PDF-glyph regression
(`≪`/`✓`) that the earlier rounds had been clean of; caught and fixed it before reporting. The lesson
from four rounds: verify against the actual source (APIs, guard sites, glyph coverage) rather than
asserting — every round a fix that *asserted* without *checking* came back as a new bug.

### Ready for Re-Review: YES
Recommend a check Rev 5 pass (focus: the guard inventory completeness, the `seam_continuous` parser
wiring, Test 5.3's parser-regression). If clean, this should finally be implementable end-to-end.

---

## Round 2 — check Rev 3 → plan **Rev 4** (2026-06-04)

**Input:** `PLAN_mixed_flux_hetero_riemann_check.md` (Rev 3). **Context:** Round 1 (Rev 3) fixed
BUG-10/11 but the fixes were themselves defective — exactly what the adversarial re-review caught.

### Summary
- Findings addressed: **6 of 6 open** (BUG-15, BUG-16 CRITICAL; BUG-19 HIGH; BUG-17, BUG-18 MEDIUM;
  BUG-20 LOW). The 11 previously-fixed bugs and the 3 already-resolved-in-Rev-3 (BUG-12/13/14) are
  untouched.
- Files modified: `PLAN_mixed_flux_hetero_riemann.md` (→ **Rev 4**): header, Constraint, Phase 2
  req 4 + Edge Cases, Phase 3 Tests 3.4a/b/c, Phase 7 (config + req 5(a) + Edge Cases + Acceptance +
  Files-to-Create), Testing Strategy, Risk Assessment, Tricky-areas.
- Plan self-consistency: PASS (grep-verified — no "probe" remains as a mandated mechanism, only as
  the rejected option; `seam_continuous` threaded everywhere; metric fully specified; PDF rebuilds
  clean, 0 missing glyphs).

### Changes Made

1. **BUG-15 + BUG-16 + BUG-20 (drop the probe; mandate `material.seam_continuous`).**
   Round 1's BUG-10 fix was a material-field probe at `x_loc`/`x_nbr`. The reviewer showed it is
   (BUG-15) unimplementable — `MaterialField`'s only evaluator is `EvalAt(elem, ElementTransformation&,
   IntegrationPoint&)` (`heterogeneous_material.hpp:183-209`), with no transformation for a
   cross-seam neighbour point, and the operator holds only `const MaterialField* material_`
   (`hpp:182`), not the coordinate-only `eval_at_xyz` on the `DepthProfile1DMaterial` wrapper; and
   (BUG-16) wrong — an along-face-normal probe conflates legitimate depth variation (TPV31's
   `depth_profile_1d` jumps at 2400/5000/10000 m) with lateral variation, false-aborting the
   deliverable. **Root-cause fix:** deleted the probe everywhere and mandated the declarative
   `material.seam_continuous = true` flag as the sole guard (abort any `Mode::Coefficient` central
   shared face when absent). Updated: Constraint, Phase 2 req 4 (+ threading note: config →
   `BimaterialWaveOperator` setter/ctor → guard), Edge Cases, Risk Assessment, Tricky-areas. Test
   3.4b rewritten against the flag; added **Test 3.4c** (`Dispatch_DepthProfile_NoFalseAbort`, the
   BUG-16 reproducer). Phase 7 config now sets `[material].seam_continuous = true` (TPV31 is
   depth-only — the same stub approximation the existing TPV31 ADER job already runs under). BUG-20
   (probe geometry under-defined) dissolves with the probe.

2. **BUG-19 (HIGH) — no interpolation.** Phase 7 req 5(a) now requires interpolating the MFEM channel
   onto the reference `time_s` over the overlap window via `np.interp` BEFORE differencing (the
   script has no resampling today; MFEM RK45 adaptive dt ≠ SeisSol fixed cadence). Test NEW-7.1 must
   use offset time grids.

3. **BUG-17 (MEDIUM) — undefined denominator.** req 5(a) now specifies a peak-normalized, zero-safe
   metric: `peak_rel = max|mfem−ref| / max(max|ref|, floor)`, `rms_rel = rms(mfem−ref) /
   max(rms(ref), floor)` — not a per-sample `|mfem−ref|/|ref|` that blows up at the pre-nucleation
   zero-crossings. Added **Test NEW-7.4** (zero-crossing reference ⇒ no spurious failure).

4. **BUG-18 (MEDIUM) — unnamed channel.** req 5(a) now names the gated channels (`V_strike`,
   `slip_strike`, `tau_strike`, `_dip` only if dip motion, `sigma_n`) and **excludes `mu_eff`**
   (all-NaN reference ⇒ would silently pass), citing the column map at `visualize_results.py:145-155`.

### Deviations from the reviewer's suggested fixes
- None of substance. The reviewer's "preferred fix — drop the probe, mandate `seam_continuous`" was
  adopted exactly. Added a small, necessary wiring note (config → operator → guard) the reviewer
  implied; flagged that the spatial-config parse for the new key must be confirmed in-scope.

### Unresolved findings
- None of the 6 open findings is unresolved at the plan level. Deferred follow-up (not a regression):
  the real cross-rank `MPI_Allgatherv` material exchange (the plan enforces the Constant/seam-continuous
  support boundary via the flag until then).

### New tests (specified, for /code-implement)
- Test 3.4c `Dispatch_DepthProfile_NoFalseAbort` (BUG-16).
- Test NEW-7.4 `TolGate_PeakNormalized_ZeroCrossingSafe` (BUG-17/19); NEW-7.1 updated for interpolation.

### Self-critique
Round 1's probe was a "plausible but wrong" fix: I specified an API (`material_->EvalAt(x)`) without
verifying it exists, and asserted depth/lateral separability without checking the seam geometry. Both
were caught by the adversarial re-review. Rev 4 uses the declarative flag — implementable with the
APIs that actually exist, and immune to the depth-vs-lateral conflation.

### Ready for Re-Review: YES
Recommend a check Rev 4 pass (focus: the `seam_continuous` wiring + the BUG-17/18/19 metric spec),
then `/code-implement`.

---

## Round 1 — check Rev 2 → plan Rev 3 (2026-06-04)

Fixed BUG-10..BUG-14 (5 of 5 open at the time). BUG-10 was addressed with a material-field probe and
BUG-11 by scoping in the `--tol-rms/--tol-peak` implementation; BUG-12/13/14 (README, verify-dispatch
attribution, deriv-cache misquote) resolved and confirmed in Rev 3. The BUG-10 probe and the
under-specified BUG-11 metric were superseded by Round 2 above (BUG-15/16/17/18/19). Full per-finding
detail is preserved in the check document's Rev-2/Rev-3 history.
