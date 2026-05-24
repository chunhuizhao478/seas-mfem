# Code Review: DEBUG_speckle_normal_velocity_jump_2026-05-23.md — analysis + next move to chase [[v_n]] leakage

**Date:** 2026-05-23
**Reviewer role:** adversarial audit of the analysis and its proposed next tests (not a code-change review).

## Review Scope
- **Document under review:** `DEBUG_speckle_normal_velocity_jump_2026-05-23.md` (the [[v_n]]-leakage root-cause claim + §6 next tests + §7 fix direction).
- **Code consulted to verify the claims:**
  - `dynamic/fault_face_flux.cpp` — `ComputeTrialTraction` (49-85), `BuildImposedState` (261-298), `WriteBackState` (300-313).
  - `dynamic/tpv205_substep_iterator.cpp` — `[SLIP]` decomposition trace (415-448; `sn_vjump`/`sn_sterm` at 427-431).
  - `dynamic/wave_operator.inl` — `EvaluateBulkAtFaultQPsCanonical` interior (1908-1965) and shared (2143-2240) loops.
  - `fault/fault_basis.hpp` — `ComputeOrientedFrame` (411-512), `ComputeQPBasis` / `ComputeQPBasisShared` / `AppendSharedFaces`.
  - `drivers/spatial_dyn_driver.cpp:814` — `pmesh.SetCurvature(cfg.mesh.order)` (mesh IS curved).
  - `tests/unit/test_fault_basis_qp_orthonormality.cpp` — existing passing orthonormality regression.
- **Domain context:** project `CLAUDE.md` (sign conventions), `seas/CLAUDE.md` (FaultBasis = `t1=dip, t2=strike`), memory `project_safs_orientation_blowup.md` (a *canonicalization/sign* frame issue — distinct from the *orthonormality/accuracy* issue this doc hypothesizes).

---

## Headline verdict

The document is **right that there is a bug** (a welded fault should not show a ~41% normal-to-tangential velocity ratio) and it has **two genuinely solid pieces** (the `sn_vjump`/`sn_sterm` decomposition; the cap-run falsification of "tension is the amplifier"). But it **mis-localizes the root cause to "curvilinear shared-face frame inaccuracy," and that hypothesis is refuted by the code from three independent directions.** Its "most direct, decisive" next test (test #1) is a **mathematical tautology** that cannot fail and will mislead the next debugging round — the same anti-pattern already flagged in `REVIEW_PLAN_speckle_slip_runaway` (the slip tripwire that could not fire).

The chase should be **re-pointed away from the frame** and toward the **raw bulk predictor traces on shared faces** (the per-substep CK/ghost path, R-1303/R-1601 family).

---

## Findings

### [R-001] CRITICAL [fault_face_flux.cpp + fault_basis.hpp] — §2.3 mechanism (`[[v_n]] = V·can_n`) is mathematically impossible; the friction slip cannot inject a normal-velocity jump

**Category:** BUG (analysis error — wrong physical mechanism)

**Description:**
§2.3 derives the leak as the friction slip-rate vector projected onto a non-orthonormal normal:
`[[v_n]] = V·can_n = V1(can_t1·can_n) + V2(can_t2·can_n)`, claiming this is non-zero because the frame is inaccurate. This is wrong **twice over**, and neither path can carry a non-zero normal jump:

1. **The imposed state welds the normal jump to exactly 0, in any frame.** `BuildImposedState` (fault_face_flux.cpp:278, 286) sets
   `Q_imp_minus[VX] = Q_minus[VX] − (1/Z_p⁻)(σ_n_corr − Q_minus[SXX])` and
   `Q_imp_plus[VX]  = Q_plus[VX]  + (1/Z_p⁺)(σ_n_corr − Q_plus[SXX])`.
   Since friction never modifies σ_n (`σ_n_corr = σ_n_trial`, confirmed WriteBackState:312 and the doc's own §2.1), substitute `σ_n_trial = η_p(Q⁻_VX − Q⁺_VX + Q⁺_SXX/Z_p⁺ + Q⁻_SXX/Z_p⁻)` with `η_p = 1/(1/Z_p⁺ + 1/Z_p⁻)` (harmonic mean, doc §2.1). Then
   `Q_imp_minus[VX] − Q_imp_plus[VX] ≡ 0` **algebraically, regardless of frame**. The friction solution sets the two tangential imposed velocities (`VY`, `VZ`); the normal channel is welded shut by construction.

2. **The frame is orthonormal by construction, so `V·can_n ≡ 0` anyway.** `ComputeOrientedFrame` (fault_basis.hpp:447-479) builds `strike = normalize(up × n)` and `dip = normalize(strike × n)`. Therefore `can_t2·can_n = (up×n)·n ≡ 0` and `can_t1·can_n = (strike×n)·n ≡ 0` to machine precision. The friction slip `V = V1·can_t1 + V2·can_t2` lies in the plane spanned by the two tangents, which is *by definition* ⊥ `can_n`. So `V·can_n = 0` exactly — there is no geometric route for tangential slip to acquire a normal component.

**What the trace actually measures:** `sn_vjump = η_p·(Q_tilde_minus[VX] − Q_tilde_plus[VX])` (tpv205_substep_iterator.cpp:427-428). `Q_tilde_±` are the **ADER predictor bulk traces** at the QP — the time-extrapolated *wavefield* velocities, **not** the friction slip vector. So `sn_vjump = 5.2 m/s` means the **bulk DG solution genuinely has a 5.2 m/s normal-velocity discontinuity in the predictor state at that QP**. That is a property of the solution/predictor, not a projection of the friction slip and not a frame defect.

**Trigger:** Any reader takes §2.3 at face value and spends the next cycle "fixing frame orthonormality" or the friction slip embedding.

**Actual behavior (of the doc):** Attributes a bulk-predictor trace discontinuity to a friction-slip→normal projection.

**Expected behavior:** State the mechanism as: *the predictor bulk-trace normal velocities `v_n⁻`, `v_n⁺` disagree by 5.2 m/s at this shared QP; since the frame is provably orthonormal and the imposed state welds the normal jump, the disagreement lives in the raw bulk traces (predictor/ghost path), not in the frame projection.*

**Suggested fix (to the analysis + a guard in code):**
```diff
- [[v_n]] = V·can_n = V1(can_t1·can_n) + V2(can_t2·can_n)   (frame leak)
+ sn_vjump = η_p·(v_n⁻ − v_n⁺) where v_n± are the PREDICTOR BULK traces
+ rotated by an orthonormal can_n; V·can_n ≡ 0 (orthonormal) and the
+ imposed [[v_n]] ≡ 0 (BuildImposedState welds it). The 5.2 m/s lives in
+ the raw bulk predictor traces, not the frame.
```
Optionally add a cheap debug assertion proving the welding identity so the claim is self-documenting:
```cpp
// In BuildImposedState, under SEAS_DIAG_FAULT_FLUX, assert the welded normal jump:
const real_t vn_imp_jump = Q_imp_minus[VX] - Q_imp_plus[VX];
MFEM_ASSERT(std::abs(vn_imp_jump) < 1e-9 * (std::abs(Q_imp_minus[VX]) + 1.0),
            "Imposed normal-velocity jump must be ~0 by construction");
```

**Test case:**
```cpp
// test_R001_imposed_normal_jump_is_zero: for random Q_plus/Q_minus and Zp,
// after Evaluate(), assert Q_imp_minus[VX] - Q_imp_plus[VX] == 0 (to 1e-12 rel),
// confirming friction slip cannot inject a normal-velocity jump.
TEST(FaultFaceFlux, R001_imposed_normal_jump_welded) {
    // build DOFData with arbitrary Zp_plus != Zp_minus, eta_p = harmonic mean
    // random Q_plus, Q_minus; run BuildImposedState via Evaluate
    EXPECT_LT(std::abs(Q_imp_minus[VX] - Q_imp_plus[VX]),
              1e-12 * (1.0 + std::abs(Q_imp_minus[VX])));
}
```

---

### [R-002] CRITICAL [DEBUG doc §6 test #1] — the "decisive" `can_n` orthogonality test is a tautology and cannot detect the hypothesized error

**Category:** BUG (proposed test is incapable of falsifying its own hypothesis)

**Description:**
Test #1 proposes logging `can_n·can_t1`, `can_n·can_t2` ("should be 0"), `|can_n|`, and `can_n` vs the geometric `CalcOrtho` normal — calling it the "most direct" test. Every one of those four quantities is fixed to its ideal value **by construction**:
- `can_n·can_t1 = n·(strike×n) ≡ 0` (scalar triple product, repeated vector) — fault_basis.hpp:461-477.
- `can_n·can_t2 = n·(up×n) ≡ 0` — fault_basis.hpp:449-457.
- `|can_n| = 1` — normalized at fault_basis.hpp:433.
- `can_n` **is** the `CalcOrtho` normal (fault_basis.hpp:96/134/166 → `ComputeOrientedFrame` → `normal = n_raw/|n_raw|`), up only to the rank-independent `sign_flipped` canonicalization.

There is **already a passing regression test** proving this: `tests/unit/test_fault_basis_qp_orthonormality.cpp:69-70,201-212` asserts `max(|n·t1|,|n·t2|,|t1·t2|) ≤ 16 ULP` at every QP. So test #1 will report ~16 ULP and "pass," and the doc will then incorrectly conclude "frame is fine ⇒ it must be a physical opening (cause 3)," which is also wrong (see R-004). A 41% leak implies `sin(θ)=0.41 ⇒ θ≈24°`; a 24° non-orthonormality is **6×10¹⁵ ULP** off — the existing test would have screamed long ago.

**Trigger:** Running test #1 on the seed QP.

**Actual behavior:** Returns ~16 ULP (tautological pass) for *any* mesh, buggy or not.

**Expected behavior:** A test that measures the **raw bulk-trace** disagreement, which is frame-independent (see Recommended Next Move).

**Suggested fix (replace test #1):**
```diff
- 1. can_n orthogonality at qp=480: log can_n·can_t1, can_n·can_t2 (should be 0)...
+ 1. RAW bulk-trace normal-velocity jump (frame-independent): at the seed QP log
+    the GLOBAL-frame bulk velocity vectors v_self[3], v_nbr[3] (before rotation)
+    and the geometric normal n_geom[3]. Compute jn = (v_self - v_nbr)·n_geom.
+    If |jn| ~= 5.2 m/s the opening is real in the bulk DG solution; the frame
+    (orthonormal by test_fault_basis_qp_orthonormality) is NOT the cause.
```

**Test case:**
```cpp
// test_R002_orthonormality_is_already_guaranteed: demonstrates test #1 is a
// tautology by sampling 1000 random normals through ComputeOrientedFrame and
// asserting |n·t1|,|n·t2| < 16 ULP -- i.e. test #1 cannot produce a nonzero.
TEST(FaultBasis, R002_test1_is_tautology) {
    for (int i = 0; i < 1000; ++i) {
        Vector n_raw = random_unit_perturbed();
        real_t n[3], t1[3], t2[3]; bool sf; real_t nl;
        FaultBasis::ComputeOrientedFrame(n_raw, 3, ref, up, n, t1, t2, sf, nl);
        EXPECT_LT(std::abs(dot3(n,t1)), 16*DBL_EPSILON);
        EXPECT_LT(std::abs(dot3(n,t2)), 16*DBL_EPSILON);
    }
}
```

---

### [R-003] MODERATE [DEBUG doc §2.3 cause (2), §7] — "± sides use different normals (ghost frame ≠ local frame)" is refuted by the rotation code

**Category:** DEVIATION (claim contradicted by source)

**Description:**
The doc lists as cause (2) that the + and − sides use **inconsistent normals** on the shared face, and §7's fix calls for making the frame "consistent across the ± (local/ghost) representations." But `EvaluateBulkAtFaultQPsCanonical` rotates **both** sides with a **single** `Tinv_can` per QP:
- Interior loop: `Q_self_can`, `Q_nbr_can` both via the same `Tinv_can` (wave_operator.inl:1939-1952).
- Shared loop: same — one `Tinv_can` for `Q_self` and `Q_nbr` (wave_operator.inl:2207-2222).
- The shared path deliberately uses `qpd.sign_flipped`, which the comment (wave_operator.inl:2190-2195) states is **rank-independent** (set at FaultBasis ctor from face geometry) precisely so both ranks agree across the seam.

So within a rank there is **no** ±-side normal mismatch, and across ranks the frame is canonicalized to agree. Any residual cross-rank difference is the ~1e-14 geometry round-off (Phase-0, job 7747304), **not** a 41% leak. Cause (2) as written cannot produce the symptom.

**Trigger:** §7 fix work that tries to "reconcile ± normals."

**Suggested fix:**
```diff
- 2. the + and − sides using different normals (ghost frame ≠ local frame)...
+ 2. [REFUTED by code] EvaluateBulkAtFaultQPsCanonical rotates BOTH sides with a
+    single Tinv_can (wave_operator.inl:2207-2222) and shared uses the rank-
+    independent sign_flipped; no ±-mismatch beyond ~1e-14. Drop this cause.
```

**Test case:**
```cpp
// test_R003_single_frame_rotates_both_sides: assert that for a shared fault QP,
// the rotation matrix applied to Q_self and Q_nbr is bit-identical (same can_n).
// (Instrument EvaluateBulkAtFaultQPsCanonical to expose Tinv_can per side and
//  compare; expect exact equality.)
```

---

### [R-004] MODERATE [DEBUG doc §2.4, TL;DR, §7] — "rate-proportional ⇒ frame, not physical opening" is a non-sequitur; the root-cause is presumed before the evidence

**Category:** ASSUMPTION (unsound inference closing the analysis prematurely)

**Description:**
§2.4 argues the leak must be a geometric frame error (not a dynamic opening) because it is "instantaneous & rate-proportional" (≈41% of V from the first sub-step). But a **DG bulk-trace consistency error** and an **under-resolved / predictor-induced opening at the rupture front** are *also* instantaneous and rate-proportional — anything driven by the local slip motion scales with V. Rate-proportionality therefore does **not** discriminate "frame" from "bulk-trace/predictor opening." Combined with R-001 (frame projection is identically zero) and R-002 (orthonormality is guaranteed), the only mechanism left standing is exactly the one §2.4 dismisses: a **real normal-velocity discontinuity in the bulk predictor traces**. Yet the TL;DR declares "root cause localized" to the frame and §7's primary fix is "frame accuracy" — presuming an outcome the doc's own standard ("must be measured, not assumed," §2.4) forbids.

**Trigger:** Committing to a frame fix on the strength of §2.4.

**Suggested fix:**
```diff
- TL;DR ... is driven by a normal-velocity jump (opening) on the curvilinear
- SHARED fault [frame inaccuracy].
+ TL;DR ... is driven by a normal-velocity jump in the PREDICTOR BULK TRACES on
+ shared faces. The frame is orthonormal+geometric by construction (R-001/R-002),
+ so the jump is a property of the bulk/predictor solution, locus UNCONFIRMED
+ pending the raw-trace measurement (Recommended Next Move).
```

**Test case:** N/A (documentation logic) — discriminator is the raw-trace diagnostic in R-002's fix and the Recommended Next Move.

---

### [R-005] MODERATE [DEBUG doc §6 test #3] — "zero the dip channel" is premised on the impossible §2.3 path and is confounded

**Category:** EDGE_CASE (test built on a refuted mechanism; result will be misread)

**Description:**
Test #3 (`SEAS_FORCE_V1_ZERO`) is justified by §2.3's `[[v_n]] = V1(can_t1·can_n) + ...` — i.e. it expects zeroing the dip slip to remove a `V1·can_t1·can_n` injection. That injection is identically zero (R-001). What zeroing V1 *actually* does is suppress dip-slip → reduce the radiated dip-tangential **bulk** velocity → indirectly reduce the bulk-trace jump *if* the predictor opening is fed by dip motion. So a positive result would be **misattributed** to "dip→normal frame leak" when it is really "less slip ⇒ less radiated bulk velocity." The test still has *some* diagnostic value (which tangential direction drives the bulk opening) but its stated rationale is wrong and the result is confounded by the slip-reduction side effect.

**Suggested fix:**
```diff
- 3. Zero the dip channel: if [[v_n]] collapses, the leak is dip-slip → normal...
+ 3. Zero the dip channel: diagnostic for WHICH tangential bulk motion feeds the
+    predictor opening. Caveat: zeroing V1 also reduces total slip/radiation, so a
+    collapse is confounded (less radiated bulk velocity, not a frame projection).
+    Prefer the raw-trace decomposition (test #1 replacement) which is unconfounded.
```

**Test case:** N/A (experiment-design caveat).

---

### [R-006] MODERATE [DEBUG doc §5, "R-007d resolution refuted"] — D_c variation does not refute mesh-h under-resolution

**Category:** ASSUMPTION (over-claimed refutation that risks closing a live lever)

**Description:**
§5 claims "cohesive-zone under-resolution refuted" because D_c=2 vs D_c=8 gave a bit-identical seed. But `D_c` is the slip-weakening **distance**, not the mesh size `h`. Increasing `d_c` *enlarges* `L_nuc` (and the cohesive zone width) at **fixed h**, so the quoted "L_nuc/h ≈ 12→48" reflects a larger *physics* length, not a finer mesh. Worse, the doc itself says the seed fires at `δ≈0.01 m ≪ d_c`, i.e. **before weakening engages** — so `d_c` *cannot* affect the seed by construction, making the "bit-identical seed" trivially expected and **not** evidence about spatial resolution. Mesh-`h` was never varied. The SeisSol/TPV14 comparison's conclusion (process-zone / `L_nuc/h` is the upstream lever; `REVIEW_seissol_drdg3d_tension_comparison`) is **not** refuted by this experiment.

**Suggested fix:**
```diff
- R-007d cohesive-zone under-resolution. D_c 2→8 (L_nuc/h ≈12→48) ... Finer mesh
- is not the lever.
+ R-007d: D_c (slip-weakening DISTANCE) 2→8 left the seed bit-identical -- but the
+ seed fires pre-weakening (δ≈0.01 m ≪ d_c) so d_c cannot affect it by
+ construction; and d_c ≠ mesh h. This refutes "d_c sensitivity of the seed," NOT
+ mesh-h under-resolution, which was not varied. R-007 (refine h) remains open.
```

**Test case:** N/A (refutation scope) — the real test is an `h`-refinement run (halve element size at fixed `d_c`).

---

### [R-007] MODERATE [DEBUG doc §1, §2.2, TL;DR "shared-only"] — the "shared-only" claim is confounded by the `V_abs > threshold` trace gate (selection bias)

**Category:** ASSUMPTION (observation cannot support the conclusion drawn)

**Description:**
The `[SLIP]` trace fires only when `s.V_abs > slip_v_thr` (tpv205_substep_iterator.cpp:415; threshold 10 m/s). So the trace only ever shows QPs that have **already** crossed into runaway — and those happen to be shared. This establishes "the runaway QPs are shared," **not** "interior QPs have no `[[v_n]]` jump." Interior QPs with a moderate (sub-threshold) bulk-trace jump are invisible to the gated trace. The conclusion "frame defect specific to shared faces" (which drives the whole §7 fix) rests on this selection-biased sample. (It may well be shared-dominated — the 1.76% shared-DOF step-forcing warning is a plausible *separate* reason — but it is not established by this trace.)

**Suggested fix:** Add an **ungated** measurement: at a fixed early time, histogram `|[[v_n]]|/V` (raw-trace, R-002 definition) over **all** interior vs all shared fault QPs. Only then is "shared-only" supportable.
```diff
+ Caveat: the [SLIP] trace gates on V_abs>10 m/s, so it samples only already-
+ runaway QPs. "shared-only" needs an ungated all-QP |[[v_n]]|/V histogram
+ (interior vs shared) at fixed early t before it can be asserted.
```

**Test case:** N/A (run-time diagnostic; add an env-gated all-QP histogram pass).

---

### [R-008] LOW [POSSIBLE] [DEBUG doc §4] — iterator-vs-output decoupling under-explains why the macro solve stays bounded under a single shared frame

**Category:** ASSUMPTION (incomplete causal story; points at the real locus)

**Description:**
§4 attributes the bounded `fault.vtkhdf` σ_n to the reconcile/macro-solve **overwriting the written field**. But the written field is separate from the σ_n the macro flux *computes*. If the cause were the frame (or any QP-static geometry), the macro-dt shared solve (`ComputeADERSharedFaceFluxRHS`) uses the **same** `FaultBasis` frame at the **same** QP, so it should compute the **same** tensile σ_n — yet §4 says its `sigma_n_trial` "stays small/bounded." The only way both are true is that the difference is **not** QP-static: the iterator collapses because it runs the **per-sub-step predictor** `Q_tilde` (CK recursion, element-local; ghost coverage via R-1303/R-1601), while the macro solve uses the base/exchanged `Q`. That squarely localizes the bulk-trace opening to the **per-sub-step predictor/ghost path on shared faces** — i.e. the actual mechanism, and *not* the frame. The doc gestures at this in §3 (loop gain) but §4's "overwrite" framing obscures it.

**Suggested fix:** Reframe §4 around the predictor: *the iterator's per-sub-step predictor `Q_tilde` on the two sides of the shared face diverges in `v_n`; the macro solve (base Q) does not — hence the iterator collapses while the output (macro) stays bounded. This implicates R-1303/R-1601 ghost coverage of the sub-step predictor, not the frame.*

**Test case:** see Recommended Next Move item (B) — compare iterator `Q_tilde_±` vs macro base `Q_±` at the same shared QP and time.

---

## Recommended Next Move (replaces the frame-chase)

Because R-001/R-002/R-003 remove the frame as a possible cause, re-point the diagnostics at the **raw bulk predictor traces**:

**(A) Frame-independent raw-trace decomposition at the seed shared QP.**
Log, in the **global** frame (before any rotation), the predictor bulk velocity vectors `v_self[3]`, `v_nbr[3]` and the geometric normal `n_geom[3]`; report `jn = (v_self − v_nbr)·n_geom` and the two tangential jumps. This is rotation-invariant, so it isolates whether the 5.2 m/s opening is real in the DG solution (it will be) and how it splits across the global axes — with **zero** dependence on the (provably orthonormal) frame.

**(B) Iterator-predictor vs macro-base Q comparison at the same QP/time (the decisive read).**
For the seed shared QP, print `Q_tilde_minus[VX]`, `Q_tilde_plus[VX]` (iterator, per-sub-step predictor) **and** the macro-solve `Q_minus[VX]`, `Q_plus[VX]` (`ComputeADERSharedFaceFluxRHS`, base/exchanged Q) at the same sub-step. If the iterator's jump is large while the macro's is small (as §4 implies), the opening is generated by the **per-sub-step predictor / ghost coverage** (R-1303/R-1601), which is the genuine shared-only mechanism — not the frame.

**(C) Ungated shared-vs-interior `|[[v_n]]|/V` histogram (settles R-007).**
At a fixed early time, over all fault QPs, to test "shared-only" without the `V_abs>10` selection bias.

**(D) Mesh-h refinement (settles R-006), and `dt` tightening (test #5).**
The two remaining *cheap* discriminators the doc has not run: halve `h` at fixed `d_c` (does the bulk-trace opening shrink ⇒ DG under-resolution), and halve the macrostep (does it shrink ⇒ the 1.76% shared-face step-forcing is the per-sub-step driver). Both are consistent with a predictor/ghost-trace origin and neither touches the frame.

---

## What is solid in the document (do NOT "fix" these)

- **The `sn_vjump`/`sn_sterm` decomposition** (tpv205_substep_iterator.cpp:427-431) is correct and well-built; it cleanly shows the collapse is in the **velocity** term and the **bulk normal-stress term ≈ 0**, which correctly retires the earlier "corrupted ghost bulk STRESS" hypothesis. Keep it.
- **The cap-run falsification of "tension is the amplifier"** (§3, job 7747835): a genuinely well-designed experiment — the compressive `qp=481` running away immune to `SEAS_NOOPENING` is a clean refutation. The conclusion "the amplifier is the sign-indefinite `[[v_n]]↔traction` feedback, not free-slide-under-tension" is sound and important. Keep it; it correctly kills the no-opening fix direction.
- **The instinct that ~41% normal/tangential is unphysical for a welded fault** is correct — there IS a real bug in the bulk-trace normal velocity. The only error is localizing it to the frame.

---

## Summary
- **Critical issues:** 2 (R-001 impossible mechanism; R-002 tautological "decisive" test)
- **Moderate issues:** 5 (R-003 refuted ±-mismatch; R-004 presumed root cause; R-005 confounded test #3; R-006 over-claimed resolution refutation; R-007 selection-biased "shared-only")
- **Low issues:** 1 (R-008 incomplete §4 story — but it points at the true locus)
- **Plan compliance (analysis ↔ code):** PARTIAL — decomposition + cap experiment are faithful and informative; the frame root-cause and test #1 are not supported (and are contradicted) by the code.
- **Verdict:** **FAIL — must re-point before proceeding.** The frame hypothesis and test #1 cannot be acted on; chasing them will burn a cycle on a provably-correct frame. Adopt the raw-trace diagnostics (A)/(B) which target the per-sub-step predictor/ghost path — the only mechanism consistent with all the evidence (welded imposed state, orthonormal frame, single rotation, shared-dominated, rate-proportional).

## Unreviewed Areas
- `ComputeADERSharedFaceFluxRHS` macro-solve internals (the "bounded" path in §4) were not line-read here beyond confirming it shares the `FaultBasis` frame; item (B) above is the recommended way to verify the predictor-vs-base divergence directly.
- The `ComputeADERSubStepStates` CK recursion and exactly when ghost `Q_tilde` is (or is not) refreshed per sub-step (R-1303/R-1601 coverage) — the likely true locus — was not audited in this pass; it is the recommended next read after (A)/(B) confirm the bulk-trace origin.
