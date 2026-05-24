# Debug findings: SAFS slip "speckle" runaway → normal-velocity-jump (opening) on the curvilinear shared fault — 2026-05-23

## TL;DR (root cause localized)

> **ROOT CAUSE (2026-05-24, confirmed; see `REVIEW_speckle_seissol_drdg3d_rootcause_2026-05-24.md`,
> the R-008 UPDATE blocks below, and §8):**
> **(R-001)** on **shared (MPI-boundary) fault faces** SAFS violates the SeisSol/drdg3d
> single-slip-rate invariant — the iterator weakens μ off the **predictor** slip while the
> applied flux is a **separate** re-solve on `I/dt` (the iterator's `I_imp` is *discarded*,
> `wave_operator.inl:4856-4895`); interior faces are correct (consume `I_imp`, :3871-3885).
> **(R-004, CONFIRMED by `test_ghost_exchange_bynodes_vs_scalar`)** the seed of the predictor
> divergence is the R-1601 `vdim=NUM_STATE` byNODES `FaceNbrData()` unpack reading **scrambled
> neighbour components**. **The frame / tension / friction-law narrative in the rest of this
> §TL;DR is SUPERSEDED** (frame refuted; tension refuted by the cap run; free-slide matches the
> references). Fix order: R-004 then R-001 (see §8.5).

The post-cross-rank-reconcile slip runaway (`max_slip → 1–3×10⁶ m` while sampled
`V_max ≈ 4 m/s`) is driven by a **normal-velocity jump (opening) on the curvilinear
SHARED fault**, which the Godunov flux turns into a tensile trial normal traction,
which the LSW **free-slide-under-tension** path then runs away.

- The σ_n collapse is **~100% the normal-velocity-jump term** `η_p·(v_n⁻−v_n⁺)` and
  **~0% the bulk normal-stress term** (decomposition, job 7747812).
- The jump is `[[v_n]] ≈ 5.2 m/s` at the seed — **~41% of the tangential slip rate**,
  present **from the first sub-step**, **on shared faces only** (`is_shared=1`),
  **resolution-independent** (D_c=2 and D_c=8 give a bit-identical seed).
- Interpretation: **the friction-imposed tangential slip is leaking ~41% into the
  normal-velocity channel** because the canonical fault frame on the curvilinear shared
  face is geometrically inaccurate (`can_n` not ⊥ the slip plane, or the ± sides use
  inconsistent normals). This is the **Part-A family** (`ComputeOrientedFrame`/
  `FaultBasis` on the curved fault); Part A fixed which side is "+", but the frame's
  orthonormality/accuracy on shared+curved faces was not addressed.

**Refuted along the way:** cohesive-zone under-resolution (R-007d), the "our LSW
free-slide deviates from SeisSol" hypothesis (R-005d as a friction-law fix), and the
"corrupted ghost bulk stress" hypothesis. See §Refuted.

**UPDATE 2026-05-24 (run 7747889, `[FRAME]`+`[MACRO]` default-ON) — R-008 CONFIRMED: the
runaway is ~100% in the per-sub-step PREDICTOR; the time-integrated macro solve is bounded
and physical.** The instruments fired at the seed (rank 104, qp 480; 8705 `[SLIP]`, 11434
`[FRAME]` on rank 104, 5717 `[MACRO]` lines; 319413/68604/34302 total) and split the two
paths cleanly over the full run (t: 0→2 s, completed, no NaN; the tail `H5Fclose` errors are
benign checkpoint-close ref-count noise):

| quantity | path | end value (t≈2 s) | verdict |
|---|---|---|---|
| `V_max` (committed) | macro | ~5–7 m/s (final 6.97) | bounded ✓ |
| `[MACRO] sigma_n_corr` (written σ_n) | macro | +46 MPa (run range 46–52) | bounded ✓ |
| macro `dv_n` (`vn_minus−vn_plus`) | macro | 0.015 m/s | bounded ✓ |
| `[MACRO] slip_rate` | macro | ~2.3 m/s | bounded ✓ |
| `V_substep_max` | predictor | 1.10e7 m/s | **RUNAWAY** |
| `[SLIP] sn_vjump ≈ sigma_n_tot` | predictor | +3.07e13 Pa | **RUNAWAY** |
| `[FRAME] dv_n` (frame-indep.) | predictor | +3.84e6 m/s | **RUNAWAY** |
| `max_slip` = ∫(predictor V)·dt | predictor | 4.62e6 m | **RUNAWAY** |

- **Predictor vs integrated at the SAME QP:** predictor `dv_n`=3.8e6 m/s vs macro
  `dv_n`=0.015 m/s — a factor ~2.6e8. The time-integrated state is smooth; the per-sub-step
  predictor *traces* are discontinuous by millions of m/s. ⇒ **fix locus = the per-sub-step
  predictor / ghost path on shared faces (R-1303/R-1601)** — NOT the frame (refuted), NOT
  tension (refuted), NOT the integrated state.
- **Symptom mechanism nailed (§4 confirmed directly):** the driver `[DIAG]` reports BOTH
  `V_max` (committed, ~5 m/s) and `V_substep_max` (predictor, 1.1e7 m/s). `max_slip` is
  accumulated from the *predictor* V (`tpv205_substep_iterator.cpp:122`), so it integrates the
  runaway → 4.6e6 m, while `V_max`/σ_n (from the macro solve) stay physical. That is exactly
  why the "slip runaway" is visible in the slip field but invisible in V/σ_n. It is
  **spreading**: `n_rupturing(V>0.5)` climbs 1037→1882 and `max_slip` is monotone (more QPs
  join the predictor runaway each step = the speckle growing).
- **Onset → blowup:** predictor `dv_n` = −5.2 m/s at t=0.477 (matches the prior `sn_vjump/η_p`
  read) → +4.4e11 Pa by t=0.547 (one downsample window, ~70 ms; sign flips compressive) →
  monotone to +3.07e13 Pa. `mu_eff` weakens 0.847→0.300 (dynamic) by t=0.547 and stays.
- **Diagnostic validated:** `sn_vjump/dv_n` = 3.0725e13 / 3.838e6 = **8.0e6 = 0.5·ρ·c_p**, so
  the `[SLIP]`/`[FRAME]` identity `sn_vjump = η_p·dv_n` holds live; `sign_flipped=0` at the
  seed throughout (no canonicalization pathology — the frame refutation stands).

**Open question driving the root-cause `/code-debug`:** why does the shared-face predictor
`Q̃⁺` (neighbor/ghost) diverge from `Q̃⁻` by millions of m/s in the normal-velocity trace when
their time-integral `I/dt` agrees to 0.015 m/s? Locus = the CK extrapolation + ghost coverage
on the shared seam (R-1303/R-1601).

### Root-cause investigation (2026-05-24, `/code-debug`) — secular feedback, NOT CK overshoot

Static code map + log analysis of run 7747889:

- **Call flow (verified):** `spatial_dyn_driver.cpp:416` `ComputeADERSubStepStates(Q,…)`
  builds the element-local CK predictor as a **point value** at each sub-step midpoint τ_o
  → `:426` `EvaluateBulkAtFaultQPsCanonical(Q_per_node[o], …)` rotates to canonical and splits
  ±, doing the R-1601 per-sub-step ghost exchange for the `+` side → `:439`
  `iterator.AdvanceWithSubStepStates(…)` friction-solves on that predictor trial traction and
  accumulates `slip += V·dt` (`tpv205_substep_iterator.cpp:125-126`). The macro corrector
  (`ComputeADERSharedFaceFluxRHS`) instead consumes the **time-integral** `I/dt`.
- **R-1601 ghost-exchange code reads correct** (`wave_operator.inl:2030-2174`): per-sub-step
  `ExchangeFaceNbrData` of the predictor, `Q_self` from local DOFs, `Q_nbr` from the face-nbr
  ghost layer with matching `nbr_idx*ndof_per_el_` indexing; frame via `sign_flipped`. No
  staleness/indexing bug found by inspection.
- **The growth is SECULAR (across macro steps), NOT a within-step CK overshoot.** At the seed,
  within one macro step the predictor `sn_vjump` is flat (o=0→o=1: −4.170e7→−4.219e7, +1%);
  **across** consecutive steps it grows geometrically −4.17e7→−5.50e7→−7.09e7→−8.96e7→−1.11e8
  →−1.36e8 = **gain ≈1.32, 1.29, 1.26, 1.24, 1.23 per macro step** (dt≈3.5e-4 s). ⇒ a
  **closed-loop instability, per-step gain ≈1.25–1.32**, not a CK large-τ artifact (that would
  give o=1≫o=0). Confirms the §3 "loop gain >1/step" hypothesis with a number.
- **Macro side individually bounded all run:** `vn_plus∈[−0.016,0.495]`, `vn_minus∈[−0.123,
  0.231]` m/s. So the predictor jump (3.8e6 m/s) and the macro trial jump (0.015 m/s) diverge
  despite both nominally deriving from the same base `Q` — they cannot be reading the same
  effective state. Two surviving sub-hypotheses, **not yet discriminable from this log**:
  - **H-A (bulk pumping):** the iterator's runaway imposed flux pumps the *bulk* velocity
    field near the fault each step (gain ~1.3); the predictor reads the growing bulk state.
    The committed friction `V_max` (≈5 m/s) would not show it — it's the bulk DG velocity at
    the fault QP, not the friction output.
  - **H-B (predictor-only feedback / exchange divergence):** the predictor/iterator path sees
    a self-or-neighbor trace that the macro `I` exchange does not.
- **Decisive next measurement (diagnostic only, no functional change):** in the `[FRAME]`
  block (`wave_operator.inl` ~2196, `Q_self`/`Q_nbr` in scope at 2158-2173) additionally print
  the predictor `Q_self[VX]` and `Q_nbr[VX]` **separately** (global frame, pre-rotation), the
  base bulk `Q_data` normal velocity at that DOF, AND the predictor + macro **tangential**
  (`VY`/`VZ`) so the tangential-pumping branch is observable. Reads: self alone grows ⇒ local
  element / bulk pumping (H-A); nbr alone grows ⇒ ghost-exchange path; base `Q_data` grows ⇒
  bulk pumping confirmed (H-A); macro tangential also runs away ⇒ the bulk IS pumped (the
  `[MACRO]` normal-only trace hid it).

### Static-read continuation + SeisSol comparison (2026-05-24, `/code-debug "also check against SeisSol"`)

CK routines read in full (`wave_operator.inl`): `ComputeADERTimeIntegrated` (1205-1284, the
integral `I`) and `ComputeADERSubStepStates` (1302-1405, the point predictor `Q̃(τ_o)`) use the
**identical** element-local recursion `D(k+1)=-Σ_d A_d ∂_xd D(k)` (same `ApplySpatialDerivative`
+ `ApplyJacobianPerDOF`, same `A_d`, separate scratch buffers, correctly phased Taylor weights —
R-1503). **Both are strictly element-local; neither reads neighbor data** (comment :1370). For
the run's order 2 both reduce to `D(0)+c·D(1)`, so point-value and `I/dt` must agree per element
to the small `D(1)` term — consistent with the measured within-step flatness. ⇒ **the bug is NOT
in the CK recursion**, and since predictor-self ≈ macro-self ≈ base-Q (same polynomial), the 10⁸
predictor-vs-macro jump divergence is in the **neighbor ghost** the predictor reads (R-1601) vs
the macro's `I`, **or** a secular **tangential** bulk pumping invisible to the normal-only
`[MACRO]` trace (predictor `dv_t1`=6.1e6 > `dv_n`=3.8e6).

SeisSol (`/Users/chunhuizhao/projects/SeisSol`) cross-check:
- **Trial-traction formula is term-for-term identical.** `FrictionSolverCommon.h:182-184`:
  `normalStress[o] = etaP·(qIMinus[U] − qIPlus[U] + qIPlus[N]·invZp + qIMinus[N]·invZpNeig)` =
  our Eq. 7a `η_p·(Δv_n + σ_n⁺/Zp⁺ + σ_n⁻/Zp⁻)`. Our physics is right; the difference is
  upstream (what feeds `Δv_n`).
- **`etaDamp` knob confirms this term is a KNOWN DR-DG instability.** SeisSol multiplies `etaP`
  by `etaPDamp` (`FrictionSolverCommon.h:150,157`), sourced from the user param `etaDamp`
  (`DRParameters.h:88`, **default 1.0 = OFF**; applied only for `t < etaDampEnd`,
  `BaseFrictionLaw.h:70-71`). It exists specifically to damp the normal radiation term `η_p·Δv_n`
  — exactly our runaway term — but defaults OFF and SeisSol is stable without it. ⇒ a candidate
  *mitigation* (an `etaDamp` analog), NOT the reason SeisSol stays bounded.
- **One-state vs our two-state — the leading STRUCTURAL difference.** SeisSol evaluates ONE
  `qInterpolated[o]` at time-quadrature intervals (`spaceTimeInterpolation`, `DynamicRupture.cpp`)
  for BOTH sides (neighbor via communicated time-derivatives, evaluated locally at the SAME
  quadrature), and uses it for the friction solve, the slip integral
  (`accumulatedSlip += timeWeight·slipRate`, `FrictionSolverCommon.h:595/603`), AND the
  time-integrated `imposedState` flux (`postcomputeImposedStateFromNewStress`, same
  `timeWeights`) — **one consistent computation**. Our path drives the friction solve + the
  `slip += V·dt_sub` accumulator (`tpv205_substep_iterator.cpp:125`) from the per-sub-step
  **point** predictor `Q_pointwise` (`spatial_dyn_driver.cpp:426,439`), while the written
  output / trial decomposition comes from the **time-integrated `I`** path
  (`ComputeADERSharedFaceFluxRHS` consumes the installed `I_imp` for the flux at
  `wave_operator.inl:3871-3885`, else re-solves `EvaluateADER_LSW` on `I` at 4886-4892). The
  slip-driving point-predictor is the state that runs away; SeisSol has no separate point-
  predictor slip path.

**Status (UPDATED 2026-05-24 — RESOLVED):** the "not yet discriminable / needs a Frontera split
diagnostic" caveat above is **superseded**. The adversarial review
`REVIEW_speckle_seissol_drdg3d_rootcause_2026-05-24.md` localized it to **R-001** (shared-face
dual-solve / single-slip-rate violation — verified in code: `wave_operator.inl:4856-4895` discards
`I_imp`, interior :3871-3885 consumes it), which also **excludes H-A** (the iterator's runaway flux
is *discarded* on shared faces, so it cannot pump the bulk). The neighbour-side seed (**R-004**) is
now **CONFIRMED by a local np=2 test** (`test_ghost_exchange_bynodes_vs_scalar`, no Frontera run
needed): the R-1601 byNODES `FaceNbrData()` unpack reads scrambled neighbour components. Fix order
**R-004 → R-001** (§8.5). `etaDamp`-analog is a known mitigation, not the fix.

---

**UPDATE 2026-05-23 (code audit `REVIEW_DEBUG_speckle_normal_velocity_jump`) — the FRAME
root cause (§2.3, §2.4, §7, test #1) is REFUTED in code; corrections below.**
Verified against source: (1) `BuildImposedState` (`fault_face_flux.cpp:278,286`) welds the
imposed normal-velocity jump to **exactly 0** — with `sigma_n_corr=sigma_n_trial`
(`:219,:813,:928`), `sigma_n_trial` form (`:62-64`), and matched SAFS material
(`eta_p=Zp/2=harmonic mean`), `Q_imp⁻[VX]−Q_imp⁺[VX] ≡ 0` algebraically, any frame.
(2) The frame is orthonormal by construction (`fault_basis.hpp:447-479`) ⇒ `V·can_n ≡ 0`;
the friction slip CANNOT inject a normal jump. (3) Both sides rotate through a single
`Tinv_can` (`wave_operator.inl:2207-2222`) ⇒ no ±-normal mismatch. **Therefore `sn_vjump`
measures the ADER PREDICTOR bulk traces `Q̃±[VX]`, not `V·can_n`; the 5.2 m/s is a genuine
normal-velocity discontinuity in the bulk PREDICTOR, locus = the per-sub-step predictor /
ghost path on shared faces (R-1303/R-1601), NOT the frame.** Also: §5 "resolution refuted"
is overclaimed — `d_c≠h` and the seed fires pre-weakening (`δ≈0.01≪d_c`) so it is
`d_c`-independent by construction; mesh-`h` was never varied (h-refine remains open).
"shared-only" (§1/§2.2) is selection-biased by the `V_abs>10` `[SLIP]` gate. The shipped
`[FRAME]` trace's `dv_n` IS the correct raw-trace probe (A) — reinterpret it, do not discard.
Decisive next read = **(B)** predictor `Q̃±` vs macro-base `Q±` at the same shared QP —
**INSTRUMENT SHIPPED**: the `[MACRO]` trace (`ComputeADERSharedFaceFluxRHS`, right after the
`EvaluateADER_LSW` dispatch — relocated there per REVIEW.md R-001, NOT the interior
`ComputeADERFaceFluxRHS`; plan `PLAN_predictor_vs_macro_diag_2026-05-23.md`; env
`SEAS_DIAG_MACRO=1`, default-ON in the Dc2 sbatch, shares the `[FRAME]` seed locator) prints
the macro solve's `sn_vjump`/`sn_sterm` on the time-integrated canonical `I_{±}_local/dt` plus
the **written** `fdata.sigma_n_corr`/`slip_rate`. Compare `[MACRO] sigma_n_corr` vs `[SLIP]
sigma_n_tot` at the same shared QP: bounded MACRO + collapsing SLIP ⇒ opening is per-sub-step
predictor/ghost (R-1303/R-1601). Byte-exact when unset (ADER interior-vs-shared rel=0);
decomposition unit-tested (`test_macro_diag_decomposition_identity`, rel≤5e-16; in `make test`
+ `test-v92-regression-gates`). Awaiting Frontera run.
What stands: the `sn_vjump`/`sn_sterm` decomposition and the cap-run amplifier refutation.

**UPDATE 2026-05-23 (cap run 7747835) — amplifier hypothesis REFUTED.** The
`SEAS_NOOPENING=1` cap does **NOT** stop the runaway. It arrests the *tensile* QP
(qp=480: max_slip 0.013 m vs 3.11 m at t=0.49, a 250× cut) but a **COMPRESSIVE**
normal-velocity-jump runaway at the neighbour qp=481 (σ_n → +TPa, `sn_vjump` **positive**)
ignites independently ~13 ms later, is immune to the cap (`sigma_n_total` never ≤ 0), and
carries the **same** runaway: max_slip 441.5 m (cap) vs 439.9 m (no-cap) at t=0.525; 5335 m
at t=0.560. ⇒ The amplifier is the **sign-indefinite `[[v_n]] ↔ traction` feedback**, not
free-slide-under-tension. Tension was one QP's symptom; no-opening is necessary-but-not-
sufficient and **cannot fix this**. See §3 (updated) and §6 test #4.

**Pending confirmations:** a `can_n` orthogonality check at the seed QP (trigger test) and
the feedback-gain / `[[v_n]]`-coupling fix. See §Decisive next tests.

Continues: `PLAN_speckle_slip_runaway_2026-05-23.md`,
`REVIEW_speckle_tension_analysis_2026-05-23.md`,
`REVIEW_seissol_drdg3d_tension_comparison_2026-05-23.md`,
debug_document/`spatial_dynamic_rupture_speckle_blowup_2026-05-22.md`.

---

## 1. Evidence chain (Frontera runs, branch `safs`)

| Job | Config | What it showed |
|---|---|---|
| 7747525 | Dc2 XRANK-diag (pre-trace) | The symptom: `max_slip → 2.7e6 m`, `V_max ~4 m/s`, speckle ≤7.7 km of the nucleus; `fault.vtkhdf` `normal_stress` **never tensile** (min 30.4 MPa). |
| 7747766 / 7747767 | Dc2 / Dc8 resolution pair (`[SLIP]` trace) | **Seed bit-identical** between D_c=2 and D_c=8 (`V_abs` 12.57 vs 12.54, σ_n +5.997 vs +6.070 MPa, δ 0.0104 both) ⇒ **resolution refuted**. Both run away (Dc8 max_slip 1.19e6 m, only ~2× smaller). |
| 7747812 | Dc2 baseline, **decomposition trace** (commit `97dd189`) | `sigma_n_trial` split into `sn_vjump` vs `sn_sterm`: the collapse is **100% the velocity-jump term**. |
| 7747835 | Dc2 **`SEAS_NOOPENING=1`** cap run | Cap arrests the tensile QP (qp=480, 250× less slip at t=0.49) but a **compressive** runaway at qp=481 (σ_n→+TPa, `sn_vjump`>0) is immune and carries max_slip (441.5 m vs 439.9 m no-cap at t=0.525). **Amplifier hypothesis refuted** — runaway is sign-indefinite. |

Instrumentation: commits `22b77c3` (Phase-1 `[SLIP]` trace + honest sub-step V_max) and
`97dd189` (decomposition `sn_vjump`/`sn_sterm`, `rank`, `is_shared`, coords; `SEAS_NOOPENING`
cap; `SEAS_NOOPENING` env on the Dc2/Dc8 sbatches). Byte-exact when `SEAS_*` unset
(TPV205 parity 14/14, friction 75/75, reconcile 4/4).

---

## 2. The mechanism, with equations

### 2.1 The trial normal traction (Godunov/Riemann), `fault_face_flux.cpp:62-66`
```
sigma_n_trial = η_p · ( v_n⁻ − v_n⁺  +  σ_n⁺/Z_p⁺ + σ_n⁻/Z_p⁻ )
              = η_p·(v_n⁻ − v_n⁺)          +  (Z_p⁻σ_n⁺ + Z_p⁺σ_n⁻)/(Z_p⁺+Z_p⁻)
              = [sn_vjump]                 +  [sn_sterm]
```
- `Z_p = ρ·c_p` (P-impedance), `c_p = sqrt((λ+2μ)/ρ)`,
  `η_p = Z_p⁺Z_p⁻/(Z_p⁺+Z_p⁻)` (harmonic mean; `= Z_p/2` for matched sides).
- `VX`(=6)/`SXX`(=0) are the **normal** velocity/stress in the canonical fault frame
  (`wave_state.hpp:27-37`).
- Physical reading: `σ_n* = (impedance-weighted avg normal stress) + (radiation damping)·(opening rate)`.
- SAFS: `ρ=2670, c_p=5996 ⇒ Z_p=1.60e7, η_p=8.0e6 Pa·s/m`.

Assembly (`tpv205_substep_iterator.cpp` StepOneQP_:88):
```
sigma_n_total = σ_n0 (+49 MPa static) + σ_n_nuc (0) + sigma_n_trial
```
Friction never modifies σ_n: `sigma_n_corr = sigma_n_trial`.

### 2.2 The decomposition at the seed (qp=480, t=0.4774 s, rank 104, is_shared=1)
```
V_abs      = 12.57 m/s        tau_abs   = 63.2 MPa     sigma_n_tot = +6.0 MPa (compressive)
sn_vjump   = −4.17e7 Pa  (−41.7 MPa)   ← entire collapse
sn_sterm   = −3.28 Pa    (≈ 0)          ← bulk normal stress unchanged
⇒ sigma_n_trial = −41.7 MPa ;  v_n⁻ − v_n⁺ = sn_vjump/η_p = −5.2 m/s
⇒ sigma_n_total = 49 − 41.7 = +6 MPa  ✓
```
This holds through the **entire** runaway: `sn_vjump` tracks `V_abs` (−4e7 → −1e9 →
−1e12), `sn_sterm` stays negligible.

### 2.3 Why a non-zero `[[v_n]]` implicates the frame
The friction slip-rate vector is `V = V1·can_t1 + V2·can_t2` — by construction in the
tangent plane. The imposed velocity jump across the fault is that vector, so
```
[[v_n]] = V·can_n = V1 (can_t1·can_n) + V2 (can_t2·can_n)
```
With an orthonormal frame (`can_t1·can_n = can_t2·can_n = 0`) this is **exactly 0** —
slip is pure shear, no opening (true to machine precision on planar TPV faults).
`[[v_n]] = 5.2 m/s ≠ 0` therefore means tangential slip is leaking into the normal
channel, via one of:
1. `can_n` not ⊥ the local fault surface / frame not orthonormal (`can_t·can_n ≠ 0`);
2. the + and − sides using **different** normals (ghost frame ≠ local frame) on the
   shared face;
3. a genuine physical opening the imposed state fails to constrain (no no-opening BC).

### 2.4 Why (1)/(2) — a geometric leak — and not (3) dynamic opening
- **Instantaneous & rate-proportional:** `[[v_n]]/V ≈ 41%` is already present at the seed
  (δ=0.01 m, essentially unslipped), proportional to the slip *rate*, not built up over
  time. (3) a dynamic opening would develop slowly; this does not.
- **Shared-only** (`is_shared=1`) and **resolution-independent** (Dc2≈Dc8 identical):
  exactly a frame/geometry defect specific to the curvilinear shared-face frame
  construction, independent of `d_c`.
- **`sn_sterm ≈ 0`:** only the velocity channel is contaminated, not the bulk stress —
  consistent with a rotation/projection error injecting velocity, not stress.

A 41% leak ⇒ a **~24° effective `can_n` error**. That is large and must be **measured**,
not assumed (test #1 below). The 41% at V=12.6 may already include a turn or two of the
free-slide feedback (the trace threshold is 10 m/s; the *initial* coupling could be
smaller — test #2 lowers the threshold).

---

## 3. The amplifier — sign-indefinite `[[v_n]] ↔ traction` feedback (NOT tension)

**The amplifier is the positive-feedback loop itself, which runs away with *either* sign
of `[[v_n]]`** (confirmed by the cap run 7747835, see below):

```
slip V ──frame leak──▶ [[v_n]] ──Godunov──▶ σ_n , τ ──friction──▶ V    (loop gain > 1/step)
```
- **`[[v_n]] < 0` (opening, qp=480):** σ_n tensile → `tau_str=0` → `V=|τ|/η_s` (free slide) → runaway.
- **`[[v_n]] > 0` (closing, qp=481):** σ_n → +TPa → `tau_str=μ_d·σ_n` huge, **but the shear
  v-jump grows `tau_abs` even faster** → `V=(tau_abs−tau_str)/η_s` still huge → runaway.

Both branches blow up; the only difference is the sign of the normal-velocity jump. This is
the closed-loop instability the memory flagged for the V1/dip channel (~1.07/step gain). The
`sn_sterm` (bulk stress) term stays ≈0 throughout in **both** signs.

**Cap experiment (7747835) — decisive refutation of "tension is the amplifier":**

| time | quantity | no-cap (7747812) | cap ON (7747835) |
|---|---|---|---|
| t=0.490 | max_slip | 3.114 m | **0.0126 m** (tensile qp=480 arrested) |
| t=0.525 | max_slip | 439.9 m | **441.5 m** (compressive qp=481 took over) |
| t=0.560 | max_slip | — | **5335 m** |

Three runaway QPs: qp=480 `(607517.9,3706359.2,−4543.0)` tensile (cap arrests it);
**qp=481 `(607701.4,3706352.7,−4743.1)` compressive — `sigma_n_total` never ≤ 0, so the cap
never fires; ignites independently even when qp=480 is suppressed, and carries the runaway**;
qp=482 compressive. ⇒ a no-opening / tension fix is **necessary-but-not-sufficient**; the fix
must target the `[[v_n]]` coupling (frame) or the loop gain.

Secondary amplifier (still relevant to *why shared-only*):
1. **Shared-fault step forcing** (the runtime warning: *"2262 of 128205 fault DOFs (1.76%)
   on shared faces see the gradual_overstress as an end-of-macrostep step rather than a
   smooth ramp"*): on shared faces the nucleation shear is applied as a large per-macrostep
   STEP; a big shear kick × the frame leak = a big normal-velocity kick each step. Interior
   faces get the smooth ramp (tiny per-step shear) → tiny leak → bounded. A second reason
   it is shared-only.

---

## 4. The iterator-vs-output decoupling (why `fault.vtkhdf` shows compressive σ_n)

The runaway slip is accumulated **only** in the iterator (`tpv205_substep_iterator.cpp:122`,
the sole `slip += V·dt` site in the code). Its `sigma_n_total` collapses tensile as above.
But the **written** `sigma_n_corr`/`slip_rate` for shared QPs is **overwritten** by the
shared-fault macro-dt solve (`wave_operator.inl:4222`) and the cross-rank reconcile, whose
`sigma_n_trial` stays small/bounded — so `fault.vtkhdf normal_stress` is a flat ~47–52 MPa
compressive (min 30.4 MPa, `n(tensile)=0`) **even at the runaway points**, while slip
explodes to 3e6 m. The output samples the bounded solve; the runaway is driven by the
collapsing one. (This is why the macro/output view never showed the tension — it lives only
in the iterator's per-sub-step trial.)

---

## 5. Refuted hypotheses (with the reasons)

- **R-007d cohesive-zone under-resolution.** D_c 2→8 (L_nuc/h ≈12→48, 4× over-resolved) left
  the runaway essentially unchanged; the **seed is bit-identical** (fires at δ≈0.01 m ≪ d_c,
  before weakening engages). Finer mesh is not the lever.
- **R-005d "our LSW free-slide deviates from SeisSol" (as a friction-law fix).** SeisSol's
  LSW free-slides under tension **identically** (`LinearSlipWeakening.h:153-160` strength
  `= -cohesion - μ·min(σ_n,0)`; `:77-79` slipRate `= max(0,(|τ|-strength)/η_s)`). SeisSol
  does not blow up because its σ_n never collapses — so changing our friction law would
  *deviate* from SeisSol, not fix toward it. The fix is the σ_n trigger, not friction.
- **"Corrupted ghost bulk stress."** The bulk normal-stress term `sn_sterm ≈ 0`, so there is
  no corrupted bulk stress to blame; and Phase 0 (job 7747304) measured the ghost path
  bit-exact cross-rank (~1e-14). The contamination is in the **velocity** channel via the
  frame, not the stress channel via the ghost.
- **"worst_rel=0 ⇒ single-rank cause."** Unsound (conceded): the reconcile *forces* the
  output fields equal, so `worst_rel=0` is guaranteed post-reconcile and cannot distinguish
  a symmetric physical cause from a copied ghost value.

---

## 6. Decisive next tests

1. **`can_n` orthogonality at the seed (the trigger test). — INSTRUMENT SHIPPED.** The
   `[FRAME]` trace (`wave_operator.inl` `EvaluateBulkAtFaultQPsCanonical` shared branch; plan
   `PLAN_frame_orthonormality_diag_2026-05-23.md`) prints, per seed shared QP:
   orthonormality (`n.t1, n.t2, t1.t2 → 0`; `|.|−1 → 0`), `sign_flipped`, `nl`, `can_n`
   (degenerate-band read: `n·ref = −can_n[1]` for SAFS `ref=(0,−1,0)`), and the **global bulk
   velocity-jump decomposition** `dv_n / dv_t1 / dv_t2`, `|dvg|`. `dv_n` ≡ the iterator's
   `sn_vjump/eta_p` (cross-check, unit-tested `test_frame_diag_projection_identity`, rel=0).
   Runtime env-gated (`SEAS_DIAG_FRAME=1`, `_XYZ`, `_R`); byte-exact when off (ADER interior-
   vs-shared equivalence rel=0). **Read:** orthonormal + `dv_n≈|dvg|` ⇒ genuine bulk opening
   (imposed-state locus, not frame); large orthonormality residual or `sign_flipped` flapping
   across seed QPs ⇒ frame/`sign_flipped` defect. Note `fault_basis.hpp` builds the frame
   per-QP & orthonormal by construction, so the residual is *expected* ~1e-15 — the live
   suspects are `sign_flipped`/degenerate-band and a real bulk opening. **Wired default-ON in
   the Dc2 sbatch; awaiting the next Frontera run.**
2. **Lower the trace threshold** (`SEAS_DIAG_SLIP_VTHR=0.1`) to catch the first slip: if
   `[[v_n]]/V ≈ 41%` from V≈0 it is a pure geometric leak; if it grows from smaller, the
   free-slide feedback is inflating it.
3. **Zero the dip channel** (`SEAS_FORCE_V1_ZERO`; dip = `can_t1`): if `[[v_n]]` collapses,
   the leak is specifically **dip-slip → normal** on the 32° fault (B.2), narrowing the fix
   to the dip tangent.
4. **`SEAS_NOOPENING=1` cap run (the amplifier test). — DONE (7747835): REFUTED.** Predicted
   arrest; instead the cap arrested only the tensile QP (qp=480) while a compressive
   runaway (qp=481) carried the same max_slip (441.5 vs 439.9 m at t=0.525). The amplifier is
   the sign-indefinite `[[v_n]]↔traction` feedback, not tension. See §3.
5. **Tighten dt** (smaller macrostep): if `[[v_n]]` shrinks, the shared-fault step forcing is
   a real amplifier (matches the warning's own advice). [still pending]
6. **Measure the per-step loop gain** at qp=481 (compressive): plot `[[v_n]]`, `tau_abs`,
   `sigma_n` vs sub-step index and fit the growth factor. Confirms gain>1 and whether the
   frame fix (test #1) brings it below 1. [new, after #1]

---

## 7. Fix direction (after the tests confirm)

- **Trigger (root):** make the canonical fault frame on the curvilinear shared face
  geometrically correct — `can_n` ⊥ the local fault surface and `{can_n, can_t1, can_t2}`
  orthonormal, consistent across the ± (local/ghost) representations — so tangential slip
  produces `[[v_n]] = 0`. Same `ComputeOrientedFrame`/`FaultBasis` locus as Part A.
- **No-opening is NOT a fix (7747835).** Capping `σ_n ≤ 0` leaves the dominant *compressive*
  runaway (qp=481) untouched — the cap is necessary-but-not-sufficient. Do not rely on it.
- **Forcing (secondary):** apply the gradual_overstress on shared faces as a smooth per-sub-
  step ramp (as on interior faces) rather than an end-of-macrostep step, removing the
  large-shear-kick amplifier.

Because the runaway is **sign-indefinite** (tensile *and* compressive QPs blow up), the only
fix that addresses both is **eliminating the `[[v_n]]` leak at the source** — i.e. the
**frame accuracy on the curvilinear shared fault** so that tangential slip produces
`[[v_n]] = 0` and the loop gain drops below 1. Tension handling is a dead end.

> **Note (2026-05-24): §7 is SUPERSEDED.** The frame-accuracy fix direction is wrong — the
> frame is built identically in the predictor and macro paths (verified: `sign_flipped` at
> `wave_operator.inl:2199` and `:4805`; orthonormal by construction) and cannot produce a
> predictor-vs-macro divergence. The actual fix direction is **R-004** (fix the R-1601 byNODES
> ghost-exchange unpack — CONFIRMED scrambled by `test_ghost_exchange_bynodes_vs_scalar`) then
> **R-001** (restore the single-slip-rate invariant on shared faces). See §8.5 and
> `REVIEW_speckle_seissol_drdg3d_rootcause_2026-05-24.md`. Do not send a fix agent at the frame.

---

## 8. Roadmap: function connections — ours vs SeisSol (2026-05-24)

This section maps how the per-macro-step functions connect in each code and where the two
architectures diverge. It is the reference companion to the R-008 confirmation (top) and the
`/code-debug` static-read + SeisSol comparison.

### 8.1 Our pipeline — `AdvanceADERWithSubStep_Spatial` (`spatial_dyn_driver.cpp:400-461`)

```
 base state Q  (committed solution; bounded — V_max ≈ 5 m/s, σ_n ≈ 46 MPa)
   │
   │ [1]  ComputeADERSubStepStates(Q, dt, O, tau_nodes → Q_per_node)   wave_operator.inl:1302
   │        element-local CK Taylor;  POINT values  Q̃(τ_o)  at sub-step midpoints
   │        ── no neighbor data ──
   ▼
 Q_per_node[o]
   │
   │ [2]  for each o:  EvaluateBulkAtFaultQPsCanonical(Q_per_node[o])  wave_operator.inl:1798
   │        eval at fault QPs; SHARED face: self = local DOFs,
   │        nbr = R-1601 per-sub-step ghost EXCHANGE of Q̃; rotate (sign_flipped), split ±
   │        ►► [FRAME] diag  (dv_n here)
   ▼
 Q_pointwise_plus/minus[o]                          ◄── the state that RUNS AWAY
   │
   │ [3]  iterator.AdvanceWithSubStepStates(Q_pointwise_±)         tpv205_substep_iterator.cpp
   │        per sub-step, per QP: ComputeTrialTraction (Eq 7a) → friction (Brent) → V
   │        slip += V · dt_sub          ◄── SLIP ACCUMULATOR (→ max_slip 4.6e6 m)   :125
   │        BuildImposedState → I_imp_±
   │        ►► [SLIP] diag  (sn_vjump here)
   ▼
 I_imp_plus/minus_flat
   │
   │ [4]  SetSubStepFaultImposedStates(I_imp_±)                    wave_operator.inl:1411
   │
   │ [5]  AdvanceADER(Q, dt, O → Q_new)
   │        ├─ ComputeADERTimeIntegrated(Q → I)   same CK, but the INTEGRAL ∫₀^dt    :1205
   │        ├─ ComputeADERVolumeUpdate(I)
   │        ├─ ComputeADERFaceFluxRHS(I)        [INTERIOR fault QPs]
   │        │     if installed: CONSUME I_imp_± as the flux (ONE solve ✓)        :3871-3885
   │        └─ ComputeADERSharedFaceFluxRHS(I)  [SHARED fault QPs]
   │              "SHARED FALLBACK" — DISCARD I_imp_± (`(void)substep_I_imp_*`),
   │              RE-SOLVE EvaluateADER_LSW on I/dt → flux + WRITTEN DOFData       :4856-4895
   │              ►► [MACRO] diag   ◄── R-001: on shared, μ←δ(predictor V) but flux←V(I/dt) — DECOUPLED
   ▼
 Q_new  =  Q + dt·(volume + fault flux)
```

**Load-bearing fact — R-001 (`REVIEW_speckle_seissol_drdg3d_rootcause_2026-05-24`): on
SHARED faces SAFS runs TWO friction solves on TWO different states, and the slip that weakens
μ comes from the solve whose flux is thrown away.** Step [3]'s iterator solves on the **point
predictor** `Q_pointwise` and is the *sole* slip-accumulation site (`slip += V·dt_sub`,
`tpv205_substep_iterator.cpp:125`). On **interior** faces step [5]'s
`ComputeADERFaceFluxRHS` *consumes* that iterator `I_imp_±` as the flux (3871-3885) → one
solve, slip↔flux coupled, invariant holds. On **shared** faces
`ComputeADERSharedFaceFluxRHS` *discards* it (the literal `(void)substep_I_imp_*` at
4893-4895; "SHARED FALLBACK… always runs the inline ADER closure", 4856-4858) and **re-solves**
`EvaluateADER_LSW` on the time-integral `I/dt` (4867-4892). So on shared QPs `μ ← δ(predictor
V)` while `flux ← V(I/dt)` — the two slip rates are unrelated. SeisSol/drdg3d make this
impossible (one slip rate drives both). The predictor jump diverges 10⁸ from the macro
(3.8e6 vs 0.015 m/s) precisely because the two solves are decoupled and nothing couples the
predictor's runaway back to the bounded applied flux. **This is the root design flaw — not a
side effect of the "two-state" framing but its concrete, shared-face-only locus.**

### 8.2 SeisSol pipeline — `BaseFrictionLaw::evaluate` (`CpuImpl/BaseFrictionLaw.h:49`)

```
 neighbor element time-DOFs ──(MPI: communicated derivatives)──┐
 local element time-DOFs ──────────────────────────────────────┤
                                                                ▼
   spaceTimeInterpolation  →  qInterpolatedPlus/Minus[o]        Kernels/DynamicRupture.cpp:51
      predictor evaluated at TIME-QUADRATURE points, BOTH sides, SAME quadrature
                                                                │
   BaseFrictionLaw::evaluate                                    │   CpuImpl/BaseFrictionLaw.h:49
     ├─ precomputeStressFromQInterpolated(qI±, etaPDamp) → faultStresses (trial, all o)   :72
     ├─ for timeIndex o = 0 … TimeSteps:
     │      updateFrictionAndSlip(faultStresses, o) → slipRate, traction                  :116
     │         slip accumulates here, SAME quadrature  (accumulatedSlip += timeWeight·slipRate)
     └─ postcomputeImposedStateFromNewStress(qI±, timeWeights)                            :158
            → imposedStatePlus/Minus   (TIME-INTEGRATED, friction-corrected)
                                                                │
            ▼  Neighbor/DR kernel:  imposedState  →  bulk flux
```

**One** `qInterpolated` → **one** friction loop → slip, imposed-state flux, *and* output, all
from the same state with the same time-quadrature weights.

### 8.3 The difference, point by point

| Aspect | **Ours** | **SeisSol** |
|---|---|---|
| State the **friction/slip** solve consumes | per-sub-step **point** predictor `Q̃(τ_o)` (`Q_pointwise`) | `qInterpolated[o]` at **time-quadrature** points |
| **Slip** accumulation | `slip += V·dt_sub` from the *point-predictor* solve (iterator) | `accumulatedSlip += timeWeight·slipRate` from the *same* solve |
| State the **flux/output** uses | **interior:** the iterator's `I_imp` (consumed, :3871-3885); **shared:** a *separate* re-solve on `I/dt` (`I_imp` discarded, :4856-4895) | the **same** `qInterpolated` → `imposedState` |
| # friction evals per face/step | **interior: one** (coupled ✓); **shared: two** (decoupled ✗ — R-001) | **one** |
| Neighbor (+) side | per-sub-step **byNODES batched exchange** of the point predictor (R-1601) — **layout bug CONFIRMED, R-004** | communicated **time-DOFs**, evaluated locally on the *same* quadrature |
| Normal radiation term `η_p·Δv_n` | undamped | optional **`etaDamp`** factor (`DRParameters.h:88`, default 1.0 = off) |

### 8.4 Why this produces our symptom and SeisSol's doesn't

- In SeisSol the slip, the flux, and the output are **the same computation** on **one** state.
  If that state's normal-velocity jump ever blew up, *everything* would blow up together —
  slip cannot run away while the flux/output stay bounded. The time-quadrature integration +
  the optional `etaDamp` keep `η_p·Δv_n` controlled.
- In our code, on **shared** faces, the **symptom-bearing slip** (from the iterator's point
  predictor) and the **bounded applied flux/output** (the separate `I/dt` re-solve) are
  computed from two different states whose `I_imp` link is *discarded* (R-001). That decoupling
  is exactly what run 7747889 shows: the point-predictor jump diverges (→ `max_slip` 4.6e6 m
  via `slip += V·dt_sub`) while the `I/dt` path keeps `V_max`/σ_n/the written output bounded and
  **hides it**. The seed of the predictor divergence is now **confirmed**: the R-1601
  `vdim=NUM_STATE` byNODES `FaceNbrData()` unpack `src[c·n_fn+j]` reads **scrambled neighbour
  components** (`test_ghost_exchange_bynodes_vs_scalar`, np=2, 918 mismatched reads/rank;
  e.g. the `(c=0,j=27)` slot returns `encode(c=1,i=108)` instead of `encode(c=0,i=135)`). The
  macro path's per-component scalar exchange is unaffected — which is why predictor-self ≈
  macro-self but predictor-neighbour diverges.

### 8.5 Fix direction (root cause located; order matters)

1. **R-004 — FIXED (2026-05-24).** The R-1601 unpack now reads the neighbour values through
   `pfes_full_state_->GetFaceNbrElementVDofs(nbr_idx, nbr_vdofs)` →
   `nbr_all[nbr_vdofs[c·ndof_per_el_ + i]]` (byNODES vdof ordering), the layout-agnostic MFEM
   pattern (cf. `elasticity_operator_debug.inl:629-638`) — replacing the wrong component-slab
   formula `src[c·n_fn+j]`. The batched R-1601 exchange is preserved (no per-component-exchange
   perf regression). Guarded by `test_ghost_exchange_bynodes_vs_scalar` (now **GREEN**, in
   `test-v92-regression-gates`): the `GetFaceNbrElementVDofs` read is bit-identical to the
   per-component scalar exchange at every (element, component, dof); the old flat-slab formula
   got ~918 reads/rank wrong. Verified locally (np=2): R-004 guard PASS, shared-fault reconcile
   cross-rank 4/4 PASS (incl. the SUBSTEP/SAFS production path, cross-rank bit-identical),
   interior-vs-shared-branch-live PASS, SAFS driver compiles. **Next: a Frontera SAFS Dc2 run to
   confirm the runaway is gone** before doing R-001. (R1600 + TPV205-crossing not run locally —
   pre-existing macOS `dyld` / missing-mesh-file env issues, unrelated to the fix.)
2. **R-001 — restore the single-slip-rate invariant on shared faces.** Either consume the
   iterator's `I_imp` on shared QPs like interior (only safe **after** R-004, else the now-
   bounded-once-fixed predictor flux can be fed back), or — the lower-risk interim —
   accumulate shared-QP slip from the *applied* macro V rather than the predictor V.
   **Order: R-004 before R-001** — fixing R-001's "consume `I_imp`" variant while the predictor
   is still scrambled would feed the runaway straight into the bulk (the np=10 `tau=4e28 →
   SIGABRT` the R-1601 fallback was added to prevent).
3. **`etaDamp` analog** on `η_p·Δv_n` — optional defense-in-depth (SeisSol's knob), not a
   substitute for 1–2.
