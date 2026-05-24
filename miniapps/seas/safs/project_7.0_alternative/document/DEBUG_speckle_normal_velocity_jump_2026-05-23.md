# Debug findings: SAFS slip "speckle" runaway → normal-velocity-jump (opening) on the curvilinear shared fault — 2026-05-23

## TL;DR (root cause localized)

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
