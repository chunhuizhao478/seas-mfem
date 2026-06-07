# Debug plan — spurious DIP-SLIP DRIFT in TPV102 & TPV31 (cross-benchmark)

Date: 2026-06-06
Author: Claude (code-plan)
Worktree/branch: `worktree-fix-dip-drift` (reset onto `b09dbfc`, `system/spatial_dyn_driver`)
Status: **DIAGNOSIS + EXPERIMENT PLAN. No code changed.** Verification on the
production mesh requires Frontera A/Bs (ask before launching, per project memory).

> **Relationship to prior docs.** This plan supersedes the root-cause *conclusion*
> of `tpv31_debug_document/tpv31_dip_slip_residual_findings_2026-06-06.md` and
> `..._code_review_2026-06-06.md` (which declared the dip path "PASS / numerical
> polish, ~1e-4 m") and **re-opens** the earlier
> `tpv104_debug_document/tpv104_sigma_n_leak_root_cause_2026-04-25.md`
> rank-2-stress-rotation hypothesis. The two prior efforts reached **opposite**
> conclusions; this plan reconciles them with new evidence (the along-strike
> antisymmetric growth they both missed) and a decisive experiment neither ran
> (freeze the per-face frame).

---

## 0. TL;DR

- The dip-slip drift is **real, structural, and at the strike edges it is large**
  (~0.075–0.117 m, i.e. **2–6 % of strike slip**), not the ~1e-4 m "polish"
  the June docs reported — they only inspected strike-center stations.
- It is **order-, flux-, and integrator-independent** ⇒ not discretization, not
  the time stepper, not the flux scheme. The defect is in code **common to every
  fault-interface path**.
- New, sharp fingerprint (this doc): the dip slip is **anti-symmetric in
  along-strike position, grows ≈ linearly with along-strike distance, peaks at
  the free surface, ≈ 0 on the symmetry plane**, and is driven by a **persistent
  (non-oscillatory) dip slip-RATE** that the dissipative friction integrates.
- Two prior root-cause theories conflict: **(A) rank-2 stress-rotation contamination
  from a per-QP fault-frame sign-flip** (April) vs **(B) mesh-asymmetry `[[v_dip]]`
  aliasing** (June). The new fingerprint **favors (A)**, but the June unit test
  does **not** actually exonerate the frame (see §4.3). One cheap experiment
  discriminates them and *is itself the candidate fix*.
- The fix, if (A), lives entirely in **editable `dynamic/` code** (precompute one
  sign-normalized rotation per face, reuse for all QPs & both sides) — `fault/
  fault_basis.hpp` (no-touch) need not change.

---

## 1. Symptom — consolidated from the gold results

Gold runs inspected (this session):
`tpv102/gold/{upwind_p1_aderO2,upwind_p2_aderO3,mixedflux_p2_rk45}_…`,
`tpv31/gold/{upwind_p1_aderO2,mixedflux_p1_rk45}_…`.
Reference overlay = DRDG3D (TPV102) / SeisSol (TPV31), both ≈ 0 dip slip.

### 1.1 Magnitude & spatial structure (the load-bearing new evidence)

TPV102 p2 upwind, **final-time dip slip vs strike slip** (m):

| station (strike, depth) km | slip_dip | slip_strike | dip/strike |
|---|---|---|---|
| (0, 3)    | −1.3e-3 | 8.21 | −1.6e-4 |
| (+12, 3)  | **−0.117** | 5.87 | −2.0e-2 |
| (−12, 3)  | **+0.117** | 5.88 | +2.0e-2 |
| (+9, 7.5) | −1.95e-2 | 6.27 | −3.1e-3 |
| (−9, 7.5) | +2.03e-2 | 6.28 | +3.2e-3 |

⇒ **anti-symmetric across the hypocenter (±x → ∓ sign, equal magnitude)** and
`|slip_dip|` grows with `|x_strike|` and toward shallow depth.

TPV31 p1 upwind, dip slip vs along-strike station (surface row, dp000):
`st000 ≈ +3.8e-4` → `st060 ≈ −3.9e-2` → `st120 ≈ −7.5e-2` (≈ linear, 2×), and
within each strike column the magnitude is **largest at the free surface
(dp000), decaying with depth**.

### 1.2 Invariances (what is exonerated)

| Varied | TPV102 result | Conclusion |
|---|---|---|
| order p1 vs p2 | dip drift ≈ −1.0e-3 vs −1.3e-3 at (0,3) — same scale | **not** resolution |
| integrator/flux: upwind ADER vs mixed-flux RK45 | both drift (−1.3e-3 vs −0.9e-3 at (0,3)) | **not** time stepper, **not** flux scheme |
| mesh: deliberately symmetric mesh | drift persists | **not** random mesh asymmetry alone |

### 1.3 Temporal character

Dip slip-rate has a small **persistent, single-signed** bias (not zero-mean
chatter): see TPV102 `flt_0_3` — `V_dip` holds O(1e-5 m/s) for seconds; `slip_dip`
ratchets monotonically. In **TPV102 (rate-state, always creeping)** it drifts
secularly to end of run; in **TPV31 (LSW, re-locks)** it **plateaus** after the
front passes. Both are the signature of a dissipative friction integrating a
biased dip traction (no reversal term).

---

## 2. Where the data flows (verified file:line)

Fault-local frame convention (CLAUDE.md, BP5/Tandem): **`t1 = dip` (component 1),
`t2 = strike` (component 2)**; for the vertical y=0 fault `can_t1=(0,0,−1)`,
`can_t2=(+1,0,0)`, `can_n=(0,−1,0)`.

1. **Per-face/per-QP frame** built in `fault/fault_basis.hpp`
   `FaultBasis::ComputeOrientedFrame` (≈411–512): `strike=normalize(up×n)`,
   `dip=normalize(strike×n)`, `t1=dip, t2=strike`, with a sign-canonicalization
   step (`NormalNeedsFlipToCanonical`, ≈362–379) that may **negate all three
   vectors**. `n` comes from `CalcOrtho(Face->Jacobian())` **per QP**. **No-touch.**
2. **Rotation matrices** `GodunovFlux::BuildRotation/Inverse`
   (`dynamic/godunov_flux.cpp:221–301`): 9×9 with rank-1 velocity block (`Q,Qᵀ`)
   and **rank-2 Voigt stress block**, from `Q=[n;t1;t2]`. Orthonormal **iff**
   `(n,t1,t2)` is. **Editable.**
3. **Trial traction** `FaultFaceFlux::ComputeTrialTraction`
   (`dynamic/fault_face_flux.cpp:49–85`): all three channels from the velocity
   jump, identical formula. Dip: `tau1_trial = eta_s·([[v_dip]] + …)`. For
   TPV102/TPV31, `tau1_0 = tau1_nuc = 0` ⇒ `tau1_total ≡ tau1_trial`. **Editable.**
4. **Friction decomposition (V ∥ tau)** `fault_face_flux.cpp:261–262`:
   `V1 = V_abs·tau1_total/(strength+eta_s·V_abs)`, likewise V2. **Correct &
   intentional** (CLAUDE.md: parallel, not anti-parallel). Identical to DRDG3D
   (`mod_wave.F90:1526/1565`) and SeisSol (`LinearSlipWeakening.h:83`,
   `RateAndState.h:229`). **Not the bug.**
5. **Closed-loop dip instability** `fault_face_flux.cpp:238–257` (in-code comment +
   `SEAS_FORCE_V1_ZERO` knob): for pure strike-slip the dip channel has no
   restoring traction; a 1-ULP `tau1_total` perturbation grows **~1.07/step**,
   reaching cm-scale dip slip in ~1 s. `tau1_corr` is fed back into the bulk via
   `BuildImposedState` (`:288–307`). **Editable. This is the amplifier, not the
   seed.**
6. **Imposed-state rotation back to global** at the `wave_operator.inl` call sites
   (interior-RK ≈2896, shared ≈3567, interior-ADER ≈2510/4071): reconstruct
   `can_n/can_t1/can_t2` from **per-QP** `FaultBasisQPData`, apply
   `T_can·Q`. **This is the rank-2 rotation April flagged.** **Editable.**
7. **Slip accumulation** `dynamic/friction_substep_iterator.cpp:143` (RS) /
   `:239` (LSW) / `tpv102_substep_iterator.cpp:182`: `slip1 += V1·dt` (plain
   forward sum, secular). **Station output** reads `slip1`(dip)/`slip2`(strike)
   in the same frame (`tpv102_stations.hpp:204`, `tpv31_stations.hpp:13`). No
   output-basis mismatch. **Editable, correct.**

**Reference invariant both DRDG3D and SeisSol satisfy** (Agent-C): the prestress
and the solver use the **same** fault basis, and the basis is built so a planar
fault has the **same orthonormal frame at every node**; true strike/dip is used
only for I/O. Because their resolved `tau_dip` is genuinely ~0, the (identical)
V∥tau decomposition produces exactly 0 dip slip and the closed loop has nothing
to amplify.

---

## 3. Hypotheses (ranked)

| ID | Hypothesis | Predicts the fingerprint? | Location | Status |
|---|---|---|---|---|
| **H1** | **Per-QP fault-frame sign/orientation inconsistency** leaks strike-shear into the dip channel through the **rank-2 stress rotation** (rank-1 velocity self-cancels, rank-2 does not). Seed sign tracks the strike-shear sign ⇒ anti-symmetric in x; amplitude tracks strike-shear ⇒ grows with rupture extent; worst where field is strong & mesh irregular ⇒ free surface. | **Yes — best match** incl. anti-symmetry + linear-in-strike + surface peak | `fault_basis.hpp` (no-touch math) + `wave_operator.inl` per-QP frame use (editable) + `godunov_flux.cpp` rank-2 block (editable) | **Prime suspect** (April-confirmed; June did not re-audit rank-2) |
| **H2** | **Closed-loop dip amplification** (`fault_face_flux.cpp:238–257`): a tiny `tau1_total` seed (from H1 or H3) is grown ~1.07/step by the friction reaction → bulk feedback. | Explains order-independence (feedback gain, not resolution) & edge magnitude (more growth steps + larger strike field) | `fault_face_flux.cpp` (editable) | **Amplifier, not seed.** References lack it only because their seed ≈ 0. |
| **H3** | **Mesh-asymmetry / near-fault upwind-dissipation `[[v_dip]]`** aliasing seed (June theory). | Explains surface peak & p-independence; **does NOT naturally explain a clean anti-symmetric linear-in-strike ramp on a symmetric mesh** | wave op / mesh / flux (editable) | **Contributing at most; demoted** by §1.1 |
| **H4** | MPI fault-plane-cut / FP-non-associativity seam (`tpv31_sigma_n_drift_findings_2026-06-02`). | random/partition-dependent sign, vanishes at np=1 — **inconsistent** with the clean structured serial pattern | driver | **Separate co-existing bug; out of scope** (confirm via np=1) |
| H5 | Setup prestress dip leak | TPV102/31 seed dip components directly as 0 (`tpv102_setup.hpp:70`) | — | **Excluded** |
| H6 | Friction regularization/floor/plate-rate injects dip bias | no per-channel floor exists; floors act on scalar magnitude only | — | **Excluded** |

---

## 4. Why the June "frame is clean" verdict is not decisive

1. **The June code review (R-1) read the trial-traction & V∥tau formulas and
   declared them correct — true, but those are not where H1 lives.** H1 is in the
   **rank-2 stress rotation of the imposed state** (`BuildRotation` Voigt block +
   per-QP `T_can` at the `wave_operator.inl` sites), which R-1 did not audit.
   April explicitly identified this as the un-audited prime suspect.
2. **"A frame bug would be depth-uniform & method-independent" is false for a
   *per-QP sign-flip* frame bug.** The leak is modulated by the rotated stress
   field, which is itself surface-peaked and method-dependent — so a frame leak
   *also* produces a surface-concentrated, method-dependent artifact. The June
   symptom-shape argument does not separate H1 from H3.
3. **The June unit test (`[[v_dip]]` survives a y-mirror mesh, ratio ~0.9) does
   not exonerate the frame.** A per-QP sign-flip frame leak is a z-tangential /
   orientation effect that is *also* immune to a y-mirror — so "survives y-mirror"
   is consistent with **both** H1 and H3. The test that *would* separate them
   (freeze the per-face frame) was never run.
4. **The June magnitude claim ("~1e-4 m, polish") came from strike-center
   stations only.** §1.1 shows ~0.1 m at the strike edges. This alone reclassifies
   the issue from "polish" to "correctness."

---

## 5. Diagnostic plan (cheap → expensive; each step is decisive)

Guiding constraints (project memory):
- **No local production-mesh TPV102/TPV31 runs.** Local **unit tests on tiny
  fixtures only**. Frontera runs need explicit approval.
- **No-touch:** `fault/`, `friction/dieterich_ruina.hpp`, `bp5/bp1/bp2/domain/
  solver`. All fixes go in **`dynamic/`** (duplicate any needed BP5 utility).
- Keep TPV102/TPV104 **byte-exact regression contract** in mind for any change.

### Phase 0 — instrument & reproduce at unit scale (LOCAL, no approval)

**P0.1 Frame-consistency probe (decisive for H1).** In
`tests/unit/test_fault_planar_serial.cpp` (already extended with `[[v_dip]]`
probes), add per-fault-QP logging of: `can_n, can_t1, can_t2`, `sign_flipped`,
`elem1_on_plus`, and the scalars `can_t1·strike_phys`, `can_t1·can_t2`,
`|can_t1|−1`. **Pass:** dip axis has zero strike projection, is identical across
all QPs of a face and (up to the documented global sign) identical on both sides.
**Fail (expected if H1):** `can_t1` carries a per-QP strike component / flips sign
across QPs. (`SEAS_DIAG_FAULT_FLUX` `[C-1s INT BASIS]` tracer at
`wave_operator.inl:~2986` already exists — turn it on.)

**P0.2 Freeze-frame discriminator (THE decisive test — H1 vs H3).** Add a
build/env-gated override at the `can_*` reconstruction sites in
`dynamic/wave_operator.inl` forcing the **exact canonical orthonormal frame**
`n=(0,−1,0), t1=(0,0,−1), t2=(1,0,0)` for **all** fault QPs & both sides (or,
more generally, **one sign-normalized `T_can` per face reused for all QPs**).
Re-run the planar serial harness, asymmetric mesh.
  - **dip leak `[[v_dip]]`/imposed-stress-dip collapses → H1 confirmed** (rank-2
    frame), April vindicated, June refuted. The override *is* the fix prototype.
  - **dip leak persists → H3** (mesh-asymmetry velocity-jump seed); pivot to §6-B.

**P0.3 Rank-2 leak microtest (confirms the mechanism of H1).** Extend
`tests/unit/test_canonical_rotation_pure_strikeslip.cpp` (currently only tests the
*ideal* hardcoded frame) to feed **per-QP sign-perturbed** frames and assert the
dip channel stays zero. Expected: a single-axis sign flip leaves the rank-1
velocity rotation clean but injects an off-diagonal dip term in the rank-2 stress
rotation — reproducing H1 at unit scale and giving a permanent regression guard.

**P0.4 Closed-loop characterization (H2).** With the `SEAS_DIAG_V1_DRIFT` tracer,
log `tau1_trial, tau1_total, V1` per substep at a `+x` and `−x` fault QP in the
harness. **Confirm** (a) `V1/V2 == tau1_total/tau2_total` (friction faithful),
(b) `tau1_trial(+x) ≈ −tau1_trial(−x)` (anti-symmetric seed), (c) step-over-step
growth ratio of `|tau1_total|` ≈ 1.07 during the active phase (the documented
instability gain).

### Phase 1 — confirm on the real benchmark (FRONTERA A/B — ask first)

Only after P0.2 identifies the seed. Minimal pair, short `tfinal`:

**P1.1 Freeze-frame A/B on TPV102** (vertical fault, cleanest): baseline vs
freeze-frame (P0.2 override) build, np=1 **and** np≥4. Station diff at the
edge stations `flt_±12_3`. **Pass:** `|slip_dip|` at edges drops by ≥1–2 orders
toward DRDG3D (~0). The np=1 leg also settles H4 (if drift persists at np=1, the
MPI seam is not the cause here).

**P1.2 Freeze-frame A/B on TPV31** (dipping/strike-slip with free surface):
same, checking the along-strike ramp `st000→st120` collapses and the surface
peak (dp000) flattens.

(If P0.2 instead pointed to **H3**, replace P1 with the §6-B cure A/Bs:
central/mixed-flux-near-fault and full over-int+resample on pure-upwind.)

### Phase 2 — fix, then regression

Apply the §6 fix matching the confirmed seed; re-run P0 unit tests (must pass +
new guards), then a Frontera validation matrix (TPV102 p1/p2, TPV31 p1/p2,
upwind + mixed-flux) confirming dip parity AND **no regression** in strike slip,
strike stress, σ_n, or the TPV102/104 byte-exact contract.

---

## 6. Candidate fixes

### 6-A. If H1 (frame) — the principled fix (preferred; matches April + references)

Precompute **one orthonormal rotation `T_can` per fault face**, after
sign-normalization, and **reuse it for all QPs of that face and both sides** —
instead of reconstructing a per-QP frame from `CalcOrtho` at every quadrature
point. This removes the per-QP sign-flip bimodality that pollutes the rank-2
stress rotation, and reproduces the reference invariant ("same orthonormal frame
at every node of a planar fault").

- **Location: editable `dynamic/` only.** Do it at the frame-reconstruction sites
  in `wave_operator.inl` and/or by caching a face-level `T_can` in the
  `FaultBasisQPData`/`precomputed_face_fluxes` path. **Do not edit
  `fault/fault_basis.hpp`** (no-touch) — consume its per-face result once,
  freeze it.
- **Risk control:** must preserve the TPV102/TPV104 byte-exact regressions on a
  mesh where the per-QP frame was already consistent (the freeze is a no-op
  there). Verify with `make test` + the byte-exact targets before any Frontera run.
- This is exactly P0.2's override promoted to production (cache, not hardcode).

### 6-B. If H3 (mesh/aliasing) — the cure levers (June's list, corrected)

Per the June §9 local proof, **over-integration ALONE worsens dip** and is
**blocked on the mixed-flux path** (R-2). Ranked:
1. **Central/mixed flux near the fault** (removes upwind-dissipation seed; mixed-
   flux TPV31 already shows smaller surface dip). Lowest risk.
2. **Full dealiasing = over-int + resample** on the **pure-upwind** path (must use
   *both* knobs; resample is identity without over-int). Needs either a Frontera
   A/B or extending the harness to route through the substep iterator so resample
   is live locally (R-4 test-gap).
3. Near-surface mesh quality (impractical to make z-symmetric with a free surface).

### 6-C. H2 amplifier (defensive, secondary)

The `SEAS_FORCE_V1_ZERO` antiplane lock breaks the dip closed loop but is **only
valid for pure strike-slip** and is wrong for dipping/oblique faults — **not** a
general fix. If the seed is killed (6-A), the loop has nothing to amplify and this
is unnecessary. Consider only a *principled* loop-stabilization (e.g. damping the
dip-channel reaction) if a residual seed is unavoidable, and gate it so TPV205/
TPV104 oblique cases are unaffected.

---

## 7. Acceptance criteria

1. **Dip parity:** at the strike-edge stations (TPV102 `flt_±12_3`,
   TPV31 `st120dp000`), `|slip_dip|` drops by ≥1–2 orders of magnitude and
   `|slip_dip|/|slip_strike| ≲ 1e-3` (toward the DRDG3D/SeisSol ~0).
2. **Anti-symmetry & ramp gone:** `slip_dip(+x) ≈ −slip_dip(−x) ≈ 0`; no
   linear-in-strike growth; no free-surface peak beyond the reference's small
   *physical* dip slip at depth.
3. **No regression:** strike slip, strike shear stress, σ_n unchanged within the
   existing tolerances; **TPV102/TPV104 byte-exact regression targets still pass**;
   `make test` green.
4. **New permanent guards:** P0.1 (frame consistency) and P0.3 (rank-2
   sign-flip → zero dip) as unit tests; the `[[v_dip]]` guard in
   `test_fault_planar_serial.cpp` tightened post-fix.

---

## 8. Open questions / risks

- **Does P0.2 actually collapse the unit-scale dip leak?** If not, H1 is wrong and
  the plan pivots to 6-B with no wasted Frontera time — that is the point of doing
  it locally first.
- **Byte-exact contract vs the freeze:** confirm the per-QP→per-face freeze is a
  no-op on the TPV102/104 regression meshes (it should be, if those frames were
  already consistent). If it changes bytes, gate the new path behind a flag and
  make the reference re-baseline an explicit, reviewed step.
- **H3 may co-exist:** even if H1 dominates the large edge drift, a smaller H3
  residual (the genuine ~1e-4 m the June docs measured at strike center) may
  remain; acceptance criterion #1 is set at the edge where H1 dominates, with #2
  catching residual structure.
- **H4 (MPI seam) is a separate bug**: the np=1 leg of P1.1 keeps it from
  confounding the dip diagnosis.

---

## 9. Evidence index

- Gold data (this session): `tpv102/gold/*`, `tpv31/gold/*` station `.dat` +
  `plots/*` (dip drift in 3rd-row-left "Slip Dip (m)" panels; anti-symmetric edge
  table §1.1).
- V∥tau decomposition + closed-loop comment + `SEAS_FORCE_V1_ZERO`:
  `dynamic/fault_face_flux.cpp:238–267`. Imposed state: `:271–308`.
- Trial traction: `dynamic/fault_face_flux.cpp:49–85`.
- Per-QP frame use / rank-2 rotation: `dynamic/wave_operator.inl` (≈349–460 ref/up
  + per-QP warning; ≈2884–2967, 3560–3568, 2510/4071 rotation sites);
  `dynamic/godunov_flux.cpp:221–325`.
- Frame math (no-touch): `fault/fault_basis.hpp:362–512`.
- Slip accumulation/output: `dynamic/friction_substep_iterator.cpp:143/239`,
  `dynamic/tpv102_stations.hpp:204`, `dynamic/tpv31_stations.hpp:13`.
- Prior docs: `tpv104_debug_document/tpv104_sigma_n_leak_root_cause_2026-04-25.md`
  (H1 origin), `..._dip_dataflow_per_step.md`, `..._bulk_asymmetry_followup_*`;
  `tpv31_debug_document/tpv31_dip_slip_residual_findings_2026-06-06.md` (+ code
  review) (H3), `tpv31_sigma_n_drift_findings_2026-06-02.md` (H4).
- Reference invariant: DRDG3D `mod_wave.F90:1526/1565`, `mod_init_fault.F90:230/263`;
  SeisSol `LinearSlipWeakening.h:83`, `RateAndState.h:229`,
  `BaseDRInitializer.cpp:249`, `ReceiverBasedOutput.cpp:220`.
- Local sources: SeisSol `/Users/chunhuizhao/projects/SeisSol`,
  DRDG3D `/Users/chunhuizhao/projects/drdg3d`.
</content>
</invoke>
