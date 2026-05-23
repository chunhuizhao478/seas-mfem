# Implementation Plan: GENERAL fault tangential-traction stability + cross-rank consistency

> ## ⚠ REVISION 2 (2026-05-23) — root cause CONFIRMED on Frontera; plan pivots
> The `SEAS_DIAG_XRANK` trace on the real Dc2 mesh (job 7747036) **refuted the
> weak-channel/covariant-decomposition hypothesis** (and H-A/H-C). Confirmed
> mechanism: a **redundant-computation + slip-onset-kink cross-rank desync**.
>
> At the rupture front (t=0.478 s) the two ranks sharing the diverging QP feed
> `EvaluateADER_LSW` **bit-identical-to-11-digits** inputs (`can_t1` identical &
> static; `tau1_0` identical; `Iself(A)==Inbr(B)` ⇒ same `(Q_plus,Q_minus)` —
> side fix correct), yet the **outputs split**: one rank slips (V1=−5e−3), the
> other stays locked (V1=0). The QP sits exactly at `|τ|≈μ_s·σ_n` (the LSW
> `V_abs=max(0,(|τ|−τ_str)/η_s)` **kink**); the ~1e−14 cross-rank difference in the
> (non-bit-exact, because the same physical QP is interpolated through two
> different element parametrizations — `shape1∘Loc1` vs `shape2∘Loc2`) bulk Q tips
> the two ranks onto opposite slip/lock branches → full desync → V1 → O(1).
>
> **Root cause:** R-701 (`wave_operator.inl:3026-3043`) deleted the v5 R-501
> owner-broadcast on the premise that the canonical frame makes both ranks'
> `Evaluate` inputs bit-identical. The *frame* is bit-identical; the *bulk Q* is
> not. So the premise is false for long LSW runs whose front crawls through the
> kink (SAFS `gradual_overstress`). Benchmarks miss it: TPV102/104 rate-state is
> smooth (no kink); TPV205's front sweeps fast.
>
> **FIX (this revision):** make the shared-fault DOFData **single-valued across
> the two ranks** — restore an owner-broadcast / reconcile in the ADER LSW
> shared-fault path (see "## REVISED FIX DESIGN" below). The covariant-
> decomposition Phase 3 (rotational covariance / rake sweep) is **DROPPED** — it
> treated a symptom. Phase 1's rake-sweep test stays as a regression oracle.
> Everything from "## Overview" down is the SUPERSEDED original plan, retained
> for the record.

## REVISED FIX DESIGN (owner-broadcast / reconcile — supersedes Phases 2–3)

### Constraints (unchanged + new)
- **Byte-exact for rate-state (TPV102/104) and TPV205**: the reconcile must be
  gated to the LSW shared-fault path used by SAFS, and must be a no-op /
  bit-identical for the rate-state path and for any run where the two ranks
  already agree. Gate: full TPV/BP5 regression `worst_rel ≤ 1e-13`.
- **MPI-collective-safe (R-1600 lesson)**: every rank with shared fault faces
  participates; ranks without must not deadlock. Mirror the existing
  point-to-point face-neighbour topology (as `ExchangeFaceNbrData` /
  `GetSharedFaceTransformations` use), not an unguarded global collective.
- **R-701 revert is justified** by the XRANK evidence (its bit-identical-inputs
  premise is false); cite job 7747036 + `wave_operator.inl:3026-3043`.

### Phase A — confirm the seed source is interpolation, not a real input gap (DIAGNOSTIC, ~done)
The XRANK trace already shows `Iself(A)==Inbr(B)` to 11 digits and the split is
at the kink. One more cheap check before coding: in the diagnostic, also dump
`I_self_can`/`I_nbr_can` to **17 digits** (`%.17e`) at the onset QP to confirm the
two ranks differ at ~1e−14 (interpolation FP), not at a larger scale (which would
indicate a real exchange bug to fix instead). If the diff is ≤ a few ULP·|I|,
proceed to Phase B. (Optional; the mechanism is already clear.)

### Phase B — owner-broadcast/reconcile of shared-fault DOFData (THE FIX)
**Goal:** after the LSW friction solve, the two ranks sharing a fault QP hold
**bit-identical** DOFData (V1,V2,tau1_corr,tau2_corr,sigma_n_corr,slip1,slip2) and
assemble their side's flux from a consistent imposed state, so the slip/lock
decision is made ONCE (by the owner) and copied — immune to the ~1e−14 input
difference.

**Owner rule:** the rank with `elem1_on_plus==true` (exactly one per shared face,
guaranteed by the side fix). Deterministic, already computed.

**Mechanism (`ComputeADERSharedFaceFluxRHS`, LSW branch only):**
1. Both ranks compute `EvaluateADER_LSW` as today (gives each rank's DOFData +
   `I_imp_plus/minus`).
2. **Exchange** the post-solve DOFData fields for every shared fault QP across the
   shared-face neighbour topology (one packed buffer per neighbour rank,
   non-blocking send/recv — same pattern/topology as the predictor-Q ghost
   exchange the substep path already does each sub-step; reuse
   `pmesh.GetSharedFace*` / the R-101 face-key pairing so each QP is matched to
   its peer).
3. **Reconcile:** the non-owner overwrites its DOFData with the owner's
   (received) values. The owner keeps its own.
4. **Imposed state:** broadcast the owner's `I_imp_plus/minus` (canonical-frame)
   alongside the DOFData; both ranks rotate the OWNER's imposed state to global
   (`T_can`, identical on both) for their own side's assembly. This makes the
   assembled fault RHS consistent too (not just the stored DOFData).
5. Gate the whole exchange on `fault_friction_law_ ∈ {LSW, LSW_ForcedRupture}`;
   rate-state skips it (byte-exact TPV102/104).

**Open implementation choices to settle in design review (do NOT guess):**
- (i) Reconcile every sub-step vs only when near the kink. Default: every
  sub-step on shared fault QPs (correctness over micro-opt; the QP set is small).
- (ii) Reuse the deleted v5 R-501 packing code (git history) vs new minimal
  pack/unpack keyed on the R-101 face-vertex-key. Prefer new + the existing
  face-key matcher (the R-101 verify already pairs shared QPs across ranks).
- (iii) Whether to broadcast `I_imp` (step 4) or have the non-owner rebuild it
  from the reconciled DOFData + its local `(Q_plus,Q_minus)`. Broadcasting is
  fully consistent; rebuilding leaves a benign ~1e−14 in the flux (but the
  DOFData — what R-101 checks — is then identical). Decide via Phase A's ULP
  measurement + a regression check.

### Phase C — verify
- The `xrank`-instrumented Dc2 run: at the onset QP, both ranks now show
  **identical** V1/V2/tau*_corr through and past t=0.478 s; R-101 no longer trips;
  the run advances past t=1.0 s (no blow-up).
- TPV102/104/205 + BP5 byte-exact regression `worst_rel ≤ 1e-13` (rate-state path
  untouched; LSW path identical when ranks already agree).
- `test_rupture_rake_sweep_cross_rank` (the existing oracle) + the local
  fault/TPV suite green.
- A NEW small unit test: a 2-tet np=2 LSW fixture driven exactly to the slip-onset
  threshold with an INJECTED 1-ULP cross-rank perturbation on one side's bulk Q;
  assert the reconciled DOFData is bit-identical across ranks (RED before Phase B,
  GREEN after). This is the local oracle the clean fixture CAN provide (injection
  makes the seed; the reconcile must absorb it).

### Risks
- **Deadlock / collective mismatch** (R-1600): the new exchange must use the same
  neighbour topology + non-blocking pattern; test at np=2,4,8 locally.
- **Perf**: per-sub-step fault exchange adds comm; bounded by #shared-fault QPs
  (small). Measure; if hot, restrict to sub-steps where any owned shared QP is
  within a small band of the kink.
- **Byte-exact regression**: gate strictly to LSW; verify rate-state unchanged.

---

# (SUPERSEDED) original plan — covariant decomposition

## Overview
The SAFS run (commit 800b328, side-fix landed) no longer blows up from the
71.6° side-flip, but still diverges (would blow up ~t=1.0 s; the R-101 guard
aborts early at t≈0.455 s) on a cross-rank divergence of `V1` (dip slip-rate) at
a shared fault QP. Investigation shows this is **not dip-specific and not
SAFS-specific**: it is the documented closed-loop instability of the **weak
tangential channel** (the slip component small relative to `|τ|`) in the LSW
friction update (`fault_face_flux.cpp:211-218`). At SAFS's rake the dip channel
is weak; under a dip-dominated load the *strike* channel would be the unstable
one. **This plan delivers a GENERAL fix: a fault tangential-traction update that
is rotationally covariant in the fault tangent plane — stable and
cross-rank-consistent for every rake angle and every (η_s, σ_n, |τ|, dt)
regime — with no privileged channel, no rake threshold, no per-config constant,
and no env-var gate.** It reduces bit-exactly to the current code on pure
strike-slip (TPV/BP5 byte-exact contract). Generality is *demonstrated*, not
assumed, by a rake-angle sweep test (the analog of the tilted-θ sweep that
fixed the side bug).

The prestress/nucleation/frame are already cross-rank-bit-identical (the
projection `τ1 = t1·S·n` is invariant under the Step-5 simultaneous negation of
`(n,t1,t2)`; the runtime flux reconstructs a rank-independent canonical frame).
So the work is the **stability** of the coupled discrete map, not a consistency
patch. **Phase 1 (rake-sweep reproducer) must confirm this and quantify the
loop gain vs rake before Phase 3 implements the fix.**

## The generality principle (the contract for "general, not problem-specific")

The fault friction update acts on the **2-vector tangential traction**
`τ_t = (τ1, τ2)` in the fault tangent plane (τ1=dip, τ2=strike, but the physics
must not privilege either). Require **rotational covariance**: for any in-plane
rotation `R_φ`, solving the friction problem with loading `R_φ·τ_t` must yield
`(R_φ·V_t, R_φ·τ_corr,t)` with **identical stability** — the spectral radius of
the per-sub-step map must be independent of φ (the rake angle). Consequences the
plan enforces as hard requirements:

- **No privileged channel.** Any code that treats τ1 differently from τ2 (other
  than the orthonormal basis vectors themselves) is a bug. Stability must be
  symmetric under the relabeling (τ1,V1)↔(τ2,V2).
- **No rake-dependent branch/threshold**, no "small-dip" special case, no
  per-config constant, no env-var (`SEAS_FORCE_V1_ZERO` must remain a no-op for
  SAFS and is rejected as a fix).
- **Scale-invariance.** Stability is a function of dimensionless groups (e.g.
  `η_s·V_c/τ_str`, a CFL-like `c_s·dt/h`), not absolute Pa/m/s. The fix and its
  proof are stated in those groups so they transfer to any mesh/material/config.
- **Exact local solve, covariant stabilization.** The local imposed-state solve
  is already the exact 2-vector Riemann solution (slip ∥ total tangential
  traction, `|τ_corr,t| = τ_str` when slipping). The instability is in the
  COUPLED ADER↔bulk map, so any stabilization must be applied **identically to
  both tangential components** (covariant), never to one channel.

## Constraints
- **Byte-exact TPV/BP5 regression (pure strike-slip) is non-negotiable.** The
  fix must be a strict no-op when one tangential channel is identically zero
  (the TPV/BP5 regime). Gate: `worst_rel ≤ 1e-13` vs the pre-change binary.
- **Files Requiring Extreme Care** (`CLAUDE.md`): `fault_face_flux.cpp/.hpp`,
  `dynamic/tpv205_friction.hpp`, `wave_operator.inl`. Full verification required.
- **Convention fixed:** `tangent1=dip`, `tangent2=strike`; canonical frame =
  ref-normal-aligned. No new sign conventions without auditing all sites
  (memory: feedback_complete_sign_sites). No hardcoded constants (memory:
  feedback_no_hardcoded_numbers).
- The fix MUST hold on **shared (cross-rank)** AND **interior (within-rank)**
  fault faces (both run the same `SolveLSW_TPV205` decomposition), and in
  **serial** (the instability has a serial manifestation with no cross-rank
  check to catch it).
- Local MPI at **np=2** (existing convention); strip the dup-`LC_RPATH` per
  binary on macOS (scripted this session).

## Background math (current scheme + where covariance breaks)

Per LSW sub-step (`tpv205_friction.hpp::SolveLSW_TPV205`,
`fault_face_flux.cpp::EvaluateADER_LSW`), all in the canonical frame:
```
τk_total = τk_0 + τk_nuc + τk_trial            (k=1 dip, 2 strike)
|τ_t|    = sqrt(τ1_total² + τ2_total²)
τ_str    = μ_eff(δ)·|σ_n_total|,   δ = |(slip1,slip2)|
V_abs    = max(0, (|τ_t| − τ_str)/η_s)
Vk       = V_abs · τk_total / |τ_t|             (slip ∥ total tangential traction)
τk_corr  = τk_trial − η_s·Vk
```
The full corrected tangential traction `τ_corr,t = τ_t − η_s·V_t = τ_t·τ_str/|τ_t|`
is exact and covariant (∥ τ_t, magnitude τ_str). **The instability is in the
discrete coupling:** the imposed state radiates `τ_corr,t` into the bulk; the
ADER predictor returns a new `τ_trial`; round-trip. For the **weak** component
(|τk| ≪ |τ_t|), the map `τk_trial(n) → τk_trial(n+1)` has gain >1 in the
strike-/dip-dominated regime (measured ~1.07 for TPV104, `fault_face_flux.cpp:
211-218`). Cross-rank, the self/nbr swap in the rotation
(`wave_operator.inl:3230-3239`) seeds a 1-ULP difference that the weak-channel
map amplifies to O(1). **A covariant scheme has the SAME gain for both
components at every rake → if pure-strike (TPV) is stable, every rake is stable.**

## Phase 1: GENERAL rake-sweep reproducer + cause/covariance discrimination (DIAGNOSTIC)

### Goal
A standalone np=2 + serial test that, **sweeping the rake angle φ ∈ {0,15,30,45,
60,75,90}° and ≥2 dimensionless-parameter scales**, reproduces the weak-channel
divergence, shows it is a function of rake (not of which physical channel), and
quantifies the per-sub-step loop gain `g(φ)` — pinning the cause and proving
the current scheme is NOT covariant.

### Files to Create
- `tests/unit/test_rupture_rake_sweep_cross_rank.cpp` — the reproducer (below).

### Files to Modify
- `Makefile` — 5 entries mirroring the tilted-test wiring (Makefile:323, 1669,
  3428): `TEST_RUPTURE_RAKE_SWEEP_SRC/_OBJ`, link target
  `seas_test_rupture_rake_sweep_cross_rank` (same object set as the tilted
  test), compile rule, run target `test-rupture-rake-sweep-cross-rank:` invoking
  `$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2`. NOT in the `make test` aggregate.
- `dynamic/fault_face_flux.cpp` — TEMPORARY, behind a new `#ifdef
  SEAS_DIAG_TANGENTIAL_LOOP`: after `SolveLSW_TPV205` returns
  (`:787`/`:904`), `fprintf` `(τ1_total, τ2_total, V_abs, V1, V2, τ1_corr,
  τ2_corr)` for the QP with `data.diag_print`. Removed at end of Phase 1.

### Detailed Requirements
1. **Rake parameterised loading.** Reuse `BuildTiltedTwoTetFaultMesh` (normal
   stays in xy-plane; tangent frame well-defined). For each rake φ, apply a
   tangential overstress of fixed magnitude `Δτ = TPV102Params::nuc_dtau`
   distributed as `tau1_nuc = Δτ·sin φ`, `tau2_nuc = Δτ·cos φ` (φ=0 → pure
   strike = TPV regime; φ=90 → pure dip; in between → oblique). Optionally also
   sweep an oblique *background* (project a rotated σ so `τ_0` carries the rake)
   to exercise `τk_0≠0`.
2. **LSW path, per-DOF params derived (no magic numbers).** Set
   `lsw_mu_s/mu_d/d_c` from `TPV205Params`/resolved ratios so the patch ruptures
   at every φ (assert `max V_abs > 0` per φ, else FAIL "did not nucleate").
3. **Two parameter scales.** Repeat the sweep at ≥2 values of the controlling
   dimensionless group (e.g. scale `η_s` and `dt` to change `c_s·dt/h` and
   `η_s·V_c/τ_str` by ~2×), to show the instability is regime- not config-bound.
4. **Run past the divergence horizon:** `kNSteps ≥ 600` (seed→O(1) at gain 1.07
   is ~530 steps).
5. **Cross-rank check each step, non-aborting.** Call
   `wave_p.VerifySharedFaultDOFDataConsistency(t, &worst_rel, &worst_field,
   /*abort_on_fail=*/false)` (new overload, Interfaces) and record
   `worst_rel`/`worst_field` per (φ, scale, step).
6. **Serial boundedness leg.** Also run serial at each φ and assert
   `|V1|,|V2|` stay bounded by a covariant cap derived from the friction
   solution (`V_abs ≤ (|τ_t|−τ_str)/η_s` with `|τ_t|` from the imposed
   loading) — catches a physically-wrong-but-cross-rank-consistent fix.
7. **Covariance/loop-gain output (the headline):**
   - Estimate `g(φ)` = geometric growth rate of `worst_rel` (or of `|V_weak|`
     in serial) over the linear-growth window, per φ and scale.
   - Print a verdict: (i) if `τk_0/τk_nuc/can_tk` differ across ranks at step 0
     → consistency gap (H-A/H-C), Phase 2; (ii) if identical at step 0 and
     `g(φ) > 1` for the weak-channel rakes but `g(0)=g(90)` are NOT equal
     → **non-covariant instability (H-B)**, Phase 3; (iii) `g(φ)` symmetric
     about φ=45° and ≤1 everywhere → already covariant/stable (unexpected;
     re-examine seed).

### Interfaces
- Add to `wave_operator.hpp/.inl`:
  `void VerifySharedFaultDOFDataConsistency(double t,
       double *worst_rel_out = nullptr, int *worst_field_out = nullptr,
       bool abort_on_fail = true) const;`
  The existing driver call (`abort_on_fail` defaulted true, out-ptrs null) is
  unchanged. Factor the field-comparison body so both modes share it.
- Reuse `TPV102Params`, `TPV205Params`, `BoundaryConfig`, `FaultFaceFlux`,
  `WaveOperator<Mesh>`/`<ParMesh>` exactly as the tilted test does.

### Edge Cases to Handle
- np≠2 → SKIPPED, return 77 (mirror tilted test).
- φ=0 (pure strike) MUST behave exactly as the existing stable TPV path
  (sanity anchor); φ=90 (pure dip) is its mirror.
- QP that never ruptures → assert nucleation occurred (Req 2).
- Degenerate normal ∥ up → excluded by the tilted-mesh geometry.

### Acceptance Criteria
- [ ] `make test-rupture-rake-sweep-cross-rank` builds + runs at np=2.
- [ ] Reproduces failure: ∃ (φ, step) with `worst_rel(V_weak) > 1e-10` (red).
- [ ] Output proves NON-covariance: `g(φ)` is NOT symmetric about 45° / exceeds
      1 for the weak-channel rakes while φ=0 is stable. (Or classifies H-A/H-C.)
- [ ] No functional source changed (only the non-aborting verify overload + a
      guarded diagnostic).

### Dependencies
- Depends on: nothing. Required by: Phase 2/3 (this is the oracle).

## Phase 2: Cross-rank consistency fix (CONDITIONAL — only if Phase 1 shows H-A/H-C)

### Goal
`τk_0`, `τk_nuc`, and the canonical frame are bit-identical across ranks at
step 0 for ALL rakes.

### Files to Modify (chosen from Phase-1 output; expected no-op given invariance)
- `drivers/spatial_dyn_driver.cpp::BuildPerDOFFaultTables` (:270-303): emit the
  **canonical** frame (`sign_flipped ? -t : t`) into `dof_basis` instead of the
  stored frame, so prestress projection and runtime flux use one identical
  frame. Only if Phase 1 shows a real step-0 inconsistency.
- `dynamic/spatial_nucleation.cpp:174`: project the nucleation amplitude onto
  the canonical frame if shown inconsistent.

### Acceptance Criteria
- [ ] Phase-1 test: step-0 per-rank `τk_0/τk_nuc/can_tk` rel_diff ≤ 16 ULP, all φ.
- [ ] TPV/BP5 byte-exact unchanged.

### Dependencies
- Depends on: Phase 1. Required by: Phase 3.

## Phase 3: GENERAL covariant tangential-traction stabilization (expected main work, H-B)

### Goal
The fault tangential-traction update is rotationally covariant: the
per-sub-step map's spectral radius is ≤1 and **independent of rake φ** across
the documented dimensionless-parameter range; `V1/V2/τ*_corr` stay bounded and
cross-rank-consistent for all φ; pure strike-slip (TPV/BP5) is byte-exact.

### Files to Modify
- `dynamic/tpv205_friction.hpp::SolveLSW_TPV205` (:93-167) and its rate-state
  analog `dynamic/fault_face_flux.cpp::CompleteFromVabs` (:190-242) — keep them
  parallel.
- Possibly the ADER imposed-state coupling in `wave_operator.inl` (the
  `Q_imp`→bulk path) if Phase 1 localizes the gain to the coupling, not the
  local solve.

### Candidate fixes (design + sign-off before implementing; chosen by Phase-1
gain localization). All MUST be covariant (applied identically to τ1 and τ2):
1. **Covariant implicit (vector) imposed-state update.** Solve the 2-vector
   fixed point `τ_corr,t = τ_trial,t − η_s·V_t`, `V_t = V_abs·τ_corr,t/|τ_corr,t|`,
   `|τ_corr,t| = τ_str` directly for the vector (one Newton/closed-form on the
   2-vector), instead of scalar `V_abs` + post-hoc split along `τ_total`. This
   is the exact local Riemann BC and, applied in the ADER corrector implicitly,
   removes the weak-component gain. (Closest to
   `PLAN_phase_R_exact_bimaterial_riemann.md`.)
2. **Covariant damping of the radiated tangential traction.** If Phase 1 shows
   the gain is in the ADER↔bulk coupling, apply the SAME stabilization
   (e.g. consistent treatment of `τ_corr,t` in the predictor and corrector) to
   the full 2-vector `τ_corr,t`, never to one component.
3. **Rejected:** any per-channel / rake-thresholded / env-gated suppression
   (violates covariance + the no-workaround constraint).

### Detailed Requirements
1. **Provably reduces to current code when one channel ≡ 0** (pure strike or
   pure dip) — bit-exact. State the algebra in the code comment; verify by
   regression. (This is what makes φ=0 byte-exact AND φ=90 equally stable.)
2. **Covariance proof:** the per-sub-step Jacobian of `τ_trial,t(n)→τ_trial,t(n+1)`
   is `≤1` in spectral radius and its eigenvalues are **independent of φ**
   (rotationally invariant). Phase 1 measures `g(φ)` before (non-flat, >1 at
   weak rakes); after the fix it must be flat in φ and ≤1.
3. **Scale-invariance:** stated/verified in dimensionless groups; the ≥2-scale
   sweep (Phase 1 Req 3) must pass unchanged by the fix.
4. No env-vars, no hardcoded constants, no rake branches.

### Edge Cases to Handle
- `|τ_t| → 0` (unloaded): `V=0`, `τ_corr=τ_trial`, no division by zero (the
  vector solve must guard `|τ_corr,t|>0`; preserve the existing `tau_abs>0`).
- Pure dip (`τ2=0`) and pure strike (`τ1=0`) are mirror images — both stable.
- Bimaterial (`Zs_plus≠Zs_minus`) out of scope (guarded elsewhere); homogeneous
  path must not regress.

### Acceptance Criteria
- [ ] Phase-1 rake sweep: `worst_rel ≤ 1e-10` for ALL φ and ALL steps and BOTH
      parameter scales; serial `|V|` bounded for all φ.
- [ ] Measured `g(φ)` is flat in φ (max/min over φ within 1% of 1.0) and ≤1.
- [ ] TPV102/104/205 + BP5 byte-exact: `worst_rel ≤ 1e-13` vs pre-Phase-3 binary.
- [ ] LSW/fault tests green: `test-ader-tpv102-smoke`,
      `test-tpv102-total-locked-fault`, `test-shared-fault-dof-data-consistency`,
      `test-rupture-multistep-serial-vs-parallel`,
      `test-rupture-tilted-fault-serial-vs-parallel`.

### Dependencies
- Depends on: Phase 1 (cause + gain), Phase 2 (consistency precondition).
  Required by: Phase 4.

## Phase 4: Full regression + Frontera re-validation

### Goal
End-to-end: SAFS Dc2 advances past t=1.0 s without R-101 abort or blow-up; the
fix is config-agnostic (also sanity-run a non-Dc2 config if quick).

### Detailed Requirements
1. Local regression battery (fault-basis trio; shared-fault role/dof-data
   consistency; interior-flux-path ×3; godunov identity; tpv102
   locked/absorbing/pepper/ader-smoke; multistep serial-vs-parallel; tilted
   serial-vs-parallel; sign-flipped truth table) — green; document unchanged
   pre-existing failures (adjacent-triangle pepper-bug, r101 missing-precondition,
   macOS MUMPS Bus error).
2. Rebuild `seas_spatial_dyn_driver`; commit; push `safs`; user re-runs the Dc2
   sbatch on Frontera. **Also** re-run (or have the user re-run) the
   physical-`D_c=1.0` baseline config to confirm the fix is not tuned to Dc2.
3. Read logs: R-101 passes through nucleation; `V_max` peaks then decreases; no
   tensile `σ_n` runaway; reaches `tfinal` (or well past t=1.0 s).

### Acceptance Criteria
- [ ] Local regression green (modulo documented pre-existing).
- [ ] Frontera Dc2 AND D_c=1.0: no R-101 abort, no blow-up.

### Dependencies
- Depends on: Phase 3 (and Phase 2 if taken).

## Testing Strategy
- **Primary general oracle:** `test_rupture_rake_sweep_cross_rank` (Phase 1) —
  red before, green after, across the full rake sweep × ≥2 parameter scales,
  np=2 AND serial. Generality is *demonstrated* by φ-independence of stability,
  not by fixing the one SAFS QP.
- **Covariance metric:** the per-φ loop gain `g(φ)` must go from non-flat/>1
  (before) to flat/≤1 (after).
- **Byte-exact regression (pure strike-slip):** TPV102/104/205 + BP5 via the
  Phase-3/4 targets, captured pre-change and diffed to ≤1e-13; use the
  stash-rebuild-baseline method to prove the fix is a strict no-op on strike-slip
  (and, by covariance, on pure dip).
- **Consistency guard:** the in-code `VerifySharedFaultDOFDataConsistency`
  (non-aborting test mode) — the same check that fails on Frontera.

## Risk Assessment
- **A fix that stabilises SAFS's rake but not all rakes** (the exact failure
  mode "general not problem-specific" guards against). Mitigation: the rake
  sweep + the φ-independence acceptance gate (Phase 3 AC 2) make a
  rake-specific fix fail the test.
- **Breaking strike-slip byte-exactness** (TPV/BP5 contract — highest impact).
  Mitigation: the algebraic reduction-to-current-code-when-one-channel≡0
  requirement (Phase 3 Req 1) + byte-exact gate.
- **Hidden non-covariance elsewhere** (e.g. the canonical-frame self/nbr
  rotation, or strike=up×n vs dip=strike×n giving the two channels different FP
  conditioning). Phase 1's φ-sweep + serial leg surface it; if found, the fix
  must remove the asymmetry, not compensate for it.
- **Gain is in the coupling, not the local solve.** Phase 1's `τk_trial` trace
  isolates intra-sub-step (local) vs across-sub-step (coupling) growth; Candidate
  2 addresses the coupling case.
- **Tricky existing code:** `wave_operator.inl:3209-3275` (canonical frame +
  self/nbr rotation = the FP-seed origin), `SolveLSW_TPV205` (the
  decomposition), `fault_face_flux.cpp:211-230` (documented instability + the
  env lock that must stay a no-op).
