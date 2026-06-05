# Code Review: Phase 2 (resample operator R) + Phase 3 (rate-state Δψ resample) + Phase 4 wiring — 2026-06-04

## Review Scope
- **Plan:** `document/fault_dealiasing_dev/seisol_overintegration_resample_speckle_2026-06-02.{md,pdf}`, §4.2, §6 (Phase 2/3/4), §Acceptance, §7 (validation).
- **Files reviewed** (Phase 1–3 worktree changes on `feat/fault-overint-resample`, commits `1cc0bbf`→`ad7f9d6`):
  - `dynamic/fault_resample.hpp` (new — builds R, applies R)
  - `dynamic/friction_substep_iterator.{hpp,cpp}` (Phase 3 Δψ resample in `RateStateSubStepIterator::Advance`)
  - `dynamic/friction_iterator.hpp` + `dynamic/friction_iterator_factory.{hpp,cpp}` (`SetFaultResample` interface + factory)
  - `drivers/spatial_dyn_driver.cpp` (Phase 1 `--fault-overint` wiring, Phase 3 `--fault-resample` build/gate)
  - `dynamic/wave_operator.{hpp,inl}` (`SetFaultOverint` / `FaultFaceQuadDegree` / `RebuildFaultQuadrature_`, fault-aware quad sites)
  - `tests/unit/test_fault_resample.cpp`, `tests/unit/test_fault_resample_apply.cpp`
- **Cross-referenced:** `dynamic/fault_face_flux.{hpp,cpp}` (`WriteBackState`/`BuildImposedState` — confirms `data.psi` is NOT touched by write-back), `VerifySharedFaultDOFDataConsistency`/`ExchangeAndPairSharedFaultQPs` (wave_operator.inl:6079,6300 — shared-fault `(face_key, qp_idx)` pairing).
- **Prior reviews consumed:** `REVIEW_phase0_code.md`, `REVIEW_phase1_code.md`, the plan `REVIEW.md`. Findings there (R-001…R-005 Phase-0/1) are NOT re-litigated.
- **Domain context:** `CLAUDE.md`, `miniapps/seas/CLAUDE.md` (sign conventions, total-Q dispatch, shared-fault redundancy, linear gmsh tet meshes).
- **Not run:** no local reproducer / TPV runs (project policy `[[feedback_no_local_reproducer]]`); this is a static audit.

## What is correct (verified, not assumed)
These were checked specifically because they are the easiest places to get the resample wrong; all hold:
- **Linear algebra of R** (`fault_resample.hpp:89-98`): `MultAtB(V,WV,M)`=VᵀWV, `Mult(V,M,A)`=V·M⁻¹, `MultABt(A,WV,R)`=V M⁻¹ (WV)ᵀ = V(VᵀWV)⁻¹VᵀW. Idempotent, rank (N+1)(N+2)/2, W-self-adjoint — matches the test.
- **Δψ apply** (`friction_substep_iterator.cpp:162-186`): snapshot `psi_before` → run sub-steps unchanged → `psi = psi_before + R·(psi_after−psi_before)` per face. Increment (not full state) is projected, matching plan §4.2 / SeisSol RS. τ_corr / I_imp are built from the **un-resampled** rate (plan §6 Phase 3 "do not rebuild τ_corr from a resampled quantity") — confirmed: resample runs *after* `RunSubSteps_` fills `I_imp`, and `WriteBackState` does not touch `psi` (`fault_face_flux.hpp:436`).
- **Per-face block layout** matches between R and `dof_data`: `BuildPerDOFFaultTables` (driver:327-362) lays out interior-faces-then-shared, each a contiguous `nbf` block in `IntRules.Get(geom, FaultFaceQuadDegree()).IntPoint(q)` order; R is built from the *same* rule (driver:1446). `ndof % nbf == 0` guard present.
- **Parallel consistency** of the per-face resample: shared-fault QPs are paired by `(face_key, qp_idx)` with the claim (wave_operator.inl:6354-6357) that `qp_idx=k` is the same physical point on both ranks. Given that invariant (which the *existing* consistency tripwire already depends on), both ranks apply R to identically-ordered Δψ ⇒ resampled `psi` stays bit-consistent. No new MPI hazard.
- **Both-knobs-off byte-exactness:** `FaultFaceQuadDegree()`=`2*order_` at k=0; resample disabled ⇒ byte-exact.
- **Factory wiring:** `MakeFrictionIterator` returns the **unified** iterators (`friction_iterator_factory.cpp:32,93,101`), all of which inherit `SubStepIteratorBase::SetFaultResample`, so `--fault-resample` is *not* silently dropped on the production path (the deprecated `*FrictionIterator` adapters that no-op `SetFaultResample` are not used).

---

## Findings

### [R-001] MODERATE [spatial_dyn_driver.cpp:2719 / friction_iterator_factory.cpp:32] — LSW resample (the plan's decisive TPV31 path) is not implemented; Phase 3 is PARTIAL

> **RESOLVED 2026-06-04** (user-approved "implement it next"). See the **Fix log** at the bottom of this file. The LSW slip-magnitude resample is now implemented in `LinearSlipWeakeningIterator::Advance` (path-length increment `Σ|V|·dt` → degree-N L2 projection per face → rescale directional slip to the dealiased magnitude); the driver abort is removed. Verified: new tests T4/T5 pass (`slip2_on = R·slip2_off` exactly, inert at R=I), and the bit-parity oracle test is **36/36** (resample-off byte-exactness preserved).

**Category:** DEVIATION / INCOMPLETE

**Description:**
Plan §6 Phase 3 specifies **two** resample targets: "**LSW** → resample the slip-rate magnitude and integrate it into the accumulated slip; **rate-state** → resample the state-variable increment Δψ." Only the rate-state branch exists. The LSW branch is a hard abort:

```cpp
// spatial_dyn_driver.cpp:2719
MFEM_VERIFY(!(cli_fault_resample && is_lsw),
            "--fault-resample is not implemented for LSW (TPV31/TPV205) ...");
```

`LinearSlipWeakeningIterator::Advance` (friction_substep_iterator.cpp:255-335) contains no resample and ignores the stored `resample_R_`.

**Why this matters:** TPV31 — the plan's explicitly **decisive** target (§7.3: "the decisive target … predict constant σ_n, dip slip → 0, completion to 15 s") — is **LSW**. Per the plan's own causal model (§3.4), the TPV31 blow-up is the *rectified secular* drift, i.e. the channel that **resample** (not over-integration) is supposed to close (§3.3 channel 2, §4.2). So the one technique the plan singles out as the secular cure for TPV31 cannot be exercised on TPV31. The plan acceptance §6 and validation §7.3 cannot be met with the current code.

**Trigger:** `seas_spatial_dyn_driver … (TPV31 LSW config) --fault-overint 2 --fault-resample` → aborts at step setup.

**Actual behavior:** loud abort (good — not silent), but Phase 3 LSW is absent.

**Expected behavior:** implement the LSW slip-magnitude resample (gather `|V|·dt`/Δslip per face, `R·Δslip`, scatter-add into accumulated slip), OR the plan/STATUS must record Phase 3 as rate-state-only and de-scope the §7.3 TPV31-with-resample acceptance until it lands.

**Suggested fix:** (scope decision required — recommend implementing, since it is the plan's headline target) add an LSW resample mirroring the RS branch, acting on the accumulated-slip increment:
```diff
// friction_substep_iterator.cpp — LinearSlipWeakeningIterator::Advance, after RunSubSteps_
+  // Phase 3 (LSW): resample the per-macro-step accumulated-slip-magnitude increment.
+  if (do_resample) {            // same gate/contiguity guards as the RS branch
+     // gather Δ|slip| = sqrt(slip1²+slip2²)_after − _before per face,
+     // proj = R·Δ|slip|, then rescale (slip1,slip2) to the projected magnitude.
+  }
```
and replace the abort at `spatial_dyn_driver.cpp:2719` with the LSW wiring. **Until then, keep the abort** — do not silently no-op.

**Test case:**
```cpp
// tests/unit/test_fault_resample_apply.cpp — add an LSW case once implemented
void test_R001_lsw_resample_changes_accumulated_slip() {
   // 1 over-integrated face (#QP > #DOF), distinct per-QP Δslip with >N content.
   // Assert: with resample ON, the per-QP accumulated slip == old_slip + R·Δ|slip|
   //         (direction preserved), and != the no-resample result (non-vacuous).
}
```

---

### [R-002] MODERATE [spatial_dyn_driver.cpp:2722-2726] — `--fault-resample` alone is NOT byte-exact at p ≥ 3; the gate keys off `nbf != ndof` instead of "over-integration on", contradicting the plan acceptance

**Category:** DEVIATION

**Description:**
Plan §Acceptance: "`--fault-resample` alone (no over-integration) ⇒ also byte-exact, because R = I (Phase 2 test iv)." Phase 2 test (iv) is "at the minimal rule (#GP = #DOF), R = I." That premise holds only on **tensor (quad) faces**; on a **triangle** the minimal `2*order` rule has #GP > #DOF for N ≥ 3 (e.g. N=3: `IntRules.Get(TRIANGLE,6)` = 12 QP > 10 DOF), so R ≠ I. The driver gates resample on exactly this geometric coincidence:

```cpp
// spatial_dyn_driver.cpp:2722-2724
const int ndof_per_face = (cfg.mesh.order + 1) * (cfg.mesh.order + 2) / 2;
const bool resample_active =
   cli_fault_resample && (nbf_per_face != ndof_per_face);
```

Consequence: at p ≥ 3, `--fault-resample` **alone** (no `--fault-overint`) ⇒ `nbf_per_face != ndof_per_face` ⇒ `resample_active == true` ⇒ the Δψ resample **acts** and changes results, while the flux is *not* over-integrated. This violates the stated regression contract and runs the resample in an off-design regime (channel-2 cleanup without the channel-1 prerequisite the plan §3.3/§4.2 says it must be "layered on top of"). It is acknowledged in the code comment ("at p>=3 … `--fault-resample` alone DOES act — review R-002") but it still breaks an explicit acceptance criterion and is untested.

**Trigger:** `--fault-resample` with no `--fault-overint`, mesh order ≥ 3 (p3 production configs exist — see `e093a67` "p2/p3 configs").

**Actual behavior:** p3 results differ from the no-knobs baseline even though only `--fault-resample` (which the plan promises is inert without over-integration) was passed.

**Expected behavior:** `--fault-resample` without `--fault-overint` is a true no-op at every order (plan: "resample is a no-op without over-integration; over-integration is the prerequisite").

**Suggested fix:** gate on the actual prerequisite — over-integration being on — not on the triangle-rule head-count:
```diff
-   const int ndof_per_face = (cfg.mesh.order + 1) * (cfg.mesh.order + 2) / 2;
-   const bool resample_active =
-      cli_fault_resample && (nbf_per_face != ndof_per_face);
+   // Plan §4.2: resample is the secular layer ON TOP OF over-integration; with
+   // over-integration off it must be a TRUE no-op at every order (incl. p>=3,
+   // where the minimal triangle rule is over-determined and R != I).
+   const bool resample_active = cli_fault_resample && (cli_fault_overint > 0);
```
(If acting at the over-determined minimal rule is genuinely wanted, then instead update the plan §Acceptance to carve out the p ≥ 3 triangle exception and add a regression test for it — but do not leave spec and code in conflict.)

**Test case:**
```cpp
void test_R002_resample_alone_is_noop_at_p3() {
   // order=3, minimal rule (no over-int). With the fix, resample_active must be
   // false, so --fault-resample alone reproduces the baseline psi bit-for-bit.
   const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2*3); // 12 QP
   DenseMatrix R; BuildFaultResampleMatrix(Geometry::TRIANGLE, 3, ir, R);
   // R != I here (12>10): prove the gate, not R, is what guarantees byte-exactness.
   assert(/* maxabs(R - I) */ > 1e-3);
   // and assert the driver gate (cli_fault_overint==0) => resample disabled.
}
```

---

### [R-003] LOW [POSSIBLE] [fault_resample.hpp:74-87] — R uses `ip.weight` only (no per-QP |J_F|); valid only for affine faces, with no guard

**Category:** ASSUMPTION

**Description:**
`BuildFaultResampleMatrix` builds `W = diag(ip.weight)` from the **reference** quadrature weights, omitting the face Jacobian `|J_F(x_q)|`. For an **affine** face `|J_F|` is constant and factors out of `R = V(VᵀWV)⁻¹VᵀW` (scale-invariant), so one reference R is exact for every straight-sided triangular face regardless of size — correct for all current targets (TPV*/SAFS use linear, straight-sided tets ⇒ flat faces, even on a faceted curved fault surface; the per-QP normal in `BuildPerDOFFaultTables` varies *between* faces, not within one). For a genuinely **curved / isoparametric** fault face, `|J_F|` varies within the face and does **not** cancel, so R is then the reference-measure projector, not the physical-L²(dA) projector the plan Eq. (4.2) specifies. There is no assert tying `--fault-resample` to affine faces.

**Trigger:** `--fault-resample` on a high-order/curved fault mesh (not used today, but `ComputeQPBasis`'s per-QP frame support invites it).

**Actual / Expected:** silently-wrong projector vs. an L²(dA) projector. Add a guard now so the latent trap fails loud if a curved mesh is ever used.

**Suggested fix:**
```diff
   const int nq   = ir.GetNPoints();
   MFEM_VERIFY(nq >= ndof, ...);
+  // Affine-face assumption: |J_F| cancels in the projector only when constant
+  // over the face. Callers that build ONE reference R for all faces must pass
+  // an affine (straight-sided) face geometry. (Document, or weight W by the
+  // per-QP |J_F| and build R per face for curved geometries.)
```
plus a driver-side `MFEM_VERIFY(fault face is affine)` when `cli_fault_resample` is set.

**Test case:** `test_R003`: build R for a deliberately curved reference map; assert it differs from the |J_F|-weighted projector, documenting the limitation (or, after a fix, assert equality of the per-face |J_F|-weighted R with the analytic L² projection of a degree-(N+1) field).

---

### [R-004] LOW [spatial_dyn_driver.cpp ~2719] — Phase-4 "disable resample on across-fault material contrast" guard is missing

**Category:** DEVIATION

**Description:**
Plan §6 Phase 4 step 2 (and §4.2 caveat) explicitly requires: "Add a guard that **disables resample on any genuine across-fault material contrast** (matching SeisSol's `BiMaterialFault` … 'resampling introduces artificial oscillations')." No such guard exists. Not currently reachable — the spatial driver's `matrix`/`BimaterialWaveOperator` path is for **depth**-heterogeneity (material identical across the fault at each depth), not a fault-normal jump — but the plan mandates the guard so a future genuine-contrast SAFS config does not silently enable an oscillation-inducing resample.

**Suggested fix:** when `cli_fault_resample`, assert (or warn-and-disable) on a detected across-fault `(ρ,cp,cs)` jump (the per-face `Zp_plus != Zp_minus` / `Zs_plus != Zs_minus` already in `DOFData`):
```diff
+  // Plan §4.2: resample must be OFF across a genuine fault-normal material
+  // contrast (SeisSol BiMaterialFault). Detect via per-QP impedance jump.
+  if (cli_fault_resample) {
+     bool contrast = /* any DOFData q with Zp_plus!=Zp_minus or Zs_plus!=Zs_minus */;
+     MFEM_VERIFY(!contrast, "--fault-resample is invalid across a fault-normal "
+                 "material contrast (artificial oscillations); disable it.");
+  }
```

**Test case:** `test_R004`: construct two `DOFData` with `Zp_plus != Zp_minus`; assert the driver guard rejects `--fault-resample`.

---

### [R-005] LOW [test_fault_resample_apply.cpp / friction_iterator_factory.hpp:26-27] — Phase-3 test gaps + stale factory doc comment

**Category:** QUALITY

**Description:**
1. `test_fault_resample_apply.cpp` exercises only `RateStateAgingIterator`. The **SRW** policy (`RateStateSlipLawSrwIterator` = TPV104, a named target) goes through the *same* templated `Advance`, but is untested for the resample. No test asserts the **conservation** property the design leans on (R is the W-orthogonal projector ⇒ it preserves the W-weighted mean of Δψ; a uniform Δψ must be unchanged). No parallel/shared-face resample-consistency test (the serial-vs-parallel invariant the plan Phase 0 calls out).
2. `friction_iterator_factory.hpp:26-27` documents the factory as returning the **deprecated** `LswFrictionIterator` / `RateStateAgingFrictionIterator`; the `.cpp` actually returns the unified `LinearSlipWeakeningIterator` / `RateStateAgingIterator` / `RateStateSlipLawSrwIterator`. The deprecated adapters *no-op* `SetFaultResample`, so the stale comment points a reader at the one wiring that would silently break the resample — fix the comment to prevent a future "why is resample dead?" regression.

**Suggested fix:** (a) add an SRW case and a "uniform Δψ is reproduced exactly / W-weighted mean preserved" assertion to `test_fault_resample_apply.cpp`; (b) update the factory header comment to name the unified iterators it returns.

**Test case:**
```cpp
void test_R005_uniform_dpsi_is_preserved() {
   // Over-integrated face; set every QP's psi increment to the same constant c.
   // R reproduces constants exactly => proj == c at every QP (mean preserved).
}
```

---

## Summary
- Critical issues: **0** (no demonstrable wrong-results bug found in the implemented rate-state + planar-fault + opt-in path; the linear algebra, increment apply, per-face layout, MPI pairing, and both-knobs-off byte-exactness were each verified).
- Moderate issues: **2 — both RESOLVED** (R-001 LSW resample implemented; R-002 gate now keys on `--fault-overint > 0`). See Fix log.
- Low issues: **3 — all RESOLVED** (R-003 affine-face guard added; R-004 material-contrast guard added; R-005 conservation + SRW tests added, factory comment fixed).
- **Plan compliance: FULL** (for the implemented opt-in dealiasing pair). Phase 1 (over-integration), Phase 2 (R operator), Phase 3 (rate-state Δψ **and** LSW slip magnitude), and Phase 4 (driver wiring + material-contrast guard) are complete. The §Acceptance "resample-alone ⇒ byte-exact at every order" clause now holds (R-002 fix). Empirical §7 validation (Frontera TPV102/31 A/B) remains as runtime work, not a code defect.
- **Verdict: PASS.** All five findings fixed and verified locally (tiny fixtures only, per `[[feedback_no_local_reproducer]]`): `seas_test_fault_resample` 56/56, `seas_test_fault_resample_apply` 14/14, `seas_test_friction_substep_iterator_parity` 36/36 (byte-exact oracle intact), driver TU + full `seas_spatial_dyn_driver` link clean.

---

## Fix log

### R-001 — LSW slip-magnitude resample (implemented 2026-06-04, user-approved)

**Design (byte-exact-preserving, plan-faithful).** This codebase drives the LSW μ(δ) from `δ = √(slip1²+slip2²)` (displacement magnitude), whereas SeisSol tracks a separate path-length accumulator. To honor the byte-exact regression contract (μ-input unchanged when the knob is off) *and* the plan §4.2/§6 ("resample the slip-rate magnitude → accumulated slip that drives μ"), the fix:
1. accumulates the per-macro-step path-length increment `Σ_substeps |V|·dt` per QP (gated on `do_resample` ⇒ zero overhead / byte-exact when off);
2. projects it per face: `proj = R · Δ|slip|`;
3. sets the dealiased magnitude `δ_new = max(0, δ_before + proj)` and **rescales** `(slip1,slip2)` to it, preserving direction (un-resampled, from the friction solve). τ_corr and the slip rate stay from the un-resampled solve (computed in `StepOneQP_` before the resample), per plan §6 Phase 3.
4. **0/0 guard:** a QP with `δ_after == 0` (front tip, no slip direction yet) is left at 0 — the small dealiased increment is dropped rather than assigning an ill-defined direction.

**Files changed:**
- `dynamic/friction_substep_iterator.cpp` — `LinearSlipWeakeningIterator::Advance`: snapshot + path-length accumulation + post-loop per-face projection/rescale.
- `drivers/spatial_dyn_driver.cpp` — removed the `MFEM_VERIFY(!(cli_fault_resample && is_lsw))` abort; updated the wiring comment. (Gate `resample_active` unchanged — R-002 deliberately **not** bundled.)
- `tests/unit/test_fault_resample_apply.cpp` — new `MakeLswDOF`/`RunLSW` (pure strike-slip, zero Q) + **T4** (R=I inert) and **T5** (`slip2_on = R·slip2_off` per face to <1e-9·‖·‖, and non-vacuous).

**Verification (local, tiny fixtures — `[[feedback_no_local_reproducer]]` respected):**
- `seas_test_fault_resample_apply`: **12/12** (incl. 4 new LSW).
- `seas_test_friction_substep_iterator_parity`: **36/36** — LSW DOFData bit-parity intact ⇒ resample-off path byte-exact.
- `drivers/spatial_dyn_driver.o` compiles clean (abort removal).

**Caveat / follow-up:** `δ_before` is the displacement magnitude while `Δ|slip|` is path length; these coincide for monotonic single-direction LSW slip (the TPV/SAFS regime) and the mixing is bounded/secular. Not byte-exact against a hypothetical separate path-length accumulator, but byte-exact against today when off — which is the binding contract. Frontera TPV31 A/B (plan §7.3) is the remaining empirical validation (not run here).

### R-002 — `--fault-resample` alone now byte-exact at every order (fixed 2026-06-04)
`spatial_dyn_driver.cpp`: `resample_active` now keys on `cli_fault_resample && (cli_fault_overint > 0)` (was `nbf_per_face != ndof_per_face`). With over-integration off the resample is a true no-op at **every** order, including p ≥ 3 where the minimal triangle rule is over-determined (R ≠ I). Removed the now-unused `ndof_per_face`; the R-build log note was retargeted to `cli_fault_overint`.

### R-003 — affine-face guard (fixed 2026-06-04)
`spatial_dyn_driver.cpp` (R-build block): when `--fault-resample` is set, samples `|J_F|` (`ftr->Face->Weight()`) at every rule QP on the probe fault face and `MFEM_VERIFY`s it is constant to 1e-10 relative. The single reference R is only the physical-L²(dA) projector for affine (straight-sided) faces; a curved face now fails loud instead of silently using the reference-measure projector. Passes trivially for the straight-sided tet meshes all current targets use.

### R-004 — material-contrast guard (fixed 2026-06-04)
`spatial_dyn_driver.cpp` (gate block): when `resample_active`, scans `dof_data` for any QP with `Zp_plus != Zp_minus` or `Zs_plus != Zs_minus` (MPI-reduced) and `MFEM_VERIFY`s none exist — matching SeisSol's `BiMaterialFault` (plan §4.2/§Phase 4). The spatial setup always sets the two sides equal per DOF (`spatial_setup.hpp:85,131`), so it never trips today; it guards a future genuine-contrast config.

### R-005 — test gaps + stale factory comment (fixed 2026-06-04)
- `tests/unit/test_fault_resample.cpp`: added a **conservation** check per `RunCase` — R reproduces a uniform field exactly and preserves the W-weighted mean of any field (no net secular drift). Now **56/56**.
- `tests/unit/test_fault_resample_apply.cpp`: added the **SRW** (slip-law, TPV104) case — `psi = psi_before + R·Δψ` exactly + non-vacuous, confirming both rate-state policies are resampled. Now **14/14**.
- `dynamic/friction_iterator_factory.hpp`: doc comment now names the **unified** iterators the `.cpp` returns (was naming the deprecated `*FrictionIterator` adapters that no-op `SetFaultResample`).

## Unreviewed Areas
- **Phase-1 ADER flux-assembly indexing under over-integration** (`ComputeADERFaceFluxRHS` / `ComputeADERSharedFaceFluxRHS`, wave_operator.inl ~2649/3437/3964): the per-QP `I_imp` ↔ `dof_data` index mapping at grown `nbf` was spot-checked (fault-aware degree ternary present) but not exhaustively traced; covered by `REVIEW_phase1_code.md` + `test_fault_overint` (k=0 byte-exact). Re-trace if over-integrated production results look wrong.
- **Numerical efficacy** (does over-int + resample actually flatten σ_n / drive [[vₙ]]→round-off): a *validation* question, not a code-defect question; the Phase-0 harness + Frontera A/B (§7) own it. Out of scope for this static audit and gated by `[[feedback_no_local_reproducer]]`.
- `Makefile` test-target plumbing and the two `ab_overint*_200m_normal.sbatch` job scripts: not reviewed (job config, not core numerics).
