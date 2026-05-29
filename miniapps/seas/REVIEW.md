# Code Review: RK45 + LSW + mixed-flux IMPLEMENTATION (adversarial audit, 2026-05-29)

> Reviews the implemented code (not the plan). The earlier plan review is archived
> at `REVIEW_rk45_lsw_plan_2026-05-29.md`; the Phase-13 review at
> `REVIEW_phase13_round2_2026-05-29.md`.

## Review Scope
- Plan: `document/mixed_flux_dev/PLAN_rk45_lsw_mixed_flux_2026-05-29.md` (Phases 1–5).
- Files reviewed (the diff vs HEAD + new files):
  - `dynamic/fault_face_flux.{hpp,cpp}` — `EvaluateLSW`
  - `dynamic/wave_operator.inl` — Mult-path LSW dispatch (interior `:~2617`, shared `:~3260`)
  - `dynamic/rk_time_stepper.hpp` — `AdvanceRKCoupledLSW_Spatial` + symmetric RS guard
  - `drivers/spatial_dyn_driver.cpp` — guard relax, gated rate_state guard, 3-way dispatch, banner
  - `tests/unit/test_lsw_rk_mixed_flux.cpp` (L1–L7), `tests/unit/test_lsw_rk_shared_fault_mpi.cpp` (L8)
  - `tpv205/configs/tpv205_spatial_rk45_mixedflux.toml`, `jobs/tpv205/tpv205_p1_rk45_mixedflux.sbatch`, `Makefile`
- Domain context: `miniapps/seas/CLAUDE.md` (sign/frame, no-local-full-mesh), the BUILD doc, the implementation report, memory (`preexisting-worktree-test-failures-2026-05`).

**Pass-1 (plan compliance):** all 5 phases implemented as specified; the R-002 (`=` not `+=`) and R-003 (symmetric guard) plan-fixes are present in code and exercised. **Pass-2 (correctness):** the algorithm is correct — verified by L1 (bit-exact `EvaluateLSW` vs `EvaluateADER_LSW`), L5 (bit-exact bulk RK4), L6 (absolute-combine vs `+=` discrimination), L8 (shared==interior==analytic), and the byte-exact guard suite (`rk_time_stepper` 22/22, `fault_face_flux_ader_equivalence` 11/11, `bimaterial_wave_operator_parity` 22/22). **No CRITICAL correctness bug found.** The findings below are real **coverage gaps** (untested production paths) and minor quality items — the kind that hide a *future* bug, not a present wrong result.

---

## Findings

### [R-001] [MODERATE] [test_lsw_rk_shared_fault_mpi.cpp] — L8 only exercises the shared-fault LSW dispatch AT REST (Q=0); the cross-rank Q± exchange is never tested with a non-trivial jump

**Category:** EDGE_CASE (the parallel path R-004-of-the-plan exists to guard is only half-covered)

**Description:**
L8 sets `Q = 0.0` on both the serial and the parallel run (lines 234 and 279). With
`Q = 0`, the trial traction is identically 0 on every fault QP, so the LSW slip-rate
`V = (τ_nuc − μ_s·σ_n)/η_s` is the SAME closed-form value at every QP — interior and
shared alike — driven purely by the per-DOF prestress in `DOFData`. The test therefore
proves the `:3260` arm *runs* and yields correct *at-rest* physics, but it does NOT
exercise the cross-rank Q± exchange: there is no velocity/stress jump for the shared-face
ghost exchange to carry, so a bug in the shared-fault Q± plumbing feeding `EvaluateLSW`
(e.g. swapped +/− sides, a missing rotation, or a stale ghost) would still pass — both
sides see Q=0. This is exactly the failure mode the plan's R-004/risk-2 wanted covered.

**Trigger:** any production run is dynamic (non-zero Q across the seam); L8 never is.

**Actual behavior:** L8 passes whether or not the shared-fault Q± pairing is correct,
because Q=0 removes the only signal that pairing affects.

**Expected behavior:** drive a one-sided Q perturbation (a VX kick on the −y side only,
as the sibling `test_interior_vs_shared_branch_live.cpp::SetOneSidedQ` does) so the
shared QPs see a real jump, then assert parallel-shared `V/τ*_corr` == serial-interior to
round-off.

**Suggested fix (extend L8):** add a one-sided Q kick before each `wave.Mult`, applied to
the element on the −y side, identical in the serial and parallel builds:
```cpp
// after `Vector Q(wave.Height()); Q = 0.0;` in BOTH the serial and parallel blocks:
const auto &fes = wave.GetFESpace();
const int ndof = fes.GetNDofs();
for (int e = 0; e < mesh.GetNE(); e++) {           // `mesh` = the (par)mesh in scope
   real_t cy = 0; Array<int> ev; mesh.GetElementVertices(e, ev);
   for (int v = 0; v < ev.Size(); v++) cy += mesh.GetVertex(ev[v])[1];
   if (cy / ev.Size() >= 0.0) continue;            // −y side only
   Array<int> ed; fes.GetElementDofs(e, ed);
   for (int j = 0; j < ed.Size(); j++) Q(VX * ndof + ed[j]) += 1.0e-5;
}
```
This makes V at the fault depend on the cross-rank jump, so the existing
`parallel == serial` assertion now actually tests the `:3260` Q± pairing.

**Test case:** the extension above is the test — with the kick, flipping the +/− side
selection in the `:3260` dispatch (or the ghost pairing) must make
`V_parallel != V_serial`, which the assertion `|V_parallel − V_serial| ≤ 1e-9·…` catches.

---

### [R-002] [MODERATE] [POSSIBLE] [drivers/spatial_dyn_driver.cpp] — the end-to-end rk45+LSW driver path past mesh-load is unverified locally

**Category:** ASSUMPTION (coverage gap on the integration path)

**Description:**
Local verification of the driver stops at the `--dry-run --verify-dispatch` banner, which
prints BEFORE the wave-operator assembly and aborts at mesh-load (the 101 MB TPV205 mesh
is absent and is a "full-mesh" op the project rule forbids locally). The unit tests (L1–L8)
drive `EvaluateLSW` / `AdvanceRKCoupledLSW_Spatial` / `Mult` DIRECTLY — they never go
through `spatial_dyn_driver`'s setup between mesh-load and the time loop. So the
driver-integration wiring for `is_rk && is_lsw` is unexercised locally:
  - `MakeFrictionIterator(cfg, fault_flux, /*rs=*/nullptr)` is still constructed for the
    LSW-RK path (`:~2515`) and then unused — if that construction asserts for a
    `time_integrator=rk45` LSW config, the run dies before the time loop, and nothing
    local catches it;
  - the 3-way dispatch actually selecting `AdvanceRKCoupledLSW_Spatial` (vs the RS arm)
    is only covered by code-reading, not execution.

**Trigger:** the first real `seas_spatial_dyn_driver … --time-integrator rk45 --mixed-flux
adjacent` run (Frontera).

**Actual behavior:** unknown past mesh-load locally; relies on the Frontera sbatch.

**Expected behavior:** at least one execution that reaches the time loop on the LSW-RK
path before production reliance.

**Suggested fix:** the `jobs/tpv205/tpv205_p1_rk45_mixedflux.sbatch` preflight already runs
`--dry-run --verify-dispatch --mesh <200m>` on Frontera, which DOES exercise mesh-load +
iterator construction + dispatch resolution. Make that the documented gate (it is), and —
to shrink the gap — add a coarse idev smoke: mesh `tpv2053d_200m.geo` at a coarse
characteristic length on an idev node (conda `pythonenv`) and run the same dry-run, OR add
a 1–2-step `--tfinal 1e-4` run on that coarse mesh to confirm the time loop enters the LSW
stepper. Document in the sbatch header that "scheme resolves locally; full driver path is
idev/Frontera-gated per the no-local-full-mesh rule."

**Test case:** N/A as a unit test (the driver is not unit-testable without a mesh); the
gate is the idev/Frontera dry-run. Acceptance: the dry-run on a real (coarse) mesh exits 0
with the `[time-integrator] … [LSW coupled-RK on (Q, slip)]` banner line printed.

---

### [R-003] [LOW] [fault_face_flux.cpp:EvaluateLSW] — `s.Theta` is computed but never used

**Category:** QUALITY

**Description:**
`EvaluateLSW` computes `s.Theta = sqrt(tau1_total² + tau2_total²)` (Step 2), but nothing
downstream reads it: `SolveLSW_TPV205` recomputes `tau_abs` from `tau{1,2}_total`
internally. It is dead arithmetic (one `sqrt` per fault QP per RK stage). It mirrors the
same dead line in `EvaluateADER_LSW`, so this is a pre-existing pattern, not a new defect —
but it is dead in the new code too.

**Trigger:** every `EvaluateLSW` call.

**Actual/Expected:** harmless; a wasted `sqrt`. Drop for clarity (or leave for parity with
`EvaluateADER_LSW`).

**Suggested fix:**
```diff
-   s.Theta         = std::sqrt(s.tau1_total * s.tau1_total
-                              + s.tau2_total * s.tau2_total);
+   // (No s.Theta needed: SolveLSW_TPV205 recomputes |tau_total| internally.)
```
(No test — dead-code removal, behaviour-neutral.)

---

### [R-004] [LOW] [fault_face_flux.cpp:EvaluateLSW] — no `SEAS_DIAG_FAULT_FLUX` diagnostic block on the LSW Mult path

**Category:** QUALITY

**Description:**
The rate-and-state `Evaluate` carries a `#ifdef SEAS_DIAG_FAULT_FLUX` hypocenter-QP trace
(the C-1/C-1n checkpoints). `EvaluateLSW` (like `EvaluateADER_LSW`) has none, so a
`SEAS_DIAG_FAULT_FLUX` build debugging an LSW-RK rupture gets zero fault-flux diagnostics
on the Mult path — the C-1 bisection tooling is silently unavailable for LSW. Consistent
with `EvaluateADER_LSW` (which also omits it), so not a regression; flagged so the gap is
on record if an LSW-RK rupture ever needs the C-1 trace.

**Trigger:** a `SEAS_DIAG_FAULT_FLUX` build running LSW + RK.

**Actual/Expected:** no diag output vs the RS path's per-QP trace. Optional: add a guarded
diag block mirroring `Evaluate`'s if/when LSW-RK needs bisection debugging.

**Suggested fix:** defer unless needed; if added, mirror `Evaluate`'s
`#ifdef SEAS_DIAG_FAULT_FLUX … if (data.diag_print) { … }` block after `ComputeTrialTraction`.
(No test — diagnostics only.)

---

### [R-005] [LOW] [test_lsw_rk_mixed_flux.cpp:L6 / DoLSWRK4Reference] — the L6 reference re-implements the stepper, so a shared conceptual error would pass

**Category:** QUALITY (test robustness)

**Description:**
L6's `DoLSWRK4Reference` reproduces the same coupled-(Q,slip) RK4 logic as
`AdvanceRKCoupledLSW_Spatial` (stage-local slip staging + absolute combine). If BOTH shared
the same conceptual error (e.g. both staged slip from the wrong stage subset), the
`stepper == reference` assertion would pass falsely. The test mitigates this with (a) the
discrimination assertion (`slip_abs != slip_bug`, which a `+=` regression fails) and (b) L1
anchoring the per-QP physics — so the residual risk is low, but the "matches reference"
check alone is not an independent oracle.

**Trigger:** a refactor that changes the staging convention in both stepper and reference.

**Actual/Expected:** the reference is faithful-by-construction, not independent. Acceptable
given the discrimination + L1 anchor, but worth noting.

**Suggested fix:** strengthen the independent anchor — assert L6's per-QP `V1/V2` after the
final stage matches the analytic LSW closed form `(|τ_total| − μ(δ)·σ_n)/η_s` at the
post-step `(Q_new, slip_new)` for at least one sliding QP (an oracle that does NOT depend
on the stepper's internal structure). (No code-under-test change; test-only hardening.)

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 shared-fault-at-rest-only; R-002 driver-path-unverified-locally)
- Low issues: 3 (R-003 dead Theta; R-004 no LSW diag block; R-005 L6 reference not independent)
- Plan compliance: **FULL** — all 5 phases implemented as specified; the plan's own
  review-fixes (`=` combine, symmetric guard, L1/L6/L8 tests) are present and pass.
- Verdict: **PASS WITH FIXES** — the implementation is correct and well-tested (L1–L8 +
  byte-exact guards green; the 5 `make test` failures are pre-existing, in files the diff
  touches 0 lines of). Apply R-001 (give L8 a real cross-rank jump) and R-002 (an
  idev/Frontera driver smoke past mesh-load) to close the two production-path coverage gaps
  before relying on the feature in production; R-003/R-004/R-005 are fix-when-convenient.

## Unreviewed Areas
- **CFL stability of RK45 + central flux at p1 on the real TPV205 mesh** — empirical,
  Frontera-only (the sbatch reads the stable dt off the dev smoke); cannot be assessed
  statically.
- **SCEC TPV205 physics correctness end-to-end** — post-merge Frontera benchmark vs the
  overlays; out of scope for a static code review.
- **`LSW_ForcedRupture` Mult abort** — reviewed (defensive, unreachable from the current
  driver since `is_lsw → LSW`); correct as a fail-loud guard.
