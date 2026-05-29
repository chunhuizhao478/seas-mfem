# Code Review: Phase 12 — Absorbing boundaries (PML) for the SAFS dynamic-rupture run (2026-05-29, fresh adversarial pass)

> Supersedes the prior 2026-05-26 `[DIAG-SIGN]` review previously at this path (recoverable via git history).
> Scope is **Phase 12 (PML)** only — not the Phase 10/11/13/14 work bundled in the same commit `30bf3e3`
> (Phase 13 has its own review at `miniapps/seas/REVIEW.md`).

## Review Scope
- Plan: `document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md` §Phase 12 (12.0–12.4)
  and the long-form reference `document/pml_dev/PLAN_pml_safs_2026-05-27.md` (Parts A–C, Phases 1–4).
- Files reviewed (read in full):
  - `dynamic/pml_layer.hpp`, `dynamic/pml_layer.cpp` (half-face mask, cubic profile, `d_max`)
  - `dynamic/wave_operator.{hpp,inl}` — `SetPML`/`GetPML` (`hpp:208–209,747`), the Mult/predictor PML branch
    (`inl:733–736`), `ApplyPMLDamping` (`inl:5606–5680`), the ADER-corrector inline PML branch (`inl:5329–5393`)
  - `drivers/spatial_dyn_driver.cpp` — CLI flags (`:539–546,652–657`), `[numerics]` merge, global bbox reduction
    (`:1111–1123`), `cp` guard (`:1124–1150`), construction block (`:1604–1752`: thickness derivation, face-mask,
    fault-clearance guard, `SetPML`, banner)
  - `spatial/code/spatial_friction.{hpp,cpp}` — `NumericsSpec` PML fields (`hpp:135–150`) + `[numerics]` parser
    (`cpp:1062,1108–1112`)
  - `tests/unit/test_pml.cpp` (6 tests), `Makefile` (`seas_test_pml` target, in default `test:`)
  - `jobs/safs/_pml_common.sh`, `spatial_dyn_pml_{on,off}_safs.sbatch`, `spatial_dyn_bigbox_ref_safs.sbatch`
  - `spatial/code/scripts/pml_reflection_metrics.py` (Phase 12.3 metric tool — existence + CLI wiring confirmed,
    not line-audited; see Unreviewed Areas)
- Domain context: `CLAUDE.md` (no-hardcoded-numbers rule, "no local full-mesh runs", don't-revert-fixes rule),
  memory notes (Phase-13 bimaterial split; pre-existing concurrent worktree edits).

## Verdict up front
Phase 12 is **wired correctly end-to-end and contains no correctness bug that this pass could find.** Several plausible
defects were specifically hunted and **verified correct** (listed below — the fix agent must NOT "fix" them, per the
don't-revert rule). The remaining findings are one MODERATE **test-coverage gap on the production (ADER-corrector) path**
and two LOW items (a documented `R_eff` relabel that is a user foot-gun, and stale TOML comments).

## Verified CORRECT (hunted, confirmed not bugs — do NOT touch)
1. **ADER-corrector damping target = `dt·Q_bg` (not `Q_bg`, not `0`).** `wave_operator.inl:5375`
   `I_qp[c] -= dt * bulk_bg_[c];` then `:5386–5387` `rhs[...] -= w * d_c * shape(i) * I_qp[c];`. The Mult/predictor
   path damps toward `Q_bg` without `dt` (`:5662` `Q_qp[c] -= bulk_bg_[c];`, residual `:5674`). Because Mult returns
   `dQ/dt` (`:713`) while the corrector operates on the time-integrated state `I≈dt·Q`, the two paths correctly differ
   by exactly the `dt` factor. This is the subtle point the plan warned must not be "simplified to the Mult form"
   (`PLAN_pml_safs_2026-05-27.md:434`) — it is implemented correctly.
2. **Global PML bounding box uses MIN for `lo`, MAX for `hi`.** `spatial_dyn_driver.cpp:1117` `MPI_Allreduce(... MPI_MIN ...)`
   for `pml_box_lo`; `:1120` `MPI_Allreduce(... MPI_MAX ...)` for `pml_box_hi`. Correct on a partitioned ParMesh.
3. **Fault→PML clearance guard uses a globally-reduced fault bbox.** `:1671` `f_lo` via `MPI_MIN`, `:1673` `f_hi` via
   `MPI_MAX`, with `±INF` sentinels (`:1657`) so empty-fault ranks are neutral; `have_fault` (`:1677`) and the
   `MFEM_VERIFY` clearance abort (`:1694`) are global/collective. Multi-rank-safe.
4. **Corrector PML branch is reached at `ader_order=2`** (the SAFS production order): `AdvanceADER` only asserts
   `order∈{2,3,4}` (`:5288`); the `if (pml_layer_)` block (`:5329`) is not order-gated, and the driver substep helper
   `AdvanceADERWithSubStep_Spatial` calls `wave.AdvanceADER(Q, dt, ader_order, Q_new)` (`:469`).
5. **Free-surface safety.** Driver builds `FaceXLo|FaceXHi|FaceYLo|FaceYHi` and adds `FaceZLo`/`FaceZHi` only when
   `pml_damp_bottom`/`pml_damp_top` (`:1641–1644`); for the SAFS box `Z[-41608,0]`, `z_max=0` is left undamped.
   `ComputeDamping` gates each half-face on its bit (`pml_layer.cpp:96–129`); `TestPMLFreeSurfaceTopUndamped` locks it.
6. **Guards present.** `MFEM_VERIFY(GetAbsorbingBackground()!=nullptr)` before `SetPML` (`:1618–1620`);
   `MFEM_VERIFY(pml_cp>0)` requiring Constant material (`:1621–1624`); abort if thickness cannot be derived
   (`mesh.lc_far_m<=0` and no explicit thickness, `:1628–1637`) — no hardcoded "3000 m" in C++; `has_bulk_bg_`
   asserted inside `ApplyPMLDamping` (`:5612`). Reflection warning suppressed when `use_pml` (`:1143`).
7. **Backward-compatible ctor / disabled path.** `half_face_mask=-1` derives the symmetric mask from `dirs`
   (`pml_layer.cpp:47–57`) so the four pre-existing tests are unchanged; with `use_pml=false` no `PMLLayer` is built and
   the run is byte-identical to pre-PML. `test_pml.cpp` compiles (named `Vector`/`Mesh` locals) and is in `make test`.

## Findings

### [R-001] [MODERATE] [tests/unit/test_pml.cpp] — the ADER-corrector PML path (the production SAFS path) has no automated test

**Category:** EDGE_CASE / coverage (regression risk; the code itself is correct-by-inspection)

**Description:**
All six tests drive the operator only through `RK4Step` → `wave.Mult(...)` (`test_pml.cpp:23–33, 134, 142, 225, 298,
350, 485`), exercising the **predictor/Mult** PML branch (`wave_operator.inl:735`). **None** advances with
`AdvanceADER`, so the **ADER-corrector PML branch** (`wave_operator.inl:5329–5393`) — the *only* PML path that runs in
SAFS production (`ader_order=2`) — is never executed by the suite. That branch is a hand-inlined re-implementation that
duplicates the elementwise loop and carries the `dt·Q_bg` subtlety (verified correct today, but a one-character edit
there, e.g. dropping the `dt` at `:5375`, would keep the whole suite green). The plan flags this exact path as its top
risk: "the SAFS run uses ADER (order 2), so the **corrector** must be the active path … if PML appears inert, the ADER
path isn't calling the damping" (`PLAN_pml_safs_2026-05-27.md:428`). The plan's own Phase-12.1 stability test
(`TestPMLFreeSurfacePulseStable`) also uses RK4, so neither the suite nor the plan AC guards the corrector.

**Trigger:** Any future edit localized to `wave_operator.inl:5329–5393` (wrong target, missing weight/`dt`, sign flip,
or accidental order-gating) — undetected by `make test`.

**Actual behavior:** Corrector PML is correct but unguarded by any local test.

**Expected behavior:** At least one test advances with the ADER corrector and asserts PML behavior (energy decays; and,
with non-zero `Q_bg`, the interior relaxes to `Q_bg`, which is what catches a dropped-`dt` regression).

**Suggested fix:** Add a corrector-path test mirroring `TestPMLEnergyDecay` but stepping via the public `AdvanceADER`
(the same entry the driver substep helper uses), with a **non-zero** background so the `dt·Q_bg` target is exercised:
```cpp
void TestPMLEnergyDecayADER()
{
   std::cout << "Test 21b: TestPMLEnergyDecayADER\n";
   // identical mesh/material/PML setup to TestPMLEnergyDecay, but a non-zero
   // uniform pre-stress background so the corrector's dt*Q_bg target matters.
   real_t bg[NUM_STATE] = {0};
   bg[SXX] = -1.0e6;                    // uniform compressive pre-stress
   wave.SetAbsorbingBackground(bg);
   // seed Q = bg + a localized P-pulse fluctuation (so the steady state is bg, not 0)
   wave.SetPML(&pml);

   const int order = 2;                 // SAFS production order
   real_t cfl = 1.0 / (3.0*(2.0*order+1));
   real_t dt = wave.ComputeMaxDt(cfl);

   real_t E_prev = ComputeEnergy(wave, Q, lambda, mu, rho);
   bool monotonic = true;
   Vector Qn(Q.Size());
   for (int step = 0; step < 50; step++)
   {
      wave.AdvanceADER(Q, dt, order, Qn);   // <-- corrector path (inl:5329)
      Q = Qn;
      real_t E = ComputeEnergy(wave, Q, lambda, mu, rho);
      if (E > E_prev * 1.001) { monotonic = false; break; }
      E_prev = E;
   }
   TEST_ASSERT(monotonic, "PML energy decays on the ADER corrector path");
   // The interior must settle toward the BACKGROUND, not toward zero:
   // a dropped dt at inl:5375 makes the PML shell relax to the wrong stress.
   // (Assert each interior QP's Q -> bg within tolerance, or that the fluctuation energy ->0.)
}
```
Register it in `main()` next to the other PML tests. (No production code change — this only adds the missing guard.)

**Test case:** the function above is the test case; it exercises the currently-untested corrector branch and would fail
if the `dt·Q_bg` target (`inl:5375`) regresses.

---

### [R-002] [LOW] [dynamic/pml_layer.cpp:59–68 + driver banner] — `target_R` is silently the *effective* reflection `R_eff = R_0^{3/4}`

**Category:** ASSUMPTION / QUALITY (documented & plan-ratified — not a correctness defect)

**Description:**
`d_max_ = (3*cp)/(2*L) * ln(1/target_R)` uses the `n=2` prefactor `3/2` for a cubic (`n=3`) profile, so the realized
reflection is `R_eff = R_0^{3/4}` (e.g. `--pml-target-R 1e-3` → `≈5.6e-3`). This is documented in the source
(`pml_layer.cpp:61–67`), the banner prints `R_eff` (`spatial_dyn_driver.cpp:1720`), and the plan explicitly chose to
"keep 3/2 and relabel" (`§12.1 req 4` / `Eq.17 caveat`). It is **not** a bug — but the user-facing knob names
(`--pml-target-R`, `[numerics].pml_target_R`) still read like `R_0`, and the relabel lives only in a comment + one
banner token, so a user can under-specify absorption by ~an order of magnitude in the exponent.

**Suggested fix (optional):** print both in the banner so the contract is unmissable:
```diff
-                   << ", R_eff=" << cfg.numerics.pml_target_R
+                   << ", R_input=" << cfg.numerics.pml_target_R
+                   << " (R_eff~=" << std::pow(cfg.numerics.pml_target_R, 0.75) << ")"
```
or switch the prefactor to the true cubic `2.0` (`pml_layer.cpp:68`) — the four existing reflection/energy asserts are
`<`-bounds, so they stay green.

**Test case:** n/a (labeling). If the prefactor is changed to `2.0`, re-run `seas_test_pml` (asserts are inequalities).

---

### [R-003] [LOW] [SAFS config TOMLs] — stale comments claim "PML is NOT wired; leave use_pml=false"

**Category:** QUALITY (stale documentation)

**Description:**
At least two SAFS production TOMLs (e.g.
`spatial_friction_rate_state_safs_projected_stress_resolution_Dc2.toml:83–90`) still carry comments stating PML "is NOT
wired … leave `use_pml=false`." Post-Phase-12 this is false and misleading — a user reading the config would conclude PML
is unavailable and never try `use_pml=true` / `--pml`.

**Trigger:** A user configuring a SAFS run from the TOML.

**Suggested fix:** Update the stale comment block in the affected TOMLs to point at Phase 12 (`--pml` /
`[numerics].use_pml` + `pml_target_R`/`pml_cells`/`pml_damp_bottom`/`pml_damp_top`) and note the
`interior_flux="scalar"` + total-Q-background requirement.

**Test case:** n/a (comment text).

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001 — ADER-corrector PML path has no automated test; production path + plan's #1 risk)
- Low issues: 2 (R-002 `R_eff` labeling [documented]; R-003 stale TOML comments)
- Plan compliance: **FULL** for the code deliverables — 12.1 half-face mask + two new tests
  (`TestPMLFreeSurfaceTopUndamped`, `TestPMLFreeSurfacePulseStable`); 12.2 driver wiring (CLI, parser, derivation,
  face-mask, fault-clearance guard, banner, total-Q + scalar-material guards, global bbox reduction); 12.3 metric script +
  three sbatch variants (`pml_on`/`pml_off`/`bigbox_ref`). No feature is missing.
- Verdict: **PASS WITH FIXES** — no correctness defect found; the implementation (including the four hunted-and-cleared
  subtleties listed above) is sound. The one fix that matters before relying on `--pml` in production is R-001 (add the
  ADER-corrector test so the production path is regression-guarded); R-002/R-003 are optional polish.

## Unreviewed Areas
- **`spatial/code/scripts/pml_reflection_metrics.py` (Phase 12.3, ~403 lines)** — existence and correct CLI flag wiring
  from the sbatch jobs were confirmed, but the metric math (reflected/incident peak ratio, interior-energy windowing,
  big-box vs PML `V_max` comparison) was not line-audited. It is a post-processing/validation deliverable, off the
  simulation-correctness path, so it was de-prioritized behind the C++ damping path.
- **True multi-rank (np>1) PML behavior** — the MIN/MAX reductions (Verified-Correct #2/#3) are correct by inspection,
  but no test runs the PML path under MPI. Per repo policy ("no local full-mesh runs") this is validated on Frontera; a
  cheap 2-rank driver-banner/bbox-diff smoke would still catch a future reduction-op regression.
- **Phases 10/11/13/14** bundled in commit `30bf3e3` were out of scope (Phase 13 has a separate review at
  `miniapps/seas/REVIEW.md`).

## Note on review reliability
The interactive shell/Read channel dropped tool output for part of this session; every Verified-Correct item and every
finding above was ultimately confirmed against the committed code (clean ASCII dumps + a verification sub-agent that
quoted exact file:line and reported working-tree == HEAD for all three core files). An earlier interim draft of this
review (`REVIEW_phase12_pml.md`, now superseded) listed three of the Verified-Correct items as `[POSSIBLE]` bugs — they
are NOT bugs; do not act on that draft.
