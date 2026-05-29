# Code Review: Phase 13 — Split scalar / matrix into separate WaveOperator classes (2026-05-28, fresh adversarial pass)

## Review Scope
- Plan: `document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md` §Phase 13
- Files reviewed:
  - `dynamic/wave_operator.hpp`, `dynamic/wave_operator.inl` (hooks, R-001 fault-site routing, bi-material extraction → scalar-only)
  - `dynamic/bimaterial_wave_operator.hpp`, `dynamic/bimaterial_wave_operator.inl` (new subclass)
  - `drivers/spatial_dyn_driver.cpp` (matrix-branch construction)
  - `tests/unit/test_bimaterial_wave_operator_parity.cpp` (new C-6), `tests/unit/test_phaseh_wave_operator_constant_parity.cpp`, `tests/unit/test_wave_operator.cpp` (repointed)
  - `Makefile`
- Domain context: `CLAUDE.md` (sign conventions, byte-exact contract, "no local full-mesh runs"), git tag `hrs-ref` (reference single-class impl), memory notes (GodunovFluxPool dedup floor; placeholder-leak class).

## What was verified correct (the fix agent should NOT touch this)
- **All 8 fault imposed-state sites** route through `FluxForElem_(elem_plus/elem_minus/e1/qa.local_elem)` and are byte-faithful to hrs-ref (`git diff hrs-ref HEAD -- wave_operator.inl` shows no `flux_.Interior(can_n` divergence). Each per-side **deposit** is materially consistent: `F_h_plus` is computed with `elem_plus`'s material AND deposited to the plus-side DOFs (symmetric for minus). No swap.
- **`(1,1,1)` sentinel containment:** the only bare-`flux_` material uses left in `wave_operator.inl` are the scalar-only hook bodies (`InteriorFaceFlux_`/`SharedInteriorFaceFlux_`/`ApplyElementJacobian_`, all overridden on the subclass), the dead `Ax_/Ay_/Az_` ctor `BuildJacobian` (referenced only in comments), and the scalar `ComputeMaxDt` return (overridden). **No virtual hook is invoked during base construction** (verified across the whole ctor, lines 24–587), so the `(1,1,1)` seed never reaches a cached/used quantity.
- **`GetFlux()`** is read only in `test_adjacent_triangle_fault_first_step_audit.cpp`, which uses a scalar `WaveOperator` — no matrix object exposes the poisoned `flux_`.
- **Acceptance grep** on `wave_operator.{hpp,inl}` returns only the mandated `UsesGodunovFluxPool()` accessor — the scalar class is bi-material-free.
- **Tripwire validated:** temporarily re-leaking one fault site to bare `flux_` made C-6a fail at field-scale rel 0.999 (finite, not NaN); reverted cleanly.

## Findings

### [R-001] [MODERATE] [Makefile / `test:` target] — C-6 test is built but never run by `make test`

**Category:** DEVIATION

**Description:**
Phase 13.2 ("Files to Modify: Makefile … link the new test; add to the aggregate `test:`") requires the C-6 parity test in the suite. The run rule `test-bimaterial-wave-operator-parity` (~Makefile:3856) and the build target (added to `SEQ_MINIAPPS`, ~864) exist, but the `make test` prerequisite chain (`test:` at line 3547) lists the sibling `test-phaseh-wave-operator-constant-parity` (line 3575) and `test-wave-operator` (3561) — it does NOT list `test-bimaterial-wave-operator-parity`. So `make all` builds C-6 but `make test` never executes it; the leak-class tripwire silently rots.

**Trigger:** `make test` (suite run) — C-6 is absent from the executed targets.

**Actual behavior:** C-6 built by `make all` (via `SEQ_MINIAPPS`) but not run by `make test`.

**Expected behavior:** `make test` runs `test-bimaterial-wave-operator-parity` with the other wave-operator tests.

**Suggested fix:** add the run target to the `test:` prerequisite list next to its sibling (in the `test:` block around line 3575):
```diff
       test-phaseh-wave-operator-constant-parity \
+      test-bimaterial-wave-operator-parity \
```

**Test case:**
```python
def test_R001_c6_in_make_test_chain():
    import re, pathlib
    mk = pathlib.Path("miniapps/seas/Makefile").read_text()
    block = re.search(r"\ntest:(.*?)(?=\n\t)", mk, re.S).group(1)   # line-continued prereq block
    assert "test-bimaterial-wave-operator-parity" in block, \
        "C-6 run target missing from the 'make test' chain"
```

---

### [R-002] [MODERATE] [POSSIBLE] [spatial_dyn_driver.cpp:~1020 reflection-time warning] — matrix path computes `cp = sqrt(0/0) = NaN`

**Category:** BUG

**Description:**
The R-107 reflection-time warning (just below the wave-operator construction) computes
`cp = sqrt((material.lambda_const + 2*material.mu_const) / material.rho_const)`.
On the **matrix** path `material` is `Mode::Coefficient` (the driver `MFEM_VERIFY`s `material.mode != Constant`). `MaterialField::MakeCoefficient` (`heterogeneous_material.cpp:21`) sets only the `*_coef` pointers and leaves `lambda_const = mu_const = rho_const = 0.0` (struct defaults, `heterogeneous_material.hpp:74-76`). So `cp = sqrt((0+0)/0) = sqrt(NaN) = NaN`; `(cp>0.0)` is false → `t_reflect = 0.0`; `0.0 < tfinal` fires the warning with a bogus `cp_max (0 s)`. The reflection warning is therefore always wrong ("0 s") on matrix runs.

**Pre-existing** (introduced with the Phase-9 matrix branch); Phase 13 only swapped the constructed type. **Diagnostic only** — no physics impact (`dt` comes from the virtual `wave.ComputeMaxDt`, which is correct). Flagged because the matrix branch is Phase 13's finalized surface and the warning is silently useless there.

**Trigger:** any `interior_flux="matrix"` run.

**Actual behavior:** prints `cp_max (0 s)` and warns nonsensically regardless of the true reflection time.

**Expected behavior:** compute `cp` from the actual matrix material, or skip the warning when `material.mode != Constant`.

**Suggested fix:** guard the warning on Constant material (minimal, safe):
```diff
   // R-107 reflection-time warning: compute min_box_dim / cp_max from
   // mesh bounding box + scalar material.
-  {
+  if (material.mode == MaterialField::Mode::Constant)
+  {
      const real_t cp = std::sqrt((material.lambda_const
                                   + 2.0 * material.mu_const)
                                   / material.rho_const);
      ...
   }
```

**Test case:**
```python
def test_R002_coefficient_material_const_fields_are_zero():
    # C++ unit (pseudocode): MakeCoefficient leaves *_const = 0, so the driver
    # must NOT derive cp from them on the matrix path.
    #   ConstantCoefficient l(3e10), m(3e10), r(2670);
    #   MaterialField f = MaterialField::MakeCoefficient(&l,&m,&r);
    #   assert f.lambda_const == 0.0 && f.rho_const == 0.0;          # current defaults
    #   real_t cp = std::sqrt((f.lambda_const+2*f.mu_const)/f.rho_const);
    #   assert std::isnan(cp);                                       # demonstrates the bug
    pass
```

---

### [R-003] [MODERATE] [test_bimaterial_wave_operator_parity.cpp:MaxRelDiff] — C-6 uses a field-scale metric, not the plan's per-component "to 1e-9"

**Category:** DEVIATION

**Description:**
Plan §Phase 13 Acceptance C-6 says "matches the scalar `WaveOperator` `Mult`/ADER **to 1e-9**", and C-5 uses a per-component metric `|a-b| / max(|a_i|,|b_i|,1.0)`. C-6 instead uses a **field-scale** metric `|a-b| / max_i(|a_i|,|b_i|)`. This is a deliberate, documented deviation: `GodunovFluxPool` stores `round_sig(·,6)` material (godunov_flux_pool.cpp:101-106), so the matrix path differs from scalar by ~1 ULP per element, which at near-equilibrium fault DOFs (`dQ/dt≈0` computed as a near-cancellation of O(1e8) terms) is a large *per-component* relative error but ~3e-10 (Mult)/~3e-8 (ADER) of the field magnitude. The per-component metric would false-fail; the field-scale metric is the correct measure for a mixed-magnitude field and still catches a `(1,1,1)` leak (field-scale rel ≈ 1, validated).

Justified and tripwire-validated, but a real departure from the plan's stated acceptance metric that must be recorded. Residual risk: the looser metric could mask a *subtle, small-magnitude* matrix-path error that only a true-bimaterial test (R-004) would expose.

**Trigger:** comparing C-6's metric to the plan's literal "to 1e-9 per component".

**Actual behavior:** field-scale-relative tolerance 1e-9 (passes at 3e-10/3e-8).

**Expected behavior (plan literal):** per-component 1e-9 (not achievable on 1e10-scale moduli due to the dedup floor — hence the justified deviation).

**Suggested fix:** No functional change required. Keep the field-scale metric; ensure the plan §Phase 13 Acceptance line is annotated that C-6 uses a field-scale metric with the round_sig rationale (the test header already documents it). Do NOT revert to the per-component metric — it would false-fail.

**Test case:**
```python
def test_R003_field_scale_metric_separates_leak_from_floor():
    # Already realized as the validated tripwire: poisoning one fault site to
    # bare flux_ gives field-scale rel ~0.999 (FAIL); correct code ~3e-10 (PASS).
    pass
```

---

### [R-004] [LOW] [test_bimaterial_wave_operator_parity.cpp] — homogeneous C-6 cannot detect a per-side / bimaterial fault-flux error

**Category:** EDGE_CASE

**Description:**
C-6 builds the bimaterial operator from a **homogeneous** Coefficient material. With equal material on both fault-adjacent elements, `FluxForElem_(elem_plus) == FluxForElem_(elem_minus)`, so the R-001 per-side selection (`elem1_on_plus ? e1 : e2`) is a **no-op** and is NOT exercised: an `e1`/`e2` swap or a genuine bimaterial asymmetry would pass C-6. Per-side correctness is currently guaranteed only by hrs-ref faithfulness (`git diff`), not by a test. This matches the plan's note that the per-side treatment is "the scaffold a bi-material fault would use" — an accepted gap, not a defect.

**Trigger:** a true bimaterial fault (different `(λ,μ,ρ)` across the fault) — no current test exercises it.

**Suggested fix (follow-up, non-blocking):** document the homogeneous scope in the C-6 header:
```diff
 // adds the fault.
+// NOTE: homogeneous material → the per-side elem_plus/elem_minus selection is a
+// no-op here; its correctness is verified by faithfulness to hrs-ref, not by C-6.
```

**Test case:** N/A (coverage/documentation note; LOW).

---

### [R-005] [LOW] [coverage] — substep fault-flux sites (4, 5) on the matrix path are not unit-tested

**Category:** EDGE_CASE

**Description:**
C-6 exercises `Mult` (fault site 1) and `AdvanceADER` (site 3) on a **serial** mesh. The substep dispatch path — site 4 (`ComputeADERFaceFluxRHS` per-QP, `elem_plus_qq`) and site 5 (`ComputeADERSharedFaceFluxRHS`, `qa.local_elem`, parallel-only) — is not exercised on the matrix operator by any unit test; site 5 also requires np>1. The R-001 routing there is verified by hrs-ref faithfulness and inspection only.

**Suggested fix (follow-up, non-blocking):** note the gap; a np>1 bimaterial substep regression case would close it.

**Test case:** N/A (coverage note; LOW).

---

## Summary
- Critical issues: 0
- Moderate issues: 3 (R-001 C-6 not run by `make test`; R-002 matrix reflection-warning NaN [pre-existing, diagnostic]; R-003 C-6 metric deviation [justified, validated])
- Low issues: 2 (R-004 homogeneous-only fault coverage; R-005 substep matrix path untested)
- Plan compliance: FULL — every Phase-13 requirement is implemented; R-003 is a justified deviation; scalar byte-exactness is construction-guaranteed (FluxForElem_==flux_, scalar hook bodies verbatim, identical `ComputeMaxDt` return) but not run-verified locally (per the "no local full-mesh runs" rule).
- Verdict: **PASS WITH FIXES** — the core split (fault-site routing, `(1,1,1)` containment, class extraction, byte-exactness) is correct and tripwire-validated. Fix R-001 (one line) so the tripwire runs in the suite; fix/guard R-002 (one line) so the matrix path stops emitting a NaN reflection warning; record R-003.

## Unreviewed Areas
- **Scalar TPV205/102/104 + BP5 byte-identical station/benchmark traces** — requires full-mesh runs (forbidden locally per CLAUDE.md / memory). Byte-exactness is preserved by construction but NOT run-verified; confirm on Frontera against a pre-Phase-13 checkpoint before production sign-off.
- **True heterogeneous (non-homogeneous-Coefficient) matrix correctness** (TPV31 / depth profile) — Phase 10 deliverable, out of Phase-13 scope.
- **np>1 matrix path** (`ExchangeBiMaterialNeighbours_` is a local-side stub, moved-not-fixed per the plan) — lateral-heterogeneous parallel runs are explicitly deferred.
