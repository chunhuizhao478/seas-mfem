# Code Review: Phase 7 — BP5 benchmark parity gate — 2026-06-04

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` §Phase 7.
- Files reviewed:
  - `drivers/spatial_seas_driver.cpp` (the bp5-native rate-state fill branch).
  - `spatial/code/spatial_friction.hpp` + `.cpp` (`RateStateBlock.bp5_analytic` flag + parse).
  - `config/safs_qd/bp5_spatial_seas_verification.toml` + `..._cgamg.toml`.
  - `tests/verification/compare_bp5_parity.py` (parity harness).
- Domain context: CLAUDE.md (4-phase init, sign rules, BP5 regression list), the plan (R-007), project memory ([[no-local-reproducer]], [[project_qd_safs_coupling_needs_clean_fault_basis]]), the implementer completion report.
- Method: three adversarial passes; verified `ComputeBP5Params` (the golden fill) field-by-field against the bp5-native fill, the `RateStatePerDOFParams` field list, and `SetRateStatePerDOF`'s reads (it does NOT read `rs.V_w`/`b`/`f_0`/`V_0`, so the unset `V_w` is harmless — no crash). The bp5-native PATH is not locally runnable (production mesh), so correctness rests on the field-by-field `ComputeBP5Params` match.

## Findings

### [R-701] MODERATE [BUG] tests/verification/compare_bp5_parity.py:compare_pair — NaN/Inf differences silently PASS

**Category:** BUG

**Description:**
`compare_pair` flags a failure with `if d > atol and r > rtol`. If either trace value is `NaN` or `Inf` (e.g., a blown-up run writes `nan`/`-inf` to `log10(V)`), `d = abs(a-b)` is `NaN`, and `NaN > atol` is **False** in Python — so the condition is never true and the row does NOT fail. A run that diverged (NaN/Inf) would therefore PASS the parity gate against a clean golden run. This is the opposite of what a validation gate must do.

**Trigger:** one run's station trace contains `nan`/`inf` (the other is finite).

**Actual behavior:** `overall_ok` stays True → exit 0 (PASS).

**Expected behavior:** any finite-vs-nonfinite mismatch (or nan in either) fails the comparison.

**Suggested fix:** treat non-finite mismatch as a hard failure in the per-element loop:
```diff
             a, b = rows_a[i][c], rows_b[i][c]
+            if math.isnan(a) != math.isnan(b) or math.isinf(a) != math.isinf(b) \
+               or (math.isnan(a) and math.isnan(b) is False):
+                ok = False
+                max_abs[c] = float("inf"); max_rel[c] = float("inf")
+                continue
+            if math.isnan(a) or math.isnan(b) or math.isinf(a) or math.isinf(b):
+                ok = False
+                max_abs[c] = float("inf"); max_rel[c] = float("inf")
+                continue
             d = abs(a - b)
```
(Simpler equivalent: at the top of the inner loop, `if not (math.isfinite(a) and math.isfinite(b)): if a != b (or nan-xor): ok=False; continue`.)

**Test case:**
```python
def test_R701_nan_fails():
    # rows_a finite, rows_b has nan in a column
    a = [[0,0,0,-9,-30,30.0,0,1.5]]
    b = [[0,0,0,float('nan'),-30,30.0,0,1.5]]
    ok, *_ = compare_pair(a, b, rtol=1e-6, atol=1e-9)
    assert ok is False   # pre-fix: passes (True)
```

---

### [R-702] MODERATE [POSSIBLE] drivers/spatial_seas_driver.cpp — no guard that the two bp5_analytic flags are set TOGETHER

**Category:** ASSUMPTION

**Description:**
BP5 parity requires BOTH `[stress].bp5_analytic` AND `[friction.rate_state].bp5_analytic`: the prestress `tau0_vec` is computed from `bp5_params.a(x2,x3)` internally, so the rate-state `a` MUST also come from `bp5_params`. If a config sets only ONE flag (e.g., `[stress].bp5_analytic=true` with a resolver-sourced uniform `a`), the prestress and friction `a` are inconsistent → `SetInitialCondition` aborts on equilibrium `> 1e-6` with a message about equilibrium, not about the flag mismatch. Nothing warns at config time.

**Trigger:** a config with exactly one of the two `bp5_analytic` flags set.

**Suggested fix:** after parsing both, warn/abort on the half-set case (rank 0):
```cpp
const bool rs_bp5 = cfg.rate_state.has_value() && cfg.rate_state->bp5_analytic;
if (rs_bp5 != cfg.stress.bp5_analytic)
{
   MFEM_ABORT("spatial_seas: [friction.rate_state].bp5_analytic ("
              << rs_bp5 << ") and [stress].bp5_analytic ("
              << cfg.stress.bp5_analytic << ") must be set TOGETHER — the BP5 "
              "prestress uses bp5_params.a(x2,x3), so the rate-state a must too "
              "(R-007).  Mixing resolver a with bp5_analytic prestress breaks "
              "equilibrium.");
}
```
(Place it near where `cfg` is finalized, before the time loop.)

**Test case:**
```
// test_R702_half_set_bp5_analytic_rejected: a config with [stress].bp5_analytic
//   = true but [friction.rate_state].bp5_analytic = false must abort at config
//   validation with the "set TOGETHER" message (not at SetInitialCondition).
```

---

### [R-703] LOW drivers/spatial_seas_driver.cpp:(bp5-native fill) — V_init direction is assumed pure +strike, not asserted

**Category:** ASSUMPTION

**Description:**
The fill stores `rs.V_init(i) = sqrt(Vi[0]^2 + Vi[1]^2)` (the magnitude) and relies on the hardcoded `init_vel_dir = (0, 1)` to reconstruct `V_init_vec = (0, +mag)` in `SetRateStatePerDOF`. This reproduces `ComputeBP5Params`'s `V_init_vec_ = (Vi[0], Vi[1])` ONLY if `Vi = (0, +Vp)` (pure positive strike — true for BP5). If `bp5_params.V_init_vec` ever returned a dip component or a negative strike, the reconstruction would silently differ (the sign/dip is dropped). Correct for BP5, but unasserted.

**Suggested fix:** assert the BP5 assumption when filling:
```cpp
real_t Vi[2];
seed.V_init_vec(cx2(i), cx3(i), Vi);
MFEM_VERIFY(std::abs(Vi[0]) <= 1e-30 && Vi[1] >= 0.0,
            "spatial_seas: bp5_analytic rate-state assumes pure +strike V_init; "
            "got Vi=(" << Vi[0] << ", " << Vi[1] << ") at DOF " << i << ".");
rs.V_init(i) = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);
```

---

### [R-704] LOW tests/verification/compare_bp5_parity.py — relative tolerance applied to log10 columns

**Category:** QUALITY

**Description:**
`compare_pair` applies the same `rtol` to every column, including `log10(V_*)` and `log10(state)`. A relative tolerance on a log-magnitude quantity (e.g., `-30`) means `30·rtol` absolute, which is an odd unit choice — for tight parity it works (the diffs are ~round-off), but a per-column scheme (absolute tol for log columns, relative for slip/tau) is more meaningful and avoids surprising sensitivity. Not a correctness bug for the intended `<1e-6` parity, but worth a per-column tolerance map.

**Suggested fix:** add an optional `--abs-cols` set (default the two `log10V` + `log10state` indices) compared with `atol` only, the rest with `rtol`.

---

### [R-705] LOW tests/verification/compare_bp5_parity.py:compare_pair — compares by row index, not by time

**Category:** EDGE_CASE

**Description:**
The comparison is row-by-row (`rows_a[i]` vs `rows_b[i]`). This is only valid if both runs wrote at the SAME times (same dt sequence). For byte-parity (same solver path) they do; but if the two runs ever sampled different times (e.g., a different output cadence), row `i` of A and row `i` of B are different physical times → a meaningless comparison. The `[ROW MISMATCH]` note only fires on a total-count difference, not on time-misalignment within the common prefix.

**Suggested fix:** verify column-0 (time) agreement per row before comparing the rest, or match rows by nearest time:
```python
if abs(rows_a[i][0] - rows_b[i][0]) > 1e-6 * (abs(rows_b[i][0]) + 1.0):
    ok = False  # time misalignment
```

---

### [R-706] LOW tests/verification/compare_bp5_parity.py:parse_station — malformed data rows silently skipped

**Category:** EDGE_CASE

**Description:**
`try: [float(t) ...] except ValueError: continue` treats ANY non-numeric line as a header to skip. A genuinely corrupt data row (e.g., `1e6 ERROR 0 ...`) is silently dropped rather than flagged, so a partially-corrupt station file could compare "clean" on its surviving rows.

**Suggested fix:** distinguish the header (only the first non-numeric line) from later malformed rows; raise/warn if a non-numeric line appears after numeric data has begun.

---

## Summary
- Critical issues: **0** — the bp5-native fill matches `ComputeBP5Params` field-by-field (a_of_x2_x3, Dc_of_x2_x3, V_init_vec, eta, uniform b/f0/V0/sigma_n) in geom's DOF order; the unset `rs.V_w` is harmless (`SetRateStatePerDOF` does not read it); build clean + config 40/40.
- Moderate issues: **2** — R-701 (harness silently passes NaN/Inf — a gate must not), R-702 (no guard that both bp5_analytic flags are set together).
- Low issues: **4** — R-703 (V_init direction unasserted), R-704 (rel-tol on log columns), R-705 (row-index vs time alignment), R-706 (malformed rows skipped).
- Plan compliance: **PARTIAL (by design)** — the configs + the bp5-native fill (the approved R-007 deviation) + the parity harness are implemented and locally verified (build, config parse, harness PASS/FAIL exit codes). The acceptance RUNS (mumps-vs-golden, cgamg-vs-mumps, np8) are Frontera (production mesh).
- Verdict: **PASS WITH FIXES** — R-701 must be fixed before the harness is trusted as a gate (a blown-up run would falsely PASS); R-702 is a cheap config guard; R-703–706 are LOW. No blocker for the local build/parse.

## Unreviewed Areas
- The actual BP5 parity (station traces vs `seas_driver` golden): Frontera — not locally runnable ([[no-local-reproducer]]).
- The SAFS-mode-vs-BP5-native **traction-formula** equivalence (the residual parity assumption): confirmed only on Frontera.
- `bp5_params.a_of_x2_x3`/`Dc_of_x2_x3`/`V_init_vec` themselves (existing, extreme-care `config/bp5_params.hpp`; out of scope — the fill consumes them).
