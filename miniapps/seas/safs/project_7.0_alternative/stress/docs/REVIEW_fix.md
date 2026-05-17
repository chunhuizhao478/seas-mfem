# Fix Report — REVIEW.md (R-901 through R-906, all-phases sweep)

## Summary
- Findings addressed: **5 of 6 fixed; 1 explicitly deferred per review (R-905)**
- Files modified: **5**
  - `miniapps/seas/io/stress_field_3d.cpp` (R-903)
  - `miniapps/seas/io/field_coefficient.{hpp,cpp}` (R-904)
  - `miniapps/seas/tests/unit/test_stress_field_3d.cpp` (R-901)
  - `miniapps/seas/safs/project_7.0_alternative/code_preprocess/test_build_stress_safs.py` (R-902)
  - `miniapps/seas/safs/project_7.0_alternative/code_preprocess/data_projection_onfaultstress/PLAN_onfaultstress.md` (R-906)
- Tests added: **4** (1 C++ + 3 Python)
- Test suite:
  - **C++ Phase 6 Tranche 1: 56 / 56 passed** (was 46; +10 R-901 assertions in `T_6_11_stress_field_coefficient_eval`)
  - **Python Phase 1-5 combined: 173 / 173 passed** (was 170; +3 R-902 cross-phase contract tests)

## Changes Made

### R-901 [MODERATE] — `StressFieldCoefficient` unit test
**File modified:** `tests/unit/test_stress_field_3d.cpp`

- Added `#include "../../io/stress_field_coefficient.hpp"` at the top.
- New test `T_6_11_stress_field_coefficient_eval`:
  - Builds a unit-hex serial mesh located inside the synthetic
    sidecar bbox.
  - Constructs `StressFieldCoefficient` over the loaded
    `StressField3D`.
  - Calls `Eval(v, T, ip)` at the element center.
  - Asserts `GetVDim() == 6` and that each of the six entries
    matches the schema-v1 canonical order (xx, yy, zz, xy, yz,
    xz) — explicit annotation "NOT Voigt order" on the v(3)
    assertion to flag the most likely silent bug.
  - Re-runs with `scale=2.0, offset=100.0` and verifies uniform
    application of scale/offset to multiple components.
- Test wired into `main()`.
- **+10 new assertions; all pass.**

### R-902 [MODERATE] — Phase 5 ↔ Phase 6 cross-phase contract pytest
**File modified:** `code_preprocess/test_build_stress_safs.py`

- New `TestPhase5Phase6Contract` class with three tests:
  - `test_R902_phase5_field_names_match_phase6_reader`: builds a
    real `stress_safs.h5` via Phase 5 on a synthetic mesh, then
    asserts the on-disk field name set is exactly Phase 6's
    expected set (xx, yy, zz, xy, yz, xz). A field-name drift
    between phases is caught here.
  - `test_R902_phase5_field_units_are_Pa`: pins every component's
    `units` attribute to `"Pa"` (R-501/R-502 sign-convention
    contract relies on Pa-valued components).
  - `test_R902_schema_v1_attrs_present`: verifies
    `schema_version`, `crs`, `z_positive` are present and
    canonical so `DataField3D` (which Phase 6 delegates to) can
    load the file.
- The class hard-codes Phase 6's expected six-field tuple as a
  module-level constant so the test fails loudly if either side
  of the contract drifts.
- **+3 new tests; all pass.**

### R-903 [LOW] — `BBox()` returns the explicit intersection
**File modified:** `io/stress_field_3d.cpp`

- Added `#include <algorithm>`.
- Ctor body rewritten to compute the strict intersection of all
  six component bboxes (per plan §1 line 1644-1645) rather than
  copying `sigma_xx_.BBox()`. The behaviour under the schema-v1
  contract is identical (all six are bit-exact equal), but the
  intersection-explicit form is robust against future schema
  relaxation that allows per-component grids.
- The `AssertConsistentGrid_()` call still runs after `bbox_` is
  computed, so any actual cross-component mismatch still aborts.

### R-904 [LOW] — `ProjectStress` accepts `scale` / `offset`
**Files modified:** `io/field_coefficient.hpp`, `io/field_coefficient.cpp`

- Header: added `real_t scale = 1.0, real_t offset = 0.0`
  parameters to `ProjectStress`. Docstring extended to note the
  uniform application.
- Implementation: forwards `scale` / `offset` to each of the six
  inner `Project(field.Field(c), target_fes, scale, offset)`
  calls.
- Default values preserve the prior pass-through behaviour
  exactly; existing call sites (the test path) need no
  changes.

### R-905 [LOW] — Deferred per review
The review explicitly flagged R-905 as defer-to-future
("Recommended: defer to a future perf-only refactor of
`FieldProjector`. Flag as LOW."). No fix applied this cycle.
The six redundant mesh-bbox `MPI_Allreduce` calls per
`ProjectStress` invocation are wasted work but not a correctness
issue.

### R-906 [LOW] — Plan sign-flip attribution corrected
**File modified:** `PLAN_onfaultstress.md`

- Plan §1 line 1601-1607: replaced the misattribution
  (`build_stress_safs.py:evaluate_stress_field_on_grid` as the
  sign-flip site) with the correct citation
  (`project_to_fault_stress.py:bulk_stress_tensor_field` —
  Phase 3 §2). Added a clarifying note that Phase 5's
  `evaluate_stress_field_on_grid` is the unit-conversion
  pass-through (MPa → Pa).
- Inline R-906 reconciliation note tags the change.
- `PHASE6_DEVIATIONS.md` already correctly attributes the flip
  to Phase 3 at line 45 (verified via grep) — no edit required.
- The source `stress_field_3d.hpp` docstring lines 10-16 was
  already correct from the original implementation — no change.

## Unresolved Findings
- **R-905** — deferred per the review's own recommendation;
  flagged for a future perf-only refactor of `FieldProjector`.

## New Tests
- `T_6_11_stress_field_coefficient_eval` (C++) — covers R-901;
  guards against Voigt-vs-schema-v1 ordering regressions in
  `StressFieldCoefficient::Eval`.
- `test_R902_phase5_field_names_match_phase6_reader` (Python) —
  pins the Phase 5 → Phase 6 field-name contract.
- `test_R902_phase5_field_units_are_Pa` (Python) — pins the
  Pa-units contract on every component.
- `test_R902_schema_v1_attrs_present` (Python) — pins the
  schema-v1 root attrs so Phase 6's `DataField3D` ctor
  succeeds.

## Verification
- [x] R-901 (MODERATE): unit test for `StressFieldCoefficient` added; 10 assertions pass.
- [x] R-902 (MODERATE): three cross-phase contract tests added; all pass.
- [x] R-903 (LOW): bbox intersection computed explicitly in `stress_field_3d.cpp` ctor.
- [x] R-904 (LOW): `scale` / `offset` parameters added to `ProjectStress` with defaults that preserve prior behaviour.
- [ ] R-905 (LOW): deferred per review.
- [x] R-906 (LOW): plan attribution corrected; deviations doc was already correct.

## Ready for Re-Review: YES

C++ suite: 56 / 56 passing. Python suite: 173 / 173 passing.
All MODERATE findings closed with regression-guard tests; all
LOW findings either fixed or explicitly deferred per the
review's own guidance.
