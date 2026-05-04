# SAFS mesh tests — legacy-test policy (P-013)

Phase 4 of `PLAN_cgal_corefine.md` migrated the production
conformalizer from the Python implementation
(`fault_intersect.py` + `conformalize_faults.py`) to the C++
CGAL tool (`tools/corefine_faults`).  The Python sources are
preserved as a regression fixture but no longer on the
production hot-path.

## Running the tests

By default, `pytest` skips the legacy Python conformalizer tests:

```bash
conda activate pythonenv
pytest miniapps/seas/safs/mesh/tests/                 # default; skips legacy
pytest miniapps/seas/safs/mesh/tests/ -m legacy       # runs legacy
pytest miniapps/seas/safs/mesh/tests/ -m "not legacy" # explicit non-legacy
```

The default behaviour is enforced by `addopts = -m "not legacy"`
in `pytest.ini`.

The legacy suite is **required to pass at every release** as long
as the retention clock below has not expired.

## 90-day retention clock

- **Clock starts**: at the Phase-4 merge date (the date this
  README first lands on the main branch).
- **Clock duration**: 90 days.
- **Reset condition**: any actual CGAL-pipeline regression that the
  legacy tests catch resets the clock to a fresh 90 days.
- **Removal eligibility**: after 90 days, if no CGAL regression has
  surfaced AND the legacy tests have not caught any
  production-relevant bug, the legacy tests + Python sources
  (`fault_intersect.py`, `conformalize_faults.py`) become removable
  in a follow-up ticket.

Removal is **not** automatic — it is a human decision triggered by
the clock expiring without incident.  If unsure, keep them.

## What lives where

| File | Status |
|---|---|
| `test_fault_intersect.py` | LEGACY (`pytest.mark.legacy`) |
| `test_conformalize_faults.py` | LEGACY (`pytest.mark.legacy`) |
| `pytest.ini` | configuration; gates legacy by default |
| `tools/tests/test_corefine_smoke.cpp` | C++ smoke (current production); run via `ctest` in `tools/build/` |

The C++ smoke tests covering the production path are NOT in this
directory (they live with the C++ source under `tools/tests/`).
