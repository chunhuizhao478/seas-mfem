"""TPV104 probe-diff tooling tests (§4.10 Step 13 gates).

Dual-mode runnable:
    python3 test_probe_diff.py             # stdlib-only (no pytest)
    python3 -m pytest test_probe_diff.py   # when pytest is available

Coverage:
    - Identical-by-construction probe file pairs diff to zero
      (every probe, every channel).
    - `tpv104_column_map.apply_reference_to_mfem` swaps strike/dip
      channels.
    - `tpv104_column_map.apply_reference_to_mfem` flips the sign of
      `sigma_n_*` channels.
    - `probe_diff.parse_probe_file` rejects files with wrong column counts.
    - `align_by_t_qp` drops unmatched rows deterministically.
"""

from __future__ import annotations

import os
import sys
import tempfile
import traceback
from typing import Dict, List, Tuple

import numpy as np

# Tests run from tpv104/scripts/tests/; put the scripts on the path.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                "..")))

import probe_diff
import tpv104_column_map as colmap


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def write_probe(path: str, probe_name: str, rows: List[List[float]]) -> None:
    column_names = probe_diff.PROBE_COLUMNS[probe_name]
    with open(path, "w") as fh:
        fh.write(f"# probe={probe_name}  code=MFEM  rank=0  nprocs=1\n")
        fh.write("# t qp_id " + " ".join(column_names) + "\n")
        for r in rows:
            fh.write(" ".join(f"{x:.16e}" for x in r) + "\n")


class _TmpDir:
    def __init__(self) -> None:
        self._td = tempfile.TemporaryDirectory()
        self.path = self._td.name

    def __call__(self, name: str) -> str:
        return os.path.join(self.path, name)

    def cleanup(self) -> None:
        self._td.cleanup()


# ---------------------------------------------------------------------------
# Test 1 — identical-by-construction pair diffs to zero (all probes).
# ---------------------------------------------------------------------------
def test_identical_pair_zero_diff_all_probes() -> None:
    tmp = _TmpDir()
    try:
        for probe_name in probe_diff.PROBE_COLUMNS.keys():
            ncols = 2 + len(probe_diff.PROBE_COLUMNS[probe_name])
            rows: List[List[float]] = []
            for k in range(5):
                row = [k * 0.01, k] + [1.0 + k * 0.1 + j * 0.5
                                       for j in range(ncols - 2)]
                rows.append(row)
            p_mfem = tmp(f"mfem_{probe_name}.txt")
            p_ref  = tmp(f"ref_{probe_name}.txt")
            write_probe(p_mfem, probe_name, rows)
            write_probe(p_ref,  probe_name, rows)

            result = probe_diff.diff_probes(probe_name, p_mfem, p_ref,
                                            apply_map=False)
            for name, r in result.items():
                assert r["count"] == 5, \
                    f"{probe_name}:{name} aligned 5 rows (got {r['count']})"
                assert r["max_abs"] == 0.0, \
                    f"{probe_name}:{name} identical max_abs != 0 " \
                    f"({r['max_abs']})"
                assert r["max_rel"] == 0.0, \
                    f"{probe_name}:{name} identical max_rel != 0 " \
                    f"({r['max_rel']})"
    finally:
        tmp.cleanup()


# ---------------------------------------------------------------------------
# Test 2 — apply_reference_to_mfem swaps tau1 ↔ tau2.
# ---------------------------------------------------------------------------
def test_column_map_strike_dip_swap() -> None:
    ref_fields = {
        "tau1_trial":    np.array([10.0, 20.0, 30.0]),
        "tau2_trial":    np.array([ 1.0,  2.0,  3.0]),
        "sigma_n_trial": np.array([-100.0, -100.0, -100.0]),
    }
    mapped = colmap.apply_reference_to_mfem("trial_traction", ref_fields)
    np.testing.assert_array_equal(mapped["tau1_trial"],
                                   ref_fields["tau2_trial"])
    np.testing.assert_array_equal(mapped["tau2_trial"],
                                   ref_fields["tau1_trial"])


# ---------------------------------------------------------------------------
# Test 3 — apply_reference_to_mfem sign-flips sigma_n_*.
# ---------------------------------------------------------------------------
def test_column_map_sign_flip() -> None:
    ref_fields = {
        "tau1_trial":    np.array([0.0, 0.0]),
        "tau2_trial":    np.array([0.0, 0.0]),
        "sigma_n_trial": np.array([-1.2e8, -1.2e8]),
    }
    mapped = colmap.apply_reference_to_mfem("trial_traction", ref_fields)
    np.testing.assert_array_equal(mapped["sigma_n_trial"],
                                   np.array([+1.2e8, +1.2e8]))


# ---------------------------------------------------------------------------
# Test 4 — parse_probe_file rejects bad column counts.
# ---------------------------------------------------------------------------
def test_parse_probe_file_rejects_bad_column_count() -> None:
    tmp = _TmpDir()
    try:
        path = tmp("bad.txt")
        with open(path, "w") as fh:
            fh.write("# probe=trial_traction  code=MFEM  rank=0  nprocs=1\n")
            # trial_traction expects 5 cols; we emit 4 to force an error.
            fh.write("0.0 0 1.0 2.0\n")
        raised = False
        try:
            probe_diff.parse_probe_file(path, "trial_traction")
        except ValueError:
            raised = True
        assert raised, "parse_probe_file must raise on column-count mismatch"
    finally:
        tmp.cleanup()


# ---------------------------------------------------------------------------
# Test 5 — align_by_t_qp drops unmatched rows deterministically.
# ---------------------------------------------------------------------------
def test_align_drops_unmatched() -> None:
    tmp = _TmpDir()
    try:
        rows_m = [[0.01, 0, 1.0, 0.0, 0.0],
                  [0.02, 0, 2.0, 0.0, 0.0],
                  [0.03, 0, 3.0, 0.0, 0.0]]
        rows_r = [[0.01, 0, 10.0, 0.0, 0.0],
                  [0.03, 0, 30.0, 0.0, 0.0]]   # no row at t=0.02
        p_m = tmp("m.txt")
        p_r = tmp("r.txt")
        write_probe(p_m, "trial_traction", rows_m)
        write_probe(p_r, "trial_traction", rows_r)

        result = probe_diff.diff_probes("trial_traction", p_m, p_r,
                                        apply_map=False)
        assert result["tau1_trial"]["count"] == 2
        assert result["sigma_n_trial"]["count"] == 2
        # |1 - 10| = 9, |3 - 30| = 27  →  max_abs = 27
        assert result["sigma_n_trial"]["max_abs"] == 27.0
    finally:
        tmp.cleanup()


# ---------------------------------------------------------------------------
# Test 6 — applying the map twice returns to the original values.
# (Sanity: the column map is self-inverse for the swap and involutive
# for the sign flip.)
# ---------------------------------------------------------------------------
def test_column_map_self_inverse_for_trial_traction() -> None:
    ref_fields = {
        "tau1_trial":    np.array([11.0, 22.0]),
        "tau2_trial":    np.array([ 3.0,  4.0]),
        "sigma_n_trial": np.array([-100.0, -50.0]),
    }
    once  = colmap.apply_reference_to_mfem("trial_traction", ref_fields)
    twice = colmap.apply_reference_to_mfem("trial_traction", once)
    np.testing.assert_array_equal(twice["tau1_trial"],
                                   ref_fields["tau1_trial"])
    np.testing.assert_array_equal(twice["tau2_trial"],
                                   ref_fields["tau2_trial"])
    np.testing.assert_array_equal(twice["sigma_n_trial"],
                                   ref_fields["sigma_n_trial"])


# ---------------------------------------------------------------------------
# Stdlib-only runner (no pytest dependency).
# ---------------------------------------------------------------------------
_TESTS = [
    test_identical_pair_zero_diff_all_probes,
    test_column_map_strike_dip_swap,
    test_column_map_sign_flip,
    test_parse_probe_file_rejects_bad_column_count,
    test_align_drops_unmatched,
    test_column_map_self_inverse_for_trial_traction,
]


def _run_stdlib() -> int:
    passed = 0
    failed = 0
    for t in _TESTS:
        name = t.__name__
        try:
            t()
            print(f"  PASSED: {name}")
            passed += 1
        except AssertionError as e:
            print(f"  FAILED: {name}: {e}")
            failed += 1
        except Exception as e:          # noqa: BLE001 — report everything
            print(f"  ERROR:  {name}: {type(e).__name__}: {e}")
            traceback.print_exc()
            failed += 1
    print()
    print("=" * 60)
    print(f"  Results: {passed} passed, {failed} failed out of "
          f"{len(_TESTS)} tests")
    print("=" * 60)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(_run_stdlib())
