"""Unit tests for check_mesh_tool_log.py.

Run with:
    conda activate pythonenv
    pytest miniapps/seas/safs/mesh/tests/test_check_mesh_tool_log.py -v
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_PARENT = _HERE.parent
sys.path.insert(0, str(_PARENT))

import check_mesh_tool_log as cmtl  # noqa: E402


# ---------------------------------------------------------------------------
# scan() — pattern-match unit tests
# ---------------------------------------------------------------------------
def test_clean_gmsh_log_returns_no_hits():
    log = (
        "Info    : Reading 'safs.geo'...\n"
        "Info    : Done reading\n"
        "Info    : Done meshing 3D (Wall 12.3s)\n"
        "Info    : 1099208 tets\n"
    )
    assert cmtl.scan(log, cmtl._PATTERN_SETS["gmsh"]) == []


def test_hxt_3d_failure_is_caught():
    log = (
        "Info    : Meshing 3D...\n"
        "Error   : HXT 3D mesh failed\n"
        "Warning : No elements in volume 1\n"
    )
    hits = cmtl.scan(log, cmtl._PATTERN_SETS["gmsh"])
    descs = [d for _, d in hits]
    assert any("HXT 3D mesh failed" in d for d in descs), descs
    assert any("no elements" in d.lower() for d in descs), descs


def test_self_intersecting_facets_caught():
    log = (
        "Info    : Done meshing 1D\n"
        "Found two exactly self-intersecting facets (dihedral 0.0).\n"
    )
    hits = cmtl.scan(log, cmtl._PATTERN_SETS["gmsh"])
    assert any("self-intersecting" in d for d in (h[1] for h in hits))


def test_failed_to_recover_constrained_caught():
    log = (
        "Info    : 3D mesh in progress\n"
        "Info    : failed to recover constrained lines/triangles\n"
    )
    hits = cmtl.scan(log, cmtl._PATTERN_SETS["gmsh"])
    assert any("recover constrained" in d for _, d in hits)


def test_mmg3d_bad_ending_caught():
    log = (
        "** MMG3DLIB BAD ENDING OF MMG3DLIB:\n"
        "##  invalid mesh\n"
    )
    hits = cmtl.scan(log, cmtl._PATTERN_SETS["mmg3d"])
    assert any("failed to terminate" in d for _, d in hits)


def test_mmg3d_error_token_caught():
    log = (
        "MMG3D - LIBRARY MODE\n"
        "## ERROR: bad geometric input\n"
    )
    hits = cmtl.scan(log, cmtl._PATTERN_SETS["mmg3d"])
    assert any("## ERROR" in d for _, d in hits)


def test_python_traceback_caught_for_safs_tools():
    log = (
        "Traceback (most recent call last):\n"
        '  File "foo.py", line 1, in <module>\n'
        "    raise ValueError('bad')\n"
        "ValueError: bad\n"
    )
    hits = cmtl.scan(log, cmtl._PATTERN_SETS["safs"])
    assert any("uncaught Python" in d for _, d in hits)


def test_assertion_error_caught():
    log = (
        "Running step 6...\n"
        "AssertionError: invariant violated\n"
    )
    hits = cmtl.scan(log, cmtl._PATTERN_SETS["safs"])
    assert any("Python assertion" in d for _, d in hits)


def test_gmsh_does_not_falseflag_error_summary_zero():
    """gmsh prints harmless 'Error : 0' summaries on clean runs.

    These should NOT trigger the gate; the pattern requires the
    canonical fatal banners ('HXT 3D mesh failed', 'Mesh generation
    error summary', 'No elements in volume') with their full text.
    """
    log = (
        "Info    : Done meshing\n"
        "Info    : 0 errors\n"
        "Error   : 0\n"   # gmsh's count-line, not a real error
    )
    hits = cmtl.scan(log, cmtl._PATTERN_SETS["gmsh"])
    descs = [d for _, d in hits]
    # Must NOT match the HXT / mesh-summary / no-elements patterns.
    fatal_descs = [
        "HXT 3D mesh failed",
        "no elements in a tagged volume",
        "gmsh terminal mesh-generation error summary",
    ]
    for fd in fatal_descs:
        assert not any(fd in d for d in descs), (
            f"false positive: {fd} matched on clean log {descs}")


# ---------------------------------------------------------------------------
# main() — exit-code integration tests
# ---------------------------------------------------------------------------
def test_main_returns_0_on_clean_log(tmp_path: Path):
    log = tmp_path / "clean.log"
    log.write_text("Info    : all done\n")
    rc = cmtl.main(["--tool", "gmsh", "--log", str(log)])
    assert rc == 0


def test_main_returns_1_on_dirty_log(tmp_path: Path, capsys):
    log = tmp_path / "dirty.log"
    log.write_text("Error   : HXT 3D mesh failed\n")
    rc = cmtl.main(["--tool", "gmsh", "--log", str(log)])
    assert rc == 1
    err = capsys.readouterr().err
    assert "HXT 3D mesh failed" in err
    assert "fatal pattern" in err


def test_main_returns_2_on_missing_log(tmp_path: Path, capsys):
    rc = cmtl.main(["--tool", "gmsh",
                     "--log", str(tmp_path / "nope.log")])
    assert rc == 2
    err = capsys.readouterr().err
    assert "log file not found" in err


def test_main_tool_choice_validation(tmp_path: Path):
    log = tmp_path / "x.log"
    log.write_text("ok\n")
    with pytest.raises(SystemExit):
        cmtl.main(["--tool", "nonsense", "--log", str(log)])


def test_main_tool_all_runs_every_pattern_set(tmp_path: Path):
    """`--tool all` must catch any of the gmsh + mmg3d + safs patterns."""
    log = tmp_path / "mixed.log"
    log.write_text(
        "Error   : HXT 3D mesh failed\n"
        "## ERROR: mmg3d aborted\n"
        "Traceback (most recent call last):\n"
    )
    rc = cmtl.main(["--tool", "all", "--log", str(log)])
    assert rc == 1
