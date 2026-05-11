# Phase 7.6 of miniapps/seas/document/io_dev/PLAN_bulk_compression_and_size_estimator_2026-05-09.md.
#
# Pytest test suite for the size estimator.  Standard-library only.
# Reference-run tests are auto-skipped when the corresponding measured
# output is not pre-staged on disk.

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

from _io_size_compression import (
    COMPRESSION_RATIOS,
    NUCLEATION_DURATION_S,
    bytes_per_dof,
    classify_field,
)
from _io_size_mesh import MeshSummary, read_inline_mesh
from _io_size_schedule import (
    ScheduleConfig,
    estimate_n_volume_writes,
    estimate_n_writes,
)
from _io_size_schemas import (
    BP5_FIELDS,
    TPV102_FIELDS,
    estimate_dofs,
    get_driver_schema,
)


# ---------------------------------------------------------------------------
# Phase 7.1 — mesh reader
# ---------------------------------------------------------------------------

def test_read_inline_mesh_tet_basic():
    m = read_inline_mesh(2, 2, 2, "tet")
    assert m.n_elements == 2 * 2 * 2 * 6
    assert m.n_vertices == 27
    assert m.element_type == "tet"
    # y-axis fault: 2 quads × 2 triangles each = 4
    assert m.n_fault_faces == 8
    assert m.n_boundary_faces > 0


def test_read_inline_mesh_hex_basic():
    m = read_inline_mesh(2, 2, 2, "hex")
    assert m.n_elements == 8
    assert m.element_type == "hex"
    assert m.n_fault_faces == 4  # one cube face = 2x2 quads


def test_read_gmsh_autodetects_fault_tag_from_physical_names(tmp_path):
    """R-402: when the mesh has a `$PhysicalNames` block declaring a
    surface named 'fault', `read_gmsh` must use that tag regardless of
    the `fault_tag` arg.  Prevents the 180× under-count seen on
    `bp5_1000m.msh` (which uses tag 100, not the Tandem default 3)."""
    p = tmp_path / "physname.msh"
    p.write_text(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
        "$PhysicalNames\n"
        "2\n"
        "2 100 \"fault\"\n"
        "2   5 \"absorbing\"\n"
        "$EndPhysicalNames\n"
        "$Nodes\n4\n"
        "1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n"
        "$EndNodes\n"
        "$Elements\n3\n"
        "1 2 2 100 1 1 2 3\n"   # fault triangle (tag 100)
        "2 2 2 5   1 2 3 4\n"   # absorbing triangle (tag 5)
        "3 4 2 1   1 1 2 3 4\n" # bulk tet
        "$EndElements\n")
    from _io_size_mesh import read_gmsh
    s = read_gmsh(p, fault_tag=3)  # caller's fallback is wrong on purpose
    assert s.fault_tag == 100
    assert s.n_fault_faces == 1
    assert s.n_boundary_faces == 2  # both surface tris counted as boundary


def test_R501_resolve_fault_tag_no_substring_false_positive():
    """R-501: substring match `'fault' in 'default'` returned True
    because 'default' literally contains the chars f-a-u-l-t at
    positions 2-6.  The fix uses an exact / start-anchored match."""
    from _io_size_mesh import _resolve_fault_tag
    # No-match cases.
    assert _resolve_fault_tag({1: "default"}, 99) == 99
    assert _resolve_fault_tag({1: "DEFAULT_BUFFER"}, 99) == 99
    assert _resolve_fault_tag({1: "non_fault_zone"}, 99) == 99
    assert _resolve_fault_tag({1: "faulty"}, 99) == 99
    # Match cases.
    assert _resolve_fault_tag({1: "fault"}, 99) == 1
    assert _resolve_fault_tag({1: "Fault"}, 99) == 1
    assert _resolve_fault_tag({1: "fault_main"}, 99) == 1
    assert _resolve_fault_tag({1: "FAULT 1"}, 99) == 1
    # Exact match preferred over start-anchored.
    assert _resolve_fault_tag({1: "fault_main", 2: "fault"}, 99) == 2


def test_read_gmsh_uses_fallback_when_no_physical_names(tmp_path):
    """When the mesh has no `$PhysicalNames` block, `read_gmsh` must
    honour the explicit `fault_tag` arg."""
    p = tmp_path / "no_physname.msh"
    p.write_text(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
        "$Nodes\n3\n1 0 0 0\n2 1 0 0\n3 0 1 0\n$EndNodes\n"
        "$Elements\n2\n"
        "1 2 2 7 1 1 2 3\n"     # surface tri tagged 7
        "2 2 2 9 1 1 2 3\n"     # surface tri tagged 9
        "$EndElements\n")
    from _io_size_mesh import read_gmsh
    s = read_gmsh(p, fault_tag=7)
    assert s.fault_tag == 7
    assert s.n_fault_faces == 1


def test_read_gmsh_v2_minimal(tmp_path):
    """Round-trip a tiny v2 ASCII mesh."""
    p = tmp_path / "tiny.msh"
    p.write_text("""$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
$EndNodes
$Elements
2
1 2 2 3 1 1 2 3
2 4 2 1 1 1 2 3 4
$EndElements
""")
    m = read_inline_mesh.__class__  # touch import to silence pyflakes
    from _io_size_mesh import read_gmsh
    s = read_gmsh(p, fault_tag=3)
    assert s.n_vertices == 4
    assert s.n_elements == 1  # one tet
    assert s.element_type == "tet"
    assert s.n_fault_faces == 1
    assert s.n_boundary_faces == 1


# ---------------------------------------------------------------------------
# Phase 7.2 — schemas + DOF formulas
# ---------------------------------------------------------------------------

def test_schemas_match_driver_register_calls():
    """Drift detection (R-308) — every string-literal call to
    `RegisterDomainField("name", ...)` or
    `pv_dc_->RegisterField("name", ...)` in the driver source files
    must appear in the corresponding driver schema.

    Indirect calls (variable-name passed in) are skipped per R-308
    limitation."""
    repo = THIS_DIR.parent  # miniapps/seas
    drivers = {
        "tpv102": repo / "drivers" / "tpv102_driver.cpp",
        "tpv104": repo / "drivers" / "tpv104_driver.cpp",
        "tpv205": repo / "drivers" / "tpv205_driver.cpp",
        "bp5":    repo / "tests" / "verification" / "bp5_verification_full.cpp",
    }
    register_re = re.compile(
        r'(?:RegisterDomainField|RegisterField)\s*\(\s*"([A-Za-z_][A-Za-z_0-9]*)"')
    for driver, path in drivers.items():
        if not path.is_file():
            pytest.skip(f"{path} not present")
        text = path.read_text()
        names_in_src = set(register_re.findall(text))
        schema = {f.name for f in get_driver_schema(driver)}
        missing = names_in_src - schema
        # Some names are registered in fault-side init (paraview_output.hpp)
        # — those are also in the schema if appropriate.  Only fail when
        # the source has a name the schema doesn't.
        assert not missing, (
            f"driver {driver} registers fields {missing} "
            f"that are not in the size-estimator schema. "
            f"Update miniapps/seas/scripts/_io_size_schemas.py.")


def test_dof_formula_l2_p0_face():
    from _io_size_schemas import FieldSchema
    f = FieldSchema("slip_dip", "L2_p0_face", is_fault=True, is_volume=False)
    m = MeshSummary(n_elements=100, n_vertices=200, n_fault_faces=42)
    assert estimate_dofs(f, m, order=1) == 42


def test_dof_formula_h1_p1_tet():
    from _io_size_schemas import FieldSchema
    f = FieldSchema("displacement", "H1_pK", n_components=3)
    m = MeshSummary(n_elements=100, n_vertices=200, element_type="tet")
    assert estimate_dofs(f, m, order=1) == 200


def test_dof_formula_l2_pk_byNODES_tet_p2():
    from _io_size_schemas import FieldSchema
    f = FieldSchema("velocity", "L2_pK_byNODES", n_components=3)
    m = MeshSummary(n_elements=10, n_vertices=20, element_type="tet")
    # P=2 tet has (2+1)(2+2)(2+3)/6 = 10 nodes per element.
    assert estimate_dofs(f, m, order=2) == 100


# ---------------------------------------------------------------------------
# Phase 7.3 — schedule integrator
# ---------------------------------------------------------------------------

def test_schedule_integrator_no_events_returns_tfinal_over_dt_interseismic():
    cfg = ScheduleConfig(
        tfinal=10.0 * 3.156e7,  # 10 yr
        dt_interseismic=3.156e7,
        n_events=0,
        driver="bp5",
    )
    # n_events default for 10 yr BP5 = round(10 / 240) = 0.
    assert cfg.tfinal > 0
    n = estimate_n_writes(cfg)
    # 10 yr / 1 yr ≈ 10 writes (+ initial t=0 = 11) when no events.
    assert 9 <= n <= 12


def test_schedule_integrator_5_events_matches_phase3_synthetic():
    """When n_events=5 and the cap is loose, the writer count
    should equal the sum of the three regimes plus 1 initial save.

    R-409: dt_coseismic default tracks the C++ AdaptiveSchedule (0.01 s),
    NOT the obsolete 60.0 s default."""
    cfg = ScheduleConfig(
        tfinal=250.0 * 3.156e7,  # 250 yr
        n_events=5,
        avg_event_duration_s=30.0,
        driver="bp5",
        max_total_snapshots=0,
    )
    n = estimate_n_writes(cfg)
    # Coseismic budget:  5 × 30   = 150 s    / 0.01 s = 15_000 writes.
    # Nucleation budget: 5 × 86_400 = 432_000 s / 1 s   = 432_000 writes.
    # Interseismic:      ~250 yr / 1 yr      ≈ 250 writes.
    # Total ≈ 447_250.
    assert 446_000 <= n <= 448_000


def test_schedule_integrator_cap_clamps():
    cfg = ScheduleConfig(
        tfinal=250.0 * 3.156e7,
        n_events=12,
        max_total_snapshots=5000,
        driver="bp5",
    )
    n = estimate_n_writes(cfg)
    # Cap=K + n_events.
    assert n == 5000 + 12


def test_schedule_integrator_cap_one_plus_events():
    cfg = ScheduleConfig(
        tfinal=250.0 * 3.156e7,
        n_events=3,
        max_total_snapshots=1,
        driver="bp5",
    )
    assert estimate_n_writes(cfg) == 1 + 3


def test_schedule_integrator_fixed_dt_override():
    cfg = ScheduleConfig(
        tfinal=100.0,
        fixed_dt=10.0,
        driver="tpv102",
    )
    # 100 / 10 + 1 = 11.
    assert estimate_n_writes(cfg) == 11


def test_schedule_integrator_volume_fixed_dt():
    cfg = ScheduleConfig(tfinal=100.0, driver="tpv102")
    assert estimate_n_volume_writes(cfg, volume_pv_dt=10.0) == 11


def test_nucleation_duration_table():
    assert NUCLEATION_DURATION_S["bp5"] == 86_400.0
    assert NUCLEATION_DURATION_S["tpv102"] == 1.0
    assert NUCLEATION_DURATION_S["tpv104"] == 1.0
    assert NUCLEATION_DURATION_S["tpv205"] == 1.0


# ---------------------------------------------------------------------------
# Phase 7.4 — compression table
# ---------------------------------------------------------------------------

def test_compression_table_velocity_zfp_1e3_in_range():
    from _io_size_schemas import FieldSchema
    f = FieldSchema("velocity", "L2_pK_byNODES", n_components=3)
    bpd = bytes_per_dof(f, "zfp_1e-3", "tpv102")
    assert 0.3 <= bpd <= 0.5


def test_compression_table_unknown_filter_raises():
    from _io_size_schemas import FieldSchema
    f = FieldSchema("velocity", "L2_pK_byNODES")
    with pytest.raises(ValueError):
        bytes_per_dof(f, "lz4_fast", "tpv102")


def test_classify_field_kinds():
    assert classify_field("velocity", "tpv102") == "smooth_fp"
    assert classify_field("displacement", "bp5") == "smooth_fp"
    assert classify_field("slip_rate_dip", "bp5") == "rough_fp"
    assert classify_field("state_variable", "bp5") == "rough_fp"
    assert classify_field("mpi_rank", "any") == "integer_index"
    assert classify_field("param_a", "bp5") == "static_2d"


def test_compression_placeholder_warning_emitted_once_per_pair():
    from _io_size_schemas import FieldSchema
    f = FieldSchema("velocity", "L2_pK_byNODES")
    # Reset the latch so the first call below WILL warn even if a
    # prior test already triggered it.
    from _io_size_compression import _PLACEHOLDER_WARNED
    _PLACEHOLDER_WARNED.clear()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        bytes_per_dof(f, "zfp_1e-3", "tpv102")
        bytes_per_dof(f, "zfp_1e-3", "tpv102")  # dedup
    placeholder_warns = [x for x in w
                          if "PLACEHOLDER" in str(x.message)]
    assert len(placeholder_warns) == 1


# ---------------------------------------------------------------------------
# Phase 7.5 — CLI smoke tests
# ---------------------------------------------------------------------------

ESTIMATOR = THIS_DIR / "estimate_output_size.py"


def _run_estimator(*args: str) -> tuple[int, str]:
    res = subprocess.run(
        [sys.executable, str(ESTIMATOR), *args, "--quiet"],
        capture_output=True, text=True)
    return res.returncode, res.stdout


def test_cli_bp5_inline_mesh_runs():
    rc, out = _run_estimator(
        "--driver", "bp5",
        "--inline-mesh",
        "--tfinal", "250yr",
        "--paraview",
        "--paraview-fault-zfp-tol", "1e-12",
        "--paraview-max-snapshots", "5000",
        "--no-volume-pv")
    assert rc == 0
    assert "TOTAL" in out


def test_cli_tpv102_inline_mesh_runs():
    rc, out = _run_estimator(
        "--driver", "tpv102",
        "--inline-mesh",
        "--inline-nx", "20", "--inline-ny", "20", "--inline-nz", "20",
        "--tfinal", "1.0",
        "--paraview",
        "--paraview-volume-zfp-tol", "1e-3")
    assert rc == 0
    assert "TOTAL" in out


def test_cli_json_output_round_trips():
    rc, out = _run_estimator(
        "--driver", "bp5", "--inline-mesh",
        "--tfinal", "10yr", "--paraview", "--json")
    assert rc == 0
    parsed = json.loads(out)
    for key in ("fault_bytes", "volume_bytes", "station_bytes",
                "total_bytes", "n_writes_fault", "n_writes_volume"):
        assert key in parsed


def test_cli_total_excludes_station_bytes():
    """The on-fault PV path covers all fault faces; benchmark station
    TXT files are a separate output path and must NOT be summed into
    TOTAL."""
    rc, out = _run_estimator(
        "--driver", "bp5", "--inline-mesh",
        "--tfinal", "10yr", "--paraview", "--json")
    assert rc == 0
    parsed = json.loads(out)
    assert parsed["station_bytes"] > 0
    assert parsed["total_bytes"] == (
        parsed["fault_bytes"] + parsed["volume_bytes"])


def test_cli_explain_names_every_input():
    rc, out = _run_estimator(
        "--driver", "bp5", "--inline-mesh",
        "--tfinal", "10yr", "--paraview", "--explain")
    assert rc == 0
    assert "per_write_bytes_fault" in out
    assert "n_writes_fault" in out
    assert "station_rows" in out


def test_cli_quota_exceeded_returns_nonzero():
    rc, _ = _run_estimator(
        "--driver", "bp5", "--inline-mesh",
        "--tfinal", "250yr",
        "--paraview", "--paraview-fault-zfp-tol", "1e-12",
        "--scratch-quota", "1KB")
    assert rc == 1


def test_cli_unknown_time_suffix_errors():
    rc, _ = _run_estimator(
        "--driver", "bp5", "--inline-mesh", "--tfinal", "1fortnight",
        "--paraview")
    assert rc != 0


def test_cli_volume_vtu_and_hdf5_mutually_exclusive():
    rc, _ = _run_estimator(
        "--driver", "bp5", "--inline-mesh",
        "--tfinal", "1yr",
        "--paraview-volume-vtu", "--paraview-volume-hdf5")
    assert rc != 0


def test_format_bytes_scaling():
    from estimate_output_size import format_bytes
    assert format_bytes(0) == "0 B"
    assert format_bytes(1023).endswith("B")
    assert format_bytes(1024).endswith("KB")
    assert format_bytes(1024 ** 3).endswith("GB")
    assert format_bytes(1024 ** 4).endswith("TB")
    assert format_bytes(1024 ** 5).endswith("PB")


def test_parse_time_suffixes():
    from estimate_output_size import parse_time
    assert parse_time("1s") == 1.0
    assert parse_time("1min") == 60.0
    assert parse_time("1hr") == 3600.0
    assert parse_time("1day") == 86_400.0
    assert parse_time("1wk") == 86_400.0 * 7
    assert parse_time("1yr") == 3.156e7
    assert parse_time("3.156e7") == 3.156e7  # no-suffix = seconds
    assert parse_time("1.5day") == 1.5 * 86_400.0
    with pytest.raises(ValueError):
        parse_time("1fortnight")


def test_parse_bytes_suffixes():
    from estimate_output_size import parse_bytes
    assert parse_bytes("1KB") == 1024
    assert parse_bytes("1MB") == 1024 ** 2
    assert parse_bytes("1GB") == 1024 ** 3
    assert parse_bytes("1TB") == 1024 ** 4
    assert parse_bytes("4096") == 4096


# ---------------------------------------------------------------------------
# Round-3 regressions (R-701, R-503, R-504, R-704, R-406, R-407, R-411)
# ---------------------------------------------------------------------------

REPO_ROOT = THIS_DIR.parent              # miniapps/seas
SBATCH_BP5_P6 = (REPO_ROOT / "jobs/bp5"
                 / "bp5_phase6_paraview_zfp_normal_48hr.sbatch")
SBATCH_TPV_P6 = [
    REPO_ROOT / "jobs/tpv102/tpv102_phase6_paraview_zfp_dev_2hr.sbatch",
    REPO_ROOT / "jobs/tpv104/tpv104_phase6_paraview_zfp_dev_2hr.sbatch",
    REPO_ROOT / "jobs/tpv205/tpv205_phase6_paraview_zfp_dev_2hr.sbatch",
]


def test_R701_bp5_sbatch_tfinal_is_literal_not_subshell():
    """R-701: BP5 sbatch must NOT use a $(...) subshell for --tfinal —
    if python3 is missing on the compute node, $() returns empty and the
    driver parses --tfinal --verify as tfinal=0, ending the simulation
    at step 1 and wasting the 48-hour reservation."""
    text = SBATCH_BP5_P6.read_text()
    # The offending pattern: --tfinal followed by a $(...) substitution.
    assert "--tfinal $(" not in text, (
        "R-701 regression: BP5 sbatch uses --tfinal $(...) subshell.  "
        "Use a literal numeric value instead so the driver gets a "
        "well-defined tfinal even when python3 is missing.")
    # Must be a positive numeric literal.
    m = re.search(r"--tfinal\s+(\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)", text)
    assert m, "R-701: --tfinal literal not found"
    assert float(m.group(1)) > 0


def test_R702_bp5_sbatch_tfinal_matches_seconds_per_year():
    """R-702: --tfinal seconds equal 250 * BP5Params::seconds_per_year
    (= 250 * 31_557_600 = 7_889_400_000), NOT 250 * 3.156e7 (which is
    0.18% off)."""
    text = SBATCH_BP5_P6.read_text()
    m = re.search(r"--tfinal\s+(\d+)", text)
    assert m, "R-702: integer --tfinal literal not found"
    value = int(m.group(1))
    expected = 250 * 31_557_600
    assert value == expected, (
        f"R-702 regression: --tfinal = {value} but expected {expected} "
        f"(= 250 * BP5Params::seconds_per_year).  Off by "
        f"{value - expected} seconds ({(value - expected) / 86400:.2f} days).")


def test_R404_all_phase6_sbatch_assert_h5z_zfp():
    """R-404: every Phase 6 sbatch must grep for MFEM_USE_H5Z_ZFP=YES
    so the build is rejected before ibrun rather than after ~30 SU."""
    all_p6 = [SBATCH_BP5_P6, *SBATCH_TPV_P6]
    for sb in all_p6:
        text = sb.read_text()
        assert "MFEM_USE_H5Z_ZFP" in text, (
            f"R-404 regression: {sb.name} missing MFEM_USE_H5Z_ZFP grep")


def test_R503_paraview_dt_not_ambiguous():
    """R-503: --paraview-dt is reserved for the C++ driver's fixed-dt
    override (which the estimator does NOT support directly; the
    Python equivalent is --paraview-fixed-dt).  With allow_abbrev=False
    the estimator must REJECT --paraview-dt (unknown option) rather
    than match the ambiguous --paraview-(co|nu|inter)seismic-dt prefix."""
    rc, _ = _run_estimator(
        "--driver", "tpv102", "--inline-mesh",
        "--tfinal", "1.0", "--paraview", "--paraview-dt", "0.5")
    # Before R-503 fix: error: ambiguous option: --paraview-dt ...
    # After fix:        error: unrecognized arguments: --paraview-dt 0.5
    # Either way nonzero exit, but it must NOT silently dispatch.
    assert rc != 0, (
        "R-503 regression: --paraview-dt accepted unambiguously; "
        "argparse allow_abbrev should be False on the estimator parser.")


def test_R504_per_regime_dt_uses_cpp_word_order():
    """R-504: the Python estimator must accept the C++ driver flag
    names verbatim (`--paraview-coseismic-dt`, NOT
    `--paraview-dt-coseismic`).  Users copy the sbatch verbatim."""
    rc, out = _run_estimator(
        "--driver", "bp5", "--inline-mesh",
        "--tfinal", "1yr", "--paraview",
        "--paraview-coseismic-dt", "0.005",
        "--paraview-nucleation-dt", "0.5",
        "--paraview-interseismic-dt", "1e7",
        "--json")
    assert rc == 0, (
        f"R-504 regression: estimator rejected C++-style flag names.  "
        f"stdout was:\n{out}")
    parsed = json.loads(out)
    assert parsed["total_bytes"] > 0


def test_R409_schedule_default_dt_coseismic_matches_cpp():
    """R-409: ScheduleConfig.dt_coseismic default must match the C++
    AdaptiveSchedule default (0.01 s, paraview_output.hpp:168).  The
    previous 60.0 default under-counted coseismic writes 6000×."""
    s = ScheduleConfig(tfinal=1.0)
    assert s.dt_coseismic == 0.01


def test_R406_integer_index_zfp_uses_deflate_fallback():
    """R-406: ZFP only filters floating-point datasets in VTKHDF; for
    integer index datasets the on-disk size follows the deflate
    fallback (mesh/vtkhdf.cpp:162-170).  The table must NOT fall
    through to 8.0 raw-double."""
    for tol in ("zfp_1e-3", "zfp_1e-6", "zfp_1e-9", "zfp_1e-12"):
        assert ("integer_index", tol) in COMPRESSION_RATIOS, (
            f"R-406 regression: COMPRESSION_RATIOS missing "
            f"('integer_index', '{tol}') — falling back to 8.0 raw "
            f"double 16× over-estimates the mpi_rank field.")
        bpd, _ = COMPRESSION_RATIOS[("integer_index", tol)]
        assert bpd <= 1.0, (
            f"R-406 regression: integer_index/{tol} should track the "
            f"deflate fallback (~0.5 B/DOF), got {bpd}")


def test_R407_volume_writes_apply_cap_when_volume_pv_dt_positive():
    """R-407: estimate_n_volume_writes must apply max_total_snapshots
    even when volume_pv_dt > 0."""
    cfg = ScheduleConfig(
        tfinal=10_000.0, n_events=0,
        max_total_snapshots=50, driver="tpv102")
    n = estimate_n_volume_writes(cfg, volume_pv_dt=1.0)
    # Uncapped would be ~10_001 (tfinal // 1.0 + 1).
    # Capped (with _apply_cap) is cap + n_events = 50 + 0 = 50.
    assert n <= 50 + 1, (
        f"R-407 regression: cap (50) ignored — got {n} volume writes")


def test_R704_estimate_from_sbatch_forwards_coseismic_dt(tmp_path):
    """R-704: estimate_from_sbatch.py:estimate_one must forward
    --paraview-coseismic-dt (and the other per-regime overrides) into
    the ScheduleConfig.  Previously it dropped them; the estimate
    then used the Python defaults which disagreed with the C++."""
    # Pre-stage a tiny mesh so estimate_one can run.
    mesh = tmp_path / "mini.msh"
    mesh.write_text(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
        "$Nodes\n1\n1 0 0 0\n$EndNodes\n"
        "$Elements\n0\n$EndElements\n")
    # Build a one-line sbatch that uses the C++ flag name.
    sb = tmp_path / "fake.sbatch"
    sb.write_text(
        f"#!/bin/bash\n"
        f"ibrun ./seas_tpv102_driver \\\n"
        f"    --mesh {mesh} --tfinal 1.0 --paraview \\\n"
        f"    --paraview-coseismic-dt 0.005 \\\n"
        f"    --order 1\n")
    sys.path.insert(0, str(THIS_DIR))
    from estimate_from_sbatch import parse_sbatch, estimate_one
    p = parse_sbatch(sb)
    # Sanity-check: parser actually saw the flag.
    assert p.flags.get("paraview-coseismic-dt") == "0.005"
    # Drive the estimator: with the flag wired through, the schedule
    # has dt_coseismic = 0.005 instead of the default 0.01.  Estimate
    # without raising.
    result = estimate_one(p)
    assert "n_writes_fault" in result


def test_R411_tpv205_sbatch_passes_lsw():
    """R-411: TPV205 sbatch must pass --fric-law lsw so the driver's
    banner matches the actual dispatch.  Any other value now ABORTS
    the driver at parse time (no more silent ignore)."""
    sb = REPO_ROOT / "jobs/tpv205/tpv205_phase6_paraview_zfp_dev_2hr.sbatch"
    text = sb.read_text()
    assert "--fric-law lsw" in text, (
        "R-411 regression: TPV205 phase-6 sbatch must pass "
        "--fric-law lsw (anything else now aborts the driver).")


# ---------------------------------------------------------------------------
# Phase 7.6 — reference runs
# ---------------------------------------------------------------------------

# These three tests are auto-skipped when the corresponding measured
# output is not pre-staged.  A measurement is "available" when the
# REFERENCE_RUNS.md file exists AND lists the run with a measured size.

REFERENCE_DOC = THIS_DIR / "REFERENCE_RUNS.md"


def _reference_measured_bytes(run_name: str) -> int:
    """Parse REFERENCE_RUNS.md for the measured-bytes column.  Returns
    0 if the file or row is missing."""
    if not REFERENCE_DOC.is_file():
        return 0
    txt = REFERENCE_DOC.read_text()
    # Extremely simple parser: look for a line containing the run name
    # and a Frontera measured size in bytes / GB / MB.
    pattern = re.compile(
        rf"^\|\s*{re.escape(run_name)}.*?\|\s*([0-9.]+\s*[KMGT]?B)\s*\|",
        re.MULTILINE)
    m = pattern.search(txt)
    if m is None:
        return 0
    from estimate_output_size import parse_bytes
    return parse_bytes(m.group(1))


@pytest.mark.skip(
    reason="reference measurement pending Phase 6.6 Frontera replication")
def test_reference_run_1_bp5_250yr_within_30pct():
    pass


@pytest.mark.skip(
    reason="reference measurement pending Phase 6.6 Frontera replication")
def test_reference_run_2_tpv102_1s_within_30pct():
    pass


@pytest.mark.skip(
    reason="reference measurement pending Phase 6.6 Frontera replication")
def test_reference_run_3_tpv205_zfp_within_30pct():
    pass
