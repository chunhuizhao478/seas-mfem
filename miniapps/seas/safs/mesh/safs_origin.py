"""SAFS mesh origin convention and CFM ↔ local-frame mapping.

Single source of truth for the local-frame definition used throughout the SAFS
mesh pipeline.  See miniapps/seas/safs/PLAN_origin.md for the design rationale.

Sign convention (memorize this):

    v_local = v_utm - origin_utm
    v_utm   = v_local + origin_utm

with the production origin at UTM Zone 11N (NAD83) (500_000, 3_765_000, 0).
"""
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

# ---------------------------------------------------------------------------
# Production origin constants.  Do NOT inline these values in any other file
# (see PLAN_origin.md acceptance criterion: grep gate).  The origin is a
# *named geographic reference frame*, not a tunable simulation parameter.
# ---------------------------------------------------------------------------
UTM_ZONE: str = "11N"
DATUM: str = "NAD83"
E0: float = 500_000.0
N0: float = 3_765_000.0
Z0: float = 0.0
ORIGIN_UTM: tuple[float, float, float] = (E0, N0, Z0)

# ---------------------------------------------------------------------------
# Canonical CFM short-name table.  Used by every script that needs to map a
# CFM filename to a stable identifier.  Keys are the CFM-ID stem (without the
# "_500m"/"_1000m"/"_2000m" suffix and without ".ts"); values are the short
# names used for STL filenames, --include-fault flag arguments, and
# fault_provenance.json keys.
# ---------------------------------------------------------------------------
FAULT_SHORT_NAMES: dict[str, str] = {
    "SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6":
        "safs_mjvs_saf",
    "ETRA-PMFZ-MULT-Pinto_Mountain_fault-CFM5":
        "safs_pmfz_pinto",
    "SAFS-SAFZ-SBMT-Mill_Creek_fault_strand-CFM4":
        "safs_sbmt_millcreek",
    "SAFS-SAFZ-SBMT-Mission_Creek_fault_strand-CFM4":
        "safs_sbmt_missioncreek",
    "SAFS-SAFZ-COAV-Mission_Creek_fault_strand-CFM4":
        "safs_coav_missioncreek",
    "SAFS-SAFZ-MULT-Banning_fault-CFM6":
        "safs_mult_banning",
    "SAFS-SAFZ-MULT-Southern_San_Andreas_fault_and_Banning-CFM6":
        "safs_mult_ssaf_banning",
    "SAFS-SAFZ-SBMT-San_Andreas_fault-CFM6":
        "safs_sbmt_saf",
    "SAFS-SAFZ-SBMT-Garnet_Hill_fault-CFM6":
        "safs_sbmt_garnethill",
}

_RES_SUFFIX_RE = re.compile(r"_(?P<res>\d+)m$")


def cfm_id_from_path(path: Path | str) -> str:
    """Strip directory, .ts extension, and resolution suffix."""
    p = Path(path)
    stem = p.stem
    m = _RES_SUFFIX_RE.search(stem)
    if m:
        stem = stem[: m.start()]
    return stem


def short_name_from_path(path: Path | str) -> str:
    cfm_id = cfm_id_from_path(path)
    if cfm_id not in FAULT_SHORT_NAMES:
        raise KeyError(
            f"unknown CFM id {cfm_id!r} for path {path}; "
            f"add it to FAULT_SHORT_NAMES in safs_origin.py"
        )
    return FAULT_SHORT_NAMES[cfm_id]


def _as_origin(o: Iterable[float] | None) -> tuple[float, float, float]:
    if o is None:
        return ORIGIN_UTM
    o = tuple(float(x) for x in o)
    if len(o) != 3:
        raise ValueError(f"origin must be a 3-tuple, got {len(o)}")
    return o  # type: ignore[return-value]


def utm_to_local(
    v_utm: np.ndarray,
    origin: Iterable[float] | None = None,
) -> np.ndarray:
    """v_local = v_utm - origin.

    Accepts (3,), (N, 3) or (..., 3) arrays.  Returns same shape, float64.
    """
    o = _as_origin(origin)
    arr = np.asarray(v_utm, dtype=np.float64)
    if arr.shape == (3,):
        return arr - np.asarray(o, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[-1] != 3:
        raise ValueError(
            f"expected last axis size 3, got shape {arr.shape}"
        )
    return arr - np.asarray(o, dtype=np.float64).reshape(
        (1,) * (arr.ndim - 1) + (3,)
    )


def local_to_utm(
    v_local: np.ndarray,
    origin: Iterable[float] | None = None,
) -> np.ndarray:
    """v_utm = v_local + origin.  Inverse of utm_to_local."""
    o = _as_origin(origin)
    arr = np.asarray(v_local, dtype=np.float64)
    if arr.shape == (3,):
        return arr + np.asarray(o, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[-1] != 3:
        raise ValueError(
            f"expected last axis size 3, got shape {arr.shape}"
        )
    return arr + np.asarray(o, dtype=np.float64).reshape(
        (1,) * (arr.ndim - 1) + (3,)
    )


# ---------------------------------------------------------------------------
# transform.json read/write with schema validation.
# ---------------------------------------------------------------------------
_SCHEMA_PATH = Path(__file__).resolve().parent / "transform.schema.json"


def _load_schema() -> dict:
    with open(_SCHEMA_PATH, "r") as fh:
        return json.load(fh)


def _validate_against_schema(data: Mapping) -> None:
    """Validate `data` against transform.schema.json.

    Uses jsonschema if available; falls back to a hand-coded check on the
    required keys and constants from the schema.
    """
    schema = _load_schema()
    try:
        import jsonschema  # type: ignore
    except ImportError:
        jsonschema = None  # type: ignore

    if jsonschema is not None:
        jsonschema.validate(instance=data, schema=schema)
        return

    # Fallback: minimal check.
    for key in schema["required"]:
        if key not in data:
            raise ValueError(
                f"transform.json missing required key {key!r}"
            )
    if data.get("version") != 1:
        raise ValueError(
            f"transform.json version must be 1, got {data.get('version')!r}"
        )
    if not re.match(r"^[0-9]{1,2}[NS]$", str(data["utm_zone"])):
        raise ValueError(
            f"utm_zone {data['utm_zone']!r} does not match pattern"
        )
    o = data["origin_utm_m"]
    if not (isinstance(o, list) and len(o) == 3):
        raise ValueError("origin_utm_m must be a 3-element list")
    for v in o:
        if not isinstance(v, (int, float)):
            raise ValueError("origin_utm_m entries must be numbers")
    conv = data["convention"]
    if conv.get("x_axis") != "UTM_easting":
        raise ValueError("convention.x_axis must equal 'UTM_easting'")
    if conv.get("y_axis") != "UTM_northing":
        raise ValueError("convention.y_axis must equal 'UTM_northing'")
    if conv.get("z_axis") != "elevation":
        raise ValueError("convention.z_axis must equal 'elevation'")
    if conv.get("z_sign") != "positive_up":
        raise ValueError("convention.z_sign must equal 'positive_up'")


def _current_git_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            cwd=Path(__file__).resolve().parent,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


def write_transform_json(
    path: Path | str,
    origin: Iterable[float] | None = None,
    extra: Mapping | None = None,
) -> None:
    o = _as_origin(origin)
    data: dict = {
        "version": 1,
        "utm_zone": UTM_ZONE,
        "datum": DATUM,
        "origin_utm_m": [float(o[0]), float(o[1]), float(o[2])],
        "convention": {
            "x_axis": "UTM_easting",
            "y_axis": "UTM_northing",
            "z_axis": "elevation",
            "z_sign": "positive_up",
        },
        "comment": (
            "v_local = v_utm - origin_utm_m. Right-handed (E,N,Up) metres. "
            "Free surface at z=0."
        ),
        "generated_by": "safs_origin.write_transform_json",
        "git_commit": _current_git_commit(),
    }
    if extra:
        for k, v in extra.items():
            data[k] = v
    _validate_against_schema(data)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, sort_keys=False)
        fh.write("\n")


def read_transform_json(path: Path | str) -> dict:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found; run ts_to_stl.py first to generate it"
        )
    with open(p, "r") as fh:
        data = json.load(fh)
    if "version" not in data:
        raise ValueError(
            f"{p} predates schema v1 (no 'version' field); "
            f"regenerate with the current ts_to_stl.py"
        )
    _validate_against_schema(data)
    return data


__all__ = [
    "UTM_ZONE",
    "DATUM",
    "E0",
    "N0",
    "Z0",
    "ORIGIN_UTM",
    "FAULT_SHORT_NAMES",
    "cfm_id_from_path",
    "short_name_from_path",
    "utm_to_local",
    "local_to_utm",
    "write_transform_json",
    "read_transform_json",
]
