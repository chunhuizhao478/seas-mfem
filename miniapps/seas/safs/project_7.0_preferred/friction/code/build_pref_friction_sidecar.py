#!/usr/bin/env python3
"""Build the PREFERRED MFEM friction sidecar (friction_safs.h5, data_projection_v1)
from the SAME ASAGI NetCDF SeisSol consumes (safs_friction_thermal_case2.nc,
compound {rs_a, rs_srW}).

Phase 3 of
  miniapps/seas/safs/project_7.0_preferred/document/
    PLAN_thermal_case2_mixedflux_port_2026-07-08.md

Why this exists
---------------
The SAFS v3_4_x THERMAL decks zone rate-and-state `a` and the strong-rate-
weakening velocity `V_w` (`rs_srW`) by TEMPERATURE, using the SCEC Community
Thermal Model (CTM, shinevar2024).  Those are genuine 3-D fields; neither MFEM's
scalar defaults, nor its 1-D `depth_profile` CSVs, nor a `boxcar_taper` rule can
express them.  SeisSol reads `rs_a` / `rs_srW` through ASAGI; MFEM reads the
SAME two baked fields through `[friction.rate_state.sidecar]`, so both codes
sample identical friction and the a(T)/V_w(T) piecewise logic stays in exactly
one place (thermal/code/build_friction_nc_thermal.py).

Decision D5 (user, 2026-07-08): MFEM never reads NetCDF.  The friction sidecar
is the same file type as the VELOCITY sidecar -- a `data_projection_v1` HDF5 --
and is written by the same shared helper, `velocity/code/sidecar.py::write_sidecar`.
This script is therefore a sibling of `build_pref_velocity_sidecar.py`, NOT of
`csm_stress_nc_to_mfem_hdf5.py` (which hand-rolls its own HDF5 writer).

No sign flip, no unit conversion
--------------------------------
`rs_a` is dimensionless and `rs_srW` is m/s in both conventions.  Unlike the
stress converter there is nothing to negate.  The only transform is the storage
transpose (z,y,x) -> (x,y,z) that the schema-v1 layout requires.

Guards
------
`write_sidecar` already enforces (and RAISES, never warns): strictly-monotone
axes, no NaN, shape == (Nx, Ny, Nz), and every cell inside `field_bounds` --
i.e. the upstream build script's G2 range gate comes for free.  This script adds
the one guard the helper cannot know about: the G4 hypocentre anchor
(`a == 0.015`, `V_w == 0.05` at the on-fault hypocentre), checked by trilinear
sampling BEFORE writing.

Usage
-----
  python3 build_pref_friction_sidecar.py \
      --in  ~/Downloads/seisol_quakeworx/safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2/safs_friction_thermal_case2.nc \
      --out <repo>/miniapps/seas/safs/project_7.0_preferred/friction/results/thermal_case2_mfem/friction_safs.h5

The emitted `.h5` is large and gitignored -- document its provenance and scp it
to Expanse alongside `velocity_safs.h5` and the stress sidecar.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# The shared schema-v1 writer, as used by build_pref_velocity_sidecar.py.
REPO = Path(__file__).resolve().parents[4]          # .../miniapps/seas
ALT = REPO / "safs" / "project_7.0_alternative"
_SIDECAR_DIR = ALT / "velocity" / "code"
if not (_SIDECAR_DIR / "sidecar.py").is_file():
    # R-004: `parents[4]` is a positional assumption about this file's depth.
    # If the script is moved or copied (which is how these builders get used),
    # sys.path.insert would silently succeed and `import sidecar` would either
    # raise a bare ModuleNotFoundError or, worse, bind to SOME OTHER sidecar.py
    # on sys.path and emit a wrong-schema file.  Fail here, naming the path.
    raise SystemExit(
        "%s: cannot find the shared schema-v1 writer at\n"
        "    %s\n"
        "This script locates it relative to its own path "
        "(parents[4] must be miniapps/seas).\n"
        "Run it from its home in the repo, or fix REPO above."
        % (__file__.rsplit("/", 1)[-1], _SIDECAR_DIR / "sidecar.py")
    )
sys.path.insert(0, str(_SIDECAR_DIR))
from sidecar import write_sidecar  # noqa: E402

# SAFS v3_4_1 PREFERRED on-fault hypocentre: the closest point on the fault
# surface at EXACTLY 10 km depth to the deck-lineage geographic target.
# Matches safs_fault.yaml's Tnuc_s LuaMap centre and the production TOML's
# [nucleation.gradual_overstress_compact_circular] centre.
HYPO_DEFAULT = (609062.8722, 3709528.1324, -10000.0)

# CASE2 anchors at the hypocentre (T = 277.8 C, mid velocity-weakening zone).
# See thermal/docs/PLAN_thermal_friction_zoning_2026-07-01.md.
A_HYPO_EXPECTED = 0.015
VW_HYPO_EXPECTED = 0.05

# G2 range gate, enforced by write_sidecar via field_bounds (it RAISES on any
# out-of-range cell).  CASE2 spans a in [0.015, 0.0332] and V_w in [0.05, 1000].
#
# On the endpoints and float32: the ASAGI members are float32, so the exact
# constants round-trip as float32(0.015) = 0.01499999966... and
# float32(0.05) = 0.05000000074...  The `a` lower bound is 0.010, comfortably
# below the former.  The `V_w` lower bound can safely be the exact 0.05 because
# float32(0.05) rounds UP, never below.  1000.0 is exactly representable.
# Keep these TIGHT: a looser bound would let a genuinely wrong field through.
A_BOUNDS = (0.010, 0.050)
VW_BOUNDS = (0.05, 1000.0)

# ASAGI compound member names, in the order build_friction_nc_thermal.py writes.
NC_A_MEMBER = "rs_a"
NC_VW_MEMBER = "rs_srW"


# ---------------------------------------------------------------------------
# netCDF read
# ---------------------------------------------------------------------------
def load_friction_nc(path):
    """Return (x, y, z, {name: (Nx,Ny,Nz) float64}).

    The ASAGI file stores ONE compound variable `data` with dims (z, y, x);
    the schema-v1 sidecar wants one dataset per member with shape (Nx, Ny, Nz).
    """
    try:
        import netCDF4
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            "build_pref_friction_sidecar: needs the `netCDF4` Python package "
            "to read the ASAGI source. Import error: %s" % exc
        )

    ds = netCDF4.Dataset(str(path))
    try:
        for axis in ("x", "y", "z"):
            if axis not in ds.variables:
                raise SystemExit(
                    "build_pref_friction_sidecar: '%s' has no '%s' coordinate "
                    "variable; is it an ASAGI grid?" % (path, axis)
                )
        if "data" not in ds.variables:
            raise SystemExit(
                "build_pref_friction_sidecar: '%s' has no 'data' variable "
                "(expected an ASAGI compound with members %s, %s)."
                % (path, NC_A_MEMBER, NC_VW_MEMBER)
            )

        x = np.asarray(ds.variables["x"][:], dtype=np.float64)
        y = np.asarray(ds.variables["y"][:], dtype=np.float64)
        z = np.asarray(ds.variables["z"][:], dtype=np.float64)
        data = ds.variables["data"][:]

        members = getattr(data.dtype, "names", None)
        if not members:
            raise SystemExit(
                "build_pref_friction_sidecar: '%s' variable 'data' is not a "
                "compound type (dtype %s)." % (path, data.dtype)
            )
        for want in (NC_A_MEMBER, NC_VW_MEMBER):
            if want not in members:
                raise SystemExit(
                    "build_pref_friction_sidecar: compound 'data' is missing "
                    "member '%s'; has %s." % (want, list(members))
                )

        def member_xyz(name):
            # (z, y, x) -> (x, y, z), promoted to float64.
            a = np.asarray(data[name], dtype=np.float64)
            if np.ma.isMaskedArray(a):
                a = a.filled(np.nan)
            return np.ascontiguousarray(np.transpose(a, (2, 1, 0)))

        fields = {
            NC_A_MEMBER: member_xyz(NC_A_MEMBER),
            NC_VW_MEMBER: member_xyz(NC_VW_MEMBER),
        }
    finally:
        ds.close()

    return x, y, z, fields


# ---------------------------------------------------------------------------
# trilinear sampling — mirrors io/data_field_3d.cpp's INTERPOLATION, so the G4
# gate below tests what MFEM will actually read, not what numpy thinks the
# nearest node is.
#
# R-006: `locate()` CLAMPS at both ends, i.e. it reproduces `OOBPolicy::Clamp`,
# NOT the schema-v1 default `OOBPolicy::Abort`.  Callers MUST verify containment
# before sampling, or a query outside the hull will silently return an edge
# value.  `check_hypocentre` does exactly that, before it calls this.
# ---------------------------------------------------------------------------
def trilinear(x, y, z, field_xyz, px, py, pz):
    def locate(ax, p):
        if p <= ax[0]:
            return 0, 0.0
        if p >= ax[-1]:
            return len(ax) - 2, 1.0
        i = int(np.searchsorted(ax, p, side="right")) - 1
        i = min(max(i, 0), len(ax) - 2)
        t = (p - ax[i]) / (ax[i + 1] - ax[i])
        return i, t

    i, tx = locate(x, px)
    j, ty = locate(y, py)
    k, tz = locate(z, pz)
    c = field_xyz
    out = 0.0
    for di, wi in ((0, 1.0 - tx), (1, tx)):
        for dj, wj in ((0, 1.0 - ty), (1, ty)):
            for dk, wk in ((0, 1.0 - tz), (1, tz)):
                out += wi * wj * wk * c[i + di, j + dj, k + dk]
    return float(out)


def check_axes(x, y, z):
    for ax, lab in ((x, "x"), (y, "y"), (z, "z")):
        if ax.size < 2:
            raise SystemExit(
                "build_pref_friction_sidecar: axis '%s' has %d node(s); "
                "trilinear interpolation needs >= 2." % (lab, ax.size)
            )
        if not np.all(np.diff(ax) > 0.0):
            raise SystemExit(
                "build_pref_friction_sidecar: axis '%s' is not strictly "
                "increasing (DataField3D requires monotone axes)." % lab
            )


def check_hypocentre(x, y, z, fields, hypo, rtol):
    """G4 gate: the CASE2 anchors at the on-fault hypocentre.

    The tolerance is RELATIVE, not absolute.  The ASAGI members are float32, so
    a faithful round-trip of `a = 0.015` lands ~3e-10 off and `V_w = 0.05`
    ~7e-10 off; float32 eps at 0.05 is ~6e-9.  An absolute 1e-9 gate would
    therefore pass or fail on storage noise rather than on physics.  `rtol`
    defaults to 1e-6 -- four orders above float32 noise, four orders below the
    0.004 a-b plateau step this gate exists to detect.
    """
    px, py, pz = hypo
    inside = (x[0] <= px <= x[-1] and y[0] <= py <= y[-1] and z[0] <= pz <= z[-1])
    if not inside:
        raise SystemExit(
            "build_pref_friction_sidecar: hypocentre (%.4f, %.4f, %.4f) is "
            "OUTSIDE the sidecar hull x[%.0f,%.0f] y[%.0f,%.0f] z[%.0f,%.0f]; "
            "the CTM grid must cover the fault." %
            (px, py, pz, x[0], x[-1], y[0], y[-1], z[0], z[-1])
        )

    a_hypo = trilinear(x, y, z, fields[NC_A_MEMBER], px, py, pz)
    vw_hypo = trilinear(x, y, z, fields[NC_VW_MEMBER], px, py, pz)

    problems = []
    for got, expected, label in ((a_hypo, A_HYPO_EXPECTED, "a"),
                                 (vw_hypo, VW_HYPO_EXPECTED, "V_w")):
        tol = rtol * abs(expected)
        if abs(got - expected) > tol:
            problems.append("%s = %.12g (expected %.12g, |err| = %.3g > tol %.3g)"
                            % (label, got, expected, abs(got - expected), tol))
    if problems:
        raise SystemExit(
            "build_pref_friction_sidecar: G4 hypocentre anchor FAILED at "
            "(%.4f, %.4f, %.4f), rtol %g: %s.  The hypocentre should sit in the "
            "velocity-weakening plateau (T = 277.8 C).  Either the nc is not "
            "CASE2, or the hypocentre moved."
            % (px, py, pz, rtol, "; ".join(problems))
        )
    return a_hypo, vw_hypo


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Convert the SeisSol ASAGI thermal friction nc into an "
                    "MFEM data_projection_v1 HDF5 sidecar.")
    ap.add_argument("--in", dest="src", required=True, type=Path,
                    help="safs_friction_thermal_case2.nc (ASAGI, compound "
                         "data{rs_a, rs_srW})")
    ap.add_argument("--out", dest="dst", required=True, type=Path,
                    help="destination friction_safs.h5")
    ap.add_argument("--hypocenter", nargs=3, type=float, metavar=("X", "Y", "Z"),
                    default=list(HYPO_DEFAULT),
                    help="on-fault hypocentre for the G4 anchor gate "
                         "(default: the v3_4_1 PREFERRED 10-km on-fault point)")
    ap.add_argument("--hypo-rtol", type=float, default=1e-6,
                    help="RELATIVE tolerance for the G4 anchor gate (default "
                         "1e-6: far above float32 storage noise ~1e-7 rel, far "
                         "below the 0.004 a-b plateau step)")
    ap.add_argument("--mesh-tag", default="safv4_deep19km_sub500m_flattop_"
                                          "nwtrim_fixed_opt",
                    help="informational `mesh_tag` root attribute")
    ap.add_argument("--case", default="THERMAL_CASE2",
                    help="informational `case` root attribute")
    args = ap.parse_args(argv)

    if not args.src.is_file():
        raise SystemExit("build_pref_friction_sidecar: --in '%s' not found"
                         % args.src)

    x, y, z, fields = load_friction_nc(args.src)
    check_axes(x, y, z)

    a = fields[NC_A_MEMBER]
    vw = fields[NC_VW_MEMBER]
    print("  grid   : %d x %d x %d  (x, y, z)" % (x.size, y.size, z.size))
    print("  x      : [%.1f, %.1f]" % (x[0], x[-1]))
    print("  y      : [%.1f, %.1f]" % (y[0], y[-1]))
    print("  z      : [%.1f, %.1f]  (elevation, positive up)" % (z[0], z[-1]))
    print("  rs_a   : [%.6g, %.6g]" % (a.min(), a.max()))
    print("  rs_srW : [%.6g, %.6g] m/s" % (vw.min(), vw.max()))

    a_hypo, vw_hypo = check_hypocentre(x, y, z, fields, args.hypocenter,
                                       args.hypo_rtol)
    print("  G4 hypocentre (%.4f, %.4f, %.4f): a = %.12g, V_w = %.12g m/s  OK"
          % (args.hypocenter[0], args.hypocenter[1], args.hypocenter[2],
             a_hypo, vw_hypo))

    args.dst.parent.mkdir(parents=True, exist_ok=True)
    write_sidecar(
        args.dst, x, y, z,
        fields={NC_A_MEMBER: a, NC_VW_MEMBER: vw},
        attrs={
            "schema_version": "data_projection_v1",
            "crs": "EPSG:32611",
            "units": "m",
            "z_positive": "elevation",
            "source": "%s (SCEC CTM shinevar2024 temperature zoning) -> "
                      "build_pref_friction_sidecar.py" % args.src.name,
            "case": args.case,
            "mesh_tag": args.mesh_tag,
        },
        field_bounds={
            NC_A_MEMBER: (A_BOUNDS[0], A_BOUNDS[1], "1"),
            NC_VW_MEMBER: (VW_BOUNDS[0], VW_BOUNDS[1], "m/s"),
        },
    )
    print("wrote %s" % args.dst)
    print("  two fields (rs_a, rs_srW); no sign flip, no unit conversion.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
