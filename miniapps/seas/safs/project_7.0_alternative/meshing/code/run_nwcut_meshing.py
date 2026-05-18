#!/usr/bin/env python3
"""
run_nwcut_meshing.py — Mesh the NW-cut alternative ALT6 fault inside a
box whose bounding box is fitted to the cut trace.

For each cut STL in ../results/stl_nwcut/, this script:

  1. Reads the STL and computes its tight (x, y, z) bounding box.
  2. Invokes gmsh on safs_fault_box_nwcut.geo, injecting the bbox
     and STL filename via -setnumber / -setstring.
  3. Converts the resulting .msh to two .vtu files
     (<base>_bulk.vtu, <base>_fault.vtu) via the sibling
     msh_to_vtu.py.
  4. Reports min element size and per-tet quality stats.

Mesh .msh files land in ../results/msh/ and the derived .vtu files
in ../results/vtu/.

Usage:
    conda activate pythonenv
    python run_nwcut_meshing.py [--res 500 1000 2000]
                                [--lc-near 1500] [--lc-far 10000]
                                [--threads N] [--quiet]

If --res is omitted, all three resolutions are meshed.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent           # meshing/code/
PILLAR = HERE.parent                              # meshing/
PROJECT = PILLAR.parent                           # project_7.0_alternative/
DATA_CUT = PILLAR / "results" / "stl_nwcut"
RESULTS_MSH = PILLAR / "results" / "msh"
RESULTS_VTU = PILLAR / "results" / "vtu"
DATA_PROJECTION = PROJECT / "velocity" / "code"   # houses the flat data-projection scripts
GEO = HERE / "safs_fault_box_nwcut.geo"
MSH_TO_VTU = HERE / "msh_to_vtu.py"

DEFAULT_RES = (2000, 1000, 500)
BASELINE_CELL_COUNTS_JSON = HERE / ".baseline_cell_counts.json"
CELL_BUDGET_FACTOR = 2.0


def _read_bulk_cell_count(stats: str) -> int | None:
    """Pull the bulk tet count out of msh_to_vtu.py's stdout block.
    Returns None if not parseable."""
    import re
    m = re.search(r"bulk\s+cells\s*[:=]\s*([0-9]+)", stats)
    if m:
        return int(m.group(1))
    return None


def _check_cell_budget(res_m: int, base: Path, stats: str) -> None:
    """G-1 mesh-budget guard: warn (don't fail) if the new mesh exceeds
    CELL_BUDGET_FACTOR x the baseline cell count for this resolution.
    Baseline counts are stored in BASELINE_CELL_COUNTS_JSON
    (populated by hand the first time from existing pre-G-1 meshes).
    """
    import json
    if not BASELINE_CELL_COUNTS_JSON.is_file():
        print(f"  NOTE: no baseline cell-count JSON at "
              f"{BASELINE_CELL_COUNTS_JSON}; cannot enforce "
              f"{CELL_BUDGET_FACTOR}x budget cap.  Run once without "
              f"--gen-size-field to populate it (or write by hand).")
        return
    try:
        baselines = json.loads(BASELINE_CELL_COUNTS_JSON.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        print(f"  WARNING: could not read {BASELINE_CELL_COUNTS_JSON}: "
              f"{exc}; budget guard skipped.")
        return
    key = str(res_m)
    if key not in baselines:
        print(f"  NOTE: no baseline entry for {res_m} m in "
              f"{BASELINE_CELL_COUNTS_JSON.name}; budget guard skipped.")
        return
    new_count = _read_bulk_cell_count(stats)
    if new_count is None:
        print(f"  WARNING: could not parse bulk cell count from "
              f"msh_to_vtu output; budget guard skipped.")
        return
    base_count = int(baselines[key])
    ratio = new_count / max(1, base_count)
    if ratio > CELL_BUDGET_FACTOR:
        print(f"  WARNING: G-1 remesh exceeds {CELL_BUDGET_FACTOR}x "
              f"cell-count budget (got {new_count}, baseline "
              f"{base_count}, ratio {ratio:.2f}); consider lowering "
              f"--alpha or raising --lc-min")
    else:
        print(f"  cell budget: {new_count} / {base_count} "
              f"= {ratio:.2f}x (within {CELL_BUDGET_FACTOR}x cap)")


def find_gmsh() -> str:
    """Locate the `gmsh` binary. Prefer the one in the running Python's
    env (so `pythonenv` users don't need to activate the env in the
    subprocess), then fall back to PATH."""
    candidate = Path(sys.executable).parent / "gmsh"
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    on_path = shutil.which("gmsh")
    if on_path:
        return on_path
    raise RuntimeError(
        f"gmsh not found in {Path(sys.executable).parent} or on PATH; "
        f"activate pythonenv (or pass an env that has gmsh) before "
        f"running this script")


def stl_bbox(path: Path) -> np.ndarray:
    """Return (2, 3) array [[xmin, ymin, zmin], [xmax, ymax, zmax]]."""
    import meshio
    m = meshio.read(str(path))
    p = np.asarray(m.points, dtype=np.float64)
    if p.size == 0:
        raise ValueError(f"empty mesh: {path}")
    return np.vstack([p.min(axis=0), p.max(axis=0)])


def stl_for_res(res_m: int) -> Path:
    return (DATA_CUT
            / f"SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_{res_m}m"
              f"_clean_clip_nwcut.stl")


def clamp_stl_top(stl_in: Path, stl_out: Path, z_top: float) -> tuple[int, float]:
    """Snap any STL vertex with z > z_top down to z = z_top (in-memory
    read, write to stl_out).  Returns (n_clamped, max_orig_z) so the
    caller can log how aggressive the clamp was.  Used to pull the
    fault top edge strictly below the bulk box's top face so the gmsh
    PLC mesher does not see a fault-embedded-surface coplanar with the
    box top (which raises 'A segment and a facet intersect at point').
    """
    import meshio
    m = meshio.read(str(stl_in))
    pts = np.asarray(m.points, dtype=np.float64)
    mask = pts[:, 2] > z_top
    n_clamped = int(mask.sum())
    max_orig_z = float(pts[:, 2].max()) if pts.size else 0.0
    if n_clamped:
        pts[mask, 2] = z_top
        m.points = pts
        meshio.write(str(stl_out), m, file_format="stl")
    else:
        if stl_in != stl_out:
            stl_out.write_bytes(stl_in.read_bytes())
    return n_clamped, max_orig_z


def out_base(res_m: int, suffix: str = "") -> Path:
    """Stem path for the .msh / .size_field.pos outputs (lands in
    results/msh/).  The .vtu pair lands in results/vtu/ under the same
    stem (see vtu_base).  ``suffix`` should include any leading
    underscore (e.g. ``"_lcfar3000"``, ``"_zgraded"``); empty for the
    baseline."""
    return RESULTS_MSH / f"safs_fault_box_nwcut_{res_m}m{suffix}"


def vtu_base(res_m: int, suffix: str = "") -> Path:
    """Stem path for the bulk/fault .vtu pair (lands in results/vtu/)."""
    return RESULTS_VTU / f"safs_fault_box_nwcut_{res_m}m{suffix}"


def run_gmsh(stl_path: Path, bbox: np.ndarray, msh_path: Path,
             lc_near: float, lc_far: float,
             dist_inner: float, dist_outer: float,
             pad_x_lo: float, pad_x_hi: float,
             pad_y_lo: float, pad_y_hi: float,
             pad_bottom: float, pad_top: float,
             threads: int, verbose: bool,
             size_field_pos: Path | None = None,
             use_z_graded: bool = False,
             lc_floor: float = 100.0,
             do_optimize: bool | None = None,
             do_optimize_netgen: bool | None = None,
             optimize_threshold: float | None = None,
             smoothing_passes: int | None = None) -> float:
    """Call gmsh -3 with the parameterised geo. Returns wall-clock seconds."""
    gmsh_bin = find_gmsh()
    cmd = [
        gmsh_bin, "-3", str(GEO),
        "-setstring", "fault_stl", str(stl_path.resolve()),
        "-setnumber", "xmin_fault", f"{bbox[0, 0]}",
        "-setnumber", "xmax_fault", f"{bbox[1, 0]}",
        "-setnumber", "ymin_fault", f"{bbox[0, 1]}",
        "-setnumber", "ymax_fault", f"{bbox[1, 1]}",
        "-setnumber", "zmin_fault", f"{bbox[0, 2]}",
        "-setnumber", "zmax_fault", f"{bbox[1, 2]}",
        "-setnumber", "lc_near", f"{lc_near}",
        "-setnumber", "lc_far", f"{lc_far}",
        "-setnumber", "dist_inner", f"{dist_inner}",
        "-setnumber", "dist_outer", f"{dist_outer}",
        "-setnumber", "pad_x_lo", f"{pad_x_lo}",
        "-setnumber", "pad_x_hi", f"{pad_x_hi}",
        "-setnumber", "pad_y_lo", f"{pad_y_lo}",
        "-setnumber", "pad_y_hi", f"{pad_y_hi}",
        "-setnumber", "pad_bottom", f"{pad_bottom}",
        "-setnumber", "pad_top",   f"{pad_top}",
        "-setnumber", "lc_min",    f"{lc_floor}",
        "-nt", str(threads),
        # MFEM's mesh/mesh_readers.cpp:ReadGmshMesh is a Gmsh v2.2-only
        # parser (the function accepts version >= 2.2 but then runs the
        # v2.2 body unconditionally — there is no v4 branch).  Emitting
        # `msh4` here would produce artifacts unreadable by every MFEM
        # driver in this tree (seas_spatial_dyn_driver,
        # seas_project_velocity_to_mesh, seas_project_stress_to_mesh,
        # any future quasi-dynamic driver) with the misleading abort
        # "Gmsh file : vertices indices are not unique" at
        # mesh_readers.cpp:1628.  See
        # docs/DEBUG_msh4_mfem_incompat.md for the full investigation.
        "-format", "msh22",
        "-o", str(msh_path),
    ]
    if size_field_pos is not None:
        cmd.extend([
            "-setnumber", "USE_SIZE_FIELD", "1",
            "-setstring", "size_field_pos",
            str(Path(size_field_pos).resolve()),
        ])
    if use_z_graded:
        cmd.extend(["-setnumber", "USE_Z_GRADED", "1"])
    # PLAN_mesh_quality.md Phase 3: optimizer overrides.  Forward only
    # when a value was explicitly provided so the .geo's defaults
    # remain authoritative for the no-override case.
    if do_optimize is not None:
        cmd.extend(["-setnumber", "do_optimize",
                    "1" if do_optimize else "0"])
    if do_optimize_netgen is not None:
        cmd.extend(["-setnumber", "do_optimize_netgen",
                    "1" if do_optimize_netgen else "0"])
    if optimize_threshold is not None:
        cmd.extend(["-setnumber", "optimize_threshold",
                    f"{optimize_threshold}"])
    if smoothing_passes is not None:
        cmd.extend(["-setnumber", "smoothing_passes",
                    f"{smoothing_passes}"])
    if not verbose:
        cmd.append("-v")
        cmd.append("3")
    print(f"  $ {' '.join(cmd)}")
    t0 = time.time()
    r = subprocess.run(cmd, cwd=str(HERE), capture_output=not verbose,
                       text=True)
    elapsed = time.time() - t0
    if r.returncode != 0:
        if not verbose:
            sys.stderr.write(r.stdout or "")
            sys.stderr.write(r.stderr or "")
        raise RuntimeError(
            f"gmsh failed (rc={r.returncode}) for {stl_path.name}")
    return elapsed


def run_msh_to_vtu(msh_path: Path, base: Path, verbose: bool) -> str:
    """Call msh_to_vtu.py and return its stdout (which contains the
    edge-length / quality stats we need to report)."""
    cmd = [sys.executable, str(MSH_TO_VTU), str(msh_path), str(base)]
    print(f"  $ {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout or "")
        sys.stderr.write(r.stderr or "")
        raise RuntimeError(f"msh_to_vtu failed (rc={r.returncode})")
    if verbose:
        print(r.stdout)
    return r.stdout


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--res", nargs="+", type=int, default=None,
                    metavar="R",
                    help="resolution(s) in metres to mesh (subset of "
                         f"{DEFAULT_RES}); default: all three")
    ap.add_argument("--lc-near", type=float, default=1500.0,
                    help="near-fault target tet edge length [m] "
                         "(default: %(default)s)")
    ap.add_argument("--lc-far", type=float, default=10000.0,
                    help="far-field target tet edge length [m] "
                         "(default: %(default)s)")
    ap.add_argument("--dist-inner", type=float, default=3000.0,
                    help="distance from fault below which size = lc_near "
                         "(default: %(default)s)")
    ap.add_argument("--dist-outer", type=float, default=40000.0,
                    help="distance from fault above which size = lc_far "
                         "(default: %(default)s)")
    # PLAN_mesh_quality.md Phase 3 -- gmsh optimizer knobs.  Defaults
    # live in safs_fault_box_nwcut.geo (do_optimize=1, do_optimize_netgen=1,
    # optimize_threshold=0.3, smoothing_passes=5); these flags forward
    # overrides via -setnumber.  Forwarding is conditional so that the
    # .geo's If(!Exists()) guards remain authoritative when nothing is
    # passed on the CLI.
    ap.add_argument("--no-optimize", action="store_true",
                    dest="no_optimize",
                    help="disable gmsh's Mesh.Optimize Laplacian-style "
                         "pass (PLAN_mesh_quality.md Phase 3; on by "
                         "default in the .geo).")
    ap.add_argument("--optimize-netgen", action="store_true",
                    dest="optimize_netgen",
                    help="ENABLE gmsh's Mesh.OptimizeNetgen tet-quality "
                         "optimizer (OFF by default in the .geo; see "
                         "PLAN_mesh_quality.md Phase 3 Risk Assessment "
                         "-- this optimizer empirically crashes with "
                         "SIGBUS on the SAFS NW-cut embedded-fault "
                         "geometry, so it is opt-in only).")
    ap.add_argument("--optimize-threshold", type=float, default=None,
                    dest="optimize_threshold",
                    help="gmsh Mesh.OptimizeThreshold (gamma metric, "
                         "0..1).  Tets with quality below this value "
                         "are attacked by the optimizer.  Default in "
                         "the .geo: 0.3 (approx. eta ~ 0.5).")
    ap.add_argument("--smoothing-passes", type=int, default=None,
                    dest="smoothing_passes",
                    help="gmsh Mesh.Smoothing pass count.  Default in "
                         "the .geo: 5.")
    ap.add_argument("--lc-floor", type=float, default=100.0,
                    dest="lc_floor",
                    help="hard lower bound on the local mesh size [m] "
                         "applied via Field[Max] in the .geo (default: "
                         "%(default)s; see PLAN_mesh_quality.md Phase 2). "
                         "Counteracts gmsh's embedded-discrete-surface "
                         "code path which can produce sub-100 m bulk "
                         "edges near the fault top despite "
                         "Mesh.CharacteristicLengthMin.  Forwarded to "
                         "the .geo as the variable `lc_min`.  Must be "
                         "< --lc-near or the driver aborts.  Set 0 to "
                         "disable the floor.  Note: distinct from "
                         "the existing --lc-min flag, which clamps the "
                         "build_pos size-field generator (used only "
                         "with --gen-size-field).")
    ap.add_argument("--pad-x", type=float, default=None,
                    help="lateral padding in +/-x [m] "
                         "(default: 50 km, or computed from --double-domain)")
    ap.add_argument("--pad-y", type=float, default=None,
                    help="lateral padding in +/-y [m] "
                         "(default: 50 km, or computed from --double-domain)")
    ap.add_argument("--pad-bottom", type=float, default=None,
                    help="padding below deepest fault vertex [m] "
                         "(default: 25 km, or 50 km with --double-domain)")
    ap.add_argument("--pad-top", type=float, default=None,
                    help="padding above shallowest (clamped) fault vertex "
                         "[m] (default: equal to --fault-top-clamp so the "
                         "box top = 0 exactly when the unclamped STL has "
                         "fault_zmax = 0).")
    ap.add_argument("--fault-top-clamp", type=float, default=100.0,
                    dest="fault_top_clamp",
                    help="[m] snap any fault vertex with z > -<this> down "
                         "to z = -<this> before passing to gmsh "
                         "(default: %(default)s).  Provides gmsh enough "
                         "vertical headroom near the free surface to "
                         "avoid slivers from the fault-corridor / "
                         "free-surface pinch (PLAN_mesh_quality.md "
                         "Phase 1).  Required because the cleaned STL "
                         "has fault_zmax = 0, which would be coplanar "
                         "with the box top face and break gmsh's "
                         "discrete-surface embedding.  The clamped STL "
                         "is written next to the original with a "
                         "'_zclamp<N>m.stl' suffix.  Physics ceiling: "
                         "250 m -- above that, SAFS co-seismic surface "
                         "rupture is suppressed (a warning is printed "
                         "but the build proceeds).  Set 0 to disable "
                         "clamping (only safe if you also set "
                         "--pad-top > 0).")
    ap.add_argument("--fault-top-clamp-policy",
                    choices=("fixed", "auto"), default="fixed",
                    dest="fault_top_clamp_policy",
                    help="'fixed' (default) uses --fault-top-clamp "
                         "verbatim. 'auto' computes "
                         "clamp = clip(surface_lc, 100, 250) where "
                         "surface_lc = z_lc_top (= 500 m, the .geo "
                         "default) in --z-graded mode, else --lc-near. "
                         "Use 'auto' to pick the smallest viable clamp "
                         "for the local surface refinement.")
    ap.add_argument("--velocity-bbox", type=float, nargs=4, default=None,
                    metavar=("LON_LO", "LAT_LO", "LON_HI", "LAT_HI"),
                    help="auto-fit the mesh into the inscribed UTM 11N "
                         "AABB of this geographic bbox.  Mutually "
                         "exclusive with --double-domain and per-axis "
                         "--pad-* flags (which take precedence if both "
                         "are passed).  pad_bottom defaults to "
                         "min(50 km, fault_zmin - velocity_zmin); "
                         "velocity_zmin is read from --velocity-z-bottom "
                         "(default -70 km, matching CVM-H archive depth).")
    ap.add_argument("--velocity-z-bottom", type=float, default=-70000.0,
                    help="elevation [m] of the velocity dataset's "
                         "lowest slice (default: %(default)s)")
    ap.add_argument("--sidecar", type=Path, default=None,
                    help="auto-fit the mesh inside the data bbox of an "
                         "existing schema-v1 HDF5 sidecar (reads "
                         "/grid/x, /grid/y, /grid/z directly so there "
                         "is no quantization gap between the mesh and "
                         "the sidecar).  Mutually exclusive with "
                         "--velocity-bbox; if both are given, --sidecar "
                         "takes precedence.  pad_top is still taken from "
                         "--pad-top (default: --fault-top-clamp, so "
                         "mesh_zmax = 0); the sidecar must have been "
                         "built with --extend-z-top >= pad_top or the "
                         "resulting mesh will fail strict interpolation "
                         "containment at the top.")
    ap.add_argument("--double-domain", action="store_true",
                    help="double the total xy box extent (fault stays "
                         "centred -> PAD_X = fault_x/2 + 100 km, "
                         "PAD_Y = fault_y/2 + 100 km), and double the "
                         "downward z padding (PAD_BOTTOM = 50 km).  "
                         "Free-surface PAD_TOP unchanged.")
    ap.add_argument("--threads", type=int,
                    default=max(1, (os.cpu_count() or 1) // 2),
                    help="gmsh thread count (default: half of cpu_count)")
    ap.add_argument("--quiet", action="store_true",
                    help="suppress gmsh verbose output")
    # ----- G-1: sidecar-driven background size field (PLAN.md Phase 2) ----
    ap.add_argument("--size-field-pos", type=Path, default=None,
                    dest="size_field_pos",
                    help="path to a pre-built gmsh PostView .pos for the "
                         "background size field (alternative to "
                         "--gen-size-field).  Activates USE_SIZE_FIELD=1 "
                         "in the .geo template.")
    ap.add_argument("--gen-size-field", action="store_true",
                    dest="gen_size_field",
                    help="auto-generate the size field from --sidecar via "
                         "build_size_field.build_pos before each gmsh "
                         "run.  Mutually exclusive with --size-field-pos.")
    ap.add_argument("--lc-min", type=float, default=500.0, dest="lc_min",
                    help="absolute floor on tet edge length applied by "
                         "--gen-size-field (default: %(default)s)")
    ap.add_argument("--alpha", type=float, default=8.0, dest="alpha",
                    help="size-field gradient sensitivity (default: "
                         "%(default)s); --alpha 6/4 if 2x cell-budget "
                         "guard fires")
    ap.add_argument("--smooth-sigma", type=float, default=1.0,
                    dest="smooth_sigma",
                    help="Gaussian smoothing sigma (voxel units) on the "
                         "sidecar gradient (default: %(default)s)")
    ap.add_argument("--voxel-stride", type=int, default=2,
                    dest="voxel_stride",
                    help="downsample factor for the PostView grid "
                         "(default: %(default)s)")
    # ----- Variant knobs (no per-variant scripts; everything is
    # parameterised on the shared safs_fault_box_nwcut.geo template) ----
    ap.add_argument("--z-graded", action="store_true", dest="z_graded",
                    help="set USE_Z_GRADED=1 in the .geo template so the "
                         "background size field is capped by a "
                         "depth-dependent step function (see the .geo "
                         "header for the four-slab schedule).  Combined "
                         "with the fault-distance field via Field[Min].")
    ap.add_argument("--suffix", type=str, default="",
                    help="suffix appended to the output mesh stem so "
                         "different variants don't overwrite each "
                         "other (e.g. '_lcfar3000', '_zgraded').  Must "
                         "start with '_' if non-empty.  Default: empty "
                         "(baseline).")
    args = ap.parse_args()
    if args.suffix and not args.suffix.startswith("_"):
        print(f"ERROR: --suffix must start with '_' or be empty; "
              f"got {args.suffix!r}", file=sys.stderr)
        return 1

    # PLAN_mesh_quality.md Phase 1: --fault-top-clamp policy resolution.
    # The .geo default for z_lc_top (zgraded basin-surface lc) is 500 m;
    # the constant here must stay in sync with safs_fault_box_nwcut.geo.
    GEO_Z_LC_TOP_DEFAULT = 500.0
    PHYSICS_CEILING_M = 250.0
    if args.fault_top_clamp_policy == "auto":
        surface_lc = (GEO_Z_LC_TOP_DEFAULT if args.z_graded
                      else float(args.lc_near))
        args.fault_top_clamp = max(100.0, min(PHYSICS_CEILING_M,
                                              float(surface_lc)))
        print(f"--fault-top-clamp-policy=auto: surface_lc={surface_lc} m "
              f"-> fault_top_clamp={args.fault_top_clamp} m")
    if args.fault_top_clamp > PHYSICS_CEILING_M:
        print(f"WARNING: --fault-top-clamp={args.fault_top_clamp} m exceeds "
              f"recommended {PHYSICS_CEILING_M:.0f} m physics ceiling -- "
              f"SAFS co-seismic surface rupture will be suppressed.",
              file=sys.stderr)

    # PLAN_mesh_quality.md Phase 2: refuse lc_floor >= lc_near (would
    # erase the fault-corridor refinement entirely).
    if args.lc_floor > 0.0 and args.lc_floor >= args.lc_near:
        print(f"ERROR: --lc-floor={args.lc_floor} >= "
              f"--lc-near={args.lc_near}: this would erase the "
              f"fault-corridor refinement.  Pick lc_floor < lc_near.",
              file=sys.stderr)
        return 1

    if not GEO.is_file():
        print(f"missing template: {GEO}", file=sys.stderr)
        return 1
    RESULTS_MSH.mkdir(parents=True, exist_ok=True)
    RESULTS_VTU.mkdir(parents=True, exist_ok=True)
    if not MSH_TO_VTU.is_file():
        print(f"missing converter: {MSH_TO_VTU}", file=sys.stderr)
        return 1
    if args.size_field_pos is not None and args.gen_size_field:
        print("ERROR: --size-field-pos and --gen-size-field are "
              "mutually exclusive", file=sys.stderr)
        return 1
    if args.gen_size_field and args.sidecar is None:
        print("ERROR: --gen-size-field requires --sidecar", file=sys.stderr)
        return 1

    # If --sidecar is supplied, read its grid bbox directly.  This is
    # the cleanest path: zero quantization gap between the mesh and the
    # sidecar (the previous --velocity-bbox approach computes the
    # inscribed AABB at runtime and can be off by up to one grid cell
    # relative to the sidecar's grid_dx-rounded axes).
    velocity_utm_aabb = None
    if args.sidecar is not None:
        try:
            import h5py
        except ImportError:
            print("--sidecar requires h5py; install via "
                  "`pip install h5py`", file=sys.stderr)
            return 1
        sc = Path(args.sidecar)
        if not sc.is_file():
            print(f"--sidecar: file not found: {sc}", file=sys.stderr)
            return 1
        with h5py.File(sc, "r") as h5:
            def _attr_str(key):
                val = h5.attrs.get(key)
                if val is None:
                    return None
                return val.decode() if isinstance(val, bytes) else val
            schema = _attr_str("schema_version")
            if schema != "data_projection_v1":
                print(f"--sidecar: unsupported schema_version='{schema}', "
                      f"expected 'data_projection_v1'", file=sys.stderr)
                return 1
            for k, expected in (("crs",        "EPSG:32611"),
                                ("units",      "m"),
                                ("z_positive", "elevation")):
                got = _attr_str(k)
                if got != expected:
                    print(f"--sidecar: attribute {k}='{got}' does not "
                          f"match canonical value '{expected}'",
                          file=sys.stderr)
                    return 1
            gx = np.asarray(h5["grid/x"][...], dtype=np.float64)
            gy = np.asarray(h5["grid/y"][...], dtype=np.float64)
            gz = np.asarray(h5["grid/z"][...], dtype=np.float64)
        v_xmin = float(gx[0]);  v_xmax = float(gx[-1])
        v_ymin = float(gy[0]);  v_ymax = float(gy[-1])
        v_zmin = float(gz[0]);  v_zmax = float(gz[-1])
        velocity_utm_aabb = (v_xmin, v_xmax, v_ymin, v_ymax,
                             v_zmin, v_zmax)
        print(f"\n=== sidecar bbox ({sc.name}) ===")
        print(f"  x = [{v_xmin:.0f}, {v_xmax:.0f}]   "
              f"({(v_xmax-v_xmin)/1000:.2f} km)")
        print(f"  y = [{v_ymin:.0f}, {v_ymax:.0f}]   "
              f"({(v_ymax-v_ymin)/1000:.2f} km)")
        print(f"  z = [{v_zmin:.0f}, {v_zmax:.0f}]   "
              f"({(v_zmax-v_zmin)/1000:.2f} km, top includes "
              f"any --extend-z-top from sidecar build)")

    # If --velocity-bbox is supplied (and --sidecar was not), compute
    # the inscribed UTM AABB up front so we can derive per-resolution
    # asymmetric pads from it.
    if velocity_utm_aabb is None and args.velocity_bbox is not None:
        lon_lo, lat_lo, lon_hi, lat_hi = args.velocity_bbox
        if not (lon_lo < lon_hi and lat_lo < lat_hi):
            print(f"--velocity-bbox: require LON_LO < LON_HI and "
                  f"LAT_LO < LAT_HI; got {args.velocity_bbox}",
                  file=sys.stderr)
            return 1
        # Lazy import so the rest of the script doesn't pull pyproj.
        sys.path.insert(0, str(DATA_PROJECTION))
        from crs import geographic_to_utm11n
        n_sample = 401
        lons = np.linspace(lon_lo, lon_hi, n_sample)
        lats = np.linspace(lat_lo, lat_hi, n_sample)
        LON, LAT = np.meshgrid(lons, lats, indexing="ij")
        X, Y = geographic_to_utm11n(LON, LAT)
        v_xmin = float(X[0, :].max())     # right-most along western lon
        v_xmax = float(X[-1, :].min())    # left-most  along eastern lon
        v_ymin = float(Y[:, 0].max())     # north-most along southern lat
        v_ymax = float(Y[:, -1].min())    # south-most along northern lat
        v_zmin = float(args.velocity_z_bottom)
        # No constraint on the velocity z_top in --velocity-bbox mode
        # (the user has not pointed us at a built sidecar yet); use
        # +infinity so the per-resolution pad_top check below trivially
        # passes.
        v_zmax = float("inf")
        velocity_utm_aabb = (v_xmin, v_xmax, v_ymin, v_ymax,
                             v_zmin, v_zmax)
        print("\n=== velocity geographic bbox ===")
        print(f"  lon ∈ [{lon_lo:.4f}, {lon_hi:.4f}], "
              f"lat ∈ [{lat_lo:.4f}, {lat_hi:.4f}]")
        print(f"  inscribed UTM 11N AABB: x=[{v_xmin:.0f}, {v_xmax:.0f}], "
              f"y=[{v_ymin:.0f}, {v_ymax:.0f}], z_lower={v_zmin:.0f}")

    res_list = args.res if args.res else list(DEFAULT_RES)
    summary = []
    rc_overall = 0
    for res in res_list:
        stl_src = stl_for_res(res)
        if not stl_src.is_file():
            print(f"\nSKIP {res} m: STL not found at {stl_src}", file=sys.stderr)
            rc_overall = 1
            continue
        base = out_base(res, args.suffix)
        msh = base.with_suffix(".msh")
        # Clamp the fault top a hair below z = 0 so the discrete fault
        # surface stays strictly interior to the bulk box (which has its
        # top face at z = 0).  Without this, gmsh's PLC mesher reports
        # 'A segment and a facet intersect at point' on the coplanar
        # top edge.
        if args.fault_top_clamp > 0.0:
            z_top_clamp = -float(args.fault_top_clamp)
            clamp_tag = f"_zclamp{int(round(args.fault_top_clamp))}m"
            stl = stl_src.with_name(stl_src.stem + clamp_tag + ".stl")
            n_clamped, max_orig_z = clamp_stl_top(stl_src, stl, z_top_clamp)
            print(f"\n=== {stl.name} ===")
            print(f"  fault-top clamp: snapped {n_clamped} vertices with "
                  f"z > {z_top_clamp:.3f} m down to z = {z_top_clamp:.3f} m "
                  f"(max original z = {max_orig_z:.4f} m)")
        else:
            stl = stl_src
            print(f"\n=== {stl.name} ===")
        bbox = stl_bbox(stl)
        print(f"  cut STL bbox:")
        print(f"    x ∈ [{bbox[0,0]:.2f}, {bbox[1,0]:.2f}]  "
              f"({(bbox[1,0]-bbox[0,0])/1000:.2f} km)")
        print(f"    y ∈ [{bbox[0,1]:.2f}, {bbox[1,1]:.2f}]  "
              f"({(bbox[1,1]-bbox[0,1])/1000:.2f} km)")
        print(f"    z ∈ [{bbox[0,2]:.2f}, {bbox[1,2]:.2f}]  "
              f"({(bbox[1,2]-bbox[0,2])/1000:.2f} km)")
        # Resolve pad values per resolution.  Order of precedence:
        #   1. --velocity-bbox: per-side pads sized to fit the velocity
        #      inscribed UTM AABB exactly (asymmetric).
        #   2. --double-domain: symmetric pads sized so domain doubles.
        #   3. Defaults: 50/50/25 km.
        # In all branches, --pad-x / --pad-y / --pad-bottom / --pad-top
        # explicit flags override the corresponding axes.
        fault_dx = float(bbox[1, 0] - bbox[0, 0])
        fault_dy = float(bbox[1, 1] - bbox[0, 1])
        fault_xmin = float(bbox[0, 0]); fault_xmax = float(bbox[1, 0])
        fault_ymin = float(bbox[0, 1]); fault_ymax = float(bbox[1, 1])
        fault_zmin = float(bbox[0, 2]); fault_zmax = float(bbox[1, 2])
        if velocity_utm_aabb is not None:
            (v_xmin, v_xmax, v_ymin, v_ymax,
             v_zmin, v_zmax) = velocity_utm_aabb
            default_pad_x_lo = fault_xmin - v_xmin
            default_pad_x_hi = v_xmax - fault_xmax
            default_pad_y_lo = fault_ymin - v_ymin
            default_pad_y_hi = v_ymax - fault_ymax
            default_pad_bottom = min(50_000.0, fault_zmin - v_zmin)
            src = "--sidecar" if args.sidecar else "--velocity-bbox"
            for label, val in (("pad_x_lo", default_pad_x_lo),
                                ("pad_x_hi", default_pad_x_hi),
                                ("pad_y_lo", default_pad_y_lo),
                                ("pad_y_hi", default_pad_y_hi),
                                ("pad_bottom", default_pad_bottom)):
                if val < 0.0:
                    print(f"  ERROR: {src} yields negative "
                          f"{label} = {val:.0f} m for {res} m fault — "
                          f"fault extends beyond velocity bbox on that "
                          f"side.  Re-cut the fault or widen velocity "
                          f"bbox.", file=sys.stderr)
                    rc_overall = 1
                    raise SystemExit(1)
            # pad_top headroom check (R-004): mesh top = fault_zmax +
            # args.pad_top must lie within the sidecar's z_top.  In
            # --velocity-bbox mode v_zmax is +inf, so this check is a
            # no-op there.
            head_room = v_zmax - fault_zmax
            effective_pad_top = (args.pad_top if args.pad_top is not None
                                  else args.fault_top_clamp)
            if effective_pad_top > head_room + 1.0e-6:
                print(f"  ERROR: pad_top={effective_pad_top:.1f} m exceeds "
                      f"available {src} headroom {head_room:.1f} m "
                      f"(z_top={v_zmax:.1f}, fault z_top="
                      f"{fault_zmax:.1f}).  Re-build the sidecar with "
                      f"--extend-z-top >= {effective_pad_top:.0f} or pass "
                      f"--pad-top <= {head_room:.0f}.",
                      file=sys.stderr)
                rc_overall = 1
                raise SystemExit(1)
        elif args.double_domain:
            default_pad_x_lo = default_pad_x_hi = fault_dx / 2 + 100_000.0
            default_pad_y_lo = default_pad_y_hi = fault_dy / 2 + 100_000.0
            default_pad_bottom = 50_000.0
        else:
            default_pad_x_lo = default_pad_x_hi = 50_000.0
            default_pad_y_lo = default_pad_y_hi = 50_000.0
            default_pad_bottom = 25_000.0

        # Explicit --pad-x / --pad-y override BOTH sides symmetrically.
        if args.pad_x is not None:
            default_pad_x_lo = default_pad_x_hi = args.pad_x
        if args.pad_y is not None:
            default_pad_y_lo = default_pad_y_hi = args.pad_y
        pad_x_lo = default_pad_x_lo
        pad_x_hi = default_pad_x_hi
        pad_y_lo = default_pad_y_lo
        pad_y_hi = default_pad_y_hi
        pad_bottom = (args.pad_bottom if args.pad_bottom is not None
                      else default_pad_bottom)
        # Default pad_top = fault_top_clamp so the box top lands at z = 0
        # exactly (clamped fault zmax = -fault_top_clamp, plus pad_top
        # = fault_top_clamp, gives mesh_zmax = 0).
        pad_top = (args.pad_top if args.pad_top is not None
                   else args.fault_top_clamp)

        domain_dx = fault_dx + pad_x_lo + pad_x_hi
        domain_dy = fault_dy + pad_y_lo + pad_y_hi
        domain_dz = (fault_zmax - fault_zmin) + pad_bottom + pad_top
        mesh_xmin = fault_xmin - pad_x_lo
        mesh_xmax = fault_xmax + pad_x_hi
        mesh_ymin = fault_ymin - pad_y_lo
        mesh_ymax = fault_ymax + pad_y_hi
        mesh_zmin = fault_zmin - pad_bottom
        mesh_zmax = fault_zmax + pad_top
        print(f"  pads: x_lo={pad_x_lo/1000:.2f} km, "
              f"x_hi={pad_x_hi/1000:.2f} km, "
              f"y_lo={pad_y_lo/1000:.2f} km, y_hi={pad_y_hi/1000:.2f} km, "
              f"bottom={pad_bottom/1000:.2f} km, "
              f"top={pad_top/1000:.3f} km")
        print(f"  mesh bbox: x=[{mesh_xmin:.0f}, {mesh_xmax:.0f}], "
              f"y=[{mesh_ymin:.0f}, {mesh_ymax:.0f}], "
              f"z=[{mesh_zmin:.0f}, {mesh_zmax:.0f}]")
        print(f"  domain: {domain_dx/1000:.1f} x {domain_dy/1000:.1f} x "
              f"{domain_dz/1000:.1f} km")
        # Per-resolution size-field generation (G-1).  If --gen-size-field,
        # build the .pos right next to the .msh so it's reproducible.
        size_field_pos: Path | None = None
        if args.gen_size_field:
            try:
                sys.path.insert(0, str(DATA_PROJECTION))
                from build_size_field import build_pos
            except ImportError as exc:
                print(f"  ERROR: --gen-size-field needs "
                      f"build_size_field: {exc}",
                      file=sys.stderr)
                rc_overall = 1
                continue
            size_field_pos = base.with_suffix(".size_field.pos")
            try:
                meta = build_pos(
                    sidecar=args.sidecar,
                    out_pos=size_field_pos,
                    field="Vs",
                    lc_near=args.lc_near,
                    lc_far=args.lc_far,
                    lc_min=args.lc_min,
                    alpha=args.alpha,
                    smooth_sigma=args.smooth_sigma,
                    voxel_stride=args.voxel_stride,
                    verbose=not args.quiet)
                print(f"  size field meta: {meta}")
            except (FileNotFoundError, KeyError, ValueError) as exc:
                print(f"  ERROR: build_pos failed: {exc}",
                      file=sys.stderr)
                rc_overall = 1
                continue
        elif args.size_field_pos is not None:
            size_field_pos = args.size_field_pos

        try:
            t_gmsh = run_gmsh(stl, bbox, msh,
                              args.lc_near, args.lc_far,
                              args.dist_inner, args.dist_outer,
                              pad_x_lo, pad_x_hi,
                              pad_y_lo, pad_y_hi,
                              pad_bottom, pad_top,
                              args.threads, verbose=not args.quiet,
                              size_field_pos=size_field_pos,
                              use_z_graded=args.z_graded,
                              lc_floor=args.lc_floor,
                              do_optimize=(False if args.no_optimize
                                           else None),
                              do_optimize_netgen=(True if
                                                  args.optimize_netgen
                                                  else None),
                              optimize_threshold=args.optimize_threshold,
                              smoothing_passes=args.smoothing_passes)
        except RuntimeError as exc:
            print(str(exc), file=sys.stderr)
            rc_overall = 1
            continue
        msh_mb = msh.stat().st_size / 1e6
        print(f"  gmsh wall-clock: {t_gmsh:.1f} s; "
              f"{msh.name} = {msh_mb:.2f} MB")
        try:
            stats = run_msh_to_vtu(msh, vtu_base(res, args.suffix),
                                   verbose=False)
        except RuntimeError as exc:
            print(str(exc), file=sys.stderr)
            rc_overall = 1
            continue
        print(stats)
        summary.append((res, msh_mb, t_gmsh, stats))

        # G-1 mesh-budget guard (PLAN.md Phase 2 R-2-G-T3).
        if size_field_pos is not None:
            _check_cell_budget(res, base, stats)

    print("\n" + "=" * 70)
    print("Summary (NW-cut meshing)")
    print("=" * 70)
    for res, mb, t_g, stats in summary:
        print(f"\n[{res} m]   {out_base(res, args.suffix).name}.msh  "
              f"({mb:.2f} MB, gmsh {t_g:.1f} s)")
        for line in stats.splitlines():
            stripped = line.strip()
            if any(k in stripped for k in
                   ("bulk  cells", "fault cells",
                    "bulk edge length", "fault edge length",
                    "tet quality", "q<0.1")):
                print("    " + stripped)
    return rc_overall


if __name__ == "__main__":
    sys.exit(main())
