#!/usr/bin/env python3
"""
TPV102 v9.2.0 plan §6.4 / §15.5.2 (rev-3b) — offline VTU outlier detector.

Post-REVIEW R-V92-C02/C05: implemented as a FREE local probe ordered
before any Frontera SU expenditure.

What it does
------------
Loads the fault-surface PVD produced by `seas_tpv102_driver --paraview`
at job 7668434 (or any equivalent TPV102 run) and, for each timestep:
- Reads per-triangle scalar fields (`normal_stress`, `traction_strike`,
  `traction_dip`, `state_variable`, `slip_rate_strike`, `slip_rate_dip`,
  `slip_strike`, `slip_dip`).
- Computes the median + MAD (median absolute deviation) of each field.
- Flags a triangle as "outlier" if ANY field deviates by more than
  k * MAD from the median (default k = 5).
- Records the (x, z) centroid of each outlier triangle.

Classifier
----------
- **H-V92-P (shared-fault DOF divergence)**: outlier count stable over
  cycles; outliers coincide with a known partition-boundary subset
  (requires §5.2 rank-id VTU patch to be applied — optional).
- **H-V92-G (bulk σ_yy amplification)**: outlier count grows
  monotonically with time.
- **H-V92-Q (friction solver per-QP instability)**: outliers scatter
  randomly regardless of partition-boundary membership.

Usage
-----
    python3 fault_vtu_outlier_detector.py \
        --pvd  tpv102/results_200m_p1_12.0s_400r_v91_job7668434/tpv102.pvd \
        --output-root r_v92_outlier_detection/ \
        --mad-k 5

Data prerequisites
------------------
The raw `FaultSurface/fault_surface_*.vtu` and the top-level `.pvd`
must be pulled from Frontera (not stored locally).  Approximate
`scp` commands:

    scp frontera:/scratch2/10024/zhaochun/seas-project/seas-mfem/\
    miniapps/seas/tpv102/results_200m_p1_12.0s_400r_v91_job7668434/\
    *.pvd ./

(+ all the VTU shards referenced by the PVD.)

Per feedback_no_local_reproducer.md this is READ-ONLY; no simulation
is run.  Safe for local analysis.

Requires: numpy, vtk (conda: `conda install -c conda-forge vtk`).
"""

import argparse
import json
import os
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

try:
    import numpy as np
except ImportError as e:
    print("ERROR: this script requires `numpy`.  Install via:")
    print("  conda install -c conda-forge numpy")
    print(f"Import error: {e}")
    sys.exit(2)

# --- Lightweight ASCII VTU reader (no vtk module dependency) ---------------
# MFEM's ParaView writer produces format="ascii" DataArrays for the TPV102
# runs, so we can parse with pure Python XML + float conversion.


# ---------------------------------------------------------------------------
FIELDS_OF_INTEREST = [
    "normal_stress",
    "traction_strike",
    "traction_dip",
    "state_variable",
    "slip_rate_strike",
    "slip_rate_dip",
    "slip_strike",
    "slip_dip",
    # Coordinates + spatially-varying parameters — READ but NOT used in
    # the outlier classifier (they vary by design across the fault).
    "fault_x2",
    "fault_x3",
    "param_a",
    "param_Dc",
    # Optional — emitted by §5.2 patch when SEAS_DIAG_VTU_RANK is set.
    "mpi_rank",
    "is_shared_face",
]

# Fields that CONTRIBUTE to the MAD-outlier classifier.  Coordinate /
# parameter fields are spatially non-uniform by design and would be
# flagged trivially, so they are excluded.
CLASSIFIER_FIELDS = {
    "normal_stress",
    "traction_strike",
    "traction_dip",
    "state_variable",
    "slip_rate_strike",
    "slip_rate_dip",
    "slip_strike",
    "slip_dip",
}


def parse_pvd(pvd_path):
    """Return list of (time, vtu_path) tuples from a ParaView .pvd file."""
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    out = []
    for ds in root.iter("DataSet"):
        t = float(ds.attrib.get("timestep", 0.0))
        f = ds.attrib.get("file", "")
        out.append((t, f))
    return out


def _parse_ascii_array(text):
    """Parse whitespace-separated float values from an ASCII DataArray body."""
    return np.fromstring(text, sep=" ", dtype=np.float64)


def _parse_ascii_int_array(text):
    return np.fromstring(text, sep=" ", dtype=np.int64)


def load_vtu_fields(vtu_path):
    """Read triangle centroids + scalar fields from an ASCII .vtu file.

    Returns (centroids, fields) where centroids is (n_cells, 3) and
    fields is {name: numpy_array[n_cells]}.

    Assumes:
      - format="ascii" DataArrays,
      - Float64 Points (3 components),
      - Int64 Cells connectivity + offsets,
      - Float64 CellData scalar fields.
    """
    tree = ET.parse(str(vtu_path))
    root = tree.getroot()
    piece = next(root.iter("Piece"))
    n_points = int(piece.attrib["NumberOfPoints"])
    n_cells  = int(piece.attrib["NumberOfCells"])

    # Empty-rank shard: no fault faces on this rank.  Return empty arrays.
    if n_cells == 0 or n_points == 0:
        return np.zeros((0, 3)), {}

    # --- Points ---
    pts_node = piece.find("Points/DataArray")
    if pts_node is None or pts_node.attrib.get("format") != "ascii":
        raise RuntimeError(f"{vtu_path}: non-ASCII or missing Points")
    pts = _parse_ascii_array(pts_node.text).reshape(n_points, 3)

    # --- Cells (connectivity + offsets) ---
    conn_arr = offsets_arr = None
    for da in piece.findall("Cells/DataArray"):
        name = da.attrib.get("Name", "")
        if name == "connectivity":
            conn_arr = _parse_ascii_int_array(da.text)
        elif name == "offsets":
            offsets_arr = _parse_ascii_int_array(da.text)
    if conn_arr is None or offsets_arr is None:
        raise RuntimeError(f"{vtu_path}: missing connectivity/offsets")

    # Triangle centroids from connectivity.
    centroids = np.zeros((n_cells, 3))
    prev = 0
    for ci in range(n_cells):
        end = int(offsets_arr[ci])
        pids = conn_arr[prev:end]
        prev = end
        centroids[ci] = pts[pids].mean(axis=0)

    # --- CellData scalar fields ---
    fields = {}
    for da in piece.findall("CellData/DataArray"):
        name = da.attrib.get("Name", "")
        if name not in FIELDS_OF_INTEREST:
            continue
        if da.attrib.get("format") != "ascii":
            continue
        vals = _parse_ascii_array(da.text)
        if vals.size == n_cells:
            fields[name] = vals

    return centroids, fields


def classify_outliers(centroids, fields, mad_k):
    """Return boolean mask of outlier triangles + per-field diagnostics."""
    n_cells = centroids.shape[0]
    outlier_mask = np.zeros(n_cells, dtype=bool)
    field_stats = {}
    for name, vals in fields.items():
        if name not in CLASSIFIER_FIELDS:
            continue
        if vals.size != n_cells:
            continue
        med = float(np.median(vals))
        mad = float(np.median(np.abs(vals - med)))
        deviations = np.abs(vals - med)
        # If MAD is too small relative to the field's magnitude, the
        # field is essentially uniform; MAD·5 ≈ 0 and every non-median
        # cell would be flagged as "outlier".  Skip such fields.
        rel_mad_threshold = 1e-6 * max(abs(med), 1.0)
        if mad < rel_mad_threshold:
            field_stats[name] = {
                "median": med, "mad": mad,
                "n_outliers": 0, "max_dev": float(np.max(deviations)),
                "skipped_uniform": True,
            }
            continue
        is_out = deviations > mad_k * mad
        outlier_mask |= is_out
        field_stats[name] = {
            "median":  med,
            "mad":     mad,
            "n_outliers": int(np.sum(is_out)),
            "max_dev": float(np.max(deviations)),
        }
    return outlier_mask, field_stats


def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pvd", required=True,
                   help="Path to the fault-surface .pvd file")
    p.add_argument("--output-root", default="r_v92_outlier_detection",
                   help="Output directory for JSON report + outlier CSVs")
    p.add_argument("--mad-k", type=float, default=5.0,
                   help="Outlier threshold in MAD multiples (default 5)")
    p.add_argument("--max-cycles", type=int, default=0,
                   help="Cap on number of timesteps processed (0 = all)")
    args = p.parse_args()

    root = Path(args.output_root)
    root.mkdir(parents=True, exist_ok=True)

    pvd_path = Path(args.pvd)
    if not pvd_path.is_file():
        print(f"ERROR: PVD not found: {pvd_path}")
        sys.exit(1)

    pvd_dir = pvd_path.parent
    ds_list = parse_pvd(pvd_path)
    if args.max_cycles > 0:
        ds_list = ds_list[:args.max_cycles]
    print(f"Found {len(ds_list)} timesteps in {pvd_path}")

    per_cycle_report = []
    for k, (t, rel_path) in enumerate(ds_list):
        vtu_path = pvd_dir / rel_path
        is_pvtu = str(rel_path).endswith(".pvtu")

        if is_pvtu:
            # PVTU references many .vtu pieces; concatenate them.
            if not vtu_path.is_file():
                print(f"  skip (missing pvtu): {vtu_path}")
                continue
            pvtu_tree = ET.parse(str(vtu_path))
            pieces = [el.attrib["Source"] for el in pvtu_tree.iter("Piece")]
            all_centroids, all_fields = [], {}
            for piece in pieces:
                piece_path = pvd_dir / piece
                if not piece_path.is_file():
                    continue
                c, f = load_vtu_fields(piece_path)
                all_centroids.append(c)
                for fname, fvals in f.items():
                    all_fields.setdefault(fname, []).append(fvals)
            if not all_centroids:
                continue
            centroids = np.concatenate(all_centroids, axis=0)
            fields = {name: np.concatenate(vs, axis=0)
                      for name, vs in all_fields.items()}
        elif vtu_path.is_file():
            centroids, fields = load_vtu_fields(vtu_path)
        else:
            print(f"  skip (missing file): {vtu_path}")
            continue

        mask, stats = classify_outliers(centroids, fields, args.mad_k)
        n_out = int(np.sum(mask))
        n_total = centroids.shape[0]

        # Partition correlation (requires §5.2 rank-id patch).
        rank_corr = None
        shared_corr = None
        if "is_shared_face" in fields and fields["is_shared_face"].size == n_total:
            is_shared = fields["is_shared_face"].astype(bool)
            n_shared = int(np.sum(is_shared))
            n_out_shared = int(np.sum(mask & is_shared))
            shared_corr = {
                "n_shared":        n_shared,
                "n_out_shared":    n_out_shared,
                "frac_out_on_shared": (n_out_shared / max(n_out, 1)),
            }
        if "mpi_rank" in fields and fields["mpi_rank"].size == n_total:
            ranks = fields["mpi_rank"]
            rank_counts = {}
            for r in np.unique(ranks):
                m = ranks == r
                rank_counts[int(r)] = {
                    "n": int(np.sum(m)),
                    "n_out": int(np.sum(mask & m)),
                }
            rank_corr = rank_counts

        per_cycle_report.append({
            "cycle":        k,
            "time":         t,
            "n_cells":      n_total,
            "n_outliers":   n_out,
            "outlier_frac": n_out / max(n_total, 1),
            "field_stats":  stats,
            "partition":    {"shared": shared_corr, "rank": rank_corr},
        })

        # Save outlier centroid CSV for ParaView re-import.
        if n_out > 0:
            out_csv = root / f"outliers_cycle_{k:04d}.csv"
            with open(out_csv, "w") as fh:
                fh.write("x,z,y\n")
                for c in centroids[mask]:
                    fh.write(f"{c[0]},{c[2]},{c[1]}\n")

        print(f"  cycle {k:4d}  t={t:6.3f}s  n_out={n_out:5d}/{n_total:5d}"
              f"  {100*n_out/max(n_total,1):.2f}%")

    report = {
        "pvd":        str(pvd_path),
        "mad_k":      args.mad_k,
        "cycles":     per_cycle_report,
    }
    report_path = root / "outlier_report.json"
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"Wrote {report_path}")

    # Trend analysis: is outlier count growing with time?
    if len(per_cycle_report) >= 3:
        n_outs = np.array([c["n_outliers"] for c in per_cycle_report])
        times  = np.array([c["time"]       for c in per_cycle_report])
        if times.max() > times.min():
            slope = np.polyfit(times, n_outs, 1)[0]
            print(f"Outlier trend: {slope:+.2f} outliers / second")
            tptp = float(np.ptp(times))
            if slope > 0.1 * n_outs.mean() / tptp:
                print("  ⇒ MONOTONIC GROWTH → H-V92-G bulk amplification")
            else:
                print("  ⇒ stable → candidate H-V92-P (shared-fault) or Q (random)")
        else:
            print("(single time step; trend not computed)")


if __name__ == "__main__":
    main()
