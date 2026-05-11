#!/usr/bin/env python3
"""sensitivity_sweep.py — Run corefine_faults.py over a parameter grid
and tabulate surface mesh size + quality, finding the best config.

For each config, runs the pipeline into a per-config directory, reads
the manifest, and produces:
  * per-config tables of (n_faces, edge_median, q_min, q_median)
  * conformality matrix (pairs conformal / total pairs)
  * shortest-edge floor violation count
  * pickled best config based on a composite quality score
"""

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


# Sweep grid (PLAN_cgal_corefine_multifault.md Phase 4, adapted to
# only the surface-mesh stage we have implemented).
GRID = [
    {"name": "baseline",       "mesh_edge_size": 1500, "min_edge": 100, "polyline_spacing":  750},
    {"name": "coarser",        "mesh_edge_size": 2000, "min_edge": 100, "polyline_spacing": 1000},
    {"name": "finer",          "mesh_edge_size": 1000, "min_edge": 100, "polyline_spacing":  500},
    {"name": "tighter_floor",  "mesh_edge_size": 1500, "min_edge": 200, "polyline_spacing":  750},
    {"name": "fine_floor",     "mesh_edge_size": 1500, "min_edge":  50, "polyline_spacing":  750},
]


def run_one(here: Path, project: Path, cfg: dict, res: int,
            sweep_root: Path, verbose: bool = False) -> dict:
    """Run one config; return summary dict."""
    cfg_dir = sweep_root / cfg["name"]
    work_dir = cfg_dir / "work"
    out_dir  = cfg_dir / "out"
    if cfg_dir.exists():
        shutil.rmtree(cfg_dir)
    cfg_dir.mkdir(parents=True)

    cmd = [sys.executable, str(here / "corefine_faults.py"),
           "--res", str(res),
           "--mesh-edge-size", str(cfg["mesh_edge_size"]),
           "--min-edge",       str(cfg["min_edge"]),
           "--polyline-spacing", str(cfg["polyline_spacing"]),
           "--out-dir", str(out_dir),
           "--workdir", str(work_dir)]
    if verbose:
        cmd.append("--verbose")
    t0 = time.time()
    res_run = subprocess.run(cmd, capture_output=not verbose)
    elapsed = time.time() - t0

    summary = {
        "config": cfg,
        "elapsed_s": elapsed,
        "exit_code": res_run.returncode,
        "out_dir": str(out_dir),
    }

    manifest_path = out_dir / "manifest.json"
    if not manifest_path.exists():
        summary["status"] = "FAILED"
        summary["error"] = "no manifest written"
        if not verbose and res_run.stderr:
            summary["stderr_tail"] = res_run.stderr.decode("utf-8")[-1000:]
        return summary

    with open(manifest_path) as f:
        manifest = json.load(f)

    n_meshes = len(manifest["meshes"])
    n_pairs  = len(manifest["pairs"])
    n_conformal = sum(1 for p in manifest["pairs"] if p["conformal"])

    n_faces = [m["n_faces"] for m in manifest["meshes"]]
    edge_min  = [m["edge_stats"]["min"]    for m in manifest["meshes"]]
    edge_med  = [m["edge_stats"]["median"] for m in manifest["meshes"]]
    edge_max  = [m["edge_stats"]["max"]    for m in manifest["meshes"]]
    q_min   = [m["tri_q_stats"]["min"]    for m in manifest["meshes"]]
    q_med   = [m["tri_q_stats"]["median"] for m in manifest["meshes"]]
    n_below_floor = [m["n_edges_below_floor"] for m in manifest["meshes"]]

    summary["status"] = "OK" if res_run.returncode in (0, 7) else "FAILED"
    summary["n_meshes"] = n_meshes
    summary["n_pairs"]  = n_pairs
    summary["n_conformal_pairs"] = n_conformal
    summary["total_faces"] = int(sum(n_faces))
    summary["faces_per_mesh_min"] = int(min(n_faces))
    summary["faces_per_mesh_max"] = int(max(n_faces))
    summary["edge_overall_min"] = float(min(edge_min))
    summary["edge_overall_med"] = float(np.median(edge_med))
    summary["edge_overall_max"] = float(max(edge_max))
    summary["q_overall_min"] = float(min(q_min))
    summary["q_overall_med"] = float(np.median(q_med))
    summary["n_edges_below_floor_total"] = int(sum(n_below_floor))
    summary["per_mesh"] = manifest["meshes"]
    return summary


def composite_score(s: dict) -> float:
    """Higher is better.  Penalizes:
    - Non-conformal pairs (weighted heavily).
    - Mesh-size deviation from target (mesh_edge_size).
    - Bad worst-case quality (q_overall_min).
    Rewards reasonable face count (not too sparse, not too dense).
    """
    if s.get("status") != "OK":
        return -1e6
    target_edge = s["config"]["mesh_edge_size"]
    edge_dev = abs(s["edge_overall_med"] - target_edge) / target_edge
    n_pairs = max(1, s["n_pairs"])
    conformality = s["n_conformal_pairs"] / n_pairs
    q_score = max(0.0, s["q_overall_min"])  # 0..1
    # Composite: 0.5 * conformality + 0.3 * (1 - edge_dev) + 0.2 * q_score
    return (0.5 * conformality
            + 0.3 * max(0.0, 1.0 - edge_dev)
            + 0.2 * q_score)


def main() -> int:
    here = Path(__file__).resolve().parent
    project = here.parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--res", type=int, default=2000)
    ap.add_argument("--sweep-root", type=Path,
                    default=here / "work" / "sweep")
    ap.add_argument("--report", type=Path,
                    default=project / "document" / "REPORT_sensitivity_2000m.md")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--configs", nargs="*", default=None,
                    help="restrict sweep to these named configs (default: all)")
    args = ap.parse_args()

    args.sweep_root.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    runs = []
    grid = [c for c in GRID if (args.configs is None or c["name"] in args.configs)]
    print(f"running sweep over {len(grid)} config(s)...")
    for cfg in grid:
        print(f"\n=== {cfg['name']} (mesh_edge={cfg['mesh_edge_size']}m,"
              f" min_edge={cfg['min_edge']}m, spacing={cfg['polyline_spacing']}m) ===")
        s = run_one(here, project, cfg, args.res, args.sweep_root, args.verbose)
        s["score"] = composite_score(s)
        runs.append(s)
        print(f"  status={s['status']} elapsed={s['elapsed_s']:.2f}s"
              f" total_F={s.get('total_faces', '?')}"
              f" conformal_pairs={s.get('n_conformal_pairs', '?')}/{s.get('n_pairs', '?')}"
              f" q_min={s.get('q_overall_min', float('nan')):.3f}"
              f" score={s['score']:.3f}")

    # Pick best.
    best = max(runs, key=lambda x: x["score"])
    print(f"\n>>> best config: {best['config']['name']} (score={best['score']:.3f}) <<<")

    # Emit markdown report.
    lines = []
    lines.append(f"# Sensitivity sweep — CGAL corefine — SAFS {args.res} m")
    lines.append(f"")
    lines.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"")
    lines.append(f"## Per-config summary")
    lines.append(f"")
    lines.append(f"| config | mesh_edge | min_edge | spacing | n_F | conformal | edge_min | edge_med | edge_max | q_min | q_med | elapsed | score |")
    lines.append(f"|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s in runs:
        c = s["config"]
        if s["status"] != "OK":
            lines.append(f"| {c['name']} | {c['mesh_edge_size']} | {c['min_edge']} | {c['polyline_spacing']} | FAILED | — | — | — | — | — | — | {s['elapsed_s']:.1f}s | — |")
            continue
        lines.append(f"| {c['name']} | {c['mesh_edge_size']} | {c['min_edge']} | {c['polyline_spacing']} | {s['total_faces']} | {s['n_conformal_pairs']}/{s['n_pairs']} | {s['edge_overall_min']:.1f} | {s['edge_overall_med']:.1f} | {s['edge_overall_max']:.1f} | {s['q_overall_min']:.4f} | {s['q_overall_med']:.4f} | {s['elapsed_s']:.1f}s | {s['score']:.3f} |")
    lines.append(f"")
    lines.append(f"**Best config:** `{best['config']['name']}` "
                 f"(`mesh_edge_size={best['config']['mesh_edge_size']}m`, "
                 f"`min_edge={best['config']['min_edge']}m`, "
                 f"`polyline_spacing={best['config']['polyline_spacing']}m`)")
    lines.append(f"")
    lines.append(f"## Per-fault detail (best config: {best['config']['name']})")
    lines.append(f"")
    lines.append(f"| fault | V | F | edge_min | edge_med | edge_max | q_min | q_med | n_edges_below_{best['config']['min_edge']}m |")
    lines.append(f"|---|---|---|---|---|---|---|---|---|")
    for m in best["per_mesh"]:
        lines.append(f"| {m['basename'][:55]} | {m['n_verts']} | {m['n_faces']} "
                     f"| {m['edge_stats']['min']:.1f} "
                     f"| {m['edge_stats']['median']:.1f} "
                     f"| {m['edge_stats']['max']:.1f} "
                     f"| {m['tri_q_stats']['min']:.4f} "
                     f"| {m['tri_q_stats']['median']:.4f} "
                     f"| {m['n_edges_below_floor']} |")
    lines.append(f"")
    lines.append(f"## Score formula")
    lines.append(f"")
    lines.append(f"```")
    lines.append(f"score = 0.5 * (n_conformal_pairs / n_pairs)")
    lines.append(f"      + 0.3 * max(0, 1 - |edge_med - mesh_edge_size| / mesh_edge_size)")
    lines.append(f"      + 0.2 * q_overall_min")
    lines.append(f"```")
    lines.append(f"")

    args.report.write_text("\n".join(lines))
    print(f"\nreport written to {args.report}")

    # Also dump runs as JSON for downstream tools.
    json_path = args.report.with_suffix(".json")
    # Strip per_mesh detail for compactness.
    runs_compact = []
    for s in runs:
        s2 = {k: v for k, v in s.items() if k != "per_mesh"}
        runs_compact.append(s2)
    json_path.write_text(json.dumps({"runs": runs_compact,
                                      "best": best["config"]["name"]},
                                     indent=2))
    print(f"json:   {json_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
