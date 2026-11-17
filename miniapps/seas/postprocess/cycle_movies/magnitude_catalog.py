#!/usr/bin/env python3
"""Moment magnitudes of the events of a 2-D fully dynamic cycle run.

    python3 magnitude_catalog.py --run RUN_FOLDER --out OUT_FOLDER --mu 32.04e9 --width-m 15e3

writes magnitude_catalog.csv (one row per catalog row) and
magnitude_catalog.png (Mw against event time, one stem per row, coloured by
the fault that carries the largest slip of the event).

A 2-D run has no rupture width, so the moment needs one from outside, W
(--width-m; for a 2.5-D run the deck's [cycle.finite_width] W_m):

    M0 = sum over fault cells of  mu * W * cell length * mean slip of the event at its 2 points
    Mw = (2/3) (log10 M0 - 9.1)                                   (M0 in N m)

mu (--mu, Pa) is the shear modulus at the fault; one value for the whole
fault (a run with a material that varies along the fault needs more than
this script).

Slip of an event = |slip(end) - slip(start)| (strike and dip components)
between two fault snapshots:
  end    the first snapshot at or after the row's t_end;
  start  the last snapshot at or before t_start when it is within
         --start-tol seconds of it (an event that starts inside a dynamic
         phase), else the first snapshot at or after t_start.  The driver
         writes no snapshot between t_start and the switch to the dynamic
         regime, so for an event that opens a phase the start is the
         snapshot at the switch and the slip before it is not counted; the
         CSV gives that gap per row (start_gap_s).

Faults: the connected groups of fault cells, numbered from the one that
reaches farthest to -x; --fault-names gives them names in that order.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cycle_run  # noqa: E402


def fault_groups(points, cells, tol_km=1e-4):
    """Group id of every cell: cells that share a point location are one fault."""
    parent = list(range(len(cells)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    seen = {}
    for c, nodes in enumerate(cells):
        for q in nodes:
            key = tuple(np.round(points[q] / tol_km).astype(np.int64))
            if key in seen:
                parent[root(c)] = root(seen[key])
            else:
                seen[key] = c
    roots = np.array([root(c) for c in range(len(cells))])
    order = sorted(set(roots.tolist()), key=lambda r: points[cells[roots == r]].reshape(-1, 2)[:, 0].min())
    return np.array([order.index(r) for r in roots]), len(order)


def events(run, mu, width_m, start_tol=0.5):
    """One dict per catalog row."""
    s = cycle_run.Series(cycle_run.find(run, "_fault_snapshots.vtkhdf"), ["slip_strike", "slip_dip"])
    cells = s.cells(2)
    length_m = 1e3 * np.hypot(*(s.points[cells[:, 1]] - s.points[cells[:, 0]]).T)
    group, ngroups = fault_groups(s.points, cells)
    point_group = np.empty(s.npts, dtype=np.int64)
    point_group[cells] = group[:, None]
    tv = s.times[:s.nsteps]
    t_on = np.array(cycle_run.switch_times(run))

    def slip(k):
        return np.stack([s.field("slip_strike", k), s.field("slip_dip", k)], axis=1)

    rows = []
    for ev, t0, t1, vpk in cycle_run.catalog(run):
        kb = int(np.searchsorted(tv, t0, side="right") - 1)
        k0 = kb if kb >= 0 and t0 - tv[kb] <= start_tol else int(np.searchsorted(tv, t0, side="left"))
        k1 = int(np.searchsorted(tv, t1, side="left"))
        if k0 >= len(tv) or k1 >= len(tv):
            continue  # the snapshots of this event are not on disk (a run still going)
        d = np.hypot(*(slip(k1) - slip(k0)).T)
        m0 = mu * width_m * float(np.sum(length_m * d[cells].mean(1)))
        j = int(np.argmax(d))
        rows.append(dict(row=int(ev), phase=int(np.searchsorted(t_on, t1, side="right")), t_start_yr=t0 / cycle_run.YEAR_S,
                         duration_s=t1 - t0, peak_vmax=vpk, max_slip_m=float(d.max()), fault=int(point_group[j]),
                         x_km=float(s.points[j, 0]), y_km=float(s.points[j, 1]), m0=m0,
                         mw=(2.0 / 3.0) * (np.log10(m0) - 9.1) if m0 > 0 else float("nan"),
                         start_gap_s=float(tv[k0] - t0)))
    lengths = [float(length_m[group == g].sum() / 1e3) for g in range(ngroups)]
    return rows, lengths


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--mu", type=float, required=True, help="shear modulus at the fault, Pa")
    p.add_argument("--width-m", type=float, required=True, help="rupture width W, m")
    p.add_argument("--fault-names", help="comma-separated names of the faults, ordered by their smallest x")
    p.add_argument("--start-tol", type=float, default=0.5, help="see above, s (default 0.5)")
    p.add_argument("--title", default=None, help="small title of the figure (default: the run folder's name)")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows, lengths = events(a.run, a.mu, a.width_m, a.start_tol)
    if not rows:
        raise SystemExit("no catalog row has its snapshots on disk")
    names = a.fault_names.split(",") if a.fault_names else ["fault %d" % (g + 1) for g in range(len(lengths))]
    if len(names) != len(lengths):
        raise SystemExit("--fault-names has %d names, the fault snapshots have %d faults" % (len(names), len(lengths)))
    csv = os.path.join(a.out, "magnitude_catalog.csv")
    with open(csv, "w") as f:
        f.write("# mu %.6e Pa, W %.6e m; M0 = sum mu W length slip; Mw = (2/3)(log10 M0 - 9.1)\n" % (a.mu, a.width_m))
        f.write("row,dynamic_phase,t_start_yr,duration_s,peak_Vmax_m_s,max_event_slip_m,fault,x_max_km,y_max_km,"
                "M0_Nm,Mw,start_gap_s\n")
        for r in rows:
            f.write("%d,%d,%.3f,%.1f,%.4g,%.4f,%s,%.2f,%.2f,%.5e,%.4f,%.2f\n"
                    % (r["row"], r["phase"], r["t_start_yr"], r["duration_s"], r["peak_vmax"], r["max_slip_m"],
                       names[r["fault"]], r["x_km"], r["y_km"], r["m0"], r["mw"], r["start_gap_s"]))
    mw = np.array([r["mw"] for r in rows])
    t = np.array([r["t_start_yr"] for r in rows])
    fl = np.array([r["fault"] for r in rows])
    print("%d rows in %d dynamic phases; Mw min %.2f median %.2f max %.2f"
          % (len(rows), len(set(r["phase"] for r in rows)), np.nanmin(mw), np.nanmedian(mw), np.nanmax(mw)))
    for g, name in enumerate(names):
        print("  %s: %.1f km, %d events" % (name, lengths[g], int((fl == g).sum())))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    unit, scale = ("kyr", 1e3) if t.max() > 5e3 else ("yr", 1.0)
    fig, ax = plt.subplots(figsize=(9, 5), dpi=200)
    colours = plt.get_cmap("tab10").colors
    for g, name in enumerate(names):
        m = fl == g
        if m.any():
            ax.vlines(t[m] / scale, 0, mw[m], color=colours[g % 10], lw=1.2)
            ax.plot(t[m] / scale, mw[m], "o", ms=4, color=colours[g % 10], label="%s (%d)" % (name, m.sum()))
    ax.set_xlabel("event start time (%s)" % unit)
    ax.set_ylabel("Mw")
    ax.set_ylim(0, np.ceil(np.nanmax(mw) * 2) / 2 + 0.25)
    ax.set_xlim(0, None)
    ax.legend(title="fault with the largest slip of the event", fontsize=8, title_fontsize=8, loc="lower left",
              bbox_to_anchor=(0.0, 1.01), ncol=min(5, len(names)), frameon=False)
    # the whole run name, small, on its own line above the legend
    fig.suptitle(a.title or os.path.basename(os.path.normpath(a.run)), fontsize=6, x=0.01, y=0.995, ha="left",
                 va="top")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    png = os.path.join(a.out, "magnitude_catalog.png")
    fig.savefig(png)
    print(csv)
    print(png)


if __name__ == "__main__":
    main()
