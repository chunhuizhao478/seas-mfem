#!/usr/bin/env python3
"""
plot_fault_normal_profile.py — off-fault fault-normal-profile ground-motion plots.

Uses the 9 receivers along the fault-NORMAL line through the epicenter
(epicenter + fault_normal_{NE,SW}_{2,5,10,20}km) from two SeisSol runs and produces:
  fig1  PGV vs fault-normal distance (log-log, NE/SW split, CVM vs constant)
  fig2  component peaks (fault-parallel / fault-normal / vertical) vs signed offset
  fig3  permanent fault-parallel displacement vs signed offset  (the static step)
  fig4  fault-normal velocity waveforms stacked across the profile (one run)

Velocity components are rotated into the LOCAL fault frame at the epicenter:
  fault-normal n0 = unit(NE_20km - SW_20km)   (the profile was built along n0)
  fault-parallel t0 = n0 rotated -90 deg in the map plane;  vertical = v3.

Run:  python3 plot_fault_normal_profile.py
      [--cvm DIR] [--const DIR] [--out DIR] [--prefix safs]
"""
import argparse
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# the 9 fault-normal-profile stations: label -> (x, y) UTM, and signed nominal offset [km]
PROFILE = {
    "epicenter":            (604801.2000, 3705030.0000,   0.0),
    "fault_normal_SW_20km": (590745.5748, 3690801.8799, -20.0),
    "fault_normal_SW_10km": (597773.3874, 3697915.9400, -10.0),
    "fault_normal_SW_5km":  (601287.2937, 3701472.9700,  -5.0),
    "fault_normal_SW_2km":  (603395.6375, 3703607.1880,  -2.0),
    "fault_normal_NE_2km":  (606206.7625, 3706452.8120,   2.0),
    "fault_normal_NE_5km":  (608315.1063, 3708587.0300,   5.0),
    "fault_normal_NE_10km": (611829.0126, 3712144.0600,  10.0),
    "fault_normal_NE_20km": (618856.8252, 3719258.1201,  20.0),
}
EPI = np.array(PROFILE["epicenter"][:2])
# local fault frame at the epicenter (profile built along n0 = NE20 - SW20)
_n0 = np.array(PROFILE["fault_normal_NE_20km"][:2]) - np.array(PROFILE["fault_normal_SW_20km"][:2])
N0 = _n0 / np.linalg.norm(_n0)               # fault-normal (points NE)
T0 = np.array([N0[1], -N0[0]])               # fault-parallel (strike), n0 rotated -90 deg

# along-strike (directivity) stations: label -> (x, y, nominal arc length s [km] NW of epicenter)
ALONGSTRIKE = {
    "epicenter":           (604801.2000, 3705030.0000,  0.0),
    "along_strike_NW_10km": (595924.3549, 3711699.7544, 10.0),
    "along_strike_NW_20km": (588836.3979, 3718720.1024, 20.0),
    "along_strike_NW_40km": (574324.2424, 3732355.1845, 40.0),
}


def load_trace(mesh_path, fault_bc=3, z_tol=1.0):
    """Surface fault-trace vertices (xy) from the PUML mesh, for local-strike estimation."""
    import h5py
    with h5py.File(mesh_path, "r") as f:
        geom = f["geometry"][:]
        connect = f["connect"][:].astype(np.int64)
        boundary = f["boundary"][:].astype(np.int64)
    face_sets = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]])
    vids = []
    for j in range(4):
        sel = np.where(((boundary >> (8 * j)) & 0xFF) == fault_bc)[0]
        if sel.size:
            vids.append(connect[sel][:, face_sets[j]].ravel())
    fv = geom[np.unique(np.concatenate(vids))]
    return fv[np.abs(fv[:, 2]) < z_tol][:, :2]


def local_strike(trace, p, radius=3000.0):
    """Unit along-strike direction from trace points within `radius` of map point p (xy)."""
    d = np.linalg.norm(trace - p, axis=1)
    nb = trace[d < radius]
    if nb.shape[0] < 2:
        nb = trace[np.argsort(d)[:5]]
    _, _, vt = np.linalg.svd(nb - nb.mean(axis=0), full_matrices=False)
    t = vt[0]
    return t / np.linalg.norm(t)


def cumtrapz0(y, t):
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))
    return out


def load_profile(run_dir, prefix):
    """Match each receiver .dat to a profile station by coordinate; return dict label -> data."""
    files = sorted(glob.glob(os.path.join(run_dir, f"{prefix}-receiver-*.dat")))
    out = {}
    for f in files:
        x = {}
        with open(f) as fh:
            for line in fh.readlines()[:8]:      # x1/x2/x3 sit after TITLE+VARIABLES
                m = re.match(r"#\s*x([123])\s+([-+0-9.eE]+)", line)
                if m:
                    x[int(m.group(1))] = float(m.group(2))
        xy = np.array([x[1], x[2]])
        # nearest profile station
        best, bestd = None, 1e30
        for label, (sx, sy, _) in PROFILE.items():
            d = np.hypot(xy[0] - sx, xy[1] - sy)
            if d < bestd:
                best, bestd = label, d
        if bestd > 50.0:            # not a profile station (e.g. along-strike receiver)
            continue
        data = np.loadtxt(f, comments="#", skiprows=2)
        t = data[:, 0]
        v = data[:, 7:10]
        par = v[:, 0] * T0[0] + v[:, 1] * T0[1]   # fault-parallel
        nor = v[:, 0] * N0[0] + v[:, 1] * N0[1]   # fault-normal
        ver = v[:, 2]                              # vertical
        out[best] = dict(
            t=t, par=par, nor=nor, ver=ver,
            offset=PROFILE[best][2],
            pgv3d=float(np.max(np.sqrt(par**2 + nor**2 + ver**2))),
            peak_par=float(np.max(np.abs(par))),
            peak_nor=float(np.max(np.abs(nor))),
            peak_ver=float(np.max(np.abs(ver))),
            perm_par=float(cumtrapz0(par, t)[-1]),
        )
    return out


def split_sides(prof, key):
    """Return (|offset|, value) arrays for NE (offset>0) and SW (offset<0), sorted by |offset|."""
    ne = sorted([(abs(d["offset"]), d[key]) for d in prof.values() if d["offset"] > 0])
    sw = sorted([(abs(d["offset"]), d[key]) for d in prof.values() if d["offset"] < 0])
    a = lambda L: (np.array([p[0] for p in L]), np.array([p[1] for p in L]))
    return a(ne), a(sw)


def signed(prof, key, drop_epi=False):
    pts = sorted([(d["offset"], d[key]) for d in prof.values()
                  if not (drop_epi and d["offset"] == 0.0)])
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def load_alongstrike(run_dir, prefix, trace):
    """Match along-strike receivers; rotate v into the LOCAL fault frame (per-station strike)."""
    files = sorted(glob.glob(os.path.join(run_dir, f"{prefix}-receiver-*.dat")))
    out = {}
    for f in files:
        x = {}
        for line in open(f).readlines()[:8]:
            m = re.match(r"#\s*x([123])\s+([-+0-9.eE]+)", line)
            if m:
                x[int(m.group(1))] = float(m.group(2))
        xy = np.array([x[1], x[2]])
        best, bestd = None, 1e30
        for label, (sx, sy, _) in ALONGSTRIKE.items():
            d = np.hypot(xy[0] - sx, xy[1] - sy)
            if d < bestd:
                best, bestd = label, d
        if bestd > 50.0:
            continue
        if trace is not None:
            t_hat = local_strike(trace, xy)
            n_hat = np.array([-t_hat[1], t_hat[0]])
        else:
            t_hat, n_hat = T0, N0
        data = np.loadtxt(f, comments="#", skiprows=2)
        t = data[:, 0]
        v = data[:, 7:10]
        nor = v[:, 0] * n_hat[0] + v[:, 1] * n_hat[1]
        par = v[:, 0] * t_hat[0] + v[:, 1] * t_hat[1]
        peak_nor = float(np.max(np.abs(nor)))
        out[best] = dict(
            t=t, nor=nor, par=par, ver=v[:, 2], s=ALONGSTRIKE[best][2],
            peak_nor=peak_nor, peak_par=float(np.max(np.abs(par))),
            pgv3d=float(np.max(np.sqrt(par**2 + nor**2 + v[:, 2]**2))),
            t_peak=float(t[np.argmax(np.abs(nor))]),
        )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cvm", default="/Users/chunhuizhao/Downloads/output_safs_v2.2.0_cvm_mat")
    ap.add_argument("--const",
                    default="/Users/chunhuizhao/Downloads/output_safs_v2.2.0_constant_mat_onfaultpoints")
    ap.add_argument("--out", default="figures_offfault")
    ap.add_argument("--prefix", default="safs")
    ap.add_argument("--mesh", default="safs_mesh.puml.h5",
                    help="PUML mesh for local-strike rotation of along-strike stations")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    runs = {}
    if os.path.isdir(a.cvm):
        runs["CVM"] = load_profile(a.cvm, a.prefix)
    if os.path.isdir(a.const):
        runs["constant"] = load_profile(a.const, a.prefix)
    if not runs:
        raise SystemExit("no run dirs found")
    style = {"CVM": dict(color="C3", marker="o"), "constant": dict(color="C0", marker="s")}

    # ---- fig1: PGV vs fault-normal distance (log-log, NE/SW split) ----
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for name, prof in runs.items():
        ne, sw = split_sides(prof, "pgv3d")
        ax.loglog(ne[0], ne[1], "-",  **style[name], label=f"{name} NE")
        ax.loglog(sw[0], sw[1], "--", color=style[name]["color"],
                  marker=style[name]["marker"], mfc="none", label=f"{name} SW")
    # reference geometric-spreading slopes, anchored at the 2 km CVM-NE point
    anchor_x, anchor_y = 2.0, None
    if "CVM" in runs:
        ne, _ = split_sides(runs["CVM"], "pgv3d")
        if len(ne[0]):
            anchor_y = ne[1][0]
    if anchor_y:
        rr = np.array([2.0, 20.0])
        for p, ls in [(-1.0, ":"), (-1.5, "--")]:
            ax.loglog(rr, anchor_y * (rr / anchor_x) ** p, ls, color="0.5", lw=1.2,
                      label=f"$R^{{{p:.1f}}}$ (ref)")
    ax.set_xlabel("fault-normal distance |offset|  [km]")
    ax.set_ylabel("PGV  (3-D peak)  [m/s]")
    ax.set_title("Off-fault PGV vs fault-normal distance")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig1_pgv_vs_distance.png"), dpi=150)
    plt.close(fig)

    # ---- fig2: component peaks vs signed offset ----
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    comps = [("peak_par", "fault-parallel"), ("peak_nor", "fault-normal"), ("peak_ver", "vertical")]
    for axc, (key, title) in zip(axes, comps):
        for name, prof in runs.items():
            o, val = signed(prof, key, drop_epi=True)
            axc.plot(o, val, "-", **style[name], label=name)
        axc.axvline(0, color="k", lw=0.8, ls=":")
        axc.set_xlabel("signed fault-normal offset  [km]  (SW<0, NE>0)")
        axc.set_title(f"peak {title} velocity")
        axc.grid(True, alpha=0.3)
    axes[0].set_ylabel("peak velocity  [m/s]")
    axes[0].legend(fontsize=8)
    fig.suptitle("Peak velocity by component vs fault-normal offset  (on-trace epicenter excluded)")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig2_components_vs_offset.png"), dpi=150)
    plt.close(fig)

    # ---- fig3: permanent fault-parallel displacement vs signed offset ----
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for name, prof in runs.items():
        o, val = signed(prof, "perm_par", drop_epi=True)
        ax.plot(o, val, "-", **style[name], label=name)
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("signed fault-normal offset  [km]  (SW<0, NE>0)")
    ax.set_ylabel("permanent fault-parallel displacement  [m]")
    ax.set_title("Static (coseismic) fault-parallel offset across the fault  (epicenter excluded)")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig3_static_offset.png"), dpi=150)
    plt.close(fig)

    # ---- fig4: fault-normal velocity waveforms stacked (CVM if present else first run) ----
    name = "CVM" if "CVM" in runs else next(iter(runs))
    prof = runs[name]
    order = ["fault_normal_SW_20km", "fault_normal_SW_10km", "fault_normal_SW_5km",
             "fault_normal_SW_2km", "epicenter", "fault_normal_NE_2km",
             "fault_normal_NE_5km", "fault_normal_NE_10km", "fault_normal_NE_20km"]
    order = [k for k in order if k in prof]
    amp = max(prof[k]["peak_nor"] for k in order) or 1.0
    fig, ax = plt.subplots(figsize=(8, 7))
    for i, k in enumerate(order):
        d = prof[k]
        ax.plot(d["t"], d["nor"] / amp * 1.2 + i, lw=0.8, color="C3")
        ax.text(d["t"][-1] * 0.98, i + 0.25, f"{k.replace('fault_normal_','')}  "
                f"(PGV_n={d['peak_nor']:.2f})", fontsize=7, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{prof[k]['offset']:+.0f} km" for k in order])
    ax.set_xlabel("time  [s]"); ax.set_ylabel("station (signed offset)")
    ax.set_title(f"Fault-normal velocity waveforms across the profile  ({name})")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig4_waveforms_fault_normal.png"), dpi=150)
    plt.close(fig)

    # ---- along-strike directivity (NW line) -----------------------------------
    trace = None
    try:
        if os.path.exists(a.mesh):
            trace = load_trace(a.mesh)
            print(f"loaded fault trace ({trace.shape[0]} surface vertices) for local-strike rotation")
        else:
            print(f"{a.mesh} not found — along-strike rotation falls back to epicenter frame")
    except Exception as e:
        print(f"trace load failed ({e}) — falling back to epicenter frame")

    as_runs = {}
    for name, d in (("CVM", a.cvm), ("constant", a.const)):
        if os.path.isdir(d):
            as_runs[name] = load_alongstrike(d, a.prefix, trace)

    if as_runs:
        # fig5: peak fault-normal velocity & pulse arrival-time vs along-strike distance
        fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.6))
        for name, prof in as_runs.items():
            pts = sorted(prof.values(), key=lambda d: d["s"])
            s = np.array([d["s"] for d in pts])
            axA.plot(s, [d["peak_nor"] for d in pts], "-", **style[name], label=name)
            # apparent velocity from peak-arrival moveout (skip s=0 on-trace anchor)
            far = [d for d in pts if d["s"] > 0]
            if len(far) >= 2:
                sf = np.array([d["s"] for d in far]); tf = np.array([d["t_peak"] for d in far])
                axB.plot(sf, tf, "-", **style[name], label=name)
                slope = np.polyfit(sf, tf, 1)[0]            # s[km]/t[s] inverse
                vapp = 1.0 / slope if slope else float("inf")
                axB.text(0.05, 0.9 if name == "CVM" else 0.8,
                         f"{name}: v_app ≈ {vapp:.1f} km/s", transform=axB.transAxes,
                         color=style[name]["color"], fontsize=9)
        axA.set_xlabel("along-strike distance NW  [km]")
        axA.set_ylabel("peak fault-normal velocity  [m/s]")
        axA.set_title("Directivity: peak fault-normal velocity along strike")
        axA.grid(True, alpha=0.3); axA.legend(fontsize=8)
        axB.set_xlabel("along-strike distance NW  [km]")
        axB.set_ylabel("fault-normal pulse peak time  [s]")
        axB.set_title("Pulse arrival moveout (slope⁻¹ = apparent velocity)")
        axB.grid(True, alpha=0.3); axB.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig5_along_strike_directivity.png"), dpi=150)
        plt.close(fig)

        # fig6: along-strike fault-normal waveform overlay (CVM if present)
        name = "CVM" if "CVM" in as_runs else next(iter(as_runs))
        prof = as_runs[name]
        order = ["epicenter", "along_strike_NW_10km", "along_strike_NW_20km", "along_strike_NW_40km"]
        order = [k for k in order if k in prof]
        amp = max(prof[k]["peak_nor"] for k in order) or 1.0
        fig, ax = plt.subplots(figsize=(8, 5.5))
        for i, k in enumerate(order):
            d = prof[k]
            ax.plot(d["t"], d["nor"] / amp * 1.2 + i, lw=0.9, color="C3")
            ax.text(d["t"][-1] * 0.98, i + 0.25,
                    f"s={d['s']:.0f} km  (PGV_n={d['peak_nor']:.2f}, t_pk={d['t_peak']:.1f}s)",
                    fontsize=8, ha="right")
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([f"{prof[k]['s']:.0f} km" for k in order])
        ax.set_xlabel("time  [s]"); ax.set_ylabel("along-strike distance NW")
        ax.set_title(f"Fault-normal velocity along the strike line  ({name})")
        ax.grid(True, axis="x", alpha=0.3)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig6_along_strike_waveforms.png"), dpi=150)
        plt.close(fig)

        print("\n== along-strike directivity ==")
        for name, prof in as_runs.items():
            print(f"-- {name}:  " + "  ".join(
                f"s={d['s']:.0f}km pk_n={d['peak_nor']:.2f} t_pk={d['t_peak']:.1f}s"
                for d in sorted(prof.values(), key=lambda d: d["s"])))

    # ---- console summary table ----
    print(f"fault frame: n0(normal)={N0.round(3)}  t0(parallel)={T0.round(3)}")
    for name, prof in runs.items():
        print(f"\n== {name} ==")
        print(f"{'offset_km':>9} {'PGV3d':>8} {'pk_par':>8} {'pk_norm':>8} {'pk_vert':>8} {'perm_par':>9}")
        for k in sorted(prof, key=lambda k: prof[k]["offset"]):
            d = prof[k]
            print(f"{d['offset']:>9.0f} {d['pgv3d']:>8.3f} {d['peak_par']:>8.3f} "
                  f"{d['peak_nor']:>8.3f} {d['peak_ver']:>8.3f} {d['perm_par']:>9.3f}")
    print(f"\nwrote 4 figures to {a.out}/")


if __name__ == "__main__":
    main()
