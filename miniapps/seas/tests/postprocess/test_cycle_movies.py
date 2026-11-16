#!/usr/bin/env python3
"""Checks of postprocess/cycle_movies (python3 tests/postprocess/test_cycle_movies.py, about a minute).

A small made-up run folder with the layout of seas_fd_cycle_driver's output
(two straight faults, three dynamic phases, known fields) is written to a
temporary folder, and the scripts are run on it:
  1. raster.locate: a linear field is reproduced at every pixel centre; a
     triangle smaller than a pixel is found; pixels outside the mesh are empty
  2. basemap: the projection keeps distance and bearing from its centre
  3. cycle_run.dynamic_frames: the frames of each phase, the cadence, and
     the last phase only with DONE
  4. cycle_movies: frame counts; 1 worker and 3 workers give the same bytes;
     the parts of --max-mb hold every frame once; the line map version runs;
     drawing only the changing layers gives the files of --full-redraw
  5. magnitude_catalog: Mw of a uniform slip against the formula
Exit 0 when every check passes.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.normpath(os.path.join(HERE, "..", "..", "postprocess", "cycle_movies"))
sys.path.insert(0, TOOL)
import basemap  # noqa: E402
import cycle_movies  # noqa: E402
import cycle_run  # noqa: E402
import magnitude_catalog  # noqa: E402
import raster  # noqa: E402

failed = 0


def check(ok, text):
    global failed
    failed += 0 if ok else 1
    print("  %s: %s" % ("ok" if ok else "FAILED", text))


def write_series(path, points, cells, fields, times, regime, phase):
    """A VTKHDF series as cycle/snapshots.hpp writes it (the parts the readers use)."""
    n = len(points)
    with h5py.File(path, "w") as f:
        w = f.create_group("VTKHDF")
        w.create_dataset("Points", data=np.column_stack([points * 1e3, np.zeros(n)]))
        w.create_dataset("Connectivity", data=cells.ravel().astype(np.int64))
        w.create_dataset("Offsets", data=np.arange(len(cells) + 1, dtype=np.int64) * cells.shape[1])
        w.create_dataset("Steps/Values", data=np.asarray(times, float))
        for name, per_step in fields.items():
            data = np.concatenate(per_step)
            w.create_dataset("PointData/" + name, data=data, chunks=(n,) + data.shape[1:])
            w.create_dataset("Steps/PointDataOffsets/" + name, data=np.arange(len(times), dtype=np.int64) * n)
        f.create_dataset("fd_cycle/regime", data=np.asarray(regime, np.int64))
        f.create_dataset("fd_cycle/phase", data=np.asarray(phase, np.int64))


def box_mesh(nx, ny, x0, x1, y0, y1):
    """Triangles of a box, every triangle with its own three points (as the volume snapshots have)."""
    xs, ys = np.linspace(x0, x1, nx + 1), np.linspace(y0, y1, ny + 1)
    pts = []
    for j in range(ny):
        for i in range(nx):
            a, b, c, d = (xs[i], ys[j]), (xs[i + 1], ys[j]), (xs[i + 1], ys[j + 1]), (xs[i], ys[j + 1])
            pts += [a, b, c, a, c, d]
    pts = np.array(pts)
    return pts, np.arange(len(pts)).reshape(-1, 3)


T_ON = [1.0e9, 2.0e9, 3.0e9]        # switch times of the three phases, s
NDYN = [30, 22, 26]                 # dynamic volume snapshots per phase, 2 s apart
MU, W_M = 30.0e9, 12.0e3
FAULT_KM = 40.0                     # length of fault 1, the one that slips
SLIP = [2.0, 0.5, 4.0]              # uniform slip of fault 1 in each event, m


def make_run(run, done):
    os.makedirs(run)
    pts, tris = box_mesh(24, 16, -60.0, 60.0, -40.0, 40.0)
    times, regime, phase, vel, disp = [0.0], [0], [0], [np.zeros((len(pts), 3))], [np.zeros((len(pts), 3))]
    for p, (t0, n) in enumerate(zip(T_ON, NDYN), 1):
        times.append(t0); regime.append(0); phase.append(p - 1)  # the QD snapshot at the switch
        vel.append(np.zeros((len(pts), 3))); disp.append(np.zeros((len(pts), 3)))
        for k in range(n):
            t = 2.0 * (k + 1)
            r = np.hypot(pts[:, 0], pts[:, 1])
            v = np.exp(-((r - 3.0 * t) / 8.0) ** 2)
            times.append(t0 + t); regime.append(1); phase.append(p)
            vel.append(np.column_stack([v, 0.5 * v, np.zeros(len(pts))])); disp.append(np.zeros((len(pts), 3)))
        times.append(t0 + 2.0 * n + 1.0); regime.append(0); phase.append(p)
        vel.append(np.zeros((len(pts), 3))); disp.append(np.zeros((len(pts), 3)))
    if not done:  # a run that stopped inside the last phase: drop the QD step after it
        times, regime, phase, vel, disp = times[:-1], regime[:-1], phase[:-1], vel[:-1], disp[:-1]
    write_series(os.path.join(run, "t_volume_snapshots.vtkhdf"), pts, tris,
                 {"velocity": vel, "displacement": disp}, times, regime, phase)
    # two faults: fault 1 along y = 0 (40 km, 40 cells), fault 2 a short one that never slips
    x = np.linspace(-20.0, 20.0, 41)
    f1 = np.array([[(x[i], 0.0), (x[i + 1], 0.0)] for i in range(40)]).reshape(-1, 2)
    f2 = np.array([[(30.0 + i, 20.0), (31.0 + i, 20.0)] for i in range(5)]).reshape(-1, 2)
    fp = np.vstack([f1, f2])
    fc = np.arange(len(fp)).reshape(-1, 2)
    on1 = np.arange(len(fp)) < len(f1)
    ft, fr, fph, ss, sd, rs, rd = [0.0], [0], [0], [np.zeros(len(fp))], [np.zeros(len(fp))], [np.zeros(len(fp))], [np.zeros(len(fp))]
    total = 0.0
    for p, (t0, n) in enumerate(zip(T_ON, NDYN), 1):
        ft.append(t0); fr.append(0); fph.append(p - 1)
        ss.append(on1 * total); sd.append(np.zeros(len(fp))); rs.append(np.zeros(len(fp))); rd.append(np.zeros(len(fp)))
        m = 4 * n  # fault snapshots 0.5 s apart
        for k in range(m):
            frac = (k + 1) / m
            ft.append(t0 + 0.5 * (k + 1)); fr.append(1); fph.append(p)
            ss.append(on1 * (total + SLIP[p - 1] * frac)); sd.append(np.zeros(len(fp)))
            rs.append(on1 * SLIP[p - 1] / (0.5 * m)); rd.append(np.zeros(len(fp)))
        total += SLIP[p - 1]
        ft.append(t0 + 0.5 * m + 1.0); fr.append(0); fph.append(p)
        ss.append(on1 * total); sd.append(np.zeros(len(fp))); rs.append(np.zeros(len(fp))); rd.append(np.zeros(len(fp)))
    if not done:
        ft, fr, fph, ss, sd, rs, rd = ft[:-1], fr[:-1], fph[:-1], ss[:-1], sd[:-1], rs[:-1], rd[:-1]
    write_series(os.path.join(run, "t_fault_snapshots.vtkhdf"), fp, fc,
                 {"slip_strike": ss, "slip_dip": sd, "slip_rate_strike": rs, "slip_rate_dip": rd}, ft, fr, fph)
    with open(os.path.join(run, "t_switch.txt"), "w") as f:
        f.write("# fd_cycle switch log\n# Columns: time(s), Vmax(m/s), from, to, reason\n")
        for t0, n in zip(T_ON, NDYN):
            f.write("%.12E 1.0E-02 QD dynamic V_on\n%.12E 1.0E-03 dynamic dynamic quench_on\n" % (t0, t0 + 5.0))
            f.write("%.12E 1.0E-04 dynamic QD end_of_phase\n" % (t0 + 2.0 * n))
    with open(os.path.join(run, "t_catalog.txt"), "w") as f:
        f.write("# fd_cycle event catalog\n# Columns: event, t_start(s), t_end(s), peak_Vmax(m/s)\n")
        for e, (t0, n) in enumerate(zip(T_ON, NDYN), 1):
            f.write("%.12E %.12E %.12E 1.0E+00\n" % (e, t0 - 30.0, t0 + 2.0 * n))
    if done:
        open(os.path.join(run, "DONE"), "w").close()


def run_tool(script, *args, quiet=False):
    r = subprocess.run([sys.executable, os.path.join(TOOL, script)] + list(args), stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True)
    if r.returncode != 0 and not quiet:
        print(r.stdout[-2000:])
    return r.returncode, r.stdout


def main():
    tmp = tempfile.mkdtemp(prefix="seas_test_cycle_movies_")
    print("1. raster.locate")
    pts, tris = box_mesh(13, 7, 0.0, 10.0, 0.0, 6.0)
    # a triangle smaller than a pixel, off to the side, around the pixel centre (24.25, 3.1)
    tiny = np.array([[24.20, 3.05], [24.30, 3.05], [24.25, 3.20]])
    pts2 = np.vstack([pts, tiny])
    tris2 = np.vstack([tris, [[len(pts), len(pts) + 1, len(pts) + 2]]])
    nx, ny, x0, y0, dx, dy = 60, 40, -1.0, -1.0, 0.5, 0.2  # pixel 0.5 x 0.2: the tiny triangle holds one centre
    tri, w = raster.locate(pts2, tris2, x0, dx, nx, y0, dy, ny)
    xc = x0 + (np.arange(nx) + 0.5) * dx
    yc = y0 + (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(xc, yc)
    inside = ((X > 0) & (X < 10) & (Y > 0) & (Y < 6)).ravel()
    field = 3.0 + 2.0 * pts2[:, 0] - 0.7 * pts2[:, 1]
    val = (field[tris2[np.where(tri >= 0, tri, 0)]] * w).sum(1)
    exact = 3.0 + 2.0 * X.ravel() - 0.7 * Y.ravel()
    found = tri >= 0
    check(np.all(found[inside]), "every pixel centre inside the box has a triangle (%d pixels)" % inside.sum())
    err = np.abs(val - exact)[found].max()
    check(err < 1e-11, "a linear field is reproduced at every found pixel: largest error %.1e" % err)
    check((tri == len(tris2) - 1).sum() == 1, "the triangle smaller than a pixel is found by the one centre inside it")
    outside = ~inside & ~(tri == len(tris2) - 1)
    check(not np.any(found[outside]), "pixels outside the mesh are empty")

    print("2. basemap projection")
    region = os.path.join(TOOL, "regions", "socal_ucerf2d.json")
    mapdir = os.path.join(tmp, "map")
    os.makedirs(mapdir)
    for name in basemap.FILES.values():  # a stand-in for the Natural Earth files: one line each
        json.dump({"features": [{"geometry": {"type": "LineString",
                                               "coordinates": [[-121.0, 33.0], [-119.0, 35.0], [-117.0, 34.0]]}}]},
                  open(os.path.join(mapdir, name), "w"))
    bm = basemap.Basemap(region, mapdir)
    lat0, lon0 = 37.25, -119.25
    lat, lon = 34.05, -118.25
    x, y = bm.project(lon, lat)
    p0, p1, dl = np.radians(lat0), np.radians(lat), np.radians(lon - lon0)
    hav = np.sin((p1 - p0) / 2) ** 2 + np.cos(p0) * np.cos(p1) * np.sin(dl / 2) ** 2
    dist = 2 * 6371.0 * np.arcsin(np.sqrt(hav))
    bearing = np.arctan2(np.sin(dl) * np.cos(p1), np.cos(p0) * np.sin(p1) - np.sin(p0) * np.cos(p1) * np.cos(dl))
    check(abs(np.hypot(x, y) - dist) < 1e-6, "distance from the centre is the great-circle distance (%.3f km)" % dist)
    check(abs(np.arctan2(x, y) - bearing) < 1e-9, "bearing from the centre is kept")
    check(np.allclose(bm.project(lon0, lat0), 0.0), "the centre maps to (0, 0)")

    print("3. frames of a run")
    run, part = os.path.join(tmp, "run"), os.path.join(tmp, "run_stopped")
    make_run(run, True)
    make_run(part, False)
    s = cycle_run.Series(cycle_run.find(run, "_volume_snapshots.vtkhdf"), ["velocity"])
    fr = cycle_run.dynamic_frames(s, run, dt=4.0)
    check([p for p, _ in fr] == [1, 2, 3] and [len(k) for _, k in fr] == [15, 11, 13],
          "volume file, 2 s cadence, dt 4: phases %s with %s frames" % ([p for p, _ in fr], [len(k) for _, k in fr]))
    f = cycle_run.Series(cycle_run.find(run, "_fault_snapshots.vtkhdf"), ["slip_strike"])
    ff = cycle_run.dynamic_frames(f, run, dt=4.0)
    check([len(k) for _, k in ff] == [15, 11, 13], "fault file, 0.5 s cadence, dt 4: the same frame counts %s"
          % [len(k) for _, k in ff])
    sp = cycle_run.Series(cycle_run.find(part, "_volume_snapshots.vtkhdf"), ["velocity"])
    check([p for p, _ in cycle_run.dynamic_frames(sp, part, dt=4.0)] == [1, 2],
          "a run without DONE that ends inside phase 3: phases 1 and 2 only")
    check([p for p, _ in cycle_run.dynamic_frames(sp, part, dt=4.0, include_last=True)] == [1, 2, 3],
          "the same with include_last: phases 1 to 3")
    check(cycle_run.switch_times(run) == T_ON, "switch log: the three QD -> dynamic times")

    print("4. cycle_movies")
    o1, o3, om = os.path.join(tmp, "w1"), os.path.join(tmp, "w3"), os.path.join(tmp, "map_out")
    rc1, out1 = run_tool("cycle_movies.py", "--run", run, "--out", o1, "--segment-frames", "7", "--preset", "veryfast")
    rc3, out3 = run_tool("cycle_movies.py", "--run", run, "--out", o3, "--segment-frames", "7", "--preset", "veryfast",
                         "--workers", "3", "--max-mb", "0.02")
    check(rc1 == 0 and rc3 == 0, "both runs exit 0")
    exe = cycle_movies.ffmpeg_exe(None)
    for name in ("particle_velocity", "fault_slip"):
        a, b = os.path.join(o1, name + ".mp4"), os.path.join(o3, name + ".mp4")
        ok = os.path.exists(a) and os.path.exists(b)
        n = cycle_movies.count_frames(exe, a) if ok else -1
        check(n == 39, "%s: %d frames (15 + 11 + 13 = 39)" % (name, n))
        check(ok and open(a, "rb").read() == open(b, "rb").read(), "%s: 1 worker and 3 workers give the same file" % name)
        parts = sorted(p for p in os.listdir(o3) if p.startswith(name + "_part"))
        total = sum(cycle_movies.count_frames(exe, os.path.join(o3, p)) for p in parts)
        check(len(parts) >= 2 and total == 39, "%s: %d parts of --max-mb hold %d frames" % (name, len(parts), total))
    check(not os.path.exists(os.path.join(o1, "particle_velocity_segments")), "the segment folder is removed")
    rcm, outm = run_tool("cycle_movies.py", "--run", run, "--out", om, "--segment-frames", "20", "--preset", "veryfast",
                         "--region", region, "--map-data", mapdir, "--workers", "2", "--last-phase", "1")
    okm = rcm == 0 and all(os.path.exists(os.path.join(om, n + "_map.mp4")) for n in ("particle_velocity", "fault_slip"))
    check(okm and cycle_movies.count_frames(exe, os.path.join(om, "fault_slip_map.mp4")) == 15,
          "line map version, phase 1 only: both videos written, 15 frames")
    for tag, base, extra in (("plain", o1, ["--segment-frames", "7"]),
                             ("line map", om, ["--segment-frames", "20", "--region", region, "--map-data", mapdir,
                                               "--last-phase", "1"])):
        of = base + "_full"
        rcf, outf = run_tool("cycle_movies.py", "--run", run, "--out", of, "--preset", "veryfast", "--full-redraw",
                             *extra)
        same = rcf == 0 and all(open(os.path.join(base, n), "rb").read() == open(os.path.join(of, n), "rb").read()
                                for n in sorted(os.listdir(of)) if n.endswith(".mp4") and "_part" not in n)
        check(same, "%s: drawing only the changing layers gives the same files as --full-redraw" % tag)
    rcx, outx = run_tool("cycle_movies.py", "--run", run, "--out", run, quiet=True)
    check(rcx != 0 and "must not be the run folder" in outx, "--out = the run folder is refused")

    print("5. magnitude_catalog")
    rows, lengths = magnitude_catalog.events(run, MU, W_M)
    check(len(rows) == 3 and abs(lengths[0] - FAULT_KM) < 1e-9 and abs(lengths[1] - 5.0) < 1e-9,
          "3 events; fault lengths %.1f and %.1f km" % tuple(lengths))
    for r, d in zip(rows, SLIP):
        mw = (2.0 / 3.0) * (np.log10(MU * W_M * FAULT_KM * 1e3 * d) - 9.1)
        check(abs(r["mw"] - mw) < 1e-9 and r["fault"] == 0 and abs(r["start_gap_s"] - 30.0) < 1e-6,
              "event %d: Mw %.4f = (2/3)(log10(mu W L D) - 9.1) with D = %g m; fault 1; start gap 30 s"
              % (r["row"], r["mw"], d))
    rcg, outg = run_tool("magnitude_catalog.py", "--run", run, "--out", os.path.join(tmp, "mag"), "--mu", str(MU),
                         "--width-m", str(W_M), "--fault-names", "main,side")
    check(rcg == 0 and os.path.exists(os.path.join(tmp, "mag", "magnitude_catalog.png"))
          and "main: 40.0 km, 3 events" in outg, "the script writes the CSV and the figure")
    print("%d check(s) failed; files kept in %s" % (failed, tmp) if failed else "all checks passed")
    if not failed:
        shutil.rmtree(tmp)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
