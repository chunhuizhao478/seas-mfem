#!/usr/bin/env python3
"""The standard videos of a fully dynamic cycle run (seas_fd_cycle_driver, 2-D map view).

    python3 cycle_movies.py --run RUN_FOLDER --out OUT_FOLDER --workers 32

writes, from the run's snapshot files (cycle_run.py),

    particle_velocity.mp4   in-plane particle speed |v| over the output box,
                            linear from 0, the fault traces on top
    fault_slip.mp4          slip rate (linear from 0) beside cumulative slip
                            on the fault traces, map view of the whole system

for the dynamic phases of every finished event, back to back, one frame per
`--dt` seconds of simulated time.  With --region and --map-data the same two
videos are drawn on a line map (basemap.py) and get the suffix _map.

Speed: the frames are cut into segments; `--workers` processes each render
their segments to H.264 files, which are then joined without re-encoding.
Every segment starts with a key frame and all use the same encoder settings,
so the joined video has the same frames whatever the number of workers.
On a login node use few workers; on a compute node one per core is fine
(jobs/postprocess/ls6_cycle_movies.sbatch).

--max-mb N also writes the video in parts of at most N MB (cut between
segments), for places with a file size limit.

Chunhui Zhao, 2026-11-16: "we need to find a way to speed this up"; "let's
develop this, later this could be a standardlized video generation process".
"""
import argparse
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cycle_run  # noqa: E402
import raster  # noqa: E402

G = {}  # what the worker processes inherit from the parent (fork): frames, grids, options


def ffmpeg_exe(given):
    if given:
        return given
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        exe = shutil.which("ffmpeg")
        if not exe:
            raise SystemExit("no ffmpeg: pip install imageio-ffmpeg, or pass --ffmpeg")
        return exe


def label(phase, t, t_on):
    return "dynamic phase %d   %s yr   t = %3.0f s" % (phase, format(int(round(t_on / cycle_run.YEAR_S)), ","),
                                                         t - t_on)


def fault_segments(run):
    """Fault cells of the fault snapshot file: end points (n, 2, 2) in km and node ids (n, 2)."""
    s = cycle_run.Series(cycle_run.find(run, "_fault_snapshots.vtkhdf"), [])
    c = s.cells(2)
    return s.points[c], c


class Blit:
    """Draw the still layers of a figure once, then per frame only the layers that change.

    changing: [(axes, zmin), ...]: in each of these axes every artist with
    zorder >= zmin is drawn again on every frame (the artist whose data
    changes and all that lie on top of it), in the order a full draw uses.
    The frames are the ones a full draw gives (checked by the test of this
    folder); --full-redraw draws everything on every frame.
    """

    def __init__(self, fig, changing, full_redraw):
        self.fig, self.full, self.bg = fig, full_redraw, None
        self.arts = []
        for ax, zmin in changing:
            arts = [c for c in ax.get_children() if c is not ax.patch and c.get_zorder() >= zmin]
            self.arts += [(ax, c) for c in sorted(arts, key=lambda c: c.get_zorder())]

    def frame(self):
        canvas = self.fig.canvas
        if self.full:
            canvas.draw()
        else:
            if self.bg is None:
                for _, c in self.arts:
                    c.set_animated(True)  # left out of the full draw below
                canvas.draw()
                self.bg = canvas.copy_from_bbox(self.fig.bbox)
            canvas.restore_region(self.bg)
            for ax, c in self.arts:
                ax.draw_artist(c)
        return np.asarray(canvas.buffer_rgba())[:, :, :3]


class VolumeMovie:
    """|v| on a pixel grid.  setup() runs once in the parent, open() once in each worker."""
    name = "particle_velocity"
    size = (960, 800)

    @staticmethod
    def setup(a, bm):
        s = cycle_run.Series(cycle_run.find(a.run, "_volume_snapshots.vtkhdf"), ["velocity"])
        tris = s.cells(3)
        if bm:
            x0, x1, y0, y1 = bm.window
        else:
            lo, hi = s.points.min(0), s.points.max(0)
            m = 0.02 * max(hi - lo)
            x0, x1, y0, y1 = lo[0] - m, hi[0] + m, lo[1] - m, hi[1] + m
        fw, fh = VolumeMovie.size
        nx = fw - 100
        ny = int(round(nx * (y1 - y0) / (x1 - x0)))
        if ny > fh - 10:
            ny = fh - 10
            nx = int(round(ny * (x1 - x0) / (y1 - y0)))
        tri, w = raster.locate(s.points, tris, x0, (x1 - x0) / nx, nx, y0, (y1 - y0) / ny, ny)
        G["vol"] = dict(window=(x0, x1, y0, y1), nx=nx, ny=ny, inside=(tri >= 0).reshape(ny, nx),
                        nodes=tris[np.where(tri >= 0, tri, 0)], w=w, segs=fault_segments(a.run)[0])
        frames = cycle_run.dynamic_frames(s, a.run, a.first_phase, a.last_phase, a.dt, a.include_last)
        return frames

    def __init__(self, a, bm):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from matplotlib.colors import Normalize
        v = G["vol"]
        self.v = v
        self.s = cycle_run.Series(cycle_run.find(a.run, "_volume_snapshots.vtkhdf"), ["velocity"])
        fw, fh = self.size
        x0, x1, y0, y1 = v["window"]
        fig = plt.figure(figsize=(fw / 100, fh / 100), dpi=100, facecolor="white")
        ax = fig.add_axes([8 / fw, (fh - v["ny"]) / 2 / fh, v["nx"] / fw, v["ny"] / fh])
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("0.6")
        if bm:
            bm.draw_under(ax)
        cm = matplotlib.colormaps["magma_r"].copy()
        cm.set_bad((0, 0, 0, 0))
        self.im = ax.imshow(np.ma.masked_all((v["ny"], v["nx"])), origin="lower", extent=(x0, x1, y0, y1),
                            cmap=cm, norm=Normalize(0, a.v_max), interpolation="nearest", aspect="auto",
                            zorder=1)
        if bm:
            bm.draw_over(ax, zorder=3)
        ax.add_collection(LineCollection(v["segs"], linewidths=0.7, colors="k", zorder=4))
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        self.lab = ax.text(0.012, 0.012, "", transform=ax.transAxes, fontsize=9, va="bottom", zorder=10,
                           bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85))
        cb = fig.colorbar(self.im, cax=fig.add_axes([(v["nx"] + 24) / fw, 0.2, 14 / fw, 0.6]), extend="max")
        cb.set_label("|v| (m/s), linear from 0", fontsize=9)
        cb.ax.tick_params(labelsize=8)
        self.fig = fig
        self.blit = Blit(fig, [(ax, 1)], a.full_redraw)

    def frame(self, phase, k, t_on):
        v = self.v
        u = self.s.field("velocity", k)[:, :2]
        ux = (u[v["nodes"], 0] * v["w"]).sum(1)
        uy = (u[v["nodes"], 1] * v["w"]).sum(1)
        self.im.set_data(np.ma.masked_array(np.hypot(ux, uy).reshape(v["ny"], v["nx"]), ~v["inside"]))
        self.lab.set_text(label(phase, self.s.times[k], t_on))
        return self.blit.frame()


class FaultMovie:
    """Slip rate | cumulative slip on the fault traces."""
    name = "fault_slip"
    size = (1280, 640)
    fields = ["slip_rate_strike", "slip_rate_dip", "slip_strike", "slip_dip"]

    @staticmethod
    def setup(a, bm):
        s = cycle_run.Series(cycle_run.find(a.run, "_fault_snapshots.vtkhdf"), FaultMovie.fields)
        frames = cycle_run.dynamic_frames(s, a.run, a.first_phase, a.last_phase, a.dt, a.include_last)
        slip_max = a.slip_max
        if slip_max is None and frames:  # the largest slip at the last frame, rounded up to 1 m
            k = frames[-1][1][-1]
            slip_max = float(np.ceil(np.hypot(s.field("slip_strike", k), s.field("slip_dip", k)).max()))
        segs, cells = fault_segments(a.run)
        if bm:
            window = bm.window
        else:
            lo, hi = s.points.min(0), s.points.max(0)
            m = 0.08 * max(hi - lo)
            window = (lo[0] - m, hi[0] + m, lo[1] - m, hi[1] + m)
        G["fault"] = dict(segs=segs, cells=cells, slip_max=slip_max or 1.0, window=window)
        return frames

    def __init__(self, a, bm):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from matplotlib.colors import Normalize
        g = G["fault"]
        self.g = g
        self.s = cycle_run.Series(cycle_run.find(a.run, "_fault_snapshots.vtkhdf"), self.fields)
        fw, fh = self.size
        x0, x1, y0, y1 = g["window"]
        fig = plt.figure(figsize=(fw / 100, fh / 100), dpi=100, facecolor="white")
        axs = [fig.add_axes([0.02, 0.16, 0.47, 0.8]), fig.add_axes([0.51, 0.16, 0.47, 0.8])]
        cms = [matplotlib.colormaps["magma_r"], matplotlib.colormaps["viridis"]]
        norms = [Normalize(0, a.rate_max), Normalize(0, g["slip_max"])]
        self.lcs = []
        for ax, cm, nm in zip(axs, cms, norms):
            if bm:
                bm.draw_under(ax)
                bm.draw_over(ax, zorder=2, fontsize=6, city_zorder=8)  # city names above the fault traces
            ax.add_collection(LineCollection(g["segs"], linewidths=5.0, colors="0.25", capstyle="round", zorder=6))
            lc = LineCollection(g["segs"], linewidths=3.6, cmap=cm, norm=nm, capstyle="round", zorder=7)
            lc.set_array(np.zeros(len(g["segs"])))
            ax.add_collection(lc)
            self.lcs.append(lc)
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color("0.6")
        c1 = fig.colorbar(self.lcs[0], cax=fig.add_axes([0.08, 0.08, 0.35, 0.025]), orientation="horizontal",
                          extend="max")
        c1.set_label("slip rate (m/s), linear 0 to %g" % a.rate_max, fontsize=9)
        c2 = fig.colorbar(self.lcs[1], cax=fig.add_axes([0.57, 0.08, 0.35, 0.025]), orientation="horizontal")
        c2.set_label("cumulative slip (m), linear 0 to %g" % g["slip_max"], fontsize=9)
        self.lab = axs[0].text(0.015, 0.015, "", transform=axs[0].transAxes, fontsize=10, zorder=10,
                               bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85))
        self.fig = fig
        self.blit = Blit(fig, [(axs[0], 7), (axs[1], 7)], a.full_redraw)

    def frame(self, phase, k, t_on):
        s, c = self.s, self.g["cells"]
        rate = np.hypot(s.field("slip_rate_strike", k), s.field("slip_rate_dip", k))
        slip = np.hypot(s.field("slip_strike", k), s.field("slip_dip", k))
        self.lcs[0].set_array(rate[c].mean(1))
        self.lcs[1].set_array(slip[c].mean(1))
        self.lab.set_text(label(phase, s.times[k], t_on))
        return self.blit.frame()


MOVIES = {"volume": VolumeMovie, "fault": FaultMovie}
_open = {}  # per process: the movie object of each kind, made on first use


def _movie(kind):
    if kind not in _open:
        _open[kind] = MOVIES[kind](G["args"], G["basemap"])
    return _open[kind]


def render_segment(job):
    """Worker: frames [i0, i1) of the flat frame list -> one H.264 file.  Returns (index, frames, seconds)."""
    kind, index, i0, i1, path = job
    a = G["args"]
    t0 = time.time()
    m = _movie(kind)
    fw, fh = m.size
    ff = subprocess.Popen(
        [G["ffmpeg"], "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", "%dx%d" % (fw, fh),
         "-r", str(a.fps), "-i", "-", "-c:v", "libx264", "-preset", a.preset, "-crf", str(a.crf),
         "-pix_fmt", "yuv420p", "-threads", "1", path], stdin=subprocess.PIPE)
    for phase, k in G["flat"][kind][i0:i1]:
        ff.stdin.write(np.ascontiguousarray(m.frame(phase, k, G["t_on"][phase - 1])).tobytes())
    ff.stdin.close()
    if ff.wait() != 0:
        raise RuntimeError("ffmpeg failed on " + path)
    return index, i1 - i0, time.time() - t0


def count_frames(ffmpeg, path):
    """Frames ffmpeg decodes from the file (every frame is read)."""
    r = subprocess.run([ffmpeg, "-loglevel", "error", "-stats", "-i", path, "-map", "0:v:0", "-f", "null", "-"],
                       stderr=subprocess.PIPE, universal_newlines=True)
    if r.returncode != 0:
        raise SystemExit("ffmpeg cannot decode %s: %s" % (path, r.stderr[-400:]))
    n = re.findall(r"frame=\s*(\d+)", r.stderr)
    return int(n[-1]) if n else -1


def join(ffmpeg, parts, out):
    lst = out + ".list.txt"
    with open(lst, "w") as f:
        for p in parts:
            f.write("file '%s'\n" % os.path.abspath(p))
    subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", lst, "-c", "copy",
                    "-movflags", "+faststart", out], check=True)
    os.remove(lst)


def make(kind, a, bm):
    cls = MOVIES[kind]
    t0 = time.time()
    frames = cls.setup(a, bm)
    flat = [(p, k) for p, ks in frames for k in ks]
    if not flat:
        print("%s: no finished dynamic phase in the range asked for; nothing written" % cls.name)
        return
    missing = [p for p, _ in frames if p > len(G["t_on"])]
    if missing:
        raise SystemExit("the switch log has %d QD -> dynamic rows but the snapshots have phase %d"
                         % (len(G["t_on"]), missing[0]))
    G["flat"][kind] = flat
    name = cls.name + ("_map" if bm else "")
    if a.frame_png is not None:
        import matplotlib.pyplot as plt
        phase, k = flat[a.frame_png]
        png = os.path.join(a.out, "%s_frame%d.png" % (name, a.frame_png))
        plt.imsave(png, _movie(kind).frame(phase, k, G["t_on"][phase - 1]))
        print(png)
        return
    tmp = os.path.join(a.out, name + "_segments")
    os.makedirs(tmp, exist_ok=True)
    nseg = max(1, int(np.ceil(len(flat) / a.segment_frames)))
    cuts = np.linspace(0, len(flat), nseg + 1).round().astype(int)
    jobs = [(kind, i, int(cuts[i]), int(cuts[i + 1]), os.path.join(tmp, "seg%05d.mp4" % i)) for i in range(nseg)]
    print("%s: %d phases (%d to %d), %d frames in %d segments; setup %.0f s"
          % (name, len(frames), frames[0][0], frames[-1][0], len(flat), nseg, time.time() - t0), flush=True)
    t1 = time.time()
    done = 0
    # the pool starts here, after setup() filled G: the workers inherit the grids (fork)
    pool = multiprocessing.get_context("fork").Pool(a.workers) if a.workers > 1 else None
    for index, n, sec in (pool.imap_unordered(render_segment, jobs) if pool else map(render_segment, jobs)):
        done += n
        if a.verbose:
            print("  segment %d: %d frames in %.1f s (%d of %d frames done)" % (index, n, sec, done, len(flat)),
                  flush=True)
    if pool:
        pool.close()
        pool.join()
    segs = [j[4] for j in jobs]
    out = os.path.join(a.out, name + ".mp4")
    join(G["ffmpeg"], segs, out)
    n = count_frames(G["ffmpeg"], out)
    if n != len(flat):
        raise SystemExit("%s has %d frames, expected %d" % (out, n, len(flat)))
    wall = time.time() - t1
    print("%s: %d frames, %.1f MB, rendered in %.0f s (%.3f s per frame wall, %d workers)"
          % (out, n, os.path.getsize(out) / 1e6, wall, wall / len(flat), a.workers), flush=True)
    if a.max_mb and os.path.getsize(out) > a.max_mb * 1e6:
        groups, size = [[]], 0.0
        for s in segs:
            b = os.path.getsize(s)
            if groups[-1] and size + b > 0.98 * a.max_mb * 1e6:
                groups.append([])
                size = 0.0
            groups[-1].append(s)
            size += b
        for i, g in enumerate(groups):
            part = os.path.join(a.out, "%s_part%d.mp4" % (name, i + 1))
            join(G["ffmpeg"], g, part)
            j0, j1 = jobs[segs.index(g[0])][2], jobs[segs.index(g[-1])][3]
            print("  %s: frames %d to %d (phase %d to %d), %.1f MB"
                  % (part, j0, j1 - 1, flat[j0][0], flat[j1 - 1][0], os.path.getsize(part) / 1e6))
    if not a.keep_segments:
        shutil.rmtree(tmp)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True, help="run folder of seas_fd_cycle_driver")
    p.add_argument("--out", required=True, help="folder for the videos (made if missing; never the run folder)")
    p.add_argument("--movies", default="volume,fault", help="which videos: volume, fault or both (default)")
    p.add_argument("--workers", type=int, default=1, help="render processes (default 1)")
    p.add_argument("--dt", type=float, default=4.0,
                   help="simulated seconds per frame (default 4): every n-th snapshot of each file, "
                        "n from the file's own dynamic cadence")
    p.add_argument("--fps", type=int, default=24, help="frames per second (default 24)")
    p.add_argument("--first-phase", type=int, default=1)
    p.add_argument("--last-phase", type=int, default=None)
    p.add_argument("--include-last", action="store_true", default=None,
                   help="take the last phase of the file although the run folder has no DONE")
    p.add_argument("--v-max", type=float, default=2.0, help="top of the |v| colour scale, m/s (default 2)")
    p.add_argument("--rate-max", type=float, default=10.0, help="top of the slip rate scale, m/s (default 10)")
    p.add_argument("--slip-max", type=float, default=None,
                   help="top of the slip scale, m (default: the largest slip at the last frame, rounded up)")
    p.add_argument("--region", help="region file for the line map (regions/*.json); with --map-data")
    p.add_argument("--map-data", help="folder with the Natural Earth files (fetch_natural_earth.sh)")
    p.add_argument("--crf", type=int, default=20, help="x264 quality (default 20; lower is better)")
    p.add_argument("--preset", default="slow", help="x264 preset (default slow)")
    p.add_argument("--segment-frames", type=int, default=120,
                   help="frames per segment (default 120); a segment is the unit of work of a worker")
    p.add_argument("--max-mb", type=float, default=None, help="also write parts of at most this many MB")
    p.add_argument("--ffmpeg", help="ffmpeg executable (default: imageio-ffmpeg's, else the one on PATH)")
    p.add_argument("--keep-segments", action="store_true")
    p.add_argument("--full-redraw", action="store_true",
                   help="draw every layer on every frame (slower; the frames are the same)")
    p.add_argument("--frame-png", type=int, default=None,
                   help="write only frame N of each video as <name>_frameN.png (a look before a long render)")
    p.add_argument("--verbose", action="store_true")
    a = p.parse_args()
    if bool(a.region) != bool(a.map_data):
        raise SystemExit("--region and --map-data go together")
    if os.path.realpath(a.out) == os.path.realpath(a.run):
        raise SystemExit("--out must not be the run folder")
    kinds = [k.strip() for k in a.movies.split(",") if k.strip()]
    for k in kinds:
        if k not in MOVIES:
            raise SystemExit("--movies: %s is not one of %s" % (k, ", ".join(MOVIES)))
    os.makedirs(a.out, exist_ok=True)
    bm = None
    if a.region:
        import basemap
        bm = basemap.Basemap(a.region, a.map_data)
    G.update(args=a, basemap=bm, ffmpeg=ffmpeg_exe(a.ffmpeg), t_on=cycle_run.switch_times(a.run), flat={})
    for kind in kinds:
        make(kind, a, bm)


if __name__ == "__main__":
    main()
