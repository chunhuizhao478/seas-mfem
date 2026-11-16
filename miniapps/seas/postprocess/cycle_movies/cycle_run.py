# Readers for the output of seas_fd_cycle_driver, shared by the movie and the
# magnitude scripts of this folder (Chunhui Zhao, 2026-11-16: "let's develop
# this, later this could be a standardlized video generation process").
#
# A run folder holds (cycle/fd_cycle_output.hpp, cycle/snapshots.hpp):
#   <prefix>_fault_snapshots.vtkhdf    fault fields, one line cell per fault
#                                      face on a 2-D mesh
#   <prefix>_volume_snapshots.vtkhdf   displacement and velocity, one triangle
#                                      per element (per sub-triangle at p > 1)
#   <prefix>_switch.txt                time, Vmax, from, to, reason
#   <prefix>_catalog.txt               event, t_start, t_end, peak_Vmax
#   DONE                               written by the job script at exit 0
# A snapshot file is a VTKHDF series: /VTKHDF/Points, Connectivity, Offsets;
# step k of field F is rows Steps/PointDataOffsets/F[k] .. + number of points
# of PointData/F; Steps/Values[k] is its time; /fd_cycle/regime[k] is 0 in
# the quasi-dynamic regime and 1 in the dynamic one, /fd_cycle/phase[k] the
# number of the dynamic phase.
#
# Everything is opened read-only with HDF5 file locking off, so a run that is
# still being written can be read; steps whose rows are not all on disk yet
# are left out.
import glob
import os

import h5py
import numpy as np

YEAR_S = 365.25 * 86400.0  # the driver's year


def find(run, suffix):
    """The one file <run>/*<suffix>."""
    hits = sorted(glob.glob(os.path.join(run, "*" + suffix)))
    if len(hits) != 1:
        raise SystemExit("expected one *%s in %s, found %d" % (suffix, run, len(hits)))
    return hits[0]


class Series:
    """One snapshot file.  Points in km (x, y)."""

    def __init__(self, path, fields):
        self.path = path
        self.f = h5py.File(path, "r", locking=False)
        w = self.f["VTKHDF"]
        self.w = w
        self.points = w["Points"][:, :2] / 1e3
        self.npts = len(self.points)
        self.conn = w["Connectivity"][:]
        self.offsets = w["Offsets"][:]
        self.times = w["Steps/Values"][:]
        self.regime = self.f["fd_cycle/regime"][:]
        self.phase = self.f["fd_cycle/phase"][:]
        self.off = {}
        n = min(len(self.times), len(self.regime), len(self.phase))
        for name in fields:
            if "PointData/" + name not in w:
                raise SystemExit("%s has no field %s" % (path, name))
            o = w["Steps/PointDataOffsets/" + name][:]
            rows = w["PointData/" + name].shape[0]
            n = min(n, len(o), int(np.searchsorted(o + self.npts, rows, side="right")))
            self.off[name] = o
        self.nsteps = n  # steps that are complete on disk for every field asked for

    def field(self, name, k):
        """Field `name` at step k: (npts,) or (npts, components)."""
        o = int(self.off[name][k])
        return self.w["PointData/" + name][o:o + self.npts]

    def cells(self, nodes_per_cell):
        """Connectivity as (ncells, nodes_per_cell); every cell must have that many nodes."""
        sizes = np.diff(self.offsets)
        if len(sizes) == 0 or not np.all(sizes == nodes_per_cell):
            raise SystemExit("%s: the cells are not all %d-node cells (triangles in the volume, "
                             "lines on the fault of a 2-D mesh)" % (self.path, nodes_per_cell))
        return self.conn.reshape(-1, nodes_per_cell)


def switch_times(run):
    """Times (s) of the QD -> dynamic switches, in order: entry p - 1 starts dynamic phase p."""
    out = []
    for line in open(find(run, "_switch.txt")):
        a = line.split()
        if line.startswith("#") or len(a) < 4:
            continue
        if a[2] == "QD" and a[3] == "dynamic":
            out.append(float(a[0]))
    return out


def catalog(run):
    """Rows event, t_start (s), t_end (s), peak_Vmax (m/s)."""
    return np.loadtxt(find(run, "_catalog.txt"), ndmin=2)


def dynamic_frames(series, run, first=1, last=None, dt=None, include_last=None):
    """[(phase, [step, ...]), ...]: the dynamic snapshots of every finished phase.

    dt (simulated seconds per frame): every n-th snapshot is taken, n = dt
    over the file's own dynamic cadence (the median step of its dynamic
    snapshots), at least 1.  None takes every snapshot.

    A phase is finished when a later step in the file is quasi-dynamic or
    belongs to a later phase.  The last phase of the file has no later step
    when the run stopped inside it or at its end: it is taken when the run
    folder has DONE (the job script's mark of exit 0), or with include_last.
    """
    n = series.nsteps
    ph, rg, tv = series.phase[:n], series.regime[:n], series.times[:n]
    if include_last is None:
        include_last = os.path.exists(os.path.join(run, "DONE"))
    stride = 1
    if dt is not None:
        d = np.diff(tv)[(rg[1:] == 1) & (rg[:-1] == 1) & (ph[1:] == ph[:-1])]
        if len(d):
            stride = max(1, int(round(dt / float(np.median(d)))))
    out = []
    for p in sorted(set(ph[rg == 1].tolist())):
        if p < first or (last is not None and p > last):
            continue
        ks = np.where((ph == p) & (rg == 1))[0]
        after = np.arange(ks.max() + 1, n)
        done = np.any(rg[after] == 0) or np.any(ph[after] > p)
        if done or include_last:
            out.append((int(p), ks[::stride].tolist()))
    return out
