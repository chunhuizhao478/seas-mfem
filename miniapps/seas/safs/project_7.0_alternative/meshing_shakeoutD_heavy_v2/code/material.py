#!/usr/bin/env python3
"""material.py -- the Vs the resolution gate is judged against.

Two cubes describe the SAME model (MUSCAL), but not equally well:

  deck    `safs_material_cvm.nc` -- what SeisSol actually reads through ASAGI.
          UTM 11N, 1500 m lateral, 250 m vertical, built from 37 UCVM
          web-viewer CSV slices at ~6 km.
  native  `MUSCAL.nc` (Zenodo 19243477) -- the model as published.
          lon/lat 0.01 deg, 210 depth levels, 50 m steps through the top 3 km.

The delivery path costs real accuracy exactly where this project's gate failures
live: measured over the cube footprint, the deck cube runs **+25.7 % fast at
50 m depth and +15.7 % at 100 m**, and agrees to <1 % below 250 m
(`Downloads/muscal_nc/compare_nc/README.md`).  A gate scored on the deck cube is
therefore too LENIENT in the shallow band -- it credits the mesh with material
that is faster than the real model's.

Neither cube bounds the other pointwise (the deck is biased fast on average but
slower at some nodes), so the gate is scored on **min(deck, native)**: the mesh
must resolve both the model the solver reads and the model the Earth has.  Where
MUSCAL is NaN -- outside its published footprint -- the deck value stands alone.

Depth convention follows the production cube's own: depth = -z, clamped at 0.
"""

import numpy as np

_UTM11N = "EPSG:32611"


class Material:
    def __init__(self, deck_nc, muscal_nc=None, box=None, margin_deg=0.1,
                 source="min", verbose=True):
        # source: 'min'    -- min(deck, native): satisfies both readings
        #         'muscal' -- the published model IS the measurement; the deck
        #                     value is used only where MUSCAL is NaN (outside
        #                     its footprint), because nothing else is available
        #         'deck'   -- legacy
        self.source = source
        import netCDF4 as ncdf
        d = ncdf.Dataset(str(deck_nc))
        self.X = d["x"][:].data; self.Y = d["y"][:].data; self.Z = d["z"][:].data
        D = d["data"][:]
        self.Vd = np.sqrt(np.maximum(D["mu"], 0) / D["rho"]).astype(np.float32)
        del D
        dz = np.diff(self.Z)
        if not np.allclose(dz, dz[0]):
            raise RuntimeError("deck cube z axis is not uniform")
        self.nz = len(self.Z)
        self.M = None
        if verbose:
            print(f"[deck] {self.Vd.shape} (z,y,x)  z {self.Z[0]:.0f}..{self.Z[-1]:.0f} "
                  f"dz {dz[0]:.0f}  Vs {np.nanmin(self.Vd):.1f}..{np.nanmax(self.Vd):.1f}")
        if muscal_nc is None:
            return

        from pyproj import Transformer
        self.fwd = Transformer.from_crs(_UTM11N, "EPSG:4326", always_xy=True)
        m = ncdf.Dataset(str(muscal_nc))
        lon = m["longitude"][:].data.astype(np.float64)
        lat = m["latitude"][:].data.astype(np.float64)
        self.dep = m["depth"][:].data.astype(np.float64)
        if box is None:
            raise ValueError("MUSCAL needs the UTM box to subset on")
        # project a DENSE sample of the box boundary -- projected edges bulge, so
        # the four corners alone under-cover the footprint.
        x0, x1, y0, y1 = box
        t = np.linspace(0, 1, 400)
        bx = np.concatenate([x0 + (x1-x0)*t, x0 + (x1-x0)*t, np.full(400, x0), np.full(400, x1)])
        by = np.concatenate([np.full(400, y0), np.full(400, y1), y0 + (y1-y0)*t, y0 + (y1-y0)*t])
        blon, blat = self.fwd.transform(bx, by)
        i0 = max(np.searchsorted(lon, blon.min() - margin_deg) - 1, 0)
        i1 = min(np.searchsorted(lon, blon.max() + margin_deg) + 1, len(lon))
        j0 = max(np.searchsorted(lat, blat.min() - margin_deg) - 1, 0)
        j1 = min(np.searchsorted(lat, blat.max() + margin_deg) + 1, len(lat))
        self.mlon = lon[i0:i1]; self.mlat = lat[j0:j1]
        self.Vm = np.asarray(m["vs"][:, j0:j1, i0:i1], np.float32)
        if verbose:
            cov = 100.0 * np.isfinite(self.Vm).mean()
            print(f"[muscal] subset {self.Vm.shape} (depth,lat,lon)  "
                  f"lon {self.mlon[0]:.2f}..{self.mlon[-1]:.2f}  "
                  f"lat {self.mlat[0]:.2f}..{self.mlat[-1]:.2f}  "
                  f"depth {self.dep[0]:.0f}..{self.dep[-1]:.0f} ({len(self.dep)} levels)  "
                  f"finite {cov:.2f} %  Vs {np.nanmin(self.Vm):.1f}..{np.nanmax(self.Vm):.1f}")
        self.M = True

    # ---- deck -------------------------------------------------------------
    def _dij(self, p):
        i = np.clip(np.rint((p[:, 0] - self.X[0]) / (self.X[1] - self.X[0])).astype(np.int32), 0, len(self.X) - 1)
        j = np.clip(np.rint((p[:, 1] - self.Y[0]) / (self.Y[1] - self.Y[0])).astype(np.int32), 0, len(self.Y) - 1)
        return i, j

    def _dk(self, z):
        return np.clip(np.rint((z - self.Z[0]) / (self.Z[1] - self.Z[0])).astype(np.int32), 0, self.nz - 1)

    # ---- muscal -----------------------------------------------------------
    @staticmethod
    def _near(axis, v):
        k = np.clip(np.searchsorted(axis, v), 0, len(axis) - 1)
        km = np.maximum(k - 1, 0)
        return np.where(np.abs(axis[km] - v) < np.abs(axis[k] - v), km, k).astype(np.int32)

    def _mij(self, p):
        lo, la = self.fwd.transform(p[:, 0], p[:, 1])
        return self._near(self.mlon, lo), self._near(self.mlat, la)

    def _mk(self, z):
        return self._near(self.dep, np.maximum(-z, 0.0))

    # ---- public -----------------------------------------------------------
    @staticmethod
    def _lat_offsets(axis, c, r):
        """Which neighbouring COLUMN indices a point within radius r can reach.

        Nearest-grid picks column i' whenever |axis[i'] - c| is minimal, so a
        descendant barycentre inside the ball of radius r can land on any column
        whose Voronoi cell the ball touches: |axis[i'] - c| <= r + step/2.
        """
        step = abs(axis[1] - axis[0])
        return int(np.floor((r.max() + 0.5 * step) / step)) if len(axis) > 1 else 0

    def at(self, p):
        """Nearest-grid Vs at points p (n,3) in UTM 11N, per `source`."""
        i, j = self._dij(p)
        out = self.Vd[self._dk(p[:, 2]), j, i]
        if self.M is None or self.source == "deck":
            return out
        mi, mj = self._mij(p)
        vm = self.Vm[self._mk(p[:, 2]), mj, mi]
        if self.source == "muscal":
            return np.where(np.isfinite(vm), vm, out)
        return np.where(np.isfinite(vm), np.minimum(out, vm), out)

    def pooled(self, v, dx=None, frac=0.25):
        """min Vs over the z window a DESCENDANT's barycentre can reach.

        The gate is scored at the barycentre, so what must be bounded is where a
        child's barycentre can move to, not the cell's whole extent.  Bisecting
        an edge (p,q) shifts the centroid by (p-q)/8, and the displacements of
        repeated splits sum to a fraction of the cell's diameter: the barycentre
        wander radius is ~frac*dx.  This is the same "min-pooled Vs over the
        wander radius" rule this project already uses for coarsening budgets.

        Pooling over the FULL vertical extent instead (the obvious reading) is a
        trap once the material is sampled finely in z: every cell touching the
        free surface then inherits the 116 m/s surface value and is required to
        reach ~175 m, and its children inherit the same demand even though they
        already pass the gate -- measured as a steady-state runaway, +505k tets
        over 50 rounds with the failure count flat.

        Pooled in z ONLY: both cubes are far coarser laterally than vertically.
        """
        b = v.mean(1)
        if dx is None:
            zlo = v[:, :, 2].min(1); zhi = v[:, :, 2].max(1)
        else:
            r = frac * dx
            zlo = np.maximum(b[:, 2] - r, v[:, :, 2].min(1))
            zhi = np.minimum(b[:, 2] + r, v[:, :, 2].max(1))
        i, j = self._dij(b)
        kb = self._dk(b[:, 2])
        k0 = np.minimum(self._dk(zlo), kb); k1 = np.maximum(self._dk(zhi), kb)
        out = self.Vd[k0, j, i].copy()
        # NOTE `np.minimum(out[m], x, out=out[m])` is a silent no-op: out[m] is a
        # fancy-index COPY on both sides.  Index explicitly.
        for o in range(1, int((k1 - k0).max()) + 1 if len(k0) else 0):
            sel = np.flatnonzero((k0 + o) <= k1)
            if not len(sel):
                continue
            out[sel] = np.minimum(out[sel], self.Vd[k0[sel] + o, j[sel], i[sel]])
        # ---- LATERAL wander -------------------------------------------------
        # Pooling in z alone leaves the barycentre free to cross a LATERAL bin
        # edge into a slower column, which is exactly what made bisection
        # oscillate on flat near-surface cells: MUSCAL's lateral step is ~1 km
        # and surface Vs contrast between adjacent columns is severe (116 vs
        # 800+ m/s at basin edges).  The wander radius is small (0.25*dx), so
        # this is at most the ONE adjacent column -- not the full 3-D min over
        # the cell, which is the documented over-conservative trap.
        if dx is not None:
            r = frac * dx
            nx = self._lat_offsets(self.X, b[:, 0], r)
            ny = self._lat_offsets(self.Y, b[:, 1], r)
            sx = abs(self.X[1] - self.X[0]); sy = abs(self.Y[1] - self.Y[0])
            for di in range(-nx, nx + 1):
                for dj in range(-ny, ny + 1):
                    if di == 0 and dj == 0:
                        continue
                    ii = np.clip(i + di, 0, len(self.X) - 1)
                    jj = np.clip(j + dj, 0, len(self.Y) - 1)
                    reach = ((np.abs(self.X[ii] - b[:, 0]) <= r + 0.5 * sx) &
                             (np.abs(self.Y[jj] - b[:, 1]) <= r + 0.5 * sy))
                    sel = np.flatnonzero(reach)
                    if not len(sel):
                        continue
                    kk0 = k0[sel]; kk1 = k1[sel]
                    acc = self.Vd[kk0, jj[sel], ii[sel]].copy()
                    for o in range(1, int((kk1 - kk0).max()) + 1 if len(sel) else 0):
                        s2 = np.flatnonzero((kk0 + o) <= kk1)
                        if not len(s2):
                            continue
                        acc[s2] = np.minimum(acc[s2], self.Vd[kk0[s2] + o, jj[sel][s2], ii[sel][s2]])
                    out[sel] = np.minimum(out[sel], acc)
        if self.M is None or self.source == "deck":
            return out
        deck_pool = out
        mi, mj = self._mij(b)
        mkb = self._mk(b[:, 2])
        m0 = np.minimum(self._mk(zhi), mkb)          # depth = -z: zhi -> shallowest
        m1 = np.maximum(self._mk(zlo), mkb)
        vm = self.Vm[m0, mj, mi].copy()
        for o in range(1, int((m1 - m0).max()) + 1 if len(m0) else 0):
            sel = np.flatnonzero((m0 + o) <= m1)
            if not len(sel):
                continue
            cand = self.Vm[m0[sel] + o, mj[sel], mi[sel]]
            vm[sel] = np.where(np.isfinite(cand), np.fmin(vm[sel], cand), vm[sel])
        if dx is not None:
            r = frac * dx
            slo = abs(self.mlon[1] - self.mlon[0]); sla = abs(self.mlat[1] - self.mlat[0])
            lo, la = self.fwd.transform(b[:, 0], b[:, 1])
            # metres -> degrees at this latitude (cheap, conservative)
            rlo = r / (111320.0 * np.cos(np.radians(la))); rla = r / 110540.0
            nlo = int(np.floor((rlo.max() + 0.5 * slo) / slo))
            nla = int(np.floor((rla.max() + 0.5 * sla) / sla))
            for di in range(-nlo, nlo + 1):
                for dj in range(-nla, nla + 1):
                    if di == 0 and dj == 0:
                        continue
                    ii = np.clip(mi + di, 0, len(self.mlon) - 1)
                    jj = np.clip(mj + dj, 0, len(self.mlat) - 1)
                    reach = ((np.abs(self.mlon[ii] - lo) <= rlo + 0.5 * slo) &
                             (np.abs(self.mlat[jj] - la) <= rla + 0.5 * sla))
                    sel = np.flatnonzero(reach)
                    if not len(sel):
                        continue
                    kk0 = m0[sel]; kk1 = m1[sel]
                    acc = self.Vm[kk0, jj[sel], ii[sel]].copy()
                    for o in range(1, int((kk1 - kk0).max()) + 1 if len(sel) else 0):
                        s2 = np.flatnonzero((kk0 + o) <= kk1)
                        if not len(s2):
                            continue
                        c2 = self.Vm[kk0[s2] + o, jj[sel][s2], ii[sel][s2]]
                        acc[s2] = np.where(np.isfinite(c2), np.fmin(acc[s2], c2), acc[s2])
                    vm[sel] = np.where(np.isfinite(acc), np.fmin(vm[sel], acc), vm[sel])
        if self.source == "muscal":
            return np.where(np.isfinite(vm), vm, deck_pool)
        return np.where(np.isfinite(vm), np.minimum(deck_pool, vm), deck_pool)
