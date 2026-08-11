#!/usr/bin/env python3
"""muscal_vs.py -- Vs sampled from MUSCAL, the SOURCE velocity model.

WHY THIS EXISTS

The gate has been judged against the deck's `safs_material_cvm.nc`, which is a
RESAMPLE of MUSCAL onto a 1500 m lateral / **250 m uniform vertical** lattice.
MUSCAL itself is 0.01 deg laterally with a **50 m** depth step through the top
500 m.  Nearest-grid on the deck file therefore hands every barycentre in the
top 125 m the z = 0 value -- Vs p50 526 m/s -- while MUSCAL says the rock at
100 m depth is p50 1,429 m/s and at 200 m is 2,091 m/s.  Sizing a mesh against
that step means paying for a VERTICAL BINNING ARTEFACT, and it is the single
thing that sets the collar's price.

Interface mirrors collar_lib.VsGrid / leb_collar.PoolVs so it can be dropped
into the census and the refiner:
    at(P)                nearest-grid Vs at points P (UTM 11N metres)
    bbox_min(lo, hi)     min Vs over each axis-aligned box  (refinement target)

Conventions:
  * mesh z is NEGATIVE down; MUSCAL `depth` is POSITIVE down.  depth = -z.
  * MUSCAL is NaN outside its valid box and over water (37 % of the surface
    level).  Those samples fall back to a supplied backup grid -- the deck nc --
    because that is exactly what ASAGI serves there at runtime.  The fallback
    count is reported; it is never silently zero-filled.
"""
import numpy as np

MUSCAL = "/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc"
UTM11N = "EPSG:32611"


# The SAFS ALT domain box in UTM 11N, plus a degree of margin.  Cropping to it
# takes MUSCAL from 210x1251x1301 (1.37 GB) to roughly 130x440x810 (185 MB),
# which matters because bbox_min() caches scipy minimum_filter windows and each
# one is a full copy of the cube.
DOMAIN_LONLAT = (-122.6, -113.4, 31.4, 36.6)
DOMAIN_DEPTH = (0.0, 45000.0)


class MuscalVs:
    def __init__(self, path=MUSCAL, backup=None, verbose=True,
                 lonlat=DOMAIN_LONLAT, depth=DOMAIN_DEPTH):
        import h5py
        from pyproj import Transformer
        with h5py.File(path, "r") as f:
            lon = f["longitude"][:].astype(np.float64)
            lat = f["latitude"][:].astype(np.float64)
            dep = f["depth"][:].astype(np.float64)
            klo = np.nonzero((lon >= lonlat[0]) & (lon <= lonlat[1]))[0]
            kla = np.nonzero((lat >= lonlat[2]) & (lat <= lonlat[3]))[0]
            kd = np.nonzero((dep >= depth[0]) & (dep <= depth[1]))[0]
            self.lon, self.lat, self.dep = lon[klo], lat[kla], dep[kd]
            self.vs = np.asarray(
                f["vs"][kd[0]:kd[-1] + 1, kla[0]:kla[-1] + 1,
                        klo[0]:klo[-1] + 1]).astype(np.float32)
        o = np.argsort(self.dep)
        self.dep, self.vs = self.dep[o], self.vs[o]
        if verbose:
            print(f"[muscal] cropped to lon {self.lon[0]:.2f}..{self.lon[-1]:.2f}, "
                  f"lat {self.lat[0]:.2f}..{self.lat[-1]:.2f}, "
                  f"depth {self.dep[0]:,.0f}..{self.dep[-1]:,.0f} m "
                  f"({self.vs.nbytes/2**20:,.0f} MB)")
        self.n_nan = int((~np.isfinite(self.vs)).sum())
        self._filled = None            # built lazily: it is another 1.4 GB and
        # only bbox_min() needs it, while the census only calls at()
        self.to_ll = Transformer.from_crs(UTM11N, "EPSG:4326", always_xy=True)
        self.backup = backup
        self.n_fallback = 0
        self.n_query = 0
        if verbose:
            print(f"[muscal] {self.vs.shape} (depth,lat,lon)  "
                  f"depth {self.dep.min():,.0f}..{self.dep.max():,.0f} m, "
                  f"top step {np.diff(self.dep)[0]:,.0f} m")
            print(f"[muscal] NaN {100*self.n_nan/self.vs.size:.1f} % of the cube "
                  f"(outside the valid box / offshore)")
        self._cache = {}
        self._max = 6

    # ---- index helpers ----------------------------------------------------
    @staticmethod
    def _near(axis, q):
        i = np.clip(np.searchsorted(axis, q), 0, len(axis) - 1)
        j = np.clip(i - 1, 0, len(axis) - 1)
        return np.where(np.abs(axis[j] - q) <= np.abs(axis[i] - q), j, i)

    def _ijk(self, P):
        lo, la = self.to_ll.transform(P[:, 0], P[:, 1])
        return (self._near(self.lon, lo).astype(np.int32),
                self._near(self.lat, la).astype(np.int32),
                self._near(self.dep, -P[:, 2]).astype(np.int32))

    # ---- the gate's sampler ----------------------------------------------
    def at(self, P):
        ilon, ilat, idep = self._ijk(P)
        v = self.vs[idep, ilat, ilon].astype(np.float64)
        bad = ~np.isfinite(v)
        self.n_query += len(v)
        if bad.any():
            self.n_fallback += int(bad.sum())
            if self.backup is not None:
                v[bad] = self.backup.at(P[bad])
            else:
                v[bad] = np.nan
        return v

    # ---- the refinement target -------------------------------------------
    def _filt(self, kd, kla, klo):
        if self._filled is None:
            # NaN -> +inf so a minimum_filter ignores it instead of poisoning
            # every window that touches the offshore mask
            self._filled = np.where(np.isfinite(self.vs), self.vs,
                                    np.float32(np.inf))
        key = (kd, kla, klo)
        if key not in self._cache:
            from scipy.ndimage import minimum_filter
            if len(self._cache) >= self._max:
                self._cache.pop(next(iter(self._cache)))
            self._cache[key] = (self._filled if (kd == 0 and kla == 0 and klo == 0)
                                else minimum_filter(self._filled,
                                                    size=(2 * kd + 1, 2 * kla + 1,
                                                          2 * klo + 1),
                                                    mode="nearest"))
        return self._cache[key]

    @staticmethod
    def _dyadic(k):
        k = np.maximum(k, 0)
        out = np.zeros_like(k)
        nz = k > 0
        out[nz] = 1 << np.ceil(np.log2(k[nz])).astype(np.int32)
        return out

    def bbox_min(self, lo, hi):
        """min Vs over each box -- a lower bound on every descendant cell."""
        i0lo, i0la, i0d = self._ijk(lo)
        i1lo, i1la, i1d = self._ijk(hi)
        # depth index runs OPPOSITE to z, so order the pairs
        d0, d1 = np.minimum(i0d, i1d), np.maximum(i0d, i1d)
        a0, a1 = np.minimum(i0la, i1la), np.maximum(i0la, i1la)
        o0, o1 = np.minimum(i0lo, i1lo), np.maximum(i0lo, i1lo)
        cd, cla, clo = (d0 + d1) // 2, (a0 + a1) // 2, (o0 + o1) // 2
        kd = self._dyadic(np.maximum(cd - d0, d1 - cd))
        kla = self._dyadic(np.maximum(cla - a0, a1 - cla))
        klo = self._dyadic(np.maximum(clo - o0, o1 - clo))
        out = np.empty(len(lo))
        key = ((kd.astype(np.int64) << np.int64(40))
               + (kla.astype(np.int64) << np.int64(20)) + klo)
        for u in np.unique(key):
            m = key == u
            out[m] = self._filt(int(kd[m][0]), int(kla[m][0]),
                                int(klo[m][0]))[cd[m], cla[m], clo[m]]
        bad = ~np.isfinite(out)
        if bad.any():
            self.n_fallback += int(bad.sum())
            if self.backup is not None:
                mid = 0.5 * (np.asarray(lo) + np.asarray(hi))
                out[bad] = self.backup.bbox_min(np.asarray(lo)[bad],
                                                np.asarray(hi)[bad]) \
                    if hasattr(self.backup, "bbox_min") else self.backup.at(mid[bad])
        return out

    def report(self):
        if self.n_query:
            print(f"[muscal] {self.n_fallback:,} of {self.n_query:,} samples fell "
                  f"back to the deck nc ({100*self.n_fallback/self.n_query:.2f} %)")
