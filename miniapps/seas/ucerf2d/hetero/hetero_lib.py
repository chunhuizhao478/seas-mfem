"""hetero_lib.py -- shared pieces of the UCERF 2.5D heterogeneity tool (plan rev2).

Frame, fault points, CVM queries and the so2 (ShakeOut-2) stress method at one depth.
Python: ~/miniforge/envs/pythonenv (numpy, scipy).  Plain-text math only.

Mesh frame: azimuthal equidistant (spherical, R = 6371 km) about (37.25 N, -119.25 E),
x east, y north, metres (ucerf-2d/analyze.py proj; build_mesh_smooth.py LAT0, LON0).
Stress: tension positive in tensors; sigma_n' > 0 is compression; tau > 0 is right-lateral
for the normal n = rot90(t) (left of the chain direction t).
"""
import math
import os
import subprocess

import numpy as np

R_EARTH = 6371.0e3
LAT0, LON0 = 37.25, -119.25
UCVM = os.path.expanduser('~/Downloads/ucvm/ucvm_install')
CVMS_BACKGROUND = (6300.0, 3637.307, 2859.770)     # cvmsi's value outside its box (queried 2026-11-20)


# ------------------------------------------------------------------ frame
def fwd(lat, lon):
    """(lat, lon) deg -> mesh-frame (x, y) m."""
    p, l = np.radians(lat), np.radians(lon)
    p0, l0 = math.radians(LAT0), math.radians(LON0)
    c = np.arccos(np.clip(np.sin(p0) * np.sin(p) + np.cos(p0) * np.cos(p) * np.cos(l - l0), -1, 1))
    k = np.where(c > 1e-12, c / np.sin(np.maximum(c, 1e-300)), 1.0)
    x = R_EARTH * k * np.cos(p) * np.sin(l - l0)
    y = R_EARTH * k * (np.cos(p0) * np.sin(p) - np.sin(p0) * np.cos(p) * np.cos(l - l0))
    return x, y


def inv(x, y):
    """mesh-frame (x, y) m -> (lat, lon) deg."""
    p0, l0 = math.radians(LAT0), math.radians(LON0)
    rho = np.hypot(x, y)
    c = rho / R_EARTH
    with np.errstate(invalid='ignore', divide='ignore'):
        lat = np.arcsin(np.cos(c) * math.sin(p0) + np.where(rho > 0, y * np.sin(c) * math.cos(p0) / rho, 0))
        lon = l0 + np.arctan2(x * np.sin(c), rho * math.cos(p0) * np.cos(c) - y * math.sin(p0) * np.sin(c))
    return np.degrees(lat), np.degrees(lon)


def convergence_deg(lat, lon):
    """Angle (deg, clockwise) from the mesh y axis to geographic north at (lat, lon).
    A direction at azimuth theta (clockwise from north) is at theta + convergence from the mesh y axis."""
    x0, y0 = fwd(lat, lon)
    xn, yn = fwd(lat + 0.01, lon)
    return np.degrees(np.arctan2(xn - x0, yn - y0))


# ------------------------------------------------------------------ fault points
def read_fault_edges(path, tag=3):
    nodes, edges, sec = {}, [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith('$'):
                sec = s
                continue
            t = s.split()
            if sec == '$Nodes' and len(t) == 4:
                nodes[int(t[0])] = (float(t[1]), float(t[2]))
            elif sec == '$Elements' and len(t) > 3 and t[1] == '1' and int(t[3]) == tag:
                nt = int(t[2])
                edges.append((int(t[3 + nt]), int(t[4 + nt])))
    return nodes, edges


def fault_chains(nodes, edges):
    """Join the fault edges into polylines, split where a node does not have two edges."""
    from collections import defaultdict
    adj = defaultdict(list)
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    used, out = set(), []
    for s in [n for n in adj if len(adj[n]) != 2]:
        for n in adj[s]:
            if frozenset((s, n)) in used:
                continue
            line = [s, n]
            used.add(frozenset((s, n)))
            while len(adj[n]) == 2:
                nx = [m for m in adj[n] if frozenset((n, m)) not in used]
                if not nx:
                    break
                used.add(frozenset((n, nx[0])))
                line.append(nx[0])
                n = nx[0]
            out.append(np.array([nodes[i] for i in line]))
    return out


def resample_chains(chains, ds=1000.0):
    """Points every ds along each chain; returns dict of arrays (x, y, tx, ty, nx, ny, chain, s)."""
    P, T, C, S = [], [], [], []
    for k, c in enumerate(chains):
        d = np.r_[0, np.cumsum(np.hypot(*np.diff(c, axis=0).T))]
        s = np.arange(0.0, d[-1], ds)
        x = np.interp(s, d, c[:, 0])
        y = np.interp(s, d, c[:, 1])
        tx, ty = np.gradient(x), np.gradient(y)
        n = np.hypot(tx, ty)
        P.append(np.c_[x, y]); T.append(np.c_[tx / n, ty / n]); C += [k] * len(s); S.append(s)
    P, T = np.vstack(P), np.vstack(T)
    return dict(x=P[:, 0], y=P[:, 1], tx=T[:, 0], ty=T[:, 1], nx=-T[:, 1], ny=T[:, 0],
                chain=np.array(C), s=np.concatenate(S))


# ------------------------------------------------------------------ CVM (UCVM cvmsi)
def ucvm_cvmsi(lon, lat, depth_m):
    """Vp, Vs, rho (float64 arrays) from UCVM cvmsi at geographic depth depth_m (scalar or array)."""
    dep = np.broadcast_to(np.asarray(depth_m, float), np.shape(lon))
    inp = ''.join('%.6f %.6f %.1f\n' % (a, b, c) for a, b, c in zip(lon, lat, dep))
    r = subprocess.run(['bash', '-c', 'source conf/ucvm_env.sh >/dev/null 2>&1; '
                        'bin/ucvm_query -f conf/ucvm.conf -m cvmsi -c gd'],
                       input=inp, capture_output=True, text=True, cwd=UCVM)
    rows = [l.split() for l in r.stdout.splitlines() if len(l.split()) >= 17 and l.split()[5] == 'cvmsi']
    if len(rows) != len(lon):
        raise RuntimeError('ucvm_query returned %d rows for %d points: %s' % (len(rows), len(lon), r.stderr[-400:]))
    a = np.array([[float(t[14]), float(t[15]), float(t[16])] for t in rows])
    return a[:, 0], a[:, 1], a[:, 2]


def in_cvms(vp, vs):
    return ~((np.abs(vp - CVMS_BACKGROUND[0]) < 1e-3) & (np.abs(vs - CVMS_BACKGROUND[1]) < 1e-2))


# ------------------------------------------------------------------ so2 stress method (one level)
def read_csm(path):
    """YHSM-2013 csv -> lon, lat, (N,3,3) tension-positive ENU tensors [MPa], A_phi (col 14)."""
    rows = [l.split(',') for l in open(path) if not l.startswith('#')][1:]
    lon = np.array([float(r[0].strip('"')) for r in rows])
    lat = np.array([float(r[1]) for r in rows])
    S = np.array([[float(r[k]) for k in range(3, 9)] for r in rows])   # See Sen Seu Snn Snu Suu
    T = np.empty((len(S), 3, 3))
    T[:, 0, 0] = S[:, 0]; T[:, 0, 1] = T[:, 1, 0] = S[:, 1]; T[:, 0, 2] = T[:, 2, 0] = S[:, 2]
    T[:, 1, 1] = S[:, 3]; T[:, 1, 2] = T[:, 2, 1] = S[:, 4]; T[:, 2, 2] = S[:, 5]
    aphi = np.array([float(r[13]) for r in rows])
    return lon, lat, T, aphi


def utm11(lat, lon):
    """(lat, lon) deg -> UTM zone 11N (EPSG:32611) x, y m, always_xy, as so2."""
    from pyproj import Transformer
    tr = Transformer.from_crs('EPSG:4326', 'EPSG:32611', always_xy=True)
    return tr.transform(np.asarray(lon), np.asarray(lat))


class CsmField:
    """so2 Steps 2 and 10 at query points (METHOD.md eqs. shape, interp).

    The Delaunay triangulation of the csv nodes is built in UTM-11, as so2 does (plan rev4, D10a); the
    query points are given as (lat, lon).  Vectorized over the query points."""

    def __init__(self, csv_path):
        from scipy.spatial import Delaunay, cKDTree
        lon, lat, T, aphi = read_csm(csv_path)
        self.x, self.y = utm11(lat, lon)
        self.T, self.aphi = T, aphi
        self.n_node = np.minimum(np.floor(aphi), 2).astype(int)          # so2 eq. (aphi)
        self.tri = Delaunay(np.c_[self.x, self.y])
        self.kd = cKDTree(np.c_[self.x, self.y])

    @staticmethod
    def _shape_R(Tq):
        w = np.linalg.eigvalsh(Tq)                                         # ascending: w0 = -sigma1
        return np.clip((w[..., 1] - w[..., 0]) / (w[..., 2] - w[..., 0]), 0.0, 1.0)

    def at(self, lat, lon):
        """Interpolated tension-positive ENU tensor, regime n, shape R, mixed and in-hull flags."""
        qx, qy = utm11(lat, lon)
        q = np.c_[qx, qy]
        nq = len(q)
        simp = self.tri.find_simplex(q)
        hull = simp >= 0
        Tq = np.empty((nq, 3, 3)); n = np.empty(nq, int); R = np.empty(nq); mixed = np.zeros(nq, bool)
        if hull.any():
            s = simp[hull]
            V = self.tri.simplices[s]                                      # (m, 3)
            Tr = self.tri.transform[s]                                     # (m, 3, 2)
            b = np.einsum('mij,mj->mi', Tr[:, :2], q[hull] - Tr[:, 2])
            lam = np.c_[b, 1 - b.sum(1)]                                   # (m, 3)
            Ti = np.einsum('mk,mkij->mij', lam, self.T[V])
            Rt = self._shape_R(Ti)
            nv = self.n_node[V]
            same = (nv[:, 0] == nv[:, 1]) & (nv[:, 0] == nv[:, 2])
            ab = np.einsum('mk,mk->m', lam, self.aphi[V])
            nm = np.minimum(np.floor(ab), 2).astype(int)
            Rm = np.where(nm == 0, 1 - ab, np.where(nm == 1, ab - 1, 3 - ab))
            Tq[hull] = Ti
            n[hull] = np.where(same, nv[:, 0], nm)
            R[hull] = np.where(same, Rt, Rm)
            mixed[hull] = ~same
        out = ~hull
        if out.any():                                                      # nearest node outside the hull
            j = self.kd.query(q[out])[1]
            Tq[out] = self.T[j]; n[out] = self.n_node[j]; R[out] = self._shape_R(self.T[j])
        return Tq, n, R, mixed, hull


def shmax_deg(T_tension):
    """Azimuth (deg E of N, [0,180)) of the most compressive eigenvector of the horizontal block
    of the compression-positive tensor (so2 toolbox shmax_from_sigma)."""
    h = -T_tension[..., :2, :2]
    w, U = np.linalg.eigh(h)
    v = U[..., :, 1]
    return np.mod(np.degrees(np.arctan2(v[..., 0], v[..., 1])), 180.0)


def magnitudes_regime(Sv, k, R, n):
    """so2 eq. (regime): (sigma_H, sigma_h) compression-positive, MPa."""
    P = (1.0 - R) * k + R
    sH = np.where(n == 0, P * Sv / k, np.where(n == 1, k * Sv / P, k * Sv))
    sh = np.where(n == 0, Sv / k, np.where(n == 1, Sv / P, P * Sv))
    return sH, sh


def horizontal_tensor(sH, sh, alpha_deg):
    """Tension-positive 2x2 tensor(s) in the mesh frame; alpha = SHmax azimuth from the mesh y axis."""
    a = np.radians(alpha_deg)
    eH = np.stack([np.sin(a), np.cos(a)], -1)
    eh = np.stack([np.cos(a), -np.sin(a)], -1)
    return -(sH[..., None, None] * eH[..., :, None] * eH[..., None, :]
             + sh[..., None, None] * eh[..., :, None] * eh[..., None, :])


def fault_traction(S2, nx, ny, tx, ty):
    """sigma_n' (compression +) and tau (right-lateral +) on the fault, same units as S2."""
    n = np.stack([nx, ny], -1); t = np.stack([tx, ty], -1)
    tr = np.einsum('...ij,...j->...i', S2, n)
    return -np.einsum('...i,...i->...', tr, n), np.einsum('...i,...i->...', tr, t)
