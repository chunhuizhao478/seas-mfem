#!/usr/bin/env python3
"""Damage geometry per ParaView frame of a pff run, same metrics as the MOOSE Vista
analysis (analyze_damage.py on Andy's Mac): A(d>=.5), A(d>=.9) with element d = mean of
its 3 nodes; wall fraction/max/mean of d on nodes with |r - rb| < 1e-5; annulus
4.15-6 mm area with d>=.5; crack count on rings r = 8, 15 mm (nodes within +-0.3 mm,
d >= 0.9, angular runs; merged if gap < 2l/R); Gamma_l = int d/(c0 l) + l|grad d|^2/c0."""
import base64, glob, math, os, re, sys
import numpy as np

RB, L, C0 = 4.15e-3, 0.6e-3, 8.0 / 3.0
DT = {'Float64': np.float64, 'Float32': np.float32, 'Int32': np.int32, 'Int64': np.int64,
      'UInt8': np.uint8, 'UInt32': np.uint32, 'UInt64': np.uint64}


def arrays(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(r'<DataArray([^>]*)>\s*(\S*)\s*</DataArray>', txt, re.S):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        assert attrs.get('format') == 'binary', attrs
        b = m.group(2)
        hdr = base64.b64decode(b[:8])
        n = int(np.frombuffer(hdr, np.uint32)[0])
        data = base64.b64decode(b[8:])[:n]
        a = np.frombuffer(data, DT[attrs['type']])
        nc = int(attrs.get('NumberOfComponents', 1))
        out[attrs.get('Name', 'Points' if nc == 3 and 'Name' not in attrs else '?')] = a.reshape(-1, nc) if nc > 1 else a
    return out


def frame(cycdir):
    P, D, T = [], [], []
    off = 0
    for f in sorted(glob.glob(os.path.join(cycdir, 'proc*.vtu'))):
        a = arrays(f)
        pts = a['Points'] if 'Points' in a else [v for v in a.values() if v.ndim == 2 and v.shape[1] == 3][0]
        conn = a['connectivity'].astype(np.int64)
        offs = a['offsets'].astype(np.int64)
        assert np.all(np.diff(np.concatenate([[0], offs])) == 3), 'triangles only'
        P.append(pts[:, :2]); D.append(a['d']); T.append(conn.reshape(-1, 3) + off)
        off += len(pts)
    return np.vstack(P), np.concatenate(D), np.vstack(T)


def rings(xy, d, R, l=L):
    r = np.hypot(xy[:, 0], xy[:, 1])
    sel = np.abs(r - R) <= 0.3e-3
    th = np.arctan2(xy[sel, 1], xy[sel, 0]); on = d[sel] >= 0.9
    o = np.argsort(th); th, on = th[o], on[o]
    # unique angles (duplicated discontinuous points)
    th, idx = np.unique(np.round(th, 9), return_index=True); on = on[idx]
    if not on.any():
        return 0, 0
    runs = []
    i, n = 0, len(on)
    while i < n:
        if on[i]:
            j = i
            while j + 1 < n and on[j + 1]: j += 1
            runs.append((th[i], th[j])); i = j + 1
        else:
            i += 1
    if len(runs) > 1 and on[0] and on[-1]:
        runs[0] = (runs[-1][0], runs[0][1] + 2 * math.pi); runs.pop()
    raw = len(runs)
    gap = 2 * l / R
    merged = 1
    for k in range(1, len(runs)):
        if runs[k][0] - runs[k - 1][1] >= gap: merged += 1
    if len(runs) > 1 and (runs[0][0] + 2 * math.pi - runs[-1][1]) < gap: merged -= 1
    return raw, max(merged, 1)


def table(xy, d, tri):
    p = xy[tri]  # (ne,3,2)
    e1, e2 = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
    det = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
    area = 0.5 * np.abs(det)
    de = d[tri]
    dm = de.mean(1)
    # P1 gradient
    gx = ((de[:, 1] - de[:, 0]) * e2[:, 1] - (de[:, 2] - de[:, 0]) * e1[:, 1]) / det
    gy = (-(de[:, 1] - de[:, 0]) * e2[:, 0] + (de[:, 2] - de[:, 0]) * e1[:, 0]) / det
    gam = np.sum((dm / (C0 * L) + L * (gx ** 2 + gy ** 2) / C0) * area)
    cen = p.mean(1); rc = np.hypot(cen[:, 0], cen[:, 1])
    r = np.hypot(xy[:, 0], xy[:, 1])
    w = np.abs(r - RB) < 1e-5
    wk = np.unique(np.round(xy[w], 9), axis=0, return_index=True)[1]
    dw = d[w][wk]
    return dict(A05=area[dm >= 0.5].sum(), A09=area[dm >= 0.9].sum(),
                wall_frac=(dw >= 0.5).mean(), wall_dmax=dw.max(), wall_dmean=dw.mean(),
                A_ann=area[(dm >= 0.5) & (rc >= RB) & (rc <= 6e-3)].sum(),
                c8=rings(xy, d, 8e-3), c15=rings(xy, d, 15e-3), gamma=gam, nwall=len(dw),
                area=area.sum())


if __name__ == '__main__':
    root = sys.argv[1]
    pvd = glob.glob(os.path.join(root, '*.pvd'))[0]
    times = re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', open(pvd).read())
    print('t_us,A05,A09,wall_frac,wall_dmax,wall_dmean,A_ann,cracks8_raw,cracks8_merged,cracks15_raw,cracks15_merged,Gamma_l')
    for t, f in times:
        xy, d, tri = frame(os.path.join(root, os.path.dirname(f)))
        s = table(xy, d, tri)
        print(f"{float(t)*1e6:.3f},{s['A05']:.6e},{s['A09']:.6e},{s['wall_frac']:.4f},{s['wall_dmax']:.4f},"
              f"{s['wall_dmean']:.4f},{s['A_ann']:.6e},{s['c8'][0]},{s['c8'][1]},{s['c15'][0]},{s['c15'][1]},{s['gamma']:.6e}",
              flush=True)
    print(f"# nwall={s['nwall']} area={s['area']:.6e}", file=sys.stderr)
