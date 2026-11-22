"""qc_hetero_grid.py -- QC of the build_hetero_grid.py file (plan rev4, phase H0).

    ~/miniforge/envs/pythonenv/bin/python qc_hetero_grid.py <file.nc> <fault_tractions.npz> <out_dir>

Q1  the file, read bilinearly at the 825 fault points, against fault_tractions.py's direct values
    (sigma_n', tau and Vs): differences are the 1 km bilinear error only.
Q2  outside the taper band every node holds the deck constants exactly, and s_zz = -S_v' everywhere.
Maps of Vs and of the SHmax azimuth / max shear stress over the box; a JSON summary.
"""
import json, os, sys
import numpy as np
import netCDF4
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hetero_lib as H
import build_hetero_grid as BG


def bilinear(x, y, F, qx, qy):
    i = np.clip(np.searchsorted(x, qx) - 1, 0, len(x) - 2); j = np.clip(np.searchsorted(y, qy) - 1, 0, len(y) - 2)
    tx = (qx - x[i]) / (x[i + 1] - x[i]); ty = (qy - y[j]) / (y[j + 1] - y[j])
    f0 = F[j, i] + tx * (F[j, i + 1] - F[j, i]); f1 = F[j + 1, i] + tx * (F[j + 1, i + 1] - F[j + 1, i])
    return f0 + ty * (f1 - f0)


def main():
    path, ftp, out = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    d = netCDF4.Dataset(path); x = d['x'][:].data; y = d['y'][:].data
    F = {k: d[k][:].data for k in ('vp', 'vs', 'rho', 'lambda', 'mu', 's_xx', 's_yy', 's_xy', 's_zz')}
    f = np.load(ftp)
    qx, qy = f['x'], f['y']
    sxx, syy, sxy = (bilinear(x, y, F[k], qx, qy) / 1e6 for k in ('s_xx', 's_yy', 's_xy'))
    S2 = np.stack([np.stack([sxx, sxy], -1), np.stack([sxy, syy], -1)], -2)
    tx = np.zeros_like(qx); ty = np.zeros_like(qy)
    for c in np.unique(f['chain']):
        m = f['chain'] == c; gx, gy = np.gradient(qx[m]), np.gradient(qy[m]); nn = np.hypot(gx, gy); tx[m], ty[m] = gx / nn, gy / nn
    sn, tau = H.fault_traction(S2, -ty, tx, tx, ty)
    vs = bilinear(x, y, F['vs'], qx, qy)
    def st(v): v = np.abs(v); return [float(np.median(v)), float(np.percentile(v, 95)), float(v.max())]
    vs_d = np.sqrt(BG.DECK['mu'] / BG.DECK['rho']); vp_d = np.sqrt((BG.DECK['lam'] + 2 * BG.DECK['mu']) / BG.DECK['rho'])
    far = (F['vs'] == vs_d) & (F['vp'] == vp_d) & (F['rho'] == BG.DECK['rho']) & (F['mu'] == BG.DECK['mu']) & (F['lambda'] == BG.DECK['lam'])
    q = dict(file=path, nx=len(x), ny=len(y), dx=float(x[1] - x[0]),
             Q1_dsigma_n_MPa=st(sn - f['sigma_n_MPa']), Q1_dtau_MPa=st(tau - f['tau_MPa']),
             Q1_dvs_rel=st(vs / f['vs'] - 1),
             Q2_nodes_exact_deck=int(far.sum()), Q2_share=float(far.mean()),
             Q2_szz_constant=bool(np.all(F['s_zz'] == F['s_zz'][0, 0])), szz_MPa=float(F['s_zz'][0, 0] / 1e6),
             vs_range=[float(F['vs'].min()), float(F['vs'].max())], rho_range=[float(F['rho'].min()), float(F['rho'].max())],
             attrs={k: str(getattr(d, k)) for k in d.ncattrs()})
    json.dump(q, open(os.path.join(out, 'qc.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in q.items() if k != 'attrs'}, indent=1))
    # maps
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    meta = json.load(open(BG.FT.MSH.replace('.msh', '.json'))); C = np.array(meta['tilt']['corners']) / 1e3
    tmax = np.sqrt(((F['s_xx'] - F['s_yy']) / 2) ** 2 + F['s_xy'] ** 2) / 1e6
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    for a, v, lab, cm in ((ax[0], F['vs'] / 1e3, 'Vs (km/s) at 10 km', 'viridis'), (ax[1], tmax, 'max horizontal shear (MPa)', 'magma')):
        s = (slice(None, None, 4), slice(None, None, 4))
        im = a.pcolormesh(x[s[1]] / 1e3, y[s[0]] / 1e3, v[s], shading='nearest', cmap=cm)
        fig.colorbar(im, ax=a, shrink=0.8, label=lab)
        a.add_patch(Polygon(C, closed=True, fill=False, ec='w', lw=1))
        for ch in range(5):
            m = f['chain'] == ch; a.plot(qx[m] / 1e3, qy[m] / 1e3, color='w', lw=1)
        a.set_aspect('equal'); a.set_xlabel('x (km)'); a.set_ylabel('y (km)')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'qc_maps.png'), dpi=150)


if __name__ == '__main__':
    main()
