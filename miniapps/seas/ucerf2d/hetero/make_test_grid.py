"""make_test_grid.py -- seas_grid2d_v1 files for the H2/H3 gates (plan rev4).

    ~/miniforge/envs/pythonenv/bin/python make_test_grid.py --deck <deck.toml> --out <file.nc>
           [--dx 2000] [--perturb 0.0] [--wavelength 20e3]

Covers the deck's mesh (bounding box + 2 cells) with the deck's own constants:
[material_constant_fallback] lambda, mu, rho and the [stress] constant_tensor (s_xx, s_yy, s_xy, s_zz).
--perturb a > 0 multiplies mu and lambda by 1 + a sin(2 pi x / L) cos(2 pi y / L) (rho kept), for a
heterogeneous smoke test; --perturb 0 gives the uniform file of gate U4 (a run with it must be bit for
bit the constant-material run).
"""
import argparse, os
import numpy as np
import netCDF4
try:
    import tomllib
except ImportError:
    import tomli as tomllib


def mesh_bbox(path):
    xs, ys, sec = [], [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith('$'):
                sec = s; continue
            t = s.split()
            if sec == '$Nodes' and len(t) == 4:
                xs.append(float(t[1])); ys.append(float(t[2]))
    return min(xs), max(xs), min(ys), max(ys)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--deck', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--dx', type=float, default=2000.0)
    ap.add_argument('--perturb', type=float, default=0.0)
    ap.add_argument('--wavelength', type=float, default=20e3)
    a = ap.parse_args()
    deck = tomllib.load(open(a.deck, 'rb'))
    seas = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    x0, x1, y0, y1 = mesh_bbox(os.path.join(seas, deck['mesh']['path']))
    x = np.arange(np.floor(x0 / a.dx) - 2, np.ceil(x1 / a.dx) + 3) * a.dx
    y = np.arange(np.floor(y0 / a.dx) - 2, np.ceil(y1 / a.dx) + 3) * a.dx
    X, Y = np.meshgrid(x, y)
    m = deck['material_constant_fallback']; s = deck['stress']
    lam = np.full(X.shape, float(m['lambda'])); mu = np.full(X.shape, float(m['mu'])); rho = np.full(X.shape, float(m['rho']))
    if a.perturb > 0:
        g = 1.0 + a.perturb * np.sin(2 * np.pi * X / a.wavelength) * np.cos(2 * np.pi * Y / a.wavelength)
        lam = lam * g; mu = mu * g
    F = dict(rho=rho, mu=mu, vp=np.sqrt((lam + 2 * mu) / rho), vs=np.sqrt(mu / rho), s_xx=np.full(X.shape, float(s['sigma_xx_pa'])),
             s_yy=np.full(X.shape, float(s['sigma_yy_pa'])), s_xy=np.full(X.shape, float(s['sigma_xy_pa'])),
             s_zz=np.full(X.shape, float(s['sigma_zz_pa'])))
    F['lambda'] = lam
    d = netCDF4.Dataset(a.out, 'w', format='NETCDF4')
    d.createDimension('y', len(y)); d.createDimension('x', len(x))
    d.createVariable('x', 'f8', ('x',))[:] = x; d.createVariable('y', 'f8', ('y',))[:] = y
    for k, v in F.items():
        d.createVariable(k, 'f8', ('y', 'x'), zlib=True)[:] = v
    d.schema = 'seas_grid2d_v1'; d.title = 'test grid from %s, perturb %g' % (a.deck, a.perturb)
    d.close(); print('wrote', a.out, len(x), 'x', len(y))


if __name__ == '__main__':
    main()
