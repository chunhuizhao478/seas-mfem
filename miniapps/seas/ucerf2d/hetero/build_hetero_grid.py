"""build_hetero_grid.py -- the 2D input file of the UCERF 2.5D heterogeneity runs (plan rev4, Step A).

    cd miniapps/seas
    ~/miniforge/envs/pythonenv/bin/python ucerf2d/hetero/build_hetero_grid.py --out <file.nc> [--dx 1000]

One NetCDF-4 file on a regular grid in the UCERF mesh frame (AEQD about 37.25 N, -119.25 E, R 6371 km),
covering the mesh box plus one cell:
  vp, vs, rho   CVM-S4.26.M01 (UCVM cvmsi) at z = 10 km; outside the model a linear taper over
                --taper_km to the deck constants (D3 = a)
  lambda, mu    rho (vp^2 - 2 vs^2), rho vs^2 [Pa], the moduli seas-mfem reads; exactly the deck's
                lambda = mu outside the taper band (so a run there is bit for bit the constant one)
  s_xx, s_yy, s_xy, s_zz   the so2 initial stress at 10 km (fault_tractions.py, plan eqs. 3.4-3.7):
                S_v' from the so2 deck, k = 1.45, CSM theta and R (UTM-11 triangulation, D10a),
                A_phi regime, theta from true north (D10b); effective stress, Pa, tension positive
Arrays are float64 (dims y, x), zlib-compressed.  The attributes record the frame, the models, the input
md5s, the recipe and the tool's commit.
"""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hetero_lib as H
import fault_tractions as FT

DECK = dict(lam=32038120320.0, mu=32038120320.0, rho=2670.0)       # [material_constant_fallback] of the UCERF decks


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def git_state():
    try:
        c = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=HERE).stdout.strip()
        d = subprocess.run(['git', 'status', '--porcelain', '--', HERE], capture_output=True, text=True, cwd=HERE).stdout.strip()
        return c, bool(d)
    except Exception:
        return 'unknown', True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--dx', type=float, default=1000.0)
    ap.add_argument('--taper_km', type=float, default=50.0)
    ap.add_argument('--k', type=float, default=FT.K_SO2)
    ap.add_argument('--chunk', type=int, default=200000)
    a = ap.parse_args()

    meta = json.load(open(FT.MSH.replace('.msh', '.json')))
    C = np.array(meta['tilt']['corners'])
    x = np.arange(np.floor(C[:, 0].min() / a.dx) - 1, np.ceil(C[:, 0].max() / a.dx) + 2) * a.dx
    y = np.arange(np.floor(C[:, 1].min() / a.dx) - 1, np.ceil(C[:, 1].max() / a.dx) + 2) * a.dx
    X, Y = np.meshgrid(x, y)
    lat, lon = H.inv(X.ravel(), Y.ravel())
    print('grid %d x %d = %d nodes' % (len(x), len(y), X.size), flush=True)

    # ---- velocity
    vp = np.empty(X.size); vs = np.empty(X.size); rho = np.empty(X.size)
    for i0 in range(0, X.size, a.chunk):
        s = slice(i0, i0 + a.chunk)
        vp[s], vs[s], rho[s] = H.ucvm_cvmsi(lon[s], lat[s], -FT.Z_S)
        print('  ucvm %d / %d' % (min(i0 + a.chunk, X.size), X.size), flush=True)
    inside = H.in_cvms(vp, vs).reshape(X.shape)
    from scipy.ndimage import distance_transform_edt
    dist = distance_transform_edt(~inside) * a.dx                         # m outside the model, 0 inside
    w = np.clip(1.0 - dist / (a.taper_km * 1e3), 0.0, 1.0).ravel()
    vs_d = np.sqrt(DECK['mu'] / DECK['rho']); vp_d = np.sqrt((DECK['lam'] + 2 * DECK['mu']) / DECK['rho'])
    vp = w * vp + (1 - w) * vp_d; vs = w * vs + (1 - w) * vs_d; rho = w * rho + (1 - w) * DECK['rho']
    # outside the taper band the values are exactly the deck's (w = 0 gives vp_d etc. up to one rounding):
    vp[w == 0] = vp_d; vs[w == 0] = vs_d; rho[w == 0] = DECK['rho']
    # the moduli seas-mfem reads (lambda, mu, rho): exactly the deck's outside the taper band
    mu = rho * vs ** 2; lam = rho * (vp ** 2 - 2 * vs ** 2)
    mu[w == 0] = DECK['mu']; lam[w == 0] = DECK['lam']

    # ---- stress (so2 method at 10 km)
    Sv = FT.sv_at(FT.Z_S)
    csm = H.CsmField(FT.CSM)
    sxx = np.empty(X.size); syy = np.empty(X.size); sxy = np.empty(X.size)
    gam = H.convergence_deg(lat, lon)
    for i0 in range(0, X.size, a.chunk):
        s = slice(i0, i0 + a.chunk)
        T, n, R, _, _ = csm.at(lat[s], lon[s])
        theta = H.shmax_deg(T)
        sH, sh = H.magnitudes_regime(Sv, a.k, R, n)
        S2 = H.horizontal_tensor(sH, sh, theta + gam[s]) * 1e6              # Pa, tension positive
        sxx[s], syy[s], sxy[s] = S2[:, 0, 0], S2[:, 1, 1], S2[:, 0, 1]
        print('  stress %d / %d' % (min(i0 + a.chunk, X.size), X.size), flush=True)
    szz = np.full(X.size, -Sv * 1e6)

    # ---- write
    import netCDF4
    commit, dirty = git_state()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    d = netCDF4.Dataset(a.out, 'w', format='NETCDF4')
    d.createDimension('y', len(y)); d.createDimension('x', len(x))
    vx = d.createVariable('x', 'f8', ('x',)); vx[:] = x; vx.units = 'm'
    vy = d.createVariable('y', 'f8', ('y',)); vy[:] = y; vy.units = 'm'
    fields = dict(vp=(vp, 'm/s'), vs=(vs, 'm/s'), rho=(rho, 'kg/m^3'), **{'lambda': (lam, 'Pa')}, mu=(mu, 'Pa'),
                  s_xx=(sxx, 'Pa'), s_yy=(syy, 'Pa'), s_xy=(sxy, 'Pa'), s_zz=(szz, 'Pa'))
    for k, (v, u) in fields.items():
        var = d.createVariable(k, 'f8', ('y', 'x'), zlib=True, complevel=4, chunksizes=(min(256, len(y)), min(256, len(x))))
        var[:] = v.reshape(X.shape); var.units = u
    att = dict(schema='seas_grid2d_v1', title='UCERF 2.5D heterogeneity input at z = 10 km (plan rev4)',
               frame='azimuthal equidistant, sphere R = 6371000 m, center lat 37.25 lon -119.25, x east, y north',
               z_m=FT.Z_S, layout='arrays (y, x); x, y ascending, uniform spacing',
               stress_convention='effective stress, Pa, tension positive (s_zz = -S_v\')',
               velocity_model='CVM-S4.26.M01 via UCVM cvmsi (ucvm_query -m cvmsi -c gd)',
               velocity_outside='linear taper over %g km to the deck constants (lambda = mu = %g Pa, rho = %g)' % (a.taper_km, DECK['mu'], DECK['rho']),
               stress_method='so2 METHOD.md Steps 1, 2, 9-11 at one depth; UTM-11 Delaunay; theta from true north; k = %g' % a.k,
               Sv_eff_MPa=Sv, csm=FT.CSM, csm_md5=md5(FT.CSM), sv_source=FT.SV_SOURCE,
               mesh=os.path.relpath(FT.MSH, FT.SEAS), mesh_md5=md5(FT.MSH), dx_m=a.dx,
               cvm_inside_share=float(inside.mean()), tool_commit=commit, tool_uncommitted=int(dirty),
               created_utc=datetime.datetime.now(datetime.UTC).strftime('%Y-%m-%dT%H:%M:%SZ'))
    for k, v in att.items():
        setattr(d, k, v)
    d.close()
    print('wrote', a.out, 'sha256', hashlib.sha256(open(a.out, 'rb').read()).hexdigest(), flush=True)


if __name__ == '__main__':
    main()
