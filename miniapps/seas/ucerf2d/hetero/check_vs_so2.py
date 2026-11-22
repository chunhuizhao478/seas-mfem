"""check_vs_so2.py -- gate S0 of plan rev2: our 10 km so2-method stress against the so2 deck file at z = -10 km.

    ~/miniforge/envs/pythonenv/bin/python check_vs_so2.py <fault_tractions.npz> <so2 stress_volume.nc> <out.json>

The so2 file is SeisSol's (tension-positive Pa, UTM-11 axes taken as east/north).  At every fault point inside
the so2 box: bilinear in x, y on the z = -10 km level, then sigma_H, sigma_h (compression +) and the SHmax
azimuth, against ours.  Differences come from (i) the Delaunay built in UTM-11 (so2) vs the mesh frame (here),
(ii) bilinear on the 1 km lattice vs exact at the point, (iii) so2's T_d0 edit near the NW fault.
"""
import json, sys
import numpy as np, netCDF4
from pyproj import Transformer

f = np.load(sys.argv[1]); d = netCDF4.Dataset(sys.argv[2])
x, y, z = d['x'][:], d['y'][:], d['z'][:]
kz = int(np.flatnonzero(np.abs(z + 10000.0) < 1e-6)[0])
lev = d['data'][kz]
X, Y = Transformer.from_crs('EPSG:4326', 'EPSG:32611', always_xy=True).transform(f['lon'], f['lat'])
inside = (X > x[0]) & (X < x[-1]) & (Y > y[0]) & (Y < y[-1])
i = np.clip(np.searchsorted(x, X) - 1, 0, len(x) - 2); j = np.clip(np.searchsorted(y, Y) - 1, 0, len(y) - 2)
tx = (X - x[i]) / (x[i + 1] - x[i]); ty = (Y - y[j]) / (y[j + 1] - y[j])
def bil(name):
    a = lev[name].astype(np.float64)
    return ((1 - tx) * (1 - ty) * a[j, i] + tx * (1 - ty) * a[j, i + 1] + (1 - tx) * ty * a[j + 1, i] + tx * ty * a[j + 1, i + 1])
sxx, syy, sxy = bil('s_xx') / 1e6, bil('s_yy') / 1e6, bil('s_xy') / 1e6
H = -np.stack([np.stack([sxx, sxy], -1), np.stack([sxy, syy], -1)], -2)     # compression +
w, U = np.linalg.eigh(H)
th_so2 = np.mod(np.degrees(np.arctan2(U[..., 0, 1], U[..., 1, 1])), 180.0)
sH_so2, sh_so2 = w[..., 1], w[..., 0]
# ours: sigma_H, sigma_h from the same eigen-split of our tensor rebuilt from the npz
from hetero_lib import magnitudes_regime
Sv = -float(np.mean(lev['s_zz'])) / 1e6
sH, sh = magnitudes_regime(Sv, 1.45, f['R'], f['regime'])
dth = np.abs(np.mod(f['shmax_deg'] - th_so2 + 90.0, 180.0) - 90.0)
m = inside
def st(v): return [float(np.percentile(v, p)) for p in (50, 95, 100)]
out = dict(points_inside_so2_box=int(m.sum()), of=int(len(m)),
           dtheta_deg_p50_p95_max=st(dth[m]),
           dsigmaH_MPa_p50_p95_max=st(np.abs(sH - sH_so2)[m]), dsigmah_MPa_p50_p95_max=st(np.abs(sh - sh_so2)[m]),
           rel_dsigmaH_p50_p95_max=st((np.abs(sH - sH_so2) / sH_so2)[m]))

# ---- tractions on our fault normals (gate S0b) ----
# so2's tensor is in UTM-11 axes.  c_u = angle of true north clockwise from the UTM y axis; gamma = the
# same from the mesh y axis.  A vector at angle a clockwise from UTM y is at a + (gamma - c_u) from mesh y,
# so the so2 tensor is rotated by beta = gamma - c_u into mesh axes.  so2 placed SHmax at theta from UTM y,
# i.e. at theta - c_u from true north; we place it at theta from true north.
import hetero_lib as HL
tr = Transformer.from_crs('EPSG:4326', 'EPSG:32611', always_xy=True)
Xn, Yn = tr.transform(f['lon'], f['lat'] + 0.01)
c_u = np.degrees(np.arctan2(Xn - X, Yn - Y))
gam = f['conv_deg']
b = np.radians(gam - c_u)
Rm = np.stack([np.stack([np.cos(b), np.sin(b)], -1), np.stack([-np.sin(b), np.cos(b)], -1)], -2)  # columns: UTM axes in mesh frame
S_utm = np.stack([np.stack([sxx, sxy], -1), np.stack([sxy, syy], -1)], -2)
S_mesh = np.einsum('...ia,...ab,...jb->...ij', Rm, S_utm, Rm)
# unit tangent per chain (no gradient across chains) and n = rot90(t)
tx_ = np.zeros_like(f['x']); ty_ = np.zeros_like(f['y'])
for c in np.unique(f['chain']):
    mc = f['chain'] == c
    gx_, gy_ = np.gradient(f['x'][mc]), np.gradient(f['y'][mc]); nn = np.hypot(gx_, gy_)
    tx_[mc], ty_[mc] = gx_ / nn, gy_ / nn
sn_so2, tau_so2 = HL.fault_traction(S_mesh, -ty_, tx_, tx_, ty_)
# our method with so2's convention (SHmax at theta - c_u from north)
sH2, sh2 = magnitudes_regime(Sv, 1.45, f['R'], f['regime'])
S_rep = HL.horizontal_tensor(sH2, sh2, f['shmax_deg'] - c_u + gam)
sn_rep, tau_rep = HL.fault_traction(S_rep, -ty_, tx_, tx_, ty_)
sn_our, tau_our = f['sigma_n_MPa'], f['tau_MPa']
out2 = dict(utm_convergence_deg=[float(np.min(c_u)), float(np.median(c_u)), float(np.max(c_u))],
            replica_vs_so2_dsn_MPa=st(np.abs(sn_rep - sn_so2)[inside]), replica_vs_so2_dtau_MPa=st(np.abs(tau_rep - tau_so2)[inside]),
            ours_vs_so2_dsn_MPa=st(np.abs(sn_our - sn_so2)[inside]), ours_vs_so2_dtau_MPa=st(np.abs(tau_our - tau_so2)[inside]),
            ours_vs_so2_dratio=st(np.abs(np.abs(tau_our) / sn_our - np.abs(tau_so2) / sn_so2)[inside]))
out.update(out2)
json.dump(out, open(sys.argv[3], 'w'), indent=1); print(json.dumps(out, indent=1))
np.savez(sys.argv[3].replace('.json', '.npz'), c_u=c_u, sn_so2=sn_so2, tau_so2=tau_so2, sn_rep=sn_rep, tau_rep=tau_rep, inside=inside)
