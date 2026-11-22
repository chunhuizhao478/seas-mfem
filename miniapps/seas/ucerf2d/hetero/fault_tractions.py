"""fault_tractions.py -- the so2 (ShakeOut-2) prestress at 10 km on the UCERF 2.5D faults,
its on-fault tractions and the nucleation length L_b (plan rev2, D4-D6).

    cd miniapps/seas
    ~/miniforge/envs/pythonenv/bin/python ucerf2d/hetero/fault_tractions.py --out <dir>

Method (so2_setup/so2_stress/docs/METHOD.md, Steps 1, 2, 9-11; read 2026-11-21):
  S_v'  = -s_zz at z = -10 km of the current so2 deck's stress_volume.nc (171.116 MPa: lateral-mean cvmsi
          density load minus hydrostatic pore pressure, rebuilt 2026-11-06; METHOD.md's older MUSCAL-density
          profile gives 161.934 MPa);
  theta, R per point from the YHSM-2013 tensor, Delaunay-linear (Step 2), rotated into the mesh frame
          by the grid convergence;
  regime n from A_phi, A_phi-bar in mixed triangles (Steps 9-10);
  (sigma_H, sigma_h) = eq. (regime) with k = sigma_1 / sigma_3 = 1.45 (Step 11), then the horizontal tensor;
  sigma_n', tau on each fault point (n = rot90(t)); tau > 0 right-lateral.
Not applied (2.5D, vertical faults): the T_d0 (dip traction) edit of Steps 5-8 and 12.
L_b = mu D_c / (b sigma_n') (the deck's convention, mu not mu'), D_c = 0.26 m, b = 0.03, mu = rho Vs^2 from
CVM-S4.26.M01 at 10 km (UCVM cvmsi).
"""
import argparse
import json
import os
import subprocess

import numpy as np

import hetero_lib as H

HERE = os.path.dirname(os.path.abspath(__file__))
SEAS = os.path.dirname(os.path.dirname(HERE))
MSH = os.path.join(SEAS, 'ucerf2d/mesh/ucerf2d_saf_sjf_500m_tilted_D415km_R70.msh')
CSM = os.path.expanduser('~/Downloads/so2/so2_setup/so2_stress/inputs/CSM_orientation.csv')
SO2_DECK = os.path.expanduser('~/Downloads/so2/so2_inputs/'
                              'alt_c200_ctm100_visco_plast_aphi_cvmsi_phi3040_dc0p2_pph1lam045_fw035h2nw_h4c005_'
                              'h6k100_h7k100_vw360n_t200_1hz')
SV_SOURCE = os.path.join(SO2_DECK, 'stress_volume.nc')   # md5 24877d65..., S_v from the cvmsi rho (2026-11-06)
Z_S = -10000.0
K_SO2 = 1.45
DC, B = 0.26, 0.03
DECK_TENSOR = np.array([[-22723874.717053052, -7382134.744727842],
                        [-7382134.744727842, -81940125.28294693]]) / 1e6      # today's [stress], MPa
CHAIN_NAMES = ['SAF Parkfield-San Gorgonio', 'SAF Coachella', 'SAF Banning', 'SJF north', 'SJF south']


def sv_at(z):
    """S_v' [MPa] of the current so2 decks at level z: -s_zz of stress_volume.nc, which is laterally constant."""
    import netCDF4
    d = netCDF4.Dataset(SV_SOURCE)
    k = np.flatnonzero(np.abs(d['z'][:] - z) < 1e-6)
    assert len(k) == 1, 'level %g not in %s' % (z, SV_SOURCE)
    szz = np.asarray(d['data'][int(k[0])]['s_zz'], np.float64)
    assert szz.max() - szz.min() < 1e-3 * abs(szz.mean()), 'S_v not laterally constant'
    return float(-szz.mean() / 1e6)


def md5(p):
    import hashlib
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--k', type=float, default=K_SO2)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    F = H.resample_chains(H.fault_chains(*H.read_fault_edges(MSH)))
    lat, lon = H.inv(F['x'], F['y'])
    gam = H.convergence_deg(lat, lon)
    Sv = sv_at(Z_S)

    csm = H.CsmField(CSM)
    T, n, R, mixed, hull = csm.at(lat, lon)
    theta = H.shmax_deg(T)                       # deg E of N
    alpha = theta + gam                          # SHmax clockwise from the mesh y axis (north is at gam from it)

    def tractions(k, sv):
        sH, sh = H.magnitudes_regime(sv, k, R, n)
        S2 = H.horizontal_tensor(sH, sh, alpha)
        return H.fault_traction(S2, F['nx'], F['ny'], F['tx'], F['ty'])

    sn, tau = tractions(a.k, Sv)
    sn0, tau0 = H.fault_traction(np.broadcast_to(DECK_TENSOR, (len(lat), 2, 2)), F['nx'], F['ny'], F['tx'], F['ty'])

    vp, vs, rho = H.ucvm_cvmsi(lon, lat, -Z_S)
    mu = rho * vs ** 2
    Lb = mu * DC / (B * sn * 1e6)
    Lb_deck = 32.03812032e9 * DC / (B * sn0 * 1e6)
    Lb_cvm_deckstress = mu * DC / (B * sn0 * 1e6)

    # k sweep: the largest |tau|/sigma_n' any point reaches, and the share above f0 = 0.6
    sweep = []
    for k in (1.45, 2.0, 2.5, 3.0, 3.5, 3.8):
        s_, t_ = tractions(k, Sv)
        r_ = np.abs(t_) / s_
        sweep.append(dict(k=k, ratio_max=float(r_.max()), ratio_median=float(np.median(r_)),
                          share_above_0p6=float(np.mean(r_ > 0.6)), share_left_lateral=float(np.mean(t_ < 0)),
                          optimal_ratio=(k - 1) / (2 * np.sqrt(k))))

    cols = dict(chain=F['chain'], s_km=F['s'] / 1e3, x=F['x'], y=F['y'], lat=lat, lon=lon, conv_deg=gam,
                shmax_deg=theta, R=R, regime=n, mixed=mixed.astype(int), in_hull=hull.astype(int),
                sigma_n_MPa=sn, tau_MPa=tau, ratio=np.abs(tau) / sn,
                sigma_n_deck_MPa=sn0, tau_deck_MPa=tau0, vp=vp, vs=vs, rho=rho,
                Lb_m=Lb, Lb_deck_m=Lb_deck, Lb_cvm_deckstress_m=Lb_cvm_deckstress)
    keys = list(cols)
    np.savetxt(os.path.join(a.out, 'fault_tractions.csv'), np.c_[tuple(np.asarray(cols[c], float) for c in keys)],
               delimiter=',', header=','.join(keys), comments='', fmt='%.10g')
    np.savez(os.path.join(a.out, 'fault_tractions.npz'), **cols)

    def st(v):
        return [float(np.min(v)), float(np.percentile(v, 5)), float(np.median(v)), float(np.percentile(v, 95)),
                float(np.max(v))]
    per_chain = {}
    for c in range(5):
        m = F['chain'] == c
        per_chain[CHAIN_NAMES[c]] = dict(
            length_km=float(m.sum()), sigma_n_MPa=st(sn[m]), tau_MPa=st(tau[m]), ratio=st(np.abs(tau[m]) / sn[m]),
            left_lateral_share=float(np.mean(tau[m] < 0)), Lb_km=st(Lb[m] / 1e3),
            regime_share=[float(np.mean(n[m] == r)) for r in (0, 1, 2)])
    try:
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=HERE).stdout.strip()
        dirty = bool(subprocess.run(['git', 'status', '--porcelain'], capture_output=True, text=True, cwd=HERE).stdout.strip())
    except Exception:
        commit, dirty = 'unknown', True
    summ = dict(
        method='so2 Steps 1, 2, 9-11 at z = -10 km; no T_d0 edit', z_m=Z_S, Sv_eff_MPa=Sv, k=a.k,
        Dc_m=DC, b=B, n_points=int(len(lat)), commit=commit, uncommitted_tool=dirty,
        inputs=dict(mesh=os.path.relpath(MSH, SEAS), mesh_md5=md5(MSH), csm=CSM, csm_md5=md5(CSM),
                    sv_source=SV_SOURCE, cvm='UCVM cvmsi (CVM-S4.26.M01)'),
        all=dict(sigma_n_MPa=st(sn), tau_MPa=st(tau), ratio=st(np.abs(tau) / sn), left_lateral_share=float(np.mean(tau < 0)),
                 Lb_km=st(Lb / 1e3), Lb_deck_km=st(Lb_deck / 1e3), Lb_cvm_deckstress_km=st(Lb_cvm_deckstress / 1e3),
                 cells_per_Lb_500m=st(Lb / 500.0), regime_share=[float(np.mean(n == r)) for r in (0, 1, 2)],
                 mixed_share=float(mixed.mean()), outside_hull_share=float((~hull).mean()),
                 conv_deg=st(gam), shmax_deg=st(theta)),
        deck_today=dict(sigma_n_MPa=st(sn0), ratio=st(np.abs(tau0) / sn0)),
        per_chain=per_chain, k_sweep=sweep)
    json.dump(summ, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    print(json.dumps(dict(all=summ['all'], k_sweep=sweep), indent=1))


if __name__ == '__main__':
    main()
