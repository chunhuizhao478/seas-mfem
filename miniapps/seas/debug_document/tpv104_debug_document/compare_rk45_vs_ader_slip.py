#!/usr/bin/env python3
"""Quantify the TPV102/TPV104 strike-slip deficit: RK45 vs gold-ADER vs SCEC benchmark.

Companion to tpv102_tpv104_rk45_slip_deficit_findings_2026-06-02.md.
No deps beyond numpy.  Reproduces the §2 tables.

Column layouts (IMPORTANT):
  * TPV104 our .dat  : t h-slip h-slip-rate ...           -> strike = col idx 1 (SCEC layout)
  * TPV102 our .dat  : time slip1(dip) slip2(strike) ...  -> strike = col idx 2 (internal layout)
  * all benchmark .txt: t h-slip ...                      -> strike = col idx 1
"""
import numpy as np, os

BASE = '/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas'
DL   = '/Users/chunhuizhao/Downloads/seas-mfem'   # user's downloaded RK45 station data

def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#') or s[0].isalpha():
                continue
            try:
                rows.append([float(x) for x in s.split()])
            except ValueError:
                continue
    return np.array(rows)

def at(a, tq, vcol):       return np.interp(tq, a[:, 0], a[:, vcol])
def tend(a):               return a[-1, 0]


def tpv104():
    g = f'{BASE}/tpv104/gold'
    sts = ['x2_0_x3_7.5', 'x2_0_x3_3', 'x2_-9_x3_7.5', 'x2_0_x3_12', 'x2_-12_x3_3']
    print('\n================= TPV104 (strike = h-slip, col idx 1) =================')
    hdr = ('station', 'tc', 'RK45', 'ADER-O2mfadj', 'ADER-O2sym', 'ADER-O3sym', 'bench')
    print('{:14s}{:>6s}{:>9s}{:>13s}{:>12s}{:>12s}{:>8s}'.format(*hdr))
    for st in sts:
        P = {
            'rk':  f'{DL}/tpv104_spatial/combined_7762891_7763034/tpv104_station_{st}.dat',
            'a2a': f'{g}/results_nonsymmetric_mesh_mixedflux_adjacent_O2_job7680545/results/tpv104_mfadj_p1_O2_station_{st}.dat',
            'a2s': f'{g}/results_symmetric_mesh_O2_job7678882/results/tpv104_hybrid_O2_station_{st}.dat',
            'a3s': f'{g}/results_symmetric_mesh_O3_job7678716/results/tpv104_hybrid_O3_station_{st}.dat',
            'dr':  f'{BASE}/tpv104/benchmark_data/DRDG3D/tpv104_drdg3d_{st}.txt',
            'se':  f'{BASE}/tpv104/benchmark_data/seisol/tpv104_seisol_{st}.txt',
        }
        d = {k: load(v) for k, v in P.items()}
        tc = min(min(tend(x) for x in d.values()), 12.0)
        bm = 0.5 * (at(d['dr'], tc, 1) + at(d['se'], tc, 1))
        pct = lambda k: (at(d[k], tc, 1) - bm) / bm * 100
        print('{:14s}{:6.2f}{:+8.1f}%{:+12.1f}%{:+11.1f}%{:+11.1f}%{:8.2f}'.format(
            st, tc, pct('rk'), pct('a2a'), pct('a2s'), pct('a3s'), bm))


def tpv102():
    g = f'{BASE}/tpv102/gold/results_mixed_flux_adjacent_p1_O2_dev_job7681589/results'
    # flt_<a>_<d>  <->  x2_<a>_x3_<d>  (n -> -)
    sts = [('flt_0_7.5', 'x2_0_x3_7.5'), ('flt_0_3', 'x2_0_x3_3'), ('flt_0_12', 'x2_0_x3_12'),
           ('flt_9_7.5', 'x2_9_x3_7.5'), ('flt_n9_7.5', 'x2_-9_x3_7.5'), ('flt_12_3', 'x2_12_x3_3'),
           ('flt_n12_3', 'x2_-12_x3_3'), ('flt_12_12', 'x2_12_x3_12'), ('flt_n12_12', 'x2_-12_x3_12')]
    print('\n================= TPV102 (our strike = slip2, col idx 2) =================')
    print('{:12s}{:>6s}{:>9s}{:>13s}{:>8s}'.format('station', 'tc', 'RK45', 'ADER-O2mfadj', 'bench'))
    rk_all, ad_all = [], []
    for flt, x2 in sts:
        rk = load(f'{DL}/tpv102_spatial/combined_7762877_7763025/tpv102_station_{flt}.dat')
        ad = load(f'{g}/tpv102_mfadj_p1_O2_station_{flt}.dat')
        dr = load(f'{BASE}/tpv102/benchmark_data/scec_drdg3d/tpv102_drdg3d_{x2}.txt')
        py = load(f'{BASE}/tpv102/benchmark_data/scec_pylith/tpv102_pylith_{x2}.txt')
        tc = min(tend(rk), tend(ad), tend(dr), tend(py), 12.0)
        bm = 0.5 * (at(dr, tc, 1) + at(py, tc, 1))          # benchmark strike = h-slip col 1
        prk = (at(rk, tc, 2) - bm) / bm * 100               # our strike = slip2 col 2
        pad = (at(ad, tc, 2) - bm) / bm * 100
        rk_all.append(prk); ad_all.append(pad)
        print('{:12s}{:6.2f}{:+8.1f}%{:+12.1f}%{:8.2f}'.format(flt, tc, prk, pad, bm))
    print('{:12s}{:>6s}{:+8.1f}%{:+12.1f}%'.format('MEAN', '', np.mean(rk_all), np.mean(ad_all)))


if __name__ == '__main__':
    if not os.path.isdir(DL):
        print(f'NOTE: {DL} not found — point DL at the downloaded RK45 station data.')
    tpv104()
    tpv102()
