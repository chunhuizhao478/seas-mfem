"""Judge one LTS Vista run against the global-stepping gold it replaces (plan item 7; criteria of the
LTS buffer-fix review, criterion 4 with the 10 Pa floor Chunhui Zhao accepted on 2026-10-20):
  1. the SCEC pass/fail (compare_nitsche_fault_stations.py items 1-3) equals the gold's;
  2. per station: arrival within 5 ms, peak slip rate within 1 %, final slip within 0.5 %
     (peak where gold peak > 1e-3 m/s, slip where gold slip > 0.1 m);
  3. noise RMS (sigma_n, tau_h, tau_v) within 10 % of the gold's at every station;
  4. period-2 part (0.5 (y_i - (y_{i-1} + y_{i+1}) / 2), RMS per 1 s window) of tau_h and sigma_n of
     the worst station is within 2x the gold's worst station, or below 10 Pa, in every window.
Also prints, for information: per station period-2 ratio, and max |run - gold| / max |gold| of slip,
V and tau over the whole trace.
Run from miniapps/seas; writes OUT_PREFIX_lts.{txt,json} and OUT_PREFIX_gold.{txt,json} (the two
compare_nitsche_fault_stations.py runs) and prints the verdict.
usage: python3 -I jobs/vista_profiling/judge_lts_gold.py tpv205|tpv104 GOLD_RESULTS_DIR RUN_DIR OUT_PREFIX [OLD_RUN_DIR]"""
import glob, json, os, subprocess, sys
import numpy as np

prob, gold, run, out = sys.argv[1:5]
old = sys.argv[5] if len(sys.argv) > 5 else None
if prob == 'tpv205':
    common = ['--stations', 'tpv205/stations/tpv205_fault_stations.csv', '--prefix', 'tpv205',
              '--ref', 'drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4', '--peak-min-stations', '14']
else:
    common = ['--stations', 'tpv104/stations/tpv104_fault_stations.csv', '--prefix', 'tpv104',
              '--ref', 'drdg3d=tpv104/benchmark_data/DRDG3D']


def compare(label, path, extra, tag):
    cmd = ['python3', '-I', 'scripts/compare_nitsche_fault_stations.py'] + common + \
          ['--run', f'{label}={path}:nitsche'] + extra + ['--json', f'{out}_{tag}.json']
    r = subprocess.run(cmd, capture_output=True, text=True)
    open(f'{out}_{tag}.txt', 'w').write(r.stdout + r.stderr)
    if r.returncode:
        print(f'compare {tag} exit {r.returncode}:\n' + (r.stdout + r.stderr)[-2000:]); sys.exit(2)
    return json.load(open(f'{out}_{tag}.json'))


extra = ['--gold', f'gold={gold}:nitsche'] + (['--gold', f'lts_r1={old}:nitsche'] if old else [])
A = compare('lts', run, extra, 'lts')
B = compare('gold', gold, [], 'gold')
ok = True
def verdict(name, good, text):
    global ok
    ok &= good
    print(f'{"PASS" if good else "FAIL"}  {name}: {text}')

pa, pb = A['pass'], B['pass']
keys = ['1', '2', '3']
same = all(pa[k] == pb[k] for k in keys)
verdict('1 SCEC pass/fail as the gold', same,
        'lts ' + ' '.join(f'{k}:{pa[k]}' for k in keys) + f' (peak count {pa.get("3_count")}); gold ' +
        ' '.join(f'{k}:{pb[k]}' for k in keys) + f' (peak count {pb.get("3_count")})')

M = {}
for k, v in A['measures'].items():
    lab, st = k.split('|'); M.setdefault(lab, {})[st] = v
L, G = M['lts'], M['gold']
rows = []
worst = {'arr': (0, ''), 'peak': (0, ''), 'slip': (0, ''), 'noise': (0, '')}
for s in sorted(G):
    a, b = L.get(s), G[s]
    if a is None: verdict('station present', False, s); continue
    da = abs(a['arrival'] - b['arrival']) if np.isfinite(a['arrival']) and np.isfinite(b['arrival']) else \
        (0.0 if np.isfinite(a['arrival']) == np.isfinite(b['arrival']) else np.inf)
    dp = abs(a['peak'] / b['peak'] - 1) * 100 if b['peak'] > 1e-3 else 0.0
    ds = abs(a['slip'] / b['slip'] - 1) * 100 if b['slip'] > 0.1 else 0.0
    nr = []
    for c in ('noise_sn', 'noise_th', 'noise_tv'):
        if np.isfinite(a[c]) and np.isfinite(b[c]) and b[c] > 0: nr.append(a[c] / b[c])
    dn = max(abs(x - 1) for x in nr) * 100 if nr else 0.0
    for key, val in (('arr', da), ('peak', dp), ('slip', ds), ('noise', dn)):
        if val > worst[key][0]: worst[key] = (val, s)
    rows.append((s, da * 1e3, dp, ds, nr))
verdict('2a arrival within 5 ms', worst['arr'][0] <= 5e-3, f'max {worst["arr"][0]*1e3:.2f} ms ({worst["arr"][1]})')
verdict('2b peak slip rate within 1 %', worst['peak'][0] <= 1.0, f'max {worst["peak"][0]:.3f} % ({worst["peak"][1]})')
verdict('2c final slip within 0.5 %', worst['slip'][0] <= 0.5, f'max {worst["slip"][0]:.3f} % ({worst["slip"][1]})')
verdict('3 noise RMS within 10 %', worst['noise'][0] <= 10.0, f'max |ratio - 1| {worst["noise"][0]:.1f} % ({worst["noise"][1]})')

# trace level
def stations(d):
    return {os.path.basename(f).split('station_')[1][:-4]: f for f in glob.glob(d + '/*_station_*.dat')
            if not f.endswith('_copies.dat')}
def load(f):
    return np.loadtxt(f, comments='#')
def nyq(y):
    return 0.5 * (y[1:-1] - 0.5 * (y[:-2] + y[2:]))
def windows(t, n, T):
    out_ = []
    for a in range(T):
        m = (t >= a) & (t < a + 1)
        out_.append(np.sqrt(np.mean(n[m] ** 2)) if m.any() else np.nan)
    return np.array(out_)
SG, SL = stations(gold), stations(run)
T = 12
WG = {6: None, 7: None}; WL = {6: None, 7: None}
per = []
tr = {2: (0, ''), 4: (0, ''), 6: (0, '')}
for s in sorted(SG):
    g, r = load(SG[s]), load(SL[s])
    n = min(len(g), len(r))
    if np.max(np.abs(g[:n, 0] - r[:n, 0])) > 2e-9:
        verdict('times equal', False, s); continue
    t = g[1:n - 1, 0]
    rat = []
    for c in (6, 7):
        wg = windows(t, nyq(g[:n, c] / 1e6), T); wl = windows(t, nyq(r[:n, c] / 1e6), T)
        WG[c] = wg if WG[c] is None else np.fmax(WG[c], wg)
        WL[c] = wl if WL[c] is None else np.fmax(WL[c], wl)
        rat.append(np.nanmax(wl / wg))
    per.append((s, rat[0], rat[1]))
    for c in (2, 4, 6):
        e = np.max(np.abs(r[:n, c] - g[:n, c])) / (np.max(np.abs(g[:n, c])) or 1.0)
        if e > tr[c][0]: tr[c] = (e, s)
for c, lab in ((6, 'tau_h'), (7, 'sigma_n')):
    q = WL[c] / WG[c]
    good = (q <= 2.0) | (WL[c] < 1e-5)
    verdict(f'4 period-2 {lab} worst station within 2x gold or below 10 Pa per 1 s window', bool(good.all()),
            'ratio per window ' + ' '.join(f'{x:.2f}' for x in q) +
            (f'; above 2x but below 10 Pa in {int(np.sum(q > 2.0))} window(s)' if np.any(q > 2.0) else ''))
    print(f'      gold  (MPa) ' + ' '.join(f'{x:.2g}' for x in WG[c]))
    print(f'      lts   (MPa) ' + ' '.join(f'{x:.2g}' for x in WL[c]))
print('info  per station max period-2 ratio over windows (tau_h, sigma_n): ' +
      ', '.join(f'{s} {a:.2f}/{b:.2f}' for s, a, b in per))
print('info  max |lts - gold| / max |gold| over 0-12 s: ' +
      ', '.join(f'{lab} {tr[c][0]:.3g} ({tr[c][1]})' for c, lab in ((2, 'slip'), (4, 'V'), (6, 'tau_h'))))
print('info  per station |d arrival| ms, |d peak| %, |d slip| %, noise ratios sn/th/tv:')
for s, da, dp, ds, nr in rows:
    print(f'      {s:16s} {da:7.2f} {dp:7.3f} {ds:7.3f}  ' + '/'.join(f'{x:.2f}' for x in nr))
print('VERDICT', 'ALL PASS' if ok else 'FAIL')
