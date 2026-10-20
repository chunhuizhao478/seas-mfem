"""Station plots of an LTS run against its global-stepping gold and the DRDG3D reference (item 7 LTS golds).
One PDF: an overview page (strike slip rate at every station), then one page per station with
(a)-(f) as in the gold plots (strike slip rate, front zoom, strike slip, strike shear stress, normal stress,
dip slip rate), (g)-(j) LTS - gold of strike slip rate, strike slip, strike shear stress and normal stress,
(k)-(l) the period-2 part 0.5 (y_i - (y_i-1 + y_i+1) / 2) of tau_h and sigma_n, RMS per 0.25 s window.
Arrival, peak and final slip in the titles come from measures() of compare_nitsche_fault_stations.py.
usage (from miniapps/seas): python3 -I jobs/vista_profiling/plot_lts_gold.py tpv205|tpv104 LTS_DIR GOLD_DIR REF_DIR
       OUT_STEM 'LTS label' 'gold label' 'ref label'
writes OUT_STEM_stations.pdf and OUT_STEM_pNN_<station>.png"""
import glob, importlib.util, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

prob, lts_dir, gold_dir, ref_dir, stem, lab_l, lab_g, lab_r = sys.argv[1:9]
spec = importlib.util.spec_from_file_location(
    'cmp', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts',
                        'compare_nitsche_fault_stations.py'))
cmp = importlib.util.module_from_spec(spec); spec.loader.exec_module(cmp)
C_G, C_L, C_R = '#1f4e8c', '#d1495b', '#222222'
T = 12.0


def files(d, ext):
    return {'x2_' + os.path.basename(f).rsplit('x2_', 1)[1][:-4]: f
            for f in glob.glob(d + '/*x2_*' + ext) if not f.endswith('_copies.dat')}


def key(s):
    w = s.split('_')
    return (float(w[3]), float(w[1]))


def nyq(y):
    return 0.5 * (y[1:-1] - 0.5 * (y[:-2] + y[2:]))


def win(t, y, w=0.25):
    c, r = [], []
    for a in np.arange(0, T, w):
        m = (t >= a) & (t < a + w)
        if m.any(): c.append(a + w / 2); r.append(np.sqrt(np.mean(y[m] ** 2)))
    return np.array(c), np.array(r)


SG, SL, SR = files(gold_dir, '.dat'), files(lts_dir, '.dat'), files(ref_dir, '.txt')
names = sorted(SG, key=key)
missing = [s for s in names if s not in SL]
if missing: sys.exit(f'LTS run lacks stations {missing}')
D = {}
for s in names:
    g, l = cmp.load(SG[s], 'nitsche'), cmp.load(SL[s], 'nitsche')
    r = cmp.load(SR[s], 'scec_mpa') if s in SR else None
    D[s] = (g, l, r)
pages = []
nc = 4; nr = int(np.ceil(len(names) / nc))
fig, axs = plt.subplots(nr, nc, figsize=(16, 3.0 * nr), squeeze=False)
for i, s in enumerate(names):
    g, l, r = D[s]; ax = axs[i // nc][i % nc]
    if r is not None: ax.plot(r['t'], r['h-slip-rate'], color=C_R, lw=0.8, label=lab_r)
    ax.plot(g['t'], g['h-slip-rate'], color=C_G, lw=1.6, label=lab_g)
    ax.plot(l['t'], l['h-slip-rate'], color=C_L, lw=0.8, ls='--', label=lab_l)
    ax.set_title(s, fontsize=10); ax.set_xlim(0, T); ax.grid(alpha=0.3)
    if i % nc == 0: ax.set_ylabel('strike slip rate (m/s)')
for j in range(len(names), nr * nc): axs[j // nc][j % nc].axis('off')
axs[0][0].legend(fontsize=8, loc='upper right')
fig.suptitle(f'SCEC {prob.upper()}: strike slip rate at every station; {lab_l} (dashed) against {lab_g} '
             f'and {lab_r}', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.97)); pages.append(('overview', fig))

for s in names:
    g, l, r = D[s]
    n = min(len(g['t']), len(l['t']))
    if np.max(np.abs(g['t'][:n] - l['t'][:n])) > 2e-9: sys.exit(f'{s}: LTS and gold times differ')
    t = g['t'][:n]
    mg, ml = cmp.measures(g, T, 0.1, 0.5), cmp.measures(l, T, 0.1, 0.5)
    ta = mg['arrival'] if np.isfinite(mg['arrival']) else 0.0
    fig, axs = plt.subplots(4, 3, figsize=(16, 15))
    A = axs.ravel()
    top = [('h-slip-rate', 1, '(a) strike slip rate (m/s)', (0, T)),
           ('h-slip-rate', 1, f'(b) strike slip rate (m/s), gold arrival {ta:.3f} s -0.5 / +2.5 s', (ta - 0.5, ta + 2.5)),
           ('h-slip', 1, '(c) strike slip (m)', (0, T)),
           ('h-shear-stress', 1, '(d) strike shear stress (MPa)', (0, T)),
           ('n-stress', 1, '(e) normal stress, compression positive (MPa)', (0, T)),
           ('v-slip-rate', 1, '(f) dip slip rate (m/s)', (0, T))]
    for ax, (k, _, title, xl) in zip(A[:6], top):
        if r is not None and k in r: ax.plot(r['t'], r[k], color=C_R, lw=0.8, label=lab_r)
        ax.plot(g['t'], g[k], color=C_G, lw=1.6, label=lab_g)
        ax.plot(l['t'], l[k], color=C_L, lw=0.8, ls='--', label=lab_l)
        ax.set_title(title, fontsize=10); ax.set_xlim(*xl); ax.grid(alpha=0.3)
    if r is None or 'n-stress' not in r: A[4].set_title('(e) normal stress, compression positive (MPa); no reference', fontsize=10)
    diffs = [('h-slip-rate', '(g) LTS - gold, strike slip rate (m/s)'), ('h-slip', '(h) LTS - gold, strike slip (m)'),
             ('h-shear-stress', '(i) LTS - gold, strike shear stress (MPa)'), ('n-stress', '(j) LTS - gold, normal stress (MPa)')]
    for ax, (k, title) in zip(A[6:10], diffs):
        ax.plot(t, l[k][:n] - g[k][:n], color=C_L, lw=0.7)
        ax.set_title(title, fontsize=10); ax.set_xlim(0, T); ax.grid(alpha=0.3)
    for ax, (k, title) in zip(A[10:12], (('h-shear-stress', '(k) period-2 part of strike shear stress (MPa)'),
                                          ('n-stress', '(l) period-2 part of normal stress (MPa)'))):
        for d, col, lab, lw in ((g, C_G, lab_g, 1.6), (l, C_L, lab_l, 0.9)):
            x, y = win(d['t'][1:-1], nyq(d[k]))
            ax.semilogy(x, np.maximum(y, 1e-12), color=col, lw=lw, label=lab)
        ax.set_title(title + ', RMS per 0.25 s', fontsize=10); ax.set_xlim(0, T); ax.grid(alpha=0.3)
    for ax in A[9:]: ax.set_xlabel('time (s)')
    A[0].legend(fontsize=8); A[10].legend(fontsize=8)
    da = (ml['arrival'] - mg['arrival']) * 1e3
    dp = (ml['peak'] / mg['peak'] - 1) * 100
    ds = (ml['slip'] / mg['slip'] - 1) * 100 if mg['slip'] > 0 else float('nan')
    line2 = f'{lab_l} vs {lab_g}: arrival {da:+.2f} ms  peak {dp:+.3f} %  final slip {ds:+.3f} %'
    if r is not None:
        mr = cmp.measures(r, T, 0.1, 0.5)
        line2 += (f';  vs {lab_r}: arrival {ml["arrival"] - mr["arrival"]:+.3f} s  peak '
                  f'{(ml["peak"] / mr["peak"] - 1) * 100:+.1f} %  slip {(ml["slip"] / mr["slip"] - 1) * 100:+.2f} %')
    w = s.split('_')
    fig.suptitle(f'SCEC {prob.upper()} station {s} (strike {w[1]} km, depth {w[3]} km)\n{line2}', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); pages.append((s, fig))
with PdfPages(stem + '_stations.pdf') as pdf:
    for i, (s, fig) in enumerate(pages):
        pdf.savefig(fig); fig.savefig(f'{stem}_p{i + 1:02d}_{s}.png', dpi=90); plt.close(fig)
print('wrote', stem + '_stations.pdf', len(pages), 'pages')
