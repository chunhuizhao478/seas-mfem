#!/usr/bin/env python3
"""Compare the station files of two Nitsche runs (rank invariance D-05/D-09, identity D-03).

  compare_station_runs.py <dirA> <dirB> <prefix> [--bytes] [--tol 1e-11]

Default: per station, per quantity group (slip = slip1/slip2, rate = V1/V2, shear = tau1/tau2,
sigma_n, state), max |a - b| over all rows divided by the group's max |value| in run A; passes when
every ratio <= tol and both files have the same times.  --bytes: the data rows (non-comment lines)
must be byte-identical; header comment lines are reported, not compared.
"""
import glob, os, sys
a, b, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
bytes_mode = '--bytes' in sys.argv
tol = float(sys.argv[sys.argv.index('--tol') + 1]) if '--tol' in sys.argv else 1e-11
fa = sorted(f for f in glob.glob(f"{a}/{prefix}_station_*.dat") if not f.endswith('_copies.dat'))
names = [os.path.basename(f) for f in fa]
fb = sorted(os.path.basename(f) for f in glob.glob(f"{b}/{prefix}_station_*.dat") if not f.endswith('_copies.dat'))
ok = names == fb and len(names) > 0
if names != fb: print("station file sets differ:", set(names) ^ set(fb))
groups = {'slip': [1, 2], 'rate': [3, 4], 'shear': [5, 6], 'sigma_n': [7], 'state': [8]}
worst = 0.0
for n in names:
    la = open(os.path.join(a, n)).read().splitlines()
    lb = open(os.path.join(b, n)).read().splitlines()
    da = [l for l in la if not l.startswith('#')]
    db = [l for l in lb if not l.startswith('#')]
    ha = [l for l in la if l.startswith('#')]
    hb = [l for l in lb if l.startswith('#')]
    if bytes_mode:
        same = da == db
        ok &= same
        print(f"{n}: {len(da)} data rows {'byte-identical' if same else 'DIFFER'}; header lines "
              f"{'identical' if ha == hb else 'differ'}")
        if not same:
            for i, (x, y) in enumerate(zip(da, db)):
                if x != y:
                    print("   first difference at data row", i, "\n   A:", x, "\n   B:", y); break
            if len(da) != len(db): print("   row counts", len(da), len(db))
        continue
    A = [[float(v) for v in l.split()] for l in da]
    B = [[float(v) for v in l.split()] for l in db]
    if len(A) != len(B) or any(x[0] != y[0] for x, y in zip(A, B)):
        print(f"{n}: times differ ({len(A)} vs {len(B)} rows)"); ok = False; continue
    out = []
    for g, cols in groups.items():
        sc = max(abs(r[c]) for r in A for c in cols)
        d = max(abs(x[c] - y[c]) for x, y in zip(A, B) for c in cols)
        rel = d / sc if sc > 0 else (0.0 if d == 0 else float('inf'))
        worst = max(worst, rel)
        out.append(f"{g} {rel:.2e}")
    print(f"{n}: " + ", ".join(out))
if not bytes_mode:
    print(f"worst relative difference {worst:.3e} (tolerance {tol:g})")
    ok &= worst <= tol
print("  PASSED" if ok else "  FAILED")
sys.exit(0 if ok else 1)
