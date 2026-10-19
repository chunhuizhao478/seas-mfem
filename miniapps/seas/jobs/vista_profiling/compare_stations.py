"""Compare the fault station files of two Nitsche runs (Chunhui Zhao, 2026-10-17).

usage: python3 -I compare_stations.py A B [--tmax T]

A and B are folders holding <prob>_station_*.dat (columns: t slip1 slip2 V1 V2 tau1 tau2
sigma_n slip_path).  Rows are matched by line; with --tmax only the rows with t <= T + 1e-9
are used, so a short run can be compared with the start of a longer one.  The times must
agree.  For each station and quantity (slip = columns 1-2, V = 3-4, tau = 5-6, sigma_n = 7)
it reports max |a - b| / max |a| over the rows, and prints the largest over all stations,
plus the number of stations whose rows are bitwise equal as printed.  Standard library only.
"""
import glob
import os
import sys

GROUPS = (('slip', (1, 2)), ('V', (3, 4)), ('tau', (5, 6)), ('sigma_n', (7,)))


def load(path, tmax):
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        vals = line.split()
        t = float(vals[0])
        if tmax is not None and t > tmax + 1e-9:
            break
        rows.append(vals)
    return rows


def main():
    args = sys.argv[1:]
    tmax = None
    if '--tmax' in args:
        i = args.index('--tmax')
        tmax = float(args[i + 1])
        del args[i:i + 2]
    a_dir, b_dir = args
    files = sorted(f for f in glob.glob(os.path.join(a_dir, '*_station_*.dat'))
                   if not f.endswith('_copies.dat'))
    if not files:
        sys.exit(f'no station files in {a_dir}')
    worst, where, same, n = 0.0, None, 0, 0
    for fa in files:
        fb = os.path.join(b_dir, os.path.basename(fa))
        if not os.path.exists(fb):
            sys.exit(f'missing {fb}')
        A, B = load(fa, tmax), load(fb, tmax)
        if len(A) != len(B):
            sys.exit(f'{os.path.basename(fa)}: {len(A)} rows against {len(B)}')
        if any(float(x[0]) != float(y[0]) for x, y in zip(A, B)):
            sys.exit(f'{os.path.basename(fa)}: times differ')
        same += A == B
        for name, cols in GROUPS:
            scale = max(abs(float(r[c])) for r in A for c in cols) or 1e-300
            d = max(abs(float(x[c]) - float(y[c])) for x, y in zip(A, B) for c in cols) / scale
            if d > worst:
                worst, where = d, f'{os.path.basename(fa)} {name}'
        n += 1
    rows = len(load(files[0], tmax))
    print(f'{b_dir} vs {a_dir}: {n} stations, {rows} rows each; {same} bitwise equal as printed; '
          f'largest difference {worst:.3e}' + (f' ({where})' if where else ''))


if __name__ == '__main__':
    main()
