#!/usr/bin/env python3
# Slow-parts study S1 (Chunhui Zhao, 2026-10-25): the ranks s1_*.sbatch records with perf, from a run's
# <prefix>_rank_balance.csv.  Prints them on stdout (one line, space-separated) and why on stderr:
#   slow       the rank with the most compute time (non-wait time less the record reduction, ph_reduce)
#   slow-other the rank with the most compute time on another node
#   typical    on each of those two nodes, the rank whose compute time is closest to the run's median
# Usage: python3 -I pick_perf_ranks.py <prefix>_rank_balance.csv
import csv
import sys


def main(path):
    with open(path) as f:
        rows = list(csv.DictReader(l for l in f if not l.startswith('#')))
    if not rows:
        sys.exit('pick_perf_ranks: no rows in %s' % path)
    red = 'ph_reduce' in rows[0]
    for r in rows:
        r['compute'] = float(r['nonwait_sec']) - (float(r['ph_reduce']) if red else 0.0)
    med = sorted(r['compute'] for r in rows)[len(rows) // 2]
    slow = max(rows, key=lambda r: r['compute'])
    other = [r for r in rows if r['node'] != slow['node']]
    picks = [('slow', slow)]
    if other:
        picks.append(('slow-other', max(other, key=lambda r: r['compute'])))
    for _, p in list(picks):
        same = [r for r in rows if r['node'] == p['node'] and r is not p and all(r is not q for _, q in picks)]
        if same:
            picks.append(('typical', min(same, key=lambda r: abs(r['compute'] - med))))
    sys.stderr.write('perf ranks from %s (compute = non-wait%s; median %.2f s):\n'
                     % (path, ' - ph_reduce' if red else '', med))
    for why, r in picks:
        sys.stderr.write('  r%s %s on %s: compute %.2f s, non-wait %s s\n'
                         % (r['rank'], why, r['node'], r['compute'], r['nonwait_sec']))
    print(' '.join(r['rank'] for _, r in picks))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('usage: pick_perf_ranks.py <prefix>_rank_balance.csv')
    main(sys.argv[1])
