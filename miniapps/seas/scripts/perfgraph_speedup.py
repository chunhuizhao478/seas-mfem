#!/usr/bin/env python3
"""perfgraph_speedup.py — diff two Caliper runtime-report "perfgraph" dumps and
report per-region speedup, for measuring the ADER hot-path optimization
(document/caliper_perfgraph_dev/ader_hotpath_optimization_plan_*.md) end-to-end
on Frontera.

The C++ micro-benchmark `seas_bench_ader_hotpath` measures the kernels locally;
this tool measures the REAL run: capture a perfgraph BEFORE and AFTER the
optimization (same job, same mesh/np/tfinal) and compare.

INPUT: the text region tree printed by
   CALI_CONFIG="runtime-report(calc.inclusive=true)"
i.e. lines of the form
   <indented region path>   <min/rank> <max/rank> <avg/rank> <time %>
(exactly what the user pastes from the *.err).  Lines that do not end in four
numeric columns (headers, "negligible", "(0.04 total)") are ignored.

USAGE:
   ./perfgraph_speedup.py BASELINE.txt OPTIMIZED.txt
   ./perfgraph_speedup.py --metric max BASELINE.txt OPTIMIZED.txt   # use max/rank
                                                                     # (wall-time proxy)

Regions are keyed by their LEAF name (last "::" component duplicates — e.g.
ApplySpatialDerivative appears under both CK recursions — are SUMMED, so the row
is the total time in that region across the step).  Pure stdlib.
"""
import argparse
import re
import sys

# A data row ends in 4 numeric columns: min, max, avg, pct.
_NUM = r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'
_ROW = re.compile(r'^(?P<name>.*\S)\s+'
                  r'(?P<min>' + _NUM + r')\s+'
                  r'(?P<max>' + _NUM + r')\s+'
                  r'(?P<avg>' + _NUM + r')\s+'
                  r'(?P<pct>' + _NUM + r')\s*$')


def parse(path, metric):
    """Return {region_name: summed_time} for the chosen metric column."""
    totals = {}
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line.strip():
                continue
            m = _ROW.match(line)
            if not m:
                continue                       # header / non-data line
            name = m.group('name').strip()
            # Drop a leading tree glyph if present, keep the full seas::... path.
            name = name.lstrip('|>+- ').strip()
            try:
                val = float(m.group(metric))
            except ValueError:
                continue
            totals[name] = totals.get(name, 0.0) + val
    if not totals:
        sys.exit('ERROR: no data rows parsed from %s '
                 '(is it a runtime-report dump?)' % path)
    return totals


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('baseline')
    ap.add_argument('optimized')
    ap.add_argument('--metric', choices=['avg', 'max', 'min'], default='avg',
                    help='which per-rank column to compare (default: avg; '
                         'use "max" as a wall-time proxy for a synchronous run)')
    ap.add_argument('--sort', choices=['baseline', 'speedup', 'delta'],
                    default='baseline', help='row ordering (default: baseline time)')
    args = ap.parse_args()

    base = parse(args.baseline, args.metric)
    opt = parse(args.optimized, args.metric)

    names = set(base) | set(opt)
    rows = []
    for n in names:
        b = base.get(n)
        o = opt.get(n)
        speedup = (b / o) if (b is not None and o not in (None, 0.0)) else None
        delta = (b - o) if (b is not None and o is not None) else None
        rows.append((n, b, o, speedup, delta))

    keymap = {
        'baseline': lambda r: -(r[1] or 0.0),
        'speedup':  lambda r: -(r[3] or 0.0),
        'delta':    lambda r: -(r[4] or 0.0),
    }
    rows.sort(key=keymap[args.sort])

    def fnum(x):
        return '       -' if x is None else '%8.2f' % x

    print('# perfgraph speedup  (metric = %s/rank)' % args.metric)
    print('# baseline : %s' % args.baseline)
    print('# optimized: %s' % args.optimized)
    print('%-52s %10s %10s %8s %10s' %
          ('region', 'baseline', 'optimized', 'speedup', 'delta(s)'))
    print('-' * 94)
    for n, b, o, s, d in rows:
        sp = '   -   ' if s is None else '%6.2fx' % s
        disp = n if len(n) <= 52 else ('…' + n[-51:])
        print('%-52s %10s %10s %8s %10s' %
              (disp, fnum(b), fnum(o), sp, fnum(d)))

    # Headline: the top-level "step" region (total wall) if present.
    step = next((n for n in names if n.endswith('::step')), None)
    if step and base.get(step) and opt.get(step):
        print('-' * 94)
        print('TOTAL step  %.2f -> %.2f s   speedup %.2fx   (-%.1f%%)' % (
            base[step], opt[step], base[step] / opt[step],
            100.0 * (1.0 - opt[step] / base[step])))


if __name__ == '__main__':
    main()
