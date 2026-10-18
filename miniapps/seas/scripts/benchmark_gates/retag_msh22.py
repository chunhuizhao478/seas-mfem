#!/usr/bin/env python3
"""Gate-only copy of a Gmsh 2.2 ASCII mesh with physical tags renamed
(e.g. 101:1,103:3,105:5) so the PF test's fixed attributes (1, 3, 5) apply.
Geometry and element order are unchanged.  Used by the B0 PF gate only."""
import sys
src, dst, spec = sys.argv[1], sys.argv[2], sys.argv[3]
m = {int(a): int(b) for a, b in (kv.split(':') for kv in spec.split(','))}
out = open(dst, 'w')
sec = None
n_changed = 0
with open(src) as f:
    for line in f:
        s = line.strip()
        if s.startswith('$'):
            sec = s
            out.write(line); continue
        if sec == '$PhysicalNames':
            p = s.split()
            if len(p) >= 3 and int(p[1]) in m:
                p[1] = str(m[int(p[1])]); line = ' '.join(p) + '\n'
        elif sec == '$Elements':
            p = s.split()
            if len(p) > 4:
                ntags = int(p[2])
                if ntags >= 1 and int(p[3]) in m:
                    p[3] = str(m[int(p[3])]); n_changed += 1
                    line = ' '.join(p) + '\n'
        out.write(line)
out.close()
print(f"{dst}: {n_changed} element physical tags renamed via {m}")
