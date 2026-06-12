#!/usr/bin/env python3
"""Report quality metrics for a tet mesh (MEDIT/.mesh, gmsh/.msh, or VTU).

Usage:
    python check_msh_quality.py <mesh_file>

Reports:
  - shortest tet edge
  - q_tet histogram (isoperimetric ratio used by msh_to_vtu.py)
  - per-region (subdomain) tet count
  - acceptance check vs Phase 3.5b: min q_tet >= 0.05, fraction(q_tet>=0.3) >= 0.995
"""
import sys
import numpy as np
import meshio


def tet_iso_q(p):
    """Isoperimetric tet quality: q = 12 * (3*V)^(2/3) / sum(face_areas)
    Equivalent to msh_to_vtu.py's metric.  q=1 is regular, q->0 degenerate."""
    a, b, c, d = p[0], p[1], p[2], p[3]
    V = abs(np.dot(b-a, np.cross(c-a, d-a))) / 6.0
    if V < 1e-30: return 0.0
    A = (np.linalg.norm(np.cross(b-a, c-a))
         + np.linalg.norm(np.cross(b-a, d-a))
         + np.linalg.norm(np.cross(c-a, d-a))
         + np.linalg.norm(np.cross(c-b, d-b))) / 2.0
    if A <= 0: return 0.0
    # Standard scaled isoperimetric: 12 * (3V)^(2/3) / A, normalized so reg tet = 1.
    return (12.0 * (3.0 * V) ** (2.0 / 3.0)) / A / (3.0 ** (1.0/3.0) * 6.0 ** (2.0/3.0)) * 6.0 ** (2.0/3.0) * 3.0 ** (1.0/3.0) / 6.0 ** (2.0/3.0)


def tet_min_edge(p):
    edges = [np.linalg.norm(p[i]-p[j]) for i in range(4) for j in range(i+1, 4)]
    return min(edges)


def main():
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <mesh_file>", file=sys.stderr)
        return 2
    path = sys.argv[1]
    print(f"loading {path}")
    msh = meshio.read(path)
    pts = msh.points
    print(f"  vertices: {len(pts)}")
    if 'tetra' not in msh.cells_dict:
        print("  no tetrahedra in mesh", file=sys.stderr)
        return 3
    tets = msh.cells_dict['tetra']
    print(f"  tetrahedra: {len(tets)}")
    # Optional region tags (gmsh:physical or medit:ref).
    tags = None
    for k in ('medit:ref', 'gmsh:physical', 'cell_tags'):
        if hasattr(msh, 'cell_data') and msh.cell_data:
            for ct, arr in msh.cell_data.items():
                if 'tetra' in msh.cells_dict and ct in (k,):
                    pass
    region_arr = None
    if 'medit:ref' in (msh.cell_data or {}):
        d = msh.cell_data['medit:ref']
        # find the entry corresponding to tetra cell-block
        for cb, ra in zip(msh.cells, d):
            if cb.type == 'tetra':
                region_arr = ra
                break
    print()
    print("computing quality + edge stats ...")
    qs = np.zeros(len(tets))
    es = np.zeros(len(tets))
    for i, t in enumerate(tets):
        p = pts[t]
        qs[i] = tet_iso_q(p)
        es[i] = tet_min_edge(p)

    print(f"  edge_min_global = {es.min():.4g} m")
    print(f"  edge_p1         = {np.percentile(es, 1):.4g} m")
    print(f"  edge_med        = {np.percentile(es, 50):.4g} m")
    print(f"  edge_max        = {es.max():.4g} m")
    print()
    print(f"  q_tet_min   = {qs.min():.4g}")
    print(f"  q_tet_p1    = {np.percentile(qs, 1):.4g}")
    print(f"  q_tet_med   = {np.percentile(qs, 50):.4g}")
    print(f"  q_tet_max   = {qs.max():.4g}")
    print(f"  fraction q_tet >= 0.30 : {(qs >= 0.30).mean():.4f}")
    print(f"  fraction q_tet >= 0.05 : {(qs >= 0.05).mean():.4f}")
    print(f"  fraction q_tet <  0.05 : {(qs <  0.05).mean():.4f}")
    print()
    if region_arr is not None:
        from collections import Counter
        c = Counter(region_arr.tolist())
        print("  per-region tet count:")
        for k, v in sorted(c.items()):
            print(f"    region {k}: {v} tets")
    print()
    print("=== Phase 3.5b acceptance ===")
    pass1 = qs.min() >= 0.05
    pass2 = (qs >= 0.30).mean() >= 0.995
    pass3 = es.min() >= 100.0
    print(f"  min q_tet >= 0.05:         {pass1}  (got {qs.min():.4g})")
    print(f"  fraction q_tet >= 0.3:     {pass2}  (got {(qs>=0.30).mean():.4f})")
    print(f"  shortest tet edge >= 100m: {pass3}  (got {es.min():.4g} m)")
    print(f"  OVERALL: {'PASS' if (pass1 and pass2 and pass3) else 'FAIL'}")
    return 0 if (pass1 and pass2 and pass3) else 1


if __name__ == "__main__":
    sys.exit(main())
