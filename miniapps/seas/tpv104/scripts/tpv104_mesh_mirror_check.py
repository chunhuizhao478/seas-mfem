#!/usr/bin/env python3
"""Check whether tpv104_1000m.msh is mirror-symmetric across x=0 and y=0.

The σ_n perturbation in the t12 dev run shows mirror-asymmetry across
x2=0 (e.g., x2=-12,x3=3 → +3.7 MPa, x2=+12,x3=3 → -4.2 MPa).  D1+D2/D3
ruled out the rotation pipeline as the source.  This script tests the
remaining hypothesis: is the mesh itself non-mirror-symmetric?

For mirror symmetry across x=0 to hold:
  - For every node (x, y, z) with |x| > tol, there must exist a mirror
    node at (-x, y, z) within position tolerance.
  - For every tet element with vertices {v1, v2, v3, v4}, there must
    exist a mirror tet whose vertices are the mirrored vertices of the
    first.

Same for y=0.

If non-symmetric, the mesh-generation step (gmsh on tpv104_*.geo) is
introducing the asymmetry that propagates into the σ_n perturbation we
observed.
"""

import argparse
import os
import sys

import numpy as np

DEFAULT_MESH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "mesh", "tpv104_1000m.msh"
)
TOL_POS = 1e-6  # meters; mesh coords are in meters

def parse_msh22(path):
    """Parse Gmsh 2.2 ASCII format. Returns (nodes_xyz, elements_tet)."""
    nodes = {}
    elements = []
    with open(path) as f:
        section = None
        for line in f:
            line = line.strip()
            if line.startswith("$"):
                if line == "$Nodes":          section = "nodes"; nstate = 0; continue
                if line == "$EndNodes":       section = None; continue
                if line == "$Elements":       section = "elements"; estate = 0; continue
                if line == "$EndElements":    section = None; continue
                continue
            if section == "nodes":
                if nstate == 0:
                    nstate = int(line)
                    continue
                parts = line.split()
                node_id = int(parts[0])
                xyz = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
                nodes[node_id] = xyz
            elif section == "elements":
                if estate == 0:
                    estate = int(line)
                    continue
                parts = line.split()
                # gmsh element line:
                # elem-id  elem-type  num-tags  tag1 ... tagN  v1 v2 ... vK
                etype = int(parts[1])
                ntags = int(parts[2])
                if etype == 4:  # 4-node tetrahedron
                    verts = list(map(int, parts[3+ntags:3+ntags+4]))
                    elements.append(verts)
                elif etype == 2:  # 3-node triangle (boundary face)
                    # Tags: physgroup, geomgroup, [partition]
                    physgroup = int(parts[3]) if ntags >= 1 else 0
                    verts = list(map(int, parts[3+ntags:3+ntags+3]))
                    if not hasattr(parse_msh22, "_tris"):
                        parse_msh22._tris = []
                    parse_msh22._tris.append((physgroup, verts))
    parse_msh22.tris = getattr(parse_msh22, "_tris", [])
    return nodes, elements


def check_mirror_nodes(nodes_xyz, axis, tol=TOL_POS):
    """For each node, check if its mirror across `axis` (0=x, 1=y, 2=z) exists."""
    coords = np.array(list(nodes_xyz.values()))
    n = len(coords)
    # Build a KDTree-free approach: sort by mirror coord, then bisect.
    # For mirror across axis a: point p mirrored is p with p[a] negated.
    mirrored = coords.copy()
    mirrored[:, axis] = -mirrored[:, axis]

    # Use a hashable rounding for fast lookup.
    hash_orig = {tuple(np.round(c / tol).astype(int)): i for i, c in enumerate(coords)}
    matches = 0
    unmatched_examples = []
    for i, m in enumerate(mirrored):
        key = tuple(np.round(m / tol).astype(int))
        if key in hash_orig:
            matches += 1
        else:
            if len(unmatched_examples) < 5:
                unmatched_examples.append(coords[i])
    return matches, n, unmatched_examples


def check_mirror_elements(nodes, elements, axis, tol=TOL_POS):
    """For each tet, check if a mirror tet exists with mirror-paired nodes."""
    coords = nodes  # dict: id → xyz
    # Build node-id mirror map: for each node, find its mirror node (or self).
    node_ids = list(coords.keys())
    node_xyz = np.array([coords[i] for i in node_ids])
    hash_orig = {tuple(np.round(c / tol).astype(int)): nid
                 for nid, c in zip(node_ids, node_xyz)}
    mirror_node = {}
    for nid in node_ids:
        c = coords[nid].copy()
        c[axis] = -c[axis]
        key = tuple(np.round(c / tol).astype(int))
        mirror_node[nid] = hash_orig.get(key, None)

    # For each element, build the mirrored vertex set (sorted) and look it up
    # in the original element set.
    elem_set = {tuple(sorted(e)) for e in elements}
    matches = 0
    unmatched_examples = []
    no_mirror_node = 0
    for e in elements:
        m = [mirror_node.get(v) for v in e]
        if any(x is None for x in m):
            no_mirror_node += 1
            continue
        key = tuple(sorted(m))
        if key in elem_set:
            matches += 1
        else:
            if len(unmatched_examples) < 3:
                unmatched_examples.append(
                    (e, m, [coords[v] for v in e]))
    return matches, len(elements), no_mirror_node, unmatched_examples


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mesh", nargs="?", default=DEFAULT_MESH,
        help=f"Path to .msh file (default: {DEFAULT_MESH})")
    args = parser.parse_args()
    print(f"Parsing {args.mesh} ...")
    nodes, elements = parse_msh22(args.mesh)
    print(f"Nodes: {len(nodes)}")
    print(f"Tetrahedral elements: {len(elements)}")

    for axis_name, axis in [("x (along-strike)", 0), ("y (fault-normal)", 1)]:
        print(f"\n=== Mirror check across {axis_name} = 0 ===")
        nmatch, ntot, unmatched = check_mirror_nodes(nodes, axis)
        print(f"  Nodes matching mirror partner: {nmatch} / {ntot} "
              f"({100.0*nmatch/ntot:.2f}%)")
        if nmatch < ntot:
            print(f"  Unmatched node samples (no mirror partner found):")
            for u in unmatched:
                print(f"    {u}")

        ematch, etot, no_mirror_node, unmatched_e = check_mirror_elements(
            nodes, elements, axis)
        print(f"  Elements matching mirror tet: {ematch} / {etot} "
              f"({100.0*ematch/etot:.2f}%)")
        if no_mirror_node:
            print(f"  Elements skipped (one or more vertices had no mirror "
                  f"node): {no_mirror_node}")
        if unmatched_e and len(unmatched_e) > 0:
            print(f"  Unmatched element samples (mirror-of-vertices set NOT "
                  f"a tet in mesh):")
            for orig, mir, xyz in unmatched_e[:2]:
                print(f"    orig vertices = {orig}")
                print(f"      coords      = {[list(c) for c in xyz]}")
                print(f"    mirror vert ids = {mir}")

    print("\n=== Fault-surface (Physical Surface 3) mirror check across x = 0 ===")
    fault_tris = [v for (g, v) in parse_msh22.tris if g == 3]
    print(f"  Fault triangles: {len(fault_tris)}")
    # Build mirror-node map across x.
    coords = nodes
    node_ids = list(coords.keys())
    hash_orig = {tuple(np.round(coords[i] / TOL_POS).astype(int)): i
                 for i in node_ids}
    mirror_node = {}
    for nid in node_ids:
        c = coords[nid].copy()
        c[0] = -c[0]
        key = tuple(np.round(c / TOL_POS).astype(int))
        mirror_node[nid] = hash_orig.get(key, None)

    fault_set = {tuple(sorted(t)) for t in fault_tris}
    matches = 0
    for t in fault_tris:
        m = [mirror_node.get(v) for v in t]
        if any(x is None for x in m):
            continue
        if tuple(sorted(m)) in fault_set:
            matches += 1
    print(f"  Fault triangles with mirror partner across x = 0: "
          f"{matches} / {len(fault_tris)} "
          f"({100.0*matches/max(len(fault_tris),1):.2f}%)")

    # Also check at the on-fault nodes specifically.
    fault_node_ids = set()
    for t in fault_tris:
        fault_node_ids.update(t)
    fault_node_ids = sorted(fault_node_ids)
    matches_n = 0
    for nid in fault_node_ids:
        c = coords[nid].copy()
        c[0] = -c[0]
        key = tuple(np.round(c / TOL_POS).astype(int))
        if key in hash_orig:
            matches_n += 1
    print(f"  Fault NODES with mirror partner across x = 0: "
          f"{matches_n} / {len(fault_node_ids)} "
          f"({100.0*matches_n/max(len(fault_node_ids),1):.2f}%)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
