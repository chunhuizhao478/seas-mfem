#!/usr/bin/env python3
"""Retag the multi-strand ``safv4_deep_500m.msh`` (gmsh v2.2) from the
preferred-project physical tags to SeisSol's gmsh boundary-condition convention,
merging all three fault strands into ONE dynamic-rupture group.

This is the multi-strand sibling of ``retag_msh_to_seissol.py``.  The original
script targets the OLD single-fault mesh (tags 101 fault / 102 top / 103 bottom /
104 sides) and would mistag this mesh (e.g. Banning->free surface); it must NOT
be used here.

New-mesh tags (RUNLOG_safv4_remesh_phases0-5_2026-06-12.md line 264):
    101 San Andreas  ┐
    102 Banning      ├-> 103  (dynamic rupture; pumgen BC 3)   [single group]
    103 Garnet Hill  ┘
    201 top / DEM    --> 101  (free surface;     pumgen BC 1)
    202 bottom       ┐
    203 sides        ┴-> 105  (absorbing;        pumgen BC 5)
      1 rock (vol)   --> 1    (material/volume group, unchanged)

Collision-safety: 101 and 103 are BOTH sources and 103 is also a target, and 201
maps to 101.  This is safe because each element's ORIGINAL first tag is read once
and the mapped value written once -- there is never a sequential re-map.

Geometry/connectivity/node ordering are byte-preserved; only physical (first)
tags of surface elements are remapped.  The input has NO $PhysicalNames block, so
a SeisSol-convention block is INSERTED after $EndMeshFormat (gmsh requires
$PhysicalNames before $Nodes).  This is NOT a mesh converter -- pumgen still does
the .msh -> .puml.h5 conversion.

Usage: retag_safv4_deep_to_seissol.py <in.msh> <out.msh>
"""
import sys

# original physical tag -> new physical tag (applied atomically in one pass)
TAG_MAP = {101: 103, 102: 103, 103: 103, 201: 101, 202: 105, 203: 105}
# 1 (rock volume) intentionally untouched.

# SeisSol-convention $PhysicalNames block (dim, tag, name)
NEW_PHYSICAL_NAMES = [
    (2, 101, "free_surface"),
    (2, 103, "fault"),
    (2, 105, "absorbing"),
    (3, 1,   "rock"),
]

GMSH_TRI = 2  # 3-node triangle


def _physical_names_block():
    block = ["$PhysicalNames\n", f"{len(NEW_PHYSICAL_NAMES)}\n"]
    for dim, tag, name in NEW_PHYSICAL_NAMES:
        block.append(f'{dim} {tag} "{name}"\n')
    block.append("$EndPhysicalNames\n")
    return block


def main(inp, outp):
    with open(inp) as f:
        lines = f.readlines()

    has_physical = any(l.strip() == "$PhysicalNames" for l in lines)

    out = []
    i = 0
    n = len(lines)
    counts = {}      # new_tag -> remapped surface-element count
    unmapped = {}    # old_tag -> count of surface elements left unchanged

    while i < n:
        line = lines[i]
        s = line.strip()

        if s == "$EndMeshFormat":
            out.append(line)
            i += 1
            if not has_physical:
                out.extend(_physical_names_block())
            continue

        if s == "$PhysicalNames":
            # replace an existing block (future inputs); skip the original
            out.extend(_physical_names_block())
            while lines[i].strip() != "$EndPhysicalNames":
                i += 1
            i += 1
            continue

        if s == "$Elements":
            out.append(line)               # $Elements
            i += 1
            out.append(lines[i])           # element count
            i += 1
            while lines[i].strip() != "$EndElements":
                tok = lines[i].split()
                # gmsh v2.2: id type ntags tag1 tag2 ... nodes...
                ntags = int(tok[2])
                if ntags >= 1:
                    old = int(tok[3])
                    if old in TAG_MAP:
                        new = TAG_MAP[old]
                        tok[3] = str(new)
                        counts[new] = counts.get(new, 0) + 1
                        out.append(" ".join(tok) + "\n")
                    else:
                        # not in map (e.g. volume tag 1) -> byte-preserve
                        if int(tok[1]) == GMSH_TRI:
                            unmapped[old] = unmapped.get(old, 0) + 1
                        out.append(lines[i])
                else:
                    out.append(lines[i])
                i += 1
            out.append(lines[i])           # $EndElements
            i += 1
            continue

        out.append(line)
        i += 1

    with open(outp, "w") as f:
        f.writelines(out)

    print(f"inserted $PhysicalNames block: {not has_physical}")
    print("retagged surface-element counts by NEW tag (pumgen BC = tag-100):")
    kind = {1: "free surface", 3: "dynamic rupture", 5: "absorbing"}
    for t in sorted(counts):
        bc = t - 100
        print(f"   tag {t} -> BC {bc} ({kind.get(bc, '?')}): "
              f"{counts[t]} triangles")
    if unmapped:
        print("WARNING: surface elements with tags not in TAG_MAP "
              "(left unchanged):", dict(sorted(unmapped.items())))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: retag_safv4_deep_to_seissol.py <in.msh> <out.msh>")
    main(sys.argv[1], sys.argv[2])
