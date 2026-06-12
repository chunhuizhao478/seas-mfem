#!/usr/bin/env python3
"""Retag a gmsh v2.2 .msh from MFEM SAFS physical tags to SeisSol's gmsh
boundary-condition convention (docs/gmsh.rst): Physical Surface 101=free
surface, 103=dynamic rupture, 105=absorbing (pumgen subtracts 100 -> BC 1/3/5).

MFEM SAFS tags -> SeisSol gmsh tags:
    101 "fault"  -> 103  (dynamic rupture)
    102 "top"    -> 101  (free surface)
    103 "bottom" -> 105  (absorbing)
    104 "sides"  -> 105  (absorbing)   [bottom+sides merge into 105]
      1 "rock"   ->   1  (volume / material group, unchanged)

Only physical (first) tags of surface/volume elements are remapped; geometry,
connectivity, node ordering are byte-preserved. This is NOT a mesh converter —
pumgen still does the .msh -> .puml.h5 conversion. It only aligns the tag
numbering to SeisSol's documented convention.
"""
import sys

# old physical tag -> new physical tag (applied atomically in one pass)
TAG_MAP = {101: 103, 102: 101, 103: 105, 104: 105}  # 1 (rock) untouched

# new $PhysicalNames block (dim, tag, name)
NEW_PHYSICAL_NAMES = [
    (2, 101, "free_surface"),
    (2, 103, "fault"),
    (2, 105, "absorbing"),
    (3, 1,   "rock"),
]


def main(inp, outp):
    with open(inp) as f:
        lines = f.readlines()

    out = []
    i = 0
    n = len(lines)
    counts = {}  # new_tag -> element count (surfaces only)
    while i < n:
        line = lines[i]
        if line.strip() == "$PhysicalNames":
            # replace the whole block
            out.append("$PhysicalNames\n")
            out.append(f"{len(NEW_PHYSICAL_NAMES)}\n")
            for dim, tag, name in NEW_PHYSICAL_NAMES:
                out.append(f'{dim} {tag} "{name}"\n')
            out.append("$EndPhysicalNames\n")
            # skip original block
            while lines[i].strip() != "$EndPhysicalNames":
                i += 1
            i += 1
            continue
        if line.strip() == "$Elements":
            out.append(line)               # $Elements
            i += 1
            count_line = lines[i]
            out.append(count_line)         # number of elements
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
                    out.append(lines[i])
                i += 1
            out.append(lines[i])           # $EndElements
            i += 1
            continue
        out.append(line)
        i += 1

    with open(outp, "w") as f:
        f.writelines(out)

    print("retagged surface-element counts by NEW tag (pumgen BC = tag-100):")
    for t in sorted(counts):
        bc = t - 100
        kind = {1: "free surface", 3: "dynamic rupture", 5: "absorbing"}.get(bc, "?")
        print(f"   tag {t} -> BC {bc} ({kind}): {counts[t]} triangles")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: retag_msh_to_seissol.py <in.msh> <out.msh>")
    main(sys.argv[1], sys.argv[2])
