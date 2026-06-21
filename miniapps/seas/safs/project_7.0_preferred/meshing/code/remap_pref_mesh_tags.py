#!/usr/bin/env python3
"""Remap safv4_deep_500m_opt.msh boundary tags to the SAFS/ALT convention so the
MFEM seas_spatial_dyn_driver hardcoded fallback (fault=101, natural={102},
absorbing={103,104}) consumes it with NO config/C++ change.

PREFERRED .msh tags (Gmsh v2.2, no $PhysicalNames):
  fault   = 101 + 102 + 103  (3 interior 2-sided sheets; SeisSol DR code 3 merged)
  free    = 201
  bottom  = 202
  sides   = 203
Remap (applied to BOTH gmsh tags of each type-2 boundary triangle):
  101,102,103 -> 101 (fault)   201 -> 102 (top)   202 -> 103 (bottom)   203 -> 104 (sides)
Tets (type 4, tag 1=rock) unchanged.  A $PhysicalNames block matching ALT is added.
"""
import sys

SRC = ("/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/"
       "safs/project_7.0_preferred/meshing/results/safv4_deep_500m_opt.msh")
DST = ("/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/"
       "safs/project_7.0_preferred/meshing/results/safv4_deep_500m_opt_safstags.msh")

TAGMAP = {101: 101, 102: 101, 103: 101, 201: 102, 202: 103, 203: 104}

PHYSNAMES = (
    "$PhysicalNames\n"
    "5\n"
    '2 101 "fault"\n'
    '2 102 "top"\n'
    '2 103 "bottom"\n'
    '2 104 "sides"\n'
    '3 1 "rock"\n'
    "$EndPhysicalNames\n"
)

tally = {}
in_elements = False
n_elem_seen = 0
declared_n_elem = None
wrote_physnames = False

with open(SRC) as fin, open(DST, "w") as fout:
    prev_was_endmeshformat = False
    for line in fin:
        s = line.rstrip("\n")

        # Insert $PhysicalNames right after $EndMeshFormat (PREFERRED has none).
        if prev_was_endmeshformat:
            fout.write(PHYSNAMES)
            wrote_physnames = True
            prev_was_endmeshformat = False

        if s == "$EndMeshFormat":
            fout.write(line)
            prev_was_endmeshformat = True
            continue

        if s == "$Elements":
            in_elements = True
            fout.write(line)
            continue
        if s == "$EndElements":
            in_elements = False
            fout.write(line)
            continue

        if in_elements:
            f = s.split()
            # First line inside $Elements is the count.
            if declared_n_elem is None and len(f) == 1:
                declared_n_elem = int(f[0])
                fout.write(line)
                continue
            n_elem_seen += 1
            etype = int(f[1])
            if etype == 2:                       # boundary triangle: remap tags
                ntags = int(f[2])
                # tags occupy f[3 : 3+ntags]; remap any that are in TAGMAP
                for k in range(3, 3 + ntags):
                    old = int(f[k])
                    f[k] = str(TAGMAP.get(old, old))
                new_phys = int(f[3])
                tally[new_phys] = tally.get(new_phys, 0) + 1
                fout.write(" ".join(f) + "\n")
            else:                                # tet (type 4) etc.: verbatim
                tally[("type", etype)] = tally.get(("type", etype), 0) + 1
                fout.write(line)
            continue

        fout.write(line)

# ---- validation ----
print("declared $Elements count :", declared_n_elem)
print("elements processed       :", n_elem_seen)
print("$PhysicalNames inserted  :", wrote_physnames)
print("post-remap boundary-tri tag tally:")
for t in sorted(k for k in tally if isinstance(k, int)):
    print(f"  tag {t}: {tally[t]}")
print("volume element tally:")
for k in tally:
    if isinstance(k, tuple):
        print(f"  type {k[1]}: {tally[k]}")

EXPECT = {101: 132333, 102: 175033, 103: 86583, 104: 36010}
ok = all(tally.get(t) == c for t, c in EXPECT.items())
ok = ok and (n_elem_seen == declared_n_elem) and wrote_physnames
print("VALIDATION:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
