#!/usr/bin/env python3
"""diagnose_plc.py -- tetgen -d as the ORACLE for PLC self-intersections.

pymeshlab's self-intersection filter false-positives on benign non-manifold welds,
so this project uses tetgen's own detector and reads the offenders it writes to
`*_skipped.face` (+ `.node` -- note the .node file is the FULL node list, a past bug
was assuming it held only the offenders).
"""
import numpy as np, os, sys, glob
from scipy.spatial import cKDTree
H=np.load("build_tmp/hull.npz"); F=np.load("build_tmp/fault_surface.npz")
HP,HT,PART=H["P"],H["T"],H["PART"]; FP,FT=F["P"],F["T"]
tree=cKDTree(HP); d,k=tree.query(FP,k=1,distance_upper_bound=1e-3); hit=np.isfinite(d)
newid=np.empty(len(FP),np.int64); newid[hit]=k[hit]
newid[~hit]=len(HP)+np.arange(int((~hit).sum()))
P=np.vstack([HP,FP[~hit]]); T=np.vstack([HT,newid[FT]])
MARK=np.concatenate([PART,np.full(len(FT),7,np.int8)])
which = sys.argv[1] if len(sys.argv)>1 else "all"
if which=="fault":
    keep=MARK==7
elif which=="hull":
    keep=MARK!=7
else:
    keep=np.ones(len(T),bool)
T=T[keep]; MARK=MARK[keep]
v=np.unique(T); rm=-np.ones(len(P),np.int64); rm[v]=np.arange(len(v)); P=P[v]; T=rm[T]
print(f"[plc:{which}] {len(T):,} facets {len(P):,} verts")
os.makedirs("build_tmp/tgd",exist_ok=True)
base="build_tmp/tgd/plc"
with open(base+".poly","w") as f:
    f.write(f"{len(P)} 3 0 0\n")
    for i,p in enumerate(P): f.write(f"{i+1} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    f.write(f"{len(T)} 1\n")
    for t,m in zip(T,MARK):
        f.write(f"1 0 {int(m)}\n3 {t[0]+1} {t[1]+1} {t[2]+1}\n")
    f.write("0\n0\n")
print(f"[poly] {base}.poly written")
