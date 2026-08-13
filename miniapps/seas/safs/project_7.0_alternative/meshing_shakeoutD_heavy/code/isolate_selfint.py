#!/usr/bin/env python3
"""isolate_selfint.py -- WHICH part of the PLC self-intersects: fault, hull, or the pair?"""
import numpy as np, tetgen, sys
from scipy.spatial import cKDTree
H=np.load("build_tmp/hull.npz"); F=np.load("build_tmp/fault_surface.npz")
HP,HT,PART=H["P"],H["T"],H["PART"]; FP,FT=F["P"],F["T"]
tree=cKDTree(HP); d,k=tree.query(FP,k=1,distance_upper_bound=1e-3); hit=np.isfinite(d)
newid=np.empty(len(FP),np.int64); newid[hit]=k[hit]
newid[~hit]=len(HP)+np.arange(int((~hit).sum()))
P=np.vstack([HP,FP[~hit]]); T=np.vstack([HT,newid[FT]])
MARK=np.concatenate([PART,np.full(len(FT),7,np.int8)])
def test(name, keep):
    t=T[keep]
    v=np.unique(t); rm=-np.ones(len(P),np.int64); rm[v]=np.arange(len(v))
    p=P[v]-P[v].mean(0); t2=rm[t]
    try:
        tg=tetgen.TetGen(np.ascontiguousarray(p),np.ascontiguousarray(t2))
        tg.tetrahedralize(switches="d")
        print(f"  {name:22s} {len(t2):>10,} facets  ->  NO self-intersection")
    except RuntimeError as e:
        print(f"  {name:22s} {len(t2):>10,} facets  ->  {str(e).strip()[:60]}")
    except Exception as e:
        print(f"  {name:22s} {len(t2):>10,} facets  ->  {type(e).__name__}: {str(e)[:50]}")
print("[tetgen -d]")
test("fault only",      MARK==7)
test("hull only",       MARK!=7)
test("top only",        MARK==1)
test("hull + fault",    np.ones(len(T),bool))
