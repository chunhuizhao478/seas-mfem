#!/usr/bin/env python3
"""Minimal: PLC -> tetgen -> counts. No analysis, so peak RSS is tetgen's own."""
import sys, time, resource
import numpy as np, tetgen
SW=sys.argv[1]; H,Z0,Z1,cx,cy=[float(x) for x in sys.argv[2:7]]
F=np.load("build_tmp/fault_surface.npz"); P,T=F["P"],F["T"]
lo=np.array([cx-H,cy-H,Z0]); hi=np.array([cx+H,cy+H,Z1])
ins=np.all((P>lo+1.)&(P<hi-1.),axis=1); Tk=T[ins[T].all(1)]
vid=np.unique(Tk); rm=-np.ones(len(P),np.int64); rm[vid]=np.arange(len(vid))
FP,FT=P[vid],rm[Tk]
X,Y,Zb=[lo[0],hi[0]],[lo[1],hi[1]],[lo[2],hi[2]]
BP=np.array([[X[i],Y[j],Zb[k]] for i in(0,1) for j in(0,1) for k in(0,1)],float)
ix=lambda i,j,k:i*4+j*2+k
BT=[]
for k in(0,1): BT+=[[ix(0,0,k),ix(1,0,k),ix(1,1,k)],[ix(0,0,k),ix(1,1,k),ix(0,1,k)]]
for j in(0,1): BT+=[[ix(0,j,0),ix(1,j,0),ix(1,j,1)],[ix(0,j,0),ix(1,j,1),ix(0,j,1)]]
for i in(0,1): BT+=[[ix(i,0,0),ix(i,1,0),ix(i,1,1)],[ix(i,0,0),ix(i,1,1),ix(i,0,1)]]
AP=np.vstack([BP,FP]); AT=np.vstack([np.array(BT,np.int64),FT+len(BP)])
base=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
t0=time.time(); shift=AP.mean(0)
tg=tetgen.TetGen(np.ascontiguousarray(AP-shift),np.ascontiguousarray(AT))
tg.tetrahedralize(switches=SW)
TT=np.asarray(tg.elem); TP=np.asarray(tg.node)
pk=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(f"{SW:<12} fault_tris {len(FT):>9,}  tets {len(TT):>10,}  verts {len(TP):>9,}  "
      f"{time.time()-t0:>7.1f}s  peakRSS {pk/2**30:>6.2f} GB  pre {base/2**30:.2f} GB  "
      f"-> {pk/max(len(TT),1):.0f} B/tet")
