#!/usr/bin/env python3
"""Sweep tetgen switches on a closed box around real ALT fault geometry.
Mirrors localbox_fill.py exactly, but parameterises switches and reports Joe-Liu eta."""
import sys, time
import numpy as np, tetgen
from scipy.spatial import cKDTree

F = np.load("build_tmp/fault_surface.npz"); P, T = F["P"], F["T"]
H = float(sys.argv[1]); Z0, Z1 = float(sys.argv[2]), float(sys.argv[3])
cx, cy = float(sys.argv[4]), float(sys.argv[5])
SWITCHES = sys.argv[6:]

lo = np.array([cx-H, cy-H, Z0]); hi = np.array([cx+H, cy+H, Z1])
ins = np.all((P > lo+1.0) & (P < hi-1.0), axis=1)
Tk = T[ins[T].all(1)]
if len(Tk) == 0: sys.exit("no fault triangles in box")
vid = np.unique(Tk); rm = -np.ones(len(P), np.int64); rm[vid] = np.arange(len(vid))
FP, FT = P[vid], rm[Tk]

V = FP[FT]
e = np.stack([np.linalg.norm(V[:,1]-V[:,0],axis=1), np.linalg.norm(V[:,2]-V[:,1],axis=1),
              np.linalg.norm(V[:,0]-V[:,2],axis=1)], 1)
S = (e**2).sum(1); A = 0.5*np.linalg.norm(np.cross(V[:,1]-V[:,0], V[:,2]-V[:,0]), axis=1)
ceil_ = 6.0*np.cbrt(A*np.sqrt(2*S)/3.0)**2/S
print(f"[patch] {len(FT):,} fault tris  min edge {e.min():.3f} m  "
      f"eta-ceiling min {ceil_.min():.5f}  tris with ceiling<0.2: {int((ceil_<0.2).sum())}")

X, Y, Zb = [lo[0],hi[0]], [lo[1],hi[1]], [lo[2],hi[2]]
BP = np.array([[X[i],Y[j],Zb[k]] for i in (0,1) for j in (0,1) for k in (0,1)], float)
idx = lambda i,j,k: i*4+j*2+k
BT = []
for k in (0,1): BT += [[idx(0,0,k),idx(1,0,k),idx(1,1,k)],[idx(0,0,k),idx(1,1,k),idx(0,1,k)]]
for j in (0,1): BT += [[idx(0,j,0),idx(1,j,0),idx(1,j,1)],[idx(0,j,0),idx(1,j,1),idx(0,j,1)]]
for i in (0,1): BT += [[idx(i,0,0),idx(i,1,0),idx(i,1,1)],[idx(i,0,0),idx(i,1,1),idx(i,0,1)]]
AP = np.vstack([BP, np.asarray(FP)]); AT = np.vstack([np.array(BT,np.int64), FT+len(BP)])
NF = len(FT)
print(f"[plc] {len(AT):,} facets, {len(AP):,} vertices\n")
print(f"{'switches':<22}{'tets':>12}{'Steiner':>10}{'eta_min':>10}{'<0.05':>9}{'<0.1':>9}"
      f"{'median':>9}{'edgemin':>9}{'faultlost':>10}{'sec':>7}")
shift = AP.mean(0)
for sw in SWITCHES:
    t0 = time.time()
    try:
        tg = tetgen.TetGen(np.ascontiguousarray(AP-shift), np.ascontiguousarray(AT))
        tg.tetrahedralize(switches=sw)
        TP = np.asarray(tg.node)+shift; TT = np.asarray(tg.elem)
    except Exception as ex:
        print(f"{sw:<22}  FAIL {str(ex).strip()[:60]}"); continue
    Vt = TP[TT]
    ee = np.stack([np.linalg.norm(Vt[:,a]-Vt[:,b],axis=1) for a,b in
                   [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]], 1)
    vol = np.abs(np.einsum('ij,ij->i', Vt[:,1]-Vt[:,0],
                 np.cross(Vt[:,2]-Vt[:,0], Vt[:,3]-Vt[:,0])))/6.0
    eta = 12.0*np.cbrt(3.0*vol)**2/np.maximum((ee**2).sum(1), 1e-300)
    _, mp = cKDTree(TP).query(AP, k=1)
    faces = np.sort(np.concatenate([TT[:,[0,1,2]],TT[:,[0,1,3]],TT[:,[0,2,3]],TT[:,[1,2,3]]]),axis=1)
    fset = set(map(tuple, np.unique(faces,axis=0).tolist()))
    ftri = np.sort(mp[AT[-NF:]], axis=1)
    lost = NF - sum(1 for t_ in map(tuple, ftri.tolist()) if t_ in fset)
    print(f"{sw:<22}{len(TT):>12,}{len(TP)-len(AP):>10,}{eta.min():>10.5f}"
          f"{int((eta<0.05).sum()):>9,}{int((eta<0.1).sum()):>9,}{np.median(eta):>9.4f}"
          f"{ee.min():>9.3f}{lost:>10}{time.time()-t0:>7.1f}")
