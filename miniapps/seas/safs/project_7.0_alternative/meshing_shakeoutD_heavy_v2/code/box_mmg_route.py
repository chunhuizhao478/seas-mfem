#!/usr/bin/env python3
"""Box: tetgen fill -> MEDIT (fault=RequiredTriangles) -> mmg3d -> eta.
Tests whether mmg can clear fault-face slivers with the fault frozen."""
import sys, time, subprocess, os
import numpy as np, tetgen
from scipy.spatial import cKDTree

SW   = sys.argv[1]
HMIN = float(sys.argv[2]); HMAX = float(sys.argv[3])
H,Z0,Z1,cx,cy = [float(x) for x in sys.argv[4:9]]
EXTRA= sys.argv[9:]
MMG="/Users/chunhuizhao/miniforge/envs/mmg/bin/mmg3d_O3"

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
NB=len(BT); NF=len(FT)
shift=AP.mean(0)
tg=tetgen.TetGen(np.ascontiguousarray(AP-shift),np.ascontiguousarray(AT))
tg.tetrahedralize(switches=SW)
TP=np.asarray(tg.node)+shift; TT=np.asarray(tg.elem)

def eta_of(Pp,Tt):
    V=Pp[Tt]
    e=np.stack([np.linalg.norm(V[:,a]-V[:,b],axis=1) for a,b in
                [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]],1)
    vol=np.abs(np.einsum('ij,ij->i',V[:,1]-V[:,0],np.cross(V[:,2]-V[:,0],V[:,3]-V[:,0])))/6.
    return 12.*np.cbrt(3.*vol)**2/np.maximum((e**2).sum(1),1e-300), e.min()

REC=np.dtype([('a',np.int32),('b',np.int32),('c',np.int32)])
rec=lambda A: np.ascontiguousarray(np.sort(A,axis=1).astype(np.int32)).view(REC).ravel()
_,mp=cKDTree(TP).query(AP,k=1)
plc=rec(mp[AT]); marks=np.concatenate([np.full(NB,104,np.int32),np.full(NF,101,np.int32)])
o=np.argsort(plc,kind='stable'); plc=plc[o]; marks=marks[o]
LOCAL=[(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
faces=[];tags=[]
for k in range(4):
    f=TT[:,list(LOCAL[k])]; fk=rec(f)
    pos=np.clip(np.searchsorted(plc,fk),0,len(plc)-1); hit=plc[pos]==fk
    if hit.any(): faces.append(f[hit]); tags.append(marks[pos[hit]])
faces=np.vstack(faces); gt=np.concatenate(tags)
keep=np.ones(len(faces),bool); isf=gt==101
fk=rec(faces[isf]); _,first=np.unique(fk,return_index=True)
idxf=np.flatnonzero(isf); keep[np.setdiff1d(idxf,idxf[first])]=False
faces,gt=faces[keep],gt[keep]
e0,_=eta_of(TP,TT)
print(f"[tetgen {SW}] {len(TT):,} tets  eta_min {e0.min():.5f}  <0.05 {int((e0<0.05).sum()):,}"
      f"  <0.1 {int((e0<0.1).sum()):,}  med {np.median(e0):.4f}   fault faces {int((gt==101).sum()):,}/{NF:,}")

with open("/tmp/box.mesh","w") as f:
    f.write("MeshVersionFormatted 2\nDimension 3\n\n")
    f.write(f"Vertices\n{len(TP)}\n")
    np.savetxt(f,np.column_stack([TP-shift,np.zeros(len(TP))]),fmt="%.15g %.15g %.15g %d")
    np.save("/tmp/shift.npy",shift)
    f.write(f"\nTriangles\n{len(faces)}\n")
    np.savetxt(f,np.column_stack([faces+1,gt]),fmt="%d")
    req=np.nonzero(gt==101)[0]+1
    f.write(f"\nRequiredTriangles\n{len(req)}\n"); np.savetxt(f,req.reshape(-1,1),fmt="%d")
    f.write(f"\nTetrahedra\n{len(TT)}\n")
    np.savetxt(f,np.column_stack([TT+1,np.ones(len(TT),int)]),fmt="%d")
    f.write("\nEnd\n")
cmd=[MMG,"-in","/tmp/box.mesh","-out","/tmp/box.o.mesh","-opnbdy","-hgrad","1.3",
     "-hausd","30","-hmin",str(HMIN),"-hmax",str(HMAX),"-v","1"]+EXTRA
t0=time.time(); r=subprocess.run(cmd,capture_output=True,text=True)
print(f"[mmg] {' '.join(cmd[5:])}   rc={r.returncode}  {time.time()-t0:.1f}s")
if r.returncode!=0:
    print(r.stdout[-600:]); sys.exit(1)
txt=open("/tmp/box.o.mesh").read().split()
def sect(name):
    i=txt.index(name); n=int(txt[i+1]); return i+2,n
i,nv=sect("Vertices"); QP=np.array(txt[i:i+4*nv],float).reshape(nv,4)[:,:3]
i,nt=sect("Tetrahedra"); QT=np.array(txt[i:i+5*nt],int).reshape(nt,5)[:,:4]-1
i,ntr=sect("Triangles"); QR=np.array(txt[i:i+4*ntr],int).reshape(ntr,4)
e1,emin=eta_of(QP,QT)
nf_out=int((QR[:,3]==101).sum())
print(f"[after mmg] {nt:,} tets  eta_min {e1.min():.5f}  <0.05 {int((e1<0.05).sum()):,}"
      f"  <0.1 {int((e1<0.1).sum()):,}  med {np.median(e1):.4f}  edge_min {emin:.3f}"
      f"  fault tris {nf_out:,} (in {int((gt==101).sum()):,})")
