#!/usr/bin/env python3
"""check_puml_faces.py -- PUML structural invariant, memory-safe (int32 keys).

Asserts, without any .msh round trip:
  P1 every face tagged BC=3 (dynamic rupture) is INTERIOR, i.e. shared by
     exactly 2 tets, and both tets carry the tag  -> count(BC3) == 2 * n_fault;
  P2 every face tagged BC=1/5 is EXTERIOR, i.e. appears exactly once;
  P3 every untagged (code 0) face is interior (shared by 2) -- no hole;
  P4 0 inverted tets.
Node ids fit int32 here (<2.1e9), so the face table is 3x cheaper than int64.
"""
import sys, numpy as np, h5py
FACE=[(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
p=sys.argv[1]
with h5py.File(p) as f:
    V=f["geometry"][:].astype(np.float64); T=f["connect"][:].astype(np.int64); B=f["boundary"][:].astype(np.int64)
n=len(T); assert T.max() < 2**31-1
print(f"tets {n:,}  nodes {len(V):,}")
Fa=np.empty((4*n,3),np.int32); code=np.empty(4*n,np.int8)
for k in range(4):
    Fa[k*n:(k+1)*n]=T[:,list(FACE[k])].astype(np.int32)
    code[k*n:(k+1)*n]=((B>>(8*k))&0xFF).astype(np.int8)
Fa.sort(axis=1)
sv=Fa.view([('a','i4'),('b','i4'),('c','i4')]).ravel()
o=np.argsort(sv,kind="stable"); sv=sv[o]; code=code[o]
first=np.ones(len(sv),bool); first[1:]=sv[1:]!=sv[:-1]
st=np.flatnonzero(first); cnt=np.diff(np.append(st,len(sv)))
print(f"distinct faces {len(st):,}  (shared x2 {int((cnt==2).sum()):,}, exposed x1 {int((cnt==1).sum()):,}, >2 {int((cnt>2).sum()):,})")
ok=True
if (cnt>2).any(): print("  !! FAIL a face is used by >2 tets"); ok=False
# per-distinct-face: the set of codes carried by its incident half-faces
mx=np.maximum.reduceat(code,st); mn=np.minimum.reduceat(code,st)
for nm,c in (("BC3 fault",3),("BC1 free-surface",1),("BC5 absorbing",5)):
    sel=(mx==c)
    n_int=int((sel&(cnt==2)).sum()); n_ext=int((sel&(cnt==1)).sum())
    if c==3:
        good=(n_ext==0)
        print(f"  P1 {nm:17s} distinct {int(sel.sum()):>9,}  interior(x2) {n_int:>9,}  exterior(x1) {n_ext:>7,}  [{'PASS' if good else 'FAIL'}]")
        # both sides tagged?
        both=int(((mn==c)&(mx==c)&(cnt==2)).sum())
        g2=(both==n_int); ok&=good and g2
        print(f"     both incident tets tagged: {both:,} of {n_int:,}  [{'PASS' if g2 else 'FAIL'}]")
    else:
        good=(n_int==0); ok&=good
        print(f"  P2 {nm:17s} distinct {int(sel.sum()):>9,}  exterior(x1) {n_ext:>9,}  interior(x2) {n_int:>7,}  [{'PASS' if good else 'FAIL'}]")
untag=(mx==0)
g3=not bool((untag&(cnt==1)).any()); ok&=g3
print(f"  P3 untagged faces exposed (holes): {int((untag&(cnt==1)).sum()):,}  [{'PASS' if g3 else 'FAIL'}]")
del Fa,sv,o,code
ninv=0
for s in range(0,n,6_000_000):
    q=V[T[s:s+6_000_000]]; d=q[:,1:]-q[:,0:1]
    ninv+=int((np.einsum('ij,ij->i',np.cross(d[:,0],d[:,1]),d[:,2])<=0).sum())
g4=(ninv==0); ok&=g4
print(f"  P4 inverted/zero-volume tets: {ninv}  [{'PASS' if g4 else 'FAIL'}]")
print(f"\nPUML STRUCTURE: {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
