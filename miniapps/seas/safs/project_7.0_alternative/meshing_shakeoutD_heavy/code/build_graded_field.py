#!/usr/bin/env python3
"""build_graded_field.py -- the g-graded 3-D target size field for LEB refinement.

This is what makes the base-plus-refine path work: gradation is a property of the
size FIELD, not of the mesher, so g = 0.15 does NOT need mmg to build the whole
mesh in one pass.  Pre-grade here, then refine the existing base with LEB
(memory-light, chunked) -- which closes the frequency gate in the SAME pass,
because the gate term is part of the field.

    h_raw = min( h_fault(d) , Vs_MUSCAL / gate , hmax )
    h     = min-plus transform  min_y [ h_raw(y) + g*|x-y| ]   (EUCLIDEAN)

Three things here are deliberate and were each measured, not assumed:

1. **Euclidean, not L1.**  A separable L1 sweep is wrong in the UNSAFE direction:
   d_L1 >= d_L2, so every term h+g*d is larger and the min is larger, i.e.
   h_L1 >= h_euclid -- L1 grades COARSER and under-refines.  A field g-Lipschitz
   in L1 is only g*(|u|_1/|u|_2)-Lipschitz in Euclidean, delivering 0.150 on axes
   but 0.212 on face diagonals and 0.260 on body diagonals: anisotropically
   failing the very spec this rebuild exists to satisfy.  Measured on a test grid:
   pure L1 max error vs exact Euclidean 61.9 %.

2. **Isotropic z for the gradation.**  With the natural non-uniform z
   (100/500/2000 m) a 26-neighbour relaxation still errs 26.1 % because the mask
   cannot represent directions when dz >> dxy.  Measured max error vs exact:

        grid            R=1 (26)   R=2 (98)   R=3 (290)
        anisotropic z     26.08 %    16.59 %    11.21 %
        isotropic z       11.66 %     4.46 %     2.23 %

   So gradation runs on a UNIFORM dz = dxy grid with an R=2 mask (4.46 %), and the
   sweep uses g/SAFETY so the DELIVERED gradation stays <= g in every direction.

3. **Each slab takes the MINIMUM of the field over its vertical extent**, sampled
   finely, so moving to 500 m dz costs no near-surface accuracy: the min is the
   conservative (finer) choice and is the same z-pooling rule the gate uses.
"""
import argparse, math, numpy as np, h5py, netCDF4 as ncdf
from scipy.spatial import cKDTree
from pyproj import Transformer

ap = argparse.ArgumentParser()
ap.add_argument("--mesh", help="PUML providing the fault (BC 3); ALT prefers --fault-npz")
ap.add_argument("--fault-npz", default="build_tmp/fault_surface.npz",
                help="deduplicated fault triangles; avoids re-reading the 5.5 GB parent "
                     "and exists before any PUML does")
ap.add_argument("--muscal", default="/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc")
ap.add_argument("--gate", type=float, default=0.8)
ap.add_argument("--hmax", type=float, default=5000.0)
ap.add_argument("--g", type=float, default=0.15)
ap.add_argument("--dxy", type=float, default=500.0)
ap.add_argument("--depth", type=float, default=80000.0)
ap.add_argument("--radius", type=int, default=2, help="relaxation mask radius")
ap.add_argument("--safety", type=float, default=1.05,
                help="sweep at g/safety to absorb the mask's residual error")
ap.add_argument("--subz", type=int, default=10, help="sub-samples per z slab")
ap.add_argument("--decim", type=int, default=21)
ap.add_argument("--no-fault", action="store_true", help="gate term only (costing)")
ap.add_argument("--out", default="build_tmp/target_g015.npz")
a = ap.parse_args()

Q = np.load("build_tmp/domain_corners_utm.npy")
X = np.arange(Q[:,0].min()-2000, Q[:,0].max()+2000+a.dxy, a.dxy)
Y = np.arange(Q[:,1].min()-2000, Q[:,1].max()+2000+a.dxy, a.dxy)
Z = -np.arange(0.0, a.depth+a.dxy, a.dxy)          # ISOTROPIC: dz = dxy
print(f"[grid] {len(X)} x {len(Y)} x {len(Z)} = {len(X)*len(Y)*len(Z):,} nodes "
      f"(isotropic {a.dxy:.0f} m)")

FACEV=[(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
tree = None
if not a.no_fault:
    # ALT takes the fault from build_tmp/fault_surface.npz, which already holds
    # exactly these triangles deduplicated.  The PUML path re-reads a 5.5 GB parent
    # (boundary alone is 1.07 GB) purely to average vertices into centroids, and on
    # this shared box that read competes with a running mmg.  The npz also exists
    # BEFORE the PUML does, which is the point: the graded field is needed for the
    # refinement pass that produces the PUML.
    if a.fault_npz:
        Fn=np.load(a.fault_npz); nfacet=len(Fn["T"])
        # decim is tuned against the PUML path, where every facet appears TWICE
        # (once per adjacent tet).  The npz is deduplicated, so halve the stride to
        # keep the same centroid density -- otherwise the distance field silently
        # thins by 2x.
        Pc=Fn["P"][Fn["T"]].mean(1)[::max(1, a.decim // 2)]
    else:
        with h5py.File(a.mesh) as f:
            V=f["geometry"][:]; C=f["connect"]; B=f["boundary"][:].astype(np.int64)
            tri=[]
            for s in range(0, C.shape[0], 8_000_000):
                c=C[s:s+8_000_000][:].astype(np.int64); b=B[s:s+8_000_000]
                for k in range(4):
                    m=((b>>(8*k))&0xFF)==3
                    if m.any(): tri.append(c[m][:,list(FACEV[k])])
            Fc=np.vstack(tri); Pc=V[Fc].mean(1)[::a.decim]
            nfacet=len(Fc)//2
    print(f"[fault] {nfacet:,} facets -> {len(Pc):,} decimated centroids")
    tree = cKDTree(Pc)
DCAP=25000.0
PD=np.array([0,125,250,500,1000,2000,4000,7500,15000],float)
PH=np.array([115,115,141,215,262,529,700,1000,1500],float)

m=ncdf.Dataset(a.muscal)
mlon=m["longitude"][:].data.astype(float); mlat=m["latitude"][:].data.astype(float)
mdep=m["depth"][:].data.astype(float); VS=np.asarray(m["vs"][:],np.float32)
fwd=Transformer.from_crs("EPSG:32611","EPSG:4326",always_xy=True)
def near(ax,v):
    i=np.clip(np.searchsorted(ax,v),0,len(ax)-1); im=np.maximum(i-1,0)
    return np.where(np.abs(ax[im]-v)<np.abs(ax[i]-v),im,i).astype(np.int32)

XX,YY=np.meshgrid(X,Y,indexing="ij"); xy=np.column_stack([XX.ravel(),YY.ravel()])
del XX,YY
lo,la=fwd.transform(xy[:,0],xy[:,1]); ii,jj=near(mlon,lo),near(mlat,la); del lo,la
H=np.empty((len(X),len(Y),len(Z)),np.float32)
half=a.dxy/2.0
for k,z in enumerate(Z):
    if tree is not None:
        d,_=tree.query(np.column_stack([xy,np.full(len(xy),z)]),k=1,
                       distance_upper_bound=DCAP,workers=-1)
        d=np.where(np.isfinite(d),d,DCAP)
        hf=np.interp(d,PD,PH,right=PH[-1])
        hf=np.where(d>15000,PH[-1]+(d-15000)*0.20,hf)
        hf=np.where(d>=DCAP-1.0,np.inf,hf)
    else:
        hf=np.full(len(xy),np.inf)
    # MINIMUM of Vs over the slab this level represents -> conservative (finer)
    vmin=None
    for zz in np.linspace(z-half, z+half, a.subz):
        v=VS[near(mdep,max(-zz,0.0)),jj,ii].astype(np.float64)
        v=np.where(np.isfinite(v),v,np.inf)
        vmin=v if vmin is None else np.minimum(vmin,v)
    H[:,:,k]=np.minimum(np.minimum(hf,vmin/a.gate),a.hmax).reshape(len(X),len(Y))
    if k%20==0: print(f"  ..z {z:8.0f} m", flush=True)
print(f"[raw] h {H.min():.0f} .. {H.max():.0f} m   median {np.median(H):.0f}")

gs = a.g / a.safety
R = a.radius
OFF=[(i,j,k) for i in range(-R,R+1) for j in range(-R,R+1) for k in range(-R,R+1)
     if (i,j,k)!=(0,0,0) and math.gcd(math.gcd(abs(i),abs(j)),abs(k))==1]
print(f"[grade] Euclidean min-plus, R={R} ({len(OFF)} offsets), "
      f"sweep g={gs:.5f} (= {a.g}/{a.safety})")
nx,ny,nz=H.shape
for it in range(400):
    before=float(H.sum())
    for di,dj,dk in OFF:
        sx=slice(max(di,0),nx+min(di,0)); dx_=slice(max(-di,0),nx+min(-di,0))
        sy=slice(max(dj,0),ny+min(dj,0)); dy_=slice(max(-dj,0),ny+min(-dj,0))
        sz=slice(max(dk,0),nz+min(dk,0)); dz_=slice(max(-dk,0),nz+min(-dk,0))
        L=np.float32(a.dxy*math.sqrt(di*di+dj*dj+dk*dk))
        np.minimum(H[dx_,dy_,dz_],H[sx,sy,sz]+np.float32(gs)*L,out=H[dx_,dy_,dz_])
    d=before-float(H.sum())
    if it%5==0 or d<=0: print(f"  it{it:03d} sum drop {d:.6e}",flush=True)
    if d<=0: break
print(f"[graded] h {H.min():.0f} .. {H.max():.0f} m   median {np.median(H):.0f}")
np.savez(a.out,X=X,Y=Y,Z=Z,H=H,gate=a.gate,g=a.g)
print(f"[out] {a.out}   {H.nbytes/2**30:.2f} GB")
