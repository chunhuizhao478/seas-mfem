#!/usr/bin/env python3
"""medit_to_puml.py -- mmg MEDIT output -> SeisSol PUML/HDF5, without a .msh hop.

Writing a 70.6 M-tet gmsh22 ASCII intermediate and reading it back with meshio
costs ~3 GB of disk and a lot of wall time for nothing, so this goes straight from
MEDIT to PUML using the conventions msh_to_puml.py reverse-engineered:

  datasets : geometry (Nnode,3) f8 | connect (Ntet,4) u8 0-based
             boundary (Ntet,) i4    | group (Ntet,) i4 = 1
  attrs    : boundary-format='i32', topology-format='geometric'
  packing  : boundary = sum_i code_i << (8*i), face i -> local verts
             f0={0,2,1} f1={0,1,3} f2={1,2,3} f3={0,3,2}
  tag->BC  : 101 fault -> 3, 102 top -> 1, 103/104 -> 5 absorbing

Tets are re-oriented to positive signed volume (SeisSol blows up on inverted ones).
Face keys are 3-field structured records, never bit-packed ints -- at 12 M vertices
three 24-bit fields overflow int64 and silently swap BC codes between faces.
"""
import numpy as np, pandas as pd, h5py
import sys
# Paths are arguments, not constants: this runs in both trees and the ALT output
# must not be written under a PREF-named file.  Defaults keep the original call
# working.
F   = sys.argv[1] if len(sys.argv) > 1 else "build_tmp/s1.mmg_out.mesh"
OUT = sys.argv[2] if len(sys.argv) > 2 else "build_tmp/shakeoutD_alt_heavy_s1.puml.h5"
# MEASURED from the file, never hardcoded. These were PREFERRED's byte layout
# (LV,LT,LS = 7, 13428231, 84344339 / NV,NT,NS = 12042113, 70641218, 3049525);
# on any other mesh those line numbers land mid-vertex-block and pandas parses
# whatever is there WITHOUT error -- coordinates read as connectivity, plausible
# mesh out the far end. medit_hdr reproduces all six of those values exactly on
# PREFERRED's file, which is how it was validated.
from medit_hdr import medit_sections
_sec = medit_sections(F)
LV, NV = _sec["Vertices"]
LT, NT = _sec["Tetrahedra"]
LS, NS = _sec["Triangles"]
print(f"[hdr] Vertices {NV:,} @ {LV:,}   Tetrahedra {NT:,} @ {LT:,}   "
      f"Triangles {NS:,} @ {LS:,}")
BC = {101: 3, 102: 1, 103: 5, 104: 5}
FACEV = [(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
REC = np.dtype([('a',np.int32),('b',np.int32),('c',np.int32)])
def rec(A):
    return np.ascontiguousarray(np.sort(A,axis=1).astype(np.int32)).view(REC).ravel()

P = pd.read_csv(F, sep=r"\s+", header=None, skiprows=LV+1, nrows=NV,
                usecols=range(3), dtype=np.float64, engine="c").to_numpy()
S = pd.read_csv(F, sep=r"\s+", header=None, skiprows=LS+1, nrows=NS,
                usecols=range(4), dtype=np.int32, engine="c").to_numpy()
tri, ref = S[:,:3]-1, S[:,3]; del S
code = np.array([BC[int(r)] for r in ref], np.int8)
key = rec(tri); o = np.argsort(key, kind="stable"); key = key[o]; code = code[o]
# mmg lists at least one fault triangle TWICE (same tag, same vertices), so any
# converter keyed on a unique facet set asserts here unless it dedupes first.
# Measured on PREFERRED (its numbers, not ALT's): 3,049,525 listed -> 3,049,524
# unique, multiplicity 2, both tag 101; deduping gives 2,761,484 unique fault
# facets, reconciling exactly with the 4 PLC facets absent from the output.
# ALT's counts differ -- everything below is derived, so nothing here is lineage-
# specific except this provenance note.
key, first = np.unique(key, return_index=True)
code = code[first]
print(f"[dedup] surface {NS:,} listed -> {len(key):,} unique")
print(f"[in] {NV:,} verts  {NT:,} tets  {NS:,} surface tris")
for t, c in ((101,3),(102,1),(104,5)):
    print(f"   tag {t} -> BC {c}: {int((ref==t).sum()):,}")
del tri, ref, o

CH = 8_000_000
T_all = np.empty((NT,4), np.int32); B_all = np.zeros(NT, np.int32)
ninv = 0; nface = np.zeros(6, np.int64)
for s in range(0, NT, CH):
    n = min(CH, NT-s)
    T = pd.read_csv(F, sep=r"\s+", header=None, skiprows=LT+1+s, nrows=n,
                    usecols=range(4), dtype=np.int32, engine="c").to_numpy()-1
    V = P[T]
    sv = np.einsum('ij,ij->i', V[:,1]-V[:,0], np.cross(V[:,2]-V[:,0], V[:,3]-V[:,0]))
    neg = sv < 0
    if neg.any():
        ninv += int(neg.sum()); T[neg] = T[neg][:, [0,1,3,2]]
    b = np.zeros(n, np.int32)
    for k in range(4):
        fk = rec(T[:, list(FACEV[k])])
        pos = np.searchsorted(key, fk); np.clip(pos, 0, len(key)-1, out=pos)
        hit = key[pos] == fk
        if hit.any():
            cc = code[pos[hit]].astype(np.int32)
            b[hit] |= cc << (8*k)
            for u in np.unique(cc): nface[u] += int((cc==u).sum())
        del fk, pos, hit
    T_all[s:s+n] = T; B_all[s:s+n] = b
    del T, V, sv, neg, b
    print(f"  ..{s+n:,}/{NT:,}", flush=True)
print(f"[orient] {ninv:,} tets re-oriented to positive volume")
print(f"[boundary] BC 1 free surface {nface[1]:,}   BC 3 fault {nface[3]:,}   "
      f"BC 5 absorbing {nface[5]:,}")
# DERIVED, not hardcoded.  The original asserted 2 x 2,761,484 -- PREFERRED's fault
# count -- which is simply wrong for any other lineage (ALT has 2,564,474) and
# would fail a correct conversion.  The invariant that actually matters is
# self-consistency: the fault is an interior crack, so every unique fault triangle
# in the surface listing must appear as a face of exactly TWO tets.
NF_SURF = int((code == 3).sum())
assert nface[3] == 2 * NF_SURF, (
    f"fault faces {nface[3]:,} != 2 x {NF_SURF:,} unique fault triangles "
    f"-> {2*NF_SURF - nface[3]:,} missing (pinholes) or single-sided")
print(f"[check] fault interior-crack invariant holds: {nface[3]:,} = 2 x {NF_SURF:,}")
with h5py.File(OUT, "w") as f:
    f.create_dataset("geometry", data=P)
    f.create_dataset("connect", data=T_all.astype(np.uint64))
    f.create_dataset("boundary", data=B_all)
    f.create_dataset("group", data=np.ones(NT, np.int32))
    f.attrs["boundary-format"] = np.bytes_("i32")
    f.attrs["topology-format"] = np.bytes_("geometric")
print(f"[out] {OUT}")
