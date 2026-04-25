# Limitation 01 — Structured tet mesh for TPV104 dynamic rupture

Date introduced: 2026-04-25
Scope: TPV104 dynamic-rupture path only (`drivers/tpv104_driver.cpp`,
`tpv104/mesh/build_symmirror_mesh.py`, `jobs/tpv104/tpv104_symmirror_*`).
Quasi-dynamic SEAS (BP5) is NOT affected.

## What we currently do

`tpv104/mesh/build_symmirror_mesh.py` generates a Y-mirror-symmetric
tetrahedral mesh by:

1. Defining a **regular Cartesian hex grid** with `(NX, NY, NZ)` cells
   over a rectangular box.
2. Splitting each hex into **6 tetrahedra** using a fixed v0–v7 diagonal
   pattern (Pattern A on the +Y half) and the Y-mirror-image of that
   pattern (Pattern B = y-bit-swap of A) on the −Y half.
3. Deriving boundary triangles topologically from the tet faces (so
   they always match the per-side splits).
4. Writing a Gmsh v2.2 ASCII file directly (NOT a `.geo` source).

Resolution is uniform: 200 m everywhere in the production run.  Domain
is a 16 × 4 × 16 km box with the fault at Y=0 and the hypocenter at
(0, 0, −7500 m).

## What this gives up vs. an unstructured graded mesh

The TPV104 reference setup (`tpv104/mesh/tpv104_200m.geo`, the original
asymmetric mesh) is built in Gmsh with:

- Refinement near the fault (lc_fault = 200 m) and coarsening far from
  it (lc = 10 km), graded smoothly.
- Half-domain extent (60 km along strike, 60 km depth) — much larger
  than the structured mesh's 16 km × 16 km.
- Quality-optimized tets (Frontal/Delaunay algorithm with mesh-quality
  passes).
- A nucleation patch with extra refinement (lc_nucl = 200 m at the
  nucleation radius).

Properties our structured mesh lacks:

| Property | Graded TPV104 mesh | Structured mirror mesh |
|---|---|---|
| Refinement near fault | yes (200 m) | yes (200 m, but uniform everywhere) |
| Domain size | 60 × 60 × 60 km | 16 × 4 × 16 km |
| Far-field dispersion control | smooth grading | absorbing BC much closer |
| Mesh quality (aspect ratio) | optimized | identical for all hexes |
| Element count at 200m | ~215k tets | 768k tets (3.6×) |
| Mirror symmetry across y=0 | none (0/215k partners) | bit-exact |
| Curved or dipping faults | supported | NOT supported (Y=0 plane only) |

## Why we accepted this for now

1. The Y-mirror symmetry was the only way to get bit-exact zero σ_n
   leak in the dynamic phase (jobs 7677822, 7677831 confirmed σ_n flat
   on the symmetric mesh; Gmsh-generated unstructured meshes had
   0.46 MPa σ_yy pollution from non-mirror tet pairs).

2. Gmsh's tet mesher is **not deterministic** on mirror-symmetric input
   geometry — even feeding it a mirrored .geo doesn't guarantee a
   mirrored mesh.  Direct construction was the only reliable way to
   produce a mirror mesh.

3. For benchmark comparison against SeisSol, the smaller 16 × 4 × 16 km
   box is sufficient for the first 2 s of rupture (front doesn't reach
   absorbing BCs in the strike or depth directions; only the Y bound-
   ary is close, and absorbing BCs handle that).

## What this prevents

- **No curved/dipping fault simulations** with this mesh path.  The
  Cartesian construction assumes the fault is on a plane perpendicular
  to one axis (Y for TPV104).
- **No graded mesh for production-quality far-field accuracy**.  The
  uniform 200 m means computational cost scales as the volume of the
  full domain; can't refine selectively.
- **Larger element count** for a given fault region.  The 768k-tet
  16 × 4 × 16 km mesh is about 3.6× the 215k-tet 60 × 60 × 60 km mesh
  in element count even though the simulation domain is much smaller.
  Cost: ~3.6× more flux computations per macro-step.

## Path to removal

Three migration paths in increasing order of effort:

### Step 1 — Gmsh post-processing for symmetry enforcement

Mesh in Gmsh as usual (graded, large domain), then:

- Run a post-processing pass that, for each tet on the +Y half, finds
  its Y-mirror partner on the −Y half (or creates one if missing).
- Reorient tet vertex orderings so the partners' face splits agree
  topologically.

Risk: Gmsh-generated tets often have no exact mirror partner.  May
require deletion + replacement of mismatched tets.

### Step 2 — Constrained Gmsh meshing

Use Gmsh's `Mesh.AnisoMax` and `Field` infrastructure to force the
mesher to produce mirror-symmetric output:

- Mesh only the +Y half-domain.
- Apply Gmsh's `Symmetry { Surface ... } { 0,1,0,0 };` operation to
  reflect into the −Y half.
- Use `Coherence` to merge coincident vertices on the y=0 plane.

Risk: Gmsh's `Coherence` may produce non-conforming faces if the +Y
half's mesh on y=0 doesn't match the mirrored mesh exactly.  Probably
needs iteration on the meshing parameters.

### Step 3 — Drop the mirror requirement entirely

Move to an **algorithmic** fix (per-face fault representation, BP5-
style) that is robust to mesh non-mirror-symmetry.  Slip is per-face,
not per-side-DOF; the fault Riemann sees a single face-integrated
traction; partition asymmetry averages out.

This is the cleanest long-term fix but requires substantial rewrite of
`dynamic/fault_face_flux.cpp` (~500–1000 lines).  It also diverges
from SeisSol's per-side reference architecture, complicating cross-
code benchmark comparison.

## When to remove

- If we add curved or dipping fault problems → Step 1 or 2 (mesh-side).
- If we want to use the original 60 × 60 × 60 km TPV104 mesh for full
  benchmark comparison → Step 1 or 2.
- If we generalize the dynamic path for production SEAS → Step 3.

For the current TPV104 verification milestone, Step 1 is sufficient
once the σ_n / dip-direction comparison is established.
