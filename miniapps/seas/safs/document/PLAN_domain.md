# Implementation Plan: Domain Bounding Box, Mesh Sizing, and BP5-Compatible Tags

## Overview
This document fixes three open design choices in the SAFS mesh pipeline:
(1) the rule that determines the bulk-volume bounding box from the union of CFM fault
extents; (2) the Gmsh size-field configuration that puts ~1 km elements on the faults
and grades smoothly out to multi-km elements at the lateral/deep boundaries; (3) the
physical-tag scheme, which adopts BP5's exact convention so the existing SEAS-MFEM
solver can consume the resulting `.msh` without changes — in particular, all eight CFM
fault surfaces share a single physical surface tag `100`. Per-fault distinction is
preserved out-of-band via a side-car JSON file.

This plan supersedes Phase 3 (domain box) and Phase 5 (size field, physical groups) of
the master `PLAN.md`. The origin convention is independently fixed in `PLAN_origin.md`
(local UTM, origin `(500 000, 3 765 000, 0)` in Zone 11N).

## Constraints

- **BP5 tag scheme is non-negotiable for solver compatibility.** From
  `miniapps/seas/bp5/mesh/bp5_v2.geo`:
  ```
  Physical Surface("xm",   1)   x = -Lx        Dirichlet
  Physical Surface("xp",   2)   x = +Lx        Dirichlet
  Physical Surface("yp",   3)   y = +Ly        Dirichlet
  Physical Surface("ym",   4)   y = -Ly        Dirichlet
  Physical Surface("ztop", 5)   z = 0          Natural (free surface)
  Physical Surface("zbot", 6)   z = -Lz        Natural (deep)
  Physical Surface("fault",100) all fault subs Internal interface
  Physical Volume ("domain",10) bulk
  ```
  And from `miniapps/seas/CLAUDE.md`: `BCMode::FarField → attrs 1-4 = Dirichlet,
  attrs 5-6 = Natural; Tandem-style consumption: PS 1=Natural, PS 3=Fault, PS 5=Dirichlet`.
  The SAFS mesh follows the **MFEM/`bp5_v2`** scheme exactly (tags 1–6 for box faces, 10
  for volume, 100 for fault).

- **Single fault tag.** All eight CFM fault surfaces (Mojave SAF, San Bernardino SAF,
  Mill Creek, Mission Creek SBMT, Mission Creek COAV, Banning, Southern SAF + Banning,
  Pinto Mountain) collapse into Physical Surface 100. Per-fault information goes to
  `output/fault_provenance.json` (Phase 3 below). This matches how BP5 expects "the
  fault" to be a single attribute that the solver iterates over.

- **Free surface at z = 0.** Per `PLAN_origin.md` and `miniapps/seas/CLAUDE.md`. ztop
  (tag 5) is the only top face; the bulk volume does not extend above z = 0.

- **No symmetry exploitation.** BP5 has a symmetric mesh because its fault is planar.
  SAFS faults bend, branch, and cross; no symmetry. The domain box is rectangular but
  not centred on origin (because the CFM faults are not centred on UTM (500 000, 3 765 000)
  — Y is well-centred, X is biased eastward).

- **Element count budget for development tier.** At `res_f = 1000 m`, the mesh must fit
  in laptop memory (≤16 GB). Hard cap: 8 M tetrahedra (≈ 8 GB peak Gmsh memory). If
  generation overshoots this, the size field must be retuned, not the cap raised.

- **No hardcoded magic numbers** (per `feedback_no_hardcoded_numbers.md`). Buffer size,
  depth, `res_f`, `res_ff`, transition width are all CLI parameters with documented
  defaults and rationale.

## Phase 1: Domain Bounding-Box Rule

### Goal
After this phase, `generate_safs_mesh.py` computes a deterministic, justified bulk box
from the per-fault local-frame bounding boxes plus CLI-controlled buffers, and passes
those numbers to `safs.geo` via `-setnumber`.

### Files to Create
None (extends scripts already specified in `PLAN.md`).

### Files to Modify
- `mesh/generate_safs_mesh.py` — compute and pass the box extents.
- `mesh/safs.geo` — accept the extents as `DefineConstant`s (already in the master plan,
  reaffirmed here).

### Detailed Requirements

1. **Inputs.** `bbox.json` (per-fault local-frame bounding boxes, written by
   `ts_to_stl.py` per `PLAN_origin.md` Phase 3) and CLI flags. Per the local-frame
   transformation (`PLAN_origin.md`), the union bbox of the eight CFM faults is
   approximately:
   ```
   x_fault_min = -134 946 m,  x_fault_max = +122 329 m   (Δx ≈ 257 km)
   y_fault_min =  -72 722 m,  y_fault_max =  +73 998 m   (Δy ≈ 147 km)
   z_fault_min =  -18 014 m,  z_fault_max =        0 m   (Δz ≈  18 km, after z>0 clamp)
   ```

2. **Box extents — exact rule.**
   ```
   x_min_box = floor((x_fault_min - buf_x) / round_to) * round_to
   x_max_box = ceil ((x_fault_max + buf_x) / round_to) * round_to
   y_min_box = floor((y_fault_min - buf_y) / round_to) * round_to
   y_max_box = ceil ((y_fault_max + buf_y) / round_to) * round_to
   z_top_box = 0                                               # free surface
   z_bot_box = -depth_box
   ```
   with parameters (and defaults):
   ```
   buf_x     = 50 000 m   (≈ 2.8 × max fault depth, ≈ 0.4 × fault X-extent)
   buf_y     = 50 000 m   (≈ 2.8 × max fault depth, ≈ 0.7 × fault Y-extent)
   depth_box = 50 000 m   (≈ 2.8 × max fault depth)
   round_to  = 10 000 m   (snap box faces to 10-km grid for sane tag IDs in v2 viewers)
   ```

3. **Why these defaults — buffer and depth.**
   - **The standard SEAS guideline** (Tandem documentation, BP5 spec, Erickson+Day SCEC
     reports) is: "lateral buffer ≥ 1× to 2× the fault depth extent W; domain depth
     ≥ 2× to 3× W". For SAFS, W ≈ 18 km (SBMT-SAF reaches z = -17 km, ETRA-Pinto reaches
     z = -17.4 km). So `buf ≥ 36 km` and `depth ≥ 36 km`. Defaults 50 km / 50 km clear
     these bounds with comfortable margin.
   - **The buffer is *not* sized as a multiple of the fault X- or Y-extent.** The
     loading mechanism in BP5/SAFS is plate-rate Dirichlet on the lateral faces, which
     drives the shear stress on the fault uniformly regardless of buffer width — the
     buffer's job is only to keep the artificial Dirichlet boundary far enough that
     stress concentrations near fault tips don't reach it. A buffer ≈ W (not ≈ L) is
     the relevant scale.
   - **`round_to = 10 km`** snaps the box to a 10-km grid. With the defaults, the SAFS
     box becomes:
     ```
     X ∈ [-190 000, +180 000] m   (Lx = 370 km)
     Y ∈ [-130 000, +130 000] m   (Ly = 260 km)
     Z ∈ [ -50 000,        0] m   (Lz =  50 km)
     ```
     Rounded numbers make boundary-face bbox queries in `safs.geo` (Phase 3 below)
     robust against floating-point drift.

4. **Validation: every fault sits strictly inside the box.** After computing the box,
   the driver asserts:
   ```python
   for fault_name, b in per_fault_bbox.items():
       assert b["xmin"] > x_min_box + 0.1 * buf_x, fault_name
       assert b["xmax"] < x_max_box - 0.1 * buf_x, fault_name
       # same for y; z handled separately because faults touch z=0
       assert b["zmin"] > z_bot_box + 0.1 * depth_box, fault_name
   ```
   The `0.1 * buf` margin catches the case where a user passes too-small `buf_x` and
   leaves a fault flush against the lateral boundary (which would couple BC and fault
   physics nonsensically). Hard fail with a message naming the offending fault.

5. **CLI surface for `generate_safs_mesh.py`.**
   ```
   --buf-x METRES       default 50000   Lateral buffer along X
   --buf-y METRES       default 50000   Lateral buffer along Y
   --depth METRES       default 50000   Domain depth (Z extent below z=0)
   --round-to METRES    default 10000   Snap box faces to this grid
   ```
   The driver also writes `output/domain_box.json` recording the resolved
   `(x_min, x_max, y_min, y_max, z_top, z_bot, buf_x, buf_y, depth, round_to)` for
   reproducibility.

6. **`safs.geo` consumption.** The driver passes:
   ```
   gmsh -3 mesh/safs.geo \
        -setnumber x_min  <x_min_box> \
        -setnumber x_max  <x_max_box> \
        -setnumber y_min  <y_min_box> \
        -setnumber y_max  <y_max_box> \
        -setnumber z_top  0 \
        -setnumber z_bot  <z_bot_box> \
        ...
   Box(1) = { x_min, y_min, z_bot,
              x_max - x_min, y_max - y_min, z_top - z_bot };
   ```
   The `Box` is created in the local frame; `Volume{1}` is later tagged
   `Physical Volume("domain", 10)` (Phase 3).

### Interfaces
- Driver writes `output/domain_box.json`:
  ```json
  {
    "x_min_m": -190000.0,
    "x_max_m":  180000.0,
    "y_min_m": -130000.0,
    "y_max_m":  130000.0,
    "z_top_m":       0.0,
    "z_bot_m":  -50000.0,
    "buf_x_m":   50000.0,
    "buf_y_m":   50000.0,
    "depth_m":   50000.0,
    "round_to_m":10000.0,
    "Lx_m":     370000.0,
    "Ly_m":     260000.0,
    "Lz_m":      50000.0
  }
  ```

### Edge Cases to Handle
- A fault's `bbox.json` entry has `zmin < z_bot_box` (fault extends below the box):
  the validation in (4) catches it; user must increase `--depth`.
- A fault's `xmin < x_min_box` (fault extends outside the buffered box laterally):
  same — validation catches it. With the default 50 km buffer no current fault
  triggers it.
- User passes `--buf-x 0` (no buffer): allowed for debugging (e.g., to mesh just the
  fault complex without far-field), but emit a `WARNING: zero lateral buffer; lateral
  Dirichlet BC will couple directly to fault stress field` to stderr.
- `round_to` does not divide the buffer cleanly: that is fine — `floor`/`ceil` snap to
  the next grid line, which only enlarges the box.

### Acceptance Criteria
- [ ] With defaults, `output/domain_box.json` matches the box `(-190, +180) × (-130, +130) ×
      (-50, 0)` km (after the 10-km snap).
- [ ] All eight CFM faults sit strictly inside the box with at least 5 km clearance from
      every lateral face and from z_bot.
- [ ] `gmsh -3 mesh/safs.geo` (with the driver-supplied `-setnumber` values, mesh-skipped
      `Mesh.Algorithm3D = 1; Mesh.MeshSizeMax = 100000;` for the smoke test) produces a
      single bulk volume with Lx × Ly × Lz matching `domain_box.json` to ±1 m.
- [ ] Passing `--buf-x 1000` (insufficient) hard-fails the driver with the message
      "fault SAFS-SAFZ-MULT-Southern_SAF_and_Banning-CFM6 too close to lateral boundary"
      (or whichever fault is closest to the X-extreme).

### Dependencies
- Depends on: `PLAN_origin.md` Phases 1–3 (origin module + `transform.json` + `bbox.json`).
- Required by: Phase 2 of this document.

---

## Phase 2: Mesh Sizing Field — Fault Refinement with Smooth Coarsening

### Goal
After this phase the SAFS mesh has `≈ res_f` (1000 m default) tetrahedra in a thin shell
around every fault surface, smoothly grades to `res_ff` (20 km default) at the lateral
and deep boundaries, and the total tet count is below the budget cap (~8 M at the
development tier).

### Files to Create
None.

### Files to Modify
- `mesh/safs.geo` — add the `Distance` + `Threshold` size-field block.
- `mesh/generate_safs_mesh.py` — accept `--res-f`, `--res-ff`, `--ramp-dist` and pass
  via `-setnumber`.

### Detailed Requirements

1. **Field topology.**
   ```
   // Sample every fault triangle for distance computation
   Field[1] = Distance;
   Field[1].SurfacesList = { fault_surfs() };  // ALL post-fragment fault patches
   Field[1].Sampling     = 100;                // points per surface in the distance octree

   // Linear ramp: res_f at distance 0 → res_ff at distance ramp_dist
   Field[2] = Threshold;
   Field[2].InField  = 1;
   Field[2].SizeMin  = res_f;
   Field[2].SizeMax  = res_ff;
   Field[2].DistMin  = 0;
   Field[2].DistMax  = ramp_dist;
   Field[2].Sigmoid  = 0;        // linear, NOT sigmoid (sigmoid creates locally fine
                                 //  bands inside the transition that bloat element count)
   Field[2].StopAtDistMax = 0;

   Background Field = 2;

   Mesh.MeshSizeExtendFromBoundary = 0;   // do not extend point sizes from box corners
   Mesh.MeshSizeFromPoints         = 0;   // ignore Point() Lc values
   Mesh.MeshSizeFromCurvature      = 0;   // no curvature-based refinement
   Mesh.MeshSizeMin                = res_f;
   Mesh.MeshSizeMax                = res_ff;
   ```

2. **Why Distance + Threshold (not BoundaryLayer, not MathEval).**
   - `BoundaryLayer` makes prismatic layers along the surface — wrong physics. We want
     isotropic tets that simply get smaller near the fault.
   - `MathEval` over a custom `f(x,y,z)` could encode region-specific behavior (e.g.,
     extra refinement at the San Gorgonio knot), but it bypasses Gmsh's distance octree
     and is far slower at SAFS scale. Distance + Threshold uses Gmsh's spatial index and
     is the canonical pattern for fault-conformal refinement (BP5 uses the same pattern
     when sizes vary).
   - `Sampling = 100` per surface gives sub-100-m distance accuracy across the eight
     CFM faults at the cost of ≈800 sample points total — negligible.

3. **Defaults and ramp width.**
   ```
   res_f      = 1000 m       (target on-fault edge length, development tier)
   res_ff     = 20000 m      (far-field edge length at the box boundary)
   ramp_dist  = 30000 m      (linear transition width)
   ```
   Justification of `ramp_dist`:
   - The grading factor is `res_ff/res_f = 20`. With a linear ramp, an element of
     edge length `h` sits at distance `(h - res_f)/(res_ff - res_f) * ramp_dist` from
     the fault. Adjacent tets along the gradient direction differ in edge length by
     `h_step ≈ res_f * ramp_dist / (h * 1.0)`. To keep neighbour-edge ratios below
     ~1.3 (Gmsh's "reasonable grading" rule of thumb), we need
     `ramp_dist ≥ 1.5 × (res_ff − res_f)`. With `res_ff = 20 km, res_f = 1 km`:
     `ramp_dist ≥ 28.5 km`. Default 30 km satisfies this.
   - `ramp_dist` cannot exceed the smaller buffer: if the lateral buffer is 50 km and
     `ramp_dist = 30 km`, the ramp completes 20 km before reaching the lateral
     boundary, so the boundary tets are fully at `res_ff`. If a user reduces
     `--buf-x` below `ramp_dist`, the driver emits a warning.

4. **Element-count estimate.** With the SAFS box `370 × 260 × 50 km` and the size field
   above, the tet count is dominated by:
   - **Near-fault shell** (thickness ≈ 5 km, where elements are < 5 km): area
     ≈ Σ (fault areas) ≈ 8 250 km² (sum of CFM `Area(km²)` from the inventory table).
     Volume ≈ 8 250 × 5 = 41 250 km³. At avg edge length ≈ 2.5 km in this shell, tet
     volume ≈ (2.5)³ / 6 ≈ 2.6 km³, so ≈ 16 000 tets.
   - **On-fault dense band** (thickness ≈ 1 km, edge length 1 km): volume ≈ 8 250 km³,
     tet volume ≈ 0.17 km³, so ≈ 49 000 tets.
   - **Transition zone** (thickness 30 km - 5 km = 25 km, mean edge ≈ 10 km):
     volume ≈ Σ areas × 25 km × 2 sides = 412 500 km³. Tet volume ≈ 167 km³ ⇒ 2 500 tets.
   - **Far-field** (rest of box, edge ≈ 20 km): volume ≈ 4.5 M km³ - above ≈ 4.0 M km³.
     Tet volume ≈ 1 333 km³ ⇒ 3 000 tets.

   **Total estimate: ≈ 70 k–100 k tets** at `res_f = 1000 m`. Comfortably under the 8 M
   cap. (The 8 M cap is reserved for `res_f = 250 m` production tier, not 1 km.)

5. **CLI surface.**
   ```
   --res-f     METRES   default 1000    Target edge length on every fault
   --res-ff    METRES   default 20000   Edge length far from any fault
   --ramp-dist METRES   default 30000   Linear ramp width
   ```
   The driver writes `output/sizing.json` with the resolved values and the *measured*
   tet count after meshing (read back from the `.msh` file via `meshio`).

6. **Resolution tiers.**
   | Tier | `--res-f` | `--res-ff` | `--ramp-dist` | Expected tets | Use |
   |------|-----------|------------|---------------|---------------|-----|
   | dev  | 1000      | 20000      | 30000         |   100 k       | smoke / iteration |
   | sci  |  500      | 15000      | 25000         |   600 k       | calibration runs |
   | prod |  250      | 10000      | 20000         |   5–8 M       | publication runs |

   The driver enforces `--res-f >= ramp_dist / 50` (so the ramp does not become
   sub-element-scale) and `--res-ff <= ramp_dist * 1.5` (so far-field elements still
   resolve the ramp). Violations emit a hard error with the offending values.

### Interfaces
- `safs.geo` reads `res_f`, `res_ff`, `ramp_dist` via `DefineConstant`. Defaults match
  the dev tier so interactive `gmsh safs.geo` produces a sane mesh without driver flags.

### Edge Cases to Handle
- `Field[1]` with `SurfacesList = {}` (the fault list resolved to empty due to a
  fragmentation failure in `PLAN.md` Phase 4): Gmsh produces a uniform mesh at
  `res_ff`. Detect this in `validate_msh.py` (mean edge ≈ res_ff with low variance) and
  fail.
- Two faults with overlapping bounding-box-defined sub-regions (e.g., Mission Creek
  and Mill Creek separated by a few hundred metres in the SBMT region): the Distance
  field automatically takes the *minimum* distance to any surface in `SurfacesList`,
  so the dense band wraps the union — correct behaviour.
- `--res-f > --res-ff` (user typo): driver hard-errors with "size-field inversion".

### Acceptance Criteria
- [ ] At default settings (`res_f=1000, res_ff=20000, ramp_dist=30000`), the produced
      mesh has `100 000 < n_tets < 1 500 000` (factor-of-15 envelope around the
      back-of-envelope 100 k estimate).
- [ ] Mean edge length within 1 km of any fault triangle is `1000 ± 300 m`. Verified
      by `validate_msh.py` sampling of fault-side tets.
- [ ] Mean edge length on the lateral box faces is `> 12 000 m`. Verified by sampling
      tets touching tag 1, 2, 3, or 4.
- [ ] Tet quality `min Gamma ≥ 0.10`, `mean Gamma ≥ 0.55`. From `gmsh -info` after meshing.
- [ ] `output/sizing.json` exists and contains the measured tet count.

### Dependencies
- Depends on: Phase 1 (box dimensions); `PLAN.md` Phase 4 (fault surfaces stitched).
- Required by: Phase 3.

---

## Phase 3: BP5-Compatible Physical Tags + Per-Fault Provenance Side-Car

### Goal
After this phase the SAFS `.msh` has exactly the tag scheme of `bp5_v2.geo` (so the
existing solver consumes it without modification), and a separate JSON file records
which CFM fault each fault triangle came from for downstream per-fault parameter
mapping.

### Files to Create
- `mesh/write_fault_provenance.py` — post-process the `.msh` to produce
  `output/fault_provenance.json`.

### Files to Modify
- `mesh/safs.geo` — add the BP5-conforming `Physical` block (replacing the per-fault
  tag scheme proposed in `PLAN.md` Phase 5).

### Detailed Requirements

1. **Physical-tag scheme — exact copy of BP5 plus a single-fault collapse.**
   ```
   // --- Boundary surfaces ---
   xm()   = Surface In BoundingBox{ x_min - eps, y_min - eps, z_bot - eps,
                                    x_min + eps, y_max + eps, z_top + eps };
   xp()   = Surface In BoundingBox{ x_max - eps, y_min - eps, z_bot - eps,
                                    x_max + eps, y_max + eps, z_top + eps };
   yp()   = Surface In BoundingBox{ x_min - eps, y_max - eps, z_bot - eps,
                                    x_max + eps, y_max + eps, z_top + eps };
   ym()   = Surface In BoundingBox{ x_min - eps, y_min - eps, z_bot - eps,
                                    x_max + eps, y_min + eps, z_top + eps };
   ztop() = Surface In BoundingBox{ x_min - eps, y_min - eps, z_top - eps,
                                    x_max + eps, y_max + eps, z_top + eps };
   zbot() = Surface In BoundingBox{ x_min - eps, y_min - eps, z_bot - eps,
                                    x_max + eps, y_max + eps, z_bot + eps };
   eps    = 1.0;   // metres; round_to=10000 keeps box faces well-separated

   Physical Surface("xm",   1) = { xm() };       // x = x_min  (Dirichlet, far-field)
   Physical Surface("xp",   2) = { xp() };       // x = x_max  (Dirichlet, far-field)
   Physical Surface("yp",   3) = { yp() };       // y = y_max  (Dirichlet, far-field)
   Physical Surface("ym",   4) = { ym() };       // y = y_min  (Dirichlet, far-field)
   Physical Surface("ztop", 5) = { ztop() };     // z = 0      (Natural, free surface)
   Physical Surface("zbot", 6) = { zbot() };     // z = z_bot  (Natural, deep boundary)

   // --- Fault surface (UNIFIED, BP5-compatible) ---
   // All eight CFM fault sub-surfaces collapse into tag 100. Tag 100 is the SAME tag
   // bp5_v2.geo uses for its fault — the SEAS-MFEM solver iterates this attribute as
   // "the fault" without caring how many physical surfaces produced it.
   Physical Surface("fault", 100) = { all_safs_fault_surfs[] };

   // --- Volume ---
   Physical Volume("domain", 10) = { Volume{:} };
   ```
   `eps = 1.0 m` works because `round_to = 10 000 m` from Phase 1 puts the box faces on
   round 10-km grid lines, far from any fault vertex.

2. **Identifying `all_safs_fault_surfs[]`.** After the `BooleanFragments` step in
   `PLAN.md` Phase 4, each pre-fragment fault STL produces a list of post-fragment
   surface tags. Concatenate those eight lists into one flat list:
   ```
   all_safs_fault_surfs[] = {
       safs_mjvs_saf_surfs[]{},
       safs_pmfz_pinto_surfs[]{},
       safs_sbmt_millcreek_surfs[]{},
       safs_sbmt_missioncreek_surfs[]{},
       safs_coav_missioncreek_surfs[]{},
       safs_mult_banning_surfs[]{},
       safs_mult_ssaf_banning_surfs[]{},
       safs_sbmt_saf_surfs[]{}
   };
   ```
   The shared 1-D edge curves at fault-fault intersections (Mission Creek meets SAF,
   Banning meets SSAF, etc.) are *not* added explicitly — they are already shared
   between the fault patches by virtue of the fragmentation, and the surface mesh on
   each side of the shared edge inherits the edge's mesh nodes. The unified tag 100
   captures every triangle of every fault.

3. **Per-fault provenance side-car — file format.**
   ```json
   {
     "schema_version": 1,
     "mesh_file": "output/safs_1000m.msh",
     "fault_tag": 100,
     "transform_json": "mesh/transform.json",
     "faults": {
       "safs_mjvs_saf": {
         "cfm_id": "SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6",
         "section_label": "Mojave",
         "avg_strike_deg": 116.0,
         "avg_dip_deg": 90.0,
         "stl_file": "mesh/stl/safs_mjvs_saf.stl",
         "n_triangles_in_msh": 14237,
         "triangle_indices_in_msh": [12001, 12002, ..., 26237]
       },
       "safs_sbmt_saf": { ... },
       ...
     }
   }
   ```
   `triangle_indices_in_msh` are 0-based indices into the `.msh` triangle element list
   filtered by tag 100. Each triangle appears in exactly one fault entry; the union
   of all `triangle_indices_in_msh` lists equals the full set of tag-100 triangles.

4. **`write_fault_provenance.py` algorithm.**
   ```
   load mesh via meshio
   load transform.json + per-fault STL files (in local frame, from PLAN_origin Phase 2)

   for each tag-100 triangle T in mesh:
       compute centroid c_T (local frame)
       for each pre-fragment STL S in {8 STLs}:
           compute squared distance from c_T to nearest triangle of S using
              KDTree of S's triangle centroids
       assign T to argmin S
       if min distance > 2 * res_f:
           mark T as "ambiguous", attach to nearest fault but log a warning

   build the JSON above and write to output/fault_provenance.json
   ```
   Implementation: use `scipy.spatial.cKDTree` over the union of pre-fragment STL
   triangle centroids, with a per-tree label = which STL the centroid came from.
   For each mesh triangle centroid, query the single nearest centroid; the label of
   that nearest centroid is the assigned fault. This is O((N_mesh + N_stl) log N_stl)
   and runs in seconds for the SAFS-scale data.

5. **Ambiguity tolerance.** "Ambiguous" triangles (assigned with > 2*res_f distance)
   should be < 0.5% of all tag-100 triangles. They occur at fault-fault intersection
   curves where the post-fragment triangulation has been re-split. Acceptable as long
   as their solver-side parameters (a, b, V_init, etc.) are continuous across faults
   (which they are for any reasonable physical setup).

6. **CLI for `write_fault_provenance.py`.**
   ```
   python write_fault_provenance.py \
       --msh output/safs_1000m.msh \
       --stl-dir mesh/stl \
       --out output/fault_provenance.json
   ```

### Interfaces
- `output/fault_provenance.json` — the contract between mesh generation and any
  future solver-side per-fault parameter mapping. Schema-versioned (`schema_version: 1`)
  so it can evolve.
- `safs.geo` Physical block — the contract between mesh generation and the existing
  SEAS-MFEM solver (BP5 attributes 1–6, 10, 100).

### Edge Cases to Handle
- A tag-100 triangle's centroid is exactly equidistant from two CFM STL centroids
  (e.g., on the San Gorgonio shared edge): KDTree breaks ties deterministically by
  index order; log this as an ambiguous case but accept the assignment.
- `--stl-dir` STLs are in a different frame than the mesh (e.g., user re-ran
  `ts_to_stl.py` with overridden origin between mesh generation and provenance):
  cross-check `transform.json` SHA hashes; refuse to proceed on mismatch.
- A pre-fragment STL has zero triangles in the final mesh (a fault was entirely
  consumed/cancelled by fragmentation): write its entry with
  `"n_triangles_in_msh": 0`, no `triangle_indices_in_msh` field, and an
  `"omitted_reason": "fragmented away"` field. Validation script flags this as a
  data-quality issue.

### Acceptance Criteria
- [ ] `output/safs_1000m.msh` contains physical groups exactly: 1, 2, 3, 4, 5, 6
      (each present, non-empty), 10 (single volume, present), 100 (single fault group,
      non-empty), and **no other physical tags**.
- [ ] `validate_msh.py` reports the same tag inventory as `bp5_v2.geo` produces (only
      tag 100 is present where 100–106 used to be in BP5; this is intentional and
      covered by the validator's "SAFS profile").
- [ ] `write_fault_provenance.py` runs in < 60 s on the dev-tier mesh and produces
      a JSON with eight fault entries summing to `n_triangles_in_msh = N_total_tag100`
      (no triangle is missed or double-counted).
- [ ] Number of "ambiguous" triangles is < 0.5% of tag-100 total at the dev tier.
- [ ] The existing SEAS-MFEM solver's mesh reader (whatever entry point loads BP5
      meshes) loads the SAFS mesh without errors. Note: this is verified manually by
      the user, since auto-running the solver is out of scope for this plan.

### Dependencies
- Depends on: Phases 1, 2; `PLAN.md` Phases 3, 4, 5.
- Required by: nothing (terminal of this sub-plan).

---

## Testing Strategy

| Test | Where | What it proves |
| ---- | ----- | -------------- |
| Domain-box rule unit test | `pytest` for `generate_safs_mesh.py:_compute_box()` | Box snaps to grid, contains all faults with margin |
| Insufficient-buffer rejection | `pytest` | `--buf-x 1000` triggers the explicit hard-fail |
| Tet count within envelope | integration | dev-tier mesh has 100 k–1.5 M tets |
| Near-fault edge length | `validate_msh.py` | tets within 1 km of fault have edge ≈ res_f |
| Far-field edge length | `validate_msh.py` | tets touching boundary have edge ≈ res_ff |
| BP5 tag inventory | `validate_msh.py` | exactly tags {1,2,3,4,5,6,10,100} present |
| Provenance partition | `write_fault_provenance.py` self-check | sum of per-fault tri counts = total tag-100 tri count |
| Solver-load smoke | manual | existing BP5 mesh-loader code reads SAFS mesh |

## Risk Assessment

1. **The SBMT-SAF surface dips at 51° and reaches the deepest point (-17 km).** Its
   western tip is also the most negative X in the union (the Mojave SAF reaches
   -134.9 km). With `--buf-x 50 km` and `--depth 50 km`, both are 35 km from the
   nearest box face — comfortable. If the user passes `--depth 25 km`, the SBMT-SAF
   pierces the box bottom; the validation in Phase 1 step 4 catches this before
   meshing.

2. **The single-fault tag complicates per-fault friction parameters.** Solvers that
   want different `a`, `b`, `Dc` per CFM fault (e.g., to reflect each fault's
   slip-rate regime) cannot read the assignment from the `.msh` alone — they must
   load `fault_provenance.json` and dispatch by triangle index. This is acceptable
   because (a) BP5 itself uses per-DOF spatially varying parameters loaded through a
   parallel side-channel (`config/bp5_params.hpp`), and (b) the SAFS side-car follows
   the same pattern. **However**, integrating `fault_provenance.json` into the solver
   is out of scope here and a known follow-up.

3. **Distance-field cost grows with fault triangle count.** At `res_f = 250 m` (prod
   tier) the eight CFM faults produce ~ 200 k surface triangles total; Gmsh's
   distance octree handles this in a few seconds. If a future SAFS variant adds many
   more faults (the SCEC CFM has ~150 fault surfaces in total), the octree build
   becomes the bottleneck — increase `Field[1].Sampling` correspondingly or split
   into multiple `Distance` fields with `Min` aggregation. Out of scope for the
   current eight-fault model.

4. **`Surface In BoundingBox` for box faces relies on the `round_to` snap.** If a
   future change reduces `round_to` below ~100 m, fault triangles near the lateral
   boundary may sneak inside the bounding-box query for `xm`/`xp`/`yp`/`ym` and get
   tagged as Dirichlet boundaries instead of remaining tag 100. The `eps = 1.0 m`
   query tolerance is calibrated for `round_to = 10 km`. Document this in the
   `safs.geo` header and assert `round_to >= 100 * eps` in `generate_safs_mesh.py`.

5. **BP5 v1 vs v2 BC mode.** `bp5.geo` (v1) uses ztop=Natural; `bp5_v2.geo` (v2) uses
   ztop=Dirichlet for the H10 fix. SAFS adopts v1's mapping (ztop = free surface,
   Natural). This must be set in the solver via `BCMode::FarField` (per
   `miniapps/seas/CLAUDE.md`). The SAFS plan does **not** require any solver-side
   changes beyond pointing it at the new mesh and selecting `BCMode::FarField`.
