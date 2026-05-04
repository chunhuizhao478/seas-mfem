# Investigation: Interior subdivision + needle collapse for fault-intersection mesh quality (2026-05-02)

## Visual evidence (user-supplied screenshots, 1.35 / 1.36 PM)

| Figure   | State                                       | Key observation                                                    |
|----------|---------------------------------------------|--------------------------------------------------------------------|
| Before   | Baseline + fix v2 step 5b (no break_wedges) | Sparse polyline, ~10 visible kink segments along the intersection. |
| After    | Fix v2 step 5 (canonical break_wedges, 30°) | **Dense zig-zag along polyline** (118 midpoints inserted), but triangles adjacent to the polyline are still elongated needles — long edges going to apexes far from the polyline. Net: γ_min improved 355× but **sliver count went UP** (78 → 129). |

The user's reading of the figures is correct: **inserting midpoints only along the polyline edges produces two needles where there was one**, because a 1→2 split halves one edge and leaves the other two unchanged.  The "narrow needle" geometry persists, and HXT's volume mesher extrudes those needles into wedge tets in 3D.

---

## Why polyline-edge-only splits are structurally insufficient

For a triangle `T = (a, b, c)` with edge `(a, b)` on the polyline and apex `c` far from the polyline:

- **1→2 split at midpoint `m`** of edge `(a,b)` → produces `(a, m, c)` and `(m, b, c)`.
- Both children inherit `c` as a vertex, so both have ONE short edge (`a→m` or `m→b`) and TWO long edges (`a→c`, `m→c` etc.).
- **Aspect ratio** of children ≈ same as parent.  Quality metric γ ≈ unchanged.

For a `(c)` apex 1 km from a 100-m polyline edge, the parent has aspect ~ 10:1.  After 1→2 split, both children have aspect ~ 10:1.  Loop subdivision (1→4 by inserting all three edge midpoints) drops aspect to ~ 5:1 because the new sub-triangles include `(m_ab, m_bc, m_ca)` in the centre — a similar triangle with all edges halved.  Still elongated, but better.

**The user's insight:** for a triangle `T` already in bad shape (γ_T < threshold) AND adjacent to the polyline, we need to insert a vertex *inside* `T` (a centroid or a Steiner point) and re-triangulate so that the new sub-triangles have controlled aspect.  This is fundamentally different from edge-midpoint splitting.

For a triangle that is a true needle (two near-coincident edges with sub-degree angle), no amount of Steiner insertion helps — the triangle should be **collapsed**: merge the two near-coincident vertices into one, eliminating the needle.

---

## The three operations needed (a complete remeshing toolkit)

### 1. Triangle interior subdivision (Steiner insertion)

For a bad-shape triangle `T = (a, b, c)` with γ_T < γ_thresh:

**Option A — centroid insertion (1→3):**
- Compute centroid `g = (a + b + c) / 3`.
- Replace `T` with `(a, b, g)`, `(b, c, g)`, `(c, a, g)`.
- All children share `g`; aspect ratio bounded by (max edge length) / (height from `g`).
- Cheap, but introduces a degree-3 interior vertex with poor regularity.

**Option B — full Loop subdivision (1→4):**
- Insert midpoints `m_ab, m_bc, m_ca` on all three edges.
- Replace `T` with `(a, m_ab, m_ca)`, `(m_ab, b, m_bc)`, `(m_ca, m_bc, c)`, `(m_ab, m_bc, m_ca)`.
- Children are similar triangles at half scale (corner three) plus one central similar triangle.
- T-junction problem: every neighbour triangle of `T` (which shares an edge with `T`) sees one of its edges suddenly split; the neighbour must do at least a 1→2 split to maintain manifoldness.  This is what `refine_fault_near_intersections.py` already implements (Loop + T-junction fix).

**Option C — constrained Delaunay re-triangulation:**
- Cluster the bad triangles into a connected patch.
- Compute a constrained Delaunay triangulation of the patch's boundary (with the polyline as a 1-D constraint).
- Replace the patch with the CDT result.
- Best quality, hardest to implement.  Triangle (Shewchuk) or `triangle` (Python wrapper) does this.

### 2. Needle / sliver triangle collapse

For a triangle that is too thin to subdivide usefully (e.g., the angle at one vertex is < 5°), we need **edge collapse**: merge two vertices into one.

For triangle `(a, b, c)` with very short edge `(a, b)` (length << others):
- **Collapse `(a, b)` to single vertex `m = (a + b) / 2`** (or one of the endpoints).
- Triangle `(a, b, c)` degenerates to a line `(m, c)` and is removed.
- Every neighbour triangle that referenced `a` or `b` now references `m`.

**Constraints:**
- If `a` or `b` is on the polyline, `m` MUST stay on the polyline (slide along the polyline rather than collapse to an endpoint).
- If `a` or `b` is on the free surface (z=0), `m` MUST stay on the free surface.
- Cross-fault conformity: collapsing in fault A must trigger the matching collapse in fault B if both share `(a, b)` as a polyline edge.
- Topology check: after collapse, verify no triangle becomes inverted (negative area) and no two triangles become identical.

This operation needs **half-edge data structures** for efficient neighbourhood traversal.  Pure NumPy gets clumsy at this scale.

### 3. Graded coarsening transition

The polyline-adjacent zone needs ~50-100 m triangles (for sliver-free meshing); the rest of the fault is fine at ~1000 m (the current `RES_F`).  A factor-10 jump in edge length would create artifact "pinch points" where the bulk mesher tries to grade quickly.

**Standard solutions:**
- **Geometric grading via mesh-size field**: `Field[2]` in `safs.geo` already does this for the bulk (Distance + Threshold from polyline), but the user is asking about the SURFACE triangulation — that's a 2-D problem before HXT runs.
- **MMG's `-hgrad` parameter**: caps the ratio of edge lengths between adjacent triangles (default `hgrad = 1.3` means ≤ 30 % size variation per step).  MMG re-triangulates while honouring this gradient.
- **Manual ring-by-ring relaxation**: insert N rings of intermediate-size triangles around the refined polyline zone before reaching the bulk size.  Tractable but tedious.

---

## Integration sketch (most surgical)

The operations have a natural ordering:

```
1. snap polyline vertices to canonical            ← already done (break_fault_wedges)
2. break_fault_wedges: split shallow-dihedral    ← already done
   polyline edges
3. for each fault, classify triangles by quality
   - keep      γ ≥ 0.3  (good)
   - subdivide γ ∈ [0.05, 0.3] AND edge-near
                polyline (1→4 Loop)
   - collapse  γ < 0.05 AND minimum angle < 5°
                (needle collapse)
4. graded smoothing: enforce hgrad ≤ 1.3 on
   triangle edge lengths
5. coarsen: remesh the rest of the fault to RES_F
6. handoff to generate_safs_mesh.py / HXT
```

Step 3-4 are exactly **`mmgs`** (the surface variant of MMG) when given a per-vertex target metric.  Step 5 is also mmgs at coarser target.  The current pipeline does step 1, 2, 6 only.

---

## Open-source tools (concrete pointers)

### Class A — surface remeshing with feature preservation (best fit)

| Tool | Repo | Why this fits |
|------|------|---------------|
| **MMG / mmgs** | https://github.com/MmgTools/mmg | Surface remesher with anisotropic metric, feature-edge protection (`-nr`), `-hgrad` gradient control, `-hausd` Hausdorff-distance bound (geometry preservation).  Reads / writes `.mesh` and `.medit` (also `.msh` via gmsh).  Splits / collapses / flips / smooths in one pass.  C++ with CMake build; conda package available (`conda install -c conda-forge mmg`).  Direct fit for steps 3 + 4 + 5. |
| **CGAL `PMP::isotropic_remeshing`** | https://github.com/CGAL/cgal | `Polygon_mesh_processing::isotropic_remeshing(mesh, target_edge_length, params...)` with `edge_is_constrained_map(ecm)` — protects the polyline edges while remeshing everything else.  Already in our CGAL 5.6.1 / 6.1 build (no new dep).  Operates on a single `Surface_mesh` per call — would need to run once per fault, with the polyline edges marked constrained on each. |
| **Geogram / vorpaline** | https://github.com/BrunoLevy/geogram | Bruno Lévy's restricted Voronoi diagram surface remeshing.  Higher quality output than `mmgs` on curved surfaces.  Feature-edge support via `vorpalite -profile=remesh`.  C++ standalone binary. |
| **OpenMesh** | https://www.graphics.rwth-aachen.de/software/openmesh/ | Half-edge data structure (no remesher; building block).  Use this if implementing custom collapse logic. |

### Class B — Steiner-augmented quality 3D tet meshing (alternative path)

If we cannot fix the surface, we can fix the volume.  These tools accept a constrained surface and produce well-shaped tets by inserting Steiner points:

| Tool | Repo | Why this fits |
|------|------|---------------|
| **TetGen with `-q1.0/15 -Y`** | https://github.com/libigl/tetgen | `-q1.0/15` requests min radius-edge ratio 1.0 AND minimum dihedral 15°.  TetGen inserts Steiner points to satisfy both.  `-Y` preserves input triangulation — we keep our conformal fault surfaces and TetGen only modifies the bulk.  Standalone binary; gmsh can hand off via `.poly` / `.smesh`. |
| **fTetWild** | https://github.com/wildmeshing/fTetWild | Robust tet meshing with envelope tolerance; tolerates non-manifold input.  Caveat: input geometry shifts by the envelope tolerance, breaking `triangle_to_fault.json` provenance.  Heavier integration. |
| **MMG3D** | https://github.com/MmgTools/mmg | The volume-mesh sister of mmgs.  Post-processes a tet mesh: edge swap, vertex collapse, point insertion to remove slivers.  Reads gmsh `.msh`.  This is the cleanest **post-HXT** sliver fixer.  `mmg3d_O3 -in safs.msh -out safs_mmg.msh -hgrad 1.3 -hausd 50.0` typically removes the slivers HXT produces. |

### Class C — bindings / wrappers (if Python integration matters)

| Tool | Repo | Notes |
|------|------|-------|
| **PyMesh** | https://github.com/PyMesh/PyMesh | Python wrappers around CGAL, MMG, Triangle, TetGen, IGL.  One-stop import.  Maintenance has slowed (last release 2020), but the existing bindings work. |
| **trimesh** | https://github.com/mikedh/trimesh | Pure-Python mesh I/O + simple operations (no remeshing).  Useful as glue. |
| **igl (libigl python bindings)** | https://github.com/libigl/libigl-python-bindings | Has a `igl.remesh_botsch` function (botsch_kobbelt isotropic remesher) that fits our use case.  pip-installable. |
| **py-mmg** | https://github.com/MmgTools/mmg/tree/develop/wrapper/python | Official MMG Python bindings (alpha as of 2024).  Reading mmg's mesh format from Python with `meshio` is the more stable path. |

---

## Recommendation in priority order

### Tier 1 — fastest path to working result (~1 day)

1. **Add a post-`break_fault_wedges` mmgs pass** that takes each fault's STL, marks the polyline edges as constrained, and remeshes with target edge ≈ 200 m near polyline → 1000 m bulk.  This implements steps 3 + 4 + 5 in one tool.

   Pipeline becomes:
   ```
   corefine_faults → break_fault_wedges → mmgs (per-fault) → generate_safs_mesh
   ```

   Sketch:
   ```python
   # mesh/mmgs_remesh_per_fault.py
   for short in include_faults:
       # Convert STL → .medit format, mark constrained edges from
       # polyline_keys (computed via collect_polyline_vertex_keys).
       run_mmgs(input_medit, output_medit,
                hmin=200, hmax=1000, hgrad=1.3, hausd=20.0,
                constrained_edges=polyline_edges)
       # Convert back to STL.
   ```

   Acceptance: sliver count drops below 30; γ_min above 0.01.  If mmgs respects polyline constraints (it does, with `-nr`), cross-fault conformity is preserved.

2. **Add a post-HXT mmg3d pass** as a safety net for residual slivers.  Even if mmgs gives perfect surface mesh, HXT still produces some slivers from the polyline geometry; mmg3d cleans them up.

   Pipeline becomes:
   ```
   ... → generate_safs_mesh → mmg3d → validate_msh
   ```

   `mmg3d_O3 -in safs.msh -out safs_mmg.msh -hgrad 1.3 -hausd 50.0 -hmin 100 -hmax 25000` typically gets γ_min above 0.05 across the whole mesh.

### Tier 2 — bespoke implementation if mmg can't be installed

Implement the three operations in pure Python on top of `OpenMesh` half-edge structures:

- `interior_subdivide_triangle(mesh, ti)` — 1→4 Loop with T-junction propagation (already mostly in `refine_fault_near_intersections.py`).
- `collapse_short_edge(mesh, ei)` — half-edge collapse with constraint preservation.
- `enforce_size_grading(mesh, target_field, hgrad=1.3)` — iterate until `max(adjacent_edge_ratio) ≤ hgrad`.

This is ~600 LoC + tests.  A focused 2-3 day effort.  Recommend ONLY if `conda install -c conda-forge mmg` is blocked.

### Tier 3 — research-grade (only if Tier 1 and 2 both fail)

Replace the surface generation entirely with **Geogram**'s RVD remesher, which is mathematically the cleanest approach for curved surfaces.  Or replace the volume meshing with **fTetWild** (accepts the messy fault soup directly, produces clean tets, with envelope tolerance documented and bounded).  These are 1-week+ integrations with implications for the SEAS provenance machinery.

---

## Cross-fault conformity considerations

The user's needle-collapse idea has a subtle conformity hazard.  When fault A collapses edge `(a, b)` to vertex `m`:

- If `(a, b)` is a polyline edge shared with fault B, fault B must collapse the same edge to the same `m` coord.  Otherwise the polyline diverges between faults and HXT rejects.
- If `(a, b)` is INTERIOR to fault A (not shared), fault A can collapse independently.

Decision rule: classify every candidate-collapse edge as polyline-shared or interior BEFORE collapsing.  For shared edges, collapse synchronously across both faults using the same canonical midpoint (just like `break_fault_wedges` already does).  For interior edges, collapse only in the owning fault.

MMG handles this correctly via its `-nr` (no edge swap on ridges) flag if the polyline edges are marked as ridges.  Custom Python implementation must replicate this gating logic.

---

## What this changes in the codebase

If we go Tier 1:

- **New file:** `mesh/mmgs_remesh_per_fault.py` (~150 LoC) — invokes mmgs binary per-fault with polyline edges marked constrained.
- **New file:** `mesh/mmg3d_post_pass.py` (~80 LoC) — invokes mmg3d binary on the gmsh output.
- **Modified:** `mesh/run_newset_step_by_step.sh` — two new optional stages:
  ```bash
  ENABLE_MMGS_REMESH=1   # before generate_safs_mesh
  ENABLE_MMG3D_POSTPASS=1 # after generate_safs_mesh
  ```
- **New tests:** `mesh/tests/test_mmgs_remesh.py`, `mesh/tests/test_mmg3d_post_pass.py` — verify constraint preservation + γ_min improvement.

The existing `refine_fault_near_intersections.py` and `break_fault_wedges.py` stay; they remain useful for non-mmg builds and as the core algorithmic backbone if Tier 2 is needed.

---

## Empirical validation criterion (re-using R-004 from prior review)

Same bit-identity guard from `test_refine_fault_near_intersections.py::assert_gamma_min_changed`:

| Run                                  | γ_min (current)         | γ_min target |
|--------------------------------------|-------------------------|--------------|
| Baseline (no break_wedges)           | 1.80e-10                | —            |
| Fix v2 (canonical break_wedges, 30°) | 6.40e-08                | (current)    |
| **Tier-1 with mmgs + mmg3d**         | (to measure)            | **> 1e-2**   |

A successful Tier-1 run improves γ_min by another 5-6 orders of magnitude AND drops the 129 sliver count below 30.

---

## Recommended next step

`conda install -c conda-forge mmg` in the `pythonenv`, then write `mmgs_remesh_per_fault.py` as a thin wrapper around the `mmgs_O3` binary.  The first end-to-end Tier-1 run should take about half a day to produce metrics.

If the user wants, I can:
1. Verify mmg installs cleanly in the current env.
2. Draft the `mmgs_remesh_per_fault.py` skeleton with constraint-edge marking.
3. Run end-to-end on Step 5 and report γ_min + sliver count.

Each step is gated; we stop and re-evaluate if mmg proves not to fit.
