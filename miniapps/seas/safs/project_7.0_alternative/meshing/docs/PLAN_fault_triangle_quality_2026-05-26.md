# Implementation Plan: Lift SAFS fault-triangle quality (min-angle ≥ 25–30° / tri_q ≥ 0.5–0.6)

> **Status: DRAFT 2026-05-26.** Awaiting review/approval.
> Target mesh:
> `meshing/results/msh/safs_fault_box_nwcut_500m_lcfar3000_z0embed.msh`.
> **Hard ground rule from the user: do NOT modify the existing mesh or the
> existing code.** Every artefact this plan produces is a *new* file
> (new remesh tool, new STL, new `.msh`, new test). The existing mesher
> `run_z0cut_meshing.py` and the existing `.msh`/`.stl` are *run* and
> *read*, never edited.

## Overview

The named mesh's **bulk** already clears the project's Q1/Q2 gates
(`tet_emin = 112 m ≥ 100 m`; Joe-Liu `η_min = 0.184 > 0.1`). The problem
is purely on the **fault surface** (Physical Surface 101): of 42,358 fault
triangles, the median quality is excellent (`tri_q = 0.9995`,
`min_angle ≈ 59°`), but a small tail is bad:

| Metric | Baseline (target mesh) | Goal |
|---|---:|---:|
| `tri_qmin` | **0.2446** | ≥ 0.5 (stretch 0.6) |
| worst `min_angle` | **9.10°** | ≥ 25° (stretch 30°) |
| tris with `tri_q < 0.5` | **45** | 0 |
| tris with `tri_q < 0.6` | **136** | ≤ a few irreducible corner tris |

This plan improves those triangles by **remeshing the input STL surface**
(holding the z = 0 trace and the NW-cut/perimeter borders fixed), then
re-running the *existing* `run_z0cut_meshing.py` against the improved STL
to emit a **new** `.msh`. It is verified with the *existing*
`check_mesh_quality.py` + `locate_fault_slivers.py` tools plus a small new
fault-triangle gate.

## Why the fix must live in the STL, not the `.msh` (root-cause diagnosis)

Measured on the target mesh (2026-05-26, `pythonenv`):

1. **Fault triangle count == STL triangle count == 42,358.** gmsh
   `merge`s the cut STL as a *discrete* surface and embeds it
   (`gmsh.model.mesh.embed(2, [fault_surf], 3, vol)`); it does **not**
   re-triangulate it. Therefore **fault-triangle quality ≡ input-STL
   triangle quality** (verified: STL `meshio` tri count equals the
   Physical-Surface-101 tri count exactly).
2. **The bad triangles sit on the z = 0 trace.** Of the 45 tris with
   `tri_q < 0.5`, **37 touch z = 0**; of the 136 with `tri_q < 0.6`,
   **124 touch z = 0**. These are the thin triangles produced by
   `ts_to_stl.py`'s *per-triangle linear clip* at z = 0 (each fault
   triangle straddling the surface is cut, leaving a needle along the cut
   line). The 3 worst (`tri_q < 0.3`, at z ≈ −900 … −4300 m) are on the
   NW-cut edge / source-`.ts` artefacts.
3. **You cannot retriangulate the fault inside the `.msh`.** Every
   Physical-Surface-101 triangle is a *shared face of two bulk tets*
   (the fault is embedded). Re-triangulating the surface alone would
   break volume conformity. The only conformity-safe options are (a) feed
   gmsh a better surface (this plan's primary path) or (b) move nodes
   without changing connectivity (TMOP surface-fitting — a limited
   *polish*, see Alternatives).

Conclusion: improve the STL surface mesh — especially the z = 0 clip
slivers — while preserving the fault geometry and the boundary loops, then
re-mesh.

### Refinement does NOT raise the minimum angle (why we remesh, not refine)

h-refinement (subdividing) and shape (angle) are **orthogonal** mesh-quality
axes. Splitting a triangle 1→4 by edge midpoints yields sub-triangles
*similar* to the parent — a 9° sliver stays a 9° sliver, you just get more
of them. Raising the minimum angle requires either node *movement*
(smoothing / TMOP) or *connectivity change* (edge swap / collapse /
remesh). This is why every path in this plan **remeshes** the surface
rather than refining it.

Classified worst-angle triangles on the target STL (min-angle < 25°,
measured 2026-05-26):

| Class | Count | Removable by remesh? | Removable by refinement? |
|---|---:|---|---|
| **Clip needles** (≤ 1 boundary edge; 171 of them touch z = 0) | **185** | ✅ collapse + retriangulate | ❌ no |
| **Geometry-forced corners** (≥ 2 boundary edges) | **2** | ❌ no | ❌ no (gets strictly worse) |

The 185 needles are `ts_to_stl.py` clip artefacts and the remesh clears
them. The **2** geometry-forced tris are both stacked at the **NW tip of
the fault footprint**, ≈ `(365012, 3838896, z ≈ −3500 … −4300 m)`, where
the cut footprint tapers to a ~9° wedge. At a sharp domain corner the apex
angle is a hard floor that *no* triangulation can beat, and refinement only
makes it worse. Options for those 2 (see Phase 2 edge case): trim/blunt the
footprint tip (geometry change, needs approval), accept a relaxed local
bar, or leave them — they sit ~4 km deep at the far NW end, far from the
seismogenic/nucleation zone, so a sliver there is likely physically
harmless (unlike the near-surface trace needles tied to the 2026-05-20
blow-up).

## What tools we can use (survey)

All paths below were verified present on this machine on 2026-05-26.

| Tool | Where (env) | Role here | Pros | Cons / risk |
|---|---|---|---|---|
| **MMG `mmgs_O3` 5.8.0** | `pythonenv` (`.../envs/pythonenv/bin/mmgs_O3`) | **PRIMARY** surface remesher | purpose-built isotropic surface remesh; `-hausd` bounds geometric drift; open-surface free-boundary edges auto-preserved as required; uniform `-hsiz`; fast; no build | Medit `.mesh` I/O (need STL↔.mesh convert); boundary nodes must be re-snapped to exact z = 0 after; needs a clean manifold input |
| **gmsh 4.15.0** | `pythonenv` (`import gmsh`) | FALLBACK A: reclassify + remesh discrete surface | same mesher already trusted by the pipeline; reuses size-field idioms; emits v2.2 directly | `classifySurfaces` can over/under-segment a curved fault into multiple patches; trickier to hold one clean z = 0 boundary curve |
| **pymeshlab** `meshing_isotropic_explicit_remeshing` | `pythonenv` | FALLBACK B: VCG remesh | one-call Python API; boundary-preserving flag | coarser quality control than mmgs/CGAL; reprojects onto self only; `OMP` init caveat (set `KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1`) — see `EXPLORE_pymeshlab_intersections.md` |
| **CGAL 6.1.1** `isotropic_remeshing` / `surface_Delaunay_remeshing` | `cgal-61` (C++) | FALLBACK C: highest-quality remesh | best quality; `edge_is_constrained_map`/`polyline_constraints` hold the z = 0 trace exactly; Delaunay-refinement quality bound | needs an ~80-line C++ build in `cgal-61`; STL float I/O can drift coords (use OFF/PLY) — see `EXPLORE_cgal_corefine.md` |
| **MFEM TMOP** (`mesh-optimizer`/`pmesh-optimizer`) | `mfem-dev` | ALTERNATIVE: in-place node polish of the existing `.msh` | operates on the named `.msh` directly; keeps conformity (moves shared fault nodes consistently in fault tris + tets); surface-fitting slides trace nodes tangentially | **no topology change** — cannot collapse a pinned-node sliver; only a polish; boundary-attribute gotcha (101–104 ∉ {1,2,3,4}) — see `EXPLORE_tmop_mesh_optimization.md` |

**Recommended primary: MMG `mmgs_O3`** — it is installed, purpose-built
for exactly this (isotropic surface remesh of an open triangulated sheet
with geometric-error control and automatic boundary preservation), and
needs no compilation. The plan tries it first (Phases 1–2) and falls back
(Phase 3) only if it misses the targets.

## Constraints

These are inherited hard invariants of the SAFS meshing pipeline
(`README.md`, `PLAN_mesh_quality.md`, project `CLAUDE.md`). Any produced
artefact must respect all of them.

- **Do not modify existing code or meshes.** New files only. The existing
  `run_z0cut_meshing.py`, `check_mesh_quality.py`, `locate_fault_slivers.py`,
  the target `.msh`, and the source `.stl` are read/run, never edited.
- **Free-surface invariant**: the produced `.msh` must have
  `mesh_zmax == 0.0` exactly (verified by the README one-liner). This
  requires the remeshed STL's `max(z) == 0` exactly — `run_z0cut_meshing.py`
  *aborts* if `abs(fault_zmax) > 1e-6` (its line ~233 guard).
- **z = 0 trace must survive as a clean boundary polyline**: the remeshed
  STL must still present its z = 0 edges as free-boundary edges, so
  `run_z0cut_meshing.py:extract_z0_trace` finds them and the seam weld
  (`removeDuplicateNodes`) succeeds. A remesh that closes/rounds the top
  edge off z = 0 breaks the embed (empty-volume Option-A failure mode —
  see `DEBUG_z0cut_option_a_failure.md`).
- **Single connected fault surface**: `run_z0cut_meshing.py` requires
  *exactly one* merged non-box surface (its line ~269 guard). The remesh
  must not split the fault into multiple shells.
- **Fault stays embedded**: every Physical-Surface-101 triangle must be a
  face of exactly two bulk tets in the final mesh.
- **Geometry fidelity**: the remeshed fault must stay within a small
  Hausdorff distance of the original Fuis-ALT6 surface (default budget
  ≤ 25 m; ≪ the ~500 m fault edge and ≪ seismogenic length). The NW-cut
  footprint (the cut perimeter) must be preserved.
- **Gmsh v2.2 output only** (`-format msh22`); the MFEM reader is v2.2-only.
  `run_z0cut_meshing.py` already writes v2.2 — do not change it.
- **Physical tags preserved**: rock=1, fault=101, top=102, bottom=103,
  sides=104 (handled by `run_z0cut_meshing.py`, unchanged).
- **No hardcoded magic numbers**: target edge size, Hausdorff budget, and
  z-snap tolerance are CLI parameters with documented defaults derived
  from the input (per project `CLAUDE.md` / `feedback_no_hardcoded_numbers`).

## Phase 1: New STL surface-remesh tool (`remesh_fault_stl.py`)

### Goal
After this phase a new, standalone CLI tool takes the existing cut fault
STL and writes a **new** STL whose triangles are near-equilateral
(min-angle ≥ target everywhere except irreducible boundary corners), whose
`max(z) == 0` exactly, and whose z = 0 trace and NW-cut/perimeter borders
are geometrically preserved. The existing pipeline is untouched.

### Files to Create
- `meshing/code/remesh_fault_stl.py` — new standalone CLI (the remesh
  driver; engine = mmgs by default).

### Files to Modify
- None. (Existing tools are invoked as subprocesses / imports only.)

### Input / output
- **Input STL** (read-only):
  `meshing/results/stl_nwcut/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m_clean_clip_nwcut.stl`
  (42,358 tris; `z ∈ [−16607.6, 0.0]`).
- **Output STL** (new):
  `meshing/results/stl_nwcut/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m_clean_clip_nwcut_triq.stl`
  (`_triq` = triangle-quality remesh).
- Intermediate Medit `.mesh`/`.sol` files go to a temp dir (or
  `--workdir`), not committed.

### Detailed Requirements

1. **CLI** (argparse; runs under `pythonenv`):
   ```
   remesh_fault_stl.py INPUT.stl OUTPUT.stl
       [--engine {mmgs,gmsh,pymeshlab}]   # default: mmgs
       [--target-edge FLOAT]              # mmgs -hsiz; default = input median 3D edge length
       [--hausd FLOAT]                    # mmgs -hausd; default 25.0 (m)
       [--hgrad FLOAT]                    # mmgs -hgrad; default 1.3
       [--ridge-angle FLOAT]             # mmgs -ar (deg); default 45
       [--z-snap-tol FLOAT]               # |z|<tol on the top boundary snapped to 0; default 1e-3 (m)
       [--mmgs-bin PATH]                  # default: shutil.which('mmgs_O3')
       [--workdir PATH] [--keep-temp] [--verbose]
   ```
   `INPUT`/`OUTPUT` are `pathlib.Path`. Abort (exit 1) if input missing.

2. **`--target-edge` default derivation** (no hardcoded size): compute the
   **median 3D edge length** of the input STL (same edge-set dedup as
   `PLAN_retriangulate.md:reference_edge_length`). For the 500 m STL this
   resolves to ≈ 500 m, matching the fault corridor `lc_near`. Print the
   resolved value.

3. **STL → Medit `.mesh`** — function
   `stl_to_medit(stl: Path, mesh: Path) -> dict` using `meshio`
   (`meshio.read(stl)` → `meshio.write(mesh, ...)`, medit format). Return a
   stats dict (n_verts, n_tris, bbox, z-range). Note the binary-STL
   `meshio` overflow `RuntimeWarning` is benign (seen on this file); the
   read is correct.

4. **Run mmgs** — function `run_mmgs(mesh_in, mesh_out, args) -> None`
   invoking, via `subprocess.run([...], check=True)`:
   ```
   mmgs_O3 -in  <mesh_in> -out <mesh_out> \
           -hsiz <target_edge> -hausd <hausd> -hgrad <hgrad> -ar <ridge_angle> \
           [-v 5 if --verbose else -v 1]
   ```
   Rationale for each knob:
   - `-hsiz`: uniform target edge → near-equilateral isotropic mesh
     (the fault here is uniform ~500 m; uniform `-hsiz` is correct).
     If a graded fault is ever needed, switch to `-hmin/-hmax` + a `.sol`
     size field — out of scope for the 500 m target.
   - `-hausd`: **must be set explicitly** — mmgs's default is
     `0.01 × bbox-diagonal`, which for a ~100 km fault is ~1 km and would
     let the surface drift badly. 25 m keeps the fault shape tight.
   - `-hgrad 1.3`: gentle size gradation.
   - `-ar 45`: ridge-detection angle; the fault is smooth so few interior
     ridges form. Free-boundary edges (z = 0 top, NW-cut, deep perimeter)
     are auto-preserved by mmgs regardless of `-ar`.
   On non-zero exit, **report mmgs stderr verbatim and stop** (do not
   silently fall back to another engine — per project `CLAUDE.md`
   debugging rules; engine choice is the user's, via `--engine`).

5. **Medit → STL with exact z = 0 snap** — function
   `medit_to_stl_snapped(mesh_out, stl_out, z_snap_tol) -> dict`:
   - Read `mesh_out` (`meshio`).
   - For every vertex with `−z_snap_tol < z ≤ +z_snap_tol`, set `z = 0.0`
     exactly. (mmgs reprojection keeps top-boundary nodes in the z = 0
     plane up to float roundoff; this restores exactness.)
   - **Assert** `max(z) == 0.0` after the snap (raise with a clear message
     if any vertex has `z > z_snap_tol` — that would mean mmgs lifted the
     surface above z = 0, which must not happen).
   - Write binary STL via `meshio`. Return stats.

6. **`main()`** orchestration: derive defaults → `stl_to_medit` →
   `run_mmgs` (or the chosen engine) → `medit_to_stl_snapped` → print a
   before/after summary (n_tris in/out, target edge, hausd, z-range out,
   and a reminder to run Phase 2 + the quality gate).

7. **`--engine gmsh` / `--engine pymeshlab`** stubs may raise
   `NotImplementedError("see Phase 3")` in this phase; they are *specified*
   in Phase 3 and only implemented if Phase 2 shows mmgs misses the target.

### Interfaces
- Module-level functions exposed for tests: `median_edge_length`,
  `stl_to_medit`, `run_mmgs`, `medit_to_stl_snapped`, `main`.
- Module constants: `DEFAULT_HAUSD = 25.0`, `DEFAULT_HGRAD = 1.3`,
  `DEFAULT_RIDGE_ANGLE = 45.0`, `DEFAULT_Z_SNAP_TOL = 1.0e-3`.

### Edge Cases to Handle
- **mmgs binary not found**: clear error naming the expected `pythonenv`
  path; exit 2.
- **Non-manifold / self-intersecting input** (mmgs aborts): surface the
  mmgs message and stop; suggest `--engine cgal` (Phase 3 C) which is more
  tolerant, or pre-repair. Do not auto-pivot.
- **mmgs splits the surface or removes the z = 0 boundary**: detected in
  Phase 2 (extract_z0_trace finds 0 segments, or >1 merged surface). The
  Phase-1 tool itself checks: output must be a single connected component
  and must still have free-boundary edges in the z = 0 plane; warn if not.
- **`max(z) > z_snap_tol` after remesh**: hard error (free-surface
  invariant would be violated downstream).
- **Output tri count wildly different** (< 0.5× or > 2× input): warn (not
  fatal) — large change hints the target edge or hausd is off.

### Acceptance Criteria
- [ ] `python remesh_fault_stl.py <500m STL> <500m _triq STL>` exits 0 and
      writes the `_triq.stl`.
- [ ] Output STL `max(z) == 0.0` exactly; `min(z)` within `hausd` of the
      input `min(z) = −16607.6`.
- [ ] Output is a single connected surface (one shell) and still has free
      edges lying in the z = 0 plane (≥ as many trace segments as the
      input had).
- [ ] Standalone tri-quality of the output STL (computed by a small inline
      check reusing the `tri_q`/`min_angle` formulas from
      `locate_fault_slivers.py:tri_quality`): worst `min_angle ≥ 25°`,
      `tri_qmin ≥ 0.5`. (This is the *surface* gate; the *embedded-mesh*
      gate is Phase 2.)
- [ ] Hausdorff distance output↔input ≤ `--hausd` (spot-check via a
      nearest-point sampling; mmgs guarantees this by construction, the
      check confirms it).

### Dependencies
- Depends on: nothing (existing STL is the input).
- Required by: Phase 2.

## Phase 2: Regenerate a NEW `.msh` via the existing mesher and verify

### Goal
After this phase a **new** `.msh` exists, built by the *existing*
`run_z0cut_meshing.py` from the `_triq` STL, in which the **fault**
triangles meet the quality target *and* the **bulk** still passes Q1/Q2
and all invariants.

### Files to Create
- `meshing/results/msh/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh`
  — new mesh (and its auto-written `_triq_{bulk,fault}.vtu` if the driver
  emits them).
- `meshing/code/test_fault_triangle_quality.py` — new pytest module that
  gates fault-triangle quality (skip-on-missing-mesh).

### Files to Modify
- None. `run_z0cut_meshing.py` is **run**, not edited.

### Detailed Requirements

1. **Re-mesh command** (reproduces the lcfar3000 500 m recipe; only `--stl`
   and `--out` differ from the canonical call — both point at new files):
   ```bash
   PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
   cd meshing/code
   $PY run_z0cut_meshing.py \
       --stl ../results/stl_nwcut/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_500m_clean_clip_nwcut_triq.stl \
       --out ../results/msh/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh
   ```
   (Defaults reproduce `lc_near=1500`, `lc_far=3000`, pads, algorithms.
   Do **not** pass overrides unless Phase 2 verification requires it; if
   any override is needed, document why here.)

2. **Bulk gate** (existing tool, unchanged):
   ```bash
   $PY check_mesh_quality.py ../results/msh/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh
   ```
   Must report `Q1 ✓` (`tet_emin ≥ 100`) and `Q2 ✓` (`η_min > 0.1`), with
   `η_min` **no worse** than the baseline 0.184 within run-to-run mesher
   noise (a small change is acceptable; a drop below 0.1 is a regression).

3. **Fault gate** (new test module + existing locator):
   - `test_fault_triangle_quality.py` parametrizes a single assertion over
     the new mesh: load Physical-Surface-101 tris via
     `check_mesh_quality._extract_blocks` + `_tri_metrics`, assert
     `tri_qmin ≥ 0.5` and (recompute min-angle) worst `min_angle ≥ 25°`.
     Skip with `pytest.skip("mesh not built")` if the `.msh` is absent.
   - Cross-check with the existing spatial locator (worst-mesh-wide block):
     ```bash
     python3 locate_fault_slivers.py \
         ../results/msh/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh \
         --target 448338 3798440 -249.7 --radius 60000
     ```
     The "5 worst min-angle fault triangles mesh-wide" block must show all
     five `≥ 25°`.

4. **Invariant checks** (must all hold on the new mesh):
   - `mesh_zmax == 0.0` exactly (README one-liner).
   - Fault embedded: every Physical-Surface-101 triangle is a face shared
     by exactly two tets (reuse the embedding snippet referenced in
     `PLAN_mesh_quality.md` testing strategy / `EXPLORE_tmop` §8 Step C).
   - `run_z0cut_meshing.py` reported `volume tets > 0` (non-empty volume)
     and a non-trivial weld count (`welded ~N` ≈ number of trace nodes).

### Interfaces
- New pytest module `test_fault_triangle_quality.py` reusing
  `check_mesh_quality` internals (no new metric code).

### Edge Cases to Handle
- **Empty volume / weld = 0** from the mesher: means the `_triq` STL's
  z = 0 boundary is not a clean free-edge loop → go back to Phase 1
  (loosen `--hausd`, raise `--ridge-angle`, or check the z-snap). Report
  the mesher's exact output; do not silently retry.
- **Bulk η_min regresses below 0.1** near the trace: the remesh changed
  the near-surface fault edge sizing; re-run Phase 1 with `--target-edge`
  closer to the original median, or apply the existing `lc_min` floor
  semantics — but note the bulk already had 0.184 margin, so this is
  unlikely.
- **The 2 geometry-forced tris at the NW footprint tip remain < 25°**
  (measured: apex ≈ 9°, at `(365012, 3838896, z ≈ −3500 … −4300 m)`).
  These are forced by the cut footprint's sharp wedge — *not* fixable by
  remesh or refinement. Present three options to the user and let them
  choose: **(a)** trim/blunt the footprint tip (a geometry change — must
  preserve the rest of the NW-cut footprint and stay within the Hausdorff
  budget), **(b)** accept a relaxed local bar (e.g. report-and-allow ≤ 2
  tris ≥ 9°), or **(c)** leave them as-is (~4 km deep, far NW, outside the
  seismogenic zone → physically harmless). Do **not** distort the footprint
  silently.

### Acceptance Criteria
- [ ] New `_triq.msh` exists; `check_mesh_quality.py` → `Q1 ✓ Q2 ✓`.
- [ ] `tri_qmin ≥ 0.5` (stretch 0.6) and worst fault `min_angle ≥ 25°`
      (stretch 30°) on the new mesh.
- [ ] `tri_q_hist` buckets `[0,0.1)` and `[0.1,0.3)` are **0**; `[0.3,0.5)`
      is 0 (or only documented irreducible-corner tris).
- [ ] `mesh_zmax == 0.0` exactly; fault embedded (all 101-tris = 2-tet
      faces); single fault surface; non-empty volume.
- [ ] `pytest -q test_fault_triangle_quality.py` passes; on the *old*
      mesh the same gate **fails** (confirms the gate has teeth).
- [ ] Bulk `η_min` not regressed below the existing 0.184 by more than
      mesher run-to-run noise (and never below 0.1).

### Dependencies
- Depends on: Phase 1.
- Required by: nothing (Phase 3 is conditional fallback).

## Phase 3: Fallback engines (only if mmgs misses the Phase-2 targets)

### Goal
If the mmgs path (Phases 1–2) cannot reach min-angle ≥ 25° / tri_q ≥ 0.5
on the embedded mesh (most likely failure: trace/corner slivers it won't
collapse, or geometry drift), implement one alternative engine in the
*same* `remesh_fault_stl.py` (`--engine`) or as a sibling tool, and
re-run Phase 2. Pick **one** based on the observed failure.

### Option A — gmsh reclassify + remesh (`--engine gmsh`, no build)
Re-triangulate the discrete STL with gmsh's own Frontal-Delaunay mesher:
```python
gmsh.merge(stl)
gmsh.model.mesh.classifySurfaces(angle=40*pi/180, boundary=True,
                                 forReparametrization=True, curveAngle=pi)  # one patch
gmsh.model.mesh.createGeometry()
# size field: constant MathEval = target_edge (no fault-distance grading needed on the sheet)
gmsh.option.setNumber("Mesh.Algorithm", 6)   # Frontal-Delaunay (matches the pipeline)
gmsh.model.mesh.generate(2)
gmsh.write(stl_out)   # then snap z=0 + ensure single surface
```
- Pro: same trusted mesher; emits clean near-equilateral tris.
- Watch: `classifySurfaces` may split a curved fault into multiple patches
  (then the z=0 boundary curve is fragmented). Tune `curveAngle=pi` to
  force a single patch; verify one surface + a clean z=0 boundary curve
  before writing. Snap z to 0 as in Phase 1 step 5.

### Option B — pymeshlab VCG remesh (`--engine pymeshlab`, no build)
```python
import os; os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE"); os.environ.setdefault("OMP_NUM_THREADS","1")
import pymeshlab as ml
ms = ml.MeshSet(); ms.load_new_mesh(str(stl_in))
ms.meshing_isotropic_explicit_remeshing(
    targetlen=ml.AbsoluteValue(target_edge), iterations=5,
    selectedonly=False, reprojectflag=True)   # boundary preserved by VCG border handling
ms.save_current_mesh(str(stl_out))
```
- Pro: one call. Con: coarser control; confirm border vertices stay in
  z = 0 (then snap). See `EXPLORE_pymeshlab_intersections.md` for the
  OMP-init caveat and filter-name notes.

### Option C — CGAL `isotropic_remeshing` with constrained z=0 trace (`cgal-61`, C++ build)
Highest quality / tightest constraint control. New dir
`meshing/code_preprocess/remesh_cgal/` with `remesh_fault.cpp` +
`CMakeLists.txt` (built in `cgal-61`):
- Read STL → `Surface_mesh`.
- Mark **constrained edges** = all border (free) edges (z = 0 trace + NW-cut
  + perimeter) via a `dynamic_edge_property_t<bool>` map.
- `PMP::isotropic_remeshing(faces(m), target_edge, m,`
  `parameters::edge_is_constrained_map(ecm).protect_constraints(false)`
  `.collapse_constraints(true).number_of_iterations(5).do_project(true))`.
- Write **OFF/PLY** (not STL — avoid float-coord drift on the z=0 trace;
  `EXPLORE_cgal_corefine.md` "Floating-point drift"), then convert to STL
  via `meshio` and snap z.
- Alternative: `surface_Delaunay_remeshing` with the z=0 trace as
  `polyline_constraints` + `protect_constraints(true)` for a Delaunay
  quality bound. Detailed API in `EXPLORE_cgal_corefine.md`.

### Acceptance Criteria (whichever option is chosen)
- [ ] Same as Phase 2, on a re-generated `_triq.msh`.
- [ ] The chosen engine is wired behind `--engine` (or a documented
      sibling tool); the default (`mmgs`) is unchanged.

### Dependencies
- Depends on: Phases 1–2 having been run and shown insufficient.
- Required by: nothing.

## Alternatives considered (documented, not part of the primary path)

- **TMOP in-place polish of the named `.msh`** (no STL rebuild). MFEM's
  `mesh-optimizer`/`pmesh-optimizer` can move the *shared* fault nodes
  (and their tet neighbours) to raise min-angle while keeping conformity,
  and surface-fitting can slide trace nodes tangentially on the fault and
  on z = 0. This directly targets the user's named file and never rebuilds
  the mesh. **But** TMOP cannot change connectivity, so it cannot collapse
  a clip needle whose vertices are pinned on the trace; it is a *polish*,
  not a cure, and carries the boundary-attribute gotcha (tags 101–104 ∉
  {1,2,3,4}). Full recipe/limits: `EXPLORE_tmop_mesh_optimization.md`.
  Could be layered *after* Phase 2 as a final smoother if a few residual
  near-target tris remain.
- **Rebuild from raw `.ts` via `PLAN_retriangulate.md`** (drop z > 0 +
  Delaunay-retriangulate the surviving cloud). This produces a clean z = 0
  boundary from the start (no clip needles), but it re-runs the *upstream*
  STL-prep pipeline (`ts_to_stl.py` → `clean_freesurface_mesh.py` →
  `nw_cut_strip.py`, all archived) and is a larger change. It conflicts
  with "do not modify existing code" more than a targeted STL remesh, so
  it is **not** the primary path — but it is the cleanest long-term fix if
  the project later regenerates the STL family.

## Testing Strategy

| Layer | Test | Where | When |
|---|---|---|---|
| Surface quality (pre-mesh) | worst min-angle ≥ 25° / tri_qmin ≥ 0.5 on the `_triq` STL | inline in `remesh_fault_stl.py` (reusing `locate_fault_slivers.tri_quality`) | Phase 1 |
| Geometry fidelity | Hausdorff(output, input) ≤ `--hausd`; `max(z)==0` | inline in `remesh_fault_stl.py` | Phase 1 |
| Bulk Q1/Q2 | `tet_emin ≥ 100`, `η_min > 0.1` | `check_mesh_quality.py` (existing) | Phase 2 |
| Fault quality (embedded) | `tri_qmin ≥ 0.5`, worst `min_angle ≥ 25°` | new `test_fault_triangle_quality.py` + `locate_fault_slivers.py` | Phase 2 |
| Free surface | `mesh_zmax == 0` exactly | README one-liner | Phase 2 |
| Embedding | every 101-tri = 2-tet face | embedding snippet | Phase 2 |
| Gate has teeth | the fault gate **fails** on the *old* mesh | `test_fault_triangle_quality.py` (parametrize old + new) | Phase 2 |

## Risk Assessment

| Risk | Likelihood | Mitigation |
|---|---|---|
| mmgs lifts/rounds the z = 0 top boundary off the plane → empty volume on re-mesh | Medium | Explicit z-snap (Phase 1 step 5) + hard `max(z)==0` assert; Phase 2 weld-count / non-empty-volume check catches it immediately. |
| mmgs drifts the fault geometry (default `-hausd` ~1 km) | Medium-High **if `-hausd` left default** | Plan mandates explicit `-hausd 25`; Phase 1 Hausdorff check confirms. |
| `classifySurfaces` (Fallback A) splits the curved fault into patches, fragmenting the z=0 curve | Medium | `curveAngle=pi` to force one patch; verify single surface + one clean z=0 boundary curve before writing. |
| 2 geometry-forced tris at the NW footprint tip (apex ≈ 9°) can't reach 25° by remesh OR refinement | **Confirmed** (measured 2026-05-26) | Not a remesh failure — a sharp-corner floor. Present the 3 options (trim tip / relaxed bar / leave) to the user; default to "leave" (deep, far NW, outside seismogenic zone). Do not silently distort the footprint. |
| Bulk η_min regresses near the re-sized trace | Low | Bulk already at 0.184 (margin to 0.1); keep `--target-edge` ≈ input median (≈500 m). |
| meshio binary-STL overflow `RuntimeWarning` mistaken for an error | Low | Documented benign on this file; the read is correct (verified 2026-05-26). |
| CGAL STL float I/O drifts the z=0 trace coords (Fallback C) | Low (only if C used) | Use OFF/PLY round-trip, convert to STL last, then snap z. |

## What is explicitly OUT of scope
- Modifying `run_z0cut_meshing.py`, `check_mesh_quality.py`,
  `locate_fault_slivers.py`, or any existing `.msh`/`.stl`.
- Regenerating the upstream STL family from raw `.ts`
  (`ts_to_stl.py`/`clean_freesurface_mesh.py`/`nw_cut_strip.py`).
- The 1000 m / 250 m / zgraded variants — the tool is parameterized to
  generalize, but acceptance is scoped to the named 500 m lcfar3000 mesh.
- Re-running the downstream velocity/stress projections onto the new mesh
  (a separate follow-up once the `_triq` mesh is accepted).
- Changing gmsh algorithm choices in the mesher (`Algorithm3D = 1` stays).
