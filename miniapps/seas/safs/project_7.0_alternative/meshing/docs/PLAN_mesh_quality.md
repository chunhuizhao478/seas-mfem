# Implementation Plan: SAFS NW-cut Mesh Quality Improvement

> **Status: ACCEPTED 2026-05-17.** Phases 1–4 implemented as designed;
> Phase 5 mesh generation completed (six target meshes regenerated;
> velocity/stress projections deferred per user direction).  Two known
> deviations from the original plan:
> 1. **Phase 2 CLI flag renamed** `--lc-min` → `--lc-floor` because
>    `--lc-min` was already taken by `build_pos`'s size-field clamp.
>    The gmsh-side variable name stays `lc_min` so the .geo is
>    unchanged from the plan.
> 2. **Phase 3 `Mesh.OptimizeNetgen` flipped to off-by-default** because
>    it SIGBUSes on this embedded-fault geometry (Risk Assessment
>    item #1 materialized).  Driver flag flipped from
>    `--no-optimize-netgen` (opt-out) to `--optimize-netgen` (opt-in).
>    `Mesh.Smoothing` default bumped from 5 → 10 to compensate.
> 3. **Q2 (η > 0.1) not universally achieved.**  Best variant
>    (2000m_zgraded) sits at η_min = 0.0985 with 1 residual sliver;
>    other Q1-passing variants have 1–2 slivers at η ≈ 0.03–0.09.
>    The 500m variants additionally miss Q1 (tet_emin = 38.72 m vs
>    100 m floor) because the cleaned 500 m STL has ~112 m fault
>    triangle edges that propagate despite the `lc_min` floor.  The
>    user accepted this state on 2026-05-17.  See README §
>    "Reproducibility — accepted mesh set" for the canonical snapshot
>    chart.

## Overview

Bring the `zgraded` and `lcfar3000` SAFS NW-cut mesh families up to two
hard quality bars: **(Q1) minimum bulk tet edge ≥ 100 m** and **(Q2)
Joe-Liu η > 0.1 for every bulk tet**.  Current state for all six meshes
(500 / 1000 / 2000 m × {zgraded, lcfar3000}) violates both bars: min
bulk edges fall to 11–64 m and 0.05–0.10 % of tets are slivers
(η < 0.1).  Spatial forensics on `1000m_zgraded` show **100 % of slivers
and 85 % of short-edge tets live in z ∈ [−700, 0] m AND distance-to-fault
< 1.5 km** — a single defect region where the fault corridor, the
basin-surface refinement step, and the 1 m fault-to-free-surface
clearance pinch the mesher.  The plan is therefore tightly scoped to
that region; we do **not** rebuild the upstream STL pipeline or change
the .geo's algorithmic backbone.

## Constraints

### Physics / geometry (cannot change)
- **Free-surface invariant**: `mesh_zmax = 0` exactly.
- **Fault fully embedded**: every Physical Surface 101 triangle must
  remain a face of exactly two bulk tets (interior fault).
- **Fault essentially flush against the surface**: SAFS is a
  surface-rupturing fault zone; the fault-tip-to-free-surface clearance
  must stay << SAFS dimensions (~100 km).  A clearance up to ~250 m is
  physically acceptable (≪ smallest seismogenic length).  A clearance
  ≥ 500 m would suppress co-seismic surface slip and is **out of bounds**.
- **Physical tags** preserved: fault=101, top=102, bottom=103,
  sides=104, rock=1.
- **NW-cut footprint** preserved: do not regenerate `stl_nwcut/` STLs.

### Pipeline (cannot change)
- Driver is `run_nwcut_meshing.py`; template is
  `safs_fault_box_nwcut.geo`; both stay in
  `miniapps/seas/safs/project_7.0_alternative/meshing/code/`.
- Mesh outputs land in `meshing/results/msh/<stem>.msh` and matching
  `meshing/results/vtu/<stem>_{bulk,fault}.vtu`.
- All commands runnable under the `pythonenv` conda environment.
- Existing CLI flags (`--lc-far`, `--z-graded`, `--res`, `--suffix`,
  `--fault-top-clamp`, etc.) keep their current semantics.

### Numerical conventions
- Quality metric is **Joe-Liu** η = 12·(3V)^(2/3) / Σℓ² normalized so
  η = 1 on a regular tet.  This is the same metric reported by the
  driver's summary block (gmsh's internal `gamma`, which is equivalent
  to Joe-Liu up to scaling).
- "Minimum edge" means the global minimum across the six pairs of
  edges of every tet — `tet_edges.min()` in numpy parlance.
- "Sliver" means η < 0.1; matches the report buckets already used by
  the driver and previous quality charts.

## Diagnosis (Pre-Phase 0; already completed, do not re-do)

Spatial breakdown of slivers in `1000m_zgraded.msh` (1.30 M tets):

| Region | Sliver count | Short-edge (<100 m) count |
|---|---:|---:|
| z ∈ [−0.2, 0] km, d_fault < 1.5 km | 762 (98.8 %) | 29 (85.3 %) |
| z ∈ [−1, −0.2] km, d_fault < 1.5 km | 9 (1.2 %) | 5 (14.7 %) |
| z < −1 km, anywhere | 0 | 0 |
| d_fault > 1.5 km, anywhere | 0 | 0 |

**Root cause**: the fault top sits at z = −1 m (=`fault_top_clamp`), the
free surface at z = 0.  In the `zgraded` family the basin-surface
step requests `lc = 500 m` for z > −500 m, yet only 1 m of vertical
headroom is available between the fault and the box top.  Gmsh's
Delaunay 3-D mesher must fit cells into a wedge ~1 m thick × ~500 m
wide × ~1500 m long, which collapses into slivers along the
fault-corridor / free-surface intersection.  In the `lcfar3000` family
the same wedge exists with `lc = lc_near = 1500 m`, still ≫ 1 m
headroom, so the same defect appears (slightly less surface-density
related, slightly more lc_near-corridor related).

Implication: changes to gmsh smoothing alone will not lift min-edge
above 100 m (Optimize moves nodes within the existing topology).  We
need a **geometric** fix (more clearance between fault and free
surface) plus a **mesh-size floor**, plus a final **optimize pass**.

## Phase 1: Geometric headroom — raise `--fault-top-clamp` default

### Goal
After this phase, the fault top sits at z = −H below the free surface
with H large enough that gmsh has a full "near"-class cell of vertical
headroom.  H is configurable per family but defaults respect the
physics ceiling (H ≤ 250 m).

### Files to Modify

- `miniapps/seas/safs/project_7.0_alternative/meshing/code/run_nwcut_meshing.py`
  - Raise the `--fault-top-clamp` **default** from `1.0` to **`100.0`**
    (metres). Update the argparse help text to explain that this is
    the gmsh-headroom budget and that the absolute physics ceiling is
    250 m.
  - Add a new CLI flag `--fault-top-clamp-policy` with choices
    `{fixed, auto}` (default `fixed`).  When `auto`, the driver
    computes `fault_top_clamp = max(100.0, min(250.0, surface_lc))`
    where `surface_lc` is:
    - For `--z-graded`: `z_lc_top` (default 500 → clamped to 250 by
      the cap).
    - Otherwise: `lc_near` (default 1500 → clamped to 250).
    The intent: in `auto` mode the clamp matches the local target LC up
    to the physics ceiling.
  - Pass the resolved value into `clean_freesurface_mesh.py`'s zclamp
    write step **unchanged** — the existing `_zclamp<N>m.stl` filename
    convention already encodes the value as `<N>m`, so larger H simply
    writes `..._zclamp100m.stl`.
- `miniapps/seas/safs/project_7.0_alternative/meshing/README.md`
  - Update the "⚠️ Free-surface invariant" callout to explain that
    `fault_top_clamp` now defaults to 100 m (not 1 m).  Document the
    physics ceiling (250 m) and the new `--fault-top-clamp-policy auto`
    flag.

### Files to Create
None.

### Detailed Requirements

1. Argparse change:
   ```python
   ap.add_argument("--fault-top-clamp", type=float, default=100.0,
                   dest="fault_top_clamp",
                   help="[m] snap any fault vertex with z > -<this> down "
                        "to z = -<this> before passing to gmsh "
                        "(default: 100). Provides gmsh enough vertical "
                        "headroom near the free surface to avoid slivers. "
                        "Physics ceiling: 250 m (above which SAFS "
                        "co-seismic surface rupture is suppressed).  "
                        "Set 0 to disable clamping (only safe if you "
                        "also set --pad-top > 0).")
   ap.add_argument("--fault-top-clamp-policy", choices=("fixed","auto"),
                   default="fixed",
                   help="'fixed' uses --fault-top-clamp verbatim; "
                        "'auto' computes "
                        "clamp = clip(surface_lc, 100, 250) where "
                        "surface_lc = z_lc_top in --z-graded mode, "
                        "else lc_near.")
   ```
2. Auto-policy resolution (call site is the body of `main()` just
   after `args = ap.parse_args()`):
   ```python
   if args.fault_top_clamp_policy == "auto":
       surface_lc = (args.z_lc_top if args.z_graded
                     else args.lc_near)
       args.fault_top_clamp = max(100.0, min(250.0, float(surface_lc)))
   ```
3. The existing `pad_top` default formula
   `pad_top = args.pad_top if not None else args.fault_top_clamp`
   continues to hold; with `fault_top_clamp = 100`, `pad_top = 100`
   and `mesh_zmax = clamped_fault_zmax + pad_top = −100 + 100 = 0`.
   **Free-surface invariant preserved.**

### Interfaces
No public Python signatures change.  CLI surface: one default flipped,
one new flag added.  The driver's call to `clean_freesurface_mesh.py`
needs no edit — that script already accepts arbitrary `--snap-strip`
(zclamp) values.

### Edge Cases to Handle
- **User passes `--fault-top-clamp 0`**: still allowed (legacy
  "no-clamp"); gmsh PLC error returns; driver must surface that error
  exactly as today.
- **User passes `--fault-top-clamp 500`**: ABOVE physics ceiling.
  Driver must `print` a warning to stderr ("clamp 500 m exceeds
  recommended 250 m ceiling — co-seismic surface rupture will be
  suppressed") but proceed.
- **`auto` with `--z-graded` and explicit `--z-lc-top 250`**: resolves
  to `clamp = 250`; legal.
- **`auto` without `--z-graded`**: uses `lc_near = 1500`, clipped to
  250.

### Acceptance Criteria
- [ ] Default-flag invocation (`python run_nwcut_meshing.py --res 1000
      --z-graded --suffix _zgraded`) produces a mesh with
      `fault_zmax = -100 m` exactly.
- [ ] `mesh_zmax = 0` exactly (free-surface invariant preserved).
- [ ] Fault still embedded (all fault triangles match a face shared by 2
      tets; verified by the existing embedding-check snippet from the
      session log).
- [ ] No regression in non-zgraded modes: passing `--res 1000` alone
      still produces a valid mesh.

### Dependencies
- Depends on: nothing.
- Required by: Phase 2 (size-floor) and Phase 3 (optimizer) depend on
  Phase 1 to remove the geometric pinch; without it, the optimizer
  cannot eliminate slivers that the topology forces.

## Phase 2: Hard size floor inside the fault corridor

### Goal
After this phase, gmsh's local mesh size at every point inside the
fault corridor and surface zone is bounded below by an explicit floor
that defaults to 100 m, eliminating the "fast taper from `lc_near` down
to embedded-surface edge length" pinch.

### Background — why gmsh's existing `Mesh.CharacteristicLengthMin`
fails

The .geo currently sets

```gmsh
Mesh.CharacteristicLengthMin = LC_NEAR;   // = 1500 m
```

but gmsh's embedded-discrete-surface code path uses the local triangle
edge length of the merged STL as the local LC, overriding the global
`CharacteristicLengthMin` near the embedded surface.  Cleaned fault STL
edges scale with `lc_near` (median ≈ 1000 m at 1000 m resolution), but
gmsh's Delaunay 3D refines aggressively where the fault top edge
intersects the free-surface plane, producing 10–60 m edges.  A hard
**local** floor field (Field[Min]'d into the background) is the only
robust countermeasure.

### Files to Modify

- `meshing/code/safs_fault_box_nwcut.geo`
  - Add a new optional parameter `lc_min` (default `100.0`).
  - Add a new `MathEval` field (Field[10]) that reports `lc_min`
    everywhere:
    ```gmsh
    If (!Exists(lc_min))  lc_min = 100.0;  EndIf
    Field[10] = MathEval;
    Field[10].F = Sprintf("%g", lc_min);
    ```
  - Where the existing `Field[9] = Min` (or `Field[4] = Min`, or
    `Background Field = bg_field_id`) is assembled, **wrap the result
    in `Field[Max]` against Field[10]**, so the final background field
    is `max(min(distance, zgraded, sidecar), lc_min)`.  Concretely:
    ```gmsh
    // existing code computes background_raw (Field[2], Field[4], or Field[9])
    Field[11] = Max;
    Field[11].FieldsList = {background_raw, 10};
    Background Field = 11;
    ```
  - Document at the field block: "Field[10] / [11] enforce a hard
    lower bound `lc_min` (default 100 m) on the local mesh size — see
    PLAN_mesh_quality.md Phase 2."
- `meshing/code/run_nwcut_meshing.py`
  - Add CLI flag `--lc-min` (default `100.0`).
  - Forward as `-setnumber lc_min` in the gmsh command list.
- `meshing/README.md`
  - Document the new `--lc-min` flag in the "Other useful flags"
    bullet list.

### Files to Create
None.

### Detailed Requirements

1. Floor field is constant `lc_min` everywhere; it does NOT make the
   mesh coarser anywhere — Max(real_lc, lc_min) only raises sub-floor
   targets.
2. The fault corridor's `lc_near = 1500 m` is unaffected (1500 > 100).
3. The `zgraded` basin-surface step at `z_lc_top = 500 m` is
   unaffected (500 > 100).
4. The fault embedded-surface edges down to ~32 m **are** affected:
   gmsh's local refinement near the embedded fault is bounded by the
   max of (its native value, lc_min).  Effectively: no tet edge may be
   shorter than 100 m.  This is exactly the goal.
5. **Do NOT** raise `Mesh.CharacteristicLengthMin` — leave at
   `LC_NEAR`.  The floor must be applied as a Field, not as a global
   gmsh option, because the latter is overridden by embedded surfaces
   (the actual root cause of the current sliver region).

### Interfaces
- New gmsh parameter `lc_min` (scalar, metres).
- New CLI flag `--lc-min` (float, default 100.0, metres).
- No Python function signatures change.

### Edge Cases to Handle
- **`--lc-min 0`**: floor disabled; same as today.  Allowed for
  comparison/debugging.
- **`--lc-min` ≥ `lc_near`**: the floor would defeat the fault-corridor
  refinement.  Driver must `print` a warning and abort:
  "lc_min=<X> ≥ lc_near=<Y>: this would erase the fault-corridor
  refinement.  Pick lc_min < lc_near."
- **lc_min between observed pinch edges and lc_near** (the intended
  regime): no special handling; this is the operating point.

### Acceptance Criteria
- [ ] `--lc-min 100` (default) on `1000m_zgraded`: post-mesh
      verification reports `min bulk edge ≥ 100 m`.
- [ ] `--lc-min 100` on every other variant (500m/2000m × zgraded,
      500/1000/2000m × lcfar3000): same.
- [ ] Fault triangulation edges unchanged in median (we only raise sub-
      floor edges; fault median is already ≥ `lc_near` ≫ 100).
- [ ] Mesh cell count grows by ≤ 25 % vs Phase 1 mesh at the same
      resolution (floor only refines NOT-yet-refined regions; should
      coarsen the over-refined wedge near the fault top, possibly
      yielding a net reduction).
- [ ] Fault still embedded (all fault triangles match 2-tet faces).

### Dependencies
- Depends on: Phase 1 (geometric headroom — without it, the floor
  alone cannot reach 100 m because gmsh has only 1 m of vertical room
  near the fault top, forcing slivers even with the floor active).
- Required by: Phase 3 (optimizer) is a polishing step; it depends on
  Phases 1–2 to have produced a near-target mesh.

## Phase 3: Gmsh post-mesh optimization passes

### Goal
After this phase, the **quality** target (Q2: η > 0.1 for all tets) is
achieved.  Phases 1+2 should already get min-edge ≥ 100 m; this phase
removes the residual slivers via gmsh's Netgen-based optimizer and a
final Laplacian smoothing pass.

### Background — gmsh's optimization knobs

| Option | Effect | Default |
|---|---|---|
| `Mesh.Optimize` | Generic Delaunay smoothing (node placement + edge swaps) | 0 |
| `Mesh.OptimizeNetgen` | Calls Netgen's tet optimizer (most effective for slivers) | 0 |
| `Mesh.OptimizeThreshold` | Quality threshold below which optimize attacks a cell | 0.3 |
| `Mesh.Smoothing` | Number of Laplacian smoothing passes on the bulk | 1 |
| `Mesh.QualityType` | 2 = Inradius/Circumradius (gamma); 4 = ICN | 2 |

None are currently enabled by `safs_fault_box_nwcut.geo`.

### Files to Modify

- `meshing/code/safs_fault_box_nwcut.geo`
  - At the end of the file (after `Mesh.Algorithm3D = 1`), add:
    ```gmsh
    // Quality optimization (PLAN_mesh_quality.md Phase 3).
    // Netgen tet optimizer is the most effective tool against slivers;
    // gmsh's Optimize is a cheaper Laplacian-style pass; both keep the
    // embedded fault topology intact (no edge-swap is allowed to break
    // a Physical Surface 101 face).
    If (!Exists(do_optimize))      do_optimize      = 1; EndIf
    If (!Exists(do_optimize_netgen)) do_optimize_netgen = 1; EndIf
    If (!Exists(optimize_threshold)) optimize_threshold = 0.3; EndIf
    If (!Exists(smoothing_passes))   smoothing_passes   = 5; EndIf
    Mesh.Optimize         = do_optimize;
    Mesh.OptimizeNetgen   = do_optimize_netgen;
    Mesh.OptimizeThreshold = optimize_threshold;
    Mesh.Smoothing        = smoothing_passes;
    Mesh.QualityType      = 2;          // Inradius/circumradius = gamma
    ```
  - The threshold `0.3` is gmsh's `gamma`, not Joe-Liu η; the two
    correlate but are not identical.  `gamma = 0.3` is roughly
    `η ≈ 0.5`, so the optimizer will attack everything below η ≈ 0.5,
    which is well above our η > 0.1 acceptance bar.

- `meshing/code/run_nwcut_meshing.py`
  - Forward the four new gmsh-side knobs as `-setnumber`
    arguments — same pattern as the existing `pad_top`, `lc_near`
    forwarding.  Default values come from the .geo's `If (!Exists)`
    guards so explicit forwarding is optional; we add it only when the
    user passes the CLI override.
  - Add CLI flags `--no-optimize`, `--no-optimize-netgen`,
    `--smoothing-passes N`, `--optimize-threshold X`.

### Files to Create
None.

### Detailed Requirements

1. `do_optimize` and `do_optimize_netgen` default to 1.  Optimization
   is on by default; users can disable via `--no-optimize` /
   `--no-optimize-netgen` (mainly for debugging or timing comparison).
2. Smoothing passes default to 5 (gmsh's default is 1; bumping to 5
   gives noticeably better quality at <10 % extra runtime for our
   mesh sizes).
3. `Mesh.OptimizeThreshold = 0.3` chosen so that optimize attacks every
   tet with gamma < 0.3, i.e. roughly η < 0.5; this guarantees η > 0.1
   slivers are visited (they have gamma << 0.1 << 0.3).
4. `Mesh.QualityType = 2` (gamma) chosen because it is gmsh's
   sliver-friendliest metric and matches OptimizeThreshold semantics.
5. **Do not change** `Mesh.HighOrderOptimize` (P1 meshes only).
6. **Do not enable** `Mesh.Algorithm3D = 10` (HXT) — the .geo's
   existing comment notes HXT fails on this geometry ("No closed
   volume"); leave at `Mesh.Algorithm3D = 1` (Delaunay).

### Interfaces
- Four new gmsh parameters (`do_optimize`, `do_optimize_netgen`,
  `optimize_threshold`, `smoothing_passes`).
- Four new CLI flags on `run_nwcut_meshing.py`.

### Edge Cases to Handle
- **Optimizer collapses an embedded fault face**: gmsh's Netgen
  optimizer respects Physical Surface constraints by default; verify
  via the embedding-check after every variant. If a regression
  appears, fall back to `Mesh.OptimizeNetgen = 0` and rely on Mesh.
  Optimize + Smoothing alone.
- **Optimizer cost explodes at 500 m resolution**: gmsh's Netgen
  optimizer is O(N log N) in tet count; expect ≤ 2× wall-time vs the
  no-optimize baseline.  If runtime exceeds 5 min, raise threshold to
  0.2 (attack only more obvious slivers).

### Acceptance Criteria
- [ ] On every variant (500 / 1000 / 2000 m × {zgraded, lcfar3000}):
      **η_min > 0.1** for all bulk tets.
- [ ] On every variant: histogram bucket `η < 0.1` is **0**.
- [ ] Mesh wall-time increases by ≤ 100 % vs the Phase 2 baseline.
- [ ] Fault still embedded (all fault triangles match 2-tet faces).
- [ ] Free-surface invariant preserved (`mesh_zmax = 0`).
- [ ] No regression in median or mean quality (medians should improve).

### Dependencies
- Depends on: Phases 1 and 2.
- Required by: Phase 4 verification.

## Phase 4: Verification harness

### Goal
After this phase, a single command produces a pass/fail report against
both Q1 and Q2 for the six target meshes, and the matching quality
chart for the README.  The harness is what the next mesh edit will
regress against.

### Files to Create

- `meshing/code/check_mesh_quality.py` — new module + CLI.
  - Public function:
    ```python
    def check_mesh_quality(
        msh_path: pathlib.Path,
        min_edge_floor_m: float = 100.0,
        min_eta_floor: float = 0.1,
    ) -> dict:
        """Return a structured report.

        Keys:
          n_tet, n_tri_fault,
          mesh_zmax,
          tet_emin, tet_emed, tet_emax,
          eta_min, eta_med, eta_mean, eta_hist (5-bucket counts),
          tri_emin, tri_emed, tri_emax,
          tri_qmin, tri_qmed, tri_qmean, tri_q_hist,
          q1_pass: bool   # tet_emin >= min_edge_floor_m
          q2_pass: bool   # eta_min  >  min_eta_floor
          violations: list[str]
        """
    ```
  - The metric formulas must exactly match the snippets used in the
    session log so charts are comparable run-to-run:
    - Joe-Liu η: `12 * (3V)**(2/3) / sum(edge_length**2)`
    - Triangle q: `4*sqrt(3)*A / sum(edge_length**2)`
    - Histogram buckets: `[0.0, 0.1, 0.3, 0.5, 0.7, 1.0+1e-9]`
  - CLI behaviour: `python check_mesh_quality.py PATH [PATH ...]`
    prints a markdown table and exits 0 iff every input mesh passes
    both Q1 and Q2; exits 1 otherwise.  Designed to be wired into the
    existing pytest suite as a smoke test.
- `meshing/code/test_mesh_quality.py` — pytest module.
  - Two test functions:
    - `test_q1_min_edge_target_meshes()`: parametrize over the six
      target meshes; assert `q1_pass`.
    - `test_q2_eta_floor_target_meshes()`: parametrize over the six
      target meshes; assert `q2_pass`.
  - Both skip with `pytest.skip("mesh not built")` if the .msh file
    is absent (so the suite passes on a fresh checkout).

### Files to Modify

- `meshing/README.md`
  - Add a "Quality gates" subsection under "Other useful flags" that
    explains Q1 / Q2 and points at `check_mesh_quality.py` as the
    canonical regression check.

### Detailed Requirements

1. The module must use only stdlib + `numpy` + `meshio` (no extra
   deps).
2. CLI exit code must be machine-readable: `0` = all-pass, `1` = any
   fail.  The markdown table is printed regardless.
3. The `q1_pass` and `q2_pass` keys are strictly defined:
   - `q1_pass = bool(tet_emin >= min_edge_floor_m)` (≥, not >).
   - `q2_pass = bool(eta_min > min_eta_floor)` (>, not ≥).
4. `violations` is a list of human-readable strings; empty if both
   pass.

### Interfaces
- One public Python function `check_mesh_quality` returning a dict.
- One CLI entry point `python check_mesh_quality.py PATH...`.

### Edge Cases to Handle
- **Mesh has no tetra block**: raise `ValueError("mesh has no tetra
  block")`; do not silently report 0/0.
- **Mesh has no Physical Surface 101**: triangle stats fall back to
  `None`; do not fail Q1/Q2 (these are bulk-only gates).
- **Mesh has duplicate Physical Surface 101 block**: take the union
  (matches what `project_to_fault_stress.py` does — see existing
  fault-VTU reader).

### Acceptance Criteria
- [ ] `python check_mesh_quality.py results/msh/safs_fault_box_nwcut_*_zgraded.msh
       results/msh/safs_fault_box_nwcut_*_lcfar3000.msh` exits 0 after
      Phases 1–3 are applied.
- [ ] `pytest -q test_mesh_quality.py` passes after Phases 1–3.
- [ ] On a deliberately-degraded mesh (e.g. an old, pre-Phase-1
      `_zgraded.msh` re-introduced), the CLI exits 1 and identifies
      the failing variant by name.

### Dependencies
- Depends on: Phases 1–3 to produce passing meshes.
- Required by: nothing.

## Phase 5: Regenerate the six target meshes and refresh downstream artefacts

### Goal
After this phase, the six target `.msh` files are up to date with all
Phase 1–3 changes; matching `.vtu` files are regenerated; the velocity
and stress projections that consume the 1000 m_zgraded variant are
re-run; the README chart is refreshed.

### Files to Create / Modify

- `meshing/results/msh/safs_fault_box_nwcut_{500,1000,2000}m_zgraded.msh`
  — regenerate via the driver.
- `meshing/results/msh/safs_fault_box_nwcut_{500,1000,2000}m_lcfar3000.msh`
  — regenerate.
- `meshing/results/vtu/safs_fault_box_nwcut_{500,1000,2000}m_{zgraded,lcfar3000}_{bulk,fault}.vtu`
  — auto-regenerated by `msh_to_vtu.py` as part of the driver call.
- `velocity/results/<version>/preview/projected_velocity_1000m_zgraded/`
  — re-project (the existing 1.2 GB .vtu file is now stale).
- `stress/results/1000m_zgraded/` — re-project (the existing
  fault/bulk-stress VTU files are stale).
- `meshing/README.md` — replace the quality-chart block with the new
  numbers.

### Detailed Requirements

1. Build order (driver invocations under `pythonenv`):
   ```bash
   cd meshing/code
   # zgraded family (uses Phase 1+2+3 defaults end-to-end)
   python run_nwcut_meshing.py --z-graded --suffix _zgraded \
       --res 500 1000 2000 --quiet
   # lcfar3000 family
   python run_nwcut_meshing.py --lc-far 3000 --suffix _lcfar3000 \
       --res 500 1000 2000 --quiet
   ```
2. Quality gate check (must pass before proceeding):
   ```bash
   python check_mesh_quality.py results/msh/safs_fault_box_nwcut_*_zgraded.msh \
                                results/msh/safs_fault_box_nwcut_*_lcfar3000.msh
   echo "exit=$?"   # must be 0
   ```
3. Mesh re-conversion for MFEM (msh4 → msh2) and re-projection of the
   `1000m_zgraded` variant for each of the three CVM versions:
   ```bash
   gmsh ../../meshing/results/msh/safs_fault_box_nwcut_1000m_zgraded.msh \
        -format msh2 -save \
        -o /tmp/safs_msh2/safs_fault_box_nwcut_1000m_zgraded.msh -v 0
   for ver in cvmh cvm_s4.26.m01 multiscale_statewise_cvm; do
       cd /…/velocity/results/$ver/preview
       …/seas_project_velocity_to_mesh \
           --mesh /tmp/safs_msh2/safs_fault_box_nwcut_1000m_zgraded.msh \
           --sidecar ../velocity_safs.h5 \
           --out projected_velocity_1000m_zgraded
   done
   ```
4. Stress re-projection (under `pythonenv`):
   ```bash
   cd stress/code
   rm -rf ../results/1000m_zgraded
   python project_to_fault_stress.py \
       ../../meshing/results/vtu/safs_fault_box_nwcut_1000m_zgraded_fault.vtu \
       safs_fault_box_nwcut_1000m_zgraded --write-bulk
   python verify_onfault_stress.py …   # the canonical incantation from stress/README.md
   ```
5. Update `meshing/README.md`'s quality-chart block with the new
   numbers from `check_mesh_quality.py`'s table output.

### Interfaces
None; this phase only consumes the interfaces created in Phases 1–4.

### Edge Cases to Handle
- **Re-projection misses a sidecar containment guard**: only relevant
  if `mesh_zmax` shifts. Phase 1 preserves `mesh_zmax = 0`; the
  existing sidecars (built with `--extend-z-top 100`) still contain
  the mesh. No re-build of sidecars is needed.
- **Stress verifier reports any non-pass field**: STOP, do not
  proceed; investigate before tagging Phase 5 complete.

### Acceptance Criteria
- [ ] All six `_zgraded` and `_lcfar3000` meshes pass `q1_pass` and
      `q2_pass`.
- [ ] `meshing/results/vtu/` contains exactly six matching bulk + fault
      VTU pairs (no stale entries).
- [ ] Velocity projection for the three CVM versions on the new
      1000m_zgraded mesh: exit 0 each; per-field reported ranges
      finite and within the schema attrs.
- [ ] Stress verifier: 12 / 12 fields pass.
- [ ] README chart reflects the new numbers; the pre-fix chart is
      moved to an `Appendix: pre-Phase 5 quality (historical)` block.

### Dependencies
- Depends on: Phases 1–4.
- Required by: nothing.

## Testing Strategy

| Layer | Test | Where | When |
|---|---|---|---|
| Q1 / Q2 gates | Per-mesh metric assertions | `check_mesh_quality.py` + `test_mesh_quality.py` (Phase 4) | After every driver call |
| Embedding | "every fault triangle shared by 2 tets" (the snippet from this session) | New `test_fault_embedded()` in `test_mesh_quality.py` | After every driver call |
| Free-surface | `mesh_zmax == 0` | New `test_freesurface_invariant()` in `test_mesh_quality.py` | After every driver call |
| Pipeline regression | Existing `test_nw_cut_strip.py` | unchanged | Standard pytest run |
| Velocity round-trip | `seas_project_velocity_to_mesh` exit 0 + finite field ranges | Manual (Phase 5) | Once per Phase-5 cycle |
| Stress round-trip | `verify_onfault_stress.py` ≥ 12/12 pass | Manual (Phase 5) | Once per Phase-5 cycle |

Convergence tests for the SEAS solver itself are out of scope here; the
mesh quality improvement should not change the physical answer beyond
the discretization-error band, and confirming that band is a separate
investigation.

## Risk Assessment

| Risk | Likelihood | Mitigation |
|---|---|---|
| Optimizer collapses an embedded fault face (loses 2-tet sharing) | Low | Phase 3 `Mesh.OptimizeNetgen` respects Physical Surface constraints by default; Phase 4's embedding test catches any regression immediately. Fallback: disable Netgen, keep only generic Optimize + Smoothing. |
| `fault_top_clamp = 100 m` suppresses near-surface co-seismic slip in benchmarks | Medium-Low | 100 m is ≪ smallest seismogenic length (~1 km); ceiling at 250 m provides margin. Document the trade-off in README. Provide `auto` policy that picks the smallest viable clamp per family. |
| Floor field interacts unexpectedly with sidecar-driven size field (USE_SIZE_FIELD=1) | Low | The sidecar pipeline currently uses lc_near = 1500 m everywhere it activates, well above the 100 m floor. The Max(...) wrapping is monotone so cannot make any region coarser than the user requested. |
| Mesh cell count grows uncomfortably (memory / runtime) | Low | Phase 2 raises only sub-floor regions; expected NET decrease at the fault top (fewer slivers replaced by a few well-shaped cells). 500 m_zgraded is the largest mesh today (1.45 M tets); Phase 5 expected to land in 1.0–1.4 M range. |
| Netgen optimizer not present in the conda mfem-dev/pythonenv gmsh build | Medium | Detect via `gmsh -info`; if absent, fall back to `Mesh.Optimize = 1` + `Mesh.Smoothing = 10` only. Document the requirement in README. |
| Pre-Phase-1 meshes still present in `results/msh/` confuse downstream consumers | Low | Phase 5 step 1 deletes all .msh in that folder before regeneration (matches user's "clear and rebuild" pattern from the session log). |
| Plan misses a sliver hot-spot outside the diagnosed wedge (e.g. NW-cut edge under coarser settings) | Medium | Phase 4's verifier reports all violations by tet index; if a new hot-spot appears at 500 m or 2000 m resolution, log it and either tighten Phase 3's `optimize_threshold` or revisit Phase 2's floor for that variant only. |

## What is explicitly OUT of scope

- Rebuilding `stl_cleaned/` or `stl_nwcut/` STLs.
- Changing the gmsh algorithm choice (`Mesh.Algorithm3D = 1` stays).
- Touching the velocity-sidecar pipeline (`build_velocity_cvmh.py`,
  HDF5 schema, etc.).
- Re-cutting the NW-cut footprint.
- Generating size-field-driven (`--gen-size-field --sidecar …`) mesh
  variants — Phase 5 only regenerates the six target meshes.
