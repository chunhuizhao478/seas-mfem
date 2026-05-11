# Implementation Plan: Sidecar→Mesh Velocity Projection Fix (v3 G-1 + G-2)

**Date:** 2026-05-10
**Triggers:** `data_projected/comparison_sidecar_vs_mesh_*_z1km_v2.png`
(NW-SE streaks in the H1-P1 panel along the cut fault trace, basin /
basement contrast).
**Supersedes / continues:** `fault_zone_projection_plan_v3.md`. v1 / v2
are kept as design history; v2 §2 (DG-L²(0)) and v2 §3 (per-side fault
projection) remain **retracted** and are NOT in scope of this PLAN.
**User decisions captured:**
- Scope: G-1 (sidecar-driven mesh refinement) **and** G-2 (H1-P2 default).
- Acceptance: per-tet integrated L² of Vs over the `[0–500 m]`
  fault-distance band ≤ **80 m/s** on the 500 m mesh.
- Mesh budget: G-1 remesh ≤ **2× current tet count**.

---

## Overview

The current `seas_project_velocity_to_mesh` driver projects the CVM-H
sidecar onto an `H1_FECollection(1)` GF. Vertex values are bit-correct,
but in-tet linear interpolation between vertices that span the
basin/basement Vs jump produces visible streaks. v3 has shown this is
representation-error proportional to `‖∇Vs‖ × h`. The fix is to (G-1)
shrink `h` where `‖∇Vs‖` is large by driving the gmsh size field from
the sidecar gradient itself, and then (G-2) cheapen the per-tet error
further by raising the projection order from P1 to P2 inside each tet.
A new Python metric runs first so we have a single number to drive
down across iterations.

## Constraints

### Interfaces that cannot change
- `io/field_coefficient.{hpp,cpp}` public API: `Project`,
  `ProjectSerial`, `ProjectVelocity` keep their existing signatures.
  G-1 needs no C++ edits; G-2 changes only a *default* in the driver,
  not the projector API.
- HDF5 sidecar v1 schema (`code_preprocess/data_projection/sidecar.py`,
  `io/data_field_3d.{hpp,cpp}`): no schema-version bump, no new
  required attributes, no new top-level groups. The size-field tool
  reads the existing `/grid/{x,y,z}` and `/fields/Vs` *only* and
  writes a *separate* gmsh `.pos` file alongside the sidecar — not into
  the `.h5`.

### Files that must NOT be touched (C2-invariant per `miniapps/seas/CLAUDE.md`)
- `bp5/`, `bp1/`, `bp2/`, `domain/`, `solver/`,
  `friction/dieterich_ruina.hpp`, `fault/fault_basis.hpp`.
- `dynamic/heterogeneous_material.hpp`: leave existing `Mode::Constant`
  and `Mode::GridFunction` paths alone. **Do not** add
  `Mode::DG_L2_GridFunction` (v2 retracted).
- All TPV-/BP-driver source files.

Per project memory ([C2] / `feedback_dynamic_folder_editable_for_tpv104.md`)
the rest of `dynamic/` is editable, but this PLAN does not touch it —
G-1 + G-2 stay entirely in `safs/project_7.0_alternative/...` and the
single driver `drivers/project_velocity_to_mesh.cpp`.

### Conventions
- Python preprocess code: lives under
  `safs/project_7.0_alternative/code_preprocess/data_projection/`.
  Same `from __future__ import annotations` + `pathlib.Path` style as
  the existing `sidecar.py`, `raw_readers.py`,
  `plot_comparison_with_dg0.py`, `build_velocity_cvmh.py`.
- Meshing: `code_meshing/safs_fault_box_nwcut.geo` is the gmsh template
  (already parameterised via `-setnumber` /
  `-setstring`); `run_nwcut_meshing.py` is the Python orchestrator.
  New CLI flags follow the existing `--lc-near`, `--lc-far`,
  `--dist-inner`, `--dist-outer` pattern.
- C++ default-argument changes: keep ABI stable; `--order P` CLI must
  continue to accept any `P ≥ 1`.
- Conda environments per `miniapps/seas/CLAUDE.md`:
  `pythonenv` for gmsh/Python preprocess; `mfem-dev` for the C++ build.

### Numerical / acceptance constraints
- The driving metric is the **per-tet integrated L²** of `(proj − f)`
  over each tet on the `z = −1 km ± 200 m` slab, with 4-point Gauss
  tet cubature (Stroud T3:5-1, equal weights ¼, barycentric
  `(α, β, β, β)` with `α = (5 − √5)/20`, `β = (5 + 3√5)/20`). This
  formula already lives at `plot_comparison_with_dg0.py:170-191`
  (`TET_CUB_BC`, `TET_CUB_W`, `tet_cell_mean_via_cubature`) and at
  `plot_comparison_with_dg0.py:226-279`
  (`stratified_metrics_h1_vs_dg0`); Phase 1 ports it into a CLI tool
  without re-deriving.
- Stratification by distance to the cut STL using
  `scipy.spatial.cKDTree` on triangle centroids, into bands
  `[0, 500), [500, 1500), [1500, 3000), [3000, ∞) m`. Same bands as
  `plot_comparison_with_dg0.py:268`.
- Acceptance: `[0, 500)` band per-tet L² (m/s) on the `Vs` field at
  `z = −1 km ± 200 m`:
  - Baseline (current H1-P1, 500 m mesh): **148** (from v3 §0).
  - **G-1 alone target: ≤ 80 m/s** on the same baseline mesh budget
    (≤ 2×) — this is the user-confirmed acceptance gate.
  - G-2 stacked (after G-1) target: ≤ 50 m/s.
- Bulk band `[3000, ∞)` per-tet L² must NOT regress beyond noise floor.
- Sidecar interpolation-only invariant
  (`io/field_coefficient.cpp:144-153, 191-196`): the post-G-1 mesh's
  bbox must remain strictly contained in the sidecar bbox, including
  any new cells the size field places in basin-edge regions. The
  remesh script must verify this before declaring success (mirrors
  the existing `--sidecar` containment guard in
  `run_nwcut_meshing.py:243-278`).

---

## Phase 1: `fault_zone_metric.py` — per-tet L² CLI

### Goal
Produce a single Python script that ingests one mesh-projected `*.vtu`
plus the sidecar `.h5` and prints the per-tet integrated L² of `Vs`
stratified by distance-to-fault, in the same units (m/s) and bands
that `plot_comparison_with_dg0.py` reports — *without* requiring
matplotlib or generating PNGs. This is the metric we drive G-1 and
G-2 against in subsequent phases.

### Files to Create
- `safs/project_7.0_alternative/code_preprocess/data_projection/fault_zone_metric.py`
  — CLI tool. ≈ 250 LOC. Importable as a module so Phase-2 / Phase-3
  acceptance tests can call it directly without a subprocess.

### Files to Modify
- *(none)*. The script lives next to the existing
  `plot_comparison_with_dg0.py` and reuses `parse_mfem_vtu.parse_vtu`
  from `/tmp/parse_mfem_vtu.py` exactly the way `plot_comparison_with_dg0.py:35`
  does (preserved to avoid silently changing behaviour during the
  metric port). If `/tmp/parse_mfem_vtu.py` is not present in CI,
  detect that case and fail with a clear "missing `parse_mfem_vtu`;
  copy it to `/tmp/` or pass `--parse-mfem-vtu PATH`" message.

### Detailed Requirements

1. **CLI surface** (`argparse`):
   ```
   python fault_zone_metric.py
       --vtu        PATH         # required — mesh-projected .vtu
       --sidecar    PATH         # required — schema-v1 .h5
       --stl        PATH         # required — cut fault STL
       --field      Vs           # default Vs (only Vs is wired today)
       --z-target   -1000.0      # default −1000 m
       --half-thick 200.0        # default 200 m (matches v3 measurement)
       --bands      0,500,1500,3000,inf   # default
       --json       PATH         # optional: write a JSON record
       --parse-mfem-vtu PATH     # optional: override the /tmp helper
   ```

2. **Public function** for in-process callers (used by Phase-2 / Phase-3
   acceptance tests):
   ```python
   def per_tet_l2_by_band(
       vtu_path: pathlib.Path,
       sidecar_path: pathlib.Path,
       stl_path: pathlib.Path,
       *,
       field: str = "Vs",
       z_target: float = -1000.0,
       half_thick: float = 200.0,
       bands: tuple[tuple[float, float], ...] = (
           (0.0, 500.0), (500.0, 1500.0),
           (1500.0, 3000.0), (3000.0, math.inf)),
   ) -> list[BandRow]
   ```
   where `BandRow` is a `dataclass(frozen=True)` with fields
   `(lo: float, hi: float, n_tets: int, l2_proj: float)`. Units of
   `l2_proj` are the same as `field` (m/s for Vs).

3. **Algorithm** (mirrors `plot_comparison_with_dg0.py:226-279`):
   - Parse `.vtu` → vertex coordinates `pts (Nv, 3)`, tet
     connectivity `conn (Ne, 4)`, vertex-DOF `Vs` array
     `vert_vals (Nv,)`. The driver currently emits H1-P1, so vertex
     values are exactly the projector's output. **In Phase 3 (G-2)**
     the driver becomes H1-P2; the script must transparently handle
     order > 1 — see Requirement 7.
   - Read `gx, gy, gz, F` from `/grid/{x,y,z}` and
     `/fields/<field>` of the sidecar (h5py, no NaN expected per
     v1 schema).
   - Slab mask: `ec = pts[conn].mean(1)`;
     `mask = (ec.z >= z_target − half_thick) & (ec.z ≤ z_target +
     half_thick)`.
   - Cubature points / weights: import the named constants
     `TET_CUB_BC`, `TET_CUB_W` from `plot_comparison_with_dg0.py`
     (single source of truth — do not redefine them in the CLI).
   - Trilinear sidecar evaluation: import `trilinear` from the same
     module (it is already pure-numpy and side-effect free).
   - For each tet on the slab:
     - 4 cubature points in physical coords:
       `qpts = einsum("qi,tij->tqj", TET_CUB_BC, V_xyz)`.
     - "Truth" values: `f_truth = trilinear(gx, gy, gz, F,
       qpts.reshape(-1,3)).reshape(n_t, 4)`.
     - Projection values at the cubature points: see Requirement 7
       (P1 closed-form; P2 via the script's own `mfem_eval`-free
       in-tet evaluator).
     - `res = proj_at_q - f_truth`;
       `l2_per_tet = sqrt((res**2 * TET_CUB_W).sum(axis=1))`.
   - STL distance: build `cKDTree` on triangle centroids
     (same as `plot_comparison_with_dg0.py:300-305, 266-267`); query
     each tet centroid for nearest distance. Bin into bands.
   - Per-band aggregate: `sqrt(mean(l2_per_tet[band]**2))`.

4. **JSON record schema** (when `--json PATH` is given):
   ```json
   {
     "vtu":          "...",
     "sidecar":      "...",
     "stl":          "...",
     "field":        "Vs",
     "z_target_m":   -1000.0,
     "half_thick_m": 200.0,
     "bands": [
       {"lo_m": 0.0,    "hi_m": 500.0,    "n_tets": 1234, "l2": 148.3},
       ...
     ],
     "global": {"n_tets": 7890, "l2": 198.7}
   }
   ```
   Stable schema — Phase 2's accept-gate Python in
   `run_nwcut_meshing.py` reads this verbatim.

5. **Stdout layout** must match the existing
   `plot_comparison_with_dg0.py:469-486` table verbatim so logs are
   diffable across the migration. Header:
   ```
   per-tet integrated L2 |proj - source|  (m/s):
       band (m)       n_tets       proj
       [0-500]            xxx     xxx.x
       [500-1500]         xxx     xxx.x
       [1500-3000]        xxx     xxx.x
       [3000-inf]         xxx     xxx.x
       GLOBAL             xxx     xxx.x
   ```

6. **No matplotlib import.** This is intentional: the script must
   import in environments that lack the GUI stack (e.g., a thin
   Frontera login node). PNGs stay in `plot_comparison_with_dg0.py`.

7. **H1-Pp evaluation at cubature points (Requirement 3 sub-bullet).**
   - **P1 case (`order = 1`, `dofs_per_elem = 4`):** vertex-only
     barycentric, identical to `plot_comparison_with_dg0.py:248-249`:
     `proj_at_q = einsum("qi,ti->tq", TET_CUB_BC, vert_vals[conn])`.
   - **P2 case (`order = 2`, `dofs_per_elem = 10`):** the P2 GF stored
     in `proc000000.vtu` is written by MFEM with `SetHighOrderOutput
     (true)` (driver line 115). Its connectivity is therefore a Lagrange
     tet of order 2 with 10 DOFs per element (4 corner + 6 edge mid).
     The in-tet evaluator computes the standard 10 nodal Lagrange
     basis at each cubature point's *barycentric* coordinates and
     dot-products with the per-element DOF vector. Closed form for
     barycentric `(λ0, λ1, λ2, λ3)`:
     ```
     phi_corner_i =      λ_i (2 λ_i − 1),    i ∈ {0,1,2,3}
     phi_edge_ij  = 4    λ_i  λ_j,           (i,j) ∈ MFEM order
     ```
     The MFEM Pᵏ tet edge ordering on the reference simplex is
     `(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)`; verify against
     a single-tet unit fixture (Test T-1-3) before trusting on the
     real mesh.
   - **Auto-detection:** read `connectivity_lengths` (or the cell
     types) from the .vtu. If `dofs_per_elem == 4`, P1 path; if
     `dofs_per_elem == 10`, P2 path; otherwise fail with a clear
     "unsupported order" message naming the observed cell-DOF count.

8. **Performance.** The 500 m mesh has ≈ 580 k tets in the full domain
   and ≈ 30 k on the `z = −1 km ± 200 m` slab. The vectorised numpy
   path in `plot_comparison_with_dg0.py` runs in < 5 s on that data;
   the CLI port must stay within 2× of that on the same hardware.

### Interfaces

- Public Python: `per_tet_l2_by_band(...)` returning
  `list[BandRow]`; called from Phase 2 / 3 tests (R-2-G-T2,
  R-3-G-T2 below).
- CLI exit code: 0 on success; 2 on bbox containment violation
  (mirrors the v1 strict-interpolation invariant —
  `io/field_coefficient.cpp:148-153`); 1 on any other error.

### Edge Cases to Handle
- Mesh has no tets on the slab → write `n_tets = 0`, `l2 = nan` for
  every band, exit 0 (so the CI step doesn't fail when an upstream
  remesh accidentally moved the slab out of the mesh — the table makes
  the cause obvious).
- Sidecar's bbox does not contain the mesh footprint at the slab z
  → exit 2, print the same "data bbox vs mesh bbox" formatted table
  the C++ side prints (`io/field_coefficient.cpp:101-115`).
- `field` requested but absent from `/fields/`: list the available
  fields and exit 1.
- P2 vtu emits collapsed coordinates for boundary-edge nodes
  (degenerate tets): warn-and-skip those tets, count the skip in the
  output, fail only if > 1 % of slab tets are skipped.

### Acceptance Criteria
- [ ] `python fault_zone_metric.py --vtu data_projected/preview/projected_velocity_500m/Cycle000000/proc000000.vtu --sidecar data_projected/velocity_safs.h5 --stl data_cutnwfault/SAFS-...-ALT6_500m_clean_clip_nwcut.stl`
  prints the same `[0–500 m]` band L² as v3 §0 reports for H1-P1
  on the 500 m mesh, **148 m/s ± 1 m/s**.
- [ ] Same command on 1000 m / 2000 m mesh prints **207 / 238 m/s ±
  1 m/s** (matching v3 §0 table).
- [ ] `--json` flag writes a parseable JSON whose `bands[0].l2` equals
  the printed `[0–500]` value.
- [ ] Importing `per_tet_l2_by_band` from another script returns
  `BandRow` objects with the same numerics.
- [ ] **T-1-1 (P1 unit fixture):** synthetic 8-tet unit cube with a
  closed-form linear field `Vs(x,y,z) = 100x + 200y + 300z`; a tiny
  fixture `.h5` written via `sidecar.write_sidecar` over a 4×4×4 grid
  containing the cube; the CLI must report L² ≤ 1e-6 m/s in every
  band (linear field, exactly captured at vertices, exactly captured
  inside each tet).
- [ ] **T-1-2 (constant fixture):** synthetic constant `Vs ≡ 2500`;
  reported L² ≡ 0 (within float epsilon).
- [ ] **T-1-3 (P2 fixture):** same linear field, projected onto a P2
  GF on the same 8-tet cube; reported L² still ≤ 1e-6 (P2 captures
  P1 fields exactly). This is the test that confirms the MFEM P2
  edge ordering in Requirement 7.
- [ ] Existing unit-test suite (`make test` in `miniapps/seas/`)
  unchanged: zero C++ files modified in this phase.

### Dependencies
- Depends on: nothing.
- Required by: Phase 2 (G-1 acceptance check), Phase 3 (G-2
  acceptance check).

---

## Phase 2: G-1 — Sidecar-driven mesh size field

### Goal
The next gmsh remesh produces tets that are smaller in regions where
`‖∇Vs‖` is large (basin/basement edges, which coincide with the
streak-prone fault corridor) **without** exceeding 2× the current tet
count. After remesh + reproject, `fault_zone_metric.py` reports
`[0–500 m]` band L² ≤ 80 m/s on the 500 m configuration.

### Files to Create
- `safs/project_7.0_alternative/code_preprocess/data_projection/build_size_field.py`
  — Python tool. Reads the sidecar, computes `‖∇Vs‖`, writes a gmsh
  PostView `.pos` file. ≈ 250 LOC.

### Files to Modify
- `safs/project_7.0_alternative/code_meshing/safs_fault_box_nwcut.geo`
  — add an *optional* `Field[PostView]` reading the new `.pos` file,
  combined with the existing `Field[Threshold]` (distance to fault) via
  `Field[Min]`. Gated on `USE_SIZE_FIELD = 0|1` so existing
  reproducibility runs are byte-identical when the flag is off.
  ≈ +35 lines.
- `safs/project_7.0_alternative/code_meshing/run_nwcut_meshing.py`
  — add `--size-field-pos PATH` CLI flag forwarded to gmsh as
  `-setstring size_field_pos PATH -setnumber USE_SIZE_FIELD 1`.
  Add a `--gen-size-field` convenience that calls
  `build_size_field.main(args)` first, drops the `.pos` next to the
  sidecar, and feeds it back. ≈ +40 lines.

### Detailed Requirements

1. **`build_size_field.py` CLI**:
   ```
   python build_size_field.py
       --sidecar       PATH                    # required
       --field         Vs                      # default Vs
       --out-pos       PATH                    # required
       --lc-near       1500                    # default (= geo default)
       --lc-far        10000                   # default (= geo default)
       --lc-min        500                     # absolute floor
       --alpha         8.0                     # gradient sensitivity
       --smooth-sigma  1.0                     # Gaussian smoothing of grad
                                                #   in voxel units
       --voxel-stride  2                       # downsample factor for the
                                                #   PostView grid (memory)
       --quiet
   ```

2. **Algorithm** (single pass, all numpy):

   Inputs (read once from `--sidecar`): `gx, gy, gz` (1-D monotone
   axes, m, UTM 11N — already enforced by `sidecar.py:24-33`) and
   `F` of shape `(Nx, Ny, Nz)`.

   - Voxel-spacing arrays: `dx, dy, dz` from
     `np.diff(g{x,y,z}).mean()` (the writer enforces strictly
     monotone but not uniform; report `min/max(dxi)/mean(dxi)` so the
     user can detect a non-uniform sidecar).
   - Gradient magnitude:
     ```python
     gFx, gFy, gFz = np.gradient(F, dx, dy, dz, edge_order=2)
     gmag = np.sqrt(gFx**2 + gFy**2 + gFz**2)              # 1/s
     if smooth_sigma > 0:
         from scipy.ndimage import gaussian_filter
         gmag = gaussian_filter(gmag, sigma=smooth_sigma)
     gmag /= max(gmag.max(), 1.0e-30)                      # → [0, 1]
     ```
   - Target tet size at every voxel:
     ```
     LC = lc_far / (1 + alpha * gmag)   # alpha controls tightness
     LC = clip(LC, lc_min, lc_far)
     ```
     This is the **exact** formula promised in v3 §3.1
     (clamp + α-knob + LC_MIN floor). Default α = 8 brings basin-edge
     voxels to LC ≈ lc_far / 9 ≈ 1.1 km (lc_far = 10 km), well
     between lc_near (1.5 km, fault corridor) and lc_far (10 km, far
     field), which is the regime where a 2× tet-count budget can
     absorb the refinement without overshoot.
   - Downsample by `--voxel-stride` so the PostView is small and
     gmsh evaluates it cheaply at every node candidate.
   - Write a gmsh PostView `.pos` (legacy ASCII format — gmsh accepts
     this directly via `Merge`) using `gmsh.view.addModelData` if
     `gmsh-sdk` Python bindings are available, else fall back to a
     hand-rolled writer. The `.pos` payload is a "ScalarPoint" view
     (`SP`) of `(x, y, z, LC)` records over the downsampled grid.
     Hand-rolled writer is ≈ 30 lines and avoids the gmsh Python
     dependency at preprocess time. **Preferred:** use the
     hand-rolled writer; gmsh's `Merge` reads ScalarPoint .pos files
     identically to `addModelData`.

3. **`safs_fault_box_nwcut.geo` modification** (the only `.geo` change):
   ```gmsh
   // After the existing Field[1] (Distance) + Field[2] (Threshold)
   // and *before* `Background Field = ...`:
   If (!Exists(USE_SIZE_FIELD)) USE_SIZE_FIELD = 0; EndIf

   If (USE_SIZE_FIELD == 1)
       If (!Exists(size_field_pos))
           Error("USE_SIZE_FIELD=1 requires -setstring size_field_pos PATH");
       EndIf
       Merge Str(size_field_pos);                   // brings in PostView 0
       Field[3] = PostView;
       Field[3].ViewIndex = 0;

       Field[4] = Min;
       Field[4].FieldsList = {2, 3};
       Background Field = 4;
   Else
       Background Field = 2;                        // unchanged default
   EndIf
   ```
   The `Else` arm reproduces the current behaviour byte-for-byte; the
   only observable difference when `USE_SIZE_FIELD == 0` is that
   gmsh now sees a guarded conditional instead of an unconditional
   `Background Field = 2`.

4. **`run_nwcut_meshing.py` modification:**
   - Add `--size-field-pos PATH` (default `None`, mutually
     compatible with `--sidecar` — passing both drives the size
     field from the same sidecar that scopes the bbox).
   - Add `--gen-size-field` boolean. When set, before each gmsh call
     synthesise `<base>.size_field.pos` next to the per-resolution
     output via:
     ```python
     from data_projection.build_size_field import build_pos
     build_pos(sidecar=args.sidecar, field="Vs",
               out_pos=base.with_suffix(".size_field.pos"),
               lc_near=args.lc_near, lc_far=args.lc_far,
               lc_min=args.lc_min, alpha=args.alpha,
               smooth_sigma=args.smooth_sigma,
               voxel_stride=args.voxel_stride)
     ```
   - When either `--size-field-pos` or `--gen-size-field` is active,
     append `-setnumber USE_SIZE_FIELD 1 -setstring size_field_pos
     PATH` to `run_gmsh`'s gmsh command list.
   - The `LC_MIN`, `α`, `smooth_sigma`, `voxel_stride` knobs surface
     as `--lc-min`, `--alpha`, `--smooth-sigma`, `--voxel-stride`
     CLI flags so a user can tune without editing scripts.
   - **Mesh-budget guard.** After `run_msh_to_vtu` reports the new
     bulk-cell count, compare against the baseline cell count for
     this resolution (read from a small JSON state file
     `code_meshing/.baseline_cell_counts.json`, populated by hand the
     first time from the existing 500/1000/2000 m vtu files).
     If `new / baseline > 2.0`, print `WARNING: G-1 remesh exceeds
     2× cell-count budget (got X, baseline Y); consider lowering
     --alpha or raising --lc-min` but **do not fail** — the user
     decides whether to ship the larger mesh.

5. **No edits to** `safs_fault_box.geo`, `safs_fault_box_buried.geo`,
   `safs_fault_box_freesurface_*.geo`. G-1 targets the NW-cut variant
   that the streak diagnosis was made on.

6. **Cell-count budget arithmetic** (sanity check the user's 2×
   ceiling answers the L² target):
   - Per-tet L² scales ~ `‖∇F‖ × h` for a smooth source
     (`‖∇F‖_basin_edge ≈ 0.5 km/s / 1 km`); halving `h` in basin-edge
     bands drops the band L² ≈ 2×, which on the 500 m mesh is
     `148 → 74 m/s`, just under the 80 m/s target. Tet count grows
     ≈ `(2)^3 = 8×` *in the refined band only*; if the band volume
     is ≤ 15 % of the domain (rough estimate from the `[0–500 m]`
     band's tet share in `plot_comparison_with_dg0.py`'s output),
     total cell count grows ≈ `1 + 0.15 × 7 = 2.05×`. Tight against
     the 2× budget; default `α = 8` and `LC_MIN = 500` are the
     starting knobs. If the count overshoots, the first lever is
     `--alpha 6`, then `--smooth-sigma 1.5`, then `--lc-min 600`.

### Interfaces

- New CLI: `build_size_field.py` (standalone) +
  `run_nwcut_meshing.py --gen-size-field --sidecar PATH`.
- New importable Python: `build_size_field.build_pos(sidecar=...,
  field='Vs', out_pos=..., lc_near=..., lc_far=..., lc_min=...,
  alpha=..., smooth_sigma=..., voxel_stride=...) -> None`. Same
  argument names as the CLI flags.
- New gmsh template variables: `USE_SIZE_FIELD` (number, default 0),
  `size_field_pos` (string, required when `USE_SIZE_FIELD == 1`).

### Edge Cases to Handle
- Sidecar has uniform `‖∇F‖ ≡ 0` (constant field): `LC ≡ lc_far`
  everywhere. The PostView is still written (1 voxel × N entries
  with constant LC); gmsh's `Field[Min]` then reduces to the
  existing distance threshold field. Acceptance: L² for the constant
  fixture stays at machine zero.
- Voxel stride yields fewer than 4 PostView records along any axis:
  abort `build_size_field` with a clear "stride X drops axis '%s'
  to %d voxels — gmsh's PostView interpolant requires ≥ 4" message.
  Mitigation: lower stride.
- Mesh-budget guard cannot find the baseline JSON: prompt the user
  to run once without `--gen-size-field` to populate it; fall back
  to "no baseline known, cannot enforce 2× cap".
- Sidecar `‖∇F‖` is concentrated in <1 % of voxels: `LC` saturates
  at `lc_min` everywhere those voxels touch, which can produce a
  narrow strip of tiny tets. The 1-voxel Gaussian smoothing default
  (`--smooth-sigma 1.0`) already mitigates this; the `--lc-min`
  floor is the second line of defence.
- `--gen-size-field` invoked without `--sidecar`: error out with a
  clear "size-field requires --sidecar" message (the bbox-only
  `--velocity-bbox` path has no Vs data to differentiate).

### Acceptance Criteria
- [ ] **R-2-G-T1 (no-op equivalence):** running
  `python run_nwcut_meshing.py --res 500` *without* `--gen-size-field`
  produces a mesh whose tet count and per-element edge-length
  histogram match the existing `safs_fault_box_nwcut_500m.msh`
  exactly. (Confirms the `.geo` `Else`-arm is byte-identical.)
- [ ] **R-2-G-T2 (acceptance gate):** after
  ```
  python run_nwcut_meshing.py --res 500 \
      --sidecar ../data_projected/velocity_safs.h5 \
      --gen-size-field
  ```
  followed by the existing `seas_project_velocity_to_mesh` driver
  (current default `--order 1` — keeps the variable change isolated
  from Phase 3), the **`fault_zone_metric.py`** CLI reports
  `[0–500 m]` band L² ≤ **80 m/s** (down from the v3-baseline 148
  m/s).
- [ ] **R-2-G-T3 (mesh budget):** the same remesh produces
  `new_cell_count ≤ 2.0 × baseline_cell_count` for the 500 m
  configuration. If overshoot, the warning fires and the run is
  marked needs-tuning (does not block CI but is loud in stdout).
- [ ] **R-2-G-T4 (bulk no-regression):** `[3000, ∞)` band L² stays
  within ± 5 m/s of the baseline (no degradation in the far field).
- [ ] **R-2-G-T5 (containment):** the new mesh's bbox is still
  contained in the sidecar bbox (the existing `--sidecar` headroom
  check in `run_nwcut_meshing.py:243-278` continues to pass).
- [ ] **R-2-G-T6 (1000 m + 2000 m):** the same flag combination on
  1000 m / 2000 m drops `[0–500 m]` L² by ≥ 30 % each (target: 207
  → ≤ 145 m/s on 1000 m; 238 → ≤ 167 m/s on 2000 m). Tighter
  acceptance for those resolutions is deferred to Phase 3.
- [ ] **R-2-G-T7 (build_size_field unit):** linear-Vs synthetic
  sidecar (constant `‖∇F‖ = c`); the resulting PostView has `LC` ≡
  `lc_far / (1 + α c / max(c))` everywhere — round-trip the .pos
  through `gmsh.merge` and assert.
- [ ] No C++ files modified; `make test` unchanged.

### Dependencies
- Depends on: Phase 1 (the metric CLI is the acceptance gate).
- Required by: Phase 3 (G-2 measures its incremental gain on the
  G-1 mesh).

---

## Phase 3: G-2 — Bump default H1 order to P2

### Goal
On the same mesh produced by Phase 2, the per-tet L² drops further by
adding within-tet curvature: `H1_FECollection(2, dim)` instead of
`H1_FECollection(1, dim)`. The change is *one default-argument bump*
in the existing driver; the rest of the projector (bbox containment,
range checks, derived `μ = ρVs²`, `λ = ρVp² − 2µ`) is order-agnostic
and untouched.

### Files to Create
- *(none)*.

### Files to Modify
- `miniapps/seas/drivers/project_velocity_to_mesh.cpp` — exactly two
  edits:
  1. Line 42: `int order = 1;` → `int order = 2;`.
  2. Line 53-54 help string: append `" (G-2 default; pass --order 1 "
     "to reproduce pre-G-2 H1-P1 output.)"`.
  Total diff: ≈ 4 lines.
  No header / signature change. The existing `pv.SetLevelsOfDetail
  (order)` and `pv.SetHighOrderOutput(order > 1)` calls (lines 114-
  115) already do the right thing for `order = 2`.

### Detailed Requirements

1. **Why P2 specifically (not P3+).** Per-tet DOFs grow 4 → 10 → 20.
   For the SAFS sidecar at 500 m mesh, P2 captures the source's
   in-tet quadratic structure (which is most of what's missing once
   G-1 has shrunk h in basin-edge tets); P3 doubles DOFs again for
   marginal gain on a `‖∇²F‖`-bounded source. Default stays P2;
   power users still get arbitrary `--order P` via the existing CLI.

2. **No edits to** `io/field_coefficient.{hpp,cpp}` — `Project /
   ProjectVelocity` already accept any `target_fes`, and the bbox
   containment + range-check logic
   (`field_coefficient.cpp:144-153, 162-169, 191-218`) is independent
   of the FEC order.

3. **No edits to** `io/data_field_3d.{hpp,cpp}` — the trilinear
   evaluator on the sidecar grid is identically called from
   `FieldCoefficient::Eval` regardless of where MFEM's `ProjectCoefficient`
   chose to evaluate it (vertices for P1; vertices + edge-mids for P2).

4. **ParaView output footprint.** P2 tets emit ~2.5× more
   `<DataArray>` entries; PNG-rendering preview pipelines downstream
   need to be aware. The existing
   `plot_comparison_with_dg0.py` slicing path (`render_h1_at_z_via_tet_slice`,
   line 114-161) currently assumes 4-vertex tets and barycentric
   linear interpolation; on a P2 vtu it will silently *under-render*
   (treat each P2 cell as if it were P1 between corner vertices,
   missing the edge-mid DOFs). This is a documentation-only issue
   for Phase 3 — **the metric CLI from Phase 1 already handles P2
   via Requirement 7** and is the authoritative measurement.
   Document the limitation in a comment at the top of
   `plot_comparison_with_dg0.py` so future readers don't trust the
   PNG over the metric.

5. **`SetHighOrderOutput(order > 1)`** at driver line 115 already
   activates ParaView's high-order rendering path; verify that the
   resulting `.pvd` opens without warnings in the user's local
   ParaView (≥ 5.13 supports tet P2 high-order natively).

### Interfaces
- No new functions, classes, or files. The driver CLI is
  back-compatible: any user who passes `--order 1` reproduces the
  pre-G-2 output exactly.

### Edge Cases to Handle
- A user with ParaView < 5.10 will see "unsupported high-order cell"
  errors. The driver's stdout already lists `order`; add a one-line
  hint when `order > 1`: "ParaView ≥ 5.10 required for high-order
  tet rendering."
- P2 increases per-rank `ParGridFunction` size by ~2.5×; verify on
  the 500 m mesh that the projection still fits in single-node
  memory before recommending Frontera for production runs. (Per
  user feedback memory `feedback_local_mpi_up_to_10.md`, the
  500 m projection runs locally with `np ≤ 10`; Phase 3 acceptance
  is locally measurable.)

### Acceptance Criteria
- [ ] **R-3-G-T1 (default landed):** invoking
  `seas_project_velocity_to_mesh --mesh ... --sidecar ... --out ...`
  with no `--order` flag produces a `.pvd` whose `<Piece>` cells are
  P2 (10-DOF) tets (verifiable by `grep` of the .vtu header).
- [ ] **R-3-G-T2 (acceptance gate, stacked):** Phase-1 metric CLI on
  the Phase-2 G-1 mesh + Phase-3 default-P2 projection reports
  `[0–500 m]` band L² ≤ **50 m/s** (≥ 30 % drop from G-1 alone, per
  v3 §3.2 prediction).
- [ ] **R-3-G-T3 (back-compat):** same driver invoked with
  `--order 1` reproduces the pre-Phase-3 output bit-identically
  (vertex values match the pre-G-2 baseline `.pvd` to machine
  precision).
- [ ] **R-3-G-T4 (no-regression on existing C++ tests):**
  `make test` in `miniapps/seas/` passes unchanged. (None of those
  tests build or run `seas_project_velocity_to_mesh`; the change is
  isolated to one driver default.)
- [ ] **R-3-G-T5 (range bounds):** the existing range check on
  derived `μ ∈ [1e7, 1e11]` Pa and `λ ∈ [0, 1e12]` Pa
  (`field_coefficient.cpp:268-277`) does not trip on the P2-
  projected output. (P2 projection can over/undershoot at edge
  nodes for fields with sharp gradients; confirm empirically. If
  it trips, the mitigation is **not** to widen bounds in C++ but
  to investigate whether G-1's size field needs more refinement
  in the offending region — this is a feedback loop into Phase 2,
  not a Phase 3 reopening.)

### Dependencies
- Depends on: Phase 2 (so we're measuring on the G-1 mesh, not the
  baseline).
- Required by: nothing in this PLAN. Downstream `WaveOperator`
  consumption of the projected `MaterialField` is the separate
  Phase-5 work in `data_projection_feature_plan_v2.md` and is
  out of scope here.

---

## Testing Strategy

| Phase | Test | Where | What it pins |
|-------|------|-------|--------------|
| 1 | T-1-1 P1 linear fixture | `code_preprocess/data_projection/test_data_projection.py` (extend) | Cubature + barycentric H1-P1 evaluation correctness. |
| 1 | T-1-2 constant fixture | same file | Sanity floor (L² ≡ 0). |
| 1 | T-1-3 P2 linear fixture | same file | MFEM P2 edge ordering + closed-form Lagrange basis. |
| 1 | T-1-4 baseline-table reproduction | same file (or a separate `test_fault_zone_metric.py`) | The CLI numbers reproduce v3 §0's 148 / 207 / 238 m/s on real data, gated on the SAFS sidecar fixture if present. |
| 2 | R-2-G-T1 no-op equivalence | manual once + JSON snapshot in CI | The `.geo` `Else`-arm is byte-identical to today. |
| 2 | R-2-G-T2 .. T6 acceptance | `fault_zone_metric.py` driven from `run_nwcut_meshing.py` post-step | The headline 80 m/s gate, plus 1000 / 2000 m and bulk no-regression. |
| 2 | R-2-G-T7 build_size_field unit | new `test_build_size_field.py` | PostView numerics + .pos round-trip. |
| 3 | R-3-G-T1 default landed | regex on a regenerated `.vtu` header in CI | Driver default actually changed. |
| 3 | R-3-G-T2 stacked gate | metric CLI | The 50 m/s post-G-2 target. |
| 3 | R-3-G-T3 back-compat | `--order 1` rerun + `diff` against snapshot | No silent semantic change for users who pin order. |
| 3 | R-3-G-T4 existing tests | `make test` | C++ regression cordon. |
| 3 | R-3-G-T5 range bounds | observation in R-3-G-T2's run logs | Catch P2 over/undershoot at nodes. |

The reference solution for every numerical check is the trilinear
evaluation of the sidecar at the cubature points (closed-form on the
sidecar grid; deterministic). No analytical closed-form is needed
beyond the unit-cube fixtures; the SAFS data is the integration test.

## Risk Assessment

- **G-1 size field shape can over-refine basin-edge strips.** The
  default `α = 8`, `LC_MIN = 500`, `smooth_sigma = 1.0` are the
  v3-recommended starting point and are tunable knobs, not invariants.
  Detection: the mesh-budget guard (R-2-G-T3) triggers a WARNING when
  cell count > 2× baseline. Mitigation order:
  `α 8 → 6 → 4`, then `smooth_sigma 1.0 → 1.5`, then
  `LC_MIN 500 → 600 → 750`. Document the chosen knobs in the .pvd's
  `mesh_tag` attribute (already supported by `sidecar.py:62`).
- **gmsh `Field[Min]` interaction with the existing
  `Field[Threshold]` may not behave as a smooth blend.** `Field[Min]`
  is a hard `min(...)`; in transition regions between fault-distance
  and gradient refinement the size field can have C⁰ kinks that gmsh
  resolves via small tets. Detection: edge-length histogram tail
  (already reported by `msh_to_vtu.py` per
  `run_nwcut_meshing.py:447`). Mitigation: replace `Field[Min]` with
  `Field[Restrict]` + `Field[Mean]`-style smoothing if the kink shows
  up in practice. Out of scope of the first G-1 landing — only invoked
  if R-2-G-T6 reveals a quality-stat regression.
- **P2 default may break a downstream consumer that hard-codes
  P1.** No such consumer exists in `miniapps/seas/` today (the only
  caller of `ProjectVelocity` is the standalone driver). External
  scripts under `safs/.../analysis_*` may do their own .vtu parsing.
  Detection: R-3-G-T3 (the back-compat snapshot) verifies the
  `--order 1` escape hatch works; document the default change in
  the project's release notes.
- **MFEM P2 tet edge ordering.** The closed-form Lagrange basis in
  Phase-1 Requirement 7 assumes the MFEM ordering `(0,1), (0,2),
  (0,3), (1,2), (1,3), (2,3)`. If MFEM's local edge ordering is
  different (e.g., due to a `permutation` map), the P2 metric will
  silently mis-evaluate. T-1-3 catches this on a unit fixture
  before any production run.
- **Sidecar `‖∇F‖` is computed in UTM-11N metres, not lon/lat.** The
  v1 schema enforces UTM (`sidecar.py:19`); `np.gradient` over
  `gx, gy, gz` (passed in metres) yields units of `1 / m` × field
  unit. This is the right thing for `α` to be dimensionless when
  `F` is normalised by `max(F)`. Document in `build_size_field.py`'s
  module docstring; the unit test `R-2-G-T7` pins the numerics.
- **Tricky existing code that interacts with this plan.** None of
  the C2-invariant files are touched. The
  `field_coefficient.cpp::ProjectVelocity` path (lines 226-291)
  derives `μ`, `λ` per-DOF from the projected `Vp`, `Vs`, `ρ`; this
  is mathematically identical for P1 vs P2 vs Pp because the algebra
  is applied at every DOF index — no order-dependent assumption is
  baked in. (Verified by reading lines 252-263.)

## Sequencing & Stop-Gates

| Order | Phase | Effort | Stop-gate |
|-------|-------|--------|-----------|
| 1 | Phase 1 metric CLI | 0.5 day | T-1-1..T-1-4 pass; baseline 148/207/238 reproduced. |
| 2 | Phase 2 G-1 size field + remesh | 1.5 day | R-2-G-T2 hits ≤ 80 m/s on 500 m within the 2× budget (R-2-G-T3). If yes → optional ship; either way → Phase 3. |
| 3 | Phase 3 G-2 P2 default | 0.5 day | R-3-G-T2 hits ≤ 50 m/s stacked. R-3-G-T3 back-compat passes. |

If Phase 2 cannot hit ≤ 80 m/s within the 2× cap after the documented
knob sequence (`α 6 → 4`, `smooth_sigma 1.5`, `LC_MIN 600 → 750`),
do **not** silently widen the cap — surface the data and re-engage
the user. The cap is a constraint, not a soft target.

---

## Files NOT to Touch (reaffirmed)

- C2-invariant set: `bp5/`, `bp1/`, `bp2/`, `domain/`,
  `friction/dieterich_ruina.hpp`, `solver/`, `fault/fault_basis.hpp`.
- `dynamic/heterogeneous_material.hpp`. Do NOT add
  `Mode::DG_L2_GridFunction` (v2 retraction stands).
- `io/field_coefficient.{hpp,cpp}`. Do NOT add `ProjectDG0` /
  `ProjectVelocityDG0` (v2 retraction stands).
- `io/data_field_3d.{hpp,cpp}`. The trilinear evaluator and bbox
  helpers are stable infrastructure.
- `code_preprocess/data_projection/sidecar.py`. The v1 schema is not
  being bumped; no new attributes, no new datasets.
- All driver source files except
  `drivers/project_velocity_to_mesh.cpp` (and that file gets exactly
  the two-line default change in Phase 3).
