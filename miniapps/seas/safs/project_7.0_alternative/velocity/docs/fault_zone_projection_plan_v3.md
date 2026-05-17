# Fault-Zone Projection-Quality Plan (v3 — retraction & redirect)

**Date:** 2026-05-09
**Status:** Supersedes the recommendation in v2 §2 (F-2) on the basis
of new per-tet integrated L² measurements.  v1 + v2 stay as design
records; the recommended next action changes.
**Triggering data:** new comparison panels at
`data_projected/comparison_sidecar_vs_mesh_*_z1km_v2.png` (sidecar
| H1-P1 | DG-0 centroid | DG-0 vol-avg | diff | per-tet L²).

---

## 0. Summary of the retraction

v2 §2 proposed **F-2 (DG-L²(0))** as the primary fix on the strength
of "the streaks disappear."  Implementation was paused before any
source edits.  A proper per-tet L² measurement now shows F-2 is
strictly worse than H1-P1 on this dataset:

| Mesh   | per-tet L² of (proj − source) over the slab,  in m/s |
|--------|------------------------------------------------------|
|        | H1-P1   |   DG-0 centroid   |   DG-0 vol-avg            |
| 500 m  | **148** |       250         |       227                 |
| 1000 m | **207** |       352         |       318                 |
| 2000 m | **238** |       395         |       353                 |

(Units: m/s; integrated over the tet via 4-point Gauss cubature, then
averaged over all tets straddling z = −1 km ± 200 m.  The same
ranking holds in every distance-to-fault band.)

**Why H1-P1 wins (and what we missed in v2):**

H1-P1 stores 1 GF DOF per vertex but every interior tet sees 4
distinct vertex values.  In-tet linear interpolation between those 4
values produces an O(h²)-accurate approximation when the source has
local linear structure — which is *most* of the SAFS sidecar volume,
even near the fault.  DG-L²(0) collapses 4 numbers per tet down to 1,
giving up the within-tet linear capture; the centroid sample (or its
volume average) is at best O(h) accurate.  v2 §0 incorrectly equated
"H1 vertex projection is bit-correct AT vertices" with "H1 has the
same per-element fidelity as DG-0"; those are different statements.

**The streaks the user originally flagged are not a fixable artefact.**
They are an honest depiction of in-tet linear interpolation on a mesh
whose tets are larger than the source's gradient length scale.
Replacing in-tet linear with in-tet constant *adds* error, even if
the visual reads more cleanly because the cell boundaries are pixel-
aligned.

---

## 1. What the v2 plan got right (kept)

- The diagnosis sequence: the H1 *vertex projection* is bit-correct;
  the streaks come from in-tet linear interpolation across the
  source's sharp basin/basement gradients.  v2 §0 still applies as
  the descriptive truth of the visualisation.
- The infrastructure inventory: every tool listed in v1 §2 still
  exists and is reusable.
- The C2-invariant boundary list (no edits to `bp5/`, `bp1/`, `bp2/`,
  `domain/`, `friction/dieterich_ruina.hpp`, `solver/`,
  `fault/fault_basis.hpp`).

## 2. What v2 got wrong (retracted)

- **v2 §2 (F-2 DG-L²(0))**: not a fix, regression on this dataset.
  Skip the implementation in §2.1–2.4 entirely.  Don't add
  `ProjectDG0`, don't add `--basis dg0` to
  `drivers/project_velocity_to_mesh.cpp`, don't add
  `Mode::DG_L2_GridFunction` *for the projection-preview path*.
  (See §5 below for one place where DG-0 is still the right
  answer — the wave-op-side material lookup, which is a different
  use case.)
- **v2 §3 (F-3 per-side fault projection)**: the rationale was that
  it pairs naturally with DG-0.  With DG-0 retracted, F-3 becomes a
  no-op on the current sidecar.  Defer F-3 entirely.

---

## 3. Recommended path forward (v3 priorities)

The dominant error is *mesh-resolution-limited representation* of a
*smooth-but-sharply-graded* source.  The only operations that move
the per-tet L² needle are:

| Phase | Lever | Mechanism |
|-------|-------|-----------|
| **G-1** | Sidecar-driven mesh refinement (= v1 Lever T5 / v1 phase F-4) | Make tets smaller where ‖∇Vs‖ is large.  Unconditional gain. |
| **G-2** | Higher-order H1 (H1-P2 or P3) | More within-tet DOFs; O(h^{p+1}) for smooth source. |
| **G-3** | (defer) Higher-order DG (DG-1 or DG-2) | Same as G-2 in fidelity for smooth source; allows face jumps if a future sidecar adds them. |

### 3.1 Phase G-1 — sidecar-driven mesh refinement

**Goal:** drive `LC_NEAR` and the `Field[Distance]` extent in
`code_meshing/safs_fault_box_nwcut.geo` from `‖∇Vs‖` of the sidecar
itself, not just from fault distance.  Refine where the source has
sharp gradients (which happen to coincide with basin/basement
boundaries — i.e., exactly the locations where the streaks appear).

**Files (offline-only; no C++ edits):**
- New: `code_preprocess/data_projection/build_size_field.py` —
  reads `velocity_safs.h5`, computes `‖∇Vs‖` at every voxel, writes
  a gmsh `.pos` file (Field[PostView]) with target tet size
  `LC = clamp(LC_FAR / (1 + α ‖∇Vs‖ / Vs_med), LC_MIN, LC_FAR)`.
  α is a CLI knob; default `α = 8`.
- Modified: `code_meshing/safs_fault_box_nwcut.geo` — add
  `Field[PostView]` reading the `.pos` and combine with the
  existing distance-to-fault field via `Field[Min]`.  Wrapped
  behind a gmsh variable `USE_SIZE_FIELD = 0|1`; default off.
- Modified: `code_meshing/run_nwcut_meshing.py` — optional
  `--size-field-pos PATH` flag forwarded to gmsh.

**Acceptance (on the F-1 metric, `[0–500 m]` band):**
- Per-tet L² ≤ **80 m/s** (≈ 50 % reduction from H1-P1 on the same
  baseline mesh) at a total tet count ≤ 2× the current 500 m mesh.
  (Arithmetic: per-tet L² scales as ~h on a sharp-gradient slab; a
  2× refinement should deliver a 2× drop, so 148 → ~80 m/s is the
  expected outcome.)

**Cost:** ~250 LOC of Python preprocessing + a mesher rerun.  No
solver edits.  Mesh size grows by some factor — measured per run.

**Risk:** if the gradient field is concentrated in too few voxels,
the size field might force unrealistically small tets in narrow
strips.  Mitigation: lower-clamp `LC_MIN` and 3-D Gaussian smoothing
of `‖∇Vs‖` before exporting the `.pos`.

### 3.2 Phase G-2 — H1-P2 projection

**Goal:** from `H1_FECollection(1, dim)` to `H1_FECollection(2, dim)`
in `drivers/project_velocity_to_mesh.cpp` (the existing `--order P`
flag already supports this — just bump the default).  Extra cost is
the larger `ParGridFunction` size and the L²-projection time at init.

**Files (the only solver-side change in this plan):**
- Modified: `drivers/project_velocity_to_mesh.cpp` — change default
  `order = 1` → `order = 2`.  Existing `--order P` CLI continues to
  work.
- (No changes to `field_coefficient.{hpp,cpp}` —
  `ProjectCoefficient` already supports any FES order.)

**Acceptance:**
- Per-tet L² should drop ≈ √(N_dof_per_elem ratio).  For tets,
  P1 → P2 is 4 → 10 DOFs (~2.5×), and the L² of a smooth field
  drops by O(h²)→O(h³) ≈ h → h² near gradients, so ≥ 30 % drop on
  the same mesh, additive with G-1.

**Cost:** ~5 LOC.  Increases preview-output size proportionally.

**Risk:** H1-P2 GFs interleave more DOFs per element; ParaView 5.13+
handles this cleanly via `SetHighOrderOutput(true)` (the driver
already does this for `order > 1`).

### 3.3 Phase G-3 (deferred) — DG path

DG basis becomes the right answer the day a future sidecar contains
a real fault-perpendicular Vs jump (e.g., basin / basement contrast
on opposite sides of the SAF).  Until then, H1-Pp matches DG-Pp on
fidelity and is cheaper to render.  Don't pre-emptively add the DG
projection path; revisit when there's a discontinuous source.

---

## 4. Recommended sequence

| Order | Phase | Effort | Decision gate |
|-------|-------|--------|---------------|
| 1 | Repurpose `fault_zone_metric.py` (F-1 from v1, never built) using the v3 per-tet integrated L² formula already implemented in `plot_comparison_with_dg0.py` | 0.5 day | Lands the metric we drive down. |
| 2 | G-1: sidecar-driven size field + remesh @ same total budget | 1.5 days | If `[0–500 m]` band L² ≤ 80 m/s → ship. |
| 3 | (assess) | — | Otherwise → G-2. |
| 4 | G-2: bump `--order` default to 2; rerun preview | 0.5 day | If L² ≤ 50 m/s → ship. |
| 5 | (assess) | — | Otherwise → revisit gradient-budget assumptions; refine size-field α. |

---

## 5. The DG-0 question — where it's still useful

DG-0 / `Mode::DG_L2_GridFunction` was never just about projection
quality; it's also the natural data layout for the wave operator's
ADER/RK4 integrand on heterogeneous material.  When (separately)
the WaveOperator gets wired to consume a `MaterialField::Mode::
GridFunction`, that material can either be H1-P1 (continuous, what
the current preview produces) or DG-0/DG-1 (per-element constants /
per-element linear).  That choice trades:

- H1-P1 material: best representation fidelity (this plan v3) but
  smooths face-located source jumps.
- DG-0 material: cheapest per-element lookup (one number, no `vdofs`
  resolution); honest about face discontinuities but worse for smooth
  data.
- DG-1 material: best of both; per-element linear; allows face
  jumps; small overhead vs DG-0.

**Recommendation for the wave-op coupling phase
(`data_projection_feature_plan_v2.md` Phase 5):** consume H1-P_p as
projected by *this* plan v3, and **only** swap to DG-1 when the
sidecar acquires a real face-located discontinuity.  Don't introduce
DG-0 for the consumer either — DG-0's per-element representation
collapse hurts the wave op the same way it hurts the preview.

---

## 6. Files NOT to touch (reaffirmed)

- C2-invariant set: `bp5/`, `bp1/`, `bp2/`, `domain/`,
  `friction/dieterich_ruina.hpp`, `solver/`, `fault/fault_basis.hpp`.
- `dynamic/heterogeneous_material.hpp`: leave the existing `Mode::
  Constant` and `Mode::GridFunction` paths alone.  Do **not** add
  `Mode::DG_L2_GridFunction` based on v2 §2.2 — that change is
  retracted.
- `io/field_coefficient.{hpp,cpp}`: leave the existing
  `Project / ProjectVelocity` API alone.  Do **not** add
  `ProjectDG0 / ProjectVelocityDG0` based on v2 §2.1 — retracted.

The only file change this plan recommends inside the C++ tree is a
**default-argument bump** for `--order` in
`drivers/project_velocity_to_mesh.cpp` (G-2), and that only after G-1
is measured.

---

## 7. Open questions (refreshed)

1. Acceptance target on `[0–500 m]` band per-tet L²: 80 m/s
   (proposed) or tighter?
2. Mesh-budget ceiling: 2× current tet count, or higher?
3. Should we still keep DG-0 plumbing in the wave op for a future
   discontinuous sidecar, even though it's not the right fix today?

Send answers and I'll write the G-1 build_size_field.py spec in the
same file-level detail as v2 §2.1.
