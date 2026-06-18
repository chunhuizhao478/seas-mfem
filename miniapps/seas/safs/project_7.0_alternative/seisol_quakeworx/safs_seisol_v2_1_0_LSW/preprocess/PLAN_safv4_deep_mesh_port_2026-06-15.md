> ## LSW twin adaptation note (2026-06-15)
> **This file is the RSSRW deep-mesh port plan, reproduced verbatim in the LSW
> case folder (`safs_seisol_v2_1_0_LSW/preprocess/`) because the LSW twin reuses
> the RSSRW deep-mesh artifacts *byte-for-byte*.** The mesh is friction-agnostic,
> so the entire mesh-processing methodology and execution log below (Phases 0–4:
> preflight → retag → pumgen → PUML verify) apply unchanged. The scripts
> (`preflight_deep_mesh.py`, `retag_safv4_deep_to_seissol.py`,
> `verify_puml_deep.py`, `build_and_run_pumgen_linux.sh`), the retagged
> `safs_mesh_deep.msh`, and the produced `safs_mesh_deep.puml.h5` are identical to
> the RSSRW ones — the LSW twin **copied the verified PUML rather than re-running
> pumgen** (PUML verify re-run here: PASS — tets 1,180,159 / nodes 279,448 / BC1
> 175,033 / BC5 122,593 / BC3 264,666=2×).
>
> **What differs in the LSW twin (friction only):**
> - `parameters.par`: `FL = 16` (linear slip-weakening), **not** `FL = 103`.
>   The RS constants block (`RS_f0/RS_b/RS_muW/...`) does not exist in the LSW
>   `parameters.par`; the LSW `OutputMask`/`EndTime` comments are LSW-specific.
> - `safs_fault.yaml`: the LSW model (`[mu_d,d_c,cohesion]` ConstantMap +
>   `[mu_s]` LuaMap deep barrier + `Tnuc_s` Gaussian). Reused verbatim from
>   `safs_seisol_v2_0_0_LSW` — coordinate-based, so it applies to every fault
>   facet automatically (same reuse logic this plan uses for the RSSRW yamls).
> - Nucleation amplitude `Tnuc_s = 10 MPa` (LSW), **not** 24 MPa (RSSW).
> - **Deep barrier vs the deep mesh (verified):** the LSW barrier locks
>   `mu_s = 1e6` for `z ∈ (−20000, −15000)`. The deep mesh fault z-extent is
>   −19,329 … +2,210 m (deepest facet −19.3 km; **0 facets below −20 km**;
>   17.3 % below −15 km), so the barrier correctly locks the deepest ~4.3 km of
>   fault and no facet falls below the band's lower edge. **No barrier edit is
>   needed.**
> - The Phase-5 run gates differ (no RS-parameter-source check; nucleation
>   10 MPa). See `README_deep_run_and_gates.md` (the LSW version in this folder)
>   for the LSW-specific gates; the RSSRW-only friction calc
>   (`CALC_rssw_friction_lengthscales`) is intentionally **not** reproduced here —
>   the LSW friction history lives inline in `safs_fault.yaml`.
>
> Everything below this box is the original RSSRW plan, unmodified.

---

# Implementation Plan: Run the RSSRW SeisSol case on the new `safv4_deep_500m` multi-strand mesh

Date: 2026-06-15
Author: planning pass (code-explore + code-plan)
Source case: `safs_seisol_v2_1_0_RSSRW/`
Target mesh: `project_7.0_preferred/meshing/results/safv4_deep_500m.msh`

## Overview

Run the exact `safs_seisol_v2_1_0_RSSRW` physics (FL=103 rate-and-state + strong
velocity weakening, CVM material via ASAGI, Gaussian over-stress nucleation) on the
new 3-strand `safv4_deep_500m.msh` instead of the original single-fault
`safs_mesh.puml.h5`. All three fault strands — San Andreas (tag 101), Banning
(tag 102), Garnet Hill (tag 103) — are merged into ONE dynamic-rupture physical group
(SeisSol boundary code 3). The friction / stress / nucleation yamls and the CVM `.nc`
are reused **verbatim** because they are coordinate-based (easi LuaMap / ASAGI), so
they apply to every fault facet automatically. The work is therefore almost entirely
mesh-processing + validation, not physics editing.

## Decision (confirmed with user, 2026-06-15)

> "All 3 rupture (shared law)" + "the fault should be treated as a whole, single
> physical group."

⇒ Retag `101, 102, 103 → 103` (one SeisSol dynamic-rupture group). There is no
per-strand distinction anywhere in the run (SeisSol's BC `3` is a single DR group;
the easi friction is purely a function of position).

## Execution log — Phases 0 & 1 (2026-06-15, PASS)

All outputs written INSIDE `safs_seisol_v2_1_0_RSSRW/` (per user constraint
"changes within this folder only"); the Phase-1 mesh output therefore lands here,
NOT in the Phase-3 sibling folder.

Created: `preflight_deep_mesh.py`, `retag_safv4_deep_to_seissol.py`,
`safs_mesh_deep.msh` (88 MB, retagged).

Phase 0 results (`python3 preflight_deep_mesh.py`, exit 0 = PASS):
- Tag inventory verified: tris `{101:100721, 102:2680, 103:28932, 201:175033,
  202:86583, 203:36010}`, tets `{1:1180159}`. bbox as in Background.
- CVM coverage: **39.2%** of nodes outside the grid xy-footprint, **28.1%** above
  z=0 (topography), 0% below z=−45000 ⇒ ASAGI clamping is real (Phase 4 confirmed).
- Hypocenter-on-fault GATE PASS: exact min point-to-triangle distance **12.6 m**
  (≪ 600 m gate). **FINDING: the nearest facet is on strand 103 = Garnet Hill, NOT
  101 = San Andreas.** Harmless because all strands are one DR group and the 6 km
  Gaussian spans nearby San-Andreas facets too — but nucleation seeds in the Garnet
  Hill region, not on the main SAF. Worth confirming the intended seed region.
- Fault orientation: median `|n_z| ≈ 0.44` ⇒ median **dip ≈ 64°** (moderately
  dipping, San Gorgonio Pass geometry), ~12% of facets dip <46°. Geologically
  realistic, but **elevates the Phase-5 per-strand strike/dip (`XRef`) risk** — one
  reference vector across three dipping strands is more likely to mirror a strand.

Phase 1 results (retag + re-parse verification, PASS):
- Counts: `103 (DR)=132,333`, `101 (free surface)=175,033`, `105 (absorbing)=122,593`.
- `$PhysicalNames` block inserted (input had none).
- Re-parse of `safs_mesh_deep.msh`: tris `{101,103,105}`, tets `{1}`, 279,448 nodes;
  node coords and triangle connectivity **byte-identical to source** (geometry
  preserved); bbox unchanged.

Folder reorg (user, 2026-06-15): final-submission files (`parameters.par`, the 3
yamls, `safs_material_cvm.nc`, `*.puml.h5`) kept at the RSSRW top level; all tooling,
the `.msh` input, docs, and `archive/` moved to `preprocess/`. Script path-defaults
fixed for the deeper location. This plan now lives in `preprocess/`.

Phases 2-4 (2026-06-15, PASS; in-place per "this folder only"):
- **Phase 2** PUML verify (`verify_puml_deep.py`, exit 0): tets 1,180,159 / nodes
  279,448; BC1 free surface 175,033; BC5 absorbing 122,593; BC3 DR **264,666 = 2x**
  the 132,333 fault triangles (this pumgen tags both adjacent tet faces - accepted);
  codes subset of {0,1,3,5}.
- **Phase 3** (deviation: IN PLACE, not a sibling folder): `parameters.par`
  `MeshFile = 'safs_mesh_deep.puml.h5'` (+ refreshed the stale retag comment). Yaml
  `!Include`/`file:` refs and the PUML all resolve at the top level. Single-fault
  `safs_mesh.puml.h5` left in place as the baseline (not in this run's upload set).
- **Phase 4** material/strike-dip gates rewritten in
  `preprocess/README_deep_run_and_gates.md` (ASAGI clamping now EXPECTED - 39.2%
  out-of-xy, 28.1% above z=0; near-fault must be un-clamped, mu~23.4 GPa spot-check;
  `.nc` not regenerable; strike/dip Gate D for the 3 dipping strands; nucleation
  seeds the Garnet Hill region).
- **Phase 5** (cluster/QuakeWorx run) prepared, NOT executed - user runs it.

## Background facts established during exploration

### Old vs new mesh tag schemes (THE central difference)

| | Old single-fault mesh | New `safv4_deep_500m.msh` (verified by parsing) |
|---|---|---|
| fault(s) | `101` (one fault) | `101` San Andreas (100,721 tri), `102` Banning (2,680), `103` Garnet Hill (28,932) |
| free surface | `102` top | `201` top / DEM topography (175,033 tri) |
| absorbing | `103` bottom, `104` sides | `202` bottom (86,583 tri), `203` sides (36,010 tri) |
| rock volume | `1` | `1` (1,180,159 tets) |
| `$PhysicalNames` block | present | **ABSENT** (only `$MeshFormat`/`$Nodes`/`$Elements`) |

Source for the new scheme: `project_7.0_preferred/document/RUNLOG_safv4_remesh_phases0-5_2026-06-12.md`
line 264 ("volume 1, faults 101/102/103 alphabetical, top 201, bottom 202, sides 203")
and the `meshing/results/gate_conformity_deep_*.json` / `vtu_safv4_deep/final_fault_*` files.

⇒ The existing `retag_msh_to_seissol.py` (`TAG_MAP = {101:103, 102:101, 103:105, 104:105}`)
would **catastrophically mistag** the new mesh: Banning→free-surface, Garnet→absorbing,
and 201/202/203 left untouched (pumgen would read them as BC 101/102/103 = garbage).
A new mapping is mandatory.

### Required SeisSol retag mapping (new mesh)

```
101 San Andreas  ┐
102 Banning      ├─→ 103   (SeisSol gmsh tag; pumgen BC 3 = dynamic rupture)   [single group]
103 Garnet Hill  ┘
201 top / DEM    ──→ 101   (pumgen BC 1 = free surface)
202 bottom       ┐
203 sides        ┴─→ 105   (pumgen BC 5 = absorbing)
  1 rock (vol)   ──→ 1     (material/volume group, unchanged)
```

Expected post-pumgen boundary-code face counts:
`BC 3 (DR) ≈ 132,333`, `BC 1 (free surface) ≈ 175,033`, `BC 5 (absorbing) ≈ 122,593`.

NOTE on collision-safety: the map sends both `101→103` and `103→103`. This is safe ONLY
if each element's ORIGINAL tag is read once and written once (no sequential re-mapping).
The original script already works this way; the new script MUST preserve that.

### Geometry vs the reused inputs

- New mesh bbox: `x 16,309 .. 711,231`, `y 3,532,686 .. 4,085,886`, `z −39,329 .. +3,006`.
- CVM `.nc` grid: `x 303,000 .. 697,500`, `y 3,612,000 .. 3,903,000`, `z −45,000 .. 0`
  (264×195×181 @ 1500 m / 250 m).
- ⇒ The mesh extends **outside** the CVM grid in the far-field (x and y) and **above
  z=0** (real topography up to +3 km). The velocity model is regional and the converter
  refuses `z_max > sidecar z[-1]` (`convert_cvm_to_asagi.py:400`), so the grid CANNOT be
  extended to cover the topography. ASAGI clamp-to-boundary will supply edge/`z=0` values
  outside the footprint. This is physically acceptable (far-field is absorbing; clamping
  near-surface velocity is reasonable) but **invalidates the old "zero `out of range.
  Fixing.`" validation gate** (`README_port.md` gate #2). That gate must be revised.
- Nucleation hypocenter `(606971, 3707270, −4965.62)` and the constant stress tensor are
  coordinate-based; they are valid only if the hypocenter still lands on a fault facet —
  must be checked (Phase 0).

## Constraints

- **Cannot edit (reuse verbatim):** `safs_fault.yaml`, `safs_initial_stress.yaml`,
  `safs_material_cvm.yaml`, `safs_material_cvm.nc`. They are coordinate-based and
  physics-defining; changing them would change the validated RSSRW physics.
- **SeisSol reads PUML only** (`.puml.h5`), never `.msh`. A `pumgen` conversion is
  mandatory and **only builds/runs on Linux** (uint64_t/size_t clash on macOS arm64 —
  see `MESH_BUILD_README.md`). All `.msh`→`.puml.h5` work happens on a Linux box /
  Expanse / Frontera / the QuakeWorx gateway, never locally on the Mac.
- **Gmsh format:** the mesh is v2.2 ASCII (`$MeshFormat 2.2`), which is what pumgen
  `-s msh2` expects. Do not convert to v4.
- **Use exactly the user-named file** `safv4_deep_500m.msh` (1,180,159 tets), NOT
  `safv4_deep_500m_opt.msh` (1,416,939 tets) which also exists in `results/`.
- **Do not run a full mesh locally** (memory `feedback-no-local-mesh-runs`): local work is
  limited to read-only mesh parsing + the retag script + unit-style count checks.
  Production/smoke runs go to a cluster or QuakeWorx.
- **Preserve the original case.** Build the new run in a NEW sibling case folder; do not
  overwrite `safs_seisol_v2_1_0_RSSRW/safs_mesh.puml.h5` (the validated single-fault PUML).

## Proposed new case folder

```
project_7.0_alternative/seisol_quakeworx/safs_seisol_v2_2_0_RSSRW_safv4deep/
```
(name adjustable — keeps the v2 RSSRW lineage, flags the safv4-deep mesh). Contents:
unchanged yamls + `.nc` (copied/hardlinked), the new retagged `.msh`, the new
`.puml.h5`, an edited `parameters.par`, the new retag script, the pumgen build script,
and a revised README with updated gates.

---

## Phase 0 — Pre-flight geometric verification (read-only, local)

### Goal
Before producing anything, prove the new mesh is compatible with the reused
coordinate-based inputs: correct tag inventory, hypocenter on the fault, CVM-clamping
extent quantified, and the fault is not in the xy-plane (SeisSol right-handed requirement).

### Files to Create
- `safs_seisol_v2_1_0_RSSRW/preflight_deep_mesh.py` — standalone Python (run with
  `/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python3`), no third-party deps beyond
  `numpy`. Parses the v2.2 `.msh` directly (the file has no `$PhysicalNames`).

### Detailed Requirements
1. `parse_msh(path) -> (nodes: dict[int, (x,y,z)], tris: list[(tag, (n1,n2,n3))], tet_tags: Counter)`.
   Read `$Nodes` then `$Elements`; collect element type 2 (triangles, gmsh code 2) with
   their first physical tag, and element type 4 (tets) tag histogram.
2. **Tag inventory check**: assert triangle tags == `{101,102,103,201,202,203}` and tet
   tags == `{1}`. Print per-tag counts; FAIL if any unexpected tag.
3. **Bounding box** of all referenced nodes; print `x/y/z` min/max.
4. **CVM-coverage report**: hard-code the grid extent
   (`x[303000,697500] y[3612000,3903000] z[-45000,0]`); report the fraction of mesh
   nodes outside the grid in xy and the fraction with `z>0`. This is informational
   (clamping is accepted), not a gate — but print it so the run README can cite real
   numbers.
5. **Hypocenter-on-fault gate**: point `H=(606971.0, 3707270.0, -4965.62)`. Over all
   triangles with tag in `{101,102,103}` compute the minimum point-to-triangle distance
   (use exact point-to-triangle, or as a cheaper proxy the min distance to any fault-
   triangle vertex). GATE: `min_dist < 600.0` m (≈ one 500 m element). Print which strand
   the nearest facet belongs to (expected: 101 San Andreas). FAIL loudly otherwise — this
   is the single most important physics gate (nucleation will not trigger off-fault).
6. **Fault-not-in-xy-plane check**: for the fault triangles, compute face normals; assert
   the median `|n_z|` is small (faults near-vertical) — i.e. the fault is NOT horizontal.
   Informational print of the normal-z distribution; warn if any large horizontal fault
   patches exist (they would clash with SeisSol's free-surface assumption).

### Edge Cases to Handle
- Mesh has no `$PhysicalNames` block — parser must not depend on it.
- Node ids may be non-contiguous — index by id, not by position.
- Large file (88 MB, 279,448 nodes, 1.6 M elements) — stream line-by-line; do not load
  the whole file into a list twice.

### Acceptance Criteria
- [ ] Tag inventory == `{101,102,103,201,202,203}` (tris) and `{1}` (tets); counts match
      the table above.
- [ ] Hypocenter-on-fault `min_dist < 600 m`, nearest strand == 101.
- [ ] CVM-coverage and normal-z reports printed (no gate, but recorded for the README).
- [ ] Script exits 0 on PASS, non-zero on any FAIL.

### Dependencies
- Depends on: nothing. Required by: Phase 1 (do not proceed if Phase 0 fails).

---

## Phase 1 — New retag script + retagged `.msh` (local)

### Goal
Produce `safs_mesh_deep.msh`, geometry-identical to `safv4_deep_500m.msh` but with
physical tags remapped to SeisSol's gmsh convention for a single DR fault group.

### Files to Create
- `safs_seisol_v2_1_0_RSSRW/retag_safv4_deep_to_seissol.py` — a new retagger modeled on
  the existing `retag_msh_to_seissol.py` but with the new map and an inserted
  `$PhysicalNames` block (the input has none).

### Detailed Requirements
1. `TAG_MAP = {101: 103, 102: 103, 103: 103, 201: 101, 202: 105, 203: 105}` (volume tag 1
   untouched). Apply by reading each element's ORIGINAL first tag once and writing the
   mapped value (collision-safe; never re-read a written tag).
2. `NEW_PHYSICAL_NAMES = [(2,101,"free_surface"), (2,103,"fault"), (2,105,"absorbing"), (3,1,"rock")]`.
   Since the input has no `$PhysicalNames`, **insert** the block immediately after
   `$EndMeshFormat` (gmsh requires `$PhysicalNames` before `$Nodes`). If a block ever IS
   present (future inputs), replace it (keep the original script's replace branch).
3. Only element type-2 (triangle) physical tags are remapped; type-4 (tet) volume tag 1 is
   passed through unchanged; geometry, connectivity, node ordering byte-preserved.
4. CLI: `retag_safv4_deep_to_seissol.py <in.msh> <out.msh>`; print the retagged
   surface-element counts by NEW tag with their pumgen BC (tag−100) and meaning, exactly
   like the original script's summary.
5. Run it:
   `python3 retag_safv4_deep_to_seissol.py \
      ../../../project_7.0_preferred/meshing/results/safv4_deep_500m.msh \
      <newcase>/safs_mesh_deep.msh`

### Edge Cases to Handle
- Input without `$PhysicalNames` (the actual case) → insert block.
- A triangle whose first tag is not in `TAG_MAP` → leave unchanged but COUNT and WARN
  (should be zero given Phase 0).
- Lines with `ntags == 0` (shouldn't occur in v2.2 with physical tags) → pass through.

### Acceptance Criteria
- [ ] Output `safs_mesh_deep.msh` byte-identical to input except `$PhysicalNames` (new)
      and the first tag of each surface element (verify with a diff that touches only
      `$PhysicalNames` + `$Elements` triangle lines).
- [ ] Printed NEW-tag counts: `103 (DR) = 132,333`, `101 (free surface) = 175,033`,
      `105 (absorbing) = 122,593`.
- [ ] Re-running Phase 0's parser on the OUTPUT shows tris `{101,103,105}`, tets `{1}`,
      same node count (279,448), same bbox.

### Dependencies
- Depends on: Phase 0 PASS. Required by: Phase 2.

---

## Phase 2 — pumgen conversion to PUML (Linux / cluster)

### Goal
Convert `safs_mesh_deep.msh → safs_mesh_deep.puml.h5` with `pumgen`, and verify the PUML
boundary codes.

### Files to Use (reuse, copy into new case folder)
- `build_and_run_pumgen_linux.sh` — unchanged; it builds the official PUMGen and runs
  `pumgen -s msh2 <msh>` producing `<msh-stem>.puml.h5`.

### Detailed Requirements
1. On the Linux machine / Expanse / Frontera (NOT the Mac):
   `bash build_and_run_pumgen_linux.sh /abs/path/safs_mesh_deep.msh`
   → produces `safs_mesh_deep.puml.h5` (+ `.xdmf`).
2. **Verify boundary codes** with the h5py snippet from `MESH_BUILD_README.md` Step 3,
   adapted to the new expected counts:
   - codes ⊆ `{0 interior, 1 free surface, 3 dynamic rupture, 5 absorbing}`
   - `3` ≈ 132,333 (the merged fault group), `1` ≈ 175,033, `5` ≈ 122,593
   - `connect.shape[0]` ≈ 1,180,159 elements, `geometry.shape[0]` ≈ 279,448 nodes.
3. **Right-handedness / negative-Jacobian**: pumgen reorders tets to positive volume;
   capture its log. If pumgen reports flipped/zero-volume elements beyond a negligible
   count, STOP and report (do not silently proceed) — the tetgen mesh passed MFEM's load
   in the RUNLOG, so this is expected clean, but confirm.

### Edge Cases to Handle
- pumgen output name defaults to `<stem>.puml.h5` → `safs_mesh_deep.puml.h5`; make sure
  `parameters.par` references that exact name (Phase 3).
- If the cluster build of pumgen is unavailable, the QuakeWorx gateway / an existing
  Expanse pumgen-build env may be reused — do not attempt to build on macOS.

### Acceptance Criteria
- [ ] `safs_mesh_deep.puml.h5` exists; h5 verification prints codes ⊆ {0,1,3,5} with the
      per-code counts above (±0 on the boundary triangle totals).
- [ ] Element/node counts match the `.msh`.
- [ ] pumgen log shows no negative-Jacobian abort.

### Dependencies
- Depends on: Phase 1. Required by: Phase 3.

---

## Phase 3 — Assemble the new case folder + edit `parameters.par`

### Goal
A self-contained, runnable SeisSol case directory pointing at the new PUML, with all
reused inputs present and the documentation/gates updated for the new mesh.

### Files to Create (new folder `safs_seisol_v2_2_0_RSSRW_safv4deep/`)
- `parameters.par` — copy of the RSSRW `parameters.par` with ONE functional edit:
  `MeshFile = 'safs_mesh_deep.puml.h5'` (line 80). Everything else (FL=103, RS constants,
  nucleation `s_0/t_0`, `XRef/YRef/ZRef`, `CFL=0.5`, `FixTimeStep=0.1`, `Format=10`,
  `OutputMask`, `EndTime=200`) stays identical.
- `safs_fault.yaml`, `safs_initial_stress.yaml`, `safs_material_cvm.yaml` — copied verbatim.
- `safs_material_cvm.nc` — copied (or hardlinked to save 112 MB; QuakeWorx upload needs a
  real file, so copy for a gateway run).
- `safs_mesh_deep.msh` (from Phase 1) and `safs_mesh_deep.puml.h5` (from Phase 2).
- `retag_safv4_deep_to_seissol.py`, `build_and_run_pumgen_linux.sh` — copied for provenance.
- `README_deep.md` — new run guide (see below).
- `output/` — created empty before running (`OutputFile = 'output/safs'` requires it).

### Detailed Requirements
1. The ONLY edit to `parameters.par` is the `MeshFile` line. Do not touch any physics or
   discretization namelist. (Mesh size ≈ same tet count as the old run, so `FixTimeStep`/
   `CFL` need no change; LTS handles the spatial variation.)
2. `README_deep.md` must record, with the real numbers from Phase 0:
   - the new single-group fault mapping and the retag command;
   - the pumgen command + boundary-code verification expectation;
   - the **revised material gate** (see Phase 4);
   - the **revised strike/dip gate** (see Phase 5);
   - a one-line provenance note: "mesh = `project_7.0_preferred/.../safv4_deep_500m.msh`,
     retagged 101/102/103→DR(single group)".
3. Confirm relative paths inside the yamls resolve from the new folder: `safs_fault.yaml`
   `!Include`s `safs_initial_stress.yaml` (relative — OK after copy); `safs_material_cvm.yaml`
   references `file: safs_material_cvm.nc` (relative — OK after copy).

### Acceptance Criteria
- [ ] `diff` of new vs old `parameters.par` shows exactly one changed line (`MeshFile`).
- [ ] All `!Include` / `file:` relative references resolve inside the new folder.
- [ ] Folder is self-contained (no path escapes except the documented mesh provenance).

### Dependencies
- Depends on: Phase 2. Required by: Phase 5.

---

## Phase 4 — Revise the material / coverage validation gate

### Goal
Replace the now-invalid "zero `out of range. Fixing.`" gate with one that reflects the
new mesh extending beyond the CVM grid.

### Files to Modify
- `README_deep.md` (created in Phase 3) — validation-gates section.

### Detailed Requirements
1. State explicitly that ASAGI `out of range. Fixing.` warnings ARE now expected, and
   why: the mesh far-field (x outside [303k,697.5k], y outside [3612k,3903k]) and the
   topography band (`z>0`) lie outside the CVM grid; ASAGI clamps to the nearest grid
   boundary there. Cite the Phase-0 coverage fractions.
2. New gate wording: clamp warnings must occur ONLY for nodes outside the documented
   footprint / above z=0. The near-fault region (which IS inside the grid) must report
   correct, non-clamped moduli. Spot-check: log the material at the hypocenter facet and
   confirm it matches the CVM value there (μ ≈ 23.4 GPa per the fault-yaml header / CALC
   doc), NOT the constant-fallback (3.2e10) and NOT a clamped edge value.
3. Keep the A/B switch note: reverting `MaterialFileName` to `safs_material.yaml`
   (constant medium) still works and is the debugging fallback.
4. Do NOT regenerate the `.nc` (impossible — the converter rejects `z_max`/`xy` beyond the
   regional sidecar footprint; `convert_cvm_to_asagi.py:400`). Record this rationale so a
   future reader doesn't try.

### Acceptance Criteria
- [ ] README documents expected clamping + the spot-check; old "zero out of range" wording
      removed.
- [ ] (At run time, Phase 5) hypocenter-facet μ ≈ 23.4 GPa, not 32 GPa fallback.

### Dependencies
- Depends on: Phase 0 (coverage numbers), Phase 3. Required by: Phase 5.

---

## Phase 5 — Smoke run + physics validation gates (cluster / QuakeWorx)

### Goal
Confirm the case initializes and produces physically correct on-fault tractions and
nucleation on the new geometry. Short `EndTime` first.

### Detailed Requirements
1. Run on the cluster / QuakeWorx with the existing build (CONVERGENCE_ORDER ≥ 2, ASAGI
   support — same build that ran the single-fault case). For a smoke, optionally lower
   `EndTime` (e.g. 10–20 s) to reach the time loop and first propagation cheaply, then
   restore `EndTime=200` for production.
2. **Init gate**: run reaches the time loop with no easi/ASAGI exception; the log line
   `RS parameter source (...): f0 1 - muW 1 - b 0` shows `b 0` (spatial `rs_b` from easi —
   a constant-b run is WRONG, per `safs_fault.yaml` and memory `safs-seissol-rssw-dc-calc`).
3. **Material gate** (Phase 4): clamp warnings only outside the footprint; hypocenter-facet
   μ ≈ 23.4 GPa.
4. **Initial-traction / stress-port gate** (`README_port.md`): at the hypocenter facet
   `Pn0` compressive with `|Pn0|` ≈ 49 MPa; `Ts0/Td0` consistent sign across the fault. If
   `|Pn0|` matches but `Ts0/Td0` flip across the fault, the strike/dip frame is mirrored →
   fix via `XRef/YRef/ZRef` (REVIEW R-001), NOT by editing stresses.
5. **Strike/dip-frame gate (NEW, multi-strand specific)**: the single reference vector
   `XRef/YRef/ZRef = (0,-1,0)`, `refPointMethod=1` was validated for the single SAF
   strand. With Banning + Garnet at different strikes, confirm each strand's `Ts0/Td0`
   signs are physical (right-lateral on the SAF). If a strand comes out mirrored, this is
   the expected place to iterate (document any change). This is the main NEW physics risk
   introduced by merging 3 differently-oriented strands into one group.
6. **Nucleation gate**: rupture initiates inside the ~6 km Gaussian patch on the San
   Andreas strand during the first ~1 s (`Tnuc_s` smoothStep), then propagates
   spontaneously. (Phase 0 already proved the patch centre is on-fault.)

### Acceptance Criteria
- [ ] Init reaches the time loop; `b 0` (easi) in the log.
- [ ] Material spot-check passes (Phase 4).
- [ ] `|Pn0|` ≈ 49 MPa, traction signs physical on all strands; any XRef fix documented.
- [ ] Nucleation triggers on the SAF patch and rupture propagates.
- [ ] No NaN / dt→0 / MFEM-abort (cf. memory `frontera-mpi-waitall-transient-fault` for
      distinguishing fabric blips from real failures).

### Dependencies
- Depends on: Phases 2–4. Required by: production run (EndTime=200).

---

## Testing Strategy

- **Phase 0/1 (local, deterministic):** the pre-flight parser and a re-parse of the
  retagged `.msh` are the unit tests — exact tag counts and bbox are the oracle.
- **Phase 2 (Linux):** the h5py boundary-code histogram is the oracle (counts derived from
  Phase 0). No solver needed.
- **Phase 5 (cluster):** physics gates compare against the single-fault RSSRW run's known
  quantities (`|Pn0|`≈49 MPa, hypocenter μ≈23.4 GPa, nucleation in ≤1 s). These are the
  same gates `README_port.md` defines, plus the new multi-strand strike/dip check.

## Risk Assessment

1. **Retag tag collision (101→103 and 103→103).** Mitigation: read original tag once,
   write mapped value; never re-process. Verified by Phase 1 re-parse showing `{101,103,105}`.
2. **Hypocenter off-fault** (different meshing pipeline than the single-fault mesh) →
   nucleation never triggers. Mitigation: Phase 0 hard gate (`min_dist < 600 m`). This is
   the highest-impact risk and is caught before any cluster time is spent.
3. **Strike/dip frame mismatch on Banning/Garnet** — one reference vector for three
   strikes may mirror a strand's traction. Mitigation: Phase 5 gate #5; lever is
   `XRef/YRef/ZRef` only. Does not block the SAF (the nucleating strand).
4. **ASAGI clamping mistaken for a bug.** Mitigation: Phase 4 rewrites the gate; clamping
   in the far-field/topography is expected, the near-fault region must be un-clamped
   (spot-check μ).
5. **pumgen only on Linux.** Mitigation: explicit in constraints; never attempt on macOS.
   Reuse an existing cluster pumgen build if available.
6. **Wrong mesh file** (`_opt` vs base; stale `gate_conformity_deep_base.json` reports
   different counts than the actual file). Mitigation: use exactly `safv4_deep_500m.msh`;
   Phase 0 prints the bbox/counts of the ACTUAL file so the operator can confirm
   (actual file = 1,180,159 tets / 279,448 nodes, NOT the stale gate JSON's 647,494).
7. **Larger absorbing far-field changes wave behavior / cost.** Informational: tet count
   is similar to the old run, so `FixTimeStep`/`CFL` are unchanged; LTS absorbs the
   spatial variation. Watch the first cluster run's reported min timestep.

## What is explicitly NOT changing

- No edits to `safs_fault.yaml`, `safs_initial_stress.yaml`, `safs_material_cvm.yaml`,
  `safs_material_cvm.nc`, or any RS/nucleation/discretization parameter.
- No `.nc` regeneration.
- No new SeisSol build (reuse the ASAGI-enabled, order ≥2 build).
- The original `safs_seisol_v2_1_0_RSSRW/` case is left intact.
```
