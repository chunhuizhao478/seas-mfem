# WORKFLOW: regional stress → on-fault traction VTU

**Target artefact (this doc explains exactly how it is produced):**

```
stress/results/safs_fault_box_nwcut_500m_lcfar3000_z0embed_hz1671_psi37_az351/
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_hz1671_psi37_az351_fault_stress.vtu
    └ sibling: ..._bulk_stress.vtu, ..._summary.json
```

Written 2026-05-28. This is a **per-run map** for the `hz1671_psi37_az351`
fault-stress projection on the 500 m `z0embed` mesh. It is the run-specific
companion to `stress/README.md` (the general pipeline reference) — read the
README §Workflow for the multi-variant / batch view; read this for *what the
filename tags mean*, *the exact command that made this file*, and *the numbers
that came out*.

The whole stress pipeline is **pure Python** (`meshio` + `numpy`) and produces
only visualization / IC-setup artefacts. No C++ runtime is involved in making
this VTU (the C++ `stress_safs.h5` sidecar is a separate, optional branch — see
README §C).

---

## TL;DR — the one command that produced the target file

```bash
conda activate pythonenv
cd stress/code

python project_to_fault_stress.py \
    ../../meshing/results/vtu/safs_fault_box_nwcut_500m_lcfar3000_z0embed_fault.vtu \
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_hz1671_psi37_az351 \
    --write-bulk \
    --SHmax-az 351
```

- **arg 1** = input fault VTU (from the meshing pipeline, see Stage 0).
- **arg 2** = `mesh_base` — the output stem *and* the output sub-directory name.
- `--write-bulk` → also emit the `_bulk_stress.vtu` (σ⁰ on every tet centroid).
- `--SHmax-az 351` is the **only** parameter override. Everything else
  (`SHmax/Shmin/Sv/P_p` magnitudes and the `--strike-hint-az 314`) is left at the
  built-in Hickman & Zoback 1671 m defaults (`DEFAULTS = _hz_demo_safod_defaults()`,
  `project_to_fault_stress.py:1480`).

> **Why the output dir is the full base name** (not `500m_lcfar3000/`):
> `_extract_lc_tag()` matches `_<N>m` followed by **at most one** variant suffix
> (`r"_(\d+m(?:_[A-Za-z0-9]+)?)$"`). The stem `..._500m_lcfar3000_z0embed_hz1671_psi37_az351`
> has many trailing segments and ends in `az351`, so the regex misses and the
> function **falls back to the full `mesh_base`**. That fallback is what created the
> long-named result directory.

---

## Filename tag decode

`safs_fault_box_nwcut_500m_lcfar3000_z0embed` `_hz1671` `_psi37` `_az351`

| Segment | Meaning | Source |
|---|---|---|
| `safs_fault_box_nwcut` | SAFS box mesh with the NW hard-cut applied | meshing pipeline |
| `500m` | 500 m near-fault target resolution | meshing |
| `lcfar3000` | far-field background size `lc = 3 km` | meshing |
| `z0embed` | fault cut **exactly at z = 0** (`run_z0cut_meshing.py`) | meshing Stage 5 |
| **`hz1671`** | **H**ickman & **Z**oback (2004) SAFOD magnitudes at **1671 m** depth: `SHmax=113, Shmin=49, Sv=45, P_p=16` MPa | `…_projection.py:328,388` |
| **`psi37`** | **Ψ ≈ 37°** — angle between SHmax and the SAF strike (`\|351 − 314\| = 37`). The *intermediate/shallow* value from H&Z's depth table, vs the deep Ψ = 69° used in `demo_safod()` | `…_projection.py:379` |
| **`az351`** | **SHmax azimuth = 351°** cw-from-N (= N9°W). Chosen so that, against the SAF strike hint 314°, Ψ = 37° | `--SHmax-az 351` |

The tags are an operator-supplied label baked into `mesh_base`; the script does
**not** parse them. They are only descriptive — the actual numeric parameters are
the ones echoed in `_summary.json["params"]`.

---

## Pipeline at a glance

```
meshing/results/vtu/…_z0embed_fault.vtu  ─┐  (triangulated SAF surface)
                          (+ _bulk.vtu)  ─┤
                                          │
   H&Z 1671 m regional params  ──────────┤
   SHmax=113 Shmin=49 Sv=45 P_p=16        │
   SHmax_az=351  fault_strike_hint=314    │
                                          ▼
           project_to_fault_stress.py  (single-file mode)
                                          │
   (1) build_bulk_stress_tensor(σ⁰)  3×3 Cauchy in (E,N,Up), compression POSITIVE
   (2) triangle_geometry             per-tri centroid + unit normal
   (3) build_fault_basis             Tandem (s = up×n̂, d = s×n̂, n̂), strike→314° hint
   (4) resolve_traction(σ⁰,n̂,s,d,P_p) σ_n, σ_n_eff, τ_strike, τ_dip, |τ|, rake, μ_app
   (5) write_fault_vtu / write_bulk_vtu / write_summary_json
                                          ▼
   results/<full mesh_base>/…_fault_stress.vtu  + _bulk_stress.vtu + _summary.json
```

---

## Stage 0 — input mesh (from the meshing pipeline)

The projector consumes the **ParaView VTU split** of the volume `.msh`
(produced by `msh_to_vtu.py`, PIPELINE Stage 7) — two siblings under
`meshing/results/vtu/`:

- `safs_fault_box_nwcut_500m_lcfar3000_z0embed_fault.vtu`
  — the 42 358-triangle SAF surface (Physical Surface 101); this is **arg 1**.
- `safs_fault_box_nwcut_500m_lcfar3000_z0embed_bulk.vtu`
  — the bulk tetra mesh; auto-found by `_derive_bulk_path()` for `--write-bulk`.

The fault triangulation is verbatim from the mesher, so the on-fault geometry
(and therefore the projected tractions) is fully determined upstream. See
`meshing/docs/PIPELINE_raw_data_to_mesh.md` for how this mesh was built.

---

## Stage 1 — regional stress parameters (H&Z 2004 SAFOD, 1671 m)

From `_summary.json["params"]` of the target run:

| Parameter | Value | Note |
|---|---|---|
| `SHmax_MPa` | 113.0 | max horizontal principal stress |
| `Shmin_MPa` | 49.0 | min horizontal principal stress |
| `Sv_MPa` | 45.0 | vertical stress at 1671 m (ρ g z) |
| `P_p_MPa` | 16.0 | hydrostatic pore pressure at 1671 m |
| `SHmax_azimuth_deg` | **351.0** | the only override (`--SHmax-az 351`) |
| `fault_strike_azimuth_hint_deg` | 314.0 | SAF strikes N46°W (default) |
| `depth_model` | `constant` | σ⁰ depth-independent; all gradients = 0 |
| `rake_sense` | right-lateral | |

Because `depth_model = constant` with zero gradients, **σ⁰ is the same tensor at
every point** — all spatial variation in the on-fault tractions comes purely from
the fault normal rotating along the curved SAF, not from depth.

---

## Stage 2 — bulk Cauchy tensor σ⁰ (`build_bulk_stress_tensor`)

Assembles the 3×3 stress in the `(east, north, up)` frame: a horizontal
`[SHmax, Shmin]` block rotated by the SHmax azimuth, with `Sv` on the vertical
axis. A **single source-site sign flip** converts the geological
compression-negative convention to the project-internal **compression POSITIVE**
(SEAS) convention. For this run, σ⁰ at z = 0 (`_summary.json`):

```
        E          N        Up
E  [  50.566   -9.889     0.0  ]
N  [  -9.889  111.434     0.0  ]   MPa  (compression POSITIVE)
Up [   0.0      0.0      45.0  ]
```

(Off-diagonal `−9.889` is the SHmax/Shmin shear coupling from the 351° azimuth;
the vertical row/col is decoupled with `Sv = 45`.)

---

## Stage 3 — per-triangle fault basis (`triangle_geometry` + `build_fault_basis`)

For each of the 42 358 triangles: centroid + **unit outward normal** `n̂`
(`CalcOrtho`-style cross product), then the **Tandem fault basis**
(CLAUDE.md "Fault-local tangent frame"):

```
s = up × n̂   (strike)      d = s × n̂   (down-dip)      n̂   (normal)
```

Strike orientation is harmonised against the **314° hint** so all `s` vectors
point consistently along the SAF (sign-flips the per-triangle `s` where its
azimuth disagrees with the hint by > 90°). `n_degenerate_basis = 0` for this run
— every triangle yielded a valid frame.

---

## Stage 4 — resolve traction (`resolve_traction`)

Traction vector `t = σ⁰ · n̂`, then project onto the fault frame:

| Output | Formula | Sign meaning |
|---|---|---|
| `sigma_n_total` | `n̂ · t` | > 0 = compression |
| `sigma_n_eff` | `sigma_n_total − P_p` | effective normal stress |
| `tau_strike` | `s · t` | **> 0 = right-lateral** |
| `tau_dip` | `d · t` | > 0 = down-dip |
| `tau_magnitude` | `√(τ_strike² + τ_dip²)` | |
| `rake_deg` | `atan2(τ_dip, τ_strike)` | |
| `mu_apparent` | `\|τ\| / sigma_n_eff` | apparent friction |

---

## Stage 5 — outputs

Three files in `results/safs_fault_box_nwcut_500m_lcfar3000_z0embed_hz1671_psi37_az351/`:

### `…_fault_stress.vtu` (24 MB) — the target artefact

ParaView triangle surface with two data layers:

- **Cell data** (per triangle, exact projection):
  `sigma_n_total_MPa_cell`, `sigma_n_eff_MPa_cell`, `tau_strike_MPa_cell`,
  `tau_dip_MPa_cell`, plus the per-cell basis vectors `strike_vec`, `dip_vec`,
  `normal_vec`, `traction_vec_MPa`, and the `gmsh:physical` / `gmsh:geometrical`
  tags.
- **Point data** (area-weighted node averages via `cell_to_node_average`):
  `sigma_n_total_MPa`, `sigma_n_eff_MPa`, `tau_strike_MPa`, `tau_dip_MPa`,
  `tau_magnitude_MPa`, `rake_deg`, `mu_apparent`.

> Node fields are the **area-weighted average of the cell projections**, *not* a
> re-projection on a node-averaged basis — the two differ on a faceted fault.
> The verifier (Stage 6) replays this exact averaging so its point-data check
> matches the writer's convention.

### `…_bulk_stress.vtu` (92 MB)

The 3 693 060-tet bulk mesh with cell-data `sigma_tensor_MPa` (full 6-component
symmetric σ⁰ at each tet centroid). For `constant` depth-model this is the same
tensor everywhere.

### `…_summary.json`

Parameter echo, σ⁰ at z = 0, cell counts, and min/median/max per field.

---

## This run's results (`_summary.json["stats"]`)

| Field | min | median | max |
|---|---:|---:|---:|
| `sigma_n_total_MPa` | 46.37 | 73.72 | 112.96 |
| `sigma_n_eff_MPa`   | 30.37 | 57.72 | 96.96 |
| `tau_strike_MPa`    | −9.57 | 26.47 | 31.49 |
| `tau_dip_MPa`       | −33.94 | −16.01 | 8.05 |
| `tau_magnitude_MPa` | 1.53 | 30.61 | 33.98 |
| `rake_deg`          | −177.92 | −29.54 | 176.91 |
| `mu_apparent`       | 0.016 | 0.562 | 0.641 |

- `fault_n_cells = 42358`, `bulk_n_cells = 3693060`, `fault_n_degenerate_basis = 0`,
  zero NaN vertices.
- Median μ_apparent ≈ 0.56 reflects Ψ = 37° (well-oriented, near-Byerlee): far
  from the Ψ = 69° deep-SAF "weak fault" case (μ_app ≈ 0.24 in `demo_safod`).
  The non-zero `tau_dip` median (−16 MPa) is the transitional SS↔reverse signature
  of the literal 1671 m magnitudes (`Sv = 45 < Shmin = 49`).

---

## Conventions (carried end-to-end)

- Frame: UTM Zone 11 N metres, `(x = east, y = north, z = up)`, `z ≤ 0` underground.
- σ⁰: **compression POSITIVE** (SEAS internal) after the single source-site flip.
- `sigma_n_eff = sigma_n_total − P_p`, `P_p > 0`.
- Tandem fault basis `s = up × n̂`, `d = s × n̂`; `tau_strike > 0` ⇒ right-lateral.
- Units: **MPa** in all VTU/JSON; Pa only in the C++ `stress_safs.h5` sidecar.

(Pinned `_CONVENTION_STRING`, echoed verbatim in `_summary.json["convention"]`.)

---

## Stage 6 — verification (optional but recommended)

`verify_onfault_stress.py` re-evaluates σ⁰ analytically at **every** cell
centroid and vertex and compares against the stored VTU values (full population,
not a sample). For the long-named `z0embed` run, point the `--*-vtu` flags at the
files in this directory:

```bash
conda activate pythonenv
cd stress
BASE=safs_fault_box_nwcut_500m_lcfar3000_z0embed_hz1671_psi37_az351
python code/verify_onfault_stress.py \
    --fault-vtu    results/$BASE/${BASE}_fault_stress.vtu \
    --bulk-vtu     results/$BASE/${BASE}_bulk_stress.vtu \
    --summary-json results/$BASE/${BASE}_summary.json \
    --report-json  results/$BASE/verify_report.json
```

Tolerances: 1 Pa for bulk + fault cell-data, 1e3 Pa for fault point-data. The
documented 12-variant set passes 12/12 fields at machine round-off (~1e-7 Pa on a
~1e8 Pa magnitude); the `z0embed` run uses the identical writer/verifier path.

---

## References

- `stress/README.md` — general (multi-variant, batch) pipeline reference + the 12-variant inventory.
- `stress/code/hickman_and_zoback_regional_stress_projection.py` — H&Z analytic core (`build_bulk_stress_tensor`, `fault_basis_vectors`, `resolve_traction`, `demo_safod`, the 1671 m docstrings).
- `stress/code/project_to_fault_stress.py` — the Phase 1–4 driver (CLI, `_extract_lc_tag`, writers).
- `stress/code/verify_onfault_stress.py` — Phase-8 analytic verifier.
- `stress/docs/PLAN_onfaultstress.md` (+ `_fix.md`) — accepted plan and decision trail.
- `meshing/docs/PIPELINE_raw_data_to_mesh.md` — how the `z0embed` input mesh is built.
- Hickman, S. & Zoback, M. (2004), *Stress orientations and magnitudes in the SAFOD pilot hole*, GRL — `stress/reference/`.
