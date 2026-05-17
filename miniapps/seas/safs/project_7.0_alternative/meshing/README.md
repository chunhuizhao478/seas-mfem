# meshing/

Build the SAFS bulk + embedded-fault mesh from the CFM ALT6 fault
representation.

```
raw/                CFM .ts fault surfaces + their direct .stl conversions
                    (input to clean_freesurface_mesh.py)

code/               ts_to_stl.py            .ts -> _unclipped.stl
                    clean_freesurface_mesh.py
                                            top-trace cleanup, isotropic remesh
                    nw_cut_strip.py         hard-cut NW tail against the
                                            project_7.0_preferred anchor
                    run_nwcut_meshing.py    bbox + Gmsh template + size-field driver
                    safs_fault_box*.geo     Gmsh templates (active: *_nwcut.geo)
                    msh_to_vtu.py           .msh -> {bulk,fault}.vtu
                    check_mesh_quality.py   Q1/Q2 verifier (PLAN_mesh_quality.md Phase 4)
                    test_mesh_quality.py    pytest gates for Q1/Q2 + embedding + zmax
                    test_nw_cut_strip.py    pytest suite (NW-cut strip)

docs/               PLAN_mesh_quality.md    accepted plan + decision trail (Phases 1-5)
                    EXPLORE_*.md, PLAN_corefine_pipeline.md, PLAN_retriangulate.md,
                    PLAN_nw_hard_cut.md, NW_HARD_CUT_*.md

results/
    stl_cleaned/    output of clean_freesurface_mesh.py (3 STLs @ 500/1000/2000 m,
                                                          each with _clip + _vdel)
    stl_nwcut/      output of nw_cut_strip.py (3 STLs, plus _zclamp100m siblings
                                                created on-the-fly by the driver)
    msh/            6 canonical .msh files:
                        safs_fault_box_nwcut_{500,1000,2000}m_{zgraded,lcfar3000}.msh
                    (size_field.pos files appear here only with --gen-size-field)
    vtu/            12 matching {bulk,fault}.vtu files, written by msh_to_vtu.py
                    immediately after each .msh
```

## ⚠️ Free-surface invariant: `mesh_zmax == 0`

The mesh's **top face is the free surface** and **must sit at z = 0
exactly** — the SEAS elasticity operator applies the Neumann (traction-
free) BC on `mesh_zmax`, and the friction law assumes the fault is
embedded in a half-space whose surface is the z = 0 plane.  A mesh
top at z > 0 (e.g. z = +100 m above the geoid) silently breaks the
physics: tractions are evaluated on an aerial slab that should not
exist, and downstream comparison against any model that treats z = 0
as the free surface (CVM-H, CSM, every SCEC benchmark) is biased.

`run_nwcut_meshing.py` defaults respect this invariant:

- `--fault-top-clamp` (default **100 m**, see PLAN_mesh_quality.md
  Phase 1) snaps any fault vertex with z > −100 m down to z = −100 m.
  This eliminates the coplanar-with-top-face PLC artefact AND gives
  gmsh enough vertical headroom near the free surface to avoid the
  sliver/short-edge pinch in the fault-corridor × free-surface
  wedge.  The 100 m clearance is well below the smallest SAFS
  seismogenic length (~1 km) so co-seismic surface rupture is
  preserved in practice.  **Physics ceiling**: 250 m — above that
  the driver prints a warning because near-surface rupture would be
  suppressed.
- `--fault-top-clamp-policy` (default `fixed`, set to `auto` for a
  per-family adaptive clamp = `clip(surface_lc, 100, 250)` where
  `surface_lc` is `z_lc_top` in `--z-graded` mode, else `lc_near`).
- `--pad-top` defaults to `--fault-top-clamp` (i.e. **100 m**), so
  the box top lands at `clamped_fault_zmax + pad_top = −100 + 100 = 0` ✓.

**Do NOT pass `--pad-top` ≠ `--fault-top-clamp`.** Doing so either
lifts `mesh_zmax` above 0 (aerial slab) or drops it below 0 (the
free-surface BC is then applied below the geoid).  If a sidecar built
with `--extend-z-top` is consumed by the projector, the extend value
must equal `pad_top` so the projector's containment check stays
satisfied.

Verify after every mesh build:

```bash
python -c "import meshio; m = meshio.read('<path>.msh'); print('zmax =', m.points[:,2].max())"
# Must print: zmax = 0.0
```

## Workflow

A single driver — `run_nwcut_meshing.py` — produces every mesh
variant.  The shared Gmsh template `safs_fault_box_nwcut.geo` exposes
`lc_far`, `USE_Z_GRADED`, `USE_SIZE_FIELD` (and the Phase 1–3 quality
knobs) via `-setnumber`, so each variant is a CLI-flag combination
rather than a separate script.

Outputs land in:

- `meshing/results/msh/<stem>.msh`
- `meshing/results/vtu/<stem>_bulk.vtu`
- `meshing/results/vtu/<stem>_fault.vtu`

where `<stem> = safs_fault_box_nwcut_<RES>m<SUFFIX>`.  By default
every invocation meshes all three resolutions (500 / 1000 / 2000 m);
restrict with `--res 500` or `--res 500 1000`.

### Canonical accepted-set commands

Two variant families, three resolutions each → six target meshes.
Every quality knob ships at the accepted default (see "Reproducibility"
below), so the only flags you must specify are the family identifier
and the suffix:

```bash
conda activate pythonenv
cd meshing/code

# zgraded family — depth-graded background size field
python run_nwcut_meshing.py --z-graded --suffix _zgraded \
    --res 500 1000 2000 --quiet

# lcfar3000 family — uniform 3 km far-field cap, no z-grading
python run_nwcut_meshing.py --lc-far 3000 --suffix _lcfar3000 \
    --res 500 1000 2000 --quiet
```

After every regeneration, gate the result against Q1/Q2:

```bash
python check_mesh_quality.py \
    ../results/msh/safs_fault_box_nwcut_500m_zgraded.msh \
    ../results/msh/safs_fault_box_nwcut_1000m_zgraded.msh \
    ../results/msh/safs_fault_box_nwcut_2000m_zgraded.msh \
    ../results/msh/safs_fault_box_nwcut_500m_lcfar3000.msh \
    ../results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh \
    ../results/msh/safs_fault_box_nwcut_2000m_lcfar3000.msh
# Markdown table to stdout; exit 0 iff all six pass Q1 and Q2.
# Accepted state (see Reproducibility table below) exits 1 because
# 500 m variants miss Q1 (38.72 m vs 100 m floor) and every variant
# has 1-2 residual slivers under Q2.
```

### Other (historical) variants

Kept for reference / cross-comparison; **not part of the accepted
mesh set** and not regenerated by default.  These use older quality
defaults and may fail Q1/Q2:

```bash
# nwcut baseline: lc_far = 10 km, no z-grading (uniform far-field)
python run_nwcut_meshing.py

# lcfar5000: lc_far = 5 km, no z-grading
python run_nwcut_meshing.py --lc-far 5000 --suffix _lcfar5000
```

### Velocity-driven size field (G-1, optional)

When you want the mesher to refine inside the fault zone proportional
to a CVM velocity contrast (rather than uniformly via `lc_near`), run
with `--gen-size-field` and one of the per-version sidecars:

```bash
python run_nwcut_meshing.py --res 500 1000 2000 --gen-size-field \
    --sidecar ../../velocity/results/cvmh/velocity_safs.h5
```

The sidecar lives under one of:

- `../../velocity/results/cvmh/velocity_safs.h5`
- `../../velocity/results/cvm_s4.26.m01/velocity_safs.h5`
- `../../velocity/results/multiscale_statewise_cvm/velocity_safs.h5`

This activates `USE_SIZE_FIELD=1` in the .geo and emits
`results/msh/<stem>.size_field.pos` alongside the .msh.  Note: this
mode has not been re-validated against Q1/Q2 since the Phase 1–3
quality changes.

### Other useful flags

- `--threads N` — gmsh thread count (default: half of cpu_count).
- `--pad-x`, `--pad-y`, `--pad-bottom`, `--pad-top` — override domain
  padding (defaults: 50 km / 50 km / 25 km / `--fault-top-clamp`
  (= 100 m), which puts `mesh_zmax` at z = 0 — see the
  free-surface-invariant note at the top of this file).
- `--velocity-bbox` / `--sidecar` / `--double-domain` — auto-fit the
  bulk box to the velocity dataset's UTM AABB instead of using fixed
  padding.
- `--lc-near`, `--dist-inner`, `--dist-outer` — tune the
  fault-distance threshold field that produces `lc_near` inside the
  fault corridor.
- `--lc-floor` (default **100 m**) — hard lower bound on the local
  mesh size applied via `Field[Max]` in the .geo (PLAN_mesh_quality.md
  Phase 2; forwarded as the gmsh variable `lc_min`).  Counteracts
  gmsh's embedded-discrete-surface code path which would otherwise
  produce sub-100 m bulk edges near the fault top.  Must be
  `< --lc-near` or the driver aborts.  Not the same as `--lc-min`
  (which clamps the `build_pos` size-field generator and is used
  only with `--gen-size-field`).
- `--quiet` — suppress gmsh stdout.

### Quality gates (Q1 / Q2)

Every produced mesh is expected to clear two hard bars (see
`docs/PLAN_mesh_quality.md`):

| Gate | Definition |
|---|---|
| **Q1** | minimum bulk tet edge length `>= 100 m` |
| **Q2** | Joe-Liu tet quality `eta > 0.1` for **all** bulk tets |

The canonical regression check is `code/check_mesh_quality.py`:

```bash
cd meshing/code
python check_mesh_quality.py ../results/msh/safs_fault_box_nwcut_*_zgraded.msh \
                              ../results/msh/safs_fault_box_nwcut_*_lcfar3000.msh
echo "exit=$?"   # 0 iff every mesh passes Q1 and Q2
```

`code/test_mesh_quality.py` wraps the same checks plus embedding +
free-surface invariants as a pytest module (skip-on-missing-mesh, so
the suite passes on a fresh checkout):

```bash
pytest -q test_mesh_quality.py
```

### Reproducibility — accepted mesh set

The six target meshes in `results/msh/` were produced from the
canonical commands in **§Workflow / Canonical accepted-set commands**
above.  Every quality knob lives in `safs_fault_box_nwcut.geo`'s
`If (!Exists()) ... EndIf` guards or in `run_nwcut_meshing.py`'s
argparse defaults, so the two driver invocations are sufficient — no
extra `-setnumber` overrides needed.  Reruns reproduce these meshes
byte-equivalent up to gmsh's parallel-mesher non-determinism.

Defaults applied (do not override unless you know why):

| Knob | Source | Default | Role |
|---|---|---|---|
| `--fault-top-clamp` | driver | **100 m** | gmsh-headroom budget; fault top at z = -100 |
| `--pad-top` | driver | `= --fault-top-clamp` | mesh top lands at z = 0 |
| `--lc-floor` | driver | **100 m** | floor on local mesh size (`lc_min` in .geo) |
| `lc_near` | .geo / `--lc-near` | 1500 m | fault-corridor target tet edge |
| `lc_far` | .geo / `--lc-far` | 10 000 m | far-field cap (zgraded family) |
| `dist_inner` / `dist_outer` | .geo | 3 km / 40 km | fault-distance threshold-field band |
| `z_lc_top` / `_basin` / `_trans` / `_deep` | .geo | 500 / 1000 / 2000 / 3000 m | zgraded depth-step caps |
| `do_optimize` | .geo | **1** | Mesh.Optimize on |
| `do_optimize_netgen` | .geo | **0** | Mesh.OptimizeNetgen off (SIGBUSes on this geometry; see plan Phase 3 Risk Assessment) |
| `optimize_threshold` | .geo | 0.3 | gmsh's gamma threshold for optimize |
| `smoothing_passes` | .geo | 10 | Mesh.Smoothing iterations |
| `Mesh.Algorithm3D` | .geo | 1 (Delaunay) | HXT = 10 reports "No closed volume" |
| `Mesh.QualityType` | .geo | 2 (Inradius/Circumradius = gamma) | matches optimize_threshold semantics |

Accepted mesh quality snapshot (Joe-Liu η; see
`docs/PLAN_mesh_quality.md` for details):

| Mesh | n_tet | tet_emin [m] | η_min | slivers (η<0.1) | Q1 ≥100 m | Q2 η>0.1 |
|---|---:|---:|---:|---:|:-:|:-:|
| **2000m_zgraded** ⭐ | 1,257,877 | 112.72 | 0.0985 | 1 | ✓ | ✗ (1 tet at η ≈ 0.0985) |
| 1000m_lcfar3000 | 963,481 | 112.18 | 0.0877 | 2 | ✓ | ✗ |
| 2000m_lcfar3000 | 917,633 | 115.55 | 0.0518 | 1 | ✓ | ✗ |
| 1000m_zgraded | 1,301,846 | 105.13 | 0.0324 | 1 | ✓ | ✗ |
| 500m_zgraded | 1,446,645 | 38.72 | 0.0693 | 17 | ✗ | ✗ |
| 500m_lcfar3000 | 1,109,906 | 38.72 | 0.0683 | 21 | ✗ | ✗ |

⭐ = best mesh (highest η_min among Q1-passing variants).

**Known residual.** No variant fully satisfies Q2 (η > 0.1) — gmsh's
only sliver-killing optimizer (`Mesh.OptimizeNetgen`) crashes with
SIGBUS on the SAFS embedded-fault geometry, so it is disabled by
default.  `Mesh.Optimize` (Laplacian) + 10 smoothing passes resolves
all but 1–2 tets per variant; the survivors sit on the fault-edge /
basin-layer intersection where smoothing alone cannot do the
edge-swap required.  The 500 m variants additionally miss the Q1
100 m floor because the cleaned 500 m STL has ~112 m fault triangle
edges that propagate into the bulk despite the `lc_min` floor.

The accepted state was discussed and approved on **2026-05-17** —
see `docs/PLAN_mesh_quality.md` for the decision trail.
