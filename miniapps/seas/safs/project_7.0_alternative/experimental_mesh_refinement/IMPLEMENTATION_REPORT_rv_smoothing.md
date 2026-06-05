# Implementation Report — reduce fault volume-ratio r_v by smoothing

**Plan:** `experimental_mesh_refinement/PLAN_rv_volume_smoothing.md`
**Date:** 2026-06-04
**Target mesh:** `safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh`
**Metric:** Eq. (18), Zhang et al. (2023) — `r_v = max{V_A/V_B, V_B/V_A}` per
fault triangle (A/B = the two tets sharing it).
**Method (user-chosen, 2026-06-04):** fully-unstructured near-fault
volume-balancing smoothing with in-plane boundary sliding (`--slide-boundary`);
**no re-mesh, no structured layer, topology unchanged.**

## Result — TARGET MET (revised A1)

| Metric | Baseline `_triq` | Smoothed `_rvsmooth` | Target | Pass |
|---|---:|---:|---:|:-:|
| r_v mean | 1.2339 | **1.0339** | < baseline | ✓ |
| r_v median | 1.1553 | **1.0188** | — | — |
| r_v p95 | 1.6379 | **1.1100** | — | — |
| r_v p99 | 2.5000 | **1.2084** | < 1.3 | ✓ |
| fault faces r_v > 1.5 | 5,746 (9.47 %) | **72 (0.12 %)** | < 0.5 % | ✓ |
| fault faces r_v > 2.0 | 1,121 (1.85 %) | **19 (0.03 %)** | — | — |
| **r_v max** | 10.7546 | **8.8855** | irreducible¹ | n/a |
| Joe-Liu eta_min (Q2) | 0.2032 | **0.2032** | ≥ baseline, > 0.1 | ✓ |
| eta_med | 0.8346 | 0.8332 | — | — |
| tet_emin (Q1) [m] | 113.4 | **122.0** | ≥ 100 | ✓ |
| inverted tets | 0 | **0** | 0 | ✓ |
| fault faces 2-tet embedded | 60,658 | **60,658** | all | ✓ |
| mesh_zmax [m] | 0.0 | **0.0** | 0 exactly | ✓ |
| fault nodes moved | — | **0** | 0 | ✓ |
| nodes moved | — | 12,297 (5.9 %) | — | — |

`pytest -q test_rv_volume_smoothing.py` → **7 passed** (gate clears; baseline
fails A1, so the gate has teeth).

¹ **The max r_v (8.89) is irreducible by node smoothing** and is *accepted +
documented* per the user's decision — see Known limitations.

## What was implemented

1. **`smooth_fault_volume_ratio.py`** — the tool (numpy only).
   - Gmsh v2.2 reader/writer that preserves node ids and the entire element
     block verbatim, rewriting only moved node coordinates (so connectivity,
     physical tags, and unmoved nodes are byte-identical).
   - `build_fault_apex_map` — pairs each fault triangle with its two owning
     tets' apex nodes; **raises** if any fault face is not 2-tet embedded.
   - Key insight: A and B share the *same* fault triangle, so `r_v = h_A/h_B`,
     the ratio of apex perpendicular distances to the fixed fault plane —
     reducing r_v ⇔ equalizing the two apex heights.
   - `build_free_axes` — fault nodes fully fixed; boundary nodes locked on the
     axis normal to each box face they lie on (so they slide only *in-plane*;
     box edges keep one axis, corners none); interior nodes fully free.
   - `smooth` — vectorized Jacobi sweeps: per fault face, pull the taller apex
     toward / push the shorter apex away from the fault to the mean height,
     scatter-averaged per node, damped by `omega`, projected onto free axes.
     **Validity guard** reverts any node whose move would invert a tet, drop its
     Joe-Liu `eta` below `min(eta0, eta_floor)`, or drop its shortest edge below
     `min(edge0, edge_floor)` — floors default to the input global minima, so
     neither Q1 nor Q2 can regress. Convergence on 99th-pct r_v.
2. **`test_rv_volume_smoothing.py`** — pytest gate (A1–A8 + teeth), skip-on-missing.
3. **`compute_eq18_volume_ratio.py`**, **`export_eq18_fault_vtu.py`** — r_v
   reporter and per-face VTU writer (shared with the TPV102 r_v work).
4. **Outputs:**
   - `safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_rvsmooth.msh` (deliverable)
   - `safs_baseline_fault_rv.vtu`, `safs_rvsmooth_fault_rv.vtu` (per-face r_v)
   - `eq18_rv_safs_before_after.png` (fault-plane before/after)

## Hard r_v ≤ 1.5 ceiling (`--rv-ceiling`) — added 2026-06-04

Zhang et al. place the SSO threshold at r_v > 1.5, so for a **benchmark** mesh we
enforce a hard guard: no fault face may exceed 1.5. `smooth_fault_volume_ratio.py
--rv-ceiling 1.5`:

- keeps iterating until **0** faces exceed the ceiling;
- when it stalls with faces still over, **adaptively relaxes** the eta/edge
  quality floors by `relax_factor` (0.85) — but **never below the Q1/Q2 project
  gates** (`eta_gate=0.105`, `edge_gate=100 m`) — trading minimum element quality
  (down to, not below, the gates) for the ceiling;
- returns the over-ceiling count; **`main()` ERRORS (exit 3) and writes nothing**
  if the ceiling cannot be met, so a violating benchmark mesh is never produced.

**TPV102 `tpv102_200m.msh` (the benchmark) — ceiling MET:**

| | mean | p99 | max | >1.5 | eta_min | tet_emin | Q1/Q2 |
|---|---:|---:|---:|---:|---:|---:|:-:|
| baseline | 1.209 | 1.750 | 2.195 | 7.68 % | 0.3552 | 131.6 | ✓ |
| `_rvsmooth` (`--rv-ceiling 1.5`) | 1.029 | 1.160 | **1.419** | **0** | 0.2567 | 106.6 | ✓ |

The adaptive guard relaxed eta only 0.355→0.302→0.257 (two steps) before clearing
the last faces; eta_median actually *rose* 0.835→0.856. Fault frozen, box
preserved, 0 inverted tets, topology byte-identical, mesh_zmax=0.
Gate: `pytest -q test_rv_ceiling_tpv102.py` → 5 passed.

**SAFS — ceiling correctly REFUSED.** With adaptive relaxation all the way to the
Q1/Q2 gates the over-count fell 70→23 but stalled; the tool errored
(`23 face(s) remain above 1.5, max 4.94 … re-triangulation required. Nothing
written`). Those 23 are the fault-spanning sliver tets — confirming the hard
guard cannot be met by node smoothing on the non-planar fault (it needs the
re-triangulation follow-up), and the guard refuses rather than ship a bad mesh.

This is why the SAFS deliverable above does **not** use `--rv-ceiling` (it would
fail), whereas the planar TPV102 benchmark does.

## Targeted Gauss-Seidel ceiling enforcement (`targeted_ceiling_pass`) — 2026-06-04

The global Jacobi sweeps can stall with faces still over the ceiling because an
apex shared by many faces gets an averaged (cancelling) displacement. Added a
**targeted Gauss-Seidel pass** that runs after the global sweeps when faces
remain over: worst-first, it moves each over-ceiling face's movable apex the
**minimal** distance (to the band edge `ceiling·margin`, damped) that an
incident-tet line search (sign + eta + edge guards, at the Q1/Q2 gates) allows,
updating immediately. Backed by a CSR node→incident-tet adjacency
(`build_node_tets_csr`). Returns the best (fewest-over) state visited. `--best-effort`
writes the best mesh with a warning when the ceiling still can't be met.

This cleanly separates the two failure modes:
- **dilution-blocked** faces (Jacobi averaging cancels the needed move): the GS
  pass clears them. Validated on **TPV102 with ceiling 1.1**: global sweeps
  stalled at 789 faces > 1.1; the GS pass drove it to **110** (max 1.33→1.21),
  Q1/Q2 still pass, 0 inverted, fault frozen.
- **inversion-blocked** faces (no valid step exists — balancing would flip a
  back tet): node motion cannot fix them; the GS pass reports them unchanged.

## TPV31 `tpv31_50m.msh` (6.3 M tets) — ceiling 1.5 NOT reachable by smoothing

| | mean | p99 | max | >1.5 | >2.0 | eta_min | tet_emin |
|---|---:|---:|---:|---:|---:|---:|---:|
| baseline | 1.238 | 1.820 | 2.292 | 11.84 % (49,246) | 0.09 % | 0.3025 | 31.8 |
| `_rvsmooth` (best-effort) | 1.050 | 1.262 | **1.868** | **0.029 % (122)** | **0** | **0.3025** | 31.8 |

The global sweeps reduced 49,246 → **122** faces > 1.5; the targeted GS pass then
moved **none** of them (122→122→122→122) — they are **inversion-blocked**
(unlike TPV102's dilution-blocked faces). So `r_v ≤ 1.5` is **not reachable by
node smoothing** for this mesh; the 122 (max 1.868, depths −350 … −14,899 m, i.e.
along the surface trace and the bottom edge where the near-fault wedge is too
tight to move an apex) require **local fault re-triangulation**.

Quality is **fully preserved** (eta_min 0.3025 unchanged — the 122 are
inversion/edge-blocked, not eta-blocked, so no relaxation room was spent;
eta_median 0.841→0.837; tet_emin 31.8 unchanged). Fault frozen, box preserved,
0 inverted, topology byte-identical, mesh_zmax=0. `Q1=False` is a **pre-existing**
property of the input (31.8 m min edge < the SAFS-derived 100 m floor, which is
inappropriate for a 50 m mesh); smoothing preserved it exactly.

**Deliverable:** `tpv31_50m_rvsmooth.msh` (best-effort, max 1.868, 0.029 % > 1.5).
A strict `r_v ≤ 1.5` benchmark mesh needs the re-triangulation follow-up.

## Local 1→4 split re-triangulation (`retriangulate_fault_rv.py`) — 2026-06-04

For faces that are **inversion-blocked** (node smoothing cannot reach the ceiling
— e.g. TPV31's 122), a conforming **local topology fix** drives them under the
ceiling with a hard guarantee.  For each over-ceiling fault face, split the
**taller** of its two tets `A = (v0,v1,v2,a)` by inserting a node `p` on the
segment `centroid(f) → a` at the **short** side's height, and replace `A` by the
four tets `(v0,v1,v2,p)`, `(v0,v1,a,p)`, `(v1,v2,a,p)`, `(v0,v2,a,p)`.  The new
fault-adjacent tet `(v0,v1,v2,p)` then has the short side's height ⇒ **r_v = 1**.

Properties (why this is safe and conforming):
- `p` lies on the centroid→apex segment, so it is strictly interior to `A` ⇒ all
  4 sub-tets are non-degenerate (orientation matched to the parent).
- All four OUTER faces of `A` are preserved ⇒ neighbours are untouched, the fault
  triangle `f` stays shared by exactly 2 tets (conforming), and **no fault or
  boundary node moves**.
- On a **planar** fault a tet borders at most one fault face (two fault triangles
  are coplanar), so the tall tets of distinct over-ceiling faces are disjoint —
  the splits are independent and collision-free.
- Cost: +1 node, +3 tets per fixed face (negligible).  `--target-rv 1.0` (place
  `p` at the short height) maximises the side sub-tet size ⇒ best η.

Validated on **TPV102 (raw, 2,914 faces > 1.5)**: split → max r_v 1.4999, **0
faces > 1.5**, 0 inverted, fault frozen, all faces still 2-tet embedded, η_min
0.14 (> Q2 floor), η_med 0.857. 25 s.

## TPV31 — r_v < 1.5 ACHIEVED (smoothing → local split)

The required pipeline for the inversion-blocked benchmark:
1. `smooth_fault_volume_ratio.py --slide-boundary --rv-ceiling 1.5 --best-effort`
   → 49,246 → 122 faces > 1.5 by node motion (η preserved at 0.3025).
2. `retriangulate_fault_rv.py --rv-ceiling 1.5 --target-rv 1.0` on the result
   → splits the 122 tall tets → **0 faces > 1.5**.

| stage | max r_v | >1.5 | η_min | tets |
|---|---:|---:|---:|---:|
| baseline | 2.292 | 11.84 % | 0.3025 | 6,304,859 |
| + smoothing | 1.868 | 122 | 0.3025 | 6,304,859 |
| **+ split (final)** | **1.4996** | **0** | **0.1836** | 6,305,225 |

**Final deliverable: `tpv31_50m_rvsmooth_split.msh` — every fault face r_v < 1.5
(max 1.4996), 0 non-positive tets, fault nodes frozen, conforming (all 416,060
fault faces 2-tet embedded), Gmsh v2.2.**  Cost: +122 nodes / +366 tets, and
η_min 0.3025 → 0.1836 (the split side sub-tets; still well above the 0.1 Q2
floor; η_median unchanged at 0.837).  This is the deliberate trade you required:
guarantee r_v < 1.5 at a modest η_min cost.

## Decisions / deviations from the original plan

- **Boundary sliding enabled (`--slide-boundary`)** — original plan invariant 2
  froze the boundary entirely. The residual r_v>1.5 tail was 96 % boundary-
  touching tets; the user chose to let boundary nodes slide *in-plane* (box
  geometry preserved exactly; only in-plane positions change), cutting the tail
  from 1.66 % to 0.12 %. Interior-only (boundary byte-frozen) is retained as the
  default mode (flag off) and also passes all gates (1.66 % > 1.5, mean 1.060).
- **A1 re-scoped** from `max < 2.0` to `r_v>1.5 fraction < 0.5 %` and `p99 < 1.3`
  — `max < 2` is structurally unreachable (see below); the user accepted this.
- **Guard is eta+edge based, not volume-fraction** (plan said `floor_frac` on
  volume). A pure volume floor let `eta_min` regress 0.203→0.117 in testing; the
  eta+edge guard keyed to the input global minima preserves both Q1 and Q2.
- **Convergence on p99, not max** — the max is pinned by the irreducible tets and
  would early-stop the bulk relaxation.

## Known limitations / not done

- **Irreducible max r_v ≈ 8.89 (~12 fault-spanning sliver tets).** These tets
  have 3 nodes on one fault triangle and their apex *on the fault* (a nearby,
  normal-offset fault triangle). The apex is a fault node we must not move, so
  no interior-node smoothing can balance them. Counts: 12 of 60,658 faces > 1.5
  carry a fault-node apex; the worst (8.89) is at depth z = −1816 m. Fixing them
  requires **local fault re-triangulation / edge-flips** (the CGAL `_triq`
  track), which changes the fault surface — out of scope here per the user
  decision (accept + document).
- **Scoped to the named 500 m lcfar3000 mesh.** The tool is parameter-free on
  geometry (reads tags 1/101/102/103/104) and generalizes, but other variants
  were not processed.
- **Downstream stress/velocity projection and the SAFS simulation were not
  re-run** (separate follow-up once the smoothed mesh is accepted).

## Files

- [created] `PLAN_rv_volume_smoothing.md`
- [created] `smooth_fault_volume_ratio.py`
- [created] `test_rv_volume_smoothing.py`
- [created] `compute_eq18_volume_ratio.py`, `export_eq18_fault_vtu.py`
- [created] `safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_rvsmooth.msh`
- [created] `safs_baseline_fault_rv.vtu`, `safs_rvsmooth_fault_rv.vtu`,
  `eq18_rv_safs_before_after.png`
- [modified] **none** (input mesh + existing code read-only; `git status` shows
  only untracked additions)

## Reproduce

```bash
conda activate pythonenv          # numpy, meshio
cd experimental_mesh_refinement

# 1. smooth (slide-boundary deliverable)
python smooth_fault_volume_ratio.py \
    --in  safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh \
    --out safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_rvsmooth.msh \
    --omega 0.6 --slide-boundary --max-sweeps 60
#   (drop --slide-boundary for the boundary-frozen interior-only variant)

# 2. gate + quality
python -m pytest -q test_rv_volume_smoothing.py
python ../meshing/code/check_mesh_quality.py \
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_rvsmooth.msh

# 3. r_v VTUs (open in ParaView, color by 'rv' / 'log10_rv')
python export_eq18_fault_vtu.py <mesh.msh> 101 1 <out.vtu>
```
