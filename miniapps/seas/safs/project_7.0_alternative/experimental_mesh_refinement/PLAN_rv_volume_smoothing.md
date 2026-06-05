# PLAN — reduce fault volume-ratio r_v by near-fault volume-balancing smoothing

**Date:** 2026-06-04
**Target mesh:** `experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh`
**Metric:** Eq. (18) of Zhang et al. (2023), JGR Solid Earth 10.1029/2022JB025817:
`r_v = max{V_A/V_B, V_B/V_A}` per fault triangle, A/B the two tets sharing it.
**Baseline (measured):** n_fault_tri = 60,658; r_v mean 1.234, median 1.155,
max **10.755**; 9.47 % of fault faces have r_v > 1.5, 1.85 % > 2.0.

## Goal

Lower the high r_v tail on this **non-planar** San-Andreas fault while keeping
the mesh a **fully unstructured Delaunay tetrahedral mesh** — no structured
boundary layer, no re-mesh, no topology change. (Method chosen by the user over
the symmetric-normal-layer alternative, 2026-06-04.)

## Key geometric fact (makes the method exact-targeting)

The two tets A and B that share a fault triangle share the **same** triangle, so
their base areas are identical. Hence

    V = (1/3)·area(face)·h   ⇒   r_v = max(h_A/h_B, h_B/h_A),

where h_A, h_B are the perpendicular distances of the two **apex** nodes (the
4th node of each tet, off the fault) from the fault triangle's plane. The fault
triangle is fixed (its 3 nodes are fault-surface nodes we never move), so its
plane is constant. Reducing r_v ⇔ equalizing the two apex heights h_A, h_B.
The only free variables are the near-fault **interior** apex nodes.

## DELIVERED MODE (revised 2026-06-04 after the user decision)

The user chose **`--slide-boundary`** (over interior-only) for the named
deliverable, and to **accept + document** the irreducible max. So the delivered
method additionally lets domain-boundary nodes slide *in-plane* on their box
face (invariant 2 below is relaxed accordingly). Interior-only remains the
default (boundary frozen) and is retained behind the off-by-default flag.

## Hard invariants (must hold after smoothing — from PLAN_mesh_quality.md / the _triq report)

1. **Fault geometry frozen.** Every fault-surface node keeps its exact
   coordinates ⇒ fault triangulation, areas, normals, and the z=0 trace are
   byte-identical; `mesh_zmax == 0.0` exactly; fault↔input Hausdorff unchanged.
2. **Domain box preserved.** Every node on Physical Surfaces 102 (top), 103
   (bottom), 104 (sides) stays exactly on its axis-aligned box face: the
   coordinate axis normal to each face it lies on is byte-unchanged (top nodes
   keep z=0, side nodes keep x/y=const, box edges keep their two locked axes,
   corners are fully fixed). In `--slide-boundary` mode the in-plane position
   may change; in the default mode boundary nodes are byte-identical. Either way
   the domain box geometry is identical.
3. **Topology frozen.** Same node count, same element count, same connectivity,
   same physical tags (rock=1, fault=101, top=102, bottom=103, sides=104). The
   fault stays embedded: all 60,658 fault triangles remain faces of exactly two
   bulk tets; the fault surface stays a single connected component.
4. **Validity.** No inverted tets (every tet keeps its original sign of signed
   volume) and no new slivers: bulk Q1 (`tet_emin ≥ 100 m`) and Q2
   (`eta_min > 0.1`) gates of `check_mesh_quality.py` must still pass, and
   `eta_min` must not regress below the baseline value.
5. **Gmsh v2.2 ASCII output** (MFEM reader requirement).
6. **No existing file modified.** New artefacts only, all under
   `experimental_mesh_refinement/`. The input `_triq` mesh is read-only.

## Algorithm (vectorized Jacobi volume-balancing smoothing)

Movable set `M` = interior nodes that are **not** fault nodes and **not**
boundary nodes (union of all nodes used by Physical Surfaces 101/102/103/104).

Precompute once (topology is constant):
- Fault faces: for each fault triangle, the two owning tets and their apex
  nodes `(apex_a, apex_b)`; face plane unit normal `n_f` and centroid `c_f`
  (constant — depends only on fault nodes).
- Per-tet original signed volume `V0` (sign + magnitude) for the validity guard.
- For each tet, its node list (for the post-move validity sweep).

Each sweep (Jacobi — all displacements computed from the current positions,
then applied together):
1. h_a = (X[apex_a] − c_f)·n_f ; h_b = (X[apex_b] − c_f)·n_f (signed; opposite
   signs since A,B are on opposite sides).
2. Per face, target apex height magnitude:
   - both apexes movable → `t_f = (|h_a| + |h_b|)/2`
   - only one movable    → `t_f = |h_of the fixed apex|` (match the fixed side)
   - neither movable      → face contributes nothing
3. Per-face desired displacement of a movable apex p along the face normal:
   `d = sign(h_p)·(t_f − |h_p|)·n_f` (pull the taller apex toward the fault,
   push the shorter one away).
4. Scatter-add each face's `d` to its movable apex/apexes; average per node over
   the faces it participates in; damp by relaxation `omega` (default 0.5).
5. Tentative `X_new[M] = X[M] + disp[M]`.
6. **Validity guard (revert-bad-nodes):** a tet is *bad* at `X_new` if it
   flipped signed-volume sign (inversion), OR its Joe-Liu `eta` dropped below
   `min(eta0_tet, eta_floor)`, OR its shortest edge dropped below
   `min(edge0_tet, edge_floor)`. `eta_floor`/`edge_floor` default to the input
   mesh's global `eta_min` / shortest edge, so neither Q2 (eta) nor Q1 (min
   edge) global minimum can regress. Revert every movable node belonging to a
   bad tet to its pre-sweep position; recompute; repeat up to `max_revert`
   (default 6) passes (reverting only relaxes constraints, so it converges).
7. **Free-axis projection:** each node's displacement is masked by its free
   axes (`disp *= free_axes`), so fault nodes never move and boundary nodes move
   only within their box face.
8. Commit `X = X_new`.

Stop after `max_sweeps` (default 60) or when the **99th-percentile** r_v
improvement between sweeps `< tol` (default 1e-3) for 3 consecutive sweeps. (p99
not max, because the max is pinned by the irreducible fault-spanning tets below
and would early-stop the bulk relaxation.)

**Hard r_v ceiling (`--rv-ceiling C`, added 2026-06-04).** When set, the stop
rule changes: iterate until **0** faces exceed `C`; if the over-ceiling count
stalls while faces remain over, **adaptively relax** the eta/edge floors by
`relax_factor` (0.85), clamped at the Q1/Q2 gates (`eta_gate=0.105`,
`edge_gate=100 m`). `smooth` returns the best (fewest-over, then lowest-p99)
state and its over-ceiling count; `main` **errors (exit 3) and writes nothing**
if that count is > 0. Used for the TPV102 SCEC benchmark (planar fault → ceiling
met, max 1.42); SAFS cannot meet it (fault-spanning tets) so the guard refuses.

Rationale: balancing raises the *flatter* (small-h, near-sliver) near-fault tet
and lowers the taller one, so it generally *improves* worst-case near-fault tet
shape; the guard blocks any move that would degrade a neighbour tet behind the
moved apex.

## Files to create (all in `experimental_mesh_refinement/`)

1. **`smooth_fault_volume_ratio.py`** — the tool. Public functions:
   - `read_gmsh22(path) -> (coords: np.ndarray[N,3], blocks: dict)` where
     `blocks` holds tetra connectivity (phys 1) and triangle blocks by physical
     tag (101/102/103/104), preserving original element order and ids.
   - `build_fault_apex_map(coords, tets, fault_tris) -> FaultTopo` with
     per-face `(apex_a, apex_b, n_f, c_f)` and the 2-tet embedding check.
   - `rv_stats(coords, fault_topo) -> dict` (n, min, mean, median, p95, p99,
     max, n_gt15, n_gt20) — same definition as `compute_eq18_volume_ratio.py`.
   - `smooth(coords, tets, fault_topo, movable_mask, *, omega, floor_frac,
     max_sweeps, tol, max_revert, verbose) -> (coords_new, history)`.
   - `write_gmsh22(path, coords, blocks)` — re-emits the mesh with new
     coordinates, identical connectivity/tags, `MeshFormat 2.2`, PhysicalNames.
   - `main(argv)` CLI: `--in --out [--omega --floor-frac --max-sweeps --tol]`;
     prints before/after r_v + Joe-Liu eta (via `check_mesh_quality`).
2. **`test_rv_volume_smoothing.py`** — pytest gate (skips if mesh absent):
   - output max r_v < `RV_MAX_TARGET` (2.0) and % > 1.5 < baseline;
   - fault + boundary nodes byte-identical to input (frozen invariant 1–2);
   - node/element counts and physical tags identical (invariant 3);
   - all fault faces still 2-tet embedded; no inverted tets (invariant 4);
   - `mesh_zmax == 0.0`;
   - **teeth:** baseline mesh has max r_v > 2.0 (so the pass is meaningful).
3. **`IMPLEMENTATION_REPORT_rv_smoothing.md`** — results, decisions, repro.

## Acceptance criteria (A1 revised 2026-06-04: max<2 is structurally unreachable)

**Why A1 changed.** Every fault face with both apexes movable is driven below
2.0, but ~12 *fault-spanning* sliver tets have a fault-node apex (a tet with 3
nodes on one fault triangle and its 4th node on a nearby, normal-offset fault
triangle). Their apex is a fault node we must not move, so `r_v ≈ 8.9` there is
irreducible by node smoothing (only fault re-triangulation fixes it; the user
chose to accept + document this). A1 is therefore re-scoped to the controllable
distribution.

| # | Criterion | Check |
|---|---|---|
| A1 | fault faces with r_v > 1.5 reduced to **< 0.5 %** (from 9.47 %); p99 r_v **< 1.3** (from 2.50) | `rv_stats` / test |
| A2 | mean r_v reduced vs 1.234 baseline | `rv_stats` |
| A3 | fault nodes byte-identical; boundary nodes stay on their box face (locked axis byte-unchanged) | test |
| A4 | node count, element count, physical tags identical | test |
| A5 | all 60,658 fault faces still embedded (2 tets each); fault 1 component | test |
| A6 | 0 inverted tets; Q1 (`tet_emin≥100`) & Q2 (`eta_min>0.1`) still pass; `eta_min` ≥ baseline; `tet_emin` ≥ 100 | `check_mesh_quality` |
| A7 | `mesh_zmax == 0.0` exactly | test |
| A8 | output is Gmsh v2.2 ASCII, loads in meshio | test |
| A9 | no existing file modified | manual / git status |
| A10 | r_v VTU written (fault triangles colored by r_v) for before & after | manual |

## Out of scope

- Re-meshing or any structured boundary layer (that is the rejected alternative).
- Second-layer / global TMOP optimization (only first-layer apex nodes move).
- Re-running downstream velocity/stress projections or the SAFS simulation.
- Regenerating the 1000 m / 250 m / zgraded variants (tool is parameterized but
  only the named 500 m mesh is processed here).
