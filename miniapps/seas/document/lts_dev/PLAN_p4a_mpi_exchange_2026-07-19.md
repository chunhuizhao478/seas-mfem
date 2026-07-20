# Phase 4a — MPI ghost-field exchange for clustered LTS (bulk): design + staging — 2026-07-19

Branch `safs-v4_0_0-alt-case1-mfem-speed` (local, unpushed). Flow: design →
implement → review → fix → verify locally (np=2==np=1 differential smoke).
Grounds §Phase 4 / requirement 1 + 3 of `PLAN_clustered_lts_ader_2026-07-18.md`.

## What is already in place (Step 0, np>1 WIRED)
- Serial-mesh clustering on rank 0 + MPI_Bcast of `serial_cluster`, Nc, λ,
  dt_base, and the LTS-aware partition `part_data` (all ranks) —
  `spatial_dyn_driver.cpp:1560-1704`.
- Serial→local cluster-id map (`lts_cluster_id`) built from `part_data` —
  `:1984-2005`.
- `BuildLtsLayout` + `LtsBulkSyncStepper` + `RunSyncInterval` tick loop; the
  bulk stepper is byte-exact-off and Nc=1==GTS at np=1.
- `ComputeADERSharedFaceFluxRHS` — the GTS seam-flux primitive: NUM_STATE
  per-component `q_gf.ExchangeFaceNbrData()` of the time-integral `I`, reads the
  ghost neighbour `I` via `nbr_idx = Elem2No - ne_`, and scatters ONLY the local
  side (`wave_operator.inl:5865-6260`). This is the reuse target.

## What blocks np>1 today (the 4a gap)
1. `BuildLtsMeshInputsFromMaterial` collapses every rank-seam face into
   `elem2 < 0` → `FaceRole::Boundary` (`spatial_dyn_driver.cpp:246`,
   `lts_layout.cpp:93-98`); the layout has no seam role and no ghost-neighbour
   cluster id.
2. The cluster corrector `AdvanceADERClusterBulk` SKIPS rank-seam faces
   (`wave_operator.inl:2232`) and never exchanges — so the seam flux is dropped.
3. `ReduceGlobalMeta(layout, nullptr)` (`:2001`) does no cross-rank reduce, so
   the global per-cluster counts that gate collectives are wrong at np>1.
4. Both LTS driver blocks `MFEM_VERIFY(nprocs == 1, ...)`.

## The seam, precisely (two cases)

A rank-seam bulk face joins a LOCAL element `e1` (cluster `c1`) to a GHOST
element (cluster `c2`, the face-neighbour). maxdiff≤1 ⇒ `|c1−c2| ∈ {0,1}`.

### Case A — same-cluster seam (`c1 == c2`, "diff-0")
Both sides step at the same rate. This is EXACTLY the GTS shared-face pattern
restricted to cluster `c`: when cluster `c` corrects, each rank exchanges its
cluster-`c` time-integral `I_c`, computes `F(I_self, I_nbr_ghost)` with its own
local orientation, and scatters ONLY its own side. Both ranks feed the SAME
two-sided data to `InteriorFaceFlux_` (its own `I` + the ghost's `I`), so the
flux is bit-identical and conservation is exact by single-evaluation.
**Reuse:** a cluster-restricted `ComputeADERSharedFaceFluxRHS`. Low risk.

### Case B — cross-cluster seam (`|c1−c2| == 1`, "diff-1")
The coarse element's in-tray (accumulate buffer) is the crux: in the local
design the FINE side deposits the coarse's share into the coarse element's
buffer. At a rank seam the coarse element (and its buffer) live on the OTHER
rank, so the fine rank cannot deposit into it. Plan requirement 3 ("no messages
needed") resolves this by making **each rank fill its OWN element's
contribution** from exchanged ghost data, using the closed-form sub-interval
`[a,b]` (rank-identical, so both agree without a message):

| Side (rank) | reads local | reads ghost | integrates | scatters into |
|---|---|---|---|---|
| fine `e1`,`c` (rank A) | fine `I` (∫ over its sub-step) | coarse `D(k)` | `IntegrateTaylor(a,b,D_coarse_ghost)` | fine `rhs` |
| coarse `e2`,`c+1` (rank B) | coarse `D(k)` | fine `I` (∫ over its sub-step) | `IntegrateTaylor(a,b,D_coarse_local)` | coarse LOCAL buffer |

Both evaluate `F(∫fine over [a,b], ∫coarse-forecast over [a,b])` from identical
data → bit-identical flux → conservation exact. The coarse rank accumulates its
share into its LOCAL buffer once per fine sub-correct (it knows the fine cluster
corrected from the rank-identical tick table), and consumes the buffer at its
own correct — no cross-rank buffer message.

**Extra payload for Case B:** the coarse's `D(k)` Taylor stack must reach the
fine rank (so the fine can integrate the forecast). This is a per-predicting-
cluster ghost exchange of the retained provider `D(k)`, NOT covered by the
current `n_collectives = num_state·|correct|` formula. The formula must gain a
provider-`D(k)` term for predicting clusters that have GLOBAL cross-rank
provider faces (mirrors the fault-predictor term already there, but for bulk).

## Collective schedule (Matched-collectives, P-007)
`RunSyncInterval` is a pure function of the rank-identical tick table, so every
rank runs `Predict(c)`/`Correct(c)` for the SAME clusters in the SAME order.
Putting the exchanges INSIDE `Predict`/`Correct` (a fixed count per call) makes
them matched by construction. Per tick:
- `Correct(c)` (bulk): `num_state` per-component `I` exchanges (the seam
  corrector) → `num_state` per correcting cluster. Case-A + Case-B fine-side +
  Case-B coarse-side all consume the SAME single `I` ghost exchange.
- `Predict(c)` (bulk, Case B only): one batched `D(k)` ghost exchange if cluster
  `c` has global cross-rank provider faces → `+1` (or `+order`) per such
  predicting cluster.
- Fault predictor exchange: dropped under D-2 (fault rank-interior).
- Per-sync (not per-tick): one V_max Allreduce + one NaN Allreduce (already in
  the driver sync block; identical on all ranks).

A per-tick debug counter (rank-local increment inside the exchange calls) must
equal `tk.n_collectives` on EVERY rank at the end of each tick (B.11).

## Ghost staleness tripwire (P-007)
Before a cluster consumes ghost data, NaN-poison the ghost slots of elements NOT
in a due cluster this tick, so a stale/wrong-tick read trips a NaN immediately
(debug builds). Implemented as an optional debug pass keyed off an env flag,
off by default (byte-exact-off).

## Staging (each stage: build + local gate before the next)

**Stage 1 — foundation + Case A (same-cluster) seams.** ReduceGlobalMeta MPI
hook; seam classification (ghost cluster id + ghost face-nbr index) in the mesh-
inputs builder; `LtsFaceSpec`/`BuildLtsLayout` seam roles; cluster corrector
processes its cluster's diff-0 seam faces via the restricted shared-face
exchange; per-tick collective counter; ghost NaN-poison; remove the bulk
`nprocs==1` guard. **FAIL LOUD on any diff-1 rank seam** (so the code is never
silently wrong before Stage 2). Gate: np=2==np=1 to 10 digits on a mesh whose
partition gives only diff-0 seams (e.g. single-cluster, or a cluster-coherent
partition), + `lts="off"` byte-exact, + the per-tick counter matches.

**Stage 2 — Case B (cross-cluster) seams.** Provider `D(k)` ghost exchange;
fine-side seam consumer (ghost `D(k)` + closed-form `[a,b]`); coarse-side seam
consumer (local buffer fill from fine ghost `I`); extend `n_collectives`;
interval-mismatch audit. Gate: B.11-style np=2==np=1 on a 3-cluster chain
crossing the rank seam; ghost-poison green.

**Out of 4a (later):** 4b EDGE flux-premultiplied payloads (Phase 5 measures
4b); Expanse physics/perf validation; the fault half at np>1 (D-2 makes fault
rank-interior, so the fault stepper's `nprocs==1` guard is a separate follow-up
once the bulk seam path is validated).

## Design fork recorded for review
The Case-B payload is the one genuine core-numerics fork the master plan leaves
open ("no messages needed" + an `I`-only `n_collectives` formula do not by
themselves fix the coarse's-buffer-is-remote problem). This doc commits to:
**exchange the coarse provider `D(k)` to ghosts on predict; each rank fills its
OWN element's seam contribution; the coarse fills its LOCAL buffer from the fine
ghost `I`; extend `n_collectives` with a bulk-provider-`D(k)` term.** The
adversarial /code-review must scrutinize this against conservation + the
matched-collective contract.
