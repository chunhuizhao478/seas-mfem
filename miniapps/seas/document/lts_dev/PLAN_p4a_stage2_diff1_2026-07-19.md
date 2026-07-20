# Phase 4a Stage 2 — cross-cluster (diff-1) rank seams: implementation plan — 2026-07-19

Branch `safs-v4_0_0-alt-case1-mfem-speed` (local). Builds on Stage 1
(`REVIEW_p4a_stage1_impl_2026-07-19.md`).  Design fork recorded in
`PLAN_p4a_mpi_exchange_2026-07-19.md`.

## The problem
A diff-1 rank seam joins a FINE local element (cluster c) on rank A to a COARSE
element (cluster c+1) on rank B (higher id = coarser = larger dt).  The coarse
element's in-tray (accumulate buffer) is on rank B, so rank A cannot deposit into
it.  Solution (plan req 3, "no messages needed" beyond the ghost field): each rank
fills its OWN element's contribution from exchanged ghost data, using the
closed-form sub-interval `[a,b]` (rank-identical).

## Mechanism (per diff-1 seam face, over one coarse step [T, T+dt_c])
The fine corrects twice (sub-intervals i=1,2). At fine correct i:
- **Rank A (fine, local):** forecast = `IntegrateTaylor(a_i,b_i, D_coarse_ghost)`;
  `F_i = InteriorFaceFlux(fine I_i, forecast)`; scatter fine's share → fine rhs.
- **Rank B (coarse, local):** forecast = `IntegrateTaylor(a_i,b_i, D_coarse_local)`;
  `F_i = InteriorFaceFlux(fine I_i ghost, forecast)`; scatter coarse's share →
  coarse element's LOCAL seam buffer.
At the coarse correct, rank B consumes its local seam buffer.  Both ranks feed
`InteriorFaceFlux_` the SAME (fine I_i, coarse forecast over [a_i,b_i]) → the flux
is bit-identical → conservation exact by single-evaluation.

`[a,b]` uses the COARSE side's expansion point (cluster c+1):
`t_origin = t_s + dt_base·2^(c+1)·floor(tick/2^(c+1))`,
`a = t(tick) − t_origin`, `b = min(t(tick)+dt_fine, t_s+T_actual) − t_origin`.

## Exchanges (matched collectives, P-007)
- **Fine I (per fine correct):** already done by Stage-1
  `ComputeADERClusterSeamFaceFluxRHS` (exchanges the correcting cluster's I).  The
  coarse rank reads the fine's ghost I from it.
- **Coarse D(k) (per coarse predict) — NEW:** when a cluster that is a GLOBAL
  cross-rank provider predicts, exchange its provider elements' Taylor stack to
  ghosts.  4a: scalar `ghost_gf_` per component ⇒ `NUM_STATE·order` exchanges per
  such predict (4b batches to `order`).  Cache `ghost_dk_[k][c]` until the next
  predict overwrites it (valid across the whole coarse step; GAP-A3 epoch).
- `n_collectives += (NUM_STATE·order) · #{predict clusters with global cross-rank
  provider faces}` — a new `LtsGlobalMeta` per-cluster flag, reduced across ranks.

## Data (operator members, all mutable, np>1 only)
- `ghost_dk_` : `std::vector<Vector>` size `order·NUM_STATE` — cached ghost Taylor
  stacks (component-major, `[k*NUM_STATE + comp]`, indexed `nbr_idx*ndof_per_el`).
- `seam_coarse_buf_` : `std::map<int,std::vector<real_t>>` (local coarse-seam elem
  → block) + `seam_coarse_fill_` : `std::map<int,int>` — the coarse-side seam
  in-tray (fill/consume/zero-at-sync analogue of LtsAccumulateBuffers).

## Methods
1. `ExchangeClusterProviderDkGhost(cluster_c, dk_data, dk_order, provider_slot_of_elem)`
   — pack cluster-c provider D(k) into `ghost_gf_`, exchange per (level,component),
   cache `ghost_dk_`.  Called from `LtsBulkSyncStepper::Predict(c)` when
   `meta.global_cross_rank_provider[c]`.
2. Extend `ComputeADERClusterSeamFaceFluxRHS(cluster_c, I, dt, rhs, t_s, dt_base,
   tick, T_actual)` — after the diff-0 case, for each shared face with local Elem1:
   - **fine-side** (local∈c, ghost∈c+1): coarse forecast from `ghost_dk_` over
     `[a,b]`; fine's share → rhs.
   - **coarse-side** (local∈c+1, ghost∈c, i.e. iterate faces whose GHOST is the
     just-corrected cluster c): coarse forecast from LOCAL D(k) over `[a,b]`;
     coarse's share → `seam_coarse_buf_[local]`.  Needs the local coarse D(k) →
     pass `dk_data`/`provider_slot_of_elem` to the seam corrector.
   - remove the diff-1 `MFEM_ABORT`.
3. Extend `AdvanceADERClusterBulk` — at cluster-c correct, CONSUME
   `seam_coarse_buf_` for cluster-c coarse-seam elements into rhs before the mass
   inverse (mirrors the local buffer consume); assert zero at sync.

## Test (gate)
`test_lts_mpi_seam_diff1` (or extend the existing): MakeCartesian3D split so
cluster ids AND the ParMesh partition BOTH split at x=mid ⇒ a GUARANTEED diff-1
rank seam at np=2 (Nc=2).  np=2==np=1 to 1e-10 over K syncs (centroid-weighted
reductions); ghost-poison; per-tick counter == n_collectives; seam-off sensitivity.

## Invariants to assert (correctness gates)
- `seam_coarse_buf_` fill count == #fine sub-steps of the closing coarse step.
- `seam_coarse_buf_` all-zero at every sync point.
- ghost D(k) read at the CURRENT coarse epoch (no stale read across a coarse step).
- fine-side and coarse-side use the SAME `[a,b]` (closed-form; interval-mismatch
  audit ≤ 1e-12·dt_c).

## Staging within Stage 2
- **2a:** D(k) ghost exchange primitive + `ghost_dk_` cache + a unit test
  (ghost D(k) == neighbour's local D(k)).  Self-contained.
- **2b:** meta global-provider flag + tick-table n_collectives extension
  (test_lts_scheduler).  Self-contained.
- **2c:** the seam corrector diff-1 paths + coarse-side seam buffer + consume;
  remove fail-loud.  Gated by the diff-1 np=2==np=1 test.
