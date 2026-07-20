# P4a Stage 1 — MPI foundation + same-cluster (diff-0) rank seams: impl + review — 2026-07-19

Branch `safs-v4_0_0-alt-case1-mfem-speed` (local, unpushed). Flow: design →
implement → review (2 adversarial forks) → fix → verify locally. Design:
`PLAN_p4a_mpi_exchange_2026-07-19.md`.

## What Stage 1 is
The fault-free BULK multi-rate stepper now runs at np>1 by exchanging ghost `I`
across rank seams and applying the seam-face flux locally.  Scope = SAME-CLUSTER
(diff-0) seams — the P-022 "LTS(1 cluster) vs GTS given the LTS partition" gate.
A cross-cluster (diff-1) seam FAILS LOUD (Stage 2).  Nothing about the fault
half or Nc>1-at-np>1 physics is claimed here.

## Implementation
- **`ComputeADERClusterSeamFaceFluxRHS(cluster_c, I_cluster, dt, rhs)`**
  (wave_operator.inl): NUM_STATE per-component `q_gf.ExchangeFaceNbrData()` of the
  cluster's `I` (packs ONLY cluster-c blocks — UB-clean, R1), then for each shared
  face whose local Elem1 ∈ cluster_c: one-sided `SharedInteriorFaceFlux_` scatter —
  a line-for-line subset of the GTS `ComputeADERSharedFaceFluxRHS`.  Diff-1 seam
  (ghost cluster ≠ cluster_c) ⇒ `MFEM_ABORT`.  Shared FAULT face ⇒ abort (D-2).
- **`EnsureGhostClusterIds_()`** — one-time ghost cluster-id exchange (ctor, all
  ranks) for diff-0/diff-1 classification.
- **ctor** stores `lts_cluster_id_`, builds ghost ids; **`AdvanceADERClusterBulk`**
  gains `cluster_id_for_seam` (-1 = np=1 no-op) + calls the seam corrector before
  the mass inverse; **`n_ghost_exchanges_`** counter + `HasSharedFaces()`.
- **Pre-existing bug fixed:** `fault_shared_faces_` wrongly tagged every seam as
  fault-shared when `fault_attr == 0` (`0 == 0` match) — broke fault-free LTS at
  np>1.  Guarded with `if (bc_.fault_attr > 0)` (mirrors fault_interior_faces_).
- **driver:** `ReduceGlobalMeta` MPI_Allreduce hook; gate uses **global** fault
  count (drops the `nprocs==1` guard + its REVIEW-A-1 divergence); per-sync
  matched-collective counter (`P-007`); `num_state==NUM_STATE` assert.
- **test** `tests/unit/test_lts_mpi_seam.cpp` (+ Makefile): Nc=1 np>1 vs np=1,
  partition-invariant reductions (energy, centroid-weighted energy, max elem norm).

## Local gates (all GREEN)
| Gate | Result |
|---|---|
| seam gate np=2 / np=3 / np=4 vs np=1 | 58/58 · 87/87 · 116/116; Δenergy=4.4e-16, Δcentroid=0.0, Δmax=0.0 |
| **gate sensitivity** (seam corrector disabled) | np=2 diverges **65%** energy, 22/58 fail — the gate is real, not false-green |
| per-sync ghost-exchange counter (P-007) | matched all syncs, all ranks |
| field-evolved negative control | passes (Q changed over K syncs) |
| LTS unit suite | layout 102 · scheduler 31 · clustering 4043 · partition 19 · predictor 31 · reorder 14 |
| fault byte-exact (tpv102 ADER) | 4/0 |
| driver build | exit 0, functional |
| pre-existing `parallel_wave` abort | confirmed identical on HEAD (stash test) — not a regression |

## Adversarial review (2 forks) — both PASS-WITH-FIXES, no CRITICAL/hang/divergence
**Reviewer A (seam corrector numerics):** bit-exactness = faithful GTS subset;
no hang (P2P exchange, unconditional per Correct(c), n_shared==0 return safe);
`has_bulk_bg_` removal safe (SharedInteriorFaceFlux_ never reads bulk_bg_); diff-1
fail-loud complete; empty-cluster `I_[c]` sizing safe.
**Reviewer B (driver + test):** gate rank-consistent (global fault count); MPI hook
types correct; counter identity holds; V2 restart enforces np match.

### Fixes applied (all verified)
- **[MODERATE] counter false-abort when a rank has `n_shared==0` at np>1** — the
  seam corrector early-returns without incrementing, so `count(0) != expected`
  would abort a correct run (isolated subdomain).  Fixed: driver expects
  `HasSharedFaces() ? sched : 0`.
- **[LOW] uninitialised `I_[c]` pack** — pack only cluster-c blocks + zero the rest
  (bit-identical; UB-clean).
- **[LOW] `num_state==NUM_STATE`** assert added where meta is built.
- **[LOW] test negative control** — assert the field evolved (>1e-6) so a no-op
  stepper can't pass; **confirmed** the gate FAILS with the seam flux off.
- **[LOW] stale np==1-only comment** (driver) updated; **[LOW] test** ref file now
  stamps + validates `(serial_ne, K)` and the message says "centroid-weighted".

## Deferred to Stage 2 (documented, not silently skipped)
- Cross-cluster (diff-1) rank seams: coarse forecast-`D(k)` ghost exchange +
  coarse-side local-buffer accumulation + extended `n_collectives` + the B.11
  np=2==np=1 multi-cluster gate + a dedicated diff-1 fail-loud abort test.
- Fault half at np>1 (fault stepper keeps its `nprocs==1` guard).
- Multi-cluster (Nc>1) at np>1 (only diff-0 exercised at Nc=1 here).
