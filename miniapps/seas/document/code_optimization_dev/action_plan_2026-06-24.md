# Action Plan — SEAS ADER spatial-dyn optimization (2026-06-24)

## 1. Goal and constraints

**Goal.** Reduce wall-clock time of the SEAS ADER spatial-dyn driver from the measured **830.65 s** baseline (run 51443560, ALT rate-state upwind-ADER, SAFS triq mesh on Expanse, 200 ranks, 1445 macro-steps) without changing the numerical result beyond documented tolerance.

**Hard constraints (binding on every item below).**
- **Numerics-preserving only.** Every adopted change is either byte-exact or within a *documented, tested* tolerance (Lever-1/3 cache: `<=1e-12` relative; Lever-2 / CK-fusion: bit-exact `0.0`). The TPV*/BP5 byte-exact regression contract must not change. No change may flip a *default* that those drivers depend on.
- **Editable surface:** everything under `dynamic/` (`wave_operator.inl`, `wave_operator.hpp`, `godunov_flux*`, `fault_face_flux*`, `elem_derivative_cache.hpp`, `friction_substep_iterator*`, `bimaterial_wave_operator.hpp`) and `drivers/spatial_dyn_driver.cpp`, plus the SAFS sbatch under `jobs/`.
- **NO-TOUCH (shared with BP5/TPV byte-exact regression):** `bp5/`, `bp1/`, `bp2/`, `domain/`, `fault/`, `solver/`, `friction/dieterich_ruina.hpp`. Reading `fault/fault_basis.hpp` is allowed; editing it is not.
- **Validation reality:** the end-to-end SAFS fingerprint is only reproducible on Frontera/Expanse. Local checks are construction-order and small-fixture unit/parallel tests only. Any non-bit-exact change requires a Frontera/Expanse fingerprint re-validation before it ships to production.

**Crux of the baseline.** The two opt-in levers that eliminate the redundant ADER work were present in the binary but **not passed by the production sbatch**: the run executed `DerivMode::OnTheFly` (no `--deriv-cache`) and the *separate* two-recursion CK path (no `--shared-ck-recursion`). So `ApplySpatialDerivative` ran **twice per step** (222.5 s + 221.4 s = 53.4% of runtime) over the full per-QP quadrature machinery on an affine P1 mesh. The largest wins are therefore *activation of already-landed, already-tested code*, not new kernels.

---

## 2. Ranked optimization backlog

Ranked by reviewed impact / effort. "Expected speedup" is the **reviewed/realistic** figure, not the proposer's optimistic one. "No-touch?" = yes means the change requires editing a no-touch file (none of the adopted items do; all adopted items are flagged No).

| ID | Title | Target region | Reviewed speedup | Effort | Risk | Recommendation | No-touch? |
|----|-------|---------------|------------------|--------|------|----------------|-----------|
| **OPT-CK-FUSE** | Activate `--shared-ck-recursion`: one CK recursion feeds substep states + integral | `AdvanceADER` 2nd `ApplySpatialDerivative` subtree (26.6%) | **~210–225 s (~25–27%)** wall; bit-exact | low (flag) | low | **adopt** | No |
| **OPT-VOL-ENABLE-LEVER3 + Lever 1** (`--deriv-cache`) | Activate cached `D_d^e` / `S_d^e` matvec; kills per-QP Jacobian re-inversion | `ApplySpatialDerivative` (53.4%) + `ComputeVolumeRHS` (10.1%) | combined w/ CK-fuse → **~1.7–2.0× overall (830→~400–480 s)**; tol `<=1e-12` | low (flag) | medium | **adopt-with-care** (needs bimaterial parity test + fingerprint check) | No |
| **OPT-FACEFLUX-ALLREDUCE-HOIST** (Lever 4) | Hoist the two per-step 1-int `has_bulk_bg_` guard Allreduces to one-time setup consensus | face-flux + shared-face guard barriers (~62.6 s wait) | **~10–30 s (1–4%)** wall; byte-exact | low | low | **adopt-with-care** (option (a) only) | No |
| **OPT-CK-FUSE-DEFAULT** | Make fused CK the driver default (keep `--no-shared-ck-recursion` opt-out) | same as OPT-CK-FUSE | 0 additional; ergonomics/robustness | low | low | **adopt-with-care** (keep opt-out; don't bundle `--deriv-cache` default) | No |
| **OPT-FUSE-3DIR** | Fuse the 3 directional derivative applies + star-matrix apply into one element pass | one CK recursion's spatial-operator cost | **~1.1–1.3× on that recursion**; bit-exact gate exists | medium | medium | **adopt-with-care** (re-baseline after levers; microbench gate `>1.1×`) | No |
| **OPT-LB-WEIGHTED-PARTITION** | Fault-weighted METIS partition to balance friction cost | friction/fault imbalance | **~5–12 s (~0.6–1.5%)**; NOT byte-exact, fingerprint risk | high | medium | **defer** (Phase 3; needs Frontera validation; bulk de-balance risk) | No |

**Defer / reject (with one-line reasons):**

| ID | Verdict | Reason |
|----|---------|--------|
| OPT-FACEFLUX-SHAPE-CACHE | defer | Sound + byte-exact, but payoff < few s at P1; sequence only after ROT-CACHE proves a residual (it does not — see below). |
| OPT-SHARED-HOIST-BULKBG-GUARD | defer | Subsumed by OPT-FACEFLUX-ALLREDUCE-HOIST done as one correctly-scoped pass over all four guard reductions; standalone double-counts. |
| OPT-DERIV-GEMM-BATCH | reject | 4×4×9 is below dgemm break-even; memory-bound; gather/scatter eats the win (~0.8–1.2×). |
| OPT-DERIV-GROUP-GEMM | reject | Flop-neutral on a bandwidth-bound kernel; needs SoA repack; ~1.0–1.1×. |
| OPT-ADER-ORDER | reject | order=2 is the hard floor (asserts {2,3,4}); dropping it changes numerics — not allowed. |
| OPT-MIXED-PRECISION | reject | fp32 injects ~1e-7 error vs `<=1e-12` gate; would change SAFS fingerprint; self-rejecting. |
| OPT-SOA-LAYOUT | reject | Locality premise false (cached `Qc` block is L1-resident, `D` dominates); ~1.0×; 60+ rewrite sites; perturbs byNODES halo. |
| OPT-VOL-AFFINE-HOIST | reject | Dominated by `--deriv-cache`; silent-wrong-answer hazard if `mesh.order>=2` (curved). |
| OPT-VOL-FUSE-CK | reject | No separate "I assembly" element loop to fuse into; self-contradictory (can't both skip the I-read and stay byte-exact). |
| OPT-VOL-GEMM-BATCH | reject | Flop count identical, `S` already L1-resident; net wash on local netlib BLAS; wrong file list. |
| OPT-FACEFLUX-ROT-CACHE | reject | Production CVM run uses `BimaterialWaveOperator` with per-face precomputed matrices — the per-QP `BuildRotation` it targets is not in the hot path (~0 s). |
| OPT-FACEFLUX-FAULT-ROT-REUSE | reject | <0.5% gain off critical path + a real MPI-determinism bug (keys cache on `elem1_on_plus` for both branches; shared branch needs `qpd.sign_flipped`). |
| OPT-SHARED-HOIST-RECONCILE-GUARD | reject | Barrier relocates to the mandatory Allgather/Allgatherv ~60 lines later; prior plan lists this region "not worth optimizing." |
| OPT-SHARED-COMBINE-GUARDS | reject | `any_shared` MAX is a deadlock gate inside the reconcile, not a foldable sibling; subsumed by the hoist. |
| OPT-SHARED-BATCH-HALO | reject | Target is 0.31 s; reintroduces the R-004 byNODES unpack that already caused a SAFS normal-velocity runaway. |
| OPT-SHARED-OVERLAP-HALO | reject | Hideable payload ~1.15 s; targets imbalance (not latency) which overlap cannot hide. |
| OPT-SHARED-ALLGATHERV-TO-NEIGHBOR | reject (revisit >1000 ranks) | Reconcile is latency-bound and sub-MB at 200 ranks; ~0 s; abort risk if a fault peer is missed. |
| OPT-LB-FAULTLOC-PARTITION | reject | Halo wait is over the *whole* shared-face graph, not fault seams; ~0 s and shifts the asymmetric-mesh sigma_n fingerprint. |
| OPT-LB-HOIST-GUARD-ALLREDUCE | reject as speedup (= Lever 4 cleanup) | A hard per-step `V_max` Allreduce remains at `spatial_dyn_driver.cpp:3359`; removing the guards relocates the wait. Land only as the byte-exact cleanup in Phase 0, not sold as 4–6%. |
| OPT-LB-OVERLAP-HALO | reject | Premise false (CK work is one fused pass *before* the exchange loop); no driver-level Isend/Irecv to split; targets imbalance. |
| OPT-AFFINE-PRECOMP-NOTE | adopt (as guidance) | Not an optimization — a "do not hand-optimize OnTheFly" note; correct, annotate the byte-exact-contract reason. |

---

## 3. Phased roadmap

> **Sequencing principle:** activate the landed levers first (Phase 0), measure, *then* re-baseline every new-code estimate against the post-activation caliper before spending effort. Several proposer estimates were computed off the *pre-activation* profile and double-count work the levers already remove.

### Phase 0 — Activation + guard hoist (quick, low-risk)

This phase is dominated by the two highest-ROI items and a byte-exact micro-cleanup. Expected: **~27% from CK-fusion + ~1–4% from the guard hoist before any kernel rewrite.**

---

#### 0A. OPT-CK-FUSE — activate `--shared-ck-recursion` (bit-exact)

- **Edit site:** `jobs/safs/safs_expanse/spatial_dyn_ratestate_v3alt_vw005_cvm_upwind_ader_caliper_triq_expanse_2N100r_compute_4hr_safs.sbatch`, the `srun` line (~L253–267).
- **Exact change:** append `--shared-ck-recursion` to the driver argument list. No source edit. This flips `use_shared_ck=true` at `drivers/spatial_dyn_driver.cpp:1316`, so `AdvanceADERWithSubStep_Spatial` calls `ComputeADERSubStepStatesAndIntegral` (L461) and `AdvanceADER(..., &shared_I)` (L513), running the CK recursion (and `ApplySpatialDerivative`) **once** instead of twice.
- **Why it is the highest-value action:** deletes the entire `AdvanceADER`-internal `ApplySpatialDerivative` subtree (221.4 s avg) — pure element-local compute on the critical path, not a barrier wait, so the saving does not merely relocate.
- **Correctness / verification:**
  - Guarded by `tests/unit/test_wave_operator_spatial_derivative.cpp` `TestSharedCKParity` (L405–475): asserts `max|merged − separate| <= 0.0` (exactly bit-exact) for both `Q_per_node` and `Q_new`, across orders {2,3,4} and **both** DerivModes. Run it (`ctest -R SpatialDerivative` or the unit binary) and confirm it passes.
  - `--shared-ck-recursion` alone is bit-exact at fixed DerivMode; do **not** bundle `--deriv-cache` here if the goal is to keep this change byte-identical to the OnTheFly baseline. (0B handles the cache separately.)
  - No new test required — the bit-exact gate already exists.
- **Measurement plan:** re-run the caliper sbatch with the flag. The `AdvanceADER → ApplySpatialDerivative` region (221.4 s) should drop to ~0; `AdvanceADERWithSubStep` should fall by ~210–225 s; `ComputeADERSubStepStates → ApplySpatialDerivative` (the friction-substep recursion, ~222 s) remains. Total wall ≈ 605–620 s.
- **Rollback:** remove the flag from the sbatch (one line).

---

#### 0B. OPT-VOL-ENABLE-LEVER3 + Lever 1 — activate `--deriv-cache` (tol `<=1e-12`)

- **Edit site:** same sbatch `srun` line. Append `--deriv-cache`.
- **Exact change:** flips `deriv_mode_` to `DerivMode::Cached` (`drivers/spatial_dyn_driver.cpp:1299–1308` → `SetDerivMode(Cached)`), which builds `elem_deriv_op_` (Lever 1) and `elem_volume_op_` (Lever 3) at `wave_operator.inl:1176–1178`. `ApplySpatialDerivative` then takes the cached matvec branch (`wave_operator.inl:1029`) and `ComputeVolumeRHS` takes the cached branch (`wave_operator.inl:869`), eliminating per-QP `CalcShape`/`CalcPhysDShape`/Jacobian-inversion and the per-element heap allocs.
- **Correctness / verification — this is the adopt-with-care gate:**
  - Existing guard `tests/parallel/test_wave_operator_cached_parallel.cpp` validates Cached==OnTheFly to `<=1.4e-15` on np2/np4 — **but only for the scalar `WaveOperator`**. The production CVM run uses `interior_flux="matrix"` → `BimaterialWaveOperator` with per-element `A_d^e` from `owned_flux_pool_`, which this test does **not** cover.
  - **New test required (editable, under `dynamic/` tests):** add a bimaterial cached-vs-OnTheFly equivalence check on a small two-material np2/np4 fixture, asserting `<=1e-12` relative on `ComputeVolumeRHS` and `ApplySpatialDerivative` outputs. Land this before trusting the flag.
  - **Frontera/Expanse fingerprint check required:** Lever-1/3 are NOT bit-exact (R-002 round-off re-association); confirm the SAFS output fingerprint is within SAFS tolerance against the reference before this becomes the production config.
  - Driver default stays `OnTheFly` — TPV*/BP5 byte-exact regression is untouched (different drivers, and the flag is opt-in).
- **Measurement plan:** with both flags, `ApplySpatialDerivative` (surviving substep recursion) should drop 4–5× (per the bench_ader_hotpath measurement of the cache at P1, where the per-QP Jacobian-inverse fixed cost is the dominant thing removed); `ComputeVolumeRHS` should drop ~1.5–2.5× (order=1-pessimistic — affine tets make OnTheFly geometry cheapest at P1, so do not expect the 3–5× microbench figure). Combined Phase-0 end-state ≈ **400–480 s wall**.
- **Rollback:** remove `--deriv-cache` (one line). Independent of 0A.

---

#### 0C. OPT-FACEFLUX-ALLREDUCE-HOIST (Lever 4) — one-time bulk-bg consensus (byte-exact)

- **Edit sites:** `dynamic/wave_operator.inl:3854` (face-flux guard) and `:4920` (shared-face guard); `dynamic/wave_operator.hpp` (new member). Also fold in the equivalent guards flagged by OPT-LB-HOIST-GUARD-ALLREDUCE — same two sites — so the four redundant per-step 1-int reductions are removed in one pass.
- **Exact change — use option (a), not (b):** add the consensus `MPI_Allreduce(MPI_MIN, has_bulk_bg_)` **once inside `SetAbsorbingBackground`'s ParMesh branch** (mirroring the existing `SetMixedFluxMode` one-time consensus at `wave_operator.inl:1639–1664`), caching the verified result in a new `bool bulk_bg_consensus_ok_` member. Replace the per-step collective at L3854/L4920 with a rank-local `MFEM_VERIFY(has_bulk_bg_)` tripwire only.
  - **Do not** use option (b) (a separate `VerifyBulkBgConsensus()` driver call): it silently drops the fail-loud guard for the ~120 unit/parallel tests and the tpv102/104/205 drivers that call `SetAbsorbingBackground` then `AdvanceADER` without any new verifier call. Placing the consensus inside `SetAbsorbingBackground` keeps every caller covered.
- **Correctness / verification:**
  - `has_bulk_bg_` has exactly one writer (`SetAbsorbingBackground`, `wave_operator.hpp:393`), called once at `drivers/spatial_dyn_driver.cpp:1999` outside the time loop, never mutated in-loop → reducing once is byte-identical. The reduced int feeds only the `MFEM_VERIFY`, never flux math → byte-exact for ADER/CK numerics.
  - Guarded by the existing parallel face-flux tests in `tests/parallel/` and the np2/np4 ParMesh suite; add a debug-build assert that `has_bulk_bg_` is unchanged across the loop. The R-1505 fail-loud deadlock contract is preserved (the same MIN now runs once, earlier, on all ranks).
- **Measurement plan:** the two per-step `MPI_Allreduce` lines (42.4 s + 20.2 s) should disappear from caliper. Expect realistic recovery of only **~10–30 s** — the waits are load-imbalance that partly resurfaces at the next collective (`ExchangeFaceNbrData`, the per-step `V_max` Allreduce at `spatial_dyn_driver.cpp:3359`). Report the *delta*, do not claim 62 s.
- **Rollback:** revert `wave_operator.inl/.hpp` (small, self-contained diff).

---

#### 0D. OPT-CK-FUSE-DEFAULT — make fused CK the driver default

- **Edit site:** `drivers/spatial_dyn_driver.cpp:1316` (invert the `HasFlag` default) plus the dispatch branch L459–468 / L513–514. Keep `--no-shared-ck-recursion` as an opt-out.
- **Exact change:** default `use_shared_ck=true`; retain the two-recursion path behind the opt-out flag for A/B parity diffing. **Do not** delete the branch (it is the reference path the parity test and on-the-fly debugging rely on). Make the existing rank-0 path-log line (L1317–1321) unconditional so job stdout records which path ran. **Do not** also flip the `--deriv-cache` default — keep the OnTheFly baseline intact so the bit-exact-at-fixed-DerivMode guarantee holds.
- **Correctness / verification:** numerically identical to 0A; same `TestSharedCKParity` coverage in both DerivModes. TPV*/BP5 binaries never call `AdvanceADERWithSubStep_Spatial`, so the regression contract is untouched.
- **Measurement plan:** confirm a fresh spatial run *without* `--shared-ck-recursion` now logs the fused path and matches 0A timing.
- **Rollback:** revert the one-line default and the branch.

---

### Phase 1 — Kernel-level fusion of `ApplySpatialDerivative` + star-matrix apply

#### 1A. OPT-FUSE-3DIR — fuse the 3 directional applies into one element pass

- **Edit sites:** `dynamic/wave_operator.inl` — the three call sites that loop `ApplySpatialDerivative` + `ApplyElementJacobian_` over `d=0,1,2`: `ComputeADERTimeIntegrated` (L1341–1348), `ComputeADERSubStepStates` (L1463–1470), `ComputeADERSubStepStatesAndIntegral` (L1580–1584); plus a new fused kernel in the Cached branch and its declaration in `dynamic/wave_operator.hpp`. Fuse **only the Cached path** (require `DerivMode::Cached`); leave the OnTheFly quadrature path alone.
- **Exact change:** replace the per-direction full-N-vector intermediate (`ck_substep_dQ_dxd_buf_`) with a single per-element loop that computes all three `D_d^e Q` blocks and immediately contracts each with its 9×9 star matrix `A_d` into `D_next`, never materializing the per-direction N-vector. Preserve the `d`-order (0,1,2), the inner `(i,c)` accumulation order, and the `ApplyJacobianPerDOF` `a==0.0` structural-zero skip.
- **Correctness / verification:**
  - Element-local, no cross-element coupling → preserves per-`(i,c)` accumulation order, so it can be **bit-exact** relative to the current cached two-pass path. Must pass `tests/unit/test_wave_operator_spatial_derivative.cpp` `TestSharedCKParity` (`dmax <= 0.0`, orders 2/3/4, both DerivModes) against the cached baseline.
  - No MPI/fault hazard: the downstream collective (`ExchangeFaceNbrData` in `EvaluateBulkAtFaultQPsCanonical`) is untouched.
- **Measurement plan / gate before committing effort:** **re-baseline against the Phase-0 caliper** (the headline 1.5–2× was computed off the pre-activation profile that runs the recursion twice). Extend `tests/unit/bench_ader_hotpath.cpp` to measure fused vs two-pass cached recursion; **adopt only if the microbench shows `>1.1×`**, otherwise defer. Realistic ~1.1–1.3× on one recursion = low single-digit % of total.
- **Rollback:** the fused kernel is additive behind the Cached branch; revert to the two-pass loop.

> **Explicitly out of Phase 1:** GEMM/batching of the matvec (OPT-DERIV-GEMM-BATCH, -GROUP-GEMM, OPT-VOL-GEMM-BATCH). At P1 ndof=4 (and even P2 10×10×9) the kernel is below the dgemm break-even and is memory-bound; the verdicts reject all three. Revisit *only* if a post-Phase-0 roofline on Expanse shows the cached matvec is compute-bound (it will not be at P1).

### Phase 2 — Face-flux compute + communication

The two compute-side face-flux items (ROT-CACHE, SHAPE-CACHE) were **rejected** for the production CVM path: it runs `BimaterialWaveOperator`, whose non-fault `InteriorFaceFlux_` already applies per-face precomputed matrices, so the per-QP `BuildRotation`/`DenseMatrix(9,9)` rebuilds those items target are not in the hot path. The halo-batch / overlap items target ~0.3–1.2 s and reintroduce the R-004 byNODES-unpack hazard. **Net: Phase 2 carries no adopted item.** If post-Phase-0 caliper still shows the 145 s face-flux compute as a top region, profile the **fault branch** rotation rebuilds and the per-QP state-assembly/matvec/scatter (not the non-fault path) before proposing anything new. Treat Phase 2 as a *measurement gate*, not a committed change.

### Phase 3 — Load-balanced partitioning

#### 3A. OPT-LB-WEIGHTED-PARTITION (defer; high effort, fingerprint risk)

- **Edit sites:** new `dynamic/fault_weighted_partition.hpp` (mirrors `fault_locality_partition.hpp`) + `drivers/spatial_dyn_driver.cpp` partition selection, behind a new `--partition-fault-weighted` flag.
- **Exact change:** build a weighted dual graph and call `METIS_PartGraphKway` directly (MFEM `Mesh::GeneratePartitioning` hardcodes `wgtflag=0`, so the weighted path is ~300 LOC of dual-graph CSR + METIS dispatch, **not** the trivial union-find of the locality helper). Weight each fault-adjacent element `base_bulk_cost + fault_qp_cost`, where `fault_qp_cost` is **calibrated from the measured friction:bulk time ratio** (not guessed), and combined with the element-count weight (not replacing it).
- **Correctness:** NOT byte-exact — changes MPI reduction order and fault-seam placement → may shift the SAFS upwind sigma_n fingerprint on the asymmetric triq mesh. Requires the full Frontera/Expanse fingerprint re-validation before adoption. `VerifyFaultLocality`-style fail-loud setup assert.
- **Measurement plan:** target the friction-overhang imbalance (`friction_substep` max−avg ≈ 10.5 s; the guard-Allreduce/halo waits are the *same* stragglers, already removed/relocated by Phase 0). Realistic recovery **~5–12 s (~0.6–1.5%)**.
- **Why deferred:** high effort + a Frontera-validation cycle for a single-digit-second gain, with a real risk of **de-balancing the already-3.5%-balanced 53% bulk `ApplySpatialDerivative`** if the fault weight is miscalibrated. Do this last, only if Phase-0 caliper confirms residual friction imbalance dominates after the levers land.
- **Rollback:** default flag off; the unweighted METIS path is unchanged.

---

## 4. Cumulative expected speedup model

Stacking the **realistic** per-phase gains against the 830.65 s baseline (wall-clock is governed by the **max** rank; the levers reduce per-rank compute, the hoists touch barrier waits):

| Stage | Lever | Realistic effect | Running wall estimate |
|-------|-------|------------------|-----------------------|
| Baseline | — | — | **830 s** |
| Phase 0A | CK-fusion (bit-exact) | −210 to −225 s (remove 2nd `ApplySpatialDerivative` subtree) | ~605–620 s |
| Phase 0B | `--deriv-cache` (Lever 1+3) | surviving `ApplySpatialDerivative` 4–5×; `ComputeVolumeRHS` ~1.5–2.5× | ~400–480 s |
| Phase 0C | Guard-Allreduce hoist | −10 to −30 s (imbalance partly resurfaces downstream) | ~380–460 s |
| Phase 1A | Fuse-3-dir (if `>1.1×` gate passes) | low single-digit % on one recursion | ~370–450 s |
| Phase 3A | Weighted partition (if validated) | −5 to −12 s (friction overhang only) | ~365–440 s |

**End-state estimate: ~370–460 s wall**, i.e. **~1.8–2.2× over baseline**, driven almost entirely by Phase 0 activation. This is below the proposer's optimistic ~350 s.

**Imbalance is the floor — be explicit.** Once the redundant flop work is removed, the residual is increasingly MPI-imbalance-bound. The per-step `V_max` Allreduce (`spatial_dyn_driver.cpp:3359`), the `ExchangeFaceNbrData` halo `MPI_Waitall` (14.8 s avg / 42.1 s max), and the macro-step boundary barrier remain hard per-step rendezvous. Hoisting the 1-int guards and even a weighted partition cannot drive wall-clock below the **slowest-rank** compute + irreducible collective latency. The collective-latency floor from the min-rank values (0.42 + 0.16 + 7.09 ≈ 7.7 s of the imbalance trio) is irreducible; the recoverable imbalance is ~60–70 s and only a *fraction* of that is reclaimable without a validated rebalance. Do not model gains below ~360 s without a measured weighted-partition result.

---

## 5. Risks and open questions

**Flagged by verdicts as risky / unsound — do not pursue:**
- **OPT-MIXED-PRECISION (unsound):** fp32 injects ~1e-7 error vs the `<=1e-12` gate; `real_t` is a global MFEM compile-time typedef, so selective fp32 is not a 3-file change. Would change the SAFS fingerprint. Rejected.
- **OPT-VOL-AFFINE-HOIST (silent-wrong-answer):** correct only for affine tets; `drivers/spatial_dyn_driver.cpp:1035` calls `SetCurvature(cfg.mesh.order)` and `mesh.order` is config-reachable `>=2`, where the "invert J once per element" fast path is **wrong** with no guard. Rejected.
- **OPT-FACEFLUX-FAULT-ROT-REUSE (MPI-determinism bug):** keys the cache on `elem1_on_plus` (rank-local) for both branches, but the shared branch must use `qpd.sign_flipped` (rank-independent) — would reintroduce the R-1601 np>1 hang. Rejected.
- **OPT-SHARED-BATCH-HALO / OPT-LB-OVERLAP-HALO (R-004 hazard):** both route through the byNODES vdim=NUM_STATE ghost unpack that already caused a SAFS normal-velocity runaway; must go through `GetFaceNbrElementVDofs`, not the slab formula. Rejected for sub-1% targets.
- **OPT-LB-FAULTLOC / OPT-LB-WEIGHTED (fingerprint shift):** both change fault-seam placement and reduction order on the **asymmetric triq mesh**, which project memory flags as the load-bearing case for the sigma_n/dip seed. NOT byte-exact for output; mandatory Frontera/Expanse fingerprint re-validation.

**Optimizations that secretly need no-touch edits:** none of the adopted/deferred items require editing `bp5/domain/fault/solver/friction`. OPT-FACEFLUX-FAULT-ROT-REUSE *reads* `fault/fault_basis.hpp` (allowed) but would cache in `wave_operator.inl`-side members (rejected anyway). Confirm during implementation that no fused-kernel change reaches into the no-touch friction iterator beyond its existing public interface.

**Open questions to resolve during implementation:**
1. **Bimaterial cache validation gap (blocks 0B):** the only Lever-3 equivalence test covers the *scalar* operator. The production run is `BimaterialWaveOperator`. The new bimaterial parity test (Section 0B) is a prerequisite, not optional.
2. **Caliper call-count discrepancy:** the profile attributes 2890 `MPI_Allreduce` calls (2/step) to `ComputeADERSharedFaceFluxRHS`, but `AdvanceADER` has a single call site (`wave_operator.inl:5741`). The second count is the nested `any_shared` reduce in `ExchangeAndPairSharedFaultQPs` (L6122). Reconcile before trusting any per-call savings magnitude for the shared-face guard.
3. **Post-Phase-0 roofline:** before any Phase-1 GEMM/layout work, profile the cached matvec on Expanse with hardware counters to confirm bound (the verdicts predict memory-bound at P1; if confirmed, all GEMM/SoA items stay rejected).
4. **Weighted-partition calibration:** the fault weight must come from the *measured* friction:bulk ratio in the Phase-0 caliper, not a guess, or Phase 3 de-balances the bulk kernel.

---

## 6. Acceptance criteria

The effort is complete and accepted when **all** of the following hold:

1. **Correctness preserved.**
   - `tests/unit/test_wave_operator_spatial_derivative.cpp::TestSharedCKParity` passes `dmax <= 0.0` (orders 2/3/4, both DerivModes) after CK-fusion and after any Fuse-3-dir change.
   - The new bimaterial cached-vs-OnTheFly equivalence test passes `<=1e-12` relative on np2/np4 before `--deriv-cache` ships.
   - The existing np2/np4 parallel suite (`tests/parallel/test_wave_operator_cached_parallel.cpp` and the face-flux parallel tests) passes unchanged.
   - The SAFS output fingerprint on Frontera/Expanse is within documented SAFS tolerance of the reference for the final production config.
   - The TPV*/BP5 byte-exact regression suite is **bit-identical** (no default flipped that those drivers depend on; `--deriv-cache` and `--shared-ck-recursion` defaults for those binaries unchanged).

2. **No-touch respected.** Zero diffs under `bp5/`, `bp1/`, `bp2/`, `domain/`, `fault/`, `solver/`, `friction/dieterich_ruina.hpp`. All edits confined to `dynamic/` + `drivers/spatial_dyn_driver.cpp` + the SAFS sbatch.

3. **Performance demonstrated by caliper, not asserted.** A re-profiled caliper run shows:
   - `AdvanceADER → ApplySpatialDerivative` region ≈ 0 (CK-fusion).
   - Surviving `ApplySpatialDerivative` and `ComputeVolumeRHS` regions reduced consistent with the cache (4–5× and ~1.5–2.5× respectively).
   - The two per-step guard `MPI_Allreduce` lines removed.
   - **Total wall-clock ≤ ~480 s** (≥1.7× over the 830 s baseline) for Phase 0 alone; any Phase-1/3 item retained only if its own caliper delta exceeds its microbench gate.

4. **Reproducibility.** The production sbatch records (via the unconditional rank-0 log line) which CK and DerivMode paths ran, so future runs are self-documenting.

5. **Rollback available.** Every adopted change is revertible independently (flag removal for 0A/0B; small self-contained diffs for 0C/0D/1A; default-off flag for 3A).
