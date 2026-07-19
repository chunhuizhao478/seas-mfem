# Code Review: LTS Phase 4 Step 0 partial (P4S0-a reorder activation + P4S0-3 partition) — 2026-07-19

## Scope
- Commits: `b137cbc` (driver reorder activation), `e6cc118` (METIS partition companion).
- Files: `drivers/spatial_dyn_driver.cpp`, `dynamic/lts_partition.{hpp,cpp}`, `tests/unit/test_lts_partition.cpp`.
- Plan: `PLAN_clustered_lts_ader_2026-07-18.md` §"Phase 4 Step 0", P-006, P-009, D-1/D-2.
- Method: one adversarial reviewer. **No CRITICAL findings.** Core logic verified sound: per-element `dt_e` binning bit-identical to the operator path for every reachable config; the checkpoint V2 perm round-trip + restart hash correct; `lts="off"` byte-identical.

## Findings + resolutions

### [A-1] MODERATE (np>1 bug) — `lts_stepping` gate omitted `nprocs==1`
`num_fault_total` is a LOCAL count, so at np>1 a fault-free-bulk + lts run had every rank enter the sync loop and abort at `MFEM_VERIFY(nprocs==1)`; a fault mesh with faults on only some ranks diverged (fault-free ranks abort, fault ranks run GTS) → MPI hang. **FIXED:** `lts_stepping = LtsEnabled() && num_fault_total==0 && nprocs==1`; the canary banner also gated on `nprocs==1`. Now all np>1 ranks take the GTS loop uniformly (the "running GTS at np>1" behavior the log promises).

### [A-2] LOW (latent footgun) — missing `T->SetIntPoint(&ip)` before `material.EvalAt`
The operator's `MaterialAtLocal_` (P2-005) sets the int point first; the pre-operator clustering eval did not. No divergence today (current sidecar coefficients read `ip` directly), but a future coefficient reading `T.GetIntPoint()` would get a stale point → different `c_p` → different `dt_e` binning than the operator. **FIXED:** added `T->SetIntPoint(&ip)` mirroring the operator.

### [A-3] LOW (unreachable) — `Mode::GridFunction` would abort deep in `EvalAt`
The spatial driver never builds GridFunction material, but the builder would hard-abort inside `EvalAt` if it did. **FIXED:** added an explicit `MFEM_VERIFY(mode != GridFunction)` with a clear message.

### [B-2] LOW (robustness) — empty cluster → zero-total METIS constraint
`BuildLtsAwarePartition` is a public entry point with no precondition check; a multi-constraint call with an empty cluster id gives METIS a zero-total constraint (version-dependent). Not reachable via `BuildLtsClustering` (contiguous ids). **FIXED:** validate cluster ids in `[0, num_clusters)` and (multi-constraint) every cluster non-empty.

### [B-6] MODERATE (test quality) — fault-lock never exercised; locality checked possibly-vacuously
T4 asserted `VerifyLtsPartitionFaultLocality==0` but never that the partition was non-trivial nor that the lock did anything (on the small grid raw METIS kept the seam intact). **FIXED:** T4 now asserts ≥2 ranks used; new T4b uses a dumbbell (two triangles bridged by the fault edge) so the balanced min-cut MUST sever the fault pair — raw cut > 0, locked cut == 0 (the lock is deterministically exercised).

### [B-7] LOW — `llround` on a value > LLONG_MAX is implementation-defined
`cell_cost = 1e30` could return a negative sentinel, defeating the clamp. **FIXED:** clamp the double to the ceiling before `llround`.

### Not changed (documented LOW)
- **A-5 fragility:** `compute_lts_layout_hash` hardcodes `rate=2` while the stepping block uses `opt.rate` (also 2). Both are 2 today; a theoretical mismatch only if `LtsClusteringOptions::rate`'s default changes. Left as-is (the plan fixes rate-2).
- **`--lts-report`:** still recomputes clustering from the operator-sourced builder (bit-identical for reachable configs; diagnostic-only).

## Verdict
PASS WITH FIXES — all applied. `test_lts_partition` 19/0; driver builds+links; LTS units green (predictor 31, reorder 14, friction 10, nucleation 62); byte-exact-off parity 36/36. np>1 parity + partition balance remain Frontera-only.
