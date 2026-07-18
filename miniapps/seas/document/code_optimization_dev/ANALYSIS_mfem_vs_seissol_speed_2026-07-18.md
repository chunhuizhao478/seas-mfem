# ANALYSIS — MFEM vs SeisSol throughput on the v4_0_0 coarse ALT mesh

**Date:** 2026-07-18
**Runs compared:** MFEM `seas_spatial_dyn_driver` jobs 52228747/52228815/52228818
(p1/p2/p3, GTS ADER pure-upwind, 2 nodes × 128 MPI ranks = 256 cores, `cfl_dg_safety=3`)
vs SeisSol job 52042189 (o4=p3, clustered rate-2 LTS, same mesh
`mesh_alt.puml.h5` ≡ `safalt_dl_safgh_2M_walls_mmg.msh`, 1,319,294 tets,
2 nodes × 8 ranks × 16 threads = 256 cores, EndTime=150 s **completed**).
**Method:** multi-agent source+log analysis, adversarially re-derived; every number
below survived independent re-derivation (closure within 8%).

---

## 1. The scoreboard

| Run | dt | steps for 150 s | sim-s/hour | 150 s ETA |
|---|---|---|---|---|
| SeisSol o4 + LTS | 41.2 µs … 42 ms (11 clusters) | — | **31.84** | **4 h 43 m** (done) |
| MFEM p1 GTS (safety 3) | 32.06 µs global | 4.68 M | 0.866 | ~7.2 days |
| MFEM p2 GTS (safety 3) | 19.24 µs global | 7.80 M | ~0.19 | ~32 days |
| MFEM p3 GTS (safety 3) | 13.74 µs global | 10.9 M | ~0.05 | ~4 months |

Observed gap (p1 basis): **36.8×** (31.84 / 0.866).

## 2. The verified factor decomposition (closes to 36.78 = observed)

```
36.8  =  38.73      ×   3      ×   3/7        ×   0.739
         (LTS)        (safety)   (order CFL,     (per-update cost —
                                  FAVORS p1)      MFEM p1 is FASTER)
```

- **38.73× — clustered LTS element-update ratio** (the structural factor).
  GTS at dt_min does 4.80×10¹² element-updates for 150 s; SeisSol's rate-2
  clusters do 1.24×10¹¹. NOTE: SeisSol's log prints "theoretical speedup 111.85"
  — that is an **arithmetic mean** of per-cell dt_cluster/dt_min
  (`ClusterLayout.cpp:74-109`), which overweights large-dt cells; the honest
  workload number is the harmonic form **38.73×**. Fault (DR-face) work benefits
  less: **22.0×** (DR faces concentrate in fine clusters 4–5).
- **3× — `cfl_dg_safety=3`.** MFEM's dt formula at safety=1 is **algebraically
  identical** to SeisSol's, verified in both sources: same length scale
  (6V/A = 2×inradius; `bimaterial_wave_operator.inl:1147-1166` vs SeisSol
  `GlobalTimestep.cpp:27-46`), same per-element cp = √((λ+2µ)/ρ), same cfl=0.5,
  same MPI-min. Numeric identity on this mesh: 41.2181 µs × 7/3 = 96.176 µs =
  MFEM's printed safety=1 dt (5+ digits; same critical element to 5.5×10⁻⁶).
  *Eliminated 2026-07-18: the speed decks now hard-wire `cfl_dg_safety=1.0`
  (commit 3b52cf9).*
- **3/7 — order CFL factor, a 2.33× ADVANTAGE for MFEM p1** vs o4 (1/(2N+1)
  rule). Beware double counting: the raw dt ratio 41.22/32.06 = 9/7 already
  contains BOTH the safety-3 and this factor.
- **0.739 — per-element-update wall cost: MFEM p1 is 1.354× FASTER per update**
  than SeisSol o4 (9.89×10⁶ vs 7.31×10⁶ updates/s on identical 256 cores).
  Independent cross-check: FLOP ratio o4/p1 ≈ 8.37 ÷ kernel-efficiency ratio
  6.18 (SeisSol 3.11 vs MFEM ~0.50 GFLOP/s/core) = 1.354 — exact match.

**Conclusion: essentially ALL of SeisSol's wall-clock advantage is stepping
policy (LTS + the safety factor), not per-update speed.** MFEM p1 updates
elements faster than SeisSol o4 does; it just updates each of them 30–1000×
more often than their own CFL requires.

## 3. Where MFEM's per-step time actually goes (p1, 256 ranks)

- **Compute/overhead-bound, NOT collective-bound** (adversarially confirmed):
  5 global collectives + 11 halo rounds per step ≈ 1–3 ms latency vs the
  observed ~133 ms/step (<3%).
- Face-corrector per-QP geometry/transform overhead ≈ 45–50% of the step
  (on-the-fly `GetFaceElementTransformations`/`CalcShape` per QP per step;
  `--face-cache` exists but was OFF).
- Kernel efficiency ~0.50 GFLOP/s/core vs SeisSol's 3.11 (generated libxsmm) —
  a ~6× efficiency gap that p1's ~8× lower FLOP/element happens to offset.
- Fault-rank load imbalance wait ~10–25% tail (prior 2026-06-24 profile:
  pure imbalance ~8.4%, guard collectives since hoisted).

## 4. SeisSol-side facts worth keeping

- Rupture window (t=0–50 s) ran 1.42× slower than post-rupture (25.7 vs 36.5
  sim-s/h) — rs-fast (FL=103) Newton iteration growth while V is high;
  `computeDynamicRupture` = **57.3% of SeisSol's kernel time** while only 26.5%
  of its FLOPs. DR is expensive for SeisSol too.
- Output overhead negligible (<0.5% blocking); setup 34 s.

## 5. Projections (what each lever buys MFEM, p1 basis)

| Configuration | sim-s/hour | 150 s ETA | vs SeisSol |
|---|---|---|---|
| GTS, safety 3 (as measured) | 0.87 | 7.2 d | 36.8× slower |
| GTS, **safety 1** (now default in decks) | 2.60 | 2.4 d | 12.3× slower |
| **LTS + safety 1** (projected, zero-overhead) | ~100 | ~1.5 h | **~3× FASTER** |

The LTS+safety projection (0.866 × 38.73 × 3 = 100.6 = 31.84 × 7/3 × 1.354) is
an upper bound: it assumes wall ∝ element-updates, zero cluster-management
overhead, and LTS-aware load balance. Even at half efficiency it overtakes
SeisSol o4 — because p1 is cheaper per update and takes 2.33× larger steps.
p3+LTS would land near SeisSol o4 (same FLOP class, minus kernel-efficiency gap).

## 6. Actions

1. **DONE (3b52cf9):** `cfl_dg_safety = 1.0` hard-wired in the three speed decks
   (formula-identical to SeisSol; stability on this mesh to be confirmed by the
   next runs — SeisSol runs the identical dt stably at o4).
2. **Clustered LTS implementation** — the structural fix; plan in
   `document/lts_dev/` (in preparation). Payoff ceiling on this mesh: 38.7×
   fewer bulk updates, 22× fewer fault updates.
3. Secondary (post-LTS) levers, in expected order of value: `--face-cache`
   (attacks the 45–50% face-corrector overhead), kernel efficiency (the ~6×
   GFLOP/s gap), LTS-aware partition weights (fault imbalance).

## 7. Caveats carried forward

- MFEM steps/s from a single early-run window (~7.5 steps/s at step 11.2k);
  rupture-phase cost growth (friction Newton) not yet measured on MFEM.
- p2/p3 rates estimated from shorter windows (4.5k / 1.6k steps).
- Caliper per-region breakdowns require a clean-exit run (short `SAFS_TFINAL`);
  the 15 h walltime-killed jobs write no region report.
- SeisSol's local checkout is newer than the job binary; formulas matched the
  printed numbers exactly, so version drift is immaterial here.
