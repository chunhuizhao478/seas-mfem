# Mechanism-level gap analysis: MFEM vs SeisSol — where the 13× actually lives

**Date:** 2026-07-24 · **SeisSol source:** local checkout `~/projects/SeisSol` (6dc8431), the v1.3.1
lineage that ran job 52365078 · **MFEM state:** λ=1 baseline, 3123 s/sim-s (52422891, reproduced
3164 = +1.3 % by 52472765).

**Why this document exists.** A0 hit (1.370×), A1 missed (0.988×). The misses share a cause: levers
were picked from *our* code's structure ("125 rounds looks like a lot") instead of from a verified
account of *which mechanisms produce SeisSol's advantage*. This is that account — every SeisSol-side
claim checked in its source, every MFEM-side cost tied to a measured number.

## The accounting (measured)

| | MFEM (λ=1) | SeisSol (as-run) | gap |
|---|---:|---:|---:|
| total | 3123 s/sim-s | 239 | 13.1× |
| compute | **1871** (59.9 %) | **111** | **16.9×** |
| exposed MPI wait | **1251** (40.1 %) | ≈ 0 (its non-compute 128 is blocking IO at our matched output cadence) | — |

Two independent problems. Fixing either alone floors at 7.8× / 5.2× respectively.

## Mechanism ledger (verified in source)

### M1 — kernel execution: generated GEMM chains vs. generic loops
**SeisSol:** every element operation is a chain of small dense GEMMs over precomputed constant
reference-element matrices, emitted by yateto → libxsmm/PSpaMM generated assembly
(`cmake/FindLibxsmm_executable.cmake`; `src/Kernels/` contains **no runtime quadrature at all**).
Sustained ~6.4 GFLOP/s/core (RESULTS_p5 §2).
**MFEM:** `DenseMatrix::Mult` + hand loops + whole-vector AXPYs; ~0.4 GFLOP/s/core equivalent.
Measured cost: `Vector::Add` 18.9 % + `operator=` 9.3 % **self**-time = 28 % pure vector traffic
(B0 perf); bench variants B/D = **2.3–8.3×** on the predictor structure under Rome contention.
**Target:** predictor 871 s/sim-s + volume share of corrector-rest 710. **This is B1/B2 — the
attack is right; only its share of compute (46.6 %, measured) was ever wrong in the plan.**

### M2 — zero runtime per-QP work: precomputed per-face flux operators
**SeisSol:** the rotated Godunov operators `n·A±·N⁻¹` are built **once at setup** per cell-side
(`src/Initializer/CellLocalMatrices.cpp:320,338` — `localKrnl.AplusT = localIntegration[cell].nApNm1[side]`)
and applied as GEMMs. Zero flux arithmetic per QP per step. Valid because linear elasticity +
per-element-constant material ⇒ the face flux is a constant linear operator.
**MFEM:** per-QP Godunov at runtime. `--face-cache` (measured **1.31×** on GTS) removes only
geometry+basis, not the flux arithmetic. **But we already own the SeisSol mechanism:**
`dynamic/precomputed_face_fluxes.cpp` is wired into the ADER face paths
(`wave_operator.inl:4372, :4849`) — and it has been **OFF in every benchmark leg**, because it is
mutually exclusive with `--face-cache` (`:741`) and the decks all pass `--face-cache`.
**Nobody has ever measured our own precomputed-flux path on this benchmark.**
(Coverage on the LTS corrector's own face sweep is unverified — check before the A/B.)

### M3 — data layout: element-local blocks vs. global strided sweeps
**SeisSol:** per-element dense blocks, cluster-contiguous; operands stay in cache across the CK
recursion. **MFEM:** component-major global vectors; the predictor sweeps the whole state per level
(the 28 % traffic above; bench A degrades 33→79 µs under contention while element-local B stays
flat ~13). Same fix as M1 (B2 tiling) — they are one program, not two.

### M4 — parallel layout: 16×15 hybrid + async ghost clusters vs. 256 pure-MPI + blocking exchange
**SeisSol:** 16 MPI ranks × 15 OMP threads (2 nodes), and communication runs as **dedicated ghost
time clusters** with `MPI_Test`-based async progression (`src/Solver/TimeStepping/
DirectGhostTimeCluster.cpp`, `GhostTimeClusterWithCopy.cpp`, `AbstractGhostTimeCluster.cpp`) —
comm is a schedulable task overlapped with interior compute, and intra-rank imbalance is absorbed
by threads.
**MFEM:** 256 single-thread ranks; every tick's exchange is a blocking `ExchangeFaceNbrData`.
Measured: wait = 40 % of wall; **A1 proved the wait is not in the message count** (halving rounds
at fixed sync count moved nothing — per-round wait exactly doubled); skew Max/Avg = 2.18.
**The untested variable is rank count.** 16× fewer ranks ⇒ 16× fewer seam partners and far less
arrival spread to expose. We have *never timed* MFEM at any count but 256 on this problem —
the one layout leg that ran (B0 Leg 3, 128 ranks/node vs 64) was never timed (harness gap).

## Why our misses happened, in one sentence each
- **A1:** optimized message count inside a mechanism (M4) whose cost driver is arrival skew.
- **Kernel plan v1:** optimized the face stage from a laptop profile that inverted on the real machine.
- Both were *within-our-structure* picks. The ledger above is *cross-implementation*, and each row
  ends in a measured number or a named experiment.

## Directions (each states what it decides and what it costs)

**D1 — the deciding experiment (no code, one job, 3 short legs + counters): rank-count sweep.**
Run the λ=1 deck at **256, 128, 64 ranks** on the same 2 nodes, with per-rank tick-arrival trace and
`perf stat` FLOP counters on one leg.
Decides: (a) does wait collapse with rank count (M4 / the SeisSol-16-rank hypothesis) — if yes, the
production fix may be as cheap as *running wider nodes with fewer ranks* or as structural as OpenMP;
(b) systematic-vs-jitter skew (sizes A6: bounds today 1.07–1.67×);
(c) memory-bound vs FLOP-bound (the old R-005 — decides whether B2's traffic mechanism can work);
(d) pays the never-timed Leg-3 debt. **This is the next Expanse job.**

**D2 — measure our own SeisSol mechanism (no code, one A/B leg): `--precomputed-face-fluxes` vs
`--face-cache`** on the λ=1 deck (GTS and LTS legs). Decides M2's remaining headroom with a
mechanism we already implemented. Prep: verify the LTS corrector honors the flag; verify memory
fits at np=256.

**D3 — B1 now (code, small, bench-backed): fuse the predictor's accumulation sweeps** — unchanged
from the plan, still the best-supported code change (28 % measured traffic target, ≥1.2× stage gate,
≤1e-12). B2 (tiling / per-element dgemm = the M1+M3 program) proceeds only after D1's counter leg
confirms memory-bound.

**D4 — deferred until D1 reports:** A6 rebalancing (only if skew is systematic), OpenMP hybrid
(only if wait collapses with rank count), seam-compute kernel (12 % of step, unclaimed).

**Sequencing:** D1 and D2 are one submission each and can share a job; D3 starts locally now.
Nothing else is authorized by current evidence.

## Honest ceiling, restated
If D3 lands at its bench central (~4× on predictor+volume) and D2 recovers ~1.3× on faces, compute
→ ~600–800 s/sim-s. With wait untouched that is ~1850–2050 total ≈ **7.7–8.6×** — consistent with
the corrected appendix. **Breaking under ~7× requires the M4 side**, which is exactly what D1 sizes.
Parity remains off the table at matched order; the far-field order drop stays a separate proposal.
