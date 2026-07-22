# Evidence record: kernel-plan round-2 (microbenchmark + roofline + comm pre-scope)

**Date:** 2026-07-22 · Produced by the evidence workflow (4 agents); feeds PLAN rev 3.
**Bench source:** `tests/bench/bench_ader_kernel_variants.cpp` (committed; Makefile untouched;
compile line in the file header comment). Run logs under the session tmp dir.

## 1. Microbenchmark (Apple M3 Max, 1 core, NE=20000 synthetic p3/O4 predictor)

| variant | µs/elem | GFLOP/s/core | vs A | verdict |
|---|---:|---:|---:|---|
| A — current-structure replica (heap DenseMatrix ops, component-major Q, whole-vector passes, transposed inner loop) | **5.21** | 13.7 | 1.00× | already 48 % of core peak — clang vectorizes the "transposed" loop fine |
| B — tiled-fused (DenseTensor packs, E=64 tiles, 3 levels fused per element, panel gather/scatter) | 4.19 | 17.0 | **1.25×** | real but small; checksum bit-identical to A |
| C — **BatchedLinAlg NATIVE** | 11.2 | 6.4 | **0.46×** | **2.2× SLOWER — eliminated as the mechanism** (serial generic kernels; the header itself warns it) |
| D — per-element LAPACK dgemm | 4.06 | 17.6 | 1.30× | best; tied with B |

**Findings that redirect the plan:**
1. **The kernel's arithmetic/data-structure layout is NOT the production bottleneck.** The
   current-structure replica runs the full O4 predictor at 5.2 µs/elem on a well-fed core —
   ~47× cheaper than the 245.9 µs Expanse per-update (which also includes the 59 % face stage
   and 4.2× hardware/contention). The recoverable *arithmetic* headroom from restructuring is
   ~25–30 %, not 3×.
2. **`BatchedLinAlg` NATIVE is affirmatively the wrong mechanism on CPU** — its CPU path is a
   serial `forall` over generic scalar kernels. If batching is done, the mechanism is
   hand-tiled fused loops (B) or small-dgemm (D).
3. **The surviving argument is memory traffic, not FLOPs:** variant A streams ~8× the DRAM
   bytes of B (~576 MB vs ~72 MB of operator traffic per rep + whole-vector state passes).
   M3's bandwidth masks this; a 128-rank Rome node at ~2 GB/s/core would not. The fusion/
   repack case must be re-made as a *bandwidth* case and proven by a Rome-side spike.
4. Caveats: single-core Apple Si; ICX/AVX2 Rome may vectorize A's strided pattern worse
   (could widen B's edge there); the replica has no virtual calls/MFEM_VERIFY/face-stage
   interleaving (the real driver's predictor is ~3× the replica locally: 15.8 vs 5.2 µs).

## 2. Roofline arithmetic (Rome, 2 GB/s/core conservative)

- **Fused design floor** (operators streamed once/elem/step): 28–36 KB/elem → **14–18 µs** →
  the ≤40 µs target has 2.2–2.8× headroom. **Level-fusion is load-bearing**: the non-fused
  floor is 45–49 µs — above the target. Arithmetic-intensity check: fused AI ≈ 5.2 FLOP/B →
  bandwidth supports 10.4 GFLOP/s/core (compute-bound, good); non-fused AI 1.65 → cap
  3.3 GFLOP/s (bandwidth-bound, target dead).
- **Rate requirement:** 40 µs at 150–160 kFLOP needs **3.75–4.0 GFLOP/s/core sustained**
  (~11 % of Rome peak; ~60 % of SeisSol's measured 6.4 on the same chip). Aggressive but not
  absurd; any FLOP-count reduction toward SeisSol's 92 k relaxes it proportionally
  (92 k/40 µs = 2.3). **The 3–3.5× fallback needs only 1.8–2.3 GFLOP/s/core — near-certain.**
- **Face-stage model** (post-fold): ~42 kFLOP/elem → 11.3 µs at 3.75 GFLOP/s; predictor+
  volume 29–32 µs; **sum 41–43 µs at uniform 3.75, ~37 µs if the GEMM-heavy predictor
  sustains ~5** → 40 µs is achievable-but-boundary.
- **Verdict: ≤40 µs is a defensible stretch target, not a safe promise; the deciding numbers
  (true HW FLOPs/update on Rome, sustainable GFLOP/s at full occupancy) are measurable in
  Phase 0 before any kernel code is written.**

## 3. End-to-end composition (the R-401 table; s/sim-s, TPV104-200m cfl=0.5)

| scenario | MFEM | vs SeisSol as-run (239) | vs async-IO SeisSol (~125) |
|---|---:|---:|---:|
| today GTS / LTS | 6465 / 4356 | 27.1× / 18.2× | — / 34.8× |
| kernel-only @6× | 2677 | **11.2× — comm caps it** | 21× |
| kernel-only @3.5× | 2917 | 12.2× (barely worse — proof comm dominates) | — |
| comm-only (/4 + overlap) | 2015–2600 | 8.4–10.9× | — |
| **both @6× + /4** | 585–921 | **2.4–3.9×** | 4.7–7.4× |
| both @3.5× fallback | 585–1161 | 2.4–4.9× | — |
| + p-drop (49 % cells p1) | ~585–760 | ~3.2× | compute-vs-compute 1.57× |

**Neither program alone approaches SeisSol. The composed end state (2.4–3.9× as-run) is
arithmetically real — and it is comm-bound (585 > 336), so the comm program is not optional.**

## 4. Comm pre-scope highlights (full detail → the comm plan document)

- Complete exchange inventory per sync at Nc=6: 63 I-exchanges + 62 forecast exchanges (both
  full-halo vdim=9), 0 from the retain-only forecast prep, 0 fault-path exchanges under D-2.
- **Merge ladder 125 → 64 → 32** (per-tick hoisting; disjoint element sets ⇒ superposition
  valid; 32 = theoretical minimum since cluster 0 predicts+corrects every tick). **Bitwise
  identical by construction** — same values, same read indices.
- Overlap: MFEM's `ExchangeFaceNbrData` is monolithic BUT the ldof tables are public → a
  ~80-line `NbrExchangerSplit` in `dynamic/` (Begin/End, persistent requests) with the entire
  predict-phase friction + correct-phase volume/interior work as the hiding window.
- **Payoff: both levers → exposed Waitall 2341 → ~150–400 s/sim-s (53 % → 7–16 % of wall),
  wall ~2230–2480 → ~1.8–2.0× LTS-leg speedup.**
- **Hard caveat: an unknown fraction of the 2341 s is late-sender skew (fault-rank friction
  imbalance), untouched by comm changes** — Phase 0 of the comm plan decomposes wire-vs-skew;
  if skew > ~40 %, the fault-weighted partition (METIS ncon=2, existing v3 plan) becomes
  co-requisite.

## 5. What this evidence changes in PLAN rev 3

1. Phase 2's mechanism: `BatchedLinAlg` NATIVE **out**; hand-tiled fused loops / small-dgemm
   **in**; the case re-founded on **traffic reduction** (8×) with a mandatory Rome-side spike
   (same bench file, run on Expanse) as the go/no-go before implementation.
2. The end-to-end composition table goes into the plan Summary (R-401) — the plan alone does
   not approach SeisSol; kernel+comm does.
3. Sequencing re-argued from the numbers (R-403): the comm merge (bitwise-identical, largest
   single lever, no unmeasured-rate dependence) is promoted to co-first-class; comm plan
   document now exists (`document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md`).
4. SeisSol reference re-measurement (async IO) + pinning added to Phase 4 (R-404).
5. Spec fixes: fault-face IntegrationRule branch (`FaultFaceQuadDegree()`), conforming-mesh
   guard, ragged-tile policy, 64-vs-128 ranks/node layout A/B, certified-flag-line pinning.
