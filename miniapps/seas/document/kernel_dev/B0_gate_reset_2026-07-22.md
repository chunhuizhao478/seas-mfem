# B0 gate-reset memo (kernel plan Phase-0 exit artifact) — **FINAL**

**Status: FINAL.** Source: Expanse job **52379131**, COMPLETED, 04:06:55 (14,815 s), all 5 legs rc=0.
Pinned protocol: TPV104-200m (2,464,689 tets), 2 nodes × 128 ranks, GTS leg
(`tpv104_200m_lts_off.toml`), flags `--deriv-cache --shared-ck-recursion --face-cache
--cfl-dg-safety 1.0`, cfl 0.5 (deck), git `a993a1b`, mesh md5 `e88ae223…`.

## One-line verdict

**B0 = 187.5 µs·core** (face-cache activation delivered **1.31×** vs the 245.8 reference — dead in
the predicted 1.2–1.4× band). **Two premise-overturning findings:** the stage split on target
hardware is **inverted** vs the local profile (predictor **56.2 %** ≫ face **23.5 %**), and the
kernel bench **reverses** its local verdict under Rome contention (**B = 2.3–6.3×**, D = 2.8–8.3×,
vs 1.25×/1.30× locally) ⇒ **Phase-2 GO**, and Phase 2 is now the *primary* phase, not Phase 1.
Fault deferral **confirmed** (0.65 % even through rupture). Comm skew **24 % < 40 %** ⇒ comm-plan
payoff holds. **Revised program target: ~2.3× vs B0 ≈ 81 µs·core ≈ 3.0× vs 245.8 — the declared
FALLBACK band. The ≤40 µs / 6× stretch is NOT reachable on these measurements.**

## Leg 1 — BLAS

`MFEM_USE_LAPACK = YES`; `dgemm_` present as a defined (`T`) symbol — PETSc's LP64 openblas is
**statically** linked from `extern/petsc/install/lib` (hence absent from `ldd`). **The 245.8 µs
was already on the optimized-dgemm LU path** ⇒ no BLAS win available; the `build_expanse.sh`
lib64 fix remains as defensive hardening only. *(Cosmetic: the leg's canned "no BLAS in ldd →
internal kernels" line mis-reads a static link — sbatch text fixed.)*

## Leg 2 — B0 (the denominator)

`per-update µs·core = wall × ranks × 1e6 / (elems × Δsteps)`; whole-run:
`8873.95 s × 256 × 1e6 / (2,464,689 × 4917)` = **187.5 µs·core** (1.805 s/step).

| | value |
|---|---|
| **B0** | **187.5 µs·core** |
| vs Phase-5 245.8 (face-cache OFF) | **1.31× activation gain** ✓ predicted band |
| windows extracted | pre-rupture steps 0→1600 (t≤0.586); rupture steps 3600→4916 (t 1.318→1.8, V_max 13.7→11.0 m/s) |
| per-window µs·core | **not separable** — per-leg wall not captured (sbatch gap, fixed). Bounded: friction is 0.65 % of the step, so B0_pre ≈ B0_rup ≈ 187.5 to <1 % |

**Stage split (Caliper avg/rank, `step` = 8874 s = 99.99 %):**

| stage | Time % | note |
|---|---:|---|
| **predictor** `ComputeADERSubStepStatesAndIntegral` | **56.23** | ← the dominant target |
|  └ `ApplySpatialDerivative` | 27.14 | |
| corrector `AdvanceADER` | 34.40 | |
|  └ face interior `ComputeADERFaceFluxRHS` | 12.96 | (`…_CachedInterior_` active; perf self 7.2 %) |
|  └ face shared `ComputeADERSharedFaceFluxRHS` | 10.56 | Min 585 / Max 1163 = **2.0× rank spread** even in GTS |
|  └ volume `ComputeVolumeRHS` | 7.88 | |
| **fault** `friction_substep` | **0.65** | Min 2.7 / Max 254 s — huge spread, tiny mean |

**⚠ PREMISE OVERTURNED:** the local macOS profile (face 59 %, predictor 27 %, face-cache OFF) is
**inverted** on target hardware with the cache on — **predictor is 2.4× the face stage.** Cause:
`--face-cache` removed the interior-face LU/CalcShape tax, and the predictor's whole-vector sweeps
are exactly what a 128-rank memory-starved node punishes. *The plan's Phase-1-first ordering was
built on the wrong split.*

perf: symbols resolved (`ApplySpatialDerivative` 23.8 % self, `Vector::Add` **18.9 %** self ← the
AXPY sweeps, `Vector::operator=` 9.3 %, `ComputeVolumeRHS` 7.3 %, `…FaceFluxRHS_CachedInterior_`
7.1 %, `ApplyElementJacobian_` 6.8 %, `DenseMatrix::Mult` 2.8 %, `GodunovFlux::Interior` 2.2 %).
`perf report` needs a perf-enabled host — `b0_gts.perf.data` is retained for a deeper pass.
**`Vector::Add` + `operator=` = 28 % of self-time in pure vector traffic** — direct confirmation of
the whole-vector-sweep diagnosis Phase 2 targets.

## Leg 3 — layout A/B (128 vs 64 ranks/node)

Ran clean (1366 steps to t=0.5, rc=0, NUMA-spread) **but per-update is NOT computable — the leg's
wall time was never printed** (sbatch gap; fixed for re-runs). The memory-starvation question is
answered indirectly and more strongly by Leg 4 (below), so this is not blocking.

## Leg 4 — kernel bench on Rome ⇒ **PHASE-2 GO**

Single-core (uncontended): A 32.98 · B 14.08 (**2.34×**) · C 40.67 (0.81×) · D 11.78 (**2.80×**) µs/elem.
Full-occupancy (256 copies, memory-contended), representative ranks:

| | A (current) | B (tiled-fused) | C (BatchedLinAlg) | D (dgemm) |
|---|---:|---:|---:|---:|
| µs/elem | 33 → **79** | **12.5–14.4 (flat!)** | 21–55 | 8.1–13.8 |
| speedup vs A | 1.00 | **2.26 – 6.34×** | 0.69–3.13× | **2.80 – 8.25×** |

**A degrades 33→79 µs under contention while B stays flat at ~13 µs.** That is the traffic
hypothesis confirmed on the target machine: A streams ~8× the DRAM bytes and collapses when 128
ranks share the memory controllers; B is cache-resident and immune. **GO/NO-GO: R = 2.3–6.3 ≥ 2 →
Phase 2 proceeds** (no full-fusion variant needed; the partial-fusion lower bound already clears).

**This reverses the local verdict** (B 1.25×, C 0.45×) exactly as flagged — Apple-silicon bandwidth
masked the effect. C stays eliminated as the *mechanism* (worst of B/C/D everywhere) though its
local "2.2× slower" was itself hardware-specific (it beats A under heavy contention). **D (dgemm)
is the best variant** and is now a first-class mechanism candidate, enabled by Leg 1's confirmed
static LAPACK.

## Leg 5 — comm wire-vs-skew

`MPI_Waitall` Min 877 / **Avg 1157** / Max 1983 s = **52.6 % of wall** (Phase-5: 52.5 % — consistent);
Max/Avg = **1.71** (Phase-5: 1.72). Injection negligible (`Isend` 0.08 %, `Irecv` 0.02 %).
Skew proxy `(Avg−Min)/Avg` = **24 % < 40 % threshold** ⇒ **the comm plan's ~1.8–2.0× payoff holds**;
~76 % is wire/sync time that merging + overlap directly attack. The fault-weighted partition stays
a secondary lever, not a co-requisite.

## Gate resets (now concrete)

| plan gate | was (provisional) | **RESET** |
|---|---|---|
| **Program target** | ≤40 µs·core (≥6× vs 245.8) | **~81 µs·core (≈2.3× vs B0, ≈3.0× vs 245.8)** — the declared fallback band; ≤40 µs unreachable |
| **Phase priority** | P1 (face 59 %) then P2 | **P2 FIRST (predictor 56.2 %)**; P1 demoted (face 23.5 %, LU already removed by face-cache) |
| Phase-2 go/no-go | contended B ≥ 2× | **PASSED (2.3–6.3×)** — proceed; mechanism = hand-tiled fused (B) and/or per-element dgemm (D) |
| Phase-2 Expanse gate | ≥1.5–2× vs B0 on pred+vol | **≥1.6× vs B0 on the predictor stage** (from 56.23 % at ~4× ⇒ 1.73× overall) |
| Phase-1 gate | face share ≤0.5× B0 share, ≥1.2× per-update | **face share ≤0.6× of 23.5 %, ≥1.10× per-update** (headroom shrank — the cache took the big bite) |
| fault work | deferred pending rupture data | **DEFER CONFIRMED — 0.65 % through rupture** (even at a 4× bulk speedup it is ~2.5 %) |
| comm plan payoff | ~1.8–2.0×, skew caveat | **HOLDS** — skew 24 % < 40 % |

**Projections from the measured split** (predictor at bench factors, face 1.5×, volume 2×):
P2 alone 1.73–1.82× · +P1 2.00× · +volume **2.31× vs B0 = 81 µs·core**.

**Composed end state (updated):** MFEM-LTS with face-cache ≈ 1538 compute + 2341 comm = 3879 s/sim-s;
kernel program → 669 + 2341 = 3010 (12.6× behind SeisSol as-run 239); **+ comm program → ~669–969
⇒ 2.8–4.1× of SeisSol as-run** (prior estimate 2.4–3.9× — holds, slightly softer).

## Harness gaps found by this run (fixed for re-runs)

1. **Per-leg wall time never printed** → Leg 3 unusable, Leg 2 windows not separable. *(fixed: per-leg `SECONDS` timing)*
2. Leg-1 interpretation line mis-reads a *static* BLAS link as "no BLAS". *(fixed)*
3. Leg-2 banner still said "EndTime 1.5s … [1.0,1.5]" after the window change to 1.8/[1.4,1.8]. *(fixed)*
4. Spurious "bench build/run failed" after Leg 4 actually succeeded (SIGPIPE from `head` in the pipeline). *(fixed)*
5. `perf report` unavailable on compute nodes — `perf.data` retained; needs a perf-enabled host for the full call tree.
