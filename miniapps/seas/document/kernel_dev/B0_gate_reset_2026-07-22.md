# B0 gate-reset memo (kernel plan Phase-0 exit artifact)

**Status: TEMPLATE — awaiting the Phase-0 Expanse job** (`jobs/kernel_dev/run_p0_baseline_expanse.sbatch`).
Fill every `‹…›` from the job output, then flip Status to FINAL and every downstream gate in
`PLAN_ader_kernel_efficiency_2026-07-21.md` becomes concrete. Until then all plan gate numbers are provisional.

Pinned protocol (so "vs B0" is unambiguous): TPV104-200m, 2 nodes × 128 ranks, GTS leg
(`tpv104_200m_lts_off.toml`), flags `--deriv-cache --shared-ck-recursion --face-cache`, cfl 0.5,
binary hash `‹git+build hash›`, mesh md5 `e88ae223…`.

## Leg 1 — BLAS

- Deployed `MFEM_USE_LAPACK` = `‹YES/NO›`  ·  `ldd` BLAS = `‹openblas path / none›`
- If NO: rebuilt after the `build_expanse.sh` lib64 fix → `‹YES›`; physics sanity vs 52344266 GTS
  (V_max(t)) = `‹PASS/FAIL, max Δ›`. **B0 is measured on the LAPACK=YES binary.**

## Leg 2 — B0 per-update µs·core (the denominator)

`per-update µs·core = wall_window[s] × ranks_in_window × 1e6 / (elems(2,464,689) × Δsteps_window)`
where **`ranks_in_window` = the ranks actually stepping the mesh** (256 for Leg 2/5, **128 for Leg 3** —
review F1; hardcoding 256 for the 64-ranks/node leg would double its number and self-corrupt the
memory-vs-kernel conclusion). Anchor check: Phase-5 GTS 2732 steps to t=1.0 in ~6465 s stepping →
`6465×256×1e6/(2.46e6×2732) = 245.8 µs·core` ✓ (matches the 245.9 reference; use STEPPING wall, not
wall-incl-init, on both sides). `Δsteps_window` = (end step index − start step index) from the two
boundary `step N/M` log lines (the driver prints every 100th step, so each is ±100 — report both
endpoints, not just the tail).

| window | start→end step | wall (s) | **µs·core/update** |
|---|---:|---:|---:|
| pre-rupture [0, 0.5] | ‹›→‹› | ‹› | **‹B0_pre›** |
| rupture [1.4, 1.8] (trails the peak) | ‹›→‹› | ‹› | **‹B0_rup›** |

Sub-function split (perf self-time, rank 0), % of step: LU/`LUFactors::Solve` ‹›, `CalcShape`/geom ‹›,
`GodunovFlux::Interior` (out-of-line) ‹›, `ApplySpatialDerivative` ‹›, `ComputeVolumeRHS` ‹›,
`ApplyMassInverse` ‹›, fault/friction **‹B0_fault_rup›** (rupture window — the fault-work decision hinges
on this; also note the PEAK-step share, not just the window mean — review F2).
Face-stage subphase share of the step = **‹B0_face_share›** (Phase-1 gate denominator).
*perf caveat (F4):* out-of-line hotspots (LU, CalcShape, GodunovFlux::Interior, ComputeVolumeRHS,
ApplyMassInverse) attribute cleanly by self-time; intra-`ProcessADERFaceToRHS_` inlined blocks
(BuildFrame/ApplySplitFlux) may fold into the parent — the local macOS `sample` profile already gave
that finer split, so cross-check against it rather than treating perf as the sole source.

## Leg 3 — layout A/B (128 vs 64 ranks/node)

- 128 ranks/node µs·core = ‹= B0_pre› · 64 ranks/node µs·core = ‹›  ·  ratio = `‹L›`
- Interpretation: L ≫ 1 ⇒ the 245.9 figure is per-core memory starvation ⇒ the Phase-2 traffic
  case strengthens and an OpenMP-hybrid layout (like SeisSol's 8×15) becomes a candidate lever.

## Leg 4 — kernel bench spike (Phase-2 go/no-go)

Single-core µs/elem: A ‹5.3 local› B ‹› C ‹elim 0.45×› D ‹›. **Full-occupancy (contended):**
A ‹› B ‹› → **contended B-vs-A = ‹R›**. Slowdown vs single-core: A ‹×› B ‹×`.

**GO/NO-GO:** kernel Phase 2 proceeds iff `R ≥ 2`. **Caveat (review F5):** bench variant B is
*partial* fusion (E=64 tile, 3 levels fused) — a LOWER BOUND on the full-level-fusion design the
roofline justifies (non-fused floor 45–49 µs vs fused 14–18 µs ≈ 2.7–3.2×). So if contended B lands
in `[1.6, 2.0)`, do NOT kill Phase 2 — first build a FULL-fusion bench variant and re-measure; only
`R < 1.6` (partial-B) is a clean no-go, deprioritizing Phase 2 below comm Phase 3 (overlap) per the
revisit trigger — the composed end state is comm-bound anyway.

## Leg 5 — comm wire-vs-skew

`MPI_Waitall` Max ‹› / Avg ‹› = `‹skew ratio›`. If the skew fraction of the 2341 s/sim-s exceeds
~40 %, the comm plan's payoff degrades toward 1.5× and the fault-weighted partition (METIS ncon=2,
existing v3 plan) is promoted to co-requisite.

## Derived gate resets (fill after the numbers above)

| plan gate | provisional | reset-from-B0 |
|---|---|---|
| kernel program target ≤40 µs·core (≥6× vs 245.9) | asserted | `‹6 × B0_pre/245.9›` × vs B0, OR the 3–3.5× fallback if Leg 4 `R<2` |
| Phase-1 face-stage share drop ≤0.5× B0 share | provisional | B0_face_share = ‹› → target ≤ ‹0.5×› |
| Phase-2 pred+vol Expanse ≥1.5–2× vs B0 | provisional | roofline-conditioned; from Leg 4 R + Leg 3 L |
| fault-work: optimize or defer | deferred | decide from B0_fault_rup (Leg 2 rupture window) |

## One-line verdict (fill last)

`‹e.g. "B0_pre = 210 µs·core; LAPACK was NO→YES (+X%); Leg-4 R=2.6 → Phase 2 GO; skew 25% → comm
plan payoff intact; kernel target reset to 5.2× vs 245.9."›`
