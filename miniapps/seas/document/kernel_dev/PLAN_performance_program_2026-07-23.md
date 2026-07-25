# Performance program — settled direction

**Date:** 2026-07-23 · **Branch:** `safs-v4_0_0-alt-case1-mfem-speed` · **Status: DIRECTION FIXED.**
Supersedes `PLAN_v2_performance_program_2026-07-22.md` and `PLAN_ader_kernel_efficiency_2026-07-21.md`.

> **Why this document is structured differently.** The two previous plans were rewritten three times
> in two days. Every rewrite chased a *number*, not a direction — the direction has not changed since
> v2, and a 283-agent adversarial review (`REVIEW.md`, 28 findings, verdict FAIL) stated explicitly
> that "the directional conclusion — communication before kernels — survives every finding and does
> not need re-arguing." The plans kept breaking because their **bodies were load-bearing on projected
> payoffs**. This one is not. The body contains only measured quantities and gates expressed as
> *measured deltas*. Every projection lives in the appendix, where being wrong costs nothing and
> forces no re-issue.

## The direction (settled — not revisited without a measurement that contradicts it)

1. **Communication first.**
2. **Kernels second**, predictor before faces, cheapest lever first.
3. **Face-cache-for-LTS** inside Track B, behind a hard correctness guard.

Three independent reasons, none of which depends on a projection:

- **It is the largest single line item.** Exposed `MPI_Waitall` is 52.5 % of the LTS wall (measured,
  job 52344266) and 52.6 % on a different job (measured, job 52379131 Leg 5).
- **It is the only bitwise-safe track.** Merging exchange rounds is identity-preserving by
  construction; kernel work needs tolerance gates.
- **It is what makes LTS worth having at all.** This is the strongest argument and it is entirely
  measured: LTS cuts compute 3.21× (87 % of its 3.68× ceiling) yet leads GTS by only **1.13×**,
  because the wait eats the win. We have already paid for LTS; comm is how we collect.

## A0 RESULT (2026-07-24, jobs 52416603 + 52422891) — first measured milestone

Full write-up: `document/comm_dev/RESULTS_a0_wiggle_2026-07-24.md`.

- **Track A GO, on measurement.** Exposed wait scales *super*-proportionally with the exchange-round
  rate: rounds ×1.5871 (predicted 1.5873 — 4-decimal agreement) → wait ×1.807.
- **`lts_wiggle="off"` is a measured 1.370× end-to-end speedup for one config line**, stable through
  rupture (V_max 17.83 vs 17.69 m/s = 0.80 %, no NaN). 4277 → **3123 s/sim-s**; 17.9× → **13.1×**
  behind SeisSol. Adopt as default for np≥256.
- **The wait is NOT bandwidth.** `Isend`+`Irecv` = 6.9 s of 4523 s = **0.15 %**. GTS pushes ~2×
  the rounds/sim-s of LTS yet waits ~0.00002 % of wall. The LTS wait is **ranks arriving at sync
  points at different times**. A1–A3 remain right, but the mechanism is *fewer exposure events*, not
  fewer bytes ⇒ **A4 de-prioritised** (it attacks 0.15 %).
- **Known-unknown #4 RESOLVED, and it points away from fault imbalance:** skew is **flat through
  rupture** (24.1 → 24.6 %, Max/Avg 1.73 → 1.72), so the imbalance is **structural** (cluster/partition
  layout), not rupture-driven. This also replaces the earlier 24 % proxy that was withdrawn for having
  been measured pre-nucleation.
- **First measured LTS stage split** (commit `884349a` made it possible): true compute is 47.1 % of
  step; **the predictor is 39.6–46.6 % of compute, not the 56.2 % the plan assumed** from the GTS
  proxy. Track B gates must be re-derived against this. The seam corrector's *compute* half (12.0 %
  of step) is a target that appears in no phase yet.
- **Independent confirmation:** measured true compute at λ=0.63 = **2015 s/sim-s**, exactly the value
  the budget below carried as a *derived* term.
- **New item A5:** `lts_clustering.cpp:compute_cost` minimises element updates with **no communication
  term** — that is why it chose a λ costing 1.37×. Hardcoding λ=1 is right here but may be wrong at
  another rank count/mesh; add a comm term to the objective.

## Where we are (measured only)

| | s/sim-s | provenance |
|---|---:|---|
| MFEM GTS, face-cache off | 6465 | job 52344266 |
| MFEM GTS, face-cache on | 4930 | job 52379131 Leg 2 |
| **MFEM LTS** | **4356** | job 52344266 |
| SeisSol o4 LTS, as-run | 239 | job 52365078 |
| SeisSol, compute phase only | 111 | job 52365078 |

**MFEM-LTS is 18.2× behind SeisSol as-run** (4356/239). Within the 4356, ~2341 s/sim-s is exposed
wait and ~2015 s/sim-s is compute — **both are derived, not directly measured**: 2015 is inferred as
6465/3.208 and 2341 is the residual. `--face-cache` does **not** apply here: `use_face_cache_` has
one functional read (`wave_operator.inl:5475`) whose sole caller is the GTS corrector, so the 1.31×
it gave GTS cannot be transferred to LTS.

## Work list

Each item states what it does and the **measured** gate that decides whether it stays. No item is
justified here by a projected payoff.

### Track A — communication

| # | work | gate (measured) |
|---|---|---|
| ~~A0~~ | **DONE 2026-07-24.** λ=1 run. | **PASSED:** rounds ×1.5871, wait ×1.807, **1.370× end-to-end**, stable through rupture. Track A GO. |
| ~~A1~~ | merge 125→64 rounds/sync | **DONE, NO GAIN (job 52472765).** Rounds fell 1.954× as designed; wait 0.971×, wall 0.988×. Bitwise-identical at np=256. Code kept, flag OFF. |
| ~~A2~~ | one wait phase per tick (64→32) | **DEVALUED by A1** — same operation A1 showed does nothing (merging rounds *within* a tick removes no rendezvous). Needs new justification or drop. |
| ~~A3~~ | split-post overlap | **DEVALUED** — hides wire time, which is 0.15 % of the wait. Nothing to hide. |
| ~~A4~~ | sparse payload | **DE-PRIORITISED by A0** — wire time is 0.15 % of the wait; bytes are not binding. |
| **A5** | comm term in `lts_clustering.cpp:compute_cost` | still valid — it generalises A0, the one lever that MEASURED a gain |
| **A6** | **rank rebalancing (NEW, now the main Track-A item)** | the wait is imbalance at tick boundaries: skew 16.5 %, **Max/Avg 2.18**. Only rebalancing or fewer tick boundaries touch it. Unsized. |

> **A1 RESULT (job 52472765) — the merge is correct and buys NOTHING.** Rounds 12,812 → 6,556
> (1.954×, exactly as designed); wait 764 → 787 s (**0.971×**); wall 1898 → 1921 s (**0.988×**).
> Bitwise identical at np=256. Per-round wait *doubled* — each surviving round absorbs the wait the
> removed ones carried, so **the wait was never in the rounds.** A0 worked because it cut the number
> of **syncs** (271→171 = fewer rendezvous in time); A1 kept every tick boundary and merged
> back-to-back exchanges *at* those boundaries, removing no rendezvous. **Corrected model: exposed
> wait scales with HOW OFTEN ranks must meet, not with how many messages they send when they do.**
> ⇒ A2/A3 devalued by the same mechanism; the Track-A projections have no measured support; the
> remaining wait is rank imbalance (Max/Avg 2.18) and only rebalancing touches it (new item A6).
> Full write-up: `document/comm_dev/RESULTS_a1_merge_2026-07-24.md`.

**A0 is DONE (2026-07-24) and it PASSED** — see the A0 RESULT section above. It confirmed the
assumption the entire track rested on: exposed wait does scale with the exchange-round rate, and
super-proportionally (1.807× per 1.587×). That converts A1–A3 from a fitted model into
measurement-backed work. It also delivered a 1.370× end-to-end speedup for one config line, and
redirected the *mechanism*: the wait is synchronisation exposure, not bandwidth (wire time is
0.15 % of it), so A1–A3 pay off by removing sync points rather than bytes.

**A3 carries a liveness gate, not just a correctness gate.** Bitwise identity cannot detect this
code's actual recorded failure mode for exchange changes — the R-1600 unmatched-collective hang
(`wave_operator.inl:3661-3672`). Deadlock produces no wrong bits; it produces no bits.

### Track B — kernels

| # | work | gate (measured) |
|---|---|---|
| **B1** | fuse the predictor's whole-vector accumulate/zero sweeps | ≥1.2× on the predictor stage; ≤1e-12 parity |
| **B2** | cache-tiled fused predictor | ≥1.6× on the predictor stage |
| **B3** | face-cache extended to the LTS corrector (see guard below) | ≥1.05× end-to-end; ≤1e-12 parity |
| **B4** | volume stage tiling | ≥1.5× on that stage |
| **B5** | face-interior tables — **re-justify before building** | drop if <5 % of the step |

**B3's guard is mandatory and is the highest-risk detail in this document.** The cache builder is
cluster-blind (`face_geom_cache.hpp:101-108` excludes only shared/boundary/fault faces), so every
`ConsumerFine` and `ProviderCoarseSkip` cluster-seam face **is** in the cache and a naive lookup
copied from the GTS site will *hit* on it and be silently wrong — wrong state vector, wrong
destination buffer, and a dropped contribution with no diagnostic. Gate on
`role == FaceRole::IntraClusterGTS`, place the lookup between the role skip and
`GetFaceElementTransformations`, and use `MFEM_VERIFY` (not `MFEM_ASSERT`, which compiles out in
production builds). **No existing test would catch a mis-guard**: the face-cache parity test drives
the GTS path only, and the multi-cluster LTS test is degenerate by construction. A new seam test is
part of B3, not optional.

**Flag semantics caveat:** `--face-cache` is already passed by in-tree LTS jobs, so extending it
changes what those recorded legs do on a re-run. Either accept that and re-baseline, or add a
separate `--lts-face-cache`.

### Prerequisites

- **Instrumentation (partly landed, commit `884349a`).** The LTS path had no Caliper step region at
  all — the existing one is GTS-gated — so no LTS stage table could ever be produced. Added: a step
  region in both LTS sync loops and a scope on the LTS corrector. **Still missing:** sub-stage scopes
  inside the corrector (interior-face, volume, friction). Until those land, the LTS split resolves
  predictor/corrector/seam but not the corrector's internals.
- **One FLOP-counter leg.** Track B's mechanism (traffic reduction) and the "hand-tuned ceiling"
  argument both assume MFEM is memory-bound. Nobody has ever counted MFEM's FLOPs per element update.
  One `perf stat` leg decides whether B1/B2 can work at all.

## Known unknowns (do not claim these are settled)

Verification of the previous plan's numbers returned **43 quantities that cannot be derived from any
existing artifact.** The ones that would change decisions:

1. ~~Whether exposed wait scales with exchange-round rate.~~ **RESOLVED by A0: yes, super-proportionally (1.807× per 1.587×).**
2. ~~The LTS stage split.~~ **MEASURED by A0** (predictor 39.6–46.6 % of compute, NOT the assumed 56.2 %; seam-corrector compute 12.0 % of step is unclaimed). **Track B gates must be re-derived against it.**
3. **MFEM's FLOP count.** → one counter leg.
4. ~~The comm skew's real size.~~ **RESOLVED by A0's 2 s through-rupture leg: skew is FLAT (24.1 →
   24.6 %, Max/Avg 1.73 → 1.72), so the imbalance is STRUCTURAL, not fault-driven.** This supersedes
   the withdrawn pre-nucleation proxy.
5. **Which rank is actually on the critical path.** The fault-imbalance argument assumes the
   max-friction rank is the straggler; nobody has checked. The shared-face stage's spread is *larger*
   in wall-seconds (225.9 s vs 196.3 s excess) and is the better suspect.
6. **Whether any of this transfers to production meshes** (34.9–52.3 M tets vs 2.46 M here).

Fault/friction is **not closed**: it is out of scope as a *kernel* target (0.65 % mean) and open as a
*balance* question inside Track A, pending items 4 and 5.

## Appendix — projections (NOT load-bearing; wrong here costs nothing)

### Speedup inventory (2026-07-24) — sized against the MEASURED post-A0 stage split

Baseline 4277 s/sim-s = 17.9× SeisSol as-run. Sizes come from the measured LTS split
(`RESULTS_a0_wiggle_2026-07-24.md`), gains are projections except where marked.

| # | item | targets | saves | running | vs SeisSol | confidence |
|---|---|---|---:|---:|---:|---|
| **A0** | `lts_wiggle="off"` | sync count | **1154** | **3123** | **13.1×** | **MEASURED — banked** |
| **A1** | merge 125→**64** rounds/sync | wait (1251) | 611 | 2512 | 10.5× | projected from A0's law |
| **A2** | single-round tick, 64→**32** | residual wait | 320 | 2192 | 9.2× | same law; **NOT a minor item** |
| A3 | split-post overlap | residual wait | 72 | 2120 | 8.9× | low; anti-synergistic with B |
| A5 | comm term in objective | — | 0 | 2078 | 8.7× | generalises A0; no new gain at np=256 |
| **B1+B2** | predictor fusion → tiling | predictor (871) | 653 | 1425 | 6.0× | bench-backed, not in-solver; B1,B2 do NOT add |
| B3 | face-cache → LTS corrector | corrector-rest | 93 | 1332 | 5.6× | ±2× error bar |
| B4 | volume tiling | corrector-rest | 143 | 1189 | 5.0× | ASSUMED — no bench of this stage anywhere |
| — | seam-corrector compute | seam (246) | 82 | 1107 | 4.6× | ASSUMED — in no phase yet |

**Endpoint if everything lands: ~1107 s/sim-s ≈ 4.6× SeisSol.** Read it as "~5×, maybe".

Three structural facts this table encodes:
1. **A1 + A2 + B1/B2 carry ~85 % of the remaining gain** (611 + 320 + 653 of ~1900). *Correction
   (2026-07-24): an earlier version of this table credited A1 with 125→32 and left A2 at 42. The comm
   plan stages it 125→**64** (Phase 1 = A1) then 64→**32** (Phase 2 = A2), so the pair splits
   611/320 — **A2 is the second-largest comm item, not a rounding error.** The A1+A2 endpoint is
   unchanged at 2192.*
2. **Neither track suffices alone.** Perfect comm / untouched kernels floors at **7.8×**; perfect
   kernels / untouched comm floors at **5.2×**. Under ~5× needs both.
3. **Confidence decreases down the table** — A0 measured, A1 one extrapolated data point, B4 and the
   seam item are unbenchmarked assumptions.

**Not in this table, and it should be:** the structural imbalance A0 found. 19–24 % of the wait is
rank spread that merging cannot remove, and it is flat through rupture. If that floor binds, A1 lands
nearer 1.3× than 1.42× and the whole A column shrinks. Only rebalancing touches it; nobody has sized
it. **Parity with SeisSol is not in this table** — the honest endpoint is ~5×.


Recorded so intent is legible, and quarantined so no correction forces a re-issue. Only the predictor
factor has any benchmark behind it; the face and volume factors are assumptions, and the predictor's
4× is a point estimate chosen from a 2.3–8.3× contended range.

- Track A, if wait scales with rounds: 2341 → ~300 s/sim-s.
- Track B, all phases: ~1016 s/sim-s off the compute term.
- B3 face-cache-for-LTS: ~1.08–1.11× end-to-end, hard ceiling 1.12×; explicitly a factor-of-two error bar.
- Composed: ~3.8–5.6× behind SeisSol as-run. **Parity is not on the table** at matched order; the
  only route past this band is the far-field order drop, which trades accuracy and is a separate proposal.
