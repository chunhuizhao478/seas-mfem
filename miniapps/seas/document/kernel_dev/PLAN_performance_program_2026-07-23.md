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
| **A0** | **λ=1 run (`lts_wiggle="off"`, one deck line, no code).** Also the decisive experiment for the whole track — see below. | exchange rounds/sim-s falls 1.587×; report the change in exposed wait |
| **A1** | merge per-correct exchanges into per-tick rounds (125 → ~32/sync) | bitwise identical output; exposed wait falls |
| **A2** | one wait phase per tick | bitwise identical; wait falls further |
| **A3** | split-post overlap (`NbrExchangerSplit`) | bitwise identical; **plus a liveness gate** |
| **A4** | sparse payload | only if A0/A1 show bytes, not round count, binding |

**A0 is the highest-value action in this document** and it is one config line. It simultaneously (a)
measures the wiggle lever and (b) tests the assumption the *entire* track rests on — **whether exposed
wait scales with exchange-round count at all.** No run has ever varied the round rate; the comm plan's
payoff model was fitted to a hypothetical. If wait does not scale with rounds, A1–A3 are worth far
less than believed and we learn it in one run instead of after weeks of engineering. λ=1 is also
SeisSol's own configuration on this benchmark, so it is not an exotic setting.

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

1. **Whether exposed wait scales with exchange-round rate.** → A0 answers it.
2. **The LTS stage split.** Every Track B per-stage number is a GTS proxy. → instrumentation + one short leg.
3. **MFEM's FLOP count.** → one counter leg.
4. **The comm skew's real size.** The 24 % skew proxy I previously reported as clearing a gate was
   measured with `--tfinal 0.5` while nucleation is at t=1.0 — it sampled a window where fault work is
   near zero *by construction*, so it cannot bound fault-related skew. Treat that gate as **not cleared**.
5. **Which rank is actually on the critical path.** The fault-imbalance argument assumes the
   max-friction rank is the straggler; nobody has checked. The shared-face stage's spread is *larger*
   in wall-seconds (225.9 s vs 196.3 s excess) and is the better suspect.
6. **Whether any of this transfers to production meshes** (34.9–52.3 M tets vs 2.46 M here).

Fault/friction is **not closed**: it is out of scope as a *kernel* target (0.65 % mean) and open as a
*balance* question inside Track A, pending items 4 and 5.

## Appendix — projections (NOT load-bearing; wrong here costs nothing)

Recorded so intent is legible, and quarantined so no correction forces a re-issue. Only the predictor
factor has any benchmark behind it; the face and volume factors are assumptions, and the predictor's
4× is a point estimate chosen from a 2.3–8.3× contended range.

- Track A, if wait scales with rounds: 2341 → ~300 s/sim-s.
- Track B, all phases: ~1016 s/sim-s off the compute term.
- B3 face-cache-for-LTS: ~1.08–1.11× end-to-end, hard ceiling 1.12×; explicitly a factor-of-two error bar.
- Composed: ~3.8–5.6× behind SeisSol as-run. **Parity is not on the table** at matched order; the
  only route past this band is the far-field order drop, which trades accuracy and is a separate proposal.
