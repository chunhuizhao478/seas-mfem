> # ⛔ SUPERSEDED (2026-07-23) — see `PLAN_performance_program_2026-07-23.md`
>
> Failed adversarial review (`REVIEW.md`, 28 findings, 5 CRITICAL). Its headline budget was wrong by
> construction: it divided the measured LTS wall by a face-cache gain that only applies to the GTS
> corrector. Its **direction — communication before kernels — survived every finding** and carries
> forward unchanged; the successor restates it at its true magnitude and, unlike this document, keeps
> no projected payoff in its body. Kept as review history.

# Performance program v2 — post-Phase-0 redesign

**Date:** 2026-07-22 · **Branch:** `safs-v4_0_0-alt-case1-mfem-speed` · **Supersedes** the phase
structure, ordering, and targets of `PLAN_ader_kernel_efficiency_2026-07-21.md` (rev 1–4).
**Evidence base (all measured, none assumed):** `B0_gate_reset_2026-07-22.md` (job 52379131),
`EVIDENCE_microbench_roofline_2026-07-22.md`, `document/lts_dev/RESULTS_p5_tpv104_200m_2026-07-21.md`.
**Governs two tracks:** this document (Track B, kernel) and
`document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md` (Track A, communication).

## Summary — read this first

**Where we actually are.** On the benchmark problem our solver takes **3879 seconds of wall time
per simulated second**; SeisSol takes 239 — we are **16× behind**. (That already includes the free
1.31× from switching on a cache we had built but never enabled.)

**What the time is actually spent on** — this is the finding that forces the redesign:

| | s/sim-s | share |
|---|---:|---:|
| **processors waiting on each other** | **2341** | **60 %** |
| computing — of which the predictor is 56 % | 1538 | 40 % |

**The redesign in one line: communication is the primary program, not the kernels.** Fixing the
waiting recovers **2041 s/sim-s**; fixing *every* kernel we have evidence for recovers **776**.
The old plan had kernels first and communication as a deferred companion — backwards by a factor
of 2.6 in payoff, and backwards again in risk (the communication fix is bit-for-bit identical by
construction; the kernel work needs tolerance gates).

**Expected outcome.** Both tracks land us at **~760–1060 s/sim-s = 3.2–4.4× behind SeisSol
as-run** (6–8.5× against its faster async-IO configuration), from 16× today. That is the honest
target; it is the fallback band the previous plan pre-committed to, and the ≤40 µs·core / ≥6×
kernel stretch goal is **withdrawn** as unreachable on the measured stage split.

**Biggest risk.** The communication saving assumes the waiting is mostly *wire and synchronisation*
time that merging and overlapping can remove. Phase 0 measured the skew fraction at **24 %** (below
the 40 % threshold that would have demoted the whole track), so the assumption is supported — but
the residual 24 % is load imbalance that no communication change touches, and it caps Track A.

**What is now closed** (measured, not deferred): fault/friction optimisation (**0.65 %** of the
step even through rupture — permanently out of scope); the BLAS question (already optimal,
statically linked dgemm — no win exists); MFEM's built-in batched linear algebra as the kernel
mechanism (**benchmarked slower than what we have**).

## What changed from v1, and why

| v1 assumed | Phase 0 measured | consequence |
|---|---|---|
| face stage 59 %, predictor 27 % (local laptop profile) | **predictor 56 %, face 23.5 %** on the cluster with the cache on | phase order **inverted**; face work demoted |
| kernels first, comm deferred | comm 60 % of wall, 2.6× the payoff, and bitwise-safe | **comm promoted to primary** |
| batched linear algebra as the mechanism | 2.2× *slower* locally; still worst variant on the cluster | mechanism replaced by hand-tiled fusion / per-element dgemm |
| tiled fusion worth ~1.25× (laptop) | **2.3–6.3× under real memory contention** | Phase-2 GO; the win is data movement, not arithmetic |
| fault deferred pending rupture data | 0.65 % through rupture | closed |
| target ≤40 µs·core (≥6×) | split doesn't support it | **withdrawn**; target = fallback band |

## Verification round (2026-07-22, post-redesign) — what was checked and what it changed

The budget above is load-bearing, so it was audited against the source before either track starts.

**CONFIRMED — the main kernel lever is real.** The concern was that the predictor's 56 % might be
*waiting* nested inside a Caliper region rather than work. It is not: the predictor entry, the LTS
cluster variant, and `ApplySpatialDerivative` contain **zero** MPI or `ExchangeFaceNbrData` call
sites. Every exchange in the wave operator lives in exactly two places — `ComputeADERSharedFaceFluxRHS`
and `ComputeADERClusterSeamFaceFluxRHS` — i.e. precisely the regions this plan assigns to Track A.
That independently confirms the shared-face reassignment above, and identifies the LTS seam
corrector as the second comm site.

**DEFECT FOUND AND FIXED — the two LTS hot paths were uninstrumented.**
`ComputeADERSubStepStatesAndIntegralCluster` (the LTS predictor — *B2's exact target*) and
`ComputeADERClusterSeamFaceFluxRHS` (the LTS seam corrector) had **no Caliper scope at all**. Both
tracks' gates are written against LTS-leg stage shares, so as written **neither track could have
measured its own gate.** Both scopes are now added at function entry, mirroring the placement of
their non-cluster siblings so region reports stay structurally comparable. Zero numerics touched.
*Open gate:* the local build has `MFEM_USE_CALIPER = NO`, so the macro expanded to nothing and the
real `CALI_CXX_MARK_SCOPE` path is **compile-verified only on the next Expanse build**, not locally.

**CAVEAT — the stage split is GTS-measured, the budget is LTS.** Because of the gap above, the
56/23/8 % split could only come from the GTS leg, while the 1538 + 2341 budget is the LTS leg,
which additionally runs a seam corrector that does not exist in GTS. Two consequences: the split
must be **re-measured on the LTS leg** now that instrumentation exists (fold this into A1/B1's
first run, no separate job needed); and the error direction is **conservative** for the predictor —
GTS's face-shared share includes its wait, so the predictor's share of true compute is if anything
slightly *higher* than 56 %. The 776 s/sim-s kernel figure already excludes face-shared entirely.

## Track A — communication (PRIMARY, 2041 s/sim-s)

Detailed phases live in `document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md`; Phase-0
measurements confirm its premises. Summary of the four phases and what Phase 0 settled:

| phase | what | Phase-0 status |
|---|---|---|
| A0 wire-vs-skew split | decompose the 2341 s | **partly done** — skew ≈ 24 % (< 40 % gate) ⇒ track proceeds; a finer per-tick split is still worth capturing during A1 |
| **A1 merge 125→64→32 rounds/sync** | hoist per-correct exchanges to per-tick; disjoint element sets superpose | **ready** — bitwise-identical by construction; the single biggest, safest win |
| A2 single-round tick exchange | one wait phase per tick (theoretical minimum) | ready |
| A3 split-post overlap | `NbrExchangerSplit` on MFEM's public tables; hide wire time behind the tick's ghost-independent compute | ready; depends on A2 |
| A4 sparse payload + GTS retrofits | ship only seam-adjacent elements | only if A0's finer split shows bytes (not skew) binding |

**Gate for the whole track:** exposed wait ≤ ~15–20 % of wall (from 53 %). **Acceptance is bitwise
identity for A1–A3** — no tolerance machinery, which is why this track is both bigger and safer.

**One item moved here from the kernel plan:** `ComputeADERSharedFaceFluxRHS` — 10.56 % of compute
with a **2× spread across ranks** even in a GTS run. That is a communication/imbalance signature,
not a kernel one; it belongs to Track A's imbalance analysis.

## Track B — kernels (SECONDARY, 776 s/sim-s), reordered cheapest-first

The mechanism is settled by measurement: **hand-tiled fusion (bench variant B, 2.3–6.3× contended)
and/or per-element dgemm (variant D, 2.8–8.3×)**, keeping each element's data in cache across all
time-expansion levels. `DenseTensor` remains as contiguous storage; `BatchedLinAlg` is out.

> **B1 and B2 savings do NOT add — B1 is a down payment on B2.** Both target the same predictor
> stage; the 776 s/sim-s figure counts the predictor **once**, at the bench's ~4× end state. B1
> exists because it reaches part of that win for a fraction of the effort and risk, not because it
> is additive. If B1 alone clears the B2 gate, B2 is re-scoped, not stacked.

### B1 — fuse the accumulation sweeps (NEW, cheapest, do first)
Perf shows **`Vector::Add` 18.9 % + `Vector::operator=` 9.3 % = 28 % of self-time in pure
whole-vector traffic** — the `O+1` full-length AXPYs and zeroing per expansion level. Fusing those
accumulations into the existing per-element loop is a **much smaller change than full tiling** and
attacks a measured 28 %. *Gate:* ≥1.2× on the predictor stage vs B0; ≤1e-12 parity (round-off
lever convention). Stop-and-reassess if <1.1×.

### B2 — tiled-fused predictor (the 56 %)
Cache-blocked tiles (E≈32–128) inside LTS clusters; operators read once per element per step,
all levels fused; D(k) still scattered to `dk_retain` bit-identically (seam contract). *Gate:*
≥1.6× vs B0 on the predictor stage (bench says 2.3–6.3× is available on the kernel itself).

### B3 — volume stage (7.9 %)
Same tiling applied to the volume update. *Gate:* ≥1.5× on that stage.

### B4 — face-interior tables (13 %, DEMOTED from v1's Phase 1)
The face cache already removed the LU/basis tax, so the remaining headroom is the fold/contraction
only. **Re-justify before building:** if B4's projected saving is <5 % of the step, drop it. The
reference-catalog table design and its corrected memory budget survive from v1 unchanged.

### Closed / withdrawn
Fault & friction (0.65 %) · BLAS relink (no win) · `BatchedLinAlg` mechanism (slower) ·
the ≤40 µs·core target · SAFS bimaterial + mixed-flux extension (still a named later phase,
unchanged, and now clearly gated behind both tracks proving out on TPV104).

## Sequencing

**Prerequisite (DONE):** the two LTS Caliper scopes, so both tracks can measure their own gates on
the LTS leg. **The first LTS run of either track must publish the re-measured LTS stage split** —
that retires the GTS-proxy caveat above and is the last unmeasured input to this plan.

**A1 first, and B1 in parallel.** They touch different code (A1: LTS steppers + seam corrector;
B1: the predictor's accumulation loops), both are small, and both have the best
payoff-per-risk in their track. Then A2/A3 and B2 in sequence. Everything stays opt-in behind
flags; defaults and the TPV/BP5 byte-exact contract are untouched throughout.

**Re-baseline rule:** B0 = 187.5 µs·core is the compute denominator; the LTS-leg budget above is
the end-to-end denominator. Because the two tracks now run concurrently, **each track's gate is
measured with the other track's flags OFF**, and the composed number is measured once at the end.

## Honest ceiling

Even with both tracks fully successful we are **~3–4× behind SeisSol at matched accuracy**, because
their element kernels are hand-tuned assembly we are not attempting to match head-on. The route to
genuine parity is the **far-field order drop** (cheap low-order math away from the fault, expensive
high-order near it) — structurally impossible for SeisSol, whose order is fixed at compile time —
but it trades accuracy where we argue it doesn't matter, which is a different and more careful
claim than "faster at identical settings." That remains a separate proposal, not part of this
program.
