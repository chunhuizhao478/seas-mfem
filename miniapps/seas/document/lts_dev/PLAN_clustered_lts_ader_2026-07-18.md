# Implementation Plan: Clustered Local Time Stepping (LTS) for the ADER path of seas_spatial_dyn_driver

**Date:** 2026-07-18 (rev 5) · **Status:** IN PROGRESS — Phase 0 DONE, Phase 1 PART-DONE + resequenced, Phase 2 IN PROGRESS (foundation + predictor byte-validated; corrector remaining)
**Rev 5:** Phase 0 + part of Phase 1 implemented on `safs-v4_0_0-alt-case1-mfem-speed`
(local, unpushed); three Phase-1 pieces (fault-QP reorder, serial-mesh clustering,
LTS-aware partition) **resequenced** to Phases 3/4 where each is first needed and
first validatable. Phase boundaries updated; Appendix A interfaces unchanged;
Appendix B phase tags corrected to match the move (B.8 split into B.8/B.8b; the
reorder canary re-tagged Phase 3). Impl review: `REVIEW_phase01_impl_2026-07-18.md`.
**Rev 4:** added **Part I — The method, explained**: a plain-language tutorial
(with the equations, the measured cluster histogram, seam/tick diagrams, and
the beat-SeisSol arithmetic) so the plan is self-contained for a reader who has
seen none of the analysis documents.
**Grounding:** two multi-agent investigations over this repo and the local SeisSol
source, adversarially verified; quantitative motivation in
`../code_optimization_dev/ANALYSIS_mfem_vs_seissol_speed_2026-07-18.md`; method
selection (all LTS/multirate families, 3-judge panel) in
`ANALYSIS_lts_method_selection_2026-07-18.md`.
**Rev 2:** adopted the EDGE-2022 upgrades (λ-wiggle, Nc-cap+auto-merge,
flux-premultiplied exchange) and the staged far-field p-drop.
**Rev 3:** full adversarial plan review (two independent reviewers; findings
P-001…P-024 in `../../REVIEW.md`) applied: the tick scheduler is corrected
(P-001), ONE accumulate-buffer design is chosen and specified (P-002), sync/
truncation semantics unified (P-003), current-step-dt threading (P-004), the
normative t_origin formula (P-005), the fault-QP reorder single-source rule
with its full consumer list (P-006), the per-tick collective count formula and
not-due-rank rule (P-007), λ re-binning made real + FP-safe integer binning
(P-008), all rev-2 amendments propagated into the phase contracts (P-009), the
in-place single-Q contract (P-010), per-kind absolute nucleation with range
overloads (P-011), a physically correct conservation harness (P-012), and two
new normative appendices: **A. Interfaces** and **B. Unit-test matrix**.

## Summary — read this first

**The problem.** Our dynamic-rupture driver advances every element of the mesh
with one shared time step, set by the single worst element. On the SAFS regional
meshes the allowed step sizes span a factor of about one thousand, so ~97% of
all element updates are wasted work. Measured on the coarse SAF-ALT
benchmark: SeisSol finishes 150 simulated seconds in under five hours while our
driver would need days — and our code is actually *faster per element update*
than SeisSol. The entire gap is scheduling.

**The fix.** Group elements into clusters by their allowed time step (each
cluster's step is a power of two times a base step), and advance each cluster at
its own rate. Neighboring clusters stay consistent because the ADER predictor
already produces, for every element, a small polynomial describing its solution
over the whole step — a coarse element hands a fine neighbor exactly the
time-integrated information the fine steps need. This is SeisSol's mechanism,
upgraded with the published EDGE improvements SeisSol ships only as
experimental, plus one axis SeisSol structurally cannot copy (per-element
polynomial order, staged later).

**Expected outcome.** ~35× realized element-update reduction on the coarse
SAF-ALT mesh (38.73× ideal × the 94–95% realization EDGE demonstrates), ~45×
with the staged far-field order drop — from ~2.4 days to hours for the 150 s
benchmark, projected 1.2–1.5× faster than SeisSol itself. The acceptance gate
claims only ≥15× end-to-end.

**Main tradeoff / biggest risk.** LTS produces a different (equally valid)
trajectory than global stepping — bit-for-bit comparisons become tolerance
comparisons. The riskiest machinery is cluster-boundary bookkeeping: get one
sub-interval, buffer fill, or truncated step wrong and the result is silent
non-conservation, not a crash. Rev 3 therefore specifies the scheduler, the
buffer lifecycle, and their invariants normatively, with dedicated unit tests
for exactly those failure modes (Appendix B).

**What this does NOT do.** No single-element-update speedup (kernel efficiency
is a separate axis); ADER-only (RK/mixed-flux path excluded by guard); zero
change to any default-GTS result — LTS lands strictly opt-in, and "default"
means default for this driver's SAFS production configs after validation, with
TPV gold decks explicitly pinned to global stepping.

## How to read this plan
- The **Summary** above is the whole idea in five paragraphs.
- **Part I — The method, explained** (next section) is the tutorial: what LTS
  is, how it works in THIS code, and why we expect to beat SeisSol — plain
  language, with the equations and diagrams. Read it once before the phases.
- **Normative sections** (Scheduling, Buffer design, Appendices A/B) are the
  implementation agent's contract — nothing there is optional.
- Skim the per-phase **"In one sentence"** lines for the arc.
- The **Glossary** defines every shorthand.

---

# Implementation status & phase resequencing (2026-07-18)

**Read this before the phases.** Phase 0 and part of Phase 1 are implemented; three
Phase-1 pieces were **resequenced** to the phase where each is first *needed* and
first *validatable*. Nothing was dropped — the work moved to its natural home, and
the moved requirements now live under Phases 3 and 4 below. Branch:
`safs-v4_0_0-alt-case1-mfem-speed` (local, unpushed). Impl review record:
`REVIEW_phase01_impl_2026-07-18.md`.

## Shipped and locally validated
- **Phase 0 (complete):** `dynamic/lts_clustering.{hpp,cpp}` (A.1) +
  `test_lts_clustering` (4043/4043); `--lts-report` hook + per-element CFL
  accessors; `lts="off"` byte-exact confirmed (TPV102 ADER smoke + wave-op suites).
- **Phase 1 (part):** `dynamic/lts_layout.{hpp,cpp}` + tick table
  `dynamic/lts_stepper.hpp` (A.2/A.3) + `test_lts_layout` (87/87); `[numerics].lts*`
  config + validators + guards; run-path clustering + layout wiring, gated on
  `lts != "off"`, built at np==1, still stepping GTS.

## The three deferred pieces are ONE coupled unit
The **fault-QP reorder**, the **serial-mesh (rank-0) clustering**, and the
**LTS-aware METIS partition** were NOT landed in Phase 1 because they share a
dependency chain:
- The reorder needs the cluster ids **at wave-operator construction** (before the
  fault-DOF layout is fixed).
- Deterministic, **rank-count-independent** ids need clustering on the SERIAL mesh
  on rank 0 (a per-rank fixpoint disagrees at seams).
- The METIS partition must choose `part_data` **before the ParMesh** — which needs
  the material sampled before the ParMesh, but today the material is built *after*
  it (`spatial_dyn_driver.cpp:~1502` vs the `ParMesh` ctor at `:1310`).

Two of the three (serial clustering, partition) are meaningful **only at np>1**.
The plan already runs np=1 through Phases 2–3 and first goes MPI in Phase 4, so
the np>1 machinery belongs in Phase 4 — not Phase 1.

## The corrected sequence
| Deferred piece | Old home | New home | Why there |
| --- | --- | --- | --- |
| Fault-QP reorder | Phase 1 | **Phase 3, step 0** | First consumer is the per-cluster friction range-sweep; its byte gate (TPV104-spatial stations identical under the permutation) is an np=1 permutation-invariance test — largely local. Needs cluster ids at operator-construction time, obtained by clustering in the np=1 window *between* material construction and operator construction (no material-before-ParMesh refactor required at np=1). |
| Serial-mesh clustering | Phase 1 | **Phase 4** | Sole purpose is rank-count-independent ids; first exercised by the np=2==np=1 gate. |
| LTS-aware METIS partition (+ material-before-ParMesh) | Phase 1 | **Phase 4** | `part_data` is an np>1 concern; the "single-cluster byte gate *given the LTS partition*" lives in Phase 4. |

## Net effect on the arc
- **Phase 2 is unblocked** — it consumes exactly the layout already shipped and
  can start now at np=1.
- **Phase 3 gains a step 0** — the reorder — plus a re-inserted **early canary**:
  the old Phase-1 gate ("lts=rate2 + still-GTS stepping, TPV104-spatial stations
  byte-identical") runs BEFORE the friction range-sweep, so a reorder bug and a
  stepper bug cannot mask each other.
- **Phase 4 gains the pre-ParMesh clustering pipeline** — serial clustering +
  partition + the material-before-ParMesh refactor + cross-rank maxdiff.
- **Estimates do not inflate** — the reorder's time moves Phase 1 → Phase 3; the
  partition/serial-clustering time moves Phase 1 → Phase 4. Total ≈ unchanged.
- **Determinism** — until serial clustering lands (Phase 4), cluster ids are
  per-rank and correct only at np=1, which is exactly the regime of Phases 2–3.

---

# Part I — The method, explained

*(Plain language. Every number below is measured or verified — sources:
`../code_optimization_dev/ANALYSIS_mfem_vs_seissol_speed_2026-07-18.md` and
`ANALYSIS_lts_method_selection_2026-07-18.md`.)*

## I.1 Why global time stepping wastes ~97% of our work

Every explicit wave solver must respect a per-element stability limit (the CFL
condition). For our ADER-DG tetrahedra, the largest stable time step of element
*e* is

```
             1        2·r_e                       6·V_e
dt_e = cfl · ────── · ──────    with   2·r_e  =  ───────   (insphere diameter)
             2N+1     c_p,e                        A_e
```

where N is the polynomial order, V_e the element volume, A_e its total face
area, and c_p,e = √((λ+2µ)/ρ) the P-wave speed at the element. Small element ⇒
small dt; slow rock ⇒ bigger dt. **This is verified to be exactly SeisSol's
formula**: at `cfl_dg_safety = 1` our driver prints dt_cfl(p3) = 4.12183e-5 s
on the benchmark mesh, and SeisSol's log for the same mesh prints "Minimum
timestep: 41.2181 µs" — agreement to five digits, same critical element, same
material pairing.

A regional fault model is *inherently* multiscale: ~100–200 m elements resolve
the fault geometry, kilometre-scale elements fill the far field, and a few
slivers from the CAD-forced fault surfaces go down to ~14 m. Here is the
benchmark mesh (1,319,294 tets), binned by allowed dt into powers of two
(measured from the SeisSol log; the dt values shown are for order 4 — our p1
values are 7/3 larger but the *ratios*, and hence the bins, are identical):

```
cluster   allowed dt      cells     share of mesh
   0        41 µs           123     ▏  0.01%   ← these 123 slivers set the global dt
   1        82 µs           822     ▏  0.06%
   2       165 µs         3,346     ▏  0.25%
   3       330 µs        13,917     ▍  1.1%
   4      0.66 ms       171,472     █████        13.0%
   5      1.32 ms       487,273     ██████████████  36.9%
   6      2.64 ms       192,506     █████▌       14.6%
   7      5.28 ms        87,246     ██▌           6.6%
   8     10.6  ms       312,888     █████████    23.7%
   9     21.1  ms        49,565     █▌            3.8%
  10     42.2  ms           136     ▏  0.01%
```

**Global time stepping (GTS, today)** advances *all* 1.32M elements at the
41 µs of those 123 slivers. Total element-updates for 150 simulated seconds:

```
GTS:   N_elem · T/dt_min  =  1,319,294 · (150 s / 41.2 µs)  ≈  4.80·10¹²
```

But each element only *needs* T/dt_e updates. If cluster c steps at
dt_min·2^c, the necessary work is

```
LTS:   Σ_c  cells_c · T/(dt_min·2^c)  =  (150/41.2 µs) · Σ_c cells_c/2^c
    =  3.64·10⁶ · 34,063  ≈  1.24·10¹¹        ⇒  4.80·10¹² / 1.24·10¹¹ = 38.7×
```

**38.7× of our element updates are pure waste** — elements updated 30–1000×
more often than their own physics requires. (SeisSol's log prints "speedup
111.85", but that is an arithmetic-mean statistic; the honest workload number
is this harmonic 38.7×.) The same sum over the fault-face column gives 22.0×
for fault work — smaller, because 97% of fault faces sit in the fine clusters
4–5.

And the punchline that makes this worth doing: we measured that our code
performs element updates **1.354× faster than SeisSol** (9.89 vs 7.31 million
updates/s on identical 256 cores). Our code was never slow — it was only ever
*scheduled* badly. LTS is a scheduler fix.

## I.2 The one lucky fact: our predictor already computes what LTS needs

The only hard problem in LTS is the *seams*: where a fast cluster touches a
slow one, the two sides live at different moments in time, yet every face flux
must be computed from both sides consistently — and conservatively.

Here our ADER discretization pays off. Each macro step begins with the
Cauchy–Kovalevskaya (CK) predictor: for the elastic system
∂_t Q + A_x∂_x Q + A_y∂_y Q + A_z∂_z Q = 0 (Q = the 9 stress+velocity fields),
each element locally converts spatial derivatives into time derivatives,

```
D(0) = Q(t_n),        D(k+1) = − Σ_d  A_d · ∂_d D(k)      (element-local recursion)
```

and thereby owns a small polynomial **in time** — a short "movie" of its own
next step:

```
Q_e(τ) = Σ_k  (τ^k / k!) · D(k),        valid for τ ∈ [0, dt_e]
```

This polynomial is *the complete answer to "what is this element doing during
its step?"* — which is exactly the question a neighbor at a different time
level needs answered. Any time-slice of the flux data can be produced by
integrating it exactly:

```
 b                              b^{k+1} − a^{k+1}
∫ Q_e(τ) dτ   =   Σ_k  D(k) · ───────────────────        (exact, closed form)
 a                                  (k+1)!
```

Today we build Q_e(τ), use it once, and **throw the coefficients away**
(they live in reusable scratch buffers). The entire LTS coupling mechanism is:
*keep them* — for the ~thin shell of elements that border a finer cluster —
and let the fine side integrate its coarse neighbor's polynomial over exactly
the sub-interval it needs. No interpolation error, no new discretization: the
seam coupling is algebraically exact in time.

## I.3 One coarse step, in pictures

Two clusters; the coarse (c=1) steps once while the fine (c=0) steps twice:

```
 time
  ↑            COARSE element                    FINE element
2dt₀ ─┤  ← both sides arrive here together (sync of this pair)
      │   correct: adds its accumulate      correct sub-step 2:
      │   buffer (= the two deposited          flux uses ∫ Q_coarse(τ)dτ over [dt₀,2dt₀]
      │   flux integrals), applies M⁻¹,        → deposits coarse's share → buffer
 dt₀ ─┤   Q += …                            correct sub-step 1:
      │        ▲                               flux uses ∫ Q_coarse(τ)dτ over [0,dt₀]
      │        │  Q_coarse(τ) = Σ τᵏ/k!·D(k)   → deposits coarse's share → buffer
  0  ─┤   predict (builds D(k), RETAINED)   predict
```

Walkthrough:
1. **Both predict at τ=0.** The coarse element's D(k) stack is retained (it is
   a "provider"); the fine one's is used as usual.
2. **Fine sub-step 1** computes the seam-face flux using its own state and the
   coarse polynomial integrated over [0, dt₀] — *both sides of the flux are at
   consistent times.* It applies its own share of the flux to itself,
   and deposits the **coarse element's share** into that element's
   *accumulate buffer* (the coarse element is asleep; the buffer is its
   in-tray).
3. **Fine sub-step 2** does the same over [dt₀, 2dt₀].
4. **Coarse correct** wakes up once: instead of visiting the seam face, it
   empties its in-tray — which now holds exactly the face-flux integral over
   its whole step — applies its mass inverse, and advances.

**Conservation is exact by construction:** each seam flux is evaluated ONCE
per fine sub-interval, and both sides consume *that same* integral — whatever
momentum leaves one side enters the other, to machine precision. (Getting the
in-tray bookkeeping right — filled *before* it is emptied, emptied exactly
once, zero at every sync — is the part our adversarial review found subtly
wrong in rev 2 and is now specified normatively with its own unit tests.)

For many clusters the pattern nests. One **sync interval** = one step of the
coarsest cluster; inside it runs a fixed **tick table** (tick = one finest
step). Three clusters:

```
tick            0        1        2        3      ← 4 ticks = one sync interval
cluster 0:    P C      P C      P C      P C      (4 steps of dt₀)
cluster 1:    P         C       P         C       (2 steps of 2·dt₀)
cluster 2:    P                           C       (1 step of 4·dt₀)

P = predict (opens a step; due when  tick mod 2^c == 0)
C = correct (closes a step; due when (tick+1) mod 2^c == 0), finest first
```

At every sync point all clusters are at the same time: outputs, V_max
reductions, checkpoints happen there and only there. Every MPI rank walks the
*same* table — even ranks with nothing to do this tick — so the number of MPI
calls per tick is identical everywhere by construction. That property is what
protects us from the matched-collective deadlocks this codebase has been
burned by before.

**The fault** needs no special coupling trick at all: both elements of every
fault face are forced into the same cluster, so the rate-and-state friction
solve simply runs at that cluster's rate over that cluster's block of fault
points — the same solver, the same sub-step structure, just per cluster. Since
97% of fault faces live in clusters 4–5, friction runs ~16–32× less often than
today while remaining exactly as resolved *relative to its own local physics*.

## I.3b Common questions about the seam

**Q1 — "So the coarse side takes results computed by the fine side?"**
Half right — the direction matters, and it is different for the two things
that cross the seam. First recall: in DG, elements never touch each other's
solutions; they interact ONLY through face fluxes. So "taking the neighbor's
result" can only ever mean "taking the flux through our shared face."

- **Coarse → fine: a FORECAST flows, and the fine side does the computing.**
  While the fine cluster sub-steps, the coarse element is asleep — it has no
  fresh results to give. What it has is its predictor polynomial Q_coarse(τ),
  its forecast over the whole big step. The fine side integrates that forecast
  over each of its own sub-intervals and evaluates the flux itself.
- **Fine → coarse: COMPUTED FLUX RESULTS flow — yes, the coarse takes them.**
  The coarse element eventually needs the seam flux integrated over its whole
  step, and the fine side has already computed exactly those fluxes (once per
  sub-interval, with time-consistent data on both sides — it is the only party
  that ever had both). Each fine sub-step deposits the coarse element's share
  into the coarse element's in-tray (accumulate buffer); the coarse correction
  just adds the in-tray. It never visits the face itself.

Slogan: **predictions flow downhill (coarse→fine); computed fluxes flow uphill
(fine→coarse); states never cross the fence at all.** Evaluating each seam
flux exactly once and letting BOTH sides consume that same evaluation is also
what makes conservation exact — two independently computed versions would
disagree slightly and leak momentum at every seam.

**Q2 — "Why does the fine side need the coarse forecast at all?"**
Because a face flux is a two-sided quantity at a single instant:
F(τ) = F(Q_fine(τ), Q_coarse(τ)). The fine side knows its own Q_fine(τ) — it
is live. But the coarse element's STORED state is frozen at the start of its
big step; halfway through, it is stale by up to 2^Δ fine steps. The coarse
element is not actually sitting still during its step — its solution evolves
the whole time, and the forecast polynomial is precisely the description of
that evolution.

Picture a P-wave traveling from the coarse side, reaching the seam midway
through the coarse element's big step:

```
                    coarse (one big step)     │     fine (4 sub-steps)
  wave →→→→→→→→→→→→→→→→ ⟍                    │
                          ⟍  arrives at the  │
                            ⟍ face at ½·dt_c  │
```

With the forecast, fine sub-step 3 evaluates Q_coarse(½dt_c) — the polynomial
contains the wave's arrival — and the wave flows into the fine region on time,
at full order. With the frozen state, the wave would not exist at the seam
until the coarse element's NEXT step: delayed, staircase-distorted, and the
scheme drops to first order exactly at every cluster boundary. There is no
third option: the fine side cannot wait for the coarse to compute more often
(that IS global stepping), and it cannot extrapolate the coarse's interior
dynamics from its own data.

The reassuring part: **this is not a new approximation.** Even under today's
global stepping, no ADER face flux is ever computed from frozen states — the
method is predictor–corrector by construction, and every flux everywhere is
already built from both elements' predicted time evolutions. LTS adds exactly
one twist: the forecast is integrated over a sub-interval [a,b] instead of the
whole step — and since it is a polynomial, that slice is algebraically exact
(the (b^{k+1}−a^{k+1})/(k+1)! formula). Same ingredient, sliced thinner; no
new error term at the seam, which is why LTS-ADER preserves the scheme's full
convergence order.

## I.4 Why this can beat SeisSol, not just match it

Wall-clock time factorizes as

```
T_wall  =  (number of element updates)  ×  (cost per update)  ÷  (parallel efficiency)
```

SeisSol's entire advantage on this benchmark is the first factor. We already
win the second (1.354× faster per update, measured). The plan attacks the
first factor with the same clustering — and then goes further on four axes:

| Lever | Gain | Why SeisSol doesn't have it |
|---|---|---|
| λ-wiggle: rescale all bin edges by λ∈(0.5,1] so a fat population just above an edge drops into the next-coarser cluster | +17.5% shown on LOH.3 (EDGE) | exists only as an experimental flag |
| Cluster cap (~5–6) + cost-model auto-merge | keeps batches big — this is also the published fix for GPU LTS (SeisSol's GPU LTS collapsed to ~1.3× without it) | experimental flag |
| Flux-premultiplied 3-buffer seam exchange (send the flux-projected payload, not raw coefficient stacks) | 1.26–1.48× over SeisSol's LTS comm; 94–95% of the theoretical speedup realized at scale (EDGE, IPDPS 2022) | EDGE-only, never merged |
| **Far-field order drop (Phase 7):** clusters ≥6 — 642,341 cells, 49% of the mesh, only 2,509 fault faces (1.6%; fault-adjacent elements are EXCLUDED from the drop and keep p3) — run p1 in the far field | ~1.16× fewer update-costs + large memory-bandwidth relief | **structurally impossible**: SeisSol's polynomial order is fixed at compile time; MFEM controls its own element layout |

Putting the measured numbers together for the p1 benchmark
(sim-seconds per wall-hour, 2 nodes / 256 cores):

```
GTS, safety=3  (measured, jobs 522287xx)      0.87   ▏
GTS, safety=1  (deployed 2026-07-18)          2.6    ▍
SeisSol o4 + LTS (measured, job 52042189)    31.8    ███████
MFEM p1 + LTS  (projected ceiling)         ~100      ██████████████████████
   = 0.87 × 3 (safety) × 38.7 (LTS)  =  31.8 × 7/3 (order dt) × 1.354 (per-update)
MFEM p1 + LTS, EDGE-realized (94–95%)       ~95      █████████████████████
+ far-field p-drop (Phase 7)               ~110      ████████████████████████
```

The two identities in the middle line are the same number computed two
independent ways — that closure (within 8%) is what makes the projection
trustworthy rather than hopeful. The formal acceptance gate claims only ≥15×
over the safety-1 GTS baseline, leaving honest room for cluster-management
overhead and load imbalance.

## I.5 Where we deliberately differ from SeisSol (engineering choices)

| Aspect | SeisSol | This plan | Why |
|---|---|---|---|
| Scheduler | asynchronous actor model, mailbox messages between cluster actors | deterministic tick table, identical on every rank | actors exist to overlap MPI; we buy correctness and reproducibility first — overlap is a v2 optimization if profiles demand it |
| Cross-rank fault | fault faces may straddle ranks (dedicated copy-layer DR machinery) | fault faces forced rank-interior (fault-locality partition implied by LTS) | deletes an entire hazard class (our historical R-1601/R-101 shared-fault defects) instead of managing it |
| Seam buffer semantics | fine element's own state-integral buffer; coarse still evaluates the face flux | flux-contribution in-tray on the coarse element; fine side does the work for both | matches the premultiplied-payload exchange; one flux evaluation per sub-interval = conservation by construction |
| Predictor variants | separate GTS/LTS kernels | one fused CK predictor (our existing `--shared-ck-recursion` path) with retained coefficients for providers | least new numerics; byte-exact GTS fallback preserved |
| Partition weights | scalar cost ~ 2^-cluster | multi-constraint (one balance constraint per cluster level), scalar fallback | published stronger scheme (Rietmann); extends our existing partition plan |
| Order | global, compile-time | per-element (two order classes, Phase 7) | the axis SeisSol cannot follow |

---


## Glossary

| Label | Plain-language meaning |
| ----- | ---------------------- |
| GTS | Global time stepping — every element uses the same dt (today's behavior). |
| LTS | Local time stepping — each cluster of elements uses its own dt. |
| dt_base | λ · min over elements of dt_e (the wiggle-scaled base step). Cluster c steps at dt_c = dt_base·2^c. |
| rate-2 binning | Element joins cluster c(λ) = max{c : dt_base·2^c ≤ dt_e}, computed by an integer comparison loop (never floor/log2 — FP determinism). Lower-edge binning ⇒ dt_c ≤ dt_e always (CFL-safe by construction). |
| maxdiff rule | Face-neighboring elements' clusters differ by ≤1; across fault faces by exactly 0. Iterative clamp to a fixed point on the serial mesh. |
| D(k) stack | Per-element Cauchy–Kovalevskaya Taylor coefficients; Q(τ)=Σ τ^k/k!·D(k), expansion point = the element's last PREDICT time. Retained (raw, unscaled) for provider elements. |
| provider element | An element with ≥1 face whose neighbor is one cluster finer. It retains its D(k) stack for consumers to integrate. |
| consumer face | A face seen from the FINE side whose neighbor is one cluster coarser. The fine side does the flux work for both sides there. |
| accumulate buffer | Per COARSE (consumer-owning) element: NUM_STATE×ndof_per_el of **pre-M⁻¹ residual** units. At each fine sub-correct, the fine side assembles the coarse element's face-flux contribution for that sub-interval — evaluated with the coarse side's own Godunov flux matrices and tested with the coarse element's basis (NOT the negation of the fine-side term; A± differ across bimaterial faces) — and adds it here. The coarse correct adds-then-zeros this buffer in place of visiting those faces, before its M⁻¹. |
| tick | One dt_base step inside a sync interval. tick times come from closed forms (`t(tick) = t_s + tick·dt_base`), never accumulation. |
| epoch counter | Per provider element: its predict count. Per accumulate buffer: its fill count. Asserted at every consumption. |
| sync point | A global time all clusters land on exactly. Outputs, checkpoints, and reductions happen only here. With `lts_sync_dt="auto"` no cluster ever truncates inside a sync interval — truncation exists only on the final interval before tfinal. |
| the consume path | The verified shared-fault substep mechanism replacing the R-1601 inline fallback (not needed under D-2, retained as background). |
| fault-locality partition | Both elements of every fault face on one rank; IMPLIED automatically by lts≠off (P-017). |
| GAP registry | The grounded-investigation gap IDs cited in this plan — defined in the table below. |

### GAP registry (P-020)

| GAP | Meaning |
|---|---|
| GAP-A1 | Stability invariant: every element steps at dt_cluster ≤ its own dt_e (lower-edge binning); asserted in production, not just tooling. |
| GAP-A3 | D(k) lifetime: a provider's stack must survive, unmodified, across all fine sub-steps of its step; epoch-asserted. |
| GAP-B1 | The accumulate-buffer graft onto the one-visit-per-face flux loop (resolved by the Buffer design section). |
| GAP-B2 | Sub-interval bookkeeping at rank seams (resolved by the t_origin formula + Phase 4). |
| GAP-B3 | Full audit of every dt read inside the corrector (resolved by current-step-dt threading, P-004). |
| GAP-C2 | The five rank-global fault structures: (1) `dof_data`, (2) `fault_coords`, (3) `Q_pointwise_plus/minus` traces, (4) `I_imp_plus/minus_flat`, (5) the iterator's deltaT/weights/tau_nodes schedule. All split per cluster in Phase 3. |
| GAP-C4 | Step-keyed gates (R-101, DIAG) rekeyed to sync points with time-based windows. |
| GAP-D1/D2 | Checkpoint partition/layout provenance + V1-under-LTS refusal (resolved in Phase 2 checkpoint spec). |
| GAP-X2 | No empirical speedup input existed → Phase 0 histogram tool is the go/no-go. |

## Technical Overview

The ADER macro step (`drivers/spatial_dyn_driver.cpp:399-515`) is predictor →
fault-friction substeps → corrector. The predictor
(`wave_operator.inl:1597-1688`) is an element-local CK recursion producing
exactly the Taylor object SeisSol's LTS couples with; it currently discards the
coefficients (ping-pong scratch `:1644-1686`) and the corrector consumes one
rank-global time-integral `I` (`:5832-5859`). LTS restructures this into
per-cluster invocations driven by a deterministic tick schedule, retains D(k)
stacks for provider elements, adds accumulate buffers on consumer-owning
elements, integrates coarse-neighbor Taylor expansions over fine sub-intervals,
and splits the fault iterator (`friction_substep_iterator.*`) into per-cluster
sweeps over cluster-contiguous fault-face blocks. Clustering runs on the serial
mesh before partitioning; partitioning gains LTS-aware weights through the
existing explicit-partition injection point (`spatial_dyn_driver.cpp:1268-1293`).

## Constraints

- **Byte-exact default:** every touched shared file recompiles into the TPV
  gold-deck binaries (Makefile `:2297-2343`); with LTS off, behavior must be
  byte-identical (repo pattern: `use_shared_ck`, `resample`).
- **ADER-only:** `lts ≠ off` rejected unless `time_integrator=ader` AND
  `mixed_flux=none` (central flux is non-dissipative under ADER — the
  dissipative upwind flux is a stability precondition for the cluster
  interfaces).
- **LTS implies fused predictor (P-016):** `lts≠off` uses shared-CK semantics
  regardless of `--shared-ck-recursion` (rank-0 log line states it).
  `--deriv-cache` / `--face-cache` compose unchanged (time-independent caches).
- **LTS implies fault locality (P-017):** `lts≠off` enables
  `BuildFaultLocalityPartitioning` automatically (log:
  `[partition] fault-locality implied by lts`); the CLI flag becomes a no-op
  alias; vacuous on fault-free meshes. Every fault face is rank-interior ⇒ the
  shared-fault×cluster hazard class (R-1601/R-101) is structurally absent.
- **Deterministic clustering:** cluster ids are a pure function of the serial
  mesh + material + config. Binning uses the integer comparison loop
  (Appendix A.1); `floor(log2(...))` is forbidden (P-008).
- **Matched collectives (P-007):** every collective is gated exclusively on the
  global tick table and serial-mesh cluster metadata (global per-cluster
  element and fault-face counts, broadcast in Phase 1). Rank-local counts never
  gate a collective. The per-tick collective count
  `n_x(tick) = Σ_{c∈due(tick)} [9 + (global fault faces in c > 0 ? O : 0)]`
  is precomputed into the tick table and asserted per tick by a debug counter
  on every rank. Under D-2, global shared-fault QPs are zero ⇒ the
  predictor-substep exchange is dropped GLOBALLY (config-level, uniform).
- **cluster-0 dt definition:** binning consumes the FULL driver product
  `dt_e = cfl · CflSafetyFactor(cfg) · h_e / cp_e` (the order/safety factor is
  driver-side — `spatial_friction.hpp:731-736` — not in ComputeMaxDt).
- **No LTS edits leak** outside `dynamic/` + the spatial driver + new files
  (memory rule [C2]).

## Design decisions (defaults chosen; sign-off requested)

| ID | Decision | Chosen default | Rejected because |
|---|---|---|---|
| D-1 | Fault faces under LTS | Per-cluster fault machinery; DR faces force both sides same-cluster in v1 (Uphoff SC'17 rule). Relaxation to maxdiff≤1 across the fault = **Phase-5 Req 5** (gated experiment, non-blocking, `lts_fault_maxdiff=1`, default 0). | All-fault-at-min-cluster costs 4.7× the whole LTS budget. |
| D-2 | Shared fault faces | `lts≠off` IMPLIES fault-locality partitioning (P-017). | Per-cluster cross-rank fault coupling doubles v1 risk. |
| D-3 | Nucleation under LTS (P-011) | ALL kinds route through their idempotent ABSOLUTE forms via `INucleationMethod`: `gradual_overstress`→`ApplyGradualOverstressAbsolute`; `gradual_overstress_compact_circular`→`ApplyGradualOverstressCompactCircularAbsolute`; `instantaneous_overstress_circular` unchanged (one-shot pre-first-sync). Both gradual appliers gain range overloads `(…, qp_begin, qp_end)`; each cluster applies its own range at its own stage times, pinned to the GTS sub-step time convention. GTS path untouched. | Per-cluster incremental telescoping duplicates state; whole-vector absolute writes stomp other clusters' τ_nuc at wrong times. |
| D-4 | Scheduler | Deterministic tick schedule (Normative Scheduling below); SeisSol's actor model deferred to v2 (MPI overlap only). | Correctness first; the tick loop preserves matched collectives trivially. |
| D-5 | Default flip | Two-stage: opt-in first; after Phase-5 acceptance, flip the parser default in the same commit that pins `lts="off"` in every TPV-spatial config + re-goldens the TPV104-spatial smoke. | One-shot flip silently changes TPV104-spatial production trajectories. |
| D-6 | λ-wiggle / auto-merge | **IN v1 from Phase 0**, inside `BuildLtsClustering` (full interface: Appendix A.1). λ-scan re-bins elements against λ-scaled edges (P-008); Nc-cap≈6 + cost-model auto-merge (merges only remove the top level, one at a time, re-evaluating cost; merge only ever reassigns an element to a SMALLER c — CFL-safe by direction). | They are the GPU story (SeisSol GPU-LTS collapsed to ~1.3× without the cap) and the realization-fraction story (EDGE 94–95%). |
| D-7 | Exchange payload | Phase 4 lands as **4a** (full ghost-field exchange per due tick; parity gate) then **4b** (EDGE fixed 3-buffer flux-premultiplied exchange; byte-compared against 4a; **Phase-5 performance is measured on 4b**). | Terminology unified (P-009); measuring Phase 5 on 4a would gate the project on machinery the plan says is not the target. |

## NORMATIVE: Scheduling (P-001, P-003, P-004, P-005)

**Sync grid.** `lts_sync_dt="auto"` ⇒ the sync time advances by
`T_s = dt_base·2^(Nc−1)` (the coarsest cluster dt). By construction no cluster
truncates inside a sync interval; truncation exists ONLY on the final interval
before `tfinal` (which may be shorter than T_s). A numeric `lts_sync_dt` is
snapped DOWN to the nearest positive multiple of the coarsest dt (rank-0 log
`sync_dt snapped X→Y`). Outputs/checkpoints/reductions are evaluated AT sync
points; the V_max-adaptive output cadence consumes the per-sync REDUCED
(global) V_max, so cadence decisions are identical on all ranks. Output
cadences are re-expressed as `ceil(dt_out/T_s)` sync intervals.

**Tick loop.** Within a sync interval of length `T_actual ≤ T_s`
(`< T_s` only on the final interval):

```
ticks_per_sync = ceil(T_actual / dt_base)
for tick in 0 .. ticks_per_sync-1:
    for c in predict_due(tick):            # any fixed order; element-local
        predict(c, dt_step(c, tick))       # CK; retain D(k) for providers; write I_c
    for c in correct_due(tick), FINE→COARSE:
        correct(c, dt_step(c, tick))       # faces + volume + buffer + M^-1 + add

predict_due(c, tick)  ⇔  tick % 2^c == 0
correct_due(c, tick)  ⇔  (tick+1) % 2^c == 0  OR  tick+1 == ticks_per_sync
dt_step(c, tick) = min(dt_base·2^c, T_actual − dt_base·(tick − tick % 2^c))
```

predict opens the step starting at the tick; correct closes the step ending at
tick+1. FINE→COARSE ordering within the correct phase guarantees the fine
cluster's final sub-interval contribution lands in the accumulate buffer
before the coarse consumes it in the same tick.

**Closed-form times (never accumulate):**
`t(tick) = t_s + tick·dt_base`;
`t_origin(c, tick) = t_s + dt_base·2^c·floor(tick/2^c)` (the cluster's last
predict time = its D(k) expansion point). A fine sub-interval on a consumer
face is `[a, b]` with `a = t(tick) − t_origin(coarse)`,
`b = min(t(tick)+dt_fine, t_s+T_actual) − t_origin(coarse)`. Because fine and
coarse truncate at the same sync time, `b ≤` the coarse's current step ≤ dt_c
— provable, and asserted (`|tracked − formula| ≤ 1e-12·dt_c`) at every
consumption.

**Residual-step merge:** a scheduled step shorter than `1e-10·dt_c` (FP
residue at tfinal) merges into the preceding step; the decision comes from the
closed-form schedule, identical on all ranks.

**Current-step dt threading (P-004 / GAP-B3):** every dt read inside the
corrector — `I/dt` (`wave_operator.inl:4319-4325, 4511`), `Q_imp·dt` (`:4741`),
PML inline (`:5907`), and the friction side-channel scale — takes the
`dt_step(c, tick)` value threaded as ONE argument from the tick loop. The
nominal dt_c may appear only as the untruncated value of dt_step. No site
recomputes dt.

**Invariants (asserted):** (i) at every `correct(c)`, each consumed accumulate
buffer's fill count equals the number of fine sub-steps of the closing step
(2^Δ untruncated; the closed-form count when ragged); (ii) all accumulate
buffers are exactly zero at every sync point; (iii) every provider read
matches the provider's current epoch (GAP-A3); (iv) `dt_cluster(e) ≤ dt_e`
for every element (GAP-A1, checked inside `BuildLtsClustering`, production
path).

## NORMATIVE: Buffer design (P-002; replaces the SeisSol "fifth rule")

One design: **flux-contribution accumulate buffers** (matches D-7's
premultiplied payloads).

- Storage: per consumer-owning (coarse) element, `NUM_STATE × ndof_per_el`,
  pre-M⁻¹ residual units, allocated only for elements in `consumer_owner_elems`
  (Phase 1 layout). NOT the global rhs vector — per-cluster corrector calls
  zero/overwrite `rhs`, so buffered contributions live in dedicated storage.
- Fill: at each fine sub-correct, the fine side evaluates the face flux once
  per QP for the sub-interval (single-flux-evaluation invariant), scatters its
  own side into its own rhs, and assembles the coarse element's contribution —
  computed with the coarse side's own Godunov flux matrices and tested with
  the coarse element's basis functions — into that element's buffer.
  "Anti-symmetric negation" is WRONG on bimaterial faces and is forbidden.
- Consume: the coarse `correct` adds the buffer into its rhs in place of
  visiting its consumer faces, then zeros it, before applying its per-element
  M⁻¹ (`ApplyMassInverse` is element-block-diagonal —
  `wave_operator.inl:5995-6021` — so per-cluster application is exact).
- Role-driven face sweep: an element may simultaneously carry (a) provider
  faces (neighbor one finer — skip: the fine side does the work), (b) GTS faces
  (equal cluster — visit normally, both-sided scatter as today), (c) consumer
  faces seen from the fine side (do the two-sided work + buffer fill). The
  sweep dispatches per face role from the Phase-1 layout tables. Per-element
  (not per-face) buffers are well-defined because maxdiff≤1 ⇒ ALL finer
  neighbors of an element are exactly c−1 (one accumulation cadence).
- The SeisSol "fifth rule" (accumulate-buffer cell must also provide
  derivatives) applies to SeisSol's state-buffer design and has NO REFERENT
  here; it is deleted. Provider-role assignment derives purely from face
  cluster differences.

## Phase 0: Cluster report & go/no-go  — ✅ DONE (2026-07-18)

**In one sentence:** Before writing stepping code, measure the cluster
histogram and predicted speedup for each target mesh with our own dt formula,
and verify against SeisSol's log on the shared mesh.

### Files to Create
- `dynamic/lts_clustering.hpp/.cpp` — pure functions (no MPI): Appendix A.1.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp` — `--lts-report` (implies `--dry-run`):
  compute per-element `dt_e = cfl·CflSafetyFactor(cfg)·h_e/cp_e` (matrix path:
  `GetPerElementCflLength()/GetPerElementMaterial()`; scalar path: promote the
  ctor loop temporaries `wave_operator.inl:65-104` to stored vectors —
  byte-exact-neutral), run clustering, print BOTH histograms (P-009):
  (a) RAW mode (`wiggle_scan=false, nc_cap≤0` ⇒ λ=1, no merge — SeisSol's
  algorithm) for the SeisSol cross-check; (b) PRODUCTION mode with the λ-scan
  curve and chosen (λ, Nc). Print element-update speedups in both SeisSol's
  arithmetic-mean form and the harmonic form.
- `dynamic/bimaterial_wave_operator.hpp` — re-label the two accessors from
  "Test-only" to production-sanctioned.

### Acceptance Criteria
- [ ] RAW histogram matches SeisSol's (123, 822, 3346, 13917, 171472, 487273, 192506, 87246, 312888, 49565, 136) on the shared mesh (per-cluster deviation ≤1% or explained by material sampling).
- [ ] `dt_cluster(e) ≤ dt_e` passes on all target meshes (both modes).
- [ ] PRODUCTION-mode predicted harmonic speedup ≥10× on the coarse ALT and ≥1 production mesh (**go/no-go**).
- [ ] `make test` unchanged; `test_lts_clustering` (Appendix B.1) passes.

**Estimate:** 3–5 days. Depends on: nothing. Required by: all later phases.

## Phase 1: Clustering wired into the run path (still GTS)  — ◑ PART-DONE + RESEQUENCED (2026-07-18)

**In one sentence:** The driver computes and carries the cluster layout on every
run (at np=1), while still stepping globally.

> **Resequenced.** The layout + tick table + config + run-path wiring shipped.
> The three np>1 / operator-construction-coupled pieces — the **serial-mesh
> clustering**, the **LTS-aware METIS partition**, and the **fault-QP reorder** —
> **moved out** (see "Implementation status & phase resequencing" above): the
> reorder to **Phase 3, step 0**; serial clustering + partition to **Phase 4**.
> Their normative requirements are struck through below and restated in their new
> phases.

### Files to Create
- `dynamic/lts_layout.hpp/.cpp` — run-side layout (Appendix A.2): per-rank
  per-cluster element index lists; per-cluster face lists tagged with roles
  (provider/GTS/consumer per the Buffer design); the derived ELEMENT sets
  `provider_elems` (coarse side of ≥1 maxdiff-1 face) and
  `consumer_owner_elems` with dense slot maps (`provider_slot_of_elem`,
  `buffer_slot_of_elem`) that Phase 2 sizes stores by (P-015); fault-face →
  cluster map; global per-cluster element/fault-face counts (broadcast — the
  collective-gating metadata); the tick-table generator (Appendix A.3).

### Files to Modify
- `spatial/code/spatial_friction.{hpp,cpp}` — parse `[numerics].lts`
  ("off"|"rate2", default "off"), `lts_nc_cap` (default 6), `lts_wiggle`
  ("scan"|numeric λ|"off", default "scan"), `lts_merge_loss_tol` (default
  0.05), `lts_sync_dt` (default "auto"), `lts_fault_maxdiff` (default 0).
  Guards: reject lts≠off with rk/mixed-flux; lts≠off implies fault-locality
  and fused-CK (Constraints).  **(DONE.)**
- `drivers/spatial_dyn_driver.cpp` — run-path clustering + `BuildLtsLayout`
  wiring, gated on `lts != "off"`, computed at **np==1**, still GTS. **(DONE.)**
- ~~`drivers/spatial_dyn_driver.cpp` — serial-mesh clustering before ParMesh;
  partition spec (P-009): rank 0 builds the element graph … `METIS_PartGraphKway`
  with `ncon = num_clusters` … companion `.cpp` (`metis.h ::real_t` collision).~~
  → **MOVED to Phase 4** (np>1 concern; needs material-before-ParMesh). Full
  spec restated there.
- ~~**Fault-QP reorder — single source of truth (P-006):** cluster-contiguous
  sort of the wave operator's fault-face list at construction, BEFORE
  `SetFaultDOFData` / station `Open` / nucleation / impedance-resolver seeding /
  first ParaView write; `fault_coords` follows; checkpoints serialize `dof_data`
  in canonical order via the permutation table; assert `n_shared_fault_qps==0`.~~
  → **MOVED to Phase 3, step 0** (first consumer is the friction range-sweep;
  gate is np=1 permutation-invariance). Full spec restated there.

### Acceptance Criteria (resequenced — what Phase 1 now delivers)
- [x] `lts="off"` byte-identical (`test-ader-tpv102-smoke` 4/4; wave-op suites 25/25 + 52/52).
- [x] `test_lts_clustering` incl. λ/merge cases (B.1, 4043/4043); `test_lts_layout` (B.2, 87/87).
- [x] `lts="rate2"` at **np==1** builds + logs the cluster histogram and layout; still steps GTS (output byte-identical to `lts="off"`).
- [ ] **→ Phase 3:** `lts="rate2"` TPV104-spatial stations byte-identical to `lts="off"` under the known permutation (requires the reorder — moved).
- [ ] **→ Phase 4:** np∈{2,10} identical cluster histograms across ranks (requires serial clustering — moved).

**Estimate:** ~1.5 weeks — **~1 wk delivered; the reorder (~0.5 wk) is now
Phase-3 step 0 and the partition/serial-clustering (~few days) is now Phase 4.**

## Phase 2: Multi-cluster stepping, bulk only (np=1)  — ◑ IN PROGRESS (foundation + predictor done)

**In one sentence:** Elements advance at their cluster's rate on one rank for a
fault-free problem, with cluster-boundary coupling by Taylor-integration, per
the Normative Scheduling and Buffer sections.

> **Progress (2026-07-18).** The STEPPING CORE is DONE + locally validated + reviewed:
> - `lts_time_basis.hpp` (`IntegrateTaylor`, A.4 — B.3 11/11).
> - `lts_stepper.{hpp,cpp}`: `LtsDkStore` / `LtsAccumulateBuffers` storage + the
>   pure tick loop `RunSyncInterval` (A.6 data + Normative Scheduling — B.4 31/31).
> - Per-cluster **predictor** `ComputeADERSubStepStatesAndIntegralCluster` (A.5):
>   element-restricted CK kernels (`ApplySpatialDerivativeElems_`,
>   `ApplyElementJacobianElems_` + BimaterialWaveOperator override) + D(k)
>   retention — **single-cluster == GTS BIT-for-bit, scalar + bimaterial**.
> - Per-cluster **corrector** `AdvanceADERClusterBulk` (A.6): element-restricted
>   volume (bit-exact) + role-driven face sweep (GTS / Boundary / ProviderSkip /
>   ConsumerFine with the coarse neighbour's `IntegrateTaylor`'d sub-interval
>   state + accumulate buffers) + restricted mass-inverse + in-place `Q +=` —
>   **single-cluster == GTS to machine epsilon (~1e-15)**, scalar + bimaterial.
> - **Multi-cluster consumer/buffer path** validated: a degenerate 2-cluster run
>   (both driven at dt, `[0,dt]` sub-intervals) reduces the coupling to GTS and
>   reproduces `AdvanceADER` to machine epsilon — 4 consumer faces + 4 provider
>   elems live, exercising D(k) retention, the ConsumerFine `IntegrateTaylor`
>   flux, buffer fill/consume, ProviderCoarseSkip, and buffers-zero-at-sync.
> - **End-to-end** `RunSyncInterval` wiring predict+correct: single-cluster LTS
>   over 6 steps == GTS. All in `test_lts_predictor` (19/19).
> - Reviewed twice (adversarial): fixed the bimaterial star-matrix override + the
>   accumulate-buffer fill-counts-sub-steps semantics (C-1).
>
> REMAINING for full Phase 2: the **driver sync-interval loop** (wire the core
> into `spatial_dyn_driver.cpp` — note the driver always carries a fault, so a
> fault-free bulk run needs a fault-free config; the fault half is Phase 3), the
> **multi-cluster conservation tests** (B.5 ragged / B.6 periodic-box / B.7
> mixed-neighbour — these exercise + validate the consumer/buffer path, which the
> single-cluster gates do not), and **Checkpoint V2**. The mathematically hard
> core (CK-recursion element-restriction + the role-driven Taylor-coupled face
> sweep) is complete and reproduces GTS.

> **Unblocked by the resequencing.** Phase 2 is bulk-only and np=1, so it needs
> NONE of the deferred trio — it consumes exactly the layout Phase 1 shipped. One
> caveat: the Checkpoint-V2 hash specifies "serial-mesh element order"; at np=1 use
> the local element order as a documented stand-in and re-base the hash on true
> serial order when serial clustering lands (Phase 4).

### Files to Create
- `dynamic/lts_stepper.hpp/.cpp` — the tick loop (Normative Scheduling;
  interface Appendix A.3).
- `dynamic/lts_time_basis.hpp` — `IntegrateTaylor` (Appendix A.4).

### Files to Modify
- `dynamic/wave_operator.{hpp,inl}` + `bimaterial_wave_operator.inl`:
  1. Per-cluster predictor `ComputeADERSubStepStatesAndIntegralCluster`
     (Appendix A.5): iterates the cluster's element index list; writes into
     the full-size `Q_per_node`/`I` vectors **zeroing ONLY the cluster's dof
     blocks** (the existing whole-vector zeroing at `:1634-1641` must not be
     reused); tau_nodes are the cluster's O midpoints on [0, dt_step];
     **RETAINS D(k) (MANDATORY, not optional) for provider elements** by
     copying out of the ping-pong scratch during the recursion (the scratch is
     shared across cluster invocations — no aliasing), raw/unscaled, epoch++.
     On the GTS path the retention branch is compiled but never taken
     (byte-exact).
  2. Per-cluster corrector `AdvanceADERCluster` (Appendix A.6): **in-place
     single-Q contract (P-010)** — ONE persistent state vector; only the
     cluster's dofs are read/written (all NUM_STATE strided blocks); the GTS
     Q→Q_new double buffer is not used on the LTS path. Order: volume + owned
     faces (role-driven sweep) + accumulate-buffer add-and-zero → per-element
     M⁻¹ → `Q +=`. All dt reads take the threaded `dt_step` (P-004).
  3. Consumer-face flux: `IntegrateTaylor(a, b, D_coarse)` per the t_origin
     formula; buffer fill per the Buffer design.
- `drivers/spatial_dyn_driver.cpp` — sync-interval outer loop replacing the
  nsteps loop when LTS is on (`while (t_s < tfinal)`); `step` becomes the sync
  counter, printed as "sync". Step-keyed consumer inventory (P-023): V_max
  print (per sync), NaN tripwire (per sync collective + optional per-tick
  rank-local isfinite spot check, collective-free), SCEC station writers (per
  sync), checkpoint-every (sync count), DIAG gates (sync count + time
  windows).
- `io/tpv104_checkpoint.hpp` — **V2 schema (P-018):** new magic tag
  `TPV104_CHECKPOINT_V2`; stores (t, sync_step, lts_mode:int, cluster-layout
  hash, partition hash). Hash = 64-bit FNV-1a over: rate (int32 LE),
  num_clusters (int32 LE), the cluster-id sequence in serial-mesh element
  order (int32 LE each), then the IEEE-754 bit patterns of dt_base and λ.
  Reader: V1 file + lts≠off → refuse (named abort); V2 + lts=off → refuse;
  V2 hash mismatch → refuse printing stored vs computed. GTS runs continue to
  write V1 byte-identically. `dof_data` on disk in canonical order — the
  permutation table for that lands with the reorder in **Phase 3 step 0**; Phase 2
  is bulk-only (no fault `dof_data`), and the hash's "serial-mesh element order"
  uses the np=1 local-order stand-in until serial clustering lands (Phase 4).

### Acceptance Criteria
- [ ] `test_lts_time_basis` (B.3), `test_lts_scheduler` (B.4 — golden tick tables + epoch/zero-at-sync property checks).
- [ ] Single-cluster LTS == GTS **byte-identical** (np=1, `lts_wiggle="off"` pinned — P-022), including tfinal=3.5·dt (truncated final step byte gate, P-004).
- [ ] Ragged-final-cycle test (B.5) green.
- [ ] Multi-cluster vs GTS on a smooth wave problem: truncation-order L2 agreement (2-resolution convergence study).
- [ ] **Conservation harness (P-012, B.6):** periodic Cartesian tet box (`Mesh::MakePeriodic`): ∫ρv_i (3) and ∫σ_ij (6) drift < 1e-12 × initial norm per sync interval; energy monotone-decreasing. (Fallback: traction-free box, ∫ρv_i only. Absorbing boundaries excluded — nothing is conserved there.)
- [ ] Mixed-neighbor 3-cluster chain test (B.7); 2×2 {deriv-cache}×{face-cache} smoke identical to 1e-15 (P-016).
- [ ] `test_lts_checkpoint_v2` refusal paths (B.8).
- [ ] `lts="off"` still byte-identical everywhere.

**Estimate:** 3 weeks (the core).

## Phase 3: Fault (dynamic rupture) under LTS (np=1)

**In one sentence:** Fault faces advance at their cluster's rate — per-cluster
friction sweeps over cluster-contiguous QP ranges, nucleation via the per-kind
absolute forms.

### Step 0 — Fault-QP reorder + early canary (MOVED from Phase 1, P-006)

Do this FIRST, in isolation, and validate it with the still-GTS canary gate
BEFORE wiring any range-sweep — so a reorder bug and a stepper bug cannot mask
each other.

- **`dynamic/wave_operator.{hpp,inl}` (+ `bimaterial_wave_operator.inl`):** apply
  the cluster-contiguous sort (by `(cluster, global face id)`, per-face QP blocks
  intact) to the wave operator's fault-face list at construction, gated on
  `lts != "off"`, BEFORE `SetFaultDOFData`. Everything downstream derives from
  this ONE ordering: station `Open`, nucleation-object construction, impedance/
  resolver seeding, `fault_coords`, and the first ParaView write. Assert
  `n_shared_fault_qps == 0` whenever the reorder is active (D-2 invariant).
- **Cluster ids at construction time (np=1):** cluster in the driver window
  BETWEEN material construction (`spatial_dyn_driver.cpp:~1502`) and operator
  construction (`~1569`), and pass the ids into the operator ctor. This needs NO
  material-before-ParMesh refactor at np=1 (that refactor is Phase 4).
- **Checkpoints:** serialize `dof_data` in CANONICAL (pre-LTS) order in both
  formats via the permutation table (the permutation's sole consumer; on-disk
  order stays layout-independent). Diag prints that name DOF indices note the ids
  now change meaning; log tooling must not assume stable ids.
- **Canary acceptance (the re-inserted Phase-1 gate):** with the reorder active
  but stepping STILL GTS, `lts="rate2"` TPV104-spatial station `.dat` files are
  **byte-identical** to `lts="off"`, and fault VTKHDF fields identical up to the
  known row permutation. This is an np=1 permutation-invariance test — runnable
  locally on a small fault fixture; it catches any missed reorder consumer (P-006)
  before the range-sweep below can hide it.

### Files to Modify
- `dynamic/friction_substep_iterator.{hpp,cpp}` — range-based
  `Advance(qp_begin, qp_end, dt_step, …)` (Appendix A.7): global indexing,
  base pointers + range; `I_imp_*_flat` outside the range untouched; per-range
  deltaT (dt_step/O each) with the Σ==dt_step check per invocation;
  `SetSubStepFaultImposedStates(qp_begin, qp_end, …)`; per-range
  ImposedGuard reset.
- `drivers/spatial_dyn_driver.cpp` — the fault half per fault-bearing due
  cluster: per-cluster tau_nodes on [0, dt_step]; per-cluster
  `EvaluateBulkAtFaultQPsCanonical` restricted to the cluster's fault faces
  (both elements same-cluster by the diff=0 clamp ⇒ time-consistent traces by
  construction; under D-2 the embedded exchange runs in `no_exchange` mode —
  P-007); per-sync reduction of `slip_rate_substep_max` feeding the output
  regime detector (GAP-C2/C4 resolved).
- `dynamic/spatial_nucleation.{hpp,cpp}` — D-3 range overloads for BOTH
  gradual kinds (new code; GTS path untouched).
- R-101: run at EVERY sync with `t_sync ≤ T_nuc_s` of the ACTIVE nucleation
  kind, plus sync 0 (no %100 throttle); under D-2 it is vacuous — retained as
  a locality tripwire (assert zero shared fault faces when lts≠off), with the
  np=2==np=1 per-sync fault-state checksum as the real cross-rank gate
  (P-019).

### Acceptance Criteria
- [ ] Single-cluster LTS == GTS byte-identical WITH fault (TPV104-spatial smoke).
- [ ] `test_lts_friction_range` (B.9): per-QP results equal a standalone-advanced reference per QP with its own dt.
- [ ] `test_lts_nucleation_absolute` (B.10): telescoping + partition-independence + range-restriction equivalence to 1e-15.
- [ ] Multi-cluster vs GTS on TPV104-spatial 200 m: station arrivals within 1%, final slip within 1%, no spurious V_max transients at fault-cluster edges.
- [ ] SAFS coarse smoke: LTS breakout time within 2% of GTS.

**Estimate:** 1.5–2 weeks.

## Phase 4: MPI

**In one sentence:** Every rank executes the same tick schedule, collectives
are counted and gated globally, and parity with single-rank LTS is the gate —
first with the simple exchange (4a), then the EDGE payloads (4b).

### Step 0 — Pre-ParMesh clustering pipeline (MOVED from Phase 1)

The np>1 machinery lands here, where rank-count-independent cluster ids first
matter (the np=2==np=1 gate) and where the LTS partition is first exercised.
Land this BEFORE 4a.

- **Material-before-ParMesh refactor:** construct (or make serially evaluable) the
  `MaterialField` before the `ParMesh` ctor (`spatial_dyn_driver.cpp:1310`), so
  rank 0 can compute per-element `dt_e = h_e / c_p,e` on the SERIAL mesh. Keep the
  existing post-ParMesh material path byte-identical for `lts="off"`.
- **Serial-mesh clustering (rank 0):** run `BuildLtsClustering` on the serial mesh
  (connected ⇒ contiguous ids), giving deterministic, rank-count-independent ids;
  broadcast/scatter to ranks; the local ids that Phase 1 computed per-rank at np=1
  are replaced by this authoritative serial result.
- **LTS-aware partition (P-009):** rank 0 builds the element graph
  (`Mesh::ElementToElementTable`), union-find-merges each fault-face element pair
  (as `fault_locality_partition.hpp`), calls `METIS_PartGraphKway` with
  `ncon = num_clusters` (`vwgt[v·ncon+c] = cellCost if cluster(v)==c else 0`),
  `ubvec = 1.05` per constraint; on METIS failure or
  `lts_partition_weights="scalar"`, single-constraint fallback
  `vwgt = cellCost·2^(maxC−c)`. Inject as `part_data` at the `ParMesh` ctor. The
  METIS call lives in a companion `.cpp` (the `metis.h ::real_t` collision — see
  the fault-weighted-partition plan).
- **Cross-rank maxdiff:** the serial fixpoint already gives global maxdiff≤1; the
  scatter preserves it, so no per-rank re-clustering (which would disagree at
  seams). Re-base the Checkpoint-V2 hash on the true serial-mesh element order
  (replacing Phase 2's np=1 local-order stand-in).

### Detailed Requirements
1. **4a:** full ghost-field exchange per due tick, per the Matched-collectives
   constraint (counts from the tick table; not-due ranks with shared faces
   still call). Debug: NaN-poison ghost slots of non-due clusters' elements
   before consumption (stale-ghost tripwire, P-007).
2. **4b (D-7 target):** per cluster-boundary ghost element, three payloads —
   (i) summed time-integral buffer, (ii) retained D(k) stack (provider cells),
   (iii) flux-premultiplied per-face block (n_seam_faces × NUM_STATE ×
   nbf_face; exact layout written into `lts_layout.hpp` before coding).
   Byte-compared against 4a on a 2-cluster np=2 box. **Phase-5 measures 4b.**
3. Rank-seam consumer faces (bulk only — fault is rank-interior): t_origin per
   the closed-form formula; no messages needed; interval-mismatch audit test.
4. **Collective audit table (required deliverable, P-007):** every collective
   in the macro step — predictor exchanges, corrector's 9 per-component
   exchanges, V_max/NaN reductions, R-101 — marked ELIMINATED (why safe) or
   MATCHED (schedule guaranteeing identical counts); per-sync debug counter
   covers ALL of them.
5. Checkpoint V2 at sync points; restart requires identical np AND matching
   layout hash.

### Acceptance Criteria
- [ ] np=2 LTS == np=1 LTS to 10 digits over ≥12 sync intervals (LSW+RS, mixed-rank fault distribution).
- [ ] np=10 symmirror gate green; locality tripwire green through nucleation.
- [ ] np=2 3-cluster rank-seam chain with ghost poisoning green (B.11).
- [ ] No hang in a 30-min np=10 SAFS coarse smoke (watchdogged).
- [ ] 4b byte-identical to 4a on the np=2 box; mid-campaign LTS restart (np=4) bit-continues.
- [ ] np>1 single-cluster byte gate: LTS(1 cluster) vs GTS **given the LTS partition** (P-022).

**Estimate:** 2.5–3 weeks.

## Phase 5: Performance validation

**In one sentence:** ≥15× end-to-end on the v4_0_0 coarse ALT case, measured on
the 4b exchange, with the achieved-vs-predicted gap explained.

1. p1 speed deck with `lts="rate2"`: **≥15×** vs the 2.60 sim-s/hour
   GTS-safety-1 baseline (ideal 38.7×; gate leaves room for overhead).
   Achieved vs Phase-0-predicted within 2× or explained (imbalance/overhead
   named).
2. Caliper `lts.cluster` attribute; per-sync aggregates.
3. LTS-weighted partition A/B (multi-constraint vs scalar vs none).
4. p3 deck rerun; document the residual kernel-efficiency gap (out of scope).
5. **D-1 experiment (gated, non-blocking):** `lts_fault_maxdiff=1` — friction
   at the finer side's rate; coarse-side fault trace integrated from retained
   D(k) per fine sub-step; TPV104-spatial tolerance gates re-run; results in
   the Phase-5 report; shipping default stays 0.

**Estimate:** ~1 week + cluster time.

## Phase 6: Default flip (D-5)

As rev 2, plus: the flip commit updates the sbatch pre-flights to assert the
lts mode in the log (the cfl_dg_safety pattern), and the flip is legal only
with the Phase-5 report accepted. `lts≠off` implies fault-locality (P-017), so
bare runs need no extra flags post-flip.

## Phase 7 (post-flip): far-field p-drop

Unchanged from rev 2 (driver-level two-order-class layout; own plan document
when Phase 6 lands).

## Appendix A — NORMATIVE interfaces

### A.1 Clustering (`dynamic/lts_clustering.hpp`)
```cpp
struct LtsClusteringOptions {
   int    rate = 2;
   int    max_clusters = 32;      // hard ceiling on raw bins (pre-merge)
   bool   wiggle_scan = true;     // lambda in (0.5,1], step 0.01; false => lambda fixed
   double lambda_fixed = 1.0;     // used when wiggle_scan == false
   int    nc_cap = 6;             // auto-merge target; <=0 disables merge (raw mode)
   double merge_loss_tol = 0.05;  // accept merge while cost <= (1+tol)*cost(uncapped)
   const std::vector<double>* cell_cost = nullptr; // nullptr => uniform 1.0
};
struct LtsClustering {
   std::vector<int> cluster;   // post-clamp, post-merge; 0-based contiguous
   int    num_clusters;
   double dt_base;             // = lambda * min(dt_e)
   double lambda;
   double modeled_cost;        // sum_e cellCost_e / (2^{c_e} * dt_base), post-clamp
};
LtsClustering BuildLtsClustering(const std::vector<double>& dt_e,
                                 const mfem::Table& elem_to_elem,
                                 const std::vector<std::pair<int,int>>& fault_face_elem_pairs,
                                 const LtsClusteringOptions& opt = {});
```
Binning (per λ candidate; NEVER floor/log2):
```cpp
int c = 0; double edge = lambda * dt_min;
while (2.0*edge <= dt_e[i] && c < opt.max_clusters-1) { edge *= 2.0; ++c; }
```
Pipeline per λ: bin → maxdiff fixpoint (diff=1; fault pairs diff=0) → cost.
Auto-merge: remove the top level one at a time while
cost ≤ (1+tol)·cost(uncapped), re-running cost each step; merge never violates
maxdiff (ids only decrease). After selection: production-path assert
`dt_base·2^{c_e} ≤ dt_e[i]` for every i, reusing the identical comparison
expression (GAP-A1).

### A.2 Layout (`dynamic/lts_layout.hpp`)
Per rank: `struct Cluster { std::vector<int> elems; std::vector<FaceRole> faces;
std::vector<int> fault_faces; };` FaceRole ∈ {IntraClusterGTS, ConsumerFine,
ProviderCoarseSkip, Boundary, Fault}. Element sets: `provider_elems`,
`consumer_owner_elems`, dense maps `provider_slot_of_elem`,
`buffer_slot_of_elem`. Global metadata (broadcast): per-cluster global element
count, global fault-face count. Built once from the serial-mesh clustering +
partition; pure function of (mesh, material, config).

### A.3 Tick table (`dynamic/lts_stepper.hpp`)
```cpp
struct LtsTick {
   std::vector<int>    predict_clusters;   // due this tick
   std::vector<int>    correct_clusters;   // due this tick, sorted FINE->COARSE
   std::vector<real_t> dt_step;            // per cluster id, current-step length
   int                 n_collectives;      // asserted by the debug counter
};
std::vector<LtsTick> BuildTickTable(int num_clusters, real_t dt_base,
                                    real_t T_actual, int ader_order,
                                    const LtsGlobalMeta& meta);
```
Regenerated per sync interval (identical on every rank; pure arithmetic).

### A.4 Taylor integration (`dynamic/lts_time_basis.hpp`)
```cpp
// a,b relative to the PROVIDER's expansion point (its last predict time).
// coeff[k] = (b^{k+1} - a^{k+1}) / (k+1)!   -- stacks stored RAW/unscaled.
void IntegrateTaylor(real_t a, real_t b,
                     const real_t* dk,     // [k][comp][i], k = 0..order-1
                     int order, int ndof,
                     real_t* out);         // [comp][i], comp = 0..NUM_STATE-1
```

### A.5 Per-cluster predictor
```cpp
void ComputeADERSubStepStatesAndIntegralCluster(
    const LtsLayout::Cluster& cl, real_t dt_step, int order,
    const Vector& Q,                       // read-only
    const std::vector<real_t>& tau_nodes,  // O midpoints on [0, dt_step]
    std::vector<Vector>& Q_per_node,       // full-size; ONLY cl dof blocks zeroed+written
    Vector& I,                             // full-size; ONLY cl dof blocks zeroed+written
    LtsDkStore& dk_store);                 // providers copied out during recursion; epoch++
```
`LtsDkStore`: raw D(k), layout `[slot][k][comp][i]`, one epoch counter per
slot; sized `order × NUM_STATE × ndof_per_el × n_provider_elems` and logged
against the deriv-cache budget.

### A.6 Per-cluster corrector
```cpp
void AdvanceADERCluster(const LtsLayout::Cluster& cl,
                        real_t dt_step, int order,
                        Vector& Q,                    // IN-PLACE; only cl's dofs written
                        const Vector& I_cluster,      // cl's dof blocks valid
                        LtsAccumulateBuffers& acc,    // consumed+zeroed for cl's owner elems
                        const LtsDkStore& dk);        // read-only neighbor stacks
```
`LtsAccumulateBuffers`: `real_t[buffer_slot][comp][i]`, pre-M⁻¹ residual
units, add-then-zero on consumption, epoch/fill counters asserted.

### A.7 Per-cluster friction
```cpp
// Range in the cluster-contiguous global QP order; global indexing throughout.
void FrictionSubStepIterator::Advance(size_t qp_begin, size_t qp_end,
                                      real_t dt_step,
                                      /* existing args: Q_pointwise_plus/minus base
                                         pointers (global-indexed), I_imp_plus/minus_flat
                                         base pointers (untouched outside range),
                                         nuc range-apply callback, ... */);
void WaveOperator::SetSubStepFaultImposedStates(size_t qp_begin, size_t qp_end,
                                                const real_t* I_plus,
                                                const real_t* I_minus);
```
Nucleation: `ApplyGradualOverstressAbsolute(dof_data, params, T_nuc, t,
qp_begin, qp_end)` and the compact-circular twin — new range overloads;
existing full-vector forms delegate with (0, N).

## Appendix B — Unit-test matrix (all new tests; file → cases)

| # | Test file | Cases | Phase |
|---|---|---|---|
| B.1 | `tests/unit/test_lts_clustering.cpp` | integer-loop binning: dt_e exactly at a lower edge joins THAT cluster (dt_cluster == dt_e, inclusive boundary); one ulp below joins c−1; single-element mesh; all-equal-dt ⇒ 1 cluster; maxdiff fixpoint convergence on a chain; fault diff=0; **λ-scan:** synthetic fat bin just above an edge ⇒ scan picks λ<1, cost strictly < cost(λ=1), CFL assert holds for every λ in the grid, cost(λ=1, nc_cap≤0) == raw baseline; **auto-merge:** num_clusters ≤ nc_cap, per-step cost growth ≤ tol, ids contiguous, maxdiff preserved, nc_cap≤0 reproduces uncapped bit-for-bit | 0 |
| B.2 | `tests/unit/test_lts_layout.cpp` | face-role assignment on a 3-cluster chain; provider/consumer element sets; slot maps dense + complete; global metadata counts | 1 |
| B.3 | `tests/unit/test_lts_time_basis.cpp` | sub-interval integrals exact for polynomials to order 4; Σ sub-intervals == whole interval to 1e-15; a=0,b=dt equals the existing whole-step integral | 2 |
| B.4 | `tests/unit/test_lts_scheduler.cpp` | golden tick tables for {2,3} clusters × ticks_per_sync {4,8} (exact predict/correct sets per tick, FINE→COARSE order, correct at (t+1)%2^c==0); properties: each cluster corrects exactly 2^(maxC−c)/interval; instrumented mock stepper: no consumer reads a stale epoch; buffers zero at sync | 2 |
| B.5 | `tests/unit/test_lts_ragged_final.cpp` | 3-cluster np=1 chain, tfinal = 3.5·dt_coarse (coarse truncates; fine's last sub-step truncates): all clusters land exactly on tfinal; consumed buffer == Σ truncated sub-interval integrals to 1e-15; buffers zero after final sync; result vs GTS at truncation order; control tfinal = 4.0·dt_coarse | 2 |
| B.6 | `tests/unit/test_lts_conservation.cpp` | periodic tet box: ∫ρv_i and ∫σ_ij drift < 1e-12·norm per sync; energy monotone decay; fallback traction-free box: ∫ρv_i only | 2 |
| B.7 | `tests/unit/test_lts_mixed_neighbor.cpp` | 3-cluster chain where one middle element carries provider+GTS+consumer faces simultaneously; conservation < 1e-12; instrumented single-flux-evaluation count (each boundary face evaluated exactly n_substeps times per coarse step, consumed once per side) | 2 |
| B.8 | `tests/unit/test_lts_checkpoint_v2.cpp` | V1 + lts=rate2 → named abort; V2 with mutated hash → abort printing stored vs computed; V2 + lts=off → abort; V2 round-trip bit-continues; GTS V1 write byte-identical pre/post change | 2 |
| B.8b | `tests/unit/test_lts_checkpoint_v2.cpp` (reorder case) | dof_data on-disk canonical order under reorder (needs the Phase-3 permutation table) | **3** |
| B.9 | `tests/unit/test_lts_friction_range.cpp` | two clusters, ranges advanced with different dt: per-QP (psi, slip, V) equal a standalone per-QP reference advanced with its own dt to 1e-15; I_imp outside range untouched; Σ deltaT == dt_step check fires on mismatch | 3 |
| B.10 | `tests/unit/test_lts_nucleation_absolute.cpp` | per kind (gradual + compact-circular): absolute form at every GTS sub-step time == telescoped incremental sum to 1e-15·amp; two different time partitions give identical τ_nuc at common times; disjoint range applies == whole-vector apply bit-for-bit | 3 |
| B.11 | `tests/unit/test_lts_mpi_seam.cpp` (np=2) | 3-cluster chain crossing the rank seam: np=2 == np=1 to 10 digits; ghost slots of non-due clusters NaN-poisoned in debug — no consumption fires; per-tick collective counter == tick-table n_collectives on both ranks | 4 |
| — | extended existing | `test-ader-tpv102-smoke` (lts=off byte gate, every phase); single-cluster byte gates with `lts_wiggle="off"` pinned incl. truncated-tfinal (2/3); 2×2 lever smoke (2); TPV104-spatial stations byte-identical under lts="rate2"+still-GTS stepping — the reorder canary **(3, step 0)** (was tagged Phase 1 pre-rev-5) | — |

## Testing Strategy (summary)
Byte-exact ladder (lts=off; single-cluster==GTS incl. truncated tfinal) →
analytic units (B.3/B.4/B.5) → conservation (B.6/B.7) → physics tolerance
(TPV104 stations 1%, SAFS breakout 2%) → parallel parity (np=2 10-digit,
np=10 watchdogged, collective counters) → performance (Phase-0 prediction vs
Phase-5 measurement within 2× or explained). Note: the two headline byte gates
CANNOT catch scheduler/buffer bugs (single cluster has no buffers) — B.4/B.5/
B.7/B.11 are the load-bearing correctness tests for the LTS core.

## Risk Assessment

| Risk (plain language) | Severity | Mitigation |
|---|---|---|
| Wrong sub-interval / buffer fill / truncated-step scaling ⇒ silent non-conservation | HIGH | Normative scheduler + buffer lifecycle with asserted invariants; B.4/B.5/B.6/B.7; current-step-dt threading (P-004) |
| Fault-QP reorder misses a consumer ⇒ silently scrambled friction | HIGH | Single-source reorder rule + 8-consumer inventory (P-006); stations byte-identical **canary gate at Phase 3 step 0** (still-GTS, np=1 permutation-invariance) — re-inserted after the reorder moved out of Phase 1 |
| Reorder moved out of Phase 1 loses its early canary ⇒ a reorder bug hides behind a stepper bug in Phase 3 | MED | Phase 3 does the reorder FIRST, in isolation, validated by the still-GTS stations gate BEFORE any range-sweep is wired |
| Collective-count mismatch ⇒ np≥10 hang | HIGH | Global gating rule + per-tick counter asserts + ghost poisoning (P-007); collective audit table deliverable |
| D(k) memory on big meshes | MED | Providers only (thin inter-cluster shell); sized+logged vs deriv-cache budget; named abort over threshold |
| LTS load imbalance eats the speedup | MED | Multi-constraint METIS weights (**Phase 4 step 0**, moved from Phase 1); Phase-5 A/B quantifies |
| Default flip changes TPV104-spatial science | MED | D-5 two-stage flip + pinning + re-goldening in one commit |
| λ-scan/merge bugs mis-cluster silently | MED | Production-path CFL assert inside BuildLtsClustering; B.1 λ/merge cases; raw mode preserved for the SeisSol cross-check |
| Checkpoint incompatibilities mid-campaign | MED | V2 magic + FNV-1a layout hash (λ, dt_base included) + refusal paths; B.8 |
| Deterministic tick loop leaves MPI idle vs actors | LOW (v1) | Accepted; actor/overlap is the named v2 axis |

**Total effort estimate (rev 3): ~9–12 weeks** (Phases 0–5; the review added
~2 weeks of specified tests and invariants — bought back many times over in
un-debugged silent-non-conservation), plus the flip; Phase 7 separately
planned. The Phase-0 go/no-go still costs only days.
