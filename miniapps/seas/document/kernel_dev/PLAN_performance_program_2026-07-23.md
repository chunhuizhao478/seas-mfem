> # ⛔ SUPERSEDED (2026-07-25) — see `PLAN_performance_2026-07-25.md`
>
> The GTS order baseline measured the gap as **order-dependent** (3.81×/7.53×/15.00×), which this
> document's single absolute endpoint could not express. Direction is unchanged and carries forward;
> the TARGET is refactored. Kept as history.

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

> **Restated 2026-07-24 after A1 (scoped update, not a direction change) — governed by
> `ANALYSIS_seissol_mechanism_gap_2026-07-24.md`:** comm remains the largest line item (40 % of
> wall), but after A1 its only live lever is unsized. Next actions, evidence-ranked:
> **D1** rank-count sweep 256/128/64 + arrival trace + FLOP counters (one job — sizes A6, tests the
> SeisSol-16-rank hypothesis, decides memory-vs-FLOP);
> **D2** `--precomputed-face-fluxes` vs `--face-cache` A/B (our own never-measured SeisSol mechanism);
> **D3 = B1 now** (best-supported code change program-wide). B2 waits on D1's counter leg.

1. ~~Communication first~~ → **D1/D2 measurements + B1 now** (see restatement above).
2. **Kernels**, predictor before faces, cheapest lever first.
3. **Face-cache-for-LTS** inside Track B, behind a hard correctness guard.

Three independent reasons, none of which depends on a projection:

- **It is the largest single line item.** Exposed `MPI_Waitall` is 52.5 % of the LTS wall (measured,
  job 52344266) and 52.6 % on a different job (measured, job 52379131 Leg 5).
- **It is the only bitwise-safe track.** Merging exchange rounds is identity-preserving by
  construction; kernel work needs tolerance gates.
- **It is what makes LTS worth having at all.** This is the strongest argument and it is entirely
  measured: LTS cuts compute 3.21× (87 % of its 3.68× ceiling) yet leads GTS by only **1.13×**,
  because the wait eats the win. We have already paid for LTS; comm is how we collect.

## Target (refactored 2026-07-25 — SCOREBOARD ONLY; the direction above is unchanged)

> **What this replaces and why.** The target used to live in the appendix as one absolute endpoint:
> *"Endpoint ~2110 s/sim-s ≈ 8.8× SeisSol on kernel items alone; [6.2, 8.8]× if A6 delivers."* It was
> (a) a single absolute number, (b) summed from bench factors and two rows marked ASSUMED, (c) silent
> about polynomial order — which the 2026-07-25 GTS baseline now measures as worth **3.81×–15.00×** —
> and (d) it folded a GTS-measured insight into an LTS budget without saying so. It was also
> arithmetically broken (Appendix A). **This refactor changes what we score and how we quote it. It
> adds exactly one work item (B0) and re-gates nothing else**: D1/D2/D3=B1 sequencing, the Track A/B
> split, and every existing stage gate stand verbatim.

### The scoreboard: an absolute, measured, monotone cost, quoted AT an order

| scoreboard quantity | today | provenance |
|---|---:|---|
| **MFEM p3 GTS — per-update compute** | **180.25 µs-core/upd** (SeisSol o4: 12.01) | jobs 52495071 / 52495072; TPV104 200 m, 2.46 M tets, 256 ranks, t=0→2 s, CFL 0.5, IO off |
| **MFEM p3 LTS λ=1 — wall** | **3123 s/sim-s** = **1871 compute + 1251 exposed wait** | job 52422891; same mesh, same rank count, `lts_wiggle="off"` |

Why this and not a ratio: it is **one leg to measure**, it is **monotone in the thing anyone cares
about** (any real speedup moves it down), it is **not gameable** (no denominator to regress), it
carries **no GTS→LTS transfer** because each row is measured in the regime it is quoted for, and it
is **order-explicit by construction** because it is quoted *at* an order. Production order is **p3**
(`jobs/lts_phase5/tpv104_200m_lts_speed_expanse/tpv104_200m_lts_rate2.toml:52`, order-matched to
SeisSol o4).

**Quotation rule (enforced at review):** a bare `s/sim-s` or `×-SeisSol` figure with no **order**,
**GTS/LTS mode**, **rank count** and **mesh** attached is rejected. No number is portable across
order — the gap alone moves 3.94× across the three orders measured.

### The measured order law — a ROUTING DIAGNOSTIC, explicitly NOT the scoreboard

Per-step cost normalised to each code's own lowest leg, step-count penalty divided out (dt and step
counts agree between the codes to 1e-16 — the baseline's fairness gate, so the penalty cancels):

| | p1↔o2 | p2↔o3 | p3↔o4 |
|---|---:|---:|---:|
| MFEM per-step vs its own p1 | 1.00 | 3.56 | **10.19** |
| SeisSol per-step vs its own o2 | 1.00 | 1.80 | **2.59** |
| **D — order-dependent excess** (ratio of the two) | **1.00** | **1.98** | **3.94** |
| measured gap | 3.81× | 7.53× | 15.00× |

`gap(p) = 3.81 × D(p)` reproduces the measured 7.53× and 15.00× to 0.1 %. **This is an algebraic
identity, not cross-validation** — `D(p) ≡ gap(p)/gap(p1)` by definition, so the relation holds for
any four positive numbers and carries no information beyond the three gap ratios. Do not quote it as
"self-validating"; the 0.1 % is rounding.

**D is a diagnostic and is NOT a progress metric.** Three reasons, all of which would be found by any
reviewer in half an hour:

1. **D is non-monotone in wall time, and B1 is the counterexample.** B1 fuses whole-vector
   accumulate/zero sweeps — traffic **linear** in modes — against a super-linear total. Removing a
   modes-linear term reduces the p1 leg proportionally *more* than the p3 leg, so **B1 landing exactly
   as designed lowers p3 wall time and RAISES D.** A scoreboard that scores the program's
   best-supported code change (D3 = B1 now) negative on day one would be a de-prioritisation of the
   settled direction disguised as a metric.
2. **D is an aggregate whose attribution is UNMEASURED.** It is the order-dependent excess of
   *whole-step* MFEM cost over *whole-step* SeisSol cost. It is **not** "kernel scaling": it provably
   contains at least three order-dependent non-kernel terms — the unconditional `SetCurvature` on a
   straight-sided mesh (B0 below), the MFEM-only per-step NaN `Q.Norml2()` + `MPI_Allreduce` (whose
   vector is 5× longer at p3 than p1), and the 256-vs-16 partition halo (face modes grow 3→6→10).
   The results doc's sentence *"Nothing else in the comparison varies with order"* is false and is
   corrected there.
3. **The confound net sign is UNKNOWN.** MFEM's order-independent per-step overheads inflate its p1
   leg and therefore **depress** D; SeisSol's o2 build pads 4 modes to vector width (NZ/HW = 55.5 %),
   inflating *its* p1 leg and therefore **inflating** D. Neither has been sized. Do not argue the
   number up or down from either.

What D *is* good for, and it is genuinely new: the **entire** order-dependence of the gap sits in
**per-step cost**, not in step count (dt and step counts are identical between the codes to 1e-16),
so it is not something the time integrator or the CFL choice can be blamed for; it **bounds** how
much of the p3 gap any order-scaling lever can ever address; and it therefore **routes work** — at p3
the levers that pay are the super-linear ones. Its consequences for wall time are arithmetic on
measured numbers and live in Appendix C, not here.

> ### ⚠ GTS-measured, LTS-produced — the caveat that travels with every use of D
> The order law was measured **GTS on both sides**. Production is **LTS**. The baseline's own README
> states *"It does not generalise to LTS."* **No LTS number may be derived from D by arithmetic** —
> that is this plan's standing rule since A1 (extrapolating one lever's law onto another). The
> transfer has already failed once in this program in exactly this shape: the GTS proxy predicted the
> predictor at 56.2 % of compute, the LTS measurement returned **39.6–46.6 %** (a 1.21–1.42× miss on a
> stage-composition quantity), and the **LTS seam corrector (12.0 % of step) has no GTS counterpart at
> all**, so D contains zero information about it. The single corroborating data point — LTS p3 13.1×
> vs GTS p3 15.0× — is one point and is not a validation.

### Target gates (measured deltas only)

- **T1 — reporting rule, applies to every Track A and Track B leg.** Report the p3 per-update compute
  cost **and** the exposed wait from the *same* Caliper split (commit `884349a`), with order, mode,
  rank count and mesh. Wins are banked in the scoreboard quantity. Recording both is also how the
  compute/wait coupling gets measured for free (T4) instead of assumed.
- **T2 — B0 geometry probe.** First action. Gate stated in the work list below.
- **T3 — attribution prerequisite, BLOCKING on any ownership claim over D.** The FLOP-counter leg
  (prerequisite #3, still unrun) must run at **p1 AND p3**, not one leg, and must report FLOPs/update
  **and achieved GFLOP/s/core**. SeisSol's rate rises 2.62 → 4.54 → 7.84 (×2.99) across the same span
  while its per-step cost rises only 2.59×, i.e. its FLOP *volume* also grows super-linearly (×7.75).
  Whether MFEM's excess is FLOP **volume** or FLOP **rate** decides whether B1/B2/B4 (traffic
  reduction) can address D at all, or whether it belongs to M1/M3 (element-local dense blocks /
  per-element GEMM), which are **in no phase**. Until this returns, D is descriptive and no share of
  it may be assigned to an owner.
- **T4 — LTS transfer.** The first change landing ≥1.2× on an LTS compute stage reports Δcompute and
  Δwait from the same leg. That measured coupling — not an assumed one — is what licenses any
  statement about what a compute win does to the LTS wall.
- **T5 — resolution and repeats.** Every leg in the order baseline is a **single run, no repeats
  anywhere in the program**; the variance of D is unbounded from the data. Repeat at least the two p3
  legs before any bound built on them is quoted. An item whose effect is smaller than the (still
  unmeasured) run-to-run spread is judged on its **stage gate alone**.
- **T6 — accuracy, the largest unaddressed threat, retired cheaply.** This is a **matched-ORDER**
  comparison; **no accuracy has been measured anywhere in this program**, so nothing here licenses a
  time-to-solution, Pareto, or "p3 is worth it" claim. Zero-compute action: post-process the six legs
  **already on disk** (jobs 52495068/69/71/72) against TPV104's SCEC on-fault metrics (rupture arrival
  time, peak slip rate, final slip at the standard stations) and record per-configuration error. This
  **gates nothing and re-opens nothing**; it is a caveat on the denominator. If it ever shows MFEM p2
  matching SeisSol o4's accuracy, the governing D row is **re-read** (1.98 at p2, not 3.94 at p3) —
  re-read from the table, not re-argued.

*Measured fact recorded because a ratio-only scoreboard hides it, and T6 is what would adjudicate it:*
**MFEM p1 GTS = 199.4 s/sim-s vs SeisSol o4 LTS as-run = 239 s/sim-s.** Cross-order, cross-mode
(GTS vs LTS), cross-accuracy, and the 239 had IO inside its stopwatch — it therefore **licenses no
claim whatsoever** about which order to run. It is recorded, not argued from.

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
one functional read (`wave_operator.inl:5651`, re-verified 2026-07-25) whose sole caller is the GTS corrector, so the 1.31×
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
| **A6** | **rank rebalancing — SIZE BEFORE BUILDING (rides D1)** | bounds from 52472765: **1.07×** (wait→Min-rank floor) to **1.67×** (wait→wire-only), a 9× spread. Sizing gate: per-rank tick-arrival trace — **systematic** lateness → proceed; **jitter** → drop. |

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

*(Historical note: after A0, this section read A0's result as "wait scales with round rate" and
declared A1–A3 measurement-backed. A1 falsified that reading — see the banner above. If A3 is ever
revived it needs a liveness gate in addition to bitwise identity: the recorded R-1600
unmatched-collective hang, `wave_operator.inl:3661-3672`, produces no wrong bits — it produces no
bits.)*

### Track B — kernels

| # | work | gate (measured) |
|---|---|---|
| **B0** | **`SetCurvature` geometry probe — FIRST ACTION, order-dependent, cheap, untested** | ≥1.05× on p3 per-update compute in **any** arm; ≤1e-12 parity. Full gate below. |
| **B1** | fuse the predictor's whole-vector accumulate/zero sweeps | ≥1.2× on the predictor stage; ≤1e-12 parity |
| **B2** | cache-tiled fused predictor | ≥1.6× on the predictor stage |
| **B3** | face-cache extended to the LTS corrector (see guard below) | ≥1.05× end-to-end; ≤1e-12 parity |
| **B4** | volume stage tiling | ≥1.5× on that stage |
| **B5** | face-interior tables — **re-justify before building** | drop if <5 % of the step |

**B0 — the `SetCurvature` probe, in full.** `pmesh.SetCurvature(cfg.mesh.order)` is applied
**unconditionally** to a **straight-sided tet mesh** at `drivers/spatial_dyn_driver.cpp:1733`, so MFEM
carries a 4/10/20-node H1 nodal transformation where SeisSol is always affine. It is
**order-dependent**, so it sits inside D; it is the cheapest falsifiable action on the board; and it
has never been measured. Make the call conditional on an actually-curved mesh and A/B it.

- **Correctness precondition (not a gate — an abort).** `MFEM_VERIFY` that the mesh is genuinely
  straight-sided before bypassing. An order-p H1 nodal transform reproduces an affine map *exactly*,
  so parity must be at machine noise. **Anything above 1e-12 means the mesh is not straight-sided and
  the bypass is invalid** — that is a bad probe, not a failed one.
- **Arms.** (i) **LTS λ=1, p3** — the deciding arm. `use_face_cache_` has exactly one functional read
  (`wave_operator.inl:5651`) and its sole caller is the **GTS** corrector, so under LTS the
  interior-face transformations *are* evaluated every step and geometry order *can* cost. (ii) **GTS
  p1 + p3 with the published baseline levers** (`--deriv-cache --shared-ck-recursion --face-cache`) —
  this arm prices the confound **as it sits inside the published 3.94×**, which is the only way to
  close the disclosure in the results doc.
- **PASS** = ≥1.05× on p3 per-update compute in any arm, at ≤1e-12 parity. Book it in the scoreboard.
- **Classification (information, never credit).** The p1 delta must be strictly smaller than the p3
  delta for the win to be order-dependent, i.e. to sit inside D rather than beside it. The win is
  banked either way — the scoreboard is absolute.
- **NULL** = p3 per-update compute within ±3 % in the GTS levers-on arm. **This outcome is predicted,
  and it must not be misread.** Under `--deriv-cache` both volume paths are quadrature-free by
  construction (`wave_operator.inl:1062-1112`, `1222-1252` — no `CalcShape`/`CalcPhysDShape`/Jacobian/
  quadrature), and `--face-cache` removes interior-face transformations from the GTS corrector. In
  that arm `SetCurvature`'s per-step footprint is therefore confined to fault/shared/boundary faces
  (fault is 0.65 % of the step). A null there means **"the levers already mitigated it"**, NOT
  "geometry is not in the gap" — which is precisely why arm (i) is mandatory and is the arm that
  decides. Only a null in **both** arms closes the item and returns its share to the rest of D.
- **Cost and its caveat.** ~145 core-h for the GTS pair at t = 0.2 s (0.1 × (38 + 686) core-h × 2
  arms), plus one LTS leg. **t = 0.2 s is entirely pre-nucleation** (T_nuc = 1.0 s), where the dynamic
  half's per-QP rate-state solves are absent — so a short leg gives a *delta*, never a baseline, and
  any short-leg instrument used later must first be validated against the 2 s legs.
- **Recurring-instrument note.** A full re-measure of D costs p1 + p3 = **724 core-h per scored item**
  against 978 core-h for the entire baseline. That is why D is not an exit gate: an unaffordable gate
  is not run, and the program silently reverts to stage gates. The scoreboard (T1) costs one leg.

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
- **FLOP-counter legs at p1 AND p3** (widened 2026-07-25 by target gate T3 — one leg is not enough).
  Track B's mechanism (traffic reduction) and the "hand-tuned ceiling" argument both assume MFEM is
  memory-bound. Nobody has ever counted MFEM's FLOPs per element update. `perf stat` at two orders
  reports FLOPs/update **and achieved GFLOP/s/core**, which decides (i) whether B1/B2 can work at all
  and (ii) whether the order-dependent excess D is FLOP **volume** (attackable by Track B) or FLOP
  **rate** (a code-generator property owned by M1/M3, which are in no phase).

## Known unknowns (do not claim these are settled)

Verification of the previous plan's numbers returned **43 quantities that cannot be derived from any
existing artifact.** The ones that would change decisions:

1. ~~Whether exposed wait scales with exchange-round rate.~~ **RESOLVED by A0+A1 jointly: it does
   NOT — wait scales with SYNC rate.** A0 (rounds/sync fixed, syncs ×0.63) → wait ×0.55; A1 (syncs
   fixed, rounds/sync ×0.51) → wait ×1.03 (job 52472765). A0 alone confounded the variables; A1
   deconfounded them.
2. ~~The LTS stage split.~~ **MEASURED by A0** (predictor 39.6–46.6 % of compute, NOT the assumed 56.2 %; seam-corrector compute 12.0 % of step is unclaimed). **Track B gates must be re-derived against it.**
3. **MFEM's FLOP count *and* achieved rate, at more than one order.** → counter legs at p1 and p3
   (target gate T3). Blocking on any claim that Track B owns a share of the order-dependent excess.
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

### A. The retired endpoint — relabelled, not deleted (so it stays falsifiable)

> ~~"Endpoint ~2110 s/sim-s ≈ 8.8× SeisSol on kernel items alone; [6.2, 8.8]× if A6 delivers its
> unsized upper bound."~~ **RETIRED 2026-07-25.** Re-expressed in the order-law frame, it was the
> assertion that B1+B2+B3+B4+seam would remove **~55 %** of the measured order-dependent excess,
> leaving D ≈ 1.8 at p3. Kept in that form so it can be *checked* later rather than forgotten. The
> pre-A1 "~5×, maybe" stays withdrawn.

It is not merely superseded — **the table it came from does not add up**, which is why that table is
replaced wholesale below rather than patched (patching would import the breakage):
- A5 was listed as saving **0** while the running column dropped **2483 → 2078** (unexplained −405).
- The running column bottomed at **1107 s/sim-s (4.6×)** while the same section's prose read
  **2110 (8.8×)** with band **[6.2, 8.8]**. Neither 4.6× nor 6.2× is reachable from the other.

### B. Confidence inventory — no running column, no sum, no endpoint

Deliberately not a ledger: summing rows of unequal confidence into one number is exactly how the
retired endpoint was produced. Each row stands alone against the measured LTS baseline
(**3123 s/sim-s = 1871 compute + 1251 exposed wait**, job 52422891).

| # | item | evidence | status |
|---|---|---|---|
| **A0** | `lts_wiggle="off"` | **MEASURED 1.370× end-to-end** (4277 → 3123) | **banked** |
| ~~A1~~ | merge 125→64 rounds/sync | **MEASURED ZERO** (52472765: rounds ×1.954, wait ×0.971) | dead |
| ~~A2~~/~~A3~~ | 64→32; split-post overlap | devalued by A1's mechanism; wire is 0.15 % of the wait | dead |
| A5 | comm term in `compute_cost` | generalises A0; **no gain claimed at np=256** | valid, unsized |
| **A6** | rank rebalancing | bounds **[1.07×, 1.67×]** on the wait — a 9× spread the plan itself calls UNSIZED | **bounds only**; D1 sizes it |
| **B1** | predictor sweep fusion | contended microbench only, **not in-solver** | best-supported code change; stage gate ≥1.2× |
| B2 | cache-tiled predictor | bench range **2.3–8.3× contended**; a point estimate from that range is not a size | mechanism blocked on T3 |
| B3 | face-cache → LTS corrector | explicit **±2×** error bar | stage gate ≥1.05× e2e |
| B4 | volume tiling | **ASSUMED — no bench of this stage anywhere** | unsized |
| — | seam-corrector compute (12.0 % of step) | **ASSUMED — in no phase, no owner** | unsized |
| **B0** | `SetCurvature` bypass | **UNTESTED** — order-dependent, ~145 core-h to decide | first action |

**Two floors, unchanged by this refactor and still the only composed statements this program will
make:** perfect comm / untouched kernels floors at **7.8×**; perfect kernels / untouched comm floors
at **5.2×**. **Under ~5× needs both tracks.** **Parity with SeisSol is not on the table at matched
order.** The one route past that band is the far-field order drop, which trades accuracy and remains
a separate proposal.

Also not in any row, and it should be: the **structural** rank imbalance A0 found — 19–24 % of the
wait is spread that merging cannot remove and it is flat through rupture. Only rebalancing touches it.

### C. What closing the order-dependent excess D would mean — GTS, definitional arithmetic

Conversion rule: `MFEM p3 GTS = 4740.2 × D/3.94`, and `× SeisSol o4 = 3.81 × D`.

| D at p3 | MFEM p3 GTS | vs SeisSol o4 GTS (315.9) |
|---:|---:|---:|
| **3.94 — today** | **4740 s/sim-s** | **15.00×** |
| 2.0 | 2406 | 7.6× |
| 1.5 | 1805 | 5.7× |
| **1.0 — full closure** | **1205** | **3.81×** |

**Four statements travel with every use of this table. Without them it is misleading.**

1. **It is a definition, not a computed prize.** `D(p) ≡ gap(p)/gap(p1)`, so "D = 1 ⇒ 3.81×" is
   literally the sentence *"the gap reverts to its p1 value."* Nothing was predicted.
2. **3.81× is NOT a floor.** Order-independent kernel wins lower the p1 leg too, and B1 — whose target
   is modes-linear vector traffic — lowers p1 proportionally *more* than p3, i.e. **B1 landing lowers
   the wall clock and raises D**. The "floor" moves under the work scheduled first. This is why the
   body scores an absolute cost, not D.
3. **No rung has an owner until T3 (the p1+p3 FLOP-counter leg) returns.** SeisSol's per-step cost
   rises 2.59× while its achieved rate rises ×2.99 — its FLOP *volume* grows ×7.75, super-linearly.
   If MFEM's rate is order-flat, most of D is SeisSol's rising arithmetic intensity, a code-generator/
   vectorisation property that **traffic-reduction levers (B1/B2/B4) cannot produce**; the levers that
   would attack it are **M1/M3 (element-local dense blocks / per-element GEMM), which are in no
   phase**. Publishing rungs before that leg is the point-estimate-ahead-of-measurement habit that
   forced three rewrites.
4. **The confound net sign is unknown** (see the body): MFEM's fixed per-step overheads depress D;
   SeisSol's 55.5 %-padded o2 build inflates it. Neither sized. Also: **single run per leg, no repeats
   anywhere** — the variance of D is unbounded from the data (T5).

### D. LTS composition — deliberately NOT computed

D is a **GTS** quantity; the production term it would act on is the LTS **compute** term (1871 of
3123). **No LTS endpoint is derived here by dividing 1871 by D**, because (i) the baseline's README
says the measurement does not generalise to LTS, (ii) the one prior GTS→LTS transfer of a
stage-composition quantity missed by 1.21–1.42× (predictor 56.2 % predicted vs 39.6–46.6 % measured),
(iii) the LTS seam corrector — 12.0 % of step — has no GTS counterpart, so D contains no information
about the largest unowned LTS stage, and (iv) the exposed-wait response to a compute cut is
**unmeasured with unknown sign**: A0+A1 established that wait scales with **sync rate** (which no
kernel change alters), while the arrival *spread in seconds* plausibly shrinks with compute. T4 is
what measures it, on a leg the program is running anyway. Until then the honest statement is the two
floors in §B.

*Campaign note, recorded once:* SAFS production is ALT/PREF × CASE1/CASE2 × k-sweep, so the binding
constraint is eventually scenarios-per-allocation, not one scenario. Nothing in this plan is sized
against that, and every number here is TPV104 200 m / 2.46 M tets / 256 ranks — production meshes are
34.9–52.3 M tets (known-unknown #6, unchanged).
