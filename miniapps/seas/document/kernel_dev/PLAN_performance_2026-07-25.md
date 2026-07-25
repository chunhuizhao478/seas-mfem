# MFEM performance program — plan (2026-07-25)

**Supersedes** `PLAN_performance_program_2026-07-23.md` (kept as history). **Direction unchanged** —
this is a *target* refactor forced by the GTS order baseline, not a re-litigation.
**Evidence:** `RESULTS_tpv104_gts_order_baseline_2026-07-25.md` (jobs 52495068/69/71/72),
`RESULTS_a0_wiggle_2026-07-24.md`, `RESULTS_a1_merge_2026-07-24.md`,
`ANALYSIS_seissol_mechanism_gap_2026-07-24.md`.
**Rule (unchanged, it has worked):** the body carries only measured quantities and gates expressed as
measured deltas; every projection lives in the quarantined appendix.

## 1. Where we are — measured, quoted at an order

| quantity | today | provenance |
|---|---:|---|
| **MFEM p3 GTS, per-update compute** | **180.25 µs-core/upd** (SeisSol o4: 12.01) | 52495071 / 52495072 |
| **MFEM p3 LTS λ=1, wall** | **3123 s/sim-s** = 1871 compute + 1251 exposed wait | 52422891 |

The gap is **not one number** — it is 3.81× (p1), 7.53× (p2), **15.00×** (p3).

**Quotation rule, enforced at review:** a bare `s/sim-s` or `×SeisSol` with no **order**, **GTS/LTS
mode**, **rank count** and **mesh** is rejected. Nothing is portable across order — the gap alone
moves 3.94× across the three orders measured.

## 2. The order law — a routing diagnostic, **not** a progress metric

Per-step cost normalised to each code's own lowest leg (step-count penalty divided out; dt and step
counts agree between codes to 1e-16, so it cancels exactly):

| | p1↔o2 | p2↔o3 | p3↔o4 |
|---|---:|---:|---:|
| MFEM per-step vs own p1 | 1.00 | 3.56 | **10.19** |
| SeisSol per-step vs own o2 | 1.00 | 1.80 | **2.59** |
| **D — order-dependent excess** | 1.00 | 1.98 | **3.94** |
| measured gap | 3.81× | 7.53× | 15.00× |

`gap(p) = 3.81 × D(p)` is an **algebraic identity, not a cross-check** — `D ≡ gap(p)/gap(p1)` by
construction. Do not quote it as self-validating.

**D must not be used as a scoreboard.** Three reasons:

1. **D is non-monotone in wall time, and B1 is the counterexample.** B1 removes traffic *linear* in
   modes from a super-linear total, so it lowers the p1 leg proportionally more than p3 — **B1
   landing exactly as designed lowers p3 wall time and RAISES D.** A metric that scores our
   best-supported change negative on day one is a de-prioritisation disguised as a number.
2. **D's attribution is UNMEASURED.** It is the order-dependent excess of *whole-step* cost and
   provably contains at least three order-dependent **non-kernel** terms: the unconditional
   `SetCurvature` (B0), the MFEM-only per-step NaN `Q.Norml2()` + `MPI_Allreduce` (vector 5× longer
   at p3), and the 256-vs-16 partition halo (face modes 3→6→10). **It may not be called "kernel
   scaling" until T3 returns.**
3. **The confounds' net sign is unknown.** MFEM's order-*independent* overheads depress D; SeisSol's
   o2 vector padding (NZ/HW 55.5 %) inflates it. Neither is sized.

What D *is* for: the entire order-dependence sits in **per-step cost**, not step count — so it cannot
be blamed on the time integrator or CFL — and it **bounds** what any order-scaling lever can address.
At p3, the levers that pay are the super-linear ones.

> **⚠ GTS-measured, LTS-produced.** The order law was measured GTS on both sides; production is LTS.
> **No LTS number may be derived from D by arithmetic** — that is A1's exact failure mode
> (extrapolating one lever's law onto another). This transfer has already failed once here: the GTS
> proxy predicted the predictor at 56.2 % of compute; the LTS measurement returned **39.6–46.6 %**,
> and the LTS seam corrector (12.0 % of step) has **no GTS counterpart at all**.

## 3. Target gates (measured deltas only)

- **T1 — reporting.** Every Track A/B leg reports p3 per-update compute **and** exposed wait from the
  *same* Caliper split (`884349a`), with order/mode/ranks/mesh. This also measures the compute–wait
  coupling for free instead of assuming it.
- **T2 — B0 geometry probe.** First action (§4).
- **T3 — attribution, BLOCKING on any ownership claim over D.** The FLOP-counter leg must run at
  **p1 AND p3** and report FLOPs/update **and** achieved GFLOP/s/core. SeisSol's rate rises ×2.99
  while its per-step cost rises ×2.59 — its FLOP *volume* grows ×7.75. Whether MFEM's excess is FLOP
  **volume** or **rate** decides between B1/B2/B4 (traffic) and M1/M3 (element-local dense blocks /
  per-element GEMM — **currently in no phase**). Until it returns, D is descriptive and no share is
  assigned to an owner.
- **T4 — LTS transfer.** The first LTS leg after any kernel change re-measures the stage split rather
  than inheriting a GTS share.

## 4. Work list

| # | item | gate |
|---|---|---|
| **B0** | **`SetCurvature` geometry probe — FIRST.** Unconditional on a straight-sided tet mesh (`spatial_dyn_driver.cpp:1733`): MFEM runs a 4/10/20-node H1 nodal transform per QP where SeisSol is affine. **Two arms:** (a) **LTS λ=1 p3** — the deciding one, since `use_face_cache_`'s only functional read is the *GTS* corrector, so LTS re-evaluates interior-face transforms every step; (b) GTS levers-on p1+p3 — prices it as it sits inside the published 3.94×. | p3 per-update compute falls. **A null in arm (b) alone means "the levers already mitigated it", NOT "geometry is not in the gap"** — under `--deriv-cache` both volume paths are quadrature-free and `--face-cache` removes interior-face transforms, leaving only fault/shared/boundary. |
| **T3 leg** | FLOP counters at p1 **and** p3 | blocks ownership claims over D |
| **D1** | rank sweep 256/128/64 + tick-arrival trace (sizes A6: bounds 1.07–1.67×) | unaffected by the baseline; still the only thing that sizes the comm side |
| **B1** | fuse predictor accumulate/zero sweeps (28 % of self-time is pure vector traffic) | ≥1.2× predictor stage; ≤1e-12 |
| B2/B4 | tiled predictor / volume | **gated on T3** — their mechanism assumes memory-bound |
| B3 | face-cache → LTS corrector (guard: `role == IntraClusterGTS`) | ≥1.05× end-to-end; ≤1e-12 |
| A6 | rebalancing | gated on D1 |
| ~~A1/A2/A3~~ | merge/overlap | **dead** — A1 measured zero; wait scales with sync rate |

## 5. Appendix — projections (not load-bearing)

**The order-scaling ceiling.** If MFEM merely *scaled* like SeisSol from its own p1:
199.4 × 2.59 × **2.333** = **1205 s/sim-s = 3.81× SeisSol** (cross-check 4740.2 / 3.934 = 1204.8).

> **RETRACTION — the superseded figure must not be reused.** An earlier revision published
> **"517 s/sim-s = 1.6× SeisSol."** That dropped MFEM's own **step-count penalty** at p3
> (dt ∝ 1/(2p+1) ⇒ ×2.333 more steps) — a factor **both codes pay identically** and no kernel change
> can remove — and it contradicted the `15.00 = 3.81 × 3.94` identity. The corrected endpoint is
> **2.33× larger**. Had the program been re-targeted on 1.6×, every future kernel result would have
> read as a shortfall.

This is a **ceiling on order-scaling work, not an attainability claim**, and it says nothing about
which lever earns it (T3 decides). Nothing beyond B0/T3/D1 is authorised by current evidence.
**Parity is not on the table** at matched order; the far-field order drop remains a separate proposal
that trades accuracy — and note MFEM is only **3.81×** behind at p1, which is what makes it
strategically interesting.
