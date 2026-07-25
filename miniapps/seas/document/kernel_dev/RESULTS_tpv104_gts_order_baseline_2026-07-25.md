# TPV104 GTS order baseline — MFEM p1/p2/p3 vs SeisSol o2/o3/o4

**Date:** 2026-07-25 · **Jobs:** 52495068 (p1), 52495069 (p2), 52495071 (p3), 52495072 (SeisSol ×3),
52495014 (SeisSol o2/o3 build) · all COMPLETED, exit 0 · **6/6 legs, 0 gate FAILs, 0 dt-check FAILs**
TPV104 200 m, 2,464,689 tets, t = 0→2 s, CFL 0.5, **GTS both codes**, IO off both sides.

## The headline: the gap is not constant — it doubles with every order

| matched order | MFEM s/sim-s | SeisSol s/sim-s | **MFEM / SeisSol** | MFEM µs-core/upd | SeisSol µs-core/upd |
|---|---:|---:|---:|---:|---:|
| p1 ↔ o2 | 199.4 | 52.3 | **3.81×** | 17.68 | 4.64 |
| p2 ↔ o3 | 1182.0 | 156.9 | **7.53×** | 62.93 | 8.35 |
| p3 ↔ o4 | 4740.2 | 315.9 | **15.00×** | 180.25 | 12.01 |

Max-vs-max and mean-vs-mean give the **same ratios to 3 s.f.** (rank imbalance Max/Avg = 1.000 on
every MFEM leg), so the headline is not an artifact of statistic choice.

## Why — within-code order scaling, step-count penalty divided out

| | per-step cost vs its own p1 | µs-core/upd vs its own p1 |
|---|---:|---:|
| **MFEM** p1→p2→p3 | 1.00 → 3.56 → **10.19×** | 1.00 → 3.56 → **10.19×** |
| **SeisSol** o2→o3→o4 | 1.00 → 1.80 → **2.59×** | 1.00 → 1.80 → **2.59×** |

Going p1→p3 multiplies modes/element by 5 (4→20). **SeisSol's per-step cost grows 2.59×** — far
*sub*-linear in modes, because its generated GEMM kernels vectorise better as the blocks get bigger
(measured GFLOP/s/core **rises** 2.62 → 4.54 → 7.84). **MFEM's grows 10.19×** — *super*-linear,
consistent with O(modes²) mode-to-mode work executed through generic `DenseMatrix`/AXPY paths.

10.19 / 2.59 = **3.94×**, which is exactly 15.00 / 3.81. **This is an algebraic identity, not a
cross-check**: 3.94 ≡ gap(p3)/gap(p1) by construction, so the relation holds for any four positive
numbers and adds no information beyond the three measured gaps.

**The entire order-dependence of the gap is in the per-step cost of the element update.**
*(Corrected 2026-07-25: an earlier revision read "…is kernel scaling. Nothing else in the comparison
varies with order." That is false. 3.94 is the order-dependent excess of **whole-step** MFEM cost over
**whole-step** SeisSol cost, and its attribution is UNMEASURED — it provably contains at least three
order-dependent non-kernel terms: the unconditional `SetCurvature` below, the MFEM-only per-step NaN
`Q.Norml2()` + `MPI_Allreduce` whose vector is 5× longer at p3 than p1, and the 256-vs-16 partition
halo whose face modes grow 3→6→10. Do not name it "kernel scaling" or assign it to Track B until the
p1+p3 FLOP-counter leg returns — see `PLAN_performance_program_2026-07-23.md` §Target, gate T3.)*

## What this confirms, and what it refutes

**Confirms M1/M3 of the mechanism ledger** (`ANALYSIS_seissol_mechanism_gap_2026-07-24.md`) — and
upgrades them from "verified in source" to "measured end-to-end". The advantage is in *how the
element kernels are executed and laid out*, and it compounds with order.

**REFUTES a carried belief.** Memory (`project_lts_plan_clustered_ader`) recorded *"MFEM p1 FASTER
per update (1.354×)"* from the Phase-0 era. This baseline measures **MFEM 3.81× SLOWER at p1**. That
figure must not be reused; whatever configuration produced it does not survive a fairness-gated
GTS comparison. The correct statement is: MFEM is slower at every order tested, least so at p1.

**Sharpens the target.** An MFEM kernel that merely scaled *like SeisSol's* from its own p1 would put
p3 at **1205 s/sim-s = 3.81× SeisSol**, not 15×. That is the arithmetic size of the order-dependent
slice — and it is a *ceiling on order-scaling work*, not an attainability claim.

> **CORRECTION (2026-07-25), load-bearing — the superseded figure must not be reused.** This
> paragraph previously read *"199.4 × 2.59 = **517 s/sim-s = 1.6× SeisSol**"*. That drops MFEM's own
> **step-count penalty** at p3 (dt ∝ 1/(2p+1) ⇒ **×2.333** more steps; measured 2.3329 MFEM /
> 2.3321 SeisSol) — a factor **both codes pay identically** and that no kernel change can remove. It
> also contradicted this document's own `15.00 = 3.81 × 3.94` identity. Correct two ways:
> 199.4 × 2.59 × 2.333 = **1204.9**, and 4740.2 / 3.934 = **1204.8**; 1204.8 / 315.9 = **3.81×**.
> The corrected endpoint is **2.33× larger** than the retracted one; had the program been re-targeted
> on 1.6×, every future kernel result would have read as a shortfall.

Which levers can claim any of that slice is **undecided**: SeisSol's per-step cost grows 2.59× while
its achieved rate grows ×2.99, so its FLOP *volume* also grows super-linearly (×7.75). Whether MFEM's
excess is FLOP **volume** or FLOP **rate** decides between B1/B2/B4 (traffic reduction) and M1/M3
(element-local dense blocks / per-element GEMM, currently in no phase). One `perf stat` leg at **p1
and p3** settles it and has not run.

## Cost (for future planning)

| leg | elapsed | core-h |
|---|---|---:|
| MFEM p1 / p2 / p3 | 8:54 / 41:22 / 2:40:41 | 38 / 176 / 686 |
| SeisSol ×3 | 18:14 | 78 |
| **total** | | **~978** |

Predicted 875–1,375 → actual 978. The estimate model (dt ∝ 1/(2p+1), cost ∝ modes^{1..2}) was sound;
p3 came in just under the low end.

## What this baseline CANNOT claim

> This is a COST baseline at NOMINALLY-matched order on ONE mesh, ONE node count, ONE problem, GTS
> only, one run per leg — it measures **NO accuracy**, so it cannot claim "time to solution at matched
> accuracy", any Pareto / "p3 is worth it" result, any scaling behaviour, and nothing about LTS or any
> other mesh.

**Disclosed MFEM-side confounds** (part of the gap, not DG arithmetic; restate wherever the ratios
are quoted):
- `pmesh.SetCurvature(order)` is applied **unconditionally to a straight-sided tet mesh**
  (`spatial_dyn_driver.cpp:1733`) — MFEM evaluates a 4/10/20-node H1 nodal transformation per QP
  where SeisSol is always affine. **This is itself order-dependent and therefore a candidate for part
  of the 3.94×** — it should be measured and, if confirmed, is a cheap fix.
- Instrumentation asymmetry: ~17 nested Caliper regions inside MFEM's timed region vs a bare
  start/pause pair in SeisSol's.
- MFEM-only per-step global NaN `MPI_Allreduce` across 256 ranks (not disableable).
- Partition asymmetry: 256-way vs 16-way (~2.5× more halo for MFEM). Largest effect at p1 — **do not
  lead with the p1 ratio.**
- SeisSol's o2 build pads 4 modes to vector width (NZ/HW = 55.5 %), so its p1 leg is its least
  efficient — another reason p1's 3.81× is the least reliable of the three.

## Next

1. **Measure the `SetCurvature` confound** — it is order-dependent, sits in the p3 ratio, and may be
   removable for a straight-sided mesh. Cheapest possible test of a real slice of the gap. **Adopted
   as plan item B0 with two arms**: an LTS λ=1 p3 arm (the deciding one — `use_face_cache_`'s only
   functional read is the *GTS* corrector, so LTS evaluates interior-face transformations every step)
   and a GTS levers-on p1+p3 arm (prices the confound as it sits inside the published 3.94×). A null
   in the GTS arm alone means "the levers already mitigated it", **not** "geometry is not in the gap":
   under `--deriv-cache` both volume paths are quadrature-free and `--face-cache` removes the
   interior-face transformations, leaving only fault/shared/boundary faces.
2. **B1/B2's ownership of the slice is UNDECIDED**, pending the p1+p3 FLOP-counter leg (above). The
   ~517 s/sim-s "prize" is retracted; the arithmetic ceiling is 1205 s/sim-s = 3.81×.
3. D1 (rank sweep) is unaffected by this result and still sizes the comm/M4 side.
