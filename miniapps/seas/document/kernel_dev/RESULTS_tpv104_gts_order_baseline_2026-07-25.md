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

10.19 / 2.59 = **3.94×**, which is exactly 15.00 / 3.81. **The entire order-dependence of the gap is
kernel scaling.** Nothing else in the comparison varies with order.

## What this confirms, and what it refutes

**Confirms M1/M3 of the mechanism ledger** (`ANALYSIS_seissol_mechanism_gap_2026-07-24.md`) — and
upgrades them from "verified in source" to "measured end-to-end". The advantage is in *how the
element kernels are executed and laid out*, and it compounds with order.

**REFUTES a carried belief.** Memory (`project_lts_plan_clustered_ader`) recorded *"MFEM p1 FASTER
per update (1.354×)"* from the Phase-0 era. This baseline measures **MFEM 3.81× SLOWER at p1**. That
figure must not be reused; whatever configuration produced it does not survive a fairness-gated
GTS comparison. The correct statement is: MFEM is slower at every order tested, least so at p1.

**Sharpens the target.** B1/B2 (predictor fusion → tiling) attack precisely the super-linear term.
An MFEM kernel that merely scaled *like SeisSol's* from p1 would put p3 at 199.4 × 2.59 = **517
s/sim-s = 1.6× SeisSol**, not 15×. That is the size of the prize, and it is a kernel prize.

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
   removable for a straight-sided mesh. Cheapest possible test of a real slice of the gap.
2. **B1/B2 now have a target and a prize** (2.59× vs 10.19× scaling; ~517 s/sim-s if closed).
3. D1 (rank sweep) is unaffected by this result and still sizes the comm/M4 side.
