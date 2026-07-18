# ANALYSIS — LTS method selection: survey of all families, MFEM fit, and the verdict

**Date:** 2026-07-18 · **Method:** 9-agent study — 4 literature sweeps (papers
read, not abstracts), repo-grounded architecture scoring, exact update-count
arithmetic on the real 1,319,294-cell histogram, and a 3-lens judge panel
(max-speedup / risk-and-correctness / architecture-and-future). All three
judges returned the SAME winner independently.

## The verdict (unanimous)

**Backbone: clustered rate-2 ADER-LTS — but the EDGE-improved package, not
plain SeisSol.** Plus one SeisSol-impossible multiplier staged behind it.

The winning package (each item is a published improvement over SeisSol's
production scheme, from Breuer & Heinecke IPDPS 2022, arXiv:2202.10313 —
EDGE beat SeisSol's own LTS by 1.26–1.48× with these):

| # | Feature | What it buys | Status vs SeisSol |
|---|---|---|---|
| 1 | λ-wiggle factor (grid-search rescale of cluster boundaries, λ∈(0.5,1]) | +17.5% algorithmic on LOH.3 by moving fat bins down a cluster | SeisSol: experimental only |
| 2 | Capped cluster count + cost-model auto-merge (Nc≈5–6) | ALSO the GPU fix: SeisSol's GPU LTS collapsed to ~1.3× (V100, 2021) because small clusters → small batches; few merged clusters preserve the big uniform batches our element-local forall GPU design needs | SeisSol: experimental |
| 3 | Fixed 3-buffer, flux-premultiplied exchange (send 9×flux-projected payloads, not raw derivative stacks) | 94–95% of theoretical LTS speedup realized at scale; static schedule = our matched-collective contract | EDGE-only |
| 4 | Per-cluster friction sweeps fed by the persisted CK coefficients | reproduces SeisSol's space-time-quadrature friction; our currently-discarded ping-pong CK buffers ARE this machinery | equal |
| 5 | Multi-constraint METIS weights (one balance constraint per cluster level, Rietmann-style) | stronger than SeisSol's scalar 2^-c weights; extends our existing ncon=2 partition plan | beyond SeisSol |
| 6 | **Phase 2: far-field p-drop (two order classes: p3 fault region / p1 far field)** | ~1.16× update-count + memory-bandwidth relief; **structurally impossible for SeisSol** (compile-time ConvergenceOrder) | MFEM-only axis |

**Projected outcome:** ~35× realized element-update reduction from the backbone
(38.73× ideal × the 94–95% realization EDGE demonstrates), ~45× with the
p-drop — a projected **1.2–1.5× faster than SeisSol on the same mesh**, on top
of our already-verified 1.354× per-update advantage.

## Why each alternative lost (arithmetic, not taste)

Exact update-count results on the real histogram (M2, all formulas shown in
the study transcript):

| Scheme | Ideal update-speedup vs GTS | Killed by |
|---|---|---|
| Clustered rate-2 | **38.73×** (fault 22.0×) | — the winner's baseline |
| Clustered rate-3 / rate-4 | 30.8× / 25.3× | strictly dominated by rate-2 on this histogram |
| Cluster caps Nc∈{4,6,8} (no merge tuning) | < rate-2 | dominated; the cap pays off only WITH auto-merge + λ |
| Elementwise (Dumbser 2007) | **55.9×** (fault 31.8×) — the highest ideal | throughput death: no batches, per-element scheduling, data-dependent async comm (exactly our historical hang class), single-dt friction contract dissolves; repo-fit score 9/25; every production code abandoned it |
| Locally implicit (clusters 0–3 implicit) | 16.2×-equivalent; break-even implicit cost **negative** for thresholds k≤5 | mathematically dead on this histogram: clusters 0–3 are 9.1% of the LTS harmonic work (Amdahl); worse, cluster 0 holds 49 DR faces ⇒ rate-and-state INSIDE the implicit solve (implicit RS documented non-convergent; all Maxwell-linear proofs void) |
| Hybrid implicit(0–3)+LTS | ≤ +6% (41.1×) | not worth a second solver framework |
| Leapfrog/Newmark-LTS (the only family with real stability proofs) | ~4× demonstrated ceiling (SPECFEM3D Tohoku 3.9×) | proofs are for mass-lumped conforming 2nd-order FEM — adopting it forfeits ADER, upwind, and order AND voids its own proofs on our discretization |
| AB-multirate / MRI-GARK | 2–3× demonstrated | replaces our time discretization entirely; no wave-DG production record; no theory for a 2^10 hierarchy |
| Tent-pitching (MTP) | conceptually ideal | immature: 2D/Maxwell scale, order-reduction fixes still landing, friction only in 2D aSDG; per-tent asynchrony is the antithesis of our friction sweep + collective contract |
| Parareal/PinT | — | provably poor for hyperbolic transport |

**The decisive structural fact:** clusters 4–5 hold **76.2% of the clustered
update cost and 96.7% of the DR faces**. Every exotic alternative attacks the
other 24%. The one real future-recovery axis is *finer-rate clustering inside
the two fat bins* (elementwise's residual 1.44× = 55.9/38.7) — reachable later
INSIDE the clustered architecture, not by switching families.

## Key literature facts worth keeping

- No GKS/energy stability proof exists for 2:1 upwind-DG subcycling in ANY
  family — clustered ADER-LTS's decade of petascale production **including
  dynamic rupture** (Uphoff SC'17, the only DR+LTS production record) is the
  strongest available evidence class. The wiggle factor doubles as the
  engineering analogue of avoiding leapfrog-LTS's discrete-dt resonances.
- Uphoff SC'17: ~10× from LTS on the Sumatra DR run; DR faces forced
  same-cluster there too — our 22×-vs-38.7× fault-face analysis appears novel;
  nobody has published per-face friction rates (a possible paper for us).
- ExaHyPE does NOT do clustered LTS (global/adaptive stepping + optimistic
  rollback); Salvus ships NO LTS (contrary to folklore); SPECFEM3D's
  Newmark-LTS peaked at 3.9×.
- Dumbser's p-adaptivity (with the max-degree zero-padded interface flux rule)
  is the published template for our two-order-class layout; ADER/CK is
  element-local so per-element order is natural.
- MFEM-library note (M1): variable-order FESpaces exist (4.9.1-dev,
  `SetElementOrder`) but require nonconforming meshes and parallel hp is
  quad/hex-scoped — **our conforming tet ParMesh is outside that envelope.**
  The driver owns its flat DOF layout, so the p-drop is implemented as a
  driver-level two-order-class layout (prefix-sum offsets, per-class kernel
  batches, max-degree interface rule) — multi-week, mechanical, and it shares
  its element-class plumbing with the LTS scheduler.

## Judge rankings (independent, all converged)

- **Max-speedup lens:** winner = phased roadmap (EDGE-package LTS now, p-drop
  Phase 2, severable) ≈ 42–43× realized; runner-up = LTS package alone (32–36×).
- **Risk lens:** winner = LTS package (only DR+LTS production evidence in
  existence; exact conservation by identical time-integrated fluxes; Nc=1
  degenerates to GTS = byte-exact gate; static exchange schedule); DR faces
  same-cluster initially, maxdiff≤1-across-fault only as a later gated
  experiment. Elementwise ranked near-bottom under this lens.
- **Architecture lens:** winner = LTS package (every requirement maps to an
  existing asset — discarded CK buffers, flat friction sweep, matched
  collectives, GPU batches; only backbone compatible with the planned
  plasticity/CDBM cross-rank collective); runner-up = p-order classes (the
  composable SeisSol-impossible axis).

## Changes applied to the implementation plan

`PLAN_clustered_lts_ader_2026-07-18.md` is amended (same file, revision 2):
1. **D-6 REVERSED:** λ-wiggle and the Nc-cap + auto-merge move INTO v1 (they
   are the GPU story and the realization-fraction story, not optimizations).
2. **Exchange design:** EDGE-style fixed 3-buffer flux-premultiplied payloads
   replace "exchange full ghost field per tick" as the Phase-4 target
   (full-field exchange stays as the v0 stepping stone).
3. **D-1 refined:** DR faces same-cluster (SeisSol/Uphoff rule) in v1;
   maxdiff≤1-across-fault becomes a named Phase-5 experiment (novel, potential
   publication).
4. **Partition:** multi-constraint (per-cluster-level) METIS weights as the
   Phase-1 target, scalar 2^-c weights as fallback.
5. **New Phase 7 (post-flip): far-field p-drop** — driver-level two-order-class
   layout, fault region p3, far field p1; gated by its own plan document.
6. Payoff restated: ~35× realized (backbone), ~45× with Phase 7; target
   1.2–1.5× faster than SeisSol o4 on the shared mesh.
