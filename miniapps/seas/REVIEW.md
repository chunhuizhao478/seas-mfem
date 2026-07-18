# Plan Review: 2026-07-18 — PLAN_clustered_lts_ader (rev 2), adversarial audit

> Reviews the LTS implementation PLAN (not code). Two independent adversarial
> reviewers (implementability lens; numerics/MPI-traps lens) + the maintainer's
> own pass. Both reviewers verified the plan's file:line claims against the repo
> (all anchors correct); the defects are in what the plan leaves implicit or
> states inconsistently. Previous review archived at
> `REVIEW_sidecar_shared_mem_2026-07-17.md`.

## Review Scope
- Plan: `document/lts_dev/PLAN_clustered_lts_ader_2026-07-18.md` (rev 2)
- Companion: `document/lts_dev/ANALYSIS_lts_method_selection_2026-07-18.md`
- Repo verification: `drivers/spatial_dyn_driver.cpp`, `dynamic/wave_operator.inl`,
  `dynamic/friction_substep_iterator.*`, `dynamic/spatial_nucleation.*`,
  `io/tpv104_checkpoint.hpp`, `spatial/code/spatial_friction.hpp`
- Convergent findings from both reviewers are merged; IDs P-001…P-024.

## Findings (consolidated; severity ▸ id ▸ title)

### CRITICAL — scheduler & coupling correctness (the interlocking core)

**[P-001] The tick-loop pseudocode is wrong in both readings.** With
`due_clusters(tick) = {c : tick % 2^c == 0}` used for BOTH predict and correct,
a coarse `correct(c)` at tick 0 consumes an accumulate buffer holding 1 of 2^Δ
fine contributions (silent non-conservation from step one); the alternative
reading never corrects the last step of any cluster. BOTH reviewers found this
independently. Neither headline gate (lts=off byte; single-cluster==GTS)
catches it. **Fix:** predict due at `tick % 2^c == 0` (opens the step); correct
due at `(tick+1) % 2^c == 0` OR `tick+1 == ticks_per_sync` (closes the step),
corrects ordered FINE→COARSE within the tick; epoch-counter assert (buffer fill
count == number of fine sub-steps of the closing coarse step); buffers exactly
zero at every sync. New `test_lts_scheduler` with golden tick tables.

**[P-002] Two incompatible accumulate-buffer designs; "anti-symmetric" is
physically wrong.** Glossary describes SeisSol's state-buffer (fine sums ITS
OWN time-integrals; coarse still visits the face); Phase 2 describes a
flux-contribution buffer (coarse skips the face). Different owner, units,
consumption. And "anti-symmetric coarse-side contribution" is false on
bimaterial faces (A± differ per side) and under different test bases. **Fix:**
choose the flux-contribution design (matches D-7 EDGE premultiplied payloads):
buffer is per CONSUMER (coarse) element, NUM_STATE×ndof_per_el, **pre-M⁻¹
residual units**, filled by the fine side evaluating the coarse side's own
Godunov flux tested with the coarse basis; consumed-then-zeroed inside the
coarse correct, before its M⁻¹. Delete "anti-symmetric". Role-driven face sweep
(an element can simultaneously carry provider/GTS/consumer faces — sweep per
face role; buffer added exactly once per coarse step; per-element storage is
well-defined because maxdiff≤1 ⇒ all finer neighbors are exactly c−1). The
transplanted SeisSol "fifth rule" has NO REFERENT in this design (our buffer is
not the fine cell's state buffer) — deleted, replaced by the role table.

**[P-003] Sync semantics self-contradictory; ragged-interval contract missing.**
Glossary says truncate-at-sync; Phase 1 default sync = coarsest dt; the risk
table says "sync at max(coarsest dt, requested cadence)" — a third rule that
would make EVERY interval ragged for every cluster (1 s cadence is not a
multiple of ~42 ms) and, with the V_max-adaptive cadence, data-dependent.
**Fix:** one rule: the sync grid advances by dt_base·2^(Nc−1) ("auto");
truncation occurs ONLY on the final interval before tfinal; outputs are
evaluated AT sync points (the adaptive cadence selects which syncs write, fed
by the per-sync REDUCED V_max so all ranks agree); a numeric `lts_sync_dt` is
snapped down to a multiple of the coarsest dt with a log line. Risk-table
sentence deleted. Ragged-final-cycle test (tfinal = 3.5·dt_coarse) with
interval-sum and fill-count asserts.

**[P-004] Truncated-step dt at the corrector's dt sites.** "takes the
per-cluster dt of the owning element's cluster" is wrong on truncated steps:
the friction side channel scales by the ACTUAL step (`accum_scale =
weight·dt_macro`, verified), so a corrector using nominal dt_c mis-scales
`I_imp/dt` and `Q_imp·dt` by up to ~2× on the last step. **Fix:** every dt
read takes the CURRENT step length `dt_step(c,tick)` threaded as one argument
from the tick loop; no site recomputes it. Cheap byte gate: single-cluster LTS
with tfinal=3.5·dt vs GTS (which already truncates at tfinal).

**[P-005] t_origin formula + FP residual steps.** "a deterministic function of
the tick index" is asserted but never given, and holds only if origins reset at
every sync and all times come from closed forms. **Fix (normative):**
`t_origin(c,tick) = t_s + dt_c·floor(tick/2^c)`; `t(tick) = t_s + tick·dt_0`
(multiplication, never `t += dt`); fine sub-interval `[a,b]` with
`b = min(t(tick)+dt_fine, t_s+T_s) − t_origin` — then `b ≤` the coarse's
current step ≤ dt_c is PROVABLE. Assert tracked-vs-formula ≤1e-12·dt_c at
every consumption. Residual steps < 1e-10·dt_c merge into the preceding step
(decided from the closed-form schedule, identical on all ranks).

**[P-006] The cluster-contiguous reorder has ≥8 consumers; the plan names 2.**
Missing: the wave operator's canonical fault-QP index `dof_idx` (THE source of
truth for `substep_I_imp_*` and `Q±` traces — reorder dof_data without it and
every friction read is scrambled silently), `fault_coords`, the
interior-then-shared split invariant, index-keyed nucleation caches, the
ParaView fault writer geometry/fields, checkpoint payload order, impedance/
resolver seeding, diag DOF ids. Also an internal contradiction: V1-under-LTS
refusal vs "permutation table for checkpoint compatibility". **Fix:** the sort
happens in ONE place (the wave operator's fault-face list at construction,
before SetFaultDOFData/stations/nucleation/PV); everything derives. Checkpoints
serialize in CANONICAL (pre-LTS) order via the permutation in BOTH formats
(this is the permutation's consumer; V1-under-LTS still refused). Assert
`n_shared_fault_qps == 0` when the reorder is active (D-2). Acceptance: with
lts="rate2" under GTS stepping, station .dat files byte-identical to lts="off".

**[P-007] Per-tick collective schedule undercounted and ungated.** Reality
today: O + 9 collectives per macro step (O predictor-substep exchanges in
`EvaluateBulkAtFaultQPsCanonical` + 9 per-component corrector exchanges), not
"one exchange". The not-due-cluster rank rule is unstated — a rank whose shared
faces all sit in not-due clusters must STILL participate; gating on rank-local
data reproduces the R-1600 hang class. **Fix (normative rule):** collectives
gate exclusively on (a) the global tick table and (b) serial-mesh cluster
metadata (global per-cluster element/fault counts). Per tick:
`n_x = Σ_{c∈due} [9 + (cluster c has fault faces GLOBALLY ? O : 0)]`,
precomputed into the tick table, asserted by a per-tick debug counter on every
rank. Under D-2, shared fault QPs are globally zero ⇒ the predictor-substep
exchange is dropped GLOBALLY (config-level fact; `no_exchange` mode asserted
safe by fault locality). Debug mode NaN-poisons ghost slots of non-due
clusters' elements before consumption; np=2 3-cluster rank-seam test.

**[P-008] λ-wiggle cost model is vacuous as written; binning FP hazard.**
Without RE-BINNING against λ-scaled edges, Σ cost/(2^c·λ·dt_min) is monotone in
λ ⇒ the scan always returns λ=1. And `floor(log2(...))` is an FP determinism
hazard violating the plan's own deterministic-clustering constraint; the
bin-edge acceptance line ("joins the LOWER cluster") contradicts the binning
formula. **Fix:** normative binning `c(λ) = max{c : λ·2^c·dt_min ≤ dt_e}` via
the integer comparison loop (never floor/log2); cost evaluated after the
maxdiff fixpoint per λ candidate; the CFL assert `dt_cluster(e) ≤ dt_e` moves
INSIDE BuildLtsClustering (production path, not just the report tool), reusing
the identical comparison expression; bin-edge test corrected (dt_e at a lower
edge joins THAT cluster, dt_cluster == dt_e, boundary inclusive); auto-merge
direction invariant (merge only reassigns to smaller c). λ/merge unit tests.

**[P-009] Phase-0 signature contradicts D-6-rev2; SeisSol cross-check
ill-defined under λ; exchange version labels clash; D-1 experiment missing
from Phase 5; the Phase-1 addendum is orphaned.** (Rev-2 amendments recorded as
decisions but never propagated into phase contracts.) **Fix:** full
`LtsClusteringOptions` interface (rate, max_clusters, wiggle_scan, nc_cap=6,
merge_loss_tol=0.05, cell_cost) + result fields (lambda, modeled_cost);
SeisSol cross-check runs in RAW mode (λ=1, no cap) with both histograms
printed — acceptance applies to raw, go/no-go to production; Phase 4 split
into 4a (full-field exchange, parity gate) / 4b (EDGE 3-buffer, byte-compared
vs 4a; **Phase-5 performance measured on 4b**); D-1 relaxation added as
Phase-5 Req 5 (gated, non-blocking, `lts_fault_maxdiff=1`); addendum folded
into Phase 1 with the concrete METIS spec (ncon=num_clusters unit-weight
constraints, ubvec 1.05, union-find fault merge composed first, fallback
scalar weights).

**[P-010] Q/Q_new double-buffer dies under LTS.** One vector holds elements at
different time levels; per-cluster correct must update only its elements
in place. **Fix:** normative single-`Q` in-place contract +
`AdvanceADERCluster` signature (cluster ref, dt_step, order, Q in-place,
I_cluster, accumulate buffers consumed+zeroed, D(k) store read-only); order:
volume+faces+buffer-add → per-element M⁻¹ (exact per cluster) → `Q +=`.

**[P-011] D-3 wires the wrong nucleation function and stomps other clusters.**
The acceptance configs use `gradual_overstress_compact_circular`, whose
absolute form is a DIFFERENT function; and the absolute appliers write the
whole fault vector. **Fix:** route ALL kinds through their absolute forms via
INucleationMethod (gradual → ApplyGradualOverstressAbsolute; compact-circular
→ ApplyGradualOverstressCompactCircularAbsolute; instantaneous already
one-shot); add range overloads `(…, qp_begin, qp_end)`; each cluster applies
its own range at its own stage times, pinned to the GTS sub-step convention;
telescoping/partition-independence unit test.

**[P-012] Conservation harness physically ill-defined.** "periodic/absorbing
box" — conservation does not hold with absorbing boundaries; energy is never
conserved under upwind flux. **Fix:** periodic Cartesian tet box
(`Mesh::MakePeriodic`): ∫ρv_i (3) and ∫σ_ij (6) are exact invariants (interior
upwind fluxes telescope), drift <1e-12×initial-norm per sync; energy monotone
decay only. Fallback: traction-free box, ∫ρv_i only.

### MODERATE

**[P-013]** Missing signatures/layouts: per-cluster predictor (writes into
full-size vectors zeroing ONLY the cluster's dof blocks — the existing
routine's whole-vector zeroing must not be reused), D(k) store (RAW unscaled
coefficients, layout `[slot][k][comp][i]`, provider_slot_of_elem built in
Phase 1, epoch counter), `IntegrateTaylor(a,b,stack)→out[NUM_STATE][ndof]`
(a,b relative to the provider's last PREDICT time — the expansion point; rename
index `c`→`comp`), `LtsTick{predict_clusters, correct_clusters, dt_actual}`,
per-cluster friction `Advance(range…)` semantics (global indexing, base
pointers + (begin,end), I_imp outside the range untouched,
`SetSubStepFaultImposedStates(range,…)`, per-range ImposedGuard).
**[P-014]** "optionally RETAIN D(k)" — retention is MANDATORY for provider
elements when lts≠off (compiled-but-untaken on GTS ⇒ byte-exact).
**[P-015]** Provider/consumer ELEMENT sets + slot maps assigned to Phase 1.
**[P-016]** LTS × levers: lts≠off implies fused (shared-CK) predictor
semantics regardless of the flag (logged); deriv-cache/face-cache compose;
2×2 lever smoke to 1e-15.
**[P-017]** Post-flip guard contradiction: `lts≠off` IMPLIES the fault-locality
partition automatically (flag becomes a no-op alias); vacuous on fault-free
meshes.
**[P-018]** Checkpoint V2: 64-bit FNV-1a over (rate, num_clusters, cluster-id
sequence in serial-mesh order, IEEE bits of dt_base and λ); GTS keeps writing
V1 byte-identically; V2+lts=off refused in v1; refusal-path unit tests.
**[P-019]** R-101 rekey: run at every sync with t_sync ≤ T_nuc of the ACTIVE
kind (drop %100); under D-2 it is vacuous — repurpose as a locality tripwire
(assert zero shared fault faces) with the np=2==np=1 per-sync fault-state
checksum as the real cross-rank gate.
**[P-020]** GAP registry + the "five rank-global fault structures" enumerated
(dof_data; fault_coords; Q_pointwise±; I_imp±_flat; the deltaT/weights/
tau_nodes schedule); "GAP: ParaView trigger" resolved (per-sync max of
slip_rate_substep_max feeds the regime detector; granularity coarsens to
sync — accepted).
**[P-021]** Unit-test matrix consolidated as a normative appendix (scheduler
goldens; ragged-final-cycle; nucleation telescoping; checkpoint refusals;
conservation harness; λ-scan/auto-merge; bin edges; truncated-tfinal byte
gate; mixed-neighbor 3-cluster; ghost-poison np=2 seam; 2×2 levers; stations
byte-identical reorder gate).

### LOW
**[P-022]** Byte gates pin `lts_wiggle="off"`; np>1 single-cluster gate
compares GTS *given the LTS partition* (legal degree of freedom — say so).
**[P-023]** Step-keyed consumer inventory (V_max print, receivers, R-101
nonfatal print, downstream log parsers) listed with new sync-keyed cadences;
NaN detection latency note + optional collective-free per-tick local isfinite.
**[P-024]** Wording: "last correction time" → "last predict time (= the
expansion point)"; header claim "Phases 1/4 amended" corrected; D(k) retention
copies out during the recursion (the ping-pong scratch is shared across
cluster invocations — no lazy aliasing).

## Summary
- Critical: 12 (P-001…P-012) — concentrated in the scheduler/coupling core and
  the rev-2 propagation gaps; both reviewers converged independently on P-001.
- Moderate: 9 (P-013…P-021) · Low: 3 (P-022…P-024)
- Plan compliance (rev-2 self-consistency): PARTIAL — amendments recorded as
  decisions but not propagated into phase contracts.
- Verdict: **FAIL as rev 2 — must fix before implementation.** All findings
  are addressed in **rev 3** (same file), which adds a normative Interfaces
  appendix and a Unit-Test Matrix appendix.

## Unreviewed Areas
- The Phase-7 p-drop stub (deliberately a pointer to its own future plan).
- Phase-6 flip mechanics beyond the guard interaction (P-017).
