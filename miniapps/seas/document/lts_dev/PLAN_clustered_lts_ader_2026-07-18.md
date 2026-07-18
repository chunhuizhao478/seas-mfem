# Implementation Plan: Clustered Local Time Stepping (LTS) for the ADER path of seas_spatial_dyn_driver

**Date:** 2026-07-18 (rev 2, same day) · **Status:** PROPOSED (grounded, not implemented)
**Grounding:** two multi-agent investigations (2026-07-18) over this repo and the
local SeisSol source (v1.3.1-2135-gdc6db6513), adversarially verified; the
quantitative motivation is `../code_optimization_dev/ANALYSIS_mfem_vs_seissol_speed_2026-07-18.md`.
**Rev 2:** the method choice was stress-tested against ALL alternative LTS/multirate
families (elementwise ADER, leapfrog/Newmark-LTS, AB-multirate, MRI-GARK,
locally-implicit/IMEX, tent-pitching, p-adaptivity) by a 9-agent literature +
repo-fit + judge-panel study — see `ANALYSIS_lts_method_selection_2026-07-18.md`
in this folder. Verdict (3 independent judges, unanimous): clustered rate-2
ADER-LTS wins, UPGRADED with the EDGE-2022 package (λ-wiggle, Nc-cap+auto-merge,
flux-premultiplied 3-buffer exchange) and a staged far-field p-drop that SeisSol
structurally cannot copy. Decisions D-1/D-6 and Phases 1/4 amended accordingly;
Phase 7 added.

## Summary — read this first

**The problem.** Our dynamic-rupture driver advances every element of the mesh
with one shared time step, set by the single worst element. On the SAFS regional
meshes the allowed step sizes of different elements span a factor of about one
thousand, so more than 99% of all element updates are wasted work. Measured on
the coarse SAF-ALT benchmark: SeisSol finishes 150 simulated seconds in under
five hours while our driver would need about seven days — and our code is
actually *faster per element update* than SeisSol. The entire gap is scheduling.

**The fix.** Group elements into clusters by their allowed time step (each
cluster's step is a power of two times the smallest), and advance each cluster
at its own rate. Neighboring clusters stay consistent because the ADER
predictor already produces, for every element, a small polynomial describing
its solution *over the whole time step* — a coarse element can therefore hand a
fine neighbor exactly the time-integrated information the fine steps need. This
is precisely SeisSol's mechanism, and our predictor is mathematically the same
object, so the design transfers rather than being invented.

**Expected outcome.** On the coarse SAF-ALT mesh, LTS removes 38.7× of the bulk
element updates (22× of the fault updates). Combined with the already-shipped
SeisSol-equivalent CFL setting, the idealized ceiling is ~100 simulated-seconds
per hour versus 2.6 today — i.e. from ~2.4 days down to ~1.5–3 hours for the
150 s benchmark, overtaking SeisSol itself. The acceptance gate claims only
≥15× end-to-end (half the ideal), to leave honest room for cluster-management
overhead and load imbalance.

**Main tradeoff / biggest risk.** LTS produces a *different* (equally valid)
numerical trajectory than global stepping — bit-for-bit comparisons with
existing results become tolerance comparisons. The riskiest machinery is the
cluster-boundary bookkeeping: a fine cluster must integrate its coarse
neighbor's polynomial over exactly the right sub-interval, and the fault's
friction solver — today a single rank-wide sweep with one dt — must be split
into per-cluster sweeps. Getting an interval or a buffer wrong produces silent
non-conservation, not a crash.

**What this does NOT do.** It does not speed up a single element's update
(kernel efficiency is a separate ~6× axis), does not apply to the RK/mixed-flux
path (ADER-only by guard), and does not change any result of the default
global-stepping mode — LTS lands strictly opt-in first, and "default" means
default for this driver's SAFS production configs after validation, with the
TPV gold decks explicitly pinned to global stepping.

## How to read this plan
- The **Summary** above is the whole idea.
- Skim the **"In one sentence"** line of each phase to see the arc.
- **Detailed Requirements** are the implementation agent's contract.
- The **Glossary** defines every shorthand used below.

## Glossary

| Label | Plain-language meaning |
| ----- | ---------------------- |
| GTS | Global time stepping — every element uses the same dt (today's behavior). |
| LTS | Local time stepping — each cluster of elements uses its own dt. |
| rate-2 clustering | Cluster c steps with dt_min·2^c; an element joins the largest c with 2^c·dt_min ≤ its own CFL dt (lower bin edge ⇒ never violates its own CFL). |
| maxdiff rule | Face-neighboring elements' clusters may differ by at most 1; across fault faces by exactly 0 (both sides same cluster). Enforced by an iterative clamp to a fixed point. |
| D(k) stack | The per-element Cauchy–Kovalevskaya Taylor coefficients the ADER predictor computes; Q(τ)=Σ τ^k/k!·D(k). Today discarded after use; LTS retains them for elements that border a finer cluster. |
| accumulate buffer | Per-element storage where a fine cluster sums its sub-step time-integrals so a coarse neighbor can consume one integral over its whole step. |
| provider / consumer roles | Per face: the coarser side provides its D(k) stack; the finer side integrates it over each sub-interval and fills the accumulate buffer the coarse side later consumes. Equal clusters exchange plain time-integrals (GTS relation). |
| sync point | A global time all clusters land on exactly (by truncating their last steps); output, checkpoints, and collectives happen only here. |
| tick | One step of the finest cluster within a sync interval; the schedule is expressed in ticks. |
| the consume path | The verified shared-fault substep mechanism (SEAS_DIAG_SHARED_SUBSTEP_CONSUME) that replaces the partition-dependent R-1601 inline fallback. |
| fault-locality partition | Existing opt-in `--partition-fault-locality`: both elements of every fault face live on one rank. |
| Phase-0 report | The cluster-histogram/speedup-predictor tool built first; its output gates the whole project. |

## Technical Overview

The ADER macro step (`drivers/spatial_dyn_driver.cpp:399-515`) is predictor →
fault-friction substeps → corrector. The predictor
(`wave_operator.inl:1597-1688`) is an element-local CK recursion producing
exactly the Taylor object SeisSol's LTS couples with; it currently discards the
coefficients (ping-pong scratch, `:1644-1686`) and the corrector consumes one
rank-global time-integral `I`. LTS restructures this into per-cluster
invocations driven by a tick schedule, adds retained D(k) stacks +
accumulate buffers on cluster-boundary faces (SeisSol `LtsSetup.cpp:96-185`
roles, including the fifth rule at `:156-162`), integrates coarse-neighbor
Taylor expansions over sub-intervals (SeisSol `TimeBasis.h:82-93`,
`TimeCluster.cpp:822-969`), and splits the fault iterator
(`friction_substep_iterator.*`) into per-cluster sweeps over
cluster-contiguous fault-face blocks. Clustering runs on the serial mesh
before partitioning (`Mesh::ElementToElementTable`, fault faces via
`FindFaultFaceIndices`), and partitioning gets LTS cost weights
(cost·2^(maxC−c), SeisSol `WeightsModels.cpp:18-35`) through the existing
explicit-partition injection point (`spatial_dyn_driver.cpp:1268-1293`).

## Constraints

- **Byte-exact default:** every touched shared file (`wave_operator.{hpp,inl}`,
  `fault_face_flux.*`, `friction_substep_iterator.*`) recompiles into the TPV
  gold-deck binaries (Makefile `:2297-2343`); with LTS off, behavior must be
  byte-identical (established repo pattern: `use_shared_ck`, `resample`).
- **ADER-only:** mixed flux requires RK; RK is a global MOL stepper. Guard:
  `--lts`/`[numerics].lts` rejected unless `time_integrator=ader` and
  `mixed_flux=none`.
- **Deterministic clustering:** cluster ids must be a pure function of the
  serial mesh + material + config (global element ids), never rank-local data —
  the R-101 cross-rank tripwire aborts otherwise.
- **Matched collectives:** the R-1600 contract (every rank with shared faces
  calls `ExchangeFaceNbrData` the same number of times) must hold ⇒ all ranks
  execute the same global tick schedule, with possibly-empty local cluster sets
  (the SeisSol model).
- **cluster-0 dt definition:** binning consumes the FULL driver dt product
  `dt_e = cfl · CflSafetyFactor(cfg) · h_e / cp_e` (the order/safety factor
  lives driver-side, `spatial_friction.hpp:731-736` — not in `ComputeMaxDt`).
- **Stability invariant:** every element steps at dt_cluster ≤ its own dt_e
  (lower bin edge); keep SeisSol's maxdiff ≤1 rule (required by the
  buffer/actor machinery; also the literature-validated regime — Dumbser &
  Käser LTS-ADER-DG with dissipative upwind flux). The pure-upwind
  (dissipative) flux is a precondition: central/mixed flux stays excluded.
- **Fault-locality prerequisite (decision D-2 below):** LTS v1 requires
  `--partition-fault-locality`, making every fault face rank-interior and
  eliminating the shared-fault-QP × cluster hazard class (R-1601/R-101)
  outright.
- **No LTS edits leak** outside `dynamic/` + the spatial driver + new files
  (memory rule [C2]: bp5/bp1/bp2/domain/fault/solver and
  `friction/dieterich_ruina.hpp` are no-touch).

## Design decisions (need sign-off; defaults chosen)

| ID | Decision | Chosen default | Alternative rejected because |
|---|---|---|---|
| D-1 | Fault faces under LTS | **Per-cluster fault machinery** (SeisSol-style: each fault face lives at its own — post-clamp — cluster rate). DR faces force BOTH sides same-cluster in v1 (the published Uphoff SC'17 rule); relaxing to maxdiff≤1 ACROSS the fault is a named Phase-5 experiment (beyond published work; potential publication). | Forcing all fault faces to the minimum fault cluster costs 5.8×10¹¹ face-updates on the ALT mesh — 4.7× MORE than the entire LTS element budget. Not viable. |
| D-2 | Shared (cross-rank) fault faces | **Require `--partition-fault-locality` for LTS v1** | Extending the consume path per-cluster across ranks couples LTS to the R-1601 promotion; sequencing both at once doubles risk. v2 may lift this. |
| D-3 | Nucleation under LTS | **Switch the LTS path to the idempotent absolute form** `ApplyGradualOverstressAbsolute` (SET τ_nuc=S(t)·amp; RK-proven, time-partition-independent) | Per-cluster incremental telescoping duplicates state and invites drift; the absolute form is gated LTS-only so GTS stays byte-identical. |
| D-4 | Scheduler | **Deterministic recursive tick schedule** (fine-to-coarse within each tick), not SeisSol's async actor model, in v1 | Actors exist to overlap MPI; v1 buys correctness first. The tick loop preserves matched collectives trivially. Actor/overlap is a v2 optimization. |
| D-5 | Default flip | **Two-stage:** (i) land opt-in (`[numerics].lts = "off"|"rate2"`, default off); (ii) after Phase-5 acceptance, set `lts="rate2"` in the SAFS production/speed configs AND flip the parser default, simultaneously pinning `lts="off"` in every TPV-spatial config + re-goldening the TPV104-spatial smoke | A one-shot default flip silently changes TPV104-spatial production trajectories (they run this driver). |
| D-6 | Wiggle factor / auto-merge | **REVERSED (rev 2): λ-wiggle grid search + Nc-cap (≈5–6) with cost-model auto-merge are IN v1** (Breuer & Heinecke IPDPS 2022; EDGE realized 94–95% of theoretical LTS speedup with them, +17.5% from λ alone) | They are not polish: the Nc cap is the published fix for GPU-LTS collapse (SeisSol GPU: ~1.3× without it) and directly serves our element-local forall GPU design; both are preprocessing-time features, cheap to carry from Phase 0 onward. |
| D-7 | Exchange payload (rev 2) | **EDGE-style fixed 3-buffer, flux-premultiplied exchange** as the Phase-4 target (send flux-projected payloads, static per-level schedule); v0 stepping stone = full ghost-field exchange per due-tick | EDGE beat SeisSol's own LTS comm by 1.26–1.48× with this; the static schedule is exactly our matched-collective contract. |

## Phase 0: Cluster report & go/no-go

**In one sentence:** Before writing any stepping code, measure — for each target
mesh — how many elements fall in each would-be cluster and what speedup that
predicts, and verify our histogram against SeisSol's own log on the shared mesh.

### Goal
Close GAP-X2: the payoff currently rests on one mesh's SeisSol histogram.
Produce the histogram + predicted speedup from OUR dt formula for: coarse ALT,
v4_0_0 ALT/PREF p3, and the fault-band graded meshes.

### Files to Create
- `dynamic/lts_clustering.hpp` — pure functions (no MPI): per-element dt →
  cluster ids → maxdiff fixpoint → histogram/speedup stats (reused by all
  later phases).

### Files to Modify
- `drivers/spatial_dyn_driver.cpp` — add `--lts-report` (implies the existing
  `--dry-run` exit): after operator construction, compute per-element
  `dt_e = cfl·CflSafetyFactor(cfg)·h_e/cp_e` (matrix path: via
  `GetPerElementCflLength()`/`GetPerElementMaterial()`; scalar path: promote
  the ctor loop temporaries `wave_operator.inl:65-104` to stored vectors —
  byte-exact-neutral), run clustering, print the SeisSol-style histogram
  (cells + fault faces per cluster), element-update speedup (harmonic form),
  and assert `dt_cluster(e) ≤ dt_e` for every element (GAP-A1).
- `dynamic/bimaterial_wave_operator.hpp` — re-label the two accessors from
  "Test-only" to production-sanctioned (one-line comment change).

### Detailed Requirements
1. `struct LtsClustering { std::vector<int> cluster; int num_clusters; double dt_base; };`
   built by `BuildLtsClustering(const std::vector<double>& dt_e, const mfem::Table& elem_to_elem, const std::vector<std::pair<int,int>>& fault_face_elem_pairs, int rate /*=2*/, int max_clusters /*=32*/)`.
   Binning: `c = floor(log2(dt_e/dt_min))` clamped to [0,max_clusters-1] with
   the SeisSol open-interval rule (element joins cluster c iff 2^c·dt_min ≤
   dt_e < 2^{c+1}·dt_min). dt_min is the global minimum of dt_e.
2. maxdiff fixpoint exactly as SeisSol (`LtsWeights.cpp:548-653`): iterate
   `cluster[i] = min(cluster[i], cluster[j] + diff)` over all element-adjacency
   pairs (diff=1) and all fault pairs (diff=0) until no change. Serial mesh
   only (runs before partitioning); deterministic by construction.
3. Report format mirrors SeisSol's (`ClusterLayout.cpp:71-119`): per-cluster
   (cells, fault faces), dt per cluster, elementwise + clustered speedup in
   BOTH SeisSol's arithmetic-mean form (for cross-checking their log) and the
   harmonic element-update form (the honest number).
4. Cross-check: on `safalt_dl_safgh_2M_walls_mmg.msh` + the CVM sidecar the
   histogram must reproduce SeisSol's (123, 822, 3346, 13917, 171472, 487273,
   192506, 87246, 312888, 49565, 136) to within the material-sampling caveat
   (document any deviation >1% per cluster).

### Acceptance Criteria
- [ ] Histogram matches SeisSol's on the shared mesh (per-cluster deviation ≤1% or explained).
- [ ] `dt_cluster(e) ≤ dt_e` assertion passes on all target meshes.
- [ ] Predicted harmonic speedup ≥10× on at least the coarse ALT and one production mesh (**go/no-go gate for the whole project**).
- [ ] `make test` unchanged (report code is dry-run-only; scalar-path storage addition is byte-exact-neutral, verified by `test-ader-tpv102-smoke`).

### Dependencies
Depends on: nothing. Required by: all later phases.
**Estimate:** 2–4 days.

## Phase 1: Clustering wired into the run path (still GTS)

**In one sentence:** The driver computes and carries the cluster layout (and the
LTS-aware partition) on every run, while still stepping globally — so all
bookkeeping can be validated with zero numerical change.

### Goal
Land `[numerics].lts = "off"|"rate2"` (+ `--lts` CLI, default off), the
serial-mesh clustering call, cluster-contiguous orderings, and the LTS-weighted
partition — with stepping unchanged.

### Files to Create
- `dynamic/lts_layout.hpp/.cpp` — run-side layout: per-rank per-cluster element
  index lists, per-cluster interior/boundary/fault face lists, cluster-boundary
  face roles (provider/consumer/GTS + the fifth rule: accumulate-buffer cell
  with any equal-cluster neighbor also provides derivatives,
  SeisSol `LtsSetup.cpp:156-162`), fault-face → cluster map, and the
  tick schedule table for one sync interval.

### Files to Modify
- `spatial/code/spatial_friction.{hpp,cpp}` — parse `[numerics].lts`
  (string, default "off"), `lts_max_clusters` (default 32), `lts_sync_dt`
  (default "auto" = coarsest cluster dt). Guards: reject lts≠off with rk/mixed
  flux; reject without `--partition-fault-locality` (D-2).
- `drivers/spatial_dyn_driver.cpp` — run clustering on the serial mesh before
  `ParMesh` construction; feed METIS an explicit partition with vertex weights
  `cellCost·2^(maxC−c)` (SeisSol ExponentialWeights, `WeightsModels.cpp:18-35`)
  through the existing explicit-partition path (`:1268-1293`), composed with
  the fault-locality union-find merge; attach cluster ids to local elements
  after partitioning (global-id keyed).
- `dynamic/fault_locality_partition.hpp` — accept optional vertex weights.

### Detailed Requirements
1. Cluster ids computed ONCE on the serial mesh (rank 0), broadcast, and
   mapped to local elements via the partition array — never recomputed from
   rank-local data (R-101 constraint).
2. Fault-face cluster = the (identical, post-clamp) cluster of its two
   elements; abort if they differ (mirror SeisSol `MeshLayout.cpp:185-189`).
3. `dof_data` fault-QP ordering becomes cluster-contiguous: fault faces sorted
   by (cluster, global face id), QP blocks per face unchanged (the resample
   contiguity contract `friction_substep_iterator.cpp:383-386` is preserved
   per face). A permutation table maps old→new for checkpoint compatibility.
4. Log line (rank 0): cluster histogram + partition imbalance per cluster.
5. With `lts="off"`: none of the above executes — byte-identical.
   With `lts="rate2"`: layout is BUILT and logged but stepping is still GTS
   (the loop ignores it) — this phase changes trajectories ONLY through the
   partition, which is already a legal degree of freedom.

### Acceptance Criteria
- [ ] `lts="off"` byte-identical on `test-ader-tpv102-smoke` + a SAFS smoke.
- [ ] `lts="rate2"` at np∈{1,2,10}: identical cluster histogram, R-101 green, run completes under GTS stepping.
- [ ] New unit test `test_lts_clustering`: binning edge cases (dt exactly at a bin edge joins the LOWER cluster; single-element mesh; all-equal dt ⇒ 1 cluster), maxdiff fixpoint convergence, fault diff=0.

### Dependencies
Depends on: Phase 0. Required by: Phases 2–6.
**Estimate:** ~1 week.

## Phase 2: Multi-cluster stepping, bulk only (np=1)

**In one sentence:** Elements actually advance at their cluster's rate — for a
fault-free (or fault-frozen) problem on one rank — with the cluster-boundary
coupling done by integrating the coarse neighbor's Taylor polynomial.

### Goal
The core numerical machinery: retained D(k), accumulate buffers, per-cluster
predictor/corrector, sub-interval integration, sync-by-truncation.

### Files to Create
- `dynamic/lts_stepper.hpp/.cpp` — the tick loop:
  ```
  for tick in 0 .. ticks_per_sync-1:
      for c in due_clusters(tick), FINE→COARSE:  predict(c)   # CK; store D(k) for provider cells; buffer I_c
      for c in due_clusters(tick), FINE→COARSE:  correct(c)   # faces + volume + M^-1 + add
  ```
  `due_clusters(tick) = { c : tick % 2^c == 0 }`. Last steps truncate to the
  sync time (`dt_step = min(dt_c, t_sync − t_c)`, SeisSol `ActorState.cpp:75-77`).
- `dynamic/lts_time_basis.hpp` — `IntegrateTaylor(a, b, D(k) stack) →
  I[c][dof]` with `coeff[k] = (b^{k+1}−a^{k+1})/(k+1)!` (SeisSol
  `TimeBasis.h:82-93`); unit-tested against analytic polynomials.

### Files to Modify
- `dynamic/wave_operator.{hpp,inl}` + `bimaterial_wave_operator.inl`:
  1. Predictor variant `ComputeADERSubStepStatesAndIntegralCluster(cluster_elems, dt_c, …)`
     iterating an element index list instead of `0..ne-1` (the kernels already
     use per-element offsets — mechanical); optionally RETAIN D(k) into
     `lts_dk_store_` for provider elements only (size: order×9×ndof_per_el ×
     n_provider_elements; sized and logged against the deriv-cache budget).
     Lifetime: a provider's D(k) is written at its cluster's predict and must
     survive until its NEXT predict (i.e. across all 2^Δ fine sub-steps) —
     enforced by an epoch counter assert (GAP-A3).
  2. Corrector variant `AdvanceADERCluster(cluster, dt_c, …)`: volume + face +
     M^-1 restricted to the cluster's elements (`ApplyMassInverse` is
     element-block-diagonal — verified `wave_operator.inl:5995-6021` — so
     per-cluster application is exact).
  3. Cluster-boundary faces (the accumulate-buffer graft, GAP-B1): the FINE
     side visits the face every fine step; its own flux uses
     `IntegrateTaylor(t_rel, t_rel+dt_fine, D_coarse)` for the neighbor state;
     it scatters its own side immediately AND adds the anti-symmetric coarse-side
     contribution into that element's accumulate buffer; the COARSE side, at
     its correction, adds its accumulate buffer INSTEAD of visiting the face.
     Conservation is preserved because both sides' totals come from the same
     per-sub-interval flux evaluations (single-flux-evaluation invariant).
  4. Full dt audit (GAP-B3): every `dt` read inside
     `ComputeADERFaceFluxRHS`/boundary branches (`I/dt` at `:4319-4325,4511`,
     `Q_imp·dt` at `:4741`, PML inline `:5907`) takes the per-cluster dt of the
     OWNING element's cluster; each site gets a code comment naming its cluster.
- `drivers/spatial_dyn_driver.cpp` — sync-interval outer loop replacing the
  nsteps loop when LTS is on: `while (t < tfinal) { advance_to(t_sync); … }`;
  outputs/V_max/NaN-tripwire/checkpoint move to sync points; `step` becomes the
  sync counter (prints labeled "sync").
- `io/tpv104_checkpoint.hpp` — **V2 schema**: new magic tag; stores
  (t, sync_step, lts mode, cluster-layout hash = hash of the serial-mesh
  cluster vector, partition hash); READER REFUSES a V1 file when LTS is on and
  refuses hash mismatches (GAP-D1/D2). V1 continues to work for GTS.

### Edge Cases to Handle
- One cluster total ⇒ the tick loop degenerates to GTS: **byte-identical**
  gate vs `lts="off"` (same partition), the strongest cheap correctness check.
- Truncated final steps: a fine cluster may need a truncated step while its
  coarse neighbor is mid-step — the sub-interval integration handles it since
  `b ≤ dt_coarse` always (tau_nodes never exceed the coarse dt; consumers
  integrate the stack rather than re-calling the predictor with shifted nodes).
- Empty clusters on a rank (after partitioning): predict/correct no-op but any
  collective in the phase still executes (matched-collective rule).

### Acceptance Criteria
- [ ] `test_lts_time_basis`: sub-interval Taylor integrals exact for polynomials up to order 4; sum of sub-interval integrals == whole-interval integral to 1e-15.
- [ ] Single-cluster LTS == GTS **byte-identical** (np=1, TPV102-style box).
- [ ] Multi-cluster LTS vs GTS on a smooth wave problem: L2 difference at t_final consistent with truncation order (convergence study at 2 resolutions), and conservation: global momentum/stress integrals drift <1e-12 per sync interval on a periodic/absorbing box.
- [ ] `lts="off"` still byte-identical everywhere (`test-ader-tpv102-smoke`).

### Dependencies
Depends on: Phase 1. Required by: Phase 3.
**Estimate:** 2–3 weeks (the core of the project).

## Phase 3: Fault (dynamic rupture) under LTS (np=1)

**In one sentence:** Fault faces advance at their own cluster's rate — the
friction iterator runs once per fault-cluster step over that cluster's
contiguous block of fault points, and nucleation switches to the
time-partition-independent absolute form.

### Goal
Split the five rank-global fault structures per cluster (GAP-C2) and validate
rupture physics LTS-vs-GTS.

### Files to Modify
- `dynamic/friction_substep_iterator.{hpp,cpp}` — `Advance(range, dt_c, …)`
  over a (begin,end) QP range (cluster-contiguous by Phase 1); per-cluster
  deltaT sequences (dt_c/O each); the Σ deltaT == dt check per invocation.
- `drivers/spatial_dyn_driver.cpp` — the fault half of predict/correct per
  fault-bearing cluster: per-cluster tau_nodes on [0,dt_c]; per-cluster
  `EvaluateBulkAtFaultQPsCanonical` restricted to the cluster's fault faces
  (both adjacent elements are same-cluster by the diff=0 clamp, so the bulk
  trace is time-consistent by construction — SeisSol's same insight,
  `LtsSetup.cpp:115-129`); per-cluster
  `SetSubStepFaultImposedStates(range)`; `slip_rate_substep_max` reduced per
  sync for output regime detection (GAP: ParaView trigger).
- `dynamic/spatial_nucleation.cpp` wiring — LTS path calls
  `ApplyGradualOverstressAbsolute(t_cluster)` (D-3); GTS path untouched.
- R-101 tripwire + DIAG gates rekeyed to sync counter (GAP-C4).

### Acceptance Criteria
- [ ] Single-cluster LTS == GTS byte-identical WITH fault (TPV104-spatial smoke config).
- [ ] Multi-cluster LTS vs GTS on TPV104-spatial 200m: rupture arrival times at the standard stations within 1%; final slip within 1%; no spurious V_max transients at cluster boundaries crossing the fault's cluster edge.
- [ ] `test_friction_substep_iterator` extended: two clusters, ranges advance with different dt, per-QP results equal a reference where each QP is advanced standalone with its own dt.
- [ ] Nucleation: LTS run reproduces the GTS breakout time on the SAFS coarse smoke within 2% (absolute-form gate).

### Dependencies
Depends on: Phase 2. Required by: Phase 4.
**Estimate:** 1–2 weeks.

## Phase 4: MPI

**In one sentence:** Every rank executes the same tick schedule (empty work
allowed), ghost data carries the right time-level per tick, and parity with the
single-rank LTS result is the gate.

### Goal
Parallel LTS preserving the matched-collective contract with fault faces kept
rank-interior (D-2).

### Detailed Requirements
1. All ranks run the identical tick loop; each tick's predictor ghost exchange
   executes on every rank (possibly empty payload) — R-1600 preserved by
   construction. v1 exchanges the full ghost field per due-tick (correct but
   unoptimized); v2 may filter by cluster.
2. Cluster-boundary faces at rank seams (bulk only — fault is rank-interior):
   the ghost side's D(k)/buffer travels with the exchange; `subTimeStart`
   bookkeeping per seam face mirrors SeisSol (`TimeCluster.cpp:731-846`):
   the fine side tracks the coarse neighbor's last correction time (a
   deterministic function of the tick index — no messages needed under the
   deterministic schedule). Interval-mismatch audit test required (GAP-B2).
3. V_max/NaN collectives per sync only. Per-sync collective count identical on
   all ranks by construction; assert with a debug counter.
4. Checkpoint V2 at sync points; restart np must equal write np (existing
   contract) AND cluster-layout hash must match.

### Acceptance Criteria
- [ ] np=2 LTS == np=1 LTS to 10 digits over ≥12 sync intervals (mirrors the consume-path validation methodology), LSW+RS, mixed-rank fault distribution.
- [ ] np=10 symmirror revalidation gate green; R-101 green through nucleation.
- [ ] No hang in a 30-min np=10 SAFS coarse smoke (the R-1600 failure mode is a 100%-CPU hang — watchdog the CI run).
- [ ] Restart mid-campaign (LTS, np=4) bit-continues (Q + dof_data + t).

### Dependencies
Depends on: Phase 3. Required by: Phase 5.
**Estimate:** 2–3 weeks.

## Phase 5: Performance validation on the benchmark

**In one sentence:** Measure end-to-end LTS speedup on the v4_0_0 coarse ALT
case and accept only if ≥15× over the GTS safety=1 baseline.

### Detailed Requirements
1. Rerun the p1 speed deck with `lts="rate2"`: target ≥15× vs the 2.60
   sim-s/hour GTS-safety-1 baseline (ideal 38.7×; the gate claims less to
   absorb overhead/imbalance). Record achieved vs Phase-0-predicted speedup.
2. Caliper: add a `lts.cluster` attribute; compare per-sync aggregates.
3. LTS-weighted partition A/B (Phase-1 weights on/off) — quantify imbalance.
4. p3 deck rerun: expect ~SeisSol-class throughput; document the residual
   kernel-efficiency gap as the next optimization axis (NOT this plan).

### Acceptance Criteria
- [ ] ≥15× end-to-end on p1 coarse ALT; fault-output fingerprint physically consistent with the GTS reference (rupture pattern, Mw within tolerance).
- [ ] A performance report doc in `document/lts_dev/`.

**Estimate:** ~1 week (mostly cluster time).

## Phase 1 addendum (rev 2)
Phase 1's clustering gains two requirements from the method study:
- the λ-wiggle grid search (λ∈(0.5,1], step 0.01, minimize the modeled update
  cost Σ cellCost/2^c/(λ·dt_min)) and the Nc-cap + auto-merge (lower the max
  cluster while cost ≤ (1+loss)·baseline) run inside `BuildLtsClustering`; the
  Phase-0 report prints the λ-scan curve and the chosen (λ, Nc);
- partition weights upgrade to **multi-constraint** (one METIS balance
  constraint per cluster level, Rietmann-style), with the scalar 2^(maxC−c)
  weight as fallback.

## Phase 6: Default flip (D-5)

**In one sentence:** After acceptance, LTS becomes the default for this driver
while every non-SAFS config that must keep its old trajectory pins `lts="off"`.

### Detailed Requirements
1. Set `lts = "rate2"` in the SAFS production + speed configs (explicit).
2. Flip the parser default `"off"` → `"rate2"` in the same commit that adds
   `lts = "off"` to every TPV-spatial config in-tree, and re-golden the
   TPV104-spatial smoke goldens under an explicit `lts="off"` pin.
3. Docs: seas CLAUDE.md section (LTS semantics of "step", checkpoint V2,
   fault-locality prerequisite); memory update.
4. The standalone tpv102/104/205 drivers never gain LTS wiring (guard stays).

### Acceptance Criteria
- [ ] A config with no `lts` key runs LTS; every gold/regression config runs identically to pre-LTS (pinned).
- [ ] `make test` green; full sbatch pre-flights updated (assert lts mode in the log like the cfl_dg_safety check).

**Estimate:** 2–3 days.

## Testing Strategy
- **Byte-exact ladder:** `lts="off"` byte-identical always; single-cluster LTS
  byte-identical to GTS; these two gates catch most wiring bugs for free.
- **Analytic:** Taylor sub-interval integration unit tests; conservation drift
  bounds on a box; convergence-order study for multi-cluster.
- **Physics tolerance:** TPV104-spatial stations (arrival <1%, slip <1%);
  SAFS breakout time <2%.
- **Parallel parity:** np=2==np=1 (10 digits), np=10 symmirror, R-101, hang
  watchdogs — reusing the shared-fault validation methodology already proven
  in this repo.
- **Performance:** Phase-0 predicted vs Phase-5 achieved speedup must agree
  within 2× or the difference must be explained (imbalance/overhead named).

## Risk Assessment

| Risk (plain language) | Severity | Mitigation |
|---|---|---|
| A wrong sub-interval or buffer produces silent non-conservation at cluster boundaries | HIGH | Single-flux-evaluation invariant; conservation drift gate; single-cluster==GTS byte gate; interval-mismatch audit test at seams |
| The fault machinery split (5 rank-global structures → per-cluster) corrupts rate-state history | HIGH | Cluster-contiguous reordering with permutation table + per-QP standalone-advance reference test; fault faces rank-interior (D-2) removes the cross-rank half of the problem |
| D(k) retention blows memory on big meshes | MED | Only provider cells store stacks (typically a thin shell between clusters); sized+logged vs the deriv-cache budget; abort with a named count if > threshold |
| LTS load imbalance eats the speedup (fine clusters concentrated on few ranks) | MED | LTS-weighted partition (Phase 1); Phase-5 A/B quantifies; SeisSol's exact weighting scheme is the template |
| Default flip silently changes TPV104-spatial science | MED | D-5 two-stage flip with in-tree pinning + re-goldening in the same commit |
| Checkpoint incompatibilities mid-campaign | MED | V2 magic tag + layout hash + refusal paths (GAP-D1/D2) |
| The deterministic tick schedule leaves MPI idle time SeisSol's actors would overlap | LOW (v1) | Accepted for v1 correctness; actor/overlap is the named v2 axis |
| Sync cadence erodes speedup (outputs force fine alignment) | LOW | Outputs already ≥1 s cadence vs coarsest dt ~42 ms; sync at max(coarsest dt, requested cadence) |

## Phase 7 (post-flip, rev 2): far-field p-drop — the SeisSol-impossible multiplier

**In one sentence:** Run the far-field clusters (≥6: 642k cells, zero fault
faces) at p1 while the fault region keeps p3, via a driver-level two-order-class
DOF layout — an axis SeisSol's compile-time fixed order structurally cannot copy.

Facts from the method study: MFEM's native variable-order FESpaces do NOT cover
our conforming tet ParMesh (nonconforming-mesh requirement; parallel hp is
quad/hex-scoped), but the driver owns its flat `[c·ndof_total + e·ndof_per_el + i]`
layout, so two order classes with prefix-sum offsets + per-class kernel batches
+ Dumbser's max-degree zero-padded interface flux rule are buildable in-driver
(multi-week, mechanical; ~40+ offset sites in `wave_operator.inl`). Ideal gain
~1.16× in update counts (clusters ≥6 are 14.7% of clustered cost) plus
memory-bandwidth relief; the friction machinery is untouched (fault region
stays p3). Gated by its own plan document when Phase 6 lands. Future
generalization: hp / damage-zone order boosting (CDBM).

**Total effort estimate: ~7–10 weeks** of focused work (Phases 0–5), plus the
flip; Phase 7 is a separately-planned follow-on. The go/no-go after Phase 0
costs only days and de-risks the rest. Realistic payoff (rev 2, from the method
study): ~35× realized element-update reduction for the backbone (94–95%
EDGE-demonstrated realization of the 38.73× ideal), ~45× with Phase 7 —
projected 1.2–1.5× faster than SeisSol o4 on the shared benchmark mesh.
