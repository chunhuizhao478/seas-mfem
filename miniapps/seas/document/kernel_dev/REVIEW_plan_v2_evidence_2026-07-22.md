# Code Review: PLAN_v2 performance program — evidence-support audit (2026-07-22)

## Review Scope

- **Plan:** `/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/document/kernel_dev/PLAN_v2_performance_program_2026-07-22.md` (163 lines, commit fb03488)

- **Files reviewed:**
  - Primary measured artifacts (the only admissible evidence base):
    - `/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/document/kernel_dev/B0_gate_reset_2026-07-22.md` (artifact **A**; Expanse job 52379131, 5 legs, FINAL)
    - `/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/document/kernel_dev/EVIDENCE_microbench_roofline_2026-07-22.md` (artifact **B**)
    - `/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/document/lts_dev/RESULTS_p5_tpv104_200m_2026-07-21.md` (artifact **C**)
    - `/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md` (artifact **D**, companion comm plan)
  - Prior review history: `document/kernel_dev/REVIEW_plan_ader_kernel_efficiency_2026-07-21.md` (**E**), `document/kernel_dev/REVIEW_p0_impl_2026-07-22.md` (**F**), `document/kernel_dev/PLAN_ader_kernel_efficiency_2026-07-21.md` (**G**, superseded v1)
  - Source: `dynamic/wave_operator.inl`, `dynamic/wave_operator.hpp`, `dynamic/lts_stepper.{hpp,cpp}`, `dynamic/lts_bulk_stepper.hpp`, `dynamic/lts_fault_stepper.hpp`, `dynamic/lts_partition.{hpp,cpp}`, `dynamic/lts_clustering.cpp`, `drivers/spatial_dyn_driver.cpp`, `tests/bench/bench_ader_kernel_variants.cpp`, `tests/bench/run_bench.sh`, `spatial/code/spatial_friction.hpp`, MFEM core `fem/pgridfunc.cpp`, `fem/pfespace.hpp`
  - Run decks: `jobs/kernel_dev/run_p0_baseline_expanse.sbatch`, `jobs/lts_phase5/tpv104_200m_lts_speed_expanse/{run_*.sbatch, tpv104_200m_lts_rate2.toml, tpv104_200m_lts_off.toml}`, `jobs/lts_phase5/tpv104_200m_seissol_expanse/{run_*.sbatch, parameters.par}`
  - Commits: `c547ac8` (2 Caliper scopes on LTS hot paths), `fb03488` (the v2 redesign)

- **Domain context consulted:** ADER-DG p3/O4 clustered-LTS structure (tick table, predict/correct ordering, seam corrector, cluster histogram at λ=0.63/Nc=6); MFEM `ParGridFunction::ExchangeFaceNbrData` semantics (blocking, 2× `MPI_Waitall` per round, buffers re-selected per call, `tag = 0` on the shared communicator); SeisSol v1.3.1 reference run configuration (16 ranks × 15 OMP, `ASYNC_MODE=SYNC`, matched dt/order); the project's recorded R-1600 matched-collective hang and the P-007 exchange audit; and the standing project rule that a previous LTS speedup projection was "far too wrong", so every claim and every directional decision must trace to a measurement.

- **Method:** Four evidence-ledger agents independently reconstructed artifacts A–D into FACT/GAP ledgers, recomputing every published number from its stated inputs and tagging each with its measurement context (leg, window, rank count, flag set, instrumentation). The plan was then audited through seven lenses — **arithmetic** (does each number reproduce), **traceability** (does each number resolve to A/B/C/D), **composition** (do the budgets and gates compose consistently), **direction** (is every directional decision measurement-backed), **code** (does the source support the plan's claims about it), **omission** (what load-bearing input is missing), **completeness** (what lever or risk was never considered). Every candidate finding was then attacked by three independent skeptics briefed to refute it on fact or on materiality; findings below survived ≥2/3, or were raised by the completeness critic and went unrefuted. Where a skeptic's dissent was substantive (severity downgrade, or a corrected fact), it is recorded in the Verification field rather than discarded.

## Findings

### [R-001] [CRITICAL] [PLAN_v2:12–22, :26] — The 1.31× face-cache gain is transferred GTS→LTS unflagged, and the only LTS measurement of it shows ~1.00×

**Category:** ASSUMPTION (lens = composition)

**Description:** `1538 = 2015 / 1.31`, where 2015 is the LTS-leg compute term (C §3, face-cache **OFF**) and 1.31 is the face-cache activation gain measured on the **GTS** leg (A Leg 2, 187.5 vs 245.8 µs·core). The plan explicitly flags the GTS→LTS transfer of the stage **split** (:76–82) but never flags this second, larger GTS→LTS transfer of the activation **factor**. The only face-cache-ON LTS datapoint in existence contradicts it: A Leg 5's implied rate is 1157/0.526 = 2199.6 s wall over 0.5 sim-s = 4399 s/sim-s against the face-cache-OFF LTS 4356 — a ~1.00× gain.

**Trigger:** The face cache does not help the LTS corrector (which uses `ComputeADERClusterSeamFaceFluxRHS`, not the cached GTS shared path).

**Actual behavior:** Budget 1538 + 2341 = 3879, comm = 60 %, kernel payoff 776, ratio 2.6×.

**Expected behavior:** Budget stated as a range (1538–2015 compute) with the transfer flagged, or the LTS face-cache-ON rate measured before the budget is quoted. If the transfer fails: budget 2015 + 2341 = 4356, comm 53.7 %, kernel figure 776 × 2015/1538 = 1017, ratio 2041/1017 = 2.0×.

**Suggested fix:** Insert after :14 — "*Transfer flagged:* the 1538 compute term is the LTS-leg 2015 s/sim-s (face-cache OFF) divided by the **GTS-leg-measured** 1.31× face-cache gain. The only face-cache-ON LTS datapoint we have — B0 Leg 5, implied 2199.6 s wall / 0.5 sim-s = **4399 s/sim-s vs 4356 face-cache-OFF, i.e. ~1.00×** — does not reproduce it. If the transfer fails, the budget is 2015+2341 = 4356 (comm 53.7 %), the kernel figure is 1017, and the payoff ratio falls from 2.6× to 2.0×. The re-measured LTS split must publish the face-cache-ON LTS compute rate."

**Evidence:** A "Composed end state" (1538 = 2015/1.31); A Leg 5 (Waitall Avg 1157 s = 52.6 % of wall, `--tfinal 0.5` per `jobs/kernel_dev/run_p0_baseline_expanse.sbatch:190-192` ⇒ wall 2199.6 s); C §Summary scoreboard (LTS 4356 s/sim-s, face-cache OFF).

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* the third skeptic conceded the arithmetic but argued materiality fails — Track B's ceiling (2015) sits below Track A's 2041 under every reading, no gate reads 1538/3879/60 %/776/2.6 (Track A gates on the measured 53 %, Track B on ratios vs a re-measured B0), and the plan already queues the LTS split for re-measurement at :80–81/:143–144. It also noted Leg 5's own wall is back-derived from one Caliper percentage over an assumed 0.5 s denominator with 4× harder init amortization, so it is a hypothesis worth a controlled leg, not a proof. **This finding is retained at CRITICAL only because R-002 supplies the source-level mechanism that converts it from "unflagged transfer" to "impossible transfer".**

---

### [R-002] [CRITICAL] [PLAN_v2:12-22] — `--face-cache` is structurally unreachable from the LTS corrector, so the headline budget 1538 + 2341 = 3879 ("16× behind") is wrong by construction

**Category:** BUG (lens = code)

**Description:** The plan's entire budget rests on LTS compute = 1538 s/sim-s, obtained as C's measured LTS compute 2015 (face-cache OFF) divided by the 1.31× face-cache gain measured on the GTS leg. The source shows this transfer is **impossible**, not merely unmeasured. `use_face_cache_` is read at exactly one place in the whole wave operator: `dynamic/wave_operator.inl:5475`, inside `ComputeADERFaceFluxRHS` (defined :5334). `ComputeADERFaceFluxRHS` has exactly one call site: `dynamic/wave_operator.inl:7288`, inside `AdvanceADER` — the **GTS** ADER corrector. The LTS corrector is `AdvanceADERClusterBulk` (:2190); its step-2 role-driven face sweep (:2229-2470) computes face geometry on the fly per quadrature point (`CalcOrtho` at :2286, `fe1->CalcShape` at :2297) and never consults `use_face_cache_` nor calls `ComputeADERFaceFluxRHS_CachedInterior_` (:909). `--face-cache` (driver:2179-2182 → `SetUseFaceCache`) therefore has **no effect on any LTS interior face**. Artifact A's own Leg 5 corroborates: LTS face-cache-ON implies 4399 s/sim-s vs the face-cache-OFF 4356 = 1.01×, and the plan does not reconcile this.

**Trigger:** Any reader taking 3879 s/sim-s, the 60/40 split, or "16× behind SeisSol" as the program's baseline; and any gate measured against a 1538 compute denominator.

**Actual behavior:** Plan states LTS compute = 1538 s/sim-s, wall 3879, comm share 60 %, 16× behind SeisSol, and calls the 1.31× "free".

**Expected behavior:** LTS compute stays at C's measured 2015 s/sim-s, wall 4356, comm share 53.7 %, **18.2× behind** SeisSol as-run 239. The kernel payoff recomputes to ~1016 s/sim-s (0.5623×2015×0.75 + 0.1296×2015×(1−1/1.5) + 0.0788×2015×0.5 = 850 + 87 + 79), not 776, so the Track-A:Track-B payoff ratio is **2.0×, not 2.6×**.

**Suggested fix:** Replace the Summary block (:12-27) with: "**Where we actually are.** On the benchmark problem our solver takes **4356 seconds of wall time per simulated second** (RESULTS_p5 §Summary, LTS leg, measured); SeisSol takes 239 — we are **18.2× behind**. The 1.31× face-cache activation measured on the GTS leg does NOT apply here: `use_face_cache_` is read only at `wave_operator.inl:5475` inside `ComputeADERFaceFluxRHS`, whose sole caller is `AdvanceADER` (:7288, GTS). The LTS corrector `AdvanceADERClusterBulk` (:2190) computes face geometry on the fly (:2286 `CalcOrtho`, :2297 `CalcShape`) and never consults the cache — consistent with B0 Leg 5's implied 4399 vs 4356 (1.01×)." Then set the table to: waiting 2341 (53.7 %), computing 2015 (46.3 %). Change :24 to "Fixing the waiting recovers **~2041 s/sim-s**; fixing every kernel we have evidence for recovers **~1016** — a factor of 2.0, not 2.6." Add a new Track-B phase **B0b**: "wire the interior-face cache into `AdvanceADERClusterBulk`'s role-driven sweep (currently GTS-only) — this is the actual free 1.31× on the LTS leg and must be measured, not assumed."

**Evidence:** `dynamic/wave_operator.inl:5475` (sole `use_face_cache_` read — verified by grep during this audit; the only other hits are the setter at :888/:905 and the guard at :871); `:7288` (sole `ComputeADERFaceFluxRHS` call, inside `AdvanceADER`); `:2286`/`:2297` (LTS sweep computes geometry on the fly); artifact A Leg 5 (implied 4399 vs C's 4356 = 1.01×); artifact C §3 (2015 + 2341 = 4356).

**Verification:** 3/3 skeptics failed to refute.

---

### [R-003] [CRITICAL] [PLAN_v2:84-98, :104-133 — the lever is absent from both] — `lts_wiggle="off"` is a zero-code lever worth ~1.59× on exchange rounds/sim-s, comparable to the entire Track B program, and neither plan mentions it

**Category:** DIRECTION (lens = completeness)

**Description:** The measured LTS leg ran at λ = 0.63 (C §1). λ is selected by `BuildLtsClustering`'s wiggle scan (`dynamic/lts_clustering.cpp:181-198`), whose objective is `compute_cost` (`:94-107`) = Σ_i w_i / (2^c_i · dt_base) — a pure **element-update count** per unit time. It contains **no sync-count, round-count, message or wait term of any kind**. The plan's own headline finding is that ~60 % of wall is communication, and the exchange-round rate is rounds/sim-s = (2·2^Nc − 3)/(2^(Nc−1)·λ·dt_min) ≈ 4/(λ·dt_min), i.e. inversely proportional to λ and near-independent of Nc. So the clustering was optimised against the 40 % term while inflating the 60 % term by 1/0.63 = **1.588×**. Arithmetic on the measured run: today 135.5 syncs/sim-s × 125 rounds = 16,938 rounds/sim-s; at λ=1, Nc=6 that is 85.35 × 125 = 10,669 (1.588× fewer); at λ=1, Nc=7 (SeisSol's clustering on the identical cross-validated histogram) 42.7 × 253 = 10,803 (1.568× fewer). Under the plan's **own** linear-in-round-count comm model — the model that produces its 2341→1200→600 payoff rows — that is 2341 × (1 − 1/1.588) = **867 s/sim-s recovered**. The compute cost of giving up the wiggle is the difference between MFEM's 3.68× ceiling at λ=0.63 and the 3.63× available at λ=1 (C §4, §Summary) = +1.4 % updates = 2015 → 2043, i.e. **+28 s/sim-s**. Net ≈ **−839 s/sim-s** for a one-line deck change: the decks never set `lts_wiggle`, so it defaults to `"scan"` (`spatial/code/spatial_friction.hpp:173`); `lts_wiggle = "off"` pins λ=1.0 (`drivers/spatial_dyn_driver.cpp:400`) and is already exercised by the P-022 byte gates. That is **108 % of Track B's entire projected 776 s/sim-s payoff**, with zero code, zero bitwise risk, and one Expanse leg. No lens in any prior review checked the clustering configuration itself; every lens audited only the numbers downstream of it.

**Trigger:** Bites immediately — every gate run, every composed measurement, and the production SAFS runs all inherit λ from the default `"scan"` and therefore pay 1.59× the exchange rounds. It also biases every Track A gate: A1/A2's round-merge payoff is measured on a baseline that is artificially round-rich, inflating the apparent win of code changes over a config change.

**Actual behavior:** The plan enumerates exactly two levers (comm code surgery 2041; kernel rewrites 776) and asserts "communication is the primary program". The clustering configuration that **sets** the communication rate is never examined; λ, wiggle, Nc and `lts_nc_cap` appear nowhere in the plan or its companion.

**Expected behavior:** A Phase-0 config sweep (λ ∈ {0.63, 0.8, 1.0} × Nc cap) measured on the LTS leg **before** either code track starts, with the clustering objective re-derived to include a communication term (cost = Σ updates + β · rounds/sim-s), since the measured budget says the comm term dominates 60/40.

**Suggested fix:** Insert as the new first item of §Sequencing, ahead of A1/B1: "**A0b — wiggle/cluster-cap sweep (config only, no code).** The measured leg ran λ=0.63, chosen by `lts_clustering.cpp:compute_cost`, an objective that counts element updates only and has no communication term. Exchange rounds/sim-s ≈ 4/(λ·dt_min), so λ=0.63 costs 1.588× the rounds of λ=1 while buying only 3.68/3.63 = 1.4 % fewer updates. Under this program's own linear-in-rounds comm model that is ~867 s/sim-s of exposed wait for ~28 s/sim-s of compute — a net ~839 s/sim-s, i.e. more than the whole of Track B, for the one-line deck change `[numerics] lts_wiggle = \"off\"`. Measure λ ∈ {0.63, 0.80, 1.00} on the LTS leg in a single job; re-baseline BOTH tracks' budgets on the winner before A1 or B1 starts. Then extend `compute_cost` with a comm term (β · rounds/sim-s, β calibrated from the same job) so the scan optimises the measured 60/40 budget rather than the 40 % half."

**Evidence:** `dynamic/lts_clustering.cpp:94-107` (`compute_cost` — updates only, no comm term); `:181-198` (scan minimises it); `drivers/spatial_dyn_driver.cpp:394-404` (`lts_wiggle = "off"` ⇒ `lambda_fixed = 1.0`); `spatial/code/spatial_friction.hpp:173` (default `"scan"`); `jobs/lts_phase5/tpv104_200m_lts_speed_expanse/tpv104_200m_lts_rate2.toml` (no `lts_wiggle` key ⇒ scan); C §1 (λ=0.63, 271 syncs, 125 rounds/sync, dt_min 3.66128e-4), §4 (3.68× ceiling), §Summary (3.63× at λ=1, histogram cross-validated identical to SeisSol's).

**Verification:** Completeness critic — unrefuted.

---

### [R-004] [CRITICAL] [PLAN_v2:39-42, :136 vs :34-37] — The plan permanently closes fault/friction, the single largest measured rank imbalance (94×), while naming that same imbalance as the cap on its primary track

**Category:** DIRECTION (lens = completeness)

**Description:** :39-40 closes fault/friction as "0.65 % of the step even through rupture — permanently out of scope". Three defects, and their interaction is the important one. (1) The 0.65 % is a Caliper **average across 256 ranks** with Min 2.7 s / Max 254 s — a **94× spread**, the largest rank imbalance recorded anywhere in artifacts A–D. The solver is bulk-synchronous, so the per-sync critical path is the **max** rank, not the mean: 254/8874 = 2.86 % of wall today, and at the plan's own 4× bulk speedup 254/(8874·0.9935/4 + 254) = **10.3 %** of the step. The closure decision uses the mean only. (2) "even through rupture" is not what was measured: A's own Harness gap #1 states the windows were "not separable", so 0.65 % is a whole-run blend in which the rupture window is 1317 of 4917 steps (27 %); **no rupture-window friction timing exists**. (3) The interaction nobody drew: :34-37 names the program's biggest risk as "the residual 24 % is load imbalance that no communication change touches, and it caps Track A". A 94× per-rank spread in a stage that only fault-adjacent ranks execute is the one measured mechanism in the evidence base that produces exactly that late-sender imbalance — and D's own Phase-0 escalation path is "if skew > ~40 %, the fault-weighted partition is promoted to co-requisite", i.e. the remedy is a fault-work-balancing change. So the plan closes, as permanently out of scope, **the only identified remedy for the cap on its primary track**. The 24 % that justified not escalating was itself measured pre-rupture (Leg 5, `--tfinal 0.5`, nucleation at t=1.0), i.e. in the window where friction work — and therefore this imbalance — is near zero by construction.

**Trigger:** Bites when Track A's later phases hit the skew floor: merging and overlap remove wire time, exposing the friction-driven late-sender imbalance, which the program has forbidden itself to fix. Bites hardest in exactly the regime the benchmark targets (through-rupture) and on the production SAFS meshes, whose fault partition imbalance is larger (C §4: fault imbalance 2.16 at np=256).

**Actual behavior:** "Fault & friction (0.65 %) · permanently out of scope", justified on a rank average from a window that cannot be separated into pre-rupture and rupture.

**Expected behavior:** Fault/friction closed for **kernel** optimisation (correct — 0.65 % of mean work) but explicitly **reopened as a load-balance item inside Track A**, gated on a rupture-window per-rank measurement.

**Suggested fix:** Replace :39-40 and the "Closed / withdrawn" entry with: "**Fault/friction: closed as a KERNEL target, reopened as a BALANCE target.** The 0.65 % is a rank average; the same Caliper table records Min 2.7 s / Max 254 s = a 94× spread, and the solver is bulk-synchronous so the critical path is the max rank — 2.86 % of wall today, ~10 % of the step after a 4× bulk speedup. The measurement also cannot be attributed to the rupture window (B0 Harness gap #1: windows not separable; rupture is 27 % of the averaged steps), so 'even through rupture' is withdrawn. Because this is the only measured mechanism that generates the late-sender skew this plan names as Track A's cap — and because the 24 % skew proxy that kept the fault-weighted partition off the critical path was measured pre-rupture at `--tfinal 0.5` — Track A's A0 must publish per-rank friction time over a rupture window, and the fault-weighted partition is promoted to co-requisite if the max-rank friction share exceeds the skew budget."

**Evidence:** A Leg 2 stage table (`friction_substep` 0.65 %, Min 2.7 s / Max 254 s, whole run to t=1.8, 4917 steps of which 1317 in the rupture window); A Harness gap #1 ("per-leg wall never captured … windows not separable"); A Leg 5 + `jobs/kernel_dev/run_p0_baseline_expanse.sbatch:190-192` (`--tfinal 0.5`, pre-nucleation); C §4 (nucleation t_0 = 1.0, V_max peak t ≈ 1.48; fault partition imbalance 2.16); D §Summary/§Phase 0 (fault-weighted partition as the >40 % escalation).

**Verification:** Completeness critic — unrefuted.

---

### [R-005] [CRITICAL] [PLAN_v2:104-108, :51, :155-163] — Mechanism choice and "honest ceiling" both presume MFEM is memory-bound rather than executing more FLOPs; neither track schedules the one cheap counter run that decides it

**Category:** ASSUMPTION (lens = completeness)

**Description:** Two load-bearing conclusions rest on an unmeasured either/or that C §2 states explicitly and refuses to resolve: MFEM "is either doing several× the FLOPs [per update] or running deeply memory/latency-bound scalar code". **No FLOP counter (PAPI, perf hardware counters) has ever been run on MFEM, on any machine, in any artifact** — the only MFEM GFLOP/s figures anywhere come from a single hard-coded `FLOPS_PER_ELEM = 71,280` in `tests/bench/bench_ader_kernel_variants.cpp:66-70`, and the "~0.4 GFLOP/s/core" in C §2 is 6.4 ÷ 17.1 computed *under the assumption of equal FLOP counts*. Yet: (a) Track B's selected mechanism is tiling/fusion, i.e. **traffic reduction**, which pays off only under the memory-bound branch (:51 asserts "the win is data movement, not arithmetic" on the strength of a synthetic contended microbench whose own control is invalid — three of four variants ran *faster* contended than uncontended); and (b) the Honest ceiling attributes the irreducible 3–4× to SeisSol's "hand-tuned assembly", which the evidence contradicts — SeisSol sustains 6.4 GFLOP/s/core = ~18 % of Rome core peak (C §2), which is ordinary well-vectorised compiled throughput, and the plan's own bench variant D reaches an implied 5.17–8.80 GFLOP/s/core contended **on the same chip**, i.e. an MFEM-side kernel in the evidence base already meets or exceeds SeisSol's sustained per-core rate. If the true cause of the 17× is a several-fold higher FLOP count per element update (algorithmic/formulation — whole-vector CK recursion, O+1 accumulators — versus SeisSol's fused ADER), then tiling cannot close it, the 3–4× ceiling is not about assembly, and the far-field order drop the plan relegates to "a separate proposal" is the **dominant remaining lever**, not an optional one. A single instrumented leg settles it and is scheduled by neither track.

**Trigger:** Bites if MFEM's measured FLOPs/update turn out to be within ~2× of SeisSol's 92 kFLOP/update — then the gap is bandwidth/latency and the plan is right; if they are 4–8×, both tracks can fully succeed and still leave a gap that neither mechanism addresses. That is the same class of error as the previous "far too wrong" LTS projection.

**Actual behavior:** Mechanism selected and ceiling explained on an unmeasured hypothesis; the decisive measurement is absent from the plan, from Track A's phases, and from Track B's phases.

**Expected behavior:** A hardware FLOP count per element update for MFEM on Rome, taken alongside B0, compared to SeisSol's counter-measured 92 kFLOP/update — **before** committing to traffic reduction as *the* mechanism and before accepting 3–4× as structural.

**Suggested fix:** Add to §Sequencing as a prerequisite alongside the Caliper scopes: "**B0b — MFEM FLOPs per element update (one leg, counters only).** RESULTS_p5 §2 states the 17× gap is 'either several× the FLOPs or deeply memory/latency-bound' and does not distinguish them; no FLOP counter has been run on MFEM anywhere in the evidence base (every MFEM GFLOP/s figure derives from the microbench's hard-coded 71,280 FLOP/elem constant). Track B's mechanism (tiling/fusion = traffic reduction) is only correct on the memory-bound branch. Run `perf stat` FP counters on one Expanse leg and compare against SeisSol's counter-measured 92 kFLOP/update (340.1 TFLOP / 3.710e9 updates). If MFEM's count is >2× SeisSol's, the mechanism is re-opened and the far-field order drop is promoted from 'a separate proposal' to the primary route." Separately, replace "their element kernels are hand-tuned assembly we are not attempting to match head-on" (:158) with the measured fact: SeisSol sustains 6.4 GFLOP/s/core ≈ 18 % of Rome core peak, and this program's own bench variant D reaches 5.2–8.8 GFLOP/s/core contended on the same chip — so per-core throughput is not the demonstrated barrier and the ceiling must be re-justified on a measured FLOP-count or bandwidth basis.

**Evidence:** C §2 ("either doing several× the FLOPs or running deeply memory/latency-bound scalar code"; SeisSol 340.1 TFLOP HW total, 92 kFLOP/update, 6.4 GFLOP/s/core ≈ 18 % of peak); `tests/bench/bench_ader_kernel_variants.cpp:66-70` (`FLOPS_PER_ELEM` hard-coded, applied identically to all variants and both machines); A Leg 4 (contended D 8.1–13.8 µs/elem ⇒ 5.17–8.80 GFLOP/s/core at that constant); B §2 (roofline built on an **assumed** 2 GB/s/core, no STREAM measurement anywhere in A–D).

**Verification:** Completeness critic — unrefuted.

---

### [R-006] [MODERATE] [PLAN_v2:12-22] — The headline budget (1538 / 3879 / "16× behind") is a GTS→LTS transfer the plan's own evidence contradicts

**Category:** UNSUPPORTED_CLAIM (lens = traceability)

**Description:** The compute term 1538 s/sim-s is not measured. It is C's LTS compute 2015 (itself a residual/ratio inference, 6465/3.208, face-cache OFF) divided by 1.31 — a face-cache activation factor measured on the GTS leg of a **different job** (52379131). The plan flags the GTS→LTS transfer only for the stage split (:76-82); it never flags that the 1538 itself is a GTS→LTS transfer. The one existing face-cache-ON LTS datapoint refutes it: A Leg 5's Waitall Avg 1157 s = 52.6 % of wall implies a wall of 2199.6 s over `--tfinal 0.5` = 4399 s/sim-s, i.e. 1.01× vs the face-cache-OFF LTS rate of 4356.

**Trigger:** Whenever the face-cache gain on the LTS leg is < 1.31× (Leg 5 indicates ~1.0×, and R-002 shows why), the budget, the 60/40 split, the 2.6× payoff ratio and the 16× gap are all wrong **in the direction that favours the plan's chosen ordering**.

**Actual behavior:** Plan asserts 3879 s/sim-s and 16× as measured fact.

**Expected behavior:** The only *measured* LTS wall on this benchmark is 4356 s/sim-s = 18.2× behind SeisSol's 239 (C §Summary); the only face-cache-ON LTS run implies ~4399.

**Suggested fix:** Replace :12-22 with: "**Where we actually are.** The last MEASURED end-to-end LTS rate on the benchmark is **4356 s/sim-s** (C §Summary, face-cache OFF) = **18.2× behind** SeisSol's as-run 239. With face-cache ON the GTS leg improved 1.31× (B0 Leg 2), but the only face-cache-ON LTS datapoint we have (B0 Leg 5 ⇒ ~4399 s/sim-s) shows **~1.0×**, i.e. no LTS gain. The projected 1538 + 2341 = 3879 (16×) therefore assumes a GTS→LTS transfer that our own LTS leg contradicts. **Retiring this is the first measurement of either track:** publish the Leg-5 wall and re-run the LTS leg with/without `--face-cache` before any budget number in this plan is quoted as fact."

**Evidence:** A Leg 5 (Waitall Min 877/Avg 1157/Max 1983 s = 52.6 % of wall, `--tfinal 0.5` ⇒ implied wall 2199.6 s ⇒ 4399 s/sim-s) vs C §3 (LTS 4356 s/sim-s, face-cache OFF); 1538 = 2015/1.31 with 1.31 from A Leg 2 (GTS leg).

**Verification:** 3/3 skeptics failed to refute.

---

### [R-007] [MODERATE] [PLAN_v2:24, :82, :104, :128] — "Fixing every kernel we have evidence for recovers 776": two of its three terms have no bench on any machine

**Category:** UNSUPPORTED_CLAIM (lens = traceability)

**Description:** 776 reconstructs exactly as 0.5623×1538×(1−1/4) + 0.1296×1538×(1−1/1.5) + 0.0788×1538×(1−1/2) = 648.6 + 66.4 + 60.6 = 775.6. The face-interior **1.5×** and volume **2×** factors exist in no measurement: the Leg-4 bench is predictor-only by construction (bench header: "Simulates the ADER-O4 p3 predictor"; no face-flux stage, no volume RHS, no corrector). B3's gate "≥1.5× on that stage" is likewise a target with no supporting bench. The phrase "every kernel we have evidence for" therefore misdescribes two of the three terms, and the 4× predictor factor is itself a chosen point inside a synthetic contended range.

**Trigger:** Comparing Track B's 776 against Track A's 2041 to justify the track ordering.

**Actual behavior:** 776 presented as the evidence-backed kernel recovery.

**Expected behavior:** Only the predictor term (648.6, and that at an unvalidated 4×) has any bench behind it; face and volume contribute 127 s/sim-s of pure assumption.

**Suggested fix:** Replace :24's second clause with: "fixing the predictor at the bench's ~4× recovers **~649** (the only kernel term with any bench behind it); adding *assumed* face-interior 1.5× and volume 2× factors — **neither stage has ever been benchmarked, on any machine** — would add 127 for a total of 776." Add to §B3: "*Gate:* ≥1.5× on that stage. **No bench supports 1.5×; B3 must be preceded by a volume-stage bench variant, or its 61 s/sim-s contribution dropped from the Track-B total.**" Also reconcile the internal inconsistency that 776 credits volume at 2× while B3's gate is ≥1.5× — a gate passing exactly at threshold delivers 40.4 s, not 60.6 s.

**Evidence:** A "Projections from the measured split" ("predictor at bench factors, face 1.5×, volume 2×"); `tests/bench/bench_ader_kernel_variants.cpp:8-30, :54-70` — predictor replica only; no face or volume kernel exists in the bench.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* the third skeptic showed 776 gates nothing (grep: it appears only at :24, :82, :104, :111 and in no gate or go/no-go), that B4's projected saving at the disputed 1.5× is 4.32 % of the step — below the plan's own "<5 % ⇒ drop it" rule at :132, so the plan already schedules that term for deletion — and that stripping both unbenched terms *strengthens* the plan's own comm-primary conclusion (ratio 2.63× → 3.14×). Retained at MODERATE as a provenance defect on a headline sentence, not as a decision defect.

---

### [R-008] [MODERATE] [PLAN_v2:142-144, :67-74] — "Prerequisite (DONE): the two LTS Caliper scopes, so both tracks can measure their own gates" is false: B3's, B4's, the LTS corrector's and the LTS friction stage's regions do not exist

**Category:** DEVIATION (lens = code)

**Description:** Commit `c547ac8` added exactly two `MFEM_PERF_SCOPE` lines (`wave_operator.inl:1929`, `:2611`). The complete scope inventory on the LTS path is therefore: `ComputeADERSubStepStatesAndIntegralCluster` (:1929), `ComputeADERClusterSeamFaceFluxRHS` (:2611), and `seas::spatial_dyn::step` (driver:4854 — **and that one is GTS-only**, being inside `for (…; !lts_stepping && !lts_fault_stepping; …)` at driver:4852, so the LTS legs have no step-region denominator at all). Everything else the plan gates against is uninstrumented: (a) `AdvanceADERClusterBulk` (:2190), the LTS corrector, whose GTS sibling `AdvanceADER` (:7239) supplies the 34.40 % parent share; (b) `ComputeVolumeRHSElems_` (:2068) — **B3's exact target**, whose GTS sibling supplies the 7.88 %; (c) the LTS interior-face sweep — **B4's target**, written inline at :2229-2470 with no function boundary at all; (d) `ApplySpatialDerivativeElems_` (:1802), the LTS analogue of the child carrying 27.14 pp of the 56.23 %; (e) the LTS friction stage — `iterator_.Advance(...)` at `lts_fault_stepper.hpp:283` has no region, so `seas::spatial_dyn::friction_substep` has no LTS counterpart. `grep MFEM_PERF lts_*.hpp lts_*.cpp` returns nothing. Consequently the promised "re-measured LTS stage split" can report at most 2 of ~6 stages, and B3's "≥1.5× on that stage", B4's "<5 % of the step" re-justification, and the closure of fault/friction at "0.65 %" all remain **unmeasurable on the leg the budget uses**.

**Trigger:** The first LTS run of A1 or B1 — the run the plan designates as retiring the GTS-proxy caveat and publishing the re-measured LTS stage split.

**Actual behavior:** Plan declares the instrumentation prerequisite DONE with two scopes.

**Expected behavior:** Six scopes are needed for the LTS split to be structurally comparable to the GTS table the budget is derived from, plus a top-level region on the LTS stepping loops so LTS shares have a step denominator rather than whole-program time.

**Suggested fix:** Add these before the next Expanse job, then change :142 to "Prerequisite (IN PROGRESS)": (1) `wave_operator.inl:2204` — `MFEM_PERF_SCOPE("seas::WaveOperator::AdvanceADERClusterBulk");` (mirrors `AdvanceADER` :7239); (2) `:2071` — `MFEM_PERF_SCOPE("seas::WaveOperator::ComputeVolumeRHSElems_");`; (3) `:1805` — `MFEM_PERF_SCOPE("seas::WaveOperator::ApplySpatialDerivativeElems_");`; (4) wrap `:2229-2470` in a braced block opened by `MFEM_PERF_SCOPE("seas::WaveOperator::ClusterInteriorFaceSweep");`; (5) bracket `lts_fault_stepper.hpp:283` with `MFEM_PERF_BEGIN/END("seas::spatial_dyn::friction_substep")`; (6) add a top-level region to the LTS stepping loops at `drivers/spatial_dyn_driver.cpp:4622` and `:4761`.

**Evidence:** `git show c547ac8 --stat` (2 insertions in `wave_operator.inl`); grep of `MFEM_PERF_SCOPE` across `dynamic/` + `drivers/` (17 hits, none on `AdvanceADERClusterBulk` / `ComputeVolumeRHSElems_` / `ApplySpatialDerivativeElems_` / any `lts_*.hpp`); `wave_operator.inl:2068, 2190, 2229-2470, 1802`; `lts_fault_stepper.hpp:283`; `drivers/spatial_dyn_driver.cpp:4852-4854`.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* the third skeptic confirmed the inventory but argued the currently-sequenced gates (B1, B2, Track A) are all measurable with what exists — B1/B2 target the instrumented predictor, Track A gates on Caliper's MPI service — and that B3/B4 are demoted and can have their scopes added alongside their own implementation at zero marginal job cost. That skeptic also supplied item (6) above by noticing that `seas::spatial_dyn::step` is GTS-only. Retained at MODERATE because the :142-144 claim, as written, is what licenses running the first gate leg.

---

### [R-009] [MODERATE] [PLAN_v2:5, :18-22] — "Evidence base (all measured, none assumed)" is false for the two headline numbers, and the contradicting Leg-5 datapoint is omitted

**Category:** UNSUPPORTED_CLAIM (lens = omission)

**Description:** The header claims the evidence base is "all measured, none assumed", and the budget table presents 1538 (compute) and 2341 (comm) as measured facts. Neither is measured. **2341** is a residual by subtraction in C §3 (4356 total minus an inferred compute term), and the same document publishes the exposed-wait share as 52.5 %, which gives 2287 — a 54 s/sim-s discrepancy never reconciled. **1538** is 2015/1.31, where 2015 is itself the GTS total 6465 divided by a 3.208 ratio (C §4 "LTS compute alone = 3.21×") and 1.31 is a face-cache activation gain measured on the **GTS** leg. Worse, the one existing face-cache-ON LTS measurement is omitted entirely: A Leg 5 implies ~4399 s/sim-s ≈ 1.0×, not 1.31×. The plan neither cites nor rebuts it.

**Trigger:** A reader takes 1538 + 2341 = 3879 and the derived 60/40 split as measured; the entire track ranking (2041 vs 776) is computed from them.

**Actual behavior:** :5 "Evidence base (all measured, none assumed)"; :20-21 table rows "2341 / 60 %" and "1538 / 40 %" with no derivation and no Leg-5 cross-check.

**Expected behavior:** The derivation shown, the transfer flagged with a direction argument, and the contradicting Leg-5 datapoint addressed.

**Suggested fix:** Change :5 to "**Evidence base:**" (drop "all measured, none assumed") and insert after the budget table (:22): "*Provenance of these two numbers (neither is a direct measurement).* **2341** is a residual: RESULTS_p5 §3 reports LTS total 4356 s/sim-s and derives compute as 6465/3.208 = 2015, leaving 2341 by subtraction; the same document's Waitall share of 52.5 % gives 2287 instead, a 54 s/sim-s (2.4 %) ambiguity in the Track-A base. **1538** = 2015/1.31, applying a face-cache activation gain measured on the **GTS** leg (B0 Leg 2) to an **LTS** compute term. **Contradicting datapoint, unresolved:** B0 Leg 5 is the only face-cache-ON LTS run in existence; it implies ~4399 s/sim-s, i.e. **~1.0×** against the 4356 face-cache-OFF baseline. Windows differ (0.5 s pre-rupture vs 2.0 s through rupture), so this is indicative, not decisive — but until B0-LTS resolves it, 1538 may be understated by up to ~30 % and the 60/40 split correspondingly overstated in Track A's favour."

**Evidence:** PLAN_v2:5, :18-22; C §3 (4356 total, 2015 compute, 2341 overhead, "52.5 % of wall") and §4 ("LTS compute alone = 3.21×"); A Leg 5 and Leg 2.

**Verification:** 3/3 skeptics failed to refute.

---

### [R-010] [MODERATE] [PLAN_v2:59-65] — The "CONFIRMED" verification claim — "every exchange in the wave operator lives in exactly two places" — is false; there are four runtime exchange sites, one firing 4 blocking full-halo rounds per GTS step inside the unaccounted 8.72 pp

**Category:** UNSUPPORTED_CLAIM (lens = completeness)

**Description:** The plan states: "Every exchange in the wave operator lives in exactly two places — `ComputeADERSharedFaceFluxRHS` and `ComputeADERClusterSeamFaceFluxRHS` — i.e. precisely the regions this plan assigns to Track A." Grep of `dynamic/wave_operator.inl` returns **six** runtime `ExchangeFaceNbrData()` call sites in **four** functions: `:2591` `EnsureGhostClusterIds_`, `:2705` and `:2729` `ComputeADERClusterSeamFaceFluxRHS`, `:3714` `EvaluateBulkAtFaultQPsCanonical`, `:4847` `ComputeSharedFaceFluxRHS`, `:6521` `ComputeADERSharedFaceFluxRHS`. The omission that matters is `:3714`: it is called from `drivers/spatial_dyn_driver.cpp:709-713` in the `for (o = 0; o < O; ++o)` loop, i.e. **O = ader_order = 4 blocking full-halo vdim=9 exchanges per GTS step**, gated on `pmesh.GetNSharedFaces() > 0` (:3677-3680) so **every** rank participates. Those four rounds sit outside the predictor's Caliper scope and outside `friction_substep` — they land in the **8.72 pp** of the GTS step that the stage split leaves unattributed (100 − 56.23 − 34.40 − 0.65). The plan never mentions that its shares do not close, and its direction argument at :80-82 accounts only for face-shared. So the audit that legitimises (a) reassigning the 10.56 % shared-face stage to Track A and (b) declaring the predictor's 56 % to be pure work was performed with an **incomplete inventory**, and 8.72 pp of the step — known to contain 4 blocking collectives per step — is silently treated as neither comm nor kernel.

**Trigger:** Bites when either track re-measures the split and finds shares that do not sum to 100, and when Track A scopes its addressable surface: on the GTS leg its true share is 10.56 % **plus an unmeasured part of 8.72 pp**. Also bites D's Phase 4, whose GTS retrofits target exactly these unlisted sites.

**Actual behavior:** Two exchange sites claimed; four exist. The 8.72 pp residual is never named, never bounded, and implicitly assumed un-speedable non-comm work.

**Expected behavior:** A complete, line-cited exchange inventory (four functions, six sites), an explicit statement that the GTS shares close only to 91.28 %, and an assignment of the residual — in particular the O=4 per-step fault-QP exchange at `:3714` — to a track.

**Suggested fix:** Replace :62-63 with: "Runtime halo exchanges live in FOUR functions, not two: `ComputeADERClusterSeamFaceFluxRHS` (`wave_operator.inl:2705, :2729` — LTS), `ComputeADERSharedFaceFluxRHS` (`:6521` — GTS ADER), `EvaluateBulkAtFaultQPsCanonical` (`:3714` — GTS ADER, O = ader_order = 4 blocking full-halo rounds per step from driver:709-713, gated on `GetNSharedFaces()>0` so all ranks participate), and `ComputeSharedFaceFluxRHS` (`:4847` — explicit-RK, 9 scalar rounds). `EnsureGhostClusterIds_` (`:2591`) adds one setup exchange. Consequence: the GTS stage shares close only to **91.28 %**; the residual 8.72 pp contains at least 4 blocking collectives per step and is NOT assumed un-speedable. Re-measure with a Caliper scope on `EvaluateBulkAtFaultQPsCanonical` and assign the residual before either track's budget is treated as final."

**Evidence:** grep of `dynamic/wave_operator.inl` (sites :2591, :2705, :2729, :3714, :4847, :6521 in four distinct functions — reproduced during this audit); `wave_operator.inl:3658-3714` (R-1600 gate on `GetNSharedFaces`, one batched collective per invocation); `drivers/spatial_dyn_driver.cpp:709-713` (O-loop) and :725-733 (`friction_substep` scope excludes it); A Leg 2 stage table (56.23 + 34.40 + 0.65 = 91.28).

**Verification:** Completeness critic — unrefuted.

---

### [R-011] [MODERATE] [PLAN_v2:29-32, :94, :151-153] — Track A and Track B are anti-synergistic: A3's payoff is hiding wire time behind compute that B1/B2 exist to delete, yet the plan composes them by simple subtraction

**Category:** ASSUMPTION (lens = completeness)

**Description:** The expected outcome 760–1060 is formed as compute 1538 − kernel 776 + residual comm (0–300), i.e. the two tracks' savings are treated as independent and additive. They are not. Track A's phase A3 is explicitly "hide wire time behind the tick's ghost-independent compute" (:94), so its recovery is **bounded by the compute available to hide behind**. Track B's entire purpose is to remove that compute — 776 s/sim-s of it, mostly from the predictor, which is precisely the ghost-independent work in the hiding window. Halving the hideable compute directly reduces A3's achievable overlap, so the residual exposed wait after **both** tracks is strictly larger than the residual measured after Track A alone. The Re-baseline rule makes this undetectable until the end: "each track's gate is measured with the other track's flags OFF, and the composed number is measured once at the end" — guaranteeing each track is gated in the configuration that maximises its own apparent payoff, and deferring the only measurement that can expose the interaction to a point where both tracks are built and their gates already recorded as PASS. Compounding it, the source shows the hiding window is already severely bounded independently of Track B: corrects run fine→coarse (`lts_stepper.hpp:130-133`, `lts_stepper.cpp:20-36`), so an `End()` before the tick's first seam sweep leaves clusters 1..Nc-1's corrects **after** the completion point, and on the 16 of 32 ticks that activate only cluster 0 the hideable work is at most ~2 × 1.09 % of the mesh's updates.

**Trigger:** Bites at the very end of the program, when the composed run is measured for the first time: both tracks can report PASS on their own gates and the composed number still miss 760–1060.

**Actual behavior:** Composition is a subtraction of two independently measured savings; no interaction term, no ordering constraint, no intermediate composed checkpoint.

**Expected behavior:** An explicit interaction statement plus at least one composed measurement before B2 is built — e.g. gate A3 with B1's flag ON, since B1 lands first and already shrinks the hiding window.

**Suggested fix:** Amend the Re-baseline rule (:151-153) to: "Each track's gate is measured with the other track's flags OFF **for attribution**, but A3 must ALSO be measured with B1 ON, because A3's payoff is overlap behind ghost-independent compute and B1/B2 delete that compute. The two tracks are anti-synergistic: the composed saving is strictly less than 776 + 2041, and the composed number must be measured at the A3/B1 boundary, not only at the end of the program. Additionally, record the hideable-compute-per-tick-class measurement A0 was supposed to produce: corrects run fine→coarse (`lts_stepper.hpp:130-133`), so clusters 1..Nc-1's corrects fall after any `End()` placed before the tick's first seam sweep, and 16 of 32 ticks activate only cluster 0 (≤1.09 % of elements) — A3's window may be small before Track B touches it."

**Evidence:** PLAN_v2:29-32, :94, :151-153; `dynamic/lts_stepper.hpp:130-133` (ascending correct order); `dynamic/lts_stepper.cpp:20-36` (predicts then corrects); C §1 cluster histogram (cluster-0 bound ≤ 26,891 / 2,464,689 = 1.09 %); D §Payoff (overlap rows are projection only).

**Verification:** Completeness critic — unrefuted.

---

### [R-012] [MODERATE] [PLAN_v2:151-153, :119, :125] — B0 = 187.5 µs·core is pinned as "the compute denominator" but is contaminated by ParaView output, Caliper and `perf record` on rank 0

**Category:** QUALITY (lens = completeness)

**Description:** The Re-baseline rule pins "B0 = 187.5 µs·core is the compute denominator", and B1/B2 are gated as ratios against it. That denominator is the whole 8873.95 s Leg-2 wall divided by element updates, and that wall carried, unbounded and unsubtracted: `--paraview --paraview-fault-vtu --paraview-max-snapshots 40` (sbatch L112), Caliper runtime-report with `calc.inclusive` + `region.count` (L124), and `perf record -F 99 -g` on rank 0 (L131-134). C §3 puts Caliper alone at ~1–3 %, and the sbatch itself (L113-114) **deliberately dropped `profile.mpi` from this leg** on the reasoning that instrumenting one rank of a bulk-synchronous solver makes it the straggler — which is exactly what `perf record` on rank 0 then does. **Direction of the error: all three overheads inflate B0**, so a gate expressed as speedup-vs-B0 is systematically *easier* to clear than the true kernel gain. Second, the denominator's identity is wrong for the target: B0's predictor is `ComputeADERSubStepStatesAndIntegral` (GTS whole-vector path) while B2's stated target is "tiles inside LTS clusters", i.e. `ComputeADERSubStepStatesAndIntegralCluster` — a structurally different routine (indexed element subsets, 20-double blocked lambdas, plus a `retain` copy of every D(k) for seam providers).

**Trigger:** Every Track B gate run. Bites hardest if the gate run uses different instrumentation than B0, in which case the measured ratio mixes an instrumentation delta with the kernel delta — the same confound that already contaminates the headline 1.31× face-cache claim (187.5 carries perf + 40 snapshots; the 245.8 reference carries `profile.mpi` + spot and 300 snapshots).

**Actual behavior:** An instrumented, IO-carrying, GTS whole-run figure used as the reference for LTS-targeted kernel gates.

**Expected behavior:** A clean B0 — one leg with ParaView off (or a bounded fixed cadence), Caliper `region.count` only, no perf — published as B0_clean alongside the instrumented figure; and a separate LTS-leg predictor denominator for B2.

**Suggested fix:** Replace the Re-baseline rule's first sentence with: "B0 = 187.5 µs·core is NOT a clean denominator — it includes `--paraview-fault-vtu` (40 snapshots), Caliper runtime-report, and `perf record -F 99` on rank 0, none of which was bounded, and rank-0 perf makes rank 0 the straggler in a bulk-synchronous solver (the same reason `profile.mpi` was dropped from this leg). All three inflate B0, so ratios against it overstate kernel gains. Before B1's gate, publish **B0_clean** (same deck, ParaView off, Caliper `region.count` only, no perf) as the denominator, and publish an **LTS predictor denominator** measured on `ComputeADERSubStepStatesAndIntegralCluster` — B2 targets that routine, not the GTS `…AndIntegral` that B0 timed. Every gate run must use identical instrumentation to its denominator."

**Evidence:** `jobs/kernel_dev/run_p0_baseline_expanse.sbatch` L112, L113-114, L124, L131-134; A Leg 2 (8873.95 s wall = the 187.5 denominator; Caliper `step` = 99.99 % of total); C §3 (Caliper ~1–3 %); `dynamic/wave_operator.inl:1955-2057` (LTS predictor is a structurally different routine).

**Verification:** Completeness critic — unrefuted.

---

### [R-013] [MODERATE] [PLAN_v2:26-27, :97-98] — "Bit-for-bit identical by construction … which is why this track is both bigger and safer" asserts safety on the wrong axis: this code's recorded failure mode for exchange changes is a production-only collective deadlock

**Category:** ASSUMPTION (lens = completeness)

**Description:** The plan justifies inverting track order partly on risk: comm is "bitwise-safe", kernels "need tolerance gates". That comparison is incomplete in a way the source contradicts. The historical, documented failure of an exchange-gating change in this exact code is not numerical drift — it is a **matched-collective deadlock**: `wave_operator.inl:3661-3672` records R-1600, where gating a per-substep exchange on the fault-only shared count instead of `pmesh.GetNSharedFaces()` produced "the production-blocking 13.5-min hang at np=10". D's own risk 3 names the same class. Track A's phases A1–A3 change exactly the property that governs collective matching — which ranks post how many exchanges on which tick — and A3 additionally introduces persistent requests (`MPI_Send_init`/`Recv_init` + `Startall`) whose correctness depends on buffer-address stability that MFEM does **not** provide (`fem/pgridfunc.cpp:243-246` re-selects the send/recv buffers at every call via `Read()`/`HostRead()`), and on tag/communicator isolation (`int tag = 0` on `pfes->GetComm()`, `:250-262`). None of that is caught by a bitwise-identity gate; a deadlock produces no output to compare. Worse, the declared bitwise gate runs at np ∈ {1,2,4,10}, and **np=1 executes zero exchanges** (both LTS sites early-return on `GetNSharedFaces()==0`), so one of the four gate configurations has no exchange coverage at all, and none reaches the np=256 partition where a rank holds multiple cluster seams simultaneously.

**Trigger:** The first np=256 Expanse gate run of A1/A2/A3: a collective-count mismatch or a persistent-request buffer invalidation hangs the job for its full walltime, consuming an allocation and producing no data — while the local np≤10 bitwise gates all passed.

**Actual behavior:** Risk characterised solely as "bitwise identity vs tolerance gates", with no deadlock/liveness gate and no acknowledgement that np=1 exercises no exchange.

**Expected behavior:** An explicit liveness gate for A1–A3 alongside the bitwise gate.

**Suggested fix:** Amend :97-98 to: "**Gate for the whole track:** exposed wait ≤ ~15–20 % of wall. **Acceptance is bitwise identity for A1–A3 PLUS a liveness gate** — bitwise identity cannot detect the failure mode this code has actually suffered from exchange changes: `wave_operator.inl:3661-3672` records the R-1600 unmatched-collective hang (13.5 min at np=10, production-blocking). Every A-phase commit must therefore (i) update `BuildTickTable.n_collectives`, the stepper and the wave counter atomically and keep the P-007 audit unrelaxed, (ii) run a timeout-guarded np=10 liveness check (note np=1 executes **zero** exchanges — both LTS sites early-return on `GetNSharedFaces()==0` — so it provides no exchange coverage), and (iii) for A3, assert buffer-address stability before `Startall`, since MFEM re-selects the send/recv buffers on every call (`fem/pgridfunc.cpp:243-246`) and uses `tag = 0` on the shared communicator (`:250-262`), so persistent requests must own their buffers and an isolated communicator."

**Evidence:** `dynamic/wave_operator.inl:3661-3672` (R-1600, 13.5-min production hang from a mis-gated exchange); `fem/pgridfunc.cpp:243-246`, `:250-262`; D §Risks item 3 and §Phase 1 gate (np ∈ {1,2,4,10}); `drivers/spatial_dyn_driver.cpp:4639-4641, 4781-4783`.

**Verification:** Completeness critic — unrefuted.

---

### [R-014] [MODERATE] [PLAN_v2:40-41, :136 vs :106-107] — "The BLAS question (already optimal, statically linked dgemm — no win exists)" is closed on a symbol grep with zero timing, and contradicts the plan's own mechanism list

**Category:** UNSUPPORTED_CLAIM (lens = completeness)

**Description:** The plan lists the BLAS question as closed — "already optimal, statically linked dgemm — no win exists". The evidence is A's Leg 1, which performed `nm`/`ldd`/`config.mk` inspection on the job-52379131 binary and took **no timing of any kind**: no A/B against an internal-LU build, no alternative BLAS, and no inspection of the Phase-5 binary that produced the 245.8/245.9 reference (its git hash and `MFEM_USE_LAPACK` setting are recorded nowhere in A or C), so Leg 1's own conclusion — "the 245.8 was already on the optimized-dgemm LU path" — is an inference about a binary that was never examined. "No win exists" is a performance claim; nothing measured supports it. Independently, the closure sits ~65 lines above the Track B mechanism paragraph, which names "per-element dgemm (variant D, 2.8–8.3×)" as a co-equal selected mechanism — the same document declares dgemm both a closed no-win and one of its two chosen levers. A further unremarked hazard: the linked OpenBLAS is PETSc's LP64 static build; variant D's per-element 20×20 dgemm calls incur per-call dispatch overhead at production shape, which the synthetic bench (NE=20,000 contiguous, no virtual calls, no `MFEM_VERIFY`) does not model.

**Trigger:** A reader takes "no win exists" at face value and skips the one measurement that would price the alternative; or B2 selects variant D and hits the dgemm path's per-call overhead at the production per-rank element count (9,628 vs the bench's 20,000), which no artifact measures.

**Actual behavior:** A build-configuration inspection with no timing reported as a closed performance question, alongside a mechanism that credits the same library with 2.8–8.3×.

**Expected behavior:** Either a measured A/B (internal LU vs dgemm on the same binary and deck), or a narrower closure: relink is not a lever; the dgemm **call shape** still is.

**Suggested fix:** Replace the closed-item text with: "the BLAS **relink** question — the deployed build already links PETSc's static LP64 OpenBLAS with a defined `dgemm_`, so swapping BLAS implementations is not a lever. Note this was established by symbol/config inspection only (B0 Leg 1 took no timing, and the binary behind the 245.8 reference was never inspected), so 'no win' is a build fact, not a measured one. The dgemm **call shape** is NOT closed: bench variant D (per-element dgemm) is one of Track B's two candidate mechanisms at 2.8–8.3× contended, measured on a synthetic replica with 20,000 contiguous elements per copy versus 9,628 per production rank and with no virtual-call or verify overhead — its per-call dispatch cost at production shape is unmeasured, and D is **not** bit-identical to the current path (dgemm reassociates), which conflicts with B2's 'D(k) still scattered to `dk_retain` bit-identically (seam contract)'."

**Evidence:** A Leg 1 (`nm`/`ldd`/`config.mk` inspection only, explicitly no timing; Phase-5 binary's LAPACK setting not recorded in A or C); PLAN_v2:106-107 vs :41 and :136; `tests/bench/bench_ader_kernel_variants.cpp:32-34` (only A and B are bit-identical by construction; D agrees "to round-off"), `:54` (NE=20000) vs 2,464,689/256 = 9,628 elems/rank.

**Verification:** Completeness critic — unrefuted.

---

### [R-015] [MODERATE] [PLAN_v2:23-27, :44-53] — The track-ordering conclusion is generalised from one 2-node / 256-rank / 9,628-elements-per-rank configuration, and the split that drives it is the quantity most sensitive to that configuration

**Category:** ASSUMPTION (lens = completeness)

**Description:** "Communication is the primary program" is stated unconditionally and is the sole justification for inverting the phase order. It is derived from exactly one configuration: TPV104-200m, 2,464,689 tets, 2 nodes × 128 ranks, **9,628 elements/rank**. The comm share is a strong function of that last number — halo surface scales as (elems/rank)^(2/3) while compute scales linearly, so the 60/40 split shifts toward compute as elements per rank grow and toward comm as ranks grow. Nothing in A, B, C or D measures the split at any other rank count, node count, or mesh. The **only** configuration variation attempted was A's Leg 3 (64 ranks/node, half occupancy), declared unusable because its wall was never printed — even though A's own job arithmetic bounds it (≤ 225.9 µs·core unconditionally, ≤ 142.2 µs·core after subtracting the implied Leg-5 wall, i.e. a possible ≥ 1.3× half-occupancy gain nobody computed). The plan's branch is `safs-v4_0_0-alt-case1-mfem-speed` and this project's production meshes are 34.9–52.3 M tets (14–21× TPV104-200m); at comparable node counts that is ~200,000 elements/rank, a regime in which the measured 60 % comm share has no evidential support and could plausibly invert the ordering this document is built on.

**Trigger:** Bites when the program is applied to production: if the comm share at production scale is materially below 60 %, Track A's 2041 shrinks proportionally while Track B's 776 does not.

**Actual behavior:** A single-configuration measurement stated as a universal property of the solver, with no scaling caveat, no scaling gate, and the one available scaling datum (Leg 3) abandoned rather than bounded.

**Expected behavior:** Either a scaling statement bounding the regime in which the 60/40 split holds, or a cheap two-point scaling measurement folded into the first gate run.

**Suggested fix:** Add to the Summary after :27: "**Scope of this conclusion.** The 60/40 split is measured at one configuration only: TPV104-200m, 2 nodes × 128 ranks, 9,628 elements/rank. Comm share grows with rank count and falls as elements/rank grows (halo ~ (elems/rank)^(2/3) against linear compute), so the track ordering is conditional on that regime. The only configuration variation attempted (B0 Leg 3, 64 ranks/node) was abandoned as unusable although it is boundable from the job wall at ≤ 225.9 µs·core — recompute it. Fold a two-point rank sweep (128 / 256 / 512 ranks, same mesh) into A1's first gate run and re-state the split; if the comm share at the production element-per-rank count of the SAFS v4_0_0 meshes falls below ~45 %, the track ordering is re-opened."

**Evidence:** A header (2 nodes × 128 ranks, 2,464,689 tets ⇒ 9,628 elems/rank), Leg 3 (declared not computable; bound derivable from the 14,815 s job wall minus Leg 2's 8873.95 s); C §1 (single mesh, single rank count for every MFEM figure); PLAN_v2:3 (branch `safs-v4_0_0-alt-case1-mfem-speed`); project record: SAFS v4_0_0 faultband meshes 34.87 M / 47.1 M / 52.3 M tets.

**Verification:** Completeness critic — unrefuted.

---

### [R-016] [LOW] [PLAN_v2:12-14] — Headline "3879 s/sim-s" is a projection, and the only face-cache-ON LTS measurement implies ~4399

**Category:** UNSUPPORTED_CLAIM (lens = arithmetic)

**Description:** 3879 = 1538 + 2341 reproduces, and 1538 = 2015/1.31 reproduces to 1538.17. But 2015 is C's LTS-leg compute (face-cache OFF, t = 0..2.0) and 1.31 is A Leg 2's **GTS**-leg activation gain. No run has ever produced 3879 s/sim-s. The single existing face-cache-ON LTS datapoint is A Leg 5: wall = 1157/0.526 = 2199.6 s over `--tfinal 0.5` = 4399 s/sim-s versus C's face-cache-OFF LTS 4356 — a gain of 1.010×, not 1.31×.

**Trigger:** Any reader or downstream gate treating 3879 as the measured "where we actually are" baseline.

**Actual behavior:** 3879 s/sim-s presented as measured fact, labelled "16× behind" (3879/239 = 16.23).

**Expected behavior:** The measured LTS wall is 4356 s/sim-s = 18.2× behind 239; 3879 is a projection resting on a GTS→LTS transfer the one available LTS datapoint puts at 1.01×.

**Suggested fix:** Replace :12-14 with: "**Where we actually are.** The last *measured* LTS wall on the benchmark is **4356 s/sim-s** (RESULTS_p5, face-cache OFF); SeisSol as-run is 239 — **18.2× behind**. Enabling `--face-cache` measured **1.31×** on the *GTS* leg (B0 Leg 2), which if it transferred would give **3879 s/sim-s (16×)**. It is not yet shown to transfer: the only face-cache-ON LTS datapoint (B0 Leg 5 ⇒ **4399 s/sim-s**) shows **~1.01×**. Every budget below is quoted against the projected 3879; re-measure the LTS wall with `--face-cache` on the first run of either track and restate if it lands near 4356."

**Evidence:** A Leg 5 ("Min 877 / Avg 1157 / Max 1983 s = 52.6 % of wall") with sbatch L190-192 `--tfinal 0.5` ⇒ 2199.6 s wall ⇒ 4399 s/sim-s; vs C §Summary ("MFEM clustered-LTS 4356 s/sim-s").

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* the third skeptic showed the 4399 is init-inclusive (Leg 5 runs the `lts_fault_stepping` block at driver:4700, which has **no** top-level Caliper scope, so its "% of wall" denominator is whole-process time) while C's 4356 is stepping-only; C's own job walls put non-stepping at 467–800 s, and 2199.6 − 1662 = 538 s sits squarely in that band, i.e. Leg 5 is *consistent* with a full 1.31×. That dissent is well-founded and is why this instance is LOW — but note it is independent of **R-002**, which kills the transfer on source-structural grounds regardless of Leg 5.

---

### [R-017] [LOW] [PLAN_v2:24 vs :34-37] — The 2041 comm lever contradicts the plan's own stated 24 % skew cap in the same Summary

**Category:** ARITHMETIC (lens = arithmetic)

**Description:** 2041 = 2341 − 300, i.e. a residual exposed wait of 300 s/sim-s = 12.8 % of 2341. Four lines later the plan states 24 % of the wait is skew that no communication change touches and that this caps Track A. 24 % of 2341 = 562, so the plan's own cap makes the maximum comm recovery 2341 − 562 = **1779** s/sim-s.

**Trigger:** Reading the Summary as a self-consistent budget; any downstream use of 2041 as Track A's payoff.

**Actual behavior:** Lever quoted as 2041 (residual 12.8 %) while the stated cap is 24 %.

**Expected behavior:** ≤ 1779 s/sim-s if the 24 % skew is genuinely untouchable, or an explicit statement that merging/overlap also removes part of the skew (which D explicitly denies: "skew is untouched by merging or overlap").

**Suggested fix:** Change :24 to: "Fixing the waiting recovers **~1780–2040 s/sim-s** (the upper end assumes the 24 % skew proxy is partly recoverable by merging; the lower end applies it as a hard floor of 562 s/sim-s)" and retitle §Track A accordingly.

**Evidence:** A Leg 5 ("skew proxy (Avg−Min)/Avg = 24 %"); 2341 × 0.24 = 561.8; D §Summary Main risk ("skew is untouched by merging or overlap").

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* 2041 appears only at :24 and :84 and gates nothing (Track A's gate is a wall fraction; Track B's are stage ratios); D:114-116 in fact says merging + overlap *do* attack soft-synchronisation skew; C §4 measures GTS element balance at Max/Avg = 1.0000, so the spread is exchange-structural, not element imbalance; and the plan already flags the cap as its "Biggest risk" and schedules a finer per-tick split during A1. Correctly a band-vs-point wording issue.

---

### [R-018] [LOW] [PLAN_v2:29-30] — Expected-outcome band 760–1060 uses two different comm residuals at its two endpoints; the low end assumes zero residual wait

**Category:** ARITHMETIC (lens = arithmetic)

**Description:** The high end reproduces as 3879 − 2041 − 776 = 1062 (residual 300). The low end reproduces only as 1538 − 776 = 762, i.e. residual **0**. Zero is outside D's own stated floor of 150–400 s/sim-s. On D's floor the band is 912–1162, i.e. 3.8–4.9× behind 239.

**Trigger:** Quoting "3.2×" as the achievable best case.

**Actual behavior:** "~760–1060 s/sim-s = 3.2–4.4×".

**Expected behavior:** "~912–1162 s/sim-s = 3.8–4.9×" on D's own residual floor.

**Suggested fix:** Replace :29-31 with: "**Expected outcome.** Both tracks land us at **~912–1162 s/sim-s = 3.8–4.9× behind SeisSol as-run**, from 16× today. The band's endpoints use the companion comm plan's own residual-exposed-wait floor of **150–400 s/sim-s**; a zero-residual best case is not supported by that plan."

**Evidence:** D §Payoff arithmetic ("merge+overlap ~150–400 exposed"); plan low end 762 = 1538 − 776 ⇒ residual 0; 762+150 = 912, 762+400 = 1162; 912/239 = 3.82, 1162/239 = 4.86.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* the band gates nothing (Track A gates on ≤15–20 % of wall, which itself encodes a nonzero residual), the 0–300 convention is inherited from A's own composed end state (669–969) and from B §3's R-401 table, and the plan is already *more* conservative than A on both endpoints. Correctly a one-line annotation.

---

### [R-019] [LOW] [PLAN_v2:49, :20] — "comm 60 % of wall, 2.6× the payoff" is listed under "Phase 0 measured" but neither number was measured

**Category:** UNSUPPORTED_CLAIM (lens = traceability)

**Description:** The column is headed "Phase 0 measured". The 60 % share is 2341/3879, whose denominator is the unmeasured 1538; every **measured** exposed-wait share in the evidence base is 52.5 % (C §3) or 52.6 % (A Leg 5). The 2.6× is 2041/776 — a projection-over-projection ratio: 2041 appears in no primary artifact, and 776 is a stage-share reconstruction using two factors that were never benchmarked. Presenting the ratio that drives the entire track re-ordering as a Phase-0 measurement is the class of error the previous LTS projection made.

**Trigger:** Any reader or future decision treating the comm-primary ordering as measurement-backed.

**Actual behavior:** Row reads "comm 60 % of wall, 2.6× the payoff, and bitwise-safe" under "Phase 0 measured".

**Expected behavior:** Measured: exposed Waitall = 52.5–52.6 % of the face-cache-OFF LTS wall. Projected: everything else.

**Suggested fix:** Rename the column to "Phase 0 measured / derived" and change the row to: "**measured:** exposed Waitall = 52.5 % of the LTS wall (C §3) / 52.6 % (A Leg 5); **projected:** comm recovery 1941–2191 vs kernel recovery ~710–776, i.e. 2.5–3.1× — *both sides are projections*." Add a footnote to the budget table's share column with the same distinction.

**Evidence:** C §3 (2341/4356 = 52.5 %); A Leg 5 (52.6 %); grep of `kernel_dev/`, `lts_dev/RESULTS_p5`, `comm_dev/` — "2041" and "776" appear only in this plan and its superseded v1.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* the third skeptic added a materially useful fact — A Leg 5 ran `CFG_RATE2` **with** `COMMON` (which includes `--face-cache`, sbatch:88), i.e. it *is* the face-cache-ON LTS configuration, and it measured 52.6 %, unchanged from face-cache-OFF; under the plan's own composition the share should have moved to 60 %. That is direct evidence against the 60 % and reinforces R-001/R-002. Still LOW because the fix is a column rename plus a footnote.

---

### [R-020] [LOW] [PLAN_v2:24, :84] — "Fixing the waiting recovers 2041 s/sim-s" contradicts the plan's own 24 % skew cap by 262 s/sim-s

**Category:** ARITHMETIC (lens = traceability)

**Description:** Same inconsistency as R-017, reached from the traceability side, plus the observation that the 300 residual's source — D's "merge+overlap → ~150–400" row — is a pure linear-in-round-count projection (2341×64/125 = 1199 ≈ D's 1200; 2341×32/125 = 599 ≈ D's 600) that assumes **zero** skew, directly contradicting D's own Summary risk statement. No post-vs-complete wire/skew decomposition exists.

**Trigger:** Any comparison of Track A against Track B payoff, and the whole "2.6× backwards" conclusion.

**Actual behavior:** Track A payoff stated as a point value 2041 s/sim-s.

**Expected behavior:** ≤ 2341 × (1 − 0.24) = 1779 under the plan's own cap; D's zero-skew model gives 1941–2191. The two are inconsistent and neither is measured.

**Suggested fix:** State Track A's payoff as a band with both bounds named and their models identified, and note that D's payoff rows are linear-in-round-count with zero skew, which its own Summary risk contradicts.

**Evidence:** PLAN_v2:37 vs :24; A Leg 5 skew proxy 24 %; D §Payoff arithmetic rows = 2341×64/125 and 2341×32/125 vs D §Summary.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* as R-017 — no gate keys off 2041; the ratio at the reviewer's own 1779 is still 2.29×, so the ordering is unchanged; and D:100-102 already pre-registers the escalation branch (fault-weighted partition) if the skew proves binding. The dissent also flagged a *separate* real issue: under a hard 562 floor, Track A's ≤15–20 %-of-wall gate is arithmetically unreachable — see R-027.

---

### [R-021] [LOW] [PLAN_v2:5, :143-144] — "Evidence base (all measured, none assumed)" and "the last unmeasured input to this plan" are both false as written

**Category:** UNSUPPORTED_CLAIM (lens = traceability)

**Description:** At least seven load-bearing inputs remain unmeasured: (1) the face-interior 1.5× and volume 2× factors inside 776 — no face or volume bench exists on any machine; (2) the face-cache gain on the LTS leg (1538); (3) the comm payoff 2341→150–400 (projection); (4) whether the Leg-4 bench ratio transfers to the in-solver predictor (no in-solver A/B); (5) rupture-phase comm skew (Leg 5 stopped at t=0.5, pre-nucleation); (6) the ~125 s/sim-s SeisSol denominator; (7) B §2's entire roofline, which rests on an **assumed** 2 GB/s/core bandwidth (no STREAM) and a hard-coded 71,280 FLOP/elem constant (no counter).

**Trigger:** A reader trusting the banner and skipping traceability checks on downstream numbers.

**Actual behavior:** ":5 **Evidence base (all measured, none assumed)**"; ":143-144 … is the last unmeasured input to this plan."

**Expected behavior:** An explicit unmeasured-inputs register.

**Suggested fix:** Change :5 to "**Evidence base:** …" (drop the parenthetical) and add a "## Unmeasured inputs (must be retired before the corresponding gate is claimed)" section listing the seven items above. Replace :143-144 with "… that retires the GTS-proxy caveat above; six further unmeasured inputs remain (see §Unmeasured inputs)."

**Evidence:** B §2 heading ("2 GB/s/core, conservative" — assumed, no STREAM anywhere in A–D); `bench_ader_kernel_variants.cpp:66-70`; C §7 open items; A "Projections from the measured split".

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* four of the seven are disclosed elsewhere (A:111 labels the face/volume factors "Projections"; PLAN_v2:76-82 caveats the GTS→LTS split; :34-37 flags the comm assumption; :124 flags the bench transfer), and item 7 underpins only the already-withdrawn ≤40 µs target. The sharper defect is :143-144's falsifiable "last unmeasured input" claim.

---

### [R-022] [LOW] [PLAN_v2:29-32] — The Expected-outcome band's low end (760) assumes zero residual exposed comm

**Category:** ARITHMETIC (lens = direction)

**Description:** 1538 − 776 = 762, so the stated band is 762 + (0 to 300). The 0 low end contradicts (a) the plan's own 2041, which embeds a 300 residual; (b) D's merge+overlap residual floor of 150–400; and (c) the plan's own :37 that ~562 s/sim-s of skew is untouchable. Under the plan's own skew cap the floor is 762 + 562 = 1328 s/sim-s = 5.6× behind SeisSol as-run.

**Trigger:** Composing the two tracks with the residuals the plan itself states.

**Actual behavior:** "~760–1060 s/sim-s = 3.2–4.4× behind SeisSol as-run".

**Expected behavior:** 762 + [150..566] = 912–1328 s/sim-s = 3.8–5.6×, with the skew-capped case as the honest floor.

**Suggested fix:** Replace :29-32 with a band whose low end uses D's own 150 floor, and state the skew-capped case explicitly: "…~910–1330 s/sim-s = 3.8–5.6× behind SeisSol as-run (239), from 16× today. A zero-residual composition (the 760 previously quoted) is not consistent with either the companion plan's floor or this plan's own skew cap. Against SeisSol's measured compute-only floor of 111 s/sim-s the same band is 8.2–12.0×; no async-IO SeisSol run exists, so no intermediate denominator is quoted."

**Evidence:** 1538−776 = 762; D §Payoff arithmetic "~150–400"; PLAN_v2:37 ⇒ 0.24×2341 = 562; C §Summary (239 as-run, 111 compute-only).

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* holding skew at 562 **absolute** while deleting 776 of the compute it disperses is itself a cross-context transfer; scaled (562 × 762/1538 = 278) the floor is 1040 = 4.35×, inside the plan's own band. The residual defect is only the zero low end.

---

### [R-023] [LOW] [PLAN_v2:5, :24] — "all measured, none assumed" is false: four of the numbers producing both headline payoffs are assumptions

**Category:** UNSUPPORTED_CLAIM (lens = direction)

**Description:** (1) face-interior 1.5× and (2) volume 2×, both inside the 776, have **no benchmark on any machine** — the bench models the predictor only. (3) predictor 4× is a chosen point inside a contended synthetic range whose control is invalid: three of four variants ran *faster* at 256-copy "full occupancy" than in the "single-core uncontended" baseline (B 12.5 vs 14.08, C 21 vs 40.67, D 8.1 vs 11.78 µs/elem), which is physically impossible and means the contention framing has no clean baseline. (4) the 300 s/sim-s residual inside 2041 is a design target.

**Trigger:** A reader treating 2041 and 776 as traceable to A/B/C/D.

**Actual behavior:** ":5 **Evidence base (all measured, none assumed):**".

**Expected behavior:** An explicit assumed-inputs list, including (v) the 1.31× face-cache factor applied to the LTS compute budget.

**Suggested fix:** Replace :5's parenthetical with an "Assumed (not measured) inputs carried by the payoff numbers" list naming (i) face 1.5×, (ii) volume 2×, (iii) predictor 4× with the invalid uncontended control spelled out, (iv) the 300 s/sim-s residual, (v) the 1.31× GTS→LTS transfer.

**Evidence:** A Leg 4 (single-core A 32.98 / B 14.08 / C 40.67 / D 11.78 vs contended B 12.5–14.4, C 21–55, D 8.1–13.8); `bench_ader_kernel_variants.cpp` header (predictor replica only); reconstruction of 776 requiring 1.5× and 2×.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* the Phase-2 GO rests on *within-condition* contended comparisons (A 79 vs B 12.5–14.4, vs D 8.1–13.8), which the single-core anomaly does not touch; the published 2.26× low end even uses the uncontended A=33 as numerator, so the range is conservative at its low end. The sharper instance is :144's "the last unmeasured input", not the :5 banner.

---

### [R-024] [LOW] [PLAN_v2:12-22, :14] — The whole budget rests on an unflagged GTS→LTS transfer of the 1.31× that the only LTS measurement contradicts

**Category:** ASSUMPTION (lens = direction)

**Description:** Same defect as R-001/R-006/R-016 reached from the direction lens, with the consequence quantified: if the true LTS factor is 1.0, the budget is 2015 + 2341 = 4356, the comm share is 53.7 % (not 60 %), and 776 rescales to 1017 — moving the headline ratio from 2.6× to 2.0×, or 1.75× once the skew cap is applied. :14 presents the 1.31× as already banked ("That already includes the free 1.31×").

**Trigger:** Face-cache providing no LTS-leg gain, as the only LTS datapoint indicates and as R-002 shows structurally.

**Actual behavior:** "3879 seconds … That already includes the free 1.31×"; budget 2341/1538 = 60/40.

**Expected behavior:** The 1.31× flagged as a GTS→LTS transfer with the contradicting LTS datapoint, and a sensitivity band on the split.

**Suggested fix:** Append to the Summary table: "*The 1538 compute term is **derived, not measured**: 2015 (C §3, LTS leg, face-cache OFF) ÷ 1.31, the activation gain measured on the **GTS** leg. If the LTS factor is 1.0 the budget is 2015 + 2341 = 4356 (comm 53.7 %) and Track B's 776 rescales to 1017, moving the A:B ratio from 2.6× to 2.0×. **Publishing the Leg-5 wall and re-running one face-cache-OFF/ON LTS pair retires this and must precede any final ordering claim.**"

**Evidence:** A §Composed end state (1538 = 2015/1.31); A Leg 5 (⇒ 4399 s/sim-s); C §Summary (LTS 4356 face-cache OFF).

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* kernels overtake comm only if compute > 3530 s/sim-s, i.e. 1.75× above the measured 2015 upper bound — so no value of the LTS face-cache factor in [1.0, 1.31] flips the ordering; and :140-144 already schedules the retiring measurement before any gate. Retained at LOW here because R-001/R-002 carry the same defect at CRITICAL with the source-level mechanism attached.

---

### [R-025] [LOW] [PLAN_v2:29-32, :155-163] — The "honest ceiling" and the "6–8.5× against async-IO SeisSol" rest on a denominator that exists in no artifact and on an accuracy match that was never measured

**Category:** UNSUPPORTED_CLAIM (lens = direction)

**Description:** (1) "6–8.5× against its faster async-IO configuration" back-implies a SeisSol denominator of ~125–127 s/sim-s. C reports only 239 (as-run) and 111 (compute-phase floor) and lists the async-IO re-run as an **open** item; the ~125 originates in B §3 as a column header. Against the measured 111, the plan's own band is 6.8–9.5×. (2) "we are ~3–4× behind SeisSol at matched accuracy" has **no accuracy measurement** behind it — C §4 and §7 list the station-level MFEM-vs-SeisSol waveform diff and the SCEC reference comparison as open. Only order (p3/o4) and dt were matched. This matters directionally: if the true comparator is 111 and neither track closes to it, the far-field order drop the plan relegates to "a separate proposal" becomes the dominant remaining lever.

**Trigger:** Any reader sizing the program against SeisSol's real production configuration.

**Actual behavior:** "(6–8.5× against its faster async-IO configuration)"; "~3–4× behind SeisSol at matched accuracy".

**Expected behavior:** Only the two measured denominators (239, 111) used, and "matched order/dt" rather than "matched accuracy".

**Suggested fix:** Replace the parenthetical on :30 with "(**8.2–12.0× against SeisSol's measured compute-only 111 s/sim-s**; no async-IO SeisSol run exists — C §7 lists that re-run as open)". Replace :156's clause with "…**at matched ORDER AND TIMESTEP** — no accuracy comparison exists (C §4, §7), so no matched-accuracy claim is admissible yet." Amend :162 to note that if the operative comparator is 111 rather than the IO-inflated 239, neither track closes the gap and the order drop becomes the dominant remaining lever.

**Evidence:** C §Summary (239 as-run, 111 compute-only); C §7 open items (async-IO SeisSol re-run NOT performed; waveform diff and SCEC comparison open); B §3 uses "~125" with no source; 760/111 = 6.85, 1060/111 = 9.55.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* E R-404 already adjudicated the 125 as a deliberate "bracket 111–239" with the async re-run scheduled, and B writes it with a tilde; the async figure appears once, parenthetically, and gates nothing. The dissent also surfaced a **separate real regression**: v2's restructure dropped R-404's scheduled async-IO SeisSol re-run, which had a home in v1's Phase 4 and now appears in no phase.

---

### [R-026] [LOW] [PLAN_v2:24, :29-32 vs :36-37] — The 24 % irreducible-skew cap makes both the 2041 comm payoff and the 760–1060 outcome band arithmetically impossible

**Category:** ARITHMETIC (lens = composition)

**Description:** 0.24 × 2341 = 562 s/sim-s is irreducible by the plan's own statement, so Track A's maximum recovery is 1779, not 2041 — the 2041 removes 262 s/sim-s of the very skew the same paragraph declares untouchable. Propagating: best composed = 1538 − 776 + 562 = 1324 s/sim-s = 5.54× SeisSol as-run, outside the stated 3.2–4.4× band at both ends.

**Trigger:** Reading the Summary payoff table together with the Biggest-risk paragraph on the same page.

**Actual behavior:** Track A payoff 2041; expected outcome 760–1060 (3.2–4.4×).

**Expected behavior:** Track A payoff ≤ 1779; composed floor stated explicitly, and the skew attacked by a partition/imbalance work item.

**Suggested fix:** State both bounds and name the fallback: "Reaching the advertised 3–4× band additionally requires removing the skew floor, which needs a fault-weighted repartition work item that is currently in neither track" (cf. R-004, which shows the plan closed the only measured remedy).

**Evidence:** PLAN_v2:36-37 × the 2341 basis (C §3); 0.24×2341 = 561.8; 1538−776+562 = 1324; 1324/239 = 5.54.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* C §4 explicitly says "`Waitall` Max/Avg 1.72 = comm skew, **not** element imbalance" with GTS element balance at Max/Avg 1.0000, so the 562 is not a measured load-imbalance floor; and holding it absolute while removing 776 of compute is itself a context transfer (scaled: 1040 = 4.35×, inside the band). D:31-33 and :99-102 already pre-register the fault-weighted-partition escalation. The surviving kernel is the band's low end (see R-018/R-022).

---

### [R-027] [LOW] [PLAN_v2:97 vs :20, :151-153] — Track A's gate ("exposed wait ≤15–20 % of wall") does not produce the 150–400 used in the composition, and its anchor contradicts the plan's own 60 %

**Category:** DEVIATION (lens = composition)

**Description:** The gate is a fraction of wall, so its absolute meaning depends on which compute term is in the denominator, and the Re-baseline rule says Track A is gated with Track B's flags OFF, i.e. at compute = 1538. Solving E/(1538+E) = 0.20 gives E = 384.5 s/sim-s; at 0.15, E = 271.4. So a **passing** Track A leaves 271–385 s/sim-s exposed, and the composed total is 762 + 271…385 = 1033–1147 = 4.3–4.8× — outside the stated 3.2–4.4× band. Under the alternative reading (gate measured post-Track-B at compute 762) the gate admits E = 134–190, i.e. 3.75–3.98× — still not the 760 low end under any reading. Separately, the gate text says "(from 53 %)", which is 2341/4356 on the face-cache-OFF wall, while the plan's headline is 60 % (2341/3879) — the two anchors differ by 7 pp and the gate has never been re-derived on the 1538 base.

**Trigger:** Track A passes its gate at the loose (20 %) end while Track B hits 776 — the plan's own definition of both tracks succeeding.

**Actual behavior:** Composition uses a residual of 0–300; the gate admits 385.

**Expected behavior:** An absolute gate consistent with the composition arithmetic, stated on the same compute base as the budget.

**Suggested fix:** Replace :97 with: "**Gate for the whole track:** exposed wait ≤ **300 s/sim-s absolute**, measured on the LTS leg with Track B's flags OFF. The percentage form is NOT the gate — it moves with the compute denominator and admits 385 s/sim-s at the 1538 base, which does not reach the composed target."

**Evidence:** PLAN_v2:20 (2341 = 60 % of 3879) vs :97 ("from 53 %" = 2341/4356); arithmetic 0.20×(1538+E) = E ⇒ E = 384.5; 762 + 384.5 = 1146.5; /239 = 4.80.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* the third skeptic showed the "53 %" anchor is the **most recent FINAL measured** exposed-wait share (A Leg 5, 52.6 %, annotated "consistent" with Phase-5's 52.5 %), not a stale one — the 60 % is the *derived* figure — so that half of the finding is inverted; and D:101-102 declares the acceptance configuration explicitly ("on the Phase-5 config"), where the gate admits E ≤ 504. The surviving defect is the gate/target coherence gap (a passing gate does not reach the published band) plus the two mismatched anchors in one document.

---

### [R-028] [LOW] [PLAN_v2:30] — The "6–8.5× against its faster async-IO configuration" claim uses a 125 s/sim-s SeisSol denominator that exists in no primary artifact

**Category:** UNSUPPORTED_CLAIM (lens = composition)

**Description:** 760/125 = 6.08 and 1060/125 = 8.48, so the parenthetical is computed against ~125 s/sim-s. C reports exactly two SeisSol rates — 239 as-run and 111 compute-phase-only — and lists the async-IO re-run as an open item in §7. Against C's actual compute floor of 111 the same band is 6.8–9.5×. The 125 is imported from B §3, where it is a column header.

**Trigger:** Any reader treating the parenthetical as a measured comparison.

**Actual behavior:** "(6–8.5× against its faster async-IO configuration)".

**Expected behavior:** Either the measured 111 s/sim-s compute floor with its basis stated, or the claim removed pending the open async-IO run.

**Suggested fix:** Replace the parenthetical with "(**6.8–9.5×** against SeisSol's measured **111 s/sim-s compute-phase floor**; there is **no measured async-IO SeisSol run** — RESULTS_p5 §7 still lists it as open, and the ~125 s/sim-s figure used in `EVIDENCE_microbench` §3 is traceable to no measurement)".

**Evidence:** C §Summary (239 as-run, 111 compute-only) and §7 open items; 760/125 = 6.08, 1060/125 = 8.48.

**Verification:** 2/3 skeptics failed to refute. *Dissent recorded:* B:53 does carry the figure as "~125" (tilde) and it interpolates between two measured endpoints (111 compute + 128 blocking IO = 239), so it is not sourceless; the plan merely dropped the tilde and the basis. Editorial precision fix.

---

## Summary

- **Critical issues:** 5
- **Moderate issues:** 10
- **Low issues:** 13
- **Plan compliance:** PARTIAL — the phase structure, per-phase gates and risk register are well formed and the two-track split is coherent, but the headline budget, the "all measured, none assumed" evidence claim, the exchange inventory, and the "prerequisite DONE" instrumentation claim are each contradicted by the source or by the primary artifacts.
- **Verdict:** **FAIL** — re-issue required before either track starts. Four things must land first: (1) the headline budget must be re-derived on the measured 4356 s/sim-s LTS wall — **R-002** shows in source that the 1.31× cannot apply to the LTS corrector, so 3879 / 60 % / "16×" / 2.6× are wrong by construction; (2) the wiggle/cluster-cap config sweep (**R-003**) must run, because a one-line deck change plausibly outweighs the entire Track B program and is currently invisible to both plans; (3) the four missing LTS Caliper scopes plus an LTS step region (**R-008**) must exist, since without them the promised "re-measured LTS stage split" cannot report 4 of its 6 stages; (4) one FLOP-counter leg (**R-005**) must decide the memory-bound-vs-more-FLOPs question that Track B's mechanism choice presumes. The plan's **directional** conclusion — communication before kernels — survives every finding and does not need re-arguing; it should simply be restated at its true magnitude (≈2.0×, not 2.6×), with the fault-balance item (**R-004**) reopened inside Track A rather than closed.

## What the evidence DOES support

These load-bearing claims are properly measured and survive the audit:

1. **The LTS leg's exposed communication is real, large, and LTS-specific.** MPI_Waitall = 52.5 % of wall (C §3, job 52344266, t = 0..2.0, 256 ranks), independently reconfirmed at 52.6 % on a different job (A Leg 5, job 52379131). The GTS leg on the same mesh shows 0.00002 %. 17.3 M waits and 585 M messages reconcile exactly with the code-derived 125 exchange rounds/sync at Nc=6 (verified against `lts_stepper.hpp:140-165` and the two `++n_ghost_exchanges_` sites at `wave_operator.inl:2706, :2730`), to 0.25 %. **Comm is the largest single line item and the merge lever's denominator is sound.** (The unreconciled tension — GTS runs 2.1× more halo rounds per sim-second yet waits ~10⁶× less — is a real open question about *mechanism*, but it does not touch the size of the LTS wait.)
2. **The GTS-leg stage split is a genuine Caliper measurement** on the pinned benchmark: predictor 56.23 % (child `ApplySpatialDerivative` 27.14 %), corrector 34.40 % (face interior 12.96 %, face shared 10.56 %, volume 7.88 %), fault 0.65 %. It closes only to 91.28 % (R-010) and it is GTS not LTS (the plan flags this), but the shares themselves are measured, not inferred. **The v1→v2 inversion of face-vs-predictor priority is correctly grounded in this.**
3. **The matched-configuration cross-code comparison is rigorous.** Same 2,464,689-tet mesh, md5-verified on both sides; both codes independently computed dt_min = 3.66128e-4 s to 6 digits at Courant 0.5; SeisSol's cluster histogram matches MFEM's raw λ=1 report exactly. **The 17.1× per-update gap is real in magnitude**, even if its normalisation is debatable (16.0× charging SeisSol the whole node; 13.0× at B0's face-cache-ON figure).
4. **The withdrawal of the ≤40 µs·core / ≥6× stretch target is correct and well-argued** — it follows from the measured stage split by Amdahl arithmetic that reproduces exactly, and every skeptic pass strengthened rather than weakened it.
5. **`BatchedLinAlg` is genuinely out.** Variant C is the worst variant on both machines (0.46× locally, 0.69–3.13× contended on Rome), measured twice independently, including a re-run during this audit that reproduced B §1 within ~6 % and showed the local verdict is not an `-O2`/`-O3` artifact.
6. **Some kernel headroom exists in the predictor.** The perf self-time attribution (`ApplySpatialDerivative` 23.8 %, `Vector::Add` 18.9 %, `Vector::operator=` 9.3 %) is a real rank-0 measurement, and the whole-vector AXPY/zero traffic it points at is confirmed present in the source (`wave_operator.inl:2031-2055`, inside the LTS predictor). **B1 is the best-supported Track B phase** — its lever is measured, its target routine is instrumented, and it is the cheapest.
7. **The companion comm plan's static analysis is accurate.** Every one of D's exchange-site claims that was checkable against source verified correct (63 I-rounds + 62 forecast = 125; `PrepareSeamCoarseForecast` retain-only; substep fault exchange contributes 0 on the LTS path; the 9-sequential shared corrector is GTS-only; MFEM's ldof tables are public so no core edit is needed). Its line citations are stale by exactly one commit (`c547ac8` inserted 2 lines), it omits one site (`ComputeSharedFaceFluxRHS:4847`) and over-counts one collective row, but its structural picture is sound.

## Unreviewed Areas

- **Runtime re-execution of any Expanse leg.** Per the standing project rule, no cluster job was submitted or inspected live during this audit; every cluster number was audited against the committed artifacts and decks only. The bound this review computes for the abandoned Leg 3 (≤ 225.9 µs·core) is arithmetic on A's own published job wall, not a new run.
- **The raw B0 profiles.** `b0_gts.perf.data` was retained but `perf report` needs a perf-enabled host, so the call tree behind the self-time list is unavailable; the Rome bench log (`$OUT_ROOT/bench.log`) is not committed, so the Leg-4 single-core and 256-copy numbers cannot be independently re-derived from the repository.
- **The Rome bench toolchain.** The compiler, flags and BLAS used to build `bench_ader_kernel_variants` on Expanse are recorded nowhere (the link line is scraped from `make -n` at run time); this could not be reconstructed.
- **Numerical correctness of the proposed mechanisms.** Whether the tiled-fused predictor (B2) can preserve the `dk_retain` bit-identity contract, and whether A1's merge is genuinely superposition-safe on a rich seam, were assessed only from the plans' descriptions and the existing source structure — no prototype exists to review.
- **The far-field order drop.** Explicitly out of scope per the plan; no p1/p2 cost data exists anywhere in A–D, so the "49 % cells p1" row in B §3 and the plan's parity route could not be evaluated on evidence. R-005 argues the *decision to defer* may need re-opening; that is not a judgement on the proposal's merits.
- **Production-mesh behaviour.** Everything here is TPV104-200m at 256 ranks. The SAFS v4_0_0 meshes (34.9–52.3 M tets) — the actual target of the branch this plan sits on — have no performance measurement in any artifact (see R-015).
