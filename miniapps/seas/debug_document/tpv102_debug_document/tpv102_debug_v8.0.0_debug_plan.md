# TPV102 v8.0.0 — Systematic Debugging Plan (Revision 7)

## Revision History
- **Rev 1–6 (deleted):** see prior REVIEW.md rounds 1–7.
- **Rev 7 (current):** addresses round-7 findings R-701 through R-705.
  Summary:
  - **R-701 bulk-DOF propagation monitors:** Phase 2 Step 5 now flags
    2–3 **BULK** DOFs at 1 / 2 / 3 km up-dip from the hypocenter (at
    `y ≈ 300 m` off the fault, into bulk).  These are the
    propagation-path witnesses — if waves are injected correctly at
    the hypocenter but die before reaching off-hypo fault stations,
    only the bulk monitors reveal it.  Each monitor is MPI-verified
    to live on exactly one rank and prints
    `max |Q[VX]| / max |Q[SXY]|` every 100 RK4 steps.  Phase 2's
    Decision Rule gains a "wave dies en route" branch (`monitor at
    1 km grows, 2 / 3 km don't` → volume-integration / mass-inverse
    bug).
  - **R-702 DIAG step throttle:** introduces a global
    `g_seas_diag_step_gate` set by the driver at each RK4 step
    (every 100 steps for routine telemetry, every 10 steps during
    the breakaway window).  `[EVAL]`, `[FLUX]`, `[BULK]` blocks
    AND-test `data.diag_print && g_seas_diag_step_gate` so log
    volume at 400 ranks / tfinal=1.2s drops from ~40 GB to ~1 GB.
  - **R-703 Tandem reference in Round 1:** Tandem moves from
    Escalation step 5 ("last resort") to a new **Round 1 Slot B**
    parallel-track ground-truth reference.  Diff Tandem's off-hypo
    station V₂(t) curves against our Slot A MFEM+DIAG output;
    branches the diagnosis cleanly (Tandem slips 9/9 + we don't
    → our code; both fail → mesh; neither reproduces → version
    mismatch).
  - **R-704 off-hypo cross-rank selection:** the previous off-hypo
    choice `flt_0_3` (4.5 km up-dip) may land on the same METIS
    partition as the hypocenter at 400 ranks, defeating the
    cross-rank propagation test.  Phase 2 Step 4 now picks
    `flt_0_3` for `tfinal=1.2 s` (arrival 1.30 s — marginal but
    accepted) AND `flt_0_12` (4.5 km **down**-dip; near-guaranteed
    different partition) for `tfinal=1.5 s` runs.  MPI_Allreduce
    verifies hypo and off-hypo are on different ranks.
  - **R-705 CFL "probe" → "log":** Phase 0's CFL section renamed
    from "probe" to "log"; gating language removed (superseded by
    R-601 retraction); clarified as informational-only.
  - **Unreviewed Areas (round 7):** `GetFaceElementTransformations`
    semantics for shared faces post-ExchangeFaceNbrData; degenerate
    h_min on anomalously thin rank elements; `diag_print` stability
    across scatter/gather; step-counter visibility in Evaluate.
- **Rev 6 (superseded):** addresses round-6 findings R-601 through R-608.
  Summary:
  - **R-601 CFL arithmetic correction:** the round-6 critique computed
    `dt_cfl = cfl · h_min / cp` with `h_min = 200 m` and concluded
    observed `dt ≈ 15 · dt_cfl`. **That arithmetic is wrong.**
    `WaveOperator::ComputeMaxDt` uses `h_min` = tet **inscribed
    diameter** = `6·V / A_total` (wave_operator.inl:55-96), which
    for a regular 200 m-edge tet is `200/√6 ≈ 81.6 m`, not 200 m.
    Correct `dt_cfl ≈ 0.0556 · 81.6 / 6000 ≈ 7.6e-4 s` on a regular
    tet, and somewhat smaller on irregular tets.  Observed
    `dt ≈ 2.85e-4 s ≈ 0.38 · dt_cfl` — **subscribing to CFL, not
    violating it**.  The CFL hypothesis is therefore ELIMINATED by
    arithmetic, not confirmed.  Phase 0 keeps a lightweight
    `dt_actual / dt_cfl` ratio probe (permanent sanity check), but
    no Phase 0A investigation is required.  The recorded-facts
    table and Phase 0 text are corrected accordingly.
  - **R-602 hypo DOF distance filter:** Phase 2 Step 4 adds a 500 m
    distance gate on the hypocenter DOF tag, symmetric to the
    existing off-hypo gate.  Without this, every rank flagged its
    local "nearest-to-hypo" DOF, producing N-rank × 2 spurious
    DIAG streams.
  - **R-603 Phase 3A propagation-corridor fixture:** Phase 3A adds
    Fixture B — a 5-tet corridor mesh (fault at leftmost face,
    probe at rightmost element) — so volume-integration bugs that
    kill bulk propagation across multiple elements become testable.
    Fixture A (existing 2-tet face-flux accumulation) remains.
  - **R-604 O(1) fault-face lookup:** Phase 1A Action 1a replaces
    the linear scan of `fault_interior_faces_` with a
    `std::unordered_set<int> fault_interior_face_set_` for O(1)
    membership.  Avoids ~10 min/job overhead at 50-rank Frontera.
  - **R-605 off-hypo arrival timing:** S-wave arrival time
    `flt_0_3` ≈ 4500/3464 ≈ 1.3 s.  At `tfinal=1.2 s` the front
    hasn't reached the station; expected-value table now notes
    this and recommends `tfinal=1.5 s` for off-hypo DIAG or an
    intermediate bulk-DOF monitor at z=-5500 m (arrival ≈ 0.58 s).
  - **R-606 banner-to-file:** Phase 1A Action 0 also writes the
    `[BUILD]` banner to `build_info.txt` so interleaved Frontera
    stderr cannot drop it.
  - **R-607 Phase 3A threshold parameterization:** Fixtures A & B
    now derive their `max_k` thresholds from the fixture's own
    `h_min` so tests don't break if fixture edge length changes.
  - **R-608 Finding Index appendix:** central map of all R-IDs
    (R-001 → R-608) to review round and plan revision addressing
    them.
- **Rev 5 (superseded):** addressed round-5 findings R-501 through R-508.
  Summary:
  - **R-501 Makefile + startup banner:** add `SEAS_EXTRA_CPPFLAGS` to
    the SEAS Makefile (bypasses MFEM's `?=` precedence) and a
    startup banner in the driver that prints whether
    `SEAS_DIAG_FAULT_FLUX` / `SEAS_DIAG_GHOST_EXCHANGE` are compiled
    in. Phase 2 grep this banner before analysing DIAG output.
  - **R-502 `tfinal`:** every diagnostic sbatch uses **`tfinal = 1.2 s`**
    (post-breakaway), never 0.3 s. The `0.3 s` wall-clock promise is
    abandoned; Tier 2 is now a 100-rank/2-node dev job at ~1.2 hr
    wall-clock.
  - **R-503 interior-branch silent-weld guard:** Phase 1A adds
    **Action 1a** — `MFEM_VERIFY` that `is_fault == (f ∈ fault_interior_faces_)`
    cross-check in `ComputeFaceFluxRHS` so misclassified fault faces
    abort rather than silently weld.
  - **R-504 off-hypocenter DIAG witness:** Phase 2 Step 5 flags TWO
    DOFs for DIAG — hypocenter AND `flt_0_3` station (4.5 km up-dip).
    The off-hypo DOF is the propagation witness.
  - **R-505 wall-clock + rank routing:** Tier 2 defaults to
    100-rank/2-node/tfinal=1.2s (~1.2 hr). 50-rank only used for
    `tfinal ≤ 0.7 s` init/partition checks, never for bug-regime
    diagnostics.
  - **R-506 success criterion:** count 0 / 1-3 / 4-9 out of the 9
    standard stations with `V2 > 1e-6` at `t=1.2s`. Distinguishes
    total failure, partial propagation, and healthy.
  - **R-507 dead-code removal:** `MFEM_ABORT`s in Phase 1A Actions 1,
    1a, 1b no longer leave `flux_.Interior(...)` as unreachable
    tail. Replaced with `// Control never reaches here.`
  - **R-508 inline unchanged sections:** Phase 1A, Phase 2 local
    track, Phase 3A/B/C, Phase 4, Phase 5 are inlined in full
    rather than "UNCHANGED from Rev N" references.
  - **Unreviewed Areas:** (a) `dt` vs `h_min/cp` startup print so
    CFL-related bulk damping can be ruled out; (b) `TEST_R007_AllInteriorFaultFacesMapped`
    must pass BEFORE the R-503 guard is added (else the guard
    fires on correct code); (c) explicit cross-check between
    `VerifySharedFaultDOFDataConsistency`'s face list and
    `ComputeSharedFaceFluxRHS`'s `sf_fault` classification at
    step 0.

## Purpose

Root-cause the v8.0.0 critical bug (rupture pinned at hypocenter;
bulk `max‖Q‖_∞ ≈ 1.66e-8` while V_max = 7.7 m/s; off-hypocenter
fault stations never perturbed) via a sequential, decision-driven
diagnostic process.

## Ground Rules

1. **One hypothesis at a time.**
2. **Every phase has a DECISION RULE** determining the next phase.
3. **Each phase produces a written OUTCOME** under "Phase Log."
4. **Phase 0 baseline → Phase 1 pre-flight + bisection → Phase 1A
   cheap hardening (both branches, BOTH is_fault classifications) →
   Phase 2 focused instrumentation (hypo + off-hypo DOFs).**
5. **Tiered reproduction [R-301, R-307, R-505]:**
   - **Tier 1 (local ≤8 ranks, cheap, iterative):** try first.
   - **Tier 2 (Frontera dev queue, 100-rank / 2 nodes, tfinal=1.2s,
     ~1.2 hr wall-clock):** if Tier 1 fails.  100-rank not 50-rank,
     because 50-rank × tfinal=1.2s exceeds the 2-hr dev-queue limit.
     Bug signature REQUIRES tfinal ≥ 1.0s (post-breakaway) per R-502.
   - **Tier 3 (Frontera dev queue, 400-rank DIAG, tfinal=1.2s,
     ~1.5 hr wall-clock):** dispositive scale.
   - **Tier 4 (normal queue):** fallback if dev can't resolve.
6. **Every fix gets a unit test added BEFORE marking the phase complete.**
7. **All diagnostic code gated on `#ifdef SEAS_DIAG_FAULT_FLUX` or
   `#ifdef SEAS_DIAG_GHOST_EXCHANGE`.**  Flag propagation via the
   `SEAS_EXTRA_CPPFLAGS` Makefile variable (R-501), NOT via
   `CXXFLAGS+=` (MFEM's `?=` precedence silently drops this).
8. **Use `MFEM_VERIFY` for always-on invariants; `MFEM_ASSERT` only
   for debug-only pre/post-conditions.**
9. **Startup banner [R-501, R-606]:** the driver MUST print the
   active DIAG flag state (`ON`/`OFF`) after MPI init, to BOTH
   stderr and the dedicated file `build_info.txt` (rank-0, CWD).
   Every Phase 2+ diagnostic analysis begins with
   `grep "\[BUILD\]" build_info.txt` — the file is deterministic
   and resistant to stderr interleaving / truncation at Frontera
   scale.  Fall back to `grep "\[BUILD\]" <job>.out` only if
   `build_info.txt` is missing.
10. **At most ONE DOF per rank per DIAG flag.**  `MFEM_VERIFY(local_diag_count <= 1)`
    prevents log-size explosion if an accidental every-QP flag.
11. **DIAG step throttle [R-702]:** every DIAG block must AND-test
    the global `g_seas_diag_step_gate`.  The driver sets the gate
    to 1 every 100 RK4 steps routinely, every 10 steps in the
    breakaway window (`step_frac ∈ [0.75 · nsteps, 1.0 · nsteps]`,
    approximately `t ∈ [0.9 s, 1.2 s]`).  At 400 ranks / tfinal=1.2 s,
    this bounds DIAG log volume to ~1 GB aggregate.  Blocks affected:
    `[EVAL]` in `FaultFaceFlux::Evaluate`, `[FLUX]` in
    `ComputeFaceFluxRHS`, `[BULK]` in the driver RK4 loop.

## What We Already Know (recorded facts)

| Fact | Source |
|---|---|
| Fault V_max saturates at 7.699 m/s at hypocenter from t≈1.2s | Frontera job 7665297 log |
| All off-hypocenter fault stations at V_ini=1e-12, tau2=75 MPa | Station files |
| `max‖Q‖_∞` stays at O(1e-8) for 1.5s | qnorm watch |
| **Breakaway at t≈1.0s** (V_max crosses 0.1 m/s) | v7.1.0 check log |
| R-101 verifier: 405 pairs, 135 shared fault faces, max_rel_diff=6.66e-16 | Step-0 log |
| `FaultFaceFlux::Evaluate` formulas match SeisSol bit-for-bit | v8.0.0 SeisSol audit |
| TPV102 1000m mesh at 2 ranks: `shared=0` | v2_fix.md:128-129 |
| 50-rank / 200m / tfinal=1.5s = ~3 hr wall-clock (sbatch header) | jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch |
| 100-rank / 200m / tfinal=1.2s ≈ ~1.2 hr wall-clock (est. 2 nodes, 2× speedup) | extrapolation |
| MFEM Makefile uses `MFEM_CPPFLAGS ?= $(CPPFLAGS)`; config.mk overrides | Makefile:384-385 |
| `WaveOperator::ComputeMaxDt(cfl)` returns `cfl · h_min / cp` | wave_operator.inl:1288-1291 |
| `h_min` = tet **inscribed diameter** = `6·V / A_total`, NOT edge length | wave_operator.inl:55-96 |
| Regular tet inscribed diameter = `L_edge / √6 ≈ 0.408 · L_edge` | geometry |
| For "200 m" mesh (edge 200 m, regular): `h_min ≈ 81.6 m`, `dt_cfl ≈ 7.6e-4 s` at p=1 | derived |
| Observed dt on Frontera 400-rank 1.5 s run: `≈ 2.85e-4 s` (5250 steps) | sbatch header |
| Observed ratio `dt_actual / dt_cfl ≈ 0.38` — below CFL limit, NOT above | derived |
| Explicit RK4 Godunov DG linear-stability envelope ≈ `2 · dt_cfl` | feedback_explicit_cfl.md |
| Job 7664250 (`--dt 5e-3 s` ≈ `6.6 · dt_cfl`) blew up with geometric V_max amplification | feedback_explicit_cfl.md |
| DG CFL scales with `(2p+1)` vs CG `p` — DG is strictly stricter at matched order | Hesthaven-Warburton Ch. 4 |
| Driver formula: `cfl = cfl_factor / (3·(2·order+1))` (3D DG upwind) | tpv102_driver.cpp:277 |

## Debugging Phases

### Phase 0 — Freeze, Label, Baseline, Allocation, CFL Probe (0.5 hr)

**Goal:** establish a reproducible baseline AND confirm Frontera
allocation AND rule out CFL-related bulk damping as a spurious cause.

**Actions:**
- `git status` clean. Working commit: `c015048` (v7.0.0 fix).
- `git tag v8.0.0-bug-state` (local only).
- Capture baseline test state:
  ```bash
  cd miniapps/seas && conda activate mfem-dev
  make -j test 2>&1 | tee /tmp/phase0_baseline_tests.log
  mpirun -np 2 ./seas_test_r101_shared_fault \
      2>&1 | tee -a /tmp/phase0_baseline_tests.log
  mpirun -np 2 ./seas_test_parallel_wave_operator \
      2>&1 | tee -a /tmp/phase0_baseline_tests.log
  ```
- **[Unreviewed Area (b)] Run TEST_R007_AllInteriorFaultFacesMapped
  BEFORE adding the R-503 guard.**  If this test fails on current
  code, the `fault_interior_faces_` array has a ctor bug and the
  R-503 guard will fire spuriously. Must fix the ctor first.
- Frontera allocation check (ssh to login node):
  ```bash
  module load TACC
  /usr/local/etc/taccinfo || gbalance  # try both, log output
  # Confirm ≥ 500 SUs remaining.  Also: df -h $SCRATCH to confirm
  # ≥ 20 GB free for DIAG output.
  ```
- **CFL / dt sanity LOG [RENAMED per R-705; R-601 retracted]:**
  This is a permanent informational log, NOT a gating probe.  The
  CFL hypothesis was retracted in Rev 6 after R-601's arithmetic
  error (`ComputeMaxDt` uses tet inscribed diameter ≈ `L/√6`, not
  edge length).  The log captures `h_min`, `cp`, `dt_cfl`, and
  `dt_actual` for the Phase Log record and future regression
  tracking.  **Proceed to Phase 1 regardless of values** — the
  only "STOP" condition is `dt_actual > 2·dt_cfl`, which triggers
  the stability-envelope check (Phase 0 Decision Rule below), and
  which the driver itself logs as a WARNING at startup.  Run locally:
  ```bash
  ./seas_tpv102_driver --mesh tpv102/mesh/tpv102_1000m.msh \
     --order 1 --cfl 0.5 --tfinal 0.01 2>&1 | \
     grep -E "dt_cfl|dt:|dt \(override|Steps|CFL" \
     > /tmp/phase0_cfl_probe.log
  cat /tmp/phase0_cfl_probe.log
  ```
  From the output, extract the two numbers the driver prints:
  - `CFL: <cfl_eff>, dt_cfl = <dt_cfl> s`
  - `dt: <dt_selected> s` (or `dt (override, --dt): <dt_selected> s`).

  **Expected (correct) arithmetic** on `tpv102_1000m.msh` (regular
  tets, L_edge ≈ 1000 m):
  - `h_min ≈ 1000 / √6 ≈ 408 m` (tet inscribed diameter)
  - `cp ≈ 6000 m/s` (TPV102 material)
  - `cfl_eff = 0.5 / (3·3) = 0.0556` (cfl_factor=0.5, order=1)
  - `dt_cfl ≈ 0.0556 · 408 / 6000 ≈ 3.78e-3 s`
  - `nsteps(1 s) ≈ 265`

  **Expected arithmetic** on `tpv102_200m.msh` (regular tets,
  L_edge ≈ 200 m, `h_min ≈ 81.6 m`):
  - `dt_cfl ≈ 0.0556 · 81.6 / 6000 ≈ 7.56e-4 s`
  - `nsteps(1.5 s) ≈ 1985` (regular-tet lower bound)

  **Observed** on Frontera 400-rank 200 m 1.5 s run: `5250 steps` →
  `dt ≈ 2.85e-4 s`.  Ratio `dt_actual / dt_cfl ≈ 0.38`.  This is
  BELOW the CFL limit — **no violation**.  (The factor ~0.38
  between `dt_actual` and my back-of-envelope `7.56e-4` is
  consistent with mesh irregularity: actual min inscribed diameter
  in the METIS-partitioned 200 m mesh is smaller than for a
  perfectly-regular tet.)

  **Decision rule (permanent sanity check):**
  | `dt_actual / dt_cfl` | Interpretation |
  |---|---|
  | `< 0.5` | under-CFL; normal, proceed |
  | `0.5 – 1.0` | near the CFL line; note but proceed |
  | `1.0 – 2.0` | above CFL-1 but within the ~2× Godunov-DG envelope; log a WARNING, proceed with care |
  | `> 2.0` | **STOP** — per `feedback_explicit_cfl.md`, this exceeds the linear-stability envelope.  Investigate `--dt` override in the sbatch, or investigate `ComputeMaxDt` before any other phase |

  **If a `--dt X` flag is in the sbatch and `X > 2 · dt_cfl`**,
  that is the fix — remove it (see `feedback_explicit_cfl.md`).

  **DG vs CG reminder (informational):** the driver's formula
  `cfl_factor / (3·(2·order+1))` is the textbook DG-upwind /
  Godunov CFL estimate for 3D, linear in `(2p+1)`.  DG is strictly
  stricter than CG at matched polynomial order (typical ratio
  `1/(2p+1)`): CG p=1 allows `dt ≲ h/cp`; DG p=1 allows
  `dt ≲ h/(cp · 3)` — ~3× stricter.  At p≥3, some analyses use
  `(2p+1)²` — our driver's linear formula could become unsafe;
  separate concern, not v8.0.0 scope.

**Decision rule:**
- All baselines + allocation + `dt_actual / dt_cfl ≤ 2.0` clean → Phase 1.
- `dt_actual / dt_cfl > 2.0` → STOP; investigate `--dt` override or
  `ComputeMaxDt` before any fault-flux instrumentation (per R-601
  correction).  Record investigation under "Phase Log: Phase 0".

---

### Phase 1 — Partition Pre-flight + BISECTION (1.5 hr)

**Goal:** identify which fault-flux branch fails to propagate.

**Actions:**
1. **Rebuild with DIAG flag guarded by startup banner [R-501]:**
   first apply Phase 1A's Action 0 (Makefile change + banner) — see
   Phase 1A below.  Build both a non-DIAG binary and a DIAG binary.

2. **Partition pre-flight (10 min) [LOCAL, 1000m + 500m ONLY, R-304]:**
   ```bash
   PRE_FLIGHT_LOG=/tmp/phase1_preflight.log
   : > "$PRE_FLIGHT_LOG"
   for np in 2 4 8; do
      for m in tpv102_1000m tpv102_500m; do
         if [ ! -f "tpv102/mesh/${m}.msh" ]; then continue; fi
         echo "=== np=$np mesh=$m ===" >> "$PRE_FLIGHT_LOG"
         set -o pipefail
         mpirun -np $np ./seas_tpv102_driver \
            --mesh tpv102/mesh/${m}.msh --tfinal 0.001 \
            --output-dir /tmp/partition_check_${np}_${m} \
            --output-prefix test 2>&1 \
          | grep -E "Fault QPs|R-101 check|MFEM abort|\[BUILD\]" \
          | tee -a "$PRE_FLIGHT_LOG"
         EXIT_CODE=$?
         [ $EXIT_CODE -ne 0 ] && \
             echo "(driver exited $EXIT_CODE — OOM or mesh failure)" \
                 >> "$PRE_FLIGHT_LOG"
      done
   done
   # If no local combo yields shared>0, log ESCALATE directive.
   if ! grep -qE "shared: [1-9]|\[R-101 check\].*[1-9][0-9]* pairs matched" \
        "$PRE_FLIGHT_LOG"; then
      echo "" >> "$PRE_FLIGHT_LOG"
      echo "PRE_FLIGHT: no local (np≤8, {1000m,500m}) yields shared>0." \
           >> "$PRE_FLIGHT_LOG"
      echo "PRE_FLIGHT: ESCALATE to Phase 1 Escalation A-Frontera." \
           >> "$PRE_FLIGHT_LOG"
   fi
   cat "$PRE_FLIGHT_LOG"
   ```

3. **Run 1a — single rank (interior-fault branch only):**
   ```bash
   ./seas_tpv102_driver --mesh tpv102/mesh/${BISECT_MESH}.msh \
      --order 1 --cfl 0.5 --tfinal 1.2 \
      --output-dir /tmp/debug_v8_1rank \
      --output-prefix test --debug-qnorm
   ```
   **`tfinal = 1.2 s` not 0.3 s (R-502):** the bug signature requires
   post-breakaway (V_max > 0.1 m/s at t ≈ 1.0 s).

4. **Run 1b — $BISECT_NP ranks (both branches active):**
   ```bash
   mpirun -np $BISECT_NP ./seas_tpv102_driver \
      --mesh tpv102/mesh/${BISECT_MESH}.msh \
      --order 1 --cfl 0.5 --tfinal 1.2 \
      --output-dir /tmp/debug_v8_Nrank \
      --output-prefix test --debug-qnorm
   ```

5. **Success criterion (propagation scoring) [R-506]:** at `t = 1.2 s`,
   count how many of the 9 SCEC fault stations (flt_0_3, flt_0_7.5,
   flt_0_12, flt_9_7.5, flt_12_3, flt_12_12, flt_n9_7.5, flt_n12_3,
   flt_n12_12) show `V2 > 1e-6 m/s`:
   - **0 of 9:** total propagation failure (v8.0.0 signature).  Case A or B.
   - **1-3 of 9:** partial failure; record which stations slip,
     which don't.  The silent stations' nearest shared fault faces
     are the focus for Phase 2.  Case B-partial.
   - **4-9 of 9:** healthy propagation.  Case C (bug only at larger scale).
   - Supplementary: `max‖Q‖_∞ > 1 Pa` at `t = 1.0 s` (qnorm watch).

**Decision rule:**

| Case | Observation | Root cause | Next phase |
|---|---|---|---|
| A | 1-rank: 0 of 9 slip; N-rank: 0 of 9 | Interior branch OR Evaluate OR accumulation | Phase 1A + Phase 2 (interior + propagation instr.) |
| B | 1-rank propagates, N-rank 0 of 9 slip | Shared-fault branch | Phase 1A + Phase 2 (shared branch) |
| B-partial | N-rank 1-3 of 9 slip | Partial partition-dependent failure | Phase 1A + Phase 2 targeted at silent stations |
| C | N-rank 4-9 of 9 slip | ≥ (N+k)-rank only | Phase 1A + Phase 2B (starting with Frontera 100-rank) |
| D | Neither completes | Harness issue | Reduce mesh/tfinal |

**Phase 1 Escalation A-Frontera (PRIMARY when pre-flight yields no local `shared > 0`):**

1. **(5 min) Sanity-submit 50-rank init-only:**
   ```bash
   cp jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch \
      jobs/tpv102/tpv102_200m_p1_0.01s_50rank_INIT.sbatch
   # Edit: --tfinal 0.01; keep other flags; do NOT add DIAG yet.
   sbatch jobs/tpv102/tpv102_200m_p1_0.01s_50rank_INIT.sbatch
   # Wait <5 min.  Confirm [R-101 check] line shows "N pairs matched" with N > 0.
   # If N == 0, escalate to 100-rank for the DIAG run.
   ```

2. **(10 min) Create the DIAG sbatch at 100 ranks / 2 nodes [R-502, R-505]:**
   ```bash
   cp jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch \
      jobs/tpv102/tpv102_200m_p1_1.2s_100rank_DIAG.sbatch
   # Edit:
   #  - #SBATCH -n 100  (was 50)
   #  - #SBATCH -N 2    (was 1)
   #  - --tfinal 1.2    (was 1.5; enough for breakaway + saturation)
   #  - In the build step: replace CXXFLAGS+=... with
   #    `make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX"`
   #  - keep --debug-qnorm
   #  - redirect stderr: add `2> diag.err` after the ibrun invocation,
   #    OR use SLURM -e filename.err.
   ```

3. **(1.2 hr wall-clock) Submit:**
   ```bash
   sbatch jobs/tpv102/tpv102_200m_p1_1.2s_100rank_DIAG.sbatch
   ```

4. **(10 min) Return to Phase 1 decision table:**
   - BEFORE analysing DIAG output, `grep "\[BUILD\]" build_info.txt`
     (preferred) or `grep "\[BUILD\]"` the `.out` file (fallback)
     to confirm `SEAS_DIAG_FAULT_FLUX = ON`.  If OFF, the Makefile
     change didn't propagate — STOP and fix the build (R-501, R-606).
   - Score the propagation: count stations slipping (R-506).
   - Fall into Case A / B / C based on the table above.

**Phase 1 Escalation A-local (no Frontera access):**
- Build a new mesh with METIS-unfriendly partitioning (structured
  grid forcing cuts across the fault).  Document and flag for the
  standard test mesh inventory.  The 2-tet fixture is NOT
  sufficient for propagation testing (no off-hypo QPs).

**Phase 1 Slot B — Tandem reference baseline [NEW per R-703, concurrent with Slot A Frontera escalation]:**

Tandem (`/Users/chunhuizhao/projects/tandem/`) has a verified
production TPV102 implementation.  A Tandem run at matched
`(200 m, tfinal = 1.2 s)` gives a reference V₂(t) time-series at
every SCEC station for direct diff against our MFEM Slot A output.
This is ORDERS OF MAGNITUDE cheaper than iterating through
Phase 2 / 3 guesswork without ground truth.

**Actions:**

1. **(10 min) Identify or create a Tandem TPV102 sbatch script on
   Frontera.** Check `/Users/chunhuizhao/projects/tandem/examples/`
   for existing SCEC TPV102 inputs; if missing, compose from the
   Tandem documentation using TPV102's SCEC spec (BP5-like layout).

2. **(submit concurrent with Slot A) Submit the Tandem reference job.**
   Same allocation EAR20006, 8 nodes / ~400 ranks, development queue.
   Mirror our sbatch's mesh + tfinal + output cadence.  Wall-clock
   ~1.5 hr (Tandem is also explicit-DG with a Godunov flux).

3. **(10 min post-completion) Diff station outputs.** For each of
   the 9 SCEC stations, compare V₂(t) and `|slip_rate|(t)` from
   Tandem's output against our Slot A output.  Use a simple script:
   ```bash
   for sta in flt_0_3 flt_0_7.5 flt_0_12 flt_9_7.5 flt_12_3 \
              flt_12_12 flt_n9_7.5 flt_n12_3 flt_n12_12; do
      python3 compare_stations.py \
         --a our/results/${sta}.dat \
         --b tandem/results/${sta}.dat \
         --tol-rel 0.1 --label "$sta" \
         >> /tmp/tandem_mfem_diff.log
   done
   ```

**Decision rule after Round 1 (Slot A + Slot B both completed):**

| Slot A (our DIAG) | Slot B (Tandem) | Conclusion | Next |
|---|---|---|---|
| Phase 1A guard tripped | (any) | Bug localized at guard site | Write test, fix, re-run Round 1 (skip Tandem this time) |
| 9/9 stations slip | 9/9 stations slip | Bug did NOT reproduce — verify git SHA matches v8.0.0 state | Re-check reproducer |
| 0/9 stations slip | 9/9 stations slip | Bug is SPECIFICALLY in our MFEM code | Phase 2 DIAG analysis on Slot A |
| 0/9 stations slip | 0/9 stations slip | Mesh file has a bug (both codes fail) | Check mesh; rebuild with alternate generator |
| 1–8/9 stations slip | 9/9 stations slip | Partial propagation bug in our code | Phase 2B (scale / partition-topology focus) + Phase 2 DIAG targeted at silent stations |
| Tandem output missing / job failed | — | Tandem comparison unavailable; proceed with Slot A DIAG only | Phase 2 (no reference) |

**Wall-clock budget:** Slot A and Slot B submitted concurrently →
single ~1.5 hr wait.  Post-job diff < 10 min.  Total Round 1 cost
~1.7 hr for TWO data points.  Compared to the Rev 6 flow (Slot A
alone → Phase 2 → Phase 3 → Escalation step 5 Tandem), this saves
4–6 hours when the bug is specifically in our propagation code.

**Phase Log:**
```
[ ] Baseline (Phase 0) passes:  __ / __
[ ] TEST_R007 passes on current code: __
[ ] Allocation SUs remaining: __
[ ] dt vs h_min/cp ratio: __
[ ] Pre-flight outcome: __  (or "ESCALATE")
[ ] 1-rank run stations slipping: __ / 9
[ ] N-rank run stations slipping: __ / 9
[ ] Case selected: A / B / B-partial / C / D
```

---

### Phase 1A — Cheap Hardening ALL silent-weld paths (45 min)
[REVISED per R-501, R-503, R-507, plus unchanged R-007, R-103, R-102]

**Goal:** catch silent-failure classes at near-zero implementation
cost.  Applies to BOTH interior and shared fault-flux branches AND
to the `is_fault` classification step that feeds both.

**Action 0 — Makefile + startup banner [NEW per R-501]:**

```diff
 # miniapps/seas/Makefile, near the top alongside other makevars
+# R-501: `make ... CXXFLAGS+=-D<FLAG>` may not reach the compile step
+# because MFEM's config.mk explicitly sets MFEM_CPPFLAGS/MFEM_CXXFLAGS.
+# Use this SEAS-specific variable for diagnostic flags:
+#   make seas_tpv102_driver \
+#        SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX -DSEAS_DIAG_GHOST_EXCHANGE"
+SEAS_EXTRA_CPPFLAGS ?=

 # ... in every compile recipe for TPV102:
-	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@
+	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI \
+		$(SEAS_EXTRA_CPPFLAGS) -c $< -o $@
```

```diff
 // drivers/tpv102_driver.cpp::main, immediately after MPI init
 int rank, nprocs;
 MPI_Comm_rank(comm, &rank);
 MPI_Comm_size(comm, &nprocs);
+if (rank == 0) {
+   // R-501 + R-606: startup banner so Phase 2 analysis can verify that
+   // the DIAG flag actually reached the binary.  Write to BOTH stderr
+   // and a dedicated file `build_info.txt` in the current working
+   // directory — stderr may be interleaved across ranks or truncated
+   // at Frontera scale, but the file is rank-0-only and deterministic.
+   // Every diagnostic run's log MUST show [BUILD] ... = ON before
+   // analysis proceeds.  Phase 2 grep is against `build_info.txt`,
+   // not stderr, for reliability.
+#ifdef SEAS_DIAG_FAULT_FLUX
+   const char *flag1 = "SEAS_DIAG_FAULT_FLUX = ON";
+#else
+   const char *flag1 = "SEAS_DIAG_FAULT_FLUX = OFF";
+#endif
+#ifdef SEAS_DIAG_GHOST_EXCHANGE
+   const char *flag2 = "SEAS_DIAG_GHOST_EXCHANGE = ON";
+#else
+   const char *flag2 = "SEAS_DIAG_GHOST_EXCHANGE = OFF";
+#endif
+   std::fprintf(stderr, "[BUILD] %s\n", flag1);
+   std::fprintf(stderr, "[BUILD] %s\n", flag2);
+
+   std::ofstream binfo("build_info.txt");
+   if (binfo.is_open()) {
+      binfo << "[BUILD] " << flag1 << "\n";
+      binfo << "[BUILD] " << flag2 << "\n";
+      binfo.close();
+   } else {
+      std::fprintf(stderr, "[WARNING] could not open build_info.txt "
+                           "for banner (Phase 2 should fall back to stderr)\n");
+   }
+}
```

**Action 1 — INTERIOR-branch silent welded-flux fallback [R-007, R-103]:**

```diff
// wave_operator.inl around line 830-844 (inside if (is_fault) / else)
 if (dof_idx >= 0 &&
     dof_idx < static_cast<int>(fault_dof_data_->size()) &&
     qpd_ptr != nullptr)
 {
    // ... canonical-frame fault flux ...
 }
 else
 {
-   // Fault face but no DOFData mapping — fall back to welded interior flux.
-   flux_.Interior(nor, Q_self, Q_nbr, F_h);
+   MFEM_ABORT("interior fault face " << f
+              << " has no FaultBasis/DOFData mapping ("
+              << "fb_idx=" << fb_idx
+              << ", dof_idx=" << dof_idx
+              << ", have_basis=" << have_basis
+              << "). fault_interior_face_to_basis_idx_ or "
+              << "fault_face_dof_offset_ is incomplete in the ctor.");
+   // Control never reaches here.  [R-507]
    ...
 }
```

**Action 1a — INTERIOR-branch `is_fault=false` misclassification guard [R-503, R-604]:**

Before the `is_fault` dispatch, cross-check that the classification
agrees with the ctor's `fault_interior_faces_` list.  To keep this
guard O(1) at the ~`(N_faces × N_fault_faces)` hot loop,
materialize the list as a `std::unordered_set<int>` in the ctor
[NEW per R-604].

```diff
 // wave_operator.hpp — new private member
+#include <unordered_set>
 ...
+// R-604: O(1) membership set for the interior fault-face
+// classification cross-check in ComputeFaceFluxRHS.  Populated
+// alongside fault_interior_faces_ in the ctor.
+std::unordered_set<int> fault_interior_face_set_;
```

```diff
 // wave_operator.inl — ctor, after fault_interior_faces_ is populated
 // (just after the face_bdr_attr_ loop near line 110-116)
+// R-604: populate the O(1)-lookup set mirroring fault_interior_faces_.
+fault_interior_face_set_.clear();
+fault_interior_face_set_.reserve(fault_interior_faces_.Size());
+for (int i = 0; i < fault_interior_faces_.Size(); i++) {
+   fault_interior_face_set_.insert(fault_interior_faces_[i]);
+}
+MFEM_VERIFY(static_cast<int>(fault_interior_face_set_.size())
+            == fault_interior_faces_.Size(),
+            "R-604: fault_interior_faces_ contains duplicate entries; "
+            "set size " << fault_interior_face_set_.size()
+            << " != array size " << fault_interior_faces_.Size());
```

```diff
// wave_operator.inl at line 695, immediately before the is_fault branch
    bool is_fault = (bdr_attr == bc_.fault_attr) && (bc_.fault_attr > 0);
+
+   // R-503 + R-604: classification cross-check.  fault_interior_faces_
+   // (and its O(1) sibling fault_interior_face_set_) is the ctor's
+   // canonical list of physical interior fault faces.  If a face
+   // appears in that list but is_fault == false, a ctor-populated
+   // face_bdr_attr_[f] is wrong and the code would silently weld this
+   // fault face.  Prerequisite: TEST_R007_AllInteriorFaultFacesMapped
+   // must pass on baseline code, OR this guard will fire on correct
+   // code that is being worked on.  See Phase 0 prerequisite check.
+   //
+   // R-604: O(1) hash-set lookup (was O(N_fault) linear scan in Rev 5).
+   const bool is_in_fault_list =
+      (fault_interior_face_set_.count(f) > 0);
+   MFEM_VERIFY(is_in_fault_list == is_fault,
+               "R-503 interior classification mismatch: face " << f
+               << " fault_interior_face_set_=" << is_in_fault_list
+               << " but is_fault=" << is_fault
+               << " (bdr_attr=" << bdr_attr
+               << ", fault_attr=" << bc_.fault_attr << "). "
+               << "face_bdr_attr_ / fault_interior_faces_ "
+               << "population inconsistent in ctor.");
```

**Action 1b — SHARED-branch silent welded-flux fallback [R-103]:**

```diff
// wave_operator.inl around line 1186-1190 (inside if (sf_fault) / else)
                if (sf_fault) {
                   ...
                   if (dof_idx >= 0 &&
                       dof_idx < static_cast<int>(fault_dof_data_->size())) {
                      ... // canonical-frame path
                      continue;
                   }
                   else
                   {
-                     flux_.Interior(nor, Q_self, Q_nbr, F_h);
+                     MFEM_ABORT("shared fault face " << sf
+                                << " has no shared_fault_dof_offset_ entry "
+                                << "(dof_idx=" << dof_idx
+                                << ", sf_to_basis_idx=" << basis_idx
+                                << "). shared_fault_dof_offset_ ctor "
+                                << "population incomplete for this face.");
+                     // Control never reaches here.  [R-507]
                   }
                }
```

**Action 2 — ghost-exchange deep-copy verification for ALL components [R-005, R-102]:**

```diff
// wave_operator.inl around line 930-936
 for (int c = 0; c < NUM_STATE; c++)
 {
    ...
    std::memcpy(nbr_data[c].GetData(), src.GetData(),
                src.Size() * sizeof(real_t));
-   if (c == 0)
-   {
-      MFEM_VERIFY(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(),
-                  "nbr_data[0] aliases q_gf.FaceNbrData() — H2 ...");
-   }
+   // R-005 / R-102: MFEM_VERIFY (always-on in Release) for EVERY c,
+   // not only c==0.  MFEM_ASSERT would compile out in Release.  An
+   // aliasing bug on c=3 (SXY) or c=5 (SXZ) would silently corrupt
+   // strike-stress ghost data — plausibly the v8.0.0 bulk-coupling
+   // failure mode.
+   MFEM_VERIFY(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(),
+               "nbr_data[" << c << "] aliases q_gf.FaceNbrData() — "
+               "silent aliasing bug would corrupt cross-rank bulk "
+               "wave propagation for this component.  Do NOT "
+               "replace with MFEM_ASSERT.");
 }
```

**Safety note on `MFEM_ABORT` inside parallel `for(sf)` loop:**
`MFEM_ABORT` → `MPI_Abort(MPI_COMM_WORLD, 1)` is unconditional; no
collective blocking, no deadlock.  Standard MFEM pattern.

**Action 3 — rebuild + rerun Phase 1:**

```bash
make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS=""   # non-DIAG binary first
# Rerun Phase 1 Runs 1a and 1b.  If any MFEM_ABORT/VERIFY fires, bug found.
make -j test 2>&1 | tee /tmp/phase1A_post_guards_tests.log
diff /tmp/phase0_baseline_tests.log /tmp/phase1A_post_guards_tests.log
```

**Decision rule:**
- `MFEM_ABORT` or `MFEM_VERIFY` fires → bug localized.  Commit the
  guard + write a unit test (Phase 4).  Skip Phase 2.
- No abort, symptom persists → Phase 2 on branch identified in Phase 1.

**Phase Log:**
```
[ ] Action 0 (Makefile + banner) applied: __
[ ] Action 1 (interior silent-weld abort): __
[ ] Action 1a (interior is_fault misclassification verify) [NEW]: __
[ ] Action 1b (shared silent-weld abort): __
[ ] Action 2 (all-components MFEM_VERIFY): __
[ ] Startup banner confirms SEAS_DIAG_FAULT_FLUX = ON/OFF: __
[ ] Rerun Phase 1: aborts fired? __
[ ] Baseline test diff: __ new failures / __ new passes
[ ] Bug caught by Phase 1A: __
```

---

### Phase 2 — Chain Instrumentation on ONE branch (1.5 hr)
[REVISED per R-104, R-105, R-106, R-501, R-504]

**Goal:** identify the specific step where the signal is lost.  Only
runs the branch identified by Phase 1.  Instruments BOTH hypocenter
AND off-hypocenter fault QPs (R-504) so the propagation failure —
not just the local injection — is observable.

**Actions:**
1. Add `bool diag_print = false` to `DOFData` (conditional on flag):
   ```cpp
   // miniapps/seas/dynamic/fault_face_flux.hpp
   struct DOFData {
      // ... existing fields ...
   #ifdef SEAS_DIAG_FAULT_FLUX
      bool diag_print = false;
   #endif
   };
   ```

**2a. Global step-gate declarations [NEW per R-702]:**

```cpp
// miniapps/seas/dynamic/fault_face_flux.hpp — near the top, inside
// the seas namespace:
#ifdef SEAS_DIAG_FAULT_FLUX
// R-702: global step-gate flags visible to every translation unit
// that contains a [EVAL]/[FLUX]/[BULK] DIAG block.  Updated by the
// driver once per RK4 step; read-only everywhere else.
extern int g_seas_diag_step;       // current RK4 step index
extern int g_seas_diag_step_gate;  // 1 iff this step should print
#endif
```

```cpp
// miniapps/seas/dynamic/fault_face_flux.cpp — any single TU defines:
#ifdef SEAS_DIAG_FAULT_FLUX
int g_seas_diag_step      = 0;
int g_seas_diag_step_gate = 0;
#endif
```

2. DIAG block inside `FaultFaceFlux::Evaluate` [throttled per R-702]:
   ```cpp
   #ifdef SEAS_DIAG_FAULT_FLUX
   if (data.diag_print && g_seas_diag_step_gate) {
      std::fprintf(stderr,
         "[EVAL] step=%d V_abs=%.3e V1=%.3e V2=%.3e\n"
         "[EVAL] tau1_corr_L=%.3e tau2_corr_L=%.3e (perturbation)\n"
         "[EVAL] data.tau1_corr=%.3e data.tau2_corr=%.3e (tau_0 + pert)\n"
         "[EVAL] Q_plus_local:  ...\n"
         "[EVAL] Q_minus_local: ...\n"
         "[EVAL] Q_imp_plus:    ...\n"
         "[EVAL] Q_imp_minus:   ...\n",
         g_seas_diag_step,
         V_abs, V1, V2, tau1_corr, tau2_corr,
         data.tau1_corr, data.tau2_corr,
         /* ... */);
   }
   #endif
   ```

3. DIAG block in `wave_operator.inl` [throttled per R-702]:
   ```cpp
   #ifdef SEAS_DIAG_FAULT_FLUX
   const bool _diag_on =
      (dof_idx >= 0 &&
       dof_idx < static_cast<int>(fault_dof_data_->size()) &&
       (*fault_dof_data_)[dof_idx].diag_print &&
       g_seas_diag_step_gate);
   if (_diag_on) {
      std::fprintf(stderr,
         "[FLUX] step=%d dof_idx=%d nor: %.3e %.3e %.3e\n"
         "[FLUX] can_n: %.3e %.3e %.3e\n"
         "[FLUX] can_t1: %.3e %.3e %.3e\n"
         "[FLUX] can_t2: %.3e %.3e %.3e\n"
         "[FLUX] Q_imp_plus_g:  ...\n"
         "[FLUX] Q_imp_minus_g: ...\n"
         "[FLUX] F_h_total:     ...\n"
         "[FLUX] w=%.3e ndof=%d\n",
         g_seas_diag_step, dof_idx, nor[0], nor[1], nor[2], ...);
      int diag_bulk_dof = VX * ndof_total_ + dof_offset1 + 0;
      std::fprintf(stderr, "[FLUX] rhs[bulk VX] pre=%.3e ", rhs[diag_bulk_dof]);
      // ... accumulation ...
      std::fprintf(stderr, "post=%.3e\n", rhs[diag_bulk_dof]);
   }
   #endif
   ```
   The `dof_idx` is in every `[FLUX]` line so hypocenter vs off-hypo
   output are distinguishable when both are instrumented.

**3a. Driver updates `g_seas_diag_step_gate` per RK4 step [NEW per R-702]:**

```cpp
// tpv102_driver.cpp — top of RK4 main loop
for (int step = 0; step < nsteps; step++) {
#ifdef SEAS_DIAG_FAULT_FLUX
   // R-702: routine cadence = 1 in 100; breakaway cadence = 1 in 10
   // during the last ~25% of the run (approximately t >= 0.9 s when
   // tfinal = 1.2 s, covering the V_max > 0.1 m/s breakaway onset
   // through saturation).
   const int break_start = (3 * nsteps) / 4;
   const int routine = 100, break_cadence = 10;
   g_seas_diag_step = step;
   if (step >= break_start) {
      g_seas_diag_step_gate = (step % break_cadence == 0) ? 1 : 0;
   } else {
      g_seas_diag_step_gate = (step % routine == 0) ? 1 : 0;
   }
#endif
   // ... RK4 stages ...
}
```

4. **Driver: flag TWO DOFs for DIAG — hypocenter AND off-hypo witness [R-504, R-602]:**

   Both DOFs MUST be gated by a distance filter so that only the rank
   that actually contains the target station sets the flag.  Without
   the gate, `FindNearestDOF` returns a valid local index on every
   rank, producing N × 2 spurious DIAG streams (R-602).

   ```cpp
   #ifdef SEAS_DIAG_FAULT_FLUX
   // R-602: distance gate threshold — one 200m mesh cell wide so
   // the filter is tight (nominally one and only one DOF will be
   // within range on any rank) but generous enough to accommodate
   // mesh-to-station coordinate mismatch.
   constexpr real_t kDiagDofMatchRadius = 500.0;  // m

   int diag_hypo_dof = -1;
   int diag_offhypo_dof = -1;

   if (num_fault_local > 0) {
      // --- Hypocenter DOF (R-602: now gated by distance) ---
      int hypo_dof = FindNearestDOF(
         {TPV102Params::hypo_along_strike,
          TPV102Params::hypo_down_dip, "hypo"},
         fault_coords, num_fault_local);
      if (hypo_dof >= 0) {
         Vector hpos = fault_coords[hypo_dof];
         real_t hdx = hpos(0) - TPV102Params::hypo_along_strike;
         real_t hdz = std::abs(hpos(2)) - TPV102Params::hypo_down_dip;
         real_t hdist = std::sqrt(hdx*hdx + hdz*hdz);
         if (hdist < kDiagDofMatchRadius) {
            dof_data[hypo_dof].diag_print = true;
            diag_hypo_dof = hypo_dof;
            std::cerr << "[diag] rank " << rank
                      << " tagging hypo DOF " << hypo_dof
                      << " (dist=" << hdist << " m)\n";
         }
      }

      // --- Off-hypocenter station DOF [R-504, R-605, R-704] ---
      // Dual selection:
      //   flt_0_3   at (along=0, down-dip=3 km)  — 4.5 km UP-dip from hypo,
      //     arrival 1.30 s.  Default for tfinal ∈ [1.0, 1.35) s.
      //     May collocate with hypocenter rank at small np (R-704).
      //   flt_0_12  at (along=0, down-dip=12 km) — 4.5 km DOWN-dip from hypo,
      //     arrival 1.30 s.  Near-guaranteed different partition at np ≥ 100
      //     (down-dip crosses Z-partitions quickly).  Use when
      //     tfinal >= 1.35 s is available OR hypo+flt_0_3 collocation needs
      //     to be broken.
      //   flt_9_7.5 at (along=9 km, down-dip=7.5 km) — arrival 2.6 s —
      //     too far for tfinal ≤ 1.5 s; NOT selected in standard DIAG.
      //
      // Selection rule: use flt_0_3 by default; if R-704 cross-rank check
      // below fails (hypo + off-hypo on same rank), the user must re-run
      // with SEAS_DIAG_OFFHYPO_FAR=1 which replaces flt_0_3 with flt_0_12.

      struct OffHypoChoice {
         const char *name;
         real_t along_strike;
         real_t down_dip;   // positive; z = -down_dip in grid
      };
      #if defined(SEAS_DIAG_OFFHYPO_FAR)
      const OffHypoChoice choice = {"flt_0_12", 0.0, 12000.0};
      #else
      const OffHypoChoice choice = {"flt_0_3",  0.0,  3000.0};
      #endif

      int offhypo_dof = FindNearestDOF(
         {choice.along_strike, choice.down_dip, choice.name},
         fault_coords, num_fault_local);
      if (offhypo_dof >= 0) {
         Vector opos = fault_coords[offhypo_dof];
         real_t odx = opos(0)        - choice.along_strike;
         real_t odz = std::abs(opos(2)) - choice.down_dip;
         real_t odist = std::sqrt(odx*odx + odz*odz);
         if (odist < kDiagDofMatchRadius) {
            dof_data[offhypo_dof].diag_print = true;
            diag_offhypo_dof = offhypo_dof;
            std::cerr << "[diag] rank " << rank
                      << " tagging off-hypo DOF " << offhypo_dof
                      << " (" << choice.name
                      << ", dist=" << odist << " m)\n";
         }
      }
   }

   // Bound total count at 2 per rank defensively:
   int local_diag_count = 0;
   for (int i = 0; i < num_fault_local; i++) {
      if (dof_data[i].diag_print) { local_diag_count++; }
   }
   MFEM_VERIFY(local_diag_count <= 2,
               "More than 2 DOFs flagged diag_print=true on rank "
               << rank << " (count=" << local_diag_count << ").  "
               << "Only hypocenter and one off-hypo witness should be flagged.");

   // R-602: globally, exactly ZERO OR ONE rank should have
   // diag_hypo_dof >= 0, and similarly for diag_offhypo_dof.  Verify
   // via Allreduce SUM.  This prevents N-rank × 2 log-noise regression.
   int hypo_flag_global = (diag_hypo_dof >= 0) ? 1 : 0;
   int offhypo_flag_global = (diag_offhypo_dof >= 0) ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &hypo_flag_global, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, &offhypo_flag_global, 1, MPI_INT, MPI_SUM, comm);
   MFEM_VERIFY(hypo_flag_global <= 1,
               "R-602: " << hypo_flag_global << " ranks flagged hypo DOF; "
               "expected 0 or 1.  Distance filter may be too wide OR "
               "the mesh has two hypocenter DOFs on different ranks.");
   MFEM_VERIFY(offhypo_flag_global <= 1,
               "R-602: " << offhypo_flag_global << " ranks flagged off-hypo "
               "DOF; expected 0 or 1.");

   // R-704: verify hypo and off-hypo are NOT on the same rank.  If
   // they are, the cross-rank propagation test collapses to an
   // intra-rank test (no ParMesh::ExchangeFaceNbrData involvement).
   int on_same_rank_local = (diag_hypo_dof >= 0 && diag_offhypo_dof >= 0) ? 1 : 0;
   int on_same_rank_global = 0;
   MPI_Allreduce(&on_same_rank_local, &on_same_rank_global, 1,
                 MPI_INT, MPI_SUM, comm);
   MFEM_VERIFY(on_same_rank_global == 0,
               "R-704: hypocenter DOF and off-hypo DOF both landed on rank "
               << rank << " (global count=" << on_same_rank_global
               << ").  Cross-rank propagation test defeated.  Pick a "
               "farther off-hypo station (e.g. flt_0_12 for tfinal=1.5s) "
               "or accept the degraded test at small rank counts.");

   if (rank == 0) {
      std::cerr << "[diag] global: hypo_flagged_ranks=" << hypo_flag_global
                << " offhypo_flagged_ranks=" << offhypo_flag_global
                << " on_same_rank=" << on_same_rank_global << "\n";
   }
   #endif
   ```

**4a. Driver: flag 2–3 BULK-DOF propagation monitors [NEW per R-701]:**

The hypo and off-hypo fault-DOF DIAG covers steps (1) and (4) of the
causal chain (hypo injects; off-hypo reads).  Steps (2) and (3) —
the bulk wave travelling from hypo's neighbourhood to off-hypo's
neighbourhood — are NOT instrumented by the fault-only DIAG.  A
volume-term, mass-inverse, or interior-flux bug that kills waves
mid-flight produces identical symptoms to "hypo injection broken":
correct hypo values, silent off-hypo.

Bulk-DOF monitors at fixed distances along the up-dip direction
break this ambiguity:

```cpp
#ifdef SEAS_DIAG_FAULT_FLUX
// R-701: flag BULK elements at 1 / 2 / 3 km up-dip from hypocenter
// as propagation-path witnesses.  Slightly off the fault (y ≈ 300 m
// into bulk) so they reside in VOLUME elements, not on the fault
// interface.  S-wave (cs=3464) arrival times: 0.29 / 0.58 / 0.87 s
// — all well within tfinal=1.2 s.
std::vector<int> diag_bulk_elem_ids;
const std::vector<real_t> diag_bulk_distances_km = {1.0, 2.0, 3.0};
const real_t bulk_y_offset = 300.0;        // m into bulk, off the fault
const real_t bulk_match_radius = 500.0;    // m

for (real_t d_km : diag_bulk_distances_km) {
   const real_t target_x = TPV102Params::hypo_along_strike;   // 0
   const real_t target_y = bulk_y_offset;                     // +300 m
   const real_t target_z = -(TPV102Params::hypo_down_dip - d_km * 1000.0);

   int best_e = -1;
   real_t best_d2 = std::numeric_limits<real_t>::max();
   Array<int> verts;
   for (int e = 0; e < wave.GetNE(); e++) {
      mesh->GetElementVertices(e, verts);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < verts.Size(); v++) {
         const real_t *p = mesh->GetVertex(verts[v]);
         cx += p[0]; cy += p[1]; cz += p[2];
      }
      cx /= verts.Size(); cy /= verts.Size(); cz /= verts.Size();
      real_t dd = (cx-target_x)*(cx-target_x) +
                  (cy-target_y)*(cy-target_y) +
                  (cz-target_z)*(cz-target_z);
      if (dd < best_d2) { best_d2 = dd; best_e = e; }
   }
   if (best_e >= 0 && std::sqrt(best_d2) < bulk_match_radius) {
      diag_bulk_elem_ids.push_back(best_e);
      std::cerr << "[diag] rank " << rank
                << " tagging bulk-prop-witness elem " << best_e
                << " at ~" << d_km << " km up-dip "
                << "(dist=" << std::sqrt(best_d2) << " m)\n";
   }
}

// R-701 + R-602 sibling: verify each bulk monitor is on exactly one
// rank globally (SUM over ranks of "I own this distance monitor").
for (size_t i = 0; i < diag_bulk_distances_km.size(); i++) {
   int mine = (i < diag_bulk_elem_ids.size()) ? 1 : 0;
   int global = 0;
   MPI_Allreduce(&mine, &global, 1, MPI_INT, MPI_SUM, comm);
   MFEM_VERIFY(global <= 1,
               "R-701: bulk-prop monitor at " << diag_bulk_distances_km[i]
               << " km flagged on " << global << " ranks (expected 0 or 1). "
               "Shrink bulk_match_radius or pick a y_offset that "
               "avoids double-occupancy across partition boundaries.");
   if (rank == 0 && global == 0) {
      std::cerr << "[diag] warning: no rank found a bulk elem within "
                << bulk_match_radius << " m of the "
                << diag_bulk_distances_km[i] << " km monitor target. "
                   "Monitor skipped.\n";
   }
}
#endif
```

**4b. Driver: emit `[BULK]` samples from the monitors [NEW per R-701]:**

At the top or bottom of each RK4 step (after the step has
completed), sample bulk `Q` at each monitored element — gated by
the R-702 step throttle:

```cpp
#ifdef SEAS_DIAG_FAULT_FLUX
// R-701 + R-702: emit at gated step cadence (every 100 steps
// normally; every 10 in breakaway window).
if (g_seas_diag_step_gate) {
   for (size_t i = 0; i < diag_bulk_elem_ids.size(); i++) {
      int e = diag_bulk_elem_ids[i];
      int dof_offset = e * ndof_per_el;
      real_t q_vx = 0.0, q_sxy = 0.0, q_sxz = 0.0;
      for (int j = 0; j < ndof_per_el; j++) {
         q_vx  = std::max(q_vx,
                          std::abs(Q[VX  * ndof_total + dof_offset + j]));
         q_sxy = std::max(q_sxy,
                          std::abs(Q[SXY * ndof_total + dof_offset + j]));
         q_sxz = std::max(q_sxz,
                          std::abs(Q[SXZ * ndof_total + dof_offset + j]));
      }
      std::fprintf(stderr,
                   "[BULK] rank=%d step=%d t=%.3e d=%gkm elem=%d "
                   "|VX|=%.3e |SXY|=%.3e |SXZ|=%.3e\n",
                   rank, step, t, diag_bulk_distances_km[i], e,
                   q_vx, q_sxy, q_sxz);
   }
}
#endif
```

5. Re-assertion invariant [R-106]:
   ```cpp
   for (int step = 0; step < nsteps; step++) {
   #ifdef SEAS_DIAG_FAULT_FLUX
      if (diag_hypo_dof >= 0) {
         MFEM_VERIFY(dof_data[diag_hypo_dof].diag_print,
                     "diag_print cleared between RK4 steps — "
                     "a code path is assigning DOFData{} wholesale.");
      }
      if (diag_offhypo_dof >= 0) {
         MFEM_VERIFY(dof_data[diag_offhypo_dof].diag_print,
                     "off-hypo diag_print cleared.");
      }
   #endif
      ...
   }
   ```

6. Rebuild **with the SEAS_EXTRA_CPPFLAGS variable [R-501]**:
   ```bash
   make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX"
   ```

7. Run the Phase 1 reproducer with `tfinal = 1.2 s` (NOT 0.3 s,
   per R-502).  Capture stderr.  FIRST: `grep "\[BUILD\]" build_info.txt`
   (preferred) or stderr fallback to confirm the flag is ON (R-606).

**Expected values at HYPOCENTER after breakaway (t ≥ 1.0 s, V_abs ≈ 7.7 m/s):**
| Quantity | Expected | Frame |
|---|---|---|
| V_abs | 7.7 m/s | scalar |
| V1, V2 | 0, 7.7 | scalar (BP5: t1=dip, t2=strike) |
| tau2_corr_L (Evaluate local var) | -36e6 Pa | perturbation |
| data.tau2_corr | +40e6 Pa (= tau2_0 + perturbation) | total |
| Q_imp_plus[VZ] (fault-local) | -3.85 m/s | can_t2=strike |
| Q_imp_plus_g[VX] (global) | -3.85 m/s | post T_can |
| Q_imp_plus_g[SXY] | +36e6 Pa | from R-003 derivation |
| F_h_total[SXY] (stress) | 0 ± 1 Pa/s | exact by stress continuity |
| F_h_total[VX] (strike velocity) | ≈ +2.2e2 m/s² | analytical |
| rhs[VX] delta per Mult call at hypo | ≈ 3.7e5 m⁴/s² per face QP | pre-mass-inverse |
| k[VX] at fault-adjacent bulk DOF | ≈ 0.3 m/s² | post-mass-inverse |

**Expected values at OFF-HYPO `flt_0_3` QP after breakaway [R-504, R-605]:**

**Wave-arrival timing** (R-605):
- Hypocenter at `(0, 0, -7500)`; station `flt_0_3` at `(0, 0, -3000)`.
  Distance = 4500 m along-fault (up-dip).
- S-wave speed `cs = √(μ/ρ) ≈ 3464 m/s` for TPV102 material.
- `t_arrive ≈ 4500 / 3464 ≈ 1.30 s`.
- **At `tfinal=1.2 s` the S-wave front is still ~150 m short of `flt_0_3`.**
  The station is therefore EXPECTED to be pristine at `tfinal=1.2 s`,
  and the "O(1e5–1e6 Pa)" entries below apply only at `tfinal ≥ 1.5 s`.

| Quantity | Expected at t=1.2s | Expected at t=1.5s | Frame |
|---|---|---|---|
| Q_plus_local[SXZ] | ≈ 0 (front not arrived) | O(1e5–1e6 Pa) — arriving shear wave | fault-local |
| tau2_trial (Evaluate local) | ≈ 0 | O(1e5–1e6 Pa) | matches Q_plus |
| tau2_total | ≈ 75e6 | ≈ 75e6 + trial (slight perturbation) | total |
| V_abs from Brent | = V_ini = 1e-12 | > V_ini by at least a few OOM | scalar |
| Q_imp_plus_g near off-hypo | ≈ 0 | bulk velocity O(mm/s to cm/s) | global |
| F_h_total[VX] at off-hypo face | ≈ 0 | O(1e1 – 1e2 m/s²) | global |

**[R-605 + R-701 bulk-DOF propagation monitors — REQUIRED for every DIAG run]:**

Bulk-DOF witnesses at 1, 2, 3 km up-dip from hypocenter (at
`y ≈ 300 m` off the fault, into bulk) break the "correct hypo /
silent off-hypo" ambiguity.  S-wave (cs=3464) arrival times:
`0.29 / 0.58 / 0.87 s` — all well within `tfinal=1.2 s`.

| Distance (up-dip from hypo) | Arrival time (cs=3464) | Expected `max|Q[VX]|` at t=1.2s | Expected `max|Q[SXY]|` at t=1.2s |
|---|---|---|---|
| 1 km | 0.29 s | O(0.1–1 m/s) | O(0.5–5 MPa) |
| 2 km | 0.58 s | O(0.05–0.5 m/s) | O(0.1–1 MPa) |
| 3 km | 0.87 s | O(0.01–0.1 m/s) | O(0.01–0.1 MPa) |

If at `tfinal=1.2 s` any of these bulk monitors is stuck at O(1e-8),
waves aren't propagating that far from the hypocenter — the bug is in
BULK PROPAGATION (volume-term / mass-inverse / inter-element flux),
NOT in fault-flux injection.  The off-hypo fault-QP DIAG (`flt_0_3`)
at 4.5 km is NOT expected to slip at `tfinal=1.2 s` (arrival 1.30 s)
on correct code — it is dispositive only at `tfinal ≥ 1.5 s`.

**Decision rule [R-701 propagation monitors + R-605 timing]:**

| 1 km | 2 km | 3 km | Hypo QP | off-hypo (`flt_0_3`) | Diagnosis | Next |
|---|---|---|---|---|---|---|
| — | — | — | wrong | — | Formula / rotation bug | Phase 3B or 3C |
| ≈ 0 | ≈ 0 | ≈ 0 | correct | — | Fault injection succeeds locally but adjacent bulk is silent → face-flux accumulation or per-component mass-inverse bug | Phase 3A Fixture A first |
| grow | ≈ 0 | ≈ 0 | correct | — | Wave dies within first 1–2 km → strong in-volume dissipation (ApplyPMLDamping on a real region? stiffness sign?) | Phase 3A Fixture B |
| grow | grow | ≈ 0 | correct | — | Wave dies 2–3 km → volume-term / mass-inverse per-component bug | Phase 3A Fixture B |
| grow | grow | grow | correct | ≈ 0 at t=1.2s (NORMAL) | Bulk propagates correctly; off-hypo hasn't arrived yet | Re-run `tfinal=1.5 s`; expect ≈ 0 bulk-monitor silence means V_ini at off-hypo |
| grow | grow | grow | correct | ≈ 0 at t=1.5s | Wave energy dies between 3 km and 4.5 km | Phase 3A + Phase 2B scale escalation |
| grow | grow | grow | correct | Q_plus grows, V_abs = V_ini at t=1.5s | Trial-traction / friction-solver reading stale bulk | Investigate `shape(xi_q)` interpolation in Evaluate caller |
| grow | grow | grow | correct | Q_plus grows, V_abs grows | Bulk + fault both working; bug is elsewhere (re-verify reproducer) | — |

**Phase Log:**
```
[ ] Phase 2 code + banner applied
[ ] Make build uses SEAS_EXTRA_CPPFLAGS=…
[ ] Banner check: SEAS_DIAG_FAULT_FLUX = __ in .out file
[ ] diag_hypo_dof set: __ (rank __)
[ ] diag_offhypo_dof set: __ (rank __)
[ ] Hypo observed vs expected — first divergence: __
[ ] Off-hypo observed vs expected — first divergence: __
[ ] Branch selected (3A/3B/3C): __
```

---

### Phase 2B — MPI / Scale Escalation (1-2 hr)
[REVISED per R-303, R-502, R-505]

**Goal:** triggered by Phase 1 Case C (local ≤ 8 propagates).
Identify the smallest rank count at which the bug manifests.

**Tier structure:**

| Tier | Scale | tfinal | Wall | Flags |
|---|---|---|---|---|
| 2B-local | 4 ranks, 1000m | 1.2 s | ~10–20 min | SEAS_DIAG_GHOST_EXCHANGE |
| 2B-100rank | 100 ranks, 2 nodes, 200m | 1.2 s | ~1.2 hr | DIAG_FAULT_FLUX + GHOST_EXCHANGE |
| 2B-400rank | 400 ranks, 200m | 1.2 s | ~1.5 hr | DIAG_FAULT_FLUX + GHOST_EXCHANGE |

50-rank intentionally skipped for 2B — at 50 ranks + tfinal=1.2s the
wall-clock is ~2.4 hr (exceeds 2-hr dev limit).  100 ranks / 2 nodes
gives ~1.2 hr wall-clock and fits.

**Actions for 2B-100rank (most common case):**
```bash
cp jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch \
   jobs/tpv102/tpv102_200m_p1_1.2s_100rank_DIAG.sbatch
# Edit:
#   #SBATCH -n 100 -N 2
#   --tfinal 1.2
#   make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX -DSEAS_DIAG_GHOST_EXCHANGE"
#   --debug-qnorm
sbatch jobs/tpv102/tpv102_200m_p1_1.2s_100rank_DIAG.sbatch
```

**Per-component ghost max expected values (same as prior revs):**
- c ∈ {SXX(0), SXY(3), SXZ(5), VX(6)}: grow from 0 to O(MPa) or
  O(m/s) on hypo-adjacent ranks by breakaway onset (t ≈ 1.0 s).
- c ∈ {SYY(1), SZZ(2), SYZ(4), VY(7), VZ(8)}: growth smaller or
  zero (TPV102 pure strike-slip).

**Decision rule:**

| Outcome | Conclusion |
|---|---|
| 2B-local shows ghost max stuck for some c | ExchangeFaceNbrData bug on component c → Phase 3A ghost-exchange path |
| 2B-local clean, 2B-100rank anomalous | 100-rank topology dependent → Phase 3A with 100-rank DIAG rerun |
| 2B-local + 2B-100rank clean, 2B-400rank anomalous | 400-rank specific → Phase 3A with 400-rank DIAG rerun |
| All tiers clean | Bug requires 400 ranks AND full 1.5s breakaway; Tandem comparison (Escalation 5) |

---

### Phase 3A — Accumulation / Mass Inverse / RK4
[REVISED per R-603, R-607]

**Goal:** verify that known-correct F_h produces bulk perturbation
AT AND AWAY FROM the fault-adjacent element.  **Highest prior** given
SeisSol audit + R-101 3-ULP pass rule out Evaluate and T_can.

Rev 6 adds a second fixture (B) to cover volume-integration bugs
that the 2-tet Fixture A cannot reach.

**Actions:**

**1. Write `tests/unit/test_fault_flux_accumulation.cpp` with TWO fixtures.**

**Fixture A (face-flux accumulation only) — existing 2-tet mesh:**
```cpp
// Fixture A: 2 tets sharing one internal fault face.
// Measures: does a known F_h produce bulk dQ/dt at the face-adjacent
// elements?
Mesh mesh = BuildTwoTetSharedFaultMeshInline();   // existing in test_r101
real_t h_A = ComputeMinInscribedDiameter(mesh);   // [R-607] param thresh.
WaveOperator<...> wave(mesh, ...);
SetDOFData_HypocenterStyle(wave, dof_data,
                           tau2_0 = tau_ini + peak_dtau,
                           psi    = psi_equilibrium);
Vector Q_zero(NUM_STATE * wave.GetNDofTotal()); Q_zero = 0.0;
Vector k(NUM_STATE * wave.GetNDofTotal());
wave.Mult(Q_zero, k);

real_t max_k = k.Normlinf();
// R-607: thresholds scale as 1/h_A because F_h is an O(1) surface
// force divided by an O(h) mass-matrix entry, so k ~ 1/h.  Fixture A
// was calibrated at h=1 m giving max_k ≈ 0.3–0.5; generalize as:
real_t k_low  = 0.05  * (1.0 / h_A);
real_t k_high = 5.0   * (1.0 / h_A);
TEST_ASSERT(max_k > k_low,
            "Fixture A: bulk dQ/dt under-scaled at h=" + std::to_string(h_A));
TEST_ASSERT(max_k < k_high,
            "Fixture A: bulk dQ/dt over-scaled / runaway at h=" + std::to_string(h_A));
```

**Fixture B (propagation corridor) — NEW per R-603:**
```cpp
// Fixture B: N=5 tets stacked along +x with fault at the leftmost
// (element 0) left-face.  DOFData hypocenter-style on that face.
// Measures: does bulk energy RADIATE AWAY from the fault over
// multiple time steps and multiple elements?
//
// A volume-integration bug that injects correctly into element 0 but
// fails to propagate to elements 1-4 CANNOT be caught by Fixture A.
// Fixture B is the canonical propagation test.
constexpr int N_elems = 5;
constexpr real_t h_B = 1.0;   // 1 m cells for easy ULP reasoning
Mesh mesh = BuildCorridorMesh(N_elems, h_B);   // new helper
WaveOperator<...> wave(mesh, ...);
SetDOFData_HypocenterStyle(wave, dof_data, /* tau2_0 = */ tau_ini + peak_dtau,
                           /* psi = */ psi_equilibrium);

Vector Q(NUM_STATE * wave.GetNDofTotal()); Q = 0.0;
Vector k(NUM_STATE * wave.GetNDofTotal());

// Integrate with explicit Euler for simplicity (200 steps is enough
// for the S-wave to cross the full 5-element corridor at dt << h/cs).
real_t cs = wave.GetCs();
real_t dt = 0.1 * h_B / cs;                // safely sub-CFL
int n_steps_cross = static_cast<int>(std::ceil(N_elems * h_B / cs / dt));
int n_steps = 2 * n_steps_cross;           // give wave time to settle

for (int step = 0; step < n_steps; step++) {
   wave.Mult(Q, k);
   Q.Add(dt, k);
}

// Probe bulk |V| at the farthest element (i = N_elems - 1).
real_t q_far_vx = MaxAbsAtElemComponent(Q, wave, /*elem*/ N_elems - 1, VX);
real_t q_mid_vx = MaxAbsAtElemComponent(Q, wave, /*elem*/ N_elems / 2, VX);

// R-607 + Rev 6 note: threshold parameterization — at h=1m with peak_dtau
// on the fault face, the radiated bulk VX at any element >0 should be
// at least a few percent of V_peak ≈ 7 m/s by the time the wave
// traverses the corridor.
TEST_ASSERT(q_far_vx > 1e-3,
            "Fixture B: bulk |VX| at farthest element did not grow;"
            " volume integration may be zeroing out bulk propagation."
            " q_far_vx=" + std::to_string(q_far_vx));
TEST_ASSERT(q_mid_vx > 1e-3,
            "Fixture B: bulk |VX| at midpoint did not grow;"
            " volume integration bug likely."
            " q_mid_vx=" + std::to_string(q_mid_vx));
// Monotone-ish check: farthest should be < midpoint (wave is spreading,
// not concentrating at the far end); but both > 0.
TEST_ASSERT(q_mid_vx >= q_far_vx * 0.5,
            "Fixture B: farthest element has more energy than midpoint — "
            "unphysical amplification / reflection");
```

**BuildCorridorMesh helper:** place into `tests/common/corridor_mesh.hpp`;
vertices at `x = 0, h, 2h, ..., Nh`, `y,z ∈ {0, h}` forming `N` cubes
each split into 6 tets (standard), with a single face at `x=0` flagged
as fault-attribute.  Total `6N` tets, 1 fault face.

**2. If Fixture A fails:** instrument `rhs` at specific DOFs
   before/after the fault-face loop to isolate where it zeroes out.
   Candidate bugs:
   - `elem1_on_plus` inverted.
   - `accum_sign` cancellation across branches.
   - Face weight anomaly.
   - Mass matrix scaling (component-specific).
   - RK4 weights.

**3. If Fixture A passes but Fixture B fails (R-603 diagnostic):**
   The bug is in volume integration (`∇·σ` term inside
   `ComputeVolumeRHS`), not in face-flux accumulation.  Instrument
   `ComputeVolumeRHS` component-by-component:
   - Apply a known `Q` with `SXY = 1 MPa`, all others 0.  Expected
     `∇·σ` contributes to `rhs[VX]` per `∂σ_xy/∂y`.
   - Candidate bugs: wrong sign of stiffness, wrong component mapping,
     missing mass-inverse on some components, wrong basis-function order.

---

### Phase 3B — `flux_.Interior` at fault discontinuity
[INLINED per R-508, corrected per R-002]

```cpp
// tests/unit/test_godunov_flux_fault_discontinuity.cpp
TEST(GodunovFlux, FaultDiscontinuityAnalytical) {
   const real_t rho = 2670.0, cs = 3464.0, cp = 6000.0;
   const real_t mu = rho * cs * cs;
   const real_t lambda = rho * cp * cp - 2.0 * mu;
   GodunovFlux flux(lambda, mu, rho);

   const real_t V = 7.7;
   const real_t tau_imp = -36e6;
   real_t nor[3] = {0, -1, 0};

   real_t Q_plus_local[NUM_STATE] = {0};
   real_t Q_minus_local[NUM_STATE] = {0};
   Q_plus_local[SXZ]  = tau_imp;
   Q_minus_local[SXZ] = tau_imp;
   Q_plus_local[VZ]  = -V / 2.0;
   Q_minus_local[VZ] = +V / 2.0;

   real_t can_n[3] = {0, -1, 0};
   real_t can_t1[3] = {0, 0, -1};
   real_t can_t2[3] = {+1, 0, 0};
   DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
   GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T);
   real_t Q_plus_g[NUM_STATE], Q_minus_g[NUM_STATE];
   T.Mult(Q_plus_local,  Q_plus_g);
   T.Mult(Q_minus_local, Q_minus_g);

   real_t F_h[NUM_STATE];
   flux.Interior(nor, Q_plus_g, Q_minus_g, F_h);

   // Closed-form checks:
   // (1) Stress flux zero when Q_plus[SXZ] == Q_minus[SXZ].
   for (int c = 0; c < 6; c++) {
      TEST_NEAR(F_h[c], 0.0, 1.0);
   }
   // (2) Velocity flux magnitude O(10^2-10^4) m/s^2.
   real_t F_v_mag = std::sqrt(F_h[VX]*F_h[VX] +
                              F_h[VY]*F_h[VY] +
                              F_h[VZ]*F_h[VZ]);
   TEST_ASSERT(F_v_mag > 10.0 && F_v_mag < 1e5,
               "velocity flux magnitude wrong: "
               + std::to_string(F_v_mag));
}
```

---

### Phase 3C — T_can Rotation (prophylactic)
[INLINED per R-508, fully specified per R-003]

```cpp
// tests/unit/test_godunov_rotation.cpp
TEST(GodunovRotation, HypocenterStrikeFault_TPV102) {
   real_t can_n[3]  = {0, -1, 0};
   real_t can_t1[3] = {0, 0, -1};
   real_t can_t2[3] = {1, 0, 0};
   DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
   GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T);
   GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv);

   // T · Tinv = I
   DenseMatrix I(NUM_STATE); Mult(T, Tinv, I);
   for (int i = 0; i < NUM_STATE; i++)
      for (int j = 0; j < NUM_STATE; j++)
         TEST_NEAR(I(i,j), (i==j) ? 1.0 : 0.0, 1e-12);

   // Velocity rotation: v_t2 = strike = -3.85 maps to v_x global
   real_t Q_local[NUM_STATE] = {0};
   Q_local[VZ] = -3.85;
   real_t Q_global[NUM_STATE];
   T.Mult(Q_local, Q_global);
   TEST_NEAR(Q_global[VX], -3.85, 1e-10);
   TEST_NEAR(Q_global[VY],  0.0, 1e-10);
   TEST_NEAR(Q_global[VZ],  0.0, 1e-10);

   // Stress rotation: SXZ_local = -36e6 maps to SXY_global = +36e6
   real_t S_local[NUM_STATE] = {0};
   S_local[SXZ] = -36e6;
   real_t S_global[NUM_STATE];
   T.Mult(S_local, S_global);

   const real_t expected_sxy_global = +36e6;
   TEST_NEAR(S_global[SXY], expected_sxy_global, 1.0);
   TEST_NEAR(S_global[SXX], 0.0, 1.0);
   TEST_NEAR(S_global[SYY], 0.0, 1.0);
   TEST_NEAR(S_global[SZZ], 0.0, 1.0);
   TEST_NEAR(S_global[SXZ], 0.0, 1.0);
   TEST_NEAR(S_global[SYZ], 0.0, 1.0);
}
```

---

### Phase 4 — Harden with Unit Tests
[REVISED per R-603, R-604, R-606; existing coverage R-501, R-503, R-504]

**Goal:** every discovered bug gets a test; plus integration coverage.

**Tests (in dependency order — read top-to-bottom):**

- `TEST(SharedFaultPropagation, BulkRadiatesFromHypocenter_2Ranks)` —
  inline 2-tet fixture, hypocenter init, 200 RK4 steps, assert
  max Q_VX > 1e-3 m/s somewhere + monotonic growth.

- `TEST_R007_AllInteriorFaultFacesMapped` — ctor's
  `fault_interior_face_to_basis_idx_` covers every fault-attr
  interior face.  **Must pass on baseline code before R-503 guard
  is added.**

- `TEST_R103_AllSharedFaultFacesMapped` — ctor's
  `shared_fault_dof_offset_` covers every shared fault face.

- `TEST_R501_DiagFlagPropagates` — build with `SEAS_EXTRA_CPPFLAGS=-DSEAS_DIAG_FAULT_FLUX`;
  run driver briefly; assert `build_info.txt` (AND stderr as fallback)
  contains `[BUILD] SEAS_DIAG_FAULT_FLUX = ON` (R-501 + R-606 file path).

- `TEST_R503_InteriorFaultMisclassification_Aborts` — corrupt
  `face_bdr_attr_[fault_face] = 0` post-ctor; assert MFEM_VERIFY trips.

- `TEST_R504_OffHypoDiagFlagged_Single_Rank` — on a mesh where both
  hypo and flt_0_3 are on rank 0, assert TWO distinct DOF IDs
  appear in `[EVAL]` stderr.

- `TEST_R006_DOFData_size_stable` — `static_assert sizeof(DOFData) <= 256`
  unless `SEAS_DIAG_FAULT_FLUX` is defined.

- `TEST_R005_R102_GhostDeepCopy_AllComponents_Release` — Release
  build; fake aliasing on component 3; assert MFEM_VERIFY fires.

**NEW for Rev 6:**

- `TEST_R601_DtRatioReportedAtStartup` — run the driver with `--cfl 0.5
  --tfinal 0.001`; grep startup stdout for `CFL:`, `dt_cfl`, and
  `dt:` / `dt (override...)`; assert:
  ```
  ratio = dt / dt_cfl
  TEST_ASSERT(ratio <= 2.0,
              "R-601: dt exceeds the explicit-RK4 stability envelope (2x dt_cfl)");
  ```
  Sanity number on 1000m regular-tet mesh: `dt_cfl ≈ 3.78e-3 s`,
  `dt ≈ 3.78e-3 s`, ratio ≈ 1.0.

- `TEST_R602_OnlyOneRankFlagsHypoDof_MultiRank` — build 4-rank
  ParMesh on `tpv102_1000m.msh`.  Build driver.  After Phase 2
  Step 4 tagging logic runs, Allreduce(`diag_hypo_dof >= 0` ? 1 : 0)
  and assert `== 1` (exactly one rank tagged hypo).  Similar for
  `diag_offhypo_dof`.

- `TEST_R603_CorridorPropagation` — the Phase 3A Fixture B above,
  assert `q_far_vx > 1e-3` after 2·N·h/(cs·dt) steps on the 5-tet
  corridor mesh.

- `TEST_R604_FaultFaceSetLookup_ConstantTime` — build WaveOperator
  with N_fault_faces ∈ {10, 100, 1000}.  For each, run a tight loop
  of `fault_interior_face_set_.count(f)` over `N_test = 10000`
  random face indices and measure wall-clock.  Assert ratio of times
  for (N=1000) / (N=10) is < 3× (near-constant), not ~100× (linear).

- `TEST_R605_OffHypoArrivalTiming_Sanity` — given the corridor
  Fixture B (R-603), verify analytically: bulk `|VX|` at element
  `k` first exceeds `1e-4 m/s` at step `n*` where
  `n* · dt ≈ k · h / cs`.  If numerical front speed differs from
  analytical `cs` by >10%, something in the volume or face
  integration is wrong.

- `TEST_R606_BuildBannerWrittenToFile` — run driver with
  `-DSEAS_DIAG_FAULT_FLUX`; assert file `build_info.txt` exists in
  CWD, contains `[BUILD] SEAS_DIAG_FAULT_FLUX = ON`, and is written
  only by rank 0 (no duplicate appends).

- `TEST_R607_Fixture_Thresholds_Scale_With_h` — run Fixture A at
  h = 1 m AND h = 0.5 m AND h = 2 m.  Assert that the k_low /
  k_high thresholds (derived from h) bracket `max_k` in all three
  cases.  Regression guard for future fixture-h changes.

**NEW for Rev 7:**

- `TEST_R701_BulkPropagationWitness_5TetCorridor` — uses the
  Fixture B corridor mesh.  Run 200 RK4 steps.  At each step,
  record `max|Q[VX]|` at element 1 (adjacent to fault) AND element 4
  (far end).  Assert both grow monotonically over the last 50
  steps; assert `Q[elem=4] > 0.1 × Q[elem=1]` (propagation, not
  decay).  If the ratio is below 0.1, volume integration damps
  bulk waves too much — exactly the symptom R-701 targets.

- `TEST_R702_DiagLog_Bounded_400Steps_Single_Rank` — compile with
  `SEAS_EXTRA_CPPFLAGS=-DSEAS_DIAG_FAULT_FLUX`; run 400 RK4 steps
  with one DIAG-flagged fault DOF; capture stderr line count.
  Assert `<= 500 lines` (with the routine-100 + breakaway-10
  throttle: ~4 routine windows × 4 stages × 10 lines/stage +
  ~100-line breakaway emission).  Without the throttle, the line
  count would be > 16000.

- `TEST_R703_TandemReference_StubbedCompare` — stub test: feed a
  synthetic 9-station V2(t) CSV and an expected-match CSV through
  `compare_stations.py`; assert the comparator returns non-zero
  exit code when stations diverge > 10% and zero when identical.
  Validates the diff script itself before Tandem runs use it.

- `TEST_R704_HypoAndOffhypo_Different_Ranks_AtProductionScale` —
  build a 4-rank ParMesh on `tpv102_1000m.msh` (small but non-
  trivial partition), run Phase 2 Step 4 tagging, then
  `MPI_Allreduce(hypo_local && offhypo_local)`.  Assert result
  `== 0` (no rank has both flagged).  Fails deterministically on
  collocation regressions.

- `TEST_R705_CflLog_IsInformational_NotGating` — run the driver
  with `--dt X` where `X = 1.5 · dt_cfl` (within the 2× envelope);
  capture startup stdout; assert the driver EMITS a WARNING line
  but DOES NOT exit.  Then run with `X = 3 · dt_cfl` (outside
  envelope); assert the driver exits OR is clearly marked unstable
  within the first 20 steps.

---

### Phase 5 — Local 2-8 Rank Verification
[INLINED per R-508]

Same as Rev 3: rerun Phase 1 reproducer after fix; check station
files for off-hypo slip; `max‖Q‖_∞ > O(Pa)`; BP5 tests still pass;
test-suite diff against `/tmp/phase0_baseline_tests.log`.

---

### Phase 6 — Frontera Re-verification
[REVISED per R-305, R-502]

**Actions:**
1. Commit fix + tests.
2. Re-pull on Frontera.
3. Submit `tpv102_1000m_p1_1.5s_4rank_dev.sbatch`.
4. Submit `tpv102_200m_p1_1.5s_50rank_dev.sbatch` (verification, no DIAG).
5. Submit `tpv102_200m_p1_1.5s_400rank_dev.sbatch` (dispositive).
6. Verify: count of stations with `V2 > 1e-6 m/s` = 9 out of 9;
   `max‖Q‖_∞ > 1e3 Pa`.
7. **Commit DIAG sbatch variants as permanent artifacts [R-305]:**
   - `jobs/tpv102/tpv102_200m_p1_1.2s_100rank_DIAG.sbatch`
   - `jobs/tpv102/tpv102_200m_p1_1.2s_400rank_DIAG.sbatch`
   All use `SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX -DSEAS_DIAG_GHOST_EXCHANGE"`
   in the build step and `tfinal=1.2s` (post-breakaway).
8. **Remove DIAG scaffolding from the non-DIAG production build
   [Phase 6 cleanup]:** the production sbatch scripts do not use
   `SEAS_EXTRA_CPPFLAGS`, so diag code is compiled out.  Scaffolding
   remains in source for future regression use.

---

## Escalation Triggers
[REORDERED per R-306; Tandem promoted to Phase 1 Slot B per R-703]

If any phase consumes > 4 hours without progress:

1. Revert to `v8.0.0-bug-state` tag.
2. Do NOT partially revert R-701 (FaultBasis reuse),
   R-801 (BP5 convention unification), R-802 (shared-flux
   conservation).  (These R-IDs are old per-revision fix tags,
   not the round-7 R-701..R-705 findings.  Kept for reversibility.)
3. **Frontera dev-job diagnostic:** if Phase 1 Slot A was run
   without DIAG (Phase 1A guards only), submit the
   `tpv102_200m_p1_1.2s_100rank_DIAG.sbatch` with Phase 2 DIAG
   instrumentation — 1.2 hr wall-clock, definitive for
   scale-dependent bugs.
4. **Finer-grained instrumentation (after step 3 signal confirmed):**
   per-QP debug file; post-process to find first anomalous QP.
5. **Re-run Phase 1 Slot B (Tandem) if not already completed in
   Round 1.** Tandem is no longer a last-resort escalation; it is
   a Round-1 parallel-track reference per R-703.  This step fires
   ONLY if Slot B was skipped due to Tandem-build problems and
   those are now fixed.

Notes:
- The per-QP debug file (step 4) is an `fprintf` to a
  rank-separated file `diag_rank_%04d_qp_%04d.txt` at the
  breakaway window cadence; retained only when step 3 analysis
  hits a dead end.

## Anti-Patterns Banned

- Speculate without running diagnostics.
- Modify `wave_operator.inl` physics beyond Phase 1A guards and
  Phase 2 DIAG printfs.
- Bump versions based on hypotheses.
- Re-submit 400-rank dispositive without Phase 1A + Phase 2 DIAG first.
- Use `MFEM_ASSERT` for Release-mode invariants (R-102).
- Use the 2-tet fixture as a propagation reproducer (R-302, R-603).
  Fixture A alone only tests face-flux accumulation, not volume
  propagation — always run both A and B (R-603).
- Include 200m in local pre-flight loop (R-304).
- Use `tfinal = 0.3 s` for diagnostic runs (R-502) — too early.
- Use `CXXFLAGS+=...` to pass DIAG flags (R-501) — use `SEAS_EXTRA_CPPFLAGS`.
- Flag only the hypocenter DOF for DIAG (R-504) — miss propagation failure.
- Tag hypocenter DOF without a distance filter (R-602) — N-rank ×
  2 spurious DIAG streams.
- Expect off-hypo (`flt_0_3`) DIAG to show arrival at `tfinal=1.2 s`
  (R-605) — S-wave arrival is 1.3 s; either use `tfinal=1.5 s` or
  add an intermediate bulk-DOF monitor at z=-5500 m.
- Override `--dt` above `2 · dt_cfl` (per `feedback_explicit_cfl.md`);
  prior job blew up with 6.6× override.  Not a coupling bug
  hypothesis (R-601 arithmetic correction ruled that out) — this
  is a stability limit.
- Use a linear scan of `fault_interior_faces_` (R-604) — use the
  `fault_interior_face_set_` hash set for O(1) lookup.
- Rely on stderr alone to verify `[BUILD]` banner (R-606) —
  `build_info.txt` is deterministic.
- Run DIAG on fault QPs only (R-701) — always include 2–3 bulk-DOF
  propagation monitors at 1, 2, 3 km up-dip from the hypocenter.
- Run DIAG without step-frequency throttling (R-702) — every
  `[EVAL]`/`[FLUX]`/`[BULK]` block gates on `g_seas_diag_step_gate`;
  log volume at 400 ranks / tfinal=1.2 s must stay below 1 GB.
- Treat Tandem comparison as a last-resort escalation (R-703) — it
  is a Round-1 Slot B parallel-track reference.
- Flag hypo + off-hypo on the same rank at production scales
  (R-704) — MPI_Allreduce(`both_flagged`) MUST be 0; if it isn't,
  rebuild with `SEAS_DIAG_OFFHYPO_FAR=1` to use `flt_0_12`.

## Escalation Priors (round-6 updates)

R-601 (CFL-violation hypothesis) was considered and **falsified by
arithmetic**: `WaveOperator::ComputeMaxDt` uses tet inscribed
diameter, not edge length, making `dt_cfl ≈ 7.6e-4 s` on the 200 m
mesh and observed `dt ≈ 2.85e-4 s` a factor 0.38 of `dt_cfl`.  Any
future reviewer proposing CFL as a root cause must first submit
corrected arithmetic OR a source-change-driven dt_cfl increase
(new mesh, new formula, new override).  Do NOT create a Phase 0A
investigation; the Phase 0 probe is sufficient as a permanent
sanity check.

## Notes on round-5 Unreviewed Areas

1. **`fault_interior_faces_` ctor completeness** — `TEST_R007` must
   pass on baseline code BEFORE R-503 guard is added.  Phase 0
   action explicit.
2. **r4 = 5e-10 non-zero interpretation** — Phase 2 with both
   hypo and off-hypo DIAG (R-504) distinguishes: if hypo's F_h
   is correct but off-hypo's Q_plus is near-zero, the hypocenter
   IS radiating but propagation dies.  If hypo's F_h is already
   tiny, the flux injection is the issue.
3. **CFL / dt sanity [CLOSED per R-601 arithmetic correction]** —
   Phase 0 prints `dt` vs `dt_cfl` (not `h_min/cp`).  Ratio ≈
   0.38, below CFL-1; no damping-masquerade.  Keep the probe
   as a permanent sanity check.
4. **Verifier ↔ `sf_fault` classification cross-check** — the
   R-103 test `TEST_R103_AllSharedFaultFacesMapped` compares the
   R-101 verifier's face list against `ComputeSharedFaceFluxRHS`'s
   runtime `sf_fault` classification.  If the verifier counts
   a face that `sf_fault` misclassifies, the test fails.

## Notes on round-6 Unreviewed Areas

1. **`wave.ComputeMaxDt(cfl)` implementation verified (R-601 source
   read):** returns exactly `cfl · h_min / cp` with `h_min` = tet
   inscribed diameter `= 6·V / A_total` (wave_operator.inl:55-96,
   1288-1291).  The `cp` is `flux_.GetCp()`, the P-wave speed.  No
   aspect-ratio / anisotropic correction — expected to be conservative
   on extreme-aspect tets.  Not a concern for v8.0.0 but noted.
2. **Whether the 400-rank sbatch uses `--dt X` override:** the
   current `tpv102_200m_p1_1.5s_400rank_dev.sbatch` does NOT set
   `--dt` (only `--cfl 0.5`), so `dt = dt_cfl` automatic.  ✓
3. **DG CFL formula order scaling:** the driver uses `(2p+1)`
   linear; some textbooks use `(2p+1)²`.  At p=1 (TPV102) the linear
   form is fine; at p≥3 separate verification needed.  Out of v8.0.0 scope.
4. **`SEAS_EXTRA_CPPFLAGS` propagation through MFEM's build tree:**
   the existing Rev 5 fix adds the flag to the per-source-file
   compile recipe directly, bypassing MFEM's config.mk; no
   interaction with `MFEM_USE_BENCHMARKS`, `MFEM_USE_ADIOS2`, etc.
   TEST_R501 verifies propagation end-to-end.

## Notes on round-7 Unreviewed Areas

1. **`GetFaceElementTransformations` semantics for shared faces
   post-`ExchangeFaceNbrData`:** the review raised a concern that
   MFEM may start reporting `Elem2No = ghost_idx` (positive) for
   shared faces after ghost exchange, possibly leading to
   double-counting between the interior and shared branches.  Mitigation:
   the existing shared-face check at `wave_operator.inl:574` uses
   `e2 < 0 && shared_mesh_face_set_.count(f) > 0`.  **Phase 1A
   Action 1a's `MFEM_VERIFY(is_in_fault_list == is_fault)`
   cross-check (Rev 5) would FIRE if a shared fault face fell
   through to the interior branch** — it would appear in
   `fault_interior_face_set_` but have `is_fault == false`.  So
   the round-7 concern is already covered.  Verified.
2. **`h_min_` computed AFTER partitioning:** `h_min_` is computed
   locally (inscribed diameter per element) then reduced via
   `MPI_Allreduce(MPI_MIN)` across ranks (`wave_operator.inl:102-108`).
   A degenerate-thin element on a single rank would pull `h_min`
   too low, making `dt_cfl` globally conservative — a performance
   concern, not a correctness one.  Log `h_min` per rank at startup
   if suspected.
3. **`dof_data.diag_print` stability across scatter/gather:** the
   Phase 2 re-assertion invariant (`MFEM_VERIFY(dof_data[diag_hypo_dof]
   .diag_print)` at the top of every RK4 step) would catch any
   code path that rebuilt `dof_data` from scratch.  No such path
   exists in the current code (`dof_data` is allocated once at
   init and never reallocated).  The invariant is a regression
   guard against future refactors.
4. **Step-counter visibility in `Evaluate`:** R-702's `g_seas_diag_step`
   + `g_seas_diag_step_gate` are file-scope externs in
   `fault_face_flux.hpp`, defined once in any TU that includes it
   with the flag active (e.g., `fault_face_flux.cpp`).  Linker-visible;
   no API change to `Evaluate`.  Alternative (pass step via data
   struct) was considered and rejected as a larger API change for
   DIAG-only benefit.

## Phase Log (fill in as execution progresses)

### Phase 0: Freeze, Label, Baseline, Allocation, CFL
Baseline pass/fail: __ / __
TEST_R007 passes: __
Allocation SUs: __
$SCRATCH free: __
dt vs h_min/cp ratio: __
Status: __

### Phase 1: Pre-flight + Bisection + Slot B (Tandem reference)
Pre-flight outcome: __
Chosen (np, mesh): __
Startup banner: SEAS_DIAG_FAULT_FLUX = __
1-rank stations slipping: __ / 9
N-rank stations slipping: __ / 9
Escalation A invoked: __
  - 50-rank INIT pair count: __
  - 100-rank DIAG submitted: __
Slot B (Tandem reference) submitted: __
  - Tandem stations slipping: __ / 9
  - Diff vs Slot A: __ matched / __ diverged
Case (A/B/B-partial/C/D): __

### Phase 1A: Hardening
Action 0 (Makefile + banner): __
Action 1 (interior silent-weld): __
Action 1a (interior classification cross-check) [NEW]: __
Action 1b (shared silent-weld): __
Action 2 (all-components MFEM_VERIFY): __
Aborts fired: __
Baseline test diff: __

### Phase 2: Instrumentation (hypo + off-hypo + bulk monitors)
Branch: __
Banner confirms DIAG=ON: __
Hypo DOF set: __ (rank __)
Off-hypo DOF set: __ (rank __, name __ [flt_0_3 | flt_0_12])
Cross-rank check (R-704): __
Bulk monitor @ 1 km up-dip set: __ (rank __, elem __)
Bulk monitor @ 2 km up-dip set: __ (rank __, elem __)
Bulk monitor @ 3 km up-dip set: __ (rank __, elem __)
Step-gate cadence (R-702): routine=__ breakaway=__
Hypo first divergence from expected: __
Off-hypo first divergence from expected: __
Bulk monitor 1 km max|VX| @ t=1.2s: __
Bulk monitor 2 km max|VX| @ t=1.2s: __
Bulk monitor 3 km max|VX| @ t=1.2s: __
Branch selected (3A/3B/3C): __

### Phase 2B: Scale escalation (Case C only)
2B-local (4 ranks): __
2B-100rank Frontera DIAG: __
2B-400rank Frontera DIAG (if needed): __
Final diagnosis: __

### Phase 3: Localised Fix
Bug description: __
File/Line: __
Fix: __

### Phase 4: Tests Added
Integration smoke: __
TEST_R007 / R-103 map completeness: __
TEST_R501 diag flag propagates: __
TEST_R503 interior classification: __
TEST_R504 off-hypo flagged: __
Other: __

### Phase 5: Local Parallel Verification
2-rank: __ / 9 stations slip
4-rank: __ / 9
8-rank: __ / 9
BP5 regression: __

### Phase 6: Frontera Verification
4-rank: __
50-rank: __
400-rank: __ / 9 stations slip
DIAG sbatch variants committed: __
v8.0.0 closed: __

---

## Finding Index (R-001 → R-608)
[NEW per R-608 — central map of all R-IDs to review round and plan revision addressing them]

Read order: review round introduced | severity | plan revision that
addressed it | where in the current plan it is applied | short gist.

| R-ID | Round | Sev | Revision | Applied in | Gist |
|---|---|---|---|---|---|
| R-001 | 1 | CRIT | Rev 1→2 | obsolete | initial Phase 2 sign-check instrumentation plan |
| R-002 | 1 | CRIT | Rev 1→2 | Phase 3B (inlined) | Godunov analytical test inputs corrected |
| R-003 | 1 | CRIT | Rev 1→2 | Phase 3C (inlined) | T_can rotation test expected values |
| R-005 | 1 | MOD | Rev 2→3 | Phase 1A Action 2 | ghost-exchange deep-copy all components |
| R-006 | 1 | MOD | Rev 3 | Phase 4 | `DOFData` size stability |
| R-007 | 1 | CRIT | Rev 2→3 | Phase 1A Action 1 | interior silent welded-flux `MFEM_ABORT` |
| R-009 | 1 | MOD | Rev 1→2 | removed | prior Phase 3A 2-tet-only proposal |
| R-101 | 2 | CRIT | Rev 2→3 | Phase 1 pre-flight | bisection mesh-partition assumption broken |
| R-102 | 2 | CRIT | Rev 2→3 | Phase 1A Ground Rule 8 | `MFEM_VERIFY` for Release-mode invariants |
| R-103 | 2 | CRIT | Rev 2→3 | Phase 1A Action 1b | shared-branch silent welded-flux `MFEM_ABORT` |
| R-104 | 2 | MOD | Rev 2→3 | Phase 2 | DIAG scope narrowed to branch identified by Phase 1 |
| R-105 | 2 | MOD | Rev 2→3 | Phase 2 | DIAG gated on `diag_print` field |
| R-106 | 2 | MOD | Rev 2→3 | Phase 2 step 5 | `diag_print` re-assertion invariant |
| R-107 | 2 | LOW | Rev 2→3 | Phase 3A | tightened accumulation assertion |
| R-108 | 2 | LOW | Rev 2→3 | Phase 4 | integration-level test of full pipeline |
| R-201 | 3 | CRIT | Rev 3→4 | Phase 1A Action 1c (see R-503) | shared-branch `is_fault` misclassification |
| R-202 | 3 | LOW | Rev 3→5→6 | Phase 3A | fixture-h parameterization (final: R-607) |
| R-203 | 3 | MOD | Rev 3→4 | Phase 4 | parallel smoke test in Phase 4 |
| R-204 | 3 | LOW | Rev 3→4 | procedural | Phase 5 terminology |
| R-205 | 3 | LOW | Rev 3→4 | Escalation section | escalation ordering |
| R-206 | 3 | LOW | Rev 3→4 | Decision rules | consistency between Phase 1 and Phase 2B |
| R-207 | 3 | LOW | Rev 3→5 | Phase 1A Actions 1, 1b | dead `flux_.Interior(...)` removal (final: R-507) |
| R-208 | 3 | LOW | Rev 3→4 | Phase Log | Phase-log structure |
| R-301 | 4 | MOD | Rev 4→5 | Ground Rule 5 | tiered reproduction — Frontera dev as first-line |
| R-302 | 4 | MOD | Rev 4→5 | Phase 1 Escalation A-local | 2-tet fixture NOT a propagation reproducer |
| R-303 | 4 | MOD | Rev 4→5 | Phase 2B | MPI / scale escalation schema |
| R-304 | 4 | LOW | Rev 4→5 | Phase 1 pre-flight | no 200m in local pre-flight loop |
| R-305 | 4 | LOW | Rev 4→5 | Phase 6 step 7 | commit DIAG sbatch variants as permanent artifacts |
| R-306 | 4 | LOW | Rev 4→5 | Escalation Triggers | reordering per user preference |
| R-307 | 4 | MOD | Rev 4→5 | Ground Rule 5 | do not limit tests all-local if Frontera dev suffices |
| R-401 | 4 | LOW | Rev 4→5 | RESULT.txt pattern | sbatch post-run pattern |
| R-402 | 4 | LOW | Rev 4→5 | sbatch recipe | dead-rank detector for `%.3e` format |
| R-403 | 4 | LOW | Rev 4→5 | sbatch recipe | submit-dir pinning for RESULT.txt |
| R-501 | 5 | CRIT | Rev 5→6 | Phase 1A Action 0 + Ground Rule 7, 9 | `SEAS_EXTRA_CPPFLAGS` Makefile variable + banner |
| R-502 | 5 | CRIT | Rev 5 | all DIAG tfinal values | `tfinal = 1.2 s` post-breakaway, not 0.3 s |
| R-503 | 5 | CRIT | Rev 5→6 | Phase 1A Action 1a | interior `is_fault=false` misclassification guard |
| R-504 | 5 | CRIT | Rev 5→6 | Phase 2 Step 4 | flag TWO DOFs (hypo + off-hypo) for DIAG |
| R-505 | 5 | MOD | Rev 5 | Ground Rule 5 | Tier 2 defaults to 100-rank / 2-node / 1.2s |
| R-506 | 5 | MOD | Rev 5 | Phase 1 Step 5 | 0/1-3/4-9 station counting |
| R-507 | 5 | LOW | Rev 5 | Phase 1A Actions 1, 1b | dead `flux_.Interior(...)` replaced with comment |
| R-508 | 5 | LOW | Rev 5 | Phases 1A, 2, 3A/B/C, 4, 5 | inlined unchanged-from-Rev-3 text |
| R-601 | 6 | CRIT (falsified) | Rev 6 | Phase 0 CFL probe + recorded facts | CFL-violation hypothesis eliminated by arithmetic correction |
| R-602 | 6 | CRIT | Rev 6 | Phase 2 Step 4 | hypo DOF tagging distance filter |
| R-603 | 6 | CRIT | Rev 6 | Phase 3A Fixture B | propagation-corridor test (5-tet) |
| R-604 | 6 | MOD | Rev 6 | Phase 1A Action 1a ctor | O(1) `fault_interior_face_set_` lookup |
| R-605 | 6 | MOD | Rev 6 | Phase 2 expected-value table | wave-arrival timing + intermediate monitor |
| R-606 | 6 | MOD | Rev 6 | Phase 1A Action 0 | banner to `build_info.txt` + stderr |
| R-607 | 6 | LOW | Rev 6 | Phase 3A Fixture A & B | threshold parameterization by h |
| R-608 | 6 | LOW | Rev 6 | this appendix | central Finding Index |
| R-701 | 7 | CRIT | Rev 7 | Phase 2 Step 4a/4b | bulk-DOF propagation monitors (1/2/3 km up-dip) |
| R-702 | 7 | MOD | Rev 7 | Phase 2 Step 2a + Ground Rule 11 | DIAG step-frequency throttle |
| R-703 | 7 | MOD | Rev 7 | Phase 1 Slot B (new subsection) | Tandem reference baseline in Round 1 |
| R-704 | 7 | MOD | Rev 7 | Phase 2 Step 4 off-hypo selection | dual flt_0_3 / flt_0_12 + MPI cross-rank verify |
| R-705 | 7 | LOW | Rev 7 | Phase 0 CFL section | rename probe → log, remove gating language |

R-IDs not listed above were either procedural-only or superseded by
later findings; consult the corresponding REVIEW.md round for details.
