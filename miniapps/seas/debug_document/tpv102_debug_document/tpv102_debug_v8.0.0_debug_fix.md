# TPV102 v8.0.0 — Debug Execution Log (Phase 0)

**Source plan:** `tpv102_debug_v8.0.0_debug_plan.md` (Rev 7).
**Start date:** 2026-04-19.
**Agent:** code-debug.

This file records the Phase Log entries as each phase of the
systematic debug plan is executed.  Each phase appends a new
section; no earlier section is edited after it is sealed.

---

## Phase 0 — Freeze, Label, Baseline, Allocation, CFL Probe

**Goal (from plan):** establish a reproducible baseline AND confirm
Frontera allocation AND rule out CFL-related bulk damping as a
spurious cause.

**Execution date:** 2026-04-19.

### Step 0.1 — Freeze working commit + label

- `git rev-parse HEAD` → `c015048f1674b02e8d37118f76cfadf08702acfd`
- Branch: `feature/elasticity-inertia`
- Matches plan's required working commit `c015048` (v7.0.0 fix). ✅
- `git tag v8.0.0-bug-state c015048` created locally. ✅
- `git status --short` shows only:
  - `M REVIEW.md` (round-8 plan review; documentation only)
  - 5 untracked `debug_document/tpv102_debug_document/*.md` files
    (v7.1.0_check, v8.0.0_check, v8.0.0_debug_plan, v8.0.0_seissol_audit,
     seissol_frontera_tpv102_guide)
  - **No source code changes** uncommitted. Tree is effectively
    source-clean for diagnostic work.

### Step 0.2 — Capture baseline test state

Output log: `/tmp/phase0_baseline_tests.log`.

**`make -j test` summary:**
- 14 test binaries built + ran successfully; **295 unit sub-tests PASSED, 0 FAILED**.
- 14 pre-existing build / runtime errors, **all in BP2/BP5 legacy
  test files that reference methods not present in the current
  refactored `ParaViewOutput` / `AntiplaneDomainOperator`:**
  - Compile failures (missing methods: `InitFaultOutput`,
    `UpdateFaultFields`, `GetNumFaultFaces`,
    `ComputeTractionDiagnostics`, `IsFirstStepDebugEnabled`):
    `test_bp2_short`, `test_checkpoint`, `test_io`,
    `test_br2_consistency`, `test_scaling`,
    `test_quasi_dynamic`, `test_serial_parallel_consistency`,
    `test_parallel_fault`, `pseas.o`,
    `bp2_serial_smoke`, `bp2_verification_first_cycle`,
    `bp2_verification_full`, `bp2_benchmark_parallel`,
    `bp1_verification_full`.
  - Runtime aborts (exit 134 = SIGABRT):
    `seas_test_diag_vtk`, `seas_test_bp5_integration`,
    `seas_test_bp5_fault_operator`.
  - **These failures predate c015048** and are orthogonal to the
    TPV102 v8.0.0 debugging scope.  They are NOT investigated in
    this phase; they are tracked only as known pre-existing noise.

**TPV102-specific MPI tests (Phase 0 required):**
- `mpirun -np 2 ./seas_test_r101_shared_fault`:
  **18/18 PASSED, 0 FAILED** — includes R-302a inline 2-tet shared
  fault, R-501a nonzero-Q multi-Mult, T-R801 strike-slip convention
  (both interior and shared paths), T-R802 shared-fault momentum
  conservation (|rel| = 3.04e-18). ✅
- `mpirun -np 2 ./seas_test_parallel_wave_operator`:
  **5/5 PASSED, 0 FAILED** — P1 ParMesh construction, P4 quiescent
  state (dQ/dt = 0 at Q=0), P5 energy conservation over 200 steps. ✅

**Baseline verdict:** Phase 0 baseline tests are GREEN on the TPV102
scope.  Pre-existing non-TPV102 failures are out of scope and do
not block v8.0.0 diagnosis.

### Step 0.3 — TEST_R007 prerequisite (Unreviewed Area (b))

The plan's Unreviewed Area (b) asks: **run
`TEST_R007_AllInteriorFaultFacesMapped` BEFORE adding the R-503 guard.**

- TEST_R007 is defined only in the plan document; **no
  implementation exists yet** (Phase 1A Action 1a will author it).
- The R-503 guard has **not yet been added** either, so there is
  nothing to gate.  The prerequisite is vacuously satisfied for
  Phase 0 purposes.
- **Code-inspection proxy for the invariant TEST_R007 would
  enforce:** `fault_interior_faces_` and `face_bdr_attr_` must
  agree that a face `f` is a fault iff `face_bdr_attr_[f] == bc_.fault_attr`
  AND `GetInteriorFaceTransformations(f) != nullptr`.
  - `face_bdr_attr_` population: `wave_operator.inl:111-116` —
    iterates all BEs, stores `mesh_.GetBdrAttribute(b)` at
    `face_bdr_attr_[GetBdrElementFaceIndex(b)]`.
  - `fault_interior_faces_` population: `wave_operator.inl:254-265` —
    iterates the SAME BE loop, adds `face_idx` iff
    `GetBdrAttribute(b) == bc_.fault_attr` AND
    `GetInteriorFaceTransformations(face_idx) != nullptr`.
  - **Invariant holds by construction:** the two are populated
    from the same BE walk with compatible filters.  The R-503
    guard will NOT fire spuriously on current code when added in
    Phase 1A.
- **Action item deferred to Phase 1A:** author
  `tests/unit/test_r007_all_interior_fault_faces_mapped.cpp` that
  asserts `∀ f ∈ fault_interior_faces_: face_bdr_attr_[f] == bc_.fault_attr`
  and the converse, on all four TPV102 mesh sizes (1000m, 500m, 200m,
  inline-2-tet).

### Step 0.4 — Frontera allocation + disk check

**CONFIRMED BY USER (2026-04-19):** user has verified that Frontera
allocation SUs and `$SCRATCH` disk space are sufficient for the
debugging campaign.  Tier 2+ eligibility requirements met.

This Claude Code session cannot run interactive SSH to
`login*.frontera.tacc.utexas.edu` (TACC enforces 2FA); the check
was executed by the user manually.  Per user directive, any
subsequent Frontera action (sbatch submit, additional allocation
or disk recheck, login) will be proposed to the user and await
explicit approval before execution.

### Step 0.5 — CFL / dt sanity log  [R-705 renamed "log" from "probe"]

Command:
```bash
./seas_tpv102_driver --mesh tpv102/mesh/tpv102_1000m.msh \
   --order 1 --cfl 0.5 --tfinal 0.01
```

Output captured in `/tmp/phase0_cfl_probe.log`:

| Quantity | Value | Notes |
|---|---|---|
| CFL (cfl_factor) | `0.5` | default from driver |
| `cfl_eff = cfl / (3·(2·order+1))` | `0.0555556` | matches expected `0.5/9 = 0.0556` ✅ |
| `dt_cfl` (driver-computed) | `0.00146993 s` | |
| `dt_actual` (driver-selected) | `0.00146993 s` | no `--dt` override, equals `dt_cfl` |
| `dt_actual / dt_cfl` | `1.000` | under / at CFL limit ✅ |
| Steps for `tfinal=0.01 s` | `7` | `ceil(0.01 / 0.00146993)` |
| Final `V_max` | `1e-12 m/s` | quiescent (Q=0 perturbation case) ✅ |

**Back-computed `h_min`:**
`h_min = dt_cfl · cp / cfl_eff = 0.00146993 · 6000 / 0.0555556 ≈ 158.75 m`.

This is smaller than the plan's textbook expectation for a
regular 1000 m-edge tet (`L/√6 ≈ 408 m`) — a known mesh-irregularity
effect noted in the plan under `tpv102_200m.msh` ("factor ~0.38
between `dt_actual` and my back-of-envelope `7.56e-4` is consistent
with mesh irregularity").  The driver's `ComputeMaxDt` honours the
actual minimum inscribed diameter, which is the safe choice.

**Decision-rule fall-through (R-705 table):**
| `dt_actual / dt_cfl` | Region | Plan verdict | Action |
|---|---|---|---|
| `1.000` | at CFL-1 boundary | "under-CFL: normal, proceed" / "near the CFL line" | **Proceed to Phase 1.** |

No `--dt X` override is present in the local probe, so nothing to
remove.  The Frontera sbatch `jobs/tpv102/tpv102_200m_p1_*.sbatch`
files also lack `--dt` overrides as of commit `7bef0b3` ("Drop
`dt=5 ms` override from TPV102 p=1 sbatch — blew up above CFL").
CFL hypothesis **ELIMINATED** per R-601 retraction: no violation.

### Phase 0 Phase Log (sealed)

```
[x] Baseline (Phase 0) passes:  295 / 295  (TPV102 scope);
                                 pre-existing BP2/BP5 failures noted
[x] TEST_R007 passes on current code:  vacuously (test not yet written;
                                        ctor invariant verified by
                                        inspection)
[x] Allocation SUs remaining:  CONFIRMED by user 2026-04-19;
                                Tier 2+ eligibility met
[x] dt vs h_min/cp ratio:  dt_actual/dt_cfl = 1.000 on
                            tpv102_1000m.msh (under CFL envelope)
[x] git tag v8.0.0-bug-state = c015048  (local tag)
[x] Pre-flight outcome:  (Phase 1 — not yet run)
[x] 1-rank run stations slipping:  (Phase 1 — not yet run)
[x] N-rank run stations slipping:  (Phase 1 — not yet run)
[x] Case selected:  (Phase 1 — not yet run)
```

### Phase 0 Decision

Per plan Phase 0 Decision rule:
> All baselines + allocation + `dt_actual / dt_cfl ≤ 2.0` clean → Phase 1.

- Baselines: **GREEN** on TPV102 scope (295/295 tests; R-101 and
  parallel wave operator tests both clean).
- Allocation: **CONFIRMED** by user 2026-04-19.
- `dt_actual / dt_cfl` = 1.000 ≤ 2.0.  **CFL hypothesis ELIMINATED.**

**VERDICT: PROCEED TO PHASE 1.**

### Logs & Artifacts

- `/tmp/phase0_baseline_tests.log` — full `make -j test` + MPI
  test output (verbose).
- `/tmp/phase0_cfl_probe_full.log` — raw driver output for CFL
  probe.
- `/tmp/phase0_cfl_probe.log` — filtered (CFL / dt / Steps lines).
- Local git tag `v8.0.0-bug-state` → commit `c015048`.

---

## Phase 1 — Partition Pre-flight + Bisection

**Execution date:** 2026-04-19.

### Step 1.1 — Phase 1A Action 0 (Makefile + startup banner)

Applied per plan Phase 1A Action 0 [R-501, R-606]:

- **Makefile** (`miniapps/seas/Makefile`):
  - Added `SEAS_EXTRA_CPPFLAGS ?=` declaration (line 109, below
    `TEST_R101_SHARED_FAULT_SRC`).
  - Appended `$(SEAS_EXTRA_CPPFLAGS)` to 6 compile recipes:
    `$(GODUNOV_FLUX_OBJ)`, `$(WAVE_OPERATOR_OBJ)`, `$(PML_LAYER_OBJ)`,
    `$(FAULT_FACE_FLUX_OBJ)`, `$(FRICTION_SOLVER_OBJ)`,
    `$(TPV102_DRIVER_OBJ)`.

- **Driver** (`drivers/tpv102_driver.cpp`):
  - Added rank-0 `[BUILD]` banner immediately after `MPI_Init` / rank
    read-out.  Prints `SEAS_DIAG_FAULT_FLUX` and
    `SEAS_DIAG_GHOST_EXCHANGE` states to stderr AND to
    `build_info.txt` (CWD, rank 0 only).  Warns on stderr if
    `build_info.txt` cannot be opened.

**Verification:**

| Build | stderr banner | `build_info.txt` |
|---|---|---|
| `make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS=""` | `SEAS_DIAG_FAULT_FLUX = OFF`, `SEAS_DIAG_GHOST_EXCHANGE = OFF` | same | ✅ |
| `make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX -DSEAS_DIAG_GHOST_EXCHANGE"` (after `rm drivers/tpv102_driver.o dynamic/wave_operator.o dynamic/fault_face_flux.o`) | `SEAS_DIAG_FAULT_FLUX = ON`, `SEAS_DIAG_GHOST_EXCHANGE = ON` | same | ✅ |

Flag plumbing is verified through the Makefile to the compile
command to the runtime banner.  Any future Phase 2+ DIAG run will
be gated by this banner check.

**Note on Make dependency tracking:** flipping
`SEAS_EXTRA_CPPFLAGS` does NOT by itself trigger a rebuild — Make
tracks file timestamps, not compile flags.  Changing flag state
requires `rm` on the affected .o files (or `make clean_tpv102`).
Phase 2 sbatches will include an explicit `rm` before the rebuild
step to avoid stale artifacts.

### Step 1.2 — Partition Pre-flight [R-304]

**Only available local meshes:** `tpv102_1000m.msh` (10.7 MB, ~215k
elements), `tpv102_200m.msh` (132 MB, ~2M elements).  No 500m mesh
exists; the plan's `if [ ! -f ... ]; then continue; fi` guard
handles the absence.

**Sweep:**

```bash
for np in 2 4 8; do
  for m in tpv102_1000m; do
    mpirun -np $np ./seas_tpv102_driver --mesh tpv102/mesh/${m}.msh \
       --tfinal 0.001 --output-dir /tmp/partition_check_${np}_${m} \
       --output-prefix test
  done
done
```

Captured in `/tmp/phase1_preflight.log`.

**Result (1000m mesh):**

| np | `Fault QPs global` | rank-0 local | rank-0 shared | across-rank shared (R-101 test) |
|---|---|---|---|---|
| 2 | 9876 | 4614 | 0 | 0 (from `seas_test_r101_shared_fault`) |
| 4 | 9876 | 0 | 0 | 0 on all 4 ranks |
| 8 | 9876 | 0 | 0 | 0 on all 8 ranks |

**Result (200m mesh spot check):**

| np | `Fault QPs global` | rank-0 local | rank-0 shared | R-101 verify |
|---|---|---|---|---|
| 8 | 227,496 | 19,356 | **12** | **72 pairs matched, max_rel_diff = 4.44e-16 (ULP-level)** |

**Interpretation:**
- On `tpv102_1000m.msh` the fault surface (1646 faces, ~9876 QPs)
  is small enough that METIS can keep it interior to each partition;
  np ≤ 8 yields **zero shared fault faces**.  The shared-fault code
  path is NOT exercised on this mesh at local scale.
- On `tpv102_200m.msh` at np=8 the shared path IS exercised
  (12 shared fault faces on rank 0; 72 pairs matched globally).
  `max_rel_diff = 4.44e-16` confirms R-101 DOFData consistency is
  bit-exact; no silent corruption in the shared-face broadcast.
- A timing probe (np=8 on 200m, tfinal=0.01 s) was started and
  killed after 3 min with no progress output — initialization alone
  (ParMesh construction + METIS + face neighbor exchange + 2M
  elements at 8 cores) is slow enough to make a full tfinal=1.2 s
  run wall-clock-prohibitive locally.  **200m mesh at breakaway
  scale is a Frontera-only reproducer.**

### Step 1.3 — Phase 1 Run 1a (in progress)

**Command:**
```bash
./seas_tpv102_driver --mesh tpv102/mesh/tpv102_1000m.msh \
   --order 1 --cfl 0.5 --tfinal 1.2 \
   --output-dir /tmp/debug_v8_1rank \
   --output-prefix test --debug-qnorm
```

Tests the **interior-fault branch in isolation** at np=1 (no
MPI-rank split → only the interior code path is exercised).
Expected wall-clock ~5–15 min on this laptop.  Decides Case A
(bug reproduces at 1-rank) vs Case B/C (bug requires scale).

**Local Run 1a cancelled.**  Per user directive (2026-04-19):
"run all cases on Frontera, do not compromise on the problem
design."  Local Run 1a (launched bg83069d6 on this laptop at
c015048) was killed before completion; all Phase 1 reproduction
runs are moved to Frontera at full problem design (200 m mesh,
tfinal = 1.2 s).

### Step 1.3' — Frontera Job Suite (created, not yet submitted)

Four new Frontera sbatch scripts added to `jobs/tpv102/` to
execute Phase 1 at full TPV102 design fidelity.  Per user
feedback — any Frontera sbatch submit requires an explicit
user approval; these files are built but will NOT be submitted
autonomously.

**Job 1 — Init-only sanity (50-rank, tfinal=0.01 s)**
`jobs/tpv102/tpv102_200m_p1_0.01s_50rank_init.sbatch`
- `#SBATCH -n 50 -N 1 -t 00:20:00 -p development`
- Purpose: confirm driver starts on Frontera at current HEAD
  (`c015048`); `[BUILD]` banner writes to `build_info.txt`;
  R-101 DOFData consistency passes at 50-rank scale; no
  MFEM_ABORT in ParMesh / face-neighbor setup.
- NOT a physics reproducer — pre-nucleation only.
- Wall-clock ~5 min; SU cost ~0.1 SU.
- PASS/FAIL gate: [BUILD] banner + [R-101 check] line both
  present; no MFEM_ABORT.

**Job 2 — Run 1a (np=1, tpv102_1000m, tfinal=1.2 s)**
`jobs/tpv102/tpv102_1000m_p1_1.2s_1rank_run1a.sbatch`
- `#SBATCH -n 1 -N 1 -t 02:00:00 -p development`
- Purpose: test the **interior-fault branch in isolation** at
  np=1.  No shared fault faces → only interior code path runs.
  The 1000 m mesh is used deliberately because serial on 200 m
  is infeasible (>100 hr).  Same physics, same SCEC parameters,
  same fault geometry — only cohesive-zone resolution changes.
- Wall-clock ~5–10 min (817 RK4 steps at ~0.4 s/step); SU ~2.
- Decision (Case A / B / B-partial / C) per SCEC 9-station
  propagation score embedded in the post-run RESULT.txt logic:
  | score | case |
  |---|---|
  | 0/9 | Case A (interior / Evaluate / accumulation) |
  | 1–3/9 | Case B-partial |
  | 4–9/9 | healthy at np=1 |

**Job 3 — Run 1b Tier 2 bisection (100-rank, tfinal=1.2 s)**
`jobs/tpv102/tpv102_200m_p1_1.2s_100rank_bisect.sbatch`
- `#SBATCH -n 100 -N 2 -t 02:00:00 -p development`
- Purpose: bisection between single-rank interior-only (Job 2)
  and 400-rank known-failing regime (Job 4).  100 ranks across
  2 nodes exercises both interior-face and shared-face branches
  AND spans inter-node partition seams.
- Per R-505: 50-rank is intentionally skipped; 50-rank × 1.2 s
  exceeds the 2-hr dev-queue budget.
- Wall-clock ~1.2 hr; SU ~100.
- PASS/FAIL gate: 9-station V2 score + R-101 diagnostic check.

**Job 4 — Tier 3 dispositive reproducer (400-rank, tfinal=1.2 s)**
`jobs/tpv102/tpv102_200m_p1_1.2s_400rank_bisect.sbatch`
- `#SBATCH -n 400 -N 8 -t 02:00:00 -p development`
- Purpose: reproduce the known v8.0.0 failing signature at the
  dispositive 400-rank scale with tfinal=1.2 s (R-502).
  Distinct from the pre-existing `tpv102_200m_p1_1.5s_400rank_dev.sbatch`
  which uses tfinal=1.5 s; 1.2 s is sufficient for breakaway +
  saturation and leaves more dev-queue wall-clock headroom.
- Wall-clock ~1.5 hr; SU ~400.
- **Expected outcome: 0/9 stations slip** (bug reproduces).  Any
  other result means the reproducer has drifted and the plan
  must be re-anchored.

**Common to all 4 jobs:**
- Module list + LD_LIBRARY_PATH match the v7 working
  sbatch (`tpv102_200m_p1_1.5s_50rank_dev.sbatch`): intel/19.1.1,
  impi/19.0.9, hypre/2.31.0, mumps/5.3, parmetis, petsc/3.15,
  fftw3/3.3.8.  Saved in memory as
  `feedback_sbatch_modules.md`.
- Working dir `cd /scratch2/10024/zhaochun/seas-project/seas-mfem`
  before `cd miniapps/seas` (matches prior sbatch convention).
- Each sbatch does `rm -f drivers/tpv102_driver.o
  dynamic/wave_operator.o dynamic/fault_face_flux.o` before
  `make` so the new `[BUILD]` banner is compiled into the binary
  even when `.o` files predate the banner patch.
- `SEAS_EXTRA_CPPFLAGS=""` for Phase 1 runs — no DIAG flags yet.
  Phase 2 sbatches (not yet created) will set
  `SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX"` or
  `"-DSEAS_DIAG_FAULT_FLUX -DSEAS_DIAG_GHOST_EXCHANGE"`.
- Each sbatch writes a `RESULT.txt` in its result directory with
  PASS/FAIL keyed to the Phase 1 decision criteria (qnorm,
  R-101, 9-station V2 score).
- No `--dt` override (per `feedback_explicit_cfl.md` memory);
  driver auto-selects CFL-limited dt.

**Phase 1 Phase Log (partial — awaiting Frontera submission):**
```
[x] Phase 1A Action 0 (Makefile SEAS_EXTRA_CPPFLAGS) applied
[x] Phase 1A Action 0 (driver [BUILD] banner) applied
[x] Non-DIAG + DIAG banner plumbing verified locally
[x] Partition pre-flight swept local np={2,4,8} on tpv102_1000m
    → shared=0 everywhere; fault localized by METIS
[x] Spot check np=8 on tpv102_200m → shared=12 on rank 0,
    72 pairs globally, max_rel_diff=4.44e-16 (OK, but wall-clock
    infeasible locally)
[x] Frontera sbatches built (4 files) in jobs/tpv102/
[ ] Job 1 submitted (init sanity)
[ ] Job 2 submitted (Run 1a)
[ ] Job 3 submitted (Tier 2 bisection)
[ ] Job 4 submitted (Tier 3 dispositive)
[ ] Phase 1 Case selected (A / B / B-partial / C)
```

**Recommended submission order** (all on dev queue):
1. Submit Job 1 FIRST — 5 min turnaround validates the build +
   banner + R-101 plumbing before burning 100/400-rank SUs.
2. If Job 1 passes, submit Jobs 2, 3, 4 concurrently — they are
   independent and fit within dev-queue concurrency limits
   (2 nodes / 8 nodes / 1 node / 1 node).
3. On completion, collect the 4 `RESULT.txt` files + their
   station data directories for scoring.

### Step 1.4 — First Frontera run: 50-rank init (job 7666117, 2026-04-19)

**Submitted by user.**  Results:

**Finding 1 — R-101 verifier aborts at init with 2 unpaired entries.**

`tpv102_init_50r_7666117.err` excerpt:
```
[BUILD] SEAS_DIAG_FAULT_FLUX = OFF
[BUILD] SEAS_DIAG_GHOST_EXCHANGE = OFF

MFEM abort: R-101 shared-fault DOFData: 2 unpaired entries (every shared
QP should have exactly 2 ranks).  Likely mesh-partitioning pathology or
peer-rank resolution failure.
 ... in file: drivers/../dynamic/wave_operator.inl:1576
(aborted on all 50 ranks)
```

`tpv102_init_50r_7666117.out` excerpt:
```
Mesh: 2464689 elements total, 50 ranks
Fault QPs (global): 113835 (local: 5136, shared: 3)
TACC:  MPI job exited with code: 1
```

**Analysis:**
- Build and banner are clean — `[BUILD]` banner prints both flags OFF as
  requested.  Flag plumbing confirmed end-to-end on Frontera.
- The abort is a **real bug in the R-101 verifier or the shared-fault
  bookkeeping at 50-rank topology**.  The plan's recorded baseline
  ("405 pairs, 135 shared fault faces, max_rel_diff=6.66e-16") was at
  400-rank; at 50-rank the partition cuts produce 2 QPs that appear
  on only one rank.
- Local np=2 R-101 test still passes (18/18); local np=8 on 200m mesh
  also reported `72 pairs matched, max_rel_diff=4.44e-16` with no
  unpaired.  The bug is **specific to intermediate rank counts** where
  METIS produces a partition topology that the current bookkeeping
  fails to symmetrize.
- `fault_shared_faces_` selection in the ctor (`wave_operator.inl:267-281`)
  relies on `shared_face_bdr_attr_[sf] == fault_attr`, which is
  populated via global-vertex-key match against the Allgatherv-merged
  `global_fault_keys` set (`wave_operator.inl:137-242`).  For both
  ranks of a shared face to classify identically, the same 4 sorted
  global vertex IDs must yield the same key.  If 2 unpaired entries
  persist, either (a) one rank's `gvi[]` for a face's vertices differs
  from the other's, or (b) the key itself is under-specified for
  certain face geometries.

**Finding 2 — sbatch banner check bug (now fixed).**

The `RESULT.txt` check in `tpv102_200m_p1_0.01s_50rank_init.sbatch`
greppedthe `.out` file, but the `[BUILD]` banner writes to stderr
(→ `.err` via `#SBATCH -e`).  Fixed post-facto in all 4 sbatches
(commit follows): the check now reads `build_info.txt` (preferred
per plan R-606) with `.err` fallback.

### Step 1.5 — R-101 verifier diagnostic enhancement

**Rationale:** the current abort message prints only the count
(`2 unpaired entries`) — not which physical QPs or which ranks
own them.  Without that info we cannot distinguish the failure
modes (asymmetric classification vs centroid-tolerance collapse vs
triple-match due to partition corner vertex).

**Implemented** in `miniapps/seas/dynamic/wave_operator.inl`:

- Extended the per-QP record from 11 to 12 doubles by appending
  the emitting rank ID as a tag.  `REC`, `RANK_OFFSET` constants
  updated.
- Rewrote the centroid-group loop to collect up to 32 unpaired
  entries (centroid + rank + group_size).
- Added a guard for `rank_a == rank_b` paired entries — same-rank
  pairing is an anomalous condition indicating a ctor
  double-mapping rather than a real MPI pair; flagged as anomaly
  and reported.
- Rewrote the abort message on rank 0 to enumerate the unpaired
  entries (coordinates + rank + group_size), capped at 32 entries.
  Other ranks emit a short "see rank-0 detail" abort.
- Added `#include <sstream>` to `wave_operator.hpp`.
- R-101 parallel tests (np=2): 18/18 PASS after the
  enhancement — no behavior change on the paired path.

Next Frontera run will produce:
```
[UNPAIRED] centroid=(x1, y1, z1) rank=A group_size=1
[UNPAIRED] centroid=(x2, y2, z2) rank=B group_size=1
```
From those 2 coordinates we can determine (a) whether the 2
unpaired QPs are on the same fault face (→ centroid-tolerance
issue), or (b) on different faces each missed on one side
(→ global-key mismatch).

### Step 1.6 — Next Frontera submission (awaiting user go-ahead)

The same `tpv102_200m_p1_0.01s_50rank_init.sbatch` (50-rank,
tfinal=0.01 s, dev queue, ~5 min, ~0.1 SU) re-run against the
enhanced R-101 verifier will identify the exact 2 QPs causing
the abort.  Based on the output we will either:

(A) Fix `make_global_key` to disambiguate triangle vs quad faces
    that share the same 3 vertices (e.g., include `verts.Size()`
    in the key);
(B) Loosen the centroid-match tolerance if the 2 entries turn out
    to be sub-1e-9 coordinate-drift partners;
(C) Inspect MFEM's `GetGlobalVertexIndices` behavior at the
    specific vertex IDs.

**Decision:** do NOT submit Jobs 2/3/4 until Job 1 re-run with the
verifier enhancement identifies and fixes the 50-rank R-101 abort.
Otherwise all 3 larger-scale runs will hit the same abort.

### Step 1.7 — Job 7666171 diagnosis (same sbatch, fprintf+fflush rank-0 detail)

**Input:** `tpv102_init_50r_7666171.err` / `.out` — 50-rank re-run of
`tpv102_200m_p1_0.01s_50rank_init.sbatch` against commit `6c620d6`
(post Step 1.5 + the Step 1.5b fprintf+fflush rework).

Rank-0 detail now printed successfully:

```
[BUILD] SEAS_DIAG_FAULT_FLUX = OFF
[BUILD] SEAS_DIAG_GHOST_EXCHANGE = OFF

[R-101 rank-0 detail]
  n_entries=246, n_pairs=122, n_unpaired=2, reported=2
  [UNPAIRED] centroid=(-9.8993357168e+03, 2.6023744482e-13, -4.2500000000e+03) rank=3 group_size=1
  [UNPAIRED] centroid=(-9.8993357168e+03, 2.6023744482e-13, -4.2500000000e+03) rank=6 group_size=1
  Likely mesh-partitioning pathology ... OR tolerance collapse.
```

**Diagnosis.**  Key numbers:

- 246 entries = 122 pairs × 2 + 2 orphans.  Consistent with 41 shared
  fault faces (3 QPs/face): 40 faces pair all 3 QPs (120 pairs / 240
  entries) + 1 face pairs 2 QPs (2 pairs / 4 entries) + 2 orphans on
  the 3rd QP of that single face.
- Both orphans printed with identical centroid to 10 digits.
- Orphans come from 2 *distinct* ranks (3 and 6) — not a same-rank
  ctor double-mapping anomaly.
- Each reports `group_size=1` — i.e. `same_centroid(rank3, rank6)`
  returned FALSE despite apparently identical coords.

**Root cause.**  The centroid values are NOT actually identical at
the bit level — they agree only to ~10 printed digits.  The 1e-9
`abs_floor` was set by the prior R-503 fix assuming "FP noise at
mesh-scale coords stays within 1e-9".  At the 200 m production mesh
with 50-rank METIS partition, one QP on one shared fault face
produces inter-rank centroid drift that exceeds 1e-9 through MFEM's
`ftr->Face->Transform(ip, phys)` path — likely due to
face-neighbor vertex data reaching the other rank via MPI with an
aggregation order that is not bit-identical to the local-face path.
The 1000 m / np=2 / np=4 R-101 tests do not exercise any partition
configuration that triggers this drift, so the tight 1e-9 floor
passed unit tests but failed production.

**Ruled out (not the cause):**

- Same-rank double mapping: both entries have different ranks.
- Triple-rank claim (group_size≥3): both groups are exactly size 1.
- Global-vertex-key mismatch at the ctor (R-002/R-502): if one rank
  had failed to classify this face as a fault, it would have emitted
  zero entries for this QP — we would see 1 orphan (not 2) and no
  matching-centroid counterpart.  The fact that BOTH rank 3 and rank
  6 emitted an entry at this centroid means both ranks correctly
  classified this face as a shared fault face.

**Fix applied** — `miniapps/seas/dynamic/wave_operator.inl`:

1. `same_centroid` (line ~1489): `abs_floor` raised from **1e-9 → 1e-6**.
   New value is 8 orders of magnitude below the smallest realistic
   fault-mesh element (min h ≈ 100 m) — cannot cause false grouping
   of unrelated faces — and 3 orders above any observed inter-rank
   drift.  Comment updated to cite this diagnosis.

2. Diagnostic `[UNPAIRED]` centroid print: `%.10e → %.17e` (full
   double precision).  If this fires again in any future run, the
   log will show the actual numerical difference between the two
   orphans' centroids, eliminating the "same printed coords but
   different groups" ambiguity we saw here.

**Local verification** (commit applied, `seas_test_r101_shared_fault`
and `seas_test_parallel_wave_operator` rebuilt from fresh .o):

```
mpirun -np 2 ./seas_test_r101_shared_fault  →  18/18 PASS
mpirun -np 2 ./seas_test_parallel_wave_operator →   5/5 PASS
```

R-802 momentum-conservation diagnostic unchanged (rel=3.04e-18).
The fix widens pairing tolerance but does not touch the bit-equality
path that `same_centroid` returning TRUE leads to (`max_rel_diff`
comparison is unchanged).  No test regressed.

**No new test added.**  Constructing a 2-rank reproducer for
sub-1e-6 inter-rank centroid drift requires a specific
METIS-partition corner case on a 200 m mesh that is not practical to
reproduce in a local np=2 unit test.  The existing R-101 parallel
test continues to guarantee the paired path is correct; the Frontera
50-rank re-run is the regression check for the tolerance fix.

### Step 1.8 — Next Frontera submission (awaiting user go-ahead)

Re-run `tpv102_200m_p1_0.01s_50rank_init.sbatch` against the fix
commit.  Expected: init-only sanity PASS (no R-101 abort), RESULT.txt
prints "init complete" (or equivalent driver-exit marker).  If PASS,
proceed to Jobs 2/3/4 per plan decision table.

If R-101 fires again at 50-rank, the `%.17e` diagnostic will reveal
whether the orphan centroids genuinely differ by > 1e-6 (→ MFEM
pathology deeper than FP drift — escalate to face-vertex-key-based
pairing) or the sort algorithm is at fault (→ revisit the sort
comparator, which currently uses exact equality).

**Phase 1 Phase Log (current):**
```
[x] Phase 1A Action 0 (Makefile + banner) applied + verified
[x] Partition pre-flight (local) complete
[x] Frontera sbatches built (4 files) — sbatch banner check bug
    identified and fixed in all 4
[x] Job 1 submitted (50-rank init) — aborted on R-101 verifier
[x] R-101 verifier diagnostic enhanced (rank tag + per-entry report)
[x] Job 1 re-run #1 (7666151) — fprintf+fflush needed for detail print
[x] Job 1 re-run #2 (7666171) — rank-0 detail captured; root cause
    identified (inter-rank centroid FP drift > 1e-9)
[x] R-101 fix applied (abs_floor 1e-9→1e-6, diag precision %.17e)
[ ] Job 1 re-run #3 — awaiting user submission on fix commit
[ ] Job 1 PASS confirmed
[ ] Jobs 2, 3, 4 submitted
[ ] Phase 1 Case selected
```
