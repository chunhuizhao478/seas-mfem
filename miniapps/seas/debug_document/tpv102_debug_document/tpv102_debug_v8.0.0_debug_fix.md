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

### Step 1.8 — Job 7666233: tolerance fix did NOT resolve the abort

**Result.**  Re-run against commit `9cfd296` (abs_floor=1e-6 +
%.17e diagnostic).  Same count: `n_entries=246, n_pairs=122,
n_unpaired=2`.  The %.17e diagnostic now reveals the actual orphan
coordinates:

```
[UNPAIRED] centroid=(-9.89933571681227295e+03, 2.60237444819293559e-13, -4.25000000000785622e+03) rank=3 group_size=1
[UNPAIRED] centroid=(-9.89933571681227113e+03, 2.60237444819293559e-13, -4.25000000000785622e+03) rank=6 group_size=1
```

Drift in x: `|(-9899.33571681227295) - (-9899.33571681227113)| =
1.82e-12`.  At magnitude 9899, that's **exactly 1 ULP** of
double precision.  y and z agree to bit-identical 17-digit
precision (drift == 0).  With `abs_floor=1e-6`, `same_centroid(A,B)`
**must** return TRUE (1.82e-12 ≪ 1e-6).

So our diagnosis from Step 1.7 (FP drift > 1e-9 breaks the
tolerance) is **incomplete**.  Something else is causing the two
entries to be assigned to separate groups of size 1 each.  Possibilities:

- (A) same_centroid returns FALSE despite coords well within tol —
  compilation / tolerance-evaluation subtlety we haven't
  understood.
- (B) another entry sorts lex-between A and B and has
  same_centroid(A, intruder) == FALSE — a sort-ordering pathology
  (unlikely given cx values are 1 ULP apart, leaving no strict-
  between double).
- (C) the grouping loop itself has a bug — e.g., the post-group
  `i = j` accidentally skips one of the pair.

None of these can be distinguished from the current log.

### Step 1.9 — Diagnostic self-check added (for job 1 re-run #4)

**Rationale.**  Submitting another Frontera run on a guess would
waste SUs.  Added a `[self-check]` block to the R-101 verifier
that runs exactly once (rank 0 only) when `unpaired.size() == 2`
and prints:

- `same_centroid(orphan_0, orphan_1) = TRUE|FALSE` — directly
  answers whether the tolerance check passes between the two
  orphans.
- `dx / dy / dz` — per-coordinate drift between them.
- `sort positions: orphan_0 at P0, orphan_1 at P1 (distance D)`
  — reveals whether they were adjacent in the sort order.
- If `D != 1`, enumerate up to 16 entries between them with
  full-precision coords + rank — shows the "intruder" records that
  broke the grouping.

**Outcome table** (what to conclude after job 1 re-run #4):

| `same_centroid` | `distance` | Diagnosis | Fix direction |
|-----------------|------------|-----------|---------------|
| FALSE | 1 | tol not applied (impossible given source?) | recheck build, `abs_floor` value at runtime |
| TRUE  | 1 | algorithm bug in grouping loop | audit `i = j` / `group_size==2` branch |
| TRUE  | >1 | sort non-determinism or NaN-like coord | audit sort lambda, check coord values |
| FALSE | >1 | intruder entry has |coord diff| > 1e-6 | intruder coords will show why; may need vertex-key pairing |

**Local verification:** `mpirun -np 2 ./seas_test_r101_shared_fault`
→ 18/18 PASS.  Self-check is inert when `unpaired.size() != 2`
and diagnostic-only when triggered.

### Step 1.10 — BP5-aligned integer-key pairing (root-cause fix)

**Comparison with BP5 (user directive):** Why does BP5 not see the same
failure?  BP5's `ElasticityOperator` (`domain/elasticity_operator_setup.inl:917-994`)
pairs shared fault faces across ranks via an integer triple of sorted
global vertex IDs (`MakeFaceKey` → `key_ranks` map).  TPV102's R-101
verifier was instead pairing by physical centroid with a FP tolerance.
The centroid approach is the source of every pathology we chased
through Steps 1.4-1.9: sub-ULP MFEM drift, ordering intruders,
tolerance tuning.  The integer-key approach is bit-exact and immune.

**Fix applied** (user directive: "reuse BP5 pattern, do not touch BP5
source, write a new file with same logic"):

1. NEW file `miniapps/seas/dynamic/shared_fault_key.hpp` — duplicates
   BP5's `ElasticityOperator::FaceVertexKey` / `MakeFaceKey` as a
   free-standing struct + inline function, bit-identical logic.  Header
   comment documents the unification plan: later BP5 will include this
   header and drop the nested copies, collapsing the duplication.

2. Rewrote `VerifySharedFaultDOFDataConsistency` in
   `dynamic/wave_operator.inl`:

   - Record now 16 doubles: `[key(3 int64 as double), qp_idx, centroid(3),
     fields(8), rank]`.
   - Sort + group by `(face_key, qp_idx)` as integers — exact, no
     tolerance.  `same_face_qp` replaces `same_centroid`.
   - Self-check block and abs_floor infrastructure removed; they
     existed only to work around the centroid-based pairing, which is
     gone.
   - Diagnostic `[UNPAIRED]` now prints face_key + qp_idx + centroid
     so any remaining orphan is traceable to a specific mesh face.

3. Zero modification of BP5 source (`domain/*`, `fault/*`, `solver/*`
   unchanged).

**Local verification:**

```
mpirun -np 2 ./seas_test_r101_shared_fault       → 18/18 PASS
mpirun -np 2 ./seas_test_parallel_wave_operator  →  5/5 PASS
```

The 1000m np=2 test continues to cover the paired path.  To exercise
the failure configuration that no unit test covers (50-rank production
topology), we ran `seas_tpv102_driver` locally:

```
mpirun --oversubscribe -np 50 ./seas_tpv102_driver \
  --mesh tpv102/mesh/tpv102_1000m.msh --mesh-scale 1 --order 1 \
  --bc-mode absorbing --cfl 0.5 --tfinal 0.0005 \
  --output-dir /tmp/tpv102_test --no-domain-pv
```

Result: **pairing succeeds** (no "unpaired entries" abort).  The
verifier now completes the pairing phase and enters the field-comparison
phase — the scenario that previously could not be reached.

### Step 1.11 — NEW finding exposed: real psi drift between rank pairs

Once pairing works, the verifier reports a **genuine** inter-rank field
inconsistency:

```
MFEM abort: R-101 shared-fault DOFData consistency FAILED.
  Field 'psi' at centroid (-1.5272e+04, 9.8e-13, -1.6038e+04)
  differs by 4.277e-3 across the two ranks sharing the face
  (scale=1.0, rel_diff=4.277e-3, rel_tol=1e-10).
```

Reproduces at np=14 (psi drift 2.19e-2), np=30 (3.12e-2), np=50
(4.28e-3) on the 1000m mesh — different faces/magnitudes, but always
the same class of defect: `psi` on the **same** shared fault QP differs
between its two rank owners after a single RK4 step.

This is NOT a new bug our fix introduced.  It is the pre-existing
v8.0.0 defect that the broken pairing was masking — the centroid
matcher failed with "unpaired entries" BEFORE the field-comparison
loop could run.  With integer-exact pairing, we now see the real
divergence.

The abort's own hint is on point: "check that the R-501 owner-broadcast
in ComputeSharedFaceFluxRHS covers every mutable field written by
FaultFaceFlux::Evaluate and the driver's RK4 averaging step".  psi
(state variable) is likely not being propagated from the owner rank
to the non-owner after the RK4 stage that updates it.

**Local reproducer established — no Frontera needed for Phase 2:**
The 1000m / np=14 case aborts in under 30 seconds on a laptop and
produces the full R-101 diagnostic (face_key + centroid + field name +
drift magnitude + rank pair).  Phase 2 of v8.0.0 (diagnose & fix the
psi drift) can proceed entirely against this local reproducer.

### Step 1.12 — Phase 2 diagnosis: "psi drift" was a verifier artifact, not a real bug

**Instrumentation:** Added `SEAS_DIAG_FAULT_FLUX` block to
`ComputeSharedFaceFluxRHS` (wave_operator.inl) that prints Evaluate
inputs (Q_plus, Q_minus hashes, psi, slip_rate, V1/V2, tau/sigma
corrections) and outputs per RK4 stage, filtered to a target
centroid.

**Ran np=14 reproducer with DIAG enabled.** Data from rank 5 and
rank 6 (the two ranks sharing the failing fault face at
`centroid=(-1.67317e+04, ~0, -1.58660e+04)`):

| Stage | psi (both ranks) | V_abs (both ranks) | h_plus / h_minus | Agreement |
|-------|------------------|---------------------|-------------------|-----------|
|  k1   | 8.221292303436545e-1 | 9.9999999999729e-13 | identical   | bit-exact |
|  k2   | 8.221292303436544e-1 | 9.9999999999731e-13 | identical   | bit-exact |
|  k3   | 8.221292303436544e-1 | 9.9999999999731e-13 | identical   | bit-exact |
|  k4   | 8.221292303436541e-1 | 9.9999999999731e-13 | identical   | bit-exact |

**Finding:** at the physically-same QP on both ranks, the
fault-state fields (psi, V_abs, V1, V2, tau_corr, sigma_n_corr) are
**bit-identical** across all 4 RK4 stages.  No drift exists.

Yet the R-101 verifier was aborting with `psi differs by 2.19e-2`.
The centroid reported by the verifier matched our DIAG target.  So
where did 2.19e-2 come from?

**Root cause (the actual one):** MFEM's shared-face orientation
differs between the two ranks.  DIAG showed:

- Rank 5: `sf=1624 q=2 cz=-15866.037060 sign_flipped=1`
- Rank 6: `sf=1262 q=0 cz=-15866.037078 sign_flipped=0`

These are the **same physical QP**, but each rank indexes it with a
different local `qp_idx` (rank 5 calls it q=2, rank 6 calls it q=0).
The Step 1.10 pairing algorithm was matching by `(face_key, qp_idx)`
— which for this face pairs rank 5's q=0 with rank 6's q=0, a
DIFFERENT physical point ~face_size/3 away.  Comparing psi at two
different physical points gave the 2.19e-2 "drift" — purely the
spatial gradient of psi across the face, not any inter-rank
inconsistency.

**Phase 2 fix (minimal, same file):**
`VerifySharedFaultDOFDataConsistency` now sorts by `(face_key
primary, physical centroid secondary)` and groups entries that
match on BOTH face_key (integer exact) AND centroid within 1e-6 m
tolerance.  Within a face, QPs are separated by ~face_size/3
(hundreds of meters on the production mesh), so the 1e-6 m floor
safely distinguishes different physical QPs on the same face while
matching the ~1 ULP FP drift between ranks' views of the same QP.

The qp_idx field stays in the record (diagnostic) but no longer
participates in grouping.  Face key remains integer-exact.

**Local verification after the Phase 2 fix:**

```
mpirun -np 2 ./seas_test_r101_shared_fault       → 18/18 PASS
mpirun -np 2 ./seas_test_parallel_wave_operator  →  5/5 PASS
```

Production-configuration runs:

| np | mesh    | result |
|----|---------|--------|
| 14 | 1000m   | R-101 OK: 27 pairs, 0 unpaired, max_rel_diff=1.11e-16 (1 ULP) |
| 30 | 1000m   | R-101 OK: 21 pairs, 0 unpaired, max_rel_diff=4.44e-16 (4 ULP) |
| 50 | 1000m   | R-101 OK: 27 pairs, 0 unpaired, max_rel_diff=1.11e-16 (1 ULP) |

All three report **bit-exact** DOFData consistency across ranks at
step 0.  The DIAG instrumentation was removed after confirming the
root cause — it remains accessible via `SEAS_DIAG_FAULT_FLUX` build
flag for future debugging.

**Propagation sanity check (np=14 / 1000m / tfinal=0.05 s, 35
steps):**

```
Step  0/35, V_max = 1.00000e-12 m/s, qnorm r0=0        (distant rank quiet)
Step 24/35, V_max = 1.00007e-12 m/s, qnorm r0=1.2e-27  (wave reached r0)
Step 34/35, V_max = 1.00248e-12 m/s, qnorm r0=7.6e-27  (wave still growing)
R-101 verifier at step 0: OK (unchanged)
```

V_max increases monotonically.  All 14 ranks have nonzero ||Q||
by step 12.  No "pinned V_max, silent stations" signature at
tfinal=0.05 s.  Breakaway (V_max > 0.1 m/s) is expected at t ≈ 1 s;
this short local run is below breakaway by design.

### Step 1.13 — Phase 1+2 complete; proposing Frontera production confirmation

The original v8.0.0 signature from Frontera job 7665297 (400-rank,
tfinal=1.5 s) was "V_max pinned at 7.699 m/s at hypocenter, off-hypo
stations at V_ini=1e-12, bulk max‖Q‖_∞ ≈ 1.66e-8".  Our theory of the
failure cascade was:

1. R-101 shared-fault DOFData drift → (hypothesized, NOW DISPROVEN).
2. Drift causes asymmetric fault flux across partition seams.
3. Seam asymmetry blocks wave propagation to off-hypo stations.

Phase 1+2 locally disprove step (1): DOFData is bit-exact across
ranks at fault QPs, including through the full RK4 step.  The
original R-101 aborts at np=50 / 200m / Frontera were
verifier-pairing artifacts (centroid tolerance + qp_idx misalignment
under MFEM face orientation flip), NOT real field drift.

**Remaining unknowns that only Frontera can confirm:**

- Does the 200m production mesh at 50/100/400 ranks pair cleanly
  with the Phase 2 verifier?  (Local 1000m mesh says yes; 200m has
  different METIS partition topology.)
- Does V_max reach breakaway at tfinal ≈ 1 s as expected by TPV102
  physics, or does the original "pinned V_max" signature persist?
  If it persists, the cause is **not** inter-rank fault-state drift
  (ruled out) and we move to Phase 3 (bulk wave propagation / ghost
  Q exchange / CFL / nucleation physics).

**Proposed Frontera jobs (USER APPROVAL REQUIRED per
`feedback_frontera_approval.md`):**

1. `jobs/tpv102/tpv102_200m_p1_0.01s_50rank_init.sbatch` — init-only
   sanity at the original failing configuration.  ~5 min dev-queue,
   ~0.1 SU.  Primary check: R-101 passes (Phase 1+2 fix holds on the
   200m mesh + Frontera Intel-MPI topology).

2. `jobs/tpv102/tpv102_200m_p1_1.2s_100rank_bisect.sbatch` — bisection
   tier 2 at 100 ranks, full tfinal=1.2s.  ~1.2 hr dev-queue, ~100 SU.
   Confirms Phase 2 fix holds past step 0, tests breakaway.

3. `jobs/tpv102/tpv102_200m_p1_1.2s_400rank_bisect.sbatch` — bisection
   tier 3 at 400 ranks, full tfinal=1.2s.  ~1.5 hr dev-queue, ~400 SU.
   This is the job 7665297 dispositive reproducer — if V_max pinning
   still happens, we know the bug is scale-dependent and NOT in R-101
   pairing.  If it resolves, the v8.0.0 bug was entirely the verifier
   artifact cascade we fixed in Phase 1+2.

Total budget: ~500 SU.  No autonomous submission — awaiting explicit
user approval.

**Phase 1+2 Phase Log (final):**
```
[x] Phase 1A Action 0 (Makefile + banner) applied + verified
[x] Partition pre-flight (local) complete
[x] Frontera sbatches built (4 files)
[x] Job 1 submitted (50-rank init) — R-101 "2 unpaired" signature
[x] Steps 1.4-1.9: iterative diagnosis of R-101 centroid pairing
[x] Step 1.10: BP5-pattern integer face_key pairing via
    `shared_fault_key.hpp` (no BP5 modification)
[x] Step 1.11: discovered apparent "psi drift" after pairing fixed
[x] Step 1.12: DIAG-instrumented np=14 reproducer; confirmed
    DOFData is bit-exact across ranks.  The "drift" was a second
    verifier bug: pairing by qp_idx mismatches under MFEM face
    orientation flip.  Added centroid tiebreaker within face group.
[x] Local verification: R-101 PASS at np=2/14/30/50.  DIAG confirms
    bit-exact fault-state consistency through all 4 RK4 stages.
[x] Propagation sanity (np=14 / 0.05 s): waves reach all ranks,
    V_max growing, no pinning signature.
[ ] Frontera confirmation (3 jobs, ~500 SU) — AWAITING USER APPROVAL
[ ] If all 3 PASS: v8.0.0 closed; mark signature resolved by Phase 1+2
[ ] If Frontera 400r still shows V_max pinning: open Phase 3 (bulk
    propagation / ghost Q exchange), with confirmed rule-out that
    fault-state drift is not the cause.
```
