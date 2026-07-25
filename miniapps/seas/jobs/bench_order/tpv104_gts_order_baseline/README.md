# TPV104 200 m — GTS order speed baseline (MFEM vs SeisSol, p1/p2/p3)

A **cost baseline at nominally-matched order**. Six legs: MFEM `p1/p2/p3` and SeisSol
`o2/o3/o4`, on the *same* SCEC TPV104 200 m production mesh (2,464,689 tets), **global time
stepping** (no LTS), CFL 0.5, run to **t = 2.0 s**, with **all output off** on both sides.

This folder is the complete artifact set: three MFEM decks, three SeisSol decks, one run
sbatch per code, the SeisSol o2/o3 build job, and a post-processor that reads both codes' logs
and emits the fairness gates and the scoreboard.

> **Nothing here has been submitted.** Every command below requires your explicit approval.

---

## 1. What it measures

| Quantity | Definition |
|---|---|
| `T_loop` (MFEM) | Caliper region `seas::spatial_dyn::step`, **Max across ranks** |
| `T_loop` (SeisSol) | terminal `Simulation time (compute):`, **max** from its `range: [min, max]` field |
| `s / sim-s` | `T_loop / 2.0` |
| **`us-core / element-update`** | `1e6 * T_loop * ncores / (steps * 2464689)` |
| `/mode` | the above divided by modes per element (4 / 10 / 20) |

The **µs-core per element-update** is the headline. `s/sim-s` alone is illegible across orders
because dt is *not* held constant (R3).

### Why `T_loop` and not the job wall

Setup — a 126 MB `.msh` read, a serial partition, assembly — is minutes, and it is neither the
same fraction across the two codes nor across the three orders. The MFEM sbatch therefore
**exits** if `libcaliper` is not linked: `MFEM_PERF_SCOPE` compiles to nothing without it
(`mfem general/annotation.hpp:31`) and the driver has no other timer. The SeisSol sbatch scrapes
only the *terminal* `Simulation time (compute):` — never the per-epoch line, whose `split()`
double-counts that epoch's IO (`SeisSol src/Solver/Simulator.cpp:130` vs `Stopwatch.cpp:35-40`).

---

## 2. Accuracy-matched pairing (R2)

SeisSol's `ConvergenceOrder` counts **basis functions**; MFEM's `[mesh].order` counts
**polynomial degree**. The matched pairing is `p_k <-> o_(k+1)`:

| leg | MFEM `[mesh].order` | MFEM `ader_order` | SeisSol binary | modes/elem | dt [s] | steps to 2.0 s | element-updates |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 2 | `seissol-elastic-o2-f64` | 4 | 8.542987e-4 | 2342 | 5.7723e9 |
| 2 | 2 | 3 | `seissol-elastic-o3-f64` | 10 | 5.125792e-4 | 3902 | 9.6172e9 |
| 3 | 3 | 4 | `seissol-elastic-o4-f64` | 20 | 3.661280e-4 | 5463 | 1.34646e10 |

At a matched pair the two codes advance **state vectors of identical size** (`4/10/20` modes ×
9 elastic quantities) on identical geometry, with the **same dt** and the **same step count**,
and both perform `o` friction sub-steps per macro step. That is the strongest fairness fact
this baseline has, and it should be stated before any timing number.

**dt equality is algebraic, not coincidental.** MFEM:
`dt = cfl * h_min / (cfl_dg_safety * (2p+1) * cp)` with `h_min = 6V/A_total`
(`dynamic/wave_operator.inl:7728` × `spatial/code/spatial_friction.hpp:744`). SeisSol:
`dt = cfl * 2*r_insphere / ((2o-1) * vp)` (`GlobalTimestep.cpp:45-46`). With `o = p+1` the order
factors coincide and `2r == 6V/A` identically; at `cfl = 0.5` **and `cfl_dg_safety = 1.0`** the
two codes take the same dt. The dt values above are *pre-registered predictions*, scaled exactly
by `7/(2p+1)` from the verified p3/o4 anchor `3.66128e-4 s` — both sbatches gate each leg on
them (0.5 % relative tolerance) and refuse to publish a leg that disagrees.

---

## 3. Fairness rules and how each is honoured

**R1 — same mesh, same window, same CFL, GTS both.**
`.msh` md5-pinned `e88ae223205530548b1d3438e6eb17c1`; `.puml.h5` md5-pinned
`5911b0b4ec29a8f61968d5b32ed0a11c` (the pumgen conversion of that same file). Both sbatches
refuse to run on a md5 mismatch. `tfinal = "2.0s"` / `EndTime = 2.0`, `cfl = 0.5` / `CFL = 0.5`,
`lts = "off"` / `ClusteredLTS = 1` — all asserted per leg before launch.

Two keys silently rescale dt and are therefore **pinned in the deck *and* gated**:
`cfl_dg_safety = 1.0` (struct default is **3.0** → dt/3, 3× the steps, MFEM inflated 3× while
every log line still looks sane) and `time_integrator = "ader"` (the RK branch uses
`RkCflFactor = 3/(2p+1)` → dt×3, and `cfl_dg_safety = 1` still prints).

**R2 — accuracy-matched pairing.** Nothing in the MFEM driver ties `ader_order` to `mesh.order`,
and `--ader-order` on the CLI *overrides the deck*. The sbatch derives `ADER = ORD+1` per leg,
asserts the deck against it, and passes that same integer.

**R3 — dt is not held equal across orders.** Each order gets its own CFL dt. Both `s/sim-s` and
`µs-core/element-update` are reported, plus `/mode`, and the within-code order ratio is split
into the exactly-known `(2p+1)` step penalty and the measured per-step work.

**R4 — IO equal-or-off on both sides.** MFEM has *no* compute/IO split, so IO cannot be
subtracted after the fact — it must be off. Four things had to be turned off, not three:

- `paraview_fault = "off"`;
- `paraview_free_surface = "off"` — **the trap**: the free-surface slice is *not* under the
  `paraview_enabled` master gate (`spatial_dyn_driver.cpp:3865-3871` vs the gate at
  `:3628-3632`), defaults to `"vtu"` @ 0.05 s, is written from inside the timed scope, and the
  driver prints `ParaView output: OFF` *while it writes*. No log gate can catch this; only the
  deck gate can;
- `[problem].tag = ""` — a `"tpv104"` tag wires the SCEC station writer, which fflushes 9 trace
  files **every step inside the timed scope** on ≤9 of 256 ranks, with a row count that grows
  with order (2342/3902/5463) — while the SeisSol legs write *nothing* inside their compute
  stopwatch (`OutputPointType = 0`, `ReceiverOutput = 0`);
- `checkpoint_every_steps = 100000` (> every leg's step count).

On the SeisSol side, `EnergyOutput = 0` is **not** an IO knob: it flips `isFrictionEnergyRequired`
and adds per-GP work *inside* the DR friction kernel. `SurfaceOutput = 0` also kills the PGV
fold, which sampled every local timestep inside the compute region. `ReceiverOutput` and
`Checkpoint` both **default to true** and are set explicitly. With all periodic output off there
is exactly **one** simulation epoch, so no sync point clips dt (the old 0.05 s cadence added ~40
truncated steps, +0.7 %, purely from IO).

*One write remains and is expected:* MFEM always writes a post-loop final checkpoint
(`:5265-5274`; the parser forbids `checkpoint_every_steps = 0`) — 256 ASCII files,
~0.7/1.8/3.6 GB at p1/p2/p3. It is **outside** the Caliper scope. The R4 audit therefore expects
*exactly* `nranks` `cp_checkpoint_r*.txt` and flags only an excess (an earlier version flagged
every leg, which would have trained the reader to ignore the check).

**R5 — each code's validated rank shape, unchanged.** MFEM 2 nodes × 128 pure-MPI ranks × 1 core
= 256 working cores. SeisSol 2 nodes × 8 ranks × 16 cores, `OMP_NUM_THREADS = 15`,
`SEISSOL_COMMTHREAD = 0` — 256 allocated, **240 active**. This is why the two codes need
separate jobs (rank shape is allocation-time, and SeisSol's prologue derives `OMP_NUM_THREADS`
from `$SLURM_CPUS_PER_TASK`). The 240/256 split is **finding M4** — disclosed, not fixed.

**R6 — per leg: wall, steps, dt, elements, µs-core/update.** All recorded, with provenance. Step
counts are taken from Caliper `region.count` (MFEM) and derived from the **measured** dt
(SeisSol) — never from the pre-registered value, so a failed dt gate cannot produce a
plausible-looking metric against a fictional step count.

### Statistic matching

MFEM's Caliper wall is a **max** over 256 ranks; SeisSol's `Simulation time (compute):` **leads
with the mean** over 16 ranks. Quoting one against the other credits SeisSol with its ~6.4 %
load imbalance. Both sbatches and the post-processor emit **max-vs-max** (headline — a
bulk-synchronous step ends when the last rank arrives) **and mean-vs-mean**.

### Core denominator

Publish **both**, never a single unqualified number:

- `C_alloc = 256` for both codes — billing-fair headline;
- `C_active = 256` MFEM / `240` SeisSol — bounds the idle-core effect (6.7 %).

Never normalise by *ranks*: 256 MFEM ranks vs 16 SeisSol ranks is a 16× error.

---

## 4. What this baseline CANNOT claim

1. **No accuracy claim of any kind.** There is no error functional in this experiment. Nominal
   order equality is not error equality — the two codes differ in flux implementation,
   quadrature, DR sub-stepping and nucleation. Do **not** say "time to solution at matched
   accuracy" or "p3 is worth it". This is a *cost* baseline at *nominally*-matched order.
2. **No scaling claim.** One mesh, one node count, one problem, three order points. No
   strong/weak scaling, no mesh-convergence study.
3. **It does not generalise to LTS.** The existing SeisSol LTS numbers (239 s as-run, 111 s
   compute) are a *different experiment* and are additionally polluted — that run had
   `EnergyOutput = 1`, PGV sampling, per-step pickpoints and 0.005 s receivers, all inside the
   compute stopwatch. The new GTS compute figure is **not** comparable to the historical 111.
4. **The MFEM number carries disclosed, unremoved overheads** that SeisSol does not pay:
   - `pmesh.SetCurvature([mesh].order)` is applied unconditionally to a *straight-sided* tet mesh
     (`spatial_dyn_driver.cpp:1733`), so MFEM evaluates a 4/10/20-node H1 nodal transformation
     per quadrature point where SeisSol is always affine. **Part of the p-scaling and of the p3
     gap is a redundant geometry order, not DG arithmetic.** `--face-cache` removes it on
     non-fault interior faces only.
   - A per-step `Q.Norml2()` + `MPI_Allreduce` NaN tripwire across all 256 ranks
     (`:5190-5199`), not disableable, with no SeisSol counterpart.
   - ~17 nested Caliper regions inside the timed scope vs SeisSol's bare start/pause stopwatch.
     (`profile.mpi` is deliberately omitted from `CALI_CONFIG` for this reason; the residual
     annotation tax is still one-sided.)
   - The order-dependent R-101 shared-fault check is removed uniformly by `SEAS_R101_SKIP=1`.
5. **Partition-count confound.** MFEM is decomposed 256 ways, SeisSol 16 — roughly **2.5× more
   partition-interface area for MFEM**, all of it MPI, worst at p1. Per-*thread* work is matched
   to 6.7 % (9,628 elem/rank vs 10,270 elem/thread); only the halo differs. Any "MFEM improves
   relative to SeisSol as order rises" is partly a halo artifact. A shape-matched control
   (SeisSol 256 pure-MPI ranks) would convert this objection into a measurement — not run here.
6. **The p1/o2 pair is the least reliable point.** At 4 modes SeisSol pads vector lanes, so its
   HW-FLOP ≫ NZ-FLOP there. Read the per-leg NZ-vs-HW ratio; do not lead with p1.
7. **The legs produce no physics artifacts** (that is the point of IO-off). A leg that diverged,
   or a p1 that failed to nucleate at 200 m, would still yield a plausible throughput number.
   Cheap sanity signals are recorded — MFEM `V_max` peak (p3 reference ≈ 17.3 m/s at t ≈ 1.48 s),
   SeisSol `Total calculated HW-FLOP` and `Load imbalance` — but **validate rupture in a separate
   IO-on run**, never in a timing job.
8. **Single run per leg.** No repeats, no steady-state isolation. Repeat at least the two p3 legs
   and report the spread; if the spread exceeds the cross-code gap, the baseline is noise-limited
   and must say so.

---

## 5. Prerequisite: the SeisSol o2/o3 binaries

**SeisSol's order is compile-time.** `o4` already exists; **`o2` and `o3` do not and must be
built.** `build_seissol_o2_o3_expanse.sbatch` clones the proven per-order recipe and changes
only `ORDER`, then asserts the result:

- source commit pinned to `dc6db65133a999eb38c578c14cc443df1000dfd0` — checked against the
  working tree **and against every binary's embedded SHA, including the pre-existing o4** (a
  clean `git rev-parse HEAD` proves nothing about a binary compiled weeks ago);
- `CMakeCache.txt` of each new build diffed against o4's, normalised — **any** non-`ORDER`
  difference exits 8 (the diff is the artifact that keeps this an order study, not a build study);
- `ORDER:STRING` asserted per build dir (exit 7);
- codegen toolchain pinned **absolutely** (numpy 2.5.1 / PSpaMM 0.3.1) *and* before-vs-after — a
  delta check alone cannot see drift that already happened, which would give o2/o3 different GEMM
  microkernels from o4;
- `PRECISION=double` passed explicitly (the recipe reads it from the environment, and
  `--export=ALL` would carry a stale `PRECISION=single` straight through);
- success is decided by **artifact checks**, never by the recipe's exit status — that script ends
  in `echo` and returns 0 even when `cmake --build` fails.

Exit codes: `0` ok · `2` root/recipe missing · `3` commit mismatch (tree or binary) · `4` exe
missing · `5` stale binary · `6` no codegen python · `7` cache ORDER mismatch · `8` cache differs
beyond ORDER · `9` toolchain drift. **First-writer-wins**, so the most fairness-relevant cause
survives to the exit status.

Cost: ~8 min/order on `shared` at 16 cores; ~8–32 SU.

---

## 6. Expected cost and walltime

Projected from the measured p3 GTS anchor (2.366 s/step, job 52344266) and the measured
p1:p2:p3 = 1 : 2.74 : 7.45 per-step ratio:

| leg | steps | ≈ s/step | stepping wall | ≈ SU (256 cores) |
|---|---|---|---|---|
| MFEM p1 | 2342 | 0.32 | ~12 min | ~55 |
| MFEM p2 | 3902 | 0.87 | ~57 min | ~240 |
| MFEM p3 | 5463 | 2.37 | ~3.6 h | ~920 |
| SeisSol o2 | 2342 | — | ~4 min | ~20 |
| SeisSol o3 | 3902 | — | ~10 min | ~45 |
| SeisSol o4 | 5463 | — | ~25 min | ~110 |

Plus one mesh read + partition per leg. **Total ≈ 1.4 k SU.**

Two scheduling notes:

- The p3 anchor was measured over `t = 0 → 1.0 s`, which is the entirely **pre-nucleation**
  second (`T_nuc = 1.0 s`). The 2.0 s window is half dynamic rupture, where the per-QP rate-state
  solves are more expensive — so treat the p3 projection as a **lower** bound. The default lever
  set includes `--face-cache` (4930 s/sim-s, not 6465), which pulls the other way.
- **A walltime kill produces no Caliper report at all** (it is written at clean exit), so an
  overrun loses the leg entirely, not just its tail. **Split the MFEM legs into three jobs** so a
  p3 overrun cannot destroy p1/p2.

**Account:** both codes are set to `usc143` (~291 k SU free). `lbl107` had only ~4.7 k SU left
on 2026-07-24, expiring 2026-08-27 — under 4× margin for this baseline, and running the two
halves of one experiment on different projects is an avoidable queue-priority asymmetry.
`ddp408` (~719 k free) is the alternative: `sbatch -A ddp408 ...`.

---

## 7. Submission sequence

### Step 0 — build the MFEM driver **with Caliper** (if not already)

```bash
ssh czhao1@login.expanse.sdsc.edu
cd <repo>/miniapps/seas/jobs/comm_dev/a0_wiggle_tpv104_200m_expanse
sbatch build_driver_expanse.sbatch          # hard-fails on a missing libcaliper
```

### Step 1 — build SeisSol o2 and o3 (~30 min, `shared`)

```bash
scp <repo>/miniapps/seas/jobs/bench_order/tpv104_gts_order_baseline/seissol/build_seissol_o2_o3_expanse.sbatch \
    czhao1@login.expanse.sdsc.edu:/expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/seisol_build/

ssh czhao1@login.expanse.sdsc.edu \
  'cd /expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/seisol_build && sbatch build_seissol_o2_o3_expanse.sbatch'

# watch:  tail -f .../seisol_build/build-o2o3-elastic.<jobid>.out
# GATE:   the job must exit 0.  Any of 3/5/7/8/9 voids the cross-order comparison.
```

### Step 2 — MFEM legs (three jobs, cheapest first)

```bash
ssh czhao1@login.expanse.sdsc.edu
cd <repo>/miniapps/seas
# confirm the 126 MB mesh is staged (it is *.msh = gitignored, it does NOT arrive via git):
ls -l tpv104/mesh/tpv104_200m.msh

B=jobs/bench_order/tpv104_gts_order_baseline
S=/expanse/lustre/scratch/$USER/temp_project/tpv104_gts     # final checkpoints land here

sbatch --export=ALL,BENCH_LEGS=1,BENCH_OUT_BASE=$S -t 1:00:00 $B/run_mfem_gts_orders.sbatch
sbatch --export=ALL,BENCH_LEGS=2,BENCH_OUT_BASE=$S -t 2:00:00 $B/run_mfem_gts_orders.sbatch
sbatch --export=ALL,BENCH_LEGS=3,BENCH_OUT_BASE=$S -t 8:00:00 $B/run_mfem_gts_orders.sbatch
```

(All three in one 8 h job is `sbatch $B/run_mfem_gts_orders.sbatch` — but then a p3 overrun
destroys p1 and p2 as well.)

### Step 3 — SeisSol legs (one job, all three orders)

Depends on step 1. If step 1 is still queued, chain it:
`sbatch --dependency=afterok:<build_jobid> ...`

```bash
# stage a FRESH directory — never reuse an existing SeisSol run dir
D=/expanse/lustre/projects/lbl107/czhao1/safs/tpv104_gts_order
rsync -avP <repo>/miniapps/seas/jobs/bench_order/tpv104_gts_order_baseline/seissol/ \
      czhao1@login.expanse.sdsc.edu:$D/
rsync -avP ~/Downloads/tpv104/tpv104.puml.h5 czhao1@login.expanse.sdsc.edu:$D/

ssh czhao1@login.expanse.sdsc.edu "cd $D && sbatch run_seissol_gts_orders.sbatch"
# single leg:  sbatch --export=ALL,BENCH_SS_LEGS=\"2\" run_seissol_gts_orders.sbatch
```

The `seissol/` folder is self-contained apart from the 103 MiB mesh — the two yaml and two dat
companions ship with it, so one rsync stages a complete run directory. Each leg runs in its own
`leg_o<N>.<jobid>/` with inputs symlinked, so legs cannot collide and no stop-file can silently
short-circuit a rerun.

### Step 4 — gather and post-process

```bash
# bring the SeisSol logs next to the MFEM ones
rsync -avP czhao1@login.expanse.sdsc.edu:$D/run_o*.log  <mfem OUT_ROOT>/

python3 <repo>/miniapps/seas/jobs/bench_order/tpv104_gts_order_baseline/postprocess_baseline.py \
        <mfem OUT_ROOT> --json baseline.json
```

---

## 8. How to read the result

The post-processor prints eight sections. Read them **in order** — the gates come before the
numbers on purpose.

- **§0 Pre-registered contract** — E, modes, fault QPs, ADER sub-steps, dt, steps, updates. All
  derived from the single p3 anchor, so the three orders cannot drift apart. Printed *before*
  any timing.
- **§1 Leg inventory** — found/missing, with the path a missing leg would have lived at.
- **§2b Cross-code dt/step equality** — the cheapest and most important fairness check. MFEM
  `[time] dt` vs SeisSol `Minimum timestep:` at each matched order. **If these disagree, stop.**
- **§2 Fairness gates** — 13 per code. Any FAIL voids that leg. Watch especially G5
  (`cfl_dg_safety == 1.0`), G8b (free-surface slice off), G9 (no in-loop checkpoint, now gated on
  the deck value because the GTS in-loop write is *silent* in the log) and G9b (`tag` empty).
- **§3 Baseline table + per-leg provenance** — where each `T_loop` came from, which statistic,
  which step source, which core denominator, and the rank imbalance.
- **§4 Cross-code ratio** — max-vs-max headline plus mean-vs-mean companion, with a
  self-consistency warning: with equal dt, steps, elements and cores the wall-ratio and the
  per-update ratio **must** coincide; a gap means one of those equalities silently broke.
- **§5 Within-code scaling** — with the exactly-known `(2p+1)` step penalty divided out. The
  `(2p+1)*modes` line (12 : 50 : 140) is printed as a **crude linear-in-DOF model, not a
  prediction** — ADER volume kernels are closer to `O(modes²)`, so a measured p3/p1 **above**
  11.67× is expected and is *not* a defect.
- **§6 SeisSol-only cross-checks** — `Load imbalance`, HW/NZ-FLOP, and the LoopStatistics
  `( per element )` regression slope. That slope is a per-rank-team (15-thread), kernel-only
  quantity: use it as an internal consistency check on SeisSol's own wall, **never** compared
  directly to an MFEM number.
- **§7 What could not be computed**, then the one-line caveat.

Exit codes: `0` clean · `1` a leg or a required R6 quantity is missing · `2` usage · `3` a gate
FAILed (or WARNed under `--strict`).

---

## 9. Files

| file | role |
|---|---|
| `tpv104_gts_p{1,2,3}.toml` | MFEM decks. **Pairwise diff is exactly 3 lines** (`order`, `ader_order`, `output_dir`) — order is the only MFEM-side variable |
| `mesh_tpv104_200m.msh` | symlink → `../../../tpv104/mesh/tpv104_200m.msh` (gitignored; the sbatch re-creates the link and hard-fails if the target is absent) |
| `run_mfem_gts_orders.sbatch` | MFEM p1/p2/p3, 2 nodes × 128 pure MPI |
| `seissol/parameters_gts_o{2,3,4}.par` | SeisSol decks. **Bodies byte-identical** (`md5 fc41f51f…`); they differ only in the `!` header, which names each leg's binary and its pre-registered dt |
| `seissol/run_seissol_gts_orders.sbatch` | SeisSol o2/o3/o4, 2 nodes × 8 ranks × 16 cores |
| `seissol/build_seissol_o2_o3_expanse.sbatch` | builds the two missing orders — **the only copy with the CMakeCache fairness diff**; run it from `$SEISSOL_ROOT` |
| `seissol/tpv104_{fault,material}.yaml`, `tpv104_{receivers,faultreceivers}.dat` | SeisSol companions (the `.dat` files are inert under IO-off but are still required by the pre-flight) |
| `postprocess_baseline.py` | stdlib-only; reads both codes' logs, applies the gates, prints the scoreboard |

Verify the deck contracts at any time:

```bash
diff tpv104_gts_p1.toml tpv104_gts_p2.toml      # must be exactly 3 lines
diff <(grep -v '^!' seissol/parameters_gts_o2.par) \
     <(grep -v '^!' seissol/parameters_gts_o4.par)   # must be empty
```
