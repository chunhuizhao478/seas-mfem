# SAFS CFL A/B (Lever A) — SeisSol-equivalent time step

**Goal.** Test whether taking the SeisSol-equivalent ADER step (dt ×3) changes the
SAFS science. SAFS carries an extra ×3 DG safety factor beyond the mandatory
Cockburn–Shu order factor `1/(2N+1)`; SeisSol does not. Lever A exposes that
factor as `[numerics].cfl_dg_safety` (default **3.0** = today, byte-exact;
**1.0** = SeisSol-equivalent). The Courant number `cfl` stays **0.5 (≤ 1)** in
both arms — we do **not** raise it. The ×3 comes only from `cfl_dg_safety`, via
the new `--cfl-dg-safety` flag.

## The two arms
Same everything (config, mesh, 2N×200r, upwind+ADER, opt stack, tfinal, snapshot
schedule); the **only** difference is `cfl_dg_safety`:

| Arm | `cfl_dg_safety` | dt vs baseline | Meaning |
|-----|-----------------|----------------|---------|
| **A** | 3.0 | ×1.0 | today (byte-exact reference) |
| **B** | 1.0 | ×3.0 | **SeisSol-equivalent** |

Script: `spatial_dyn_cflAB_ratestate_v3alt_vw005_cvm_upwind_ader_triq_expanse_2N200r_compute_2hr_safs.sbatch`

## Submit (on Expanse — do NOT run the production mesh locally)
```bash
cd <root>/miniapps/seas/jobs/safs/safs_expanse
SAFS_CFL_DG_SAFETY=3.0 sbatch spatial_dyn_cflAB_ratestate_v3alt_vw005_cvm_upwind_ader_triq_expanse_2N200r_compute_2hr_safs.sbatch   # arm A
SAFS_CFL_DG_SAFETY=1.0 sbatch spatial_dyn_cflAB_ratestate_v3alt_vw005_cvm_upwind_ader_triq_expanse_2N200r_compute_2hr_safs.sbatch   # arm B
```
Outputs/logs are job-ID + `cfldg{3p0,1p0}` keyed, so the arms never collide.

## ⚠️ MUST rebuild the binary on Expanse first
`--cfl-dg-safety` is new (2026-07-01). A pre-existing binary **silently ignores**
it (argv-scan) and runs `cfl_dg_safety=3.0` for **both** arms → arm B == arm A
with no error. Rebuild from this branch:
```bash
bash <root>/build_expanse.sh
```
**Detector:** the driver now logs `[time] cfl = 0.5 ... cfl_dg_safety = <v>`. The
script greps for it and WARNS if absent (⇒ stale binary). Arm B's log must show
`cfl_dg_safety = 1` and `dt_cfl` ≈ 3× arm A's `dt_cfl`.

## What to compare
1. **Speedup** — from each log: `[time] nsteps = N` and `[time] dt_cfl = …`.
   Expect `nsteps_A / nsteps_B ≈ 3` and `dt_cfl_B ≈ 3·dt_cfl_A`. Per-macro-step
   cost is dt-independent, so throughput (sim-time reached in the fixed 2 h wall,
   or wall-to-tfinal if completed) scales the same ≈3×.
2. **Fingerprint (C9 — the sole science arbiter)** — diff the two runs' ParaView
   fault output at **matched sim-times** (both arms share tfinal + snapshot
   schedule → snapshots land at identical sim-times over the overlap). Compare
   peak slip-rate, rupture-front arrival, dip-slip. **Bank `cfl_dg_safety=1.0`
   only if the fingerprint is clean;** otherwise the largest clean value is the
   ceiling (sweep 3.0 → 2.0 → 1.5 → 1.0 by re-submitting with those
   `SAFS_CFL_DG_SAFETY` values).

## Wall vs sim (pick one)
- **Throughput mode (default, `SAFS_TFINAL=100s`):** neither arm completes in 2 h;
  each checkpoints. Arm B reaches ≈3× arm A's sim-time — that ratio *is* the
  speedup. Fingerprint on the `[0, t_A]` overlap (covers nucleation + breakout).
- **Completed-wall mode:** submit both with `SAFS_TFINAL=10s` (small enough that
  the slow arm A finishes in 2 h); compare total wall directly. Recommended if
  you want a single clean speedup number.

## Caveats
- **Scheme scope.** Lever A only speeds the **ADER** path. This A/B forces
  `--mixed-flux none --time-integrator ader` (the SeisSol-faithful upwind+ADER
  scheme = the caliper baseline 51443560), so the lever is active. The config
  file's *stated* default is `rk4 + adjacent` (central flux), where
  `cfl_dg_safety` is **inert** (`RkCflFactor` governs). If the production science
  target is the rk4+adjacent variant, Lever A does not apply to it — confirm
  which scheme is the target before banking.
- **ALT nucleation.** The v3 ALT deck's `Tnuc` center was re-snapped on-fault
  (memory `project_safs_alt_nucleation_offfault_recenter`); confirm the log shows
  `V_max` rising (a non-nucleating run has no fingerprint to compare).
- **Approval.** Per project policy, cluster submission waits for explicit user
  approval; this directory only *prepares* the job.

## Next: combined CFL + balanced-MPI A/B
Lever B (the fault-weighted `ncon=2` partition) is **not implemented yet** (the
`--partition-fault-weighted` flag / companion `.cpp` do not exist). The combined
"same-CFL + balanced-MPI" A/B needs Lever B built first; then it becomes a third
arm (B + partition) on top of this CFL A/B — evaluated against the **post-CFL**
wall (Lever A shrinks the partition's absolute imbalance-seconds ~3×).
