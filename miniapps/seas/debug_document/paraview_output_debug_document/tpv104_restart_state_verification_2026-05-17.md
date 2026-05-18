# TPV104 Restart — Complete State Verification (Phase C + Seam + Reference)

**Date:** 2026-05-17
**Driver:** user request "let's add all to complete the restart check"
**Scope:** extend `jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch`
from 7 file-existence / file-isolation validations to 11 validations
including state-trajectory comparison against a no-restart reference run.

## Gap closed

Prior sbatch validations #1–#7 verified:
- ParaView + station outputs are produced (#1, #2, #6).
- Phase A's outputs are not clobbered by Phase B (#3 fault.vtkhdf size,
  #7 station md5).
- Restart code path fired and time scalar was loaded (#4 log line, #5
  first-t threshold).

What they did NOT verify:
1. **State correctness.**  A hypothetical field-order swap in
   `WriteTpv104Checkpoint` (e.g., psi ↔ slip_rate columns) would still
   produce files of the right size in the right directories, log
   "TPV104 restart loaded:" correctly, and pass first-t ≥ threshold.
   The simulation would continue with corrupted state but all 7
   validations would PASS.
2. **Trajectory continuity vs no-restart reference.**  No single-shot
   run existed to compare against.
3. **Final-time correctness.**  Phase B's last station entry was not
   checked against tfinal.

This update adds:
- A third **reference Phase C** that runs 0 → tfinal_C in one shot
  (no restart).  Same parameters as A+B, same output directory layout
  under `segment_ref/`.
- **Validation #8**: Phase B's last station t ≥ tfinal − ε.
- **Validation #9**: Phase C's fault.vtkhdf landed > 1 KiB.
- **Validation #10 SEAM**: at the central nucleation station
  `x2_0_x3_7.5`, Phase A's LAST row's `h-slip` and `h-slip-rate`
  agree with Phase B's FIRST row within tolerance.  Catches gross
  state-load failures (e.g., field-order swap that wipes slip).
- **Validation #11 REFERENCE**: at the same station, Phase A+B at
  t=tfinal_B agrees with Phase C at t=tfinal_C within tighter
  tolerance.  Compares `h-slip`, `h-slip-rate`, and `psi`.  This is
  the gold-standard restart-correctness check.

## Time budget

| Phase | Original | New | Wall (@ 2.5s/hr) |
|-------|----------|-----|------------------|
| A     | 0 → 2.5  | 0 → 1.0   | ~24 min |
| B     | 2.5 → 5.0 (restart) | 1.0 → 2.0 (restart) | ~24 min |
| C (NEW) | —       | 0 → 2.0   | ~48 min |
| Total | ~120 min | ~96 min | 24 min margin under 2h dev queue |

A's tfinal=1.0 and B's tfinal=2.0 (absolute) chosen so A+B endpoint
matches C's endpoint exactly at t=2.0 for the reference comparison.

## Station-file column reference

From `dynamic/tpv104_setup.hpp:355` (TPV104StationWriter::WriteStep):

```
# Columns (SCEC TPV104 trace layout):
# t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress  psi
#  $1    $2        $3            $4         $5      $6           $7              $8       $9
```

Restart-correctness validation uses:
- `h-slip` ($2) — integrated quantity, smooth across seam.
- `h-slip-rate` ($3) — evolves rapidly during rupture; sensitive to
  state-load errors.
- `psi` ($9) — directly checkpointed; perfect round-trip canary
  (any psi corruption immediately shows here).

## Tolerance design

**Seam (#10)** — compares t=tfinal_A row of Phase A vs first post-restart
row of Phase B at the same station.  Times differ by O(output_dt) — a
few dt steps.  During rupture, slip-rate can change rapidly across
these few steps.  Loose tolerance just catches gross failures:
- `h-slip`:        abs_tol = 0.1 m  (rupture moves at ~1 m/s × output_dt = 0.05 m)
- `h-slip-rate`:   abs_tol = 5.0 m/s OR rel_tol = 50% (loose; just sanity)

**Reference (#11)** — compares A+B at t=tfinal_B vs C at t=tfinal_C
(same simulation time, different code path).  Checkpoint format is
text 17-digit scientific (round-trips IEEE 754 doubles exactly per
`tpv104_checkpoint.hpp:70`).  Divergence sources after restart:
- MPI reduction-order non-determinism: ~1e-12 per step.
- ADER substep coefficient summation order: ~1e-14 per step.
- Over ~3500 Phase B steps: cumulative ~1e-9 to 1e-7.

Initial tight tolerance:
- `h-slip`:        rel_tol = 1%, abs_tol = 1e-3 m
- `h-slip-rate`:   rel_tol = 5%, abs_tol = 1e-3 m/s
- `psi`:           rel_tol = 0.1%, abs_tol = 1e-6 (psi is dimensionless O(0.5))

If real Frontera output shows wider drift, document and loosen here.

## Implementation

- `jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch` — add Phase C
  invocation; shrink A/B tfinal; add 4 new validations (#8 final-t,
  #9 reference-fault.vtkhdf, #10 seam, #11 reference).  Validation
  helpers are inline awk/bash; no external deps.
- `tests/unit/test_tpv104_checkpoint.cpp` Sub-test 6 — add grep
  assertions that the sbatch contains the new sections so a future
  regression that drops them is caught locally.

## Future tightening

If Validation #11 passes with 1% but the user wants bit-identity
verification:
- Replace text-format checkpoint with binary (eliminates 17-digit
  scientific round-trip noise even though it shouldn't matter).
- Force `MPI_Allreduce` reduction order via custom op (eliminates
  reduction-order non-determinism).  Costs perf.
- Add a Validation #12 that fault.vtkhdf snapshots at t=2.0 from B
  and C are byte-identical (h5diff -c).  Strongest possible check.

These are deferred until the looser tolerances expose real drift on
Frontera.
