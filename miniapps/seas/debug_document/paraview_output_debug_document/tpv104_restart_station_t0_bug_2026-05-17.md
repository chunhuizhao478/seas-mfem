# TPV104 Restart — Station File Starts at t=0 Bug (Frontera job 7729560)

**Date:** 2026-05-17
**Symptom:** Validation #5 of `jobs/tpv104/tpv104_restart_test_v1_dev_2hr.sbatch`
on Frontera FAILed with:

```
Validation #5: Phase B station files start near t = 2.5 s
  First t in .../tpv104_restart_seg002_jobXXX_station_surf_0_9.dat:
    0.0000000000e+00
  FAIL — Phase B station file starts at t=0.0000000000e+00,
         expected >= 2.4 (Phase A tfinal was 2.5).
```

Validations #1–#4 all PASSed; restart-load fired correctly per Validation #4:

```
TPV104 restart loaded: t=2.5 s, step=8945, dt(saved)=0.000279492.
Resuming time loop.
```

So the V1 checkpoint round-trip is fine; the bug is in the **station-writer
initialisation sequence**.

## Root cause

`drivers/tpv104_driver.cpp:1902` and `:1908` unconditionally write a t=0 row
to both station files **before** the restart-load block runs (which is at
line 2349).  Sequence:

```
Line 1894  TPV104StationWriter station_writer;
Line 1896  station_writer.Open(output_dir, output_prefix, ...);     // creates file in segment_002
Line 1902  station_writer.WriteStep(0.0, dof_data);                 // writes t=0 row (BUG)
Line 1906  surface_writer.Open(output_dir, output_prefix, ...);     // creates file in segment_002
Line 1908  surface_writer.WriteStep(0.0, Q);                        // writes t=0 row + zero-Q (BUG)
...
Line 2333  real_t t = 0.0;
Line 2349  if (!restart_prefix.empty()) {                           // restart-load block
              ReadTpv104Checkpoint(..., t, ...);                    // t becomes 2.5
              ...
           }
Line 2387  for (int step = restart_step; step < nsteps; ++step) {
Line 2592    station_writer.WriteStep(t, dof_data);                 // appends t=2.5+, t=2.5+output_dt, ...
Line 2593    surface_writer.WriteStep(t, Q);
              ...
           }
```

So Phase B's station file ends up with:
- Row 1: `t = 0.0` (from line 1902, written before restart load)
- Row 2+: `t = 2.5...` (from in-loop writes after restart)

Validation #5 reads row 1 → reads 0.0 → FAIL.

The surface_writer at line 1908 has the SAME bug and an additional one: at
the time of the t=0 write, `Q` is still the uninitialised zero vector
(restart overwrites it at line 2354 via `ReadTpv104Checkpoint`).  So that
initial entry is wrong on BOTH axes.

## Why Validation #5 caught it on Frontera but not locally

The local unit test (`test_tpv104_checkpoint.cpp`) does not spin up the full
driver — it only round-trips the checkpoint format.  The driver's station
writer is exercised only by the cluster sbatch.  Validation #5 reads row 1
of the .dat file; this is the FIRST FAILURE OBSERVED IN THE FULL DRIVER for
this bug because the BP5 driver's analogous restart path (which DOES skip
the initial probe write on restart per its own pattern) doesn't share code
with the TPV104 driver.

## Fix

Guard both initial-write calls with `restart_prefix.empty()`:

```diff
@@ drivers/tpv104_driver.cpp:1900-1908
    station_writer.Open(output_dir, output_prefix, stations,
                        fault_coords, num_fault_local);
 #endif
-   station_writer.WriteStep(0.0, dof_data);
+   if (restart_prefix.empty())
+   {
+      // Write the pre-evolution t=0 state only on a fresh run.
+      // On restart, the first in-loop WriteStep at step==restart_step
+      // (or the next multiple of output_interval) writes the first
+      // row at t = restart_t (= 2.5 in the sbatch test).  Without
+      // this gate, Phase B's station file starts at t=0 even though
+      // restart loaded t=2.5 (Frontera job 7729560 Validation #5).
+   }

    auto surface_stations = DefaultSurfaceStations_TPV104();
    TPV104SurfaceStationWriter surface_writer;
    surface_writer.Open(output_dir, output_prefix, surface_stations,
                        pmesh, fes);
-   surface_writer.WriteStep(0.0, Q);
+   if (restart_prefix.empty())
+   {
+      // Same gate as station_writer: on restart, Q is still zero at
+      // this point — restart load at line 2354 overwrites it.  The
+      // first in-loop WriteStep writes the LOADED Q at t = restart_t.
+      surface_writer.WriteStep(0.0, Q);
+   }
```

Wait — the second hunk's `WriteStep(0.0, Q);` is INSIDE the `if (restart_prefix.empty())`
guard, so it only runs on fresh runs.  On restart it is skipped — correct.

Actually let me re-write the fix more cleanly:

```diff
-   station_writer.WriteStep(0.0, dof_data);
+   // Write the pre-evolution t=0 state ONLY on a fresh run.  On
+   // restart, the first in-loop WriteStep at step >= restart_step
+   // writes the first row at t = restart_t (= 2.5 in the sbatch
+   // test).  Without this gate, Phase B's station file starts at
+   // t=0 even though restart loaded t=2.5 (Frontera job 7729560
+   // Validation #5).
+   if (restart_prefix.empty())
+   {
+      station_writer.WriteStep(0.0, dof_data);
+   }
```

Same for surface_writer.

## Validation expectation after fix

With the fix:
- Phase A (fresh): row 1 of station file is t=0 (initial), row 2+ at output cadence.  Unchanged.
- Phase B (restart): row 1 of station file is t ≈ 2.5+δ (first in-loop write at
  `step == restart_step` or the next multiple of `output_interval` thereafter).
  Validation #5's threshold `>= 2.4` will PASS.

Local unit test impact: NONE — `test_tpv104_checkpoint.cpp` does not touch
the driver's station writer init path.  Cluster re-submission is required to
confirm.

## Secondary issue (not blocking validation #5)

The Phase B run also emitted HDF5 cleanup noise at end-of-run:

```
HDF5-DIAG: Error detected in HDF5 (1.14.6):
  #000: ../../src/H5P.c line 1468 in H5Pclose(): can't close
  #003: ../../src/H5Iint.c line 948 in H5I__dec_ref(): can't locate ID
        (already closed?)
```

Pattern: "decrementing ID failed / already closed" suggests double-close of
HDF5 property lists / groups / files at ParaView shutdown.  This is
cosmetic — Phase B's fault.vtkhdf was correctly written (237 MB per
Validation #2) and Validation #3 confirmed Phase A unchanged.  Likely
related to `ParaViewHDFDataCollection` close path being called twice (once
explicitly in destructor cleanup, once implicitly when HDF5 plugin
finalises).  **Out of scope for this bug; track separately if needed.**

## Subsequent test improvement (R-109 candidate)

The cluster validation #5 caught this; the unit test missed it.  An
opt-in driver-level smoke test that spins up a small mesh, writes a
checkpoint, restarts from it, then greps the station file for the FIRST
t value would catch this class of bug locally.  Out of scope for the
immediate fix; tracked here as a follow-up.
