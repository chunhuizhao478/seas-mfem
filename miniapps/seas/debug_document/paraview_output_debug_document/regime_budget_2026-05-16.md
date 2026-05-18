# ParaView regime budget — empirical measurement from BP5 production run

**Date:** 2026-05-16
**Source run:** `~/Downloads/seas-mfem/results_phase7_production_job7650438/`
(SLURM job 7650438, BP5-QD, 250 yr target, completed to ~1800 yr)
**Station analyzed:** `bp5_phase7_prod_job7650438_fltst_strk+00dp+10.txt`
(in the BP5 nucleation patch; x2 = 0 km along strike, x3 = 10 km depth)

## What this document is

Concrete numbers for how much simulated time BP5-QD spends in each ParaView
adaptive-schedule regime over a multi-event run. Use this to size
`--paraview-max-snapshots`, set `--paraview-dt-*` overrides, or budget disk
quota before submitting a new production sbatch.

## Regime classification (current defaults)

Decided by global `V_max` (m/s) at each ODE step:

| Regime         | V_max range            | `dt_*` default               | Cadence |
|----------------|------------------------|------------------------------|---------|
| Interseismic   | V ≤ 1e-6               | `1.0 * seconds_per_year`     | 1 yr    |
| Nucleation     | 1e-6 < V ≤ 1e-3        | `1.0`                        | 1 s     |
| Coseismic      | V > 1e-3               | `0.01`                       | 10 ms   |

Defined in `miniapps/seas/io/paraview_output.hpp:162-170`. NaN/Inf V is
classified as coseismic by `NextRegime` (line 219) so output stays dense
during a blowup.

## Measured elapsed time per regime (this run)

Total simulated time: **1775.58 yr** (5.603 × 10¹⁰ s)

| Regime         | Elapsed time            | % of run     | Rows written |
|----------------|-------------------------|--------------|--------------|
| Interseismic   | 1775.22 yr (≈ 5.60e10 s) | 99.98%       | 283,290      |
| Nucleation     | 134 days (≈ 1.16e7 s)   | 0.0206%      | 38,990       |
| Coseismic      | 51 min (3,080 s)        | 0.0000055%   | 45,893       |

## Per-event breakdown (≈ 7 events in 1775 yr at ~240 yr recurrence)

| Phase          | Per-event duration      | Per-event writes @ default cadence |
|----------------|-------------------------|------------------------------------|
| Nucleation     | ~19 days                | ~1.6 million (at dt = 1 s)         |
| Coseismic      | ~7 min                  | ~42,000 (at dt = 0.01 s)           |

## What this means for `--paraview-max-snapshots`

The Phase 6 sbatch
(`miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch`)
sets `--paraview-max-snapshots 5000`. The Python estimator
(`scripts/_io_size_schedule.py`) projects ~86k uncapped writes
*per event*, dominated by nucleation, assuming
`NUCLEATION_DURATION_S["bp5"] = 86,400 s` (1 day).

Empirically the nucleation phase lasts **~19 days per event**, not 1 day.
Real per-event nucleation writes at 1 s cadence are therefore **~1.6 M**,
not 86k — a **19× underestimate** in the planner.

Combined with R-002 (current `SnapshotCapAwareInterval` only caps the
interseismic regime), the 5000 cap is unreachable on this workload by
two orders of magnitude. The user observed 11,408 writes at simulated
time ≈ 114 s in the failing Phase 6 sbatch — consistent with the system
entering coseismic immediately and writing every 10 ms for ~114 s.

## Recommendations

For a 250 yr × ~1 event production run, the conservative numbers based
on the measured per-event durations above:

| Goal                                     | Suggested CLI flags                                              |
|------------------------------------------|------------------------------------------------------------------|
| Capture full resolution                  | (no overrides) → expect ~1.6 M nucleation writes per event       |
| Drop nucleation cadence 100×             | `--paraview-dt-nu 100`  → ~16k nucleation writes per event       |
| Drop coseismic cadence 10×               | `--paraview-dt-co 0.1`  → ~4,200 coseismic writes per event      |
| Fit under 5000 total snapshots / event   | `--paraview-dt-nu 1000 --paraview-dt-co 1.0`                     |
| Stretch interseismic 10×                 | `--paraview-dt-inter-yr 10`                                      |

These do NOT replace the R-002 fix — until the cap is hard, any of the
above can be silently overshot if a regime stays active longer than
projected. They reduce the *projected* count, not the *bounded* count.

## How to reproduce

```python
# Run from any directory with the station files in it.
import glob
f = "bp5_phase7_prod_job7650438_fltst_strk+00dp+10.txt"
t_inter = t_nucl = t_coseis = 0.0
prev_t = None
prev_v = None
with open(f) as fh:
    for line in fh:
        if line.startswith("#"): continue
        parts = line.split()
        if len(parts) < 5: continue
        t = float(parts[0])
        # Columns: time, slip_s, slip_d, log10(V_s), log10(V_d), tau_s, tau_d, log10(state)
        lvs, lvd = float(parts[3]), float(parts[4])
        v = max(10**lvs if lvs > -290 else 0.0,
                10**lvd if lvd > -290 else 0.0)
        if prev_t is not None:
            dt = t - prev_t
            if   prev_v > 1e-3: t_coseis += dt
            elif prev_v > 1e-6: t_nucl   += dt
            else:               t_inter  += dt
        prev_t, prev_v = t, v
print(f"interseismic = {t_inter:.3e} s ({t_inter/3.15576e7:.2f} yr)")
print(f"nucleation   = {t_nucl:.3e} s ({t_nucl/86400:.2f} days)")
print(f"coseismic    = {t_coseis:.3e} s ({t_coseis/60:.2f} min)")
```

## Caveats

1. **Single-station approximation.** This uses local V at one station
   in the nucleation patch. Global V_max (which the scheduler actually
   uses) is the max across the whole fault and can be slightly larger,
   so the real nucleation / coseismic durations are likely a bit
   longer than reported here. Order of magnitude is reliable.
2. **Files have unequal row counts** across stations
   (366,819 to 372,912). All stations share the same time stamps for
   the rows that exist on both; this run was cut off mid-write so
   later-flushed station files reach further in time. Cross-station
   max-V analysis was attempted but failed the equal-length assertion
   and was abandoned in favour of single-station numbers.
3. **The bench_out probe file uses the same adaptive schedule as
   ParaView** (`BP5BenchmarkOutput::OutputInterval` constructs a
   default `AdaptiveSchedule` and queries `Interval`). So this regime
   count is also the bench_out write count, NOT just a sampled subset.
4. **Recurrence** of ~240 yr per CLAUDE.md "What Constitutes a
   Regression" — 1775 yr / 240 yr ≈ 7 events used for per-event
   averaging above. If the actual event count differs (e.g., system
   accelerated during the run), per-event numbers scale accordingly.

## See also

- `REVIEW.md` R-002 — `SnapshotCapAwareInterval` cap only applies in
  interseismic regime; coseismic / nucleation are uncapped.
- `miniapps/seas/io/paraview_output.hpp:1807-1829` — the
  `SnapshotCapAwareInterval` implementation.
- `miniapps/seas/scripts/_io_size_schedule.py` — Python estimator;
  `NUCLEATION_DURATION_S["bp5"]` value is the lever to correct to the
  measured ~19 days per event.
- `miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch`
  — the failing sbatch whose snapshot-cap assumptions this document
  invalidates.
