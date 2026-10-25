# Horizon GPU job 69011 against the CPU gold 1074099 and Tandem p=1

These files compare the GPU full-cycle run in the folder above (`../bp5_cycle_p1_*.txt`) with the
Vista CPU gold (`../../golden_vista_144r_p1_1000m/`, job 1074099) and Tandem's p=1, 1000 m reference
(`../../golden_tandem_p1_1000m/`), in three pairs: GPU - CPU, GPU - Tandem and CPU - Tandem. They
were made on 2026-10-25 after the run ended and are copied here unchanged.

## Files

- `cmp_like_tandem.py`: makes the tables. `<pair>_events.csv` has one row per event (1-7) and
  station: onset difference, peak log10 V_strike of each code and their difference, and the strike
  slip one year after the onset. `<pair>_interseismic.csv` has one row per window between events and
  station: the largest |difference| in tau_strike, slip_strike and log10 V_strike, and when it occurs.
  `summary.json` holds all of it with the event onsets.
- `plot_like_tandem.py`: makes `bp5_gpu69011_vs_cpu1074099_vs_tandem_p1.pdf` from the tables and the
  station files. Tandem is drawn first, the CPU next and the GPU last, on top.

From `miniapps/seas/bp5/benchmark_data`:

```
H=golden_horizon_gb200_p1_1000m; C=$H/compare_tandem_p1
python3 $C/cmp_like_tandem.py $H golden_vista_144r_p1_1000m golden_tandem_p1_1000m OUT
python3 $C/plot_like_tandem.py $H golden_vista_144r_p1_1000m golden_tandem_p1_1000m OUT OUT/cmp.pdf OUT/png
```

The first command gives the six CSV files and `summary.json` here byte for byte (checked 2026-10-25).

## Method

The method of `../../golden_vista_144r_p1_1000m/compare_tandem_p1/`; on the CPU - Tandem pair it gives
the same values as the tables there (the interseismic table is identical byte for byte).

- Times are in years of 365.25 days.
- A station's onset is its first upward crossing of log10 V_strike = -3, interpolated linearly. An
  event's onset is the earliest station onset.
- The peak is searched from 1 h before to 0.5 yr after the event onset.
- Between events, from 1 yr after an onset to 1 day before the next, the differences interpolate one
  code onto the other's sample times, in both directions, and keep each point only where the
  interpolated code is sampled at least as densely there.

## Results

- All three have 7 events in 1500 yr (Tandem's run ends at 1450.6 yr). The GPU's event onsets are at
  0, 225.289, 450.991, 676.705, 902.432, 1128.161 and 1353.891 yr.
- Onset differences (s), events 1-7. GPU - CPU: 0, +1.09, +3.02, -1.79, -1.80, -0.59, -3.60.
  GPU - Tandem: 0, +1.61, +3.11, +0.91, +2.70, +5.57, +4.58. CPU - Tandem: 0, +0.52, +0.09, +2.70,
  +4.50, +6.16, +8.18.
- At +00dp+00, +00dp+10, +00dp+22, -16dp+10, +16dp+10 and -24dp+10, the GPU differs from Tandem by at
  most 9.6e-4 in peak log10 V and 6.8e-5 m in slip one year after onset. Between events, by at most
  7.0e-3 MPa in tau, 1.6e-5 m in slip and 1.0e-2 in log10 V. The CPU's numbers against Tandem are
  9.6e-4, 8.5e-5 m, 7.5e-3 MPa, 2.1e-5 m and 1.2e-2.
- The four surface stations off the center (-36dp+00, -16dp+00, +16dp+00, +36dp+00) have the same
  offset from Tandem in both runs: tau up to 0.128 MPa between events and peak log10 V up to 2.0e-2.
  This is the open item of the CPU gold; the GPU run does not change it.
- GPU - CPU: at the 10 stations, the peak log10 V differs by at most 1.3e-4 and the slip one year after
  onset by at most 3.7e-5 m.
- GPU - CPU between events: the largest differences (0.17 MPa, 4.0e-2 m, 0.19 in log10 V, at
  +00dp+10 and +16dp+10, 12.5 days before event 2 and 2.6 days before event 7) come from the
  interpolation. There, in the last days of nucleation, the runs write station rows 1.6 to 36 days
  apart while the slip rate grows about tenfold, and the interpolated run is linear across that gap.
  The rows themselves agree: at +16dp+10 the GPU's row at 225.2547 yr has 7.49126 m of slip, between
  the CPU's rows at 225.2439 yr (7.46514 m) and 225.2859 yr (7.70765 m).
