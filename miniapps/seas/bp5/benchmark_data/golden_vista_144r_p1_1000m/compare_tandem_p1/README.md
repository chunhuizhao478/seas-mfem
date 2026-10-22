# Vista job 1074099 against Tandem p=1

These tables compare the full-cycle run in the folder above (`../bp5_cycle_p1_*.txt`)
with Tandem's p=1, 1000 m reference (`../../golden_tandem_p1_1000m/`). They were made
after the run ended and are copied here unchanged.

## Files

- `cycle_events_vs_tandem.csv`: one row per event (1-7) and station. For each code it
  gives the onset time, the peak log10 V_strike and its time, and the strike slip one
  year after the onset, then the run minus Tandem difference of each.
- `cycle_interseismic_vs_tandem.csv`: one row per interseismic window (event k to
  k+1) and station: the largest |difference| in tau_strike, slip_strike and
  log10 V_strike over the window, and when it occurs.
- `cycle_vs_tandem_summary.json`: both tables, the event onsets and the method notes.

## Method

- Times are in years of 365.25 days.
- A station's onset is its first upward crossing of log10 V_strike = -3, interpolated
  linearly. An event's onset is the earliest station onset.
- The peak is searched from 1 h before to 0.5 yr after the event onset.
- Interseismic differences interpolate one code onto the other's sample times, in
  both directions, and keep each point only where the interpolated code is sampled
  at least as densely there.

## Results

- Both codes have 7 events in 1500 yr. The run's event onsets are at 0, 225.289,
  450.991, 676.705, 902.432, 1128.161 and 1353.891 yr. Run minus Tandem: 0, +0.52,
  +0.09, +2.70, +4.50, +6.16 and +8.18 s.
- At +00dp+00, +00dp+10, +00dp+22, -16dp+10, +16dp+10 and -24dp+10, the peak
  log10 V differs by at most 1.0e-3 and the slip one year after onset by at most
  8.5e-5 m. Between events, tau differs by at most 7.6e-3 MPa and slip by at most
  2.1e-5 m.
- Open: the four surface stations off the center (-36dp+00, -16dp+00, +16dp+00,
  +36dp+00) keep a steady tau offset between events of 0.026-0.030, 0.044-0.047,
  0.027-0.032 and 0.126-0.128 MPa. At -36dp+00 and +36dp+00 the peak log10 V is
  0.011-0.020 higher than Tandem's and comes 710-2,480 s earlier. Every event shows
  the same offsets. The cause has not been found.
