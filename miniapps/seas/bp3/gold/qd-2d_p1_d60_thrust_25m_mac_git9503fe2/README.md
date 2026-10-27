# BP3-QD gold: dip 60 thrust, quasi-dynamic cycles on a 2-D mesh, p1, 25 m, 0-1500 yr (Mac run at 9503fe25)

- Code: seas-mfem branch 2d-inplane, commit 9503fe25ba84a4a3a999d701760f82b9a0e800d5 (9503fe25, on GitHub),
  built and run on Chunhui Zhao's Mac (macOS arm64, conda env mumps-petsc) in the clone
  /Users/chunhuizhao/projects/seas-mfem-2d-inplane at HEAD 9503fe25 with no tracked file modified (RUN_INFO.txt:
  "dirty 0"; only untracked reference data and binaries). jobs/run.sh is the script as run: it records the commit,
  the tracked-tree state and the binary sha256, then runs the driver through jobs/locked_exec.sh (the Mac's
  one-run-at-a-time lock). The driver was built 15 s before the run by
  `make -f cycle/cycle.mk -j6 SEAS_FS_LIB= seas_fd_cycle_driver` (jobs/make_cycle.log, rc 0). Binary sha256
  4f18eb090554757b9a460d6776de85178f88ff805f13f7ef4f1e16f7a5f4adaf. This is a workstation run, not a cluster job:
  the CLAUDE.md guard script was not used, and the build log does not record `git status` separately.
- Deck: config/fd_cycle/bp3_qd_d60_thrust.toml at 9503fe25, sha256
  5857f2f9e48b0832a3571c84fa200d2b3f6670db852e83f9ebadc554b00da693. SCEC BP3-QD, dip 60 deg, thrust: plane strain,
  free surface, rate and state with the aging law, a from 0.010 to 0.025 between 15 and 18 km down-dip, creep at
  V_p below 40 km (the deck header has the frame and sign conventions). The QD regime only, RK 5(4) Dormand-Prince
  with atol 1e-7 and the infinity norm (tests/verification/petsc_ts_rk45_tandem.cfg, the options of Tandem's
  rk45.cfg). The QD stages go through the 2-D discrete Green's function (phase F); its 466 checks against the solve
  path over 4,662,572 stages differ by at most 1.1e-9 (shear) and 4.6e-10 (normal), relative.
- Mesh: bp3/mesh/bp3_d60_25m.msh (in the repo), sha256 fd08af92b3b48b78770ecfacee64c7b0b7ae93f3607fc079eb714bccff2c6d3c;
  106,448 vertices, 212,659 triangles, 1,600 fault segments of 25 m. Tandem p = 1 ran on the same file.
- Run: `mpirun -np 4`, OMP_NUM_THREADS=1, 2026-10-27 04:30:42 to 09:52:38 PDT, wall 19,316 s (5 h 22 min), peak
  resident memory 8.2 GB; 0 to 1500 yr (4.7304e10 s, 365-day years), 743,001 accepted steps and 33,966 rejected,
  16 events, exit 0 ("status: done").
- results/: what the driver wrote, all stamped with the commit in their headers: the 12 fault and 8 surface
  station files (28,004 rows each), the five dip profiles (1,993 rows of 81 points), the global file (one row per
  step), the source time series, the event catalog and event list, and the rupture times of each event;
  outputs.sha256 holds their sha256. RUN_INFO.txt, run_stats.log, run.out and t0.txt are the run script's records.
  results/scec_export/: the 23 SCEC upload files made from them by cycle/scec_bp3_export.py (command below;
  scec_export.sha256, scec_export_stdout.txt). Not kept: the 4 checkpoint files and the 22 MB driver log (on the
  Mac in seas-mfem-2d-inplane-results/g3-bp3/seas_p1_d60_thrust/; every file here was moved off the Mac gzipped and
  checked against the sha256 of the Mac's copy).
- plots/: cycle_compare_ours_vs_tandem_p1_1500yr.pdf (24 pages: event table, onsets against event number, one page
  per fault and surface station, notes), event_table.txt, largest_differences.txt, front_offsets.txt and
  page_order.txt, made by scripts/plot_bp3_qd_tandem_compare.py (command below) from results/scec_export and
  bp3/benchmark_data/TANDEM2D_p1_25m. Tandem p = 1 is lighter underneath and ours red on top; each station page
  also shows ours minus Tandem between events and in the event-16 window.

## Against Tandem p = 1 (bp3/benchmark_data/TANDEM2D_p1_25m: same mesh, same order, run on the same Mac)

Both runs are cut at 1499.9083 yr, the end of Tandem's output. Onsets of both codes are taken from the station data
with one definition (the first time max V over the 12 fault stations exceeds 1e-3 m/s after more than 1 yr below).

- Events: 16 in each code. Onsets agree within 1.04e-4 yr (largest at event 16, about 55 minutes); the mean
  recurrence is 85.288 yr from event 2 and 87.628 yr from event 3 in both codes. The SCEC submissions (Tandem p8
  250 m, fdra and hbi 25 m) have about 85.0 yr and the first event near 177.2 yr: at p = 1 both codes run about 3 %
  long, and their first event comes 0.5 yr early.
- Slip one year after each onset: the same within 0.5 mm at every fault station and event. Event durations within
  0.4 s. In event 16 our front reaches each fault station 0.10 to 0.16 s earlier (front_offsets.txt; dp050 is
  explained below).
- Between events (onset - 1 day to onset + 30 days of every event left out): shear stress within 0.016 MPa, log10 V
  within 0.009 and slip within 2e-4 m at every fault station, except dp150 before 70 yr (below).
- **Differences in peak magnitude come from the output sampling, adaptive (Tandem) against fixed (ours), not from
  the solution** (Chunhui Zhao, 2026-10-27: "indicate that the difference in magnitude is due to adaptive/fix
  sampling"). The two codes write rows on different schedules:
  - Tandem samples adaptively (app/tandem/AdaptiveOutputStrategy.h): a row when the fault-wide max slip rate has
    changed by 1 % since its last row (rtol 0.01, atol 1e-50), or after t_max = 0.3 yr (bp3.toml). Where max V
    changes fast it writes every few ms; on the mid-rupture plateau (about 17 to 22 s after onset) it writes nothing
    for up to about 2 s.
  - seas-mfem samples on a fixed clock ([cycle.output]): a row every 0.1 s while max V > 1e-3 m/s, else every
    0.1 yr.
  - Consequences: at dp000 and the surface stations Tandem's ms rows catch spikes of a few hundredths of a second
    that our 0.1 s rows miss, so Tandem's peaks are larger (dp000 peak V, mean over events 2 to 16: Tandem 5.64 m/s,
    ours 3.67 m/s, Tandem read at our row times 3.47 m/s; peak shear stress 141.7, 115.8 and 121.6 MPa). Sampled the
    same way the peaks agree within about 5 %. At dp050, from event 3 on, the front passes inside Tandem's
    mid-rupture gap, so Tandem's file holds neither the arrival nor the peak (its front reads 0.3 to 1.4 s late and
    its peak V is 0.27 to 0.47 m/s against our 0.58 to 0.84 m/s); in events 1 and 2, where the front falls on Tandem
    rows, the codes agree (peak V 0.227 and 0.226 m/s, 1.01 and 0.98 m/s; front within 0.16 s). Use the onset and
    slip columns of the event table, not peak V.
- Open: dp150 (15 km down-dip, the kink of a) before event 1, between about 30 and 70 yr: our shear stress is up to
  0.376 MPa lower and log10 V up to 0.325 lower (at 63.5 yr); after 100 yr it agrees within 0.0055 MPa. Both codes
  write densely there, so this is not sampling. The station sits on a mesh vertex where two fault faces meet and
  each keeps its own value; ours reads the up-dip face, and Tandem's values lie on the deeper, higher-stress side.
  That the codes read different faces is the likely cause but is not confirmed (it needs stations a few metres to
  either side). dp125 and dp175 agree within 4e-4 MPa over the same years.
- fn-00 and fn+00: ours evaluates both on the trace (x1 = 0) as the limit from each side; Tandem's sit at x1 = -1 m
  and +1 m. Their values differ by up to 3e-4 m between events for that reason.
- Tandem's surface velocity is du/dt by finite differences (its probes write displacement only).
- Our post-event output: 10 s after max V falls below 1e-3 m/s the rows fall back to every 0.1 yr, so afterslip
  peaks 15 to 30 minutes after an event (dp025, dp175, dp200) are not sampled. Up to the end of the 0.1 s rows our
  slip rate matches Tandem's.

- Status: BP3-QD dip 60 thrust at p = 1 against Tandem p = 1 is completed for the 2-D in-plane design (gate G3;
  Chunhui Zhao, 2026-10-27: "gold the results p=1 comparison with tandem p = 1 and move on"); see
  document/2d_inplane/BENCHMARK_STATUS.md.

## Reproduce

In a checkout of 9503fe25 with MFEM built (MPI, MUMPS, METIS) and PETSc, from miniapps/seas:

    make -f cycle/cycle.mk seas_fd_cycle_driver SEAS_FS_LIB=
    mpirun -np 4 ./seas_fd_cycle_driver --deck config/fd_cycle/bp3_qd_d60_thrust.toml --output-dir <out>
    python3 cycle/scec_bp3_export.py <out> bp3_qd_d60_thrust <scec> --modeler "Chunhui Zhao" \
        --element-size "25 m" --order 1 --scenario "dip 60, thrust"

On another machine or rank count the station data are expected to agree to round-off, not byte for byte (MUMPS;
not tested for this deck). Plots (from
miniapps/seas, D = this folder):

    python3 -I scripts/plot_bp3_qd_tandem_compare.py $D/results/scec_export \
        $D/results/bp3_qd_d60_thrust_events.txt bp3/benchmark_data/TANDEM2D_p1_25m $D/plots

- Requested by Chunhui Zhao (2-D in-plane design, gate G3, 2026-10-27).
