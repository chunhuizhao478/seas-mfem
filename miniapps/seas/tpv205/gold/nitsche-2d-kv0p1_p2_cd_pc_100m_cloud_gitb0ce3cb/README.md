# TPV205-2D gold: Nitsche dynamic rupture on a 2-D mesh, Kelvin-Voigt damping, p2, 100 m (cloud run at b0ce3cb6)

- Code: seas-mfem branch 2d-inplane, commit b0ce3cb6d0b0a437ac773d2f61f27241b1771838 (b0ce3cb6, on GitHub),
  built and run in a Claude Code cloud container by jobs/run_gold_cloud.sh: the CLAUDE.md guard (git fetch,
  clean tracked tree, HEAD on GitHub), every `*.o` deleted, seas_nitsche_dyn_driver built at HEAD with
  jobs/mk.sh (the container's make wrapper around cycle/cycle.mk; MFEM built from the same checkout, whose
  files outside miniapps/seas are those of the branch point 3926e9c8), the tree checked unchanged after
  the build. Binary sha256 e7fd5b6f7cdf248bf9064d3d895b88858d31dff7e423b008d4b0f20564d171bf.
- Deck: jobs/tpv205_2d_nitsche.toml (= tpv205/configs/tpv205_2d_nitsche.toml at b0ce3cb6). Plane strain in the
  (x, y) plane, fault on y = 0 from -15 to +15 km; linear slip weakening (mu_s 0.677, mu_d 0.525, d_c 0.40 m),
  sigma_n 120 MPa, tau 70 MPa with the three patches of the 7.5 km depth slice (81.6 MPa on |x| <= 1.5 km,
  78.0 MPa on x in [-9, -6] km, 62.0 MPa on x in [6, 9] km); Prakash-Clifton v* 1 m/s, L 0.01 m;
  [numerics] kv_eta 0.1, kv_scale "node", kv_length "inscribed"; cfl 0.5 of the 2-D C(2) = 0.365;
  absorbing outer boundary; `up = [0, 0, 1]` puts the in-plane slip in the strike slot (slip2 = SCEC h-slip).
- Mesh: tpv205/mesh/tpv205_2d_100m.msh (in the repo, committed in phase B2), 40,526 triangles, 300 fault
  edges, sha256 b1b9467a67d7c53ea22cbd58fd552d0910cd8535f11a917892223fdf5e3d6909.
- Run: cloud container (Intel Xeon @ 2.10 GHz, 4 cores), `mpirun -np 4 --bind-to none`, OMP_NUM_THREADS=1,
  2026-10-26 16:33-16:42 UTC; `--order 2 --log-every 2000`; 0-12 s, 39,136 steps, dt 3.066261374e-4 s,
  central difference ("cd"), driver wall 532.7 s (13.1 ms/step on the slowest rank), max|V| 7.97 m/s, no NaN;
  sigma_n_eff stays in [111.6, 124.7] MPa.
- results/: the 7 on-fault station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip,
  2 = strike; SI; sigma_n compression positive; tpv205_2d_station_conventions.txt), the station snap list
  (every station within 2.2 mm of a fault node), driver.log and RUN_INFO.txt; stations.sha256 holds the
  sha256 of each station file. The `*_copies.dat` files (the two DG copies of each station node) and the
  initial-state CSVs the run also wrote are not kept. The 2-D driver writes no fault ParaView output
  (phase F).
- Reproducibility: the station data (all lines but the commit-stamped headers) are byte for byte those of
  the same deck run at 3bd1bc23 (review 6) in the same container, which matched the Mac's gate G1 run
  (mumps-petsc, 4 ranks, 39,136 steps, same done line). Reviews 7 and 8 leave this deck's path unchanged.
- plots/: tpv205_2d_stations.pdf and its page PNGs (misfit overview, then per station slip rate, front zoom,
  slip, shear stress, normal stress, shear / normal; PyLith 200 m black, this run red), compare_vs_ref.txt
  and summary.json, made by scripts/plot_tpv205_2d_stations.py (command below).
- vs PyLith 200 m (tpv205/benchmark_data/PYLITH2D_200m, the only TPV205-2D code available), window 0-12 s:
  - final slip within 1.4 % from -12 to +4.5 km (+1.40, +0.70, +0.18, +0.86, +0.75 %), +3.06 % at +7.5 km
    and +3.29 % at +12 km;
  - arrival (first |V| > 1 mm/s): PyLith saved one sample every 0.06 s, which delays its arrival by up to
    0.06 s; read at PyLith's sample times, ours is within one sample at all 7 stations (-0.06 to +0.06 s);
    on our own time grid the differences are -0.10 to +0.04 s;
  - peak slip rate: -4.1 % at 0, +5.9 / +7.0 % at -12 / -7.5 km, about +29 % at +-4.5 and +7.5 km and
    +70 % at +12 km. The peaks are not converged at 100 m p2: in the resolution runs of review 6 and G1
    they keep rising with refinement (+12 km at p1: 3.6, 4.3, 5.0 m/s at 200, 100, 50 m), and PyLith's
    0.06 s sampling and 200 m elements lower its peaks too.
  - Open: at +7.5 km (the low-stress patch) slip stays about +3 % and the peak at least +22 % above PyLith
    in every G1 run (p1 at 200, 100 and 50 m and p2 at 100 m; +29.9 % here). PyLith's header places this station at 7.6 km (its
    200 m node), which was not tested as the cause; a second TPV205-2D code would settle it.
- Status: TPV205-2D is completed for the 2-D in-plane design (gate G1; Chunhui Zhao 2026-10-26: "gold
  tpv205-2d, tpv105-2d results p=2, and mark as completed"); see document/2d_inplane/BENCHMARK_STATUS.md.

## Reproduce

In a checkout of b0ce3cb6 with MFEM built (MPI, MUMPS, METIS) and the build variables of jobs/mk.sh adapted
to the machine, from miniapps/seas:

    make -f cycle/cycle.mk seas_nitsche_dyn_driver        # plus the MFEM_BUILD_DIR etc. of jobs/mk.sh
    mpirun -np 4 ./seas_nitsche_dyn_driver --config tpv205/configs/tpv205_2d_nitsche.toml \
        --order 2 --log-every 2000 --output-dir <out>

jobs/run_gold_cloud.sh is the script as run (container paths). Plots (from miniapps/seas, D = this folder):

    python3 -I scripts/plot_tpv205_2d_stations.py $PWD PyLith2D_200m=tpv205/benchmark_data/PYLITH2D_200m \
        $D/plots 12 gold_p2_100m=$D/results

- Requested by Chunhui Zhao (2-D in-plane design, phase C gate G1, 2026-10-26).
