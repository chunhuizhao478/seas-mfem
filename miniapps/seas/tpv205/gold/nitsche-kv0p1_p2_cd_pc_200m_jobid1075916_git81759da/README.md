# TPV205 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p2, 200 m (Vista job 1075916)

- Code: seas-mfem branch nitsche-benchmarks, commit 81759da2fb36f180fbf3bdd087318b6d680b732f (81759da, on GitHub),
  built on TACC Vista by
  jobs/vista_nitsche/build_verify_vista.sbatch (clean checkout, BUILD_COMMIT_VISTA = 81759da).
  Driver seas_nitsche_dyn_driver. Same commit, deck and mesh as the 200 m p1 gold
  (nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da); only the order differs.
- Deck and job script: jobs/ (= tpv205/configs/tpv205_nitsche.toml and jobs/vista_nitsche/{tpv205_200m_p2.sbatch,
  run_nitsche_vista.sh} at 81759da; the sbatch and deck are byte-identical to the copies the job saved). Linear
  slip weakening (mu_s 0.677, mu_d 0.525, d_c 0.40 m), three prestress patches at 7.5 km depth; Prakash-Clifton
  v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale "node", kv_length "inscribed"; cfl 0.5.
- Mesh: tpv205/mesh/tpv2053d_200m.msh, made on Vista by build_verify_vista.sbatch (gmsh 4.15.0, msh 2.2):
  1,894,583 tets, 26,228 fault faces, sha256 30c63b923781177830330a36aebcc4159cfc8043374ea4dcd883a9b2eec13538
  (not in the repo). The mesh is not mirror-symmetric about the fault (see the normal-stress note below).
- Run: Vista gg, 16 nodes / 2304 ranks (144 per node), account EAR26024, 2026-10-17 3:05-3:50 am PDT
  (44:54 job wall, about 4.0 SU), 0-12 s, 88,725 steps, dt 1.352508279e-4 s, central difference ("cd"),
  2670 s driver wall (293 s of it fault output), max|V| 10.20 m/s, no NaN.
- results/: the 16 SCEC station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 = strike;
  SI; sigma_n compression positive; tpv205_station_conventions.txt), the station snap list, driver.log and
  RUN_INFO.txt. stations.sha256 holds the sha256 of each station file as written on Vista (checked after the
  transfer). Read with `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (Vista scratch) is not here.
- plots/: tpv205_200m_p2_vista1075916_stations.pdf and its page PNGs (slip-rate overview of the 16 stations,
  then per station slip rate, front zoom, slip, shear stress, normal stress, dip slip rate; this run vs
  DRDG3D 100 m and 200 m, the mixed-flux p2 gold and the 200 m and 100 m Nitsche p1 golds); compare.txt /
  compare_summary.json from compare_nitsche_fault_stations.py vs DRDG3D 100 m, DRDG3D 200 m, both mixed-flux
  golds and both Nitsche p1 golds (command below).
- vs DRDG3D 100 m (median |misfit| of 16 stations): arrival 0.023 s (late at 11 of 16), final slip 0.47 %
  (low at all 16), peak h-slip-rate 1.8 %; noise RMS sigma_n / tau_h / tau_v 0.021 / 0.018 / 0.023 MPa
  (DRDG3D 100 m tau_h / tau_v 0.017 / 0.021). PASSES all three pass tests, the first Nitsche run to do so:
  arrival within 0.1 s everywhere (worst +0.076 s at x2 = 7.5 km, 7.5 km depth); slip within 2 % everywhere
  (worst -1.53 % at x2 = -12 km on the surface); peak within 10 % at 14 of 16 (need 14). The two misses are
  x2 = -12 km on the surface, 28.1 % low (mixed-flux p2 16.4 % low there, DRDG3D 200 m 2.3 % low), and
  x2 = 12 km at 7.5 km depth, 11.4 % low (DRDG3D 200 m 10.7 % low). Every median is far below the p1 golds'
  (200 m: 0.043 s, 1.46 %, 19.3 %; 100 m: 0.022 s, 1.15 %, 8.6 %).
- Fault-wide normal stress (driver.log): 58.8-183 MPa after the rupture reaches the fault edges, against
  120 MPa at rest; the stations see 0.021 MPa RMS. Cause: the tets on the two sides of each fault face differ
  in shape on this gmsh mesh; a mesh mirrored across the fault keeps 120 MPa exactly (1000 m tests,
  2026-10-17). Mirror-mesh reruns use the scripts of commit 41879e8.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-benchmarks
    cd seas-mfem-benchmarks && git checkout 81759da2fb36f180fbf3bdd087318b6d680b732f
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv2053d_200m.msh, gate tests
    sbatch tpv205_200m_p2.sbatch

build_verify_vista.sbatch generates tpv205/mesh/tpv2053d_200m.msh from the committed tpv205/mesh/tpv2053d_200m.geo with Vista's gmsh 4.15.0;
gmsh output depends on the platform, so compare its sha256 with the one above. The job refuses to start unless
the checkout is clean, HEAD is on GitHub and the driver was built at HEAD. Compare (from miniapps/seas, with
D = this folder; rerunning it on these results/ reproduces compare.txt apart from its last line, the path it
wrote the JSON to, and compare_summary.json to round-off):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv205/stations/tpv205_fault_stations.csv \
        --prefix tpv205 --run vista_p2_200m=$D/results \
        --ref drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4 --ref drdg3d_200m=tpv205/benchmark_data/DRDG3D_200m_O4 \
        --gold mixedflux_p1=tpv205/gold/mixedflux_p1_rk45_jobidcombined_7762879_7763019_git04674fc/results \
        --gold mixedflux_p2=tpv205/gold/mixedflux_p2_rk45_jobid7762880_git04674fc/results \
        --gold vista_p1_200m=tpv205/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da/results:nitsche \
        --gold vista_p1_100m=tpv205/gold/nitsche-kv0p1_p1_cd_pc_100m_jobid1076419_git4803a90/results:nitsche \
        --peak-min-stations 14 --json compare_summary.json > compare.txt

- Requested by Chunhui Zhao (TPV205/TPV104 benchmarks, 2026-10-17).
