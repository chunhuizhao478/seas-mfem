# TPV205 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p1, 100 m (Vista job 1076419)

- Code: seas-mfem branch nitsche-benchmarks, commit 4803a9000b6dfca1be64a42372517c30e6e7d36d (4803a90, on
  GitHub), built on TACC Vista by jobs/vista_nitsche/build_verify_vista.sbatch (clean checkout,
  BUILD_COMMIT_VISTA = 4803a90; build_vista.sh with gcc/14.2.0, openmpi/5.0.5, hypre/2.31.0, metis/5.1.0.3,
  phdf5/1.14.6). Driver seas_nitsche_dyn_driver.
- Deck and job script: jobs/ (= tpv205/configs/tpv205_nitsche.toml and jobs/vista_nitsche/{tpv205_100m_p1.sbatch,
  run_nitsche_vista.sh} at 4803a90; the deck is the same file as in the 200 m gold at 81759da). Linear slip
  weakening (mu_s 0.677, mu_d 0.525, d_c 0.40 m), three prestress patches at 7.5 km depth; Prakash-Clifton
  v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale "node", kv_length "inscribed"; cfl 0.5.
- Mesh: tpv205/mesh/tpv2053d_100m.msh, gmsh 4.15.0 output (msh 2.2) of tpv205/mesh/tpv2053d_100m.geo made on the
  Mac (the mesh of condo job 6490984), copied into the Vista checkout: 5,438,053 tets, 98,984 fault faces,
  sha256 2f5662fb7267824a3e1dbace6564d62508c16aaed506962ad4b0affb63640495 (not in the repo; listed in
  jobs/vista_nitsche/meshes.sha256 at 4803a90, which the run checks). gmsh output depends on the platform, so
  a regenerated mesh may differ in bytes; the sha256 identifies the one used here.
- Run: Vista gg, 8 nodes / 512 ranks (64 per node), account EAR26024, 2026-10-16 9:33-11:01 pm PDT, 0-12 s,
  110,004 steps, dt 1.09087e-4 s, central difference ("cd"), 5276 s wall (958 s of it fault output),
  max|V| 9.05 m/s, no NaN.
- results/: the 16 SCEC station files (t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 = strike;
  SI; sigma_n compression positive; tpv205_station_conventions.txt), the station snap list, driver.log and
  RUN_INFO.txt. stations.sha256 holds the sha256 of each station file as written on Vista (the files were
  moved in two parts each and checked against it after joining). Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (Vista scratch) is not here.
- plots/: tpv205_100m_p1_vista1076419_stations.pdf and its page PNGs (slip-rate overview of the 16 stations,
  then per station slip rate, front zoom, slip, shear stress, normal stress, dip slip rate; this run vs
  DRDG3D 100 m and 200 m, the mixed-flux p2 gold and the 200 m p1 gold); compare.txt / compare_summary.json
  from compare_nitsche_fault_stations.py vs DRDG3D 100 m, DRDG3D 200 m, both mixed-flux golds and the 200 m
  p1 gold (command below; rerunning it on these results/ reproduces compare.txt).
- vs DRDG3D 100 m (median |misfit| of 16 stations): arrival 0.022 s (early at all 16 but the hypocentre),
  final slip 1.15 %, peak h-slip-rate 8.6 % (low at 15 of 16); noise RMS sigma_n / tau_h / tau_v
  0.020 / 0.017 / 0.020 MPa (DRDG3D 100 m tau_h / tau_v 0.017 / 0.021). Pass test: arrival within 0.1 s
  everywhere PASS; slip within 2 % FAIL at 2 stations (x2 = -12 km 2.7 % low and +12 km 2.2 % low, both at
  the surface); peak within 10 % at 9 of 16 (need 14). Every median is smaller than in the 200 m p1 gold
  (arrival 0.043 s, slip 1.46 %, peak 19.3 %).
- Fault-wide normal stress (driver.log): 38-231 MPa after the rupture reaches the fault edges, against
  120 MPa at rest, and it stays there after slip stops; the stations see 0.02 MPa RMS. Under investigation
  (TPV104 thread notes, 2026-10-17).

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-benchmarks
    cd seas-mfem-benchmarks && git checkout 4803a9000b6dfca1be64a42372517c30e6e7d36d
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build + gate tests
    # copy the 100 m mesh (sha256 above) to miniapps/seas/tpv205/mesh/tpv2053d_100m.msh, then
    sbatch tpv205_100m_p1.sbatch

The job refuses to start unless the checkout is clean, HEAD is on GitHub, the driver was built at HEAD and the
mesh matches meshes.sha256. Compare (from miniapps/seas, with D = this folder):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv205/stations/tpv205_fault_stations.csv \
        --prefix tpv205 --run vista_p1_100m=$D/results \
        --ref drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4 --ref drdg3d_200m=tpv205/benchmark_data/DRDG3D_200m_O4 \
        --gold mixedflux_p1=tpv205/gold/mixedflux_p1_rk45_jobidcombined_7762879_7763019_git04674fc/results \
        --gold mixedflux_p2=tpv205/gold/mixedflux_p2_rk45_jobid7762880_git04674fc/results \
        --gold vista_p1_200m=tpv205/gold/nitsche-kv0p1_p1_cd_pc_200m_jobid1076378_git81759da/results:nitsche \
        --peak-min-stations 14 --json compare_summary.json > compare.txt

- Requested by Chunhui Zhao (TPV205/TPV104 benchmarks, 2026-10-16).
