# TPV31 gold: Nitsche dynamic rupture, Kelvin-Voigt damping, p1, 200 m (Vista job 1085893)

- Code: seas-mfem branch nitsche-benchmarks, commit 960771241718762466c6802a581184fbea37b7dc (9607712, on GitHub),
  built on TACC Vista by jobs/vista_nitsche/build_verify_vista.sbatch (job 1084231, clean checkout,
  BUILD_COMMIT_VISTA = 9607712). 9607712 is trunk a4c9dee plus the p3 jobs and golds (b50741b, daae761), the TPV6/TPV7/TPV31
  setup (652e976), the gmsh built-in kernel meshes (d910b2f) and the run-layout memory checks; the Nitsche speedups are on by default (the
  run logs "fast fault on, fused K/K_eta off, dense blocks on"). Driver seas_nitsche_dyn_driver.
- Deck and job script: jobs/ (= tpv31/configs/tpv31_nitsche.toml and jobs/vista_nitsche/{tpv31_200m_p1.sbatch,
  run_nitsche_vista.sh} at 9607712). SCEC TPV31: layered 1D medium ([material] depth_profile_1d, jumps at 2.4, 5
  and 10 km depth, linear between 2.4 and 5 km, sampled once per element at its centroid); prestress scaled by
  mu/mu_0 ([stress_scale_by_mu]; sigma_n 60 MPa and along-strike shear 30 MPa at mu_0 = 32.04 GPa); linear slip
  weakening (mu_s 0.58, mu_d 0.45, d_c 0.18 m) with cohesion 0.000425 MPa/m above 2.4 km depth; instantaneous
  circular overstress 4.95 MPa mu/mu_0 within 1.4 km of the hypocentre at 7.5 km depth, cosine taper to 2 km;
  Prakash-Clifton v* 1 m/s, L 0.01 m; [numerics] kv_eta 0.1, kv_scale "node", kv_length "inscribed"; cfl 0.5.
- Mesh: tpv31/mesh/tpv31_nitsche_200m.msh, made on Vista by build_verify_vista.sbatch from the committed
  tpv31/mesh/tpv31_nitsche_200m.geo (gmsh 4.15.0, built-in geometry kernel, msh 2.2): 2,029,790 tets, 26,248 fault
  faces, sha256 5bd0b52ddd5080c531e2fc1e96b7b5d86edee0b8957c55f11e7d5294c02b9a71 (not in the repo). The layer
  depths are mesh surfaces. The mesh is not mirror-symmetric about the fault.
- Run: Vista gg, 4 nodes / 576 ranks (144 per node), account EAR24034, 2026-10-21, 9:41 job wall (about 0.2 SU),
  0-15 s, 58,551 steps, dt 2.561883944e-4 s, central difference ("cd"), 559.55 s driver wall (62.0 s of it fault
  output), max|V| 2.256 m/s, no NaN.
- results/: the 30 SCEC station files (strike 0, 6, 12 km x depths 0, 0.2, 0.5, 1, 2.4, 3, 5, 7.5, 10, 12 km;
  t slip1 slip2 V1 V2 tau1 tau2 sigma_n slip_path; 1 = dip, 2 = strike; SI; sigma_n compression positive;
  tpv31_station_conventions.txt), the station snap list (stations snap to the nearest fault node, up to 114 m
  from the requested point), driver.log, RUN_INFO.txt and the job's stdout tpv31_200m_p1_1085893.out, the source
  of the speed numbers below. stations.sha256 holds the sha256 of each station file as written on Vista; the
  copies here were checked against it. Read with
  `scripts/compare_nitsche_fault_stations.py --gold LABEL=<this folder>/results:nitsche`.
  fault.vtkhdf (300 frames, Vista scratch) is not here.
- plots/: tpv31_200m_p1_vista1085893_stations.pdf and its page PNGs (misfit maps of arrival, final slip and peak
  slip rate over the 30 stations, then per station slip rate, front zoom, slip, shear stress, normal stress and
  shear/normal ratio; this run in red, drawn last, vs SeisSol in black and EQdyna in gray); compare.txt /
  compare_summary.json from compare_nitsche_fault_stations.py vs SeisSol and EQdyna (command below). SeisSol ends
  at 14.9964 s, so every source is measured to there. compare.txt equals the job's compare_vs_scec.txt on Vista
  except for the pass-3 line, which now gives the 24-station requirement, and the last line (the JSON path).
- vs SeisSol (median |misfit| of 30 stations): arrival 0.113 s (early at all 30), final slip 1.73 %, peak
  h-slip-rate 21.8 % (low at 28 of 30); noise RMS sigma_n / tau_h / tau_v 0.006 / 0.003 / 0.005 MPa (SeisSol
  tau_h / tau_v 0.004 / 0.005). Fails pass tests 1-3: arrival > 0.1 s at 17 stations (all 10 at 12 km along strike,
  7 at 6 km; worst -0.314 s at x2 = 12 km, 2.4 km depth), slip > 2 % at 14 (2 to 5 % low: 9 of the 10 stations at 12 km
  along strike, and the deepest stations at 0 and 6 km), peak within 10 % at 5 of 30 (need 24; worst -46.1 % at x2 = 12 km, 2.4 km depth).
  The error grows along strike: median arrival -0.035 / -0.113 / -0.204 s and peak -12 / -22 / -32 % at 0 / 6 /
  12 km. For scale, EQdyna vs SeisSol: max |dt| 0.030 s, max |slip| 2.02 %, peak within 10 % at 25 of 30.
- Initial stresses match EQdyna within 1.6 % at every station. SeisSol's initial shear and normal stress differ
  from both by up to 15 % at the stations on the material interfaces (2.4, 5 and 10 km depth), higher at some
  and lower at others, with the same shear/normal ratio: read here as the side of the interface SeisSol samples
  (inferred from the station files, not from SeisSol's setup).
- Fault-wide normal stress (driver.log): 24.5-81.1 MPa at rest (the mu/mu_0 scaling), 14.5-100.1 MPa over the
  run, reached by t = 7.7 s and held to 15 s.

## Speed

From results/tpv31_200m_p1_1085893.out and driver.log; job wall from sacct; Vista gg charges 0.33 SU per
node-hour.

| | p1, this run (9607712) |
|---|---|
| steps to 15 s, dt | 58,551, 2.562e-4 s |
| job wall (sacct) | 9:41 |
| driver wall | 559.55 s (9.56 ms per step) |
| time in ode.Step per rank, slowest | 428.7 s = 7.32 ms per step |
| time in ode.Step per rank, min / mean | 5.80 / 6.80 ms per step (max / mean 1.08) |
| step 1 (with operator assembly) | 0.92 s |
| setup plus the first 1000 steps | 65.5 s (driver.log wall at step 1000) |
| fault work, busiest rank | 1.06 ms per step |
| fault output (rank 0) | 300 frames, 62.0 s (0.207 s per frame) |
| node-hours, SU | 0.65, 0.21 |
| first-node memory peak | 141.6 of 237.6 GB |
| output folder on Vista | 2.7 GB |

- The driver's "time in ode.Step" excludes setup; the slowest rank sets the pace.
- For scale: the TPV205 200 m p1 gold (1076378, 1.89 M tets, same 4 x 144 layout, before the speedups) took
  26.5 ms of driver wall per step.

## Reproduce

On Vista (login node for the network steps):

    cd $WORK && git clone --recurse-submodules -b nitsche-benchmarks \
        https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-tpv6731
    cd seas-mfem-tpv6731 && git checkout 960771241718762466c6802a581184fbea37b7dc
    cd miniapps/seas/jobs/vista_nitsche && sbatch build_verify_vista.sbatch   # build, tpv31_nitsche_200m.msh, gate tests
    sbatch tpv31_200m_p1.sbatch

build_verify_vista.sbatch generates tpv31/mesh/tpv31_nitsche_200m.msh from the committed .geo with Vista's gmsh
4.15.0; gmsh output depends on the platform, so compare its sha256 with the one above (the run script refuses a
mesh whose sha256 is not in MESH_SHA256_VISTA). The job refuses to start unless the checkout is clean, HEAD is on
GitHub and the driver was built at HEAD. Compare (from miniapps/seas, with D = this folder; rerunning it on these
results/ reproduces compare.txt and compare_summary.json; checked 2026-10-21):

    python3 -I scripts/compare_nitsche_fault_stations.py --stations tpv31/stations/tpv31_fault_stations.csv \
        --prefix tpv31 --run p1_200m=$D/results \
        --ref seisol=tpv31/benchmark_data/scec_seisol --ref eqdyna=tpv31/benchmark_data/scec_eqdyna \
        --peak-min-stations 24 --t-end 15 --json $D/plots/compare_summary.json > $D/plots/compare.txt

- Requested by Chunhui Zhao (TPV6/TPV7/TPV31 benchmarks, 2026-10-21).
