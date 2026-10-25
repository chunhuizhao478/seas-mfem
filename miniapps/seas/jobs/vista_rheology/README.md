# Vista jobs for TPV26/TPV27 and rheology/

Branch `nitsche-dynamic-rupture-nonlinearity` (Chunhui Zhao, 2026-10-20). The jobs carry out the Vista
gates of the rheology plan rev 1 (`/mnt/project-files/nonlinearity/rheology-plan-rev1.pdf`, sections 10.3
and 11) on account EAR26024, partition `gg`, 144 cores per node. Phase R1: G1 = U8, U11, B0, B1, B2.
Phase R3: G3 = O1-O7 (here and in the cloud) and B3 (cloud), plus the first run of the rheology output,
which needs parallel HDF5.

## Files

| File | What it does |
|---|---|
| `build_tpv2627_vista.sbatch` | 1 node. The CLAUDE.md guard; `build_vista.sh` (writes `BUILD_COMMIT_VISTA`), then `seas_test_nitsche_lts` and the rheology unit tests; runs U1-U11 (`make test-rheology`), the six Nitsche gate test sets, the LTS suite, the operator tests O1-O7 (`make test-nitsche-rheology`, 1, 2 and 4 ranks) and the station output and mesh parts tests (`make test-station-output`, `make test-mesh-parts`); makes the TPV26/27 1000 m and 200 m mirror meshes (`tpv26/mesh/make_tpv2627_meshes.sh`, sha256 to `MESH_SHA256_TPV2627_VISTA`); copies the TPV205 and TPV104 200 m gold meshes from `$WORK/seas-mfem-benchmarks` and checks them against `jobs/vista_profiling/meshes.sha256`; the mesh parts test on the TPV26/27 200 m mirror mesh (16 ranks, fault attribute 103); B1 (TPV26 1000 m p=1, 0 to 2.5 s, 8 ranks: no abort and rupture beyond r_crit = 4 km); the R3 output check (the TPV27 deck, 1000 m p=1, 0 to 1.0 s, 8 ranks: `rheology.vtkhdf` with 3 snapshots); a 200 m p=2 load check on 144 ranks (memory in use and committed, tets, dt; informational, it does not gate). |
| `b0_tpv205_200m_p1.sbatch`, `b0_tpv205_200m_p2.sbatch`, `b0_tpv104_200m_p1.sbatch`, `b0_tpv104_200m_p2.sbatch` | B0: the gold's deck, mesh and rank layout (1076378, 1075916, 1075919, 1075917), the defaults except `--mesh-setup every-rank` (the LTS gold jobs' mesh setup), 0 to 1.0 s. Pass when the run reproduces the trunk code's 1 s run of the LTS gold job beside that gold: the same `=== done` line as its `gts_driver.log` and the same largest station difference from the gold (`compare_stations.py`). |
| `run_b0_vista.sh` | Shared body of the B0 jobs. |
| `tpv26_200m_p2.sbatch` | B2: TPV26 on the 200 m mirror mesh, p = 2, 0 to 13.0 s, 24 nodes x 128 ranks, rank-0 mesh setup (the default; logs the first node's committed memory), through `jobs/vista_nitsche/run_nitsche_vista.sh` (output `$SCRATCH/nitsche_benchmarks/tpv26_mirror200m_p2_<job>`). Compared with DRDG3D TPV26 over 0 to 12 s off Vista; the rupture-time contour file comes from `scripts/make_cplot.py`. |

TPV102 has no Vista gold (its golds ran on Frontera), so its B0 is the cloud check only: the TPV205,
TPV104 and TPV102 decks on 1000 m meshes, 3 ranks, 0 to 0.6 s, every output file byte-identical to the
trunk code's.

## Submit

Each job is submitted only on Chunhui Zhao's word for it. On the login node:

```
cd $WORK && git clone --recurse-submodules -b nitsche-dynamic-rupture-nonlinearity \
    https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-nonlinearity
cd seas-mfem-nonlinearity/miniapps/seas/jobs/vista_rheology
sbatch build_tpv2627_vista.sbatch
# after it prints "all checks passed":
sbatch b0_tpv205_200m_p1.sbatch; sbatch b0_tpv104_200m_p1.sbatch
sbatch b0_tpv205_200m_p2.sbatch; sbatch b0_tpv104_200m_p2.sbatch
sbatch tpv26_200m_p2.sbatch
```

Cost at most: build 0.66 SU, B0 0.66 + 0.66 + 2.64 + 2.64 SU, B2 15.8 SU (estimate 7 SU).

Memory: Vista compute nodes use strict overcommit (mode 2, ratio 50, no swap; read on i617-101 by the
benchmarks thread, 2026-10-21), so a node can commit at most 118.8 GiB, half its RAM, and past that every
allocation fails without a process being killed.  Up to 2f4ede1 every rank of the Nitsche driver read the
whole serial mesh and partitioned it.  Since the merge of the Nitsche profiling branch (bc91ac9) the
default `--mesh-setup rank0` has rank 0 read and partition the mesh and send each rank its part
(`--mesh-setup every-rank` is the old path); on TPV205 200 m p=2, 16 x 144, the first node's peak
committed memory went from 114.1 GiB to 37.6 GiB (profiling job 1090566).

## Runs

- Build 1085864 (babe089, 2026-10-21, node i617-051, 7 min 46 s): the build and every test passed
  (rheology U1-U11 24/13/36/25, the six Nitsche gate sets, LTS 24, O1-O7 108 checks), then the job stopped
  at the mesh step: the TACC gmsh module (4.15.0-git) has no OpenCASCADE, which the first
  `tpv2627_half.geo` used, so gmsh wrote an empty mesh.  The .geo now uses gmsh's built-in kernel.  The
  four B0 jobs 1085865-1085868 (afterok on it) never started; they were cancelled.
- Build 1086883 (2f4ede1, 2026-10-21, node i618-023, 11 min 51 s): all tests passed (as in 1085864); Vista
  gmsh made the meshes (1000 m: 242,154 tets, 1,858 fault faces; 200 m: 2,995,976 tets, 46,328 fault faces,
  sha256 8a72ea76...); the TPV205 and TPV104 gold meshes matched; B1 passed (54 of 54 nodes inside r_crit and
  95 of 936 outside ruptured by 2.5 s, the farthest 6.8 km away); the R3 output check passed (rheology.vtkhdf
  with 3 snapshots).  The 200 m load check printed the mesh and dt (1.330e-4 s, 97,755 steps to 13 s) and
  was then killed (signal 9, 202.5 GiB in use of 237.6 GiB): the whole problem on one node does not fit.
- B0 1086884-1086887 (2f4ede1): all four PASS (the same done line and the same largest station difference
  from the gold as the trunk code's 1 s run: TPV205 p2 2.691e-11, TPV104 p2 1.225e-11).
- TPV26 1087808 (2f4ede1, 24 nodes x 128 ranks): stopped at setup after 59 s, before the mesh summary,
  with std::bad_alloc on 21 ranks, all on node i617-101; driver.log in
  $SCRATCH/nitsche_benchmarks/tpv26_mirror200m_p2_1087808.  The node commit limit above, with every rank
  holding the serial mesh, is the likely cause (inferred); the rerun uses rank-0 mesh setup.
- Build 1094996 and TPV26 1094998 (b618bbb5, 2026-10-24): the build passed; TPV26 ran 0 to 13 s in 22 min
  (12.8 ms/step), the first node's peak committed memory 56.6 GiB of 118.8.  Against DRDG3D (station PDF in
  /mnt/project-files/nonlinearity/tpv26-200m-p2-stations-2026-10-24/): median arrival difference 0.023 s,
  slip 0.30 %, peak slip rate 4.6 %; sigma_n 0.3 to 1 % above DRDG3D at five stations, whose station files
  held a fault node 50 m deeper than the SCEC point (sigma_n grows 17.6 MPa per km).  The .geo now makes
  the 12 stations mesh vertices.
