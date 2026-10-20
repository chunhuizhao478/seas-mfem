# Vista jobs for TPV26/TPV27 and rheology/

Branch `nitsche-dynamic-rupture-nonlinearity` (Chunhui Zhao, 2026-10-20). The jobs carry out the Vista
gates of the rheology plan rev 1 (`/mnt/project-files/nonlinearity/rheology-plan-rev1.pdf`, sections 10.3
and 11) on account EAR26024, partition `gg`, 144 cores per node. Phase R1 (this folder so far): G1 = U8,
U11, B0, B1, B2.

## Files

| File | What it does |
|---|---|
| `build_tpv2627_vista.sbatch` | 1 node. The CLAUDE.md guard; `build_vista.sh` (writes `BUILD_COMMIT_VISTA`), then `seas_test_nitsche_lts` and the rheology unit tests; runs U1-U11 (`make test-rheology`), the six Nitsche gate test sets and the LTS suite; makes the TPV26/27 1000 m and 200 m mirror meshes (`tpv26/mesh/make_tpv2627_meshes.sh`, sha256 to `MESH_SHA256_TPV2627_VISTA`); copies the TPV205 and TPV104 200 m gold meshes from `$WORK/seas-mfem-benchmarks` and checks them against `jobs/vista_profiling/meshes.sha256`; B1 (TPV26 1000 m p=1, 0 to 2.5 s, 8 ranks: no abort and rupture beyond r_crit = 4 km); a 200 m p=2 load check on 144 ranks (memory, tets, dt). |
| `b0_tpv205_200m_p1.sbatch`, `b0_tpv205_200m_p2.sbatch`, `b0_tpv104_200m_p1.sbatch`, `b0_tpv104_200m_p2.sbatch` | B0: the gold's deck, mesh and rank layout (1076378, 1075916, 1075919, 1075917), the defaults, 0 to 1.0 s. Pass when the run reproduces the trunk code's 1 s run of the LTS gold job beside that gold: the same `=== done` line as its `gts_driver.log` and the same largest station difference from the gold (`compare_stations.py`). |
| `run_b0_vista.sh` | Shared body of the B0 jobs. |
| `tpv26_200m_p2.sbatch` | B2: TPV26 on the 200 m mirror mesh, p = 2, 0 to 13.0 s, 24 nodes x 128 ranks, through `jobs/vista_nitsche/run_nitsche_vista.sh` (output `$SCRATCH/nitsche_benchmarks/tpv26_mirror200m_p2_<job>`). Compared with DRDG3D TPV26 over 0 to 12 s off Vista; the rupture-time contour file comes from `scripts/make_cplot.py`. |

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

## Runs

None yet.
