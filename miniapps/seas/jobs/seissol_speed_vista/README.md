# SeisSol and MFEM speed test on Vista

Plan: `/mnt/project-files/profiling/seissol-speed-test-2026-10-30/seissol-speed-test-plan-rev1.pdf`
(Chunhui Zhao, 2026-10-30; account EAR24034, SeisSol built in `$WORK/seissol`).

The test times TPV205, TPV104 and TPV102 on the same 200 m Gmsh meshes and the same Vista nodes:

- **MFEM:** the Nitsche solver at p = 1, 2, 3, with central difference and `--lts rate2`.
- **SeisSol:** o2, o3 and o4, with `ClusteredLTS = 2` and CFL 0.5.

Each pair runs in one allocation, MFEM first, then SeisSol, and both write stations only. TPV102 runs on the TPV104 mesh, as its p = 3 gold did.

## Files

| File | What it is |
|---|---|
| `build_speed_test_vista.sbatch` | Step B, one node. It builds MFEM (`build_vista.sh`) and checks the two gold meshes. It builds SeisSol elastic o2, o3 and o4 (sve128, PSpaMM, fp64) with the kit's environment `$WORK/seissol/kit/seissol-env-vista.sh`, pinned by sha256, from source `a8e7f04f`. It builds PUMGen v1.1.3. It converts the meshes to `$WORK/seissol/meshes/*.puml.h5` and writes `MESHES.sha256`. Last, it runs `check_parity.py` for all three problems. |
| `run_pair_vista.sh` | The body of every run job: the guard and input checks (driver commit, mesh and puml sha256, the SeisSol binary in the build manifest, deck parity), then MFEM, then SeisSol in each `SS_LAYOUTS` layout, then `to_scec.py`. |
| `probe_layout.sbatch` | TPV205 o3 on 4 nodes, 0 to 2 s, SeisSol only, one run per entry of `SS_LAYOUTS` (now 4, 4, 8, 16 ranks per node). It picked `SS_LAYOUTS = 8` for the run jobs. |
| `run_<prob>_p<N>.sbatch` | The nine paired runs, 0 to 12 s: p1 on 4 nodes, p2 and p3 on 16. MFEM uses 144 ranks per node, except 128 at p = 3 for TPV104 and TPV102, as their golds. |
| `mfem_decks/<prob>_nitsche_sst.toml` | The gold decks with `paraview_enabled = false`: `tpv205_nitsche.toml`, `tpv104_nitsche.toml` and `tpv102_nitsche_pc_p2_kv0p1.toml`. |
| `decks/<prob>/` | The SeisSol decks: `parameters.par`, fault and material yaml, on-fault and off-fault receivers. |
| `check_parity.py` | Compares each MFEM deck with its SeisSol deck: friction, stresses, nucleation, material, stations, window, output. Fields are compared at about 1,600 to 2,400 fault points with a relative tolerance of 1e-9. It exits 1 on a difference. |
| `to_scec.py` | Converts SeisSol receiver output to our station-file format. |
| `retag_msh2.py` | Renames the TPV205 mesh's boundary tags 101/103/105 to 1/3/5 in a copy for pumgen. SeisSol reads 101 and 105 as rupture faces. |
| `NOTES.md` | The parity tables with file:line on both sides, and the open items. |

## Submit (in `$WORK/seas-mfem-profiling`, at the pushed commit, from this folder)

```
sbatch build_speed_test_vista.sbatch                                   # step B, then:
sbatch --dependency=afterok:<build> probe_layout.sbatch
# after the probe: set SS_LAYOUTS in run_*.sbatch (a reviewed commit), then step C:
for f in run_tpv205_p{1,2,3}.sbatch run_tpv104_p{1,2,3}.sbatch run_tpv102_p{1,2,3}.sbatch; do sbatch $f; done
```

## What the physics does not share (NOTES.md, Open)

- **Damping and Prakash-Clifton.** The MFEM decks are the damped golds' decks: Kelvin-Voigt damping (`kv_eta_bar = 0.1`) and Prakash-Clifton (v* = 1 m/s, L = 0.01 m). SeisSol at `a8e7f04f` has no Kelvin-Voigt term. It has Prakash-Clifton only on the strength in FL = 6, a different law. Each code runs as it is used in production, and the station error against the references shows what that costs.
- **CFL.** The two codes define CFL differently, so their time steps differ.
- **Off-fault stations.** TPV205 has no off-fault stations in MFEM, so `ReceiverOutput = 0` in its SeisSol deck. The TPV104 and TPV102 off-fault receivers have no MFEM Nitsche counterpart.

## Results

| Job | Commit | What | Result |
|---|---|---|---|
| 1112276 | 61484459 | build | FAILED at the mesh step: pumgen writes under the exact name given (fixed in 1831dfa5). |
| 1112915 | 1831dfa5 | build | Passed. Tet counts equal MFEM's (TPV205 1,894,583; TPV104 2,466,549); boundary tags TPV205 1:22648 3:52456 5:7256, TPV104 1:25142 3:75856 5:8006; parity passes for all three. |
| 1112916 | 1831dfa5 | probe, 2/4/8 per node | 2 per node 22.68 s, 8 per node 15.45 s; 4 per node crashed in its first steps (UCX "Local length error on mlx5_0:1/IB" under MPI_Test on the communication thread). |
| 1114522 | 060551aa | probe, 4/4/8/16 per node | All exited 0, no crash: 4 per node 17.77 s and 16.80 s, 8 per node 15.25 s, 16 per node 15.18 s (load imbalance 16.3, 7.9, 6.4, 7.5 %). Binary loads Open MPI 5.0.5 (gcc14) and UCX 1.17.0. 8 and 16 tie, so the runs use 8. Cluster histogram identical in every run; 26,228 rupture faces, as MFEM. |
