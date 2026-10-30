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
| `probe_layout.sbatch` | TPV205 o3 on 4 nodes, 0 to 2 s, SeisSol only, at 2, 4 and 8 ranks per node. It picks `SS_LAYOUTS` for the run jobs (a follow-up commit). |
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
