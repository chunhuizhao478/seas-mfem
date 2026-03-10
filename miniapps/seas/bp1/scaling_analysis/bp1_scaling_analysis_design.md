# BP1 Scaling Analysis Design

## Context

The BP1 SEAS simulation (`seas_bp1_bdrload`) lacks wall-clock timing output and has no scaling test infrastructure. We need to instrument the driver with `MPI_Wtime()`, create SLURM job scripts for strong/weak scaling on TACC, and build a Python postprocessing pipeline to generate scaling plots.

All scaling tests run the full 3000-year BP1 benchmark to capture complete earthquake cycle performance.

## Design Decisions

- **tfinal**: 3000 years (full BP1 benchmark, all earthquake cycles) — production-level scaling test
- **Timing output**: Single parseable line `SCALING_DATA: np=N elements=E steps=S elapsed=T mesh=M solver=X dg=Y`
- **Strong scaling mesh**: 50m (`bp1_ss_50m.msh`, ~290k elements), fixed across all core counts
- **Weak scaling**: Pre-generate meshes locally via gmsh with `h = 0.30 * sqrt(4/np)`, scp to cluster. ~3k elements/core
- **Unit test**: Python-based validation in the postprocessing script (verify parseable output, elapsed > 0)

## Weak Scaling Mesh Formula

For a 2D triangular mesh, element count scales as `N ~ C / h^2`. For weak scaling we need `N/np = const`, so `N ~ np`, giving `h ~ 1/sqrt(np)`.

Calibrated with `np=4` as base case (h=0.300 km, ~12k elements):

```
h(np) = 0.300 * sqrt(4 / np)
```

## Files to Modify

### 1. `tests/verification/bp1_bdrload.cpp` — Add timing instrumentation

- Add `MPI_Wtime()` before/after the time-stepping loop (lines ~763 and ~888)
- Add `--max-steps` CLI flag (reuse existing `max_steps` variable)
- Print `SCALING_DATA:` structured line in the summary section (~line 910):

```cpp
// Before time-stepping loop:
mpi.Barrier();
double wall_start = MPI_Wtime();

// After time-stepping loop:
mpi.Barrier();
double wall_end = MPI_Wtime();
double elapsed = wall_end - wall_start;

// In summary section:
if (mpi.IsRoot())
{
   std::cout << "\nSCALING_DATA:"
             << " np=" << mpi.Size()
             << " elements=" << global_ne
             << " steps=" << step
             << " elapsed=" << std::fixed << std::setprecision(2) << elapsed
             << " mesh=" << mesh_file
             << " solver="
#ifdef MFEM_USE_MUMPS
             << "MUMPS"
#else
             << "CG_ILU"
#endif
             << " dg=" << dg_method_str
             << "\n";
}
```

## Files to Create

### 2. Strong Scaling Jobs: `jobs/bp1_strongscaling/`

8 sbatch files for np = 1, 2, 4, 8, 16, 32, 40, 80:

- Fixed mesh: `bp1/mesh/bp1_ss_50m.msh`
- Full 3000-year simulation (default tfinal), `--checkpoint-interval 5000`
- np<=40: 1 node; np=80: 2 nodes
- Wall time: 48h (full earthquake cycle simulation)
- Output to `bp1/scaling_analysis/strong_scaling/npN/`
- Plus `submit_all.sh` launcher

**Note:** With 3000-year tfinal on the 50m mesh, np=1 and np=2 may need multi-day runs. Consider starting with np>=4 if allocation is limited.

### 3. Weak Scaling Jobs: `jobs/bp1_weakscaling/`

6 sbatch files for np = 4, 8, 16, 32, 40, 80:

| np | Nodes | h (km) | ~Elements | Elem/core |
|----|-------|--------|-----------|-----------|
| 4  | 1     | 0.300  | ~12k      | ~3k       |
| 8  | 1     | 0.212  | ~24k      | ~3k       |
| 16 | 1     | 0.150  | ~48k      | ~3k       |
| 32 | 1     | 0.106  | ~96k      | ~3k       |
| 40 | 1     | 0.095  | ~120k     | ~3k       |
| 80 | 2     | 0.067  | ~240k     | ~3k       |

- Meshes pre-generated locally and scp'd to cluster (no gmsh on compute nodes)
- Full 3000-year simulation (default tfinal), `--checkpoint-interval 5000`
- Output to `bp1/scaling_analysis/weak_scaling/npN/`
- Plus `submit_all.sh` launcher

### 4. Local Mesh Generation: `jobs/bp1_weakscaling/generate_weak_meshes.sh`

- Runs locally (not on cluster) with `conda activate pythonenv`
- Generates meshes: `bp1_ss_weak_np{4,8,16,32,40,80}.msh` in `bp1/mesh/`
- Uses `gmsh -2 bp1/mesh/bp1_selfsimilar.geo -o OUTFILE -setnumber h H` for each h value
- User scp's the generated meshes to cluster before submitting jobs

### 5. Postprocessing: `bp1/scaling_analysis/scaling_analysis.py`

- Parse `SCALING_DATA:` from SLURM `.out` files via regex
- Compute strong scaling: speedup = T1/Tn, efficiency = speedup/(n/n0)
- Compute weak scaling: efficiency = T_base/Tn
- Generate plots:
  - **Strong scaling** (3-panel): wall time vs cores, speedup vs cores (+ ideal line), parallel efficiency vs cores
  - **Weak scaling** (2-panel): wall time vs cores (+ ideal flat line), weak efficiency vs cores
- Save CSV summaries alongside plots
- Include `--validate` flag that checks all outputs are parseable and elapsed > 0 (unit test)

## Verification

1. Build: `make seas_bp1_bdrload` — ensure compilation succeeds with timing changes
2. Local test: `mpirun -np 2 ./seas_bp1_bdrload --mesh bp1/mesh/bp1_ss_50m.msh --tfinal 1e5 --checkpoint-interval 0` — quick run to verify `SCALING_DATA:` line appears in stdout
3. Generate weak scaling meshes locally: `bash jobs/bp1_weakscaling/generate_weak_meshes.sh`, then scp to cluster
4. Submit scaling jobs on TACC, collect `.out` files
5. Run `python bp1/scaling_analysis/scaling_analysis.py --strong-dir DIR --weak-dir DIR --output-dir plots`
6. Run `python bp1/scaling_analysis/scaling_analysis.py --validate --strong-dir DIR --weak-dir DIR` to verify all outputs
