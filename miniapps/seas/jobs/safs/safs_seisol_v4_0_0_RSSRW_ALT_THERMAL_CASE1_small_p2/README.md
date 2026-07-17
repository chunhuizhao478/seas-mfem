# safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small_p2

MFEM port of the SeisSol deck **safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small**,
run at **DG order p2** (ADER temporal order 3) with **pure upwind (Godunov)
flux + ADER** time integration. One of the p1/p2/p3 speed-benchmark triplet.

See the full writeup (physics, input mapping, decisions, sidecar build commands,
speed methodology) in:
`../README_safs_seisol_v4_0_0_alt_case1_mfem_speed.md`

## Contents
- `spatial_friction_safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_p2.toml` — config ([mesh].order = 2)
- `run_safs_v4_0_0_alt_case1_p2_upwind_ader_expanse.sbatch` — Expanse submit script
- `mesh_alt.msh`, `velocity_safs.h5`, `stress_andersonian_k1.65.h5`, `friction_thermal_case1.h5` —
  symlinks to the shared inputs under `safs/project_7.0_alternative/` (gitignored; rebuild via the parent README).

## Submit (from this folder, after `build_expanse.sh`)
```bash
sbatch run_safs_v4_0_0_alt_case1_p2_upwind_ader_expanse.sbatch
```
Numerics: `--mixed-flux none --time-integrator ader --ader-order 3` (+ `--deriv-cache --shared-ck-recursion`).
The sbatch asserts the config's `[mesh].order`==2 and `ader_order`==3, and fails loudly on any dangling input symlink.
