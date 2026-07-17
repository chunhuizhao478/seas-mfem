# MFEM port of `safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small` — p1 / p2 / p3

**Purpose:** measure MFEM `seas_spatial_dyn_driver` wall-clock **against SeisSol**
on the *same* coarse SAF-ALT mesh, using **pure upwind (Godunov) flux + ADER**
time integration, at DG orders p1, p2, p3.

SeisSol runs this deck at **o4 = p3** (pure Godunov, clustered LTS). The p1/p2/p3
MFEM triplet brackets that on the byte-identical mesh, sidecars, and physics — so
the **only** experimental variable is the discretization (p-order, and MFEM
global-dt ADER-DG vs SeisSol LTS ADER-DG).

## Deliverable layout

Three self-contained folders under `miniapps/seas/jobs/safs/`:

```
safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small_p1/   (order 1, ADER order 2)
safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small_p2/   (order 2, ADER order 3)
safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small_p3/   (order 3, ADER order 4)
```

Each folder contains:

| File | What |
|------|------|
| `spatial_friction_..._CASE1_p{N}.toml` | the MFEM config (order N baked in) |
| `run_safs_v4_0_0_alt_case1_p{N}_upwind_ader_expanse.sbatch` | Expanse submit script |
| `README.md` | one-page pointer to this doc + submit command |
| `mesh_alt.msh` → | symlink to the shared Gmsh v2.2 mesh |
| `velocity_safs.h5` → | symlink to the shared CVM material sidecar |
| `stress_andersonian_k1.65.h5` → | symlink to the shared k=1.65 stress sidecar |
| `friction_thermal_case1.h5` → | symlink to the shared CASE1 friction sidecar |

The four inputs are **identical across the three folders** (the p-order lives only
in the `.toml`), so they are single physical files shared by symlink.

## Input mapping: SeisSol → MFEM

| SeisSol input | MFEM input (this port) | Notes |
|---------------|------------------------|-------|
| `mesh_alt.puml.h5` (1,319,294 tets) | `safalt_dl_safgh_2M_walls_mmg.msh` | Gmsh v2.2, **same** 262,221 nodes / 1,319,294 tets; safstags fault=101 top=102 bottom=103 sides=104 |
| `safs_material_cvm.nc` (ASAGI) | `velocity_safs.h5` (`data_projection_v1`) | Vp/Vs/ρ from ρ,μ,λ; 452×361×194 |
| `safs_stress_andersonian_k1.65.nc` | `stress_andersonian_k1.65.h5` | sign-flipped to compression-**positive**, `--pad-m 3000`; effective (P_p already removed → `P_p_pa=0`) |
| `safs_friction_thermal_case1.nc` | `friction_thermal_case1.h5` | `rs_a`, `rs_srW` (CASE1 zoning) transposed; G4 anchor a=0.015/V_w=0.05 at hypo OK |
| `parameters.par` + `safs_fault.yaml` | `spatial_friction_..._p{N}.toml` | see the config header for the key-by-key mapping |

**Physics is verbatim from the deck:** FL=103 rate-state + strong rate weakening,
f₀=0.6, b=0.019, Dc=0.10, V₀=1e-6, V_init=1e-12, **f_w=0**; SCEC compact-bell
nucleation **75 MPa, R=2000 m, T_nuc=1 s** at the ALT on-fault 10 km hypocentre
`(604446.944, 3704576.3853, -10067.9819)`; CFL 0.5; tfinal 150 s.

## Numerics (the requested combo)

```
mixed_flux      = "none"    # pure Godunov upwind (no central-flux corridor)
time_integrator = "ader"    # ADER-DG (Godunov upwind is stable under ADER — no RK)
interior_flux   = "matrix"  # heterogeneous CVM Riemann solver
ader_order      = order + 1 # p1→2, p2→3, p3→4  (SeisSol ORDER = poly-degree + 1)
```

There is **no `--order` CLI flag** — the polynomial order is the TOML key
`[mesh].order`, which is why there is one config per order. The sbatch asserts
`[mesh].order`/`ader_order` match the file it belongs to.

ADER perf levers `--deriv-cache --shared-ck-recursion` are ON by default (round-off /
byte-identical; the biggest wins in the SAFS caliper study). Set `SAFS_DERIV_OPT=""`
to time the un-optimised ADER path.

## Decisions / deviations from a bit-for-bit SeisSol reproduction (flagged)

1. **Cluster = Expanse (SDSC), 2 exclusive lustre compute nodes, `-t 15:00:00`.**
   Same hardware as the SeisSol baseline (AMD EPYC 7742, `--constraint=lustre`,
   `--account=lbl107`, `--export=ALL`), required for a fair wall-clock comparison.
   SeisSol uses 2N × 8 MPI × 16 OMP = 256 cores (**hybrid**); MFEM is **pure-MPI**
   (no OpenMP in the driver hot path) and uses 2N × 128 ranks = 256 ranks (**full
   core parity**, `--cpus-per-task=1`). Because the CVM/stress/friction sidecars
   are replicated **per MPI rank** (OpenMP threads would share one copy; separate
   ranks each hold their own), pure-MPI at 128 ranks/node needs ~135 GB of sidecar
   copies — so `--mem=249000M` (the exclusive 256 GB node, free) instead of the
   SeisSol deck's 140 GB (which sufficed for its 8 ranks/node). For headroom submit
   with `--ntasks-per-node=64 --mem=249000M` (uses 64/128 cores; the sbatch memory
   guard also suggests this on an OOM-risk sizing).
2. **10 MPa σ_n strength floor (`sigma_n_strength_floor_pa`).** MFEM-only stabiliser
   with **no SeisSol analogue**. The ALT k=1.65 stress is *not* a freeze500m field,
   so shallow daylighting facets can reach low σ_n; under a pure Godunov flux those
   cells can free-slip ("trace-cell runaway": σ_n 0.1–5 MPa cells sliding ~4 m/s,
   which collapses dt and ruins the timing). Matches the established ALT
   pure-upwind+ADER config (`v3_0_0_RSSRW_ALT`). **For strict SeisSol parity, delete
   this one key** from the `[friction]` block.
3. **75 MPa nucleation (deck value).** MFEM's `--print-derived` nucleation heuristic
   is stricter than SeisSol's (it adds the ~8 MPa direct effect at V_init=1e-12).
   The margins clear at the centre (S_E ×1.137, static ×1.171), but if a run only
   **creeps** (no breakout), raise `delta_tau_pa` to ~85e6 before touching radius or
   centre. `--print-derived` is a *warning*, not a hard abort.
4. **Output is light** (fault + free-surface **VTU**, 4 s cadence, no volume/bulk),
   so I/O does not pollute the timing and the run does not depend on `MFEM_USE_HDF5`.

## Build the sidecars (gitignored / large — do this before submitting)

The three `.h5` sidecars are large (velocity 760 MB, stress 256 MB, friction 87 MB)
and gitignored. The symlinks in each folder point to their canonical targets under
`safs/project_7.0_alternative/`. Rebuild them from the SeisSol `.nc` files (run from
`miniapps/seas`, `conda activate pythonenv`):

```bash
SRC=~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small

# velocity/material (452×361×194 CVM → Vp/Vs/ρ)
python3 safs/project_7.0_preferred/velocity/code/build_pref_velocity_sidecar.py \
  --in  "$SRC/safs_material_cvm.nc" \
  --out safs/project_7.0_alternative/velocity/results/multiscale_statewise_cvm_v4_0_0/velocity_safs.h5 \
  --mesh-tag safalt_dl_safgh_2M_walls_mmg

# stress (Andersonian k=1.65; sign-flip + pad 3 km)
python3 safs/seisol_quakeworx/toolbox/on_fault_stress_projection_csm/csm_stress_nc_to_mfem_hdf5.py \
  --in  "$SRC/safs_stress_andersonian_k1.65.nc" \
  --out safs/project_7.0_alternative/stress/results/andersonian_k1p65_mfem/safs_stress_andersonian_k1.65_mfem.h5 \
  --pad-m 3000 --hypo 604446.944 3704576.3853 -10067.9819

# friction (THERMAL CASE1; rs_a / rs_srW; G4 anchor at the ALT hypo)
python3 safs/project_7.0_preferred/friction/code/build_pref_friction_sidecar.py \
  --in  "$SRC/safs_friction_thermal_case1.nc" \
  --out safs/project_7.0_alternative/friction/results/thermal_case1_mfem/friction_safs.h5 \
  --case THERMAL_CASE1 --hypocenter 604446.944 3704576.3853 -10067.9819 \
  --mesh-tag safalt_dl_safgh_2M_walls_mmg
```

On Expanse, build them locally the same way (or scp the local copies) to those
canonical paths inside the repo checkout, then the folder symlinks resolve.

## Submit

```bash
# after `bash build_expanse.sh` has produced ./seas_spatial_dyn_driver
cd miniapps/seas/jobs/safs/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small_p1
sbatch run_safs_v4_0_0_alt_case1_p1_upwind_ader_expanse.sbatch     # then _p2, _p3
```

**Validate first (on Expanse, cheap):** append `--dry-run` (or `--print-derived`,
which continues) to the driver line, or run interactively — it builds the mesh +
sidecars + operators and prints the derived numbers (L_nuc, h_min, dt, on-fault
σ_n) without stepping. This is the real config gate; it was **not** run locally
(the 1.3 M-tet production mesh is cluster-only per project policy). Local checks
done here: TOML parse, all input paths resolve, mesh is Gmsh v2.2 with the
101/102/103/104 safstags, friction grid within bounds + hypo anchor.

## Speed metric

This is a **throughput** comparison, not a bit-for-bit match. From each MFEM log
read the per-step wall-time (steps/s, simulated-s per wall-hour) and compare to the
SeisSol o4 log for the same mesh. The runs checkpoint (`--checkpoint-every 2000`);
they need not reach tfinal=150 s to yield a rate — but for a like-for-like number,
compare over a phase both codes reach (e.g. through nucleation + first ~seconds of
rupture). Expect cost to rise steeply p1→p2→p3 (DOFs/element and the ADER predictor
both grow with order); the **p3** number is the apples-to-apples comparison with
SeisSol's o4.
