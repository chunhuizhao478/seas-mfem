# Implementation Plan: the SeisSol v3.4.1 RSSRW PREFERRED THERMAL CASE2 problem in SEAS-MFEM, with mixed flux

Date: 2026-07-08
Author: exploration + plan pass over `seas-mfem-spatial-dyn-driver`, branch `system/spatial_dyn_driver`

Successor to the v3.0.0 port plan (2026-06-19), in `project_7.0_alternative/document/`.

---

## Overview

Run the **same physical problem** as the SeisSol deck

```
~/Downloads/seisol_quakeworx/safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2/
```

through `seas_spatial_dyn_driver` with MFEM's **mixed (adjacent-central) flux + explicit RK**, which
SeisSol does not have. Same mesh, same CVM material, same Andersonian k=1.8 initial stress, same
temperature-zoned rate-and-state + strong-rate-weakening friction, same 75 MPa / R=2000 m compact
nucleation bump, same 150 s of simulated time.

This is **not** a bit-for-bit reproduction: SeisSol runs ADER-DG order 4 (p=3) with clustered LTS and a
pure Godunov flux; MFEM runs p=1 global time stepping with a central flux on fault-adjacent faces. The
*inputs* transfer exactly; the *discretization* is the experimental variable.

### Headline

**Four of the six input channels already work.** The mesh is already in the repo in MFEM-readable form and
is provably the same mesh SeisSol ran. Stress and material have existing `.nc → data_projection_v1 HDF5`
converters that need a re-run, not a rewrite. The flux/integrator combination, the FL=103 friction law, and
the compact-bell nucleation are all implemented and unit-tested.

**Two things genuinely block the run, and both are small:**

| # | Blocker | Fix size |
|:--|:------------------------------------------------------------|:----------------------|
| **B1** | No way to feed a 3-D `rs_a(x,y,z)` / `rs_srW(x,y,z)` field to the friction resolver. The thermal deck's entire novelty (CTM temperature zoning) lives in that field. | ~120 lines C++ + ~150 lines Python + 3 tests |
| **B2** | `f_w_default` parser guard rejects `f_w = 0.0`, which is exactly what `RS_muW = 0.0` in the deck asks for. | 1 line + 1 test |

Everything else is data preparation, config authoring, and one optional new output writer.

---

## Background — verified findings

Every claim below was checked against the source or the files themselves during this exploration.

### The deck (SeisSol side)

`parameters.par` + `safs_fault.yaml` + `safs_initial_stress.yaml` + `safs_material_cvm.yaml`:

| Quantity | Value |
|---|---|
| Friction law | `FL = 103` (rate-and-state slip law + strong rate weakening) |
| `RS_f0` / `RS_sr0` | `0.6` / `1e-6` |
| `RS_b` | `0.019` (constant — v1.1.3 does not read spatial `rs_b`) |
| `RS_muW` (`f_w`) | **`0.0`** |
| `rs_sl0` (`Dc`) | `0.10` m (constant) |
| `RS_iniSlipRate1/2` | `1e-12` / `0.0` |
| `rs_a`, `rs_srW` | **spatial**, from `safs_friction_thermal_case2.nc` (ASAGI) |
| Nucleation | `Tnuc_s = 75e6 · exp(r²/(r²−R²))`, `R = 2000 m`, centre `(609062.8722, 3709528.1324, −10000.0)`, `Tnuc_n = Tnuc_d = 0` |
| Nucleation ramp | `s_0 = 0`, `t_0 = 1.0 s` (SeisSol smoothStep) |
| Fault reference | `XRef,YRef,ZRef = (0,−1,0)`, `refPointMethod = 1` (reference **vector**) |
| Initial stress | `safs_stress_andersonian_k1.8.nc`, Andersonian C1 k=1.8, hydrostatic Pp **already removed** (effective), compression **NEGATIVE** |
| Material | `safs_material_cvm.nc`, compound `{rho, mu, lambda}` |
| Mesh | `mesh_preferred.puml.h5`, 1,241,419 tets / 250,588 vertices, 172,593 fault facets |
| Discretization | `CFL = 0.5`, `ClusteredLTS = 2`, order 4 (compile-time: `SeisSol_Release_drome_4_elastic`) |
| Flux | `numflux = numfluxnearfault = 'godunov'` (pure upwind everywhere) |
| `EndTime` | `150.0` s |
| Plasticity | off |

Sidecar grids (all ASAGI: NETCDF4, one compound var `data`, dims `(z,y,x)`, float32 members, float64
strictly-increasing equidistant axes, z = elevation):

| file | dims (z,y,x) | Δx, Δy, Δz [m] | extent [m] | members |
|---|---|---|---|---|
| friction (CASE2) | 107 × 194 × 263 | 1500, 1500, 200 | x[303,696]k y[3612,3901.5]k z[−21000,+200] | `rs_a` ∈ [0.015, 0.0332]; `rs_srW` ∈ [0.05, 1000] |
| material (CVM) | 194 × 361 × 452 | 1500, 1500, 250 | x[28.5,705]k y[3543,4083]k z[−45000,+3250] | `rho`, `mu`, `lambda` |
| stress (Anders. k=1.8) | 93 × 281 × 434 | 1000, 1000, 250 | x[195,628]k y[3685,3965]k z[−20000,+3000] | `s_xx..s_xz`; **`s_yz ≡ s_xz ≡ 0`** |

CASE2 friction profile (b = 0.019 constant, slope m = 8.0e−5 /°C):

```
a − b  =  −0.004                       for T ≤ 300 °C     (a = 0.015; VW everywhere shallow)
a − b  =  −0.004 + m·(T − 300)         for T ≥ 300 °C     (zero-crossing exactly at 350 °C)
V_w    =  0.05 m/s                     for T ≤ 350 °C
V_w    =  0.05 + (1000−0.05)·(T−350)/50   for 350 < T < 400 °C
V_w    =  1000 m/s                     for T ≥ 400 °C
```

CASE1 (the sibling deck) adds a shallow VS cap at the 100 °C isotherm; **CASE2 has no shallow VS zone** —
the shallow fault keeps `V_w = 0.05` and can weaken dynamically. At the hypocentre column: 100 °C at
2.93 km, 350 °C at 13.19 km depth. `T(hypocentre) = 277.8 °C` ⇒ `a = 0.015`, `V_w = 0.05`.

### The mesh is already here, and it is the same mesh

```
$ shasum -a 256 .../meshing_deep19km/results/safv4_deep19km_sub500m_flattop_nwtrim_fixed_opt.puml.h5
                ~/Downloads/seisol_quakeworx/.../mesh_preferred.puml.h5
f08aef34...  (identical)
```

The PUML the SeisSol job read was produced by `msh_to_puml.py` from

```
safv4_deep19km_sub500m_flattop_nwtrim_fixed_opt_safstags.msh
```

whose default map is `101 fault → 3`, `102 top → 1`, `103 bottom → 5`, `104 sides → 5`. That `.msh` is
Gmsh **v2.2 ASCII** (MFEM's only supported flavour), 250,588 nodes, 1,543,222 elements (1,241,419 tets
plus 301,803 boundary triangles), physical tags `101/102/103/104` + volume `1`.

`seas_spatial_dyn_driver` with **no** `[boundary]` block falls back to exactly
`fault_attr=101, natural={102}, absorbing={103,104}` (`spatial_dyn_driver.cpp:1070`). So the
mesh channel is **zero work**.

For scale: the MFEM SAFS jobs already in production run `safv4_deep_500m_opt_safstags.msh` — **314,284
nodes / 1,846,898 elements**, i.e. *larger* than the deep19km mesh — on 2 nodes × 50 ranks.

### What MFEM already implements

- **`interior_flux="matrix"` + `mixed_flux="adjacent"` + `time_integrator="rk4"|"rk45"`** is the intended
  heterogeneous-CVM mixed-flux mode. `BimaterialWaveOperator::SetMixedFluxMode`
  (`bimaterial_wave_operator.inl:1170`) explicitly *lifts* the old R-003 matrix×mixed-flux
  exclusion. Central flux is non-dissipative ⇒ **unstable under ADER**; both `ComputeMaxDt` overloads and
  a driver-level guard (`MatrixMixedFluxUnderAder`) abort on that combination. The comment at
  `wave_operator.hpp:335` claiming "mixed flux is scalar-only" is **stale**.
- **FL=103 equivalent**: `state_evolution = "slip_law_strong_rate_weakening"` (alias `"slip_law_srw"`),
  `slip_law_srw_psi.hpp`:
  `f_ss(V) = f_w + (f_LV − f_w)/(1+(V/V_w)⁸)^{1/8}`, `f_LV = max(0, f_0 − (b−a)·ln(V/V_0))`,
  `ψ_ss = a·ln((2V_0/V)·sinh(f_ss/a))`. Per-QP `a` and per-QP `V_w` are threaded through
  `Rate_SRW` / `SteadyState_SRW` / `PsiSS_SRW`.
- **Compact-bell nucleation**: `[nucleation] kind = "gradual_overstress_compact_circular"` →
  `CompactBellFactor(r,R) = exp(r²/(r²−R²))`, strike-only (`tau2_nuc`), SCEC smoothStep temporal ramp over
  `[0, T_nuc_s]` — the same function SeisSol uses for `t_0`. Parser at
  `spatial_friction.cpp:1363`, driver dispatch at `spatial_dyn_driver.cpp:1911`.
  **Both `CLAUDE.md` (§ *Spatial driver nucleation mechanism*) and
  `spatial_friction_config_schema.md` claim the spatial driver has only
  `gradual_overstress`, and that rate-state ignores `[nucleation]`.
  Both statements are stale** — `tpv104/configs/tpv104_spatial.toml`
  is `law="rate_state"` + `kind="gradual_overstress_compact_circular"` and is a validated benchmark.
  Fix both docs as part of this work.
- **Stress sidecar**: `[stress] kind="sidecar_hdf5"` → `StressField3D` (`stress_field_3d.hpp`), six
  `DataField3D` readers over a `data_projection_v1` HDF5 (`grid/{x,y,z}` float64 ascending,
  `fields/sigma_{xx,yy,zz,xy,yz,xz}` shaped `(Nx,Ny,Nz)`, attrs `schema_version`,
  `crs="EPSG:32611"`, `z_positive="elevation"`, per-dataset `min_value`/`max_value`).
  **Strict pass-through — compression POSITIVE.** The `σ_xy` negation at
  `spatial_dyn_driver.cpp:1716` is in the `ConstantTensor` branch **only**; the sidecar branch
  does not touch signs.
- **Material sidecar**: `[velocity]` bundle reads `fields/{Vp, Vs, density}` from the same schema;
  `far_field_clamp = true` ⇒ `OOBPolicy::Clamp`, the ASAGI nearest-edge hold. `μ = ρVs²`,
  `λ = ρVp² − 2μ`. Evaluated per quadrature point (`Mode::Coefficient`), no P1 smearing.
- **Converters already written**:
  - `seisol_quakeworx/toolbox/on_fault_stress_projection_csm/` `csm_stress_nc_to_mfem_hdf5.py`
    — ASAGI stress `.nc` → schema-v1 HDF5, applies the **global −1 sign flip** on all six components and
    already nudges `max_value` for identically-zero components (which `s_yz`/`s_xz` are here).
  - `safs/project_7.0_preferred/velocity/code/build_pref_velocity_sidecar.py`
    — `safs_material_cvm.nc {rho,mu,lambda}` → schema-v1 HDF5 `{Vp,Vs,density}`.
- **Checkpoint/restart**, ParaView fault output (7 fields: `slip_{dip,strike}`,
  `slip_rate_{dip,strike}`, `traction_{dip,strike}`, `state_variable`), free-surface slice.

### Sign / frame conventions (established, audited 2026-07-07)

| | MFEM | SeisSol | handoff |
|---|---|---|---|
| Frame | UTM 11N, x=E y=N z=Up | same | 1:1 |
| Stress | compression **positive** | compression **negative** | global `×(−1)` at the converter, once |
| Fault tangent frame | `t1 = dip`, `t2 = strike` (FaultBasis/Tandem) | `1 = strike`, `2 = dip` | `tau_strike_MFEM ≡ Ts_SeisSol` exactly (two flips cancel); `tau_dip_MFEM = −Td_SeisSol` |
| `ref_normal` | `[0,−1,0]` (default) | `XRef=(0,−1,0)`, `refPointMethod=1` | identical |
| σ_n | `> 0` compression | `min(σ_n, 0)` | — |
| Pore pressure | `σ_n_eff = σ_n − P_p`, applied unconditionally | — | nc is **already effective** ⇒ **`P_p_pa = 0.0`** (else double-subtract) |

The Andersonian nc's own global attribute states it: `convention = "SeisSol compression-NEGATIVE;
effective (hydrostatic P_p removed)"`.

---

## Gap table

| ID | Gap | Severity | Phase |
|:---|:----------------------------------------------------------------------------|:---------|:------|
| **G1** | No 3-D field ingestion for `rs_a` / `rs_srW`. Resolver seeds `a`,`b` from scalars or a 1-D depth CSV; `V_w` from scalars, box rules, or a boxcar taper. `spatial_friction.cpp:2031`. | **Blocker** | 3 |
| **G2** | `MFEM_VERIFY(f_w_default > 0.0 && f_w_default < 1.0)` at `spatial_friction.cpp:713` rejects `f_w = 0`. | **Blocker** | 2 |
| **G3** | `velocity_safs.h5` was built from the *v3.0.0* CVM nc (297×196×194); the deck ships the statewide nc (452×361×194). `build_pref_velocity_sidecar.py` also hard-codes a path under the old `seisol_quakeworx/` location, **which no longer exists** (relocated 2026-06-26). | Data | 4 |
| **G4** | No MFEM stress sidecar for `safs_stress_andersonian_k1.8.nc`. | Data | 5 |
| **G5** | No on-fault pickpoint / off-fault receiver time-series writer for SAFS. Station writers are hard-coded SCEC grids keyed on five fixed `[problem].tag` values. The deck ships 10 pickpoints + 10 receivers. | Feature | 7 (opt.) |
| **G6** | No energy output (SeisSol `EnergyOutput=1`). | Cosmetic | not planned |
| **G7** | MFEM has **no LTS**. SeisSol's clustered LTS gave a 30.5× speedup over GTS on this mesh. | Accept | 1 (measure) |
| **G8** | MFEM `[mesh].order = 1` (p1) vs SeisSol order 4 (p3). Not the same accuracy per DOF. | Decision | 6 |
| **G9** | `sigma_n_strength_floor_pa = 10e6` in every MFEM SAFS config. The Andersonian field has min σ_n = 2.15 MPa on this fault, so **the floor is active** on shallow facets and makes them stronger than SeisSol. No SeisSol analog. | Decision | 6 |
| **G10** | Compact-bell radius: MFEM measures `r` **in the fault plane** (projection onto the per-DOF dip/strike basis, `spatial_nucleation.hpp:231`); SeisSol's Lua uses the full 3-D Euclidean distance. On a curved fault the two differ by the normal offset; MFEM's `r` ≤ SeisSol's `r`, so MFEM's bump is marginally wider. | Quantify | 6 |
| **G11** | `spatial_dyn_driver.cpp:1418` warns: `rate_state` at np>1 uses **end-of-step ψ** on rank-seam fault QPs (1st order). Fires on every production SAFS run today. | Known | 8 |
| **G12** | *(found during implementation)* `io/data_field_3d.cpp` holds each field's **entire** array on **every rank** (`data_.assign(Nx*Ny*Nz,0)` + full `H5Dread`). v3_4_1's three sidecars total **1329 MB/rank** (v3.0.0: 779 MB) because the statewide CVM cube replaced the smaller one and friction is new. The stock 2N×50r / `--mem=100000M` sizing gives 2000 MB/rank — the sidecars alone are 66% of it. | **Blocker (runtime)** | 8 |
| **G13** | *(found during implementation)* `SEAS_HEADERS` omitted `spatial/code/*.hpp`, and this Makefile has no auto-dependency generation. Adding a member to `RateStateBlock` therefore did **not** rebuild `nucleation_factory.o` / `fault_face_flux.o` on an incremental `make`, silently linking mixed-ABI objects (observed: `seas_test_nucleation_factory` segfaulting). | **Blocker (build)** | 3 |

---

## Constraints

1. **Byte-exact regression contract.** TPV102 / TPV104 / TPV205 / BP5 outputs must not change. Every
   addition below defaults to "absent ⇒ old behaviour". `make test` must stay green.
2. **Do not edit** the files `miniapps/seas/CLAUDE.md` lists as requiring extreme care:

   ```
   friction/dieterich_ruina.hpp      friction/state_evolution.hpp
   domain/elasticity_operator.hpp    fault/fault_basis.hpp
   fault/rate_state_fault.hpp        solver/seas_operator.hpp
   solver/time_stepper.hpp           config/bp5_params.hpp
   ```

   Nothing in this plan needs to.
3. **No NetCDF in C++.** The driver links HDF5 only (`data_field_3d.cpp` uses the raw HDF5 C API). All
   `.nc → .h5` conversion stays in offline Python. Do not add a netcdf-c dependency.
4. **Gmsh v2.2 ASCII only** for `[mesh].path`.
5. **Frontera/Expanse runs need explicit user approval** before submission.
6. Mixed flux ⇒ RK. `mixed_flux != "none"` under ADER aborts by design; do not "fix" that.

---

## Phase 1 — Baseline measurement (no code)

### Goal
Know the actual time-step and per-step cost on the deep19km mesh before committing to a wall-time budget.

### Steps
1. Build the driver locally (`conda activate mfem-dev; make seas_spatial_dyn_driver`).
2. Author a throwaway config that points at the deep19km `.msh` with `--no-sidecar-material`
   (constant material), `interior_flux="scalar"`, `mixed_flux="none"`, `time_integrator="ader"`, and run
   `mpirun -np 8 ./seas_spatial_dyn_driver --config <tmp>.toml --dry-run --print-derived`.
   Record `h_min`, `dt_cfl`.
3. Repeat with `--time-integrator rk4 --mixed-flux adjacent` (still constant material, scalar flux) and
   record `dt_cfl`.

### Expected numbers (derive, then check)
`h` is the **inscribed diameter** `6·Vol/A_total` (`wave_operator.inl:96`) — the same
`2·r_insphere` SeisSol uses. SeisSol reported `dt_min = 5.13722e-5 s` at order 4 (N=3), CFL 0.5, from
`dt = cfl·2r/((2N+1)·c_p)`, hence `min_e(h_e/c_p,e) = 7 · 5.13722e-5 / 0.5 = 7.192e-4 s`.

MFEM: `dt = MixedFluxCflFactor · (cfl · RkCflFactor) · min_e(h_e/c_p,e)` where
`RkCflFactor = 3/(2p+1) = 1.0` at p=1 and `MixedFluxCflFactor(Adjacent, rk_aware) = 0.6`
(`wave_operator.hpp:870`; `spatial_friction.hpp:711`;
`spatial_dyn_driver.cpp:2221`):

```
dt(p=1, rk4, adjacent) = 0.6 · 0.5 · 7.192e-4  =  2.16e-4 s   →  ~695,000 steps × 4 RHS for 150 s
dt(p=1, ader, upwind, cfl_dg_safety=3) = (0.5/9) · 7.192e-4 = 4.00e-5 s  →  ~3.75M steps
dt(p=1, ader, upwind, cfl_dg_safety=1) = (0.5/3) · 7.192e-4 = 1.20e-4 s  →  ~1.25M steps
```

Note `--cfl-dg-safety` is **inert on the RK path** — it only scales `CflSafetyFactor`, which the ADER
branch uses. The mixed-flux+RK production target gets no benefit from Lever A.

For calibration: SeisSol reached only `t ≈ 98 s` of 150 s in 300 min on 160 ranks (job hit the Tapis time
limit), *with* a 30.5× LTS speedup. Expect MFEM GTS at p1 to need a **checkpoint/restart chain**, not a
single 2 h job.

### Acceptance
`dt_cfl` printed and within 2× of the estimate above. If it is wildly smaller, a sliver tet is driving
`h_min` and that must be investigated before anything else (see `project_flattop_badtet_fix` /
`QUALITY_REPORT_deep19km.md`, which reports `eta_min` PASS).

### RESULT (2026-07-09, np=1 `--dry-run --print-derived`, production parity config)

```
[time]    dt_cfl   = 0.000215799 s          predicted 2.16e-4 s   -- MATCHES
[derived] nsteps   = 695,093                predicted ~695,000    -- MATCHES
[derived] mesh h_min = 70.5003 m
[mixed-flux] mode = adjacent, central-flux faces (rank-summed) = 708,810
```

`dt_cfl` lands on the prediction to 3 significant figures. That is a genuine cross-check, not a
coincidence: the estimate was derived from SeisSol's own reported `dt_min = 5.13722e-5 s` at order 4,
inverted to `min_e(h_e/c_p,e) = 7.192e-4 s`. MFEM independently measures the same quantity from the same
mesh through a different code path (`6·Vol/A_total` vs SeisSol's inscribed-sphere diameter) and agrees.
No sliver tet is driving `h_min`.

Setup cost at np=1 was 851 s wall, 15.7 GB peak RSS (serial mesh + all three sidecars). Per-rank at
np=100 the mesh share collapses but the 1329 MB of sidecars does **not** — see G12.

---

## Phase 2 — Allow `f_w = 0` (G2)

### Goal
`RS_muW = 0.0` (zero residual friction — the rate-state analogue of LSW `mu_d = 0`) must parse.

### Why it is safe
`f_w` enters only additively: `f_ss = f_w + (f_LV − f_w)/(1+(V/V_w)⁸)^{1/8}`. There is no division by
`f_w` and no `log(f_w)`. `ψ_ss = a·ln((2V_0/V)·sinh(f_ss/a))` needs `f_ss > 0`; with `f_w = 0`,
`f_ss = f_LV/(1+(V/V_w)⁸)^{1/8}` and `f_LV = max(0, f_0 − (b−a)·ln(V/V_0))` only reaches 0 at
`V = V_0·exp(f_0/(b−a)) = 1e-6·e^{150}`, which is unreachable. `f_ss > 0` strictly ⇒ `ψ_ss` finite.
SeisSol runs this deck with `RS_muW = 0` on v1.1.3. The deck header flags it as an
"energy-max idealization below Lachenbruch & Sass `mu_d ≈ 0.10`".

### Files to change
- `spatial_friction.cpp:713` — relax to `f_w_default >= 0.0 && f_w_default < 1.0`; update the
  R-004 comment to record *why* `0` is admissible and that `V_w > 0` is the load-bearing guard (never
  `V_w = 0` under FL=103).

### Files to create
- `seas_test_slip_law_srw_fw_zero/` (or extend `seas_test_slip_law_srw_psi`): assert
  `PsiSS_SRW(V, V_w, a)` is finite and monotone for `f_w = 0` across `V ∈ [1e-12, 10]` m/s with
  `a=0.015, b=0.019, f_0=0.6, V_0=1e-6, V_w=0.05`; assert the parser accepts `f_w_default = 0.0` and still
  rejects `-0.1` and `1.0`.

### Acceptance
- `make test` green.
- `ParseSpatialFrictionConfigString` with `f_w_default = 0.0` + `state_evolution = "slip_law_srw"` returns
  without abort.

---

## Phase 3 — 3-D friction-field sidecar (G1) — **the only real feature work**

### Goal
Let `[friction.rate_state]` seed per-DOF `a` and `V_w` from a `data_projection_v1` HDF5, exactly as
`[stress]` and `[velocity]` already do. This is what carries the CTM temperature zoning.

### Design
Mirror the stress/velocity pattern. Do **not** invent a new file format; reuse `DataField3D` (trilinear,
`OOBPolicy::Clamp` = ASAGI edge-hold).

**Why baked `a`/`V_w` and not a temperature field + inline `a(T)` formula:** SeisSol itself reads baked
`rs_a`/`rs_srW` through ASAGI. Reading the same two fields means the two codes sample *identical* friction,
and the `a(T)`/`V_w(T)` piecewise logic stays in one place (`build_friction_nc_thermal.py`) instead of
being duplicated in TOML.

#### Python: `safs/project_7.0_preferred/friction/code/build_pref_friction_sidecar.py` (new)

**Decision D5 (user, 2026-07-08): MFEM never reads `.nc`. The friction sidecar is the same file type as the
velocity sidecar** — a `data_projection_v1` HDF5, written by the **same shared helper**
`safs/project_7.0_alternative/velocity/code/sidecar.py::write_sidecar` that
`build_pref_velocity_sidecar.py` already uses. Do **not** hand-roll an HDF5 writer (the stress converter
has its own inline copy; that is technical debt, not a template).

So this script is a sibling of `build_pref_velocity_sidecar.py`, not of `csm_stress_nc_to_mfem_hdf5.py`:

```python
# in : safs_friction_thermal_case2.nc   compound data{rs_a, rs_srW}, dims (z,y,x), f4
# out: .../friction/results/thermal_case2_mfem/friction_safs.h5
sys.path.insert(0, f"{ALT}/velocity/code")
from sidecar import write_sidecar

data = netCDF4.Dataset(NC).variables["data"][:]            # (z, y, x) compound
rs_a   = np.transpose(np.asarray(data["rs_a"],   np.float64), (2, 1, 0))   # -> (x, y, z)
rs_srW = np.transpose(np.asarray(data["rs_srW"], np.float64), (2, 1, 0))

write_sidecar(
    OUT, x, y, z,
    fields       = {"rs_a": rs_a, "rs_srW": rs_srW},
    attrs        = {"schema_version": "data_projection_v1", "crs": "EPSG:32611",
                    "units": "m", "z_positive": "elevation",
                    "source": "safs_friction_thermal_case2.nc (CTM shinevar2024, CASE2 "
                              "a(T)/V_w(T)) -> build_pref_friction_sidecar.py",
                    "mesh_tag": "safv4_deep19km_sub500m_flattop_nwtrim_fixed_opt"},
    field_bounds = {"rs_a":   (0.010, 0.050,  "1"),
                    "rs_srW": (0.05,  1000.0, "m/s")},
)
```

`write_sidecar` already enforces: strictly-monotone axes, no NaN, shape `(Nx,Ny,Nz)`, and every cell inside
`field_bounds` — i.e. it gives the build script's G2 gate (`rs_a > 0`, `rs_srW ∈ [0.05, 1000]`) for free and
raises rather than warns.

Add one guard the helper cannot know about: assert `a(hypocentre) == 0.015` and `V_w(hypocentre) == 0.05`
to 1e-9 by trilinear sampling before writing (the deck's own G4 gate).

Resulting layout, byte-comparable to `velocity_safs.h5`:

```
grid/x (263,) f8   grid/y (194,) f8   grid/z (107,) f8      # ascending, equidistant
fields/rs_a   (263,194,107) f8   attrs units="1",   min_value, max_value
fields/rs_srW (263,194,107) f8   attrs units="m/s", min_value, max_value
```

#### C++: schema

`spatial_friction.hpp`:

```cpp
/// [friction.rate_state.sidecar] — per-DOF a / V_w (and optionally b, Dc) from a
/// data_projection_v1 HDF5.  Absent ⇒ enabled=false ⇒ resolver behaviour byte-identical.
struct FrictionSidecarSpec
{
   bool        enabled          = false;
   std::string path;                       // required when enabled
   std::string a_field          = "rs_a";
   std::string V_w_field        = "rs_srW";
   std::string b_field;                    // "" ⇒ keep b_default
   std::string Dc_field;                   // "" ⇒ keep Dc_default
   bool        far_field_clamp  = true;    // OOBPolicy::Clamp (ASAGI edge hold)
};
```

Add `FrictionSidecarSpec sidecar;` to `RateStateBlock`.

Loaded fields (owned by the driver, handed to the resolver):

```cpp
struct RateStateSidecarFields
{
   std::unique_ptr<DataField3D> a, V_w, b, Dc;   // b/Dc may be null
};
```

`SpatialFrictionResolver` gains an optional ctor argument
`std::shared_ptr<const RateStateSidecarFields>`; the default ctor leaves it null, so every existing call
site (and every TPV/BP5 test) is byte-unchanged.

#### C++: parser (`spatial_friction.cpp`, in `parse_rate_state`)

- Read `[friction.rate_state.sidecar]`; `path` required and non-empty when the table is present.
- **`MFEM_VERIFY(!(sidecar.enabled && depth_profile.enabled))`** — both seed `a`/`b`; a config that sets
  both is ambiguous and must abort with a message naming both keys.
- Keep the existing `a_default < b_default` guard satisfied by the config (`0.015 < 0.019`); `a_default`
  becomes a pure fallback.

#### C++: resolver (`resolve_rs_impl`, around `spatial_friction.cpp:2031`)

Insert **after** the depth-profile/scalar seed and **before** the spatial-rule loop, so `box` /
`boxcar_taper` rules can still override a sidecar value:

```cpp
if (sidecar_)                       // null ⇒ untouched legacy path
{
   a_i   = sidecar_->a  ->Evaluate(x, y, z);
   V_w_i = sidecar_->V_w->Evaluate(x, y, z);
   if (sidecar_->b)  { b_i  = sidecar_->b ->Evaluate(x, y, z); }
   if (sidecar_->Dc) { Dc_i = sidecar_->Dc->Evaluate(x, y, z); }
}
```

The existing `V_w_i > 0` guard (`:2166`) and the per-DOF `a > b` allowance (R-011, `:2153`) already cover
the VS zone (`a` up to 0.0332 > `b` = 0.019).

#### C++: driver (`spatial_dyn_driver.cpp`)

- Construct `RateStateSidecarFields` right before `ResolveRateState` (after the mesh + material exist).
- `OOBPolicy` = `Clamp` when `far_field_clamp`, else `Abort`.
- Extend `--print-derived` to report, over fault DOFs: `min/max/mean` of `a` and `V_w`, the count of
  DOFs with `a − b < 0` (the VW area), and the values at `[hypocenter]`. This is the offline↔online
  cross-check against the SeisSol build script's own G4 gate.

### Files to create
- `safs/project_7.0_preferred/friction/code/build_pref_friction_sidecar.py`
- `seas_test_spatial_friction_sidecar/` — three tests:
  1. **Round-trip**: write a synthetic 3×3×3 sidecar with a known linear `a(x,y,z)`; resolve at 8 corners
     and the centre; assert trilinear values to 1e-12.
  2. **Clamp**: query outside the hull with `far_field_clamp = true` ⇒ edge value; with `false` ⇒ abort.
  3. **Rule precedence**: sidecar + a `box` rule setting `V_w` ⇒ rule wins inside the box, sidecar outside.
- Parser tests: `sidecar` + `depth_profile` together ⇒ abort; missing `path` ⇒ abort.

### Files to change
- `spatial_friction.hpp` (schema + resolver ctor)
- `spatial_friction.cpp` (parse + resolve)
- `spatial_dyn_driver.cpp` (load + `--print-derived`)
- `safs/project_7.0_alternative/document/spatial_friction_config_schema.md` (document the new block; also
  fix the stale "rate-state ignores `[nucleation]`" claim at lines 256-258 / 376-377)
- `miniapps/seas/CLAUDE.md` (fix the stale "single nucleation kind" paragraph)

### Acceptance
- `make test` green; `seas_test_spatial_friction_resolver` byte-unchanged (no sidecar ⇒ no behaviour change).
- `--print-derived` on the real sidecar reports `a(hypo) = 0.015`, `V_w(hypo) = 0.05`, `min a = 0.015`,
  `max a ≈ 0.0332`, `min V_w = 0.05`, `max V_w = 1000`.
- Fault-DOF VW fraction is consistent with "350 °C isotherm at ~13.2 km at the hypocentre column".

### RESULT (2026-07-09, np=1 `--dry-run --print-derived`, production parity config)

```
[print-derived] friction sidecar (a, V_w)
  fault DOFs (global) : 517779
  a      min/mean/max : 0.015 / 0.0177678 / 0.0298597
  V_w    min/mean/max : 0.05 / 225.437 / 1000 m/s
  velocity-weakening  : 371373 / 517779 DOFs with a - b < 0  (71.72 %)
  sidecar at [nucleation.gradual_overstress_compact_circular] centre (609063, 3.70953e+06, -10000):
      a   = 0.015          <- the CASE2 VW-plateau anchor
      V_w = 0.05 m/s       <- ditto
```

All gates met, with one refinement to the wording above: `max a = 0.0299` **on the fault**, not 0.0332.
0.0332 is the maximum over the whole CTM *grid*; the fault bottom (−19.3 km) never reaches the grid's
hottest cells. The two numbers are consistent, and the offline builder's own G2 gate still checks the
grid-wide `0.0332263`.

The VW fraction (71.7 %) matches the geometry: the 350 °C isotherm sits at ≈13.2 km at the hypocentre
column and the fault runs to −19.33 km, so ≈13.2/19.33 ≈ 68 % of a uniformly-sampled fault should be VW;
71.7 % reflects the shallower isotherm away from the hypocentre column.

`seas_test_spatial_friction_resolver` stayed at 127/127 (byte-unchanged), and 24 affected test binaries
pass from a clean build.

---

## Phase 4 — Rebuild the CVM material sidecar (G3)

### Goal
`velocity_safs.h5` sampled from the **same** `safs_material_cvm.nc` the SeisSol job read.

### Files to change
- `safs/project_7.0_preferred/velocity/code/build_pref_velocity_sidecar.py`
  - Replace the hard-coded source path

    ```python
    NC = f"{ALT}/seisol_quakeworx/safs_seisol_v3_0_0_LSW_PREFERRED/safs_material_cvm.nc"
    ```

    (**this path no longer exists** — `seisol_quakeworx` moved to `safs/seisol_quakeworx/` on 2026-06-26)
    with `--in` / `--out` argparse options. Keep the old defaults out; make both required.

### Output
`project_7.0_preferred/velocity/results/` `multiscale_statewise_cvm_v3_4_1/velocity_safs.h5`
(new directory so the v3.0.0 sidecar is not clobbered; both are gitignored, scp to Expanse).

Source nc: `safs_material_cvm.nc` from the deck directory (452 × 361 × 194; the existing h5 is
297 × 196 × 194 — a **different, smaller** grid).

### Acceptance
- `Vp = sqrt((λ+2μ)/ρ)`, `Vs = sqrt(μ/ρ)`, `density = ρ`, all finite, `ρ > 0`, `μ ≥ 0`.
- Grid attrs: `crs="EPSG:32611"`, `z_positive="elevation"`, axes ascending.
- Mesh bbox `x[253447,711231] y[3532686,3956757] z[−39329,+26]` pokes outside the CVM hull only at the
  far-field absorbing corners (`x` high ≤ 6.2 km, `y` low ≤ 10.3 km) ⇒ **`far_field_clamp = true` is
  mandatory**, exactly as SeisSol relies on ASAGI's edge clamp + the yaml `!ConstantMap`.

---

## Phase 5 — Convert the Andersonian stress field (G4)

### Goal
`safs_stress_andersonian_k1.8.nc` → schema-v1 HDF5, compression-**positive**, effective Pa.

### Command
```
python3 safs/seisol_quakeworx/toolbox/on_fault_stress_projection_csm/csm_stress_nc_to_mfem_hdf5.py \
  --in  ~/Downloads/seisol_quakeworx/safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2/safs_stress_andersonian_k1.8.nc \
  --out safs/project_7.0_preferred/stress/results/andersonian_k1p8_mfem/safs_stress_andersonian_k1.8_mfem.h5
```

The script needs no modification: the member names are already `s_xx..s_xz`, and its `min<max` nudge
already covers `s_yz ≡ s_xz ≡ 0`.

### Cautions
- The script's `HYPO_DEFAULT = (606971.0, 3707270.0, -4965.62)` is the **v3.0.0** hypocentre used only for
  a diagnostic print. Pass the v3.4.1 on-fault hypocentre `(609062.8722, 3709528.1324, -10000.0)` if the
  script exposes a flag; otherwise ignore the print (it is not used in the output).
- **`[pore_pressure] P_p_pa = 0.0`.** The nc is already effective (its own `convention` attribute says so).
  A non-zero `P_p_pa` would double-subtract in `ProjectFaultPreStress`.
- Grid `z[−20000,+3000]` and `x[195,628]k / y[3685,3965]k` contain the fault bbox
  `x[363.5,626.7]k y[3692.0,3839.0]k z[−19329,+26]`, so the default `OOBPolicy::Abort` in
  `ApplyCsmStressSidecar` is fine — stress is only projected onto **fault** DOFs.

### Acceptance
- Sampled at the hypocentre facet: `σ_n = 169.83 MPa`, `mu_app ≤ 0.2931`, `τ_0 ≥ 0.006 MPa` — the numbers
  the deck header records for this fault.
- `min σ_n` over fault DOFs `≈ 2.15 MPa`, 0 non-finite.

---

## Phase 6 — Production TOML

### File to create

```
safs/project_7.0_preferred/config/
    spatial_friction_safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2.toml
```

```
[meta]
schema_version = 1
law            = "rate_state"
description    = "MFEM port of seisol_quakeworx/safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2 (FL=103, CTM temperature-zoned rs_a/rs_srW via sidecar, f_w=0, Andersonian k=1.8 effective stress, statewide CVM) with MIXED FLUX (adjacent) + RK. 'Do what SeisSol did, but with mixed flux.'"

[material_constant_fallback]              # only reachable via --no-sidecar-material
lambda = 32.0e9
mu     = 32.0e9
rho    = 2670.0

[pore_pressure]
P_p_pa            = 0.0                   # nc is EFFECTIVE (hydrostatic Pp already removed)
P_p_grad_pa_per_m = 0.0
min_sigma_n_pa    = 1.0e6                 # inert: min sigma_n on this fault is 2.15 MPa

[mesh]
path  = "safs/project_7.0_preferred/meshing_deep19km/results/safv4_deep19km_sub500m_flattop_nwtrim_fixed_opt_safstags.msh"
order = 1                                 # DECISION D2 (SeisSol runs p3)

[velocity]
use_sidecar     = true
model           = "multiscale_statewise"
dataset_root    = "safs/project_7.0_preferred"
override_path   = "safs/project_7.0_preferred/velocity/results/multiscale_statewise_cvm_v3_4_1/velocity_safs.h5"
far_field_clamp = true                    # MANDATORY: mesh far field exceeds CVM hull

[material]
kind         = "sidecar_hdf5"             # required: interior_flux="matrix" forbids kind="constant"
sidecar_path = "safs/project_7.0_preferred/velocity/results/multiscale_statewise_cvm_v3_4_1/velocity_safs.h5"

[stress]
kind         = "sidecar_hdf5"
sidecar_path = "safs/project_7.0_preferred/stress/results/andersonian_k1p8_mfem/safs_stress_andersonian_k1.8_mfem.h5"

[numerics]
ader_order      = 2                       # ignored on the RK path
mixed_flux      = "adjacent"
interior_flux   = "matrix"                # heterogeneous CVM ⇒ bimaterial Riemann
time_integrator = "rk4"                   # central flux is UNSTABLE under ADER
cfl             = 0.5                     # = SeisSol CFL
use_pml         = false

[time]
tfinal     = "150s"                       # = SeisSol EndTime
t_initial  = 0.0
dt_initial = "auto"
dt_max     = "0.1s"                       # = SeisSol FixTimeStep

[output]
output_dir               = "output_safs_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2"
restart_prefix           = "cp"
paraview_volume          = "off"
paraview_bulk            = "off"
paraview_fault           = "hdf5"
paraview_fault_dt        = "4.0s"         # = SeisSol Elementwise printtimeinterval_sec
paraview_fault_zfp_tol   = 1.0e-12
paraview_free_surface    = "vtu"
paraview_free_surface_dt = "4.0s"         # = SeisSol SurfaceOutputInterval
max_snapshots            = 5000
checkpoint_every_steps   = 10000

[nucleation]
kind = "gradual_overstress_compact_circular"

[nucleation.gradual_overstress_compact_circular]
# = SeisSol Tnuc_s = 75e6 * exp(r^2/(r^2-R^2)), R=2000, smoothStep over [s_0, s_0+t_0] = [0,1] s.
# Centre = the on-fault point at EXACTLY 10 km depth closest to the deck-lineage geographic target.
center_x_m   =  609062.8722
center_y_m   = 3709528.1324
center_z_m   =  -10000.0
radius_m     =    2000.0
delta_tau_pa =      75.0e6
T_nuc_s      = "1.0s"

[friction]
sigma_n_strength_floor_pa = -1.0          # DECISION D1: -1 = disabled = SeisSol parity

[friction.rate_state]
state_evolution = "slip_law_strong_rate_weakening"   # = SeisSol FL=103
f_w_default     = 0.0        # = RS_muW   (requires Phase 2)
V_w_default     = 0.05       # FALLBACK ONLY (sidecar governs)
a_default       = 0.015      # FALLBACK ONLY (sidecar governs); must stay < b_default
b_default       = 0.019      # = RS_b  (constant; v1.1.3 reads no spatial rs_b)
Dc_default      = 0.10       # = rs_sl0 (constant)
V_0_default     = 1.0e-6     # = RS_sr0
f_0_default     = 0.6        # = RS_f0
V_init_default  = 1.0e-12    # = RS_iniSlipRate1
sigma_n_default = 169.83e6   # FALLBACK ONLY (live sigma_n is per-DOF from the stress sidecar)
eta             = "auto"     # 0.5*sqrt(mu*rho) per DOF  == SeisSol radiation damping

[friction.rate_state.sidecar]             # NEW (Phase 3) — same file type as [velocity]
path            = "safs/project_7.0_preferred/friction/results/thermal_case2_mfem/friction_safs.h5"
a_field         = "rs_a"
V_w_field       = "rs_srW"
far_field_clamp = true                    # CTM hull covers the fault, not the far field

# No [boundary] block: the driver falls back to fault=101, natural={102}, absorbing={103,104},
#   which is exactly what the safstags mesh uses.
# No [fault_geometry] block: defaults ref_normal=[0,-1,0], up=[0,0,1] == XRef/refPointMethod=1.
# No [friction.rate_state.depth_profile]: mutually exclusive with [.sidecar].
```

### Decisions — RESOLVED (user, 2026-07-08)

**D1 — `sigma_n_strength_floor_pa`: two arms.** Every existing MFEM SAFS RS config sets `10.0e6`. On this
fault `min σ_n = 2.15 MPa`, so the floor **is active** on shallow facets and makes them stronger than
SeisSol's. Author **both** configs:

| config | `sigma_n_strength_floor_pa` | role |
|---|---|---|
| `..._THERMAL_CASE2.toml` | `-1.0` (disabled) | **parity arm — run this first** |
| `..._THERMAL_CASE2_snfloor10.toml` | `10.0e6` | fallback arm, only if the parity arm shows trace-cell runaway |

The two files differ in **that one key and `output_dir`** and nothing else. Runaway signature (from
`project_pre_k25_case2_runaway_slip_not_mesh`): isolated shallow cells at σ_n 0.1–5 MPa sliding ~4 m/s
indefinitely, visible in the first few fault snapshots. `min_sigma_n_pa = 1.0e6` is a *different* knob
(it clamps σ_n_eff, not the strength) and is inert here — keep it in both.

**D2 — `[mesh].order`: p1 production + p2 convergence window.** `order = 1` for the 150 s run (matches every
existing SAFS production run and the Phase-1 cost estimate). Separately run `order = 2` with
`tfinal = "20s"` as a convergence check: if the p2 rupture front and peak slip rate agree with p1 to within
the snapshot cadence, p1 resolves the cohesive zone well enough; if not, the p1 result is quantitatively
suspect and must be reported as such. Under mixed-flux+RK, p2's `dt` penalty is only
`RkCflFactor(2)/RkCflFactor(1) = 0.6`, but DOFs rise ~2.5× ⇒ ~4× total cost. Add this as gate 6b in Phase 8.

**D3 — `time_integrator`: `rk4`.** Fewer RHS evaluations per step than `rk45` at the same fixed `dt`. The
existing sbatch scripts pass `--time-integrator rk45`; the sbatch override wins, so set the default in the
sbatch to `rk4` for this case.

**D4 — `tfinal`: 150 s** (= `EndTime`). But the reference SeisSol run only reached **t ≈ 98 s** before its
Tapis time limit, so **0–98 s is the entire comparison window** the reference data supports.

**D5 — friction sidecar file type: `data_projection_v1` HDF5**, identical in kind to the velocity sidecar,
written by the shared `sidecar.py::write_sidecar`. MFEM reads no NetCDF anywhere. See Phase 3.

**D6 — receivers/pickpoints (Phase 7): deferred** to a follow-up. Phases 1–6 + 8 give a runnable,
comparable simulation on fault snapshots and free-surface slices; the per-point `.dat` writers upgrade that
to waveform misfit once there is a run worth comparing.

---

## Phase 7 — SAFS receiver / pickpoint output (G5) — **optional, scope decision**

Without this, the only comparison surface is ParaView fault snapshots at 4 s cadence and free-surface
slices. The deck's `safs-faultreceiver-*.dat` (10 on-fault, `printtimeinterval = 10` steps) and
`safs-receiver-*.dat` (10 off-fault, `pickdt = 0.005` s = 100 Hz) are the highest-signal comparison data
and they already exist for the SeisSol side.

### Design
Generalize the existing `wire_stations` block (`spatial_dyn_driver.cpp:3297`) rather than
adding a fifth hard-coded SCEC grid:

- New TOML: `[stations] on_fault_file = "..."`, `off_fault_file = "..."`, `on_fault_dt`, `off_fault_dt`.
- Reuse the `.dat` format the deck ships (whitespace-separated `x y z`, one point per line — see
  `safs_pickpoints.dat` / `safs_receivers.dat`, both 10 lines).
- On-fault: nearest fault DOF per point (assert the residual is `< 0.5 · median fault edge = 200 m`);
  write `t, slip_strike, slip_dip, slip_rate_strike, slip_rate_dip, traction_strike, traction_dip, sigma_n, psi`.
- Off-fault: MFEM point-location on the volume mesh; write `t, vx, vy, vz`. Receivers sit at `z = −1 m`
  precisely because SeisSol v1.1.3 drops `z = 0` points; MFEM has no such restriction, but keep `z = −1`
  so both codes sample the same location.

### Recommendation
**Defer to a follow-up.** Phases 1–6 + 8 produce a runnable, comparable simulation. Phase 7 upgrades the
comparison from qualitative (rupture front, slip map) to quantitative (waveform misfit). Do it once the
first run has produced something worth comparing.

---

## Phase 8 — Job scripts and validation

### Files to create
```
jobs/safs/safs_expanse/
  spatial_dyn_ratestate_v3_4_1_pref_thermal_case2_cvm_mixedflux_rk_deep19km_expanse_2N100r_compute_2hr_safs.sbatch
  ..._2hr_restart_safs.sbatch          # the --restart sibling
```

Copy the existing `..._v3pref_vw002_cvm_mixedflux_rk_safv4deep_..._2hr_safs.sbatch`
verbatim and change only `SAFS_CONFIG`, `SAFS_MESH`, `SAFS_OUT`, and the job name. Keep the existing
`mixed_flux != none && integrator == ader` abort guard.

**Promote the validated perf levers**: `--deriv-cache` is `≤1e-12`-parity and works on the RK path
(`--face-cache` and `--shared-ck-recursion` are ADER-only; `--cfl-dg-safety` is inert on RK). Today's
production scripts pass none of them.

Sidecars to `scp` to Expanse alongside the mesh (all gitignored, ~700 MB total):
`velocity_safs.h5` (from Phase 4), `safs_stress_andersonian_k1.8_mfem.h5` (Phase 5),
`safs_friction_thermal_case2_mfem.h5` (Phase 3).

### Validation gates, in order

1. **`--dry-run --print-derived`, np=1, local.** Confirms mesh reads, all three sidecars load, boundary
   attrs resolve to 101/102/103/104, `dt_cfl` matches Phase 1.
2. **`--print-derived` friction gate.** `a(hypo)=0.015`, `V_w(hypo)=0.05`, `min/max a = 0.015/0.0332`,
   `min/max V_w = 0.05/1000`. Compare against `build_friction_nc_thermal.py`'s G2/G4 gates.
3. **Stress gate.** `σ_n(hypo facet) = 169.83 MPa`; `min σ_n` over fault DOFs ≈ 2.15 MPa; `max mu_app`
   ≈ 0.2931 `< f_0` ⇒ **no t=0 pre-slip**.
4. **Nucleation gate. — FAILED at 75 MPa (2026-07-09). DECISION REQUIRED.**

   The deck's own margins are thin: barrier `S_E = a·σ_n·ln(V_dyn/V_init) = 0.015 · 169.83e6 · 27.631
   = 70.39 MPa` (75/70.39 = ×1.065) and static excess `(f_0 − mu_app)·σ_n = 72.11 MPa` (×1.040, binding).

   MFEM's `--print-derived` applies a **stricter** criterion — it compares `|τ_pre + Δτ|` against the
   *steady-state* strength `f_ss·σ_n` rather than against `f_0·σ_n` — and on the production config it
   reports:

   ```
   [derived] nucleation peak |F(r) * delta_tau|       = 7.49448e+07 Pa
   [derived] velocity-weakening DOFs inside patch      = 851
   [derived] most-overstressed VW DOF at (608882, 3.70954e+06, -9810.02)
   [derived] amplitude |delta_tau| at that DOF         = 7.36938e+07 Pa
   [derived] nucleation overshoot                      = -6.85895 MPa   (BELOW steady-state strength)
   [derived] WARNING: ... the patch only creeps faster and may not nucleate within tfinal.
                      Increase delta_tau_*_pa.
   ```

   Scaling the amplitude (only `Δτ` moves; `τ_pre` and `f_ss·σ_n` are fixed, and the bell factor at that
   DOF is `F = 0.9826`):

   | `delta_tau_pa` | overshoot | verdict |
   |---|---|---|
   | 75 MPa (deck value) | **−6.86 MPa** | will not nucleate |
   | 80 MPa | −1.95 MPa | still below strength |
   | **85 MPa** (the plan's escalation) | **+2.97 MPa** | nucleates |
   | 90 MPa | +7.88 MPa | nucleates |

   **This is a physics decision, not a bug.** Two options, and they answer different questions:
   - **Keep 75 MPa** — a *faithful* port. SeisSol's own criterion says it nucleates (×1.040 static
     excess); MFEM's stricter `f_ss` criterion says it does not. Running it settles which criterion is
     right for this fault, at the cost of possibly burning the job.
   - **Raise to 85 MPa** — the plan's pre-authorized escalation. Guarantees nucleation on both criteria,
     but the run is then no longer the *same* forcing SeisSol used, so rupture-front timing is not
     directly comparable.

   Do **not** touch `radius_m` or the hypocentre either way.

   Also quantify G10 here: report `max |r_MFEM − r_3D|` over bump-support DOFs; if it exceeds ~1% of `R`,
   the effective bump is measurably wider than SeisSol's.
5. **np=4 smoke, 200 steps, local.** Watch the R-020 `rate_state at np>1` warning (G11) — it fires on
   every SAFS RS run today; confirm it is the ADER-substep path and not the RK path, or accept it.
6. **Short Expanse pilot (user approval required).** `--tfinal 5s`, 2N×100r, `--paraview-fault-dt 0.5s`.
   Confirm: rupture nucleates in the bump, no NaN, no free-surface daylighting blow-up, `V_max` peaks and
   decays. Wall-time per simulated second × 150 gives the real budget.
6b. **p2 convergence window (D2).** Same config with `[mesh].order = 2`, `tfinal = "20s"`. Compare the
   rupture-front position and peak slip rate against the p1 run at `t ∈ {4, 8, 12, 16, 20}` s. Agreement
   within the snapshot cadence ⇒ p1 resolves the cohesive zone; disagreement ⇒ report the p1 150 s result
   as qualitative only. ~4× the p1 cost per simulated second.
7. **Production, checkpoint-chained.** Run the **parity arm** (`sigma_n_strength_floor_pa = -1.0`) first.
   Inspect the first two fault snapshots for the trace-cell runaway signature before committing the full
   wall-time; switch to the `snfloor10` arm only if it appears. Compare against
   `~/Downloads/SAFS_V3.4.1_..._outputs/`:
   fault-surface `SRs/SRd/Sls/Sld/Ts0/Td0/Pn0` at matched `t ∈ {4, 8, ..., 96}` s. Reference data stops at
   t ≈ 98 s. Key metrics: rupture-front arrival along the corridor arclength `s`, peak slip rate,
   whether the front crosses the San Gorgonio gate (PREFERRED corridor `kbar_min = 0.98–0.99`, barely above
   the 0.9 floor — the run decides), final slip, `Mw`.

Use `safs/seisol_quakeworx/toolbox/combined_workflow/postprocess.ipynb` for the SeisSol side — it already
reads any SeisSol output directory and emits `Mw` / slip / stress-drop / finished?.

---

## Risks

| Risk | Mitigation |
|:----------------------------------------|:----------------------------------------|
| **Nucleation margins are ×1.04 (binding).** MFEM's p1 discretization resolves the 583 m supra-`S_E` core with ~1.5 median fault edges; the bump may fail to nucleate where SeisSol's p3 succeeded. | Gate 4. Escalate `delta_tau_pa → 85e6` per the deck's own guidance. |
| **Trace-cell runaway with `sigma_n_strength_floor_pa = -1`** (D1). | Run the `10e6` arm in parallel; `project_pre_k25_case2_runaway_slip_not_mesh` documents the signature (isolated cells sliding ~4 m/s forever at σ_n 0.1–5 MPa). |
| **GTS cost.** MFEM has no LTS; SeisSol's 30.5× LTS speedup is unavailable. | Phase 1 measures it. Checkpoint-chain the production run. The mesh is *smaller* than what MFEM already runs at 2N×100r, so this is a budget problem, not a feasibility one. |
| **Mixed flux is not SeisSol's flux.** Central-on-fault-adjacent faces vs pure Godunov. Differences in the rupture front are the *point* of the experiment, not a bug. | Also run the `--mixed-flux none --time-integrator ader` arm (the upwind sibling sbatch) so the flux is the only variable. |
| **p1 vs p3.** Under-resolved cohesive zone ⇒ artificially high peak slip rate / delayed front. | D2; p2 convergence check on a 20 s window. |
| **Doc rot.** `CLAUDE.md` and `spatial_friction_config_schema.md` both make stale claims about spatial-driver nucleation. An implementer following them will do the wrong thing. | Phase 3 fixes both. |
| **Sidecar memory (G12).** The three sidecars are replicated per rank and total 1329 MB. At 2 GB/rank the job likely OOMs after ~20 min of setup. | The new sbatch defaults to `--mem=200000M` (4 GB/rank; `compute` nodes are exclusive and charged whole-node, so this is free) and carries a pre-flight that recomputes the ratio from `SLURM_MEM_PER_NODE / SLURM_NTASKS_PER_NODE` and errors above 90%. If memory ever gets tight, crop the CVM sidecar to the mesh bbox — that alone removes ~half of `velocity_safs.h5`. |
| **Stale-ABI builds (G13).** Any future change to a struct in `spatial_friction.hpp` silently corrupts incrementally-built binaries. | Fixed: `SPATIAL_HEADERS` added to `SEAS_HEADERS`, plus an explicit prerequisite on the `fault_face_flux.o` rule (which uses `DYNAMIC_HEADERS`). Verified by touching the header and confirming all seven dependent TUs recompile. |

---

## Summary of files

**Create (6):**

```
safs/project_7.0_preferred/friction/code/build_pref_friction_sidecar.py
    uses the shared velocity/code/sidecar.py::write_sidecar; no hand-rolled HDF5

safs/project_7.0_preferred/config/
    spatial_friction_safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2.toml
        parity arm, sigma_n_strength_floor_pa = -1.0
    spatial_friction_safs_seisol_v3_4_1_RSSRW_PREFERRED_THERMAL_CASE2_snfloor10.toml
        fallback arm; differs only in that key + output_dir

seas_test_spatial_friction_sidecar/          3 tests
                                             + the f_w = 0 parser/psi tests from Phase 2

jobs/safs/safs_expanse/
    spatial_dyn_ratestate_v3_4_1_pref_thermal_case2_..._2hr_safs.sbatch  (+ restart sibling)
```

**Change (6):**

```
spatial/code/spatial_friction.hpp     FrictionSidecarSpec, RateStateSidecarFields,
                                      resolver ctor
spatial/code/spatial_friction.cpp     :713 f_w >= 0 guard; parse + resolve the sidecar
drivers/spatial_dyn_driver.cpp        load sidecar fields, extend --print-derived

safs/project_7.0_preferred/velocity/code/build_pref_velocity_sidecar.py
                                      argparse; stale hard-coded path
safs/project_7.0_alternative/document/spatial_friction_config_schema.md
                                      new block; fix stale nucleation claims
miniapps/seas/CLAUDE.md               fix stale "single nucleation kind" paragraph
```

**Re-run, no change (2):**
- `csm_stress_nc_to_mfem_hdf5.py` → the Andersonian k=1.8 sidecar
- `msh_to_puml.py` — not needed; the mesh already exists on both sides

**Zero work:** mesh, boundary attrs, flux/integrator wiring, FL=103 friction law, compact-bell nucleation,
stress reader, material reader, checkpoint/restart, ParaView fault + free-surface output.
