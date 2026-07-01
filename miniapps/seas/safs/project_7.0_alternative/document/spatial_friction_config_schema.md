# Spatial-friction config TOML schema (Phase 0, rev-3)

**Status:** schema spec only.  Phase 1 implements the parser; Phase 3b
adds the constant-tensor stress source; Phase 4 wires the schema into
the `seas_spatial_dyn_driver` CLI.

This document is the single source of truth for the TOML grammar that
`seas_spatial_dyn_driver --config <PATH.toml>` consumes.  Phase 1's
`LoadSpatialFrictionConfig(...)` parser MFEM_VERIFY's every validation
rule below.  When the schema is bumped (e.g. `schema_version = 2`), this
document must be edited in the same commit.

---

## Top-level layout

```toml
[meta]                         # mandatory
[material_constant_fallback]   # mandatory
[pore_pressure]                # mandatory
[mesh]                         # mandatory
[velocity]                     # mandatory
[stress]                       # mandatory
[numerics]                     # mandatory
[time]                         # mandatory
[output]                       # mandatory
[nucleation]                   # optional (when absent, no nucleation perturbation)
[nucleation.gradual_overstress]       # required when [nucleation].kind = "gradual_overstress"
[friction]                     # optional table; holds the law-agnostic sigma_n strength floor
[friction.slip_weakening]      # required when [meta].law = "slip_weakening"
[[friction.slip_weakening.spatial]]   # zero or more
[friction.rate_state]          # required when [meta].law = "rate_state"
[[friction.rate_state.spatial]]       # zero or more
```

Either `[friction.slip_weakening]` or `[friction.rate_state]` is
populated, never both.  The non-matching block MUST be absent; the
parser aborts on extra keys.

---

## `[meta]`

| Key              | Type    | Default       | Validation              |
|------------------|---------|---------------|-------------------------|
| `schema_version` | int     | —             | must equal `1`          |
| `law`            | string  | —             | `"slip_weakening"` or `"rate_state"` |
| `description`    | string  | `""`          | free-form               |

---

## `[material_constant_fallback]`

Used when `--no-sidecar-material` is on or the velocity sidecar load
fails.  D-3 does NOT relax these checks:

| Key      | Unit  | Default  | Validation              |
|----------|-------|----------|-------------------------|
| `lambda` | Pa    | `32.0e9` | `lambda + 2*mu > 0`     |
| `mu`     | Pa    | `32.0e9` | `mu > 0`                |
| `rho`    | kg/m³ | `2670.0` | `rho > 0`               |

Equal `lambda = mu` produces Poisson's ratio `ν = 0.25` (typical crust).

---

## `[pore_pressure]`

Effective normal stress is `σ_n_eff = σ_n_total − P_p` with the depth
convention `z < 0` is below the free surface (CLAUDE.md "Depth
coordinate").  Pore pressure is

```
P_p(depth) = P_p_pa + P_p_grad_pa_per_m * max(0, -z_dof)
```

| Key                 | Unit  | Default | Validation                                |
|---------------------|-------|---------|-------------------------------------------|
| `P_p_pa`            | Pa    | `0.0`   | none                                      |
| `P_p_grad_pa_per_m` | Pa/m  | `0.0`   | warn-only when negative                   |
| `min_sigma_n_pa`    | Pa    | `0.0`   | `>= 0`; `0` means "no floor"; otherwise σ_n_eff is clamped from below |

Hydrostatic gradient ≈ `1.0e4` Pa/m.

---

## `[mesh]`

| Key    | Type   | Default | Validation                                |
|--------|--------|---------|-------------------------------------------|
| `path` | string | —       | non-empty; file existence is checked by the mesh loader, not the parser |
| `order`| int    | `1`     | `>= 1`                                    |

---

## `[velocity]`

| Key             | Type   | Default       | Validation                                       |
|-----------------|--------|---------------|--------------------------------------------------|
| `model`         | string | `"cvmh"`      | `"cvmh"`, `"cvm_s4.26.m01"`, `"multiscale_statewise"` |
| `dataset_root`  | string | `""`          | non-empty when `use_sidecar` is true AND `override_path` is empty |
| `override_path` | string | `""`          | when non-empty, bypasses model-based resolution  |
| `use_sidecar`   | bool   | `true`        | when `false`, the driver skips the sidecar load and uses `[material_constant_fallback]` — `model` / `dataset_root` / `override_path` are then ignored and may be empty |

Path resolution (rev-3, R-003 corrected; applies only when `use_sidecar = true`):

```
cvmh                 → <dataset_root>/velocity/results/cvmh/velocity_safs.h5
cvm_s4.26.m01        → <dataset_root>/velocity/results/cvm_s4.26.m01/velocity_safs.h5
multiscale_statewise → <dataset_root>/velocity/results/multiscale_statewise_cvm/velocity_safs.h5
```

There is **no `mesh_tag`** in the filename — one sidecar per CVM model.

**CLI interaction.** The CLI flag `--no-sidecar-material` is honoured as an
override: if either `use_sidecar = false` (TOML) or `--no-sidecar-material`
(CLI) is set, the constant-fallback path is taken.  The rank-0 startup log
prints `gated by TOML` / `CLI` / `TOML+CLI` to indicate which source
forced the constant path.  This combination is intentional: TOMLs that
require Phase H heterogeneous material (`use_sidecar = true`) can be
temporarily forced through the scalar path from the CLI without editing
the file.

---

## `[stress]` (D-1)

Two modes; exactly one is populated.

| Key                          | Mode             | Unit | Default | Validation                       |
|------------------------------|------------------|------|---------|----------------------------------|
| `kind`                       | always           | str  | —       | `"constant_tensor"` or `"sidecar_hdf5"` |
| `sigma_xx_pa`                | constant_tensor  | Pa   | `0.0`   | required; absent in sidecar mode |
| `sigma_yy_pa`                | constant_tensor  | Pa   | `0.0`   | required; absent in sidecar mode |
| `sigma_zz_pa`                | constant_tensor  | Pa   | `0.0`   | required; absent in sidecar mode |
| `sigma_xy_pa`                | constant_tensor  | Pa   | `0.0`   | required; absent in sidecar mode |
| `sigma_yz_pa`                | constant_tensor  | Pa   | `0.0`   | required; absent in sidecar mode |
| `sigma_xz_pa`                | constant_tensor  | Pa   | `0.0`   | required; absent in sidecar mode |
| `sidecar_path`               | sidecar_hdf5     | str  | `""`    | required; absent in constant mode |

Cauchy tensor is in EAST-NORTH-UP frame, **compression positive** (SEAS
convention, matches `stress/code/hickman_and_zoback_*`).

---

## `[numerics]`

| Key          | Type | Default | Validation                                   |
|--------------|------|---------|----------------------------------------------|
| `ader_order` | int  | `2`     | `>= 1`                                       |
| `mixed_flux` | str  | `"none"`| `"none"`, `"adjacent"`, `"all_continuous"`   |
| `cfl`        | real | `0.5`   | `> 0`, `< 1` (Courant number — never raise to get bigger steps) |
| `cfl_dg_safety` | real | `3.0` | `>= 1.0` (extra DG margin; `1.0` = SeisSol-equivalent, `dt×3`) |
| `use_pml`    | bool | `false` | none                                         |

---

## `[time]`

Time strings accept the `s` suffix (only): `"12s"`, `"0.001s"`,
`"1.5e-2s"`.  `"auto"` for `dt_initial` returns sentinel `-1`.

| Key          | Type        | Default  | Validation                          |
|--------------|-------------|----------|-------------------------------------|
| `tfinal`     | str or real | —        | `> 0`                               |
| `t_initial`  | real        | `0.0`    | `>= 0`                              |
| `dt_initial` | str or real | `"auto"` | `> 0` OR string `"auto"`            |
| `dt_max`     | str or real | `"0.1s"` | `> 0`                               |

---

## `[output]`

| Key                          | Type | Default                  | Validation       |
|------------------------------|------|--------------------------|------------------|
| `output_dir`                 | str  | —                        | non-empty        |
| `restart_prefix`             | str  | `"cp"`                   | non-empty        |
| `paraview_volume`            | str  | `"hdf5"`                 | `hdf5`/`vtu`/`off` |
| `paraview_bulk`              | str  | `"hdf5"`                 | `hdf5`/`vtu`/`off` |
| `paraview_fault`             | str  | `"hdf5"`                 | `hdf5`/`vtu`/`off` |
| `paraview_volume_dt`         | str  | `"0.05s"`                | `> 0`            |
| `paraview_bulk_dt`           | str  | `"0.05s"`                | `> 0`            |
| `paraview_fault_dt`          | str  | `"0.001s"`               | `> 0`            |
| `paraview_volume_zfp_tol`    | real | `1e-3`                   | `>= 0`           |
| `paraview_bulk_zfp_tol`      | real | `1e-3`                   | `>= 0`           |
| `paraview_fault_zfp_tol`     | real | `1e-12`                  | `>= 0`           |
| `max_snapshots`              | int  | `5000`                   | `>= 1`           |
| `checkpoint_every_steps`     | int  | `10000`                  | `>= 1`           |

---

## `[nucleation]` (optional)

Time-domain nucleation perturbation.  When the entire block is absent,
the parser sets `NucleationSpec::enabled = false` and the driver runs
without any nucleation perturbation (the rupture is driven entirely by
the static initial stress field).  When the block is present, the
parser sets `enabled = true`, validates `kind`, and reads the matching
sub-block.

| Key    | Type   | Default | Validation                  |
|--------|--------|---------|-----------------------------|
| `kind` | string | —       | `"gradual_overstress"` (only supported kind) |

### `[nucleation.gradual_overstress]` (required when `kind = "gradual_overstress"`)

Per-DOF shear-stress perturbation that ramps smoothly from 0 to the
full amplitude `Δτ · F(r)` over the interval `[0, T_nuc_s]`, where:

- `F(r)` is a Gaussian spatial factor centred on `(center_x_m,
  center_y_m, center_z_m)` with e-fold radii `radius_dip_m` (down-dip
  direction) and `radius_strike_m` (along-strike direction);
- the temporal ramp is the SCEC smoothStep function
  `exp(τ²/(t·(t − 2·T_nuc_s)))` on `(0, T_nuc_s)`, 0 below, 1 above —
  `C∞` everywhere except at `t = 0`.  `T_nuc_s` plays both roles of
  "ramp duration" and smoothStep `t0` parameter (R-007);
- per sub-step the resolver writes the increment
  `ΔS(t, Δt) · F(r) · Δτ` into `DOFData[i].tau1_nuc` (dip component)
  and `DOFData[i].tau2_nuc` (strike component), so that summed over
  `[0, T_nuc_s]` the accumulated perturbation telescopes to the full
  `Δτ · F(r)`.

| Key                   | Unit | Default | Validation                                |
|-----------------------|------|---------|-------------------------------------------|
| `center_x_m`          | m    | —       | required; in mesh bbox                    |
| `center_y_m`          | m    | —       | required; in mesh bbox                    |
| `center_z_m`          | m    | —       | required; typically `< 0` (depth)         |
| `radius_dip_m`        | m    | —       | required; `> 0`; recommended `≥ 3 · h_min`|
| `radius_strike_m`     | m    | —       | required; `> 0`; recommended `≥ 3 · h_min`|
| `delta_tau_dip_pa`    | Pa   | `0.0`   | none                                      |
| `delta_tau_strike_pa` | Pa   | `0.0`   | none                                      |
| `T_nuc_s`             | s    | —       | required; `> 0`; recommended `≤ tfinal/10`.  Plays both roles of "ramp duration" and smoothStep `t0` parameter (the redundant `t0_smooth_s` field was dropped per PLAN_first_safs_run R-007). |

The Gaussian spatial factor is

```
F(r) = exp( - ((dx/radius_dip_m)^2 + (ds/radius_strike_m)^2) )
```

where `dx` and `ds` are the DOF-to-centre offsets resolved in the
fault-local dip and strike directions (using each DOF's local
`FaultBasis`).  At `r = 0`, `F = 1`; outside ~`3 · radius_*`, `F` is
numerically zero and the accumulator is a no-op on those DOFs.

The per-sub-step smoothStep increment is

```
ΔS(t, Δt) = smoothStep(t, t₀_smooth) − smoothStep(t − Δt, t₀_smooth)
```

with `smoothStep(t, t₀) = 0` for `t ≤ 0`, `exp(τ²/(t·(t − 2·t₀)))`
for `0 < t < t₀`, and `1` for `t ≥ t₀` (where `τ = t − t₀`).  At any
`t ≥ T_nuc_s` the perturbation has telescoped to its full value
`F(r) · |Δτ|` and the accumulator is a no-op for the rest of the run.

Rate-and-state runs (`law = "rate_state"`) ignore the `[nucleation]`
block entirely — nucleation in the rate-state path is achieved via
the `V_init` field, not via a stress accumulator.

---

## `[friction]` (optional; law-agnostic strength floor)

A bare `[friction]` table (declared **before** the law sub-block) carries
one optional key that applies to **both** friction laws:

| Key                         | Unit | Default | Validation                                        |
|-----------------------------|------|---------|---------------------------------------------------|
| `sigma_n_strength_floor_pa` | Pa   | `-1.0`  | when present: finite **and** `>= 0` (a negative value is the disabled sentinel and must be expressed by OMITTING the key) |

The **compressive normal-stress strength floor** (sliver-blowup plan
2026-05-26).  The normal stress entering the **shear strength** is
`max(sigma_n_total, sigma_n_strength_floor_pa)`; above the floor the
strength is the usual `mu * sigma_n` (proportional), and below the floor
it saturates at the cohesion-like constant `mu * floor` instead of the
spurious tensile `0` (LSW `max(sigma_n, 0)`) or `|sigma_n|` (RS).  This
breaks the `sigma_n -> strength -> radiation` feedback that drives the
tensile free-slip runaway (job 7748818, debug doc §2c).

- **Disabled by default** (`< 0` sentinel): each law keeps its exact
  current strength expression, so the TPV205/102/104 and BP5 byte-exact
  regressions are untouched.
- Only the strength's `sigma_n` is floored; the written-back
  `normal_stress` output channel and the friction-solver `sigma_n`
  argument are NOT floored (v1; plan §Phase 3 decision 2).
- The SAFS configs set `sigma_n_strength_floor_pa = 10.0e6` (10 MPa).
- **The floor must be a TOP-LEVEL `[friction]` key.**  Nesting it under
  `[friction.slip_weakening]` / `[friction.rate_state]` is a hard error
  (the parser aborts), because the per-law sub-block parsers would
  otherwise silently ignore it and leave the floor disabled.
- **Effect differs by law.**  On the **LSW** path the floor directly
  bounds `tau_strength = mu_eff * max(sigma_n, floor)`, so it actively
  bounds the slip rate `V_abs = (|tau| - tau_strength)/eta_s` under a
  tensile excursion — this is the blow-up fix.  On the **rate-and-state**
  path the floor is applied only to the post-solve slip-rate
  *decomposition* strength; the Brent solver that sets the slip-rate
  *magnitude* `V_abs` still uses the unfloored `|sigma_n|` (decision 2).
  Because a tensile `|sigma_n|` makes RS friction *stronger* (so the
  solver locks `V_abs ≈ 0` and the decomposition strength is never
  reached), the floor is largely inert on the RS path — RS does not have
  the LSW tensile free-slip failure mode in the first place.  Enabling it
  on RS configs is harmless but mostly a no-op; do not treat it as an
  active tensile-slip safeguard for RS.

---

## `[friction.slip_weakening]` (required when `law = "slip_weakening"`)

| Key                | Unit | Default | Validation                  |
|--------------------|------|---------|-----------------------------|
| `mu_s_default`     | —    | `1.1`   | `> 0`; **no upper bound** (D-3) |
| `mu_d_default`     | —    | `0.5`   | `> 0`; `mu_d_default < mu_s_default` |
| `d_c_default`      | m    | `0.5`   | `> 0`; alias `d_o_default` accepted (D-3, deprecation notice on first use) |
| `cohesion_default` | Pa   | `0.0`   | `>= 0`                      |

Defaults match `friction/slip-weakening/geoffrey2010.md` verbatim
(R-113).

### `[[friction.slip_weakening.spatial]]` (zero or more rules; document
order; last-match wins per key)

Common keys for all spatial rules:

| Key            | Type | Default | Notes                                     |
|----------------|------|---------|-------------------------------------------|
| `kind`         | str  | `"depth"` | `"depth"`, `"box"`, `"region_attribute"`, `"barrier"` |
| `x_min_m`      | real | `-inf`  | used by `box` only                        |
| `x_max_m`      | real | `+inf`  | used by `box` only                        |
| `y_min_m`      | real | `-inf`  | used by `box` only                        |
| `y_max_m`      | real | `+inf`  | used by `box` only                        |
| `z_min_m`      | real | `-inf`  | used by `depth`/`box`/`barrier`           |
| `z_max_m`      | real | `+inf`  | used by `depth`/`box`/`barrier`           |
| `region_attr`  | int  | `-1`    | used by `region_attribute` only           |

Per-key overrides (NaN sentinel = "do not override"):

- `mu_s`, `mu_d`, `d_c` (alias `d_o`), `cohesion`

Nucleation is entirely time-domain (`[nucleation]` block, kind =
`gradual_overstress`); friction spatial rules NEVER override
`tau_pre_*` or `sigma_n`.  (The `nucleation_box` kind that existed in
rev-1 and rev-2 was removed in rev-3.)

**Barrier (R-114):** `kind = "barrier"` locks the DOF.  The resolver
internally sets `mu_s = 1.0e6` (the existing `fault_face_flux.hpp`
sentinel).  Users NEVER type `1.0e6` themselves: the parser REJECTS
any user-supplied `mu_s > 1.0e5` with the message "to mark a barrier
region, use `kind = \"barrier\"` instead of a large mu_s".

---

## `[friction.rate_state]` (required when `law = "rate_state"`)

| Key                | Unit | Default  | Validation               |
|--------------------|------|----------|--------------------------|
| `f_0_default`      | —    | `0.6`    | `0 < f_0 < 1`            |
| `V_0_default`      | m/s  | `1.0e-6` | `> 0`                    |
| `eta`              | —/str| `"auto"` | `"auto"` (per-DOF eta = `0.5 * sqrt(mu * rho)`) OR a positive real; R-011 |
| `a_default`        | —    | `0.010`  | `> 0`                    |
| `b_default`        | —    | `0.015`  | `> 0`; `a_default < b_default` **only when `[friction.rate_state.depth_profile]` is absent** (Phase 11b) |
| `Dc_default`       | m    | `0.004`  | `> 0`                    |
| `V_init_default`   | m/s  | `1.0e-9` | `> 0`                    |
| `sigma_n_default`  | Pa   | `50.0e6` | `> 0`                    |

### `[[friction.rate_state.spatial]]` (same kind set as LSW; no
`nucleation_box`)

Per-key overrides (NaN sentinel = "do not override"):

- `a`, `b`, `Dc`, `V_init`, `eta`, `sigma_n` — **allowed** (per-DOF). `b` was
  rejected pre-Phase-11a; it is now plumbed through `DOFData.b` + the aging-law
  ψ-update + the equilibrium seed.
- `f_0`, `V_0` — **rejected** (scalar aging-law globals; `V_0` is pinned to
  `FrictionSolver::V0` by the R-009 guard). Set them only in the defaults above.

Rate-and-state nucleation is achieved through the `V_init` field, NOT
through `[nucleation]` (the `[nucleation]` block is for LSW
`gradual_overstress` only and is ignored when `law = "rate_state"`).

### `[friction.rate_state.depth_profile]` (optional; Phase 11b)

Depth-varying `a(z)` and `b(z)` supplied as **two CSV files**, linearly
interpolated onto each fault DOF at `depth = max(0, -z)` (metres). When present,
the resolver seeds per-DOF `a`/`b` from the profile instead of the scalar
`a_default`/`b_default` (those become an unused fallback, and the
`a_default < b_default` default check is skipped). Per-DOF `[[…spatial]]`
overrides still apply on top (last-match-wins).

| Key                   | Type | Default | Validation                                  |
|-----------------------|------|---------|---------------------------------------------|
| `param_a_csv`         | str  | —       | required, non-empty; file readable at parse time |
| `param_a_minus_b_csv` | str  | —       | required, non-empty; file readable at parse time |
| `depth_units`         | str  | `"km"`  | `"km"` (depth × 1000 → m) or `"m"`          |

Each CSV: no header, comma- OR whitespace-separated, **two columns per row
`value, depth` (value FIRST, depth SECOND)**; `#` comments and blank lines are
skipped; ≥ 2 rows; duplicate depths abort. `param_a.csv` holds `a` (must be
`> 0`); `param_a_minus_b.csv` holds `a − b` (may be negative — velocity-
weakening). The two files are interpolated **independently** (they may use
different depth grids), then combined:

```
a(z) = interp(param_a)
b(z) = a(z) − interp(param_a_minus_b)
```

Both curves are **flat-clamped** (constant) outside their sampled depth range —
e.g. a fault deeper than the CSVs extends the deepest sampled value (so a future
deeper re-mesh needs no profile change). The resulting per-DOF `b > 0` is
enforced by the resolver (`a − (a−b) ≤ 0` aborts, naming `param_a_minus_b.csv`).
`--print-derived` echoes the CSV paths, the sampled depth ranges, the mesh fault
max depth, and the VW↔VS transition depth (where `a − b = 0`).

---

## Validation rules (parser MFEM_VERIFYs each; complete enumeration)

1. `meta.schema_version == 1`; abort with the offending value.
2. `meta.law ∈ {"slip_weakening", "rate_state"}`.
3. Either `[friction.slip_weakening]` or `[friction.rate_state]` is
   present (not both).  Extra keys in the non-matching block →
   MFEM_ABORT.
4. **Friction validator (D-3 relaxed):**
   - LSW: `0 < mu_d`, `mu_d < mu_s`, `d_c > 0`, `cohesion >= 0` at
     every DOF after spatial rules.  **No upper bound on `mu_s`.**
   - Rate-state: `0 < a`, `0 < b`, `Dc > 0`, `V_0 > 0`, `sigma_n_eff > 0`,
     `0 < f_0 < 1` at every DOF.  **`a < b` is NOT required per-DOF** (R-011):
     velocity-strengthening regions (`a > b`) are allowed (fault-edge / deep
     arrest, depth-profile VS tapers).  The scalar `a_default < b_default` check
     applies only to the fallback when no depth profile is configured.
5. **`d_o` alias (D-3):** accepts `d_o_default` as a synonym for
   `d_c_default`, and `d_o` as a synonym for `d_c` inside spatial
   rules.  Emits a one-line `mfem::out` deprecation notice on first
   use.
6. **Barrier-sentinel guard (R-114):** if a user supplies `mu_s > 1.0e5`
   in any spatial rule's per-DOF override, abort with the message
   above.
7. **Pore pressure:** `min_sigma_n_pa >= 0`; `P_p_grad_pa_per_m >= 0`
   is expected (warn but do not abort on negative).
8. **Material fallback:** `mu > 0`, `rho > 0`, `(lambda + 2*mu) > 0`.
9. **Stress mode (D-1):** if `kind = "constant_tensor"`, the six
   `sigma_*_pa` keys must be present and `sidecar_path` must be
   absent (or empty).  If `kind = "sidecar_hdf5"`, `sidecar_path`
   must be a non-empty string AND the six `sigma_*_pa` keys must be
   absent.
10. **Time:** `tfinal > 0`, `dt_max > 0`.  `tfinal` parses successfully
    via `spatial_time_parser` (otherwise MFEM_ABORT with the offending
    string).
11. **Numerics:** `cfl > 0`, `cfl < 1`; `cfl_dg_safety >= 1.0` (the extra
    DG safety margin beyond the mandatory order factor `1/(2N+1)`; `1.0` =
    SeisSol's validated step); `ader_order >= 1`; `mixed_flux` is one of the
    three accepted strings.
12. **Output:** `output_dir` non-empty; `paraview_*` modes are one of
    `hdf5`/`vtu`/`off`; ZFP tolerances `>= 0`; `max_snapshots >= 1`;
    `checkpoint_every_steps >= 1`.
13. **Nucleation:** when the `[nucleation]` block is present,
    `kind = "gradual_overstress"` (the only supported kind) AND the
    `[nucleation.gradual_overstress]` sub-block is populated with
    `center_*_m`, `radius_*_m > 0`, and `T_nuc_s > 0`.  When absent,
    `NucleationSpec::enabled = false` and the driver runs without any
    nucleation perturbation.
14. **Sigma_n strength floor:** when `[friction].sigma_n_strength_floor_pa`
    is present it must be finite and `>= 0`; an explicit negative value
    aborts (the disabled state is expressed by OMITTING the key, which
    defaults the field to `-1.0`).

---

## EXAMPLE files

Three `EXAMPLE_*.toml` files are committed under
`safs/project_7.0_alternative/friction/`:

1. `EXAMPLE_spatial_friction_slip_weakening_safs.toml` — LSW with
   `geoffrey2010.md` defaults (`mu_s=1.1, mu_d=0.5, d_c=0.5`), the
   `[nucleation]` block populated with `kind = "gradual_overstress"`,
   and a single barrier rule.
2. `EXAMPLE_spatial_friction_rate_state_safs.toml` — rate-and-state
   (Aging law) with BP5-style scalars + a single VW box rule.
3. `EXAMPLE_spatial_stress_constant_tensor_safs.toml` — minimum-viable
   constant-tensor `[stress]` block with a representative
   SAFOD-style tensor (SHmax ≈ 100 MPa NE-ward).

All three load through `LoadSpatialFrictionConfig(...)` without error
(Phase 1 unit test `test_spatial_friction_config`).
