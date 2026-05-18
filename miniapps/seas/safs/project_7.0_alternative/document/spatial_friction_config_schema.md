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
[nucleation]                   # optional (when absent, no forced rupture)
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
| `model`         | string | —             | `"cvmh"`, `"cvm_s4.26.m01"`, `"multiscale_statewise"` |
| `dataset_root`  | string | —             | non-empty when `override_path` is empty          |
| `override_path` | string | `""`          | when non-empty, bypasses model-based resolution  |

Path resolution (rev-3, R-003 corrected):

```
cvmh                 → <dataset_root>/velocity/results/cvmh/velocity_safs.h5
cvm_s4.26.m01        → <dataset_root>/velocity/results/cvm_s4.26.m01/velocity_safs.h5
multiscale_statewise → <dataset_root>/velocity/results/multiscale_statewise_cvm/velocity_safs.h5
```

There is **no `mesh_tag`** in the filename — one sidecar per CVM model.

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
| `cfl`        | real | `0.5`   | `> 0`, `< 1`                                 |
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

## `[nucleation]` (D-4, optional)

TPV26/27-style gradual forced rupture.  When the entire block is
absent, the parser sets `NucleationSpec::enabled = false` and the
driver runs WITHOUT forced rupture.  When the block is present, the
parser sets `enabled = true` and the keys below are required.

| Key                | Unit | Default  | Validation               |
|--------------------|------|----------|--------------------------|
| `hypocenter_x_m`   | m    | —        | none                     |
| `hypocenter_y_m`   | m    | —        | none                     |
| `hypocenter_z_m`   | m    | —        | typically `< 0` (depth)  |
| `r_crit_m`         | m    | `4000.0` | `> 0`                    |
| `t0_decay_s`       | s    | `0.5`    | `> 0`                    |

The per-DOF time-of-forced-rupture `T(r)` is computed in Phase 1
`ResolveForcedRupture` from the local shear-wave speed `Vs(r)` per
TPV26/27 §Part 5 (variable-speed rupture front):

```
T(r) = r / (0.7 Vs) + 0.081 r_crit / (0.7 Vs) * (1 / (1 - (r/r_crit)^2) - 1)
                                                     for r < r_crit
T(r) = 1.0e9                                          for r >= r_crit
```

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

**Rev-3 D-4 removes** the `nucleation_box` kind that existed in rev-1
and rev-2.  Nucleation is now entirely time-domain (`[nucleation]`
block); spatial rules NEVER override `tau_pre_*` or `sigma_n`.

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
| `b_default`        | —    | `0.015`  | `> 0`; `a_default < b_default` |
| `Dc_default`       | m    | `0.004`  | `> 0`                    |
| `V_init_default`   | m/s  | `1.0e-9` | `> 0`                    |
| `sigma_n_default`  | Pa   | `50.0e6` | `> 0`                    |

### `[[friction.rate_state.spatial]]` (same kind set as LSW; no
`nucleation_box`)

Per-key overrides (NaN sentinel = "do not override"):

- `a`, `b`, `Dc`, `V_init`, `f_0`, `V_0`, `eta`, `sigma_n`

Rate-and-state nucleation is achieved through the `V_init` field, NOT
through `[nucleation]` (the `[nucleation]` block is for LSW forced
rupture only and is ignored when `law = "rate_state"`).

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
   - Rate-state: `0 < a`, `0 < b`, `a < b`, `Dc > 0`, `V_0 > 0`,
     `sigma_n_eff > 0`, `0 < f_0 < 1` at every DOF.
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
11. **Numerics:** `cfl > 0`, `cfl < 1`; `ader_order >= 1`; `mixed_flux`
    is one of the three accepted strings.
12. **Output:** `output_dir` non-empty; `paraview_*` modes are one of
    `hdf5`/`vtu`/`off`; ZFP tolerances `>= 0`; `max_snapshots >= 1`;
    `checkpoint_every_steps >= 1`.
13. **Nucleation:** when the `[nucleation]` block is present,
    `r_crit_m > 0` and `t0_decay_s > 0`.  When absent,
    `NucleationSpec::enabled = false`.

---

## EXAMPLE files

Three `EXAMPLE_*.toml` files are committed under
`safs/project_7.0_alternative/friction/`:

1. `EXAMPLE_spatial_friction_slip_weakening_safs.toml` — LSW with
   `geoffrey2010.md` defaults (`mu_s=1.1, mu_d=0.5, d_c=0.5`), the
   `[nucleation]` block populated, and a single barrier rule.
2. `EXAMPLE_spatial_friction_rate_state_safs.toml` — rate-and-state
   (Aging law) with BP5-style scalars + a single VW box rule.
3. `EXAMPLE_spatial_stress_constant_tensor_safs.toml` — minimum-viable
   constant-tensor `[stress]` block with a representative
   SAFOD-style tensor (SHmax ≈ 100 MPa NE-ward).

All three load through `LoadSpatialFrictionConfig(...)` without error
(Phase 1 unit test `test_spatial_friction_config`).
