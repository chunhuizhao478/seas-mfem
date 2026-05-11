# Data Projection Schema v1

**Schema id:** `data_projection_v1`
**Status:** Stable. Any breaking change bumps the schema version (`v2`, `v3`,
…) and the runtime `DataField3D` loader rejects mismatched versions.

This document is the single source of truth for the HDF5 sidecar layout
consumed by the data-projection feature (see
`data_projection_feature_plan_v2.md`). The Phase 2 writer
(`code_preprocess/data_projection/sidecar.py`) and the Phase 3 reader
(`miniapps/seas/io/data_field_3d.hpp`) cite this document.

---

## 1. Layout

```
<sidecar>.h5
├── attrs (root):
│     schema_version = "data_projection_v1"          (string, required)
│     crs            = "EPSG:32611"                   (string, required)
│     units          = "m"                            (string, required)
│     z_positive     = "elevation"                    (string, required)
│     created_at     = ISO-8601 timestamp             (string, required)
│     source         = "<comma-separated filenames>"  (string, optional)
│     source_crs     = "EPSG:4326"                    (string, optional)
│     mesh_tag       = "<informational tag>"          (string, optional)
├── /grid
│     ├── x  [Nx]  float64, strictly monotone increasing
│     ├── y  [Ny]  float64, strictly monotone increasing
│     └── z  [Nz]  float64, strictly monotone increasing
└── /fields
      ├── <name1>  [Nx, Ny, Nz]  float64
      │     attrs:  units      = "<unit>"   (string, required)
      │             min_value  = <float>    (required)
      │             max_value  = <float>    (required)
      ├── <name2>  ...
      └── ...
```

## 2. Coordinate convention

- **CRS:** UTM zone 11 N, WGS-84 datum (EPSG:32611). All `(x, y)`
  coordinates are in metres east / north. The reserved `"local-tangent"`
  value is not implemented in v1.
- **Units:** every spatial coordinate is in metres.
- **Z direction:** `z_positive = "elevation"`. Values increase upward.
  `z = 0` is the free surface; `z < 0` is at depth. Source datasets that
  store "depth, positive down" must be flipped to elevation by the
  writer (`z_canon = -depth_m`).

## 3. Index order and shape

- All field arrays are 3-D, row-major (HDF5 C order), with shape
  exactly `(len(/grid/x), len(/grid/y), len(/grid/z))`.
- The runtime loader checks shape match and the strict-monotone
  property of each axis, and aborts otherwise.
- Single-cell axes (e.g. `Nx == 1`) are permitted by the schema but
  cause the runtime to treat any off-axis query as out-of-bbox.

## 4. Value constraints (writer-side AND reader-side)

The schema **forbids**:

1. **NaN** anywhere in any `/fields/<name>` dataset.
2. Any cell value outside the per-field declared `[min_value,
   max_value]` interval.

The Phase 2 writer enforces both via guards `G-1` and `G-2`; the
Phase 3 reader independently re-checks both at load time. Both
violations are hard aborts with a precise error message.

### 4.1 Out-of-bounds (OOB) policy at evaluation time

The runtime is hard-wired to `OOBPolicy::Abort`. There is no clamp
or extrapolation in v1. A `(x, y, z)` query that falls outside
`[/grid/x[0], /grid/x[Nx-1]] × [...] × [...]` causes
`DataField3D::Evaluate` to abort with the offending coordinate.

A future schema bump (`data_projection_v2`) may re-introduce a
`clamp` policy with an explicit attribute and runtime opt-in.

## 5. Field-set conventions

The following named sets are recognised as canonical bundles. Adding a
new field is just a matter of writing a new `/fields/<name>` dataset;
no schema change is required.

| Field set       | Field names                                                       | Default units |
|-----------------|-------------------------------------------------------------------|---------------|
| `velocity`      | `Vp`, `Vs`, `density`                                             | m/s, m/s, kg/m³ |
| `stress`        | `sigma_xx`, `sigma_yy`, `sigma_zz`, `sigma_xy`, `sigma_yz`, `sigma_xz` | Pa each |
| `pore_pressure` | `p`                                                               | Pa            |

### 5.1 Default sanity bounds for the velocity set

The defaults below are **permissive lower bounds** that accept the
minima observed in the SCEC CVM-H 15.1.1 archive (Vs ~ 120 m/s,
Vp ~ 1200 m/s, ρ ~ 1420 kg/m³ in water-saturated near-coast cells).
Upper bounds remain at the competent-rock ceiling.  Callers may
tighten the lower bound per-field via the writer CLI flags
(`--vp-min-mps`, `--vs-min-mps`, `--rho-min-kgm3`) or by passing
explicit `field_bounds` to `write_sidecar`.

| Field    | min_value (default) | max_value (default) | Notes                                                                    |
|----------|---------------------|---------------------|--------------------------------------------------------------------------|
| `Vp`     | 1000 m/s            | 9000 m/s            | Permissive lower bound; tighten to 2000 m/s to reject sub-seafloor cells.|
| `Vs`     | 100 m/s             | 5000 m/s            | Permissive lower bound; tighten to 1500 m/s to reject sub-seafloor cells.|
| `density`| 1000 kg/m³          | 3500 kg/m³          | Permissive lower bound; tighten to 2000 kg/m³ for competent-rock-only.   |

## 6. Required attributes

| Path / attribute            | Type   | Required | Notes                                                                |
|-----------------------------|--------|----------|----------------------------------------------------------------------|
| root: `schema_version`      | string | yes      | Must equal `"data_projection_v1"`.                                   |
| root: `crs`                 | string | yes      | Must equal `"EPSG:32611"`.                                           |
| root: `units`               | string | yes      | Must equal `"m"`.                                                    |
| root: `z_positive`          | string | yes      | Must equal `"elevation"`.                                            |
| root: `created_at`          | string | yes      | ISO-8601 (`YYYY-MM-DDTHH:MM:SSZ`).                                   |
| root: `source`              | string | no       | Comma-separated source filenames for provenance.                      |
| root: `source_crs`          | string | no       | CRS of the source data before re-projection.                          |
| root: `mesh_tag`            | string | no       | Informational mesh stem (e.g. `safs_fault_box_nwcut_500m`).          |
| `/grid/x`, `/grid/y`, `/grid/z` | dataset (float64, 1-D) | yes      | Strictly monotone increasing.                                         |
| `/fields/<name>`            | dataset (float64, 3-D) | ≥ 1 required | Shape exactly `(Nx, Ny, Nz)`.                                       |
| `/fields/<name>/units`      | string | yes      | Per-field unit string.                                                |
| `/fields/<name>/min_value`  | float  | yes      | Sanity lower bound enforced at writer + reader.                      |
| `/fields/<name>/max_value`  | float  | yes      | Sanity upper bound enforced at writer + reader.                      |

## 7. Schema-version policy

The runtime `DataField3D` constructor reads the root `schema_version`
attribute and aborts unless it equals `"data_projection_v1"` exactly.
A future change to any of:

- the OOB policy (e.g. opt-in `clamp`)
- the index order (e.g. switch to `(Nz, Ny, Nx)`)
- the value-constraint contract (e.g. allow NaN with a fill attribute)
- removal of a previously-required attribute

bumps the schema id to `data_projection_v2`, etc. Loaders may support
multiple historical schemas via dispatch on the version attribute.

## 8. Provenance traceability

The optional `source` and `source_crs` attributes are recommended for
real-data sidecars. They are not consumed by the runtime; their
purpose is to make a sidecar self-describing for humans inspecting it
with `h5dump -A`.
