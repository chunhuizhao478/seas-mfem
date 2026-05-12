# Plan: Split bulk solutions into independent kinematics + stress files

**Date:** 2026-05-12
**Status:** active (companion to PLAN_paraview_compaction_2026-04-28.md + PLAN_bulk_compression_and_size_estimator_2026-05-09.md)
**Authors:** Phase 6 follow-up.

---

## Motivation

After Phase 6 landed, every TPV* run writes three `.vtkhdf` files:

- `volume.vtkhdf` (primary) — velocity + mpi_rank + L2-p0 fault projections (slip_dip, slip_rate_*, traction_*, state, normal_stress, friction params a/Dc, x2, x3)
- `wave_bulk.vtkhdf` (secondary) — velocity + mpi_rank + sigma_yy + sigma_xy + sigma_xz
- `fault_surface.vtkhdf` — full fault DOFData on a 2-D triangulated surface

This produces two structural problems:

1. **`velocity` and `mpi_rank` are stored twice.** Both `pv_out` and `pv_bulk_out` register the same GridFunction pointers (`tpv102_driver.cpp:1921,1928` and `:2116,2117`). When `--paraview-dt == --paraview-bulk-dt` and the two ZFP tolerances match — the Phase 6 sbatch default — the velocity bytes in `wave_bulk.vtkhdf` are byte-for-byte identical to those in `volume.vtkhdf`. Roughly half of `wave_bulk.vtkhdf` is wasted.

2. **`volume.vtkhdf` carries fault-projected L2-p0 fields that are already in `fault_surface.vtkhdf`.** `UpdateFaultFieldsBP5` populates 8+ scalar GFs on the volume mesh restricted to elements adjacent to a fault face; the same per-DOF data lands in `fault_surface.vtkhdf`'s triangulated representation. The volume-side copy is a debug rendering aid that few users actually look at; in production it's a fixed-cost duplication.

3. **Stress fields are incomplete.** `pv_bulk_out` writes only 3 of the 6 independent components of the symmetric stress tensor (yy/xy/xz). The full tensor (xx,yy,zz,xy,xz,yz) is in the solver state vector at zero additional compute cost; only the registration was incomplete.

4. **Names don't say what's inside.** "volume" and "wave_bulk" are physically vague; "fault_surface" is OK but inconsistent with the proposed rename of the other two.

---

## Design

### File names

| Old name | New name | Path |
|---|---|---|
| `volume.vtkhdf` | `kinematics.vtkhdf` | `<output>/ParaView_kinematics/kinematics.vtkhdf` (was `<output>/ParaView/`) |
| `wave_bulk.vtkhdf` | `stress.vtkhdf` | `<output>/ParaView_stress/stress.vtkhdf` (was `<output>/ParaView_bulk/`) |
| `fault_surface.vtkhdf` | `fault.vtkhdf` | `<output>/fault.vtkhdf` (top-level, unchanged dir) |

### Default field contents (TPV*)

| Collection | Fields |
|---|---|
| `kinematics.vtkhdf` | `velocity` (3-vec), `mpi_rank` (cell scalar) |
| `stress.vtkhdf` | `sigma_xx`, `sigma_yy`, `sigma_zz`, `sigma_xy`, `sigma_xz`, `sigma_yz` (full symmetric tensor; 6 cell scalars) |
| `fault.vtkhdf` | unchanged — emitted by `WriteFaultSurfaceVTU` with the existing field set |

The fault L2-p0 projections on the volume mesh (slip_dip, slip_rate_*, traction_*, state, friction-params) are dropped from the primary collection. They were debug renderings of data already in `fault.vtkhdf`.

### BP5

BP5 has no `pv_bulk_out`; its single volume collection currently writes `displacement` + fault L2-p0 projections.

For consistency: rename the BP5 collection to `kinematics` too, drop the L2-p0 projections (they're in `fault.vtkhdf`), keep `displacement` as the BP5 "kinematic" field. No stress file.

### CLI flag renames (tpv* drivers)

| Old flag | New flag |
|---|---|
| `--paraview-volume-vtu` / `--paraview-volume-hdf5` | `--paraview-kinematics-vtu` / `--paraview-kinematics-hdf5` |
| `--paraview-volume-zfp-tol` / `--paraview-volume-deflate-level` | `--paraview-kinematics-zfp-tol` / `--paraview-kinematics-deflate-level` |
| `--paraview-bulk-dt` | `--paraview-stress-dt` |
| `--paraview-bulk-zfp-tol` / `--paraview-bulk-deflate-level` | `--paraview-stress-zfp-tol` / `--paraview-stress-deflate-level` |
| `--no-volume-pv` | `--no-kinematics-pv` (with `--no-volume-pv` kept as deprecated alias) |

Per-collection field selectors (new):

- `--paraview-kinematics-fields=velocity,mpi_rank` — comma-separated names; empty = all defaults.
- `--paraview-stress-fields=sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz` — same.
- `--paraview-fault-fields=...` — already exists via `SetFaultVTUFields`; no change.

A flag mismatch (unknown field name) aborts with the available-name list.

### TOML schema (future seas_driver wiring; deferred per R-313)

```toml
[paraview]
enabled = true

[paraview.kinematics]
fields    = ["velocity", "mpi_rank"]
zfp_tol   = 1e-3              # absolute-error tolerance, m/s for velocity
# deflate_level = 6           # alternative to zfp_tol; mutually exclusive
dt        = 0.5               # seconds; 0 = adaptive
backend   = "hdf5"            # or "vtu"

[paraview.stress]
fields        = ["sigma_xx", "sigma_yy", "sigma_zz",
                 "sigma_xy", "sigma_xz", "sigma_yz"]
zfp_tol       = 1e-3          # Pa; stresses span ~1e6..1e8 so 1e-3 is overkill but safe
dt            = 0.5
# enabled     = false         # drop the file entirely

[paraview.fault]
zfp_tol       = 1e-12         # m/s; slip-rate floor preserved
deflate_level = -1            # -1 = use ZFP, else deflate level 0..9
backend       = "hdf5"
fields        = []            # empty = all available; subset = filter
```

`seas_driver.cpp` currently doesn't construct `seas::ParaViewOutput` (Phase 4 deferred deviation R-313 documented in `miniapps/seas/CLAUDE.md`). The TOML schema lands here as a future-use spec; the parsing wiring follows when BP5 production runs through `seas_driver` actually need a `kinematics.vtkhdf`.

---

## Implementation phases

| # | Phase | Files touched |
|---|---|---|
| 1 | tpv102 canonical refactor | `drivers/tpv102_driver.cpp` |
| 2 | Mirror to tpv104, tpv205 | `drivers/tpv104_driver.cpp`, `drivers/tpv205_driver.cpp` |
| 3 | Rename fault file in writer | `io/fault_vtkhdf_writer.hpp` (one line) + any test that hardcodes the filename |
| 4 | BP5 driver rename | `tests/verification/bp5_verification_full.cpp` (one line) |
| 5 | Sbatch updates | `jobs/{bp5,tpv102,tpv104,tpv205}/*phase6*.sbatch` (CLI flag names) |
| 6 | TOML schema (docs only) | new `config/paraview_example.toml` |
| 7 | CLAUDE.md + walkthrough doc | `CLAUDE.md`, `document/codebase_explorer/paraview_io_walkthrough_2026-05-11.{md,pdf}` |
| 8 | Local build verification | all four drivers |

No backwards-compat shims for the volume/bulk flag names — per CLAUDE.md feedback, just rename. Sbatch files in this repo are updated in the same commit; out-of-tree sbatch will break loudly at parse time (better than silent semantic change).

`--no-volume-pv` keeps an alias because it's documented in `miniapps/seas/CLAUDE.md` and likely lives in user sbatch outside this repo.

---

## Acceptance criteria

- [ ] All four drivers build locally under mfem-dev conda env (no warnings about renamed flags).
- [ ] A TPV102 run on a small mesh emits exactly: `ParaView_kinematics/kinematics.vtkhdf`, `ParaView_stress/stress.vtkhdf`, `fault.vtkhdf` (and station `.dat`).
- [ ] `kinematics.vtkhdf` contains only `velocity` + `mpi_rank` (verifiable via `h5ls -r`).
- [ ] `stress.vtkhdf` contains all 6 sigma components.
- [ ] No field-name overlap between `kinematics.vtkhdf` and `stress.vtkhdf`.
- [ ] `--paraview-volume-zfp-tol` aborts at parse time with "did you mean --paraview-kinematics-zfp-tol?".
- [ ] BP5 driver emits `kinematics.vtkhdf` (single collection, no stress).
- [ ] All four sbatch files updated and pass `bash -n`.
- [ ] CLAUDE.md and walkthrough doc reflect new names.

---

## Out of scope (for this plan)

- Wiring TOML parsing into `seas_driver` — deferred per R-313 until BP5 production through seas_driver needs PV.
- Backwards-compat alias for `--paraview-volume-*` and `--paraview-bulk-*` flag names — just rename.
- Migration tool to rename existing `.vtkhdf` files on disk — they're regenerable.
