# Implementation Plan — Spatial-Driver ParaView Output Parity with TPV104 Phase-6 Reference (2026-05-18)

## Overview

The `seas_spatial_dyn_driver` ParaView wiring (Phase 4 round-7 state) is a
**proper subset** of the Phase-6 paraview-compaction feature set that the
TPV* drivers consume via `jobs/tpv104/tpv104_phase6_paraview_zfp_dev_2hr.sbatch`.
This plan brings the spatial driver to **functional parity** with that sbatch's
flag surface, then extends fault-output to expose the SAFS spatial-friction
heterogeneity (per-DOF `mu_s`, `mu_d`, `d_c`, `T_forced`, initial stresses) so a
production rupture run can be visualised end-to-end without bouncing back to
ASCII dumps.

The work is **purely additive** to `drivers/spatial_dyn_driver.cpp`,
`spatial/code/spatial_friction.{hpp,cpp}`,
`safs/project_7.0_alternative/document/spatial_friction_config_schema.md`,
and (informationally) `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md`.
No `io/paraview_output.hpp` change — the API consumed already exists on this
branch (merged from `feature/paraview-compaction`, commits `240581e` and `0ce73a4`).
No TPV* / BP5 driver change — byte-exact contracts preserved.

## Why this matters now

The round-7 review (`spatial_dynamic_rupture_review_round7_2026-05-18.md`)
confirmed the dispatch chain (LSW_ForcedRupture → iterator → mu_eff) is
correctly wired on the production strength-reduction path.  The next user-
visible defect on the SAFS production path is **observability**: today, even
when the driver runs end-to-end, the operator gets `fault.vtkhdf` with only the
five BP5-style fault projections (slip, slip_rate, traction, mu_eff, sigma_n).
There is **no** volume velocity field, **no** secondary stress collection,
**no** regime-adaptive cadence, and **no** spatial-friction parameter map.
Compared to the TPV104 sbatch's flag set, the spatial driver is missing 14 CLI
flags, 6 HDF5/ZFP build guards, and an entire secondary `pv_bulk_out`
collection.

## What the TPV104 sbatch uses (the reference contract)

From `jobs/tpv104/tpv104_phase6_paraview_zfp_dev_2hr.sbatch:125-146`:

```
--paraview                            ← boolean enable
--paraview-dt 0.5                     ← fault snapshot cadence (umbrella)
--paraview-bulk-dt 0.5                ← secondary collection cadence
--paraview-volume-zfp-tol 1e-3        ← primary HDF5 ZFP on velocity
--paraview-bulk-zfp-tol 1e-3          ← secondary HDF5 ZFP on stress
--paraview-fault-zfp-tol 1e-12        ← fault HDF5 ZFP on slip-rate
--paraview-max-snapshots 200          ← AdaptiveSchedule cap
```

The driver behind those flags (`drivers/tpv104_driver.cpp:597-770` for parsing,
`drivers/tpv104_driver.cpp:1955-2330` for wiring) does substantially more:

1. **CLI surface** (12 additional knobs not in the sbatch but accepted by the
   parser): `--paraview-every N`, `--paraview-fault-{vtu,hdf5,legacy-ascii}`,
   `--paraview-{volume,bulk,fault}-deflate-level N`,
   `--paraview-volume-{vtu,hdf5}`,
   `--paraview-{coseismic,nucleation,interseismic}-dt X`.
2. **Build guards** — every `--paraview-*-{hdf5,zfp-tol,deflate-level}` flag
   verifies `MFEM_USE_HDF5`/`MFEM_USE_H5Z_ZFP` and aborts with a precise
   message instead of crashing in the writer.
3. **Two collections**: primary `kinematics.vtkhdf` (velocity + mpi_rank +
   BP5 fault projections) and secondary `stress.vtkhdf` (6 sigma components).
4. **Static fault parameter arrays**: `pv_out->SetFaultParamsBP5(a, Dc, x2, x3)`
   so ParaView users can colour by `a`/`Dc` directly.  For SAFS this list MUST
   extend to the spatial-friction quantities (which TPV104 keeps constant).
5. **`SetTotalRunTime(tfinal)`** so the AdaptiveSchedule's interseismic-budget
   rationing has a horizon to plan against.
6. **Rank-0 status banner** printing back the active mode + cadence.

## Constraints

- **`io/paraview_output.hpp` is read-only.**  This plan only consumes the
  paraview-compaction API surface (`SetVolumeHDFCompression`,
  `SetFaultHDFCompression`, `SetFaultParamsBP5`, `SetVolumeSaveEnabled`,
  `SetTotalRunTime`, `GetSchedule()`, `output_every_n_steps`, `fixed_dt`,
  `PeekShouldWrite`, `CommitSchedule`).
- **TPV* and BP5 drivers are read-only.**  Only `drivers/spatial_dyn_driver.cpp`
  changes.  The PV regression gates from §"Files this plan does NOT touch"
  in `spatial_dynamic_rupture_plan.md` apply unchanged.
- **TOML is the source of truth; CLI overrides win.**  Adding a CLI flag MUST
  also add a TOML field with the same default, mirrored by an entry in
  `safs/project_7.0_alternative/document/spatial_friction_config_schema.md`.
- **Default is OFF for volume / bulk** (matches BP5 + R-313 in `CLAUDE.md`):
  spatial-driver defaults should be `paraview_volume = "off"`,
  `paraview_bulk = "off"`, `paraview_fault = "hdf5"`.  The current schema
  defaults (`hdf5` for all three) are too aggressive for SAFS production
  meshes (~100M elements × 0.001s fault dt → multi-TB output).  This is a
  **default-flip**, explicitly called out in the schema doc.
- **All physical units SI**; tfinal in TOML accepts `"12s"` suffix only
  (matches the existing schema in `spatial_friction.cpp`'s `toml_time_seconds`).
- **No `safs_*` filenames**.  All work lives under
  `drivers/spatial_dyn_driver.cpp` and `spatial/code/spatial_friction.{hpp,cpp}`,
  per the naming convention in `spatial_dynamic_rupture_plan.md:728`.
- **Spatial-friction static fields are heterogeneous per DOF** (unlike TPV104),
  so the fault-static-parameter array set published via `SetFaultParamsBP5`
  MUST be a SAFS-specific overload that publishes `lsw_mu_s`, `lsw_mu_d`,
  `lsw_d_c`, `T_forced_rupture`, `sigma_n_corr` (initial), `tau1_corr` (initial),
  `tau2_corr` (initial) — NOT the BP5 `a, Dc, x2, x3`.

## Phase 1: CLI flag surface + build guards

### Goal
After Phase 1, `seas_spatial_dyn_driver --help` accepts every paraview flag
the TPV104 sbatch uses (and the 12 additional knobs from the TPV104 driver),
each guarded against MFEM build-flag absence.  No on-disk output change yet.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp` — extend CLI parser at lines 497-506 and
  the wiring block at lines 1045-1108.

### Detailed Requirements

1. **Add CLI parsing for the missing flags.**  Insert immediately after line
   506 (after the existing `--paraview-max-snapshots` line):
   ```cpp
   const bool   cli_pv                       = HasFlag(argc, argv, "--paraview");
   const real_t cli_pv_dt                    = GetRealArg(argc, argv, "--paraview-dt", -1.0);
   const int    cli_pv_every                 = GetIntArg (argc, argv, "--paraview-every", 0);
   const bool   cli_pv_force_vtu             = HasFlag(argc, argv, "--paraview-fault-vtu");
   const bool   cli_pv_force_hdf5            = HasFlag(argc, argv, "--paraview-fault-hdf5");
   const bool   cli_pv_legacy_ascii          = HasFlag(argc, argv, "--paraview-fault-legacy-ascii");
   const int    cli_pv_fault_deflate         = GetIntArg (argc, argv, "--paraview-fault-deflate-level", -1);
   const int    cli_pv_bulk_deflate          = GetIntArg (argc, argv, "--paraview-bulk-deflate-level",  -1);
   const bool   cli_pv_volume_force_vtu      = HasFlag(argc, argv, "--paraview-volume-vtu");
   const bool   cli_pv_volume_force_hdf5     = HasFlag(argc, argv, "--paraview-volume-hdf5");
   const int    cli_pv_volume_deflate        = GetIntArg (argc, argv, "--paraview-volume-deflate-level", -1);
   const real_t cli_pv_coseismic_dt          = GetRealArg(argc, argv, "--paraview-coseismic-dt",   -1.0);
   const real_t cli_pv_nucleation_dt         = GetRealArg(argc, argv, "--paraview-nucleation-dt",  -1.0);
   const real_t cli_pv_interseismic_dt       = GetRealArg(argc, argv, "--paraview-interseismic-dt",-1.0);
   ```

2. **Add config-merge logic** at lines 571-579 (after existing `cfg.output.*`
   overrides), mapping each new CLI value into the corresponding new
   `cfg.output.*` field added in Phase 2:
   ```cpp
   if (cli_pv)                            { cfg.output.paraview_enabled       = true; }
   if (cli_pv_dt              > 0.0)      { cfg.output.paraview_fault_dt      = cli_pv_dt; }
   if (cli_pv_every           > 0)        { cfg.output.paraview_every_steps   = cli_pv_every; }
   if (cli_pv_fault_deflate  >= 0)        { cfg.output.paraview_fault_deflate_level   = cli_pv_fault_deflate; }
   if (cli_pv_bulk_deflate   >= 0)        { cfg.output.paraview_bulk_deflate_level    = cli_pv_bulk_deflate;  }
   if (cli_pv_volume_deflate >= 0)        { cfg.output.paraview_volume_deflate_level  = cli_pv_volume_deflate; }
   if (cli_pv_coseismic_dt    > 0.0)      { cfg.output.paraview_coseismic_dt   = cli_pv_coseismic_dt; }
   if (cli_pv_nucleation_dt   > 0.0)      { cfg.output.paraview_nucleation_dt  = cli_pv_nucleation_dt; }
   if (cli_pv_interseismic_dt > 0.0)      { cfg.output.paraview_interseismic_dt = cli_pv_interseismic_dt; }
   // Mode-string flags override the [output].paraview_{fault,volume} TOML.
   if (cli_pv_force_vtu)                  { cfg.output.paraview_fault  = "vtu";  }
   if (cli_pv_force_hdf5)                 { cfg.output.paraview_fault  = "hdf5"; }
   if (cli_pv_volume_force_vtu)           { cfg.output.paraview_volume = "vtu";  }
   if (cli_pv_volume_force_hdf5)          { cfg.output.paraview_volume = "hdf5"; }
   if (cli_pv_legacy_ascii)               { cfg.output.paraview_fault_legacy_ascii = true;
                                            cfg.output.paraview_fault = "vtu"; }
   ```

3. **Add 6 build-flag-guard `MFEM_ABORT`s** before the `pv_out` construction
   block, mirroring `tpv104_driver.cpp:666-757`:
   - `--paraview-fault-vtu` + `--paraview-fault-hdf5` → abort (incompatible).
   - `--paraview-fault-legacy-ascii` + `--paraview-fault-hdf5` → abort.
   - `--paraview-{fault,bulk,volume}-zfp-tol` + matching deflate-level → abort.
   - `--paraview-volume-vtu` + `--paraview-volume-hdf5` → abort.
   - On non-`MFEM_USE_HDF5` builds: any `*-hdf5` opt-in or any `*-zfp-tol > 0`
     or any `*-deflate-level >= 0` → abort with build instruction.
   - On non-`MFEM_USE_H5Z_ZFP` builds: any `*-zfp-tol > 0` → abort.

   These guards MUST fire **before** any `pv_out` constructor invocation so a
   bad flag combination fails at parse time, not after the run has spun up.

4. **Honor `cfg.output.paraview_enabled` as the master gate.**  Set
   `volume_pv_enabled` / `fault_pv_enabled` / `bulk_pv_enabled` to `false` when
   `cfg.output.paraview_enabled == false`, regardless of the per-collection
   mode strings.  This lets the SAFS production default (Phase 3 below) be
   `paraview_enabled = false` even though `paraview_fault = "hdf5"` is in the
   schema for documentation purposes.

### Interfaces
- No new public/private member function signatures introduced; all changes are
  local to the spatial driver and the `OutputSpec` struct (Phase 2).

### Edge Cases to Handle
- `--paraview` alone (no `--paraview-dt`, no `--paraview-bulk-dt`) → enables
  fault output at the TOML `paraview_fault_dt` cadence (default `0.001 s`);
  volume + bulk stay off unless their per-collection mode was set explicitly.
- `--paraview-volume-zfp-tol 1e-3` on a `MFEM_USE_HDF5=NO` build → abort at
  parse time, not at first ForceSave.
- `--paraview-coseismic-dt 0.01 --paraview-interseismic-dt 1.0` without
  `--paraview-bulk-dt`/`--paraview-volume-dt` → applies only to the primary
  `pv_out` schedule; secondary `pv_bulk_out` (if Phase 4 enables it) uses its
  own fixed_dt.  Tested in P-1 below.
- `paraview_enabled = false` AND `cli_pv_dt > 0` → CLI wins, `paraview_enabled`
  flips to `true` (per item #2 above).
- `--paraview-fault-legacy-ascii` on a `MFEM_USE_HDF5=YES` build but with
  `paraview_fault = "hdf5"` in TOML → CLI override should switch the mode to
  `"vtu"` AND set the legacy flag (per item #2).

### Acceptance Criteria
- [ ] `seas_spatial_dyn_driver --paraview --paraview-dt 0.5` exits the parse
      phase without warnings on a build with `MFEM_USE_HDF5=YES`.
- [ ] Same invocation on `MFEM_USE_HDF5=NO` aborts at parse time with the
      message `"--paraview-fault-hdf5 requires the seas-mfem build to ..."`
      (the abort fires because the default `cfg.output.paraview_fault = "hdf5"`
      survives in the schema; flip-to-vtu is the user's job).
- [ ] `--paraview-fault-vtu --paraview-fault-hdf5` aborts with the exact
      message `"--paraview-fault-vtu and --paraview-fault-hdf5 are ..."`.
- [ ] Existing TPV104/205 byte-exact regression suite remains green
      (`make test-tpv104-checkpoint`, `make test-tpv205-byteexact`).

### Dependencies
- Depends on: nothing (all new code is in `spatial_dyn_driver.cpp`).
- Required by: Phase 2 (config-schema extension), Phase 4 (secondary `pv_bulk_out`
  wiring).

---

## Phase 2: OutputSpec schema extension + TOML parser + docs

### Goal
After Phase 2, the `OutputSpec` struct in `spatial_friction.hpp` holds every
new field referenced by Phase 1's CLI merge logic, the TOML parser reads them
from `[output]`, the validator rejects out-of-range values, and the schema doc
documents the new fields.

### Files to Modify
- `spatial/code/spatial_friction.hpp` — extend `struct OutputSpec`.
- `spatial/code/spatial_friction.cpp` — extend `toml::*` parser block at
  lines 655-696.
- `safs/project_7.0_alternative/document/spatial_friction_config_schema.md` —
  document new keys + the flipped defaults.

### Detailed Requirements

1. **Extend `OutputSpec`** (insert after `int max_snapshots = 5000;` at
   `spatial_friction.hpp:120`):
   ```cpp
   bool        paraview_enabled              = false;  // master gate
   int         paraview_every_steps          = 0;      // 0 = use *_dt
   bool        paraview_fault_legacy_ascii   = false;  // debug only
   int         paraview_volume_deflate_level = -1;     // -1 = none, [0..9]
   int         paraview_bulk_deflate_level   = -1;
   int         paraview_fault_deflate_level  = -1;
   real_t      paraview_coseismic_dt         = -1.0;   // -1 = unset (use AdaptiveSchedule default)
   real_t      paraview_nucleation_dt        = -1.0;
   real_t      paraview_interseismic_dt      = -1.0;
   ```
   **Flip the defaults** for the three per-collection mode strings (replacing
   the existing default of `"hdf5"` for all three):
   ```cpp
   std::string paraview_volume        = "off";   // was "hdf5"
   std::string paraview_bulk          = "off";   // was "hdf5"
   std::string paraview_fault         = "hdf5";  // unchanged
   ```

2. **Extend the TOML parser** (insert after the existing `max_snapshots`
   parse at `spatial_friction.cpp:664`):
   ```cpp
   cfg.output.paraview_enabled              = toml_bool(o, "paraview_enabled",              false);
   cfg.output.paraview_every_steps          = toml_int (o, "paraview_every_steps",          0);
   cfg.output.paraview_fault_legacy_ascii   = toml_bool(o, "paraview_fault_legacy_ascii",   false);
   cfg.output.paraview_volume_deflate_level = toml_int (o, "paraview_volume_deflate_level", -1);
   cfg.output.paraview_bulk_deflate_level   = toml_int (o, "paraview_bulk_deflate_level",   -1);
   cfg.output.paraview_fault_deflate_level  = toml_int (o, "paraview_fault_deflate_level",  -1);
   cfg.output.paraview_coseismic_dt         = toml_time_seconds(o, "paraview_coseismic_dt",    -1.0);
   cfg.output.paraview_nucleation_dt        = toml_time_seconds(o, "paraview_nucleation_dt",   -1.0);
   cfg.output.paraview_interseismic_dt      = toml_time_seconds(o, "paraview_interseismic_dt", -1.0);
   ```

3. **Extend the validator** (insert after the existing
   `max_snapshots >= 1` check at `spatial_friction.cpp:695`):
   ```cpp
   MFEM_VERIFY(cfg.output.paraview_every_steps >= 0,
               "[output].paraview_every_steps must be >= 0; got "
               << cfg.output.paraview_every_steps);
   for (auto kv : { std::pair{"paraview_volume_deflate_level", cfg.output.paraview_volume_deflate_level},
                    std::pair{"paraview_bulk_deflate_level",   cfg.output.paraview_bulk_deflate_level},
                    std::pair{"paraview_fault_deflate_level",  cfg.output.paraview_fault_deflate_level} })
   {
      MFEM_VERIFY(kv.second == -1 || (kv.second >= 0 && kv.second <= 9),
                  "[output]." << kv.first << " must be -1 (none) or in [0..9]; got "
                  << kv.second);
   }
   for (auto kv : { std::pair{"paraview_coseismic_dt",    cfg.output.paraview_coseismic_dt},
                    std::pair{"paraview_nucleation_dt",   cfg.output.paraview_nucleation_dt},
                    std::pair{"paraview_interseismic_dt", cfg.output.paraview_interseismic_dt} })
   {
      MFEM_VERIFY(kv.second < 0.0 || kv.second > 0.0,
                  "[output]." << kv.first << " must be unset (negative) or > 0; got "
                  << kv.second);
   }
   ```

4. **Document new keys in
   `safs/project_7.0_alternative/document/spatial_friction_config_schema.md`.**
   Insert the 9 new rows into the existing `[output]` table at line 174 (after
   `max_snapshots`).  Update the **default** for `paraview_volume` and
   `paraview_bulk` to `"off"` and add a footnote: "SAFS production default is
   OFF for volume+bulk because the SAFS mesh (~100M elements × 0.001 s fault
   cadence) would emit multi-TB output if all three streams were on; opt in
   per collection."

### Interfaces
- `OutputSpec` gains 9 new public fields.  No method signatures change.
- TOML parser is purely additive; all new fields have defaults so existing
  `.toml` files continue to parse byte-equivalently (verified by Phase 5
  test T-25).

### Edge Cases to Handle
- A user TOML with `paraview_volume = "off"` but `paraview_volume_zfp_tol = 1e-3` →
  validator accepts (the ZFP tol is dormant); driver emits a one-time rank-0
  info line "paraview_volume = off; paraview_volume_zfp_tol ignored".
- TOML with `paraview_enabled = false` but CLI `--paraview` → CLI override
  flips to true (Phase 1 item #2); validator never sees the conflict.
- `paraview_every_steps = 5 AND paraview_fault_dt = 0.001` → CLI/TOML
  priority must be defined: `paraview_every_steps` wins (matches TPV104
  `tpv104_driver.cpp:2052`).  Documented in the schema doc.

### Acceptance Criteria
- [ ] Existing `test_spatial_friction_config` suite (28 tests) still passes
      (additive change, no default-string churn outside the volume/bulk flip).
- [ ] New tests T-24 (parser accepts the 9 new fields), T-25 (defaults match
      Phase 1 contract), T-26 (validator rejects deflate=10), T-27 (validator
      rejects zero-or-negative cadence), T-28 (default flip: parsing
      `paraview_volume` without a value returns `"off"`) all pass.

### Dependencies
- Depends on: Phase 1 (the CLI-merge logic references these fields).
- Required by: Phase 3, Phase 4.

---

## Phase 3: Default-off enforcement + master enable gate

### Goal
After Phase 3, a spatial-driver run with no `--paraview*` flags and a TOML
that does not set `[output].paraview_enabled = true` writes ZERO ParaView
files (matches the BP5 production posture per `CLAUDE.md` R-313).

### Files to Modify
- `drivers/spatial_dyn_driver.cpp` — modify the `pv_out` construction
  guard at lines 1073-1108.

### Detailed Requirements

1. **Replace the current `volume_pv_enabled || fault_pv_enabled` guard**
   (line 1073) with an additional check on the master gate:
   ```cpp
   const bool any_pv_requested =
      cfg.output.paraview_enabled
      || volume_pv_enabled
      || fault_pv_enabled
      || bulk_pv_enabled;
   if (any_pv_requested && (volume_pv_enabled || fault_pv_enabled))
   {
      // ... existing pv_out construction
   }
   ```
   The `any_pv_requested` shortcut keeps the construction skipped when ALL
   three modes are `"off"` AND the master gate is `false`.  The
   `volume_pv_enabled || fault_pv_enabled` half ensures we still construct
   `pv_out` when either of those is opted in (the bulk-only path goes to
   `pv_bulk_out` in Phase 4).

2. **Emit a clean rank-0 banner** when `pv_out` IS constructed, mirroring
   `tpv104_driver.cpp:2117-2148`:
   ```cpp
   if (rank == 0)
   {
      std::cout << "ParaView output: ON (prefix=" << cfg.output.output_dir
                << ")\n";
      std::cout << "  Volume:   " << cfg.output.paraview_volume
                << " (dt=" << cfg.output.paraview_volume_dt << " s)\n";
      std::cout << "  Bulk:     " << cfg.output.paraview_bulk
                << " (dt=" << cfg.output.paraview_bulk_dt   << " s)\n";
      std::cout << "  Fault:    " << cfg.output.paraview_fault
                << " (dt=" << cfg.output.paraview_fault_dt  << " s)\n";
      if (cfg.output.paraview_volume_zfp_tol > 0.0)
      { std::cout << "  Volume ZFP tol: " << cfg.output.paraview_volume_zfp_tol << "\n"; }
      // ... bulk, fault similar
      if (cfg.output.paraview_coseismic_dt > 0.0
          || cfg.output.paraview_nucleation_dt > 0.0
          || cfg.output.paraview_interseismic_dt > 0.0)
      { std::cout << "  Regime-adaptive cadence active.\n"; }
   }
   ```

3. **Symmetric quiet exit when PV is fully off** — emit a one-line rank-0
   note `"ParaView output: OFF"`.  Prevents user confusion about a silent
   driver.

### Interfaces
- No new functions.  Behavior change only.

### Edge Cases to Handle
- TOML with `paraview_enabled = true` but every per-collection mode is
  `"off"` → emit a rank-0 warning `"WARNING: [output].paraview_enabled = true
  but every per-collection mode is 'off' — no ParaView files will be
  written."` and proceed.  Do not abort (user might be sweeping configs).
- TOML with `paraview_enabled = false` AND no CLI flags → no `pv_out`
  constructed, no files in `<output_dir>/ParaView*`, banner says OFF.

### Acceptance Criteria
- [ ] Smoke test: `seas_spatial_dyn_driver --config <minimal>.toml --tfinal 0.1s`
      with no `paraview_enabled` and no `--paraview*` produces ZERO files
      matching `<output_dir>/*.vtkhdf` or `<output_dir>/ParaView/`.
- [ ] Smoke test: same config + `--paraview --paraview-dt 0.05` produces
      exactly one `<output_dir>/volume.vtkhdf` (fault projections only since
      `paraview_volume = "off"` is the new default).

### Dependencies
- Depends on: Phase 2.
- Required by: Phase 4 (secondary collection only constructed when the master
  gate is on AND `paraview_bulk != "off"`).

---

## Phase 4: Secondary `pv_bulk_out` collection (6 stress components)

### Goal
After Phase 4, when `paraview_bulk != "off"`, the driver constructs a
secondary `seas::ParaViewOutput<ParMesh>` named `"stress"` under
`<output_dir>/ParaView_bulk/`, registers six L2 stress component GFs
(σ_xx, σ_yy, σ_zz, σ_xy, σ_xz, σ_yz), and writes them at the
`paraview_bulk_dt` cadence with `paraview_bulk_zfp_tol` ZFP compression.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp` — extend the section "17. ParaView
  output wiring" (currently lines 1045-1108) and "paraview_write" lambda
  (currently lines 1199-1239).

### Detailed Requirements

1. **Declare bulk-collection state** alongside the existing `pv_out`
   declaration (line 1071):
   ```cpp
   std::unique_ptr<seas::ParaViewOutput<ParMesh>> pv_bulk_out;
   std::unique_ptr<L2_FECollection>               pv_bulk_sigma_fec;
   std::unique_ptr<typename seas::GFType<ParMesh>::FESType> pv_bulk_sigma_fes;
   std::unique_ptr<typename seas::GFType<ParMesh>::type>
      pv_bulk_sxx_gf, pv_bulk_syy_gf, pv_bulk_szz_gf,
      pv_bulk_sxy_gf, pv_bulk_sxz_gf, pv_bulk_syz_gf;
   ```

2. **Construct the secondary collection** under
   `<output_dir>/ParaView_bulk/` with collection name `"stress"` when
   `bulk_pv_enabled && cfg.output.paraview_enabled` (mirroring
   `tpv104_driver.cpp:2191-2228`):
   ```cpp
   if (bulk_pv_enabled && cfg.output.paraview_enabled)
   {
      const std::string bulk_dir = cfg.output.output_dir + "/ParaView_bulk";
      if (rank == 0) { std::filesystem::create_directories(bulk_dir); }
   #ifdef MFEM_USE_MPI
      MPI_Barrier(comm);
   #endif
      auto bulk_mode = ParseVolumeMode(cfg.output.paraview_bulk, /*&out*/bulk_pv_enabled);
      pv_bulk_out = std::make_unique<seas::ParaViewOutput<ParMesh>>(
         bulk_dir, pmesh, cfg.mesh.order, "stress", bulk_mode);

      pv_bulk_sigma_fec = std::make_unique<L2_FECollection>(
         cfg.mesh.order, 3, BasisType::GaussLobatto);
      pv_bulk_sigma_fes = std::make_unique<...FESType>(&pmesh, pv_bulk_sigma_fec.get());
      pv_bulk_sxx_gf = std::make_unique<...GF>(pv_bulk_sigma_fes.get());
      // ... 5 more
      *pv_bulk_sxx_gf = 0.0;  // ... 5 more
      pv_bulk_out->RegisterDomainField("sigma_xx", pv_bulk_sxx_gf.get());
      // ... 5 more

      pv_bulk_out->fixed_dt = cfg.output.paraview_bulk_dt;

   #ifdef MFEM_USE_HDF5
      if (cfg.output.paraview_bulk_zfp_tol > 0.0
          && bulk_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_bulk_out->SetVolumeHDFCompression(
            ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            cfg.output.paraview_bulk_zfp_tol);
      }
      else if (cfg.output.paraview_bulk_deflate_level >= 0
               && bulk_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_bulk_out->SetVolumeHDFCompression(
            ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(cfg.output.paraview_bulk_deflate_level));
      }
   #endif
   }
   ```

3. **Extend the `paraview_write` lambda** at line 1199 to publish the
   sigma fields to `pv_bulk_out` when the lambda is called and
   `pv_bulk_out` is alive:
   ```cpp
   auto paraview_write = [&](int step_num, real_t time, real_t V_max)
   {
      const bool fault_wants = pv_out      && pv_out     ->PeekShouldWrite(step_num, time, V_max);
      const bool bulk_wants  = pv_bulk_out && pv_bulk_out->PeekShouldWrite(step_num, time, V_max);
      if (!fault_wants && !bulk_wants) { return; }

      if (bulk_wants)
      {
         // Q ordering from dynamic/wave_state.hpp:
         //   SXX=0, SYY=1, SZZ=2, SXY=3, SYZ=4, SXZ=5.
         std::memcpy(pv_bulk_sxx_gf->GetData(), Q.GetData() + SXX * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_syy_gf->GetData(), Q.GetData() + SYY * ndof_total,
                     ndof_total * sizeof(real_t));
         // ... SZZ, SXY, SYZ, SXZ
         pv_bulk_out->ForceSave(step_num, time);
      }
      if (!fault_wants) { return; }
      // ... existing pv_local_slip / pv_local_state population
      // ... existing pv_out->UpdateFaultFieldsBP5 + ForceSave + CommitSchedule
   };
   ```

### Interfaces
- No new functions exposed.  All state is local to the driver.

### Edge Cases to Handle
- `paraview_bulk = "vtu"` AND `paraview_bulk_zfp_tol > 0` → the ZFP setter is
  guarded by `bulk_mode == Hdf5`, so it's a no-op; emit a one-time rank-0
  warning (matches Phase 1 build-guard policy).
- `paraview_bulk_dt = 0` would be rejected by the validator (Phase 2 item #3).
- `pv_bulk_out` constructed but `bulk_wants == false` at every step (cadence
  too coarse) → no files written, no abort.
- Restart: the secondary `pv_bulk_out` collection re-initialises in its
  default state at restart time, emits a frame at `t = restart_t`, and starts
  its schedule fresh.  This matches TPV104's known limitation
  (`tpv104_driver.cpp:2403-2412`).  Plan documents this in the spatial
  driver's restart section.

### Acceptance Criteria
- [ ] Smoke test: `seas_spatial_dyn_driver --config <bp5-style>.toml
      --paraview --paraview-bulk-dt 0.5 --paraview-bulk-zfp-tol 1e-3 --tfinal 1.0s`
      on a `MFEM_USE_HDF5=YES,MFEM_USE_H5Z_ZFP=YES` build produces exactly
      one `<output_dir>/ParaView_bulk/stress.vtkhdf` file containing the 6
      sigma components at every committed cycle.
- [ ] On `MFEM_USE_HDF5=NO` build, same flags abort at Phase 1's build guard
      BEFORE the secondary collection is constructed.

### Dependencies
- Depends on: Phase 1 (build guards), Phase 2 (config fields), Phase 3
  (master gate).
- Required by: nothing.

---

## Phase 5: Volume domain field registration (velocity + mpi_rank)

### Goal
After Phase 5, when `paraview_volume != "off"`, the primary `pv_out`
collection contains a `velocity` vector field (3 components, L2 p=order) and
an `mpi_rank` scalar field (L2 p=0).  ParaView users can see the actual
wavefield, not just the fault projections.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp` — extend the `pv_out` construction
  block (Phase 1/3 location) and the `paraview_write` lambda (Phase 4
  extension).

### Detailed Requirements

1. **Declare and register the volume GFs** inside the existing
   `if (volume_pv_enabled || fault_pv_enabled) { ... }` block, conditionally
   on `volume_pv_enabled`:
   ```cpp
   std::unique_ptr<L2_FECollection> pv_vel_fec, pv_rank_fec;
   std::unique_ptr<typename seas::GFType<ParMesh>::FESType> pv_vel_fes, pv_rank_fes;
   std::unique_ptr<typename seas::GFType<ParMesh>::type>    pv_vel_gf,  pv_rank_gf;

   if (volume_pv_enabled)
   {
      pv_vel_fec  = std::make_unique<L2_FECollection>(
         cfg.mesh.order, 3, BasisType::GaussLobatto);
      pv_vel_fes  = std::make_unique<...FESType>(
         &pmesh, pv_vel_fec.get(), 3, Ordering::byNODES);
      pv_vel_gf   = std::make_unique<...GF>(pv_vel_fes.get());
      *pv_vel_gf  = 0.0;
      pv_out->RegisterDomainField("velocity", pv_vel_gf.get());

      pv_rank_fec = std::make_unique<L2_FECollection>(0, 3);
      pv_rank_fes = std::make_unique<...FESType>(&pmesh, pv_rank_fec.get());
      pv_rank_gf  = std::make_unique<...GF>(pv_rank_fes.get());
      *pv_rank_gf = static_cast<real_t>(rank);
      pv_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());
   }
   ```

2. **Update the `paraview_write` lambda** to memcpy the velocity slice
   from `Q` into `pv_vel_gf` BEFORE `pv_out->ForceSave(step_num, time)` —
   immediately after the `bulk_wants` branch from Phase 4:
   ```cpp
   if (fault_wants && pv_out->GetVolumeSaveEnabled() && pv_vel_gf)
   {
      std::memcpy(pv_vel_gf->GetData(),
                  Q.GetData() + VX * ndof_total,
                  3 * ndof_total * sizeof(real_t));
   }
   ```
   The `VX` / `VY` / `VZ` slots are CONTIGUOUS in the byNODES ordering, so
   a single memcpy of `3 * ndof_total` doubles is correct (per
   `dynamic/wave_state.hpp:35-37` and the TPV104 reference at
   `tpv104_driver.cpp:2253-2257`).

### Interfaces
- No new functions.  `pv_out->RegisterDomainField(...)` is the existing
  paraview-compaction API.

### Edge Cases to Handle
- `paraview_volume = "off"` AND `paraview_fault = "hdf5"` → `pv_vel_gf` is
  null, the lambda skips the memcpy, ForceSave writes the fault projections
  only.  No-op on the volume side.
- `Q` re-allocated mid-loop (it isn't in the current driver, but the lambda
  captures `Q.GetData()` by `&Q`, so a reallocation would invalidate the
  pointer — the lambda re-fetches `Q.GetData()` on every call, so this is
  safe by construction).
- Order-0 mesh (`cfg.mesh.order == 0`) → `L2_FECollection(0, 3)` is valid
  (one DOF per element), velocity becomes piecewise constant.  Acceptable
  for visualization.

### Acceptance Criteria
- [ ] Smoke test: same as Phase 4, but with `--paraview-volume hdf5
      --paraview-volume-zfp-tol 1e-3`.  The primary `volume.vtkhdf` contains
      `velocity` and `mpi_rank` fields readable by ParaView 5.11+.
- [ ] When `paraview_volume = "off"`, the primary `volume.vtkhdf` contains
      ONLY the BP5 fault projection fields (no `velocity`, no `mpi_rank`) —
      preserving the round-7 file size.

### Dependencies
- Depends on: Phase 3 (the gate that decides when `pv_out` is constructed).
- Required by: nothing.

---

## Phase 6: SAFS-aware fault static-parameter publishing

### Goal
After Phase 6, the fault collection includes per-DOF heterogeneous static
arrays: `lsw_mu_s`, `lsw_mu_d`, `lsw_d_c`, `T_forced_rupture`,
`sigma_n_init`, `tau1_init`, `tau2_init`.  This replaces TPV104's
`SetFaultParamsBP5(a, Dc, x2, x3)` (which assumes homogeneous a, Dc) with a
SAFS-specific call that publishes the spatial heterogeneity the driver
already has in `dof_data`.

### Files to Create
- None.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp` — extend the `pv_out` setup block.
- `io/paraview_output.hpp` — **READ-ONLY exception**: the existing
  `SetFaultParamsBP5` API publishes exactly 4 named arrays
  (`a`, `Dc`, `x2`, `x3`).  Adding a SAFS-specific overload requires either:
   (a) extending the existing function to take an `std::map<std::string, const
       Vector*>` of arbitrary named arrays, or
   (b) adding a new method `SetFaultParamsSpatial(const std::vector<
       std::pair<std::string, const Vector*>>& named_arrays)`.
  Option (b) is preferred (additive, no signature change to the existing
  TPV104/205 call sites).  This is the ONLY io/paraview_output.hpp change in
  this plan and MUST be guarded by a "TPV104/205 byte-exact regression
  still passes" gate.

### Detailed Requirements

1. **Add `SetFaultParamsSpatial` to `seas::ParaViewOutput`** with this
   signature:
   ```cpp
   /// @brief SAFS analog of SetFaultParamsBP5 — publishes named static
   /// per-DOF arrays as fault-surface fields.  Each Vector MUST have
   /// size == num_fault_total (the same sizing rule as the dynamic
   /// fields published via UpdateFaultFieldsBP5).
   void SetFaultParamsSpatial(
      const std::vector<std::pair<std::string, const Vector*>>& named_arrays);
   ```
   The body stores the arrays in a new member `fault_params_spatial_` and
   writes them to the fault VTU/VTKHDF in `WriteFaultSurfaceVTU` (and the
   VTKHDF mirror) the same way the existing `SetFaultParamsBP5` arrays are
   written.

2. **Spatial driver pushes the 7 SAFS arrays** at `pv_out` setup time
   (insert after the existing `pv_out->InitFaultOutputBP5(...)` call at
   `spatial_dyn_driver.cpp:1106`):
   ```cpp
   Vector pv_lsw_mu_s   (num_fault_total),
          pv_lsw_mu_d   (num_fault_total),
          pv_lsw_d_c    (num_fault_total),
          pv_T_forced   (num_fault_total),
          pv_sig_n_init (num_fault_total),
          pv_tau1_init  (num_fault_total),
          pv_tau2_init  (num_fault_total);
   for (int i = 0; i < num_fault_total; ++i)
   {
      const DOFData &d = dof_data[i];
      pv_lsw_mu_s  (i) = d.lsw_mu_s;
      pv_lsw_mu_d  (i) = d.lsw_mu_d;
      pv_lsw_d_c   (i) = d.lsw_d_c;
      pv_T_forced  (i) = d.T_forced_rupture;
      pv_sig_n_init(i) = d.sigma_n_corr;   // captured at init time
      pv_tau1_init (i) = d.tau1_corr;
      pv_tau2_init (i) = d.tau2_corr;
   }
   pv_out->SetFaultParamsSpatial({
      {"lsw_mu_s",         &pv_lsw_mu_s},
      {"lsw_mu_d",         &pv_lsw_mu_d},
      {"lsw_d_c",          &pv_lsw_d_c},
      {"T_forced_rupture", &pv_T_forced},
      {"sigma_n_init",     &pv_sig_n_init},
      {"tau1_init",        &pv_tau1_init},
      {"tau2_init",        &pv_tau2_init},
   });
   ```
   The 7 arrays MUST live for the lifetime of `pv_out` — store them in the
   enclosing scope, not inside the construction block.

3. **Set `pv_out->SetTotalRunTime(cfg.time.tfinal)`** so the AdaptiveSchedule
   has a horizon (mirroring `tpv104_driver.cpp:2113`).  Insert immediately
   after the `GetSchedule().max_total_snapshots = ...` line at
   `spatial_dyn_driver.cpp:1082`:
   ```cpp
   pv_out->SetTotalRunTime(cfg.time.tfinal);
   if (cfg.output.paraview_coseismic_dt    > 0.0)
   { pv_out->GetSchedule().dt_coseismic    = cfg.output.paraview_coseismic_dt;    }
   if (cfg.output.paraview_nucleation_dt   > 0.0)
   { pv_out->GetSchedule().dt_nucleation   = cfg.output.paraview_nucleation_dt;   }
   if (cfg.output.paraview_interseismic_dt > 0.0)
   { pv_out->GetSchedule().dt_interseismic = cfg.output.paraview_interseismic_dt; }
   pv_out->GetSchedule().Validate();
   ```
   The regime-adaptive schedule is dormant by default (`dt_coseismic = 0.01 s`,
   `dt_nucleation = 1 s`, `dt_interseismic = 1 yr`) — for a 12-s dynamic
   rupture, the default `dt_coseismic = 0.01 s` is the operative cadence
   once V_max exceeds the regime threshold.  The new TOML/CLI fields let
   the user tune this without recompiling.

### Interfaces
- New method `seas::ParaViewOutput::SetFaultParamsSpatial(...)` (additive).
- All other interfaces unchanged.

### Edge Cases to Handle
- A SAFS DOF with `T_forced_rupture = 1e9` (the "never forced" sentinel) —
  the array publishes the sentinel value; ParaView users see `1e9` as a clear
  "outside nucleation zone" marker.
- `num_fault_total == 0` on a rank with no fault DOFs — the 7 arrays are
  empty Vectors; `SetFaultParamsSpatial` accepts zero-size and is a no-op
  on that rank.
- Restart: the static arrays are re-published from the post-init `dof_data`
  on every run.  Since they don't evolve in time, this is correct.

### Acceptance Criteria
- [ ] TPV104 + TPV205 byte-exact regression gates still pass (the new method
      is additive; existing `SetFaultParamsBP5` call sites untouched).
- [ ] Smoke test: `seas_spatial_dyn_driver --paraview --paraview-fault hdf5`
      produces a `fault.vtkhdf` with 5 dynamic fields (slip, slip_rate,
      traction, mu_eff, sigma_n) + 7 static fields (lsw_mu_s, lsw_mu_d,
      lsw_d_c, T_forced_rupture, sigma_n_init, tau1_init, tau2_init).
- [ ] On a heterogeneous TOML (e.g. spatial rule sets `lsw_mu_s = 1.2` inside
      a depth band), the fault VTKHDF shows the depth-band step in the
      `lsw_mu_s` field.

### Dependencies
- Depends on: Phase 1-5.
- Required by: nothing.

---

## Testing Strategy

Each phase has acceptance criteria as ad-hoc smoke tests.  In addition, the
following unit tests are required:

### Phase 1 — CLI + build guards
- **T-30** (new file `tests/unit/test_spatial_paraview_cli.cpp`):
  Drive `argv = {"prog", "--paraview", "--paraview-dt", "0.5"}` through
  the spatial driver's parser (extracted into a small testable helper if
  necessary), assert `cfg.output.paraview_enabled == true` and
  `cfg.output.paraview_fault_dt == 0.5`.
- **T-31**: Same, with `argv = {..., "--paraview-fault-vtu",
  "--paraview-fault-hdf5"}` — expect `MFEM_ABORT` (use the existing
  `RunInChild` pattern from `test_spatial_friction_resolver.cpp`).
- **T-32**: Same, but on a build with `#ifndef MFEM_USE_HDF5` mocked via
  a compile-time guard — expect abort on `--paraview-fault-hdf5`.

### Phase 2 — TOML parser
- **T-33** (extend `test_spatial_friction_config.cpp`): parse a TOML with
  all 9 new `[output]` fields set explicitly, assert each round-trips.
- **T-34**: parse a TOML omitting all new fields, assert defaults match
  Phase 2 item #1 verbatim.
- **T-35**: parse a TOML with `paraview_volume_deflate_level = 10`, expect
  `MFEM_VERIFY` abort.
- **T-36**: parse a TOML with `paraview_coseismic_dt = 0`, expect
  `MFEM_VERIFY` abort.
- **T-37**: parse a TOML with `paraview_volume` omitted, assert
  `cfg.output.paraview_volume == "off"` (the new default).

### Phase 3 — Default-off
- **T-38** (new file `tests/integration/test_spatial_paraview_default_off.cpp`):
  Run the driver with `--config <minimal>.toml --tfinal 0.1s`, no
  paraview flags.  Assert ZERO files matching `<output_dir>/*.vtkhdf` or
  `<output_dir>/ParaView*/`.

### Phase 4 — Secondary collection
- **T-39**: Run with `--paraview --paraview-bulk hdf5 --paraview-bulk-dt 0.5
  --tfinal 1.0s`.  Assert `<output_dir>/ParaView_bulk/stress.vtkhdf` exists
  and is non-empty.  Assert it contains all 6 sigma_* fields (use the
  existing `test_vtkhdf_zfp.cpp` pattern for HDF5 field enumeration).

### Phase 5 — Volume velocity registration
- **T-40**: Run with `--paraview --paraview-volume hdf5 --paraview-volume-dt
  0.5 --tfinal 1.0s`.  Assert `<output_dir>/volume.vtkhdf` contains
  `velocity` and `mpi_rank` fields.

### Phase 6 — SAFS fault static fields
- **T-41**: Run with `--paraview --paraview-fault hdf5 --tfinal 0.5s` on a
  spatially heterogeneous friction TOML.  Assert `<output_dir>/fault.vtkhdf`
  contains the 7 new static fields AND that the per-DOF values match
  `dof_data` (read back via the existing HDF5 mirror harness).
- **T-42**: TPV104 byte-exact regression (`make test-tpv104-checkpoint`) —
  the new `SetFaultParamsSpatial` method is additive, so TPV104's existing
  `SetFaultParamsBP5(a, Dc, x2, x3)` output must remain byte-equal.

### Regression gates that MUST stay green throughout
- `make test-tpv104-checkpoint` (164 / 164)
- `make test-tpv205-byteexact`
- `make test-bp5-smoke` (≥ 4 ranks)
- `seas_test_spatial_friction_config` (currently 28 / 28; will be 33 / 33
  after Phase 2)
- `seas_test_spatial_friction_resolver` (currently 110 / 110)
- All existing `tests/unit/test_paraview_*` + `test_fault_surface_*`
  + `test_vtkhdf_zfp` tests.

## Risk Assessment

### High-risk areas

1. **`SetFaultParamsSpatial` is the only `io/paraview_output.hpp` change in
   this plan.**  Mis-design here could ripple into TPV104/205 if the new
   method shares state with `SetFaultParamsBP5`.  Mitigation: implement the
   new method to write into a SEPARATE member container
   (`fault_params_spatial_`, NOT `fault_params_bp5_`), so the existing
   call sites see byte-equivalent output.  Verified by T-42.

2. **`Q` ordering assumption in Phase 4/5 (`SXX..VZ`)**: a future
   `dynamic/wave_state.hpp` reordering would silently corrupt the volume
   stress + velocity fields.  Mitigation: add a `static_assert` in the
   driver immediately above the memcpy block:
   ```cpp
   static_assert(VX == 6 && VY == 7 && VZ == 8,
                 "Phase 5 velocity memcpy assumes contiguous VX/VY/VZ; "
                 "wave_state.hpp ordering changed.");
   static_assert(SXX == 0 && SYY == 1 && SZZ == 2 && SXY == 3 && SYZ == 4 && SXZ == 5,
                 "Phase 4 sigma memcpy assumes the wave_state.hpp ordering.");
   ```

3. **Default flip from `hdf5` to `off` for `paraview_volume` and
   `paraview_bulk`** is user-visible: any existing SAFS TOML that depended
   on the old default would silently lose volume output.  Mitigation:
   - Search `safs/` and `jobs/safs/` for `.toml` files that omit
     `paraview_volume`/`paraview_bulk`; before the flip, update each one
     to spell out the previous default `"hdf5"` explicitly.
   - Document the flip prominently in `CLAUDE.md` "ParaView output mode"
     section and in the schema doc's footnote.

4. **Regime-adaptive schedule interaction with `paraview_fault_dt`**:
   When `paraview_coseismic_dt > 0` AND `paraview_fault_dt > 0`, the
   `AdaptiveSchedule::Interval()` returns the regime cadence and the
   `fixed_dt` member is ignored.  Document this priority in the schema
   doc; assert in T-39 that the regime cadence wins.

### Medium-risk areas

5. **`tfinal` units**: `cfg.time.tfinal` is already parsed in seconds by
   `toml_time_seconds`.  `SetTotalRunTime(cfg.time.tfinal)` consumes
   seconds.  Match.  No risk.

6. **`SetFaultParamsSpatial` empty-array edge case**: on a rank with
   `num_fault_total == 0`, the 7 vectors are empty; the new method must
   handle a zero-size input.  Covered by T-41 with a 1-rank setup that
   has all fault DOFs on one rank.

### Out of scope (deliberate)

- Restart of the SECONDARY `pv_bulk_out` collection's schedule state.
  TPV104 has the same V2 limitation (`tpv104_driver.cpp:2403-2412`); a V3
  checkpoint schema with two collections' state is future work.
- Per-fault-DOF "label" (which spatial rule applied to this DOF).  Useful
  for debugging spatial rule precedence but not in the user's current
  request.

## Files to Update — Summary Table

| File                                                          | Phase    | Change type                                            |
|---------------------------------------------------------------|----------|--------------------------------------------------------|
| `drivers/spatial_dyn_driver.cpp`                              | 1,3,4,5,6 | Additive — CLI flags, guards, banner, two collections, volume GFs, SAFS params |
| `spatial/code/spatial_friction.hpp`                           | 2         | Additive — 9 new `OutputSpec` fields, 2 default flips   |
| `spatial/code/spatial_friction.cpp`                           | 2         | Additive — TOML parse + validate the 9 new fields       |
| `safs/project_7.0_alternative/document/spatial_friction_config_schema.md` | 2 | Additive — document 9 new keys + the flipped defaults |
| `safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md` | (informational only) | Add a note pointing at this plan; do NOT mutate Phase 4 step list in the parent plan |
| `io/paraview_output.hpp`                                      | 6         | Additive — `SetFaultParamsSpatial` method only          |
| `tests/unit/test_spatial_paraview_cli.cpp`                    | 1         | NEW                                                    |
| `tests/unit/test_spatial_friction_config.cpp`                 | 2         | EXTEND — add T-33..T-37                                 |
| `tests/integration/test_spatial_paraview_default_off.cpp`     | 3,4,5,6   | NEW                                                    |
| `jobs/safs/safs_dyn_paraview_zfp_dev.sbatch`                  | Phase 6 follow-up | NEW (mirror of `tpv104_phase6_paraview_zfp_dev_2hr.sbatch`) |
