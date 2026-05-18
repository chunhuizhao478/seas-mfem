# ParaView I/O — Codebase Walkthrough

**Date:** 2026-05-11 (originally) — last expanded 2026-05-17
**Scope:** ParaView visualization output for the seas miniapp (BP5 + TPV102/104/205 drivers), including the Phase 6 VTKHDF+ZFP compaction work, the 2026-05-12 split-bulk-solutions refactor, the 2026-05-16 schedule-cap rework / zero-fault-rank crash fix, and the 2026-05-17 V2 restart support. Covers the wrapper `seas::ParaViewOutput<MeshType>`, its three underlying collections, the MFEM-level VTKHDF substrate, the ZFP compression layer, the new `petsc_ts_checkpoint.hpp` / `tpv104_checkpoint.hpp` machinery, and the driver-side glue.
**Audience:** Computational scientist who has never read this code before. The math is elsewhere — this document is exclusively about I/O.

---

## READ THIS FIRST: Opening Phase 6 Output on a macOS ParaView

The `.vtkhdf` files this codebase produces are ZFP-compressed. ParaView 5.11+ knows how to *render* them but does NOT bundle the H5Z-ZFP plugin needed to *decompress* them. Out of the box, opening one of our files in ParaView fails with:

```
HDF5: object 'Points' doesn't exist
HDF5: object 'Values' doesn't exist
Cannot read the Points_0_ array from file
ERROR: vtkHDFReader: Cannot read the Points array
```

That's not a file corruption — it's the HDF5 read path silently failing because the ZFP filter isn't registered. The plugin (`libh5zzfp.dylib`, filter id 32013) must be **discoverable at runtime via `HDF5_PLUGIN_PATH`** AND its dependencies must resolve without conflicting with ParaView's bundled HDF5.

### macOS-specific gotcha: `DYLD_LIBRARY_PATH` causes a segfault

The plugin from the `mfem-dev` conda env has `@rpath/libzfp.1.dylib` and `@rpath/libhdf5.310.dylib` as dependencies. Setting `DYLD_LIBRARY_PATH=<conda-env>/lib` to satisfy those drags conda's `libhdf5.310.dylib` into the ParaView process *alongside* ParaView's own bundled `libhdf5.310.dylib` — two copies of HDF5 in one process -> instant segfault on launch.

`DYLD_FALLBACK_LIBRARY_PATH` works only sometimes (macOS SIP / amfi may strip it). The reliable solution is to `install_name_tool` the plugin so its dependencies point at **absolute paths** baked in, then re-sign with an ad-hoc signature so Gatekeeper accepts it.

### One-time setup (run this once, then never again)

```bash
# 1. Build a private, self-contained copy of the plugin with absolute
#    dependency paths.
mkdir -p ~/.h5z-zfp-deps/plugin
cp /Users/chunhuizhao/miniforge/envs/mfem-dev/plugin/libh5zzfp.dylib \
   ~/.h5z-zfp-deps/plugin/
install_name_tool \
  -change @rpath/libzfp.1.dylib \
          /Users/chunhuizhao/miniforge/envs/mfem-dev/lib/libzfp.1.dylib \
  -change @rpath/libhdf5.310.dylib \
          /Applications/ParaView-6.0.0-RC3.app/Contents/Libraries/libhdf5.310.dylib \
  ~/.h5z-zfp-deps/plugin/libh5zzfp.dylib

# 2. Re-sign the modified dylib — install_name_tool invalidates the
#    code signature and macOS kills any process trying to load it.
codesign --force --sign - ~/.h5z-zfp-deps/plugin/libh5zzfp.dylib

# 3. Make HDF5_PLUGIN_PATH and a `paraview` alias permanent in zsh.
echo 'export HDF5_PLUGIN_PATH=$HOME/.h5z-zfp-deps/plugin' >> ~/.zshrc
echo 'alias paraview="/Applications/ParaView-6.0.0-RC3.app/Contents/MacOS/paraview"' >> ~/.zshrc
source ~/.zshrc
```

After that, opening any of our `.vtkhdf` files is a one-liner:

```bash
paraview fault_surface.vtkhdf &
```

**Note:** `~/.zshrc` covers only the terminal. To also open `.vtkhdf` by double-clicking in Finder (or via Spotlight), install a LaunchAgent that sets `HDF5_PLUGIN_PATH` system-wide:

```bash
mkdir -p ~/Library/LaunchAgents
cat > ~/Library/LaunchAgents/com.user.hdf5pluginpath.plist <<'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN"
   "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key><string>com.user.hdf5pluginpath</string>
  <key>ProgramArguments</key>
  <array>
    <string>launchctl</string><string>setenv</string>
    <string>HDF5_PLUGIN_PATH</string>
    <string>/Users/chunhuizhao/.h5z-zfp-deps/plugin</string>
  </array>
  <key>RunAtLoad</key><true/>
</dict>
</plist>
EOF
launchctl load ~/Library/LaunchAgents/com.user.hdf5pluginpath.plist
```

### Why this song and dance is necessary (3-step root cause)

1. **ZFP compression is non-optional.** Without it, BP5 production output is hundreds of GB; the whole point of Phase 6 was to bring it under a GB.
2. **HDF5 filters are dlopen'd at decompression time** from `HDF5_PLUGIN_PATH`. The plugin is not linked into ParaView.
3. **macOS Gatekeeper + the two-HDF5-in-one-process trap** make the obvious environment-variable solutions fail (segfault) or get blocked (kill). The `install_name_tool` + `codesign --sign -` recipe is the only reliable workaround for an unmodified ParaView .app.

If you ever reinstall ParaView or update conda's `libhdf5`, the absolute paths inside `~/.h5z-zfp-deps/plugin/libh5zzfp.dylib` will become stale — re-run the one-time setup. The plugin's `otool -L` output should show no `@rpath/...` entries when correctly set up.

---

## Overview

A single SEAS run on Frontera produces three ParaView artefacts (renamed
per `document/io_dev/PLAN_split_bulk_solutions_2026-05-12.md` — older
runs and out-of-tree sbatch may still show the legacy names in
parentheses):

| File | Contents | Purpose |
|---|---|---|
| `<output>/ParaView/kinematics.vtkhdf` (was `volume.vtkhdf`) | **Kinematics** — velocity (TPV*) or displacement (BP5) + mpi_rank | Animation in ParaView 5.11+ |
| `<output>/ParaView_bulk/stress.vtkhdf` (was `wave_bulk.vtkhdf`) | **Stress** — full symmetric tensor sigma_xx, _yy, _zz, _xy, _xz, _yz | Wavefield diagnostics; TPV* only (BP5 has no stress file — quasi-dynamic) |
| `<output>/fault.vtkhdf` (was `fault_surface.vtkhdf`) | **Fault interface** — slip, slip-rate, traction, state on the 2D fault | The science output |

Each collection has its own writer (back-end), cadence, compression tolerance, and ParaView-side time series. The wrapper `seas::ParaViewOutput<MeshType>` owns the kinematics + fault collections; the stress collection is a second instance of `seas::ParaViewOutput` created next to it in the TPV drivers.

The historical context (per `document/io_dev/PLAN_paraview_compaction_2026-04-28.md`): a 250-yr BP5 production run on Frontera at np=800 used to emit **3.9 million** per-rank-per-cycle ASCII VTU files into one Lustre directory, breaking `du`, `ls`, `tar`, and `rsync`. The compaction work replaced that with single-file VTKHDF (HDF5 + chunked + ZFP-compressed), reducing both file count (millions -> 3) and bytes (~25 GB -> ~200–400 MB) per run.

---

## How the split-bulk-solutions refactor works, and how to adjust outputs

Before the 2026-05-12 refactor, the TPV drivers wrote two volume files (`volume.vtkhdf` + `wave_bulk.vtkhdf`) that BOTH carried `velocity` + `mpi_rank` — a byte-for-byte duplication whenever the cadences matched (the default). The primary file ALSO carried 12 fault projections that duplicated data already in `fault_surface.vtkhdf`. Three problems:

1. ~50 % of `wave_bulk.vtkhdf` was a duplicate of `volume.vtkhdf`.
2. The 12 fault projections in `volume.vtkhdf` duplicated `fault_surface.vtkhdf`.
3. The names ("volume", "wave_bulk", "fault_surface") were physically vague.

The refactor (`document/io_dev/PLAN_split_bulk_solutions_2026-05-12.md`) does three things:

| Aspect | Old | New |
|---|---|---|
| Names | `volume.vtkhdf` / `wave_bulk.vtkhdf` / `fault_surface.vtkhdf` | `kinematics.vtkhdf` / `stress.vtkhdf` / `fault.vtkhdf` |
| Stress fields | 3 components (sigma_yy, sigma_xy, sigma_xz) | **6** — full symmetric tensor (sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz) |
| Velocity / mpi_rank | In BOTH volume + bulk | In `kinematics.vtkhdf` only |
| 12 L2-p0 fault projections (slip_dip, slip_strike, slip_rate_*, traction_*, state_variable, normal_stress, param_a, param_Dc, fault_x2, fault_x3) | Always in `volume.vtkhdf` | Default OFF; opt back in via `SetRegisterFaultProjectionsInVolumePV(true)` (see "Layer 1" below) |

### Three layers a user can adjust

Output content + cadence + compression are controlled at three layers, from least to most invasive:

#### Layer 1 — driver code, before the time loop (compile-time)

For driver authors. Each TPV driver constructs two `seas::ParaViewOutput<MeshT>` instances:

```cpp
// drivers/tpv102_driver.cpp ~ line 1904
pv_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
   output_dir + "/ParaView", pmesh, order,
   /*collection_name=*/"kinematics", volume_mode);
pv_out->RegisterDomainField("velocity",   pv_vel_gf.get());
pv_out->RegisterDomainField("mpi_rank",   pv_rank_gf.get());
pv_out->InitFaultOutputBP5(...);              // allocates the 12 GFs but
                                              // does NOT register them
                                              // by default (R-001 fix)

// drivers/tpv102_driver.cpp ~ line 2120 (only when --paraview-bulk-dt > 0)
pv_bulk_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
   output_dir + "/ParaView_bulk", pmesh, order,
   /*collection_name=*/"stress", volume_mode);
pv_bulk_out->RegisterDomainField("sigma_xx", pv_bulk_sxx_gf.get());
pv_bulk_out->RegisterDomainField("sigma_yy", pv_bulk_syy_gf.get());
// ... 4 more sigma components
```

To add a field, add a GridFunction and call `RegisterDomainField("name", gf)` on the appropriate collection before any `Save()`. To remove one, comment out its `RegisterDomainField` call. To put the 12 fault projections back into `kinematics.vtkhdf`, call `pv_out->SetRegisterFaultProjectionsInVolumePV(true)` **BEFORE** `InitFaultOutputBP5(...)`. The default is OFF.

The unit test `seas_test_kinematics_field_set` (`tests/unit/test_kinematics_field_set.cpp`) guards both scenarios — flipping the default to true makes the test fail 13/31.

#### Layer 2 — sbatch CLI flags (per-run)

For job submitters. Per-run cadence and ZFP tolerance go on the `ibrun` line in the sbatch:

```bash
ibrun ./seas_tpv102_driver \
    ...
    --paraview                          # enable PV output
    --paraview-dt 0.5                   # kinematics + fault cadence (s)
    --paraview-bulk-dt 0.5              # stress cadence (s; 0 = no stress file)
    --paraview-volume-zfp-tol 1e-3      # ZFP abs-error on kinematics.vtkhdf
    --paraview-bulk-zfp-tol   1e-3      # ZFP abs-error on stress.vtkhdf
    --paraview-fault-zfp-tol  1e-12     # ZFP abs-error on fault.vtkhdf
    --paraview-max-snapshots  200       # soft cap on total fault writes
    --no-volume-pv                      # BP5 production: suppress kinematics
```

The volume flag prefix is still `--paraview-volume-*` and `--paraview-bulk-*` after the refactor — the driver-side rename to `--paraview-kinematics-*` / `--paraview-stress-*` is queued as a follow-up. Same flag semantics; only the collection name on disk changed.

To **drop the stress file entirely**, omit `--paraview-bulk-dt` (the secondary collection isn't constructed). To **drop the kinematics file**, pass `--no-volume-pv`. To skip `fault.vtkhdf`, well… don't pass `--paraview` (everything turns off).

Independent ZFP tolerances are the architectural reason both volume files exist:

| Field | Dynamic range | Recommended tol | Why |
|---|---|---|---|
| velocity (kinematics) | ~1e-6..10 m/s, mostly O(1) | `1e-3` | ±1 mm/s invisible at coseismic peak |
| stress components | ~1e6..1e8 Pa | `1e-3` (Pa) | ~10 decades below peak; safe |
| slip rate (fault) | 1e-9..1e0 m/s during nucleation->event | `1e-12` | 3 decades below smallest interesting value |

Picking tolerances bigger than these will quantize interesting features. Picking smaller wastes storage with no science benefit.

#### Layer 3 — TOML config (future, deferred per R-313)

For `seas_driver` users, once PV wiring lands there. The planned schema:

```toml
[paraview]
enabled = true

[paraview.kinematics]
fields    = ["velocity", "mpi_rank"]
zfp_tol   = 1e-3
dt        = 0.5
backend   = "hdf5"

[paraview.stress]
fields  = ["sigma_xx", "sigma_yy", "sigma_zz",
           "sigma_xy", "sigma_xz", "sigma_yz"]
zfp_tol = 1e-3
dt      = 0.5

[paraview.fault]
zfp_tol = 1e-12
backend = "hdf5"
```

Currently `seas_driver.cpp` doesn't construct a `seas::ParaViewOutput` (R-313 deferred deviation; see `CLAUDE.md` § "Phase 4 deferred deviation"). The TOML schema lands here as a future-use spec; the parsing wiring follows when BP5 production runs through `seas_driver` actually need a `kinematics.vtkhdf`.

### Quick recipes

- **"I want stress data but no velocity"** -> add `--paraview-bulk-dt 0.5` to the sbatch, ALSO pass `--no-volume-pv` to suppress `kinematics.vtkhdf`. `stress.vtkhdf` will be the only volume file.
- **"I want lossless ZFP"** -> pass `--paraview-{volume,bulk,fault}-deflate-level 6` (omit the `-zfp-tol` flags). Output is ~2-3× larger.
- **"I want the fault projections in the kinematics file for visualization"** -> in your driver code, call `pv_out->SetRegisterFaultProjectionsInVolumePV(true)` before `InitFaultOutputBP5(...)`. No sbatch change needed. Note this re-introduces ~12× duplicate data with `fault.vtkhdf`.
- **"I want a different sigma component subset"** -> patch the driver: comment out the unwanted `pv_bulk_out->RegisterDomainField("sigma_XX", ...)` calls in `tpv*_driver.cpp` (around line 2141–2148 for tpv102). Then the corresponding memcpy in `paraview_write` is dead code but harmless.

---

## Architecture

### Three back ends, three switches

```
Build flag            Default fault back end       Default volume back end
-----------------     --------------------------    -------------------------
MFEM_USE_HDF5=YES     FaultOutputMode::Hdf5         VolumeOutputMode::Hdf5
                       -> mfem::ParaViewHDFData       -> mfem::ParaViewHDFData
                         Collection                    Collection
                       -> fault_surface.vtkhdf        -> volume.vtkhdf
MFEM_USE_HDF5=NO      FaultOutputMode::Vtu          VolumeOutputMode::Vtu
                       -> single binary VTU/cycle     -> mfem::ParaViewData
                         (gather-to-rank-0)            Collection (per-rank VTU)
```

CLI overrides on every driver (`tpv102`, `tpv104`, `tpv205`):
- `--paraview-fault-vtu` / `--paraview-fault-hdf5` — pick fault back end
- `--paraview-volume-vtu` / `--paraview-volume-hdf5` — pick volume back end
- `--paraview-fault-legacy-ascii` — restore the pre-Phase-1 per-rank ASCII writer (debugging only)
- `--no-volume-pv` — suppress the primary volume file entirely (BP5 production default)

### Key files

```
miniapps/seas/io/
├── paraview_output.hpp       ~2280 lines — main wrapper template
│                              `seas::ParaViewOutput<MeshType>`.
│                              2026-05-12: SetRegisterFaultProjectionsInVolumePV
│                              + GetVolumeDataCollection.
│                              2026-05-17: Set/Get accessors for V2
│                              restart (TotalSnapshotsWritten,
│                              RestoreScheduleState, LastCommittedCycle,
│                              LastVolumeWriteTime, ...).
├── fault_vtu_binary.hpp        735 lines — Phase 1 gather-to-rank-0
│                              binary VTU + LocalFaultPack /
│                              GatheredFaultPack data structures.
├── fault_vtkhdf_writer.hpp     304 lines — Phase 2b rank-0 VTKHDF
│                              writer for the fault surface.
├── petsc_ts_checkpoint.hpp     ~220 lines (new 2026-05-17) — V2
│                              `PETSC_TS_V2` trailing block reader /
│                              writer.  Appends to the per-rank V1
│                              checkpoint produced by checkpoint.hpp's
│                              WriteCheckpoint.  Used by
│                              bp5_verification_full.cpp on the
│                              `--restart` + `--petsc-ts` path.
└── tpv104_checkpoint.hpp       ~370 lines (new 2026-05-17) — TPV104
                                V1 checkpoint format (distinct magic
                                tag TPV104_CHECKPOINT_V1).  Per-rank
                                file, round-trips wave-field Q + 9
                                dynamic DOFData fields.  No ParaView
                                schedule state yet (deferred V3).

mesh/vtkhdf.{hpp,cpp}           Low-level VTKHDF (HDF5) writer.
                                Used by ParaViewHDFDataCollection.
                                MFEM-level; not seas-specific.

fem/datacollection.hpp/cpp      mfem::ParaViewDataCollection (legacy
                                per-rank VTU) and
                                mfem::ParaViewHDFDataCollection
                                (single-file VTKHDF).

miniapps/seas/drivers/
├── tpv102_driver.cpp           Constructs `pv_out` and (conditionally)
├── tpv104_driver.cpp           `pv_bulk_out`; calls into the wrapper
├── tpv205_driver.cpp           inside the time loop.
└── seas_bp5_full doesn't construct ParaViewOutput directly; the
    BP5 verification driver (tests/verification/bp5_verification_full.cpp)
    constructs it and calls Save/WriteFaultSurfaceVTU itself.
```

### Data flow per write cycle (parallel build, HDF5 mode)

```
                                 ┌─────────────────────────────────┐
                                 │  PRIMARY (volume.vtkhdf)        │
                                 │  ─────────────────────────────  │
displacement /  ──────────────->  │  pv_dc_->Save()                 │
velocity GF                      │  -> ParaViewHDFDataCollection    │
                                 │    -> VTKHDF::SaveMesh + fields  │
                                 │    -> parallel collective HDF5   │
                                 │      writes into one .vtkhdf    │
                                 │      (HDF5 dataspace chunked,   │
                                 │       ZFP filter applied to fp) │
                                 └─────────────────────────────────┘

per-DOFData fault    ┌─────────────────────────────────┐
fields (slip, V,     │  FAULT (fault_surface.vtkhdf)   │
traction, state)     │  ─────────────────────────────  │
                     │  GatherFaultPackToRoot(comm)    │
all ranks -> rank 0   │  ─────────────────────────────  │
                     │  rank 0 only:                   │
                     │  state.dc->Save()               │
                     │  -> ParaViewHDFDataCollection    │
                     │    (built on MPI_COMM_SELF)     │
                     │    -> serial HDF5 write          │
                     │  non-root ranks return after    │
                     │  the gather.                    │
                     └─────────────────────────────────┘

σyy, σxy, σxz GFs    ┌─────────────────────────────────┐
(TPV* only)          │  SECONDARY (wave_bulk.vtkhdf)   │
                     │  ─────────────────────────────  │
                     │  pv_bulk_out->ForceSave()       │
                     │  -> independent                  │
                     │  ParaViewHDFDataCollection      │
                     │  built on the SAME ParMesh as   │
                     │  pv_out but a different file    │
                     │  + different ZFP tolerance      │
                     └─────────────────────────────────┘
```

### External dependencies

| Library | Used for |
|---|---|
| HDF5 (parallel) | VTKHDF file substrate. Required: `H5_HAVE_PARALLEL` => `MFEM_PARALLEL_HDF5`. ParMesh runs `static_assert` if missing (`paraview_output.hpp:126`). |
| ZFP | Lossy floating-point compression. `extern/zfp/install/lib64/libzfp.so`. |
| H5Z-ZFP | HDF5 filter plugin id 32013 that wraps ZFP. Loaded at runtime via `HDF5_PLUGIN_PATH`; **not** linked into MFEM. |
| zlib | Lossless fallback for integer datasets (connectivity, offsets, cell types). Available unconditionally with HDF5. |

The `extern/` directory is populated by `build_frontera.sh`; the runtime plugin path is set in each Phase 6 sbatch under `miniapps/seas/jobs/{bp5,tpv*}/`.

---

## Call Graph (critical path, TPV* driver)

```
tpv*_driver.cpp::main()
  └-> argv parsing (--paraview-* family)
  └-> ParaViewOutput<ParMesh> pv_out(prefix, pmesh, order,
                                     "volume", volume_mode)
  └-> pv_out->RegisterDomainField("velocity", ...)
  └-> pv_out->InitFaultOutputBP5(interior_faces, shared_faces, nbf)
  └-> pv_out->SetFaultHDFCompression(ZfpAccuracy, fault_zfp_tol)
  └-> pv_out->SetVolumeHDFCompression(ZfpAccuracy, volume_zfp_tol)
  └-> pv_out->SetVolumeSaveEnabled(use_volume_pv)
  └-> pv_out->SetVolumePVDt(volume_pv_dt)
  └-> pv_out->GetSchedule().{v_coseismic, dt_coseismic, ...} = CLI overrides
  └-> pv_out->GetSchedule().Validate()
  └-> pv_out->SetTotalRunTime(tfinal)         // required for max_total_snapshots
  └-> [optional] pv_bulk_out = ParaViewOutput<ParMesh>(prefix+"/ParaView_bulk",
                                                       pmesh, order, "wave_bulk")
  └-> pv_bulk_out->fixed_dt = bulk_dt
  └-> pv_bulk_out->RegisterDomainField("velocity"/"sigma_yy"/...)
  └-> pv_bulk_out->SetVolumeHDFCompression(ZfpAccuracy, bulk_zfp_tol)
  │
  └-> time loop:
       └-> paraview_write(step_num, time, V_max):
            ├-> pv_out->PeekShouldWrite(...)         // schedule probe, NO mutation
            ├-> pv_bulk_out && pv_bulk_out->PeekShouldWrite(...)
            ├-> pack pv_vel_gf, pv_bulk_*_gf from solver state
            ├-> pv_bulk_out->ForceSave(step_num, time)
            │   └-> pv_dc_->Save()  (ParaViewHDFDataCollection)
            │        └-> VTKHDF::SaveMesh + AppendField
            │             └-> collective HDF5 writes (chunked + filtered)
            ├-> pv_out->UpdateFaultFieldsBP5(...) -> fills L2-p0 fault GFs
            ├-> pv_out->ForceSave(step_num, time)
            │   └-> ForceSaveImpl: if (emit_volume_save_) pv_dc_->Save()
            ├-> pv_out->CommitSchedule(time, V_max)  // advance regime FSM
            └-> pv_out->WriteFaultSurfaceVTU(...):
                  ├-> build LocalFaultPack from per-DOFData
                  ├-> vtu::GatherFaultPackToRoot(local, rank, nranks, comm)
                  ├-> if FaultOutputMode::Hdf5:
                  │     vtkhdf::WriteFaultPackHdf(state, prefix,
                  │                                cycle, time, ...)
                  │       └-> rank 0: state.dc->Save() (serial HDF5)
                  │       └-> non-root: early return
                  └-> if FaultOutputMode::Vtu:
                        vtu::WriteFaultPackVTU(path, gathered, BINARY)
                          └-> rank 0 writes one .vtu, appends PVD entry
```

### Call graph (BP5)

`seas_bp5_full` (= `tests/verification/bp5_verification_full.cpp`) is similar but:

- No secondary `pv_bulk_out`.
- Adaptive-schedule (V_max-driven) is the default cadence; `--paraview-dt`/`--paraview-every` override it. Per-regime overrides use **BP5-specific** flag names `--paraview-dt-co / -nu / -inter-yr` (TPV* uses `--paraview-{co,nucleation,inter}seismic-dt`).
- Production runs default to `--no-volume-pv` (or `--paraview-fault-only`) — only `fault.vtkhdf` and station-CSV are emitted.
- **V2 restart hook** (new 2026-05-17): when `--restart PREFIX` AND `--petsc-ts` AND the V2 trailing block is present, the driver's restart block at `bp5_verification_full.cpp:2583–2691` calls `ReadPetscTSCheckpoint` then (inside `if (pv_out)`) the three setters `pv_out->SetTotalSnapshotsWritten(...)`, `pv_out->RestoreScheduleState(...)`, `pv_out->SetLastCommittedCycle(...)` BEFORE the time loop starts. The monitor callback at lines 643–679 and the end-of-run write at lines 3020–3058 call `WritePetscTSCheckpoint` with the five `ParaViewOutput` getters (`GetTotalSnapshotsWritten / GetLastWriteTime / GetLastVMax / GetCurrentRegime / GetLastCommittedCycle / GetLastVolumeWriteTime`).

---

## Detailed Walkthrough

### `seas::ParaViewOutput<MeshType>` (paraview_output.hpp:69 — 2093 lines)

Template parameter `MeshType ∈ {Mesh, ParMesh}` selects the corresponding `GridFunction` / `FiniteElementSpace` types via the `GFType` trait at the top of the file.

#### Construction

```cpp
ParaViewOutput(prefix, mesh, order,
               collection_name = "volume",
               mode = DefaultVolumeOutputMode())
```

- `prefix` — output directory.
- `mesh` — domain mesh (held by reference; not owned).
- `order` — polynomial order, passed to `SetLevelsOfDetail`.
- `collection_name` — base filename for HDF5 mode (`<prefix>/<name>.vtkhdf`) or subdir for VTU mode.
- `mode` — locks the volume back end for the lifetime of the wrapper.

`DefaultVolumeOutputMode()` is a `static constexpr` function that returns `Hdf5` whenever the build supports it for the current `MeshType`, otherwise `Vtu`. For `ParMesh` on a build without `MFEM_PARALLEL_HDF5`, it triggers a `static_assert` rather than silently falling back to VTU (R-401, paraview_output.hpp:126–131). This is intentional: the silent fallback used to mask every Frontera ParMesh run dropping into the legacy per-rank VTU path.

#### Three writers, three responsibilities

| Member | Type | When written | Lifecycle |
|---|---|---|---|
| `pv_dc_` | `unique_ptr<ParaViewDataCollectionBase>` (pointing at either `ParaViewDataCollection` or `ParaViewHDFDataCollection`) | every `ForceSaveImpl` when `emit_volume_save_` AND the user dt elapsed | constructed in the ctor; survives the run |
| `fault_hdf_state_` | `vtkhdf::FaultHDFState` (rank-0 only) | every `WriteFaultSurfaceVTU` when `FaultOutputMode::Hdf5` | constructed lazily on first write |
| (binary VTU path) | none — stateless emit via `vtu::WriteFaultPackVTU` | every `WriteFaultSurfaceVTU` when `FaultOutputMode::Vtu` | n/a |

The fault-surface and primary-volume writers DO NOT share files. They share the schedule (`adaptive_`, `last_write_time_`, `current_regime_`) but their on-disk artefacts are independent.

#### Schedule state machine — `AdaptiveSchedule`

A 3-regime finite state machine keyed on the global max slip rate `V_max`:

```
regime 0 (interseismic)   V ≤ v_nucleation       interval = dt_interseismic (default 1 yr)
regime 1 (nucleation)     v_nucleation < V       interval = dt_nucleation   (default 1 s)
                          ≤ v_coseismic
regime 2 (coseismic)      V > v_coseismic        interval = dt_coseismic    (default 0.01 s)
```

State transitions in `NextRegime(V, prev_regime)`:
- Entries use strict `>` comparisons.
- Exits use `V / hysteresis_factor` (≥ 1.0). `hysteresis_factor = 1.0` (default) disables hysteresis.
- Non-finite `V` -> regime 2 (so we keep dense output during a blowup).
- "Drop-past-threshold": a coseismic->V drop that lands in the nucleation band skips regime 1 and goes directly to interseismic. Rationale at `paraview_output.hpp:232–243`: avoids over-sampling post-rupture ring-down.

#### Three save methods, three semantics

| Method | Mutates schedule | Calls `pv_dc_->Save()` | When to use |
|---|---|---|---|
| `Save(cycle, time, V_max)` | yes (regime + last_write_time + last_v_max) | yes (subject to `emit_volume_save_`) | Standard one-shot "schedule decides + write everything" |
| `ShouldWrite(cycle, time, V_max)` | yes | no | "Schedule decides but skip the volume" — fault-only with the snapshot counter still bumping |
| `PeekShouldWrite(cycle, time, V_max) const` | no | no | Schedule probe before doing expensive packing work. Must be followed by `CommitSchedule(time, V_max)` to advance state |
| `ForceSave(cycle, time)` | bumps `last_write_time_` + snapshot counter | yes (subject to `emit_volume_save_`) | Caller decides ("write now") — bypass schedule |
| `CommitSchedule(time, V_max)` | yes | no | Advance schedule state without writing volume (used after `PeekShouldWrite` + `WriteFaultSurfaceVTU` only) |

The TPV* drivers use `PeekShouldWrite -> ForceSave/CommitSchedule -> WriteFaultSurfaceVTU` so they can decide whether to copy fault data into the registered GFs (an O(num_dofs) memcpy) on each step.

#### Snapshot cap — `max_total_snapshots`

If set (default 0 = uncapped), `SnapshotCapAwareInterval` inflates `dt_interseismic` to keep the projected total snapshot count under the cap:

```
remaining_budget = K - snapshots_so_far
new_dt_inter     = max(dt_inter, (T - t) / remaining_budget)
```

Coseismic and nucleation cadences are never inflated; they're physically meaningful. The cap is "soft" — a long-running event may push the final count slightly over K (≤ 1 extra per event, per design).

The cap requires `SetTotalRunTime(tfinal)` — without it the cap is a no-op and the constructor asserts loudly rather than silently degrading.

#### Fault-side state — `InitFaultOutputBP5` and `UpdateFaultFieldsBP5`

`InitFaultOutputBP5(fault_interior_faces, fault_shared_faces, nbf_per_face)` constructs the L2-p0 fault GridFunctions on the **domain** mesh:

- 8 scalar GFs: `fault_slip_{dip,strike}`, `fault_slip_rate_{dip,strike}`, `fault_trac_{dip,strike}`, `fault_state`, `fault_normal_stress`.
- One DOF per element (L2 p=0).
- Fault values are scattered from the per-face DOFData layout to the two volume elements adjacent to each fault face — only elements touching a fault face carry nonzero values.
- For shared faces, only the local-side `Elem1` carries the value (`Elem2` lives on a neighbour rank).

These GFs are registered with the primary `pv_dc_` so they appear in `volume.vtkhdf` alongside displacement/velocity. They are NOT the same data the fault-surface writer emits — see "Gotchas" below.

`UpdateFaultFieldsBP5` is called every output cycle to refresh the GF values from the latest DOFData snapshot.

#### Fault-surface writer — `WriteFaultSurfaceVTU`

Five branches:

1. **Legacy ASCII** (`legacy_ascii_ = true`): one ASCII VTU per rank per cycle + per-cycle PVTU + cumulative PVD. Pre-Phase-1 layout. Kept verbatim because four unit tests depend on the ASCII regex shape (`test_fault_surface_vtu_continuity` and friends).

2. **Phase 1 binary VTU** (`output_mode_ = Vtu`, default on `MFEM_USE_HDF5=NO`): every rank packs a `vtu::LocalFaultPack`, MPI-gathers to rank 0 via `GatherFaultPackToRoot`, rank 0 writes ONE `<prefix>/FaultSurface/fault_surface_c<cycle>.vtu` (binary, base64-inline), appends a PVD entry, writes PVD index. File count drops from `nranks × ncycles` to `1 × ncycles`.

3. **Phase 2b VTKHDF** (`output_mode_ = Hdf5`, default on HDF5 builds): same gather, but rank 0 calls `vtkhdf::WriteFaultPackHdf` which feeds the gathered pack into a `ParaViewHDFDataCollection` constructed against a serial in-memory `mfem::Mesh` (one triangle per fault face). One `.vtkhdf` file for the entire run.

4. **`has_fault_output_ == false`** (collective early return): all ranks return — the wrapper was constructed but `InitFaultOutputBP5` was never called.

5. **Compile-time guards**: on `MFEM_USE_HDF5=NO` builds, `output_mode_ == Hdf5` aborts at runtime; CLI arg parsing in the drivers rejects `--paraview-fault-hdf5` at parse time so the abort is unreachable from external callers.

#### Uniformity check — paraview_output.hpp:1076–1095

The first `WriteFaultSurfaceVTU` call on a parallel build broadcasts `(legacy_ascii_, output_mode_)` from rank 0 and asserts every rank matches. If `SetFaultOutputMode`/`SetLegacyAsciiVTU` are called from a single rank only, the gather collective would deadlock; this check converts the deadlock into a loud abort. The check is cached (`output_mode_uniformity_checked_`) so the broadcast doesn't fire 4M times on a 5000-cycle BP5 run.

---

### `vtkhdf::FaultHDFState` and `WriteFaultPackHdf` (fault_vtkhdf_writer.hpp)

The Phase 2b fault writer **deliberately deviates from the plan**. The plan specified a `ParSubMesh`-based per-rank fault submesh + collective parallel HDF5 write. The deviation note at fault_vtkhdf_writer.hpp:22–47 explains:

- `ParSubMesh::CreateFromBoundary` requires the submesh source to be a boundary attribute; fault faces are *interior*. Promoting them to boundary would invalidate the elasticity operator's attribute contract.
- The fault surface is at most a few thousand triangles even at BP5 scale; parallel-HDF5 throughput gain is negligible vs reusing the already-tested Phase 1 gather.
- The on-disk result is identical: one `.vtkhdf` file per run.

Implementation:

- All ranks call `GatherFaultPackToRoot(local, rank, nranks, comm)` collectively. Non-root ranks then return.
- Rank 0 constructs a serial 2D-in-3D `mfem::Mesh` of triangles on first call (`Init(prefix, g)`), with one triangle per fault face, L2-p0 finite element space, one GridFunction per emitted scalar field, and a `ParaViewHDFDataCollection` on the serial mesh.
- Subsequent calls reuse the mesh + GFs (geometry is static across cycles, asserted via `MFEM_VERIFY`).
- The chunk filter (Deflate level OR ZFP tolerance) is set on the data collection at Init time; cannot change mid-run.
- Each cycle: copy gathered field values into the L2-p0 GFs, `dc->SetCycle(cycle); dc->SetTime(time); dc->Save();`. Inside MFEM, `ParaViewHDFDataCollection::EnsureVTKHDF` dynamic_casts the mesh as serial, so HDF5 calls are non-collective (`fem/datacollection.cpp:1374–1402`).

### `vtu::LocalFaultPack` / `GatherFaultPackToRoot` (fault_vtu_binary.hpp)

Three data structures + three helpers:

| Type | Layout |
|---|---|
| `LocalFaultPack` | `vertices[3*ncells]`, `triangles[ncells]={3i,3i+1,3i+2}`, `field_arrays[k][ncells]`, `field_names[k]` |
| `GatheredFaultPack` | Concatenations of per-rank LocalFaultPacks with triangle indices remapped to point at the global vertex array |

Invariant: `triangles[i] == {3i, 3i+1, 3i+2}` — each triangle owns three unique vertices, no inter-cell sharing. This matches the legacy per-face emit and makes the gather's index remap a constant offset add.

The gather is collective on the supplied communicator (or a no-op on serial builds). The field-name list is OR-reduced across ranks so every rank packs the same fields — the `_k4` diagnostic fields are present iff at least one rank has them.

`WriteFaultPackVTU(path, g, format, compression_level)` writes the gathered pack as a binary-with-base64-inline VTU. The format deviation (inline-base64 vs raw-appended) is documented at fault_vtu_binary.hpp:26–41 — chose inline-base64 because (a) it reuses MFEM's already-tested writers and (b) `MFEM_USE_ZLIB=YES` recovers the size delta vs raw-appended.

---

### MFEM-level `ParaViewHDFDataCollection` (fem/datacollection.{hpp,cpp})

The MFEM-side counterpart of `ParaViewDataCollection`. The class:

```cpp
class ParaViewHDFDataCollection : public ParaViewDataCollectionBase {
   std::unique_ptr<VTKHDF> vtkhdf;
   void EnsureVTKHDF();
   template <typename FP_T> void TSave();
public:
   ParaViewHDFDataCollection(name, mesh = nullptr);
   void SetCompression(bool) override;
   enum class HDFCompression { Deflate, ZfpAccuracy };
   void SetHDFCompression(HDFCompression alg, double param);
   void Save() override;
};
```

- One `.vtkhdf` file per collection, lazily constructed on first `Save()` via `EnsureVTKHDF()`.
- Two compression modes:
  - `Deflate` (default at level 6): lossless zlib (preceded by HDF5 byte shuffle).
  - `ZfpAccuracy`: LLNL ZFP filter id 32013 at absolute-error tolerance `param`. Applied **only** to FP32/FP64 datasets; integer connectivity falls back to lossless deflate via `VTKHDF::EnsureDataset`.
- Filter selection is captured on the collection and applied to the underlying `VTKHDF` on the next `Save()`. It is **locked at first save** — calling `SetHDFCompression` after the first save updates the collection state but does not retroactively re-filter previously-written timesteps.

### MFEM-level `VTKHDF` (mesh/vtkhdf.{hpp,cpp})

Low-level VTKHDF writer. The seas miniapp doesn't call it directly; everything goes through `ParaViewHDFDataCollection`. Useful facts for debugging:

- One HDF5 file per `VTKHDF` instance. Opens with `H5F_ACC_TRUNC` by default (no restart yet; planned follow-up).
- Hierarchy: `/VTKHDF/{Points, Connectivity, Offsets, Types, NumberOfPoints, NumberOfCells, NumberOfConnectivityIds, PointData/<name>, CellData/<name>, Steps, Steps/PartOffsets, ...}`.
- Datasets grow on the first dimension via `H5S_UNLIMITED`. Each `Save()` extends and appends.
- Chunked storage with chunk size targeting ~0.5 MB. Filter chain attached at `EnsureDataset` time (lossless for integers, ZFP-or-deflate for floats).
- Phase 2d requirements pulled into MFEM proper (`H5S_BLOCK`, `unsigned long long` hsize_t) bumped the minimum HDF5 to **1.14.x**. On Frontera, this forced a from-source build at `extern/hdf5/install` (the intel-19-compatible module caps at 1.12.2).
- A type-id template `TypeID<T>` maps C++ types to HDF5 type ids. `hsize_t`'s underlying type changed from `uint64_t` (= `unsigned long` on Linux x86_64) to `unsigned long long` in 1.14, requiring exactly that specialization (mesh/vtkhdf.cpp:32–36).

---

### `build_frontera.sh` + Phase 6 sbatch glue

The build- and run-time substrate that makes the Phase 6 path work on Frontera. None of this is in the C++ source — but the C++ won't link or run without it.

**`build_frontera.sh`** (project root) builds three from-source extern installs:

1. **HDF5 1.14.6** -> `extern/hdf5/install` (parallel build against intel/19 + impi/19; Frontera's phdf5 module caps at 1.12.2 which is too old).
2. **ZFP 1.0.1** -> `extern/zfp/install` (LLNL's reference encoder/decoder; built with CMake, mpicc, intel/19).
3. **H5Z-ZFP v1.1.1** -> `extern/h5z-zfp/install/plugin/libh5zzfp.so` (HDF5 filter plugin id 32013 that wraps ZFP; built with `make HDF5_HOME=... ZFP_HOME=...`).
4. **gmsh 4.13.1** -> `extern/gmsh/bin/gmsh` (built from source with `gcc/g++`; needed only for `.geo -> .msh` regen at sbatch submit time, not for MFEM/seas itself).

MFEM is then configured with:

```
MFEM_USE_HDF5=YES
MFEM_USE_H5Z_ZFP=YES        # required so #define lands in _config.hpp
HDF5_OPT=-I${HDF5_INSTALL}/include
HDF5_LIB=-L${HDF5_INSTALL}/lib -lhdf5_hl -lhdf5 -lz
```

The H5Z-ZFP plugin is **not** linked into MFEM. HDF5 discovers it at runtime via `HDF5_PLUGIN_PATH`. Each Phase 6 sbatch sets:

```bash
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"
```

If the plugin isn't on `HDF5_PLUGIN_PATH` at run time, HDF5 fails the first ZFP write with "required filter is not registered". `seas::ParaViewOutput::ProbeH5ZZfpPluginOrAbort()` (paraview_output.hpp:503–546) probes the env-var before the first save and aborts with an actionable message instead.

---

## Phase chronology (from `git log`, `document/io_dev/PLAN_*.md`)

| Phase | Commit | What landed |
|---|---|---|
| Pre-Phase-1 | 297d249 | Initial per-rank ASCII VTU + cumulative PVD for BP5. |
| | 1cf3291–c409afb | Misc tightening (fault PVD time series, friction parameter fields, per-DOF mapping audit). |
| | 079e348 | `--paraview-dt X` (fixed time interval, overrides step interval). |
| | a78c623, b2b9531 | Adaptive V_max schedule with hysteresis (BP5 only). |
| Phase 1 | 4f856e7 | Gather-to-rank-0 binary appended VTU writer. One file per cycle instead of `nranks` files per cycle. |
| Phase 2b | (same WIP) | Single-file VTKHDF for the fault. Rank-0-only serial HDF5 inside a `ParaViewHDFDataCollection` on `MPI_COMM_SELF`. |
| Phase 2d | de73e7d (folded in) | H5Z-ZFP lossy compression on the fault path. Adds `--paraview-fault-zfp-tol`. |
| Phase 3 | (folded in) | `max_total_snapshots` cap + `--paraview-{co,nucleation,inter}seismic-dt` CLI overrides. |
| Phase 4 | (folded in) | `--no-volume-pv` / `--paraview-fault-only`, decoupling fault from volume PV. |
| Phase 6 | de73e7d | **Primary volume goes VTKHDF by default** on `MFEM_USE_HDF5=YES` builds. `--paraview-volume-hdf5/--paraview-volume-vtu` CLI selector. `SetVolumeHDFCompression(...)` actually applies. Secondary `pv_bulk_out` reroutes `--paraview-bulk-*` from `pv_out` to a new collection in `<prefix>/ParaView_bulk/wave_bulk.vtkhdf`. |
| Phase 6 runtime fixes | (multiple commits) | HDF5 1.14.6 from-source on Frontera; gmsh from-source; `module load python3`; SEAS_MFEM_ROOT auto-detect in sbatch; `MFEM_USE_H5Z_ZFP` placeholder in `config/config.mk.in`. |
| Split-bulk refactor | 7c16e62 / 64ea4cf / 78fe29d / b69f7ca | File renames `volume`→`kinematics`, `wave_bulk`→`stress`, `fault_surface`→`fault`. Stress collection now carries the FULL 6-component symmetric tensor. Velocity/mpi_rank deduplicated. 12 L2-p0 fault projections OFF-by-default on kinematics; opt back in via `SetRegisterFaultProjectionsInVolumePV(true)`. New `seas_test_kinematics_field_set` unit test (31/31). New `GetVolumeDataCollection()` read-only accessor. |
| Phase 6 ParaView hardening | 407f456 (2026-05-16) | **Hard cap reverted to SOFT (interseismic-only)** — coseismic/nucleation regimes always use natural cadence. R-001 rank-0 gate on cap-exhausted warning (was N-rank flood). R-005 honest warning text (geometric-halving cadence). ResetMemory crash on zero-fault-DOF ranks fixed via driver-side padded-shrink. Trace + fault_dof_coords debug outputs gated behind `SEAS_DEBUG_FACE_TRACE` / `SEAS_DEBUG_FAULT_DOF_COORDS` env vars (default OFF). Production sbatch passes `--paraview-dt-co/-nu/-inter-yr` explicitly. New tests: `seas_test_paraview_schedule_cap` (23/23), `seas_test_bp5_petsc_ts_zero_fault_rank` (8/8), `seas_test_paraview_rank0_warning_gate` (MPI-only). |
| BP5 V2 restart + PETSc TS | 19e128e + a921246 (2026-05-17) | `--restart` + `--petsc-ts` now works end-to-end. New `io/petsc_ts_checkpoint.hpp` adds a `PETSC_TS_V2` trailing block (10 fields: 4 PETSc TS + 5 ParaView schedule + cumulative rejections) appended to the existing per-rank V1 checkpoint. New `ParaViewOutput::{SetTotalSnapshotsWritten, RestoreScheduleState, SetLastCommittedCycle}` setters + 5 getters thread schedule state across the restart seam (R-304 / R-004 / R-006). New unit test `seas_test_bp5_petsc_ts_restart` (57/57) + cluster sbatch `bp5_restart_test_v2_dev_2hr.sbatch` (two-phase fresh→restart inside one reservation). |
| TPV104 V1 checkpoint | 63373f8 + 6cbbb9d (2026-05-17) | New `io/tpv104_checkpoint.hpp` (TPV104-distinct magic tag `TPV104_CHECKPOINT_V1`, per-rank file like BP5). Round-trips wave-field Q + 9 dynamic DOFData fields. No ParaView-state extension yet (Phase-4b follow-up). New `seas_test_tpv104_checkpoint`. Makefile auto-detects `-lstdc++fs` for the `<filesystem>` use in `--restart` / `--output-dir` safety check. |

The bulk-compression plan (`PLAN_bulk_compression_and_size_estimator_2026-05-09.md`) is the Phase 6 spec; the older `PLAN_paraview_compaction_2026-04-28.md` carries Phases 0–5. The 2026-05-17 restart work is specified in `debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md` (Phase 1 is the production deliverable; Phases 2–3 are documented follow-ups for bit-exact restart). Cumulative state snapshot: `document/io_dev/STATUS_2026-05-12.md`.

---

## Updates since 2026-05-12 — Restart, cap rework, hardening

This section catalogues the changes that landed between 2026-05-13 and 2026-05-17, in commit order. The earlier walkthrough sections describe steady-state behaviour; this section documents new mechanisms and their justifications.

### Snapshot cap reverted to SOFT (interseismic-only) — 407f456

**What changed.** The R-002 "hard cap" attempt that stretched every regime's interval to fit a fixed snapshot budget was reverted after the round-2 review of `REVIEW.md 2026-05-16` showed it erased coseismic / nucleation detail for any plausible K. The cap is back to the pre-existing soft semantics: `SnapshotCapAwareInterval` returns the **base** (uncapped) cadence whenever `regime != 0`, so only the interseismic regime is throttled.

**How sbatch should use this.** Production runs should size cadences directly using per-regime flags rather than relying on the cap to clip down:

| BP5 driver flag | TPV* driver flag | Default | Recommended |
|---|---|---|---|
| `--paraview-dt-co <s>` | `--paraview-coseismic-dt <s>` | 0.01 s | 1.0 s |
| `--paraview-dt-nu <s>` | `--paraview-nucleation-dt <s>` | 1.0 s | 43200 s (12 h) |
| `--paraview-dt-inter-yr <yr>` | `--paraview-interseismic-dt <s>` | 1.0 yr | 1.0 yr |

`bp5_phase6_paraview_zfp_normal_48hr.sbatch` now passes all three explicitly (~708 frames per 250 yr / 1-event run; ~5,123 per 1775 yr / 7-event run) and does NOT pass `--paraview-max-snapshots` — the natural cadence is trusted.

`bp5_verification_full.cpp` and `tpv*_driver.cpp` use DIFFERENT flag names (BP5 chose `-dt-co/-nu/-inter-yr`, TPV* kept the older long-form names). Both wire through to the same `AdaptiveSchedule::Set{Coseismic,Nucleation,Interseismic}Dt` setters. The follow-up to unify the flag names is queued but not blocking.

**Warning text + rank-0 gate (R-001, R-005, R-203).** When the interseismic budget exhausts, the diagnostic line:

- Now prints **once on rank 0 only** (was N copies, one per rank, torn at character boundaries at np=400).
- Now correctly describes the **geometric-halving** cadence (`t/2, 3t/4, 7t/8, ...`) that the soft cap actually produces, instead of the pre-revert "one final write" phrasing that was inherited from the hard-cap branch.
- Hard-codes "interseismic regime" in the text because the soft cap is structurally only reachable from regime 0 (the `regime != 0` early-return in `SnapshotCapAwareInterval` makes the coseismic / nucleation cases dead code).

The implementation uses `mfem::GetGlobalMPI_Comm()` for the rank query because `mfem::out` wraps `std::cout` unconditionally on parallel builds — the older comment claiming "rank-0-only" was wrong.

### ResetMemory crash on zero-fault-DOF ranks — 407f456 (R-003 / R-204)

On the 8N × 400-rank Frontera BP5 run, ranks whose subdomain contains no fault DOFs construct `Vector state(0)` whose underlying `Memory` is NULL-backed. `PetscODESolver::Run` calls `PetscParVector::ResetMemory` which cannot tolerate that, killing the job mid-loop.

**Driver-side workaround** at `bp5_verification_full.cpp:1740–1747`:

```cpp
const int actual = state.Size();
const int padded = std::max(actual, 1);
state.SetSize(padded);     // ensures a real allocation
state.SetSize(actual);     // shrink back; allocation is preserved
```

The two-step `SetSize` is **load-bearing**: the second call shrinks the logical size but the underlying allocation survives, so `PlaceMemory` aliases a valid pointer and `ResetMemory` finds something to release.

A dedicated unit test `seas_test_bp5_petsc_ts_zero_fault_rank` (8/8 PASS) covers:
- The padded-shrink invariant (R-003).
- A driver-source grep guard (R-204) that fails if the two-call pattern is ever refactored away.

The proper fix lives in MFEM's `PetscParVector::ResetMemory` itself; documented as a deferred upstream PR.

### V2 PETSc-TS restart format — 19e128e (`io/petsc_ts_checkpoint.hpp`)

**Old behaviour.** The driver hard-aborted at `bp5_verification_full.cpp:1091–1099` if you combined `--restart` with `--petsc-ts`. The V1 checkpoint carried MFEM-side state only (slip+psi vector, displacement, traction, slip_rate, FSAL k0) — PETSc TS's internal `t / dt / step / rejections` were absent, so a naive resume would silently lose adaptive-dt history.

**New behaviour.** A `PETSC_TS_V2` block is **appended** to the existing V1 file (no V1 reformat — V1-only files still parse cleanly). Ten fields in fixed order:

```
PETSC_TS_V2
petsc_ts_time                    <real_t>
petsc_ts_dt_next                 <real_t>
petsc_ts_step                    <int>
petsc_ts_rejections              <int>   (cumulative across restart chain — R-005)
paraview_snapshots               <int>   (R-304)
paraview_last_write_time         <real_t>  (R-304)
paraview_last_v_max              <real_t>  (R-304)
paraview_current_regime          <int>   (R-304; clamped to [0,2] on read — R-007)
paraview_last_committed_cycle    <int>   (R-004 dedup key; negative → INT_MIN sentinel)
paraview_last_volume_write_time  <real_t>  (R-006 independent volume-PV cadence)
```

The five `paraview_*` fields are what makes restart work correctly with ParaView output enabled — without them the first `ShouldWrite` after restart fires unconditionally (because the default `last_write_time_ = -1e30` makes `time - last_write_time_` look like +infinity), and the regime FSM resets to interseismic regardless of where the pre-checkpoint trajectory was. R-006 specifically protects the independent `--volume-pv-dt` cadence.

**Driver wiring sites** (`bp5_verification_full.cpp`):

1. **Read at restart** (lines 2583–2691): the V2 block is read AFTER the V1 block — ordering is load-bearing because the V1→V2 cross-check `std::abs(t - ts_t) < 1e-12 * max(|t|,1)` (R-003 round-5 scale-aware tolerance) is only meaningful once `t` has been populated by V1. Updates `t`, `current_dt` BEFORE `PetscODESolver::Run()` because `Run` overwrites `TSSetTime`/`TSSetTimeStep` internally (R-002); `TSSetStepNumber(ts, ts_step)` survives `Run` and is set explicitly. R-008 dt-survival guard: `v2_authoritative_dt` sentinel + bit-exact `MFEM_VERIFY` at the `Run()` call site catches any future CFL clamp that would silently re-introduce the R-002 failure mode.

2. **Write at checkpoint sites**: two sites, both calling `WritePetscTSCheckpoint` immediately after `WriteCheckpoint` (V1):
   - Monitor callback (lines 643–679) — installed only when `--petsc-ts` is active. Sees only `void* ctx → BP5MonitorCtx*`, so `BP5MonitorCtx` was extended (lines 519–523) with two members `pv_out*` and `restart_rejections_carryover` that are wired up in `main()` before `TSMonitorSet`.
   - Final checkpoint in `main` (lines 3020–3058) — same fields, sourced directly from in-scope locals.

3. **Cumulative rejection accumulator (R-005)**: `restart_rejections_carryover` plumbs the prior-chain rejection count into every `WritePetscTSCheckpoint` call site (`cum_rejects = carryover + TSGetStepRejections(ts)`), so the count survives an arbitrary-length restart chain.

4. **`dt_next <= 0` fallback (R-003)**: at seam edges PETSc can return 0 from `TSGetTimeStep`; the driver falls back to `dt_init` rather than passing 0 into `TSSetTimeStep`.

**Caller contract on `ReadPetscTSCheckpoint`.** Returns `true` when a V2 block was parsed; returns `false` (out-params untouched) when the file is V1-only. Driver MUST treat the V1-only case as a hard error when `--petsc-ts` is active — V1-only fallback ("restart MFEM state but let PETSc TS start fresh from dt_init") produces a misleading trajectory and is rejected via `MFEM_VERIFY`.

**`ParaViewOutput` accessors used by the restart machinery** (paraview_output.hpp:418–501):

| Setter (restart-time only, before any Save/ShouldWrite) | Getter (checkpoint-time) |
|---|---|
| `SetTotalSnapshotsWritten(int n)` — clamps negative → 0 | `GetTotalSnapshotsWritten()` |
| `RestoreScheduleState(t_last, V_max_last, regime, t_last_vol)` — clamps regime to [0,2] | `GetLastWriteTime() / GetLastVMax() / GetCurrentRegime() / GetLastVolumeWriteTime()` |
| `SetLastCommittedCycle(int)` — negative input clamps to `INT_MIN` sentinel | `GetLastCommittedCycle()` |

Cluster validation: `bp5_restart_test_v2_dev_2hr.sbatch` runs a two-phase fresh→restart inside a single 2 h development reservation (Phase A ~12.7 yr → V1+V2 checkpoint → Phase B restart to ~25.4 yr), with five built-in acceptance checks (V1+V2 coverage, V2-block log line, t advanced, R-005 monotonic, chainable Phase-B checkpoint) and distinct exit codes 1–7 for diagnosis. ParaView is enabled at production cadence on purpose — the V2 schedule-state fields exist precisely because restart breaks the adaptive schedule without them, and `--no-paraview` would skip every `RestoreScheduleState` / `SetTotalSnapshotsWritten` / `SetLastCommittedCycle` call site behind the `if (pv_out)` guard.

### TPV104 V1 checkpoint — 63373f8 (`io/tpv104_checkpoint.hpp`)

A separate, TPV104-only checkpoint format. Distinct magic tag `TPV104_CHECKPOINT_V1` so cross-driver restarts fail loudly at the header check instead of corrupting state. Per-rank files, one file per MPI rank, reusing BP5's `CheckpointFilename` helper for the path scheme.

Round-trips:
- Wave-field state `Q` (bulk fluctuation, `NUM_STATE * ndof_total`).
- 9 dynamic DOFData fields per fault DOF: `psi, slip_rate, V1, V2, slip1, slip2, tau1_nuc, tau2_nuc, sigma_n_nuc`.
- Scalars: `time, dt, step`.

Two public overloads on each of read/write — one takes `MPIContext*` (BP5 pattern), the other takes raw `(rank, size, MPI_Comm)` (TPV104's pattern, because the TPV104 driver doesn't define `SEAS_USE_MPI` globally). Both forward into an `internal::WriteImpl` / `ReadImpl` body so R-002 + R-104's `expected_Q_size` enforcement lives in exactly one place — the sentinel "skip the check" mode was removed in R-104 to close the silent-bypass hole.

**Out of scope for V1 (deferred to V3):**
- ParaView schedule state for BOTH `pv_out` AND `pv_bulk_out` collections (V2 carries the PRIMARY collection only; TPV104's V1 carries no PV state yet).
- DOFData STATIC fields (impedances, a, Dc, prestress, LSW params) — re-initialised by `InitializeFaultDOFs_TPV104` and don't need to round-trip.
- DOFData corrected-traction fields (`tau1_corr, tau2_corr, sigma_n_corr`) — recomputed by `FaultFaceFlux::Evaluate` on the first post-restart step.

Driver wiring: `tpv104_driver.cpp` reads at line 2355 and writes at lines 2647 / 2671 (intermediate via `--checkpoint-interval` and final-of-run). Restart is mutually exclusive with `--paraview`'s schedule replay — a `WARNING: --restart resumes the PRIMARY ParaView collection's ...` line at line 2376 alerts the operator that the TPV104 V1 checkpoint does not carry PV schedule state.

Unit test: `seas_test_tpv104_checkpoint` (699 lines) covers magic-tag mismatch, rank-count mismatch, Q-size mismatch (R-002 / R-104), and DOFData round-trip.

### Per-collection field-set introspection — `GetVolumeDataCollection()`

A new read-only accessor on `seas::ParaViewOutput`:

```cpp
const ParaViewDataCollectionBase *GetVolumeDataCollection() const;
```

Returns the underlying `pv_dc_` (either `ParaViewDataCollection` or `ParaViewHDFDataCollection`). The motivating user is `seas_test_kinematics_field_set`, which calls `HasField("name")` on the returned handle to verify that the 12 L2-p0 fault projections are NOT registered with the kinematics collection by default and ARE registered when `SetRegisterFaultProjectionsInVolumePV(true)` is set before `InitFaultOutputBP5`. Returns `nullptr` only in the currently-unreachable case where the constructor failed to allocate. Not intended for external write access.

### Debug diagnostic outputs gated behind env vars

The face-trace and fault-DOF-coords diagnostics that production sbatch was leaving on the scratch dir now require explicit env-var opt-in:

| Env var | What it gates | Default |
|---|---|---|
| `SEAS_DEBUG_FACE_TRACE` | The per-rank `[face-trace] ...` log lines + the per-cycle trace CSV | OFF |
| `SEAS_DEBUG_FAULT_DOF_COORDS` | The one-shot `fault_dof_coords_r<rank>.vtp` dump | OFF |

Both are independent (a user who needs the VTP for visualization sanity-check can enable only that one). Earlier revisions had `fault_dof_coords` piggy-back on either gate, which meant the trace-debug bundle silently emitted a 50 MB-per-rank VTP that nobody read. The two-env-var design is documented at `bp5_verification_full.cpp:2153–2174`.

### Build glue — `-lstdc++fs` auto-detect (6cbbb9d)

The BP5 + TPV104 driver `--restart` / `--output-dir` safety check includes `<filesystem>` and calls `std::filesystem::weakly_canonical`. On Frontera's gcc-8.3 + libstdc++ combination, `weakly_canonical` lives in the separate library `-lstdc++fs`; on gcc-9+ it's in the main library and `-lstdc++fs` is a no-op. The Makefile (around line 360) now runs a tiny compile probe at make time and conditionally adds `-lstdc++fs` to `SEAS_FS_LIB`. Operators can override with `make seas_tpv104_driver SEAS_FS_LIB=-lstdc++fs`.

---

## Conventions

- **MeshType templating.** Every method that touches mesh data is a template on `MeshType` (Mesh / ParMesh) via the `GFType` trait. Serial path is the parallel path with all MPI calls compiled away.
- **Collective vs rank-0.** `Save()`, `ShouldWrite()`, `CommitSchedule()`, `WriteFaultSurfaceVTU()` are all collective: every rank must call. The fault file content is rank-0-only but reaching the rank-0 write requires a collective gather.
- **No owning pointers.** Mesh, fault face arrays, DOFData vectors are all passed by reference and not owned. The wrapper owns only its own GFs (`fault_*_`), the underlying data collections (`pv_dc_`), and the H5Z-ZFP probe state.
- **CLI parsing.** Drivers use ad-hoc `HasFlag`/`GetIntArg`/`GetRealArg`/`GetStringArg` helpers (defined in each `tpv*_driver.cpp`). No `argparse`-style library; the flag names are kept identical across all four drivers by manual mirroring.
- **Error handling.** `MFEM_VERIFY` (always-on) for invariants; `MFEM_ABORT` for unreachable branches. `mfem::out` for rank-0 prints (no-op on non-root); `std::cout` for per-rank.
- **Filter locking.** `SetFaultHDFCompression` / `SetVolumeHDFCompression` MUST be called before the first `Save()` on the corresponding collection — the filter selection is read by `VTKHDF::EnsureDataset` exactly once when the first dataset is created. Later calls are accepted but only take effect on **new** datasets (which there generally aren't, since the file is built up entirely on the first save).

---

## Gotchas

- **Two fault representations coexist.** `WriteFaultSurfaceVTU` emits **one triangle per fault face** in `fault_surface.vtkhdf`. `UpdateFaultFieldsBP5` + the primary volume save emit **L2-p0 cell data on the volume mesh** (one DOF per *volume* element, scattered to elements adjacent to fault faces). The fault file is the science output; the volume's `fault_*` fields are debug rendering. They will differ near shared-face boundaries because L2-p0 averages across faces.

- **(RESOLVED 2026-05-12)** Earlier versions of the codebase wrote `velocity` and `mpi_rank` into BOTH `volume.vtkhdf` and `wave_bulk.vtkhdf` — byte-for-byte duplication when cadences/tolerances matched (which they did in the Phase 6 sbatch defaults). The split-bulk-solutions refactor (see `document/io_dev/PLAN_split_bulk_solutions_2026-05-12.md`) fixed this:
  - `kinematics.vtkhdf` (renamed from `volume.vtkhdf`) holds velocity + mpi_rank (no stresses, no fault projections).
  - `stress.vtkhdf` (renamed from `wave_bulk.vtkhdf`) holds the full symmetric stress tensor (6 components: sigma_xx, _yy, _zz, _xy, _xz, _yz) and **no** velocity / mpi_rank.

  Consequence for downstream users: if you open just `stress.vtkhdf` in ParaView to view a stress wavefield, there's no longer a velocity field in the same view — you'd need to also open `kinematics.vtkhdf` and use *Group Datasets* if you want both. This is by design: the duplication cost wasn't worth the convenience.

- **L2-p0 fault GFs need `nbf_per_face > 0`.** `InitFaultOutputBP5` writes per-face averages `inv_nbf = 1.0 / nbf_per_face_`; `nbf=0` would silently produce NaN CellData. R-007 added a hard `MFEM_VERIFY` at the call boundary.

- **`fix_orientation=true` in the rank-0 fault mesh** (fault_vtkhdf_writer.hpp:169) — gmsh may reorder triangle vertices to enforce a positive Jacobian. This is safe for **scalar** cell data (the only thing this writer currently emits). If anyone adds vector cell data or per-vertex data, the gathered `(3i, 3i+1, 3i+2)` ordering will no longer match the rendered vertex order (R-111 note inline).

- **`PeekShouldWrite` then `CommitSchedule(time)` (single-arg) is only correct with `hysteresis_factor = 1.0`.** If hysteresis is enabled the regime FSM is not memoryless in V, so the single-arg shim uses the **last** V_max recorded by `Save`/`ShouldWrite` rather than the *current* V_max. Hysteresis-enabled callers MUST use `CommitSchedule(time, V_max)`.

- **`current_regime_` advancing only when something writes.** If a fault-only path skips `pv_->Save()` and forgets `CommitSchedule`, `current_regime_` pins at 0 and hysteresis silently no-ops. The two-arg `CommitSchedule(time, V_max)` is the cure.

- **`max_total_snapshots` requires `SetTotalRunTime(tfinal)`.** Without it, the cap is silently a no-op. `SnapshotCapAwareInterval` aborts with a readable message rather than degrading.

- **`--paraview-bulk-*` flag scope changed in Phase 6.4.** Pre-Phase-6 those flags applied to the primary `pv_out` (with a one-time rank-0 warning that they didn't work on the VTU path). Post-Phase-6 they route to the **secondary** `pv_bulk_out` collection. BP5 has no secondary, so `--paraview-bulk-*` on BP5 is now a one-shot warning that the flag has no effect.

- **Default volume back end is HDF5, even when you didn't ask.** On any build with `MFEM_USE_HDF5=YES` + (`MFEM_PARALLEL_HDF5=YES` for ParMesh), the volume collection writes `volume.vtkhdf`, not a per-rank VTU. Downstream scripts that grep `<prefix>/ParaView/<basename>_<rank>_<cycle>.vtu` paths break silently. Pass `--paraview-volume-vtu` to recover the legacy layout.

- **`MFEM_USE_H5Z_ZFP=YES` was missing from `config/config.mk.in`** (now fixed). MFEM's build did `#define MFEM_USE_H5Z_ZFP` in `_config.hpp` but `config.mk` was missing the corresponding `MFEM_USE_H5Z_ZFP = YES` line. The Phase 6 sbatch had a `grep -q '^MFEM_USE_H5Z_ZFP *= YES' config.mk || exit` gate that therefore failed on every run. Fixed in `config.mk.in` to add the placeholder.

- **`HDF5_PLUGIN_PATH` is required at run time.** The H5Z-ZFP plugin is not linked into MFEM; HDF5 loads it via this env var. If it's unset or wrong, the first ZFP write fails inside HDF5's filter pipeline (cryptic "required filter is not registered"). The wrapper probes early and aborts with an actionable message.

- **Filter selection is locked at first save.** Calling `SetVolumeHDFCompression(ZfpAccuracy, 1e-3)` after `pv_out->Save()` has fired once is silently a no-op for the existing datasets. The wrapper makes the call eagerly in the constructor path; if you ever add a "change compression mid-run" flow, the only way to honour it is to close + reopen the file with a new MFEM `ParaViewHDFDataCollection`.

- **Mid-run mesh swap is unsupported on the Hdf5 path.** The `.vtkhdf` file is keyed on the first mesh handed to the writer; rewiring `pv_dc_->SetMesh(new_mesh)` would invalidate the dataset shapes. Document at paraview_output.hpp:384–386.

- **Reruns overwrite `.vtkhdf` files.** The V2 restart machinery (2026-05-17) preserves the ParaView **schedule state** across restart (snapshot count, last-write times, regime, dedup cycle) so the cap budget and adaptive cadence survive, but the `.vtkhdf` **file itself** is still opened with `H5F_ACC_TRUNC` on every run — the post-restart write starts a fresh time series, it does not append to the pre-checkpoint one. Operators currently rename the pre-restart `kinematics.vtkhdf` / `stress.vtkhdf` / `fault.vtkhdf` files before restart and concatenate the two series in post. Single-file append-on-restart for VTKHDF is a planned follow-up.

- **V2 restart accessors must be called BEFORE any Save / ShouldWrite / ForceSave.** `SetTotalSnapshotsWritten`, `RestoreScheduleState`, and `SetLastCommittedCycle` set the internal `last_write_time_ / last_v_max_ / current_regime_ / last_committed_cycle_ / last_volume_write_time_` directly; calling them after the schedule has already advanced silently discards the just-computed state and re-runs the regime FSM from a stale baseline. The driver wires them at the V2-restart block immediately after `ReadPetscTSCheckpoint` succeeds, before entering the time loop.

- **V1-only file + `--petsc-ts` + `--restart` is now a hard error**, not a silent degradation. `ReadPetscTSCheckpoint` returns `false` and the driver's `MFEM_VERIFY` aborts with a message pointing the operator at either (a) re-running fresh, or (b) dropping `--petsc-ts` to use MFEM's RK45 (which doesn't need TS-internal state). Don't try to "salvage" by editing the checkpoint by hand — the V2 block expects 10 fields in fixed order.

- **TPV104's V1 checkpoint format is distinct from BP5's V1.** Magic tags are `TPV104_CHECKPOINT_V1` vs `SEAS_CHECKPOINT_V1`. Reading a BP5 checkpoint with the TPV104 reader (or vice versa) aborts at the first `read_tag` line — the failure is loud, not silent, but the error message points at the magic-tag mismatch rather than at the operator's cross-driver `--restart` invocation. Check the `prefix` you passed.

- **Hard cap is gone — soft cap throttles only the interseismic regime.** `--paraview-max-snapshots K` used to (briefly, in an interim R-002 attempt) stretch every regime's cadence to fit K. That was reverted because it sub-sampled the coseismic frames that earthquake animations exist to show. Today K only throttles interseismic; coseismic + nucleation always use their natural cadence. **Operators must size the natural cadences directly via `--paraview-dt-co / -dt-nu / -dt-inter-yr` (BP5) or `--paraview-{co,nucleation,inter}seismic-dt` (TPV*).** If the natural-cadence count is bigger than K, fix the cadences — the cap will NOT save you.

- **The BP5 driver `--paraview-dt-co / -nu / -inter-yr` flags differ in name from the TPV* `--paraview-{co,nucleation,inter}seismic-dt` flags.** They wire through to the same `AdaptiveSchedule` setters; only the spelling differs. A copy-paste between BP5 and TPV* sbatch will silently parse as "unknown flag" and the schedule will keep its defaults. A flag-name unification is queued but not blocking.

- **Zero-fault-DOF ranks need the padded-shrink workaround.** `bp5_verification_full.cpp:1740–1747` does `state.SetSize(max(actual,1))` then `state.SetSize(actual)` so the allocation survives the second shrink. Don't refactor this away — `seas_test_bp5_petsc_ts_zero_fault_rank` greps the driver source for the two-call pattern and fails if it disappears. The proper fix lives in MFEM's `PetscParVector::ResetMemory` upstream.

- **Cap-exhausted warning is silent on non-zero ranks.** R-001's rank-0 gate means only rank 0 prints the "max_total_snapshots= K exhausted" line. On a 400-rank run that's the right behaviour; on a 1-rank serial run nothing changes. But if your grep filter only matches the line text and you forget that it's per-run not per-rank, you'll think it was emitted "once" when actually the latch (`cap_exhausted_warned_` on each rank) guarantees it's emitted at most once per run **per rank**, and rank 0 is the only rank that actually writes.

- **Debug-CSV outputs are now opt-in.** `[face-trace]` log lines and `fault_dof_coords_r<rank>.vtp` files no longer appear unless `SEAS_DEBUG_FACE_TRACE=1` / `SEAS_DEBUG_FAULT_DOF_COORDS=1` are set in the sbatch's `export` block. Old sbatch that relied on the pre-2026-05-16 default-ON behaviour will produce a cleaner results dir than they did before — which is correct, but unexpected.

---

## Open Questions

- **`pv_bulk_out` calls `Save()`-like methods directly via `ForceSave`**, bypassing `WriteFaultSurfaceVTU`. There is no equivalent `SetVolumeSaveEnabled(false)` path for the secondary collection — you'd have to not construct it. Is the asymmetry intentional? (My guess: yes, because the bulk is "volume only" by definition. Worth a comment in the driver.)

- **Phase 5 (post-run packaging) is not yet shipped** per the plan doc. There's no tool that takes an already-emitted multi-million-file directory and re-emits it as VTKHDF. Users with pre-Phase-2 output have to re-run.

- **`SetTotalRunTime(tfinal)` requirement isn't enforced at construction time.** It only matters once `max_total_snapshots > 0`, and only at the moment `SnapshotCapAwareInterval` first runs. If you set `max_total_snapshots` AND forget `SetTotalRunTime`, the abort fires inside the time loop, not at setup. A constructor-time invariant check would be friendlier.

- **The bulk collection inherits the same `VolumeOutputMode` as the primary.** The TPV drivers pin `pv_bulk_out`'s mode to whatever `pv_out` chose. If you want a per-collection back-end choice (e.g. lossless deflate on bulk velocity, lossy ZFP on bulk stresses), the API doesn't currently allow that — you'd need to split the bulk into two more collections.

- **`output_mode_uniformity_checked_` only covers `legacy_ascii_ + output_mode_`.** Other state that must be uniform across ranks (filter algorithm, ZFP tolerance, snapshot cap, the new `register_fault_projections_in_volume_pv_` toggle, the V2 restart accessors) is not checked. In practice these are set from CLI / V2-block parsing which is uniform, but the assumption isn't enforced. A rank-divergent `SetRegisterFaultProjectionsInVolumePV` would deadlock on the gather collective inside `WriteFaultSurfaceVTU`.

- **`MFEM_USE_H5Z_ZFP` doesn't actually gate any link-time code** — it's pure preprocessor, used by `seas::ParaViewOutput::ProbeH5ZZfpPluginOrAbort` to know whether to probe. The flag exists mostly so the sbatch grep gate can fail fast. Could be inlined.

### Restart-specific open questions

- **VTKHDF append-on-restart is not implemented.** The V2 checkpoint preserves the ParaView **schedule** state but the `.vtkhdf` files themselves are truncated on every fresh run; restart writes a new series rather than appending to the pre-checkpoint one. Operators currently post-concatenate. A clean fix would be to detect "restart mode" in `EnsureVTKHDF` and open with `H5F_ACC_RDWR` + seek-to-end on the unlimited datasets, but `Steps/PartOffsets` and the global `NumberOfPoints / NumberOfCells / NumberOfConnectivityIds` accumulators would all need new "resume" code paths.

- **TPV104 V2 (PV-state extension) is not designed yet.** TPV104's V1 carries no `ParaViewOutput` schedule state — restart loses both the snapshot count and the regime FSM. BP5's V2 trailing block carries only the PRIMARY collection's state; even when ported to TPV104, the SECONDARY `pv_bulk_out` collection's `last_volume_write_time_` for stress fields would need a parallel set of fields (or a single block that loops over collections). Plan: V3 layout to be designed.

- **The `restart_rejections_carryover` accumulator (R-005) is BP5-only.** TPV104 uses ADER-DG explicit substepping (`tpv104_substep_iterator.cpp`) which has no rejection concept, so this is fine for TPV104 specifically. But if a future PETSc-TS-based explicit-RK TPV driver appears, it'll need the same accumulator pattern.

- **The V1 ↔ V2 cross-check uses `1e-12 * max(|t|, 1)` (R-003 round 5).** For runs near `t = 0` (a fresh restart of an early-aborted run, say) the `max(|t|, 1)` floor gives a fixed tolerance of `1e-12 s`, which is below dt for any realistic SEAS problem. For runs at `t ~ 1e12 s` (millions of years), the relative tolerance scales appropriately. But the analysis assumes the V1 write and V2 write happen in the **same** SLURM step — if they're separated by a `Barrier` + collective write that touches `t` between them, the check could spuriously fail. Currently the two writes are bracketed tightly enough that this hasn't been observed, but it's not enforced.

- **`bp5_restart_test_v2_dev_2hr.sbatch` has not yet been ported to a production-queue 48h reservation.** The dev-queue 2h reservation is sized for the two-phase fresh→restart proof of concept (12.7 yr + 12.7 yr); a real long-run scenario (250 yr → checkpoint → 1775 yr) needs a `normal`-queue analogue with the V1+V2 file management (the V2 trailing block from Phase A becomes Phase B's restart prefix, etc.).
