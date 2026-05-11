# Implementation Plan: Bulk Volume Compression and Output Size Estimator

**Date:** 2026-05-09
**Status:** Draft for review.  No code changes proposed in this document.
**Relationship to prior plan:** This document extends
`PLAN_paraview_compaction_2026-04-28.md` with two phases — **Phase 6**
(bulk volume PV → VTKHDF + ZFP) and **Phase 7** (pre-submit output size
estimator).  The earlier plan's Phases 0–5 are presumed complete and
its conventions (file naming, MFEM_USE_HDF5 / MFEM_USE_H5Z_ZFP guards,
driver wiring style) are inherited verbatim.

**Owners:**
- Phase 6 — `miniapps/seas/io/paraview_output.hpp` (primary);
  drivers `tpv102_driver.cpp`, `tpv104_driver.cpp`, `tpv205_driver.cpp`,
  `tests/verification/bp5_verification_full.cpp` (CLI wiring); selected
  unit tests under `tests/unit/`.
- Phase 7 — `miniapps/seas/scripts/estimate_output_size.py` (new file);
  optional helper data in `miniapps/seas/scripts/_io_size_*.py`.

---

## Overview

Two follow-ups to the io_dev compaction work:

1. **Phase 6** — the *fault* path is single-file VTKHDF with optional
   ZFP after Phase 2/2d.  The *volume* path is still the legacy
   `mfem::ParaViewDataCollection` (per-rank ASCII/binary VTU + PVTU +
   PVD) — exactly the layout Phase 1 was written to eliminate.  For
   any run that wants both fault and volume output, the volume side
   reproduces the original Frontera pain (millions of files in
   `output_dir/ParaView/` and `output_dir/ParaView_bulk/` for tpv*
   drivers).  Phase 6 promotes the volume path to
   `ParaViewHDFDataCollection`, exposes the bulk-side ZFP / deflate
   knobs that Phase 2d.3 already parses-but-warns, and keeps the
   legacy VTU path as an opt-in fallback for ParaView < 5.11
   compatibility.

2. **Phase 7** — even after Phase 6, a 250-yr BP5 run on 800 ranks at
   1000 m mesh resolution produces hundreds of GB of compressed HDF.
   Today the only way to know whether a configuration fits in `$SCRATCH`
   is to submit and wait.  Phase 7 ships a pre-submit Python tool that
   reads the mesh + driver parameters, integrates the adaptive schedule
   over `tfinal`, and prints an estimated total output size with a
   per-field breakdown.  No MFEM dependency — runs on a laptop in <5 s.

Quantitative targets (validated by acceptance tests in §Phase 6 / §Phase 7):

| Metric                               | Today          | After Phase 6 (deflate) | After Phase 6 (ZFP)    |
|--------------------------------------|----------------|-------------------------|------------------------|
| Volume PV file count, BP5 250 yr     | ~5e6           | **1**                   | **1**                  |
| Bulk PV bytes/cycle (TPV102, 200 m)  | ~3 GB ASCII    | ~0.6–1 GB               | ~0.05–0.15 GB          |
| `--paraview-bulk-zfp-tol 1e-3` semantics | `warn-and-ignore` | `warn-and-ignore` | **applied** |

For Phase 7:

| Metric                             | Today           | After Phase 7   |
|------------------------------------|-----------------|-----------------|
| Pre-submit size estimate accuracy  | n/a (no tool)   | ≤ 30 % of measured |
| Time to estimate                   | n/a             | < 5 s on laptop |
| User-visible fields                | n/a             | per-field bytes, total, fits-in-$SCRATCH yes/no |

---

## Constraints (inherited from PLAN_paraview_compaction_2026-04-28)

### Interface constraints (hard)

- `seas::ParaViewOutput<MeshType>` public API may extend but **must
  not break**.  Specifically:
  - `RegisterDomainField(name, gf)` keeps its signature; the call
    forwards through the new `pv_dc_` (a
    `ParaViewDataCollectionBase*`) via the virtual
    `DataCollection::RegisterField` (R-311: `RegisterField` is
    declared virtual on the GRANDPARENT `DataCollection` at
    `fem/datacollection.hpp:245`, NOT on `ParaViewDataCollectionBase`).
    Both `ParaViewDataCollection` and `ParaViewHDFDataCollection`
    inherit the `DataCollection::RegisterField` body
    (direct-into-`field_map`); neither overrides it.  Virtual
    dispatch through the base pointer therefore reaches the same
    implementation regardless of which concrete subclass is active.
  - `Save(cycle, time, V_max)`, `ForceSave`, `ShouldWrite`,
    `PeekShouldWrite`, `CommitSchedule(...)` — semantics unchanged.
  - `SetVolumeSaveEnabled`, `SetVolumePVDt`, `SetFaultOutputMode`,
    `SetFaultHDFCompression` — semantics unchanged.
  - `SetVolumeHDFCompression(alg, param)` — currently a one-time
    rank-0 warning; Phase 6 makes it a no-op when the active volume
    writer is VTU AND a real apply when the writer is HDF.

- The on-disk **logical layout** that ParaView sees changes from
  `output_dir/ParaView/fault_surface_*.{vtu,pvtu}` + a PVD to a single
  `output_dir/ParaView/<collection_name>.vtkhdf` (per
  `mfem::ParaViewHDFDataCollection`'s naming convention).  Existing
  workflows that opened the legacy PVD must be migrated; legacy paths
  are preserved via an explicit CLI flag (see §Phase 6.3).

### Dependency constraints

- **Phase 6 requires `MFEM_USE_HDF5 = YES` AND
  `MFEM_PARALLEL_HDF5 = YES`** for the parallel default.  The build
  must use the OpenMPI- or MPICH-flavoured parallel HDF5
  (`mpi_openmpi_*` / `mpi_mpich_*` conda-forge selector, or a TACC
  `phdf5/<version>` module), NOT the serial conda-forge `hdf5`
  package.  See `PLAN_paraview_compaction_2026-04-28.md` §Phase 2a.2
  for the install recipe.
  - **Serial-HDF5 builds** (`MFEM_USE_HDF5=YES, MFEM_PARALLEL_HDF5=NO`):
    `ParaViewHDFDataCollection` aborts at construction time on a
    `ParMesh` per `fem/datacollection.cpp:1414` (`MFEM_ABORT("Requires
    HDF5 library with parallel support enabled")`).  Therefore the
    `DefaultVolumeOutputMode()` helper in §Phase 6.1 falls back to
    `Vtu` when the build is `ParMesh + MFEM_USE_HDF5 +
    !MFEM_PARALLEL_HDF5`.
  - **Non-HDF5 builds** (`MFEM_USE_HDF5=NO`): default is Vtu;
    `--paraview-volume-hdf5` is rejected at parse time.
  - All Phase 6 code is `#ifdef MFEM_USE_HDF5`-guarded so non-HDF5
    builds compile and run.
- **Phase 6 ZFP applied to bulk** requires `MFEM_USE_H5Z_ZFP = YES`
  AND the H5Z-ZFP plugin discoverable via `HDF5_PLUGIN_PATH` at
  runtime.  Same conditions as Phase 2d.3 fault-side ZFP.
- **Phase 7 has zero MFEM dependency** — it is a pure Python tool.
  Required: Python 3.8+, standard library only.  Optional: `pyvista`
  for direct comparison against an existing run's actual output (the
  sanity-check mode in §Phase 7.5).

### Convention constraints

- Follow the io_dev prior-plan conventions: `--paraview-volume-*` flags
  for the new volume controls, mirroring the `--paraview-fault-*`
  flags introduced by Phase 2b/2d.3.  Mutually exclusive validation
  identical to the fault side: `--paraview-volume-zfp-tol` and
  `--paraview-volume-deflate-level` cannot both be set.
- New Python helper modules under `miniapps/seas/scripts/_io_size_*.py`
  use the leading-underscore convention to mark them private to the
  estimator (not user-facing entry points).
- Do not hardcode numerical constants (per `feedback_no_hardcoded_numbers.md`):
  the compression-ratio table, regime cadence defaults, and field
  shapes must all be derived from project headers (`bp5_params.hpp`)
  or driver inspection — not magic numbers in the estimator.

### Numerical / performance constraints

- Phase 6 promotion to HDF must not regress fault-side throughput.
  Volume Save() is the dominant cost on TPV runs; the move to
  parallel HDF5 collective writes should be **at most 10 % slower**
  than the current per-rank VTU writes on np ≥ 100, and faster on
  np ≥ 400 (where Lustre per-file metadata dominates).
- Phase 7 estimator accuracy target: ≤ 30 % overestimate, ≤ 50 %
  underestimate of the actual on-disk total bytes for any of the
  three reference scenarios in §Phase 7.6.  Overestimate-bias is
  preferred since the user planning around `$SCRATCH` quotas would
  rather provision too much.

---

## Phase 6: Bulk Volume PV via VTKHDF + ZFP (production target)

### Goal

After this phase, the volume `pv_` member of
`seas::ParaViewOutput<MeshType>` writes a single VTKHDF file per run
(matching the fault-side Phase 2 layout) on builds where
`MFEM_USE_HDF5=YES`, and `--paraview-volume-{zfp-tol,deflate-level}`
are no longer warn-and-ignore — they apply real chunk filters to bulk
velocity / stress / displacement fields.  The legacy per-rank VTU
path remains accessible via `--paraview-volume-vtu` for ParaView <
5.11 compatibility and for debugging.

### Sub-phase ordering

| Sub-phase | What | Where | Acceptance gate |
|---|---|---|---|
| **6a** | Library refactor — replace concrete `ParaViewDataCollection pv_` member with abstract `std::unique_ptr<ParaViewDataCollectionBase> pv_dc_`; add `VolumeOutputMode` enum; default selection from `MFEM_USE_HDF5`. | `miniapps/seas/io/paraview_output.hpp` | Existing tests pass on `MFEM_USE_HDF5=YES` build with `--paraview-volume-vtu`; same tests pass with default (HDF5) on a 1-yr smoke run. |
| **6b** | Make `SetVolumeHDFCompression` a real apply when active writer is HDF; warn-and-ignore otherwise. | `paraview_output.hpp` | New unit test `test_volume_hdf_compression.cpp` passes both `Deflate` and `ZfpAccuracy` paths; `test_paraview_tolerance_isolation.cpp` from Phase 2d.3 newly enables the bulk half. |
| **6c** | Driver CLI `--paraview-volume-{vtu,hdf5}` selectors. | `tpv102/104/205_driver.cpp`, `bp5_verification_full.cpp` | Each driver's `--help` lists the new flag; parse-time rejection on non-HDF5 builds for `--paraview-volume-hdf5`. |
| **6c+** (a.k.a. **6.3a**, R-301) | Driver CLI `--paraview-volume-{zfp-tol,deflate-level}` for the PRIMARY collection's volume PV; mutual-exclusion + build-flag rejection identical to fault side; routes to `pv_out->SetVolumeHDFCompression`. | `tpv102/104/205_driver.cpp`, `bp5_verification_full.cpp` | Driver shell test: BP5 run with `--paraview-volume-zfp-tol 1e-3` produces a `<output>/ParaView/volume.vtkhdf` whose displacement dataset carries filter id 32013. |
| **6d** | `pv_bulk_out` (secondary `ParaViewOutput` instance in tpv* drivers) inherits the same VolumeOutputMode default; **`--paraview-bulk-*` is RE-ROUTED from `pv_out` (Phase 2d.3 warn-and-ignore) to `pv_bulk_out` (real apply)** — see R-310 migration note. | `tpv102/104/205_driver.cpp` | New test `test_pv_bulk_out_zfp.cpp` writes synthetic velocity to a `pv_bulk_out` with `--paraview-bulk-zfp-tol 1e-3` and asserts filter id 32013 on the velocity dataset of `<output>/ParaView_bulk/wave_bulk.vtkhdf`. |
| **6e** | Test migration — existing tests that opened the legacy `<prefix>/ParaView/<name>_*.vtu` must either opt into Vtu via `--paraview-volume-vtu` OR migrate to read VTKHDF (per the R-306 decision rule). | `tests/unit/test_*paraview*.cpp` | All existing tests pass; no spurious skips. |
| **6f** | Frontera replication; sbatch updates for the production templates. | Frontera-side build + run | One TPV102 microrun (1 s, np=400) produces a single `<output>/ParaView/wave.vtkhdf` < 1 GB at default deflate; opens in ParaView 5.13+. |

---

### Phase 6.1: Library refactor (sub-phase 6a)

#### Files to Create

None.  Pure refactor of an existing file.

#### Files to Modify

- `miniapps/seas/io/paraview_output.hpp`:
  - Add `enum class VolumeOutputMode { Vtu, Hdf5 };` directly under the
    existing `enum class FaultOutputMode { Vtu, Hdf5 };` (line ~80 of
    the current file).
  - Replace the concrete member
    `ParaViewDataCollection pv_;`
    (currently at line 1422) with an abstract handle:
    ```cpp
    std::unique_ptr<ParaViewDataCollectionBase> pv_dc_;
    VolumeOutputMode volume_output_mode_ =
    #ifdef MFEM_USE_HDF5
        VolumeOutputMode::Hdf5;     // default when HDF5 is built in
    #else
        VolumeOutputMode::Vtu;
    #endif
    ```
  - In the constructor (line 313-326), construct the right concrete
    type into `pv_dc_`:
    ```cpp
    if (volume_output_mode_ == VolumeOutputMode::Hdf5)
    {
    #ifdef MFEM_USE_HDF5
       pv_dc_ = std::make_unique<ParaViewHDFDataCollection>(
                  GenerateCollectionName(prefix), &mesh);
    #else
       MFEM_ABORT("VolumeOutputMode::Hdf5 selected on a build without "
                  "MFEM_USE_HDF5=YES");
    #endif
    }
    else
    {
       pv_dc_ = std::make_unique<ParaViewDataCollection>(
                  GenerateCollectionName(prefix), &mesh);
    }
    pv_dc_->SetPrefixPath(prefix);
    pv_dc_->SetDataFormat(VTKFormat::BINARY);
    pv_dc_->SetHighOrderOutput(true);
    pv_dc_->SetLevelsOfDetail(order);
    ```
  - The collection name is supplied as a constructor parameter (NOT
    derived from `prefix`).  R-305 resolves the default to
    `"volume"` (cross-driver uniform, matching the fault path's
    `"fault_surface"` convention).  No `GenerateCollectionName`
    helper is needed; the constructor parameter has a default.
    The final constructor signature (R-303 / R-305): `collection_name`
    AND `mode` appended at the end with defaults.  See the §Interfaces
    block below for the canonical declaration.  No back-compat shim
    is needed because default arguments handle the existing 3-arg
    `(prefix, mesh, order)` call sites without source changes.
    Existing callers get `"volume"` as the file basename automatically.
  - Replace every direct `pv_.X(...)` call site in this file with
    `pv_dc_->X(...)`:
    - `pv_.SetPrefixPath` (line 322), `SetDataFormat` (323),
      `SetHighOrderOutput` (324), `SetLevelsOfDetail` (325)
    - `pv_.RegisterField` (line 335; called from `RegisterDomainField`)
    - All 17 `pv_.RegisterField(...)` calls in `InitFaultOutputBP5`
      (lines ~643-660)
    - `pv_.SetCycle`, `pv_.SetTime`, `pv_.Save()` in `ForceSaveImpl`
      (lines ~1507-1510)
  - Add a public **getter only** (no setter) for the new mode:
    ```cpp
    VolumeOutputMode GetVolumeOutputMode() const { return volume_output_mode_; }
    ```
    R-302: the mode is set EXCLUSIVELY through the constructor
    parameter — no setter exists.  Once `pv_dc_` is built, the mode
    is locked.  This mirrors the immutability of `mesh_` and
    `order_` already in `ParaViewOutput`.  A public setter would be
    unreachable from external callers (by the time a caller has a
    `ParaViewOutput*`, construction is by definition complete) and
    would therefore be dead code.
  - R-303: extend the constructor by APPENDING the new parameters
    at the end (with defaults) so the existing 3-arg call sites in
    8 drivers + 6 unit tests keep compiling without source edits.
    R-304: the default selector consults BOTH `MFEM_USE_HDF5` AND
    `MFEM_PARALLEL_HDF5` for `ParMesh` builds — a serial-HDF5 build
    cannot drive `ParaViewHDFDataCollection` on a `ParMesh` (see
    §Dependency constraints), so it must default to Vtu in that
    case.
    ```cpp
    ParaViewOutput(const std::string &prefix,
                   MeshType &mesh,
                   int order,
                   const std::string &collection_name = "volume",
                   VolumeOutputMode mode = DefaultVolumeOutputMode());

    // Compile-time helper that picks Hdf5 only when the build can
    // actually use it (HDF5 + parallel-HDF5 for ParMesh; HDF5 alone
    // for serial Mesh).  Drops to Vtu otherwise.
    static constexpr VolumeOutputMode DefaultVolumeOutputMode()
    {
    #if defined(MFEM_USE_HDF5)
       if constexpr (std::is_same_v<MeshType, ParMesh>)
       {
    #  if defined(MFEM_PARALLEL_HDF5)
          return VolumeOutputMode::Hdf5;
    #  else
          return VolumeOutputMode::Vtu;
    #  endif
       }
       else
       {
          return VolumeOutputMode::Hdf5;
       }
    #else
       return VolumeOutputMode::Vtu;
    #endif
    }
    ```
    The 3-arg shim is no longer needed (default arguments handle
    both legacy `(prefix, mesh, order)` calls AND new
    `(prefix, mesh, order, "wave", Hdf5)` calls).

#### Detailed Requirements

1. **Casting helper** — many call sites need
   `ParaViewHDFDataCollection`-only methods (`SetCompression(true)`,
   `SetHDFCompression(...)`).  Add a private helper:
   ```cpp
   ParaViewHDFDataCollection* GetHDFDC()
   {
   #ifdef MFEM_USE_HDF5
      return dynamic_cast<ParaViewHDFDataCollection*>(pv_dc_.get());
   #else
      return nullptr;
   #endif
   }
   ```
   Returns `nullptr` if the active collection is VTU; the call site
   checks and either applies the HDF-specific call or skips it (with
   the existing one-time-warning for `SetVolumeHDFCompression`).

2. **`ForceSaveImpl` cadence gate** (currently at line 1482) preserves
   its Phase 4 logic verbatim, except the calls become:
   ```cpp
   pv_dc_->SetCycle(cycle);
   pv_dc_->SetTime(time);
   pv_dc_->Save();
   ```
   The HDF path's `Save()` is collective on `mesh_.GetComm()` for
   `ParMesh`; the VTU path's `Save()` is per-rank with rank-0 PVD
   write.  This change is transparent — both inherit `Save()` from
   `ParaViewDataCollectionBase` (overridden in each subclass).

3. **`SetVolumeSaveEnabled(false)` + cadence skip** (Phase 4 R-001
   logic) continues to work without change because the cadence check
   gates the entire `pv_dc_->Save()` call.

4. **No-op-when-empty invariant** — if the user calls
   `RegisterDomainField` with zero fields and then `Save()`, MFEM's
   HDF collection writes only the mesh and an empty `Steps/` group.
   That's accepted; current VTU behavior is similar (empty CellData
   + PointData blocks).  No new code needed.

#### Interfaces

```cpp
// New in seas::ParaViewOutput<MeshType>
public:
   enum class VolumeOutputMode { Vtu, Hdf5 };

   // R-303: new params APPENDED to the existing 3-arg signature so
   // every existing 3-arg call site (8 drivers + 6 unit tests) keeps
   // compiling unchanged.  No back-compat shim needed.
   ParaViewOutput(const std::string &prefix,
                  MeshType &mesh,
                  int order,
                  const std::string &collection_name = "volume",
                  VolumeOutputMode mode = DefaultVolumeOutputMode());

   // R-302: getter only — no setter.  Mode is locked at construction.
   VolumeOutputMode GetVolumeOutputMode() const;

   // R-304: compile-time default that respects both MFEM_USE_HDF5
   // and MFEM_PARALLEL_HDF5 (for ParMesh).
   static constexpr VolumeOutputMode DefaultVolumeOutputMode();

private:
   std::unique_ptr<ParaViewDataCollectionBase> pv_dc_;
   VolumeOutputMode volume_output_mode_;
   ParaViewHDFDataCollection *GetHDFDC();   // nullptr if active is VTU
```

#### Edge Cases to Handle

- **Build without HDF5**: `VolumeOutputMode::Hdf5` is the default per
  the `#ifdef`, but we never select it on a non-HDF5 build because
  the `#else` flips the default to Vtu.  Defense-in-depth:
  `MFEM_ABORT` in the constructor if Hdf5 is requested without
  HDF5 — should be unreachable but clearly diagnosed if a downstream
  user passes `Hdf5` explicitly.
- **Existing 3-arg constructor**: must continue to compile and link
  without source changes at every call site.  All current call sites
  (8 across drivers + tests, plus 2 in unit tests) use the 3-arg
  form; the shim preserves binary compatibility.
- **Mesh swap mid-run**: `pv_dc_->SetMesh(...)` is inherited from
  `DataCollection`; `ParaViewHDFDataCollection` does not support a
  mid-run mesh swap (the HDF file is keyed on the first mesh).
  Document this in `paraview_output.hpp` as a known limitation; no
  driver currently swaps the mesh, so no test is needed.
- **Empty fault but volume on**: `pv_out->ForceSave(0, 0.0)` must
  still emit a valid VTKHDF file with no fault datasets — only
  domain fields registered via `RegisterDomainField`.
  `ParaViewHDFDataCollection` handles this natively.
- **Re-run with same `output_dir`**: ParaViewHDFDataCollection
  truncates the file by default unless `UseRestartMode(true)` is
  set.  Same as the fault path; document in `--help`.

#### Acceptance Criteria

- [ ] `make seas_test_paraview_volume_mode_smoke` (new test, see
      §Phase 6.5) passes serial AND `mpirun -np 4` on an
      `MFEM_USE_HDF5=YES` build.
- [ ] All existing `tests/unit/test_*` targets still build and pass.
      Specifically `test_io.cpp` (uses the legacy 3-arg constructor)
      passes without modification.
- [ ] On a default `MFEM_USE_HDF5=YES` build, the existing TPV102
      smoke test (`make test-tpv102`) emits ONE `volume.vtkhdf` (or
      whatever the collection name is) instead of N+1 per-cycle VTU
      files.  File-count check: `ls output_dir/ParaView/ | wc -l`
      returns ≤ 2 (the `.vtkhdf` plus optionally a sidecar PVD-like
      JSON).
- [ ] On the SAME build with `--paraview-volume-vtu`, the legacy
      per-rank-VTU layout is emitted (regression gate; > 2 files).

#### Dependencies

- Depends on: Phase 2 (existing) — `ParaViewHDFDataCollection` API.
- Required by: Phase 6.2 (compression wiring), Phase 6.3 (driver CLI).

---

### Phase 6.2: SetVolumeHDFCompression actually applies (sub-phase 6b)

#### Files to Create

- `miniapps/seas/tests/unit/test_volume_hdf_compression.cpp` — verifies
  filter id 32013 (ZFP) and filter id 1 (deflate) actually attach to
  the bulk velocity dataset.  Mirrors `test_vtkhdf_zfp.cpp` from Phase
  2d.2 but exercises the public seas API
  `SetVolumeHDFCompression`.

#### Files to Modify

- `miniapps/seas/io/paraview_output.hpp`:
  - Replace the existing `SetVolumeHDFCompression` body (currently
    lines ~480-525, the warn-and-ignore stub) with a real apply when
    `pv_dc_` is HDF:
    ```cpp
    void SetVolumeHDFCompression(
        ParaViewHDFDataCollection::HDFCompression alg,
        double param)
    {
       if (alg == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy)
       {
          MFEM_VERIFY(param > 0.0, "...");
       }
       auto *hdf_dc = GetHDFDC();
       if (hdf_dc)
       {
          hdf_dc->SetHDFCompression(alg, param);
          if (alg == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy)
          {
             ProbeH5ZZfpPluginOrAbort();   // existing R-307 helper
          }
       }
       else
       {
          // Fallback: VTU active — emit the existing warn-and-ignore
          // message so the user knows.
          if (!volume_hdf_compression_warned_) { mfem::out << "..."; }
          volume_hdf_compression_warned_ = true;
       }
    }
    ```
  - The `volume_hdf_compression_warned_` member (currently `mutable bool`)
    is unchanged.

#### Detailed Requirements

1. **Per-dataset filter check** — the new test must verify that
   filter 32013 (ZFP) actually attaches to the velocity dataset in
   the volume `.vtkhdf`, NOT just to mesh / connectivity (which is
   integer and falls back to deflate per the
   `mesh/vtkhdf.cpp:EnsureDataset` dispatch).  Use the same
   `ReadFilterIds(file, "/VTKHDF/PointData/velocity")` helper
   pattern from `test_vtkhdf_zfp.cpp`.

2. **Per-dataset accuracy check** — write a smooth analytic velocity
   field (e.g. `v(x,y,z) = sin(2πx)cos(πy)e^(-z)`), Save with
   `--paraview-volume-zfp-tol 1e-3`, read back with raw HDF5 / pyvista,
   assert `max_abs_err ≤ 1e-3` AND compression-ratio `≥ 4×` vs the
   same field with default deflate.

3. **Driver-side acceptance** — once the bulk side actually applies,
   the driver-side parse-time validation (already present in
   tpv102/104/205) starts being meaningful.  Verify on a
   `MFEM_USE_HDF5=NO` build that `--paraview-volume-zfp-tol 1e-3`
   aborts at parse time with the documented error string.

#### Interfaces

```cpp
// In seas::ParaViewOutput<MeshType> — body changes only; signature
// is unchanged from Phase 2d.3.
void SetVolumeHDFCompression(
    ParaViewHDFDataCollection::HDFCompression alg,
    double param);
```

#### Edge Cases to Handle

- **Active writer is VTU but the user asked for ZFP**: print the
  existing one-time warning AND keep `volume_hdf_compression_warned_`
  latched.  Bulk output proceeds as plain VTU.  No abort.
- **Active writer is HDF but the user asked for `ZfpAccuracy`**
  AND `MFEM_USE_H5Z_ZFP=NO`: existing `MFEM_VERIFY` in
  `VTKHDF::SetZfpAccuracy` aborts with a clear message; driver-side
  parse-time validation should have caught this earlier, but
  defense-in-depth.
- **`SetVolumeHDFCompression` called after the first `Save()`**:
  the Phase 2d.2 MFEM patch only consults `algorithm` /
  `compression_level` at `EnsureDataset` (i.e. the FIRST write of a
  given dataset).  Subsequent writes append to the existing dataset's
  filter chain.  Net effect: changing the algorithm mid-run only
  affects newly-created datasets (e.g. a field registered after the
  first Save).  Document this in the doc-comment.
- **VTU active + `ZfpAccuracy` requested** on a build with
  `MFEM_USE_H5Z_ZFP=YES`: still warn-and-ignore.  The user's intent
  was to use ZFP; the active writer can't.  No abort.

#### Acceptance Criteria

- [ ] `make test-volume-hdf-compression` passes (new target):
      - filter id 32013 attached to `velocity` dataset under ZFP mode;
      - max_abs_err ≤ 1e-3 vs ground truth at `tol=1e-3`;
      - per-dataset compression ratio ≥ 4× vs default deflate at
        same input.
- [ ] `make test-paraview-tolerance-isolation` (Phase 2d.3 acceptance,
      previously bulk-half-deferred) passes both halves now: bulk
      `tol=1e-3` applies to velocity, fault `tol=1e-12` applies to
      slip_rate_strike, ratios differ by > 10× (proves no
      cross-talk).
- [ ] `make test-vtkhdf-zfp test-fault-surface-vtkhdf-zfp` still pass
      (regression — Phase 2d.2/2d.3 paths unchanged by Phase 6.2).
- [ ] Driver shell test on a non-H5Z-ZFP build: `--paraview-volume-zfp-tol
      1e-3` aborts at parse time.

#### Dependencies

- Depends on: Phase 6.1 (HDF-active volume writer).
- Required by: Phase 6.6 (Frontera deploy).

---

### Phase 6.3: Driver CLI volume-mode selectors (sub-phase 6c)

#### Files to Create

None.

#### Files to Modify

- `miniapps/seas/drivers/tpv102_driver.cpp`:
  - Add CLI flags parsed alongside the existing `--paraview-fault-vtu` /
    `--paraview-fault-hdf5` block (currently lines ~519-545 of the
    file):
    - `--paraview-volume-vtu` — force `VolumeOutputMode::Vtu`
    - `--paraview-volume-hdf5` — force `VolumeOutputMode::Hdf5` (the
      new default; explicit version for sbatch readability)
  - Parse-time validation:
    - `--paraview-volume-vtu` AND `--paraview-volume-hdf5` together
      → `MFEM_ABORT("--paraview-volume-vtu and --paraview-volume-hdf5
      are mutually exclusive.")`
    - `--paraview-volume-hdf5` on `#ifndef MFEM_USE_HDF5` →
      `MFEM_ABORT("--paraview-volume-hdf5 requires MFEM_USE_HDF5=YES;
      current build has it disabled.")`
  - Wire into the `ParaViewOutput` constructor:
    ```cpp
    auto volume_mode =
    #ifdef MFEM_USE_HDF5
        seas::ParaViewOutput<MeshT>::VolumeOutputMode::Hdf5;
    #else
        seas::ParaViewOutput<MeshT>::VolumeOutputMode::Vtu;
    #endif
    if (paraview_volume_force_vtu)  { volume_mode = ...::Vtu; }
    if (paraview_volume_force_hdf5) { volume_mode = ...::Hdf5; }

    // R-303 / R-305: append-at-end constructor signature with
    // collection_name = "volume" (the cross-driver uniform default).
    pv_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
       output_dir + "/ParaView", pmesh, order,
       /*collection_name=*/"volume", volume_mode);
    ```
  - Update the `--help` block (the comment-block at line ~480 that
    lists CLI flags) to include the new flags.
- `tpv104_driver.cpp`, `tpv205_driver.cpp`: identical changes —
  cross-driver-uniform collection name `"volume"`.
- `tests/verification/bp5_verification_full.cpp`: same treatment with
  collection name `"volume"`.  The driver's existing hand-rolled CLI
  parser (`for (int i = 1; i < argc; ...)`) gets the new flags at the
  same site as the Phase-2b/3/4 flags wired earlier.

#### Detailed Requirements

1. **Default is HDF5** when the build defines `MFEM_USE_HDF5`.  Users
   on Frontera with the post-Phase-6f sbatch templates get VTKHDF
   automatically; they don't need to pass `--paraview-volume-hdf5`.
2. **No silent migration** — when the default flips, the run banner
   prints the active mode:
   ```
   Volume PV: ON (HDF5; output_dir/ParaView/wave.vtkhdf)
   ```
   or:
   ```
   Volume PV: ON (VTU; output_dir/ParaView/wave_*.vtu)
   ```
3. **The `pv_no_domain` / `--no-volume-pv` Phase-4 path is
   independent** — `--no-volume-pv` suppresses ALL volume saves
   regardless of mode.  The new `--paraview-volume-{vtu,hdf5}` only
   matters when volume save is enabled.

#### Interfaces

(No new C++ interfaces; this is driver CLI wiring only.)

```bash
# tpv102 / tpv104 / tpv205 / seas_bp5_full now accept:
--paraview-volume-vtu           # legacy per-rank VTU
--paraview-volume-hdf5          # single VTKHDF (default on HDF5 builds)
```

#### Edge Cases to Handle

- **`pv_bulk_out` (secondary collection in tpv*) inherits the
  selector** — see §Phase 6.4.  This phase only addresses the
  primary `pv_out`.
- **Existing sbatch scripts that DO NOT pass `--paraview-volume-*`**
  silently switch to HDF5 on `MFEM_USE_HDF5=YES` builds.  Plan §Phase
  2b's exact same rationale applies: this is the production
  default, but legacy workflows opening per-rank VTU directly will
  break.  Documented in the run banner (see step 2 above).

#### Acceptance Criteria

- [ ] `--paraview-volume-vtu`, `--paraview-volume-hdf5`, and the
      mutual-exclusion abort all live in tpv102/104/205 + bp5_full
      drivers.
- [ ] Each driver's `--help` output (the comment-block at the top of
      the parsing section) lists the new flags.
- [ ] On `MFEM_USE_HDF5=NO` build, `--paraview-volume-hdf5` aborts
      at parse time with the documented message.
- [ ] Banner-output regression: existing `make test-tpv102-banner`
      (if such a test exists; otherwise add a new shell test) shows
      `Volume PV: ON (HDF5; ...)` on a default HDF5 build.

#### Dependencies

- Depends on: Phase 6.1 (constructor with mode parameter).
- Required by: Phase 6.3a, Phase 6.4 (`pv_bulk_out` parallel wiring).

---

### Phase 6.3a: Volume compression CLI for the PRIMARY collection (R-301)

#### Goal

After this sub-phase, users can compress the volume PV output of the
PRIMARY `pv_out` collection (the one written to `<output>/ParaView/`)
via the new `--paraview-volume-{zfp-tol,deflate-level}` flags.  This
closes the R-301 design gap: without these flags, BP5 users (no
`pv_bulk_out`) and tpv* users with `--paraview-bulk-dt` unset would
have NO way to compress volume PV from the CLI after Phase 6.4
re-routes `--paraview-bulk-*` to the secondary collection.

#### Files to Create

None.

#### Files to Modify

- `miniapps/seas/drivers/tpv102_driver.cpp`,
  `miniapps/seas/drivers/tpv104_driver.cpp`,
  `miniapps/seas/drivers/tpv205_driver.cpp`,
  `miniapps/seas/tests/verification/bp5_verification_full.cpp`:
  - Add CLI flags parsed alongside the Phase 6.3 volume-mode block:
    - `--paraview-volume-zfp-tol X` (real, default 0 = disabled)
    - `--paraview-volume-deflate-level N` (int 0..9, default -1 = MFEM default)
  - Parse-time validation (mirroring `--paraview-fault-*`):
    - `MFEM_USE_HDF5=NO` rejects `--paraview-volume-deflate-level` AND
      `--paraview-volume-zfp-tol` with the documented message
      (`"requires the seas-mfem build to define MFEM_USE_HDF5=YES"`).
    - `MFEM_USE_H5Z_ZFP=NO` rejects `--paraview-volume-zfp-tol` with
      `"requires MFEM_USE_H5Z_ZFP=YES"`.
    - Mutually exclusive: `--paraview-volume-zfp-tol > 0` AND
      `--paraview-volume-deflate-level >= 0` is rejected with
      `"choose either ZFP-accuracy OR deflate, not both"`.
  - Wire through to the PRIMARY `pv_out`:
    ```cpp
    #ifdef MFEM_USE_HDF5
    if (paraview_volume_zfp_tol > 0.0)
    {
       pv_out->SetVolumeHDFCompression(
          mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
          paraview_volume_zfp_tol);
    }
    else if (paraview_volume_deflate >= 0)
    {
       pv_out->SetVolumeHDFCompression(
          mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
          static_cast<double>(paraview_volume_deflate));
    }
    #endif
    ```

#### Detailed Requirements

1. **Naming separation** — `--paraview-volume-*` ALWAYS controls
   `pv_out` (the primary collection's volume PV).
   `--paraview-bulk-*` ALWAYS controls `pv_bulk_out` (the secondary
   wavefield collection in tpv* drivers; see §Phase 6.4).  Document
   in `CLAUDE.md` "ZFP lossy output" section per R-310.

2. **BP5 driver** has no secondary collection.  `--paraview-bulk-*`
   on `seas_bp5_full` therefore emits a one-time "no effect" warning
   per §Phase 6.4 step 3.  `--paraview-volume-*` on `seas_bp5_full`
   is the BP5-recommended way to compress the displacement /
   velocity volume PV.

3. **Default values** — `--paraview-volume-zfp-tol 0.0` (= disabled,
   default).  `--paraview-volume-deflate-level -1` (= use MFEM
   default, currently 6).

4. **Run-banner output** — when either flag is set, the run banner
   shows the active filter:
   ```
   Volume PV compression: ZFP @ 1e-3 (--paraview-volume-zfp-tol)
   ```
   or:
   ```
   Volume PV compression: deflate level 9 (--paraview-volume-deflate-level)
   ```

#### Interfaces

(Driver-only; no new C++ public API.  `SetVolumeHDFCompression` already
exists post-Phase 6.2.)

```bash
# tpv102 / tpv104 / tpv205 / seas_bp5_full now accept:
--paraview-volume-zfp-tol X            # real; ZFP accuracy on pv_out
--paraview-volume-deflate-level N      # int 0..9; deflate level on pv_out
```

#### Edge Cases to Handle

- **Active volume mode is VTU** (user passed `--paraview-volume-vtu`
  on a build that supports HDF, OR the build is non-HDF5):
  `SetVolumeHDFCompression` falls through to the existing warn-and-
  ignore branch from Phase 6.2.  No abort.
- **Both `--paraview-bulk-*` and `--paraview-volume-*` set on a
  tpv* driver with `--paraview-bulk-dt > 0`**: independent
  collections; both apply.  Document in `--help` that the two flag
  sets target different collections.
- **`--paraview-volume-zfp-tol` requested with VTU active**:
  per Phase 6.2 the warn-and-ignore message is emitted.  ZFP is
  still rejected at parse time on non-H5Z-ZFP builds (defense in
  depth).

#### Acceptance Criteria

- [ ] BP5 driver shell test: `seas_bp5_full --paraview
      --paraview-volume-zfp-tol 1e-3 --tfinal 1` produces a
      `<output>/ParaView/volume.vtkhdf` whose `displacement` dataset
      carries filter id 32013 (verified via `h5dump -p`).
- [ ] tpv102 driver shell test with both volume and bulk ZFP at
      different tolerances: `--paraview-volume-zfp-tol 1e-3
      --paraview-bulk-dt 0.05 --paraview-bulk-zfp-tol 1e-6` produces
      `volume.vtkhdf` with filter 32013 at tol=1e-3 AND
      `wave_bulk.vtkhdf` with filter 32013 at tol=1e-6.
      Compression ratios DIFFER between the two files.
- [ ] On `MFEM_USE_HDF5=NO` build: `--paraview-volume-zfp-tol 1e-3`
      AND `--paraview-volume-deflate-level 9` BOTH abort at parse
      time with the documented message.

#### Dependencies

- Depends on: Phase 6.2 (`SetVolumeHDFCompression` real-apply on HDF
  pv_dc_), Phase 6.3 (volume-mode CLI in place).
- Required by: Phase 6.4 (depends on Phase 6.3a so the two flag sets
  do not collide).

---

### Phase 6.4: pv_bulk_out independent ZFP (sub-phase 6d)

#### Files to Create

- `miniapps/seas/tests/unit/test_pv_bulk_out_zfp.cpp` — verifies that
  the secondary `pv_bulk_out` instance in tpv* drivers accepts an
  independent ZFP tolerance and the resulting VTKHDF carries filter
  32013 on velocity / sigma_yy / sigma_xy / sigma_xz.

#### Files to Modify

- `tpv102_driver.cpp`, `tpv104_driver.cpp`, `tpv205_driver.cpp`:
  - The existing `pv_bulk_out` construction (currently around line
    1996 in tpv102) gains the same `volume_mode` and ZFP wiring as
    `pv_out`:
    ```cpp
    // R-303: append-at-end constructor.  The secondary collection
    // uses a DISTINCT name "wave_bulk" because it lives in a
    // sibling directory ParaView_bulk/ and stores a different
    // dataset (wavefield only).  The primary collection uses the
    // cross-driver default "volume" per R-305.
    pv_bulk_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
       output_dir + "/ParaView_bulk", pmesh, order,
       /*collection_name=*/"wave_bulk", volume_mode);
    if (paraview_bulk_zfp_tol > 0.0)
    {
       pv_bulk_out->SetVolumeHDFCompression(
          ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
          paraview_bulk_zfp_tol);
    }
    else if (paraview_bulk_deflate >= 0)
    {
       pv_bulk_out->SetVolumeHDFCompression(
          ParaViewHDFDataCollection::HDFCompression::Deflate,
          static_cast<double>(paraview_bulk_deflate));
    }
    ```
  - Note: the existing `--paraview-bulk-zfp-tol` / `--paraview-bulk-deflate-level`
    flags (added in Phase 2d.3) currently route to `pv_out`'s
    `SetVolumeHDFCompression` — the warn-and-ignore stub.  After
    Phase 6.4 they route to `pv_bulk_out`'s
    `SetVolumeHDFCompression` (a real apply).  This is a semantic
    flip: `--paraview-bulk-*` historically affected the SECONDARY
    collection in dynamic-rupture studies (the wavefield-only
    one).  Phase 6.4 makes the flag name match its old intent.
  - Driver banner adds:
    ```
    Bulk collection: ON (ParaView_bulk/wave_bulk.vtkhdf, every 0.05 s, ZFP @ 1e-3)
    ```

#### Detailed Requirements

1. **`--paraview-bulk-*` semantics** are RE-routed to `pv_bulk_out`
   from `pv_out`.  R-310: this is a **SEMANTIC FLIP** from the
   Phase 2d.3 implementation, which routed `--paraview-bulk-*` to
   `pv_out->SetVolumeHDFCompression` (a warn-and-ignore stub
   because `pv_out`'s writer was VTU).  After Phase 6.4 the flag
   routes to `pv_bulk_out->SetVolumeHDFCompression` (a real apply,
   on the secondary wavefield collection that is only constructed
   when `--paraview-bulk-dt > 0` is also set).
   - **Migration impact**: users with Phase 2d.3 sbatch scripts
     that pass `--paraview-bulk-zfp-tol 1e-3` saw the warn-and-
     ignore message before; after Phase 6.4 the flag silently
     STARTS APPLYING — but to a different dataset (the secondary
     wavefield, not the primary volume) than the user might have
     expected.  The companion Phase 6.3a flag
     `--paraview-volume-zfp-tol` is the new way to compress the
     PRIMARY collection's volume PV.
   - `CLAUDE.md` "ZFP lossy output" section MUST be updated by
     §Phase 6.6 to document this flip and the corresponding
     `--paraview-volume-*` vs `--paraview-bulk-*` scope mapping.
     This is a hard documentation gate; an implementation report
     that lands Phase 6.4 without the CLAUDE.md update fails the
     §Phase 6 acceptance.
2. **`--paraview-volume-*` semantics** apply to BOTH `pv_out` and
   `pv_bulk_out` simultaneously (the binary-back-end choice is global
   to the run).  Allowing different modes per collection is
   unjustifiably complex.
3. **`pv_bulk_out` is only constructed when `--paraview-bulk-dt X`**
   (existing flag) is positive; if not, the bulk-side ZFP flags are
   silently irrelevant.  Add a one-time rank-0 warning if the user
   sets `--paraview-bulk-zfp-tol > 0` without `--paraview-bulk-dt`:
   ```
   warning: --paraview-bulk-zfp-tol set but --paraview-bulk-dt not
   provided; the secondary bulk collection is disabled, the flag has
   no effect.
   ```

#### Interfaces

(Driver-only; no new C++ public API.)

#### Edge Cases to Handle

- **BP5 driver does not construct a `pv_bulk_out`** (the secondary
  collection is a TPV-only construct).  `--paraview-bulk-*` flags on
  `seas_bp5_full` therefore emit the "no effect" warning and proceed.
  Document in the BP5 driver's `--help`.
- **TPV102 / 104 / 205 with `--paraview-bulk-dt 0.05` AND
  `--paraview-volume-vtu`**: both `pv_out` and `pv_bulk_out` are VTU
  writers; bulk ZFP flags emit the existing volume-side warn-and-ignore.

#### Acceptance Criteria

- [ ] `make test-pv-bulk-out-zfp` passes:
      - synthetic velocity in `pv_bulk_out` written with `tol=1e-3`;
      - filter id 32013 on the velocity dataset of the bulk
        `.vtkhdf`;
      - max-abs-err ≤ 1e-3 vs analytic field;
      - per-dataset compression ratio ≥ 4× vs deflate.
- [ ] tpv102 driver's banner shows the bulk-side ZFP setting when
      `--paraview-bulk-zfp-tol 1e-3` is passed.
- [ ] BP5 driver run with `--paraview-bulk-zfp-tol 1e-3` (without
      `--paraview-bulk-dt`) emits the "no effect" warning and exits 0.

#### Dependencies

- Depends on: Phase 6.3 (volume-mode CLI wiring); Phase 6.2 (real
  SetVolumeHDFCompression).
- Required by: Phase 6.6 (Frontera).

---

### Phase 6.5: Existing-test migration (sub-phase 6e)

#### Files to Create

None.

#### Files to Modify

- `miniapps/seas/tests/unit/test_io.cpp` — currently constructs a
  `ParaViewOutput` with the 3-arg constructor; verify the back-compat
  shim works.  No changes expected to the test source if the shim
  is correct; this is a regression-gate review.
- Any unit test that opens `<prefix>/ParaView/<basename>_*.vtu` files
  directly (rather than through `ParaViewOutput` API) must:
  1. Either pass `--paraview-volume-vtu` to its driver (preferred for
     existing tests that depend on the per-rank-VTU layout).
  2. OR migrate to read VTKHDF via h5py/pyvista (preferred for new
     tests).

#### Detailed Requirements

1. **Inventory pass** (R-306).  Run the following grep to enumerate
   every test file that hard-codes a legacy per-rank VTU / PVTU
   filename pattern.  The regex covers MFEM's actual naming
   conventions (basename + cycle / rank suffixes) and BOTH
   `ParaView/` and `ParaView_bulk/` directories:

   ```bash
   grep -rEn '(ParaView|ParaView_bulk)/[A-Za-z_0-9]+\.(vtu|pvtu)' \
       miniapps/seas/tests/
   ```

   For every match `(file, lineno, name)`, classify and resolve
   per the rule below.  Record each resolution `(a / b / c)` in the
   implementation report.

   **Decision rule:**
   - **(a) Opt into VTU**: the test PARSES the VTU file's content
     (via `ParseFaultVTUCellData`, regex on `<DataArray>`, or any
     XML walker) AND the parser is tied to a specific VTU layout.
     Add `pv->SetFaultOutputMode(...::Vtu)` (for fault-side
     parsers) or driver flag `--paraview-volume-vtu` (for
     volume-side parsers) so the test's parser keeps working.
   - **(b) File-existence check only**: the test only checks that a
     file exists / counts files.  MIGRATE to check
     `<output>/ParaView/volume.vtkhdf` (or
     `<output>/ParaView_bulk/wave_bulk.vtkhdf`) instead.
   - **(c) Time-series PVD check**: the test depends on the
     ParaView PVD time-series output.  MIGRATE to read the VTKHDF
     `Steps/` group via h5py.

2. **No silent breakage** — any test that fails after Phase 6.1 must
   either be fixed (preferred) or its breakage explicitly documented
   in the implementation report with a follow-up issue.

#### Acceptance Criteria

- [ ] `make test` (all unit tests) passes on a default `MFEM_USE_HDF5=YES`
      build.
- [ ] `make test` passes on `MFEM_USE_HDF5=NO` build (defaults to VTU
      automatically; nothing to migrate).

#### Dependencies

- Depends on: Phase 6.1.

---

### Phase 6.6: Frontera replication (sub-phase 6f)

#### Files to Create

None.

#### Files to Modify

- `miniapps/seas/jobs/bp5/bp5_1000m.sbatch` — already has the
  `module load phdf5 || hdf5/1.14.4 || true` line from Phase 2c.
  Add a banner-check comment showing the new flags:
  ```bash
  # Phase 6 — single-file VTKHDF for volume PV (default on
  # MFEM_USE_HDF5=YES builds).  To revert to legacy per-rank VTU:
  #   --paraview-volume-vtu
  #
  # To enable bulk-side ZFP (after Phase 2d.4 / Phase 6 are verified
  # on this cluster):
  #   --paraview-volume-zfp-tol 1e-3
  ```
- Production sbatch templates under `jobs/tpv102/`, `jobs/tpv104/`,
  `jobs/tpv205/` — same treatment.

#### Detailed Requirements

(Operational; no code changes beyond the sbatch comment block.)

1. Re-run a small TPV102 (1 s, np=400) on Frontera with the default
   HDF5 build.  Verify the output directory contains a single
   `wave.vtkhdf` instead of N+1 per-rank VTU files.
2. Re-run with `--paraview-volume-zfp-tol 1e-3` (assuming Phase 2d.4
   has built H5Z-ZFP under `$WORK`).  Verify the output is < 30 % the
   size of the deflate-only run.

#### Acceptance Criteria

- [ ] Frontera TPV102 microrun produces 1 file in `output_dir/ParaView/`
      and 1 file in `output_dir/ParaView_bulk/` (down from
      `(n_cycles * (nranks + 1)) + 1` per directory).
- [ ] BP5 1-yr microrun on Frontera with `--paraview` produces 1
      `volume.vtkhdf` (instead of the previous per-rank-VTU stream).
- [ ] `du -sh output_dir/` shows the targeted >5× reduction at default
      deflate; >20× at `--paraview-volume-zfp-tol 1e-3`.

#### Dependencies

- Depends on: Phase 6.1–6.5 complete and green locally.
- Required by: nothing.

---

## Phase 7: Output Size Estimator

### Goal

After this phase, a Python script
`miniapps/seas/scripts/estimate_output_size.py` reads the planned
run configuration (mesh path, polynomial order, driver name, schedule
parameters, compression settings) and prints an estimated total
output size (split into "fault" + "volume" + "stations" bytes), a
"fits-in-$SCRATCH" indicator vs a user-supplied quota, and a
per-field breakdown.  No MFEM dependency — runs in `python3
estimate_output_size.py ...` on a laptop in < 5 s.

### Sub-phase ordering

| Sub-phase | What | Where | Acceptance gate |
|---|---|---|---|
| **7a** | Mesh reader: count elements, count surface elements, count fault triangles given mesh + boundary tag list. | `_io_size_mesh.py` | `pytest test_estimate_output_size.py::test_mesh_reader_bp5_1000m` passes. |
| **7b** | Driver schemas: per-driver field list with shapes. | `_io_size_schemas.py` | `pytest ::test_schemas_match_driver_register_calls` passes. |
| **7c** | Schedule integrator: integrate `AdaptiveSchedule` over `tfinal` to estimate total cycles. | `_io_size_schedule.py` | unit test on synthetic V_max sequences (mirrors Phase 3 test). |
| **7d** | Compression-ratio table: empirical bytes-per-element for each (field, filter) combination. | `_io_size_compression.py` | `pytest ::test_compression_table_matches_known_runs` passes against three reference outputs. |
| **7e** | CLI + reporting: argparse, table-output, fits-in-quota check. | `estimate_output_size.py` | end-to-end smoke against three driver configurations. |
| **7f** | Sanity check vs real output: compare estimate against actual `du -sb`. | (manual / Frontera) | < 30 % overestimate, < 50 % underestimate. |

---

### Phase 7.1: Mesh reader (sub-phase 7a)

#### Files to Create

- `miniapps/seas/scripts/_io_size_mesh.py` — minimal Gmsh `.msh` v2 /
  v4 reader.  Returns:
  ```python
  @dataclass
  class MeshSummary:
      n_elements: int        # total volume elements (tets/hexes)
      n_vertices: int
      n_fault_faces: int     # surface elements with the fault tag
      n_boundary_faces: int  # surface elements with non-fault tags
      element_type: str      # "tet" / "hex" / ...
      element_order: int     # geometric order (1 for linear)
      fault_tag: int         # which Physical Surface ID is the fault
  ```

#### Files to Modify

None.

#### Detailed Requirements

1. **Format support** — Gmsh ASCII v2 and v4 (the formats currently
   produced by `gmsh -3 bp5/mesh/bp5.geo`).  Reject `.med`, `.cgns`
   with a clear error.
2. **Lazy parsing** — only read the `$Nodes` and `$Elements` blocks;
   skip everything else.  A 1.5 GB BP5 1000m mesh must parse in < 3 s
   on a laptop.
3. **Tag inference** — caller passes the fault tag (default 3 for
   BP5 per `miniapps/seas/CLAUDE.md` "Tandem mesh tags").  Override via
   the `--fault-tag N` CLI flag in §Phase 7.5.
4. **Error reporting** — every parse error includes the file path,
   line number, and the problematic token.

#### Interfaces

```python
def read_gmsh(path: pathlib.Path,
              fault_tag: int = 3) -> MeshSummary: ...

def read_inline_mesh(nx: int, ny: int, nz: int,
                     element_type: str = "tet") -> MeshSummary: ...
   # For drivers that use `--inline-mesh` (BP5 has this).  No file IO.
```

#### Edge Cases to Handle

- **Mesh file does not exist**: raise `FileNotFoundError` with the
  expanded path.
- **Mesh has no Physical Surface tags**: fall back to "all surface
  elements" → `n_fault_faces = n_boundary_faces`, with a stderr
  warning.
- **Higher-order curved-mesh elements** (Gmsh order > 1): count them
  as 1 element each but multiply n_dofs (in §Phase 7.2) by the order's
  DOF inflation factor.

#### Acceptance Criteria

- [ ] `read_gmsh("bp5/mesh/bp5_1000m.msh", fault_tag=3)` returns
      `n_elements ≈ 80,000`, `n_fault_faces ≈ 5,000`,
      `element_type == "tet"`.
- [ ] Parse time on the 1000m BP5 mesh < 3 s on a laptop.
- [ ] `read_inline_mesh(2, 2, 2, "tet")` returns `n_elements == 48`
      (8 cubes × 6 tets/cube via Gmsh's default tet refinement),
      `n_fault_faces == 0`.

#### Dependencies

- Depends on: nothing (pure Python).
- Required by: Phase 7.2.

---

### Phase 7.2: Driver field schemas (sub-phase 7b)

#### Files to Create

- `miniapps/seas/scripts/_io_size_schemas.py` — per-driver registered
  field tables:
  ```python
  @dataclass
  class FieldSchema:
      name: str
      dof_basis: str    # "L2_p0" / "H1_pK" / "L2_pK"
      n_components: int
      sizeof_dtype: int = 8  # double
      is_fault: bool = False
      is_volume: bool = True
      is_static: bool = False  # only emitted at first save
  ```
  Per-driver schemas:
  ```python
  TPV102_FIELDS = [
      FieldSchema("velocity", dof_basis="L2_pK_byNODES",
                  n_components=3),
      FieldSchema("mpi_rank", dof_basis="L2_p0", n_components=1,
                  is_static=True),
      FieldSchema("sigma_yy", ...), FieldSchema("sigma_xy", ...),
      FieldSchema("sigma_xz", ...),
      # 17 fault fields (12 standard + 5 _k4):
      FieldSchema("slip_dip", dof_basis="L2_p0_face",
                  n_components=1, is_fault=True, is_volume=False),
      ...
  ]
  TPV104_FIELDS = [...]
  TPV205_FIELDS = [...]
  BP5_FIELDS = [...]
  ```

#### Files to Modify

None (data file only).

#### Detailed Requirements

1. **Schema-source pinning** (R-308).  The schemas must match the
   actual `RegisterDomainField(...)` and `RegisterField(...)` calls
   in each driver / `paraview_output.hpp`.  Drift detection is
   implemented via the following regex:

   ```python
   _REGISTER_RE = re.compile(
       r'(?:RegisterDomainField|RegisterField)'
       r'\s*\(\s*'
       r'"([A-Za-z_][A-Za-z_0-9]*)"'  # field-name string literal
       r'\s*,'
   )
   ```

   The drift test scans these files (relative to project root):
   - `miniapps/seas/io/paraview_output.hpp`
   - `miniapps/seas/drivers/{tpv102,tpv104,tpv205}_driver.cpp`
   - `miniapps/seas/tests/verification/bp5_verification_full.cpp`

   For every match `(file, lineno, name)`, the test asserts that
   `name` appears in `get_driver_schema(driver)` for the driver
   inferred from the file path.  On failure, the test prints:

   ```
   Schema drift detected:
     Field "<name>" registered at <file>:<lineno>
     but not in <DRIVER>_FIELDS (in _io_size_schemas.py).
   To fix: add `FieldSchema("<name>", dof_basis=..., n_components=...)`
   to <DRIVER>_FIELDS, OR remove the registration if the field
   is no longer needed.
   ```

   The regex deliberately requires a STRING-LITERAL field name; it
   does NOT match indirect calls such as
   `RegisterField(get_name(...), ...)`.  Such a future driver pattern
   requires a manual schema entry; document the limitation in
   `_io_size_schemas.py`'s module docstring.
2. **DOF count formula**:
   - `L2_p0` element-wise: `n_dofs = n_elements`
   - `L2_p0_face` (fault): `n_dofs = n_fault_faces`
   - `L2_pK_byNODES` for tet of polynomial order K:
     `n_dofs = n_elements × (K+1)(K+2)(K+3)/6`
   - `H1_pK` for tet of polynomial order K:
     `n_dofs = n_vertices` for K=1; for K≥2 use the standard
     simplex DOF count `binomial(K+3, 3)` per element minus shared
     DOFs (estimator approximates by `n_elements × binomial(K+3, 3)
     / shared_factor`, where `shared_factor ≈ 4` for tets at K≤4 —
     pinned via the unit test against MFEM-actual-DOF counts on a
     small test mesh).
3. **Static-field handling** — fields with `is_static=True`
   (e.g. `mpi_rank`, `param_a`, `param_Dc`) are emitted once and
   reused across cycles (per Phase 2b §2 of the prior plan).
   Estimator counts them once, not per-cycle.

#### Interfaces

```python
def get_driver_schema(driver_name: str) -> List[FieldSchema]: ...
   # Returns the registered-field list for one of:
   # "tpv102", "tpv104", "tpv205", "bp5".

def estimate_dofs(field: FieldSchema,
                  mesh: MeshSummary,
                  order: int) -> int: ...
   # Returns the DOF count for `field` on `mesh` at polynomial order.
```

#### Edge Cases to Handle

- **Driver name not in the schema table**: raise `ValueError` listing
  the supported drivers.
- **`-pv-low-order` flag** (tpv*): forces `H1_p1` for displacement and
  `L2_p0` for everything else.  Estimator must accept an `order_override`
  parameter that defaults to the driver's natural order.

#### Acceptance Criteria

- [ ] `pytest test_estimate_output_size.py::test_schemas_match_driver_register_calls`
      passes: every `pv_out->RegisterDomainField("name", ...)` and
      `pv_.RegisterField("name", ...)` in the four drivers appears in
      the corresponding schema list.
- [ ] `estimate_dofs(velocity_field, bp5_1000m_mesh, order=1)` is
      within 5 % of the actual MFEM-reported DOF count for that
      mesh / order.

#### Dependencies

- Depends on: Phase 7.1.
- Required by: Phase 7.4.

---

### Phase 7.3: Schedule integrator (sub-phase 7c)

#### Files to Create

- `miniapps/seas/scripts/_io_size_schedule.py` — pure-Python
  `AdaptiveSchedule` mirror that integrates over `tfinal`:
  ```python
  @dataclass
  class ScheduleConfig:
      tfinal: float
      v_coseismic: float = 1e-3
      v_nucleation: float = 1e-7
      dt_coseismic: float = 60.0          # seconds
      dt_nucleation: float = 1.0
      dt_interseismic: float = 3.156e7    # 1 yr
      hysteresis_factor: float = 1.0
      max_total_snapshots: int = 0
      output_every_n_steps: int = 0       # step-based override
      fixed_dt: float = 0.0                # fixed-dt override
      n_events: int = 0                    # estimated coseismic events
      avg_event_duration_s: float = 30.0   # how long each event lasts

  def estimate_n_writes(cfg: ScheduleConfig) -> int: ...
  ```

#### Files to Modify

None.

#### Detailed Requirements

1. **Three regime time-budget**:
   - coseismic budget = `n_events × avg_event_duration_s`
   - nucleation budget = `n_events × NUCLEATION_DURATION_S[driver]`,
     where `NUCLEATION_DURATION_S` is a per-driver constant in
     `_io_size_compression.py` (R-307).  Initial values, calibrated
     against post-Phase-2 Frontera reference runs:
     - `bp5`: `86_400.0`  (~1 day; Dieterich-Ruina nucleation phase
       takes hours-to-days for the V threshold to climb from
       `v_nucleation = 1e-7 m/s` to `v_coseismic = 1e-3 m/s`
       per `miniapps/seas/CLAUDE.md` "What Constitutes a Regression")
     - `tpv102` / `tpv104` / `tpv205`: `1.0`  (dynamic-rupture
       nucleation completes in a single time-step window)
     The 100 s figure used in the prior plan draft was wrong by
     ~3 orders of magnitude for BP5 (R-307); the Phase 7.6 reference
     run #1 (BP5 250 yr / 5 events) must predict ~432,000 nucleation
     writes when uncapped (5 × 86,400 s ÷ 1 s `dt_nucleation`).
   - interseismic budget = `tfinal − coseismic − nucleation`
   - n_writes_co     = `coseismic_budget / dt_coseismic`
   - n_writes_nu     = `nucleation_budget / dt_nucleation`
   - n_writes_inter  = `interseismic_budget / dt_interseismic`
   - total = sum of three
2. **Cap interaction** — if `max_total_snapshots > 0` AND total > cap,
   total is clamped to `cap + n_events` (per Phase 3 acceptance:
   "K + n_events").
   - Note: BP5 production runs typically set
     `--paraview-max-snapshots 5000`, which clamps the total well
     below the uncapped ~430k.  The nucleation-budget calibration
     therefore matters most for **uncapped** runs and for accuracy
     when the cap is loose (`max_total_snapshots > 100k`).
3. **Step-based / fixed-dt overrides** take precedence per Phase 3
   §Edge Cases.
4. **Events default**:
   - For BP5: `n_events = round(tfinal / 240yr)` (recurrence default
     per `miniapps/seas/CLAUDE.md` "What Constitutes a Regression").
   - For TPV: `n_events = 1` (one nucleation event per dynamic-rupture
     run).

#### Interfaces

```python
@dataclass
class ScheduleConfig: ...

def estimate_n_writes(cfg: ScheduleConfig) -> int:
    """Return total number of fault-side write events the schedule
    will trigger over `cfg.tfinal`.  Caller multiplies this by
    per-write bytes (from §7.4) to get total fault-output bytes.
    Volume-side n_writes is the same unless `volume_pv_dt > 0`,
    in which case use `estimate_n_volume_writes(cfg)`."""

def estimate_n_volume_writes(cfg: ScheduleConfig,
                             volume_pv_dt: float) -> int: ...
```

#### Edge Cases to Handle

- **`max_total_snapshots = 1`**: returns 1 + n_events (Phase 3
  acceptance edge case).
- **`tfinal < dt_interseismic`**: returns at least 1 (the t=0
  initial snapshot all drivers emit unconditionally).
- **`fixed_dt > 0`**: returns `tfinal / fixed_dt`, ignoring the
  three-regime breakdown.

#### Acceptance Criteria

- [ ] `pytest ::test_schedule_integrator_no_events_returns_tfinal_over_dt_interseismic`
      passes.
- [ ] `pytest ::test_schedule_integrator_5_events_matches_phase3_synthetic`
      passes — same V_max sequence as `test_paraview_schedule_cap.cpp`
      yields ≈ same write count (± 1 per event).
- [ ] `pytest ::test_schedule_integrator_cap_clamps`
      passes — cap=K returns ≤ K + n_events.

#### Dependencies

- Depends on: nothing.
- Required by: Phase 7.5.

---

### Phase 7.4: Compression-ratio table (sub-phase 7d)

#### Files to Create

- `miniapps/seas/scripts/_io_size_compression.py` — table mapping
  `(field_kind, filter)` to bytes-per-DOF.  R-309: every numerical
  entry whose source is the literal string `"PLACEHOLDER"` is an
  initial guess pending §Phase 7.6 calibration against on-disk
  bytes; the runtime warns when a placeholder is consulted (see
  Detailed Requirements §1).  Anchors are taken from the prior
  Phase 1 plan's quantitative-targets table (deflate@VTU vs raw
  was ~3-5×; ZFP@1e-3 was ~10-20× on smooth-data fields).
  ```python
  COMPRESSION_RATIOS = {
      # (field_kind, filter) -> (bytes_per_dof, source_or_PLACEHOLDER)
      ("any", "none"):                 (8.0, "raw IEEE-754 double"),
      ("any", "vtu_ascii"):            (25.0, "legacy ASCII baseline"),
      ("any", "vtu_binary"):           (8.0, "raw double"),
      ("smooth_fp", "deflate_6"):      (3.2, "PLACEHOLDER"),
      ("smooth_fp", "zfp_1e-3"):       (0.4, "PLACEHOLDER"),
      ("smooth_fp", "zfp_1e-12"):      (4.8, "PLACEHOLDER"),
      ("rough_fp", "deflate_6"):       (5.5, "PLACEHOLDER"),
      ("rough_fp", "zfp_1e-12"):       (5.0, "PLACEHOLDER"),
      ("integer_index", "deflate_6"):  (0.5, "PLACEHOLDER"),
      ...
  }
  ```
  Plus a classifier `classify_field(name, driver) -> str` that maps
  a field name to its compression "kind" (`smooth_fp` /
  `rough_fp` / `integer_index` / `static_2d`).

#### Files to Modify

None.

#### Detailed Requirements

1. **Calibration source + runtime warning** (R-309) — every entry's
   second tuple element is either a concrete Frontera run name
   ("calibrated") or the literal string `"PLACEHOLDER"` (initial
   guess).  `bytes_per_dof(field, filter, driver)` emits a
   one-time `warnings.warn(...)` per `(field_kind, filter)` pair
   whose source equals `"PLACEHOLDER"`:
   ```
   warning: compression ratio for (field_kind=smooth_fp,
   filter=zfp_1e-3) is a PLACEHOLDER pending Phase 7.6 calibration;
   the size estimate may be off by >2× until calibration completes.
   ```
   The `--quiet` flag on the estimator suppresses these warnings;
   absence of `--quiet` lets the user see which entries are not
   yet calibrated.  When §Phase 7.6 runs the calibration, every
   `"PLACEHOLDER"` is replaced with a concrete run name; the
   §Phase 7.4 acceptance test asserts the table is within ±15 %
   of the on-disk bytes for any entry whose source is NOT
   `"PLACEHOLDER"`.
2. **Out-of-table fallback** — `(field_kind, filter)` not in the
   table → assume `8.0` bytes/DOF (raw double) and emit a stderr
   warning.

#### Interfaces

```python
def bytes_per_dof(field: FieldSchema,
                  filter_name: str,
                  driver: str) -> float: ...

def classify_field(field_name: str, driver: str) -> str:
   """Return one of: 'smooth_fp', 'rough_fp', 'integer_index',
   'static_2d', 'unknown'."""
```

#### Edge Cases to Handle

- **Filter name unrecognised**: raise `ValueError` listing supported
  filters: `none`, `vtu_ascii`, `vtu_binary`, `deflate_<N>` for
  N=0..9, `zfp_<tol>` for tol in {1e-3, 1e-6, 1e-9, 1e-12}.

#### Acceptance Criteria

- [ ] `pytest ::test_compression_table_matches_known_runs` passes
      against three downloaded reference HDF files (or skips if not
      present).
- [ ] `bytes_per_dof(velocity_field, "zfp_1e-3", "tpv102")` returns
      a value within `[0.3, 0.5]` (matches the table entry calibrated
      on TPV102 wave runs).

#### Dependencies

- Depends on: nothing.
- Required by: Phase 7.5.

---

### Phase 7.5: CLI and reporting (sub-phase 7e)

#### Files to Create

- `miniapps/seas/scripts/estimate_output_size.py` — main CLI entry
  point:
  ```bash
  python3 estimate_output_size.py \
      --driver bp5 \
      --mesh bp5/mesh/bp5_1000m.msh \
      --tfinal 250yr \
      --paraview \
      --paraview-fault-zfp-tol 1e-12 \
      --paraview-max-snapshots 5000 \
      --no-volume-pv \
      --np 800 \
      --scratch-quota 1TB
  ```
  Output:
  ```
  Mesh: bp5_1000m.msh — 80,432 tets, 5,120 fault triangles
  Driver: bp5 (P1 elasticity, 17 fault fields, 1 volume field)
  Schedule: adaptive cap=5000 over 250yr → 5,012 writes
            (12 events: 12 coseismic, 12 nucleation, 4988 interseismic)

  Per-write bytes:
    fault    (zfp_1e-12, 12 fields × 5,120 cells × 8 B × 5.0)
              = 2,457,600 B
    volume   (suppressed via --no-volume-pv)
    stations (BP5 7 stations × 8 B × 8 cols)        = 448 B

  Total run output:
    fault    : 12.3 GB
    volume   : 0
    stations : 14 MB
    -----------------
    TOTAL    : 12.3 GB

  Quota check ($SCRATCH = 1.0 TB): 1.2 % of quota — OK
  ```

#### Files to Modify

None.

#### Detailed Requirements

1. **Time-unit suffixes** (R-313) — `--tfinal` accepts the suffixes
   `s` (seconds, ALSO THE DEFAULT WHEN NO SUFFIX), `min`, `hr` /
   `h`, `day` / `d`, `wk` / `w`, `yr` / `y`.  Whitespace between the
   number and the suffix is permitted (`1 day` and `1day` both
   valid).  Decimal values are permitted (`1.5day`, `0.5yr`).
   Examples: `250yr`, `3.156e7` (= 3.156×10⁷ s, no suffix → seconds),
   `1day`, `36hr`, `2.5wk`.  Unknown suffixes raise `ValueError`
   listing the supported set.
2. **Driver-flag passthrough** — every `--paraview-*` / `--no-volume-pv` /
   `--volume-pv-dt` / `--paraview-bulk-*` / `--paraview-volume-*` flag
   accepted by the actual drivers must be parsed by the estimator AND
   route through to the schedule + compression configuration.  Drift
   detection: a unit test reads the four drivers' CLI parsers and
   asserts the estimator parses every relevant flag.
3. **Quota check** — `--scratch-quota` accepts `1TB`, `500GB`, etc.
   Output color codes (when stdout is a tty) use red for >100 %, yellow
   for 80–100 %, green otherwise.  R-312: byte-size formatting uses a
   hand-rolled `format_bytes(n)` helper (no third-party dependency)
   that returns `"12.3 GB"` / `"1.2 TB"` / `"456 MB"` with consistent
   3-significant-digit precision; column alignment uses Python's
   built-in f-string `>12s` width specifier rather than `tabulate`,
   so the script remains standard-library-only.  This handles the
   "future bigger run" concern — the byte-format helper auto-scales
   from B → KB → MB → GB → TB → PB so the column never overflows.
4. **JSON output** — `--json` flag prints the result as machine-readable
   JSON for scripting.
5. **`--explain` mode** — prints the formula chain that led to the
   estimate (per-field DOFs, per-cycle bytes, n_writes, total) so
   the user can audit any number that looks wrong.

#### Interfaces

```python
def estimate(driver: str,
             mesh: MeshSummary,
             order: int,
             schedule: ScheduleConfig,
             filters: Dict[str, str],   # {"fault": "zfp_1e-12", ...}
             scratch_quota_bytes: Optional[int] = None) -> EstimateReport:
    ...

@dataclass
class EstimateReport:
    fault_bytes: int
    volume_bytes: int
    station_bytes: int
    total_bytes: int
    per_write_bytes_fault: int
    per_write_bytes_volume: int
    n_writes_fault: int
    n_writes_volume: int
    quota_fraction: Optional[float]
```

#### Edge Cases to Handle

- **No mesh provided AND `--inline-mesh`**: use the inline-mesh
  reader (Phase 7.1).
- **`--paraview` not set**: fault and volume bytes both report 0;
  station bytes still computed.
- **Unknown driver**: `argparse` `choices` constraint catches before
  the estimator runs.

#### Acceptance Criteria

- [ ] `make test-estimate-output-size` (new Makefile target)
      runs the script against three pre-staged scenarios:
      1. BP5 250 yr / 800 ranks / cap=5000 / fault-only / ZFP 1e-12
      2. TPV102 1 s / 400 ranks / volume-on / deflate-6
      3. TPV205 1 s / 200 ranks / volume-on / ZFP 1e-3
      Each prints to stdout, exits 0, and the printed total is
      within the per-scenario tolerance defined in §Phase 7.6.
- [ ] `--json` output round-trips through `json.loads` cleanly.
- [ ] `--explain` output names every formula input.

#### Dependencies

- Depends on: Phase 7.1, 7.2, 7.3, 7.4.
- Required by: Phase 7.6 (sanity-check vs real run).

---

### Phase 7.6: Sanity check vs real output (sub-phase 7f)

#### Files to Create

- `miniapps/seas/scripts/test_estimate_output_size.py` — pytest
  test suite covering 7.1–7.5 (referenced throughout above).
- `miniapps/seas/scripts/REFERENCE_RUNS.md` — short table listing
  three reference runs and their measured output sizes.  Used as the
  ground truth for the ±30 % accuracy gate.

#### Files to Modify

None.

#### Detailed Requirements

1. **Reference run #1: BP5 250 yr / 800 ranks / cap=5000 / fault-only /
   deflate-6** — measured ~12 GB on Frontera (post-Phase 2 deploy).
   Estimate must be within `[8.4 GB, 15.6 GB]` (i.e. ±30 %).
2. **Reference run #2: TPV102 1 s / 400 ranks / volume-on / deflate-6** —
   measured ~3.2 GB on Frontera.  Estimate ∈ `[2.24 GB, 4.16 GB]`.
3. **Reference run #3: TPV205 1 s / 200 ranks / volume-on / ZFP 1e-3** —
   measured ~0.4 GB on Frontera.  Estimate ∈ `[0.28 GB, 0.52 GB]`.
4. **Underestimate-bias** — if the estimator UNDERESTIMATES any
   reference run by > 50 %, that's a hard failure (the user would
   provision insufficient $SCRATCH).  Overestimate by up to 100 %
   is acceptable (over-provision is fine).

#### Interfaces

(Pytest only.)

#### Edge Cases to Handle

- **Reference run not yet measured** (Frontera deploy pending): the
  test file marks the corresponding test cases `@pytest.mark.skip(reason=
  "reference measurement pending Phase 6.6 Frontera replication")`.

#### Acceptance Criteria

- [ ] `pytest test_estimate_output_size.py` exit code 0 with the
      reference runs available.
- [ ] All three reference scenarios estimate within ±30 % of measured.
- [ ] No reference scenario underestimates by > 50 %.

#### Dependencies

- Depends on: Phase 7.1–7.5; the three reference runs being available
  on disk.
- Required by: nothing.

---

## Testing Strategy

### Per-phase

- **Phase 6.1**: unit test verifying both VolumeOutputMode paths
  produce ParaView-readable output for a 1-element mesh + 1 GF.
- **Phase 6.2**: `test_volume_hdf_compression.cpp` — filter id +
  per-dataset accuracy + compression-ratio assertions, mirroring
  `test_vtkhdf_zfp.cpp` style.
- **Phase 6.3**: shell test verifying the parse-time abort message
  on a non-HDF5 build for `--paraview-volume-hdf5`.
- **Phase 6.4**: `test_pv_bulk_out_zfp.cpp` — independent ZFP on the
  secondary collection.
- **Phase 6.5**: re-run the existing `make test` matrix on both HDF5
  and non-HDF5 builds.
- **Phase 6.6**: Frontera microrun acceptance (operational).
- **Phase 7.1–7.4**: pytest unit tests; no MFEM dependency.
- **Phase 7.5**: end-to-end CLI smoke against three pre-staged
  scenarios.
- **Phase 7.6**: sanity-check vs measured Frontera output.

### Cross-phase

- **Phase 6 + Phase 7 integration test**: run the estimator with
  the same flags as a real BP5 microrun, then run the actual driver,
  compare estimate vs `du -sb output_dir/`.  Acceptance: ±30 %.
  This is the same gate as Phase 7.6 acceptance #1, just exercised
  end-to-end.

### Regression gates

- `make test test-bp5-smoke test-tpv102` continue to pass on default
  HDF5 builds.
- `make test` continues to pass on non-HDF5 builds (default falls
  back to VTU automatically).
- `make test-paraview-tolerance-isolation` (Phase 2d.3 deferred half)
  passes both halves now (bulk-side ZFP applies after Phase 6.2).

---

## Risk Assessment

### What could go wrong (Phase 6)

1. **Existing test files that hard-code the legacy
   `<output>/ParaView/<basename>_<rank>_<cycle>.vtu` filename**:
   these will break silently (no file at the expected path) on the
   new HDF5 default.  **Mitigation**: §Phase 6.5 inventory pass +
   per-test resolution.  All resolutions documented in the
   implementation report.

2. **Parallel HDF5 collective-write deadlock** on a rank that
   returns early from the `pv_dc_->Save()` path (e.g. zero-DOF
   partition).  MFEM's `ParaViewHDFDataCollection::Save` is
   collective on the mesh's communicator; an empty partition is
   well-defined per MFEM's tests.  **Mitigation**: §Phase 6.1
   acceptance test on `np=4` with one rank holding zero DOFs.

3. **`pv_bulk_out` constructor with `mesh = pmesh`**: TPV
   drivers construct a SECOND ParaViewOutput on the SAME mesh.
   `ParaViewHDFDataCollection` opens its own HDF5 file (different
   collection name) so there's no file-handle conflict — but the
   underlying VTKHDF object is per-collection, not per-mesh.
   Verified by §Phase 6.4 acceptance test.

4. **VTKHDF restart-mode** is required for any sbatch script that
   uses `requeue` (e.g. Frontera 48 h walls with checkpointing).
   Currently `UseRestartMode(true)` is not wired through Phase 6;
   document as a Phase 6 follow-up — running without restart-mode
   means a requeued job's later writes overwrite prior cycles.
   **Mitigation**: explicit warning in the run banner when
   `--restart` is set: "VTKHDF restart-mode not yet supported;
   prior cycles will be overwritten on requeue."

5. **Driver-side default-flip surprise**: long-running production
   sbatch scripts (BP5 250 yr) that do not explicitly set
   `--paraview-volume-vtu` switch to HDF5 silently.  This is the
   intended behaviour, but a user re-running an old script and
   expecting per-rank VTU will be surprised.  **Mitigation**: bold
   call-out in the §Phase 6.6 sbatch template comment block AND
   in `miniapps/seas/CLAUDE.md` "ParaView output mode" section.

### What could go wrong (Phase 7)

1. **Compression-ratio table drift over time** — as the seas
   numerics evolve (different fields, different solver tolerances),
   the smooth/rough field classification can shift.  **Mitigation**:
   §Phase 7.4 acceptance test re-calibrates against three measured
   reference runs; re-run when any of the three drivers ships a
   semantically new field.

2. **Schedule integrator over-counts when a real run hits the
   schedule cap mid-run**: my Phase 3 cap math assumes a uniform
   distribution of events through `tfinal`; a clustered-event run
   (e.g. all 12 events in the first half of a 500 yr BP5 run)
   triggers the cap-exhausted branch earlier than estimated.
   **Mitigation**: estimator prints the cap-exhausted-time estimate;
   user can spot-check vs their event history.

3. **Mesh reader fails on Gmsh 5.x format**: format is binary-mostly
   and includes new tags.  **Mitigation**: detect the version line
   and abort with a clear "supported: Gmsh ASCII v2 / v4" message;
   defer v5 support to a follow-up.

4. **User runs estimator on a laptop without the mesh file**: the
   tool crashes with `FileNotFoundError`.  **Mitigation**: doc the
   workflow as "scp the mesh to your laptop OR pass `--mesh-summary`
   (a JSON pre-computed from a prior `--dump-mesh-summary` run on
   Frontera)".  Add `--dump-mesh-summary` flag in §Phase 7.5 if
   user feedback requests it.

### Known tricky areas in the existing code

- `paraview_output.hpp:RegisterDomainField` (line 333) currently
  forwards to the concrete `pv_.RegisterField`.  After Phase 6.1 it
  forwards to `pv_dc_->RegisterField` (virtual via
  `ParaViewDataCollectionBase`).  The base class declares
  `RegisterField` as a virtual override of `DataCollection`'s
  non-virtual `RegisterField`, so existing call sites keep working
  but a `dynamic_cast`-based dispatch would be wrong (no
  `dynamic_cast` needed; virtual dispatch handles it).

- `bp5_verification_full.cpp:1618` registers
  `seas_op.GetDisplacement()` — a `ParGridFunction` reference cast
  to a non-const pointer.  After Phase 6.1 this still works because
  `RegisterField` accepts `GridFunction*`.

- `tpv102_driver.cpp:1996` constructs `pv_bulk_out` with the 3-arg
  back-compat shim.  Phase 6.4 changes this to the 5-arg form so
  the bulk side gets the volume-mode parameter.  Verify on rebuild.

---

## Open questions to resolve before implementation

1. **(RESOLVED 2026-05-09 / R-305)** **Default collection name** —
   `"volume"` for cross-driver uniformity with the fault path's
   `"fault_surface"`.  All four drivers (tpv102, tpv104, tpv205,
   bp5_full) use the default; none passes a custom collection name.
   The constructor parameter exists for future flexibility but is
   not exercised by any driver in this plan.

2. **VTKHDF restart-mode** for requeue-aware runs (Risk #4) —
   detect-and-abort vs detect-and-warn-and-overwrite.  Plan §Phase
   6.6 proposes warn-and-overwrite + a Phase 6 follow-up to
   actually wire `UseRestartMode(true)`.  Confirm before
   implementation; if the user wants abort-instead, change the
   §Phase 6.6 sbatch comment block accordingly.

3. **Phase 7 estimator mesh-summary JSON** — should the estimator
   support a pre-computed mesh-summary JSON so the user runs the
   estimator on a laptop without scp'ing the multi-GB mesh?
   Plan §Risk 4 of Phase 7 proposes adding `--dump-mesh-summary`
   if user feedback requests it; defer to user.

The implementation agent should resolve open questions 1–3 with
the user BEFORE starting Phase 6.1 / Phase 7.1.

---

## Acceptance for the plan as a whole

- [ ] Phases 6.1–6.5 complete locally; `make test` green on both
      HDF5 and non-HDF5 builds.
- [ ] Phase 6.6 Frontera microrun shows expected file-count + size
      reductions.
- [ ] Phase 7.1–7.5 complete locally; `pytest` green.
- [ ] Phase 7.6 sanity check passes against three reference runs
      within ±30 % accuracy.
- [ ] CLAUDE.md "ParaView output mode" and "ZFP lossy output"
      sections updated to describe the new defaults and the
      `--paraview-volume-*` flags.
- [ ] One Frontera production sbatch template (`bp5_1000m.sbatch`
      or similar) verified end-to-end with the new defaults.
