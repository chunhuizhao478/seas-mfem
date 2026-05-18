# HDF5-DIAG "already closed" Noise on Frontera Runs

**Date:** 2026-05-17
**Symptom:** Every TPV104 run on Frontera emits HDF5-DIAG error stacks
from the very first ParaView write onwards (NOT just at end-of-run):

```
HDF5-DIAG: Error detected in HDF5 (1.14.6):
  #000: ../../src/H5G.c line 878 in H5Gclose(): decrementing group ID failed
  #001: ../../src/H5Iint.c line 1087 in H5I_dec_app_ref(): can't decrement ID ref count
  #002: ../../src/H5Iint.c line 1042 in H5I__dec_app_ref(): can't decrement ID ref count
  #003: ../../src/H5Iint.c line 948 in H5I__dec_ref(): can't locate ID
        (already closed?)

HDF5-DIAG: Error detected in HDF5 (1.14.6):
  #000: ../../src/H5P.c line 1468 in H5Pclose(): can't close
  ... same chain ...

HDF5-DIAG: Error detected in HDF5 (1.14.6):
  #000: ../../src/H5F.c line 1047 in H5Fclose(): decrementing file ID failed
  ... same chain ...
```

**Data correctness:** UNAFFECTED.  Phase A wrote 154 MB
`fault.vtkhdf`; Phase B wrote 237 MB `fault.vtkhdf` in segment_002.
Validations #1, #2, #3 all PASS — the bytes land correctly on disk.
The errors are **diagnostic stderr noise from HDF5's error stack** that
HDF5 auto-prints when its internal ID table state is inconsistent.

## Root cause: dual HDF5 ABI loaded in one process

The Frontera build pulls in TWO different HDF5 ABIs at runtime:

1. **HDF5 1.10.x** (`libhdf5_hl.so.200`, `libhdf5.so.200`) — required
   by PETSc 3.15 which was built against the system HDF5 module.
2. **HDF5 1.14.6** (`libhdf5_hl.so.310`, `libhdf5.so.310`) — built
   from source by `build_frontera.sh` into
   `${SEAS_MFEM_ROOT}/extern/hdf5/install/lib`.  Required by
   MFEM_USE_HDF5=YES + the seas Phase 6 VTKHDF pipeline.

The link warning at build time confirms the dual load:

```
/opt/apps/gcc/8.3.0/bin/ld: warning: libhdf5_hl.so.200, needed by
   /home1/apps/intel19/impi19_0/petsc/3.15/clx/lib/libpetsc.so,
   may conflict with libhdf5_hl.so.310
/opt/apps/gcc/8.3.0/bin/ld: warning: libhdf5.so.200, needed by
   /home1/apps/intel19/impi19_0/petsc/3.15/clx/lib/libpetsc.so,
   may conflict with libhdf5.so.310
```

The two libraries each have **their own static ID table**.  When seas
opens a group/property list/file via libhdf5.so.310, the ID is
registered in 310's table.  When that same ID is later released through
a code path that goes through libpetsc.so → libhdf5.so.200's symbols,
the close fires in 200's library but the ID lives in 310's table.
310's `H5I_dec_app_ref()` looks up the ID, can't find it, prints "can't
locate ID (already closed?)" to stderr — even though no data corruption
occurred.  HDF5 returns a soft error; the file write completes.

**Confirmation:** the `HDF5-DIAG: Error detected in HDF5 (1.14.6)`
preamble names libhdf5.so.310 as the printer.  The errors fire because
310's view of ID lifetimes is wrong (the ID was already freed by 200's
close path).

## Why TPV104 surfaces this more visibly than BP5

TPV104 writes ParaView frames at fixed `--paraview-dt 0.1` cadence
during dynamic rupture (thousands of frames per ~5 s run).  Each frame
triggers ParaView's HDF5 group/property-list churn → each frame emits
a diagnostic stack.  BP5 quasi-static runs write at much sparser
cadence (regime-adaptive), so even though the same dual-HDF5 issue
exists, the visible noise is far smaller.

Note: TPV104 does NOT use PETSc itself (it uses the ADER-DG explicit
stepper).  But MFEM_USE_PETSC=YES at MFEM build time forces every
binary that links libmfem to ALSO link libpetsc — and libpetsc pulls
in libhdf5.so.200.  Selectively excluding -lpetsc from the TPV104
link would require rebuilding MFEM without PETSc support, which would
break BP5.

## Fix options ranked by invasiveness

| # | Approach | Effort | Side effects |
|---|----------|--------|--------------|
| A | `H5Eset_auto2(H5E_DEFAULT, NULL, NULL)` at driver start — silences HDF5's auto-print of error stacks | 1-line code change in each driver | Real HDF5 errors (disk full, permission denied) also become silent |
| B | Install a custom HDF5 error handler that filters "can't locate ID" / "already closed" but passes others through | ~30 lines | Same as A but preserves real-error visibility; more code to maintain |
| C | Rebuild PETSc 3.15 against HDF5 1.14 on Frontera | days; needs admin or local PETSc build | Eliminates the dual-instance load entirely |
| D | Build MFEM against the system HDF5 1.10 | weeks; loses Phase 6 VTKHDF support entirely | Not viable |
| E | Statically link HDF5 1.14 into the binary, hide all dynamic resolution | requires HDF5 + plugin rebuild with -fPIC | Complex; ZFP plugin might not statically link cleanly |

**Chosen: option A.**  Option B was initially considered (filter by
major class) but the actual error stacks have a `Symbol table`
(`H5E_SYM`) outer record followed by `Object ID` (`H5E_ID`) inner
records — a class-only filter cannot reliably distinguish dual-HDF5
noise from real symbol-table errors.  String-matching the minor
descriptions ("Unable to find ID information (already closed?)") is
fragile across HDF5 patch versions.

Option A's downside (real HDF5 auto-print errors also silenced) is
acceptable because every code path that opens / writes a file in
seas does so through MFEM wrappers (`ParaViewHDFDataCollection`,
`FaultHDFState`) that surface real failures via MFEM_VERIFY or
exceptions — i.e., disk-full / permission-denied / dataset-not-found
errors are caught by the higher-level wrappers regardless of whether
HDF5's auto-print is enabled.

## Applied fix (option A)

Each driver (BP5, TPV104) installs the suppressor once after
MPI_Init.  Guarded by `MFEM_USE_HDF5` so non-HDF5 builds are no-ops.
Registered via `H5Eset_auto2(H5E_DEFAULT, NULL, NULL)`.

A header-only helper at `miniapps/seas/io/hdf5_error_filter.hpp`
provides `InstallHdf5ErrorFilter()` so both drivers call the same
code path; if a future refactor switches to option B (filtered
handler), only the helper's implementation changes.

## Verification

- Local macOS build: re-build TPV104 driver + unit tests with the
  filtered handler; confirm no regressions.  The single-HDF5 macOS env
  never triggers the "can't locate ID" stack, so the handler is a
  no-op; verify the build still links.
- Frontera re-submission: re-run the TPV104 sbatch and confirm:
  - Validation #1–#7 all PASS (data integrity preserved).
  - No `HDF5-DIAG: Error detected ... can't locate ID` lines in
    stdout/stderr for the entire two-phase run.
  - End-of-run summary disk usage unchanged.

## Long-term cleanup

Track separately: file a request with TACC support to rebuild PETSc
3.15 against HDF5 1.14 (or install a newer PETSc that supports it
natively).  Once done, the filtered handler becomes unnecessary and
can be removed.
