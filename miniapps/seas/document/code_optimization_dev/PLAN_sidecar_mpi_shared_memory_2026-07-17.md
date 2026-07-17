# PLAN — Per-node MPI-3 shared-memory windows for the gridded sidecars

**Date:** 2026-07-17
**Status:** Phases 0–2 + §6 tests IMPLEMENTED (2026-07-17); Phase 3 (flip the
Expanse benchmark to `--sidecar-shared-mem` + re-measure) pending on-cluster.
New test gate: `make test-data-field-3d-shared-mem` (np=2 + np=4).
**Owner folder:** `document/code_optimization_dev/`
**Motivating context:** the v4_0_0 ALT CASE1 MFEM-vs-SeisSol speed benchmark
(`jobs/safs/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_small_p{1,2,3}`) had to
raise `--mem` from the SeisSol deck's `140000M` to `249000M` because MFEM is
pure-MPI and replicates the sidecars **per rank**.

---

## 1. Problem

`seas_spatial_dyn_driver` reads three gridded sidecars — material (CVM velocity),
on-fault stress, and rate-state friction — as `data_projection_v1` HDF5. Each is
loaded by `DataField3D`, which holds the full grid in a **per-rank** member
`std::vector<real_t> data_` (`io/data_field_3d.hpp:151`), filled once in the ctor
(`io/data_field_3d.cpp:325-340`) and **read-only thereafter**.

Because the fill is per-*rank* (each rank does its own `H5Fopen`/`H5Dread`), a node
running `R` MPI ranks holds `R` identical copies. For the SAFS production run:

| Sidecar | `DataField3D` grids | Size (this deck) | Residency |
|---|---|---|---|
| Material (Vp/Vs/ρ) | 3 | ~760 MB | **whole run** (`vel_bundle`, `drivers/spatial_dyn_driver.cpp:1332`) |
| Stress (σ_xx…σ_xz) | 6 (via `StressField3D`) | ~256 MB | **transient** (freed after `ComputeParams`) |
| Friction (rs_a/rs_srW) | 2 | ~87 MB | setup → held through resolver |

- Whole-run resident ≈ **847 MB/rank** (material + friction); setup **peak** ≈ 1.1 GB/rank.
- On a 128-rank node: **~108 GB** just for whole-run grids (≈141 GB at setup peak).
- That is what forced `--mem=249000M`; it also caps how many ranks fit per node on
  memory-tighter partitions.

This is pure duplication of **read-only** data. SeisSol solves the identical problem
with ASAGI's per-node MPI-shared mode (`SEISSOL_ASAGI_MPI_MODE`); MFEM has no analog
today (repo-wide grep for `MPI_Win`/`MPI_COMM_TYPE_SHARED` → **zero matches**).

## 2. Goal / non-goals

**Goal.** Hold **one physical copy of each grid per node**, shared read-only by all
node-local ranks via an MPI-3 shared-memory window, cutting per-node sidecar memory
from `R × 847 MB` to `~847 MB` — so the benchmark returns to `--mem=140000M` at full
128 ranks/node, and memory stops bounding ranks-per-node.

**Non-goals.** No change to numerics, interpolation, values, or the schema. No change
to the compute path (this is I/O-layer only). Not about read *speed* — see §3.

## 3. Why this is the right, low-risk fix (and OpenMP is not needed for memory)

`DataField3D::Evaluate` is called **once per DOF during one-time setup projections,
never inside the time loop** (grounded):

- stress → per-fault-DOF loop in `FaultGeometry::ComputeParams`
  (`fault/fault_geometry_safs_templated.inl:88-95`), once from `ApplyCsmStressSidecar`;
- material → wrapped in MFEM `Coefficient`s (`io/material_coefficients.cpp:21-43`)
  consumed at **operator-build** time (the stiffness/flux pool is "assembled once and
  reused", per `CLAUDE.md`); never touched in `wave_operator.inl`'s ADER loop;
- friction → per-fault-DOF `SpatialFrictionResolver::ResolveRateState`
  (`drivers/spatial_dyn_driver.cpp:2097`).

Consequences: (a) the win is **pure footprint**, so shared-window NUMA/read-latency is
irrelevant — reads happen once at setup; (b) the data being **write-once/read-only**
means **no race conditions** and no synchronization in the hot path; (c) the compute
model stays pure-MPI (unchanged), so this is orthogonal to — and far cheaper than —
the OpenMP-hybrid effort (see `document/gpu_dev/PLAN_openmp_hybrid_2026-07-17.md`).
OpenMP would reduce this memory only as a side effect of running fewer ranks; this
plan removes the duplication directly.

## 4. Design

Replace `DataField3D`'s owning `std::vector<real_t> data_` with a **pointer + length
into an MPI-3 shared-memory window** that is allocated once per node:

- Split the caller's communicator into node-local groups:
  `MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm)`.
- **Node-local rank 0** calls `MPI_Win_allocate_shared(Nx*Ny*Nz*sizeof(real_t), …,
  node_comm, &baseptr, &win)`; every other node-local rank allocates size 0 and
  obtains the shared base with `MPI_Win_shared_query(win, 0, &sz, &disp, &baseptr)`.
- **Only node-local rank 0** does the `H5Fopen`/`H5Dread` fill into `baseptr` (with the
  same float64→`real_t` conversion currently at `io/data_field_3d.cpp:337-340`), then
  an `MPI_Barrier(node_comm)` so readers see a completed buffer before first use.
- `data_` becomes `const real_t* data_; size_t data_len_;` (or a thin non-owning span).
  `flat_index_`/`evaluate_trilinear_`/`evaluate_catmull_rom_` already index `data_[…]`
  uniformly (`io/data_field_3d.cpp:424-459,524-618`) — the read code is **unchanged**.
- The window handle is owned by the `DataField3D` instance; its dtor does
  `MPI_Win_free(&win)` (frees the shared mapping on all node-local ranks; MFEM tears
  these down before `MPI_Finalize`).
- The tiny per-rank members (`x_,y_,z_` axes at `io/data_field_3d.hpp:150`, `bbox_`,
  scalars) **stay replicated** — they are kilobytes and every rank needs them for
  `find_index_`/`ContainsBBox`; only the big `data_` payload is shared.

### Value-range + NaN checks (defense-in-depth) stay
The post-load sanity loop (`io/data_field_3d.cpp:347-372`) must run on **node-local
rank 0 only** (it owns the fill), before the barrier; a failure there must
`MPI_Abort` so all ranks stop (an `MFEM_ABORT` on one rank alone would hang the
others at the barrier). This is the one new correctness subtlety.

### Backward-compatible constructor
Add an **optional** node-local communicator to the ctor
(`io/data_field_3d.hpp:85-87`), defaulting to `MPI_COMM_NULL` = **the current
serial per-rank path** (one owning `std::vector`, no window). This preserves every
unit test verbatim (they construct readers with no comm — `tests/unit/
test_data_field_3d.cpp`) and the standalone projector tools. `StressField3D`'s ctor
(`io/stress_field_3d.hpp:67`) forwards the comm to its six `DataField3D`.

## 5. Phased implementation

**Phase 0 — feature flag + comm plumbing (no behavior change).**
Add `--sidecar-shared-mem` CLI flag (default OFF) and, in the driver, build a
node-local comm once: `MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, …)`
near the comm setup (`drivers/spatial_dyn_driver.cpp:750-754`). Thread it (optional
arg) into the three loader sites:
- material `spatial::LoadSpatialVelocityBundle` → `spatial/code/spatial_velocity.cpp:75-77`
  (already takes `ParMesh&`; add the node-comm arg);
- friction `MakeRateStateSidecarFields` → `drivers/spatial_dyn_driver.cpp:549`
  (free function; `comm` in local scope);
- stress `ApplyCsmStressSidecar`/`StressField3D` → `spatial/code/spatial_stress.cpp:58`
  (add a comm param to the call chain).
With the flag OFF, everything constructs `MPI_COMM_NULL` → identical to today.

**Phase 1 — `DataField3D` dual-mode storage.**
Implement the window path behind the ctor comm arg: `MPI_COMM_NULL` → today's
`std::vector` owner; a real node-comm → shared-window path (§4). Unify reads on a
`const real_t*` so `Evaluate` is storage-agnostic. Gate with `test-data-field-3d`
run at `np=2` and `np=4` on one node (see §6).

**Phase 2 — `StressField3D` forwarding + transient free.**
Forward the comm through the six sub-fields (`io/stress_field_3d.cpp:45-53`). Because
stress is transient (freed after `ComputeParams`), its six windows can be freed early
— confirm the `StressField3D` local in `apply_csm_impl` (`spatial/code/
spatial_stress.cpp`) destructs (and thus `MPI_Win_free`s) after the projection loop.

**Phase 3 — flip the benchmark + measure.**
Set `--sidecar-shared-mem` in the three v4_0_0 sbatch and drop `--mem` back to
`140000M` (or lower). Confirm: identical `[derived]`/fault-projection output vs the
non-shared run (byte-for-byte — same reads, same values), and per-node RSS drops from
~108 GB to ~1 GB of sidecars. Then make the flag the default once proven.

## 6. Testing / regression gate

The values are unchanged, so the gate is **equality vs the current path** + a
multi-rank memory/consistency check. Build with `conda activate mfem-dev`.

- `make test-data-projection` (`Makefile:4105`) → `test-data-field-3d`,
  `test-stress-field-3d`, `test-field-projector`, `test-heterogeneous-material`.
- `make test-project-fault-prestress` (`Makefile:4180`),
  `make test-spatial-friction-sidecar` (`4330`),
  `make test-spatial-velocity-bundle` (`4336`).
- **New tests to add:** (a) a `np≥2` single-node test that constructs a `DataField3D`
  with a node-comm and asserts `Evaluate` at a grid of points is **bit-identical** to
  the `MPI_COMM_NULL` path (same file); (b) an assertion that non-zero-local-rank
  windows report the same `baseptr` contents (shared, not copied) — e.g. rank 0 writes
  a sentinel pre-barrier is *not* possible (read-only), so instead compare a checksum
  of `data_` across node-local ranks == equal.
- Note `test-spatial-stress-bundle` is **currently a no-op skip**
  (`Makefile:4339-4346`, depends on a segfaulting BP5 fixture) — do **not** rely on it;
  cover stress via `test-stress-field-3d` + `test-project-fault-prestress`.

Because Phase 3 asserts byte-for-byte equality of the projected fields, this change is
**not** subject to the float-reassociation tolerance problems that gate the OpenMP/
cached compute paths — the arithmetic is identical; only the buffer's owner changed.

## 7. Risks & mitigations

| Risk | Mitigation |
|---|---|
| One rank `MFEM_ABORT`s in the load sanity check → others hang at the node barrier | Load/validate on node-rank-0 only; on failure `MPI_Abort(MPI_COMM_WORLD)`, never a single-rank abort. |
| MPI-3 shared-mem unsupported / heterogeneous-memory node | Feature-flagged (default OFF initially); `MPI_Win_allocate_shared` is MPI-3.0 (OpenMPI 4.1.x on Expanse supports it). Fallback = the `MPI_COMM_NULL` per-rank path. |
| Window lifetime vs `MPI_Finalize` | (Reworded per REVIEW.md R-001.)  The LONG-LIVED readers (`vel_bundle`'s 3 windows, the friction readers) CANNOT be destroyed before `MPI_Finalize`: they are main-scope locals, the driver finalizes before main unwinds, and the R-008 ordering contract (`vel_bundle` outlives `wave_ptr`) forbids an early reset.  Actual mechanism: the dtor checks `MPI_Finalized()` and skips the collective `MPI_Win_unlock_all`/`MPI_Win_free` after finalize (the OS reclaims the shm segment at process exit; SPMD flow keeps the decision consistent across node ranks).  Formally the MPI standard wants windows freed pre-finalize; OpenMPI tolerates this, and the path is covered by the permanent np=1 regression `seas_test_data_field_3d_shared_mem_finalize`.  The TRANSIENT stress windows (6) free normally inside `apply_csm_impl`. |
| `real_t` size vs on-disk `float64` | Window sized in `real_t` (matches today's converted `data_`), fill converts as now (`:337-340`). |
| Serial/test path regressions | Default `MPI_COMM_NULL` preserves the exact current path; unit tests unchanged. |

## 8. Effort & decision

**Scope:** ~1 core file (`io/data_field_3d.{hpp,cpp}`), `StressField3D` forwarding,
3 loader call-site signature tweaks, 1 driver flag + node-comm, 1–2 new tests.
**Estimate:** ~1–2 focused days incl. tests and the benchmark re-run. **Risk: low**
(read-only data, feature-flagged, equality-gated).

**Recommendation:** do this. It is the correct, targeted fix for the memory problem
that motivated the `--mem` bump, keeps the fast pure-MPI compute untouched, and is the
MFEM analog of SeisSol's own ASAGI-MPI-shared mode. It does **not** require or block
the OpenMP-hybrid work — if the goal is memory, this alone suffices.
