# Code Review: LTS Phase 4 Step 0 — serial clustering + material-before-ParMesh (P4S0-2, commit 300649a) — 2026-07-19

## Scope
- File: `drivers/spatial_dyn_driver.cpp` (serial clustering block + serial→local reconstruction).
- Plan: §"Phase 4 Step 0", P-009, D-1, D-2.
- Method: one adversarial reviewer, focused on the two load-bearing np>1 assumptions (serial→local element order; rank-0 sidecar-read collectiveness).

**Verified CORRECT:** the serial→local map (MFEM `ParMesh` ctor assigns local elements in ascending serial order among `part[e]==rank` — `pmesh.cpp` `BuildLocalElements`, no post-ctor reorder; conforming meshes); the rank-0 sidecar read is non-collective on the happy path (2-arg `Mesh&` overload → serial HDF5, `node_comm==MPI_COMM_NULL`); the MPI_Bcast sequence is rank-uniform; serial-vs-local id usage (hash=serial, stepper/layout/reorder=local); byte-exact-off gating; np=1 mechanics; lifetime (no read of cleared `smesh`).

## Findings + resolutions

### [F-1] CRITICAL — serial clustering + partition used the RAW `cfg.boundary.fault_attr` (−1 on SAFS) → fault_pairs empty → D-1 & D-2 broken
SAFS decks omit `[boundary]`, so `cfg.boundary.fault_attr = -1`; the driver's fallback to 101 lands in `bc.fault_attr` only *after* the serial block. Passing −1 left `in.fault_pairs` empty, so `BuildLtsClustering` never forced a fault face's two elements into one cluster (D-1) and `BuildLtsAwarePartition`'s fault-lock was a no-op (D-2 — fault faces cut across ranks, resurrecting the shared-fault hazard). Invisible to TPV tests (which set `fault=3`); bites every SAFS mesh (Frontera-only). **FIXED:** compute `fault_attr_eff = (cfg.boundary.fault_attr > 0) ? cfg.boundary.fault_attr : 101` at the top of the serial block (mirroring the fault-locality block's idiom) and use it for the serial clustering + partition; added the nonconforming guard (F-1b).

### [F-4] MODERATE — METIS-failure fallback aborted with a misleading message
On `BuildLtsAwarePartition` throw, `lts_part_vec` was cleared → `part_data` fell to `nullptr` → the ParMesh used MFEM's *internal* partition, which the serial→local map can't reproduce → the count `MFEM_VERIFY` aborted ~280 lines later with a "partition/element-ordering mismatch" message that never mentioned METIS. **FIXED:** the fallback now builds a RECONSTRUCTABLE explicit partition (`FindFaultFaceIndices` + `BuildFaultLocalityPartitioning` = GeneratePartitioning + fault-lock), preserving D-2 and letting the reconstruction work.

### [F-2] MODERATE — uncaught exception in the rank-0 block → hang
A throw in the rank-0 serial block (realistically `bad_alloc` from the ~1.3 GB sidecar) would terminate rank 0 while peers waited at `MPI_Bcast` → deadlock. **FIXED:** wrapped the rank-0 block in `try/catch` → `MPI_Abort(comm, 1)` (mirrors the post-ParMesh material load).

### [F-6] LOW–MODERATE — stale np>1 log
The Phase-1 log said np>1 lts=rate2 runs GTS "byte-identical to lts=off" — no longer true (the reorder + LTS partition are now active at np>1). **FIXED:** the log now fires on `lts_layout_ready` (any np), states the Step-0 machinery is active at np>1, that stepping stays GTS (MPI stepping = Phase 4), and that np>1 parity is Frontera-validated.

### [F-1b] LOW — nonconforming-mesh guard
The serial→local map is a conforming-mesh guarantee. **FIXED:** `MFEM_VERIFY(!smesh.Nonconforming(), ...)`.

## Verdict
PASS WITH FIXES — all applied. Driver builds+links; LTS units green (predictor 31, reorder 14, friction 10, nucleation 62, partition 19); byte-exact-off parity 36/36, tpv102 4/4. np>1 parity (np=2==np=1, symmirror) + partition balance + the rank-0 sidecar-read collectiveness under production HDF5 remain Frontera-only.
