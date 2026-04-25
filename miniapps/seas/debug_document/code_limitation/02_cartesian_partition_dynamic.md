# Limitation 02 — A1 Cartesian X+Z partition for TPV104 dynamic

Date introduced: 2026-04-25
Scope: TPV104 dynamic-rupture path only.  Triggered by
`--partition-file <path>` flag in `drivers/tpv104_driver.cpp`; partition
files are produced by `tpv104/mesh/build_symmirror_mesh.py
--emit-partition NP`.  BP5 path uses default ParMETIS (not affected).

## What we currently do

When the dynamic driver is launched with `--partition-file <path>`, the
partition is loaded from a sidecar file with format:

```
np <N>
ne <M>
<rank_0>
<rank_1>
...
<rank_{M-1}>
```

The partition is produced by Cartesian tiling along X and Z only, with
the full Y extent assigned to each rank:

```
For each hex (i, j, k):
   nx_ranks, nz_ranks = factor_2d(NP)        # closest 2D factorization
   rx = (i * nx_ranks) // NX
   rz = (k * nz_ranks) // NZ
   rank = rx + rz * nx_ranks                 # j (Y index) ignored
```

Each rank owns a (X, Z)-tile that spans the full Y extent of the
domain.  The y=0 plane is never crossed by a partition boundary.  Each
rank's local domain is internally y-mirror-symmetric.

## Why this fixes the dip-direction MPI pollution

Documented in detail in `tpv104_bulk_asymmetry_followup_2026-04-25_pm.md`
§13.  Summary:

- Default ParMETIS at np ≥ 4 splits the mesh along Y as well as X, Z.
- The y=0 fault face and nearby non-fault interior faces get cut by
  partition boundaries.
- Ghost-cell exchanges between +y and −y ranks introduce per-rank FP
  ordering differences in flux accumulation (non-associativity).
- These compound through every macro-step and get exponentially
  amplified by the rupture front.
- Result: spurious dip-direction slip (~7 mm at t=12 s in production).

Cartesian X+Z partition eliminates this by construction: no Y cut →
no ghost exchange across y=0 → per-rank flux accumulations are
internally mirror-symmetric.  Local verification confirms bit-exact
identical state at np = {1, 2, 4, 8} on the 1000m mirror mesh.

## What this gives up vs. ParMETIS-default

| Property | ParMETIS default | A1 Cartesian X+Z |
|---|---|---|
| Bit-exact mirror at np ≥ 4 | NO | YES |
| 3D load balancing | full (X, Y, Z) | restricted to (X, Z) |
| Optimal cut-edge minimization | yes | no — ignores Y |
| Adaptable to AMR | yes (re-partition on refinement) | requires re-build |
| Works for non-Cartesian meshes | yes | requires structured cell layout |
| Works for non-y=0 faults | yes | NO (TPV104 specific) |

## Specific restrictions imposed by A1

1. **NP must be expressible as nx_ranks × nz_ranks**.  Any composite
   number works (4 = 2×2, 6 = 2×3, 8 = 2×4 or 4×2, etc.); primes are
   one-axis only (e.g., NP=7 → 1×7).  At NP=400, factorization is
   20×20.

2. **Memory per rank scales with full Y-extent**.  With Y = 4 km / dx
   = 200 m → 20 cells in Y per rank.  Total per-rank cells:
   `NX/nx_ranks × 20 × NZ/nz_ranks`.  At NP=400 with our 80×20×80
   grid: `4 × 20 × 4 = 320` hexes = 1920 tets per rank — same as
   default but without 3D balance flexibility.

3. **No re-partition between SEAS phases**.  This is intentional and
   desirable per user directive — partition file is generated once at
   mesh-build time and used for all phases.  But it also means we
   cannot dynamically rebalance load if dynamic phase becomes much
   more expensive than quasi-dynamic.

4. **Limited to flat Y=0 faults**.  For a curved or dipping fault,
   Cartesian X+Z tiling does not respect the fault plane and the
   bit-exact-mirror property is lost.

## What this prevents

- **General fault-geometry simulations** with this partition path.
  Curved faults, dipping faults, branched faults all need a different
  partitioning approach (the G1 Union-Find approach in
  `dynamic/fault_locality_partition.hpp` is partially compatible but
  was shown to be insufficient — see `tpv104_bulk_asymmetry_followup_*.md` §13).

- **AMR with mid-simulation re-partitioning**.  Any refinement that
  changes the cell count would invalidate the sidecar partition file.

- **Heterogeneous load balancing**.  PETSc/HYPRE for the BP5 quasi-
  dynamic phase usually benefits from 3D-balanced partitions for
  optimal preconditioner convergence; we trade some BP5 efficiency
  for dynamic-phase correctness.  (BP5 path is not currently using
  this partition file, so this is a future concern.)

## Path to removal

### Step 1 — generalize G1 to enforce no-cut on a user-specified plane

Currently G1 (Union-Find on fault adjacencies) is implemented but
insufficient — non-fault interior faces near the fault still get cut.
Extend it to add a "no-cut plane" constraint: any interior face whose
two adjacent elements lie on opposite sides of the plane (regardless
of being a fault face) gets an "infinite" merge weight.

For TPV104 with y=0 fault: this would merge BOTH fault and non-fault
faces crossing y=0, replicating A1's effect via a more general API.

### Step 2 — switch to per-face fault interface (BP5-style) for dynamic

Per-face friction representation makes the fault Riemann partition-
robust by averaging.  Partition asymmetry no longer directly drives
V_dip via per-side Riemann formula; instead it gets integrated into
per-face traction noise that doesn't compound exponentially.

This eliminates the requirement for any specific partition geometry —
ParMETIS default works fine.

### Step 3 — symmetric ParMETIS via edge weights + quality

Implement true mirror-aware partitioning by setting infinite edge
weights on fault adjacencies AND on edges crossing the mirror plane.
Pass this weighted graph to METIS_PartGraphKway directly.

## When to remove

- For curved/dipping fault problems: Step 1.
- For SEAS unification with BP5: Step 2.
- For optimal load balance with arbitrary fault geometry: Step 3.

For the current TPV104 verification milestone, A1 is sufficient and
the simplest implementation.

## Verification record (local, 2026-04-25)

4 × 4 × 8 km / 1000 m mirror mesh, tfinal = 2 s, ADER-O2:

```
WITHOUT --partition-file (default ParMETIS):
  np=1:  slip_dip = -3.6395e-02   slip_strike = 4.9436
  np=2:  slip_dip = -3.6395e-02   slip_strike = 4.9436   (matches np=1)
  np=4:  slip_dip = +1.1920e-01   slip_strike = 4.7393   (DIFFERS)
  np=8:  slip_dip = +1.1920e-01   slip_strike = 4.7393   (DIFFERS)

WITH --partition-file (Cartesian X+Z, full Y per rank):
  np=1:  slip_dip = -3.63949181830000e-02   slip_strike = 4.94355048
  np=2:  slip_dip = -3.63949181830000e-02   slip_strike = 4.94355048
  np=4:  slip_dip = -3.63949181830000e-02   slip_strike = 4.94355048
  np=8:  slip_dip = -3.63949181830000e-02   slip_strike = 4.94355048

  V_max identical at every macro-step: 12.6092 m/s.
```

Bit-exact identical to 14-digit display precision.

## Frontera production sbatch reference

See `jobs/tpv104/tpv104_symmirror_200m_partfile.sbatch` (added in same
commit as this doc).  Generates the partition file in-job alongside
the mesh and uses it via `--partition-file`.
