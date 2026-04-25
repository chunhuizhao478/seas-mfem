# Code limitations & accepted compromises

This folder records temporary compromises in the SEAS-MFEM dynamic-
rupture (TPV104) implementation that were necessary to reach a working
benchmark comparison against SeisSol but are NOT general-purpose.
Each compromise has a clear path to removal, captured below.

Index:

- [01_structured_mesh_dynamic.md](01_structured_mesh_dynamic.md) — TPV104 dynamic
  uses a structured Cartesian tetrahedral mesh, not an unstructured
  graded mesh.  Limits production-mesh generation flexibility.
- [02_cartesian_partition_dynamic.md](02_cartesian_partition_dynamic.md) — TPV104
  dynamic uses an explicit Cartesian X+Z partition (full Y per rank)
  instead of ParMETIS-default partitioning.  Required for bit-exact
  y-mirror; limits load-balancing flexibility.
