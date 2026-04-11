# BP5 Debug v57/v58: MPI Core-Count Dependence — RESOLVED

**Date:** 2026-04-01 to 2026-04-02
**Issue:** Simulation results depend on number of MPI ranks (100 vs 400)
**Status:** FIXED (commit 932b797)
**Verified:** 100-rank and 400-rank V_max trajectories match; results agree with Tandem

## 1. Problem Statement

Two runs of the same BP5 v57 configuration (1000m, IP, p=1, PETSc TS RK45, MUMPS)
produced different results depending on MPI rank count. Divergence grew from ~1% at
step 1 to ~7% by step 410.

## 2. Root Cause

**One-sided shared fault face detection in tag-based path.**

Internal boundary elements (attr=3 in Gmsh) only exist on ONE rank per shared face
in MFEM's ParMesh. `BuildFaultTaggedFaces()` iterates boundary elements to detect
fault faces, so only the rank owning the boundary element detects the shared face.
The other rank never processes it in `AssembleSlipContributionIPShared`, silently
dropping that element's RHS contribution (`elvec2`).

**Code location:** `elasticity_operator.hpp:508-539` (BuildFaultTaggedFaces)

**Mechanism:**
- Rank A has the boundary element → detects shared fault face → contributes `elvec1`
- Rank B does NOT have the boundary element → never sees the face → `elvec2` is lost
- With more MPI ranks → more shared fault faces → more missing contributions
- 100 ranks: 4 one-sided shared fault faces
- 400 ranks: 55 one-sided shared fault faces
- Missing contributions cause 0.037% ||b||_2 mismatch → amplified by κ≈160 → 6.4% ||u|| error

## 3. Fix (commit 932b797)

After building `fault_shared_tagged_` from local boundary elements, Allgather the
global face IDs of all detected shared fault faces across all ranks. Each rank then
checks its own shared faces against this global set and adds any it missed.

```cpp
// v58 fix: Exchange shared fault face tags between neighboring ranks.
// Each rank broadcasts its shared fault face GIDs, then all ranks
// check their own shared faces against the global set.
```

Both ranks now detect every shared fault face, ensuring both `elvec1` and `elvec2`
contributions are assembled in `AssembleSlipContributionIPShared`.

## 4. Verification

- 100-rank and 400-rank V_max trajectories match after fix
- Results agree with Tandem reference (1000m p=1)
- Micro-test confirmed raw face integrator is bit-identical for interior vs shared

## 5. Discarded Hypotheses

| Hypothesis | Fix Attempted | Result | Why Discarded |
|-----------|---------------|--------|---------------|
| Shared fault DOF duplication | Owned-fault layout (f94d000) | No effect | Only 4 shared fault face instances |
| DOF-within-face permutation | Canonical DOF perm (978155d) | No effect | Permutation doesn't affect assembly |
| Shared-face parameterization | KeepNbrBlock, skip_zeros | No effect | K matrix values identical (5e-15 rel) |
| Face integrator orientation | Micro-test | Bit-identical | Raw integrator is correct |
| FaultBasis tangent mismatch | SLIP-EMBED diagnostic | Consistent | Same tangent frame for interior/shared |

## 6. Diagnostic Data That Led to Root Cause

| Quantity | 100 ranks | 400 ranks |
|----------|-----------|-----------|
| `||slip||_inf` | 2.00e-05 | 2.00e-05 |
| `sum|K_ij|` | 8.1196e+21 | 8.1196e+21 |
| `||b||_2` | 2.3240e+11 | **2.3232e+11** |
| `||b_interior||_2` | 2.3221e+11 | 2.3204e+11 |
| `||b_shared||_2` | 9.42e+09 | 1.14e+10 |
| `||u||_inf` | 6.85e-05 | **6.41e-05** |
| shared fault faces | 4 | 55 |

Key insight: `||b||` differs pre-solve → not MUMPS, not K matrix. The slip is
identical but the RHS differs → the manual fault RHS scatter path is the suspect.
The shared path only scatters `elvec1` (line 2125), relying on the other rank to
contribute `elvec2`. But the other rank doesn't detect the face.

## 7. References

- Fix commit: 932b797
- Bug location: `elasticity_operator.hpp:BuildFaultTaggedFaces()`
- Scatter path: `elasticity_operator.hpp:AssembleSlipContributionIPShared()`
- MFEM ParMesh boundary element distribution: one rank per internal boundary element
