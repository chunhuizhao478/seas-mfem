# BP5 Debug v52: Zero Dip Traction — Isolating Cross-Component DG Contamination

**Date**: 2026-03-24
**Status**: PRODUCTION RUN SUBMITTED
**Previous**: v51 (dip offset root cause investigation, elastic σ_n, K-matrix coupling diagnostics)
**Branch**: `feature/elasticity`

---

## 1. Motivation

v51 established that the DG discretization produces ~21% cross-component dip/strike traction
contamination (Section 19), and that elastic σ_n feedback partially suppresses it (Section 25).
However, elastic σ_n is an indirect fix — it adjusts friction strength to absorb spurious dip,
rather than eliminating the spurious dip at its source.

This version takes the direct approach: **zero τ_dip after ComputeTraction**.

### 1.1 Why This Is Benchmark-Correct

For BP5 on a planar vertical fault with pure strike-slip loading:
- SCEC spec Eq. 15b: V_3 = 0 (zero dip velocity) outside Ω_f
- SCEC spec Eq. 16: V_3 = V_zero = 1e-20 m/s (negligibly small)
- τ_dip = 0 analytically at all times, everywhere on the fault
- Even Tandem's 0.2-0.7% dip/strike ratio is a discretization artifact

The `--zero-dip-traction` flag is not a workaround — it enforces the analytically correct
value that the DG discretization fails to preserve.

### 1.2 What This Tests

By zeroing only the elastic dip traction while leaving everything else unchanged, we isolate
whether the DG cross-component contamination is responsible for:

1. **Dip offset** — should go to exactly zero (no source of spurious dip)
2. **Interseismic strike deviations** — τ_strike was ~1 MPa low; the cascade
   (Section 4 of v51) predicts this is caused by dip contamination weakening the fault
3. **Earthquake timing drift** — events 2+ diverge from Tandem

If results improve on all three, the root cause is confirmed and we can pursue a proper
fix (e.g., reformulating the DG traction recovery to preserve component decoupling).

---

## 2. Code Changes (v51 → v52)

### 2.1 Removed: Single Elastic Solve Diagnostic

The `--single-solve-uz` flag and its entire implementation block (~100 lines) were removed
from `bp5_verification_full.cpp`. This was a one-time diagnostic that prescribed a slip
pattern, performed a single elastic solve, and dumped u_z/u_x statistics. The diagnostic
served its purpose in v51 Section 19 (confirming u_z contamination) and is no longer needed.

Removed:
- Variable declaration: `bool single_solve_uz`
- CLI parsing: `--single-solve-uz` flag
- Implementation block: slip prescription, elastic solve, MPI reduction, statistics output

### 2.2 Removed: K-Matrix Coupling Test

The `test_k_matrix_coupling.cpp` Makefile entries were removed. This unit test was a v51
diagnostic for verifying cross-component coupling in the stiffness matrix. The test confirmed
that K-matrix components match Tandem (v51 Section 22), so the build rules are no longer
needed.

Removed from Makefile:
- `TEST_K_COUPLING_SRC` / `TEST_K_COUPLING_OBJ` variable definitions
- `seas_test_k_coupling` build target
- Compilation rule for `TEST_K_COUPLING_OBJ`

### 2.3 Retained: `--zero-dip-traction` Implementation

The zero-dip-traction logic in `seas_operator.hpp` (lines 255-262) is unchanged:

```cpp
// v51: Zero dip traction component (index 0 of each DOF's [dip, strike] pair)
if (zero_dip_traction_)
{
   for (int i = 0; i < traction_.Size() / 2; i++)
   {
      traction_(2 * i) = 0.0;  // tau_dip = 0
   }
}
```

This executes after `ComputeTraction()` and before `ComputeRHS()`, so:
- The elastic solve runs normally (3D coupling intact)
- Only the dip component of the resulting traction is zeroed
- `tau_total_dip = tau_pre_dip + 0 ≈ 0` (tau_pre_dip is O(1e-11) Pa, negligible)
- The friction law then sees effectively zero dip stress → V_dip ≈ V_zero

---

## 3. Production Run Configuration

**Job**: `bp5_v52_zero_dip_prod.sbatch`

```
Mesh:       bp5_tandem.msh (1000m reference mesh)
Order:      p=2
DG method:  IP (SIPG)
Solver:     MUMPS-BLR (tol=1e-12)
BC mode:    FarField (default)
t_final:    5.68e10 s (~1800 years)
Nodes:      8 (400 MPI ranks)
Queue:      normal (48hr)
Flags:      --zero-dip-traction --check-residual --write-every-step
Output:     bp5/results_v52_zero_dip_prod/
```

This matches the v50 production baseline exactly, with only `--zero-dip-traction` added.
No `--elastic-sigma-n` — we test each fix independently.

---

## 4. Expected Outcomes

### 4.1 If Zero Dip Traction Fixes All Issues

- Dip slip offset → 0 at all stations
- Dip V → O(1e-20) (matches V_zero)
- τ_strike recovers ~1 MPa deficit
- Earthquake timing matches Tandem through multiple cycles
- State variable ψ matches (no longer depressed by spurious |V| from dip)
- **Conclusion**: DG cross-component contamination was the sole root cause

### 4.2 If Zero Dip Traction Fixes Dip But Not Strike

- Dip quantities correct
- τ_strike still ~1 MPa low, timing still drifts
- **Conclusion**: Strike deficit has a separate source (e.g., traction magnitude error,
  boundary loading, or penalty-related stress redistribution)

### 4.3 If Zero Dip Traction Has Minimal Effect

- Dip should still improve (we're removing the source directly)
- If it doesn't, something else feeds dip (e.g., slip evolution, state feedback)
- **Conclusion**: Unlikely — would indicate a bug in the zeroing implementation

---

## 5. Comparison Plan

When the production run completes, compare against:

1. **Tandem reference**: `bp5/benchmark_data/bp5-qd-*.txt` (8-column SCEC format)
2. **v50 production**: `bp5/results_v50_prod/` (baseline without any dip fix)
3. **v51d elastic σ_n**: if available, to compare direct zeroing vs. indirect feedback

Key metrics per station:
- Dip slip offset (should be ~0)
- Dip/strike traction ratio (should be ~0)
- τ_strike absolute value vs. Tandem
- Earthquake recurrence interval
- State variable ψ trajectory

---

## 6. Revision History

| Date | Section | Change |
|------|---------|--------|
| 2026-03-24 | All | Initial v52 document. Removed single-solve and K-coupling diagnostics. Production run submitted with --zero-dip-traction. |
