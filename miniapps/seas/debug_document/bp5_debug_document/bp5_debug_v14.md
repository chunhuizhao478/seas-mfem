# BP5 Debug v14: Boundary Condition + Output Floor Fix

**Date**: 2026-03-13
**Status**: Analysis complete, fixes designed
**Previous**: v13 (H1-H4, H7 implemented)

---

## 1. Executive Summary

Analysis of BR2 1000m results (post-v13 fixes) reveals two remaining issues:

1. **Wrong boundary conditions** (H8): Dirichlet BC applied on ALL non-free-surface boundaries (attrs 1-4,6), but only the bottom (attr 6) should be Dirichlet. This over-constrains the domain, reducing effective loading to ~36% of Vp and producing a 3.7× longer recurrence interval.

2. **Output floor at 1e-30** (H9): `bp5_benchmark_output.hpp` clamps V to 1e-30 before log10, creating a sharp artificial floor at dp+10 stations where V genuinely drops below 1e-30.

---

## 2. BR2 1000m Results Analysis

### 2.1 MFEM Results (post-v13)

| Metric | MFEM BR2 1000m | Tandem ref | Ratio |
|--------|---------------|------------|-------|
| Events in 1800 yr | 3 | 8 | 0.38 |
| Recurrence interval | ~749 yr | ~240 yr | 3.1× |
| Long-term slip rate | ~11.5 mm/yr | ~31.6 mm/yr (Vp) | 0.36 |
| Peak V_strike | ~0.03 m/s | ~0.65 m/s | 0.05 |
| V floor at dp+10 | -30 (artificial) | no floor | — |

### 2.2 Station-Level V Floor Analysis

Minimum log10(V) across all stations:

| Station | min log10(V_strike) | min log10(V_dip) | Hits floor? |
|---------|--------------------|--------------------|-------------|
| strk+00dp+00 | -10.68 | -20.00 | No (V_dip = V_zero) |
| strk+16dp+00 | -10.61 | -20.00 | No |
| strk-16dp+00 | -10.60 | -20.00 | No |
| strk+36dp+00 | -10.14 | -20.00 | No |
| strk-36dp+00 | -10.14 | -20.00 | No |
| strk+00dp+22 | -10.04 | -20.00 | No |
| **strk+00dp+10** | **-30.00** | **-30.00** | **Yes** |
| **strk+16dp+10** | **-30.00** | **-30.00** | **Yes** |
| **strk-16dp+10** | **-30.00** | **-30.00** | **Yes** |
| **strk-24dp+10** | **-30.00** | **-30.00** | **Yes** |

All dp+10 stations (10km depth, deep VW zone) hit the 1e-30 floor. dp+00 stations never hit it.

---

## 3. Root Cause: Boundary Conditions (H8)

### 3.1 Current Code (Buggy)

`elasticity_operator.hpp:SetupBoundaryMarkers()` lines 208-228:

```cpp
// SCEC BP5-QD: all non-free-surface boundaries get Dirichlet
for (int be = 0; be < mesh_.GetNBE(); be++)
{
   int attr = mesh_.GetBdrAttribute(be);
   if (attr == 1 || attr == 2 || attr == 3 || attr == 4 || attr == 6)  // WRONG
   {
      dirichlet_bdr_marker_[attr - 1] = 1;
   }
}
```

This applies Dirichlet u = (0, sgn(x)·Vp·t/2, 0) on:
- attr 1: x = -Lx (far-field left)
- attr 2: x = +Lx (far-field right)
- attr 3: y = +Ly (along-strike +)
- attr 4: y = -Ly (along-strike -)
- attr 6: z = Lz (bottom)

### 3.2 Correct Boundary Conditions

Per SCEC BP5-QD spec, Tandem implementation, and the project's own
`tandem_boundary_condition_analysis.md` (Section 7.2, Full Domain table):

| Boundary | Attr | Correct BC |
|----------|------|-----------|
| x = -Lx (left) | 1 | **Natural** (zero traction) |
| x = +Lx (right) | 2 | **Natural** (zero traction) |
| y = +Ly | 3 | **Natural** (zero traction) |
| y = -Ly | 4 | **Natural** (zero traction) |
| z = 0 (free surface) | 5 | **Natural** (zero traction) |
| z = Lz (bottom) | 6 | **Dirichlet** u = (0, sgn(x)·Vp·t/2, 0) |

Only the bottom boundary should drive plate loading. The x/y boundaries should be
free (zero traction), allowing the domain to deform naturally under the bottom loading.

### 3.3 Why This Explains the Results

With Dirichlet on all 5 boundaries, the displacement field is pinned everywhere.
The only "free" boundary is the top (z=0, free surface). The fault sees much less
differential loading because the x-boundaries force displacement = sgn(x)·Vp·t/2
throughout the domain, leaving no stress concentration at the fault.

The effective loading rate reduces to ~36% of Vp (matching the observed
11.5/31.6 = 0.36 ratio), directly explaining:
- 3.1× longer recurrence (749/240 ≈ 3.1)
- 20× lower peak V (less energy release per event)
- Fewer total events (3 vs 8 in 1800 yr)

### 3.4 Tandem Confirmation

Tandem's `bp5.lua` boundary function:
```lua
function BP5:boundary(x, y, t)
   local Vh = self.Vp / 2.0 * t
   if y > 1.0 then return {Vh, 0, 0}
   elseif y < -1.0 then return {-Vh, 0, 0}
   else return {self.Vp * t, 0, 0} end
end
```

This is only applied to faces marked `boundary_linear = true`, which in Tandem's
BP5 mesh corresponds to the **bottom boundary only**. All other boundaries
(x±, y±) use Natural BC.

---

## 4. Root Cause: Output Floor (H9)

### 4.1 Current Code

`bp5_benchmark_output.hpp` WriteFromGlobalData (lines 253-254, 261):
```cpp
real_t V_dip    = std::max(std::abs(global_V_dip(dof)), 1e-30);
real_t V_strike = std::max(std::abs(global_V_strike(dof)), 1e-30);
real_t th = std::max(global_theta(dof), 1e-30);
```

WriteRow path (lines 358-359, 367): same pattern.

### 4.2 Why dp+10 Stations Hit the Floor

At 10km depth in the VW zone, after an earthquake the fault is deeply locked.
The state variable psi increases (strengthening), driving V to very small values.

For psi/a ≈ 60 (typical post-earthquake):
- V ≈ 2V₀ · (τ/σ_n) / (a · exp(psi/a)) ≈ 2×10⁻⁶ × 0.28 / (0.004 × 10²⁶) ≈ 10⁻³¹

So V genuinely drops below 1e-30 at locked stations. The output clamping clips
these to exactly -30 on the plot.

### 4.3 Tandem Comparison

Tandem writes **raw V values** to output files (no log10, no clamping).
The 1e-30 floor is entirely MFEM-specific. The log10 transformation and any
floor handling is done by post-processing scripts.

### 4.4 Relationship to H7

The H7 fix (changing `tau_abs < 1e-30` to `tau_abs <= 0.0` in the friction solver)
was correct but insufficient. It removed the floor in the Brent solver's
zero-traction guard, but the output code still clamps independently. The user
confirmed H7 alone did not fix the sharp cutoff — the output floor is the cause.

---

## 5. Fixes

### H8 — Boundary condition: Dirichlet only at bottom (z=Lz)

**File**: `miniapps/seas/domain/elasticity_operator.hpp`

**Change 1** — `SetupBoundaryMarkers()` (line 223):
```cpp
// BEFORE (buggy):
if (attr == 1 || attr == 2 || attr == 3 || attr == 4 || attr == 6)

// AFTER (correct):
if (attr == 6)
```

**Change 2** — Update comments at lines 45, 210-218, and 1340-1341 to reflect
that only the bottom boundary is Dirichlet.

No changes needed to `AssembleDirichletLoading()` logic — it already uses
`dirichlet_bdr_marker_` to filter boundary elements.

### H9 — Remove 1e-30 output floor

**File**: `miniapps/seas/io/bp5_benchmark_output.hpp`

Replace `std::max(V, 1e-30)` with proper zero-handling in both WriteFromGlobalData
and WriteRow:

```cpp
// BEFORE:
real_t V_dip    = std::max(std::abs(global_V_dip(dof)), 1e-30);
real_t V_strike = std::max(std::abs(global_V_strike(dof)), 1e-30);
real_t th = std::max(global_theta(dof), 1e-30);

// AFTER:
real_t V_dip    = std::abs(global_V_dip(dof));
real_t V_strike = std::abs(global_V_strike(dof));
real_t th = global_theta(dof);
```

Then in row construction, handle zero before log10:
```cpp
V_strike > 0.0 ? std::log10(V_strike) : -300.0,
V_dip > 0.0 ? std::log10(V_dip) : -300.0,
th > 0.0 ? std::log10(th) : -300.0
```

**File**: `miniapps/seas/tests/verification/bp5_verification_full.cpp`

Lines 557, 727: Replace `std::max(V_max, 1e-30)` with same pattern.

---

## 6. Implementation Status

| Fix | Issue | Status | Files |
|-----|-------|--------|-------|
| H1 | PsiToTheta per-DOF Dc | ✅ DONE (v13) | `dieterich_ruina.hpp`, `rate_state_fault.hpp` |
| H2 | Output intervals to SCEC spec | ✅ DONE (v13) | `bp5_benchmark_output.hpp` |
| H3 | IP tensor-coupled penalty | ✅ DONE (v13) | `elasticity_operator.hpp` |
| H4 | dt_max 0.5yr → 0.1yr | ✅ DONE (v13) | `bp5_verification_full.cpp`, `time_stepper.hpp` |
| H5 | V_nuc mismatch | Info only | — |
| H6 | Pre-stress formula | Info only | — |
| H7 | tau_abs floor in solver | ✅ DONE (v13) | `dieterich_ruina.hpp` |
| **H8** | **BC: Dirichlet only at bottom** | **TODO** | `elasticity_operator.hpp` |
| **H9** | **Output 1e-30 floor** | **TODO** | `bp5_benchmark_output.hpp`, `bp5_verification_full.cpp` |

---

## 7. Expected Impact

### After H8 (boundary fix):

| Metric | Current | Expected | Tandem ref |
|--------|---------|----------|------------|
| Recurrence | ~749 yr | ~240 yr | ~240 yr |
| Long-term slip rate | 11.5 mm/yr | ~31.6 mm/yr | Vp |
| Peak V | 0.03 m/s | ~0.6 m/s | 0.65 m/s |
| Events in 1800 yr | 3 | ~7-8 | 8 |

### After H9 (output fix):

dp+10 station V values will show true slip rates (potentially -40 to -60 in log10)
instead of clipping at exactly -30.

---

## 8. Verification Plan

```bash
conda activate mfem-dev
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas
make -j seas_bp5_full
mpirun -np 4 tests/verification/seas_bp5_full --inline-mesh --tfinal 31557600 --write-every-step
```

1. Compile succeeds
2. First event at ~240 yr (not ~749 yr)
3. Long-term slip rate ≈ Vp = 31.6 mm/yr
4. Peak V ~0.5-0.7 m/s
5. dp+10 stations: no flat line at -30
6. Station time series qualitatively match Tandem
