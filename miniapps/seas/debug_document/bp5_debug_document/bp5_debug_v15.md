# BP5 Debug v15: All-Boundary Dirichlet Loading + Configurable Nucleation

**Date**: 2026-03-13
**Status**: Analysis complete, fixes designed
**Previous**: v14 (H8 bottom-only Dirichlet, H9 output floor removal)

---

## 1. Executive Summary

Post-v14 results (H8 bottom-only Dirichlet + H9 output floor fix) show good surface behavior (~250yr recurrence) but the VW zone center (dp+10) barely participates:

- dp+10 total slip: 7m over 1800yr (Tandem: 55m)
- dp+10 interseismic V: 10^(-57) (effectively locked)
- dp+10 loading rate: ~6.3 kPa/yr (Tandem: ~50 kPa/yr)

**Root cause**: H8 (bottom-only Dirichlet) was based on an incorrect interpretation of Tandem's boundary conditions. Re-analysis of Tandem's actual source code (`bp5.lua`, `bp5.geo`, `bp5.toml`) reveals that Tandem applies Dirichlet loading on **ALL non-fault boundaries** (top, bottom, and all four sides), not just the bottom.

**Secondary issue**: SCEC spec uses V_nuc=0.03 + delta_tau=eta*V_nuc in the nucleation zone, while Tandem uses V_nuc=0.01 with no delta_tau. This produces an oversized first event (6.37m vs Tandem's ~3m), creating a 201yr stress shadow that delays VW zone recovery.

---

## 2. Post-v14 Results Analysis

### 2.1 Surface Stations (dp+00) - Reasonable

| Metric | MFEM (post-v14) | Tandem |
|--------|-----------------|--------|
| Recurrence | ~250 yr | ~240 yr |
| Number of events | ~7 | ~8 |
| Peak V_strike | ~0.03 m/s | ~0.65 m/s |

Surface behavior is qualitatively correct but dp+10 is disconnected.

### 2.2 Depth Stations (dp+10) - Problem

| Metric | MFEM (post-v14) | Tandem |
|--------|-----------------|--------|
| Total slip (1800yr) | 7m | ~55m |
| Interseismic V | 10^(-57) | 10^(-9) |
| Event participation | Barely | Full |
| Loading rate | ~6.3 kPa/yr | ~50 kPa/yr |
| Psi (interseismic) | ~1.0 | ~0.6 |
| Friction strength | 24.3 MPa | — |
| Actual stress | 11.46 MPa | — |

The friction strength (24.3 MPa) far exceeds the actual stress (11.46 MPa) at dp+10, meaning the VW zone center would need ~2000yr of loading to reach failure under the current loading rate.

---

## 3. Root Cause: H8 Was Wrong

### 3.1 What H8 Did

v14 analysis incorrectly concluded that Tandem applies Dirichlet BC only on the bottom boundary. `SetupBoundaryMarkers()` was changed to mark only attr 6 (bottom).

### 3.2 What Tandem Actually Does

Re-analysis of Tandem's source code:

**bp5.toml** (line 11):
```toml
boundary_linear = true
```
This flag applies boundary loading to **all marked boundary surfaces**.

**bp5.geo** (lines 46-56):
```
Physical Surface(1) = {bottom(), top()};     // bottom AND top
Physical Surface(5) = {diri()};               // all four sides
```
Both surface groups (1 and 5) receive boundary loading via `boundary_linear = true`.

**bp5.lua** (lines 27-34):
```lua
function BP5:boundary(x, y, z, t)
    local Vh = self.Vp * t
    if y > 1 then Vh = Vh / 2.0
    elseif y < -1 then Vh = -Vh / 2.0 end
    return Vh, 0, 0   -- u_x = sgn(y)*Vp*t/2
end
```
Applied on ALL non-fault boundaries (top, bottom, and all four sides).

### 3.3 Boundary Attribute Comparison

| Boundary | Attr | Pre-H8 | H8 (v14) | H10 (Tandem-correct) |
|----------|------|--------|----------|---------------------|
| x=-Lx | 1 | Dirichlet | Natural | **Dirichlet** |
| x=+Lx | 2 | Dirichlet | Natural | **Dirichlet** |
| y=+Ly | 3 | Dirichlet | Natural | **Dirichlet** |
| y=-Ly | 4 | Dirichlet | Natural | **Dirichlet** |
| z=0 (top) | 5 | Natural | Natural | **Dirichlet** |
| z=Lz (bottom) | 6 | Dirichlet | Dirichlet | **Dirichlet** |

Note: Pre-H8 code had attrs 1-4,6 as Dirichlet (missing top surface attr 5). Pre-H8 results were poor (749yr recurrence) due to other bugs (H1-H7), not the boundary approach. With H1-H7 now fixed, all-boundary loading should work correctly.

### 3.4 Why Bottom-Only Loading Starves dp+10

With only the bottom boundary driving plate motion, stress must propagate upward through 100km of elastic material. The VW zone center at 10km depth receives drastically reduced loading:

- Loading rate at dp+10: ~6.3 kPa/yr (only 12% of Tandem's ~50 kPa/yr)
- Far-field boundaries (x,y) are traction-free, so no lateral plate motion is imposed
- The fault essentially only sees the distant bottom drive, filtered by elastic compliance

With all boundaries loaded, the displacement field `u_y = sgn(x)*Vp*t/2` is imposed everywhere, creating a uniform far-field plate motion that directly loads the entire fault.

---

## 4. Coordinate Conversion Verification (Tandem <-> MFEM)

Complete verification of coordinate mapping between Tandem and MFEM:

| Aspect | Tandem | MFEM | Status |
|--------|--------|------|--------|
| Along-strike axis | x | y (x2) | Correct |
| Fault-normal axis | y | x (x1) | Correct |
| Depth convention | z < 0 | z > 0 | Handled consistently |
| Loading displacement | u_x = sgn(y)*Vp*t/2 | u_y = sgn(x)*Vp*t/2 | Correct |
| Fault normal | [0,-1,0] | [1,0,0] | Correct (rotated coords) |
| Strike direction | up x n | up x n | Correct |
| Dip direction | s x n | -(s x n) | Correct (negated for z>0 depth) |
| tau_pre ordering | (dip, strike) | (dip, strike) | Correct |
| V_init ordering | (dip, strike) | (dip, strike) | Correct |
| V_nuc value | 0.01 m/s | 0.03 m/s | Intentional (MFEM follows SCEC) |

**Conclusion**: No coordinate conversion bugs found. The `AssembleDirichletLoading()` sign logic (`sgn(centroid(0))`) correctly maps Tandem's `sgn(y)` to MFEM's `sgn(x)`.

---

## 5. Oversized First Event Analysis

### 5.1 The Problem

First event at dp+10 produces 6.37m slip (Tandem: ~3m). The below-Wf region then takes 201yr to "catch up" (6.37m / 31.6 mm/yr), during which the VW zone has a negative stress shadow.

### 5.2 Cause: SCEC vs Tandem Nucleation Parameters

| Parameter | SCEC BP5-QD | Tandem |
|-----------|-------------|--------|
| V_nuc | 0.03 m/s | 0.01 m/s |
| delta_tau | eta * V_nuc | 0 |
| Effect | Stronger nucleation push | Gentler nucleation |

The SCEC-standard parameters create a stronger initial perturbation, leading to a larger first event. Making these configurable allows comparison against Tandem reference data.

---

## 6. Fixes

### H10 — All-Boundary Dirichlet Loading (like Tandem)

**File**: `miniapps/seas/domain/elasticity_operator.hpp`

**Change** — `SetupBoundaryMarkers()` (line 208-227):
```cpp
void SetupBoundaryMarkers()
{
   // BP5: All-boundary Dirichlet loading (matching Tandem bp5.lua)
   //   All attrs 1-6: u = (0, sgn(x)*Vp*t/2, 0)
   //   No free surface -- consistent with Tandem's boundary_linear=true
   int num_bdr = mesh_.bdr_attributes.Size() > 0 ? mesh_.bdr_attributes.Max() : 0;
   dirichlet_bdr_marker_.SetSize(num_bdr);
   dirichlet_bdr_marker_ = 0;

   for (int be = 0; be < mesh_.GetNBE(); be++)
   {
      int attr = mesh_.GetBdrAttribute(be);
      if (attr >= 1 && attr <= num_bdr)
      {
         dirichlet_bdr_marker_[attr - 1] = 1;
      }
   }
}
```

No changes to `AssembleDirichletLoading()` logic -- the existing sign computation from centroid(0) works for all boundary faces. The `dirichlet_bdr_marker_` update propagates automatically to both stiffness and loading assembly.

Update comments in class docstring and `AssembleDirichletLoading()` to reflect all-boundary approach.

### H11 — Configurable Nucleation Parameters

**File**: `miniapps/seas/config/bp5_params.hpp`

Add field:
```cpp
/// Delta-tau multiplier for nucleation zone pre-stress.
/// SCEC BP5-QD: 1.0 (delta_tau = eta * V_nuc)
/// Tandem: 0.0
real_t delta_tau_factor = 1.0;
```

In `tau0_vec()` (line 298-302):
```cpp
// BEFORE:
tau0_scalar += eta_val * Vi_abs;

// AFTER:
tau0_scalar += delta_tau_factor * eta_val * Vi_abs;
```

**File**: `miniapps/seas/tests/verification/bp5_verification_full.cpp`

Add CLI options:
```
--V-nuc <val>              Nucleation slip rate [m/s] (default: 0.03)
--delta-tau-factor <val>   Delta-tau multiplier (default: 1.0, Tandem: 0.0)
```

### VTK Boundary Diagnostic

**File**: `miniapps/seas/tests/verification/bp5_verification_full.cpp`

Add `--dump-bdr-vtk` flag using MFEM's built-in `PrintBdrVTU()`:
```cpp
if (dump_bdr_vtk)
{
   pmesh.PrintBdrVTU(output_dir + "/boundary_attributes");
}
```

### H12 — Update `miniapps/seas/CLAUDE.md`

Current text says "Plate loading is applied at **bottom boundary (z=Lz)**, NOT at far-field x-boundaries". Update to reflect all-boundary Dirichlet approach confirmed from Tandem's actual code.

---

## 7. Files to Modify

| File | Change | Fix |
|------|--------|-----|
| `domain/elasticity_operator.hpp` | `SetupBoundaryMarkers`: all attrs Dirichlet; update comments | H10 |
| `config/bp5_params.hpp` | Add `delta_tau_factor`, use in `tau0_vec()`, add to `Print()` | H11 |
| `tests/verification/bp5_verification_full.cpp` | Add `--V-nuc`, `--delta-tau-factor`, `--dump-bdr-vtk` CLI options | H11, VTK |
| `tests/unit/test_bp5_integration.cpp` | Update if boundary marker changes affect test expectations | H10 |
| `miniapps/seas/CLAUDE.md` | Update BC documentation | H12 |

---

## 8. Implementation Status

| Fix | Issue | Status | Files |
|-----|-------|--------|-------|
| H1 | PsiToTheta per-DOF Dc | Done (v13) | `dieterich_ruina.hpp`, `rate_state_fault.hpp` |
| H2 | Output intervals to SCEC spec | Done (v13) | `bp5_benchmark_output.hpp` |
| H3 | IP tensor-coupled penalty | Done (v13) | `elasticity_operator.hpp` |
| H4 | dt_max 0.5yr -> 0.1yr | Done (v13) | `bp5_verification_full.cpp`, `time_stepper.hpp` |
| H5 | V_nuc mismatch | Info only | -- |
| H6 | Pre-stress formula | Info only | -- |
| H7 | tau_abs floor in solver | Done (v13) | `dieterich_ruina.hpp` |
| H8 | BC: bottom-only Dirichlet | Done (v14) | `elasticity_operator.hpp` |
| H9 | Output 1e-30 floor | Done (v14) | `bp5_benchmark_output.hpp` |
| **H10** | **BC: all-boundary Dirichlet** | **TODO** | `elasticity_operator.hpp` |
| **H11** | **Configurable nucleation** | **TODO** | `bp5_params.hpp`, `bp5_verification_full.cpp` |
| **H12** | **Update CLAUDE.md** | **TODO** | `miniapps/seas/CLAUDE.md` |

---

## 9. Expected Impact

With all-boundary Dirichlet (matching Tandem) + configurable nucleation:

| Metric | Current (bottom-only) | Expected (all-BC) | Tandem |
|--------|----------------------|-------------------|--------|
| dp+10 loading rate | 6.3 kPa/yr | ~40-50 kPa/yr | ~50 kPa/yr |
| dp+10 total slip (1800yr) | 7m | ~40-55m | ~55m |
| dp+10 event V_peak | 10^(-3) | ~0.1-1 m/s | ~1 m/s |
| Recurrence (surface) | ~250yr | ~240yr | ~240yr |
| Recurrence (depth) | >1000yr | ~240yr | ~240yr |

---

## 10. Verification Plan

### Step 1: Build and run unit tests
```bash
conda activate mfem-dev
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas
make -j test-bp5-integration test-bp5-output
./test-bp5-integration && ./test-bp5-output
```

### Step 2: Verify boundary attributes visually
```bash
make -j seas_bp5_full
mpirun -np 4 tests/verification/seas_bp5_full --inline-mesh --nx 4 --ny 4 --nz 2 --tfinal 0 --dump-bdr-vtk
# Open output/boundary_attributes.vtu in ParaView -> color by Attribute
# Verify all 6 boundary faces are present with distinct attribute numbers
```

### Step 3: Quick smoke test
```bash
mpirun -np 4 tests/verification/seas_bp5_full --inline-mesh --nx 4 --ny 4 --nz 2 --tfinal 31557600
```

### Step 4: Production run (SCEC ICs, all-boundary BC)
```bash
ibrun ./seas_bp5_full --mesh bp5/mesh/bp5_1000m.msh \
   --ref-dir bp5/benchmark_data --output-dir bp5/results_1000m_allbc
```

### Step 5: Production run (Tandem ICs + all-boundary BC)
```bash
ibrun ./seas_bp5_full --mesh bp5/mesh/bp5_1000m.msh \
   --V-nuc 0.01 --delta-tau-factor 0.0 \
   --ref-dir bp5/benchmark_data --output-dir bp5/results_1000m_allbc_tandem_ic
```

### Checks:
1. dp+10 participates in events with V_peak > 0.01 m/s
2. dp+10 total slip over 1800yr > 30m
3. Surface recurrence ~240yr maintained
4. No displacement blowup or NaN
