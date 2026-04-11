# BP5 Debug v9: Missing Dirichlet BCs and RHS BLOWUP False Alarm

## Problem

After the v8 sign convention fix, BP5 successfully captures earthquake cycles. However, comparing MFEM results against Tandem over 1800 years reveals progressive divergence:

### Observations from Station Plots (1000m mesh, CG+AMG run)

1. **Spurious dip-slip** grows monotonically at all stations:
   - strk+00dp+10: dip slip → -0.08 m (should be ~0)
   - strk+16dp+00: dip slip → -0.4 m
   - strk+36dp+00: dip slip → -0.2 m
   - The dip slip grows steadily between earthquakes (not just coseismically)

2. **Dip-direction shear stress** oscillates with growing amplitude:
   - ±0.3–0.5 MPa at center stations (Tandem shows ~0)
   - This drives V_dip ~ 1e-11 m/s during interseismic (should be ~0)

3. **Earthquake timing drifts** — MFEM cycles get progressively earlier than Tandem:
   - First 2-3 cycles match well
   - By 8th cycle, MFEM is ~50 years ahead

4. **Nucleation zone (strk-24dp+10) barely participates**:
   - Initial V_strike = 0.03 m/s (correct nucleation velocity)
   - But subsequent earthquake rupture barely reaches this station
   - V_strike drops to ~1e-14 between events (vs Tandem's ~1e-9)

5. **Late-time instability**: V_max eventually reaches ~3748 m/s, suggesting accumulated error leads to blowup

### Observations from Frontera MUMPS Run (Job 7594911, 4 nodes/224 cores)

1. **MUMPS INFOG(1)=-9** at startup: Insufficient memory for factorization on some ranks
   - MFEM auto-retries with ICNTL(14) += 20 (benign, simulation continues)

2. **69,018 "RHS BLOWUP" false-alarm messages** starting at ~6.4 years:
   ```
   [Rank 208] RHS BLOWUP: before_slip=0 after_slip=0 final=1.00002e+15 slip_max=0 time=2.01487e+08
   ```
   - Rank 208 has **no fault DOFs** (slip_max=0, before_slip=0, after_slip=0)
   - The entire `final` value comes from **Dirichlet loading**: u_D = Vp·t/2 grows linearly with time
   - At t ≈ 6.4 yr (2e8 s), the RHS norm crosses the hardcoded 1e15 threshold
   - By t ≈ 300 yr (9.5e9 s), `final` reaches 4.7e16
   - **This is expected behavior, not a blowup** — the RHS naturally scales with displacement

3. **Only 300 years completed in 48 hours** (should be 1800 yr):
   - I/O bottleneck from 69K warning messages (~10 messages per RK stage)
   - Only 1 earthquake detected (the initial nucleation at t=0)
   - Interseismic V_max = 8e-11 m/s at t=300 yr (reasonable)

## Root Cause: Missing Dirichlet Boundary Conditions

### SCEC BP5-QD Specification (Jiang et al., 2022)

All non-free-surface external boundaries require:

**u(x, t) = (Vp/2) · t · sgn(x₁) · ê₂**

where x₁ is the fault-normal coordinate and ê₂ is the along-strike direction.

### Current MFEM vs. Required BCs

| Boundary | Attr | SCEC Spec | MFEM Current | Status |
|----------|------|-----------|-------------- |--------|
| x = -Lx (fault-normal, -100 km) | 1 | Dirichlet: u_y = -Vp·t/2 | Dirichlet: u_y = -Vp·t/2 | ✅ |
| x = +Lx (fault-normal, +100 km) | 2 | Dirichlet: u_y = +Vp·t/2 | Dirichlet: u_y = +Vp·t/2 | ✅ |
| y = +Ly (along-strike, +200 km) | 3 | Dirichlet: u_y = sgn(x)·Vp·t/2 | **Natural BC (zero traction)** | ❌ |
| y = -Ly (along-strike, -200 km) | 4 | Dirichlet: u_y = sgn(x)·Vp·t/2 | **Natural BC (zero traction)** | ❌ |
| z = 0 (free surface) | 5 | Zero traction | Zero traction | ✅ |
| z = Lz (bottom, 100 km depth) | 6 | Dirichlet: u_y = sgn(x)·Vp·t/2 | **Natural BC (zero traction)** | ❌ |

**Three of five non-free-surface boundaries are missing Dirichlet conditions.**

### Evidence

1. **Gmsh geo file** (`bp5/mesh/bp5.geo`) comments correctly label attrs 3,4 as "Dirichlet" — the intent was always to have them, but the code at `elasticity_operator.hpp:221` only marks attrs 1,2.

2. **Tandem comparison**: Tandem's `BP5:boundary(x, y, z, t)` applies Dirichlet BCs at ALL external boundaries with position-dependent sign: `u_x = ±Vp·t/2` depending on `y` (Tandem's fault-normal coordinate).

3. **Code location**: `SetupBoundaryMarkers()` at `elasticity_operator.hpp:218-225`:
   ```cpp
   // Only marks attrs 1 and 2!
   if (attr == 1 || attr == 2) { dirichlet_bdr_marker_[attr - 1] = 1; }
   ```

### Why Missing BCs Cause Spurious Dip-Slip

1. **Free bottom boundary (z=100 km)**: Without Dirichlet constraint, the bottom boundary deforms freely. Strike-slip loading creates vertical displacement relaxation that propagates back to the fault as dip-direction stress.

2. **Free along-strike boundaries (y=±200 km)**: Without constraint, the material at the along-strike edges relaxes freely, creating edge effects that affect stress transfer rates.

3. **Accumulation**: Each earthquake cycle adds slip that interacts with the unconstrained boundaries. The boundary effects are small per cycle but compound over 8 cycles → progressive divergence.

4. **Dip stress magnitude**: The observed ~0.3 MPa dip traction is ~1.5% of the strike traction (~20 MPa). This is consistent with the free boundary at z=100 km (60 km below fault tip) creating O(1%) stress leakage.

### Why BCs at attrs 3, 4, 6 need sgn(x)

The mesh is split into two volumes at x=0 (fault plane). On the y=±Ly and z=Lz boundaries, mesh faces exist on BOTH sides of x=0:
- Faces with centroid x > 0 → u_y = +Vp·t/2 (positive plate side)
- Faces with centroid x < 0 → u_y = -Vp·t/2 (negative plate side)

This is straightforward to implement since each face is entirely on one side of x=0 (guaranteed by the two-volume mesh split).

## Fix

### Change 1: Mark Additional Dirichlet Boundaries

**File**: `domain/elasticity_operator.hpp`, `SetupBoundaryMarkers()` (~line 221)

```cpp
// BEFORE:
if (attr == 1 || attr == 2) { dirichlet_bdr_marker_[attr - 1] = 1; }

// AFTER:
if (attr == 1 || attr == 2 || attr == 3 || attr == 4 || attr == 6)
{
   dirichlet_bdr_marker_[attr - 1] = 1;
}
```

This automatically includes the new boundaries in:
- The **stiffness matrix** (boundary face integrators use `dirichlet_bdr_marker_`)
- The **Dirichlet loading RHS** (AssembleDirichletLoading loops over marked boundaries)

### Change 2: Position-Based Dirichlet Sign

**File**: `domain/elasticity_operator.hpp`, `AssembleDirichletLoading()` (~line 1351-1356)

Replace attribute-based sign logic with position-based sgn(x) logic:

```cpp
// BEFORE:
if (attr != 1 && attr != 2) { continue; }
real_t u_D[3] = {0.0, 0.0, 0.0};
if (attr == 2) { u_D[1] = Vp_ * time / 2.0; }   // +x1
else           { u_D[1] = -Vp_ * time / 2.0; }   // -x1

// AFTER: Use face centroid x-coordinate to determine sign
if (dirichlet_bdr_marker_[attr - 1] != 1) { continue; }

// Compute face centroid to determine which side of fault
Vector centroid(3);
centroid = 0.0;
{
   ElementTransformation *eltransf = mesh_.GetBdrElementTransformation(be);
   const IntegrationRule &ir_c = IntRules.Get(eltransf->GetGeometryType(), 1);
   for (int p = 0; p < ir_c.GetNPoints(); p++)
   {
      eltransf->SetIntPoint(&ir_c.IntPoint(p));
      Vector phys(3);
      eltransf->Transform(ir_c.IntPoint(p), phys);
      centroid.Add(1.0 / ir_c.GetNPoints(), phys);
   }
}
real_t sign = (centroid(0) > 0.0) ? 1.0 : -1.0;

real_t u_D[3] = {0.0, 0.0, 0.0};
u_D[1] = sign * Vp_ * time / 2.0;
```

This correctly handles ALL Dirichlet boundaries:
- x=+Lx: centroid.x = +100 km > 0 → sign = +1 → u_y = +Vp·t/2 ✓
- x=-Lx: centroid.x = -100 km < 0 → sign = -1 → u_y = -Vp·t/2 ✓
- y=±Ly: centroid.x varies per face → position-dependent ✓
- z=Lz: centroid.x varies per face → position-dependent ✓

### Change 3: Fix RHS BLOWUP False Alarm

**File**: `domain/elasticity_operator.hpp` (~line 1717)

The Dirichlet RHS grows linearly: ||RHS|| ∝ μ · |u_D| / h ∝ μ · Vp · t / h. The hardcoded 1e15 threshold is crossed at ~6.4 years and triggers 69K false alarms.

Options (pick one):
- **(A) Remove the non-NaN check entirely** — the diagnostic has served its purpose
- **(B) Scale threshold with time**: `max(1e15, C * |time * Vp| * mu / h)`
- **(C) Rate-limit the output**: only print once per 1000 triggers

Recommended: **(A)** — simplest. Keep only the NaN check:

```cpp
// BEFORE:
if (rhs_final > 1e15 || std::isnan(rhs_final))

// AFTER:
if (std::isnan(rhs_final))
```

### Change 4: Update Comments

- `bp5.geo` lines 13-14: Fix attrs 1,2 from "natural BC" to "Dirichlet"
- `elasticity_operator.hpp` line 45: Update to mention all Dirichlet boundaries
- `elasticity_operator.hpp` lines 1334-1342: Update Dirichlet loading comments

## What Does NOT Change

| Component | Reason |
|-----------|--------|
| Stiffness matrix code | `dirichlet_bdr_marker_` already passed to `AddBdrFaceIntegrator` — marking new attrs automatically includes them |
| Friction solver | Sign convention unchanged from v8 |
| Slip assembly | `sign * EmbedSlip` convention unchanged |
| Traction computation | Same `sign` convention |
| Below-fault rate | Already correct: (0, Vp, 0) |
| MUMPS configuration | Auto-retry handles INFOG(1)=-9 |

## Impact of the Stiffness Matrix Change

Adding attrs 3, 4, 6 to the Dirichlet marker adds boundary face integrals to the stiffness matrix. Since the matrix is assembled once and cached:
- **One-time cost**: Slightly larger matrix (more boundary face contributions)
- **MUMPS**: Single re-factorization (may need slightly more memory)
- **CG+AMG**: Setup unchanged, iteration count may decrease (better-conditioned problem)

## Verification

1. **Build**: `conda activate mfem-dev && make seas_bp5_full`
2. **Unit tests**: All BP5 tests should still pass (friction, params, fault operator, integration, fault basis)
3. **Local smoke test**: Run small BP5 mesh — verify no BLOWUP spam, MUMPS succeeds
4. **Frontera** (1000m mesh): Expected improvements:
   - No BLOWUP message flood (I/O bottleneck removed)
   - Dip slip should be negligible (~0)
   - Earthquake timing should match Tandem more closely
   - Nucleation zone should participate in earthquake cycles
   - Simulation should reach 1800 years within wall time
