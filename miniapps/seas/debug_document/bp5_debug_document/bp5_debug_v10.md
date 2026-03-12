# BP5 Debug v10: Dip Sign Convention, Mesh Refinement, and DG Traction Noise

## Problem

After the v9 fix (Dirichlet BCs on attrs 3,4,6), the BP5 simulation shows:

1. **Slip_dip sign reversal**: MFEM slip_dip grows positive (+0.03m at 190yr), Tandem is negative (-0.15m at 1800yr). Same physical motion, opposite sign convention.
2. **V_dip spikes**: Single-timestep velocity spikes during seismic events at 6/10 stations.
3. **Sharp mesh transition**: 1km (fault) → 40km (far-field) with no intermediate grading, causing DG traction noise.

### Evidence from New Run Plots (with v9 BC fix)

| Station | slip_dip (MFEM) | slip_dip (Tandem) | Sign match? |
|---------|-----------------|-------------------|-------------|
| (0, 0)  | +0.03m (190yr)  | -0.15m (1800yr)   | REVERSED    |
| (16, 0) | +0.15m (190yr)  | negative           | REVERSED    |
| (-24, 10) | +0.0008m       | negative           | REVERSED    |

V_dip shows a large spike at the first seismic event (~10^-3 m/s), then settles to ~10^-12 m/s.

## Root Cause 1: Dip Direction Sign Convention Mismatch

### MFEM Fault Basis (fault_basis.hpp:112-136)

```
MFEM coordinates: x = fault-normal, y = along-strike, z = depth (positive down)
ref_normal = (1, 0, 0)
up = (0, 0, 1)  ← points downward (into the earth) in physical space

strike = up × n = (0,0,1) × (1,0,0) = (0, 1, 0)   → along-strike (+y)
dip = strike × n = (0,1,0) × (1,0,0) = (0, 0, -1)  → UPWARD (toward surface)
```

### Tandem Fault Basis (Curvilinear.cpp:281-291)

```
Tandem coordinates: X = along-strike, Y = fault-normal, Z = elevation (positive UP)
ref_normal = (0, -1, 0)
up = (0, 0, 1)  ← points upward in physical space

strike = up × n = (0,0,1) × (0,-1,0) = (1, 0, 0)   → along-strike
dip = strike × n = (1,0,0) × (0,-1,0) = (0, 0, -1)  → DOWNWARD (into earth)
```

### Coordinate Mapping

From `bp5.geo`: MFEM_x = Tandem_Y, MFEM_y = Tandem_X, MFEM_z = -Tandem_Z

Converting Tandem dip (0,0,-1)_Tandem to MFEM coords:
- MFEM_x = Tandem_Y component = 0
- MFEM_y = Tandem_X component = 0
- MFEM_z = -Tandem_Z component = -(-1) = **+1**

**Tandem dip in MFEM coords = (0, 0, +1) → physically DOWNWARD (into the earth)**
**MFEM dip = (0, 0, -1) → physically UPWARD (toward surface)**

The dip basis vectors are exactly **opposite in physical direction**.

### Impact

For the same physical motion (upward dip slip at station (0,0)):
- Tandem: slip_dip = **negative** (slip in -dip direction, where dip=down → -dip=up)
- MFEM: slip_dip = **positive** (slip in +dip direction, where dip=up)

All dip quantities (slip_dip, V_dip, tau_dip) have opposite sign convention.

## Root Cause 2: DG-BR2 Traction Noise (Sharp Mesh Transition)

### Current Mesh Sizing (bp5.geo:97-99)

```gmsh
MeshSize{ PointsOf{Volume{:};} } = res;         // 40 km far-field
MeshSize{ PointsOf{Surface{fault_surfs()};} } = res_f;  // 1 km on fault
// NO sizing on nuc1/nuc2/nuc3 refinement zones
// NO volume mesh grading field
```

The element size jumps from 1 km on the fault surface to 40 km in the volume — a **40:1 ratio**. Gmsh's automatic gradation helps somewhat, but the transition is still very abrupt compared to Tandem's structured hex mesh.

### Why This Causes V_dip Spikes

During seismic events, `V_abs` reaches ~1-5 m/s. The DG-BR2 traction computation has numerical noise proportional to element size mismatch. For pure strike-slip loading, `tau_dip` should be exactly zero, but DG noise produces O(0.01-0.1 MPa) in the dip direction.

The friction solver computes:
```
V_dip = (V_abs / tau_abs) * tau_dip
```

With V_abs ~ 1 m/s, tau_abs ~ 20 MPa, tau_dip noise ~ 0.01 MPa:
**V_dip ~ 5e-4 m/s** — appearing as a spike from 10^-12 to 10^-3 m/s.

### Why Tandem Doesn't Need a Guard

Tandem uses **SBP-SAT on structured hexahedral meshes**, which produces inherently cleaner off-diagonal tractions. Tandem's friction solver uses the same unguarded formula:
```cpp
return -(V / tauAbs) * tauAbsVec;  // DieterichRuinaAgeing.h:119
```
No directional guard, no traction filter. The noise level is simply lower.

## Fix 1: Dip Sign Convention

**File**: `miniapps/seas/fault/fault_basis.hpp` (lines 125-136)

Negate the dip vector after computing `d_vec = strike × n`:

```cpp
// BEFORE:
d_vec[0] = s[1] * n_raw(2) - s[2] * n_raw(1);
d_vec[1] = s[2] * n_raw(0) - s[0] * n_raw(2);
d_vec[2] = s[0] * n_raw(1) - s[1] * n_raw(0);

// AFTER:
// Negate dip to match Tandem/SCEC convention (positive dip = physically downward).
// In MFEM coords (z = depth positive down), strike × n gives (0,0,-1) = upward.
// Tandem's dip convention points downward (into the earth), so we negate.
d_vec[0] = -(s[1] * n_raw(2) - s[2] * n_raw(1));
d_vec[1] = -(s[2] * n_raw(0) - s[0] * n_raw(2));
d_vec[2] = -(s[0] * n_raw(1) - s[1] * n_raw(0));
```

### Self-Consistency Check

tangent1 (dip) is used in two places:

1. **ProjectTraction** (`fault_basis.hpp:271-291`):
   `tau_dip = T_global · tangent1` → sign of tau_dip flips ✓

2. **EmbedSlip** (`fault_basis.hpp:194-215`):
   `delta_u += slip_dip * tangent1` → sign of embedded displacement flips ✓

The friction solver sees `tau_vec[0] = tau_pre_dip + traction_dip`. Both `tau_pre_dip` and `traction_dip` change sign because they're both projected through tangent1. The scalar V_abs is unchanged (depends on |tau|). V_dip changes sign. slip_dip changes sign.

**All dip quantities are self-consistently negated. Physics is identical.**

### Impact on tau_pre_dip

In `bp5_params.hpp:tau0_vec()`, the pre-stress is constructed as `tau = tau0_scalar * V_init / |V_init|`. Since V_init has zero dip component for all stations (BP5 is pure strike-slip), tau_pre_dip = 0 regardless of the dip sign convention. No change needed in bp5_params.hpp.

## Fix 2: Mesh Refinement Around Fault Zone

**File**: `miniapps/seas/bp5/mesh/bp5.geo` (lines 97-99)

### Strategy

Use Gmsh's Distance + Threshold fields to create smooth element size grading from the fault surface into the volume. Also add explicit intermediate sizing on the nucleation refinement zones.

### Changes

Replace the current mesh sizing section with:

```gmsh
// --- Mesh sizing ---
// Far-field volume: coarse
MeshSize{ PointsOf{Volume{:};} } = res;

// Fault surface: benchmark resolution
MeshSize{ PointsOf{Surface{fault_surfs()};} } = res_f;

// Smooth mesh grading from fault into volume.
// Prevents sharp element size jumps that cause DG traction noise.
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfs()};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = res_f;      // At fault: res_f (1 km for benchmark)
Field[2].SizeMax = res;         // Far from fault: res (40 km)
Field[2].DistMin = 0;           // Start grading at fault surface
Field[2].DistMax = 80;          // Reach far-field size at 80 km distance

Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;
```

### Expected Effect

| Distance from fault | Element size |
|--------------------|-------------|
| 0 km (on fault)    | 1 km        |
| 5 km               | ~3.4 km     |
| 10 km              | ~5.9 km     |
| 20 km              | ~10.8 km    |
| 40 km              | ~20.5 km    |
| 80+ km             | 40 km       |

Adjacent element size ratio is controlled by Gmsh's mesh gradation algorithm (typically ≤ 1.3-1.5), ensuring no sharp jumps.

### Impact on Element Count

The smooth grading will significantly increase the total element count compared to the current mesh. Rough estimate:
- Current (sharp transition): ~100K-200K tets
- With grading: ~500K-1M tets

This increases memory and solve time but should dramatically improve DG traction quality near the fault.

## Fix 3 (CONDITIONAL): Traction Noise Filter

**Only apply if mesh refinement (Fix 2) alone does not resolve V_dip spikes.**

**File**: `miniapps/seas/fault/rate_state_fault.hpp` (lines 415-416)

After assembling `tau_vec`, suppress traction components below a relative noise floor:

```cpp
real_t tau_vec[2] = {tau_pre_(2*i) + traction(2*i),
                     tau_pre_(2*i+1) + traction(2*i+1)};

// DG traction noise floor: suppress components below relative threshold.
// DG-BR2 on unstructured tets produces O(h) off-diagonal traction noise
// that Tandem's SBP-SAT on structured hex meshes does not.
real_t tau_abs = std::sqrt(tau_vec[0]*tau_vec[0] + tau_vec[1]*tau_vec[1]);
if (tau_abs > 0.0)
{
   constexpr real_t r_thresh = 1e-4;  // relative noise floor
   for (int c = 0; c < 2; ++c)
   {
      if (std::abs(tau_vec[c]) < r_thresh * tau_abs)
      {
         tau_vec[c] = 0.0;
      }
   }
}
```

**Threshold choice**: `r_thresh = 1e-4`. For BP5 with tau_strike ~ 20 MPa, this suppresses dip traction below 2 kPa — well below any physically meaningful signal.

**Why filter traction, not velocity**: Acts at the source (noisy DG traction), not the symptom (amplified velocity). Keeps the friction solver (`SolveSlipRateVectorPsi`) mathematically identical to Tandem's implementation.

## What Does NOT Change

| Component | Reason |
|-----------|--------|
| Friction solver (`dieterich_ruina.hpp`) | Mathematically correct, identical to Tandem |
| Pre-stress formula (`bp5_params.hpp`) | Follows SCEC BP5 spec (V_nuc=0.03, delta_tau=eta*Vi) |
| Stiffness matrix assembly | dirichlet_bdr_marker_ already includes attrs 1,2,3,4,6 (v9 fix) |
| Dirichlet loading (`AssembleDirichletLoading`) | Centroid-based sign already correct (v9 fix) |
| RHS blowup check | Already NaN-only (v9 fix) |
| Slip assembly / traction computation | Sign convention unchanged; tangent1 negation is self-consistent |

## Files Modified

| File | Change |
|------|--------|
| `fault/fault_basis.hpp:125-136` | Negate dip direction to match Tandem/SCEC convention |
| `bp5/mesh/bp5.geo:97-99` | Add Distance/Threshold fields for smooth mesh grading |
| `fault/rate_state_fault.hpp:415-416` | (conditional) Traction noise filter |
| `tests/unit/test_fault_basis.cpp` | Update expected dip direction |

## Verification

1. **Build**: `conda activate mfem-dev && make seas_bp5_full`
2. **Unit tests**: Run `test_fault_basis` — verify dip = (0,0,+1) in MFEM coords
3. **Mesh generation**: `conda activate pythonenv && gmsh -3 bp5.geo -setnumber res_f 1 -o bp5_1000m.msh`
   - Verify gradual element size transition (no sharp jumps near fault)
   - Check total element count
4. **Local smoke test**: Run small BP5 mesh — no crashes, correct dip sign at t=0
5. **Cluster run (1000m)**: Submit with new mesh and check:
   - slip_dip sign matches Tandem (both negative for physical upward motion)
   - V_dip spikes reduced or eliminated
   - tau_dip oscillation amplitude reduced
   - Strike components unchanged (recurrence ~200yr, amplitudes match)
   - Nucleation zone participates in events

## Commit Strategy

1. **Commit A**: Dip sign convention (`fault_basis.hpp` + `test_fault_basis.cpp`)
2. **Commit B**: Mesh refinement (`bp5.geo`) — requires regenerating mesh
3. **Commit C** (if needed after B results): Traction noise filter (`rate_state_fault.hpp`)
