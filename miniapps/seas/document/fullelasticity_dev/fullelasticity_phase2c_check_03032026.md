# Phase 2c Implementation Review: FaultGeometry 3D Extension

**Date**: 2026-03-03
**Last updated**: 2026-03-03 (revision 2 — all bugs resolved, all test gaps covered)
**Reviewing**: `fullelasticity_phase2c_plan_03022026.md` vs. current implementation
**Files checked**: `fault/fault_geometry.hpp`, `config/bp5_params.hpp`,
`tests/unit/test_elasticity_operator.cpp`, `tests/unit/test_bp5_params.cpp`,
`domain/elasticity_operator.hpp` (GetFaultCoords2D)
**References**: Tandem `examples/tandem/3d/bp5.lua`, SCEC BP5-QD spec

---

## 1. Overall Assessment

Phase 2c is functionally complete — the BP5 FaultGeometry constructor,
ComputeBP5Params(), new accessors, and backward compatibility with BP2 are
all implemented. BP5Params has excellent independent test coverage (12 test
functions in test\_bp5\_params.cpp).

**Revision 2 status:**
- **0 open bugs** (all 4 from revision 1 are COMPLETED)
- **0 missing tests** (all 6 test gaps now covered by TestFaultGeometry3DValues)
- **3 INFO items** (dedup, asinh, nucleation tolerance)
- BP5Params unit tests are thorough and self-consistent
- Implementation is ready for integration testing with the SEAS time-stepper

**Change summary since revision 1:**
- §2.1 (Print uses BP2 b): COMPLETED — now uses `is_bp5_ ? bp5_params_.b : params_.b`
- §2.2 (IsVelocityWeakening uses BP2): COMPLETED — signature changed to `int dof_idx`,
  uses precomputed `a_values_` and correct `b` value
- §2.3 (GetVWDepth/GetVSDepth): COMPLETED — `MFEM_VERIFY(!is_bp5_, ...)` guards added
- §2.4 (V\_init\_vec comment labels): COMPLETED — comments now correctly say
  V[0]=along-dip, V[1]=along-strike
- TestFaultGeometry3DValues added (lines 358-431): 13 new assertions covering all 6 gaps

---

## 2. Bugs / Potential Bugs

### 2.1 [COMPLETED] Print() uses BP2 params\_.b for BP5 FaultGeometry

Print() now uses the correct `b` value (line 393):
```cpp
real_t b_val = is_bp5_ ? bp5_params_.b : params_.b;
```
Verified at line 397: `if (a_values_(i) < b_val) { vw_count++; }`. ✓
New TestFaultGeometry3DValues exercises Print() for BP5 path (line 426). ✓

### 2.2 [COMPLETED] IsVelocityWeakening/GetVWDOFs/GetVSDOFs use BP2 params\_

Signature changed from `IsVelocityWeakening(real_t z)` to
`IsVelocityWeakening(int dof_idx)` (line 316). Now uses precomputed
`a_values_` (correct for both BP2 and BP5) and the appropriate `b`:
```cpp
bool IsVelocityWeakening(int dof_idx) const
{
   real_t b_val = is_bp5_ ? bp5_params_.b : params_.b;
   return a_values_(dof_idx) < b_val;
}
```
GetVWDOFs (line 350) and GetVSDOFs (line 363) now call the fixed method. ✓

### 2.3 [COMPLETED] GetVWDepth/GetVSDepth use BP2 params\_ for BP5

Both methods now have `MFEM_VERIFY(!is_bp5_, ...)` guards (lines 328-330, 339-341):
```cpp
real_t GetVWDepth() const
{
   MFEM_VERIFY(!is_bp5_, "GetVWDepth() not applicable for BP5 (2D VW zone)");
   return -params_.H;
}
```
Calling these for a BP5 geometry will produce a clear error message. ✓

### 2.4 [COMPLETED] V\_init\_vec/tau0\_vec comment labels were SWAPPED

Comments in `bp5_params.hpp` now correctly read (line 236-237):
```cpp
///   V[0] = along-dip ≈ 0 (V_zero placeholder)
///   V[1] = along-strike = dominant slip component
```
This matches Tandem's convention (Vi1=dip≈0, Vi2=strike=dominant) and the
accessor docstring layout `[V_dip_0, V_strike_0, ...]`. ✓

---

## 3. Deviations from Plan

### 3.1 Plan lists tau\_pre/V\_init\_vec layout as (dip, strike)

**Plan**: `[τ_dip_0, τ_strike_0, ...]` and `[V_dip_0, V_strike_0, ...]`
**Code**: Stores `Vi[0]` at `2*i` and `Vi[1]` at `2*i+1`

Since Vi[0] = dip and Vi[1] = strike (see §2.4), the actual storage is
(dip, strike) per DOF. Plan and code match. ✓

### 3.2 BP5 constructor parallel handling

**Plan**: "Each rank computes parameters for its local fault DOFs only"
**Code**: BP5 constructor properly handles parallel case with
`num_global_fault_dofs_` computed via `GlobalSumInt()`. ✓

### 3.3 Missing early return for zero global DOFs

The BP2 constructor checks `num_global_fault_dofs_ == 0` and warns:
```cpp
if (num_global_fault_dofs_ == 0)
{
   MFEM_WARNING("FaultGeometry: No fault DOFs found on any rank");
   return;
}
```

The BP5 constructor skips this check, going directly to:
```cpp
if (num_fault_dofs_ == 0) { return; }
```

In serial, `num_fault_dofs_ == 0` implies `num_global_fault_dofs_ == 0`,
so this is equivalent. In parallel, a rank with local DOFs=0 returns early
regardless. Not a bug, but the warning message is lost for BP5.

---

## 4. Correctness Verification Against Tandem

### 4.1 a(x2, x3): CORRECT ✓

**Tandem** (`bp5.lua:77-87`):
```lua
local d = -z
local s = math.abs(x)
-- VW core: h_s+h_t <= d <= h_s+h_t+H and s <= l/2
-- VS zone: d <= h_s or d >= h_s+2*h_t+H or s >= l/2+h_t
-- Transition: r = max(|d-h_s-h_t-H/2|-H/2, s-l/2)/h_t
```

**MFEM** (`bp5_params.hpp:162-183`): Uses `x3` directly (positive downward,
equivalent to `d = -z` in Tandem). Formula, zone boundaries, and transition
`r` computation match exactly. Tested in test\_bp5\_params.cpp with ~15
assertions (VW core, VS zones, transition, symmetry, edge cases). ✓

### 4.2 IsNucleationZone: CORRECT ✓

**Tandem** (`bp5.lua:49-57`):
```lua
local d = -z; local s = x  -- Note: signed x, not |x|
-- h_s+h_t <= d+eps and d-eps <= h_s+h_t+H and -l/2 <= s+eps and s-eps <= -l/2+w
```

**MFEM** (`bp5_params.hpp:192-197`): Same zone bounds. No `eps` tolerance
(equivalent to Tandem's `bp5_exact = BP5.new({eps=0.0})`). See §5.1 for
tolerance discussion. Tested with 9 assertions in test\_bp5\_params.cpp. ✓

### 4.3 L\_of\_x2\_x3 / Dc\_of\_x2\_x3: CORRECT ✓

**Tandem**: L=0.13 in nucleation, 0.14 elsewhere.
**MFEM**: `L_nuc=0.13`, `L0=0.14`. `Dc_of_x2_x3` is alias for `L_of_x2_x3`. ✓

### 4.4 V\_init\_vec: CORRECT ✓

**Tandem** (`bp5.lua:70-74`):
```lua
return self.Vzero, 0.01   -- nucleation
return self.Vzero, self.Vp  -- elsewhere
```

**MFEM** (`bp5_params.hpp:247-258`):
```cpp
V[0] = V_zero;
V[1] = IsNucleationZone(...) ? V_nuc : V_init;
```

Values match: V\_zero=1e-20, V\_nuc=0.01, V\_init=1e-9 = Vp. ✓
Component labels now correct (§2.4 fixed). ✓

### 4.5 tau0\_vec: CORRECT (with INFO-level difference) ✓

**Tandem** (`bp5.lua:94-103`):
```lua
local tau0 = sn * ax * math.asinh((Vi2 / (2.0 * self.V0)) * e) + eta * Vi2
return -tau0 * Vi1 / Vi, -tau0 * Vi2 / Vi
```

**MFEM** (`bp5_params.hpp:274-296`):
```cpp
real_t tau0_scalar = sigma_n * a * asinh((Vi_abs / (2*V0)) * e) + eta * Vi_abs;
tau[0] = -tau0_scalar * Vi[0] / Vi_abs;
tau[1] = -tau0_scalar * Vi[1] / Vi_abs;
```

Tandem uses `Vi2` (second component only) in the asinh and damping terms;
MFEM uses `Vi_abs` (full magnitude). Since `Vi[0] = V_zero = 1e-20` and
`Vi[1] = O(1e-9)` to `O(1e-2)`, the relative difference is ~1e-18, which
is below machine epsilon. Direction computation is identical. ✓

Self-consistency verified in test\_bp5\_params.cpp TestTau0Vec: friction
solver recovers |V\_init| from tau0 both outside and inside the nucleation
zone, within 1e-6 relative tolerance. ✓

### 4.6 eta: CORRECT ✓

**Tandem**: `eta = cs * rho / 2` = 3464 × 2670 / 2 = 4,624,440 Pa·s/m
**MFEM**: `eta = mu / (2*cs)` = `rho*cs^2 / (2*cs)` = `rho*cs/2` = same. ✓

### 4.7 Material parameters: CORRECT ✓

All constants match (verified in test\_bp5\_params.cpp TestMaterialProperties):

| Parameter | Tandem | MFEM | Match |
|-----------|--------|------|-------|
| rho | 2.670 g/cm³ | 2670.0 kg/m³ | ✓ |
| cs | 3.464 km/s | 3464.0 m/s | ✓ |
| nu | 0.25 | 0.25 | ✓ |
| a0 | 0.004 | 0.004 | ✓ |
| amax | 0.04 | 0.04 | ✓ |
| b | 0.03 | 0.03 | ✓ |
| f0 | 0.6 | 0.6 | ✓ |
| V0 | 1e-6 | 1e-6 | ✓ |
| Vp | 1e-9 | 1e-9 | ✓ |
| V\_zero | 1e-20 | 1e-20 | ✓ |
| V\_nuc | 0.01 | 0.01 | ✓ |
| sigma\_n | 25 MPa | 25e6 Pa | ✓ |
| h\_s | 2 km | 2e3 m | ✓ |
| h\_t | 2 km | 2e3 m | ✓ |
| H | 12 km | 12e3 m | ✓ |
| l | 60 km | 60e3 m (l\_vw) | ✓ |
| w | 12 km | 12e3 m (w\_nuc) | ✓ |
| L0 | 0.14 m | 0.14 m | ✓ |
| L\_nuc | 0.13 m | 0.13 m | ✓ |

### 4.8 ComputeBP5Params: CORRECT ✓

**File**: `fault_geometry.hpp` lines 565-593

Correctly iterates over all local fault DOFs, reads (x2, x3) coordinates,
calls BP5Params parameter functions, and stores results in correctly-sized
vectors. Eta is computed once outside the loop (constant). ✓

### 4.9 GetFaultCoords2D: CORRECT ✓

**File**: `elasticity_operator.hpp` lines 1135-1166

Returns mesh y-coordinate as x2 (along-strike) and z-coordinate as x3
(depth) from fault face centroids. Consistent with the SCEC convention
where mesh (x, y, z) maps to (x1, x2, x3). ✓

---

## 5. Design / Info Items

### 5.1 [INFO] IsNucleationZone uses exact boundaries (no eps tolerance)

Tandem uses a configurable `eps` tolerance in `in_nucleation`:
```lua
self.h_s + self.h_t <= d+eps and d-eps <= self.h_s + self.h_t + self.H
```

with `eps=0.0` for exact, `±1e-3` for expanded/contracted zones (in km).

MFEM uses exact comparison (equivalent to `eps=0.0`). Since fault DOF
coordinates come from face centroids computed with double precision, and
zone boundaries are at km scale, floating-point edge cases are extremely
unlikely. No action needed.

### 5.2 [INFO] GatherToRootDedup uses depth-only merging for BP5

`GatherToRootDedup` and `GatherFieldsToRootDedup` sort and merge duplicates
by depth (z-coordinate) only. For BP5's 2D fault surface, DOFs at the same
depth but different along-strike coordinates would be incorrectly merged.

Not a bug for serial (no partition-boundary duplicates) or for Cartesian
partitioning where each partition boundary has unique depths. Could be
an issue for non-Cartesian partitioning in parallel.

This is acknowledged in the plan §Parallel Considerations:
"Future: Full 2D (x2, x3) deduplication may be needed for non-Cartesian
partitioning."

### 5.3 [INFO] tau0\_vec uses Vi\_abs vs Tandem's Vi2

See §4.5. MFEM uses the full velocity magnitude; Tandem uses the second
component only. Difference is negligible (< machine epsilon) since
V[0] = 1e-20 << V[1].

---

## 6. Test Coverage Assessment

### BP5Params unit tests (test\_bp5\_params.cpp): EXCELLENT ✓

**12 test functions**, ~35+ assertions covering:
- [x] Material properties (mu, lambda, eta) — TestMaterialProperties
- [x] a(x2,x3): VW core, VS zones, transition, edges — TestAFunction
- [x] Nucleation zone: inside/outside, corners — TestNucleationAndL
- [x] V\_init\_vec: nucleation vs outside, magnitude — TestVinitVec
- [x] tau0\_vec: magnitude, direction, self-consistency with friction
      solver — TestTau0Vec
- [x] Print output — TestPrint
- [x] a() depth symmetry within VW zone — TestADepthSymmetry
- [x] Zone boundary edge cases (exact VW/VS boundaries) — TestZoneBoundaryEdgeCases
- [x] tau0\_vec in VS zone with self-consistency — TestTau0VecVSZone
- [x] tau0\_vec direction (anti-parallel to V\_init) — TestTau0VecDirection
- [x] psi\_init() steady-state value — TestPsiInit
- [x] Validate() on default params — TestValidate

### FaultGeometry 3D tests: COMPLETE ✓

**2 test functions** in test\_elasticity\_operator.cpp:

**TestFaultGeometry3D** (lines 310-353, 6 assertions):
- [x] DOF count matches domain operator
- [x] IsBP5() returns true
- [x] GetAValues() size == N
- [x] GetEtaValues() size == N, values match eta()
- [x] GetDcValues() size == N
- [x] GetTauPre() size == 2N

**TestFaultGeometry3DValues** (lines 358-431, ~13 assertions) *(NEW)*:
- [x] GetVInit() size == 2N and non-zero
- [x] GetCoordsX2() size == N, values within mesh bounds
- [x] GetCoordsX3() size == N, values within mesh bounds
- [x] a\_values all in [a0, amax]
- [x] dc\_values all equal to L0 or L\_nuc
- [x] tau\_pre non-zero, physically reasonable magnitude
- [x] GetBP5Params() returns correct b value
- [x] Print() for BP5 path outputs "VW DOFs" (exercises §2.1 fix)

### Previous test gaps — all resolved:

1. **[COMPLETED] GetVInit()**: Size and non-zero norm verified (lines 381-383)
2. **[COMPLETED] GetCoordsX2/X3**: Size and bounds verified (lines 386-394)
3. **[COMPLETED] a\_values range**: Per-DOF `[a0, amax]` check (lines 397-402)
4. **[COMPLETED] dc\_values range**: Per-DOF L0 or L\_nuc check (lines 405-411)
5. **[COMPLETED] tau\_pre values**: Non-zero and magnitude check (lines 414-418)
6. **[COMPLETED] Print() BP5 path**: Output contains "VW DOFs" (lines 425-429)

---

## 7. Backward Compatibility

### 7.1 BP2 constructor: Unchanged ✓

The existing BP2Params constructor (lines 49-93) is unmodified. All existing
BP2 tests pass without changes.

### 7.2 Data member initialization: CORRECT ✓

- `is_bp5_` defaults to `false` (line 416), set to `true` only in BP5 constructor
- `bp5_params_` is default-constructed (empty struct) for BP2 path
- BP5-only vectors (`dc_values_`, `tau_pre_`, `V_init_vec_`, `coords_x2_`,
  `coords_x3_`) are only populated in ComputeBP5Params()

### 7.3 Shared accessors: CORRECT ✓

- `GetAValues()`, `GetEtaValues()`, `GetDepths()` work for both BP2 and BP5
- `NumFaultDOFs()`, `NumLocalFaultDOFs()`, `NumGlobalFaultDOFs()` work for both
- `FindNearestDOF()` uses `depths_` which is populated for both paths

---

## 8. Summary Table

| Item | Status | Severity | Section |
|------|--------|----------|---------|
| Print() uses BP2 params\_.b for BP5 | **COMPLETED** | ~~MEDIUM~~ | §2.1 |
| IsVelocityWeakening/GetVWDOFs/GetVSDOFs | **COMPLETED** | ~~MEDIUM~~ | §2.2 |
| GetVWDepth/GetVSDepth guards | **COMPLETED** | ~~LOW~~ | §2.3 |
| V\_init\_vec/tau0\_vec comment labels | **COMPLETED** | ~~MEDIUM-DOCS~~ | §2.4 |
| Nucleation zone no eps tolerance | — | INFO | §5.1 |
| Depth-only dedup for 2D fault | — | INFO | §5.2 |
| tau0\_vec Vi\_abs vs Vi2 | — | INFO | §5.3 |
| a(x2,x3) formula | Correct | — | §4.1 |
| IsNucleationZone | Correct | — | §4.2 |
| L\_of\_x2\_x3 / Dc\_of\_x2\_x3 | Correct | — | §4.3 |
| V\_init\_vec values | Correct | — | §4.4 |
| tau0\_vec values | Correct | — | §4.5 |
| eta computation | Correct | — | §4.6 |
| Material parameters (all 17) | Correct | — | §4.7 |
| ComputeBP5Params | Correct | — | §4.8 |
| GetFaultCoords2D | Correct | — | §4.9 |
| BP5Params unit tests (12 functions) | Excellent | — | §6 |
| FaultGeometry 3D tests (2 functions) | **Complete** | — | §6 |
| GetVInit tested | **COMPLETED** | ~~MEDIUM~~ | §6 |
| GetCoordsX2/X3 tested | **COMPLETED** | ~~MEDIUM~~ | §6 |
| a\_values range verified | **COMPLETED** | ~~LOW~~ | §6 |
| dc\_values range verified | **COMPLETED** | ~~LOW~~ | §6 |
| tau\_pre values verified | **COMPLETED** | ~~LOW~~ | §6 |
| Print BP5 path tested | **COMPLETED** | ~~LOW~~ | §6 |
| BP2 backward compatibility | Correct | — | §7 |

---

## Remaining items summary

| Category | Count | Details |
|----------|-------|---------|
| Open code bugs | **0** | All 4 resolved (§2.1–§2.4) |
| Info items | **3** (§5.1, §5.2, §5.3) |
| Missing tests | **0** | All 6 gaps covered by TestFaultGeometry3DValues |
| Correctness checks | **All 9 passed** (§4.1–§4.9) |
| Backward compatibility | **Verified** (§7) |

**Compared to revision 1**: All 4 bugs moved to COMPLETED, all 6 test gaps
covered. Phase 2c implementation is fully complete with no open bugs and
comprehensive test coverage (12 BP5Params tests + 2 FaultGeometry 3D tests
with ~19 total assertions). Ready for integration testing with the SEAS
time-stepper.
