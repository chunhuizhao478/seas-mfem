# BP5 Debug v30: Adopt Tandem Coordinate System Exclusively

**Date**: 2026-03-14
**Status**: Implemented, all unit tests pass. Ready for IP simulation run.
**Previous**: v29 (BR2+new mesh still shows VS lockup, IP still blows up)

---

## 1. Motivation

The dual MFEM/Tandem coordinate system created subtle sign bugs and dead-code branches.
After v24-v29 showed persistent issues (IP blowup, VS lockup), we eliminate the dual
system entirely. Use Tandem's mesh directly. No more coordinate transformations.

## 2. Coordinate System (Tandem, now exclusive)

| Axis | SCEC Abstract | Tandem Mesh | Old MFEM (removed) |
|------|---------------|-------------|---------------------|
| Along-strike | x2 | X | y |
| Fault-normal | x1 | Y | x |
| Depth | x3 (positive down) | −Z (Z negative down) | z (positive down) |
| Fault plane | x1 = 0 | Y = 0 | x = 0 |
| Free surface | x3 = 0 | Z = 0 | z = 0 |
| Domain depth | x3 = 100 km | Z = −100 km | z = 100 km |

**Boundary tags** (Tandem Physical Surface):
- 1 = Natural (top Z=0 + bottom Z=Z₀), zero traction
- 3 = Fault (Y=0 interior), handled by DG slip BC
- 5 = Dirichlet (far-field X=±X₁, Y=±Y₁), plate loading

**Loading**: u_D = (sgn(Y)·Vₚ·t/2, 0, 0) — along-strike (X) component only.

## 3. Changes Made

### 3.1 `fault/fault_basis.hpp` — Remove dip negation

**Old**: `dip = −(strike × n)` (negated for MFEM z-positive-down convention)
**New**: `dip = strike × n` (no negation, matches Tandem `Curvilinear.cpp:facetBasis()`)

With ref_normal=(0,−1,0), up=(0,0,1):
- n = (0, −1, 0)
- strike = up × n = (1, 0, 0)
- dip = strike × n = (0, 0, −1) — points downward ✓

### 3.2 `domain/seas_boundary_tags.hpp` — Update tag constants

Old: `FARFIELD_LEFT=1, FARFIELD_RIGHT=2, FREE_SURFACE=3, BOTTOM=4, FAULT=5`
New: `NATURAL=1, FAULT=3, DIRICHLET=5`

Fixed `FindFaultSharedFaces()`: `coords(0)` → `coords(1)` (fault at Y=0).

### 3.3 `domain/elasticity_operator.hpp` — Major cleanup (13 changes)

| Change | Detail |
|--------|--------|
| Remove `fault_normal_axis_` | Eliminated dual-system variable (was hardcoded to 0) |
| ref_normal | `(0,1,0)` → `(0,−1,0)` — matches Tandem |
| `SetupBoundaryMarkers` | Simplified: only mark attr 5 as Dirichlet |
| `BuildFaultTaggedFaces` | Only support attr 3 (removed attr 100) |
| `IsFaultFace3D` | Removed `mfem_check` branch, Tandem-only: Y≈0, X∈[−lf/2,lf/2], Z∈[−Wf,0] |
| `IsFaultFace3DShared` | Same simplification |
| Sign computation (×4) | `nor(fault_normal_axis_)` → `nor(1)` |
| BR2 sign bug (×2) | `nor_q(0)` → `nor_q(1)` — was hardcoded to wrong axis |
| `AssembleDirichletLoading` | Removed if/else: always `u_D[0] = sgn(Y)·Vₚ·t/2` |
| `GetFaultDepths` | Always `depth = −coords(2)`, fixed shared-faces bug (`coords(2)` → `−coords(2)`) |
| `GetFaultCoords2D` | Always `x2=coords(0)`, `x3=−coords(2)`, fixed shared-faces bug |

### 3.4 `config/bp5_params.hpp` — tau_pre sign fix (CRITICAL)

**Bug found**: `tau0_vec()` returned positive tau_pre, but Tandem returns negative.

**Physics**: With ref_normal=(0,−1,0), the elastic traction for right-lateral loading is
negative in the strike direction (T = σ·n, T_x = −σ_xy < 0). The pre-stress `tau_pre` must
have the **same sign** so that `tau_total = tau_pre + tau_elastic` reinforces:

| | `tau_pre` | `tau_elastic` (right-lateral) | `tau_total` behavior |
|---|---|---|---|
| Tandem | −τ₀ | −τ_el | Both negative → |τ_total| grows → nucleation ✓ |
| MFEM (old, bug) | +τ₀ | −τ_el | Opposite signs → |τ_total| initially shrinks ✗ |
| MFEM (fixed) | −τ₀ | −τ_el | Matches Tandem ✓ |

**Fix**: `tau[i] = -tau0_scalar * Vi[i] / Vi_abs` (negate, matching Tandem's `bp5.lua`).

This was a **pre-existing bug** (not introduced by the coordinate change). It would cause
delayed earthquake nucleation and wrong traction balance at t=0.

### 3.5 Test file updates

| File | Changes |
|------|---------|
| `test_fault_basis.cpp` | Updated expected dip from (0,0,−1) to (0,0,+1) for ref_normal=(1,0,0) test |
| `test_elasticity_operator.cpp` | Mesh Z-range [−Lz,0], Tandem boundary tags (1,5), lf=2Lx, traction sign assertions |
| `test_fault_detection.cpp` | Updated `SEASBoundaryTags` constant names |
| `test_bp5_params.cpp` | Updated tau_pre sign assertions |
| `bp5_verification_full.cpp` | Inline mesh Z-range [−Lz,0], Tandem boundary tags |

### 3.6 Other files

- `bp5/mesh/bp5.geo` — Added deprecation notice (uses old MFEM convention)
- `config/bp5_params.hpp` — Updated coordinate mapping comments

## 4. Term-by-Term Verification vs Tandem

### 4.1 Fault Basis

| | Tandem (`Curvilinear.cpp:facetBasis`) | SEAS-MFEM |
|---|---|---|
| ref_normal | (0, −1, 0) | (0, −1, 0) ✓ |
| up | (0, 0, 1) | (0, 0, 1) ✓ |
| strike | `u.cross(n).normalized()` = (1,0,0) | `up × n` = (1,0,0) ✓ |
| dip | `s.cross(n).normalized()` = (0,0,−1) | `strike × n` = (0,0,−1) ✓ |

### 4.2 DG Weak Form (IP method, SIPG ε=−1)

**Bilinear form** on interior fault face F with prescribed jump g = sign·δu:

| Term | Formula | Tandem | SEAS-MFEM |
|------|---------|--------|-----------|
| Symmetry | ε·½·∫_F {σ(v)·n}·g ds | `c1 = ε*0.5` | `ε * ip.w / (2*detJ)` ✓ |
| Penalty (IP) | κ·∫_F g·[[v]] ds | `(p0+p1)/4` | Same formula ✓ |
| Penalty (BR2) | η·∫_F C:r(g)⊗n·[[v]] ds | `NumFacets * M⁻¹ lifting` | Same ✓ |

**IP penalty formula** (both codes):
```
penalty = (p0 + p1) / 4
p_e = (D+1) · c_{N,1} · (face_area / elem_vol) · (c1² / c0)
c0 = 2μ, c1 = Dλ + 2μ
```

### 4.3 Sign Convention

| | Tandem | SEAS-MFEM |
|---|---|---|
| Slip sign | `side==1 ? −1 : +1` (per element) | `nor(1)>0 ? −1 : +1` (per face) |
| Elem1 penalty | `+penalty · sign · f_q` | `+penalty · sign · δu` ✓ |
| Elem2 penalty | `−penalty · sign · f_q` | `−penalty · sign · δu` ✓ |

Both are equivalent: the per-face sign with ± application to elem1/elem2 produces the
same result as per-element side-based sign.

### 4.4 Traction Computation

Both: `T = {σ}·n − penalty·([[u]] − g)` with scalar penalty.

### 4.5 Boundary Loading

Tandem (`bp5.lua`): `return sgn(Y) * Vp*t/2, 0, 0`
SEAS-MFEM: `u_D[0] = sgn(centroid(1)) * Vp * time / 2` ✓

### 4.6 tau_pre (fixed)

Both: `tau_total = tau_pre + tau_elastic`, with `tau_pre = −τ₀ · V/|V|` ✓

## 5. Test Results

| Test Suite | Tests | Result |
|------------|-------|--------|
| `seas_test_fault_basis` | 187 | All pass ✓ |
| `seas_test_elasticity_operator` | 84 | All pass ✓ |
| `seas_test_elasticity_br2` | 46 | All pass ✓ |
| `seas_test_fault_detection` | 7 | All pass ✓ |
| `seas_test_bp5_params` | 96 | All pass ✓ |
| `seas_test_domain_interface` | 32 | All pass ✓ |

## 6. Breaking Changes

- Old MFEM-convention meshes (attrs 1-6, 100) **no longer work**
- `bp5.geo` is deprecated; use Tandem's `bp5_tandem.geo` directly
- `BCMode::XOnly` now behaves same as `FarField` (Tandem lumps all far-field into tag 5)

## 7. Next: IP Simulation Run

Test the coordinate-system fix with an IP run on Tandem mesh:

```bash
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas
conda activate mfem-dev

# IP method, Tandem mesh, short run for validation
mpirun -np 8 ./seas_bp5_full \
  --mesh bp5/mesh/reference/bp5_tandem.msh \
  --mesh-scale 1000 \
  --dg-method IP \
  --solver mumps-blr \
  --tfinal 1e8 \
  --diag-vtk \
  --output-dir output_v30_ip \
  --output-prefix bp5_v30_ip \
  --check-residual
```

**What to look for in diagnostic VTK**:
1. Fault basis vectors: dip should point −Z (downward), strike should point +X
2. Displacement field: X-component loading, sgn(Y) pattern
3. Boundary attrs: only 1 (Natural) and 5 (Dirichlet) visible
4. tau_pre: negative in strike direction (right-lateral convention)
5. a-parameter distribution: VW core at correct depth/along-strike location

**Success criteria**:
- IP run survives beyond t > 1e8 s (~3 years) without blowing up
- Traction values physically reasonable (order of sigma_n * f ≈ 15 MPa)
- If IP still blows up, try BR2 with same mesh to isolate cause
