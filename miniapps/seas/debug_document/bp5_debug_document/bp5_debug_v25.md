# BP5 Debug v25: Mesh Quality Fix + Traction Formula Fix

**Date**: 2026-03-14
**Status**: Two fixes applied, test runs designed
**Previous**: v24 (traction formula mismatch identified, mesh quality issue discovered)

---

## 1. Summary of Two Fixes

### Fix H26: Traction formula — IP-style penalty matching Tandem

**Problem**: MFEM's `ComputeTraction` used BR2 lifting for the penalty correction,
amplifying the DG residual by μ/h ≈ 3.2×10⁷. The benchmark document (Algorithm
line 15) specifies traction = `{{C:∇u}}·n̂` (average stress only). Tandem uses
an IP-style penalty with `penalty = NumFacets = 4` (dimensionless), which is
numerically negligible.

**Fix**: Replace the BR2 lifting in `ComputeTraction` with Tandem's IP-style:
```cpp
correction[c] = penalty_val * jump[c];  // penalty_val = 4 for tet
```

**Files**: `domain/elasticity_operator.hpp`, both interior and shared fault face
traction paths.

### Fix H27: Mesh quality — uniform fault resolution matching Tandem

**Problem**: MFEM's `bp5.geo` applied different MeshSize values to different
nucleation sub-zones (fault: res_f, nuc1: 2×res_f, nuc2: 1.5×res_f, nuc3: res_f)
plus a Distance/Threshold background field. This created non-uniform elements on
the fault surface with size transitions at zone boundaries (visible in ParaView).

Five fault faces near the nuc2 boundary (y≈28-29 km, z≈5-6 km) had pathological
average stress (~1 GPa) due to poor element quality at the size transition. With
the old BR2 traction correction, these were masked (correction subtracted ~2.4 GPa).
With fix H26 (IP-style), the raw average stress was exposed, causing traction blowup
and single-DOF runaway (slip reaching 98m in 2 years).

**Root cause**: The blowup DOFs were all clustered at y=28-29 km, 1-2 km inside
the VW/nuc2 boundary at y=30 km, exactly where the mesh sizing transitioned from
res_f to 1.5×res_f.

**Fix**: Rewrote `bp5.geo` to match Tandem's mesh sizing approach:
- **Uniform** `res_f` on ALL fault surface points (no per-zone sizing)
- **No background field** (no Distance/Threshold grading)
- Gmsh handles volume grading naturally
- Nuc zones kept for BooleanFragments geometry only, not mesh sizing

**Comparison**:

| Aspect | Old MFEM mesh | New MFEM mesh | Tandem mesh |
|--------|---------------|---------------|-------------|
| Fault sizing | res_f + 1.5× + 2× per zone | Uniform res_f | Uniform res_f |
| Volume grading | Background field | Natural Gmsh | Natural Gmsh |
| Elements q < 0.10 | Some | **0** | ~0 |
| Elements q < 0.20 | Some | **0** | ~0 |
| Volume elements | 66,220 | 62,992 | 63,451 |
| Nodes | 11,677 | 11,220 | 11,306 |

The new mesh is very close to Tandem's in element count and quality.

**File**: `bp5/mesh/bp5.geo`

---

## 2. Mesh Generation

The corrected generation command:
```bash
gmsh -3 bp5/mesh/bp5.geo -setnumber res_f 1 -setnumber res 40 -o bp5/mesh/bp5_1000m.msh
```

**Note**: Previous sbatch files used `-setnumber h 1.0`, which had no effect because
the geo file uses `res_f`, not `h`. The default `res_f` was 10 (10 km) in the old
geo file. The new geo file defaults to `res_f = 1` (1 km).

All sbatch files must be updated to use `-setnumber res_f 1` (or rely on the new
default).

---

## 3. Blowup Location Analysis (from v24)

The five blowup faces were all at:
```
(0, 28060, 6200)   y=28.1 km, z=6.2 km   VW core, 1.9 km from nuc2 edge
(0, 28102, 5574)   y=28.1 km, z=5.6 km   VW core, 1.9 km from nuc2 edge
(0, 28631, 5465)   y=28.6 km, z=5.5 km   VW core, 1.4 km from nuc2 edge
(0, 28894, 5055)   y=28.9 km, z=5.1 km   VW core, 1.1 km from nuc2 edge
(0, 28925, 5958)   y=28.9 km, z=6.0 km   VW core, 1.1 km from nuc2 edge
```

All within 1-2 km of the nuc2 boundary at y=30 km — exactly where the old mesh
had a sizing transition from res_f to 1.5×res_f. The new uniform mesh eliminates
this transition.

---

## 4. Test Plan

### Test 1: Uniform fault with both fixes (mesh + traction)

Verifies the VS zone maintains V ≈ Vp without earthquake, with both fixes active.

```
sbatch: bp5_v25_test1_uniform.sbatch
Flags:  --delta-tau-factor 0 --V-nuc 1e-9
Mesh:   bp5_1000m.msh (regenerated with uniform res_f=1)
```

**Pass**: V/Vp stays in [0.5, 2.0] at z=22 km for 300 years. No traction blowup.
**Fail**: VS zone still drifts or blows up.

### Test 2: Full BP5 with Tandem-matching initialization

Full earthquake simulation with both fixes.

```
sbatch: bp5_v25_test2_tandem.sbatch
Flags:  --psi-init-mode tandem --V-nuc 0.01
Mesh:   bp5_1000m.msh (regenerated)
tfinal: 600 years
```

**Expected**: First event at ~150 yr, deep VS maintains V ≈ Vp, subsequent cycling.

### Test 3: Full BP5 with SCEC-correct initialization

Test whether SCEC-correct ψ + δτ works now that the mesh and traction are fixed.

```
sbatch: bp5_v25_test3_scec.sbatch
Flags:  (defaults — SCEC psi, V_nuc=0.03, delta_tau_factor=1)
Mesh:   bp5_1000m.msh (regenerated)
tfinal: 600 years
```

**Expected**: Immediate earthquake (from δτ), then VS zone recovers and cycles.

---

## 5. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H25 | Previous fixes (see v21-v23 docs) | Done |
| **H26** | **Traction: IP-style penalty matching Tandem (no BR2 lifting)** | **Done (v24)** |
| **H27** | **Mesh: uniform fault resolution, no per-zone sizing** | **Done (v25)** |

---

## 6. Key Lessons

1. **The BR2 traction correction was masking two separate issues**: a traction formula
   error (BR2 amplification by μ/h) AND poor mesh quality at zone boundaries.

2. **Mesh quality matters for DG**: Non-uniform element sizes on the fault surface
   create pathological average stress at transition faces. Tandem avoids this with
   uniform fault resolution.

3. **The VS zone lockup** (the original problem) was caused by the BR2 traction
   amplification biasing stress transfer between VS and VW zones. The mesh quality
   issue was secondary but prevented the traction fix from working alone.

4. **Both fixes are needed together**: The traction fix alone blows up on the old mesh.
   The mesh fix alone wouldn't resolve the VS lockup (BR2 amplification still present).

---

## 7. Files Modified

| File | Change |
|------|--------|
| `domain/elasticity_operator.hpp` | H26: BR2 traction → IP-style penalty (both interior and shared) |
| `bp5/mesh/bp5.geo` | H27: Uniform fault sizing, removed per-zone MeshSize and background field |
| `bp5/mesh/bp5_1000m.msh` | Regenerated with `res_f=1` using new geo |
| `config/bp5_params.hpp` | V_nuc >= V_init validation (for uniform test) |
| `fault/fault_geometry.hpp` | Print() shows local rank stats (global reduction reverted — deadlock) |
| `fault/rate_state_fault.hpp` | SetScecPsiInit() flag, Init() branching |
| `tests/verification/bp5_verification_full.cpp` | --psi-init-mode, --diag-vtk traction output |
