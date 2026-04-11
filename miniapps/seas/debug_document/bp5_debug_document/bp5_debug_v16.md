# BP5 Debug v16: Enhanced Diagnostic VTK + Region-Conformal Mesh

## Context

After implementing H10 (all-boundary Dirichlet) and H11 (configurable nucleation), we need:
1. **Enhanced diagnostic VTK** — Include on-fault slip rate, slip, state variable, and traction (not just "a" and tau_pre)
2. **Region-conformal Gmsh mesh** — Pre-create surfaces for VW, VS, and transition boundaries so element faces align with region boundaries and parameters are interpolated correctly
3. **Unit tests** for all new capabilities

## Part 1: Enhanced Diagnostic VTK Output (PVD format)

### Current State

The `--diag-vtk` flag in `bp5_verification_full.cpp` currently outputs:
- `displacement` — Volume-mesh displacement from BC solve at t=1yr
- `fault_a` — Rate-and-state "a" mapped to fault-adjacent elements (from analytic formula at centroids)
- `tau_pre_magnitude` — Pre-stress magnitude (from analytic formula)
- `V_init_magnitude` — Initial velocity magnitude (from analytic formula)
- `boundary_attributes` — Separate VTU with boundary attribute IDs

### Problem

These fields come from **analytic evaluation at element centroids**, not from the actual simulation state. Missing: slip, slip rate, state variable (psi), and traction — the quantities computed during `SetInitialCondition()`.

### Design

Enhance `--diag-vtk` to also output fault quantities from the **actual simulation state** after initialization. Use the `ParaViewOutput::InitFaultOutput()` pattern (fault face → element mapping via `GetFaultInteriorFaces()`).

**New fields** (all mapped to fault-adjacent volume elements via L2 p=0):

| Field Name | Source | Size per DOF | Notes |
|---|---|---|---|
| `slip_dip` | `fault_op.GetSlip(state, slip)` → `slip(2*i)` | scalar | Should be 0 at init |
| `slip_strike` | `slip(2*i+1)` | scalar | Should be 0 at init |
| `V_dip` | `fault_op.GetSlipRate()` → `V(2*i)` | scalar | ≈ V_zero at init |
| `V_strike` | `V(2*i+1)` | scalar | V_init or V_nuc |
| `V_magnitude` | `sqrt(V_dip² + V_strike²)` | scalar | Key: shows nucleation zone |
| `tau_dip` | `seas_op.GetTraction()` + `fault_geom.GetTauPre()` | scalar | Total stress = tau_pre + elastic |
| `tau_strike` | same | scalar | Dominant component |
| `tau_magnitude` | magnitude of total stress | scalar | |
| `state_psi` | `fault_op.GetTheta(state, theta)` → convert to psi | scalar | Initial psi_ss everywhere |
| `fault_a` | `fault_geom.GetAValues()` | scalar | Already exists, keep |
| `fault_dc` | `fault_geom.GetDcValues()` | scalar | L0 vs L_nuc |

### Implementation

**File: `tests/verification/bp5_verification_full.cpp`** — Enhance the `if (diag_vtk)` block.

Instead of computing parameters from analytic formulas at element centroids (which is approximate), use the **fault face → element mapping** from `domain.GetFaultInteriorFaces()`:

```cpp
// Build fault face → element mapping (same pattern as ParaViewOutput)
const Array<int> &fault_int_faces = domain.GetFaultInteriorFaces();
int nf_int = fault_int_faces.Size();
std::vector<int> face_elem1(nf_int), face_elem2(nf_int);
for (int i = 0; i < nf_int; i++)
{
   FaceElementTransformations *FTr =
      pmesh.GetInteriorFaceTransformations(fault_int_faces[i]);
   face_elem1[i] = FTr->Elem1No;
   face_elem2[i] = FTr->Elem2No;
}

// Extract fault quantities from state
Vector slip, theta;
fault_op.GetSlip(state, slip);       // 2*N interleaved [dip, strike]
fault_op.GetTheta(state, theta);     // N (converted from psi)
const Vector &V = fault_op.GetSlipRate();  // 2*N interleaved
const Vector &traction = seas_op.GetTraction();  // 2*N interleaved
const Vector &tau_pre = fault_geom.GetTauPre();  // 2*N interleaved
const Vector &a_vals = fault_geom.GetAValues();  // N
const Vector &dc_vals = fault_geom.GetDcValues(); // N

// Create L2 p=0 fields and map fault DOFs to adjacent elements
// For each fault interior face i:
//   Both face_elem1[i] and face_elem2[i] get the same fault DOF value
//   DOF index i maps to slip(2*i), slip(2*i+1), theta(i), etc.
```

**Key detail — DOF ordering**: Fault DOFs are ordered as interior faces first, then shared faces. For the diagnostic VTK, we only map interior faces (index 0..nf_int-1). Shared face DOFs (parallel partition boundaries) won't appear in the VTK — acceptable for diagnostic purposes since each shared face exists as an interior face on one rank.

**Keep existing centroid-based fields** (fault_a, tau_pre_magnitude, V_init_magnitude) as they provide useful reference, and **add** the new state-based fields alongside.

### Parallel Considerations

- `ParGridFunction` on `ParFiniteElementSpace(L2 p=0)` handles MPI automatically
- Each rank maps only its local interior faces → no communication needed
- `ParaViewDataCollection` on `ParMesh` writes per-rank VTU files + PVTU metadata
- Output format: PVD (`.pvd` + `.vtu` files), viewable in ParaView

## Part 2: Region-Conformal Gmsh Mesh

### Current State

The existing `bp5.geo` already embeds 4 surfaces via BooleanFragments:
- `fault` — Full fault plane: y ∈ [-50, 50] km, z ∈ [0, 40] km
- `nuc1` — Outer transition: y ∈ [-32, 32], z ∈ [2, 18] km
- `nuc2` — VW core: y ∈ [-30, 30], z ∈ [4, 16] km
- `nuc3` — Nucleation patch: y ∈ [-30, -18], z ∈ [4, 16] km

These create conformal boundaries at VW/VS/transition edges. However:
1. The surfaces are used only for mesh refinement, not region identification
2. No Physical Surface tags for individual regions
3. Missing surface for the shallow VS strip (z ∈ [0, 2]) and deep VS strip (z ∈ [18, 40])

### Design: `bp5_v2.geo`

Create a new Gmsh file `bp5/mesh/bp5_v2.geo` that:

1. **Preserves** the current domain geometry and boundary attributes (1-6)
2. **Adds** Physical Surface tags for fault sub-regions (attrs 101-106)
3. **Adds** separate surfaces for shallow VS and deep VS zones

**Region Physical Surface tags:**

| Tag | Region | y range (km) | z range (km) | a value |
|-----|--------|---------|---------|---------|
| 100 | Full fault (all) | [-50, 50] | [0, 40] | varies |
| 101 | Shallow VS | [-50, 50] | [0, 2] | amax |
| 102 | Shallow transition | [-32, 32] | [2, 4] | varies |
| 103 | VW core (non-nuc) | [-30, 30] | [4, 16] | a0 |
| 104 | Deep transition | [-32, 32] | [16, 18] | varies |
| 105 | Deep VS | [-50, 50] | [18, 40] | amax |
| 106 | Nucleation zone | [-30, -18] | [4, 16] | a0, L_nuc |

Lateral VS (|y| > 32, z ∈ [2, 18]) and lateral transition (|y| ∈ [30, 32]) are covered by the BooleanFragments partition of the fault surface.

**Key geometry dimensions (from BP5 params):**
```
h_s = 2 km      → shallow VS width
h_t = 2 km      → transition width
H   = 12 km     → VW core height
l   = 60 km     → VW along-strike extent (l_vw)
w   = 12 km     → nucleation patch width (w_nuc)
W_f = 40 km     → fault depth extent
l_f = 100 km    → fault along-strike extent
```

**Critical boundaries on fault plane:**
- z = 0: top surface
- z = 2 km (h_s): VS ↔ transition
- z = 4 km (h_s+h_t): transition ↔ VW
- z = 16 km (h_s+h_t+H): VW ↔ transition
- z = 18 km (h_s+2h_t+H): transition ↔ VS
- z = 40 km (W_f): bottom of fault
- y = ±30 km (l/2): VW ↔ transition along-strike
- y = ±32 km (l/2+h_t): transition ↔ VS along-strike
- y = -18 km (-l/2+w): nucleation zone right boundary

**Implementation:**

```
// Additional surfaces for complete fault partitioning
// Shallow VS strip: full width, above nuc1
vs_shallow = news;
Rectangle(vs_shallow) = {0, -l_f/2, 0, h_s, l_f};
Rotate{ {0, 1, 0}, {0, 0, 0}, -Pi/2} { Surface{vs_shallow}; }

// Deep VS strip: full width, below nuc1
vs_deep = news;
Rectangle(vs_deep) = {h_s+2*h_t+H, -l_f/2, 0, W_f-(h_s+2*h_t+H), l_f};
Rotate{ {0, 1, 0}, {0, 0, 0}, -Pi/2} { Surface{vs_deep}; }

// Include all surfaces in BooleanFragments
BooleanFragments{ Volume{1,2}; Delete; }
   { Surface{fault,nuc1,nuc2,nuc3,vs_shallow,vs_deep}; Delete; }

// After BooleanFragments, identify sub-regions by bounding box
// and assign Physical Surface tags 101-106
```

**Mesh sizing** (same as current, with explicit control per region):
- VW core (nuc2): `res_f` (1 km for benchmark)
- Transition (nuc1 \ nuc2): `res_f * 1.5`
- Nucleation (nuc3): `res_f` (finest)
- VS zones: `res_f * 2` near fault, grading to `res` far from fault
- Volume: smooth grading via Distance + Threshold fields

**Python generation script** `bp5/mesh/generate_bp5_mesh.py`:
```python
#!/usr/bin/env python3
"""Generate BP5 mesh with region-conformal surfaces."""
import subprocess, argparse

def generate(res_f_km, output):
    cmd = f"gmsh -3 bp5_v2.geo -setnumber res_f {res_f_km} -o {output}"
    subprocess.run(cmd.split(), check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--res", type=float, default=1.0, help="Fault resolution (km)")
    parser.add_argument("-o", "--output", default="bp5_v2_1000m.msh")
    args = parser.parse_args()
    generate(args.res, args.output)
```

### Backward Compatibility

The new mesh is fully backward-compatible:
- Boundary attributes 1-6 unchanged
- Volume attribute 10 unchanged
- Internal fault surface (attr 100) unchanged
- New region tags (101-106) are additional Physical Surface groups — ignored by code that doesn't query them
- The existing `bp5.geo` remains for reference

## Part 3: Unit Tests

### Test 1: `test_diag_vtk.cpp` — Diagnostic VTK Output

Tests that the `--diag-vtk` code path produces correct output:

```
TestDiagVTKFileCreation:
  - Create inline mesh (nx=2, ny=2, nz=1)
  - Build domain + fault + SEAS operator stack
  - Call SetInitialCondition
  - Run diagnostic VTK output code
  - Verify: output directory exists, .pvd file exists, .vtu files exist

TestDiagVTKFaultMapping:
  - Create inline mesh with known fault face count
  - Build fault face → element mapping
  - Verify: nf_int interior faces found, each has valid elem1/elem2
  - Verify: L2 p=0 field has correct size (NE elements)
  - Map a test value (1.0) to fault-adjacent elements
  - Verify: fault-adjacent elements = 1.0, others = 0.0

TestDiagVTKFieldValues:
  - Create full BP5 stack, init
  - Extract slip, V, theta, traction from state
  - Map to L2 p=0 fields
  - Verify: initial slip = 0
  - Verify: V_magnitude > 0 on fault elements
  - Verify: state_psi > 0 on fault elements
  - Verify: tau fields are finite and physically reasonable (5-25 MPa)
  - Verify: a values match bp5_params.a_of_x2_x3() at face centroids
```

### Test 2: `test_bp5_mesh_regions.cpp` — Mesh Region Verification

Tests that mesh generation produces correct region tagging (also verifiable with inline mesh + analytic check):

```
TestInlineMeshBoundaryAttributes:
  - Create inline mesh
  - Verify all 6 boundary attributes present
  - Verify face counts per attribute

TestFaultRegionParameters:
  - Create inline mesh with sufficient resolution (nx=4, ny=8, nz=4)
  - Build FaultGeometry
  - For each fault DOF, verify:
    - a(x2, x3) matches expected region (VW/VS/transition)
    - Dc matches expected value (L_nuc in nucleation zone, L0 elsewhere)
    - tau_pre is finite and positive
  - Verify: at least one VW DOF exists (a = a0 = 0.004)
  - Verify: at least one VS DOF exists (a = amax = 0.04)

TestFaultFaceElementMapping:
  - Create inline mesh, build domain operator
  - Get fault interior faces
  - For each, verify elem1 and elem2 are valid element indices
  - Verify elem1 != elem2
  - Verify face centroid has x ≈ 0
```

### Makefile Changes

Add to `Makefile`:
```makefile
# New test sources
TEST_DIAG_VTK_SRC = tests/unit/test_diag_vtk.cpp
TEST_DIAG_VTK_OBJ = $(TEST_DIAG_VTK_SRC:.cpp=.o)
TEST_MESH_REGIONS_SRC = tests/unit/test_bp5_mesh_regions.cpp
TEST_MESH_REGIONS_OBJ = $(TEST_MESH_REGIONS_SRC:.cpp=.o)

# Add to SEQ_MINIAPPS
seas_test_diag_vtk seas_test_mesh_regions

# Executable + compile rules (follow existing pattern)
# Test targets:
test-diag-vtk: seas_test_diag_vtk
	./seas_test_diag_vtk
test-mesh-regions: seas_test_mesh_regions
	./seas_test_mesh_regions
```

## Part 4: Update Existing Files

### `bp5.geo` Header Comment
Update boundary attribute comments to reflect H10 (all Dirichlet):
```
//   5 = z = 0    (Dirichlet: u_y = sgn(x)*Vp*t/2, matching Tandem)
```

## Files to Create/Modify

| File | Action | Description |
|------|--------|-------------|
| `tests/verification/bp5_verification_full.cpp` | **Modify** | Enhance `--diag-vtk`: add fault state fields via face→element mapping |
| `bp5/mesh/bp5_v2.geo` | **Create** | New mesh with region-conformal surfaces + Physical Surface tags |
| `bp5/mesh/generate_bp5_mesh.py` | **Create** | Python wrapper for mesh generation |
| `tests/unit/test_diag_vtk.cpp` | **Create** | Unit tests for diagnostic VTK output |
| `tests/unit/test_bp5_mesh_regions.cpp` | **Create** | Unit tests for mesh regions and parameter mapping |
| `Makefile` | **Modify** | Add new test targets |
| `bp5/mesh/bp5.geo` | **Modify** | Update BC comment to match H10 |
| `debug_document/bp5_debug_document/bp5_debug_v16.md` | **Create** | This document |

## Verification

### Step 1: Build and run ALL unit tests
```bash
conda activate mfem-dev
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas
make -j test-bp5-integration test-bp5-output test-diag-vtk test-mesh-regions
./test-bp5-integration && ./test-bp5-output && ./test-diag-vtk && ./test-mesh-regions
```

### Step 2: Generate new mesh
```bash
conda activate pythonenv
cd bp5/mesh
gmsh -3 bp5_v2.geo -setnumber res_f 1 -o bp5_v2_1000m.msh
```

### Step 3: Run diagnostic VTK on new mesh
```bash
conda activate mfem-dev
make -j seas_bp5_full
mpirun -np 4 ./seas_bp5_full --mesh bp5/mesh/bp5_v2_1000m.msh \
   --tfinal 0 --diag-vtk --output-dir bp5/diag_v2
```

### Step 4: Visual verification in ParaView
Open `bp5/diag_v2/bp5_diag/bp5_diag.pvd`:
1. **displacement**: Threshold on magnitude > 0 → verify all-boundary Dirichlet pattern
2. **fault_a**: Threshold > 0 → verify VW core (a=0.004), VS (a=0.04), smooth transition
3. **V_magnitude**: Threshold > 0 → verify nucleation zone has V ≈ 0.03, elsewhere V ≈ 1e-9
4. **fault_dc**: Threshold > 0 → verify L_nuc = 0.13 in nucleation zone, L0 = 0.14 elsewhere
5. **tau_magnitude**: Threshold > 0 → verify physically reasonable (15-20 MPa range)
6. **state_psi**: Threshold > 0 → verify uniform initial psi

Open `bp5/diag_v2/boundary_attributes/boundary_attributes.pvtu`:
- Color by Attribute → verify 6 distinct boundary faces

## Previously Implemented Fixes

- H1-H7: Various friction, output, penalty fixes ✅
- H8: Bottom-only Dirichlet BC ✅ (superseded by H10)
- H9: Output floor removed ✅
- H10: All-boundary Dirichlet loading ✅
- H11: Configurable nucleation parameters ✅
- H12: Updated CLAUDE.md BC documentation ✅
- Diagnostic VTK (basic): displacement, fault_a, tau_pre, V_init, boundary_attrs ✅
- `--tfinal 0` fix ✅

## Implementation Results (2026-03-13)

### All Changes Applied

1. **`bp5_verification_full.cpp`** — Enhanced `--diag-vtk` block with 11 new state-based fields via fault face→element mapping. Kept original centroid-based analytic fields.

2. **`bp5/mesh/bp5_v2.geo`** — New Gmsh mesh with:
   - 2 additional surfaces (`vs_shallow`, `vs_deep`) in BooleanFragments
   - 6 Physical Surface tags (101-106) for fault sub-regions
   - Updated header for H10 (all-boundary Dirichlet)

3. **`bp5/mesh/generate_bp5_mesh.py`** — Python wrapper for mesh generation.

4. **`tests/unit/test_diag_vtk.cpp`** — 3 tests, 22 assertions, all passing:
   - `TestFaultFaceElementMapping`: Validates face→element mapping
   - `TestDiagVTKFieldValues`: Validates state extraction after init
   - `TestL2FieldMapping`: Validates L2 p=0 consistency

5. **`tests/unit/test_bp5_mesh_regions.cpp`** — 3 tests, 19 assertions, all passing:
   - `TestBoundaryAttributes`: All 6 boundary attributes present
   - `TestFaultRegionParameters`: a ∈ [a0, amax], Dc > 0, tau_pre finite
   - `TestFaultFaceElementMapping`: Valid indices, distinct elements, x≈0

6. **`Makefile`** — Added `seas_test_diag_vtk`, `seas_test_mesh_regions` targets.

7. **`bp5/mesh/bp5.geo`** — Updated BC comment for z=0 (Dirichlet, H10 fix).

### VTK Fields in Output

```
displacement, fault_a, tau_pre_magnitude, V_init_magnitude,
slip_dip, slip_strike, V_dip, V_strike, V_magnitude,
tau_dip, tau_strike, tau_magnitude, state_psi, fault_a_state, fault_dc
```

### Test Results

- `seas_test_diag_vtk`: 22/22 passed
- `seas_test_mesh_regions`: 19/19 passed
- `seas_test_bp5_integration`: 33/33 passed (no regression)
- `seas_test_bp5_output`: 45/45 passed (no regression)
- `seas_test_bp5_mesh`: 25/25 passed (no regression)
- `seas_bp5_full --diag-vtk`: Runs successfully on 2 ranks, all fields present in VTK
