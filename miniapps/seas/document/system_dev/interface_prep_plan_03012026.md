# Phase 1a + 1b: Fault Physical Group Detection & Interface Cleanup

**Prepared:** 2026-03-01
**Branch:** `feature/interface-prep`
**Scope:** Preparatory refactoring for BP5 — add Gmsh Physical Curve tag-based fault detection and generalize DomainOperator interface for multi-component slip.

## Context

The SEAS miniapp currently identifies fault faces by coordinate checking (`|x| < tol`). This works for BP1/BP2 (vertical fault at x=0) but fails for BP5 (3D fault surface), BP3 (dipping fault), or any non-axis-aligned fault. Phase 1a adds Gmsh Physical Curve tag-based fault detection. Phase 1b adds interface methods to `DomainOperator` needed by BP5 (multi-component slip, 2D fault coords) and generalizes `SEASQuasiDynamicOperator`/`SEASBdrLoadOperator` to accept any domain operator type.

All changes are backwards-compatible. Existing code works unchanged with `fault_tag=-1` (default).

---

## Phase 1a: Fault Physical Group Detection

### Step 1a.1: Create `domain/seas_boundary_tags.hpp`

New file with:
```cpp
struct SEASBoundaryTags {
    static constexpr int FARFIELD_LEFT  = 1;
    static constexpr int FARFIELD_RIGHT = 2;
    static constexpr int FREE_SURFACE   = 3;
    static constexpr int BOTTOM         = 4;
    static constexpr int FAULT          = 5;  // NEW
};
```

Plus a static utility class `FaultBoundaryData`:
```cpp
class FaultBoundaryData {
public:
    // Scan boundary elements for fault_tag, return interior face indices
    static Array<int> FindFaultInteriorFaces(const Mesh &mesh, int fault_tag);
    // Parallel: scan shared faces (coordinate fallback for now)
    static Array<int> FindFaultSharedFaces(const ParMesh &pmesh, int fault_tag,
                                           real_t coord_tol = 1e-10);
};
```

`FindFaultInteriorFaces` implementation:
```cpp
for (int be = 0; be < mesh.GetNBE(); be++) {
    if (mesh.GetBdrAttribute(be) != fault_tag) continue;
    int face, ori;
    mesh.GetBdrElementFace(be, &face, &ori);
    if (mesh.GetInteriorFaceTransformations(face) != nullptr)
        result.Append(face);
}
```

For shared faces in parallel: MFEM doesn't expose boundary attributes on shared faces directly. Use a hybrid approach — tag-based for interior faces, coordinate-based fallback for shared faces (same as current). This is acceptable since shared faces at partition boundaries are a parallel artifact; the fault geometry is the same.

### Step 1a.2: Update `.geo` files

**`bp1/mesh/bp1.geo`** — add after line 110:
```geo
// FAULT interior interface at x=0 (entire fault + below-Wf)
Physical Curve(5) = {7, 8, 9, 10, 11, 12};
```
Lines 7-10 = fault (z=0 to -Wf), Lines 11-12 = below-Wf to bottom.

**`bp1/mesh/bp1_selfsimilar.geo`** — add after line 193:
```geo
// FAULT interior interface at x=0 (z=0 to z=-Zf, z=-Zf to z=-D)
Physical Curve(5) = {17, 18};
```
Lines 17-18 = vertical lines at x=0 (upper + lower blocks).

**`bp2/mesh/bp2.geo`** — add after line 93:
```geo
// FAULT interior interface at x=0
Physical Curve(5) = {7, 8, 9, 10, 11};
```
(BP2 has lines 7-10 for fault zone, line 11 for below-Wf to bottom.)

### Step 1a.3: Update `AntiplaneDomainOperator` constructor + `SetupFaultInfo`

**File:** `domain/antiplane_operator.hpp`

Add `fault_tag` parameter to constructor (default -1 for backward compat):
```cpp
AntiplaneDomainOperator(MeshType &mesh, int order, real_t mu, real_t Vp,
                        real_t Wf = 40.0e3, DGMethod method = DGMethod::IP,
                        int fault_tag = -1);
```

Add private member `int fault_tag_;`.

Modify `SetupFaultInfo()` to dual-path:
```cpp
void SetupFaultInfo() {
    fault_interior_faces_.SetSize(0);
    if (fault_tag_ >= 1) {
        // Tag-based detection
        fault_interior_faces_ = FaultBoundaryData::FindFaultInteriorFaces(
            mesh_, fault_tag_);
    } else {
        // Legacy coordinate-based detection
        for (int f = 0; f < mesh_.GetNumFaces(); f++) {
            if (mesh_.GetInteriorFaceTransformations(f) && IsFaultFace(f))
                fault_interior_faces_.Append(f);
        }
    }
    // Shared faces: same as current (coordinate-based) for both paths
    // ... existing shared face code unchanged ...
}
```

Add public accessor: `int GetFaultTag() const { return fault_tag_; }`

### Step 1a.4: Update `AntiplaneBdrLoadOperator` — same changes

**File:** `domain/antiplane_bdrload_operator.hpp`

Same constructor change (add `int fault_tag = -1`), same `SetupFaultInfo()` dual-path logic, same accessor.

### Step 1a.5: Update `.geo` comment headers

Update the comment at the top of each `.geo` file to document tag 5 = FAULT.

### Step 1a.6: Create unit test `tests/unit/test_fault_detection.cpp`

Test cases:
1. **`test_tag_based_detection`**: Load a small mesh with Physical Curve(5), create `AntiplaneDomainOperator` with `fault_tag=5`, verify correct number of fault faces found.
2. **`test_coordinate_based_fallback`**: Same mesh, `fault_tag=-1`, verify same fault faces found via coordinate check.
3. **`test_tag_vs_coordinate_match`**: Verify both methods produce identical fault face lists on BP1 mesh.
4. **`test_no_tag_in_mesh`**: Load old mesh without tag 5, `fault_tag=5` → 0 fault faces found (graceful fallback).
5. **`test_boundary_attributes_max`**: Verify `mesh.bdr_attributes.Max() >= 5` after loading tagged mesh.

### Step 1a.7: Update verification drivers (minimal)

**`tests/verification/bp1_verification_full.cpp`** and **`tests/verification/bp1_bdrload.cpp`**: No changes needed — they use default `fault_tag=-1`, so coordinate-based detection continues to work. Future runs can opt into tag-based by passing `fault_tag=5` (but this is not required for Phase 1a).

### Step 1a.8: Verification

- Regenerate `.msh` files from updated `.geo` files
- Run BP1 verification with both `fault_tag=-1` (legacy) and `fault_tag=5` (new)
- Confirm identical fault face counts, identical traction, identical time series output

---

## Phase 1b: DomainOperator Interface Cleanup

### Step 1b.1: Add virtual methods to `DomainOperator`

**File:** `domain/domain_operator.hpp`

Add three virtual methods with default implementations (non-breaking):

```cpp
/// Number of slip components per fault DOF (1=antiplane, 2=BP5 3D)
virtual int NumSlipComponents() const { return 1; }

/// Get 2D fault coordinates. Default: x2=0, x3=depths.
virtual void GetFaultCoords2D(Vector &coords_x2, Vector &coords_x3) const
{
    Vector depths;
    GetFaultDepths(depths);
    coords_x2.SetSize(depths.Size());
    coords_x2 = 0.0;
    coords_x3 = depths;
}

/// Get off-fault displacement at given points. Default: no-op.
virtual void GetOffFaultDisplacement(
    const std::vector<Vector> &points, Vector &displacements) const
{
    displacements.SetSize(0);
}
```

### Step 1b.2: Generalize `SEASQuasiDynamicOperator`

**File:** `solver/seas_operator.hpp`

Change from:
```cpp
template <typename MeshType = Mesh>
class SEASQuasiDynamicOperator : public TimeDependentOperator
{
public:
   using DomainOpType = AntiplaneDomainOperator<MeshType>;
```

To:
```cpp
template <typename MeshType = Mesh,
          typename DomainOpType = AntiplaneDomainOperator<MeshType>>
class SEASQuasiDynamicOperator : public TimeDependentOperator
{
public:
```

The rest of the class stays the same — `DomainOpType*` member, constructor, `Mult()`, etc. All existing instantiations (`SEASQuasiDynamicOperator<ParMesh>`) continue to work because the default template argument provides `AntiplaneDomainOperator<ParMesh>`.

Also add `#include "../domain/domain_operator.hpp"` (already included transitively but be explicit).

### Step 1b.3: Generalize `SEASBdrLoadOperator`

**File:** `solver/seas_bdrload_operator.hpp`

Same pattern:
```cpp
template <typename MeshType = Mesh,
          typename DomainOpType = AntiplaneBdrLoadOperator<MeshType>>
class SEASBdrLoadOperator : public TimeDependentOperator
{
public:
```

### Step 1b.4: Verify compilation

All existing code must compile without changes:
- `bp1_verification_full.cpp` uses `SEASQuasiDynamicOperator<ParMesh>` → default arg fills in `AntiplaneDomainOperator<ParMesh>` ✓
- `bp1_bdrload.cpp` uses `SEASBdrLoadOperator<ParMesh>` → default arg fills in `AntiplaneBdrLoadOperator<ParMesh>` ✓

### Step 1b.5: Create `tests/unit/test_domain_operator_interface.cpp`

Test cases:
1. **`test_num_slip_components_antiplane`**: Verify `NumSlipComponents() == 1` on `AntiplaneDomainOperator`.
2. **`test_get_fault_coords_2d_default`**: Call `GetFaultCoords2D()`, verify `x2` is all zeros, `x3` matches `GetFaultDepths()`.
3. **`test_polymorphic_dispatch`**: Create operator, store in `DomainOperator<Mesh>*` base pointer, verify all interface methods work through base pointer.
4. **`test_off_fault_displacement_default`**: Call `GetOffFaultDisplacement()`, verify empty result.

### Step 1b.6: Update build system

Add new test targets to both `Makefile` and `CMakeLists.txt`:
- `seas_test_fault_detection` (Phase 1a)
- `seas_test_domain_operator_interface` (Phase 1b)

Add `domain/seas_boundary_tags.hpp` to `SEAS_HEADERS`.

---

## Files Modified (Summary)

| File | Change |
|------|--------|
| `domain/seas_boundary_tags.hpp` | **NEW** — SEASBoundaryTags constants + FaultBoundaryData utility |
| `domain/domain_operator.hpp` | Add 3 virtual methods with defaults |
| `domain/antiplane_operator.hpp` | Add `fault_tag` constructor param, dual-path `SetupFaultInfo()` |
| `domain/antiplane_bdrload_operator.hpp` | Same as antiplane_operator.hpp |
| `solver/seas_operator.hpp` | Add `DomainOpType` template param with default |
| `solver/seas_bdrload_operator.hpp` | Add `DomainOpType` template param with default |
| `bp1/mesh/bp1.geo` | Add `Physical Curve(5) = {7,8,9,10,11,12};` |
| `bp1/mesh/bp1_selfsimilar.geo` | Add `Physical Curve(5) = {17,18};` |
| `bp2/mesh/bp2.geo` | Add `Physical Curve(5) = {7,8,9,10,11};` |
| `tests/unit/test_fault_detection.cpp` | **NEW** — tag-based vs coordinate-based tests |
| `tests/unit/test_domain_operator_interface.cpp` | **NEW** — interface method tests |
| `Makefile` + `CMakeLists.txt` | Add new test targets + header |

## Verification Plan

1. Build all existing tests — must compile with zero changes to test source
2. Run `seas_test_fault_detection` — verify tag-based and coordinate-based produce identical results
3. Run `seas_test_domain_operator_interface` — verify new interface methods
4. Run existing `seas_test_antiplane`, `seas_test_fault_operator`, `seas_test_quasi_dynamic` — must pass unchanged
5. Regenerate BP1 mesh from updated `.geo`, run short BP1 smoke test with `fault_tag=5` — verify identical output to `fault_tag=-1`
