# Bug Fix Plan: Traction Extraction — Normal Direction & Quadrature Integration

**Date:** 2026-02-28
**File to modify:** `miniapps/seas/domain/antiplane_operator.hpp`
**Related analysis:** `antiplane_implementation_check_02282026.md` (Issues 2.3 & 2.4)

---

## Context

`ComputeTraction` in `domain/antiplane_operator.hpp` has two related deficiencies
identified in `antiplane_implementation_check_02282026.md` (Issues 2.3 & 2.4):

- **Issue #3 (normal direction):** The traction uses `grad1(0)` — the raw x-component of
  the gradient — instead of dotting with the actual face normal `{μ∇u · n̂}`. This hardcodes
  the assumption that the fault normal is exactly `ê_x`. For the current vertical fault (x=0)
  this is numerically correct only if `CalcOrtho` always returns a `+x` normal, but MFEM's
  `CalcOrtho` may orient face normals in either `±x` direction depending on mesh ordering,
  and the formula fails for any non-axis-aligned fault face.

- **Issue #4 (single midpoint):** The gradient is evaluated at a single face midpoint
  (`ip.x = 0.5`). For polynomial order `p ≥ 2`, the gradient is a degree-`p−1` polynomial
  over the face and cannot be exactly represented by a single point. The slip assembly already
  uses a proper quadrature rule (`IntRules.Get(FTr->FaceGeom, 2*p+1)`) — the traction
  extraction must be consistent with this.

Both issues share a root cause: the traction loop was not written with the same rigour as the
slip assembly. The fix brings `ComputeTraction` up to the same standard.

---

## Correct Formula

### What the document specifies

Algorithm Step 15 (traction extraction):
```
τ_qs = {μ ∇u_h · n̂}|_Γᶠ
```
where `{·}` is the average operator across the fault face.

### Sign convention analysis

`CalcOrtho(FTr->Jacobian(), nor)` returns a normal whose x-component `nor(0)` can be positive
or negative depending on which element MFEM labels as Elem1. For the fault at x=0:

| MFEM orientation | `nor` direction | Effect on `grad1(0)` |
|---|---|---|
| Elem1=left, Elem2=right | `nor(0) > 0` → `n̂ = +ê_x` | `{∂u/∂x}` (correct sign) |
| Elem1=right, Elem2=left | `nor(0) < 0` → `n̂ = -ê_x` | `-{∂u/∂x}` (wrong sign) |

The slip assembly already handles this via:
```cpp
real_t slip_imposed = (nor(0) > 0) ? -slip_phys : slip_phys;
```
The traction fix mirrors this convention by multiplying `{∇u · n̂}` by `sgn(nor(0))`.

**Combined formula** (both fixes together):
```
tau_face = μ * sgn(nor(0)) * (1/|e|) * ∫_e { ∇u_h · n̂ } ds
```

For a vertical fault with `n̂ = ±ê_x`:
```
sgn(nor(0)) * {∇u · n̂} = sgn(nor(0)) * {∂u/∂x} * sgn(nor(0)) = {∂u/∂x}
```
This always equals `{∂u/∂x}` regardless of mesh orientation. For `p=1`, the gradient is
constant and the face average equals the midpoint value — **results are numerically identical
to the current implementation for all existing BP1 tests.**

---

## Two Code Paths

`ComputeTraction` has two paths that both need to be fixed:

1. **Serial interior face path** (lines ~875–915): uses `GradientGridFunctionCoefficient`
2. **Parallel shared face path** (lines ~929–1011): uses manual `dshape * Jinv * u_dofs`

---

## Implementation Plan

### Serial interior face path

**Before (current code):**
```cpp
IntegrationPoint ip;
ip.x = 0.5;
FTr->SetAllIntPoints(&ip);
const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();
Vector grad1(mesh_.Dimension()), grad2(mesh_.Dimension());
grad_u.Eval(grad1, *FTr->Elem1, eip1);
grad_u.Eval(grad2, *FTr->Elem2, eip2);
real_t avg_dudx = 0.5 * (grad1(0) + grad2(0));
real_t tau_face = mu_ * avg_dudx;
```

**After (fixed code):**
```cpp
const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
const FiniteElement *fe2 = fes_->GetFE(FTr->Elem2No);
int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*face_order + 1);

int dim = mesh_.Dimension();
real_t flux_sum = 0.0;
real_t face_len = 0.0;
real_t traction_sign = 0.0;   // set on first quadrature point

for (int q = 0; q < ir.GetNPoints(); q++)
{
    const IntegrationPoint &ip = ir.IntPoint(q);
    FTr->SetAllIntPoints(&ip);
    const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
    const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

    // Face normal (magnitude = face measure ds contribution)
    Vector nor(dim);
    CalcOrtho(FTr->Jacobian(), nor);
    real_t nor_len = nor.Norml2();

    // Consistent sign: same convention as AssembleSlipContributionIP
    if (q == 0) { traction_sign = (nor(0) >= 0.0) ? 1.0 : -1.0; }

    // Unit normal
    Vector n_hat(nor);
    n_hat /= nor_len;

    // Gradients from both elements via GradientGridFunctionCoefficient
    Vector g1(dim), g2(dim);
    grad_u.Eval(g1, *FTr->Elem1, eip1);
    grad_u.Eval(g2, *FTr->Elem2, eip2);

    // Average normal flux: {∇u · n̂}
    real_t avg_flux = 0.5 * ((g1 * n_hat) + (g2 * n_hat));

    // Accumulate face-integrated flux (ds = nor_len * ip.weight)
    flux_sum += ip.weight * nor_len * avg_flux;
    face_len  += ip.weight * nor_len;
}

real_t tau_face = mu_ * traction_sign * flux_sum / face_len;
```

Note: `GradientGridFunctionCoefficient grad_u(&displacement)` is created once before the
outer face loop, not per face (same as current code).

---

### Parallel shared face path

The parallel path computes gradients manually because Elem2 is a face-neighbor whose DOFs
come from `pfes->GetFaceNbrElementVDofs`. DOF values are fetched once per face; shape
derivatives are recomputed at each quadrature point.

**Before (current code — schematic):**
```cpp
// Single midpoint
IntegrationPoint ip; ip.x = 0.5;
FTr->SetAllIntPoints(&ip);
// ... compute grad1 and grad2 at midpoint ...
real_t avg_dudx = 0.5 * (grad1(0) + grad2(0));
real_t tau_face = mu_ * avg_dudx;
```

**After (fixed code):**
```cpp
// --- Fetch DOF values once, outside quadrature loop ---
const FiniteElement *fe1 = fes_->GetFE(FTr->Elem1No);
int nbr_idx = FTr->Elem2No - mesh_.GetNE();
const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);

Array<int> dofs1;
fes_->GetElementDofs(FTr->Elem1No, dofs1);
Vector u1(fe1->GetDof());
par_u.GetSubVector(dofs1, u1);

Array<int> dofs2;
pfes->GetFaceNbrElementVDofs(nbr_idx, dofs2);
const Vector &nbr_data = par_u.FaceNbrData();
Vector u2(fe2->GetDof());
for (int k = 0; k < fe2->GetDof(); k++) { u2(k) = nbr_data(dofs2[k]); }

// --- Quadrature loop ---
int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*face_order + 1);

real_t flux_sum = 0.0;
real_t face_len = 0.0;
real_t traction_sign = 0.0;

for (int q = 0; q < ir.GetNPoints(); q++)
{
    const IntegrationPoint &ip = ir.IntPoint(q);
    FTr->SetAllIntPoints(&ip);
    const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
    const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

    // Face normal
    Vector nor(dim);
    CalcOrtho(FTr->Jacobian(), nor);
    real_t nor_len = nor.Norml2();
    if (q == 0) { traction_sign = (nor(0) >= 0.0) ? 1.0 : -1.0; }
    Vector n_hat(nor);
    n_hat /= nor_len;

    // Elem1 gradient: dshape_ref * Jinv * u1
    DenseMatrix dshape1(fe1->GetDof(), dim);
    fe1->CalcDShape(eip1, dshape1);
    DenseMatrix Jinv1(dim);
    CalcInverse(FTr->Elem1->Jacobian(), Jinv1);
    DenseMatrix dshape1_phys(fe1->GetDof(), dim);
    Mult(dshape1, Jinv1, dshape1_phys);
    Vector grad1(dim);
    grad1 = 0.0;
    for (int k = 0; k < fe1->GetDof(); k++)
        for (int d = 0; d < dim; d++)
            grad1(d) += dshape1_phys(k, d) * u1(k);

    // Elem2 gradient: dshape_ref * Jinv * u2
    DenseMatrix dshape2(fe2->GetDof(), dim);
    fe2->CalcDShape(eip2, dshape2);
    DenseMatrix Jinv2(dim);
    CalcInverse(FTr->Elem2->Jacobian(), Jinv2);
    DenseMatrix dshape2_phys(fe2->GetDof(), dim);
    Mult(dshape2, Jinv2, dshape2_phys);
    Vector grad2(dim);
    grad2 = 0.0;
    for (int k = 0; k < fe2->GetDof(); k++)
        for (int d = 0; d < dim; d++)
            grad2(d) += dshape2_phys(k, d) * u2(k);

    // Average normal flux and accumulate
    real_t avg_flux = 0.5 * ((grad1 * n_hat) + (grad2 * n_hat));
    flux_sum += ip.weight * nor_len * avg_flux;
    face_len  += ip.weight * nor_len;
}

real_t tau_face = mu_ * traction_sign * flux_sum / face_len;
```

---

## Key Design Decisions

| Decision | Choice | Reason |
|---|---|---|
| Sign correction | `sgn(nor(0))` captured at first quadrature point | Matches slip assembly convention; stable even if `CalcOrtho` flips orientation |
| Integration formula | Face-averaged: `flux_sum / face_len` | Each DG face = one fault "patch" with uniform δ; representative traction is the face average |
| Quadrature order | `2*face_order + 1` | Exact for `p−1` polynomial gradient; matches slip assembly |
| Normal computation | `CalcOrtho(FTr->Jacobian(), nor)` inside loop | `FTr->Jacobian()` is updated by `SetAllIntPoints`; must call per quadrature point |
| Serial gradient eval | `GradientGridFunctionCoefficient` in loop | Idiomatic MFEM; handles element-to-physical mapping automatically |
| Parallel gradient | DOF fetch outside loop, `dshape*Jinv*u` inside loop | DOFs don't change with quadrature point; shape derivatives do |

---

## Backward Compatibility

For `order=1` on a uniform straight mesh (all current BP1 tests):
- The quadrature rule has a single Gauss point ≈ face midpoint
- The gradient is constant (degree 0 polynomial) within each element
- `face_len = |e|` (face length) and `flux_sum = |e| * avg_flux`
- Therefore `tau_face = mu_ * traction_sign * avg_flux`, which equals the current value

**Results are numerically identical to the pre-fix code for all p=1 tests.**

---

## Verification Steps

1. **Existing unit tests** — `tests/unit/test_antiplane.cpp`, `debug_ip_consistency.cpp`,
   `test_fault_operator.cpp`: run without regression.

2. **MMS convergence test** (`tests/parallel/mms_antiplane_parallel.cpp`): convergence
   rate should be unchanged for `order=1,2,3`.

3. **BP1 mesh convergence** (`miniapps/seas/bp1/mesh_convergence.py`): traction magnitudes
   and slip rates for the 100m, 50m, 25m meshes should be identical (within floating-point)
   to pre-fix results for `order=1`.

4. **Sign check**: for a uniform right-lateral slip δ > 0 on all fault faces, verify
   `traction(i) > 0` for all i regardless of mesh ordering.

5. **Higher-order accuracy check**: for `order=2`, run a single quasi-static solve with
   known analytical solution and verify that the traction error converges at order ≥ 1.
