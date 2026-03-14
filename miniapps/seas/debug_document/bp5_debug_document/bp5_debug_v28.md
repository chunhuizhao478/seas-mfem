# BP5 Debug v28: IP Bilinear Form Penalty Mismatch — |nor|² vs |nor|

**Date**: 2026-03-14
**Status**: Bilinear form discrepancy identified, custom integrator needed
**Previous**: v27 (traction scalar penalty, still blows up)

---

## 1. The Discovery

The IP bilinear form penalty in MFEM's `DGElasticityIntegrator` uses a
**different geometric scaling** from Tandem's IP formulation.

### MFEM (bilininteg.cpp line 4161)

```cpp
const real_t jmatcoef = kappa * (nor*nor) * wLM;
//                              ^^^^^^^^
//                              |nor|² = face_area²
```

The penalty per quadrature point:
```
kappa * |nor_q|² * (ip.weight/2) * (λ+2μ) * (1/detJ₁ + 1/detJ₂)
```

### Tandem (elasticity.py lines 98-100, 123)

```python
# surfaceOp apply (line 100):
c2[x] * w['q'] * E_q[x]['kq'] * u_jump['pq'] * nl_q['q']

# lift_ip (line 123):
L_q[x]['lpuq'] <= E_q[x]['lq'] * delta['pu'] * nl_q['q']
```

where `nl_q = |nor_q|` = NormalLength (DGCurvilinearCommon.cpp line 84:
`length[i] = norm(normal[i])`)

The penalty per quadrature point:
```
penalty * ip.weight * |nor_q|
```

### The difference

| | MFEM | Tandem |
|---|------|--------|
| Normal scaling | **|nor|²** (squared) | **|nor|** (linear) |
| Penalty coefficient | `kappa * (λ+2μ)/detJ` | `(D+1)*c_N_1*(A/V)*(c₁²/c₀)/4` |
| Material coupling | `(λ+2μ)` only | `c₁²/c₀ = (3λ+2μ)²/(2μ)` |
| Structure | Same-component only (diagonal) | Same-component only (diagonal) |

The `|nor|²` vs `|nor|` difference means the penalty has a different
dependence on face geometry. For a face with area A:
- MFEM scales as A²
- Tandem scales as A × (material-dependent scalar)

This produces different effective penalty magnitudes, different DG solutions,
different displacement jumps, and different tractions.

---

## 2. Impact

The MFEM IP penalty being differently scaled from Tandem's means:

1. The DG stiffness matrix is different → DG solution is different
2. The displacement jumps [[u]] at fault faces are different
3. The traction (which depends on the DG solution) is different
4. Even with Tandem-matching traction formula, the traction values differ
   because the SOLUTION is different

This is why fixing only the traction (v24-v27) never worked — the
underlying DG solution was already wrong because of the bilinear form
penalty mismatch.

---

## 3. Tandem's Full IP Formulation

### Bilinear form assembly (assembleSurface, line 141-145)

```python
a[x][y]['kplu'] <=
    c0[x] * E_q[x]['kq'] * w['q'] * traction_op_q[y]['lupq'] +    # consistency
    c1[y] * E_q[y]['lq'] * w['q'] * traction_op_q[x]['kpuq'] +    # symmetry
    c2[abs(y-x)] * w['q'] * E_q[x]['kq'] * L_q[y]['lpuq']         # penalty
```

### Coefficients (from Elasticity.cpp)

Interior face:
```
c0[0] = -0.5,  c0[1] = +0.5
c1[0] = ε*0.5 = -0.5,  c1[1] = +0.5
c2[0] = +penalty,  c2[1] = -penalty
```

Boundary face:
```
c0[0] = -1.0
c1[0] = ε = -1.0
c2[0] = +penalty
```

### Penalty (Elasticity.h line 131-136, Elasticity.cpp line 278-290)

```cpp
double penalty(std::size_t fctNo) const {
    if (method_ == DGMethod::BR2) {
        return NumFacets;  // = D+1 = 4
    }
    return penalty_[fctNo];  // IP: material-dependent
}

// IP penalty:
auto const p = [&](int side) {
    const auto [c0, c1] = stiffness_tensor_bounds(info.up[side]);
    constexpr double c_N_1 = InverseInequality<Dim>::trace_constant(PolynomialDegree - 1);
    return (Dim + 1) * c_N_1 * (area_[fctNo] / volume_[info.up[side]]) * (c1 * c1 / c0);
};
penalty_[fctNo] = (p(0) + p(1)) / 4.0;  // interior
// or: penalty_[fctNo] = p(0);            // boundary
```

Material bounds:
```
c0 = 2μ                    (min eigenvalue of C)
c1 = Dim*λ + 2μ = 3λ+2μ   (max eigenvalue of C)
c1²/c0 = (3λ+2μ)²/(2μ)
```

### Consistency/symmetry term (traction_op_q, line 137-139)

```python
traction_op_q[x]['kpuq'] <=
    lam_q[x] * Dx_q[x]['kpq'] * n_q['uq'] +
    mu_q[x] * (n_q['jq'] * Dx_q[x]['kjq'] * delta['pu'] + Dx_q[x]['kuq'] * n_q['pq'])
```

This is: `[C : ∇φ_k e_p]_u · n = λ (∂φ_k/∂x_p) n_u + μ(∂φ_k/∂x_j n_j δ_{pu} + ∂φ_k/∂x_u n_p)`

Uses the **unnormalized** normal `n_q` (from face Jacobian, magnitude ≈ face area).
Gradients `Dx_q` are in **physical** coordinates (using J⁻¹).
The division by detJ is embedded in the gradient computation.

### IP penalty lift (lift_ip, line 122-123)

```python
L_q[x]['lpuq'] <= E_q[x]['lq'] * delta['pu'] * nl_q['q']
```

This is: `L_q_{lpu}(q) = φ_l(q) * δ_{pu} * |nor_q|`

- Same-component only (δ_{pu})
- Scaled by `|nor_q|` = NormalLength (linear, NOT squared)
- The basis function φ_l is in reference space

### Traction (compute_traction, line 242-244)

```python
traction_q['pq'] <= 0.5 * (traction(0, n_unit_q) + traction(1, n_unit_q)) +
                    c0[0] * (E_q[0] * u[0] - E_q[1] * u[1] - f_q)
```

- Average stress uses **unit** normal `n_unit_q`
- Penalty: `c0 = -penalty` (scalar) × displacement jump
- For BR2: penalty = 4 (negligible)
- For IP: penalty = material-dependent (significant)

---

## 4. Custom IP Integrator Design

Replace MFEM's `DGElasticityIntegrator` with a custom integrator that
matches Tandem's formula exactly.

### Class: `DGElasticityTandemIPIntegrator`

```cpp
class DGElasticityTandemIPIntegrator : public BilinearFormIntegrator
{
public:
    DGElasticityTandemIPIntegrator(Coefficient &lambda, Coefficient &mu,
                                    real_t epsilon, int dim = 3);

    void AssembleFaceMatrix(const FiniteElement &el1,
                            const FiniteElement &el2,
                            FaceElementTransformations &Trans,
                            DenseMatrix &elmat) override;

private:
    Coefficient &lambda_, &mu_;
    real_t epsilon_;  // -1 for SIPG
    int dim_;

    // Compute Tandem-style penalty for a face
    real_t ComputePenalty(FaceElementTransformations &Trans,
                          real_t lam, real_t mu) const;
};
```

### Penalty computation

```cpp
real_t ComputePenalty(FaceElementTransformations &Trans,
                      real_t lam, real_t mu) const
{
    real_t c0 = 2.0 * mu;                    // min eigenvalue of C
    real_t c1 = dim_ * lam + 2.0 * mu;       // max eigenvalue of C

    // Inverse inequality trace constant for p=0
    // (PolynomialDegree - 1 = 0 for order-1 DG)
    real_t c_N_1 = 1.0;

    // Face area (from unnormalized normal at centroid)
    const IntegrationPoint &ip = Geometries.GetCenter(Trans.GetGeometryType());
    Trans.SetAllIntPoints(&ip);
    Vector nor(dim_);
    CalcOrtho(Trans.Jacobian(), nor);
    real_t face_area = nor.Norml2();

    // Element volumes
    real_t vol1 = Trans.Elem1->Weight();
    real_t vol2 = (Trans.Elem2No >= 0) ? Trans.Elem2->Weight() : vol1;

    real_t p0 = (dim_ + 1) * c_N_1 * (face_area / vol1) * (c1 * c1 / c0);
    real_t p1 = (dim_ + 1) * c_N_1 * (face_area / vol2) * (c1 * c1 / c0);

    if (Trans.Elem2No >= 0) {
        return (p0 + p1) / 4.0;  // interior face
    } else {
        return p0;                // boundary face
    }
}
```

### Assembly structure

For each quadrature point q on the face:

**Consistency** (same as MFEM, uses unnormalized nor and adjugate gradients):
```
elmat contribution:
  -0.5 * [[φ_k e_i]] · {C : ∇(φ_l e_u)} · n_q * w_q    (for interior)
  -1.0 * φ_k e_i · (C : ∇(φ_l e_u)) · n_q * w_q         (for boundary)
```

**Symmetry** (same as MFEM, with ε = -1):
```
  ε * 0.5 * {C : ∇(φ_k e_i)} · n_q * [[φ_l e_u]] * w_q  (for interior)
  ε * 1.0 * (C : ∇(φ_k e_i)) · n_q * φ_l e_u * w_q       (for boundary)
```

**Penalty** (Tandem-style: penalty × |nor| × δ_{iu} × shape × shape):
```
  +penalty * |nor_q| * w_q * φ_k * φ_l * δ_{iu}    (for interior, same side)
  -penalty * |nor_q| * w_q * φ_k * φ_l * δ_{iu}    (for interior, cross side)
  +penalty * |nor_q| * w_q * φ_k * φ_l * δ_{iu}    (for boundary)
```

Note: `|nor_q|` is the NormalLength (linear), NOT `|nor_q|²` (squared).

---

## 5. Files to Create/Modify

| File | Action |
|------|--------|
| `integrator/dg_elasticity_tandem_ip_integrator.hpp` | **NEW**: Custom IP integrator matching Tandem |
| `domain/elasticity_operator.hpp` | Use new integrator for IP method |
| `domain/elasticity_operator.hpp` | IP traction: Tandem scalar penalty (already done in H28) |

---

## 6. Slip RHS for IP

The slip RHS also needs to match Tandem. Currently MFEM's IP slip path
(`AssembleSlipContributionIP`) uses the same `kappa * |nor|²` scaling.
This must be updated to use Tandem's scalar penalty formula.

The Tandem IP RHS formula (from elasticity.py line 174-178):
```
b['kp'] += c1 * tractionTest(0, f_q) * w['q'] +      # consistency
           c2 * w['q'] * E_q['kq'] * f_lifted_q['pq']  # penalty
```

where for IP: `f_lifted_q = f_q * nl_q` (line 156: `rhs_lift_ip`)

So the slip RHS penalty is: `penalty * w_q * φ_k * slip_p * |nor_q|`

This matches the bilinear form penalty structure (penalty × |nor| × shape × data).

---

## 7. Dirichlet Loading RHS for IP

Similarly, the Dirichlet loading RHS needs the same penalty formula.
Currently `AssembleDirichletLoading` IP path uses `kappa * |nor|²`.

---

## 8. Summary: What Needs to Match

| Component | Tandem formula | MFEM current | Match? |
|-----------|---------------|--------------|--------|
| IP bilinear consistency | `0.5 * [[u]]·{C:∇v}·n` | Same | ✓ |
| IP bilinear symmetry | `ε*0.5 * {C:∇u}·n·[[v]]` | Same | ✓ |
| **IP bilinear penalty** | **penalty × |nor| × [[u]]·[[v]]** | **kappa × |nor|² × (L+2M)/detJ × [[u]]·[[v]]** | **✗** |
| **IP slip RHS penalty** | **penalty × |nor| × slip·[[v]]** | **kappa × |nor|² × ... × slip·[[v]]** | **✗** |
| **IP Dirichlet RHS penalty** | **penalty × |nor| × g^D·v** | **kappa × |nor|² × ... × g^D·v** | **✗** |
| IP traction | penalty × jump (scalar) | penalty × jump (scalar, H28) | ✓ |

Three penalty terms need fixing: bilinear form, slip RHS, and Dirichlet RHS.
All need to use Tandem's `penalty * |nor_q|` instead of MFEM's `kappa * |nor_q|²`.

---

## 9. Implementation Complete

### New files
- `integrator/dg_elasticity_ip_penalty_integrator.hpp`: Custom IP penalty integrator
  using `penalty * |nor| * [[u]]·[[v]]` (isotropic scalar penalty)

### Modified files
- `domain/elasticity_operator.hpp`:
  - IP bilinear form: split into `DGElasticityIntegrator(kappa=0)` +
    `DGElasticityIPPenaltyIntegrator` (both interior and boundary faces)
  - IP slip RHS: `wq_penalty = penalty_ip * ip.weight * nl_q`
    (interior faces at line ~809, shared faces at line ~1239)
  - IP Dirichlet RHS: `wq_penalty = p0 * ip.weight * nl_q`
    (boundary faces at line ~1588)
  - IP traction: `correction = penalty_ip * jump` (scalar × vector)
    (interior at line ~2207, shared at line ~2415)

### Unit tests: 675 pass, 0 fail
5 new tests added to `test_elasticity_br2.cpp`:
- TestIPPenaltySymmetry: face matrix symmetric on hex mesh ✓
- TestIPPenaltyPSD: positive semi-definite on tet mesh ✓
- TestIPPenaltyValue: non-zero with positive diagonal ✓
- TestSplitIPConsistency: split IP same sparsity as full IP ✓
- TestIPPenaltyTet: full assembly on tets symmetric and positive ✓

---

## 10. Test Plan

### Test 1: Uniform fault (no earthquake) — stability check

```
sbatch: bp5_v28_test1_ip_uniform.sbatch
--dg-method IP --delta-tau-factor 0 --V-nuc 1e-9
```

PASS: V/Vp in [0.5, 2.0] at z=22km for 300yr, no blowup, no VS lockup.

### Test 2: Tandem initialization — cycling check

```
sbatch: bp5_v28_test2_ip_tandem.sbatch
--dg-method IP --psi-init-mode tandem --V-nuc 0.01
```

PASS: First event ~150yr, VS maintains V≈Vp, subsequent cycling.

### Test 3: SCEC initialization — full benchmark

```
sbatch: bp5_v28_test3_ip_scec.sbatch
--dg-method IP (defaults: V_nuc=0.03, SCEC psi, delta_tau=1)
```

PASS: Immediate earthquake, VS recovers, cycling, matches Tandem.

---

## 11. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H25 | Previous fixes | Done |
| H26 | Traction: scalar penalty for BR2 | Done |
| H27 | Mesh: uniform fault resolution | Done |
| H28 | IP traction: scalar penalty (isotropic) | Done |
| H29 | IP bilinear form: custom penalty integrator (`penalty * |nor|`) | Done |
| H30 | IP slip RHS: same penalty formula | Done |
| H31 | IP Dirichlet RHS: same penalty formula | Done |
| — | Unit tests: 675 pass, 0 fail | Done |
