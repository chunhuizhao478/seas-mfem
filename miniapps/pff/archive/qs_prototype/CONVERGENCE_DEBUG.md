# Phase Field Fracture Convergence Issues: Debug Analysis

## Problem Description
The MFEM PFF solver experiences convergence difficulties when `d_max` approaches 0.6. This document analyzes potential bugs and formulation differences compared to MOOSE/Raccoon.

---

## 1. Diffusion Coefficient Mismatch (CRITICAL BUG)

### MFEM Implementation (`pff_solver.hpp`, line 1387-1388):
```cpp
// Diffusion coefficient: Gc * l / c0
real_t diff_coeff = mat_.Gc * mat_.l / mat_.c0;
```

### Raccoon Implementation (`ADPFFDiffusion.C`, line 40):
```cpp
return 2 * _Gc[_qp] * _l[_qp] / _c0[_qp] * value;
```

### Analysis:
The AT2 phase-field energy functional is:
```
Γ = ∫ [Gc/c0 * (α(d)/l + l|∇d|²)] dΩ
```
where α(d) = d² for AT2.

Taking the functional derivative with respect to d:
```
δΓ/δd = Gc/c0 * [dα/dd / l - 2l Δd]
      = Gc/c0 * [2d/l - 2l Δd]
      = 2Gc/(c0*l) * d - 2Gc*l/c0 * Δd
```

In weak form (testing with v):
```
∫ 2Gc*l/c0 * ∇d·∇v dΩ + ∫ 2Gc/(c0*l) * d * v dΩ
```

**Conclusion:**
- **Raccoon's diffusion coefficient: `2*Gc*l/c0`** ✓
- **MFEM's diffusion coefficient: `Gc*l/c0`** ✗ (MISSING FACTOR OF 2!)

This 2x error in diffusion affects:
- The balance between diffusion and reaction terms
- The effective regularization length
- Can cause oscillations and non-physical solutions at higher damage values

---

## 2. Damage Source Term: Linear vs. Nonlinear Treatment

### MFEM Implementation (`pff_solver.hpp`, lines 1400-1405):
```cpp
// Source term: -g'(d)*H = 2*(1-d)*(1-η)*H
DamageSourceCoefficient src_coeff(d_gf_.get(),
                                  &damage_op_->GetStrainEnergyHistory(),
                                  mat_.eta);
b.AddDomainIntegrator(new DomainLFIntegrator(src_coeff));
```

The source is treated as a **right-hand side** that depends on the current d.

### Raccoon Implementation (`ADPFFSource.C`, line 32):
```cpp
return _dpsi_dd[_qp];
```

Where `psi = alpha*Gc/c0/l + g*psie_active`, so:
```
dpsi/dd = dalpha/dd * Gc/(c0*l) + dg/dd * psie_active
        = 2d * Gc/(c0*l) - 2(1-d)*(1-η) * psie_active
```

### Analysis:
MFEM's formulation separates the equation as:
```
A*d = B
where A = diffusion + reaction (linear in d)
      B = source (depends on current d)
```

This is a **Picard/fixed-point linearization**. The problem is:
1. The source `-g'(d)*H` is evaluated at the **current** d, not the **new** d
2. This creates a lag that can cause oscillations when d changes rapidly
3. Near d=0.6, the nonlinearity is strong (g'(d) varies significantly)

Raccoon uses **full Newton iteration** where dpsi/dd is:
```
2d*Gc/(c0*l) - 2(1-d)*(1-η)*H
```
This is linearized consistently at each Newton step with the Jacobian:
```
d(dpsi/dd)/dd = 2*Gc/(c0*l) + 2*(1-η)*H
```

**Potential Bug:** MFEM's approach doesn't account for the Jacobian contribution from the source term.

---

## 3. VI Solver: Simplified vs. Full Complementarity

### MFEM Implementation (`pff_solver.hpp`, lines 1459-1483):
```cpp
// Solve for new damage: A * d_new = B
// (This is a direct solve, not incremental)
...
cg.Mult(B, d_new);

// Project onto bounds [d_lower, d_upper]
for (int i = 0; i < d_new.Size(); i++)
{
   d_new(i) = std::max(d_lower(i), std::min(d_upper(i), d_new(i)));
}
```

### Raccoon Implementation:
Uses PETSc's `vinewtonrsls` (Variational Inequality Newton with Reduced Linear System):
```
petsc_options_iname = '-snes_type'
petsc_options_value = 'vinewtonrsls'
```

### Analysis:
MFEM uses a **projected iteration** approach:
1. Solve unconstrained system
2. Project onto bounds

This is simple but has issues:
- **No active set identification**: The solver doesn't know which DOFs are at bounds until after projection
- **Modified system not assembled**: For DOFs at bounds, the system should be modified (row/column elimination)
- **Oscillation risk**: Solution can oscillate between projection and solve

The `vinewtonrsls` solver properly:
1. Identifies active set (DOFs at bounds)
2. Eliminates active DOFs from the reduced linear system
3. Only solves for free DOFs
4. Updates active set iteratively

---

## 4. Nonlinear vs. Linear Damage Solve

### MFEM Implementation:
Damage equation is assembled as a **linear system** A*d = B, where:
- A is constant (diffusion + reaction)
- B depends on current d

This is a **quasi-linear** approach that requires outer iterations.

### Raccoon Implementation:
Uses full Newton on the nonlinear residual:
```
R(d) = ∫ 2Gc*l/c0 * ∇d·∇v dΩ + ∫ dpsi/dd * v dΩ = 0
```

The Jacobian is:
```
J = ∫ 2Gc*l/c0 * ∇φ·∇v dΩ + ∫ d²psi/dd² * φ * v dΩ
```

where:
```
d²psi/dd² = d²alpha/dd² * Gc/(c0*l) + d²g/dd² * psie_active
          = 2 * Gc/(c0*l) + 2*(1-η) * psie_active
```

**Potential Bug:** MFEM's Jacobian is missing the `2*(1-η)*H` contribution from the degradation function.

---

## 5. Convergence Check Location

### MFEM Implementation (`pff_solver.hpp`, lines 1450-1457):
```cpp
// Check convergence
if (res_norm < vi_abs_tol ||
    (res_norm_0 > 0 && res_norm / res_norm_0 < vi_rel_tol))
{
   break;
}
```

Convergence is checked **before** the linear solve, based on the previous iteration's residual.

### Issue:
When damage is near the transition region (d ≈ 0.5-0.7):
- The residual might appear small after complementarity zeroing
- But the actual solution change could still be significant
- Early termination can leave the solution in an inconsistent state

---

## 6. Staggered Iteration Instability at High Damage

### Issue Description:
At d ≈ 0.6, the coupling between elasticity and damage becomes highly nonlinear:

1. **Degradation function sensitivity**:
   - g(0.6) = 0.4² * (1-η) + η ≈ 0.16
   - g'(0.6) = -2 * 0.4 * (1-η) ≈ -0.8

2. **Strain energy feedback**:
   - High d → low stiffness → higher strain → higher H
   - Higher H → drives more damage → potential runaway

3. **Fixed-point convergence radius**:
   - The spectral radius of the staggered iteration increases with d
   - May exceed 1 near critical damage levels

### Raccoon's Approach:
- Uses tighter tolerances: `nl_rel_tol = 1e-8`, `nl_abs_tol = 1e-10`
- `automatic_scaling = true` helps condition the system
- Full Newton captures the correct sensitivities

---

## 7. Summary of Potential Bugs

| # | Issue | Severity | Location |
|---|-------|----------|----------|
| 1 | Diffusion coefficient missing factor of 2 | **CRITICAL** | `pff_solver.hpp:1387-1388` |
| 2 | Source term Jacobian not included | HIGH | `pff_solver.hpp:1400-1405` |
| 3 | Simplified VI projection instead of reduced system | MEDIUM | `pff_solver.hpp:1476-1480` |
| 4 | Linear damage solve instead of Newton | MEDIUM | `pff_solver.hpp:1459-1474` |
| 5 | Convergence check before solve | LOW | `pff_solver.hpp:1450-1457` |

---

## 8. Recommended Fixes (Priority Order)

1. **Fix diffusion coefficient**: Change `Gc*l/c0` to `2*Gc*l/c0`

2. **Implement full Newton for damage**: Include the Jacobian contribution from the source term:
   ```cpp
   // Full Jacobian: diffusion + reaction + source derivative
   // d(dpsi/dd)/dd = 2*Gc/(c0*l) + 2*(1-η)*H
   // Add mass integrator with coefficient 2*(1-η)*H
   ```

3. **Proper VI active set method**: Use reduced system approach that eliminates active DOFs from the linear system

4. **Add relaxation/line search**: For staggered iteration stability at high damage:
   ```cpp
   d_new = (1-omega)*d_old + omega*d_solved;
   ```
   where omega < 1 for damping

---

## 9. Verification Steps

1. Run with artificially low Gc to delay fracture, verify elastic response
2. Compare damage profile width with analytical: w ≈ 4*l for AT2
3. Check energy dissipation: ∫ Gc*(alpha/l + l|∇d|²)/c0 dΩ ≈ Gc * crack_length
4. Compare force-displacement curve with Raccoon reference

---

## References

- Raccoon source: `/Users/chunhuizhao/projects/farms_cdms/raccoon/`
- MOOSE framework: https://mooseframework.org/
- Bourdin et al., "Numerical experiments in revisited brittle fracture" (2000)
- Miehe et al., "A phase field model for rate-independent crack propagation" (2010)
