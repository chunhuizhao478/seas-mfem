# BP5 Debug v6: Revert Wrong BR2 Normal Fix + Deep Analysis

## Summary

The Phase 2 "BR2 normal convention fix" applied in v5 is **dimensionally wrong** and must be reverted. It changed the BR2 traction correction in `ComputeTraction()` from using `basis.normal` (unit normal) to `nor` (unnormalized normal) with `/ nor_sq` division, producing units of [Pa/m^2] instead of [Pa]. This effectively zeroes the BR2 penalty at typical mesh scales (h=1000m), removing DG stabilization and causing V_max to grow unboundedly.

## Frontera v5 Results (100 ranks, 1000m Gmsh mesh)

```
- 40650 RHS BLOWUP messages
- 82732 TRACTION BLOWUP messages
- 0 DISPLACEMENT BLOWUP messages
- 0 zeroIn errors
- V_max: 0.03 -> 50 (step 143) -> 1178 m/s (step 1257)
- Affected ranks: 6, 18, 25, 26, 27, 29, 30, 34
```

V_max grows without saturation. This is fundamentally different from the v4 behavior (localized zeroIn errors at specific ranks).

## Root Cause Analysis

### The BR2 Traction Correction Formula

The DG traction at a fault face is:

```
T = {sigma(u)} . n_hat - eta * (C : r_e([[u]] - delta)) . n_hat
```

Where:
- `{sigma(u)}` = average stress from both sides (consistency term)
- `eta * (C : r_e(...)) . n_hat` = BR2 penalty correction
- `n_hat` = unit outward normal
- `r_e` = BR2 lifting operator

The BR2 lifting `r_e` maps a jump function to a tensor field via:

```
integral_Omega r_e : tau dOmega = integral_e {psi x n_hat} dS
```

### Discrete Implementation

The discrete lifting at the face centroid:

```
face_int[us, m] = shape_m * jump[u] * nor[s]     (nor = n_hat * |J_face|)
f_lifted = 0.5 * face_int * Minv^T
eval = sum_m shape_m * f_lifted[us, m]
```

The `face_int` correctly uses `nor` (unnormalized, from CalcOrtho) because this is part of the face integral, encoding the surface area element.

The final traction contraction uses:

```
T_BR2[i] = eta * sum_{u,s} C_{ius}(n) * eval[us]
```

Where `C_{ius}(n)` must use the **unit normal** `n_hat` because traction = force per unit area evaluated at a point.

### Dimensional Analysis

**Original code (CORRECT):**
- `face_int` ~ [m] x [m^2] = [m^3]  (jump [m] x unnorm-normal [m^2])
- `Minv` ~ [m^-3]
- `f_lifted` ~ dimensionless
- `eval` ~ dimensionless
- `tn = C * basis.normal` ~ [Pa] x [1] = [Pa]  (C is [Pa], unit normal is dimensionless)
- `correction = sum` ~ **[Pa]**  CORRECT

**Phase 2 fix (WRONG):**
- `tn = C * nor` ~ [Pa] x [m^2] = [Pa.m^2]  (nor has units [m^2])
- `sum` ~ [Pa.m^2]
- `nor_sq = nor . nor` ~ [m^4]
- `correction = sum / nor_sq` ~ [Pa.m^2] / [m^4] = **[Pa/m^2]**  WRONG

For h=1000m mesh: `nor_sq ~ h^4 = 10^12`. The correction is divided by 10^12, making it ~10^6x too small. The BR2 penalty is effectively zero, removing DG stabilization entirely.

## Deep Comparison with Tandem

### Tandem's Traction Formula (Simpler, Dimensionally Inconsistent for BR2)

Tandem uses the SAME formula for both IP and BR2 traction:

```python
# app/kernels/elasticity.py:242-244
traction_q = 0.5 * (sigma_1.n_hat + sigma_2.n_hat)
           + c0 * (E_q[0]*u[0] - E_q[1]*u[1] - f_q)
```

Where:
- `c0 = -penalty(fctNo)`
- For IP: `penalty = kappa * {mu/h}` → has dimensions [Pa/m]
- For BR2: `penalty = NumFacets` (= 4 for tets, dimensionless)
- `f_q` = prescribed slip in physical 3D coordinates (units [m])
- `E_q[x]*u[x]` = displacement at face quad points (units [m])

For BR2, this gives: `T = [Pa] - 4 * [m]` — **dimensionally inconsistent!**

However, this works in practice because `[[u]] - delta ≈ 0` in a well-converged DG solution, so the penalty term is negligible regardless of its dimensions.

### Our Implementation is MORE Correct

Our implementation computes the FULL BR2 lifting for the traction correction:
```
T -= br2_penalty * 0.5 * C : r_e([[u]] - delta) . n_hat
```

This involves mass matrix inverse, shape function evaluation, and elasticity tensor coupling — giving proper [Pa] dimensions. Both approaches produce similar results since [[u]] - δ ≈ 0.

### Key Tandem Details

- **`f_q` computation** (`app/localoperator/ElasticityAdapter.cpp:37-49`): The slip is transformed from 2D fault coordinates to 3D physical coordinates using the fault basis vectors. No BR2 lifting or mass inverse is applied to the slip before it enters the traction formula.
- **Unit normal** (`app/localoperator/Elasticity.cpp:982`): `krnl.n_unit_q = fct[fctNo].get<UnitNormal>().data()->data()`
- **Penalty values** match: Tandem NumFacets=4 (tets) corresponds to our `br2_penalty = dim+1 = 4` (tets) or `2*dim = 6` (hex)

## Deep Analysis of Our Implementation

### Consistency Traction Term: CORRECT

Both interior and shared faces use identical structure:
1. `dshape_phys = dshape_ref * Jinv` (physical gradients via Jacobian inverse)
2. `avg_grad = 0.5 * (grad1 + grad2)` (average gradient)
3. `strain = 0.5 * (avg_grad + avg_grad^T)` (symmetric strain)
4. `stress = lambda * tr(eps) * I + 2*mu*eps` (isotropic elasticity)
5. `T = stress * basis.normal` (traction with unit normal)

No bugs found. Interior and shared faces are structurally identical.

### BR2 Quadrature: CONSISTENT for BP5

| Aspect | Traction | RHS | Matrix |
|--------|----------|-----|--------|
| Quad points | 1 (centroid) | Full (nqp) | Full (nqp) |
| Weights in face_int | None | YES (w_q) | YES (w_q) |
| Normal | Unnormalized | Unnormalized | Unnormalized |

For p=1 on flat fault faces (BP5 uses x=0 plane), single-point centroid evaluation is equivalent to full quadrature because:
- Reference face weight for 1-point rule = 1.0
- Shapes are linear → centroid value = average value on flat face
- Normal is constant on flat faces

This is NOT a bug for BP5.

### Data Exchange in ComputeTraction: CORRECT

`ExchangeFaceNbrData()` is called at lines 2023-2024:
```cpp
pfes->ExchangeFaceNbrData();       // FES exchange
par_u.ExchangeFaceNbrData();       // Solution vector exchange
```
Both called BEFORE the shared face loop. Neighbor DOF values accessed via `par_u.FaceNbrData()` at line 2058.

### BR2 Correction Magnitude (After Revert)

For [[u]] - delta = 0.001 m (1mm DG error), h = 1000m:
- face_int ~ 0.001 * h^2 = 10^3
- Minv ~ 1/h^3 = 10^-9
- f_lifted ~ 10^-6
- tn ~ mu = 32 GPa
- correction ~ 6 * 0.25 * 32e9 * 10^-6 ≈ 48 kPa

This is small compared to the initial shear stress (~15 MPa), confirming the correction is well-behaved when the solver produces a good solution.

## What Was Ruled Out

1. **`scalar_fes_` being serial** — only used for `GetFE(Elem1No)` (local elements)
2. **Mass inverse indexing** — `elem_mass_inv_[FTr->Elem2No]` correct (`ne + i` for shared faces)
3. **Missing ExchangeFaceNbrData in ComputeTraction** — confirmed called before shared face loop
4. **Double-counting shared faces** — interior and shared are disjoint in MFEM
5. **Sign convention** — matches interior face pattern and antiplane reference
6. **MUMPS INFOG(1)=-9** — auto-retry succeeds, not root cause
7. **RHS magnitude ~1e15** — expected for ~30m accumulated slip
8. **Consistency traction term** — structurally identical for interior and shared faces
9. **Quadrature weight mismatch** — irrelevant for p=1 on flat faces
10. **BR2 penalty value** — matches Tandem (number of faces per element)

## Fix

### Revert Phase 2 in `ComputeTraction()` (elasticity_operator.hpp)

#### Interior faces (~lines 1974-1999):

```cpp
// CURRENT (WRONG):
real_t nor_sq = nor * nor;
// ...
real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * nor(i)
   + mu_val_ * ((i == u ? 1.0 : 0.0) * nor(s)
               + (i == s ? 1.0 : 0.0) * nor(u));
// ...
correction[i] = br2_penalty * 0.5 * sum / nor_sq;

// REVERT TO (CORRECT):
real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * basis.normal[i]
   + mu_val_ * ((i == u ? 1.0 : 0.0) * basis.normal[s]
               + (i == s ? 1.0 : 0.0) * basis.normal[u]);
// ...
correction[i] = br2_penalty * 0.5 * sum;
```

#### Shared faces (~lines 2203-2224):

Same change: replace `nor(ci)` with `basis.normal[ci]`, remove `nor_sq` line and `/ nor_sq` division.

`basis.normal` is available in both sections via `fault_basis_.GetBasis()` (line 1860 for interior, line 2121 for shared).

### Keep Phase 1 Diagnostics

All Phase 1 diagnostics (traction bounds, solver convergence, RHS assembly) remain in place. They work correctly and will help diagnose the separate v4 bug (zeroIn at scale).

### Update Comments

Remove misleading comments about "unnormalized normal for consistency with BR2 matrix/RHS" (lines 1974-1975 and 2203-2204). The traction is a point evaluation requiring the unit normal, which is different from the volume integral convention used in the matrix/RHS.

## Implementation Results

### Changes Applied

1. **BR2 normal revert** — Reverted `nor(i)` → `basis.normal[i]` and removed `/ nor_sq` in both interior and shared face sections of `ComputeTraction()`.

2. **Solver: MUMPS → CG+BoomerAMG** — Replaced MUMPS direct solver with `CGSolver` + `HypreBoomerAMG` preconditioner in `AssembleStiffness()`. This matches Tandem's approach (iterative solver for 3D). MUMPS was causing hangs on macOS due to slow factorization, even for small meshes.

3. **Smoke test: `GetGlobalNE()` deadlock fix** — `ParMesh::GetGlobalNE()` is a **collective MPI call** (internally uses `MPI_Allreduce`). It was called only on rank 0 inside `if (mpi.IsRoot())`, causing rank 1+ to never participate → MPI deadlock. Fixed by calling on all ranks before the root-only print.

4. **Smoke test: mesh resolution fix** — The original BP5-sized mesh (4×4×1 = 16 hex over 400×200×100 km) had face centroids at z=50km > Wf=40km, so zero fault DOFs were detected. Changed to small mesh (4×2×1 = 8 hex, Lx=4, Ly=2, Lz=2) with Wf=Lz, lf=2*Ly to reliably capture fault faces.

### Test Results

| Test | Result |
|------|--------|
| Serial elasticity (`seas_test_elasticity_operator`) | **77/77 passed** |
| Parallel 2-rank (`seas_test_parallel_elasticity`) | **6/6 passed** |
| Parallel 4-rank (`seas_test_parallel_elasticity`) | **6/6 passed** |
| BP5 smoke np=1 (`seas_test_bp5_parallel_smoke`) | **13/13 passed** |
| BP5 smoke np=2 (`seas_test_bp5_parallel_smoke`) | **13/13 passed** |

## Lessons Learned

### 1. `ParMesh::GetGlobalNE()` is Collective

**Rule**: NEVER call `GetGlobalNE()` inside `if (mpi.IsRoot())`. It uses `MPI_Allreduce` internally and all ranks must participate.

**Pattern** (correct):
```cpp
long long global_ne = pmesh.GetGlobalNE();  // ALL ranks call
if (mpi.IsRoot()) {
   std::cout << "Elements: " << global_ne << "\n";  // Only root prints
}
```

**Anti-pattern** (deadlock):
```cpp
if (mpi.IsRoot()) {
   std::cout << "Elements: " << pmesh.GetGlobalNE() << "\n";  // DEADLOCK!
}
```

This same pattern applies to any MFEM method that may internally call MPI collectives. When in doubt, call the method on all ranks and only guard the I/O.

### 2. Tandem Uses Iterative Solvers for 3D

Tandem's BP5 config uses:
- **Matrix-free mode** (`matrix_free = true`)
- **p-Multigrid** with logarithmic strategy
- **PETSc CG/FGMRES** (not MUMPS direct solve)

MUMPS is only used for 2D problems (`lu_mumps.cfg`). For 3D DG elasticity, MUMPS is too slow even for tiny meshes (16 hex elements hang on macOS). `CG + HypreBoomerAMG` works well as a drop-in replacement.

### 3. Unit Normal vs Unnormalized Normal in DG

The BR2 traction computation has two conceptually distinct parts:

| Part | What it does | Which normal to use |
|------|-------------|-------------------|
| `face_int` (lifting) | Face integral ∫ ψ⊗n dS | `nor` (unnormalized, includes |J_face|) |
| `tn` (elasticity coupling) | Point evaluation T = C:ε·n̂ | `basis.normal` (unit normal) |

Using `nor` for the elasticity tensor coupling produces units of [Pa·m²] instead of [Pa], which when divided by `nor_sq = ||nor||² ~ h⁴` gives [Pa/m²] — wrong by a factor of ~h² (10⁶ for h=1000m).

### 4. Smoke Test Mesh Must Match Fault Bounds

When using `params.Wf` and `params.lf` for fault detection, the mesh element size must place face centroids within the fault zone:
- Face centroid z < Wf (fault depth limit)
- Face centroid |y| < lf/2 (fault length limit)

For BP5 (`Wf=40km, lf=100km`), this requires element height < 80km in z and element width < 100km in y. A simpler approach is to use a unit-scale mesh with `Wf=Lz, lf=2*Ly`.

## Remaining Issue: v4 Bug (Separate from Phase 2)

After reverting Phase 2, the code returns to v4 behavior. The v4 bug (zeroIn errors at 56+ ranks) is a separate issue:

- v4 on Frontera: `zeroIn: F(a)=2.16e10` at Ranks 25, 34
- Only manifests at scale (56+ ranks), passes locally with 2-4 ranks
- tau_abs ~ 21.6 GPa is physically impossible (mu=32 GPa would require strain ~0.67)

Possible v4 root causes (to investigate after revert):
1. **Solver accuracy at scale** — CG+AMG may produce different convergence behavior than MUMPS at high rank counts → large [[u]] - delta → large BR2 correction → traction blowup
2. Mesh partitioning edge cases creating problematic face configurations
3. Something else that only appears with many ranks

The Phase 1 diagnostics will identify which subsystem (RHS, displacement, or traction) blows up first in the v4 case.
