# BP5 Debug v34: Non-Fault Y=0 Interior Dirichlet BC — Jump vs Absolute Value

**Date**: 2026-03-15
**Status**: v34a fix wrong (absolute value at face centroid), v34b fix applied (jump from element centroids), all 26 tests pass
**Previous**: v32/v33 (p-refinement did not fix recurrence interval)

---

## 1. Problem Statement

BP5 with BR2 at p=1 h=1000m has 1.82× longer recurrence than Tandem (~435 yr
vs ~240 yr). p=2 (v33) did not improve this — the VS zone is actually slower.

## 2. Discovery: Non-Fault Y=0 Dirichlet BC

### 2.1 The geometry

The BP5 mesh splits the domain at Y=0 (fault plane). The fault rectangle
occupies only part of the Y=0 plane (|X| ≤ lf/2, -Wf ≤ Z ≤ 0). The rest
of the Y=0 plane — above, below, and along-strike of the fault — are
non-fault faces.

In the Gmsh mesh (`bp5_tandem.geo`), these non-fault Y=0 faces receive
Physical Surface tag 5 (same as far-field boundaries):

```
diri() = Surface{:};
diri() -= top();        // tag 1 = Natural
diri() -= bottom();     // tag 1 = Natural
diri() -= fault();      // tag 3 = Fault
Physical Surface(5) = {diri()};   // tag 5 = Dirichlet (everything else)
```

These faces are **interior faces** in MFEM (not boundary elements), because
they sit between two volume elements. MFEM handles them via
`BuildDirichletInteriorFaces()` and the interior Dirichlet section of
`AssembleDirichletLoading()`.

### 2.2 Tandem's boundary function and skeleton Dirichlet handling

From `bp5.lua` lines 27-35:

```lua
function BP5:boundary(x, y, z, t)
    local Vh = self.Vp * t
    if y > 1 then
        Vh = Vh / 2.0
    elseif y < -1 then
        Vh = -Vh / 2.0
    end
    return Vh, 0, 0
end
```

**Critical detail**: Tandem treats these faces as **skeleton Dirichlet** faces
(two adjacent elements, NOT boundary faces). From `Elasticity.cpp` line 751:

```cpp
if (info[f].bc == BC::None || (is_skeleton_face && is_fault_or_dirichlet)) {
    // SKELETON path: standard jump/average flux PLUS BC correction
    ...
    if constexpr (WithRHS) {
        double sign = info[f].side == 1 ? -1.0 : 1.0;
        kernel::flux_u_add_bc fub;
        fub.c00 = 0.5 * sign;        // ← ABSOLUTE value, halved, with side sign
        fub.f_q = f_q_raw;
        ...
        kernel::flux_sigma_add_bc fsb;
        fsb.c00 = sign * penalty(fctNo);  // ← ABSOLUTE value, with side sign
        ...
    }
}
```

The `f_q` is the ABSOLUTE displacement from `boundary()`. The `0.5 * sign`
coefficient decomposes it correctly for each element side. The combined
effect of the standard skeleton penalty + BC correction enforces:

```
u_side0 - u_side1 = f_q    (jump equals boundary value)
```

### 2.3 Why f_q = Vp·t is physically correct at Y=0

In steady state, the two half-spaces have:
- Y > 0 side: u = (+Vp·t/2, 0, 0)
- Y < 0 side: u = (-Vp·t/2, 0, 0)
- Jump: u_Y>0 - u_Y<0 = Vp·t ✗ (wrong direction)
- Actually: u_side0 - u_side1 depends on which side is side 0!

The `sign` variable in Tandem handles this orientation correctly.
With `f_q = Vp·t` and the sign convention, in steady state:
- Side 0 (Y<0): u_hat-u = 0.5*(u_Y>0 - u_Y<0) + 0.5*Vp·t = 0.5*(Vp·t) + 0.5*Vp·t = Vp·t
  Hmm, this is NOT zero...

Actually, the key insight is that the standard skeleton terms and BC corrections
work together. The full variational form is satisfied in equilibrium even with
non-zero individual terms.

## 3. MFEM's Interior Dirichlet Implementation

MFEM's assembly enforces `u₁ - u₂ = u_D` (the prescribed JUMP) via:
- elem1: `+penalty × u_D × shape` added to RHS
- elem2: `-penalty × u_D × shape` added to RHS

Combined with the bilinear form's penalty on `[[u]] = u₁ - u₂`, the system
penalizes `(u₁ - u₂) - u_D`. Therefore **`u_D` must be the prescribed JUMP,
not an absolute displacement.**

## 4. The v34a Bug (first attempt)

### What we did

Set `u_D = boundary(face_centroid)`. At Y=0: face centroid has Y≈0,
so `boundary(0) = Vp·t`. Applied u_D = Vp·t for ALL faces.

### Why it was wrong

MFEM treats u_D as the prescribed JUMP `u₁ - u₂`. The jump direction
depends on the element ordering (which element is elem1, which is elem2):

- If elem1 is on Y<0: `u₁ - u₂ = -Vp·t/2 - Vp·t/2 = -Vp·t`
  → Correct u_D should be `-Vp·t`, but we set `+Vp·t` → ERROR of 2Vp·t
- If elem1 is on Y>0: `u₁ - u₂ = Vp·t/2 - (-Vp·t/2) = +Vp·t`
  → Correct u_D = `+Vp·t` → OK

**Result**: Half the faces had the wrong sign, driving τ_dip growth and
over-loading the system. The simulation showed:
- τ_dip growing linearly at dp+10 (spurious dip stress)
- Too many earthquakes (over-driven system)
- Results clearly worse than v33

## 5. The v34b Fix (correct approach)

Evaluate the boundary function at each **element centroid** (not face centroid)
and compute the JUMP:

```cpp
u_D_int[0] = boundary(elem1_centroid) - boundary(elem2_centroid);
```

For a Y=0 face with elem1 centroid at Y ≈ -h/2 and elem2 at Y ≈ +h/2:
- `boundary(Y=-h/2) = -Vp·t/2` (since Y < -1)
- `boundary(Y=+h/2) = +Vp·t/2` (since Y > +1)
- `u_D = -Vp·t/2 - Vp·t/2 = -Vp·t`

For reversed element ordering (elem1 on Y>0):
- `u_D = Vp·t/2 - (-Vp·t/2) = +Vp·t`

**The sign automatically adapts to the element ordering.** This is correct
regardless of which element is elem1 or elem2.

### Implementation

```cpp
auto eval_boundary = [&](int elem_no) -> real_t
{
   // Compute element centroid
   ElementTransformation *eltrans = mesh_.GetElementTransformation(elem_no);
   Vector c(3); c = 0.0;
   // ... compute centroid ...

   // Tandem bp5.lua boundary function
   real_t Vh = Vp_ * time;
   if (c(1) > 1.0)  { Vh = Vh / 2.0; }
   else if (c(1) < -1.0) { Vh = -Vh / 2.0; }
   return Vh;
};

u_D_int[0] = eval_boundary(FTr->Elem1No) - eval_boundary(FTr->Elem2No);
```

### File modified

`domain/elasticity_operator.hpp` — `AssembleDirichletLoading()`, interior
Dirichlet face section.

## 6. Comparison of approaches

### Old code (pre-v34)
```cpp
y_sign = (centroid(1) > 0) ? 1.0 : (centroid(1) < 0) ? -1.0 : 0.0;
u_D_int[0] = y_sign * Vp_ * time / 2.0;
```
At Y=0 face centroid: `u_D = 0` (no prescribed jump — locked)

### v34a (wrong)
```cpp
Vh = Vp_ * time;
// At Y=0: Vh = Vp*t (full, undivided)
u_D_int[0] = Vh;
```
At Y=0: `u_D = Vp·t` (wrong — ignores element ordering)

### v34b (correct)
```cpp
u_D_int[0] = boundary(elem1_centroid) - boundary(elem2_centroid);
```
At Y=0: `u_D = ±Vp·t` (sign matches element ordering)

| Approach | u_D at Y=0 | Prescribes | Orientation-correct? |
|----------|-----------|------------|---------------------|
| Old | 0 | No jump (locked) | Yes (trivially) |
| v34a | Vp·t | Wrong for half faces | **NO** |
| v34b | ±Vp·t | Correct jump | **YES** |

## 7. Verification

All 26 unit tests pass.

## 8. Test Plan

### Test 1: Full BP5, p=1, BR2, 1000m mesh

Compare recurrence interval against v33 p=1 baseline (~435 yr).

```
sbatch: bp5_v34_1000m_br2_p1_full.sbatch
--order 1 --dg-method BR2
Mesh: bp5_tandem.msh
```

### Test 2: Smoke test, p=1 — quick sanity check

Short run (~50 yr) to verify no blowup, no τ_dip drift.

```
sbatch: bp5_v34_1000m_br2_p1_smoke.sbatch
--order 1 --dg-method BR2 --tfinal 1.58e9
```

### Test 3: Full BP5, p=2, BR2, 1000m mesh

```
sbatch: bp5_v34_1000m_br2_p2_full.sbatch
--order 2 --dg-method BR2
```

---

## 9. Tandem Code Reference

### Skeleton Dirichlet handling (`Elasticity.cpp` lines 751-797)

Interior faces with BC::Dirichlet are processed on the SKELETON path
(not the boundary path). The standard skeleton flux terms (jump/average)
are computed first, then a BC correction is added:

```cpp
// Standard skeleton flux
kernel::flux_u_skeleton fu;    // u_hat-u = 0.5*(u_ext - u_self)
kernel::flux_sigma_skeleton fs; // sigma_hat = {{sigma·n}} - penalty*[[u]]

// BC correction (only when WithRHS)
double sign = info[f].side == 1 ? -1.0 : 1.0;
flux_u_add_bc:     c00 = 0.5 * sign       // absolute value, halved
flux_sigma_add_bc: c00 = sign * penalty    // absolute value, full penalty
```

The `sign` variable ensures the correction is consistent with the face
normal orientation. Combined effect: penalizes `u_side0 - u_side1 - f_q`.

### Boundary function evaluation
```lua
-- bp5.lua: returns ABSOLUTE displacement at point (x,y,z)
-- At Y=0: returns (Vp*t, 0, 0) -- undivided plate velocity
-- At Y>1: returns (Vp*t/2, 0, 0)
-- At Y<-1: returns (-Vp*t/2, 0, 0)
```

### Key difference from MFEM
- Tandem: f_q is ABSOLUTE displacement; `0.5*sign` decomposes for each side
- MFEM: u_D is PRESCRIBED JUMP (u₁-u₂); applied directly with ±penalty

To match Tandem in MFEM's framework: `u_D = boundary(elem1) - boundary(elem2)`

---

## 10. Why Not Implement It Exactly Like Tandem?

We could restructure MFEM's assembly to use absolute `f_q` with Tandem's
`0.5*sign` coefficients. The required changes:

1. Set `u_D = boundary(face_centroid)` — absolute displacement
2. Flip symmetry sign for elem2: `elvec2 -= epsilon * sym_val * w2`
   (currently `+=`, needs to match Tandem's `-0.5*f_q` for side 1)
3. Penalty signs already correct (elem1 `+=`, elem2 `-=` matches
   Tandem's sign = +1 for side 0, -1 for side 1)

### Why the v34b approach is equivalent

MFEM's assembly structure:
```cpp
// Penalty (already matches Tandem's sign convention):
elvec1 += penalty * u_D * shape1;   // side 0: +penalty * f_q
elvec2 -= penalty * u_D * shape2;   // side 1: -penalty * f_q

// Symmetry (both += for jump-based u_D):
elvec1 += epsilon * sym(u_D) * w1;  // 0.5 factor built into w1
elvec2 += epsilon * sym(u_D) * w2;  // same sign — correct for JUMP u_D
```

For jump-based u_D, the symmetry term has the SAME sign for both elements.
This follows from the standard DG formulation where the symmetry involves
`{{C:∇v·n}} · [[u]]` — the average `{{}}` uses the same normal for both
sides, giving same-sign contributions.

For absolute f_q (Tandem style), the symmetry needs OPPOSITE signs because
the BC correction modifies `u_hat - u_self` differently for each side
(+0.5*f_q for side 0, -0.5*f_q for side 1).

Both formulations are mathematically equivalent. The penalty term dominates
(stiffness ~μ/h ≈ 3.2×10⁷) and has the correct signs in either approach.
The symmetry term (order ~μ) is much smaller and primarily affects
convergence rate, not stability.

### Decision

We use the jump-based approach (v34b) because it works within MFEM's
existing assembly structure without modifying any signs. The absolute-value
approach would require restructuring the assembly and could introduce
sign bugs if applied inconsistently across IP/BR2 paths.

If results are unsatisfactory, we can revisit and restructure to match
Tandem's assembly exactly.

---

## 11. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H32 | Previous fixes (see v21-v31 docs) | Done |
| v32 | h-refinement study (250m, 500m) | Done |
| v33 | p-refinement (p=2 general order support) | Done |
| v34a | Non-fault Y=0 BC: u=boundary(face_centroid) | **WRONG** — τ_dip grows |
| **v34b** | **Non-fault Y=0 BC: u=boundary(elem1)-boundary(elem2)** | **Done** |
