# BP5 Debug v34: Non-Fault Y=0 Dirichlet BC Mismatch

**Date**: 2026-03-15
**Status**: Fix implemented, all 32 tests pass, sbatch files ready
**Previous**: v32/v33 (p-refinement did not fix recurrence interval)

---

## 1. Problem Statement

BP5 with BR2 at p=1 h=1000m has 1.82× longer recurrence than Tandem (~435 yr
vs ~240 yr). p=2 (v33) did not improve this — the VS zone is actually slower.

## 2. Discovery: Non-Fault Y=0 Dirichlet BC Value

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

### 2.2 Tandem's boundary function

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

This function is evaluated at ALL Dirichlet faces (tag 5). At Y=0, neither
`y > 1` nor `y < -1` triggers, so `Vh = Vp * t` (undivided):

| Location | Y coord | Tandem u_D |
|----------|---------|-----------|
| Far-field Y > 0 | +100 km | (Vp·t/2, 0, 0) |
| Far-field Y < 0 | -100 km | (-Vp·t/2, 0, 0) |
| **Non-fault Y=0** | **0** | **(Vp·t, 0, 0)** |

The non-fault Y=0 value represents a **locked fault**: both sides of the
plane move together at the full plate velocity. No relative slip.

### 2.3 MFEM's code (before fix)

From `elasticity_operator.hpp` lines 1983-1988:

```cpp
real_t y_sign = (centroid(1) > 0.0) ? 1.0
              : (centroid(1) < 0.0) ? -1.0 : 0.0;
u_D_int[0] = y_sign * Vp_ * time / 2.0;
```

At Y=0, `y_sign = 0`, so `u_D = (0, 0, 0)`.

### 2.4 The discrepancy

| | Tandem | MFEM (before fix) |
|---|--------|-------------------|
| Non-fault Y=0 | **u = (Vp·t, 0, 0)** | **u = (0, 0, 0)** |

MFEM was pinning the non-fault Y=0 plane to zero displacement. Tandem
prescribes the full plate velocity there.

## 3. Physical Impact

The non-fault Y=0 faces surround the fault rectangle on three sides:
- Along-strike: |X| > lf/2 (beyond the fault length)
- Depth: Z < -Wf (below the fault)
- Surface: Z > 0 (above the fault, if any)

Pinning these faces to u=0 creates artificial resistance to plate loading:
- The far-field drives displacement at ±Vp·t/2, but the non-fault Y=0
  boundary holds at zero → stress concentration at the fault edges
- This artificial constraint reduces the effective loading rate on the
  VS zone, slowing the stress accumulation needed for earthquake nucleation
- The effect is strongest at the VS/VW transition and the fault base,
  exactly where the recurrence interval is controlled

## 4. Fix

Replace the `y_sign` formula with Tandem's exact `boundary()` logic:

```cpp
real_t Vh = Vp_ * time;
if (centroid(1) > 1.0)
{
   Vh = Vh / 2.0;
}
else if (centroid(1) < -1.0)
{
   Vh = -Vh / 2.0;
}
u_D_int[0] = Vh;
```

At Y=0 (all non-fault interior faces): `u_D = (Vp·t, 0, 0)`.

This matches Tandem's `bp5.lua` boundary function exactly. The 1-meter
threshold between Y>1 and Y<-1 matches Tandem's Lua code.

### File modified

`domain/elasticity_operator.hpp` — `AssembleDirichletLoading()`, interior
Dirichlet face section (lines ~1983-1999).

## 5. Verification

All 32 unit tests pass. The far-field boundary faces (actual boundary
elements) are unaffected — they use a separate code path at line 1739
with `sign = (centroid(1) > 0.0) ? 1.0 : -1.0` (no Y=0 case needed
since boundary faces are at Y=±100 km).

## 6. Test Plan

### Test 1: Full BP5, p=1, BR2, 1000m mesh

Compare recurrence interval against v33 p=1 baseline (~435 yr).
If the BC fix improves the loading rate, the recurrence should decrease
toward Tandem's ~240 yr.

```
sbatch: bp5_v34_1000m_br2_p1_full.sbatch
--order 1 --dg-method BR2
Mesh: bp5_tandem.msh (same as v33)
```

### Test 2: Full BP5, p=2, BR2, 1000m mesh

Same fix at p=2. Compare against v33 p=2 result.

```
sbatch: bp5_v34_1000m_br2_p2_full.sbatch
--order 2 --dg-method BR2
Mesh: bp5_tandem.msh (same as v33)
```

### Test 3: Smoke test, p=1, BR2 — quick sanity check

Short run (~50 yr) to verify no blowup and correct initial behavior.

```
sbatch: bp5_v34_1000m_br2_p1_smoke.sbatch
--order 1 --dg-method BR2 --tfinal 1.58e9
```

---

## 7. Tandem Code Reference

### BC enum (`src/form/BC.h`)
```cpp
enum class BC : int { None = 0, Natural = 1, Fault = 3, Dirichlet = 5 };
```

### Tag → BC mapping (`src/io/GlobalSimplexMeshBuilder.cpp` lines 79-97)
```cpp
switch (tag) {
case static_cast<long>(BC::Dirichlet): bc = BC::Dirichlet; break;  // tag 5
case static_cast<long>(BC::Fault):     bc = BC::Fault;     break;  // tag 3
case static_cast<long>(BC::Natural):   bc = BC::Natural;   break;  // tag 1
}
```

### Boundary enforcement (`app/localoperator/Elasticity.cpp` line 448)
```cpp
bool Elasticity::assemble_boundary(...) {
    if (info.bc == BC::Natural) { return false; }  // traction-free
    // Dirichlet/Fault: assemble penalty + consistency + symmetry
}
```

### Dirichlet value source (`app/form/SeasQDOperator.cpp` line 57)
```cpp
dgop_->set_dirichlet((*fun_boundary_)(time));  // calls bp5.lua boundary()
```

---

## 8. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H32 | Previous fixes (see v21-v31 docs) | Done |
| v32 | h-refinement study (250m, 500m) | Done |
| v33 | p-refinement (p=2 general order support) | Done |
| **v34** | **Non-fault Y=0 Dirichlet BC: u=0 → u=Vp·t (match Tandem)** | **Done** |
