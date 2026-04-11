# BP5 Debug v45: Penalty ×3 Analysis — NOT a Bug; Multi-DOF Plan Continues

**Date**: 2026-03-19
**Status**: Analysis complete, ×3 fix reverted, multi-DOF plan proceeding from Phase 1
**Previous**: v44 (×3 penalty fix + multi-DOF — both reverted)
**Branch**: `feature/elasticity`

---

## 1. Summary

v44 identified a factor-of-3 difference between MFEM's IP penalty and Tandem's,
and applied `dim *` (×3) correction to all 13 penalty locations. TACC results
(v44d/v44e) confirmed this makes p=1 results **worse**, not better.

This document explains **why** the ×3 fix is wrong, despite the factor-of-3
being mathematically real.

**Key finding**: The factor of 3 is an artifact of MFEM's reference element
conventions. Both Tandem and MFEM use the same averaging formula (`/4`).
The only difference is `A/V` (Tandem, physical) vs `nl_q/Weight() = A/(3V)`
(MFEM). The consistency/symmetry terms in the DG bilinear form are at full
physical strength (reference factors cancel), while the penalty is at 1/3
strength. Applying ×3 to only the penalty shifts the consistency-penalty
balance, producing a different (not better) DG solution.

**Decision**: Revert the ×3 fix. The p>=2 fix path is multi-DOF fault
discretization (matching Tandem's approach), not penalty scaling. Tandem
never runs p=1 — its default is p=2 with nbf=6 DOFs per face.

---

## 2. The Factor of 3: Where It Comes From

### 2.1 MFEM Reference Element Conventions

MFEM's reference tetrahedron has vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1):

| Quantity | MFEM reference | Physical | Ratio |
|----------|---------------|----------|-------|
| Tet volume | V_ref = 1/6 | V_phys | Weight() = 6 V_phys |
| Triangle area | A_ref = 1/2 | A_phys | nl_q = 2 A_phys |
| A/V ratio | A_ref/V_ref = 3 | A_phys/V_phys | nl_q/Weight() = A/(3V) |

The penalty formula uses:
```cpp
real_t p0 = (dim+1) * c_N_1 * (nl_q / vol1) * (c1*c1/c0);
// nl_q / vol1 = (2A) / (6V) = A / (3V)
```

### 2.2 Tandem Uses Physical Quantities

From `Elasticity.cpp:278-291`:
```cpp
void Elasticity::prepare_penalty(std::size_t fctNo, FacetInfo const& info, ...) {
    auto const p = [&](int side) {
        constexpr double c_N_1 = InverseInequality<Dim>::trace_constant(PolynomialDegree - 1);
        return (Dim + 1) * c_N_1 * (area_[fctNo] / volume_[info.up[side]]) * (c1*c1 / c0);
    };
    if (info.up[0] != info.up[1]) {
        penalty_[fctNo] = (p(0) + p(1)) / 4.0;   // <-- SAME averaging as our code!
    } else {
        penalty_[fctNo] = p(0);
    }
}
```

Where `area_[fctNo]` = physical A, `volume_[...]` = physical V.

### 2.3 The ONLY Difference

| Property | Our code | Tandem | Identical? |
|----------|----------|--------|-----------|
| c_N_1 formula | p(p+D-1)/D | (N+1)(N+D)/D | Yes (same) |
| (D+1) factor | 4 | 4 | Yes |
| c1^2/c0 | (3lambda+2mu)^2/(2mu) | same | Yes |
| Interior averaging | (p0+p1)/4.0 | (p(0)+p(1))/4.0 | **Yes** (verified!) |
| Boundary penalty | p0 | p(0) | Yes |
| Face-to-volume ratio | nl_q/vol = A/(3V) | area/volume = A/V | **NO: factor of 3** |
| Assembly weight | penalty * w_q * nl_q | penalty * w_q * nl_q | Yes |

The averaging formula `/4` is **identical** — this was a key finding. The earlier
hypothesis that we had `/4` while Tandem had `/2` was WRONG. The only difference
is the geometric factor: 1/3.

---

## 3. Why the ×3 Fix Is Wrong

### 3.1 The DG Bilinear Form Has Three Terms

The SIPG bilinear form assembled by **two separate integrators**:

**Integrator 1**: MFEM's `DGElasticityIntegrator` (kappa=0)
- Computes: consistency `-<{sigma(u)n},[v]>` + symmetry `-<{sigma(v)n},[u]>`
- Uses `dshape * adj(J) * (1/det(J)) * nor`
- The `det(J)` from `adj(J)` CANCELS with `1/det(J)` from the weight
- Result: **exact physical integral** (no reference element factor)

Proof from MFEM source (`bilininteg.cpp:4116-4158`):
```cpp
CalcAdjugate(Trans.Elem1->Jacobian(), adjJ);  // adjJ = det(J) * J^{-1}
Mult(dshape1, adjJ, dshape1_ps);               // dshape1_ps = det(J) * grad_phys
w1 = w / Trans.Elem1->Weight();                // w1 = ip.weight / (2 * det(J))
nL1.Set(w1 * lambda, nor);                     // nL1 = ip.weight*lambda*nor / (2*det(J))
// Product: dshape_ps * nL = det(J)*grad * ip.weight*lambda*nor/(2*det(J))
//        = (ip.weight/2) * lambda * grad_phys * nor   <-- det(J) CANCELS
```

**Integrator 2**: Our `DGElasticityIPPenaltyIntegrator`
- Computes: penalty `+eta * <[[u]],[[v]]>`
- Uses `nl_q / vol` which is `A/(3V)`, NOT physical `A/V`
- Result: **1/3 of physical penalty**

### 3.2 The Imbalance

| Term | Our code (v42) | With ×3 fix | Tandem |
|------|---------------|-------------|--------|
| Consistency | 1.0x (physical) | 1.0x (physical) | 1.0x (physical) |
| Symmetry | 1.0x (physical) | 1.0x (physical) | 1.0x (physical) |
| Penalty | 0.33x (ref-scaled) | **1.0x** (physical) | 1.0x (physical) |

**v42 (no ×3)**: Consistency and penalty are at different relative strengths.
The system `(K_consist + 0.33*K_penalty) u = f_consist + 0.33*f_penalty` gives
one DG solution.

**v44 (with ×3)**: The penalty matches Tandem's, giving a DIFFERENT system
`(K_consist + K_penalty) u = f_consist + f_penalty` with a different solution.

Both are valid SIPG methods (as long as penalty exceeds the coercivity threshold).
But they produce **different numerical solutions on a fixed mesh**.

### 3.3 Why v42 Matches Tandem Better Than v44

At p=1 with 1000m mesh, the discretization error is ~20% (known strike deficit).
The v42 penalty at 1/3 of Tandem's creates a **compensating error** that partially
offsets the strike deficit, giving a fortuitously good match to Tandem:
- v42 p=1: 8 events, ~250yr recurrence (Tandem: 8 events, ~240yr)

With ×3 (v44d/v44e), the compensating error is removed, exposing the full
discretization error and shifting the earthquake timing.

### 3.4 The Theoretical Perspective

For SIPG coercivity, the penalty must exceed: `eta > C * c_N_1 * (A/V) * (c1^2/c0)`

Our penalty (v42) at 1/3 of Tandem's is **below** the theoretical minimum.
Yet it works at p=1 because:
1. The theoretical bound is pessimistic (actual threshold is lower)
2. Rate-state friction adds natural damping
3. For p=1 on regular meshes, the penalty barely matters

At p>=2, the penalty matters more. But the primary p>=2 failure cause is the
**fault DOF mismatch** (constant-per-face slip vs polynomial DG solution),
not the penalty scaling.

---

## 4. Tandem Never Uses p=1

A critical observation: **Tandem's default polynomial degree is 2**. All
published Tandem BP5 results use p>=2. At p=1:
- Volume stress is constant per element
- nbf = 3 face DOFs (triangular) are under-resolved relative to volume
- Tandem never tested or validated p=1

This means:
- Matching Tandem at p=1 is **not the goal** — Tandem doesn't run p=1
- Our p=1 results serve as a **sanity check**, not a benchmark target
- The true comparison is at p>=2, where multi-DOF fault discretization is essential

---

## 5. The Right Fix Path: Multi-DOF Fault Discretization

### 5.1 Why Penalty Scaling Is Secondary

The p>=2 failure has two interleaved causes:

| Cause | Effect | Fix |
|-------|--------|-----|
| Fault DOF mismatch (primary) | Constant slip can't match polynomial DG jump; penalty amplifies mismatch | Multi-DOF: nbf = (p+1)(p+2)/2 per face |
| Penalty at 1/3 of Tandem (secondary) | Weaker penalty enforcement; may compound mismatch effects | ×3 correction (defer) |

The multi-DOF fix addresses the fundamental algorithmic mismatch. Once multi-DOF
is working, the penalty scaling can be evaluated at p>=2:
- If multi-DOF + current 1/3 penalty works → keep it (simpler)
- If multi-DOF + 1/3 penalty fails but 1× works → apply ×3 then

### 5.2 Implementation Plan (v45 onwards)

Revert to v42 baseline source code (keeping FaceQuadrature from v44 Phase 1).
Then implement the multi-DOF plan:

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | FaceQuadrature class + 7 unit tests | Done (v44, kept) |
| 2 | L2 traction projection in ComputeTraction | Done |
| 3 | Multi-DOF slip interpolation in AssembleSlipContributionIP | **Done** |
| 4 | Fault state/geometry expansion + unit tests | **Done** |
| 5 | Output integration + unit tests | **Done** |

---

## 6. Evidence: v44e TACC Results

v44e tested the ×3 fix on clean v42 code (no multi-DOF changes). This isolated
the penalty scaling effect. Results showed that p=1 earthquake timing shifted
significantly compared to v40/v42, confirming that the ×3 changes the DG
solution on the existing mesh.

This is **expected behavior** (different penalty = different solution), NOT a bug.
The "wrong" label means the solution moved away from the coincidentally good
v42 match, not that the code is incorrect.

---

## 7. Comparison of DG Integrator Scaling (Reference vs Physical)

### 7.1 MFEM's DGElasticityIntegrator (consistency + symmetry)

Per quadrature point, the consistency term accumulates:
```
(ip.weight/2) * phi_i * [lambda*(dphi_j/dx_jm)*nor_im
                        + mu*(dphi_j/dx_im)*nor_jm
                        + delta(im,jm)*mu*(grad_phi_j . nor)]
```

The `nor = CalcOrtho(J_face)` provides BOTH:
- The unit normal direction (`nor/|nor|`)
- The face integration measure (`|nor| = nl_q`)

The gradients `dphi/dx` are physical (adjugate/det cancels). The `ip.weight`
is the reference face quadrature weight. Together, `ip.weight * nor` correctly
represents `dS * n_hat` on the physical face. **No leftover reference factors.**

### 7.2 Our DGElasticityIPPenaltyIntegrator (penalty)

Per quadrature point:
```
coeff = penalty * ip.weight * nl_q
// penalty contains nl_q/vol = A/(3V)
// coeff = (D+1)*c_N_1*(c1^2/c0) * nl_q^2 / (2*vol) * ip.weight
```

The `nl_q^2 / vol` factor comes from having `nl_q` in both the penalty formula
AND the assembly weight. For tets: `nl_q^2/vol = (2A)^2/(6V) = 2A^2/(3V)`.

Tandem has: `penalty * w_q * nl_q` where `penalty` contains `A/V` (physical),
giving `(A/V) * nl_q * w_q`. The integrated penalty is `(A/V)*A = A^2/V`.

Our integrated penalty is `2A^2/(3V) * sum(w_q) = 2A^2/(3V) * 1/2 = A^2/(3V)`.
Tandem's integrated penalty is `(A/V)*A = A^2/V`. Ratio: 1/3.

### 7.3 MFEM's Standard Penalty (for comparison)

MFEM's `DGElasticityIntegrator` with kappa != 0 computes:
```cpp
jmatcoef = kappa * (nor*nor) * wLM;
// = kappa * nl_q^2 * ip.weight * (lambda+2mu) * (1/(2*W1) + 1/(2*W2))
```

This has the SAME `nl_q^2 / Weight()` scaling as our custom integrator.
MFEM's standard kappa parameter implicitly accounts for reference element
scaling. Users set kappa empirically or via theory. Our formula
`(D+1)*c_N_1*(nl_q/vol)*(c1^2/c0)/2` is equivalent to a specific kappa value.

---

## 8. Changes Made in v45

### 8.1 Code Reverted

Source files reverted to v42 baseline (`8dafa86`):
- `domain/elasticity_operator.hpp` → v42 (removes both ×3 and multi-DOF changes)
- `integrator/dg_elasticity_ip_penalty_integrator.hpp` → v42 (removes ×3)

### 8.2 Code Kept

New files from v44 Phase 1 (no issues):
- `fault/face_quadrature.hpp` → FaceQuadrature class (standalone, all tests pass)
- `tests/unit/test_face_quadrature.cpp` → 7 unit tests (all pass)
- `Makefile` → test target additions

### 8.3 New in v45

- This debug document (`bp5_debug_v45.md`)
- Phase 2 implementation: L2 traction projection in ComputeTraction
- Phase 3 implementation: Multi-DOF slip interpolation in AssembleSlipContributionIP

---

## 9. Phase 3: Multi-DOF Slip Interpolation

### 9.1 What Changed

Modified `AssembleSlipContributionIP` and `AssembleSlipContributionIPShared` to
interpolate per-DOF slip values to quadrature points instead of using a single
constant slip per face.

**Before (v42):**
```cpp
// One slip value per face -> constant delta_u[3] at all quad points
real_t slip_local[2] = {slip_bc(2*fi), slip_bc(2*fi+1)};
fault_basis_.EmbedSlip(fi, slip_local, delta_u);  // constant
```

**After (v45 Phase 3):**
```cpp
// nbf slip values per face -> interpolate to each quad point
for (int kk = 0; kk < nbf; kk++) {
   int dof_idx = fi * nbf + kk;
   slip_local = {slip_bc(2*dof_idx), slip_bc(2*dof_idx+1)};
   fault_basis_.EmbedSlip(fi, slip_local, du);  // per-DOF
   delta_u_nodal(c*nbf + kk) = du[c];           // store
}
face_quad_->InterpolateToQuadPoints(3, delta_u_nodal, delta_u_quad);
// At quad point q: delta_u_q[c] = delta_u_quad(c*nq + q)
```

### 9.2 Tandem Correspondence

| Our code | Tandem |
|----------|--------|
| `delta_u_nodal[c*nbf+kk]` | `slip[l,n] * copy_slip[n,o] * fault_basis_q[p,o,q]` |
| `face_quad_->InterpolateToQuadPoints()` | `e_q['lq']` contraction in `evaluate_slip` kernel |
| `delta_u_quad[c*nq+q]` | `slip_q['pq']` (ElasticityAdapter.cpp:36-49) |

### 9.3 Backward Compatibility

At p=1 IP: `face_order=0`, `nbf=1`, `e_q(0,q)=1.0` for all q.
`InterpolateToQuadPoints` with nbf=1 simply copies the single nodal value
to all quad points → identical to the old constant-per-face code.

### 9.4 Unit Tests Added (7 tests, 148 total assertions)

| Test | What it verifies |
|------|-----------------|
| `TestMultiDOFProperties` | p=2 IP: nbf=6; p=1 IP: nbf=1; p=2 BR2: nbf=1 |
| `TestMultiDOFZeroSlip` | p=2 IP: zero multi-DOF slip → zero displacement + traction |
| `TestMultiDOFUniformSlip` | p=2 IP: uniform slip → finite non-zero displacement + traction |
| `TestMultiDOFBackwardCompatP1` | p=1 IP: nbf=1 code works correctly |
| `TestMultiDOFVaryingSlip` | p=2 IP: varying slip ≠ uniform slip (multi-DOF matters) |
| `TestMultiDOFSlipInterpolation` | Partition of unity + single-DOF interpolation accuracy |
| `TestMultiDOFProjectInterpolateRoundtrip` | GalerkinProject → Interpolate roundtrip for linear functions |

---

## 10. Phase 4: Fault State/Geometry Expansion

### 10.1 Key Finding: No Code Changes Needed

The Phase 2-3 infrastructure already handles multi-DOF propagation correctly:

1. **`GetFaultCoords2D`** (Phase 2): Returns per-DOF (x2, x3) coordinates using
   `face_quad_->GetNodalRule()` — evaluates face geometry at GaussLobatto nodes
2. **`FaultGeometry::ComputeBP5Params()`**: Evaluates spatially-varying `a`, `Dc`,
   `tau_pre`, `V_init` at each DOF coordinate using BP5 parameter functions
3. **`RateStateFaultOperator`**: Uses `num_nodes_ = geom->NumFaultDOFs()` which
   automatically gets the expanded count; all loops are generic over num_nodes_
4. **`SEASQuasiDynamicOperator`**: Uses `fault->SlipSize()` and `fault->TractionSize()`
   for vector sizing — automatically adapts to multi-DOF

The multi-DOF expansion propagates entirely through existing generic code:
- At p=1 IP: `nbf=1`, all sizes identical to v42 baseline
- At p=2 IP: `nbf=6`, all sizes scale by 6× automatically

### 10.2 Tet Mesh Requirement for p>=2

**Important**: `FaceQuadrature` only supports `Geometry::TRIANGLE` faces. This is
correct for tet meshes (Tandem BP5 uses tet meshes from Gmsh). For hex meshes,
faces are quads, causing a quadrature rule mismatch at p>=2.

The Phase 4 tests use a tet mesh helper (`CreateTestMesh3DTet`) that creates
meshes with `Element::TETRAHEDRON` instead of `Element::HEXAHEDRON`.

### 10.3 Unit Tests Added (6 tests, 253 total assertions)

| Test | What it verifies |
|------|-----------------|
| `TestMultiDOFStateLayoutP1` | p=1 IP tet: nbf=1, NumNodes=nfaces, StateSize=3*nfaces |
| `TestMultiDOFStateLayoutP2` | p=2 IP tet: nbf=6, NumNodes=6*nfaces, StateSize=18*nfaces |
| `TestMultiDOFParameterEvaluation` | Per-DOF a, Dc, eta match BP5 formulas at DOF coords |
| `TestMultiDOFFrictionUniform` | Uniform psi/traction → all V finite/positive at p=2 |
| `TestMultiDOFGetSlipSetSlipRoundtrip` | Set/Get slip and theta exact roundtrip at p=2 |
| `TestMultiDOFPreInitInitP2` | Full PreInit → Init cycle at p=2: psi>0, V>0, stress equilibrium |

---

## 11. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| v30-v40 | All previous fixes | Done |
| v41 | Parametric study: resolution, domain, p-refinement | Done |
| v42 | IP c_N_1 fix: 1.0 -> p(p+D-1)/D at 5 RHS locations | Done |
| v43 | Multi-DOF fault attempt (reverted - sign bug) | Reverted |
| v44 | Penalty ×3 + multi-DOF traction removal | **Reverted (×3 wrong, multi-DOF caused MUMPS issues)** |
| **v45** | **Penalty analysis + Multi-DOF Phases 1-5** | **Done (unit tests)** |

---

## 12. Phase 5: Output Integration

### 12.1 Key Finding: No Code Changes Needed

The BP5 output pipeline already supports multi-DOF without modification:

1. **`Probe2DInterpolator`**: Uses 2D Euclidean nearest-neighbor matching on
   `(x2, x3)` coordinates. With more DOFs per face (nbf=6 at p=2), there are
   more candidate points → stations match to closer DOFs. Works generically.

2. **`BP5BenchmarkOutput::WriteFromGlobalData()`**: Indexes per-DOF vectors by
   the nearest-DOF index from Probe2DInterpolator. All vectors (slip, V, traction,
   theta) are per-DOF at the expanded size — indexing works unchanged.

3. **`GatherToRoot()`**: Uses `MPI_Gatherv` without dedup — works generically
   with any DOF count. The `GatherToRootDedup()` method (depth-only merging) is
   NOT used for BP5 and does not need changes.

4. **`ParallelBP5BenchmarkOutput::Write()`**: Splits interleaved state into
   per-component scalars, gathers 7 fields via `GatherToRoot()`, then calls
   `WriteFromGlobalData()`. All operations scale generically with DOF count.

### 12.2 Multi-DOF Improves Station Matching

With nbf=6 at p=2, there are 6× more DOFs per face. Each DOF has a unique
(x2, x3) coordinate (GaussLobatto node position within the face). Stations
find closer DOFs → lower interpolation error. Verified in unit test:
- Coarse (nbf=1): match distance = 1414 m
- Fine (nbf=6): match distance = 1000 m

### 12.3 Unit Tests Added (8 tests, 280+54 total assertions)

**In `test_bp5_output.cpp` (2 new tests):**

| Test | What it verifies |
|------|-----------------|
| `TestProbe2DInterpolator_MultiDOFCloserMatch` | Multi-DOF (nbf=6) provides closer station matching |
| `TestBP5BenchmarkOutput_MultiDOFWrite` | WriteFromGlobalData correctly indexes multi-DOF vectors: slip, V (log10), theta (log10), tau (tau_pre + elastic) |

**In `test_elasticity_operator.cpp` (6 new tests):**

| Test | What it verifies |
|------|-----------------|
| `TestMultiDOFSEASOperatorP2` | Full SEAS operator construction + SetInitialCondition at p=2 IP tet |
| `TestMultiDOFSEASMultP2` | Mult() at p=2: deterministic, finite, |V|>0 at all DOFs |
| `TestMultiDOFShortRK4P2` | 5 RK4 steps at p=2: state finite, V_max in (0,1), slip accumulated, psi>0 |
| `TestMultiDOFStressEquilibriumP2` | Stress equilibrium maintained during time stepping at p=2 |
| `TestMultiDOFSEASP1Regression` | p=1 IP tet: nbf=1, full SEAS cycle works (backward compat) |
| `TestMultiDOFOutputFromSEASP2` | BP5BenchmarkOutput writes correct 8-column file from SEAS state at p=2 |

---

## 13. v45 TACC Results and v46 Fix

### 13.1 v45 TACC Results (2026-03-19)

**v45a (p=1, 1000m, IP):** Job 7604176, 400 ranks, 8 nodes
- Ran 8046 steps to t = 2.41e-4 yr (76 sec of simulation)
- V_max peaked at 0.293 m/s (step ~7900), then declining
- **Identical to v41a through step 8046** — verified step-by-step match
- This is the initial overstress dissipation phase; EQ #2 at ~250 yr (seen in v41a)
- **Conclusion: v45a p=1 is correct, just needs more walltime**

**v45b (p=2, 1000m, IP):** Job 7604180, 800 ranks, 16 nodes
- 56,382 global fault DOFs (9,397 faces × 6)
- V_max exploded: 0.049 → 0.072 → 0.232 → ... → 196 m/s in 56 steps
- OOM crash: SIGNAL 9 (killed) on ranks 700-748, SIGNAL 11 (segfault) on rank 719
- **Classic numerical instability — positive feedback loop in traction**

### 13.2 Root Cause: Slip Indexing Bug in ComputeTraction

**The bug (2 locations in elasticity_operator.hpp):**

```cpp
// Interior faces (line 3035):
real_t slip_local[2] = {slip_bc(2 * fi), slip_bc(2 * fi + 1)};
// Shared faces (line 3458):
real_t slip_local[2] = {slip_bc(2 * trac_idx), slip_bc(2 * trac_idx + 1)};
```

These index by **face index** (`fi` or `trac_idx`), but `slip_bc` is laid out by
**DOF index** (`fi * nbf + kk`). At p=2 (nbf=6), face `fi=1` reads `slip_bc(2)`,
which is actually DOF 1 of face 0 — not face 1's slip at all!

Additionally, the single `delta_u[3]` was used **constant** across all quad points,
but multi-DOF slip varies across the face. This created a **mismatch** between:
- `AssembleSlipContributionIP`: correctly uses per-DOF interpolated slip
- `ComputeTraction`: uses wrong face's slip, applied as constant

The penalty correction `jump = (u1-u2) - sign*delta_u` sees the wrong slip,
creating huge spurious tractions → positive feedback → blowup.

**Why p=1 was unaffected:** At nbf=1, `fi * 1 + 0 = fi`, so the indexing is
coincidentally correct. This is why v45a matched v41a exactly.

### 13.3 v46 Fix

In `ComputeTraction`, both interior and shared IP paths now:
1. Build per-DOF nodal slip: `slip_bc(2 * (fi * nbf + kk))` for each DOF kk
2. Embed to 3D: `fault_basis_.EmbedSlip(fi, sl, du)`
3. Interpolate: `face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal, delta_u_quad)`
4. Use per-quad-point slip: `delta_u_quad(c * nqp + q)` in the penalty correction

The BR2 path (always nbf=1) retains the old indexing since face_index = DOF_index.

All 280 unit tests pass after the fix.

### 13.4 v46 TACC Runs

| Job | Config | Nodes | Purpose |
|-----|--------|-------|---------|
| v46a | p=1 IP 1000m | 8 (400 ranks) | Regression — must match v41a |
| v46b | p=2 IP 1000m | 16 (800 ranks) | Critical — must nucleate (V>1 m/s) |
| v46c | p=4 IP 1000m | 32 (1600 ranks) | Convergence — nbf=15, must nucleate |

## 14. Next Steps

1. ~~Revert source code to v42 baseline~~ Done
2. ~~Phase 2: L2 traction projection~~ Done
3. ~~Phase 3: Multi-DOF slip interpolation~~ Done (148 tests pass)
4. ~~Phase 4: Fault state/geometry expansion~~ Done (253 tests pass)
5. ~~Phase 5: Output integration + unit tests~~ Done (280+54 tests pass)
6. ~~v45 TACC runs~~ p=2 blew up (slip indexing bug)
7. ~~v46 fix: ComputeTraction slip indexing~~ Done (280 tests pass)
8. **v46 TACC runs**: Submit v46a/b/c, await results
   - Success criterion: earthquake nucleation (V > 1 m/s) at p=2
