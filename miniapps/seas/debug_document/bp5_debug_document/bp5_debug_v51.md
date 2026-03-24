# BP5 Debug v51: Dip-Component Offset — Systematic Cross-Component Coupling Error

**Date**: 2026-03-24
**Status**: u_z = 18.5% of u_x confirmed in elastic solution. Convergence hypothesis DISPROVED: our p=4 on 4000m has 10× worse dip/strike ratio (0.4-8.1%) than Tandem p=4 on the SAME mesh (0.0-0.7%). DG face integrators confirmed identical — the difference must be in volume integrator, boundary face treatment, or quadrature. Investigation ongoing.
**Previous**: v50 (production defaults, first earthquake, ClosedUniform nodes, p=4/p=6 working)
**Branch**: `feature/elasticity`

---

## 1. Problem Statement

The v50 production run (1000m, p=2, 400 ranks) completed the full 1800-year simulation
with 8+ earthquake cycles. While the **strike-slip component** matches Tandem well for the
first 1-2 events, the **dip component** shows systematic errors from t=0:

- Dip slip has a **constant offset** (not accumulating — appears immediately)
- Dip slip rate is **3-5 orders of magnitude too high** (1e-7 vs 1e-11)
- Dip shear stress has a **systematic offset** (varies by station)
- State variable (log10 ψ) is **systematically lower** (~0.5-1 units)
- Strike components start well but **diverge by event 2-3** due to dip contamination

This is the same pattern seen in v46, confirming the issue was never resolved — only
masked by the focus on CFL stability (v49) and node conditioning (v50).

**Both codes use SIPG (identical formulation).** The difference is NOT in the DG method.

---

## 2. Station-by-Station Analysis (Closeup 0-150 years)

### 2.1 Nucleation Station: strk-24dp+10 (x2=-24km, x3=10km)

**Strike (left column):**
- Slip: Good match through first event, slight timing drift by event 2
- V_strike: Matches Tandem envelope, earthquake timing ~5-8s late
- τ_strike: Initial match (21.15 MPa), ~0.5 MPa drift by t=50yr

**Dip (right column):**
- **Slip_dip: CONSTANT -0.015m offset from t=0** — flat line, no accumulation
- **V_dip: ~1e-7 m/s vs Tandem's ~1e-11** — 4 orders of magnitude too high
- **τ_dip: -0.4 MPa offset** (our code more negative)
- **State: ~0.5 units lower in log10(ψ)** — weaker fault

### 2.2 Center Depth: strk+00dp+10 (x2=0km, x3=10km)

**Strike:** Good initial match, τ_strike ~1 MPa low from early interseismic.
**Dip:**
- **Slip_dip: CONSTANT -0.025m offset from t=0**
- **V_dip: ~1e-10 vs ~1e-14** — 4 orders too high
- **τ_dip: -0.15 MPa offset**
- **State: systematically lower**

### 2.3 Near-Nucleation Depth: strk-16dp+10 (x2=-16km, x3=10km)

**Strike:** Good match, slight timing drift.
**Dip:**
- **Slip_dip: CONSTANT -0.015m offset from t=0**
- **V_dip: elevated by ~4 orders**
- **τ_dip: +0.1 MPa offset** (opposite sign from center!)
- **State: lower**

### 2.4 Far Depth: strk+16dp+10 (x2=+16km, x3=10km)

**Strike:** Good match initially, τ_strike ~1 MPa low.
**Dip:**
- **Slip_dip: CONSTANT -0.046m offset from t=0** (larger than center)
- **V_dip: elevated**
- **τ_dip: +0.1 MPa offset**

### 2.5 Center Surface: strk+00dp+00 (x2=0km, x3=0km)

**Strike:** Reasonable match (surface VS zone, slow dynamics).
**Dip:**
- **Slip_dip: CONSTANT -0.05m offset from t=0**
- **V_dip: ~1e-10 vs ~1e-13** — 3 orders too high
- **τ_dip: -0.6 MPa offset**

### 2.6 Far Surface: strk+36dp+00 (x2=+36km, x3=0km)

**Strike:** Good match.
**Dip:**
- **Slip_dip: -0.12m offset, GROWING over time** (not constant like depth stations)
- **V_dip: ~1e-10** — elevated
- **τ_dip: -1.2 MPa offset** (largest offset of any station)

### 2.7 Far Surface Opposite: strk-36dp+00 (x2=-36km, x3=0km)

**Strike:** Reasonable match.
**Dip:**
- **Slip_dip: +0.2m offset — OPPOSITE SIGN from all other stations!**
- **V_dip: ~1e-3** — orders of magnitude too high
- **τ_dip: +1.0 MPa offset** (also opposite sign)
- This station has the WORST dip error

---

## 3. Spatial Pattern of the Dip Offset

| Station | Dip slip offset | Dip τ offset | Pattern |
|---------|----------------|-------------|---------|
| strk-36dp+00 | **+0.20 m** | **+1.0 MPa** | Positive (anomalous) |
| strk-24dp+10 | -0.015 m | -0.4 MPa | Negative, small |
| strk-16dp+10 | -0.015 m | +0.1 MPa | Mixed sign |
| strk+00dp+10 | -0.025 m | -0.15 MPa | Negative, moderate |
| strk+00dp+00 | -0.05 m | -0.6 MPa | Negative, larger at surface |
| strk+16dp+10 | -0.046 m | +0.1 MPa | Mixed sign |
| strk+36dp+00 | -0.12 m | -1.2 MPa | Negative, largest |

Key observations:
1. **Offset INCREASES with distance from nucleation zone** — small at x2=-24, large at x2=±36
2. **Offset is LARGER at surface (x3=0) than at depth (x3=10)**
3. **strk-36dp+00 has OPPOSITE sign** — suggesting an asymmetry in the fault geometry or BC coupling
4. **The offset is IMMEDIATE (t=0)** — not accumulated over time

---

## 4. Cascade of Errors

The dip offset creates a cascade that contaminates all components:

```
Spurious dip displacement (from elastic solve at t=0)
    ↓
Non-zero dip traction (from ComputeTraction)
    ↓
Elevated dip V (friction law balances η·V = τ_dip)
    ↓
Lower state variable ψ (aging law: dψ/dt depends on |V|, not just V_strike)
    ↓
Weaker fault (lower ψ → lower friction strength)
    ↓
Strike τ deficit (stress redistributes across weaker fault)
    ↓
Earlier/different earthquake timing (events 2+ diverge)
```

---

## 5. Root Cause Analysis

### 5.1 What We Know

1. **Both codes use SIPG** — the DG formulation is identical
2. **Initial conditions match to 6 digits** (tau, V, psi at t=0)
3. **The dip offset is IMMEDIATE** — present from the first time step
4. **The offset is SPATIALLY varying** — not a uniform bias
5. **The offset depends on position relative to fault boundaries** — largest at edges
6. **This was also seen in v46** — not a new issue

### 5.2 Hypothesis: Fault Basis / Projection Error

The most likely root cause is in the **FaultBasis** class, which handles:
- `EmbedSlip(face_idx, slip_local[2], du[3])`: converts (slip_dip, slip_strike) → (Δu_x, Δu_y, Δu_z)
- `ProjectTraction(face_idx, T_global[3], tau_local[2])`: converts (T_x, T_y, T_z) → (τ_dip, τ_strike)

If the rotation between local (dip, strike) and global (x, y, z) frames has even a small
error, it would:
- Mix strike slip into the dip component (creating spurious dip displacement)
- Mix global traction components into the dip traction
- The error would be **spatially varying** because the fault geometry varies

Tandem constructs its tangent/normal basis differently. Any difference in:
- Normal vector direction or sign
- Tangent vector definition (which direction is "dip" vs "strike")
- Orthogonalization method
...would produce exactly this pattern.

### 5.3 Hypothesis: Far-Field BC Coupling into Dip

The spatial pattern (larger offset at edges, opposite sign at strk-36) suggests
the **far-field Dirichlet BCs** may be coupling into the dip component through the
3D elasticity tensor. If the prescribed far-field displacement has a non-zero dip
component (even a tiny one from numerical imprecision), the penalty enforcement
at far-field boundary faces would create a dip traction that propagates to the fault.

### 5.4 Hypothesis: Initial Elastic Solve Residual

At t=0, slip=0 everywhere, so the elastic solve should give u=0 and traction=0.
But if there's any numerical residual in the first solve (from the DG formulation,
solver tolerance, or BC enforcement), it could produce a non-zero dip displacement
at fault faces. The penalty term in ComputeTraction would then amplify this residual
into a spurious dip traction.

---

## 6. Investigation Plan

### Priority 1: Verify t=0 traction

Dump the per-DOF traction at the first ODE evaluation (slip=0). The dip traction
should be EXACTLY zero. If it's not, the error originates in the elastic solve or
traction computation.

### Priority 2: Compare FaultBasis with Tandem

Line-by-line comparison of:
- Normal vector construction
- Tangent1 (dip) and tangent2 (strike) definition
- EmbedSlip rotation matrix
- ProjectTraction rotation matrix

Files to compare:
- Our: `miniapps/seas/fault/fault_basis.hpp`
- Tandem: `app/form/FacetFunctionalInterior.h` or similar

### Priority 3: Compare far-field BC implementation

Check whether our Dirichlet BC for far-field faces:
- Prescribes only strike displacement (correct)
- Or prescribes all 3 components (potentially introducing dip coupling)

Compare with Tandem's `boundary_linear` implementation.

### Priority 4: Dump displacement at t=0

After the first elastic solve (slip=0), dump the displacement field at fault faces.
Check if there's a non-zero dip component. If so, trace it back to the BC or DG assembly.

---

## 7. FaultBasis Comparison: IDENTICAL — NOT the Root Cause

### 7.1 Line-by-Line Comparison

Comprehensive comparison of fault tangent/normal basis construction between our code
(`fault_basis.hpp`) and Tandem (`Curvilinear.cpp:facetBasis()`):

**Both codes compute for BP5 (ref_normal=(0,-1,0), up=(0,0,1)):**

| Vector | Formula | Tandem | SEAS-MFEM | Match? |
|--------|---------|--------|-----------|--------|
| strike | normalize(up × n) | (1, 0, 0) | (1, 0, 0) | ✓ Identical |
| dip | strike × n | (0, 0, -1) | (0, 0, -1) | ✓ Identical |
| normal | ref_normal | (0, -1, 0) | (0, -1, 0) | ✓ Identical |

**Component ordering:**

| Index | Tandem | SEAS-MFEM | Match? |
|-------|--------|-----------|--------|
| 0 | dip | tangent1 = dip | ✓ |
| 1 | strike | tangent2 = strike | ✓ |

**EmbedSlip (both codes):**
```
delta_u[d] = slip_local[0] * dip[d] + slip_local[1] * strike[d]
           = slip_dip * (0,0,-1) + slip_strike * (1,0,0)
           = (slip_strike, 0, -slip_dip)
```

**ProjectTraction (both codes):**
```
tau_local[0] = T · dip    = T · (0,0,-1) = -T_z     (dip traction)
tau_local[1] = T · strike = T · (1,0,0)  = T_x      (strike traction)
```

**Rotation matrix R (columns = dip, strike):**
```
R = [ 0  1 ]
    [ 0  0 ]
    [-1  0 ]
```
R^T × R = I (orthonormal). EmbedSlip uses R, ProjectTraction uses R^T. **No transposition error.**

### 7.2 Only Structural Difference: Sign-Flip Post-Processing

Tandem's `AdapterBase::prepare()` has an additional step: when the raw mesh normal
disagrees with ref_normal (i.e., `sign_flipped = true`), Tandem flips the entire
fault_basis_q (all 3 columns negated). This is self-canceling: the negation in
ProjectTraction cancels with the negation in EmbedSlip, so the friction law sees
the same physics regardless of face orientation.

SEAS-MFEM handles this differently — it computes tangent vectors from the already-
corrected normal, so tangent vectors are always consistent. On a planar fault (BP5),
both approaches give identical results.

### 7.3 Conclusion

**The FaultBasis is NOT the source of the dip offset.** The rotation matrices, component
ordering, and sign conventions are identical between the two codes. Cross-component
coupling does not originate from the fault basis projection.

---

## 8. History of Dip-Related Findings in Previous Debug Documents

| Version | What was found | Conclusion at the time | Still valid? |
|---------|---------------|----------------------|-------------|
| **v30** | Dip negation bug: `dip = -(s×n)` | **FIXED** to `dip = s×n` | ✓ Fix confirmed correct |
| **v31** | Interior Dirichlet sign bug → spurious τ_dip | **FIXED** | ✓ Fix confirmed correct |
| **v34** | τ_dip grows linearly, antisymmetric at edges | "Discretization artifact at p=1" | **✗ Still present at p=2** |
| **v39/v40** | IP traction sign convention with flipped normal | **FIXED** | ✓ Fix confirmed correct |
| **v45** | Penalty ×3 from reference element conventions | Analyzed, applied in v47 | ✓ Separate issue |
| **v46** | Multi-DOF slip indexing bug at p≥2 | **FIXED** | ✓ Fix confirmed correct |
| **v48** | Shared face normal flips between ranks | "By design, NOT a bug" | ✓ Confirmed correct |
| **v49** | CalcOrtho=(0,+1,0) vs basis.normal=(0,-1,0) | "Correct by design" | ✓ Confirmed correct |

**Critical re-evaluation of v34:** The v34 conclusion that τ_dip growth is a "discretization
artifact at p=1, resolved at p=6" was never re-tested at p=2. Our v50 production run IS
at p=2, and the dip offset is still present with the same spatial pattern (antisymmetric
at edges, growing with distance from fault center). The v34 explanation is **insufficient**
— the dip offset is NOT purely a p-refinement issue.

---

## 9. Far-Field BC Comparison: IDENTICAL — NOT the Root Cause

### 9.1 Prescribed Displacement

Both codes prescribe identical far-field displacement:

| Face | Our code | Tandem | Match? |
|------|----------|--------|--------|
| Far-field (exterior) | `u_D = (sgn(Y)*Vp*t/2, 0, 0)` | `return Vh, 0, 0` with `Vh = sgn(y)*Vp*t/2` | ✓ |
| Y=0 non-fault (interior) | `u_D_jump = (Vp*t, 0, 0)` | Same jump via sign flip | ✓ |
| **Dip component** | **Explicitly zero** | **Explicitly zero** | ✓ |

**The prescribed BC has NO dip (z) component in either code.**

### 9.2 Boundary Surfaces

Both codes apply Dirichlet to the same surfaces:
- x = ±Lx faces → Dirichlet
- y = ±Ly faces → Dirichlet
- Non-fault Y=0 interior faces → Interior Dirichlet (jump-based)
- z = 0 (top) → Natural (free surface)
- z = -Lz (bottom) → Natural

Tandem Gmsh: Physical Surface(1)={top,bottom}→Natural, Surface(5)={far-field}→Dirichlet
Our code: BCMode::FarField matches correctly for Tandem's Gmsh mesh.

### 9.3 DG Penalty Structure

Both codes enforce all 3 components simultaneously:
- Tandem: `sigma_hat = C:grad(u) + penalty * (u - u_D) * n_unit`
- Ours: `elvec(idx) += wq_penalty * u_D[i] * shape(k)` for all `i` in dim

Same penalty formula (after v47 ×3 correction).

### 9.4 `boundary_linear` Flag

Tandem's `boundary_linear = true` (bp5.toml) is a mode guard for discrete Green's function
optimization. It does NOT change the DG assembly. The standard QD mode evaluates the full
boundary function at each solve.

### 9.5 Conclusion: Far-field BC is NOT the source of dip offset.

---

## 10. Revised Understanding: Dip Slip Accumulates During Earthquakes

### 10.1 Re-interpretation of the "Constant Offset"

The closeup plots (0-150 years) show dip slip that APPEARS constant. But the time
resolution is years — the dip slip actually accumulates during each **coseismic event**
(~30 seconds) and then stays constant during the interseismic period (~100 years).

The "constant" offset at the nucleation station (-0.015m) is actually -0.015m of dip
slip accumulated during the FIRST earthquake.

### 10.2 Mechanism: Vector Friction Law Couples Dip and Strike

In BP5's vector rate-state friction, the slip velocity direction follows the traction:
```
V_vec = V_scalar * (τ_vec / |τ_vec|)
```
So `V_dip / V_strike = τ_dip / τ_strike`.

If our code computes a spurious τ_dip = 0.15 MPa while τ_strike = 20 MPa:
```
V_dip / V_strike = 0.15 / 20 = 0.0075
```
During a 30-second earthquake with V_strike ≈ 0.5 m/s:
```
V_dip ≈ 0.004 m/s
dip_slip ≈ 0.004 × 30 ≈ 0.12 m
```
This matches the observed offset at far-field stations!

At the nucleation station (smaller τ_dip offset): dip_slip ≈ 0.015m — also matches.

### 10.3 The Root Question

**Where does the ~0.15 MPa spurious τ_dip come from?**

For a pure strike-slip loading on a planar fault with normal n=(0,-1,0):
- Strike loading creates σ_xy, which produces T_x (strike traction) ✓
- Strike loading does NOT create σ_zy, so T_z (dip traction) should be ZERO
- Poisson coupling creates σ_zz from ε_xx, but σ_zz does NOT produce traction on a y-normal face

Yet our ComputeTraction reports non-zero T_z. This can only come from:
1. **Mesh geometry**: tetrahedral elements are never perfectly aligned — numerical integration
   produces small cross-component stress gradients
2. **DG formulation**: the penalty/consistency/symmetry terms in the traction computation
   amplify these small cross-component errors

Tandem uses the SAME DG formulation on the SAME mesh. If mesh geometry were the sole cause,
Tandem would have the same τ_dip. **The difference must be in HOW the traction is computed.**

---

## 11. Remaining Hypotheses (FaultBasis + BC Both Ruled Out)

### 11.1 Hypothesis A: Traction Post-Processing Cross-Component Error (MOST LIKELY)

Both codes solve the same K*u = f. The displacement u is identical (same formulation,
same mesh, same solver). But the TRACTION is computed differently:

**Tandem**: Computes traction WITHIN the DG operator application. The consistency term
`{σ·n}` and penalty correction are evaluated together as part of the same kernel. The
cross-component coupling in the elasticity tensor is handled consistently.

**Our code**: Post-processes traction SEPARATELY in ComputeTraction. The stress average
`{σ·n}` is computed from the displacement gradient, and the penalty correction is subtracted.
These are separate evaluations that may handle the λ (cross-coupling) term differently.

For a tet mesh with non-aligned faces, the stress gradient ∇u has all 9 components
non-zero. The traction `{σ·n}` = `{(λ tr(ε)I + 2με)·n}` includes:
```
T_z = λ*(ε_xx + ε_yy + ε_zz) * n_z + 2μ * (ε_zx*n_x + ε_zy*n_y + ε_zz*n_z)
```
On a y-normal face: n = (0, -1, 0), so:
```
T_z = -2μ * ε_zy = -μ * (∂u_z/∂y + ∂u_y/∂z)
```
The λ term DROPS OUT for dip traction on a y-normal face. The dip traction depends
ONLY on the off-diagonal strain ε_zy.

For a perfect strike-slip solution, ∂u_z/∂y = ∂u_y/∂z = 0, so T_z = 0. But numerical
errors in ∂u_z/∂y (from DG element-wise polynomial approximation of the displacement
field) produce non-zero T_z.

**Critical question**: Does Tandem's operator-internal traction computation cancel this
numerical error more effectively than our separate post-processing?

### 11.2 Hypothesis B: Initial Condition — tau_pre Dip Component

In BP5, the pre-stress τ_pre is computed from the steady-state friction law:
```
τ_pre = σ_n * [f0 + (a-b) * ln(V_init/V0)]
```
This is a SCALAR. The direction is V_init / |V_init|.

V_init = (V_zero, V_init_strike) where V_zero ≈ 0.

If V_zero is EXACTLY 0: τ_pre_dip = 0 (correct)
If V_zero is tiny but non-zero: τ_pre_dip = τ_pre * V_zero / |V_init| (negligible)

**Need to verify**: What is V_zero set to in our code? If it's 0, this is not the cause.

---

## 12. Updated Investigation Plan

### Priority 1: Dump per-component traction at first coseismic step

At the first ODE evaluation where slip is non-zero (first RK stage of step 1),
dump for each fault face DOF:
- `T_global = (T_x, T_y, T_z)` from ComputeTraction
- `T_stress = {σ·n}` component (stress average only)
- `T_penalty = penalty * correction` component
- `tau_dip = T · (0,0,-1)` and `tau_strike = T · (1,0,0)`

If T_z is non-zero: is it from T_stress or T_penalty?
If T_stress: the displacement gradient has spurious ∂u_z/∂y — DG discretization error
If T_penalty: the penalty correction has cross-component leakage

### Priority 2: Compare with Tandem's traction at same time

Run Tandem on the same mesh with the same parameters, dump traction components at the
same time step. If Tandem's T_z = 0 and ours ≠ 0, the difference is in the traction
computation method. If both have non-zero T_z, the mesh geometry is the cause.

### Priority 3: Check V_zero value

Verify that the dip component of V_init is exactly 0.0 in our initialization.

---

## 13. v51 Diagnostic Results: Per-Component Traction Dump

### 13.1 Test Configuration

Serial (1 rank), 1000m mesh, p=2, IP method. The `--diag-dip-traction` flag dumps
per-component traction (global xyz + local dip/strike) at the first non-zero-slip
ComputeTraction call. 102,520 fault face DOFs reported.

### 13.2 SMOKING GUN: 21% Cross-Component Contamination in {σ·n}

**Dip/Strike ratio in T_stress = {σ·n} across 84,322 fault faces:**

| Statistic | |Tz/Tx| ratio |
|-----------|-------------|
| Median | **20.6%** |
| Mean | **21.3%** |
| Max | 35.8% |
| Min | 0.0% |

For every 1 Pa of strike traction from {σ·n}, there is **0.21 Pa of spurious dip
traction**. This is a 21% cross-component contamination from the DG stress average
on tetrahedral elements.

### 13.3 Penalty Correction is NOT the Source

| Comparison | Median ratio |
|-----------|-------------|
| |Tc_x / Ts_x| (penalty / stress) | 1.35% |

The penalty correction is **negligible** at most faces. The dip contamination comes
almost entirely from **{σ·n}** (the stress average), NOT from the penalty term.

Exception: nucleation boundary face (fi=16, x=-30.7km) where penalty correction
is 270% of stress — catastrophic cancellation at that specific face.

### 13.4 Sample Data (Faces with Highest Dip/Strike Ratio)

```
fi    x(km)    z(km)  Tstress_x(Pa)  Tstress_z(Pa)  |Tz/Tx|
 2    -14.5    -22.9       -16.11          5.78      35.8%
14    -14.5    -21.2       -23.95          8.35      34.9%
 3    -15.5    -24.6       -12.12          4.17      34.4%
 6    -34.5    -24.6       -10.63         -3.65      34.3%
 4    -15.0     -6.8      -190.50        -28.26      14.8%
 0    -12.5    -14.0       -64.42          4.95       7.7%
```

The highest ratios (30-36%) are at faces near the VW zone boundary at depth 20-27km.
The absolute magnitudes are small (Pa range at the first step), but the RATIO is large.

### 13.5 Why This Causes Dip Slip During Earthquakes

At the first time step, the absolute dip traction (~5 Pa) is negligible vs tau_pre
(~20 MPa). But during earthquakes:

1. tau_strike drops from 20 → 8 MPa (stress drop = 12 MPa)
2. The elastic traction T_elastic grows with slip magnitude
3. The 21% contamination ratio means T_dip ≈ 0.21 × T_elastic_strike
4. Even if the ratio drops to ~1% during coseismic (better cancellation at large slip):
   ```
   tau_dip_contam ≈ 0.01 × 12 MPa = 0.12 MPa
   V_dip / V_strike = 0.12 / 8 = 1.5%
   dip_slip_per_event ≈ 0.015 × 0.5 m/s × 30s ≈ 0.23 m
   ```
5. This matches the observed 0.05-0.12 m per-event dip offset

### 13.6 V_zero Check

`V_zero = 1e-20` in bp5_params.hpp. Produces `τ_pre_dip = 2e-4 Pa`, 750,000× too small
to explain the offset. **NOT the cause.**

---

## 14. Tandem Comparison: Same Formula, Same Contamination?

### 14.1 Tandem's Traction Formula

**File:** `tandem/app/kernels/elasticity.py`, line 242-244

```python
traction_q = 0.5 * (sigma_0.n + sigma_1.n) + c0 * (u_0 - u_1 - f_q)
```

where `c0 = -penalty`. This is **IDENTICAL** to our formula:
`T = {σ·n} - penalty*(jump - slip)`.

**Tandem uses the same DG traction formula as our code.** The {σ·n} stress average
on tet meshes produces the same cross-component contamination in both codes.

### 14.2 Tandem's Traction Projection

**File:** `tandem/app/kernels/elasticity_adapter.py`, lines 24-27

```python
traction[k,p] = M_inv[k,l] * Σ_q w[q] * |n|[q] * e[l,q] * Σ_o traction_q[o,q] * fault_basis[o,p,q]
```

This is an **L2 projection combined with local-frame projection** in a single kernel.
The contraction `traction_q[o,q] * fault_basis[o,p,q]` projects global → local at
each quadrature point, then L2-averages.

**Our code does the same** but in two separate steps:
1. L2-project global T_q → per-DOF T_global (via `GalerkinProject`)
2. Project T_global → local (τ_dip, τ_strike) per DOF (via `ProjectTraction`)

For a **flat fault** (BP5: planar y=0), the fault basis is CONSTANT over the face.
The ordering (project-then-average vs average-then-project) is mathematically
equivalent. This difference does NOT explain the dip offset.

### 14.3 Critical Question: Does Tandem Also Have 21% Contamination?

**Yes, Tandem MUST have the same {σ·n} contamination** because:
- Same DG formulation (SIPG)
- Same mesh (tetrahedral, non-aligned)
- Same stress average formula
- Same penalty formula

Yet Tandem's dip slip is essentially zero. **Why?**

### 14.4 Three Possible Explanations

**A. Tandem's contamination magnitude differs at coseismic scale.**

During the earthquake, the contamination ratio may differ between the two codes
due to subtle differences in penalty magnitude, solver accuracy, or element-level
displacement field. Even a factor-of-2 difference in coseismic contamination ratio
(e.g., 0.5% vs 1%) would halve the dip slip accumulation.

**B. Tandem's friction law handles small dip traction differently.**

The vector RSF law computes V_vec = V × τ_vec / |τ_vec|. If |τ_dip| << |τ_strike|,
then V_dip = V × τ_dip / |τ|. The RELATIVE magnitude matters. Tandem may:
- Normalize differently (using |τ| = sqrt(τ_d² + τ_s²) vs τ_s alone)
- Have a threshold for small-component velocity

**C. The contamination accumulates differently due to state variable coupling.**

The state variable ψ evolves as dψ/dt = f(V, ψ, L), where V = |V_vec| =
sqrt(V_dip² + V_strike²). A non-zero V_dip increases |V|, which changes ψ evolution,
which changes τ_friction, which changes V_strike. This coupling could make the
dip contamination self-correcting in Tandem's implementation but self-reinforcing
in ours, depending on subtle numerical differences.

---

## 15. Friction Law Comparison: Three Differences Found

### 15.1 Velocity Direction Computation: IDENTICAL

Both codes compute:
```
V_dip = V_scalar × τ_dip / |τ|
V_strike = V_scalar × τ_strike / |τ|
```
where `|τ| = sqrt(τ_dip² + τ_strike²)`.

No thresholds, no clipping, no regularization for small dip/strike ratios.
The Brent solver for V_scalar uses the same residual equation and bracketing.

**Hypothesis B (friction law) is RULED OUT** — the velocity direction logic is identical.

### 15.2 State Variable Evolution: Nearly Identical

Both codes: `dψ/dt = b·V₀/L × (exp((f₀-ψ)/b) - V/V₀)` where `V = |V_vec| = sqrt(V_dip² + V_strike²)`.

**One minor difference**: SEAS-MFEM caps `exp((f₀-ψ)/b)` at `exp(20)` to prevent overflow.
Tandem has no cap. This affects post-earthquake healing rate when ψ << f₀, but is
unlikely to cause the dip offset during interseismic or coseismic.

### 15.3 Normal Stress: CRITICAL DIFFERENCE

| | Our Code | Tandem |
|--|----------|--------|
| **σ_n** | **Constant 25 MPa** | **-sn_elastic + sn_pre** (varies) |
| Source | `params.sigma_n` (line 130) | `traction(node, 0)` from DG solver |

**Tandem uses the normal traction perturbation from the elasticity solver.**
Our code discards it — we compute the normal traction in `ProjectTraction` but only
pass the tangential (dip, strike) components to the friction law.

Impact: The DG stress average has cross-component contamination in the NORMAL direction
too (σ_yy perturbation from pure x-loading, via the λ·tr(ε)·I term). In Tandem, this
perturbation modifies σ_n, which changes friction strength:
- Positive σ_n perturbation → STRONGER fault → LESS slip → dampens dip motion
- Negative σ_n perturbation → WEAKER fault → MORE slip

This normal stress feedback provides a natural damping mechanism that our code lacks.
However, the perturbation magnitude (~10-100 Pa vs 25 MPa baseline = 0.0004%) is very
small, so this alone likely does NOT explain the full dip offset.

### 15.4 Radiation Damping: IDENTICAL

Both codes: `η·V` where V = V_scalar (from Brent solver). Same formula, same convention.

### 15.5 Summary Table

| Component | Compared? | Result |
|-----------|----------|--------|
| Velocity direction (V_dip/V_strike) | ✓ Line-by-line | **IDENTICAL** |
| Scalar RSF equation | ✓ | **IDENTICAL** (same Brent solver) |
| |τ| normalization | ✓ | **IDENTICAL** (L2 norm) |
| State evolution dψ/dt | ✓ | Nearly identical (exp cap difference) |
| **Normal stress σ_n** | ✓ | **DIFFERENT** — ours constant, Tandem varies |
| Radiation damping η·V | ✓ | **IDENTICAL** |
| Thresholds/clipping | ✓ | Both have none for dip motion |

### 15.6 Conclusion

The friction law is NOT the primary cause of the dip offset. Both codes produce
`V_dip ∝ τ_dip / |τ|` identically. The root cause remains the 21% cross-component
contamination in {σ·n} from the DG stress average on tet meshes.

The normal stress difference (Section 15.3) is a secondary finding that should be
addressed separately — it could improve overall benchmark agreement.

---

## 16. What Has Been Ruled Out

| Component | Compared? | Result |
|-----------|----------|--------|
| FaultBasis (normal, dip, strike vectors) | ✓ Line-by-line | **IDENTICAL** |
| EmbedSlip rotation matrix | ✓ | **IDENTICAL** |
| ProjectTraction rotation matrix | ✓ | **IDENTICAL** |
| Component ordering (index 0=dip, 1=strike) | ✓ | **IDENTICAL** |
| Far-field BC prescribed displacement | ✓ | **IDENTICAL** (zero dip) |
| Boundary surface assignment | ✓ | **IDENTICAL** |
| DG penalty structure | ✓ | **IDENTICAL** |
| DG traction formula | ✓ | **IDENTICAL** ({σ·n} - penalty*(jump-slip)) |
| Traction projection ordering | ✓ | Equivalent for flat faces |
| `boundary_linear` flag effect | ✓ | No effect on assembly |
| Sign fixes (v30, v31, v39, v40) | ✓ | All still correct |
| Multi-DOF indexing (v46) | ✓ | Fixed, still correct |
| V_zero initialization | ✓ | 1e-20, negligible |
| Velocity direction V_dip/V_strike | ✓ | **IDENTICAL** |
| Scalar RSF solver (Brent) | ✓ | **IDENTICAL** |
| State evolution dψ/dt | ✓ | Nearly identical (exp cap) |
| Radiation damping η·V | ✓ | **IDENTICAL** |

---

## 17. Tests in Progress

### v51a: Zero dip traction (production, 48hr, 400 ranks)

`--zero-dip-traction --diag-coseismic-dip`

Sets τ_dip = 0 after ComputeTraction, before friction law. If this eliminates
the dip offset AND improves strike match → confirms the mechanism and provides
a workaround for production runs.

### v51b: Coseismic dip diagnostic (4hr, 400 ranks)

`--diag-coseismic-dip`

Baseline run that dumps the mean and max |τ_dip/τ_strike| ratio when V_max > 0.1.
Quantifies the actual contamination during the earthquake that drives dip slip.

---

## 18. Updated Investigation Plan

### Completed
- ✅ FaultBasis comparison (IDENTICAL)
- ✅ Far-field BC comparison (IDENTICAL)
- ✅ Tandem traction formula comparison (IDENTICAL)
- ✅ Per-component traction diagnostic (21% contamination found)
- ✅ Friction law comparison (velocity direction IDENTICAL, normal stress DIFFERENT)
- ✅ V_zero check (negligible)

### Waiting for results
- ⏳ v51a: zero-dip-traction production run
- ⏳ v51b: coseismic dip/strike ratio diagnostic
- ⏳ v51d: p=1 1000m BLR 1e-12 (missed test from v47)

### Remaining
1. **If v51a confirms**: implement `--zero-dip-traction` as default for BP5 production
2. **Normal stress coupling**: add elastic normal stress perturbation to friction law (match Tandem)
3. **Root cause of 10× worse dip ratio vs Tandem** (Section 19.5): compare volume integrator,
   boundary face treatment, quadrature orders
4. **If v51d nucleates**: p=1 with ×3 penalty is viable — revisit dip contamination at p=1

---

## 17. Revision History

| Version | Change | Status |
|---------|--------|--------|
| v46 | First observation of dip offset and timing drift | Observed |
| v47 | Penalty ×3 correction — improved strike, did NOT fix dip | Partial fix |
| v49 | CFL + V-guard — stabilized parallel runs | Done |
| v50 | ClosedUniform nodes — fixed p=4 conditioning | Done |
| v50 prod | Full 1800yr run: 8 earthquakes, good strike match, **dip offset persists** | Confirmed |
| **v51** | **Dip offset analysis: systematic, spatially varying, accumulates during EQs** | **INVESTIGATING** |
| v51 | FaultBasis comparison: IDENTICAL to Tandem — NOT the root cause | **RULED OUT** |
| v51 | Far-field BC comparison: IDENTICAL to Tandem — NOT the root cause | **RULED OUT** |
| v51 | Tandem traction formula comparison: IDENTICAL (same {σ·n} - penalty*(jump-slip)) | **CONFIRMED** |
| v51 | Previous debug doc review: v30-v49 dip fixes all confirmed correct | **REVIEWED** |
| v51 | v34 "p=1 artifact" conclusion invalidated — offset persists at p=2 | **RE-EVALUATED** |
| v51 | Revised understanding: dip slip accumulates during EQs via vector friction | **KEY INSIGHT** |
| **v51 diag** | **SMOKING GUN: {σ·n} has 21% dip/strike cross-component contamination** | **FOUND** |
| v51 diag | Penalty correction is negligible (1.35% of stress) — NOT the source | **CONFIRMED** |
| v51 diag | V_zero = 1e-20 — negligible, NOT the cause | **CONFIRMED** |
| v51 | Tandem has SAME contamination from same DG formula on same mesh | **KEY QUESTION** |
| v51 | Three hypotheses for why Tandem doesn't accumulate: (A) magnitude, (B) friction, (C) state coupling | HYPOTHESIZED |
| **v51** | **Friction law comparison: velocity direction IDENTICAL, normal stress DIFFERENT** | **COMPLETED** |
| v51 | Hypothesis B (friction law) RULED OUT — velocity direction logic identical | **RULED OUT** |
| v51 | **Found: Tandem uses elastic σ_n, we use constant 25 MPa** | **SECONDARY FINDING** |
| v51a | Zero-dip-traction production test (48hr, 400 ranks) | **SUBMITTED** |
| v51b | Coseismic dip/strike ratio diagnostic (4hr, 400 ranks) | **SUBMITTED** |
| **v51d** | **p=1 1000m BLR 1e-12: disambiguate v47a V-decay (BLR artifact vs penalty stiffness)** | **PLANNED** |
| v51f/f2 | Strategy 1 (stress-only traction): BLOWUP at p=2 and p=4 | **DISPROVED** |
| **v51c** | **DEFINITIVE: u_z = 18.5% of u_x in the elastic solution itself** | **ROOT CAUSE** |
| v51c | Source: K-matrix cross-coupling from non-fault interior face DG terms | **ANALYZED** |
| v51c | f_z = 0 confirmed on fault faces → u_z comes from K_zx coupling | **CONFIRMED** |
| v51c | Initial conclusion: DG discretization error scaling with p | HYPOTHESIZED |
| **v51c+** | **DISPROVED: p=4 on 4000m has 10× worse ratio than Tandem p=4 on same mesh** | **CODE DIFF EXISTS** |
| v51c+ | Tandem ALSO has non-zero dip slip — it's a real 3D effect, but our ratio is 10× too high | **KEY FINDING** |
| v51c+ | DG face integrators confirmed identical → difference must be in volume/BC/stabilization | **NARROWED** |
| v51c+ | Next: compare volume integrator, boundary face treatment, quadrature orders | **PLANNED** |

---

## 19. DEFINITIVE FINDING: u_z Contamination in Elastic Solution (v51c)

### 19.1 Setup

Added `--diag-uz-fault` flag: dumps u_x, u_y, u_z at all fault face centroids after
the first elastic solve with non-zero slip. For pure strike-slip, u_z should be exactly 0.

### 19.2 Results: u_z = 18.5% of u_x

| Component | Max value | Role | Expected |
|-----------|----------|------|----------|
| u_y (normal) | 9.24e-5 m | Fault-normal (Poisson) | Non-zero ✓ |
| u_x (strike) | 2.92e-5 m | Strike displacement | Non-zero ✓ |
| **u_z (dip)** | **5.39e-6 m** | **Dip displacement** | **Should be 0** ✗ |

**max |u_z / u_x| = 18.5%** — the displacement solution ITSELF has massive dip contamination.

### 19.3 Verified: f_z = 0 on Fault Faces

The RHS has zero z-component on fault faces:
- **Penalty term**: `penalty * sign * delta_u_q[i=z] * shape(k)` — delta_u_q[z] = 0 for pure x-slip ✓
- **Symmetry term**: `trac_20 = λ * dshape(k,2) * nor(0) + μ * dshape(k,0) * nor(2)` — both nor(0) = 0 and nor(2) = 0 for fault normal (0, ±n_y, 0) ✓

### 19.4 Source: K-Matrix Cross-Coupling from Interior Faces

The DG bilinear form K includes face integrals on ALL interior faces (not just fault).
On non-fault interior faces, the face normal n = (n_x, n_y, n_z) has non-zero x,z components.
The SIPG consistency/symmetry terms:

```
-∫_F {σ(u)·n} · [v] dS        (consistency)
-ε ∫_F {σ(v)·n} · [u] dS      (symmetry)
```

create cross-component coupling K_xz ≠ 0 through the elasticity tensor (λ term couples
div(u) = ∂u_x/∂x + ∂u_y/∂y + ∂u_z/∂z to all components). With f_z = 0 but K_xz ≠ 0:

```
K_zz * u_z + K_zx * u_x = 0
→ u_z = -K_zz⁻¹ * K_zx * u_x ≠ 0
```

### 19.5 Convergence Hypothesis — DISPROVED by p=4 4000m Results

Initial hypothesis: the cross-coupling is a discretization error scaling with p, and higher p
would reduce the contamination. This was **DISPROVED** by the p=4 4000m results:

**Same mesh (4000m), same order (p=4) — 10× worse than Tandem:**

| Run | Mesh | Order | Dip/Strike ratio (accumulated) |
|-----|------|-------|-------------------------------|
| **Our p=4** | **4000m** | **p=4** | **0.4-8.1%** |
| **Tandem p=4** | **4000m** | **p=4** | **0.0-0.7%** |
| **Tandem p=6** | **4000m** | **p=6** | **0.0-0.6%** |
| Our p=2 | 1000m | p=2 | 0.3-5.7% |

**Station-by-station dip/strike ratio comparison:**

| Station | Our p=2 | Our p=4 | Tandem p=4 | Tandem p=6 |
|---------|---------|---------|------------|------------|
| Nucleation (x2=-24, x3=10) | 0.5% | 0.4% | 0.2% | 0.2% |
| Near-nuc depth (x2=-16, x3=10) | 0.3% | 0.5% | 0.0% | 0.0% |
| Center depth (x2=0, x3=10) | 0.4% | 2.5% | 0.1% | 0.1% |
| Center surface (x2=0, x3=0) | 0.8% | **6.8%** | 0.2% | 0.2% |
| Far depth (x2=16, x3=10) | 0.6% | **2.4%** | 0.2% | 0.3% |
| Far surface (x2=-36, x3=0) | 5.7% | **8.1%** | 0.7% | 0.6% |

**Critical finding**: Our p=4 on 4000m (0.4-8.1%) is **WORSE** than our p=2 on 1000m
(0.3-5.7%), and both are much worse than Tandem at the SAME mesh and order (0.0-0.7%).

**Important**: Tandem also has non-zero dip slip — it's a real 3D effect. But Tandem's
dip/strike ratio is ~10× smaller than ours at the same configuration.

### 19.6 Revised Understanding: NOT Just Convergence — Code Difference Exists

Since Tandem and our code use the same mesh, same order, and the DG face integrators
have been confirmed mathematically identical (Section 14), there MUST be a code-level
difference that we have NOT yet identified.

**What has been ruled out:**
- DG face integrators (consistency, symmetry, penalty) — IDENTICAL (agent comparison)
- Penalty formula and c_N_1 — IDENTICAL
- Normal computation — IDENTICAL (CalcOrtho vs cofactor)
- Fault basis (EmbedSlip, ProjectTraction) — IDENTICAL (Section 7)
- Far-field BC formula — IDENTICAL (Section 9)
- Friction law velocity direction — IDENTICAL (Section 15)

**What has NOT been checked:**
1. **Volume integrator** — MFEM's `ElasticityIntegrator` vs Tandem's `assembleVolume` kernel.
   The volume integral ∫ σ(u):ε(v) dx creates K_xz coupling through the λ term. If the
   quadrature order or Jacobian handling differs, K_vol differs → different u_z.

2. **Boundary face integrator** — MFEM uses `w = ip.weight` (no 1/2 factor) on boundary
   faces. If Tandem uses a different weighting, the boundary stiffness contribution differs.

3. **Quadrature order for volume terms** — MFEM defaults to 2p for volume integrals.
   Tandem may use a different rule (e.g., 2p+1 or exact integration).

4. **Additional stabilization** — Tandem may have a cross-component penalty or stabilization
   that we lack, reducing the K_xz / K_zz ratio.

### 19.7 Additional Finding: p=4 Earthquake Did Not Fully Propagate

The p=4 4000m run shows:
- Nucleation (x2=-24): slip_s = 4.59 m — full earthquake ✓
- Near-nuc (x2=-16): slip_s = 4.92 m — full earthquake ✓
- Center depth (x2=0): slip_s = **0.14 m** — earthquake did NOT reach center ✗
- Far depth (x2=16): slip_s = 2.98 m — partial earthquake

The earthquake propagated from the nucleation zone but did NOT reach the center of the
fault. This is a separate issue from the dip offset — possibly related to the coarser
mesh (4000m vs 1000m) providing insufficient resolution to sustain rupture propagation.

### 19.8 Next Steps: Find the Missing Difference

**Priority 1**: Compare volume integrator (MFEM vs Tandem) — quadrature order, Jacobian
handling, cross-component structure.

**Priority 2**: Compare boundary face integrator — weighting, penalty at boundaries.

**Priority 3**: Check if Tandem has any additional stabilization or post-processing that
reduces cross-component coupling.

## 20. Missed Test from v47: p=1 1000m with Tight BLR Tolerance

### 20.1 Gap Identified

v47a ran p=1 on 1000m with ×3 penalty and MUMPS-BLR at default tolerance (1e-10).
Result: V_nuc decayed from 0.01 → 5e-9, no earthquake. This was attributed to the
×3 penalty making the effective elastic stiffness exceed k_critical.

However, the v47 investigation plan (Section 8.1.2) explicitly called for:

> "p=1 with ×3 penalty + exact MUMPS: Run p=1 (not p=2) with exact MUMPS
> (smaller system, should fit in memory). If V_nuc still decays → confirms the
> issue is formulation, not solver. If V_nuc grows → BLR IS the issue at p=1."

**This test was never executed.** The tighter-tolerance tests (v47c at BLR 1e-14,
v47d exact MUMPS) were done for **p=2 only**. The investigation pivoted to the
multi-DOF blowup pattern (nbf > 1 at p≥2) and moved on to v48+.

### 20.2 Why This Matters

At p=2, tighter BLR tolerance made the instability **worse** (v47d: BLR 1e-14
produced larger τ_max than v47b: BLR 1e-10). But p=1 has a fundamentally different
failure mode (V decay, not blowup) — the BLR interaction could go the other way.

The v47a V-decay diagnosis assumed the cause was "×3 penalty → k_elastic > k_critical"
(Section 7.4). But BLR at 1e-10 introduces solver error that contaminates the
displacement field. At p=1 with nbf=1 (no multi-DOF issue), the only variables are:
1. Penalty magnitude (×3 — correct, matching Tandem)
2. Solver accuracy (BLR 1e-10 — never tested tighter for p=1)

If p=1 with BLR 1e-12 nucleates successfully, it would mean:
- The v47a V-decay was a BLR artifact, NOT a penalty stiffness issue
- p=1 with ×3 penalty is viable at tight tolerance
- The dip contamination analysis (Section 19) needs revisiting at p=1

### 20.3 Test: v51d — p=1, 1000m, BLR 1e-12

**Configuration**: p=1, 1000m mesh (bp5_tandem.msh), IP method, MUMPS-BLR 1e-12,
400 ranks, 48hr. All v50 defaults (CFL-aware dt, V-guard ON).

**Job script**: `jobs/bp5/bp5_v51d_1000m_ip_p1_blr12.sbatch`

**Success criterion**: V_nuc grows past 0.01 and earthquake nucleates.
**Failure criterion**: V_nuc decays as in v47a → confirms penalty stiffness is the cause.

---

## 21. Strategy 1 Result: Stress-Only Traction (v51f/f2 — FAILED)

Tested removing the penalty correction from traction recovery:
- **v51f (p=4)**: BLOWUP — V = 10 m/s in 35 steps
- **v51f2 (p=2)**: BLOWUP — V = 46 m/s in 397 steps

The penalty correction is ESSENTIAL for stability. Cannot be removed.
The dip contamination is in {σ·n} (from u_z), not in the penalty correction.
