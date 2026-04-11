# BP5 Debug v51: Dip-Component Offset — Systematic Cross-Component Coupling Error

**Date**: 2026-03-24
**Status**: ELASTIC σ_n DISPROVED AS FIX. Fair same-time comparison shows elastic σ_n makes dip 2-3× WORSE (not better). Initial "improvement" was artifact of comparing different time points. All K/f components confirmed identical to Tandem. Root cause of 10-30× excess dip at p=2 remains unidentified.
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
- ✅ v51d: p=1 1000m BLR 1e-12 → V decays, penalty stiffness confirmed
- ✅ v51e: p=1 500m BLR 1e-12 → V decays, p=1 dead end (polynomial order, not h)

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
| **v51d** | **p=1 1000m BLR 1e-12: V decays — penalty stiffness CONFIRMED (Section 20.4)** | **COMPLETED** |
| **v51e** | **p=1 500m BLR 1e-12: V decays — p=1 dead end, p≥2 required (Section 20.5)** | **COMPLETED** |
| v51 | **IP penalty ×3 first-principles verification: proven correct (Section 21)** | **CONFIRMED** |
| **v51** | **Benchmark conformance audit: exp cap removed, σ_n constant confirmed (Section 24)** | **CODE CHANGE** |
| v51a | Zero-dip-traction (2hr dev queue) — failed: walltime + V-guard rejection flood | **FAILED (infra)** |
| **v51f** | **Zero-dip-traction PRODUCTION (48hr, 400 ranks) — benchmark-correct per Eq. 15b** | **PLANNED** |
| v51f/f2 | Strategy 1 (stress-only traction): BLOWUP at p=2 and p=4 | **DISPROVED** |
| **v51c** | **DEFINITIVE: u_z = 18.5% of u_x in the elastic solution itself** | **ROOT CAUSE** |
| v51c | Source: K-matrix cross-coupling from non-fault interior face DG terms | **ANALYZED** |
| v51c | f_z = 0 confirmed on fault faces → u_z comes from K_zx coupling | **CONFIRMED** |
| v51c | Initial conclusion: DG discretization error scaling with p | HYPOTHESIZED |
| **v51c+** | **DISPROVED: p=4 on 4000m has 10× worse ratio than Tandem p=4 on same mesh** | **CODE DIFF EXISTS** |
| v51c+ | Tandem ALSO has non-zero dip slip — it's a real 3D effect, but our ratio is 10× too high | **KEY FINDING** |
| v51c+ | DG face integrators confirmed identical → difference must be in volume/BC/stabilization | **NARROWED** |
| v51c+ | Next: compare volume integrator, boundary face treatment, quadrature orders | **PLANNED** |
| **v51c++** | **Volume integrator + boundary face: ALL MATCH. Full K audit complete.** | **CONFIRMED** |
| v51c++ | Every K component, every f component, every pipeline element checked — ALL IDENTICAL | **PARADOX** |
| v51c++ | K identical + f identical yet u_z differs 10× → remaining: σ_n, Tandem pipeline, assembly bug | **NARROWED** |
| **v51c+++** | **RESOLUTION: Tandem uses elastic σ_n = σ_n_pre - T_n_elastic; we use constant 25 MPa** | **ROOT CAUSE** |
| v51c+++ | Tandem ODE pipeline: clean, no filtering — NOT the cause | **RULED OUT** |
| v51c+++ | Assembly pipeline: clean, no DOF bugs — NOT the cause | **RULED OUT** |
| v51c+++ | Elastic σ_n provides self-consistent feedback suppressing spurious dip accumulation | **MECHANISM** |
| v51c+++ | FaultBasis::NormalStress() EXISTS in our code but is NEVER CALLED | **FOUND** |
| v51c+++ | Implementation plan: expand traction to 3 components, pass σ_n_eff to friction | **PLANNED** |
| **v51d** | **p=1 1000m BLR 1e-12: V decays 0.010→0.0003 — confirms ×3 penalty too stiff at p=1** | **CONFIRMED** |
| v51d | Tighter BLR (1e-12 vs 1e-10) does NOT help — problem is penalty magnitude, not solver | **RULED OUT** |
| v51d_esn | p=2 1000m elastic σ_n: initial comparison INVALID (different time points) | **RETRACTED** |
| v51e_esn | p=4 4000m elastic σ_n: initial comparison INVALID (different time points) | **RETRACTED** |
| **v51d/e+** | **FAIR same-time comparison: elastic σ_n is 2-3× WORSE than constant σ_n** | **✗ DISPROVED** |
| v51d/e+ | Elastic T_n has same DG contamination as T_dip → amplifies error, not suppresses | **MECHANISM** |

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

### 20.4 Result: V Decays — Penalty Stiffness CONFIRMED

**Job**: `bp5_v51d_ip_p1_blr12_7610045.out`

```
Step   1: V = 0.01015  ↑
Step  15: V = 0.01058  ← peak
Step  50: V = 0.00580  ↓
Step 150: V = 0.00160  ↓
Step 300: V = 0.00033  ↓  (dt = 1.16s — entered interseismic regime)
```

V peaked at step 15 (V=0.0106) then monotonically decayed to 0.0003 by step 300.
Same pattern as v47a. **Tighter BLR tolerance (1e-12 vs 1e-10) makes no difference.**

**Conclusion**: The v47a nucleation failure at p=1 is caused by the ×3 penalty being
too stiff, NOT by BLR solver accuracy. The corrected penalty (which IS mathematically
correct — see Section 21) creates an effective elastic stiffness that overwhelms the
friction weakening rate at p=1 with DOF spacing = 1000m/1 = 1000m.

**Why p=2 works**: At p=2, DOF spacing = 1000m/2 = 500m. More DOFs per nucleation zone
(24 vs 12 at p=1) provide better resolution of the nucleation instability. The penalty
per face is stronger at p=2 (c_N_1=8/3 vs 1.0), but the finer DOF spacing compensates.

**Resolution**: p=1 on 1000m mesh with the correct ×3 penalty is **under-resolved for
nucleation**. This is NOT a bug — it's a resolution limitation.

This test closes the gap identified in v47's investigation plan.

### 20.5 Result v51e: p=1 500m — ALSO Cannot Nucleate

**Job**: `bp5_v51e_500m_ip_p1_7610108.out`

```
Step   1: V = 0.01011  ↑
Step  50: V = 0.01095  ← peak
Step 100: V = 0.01006  ↓
Step 150: V = 0.00865  ↓
Step 188: V = 0.00763  ↓  (walltime, still decaying)
```

Same pattern as v51d (1000m). V peaked at step 50 (V=0.0110) then monotonically
decayed. Slightly better than 1000m (peak 0.0110 vs 0.0106, slower decay rate)
but outcome is identical: **nucleation fails**.

**Key finding: halving h does NOT fix p=1 nucleation.** Both p=1 h=500m and p=2
h=1000m have the same DOF spacing (500m), but only p=2 nucleates. This proves the
issue is **polynomial order**, not DOF spacing or mesh resolution.

The linear polynomial (p=1) cannot represent the displacement field accurately
enough on tet meshes — the DG jump residual remains large regardless of h, keeping
the effective elastic stiffness above k_crit = 5.77 MPa/m.

**Comparison of p=1 attempts:**

| Run | h | DOF spacing | BLR tol | V peak | Outcome |
|-----|---|-------------|---------|--------|---------|
| v47a | 1000m | 1000m | 1e-10 | 0.0106 | Decays |
| v51d | 1000m | 1000m | 1e-12 | 0.0106 | Decays |
| **v51e** | **500m** | **500m** | **1e-12** | **0.0110** | **Decays** |
| v50 (p=2) | 1000m | 500m | 1e-12 | grows | **Nucleates** ✓ |

**Conclusion: p=1 with correct ×3 penalty is a dead end for BP5 nucleation.**
The minimum viable configuration is p=2 on 1000m. This is consistent with the
SCEC benchmark suggestion of Δz = Δh/N = 1000m (implying N≥2 for DG methods).

---

## 21. IP Penalty ×3: First-Principles Verification

### 21.1 Motivation

The ×3 penalty correction (v47) is load-bearing for matching Tandem's formulation.
v45 previously reverted it, claiming it was wrong. v47 re-applied it and showed it
breaks p=1 nucleation. Given its importance — and the fact that it changes nucleation
behavior — this section provides a complete first-principles derivation.

### 21.2 MFEM Reference Element Measures

**Reference tetrahedron** (vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1)):
- V_ref = 1/6 = 1/D!

**Reference triangle** (vertices (0,0), (1,0), (0,1)):
- A_ref = 1/2 = 1/(D-1)!

### 21.3 What MFEM Functions Return

**`Trans.Elem1->Weight()`** — `densemat.cpp:553-558`: returns `det(J_elem)` for
square Jacobians. For a flat tet:

```
V_phys = ∫_ref |det(J)| dξ = |det(J)| × V_ref = det(J) / 6
→  Weight() = det(J_elem) = 6 × V_phys
```

**`CalcOrtho(Trans.Jacobian(), nor)`** — `densemat.cpp:2692-2716`: computes cross
product of the two columns of the 3×2 face Jacobian. For a flat triangle:

```
A_phys = ∫_ref |nor| dξ = |nor| × A_ref = |nor| / 2
→  nl_q = |nor| = 2 × A_phys
```

### 21.4 Concrete Verification: Right-Angle Tet

Use a right-angle tet: v0=(0,0,0), v1=(a,0,0), v2=(0,b,0), v3=(0,0,c).

| Quantity | Formula | Value |
|----------|---------|-------|
| V_phys | abc/6 | abc/6 |
| A_phys (base face z=0) | ab/2 | ab/2 |
| Weight() = det(J) | abc | = 6 × V_phys ✓ |
| \|CalcOrtho\| = nl_q | \|e1 × e2\| = ab | = 2 × A_phys ✓ |

**The ratio nl_q / Weight():**

```
nl_q / Weight() = ab / abc = 1/c
Physical A/V = (ab/2) / (abc/6) = 3/c

→  nl_q / Weight() = (1/3) × (A/V)     ← WRONG by factor 3
→  dim × nl_q / Weight() = 3/c = A/V   ← CORRECT
```

The factor of 3 = D arises from the ratio of reference measures:

```
nl_q / Weight() = (2 × A_phys) / (6 × V_phys) = A / (3V)
The "3" = D!/(D-1)! = D = 3
```

### 21.5 Tandem Verification

Tandem precomputes physical area and volume directly:

**Volume** (`DGCurvilinearCommon.cpp:55-59`):
```cpp
volume = Σ_q w_q × |det(J_q)|    // w_q sums to V_ref = 1/6
       = (1/6) × |det(J)| = V_phys     ✓
```

**Area** (`DGCurvilinearCommon.cpp:90-94`):
```cpp
area = Σ_q w_q × |normal_q|      // w_q sums to A_ref = 1/2
     = (1/2) × 2A = A_phys             ✓
```

Where `normal_q` is computed via the cofactor formula (`Curvilinear.cpp:234-244`):
```cpp
normal = |det(J_elem)| × J_elem^{-T} × N_ref
```

Tandem's penalty uses `area/volume = A_phys/V_phys` directly — no reference scaling.

### 21.6 Penalty Formula Comparison

**Tandem** (`Elasticity.cpp:188-191`):
```
p(side) = (D+1) × c_N_1 × (area / volume) × (c₁²/c₀)
        = (D+1) × c_N_1 × (A/V) × (c₁²/c₀)
```

**Our code WITH ×3** (`dg_elasticity_ip_penalty_integrator.hpp:122`):
```
p0 = (D+1) × c_N_1 × (dim × nl_q / Weight()) × (c₁²/c₀)
   = (D+1) × c_N_1 × (3 × 2A / 6V) × (c₁²/c₀)
   = (D+1) × c_N_1 × (A/V) × (c₁²/c₀)        ← MATCHES TANDEM ✓
```

**Our code WITHOUT ×3** (pre-v47):
```
p0 = (D+1) × c_N_1 × (nl_q / Weight()) × (c₁²/c₀)
   = (D+1) × c_N_1 × (A/(3V)) × (c₁²/c₀)     ← 1/3 OF TANDEM ✗
```

### 21.7 Assembly-Level Comparison

Both codes integrate: `∫_F η_F × [u]·[v] dS ≈ Σ_q η_F × w_q × nl_q × φ_i × φ_j`

**Tandem** (per face, flat):
```
penalty × Σ_q w_q × nl_q × φ_i × φ_j
= (D+1)·c_N_1·(A/V)·(c₁²/c₀) × A_phys × ⟨φ_i,φ_j⟩_ref
```

**Our code WITH ×3** (per quad point, flat face — nl_q = const):
```
Σ_q [(D+1)·c_N_1·(3·nl_q/W)·(c₁²/c₀)] × w_q × nl_q × φ_i × φ_j
= (D+1)·c_N_1·(A/V)·(c₁²/c₀) × A_phys × ⟨φ_i,φ_j⟩_ref     ← IDENTICAL ✓
```

### 21.8 Why the Consistency Term Does NOT Need ×3

In MFEM's `DGElasticityIntegrator::AssembleBlock` (`bilininteg.cpp:4002-4014`):
```cpp
w1 = ip.weight / (2 × Weight())        // = w / (2·det(J))
nL1 = w1 × λ × nor                      // = λ·w/(2·det(J)) × nor
dshape_ps = dshape_ref × adj(J)          // = det(J) × ∂φ/∂x_phys
```

Product: `dshape_ps(j,m) × nL1(i)`:
```
= [det(J) × ∂φ_j/∂x_m] × [λ × w / (2·det(J)) × nor_i]
= λ × (w/2) × (∂φ_j/∂x_m) × nor_i
```

**det(J) cancels between adj(J) and 1/Weight().** The consistency term naturally
produces the correct physical traction — no reference element scaling remains.

Only the penalty has the ×3 issue because it explicitly uses the A/V ratio (nl_q/Weight),
while the consistency term uses the gradient×normal product where det(J) self-cancels.

### 21.9 Numerical Values for BP5 at p=1, h=1000m

Using λ = μ = 32.04 GPa, regular tet h=1000m:

```
c₀ = 2μ = 64.08 GPa
c₁ = 3λ + 2μ = 160.2 GPa
c₁²/c₀ = 400.5 GPa
c_N_1 = p(p+D-1)/D = 1×3/3 = 1.0
(D+1) = 4
A/V ≈ 3.674e-3 m⁻¹
```

| Penalty version | p_side (GPa/m) | penalty_ip (GPa/m) | k_spring/k_phys |
|----------------|----------------|---------------------|------------------|
| WITH ×3 (correct) | 5.87 | 2.94 | **37×** |
| WITHOUT ×3 (1/3) | 1.96 | 0.98 | 12× |
| Tandem | 5.87 | 2.94 | 37× |

### 21.10 Why ×3 Breaks p=1 Nucleation (Despite Being Correct)

The nucleation criterion (spring-slider analog): `k_elastic < k_crit = σ_n·b/D_c = 5.77 MPa/m`

v47 Section 7.4 measured k_elastic at the nucleation station (x2=-24km, x3=10km):

| Time | Slip (mm) | k_elastic (MPa/m) | vs k_crit | V trend |
|------|-----------|-------------------|-----------|---------|
| 0.44s | 4.6 | 3.9 | BELOW | Growing |
| 0.68s | 7.0 | 5.0 | BELOW | Growing |
| 0.90s | 9.3 | **5.9** | **ABOVE** | **Decaying** |
| 1.83s | 17.4 | 7.9 | ABOVE | Decaying |

At ~9mm slip (t≈0.9s), k_elastic crosses k_crit and V decays irreversibly.

Tandem runs at p=4 where the polynomial better resolves the displacement field —
slip spreads across a wider area, keeping the effective stiffness below k_crit.
At p=1, the DG solution is too coarse: the correct penalty creates a fault coupling
that is too stiff for nucleation at h=1000m.

The 1/3 penalty (pre-v47) was a **compensating error**: wrong penalty value, but
accidentally soft enough for p=1 nucleation. This gave fortuitously good results
at p=1 that would not converge correctly at higher p.

### 21.11 Conclusion

The ×3 factor is **mathematically proven correct**:
1. `Weight() = 6V`, `|CalcOrtho| = 2A` → `nl/W = A/(3V)` → need ×3 for A/V
2. Tandem uses physical A/V directly → our ×3 matches exactly
3. Assembly-level total penalty is identical between codes
4. Consistency term self-cancels det(J) → no ×3 needed there
5. v45's claim that ×3 was wrong was **incorrect** — the consistency-penalty balance
   is restored (not broken) by the ×3

The v47a p=1 nucleation failure with ×3 is a **resolution issue** (p=1 too coarse),
not a penalty error. The v51d test (BLR 1e-12) will determine whether solver accuracy
is a confounding factor.

---

## 22. Strategy 1 Result: Stress-Only Traction (v51f/f2 — FAILED)

Tested removing the penalty correction from traction recovery:
- **v51f (p=4)**: BLOWUP — V = 10 m/s in 35 steps
- **v51f2 (p=2)**: BLOWUP — V = 46 m/s in 397 steps

The penalty correction is ESSENTIAL for stability. Cannot be removed.
The dip contamination is in {σ·n} (from u_z), not in the penalty correction.

## 22. Volume Integrator + Boundary Face Comparison: ALL MATCH

### 22.1 Volume Integrator Comparison

Exhaustive comparison of MFEM's `ElasticityIntegrator::AssembleElementMatrix`
(bilininteg.cpp) vs Tandem's `assembleVolume` kernel (elasticity.py):

| Aspect | MFEM | Tandem | Match? |
|--------|------|--------|--------|
| Formula | ∫ (λ div(u)·div(v) + 2μ ε:ε) dx | Same | ✓ |
| Quadrature order | 2p-2 | 2p+1 | Both exact for linear tets ✓ |
| Jacobian handling | adj(J)/det(J) cancellation | J⁻¹ directly | Equivalent ✓ |
| Material evaluation | Direct at quad points | L2 projection to quad points | Same for constant material ✓ |
| Cross-component coupling | λ creates K_xz through div(u)·div(v) | Same | ✓ |

For linear tets: both quadrature orders (2p-2 and 2p+1) give EXACT integration of the
degree 2p-2 integrand. No numerical difference.

### 22.2 Boundary Face Comparison

| Aspect | MFEM | Tandem | Match? |
|--------|------|--------|--------|
| Consistency term | w (no 1/2 factor) = double interior | c00=-1.0 (double of -0.5) | ✓ |
| Symmetry term | α * w (no 1/2 factor) | c10=ε*1.0 (double of ε*0.5) | ✓ |
| Penalty | kappa * |nor|² * wLM (no 1/2) | 2 * penalty (double of interior) | ✓ |
| Quadrature order | 2p | 2p+1 | Both exact for linear tets ✓ |

### 22.3 Complete K-Matrix Audit Summary

**Every component of K has been checked and confirmed IDENTICAL:**

| K component | Section | Result |
|------------|---------|--------|
| Interior face consistency+symmetry | §14 | IDENTICAL |
| Interior face penalty | §14 | IDENTICAL |
| Volume integrator | §22 | IDENTICAL |
| Boundary face consistency+symmetry | §22 | IDENTICAL |
| Boundary face penalty | §22 | IDENTICAL |

**Every component of f has been checked:**

| f component | Section | Result |
|------------|---------|--------|
| Fault face slip RHS (penalty term) | §19 | f_z = 0 CONFIRMED |
| Fault face slip RHS (symmetry term) | §19 | f_z = 0 CONFIRMED |
| Fault basis EmbedSlip | §7 | IDENTICAL to Tandem |
| Far-field Dirichlet BCs | §9 | IDENTICAL |

**Every other component checked:**

| Component | Section | Result |
|----------|---------|--------|
| Friction law velocity direction | §15 | IDENTICAL |
| FaultBasis (ProjectTraction) | §7 | IDENTICAL |
| Normal computation (CalcOrtho) | §14 | IDENTICAL |
| Penalty coefficient c_N_1 | §14 | IDENTICAL |

### 22.4 The Paradox

K is identical. f is identical. Yet u_z differs by 10× between our code and Tandem
at the same mesh and polynomial order. The remaining unexplored possibilities:

1. **Something in Tandem we haven't read** — an additional post-solve processing step,
   a filter on the traction, or a different ODE coupling strategy that we missed.

2. **The elastic σ_n difference** (Section 15) — Tandem uses elastic normal stress
   from the displacement field; we use constant σ_n = 25 MPa. This is the ONLY
   confirmed implementation difference. While the direct effect seems small (~4%
   change in σ_n during coseismic), it could affect how the friction law responds
   to small dip traction perturbations over many earthquake cycles.

3. **Tandem's matrix-free operator application** vs our assembled sparse matrix +
   MUMPS direct solve — in principle equivalent, but floating-point accumulation
   order differs. However, this should only cause machine-epsilon differences.

4. **A subtle assembly or pipeline bug** in our code that doesn't show up when
   comparing individual integrators — e.g., DOF indexing, element connectivity
   ordering, or how shared face contributions are accumulated across MPI ranks.

### 22.5 Recommended Next Steps

**Option A**: Dig deeper into Tandem's ODE evaluation pipeline — read how Tandem's
`SeasQDOperator::evaluate()` calls the traction computation and whether any
post-processing happens between the elastic solve and the friction law.

**Option B**: Implement the elastic σ_n — the ONLY confirmed difference. Test whether
it reduces the dip/strike ratio.

**Option C**: Direct numerical test — on a small mesh (4-8 elements), compute K*e_x
(K applied to a pure x-displacement vector) and check the z-component. Compare with
Tandem on the same small mesh. This bypasses the integrator comparison and directly
tests the assembled K matrix.

---

## 23. RESOLUTION: Elastic Normal Stress σ_n — The Missing Implementation

### 23.1 The Investigation

All three options from Section 22.5 were investigated simultaneously:

| Investigation | Finding | Explains 10×? |
|--------------|---------|---------------|
| **A: Tandem ODE pipeline** | Clean. No filtering/projection/post-processing. Same structure as ours. | No |
| **B: Elastic σ_n** | **CRITICAL DIFFERENCE CONFIRMED** — Tandem uses elastic σ_n, we use constant | **YES** |
| **C: Assembly pipeline** | Clean. No DOF bugs, no BC contamination, no residual between steps. | No |

### 23.2 What Tandem Does (That We Don't)

**Tandem passes 3 traction components to the friction law:**

```
traction(node, 0) = T · n_fault    (normal traction)
traction(node, 1) = T · dip_vec    (dip traction)
traction(node, 2) = T · strike_vec (strike traction)
```

**Tandem uses elastic σ_n in the friction law:**

From `RateAndState.h:159` and `DieterichRuinaAgeing.h:86`:
```cpp
auto sn = t_mat(node, 0);                     // elastic normal traction
double snAbs = -sn + p_[index].get<SnPre>();   // σ_n_eff = σ_n_pre - T_n_elastic
```

**Our code discards the normal traction entirely:**

From `rate_state_fault.hpp:436`:
```cpp
dr_friction_->SolveSlipRateVectorPsi(
    tau_vec, psi, sigma_n_bp5_, eta, a, V_vec);  // sigma_n = constant 25 MPa
```

Our `FaultBasis::NormalStress()` method EXISTS (fault_basis.hpp line 323) but is
**NEVER CALLED** in the friction pipeline.

### 23.3 How This Explains the 10× Dip Ratio Difference

**The mechanism:**

The DG elastic solution on tet mesh produces spurious cross-component displacements
(u_z ≠ 0 for pure strike-slip, measured at 18.5% of u_x). This creates non-zero
traction in ALL three directions: T_strike, T_dip, AND T_normal. All three components
are correlated — they arise from the same discretization error on tet faces.

**In Tandem** (3-component traction + elastic σ_n):
1. T_n_elastic is computed from the SAME stress field as T_dip and T_strike
2. σ_n_eff = 25 MPa - T_n_elastic varies per node, per time step
3. The friction strength F = σ_n_eff * [f₀ + a·ln(V/V₀) + b·ln(V₀ψ/L)] adjusts accordingly
4. Through the nonlinear RSF law (V ~ sinh(τ/(a·σ_n))), the σ_n correction modifies
   the velocity magnitude and direction
5. Because T_n is correlated with T_dip, the σ_n feedback provides a **self-consistent
   correction** that suppresses spurious dip accumulation
6. The result: dip/strike ratio 0.0-0.7% (small, dominated by real 3D effects)

**In our code** (2-component traction + constant σ_n):
1. T_n_elastic is NEVER COMPUTED (discarded in ProjectTraction)
2. σ_n = 25 MPa always — no feedback from the elastic solution
3. The friction strength is computed with the wrong normal stress
4. The spurious T_dip propagates unchecked through the friction law
5. No self-consistent feedback to suppress dip accumulation
6. The result: dip/strike ratio 0.4-8.1% (10× worse than Tandem)

**Why even a small σ_n correction matters:**

The RSF law is exponentially sensitive to stress/strength ratio:
```
V = 2·V₀ · sinh(|τ|/(a·σ_n)) · exp(-ψ/a)
```

A 0.1% change in σ_n changes the argument of sinh by 0.1%, which during coseismic
(where the argument is large) shifts V by several percent. Over 30 seconds of
coseismic slip at V ~ 0.5 m/s, this accumulates to measurable dip slip differences.

More importantly, the σ_n correction is **correlated** with the dip contamination
(both come from the same 3D stress field). So the correction systematically opposes
the dip error, not randomly.

### 23.4 Comparison Table

| Aspect | Our Code | Tandem |
|--------|----------|--------|
| Traction output components | 2 (dip, strike) | 3 (normal, dip, strike) |
| σ_n in friction law | Constant 25 MPa | σ_n_pre - T_n_elastic (per node, per step) |
| Normal traction computation | Discarded | L2-projected same as tangential |
| FaultBasis::NormalStress() | EXISTS but never called | Equivalent: fault_basis_q[:,0,:] |
| Self-consistent feedback | ✗ None | ✓ T_n correlated with T_dip |
| Dip/strike ratio (p=4, 4000m) | 0.4-8.1% | 0.0-0.7% |

### 23.5 Tandem Code References

**Traction evaluation** — `elasticity_adapter.py:26-27`:
```python
traction['kp'] <= minv['lk'] * e_q_T['ql'] * w['q'] * nl_q['q'] *
                  traction_q['oq'] * fault_basis_q['opq']
```
Output shape: `(nbf, 3)` — columns are (normal, tangent1, tangent2) = (n, dip, strike).

**Friction law σ_n usage** — `DieterichRuinaAgeing.h:86`:
```cpp
double snAbs = -sn + p_[index].get<SnPre>();
// sn = traction(node, 0) = elastic T · n_fault
// SnPre = 25.0 MPa (from bp5.lua)
// snAbs = σ_n_effective (positive in compression)
```

**Our unused method** — `fault_basis.hpp:323`:
```cpp
real_t NormalStress(int face_idx, const real_t T_global[3]) const
{
    const auto &basis = face_bases_[face_idx];
    return -(T_global[0]*basis.normal[0] +
             T_global[1]*basis.normal[1] +
             T_global[2]*basis.normal[2]);
}
```
Sign convention: returns positive value for compressive normal stress (correct for RSF).

### 23.6 Implementation Plan

**Step 1**: Expand traction vector from `2 * num_fault_dofs_` to `3 * num_fault_dofs_`
- Layout: `[T_n_0, T_dip_0, T_strike_0, T_n_1, T_dip_1, T_strike_1, ...]`
- Modify ComputeTraction to also compute `T_n = NormalStress(fi, T_global)` per DOF

**Step 2**: Pass elastic σ_n to friction law
- In `rate_state_fault.hpp` ComputeRHS: extract T_n from traction vector
- Compute `sigma_n_eff = sigma_n_pre + T_n` (NormalStress returns positive for compression)
- Pass `sigma_n_eff` to `SolveSlipRateVectorPsi` instead of constant `sigma_n_bp5_`

**Step 3**: Add `--elastic-sigma-n` flag for A/B testing
- Default: OFF (constant σ_n, matching current behavior)
- ON: use elastic σ_n (matching Tandem)
- Run both on same mesh, compare dip/strike ratio

**Step 4**: Verify
- p=4 on 4000m with elastic σ_n: expect dip/strike ratio to drop from 0.4-8.1% to ~0.0-0.7%
- p=2 on 1000m with elastic σ_n: expect improvement but may still be higher than Tandem
  (because u_z contamination is larger at p=2)

### 23.7 What Was Ruled Out (Complete Audit)

After exhaustive investigation spanning Sections 7-22, every other hypothesis has been
eliminated:

- FaultBasis (EmbedSlip, ProjectTraction) — IDENTICAL (§7)
- Far-field BC — IDENTICAL (§9)
- Friction law velocity direction — IDENTICAL (§15)
- DG face integrators (consistency, symmetry, penalty) — IDENTICAL (§14, §22)
- Volume integrator — IDENTICAL (§22)
- Boundary face treatment — IDENTICAL (§22)
- Assembly pipeline, DOF ordering — CLEAN (§22)
- Tandem ODE pipeline — CLEAN, no filtering (§23)

**The elastic σ_n is the ONLY remaining implementation difference.**

---

## 24. Rate-and-State Benchmark Conformance Audit

### 24.1 Motivation

Line-by-line comparison of SCEC BP5 spec (SEAS_BP5.pdf) equations against our code
and Tandem, to verify all three are consistent (or document where they diverge).

### 24.2 Equation-by-Equation Comparison

| # | Benchmark | Our Code | Tandem | Status |
|---|-----------|----------|--------|--------|
| **Eq. 13** | f = a·asinh[(V/2V₀)·exp((f₀+b·ln(V₀θ/L))/a)] | f = a·asinh[(V/2V₀)·exp(ψ/a)] | Same | ✓ Identical (ψ = f₀+b·ln(V₀θ/L)) |
| **Eq. 12** | dθ/dt = 1 − Vθ/L | dpsi/dt = (bV₀/L)·[exp((f₀−ψ)/b) − V/V₀] | Same, **no exp cap** | **✓ Fixed** (exp cap removed) |
| **Eq. 10** | \|τ⁰+Δτ\| = σ_n·f + η·V (QD) | Brent solve on same equation | Same | ✓ Identical |
| **Eq. 11** | V direction = τ̂ | V = +(V/\|τ\|)·τ | V = **−**(V/\|τ\|)·τ | ✓ Self-consistent (sign conventions differ) |
| **Sect. 3** | **σ_n = constant** | **Constant 25 MPa** (default) | **−sn + sn_pre** (variable) | **Tandem deviates from spec** |
| **Eq. 16** | V = [V_init, V_zero] | V = [V_zero, V_init] (dip, strike) | Same | ✓ (component order matches geometry) |
| **Eq. 18** | θ(0) = L/V_init | ψ = f₀+b·ln(V₀/Vp) = f₀+b·ln(V₀/V_init) | Same | ✓ (Vp = V_init for BP5) |
| **Eq. 20** | τ⁰ = σ_n·f(V_init,θ) + η·V_init | Same | Same | ✓ Identical |
| **Eq. 23** | τ⁰_nuc = σ_n·f(V_i,θ) + δτ, δτ=ηV_i | sigma_n·f + η·V_nuc (dtf=0) | Same (dtf=0) | ✓ (δτ IS the η·V term, not added on top) |
| **Eq. 15** | Outside Ω_f: V₂=Vp, V₃=0 | dslip = (0, Vp, 0) | Same | ✓ Identical |
| **Table 1** | V_i = 0.03 m/s | V_nuc = 0.01 (Tandem default) | 0.01 | **Both deviate from SCEC** |

### 24.3 Findings and Code Changes

**Finding 1: AgingLawPsi exp cap removed** (this section)

The `AgingLawPsi::Rate()` and `RateDerivativeTheta()` had an `exp_arg_max_ = 20`
cap that limited `exp((f₀−ψ)/b)` to `exp(20) ≈ 4.85×10⁸`. This deviated from both
the benchmark (Eq. 12, no cap) and Tandem (no cap).

**Change**: Removed the cap from both methods in `friction/state_evolution.hpp`.
The adaptive RK45 stepper handles post-earthquake stiffness by reducing dt.

**Affected code**: `AgingLawPsi::Rate()` (line 191) and `RateDerivativeTheta()` (line 214).

**Risk**: During post-earthquake healing, `exp((f₀−ψ)/b)` can reach `exp(100+) ~ 10⁴³`,
making dpsi/dt ~ 10³⁵. The RK45 stepper must take very small steps (dt ~ 10⁻³⁵ s)
to resolve this. If the stepper's dt_min is too large, it may fail.

Current `dt_min = 1e-6` in `bp5_verification_full.cpp:1101`. This may need to be
reduced for the uncapped formulation. Monitor for "dt below minimum" errors.

**Finding 2: σ_n is correctly constant (matching benchmark)**

The SCEC spec explicitly states σ_n remains constant for identical materials.
Our default (constant 25 MPa) is more faithful to the spec than Tandem's elastic
σ_n feedback. The `--elastic-sigma-n` flag should NOT be used for benchmark-conforming
runs. (It may be useful for non-planar or bimaterial extensions.)

**Finding 3: Tandem slip rate sign is negated**

Tandem returns `V = −(V/|τ|)·τ` while the benchmark defines V parallel to τ (Eq. 11).
This is a self-consistent convention within Tandem (EmbedSlip and traction extraction
have matching signs). When comparing output plots, Tandem's slip has the opposite sign
from the benchmark's s_j = u_j(0⁺) − u_j(0⁻) definition.

**Finding 4: delta_tau_factor = 0 is correct**

The benchmark Eq. 23 defines δτ = ηV_i for QD. This δτ IS the radiation damping
contribution (replacing the +ηV term in Eq. 20), not an overstress on top of it.
With our code structure (general formula already includes +ηV), setting
`delta_tau_factor = 0` gives the correct equilibrium initialization.
Setting `delta_tau_factor = 1` would DOUBLE the radiation damping → incorrect.

---

## 25. Elastic σ_n Test Results: DRAMATIC Improvement (v51d_esn, v51e_esn)

### 25.1 Test Configuration

| Run | Mesh | Order | σ_n mode | Flag | Status |
|-----|------|-------|----------|------|--------|
| v51d_esn | 1000m | p=2 | **Elastic** | `--elastic-sigma-n` | Running (48hr) |
| v51e_esn | 4000m | p=4 | **Elastic** | `--elastic-sigma-n` | Running (48hr) |
| v50 prod | 1000m | p=2 | Constant 25 MPa | (default) | Complete (223 yr) |
| v50h | 4000m | p=4 | Constant 25 MPa | (default) | Complete (0.43 yr) |

### 25.2 CORRECTION: Initial Comparison Was INVALID

The initial comparison (below, struck through) compared v50 at t=223 years with v51 at
t=23 seconds. This was meaningless — dip accumulates during earthquakes (t≈50-80s), so
any run at t=23s trivially has near-zero dip ratio regardless of σ_n treatment.

~~p=2: dip ratio drops 0.3-5.4% → 0.00-0.03%~~ ← WRONG: different time points
~~p=4: dip ratio drops 0.4-8.1% → 0.00-0.08%~~ ← WRONG: different time points

### 25.3 FAIR Same-Time Comparison: Elastic σ_n is WORSE

When compared at the SAME time (t=30s for p=2, t=35s for p=4), the elastic σ_n
correction makes the dip deviation **2-3× worse**, not better.

**p=2 Nucleation station (strk-24dp+10) — same-time comparison:**

| t (s) | v50 slip_d | v51d_esn slip_d | Tandem p6 slip_d | v50 − Tandem | v51d − Tandem | Result |
|-------|-----------|----------------|-----------------|-------------|--------------|--------|
| 10 | 4.26e-6 | 1.29e-5 | -4.77e-7 | 4.74e-6 | 1.33e-5 | **✗ 2.8× worse** |
| 20 | 3.55e-5 | 9.67e-5 | 9.87e-7 | 3.45e-5 | 9.57e-5 | **✗ 2.8× worse** |
| 30 | 1.34e-4 | 3.29e-4 | 9.99e-6 | 1.24e-4 | 3.19e-4 | **✗ 2.6× worse** |

**p=4 Nucleation station — same-time comparison:**

| t (s) | v50h slip_d | v51e_esn slip_d | Tandem p4 slip_d | v50h − Tandem | v51e − Tandem | Result |
|-------|-----------|----------------|-----------------|-------------|--------------|--------|
| 20 | -2.85e-6 | -1.20e-5 | 3.88e-6 | -6.74e-6 | -1.59e-5 | **✗ 2.4× worse** |
| 30 | -1.18e-5 | -3.33e-5 | 2.51e-5 | -3.69e-5 | -5.84e-5 | **✗ 1.6× worse** |
| 35 | 1.59e-5 | -1.78e-5 | 6.67e-5 | -5.09e-5 | -8.45e-5 | **✗ 1.7× worse** |

### 25.4 Why Elastic σ_n Makes It Worse

The elastic T_n has the **same DG cross-component contamination** as T_dip (both
originate from the same tet mesh discretization error in {σ·n}). Adding contaminated
T_n into σ_n_eff introduces an ADDITIONAL error source into the friction law.

The "self-consistent feedback" hypothesis was wrong:
- Hypothesis: T_n is correlated with T_dip → σ_n correction opposes T_dip error
- Reality: T_n adds contaminated noise that AMPLIFIES the dip deviation

The feedback is self-consistent with the WRONG stress field (contaminated by DG
error on tets), so it reinforces the error rather than canceling it.

### 25.5 Implications

1. **Elastic σ_n is NOT the fix** for the dip offset — it makes things worse
2. The dip deviation exists at BOTH p=2 and p=4, at 10-30× above Tandem
3. The deviation starts during nucleation (t≈12s), well before the earthquake
4. The root cause remains unidentified after exhaustive audit of all K/f components

### 25.6 What's Left to Investigate

After ruling out:
- Face integrators (identical) ✓
- Volume integrator (identical) ✓
- Boundary face treatment (identical) ✓
- FaultBasis (identical) ✓
- Friction law (identical) ✓
- Elastic σ_n (makes it worse) ✓
- Assembly pipeline (clean) ✓
- Tandem ODE pipeline (clean) ✓

Remaining possibilities:
1. **Tandem's matrix-free evaluation** vs our assembled sparse matrix — floating-point
   accumulation order differs, could systematically reduce cross-component error
2. **Tandem's element Jacobian caching** — Tandem precomputes and caches Jacobians;
   our code re-evaluates per quadrature point. Numerical differences accumulate.
3. **A subtle numerical difference** in how MFEM computes CalcAdjugate/CalcInverse
   vs Tandem's Jacobian routines — tiny per-element differences that accumulate
   across thousands of elements into a measurable u_z
4. **The dip deviation IS inherent** to our DG implementation at p=2 on 1000m tets,
   and Tandem at p=4/p=6 on 4000m simply has less because the discretization error
   is smaller at higher effective resolution
| Far-nuc surface | 5.42% | **0.01%** | 0.61% | **540×** better |

**p=4 (4000m): v50h (constant σ_n) vs v51e_esn (elastic σ_n)**

| Station | v50h |d/s| | v51e |d/s| | Tandem p4 |d/s| | Improvement |
|---------|-------------|-------------|-----------------|-------------|
| Nucleation | 0.39% | **0.00%** | 0.19% | **eliminated** |
| Near-nuc depth | 0.52% | **0.08%** | 0.00% | **6.5×** better |
| Center depth | 2.51% | **0.01%** | 0.14% | **250×** better |
| Center surface | 6.78% | **0.00%** | 0.23% | **eliminated** |
| Far depth | 2.36% | **0.00%** | 0.24% | **eliminated** |
| Far surface | 3.18% | **0.00%** | 0.66% | **eliminated** |
| Far-nuc surface | 8.09% | **0.00%** | 0.67% | **eliminated** |

### 25.3 Assessment

The elastic σ_n correction produces **dramatic improvement** at both polynomial orders:

- p=2: dip/strike ratio drops from 0.3-5.4% to **0.00-0.43%** — matching or BETTER than Tandem's 0.02-0.67%
- p=4: dip/strike ratio drops from 0.4-8.1% to **0.00-0.08%** — ALL stations near zero

**Every station at both orders shows improvement** (except Near-nuc depth at p=2, which shows slight noise at the early time point — likely meaningless).

### 25.4 Caveats

Both v51 runs are still in the **early nucleation phase** (t=23-31s). The dip/strike ratios at this stage are near zero because total slip is tiny. The critical test is whether the improvement **persists through the earthquake and into interseismic**.

However, the comparison is meaningful because:
1. v50 at the same early times ALSO had tiny slip, yet already showed measurable dip ratios
2. The absolute dip slip in v51 is orders of magnitude smaller than v50 at comparable times
3. The elastic σ_n eliminates the **seed** of dip contamination, which should prevent long-term accumulation

The runs are continuing on the 48-hour normal queue. Full earthquake cycle comparison will be available when they complete.

### 25.5 Mechanism Confirmed

The elastic σ_n provides **self-consistent feedback** (Section 23.3):
- The DG stress field produces correlated T_n, T_dip, T_strike perturbations
- With constant σ_n: T_dip drives V_dip unchecked → accumulates dip slip
- With elastic σ_n: T_n modifies σ_n_eff → friction strength adjusts → V_dip suppressed

This is NOT a workaround or filter — it's the **physically correct** implementation that Tandem uses. The SCEC spec assumes constant σ_n for the analytical initial condition (Eq. 20), but the elastic solve naturally produces σ_n perturbations that must be fed back for self-consistency.

### 25.6 Note on Benchmark Conformance

Section 24.2 noted that the SCEC spec defines σ_n as constant. Tandem's elastic σ_n technically deviates from the literal spec. However:
1. Tandem IS the reference implementation for BP5
2. The elastic σ_n perturbation is small (~0.1% of 25 MPa during coseismic)
3. Without elastic σ_n, the DG discretization error produces 10× worse dip contamination
4. The spec's constant-σ_n assumption is for identical materials on a planar fault — the DG cross-component coupling is an implementation artifact, not a physical feature

**Recommendation**: Use `--elastic-sigma-n` as the **default** for all production runs. The constant-σ_n mode should only be used for strict benchmark conformance testing.
