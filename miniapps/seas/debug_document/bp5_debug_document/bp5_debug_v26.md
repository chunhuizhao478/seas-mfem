# BP5 Debug v26: Traction Instability at Specific DOFs — Open Investigation

**Date**: 2026-03-14
**Status**: Instability persists on new mesh, root cause at specific DOFs unknown
**Previous**: v25 (mesh quality fix + IP-style traction)

---

## 1. Current State

We have two fixes:
- **H26**: Traction formula — IP-style penalty matching Tandem (no BR2 lifting)
- **H27**: Mesh quality — uniform fault resolution matching Tandem

**Result on v25 Test 1 (uniform fault, no earthquake, new mesh + IP traction):**
- Station DOFs: all healthy (V/Vp=0.987 at z=22km, t=1.95yr)
- But Vmax=568 m/s at a non-station DOF — single-DOF runaway
- ODE integrator stuck at t=1.95yr rejecting steps
- Same instability as on old mesh, same time (~2yr)

---

## 2. The Core Problem

A small number of fault DOFs develop unbounded average stress `{{σ·n̂}}`
over time, even in a uniform test (V=Vp everywhere, no earthquake).

| Configuration | Stable? | VS lockup? | Notes |
|---------------|---------|-----------|-------|
| BR2 traction + old mesh | ✓ | YES (V/Vp→0.078) | BR2 correction masks instability + biases traction |
| IP traction + old mesh | ✗ blowup | N/A | 5 DOFs at y≈28-29km, z≈5-6km hit 1 GPa |
| IP traction + new mesh | ✗ blowup | N/A | Different DOF(s), Vmax→568 at t≈1.95yr |
| **BR2 traction + new mesh** | **? untested** | **?** | **Key test — may reduce VS lockup bias** |

The instability is:
- **Independent of mesh quality** (occurs on both old and new mesh)
- **Independent of mesh resolution** (both meshes are ~1km on fault)
- **Localized to specific DOFs** (all station DOFs remain healthy)
- **Time-dependent** (develops over ~2 years from perfectly uniform initial conditions)
- **Suppressed by BR2 correction** (which provides negative feedback on DG jumps)

---

## 3. What We Don't Know

### 3.1 Why these specific DOFs?

On the old mesh, the 5 blowup DOFs were at y≈28-29 km, z≈5-6 km (near nuc2
boundary). On the new mesh, the blowup DOF location is unknown (the TRACTION
BLOWUP diagnostic with coordinates was not in this build — it was added to
the old mesh code only).

**Open question**: Are the unstable DOFs always at the same geometric location?
Or do they move with mesh changes? We need to add coordinate output to identify
the blowup DOF on the new mesh.

### 3.2 What makes a fault face unstable?

Possible causes (not yet verified):
1. **Element geometry**: Tet faces with high aspect ratio or near-degenerate shapes
   produce poor gradient approximation → anomalous average stress
2. **Mesh topology**: Faces where the two adjacent elements have very different sizes
   or orientations → inconsistent gradients from each side
3. **Fault edge effects**: Faces near the fault boundary (y=±50km, z=0, z=40km)
   where the fault surface terminates in the mesh
4. **BooleanFragments artifacts**: The Gmsh BooleanFragments operation creates
   sub-surfaces at nuc zone boundaries — faces near these geometric edges may
   have special topology
5. **DG basis function alignment**: On certain tet orientations, the linear basis
   may produce gradients that systematically amplify in the average stress formula

### 3.3 Does Tandem have this instability?

Tandem uses penalty=4 (dimensionless) for BR2 traction — effectively no correction.
If Tandem is stable, either:
- Tandem's mesh doesn't have the pathological faces
- Tandem's DG implementation differs in some way we haven't identified
- Tandem's time stepper handles the instability differently

We cannot test Tandem to verify.

---

## 4. The Dilemma

| Approach | Stability | VS zone | Problem |
|----------|-----------|---------|---------|
| Full BR2 correction | ✓ stable | ✗ locks up | Amplifies DG residual by μ/h, biases traction |
| No correction (IP-style) | ✗ blows up | ✓ healthy | ~1 DOF becomes unstable |
| ??? middle ground | ? | ? | Need to find |

The BR2 correction serves dual roles:
1. **Stabilization** at pathological faces (needed)
2. **Traction bias** at all faces (harmful)

We need to decouple these.

---

## 5. Proposed Approaches

### Approach A: BR2 traction + new mesh (quick test)

Revert H26 (restore original BR2 traction), keep H27 (new mesh).
The better mesh may reduce the BR2 correction magnitude at healthy faces,
reducing the VS lockup bias while maintaining stability.

**Hypothesis**: On the old mesh, the BR2 correction was ~2.4 GPa at pathological
faces and ~0.06 MPa at healthy faces. On the new mesh, the correction at healthy
faces may be smaller, reducing the VS lockup.

**Effort**: Revert one code change, resubmit.

### Approach B: Identify and exclude pathological DOFs

Add diagnostic output to find the exact coordinates of the unstable DOF(s) on the
new mesh. Then either:
- Exclude those faces from traction computation (set traction=0, let interpolation handle it)
- Use BR2 correction only at those faces, IP-style everywhere else
- Investigate the mesh/element quality at those specific faces

**Effort**: Add diagnostics, resubmit, analyze, then fix.

### Approach C: Adaptive correction strength

Use BR2 correction but with a reduced stabilization parameter. Instead of
`σ_BR2 = 4`, use `σ_BR2 = α` where α < 4. Smaller α reduces the traction bias
but may still provide enough stabilization.

Or: scale the correction by the relative jump magnitude:
```cpp
real_t jump_mag = sqrt(jump[0]*jump[0] + jump[1]*jump[1] + jump[2]*jump[2]);
real_t slip_mag = sqrt(delta_u[0]*delta_u[0] + delta_u[1]*delta_u[1] + delta_u[2]*delta_u[2]);
real_t scale = (slip_mag > 0) ? min(jump_mag / slip_mag, 1.0) : 1.0;
// Only apply correction when jump is large relative to slip
correction[i] = scale * br2_penalty * 0.5 * sum;
```

**Effort**: Implement + test.

### Approach D: Deep investigation of DG traction formula

The average stress `{{σ·n̂}}` should be bounded for a well-posed DG solution.
If it grows unboundedly, either:
- The DG solution itself is unstable (solver issue)
- The average stress formula has an error
- The feedback loop (traction → V → slip → new solve → traction) is unstable

Investigate by monitoring the displacement solution `u` at the pathological face:
print `u1`, `u2`, `grad_u1`, `grad_u2`, `[[u]]`, `slip` at each time step.
This would reveal whether `u` is growing or `grad_u` is oscillating.

**Effort**: Significant diagnostic code + analysis.

---

## 6. New Direction: Switch to IP Method

### Why IP instead of BR2

The BR2 traction issue is fundamentally about the lifting operator amplifying
the DG residual by μ/h in the traction post-processing. The correct traction
formula for BR2 requires careful derivation from the numerical flux, and the
relationship between the BR2 bilinear form penalty and the traction penalty
is subtle (see Section 3.2 above).

The **IP (Interior Penalty / SIPG) method** avoids this entirely:
- The penalty is a **scalar** with units Pa/m (properly scaled)
- The **same penalty** appears in the bilinear form AND the traction
- No lifting operator, no μ/h amplification
- Dimensionally consistent throughout
- The Tandem paper (Uphoff et al. 2023, GJI 233(1), 586-626) primarily
  develops and analyzes the IP/SIPG method, with sharp penalty bounds

MFEM already supports IP via `--dg-method IP`. The IP traction path in
`ComputeTraction` (lines 2141-2169) uses the standard IP penalty formula.

### Test plan: IP + new mesh

This sidesteps the BR2 traction problem entirely. If IP + new mesh produces
correct VS zone behavior and earthquake cycling, the BR2 traction issue
becomes a separate research question rather than a blocking bug.

### References

- Uphoff, May & Gabriel (2023). "A discontinuous Galerkin method for sequences
  of earthquakes and aseismic slip on multiple faults using unstructured
  curvilinear grids." Geophysical Journal International, 233(1), 586-626.
  https://doi.org/10.1093/gji/ggac467
- Tandem implementation: https://github.com/TEAR-ERC/tandem

---

## 7. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H25 | Previous fixes | Done |
| H26 | Traction: IP-style penalty (no BR2 lifting) | Done but causes instability on ~1 DOF |
| H27 | Mesh: uniform fault resolution (Tandem-style) | Done, committed to git |
| **H28** | **Switch DG method from BR2 to IP (`--dg-method IP`)** | **NEXT** |

---

## 8. Key Data Points

### v25 Test 1 station data (IP traction + new mesh, t=1.95yr)

| Station | V/Vp | tau_s (MPa) | Status |
|---------|------|-------------|--------|
| strk+00dp+22 (z=22km, VS) | 0.987 | 13.264 | Healthy |
| strk+00dp+10 (z=10km, VW) | 0.857 | 19.515 | Locking normally |

All 10 station DOFs healthy. The instability is at a non-station DOF.

### Vmax evolution (from .out file)

| Step | Time (yr) | Vmax (m/s) | Status |
|------|-----------|-----------|--------|
| 10 | 0.64 | 1.09e-9 | Normal (≈Vp) |
| 114 | 1.95 | 5.1e-3 | Unstable |
| 200 | 1.95 | 9.1 | Growing |
| 400 | 1.95 | 90 | Growing |
| 1048 | 1.95 | 568 | Growing (integrator stuck) |

Time is frozen at 1.95 yr — ODE integrator rejects all steps due to the
runaway DOF.
