# BP5 Debug v59: Tandem Comparison Proves MFEM-Specific Normal Traction Bug

**Date:** 2026-04-06
**Status:** Active — root cause narrowed to assembly/solve path, not friction or traction formula
**Scope:** BP5, 1000 m mesh, IP, p=1, same mesh for both codes

---

## 1. Problem Statement

MFEM's BP5 run blows up at ~25.27 yr due to localized effective normal stress
collapse at the fault-tip DOF (depth = Wf = 40 km). The elastic normal traction
grows linearly at -0.98 MPa/yr until it cancels the 25 MPa background compression.

This document records the decisive Tandem comparison that proves the linear
normal-stress erosion is an MFEM-specific bug, not expected crack-tip physics.

---

## 2. Tandem Comparison Setup

### 2.1 Config

Tandem BP5 config with fault probes at the MFEM blowup coordinates:

- `examples/tandem/3d/bp5_normal_diag.toml`
- Same mesh: `bp5_tandem_1000m.msh`
- Same solver: MUMPS direct, RK45, PETSc TS
- Added probes:
  - `tip_strk-37dp+40`: x=-37 km, depth=40 km (matches MFEM r96 fi=28 dof=0)
  - `tip_strk-36dp+40`: x=-36 km, depth=40 km
  - `tip_strk+00dp+40`: x=0 km, depth=40 km (center of fault bottom edge)
  - `int_strk-37dp+39`: x=-37 km, depth=39 km (just above tip)
  - `int_strk-36dp+39`: x=-36 km, depth=39 km
  - `int_strk+00dp+39`: x=0 km, depth=39 km
  - `ctl_strk-37dp+20`: x=-37 km, depth=20 km (shallow control)
  - `ctl_strk+00dp+20`: x=0 km, depth=20 km (shallow control)
- Output: `normal-stress` column = `snAbs = -sn + SnPre` (compression-positive,
  equivalent to MFEM's `sigma_n_eff`)

### 2.2 Tandem run

- Frontera job 7634234 (later resubmitted as 7634266)
- 8 nodes, 400 MPI ranks
- Ran successfully to 225.44 yr with no instability

### 2.3 MFEM run

- Frontera job 7633584
- 8 nodes, 400 MPI ranks
- Aborted at 25.29 yr due to NaN psi at the fault-tip DOF

---

## 3. Decisive Result: Normal Stress at Fault Tip (x=-37 km, depth=40 km)

### 3.1 Side-by-side comparison

| Time (yr) | Tandem sigma_n_eff (MPa) | MFEM sigma_n_eff (MPa) | Difference (MPa) |
|-----------|--------------------------|------------------------|-------------------|
| 0         | 25.000                   | 25.000                 | 0.000             |
| 10        | 24.641                   | 15.040                 | 9.601             |
| 20        | 24.669                   | 5.230                  | 19.439            |
| 24        | 24.681                   | 1.250                  | 23.431            |
| 25.27     | 24.685                   | 0.019                  | 24.666            |

### 3.2 Growth rates

- **Tandem**: slope = +0.001 MPa/yr (essentially flat, < 1.3% erosion at 25 yr)
- **MFEM**: slope = -0.98 MPa/yr (linear collapse, cancels 25 MPa at 25.27 yr)

### 3.3 Tandem at all probes (t ≈ 25 yr)

| Probe                 | sigma_n_eff (MPa) | Erosion (%) |
|-----------------------|-------------------|-------------|
| TIP  x=-37 z=-40      | 24.684            | 1.26        |
| TIP  x=-36 z=-40      | 24.603            | 1.59        |
| TIP  x=  0 z=-40      | 25.030            | -0.12       |
| INT  x=-37 z=-39      | 25.005            | -0.02       |
| INT  x=-36 z=-39      | 25.039            | -0.15       |
| INT  x=  0 z=-39      | 25.164            | -0.66       |
| CTL  x=-37 z=-20      | 25.053            | -0.21       |
| CTL  x=  0 z=-20      | 25.003            | -0.01       |

Tandem shows < 2% erosion everywhere, even at the fault-tip corners. The
interior and control probes are essentially at the baseline 25 MPa.

---

## 4. What This Proves

### 4.1 The -0.98 MPa/yr erosion is NOT expected physics

Both codes use:
- Same mesh (`bp5_tandem_1000m.msh`)
- Same DG IP formulation
- Same penalty formula: `(p0+p1)/4`, trace inequality `p*(p+D-1)/D`
- Same Dirichlet loading: `boundary(x,y,z,t)` from `bp5.lua`
- Same friction law: Dieterich-Ruina aging, same parameters
- Same solver: MUMPS direct, RK45

Tandem does not see the erosion. Therefore the erosion is an MFEM implementation
artifact, not a physical consequence of the DG formulation on this mesh.

### 4.2 The bug is NOT in

The following have been verified through line-by-line code comparison and/or
unit tests, and are confirmed consistent between MFEM and Tandem:

1. **Traction formula** — stress average `{sigma.n_hat}` and penalty correction
   `-eta*(u1-u2-slip)` match exactly (v58, Section 3.5)
2. **Penalty scaling** — trace inequality constant, area/volume ratio, material
   handling all match (v58, Sections 3.5.4–3.5.7)
3. **Sign conventions** — normal convention, sign_flipped handling, dir_sign for
   Dirichlet loading all match Tandem's facet-functional pattern
4. **L2 projection** — same per-face mass matrix, same quadrature rule (2p+1),
   same nl_q weighting, same basis-vector rotation
5. **Normal extraction** — both use `-(T.n_hat)` for compression-positive
6. **Face classification** — no face gets both fault and Dirichlet loading; strict
   partition with validated tag recovery
7. **Stiffness matrix K** — assembled identically on all interior faces regardless
   of face BC type
8. **RHS partitioning** — fault faces get slip, Dirichlet faces get prescribed
   displacement; no overlap
9. **Friction solver** — same sigma_n <= 0 guard (V = tau/eta), same state
   evolution law, same dpsi/dt formula
10. **Owned DOF ordering** — FaultGeometry uses owned-restricted coordinates
    correctly

### 4.3 The bug IS in

The displacement field `u` that MFEM computes from `K * u = b_slip + b_dirichlet`
must be different from Tandem's, because the traction formula `T(u)` is verified
identical but the output `T` is drastically different.

Specifically, the normal component of traction at the fault-tip DOF diverges
linearly with time in MFEM but stays flat in Tandem. Since `T = f(u)` and `f` is
the same, `u` must be different.

The remaining suspects are in the **assembly of `b`** (the RHS vector):

- `AssembleSlipContributionIP()` — how fault slip enters the RHS
- `AssembleDirichletLoading()` — how the Dirichlet prescribed displacement
  enters the RHS on interior skeleton faces
- The interaction between these two at the fault/Dirichlet junction

Or in the **stiffness matrix K assembly** itself, despite the formula matching:

- Quadrature rule differences in the actual `BilinearForm` assembly path
- Element-side Jacobian handling in MFEM's `FaceElementTransformations`
- Any MFEM-internal transform caching that could give wrong values on specific
  face geometries

---

## 5. Evidence Supporting Assembly/Solve as Root Cause

### 5.1 v58 first-solve jump_y mismatch

The v58 investigation (Section 4.4) found that at the very first accepted step:

- MFEM `jump_y` ~ 1e-11 to 1e-10
- Tandem `jump_y` ~ 1e-08

A 100–1000x difference in the displacement jump at the first solve means the
displacement field is already different. This is upstream of traction recovery.

### 5.2 Linear growth rate

The -0.98 MPa/yr rate is consistent with a systematic, time-proportional error
in the RHS. Since the Dirichlet loading grows as `Vp * t` and the slip also
grows approximately linearly, a wrong coefficient or sign in one RHS term would
produce exactly this linear growth pattern.

### 5.3 Localization to fault tip

The erosion is worst at the fault-tip DOFs (depth = 40 km) and negligible at
interior DOFs (depth = 39 km) and controls (depth = 20 km). This pattern
matches a boundary-effect or junction-effect in the RHS assembly, not a global
formula error.

---

## 6. Potential Issue 1: IP Dirichlet RHS Still Uses `Face->Transform`

### 6.1 Why this is suspicious

In the MFEM IP Dirichlet RHS path, the BP5 prescribed displacement is chosen by
evaluating physical coordinates at each quadrature point and then applying the
piecewise `boundary(x,y,z,t)` logic from `bp5.lua`:

- `Vh = Vp * t / 2` if `y > 1 km`
- `Vh = -Vp * t / 2` if `y < -1 km`
- `Vh = Vp * t` if `|y| <= 1 km`

The current code computes `y` using `Face->Transform(...)` in all three IP
Dirichlet RHS sites:

- true boundary Dirichlet faces
- interior `y=0` Dirichlet skeleton faces
- shared `y=0` Dirichlet skeleton faces

This is suspicious because the earlier geometry-getter bug already showed that
`Face->Transform` can disagree with the fully initialized element transform on
the BP5/tet path. The getters were only made consistent after switching to:

- `FTr->SetAllIntPoints(&ipq)`
- `FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords)`

So the Dirichlet RHS is still using the same transform path that was already
proven unreliable elsewhere in this code.

### 6.2 Why this could cause the observed linear erosion

This is not just a diagnostic issue. The returned physical coordinate directly
controls which BP5 loading branch is applied at that quadrature point.

If `Face->Transform(...)` gives a slightly wrong `y` on a subset of faces, MFEM
can apply the wrong prescribed displacement:

- `Vp*t` instead of `Vp*t/2`
- the wrong sign on one side of the `y = +/- 1 km` threshold
- or the wrong continuation value on the interior/shared `y=0` faces adjacent
  to the fault

Because the imposed loading scales linearly with time, any systematic branch
selection error in the Dirichlet RHS would also produce a linearly growing
error in the displacement field and therefore a linearly growing traction
error. That is consistent with the observed `-0.98 MPa/yr` normal-stress
erosion.

### 6.3 Why this fits the localization

The most suspicious locations are the interior/shared `y=0` Dirichlet faces near
the fault/Dirichlet junction at the lower fault tip. That is exactly where the
MFEM erosion is strongest:

- localized to the bottom patch
- strongest on the fault-tip DOFs
- negligible at the interior and shallow control probes

This is the spatial pattern expected from a local RHS loading inconsistency,
not from a global traction formula error.

### 6.4 Recommended check/fix

As a narrow IP-only test, replace the coordinate evaluation in the Dirichlet
RHS paths with the same robust element-transform pattern already used in the
fixed geometry getters:

```cpp
FTr->SetAllIntPoints(&ipq);
Vector phys(dim);
FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
```

Use `phys(1)` from that element-transform path when applying the BP5
`boundary(x,y,z,t)` piecewise logic.

If the linear tip erosion changes materially after this replacement, the bug is
in the IP Dirichlet coordinate-evaluation path rather than in the traction
formula itself.

---

## 7. Diagnostic Bug Fixed in This Investigation

### 7.1 Normal decomposition projection bug

In `dg_elasticity_ip_combined_integrator.hpp`, the `ProjectTractionToFaultDOFs`
function had a bug when called with `ncomp_local=1` and per-QP `qp_data` provided.

**Bug**: `basis_vec[0]` was set to `qd.tangent1` instead of `qd.normal` when
`ncomp_local < 3`. This caused the normal decomposition outputs (`normal_stress`,
`normal_correction`) to project onto the wrong vector (tangent1 instead of normal).

**Impact**: The diagnostic `normal_stress` and `normal_corr` columns in
`trace_timeseries_r*.csv` were wrong. The decomposition did not close:
`normal_stress + normal_corr != normal_traction`. This prevented correct
diagnosis of whether the stress or correction term drives the sigma_n_eff collapse.

**Fix**: Changed `basis_vec[0]` selection so that for `ncomp_local=1`, it uses
`n0` (which is `qd.normal` in the qpd path, or `tangents[0]` in the non-qpd path).

**Production impact**: None. The production `normal_traction_` used for
`sigma_n_eff` in friction comes from the `ncomp=3` unified projection, which was
always correct. This bug only affected diagnostic CSV outputs.

**Verification**: After the fix, the TestNormalTractionLeakageIP test shows
decomposition closure to machine precision (`jumpres ~ 1e-14`).

---

## 8. Immediate Next Debug Steps

### 8.1 Compare the displacement field directly

The most direct way to isolate the RHS bug is to compare MFEM's and Tandem's
displacement `u` at the first accepted step:

- Dump MFEM's `u` at element DOFs adjacent to the fault tip
- Dump Tandem's `u` at the equivalent element DOFs
- Compare `u1 - u2` (the DG jump) on the last fault face

If `u` differs at the first step, the bug is in the RHS assembly (since K is
the same and both use direct solvers).

### 8.2 Compare the RHS vector b

Dump the assembled RHS vector `b` at the first step in both codes:

- `b_slip` contribution from `AssembleSlipContributionIP()`
- `b_dirichlet` contribution from `AssembleDirichletLoading()`
- Total `b = b_slip + b_dirichlet`

Compare element-by-element on the elements adjacent to the fault tip.

### 8.3 Focused unit test

Build a minimal test case:

- Small mesh with a fault rectangle and Dirichlet continuation on y=0
- Apply Vp*t loading for one step
- Compare the resulting displacement and traction against Tandem's output
- Focus on the last fault face at the fault/Dirichlet junction

---

## 9. How to Interpret the First-Step Dumps

The first-step dumps are not just a "same or different" check. Their value is
that they identify the earliest layer where MFEM diverges from Tandem.

### 9.1 If `b_slip` differs first

If the fault-slip RHS contribution already differs on the bottom-tip fault face,
then the bug is in the **fault slip RHS path**, not in the Dirichlet loading.

This implicates:

- `AssembleSlipContributionIP()`
- `AssembleSlipContributionIPShared()`
- `DGElasticityIPCombinedIntegrator::AssembleSlipFaceRHS()`

Likely sub-causes:

- slip orientation or sign handling
- shared-face local/neighbor handling
- per-QP fault basis / sign mismatch

### 9.2 If `b_dirichlet` differs first

If `b_slip` matches but the Dirichlet contribution differs, then the bug is in
the **Dirichlet RHS path**.

This implicates:

- `AssembleDirichletLoading()`
- `DGElasticityIPCombinedIntegrator::AssembleBoundaryFaceRHS()`
- the skeleton Dirichlet reuse of `AssembleSlipFaceRHS()`

Likely sub-causes:

- wrong physical coordinates used for `boundary(x,y,z,t)`
- wrong interior/shared Dirichlet face orientation
- wrong shared-face handling on `y=0`

### 9.3 If `b_slip` and `b_dirichlet` match, but total `b` differs

Then the local formulas are likely correct and the problem is in **assembly
combination or scatter**, not in the face kernels themselves.

This points to:

- global scatter into `rhs`
- duplicated or omitted face contributions
- shared-face ownership logic
- overlap between interior/shared/boundary processing

### 9.4 If total `b` matches, but first-step `u` differs

Then the issue is in **K or the solve path**, not in the RHS.

This implicates:

- `AssembleStiffness()`
- shared-face matrix assembly
- neighbor-block handling
- transform use inside the face-matrix assembly path

With direct solvers, a first-step mismatch here is strong evidence that the
stiffness assembly or matrix export path differs.

### 9.5 If `u` matches, but traction differs

Then the bug is in **traction recovery**, not in assembly/solve.

This would implicate:

- `ComputeTractionImpl()`
- quadrature-point traction computation
- projection to fault DOFs

Current evidence argues this is less likely, but the interpretation rule is
useful to keep explicit.

### 9.6 MPI-specific interpretation

The same first-step dump can also tell whether the bug is local or shared-face
parallel-specific:

- If the mismatch appears only on shared faces and not on purely interior faces,
  suspect a **parallel shared-face bug**.
- If the mismatch appears on both interior and shared faces of the same class,
  suspect a **common local kernel issue**.

The shared-face-only case would point most strongly to:

- `AssembleSlipContributionIPShared()`
- the shared Dirichlet IP path in `AssembleDirichletLoading()`
- local `Elem1` vs neighbor `Elem2` handling

### 9.7 Summary decision table

- `b_slip` differs first → fault slip RHS path
- `b_dirichlet` differs first → Dirichlet RHS path
- only total `b` differs → scatter / duplication / omission
- `b` matches, `u` differs → stiffness matrix / solve path
- `u` matches, traction differs → traction recovery path

This is why the first-step dump is the highest-signal next diagnostic. It tells
which layer first diverges, instead of inferring root cause indirectly from the
late-time fault-tip blowup.

### 9.8 Current MFEM first-step CSV outputs

The current IP-only first-step debug mode writes four rank-local CSVs into the
BP5 `--output-dir` on the selected debug rank.

#### `first_step_face_rhs_r<rank>.csv`

Purpose:

- face-local RHS contributions before global scatter
- both fault-slip and Dirichlet face assembly paths

Columns:

- `time`
- `phase`
  - `slip` or `dirichlet`
- `face_kind`
  - `fault_interior`
  - `fault_shared`
  - `dirichlet_boundary`
  - `dirichlet_interior`
  - `dirichlet_shared`
- `fault_idx`
  - local fault-face index for fault faces
  - `-1` for Dirichlet-only faces
- `mesh_face`
- `elem1`, `elem2`
- `record`
  - `phys_y`
  - `input_qp`
  - `elvec1`
  - `elvec2`
- `index0`, `index1`
- `value`

Interpretation:

- `phys_y` checks the physical coordinate used in the BP5 `boundary(x,y,z,t)`
  branch logic
- `input_qp` is the imposed vector at each face quadrature point
  - `u_D_3d` for Dirichlet faces
  - embedded slip for fault faces
- `elvec1` / `elvec2` are the element-local face RHS vectors before scatter

#### `first_step_elem_data_r<rank>.csv`

Purpose:

- assembled element-restricted vectors after face scatter
- first-step solved displacement on the same elements

Columns:

- `time`
- `quantity`
  - `rhs_slip`
  - `rhs_dirichlet`
  - `rhs_total`
  - `u`
- `elem`
- `component`
- `local_dof`
- `vdof`
- `value`

Interpretation:

- use this to compare the assembled algebraic vectors on the elements adjacent
  to the target fault patch
- if face-local terms match but `rhs_total` differs here, the problem is in
  scatter/ownership rather than the local kernels

#### `first_step_face_jump_r<rank>.csv`

Purpose:

- post-solve fault jump check at face quadrature points

Columns:

- `time`
- `face_kind`
  - `fault_interior`
  - `fault_shared`
- `fault_idx`
- `mesh_face`
- `elem1`, `elem2`
- `record`
  - `u1_q`
  - `u2_q`
  - `slip_q`
  - `jump_minus_slip`
- `q`
- `component`
- `value`

Interpretation:

- `jump_minus_slip = (u1 - u2) - slip`
- this is the most direct first-step check for whether the DG interface state is
  already inconsistent on the bad fault patch

#### `first_step_face_trac_r<rank>.csv`

Purpose:

- post-solve IP traction evaluation on the same target fault faces

Columns:

- `time`
- `face_kind`
  - `fault_interior`
  - `fault_shared`
- `fault_idx`
- `mesh_face`
- `elem1`, `elem2`
- `record`
  - `nl_q`
  - `traction_q`
  - `traction_dip`
  - `traction_strike`
  - `normal_traction`
- `index0`, `index1`
- `value`

Interpretation:

- `traction_q` is the quadrature-point traction vector from the production IP
  traction path
- `traction_dip`, `traction_strike`, and `normal_traction` are the projected
  fault-DOF outputs derived from the same solve result
- this file should only be interpreted after `b`, `u`, and `jump` are compared

---

## 10. Constant-σn Test: MFEM Stress Buildup Is 5× Too Slow

### 10.1 Test setup

An MFEM BP5 run was performed with **forced constant normal stress** (elastic
σn feedback disabled, `sigma_n_eff = sigma_n_bp5_ = 25 MPa` at all times).
This eliminates the normal traction erosion and allows the run to proceed
indefinitely.

- Frontera job: 7634218
- Config: same mesh, IP p=1, 1000 m, MUMPS, RK45
- Only change: elastic normal stress disabled (`elastic_sigma_n_ = false`)
- Ran to 1800 yr without blowup

### 8.2 Earthquake recurrence

| Code               | First earthquake (yr) | Recurrence |
|--------------------|-----------------------|------------|
| Tandem (elastic σn)| ~225                  | Consistent with SCEC BP5 |
| MFEM (const σn)    | ~1191                 | 5.3× too slow |

MFEM's first earthquake is delayed by a factor of 5.3 compared to Tandem.

### 8.3 Shear stress buildup comparison

At center station (strk+00dp+00, surface, VW region):

| Time (yr) | MFEM |τ| (MPa) | Tandem |τ| (MPa) | Deficit (MPa) |
|-----------|----------------|------------------|---------------|
| 1         | 12.41          | 12.35            | +0.06         |
| 10        | 11.90          | 11.93            | -0.04         |
| 50        | 11.56          | 12.01            | -0.45         |
| 100       | 11.50          | 12.55            | -1.05         |
| 150       | 11.52          | 13.11            | -1.59         |
| 200       | 11.58          | 13.22            | -1.64         |

MFEM's shear stress is nearly flat (11.5 MPa), while Tandem's builds up from
12.4 to 13.2 MPa over 200 yr. The stress accumulation rate is drastically
different.

At depth=10 km station (strk+00dp+10, VW region), the deficit is even larger:

| Time (yr) | MFEM |τ| (MPa) | Tandem |τ| (MPa) | Deficit (MPa) |
|-----------|----------------|------------------|---------------|
| 10        | 9.75           | 9.89             | -0.13         |
| 100       | 10.38          | 12.34            | -1.97         |
| 200       | 10.95          | 16.46            | -5.51         |

### 8.4 Slip accumulation comparison

At center station (strk+00dp+00):

| Time (yr) | MFEM slip (m) | Tandem slip (m) | Deficit (m) |
|-----------|---------------|-----------------|-------------|
| 1         | 5.327         | 5.328           | -0.001      |
| 10        | 5.450         | 5.457           | -0.008      |
| 50        | 5.542         | 5.582           | -0.040      |
| 100       | 5.585         | 5.698           | -0.113      |
| 200       | 5.639         | 6.362           | -0.723      |

MFEM accumulates 0.72 m less slip over 200 yr than Tandem. The deficit grows
with time, consistent with a loading rate error.

At depth=10 km station (strk+00dp+10), MFEM shows essentially **zero slip
accumulation** after t ≈ 1 yr (stuck at 6.5705 m), while Tandem continues.

### 8.5 Interpretation

The constant-σn test proves two things:

1. **The normal traction erosion is a symptom, not the root cause.** Even with
   normal stress forced constant, MFEM's stress buildup is 5× too slow. The
   earthquake cycle is fundamentally delayed.

2. **The bug is in the elastic loading path.** The Dirichlet loading `Vp·t` is
   not producing the correct stress accumulation rate on the fault. The same
   error that makes shear stress build up too slowly also makes normal traction
   grow spuriously — both are consequences of a wrong displacement field.

This narrows the root cause further: the assembled RHS vector `b` or the
stiffness matrix `K` must produce a displacement field that evolves too slowly
under the applied far-field loading. Since the formulas for K and b have been
verified line-by-line, the error is likely a coefficient, sign, or scaling
issue in one specific assembly path — most likely in `AssembleDirichletLoading()`
where the `Vp·t` boundary condition enters the DG skeleton RHS.

### 8.6 Updated root cause ranking

Given both the normal traction comparison (Section 3) and the constant-σn test:

1. **Highest**: `AssembleDirichletLoading()` — the Dirichlet skeleton RHS
   produces wrong loading on interior y=0 faces. A coefficient or scaling error
   here would directly explain both the slow stress buildup and the spurious
   normal traction growth.

2. **High**: `AssembleSlipContributionIP()` — the fault slip RHS has a wrong
   coefficient. Less likely because the slip itself matches Tandem at early
   times (deficit grows slowly).

3. **Medium**: Stiffness matrix K assembly — a quadrature or Jacobian error in
   the `BilinearForm` assembly that makes K too stiff. Both codes use the same
   formulas, but MFEM's internal implementation may differ.

4. **Low**: Friction/state evolution — the const-σn test uses the same friction
   law and still shows the delay, so friction is not the primary cause.

---

## 11. Corrected Normal Decomposition: Penalty Correction Drives 84% of the Collapse

### 11.1 Data source

After fixing the `ProjectTractionToFaultDOFs` basis-vector bug (Section 6), a
new BP5 run with the corrected tracer produced the first trustworthy normal
decomposition at the fault-tip DOF.

- Trace data: `Archive_2026-04-06T16_44_22/trace_*_r96.csv`
- Worst DOF: fi=28, dof=0, depth=40 km, x2=-37 km (same DOF as all previous runs)
- Decomposition closure: `|Tn_total - (Tn_stress + Tn_corr)| < 1e-10 Pa` (machine precision)

### 11.2 Normal traction decomposition at the worst DOF

| Time (yr) | Tn_total (Pa) | Tn_stress (Pa) | Tn_corr (Pa) | corr/stress | % from corr |
|-----------|--------------|----------------|--------------|-------------|-------------|
| 0.1       | -1.75e5      | -5.26e4        | -1.22e5      | 2.3         | 69.8%       |
| 1         | -1.11e6      | -2.15e5        | -8.97e5      | 4.2         | 80.6%       |
| 5         | -5.16e6      | -8.78e5        | -4.28e6      | 4.9         | 83.0%       |
| 10        | -1.00e7      | -1.66e6        | -8.38e6      | 5.1         | 83.5%       |
| 20        | -1.99e7      | -3.16e6        | -1.67e7      | 5.3         | 84.1%       |
| 25.27     | -2.50e7      | -3.93e6        | -2.10e7      | 5.4         | 84.3%       |

### 11.3 Growth rates

Linear fits for t > 1 yr:

| Component       | Slope (Pa/yr) | Fraction of total |
|-----------------|--------------|-------------------|
| Tn_corr         | -828,440     | 84.4%             |
| Tn_stress       | -152,780     | 15.6%             |
| Tn_total        | -981,220     | 100%              |

Both stress and correction grow linearly. The correction term grows 5.4× faster
than the stress term and is responsible for 84% of the sigma_n_eff erosion.

### 11.4 Implied normal opening

The correction term `Tn_corr = -penalty * (u1_y - u2_y)` directly measures the
normal displacement jump across the fault face. At t=25 yr with
penalty ≈ 2.183e9 Pa:

| Code   | Normal corr (Pa) | Implied normal jump | Notes |
|--------|------------------|--------------------:|-------|
| MFEM   | -2.10e7           | 9.5 mm             | From trace data |
| Tandem | ~-3.2e5           | 0.15 mm            | Estimated from 1.3% sigma_n erosion |

MFEM's displacement field develops a **65× larger spurious normal jump** at the
fault-tip face than Tandem's. Both codes use the same penalty, the same traction
formula, and the same mesh. The jump difference must come from the displacement
solution `u`.

### 11.5 Interpretation

The blowup is NOT caused by the elastic stress average `{σ·n̂}` — that term
contributes only 3.9 MPa at 25 yr, which is small relative to the 25 MPa
background compression.

The blowup IS caused by the DG penalty correction `-η * (u1_y - u2_y)`, which
reaches 21 MPa by 25 yr. This term measures how much the displacement solution
violates the no-opening condition on the fault face. MFEM's solution violates
it 65× more than Tandem's.

Since the penalty formula, traction formula, sign conventions, and L2 projection
all match Tandem term-by-term (Sections 4.2, 11.2), the 65× difference in the
normal jump must originate upstream: in the assembled RHS `b` or the stiffness
matrix `K` that produce the displacement field `u = K^{-1} b`.

### 11.6 Updated root cause ranking

The decomposition narrows the suspect further:

1. **Highest**: The RHS vector `b` produces a displacement `u` with too much
   normal discontinuity at fault-tip faces. The most likely sub-cause is a subtle
   issue in how `b_slip` or `b_dirichlet` load the elements adjacent to the
   fault tip — not a formula error (verified), but a numerical value difference
   from MFEM's internal mesh/transform machinery.

2. **High**: The stiffness matrix `K` is assembled with slightly different
   effective stiffness on fault-tip elements, making the penalty less effective
   at suppressing the normal jump.

3. **Medium**: An MFEM-internal issue in `FaceElementTransformations` or
   `CalcOrtho` that produces subtly different Jacobians/normals on specific
   face geometries, affecting both K and b consistently but differently from
   Tandem.

4. **Low**: The traction formula or projection itself — ruled out by
   line-by-line code comparison and decomposition closure.

---

## 13. First-Step Cross-Code Comparison: RHS Matches, Displacement Differs

### 13.1 Setup

Both codes were run with first-step dump diagnostics on the same mesh
(`bp5_tandem_1000m.msh`) to compare the assembled RHS element vectors,
displacement jumps, and traction at fault-tip faces.

- **MFEM**: Frontera job 7635265, rank 96, fault faces fi=4/28/29/33
  - Dump time: t = 0.002 s (first RK45 stage with t > 0, dt_init=0.01)
  - Units: meters, Pa

- **Tandem**: Frontera job 7635260, rank 40, matched faces fct=31/38/66/88/146/160/182
  - Dump time: t = 0.02 s (first accepted step, dt_init=0.01)
  - Units: km (coordinates), meters (displacement)

Face matching by physical centroid coordinates:

| MFEM face | Location | Tandem face |
|-----------|----------|-------------|
| fi=29     | x=-36.5 km, depth=39.8 km | fct=31 |
| fi=33     | x=-37.0 km, depth=39.5 km | fct=38 |
| fi=4      | x=-36.0 km, depth=39.1 km | fct=146 |

### 13.2 RHS element vectors match exactly

The per-face RHS element vectors (`elvec1`, `elvec2`) match between the two
codes to **7 significant figures**, with a constant ratio of 1.00e+11 that is
the expected unit conversion factor between MFEM (meters) and Tandem (km).

| Face pair | MFEM elvec1 | Tandem elvec2 | Ratio |
|-----------|-------------|---------------|-------|
| fi=29 ↔ fct=31 | 3.224545e+02 | 3.224545e-09 | 1.00e+11 |
| fi=33 ↔ fct=38 | 3.940957e+02 | 3.940957e-09 | 1.00e+11 |
| fi=4 ↔ fct=146 | 3.633561e+02 | 3.633561e-09 | 1.00e+11 |

The elem1/elem2 swap between some pairs is from different mesh element
numbering. The digit-level agreement proves:

- `AssembleSlipFaceRHS` produces identical values to Tandem's `rhs_skeleton`
- `AssembleDirichletLoading` produces identical values (same code path)
- The penalty, quadrature, shape functions, and Jacobians all agree
- The sign/orientation conventions are correct

**The RHS assembly is not the bug.**

### 13.3 Displacement field differs at first step

Despite identical RHS, the displacement solution `u` differs:

| Quantity | MFEM (t=0.002s) | Tandem (t=0.02s) | Expected ratio | Actual ratio |
|----------|-----------------|-------------------|----------------|-------------|
| u_y (normal disp) | ~4e-8 m | ~2e-7 m | 10x (time) | 5x |
| jump_y (normal jump) | ~1e-11 m | ~1e-8 m | 10x (time) | ~1000x |

After accounting for the 10x time difference:
- MFEM u_y is **~2× too large** (extrapolated)
- MFEM jump_y is **~100× too small**

The normal displacement is too large but the normal jump (u1−u2) is too small.
This means MFEM's solution is more **continuous** across the fault face — both
elements agree more on u_y — while Tandem's solution has a larger physical
discontinuity, as expected from the DG formulation.

### 13.4 RHS Dirichlet on fault-tip elements is zero

The MFEM elem_data dump shows:

```
rhs_dirichlet: 96 values, max|v|=0.000000e+00
```

This is expected: the target elements are fault-face elements, and Dirichlet
loading is applied to Dirichlet faces (not fault faces). The Dirichlet loading
affects fault-tip elements only through the global K matrix coupling.

### 13.5 Conclusion: bug is in the stiffness matrix K

Since `u = K⁻¹ · b` and `b` is proven identical, the difference in `u` must
come from the stiffness matrix `K` or the linear solver.

The stiffness matrix K is assembled by MFEM's `ParBilinearForm` using
`DGElasticityIPCombinedIntegrator::AssembleFaceMatrix`. The formula in
`AssembleFaceMatrix` has been verified line-by-line against Tandem's
`assembleSurface` kernel. But the global assembly process — how MFEM's
`BilinearForm::Assemble()` and `ParBilinearForm::AssembleSharedFaces()`
collect, scatter, and store face matrix contributions — has NOT been directly
verified against Tandem.

The most likely sub-causes:

1. **A face processed by the wrong integrator path**: if an interior face
   is accidentally processed as a boundary face (or vice versa) in K but not
   in b, the K/b consistency would be broken.

2. **A shared face with wrong K contribution**: if the shared-face K assembly
   scatters entries to the wrong DOFs or with the wrong sign.

3. **The volume integrator** (`ElasticityIntegrator`): not directly verified
   against Tandem's volume assembly. A difference in the volume K terms would
   affect the displacement but not the face RHS.

4. **The MUMPS interface**: if `HypreParMatrix` drops or misindexes entries
   when converting to MUMPS format. Unlikely but not excluded.

### 13.6 Next step

Compare K matrix entries directly:

- Dump specific rows/columns of K from both codes for elements adjacent to
  the fault tip
- Or: apply K to known unit vectors and compare the output on ALL ranks
  simultaneously (MPI-safe)
- Focus on the volume integrator contribution, which was never directly
  compared

---

## 14. First-Step Dump Data References

### MFEM first-step dump

- Frontera job: 7635265 (development, 2hr, 8 nodes, 400 ranks)
- Output: `results_v59_first_step_job7635265/first_step_*_r96.csv`
- Dump time: t = 0.002 s (first RK45 stage with t > 0)
- Target rank: 96, fault faces: fi=4/28/29/33

### Tandem first-step dump

- Frontera job: 7635260 (development, 2hr, 8 nodes, 400 ranks)
- Output: `results_tandem_1000m_p1_job7635260/first_step_*_r40.csv`
- Dump time: t = 0.02 s (first accepted step after initial)
- Matched faces on rank 40: fct=31/38/66/88/146/160/182
- Also rank 17: 2 Dirichlet skeleton faces

---

## 15. Per-Element K Matrix Comparison: Integrators Are Identical

### 15.1 Setup

Dumped per-element volume K (`ElasticityIntegrator::AssembleElementMatrix`) and
per-face skeleton K (`DGElasticityIPCombinedIntegrator::AssembleFaceMatrix`)
from both codes on the production 1000 m mesh.

- MFEM: Frontera job 7636174, rank 96, target elements adjacent to fi=4/28/29/33
- Tandem: Frontera job 7636202, rank 40, coordinate-matched faces + all local volumes

### 15.2 Face K comparison (sorted absolute values, DOF-ordering independent)

Three matched fault-tip faces compared, all 4 blocks (A00, A01, A10, A11):

| Face | MFEM face | Tandem fct | Entries compared | Ratio MFEM/Tandem | Max deviation |
|------|-----------|------------|------------------|-------------------|---------------|
| x=-36.5km, d=39.8km | 144 | 31 | 12 per block × 4 blocks | 1.0000000000e+12 | 1.2e-15 |
| x=-37.0km, d=39.5km | 159 | 38 | 12 per block × 4 blocks | 1.0000000000e+12 | 1.1e-15 |
| x=-36.0km, d=39.1km | 17  | 146 | 12 per block × 4 blocks | 1.0000000000e+12 | 5.6e-16 |

The 1e12 ratio is the expected unit scaling (meters vs km coordinate system).
Every sorted absolute value matches to **10 significant figures**. Max deviation
from the mean ratio is at **machine epsilon** (~1e-15).

Full-matrix norms confirm:

| Face | Frobenius ratio/1e12 | L1 ratio/1e12 | Max ratio/1e12 |
|------|---------------------|---------------|----------------|
| 144↔31  | 1.0000000000 | 1.0000000000 | 1.0000000000 |
| 159↔38  | 1.0000000000 | 1.0000000000 | 1.0000000000 |
| 17↔146  | 1.0000000000 | 1.0000000000 | 1.0000000000 |

### 15.3 Volume K comparison (sorted absolute values)

Four matched element pairs compared (parents of the matched faces):

| MFEM elem | Tandem elem | Entries | Ratio | Max deviation |
|-----------|-------------|---------|-------|---------------|
| 48 | 129 | 128 | 1.0000000000e+12 | 7.1e-15 |
| 85 | 8   | 128 | 1.0000000000e+12 | 4.0e-15 |
| 56 | 10  | 144 | 1.0000000000e+12 | 7.6e-15 |
| 105 | 74 | 144 | 1.0000000000e+12 | 9.7e-15 |

Volume K matches to the same precision as face K.

### 15.4 CG vs MUMPS solver comparison

Both solvers produce the same displacement to machine precision:

- max|u_mumps - u_cg| = 2.26e-14
- relative diff = 5.17e-07

RHS is bit-identical between MUMPS and CG runs (zero diff).

The solver is not the issue.

### 15.5 Summary of what matches and what doesn't

| Quantity | MFEM vs Tandem | Status |
|----------|---------------|--------|
| Face K (per-face local matrix) | ratio = 1e12 exact | ✓ identical |
| Volume K (per-element local matrix) | ratio = 1e12 exact | ✓ identical |
| Face RHS elvec (per-face) | ratio = 1e11 exact | ✓ identical |
| Input f_q | ratio = 0.1 exact | ✓ identical (10x time) |
| CG vs MUMPS displacement | diff = 2e-14 | ✓ identical |
| **Displacement u** | **5-8x different** | **✗ different** |
| **Normal jump** | **100-1000x different** | **✗ different** |

### 15.6 Conclusion

The integrators produce bit-identical local matrices. The solver produces
the correct inverse. Yet the global displacement differs.

The bug must be in the **global assembly** step — how the local matrices
are scattered into the global sparse matrix via
`BilinearForm::Assemble() → mat->AddSubMatrix(vdofs, vdofs, elemmat)`.

Possible failure modes:

1. **Wrong vdofs**: if `GetElementVDofs` returns different DOF indices for
   some elements compared to Tandem's element-to-global DOF mapping, the
   local matrix entries would be placed at wrong positions in the global
   matrix.

2. **Face processed wrong number of times**: if an interior face is
   processed twice (by both the interior loop and the boundary loop) or
   skipped entirely, the global matrix would have extra or missing
   contributions even though each individual contribution is correct.

3. **Shared face assembly mismatch**: if `AssembleSharedFaces` scatters
   entries to wrong DOF positions for some shared faces, the global K
   would differ from Tandem's even though the per-face K is correct.

4. **ParallelAssemble DOF mapping**: if the local-to-global DOF mapping
   in `ParallelAssemble` (via `face_nbr_glob_ldof`) has an error for
   specific faces, the HypreParMatrix would have misplaced entries.

### 15.7 K dump data references

- MFEM: `results_v59_fs_revert_job7636174/first_step_K_volume_r96.csv`,
  `first_step_K_skeleton_r96.csv`
- Tandem: `results_tandem_1000m_p1_job7636202/first_step_K_volume_r40.csv`,
  `first_step_K_skeleton_r40.csv`

---

## 16. K Assembly Face Count Comparison: Boundary Face Mismatch

### 16.1 Global assembly counts (corrected)

| Quantity | MFEM | Tandem | Match? |
|----------|------|--------|--------|
| Volume elements | 62874 | 62874 | ✓ |
| Interior faces (both-inside) | 106618 | 106031 | ✗ (+587 in MFEM) |
| Shared faces (one-inside) | 36006 | 37180 | ✗ (-1174 in MFEM) |
| Interior+Shared total | 142624 | 143211 | ✗ (-587 in MFEM) |
| Bdr faces (Dirichlet K) | **1628** | — | — |
| Bdr elements (all attrs) | 12954 | — | — |
| Tandem boundary faces (one-sided) | — | **2254** | — |

### 16.2 MFEM per-attr boundary element breakdown

| Attr | Boundary elements | K integrator? | Tandem BC |
|------|-------------------|---------------|-----------|
| 1 | 2014 | **NO** | Natural (skipped) |
| 3 | 9312 | **NO** | Fault (K added!) |
| 5 | 1628 | **YES** (Dirichlet) | Dirichlet (K added) |
| Total | 12954 | 1628 | — |

Note: 12954 is ALL boundary elements. Many are interior faces in MFEM
(both parents local, `GetBdrFaceTransformations` returns null). The count
of TRUE one-sided boundary faces is ~2254 (matching Tandem).

### 16.3 The discrepancy

MFEM's `AddBdrFaceIntegrator` with `dirichlet_bdr_marker_` only adds boundary
K terms for **attr 5** (1628 faces). Tandem's `assemble_boundary` adds K terms
for ALL non-Natural boundary faces — both `BC::Dirichlet` (attr 5) AND
`BC::Fault` (attr 3).

The difference: **attr-3 boundary faces get K in Tandem but not in MFEM.**

But the critical question is: how many of the 9312 attr-3 boundary ELEMENTS
are actually true one-sided boundary faces (where `assemble_boundary` fires)?

The BP5 fault is entirely inside the domain (y=0 plane, x∈[-50,50], z∈[-40,0]).
So most attr-3 faces should be interior faces (two parent elements), not true
boundary faces. The number of attr-3 TRUE boundary faces is unknown — it could
be 0 (all interior) or nonzero (if the mesh topology creates some one-sided
fault faces at edges).

### 16.4 Tandem's perspective

Tandem classifies faces as boundary (`info.up[0] == info.up[1]`) based on its
own mesh topology, which may differ from MFEM's. In Tandem:

- `assemble_boundary` returns false for `BC::Natural` (attr 1) → skipped
- `assemble_boundary` returns true for `BC::Fault` (attr 3) → **K ADDED**
- `assemble_boundary` returns true for `BC::Dirichlet` (attr 5) → K added

If Tandem has `N` attr-3 boundary faces with K, then MFEM is missing `N`
boundary K contributions. The per-BC breakdown from Tandem (pending) will
give the exact number.

### 16.5 Interior/shared face discrepancy

587 faces that MFEM classifies as interior (both elements local) are classified
as shared (one element ghost) in Tandem. This is from different mesh
partitioners and is expected. The total skeleton face count differs by 587,
which could also be related to the boundary face classification difference.

### 16.6 Current assessment

**Confirmed so far:**
- Per-element volume K matches exactly (Section 15)
- Per-face skeleton K matches exactly (Section 15)
- RHS matches exactly (Section 13)
- Solver is correct (CG = MUMPS, Section 15.4)
- Volume element count matches (62874)

**Pending confirmation:**
- How many attr-3 or attr-1 faces does Tandem treat as boundary with K?
  (per-BC breakdown submitted, waiting for output)
- Are those the same faces MFEM misses?

**If Tandem's `BC::Fault` boundary count is > 0:** MFEM is missing boundary K
contributions on fault boundary faces. The fix is to add those faces to the
boundary face integrator marker. This would be the root cause of the
displacement difference that causes the 25 yr blowup.

**If Tandem's `BC::Fault` boundary count is 0:** The discrepancy is only in
attr-1 (Natural) faces, which both codes skip. Then the root cause is still
unknown and likely in the interior/shared face partitioning difference.

---

## 17. Files and Commits

### MFEM

- `miniapps/seas/domain/elasticity_operator.hpp` — RHS assembly
- `miniapps/seas/integrator/dg_elasticity_ip_combined_integrator.hpp` — traction + projection
- `miniapps/seas/solver/seas_operator.hpp` — Mult() coupling
- `miniapps/seas/trace/face_trace_logger.hpp` — tracer
- `miniapps/seas/tests/unit/test_cross_verify_tandem.cpp` — cross-verify tests
- `miniapps/seas/tests/unit/test_elasticity_operator.cpp` — normal leakage test

### Tandem

- `examples/tandem/3d/bp5_normal_diag.toml` — diagnostic probe config
- `jobs/run_bp5_1000m_p1_normal_diag.sbatch` — Frontera job script
- `app/localoperator/Elasticity.cpp` — DG assembly
- `app/localoperator/ElasticityAdapter.cpp` — traction projection
- `app/localoperator/DieterichRuinaBase.h` — friction solver

### Tandem comparison data

- Frontera job: 7634266 (8 nodes, 400 ranks, ran to 225 yr)
- Output: `/scratch2/10024/zhaochun/tandem/bp5_1000m_p1_normal_diag/fltst_*.csv`
- Key file: `fltst_tip_strk-37dp+40.csv` (normal-stress column)

### MFEM trace data (elastic σn)

- Frontera job: 7633584 (8 nodes, 400 ranks, aborted at 25.29 yr)
- Output: `results_v58_job7633584/trace_*.csv`
- Key DOF: r96 fi=28 dof=0 (owned_dof_idx=84, depth=40 km, x2=-37 km)

### MFEM corrected-decomposition trace data

- Trace data: `Archive_2026-04-06T16_44_22/trace_*_r96.csv`
- Run with fixed `ProjectTractionToFaultDOFs` basis-vector selection
- Decomposition closure verified to machine precision
- Key result: correction term is 84% of total normal traction at fault tip

### MFEM const-σn data

- Frontera job: 7634218 (8 nodes, 400 ranks, ran to 1800 yr)
- Output: `results_v58_const_sn_job7634218/bp5_v58_csn_job7634218_fltst_*.txt`
- Config: elastic_sigma_n disabled, everything else identical
- First earthquake: ~1191 yr (5.3× delayed vs Tandem's ~225 yr)

### Relevant commits

- `a329463` — fix double-loading of shared Dirichlet faces (fixed 0.5 yr blowup)
- `65d2e9d` — tag-only facet BC classification for BP5
- `8cb7c00` — accepted-step face tracing
- `53e5926` — per-DOF BP5 tracing diagnostics
- `da07d15` — BP5 geometry and traction regression tests
- `e7d6fa9` — job-specific output paths
- `2b3fbd1` — fix trace window depth sign
- Normal decomposition fix — `ProjectTractionToFaultDOFs` basis_vec selection (local, not yet committed)

---

## 16. K Matrix Investigation Roadmap

### 16.1 Established facts

- `b` (RHS) matches Tandem to 7 significant figures (Section 13.2)
- `u = K⁻¹·b` differs → K must differ (Section 13.5)
- Face integrator formula (`AssembleFaceMatrix`) verified line-by-line against
  Tandem's `assembleSurface` — no formula-level difference found
- The stress buildup is 5× too slow even with constant σ_n (Section 10)

### 16.2 K assembly call chain

```
AssembleStiffness() [elasticity_operator.hpp:2345]
│
├─ Volume: ElasticityIntegrator::AssembleElementMatrix()     [bilininteg.cpp:3208]
│   └─ ∫_Ω σ(u):ε(v) dV  — λ(div u)(div v) + 2μ ε(u):ε(v)
│
├─ Interior faces: DGElasticityIPCombinedIntegrator::AssembleFaceMatrix()
│                                    [dg_elasticity_ip_combined_integrator.hpp:45]
│   ├─ ComputePenalty()              [line 829]
│   └─ AssembleFaceBlock() × 4      [line 857]  — blocks (0,0),(0,1),(1,0),(1,1)
│       └─ 3 terms: consistency (c0), symmetry (c1), penalty (c2)
│
├─ Boundary faces: same DGElasticityIPCombinedIntegrator (1×1 block only)
│   └─ filtered by dirichlet_bdr_marker_
│
├─ BilinearForm::Assemble(0)        [bilinearform.cpp:456]
│   ├─ Domain loop (507-602)         — AssembleElementMatrix per element
│   ├─ Interior face loop (678-695)  — AssembleFaceMatrix per interior face
│   └─ Boundary face loop (724-749)  — AssembleFaceMatrix per bdr face
│
├─ ParBilinearForm::AssembleSharedFaces(0)  [pbilinearform.cpp:229]
│   └─ For each MPI-shared face → interior_face_integs->AssembleFaceMatrix()
│
├─ Finalize()
└─ ParallelAssemble() → HypreParMatrix
```

### 16.3 Verification status of each layer

| Layer | Component | Verified? | Method | Notes |
|-------|-----------|-----------|--------|-------|
| 1 | Face integrator formula | ✅ Yes | Line-by-line code comparison | Matches Tandem `assembleSurface` |
| 2 | Face integrator coefficients | ✅ Yes | Code comparison | c0/c1/c2 per block match Tandem |
| 3 | Penalty formula | ✅ Yes | Code comparison + unit test | `(D+1)*c_{N-1}*(D*nl_q/detJ)*(c1²/c0)` |
| 4 | DOF ordering (byNODES) | ✅ Yes | Code comparison | `row = p*ndof + k` consistent |
| 5 | Volume integrator | ❌ **No** | Never compared | MFEM's `ElasticityIntegrator` vs Tandem |
| 6 | Raw per-face elemmat values | ❌ **No** | Not dumped | Formulas match but numerics not checked |
| 7 | Face classification in K | ❌ **No** | Not verified | Interior vs boundary assignment |
| 8 | Shared-face K scatter | ❌ **No** | Not verified | DOF indexing for cross-partition faces |
| 9 | ParallelAssemble (P^T A P) | ❌ **No** | Not verified | HypreParMatrix construction |
| 10 | MUMPS interface | ❌ **No** | Assumed correct | Low suspicion |

### 16.4 Investigation steps (priority order)

#### Step 1: Volume integrator comparison (HIGHEST priority)

**Why first:** The volume integrator is the only K component never verified
against Tandem. It contributes the bulk of K. A wrong volume term would produce
the observed "too stiff / too compliant" behavior globally — consistent with
5× slow stress buildup.

**What to do:**
- Read MFEM's `ElasticityIntegrator::AssembleElementMatrix()` in
  `fem/bilininteg.cpp:3208`
- Read Tandem's volume assembly in `app/localoperator/Elasticity.cpp`
  (`assemble_volume` or equivalent)
- Compare: Lamé formulation, Jacobian handling, quadrature order, DOF layout
- Check whether Tandem uses Voigt notation vs component-wise and whether
  the resulting matrix is algebraically identical

**Key questions:**
- Does MFEM use `λ(div u)(div v) + 2μ ε:ε` or `λ(div u)(div v) + μ(∇u:∇v + ∇u^T:∇v)`?
  These are equivalent only for symmetric gradient.
- Does the quadrature order match? MFEM uses `2*OrderGrad(fe)`, Tandem uses
  `MinQuadOrder()`.
- Is the Jacobian determinant handled identically (reference → physical)?

#### Step 2: Dump raw per-face elemmat (HIGH priority)

**Why:** Even though formulas match, the actual numerical values could differ
due to MFEM-internal transform/Jacobian machinery. Dumping the raw face element
matrix before global scatter bypasses all parallel assembly issues.

**What to do:**
- In `DGElasticityIPCombinedIntegrator::AssembleFaceMatrix`, add a conditional
  dump of `elmat` for target faces (identified by element numbers or centroids)
- Dump as CSV: face_id, row, col, value
- Compare against Tandem's `A00/A01/A10/A11` blocks for matched faces
- This is purely local — no MPI issues

**What this proves:**
- If elemmat matches → bug is in assembly scatter, ParallelAssemble, or solve
- If elemmat differs → bug is in how MFEM evaluates transforms/Jacobians at
  face quadrature points, despite the formula being the same

#### Step 3: Dump raw per-element volume matrix (MEDIUM priority)

**Why:** If Step 1 finds a formula difference, this confirms it numerically.
If Step 1 finds no formula difference (like the face integrator), this catches
transform/Jacobian-level issues in the volume path.

**What to do:**
- In `ElasticityIntegrator::AssembleElementMatrix`, add conditional dump of
  `elmat` for target elements (those adjacent to fault-tip faces)
- Compare against Tandem's volume element matrix for matched elements

#### Step 4: Verify face classification consistency between K and b (MEDIUM)

**Why:** If a face is classified as interior in K assembly but as boundary
(or skipped) in RHS assembly (or vice versa), the K/b consistency breaks.
The formulas can be correct but K and b see different face sets.

**What to do:**
- During `BilinearForm::Assemble()`, log which faces are visited by each loop
  (domain, interior face, boundary face)
- Cross-reference against the face sets used in `AssembleSlipContributionIP()`
  and `AssembleDirichletLoading()`
- Check: every face that contributes to b_slip or b_dirichlet must also appear
  in the interior face or boundary face K loop with the matching integrator type

#### Step 5: Verify shared-face K assembly (MEDIUM-LOW)

**Why:** Shared faces are assembled by `ParBilinearForm::AssembleSharedFaces()`
using the interior face integrator. If the elem1/elem2 assignment or DOF
indexing differs from what the RHS path uses, K would be inconsistent.

**What to do:**
- Log which shared faces are visited during K assembly
- Verify elem1/elem2 assignment matches the RHS shared-face paths
- Check that `keep_nbr_block=true` is set (it is, per previous review)

#### Step 6: Serial reproduce (LOW priority, high diagnostic value)

**Why:** Running on 1 MPI rank eliminates all shared-face, ParallelAssemble,
and cross-partition issues. If the bug persists on 1 rank, it is in the
local assembly (volume or interior face integrator). If it disappears,
it is in the parallel assembly path.

**What to do:**
- Build a small mesh (e.g., 5 km × 5 km × 5 km) that fits on 1 rank
- Run both MFEM and Tandem single-rank
- Compare K·e_hat for a few unit vectors

### 16.5 Decision tree

```
Step 1: Volume integrator formula comparison
  ├─ Formula differs → FIX IT → likely root cause
  └─ Formula matches →
      Step 2: Dump raw per-face elemmat
        ├─ Face elemmat differs → MFEM transform/Jacobian bug
        └─ Face elemmat matches →
            Step 3: Dump raw volume elemmat
              ├─ Volume elemmat differs → MFEM volume transform bug
              └─ Volume elemmat matches →
                  Step 4: Face classification audit
                    ├─ Inconsistency found → FIX face set mismatch
                    └─ All consistent →
                        Step 5: Shared-face K audit
                          ├─ Issue found → FIX shared-face DOF scatter
                          └─ All correct →
                              Step 6: Serial vs parallel comparison
                                ├─ Serial matches Tandem → parallel assembly bug
                                └─ Serial also wrong → deeper MFEM framework issue
```

### 16.6 Expected outcome

The most likely outcome based on current evidence:

1. **Volume integrator** (40% probability): A subtle difference in how MFEM's
   `ElasticityIntegrator` handles the Lamé formulation, Jacobian, or quadrature
   compared to Tandem. This would explain why the displacement field is globally
   wrong (5× slow stress buildup) rather than locally wrong at specific faces.

2. **Face elemmat numerical difference** (30%): The formula matches but MFEM's
   `FaceElementTransformations` produces slightly different Jacobians or normals
   on specific face geometries (e.g., degenerate tets at the fault tip). This
   would explain the localization to fault-tip DOFs.

3. **Face classification mismatch** (15%): A face is processed by interior-face
   integrator in K but skipped (or processed as boundary) in b, or vice versa.
   This would break K/b consistency on specific faces.

4. **Parallel assembly issue** (10%): Shared-face DOF scatter or
   ParallelAssemble produces wrong entries for cross-partition faces.

5. **MUMPS interface** (5%): Matrix conversion issue. Very unlikely given that
   MUMPS is widely used and tested.

---

## 17. Step 1 Result: Volume Integrator Comparison — Formulas Match

### 17.1 MFEM volume integrator

**File:** `fem/bilininteg.cpp:3208-3288`
**Function:** `ElasticityIntegrator::AssembleElementMatrix()`

```cpp
// Per quadrature point:
el.CalcDShape(ip, dshape);                          // dshape[k,i] = ∂N_k/∂ξ_i
w = ip.weight * Trans.Weight();                      // w = w_q * |det(J)|
Mult(dshape, Trans.InverseJacobian(), gshape);       // gshape[k,i] = ∂N_k/∂x_i
MultAAt(gshape, pelmat);                             // pelmat[k,l] = Σ_i (∂N_k/∂x_i)(∂N_l/∂x_i)
gshape.GradToDiv(divshape);                          // divshape[j*dof+k] = gshape[k,j] = ∂N_k/∂x_j

// Term 1 — λ(div u)(div v):
AddMult_a_VVt(L*w, divshape, elmat);
// → elmat(p*dof+k, u*dof+l) += λ*w * (∂N_k/∂x_p)(∂N_l/∂x_u)

// Term 2 — μ(∇u·∇v) on diagonal blocks:
elmat(dof*d+k, dof*d+l) += M*w * pelmat(k,l);
// → elmat(p*dof+k, p*dof+l) += μ*w * Σ_j (∂N_k/∂x_j)(∂N_l/∂x_j) * δ_{pu}

// Term 3 — μ(∂u_i/∂x_j)(∂v_j/∂x_i) cross terms:
elmat(dof*ii+kk, dof*jj+ll) += M*w * gshape(kk,jj) * gshape(ll,ii);
// → elmat(p*dof+k, u*dof+l) += μ*w * (∂N_k/∂x_u)(∂N_l/∂x_p)
```

Combined, the MFEM element matrix entry at (p·dof+k, u·dof+l) is:

```
A[p,k,u,l] = Σ_q w_q |J_q| * [
    λ * (∂N_k/∂x_p)(∂N_l/∂x_u)                        // div-div
  + μ * Σ_j (∂N_k/∂x_j)(∂N_l/∂x_j) * δ_{pu}          // grad·grad (diag)
  + μ * (∂N_k/∂x_u)(∂N_l/∂x_p)                        // transpose coupling
]
```

DOF layout: **byNODES** — DOF index = component × ndof + basis_function.

### 17.2 Tandem volume integrator

**File:** `app/kernels/elasticity.py:55-56`
**Kernel:** `assembleVolume`

```python
A['kplu'] <= lam_W_J_Q['q'] * Dx_Q['luq'] * Dx_Q['kpq']
           + mu_W_J_Q['q'] * Dx_Q['kjq'] * (Dx_Q['ljq'] * delta['pu']
                                            + Dx_Q['lpq'] * delta['ju'])
```

Where:
- `Dx_Q[k,i,q] = Σ_e JInv[e,i,q] * Dxi_Q[k,e,q]` = ∂N_k/∂x_i at quad point q
- `lam_W_J_Q[q] = λ * w_q * |det(J)|_q` (precomputed)
- `mu_W_J_Q[q] = μ * w_q * |det(J)|_q` (precomputed)
- `delta[p,u] = δ_{pu}` (Kronecker delta)

Expanding term by term:

```
A[k,p,l,u] = Σ_q [
    λ*w*|J| * (∂N_l/∂x_u)(∂N_k/∂x_p)                  // div-div
  + μ*w*|J| * Σ_j (∂N_k/∂x_j)(∂N_l/∂x_j) * δ_{pu}    // grad·grad (diag)
  + μ*w*|J| * (∂N_k/∂x_u)(∂N_l/∂x_p)                  // transpose coupling
              ↑ from Σ_j Dx_Q[k,j]*Dx_Q[l,p]*δ_{ju} = Dx_Q[k,u]*Dx_Q[l,p]
]
```

DOF layout: tensor `A[k,p,l,u]` — DOF index = basis_function + Nbf × component.

**Called from:** `Elasticity::assemble_volume()` in `app/localoperator/Elasticity.cpp:525-573`

### 17.3 Term-by-term comparison

| Term | Physics | MFEM | Tandem | Match? |
|------|---------|------|--------|--------|
| λ | div-div | λ·w·\|J\|·(∂N_k/∂x_p)(∂N_l/∂x_u) | λ·w·\|J\|·(∂N_l/∂x_u)(∂N_k/∂x_p) | ✅ Same (commutative) |
| μ diag | grad·grad | μ·w·\|J\|·Σ_j(∂N_k/∂x_j)(∂N_l/∂x_j)·δ_{pu} | μ·w·\|J\|·Σ_j(∂N_k/∂x_j)(∂N_l/∂x_j)·δ_{pu} | ✅ Identical |
| μ cross | transpose | μ·w·\|J\|·(∂N_k/∂x_u)(∂N_l/∂x_p) | μ·w·\|J\|·(∂N_k/∂x_u)(∂N_l/∂x_p) | ✅ Identical |
| Weight | — | ip.weight × Trans.Weight() | W[q] × J[q] (precomputed) | ✅ Same |
| Gradient | — | dshape × J⁻¹ | JInv × Dxi_Q | ✅ Same |

All three terms are algebraically identical.

### 17.4 DOF layout difference (irrelevant)

- MFEM: `elmat(p*dof+k, u*dof+l)` → DOF = component × Ndof + basis
- Tandem: `A[k,p,l,u]` → DOF = basis + Nbf × component

This is a permutation of rows/columns, not a formula difference. The global
assembly uses DOF maps (vdofs) that account for the local ordering, so the
physics is identical regardless of local DOF numbering.

### 17.5 Quadrature order difference (irrelevant for p=1 affine tets)

**MFEM:** `order = 2 * Trans.OrderGrad(&el)`

For `IsoparametricTransformation::OrderGrad` with `FunctionSpace::Pk`:
```cpp
return (k-1)*(d-1) + (l-1);
// k = mesh order = 1, d = dim = 3, l = FE order = 1
// → (1-1)*(3-1) + (1-1) = 0
```
So MFEM uses `order = 2*0 = 0` → 1-point rule on tet.

**Tandem:** `MinQuadOrder() = 2p+1 = 3` → 4- or 5-point rule on tet.

**Why it doesn't matter:** For p=1 on affine (straight-sided) tets:
- Shape function gradients are **constant** (degree 0)
- Jacobian determinant is **constant** (affine map)
- Lamé parameters are **constant** (uniform material in BP5)
- The integrand is a **constant** → any quadrature rule with ≥ 1 point is exact

Both 1-point and multi-point rules integrate constants exactly on simplices.

### 17.6 Conclusion

**The volume integrator formulas match exactly.** The volume integrator is
**not the cause** of the K matrix mismatch.

This eliminates suspect #1 from the roadmap (Section 16.4, Step 1). The
investigation should proceed to **Step 2: dump raw per-face elemmat** to
check whether the face integrator produces numerically identical element
matrices despite having identical formulas.

### 17.7 Updated probability ranking

With the volume integrator eliminated:

1. **Highest (40%)**: Face elemmat numerical difference — same formulas but
   MFEM's `FaceElementTransformations` produces different Jacobians/normals
   on specific face geometries
2. **High (25%)**: Face classification mismatch — a face processed by wrong
   integrator path in K vs b
3. **Medium (20%)**: Shared-face K assembly — DOF scatter or elem1/elem2
   assignment error for cross-partition faces
4. **Medium-Low (10%)**: ParallelAssemble / HypreParMatrix construction
5. **Low (5%)**: MUMPS interface

---

## 18. Steps 2-3 Result: K Matrix Numerical Match — K Is Identical

### 18.1 Setup

Per-element and per-face K matrix entries were dumped in both codes and compared
numerically, using the face matching from Section 13.1.

**Data:**
- MFEM: Frontera job 7636174, `results_v59_fs_revert_job7636174/first_step_K_*.csv`
- Tandem: Frontera job 7636202, `results_tandem_1000m_p1_job7636202/first_step_K_*.csv`

**Element matching** (established via face geometry + volume K diagonal):
- MFEM elem 4 ↔ Tandem elem 73 (basis perm: {0:2, 1:3, 2:1, 3:0})
- MFEM elem 94 ↔ Tandem elem 45 (basis perm: {0:3, 1:2, 2:1, 3:0})
- MFEM elem 48 ↔ Tandem elem 129 (basis perm: {0:3, 1:0, 2:1, 3:2})
- MFEM elem 85 ↔ Tandem elem 8 (basis perm: {0:0, 1:3, 2:1, 3:2})

**Face matching** (from Section 13.1):
- MFEM face 17 (e1=4, e2=94) ↔ Tandem fct 146 (e0=45, e1=73)
- MFEM face 144 (e1=48, e2=85) ↔ Tandem fct 31 (e0=8, e1=129)

**Element-side mapping** (MFEM elem1 ↔ Tandem elem1, NOT elem0):
- MFEM A00 (elem1 side) ↔ Tandem A11 (elem1 side)
- MFEM A11 (elem2 side) ↔ Tandem A00 (elem0 side)
- MFEM A01 ↔ Tandem A10, MFEM A10 ↔ Tandem A01

**Unit scaling:** K_MFEM / K_Tandem = 1e12 (MFEM in meters, Tandem in km).

### 18.2 Volume K results

| Element pair | Entries compared | Matched | Max reldiff |
|-------------|-----------------|---------|-------------|
| MFEM 4 ↔ Tandem 73 | 126 | 126 | 7.8e-14 |
| MFEM 94 ↔ Tandem 45 | 126 | 126 | < 1e-12 |

All 252 volume K entries match to **machine precision** with a uniform
scaling factor of 1e12.

### 18.3 Face K results

| Face pair | Block | Exact matches | Machine-eps entries | Significant mismatches |
|-----------|-------|---------------|--------------------|-----------------------|
| face17 ↔ fct146 | A00→A11 | 81 | 0 | 0 |
| face17 ↔ fct146 | A11→A00 | 81 | 0 | 0 |
| face17 ↔ fct146 | A01→A10 | 69 | 6 (max ~7e-4) | 0 |
| face17 ↔ fct146 | A10→A01 | 69 | 6 (max ~7e-4) | 0 |
| face144 ↔ fct31 | A00→A11 | 87 | 20 (max ~2e-2) | 0 |
| face144 ↔ fct31 | A11→A00 | 87 | 16 (max ~1e-2) | 0 |
| face144 ↔ fct31 | A01→A10 | 81 | 25 (max ~5e-2) | 0 |
| face144 ↔ fct31 | A10→A01 | 81 | 25 (max ~5e-2) | 0 |

All "machine-eps entries" have values 14+ orders of magnitude smaller than the
matrix scale (~1e12). They represent off-face basis function interactions at
floating-point noise level. **Zero significant mismatches.**

### 18.4 Conclusion: K is identical

The stiffness matrix K is **numerically identical** between MFEM and Tandem to
double-precision machine epsilon. This holds for both:
- Volume element matrices (∫σ:ε dΩ)
- Face element matrices (consistency + symmetry + penalty terms)

### 18.5 Implications: the paradox

Combined with the prior result that **b matches to 7 significant figures**
(Section 13.2), we now have:

```
K matches ✅  +  b matches ✅  +  MUMPS direct solver  →  u MUST match
```

But the simulation produces a 5× slow stress buildup (Section 10). Therefore:

1. **The first-step u comparison (Section 13.3) may have been misleading.**
   It compared MFEM at t=0.002s with Tandem at t=0.02s (10× time difference)
   and concluded u_y was "~2× too large" after extrapolation. But if K and b
   match globally, u = K⁻¹·b must also match. The apparent u discrepancy may
   reflect comparison methodology (different times, unit conversion, face
   matching) rather than a real difference.

2. **The bug is NOT in the linear solve (K·u = b).** Since K, b, and the
   solver are all verified, the displacement u at any given time step is
   correct given the inputs (slip, Dirichlet loading) at that step.

3. **The bug must be in the coupling loop** — how the solved u feeds back
   through traction → friction → slip → next b:
   - Traction computation: T(u) on fault faces
   - L2 projection: QP traction → fault DOF traction
   - Friction law: V(T, ψ) → slip rate
   - State evolution: dψ/dt
   - Slip integration: how V·dt accumulates into slip
   - Slip → RHS: how cumulative slip enters b_slip at the next step

4. **Or in the face/element counting** — if MFEM assembles K and b over a
   different set of faces/elements than Tandem (e.g., extra boundary faces,
   missing interior faces), the global matrix could differ even though
   individual face contributions match.

### 18.6 Revised investigation direction

The K investigation roadmap (Section 16) is now **complete through Step 3**.
Steps 1-3 all show matches. The remaining roadmap steps (4-6) address assembly
scatter and parallel issues, but these are less likely given that the K formula,
volume elemmat, and face elemmat all match.

The investigation must shift to the **coupling loop**:

**Priority 1: Verify face/element counts**
- How many volume elements, interior faces, shared faces, and boundary faces
  does each code process during K assembly?
- Are these counts identical? (Commit 2c39032 added diagnostics for this.)
- A count mismatch would mean K or b is assembled over different face sets.

**Priority 2: Traction computation comparison**
- Dump T(u) at fault faces from both codes at the first step
- Compare against Tandem's traction output
- A mismatch here (with matching u) would point to the traction formula or
  its inputs (normals, Jacobians, penalty)

**Priority 3: Friction coupling audit**
- How does the traction from the elastic solve enter the friction law?
- Are the sign conventions consistent? (compression-positive σ_n)
- Is the slip deficit computed correctly?
- How does the ODE integrator advance slip and state?

**Priority 4: Time-stepping feedback test**
- Fix the slip to a known analytical profile (no friction feedback)
- Compare u, traction, and normal stress between the two codes
- If they match with fixed slip, the bug is in the friction coupling

---

## 19. ~~ROOT CAUSE FOUND~~: INVALIDATED — Tandem Per-BC Data Shows No Missing Fault Faces

### 19.1 Face count discrepancy

K assembly diagnostic counts (global, MPI_Allreduce):

| Quantity | MFEM | Tandem | Diff |
|----------|------|--------|------|
| Volume elements | 62,874 | 62,874 | 0 |
| Interior + Shared faces | 142,624 | 143,211 | -587 |
| Boundary faces with K | 1,628 | 2,254 | -626 |
| **Total K contributions** | **207,126** | **208,339** | **-1,213** |

MFEM is missing **~1,213 K contributions** compared to Tandem.

### 19.2 The bug

**Tandem** (`Elasticity::assemble_boundary`, line 724):
```cpp
if (info.bc == BC::Natural) { return false; }
// Adds boundary K (c0=-1, c1=ε, c2=penalty) for ALL non-Natural boundary faces
// This includes both Dirichlet (attr 5) AND Fault (attr 3) boundary faces
```

**MFEM** (`SetupBoundaryMarkers`, line 1350):
```cpp
if (attr == 5) { dirichlet_bdr_marker_[i] = 1; }
// Only attr 5 → boundary K only for Dirichlet faces
```

**MFEM** (`AssembleStiffness`, line 2528):
```cpp
cached_a_->AddBdrFaceIntegrator(
    new DGElasticityIPCombinedIntegrator(...),
    dirichlet_bdr_marker_);  // Only fires for attr 5
```

One-sided fault faces (attr 3 on the mesh surface — e.g., where the fault
plane meets z=0 top surface) get **no boundary K contribution** in MFEM:
- No penalty enforcement
- No consistency term
- No symmetry term

Without these terms, the DG formulation is incomplete on these faces. The
displacement is effectively under-constrained on one-sided fault faces,
allowing spurious discontinuity to develop.

### 19.3 Why this explains the symptoms

1. **5× slow stress buildup**: The missing boundary K changes the global
   stiffness matrix. Since K is used for ALL time steps (cached), every
   elastic solve produces a slightly wrong displacement field. The error
   accumulates through the coupling loop.

2. **Fault-tip localization**: The one-sided fault faces are at the edges of
   the fault rectangle — exactly where the fault meets the mesh surface at
   z=0 and possibly at x=±Lf. The fault-tip DOFs at depth=40 km are coupled
   to these surface fault faces through the global K matrix.

3. **Linear growth rate**: The Dirichlet loading grows as V_p·t. Each time
   step, the slightly wrong K produces a displacement with slightly wrong
   stress. Since the loading is linear in t, the error is linear in t.

### 19.4 Proposed fix

Add a `boundary_k_marker_` that includes both attr 3 (Fault) and attr 5
(Dirichlet). Use it for the boundary face K integrator. Keep
`dirichlet_bdr_marker_` (attr 5 only) for the Dirichlet RHS loading.

```cpp
// New member
Array<int> boundary_k_marker_;  // attrs 3 + 5 for boundary K

// In SetupBoundaryMarkers:
boundary_k_marker_.SetSize(num_bdr);
boundary_k_marker_ = 0;
for (int i = 0; i < num_bdr; i++) {
    int attr = i + 1;
    if (attr == 3 || attr == 5) { boundary_k_marker_[i] = 1; }
}

// In AssembleStiffness:
cached_a_->AddBdrFaceIntegrator(
    new DGElasticityIPCombinedIntegrator(...),
    boundary_k_marker_);  // Now fires for attr 3 AND attr 5
```

### 19.5 Related concern: boundary slip RHS

One-sided fault boundary faces also need their slip RHS contribution.
Currently, `AssembleSlipContributionIP` processes fault faces from
`fault_faces_` (interior) and `shared_fault_faces_` (shared). If one-sided
boundary fault faces are not in either list, they are missing their slip RHS
contribution as well.

For boundary fault faces, the slip RHS should use **boundary coefficients**
(c1=ε, c2=penalty_single), not interior coefficients (c1=0.5ε, c2=±penalty).
This matches Tandem's `rhs_boundary` which uses the same coefficients for
both Dirichlet and Fault boundary faces.

This requires a separate check: are one-sided fault faces in `fault_faces_`
or `shared_fault_faces_`? If not, a corresponding RHS fix is also needed.

### 19.6 INVALIDATED by Tandem per-BC breakdown (job 7636327)

Tandem's per-BC boundary face breakdown:

```
BC::None:                0
BC::Natural (skipped):   2014
BC::Fault (K added):     0
BC::Dirichlet (K added): 240
Bdr faces with K:        240
```

**BC::Fault = 0.** There are NO one-sided fault boundary faces in Tandem's
topology. All 9312 attr-3 boundary elements are interior (two-sided) faces.

**The actual comparison is reversed:**

| Quantity | MFEM | Tandem | Diff |
|----------|------|--------|------|
| Bdr faces with K | **1628** | **240** | **MFEM +1388** |
| Interior+Shared | 142624 | 143211 | MFEM -587 |
| Total K contributions | **207126** | **206325** | **MFEM +801** |

MFEM has **801 MORE** total K contributions than Tandem, not fewer.

The 1628 vs 240 boundary K discrepancy: MFEM processes 1628 attr-5 boundary
faces with K, while Tandem only processes 240 `BC::Dirichlet` faces. The
difference of 1388 could be attr-5 boundary elements that MFEM treats as
true boundary faces but Tandem treats as interior (skeleton) faces.

**The Section 19.1-19.5 root cause hypothesis is WRONG.** The bug is not
"missing boundary K on fault faces." The actual K discrepancy is:
- MFEM has more boundary K (1388 extra attr-5 faces)
- MFEM has fewer skeleton K (587 fewer interior faces)
- Net: MFEM has 801 more K contributions

**The root cause remains unknown.** The per-element K matrices match exactly
(Section 15), the RHS matches exactly (Section 13), but the displacement
differs. The K contribution COUNTS differ by ~801, suggesting some faces
are classified differently between the two codes' mesh topologies.

### 19.7 Revised next steps

1. Understand why MFEM has 1388 more attr-5 boundary K faces than Tandem
   (1628 vs 240). Are 1388 faces treated as boundary in MFEM but skeleton
   in Tandem?

2. Understand the 587 interior face difference. Combined with the boundary
   difference: if 1388 faces move from skeleton→boundary in MFEM, and
   587 faces move from shared→interior, the topology is significantly
   different.

3. Consider whether the MFEM boundary face integrator formula (c1=ε,
   c2=penalty_single) vs the skeleton formula (c1=ε/2, c2=penalty_avg)
   matters for the 1388 faces that have different classification.
