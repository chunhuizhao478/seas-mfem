# BP5 Debug v58: First-Solve MFEM vs Tandem Status and Ranked Bug Candidates

**Date:** 2026-04-04  
**Status:** Historical first-solve note. Later work fixed the earlier `~0.5 yr` blow-up via shared-Dirichlet RHS de-duplication; current remaining failure is a later `~2.5 yr` blow-up.  
**Scope:** BP5, 1000 m mesh, IP, p=1, PETSc TS RK45, exact Tandem comparison case

---

## 1. Problem Statement

Tandem runs the exact BP5 `p=1`, `1000 m` case without the traction blow-up seen in MFEM
around `t ≈ 0.5 yr`. Therefore MFEM still has a real bug or inconsistency in the BP5
coupling path, even though several earlier suspects have now been eliminated.

The current goal of the debug work is narrower than “explain all late-time divergence”:

1. compare MFEM and Tandem on the same physical fault face during the first solve
2. identify what is already matched
3. identify what is already different
4. rank the remaining code locations that could still create the later traction blow-up

---

## 2. What Data We Have From the First-Solve Comparison

### 2.1 Tandem exact-face dump

Tandem exact matched face:

- face key: `key=(63,3300,3910)`
- dump time: `t = 2.0000000000e-02`
- penalty: `2.1833007129e+03`
- area: `4.0958328773e-01`
- volumes: `vol0=1.5081597266e-01`, `vol1=1.4970283849e-01`

Per-QP Tandem data:

| q | xyz (km) | ny | u0_y | u1_y | slip_y | jump_y | Ty_stress | Ty_penalty | Ty |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | `(-49.2162, 0, -2.1949)` | `-1` | `-7.5228e-07` | `-7.6581e-07` | `-5.2923e-27` | `1.3530e-08` | `2.9987e-05` | `-2.9540e-05` | `4.4686e-07` |
| 1 | `(-49.3845, 0, -2.2913)` | `-1` | `-1.3932e-06` | `-1.4077e-06` | `-2.0000e-22` | `1.4416e-08` | `3.1507e-05` | `-3.1475e-05` | `3.2752e-08` |
| 2 | `(-49.3199, 0, -1.7653)` | `-1` | `-1.2949e-06` | `-1.3477e-06` | `2.0000e-11` | `5.2834e-08` | `1.1514e-04` | `-1.1535e-04` | `-2.1136e-07` |
| 3 | `(-49.7203, 0, -1.9946)` | `-1` | `0.0` | `0.0` | `-5.2923e-27` | `5.2923e-27` | `-1.1968e-06` | `-1.1555e-23` | `-1.1968e-06` |

### 2.2 MFEM exact-face dumps

MFEM matched face for the same physical location:

- face key: `key=(551,555,613)`
- centroid: `cx=-4.9440115035e+04`, `cz=2.0121236602e+03`
- penalty: `2.1833007129e+09`
- area: `4.0958328773e+05`
- volumes: `vol0=8.9821703093e+08`, `vol1=9.0489583596e+08`

We now have three MFEM exact-face snapshots for this same face:

1. **Old PETSc stage-state dump**  
   - `step=1`, `t=2.0e-03`, `dt=1.0e-02`  
   - this was the bad diagnostic path later fixed in `e3d4d31`

2. **Accepted first-step MFEM dump**  
   - `step=1`, `t=1.0e-02`, `dt=1.0e-02`

3. **Accepted post-crossing MFEM dump**  
   - `step=2`, `t=2.8748438534e-02`, `dt=1.8748438534e-02`

The accepted-state dumps are the important ones. They show that the exact-face mismatch
survives across the early time window that brackets Tandem’s `t=2.0e-02` dump.

#### 2.2.1 Accepted first-step MFEM dump (`t=1.0e-02`)

| q | xyz (m) | ny | u0_y | u1_y | slip_y | jump_y | Ty_stress | Ty_penalty | Ty |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | `(-4.9440e+04, 0, -2.0121e+03)` | `+1` | `-9.0916e-07` | `-9.0912e-07` | `2.6462e-27` | `-4.0484e-11` | `9.5424e-02` | `8.8388e-02` | `1.8381e-01` |
| 1 | `(-4.9321e+04, 0, -2.2098e+03)` | `+1` | `-9.2230e-07` | `-9.2238e-07` | `2.6462e-27` | `7.6294e-11` | `9.5424e-02` | `-1.6657e-01` | `-7.1149e-02` |
| 2 | `(-4.9336e+04, 0, -1.8193e+03)` | `+1` | `-9.0978e-07` | `-9.0977e-07` | `2.6462e-27` | `-5.7769e-12` | `9.5424e-02` | `1.2613e-02` | `1.0804e-01` |
| 3 | `(-4.9664e+04, 0, -2.0073e+03)` | `+1` | `-8.9539e-07` | `-8.9520e-07` | `2.6462e-27` | `-1.9197e-10` | `9.5424e-02` | `4.1913e-01` | `5.1455e-01` |

#### 2.2.2 Accepted second-step-crossing MFEM dump (`t=2.8748438534e-02`)

| q | xyz (m) | ny | u0_y | u1_y | slip_y | jump_y | Ty_stress | Ty_penalty | Ty |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | `(-4.9440e+04, 0, -2.0121e+03)` | `+1` | `-2.5940e-06` | `-2.5939e-06` | `7.6074e-27` | `-1.1528e-10` | `2.7503e-01` | `2.5168e-01` | `5.2672e-01` |
| 1 | `(-4.9321e+04, 0, -2.2098e+03)` | `+1` | `-2.6315e-06` | `-2.6317e-06` | `7.6074e-27` | `2.1809e-10` | `2.7503e-01` | `-4.7616e-01` | `-2.0113e-01` |
| 2 | `(-4.9336e+04, 0, -1.8193e+03)` | `+1` | `-2.5958e-06` | `-2.5957e-06` | `7.6074e-27` | `-1.6328e-11` | `2.7503e-01` | `3.5650e-02` | `3.1068e-01` |
| 3 | `(-4.9664e+04, 0, -2.0073e+03)` | `+1` | `-2.5547e-06` | `-2.5542e-06` | `7.6074e-27` | `-5.4759e-10` | `2.7503e-01` | `1.1956e+00` | `1.4706e+00` |

---

## 3. What Is Matched

The following items are now matched strongly enough that they are no longer the top suspects.

### 3.1 Same physical face

The earlier “wrong face” problem is solved.

- Tandem exact face: `key=(63,3300,3910)`
- MFEM exact face: `key=(551,555,613)`

These are different numbering systems, but they identify the same physical face.

### 3.2 Geometry and penalty scaling

After unit conversion, the face geometry is consistent:

- face area matches
- face volumes are consistent
- penalty matches exactly up to the expected `Pa` vs `MPa` scaling:
  - Tandem: `2.1833007129e+03`
  - MFEM: `2.1833007129e+09`

### 3.3 D8 slip-sign chain

The explicit D8 sign-chain unit test now verifies that MFEM and Tandem produce the same
prescribed DG jump sign for both `sign_flipped=false` and `sign_flipped=true`.

This is covered by the added assertions in:

- `miniapps/seas/tests/unit/test_cross_verify_tandem.cpp`

and the test passes:

- `seas_test_cross_verify`: `119 passed, 0 failed`

This lowers the probability that the blow-up comes from the BP5 slip sign convention itself.

### 3.4 Shared fault face MPI ownership bug from v57

The old one-sided shared fault face detection bug was already fixed in v57 (`932b797`).
That bug explained MPI core-count dependence, but it does not explain the remaining
MFEM-vs-Tandem difference for the exact `p=1`, `1000 m` case.

### 3.5 Production jump formula: line-by-line MFEM-vs-Tandem match

A line-by-line code comparison of the traction-extraction jump path confirms the production
formula is identical between MFEM and Tandem. Examined locations:

- MFEM: `dg_elasticity_ip_combined_integrator.hpp:653-667`
- Tandem kernel: `app/kernels/elasticity.py:242-244`
- Tandem C++ debug: `app/localoperator/Elasticity.cpp:1088-1097`

#### 3.5.1 Jump formula: exact match

Both codes compute: `jump = u⁻ - u⁺ - f_q`

| | MFEM | Tandem |
|---|---|---|
| u⁻ reconstruction | `Σ_k shape1(k) * u1_dofs(c*ndof1+k)` | `Σ_l E_q[0](l,q) * u[0](l,p)` |
| u⁺ reconstruction | `Σ_k shape2(k) * u2_dofs(c*ndof2+k)` | `Σ_l E_q[1](l,q) * u[1](l,p)` |
| slip | `slip_3d(c*nq+q)` | `f_q(p,q)` |
| stored jump | `u1q - u2q - f_q` | `u0_y - u1_y - slip_y` |

Storage layouts differ (MFEM: component-major `u[c*ndof+k]`; Tandem: basis-major `u[l*Dim+c]`)
but dot products are algebraically identical.

#### 3.5.2 Stress average: exact match

Both compute `T_p = Σ_j 0.5*(σ1_pj + σ2_pj) * n̂_j` using:

- Per-side material properties (MFEM: `lam1t/mu1t`, `lam2t/mu2t`; Tandem: `lam_q[0]`, `lam_q[1]`)
- Unit normal (MFEM: `n_hat = nor/nl`; Tandem: `n_unit_q`)
- Same Cauchy traction: `σ·n̂ = λ·div(u)·n̂ + μ·(∇u+∇uᵀ)·n̂`

MFEM code: lines 630-641.  Tandem Python: `traction(x, normal)` at line 89-91.

#### 3.5.3 Penalty correction sign: exact match

- MFEM: `corr = (-penalty) * (u1q - u2q - f_q)` → added to traction (line 662)
- Tandem: `krnl.c00 = -penalty(fctNo)` (line 977), kernel adds `c00 * jump` (line 244)

Both: `T += -penalty * (u⁻ - u⁺ - f_q)`.

#### 3.5.4 Trace inverse inequality constant: exact match

- Tandem: `InverseInequality<D>::trace_constant(PolynomialDegree-1)` = `(N+1)(N+D)/D` with `N=p-1`
  → `p*(p+D-1)/D`  (`src/form/InverseInequality.h:27-28`)
- MFEM: `p*(p+dim_-1)/dim_` (line 830)

Identical formula.

#### 3.5.5 Normal convention: consistent

- MFEM: `CalcOrtho(Trans.Jacobian(), nor)` → outward from elem1; jump = elem1 − elem2
- Tandem: `cl_->normal(info.localNo[0], ...)` → outward from `up[0]`; jump = up[0] − up[1]

Both: normal points from the "minus" side (side 0) outward.

#### 3.5.6 Area/volume ratio in penalty: matches for affine tets

- Tandem: `area_[fctNo] / volume_[info.up[side]]` (precomputed integrals)
- MFEM: `dim_ * nl_q / detJ` (per-quad-point Jacobian ratio)

For affine tets with MFEM reference geometry (`V_ref=1/6`, face `S_ref=1/2`):
`dim*nl_q/detJ = 3*(area/S_ref)/(vol/V_ref) = area/volume`. Identical for straight-sided
elements. Would diverge for curved elements (MFEM varies per QP, Tandem uses a single
precomputed value), but BP5 uses affine tets only.

#### 3.5.7 Penalty material mismatch (latent, harmless for BP5)

One real code-level mismatch found in `ComputePenalty`:

**MFEM** (`dg_elasticity_ip_combined_integrator.hpp:647`):
```
penalty = ComputePenalty(fe1, fe2, detJ1, detJ2, lam1t, mu1t, nl, true);
//                                                 ^^^^  ^^^^
//                                                 element-1 material only
```

Both `p0` and `p1` inside `ComputePenalty` use the same `ratio = c1²/c0` derived exclusively
from element 1's λ, μ (lines 831-833).

**Tandem** (`Elasticity.cpp:285-288`):
```
p(side) = (D+1) * c_N_1 * (area/volume[side]) * (c1[side]²/c0[side])
```

Each side uses `stiffness_tensor_bounds(info.up[side])` — the material properties of its own
element.

**Impact for BP5**: None. BP5 is homogeneous (same λ, μ everywhere), so both sides always
have identical material. The penalty values match exactly, as confirmed by the first-solve
dump comparison (Section 2, penalty = `2.1833007129e+03` / `2.1833007129e+09`).

**Impact for heterogeneous-material problems**: Incorrect. MFEM would under- or over-penalize
on faces between materials with different stiffness. This should be fixed if MFEM is later
used for heterogeneous elasticity.

### 3.6 Sign-chain cleanup: state → slip → DG f_q

The full state-to-DG-slip path was audited and cleaned up. The MFEM sign chain now
matches Tandem step-for-step, with one intentional implementation difference that is
algebraically equivalent.

#### 3.6.1 Steps that match exactly (same sign, same term)

| Step | Tandem | Line | MFEM | Line |
|------|--------|------|------|------|
| V direction | `-(V/\|τ\|)·τ` | `DieterichRuinaBase.h:174` | `-(τ·inv_tau)·V_abs` | `dieterich_ruina.hpp:458` |
| ODE RHS | `r_mat(node,t) = Vi[t]` | `RateAndState.h:234` | `rate(i*SPN+c) = V_vec[c]` | `rate_state_fault.hpp:530` |
| State → slip | `krnl.slip = state.data()` | `ElasticityAdapter.cpp:46` | `slip = state` (via GetSlip) | `rate_state_fault.hpp:631` |
| DG formula | `c0*(E_q[0]*u[0]-E_q[1]*u[1]-f_q)` | `elasticity.py:243` | `(-pen)*(u1q-u2q-f_q)` | `combined_integrator.hpp:662` |
| Traction proj | `T·fault_basis_q` | `elasticity_adapter.py:26` | `T·basis_vec·sf` | `combined_integrator.hpp:778` |

Each of these steps can be read side-by-side with the same sign.

#### 3.6.2 Intentional implementation difference: sign_flipped handling

This is the one step where the code is **not** line-for-line identical, by design.

**The problem both codes solve:** The mesh may orient a face normal opposite to the
physical reference normal. When this happens, the DG formula's `u1-u2` reverses
(elem1/elem2 swap), so the prescribed slip `f_q` must also reverse to keep the
physical constraint `u⁻-u⁺ = f_q` correct.

**Tandem's approach** — mutate the stored basis (`AdapterBase.cpp:71-84`):

```
sign_flipped = n_ref_dot_n < 0;
if (sign_flipped) { normal_i = -1.0 * normal_i; }  // flip normal
cl_->facetBasis(up_, normal, fault_basis_q);         // compute basis from flipped normal
if (sign_flipped) { fault_basis_q *= -1.0; }         // negate entire basis
```

After this, `fault_basis_q` is permanently modified. Both `evaluate_slip` and
`evaluate_traction` use the (possibly negated) basis with no runtime sign:

```python
slip_q = e_q * fault_basis_q * slip * copy_slip      # adapter.py:21
traction = minv * e_q_T * ... * traction_q * fault_basis_q  # adapter.py:26
```

**MFEM's approach** — runtime scalar (`elasticity_operator.hpp`):

```cpp
real_t sign = basis.sign_flipped ? -1.0 : 1.0;  // lines 1761, 4335, 5147, 2177
...
delta_u_quad(c * nqp + q) = sign * du[c];        // slip embedding
...
sf = sign_flipped ? -1.0 : 1.0;                  // combined_integrator.hpp:755,765
T_local += traction_q(...) * basis_vec[t][p] * sf;  // traction projection
```

The tangent vectors `(t1, t2)` are stored as the true geometric frame (never negated).
The `sign` scalar applies the same `-1` at runtime instead.

**Why MFEM does not mutate the basis:**

MFEM's `FaultBasis` stores `(normal, tangent1, tangent2)` and is read by multiple
consumers beyond the DG path: diagnostic dumps (`[SLIP-EMBED]`, `[TIP-FACE]` prints),
per-QP coherence checks, and any future code that needs the physical tangent frame.
Negating the stored vectors would make them geometrically wrong (tangents pointing
opposite to dip/strike), which would corrupt diagnostics and violate the principle
that stored geometry should describe the physical face.

Tandem can mutate `fault_basis_q` because it is a private adapter tensor consumed only
by the two kernels. MFEM's `FaultBasis` is a shared data structure, so mutation is not safe.

**Algebraic equivalence proof:**

Let `B = [t1 | t2]` be the tangent matrix (columns = dip, strike in 3D).
Let `s = (s_dip, s_strike)` be the fault-local slip.

For slip embedding when `sign_flipped = true`:

- Tandem: `slip_q = (-B) · s = -(B · s)`
- MFEM: `delta_u = (-1) · (B · s) = -(B · s)`

For traction projection when `sign_flipped = true`:

- Tandem: `τ_local = T · (-B) = -(T · B)`
- MFEM: `τ_local = (T · B) · (-1) = -(T · B)`

Both use `(-B)·s = (-1)·(B·s)` (linearity of matrix multiplication). The `-1` lives
in different places (inside `B` vs outside as a scalar), but the product is identical.

Both directions use the same sign formula `sign_flipped ? -1 : +1`, matching Tandem's
single "negate the entire basis" operation.

#### 3.6.3 Output convention boundary

The internal API (`GetSlip()`) returns Tandem-internal sign (S anti-parallel to τ
for BP5). All SCEC output paths negate at the file-write boundary:

- `bp5_parallel_output.hpp:268-269`: `sd = -interp(slip_dip)`
- `bp5_benchmark_output.hpp:630-633`: `slip_dip = -interp(global_slip_dip)`
- `bp5_benchmark_output.hpp:850-853`: `slip_dip = -interp_interleaved(slip, 0)`

This converts to SCEC physical slip: `δ_SCEC = u⁺ − u⁻ = -S` (parallel to τ).

BP1/BP2 output paths are unaffected (scalar state is already parallel to τ).

#### 3.6.4 Test coverage for the sign chain

- `seas_test_cross_verify`: 119 passed — validates D8 sign chain end-to-end
- `seas_test_bp5_output`: 78 passed — validates SCEC output signs
- `seas_test_elasticity_operator`: 431 passed — validates explicit IP slip RHS matches
  production assembly (the explicit form uses the same `sign_flipped ? -1 : +1` convention)

### 3.7 Slip RHS raw-vs-lifted IP form

There is no bug in the specific “raw `f_q` vs lifted `f_lifted_q`” suspicion for IP.

Tandem IP kernel:

- `rhs_lift_ip: f_lifted_q = f_q * nl_q`
- `rhsFacet: ... + c2 * w * E_q * f_lifted_q`

MFEM IP combined integrator:

- `... + penalty * w_q * nl_q * shape * f_q`

Those are algebraically the same for the IP case.

---

## 4. What Is Different

The first-solve comparison still contains real differences. Some are expected, some are not.

### 4.1 The dump times are still not exactly matched

This is the single biggest limitation of the current first-solve comparison.

- MFEM accepted dump 1: `t=1.0e-02`
- MFEM accepted dump 2: `t=2.8748438534e-02`
- Tandem dump: `t=2.0e-02`

So the two exact-face dumps are **not from exactly the same physical time**.

Important note:

- the old MFEM exact-face dump was taken from a PETSc RK stage state, not the accepted step state
- this was a real diagnostic bug
- it was fixed in commit `e3d4d31`

This means the exact-face numbers are not the final apples-to-apples physics comparison.
However, they are now sufficient to answer a more limited question:

- does the MFEM-vs-Tandem early `jump_y` mismatch persist across the early accepted-step window?

The answer is yes.

### 4.2 QP ordering differs

MFEM and Tandem do not enumerate the quadrature points in the same order, and the face-local
reference mappings are not identical. Therefore `q=0` in one code should not be assumed to be
the same physical quadrature point as `q=0` in the other code.

The comparison must be done by physical `xyz`, not by raw `q` index.

### 4.3 Normal convention differs

- MFEM: `ny=+1`
- Tandem: `ny=-1`

This by itself is not a bug. It means signs in per-QP traction/jump fields must be interpreted
through each code’s convention, not compared naively.

### 4.4 `jump_y` is qualitatively different

This is the strongest remaining physical mismatch in the current first-solve data:

- MFEM accepted `t=1.0e-02`: `jump_y ~ 1e-11 .. 1e-10`
- MFEM accepted `t=2.8748e-02`: `jump_y ~ 1e-11 .. 1e-10`
- Tandem: `jump_y ~ 1e-08`

Even allowing for the unmatched times, this is a stable early-time difference, not a one-off
artifact of a single dump. This keeps the face-jump path near the top of the suspect list.

### 4.5 Stress/correction split is qualitatively different

MFEM exact-face dump:

- `Ty_stress` is nearly constant across the QPs
- `Ty_penalty` varies strongly

Tandem exact-face dump:

- `Ty_stress` varies more across QPs
- `Ty_penalty` tends to nearly cancel `Ty_stress`

Because the times are not yet matched, this difference cannot yet be called the root cause.
But it is still a strong warning sign that the two codes are not carrying the same face state
into traction recovery during the first solve.

---

## 5. What We Know the Current Data Can and Cannot Prove

### 5.1 What it can prove

The current debug data proves:

1. MFEM and Tandem are looking at the same physical face.
2. The face geometry and penalty scaling are consistent.
3. The BP5 slip sign chain is not obviously wrong.
4. The old MPI shared-face detection bug is not the current explanation.
5. The old MFEM exact-face diagnostic timing was wrong under PETSc TS.

### 5.2 What it cannot yet prove

The current data does **not** yet prove:

1. the exact same-time MFEM state at `t=2.0e-02`
2. that the exact-face traction kernel itself is wrong
3. whether the entire later 0.5 yr blow-up can be reduced to this early jump mismatch alone

So the current status is:

- there is still a real MFEM bug or inconsistency
- the exact same-time comparison would still be nice to have
- but it is no longer necessary to justify moving upstream into the jump path, because the
  early mismatch is already persistent across the accepted-step window bracketing `t=2.0e-02`

---

## 6. Ranked Potential Bug Locations

The following list is ordered from **highest current probability** to **lowest** for causing the
later traction blow-up around `0.5 yr`.

### 6.1 Highest probability: early fault-state / jump path before friction update

**Why it is still high probability**

The strongest remaining mismatch is still the tiny MFEM `jump_y` on the exact face. If that
tiny jump survives the accepted-step comparison, then MFEM is feeding the friction law a fault
state that is already different from Tandem near the very first solve.

**Code areas**

- `miniapps/seas/solver/seas_operator.hpp`
  - `SetInitialCondition()`
  - `Mult()`
- `miniapps/seas/domain/elasticity_operator.hpp`
  - `ExpandOwnedToLocalFault()`
  - `Solve()`
  - `ComputeTractionImpl()`

**Failure mode**

The wrong local fault slip or wrong side trace arrives at the domain solve / traction recovery,
so MFEM sees `[[u]] - slip` much closer to zero than Tandem does. That then perturbs the
quasi-dynamic friction update and grows into the 0.5 yr blow-up.

### 6.2 High probability: face trace / jump construction on the exact fault face

**Why it is still high probability**

The direct first-solve quantity that still looks wrong is `jump_y`.

**Code areas**

- `miniapps/seas/domain/elasticity_operator.hpp`
  - `ComputeTractionImpl()`
  - exact-face diagnostic block around the `jump_y = u0y - u1y - slip_y` print
- `miniapps/seas/integrator/dg_elasticity_ip_combined_integrator.hpp`
  - `ComputeTractionAtQuadPointsDecomposed()`
  - `AssembleSlipFaceRHS()`

**Failure mode**

MFEM may be reconstructing `u^-`, `u^+`, or the sign of the embedded slip on the exact face
in a way that is internally consistent in unit tests but still different from Tandem in the
real BP5 coupled solve.

### 6.3 Medium probability: accepted-step PETSc coupling path

**Why it matters**

This is already a proven diagnostic bug:

- MFEM exact-face dumps were previously taken from PETSc RK stage states, not accepted step states

That by itself does not explain the 0.5 yr blow-up, but it means the current exact-face
comparison is still not final.

**Code areas**

- `miniapps/seas/tests/verification/bp5_verification_full.cpp`
- PETSc stepping / post-step re-evaluation path

**Failure mode**

Not the physics bug itself, but this can hide or mislocate the physics bug by comparing the
wrong state.

### 6.4 Medium probability: full-system IP assembly vs slip RHS consistency in the real BP5 path

**Why it is lower than the jump path but still relevant**

The local unit tests strongly suggest the IP combined integrator is internally consistent.
However, the full BP5 solve still couples:

- global matrix assembly
- slip RHS assembly
- shared-face communication
- domain solve

in a way that is more complex than the toy tests.

**Code areas**

- `miniapps/seas/domain/elasticity_operator.hpp`
  - `AssembleStiffness()`
  - `AssembleSlipContributionIP()`
  - shared-face slip scatter path
- `miniapps/seas/integrator/dg_elasticity_ip_combined_integrator.hpp`

**Failure mode**

The BP5 production solve may still be using a slightly different effective jump enforcement than
the local cross-verification tests expose, especially across the full assembled system.

### 6.5 Medium-to-low probability: normal stress / friction feedback amplifies a small early mismatch

**Why it is plausible**

BP5 is highly sensitive. A small early traction or jump mismatch can be magnified by the
rate-state friction law until the divergence becomes visible at `0.5 yr`.

**Code areas**

- `miniapps/seas/fault/rate_state_fault.hpp`
- `miniapps/seas/friction/dieterich_ruina.hpp`
- `miniapps/seas/solver/seas_operator.hpp`

**Failure mode**

The root cause may be a small early elastic mismatch, but the observed “traction blow-up at
0.5 yr” is the friction system amplifying it rather than the original local elastic error.

### 6.6 Low probability: D8 slip sign convention

This was an important suspect earlier, but it is now lower probability because the strengthened
cross-verification test explicitly checks the prescribed-jump sign and passes for both
`sign_flipped` branches.

### 6.7 Low probability: exact-face identification, face area, penalty scaling

These are no longer serious suspects.

- face identification is matched
- area is matched
- penalty scaling is matched

### 6.8 Lowest probability: final traction formula itself

Line-by-line code inspection (Section 3.5) now **confirms** the traction formulas are
identical, not just structurally aligned:

- jump formula: `u⁻ - u⁺ - f_q` — exact match
- stress average: `0.5*(σ0+σ1)·n̂` with per-side material — exact match
- penalty correction: `-penalty * jump` — exact match
- trace inequality constant: `p*(p+D-1)/D` — exact match
- normal convention: outward from minus side — consistent
- area/volume in penalty: equivalent for affine tets — verified
- IP raw-vs-lifted form — equivalent (Section 3.6)

One latent mismatch found: `ComputePenalty` uses element-1 material for both sides
(Section 3.5.7). Harmless for BP5 (homogeneous material), confirmed by matching penalty
values in the first-solve dump.

This is now the **most thoroughly ruled out** location. The blow-up bug is upstream of the
local traction extraction kernel.

---

## 7. Current Working Hypothesis

The most likely explanation at this stage is:

1. MFEM carries a wrong fault jump / face state into the early accepted solves,
   probably before or during the local `[[u]] - slip` construction.
2. This error is not due to face identification, penalty scaling, or the D8 sign chain.
3. The quasi-dynamic BP5 friction update amplifies that early mismatch until it becomes
   visible as the traction blow-up around `0.5 yr`.

This is a narrower and more defensible statement than “the traction kernel is wrong.”

---

## 8. Immediate Next Debug Step

The next decisive step is no longer another timing refinement. It is code-level
investigation of the early fault-state / jump path, because the accepted-step MFEM dumps
already show that the tiny `jump_y` persists on both sides of Tandem’s `t=2.0e-02` dump.

Useful interpretation from the current data:

- accepted-step MFEM still has `jump_y ~ 1e-11 .. 1e-10`
- Tandem has `jump_y ~ 1e-08`
- therefore the bug is almost certainly upstream of the final exact-face traction readout

---

## 9. Retrospective Update After v58

This note originally treated the `~0.5 yr` failure as still unresolved. That is no longer the
current state.

### 9.1 What fixed the earlier `~0.5 yr` blow-up

The key later fix was commit `a32946316883ae42bc7115d3c0c71a499074b230`:

- `Fix double-loading of shared Dirichlet faces in boundary RHS`

The bug was in `AssembleDirichletLoading()` in
`miniapps/seas/domain/elasticity_operator.hpp`.

Before that fix, the exterior attr-`5` boundary loop processed all attr-`5` boundary elements,
including faces that were already being handled through the interior/shared Dirichlet skeleton
paths. As a result, shared Dirichlet faces could receive both:

- boundary-face RHS loading
- shared/interior skeleton RHS loading

Those two DG paths are not equivalent copies of the same operation. They use different formulas
and penalty roles, so this produced an incorrect effective RHS on shared Dirichlet faces.

Commit `a329463` fixed this by building exclusion sets from:

- `dirichlet_interior_faces_`
- `dirichlet_shared_faces_`

and skipping those face indices in the exterior boundary attr-`5` loop. The same exclusion logic
was also applied to the startup Dirichlet-face diagnostic so the reported counts matched the
actual assembled loading path.

In short:

- earlier `~0.5 yr` blow-up: strongly tied to double-loading of shared Dirichlet faces
- `a329463`: removed that double-loading and allowed the run to progress past `0.5 yr`

### 9.2 What changed after that

Later cleanup and hardening work replaced the temporary/gap-fix style BP5 face recovery with a
tag-only startup path:

- `65d2e9d` — `Use tag-only facet BC classification for BP5`
- `8cb7c00` — `Add accepted-step face tracing for BP5 diagnostics`

That later work:

- removed coordinate-based production classification for BP5
- removed the old `Y=0` gap-fix behavior
- recovered fault and Dirichlet faces from attrs `3` and `5` only
- propagated shared tagged faces by canonical face key
- added exact startup validation and shared-face agreement audit
- added accepted-step runtime tracing for the current blow-up investigation

### 9.3 Current interpretation

The original v58 analysis is still useful as a record of the early first-solve comparison, but it
should no longer be read as the latest statement of the BP5 failure mode.

Current status:

- the earlier `~0.5 yr` instability was successfully bypassed by fixing shared Dirichlet
  double-loading
- startup boundary/fault-face classification is now tag-only and validated
- the remaining active problem is a later `~2.5 yr` blow-up, which is more likely in the runtime
  traction / normal-stress / friction evolution than in startup face classification

---

## 10. References

### MFEM

- `miniapps/seas/solver/seas_operator.hpp`
- `miniapps/seas/domain/elasticity_operator.hpp`
- `miniapps/seas/integrator/dg_elasticity_ip_combined_integrator.hpp`
- `miniapps/seas/fault/rate_state_fault.hpp`
- `miniapps/seas/tests/unit/test_cross_verify_tandem.cpp`
- `miniapps/seas/tests/verification/bp5_verification_full.cpp`

### Tandem

- `app/form/SeasQDOperator.cpp`
- `app/form/AdapterOperator.h`
- `app/localoperator/ElasticityAdapter.cpp`
- `app/localoperator/AdapterBase.cpp`
- `app/localoperator/Elasticity.cpp`
- `app/kernels/elasticity.py`

### Relevant commits

- `932b797` — v57 shared fault face MPI detection fix
- `e3d4d31` — dump MFEM exact-face data at accepted PETSc step state
- `fef1e9f` — replay MFEM TQ at second accepted PETSc step
- `85525e2` — trigger MFEM TQ replay by accepted PETSc time crossing
- `a329463` — fix double-loading of shared Dirichlet faces in boundary RHS
- `65d2e9d` — use tag-only facet BC classification for BP5
- `8cb7c00` — add accepted-step face tracing for BP5 diagnostics
