# BP5 Debug v61: TSSolve Ineffective — Parallel Bug Hunt

**Date:** 2026-04-07
**Status:** All parallel shared-face suspects ruled out by diagnostic tests — awaiting serial (1-rank) test and partition-dependence results
**Scope:** BP5, bp5_tandem_exact.msh, IP, p=1, mesh-scale=1000, TSSolve

---

## 1. Critical Finding: TSSolve Has Zero Effect on Blowup

The v59 investigation concluded that the manual `TSStep` loop (breaking PETSc's
FSAL optimization) was the root cause of the ~25 yr blowup. Commit 51da03b
replaced the TSStep loop with a single `TSSolve` call matching Tandem's
architecture.

**Result:** The blowup still occurs at ~25 yr with TSSolve. The timing is
unchanged. The v59 root cause diagnosis was wrong.

This means the bug is in the **per-evaluation physics**, not the time stepping
framework. Something in each `Mult()` call produces a small but consistent
error that accumulates over ~25 yr regardless of whether PETSc manages the
steps via TSStep or TSSolve.

---

## 2. Observations from ParaView (TSSolve run, t ≈ 4.7 yr)

### 2.1 Displacement magnitude

- Max displacement: **3.8 m** at t = 1.47e8 s (4.68 yr)
- Expected: Vp*t = 1e-9 × 1.47e8 = **0.15 m**
- Ratio: **~25× too large**
- The displacement concentration has a **rectangular footprint matching the
  fault rectangle** (x ∈ [-50, 50] km, z ∈ [-40, 0] km on y=0)
- Secondary displacement concentrations below the fault (z < -40 km)
  at the fault-tip corners

### 2.2 Traction (dip component)

- **±5 MPa oscillatory pattern** at the boundary of the fault rectangle
- The oscillation follows the rectangular outline of the fault region:
  all four edges (top z=0, bottom z=-40 km, left x=-50 km, right x=+50 km)
- Interior of the fault and exterior Dirichlet region are smooth

### 2.3 Normal stress

- Normal stress goes **tensile** between the fault region (attr=3) and
  the Dirichlet continuation region (attr=5) on y=0
- The tensile zone follows the fault rectangle boundary

### 2.4 Displacement field details

- Displacement 0 (x-component): -1.7 to 17 m. Isolated hot spots (up to 17 m)
  at fault rectangle corners. Bulk fault region ~0-4 m.
- Displacement 1 (y-component): -1.0 to 1.2 m on the fault surface.
- The displacement is unphysical — no need to compare with Tandem. The Dirichlet
  continuation region (creeping at Vp) should have displacement O(Vp*t) = O(0.15 m).

---

## 3. What Has Been Verified Correct

### 3.1 Stiffness matrix K (from v59, still valid)

- Global K·1 norms match Tandem to **15 significant figures** (v59 Section 21)
- Per-element volume K matches to machine precision (v59 Section 18.2)
- Per-face skeleton K matches to machine precision (v59 Section 18.3)
- Volume integrator formula matches Tandem term-by-term (v59 Section 17)
- Face integrator formula matches Tandem term-by-term (v59 Section 17)
- CG vs MUMPS displacement matches to 2e-14 (v59 Section 15.4)
- No double-counting of boundary K on interior Dirichlet faces:
  `ExchangeFaceNbrData()` is called at line 2617 BEFORE `Assemble(0)` at
  line 2685, so `FaceIsTrueInterior` correctly returns true for shared
  faces, preventing `GetBdrFaceTransformations` from processing them

### 3.2 RHS vector b (from v59, first step only)

- Per-face RHS elvec matches Tandem to **7 significant figures** (v59 Section 13.2)
- Global b norms match at Stage 1 (pure initial conditions) with exact
  time/unit scaling (v59 Section 23.2)
- **Caveat:** The first-step comparison was at t=0.02s where Dirichlet
  loading ≈ Vp*0.004 = 4e-12 m (essentially zero). A bug proportional
  to the loading magnitude would be invisible at this time.

### 3.3 Coupling loop (from v59, first step only)

- Traction: 8-11 digit match per RK45 stage (v59 Section 25.2)
- Slip rate: 10-11 digit match per RK45 stage (v59 Section 25.3)
- Displacement: 6-7 digit match per RK45 stage (v59 Section 25.4)
- **Caveat:** Same as 3.2 — all values are near-zero at the first step.
  A bug that grows with slip or displacement would be invisible.

### 3.4 Mesh face classification (this investigation)

Y=0 Face Classification Audit output:
```
Fault (attr=3):      9257
Dirichlet (attr=5):  1377
NONE (interior):     0
NONE (shared):       0
OK: all y=0 interior faces classified.
```

No unclassified faces on y=0. All interior and shared faces have either
fault or Dirichlet BC. The mesh file has:
- 9,312 attr=3 triangles (fault on y=0)
- 1,388 attr=5 triangles (Dirichlet continuation on y=0)
- 240 attr=5 triangles (far-field true boundary)
- 2,014 attr=1 triangles (natural: top/bottom)
- Total: 12,954 boundary elements, 62,874 tetrahedra

The 66-face difference between mesh file counts (10,700) and audit counts
(10,634) reflects shared faces counted by one rank's `shared_face_bc_` but
whose centroid check ran on the other rank.

### 3.5 Friction parameter a

User verified via ParaView that the spatially varying friction parameter
`a(x2, x3)` is correctly assigned to all fault DOFs. The `a` field
matches the BP5 specification (a=0.004 in VW zone, a=0.04 in VS zone,
smooth transition in the 2 km ht zone).

### 3.6 Below-Wf hardcoded branch (dead code)

The `depths(i) > Wf_bp5_ + 1.0` check in `rate_state_fault.hpp` was
dead code: all fault DOFs have depth ≤ Wf = 40,000 m (the fault rectangle
stops at z = -40 km). The mesh has no fault faces below z = -40 km.
The removal of this branch (current diff on `test/revert-transform`)
has no effect on the simulation.

### 3.7 Coordinate systems

- Mesh file: coordinates in km (Tandem convention)
- MFEM: `--mesh-scale 1000` converts km → m
- Dirichlet loading threshold: |y| > 1000 m (= 1 km) — correct after scaling
- Fault depths: computed from mesh coordinates, stored as positive values
- All coordinate evaluations use reference (undisplaced) coordinates —
  correct for small-deformation linear elasticity

### 3.8 PETSc configuration

```
-ts_type rk
-ts_rk_type 5dp
-ts_rtol 1e-50
-ts_atol 1e-7
-ts_adapt_wnormtype infinity
-ts_adapt_dt_max 3.15576e6   (= 0.1 yr)
```

dt_max is NOT the binding constraint — adaptive dt is mostly below 0.1 yr.
Removing dt_max would not change behavior.

### 3.9 RHS reset

The `Solve()` function creates a fresh `LinFormType b` each call (line 4727),
assembles it to zero, then adds slip and Dirichlet contributions. No
accumulation of RHS across time steps.

---

## 4. Code Fix Applied (Non-causal)

### 4.1 Face->Transform → Elem1->Transform in fault coordinate getters

`GetFaultDepths()` and `GetFaultCoords2D()` in `elasticity_operator.hpp`
were using `FTr->Face->Transform(nip, coords)` — the same unreliable
path identified in v59 Section 6 for Dirichlet loading.

**Fix:** Replaced with the robust element-transform pattern:
```cpp
FTr->SetAllIntPoints(&nip);
FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords);
```

**Impact on blowup:** None. User verified `a` values are correct via
ParaView, meaning `Face->Transform` coordinates were close enough for
correct friction parameter assignment. This fix is a correctness
improvement but does not address the blowup.

---

## 5. Remaining Suspects

### 5.1 Shared-face RHS assembly (HIGHEST priority)

The blowup originates at the fault/Dirichlet boundary on y=0 — exactly
where faces are most likely to become shared faces after MPI partitioning.

The v59 first-step comparison verified b at t=0.02s where loading ≈ 0.
A bug in the shared Dirichlet RHS that scales with Vp*t would be
invisible at the first step but grow linearly to cause the blowup.

Specific suspects in `AssembleDirichletLoading()` shared path (lines 4217-4279):
- `ComputeSkeletonDirichletSign(FTr)` for shared faces: the Jacobian
  from `GetSharedFaceTransformations` might differ from the interior
  version, producing wrong `dir_sign` on some faces
- elem1/elem2 assignment: for shared faces, only elem1 is local.
  The RHS assembles only ev1. If elem1/elem2 are swapped relative to
  the interior case, the sign and magnitude of ev1 would be wrong.
- Missing or double contributions: if some shared Dirichlet faces
  are processed by both ranks instead of one, or skipped by both.

### 5.2 Shared-face slip RHS assembly (HIGH priority)

Same logic as 5.1 but for `AssembleSlipContributionIPShared()`. A
bug in the shared fault slip RHS at the fault/Dirichlet junction
faces could produce wrong displacement at the boundary.

### 5.3 Some other parallel-specific path (MEDIUM priority)

- `ExpandOwnedToLocalFault` / `RestrictToOwnedFault` — ghost DOF
  communication for slip and traction. If wrong values are sent to
  ghost DOFs on shared faces, the per-face RHS would be wrong.
- `canonical_to_local_perm_` — DOF permutation between canonical
  (sorted global vertex ID) and MFEM-local ordering. If wrong on
  some faces, slip and traction DOFs would be misaligned.

### 5.4 Local kernel bug masked by first-step comparison (LOW priority)

A bug in the per-face RHS formula that is proportional to the input
value (slip or u_D). At the first step, inputs are ~0, so the error
is invisible. At finite time, the error grows with the input.

Example: a wrong coefficient (0.5 instead of 1.0) in one of the
DG terms. The first-step comparison would show the ratio is correct
(because 0 * 0.5 = 0 * 1.0 = 0), but at finite time the error
would be proportional to the loading.

---

## 6. New Evidence: Displacement Hot Spots Correlate with MPI Partition Boundaries

### 6.1 Observation

ParaView visualization of displacement magnitude overlaid with mesh wireframe
shows displacement hot spots (up to 16 m) at specific isolated locations:

- Below the fault rectangle (z < -40 km) at ~2 locations
- At the far-left domain corner (x ≈ -200 km)

A separate ParaView visualization of the `mpi_rank` field shows partition
boundaries (rank color transitions) at approximately the same spatial locations.

### 6.2 Significance

A displacement of 16 m at x ≈ -200 km is **~200× too large** (expected
Vp*t/2 ≈ 0.07 m). This is in the far-field, far from the fault, where
the displacement should be smooth and small. The only special thing about
these locations is that they sit on **MPI partition boundaries**.

This is strong evidence that the blowup is caused by a **parallel
shared-face assembly bug** — not a formula error in the local kernels.

### 6.3 Mechanism

When interior Dirichlet faces at y=0 (attr=5) are split across MPI
partitions, they become shared faces. The RHS assembly for these faces
goes through `AssembleDirichletLoading` shared path (lines 4217-4279)
and `AssembleSlipContributionIPShared`. If these paths have a sign,
coefficient, or double-contribution error, the wrong RHS on shared
faces would:

1. Produce wrong displacement at elements adjacent to the shared face
2. The wrong displacement would propagate through the traction computation
3. The wrong traction would feed back into the friction law
4. The error would accumulate over time proportional to Vp*t

The localization to MPI boundaries explains why:
- The first-step comparison (v59) did not catch it: loading ≈ 0 at t=0.02s
- The per-element K comparison matched: K is the same on all faces
- The global K·1 norm matched: K assembly uses ExchangeFaceNbrData correctly
- The blowup is independent of time stepping method: the per-step error
  is in the RHS assembly, not the ODE integration

---

## 7. INVALIDATED: Ghost DOF Permutation Hypothesis

The initial v61 analysis hypothesized that `ExpandOwnedToLocalFault`
sends DOFs in MFEM-local order while the receiver expects canonical
order. **This was wrong.** Tracing the data flow:

1. `owned_fault_dof_to_local_dof_` (line 2459-2464) is built with the
   permutation baked in: position `kk` in owned_data maps to canonical
   DOF `kk`.
2. `RestrictToOwnedFault` uses this mapping, so `owned_data` is already
   in canonical order.
3. The old send code `buffer[kk] = owned_data[owned_face*nbf + kk]`
   sends canonical DOF `kk` at buffer position `kk` — **correct**.
4. The receive code reads `buffer[canonical_kk]` and applies
   `canonical_to_local_perm_` — **also correct**.

The proposed send-path "fix" was reverted as it would introduce a
double-permutation bug.

---

## 8. Things Confirmed NOT the Bug

### 8.1 Ghost DOF permutation in ExpandOwnedToLocalFault — NOT the bug

`owned_data` is in canonical order (permutation applied during
`RestrictToOwnedFault` via `owned_fault_dof_to_local_dof_`). The send
path sends in canonical order. The receive path applies
`canonical_to_local_perm_` to place at MFEM-local positions. Both
directions are correct.

### 8.2 K symmetry term sign error on shared faces — NOT the bug

Analysis of `AssembleFaceBlock` shows the off-diagonal blocks use
different c1 coefficients that compensate for the normal reversal:

- Block (1,0): `c1 = 0.5*ε = -0.5` (SIPG)
- Block (0,1): `c1 = -0.5*ε = +0.5` (SIPG)

When Rank B processes its (0,1) block with reversed normal (-n):
- Symmetry = `+0.5 * shape_A * T_B(-n)` = `+0.5 * (-T_B(+n))` = `-0.5 * T_B(+n)`
- Matches serial (1,0): `-0.5 * T_B(+n)` ✓

The sign flip in c1 exactly compensates the sign flip in the traction
operator. The explicit symmetrization (lines 173-181) then enforces
exact machine-precision symmetry. K·1 norms matching to 15 digits
confirms this.

### 8.3 Double contribution on shared faces — NOT the bug

`fault_interior_faces_` and `fault_shared_faces_` are disjoint
(populated from separate branches in `BuildFacetBCTables`, lines
1553-1571). Same for `dirichlet_interior_faces_` / `dirichlet_shared_faces_`.
No face appears in both lists. Each shared face is processed once per
rank, contributing only ev1 (local element).

### 8.4 DOF permutation in BuildSlipAtQuadPoints — NOT the bug

`local_slip_` is in MFEM-local order (permutation correctly applied
during `ExpandOwnedToLocalFault`). `InterpolateToQuadPoints` uses
`e_q_(k,q)` — reference shape functions at quadrature points. MFEM-local
DOF k corresponds to the k-th reference vertex on each rank's face
parametrization, so the interpolation is self-consistent.

### 8.5 Double K on interior Dirichlet faces — NOT the bug

`ExchangeFaceNbrData()` is called at line 2617 BEFORE `Assemble(0)` at
line 2685. So `FaceIsTrueInterior` correctly returns true for shared
faces (Elem2Inf >= 0), and `GetBdrFaceTransformations` returns null,
preventing boundary K on interior/shared faces.

### 8.6 Below-Wf hardcoded branch removal — NOT the bug (dead code)

All fault DOFs have depth ≤ Wf = 40,000 m. The mesh has no fault faces
below z = -40 km. The `depths(i) > Wf_bp5_ + 1.0` check never triggered.
Removal has no effect.

### 8.7 Friction parameter a — NOT the bug

User verified via ParaView that spatially varying `a(x2, x3)` is correct
on all fault DOFs.

### 8.8 Face normal MUST flip on shared faces (structural proof from K)

**Proof by contradiction from K correctness (v59 Section 21).**

The DG bilinear form K_11 (elem1 test × elem1 trial) and K_22 (elem2
test × elem2 trial) have opposite-sign consistency terms:

```
K_11 = -ε/2 · ∫ σ(u1)·n · v1  -  ε · ε/2 · ∫ σ(v1)·n · u1  +  κ · ∫ u1·v1/h
K_22 = +ε/2 · ∫ σ(u2)·n · v2  +  ε · ε/2 · ∫ σ(v2)·n · u2  +  κ · ∫ u2·v2/h
        ^^^^                        ^^^^
       opposite sign               opposite sign
```

For MFEM's `ParBilinearForm::AssembleSharedFaces`, rank B's elem1 is
physically the same element as rank A's elem2. Rank B computes K_11'
using its face normal `nor_B`. For the global K to be correct:

- If `nor_B = nor_A`: K_11'_B has `-ε/2` sign → but K_22 needs `+ε/2` → **WRONG**
- If `nor_B = -nor_A`: K_11'_B has `-ε/2 · (-1) = +ε/2` → matches K_22 → **CORRECT**

Since K matches Tandem to **15 significant figures** (v59 Section 21),
the face normal **must flip**: `nor_B = -nor_A` on shared faces. This is
a structural invariant enforced by MFEM's `GetSharedFaceTransformations`.

All subsequent shared-face sign analyses use this proven invariant.

### 8.9 Shared-face slip RHS (AssembleSlipContributionIPShared) — NOT the bug

**Given** `nor_B = -nor_A` (Section 8.8), the Tandem fault basis sign
convention ensures correctness.

`ComputeOrientedFrame` (fault_basis.hpp:366-438) orients the raw face
normal to `ref_normal = (0,-1,0)` and then, if the raw normal was
flipped, **negates ALL output vectors** (normal, tangent1, tangent2) at
line 430-438. Trace for y=0 face:

| Rank | Raw normal   | dot < 0? | sign_flipped | Final normal | tangent1 (dip)| tangent2 (strike)|
|------|-------------|----------|-------------|-------------|---------------|-----------------|
| A    | (0,+H,0)   | Yes      | true        | (0,+1,0)    | (0,0,+1)      | (-1,0,0)        |
| B    | (0,-H,0)   | No       | false       | (0,-1,0)    | (0,0,-1)      | (1,0,0)         |

All basis vectors are negated: `tangent_A = -tangent_B`. Therefore the
embedded 3D slip also flips: `f_q_A = -f_q_B`.

In `AssembleSlipFaceRHS` (line 267-313), with `nor_B = -nor_A` and
`f_q_B = -f_q_A`:

```
trac_test_B = λ·Dx2·(f_B·nor_B) + μ·(Dx2·nor_B · f_B + Dx2·f_B · nor_B)
            = λ·Dx2·((-f_A)·(-n_A)) + μ·(Dx2·(-n_A)·(-f_A) + Dx2·(-f_A)·(-n_A))
            = λ·Dx2·(f_A·n_A) + μ·(Dx2·n_A·f_A + Dx2·f_A·n_A)
            = trac_test_2_A  ← same as rank A's elem2 traction test

penalty_B = +penalty · shape2 · (-f_A) · nl = -penalty · shape2 · f_A · nl
          = penalty_2_A     ← matches rank A's elem2 penalty (-penalty)
```

**Result:** `elvec1_B = elvec2_A`. Each physical element gets its correct
contribution. ✓

### 8.10 Shared-face Dirichlet RHS (AssembleDirichletLoading shared) — NOT the bug

**Given** `nor_B = -nor_A`, `ComputeSkeletonDirichletSign(FTr)` (line
1300-1311) returns opposite signs on the two ranks:

```
dir_sign_A = sign(dot(nor_A, ref_normal))
dir_sign_B = sign(dot(-nor_A, ref_normal)) = -dir_sign_A
```

So `u_D_B = -u_D_A`. This is structurally identical to the fault slip
case (f_q flips). The same traction-test + penalty analysis applies:

```
elvec1_B(nor_B=-n_A, u_D_B=-u_D_A) = elvec2_A(nor_A, u_D_A) ✓
```

The Dirichlet loading for both the IP path (line 4255-4279) and the BR2
path (line 4281-4435) use this same sign mechanism.

### 8.11 Shared-face traction computation (ComputeTractionImpl) — NOT the bug

The traction formula at each quad point (dg_elasticity_ip_combined_
integrator.hpp:500-539) is:

```
T = 0.5·(σ1 + σ2)·n_hat + (-penalty)·(u1 - u2 - f_q)
```

With `nor_B = -nor_A`, `u1_B = u2_A`, `u2_B = u1_A`, `f_q_B = -f_q_A`:

- Stress: `0.5·(σ2_A + σ1_A)·(-n_hat_A) = -T_stress_A`
- Penalty: `(-penalty)·(u2_A - u1_A - (-f_q_A)) = +penalty·(u1_A - u2_A - f_q_A) = -penalty_corr_A`
- Total: `T_B = -(T_stress_A + penalty_corr_A) = -T_A`

The 3D traction flips: `T_B = -T_A`. But the fault basis vectors also
flip (Section 8.9), so the projected tangential traction is the same:

```
tau_dip_B  = T_B · tangent1_B = (-T_A)·(-tangent1_A) = T_A·tangent1_A = tau_dip_A  ✓
tau_strike_B = T_B · tangent2_B = (-T_A)·(-tangent2_A) = T_A·tangent2_A = tau_strike_A ✓
sigma_n_B    = T_B · normal_B   = (-T_A)·(-normal_A)   = T_A·normal_A   = sigma_n_A   ✓
```

For owned shared faces, only the owner's traction is kept (via
`RestrictToOwnedFault`), and the owner's projected traction is correct.

### 8.12 Face-neighbor displacement exchange — NOT the bug

`ComputeTractionImpl` calls `pfes->ExchangeFaceNbrData()` and
`par_u.ExchangeFaceNbrData()` at lines 5517-5518, immediately before the
shared fault face traction loop. The face-neighbor displacement `u2` is
current for each solve.

### 8.13 No double-counting between boundary and skeleton Dirichlet — NOT the bug

`AssembleDirichletLoading` boundary loop (line 3674) builds
`dir_interior_set` and `dir_shared_set` (lines 3658-3672) and skips any
face in these sets (lines 3684-3685). Interior Dirichlet faces go through
the skeleton path (line 3896). Shared Dirichlet faces go through the
shared path (line 4217). No face is processed twice.

---

## 9. RESOLVED: Shared-Face RHS Sign Error Hypothesis

### 9.1 Original hypothesis

If the face normal does NOT flip between ranks, `AssembleSlipFaceRHS`
would produce `elvec1_B` with the wrong penalty sign, since elem1
always gets `+penalty` but rank B's elem1 is physically elem2 (which
needs `-penalty`).

### 9.2 Resolution

**The normal DOES flip** — proved structurally from K correctness
(Section 8.8). Given the flip:

- Fault slip: `f_q` flips via Tandem basis convention → `elvec1_B = elvec2_A` (Section 8.9)
- Dirichlet: `u_D` flips via `ComputeSkeletonDirichletSign` → `elvec1_B = elvec2_A` (Section 8.10)
- Traction: `T_B = -T_A`, projected traction same due to basis negation (Section 8.11)

**All shared-face sign conventions are correct.** The hypothesis is ruled out.

### 9.3 What this means

The MPI-boundary-correlated hot spots (Section 6) are NOT caused by a
sign error in the shared-face RHS or traction assembly. The root cause
must be a more subtle mechanism — possibly:

- A bug in a code path not audited (e.g., MFEM internal face orientation
  edge case, specific face geometries)
- A local kernel bug that merely correlates with MPI boundaries due to
  mesh density patterns near partition edges
- A numerical issue (penalty scaling, element quality) amplified at
  partition boundaries

The **serial (1-rank) test** (Section 11) remains the decisive experiment.

---

## 10. Partition-Dependence Test (submitted)

Submitted runs at **200, 400, and 600 MPI ranks** on the same mesh.

**If hot spots move with partitioning:** Definitive proof the bug is
MPI-partition-dependent. The hot spots follow the shared faces, not the
physics. This confirms the shared-face RHS assembly suspect (Section 9).

**If hot spots stay at the same locations:** The bug is in a local kernel
and merely correlates with mesh features near partition boundaries (e.g.,
mesh refinement transition at fault edges). This would shift suspicion to
Section 5.4 (local kernel bug proportional to input).

---

## 11. Proposed Additional Experiment

### 11.1 Serial (1-rank) test on the same mesh

Run the production mesh (`bp5_tandem_exact.msh`) on **1 MPI rank** with
MUMPS direct solver. The mesh has 62,874 tets → ~750K DOFs, which
should fit on a single Frontera node (96 GB).

**If the blowup disappears at 1 rank:** the bug is in the parallel
shared-face path (Section 5.1-5.3). This would narrow the investigation
to `AssembleDirichletLoading` shared path, `AssembleSlipContributionIPShared`,
or the ghost DOF communication.

**If the blowup persists at 1 rank:** the bug is in the local kernel
(Section 5.4). There are no shared faces at 1 rank, so the investigation
shifts to the interior-face RHS assembly, traction computation, or
friction coupling.

### 6.2 Why this is decisive

At 1 rank:
- No shared faces → no shared-face K, no shared-face RHS
- No ghost DOF communication
- No `ExpandOwnedToLocalFault` / `RestrictToOwnedFault` communication
- No `canonical_to_local_perm_` permutation issues

This eliminates ALL parallel-specific suspects in a single test.

### 6.3 Job configuration

```bash
ibrun -n 1 ./seas_bp5_full \
      --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
      --mesh-scale 1000 \
      --order 1 \
      --dg-method IP \
      --tandem-time-stepping \
      --solver mumps \
      --tfinal 1e9 \
      --ref-dir bp5/benchmark_data \
      --output-dir "${RESULT_DIR}" \
      --output-prefix "${OUTPUT_PREFIX}" \
      --petsc-ts \
      --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
      --paraview-dt 1577880
```

Note: K factorization will be slow (~minutes to an hour for 750K DOFs
with MUMPS), but only happens once. The simulation should reach a few
years of simulated time within a 2-hour development job. Check
displacement at t ≈ 5 yr — if it's O(0.15 m), the serial case is
correct and the bug is parallel-specific.

---

## 12. Parallel Diagnostic Tests (v61)

Three targeted tests added to `test_parallel_elasticity.cpp` to rule out
the six shared-face suspects. All tests use forced y-partitioning on a
`CreateTestMesh3DTet(2,1,1, ...)` mesh to guarantee ALL fault faces
become shared (zero interior fault faces).

### 12.1 Test A: Coordinate consistency (Face→Transform vs Elem1→Transform)

**Test:** At each DOF point on every shared fault face, compare physical
coordinates from `Face->Transform(ip)` vs `Elem1->Transform(GetElement1IntPoint())`.

**Result:** PASSED (`max_diff=8.88e-16`, 16 shared faces, 0 QP disagreements).
Face and Elem1 transforms agree to machine precision at shared faces.

**Rules out:** Suspects #3 (wrong coordinate transform at shared faces)
and #5 (Dirichlet loading mis-coords).

### 12.2 Test B: Fault basis sign-baking consistency

**Test:** Both ranks sharing a fault face compute CalcOrtho independently.
Verify that basis vectors differ by exactly a global sign factor (±1):
`basis_A = ±basis_B` with all 3 vectors (normal, tangent1, tangent2)
sharing the same sign.

**Result:** PASSED (`max_sign_err=0.0`, 8 pairs, all flipped, 0 inconsistent).
Sign-baking is consistent: all shared faces have `basis_A = -basis_B`,
correctly compensated by `nor_B = -nor_A`.

**Rules out:** Suspect #1 (inconsistent basis vectors at shared faces).
The sign-baking is a designed behavior (Tandem convention), not a bug.

### 12.3 Test C: Ghost DOF expand/restrict roundtrip

**Test:** Set owned fault DOFs to `f(x2, x3) = 7·sin(x2) + 3·cos(x3)`
(a smooth function of coordinates). ExpandOwnedToLocalFault distributes
to ghost DOFs via MPI. Verify ALL local DOFs (including ghosts) match
expected values from independently-computed local coordinates.

**Result:** PASSED (`max_err=5.55e-15`, `rt_err=0.0`, 0 mismatches, 8 owned + 16 total).
Ghost communication, face matching, and canonical→local permutation are
all correct.

**Rules out:** Suspects #6 (ghost comm), #6b (face matching in
`shared_fault_comm_blocks_`), and canonical→local permutation errors.

### 12.4 Summary

| Suspect | Test | Result | Verdict |
|---------|------|--------|---------|
| #1 Basis vector disagreement | Test B | sign-baking consistent | Ruled out (by design) |
| #2 CalcOrtho normal flip | Test B | flips correctly | Ruled out |
| #3 Wrong coordinate transform | Test A | 8.88e-16 agreement | Ruled out |
| #4 Dirichlet sign computation | Section 8.10 | structurally correct | Ruled out |
| #5 Dirichlet loading mis-coords | Test A | 8.88e-16 agreement | Ruled out |
| #6 Ghost DOF communication | Test C | 5.55e-15 roundtrip | Ruled out |
| #6b Face matching | Test C | 0 mismatches | Ruled out |

All six shared-face suspects have been ruled out on the small test mesh.
The existing `test_serial_parallel_displacement_match` also passes with
`rel_err=7.38e-13` (serial and parallel displacement identical).

**Conclusion:** The blowup is NOT caused by any of the custom shared-face
assembly paths. The bug is either in a code path not yet tested, or
manifests only on the production mesh / after extended time stepping.

### 12.5 Code changes in this section

- `fault_basis.hpp` `ComputeOrientedFrame`: Added documentation clarifying
  the sign-baking convention and its correctness for shared faces
- `elasticity_operator.hpp`: Wrapped `GetComm()` debug code in
  `if constexpr (IsParallelMesh)` to fix compilation with serial Mesh
- `test_parallel_elasticity.cpp`: Added 3 new diagnostic tests

---

## 13. Summary of v59 Conclusions That Remain Valid

| v59 Finding | Still valid? | Notes |
|-------------|-------------|-------|
| K matches Tandem (global, per-element, per-face) | Yes | Independent of time stepping |
| b matches at first step | Yes | But first step has ~0 loading |
| Coupling matches stage-by-stage | Yes | But only at ~0 loading |
| Volume integrator formula matches | Yes | Independent of time stepping |
| Face integrator formula matches | Yes | Independent of time stepping |
| TSSolve is root cause | **NO** | Blowup unchanged with TSSolve |
| Per-element K matrices identical | Yes | Independent of time stepping |
| Solver (MUMPS) is correct | Yes | CG matches MUMPS |

---

## 13. Files and Commits

### Modified in this investigation

- `miniapps/seas/domain/elasticity_operator.hpp`
  - `GetFaultDepths()`: Face->Transform → Elem1->Transform
  - `GetFaultCoords2D()`: Face->Transform → Elem1->Transform

### Run configuration

- Mesh: `bp5/mesh/reference/bp5_tandem_exact.msh`
- Mesh scale: 1000 (km → m)
- Job script: `jobs/bp5/bp5_v60_paraview_dev_2hr.sbatch`
- Flags: `--tandem-time-stepping --petsc-ts`
- PETSc config: `tests/verification/petsc_ts_rk45_tandem.cfg`
- 8 nodes, 400 MPI ranks on Frontera
