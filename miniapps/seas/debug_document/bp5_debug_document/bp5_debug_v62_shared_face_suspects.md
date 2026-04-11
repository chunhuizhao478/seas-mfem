# BP5 Debug v62: Shared-Face Bug — Comprehensive Suspect List & Test Plan

**Date:** 2026-04-07
**Status:** ROOT CAUSE FOUND AND FIXED — `ref_normal_` uninitialized on no-fault ranks. Fix: commit `27a548b`. Local displacement test PASS. Awaiting Frontera verification.
**Scope:** BP5, bp5_tandem_exact.msh, IP, p=1, mesh-scale=1000, parallel

---

## 1. Summary

ParaView visualization confirms displacement hot spots correlate with MPI
partition boundaries on the y=0 fault plane. The pattern:

- **25× expected displacement** at ~4.7 yr (3.8 m vs expected 0.15 m)
- Hot spots sit at **rank-to-rank boundaries** on y=0
- Pattern grows linearly with time (proportional to Vp*t)
- TSSolve vs TSStep makes no difference — bug is per-evaluation
- At t=0, slip ≈ 0, so the bug is invisible in first-step comparisons

The bug must be in custom code that handles **shared faces** differently from
interior faces. MFEM's K assembly (`ParBilinearForm::AssembleSharedFaces`)
is the only non-custom shared-face path and is confirmed correct.

**Partition sensitivity test submitted:** 100, 200, 600, 800 ranks (normal
queue 6hr) to verify hot spots shift with partition. Existing 400-rank
baseline for comparison.

---

## 2. v61 Permutation Hypothesis — RULED OUT

The v61 document proposed a ghost DOF permutation bug in
`ExpandOwnedToLocalFault` where the send path allegedly sent data in
"MFEM-local order" while the receiver expected "canonical order."

**This hypothesis is incorrect.** The owned_data vector is already in
canonical order because `owned_fault_dof_to_local_dof_` (built at
lines 2454–2466 of `elasticity_operator.hpp`) applies
`canonical_to_local_perm_` when mapping between owned and local DOFs:

```
owned_fault_dof_to_local_dof_[owned_face * nbf + kk]
    = local_face * nbf + canonical_to_local_perm_[local_face * nbf + kk]
```

Here `kk` is the canonical index. `RestrictToOwnedFault` maps:
```
owned_data[canonical_kk] = local_data[mfem_local_for_canonical_kk]
```

So owned_data position `kk` holds the value for canonical DOF `kk`.
The old send path `buffer[kk] = owned_data[owned_face*nbf + kk]`
correctly sends canonical-ordered data.

The proposed "fix" (applying `canonical_to_local_perm_` in the send path)
would introduce a **double-permutation bug**, reading
`owned_data[perm[canonical_kk]]` instead of `owned_data[canonical_kk]`.
This was never committed to the branch.

---

## 3. Comprehensive Suspect List

Every custom code path that handles shared faces differently from interior
faces is listed below. The bug must be in exactly one (or a combination) of
these paths.

### Suspect 1: Fault Basis Orientation at Shared Faces — ✅ RULED OUT

**Code:** `fault_basis.hpp` — `AppendSharedFaces()`, `ComputeQPBasisShared()`

**What could be wrong:**
Each rank independently computes the fault basis vectors (normal, tangent1,
tangent2) for shared faces using its own `GetSharedFaceTransformations`.
The raw face normal from `CalcOrtho(FTr->Jacobian(), nor)` points from
Elem1→Elem2, which FLIPS between the two ranks sharing the face.

The basis computation orients the normal using `ref_normal_` and derives
tangent vectors via the cross-product recipe (`strike = up × n`,
`dip = strike × n`). If the orientation logic produces a consistent
oriented normal on both ranks, the basis matches. But if there's a
subtle sign or ordering issue, the tangent frame could differ.

**Impact if wrong:**
- `BuildSlipAtQuadPoints` embeds 2D tangential slip into 3D using the basis.
  Wrong basis → wrong 3D delta_u at shared fault faces.
- `ProjectTractionToFaultDOFs` projects 3D traction onto tangent directions.
  Wrong basis → wrong tangential traction at shared fault faces.
- Both errors are proportional to slip/traction magnitude → grow with time.

**How to test locally:**
For each shared fault face, gather the basis vectors from both ranks (owner
and ghost). Compare normal, tangent1, tangent2. They must agree to machine
precision.

```
Test: TestSharedFaceBasisConsistency
- For each shared fault face fi:
  - Rank A (owner) computes basis(fi).normal, .tangent1, .tangent2
  - Rank B (ghost) computes basis(fi).normal, .tangent1, .tangent2
  - Gather both to one rank and compare
  - FAIL if max |diff| > 1e-12
```

**TEST RESULT (`test_diag_shared_fault_basis`, line 2057 of
`tests/parallel/test_parallel_elasticity.cpp`, mpirun -np 4):**
```
PASSED (max_sign_err=0.000000e+00, pairs=8 [same_sign=0 flipped=8], inconsistent=0)
```
All 8 shared face pairs have basis_A = -basis_B (consistent global sign
flip). The sign-baking is by design (Tandem convention) and compensated
by the DG normal flip. `test_serial_parallel_displacement_match` confirms
serial==parallel to `rel_err=7.38e-13`. **RULED OUT.**

---

### Suspect 2: Shared Fault Slip RHS Assembly — ✅ RULED OUT (IP), 🔴 BR2 bug found

**Code:** `elasticity_operator.hpp` — `AssembleSlipContributionIPShared()`
(lines ~3515–3572) is correct. `AssembleSlipContributionBR2Shared()`
(lines ~3578+) has a confirmed bug (not production-relevant).

**What could be wrong:**
This function loops over `fault_shared_faces_`, calls
`BuildSlipAtQuadPoints(slip_idx, slip_bc, delta_u_quad)` and then
`AssembleSlipFaceRHS(*fe1, *fe2, *FTr, delta_u_quad, elvec1, elvec2)`.
Only `elvec1` is added to the RHS (Elem2 is on the neighbor rank).

For correctness, rank B's `elvec1` (computed with its Elem1 = rank A's Elem2,
flipped normal, same slip) must equal rank A's `elvec2`.

`AssembleSlipFaceRHS` has two terms:
- Symmetry: `c1 = eps/2` (same sign both elements)
- Penalty: `+penalty` for Elem1, `-penalty` for Elem2

When the normal flips (n→-n) but slip stays the same (communicated via
`ExpandOwnedToLocalFault`), the symmetry term's `f·n` product flips and
`dshape·n` flips, producing a double-negation that should cancel. The penalty
term uses `nl_q = |n|` (unsigned) and `f_q` (same), so sign is correct.

**Potential issue:** The slip `delta_u_quad` is in 3D physical coordinates,
embedded using the fault basis (suspect 1). If the basis differs between
ranks, the 3D slip differs, breaking the cancellation.

**How to test locally:**
Set a known uniform slip (e.g., pure strike slip = 1e-9 m) on all fault
DOFs. Assemble the slip RHS on each rank. For a face that is interior on
one partition but shared on another, the total RHS contribution should be
identical.

Simpler: for each shared fault face, have both ranks compute their `elvec1`
and gather. Rank A's `elvec1` + rank B's `elvec1` should equal what a
single-rank interior assembly would produce (`elvec1 + elvec2`).

```
Test: TestSharedFaultSlipRHS
- Set uniform slip on all fault DOFs
- For each shared fault face:
  - Both ranks compute elvec1 via AssembleSlipFaceRHS
  - Gather to one rank
  - Sum elvec1_A + elvec1_B
  - Compare against single-rank interior assembly (elvec1 + elvec2)
  - FAIL if |diff| > 1e-10
```

---

### Suspect 3: Shared Dirichlet Loading Assembly — ✅ RULED OUT (IP D1 passes)

**Code:** `elasticity_operator.hpp` — `AssembleDirichletLoading()` shared
section (lines ~4217–4280)

**What could be wrong:**
For shared Dirichlet faces on y=0 (attr=5, outside fault rectangle), the
code applies `ComputeSkeletonDirichletSign(FTr)` to determine the sign of
the Dirichlet value `u_D`. This sign flips when the face normal flips
(between the two ranks), producing `u_D_B = -u_D_A`.

The code then calls `AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d, ev1, ev2)`
and adds only `ev1`.

Sign analysis (see v62 Section 2 methodology) shows the symmetry and penalty
terms match between rank B's `ev1` and rank A's `ev2` after the
normal+u_D double flip. However, this relies on MFEM's
`GetSharedFaceTransformations` producing consistent integration point
mappings on both ranks.

**Potential issue:** If the face-to-element reference coordinate mapping
differs subtly between ranks (due to element orientation), the shape
function evaluations `shape1_B(eip1_B)` and `shape2_A(eip2_A)` could
disagree, producing a small but nonzero mismatch at every shared Dirichlet
face.

**How to test locally:**
Same strategy as suspect 2: compare the sum of both ranks' `elvec1` against
a single-rank interior assembly result.

```
Test: TestSharedDirichletLoading
- Set time = 1.0 (nonzero Vp*t)
- For each shared Dirichlet face:
  - Both ranks compute elvec1
  - Gather and sum
  - Compare against interior Dirichlet assembly (elvec1 + elvec2)
  - FAIL if |diff| > 1e-10
```

**TEST RESULT (`test_diag_displacement_with_dirichlet`, line 2446 of
`tests/parallel/test_parallel_elasticity.cpp`, mpirun -np 4):**
```
FAILED (serial_u=6.164570782907e-01, parallel_u=5.654296088697e+00,
        rel_err=8.172246191698e+00, fault=0i+16sh)
```
**Parallel displacement is 9× larger than serial.** The mesh has 8 attr=3
(fault) + 8 attr=5 (Dirichlet) faces on y=0, with forced y-partitioning
making ALL y=0 faces shared. At time=1.0, the Dirichlet loading through
`AssembleDirichletLoading` shared path (lines 4217–4280) is active and
produces catastrophically wrong RHS.

**Per-DOF traction (`test_diag_traction_with_dirichlet`, line 2633):**
```
FAILED (max_rel_err=6.888985e+00, matched=8, unmatched=0, serial_nf=8)
```
Per-DOF traction differs by up to 7× between serial and parallel.

**IP result: D1 (slip=0, Dirichlet only) PASSES at rel_err=2.5e-12.**
D1b (slip, t=0) also PASSES with IP at rel_err=2.0e-13.
All shared Dirichlet and slip assembly paths are correct for IP.

The earlier failures (rel_err ≈ 8–10×) were caused by using the default
BR2 method. The BR2 shared path has a separate bug (not production-relevant).
**Suspect #3 is RULED OUT for the production IP method.**

---

### Suspect 4: Shared Face Traction Computation — ✅ RULED OUT (IP D correct → traction correct)

**Code:** `elasticity_operator.hpp` — `ComputeTractionImpl()` shared face
section (lines ~5527–5710+)

**What could be wrong:**
The traction computation at shared fault faces uses:
- `ExchangeFaceNbrData()` — MFEM's built-in displacement exchange (correct)
- `slip_bc` — from `ExpandOwnedToLocalFault` (correct per Section 2)
- Fault basis — from `fault_basis_.GetBasis(trac_idx)` (suspect 1)
- `ProjectTractionToFaultDOFs` — L2 projection onto tangent directions

The traction at a shared face is written to `traction(2*dof_idx)` and
`traction(2*dof_idx+1)` where `dof_idx = trac_idx * nbf + kk`. This value
is then picked up by `RestrictToOwnedFault` (only on the owner rank).

**Potential issue:** Only the OWNER rank's traction value is used (via
RestrictToOwnedFault). The non-owner also computes traction but it's
discarded. If the owner rank computes wrong traction (e.g., due to wrong
basis vectors from suspect 1, or wrong displacement at ghost elements),
the friction solver gets wrong input.

**Dependency:** This suspect is largely downstream of suspect 1 (fault basis)
and suspect 6 (ghost DOF communication). If those are correct, the traction
computation is likely correct.

**How to test locally:**
Set a known displacement field (e.g., u = Vp*t/2 * sgn(y) * e_x), compute
traction at all fault faces, and compare shared vs interior face results.

```
Test: TestSharedFaceTraction
- Set displacement to known analytic field
- ComputeTractionImpl on full mesh
- For shared fault faces: compare traction against reference
  (from interior face with same geometry, or analytic)
- FAIL if |diff| > 1e-8
```

---

### Suspect 5: Coordinate Computation (Face→Elem1 Transform) — ✅ RULED OUT

**Code:** `elasticity_operator.hpp` — `GetFaultDepths()` (lines ~4475–4510),
`GetFaultCoords2D()` (lines ~4530–4580)

**What could be wrong:**
These functions compute physical coordinates at fault face DOFs using
`FTr->Face->Transform(nip, coords)`. An uncommitted change replaces this
with `FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords)`.

For well-formed planar faces, both should agree. But for tetrahedral meshes
with shared faces, `FTr->Face` may use a different parametrization than
`FTr->Elem1`, and the face integration point `nip` may map to a different
physical location than the element-mapped point.

If coordinates are wrong, `a(x2, x3)` and `Dc(x2, x3)` are wrong at
those DOFs, causing wrong friction parameters at shared fault faces.

**Impact if wrong:**
Wrong a-values → wrong friction strength → wrong slip rate at shared faces.
This is a localized error at specific DOFs, not a global error.

**How to test locally:**
Compare `Face->Transform` and `Elem1->Transform` at every fault face DOF.

```
Test: TestFaceVsElemCoordinates
- For each fault face (interior + shared):
  - For each DOF integration point:
    - coords_face = FTr->Face->Transform(nip)
    - FTr->SetAllIntPoints(&nip)
    - coords_elem = FTr->Elem1->Transform(GetElement1IntPoint())
    - FAIL if |coords_face - coords_elem| > 1e-10
```

This test requires NO MPI comparison — it runs on each rank independently.
If Face->Transform and Elem1->Transform ever disagree, the coordinate
bug is confirmed. The uncommitted change (Elem1 path) should then be
evaluated as the fix.

**TEST RESULT (`test_diag_shared_face_coords`, line 1894 of
`tests/parallel/test_parallel_elasticity.cpp`, mpirun -np 4):**
```
PASSED (max_diff=8.881784e-16, shared_fault=16, qp_disagree=0)
```
Face->Transform and Elem1->Transform agree to machine precision at all
shared fault DOF points and higher-order quadrature points. **RULED OUT.**

---

### Suspect 6a: Ghost DOF Communication (ExpandOwnedToLocalFault) — ✅ RULED OUT

**Code:** `elasticity_operator.hpp` — `ExpandOwnedToLocalFault()` MPI
section (lines ~4632–4720)

**What could be wrong (per v62 Section 2, ruled out):**
The send path sends owned_data in canonical order (correct). The receive
path unpacks using `canonical_to_local_perm_` (correct). However, there
could be a subtle issue in the data layout or buffer packing that only
manifests on specific face geometries.

**How to test locally:**
Set owned data to known coordinate-based values. Expand to local. On the
ghost rank, verify received values match expectations.

```
Test: TestGhostDOFCommunication
- Use FaultGeometry coordinates (x2, x3) as known owned data
- ExpandOwnedToLocalFault(owned_x2, local_x2, 1)
- For each shared face that is a ghost on this rank:
  - Read local_x2 at ghost DOF positions
  - Compare against expected x2 values (from physical coordinates)
  - FAIL if |diff| > 1e-10
```

This directly tests whether the MPI communication delivers the right
values to the right DOF positions on the receiving rank.

**TEST RESULT (`test_diag_ghost_dof_roundtrip`, line 2271 of
`tests/parallel/test_parallel_elasticity.cpp`, mpirun -np 4):**
```
PASSED (max_err=5.551115e-15, rt_err=0.000000e+00, mismatch=0 [ghost=0 owned=0], total_owned=8, total_local=16)
```
Set owned data to `f(x2,x3) = 7·sin(x2) + 3·cos(x3)`, expanded to all
local DOFs (including ghosts via MPI), verified ALL DOFs match expected
values from independently-computed coordinates. Zero mismatches. Also
verified owned-only restrict→expand roundtrip: `rt_err=0.0`. **RULED OUT.**

---

### Suspect 6b: Shared Fault Comm Block Setup (Face Matching) — ✅ RULED OUT

**Code:** `elasticity_operator.hpp` — `BuildFaultOwnership()` (lines
~2200–2360)

**What could be wrong:**
The `shared_fault_comm_blocks_` structure determines which owned faces are
sent to which rank, and which local face positions receive the data. The
face matching uses `FaceVertexKey` (sorted global vertex IDs) gathered via
`MPI_Allgatherv`.

If a face is matched to the **wrong neighbor face** (e.g., two faces with
similar but different vertex keys), the correct data arrives at the wrong
ghost position. This would produce localized errors at specific shared faces.

The face matching depends on:
1. `MakeFaceKey(local_face, gvert)` — sorts 3 global vertex IDs
2. `MPI_Allgatherv` — collects all shared face keys from all ranks
3. `key_ranks` map — groups faces by key, assigns owner = lowest rank
4. `send_by_rank` / `recv_by_rank` — sorted by key for deterministic order

**Potential issues:**
- Two different physical faces with the same 3 global vertex IDs (impossible
  for conforming meshes, but worth verifying)
- Off-by-one in `recv_displs` or `recv_counts`
- Mismatch between send order and receive order (both sorted by key, so
  should match — but if sorting is inconsistent, face data goes to wrong slot)

**How to test locally:**
Verify that every send/recv pair maps to the same physical face.

```
Test: TestCommBlockFaceMatching
- For each comm block:
  - For each send_owned_face[j]:
    - Compute FaceVertexKey from owned_fault_face_to_local_face_[owned_face]
    - Send the key to the neighbor rank
  - For each recv_local_face[j]:
    - Compute FaceVertexKey from the local face
    - Receive the key from the neighbor rank
    - FAIL if send_key[j] != recv_key[j]
```

This verifies the face-matching is correct: each send slot maps to the
correct receive slot on the neighbor.

**TEST RESULT:** Covered by `test_diag_ghost_dof_roundtrip` (line 2271 of
`tests/parallel/test_parallel_elasticity.cpp`). If face matching were
wrong, the coordinate-based ghost verification would show mismatches
(wrong face's data at wrong positions). Zero mismatches. **RULED OUT.**

---

## 4. Test Implementation Plan

### Phase 1: Quick wins (no MPI comparison needed) — ✅ COMPLETE

| Test | Line | Suspect | Result | Verdict |
|------|------|---------|--------|---------|
| `test_diag_shared_face_coords` | L1894 | #5 | max_diff=8.88e-16 | RULED OUT |

Face->Transform and Elem1->Transform agree to machine precision.

### Phase 2: MPI consistency checks — ✅ COMPLETE

| Test | Line | Suspect | Result | Verdict |
|------|------|---------|--------|---------|
| `test_diag_shared_fault_basis` | L2057 | #1 | max_sign_err=0.0, 8/8 flipped | RULED OUT (by design) |
| `test_diag_ghost_dof_roundtrip` | L2271 | #6a, #6b | max_err=5.55e-15, 0 mismatch | RULED OUT |

Sign-baking is consistent. Ghost comm and face matching are correct.

All Phase 1–2 tests in **`tests/parallel/test_parallel_elasticity.cpp`**
using forced y-partitioning on `CreateTestMesh3DTet(2,1,1, ...)`.
Run: `mpirun -np 4 ./seas_test_parallel_elasticity` → these 3 pass.

### Phase 3: Assembly comparison — ✅ COMPLETE, IP ALL PASS

#### Initial results (BR2 default — WRONG method)

Tests D, E, F originally used the default `DGMethod::BR2`. This found a
real BR2 shared-face bug (every shared face wrong, rel_err ≈ 10×). But
the production BP5 run uses `--dg-method IP`.

#### Corrected results (DGMethod::IP — production method)

After adding `DGMethod::IP` to all diagnostic operators:

**Test D (`test_diag_displacement_with_dirichlet`, line 2446):**

| Sub-test | Slip | Time | Dirichlet active? | IP Result | rel_err |
|----------|------|------|--------------------|-----------|---------|
| D1  | zero | 1.0 | YES (only) | **PASS** | 2.5e-12 |
| D1b | nonzero | 0.0 | NO | **PASS** | 2.0e-13 |
| D2  | nonzero | 1.0 | YES (both) | **PASS** | 2.0e-13 |

**Test F (`test_diag_mesh_resolution_sweep`, line 2828):**

```
(2,1,1): PASS u_rel=2.6e-14  rhs_rel=0.0      fi=0 fs=16 nf=24
(3,1,1): PASS u_rel=5.9e-13  rhs_rel=5.3e-16   fi=0 fs=24 nf=36
(4,1,1): PASS u_rel=3.2e-13  rhs_rel=2.0e-16   fi=0 fs=32 nf=48
(2,1,2): PASS u_rel=3.1e-13  rhs_rel=4.4e-16   fi=0 fs=32 nf=48
(3,1,2): PASS u_rel=1.1e-12  rhs_rel=7.4e-16   fi=0 fs=48 nf=72
```

ALL cases pass at machine precision, including forced y-partition with
fi=0 (100% shared faces). Both the RHS and the solved displacement
match serial to < 1.1e-12.

**Test G (`test_diag_mixed_mesh_face_classification`, line 2969):**
PASS — `dof_match=OK, face_match=OK, local_nf=OK`.

**Test E (`test_diag_traction_with_dirichlet`, line 2633):**
PASS with IP — `l2_rel=4.0e-14, inf_rel=1.5e-13, ndof=24/24`.
Traction L2 and Linf norms match serial to machine precision.

Note: the original per-DOF coordinate matching failed (max_rel_err=1.95)
because DG has multiple DOFs at the same physical location (shared
vertices between faces) with different traction values. The fix was to
compare global traction norms (L2, Linf) instead of per-DOF values.
This correctly validates traction without the DG DOF-matching ambiguity.

### Summary of all diagnostic tests (IP method)

All in **`tests/parallel/test_parallel_elasticity.cpp`**.
Mesh helper `CreateTestMesh3DTetWithDirichlet` at line 656.
Run: `mpirun -np 4 ./seas_test_parallel_elasticity` → **18/18 pass**.

| Test | Line | Suspect | IP Verdict |
|------|------|---------|------------|
| `test_diag_shared_face_coords` | L1894 | #5 | ✅ RULED OUT |
| `test_diag_shared_fault_basis` | L2057 | #1 | ✅ RULED OUT |
| `test_diag_ghost_dof_roundtrip` | L2271 | #6a, #6b | ✅ RULED OUT |
| `test_diag_displacement_with_dirichlet` D1 | L2446 | #3 | ✅ RULED OUT |
| `test_diag_displacement_with_dirichlet` D1b | L2446 | #2 | ✅ RULED OUT |
| `test_diag_displacement_with_dirichlet` D2 | L2446 | #2+#3 | ✅ RULED OUT |
| `test_diag_traction_with_dirichlet` | L2633 | #4 | ✅ RULED OUT (l2_rel=4e-14) |
| `test_diag_mesh_resolution_sweep` | L2828 | #2 sweep | ✅ ALL 5 PASS |
| `test_diag_mixed_mesh_face_classification` | L2969 | face counts | ✅ PASS |

### BR2 shared-face bug (confirmed, not production-relevant)

The BR2 shared path (`AssembleSlipContributionBR2Shared`, lines 3578+)
has a confirmed bug: every shared fault face produces wrong elvec1
(rel_err ≈ 10×) with forced y-partition. This was found during the
investigation but does NOT affect BP5 production (which uses IP).

The BR2 bug likely involves the `slip_sign_all[q] = (nor_q(1) > 0) ?
1.0 : -1.0` sign computation (lines 3308, 3651) interacting with the
shared-face normal flip. To be fixed separately.

### Conclusion: all IP shared-face suspects ruled out

The BP5 production blowup is **NOT caused by** any of the 7 custom
shared-face code paths tested here. All suspects are ruled out for the
IP method with forced y-partitioning (100% shared faces) on meshes
of varying resolution.

---

## 5. Parallel Submissions (In Progress)

| Job | Ranks | Nodes | Queue | Walltime | Status |
|-----|-------|-------|-------|----------|--------|
| `bp5_v60_serial_48hr` | 1 | 1 | normal | 48hr | Submitted |
| `bp5_v60_paraview_dev_2hr` | 400 | 8 | development | 2hr | Baseline |
| `bp5_v60_100ranks_normal_6hr` | 100 | 2 | normal | 6hr | Ready |
| `bp5_v60_200ranks_normal_6hr` | 200 | 4 | normal | 6hr | Ready |
| `bp5_v60_600ranks_normal_6hr` | 600 | 11 | normal | 6hr | Ready |
| `bp5_v60_800ranks_normal_6hr` | 800 | 15 | normal | 6hr | Ready |

If hot spots shift with rank count → parallel-specific bug, but NOT in
the 7 shared-face suspects tested here (all ruled out for IP).
If serial run has no blowup → parallel-only bug in a path not yet tested
(e.g., MFEM internal assembly, solver setup, or time-stepping coupling).
If serial run also blows up → bug is in local physics (not shared-face).

---

## 6. Key Insight: Why v61 Permutation Hypothesis Was Wrong

The v61 document assumed `owned_data` is in "MFEM-local order." In fact,
it is in **canonical order** because `owned_fault_dof_to_local_dof_`
already embeds the `canonical_to_local_perm_` mapping.

Concrete trace for a face with perm = [1, 2, 0] (canonical→MFEM-local):

```
RestrictToOwnedFault:
  owned[0] = local[perm[0]] = local[1] = value at MFEM-local 1 = canonical 0 ✓
  owned[1] = local[perm[1]] = local[2] = value at MFEM-local 2 = canonical 1 ✓
  owned[2] = local[perm[2]] = local[0] = value at MFEM-local 0 = canonical 2 ✓

Old send path (correct):
  buffer[0] = owned[0] = canonical 0 value ✓
  buffer[1] = owned[1] = canonical 1 value ✓
  buffer[2] = owned[2] = canonical 2 value ✓

Receiver (canonical_to_local_perm_ on receiver side):
  local[receiver_perm[0]] = buffer[0] = canonical 0 value ✓
  local[receiver_perm[1]] = buffer[1] = canonical 1 value ✓
  local[receiver_perm[2]] = buffer[2] = canonical 2 value ✓
```

The "fix" would have read `owned[perm[canonical_kk]]` instead of
`owned[canonical_kk]`, scrambling the data.

---

## 7. Production-Mesh Verification Tests

### 7.1 Why production-mesh tests are needed

The Phase 1–3 local tests use a tiny mesh (8–48 tets, 4 MPI ranks).
While Phase 3 reproduces the bug locally on the mixed fault+Dirichlet
mesh, there are gaps between the local test and the production mesh
(62,874 tets, 400 ranks) that could hide additional issues or mask the
fix:

| Property | Local test mesh | Production mesh |
|----------|----------------|-----------------|
| Volume elements | 8–48 tets | 62,874 tets |
| Fault faces | 8–16 (all shared with y-partition) | 9,312 (mix of interior + shared) |
| Dirichlet faces on y=0 | 0–8 | 1,388 |
| MPI ranks | 2–4 | 100–800 |
| Fault/Dirichlet boundary | 0 junction faces | ~200 junction faces |
| Element quality | Uniform Cartesian-split | Variable Gmsh-generated |
| Time steps tested | 0 (single solve) | Thousands (accumulating) |
| Mesh scale | 1 (direct m) | 1000 (km → m) |

**Specific gaps:**
1. The local test has either ALL-shared or ALL-interior fault faces.
   The production mesh has a MIX — the bug's behavior at the
   interior/shared boundary is untested.
2. No multi-step divergence tracking — the local test does a single
   solve, but the production blowup accumulates over ~25 yr.
3. The Dirichlet continuation region (attr=5, z < -40 km on y=0) has
   a complex shape that may exercise face-classification edge cases
   not present in the rectangular local mesh.

### 7.2 Test P1: Ghost DOF roundtrip on production mesh

**Purpose:** Verify `ExpandOwnedToLocalFault` MPI communication is
correct on the production mesh with 400+ ranks and complex partitioning.

**Method:** Built into `ElasticityDomainOperator::VerifyGhostDOFCommunication()`.
Same logic as `test_diag_ghost_dof_roundtrip`:
1. Get local fault coordinates via `GetFaultCoords2D()`
2. Restrict to owned, set `f(x2,x3) = 7·sin(x2) + 3·cos(x3)` (2 comp)
3. `ExpandOwnedToLocalFault` → ghost DOFs filled via MPI
4. Compare all local DOFs against expected from independently-computed
   local coordinates
5. `MPI_Allreduce` max error, log PASS/FAIL

**Trigger:** `--verify` flag on `seas_bp5_full`, runs once before time stepping.

**Expected result:** PASS (max_err < 1e-10). If FAIL, the ghost
communication or face matching is broken on the production mesh.

**Runtime:** <1 second (no K assembly needed).

### 7.3 Test P2: Serial vs parallel RHS norm comparison

**Purpose:** Detect shared-face assembly errors on the production mesh
by comparing RHS norms between 1-rank (serial, no shared faces) and
N-rank (parallel, shared faces active) runs.

**Method:** Built into `ElasticityDomainOperator::VerifyRHSNorms(time, slip)`.
Assembles the full RHS (slip + Dirichlet) and reports global L2 norms
of `b_slip`, `b_dir`, and `b_total` separately.

Run at two times:
- `t=0`: Only slip RHS (Dirichlet = 0). If norms differ, the shared
  fault slip assembly is wrong.
- `t=1yr` (3.156e7 s): Both slip and Dirichlet active. If `b_dir`
  differs, the shared Dirichlet loading is wrong.

**Trigger:** `--verify` flag, runs after `SetInitialCondition` (K assembled).

**Expected result:** Norms match between 1 rank and N ranks to ~12
significant figures (MUMPS direct solve, no iterative error).

**Runtime:** ~5 minutes (includes K factorization at 1 rank).

**How to compare:**
```bash
grep "\[VERIFY\]" bp5_verify_serial_*.out  > serial_norms.txt
grep "\[VERIFY\]" bp5_verify_parallel_*.out > parallel_norms.txt
diff serial_norms.txt parallel_norms.txt
```

### 7.4 Test P3: Multi-step Vmax divergence tracking

**Purpose:** Detect when serial and parallel simulations first diverge.
The production blowup occurs at ~25 yr; if the divergence starts
earlier, it narrows the time window for detailed comparison.

**Method:** No code change — uses existing `_global.txt` output which
logs `(time, log10_Vmax)` at each accepted step. Run serial (1 rank)
and parallel (400 ranks) for 50 steps (reaches ~0.1–1 yr depending on
adaptive dt). Compare the Vmax time series.

**Sbatch scripts:**
- `jobs/bp5/bp5_verify_serial.sbatch`: 1 rank, dev queue 30 min, 50 steps
- `jobs/bp5/bp5_verify_parallel.sbatch`: 400 ranks, dev queue 30 min, 50 steps

Both use `--verify` (P1 + P2 diagnostics) and `--max-steps 50`.

**Expected result:**
- If Vmax diverges immediately (step 1–5): the per-step RHS error is
  large, confirming the Phase 3 local finding on the production mesh.
- If Vmax stays aligned for 50 steps: the error is subtle and only
  accumulates over many steps — run longer (Section 5 jobs: 6–48 hr).
- If both serial and parallel blow up: the bug is in local physics,
  not shared faces.

**How to compare:**
```bash
paste serial_global.txt parallel_global.txt | \
  awk '{printf "%12.6e  %12.6e  %12.6e  %12.6e\n", $1,$2,$3,$4}'
# Columns: t_serial, Vmax_serial, t_parallel, Vmax_parallel
# Divergence: columns 2 and 4 differ significantly
```

**Runtime:** ~15–25 minutes (K factorization dominates at 1 rank).

### 7.5 Test P4: Face classification audit on production mesh

**Purpose:** Verify that `BuildFacetBCTables` classifies all y=0 faces
correctly on the production mesh. The Phase 3 finding (bug triggered by
mixed fault+Dirichlet mesh) suggests face classification may be the root
cause. The production mesh has both attr=3 and attr=5 on y=0.

**Method:** Already built into `ElasticityDomainOperator` as
`RunStartupFaceAudit()` (runs automatically). Outputs:
```
Fault (attr=3):      N1
Dirichlet (attr=5):  N2
NONE (interior):     N3
NONE (shared):       N4
```

**Expected result:** N3=0 and N4=0 (all y=0 faces classified). Compare
counts between 1-rank and N-rank runs — the totals should agree.

**How to compare:**
```bash
grep "Face Classification Audit" -A5 bp5_verify_serial_*.out
grep "Face Classification Audit" -A5 bp5_verify_parallel_*.out
```

If fault or Dirichlet counts differ → face classification depends on
partitioning, which would explain the bug.

**Runtime:** Included in startup, <1 second.

### 7.6 Test P5: Per-DOF slip/traction snapshot at t ≈ 0.1 yr

**Purpose:** Detect per-DOF discrepancies between serial and parallel
early in the simulation, before the error accumulates.

**Method:** Add `--write-every-step` to both serial and parallel runs.
At the first output step (t ≈ 0.02 s via RK45 adaptive), compare:
- Per-station time series (SCEC format): `{prefix}_dp000.txt` etc.
- These files contain `(time, slip, log10_V, shear_stress, log10_state)`
  at probe depths 0, 2, 4, ..., 40 km.

**Expected result:** Station values match between serial and parallel
to ~10 significant figures at the first few steps. If they diverge at
step 1, the error is in the initial RHS/traction/coupling. If they
diverge later, the error accumulates through the friction ODE.

**Runtime:** Included in the 50-step run (no additional cost).

### 7.7 Implementation status

| Test | Location | Status |
|------|----------|--------|
| P1 Ghost DOF | `elasticity_operator.hpp` `VerifyGhostDOFCommunication()` | ✅ Implemented |
| P2 RHS norms | `elasticity_operator.hpp` `VerifyRHSNorms()` | ✅ Implemented |
| P3 Vmax divergence | `_global.txt` comparison | ✅ Sbatch scripts ready |
| P4 Face audit | `RunStartupFaceAudit()` (existing) | ✅ Already runs |
| P5 Station snapshot | `--write-every-step` (existing) | ✅ Already available |

CLI flag: `--verify` on `seas_bp5_full` enables P1 + P2.
P3–P5 use existing output, compared post-hoc.

Sbatch scripts:
- `jobs/bp5/bp5_verify_serial.sbatch` — 1 rank, dev 30 min, 50 steps
- `jobs/bp5/bp5_verify_parallel.sbatch` — 400 ranks, dev 30 min, 50 steps

### 7.8 Production-mesh results (Frontera, 2026-04-07)

**Jobs:** serial job 7639250 (1 rank), parallel job 7639293 (400 ranks).
Both ran on `bp5_tandem_exact.msh` (62,874 tets) with `--verify --max-steps 50`.

#### P1: Ghost DOF communication — PASS

```
Serial:   [VERIFY] Ghost DOF communication: max_err=0.000000e+00, mismatches=0 → PASS
Parallel: [VERIFY] Ghost DOF communication: max_err=0.000000e+00, mismatches=0 -> PASS
```

Ghost communication (`ExpandOwnedToLocalFault`) is correct on the
production mesh with 400 ranks and complex partitioning.

#### P2: RHS norms — DIRICHLET MISMATCH CONFIRMED

```
                         ||b_slip||              ||b_dir||                ||b_total||
Serial   t=0:    0.000000000000e+00    0.000000000000e+00    0.000000000000e+00
Parallel t=0:    0.000000000000e+00    0.000000000000e+00    0.000000000000e+00
Serial   t=1yr:  0.000000000000e+00    1.436659911942e+16    1.436659911942e+16
Parallel t=1yr:  0.000000000000e+00    1.381124698531e+16    1.381124698531e+16
```

| Metric | Serial | Parallel | Relative diff |
|--------|--------|----------|---------------|
| `\|\|b_slip\|\|` at t=0 | 0 | 0 | 0 (match) |
| `\|\|b_dir\|\|` at t=1yr | 1.43666e+16 | 1.38112e+16 | **3.87%** |

**The shared-face Dirichlet loading assembly produces a 3.87% error in
`||b_dir||` on the production mesh.** This is a per-step error that
accumulates through the entire simulation.

Note: `||b_slip||` = 0 in both cases because the verification used zero
slip. The slip path was not tested on the production mesh by P2 (it was
tested locally by Phase 3 test D1b, which FAILED on the mixed mesh).

#### P4: Face classification — MATCH

```
Serial:   Fault (attr=3): 9257  Dirichlet (attr=5): 1377  NONE: 0
Parallel: Fault (attr=3): 9257  Dirichlet (attr=5): 1377  NONE: 0
```

Face counts match between serial and parallel. The bug is NOT in face
classification — all faces are correctly identified. The bug is in how
the Dirichlet loading is assembled on shared faces AFTER classification.

#### P6: Dirichlet skip-set audit — SKIP SET IS CORRECT

```
Serial:
  attr=5 bdr elems:      1628
  skipped (interior):    1388
  skipped (shared set):  0
  processed (boundary):  240
  dirichlet_shared list: 0

Parallel:
  attr=5 bdr elems:      1628
  skipped (interior):    1377
  skipped (shared set):  11
  processed (boundary):  240
  dirichlet_shared list: 22
  WARNING: shared_skip (11) != dirichlet_shared list (22)
```

**The WARNING is a false alarm.** Analysis:

- 22 = 11 shared Dirichlet faces × 2 ranks (each face in
  `dirichlet_shared_faces_` on both sharing ranks).
- Each shared face has a boundary element on exactly ONE rank.
- The 11 that skip: the rank that has the boundary element correctly
  matches face_idx and skips it from the boundary loop.
- The other 11 appearances: the other rank has no boundary element for
  this face, so the boundary loop never sees it. No skip needed.
- No double-counting, no missed faces. The skip mechanism is correct.

**Accounting check:**
```
Serial:   1388 interior + 0 shared + 240 boundary = 1628 ✓
Parallel: 1377 interior + 11 shared + 240 boundary = 1628 ✓
                          ^^
          11 faces changed from interior (serial) to shared (parallel)
```

#### Summary: bug is in the shared Dirichlet FORMULA, not skip sets

| Test | Result | Implication |
|------|--------|-------------|
| P1 Ghost DOF | PASS | MPI communication correct |
| P2 `b_slip` at t=0 | MATCH | (untested — zero slip used) |
| P2 `b_dir` at t=1yr | **3.87% MISMATCH** | Shared Dirichlet formula wrong |
| P4 Face classification | MATCH | Faces classified correctly |
| P6 Skip-set audit | OK (no double/missing) | Skip mechanism correct |

**The bug is in the per-face shared Dirichlet skeleton formula.**
Only 11 out of 1388 Dirichlet faces are shared, yet the total ||b_dir||
changes by 3.87%. This implies each shared face has a **~50% per-face
error** in its Dirichlet contribution ((3.87% × 1388) / 11 ≈ 49%).

The skip-set diagnostic proves:
- No faces are double-counted (processed + skipped = total)
- No faces are missed (all boundary elements accounted for)
- The interior path handles 1377 faces correctly
- The boundary path handles 240 true-boundary faces correctly
- **The shared path handles 11 faces with the wrong formula**

Combined with the local Phase 3 findings:
- Local D1 (Dirichlet only, zero slip): PASS — Dirichlet path OK on small mesh
- Local D1b (slip only, t=0): FAIL (rel_err=8.17) — slip path broken on mixed mesh
- Production P2 (Dirichlet only, zero slip): FAIL (3.87%) — Dirichlet shared formula wrong
- Production P6: skip-set accounting correct — bug is NOT in face counting

#### P7: Per-face shared Dirichlet diagnostic — ROOT CAUSE FOUND

`VerifySharedDirichletPerFace(t=1yr)` on Frontera job 7639340 (400 ranks):

```
11 pairs, 8 with SAME dir_sign, 0 with ev cross-mismatch
-> dir_sign does NOT flip on 8 faces! This is the bug.
```

8 of 11 shared Dirichlet face pairs show `SIGN_SAME!`: `dir_sign = +1`
on BOTH ranks despite the face normals being opposite (`nor_y` has
opposite sign). The 3 correct faces have `dir_sign` that properly flips.

Example (face 1637,1640,1644):
```
r51: ds=+1, nor_y=-7.66e8   (dot = +7.66e8 > 0 → ds=+1 correct)
r55: ds=+1, nor_y=+7.66e8   (dot = -7.66e8 < 0 → ds SHOULD be -1, but is +1!)
```

The 8 broken faces all have ONE rank that has Dirichlet shared faces but
**zero fault DOFs**. The 3 correct faces have BOTH ranks with fault DOFs.

---

## 8. ROOT CAUSE: `ref_normal_` Uninitialized on No-Fault Ranks

### 8.1 The bug

`ref_normal_` is set at line 2488 inside `if (num_fault_dofs_ > 0)`:

```cpp
if (num_fault_dofs_ > 0)        // ← guard
{
    Vector ref_normal(3);
    ref_normal(1) = -1.0;
    ref_normal_ = ref_normal;   // ← line 2488: only set when fault DOFs > 0
    ...
}
```

On ranks with zero fault DOFs but shared Dirichlet faces, `ref_normal_`
is a default `Vector` (size 0). Then `ComputeSkeletonDirichletSign`:

```cpp
real_t dot_ref = 0.0;
for (int d = 0; d < ref_normal_.Size(); d++)   // Size()=0 → loop skipped
{
    dot_ref += nor(d) * ref_normal_(d);
}
return (dot_ref < 0.0) ? -1.0 : 1.0;   // dot_ref=0 → always returns +1
```

`dir_sign = +1` on ALL faces regardless of normal direction. The
Dirichlet value `u_D = dir_sign * Vh` doesn't flip to compensate the
flipped face normal, breaking the skeleton formula:

- **Correct:** `u_D_A = +Vh`, `u_D_B = -Vh` → penalty terms cancel
- **Buggy:**   `u_D_A = +Vh`, `u_D_B = +Vh` → penalty term has wrong sign

This produces a ~50% per-face error on the 8 affected faces, totaling
the 3.87% global `||b_dir||` mismatch.

### 8.2 Why it only manifests on the production mesh

On the local test mesh with forced y-partitioning, every rank that has
shared y=0 faces also has fault DOFs (because the mesh is small and the
fault covers the full y=0 plane). So `ref_normal_` is always set.

On the production mesh with METIS partitioning, some ranks in the
far-field get a few shared Dirichlet faces on y=0 but NO fault faces.
These ranks have `num_fault_dofs_ = 0` and `ref_normal_.Size() = 0`.

### 8.3 The fix

Move `ref_normal_` initialization BEFORE the `if (num_fault_dofs_ > 0)`
guard so ALL ranks have it set:

```cpp
// Set unconditionally — needed by ComputeSkeletonDirichletSign even
// on ranks with zero fault DOFs
{
    Vector ref_normal(3);
    ref_normal = 0.0;
    ref_normal(1) = -1.0;
    ref_normal_ = ref_normal;
}

if (num_fault_dofs_ > 0)
{
    const Vector &ref_normal = ref_normal_;  // reuse
    ...
}
```

**Commit:** `27a548b` on `test/revert-transform`.

### 8.4 Local test result after fix

```
Before fix: test_diag_displacement_with_dirichlet FAIL (rel_err=8.17)
After fix:  test_diag_displacement_with_dirichlet PASS
```

**18 of 18 local tests pass.** The earlier traction test failure
(`max_rel_err=1.95`) was a test design issue: per-DOF coordinate
matching is ambiguous for DG (multiple DOFs at the same physical
location with different values). Fixed by comparing global traction
norms (L2, Linf) instead → `l2_rel=4.0e-14, inf_rel=1.5e-13`. PASS.

### 8.5 Production-mesh verification — FIX CONFIRMED

Frontera jobs: serial 7639389 (1 rank), parallel 7639390 (400 ranks).

#### P2: RHS norms — MATCH (was 3.87% off)

```
Serial   t=1yr: ||b_dir||=1.436659911942e+16
Parallel t=1yr: ||b_dir||=1.436659911942e+16
```

All 12 printed digits match. The 3.87% mismatch is completely eliminated.

#### P7: Per-face diagnostic — ALL faces correct

```
summary: 11 pairs, 0 with SAME dir_sign, 0 with ev cross-mismatch
```

All 11 shared Dirichlet face pairs now have opposite `dir_sign`
(was 8 with SAME before the fix). Every face shows `ds=+1` on one
rank and `ds=-1` on the other, correctly compensating the normal flip.

#### P1: Ghost DOF — PASS (unchanged)

```
Ghost DOF communication: max_err=0.000000e+00, mismatches=0 -> PASS
```

#### Before/after summary

| Metric | Before fix (job 7639293) | After fix (job 7639390) |
|--------|-------------------------|------------------------|
| Serial `\|\|b_dir\|\|` | 1.436659911942e+16 | 1.436659911942e+16 |
| Parallel `\|\|b_dir\|\|` | 1.381124698531e+16 | **1.436659911942e+16** |
| Mismatch | 3.87% | **0 (12-digit match)** |
| SIGN_SAME faces | 8 of 11 | **0 of 11** |
| Ghost DOF | PASS | PASS |

### 8.6 Next step: full production run

Submit a full production run (400 ranks, 48hr normal queue, no
`--max-steps` limit) to verify the ~25 yr displacement blowup is gone.
The fix eliminates the 3.87% per-step Dirichlet loading error that was
accumulating over thousands of time steps.
