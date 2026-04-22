# Implementation Plan: Local Unit Tests to Constrain Pepper-Bug Hypotheses H1–H5

## Overview
The `tpv102_mfem_vs_seissol_codeflow_pepper_report_2026-04-22.md` document lists five plausible
mechanisms (H1–H5) for the triangle-local "pepper" pattern observed on TPV102 fault outputs.
Per user directive *no TPV102 reproducer runs locally, no Frontera sbatch until asked*, this
plan specifies unit tests that can falsify each hypothesis on a laptop (≤ 8-element fixtures,
MPI np ≤ 4) without any production-scale simulation. All tests produce a binary PASS/FAIL
decision and, on FAIL, print the exact field/component that diverged.

Target files: `miniapps/seas/tests/unit/` and `miniapps/seas/tests/parallel/`, wired into
`miniapps/seas/Makefile` in the same style as the existing `seas_test_interior_fault_flux_path`
and `seas_test_shared_fault_dof_data_consistency` targets.

---

## Local Testability Summary (scorecard)

| Hyp | Existing tests that cover it                                                                                                      | Gap                                                                                                    | Local-testable? | New test(s) in this plan  |
|-----|-----------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------|:---------------:|---------------------------|
| H1  | `test_fault_flux_interior_vs_shared_branch_equivalence.cpp` (header admits TAUTOLOGY post-v9.0.0); `test_r101_shared_fault.cpp`   | No live cross-check of interior-branch vs shared-branch rhs on the **same physical face**              |       YES       | Phase 1 (P1-T1)           |
| H2  | `test_canonical_rotation_pure_strikeslip.cpp`, `test_per_qp_vs_centroid_basis_planar_1el.cpp`, `test_godunov_rotation_identity`   | No truth-table test that `sign_flipped` **and** `elem1_on_plus` are **jointly consistent** per branch  |       YES       | Phase 2 (P2-T1, P2-T2)    |
| H3  | `test_interior_fault_flux_path.cpp` (only anti-symmetric Q); `test_fault_face_flux*.cpp` (skips Q-interpolation)                  | No **asymmetric** Q probe; no shared-fault analog of interior flux-path test                           |       YES       | Phase 3 (P3-T1, P3-T2)    |
| H4  | `test_shared_fault_dof_data_consistency.cpp` Phase D (τ_ini·tanh(y/L0)), `test_shared_fault_role_consistency.cpp`                 | No unit test of `MakeFaceKey` uniqueness across ranks; no DOF-offset overlap assertion                 |       YES       | Phase 4 (P4-T1, P4-T2)    |
| H5  | `test_no_penalty_dynamic_rupture.cpp`, `test_volume_jacobian_single_channel.cpp`, `test_wave_operator_spatial_derivative.cpp`     | No analytic plane-wave pattern-recovery test; no RK4 stability probe for per-DOF noise                 |       YES       | Phase 5 (P5-T1, P5-T2)    |

All five hypotheses **can** be constrained by local unit tests. The plan proceeds phase-by-phase
from the hypothesis closest to the symptom (H3, fault-side Q trace) to the most upstream (H5, bulk).

---

## Constraints

- **Fixture size**: ≤ 8 hex/tet elements per test, ≤ 3 fault QPs. Total runtime per test ≤ 2 s.
- **MPI ranks**: serial for interior-branch tests; `mpirun -np 2` for shared-branch tests; one
  `np=4` test for DOFData consistency. No `-np > 4` anywhere — matches user's local laptop cap.
- **No TPV102 production mesh** — per user memory, don't run even 1000 m reproducer locally.
- **Build system**: every new `.cpp` gets (a) a `*_SRC` / `*_OBJ` block in the Makefile, (b) a
  `seas_test_*` link rule, (c) a `test-*` target that executes it, (d) inclusion in the
  `test` / `test-tpv102-pepper` aggregate target (new — see Phase 6).
- **Interface preservation**: no changes to `WaveOperator`, `FaultFaceFlux`, `FaultBasis`, or
  `GodunovFlux` **public** signatures. Where a test needs access to internal state
  (`fault_interior_face_to_basis_idx_`, `shared_fault_elem1_on_plus_`,
  `fault_face_dof_offset_`), either (i) add a `const &` getter guarded by
  `#ifdef SEAS_TEST_INTERNAL` (preferred, keeps production header untouched in non-test
  builds), or (ii) have the test call `WaveOperator::Mult` and probe indirectly via DOFData.
  Default choice: use getters with `#ifdef SEAS_TEST_INTERNAL`; `-DSEAS_TEST_INTERNAL` is added
  to every new test's compile flags (`CCC` in the Makefile) but **not** to production
  drivers. This matches the pattern already used by `SEAS_DIAG_FAULT_FLUX`.
- **Diagnostic output**: every test prints `{PASSED, FAILED}` per assertion and the file/line
  of the failing check; non-zero exit code on any FAIL.
- **No modifications to BP5 source** — per user memory, pepper fixes live in `dynamic/` only.
  If any production code changes are needed (e.g. to expose a getter), only
  `dynamic/wave_operator.{hpp,inl}` and `dynamic/fault_face_flux.{hpp,cpp}` are touched.
- **Numerical precision**: use `real_t` (= `double`) throughout. Tolerances specified
  per-test; prefer bit-exact (tol = 0.0) where the algebra permits (signed-permutation
  rotations on analytic inputs).

---

## Phase 1: Live interior-vs-shared branch cross-check (H1)

### Goal
After this phase: one new test drives a bulk state through BOTH the interior-fault branch
(serial, 2-tet fixture with one interior face at y=0) AND the shared-fault branch (np=2,
same 2-tet fixture partitioned to one tet per rank) and compares `rhs[Elem1]` contributions
at the fault QP bit-exactly (tol ≤ 1e-12 × |rhs|). The existing
`test_fault_flux_interior_vs_shared_branch_equivalence.cpp` is replaced per its own header
directive.

### Files to Create
- `tests/unit/test_interior_vs_shared_branch_live.cpp` — serial + parallel section, the
  MPI section compiles out in serial builds (`#ifdef MFEM_USE_MPI`). Serial section builds
  the serial `Mesh`; parallel section builds a `ParMesh` on `MPI_COMM_WORLD` with
  `Mesh::Partition` assigning tet 0 to rank 0, tet 1 to rank 1. Both sections construct the
  same `WaveOperator`, set the same `Q`, call `Mult(Q, k)`, and extract
  `k[VX · ndof_total + dof_offset_elem1 + 0]`. The shared branch is exercised on rank 0 only;
  rank 1 is a mirror.

### Files to Modify
- `dynamic/wave_operator.hpp` — add `const std::unordered_map<int,int>& GetFaultInteriorFaceToBasisIdx() const`,
  `const std::vector<bool>& GetSharedFaultElem1OnPlus() const`,
  `const std::unordered_map<int,int>& GetFaultFaceDofOffset() const`,
  `const std::unordered_map<int,int>& GetSharedFaultDofOffset() const`. All four under
  `#ifdef SEAS_TEST_INTERNAL`. No change to member visibility in non-test builds.
- `Makefile` — add `TEST_INTERIOR_VS_SHARED_LIVE_SRC`, `_OBJ`, link rule
  `seas_test_interior_vs_shared_branch_live`, and `test-interior-vs-shared-branch-live`
  target. Add `-DSEAS_TEST_INTERNAL` to this test's compile flags via
  `CXXFLAGS += -DSEAS_TEST_INTERNAL` in a per-target override block.
- `Makefile` — deprecate `test-fault-flux-interior-vs-shared-branch-equivalence` target
  (rename to `test-fault-flux-interior-vs-shared-branch-equivalence-LEGACY` and stop including
  it in `test` aggregate); the existing test's own header already directs this
  ("REPLACE this test with a live cross-check").

### Detailed Requirements
1. **Fixture (2-tet, one interior fault face)**:
   - Vertices: V0=(0,-1,0), V1=(0,0,0), V2=(1,0,0), V3=(0,0,1), V4=(0,+1,0).
   - Tet 0 = {V0, V1, V2, V3} (centroid y=-0.25, the − side).
   - Tet 1 = {V4, V1, V2, V3} (centroid y=+0.25, the + side).
   - Interior face = {V1, V2, V3} at y=0; **attribute 3** (fault); all other triangles get
     attribute 1 (free-surface).
   - Builder: `Mesh BuildTwoTetFault()` returns the mesh; identical code path in both sections.
2. **Material**: homogeneous TPV102 (`TPV102Params::lambda, mu, rho`). `FaultFaceFlux`
   constructed with those values.
3. **Bulk state**: `Q` zero EXCEPT `Q[SXY]` on tet 0 = +A, `Q[SXY]` on tet 1 = -A, where
   A = 1.0e-3 · τ_ini. This is the **asymmetric** pattern (NOT anti-symmetric in the
   velocity sense — the stress itself has a true jump across the fault). A is chosen small
   enough that the friction solver stays at pre-nucleation equilibrium but large enough that
   the trial traction is machine-distinguishable from pre-stress.
4. **Serial branch** (rank 0 only when parallel):
   - Construct `WaveOperator<Mesh> wave(mesh, fes, bc, ...)`.
   - Populate `DOFData` via `TPV102Setup::InitializeFaultDOFs(...)` — same call the driver uses.
   - Call `wave.Mult(Q, k)`. Record `k[c * ndof_total + dof_offset_tet0 + 0]` for
     c ∈ {SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ} into `k_serial[9]`.
5. **Parallel branch** (np=2):
   - Partition so rank 0 owns tet 0, rank 1 owns tet 1. The face {V1,V2,V3} becomes a shared
     fault face.
   - Construct `WaveOperator<ParMesh> wave_par(pmesh, pfes, bc, ...)`.
   - Populate DOFData identically.
   - Call `wave_par.Mult(Q_local, k_local)` on each rank. Record
     `k_local[c * ndof_total + dof_offset_tet0_rank0 + 0]` from rank 0 into `k_parallel[9]`.
6. **Assertion**: `|k_serial[c] - k_parallel[c]| ≤ 1e-12 * max(|k_serial[c]|, |k_parallel[c]|, 1.0)`
   for every c. If any component fails, the test name identifies the c and dumps both values.

### Interfaces
- `const std::unordered_map<int,int>& WaveOperator::GetFaultInteriorFaceToBasisIdx() const`
  (header-only, `#ifdef SEAS_TEST_INTERNAL`, returns `fault_interior_face_to_basis_idx_`).
- `const std::vector<bool>& WaveOperator::GetSharedFaultElem1OnPlus() const` (same pattern).
- `const std::unordered_map<int,int>& WaveOperator::GetFaultFaceDofOffset() const`.
- `const std::unordered_map<int,int>& WaveOperator::GetSharedFaultDofOffset() const`.

### Edge Cases to Handle
- **MPI-less build**: the parallel section is wrapped in `#ifdef MFEM_USE_MPI`. When built
  without MPI the test reports "skipped: MPI not built" and exits 0. This matches the pattern
  in `test_r101_shared_fault.cpp`.
- **Partitioner places both tets on one rank**: MFEM's default partitioner can collapse a
  2-element mesh onto rank 0 when np=2. Force `ParMesh` constructor with an explicit
  `partitioning[]` array: `{0, 1}`, guaranteeing one tet per rank. If this still yields zero
  shared fault faces on one rank, `MFEM_ABORT` with a clear message — this is a fixture bug,
  not a production-code bug.
- **Rank ordering of (+, −) side**: the serial branch's `elem1_on_plus` classification
  depends on geometric projection (interior branch uses `LookupInteriorFaultBasisIndex`);
  the parallel branch's depends on `shared_fault_elem1_on_plus_` (centroid projection in
  `wave_operator.inl:403-446`). Both must agree on "tet 0 is the − side". Assert this
  explicitly in the test before comparing `k_*`.
- **DOFData populated differently on the two branches**: the test populates DOFData once per
  `WaveOperator`, so the serial and parallel instances have independent copies. Both branches
  must observe the **same** DOFData immediately after `InitializeFaultDOFs`; assert this via
  `sigma_n_corr`, `tau2_corr` equality before calling `Mult`.

### Acceptance Criteria
- [ ] Serial section runs in < 200 ms.
- [ ] Parallel section runs in < 500 ms on np=2.
- [ ] All 9 components of `k_serial[c] == k_parallel[c]` to 1e-12 relative.
- [ ] Test links against `WAVE_OPERATOR_OBJ, GODUNOV_FLUX_OBJ, PML_LAYER_OBJ, FAULT_FACE_FLUX_OBJ, FRICTION_SOLVER_OBJ`.
- [ ] `make test-interior-vs-shared-branch-live` succeeds.
- [ ] Test's printed failure message names the failing component (SXY, VX, etc.) and
      numeric values on both sides when the assertion fires.
- [ ] Tampering demonstration: manually flip the `elem1_on_plus` classification in the
      parallel section (patch `shared_fault_elem1_on_plus_[0] = !shared_fault_elem1_on_plus_[0]`
      before the `Mult` call — gated by env var `SEAS_TEST_TAMPER=1`) — the test MUST FAIL.
      This replaces the tautology of the existing equivalence test.

### Dependencies
- Depends on: none (pure unit test).
- Required by: Phase 3 P3-T2 (which reuses the 2-tet fixture).

---

## Phase 2: sign_flipped / elem1_on_plus truth-table (H2)

### Goal
After this phase: two new tests verify the `sign_flipped → elem1_on_plus` logic in both
branches against a hand-built truth table. Specifically, for the interior branch
(`wave_operator.inl:1117`) `elem1_on_plus = !qpd.sign_flipped` must hold iff MFEM's
`CalcOrtho` on the same face gives a normal anti-aligned with `ref_normal = (0,-1,0)`
when the canonical + side sits on tet 0. For the shared branch
(`wave_operator.inl:1567-1568`) the same relation must hold via
`shared_fault_elem1_on_plus_`.

### Files to Create
- `tests/unit/test_sign_flipped_elem1_on_plus_truth_table.cpp` — serial, 4 fixture variants.
- `tests/parallel/test_sign_flipped_elem1_on_plus_truth_table_shared.cpp` — np=2, same
  4 fixture variants with shared-face geometry.

### Files to Modify
- `dynamic/fault/fault_basis.hpp` — no changes; expose `sign_flipped` via existing
  `FaultBasisQPData::sign_flipped` (already public).
- `Makefile` — add 2 link rules + test targets.

### Detailed Requirements
1. **Four fixture variants**, each exercising a different geometric combination:
   - **A**: 2-tet with interior face normal aligned with `ref_normal=(0,-1,0)` (tet 1 above).
     Expected: `sign_flipped = false` on at least one QP, `elem1_on_plus = true` for the upper tet.
   - **B**: Same mesh with vertex order of tets swapped (tet 0 ↔ tet 1 in the element table).
     Expected: `sign_flipped` flips accordingly; `elem1_on_plus` follows.
   - **C**: Rotated mesh with fault at x=0 (`ref_normal=(+1,0,0)`) — tests that the logic
     is not hard-coded to y=0.
   - **D**: Degenerate case: fault face exactly at axis-aligned orientation but with a
     near-zero-volume neighbor (one tet height = 1e-8). `sign_flipped` must still be a
     pure sign comparison, not a tolerance-sensitive decision.
2. **Assertion table** for each variant:
   - `FaultBasis::ComputeQPBasis(...)` produces `qp_data[q].sign_flipped` consistent with
     `dot(nor_mfem, ref_normal) < 0` exactly (`<` not `<=` — no ties in A/B/C;
     D's tie-breaking is by the existing library code and we accept whatever sign it gives
     as long as it is consistent across QPs on the same face).
   - `elem1_on_plus_reconstructed = !qpd.sign_flipped` (interior branch) matches the
     geometric truth computed from tet centroids projected onto `ref_normal`: `elem1_on_plus
     = (projection(centroid_elem1) > projection(centroid_elem2))`.
3. **Shared-branch companion** (`np=2`): same four fixtures partitioned across ranks.
   Assert `shared_fault_elem1_on_plus_[0]` (from the cached vector) matches the interior
   branch's `!qpd.sign_flipped` on the same physical face — i.e., that the two branches use
   the SAME definition of "+ side" for the same face.

### Interfaces
- No new public interface. Uses existing:
  - `FaultBasis::ComputeQPBasis`, `FaultBasisQPData::sign_flipped`, `normal`, `tangent1`,
    `tangent2` (all public in `fault_basis.hpp`).
  - `WaveOperator::GetSharedFaultElem1OnPlus()` (added in Phase 1).

### Edge Cases to Handle
- **Fault exactly parallel to ref_normal (degenerate)**: variant D is fragile; if the test
  fails on D, report as "degenerate geometry — skipped" rather than FAILED (set a
  `SKIPPED` counter). This avoids false positives on implementation-defined tie-breaking.
- **MFEM builds with single precision**: `sign_flipped` is computed from a dot-product
  whose sign under float (vs double) may differ near zero. The fixture A/B/C geometries
  give dot-products of magnitude ≥ 0.5, so this is not an issue in practice; assert
  `|dot(nor_mfem, ref_normal)| > 0.1` before using `sign_flipped`.
- **Same fixture, multiple QPs**: all QPs on a planar face must share the same
  `sign_flipped` (follows from constant Jacobian). Assert this per-fixture.

### Acceptance Criteria
- [ ] All 4 variants PASS on both serial and np=2 builds.
- [ ] The test prints a truth table: `(variant, qp_idx, sign_flipped, elem1_on_plus, agree?)`.
- [ ] Tampering: manually patch `qp_data[q].sign_flipped = !sign_flipped` for one QP via
      env var `SEAS_TEST_TAMPER=1` — the test MUST FAIL and name the failing QP.

### Dependencies
- Depends on: Phase 1 (for `GetSharedFaultElem1OnPlus` getter).
- Required by: nothing downstream; completes H2 coverage.

---

## Phase 3: Fault-side Q trace reconstruction (H3)

### Goal
After this phase: extend `test_interior_fault_flux_path.cpp` with **asymmetric** Q test
cases (current T1–T3 use anti-symmetric) and add a shared-fault analog. These two tests
together probe whether the per-side Q interpolation + rotation produces the canonical
Pelties eq. (7) trial traction when Q_plus ≠ −Q_minus.

### Files to Create
- `tests/unit/test_interior_fault_flux_path_asymmetric.cpp` — NEW, siblings of
  `test_interior_fault_flux_path.cpp`, same 2-tet fixture, different Q patterns.
- `tests/parallel/test_shared_fault_flux_path.cpp` — NEW, np=2 analog.

### Files to Modify
- `Makefile` — add 2 link rules + test targets.

### Detailed Requirements

1. **Interior asymmetric cases** (T1a–T3a, complementing the existing T1–T3):
   - **T1a (asymmetric v_y)**: Q[VY] = +V_test on tet 0 (minus side), Q[VY] = 0 on tet 1
     (plus side). Pelties 7a: `sigma_n_trial = eta_p · (v_n⁻ - v_n⁺) = eta_p · (+V_test - 0)
     = eta_p · V_test = 0.5 · Zp · V_test`. Expected:
     `data.sigma_n_corr - sigma_n0 ≈ +0.5·Zp·V_test`. V_test = 1e-6 m/s (30 Pa on 120 MPa
     pre-stress, stays at friction equilibrium).
   - **T2a (asymmetric v_x)**: Q[VX] = +V_test on tet 0, Q[VX] = 0 on tet 1. Local:
     `v_t2⁺ = 0`, `v_t2⁻ = +V_test`, `v_t2⁻ - v_t2⁺ = +V_test`. Pelties 7c:
     `tau_2_trial = eta_s · (-V_test) = -0.5·Zs·V_test`. Sign opposite to T2 of the
     symmetric case, confirming the sign map.
   - **T3a (asymmetric v_z)**: Q[VZ] = +V_test on tet 0, Q[VZ] = 0 on tet 1. Expected:
     `tau_1_trial = eta_s · (+V_test) = +0.5·Zs·V_test`.
   - Tolerance: 1% of expected perturbation (matches the existing test's convention).

2. **Interior stress-jump case** (T4a):
   - Q[SXY] = +A on tet 0, Q[SXY] = 0 on tet 1 (A = 1 kPa). Pelties 7b:
     `tau_1_trial` uses `tau_1⁺/Zs⁺ + tau_1⁻/Zs⁻`, but in fault-local coordinates
     `tau_1 = σ_{n,t1} = -σ_yz` (since can_n = -ŷ, can_t1 = -ẑ: n⊗t1 = (-ŷ)⊗(-ẑ) = ŷẑ, so
     σ_{local,SXY} = +σ_yz), while `tau_2 = σ_{n,t2} = -σ_yx = -σ_xy`. So Q[SXY] (global)
     → local `tau_2` (not tau_1). With A on − side only:
     `tau_2_trial = eta_s · (-A/Zs_plus + 0/Zs_minus)` ... **this gets tricky** — the
     stress-jump formula requires careful unpacking. Implementation note: compute the
     expected value by running the same rotation (`Tinv_can`) on both sides in the test
     harness, then calling `FaultFaceFlux::ComputeTrialTraction` directly with those
     rotated inputs. Compare the `Mult()`-driven `data.tau2_corr - tau2_0` against that
     reference. Tolerance 1%. This form — "test harness reproduces the whole pipeline
     except the Mult-driven interpolation" — isolates the per-side Q interpolation step
     from the rotation + trial-traction step.

3. **Shared-fault path** (`tests/parallel/test_shared_fault_flux_path.cpp`, np=2):
   - Same 2-tet fixture partitioned one-tet-per-rank.
   - Run the same T1a–T4a patterns. Each rank reads its local DOFData after `Mult`; since
     the fault QP lives on exactly one DOF in the interior and appears twice under shared
     (once per rank, independent copies), assert that BOTH ranks' DOFData.{tau1_corr,
     tau2_corr, sigma_n_corr} match the analytic target AND match each other bit-exactly
     (relative 1e-12). This is the operational probe for H3 under the shared branch.

### Interfaces
- No new public interfaces. Reuses `FaultFaceFlux::ComputeTrialTraction`,
  `GodunovFlux::BuildRotationInverse`, `WaveOperator::Mult`, `WaveOperator::SetFaultDOFData`.

### Edge Cases to Handle
- **Tolerance tightness vs. friction-solver drift**: the existing T1–T3 use a 10% tolerance
  on top of expected `Zp*V_test` due to Brent-solver corrections at V ≈ V_ini. Use the same
  convention here. If T1a–T3a require tighter than 5%, it is a signal that the interior
  branch has extra coupling beyond Pelties 7 — FAIL and report the discrepancy.
- **np=2 partitioner non-determinism**: reuse the explicit `partitioning[]={0,1}` pattern
  from Phase 1. Skip (not fail) if the partitioning degenerates.

### Acceptance Criteria
- [ ] Interior asymmetric T1a passes (check σ_n coupling).
- [ ] Interior asymmetric T2a passes (strike, the dominant TPV102 load direction).
- [ ] Interior asymmetric T3a passes (dip).
- [ ] Interior asymmetric T4a passes (stress jump).
- [ ] Shared-fault np=2 T1a–T4a all pass AND both ranks agree bit-for-bit.
- [ ] Tampering: if T2a FAILS on interior but T2 (symmetric) PASSES on interior, the
      pepper source is indeed in the per-side Q interpolation path.

### Dependencies
- Depends on: Phase 1 (for `GetFaultFaceDofOffset` to confirm DOF placement before Mult).
- Required by: nothing.

---

## Phase 4: Shared-face DOF pairing (H4)

### Goal
After this phase: two new tests verify the `MakeFaceKey` integer-exact pairing invariant
and the `fault_face_dof_offset_` / `shared_fault_dof_offset_` non-overlap invariant, both
at the unit level on a 4-tet fixture across np=2.

### Files to Create
- `tests/parallel/test_make_face_key_uniqueness.cpp` — np=2.
- `tests/unit/test_fault_dof_offset_invariants.cpp` — serial.

### Files to Modify
- `dynamic/shared_fault_key.hpp` — no changes; `MakeFaceKey` is already public.
- `Makefile` — 2 link rules + test targets.

### Detailed Requirements

1. **MakeFaceKey uniqueness** (np=2):
   - Fixture: 4-tet mesh with two interior fault faces at y=0, chosen so both faces span
     the rank boundary. (Vertices: V0=(0,-1,0), V1=(0,0,0), V2=(1,0,0), V3=(0,0,1),
     V4=(0,+1,0), V5=(2,0,0), V6=(0,0,2); 4 tets wrapping two shared fault faces.)
   - On each rank, enumerate `fault_shared_faces_`, compute `MakeFaceKey(local_face, gvi,
     pmesh)` for each via existing `wave_operator.inl:2733` logic.
   - MPI_Allgather all keys.
   - Assert: each `{key.v[0], key.v[1], key.v[2]}` triple appears on **exactly 2 ranks**
     (once from rank A, once from rank B) and the matching pair corresponds to the same
     physical face (verified by centroid coordinates agreeing to 1e-10 m).
2. **DOF offset invariants** (serial, 2-tet fixture from Phase 1):
   - After `WaveOperator` ctor + `SetFaultDOFData`, read
     `fault_face_dof_offset_` and `shared_fault_dof_offset_` (via getters from Phase 1).
   - Assert: the set `{dof_offset + q : (face_idx, dof_offset) ∈ fault_face_dof_offset_,
     q ∈ [0, nqp_per_face)}` is disjoint from the analogous set for shared-face offsets.
   - Assert: the union is a contiguous range `[0, fault_dof_data_->size())`.
3. **Repeat serial test on 4-tet fixture with two fault faces** to confirm the offset
   invariant holds for multi-face cases.

### Interfaces
- No new public interfaces. Phase-1 getters are sufficient.

### Edge Cases to Handle
- **Zero shared fault faces** (np=1 build): np=2 test reports "skipped: no shared faces
  produced by partitioner" on np=1. Exit 0.
- **MakeFaceKey hash collision across distinct faces**: cannot occur in principle — the key
  is the sorted triple of global vertex indices, so distinct triangular faces have distinct
  keys. Test asserts uniqueness within each rank as a sanity check.

### Acceptance Criteria
- [ ] MakeFaceKey test on np=2 produces a perfect 2-per-key histogram.
- [ ] DOF-offset test produces 0 overlaps and a contiguous union.
- [ ] Tampering: corrupt one entry of `fault_face_dof_offset_` via env var
      `SEAS_TEST_TAMPER=1` — test MUST FAIL with offset-overlap message.

### Dependencies
- Depends on: Phase 1 (getters).
- Required by: nothing.

---

## Phase 5: Bulk solver pattern recovery (H5)

### Goal
After this phase: two new tests on a fault-less cubic fixture probe whether the volume +
Godunov-interior-flux + RK4 stack can reproduce an analytic plane wave, and whether the
stack is stable under per-DOF O(1e-10 × τ_ini) noise for 100 RK4 steps.

### Files to Create
- `tests/unit/test_bulk_plane_wave_recovery.cpp` — serial, 2×2×2 hex.
- `tests/unit/test_bulk_rk4_noise_stability.cpp` — serial, 2×2×2 hex.

### Files to Modify
- `Makefile` — 2 link rules + test targets.

### Detailed Requirements

1. **Plane-wave recovery** (T5-PW):
   - Fixture: 2×2×2 hex mesh, [0,1]³, absorbing BCs (attr 5) on all external faces. No fault.
   - Initial condition: a P-wave propagating along +x, analytic solution
     `Q[SXX](x,0) = A · sin(k·x)`, `Q[VX](x,0) = -(A/(ρ·cp)) · sin(k·x)`, others zero.
     k = 2π (one wavelength per box). A = 1 MPa.
   - Advance one RK4 step with dt = 0.1 · (1/cp) · (1/2) [CFL safe].
   - Analytic: after dt, `Q[SXX](x, dt) = A · sin(k·(x - cp·dt))`.
   - Assertion: max-DOF error in `Q[SXX]` at t=dt is bounded by
     `(cp·dt·k·A)³/6 · ndof` (pure Taylor truncation of the analytic solution); this
     is the RK4 local truncation error. Tolerance ≈ 1e-4 × A for this fixture. If violated,
     the volume + Godunov path is misconstructed.
2. **RK4 noise stability** (T5-N):
   - Fixture: same 2×2×2 hex, absorbing BCs.
   - Initial condition: `Q` = 0 + uniform random per-DOF perturbation of amplitude
     1e-10 × τ_ini (seed fixed). All 9 components perturbed.
   - Advance 100 RK4 steps with dt = 0.1 · (h/cp) where h = 0.5 (one hex).
   - Assertion: `max|Q|_final ≤ 10 × max|Q|_initial`. Absorbing BCs should radiate noise
     out, so the final amplitude should actually DECREASE. A growth factor > 10 indicates
     the bulk RK4 stack is amplifying per-DOF noise — that IS the pepper signature at the
     bulk level. If this test fails, H5 is the primary suspect, not H1–H4.
   - Secondary assertion: per-element variance of `Q[SYY]` at t=final does not exceed the
     per-element variance of `Q[SXX]` by more than 10%. Symmetric fixture + symmetric IC
     should yield symmetric results; directional bias in the bulk noise = directional
     instability.

### Interfaces
- No new interfaces.

### Edge Cases to Handle
- **CFL violation**: dt chosen conservatively (factor 0.1 below the nominal cp·Δt/h = 1
  limit). Test aborts with a CFL message if dt > 0.5·h/cp.
- **Random seed determinism**: use `std::mt19937` with seed 42; never `rand()`. Asserts are
  deterministic across runs on the same machine.

### Acceptance Criteria
- [ ] T5-PW: max-DOF error < 1e-4 × A for order=1 DG, stricter at order=2 (order-2 is
      optional and kept as a secondary case).
- [ ] T5-N: max|Q|_final / max|Q|_initial ≤ 10 AND per-element SYY/SXX variance ratio
      within 10%.
- [ ] Both tests run in < 2 s.

### Dependencies
- Depends on: nothing.
- Required by: nothing.

---

## Phase 6: Aggregator target (`test-tpv102-pepper`)

### Goal
After this phase: a single `make test-tpv102-pepper` runs every new test in phases 1–5,
plus the five pre-existing tests that already constrain H1–H5, and reports a combined
PASS/FAIL summary. Gives the user a one-command check after any `dynamic/` change.

### Files to Modify
- `Makefile` — add:
  ```
  test-tpv102-pepper: \
      test-interior-vs-shared-branch-live \
      test-sign-flipped-elem1-on-plus-truth-table \
      test-interior-fault-flux-path \
      test-interior-fault-flux-path-asymmetric \
      test-shared-fault-flux-path \
      test-canonical-rotation-pure-strikeslip \
      test-per-qp-vs-centroid-basis-planar-1el \
      test-no-penalty-dynamic-rupture \
      test-volume-jacobian-single-channel \
      test-make-face-key-uniqueness \
      test-fault-dof-offset-invariants \
      test-bulk-plane-wave-recovery \
      test-bulk-rk4-noise-stability \
      test-shared-fault-role-consistency \
      test-shared-fault-dof-data-consistency
  	@echo "=========================================="
  	@echo "  TPV102 pepper-bug test suite: COMPLETE"
  	@echo "=========================================="
  ```
- Document at the top of `PLAN.md`: "Run `make test-tpv102-pepper` after any change to
  `dynamic/wave_operator.{hpp,inl}` or `dynamic/fault_face_flux.{hpp,cpp}`."

### Acceptance Criteria
- [ ] `make test-tpv102-pepper` sequences all 15 tests.
- [ ] Total wall-time ≤ 60 s on a laptop.
- [ ] Any single test failure stops the suite and exits non-zero.

### Dependencies
- Depends on: phases 1–5 complete.

---

## Testing Strategy

- **Unit-level only**: no test invokes the TPV102 driver, full mesh, or Frontera job.
- **Tampering demonstrations**: each phase includes a `SEAS_TEST_TAMPER` env-var path that
  introduces a known wrong value; the test MUST fail. This replaces tautological checks.
- **Diagnostic output**: on failure, print (test name, failing assertion line, expected,
  actual, |difference|) — so a CI log is self-contained.
- **Numerical tolerances**: bit-exact (tol = 0) where the algebra is signed-permutation
  (Phase 2); 1e-12 relative for floating-point rotation compositions (Phases 1, 3); 1%
  absolute for `Mult`-driven Pelties-7 checks (Phase 3); RK4-truncation-based for Phase 5.
- **Coverage check**: after all phases are in place, run `gcov` / `llvm-cov` against the
  hotspots `wave_operator.inl:1055-1188, 1608-1665, 258-320, 2680-2775, 483-499, 831-925`
  and `fault_face_flux.cpp:111-170`. Each line in these ranges should be exercised by at
  least one of the new tests (target: ≥ 90% line coverage of the hotspots).

---

## Risk Assessment

- **Risk 1 — Internal-state exposure**: Phase 1's getters under `#ifdef SEAS_TEST_INTERNAL`
  add a test-only surface. If the project accidentally builds production code with
  `-DSEAS_TEST_INTERNAL`, the getters become accessible in production. Mitigation: the
  Makefile adds this flag only via per-target overrides on the new test link rules; no
  change to the top-level `CXXFLAGS`. A grep guard in CI (`grep -r SEAS_TEST_INTERNAL
  drivers/` returns 0 hits) closes the loop.
- **Risk 2 — Partitioner dependence**: Phase 1 and Phase 3 rely on np=2 putting one tet per
  rank. Mitigation: explicit `partitioning[]` array passed to `ParMesh`. The test aborts
  with a fixture-bug message if the partitioner overrides this (it should not, but the
  assertion makes the test robust).
- **Risk 3 — Test 5 false positive**: RK4 noise stability can fail for physical reasons
  (e.g., absorbing BCs insufficient to damp). Mitigation: the 10× growth factor is
  generous; if it fires, it's a real signal, not implementation drift. Cross-check against
  the existing `test_absorbing_bc_energy_decay`.
- **Risk 4 — Phase 3 stress-jump algebra error**: T4a's expected value depends on a
  careful bookkeeping of Tinv_can sign conventions. Mitigation: the test harness recomputes
  the expected value by calling `FaultFaceFlux::ComputeTrialTraction` on the rotated
  states directly; the `Mult`-driven value is compared against this reference, not against
  a paper-derived formula. This removes the author-algebra risk.
- **Risk 5 — Makefile maintenance burden**: 10+ new test targets. Mitigation: the
  Phase-6 aggregator gives a single command; individual targets are only used for
  debugging specific failures.
- **Known tricky area — DOFData identity across the two branches**: `test_shared_fault_dof_data_consistency`
  Phase D already handles this at the integration level. The new Phase 3 shared test
  duplicates some coverage but at a smaller fixture — intentional, to isolate failure
  modes.

---

## Appendix: Mapping back to the pepper report hotspots

| Hotspot in pepper report                                    | Covered by            |
|-------------------------------------------------------------|-----------------------|
| `dynamic/wave_operator.inl:1055-1188` (interior branch)     | Phase 1, Phase 3 T1a–T4a, Phase 2 interior |
| `dynamic/wave_operator.inl:1608-1665` (shared branch)       | Phase 1, Phase 3 shared, Phase 2 shared   |
| `dynamic/wave_operator.inl:1092-1134` (can_*, elem1_on_plus, Tinv) | Phase 2, Phase 3, existing `test_canonical_rotation_pure_strikeslip` |
| `dynamic/wave_operator.inl:887-913` (Q_self interpolation)  | Phase 3 T4a, existing `test_wave_operator_spatial_derivative` |
| `dynamic/wave_operator.inl:258-320` (fault face lists + basis idx) | Phase 4 P4-T2           |
| `dynamic/wave_operator.inl:2680-2775` (MakeFaceKey + shared output) | Phase 4 P4-T1     |
| `dynamic/wave_operator.inl:483-499, 831-925` (Mult + bulk)  | Phase 5, existing `test_no_penalty_dynamic_rupture`, `test_volume_jacobian_single_channel` |
| `dynamic/fault_face_flux.cpp:111-170` (Evaluate)            | existing `test_fault_face_flux.cpp`, strengthened via Phase 3 call-path tests |

Every hotspot is hit by at least one test from this plan plus existing coverage; H1–H5
can all be falsified locally.
