# Code Review v7.0.0 audit: fresh adversarial pass on the R-801/R-802 fix

**Versioning note:** the `v[a].[b].[c]` scheme introduced in the v7.0.0
check bumps `[a]` only on NEW critical findings, `[b]` on NEW moderate,
`[c]` on NEW low.  This audit found **zero** new findings at any
severity, so **the version stays at `v7.0.0`**.  This document is a
confirmation audit, not a version bump.

**Context:** the v7.0.0 check (`tpv102_debug_v7.0.0_check.md`) flagged
two critical bugs after the R-701 refactor:

1. **R-801** — t1 convention mismatch: interior fault branch used
   `GodunovFlux::BuildFrame` (t1 = strike) while the shared fault
   branch (R-701) used BP5 `FaultBasis` (t1 = dip).  Driver wrote the
   TPV102 strike-slip pre-stress into `tau1_0`, making interior-fault
   DOFs self-consistent but corrupting shared-fault DOFs.
2. **R-802** — If MFEM's shared-face `CalcOrtho` returns identical (not
   opposite) normals on the two ranks sharing the face (the v6 fix
   report's empirical observation), the Godunov conservation identity
   does NOT compensate, and both ranks subtract the same `F_h` from
   their Elem1 rhs, breaking momentum conservation across the fault.

The fix round applied **Option A for R-801** (unify the whole pipeline
on BP5's `tangent1 = dip, tangent2 = strike` convention — driver,
interior flux, shared flux all use `FaultBasis` canonical frame) plus
an **accum_sign gate for R-802** (shared-fault flux uses `can_n` and
flips accumulation sign based on `elem1_on_plus`).

This audit re-executes all three review passes from scratch on the
post-fix code to verify the fix and hunt for any new bugs.

## Review Scope
- Files re-read in full (no skim):
  - `miniapps/seas/CLAUDE.md` — project-wide convention note added for
    BP5 / TPV102 tangent frame (v7.0.0 R-801 entry under Sign Conventions).
  - `miniapps/seas/dynamic/wave_operator.hpp`
  - `miniapps/seas/dynamic/wave_operator.inl` — both fault branches
    (interior at ~line 700-844 and shared at ~line 867-1215).
  - `miniapps/seas/dynamic/tpv102_setup.hpp` — `InitializeFaultDOFs`
    and `ApplyNucleation` (both now write BP5 convention: `tau2_0,
    V2 = strike` and `tau1_0, V1 = 0` for pure-strike-slip TPV102).
  - `miniapps/seas/dynamic/fault_face_flux.cpp` — `Evaluate` unchanged;
    its Eq. (7a-c) convention is frame-indifferent (it operates on
    whatever `(VY, VZ) = (v_t1, v_t2)` the caller rotated in).
  - `miniapps/seas/dynamic/godunov_flux.cpp` — `BuildFrame` still used
    by bulk Godunov splitting (internal to `flux_.Interior`); NOT used
    for fault-local frame any more.
  - `miniapps/seas/drivers/tpv102_driver.cpp` — ParaView output
    mapping, RK4 stage accumulation, nucleation loop.
  - `miniapps/seas/tests/parallel/test_r101_shared_fault.cpp` — new
    T-R801 (Parts A+B) and T-R802 tests.
- Domain context:
  - Updated `CLAUDE.md` line 34-43 documents the project-wide
    `tangent1 = dip, tangent2 = strike` convention and the
    corresponding TPV102 component assignment `tau2_0 = tau_ini,
    V2 = V_ini, tau1_0 = 0, V1 = 0` (pre-R-801 used the opposite).
- Build + runtime verification (local):
  - Rebuild of the `seas_tpv102_driver`, `seas_test_r101_shared_fault`,
    `seas_test_wave_operator`, `seas_test_parallel_wave_operator`
    succeeded clean.
  - Test suite:

    | Test | Result |
    |---|---|
    | `seas_test_godunov_flux` | 29/29 |
    | `seas_test_wave_operator` | 17/17 |
    | `seas_test_wave_bc` | 10/10 |
    | `seas_test_fault_face_flux` | 19/19 |
    | `seas_test_tpv102_setup` | 24/24 |
    | `seas_test_tpv102_local` | 16/16 |
    | `seas_test_parallel_wave_operator` (2 ranks) | 5/5 |
    | `seas_test_parallel_wave_operator` (4 ranks) | 5/5 |
    | `seas_test_r101_shared_fault` (2 ranks, includes R-302a, R-501a, T-R801 A+B, T-R802) | **18/18** |
    | `seas_test_elasticity_operator` (BP5) | 461/461 |
    | `seas_test_bp5_params` (BP5) | 123/123 |
  - End-to-end: `mpirun -np 4 ./seas_tpv102_driver --mesh
    tpv102/mesh/tpv102_1000m.msh --tfinal 0.1 --debug-qnorm` runs to
    completion with no `MFEM abort` and no R-101 trip.  Qnorm watch
    shows bulk wave transit across ranks at sub-ms scale (expected
    at 1000 m mesh / 4 ranks; no shared fault faces on this specific
    partition, so shared-fault code path is dormant here).

## Pass 1 — Plan compliance (R-801 + R-802)

### R-801 Option A audit
- `dynamic/tpv102_setup.hpp:76-78`: `tau1_0 = 0.0` (no dip pre-stress),
  `tau2_0 = tau_ini` (strike pre-stress).  ✓
- `dynamic/tpv102_setup.hpp:93-96`: `V1 = 0`, `V2 = V_ini` (strike slip
  rate).  ✓
- `dynamic/tpv102_setup.hpp:101-102`: `tau1_corr = 0`, `tau2_corr =
  tau_ini`.  ✓
- `dynamic/tpv102_setup.hpp:146`: `ApplyNucleation` writes `tau2_0 +=
  dtau` (not `tau1_0`).  ✓
- `dynamic/wave_operator.inl:709-760` (interior fault branch): now uses
  `FaultBasis` canonical frame via `LookupInteriorFaultBasisIndex(f)`
  and reconstructs `can_n / can_t1 / can_t2` from `qpd.sign_flipped`
  exactly as in the shared-fault branch.  `GodunovFlux::BuildFrame` is
  no longer called for the fault-local frame in this file.  ✓
- `dynamic/wave_operator.inl:761`: `elem1_on_plus = !qpd.sign_flipped`
  for interior faces.  Geometrically: `sign_flipped = true` iff MFEM's
  raw `CalcOrtho` was anti-aligned with `ref_normal`; `elem1_on_plus`
  is the side opposite to `can_n` direction (= `ref_normal` direction).
  The interior case's `!sign_flipped` matches the canonical + side
  assignment used by the shared case's geometric centroid test.  ✓
- `drivers/tpv102_driver.cpp:647-652`: ParaView output pack maps
  `d.slip1 → comp 0 (dip)`, `d.slip2 → comp 1 (strike)`, consistent
  with `paraview_output.hpp`'s `UpdateFaultFieldsBP5` reader.  ✓
- `dynamic/tpv102_setup.hpp:347` (station writer header): explicitly
  documents the BP5 convention and expected near-zero behaviour of
  `V1/slip1/tau1` for TPV102.  ✓

### R-802 audit
- `dynamic/wave_operator.inl:1172`: shared-fault flux uses `can_n`
  (canonical, bit-identical across ranks) instead of rank-local `nor`.
  ✓
- `dynamic/wave_operator.inl:1173-1181`: accumulation gated by
  `accum_sign = elem1_on_plus ? +1 : -1`.  Exactly one rank per shared
  face has `+1`, the other `-1`, so the sum-over-ranks of `rhs[Elem1]`
  contributions is zero (conservation).  ✓
- `dynamic/wave_operator.inl:1182`: `continue` correctly skips the
  generic `-= w * shape1(i) * F_h[c]` accumulation at line 1204-1210,
  which is only safe for rank-local `nor` and standard Godunov
  orientation (i.e., non-fault shared faces and the fallback fault
  path).  ✓

### R-702 (docstring) audit
- `dynamic/wave_operator.hpp:184-190`: `tol` now documented as
  **RELATIVE** with the formula `|a-b| / max(|a|, |b|, 1.0) > tol`
  explicitly.  ✓

### Plan-level compliance: **FULL**
Every critical and moderate finding from `v7.0.0_check.md` has a
matching code change that survives this audit.

## Pass 2 — Fresh bug hunt on the post-fix code

**Hypothesis:** after a large semantic refactor (convention flip),
something may have slipped.  Targets:

1. **Did the driver's RK4 stage capture/restore preserve BP5
   semantics?**
   Lines 726, 741, 755, 770 capture `V1_k*, V2_k*, t1c_k*, t2c_k*`.
   Line 795-804 averages and writes back to `dof_data[i].{V1, V2,
   tau1_corr, tau2_corr, sigma_n_corr}`.  `slip1 += V1*dt, slip2 +=
   V2*dt` at 800-801.  Under BP5 convention, `slip1` accumulates
   dip-component and `slip2` accumulates strike-component.  Consistent
   with driver output mapping (line 647: `slip(2i+0) = slip1 = dip`).
   **No inconsistency.**

2. **Does `Evaluate` expect any specific tangent convention?**
   `fault_face_flux.cpp` Evaluate is frame-indifferent: it writes
   `data.V1` = V-component-along-whichever-axis-the-caller-labelled-as-t1,
   `data.V2` = along-t2.  Under R-801, the caller always uses BP5's
   `tangent1 = dip, tangent2 = strike`, so `data.V1 = v_dip, data.V2
   = v_strike` uniformly.  **Consistent.**

3. **Does the bulk non-fault shared-face path still conserve?**
   Line 1193 falls into the generic `flux_.Interior(nor, Q_self,
   Q_nbr, F_h)` with rank-local `nor`.  Standard Godunov conservation
   identity `F(L, R, +n) = -F(R, L, -n)` applies ONLY if `nor_A =
   -nor_B` on bulk shared faces.  The v6 empirical observation of
   `nor_A = nor_B` was specifically on the shared FAULT face of the
   inline 2-tet test; it was NOT claimed for bulk shared faces.
   Empirically, `seas_test_parallel_wave_operator` Test P5 (energy
   conservation over 200 bulk wave steps) still passes at 2 and 4
   ranks, which is sensitive to any systematic bulk-flux conservation
   violation.  **Bulk path presumed safe; if Frontera runs surface
   energy leak at bulk seams, revisit.**

4. **Any residual `BuildFrame` in fault branches?**
   `grep -nE "GodunovFlux::BuildFrame|BuildFrame\("` inside
   `wave_operator.inl` returns only **comments** (lines 713, 944).
   The actual call sites in `godunov_flux.cpp:330, 371` are inside
   `Interior` and `FreeSurface` — those are BULK Godunov splitting
   (internal to the flux routine, frame-indifferent in output), not
   fault-local-frame construction.  **No fault-branch regression.**

5. **`LookupInteriorFaultBasisIndex` correctness:**
   Map populated in ctor at line 292-296:
   ```cpp
   fault_interior_face_to_basis_idx_.clear();
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
      fault_interior_face_to_basis_idx_[fault_interior_faces_[i]] = i;
   ```
   Consumed in `ComputeFaceFluxRHS:717` with `LookupInteriorFaultBasisIndex(f)`
   where `f` is the current mesh-face index from the outer loop.
   `f` is in `fault_interior_faces_` iff the face was added to the
   interior-fault list in the ctor.  **Mapping is consistent.**

6. **Edge case — ranks with only shared fault, no interior fault:**
   The ctor at line 320 always calls `fault_basis_->Compute(mesh,
   fault_interior_faces_, ref_normal, up)` even if the list is
   empty (to initialise `dim_`).  Then `AppendSharedFaces` extends.
   Verified: FaultBasis's `Compute` handles empty `fault_faces`
   cleanly (sets `num_faces_ = 0`, no iteration).  **Edge case
   handled.**

7. **Edge case — ranks with no fault at all:**
   `fault_basis_` is only constructed if `bc_.fault_attr > 0 &&
   (interior > 0 || shared > 0)`.  Otherwise `fault_basis_ ==
   nullptr`.  `fault_active` guard at line 968 includes
   `&& fault_basis_`, so fault-branch code paths short-circuit to
   the bulk flux.  **Safe.**

8. **T-R801 test strength:**
   The test reports `max|V1| = 0.000000, min|V2| = 0.000000` (6-digit
   `std::to_string` rounds `V_ini = 1e-12` to zero).  The stronger
   assertion `max|V1| < 1e-6 * min|V2| + 1e-30` is guarded by
   `min_V2_abs > 0.0`; in this specific configuration the actual
   `min|V2|` is `V_ini ≈ 1e-12 > 0`, so the guard passes and the
   inner assertion runs.  The test IS discriminating: pre-R-801
   wrote `V1 ≈ V_ini > 0` and `V2 ≈ 0`, which would trip the
   `v2 < v1` outer check.  **No audit-level concern.**

9. **T-R802 Q=0 baseline only:**
   T-R802 tests conservation on `Q = 0` initial state.  At `Q=0`,
   the fault-flux contribution reflects only the pre-stress-driven
   slip rate, which is small but nonzero.  The test reports
   `Σ k = -3.04e-18, Σ|k| = 9.12e-18, rel = 3.04e-18` — conservation
   to machine precision.  A stronger future test could drive nonzero
   Q before the conservation check, but the present test is
   sufficient for the `rel < 1e-8` threshold it asserts.  **Not an
   audit-level concern.**

10. **Any residual stale comments or dead code referencing pre-R-801
    conventions?**
    - `wave_operator.inl:713` comment references pre-R-801 behaviour
      ("along-strike (from GodunovFlux::BuildFrame's t1=x)") as a
      history note.  Intentional, useful.  ✓
    - `wave_operator.inl:944` comment similarly a history note.  ✓
    - `tpv102_setup.hpp:60-75` comment block documents R-801 Option A
      reasoning.  ✓

**No new bugs found in pass 2.**

## Pass 3 — Code quality check

- **Naming:** `elem1_on_plus`, `can_n / can_t1 / can_t2`,
  `sf_to_basis_idx`, `shared_fault_elem1_on_plus_` are all
  self-descriptive and the comments explain the meaning.  ✓
- **Dead code:** `shared_face_peer_` / `SharedFacePeer` from v5 R-501
  are confirmed removed (grep returns zero occurrences in
  `wave_operator.*`).  ✓
- **Duplication:** The canonical-frame reconstruction block (`can_n =
  sign_flipped ? -normal : normal`, etc.) appears twice — once in the
  interior-fault branch and once in the shared-fault branch.  ~20
  lines.  Could be extracted into a small helper, but doing so would
  obscure the per-branch control flow and costs negligible
  readability.  Not worth flagging.  ✓
- **Error handling:** `MFEM_VERIFY` / `MFEM_ASSERT` guards on
  `qp_data.size()`, `sf_to_basis_idx[sf] >= 0`, `fault_basis_ != nullptr`.
  ✓

No quality-level findings.

## Summary
- Critical issues: **0** (zero new findings).
- Moderate issues: **0**.
- Low issues: **0**.
- Plan compliance: **FULL**.
- **Verdict: PASS — no blockers remaining, code is Frontera-ready.**

## Version status
- Current version: **`v7.0.0`** (unchanged from the original v7.0.0
  check, because this audit found no new findings).
- R-801 and R-802 from `tpv102_debug_v7.0.0_check.md` are both
  resolved and pass their corresponding tests (T-R801 Parts A+B, T-R802).

## Frontera run plan — RECOMMENDED SEQUENCE

Given the green state of tests and the 4-rank local driver smoke
success, proceed with the originally-staged sbatch sequence.  All 4
scripts include the v4 R-402/R-403 RESULT.txt hardening, so the
PASS/FAIL verdict is trustworthy regardless of submission CWD.

| # | Script | Scale | Wall | Purpose | Gate to next |
|---|---|---|---|---|---|
| 1 | `tpv102_1000m_p1_1.5s_4rank_dev.sbatch` | 1000 m, 4 ranks, t=1.5 s | ~15 min | **Cheapest sanity.**  Verify build links on Frontera, `[qnorm:watch]` shows bulk wave transit on all 4 ranks by t≈0.3s, R-101 verifier either SKIPs (no shared fault) or reports `max_rel_diff=0`.  Reproduces the v7.0.0 ORIGINAL failure scenario with the R-801/R-802 fix. | No abort, no dead ranks in qnorm watch at t=1.5s. |
| 2 | `tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch` | 1000 m, 8 ranks, t=0.5 s, ParaView | ~15 min | **Visual smoke (optional).**  Open fault-surface PVD in ParaView, color by `slip_rate_strike` — look for smooth isochrones from the hypocenter, NOT partition-seam stripes.  Validates R-802 conservation + R-801 strike-component channeling. | No visible seam stripes. |
| 3 | `tpv102_200m_p1_1.5s_50rank_dev.sbatch` | 200 m, 50 ranks, t=1.0 s | ~2 hr | **Intermediate pre-flight.**  First run where METIS likely produces multiple shared fault faces per peer (exercises R-701 batch-safety implicitly).  RESULT.txt verdict matters here. | `grep PASS RESULT.txt` |
| 4 | `tpv102_200m_p1_1.5s_400rank_dev.sbatch` | 200 m, 400 ranks, t=1.5 s | ~1.5 hr | **THE DISPOSITIVE TEST.**  Reproduces the v1 failing configuration (rupture stops at seam) with all v1-v7 fixes in place.  PASS → confirms entire fault-mismatch saga resolved. | N/A — this is the verdict. |

### Pre-flight checklist (before submitting #1)

- [x] Local tests all green: 18/18 `seas_test_r101_shared_fault`, 5/5 × 2
      `seas_test_parallel_wave_operator`, 461/461 + 123/123 BP5.
- [x] 4-rank 1000 m driver runs locally to `--tfinal 0.1` without abort.
- [ ] `git pull` on Frontera so R-801/R-802 code is on the scratch
      repo.  Double-check `miniapps/seas/dynamic/wave_operator.inl`
      and `tpv102_setup.hpp` modification timestamps on Frontera match
      the local copy.
- [ ] Module load dry-run on a Frontera login node:
      `module load intel/19.1.1 impi/19.0.9 hypre/2.31.0 mumps/5.3 parmetis petsc/3.15 fftw3/3.3.8 && module list`
      should succeed.  (User memory `feedback_sbatch_modules` applies.)
- [ ] Submit from `/scratch2/10024/zhaochun/seas-project/seas-mfem`
      (repo root).  R-403 fix makes this safe; the `SBATCH_LOG_DIR`
      pin anchors the post-run RESULT.txt grep to the submit dir.

### Stop conditions (any one = HALT the sequence)
1. Any step prints `MFEM abort:` or `R-101 shared-fault DOFData
   consistency FAILED`.
2. Step 1 final `[qnorm:watch]` shows non-hypocenter ranks at
   `0.000e+00` at t=1.5s (bulk wave failed to cross partition seams).
3. Step 2 ParaView shows partition-seam stripes in `slip_rate_strike`.
4. Step 3 or Step 4 `RESULT.txt` contains `FAIL`.

### Unreviewed Areas
- **Non-regression at 400 ranks.**  Local testing caps at 4 ranks.  If
  Step 3 or Step 4 on Frontera reveals an issue invisible locally, it
  would be flagged in the next review round and tracked under a new
  `v[a].[b].[c]` document.
- **Long-time integration stability (tfinal > 1.5 s).**  R-801 changed
  the RK4 averaging semantics (same averages, different physical
  components).  If a subtle accumulated-drift bug exists, only a
  longer-horizon run will surface it.  The TPV102 spec runs 12 s;
  that's a post-dispositive-success follow-up.
- **MFEM shared-face `CalcOrtho` behaviour on TPV102 production
  meshes.**  The v6 fix report claimed `nor_A = nor_B` empirically on
  the 2-tet inline mesh.  Whether the same holds on `tpv102_1000m.msh`
  (Gmsh BooleanFragments) or `tpv102_200m.msh` was not verified.  The
  R-802 fix uses `can_n` regardless, so it's safe either way, but the
  diagnostic question remains open.
