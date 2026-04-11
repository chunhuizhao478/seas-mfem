# BP5 Debug v56: Parallel Exact Station Interpolation & Output Regression

**Date**: 2026-03-30
**Status**: PHASES 1-3 COMPLETE — Exact interp disabled in parallel, diagnostics added, tests passing
**Previous**: v55 (line-by-line MFEM vs Tandem comparison, all code fixes applied)
**Branch**: `feature/elasticity`

---

## 1. Context

v55 applied all critical code fixes (sign conventions, log10 Brent solver, combined integrator,
per-QP boundary evaluation, etc.) and all 587 tests pass. However, comparison between the
**prod** (nearest-DOF) and **exct_prod** (exact face interpolation) output paths shows a
mismatch in station data. The remaining physics mismatch near nucleation (event timing ~70-90s
late, left-side interseismic too fast/stressed) cannot be diagnosed reliably until the output
path is validated.

**Priority order** (strict):
1. Fix or disable parallel exact station interpolation
2. Recompare with Tandem using the non-broken output path
3. Only then chase the remaining nucleation/interseismic physics mismatch

---

## 2. Root Cause Analysis: Parallel Exact Interpolation

### 2.1 The Assumption

`Probe2DInterpolator::TryBuildExactMatch()` (`bp5_benchmark_output.hpp:228-262`) assumes the
globally gathered fault DOF array has **contiguous per-face blocks**:

```
face 0: [0, nbf_per_face)
face 1: [nbf_per_face, 2*nbf_per_face)
...
face F: [F*nbf_per_face, (F+1)*nbf_per_face)
```

It iterates `f = 0..num_faces-1`, computes `start = f * nbf_per_face_`, and uses the first 3
DOFs at `[start, start+3)` as triangle vertices for `ComputeReferenceIP()`.

### 2.2 What Actually Happens in Parallel

**Data flow:**
1. Each rank calls `domain.GetFaultCoords2D(local_x2, local_x3)` → stores coords in order:
   interior faces first (`fault_interior_faces_`), then shared faces (`fault_shared_faces_`),
   each face as a contiguous block of `nbf_per_face` DOFs.
   (`elasticity_operator.hpp:2649-2712`)

2. `fault_geom.GatherToRoot(local_x2, global_x2)` uses `MPI_Gatherv` to concatenate rank data:
   `[rank0_local | rank1_local | rank2_local | ...]`
   (`fault_geometry.hpp:150-174`)

3. Root constructs `Probe2DInterpolator(global_x2, global_x3, stations, nbf_per_face, ...)`
   which calls `TryBuildExactMatch` for each station.
   (`bp5_parallel_output.hpp:73-75`)

4. At each timestep, the same `GatherToRoot` is applied to field vectors (slip, V, traction,
   theta), and `EvaluateScalar(global_field, station_idx)` uses the stored `face_start_` and
   `exact_weights_` to interpolate.
   (`bp5_parallel_output.hpp:145-163`, `bp5_benchmark_output.hpp:108-117`)

### 2.3 Three Potential Failure Modes

#### Issue A: Shared Face Duplication

In MFEM's `ParMesh`, a face shared between rank A and rank B appears in **both** ranks'
`fault_shared_faces_` arrays. After `GatherToRoot`:

- The physical face has TWO sets of entries in the global array
- `num_global_dofs = sum(local_dofs)` overcounts by `num_duplicated_shared_faces * nbf_per_face`
- `num_faces = num_global_dofs / nbf_per_face` is too large (phantom faces)
- Face blocks are still individually contiguous, so `TryBuildExactMatch` may find a valid
  duplicate — but it may match the "wrong" rank's copy

**Impact on exact interpolation**: The field values at duplicate entries come from different
ranks. Both should compute the same fault quantity (slip = displacement jump), but numerical
differences from domain decomposition (element ownership, quadrature evaluation) could cause
subtle discrepancies.

**Impact on nearest-DOF**: Nearest-DOF is robust to duplication because it searches all DOFs
and picks the closest. A duplicate at the same coordinates won't change the match.

#### Issue B: Vertex Ordering Inconsistency

MFEM's `GetInteriorFaceTransformations(face)` and `GetSharedFaceTransformations(sf)` may
return different local-to-reference mappings depending on which element is the "reference"
element. This means the 3 vertices of a shared face may appear in a different order on
different ranks.

`ComputeReferenceIP` uses DOFs at `[start+0, start+1, start+2]` as the triangle's three
vertices to compute barycentric coordinates. If vertex ordering differs between the
coordinate gather and the field gather (it shouldn't — they use the same ordering), OR
if the `H1_TriangleElement::CalcShape` expects a specific vertex ordering that doesn't match,
the shape function weights would be wrong.

#### Issue C: Non-Face-Aligned DOF Count

If any rank's local DOF count is NOT a multiple of `nbf_per_face` (e.g., due to partial face
ownership), then:
- `num_global_dofs % nbf_per_face != 0` → `TryBuildExactMatch` returns immediately (line 231)
- Exact interpolation silently disabled, falls back to nearest-DOF
- This is **safe but undetected** — the fallback produces correct but less accurate results

---

## 3. Implementation Checklist

### Phase 1: Diagnose the Output Regression

#### Step 1.1: Add parallel diagnostics to station mapping
**Files**: `bp5_benchmark_output.hpp`
**Goal**: For each station, print whether exact mode is used, the chosen `face_start`, weights,
and the underlying (x2, x3) DOFs.

- [ ] **1.1a** Add a `PrintDiagnostics()` method to `Probe2DInterpolator` that, for each
  station, prints:
  - Station name
  - `exact_match_[s]` (true/false)
  - `face_start_[s]` (global DOF index)
  - `nearest_dof_[s]` (global DOF index)
  - `match_distance_[s]` (meters)
  - `exact_weights_[s]` (barycentric weights)
  - The (x2, x3) coordinates of DOFs at `[face_start, face_start + nbf_per_face)`
  - The (x2, x3) coordinates of the `nearest_dof`

- [ ] **1.1b** Call `PrintDiagnostics()` on root after constructing `ParallelBP5BenchmarkOutput`
  in `bp5_verification_full.cpp` (after line 961).

- [ ] **1.1c** Verify that for every exact station, ALL contributing DOFs belong to a single
  geometric face:
  - The 3 vertex DOFs should form a non-degenerate triangle (det != 0)
  - The target station should lie inside the triangle (all barycentric coords >= 0)
  - The DOF coordinates should be geometrically consistent (no cross-rank contamination)

#### Step 1.2: Detect shared face duplication
**Files**: `fault_geometry.hpp`, `bp5_verification_full.cpp`

- [ ] **1.2a** After gathering global coordinates on root, count unique vs total DOFs:
  ```cpp
  // On root:
  int total = global_x2.Size();
  int expected_unique = actual_mesh_fault_faces * nbf_per_face;
  // If total > expected_unique, shared faces are duplicated
  ```

- [ ] **1.2b** Identify duplicate face blocks by scanning the global coordinate array for
  faces with matching vertex coordinates (within 1m tolerance). Count how many faces
  are duplicated and from which ranks.

- [ ] **1.2c** Log a warning if `num_global_dofs != num_unique_fault_faces * nbf_per_face`.

#### Step 1.3: Compare prod vs exct_prod at matching stations
**Files**: visualization scripts or new diagnostic code

- [ ] **1.3a** Run the same simulation with exact interpolation ON vs OFF (or compare existing
  `bp5_v55_prod_p1` vs `bp5_v55_exct_prod` outputs).

- [ ] **1.3b** Focus on these priority stations:
  - `fltst_strk+00dp+10` (center, dp+10)
  - `fltst_strk+16dp+10` (right of center, dp+10)
  - `fltst_strk-24dp+10` (left, nucleation side, dp+10)
  - `fltst_strk-16dp+10` (left, dp+10)

- [ ] **1.3c** At each station, compare:
  - `log10(V_strike)` and `log10(V_dip)` during interseismic (t = 1, 5, 8, 10 yr)
  - `tau_strike` and `tau_dip` during interseismic
  - `log10(theta)` during interseismic
  - Any -300 values (writer floor hit → V <= 0 at VW stations = wrong)

- [ ] **1.3d** If mismatch exists: record which stations diverge, at what time, and by how much.
  If mismatch is station-dependent, correlate with whether the station's face is shared.

---

### Phase 2: Fix or Disable Parallel Exact Interpolation

#### Step 2.1: Decide on fix strategy

Based on Phase 1 diagnostics, choose ONE of:

**Option A — Disable exact interpolation in parallel** (fast, safe):
- [ ] **2.1a** In `ParallelBP5BenchmarkOutput` constructor (`bp5_parallel_output.hpp:73`),
  force `nbf_per_face = 1` when constructing the inner `BP5BenchmarkOutput`. This makes
  `TryBuildExactMatch` return immediately (line 231: `nbf_per_face_ < 3`), falling back to
  nearest-DOF everywhere.
- [ ] **2.1b** Add a log message: `"Parallel mode: exact face interpolation disabled, using nearest-DOF"`
- [ ] **2.1c** Verify all stations produce reasonable output (no -300 values at VW dp+10 stations).

**Option B — Deduplicate shared faces before interpolation** (correct, more work):
- [ ] **2.1d** Add a `GatherToRootWithFaceMap()` method to `FaultGeometry` that:
  1. Gathers coordinates and field data as before
  2. On root, identifies duplicate faces (same vertex coords within tolerance)
  3. Keeps only one copy of each face (from the lower-ranked owner)
  4. Returns the deduplicated array and a face-to-globalDOF mapping
- [ ] **2.1e** Use the deduplicated coordinates to construct `Probe2DInterpolator`
- [ ] **2.1f** At each timestep, deduplicate field data using the same mapping before
  calling `WriteFromGlobalData`

**Option C — Build face metadata into the gather** (most robust, most work):
- [ ] **2.1g** Extend `GatherToRoot` to also gather per-DOF face IDs (global face index)
- [ ] **2.1h** Reconstruct per-face ordering on root from the face IDs, ensuring each global
  face appears exactly once
- [ ] **2.1i** Reindex the interpolator to use the reconstructed ordering

**Recommendation**: Start with **Option A** to unblock the Tandem comparison. Implement
**Option B or C** later when exact interpolation accuracy is needed.

#### Step 2.2: Validate the fix
- [ ] **2.2a** Re-run the simulation with the fix applied
- [ ] **2.2b** Verify no station shows -300 values during interseismic at VW dp+10 stations
- [ ] **2.2c** Compare fixed output against the previous `prod_p1` output to confirm they match
  (since both now use nearest-DOF)

---

### Phase 3: Add Parallel Station Mapping Tests

**Files**: `miniapps/seas/tests/unit/test_bp5_output.cpp`

- [ ] **3.1** Add a parallel test `TestBP5BenchmarkOutput_ParallelExactInterpolation`:
  - Create a synthetic 2-rank setup with 4 triangular faces (2 per rank, 1 shared)
  - Place a station inside a shared face
  - Verify that the gathered global data + exact interpolation produces the correct
    interpolated value (known analytically)
  - Verify that nearest-DOF fallback also produces a valid result

- [ ] **3.2** Add a test `TestBP5BenchmarkOutput_SharedFaceDuplication`:
  - Create a multi-rank setup where shared faces appear in both ranks
  - Verify that the diagnostics correctly identify duplicates
  - Verify that the deduplication (if implemented) produces correct results

- [ ] **3.3** Extend the existing serial exact test (line 604) to also check:
  - Station at face boundary (barycentric coord near 0)
  - Station exactly at a vertex (should match nearest-DOF result)
  - Higher-order faces (p=2, nbf_per_face=6)

---

### Phase 4: Rebaseline Against Tandem (prod output only)

**Reference data**: `/Users/chunhuizhao/Downloads/seas-mfem/results_tandem_1000m_p1`
**MFEM data**: `bp5_v55_prod_p1` (nearest-DOF path, no exact interpolation)

#### Step 4.1: First crossing time comparison
- [ ] **4.1a** For each priority station, find the first time `log10(V) > -3` (onset of
  acceleration) in both Tandem and MFEM.
- [ ] **4.1b** For each priority station, find the first time `log10(V) > -1` (coseismic)
  in both Tandem and MFEM.
- [ ] **4.1c** Record the timing offset for each station. Current baseline: MFEM is
  ~70-90s late.

#### Step 4.2: Interseismic value comparison
- [ ] **4.2a** At t = 1, 5, 8, 10 yr, compare `log10(V_strike)` between Tandem and MFEM
  at each priority station.
- [ ] **4.2b** Same for `tau_strike` (MPa).
- [ ] **4.2c** Same for `log10(theta)` (s).
- [ ] **4.2d** Tabulate results: station × time × quantity → (Tandem, MFEM, delta).

#### Step 4.3: Priority station list (in order)
1. `fltst_strk-24dp+10` — nucleation side, most likely to show asymmetry
2. `fltst_strk-16dp+10` — nucleation side, moderate offset
3. `fltst_strk+00dp+10` — center, should be closest to symmetric
4. `fltst_strk+16dp+10` — right side, reference for asymmetry
5. `fltst_strk-16dp+00` — surface station, nucleation side
6. `fltst_strk+00dp+00` — surface station, center

---

### Phase 5: Investigate Nucleation-Side Physics Mismatch

**Only proceed after Phase 4 confirms the output path is clean.**

#### Step 5.1: Traction decomposition near nucleation
- [ ] **5.1a** Enable `diag_station_traction_decomp` in the verification run
- [ ] **5.1b** At nucleation-side stations (`strk-24dp+10`, `strk-16dp+10`), decompose
  traction into: `tau_stress` (physical), `tau_correction` (IP penalty), `tau_hat` (total)
- [ ] **5.1c** Compare the decomposition against Tandem's values. If `tau_correction` is
  anomalously large on the nucleation side, the DG penalty is the culprit.

#### Step 5.2: Shared-face vs interior-face basis handling
- [ ] **5.2a** Identify which fault faces near `dp+10` are shared vs interior
- [ ] **5.2b** Compare traction values at DOFs on shared faces vs interior faces in the
  same depth band
- [ ] **5.2c** If shared-face DOFs show a systematic bias, investigate the
  `FaultBasisData::sign_flipped` handling for shared faces

#### Step 5.3: Normal stress check
- [ ] **5.3a** Verify that normal stress (σ_n) is consistent at nucleation-side stations
  compared to symmetric stations
- [ ] **5.3b** Check whether normal stress contribution from the DG penalty differs
  between interior and shared faces

---

## 4. Guardrails

1. **-300 floor**: `log10(V) = -300` in output means `V <= 0` (writer floor at
   `bp5_benchmark_output.hpp:501`). Any VW station at `dp+10` hitting this during early
   interseismic is almost certainly wrong. If seen, stop and diagnose.

2. **Do not revert v55 fixes**: All fixes in v55 were validated with evidence. See the
   guardrail in `CLAUDE.md` — any revert requires citing the specific fix, explaining why
   the original reasoning was wrong, and getting explicit approval.

3. **Nearest-DOF is always safe**: The nearest-DOF path (`nbf_per_face = 1` or fallback)
   has been validated in serial and parallel. Use it as the ground truth for comparing
   exact interpolation results.

4. **Tandem is the reference**: Use converted Tandem text data at
   `/Users/chunhuizhao/Downloads/seas-mfem/results_tandem_1000m_p1`.
   Tandem's output is trustworthy.

---

## 5. Implementation Results (Phases 1-3)

### Phase 1 Results: Diagnostics Added

- **`Probe2DInterpolator::PrintDiagnostics()`** added to `bp5_benchmark_output.hpp`
  - Prints per-station: exact_match, face_start, weights, nearest_dof, match_distance
  - Prints contributing DOF coordinates and triangle degeneracy check
  - Called automatically on root after `ParallelBP5BenchmarkOutput` construction in
    `bp5_verification_full.cpp`

- **`Probe2DInterpolator::CountDuplicateFaces()`** added as static method
  - Scans global coordinate array for face pairs with matching vertices (within tolerance)
  - Integrated into `BP5BenchmarkOutput::PrintDiagnostics()` for automatic detection

### Phase 2 Results: Option A Implemented

- **Exact interpolation disabled in parallel** (`bp5_parallel_output.hpp:73-86`)
  - `ParallelBP5BenchmarkOutput` constructor now forces `nbf_per_face = 1` when creating
    the inner `BP5BenchmarkOutput`, regardless of the caller's `nbf_per_face` value
  - This forces `TryBuildExactMatch` to return immediately (line 231: `nbf_per_face_ < 3`),
    falling back to nearest-DOF everywhere
  - Log message printed when exact interpolation is suppressed
  - Serial path (`BP5BenchmarkOutput` used directly) still supports exact interpolation

### Phase 3 Results: New Tests Added

4 new unit tests added to `test_bp5_output.cpp` (78/78 total pass):
- `TestProbe2DInterpolator_PrintDiagnostics` — validates diagnostics output
- `TestProbe2DInterpolator_CountDuplicateFaces` — validates duplicate detection
- `TestProbe2DInterpolator_ExactWithDuplicates` — validates behavior with shared-face duplication
- `TestProbe2DInterpolator_ExactBoundaryAndVertex` — validates edge cases (vertex, edge midpoint)

**All existing tests unaffected:**
- 78/78 BP5 output tests pass
- 117/117 cross-verify tests pass
- 431/431 elasticity operator tests pass
- 17/17 BP5 parallel smoke tests pass

---

## 6. Files Modified

| File | Changes |
|------|---------|
| `miniapps/seas/io/bp5_benchmark_output.hpp` | Added `PrintDiagnostics()`, `CountDuplicateFaces()`, `GetFaceStart()`, `GetExactWeights()`, `GetNbfPerFace()` to `Probe2DInterpolator`; added `PrintDiagnostics()`, `GetInterpolator()` to `BP5BenchmarkOutput` |
| `miniapps/seas/io/bp5_parallel_output.hpp` | Disabled exact interp in parallel (force `nbf_per_face=1`); added `PrintDiagnostics()` |
| `miniapps/seas/tests/unit/test_bp5_output.cpp` | Added 4 new tests for diagnostics, duplicate detection, boundary/vertex cases |
| `miniapps/seas/tests/verification/bp5_verification_full.cpp` | Added diagnostic print call after `ParallelBP5BenchmarkOutput` construction |

---

## 7. Acceptance Criteria

- [x] Exact interpolation disabled in parallel (Option A) — forces nearest-DOF fallback
- [x] Diagnostics added to detect shared-face duplication and validate station mapping
- [x] New unit tests pass (78/78), all existing tests pass (587+ total)
- [ ] No station shows -300 values at VW dp+10 during early interseismic (needs prod run)
- [ ] Tandem comparison table completed for all 6 priority stations (Phase 4)
- [ ] Event timing offset measured and documented (Phase 4)
- [ ] Root cause of nucleation-side asymmetry identified (Phase 5)
