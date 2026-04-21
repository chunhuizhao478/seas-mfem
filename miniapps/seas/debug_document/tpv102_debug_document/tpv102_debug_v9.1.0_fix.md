# TPV102 Debug v9.1.0 Fix Log

> Companion to `tpv102_debug_v9.1.0_debug_plan.md` and `REVIEW.md`
> (rev 3).  One entry per substantive change applied on the laptop side
> of the v9.1.0 plan (Pelties-9 follow-through for R-V91-A speckle and
> R-V91-B dip-transient).  No Frontera runs are recorded here; those
> wait on explicit user approval per `feedback_frontera_approval.md`.
>
> **Status:** C0 + C1 + C2 applied.  530/530 laptop regression assertions
> green (13 test binaries).  Plan §8 C3 (verification helper script),
> C5 (sbatch), and C6 (post-Frontera closure) still pending.

---

## Background — why a fix doc is needed (REVIEW R-009)

`WriteFaultSurfaceVTU` in `io/paraview_output.hpp` has bounced between
CellData and PointData output across the project's history:

| Commit / version | Output form | Why | Consequence |
|---|---|---|---|
| pre-0756b80 (pre-2026-04-07) | CellData (one value per cell) | Original ParaView integration | BP5 worked; TPV102 not yet live |
| **0756b80 (2026-04-07)** | **PointData (one value per vertex)** | Commit message: "CellData→PointData for smooth interpolation" | Introduced **R-001** (H-V91-A4): QP DOFs placed at reference-triangle corners, visible as speckle in TPV102 job 7667881 slip-rate VTU |
| **v9.1.0 C0 (this change)** | **CellData (one value per cell, averaged over nbf DOFs)** | REVIEW R-001 (2026-04-20) identifies QP/vertex-position mismatch; Pelties-9 fix is code-path-correct so the artifact is purely in the output writer | Speckle resolved; BP5 BR2 latent bug (R-005 cross-face read + OOB) resolved as free side effect |

Without this doc, the round-trip (CellData → PointData → CellData)
looks like an undocumented revert.  Future readers grep-ing for
`<PointData>` in git history would find 0756b80's justification
("smooth interpolation") and re-introduce the bug.  Per REVIEW.md §11
R-009: **this doc's existence is itself the regression gate.**

---

## C0 — source fixes (single commit)

### C0.1  `io/paraview_output.hpp::WriteFaultSurfaceVTU` — PointData → CellData rewrite

**Review IDs resolved:** R-001 (CRITICAL), R-002 (MODERATE, subsumed),
R-004 (LOW, docstring), R-005 (LOW, BP5 BR2 latent).

**What changed.**

1. Accumulation loop inside `process_face` (lines ~540-605 post-fix):
   - Replaced 12 per-vertex std::vector<double> `f_*` with per-cell
     `c_*`.  Each `c_*.push_back(...)` runs once per face rather than
     once per vertex.
   - Inner loop bound `for (int k = 0; k < 3; k++)` → `for (int k = 0;
     k < nbf; k++)` with an arithmetic-mean accumulator (`a_* / nbf`).
     Fixes R-002 (k<3 hardcode dropped QPs 3..nbf-1 for order>=2) and
     R-005 (nbf=1 BR2 loop cross-read faces `i`, `i+1`, `i+2` and
     OOB-read the last 1-2 faces per rank).

2. Per-rank VTU XML block (lines ~675-701 post-fix):
   - `<PointData>` → `<CellData>`.  DataArrays now emit 1 value per
     cell rather than 3 values per vertex.

3. PVTU index block (lines ~717-730 post-fix):
   - `<PPointData>` → `<PCellData>`.

4. MPI barrier guard (lines ~668-673 post-fix):
   - `#ifdef MFEM_USE_MPI { MPI_Barrier(mesh_.GetComm()); } #endif`
     wrapped in `if constexpr (std::is_same_v<MeshType, ParMesh>)` so
     serial `Mesh` specialisation compiles.  Latent bug exposed when
     `test_fault_surface_vtu_continuity` became the first serial
     caller of the writer.

5. Docstring (lines ~505-518):
   - Rewritten to name the CellData contract and explicitly reference
     R-001 / H-V91-A4 so a future auditor who sees the round-trip in
     git history has an in-source anchor.

### C0.2  `fault/fault_basis.hpp::ComputeOrientedFrame` — explicit dip normalize

**Review IDs resolved:** R-003 (LOW).

**What changed.** Lines ~418-431: between the `dip = strike × n_raw`
cross product and the `tangent1[d] = dv[d]` assignment, inserted:

```cpp
const real_t d_len = std::sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
MFEM_VERIFY(d_len > 1e-12, "ComputeOrientedFrame: degenerate dip ...");
const real_t d_inv = 1.0 / d_len;
dv[0] *= d_inv; dv[1] *= d_inv; dv[2] *= d_inv;
```

Restores orthonormality to ~1 ULP on the dip vector (pre-fix: ~5 ULP
drift per R-003 numerical analysis).  Plan §3.4 H-V91-B1 explicitly
notes this is a defensive latent-class removal, NOT an explanation
for R-V91-B's 0.3 MPa hypocenter dip spike.

### C0.3  `fault/fault_basis.hpp::FaultBasis` — `ComputeOrientedFrame` visibility

**Decision:** promoted `static ComputeOrientedFrame` from `private:`
to a dedicated `public:` block (line ~336).  Required by
§4.2 test which needs to call the routine on 1000 perturbed normals
without constructing a mesh per trial.  Method is static with no
class-state dependency; promotion is API-only.

**Review IDs resolved:** prerequisite for R-003 test (§4.2).

### C0.4  Docstring-only rewrite in `paraview_output.hpp`

**Review IDs resolved:** R-004 (LOW).

**What changed.** Lines 505-518 docstring: removed the stale
"smooth interpolation within each face in ParaView" claim and added
explicit reference to R-001 / H-V91-A4 speckle and R-005 BR2 BP5
side-effect resolution.

---

## C1 — `test_fault_surface_vtu_continuity` (plan §4.1)

**File:** `tests/unit/test_fault_surface_vtu_continuity.cpp` (NEW).
**Makefile:** new SRC/OBJ variables, link rule, obj-rule, phony target
`test-fault-surface-vtu-continuity`.

**What it tests.**

- **Test 1 (nbf=3, linear field).**  Builds a 2×2×1 hex-tet cartesian
  mesh, collects up to 6 interior faces, feeds the VTU writer a
  linear field `f(x2,x3) = 0.5 + 0.1*x2 + 0.2*x3` evaluated at each
  face's three degree-2 quadrature points, writes to a temp dir,
  parses the VTU XML, and asserts each cell value equals
  `f(face_centroid)` to `5e-9 * max(1, |f|)` relative across 12 output
  fields.  Tolerance is loose enough to absorb the ASCII printed-
  precision roundoff (`paraview_output.hpp:666` uses
  `setprecision(10)`) but tight enough to catch the pre-R-001 writer
  (residue `O(h * |grad f|) ~ 1e-2`, 7 OOM above the threshold).

- **Test 2 (nbf=6, R-002 regression).**  Sets `local_slip_rate[2*d]=d`
  for globally-unique per-QP IDs, asserts each face's cell value
  equals the mean of its six per-QP IDs (=`fi*6+2.5`).  A pre-fix
  writer reading only the first 3 QPs would emit `fi*6+1.0`; test
  includes a defensive `!=` assertion on that pre-fix value.

**Result:** 63/63 assertions PASS on the post-C0 writer.

---

## C2 — `test_fault_basis_dip_strike_symmetry` (plan §4.2)

**File:** `tests/unit/test_fault_basis_dip_strike_symmetry.cpp` (NEW).
**Makefile:** new SRC/OBJ variables, link rule, obj-rule, phony target
`test-fault-basis-dip-strike-symmetry`.

**What it tests.**

- **Canonical frame.**  Feeds `ref_normal=(0,-1,0)`, `up=(0,0,1)`,
  asserts the expected exact IEEE-754 canonical strike=(1,0,0) and
  dip=(0,0,-1).  9 TEST_NEAR assertions at tol `1e-15`.

- **1000 perturbed-normal trials.**  xorshift64 PRNG draws
  `delta ~ U(-0.1/sqrt(3), +0.1/sqrt(3))^3`, constructs
  `n_raw = (ref_normal + delta) * 2.7` (unnormalized; writer
  normalizes), calls `ComputeOrientedFrame`, checks
  `|normal|, |dip|, |strike| == 1` and all three pairwise dots = 0.
  Aggregates into 2 PASS/FAIL counters (6000 underlying checks).

- **Tolerance (R-006, rev 3).**  `tol_len = 2.0 * eps`, `tol_dot =
  10.0 * eps`.  Pre-C0.2 (R-003 revert) the dip-length drift is
  1-4 ULP which would fail the 2·eps length check — so the test
  functions as a regression gate for accidental revert of the
  explicit dip-normalize block.

**Result:** 17/17 top-level assertions (6000 aggregate checks) PASS
on post-C0.2 code.

---

## Rev-3 follow-ups (REVIEW R-006 / R-007 / R-008 / R-009)

Applied 2026-04-20 post-rev-3 review:

- **R-006.**  Tightened `test_fault_basis_dip_strike_symmetry.cpp`
  `tol_len` from `10 * eps` to `2 * eps`.  Reverting the R-003 normalize
  block now fails the length check instead of silently passing.  Log
  strings updated to cite the new tolerance.
- **R-007.**  Replaced the inside-lambda `MFEM_ASSERT(nbf > 0, ...)`
  (no-op in Release builds) with `MFEM_VERIFY(nbf > 0, ...)` at the
  `WriteFaultSurfaceVTU` function entry.  `nbf=0` now fails loudly in
  every build rather than propagating `1.0/0.0 = +Inf` through the
  averaging and emitting NaN CellData.
- **R-008.**  Added a trailing empty `private:` block in
  `fault_basis.hpp` after the `ComputeOrientedFrame` method.  The
  mid-class `public:` inserted by R-003 is now explicitly re-closed;
  future members appended to the class cannot silently inherit
  `public:` visibility.
- **R-009.**  This document.

---

## Pending (not in this commit)

- **C3** — `tpv102/verification/vtu_speckle_audit.py` post-processing
  helper (plan §5.1/§5.2).  Not on the critical path for R-V91-A
  closure.
- **C5** — `jobs/tpv102/tpv102_200m_p1_12.0s_400rank_v91.sbatch`.
  Gated on user approval per `feedback_frontera_approval.md`; draft
  parameters already specified in plan §6.2.
- **C6** — R-V91-A / R-V91-B closure documentation after the 12 s
  Frontera re-run lands.  Dependent on C5.
- **Separate ticket** — 2026-04-14 RK4/station-output consistency
  finding (plan §10 carry-forward).  Independent of R-V91-A/B; has
  its own proposed diff.

---

## Test summary post-C0+C1+C2 (laptop, with R-006/R-007/R-008 cleanup)

| Test binary | Assertions | Failed |
|---|---|---|
| `seas_test_wave_operator` | 17 | 0 |
| `seas_test_wave_bc` | 10 | 0 |
| `seas_test_pml` | 5 | 0 |
| `seas_test_fault_basis` | 219 | 0 |
| `seas_test_fault_face_flux` | 19 | 0 |
| `seas_test_fault_face_flux_frame` | 12 | 0 |
| `seas_test_fault_face_flux_frame_and_flux` | 14 | 0 |
| `seas_test_godunov_identity_normal_reversal` | 5 | 0 |
| `seas_test_fault_flux_interior_vs_shared_branch_equivalence` | 4 | 0 |
| `seas_test_godunov_interior_equal_sides_identity` | 135 | 0 |
| `seas_test_fault_face_flux_per_side_assembly` | 10 | 0 |
| **`seas_test_fault_surface_vtu_continuity` (NEW)** | **63** | **0** |
| **`seas_test_fault_basis_dip_strike_symmetry` (NEW)** | **17** | **0** |
| **TOTAL** | **530** | **0** |

Both production drivers (`seas_tpv102_driver`, `seas_bp5_full`) re-link
clean.  DIAG build smoke-test still clean under
`SEAS_EXTRA_CPPFLAGS=-DSEAS_DIAG_FAULT_FLUX` (inherited from v9.0.0
§14.4).

---

## References

- `tpv102_debug_v9.1.0_debug_plan.md` (plan).
- `REVIEW.md` rev 1, rev 2, rev 3 — the authoritative finding list
  this doc closes out.
- `tpv102_debug_v9.0.0_seissol_flux_comparison.md` — the upstream
  Pelties-9 per-side flux fix.  Orthogonal to v9.1.0; neither R-V91-A
  nor R-V91-B touches the solver.
- de la Puente et al. 2009 · Pelties et al. 2012 — the per-side flux
  formulation that v9.0.0 ported.
- Dumbser & Käser 2006 — the welded-face upwind flux used at non-fault
  interior faces (unchanged by v9.0.0 / v9.1.0).
