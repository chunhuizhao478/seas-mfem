# Phase 3 STOP Findings — 2026-04-23

Status: **HARD STOP on the flux-layer hypothesis.** Phase 3 §8.4 STOP condition
("hypothesis wrong") fired. Flux-layer refactor (Phases 1–2–2b) is
mechanically correct but is a no-op on the actual bug.

## Run configuration

- Binary: `seas_test_adjacent_triangle_fault_first_step_audit`
- Fixture: M0 (Cartesian-fault 2×2×2 tet mesh, Kuhn split)
- Mode: `SEAS_TEST_USE_PRECOMPUTED_FLUX=1`  → `wave.UsePrecomputedFaceFluxes(true)`
- Rank count: np=1 serial
- Build: `mfem-dev` conda env, clean build, post-Phase-1/2/2b fixes
- Totals: **46 passed / 43 failed out of 89**

Phase 2b MPI test (`test_precomputed_fluxes_phase2b_mpi_switch`) remains GREEN
at np=1/2/4: 8/16/32 of 8/16/32 assertions pass. Phase 1 remains GREEN
(`test_precomputed_fluxes_interior_identity` 78/78 +
`test_precomputed_fluxes_boundary_identity` 9/9). The failing gates are
specifically the orbit-drift / constant-state lift checks in the audit —
none of the precomputed-path unit-correctness gates fail.

## Failing-assertion inventory

### A — Runtime-path pre-existing pepper (flag-OFF baseline)

| Gate | Line | Metric | Drift |
|---|---|---|---|
| 4 bface SXY L/U, SXZ L/U | 1825–1828 | sorted signature orbit mismatch | 1.000 × 4 |
| 4 iface SXY L/U | 1829–1830 | sorted signature orbit mismatch | 1.479 × 2 |
| 4 iface SXZ L/U | 1831–1832 | sorted signature orbit mismatch | 0.814 × 2 |
| 4 manual nonfault-face SXY L/U | 1833–1834 | sorted signature orbit mismatch | 1.676 / 1.479 |
| 4 manual nonfault-face SXZ L/U | 1835–1836 | sorted signature orbit mismatch | 1.000 / 0.770 |
| 4 nonfault residual SXY L/U | 1837–1838 | sorted signature orbit mismatch | 1.676 / 1.479 |
| 4 nonfault residual SXZ L/U | 1839–1840 | sorted signature orbit mismatch | 1.000 / 0.770 |
| 4 step-1 SXY L/U | 1842–1843 | after full ADER step | 9.7e−4 / 5.0e−4 |
| 4 step-1 SXZ L/U | 1844–1845 | after full ADER step | 1.000 / 0.893 |
| 5a | 1881 | per-face F_h orbit drift (raw Interior) | 2.000 |
| 5c | 1885 | per-face F_h orbit drift (normal canonicalized) | 2.000 |
| 6a/6b/6c | 1922/1924/1926 | per-face F_h orbit drift (Gamma / Godunov / sign-symm) | 1.876 × 3 |
| 7a/7b/7c | 1991/1992/1993 | element-level lifted drift after branch fixes | 1.000 / 1.479 / 1.000 |
| 14b | 2454 | constant-I lifted boundary rhs | 1.000 |

### B — Precomputed-path (§8.3 P3 gates)

| Gate | Plan | Drift | Status |
|---|---|---|---|
| 7′ iface SXY L/U — **P3.3** | §8.3 | 1.479 / 1.479 | FAIL |
| 7′ iface SXZ L/U — **P3.3** | §8.3 | 0.814 / 0.814 | FAIL |
| 7′ bface SXY L/U — **P3.4 CRITICAL** | §8.3 | 1.000 / 1.000 | FAIL |
| 7′ bface SXZ L/U — **P3.4 CRITICAL** | §8.3 | 1.000 / 1.000 | FAIL |
| 14′ iface — **P3.5** | §8.3 | **1.999** | **FAIL (decisive)** |
| 14′ bface — **P3.5** | §8.3 | 1.000 | FAIL |

## Decisive signal: Gate 14a vs Gate 14′ interior

- **Gate 14a** (runtime path, constant-I input, interior lifted rhs): **0.000** — MFEM's runtime interior lift is orbit-covariant on constant input.
- **Gate 14′ iface** (precomputed path, constant-I input, interior lifted rhs): **1.999**.

Same input, same mesh, same element-level accumulation loop — differs only
in whether the 9×9 `nApNm1` / `nAmNm1` is built from the per-cell topology
frame (precomputed) or the per-QP runtime `(BuildFrame(CalcOrtho))` frame.

For constant Q the physical flux is identically zero (∇·F = 0). Any nonzero
`rhs` on a constant input is a basis / lift artefact of the numerical
assembly, not a flux defect. The **precomputed path produces more drift on
a ZERO-physical-flux input than the runtime path does on the same input**.
That cannot be a flux-layer bug.

## Diagnosis — where the bug is NOT

1. NOT in `PrecomputedFaceFluxes` 9×9 composition — Phase 1 P_MATRIX_IDENTITY
   is green (0.000 worst error) and the eigenvalue-invariant receipts pass
   (trace, projector orthogonality).
2. NOT in §6.2 dispatch — R5-001 hoist compile test passes; R008 fault-
   branch preservation passes; Phase 2a linear-wave equivalence is 1.1e-16.
3. NOT in Phase 2b MPI shared-face — 8/16/32 assertions green at np=1/2/4.

## Diagnosis — where the bug IS (candidates, ranked by evidence)

1. **Element-local nodal basis is not D4-covariant on the Kuhn split.**
   MFEM's L2_FECollection(GaussLobatto) gives DOF indices tied to per-
   element vertex ordering. Orbit-related elements have different local
   vertex orderings on Kuhn tets, so geometrically-coincident physical
   points do NOT map to the same local DOF index. SeisSol's modal ADER-
   DG basis is orbit-invariant by construction — that's why SeisSol's
   production TPV102 on the same Kuhn split is pepper-free while ours
   is not.
2. **Boundary-face `Loc1.Transform` + `Geometry::...::Orient[...]` tables**
   produce orbit-dependent QP-to-element-local maps on boundary faces.
3. **`ApplyMassInverse` compounds the basis non-covariance** — per-element
   M⁻¹ is applied with the element's local-DOF permutation, so even if
   raw rhs were orbit-covariant at the physical-component level, the
   M⁻¹ multiply breaks it.

## Why Phase 1/2/2b still pass despite the real bug

All three Phase 1/2/2b unit gates assert matrix-level or flux-value
equivalence for a FIXED normal. They compare `precomputed` vs `runtime`
using the SAME element and SAME face frame, so any basis/lift artefact
cancels. The orbit-symmetry gates (Gate 7 / 7′ / 14 / 14′) compare ACROSS
orbit-related elements — that's where MFEM's element-local basis
asymmetry surfaces. The unit-test harness cannot see it because unit
tests are within-element, not cross-orbit.

## Proposed next direction (three-arm triage)

### Arm 1 — Localization probes (1–2 days)

Four constant-state probes, no production-code changes:
- `G_CONST_VOL` — `ComputeVolumeRHS(Q_const)` → cross-orbit drift.
- `G_CONST_VOL_MINV` — same after `ApplyMassInverse`.
- `G_CONST_BFACE_LIFT_PROBE` — boundary-face-only lift on constant I,
  runtime vs precomputed, cross-orbit.
- `G_ORBIT_DOF_MAP_PROBE` — enumerate each orbit-pair's local DOFs at
  geometrically-coincident physical points; assert matching indices.

Decision: first probe to fail pinpoints the fix scope.

### Arm 2 — D4-equivariant fixture swap (1 day if DOF-ordering is the cause)

Plan §8.5 Option F2: replace `Mesh::MakeCartesian3D(..., TETRAHEDRON, ...)`
with a hand-built D4-equivariant tet decomposition. Rerun the full audit.
- Gate 14′ interior drops to 0 → fixture was the bug; close investigation.
- Gate 14′ interior still drifts → MFEM nodal basis itself is the bug.

### Arm 3 — Basis replacement or SeisSol bridge (2–4 days)

If Arm 2 rules out the fixture:
- **3a**: modal L2 basis with explicit D4-covariant ordering, OR
- **3b**: bridge to SeisSol's ADER kernel (MFEM handles mesh/IO, SeisSol
  handles numerics).

## Freeze policy

Until Arm 1 produces localization evidence:
- No flux-layer code changes.
- No Phase 4 (RK path) work.
- No Phase 5 (dedup, default-on, Frontera) work.
- Phase 3 audit stays wired for regression detection on any change.

## References

- Plan §8 (Phase 3 gates) — `tpv102_seissol_aligned_flux_plan_2026-04-23.md:1105-1143`
- Plan §8.4 STOP criteria — "hypothesis wrong" branch fired on this run
- Plan §8.5 Option F1/F2 — F1 (threshold relaxation) INSUFFICIENT per Gate 14′
  interior = 1.999 on constant input; F2 (fixture replacement) is Arm 2
- Audit source — `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
- Run log — this document
