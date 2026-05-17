# Phase 6 — Implementation Deviations and Deferred Subscopes

Plan reference: `PLAN_onfaultstress.md` Phase 6 (lines 1543–1977),
including Phase 6.A (lines 1704–1789).

This document records the parts of Phase 6 implemented in this
cycle vs. deferred to a follow-up cycle, and the rationale.

## Implemented in this cycle

### A. Phase 6 §1 — `StressField3D` (NEW)
- **Files created:** `io/stress_field_3d.hpp`, `io/stress_field_3d.cpp`
- Six-component symmetric Cauchy stress reader; wraps six
  `DataField3D` instances in schema-v1 canonical order
  (sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz).
- `Evaluate(x, y, z)` returns a 3x3 `mfem::DenseMatrix` in Pa.
- `BBox()` returns the (asserted-consistent) shared grid bbox.
- `ContainsBBox(...)` thin gate.
- `SetInterpMode(mode)` propagates to all six readers.
- Static `ComponentName(idx)` lookup.
- Ctor asserts that all six components share the same grid axes
  / bbox (defense-in-depth for hand-edited sidecars; the standard
  Phase 5 writer produces a single `/grid/*` shared by all six
  components).

### B. Phase 6 §2 — `StressFieldCoefficient` (NEW, header-only)
- **File created:** `io/stress_field_coefficient.hpp`
- `mfem::VectorCoefficient` of dimension 6.
- `Eval(v, T, ip)` returns the six components in **schema-v1
  canonical order** (xx, yy, zz, xy, yz, xz) — NOT Voigt — with
  uniform `scale` / `offset` applied.

### D. Phase 6 §3 — `FieldProjector::ProjectStress` (additive)
- **Files modified:** `io/field_coefficient.hpp`,
  `io/field_coefficient.cpp` (additive only — existing
  `Project`, `ProjectSerial`, `ProjectVelocity` API unchanged).
- New `StressFields` POD aggregating six `ParGridFunction` shared
  pointers plus per-field (min, max) observed values.
- `static StressFields ProjectStress(sidecar_path, target_fes,
  interp)` constructs a single `StressField3D` and calls
  `Project` six times. Each call increments `call_count_` (so
  `ProjectStress` advances it by six per invocation —
  documented in the header).
- Sign convention: pure pass-through (R-501/R-502). The
  source-site flip lives in Phase 3's `bulk_stress_tensor_field`;
  no sign manipulation here.

### Tests added
- **File created:** `tests/unit/test_stress_field_3d.cpp` —
  10 tests covering construct + Evaluate, Field accessor / names,
  BBox / ContainsBBox, SetInterpMode propagation, pass-through
  (no sign flip) on negative on-disk values, missing-component
  abort, mismatched-grid abort, DenseMatrix shape and symmetry,
  OOB Evaluate abort, out-of-range Field idx abort.
  **Final run: 46 / 46 passed.**
- **Makefile updates:** `STRESS_FIELD_3D_{HEADERS,SRC,OBJ}`,
  `TEST_STRESS_FIELD_3D_{SRC,OBJ}`, build target
  `seas_test_stress_field_3d`, object rules for
  `STRESS_FIELD_3D_OBJ` and `TEST_STRESS_FIELD_3D_OBJ`,
  test target `test-stress-field-3d`, and integration into
  `test-data-projection`.
- Updated `FIELD_COEFFICIENT_HEADERS` to include the two new
  headers so the existing `field_coefficient.o` rebuilds when
  they change.

### Build verification
Building from `/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas`
via the in-tree `Makefile` requires a built MFEM library
(`libmfem.a` in `MFEM_BUILD_DIR`).  At the time of this cycle the
SAFS repo's own MFEM tree had not been built (`libmfem.a`
absent under `seas-mfem-safs/`).  As a workaround, the new files
were verified to compile and link against a sibling MFEM build
(`/Users/chunhuizhao/projects/seas-mfem/libmfem.a`) via direct
`mpicxx` invocation.  All three TUs (`stress_field_3d.cpp`,
`field_coefficient.cpp` [post-modification], `test_stress_field_3d.cpp`)
compile clean; the linked test binary runs and all 46 tests pass.
Once the SAFS-tree MFEM is built, the `make seas_test_stress_field_3d`
target works from `miniapps/seas/` without any additional changes.

## Deferred to a follow-up cycle (Tranches 2 and 3)

### C. Phase 6.A — `FaultGeometry` per-DOF accessors (DEFERRED)
- **Files to modify:** `fault/fault_geometry.hpp`,
  `domain/domain_operator.hpp` (likely needs new virtual to
  expose per-DOF reference coords and per-face -> per-DOF basis
  evaluator).
- **Reason for deferral:** the templated `FaultGeometry<MeshType>`
  is a "Files Requiring Extreme Care" header (per
  `miniapps/seas/CLAUDE.md`).  Adding per-DOF coords + basis
  requires either:
    1. extending `DomainOperator` with new virtual methods
       (`GetFaultDOFCoords3D` and `GetFaultDOFBasis`), every
       concrete implementation needing the override, OR
    2. driving the per-DOF iteration inside `FaultGeometry`
       using the existing fault DG space (`fault_basis_q_` is
       per-quadrature-point, not per-DOF, so this is a new
       traversal).
  Either approach merits a focused implementation cycle with a
  dedicated `test_fault_dof_basis` regression-gating the BP5
  vector path (and ideally a parity check against the existing
  per-quadrature-point basis).  Pushing it into the same cycle
  as the additive io/ deliverables would risk a silent breakage
  of the BP5 hot path.

### E. Phase 6 §4 — `FieldProjector::ProjectFaultPreStress` (DEFERRED)
- Depends on Phase 6.A.  Pure rotation of `σ_seas(x_i)` onto the
  per-DOF Tandem fault basis.  Cannot be implemented without
  the per-DOF coords / basis from §C.

### F. Phase 6 §5 — `FaultGeometry::ComputeSAFSParams` (DEFERRED)
- Adds a new third construction path to `FaultGeometry`
  (currently has BP2 / BP5 paths via separate ctors).  Depends
  on Phase 6.A for the per-DOF coords; depends on Phase 6 §4 for
  the sidecar-driven `tau_pre_` / `sigma_n_per_dof_`.

### G. Phase 6 §6 — `RateStateFaultOperator` SAFS-mode wiring (DEFERRED)
- **The highest-risk subscope.**  Modifies
  `fault/rate_state_fault.hpp` — on the "Files Requiring Extreme
  Care" list — at every site where the BP5 vector path reads
  `params_.sigma_n` or `sigma_n_bp5_` (Init line 303-304,
  ComputeRHS line 825-826, SolveSlipRateVectorPsi at
  rate_state_fault.hpp:831 / 838, equilibrium verification at
  line 856, plus any other site identified by the audit per
  the user-feedback memory `feedback_complete_sign_sites.md`).
- The BP5 bit-exact requirement (Phase 6 acceptance line
  1964-1968) demands that every modified site preserve the
  existing branch byte-for-byte under `safs_mode_ == false`.
- **Verification:** `make test-bp5-integration` must pass
  bit-exact relative to the pre-Phase-6 baseline.
- This subscope alone deserves its own implement -> review ->
  fix cycle, with the bit-exact baseline captured before any
  change and compared after.

### H. Phase 6 §7 — TOML schema extension (DEFERRED)
- Adds `StressConfig` to `SeasConfig`, parser/bridge updates,
  driver integration.  Depends on §G being in place.

## Why the partial implementation is still useful

The implemented tranche (StressField3D + StressFieldCoefficient
+ ProjectStress) is the **read-only, mesh-level** consumer of the
Phase 5 sidecar.  A downstream caller (e.g. the Phase 7 driver
`seas_project_stress_to_mesh`) can already:

1. Load `stress_safs.h5` via `StressField3D` (asserts schema-v1
   invariants on all six components).
2. Project the six components onto a `ParFiniteElementSpace` via
   `FieldProjector::ProjectStress` (returns six bulk
   `ParGridFunction`s in compression-positive Pa for VTK
   visualisation or downstream use).
3. Evaluate the symmetric Cauchy tensor at any mesh point via
   `StressField3D::Evaluate(x, y, z)`.

The deferred subscopes (Phase 6.A + §4-§7) are the **fault-side
plumbing** that wires the sidecar into the rate-state friction
solver's per-DOF pre-stress.  None of those are required for the
Phase 7 verification driver (which uses bulk projection, not
fault-mode), so Phase 7 can land before the deferred subscopes.

## Next steps (recommended order)

1. **Code-review** the implemented tranche (this cycle's work).
2. **Code-fix** any findings.
3. Land Phase 6.A in a focused cycle (Tranche 2 of Phase 6).
4. Land §4 + §5 (depends on 6.A).
5. Land §6 (BP5-bit-exact wiring) in its own cycle with the
   baseline-comparison protocol.
6. Land §7 (TOML schema) after §6.

Phases 7 and 8 of `PLAN_onfaultstress.md` can be implemented
independently of the deferred §4-§7, since they consume the
already-implemented tranche.
