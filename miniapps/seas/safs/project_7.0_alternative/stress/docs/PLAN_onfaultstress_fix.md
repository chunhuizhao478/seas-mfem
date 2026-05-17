# Fix Report — PLAN_onfaultstress.md

## Round 2 (R-101 through R-107)

### Summary
- Findings addressed: **7 of 7** (R-101 through R-107)
- Files modified: 1 (`PLAN_onfaultstress.md`; line count 2691 → 2761)
- Tests added: **0** (per input contract — the plan describes future code)
- Verification grep results (all passing):
  - "Compression negative" — already audited round 1, no regressions.
  - "sigma_n = np.einsum" — **0 occurrences** (R-103 fix applied; now `sigma_n_total = -np.einsum`).
  - "(recomputed at the node)" — **0 occurrences** (R-104 fix applied).
  - "via the H1(p) shape functions" — **0 occurrences** (R-107 fix applied).
  - "σ_seas(p_v)" (verifier) — **0 occurrences** (R-101 fix applied; verifier now uses σ_HZ for shear and -σ_HZ for σ_n).
  - "d_hat = np.cross(n_hat, s_hat)" — 1 occurrence at line 534, which is a *citation* of H&Z's source line numbers inside the math justification (not a recipe for our Python). Legitimate.

## Changes Made (per finding)

### R-101 [CRITICAL] — Phase 8 verifier prediction must match Phase 3's split emission
**Where:** PLAN_onfaultstress.md:2086-2108 (Phase 8 §5, fault-point and fault-cell predictor blocks).

**Change:** Replaced the symmetric-flip prediction
`σ_n_v = n_v · σ_seas · n_v, τ_s_v = s_v · σ · n_v, τ_d_v = d_v · σ · n_v`
with the asymmetric split-flip predictor:

```
σ_n_v        = -n_v · σ_HZ(p_v) · n_v    (== n_v · σ_seas · n_v)
τ_strike_v   =  s_v · σ_HZ(p_v) · n_v    (unflipped; H&Z sign)
τ_dip_v      =  d_v · σ_HZ(p_v) · n_v    (unflipped; H&Z sign)
```

Also updated the fault-cell block to explicitly say "same split-flip
predictor (σ_n flipped, shear unflipped)". This makes the verifier's
prediction algebra identical to the Phase 3 §1 emission contract, so
a correctly-written fault VTU now passes the verifier.

### R-102 [CRITICAL] — Python dip definition aligned with Tandem/SEAS (down-dip)
**Where:** PLAN_onfaultstress.md:500-540 (Phase 2 §1 `per_triangle_basis_raw`)
and PLAN_onfaultstress.md:803-816 (Phase 3 §1 cross-check stanza).

**Change in Phase 2 §1:**
- Recipe: `dips[i] = n_i × strikes[i]` → `dips[i] = strikes[i] × n_i`
- Math justification: `d = n × s` (up-dip, H&Z) → `d = s × n` (down-dip, Tandem/SEAS)
- Added a long inline comment explaining the deliberate divergence from
  H&Z and citing CLAUDE.md's `can_t1 = (0, 0, -1)` rule.

**Change in Phase 3 §1 cross-check stanza:** updated to acknowledge that
`τ_dip_seas = -τ_dip_HZ` (sign-flipped), while `τ_strike` matches
H&Z exactly. The test `test_resolve_traction_against_HZ_safod` now
asserts:
- `σ_n_seas       = -σ_n_HZ`         to 1e-12
- `τ_strike_seas  =  τ_strike_HZ`    to 1e-12
- `τ_dip_seas     = -τ_dip_HZ`       to 1e-12

This makes the Python QA τ_dip plot match the C++ production
friction solver's τ_dip sign, eliminating the parity-test failure.

### R-103 [MODERATE] — Vectorised einsum has the σ_n negation
**Where:** PLAN_onfaultstress.md:813-822 (Phase 3 §1 vectorised math).

**Change:** Replaced

```python
sigma_n = np.einsum('ki,ki->k', normals, t)
tau_s   = np.einsum('ki,ki->k', strikes, t)
tau_d   = np.einsum('ki,ki->k', dips,    t)
```

with

```python
# sigma_global is σ_HZ (compression negative). The split-flip
# contract emits sigma_n_total in compression-POSITIVE
# convention; the negation is applied to the rotated scalar.
sigma_n_total = -np.einsum('ki,ki->k', normals, t)
# tau_strike and tau_dip retain the H&Z sign (no flip).
tau_strike    =  np.einsum('ki,ki->k', strikes, t)
tau_dip       =  np.einsum('ki,ki->k', dips,    t)
```

The variable names now match the dict-key names, so an implementer
who lifts this block produces the contract-correct output.

### R-104 [MODERATE] — `traction_vec_MPa` emitted in SEAS convention
**Where:** PLAN_onfaultstress.md:757-768 (Phase 3 §1 dict entry),
PLAN_onfaultstress.md:773-789 (Phase 3 §1 math block), and
PLAN_onfaultstress.md:1060 (Phase 4 §1 fault-VTU table row).

**Change in Phase 3 §1:**
- Dict entry: rewrote `"traction_vec"` description to say
  `t = -σ_HZ · n (= σ_seas · n)` (compression-positive SEAS),
  with a usage hint noting that
  `n · traction_vec ≈ sigma_n_total` (consistent), while
  `s · traction_vec ≈ -tau_strike` (because shear keeps H&Z sign).
- Math block: added a final line `traction_vec = -t # 3-vector
  emitted to VTU; SEAS sign so n · traction_vec = σ_n_total`.

**Change in Phase 4 §1:** replaced the misleading
"global-frame traction `σ·n_node` (recomputed at the node)" with
"`σ_seas · n_cell` (area-weighted from the cell-data via
`cell_to_node_average`; **NOT** recomputed from `σ·n_node` at the
vertex). Convention: compression POSITIVE (consistent with
`sigma_n_total_MPa`)."

### R-105 [MODERATE] — `write_bulk_vtu` documents its bulk-path flip
**Where:** PLAN_onfaultstress.md:1113-1160 (Phase 4 §2 bulk VTU writer).

**Change:**
- Signature comment: `sigma_field: np.ndarray, # (N_tet, 3, 3) — H&Z convention input`.
- Cell-data table: added "compression POSITIVE (SEAS)" note to all
  six `sigma_*_MPa` rows and to `sigma_tensor_MPa`.
- Added an explicit "Bulk-path sign-flip site" paragraph stating
  that this writer is one of the two bulk-path flip sites (the
  other being Phase 5 §4), receives H&Z input, and emits
  compression-positive via `sigma_seas = -sigma_field`.

### R-106 [MODERATE] — Basis orthonormality test covers all six relations
**Where:** PLAN_onfaultstress.md:2399-2412 (Phase 4 acceptance criterion).

**Change:** Expanded
`test_fault_vtu_basis_orthonormal_at_node` from three assertions
(`s·n < 1e-9`, `|s| = 1`, `|d| = 1`) to six:
- `|s·n| < 1e-9`, `|d·n| < 1e-9`, `|s·d| < 1e-9`
- `|s| = 1 ± 1e-9`, `|d| = 1 ± 1e-9`, `|n| = 1 ± 1e-9`

Added a comment explaining that partial coverage masks Gram-Schmidt
regressions like the R-102 dip-sign bug.

### R-107 [LOW] — Phase 6.A direct evaluation (not H1 projection)
**Where:** PLAN_onfaultstress.md:1721-1729 (Phase 6.A §2 initialiser
step).

**Change:** Replaced
"Project the per-quadrature-point `fault_basis_q_` data onto the
DOF via the H1(p) shape functions of the fault DG space..." with
"Evaluate the fault-face basis at each DOF's reference coordinate
directly: invoke the `FaultBasis` evaluator at the per-DOF
reference coordinate `ip_dof` on the parent face. The basis is
per-face and discontinuous across face boundaries; this is correct
and matches the per-face storage of `tau_pre_`. **Do NOT** project
the quadrature-point values via H1 shape functions — the fault
basis is not in H1."

## Verification
- [x] R-101 (CRITICAL): Phase 8 verifier predicts σ_n via -σ_HZ flip and shear via σ_HZ (unflipped), matching Phase 3 §1 emission.
- [x] R-102 (CRITICAL): Python dip recipe and math justification both use `strikes × n` (down-dip, Tandem/SEAS); Phase 3 cross-check stanza acknowledges `τ_dip_seas = -τ_dip_HZ`.
- [x] R-103 (MODERATE): Vectorised einsum has `sigma_n_total = -np.einsum(...)`; renamed variables to match the dict keys.
- [x] R-104 (MODERATE): `traction_vec` emitted in SEAS convention (`-σ_HZ · n`); Phase 4 table row updated.
- [x] R-105 (MODERATE): `write_bulk_vtu` cell-data table marked compression POSITIVE; bulk-path flip explicitly documented.
- [x] R-106 (MODERATE): Basis-orthonormality test expanded to all six relations.
- [x] R-107 (LOW): Phase 6.A §2 says "Evaluate at the DOF reference coordinate directly" (no H1 projection).

## Unresolved Findings
- None. All seven round-2 findings applied.

## New Tests
- None added. The plan describes future code; the test specs attached
  to each REVIEW.md finding are the implementer's contract and will
  be written when the corresponding implementation lands.

## Ready for Re-Review: **YES**
The split-site sign-flip contract is now consistent end-to-end:
- **Fault path** (Phase 3 §1, Phase 4 fault VTU, Phase 8 verifier):
  internally compute on σ_HZ; emit σ_n flipped, τ_strike unflipped,
  τ_dip unflipped (Python uses Tandem-SEAS down-dip basis from
  R-102, so τ_dip carries the production sign).
- **Bulk path** (Phase 4 §2 write_bulk_vtu, Phase 5 §4
  evaluate_stress_field_on_grid, Phase 6 StressField3D consumer):
  full-tensor flip at the write boundary; downstream reads
  compression-positive Pa.
- **C++ rotation** (Phase 6 ProjectFaultPreStress, Phase 6.A
  FaultGeometry::fault_dof_basis): the C++ basis is Tandem's
  facetBasis (down-dip), matching the Python after R-102.

Recommended next pass: spot-check the Phase 0 → Phase 3 dump-file
loading path (does the round-1 R-004 `--hz-dump-file` flag still
work after R-103's variable-name changes?), and re-read the round-1
sign-convention sections of the Constraints block to confirm no
internal contradiction with the new R-102 down-dip statement.

---

## Round 1 (R-001 through R-010)

[Original round-1 fix report retained below; see git history for the
content prior to round 2.]

### Summary
- Findings addressed: **10 of 10** (R-001 through R-010, plus 3 downstream consistency edits triggered by R-001)
- Files modified: 1 (`PLAN_onfaultstress.md`; line count 2498 → 2691)
- Tests added: **0** (per input contract — the plan describes future code; test specs in the review are the implementer's contract, not artefacts to land now)
- All 5 round-1 verification anchor-greps passed.

### Critical fixes
- **R-001** (sign-flip inverting τ_strike): split into a per-component
  flip — σ_n flipped at the resolver, τ_strike/τ_dip preserved.
- **R-002** (wrong dispatch site): Phase 6 §6 targets only the BP5
  vector path (`rate_state_fault.hpp:303-304, 825-826, 831, 838, 856`)
  with correct `(2*i)` / `(2*i+1)` indexing.
- **R-003** (StressField3D docstring): replaced "Compression negative"
  with "compression POSITIVE (SEAS internal convention)".

### Moderate fixes
- **R-004**: Phase 0 dependency declaration corrected; `--hz-dump-file`
  flag added to Phase 4 CLI.
- **R-005**: New Phase 6.A inserted (`fault_geometry.hpp` per-DOF 3-D
  coords + Gram-Schmidt-orthonormalised basis accessors).
- **R-006**: μ_apparent NaN guard now triggers on `σ_n_eff ≤ 0`.
- **R-010**: Exact `np.reshape` / `np.repeat` / `np.add.at` broadcast
  recipe added.

### Low fixes
- **R-007**: "Voigt" → "schema-v1 canonical order".
- **R-008**: "ENE half-space" → "SW half-space".
- **R-009**: `DEFAULT_DUMP_PATH` declared.

### Cascading edits (round 1)
JSON summary stats and convention string, Phase 6 §4 sign-convention
contract block, and the obsolete "Sign / compression convention is
geomechanics" Constraints stanza.
