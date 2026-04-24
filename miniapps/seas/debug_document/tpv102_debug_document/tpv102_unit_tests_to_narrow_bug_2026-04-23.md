# TPV102 — Unit Tests to Narrow the Adjacent-Triangle Bug

Date: 2026-04-23 (revised — audit has isolated the failure to the second-step **non-fault** face contribution)

Companion to [tpv102_fresh_review_recent_changes_2026-04-23.md](./tpv102_fresh_review_recent_changes_2026-04-23.md).

## 1. Symptom Recap (post-audit, revision 3)

The first-step audit (`tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp:732`,
Makefile target at `Makefile:149`) and the mesh-configurable existing test
(`tests/unit/test_adjacent_triangle_fault_uniformity.cpp:147`,
`SEAS_TEST_FAULT_NX/NY/NZ`) have narrowed the bug through to the second step. Current
state:

| step-0 gate | result |
|---|---|
| fault-local solve uniform | PASS |
| first fault-face lift uniform pre- and post-`M^{-1}` | PASS |
| manual chain = production `AdvanceADER(Q=0)` | PASS (bit-exact) |

| step-1 gate on refined M_ref 4×2×4 (serial, rerun 2026-04-23) | result |
|---|---|
| `I_plus / I_minus` at fault QPs uniform across orbit | PASS |
| `Q_imp_plus / Q_imp_minus` uniform across orbit | PASS |
| `F_h` (fault) uniform across orbit | PASS |
| fault DOF state uniform across orbit | PASS |
| fault-face-only lift uniform across orbit | PASS (0.0) |
| volume-only update uniform across orbit | PASS (~1e-15) |
| boundary-face-only update uniform across orbit | **FAIL** (SXY 1.00 / SXZ 1.00 both sides) |
| interior-non-fault-only update uniform across orbit | **FAIL** (SXY 1.749 / SXZ 0.649 both sides) |
| manual (boundary + interior-non-fault) uniform across orbit | **FAIL** (SXY 1.749 / SXZ 1.00 both sides) |
| closure: `(manual nonfault) − (full − fault − volume)` | PASS (7.94e-13 ≪ 1e-11) |
| full step-1 bulk signature uniform across orbit | **FAIL** (SXY 1.70e-3 / 9.55e-4; SXZ 1.00 / 1.00) |

Failing full-bulk signature on M_ref 4×2×4 serial (sorted-signature orbit drift):

| channel | lower | upper |
|---|---|---|
| SXY | 1.70e-3 | 9.55e-4 |
| SXZ | 1.00 | 1.00 |

M_ref 4×2×4 `np=4`, existing 20-step uniformity test, worst step-19 spreads:

- slip_rate 9.21e-5
- **tau1_corr 1.03e+01**
- tau2_corr 3.60e-6
- sigma_n_corr 4.07e-6

**What is now established**:

- The fault path — `EvaluateADER`, canonical rotation (`can_n`, `can_t1`, `can_t2`, `T_can`, `Tinv_can`), fault-face lift + `M^{-1}` — is clean on the first two steps of a planar fault with uniform material + uniform persistent nucleation.
- The volume-only update is clean.
- **Both** non-fault face branches break at step 2 on orbit-uniform nonzero `I`:
  - interior-non-fault contribution is orbit-asymmetric (SXY 1.749 / SXZ 0.649).
  - boundary contribution is orbit-asymmetric (SXY 1.00 / SXZ 1.00).
- Decomposition arithmetic is clean: manual (boundary + interior-non-fault) agrees with `full − fault − volume` to 7.94e-13, so H_VOLUME_COUPLING is **eliminated**.
- **Channel-to-branch attribution** (from matching drift magnitudes):
  - **SXY drift is driven by interior-non-fault** (1.749 on interior ≈ 1.749 on manual).
  - **SXZ drift is driven by boundary** (boundary and manual both 1.00; interior alone is 0.649).
  - The screenshots' sigma_xz lobe/quadrant pattern therefore traces to the boundary branch; sigma_xy's weaker but real pattern traces to the interior non-fault branch.

## 2. Hypotheses (revision 4, post-run)

H_VOLUME_COUPLING is **eliminated** by the closure gate (7.94e-13). H_IFACE and H_BFACE
are **both confirmed**. The new live hypotheses refine each branch's internal mechanism.

| # | hypothesis | post-audit |
|---|---|---|
| H_IFACE | Interior non-fault branch is not orbit-uniform — CONFIRMED by Gate 3A failure (SXY 1.749, SXZ 0.649 on M_ref serial). | confirmed |
| H_BFACE | Boundary branch is not orbit-uniform — CONFIRMED by Gate 3B failure (SXY 1.00, SXZ 1.00 on M_ref serial). | confirmed |
| H_IFACE_NOR | Inside H_IFACE: the asymmetry is induced by MFEM's uncanonicalized `nor = CalcOrtho(J_F)` in the interior branch. Two orbit-equivalent faces can have `nor` in opposite directions depending on the Kuhn split's Elem1/Elem2 ordering; `flux_.Interior(nor, I_self, I_nbr, F_h)` must be symmetric under `(nor → −nor, swap I_self ↔ I_nbr)` for the lift to be orbit-invariant. | new, primary |
| H_BFACE_NOR | Inside H_BFACE: the absorbing (or natural) BC flux `AbsorbingTotal(nor, I_self, Q_bg, F_h)` inherits the outward-normal sign, but the local volume-element on the other side of the cube (y=0 vs y=L) has a mirror-image orbit expectation that only holds if `AbsorbingTotal` is sign-anti-invariant for the state-dissipation term. A residual in the Rusanov dissipation term that does not flip sign under `nor → −nor` would produce exactly the 1.00 symmetric drift observed. | new, primary |
| H_IFACE_ACCUM | Inside H_IFACE: even with consistent `nor`, per-DOF accumulation `rhs ±= w·shape(i)·F_h[c]` with `shape1(i)` indexed by element-local DOF order could mis-match across orbit because the Kuhn split gives two orbit-equivalent tets different element-local DOF orderings. This would NOT cancel at the element-sorted-signature level, so it cannot explain a nonzero sorted-signature drift. Flagged for completeness only. | retired |
| H_SFACE | Shared non-fault face branch is orbit-asymmetric at `np ≥ 2` — not yet tested serially/parallelly; plan keeps this as a residual gate. | untested (parallel only) |

Retired by previous audit:

- H_CK (CK time predictor on nonzero Q) — step-1 I_plus/I_minus at the fault were uniform, i.e. the CK predictor ran fine on the orbit-uniform step-0 Q.
- H_TRACE (face-trace reconstruction at fault) — step-1 I_plus/I_minus uniform.
- H_LIFT2 (fault lift with real imposed states) — step-1 fault-face-only lift uniform.
- H_KUHN_FIXTURE, H_DIP_NULL_MODE — refinement strengthens the failure.
- H_FAULT_BASIS — Gate 1–3 clean.

## 3. Design: Non-Fault Face Split and Boundary vs Interior Isolation

### Mesh tiers (unchanged, for continuity)

| tier | NX×NY×NZ | fault tris | ranks used |
|---|---|---|---|
| M0 | 2×2×2 | 8 | 1 |
| M_ref | 4×2×4 | 32 | 1, 4 |
| M_big | 8×4×8 | 128 | 8, 10 |

All exposed through `SEAS_TEST_FAULT_NX/NY/NZ`. MPI budget ≤ 10. No Frontera.

### Test 1 — Branch-level 4-way split  *(shipped and run)*

Already present in the shipped audit at
`test_adjacent_triangle_fault_first_step_audit.cpp:1119–1134`. Gate 3A (interior
non-fault) and Gate 3B (boundary) both fail on M_ref serial. Gate 3C (closure) passes.
No further work required on Test 1 itself; the per-branch drill-downs below consume
its attribution.

### Test 2 — `test_interior_nonfault_flux_face_orbit_identity`  *(follow-up to Gate 3A, primary)*

Extend the audit's `RunADERNonFaultFaceAudit` helper with **per-face orbit recording**
on the interior-non-fault branch only.

For every interior non-fault face `f` adjacent to the fault-adjacent element set:

1. Record the mesh face index, its two elements `(e1, e2)`, the per-QP `nor` (post-CalcOrtho, post-normalization), `w = ip.weight * nor_len`, `I_self[c]`, `I_nbr[c]`, and the per-QP `F_h[c] = flux_.Interior(nor, I_self, I_nbr)`.
2. Record the pre-`M^{-1}` contribution `shape1(i) * F_h[c]` and `shape2(i) * F_h[c]` summed per component per DOF index `i` (an `ndof × NUM_STATE` array per face).
3. Classify each face into orbit under the y=L/2 reflection (key: `(sign(cy − L/2), round(|x − L/2|·1e6), round(|z − L/2|·1e6))`). Within each orbit, compare the **sorted** per-component `F_h` tuples across faces.

**Gates 2a–2d** added at the audit level:

- **2a**: for every orbit of interior non-fault faces, the sorted tuple of `F_h[c]` across orbit members agrees to 1e-12 componentwise.
- **2b**: with `nor` forcibly replaced by a canonicalized, always-outward-positive per-face normal (via the fault basis's `can_n` pattern extended to non-fault faces), repeat 2a. If 2a fails and 2b passes, **H_IFACE_NOR confirmed**.
- **2c**: with `nor` left as-is but `flux_.Interior(nor, I_self, I_nbr)` replaced by a symmetrized variant `0.5 * (flux_.Interior(nor, I_self, I_nbr) + flux_.Interior(−nor, I_nbr, I_self))`, repeat 2a. A PASS isolates the bug to the non-symmetrization of `flux_.Interior` under `(nor, L↔R)` flip.
- **2d**: report the per-orbit max pairwise drift in the Rusanov dissipation term alone (`α/2 · (I_nbr − I_self)` component) vs the central-flux term (`(A(nor)·I_self + A(nor)·I_nbr)/2`). Identifies which of the two Godunov-flux pieces is the one that breaks.

**Fixture**: M_ref serial. M0 is too coarse for the orbit attack.

### Test 3 — `test_boundary_face_orbit_identity`  *(follow-up to Gate 3B, primary)*

Analogous extension for the boundary branch.

For every boundary face `f` adjacent to the fault-adjacent element set:

1. Record `(f, e1, bdr_attr, FaceBC classification, nor, w, I_self[c], F_h[c])`.
2. Orbit-classify faces by outward-normal direction and mirror position across y=L/2.
3. Gates:
   - **3a**: sorted per-orbit `F_h` agrees to 1e-12.
   - **3b**: replace `AbsorbingTotal(nor, I_self, bulk_bg_scaled, F_h)` with its sign-symmetrized counterpart `0.5 * (AbsorbingTotal(nor, I_self, bg) + S · AbsorbingTotal(−nor, I_self, bg))` where `S` is the stress-component sign-flip that a reflection imposes (the fault path's `Tinv_can` then `T_can` round-trip is the correct analog). A PASS isolates H_BFACE_NOR to the non-anti-symmetry of the absorbing dissipation.
   - **3c**: toggle `SetAbsorbingBackground(Q_bg = 0)` vs `Q_bg = a nonzero orbit-symmetric background`. If 3a fails at `Q_bg = 0` but passes at a nonzero orbit-symmetric `Q_bg`, the dead-code branch under zero background is the culprit (a long-standing risk per CLAUDE.md).
   - **3d**: split natural vs free-surface vs dead absorbing by constructing three mini-fixtures that set `bc.natural_attrs`, `bc.free_surface_attrs`, and default (absorbing) respectively on the same exterior attribute. The first fixture in which 3a passes identifies the only-correct BC mode and names the wrong ones.

**Fixture**: M_ref serial.

### Test 4 — `test_shared_nonfault_flux_orbit_identity`  *(Gate 3D follow-up, parallel)*

Same structure as Test 2 but targets `ComputeADERSharedFaceFluxRHS`'s non-fault branch.
Required to confirm the step-19 parallel drift (`tau1_corr = 1.03e+01` at np=4) is
consistent with the serial interior-non-fault + boundary attribution plus the shared
seam, not a new mechanism.

**Fixture**: M_ref `np=4`, M_big `np=8`, M_big `np=10`.

### Test 5 — `test_multistep_nonfault_feedback_scaling`  *(diagnostic, unchanged)*

Non-gating. Runs steps 0–6 and reports per-step orbit spread of SXY and SXZ on
fault-adjacent tets to scope the regression-gate step count.

**Fixture**: M_ref serial + `np=4`.

## 4. Run Order and Rank Budget

```
shipped: first-step audit (step-0, step-1 fault, step-2 3-way split) → non-fault residual fails
                        │
                        ▼
Test 1 — extend to 4-way split (Gates 3A / 3B / 3C / 3D)
                        │
       ┌────────────────┼────────────────┐
       │                │                │
   3A fails          3B fails        3D fails (np≥2 only)
       │                │                │
       ▼                ▼                ▼
   Test 2           Test 3           Test 4
  (iface          (bface         (shared iface
   nonfault)       branch)        nonfault)

                          (orthogonal, diagnostic)
                          Test 5 — multistep feedback scaling
```

Per-test MPI rank plan (all ≤ 10):

| Test | Serial | Parallel |
|---|---|---|
| 1 | M0, M_ref | M_ref `np=4`, M_big `np=8`, M_big `np=10` |
| 2 | M_ref | M_ref `np=4` |
| 3 | M_ref | M_ref `np=4` |
| 4 | — | M_ref `np=4`, M_big `np=8`, M_big `np=10` |
| 5 | M_ref | M_ref `np=4` |

Commands (matching the form the user already ran):

```
make -C miniapps/seas \
  seas_test_adjacent_triangle_fault_first_step_audit

./miniapps/seas/seas_test_adjacent_triangle_fault_first_step_audit

SEAS_TEST_FAULT_NX=4 SEAS_TEST_FAULT_NY=2 SEAS_TEST_FAULT_NZ=4 \
  ./miniapps/seas/seas_test_adjacent_triangle_fault_first_step_audit

SEAS_TEST_FAULT_NX=4 SEAS_TEST_FAULT_NY=2 SEAS_TEST_FAULT_NZ=4 \
  mpiexec -n 4 \
  ./miniapps/seas/seas_test_adjacent_triangle_fault_first_step_audit

SEAS_TEST_FAULT_NX=8 SEAS_TEST_FAULT_NY=4 SEAS_TEST_FAULT_NZ=8 \
  mpiexec -n 8 \
  ./miniapps/seas/seas_test_adjacent_triangle_fault_first_step_audit

SEAS_TEST_FAULT_NX=8 SEAS_TEST_FAULT_NY=4 SEAS_TEST_FAULT_NZ=8 \
  mpiexec -n 10 \
  ./miniapps/seas/seas_test_adjacent_triangle_fault_first_step_audit
```

## 5. What Each Outcome Tells Us

| Outcome | Conclusion |
|---|---|
| 3A fails, 3B/3C/3D pass | **Primary-suspect confirmed**: interior non-fault DG face flux is not orbit-uniform on the step-2 Q. Focus on `nor = CalcOrtho(J_F)` routing vs canonical-normal handling in the non-fault interior branch of `ComputeADERFaceFluxRHS`. Run Test 2. |
| 3B fails, 3A pass | Boundary branch is the culprit. Run Test 3 to split natural / free-surface / absorbing. |
| 3C fails | Decomposition arithmetic mismatch (H_VOLUME_COUPLING). Rewrite the audit's step-2 split to re-accumulate contributions non-destructively; the underlying code may actually be correct. |
| 3A pass, 3B pass, 3C pass, 3D fails | MPI shared non-fault branch. Run Test 4. Parallel-only code path. |
| 3A pass, 3B pass, 3C pass, 3D not exercised, but `np ≥ 2` still fails | The parallel extra drift comes from a combination of shared-non-fault and rank-local reduction ordering; cover with Test 4 and then a small reduction-order consistency probe. |
| All pass | Re-examine the subtraction identity in the audit; the nonfault-residual drift may be recomputed from instrumentation that does not match the production accumulation order. |

## 6. What I Would Not Do

- Reopen any fault-path tests — the audit has cleanly ruled out `EvaluateADER`, canonical fault rotation, and the fault lift + `M^{-1}`.
- Revisit the split-prestress contract — both `test_persistent_nuc_prestress_channel` and Gate 1 stay clean.
- Exceed 10 ranks locally.
- Pivot to a Frontera run before Test 1 has produced a concrete failing gate. Per the project's Frontera-run policy, Frontera is reserved for scale-only reproducers we cannot show locally; every hypothesis above is reproducible within M_big `np=10`.

## 7. Final debug run (2026-04-23, Gates 5–14 all shipped — errors solid)

Added gates:

- **Gate 5** — per-face F_h orbit drift in interior non-fault branch (raw / symmetrized / canonical-normal).
- **Gate 6** — per-face F_h orbit drift in boundary branch (baseline / Godunov / sign-symmetrized; kept as diagnostic only — per-face drift is not a clean test for boundaries because `I_self` depends on per-face shape ordering).
- **Gate 7** — element-level lifted orbit drift after applying each symmetrization in isolation and together. **Definitive** diagnostic for branch-level fixes.
- **Gate 8** — within-outer-side boundary F_h orbit breakdown (by the six cube sides).
- **Gate 9** — Gamma vs Godunov lifted bface drift comparison.
- **Gate 10** — zero-input sanity (with I=0, interior and boundary rhs should be zero bit-exactly). PASSES.
- **Gate 11** — constant-input boundary D4 covariance probe (constant I so `I_self = val` at every face QP regardless of shape ordering). Measures relative F_h drift.
- **Gate 12** — boundary-face `nor` orientation audit (does `CalcOrtho(J_F)` return outward?). PASSES: all 128 boundary faces have outward `nor`.
- **Gate 13** — raw per-face dump on x=0 side (six sample faces; confirms F_h bit-identical across them).
- **Gate 14** — constant-I lifted element rhs orbit drift. **Key isolator**: removes all `I_self`-through-shape variation; any nonzero bface drift here is purely a mesh-topology / DOF-attachment issue.

### Final numbers (M_ref 4×2×4 serial, confirmed same on 8×4×8 serial and 4×2×4 np=4)

| gate | metric | value | verdict |
|---|---|---|---|
| 5a | per-face F_h drift, interior raw `Interior(n,L,R)` | 2.000 | FAIL |
| 5b | per-face F_h drift, symmetrized `0.5*(Interior(n,L,R)+Interior(-n,R,L))` | **5.68e-14** | **PASS** |
| 5c | per-face F_h drift, canonical normal only (no L↔R swap) | 2.000 | FAIL |
| 7a | element-level lifted drift, iface symmetrized only | iface **5.98e-16** / bface 1.00 | interior fixed |
| 7b | element-level lifted drift, bface symmetrized only | iface 1.749 / bface SXZ **8.30e-19**, bface SXY **1.00** | SXZ fixed; SXY residual |
| 7c | element-level lifted drift, both symmetrized | iface **0** / bface SXY **1.00**, bface SXZ **8.30e-19** | SXY residual remains |
| 10 | zero-input interior and boundary rhs max | 0 / 0 | PASS (audit is clean) |
| 11 | relative F_h drift on constant I | stress 1.4e-16–1.9e-16 (ULP); VX 1.0e-6 (non-ULP signal) | flux is ULP-D4-covariant |
| 12 | inward-pointing boundary normals | 0 of 128 | PASS |
| 13 | six x=0 faces under constant I | bit-identical F_h across all six | direct proof Gate 11 is ULP |
| 14a | constant-I lifted iface orbit drift | 0 | PASS |
| 14b | constant-I lifted bface orbit drift | **1.000** | **FAIL — fixture-level** |
| — | scale invariance | bface drift = 1.000 on 2×2×2, 4×2×4, 8×4×8 | Kuhn artefact is h-independent |

### Three solid mechanisms

1. **H_IFACE_NOR — CONFIRMED CODE BUG** (primary).
   `flux_.Interior(nor, I_self, I_nbr)` is not invariant under the simultaneous `(nor → −nor, I_self ↔ I_nbr)` swap on orbit-mirror faces. MFEM's `nor = CalcOrtho(J_F)` orientation on an interior face depends on Kuhn Elem1/Elem2 traversal; orbit-equivalent mirror faces receive opposite `nor` directions with L/R unswapped. Fix: replace the call with its `(n, L↔R)` symmetrization. Gate 5b + Gate 7a (iface residual 5.98e-16 → 0 on 8×4×8 scale too) are ironclad.

2. **H_BFACE_NOR (SXZ) — CONFIRMED CODE BUG** (secondary).
   `FreeSurfaceTotal(nor, I_self, bg, F_h)` is not anti-symmetric in its SXZ response under `nor → −nor`; the Gamma and Godunov variants both show it. Fix: sign-symmetrize the BC flux call. Gate 7b SXZ residual 8.30e-19 is ULP.

3. **H_KUHN_FIXTURE — FIXTURE ARTEFACT, NOT A CODE BUG** (reopened).
   Gate 14b proves under a literally constant input `I` the boundary contribution's per-element sorted signature has 100% relative drift. Per-face F_h is ULP-identical (Gates 11 / 13), so the only remaining degrees of freedom are the per-tet `{(w, shape1(i), dof_i)}` attachments — which are NOT D4-orbit-equivalent under `MakeCartesian3D(..., TETRAHEDRON, ...)`'s per-hex Kuhn diagonal. This defect is scale-invariant (confirmed 2×2×2, 4×2×4, 8×4×8 all give 1.000). **Refinement does not dilute it** — the earlier elimination of H_KUHN_FIXTURE was based on the wrong assumption that it would.

### Auxiliary finding (flag for later)

Gate 11 shows **VX relative drift ~1e-6** on every outer side under constant I (non-ULP). Godunov mode has it asymmetric between y=0 (2.94e-07) and y=L (1.61e-06). This is a separate BC-implementation asymmetry, much smaller than the primary bugs. Lower priority; flag for follow-up.

## 8. Production fix decisions (waiting on user authorization before touching `godunov_flux` or `wave_operator.inl`)

### Fix A — interior non-fault flux symmetrization  *(high value, production-relevant)*

Two equivalent options:

**Option A1 — symmetrize inside `GodunovFlux::Interior`**: wrap the existing body with
```
Interior(n, L, R, F_h) {
  InteriorRaw(n, L, R, F_fwd);
  InteriorRaw(-n, R, L, F_rev);
  F_h = 0.5 * (F_fwd + F_rev);
}
```
Smallest surface (one file). Matches SeisSol's precomputed-face-Riemann symmetry.

**Option A2 — canonicalize `nor` + L/R swap at the call site**: change the interior branch in `ComputeADERFaceFluxRHS` (and the RK4 `ComputeFaceFluxRHS`) to pick an orbit-invariant `nor` (e.g. `(nor.dot(c_e2 − c_e1) > 0) ? nor : −nor`) and swap `I_self ↔ I_nbr` to match. Larger surface; touches two production paths.

Recommend A1.

### Fix B — boundary-face sign-symmetrization  *(secondary; SXZ-like components)*

Wrap `FreeSurfaceTotal` / `FreeSurfaceGodunovTotal` similarly:
```
FreeSurfaceTotalSymmetric(n, I_self, bg, F_h) {
  FreeSurfaceTotal( n, I_self, bg, F_fwd);
  FreeSurfaceTotal(-n, I_self, bg, F_rev);
  F_h = 0.5 * (F_fwd + F_rev);
}
```
Called from the boundary branches. Small surface.

### Fix C — test fixture D4 symmetry  *(test-side; required to remove the stubborn 1.00 drift from `test_adjacent_triangle_fault_uniformity`)*

Replace `BuildCartesianFaultMesh()`'s use of `MakeCartesian3D(..., TETRAHEDRON, ...)` with a hand-built tet mesh whose per-hex Kuhn diagonal orientations are D4-equivariant. Alternatively: relax the uniformity test's 1e-10 threshold on boundary-adjacent fault-adjacent tets, or gate only on the interior subset of fault-adjacent tets that have no outer-boundary face.

### Expected outcome after Fix A + Fix B on the Kuhn fixture

After both code fixes, the step-1 and step-2 drifts in `test_adjacent_triangle_fault_uniformity` should drop from their current magnitudes down to the Kuhn-fixture floor (~1.7e-3 to 1e-3 on SXY, essentially zero on SXZ after the boundary SXZ fix, unchanged SXY due to the fixture residual).

If the user needs the full 1e-10 pass, Fix C (or a relaxed threshold) is required.

## 9. Handoff

Errors are pinned:
- `H_IFACE_NOR` — real, fixable in `GodunovFlux::Interior` via symmetrization.
- `H_BFACE_NOR (SXZ)` — real, fixable analogously for the boundary BC calls.
- `H_KUHN_FIXTURE` — fixture-level, not a code bug; requires Fix C or threshold relaxation.

Recommended next decision: authorize Fix A + Fix B in production (subject to CLAUDE.md's "Extreme Care" policy on `godunov_flux` / `wave_operator.inl`; propose to land them behind a guarded toggle first), and separately decide whether Fix C is worth the effort.

Added Gates 5, 6, 7 to `test_adjacent_triangle_fault_first_step_audit.cpp`:

- **Gate 5** — per-face F_h orbit drift in the interior non-fault branch.
  - 5a raw `flux_.Interior(nor, L, R)`: max orbit F_h drift 2.000 on SXX (FAIL).
  - 5b symmetrized `0.5*(Interior(n, L, R) + Interior(−n, R, L))`: drift **5.68e-14** (PASS).
  - 5c canonical-normal only (no L↔R swap): drift 2.000 (FAIL; canonical normal alone is insufficient).
- **Gate 6** — per-face F_h orbit drift in the boundary branch (not physically meaningful as written; kept as a diagnostic only. Per-face `I_self` is orbit-sensitive because it depends on per-face shape-ordering, so raw F_h drift is not a clean test of branch correctness).
- **Gate 7** — element-level lifted orbit drift after M^{-1} on the fault-adjacent tets, under four modes:

| mode | iface SXY | iface SXZ | bface SXY | bface SXZ |
|---|---|---|---|---|
| raw baseline | 1.749 | 0.649 | 1.00 | 1.00 |
| interior symmetrized only | **5.98e-16** | **0** | 1.00 | 1.00 |
| boundary sign-symmetrized only | 1.749 | 0.649 | 1.00 | **8.30e-19** |
| both together | **5.98e-16** | **0** | 1.00 | **8.30e-19** |

**Findings**:

- **H_IFACE_NOR CONFIRMED.** The interior non-fault branch's orbit drift is fully eliminated by replacing `flux_.Interior(nor, I_self, I_nbr)` with its `(n, L↔R)` symmetrization in a test-side clone of the branch. Gate 7a's residual is 5.98e-16 — floating-point noise. Mechanism: MFEM's `nor = CalcOrtho(J_F)` orientation for an interior face depends on the Elem1 / Elem2 traversal order from the Kuhn split, and `flux_.Interior(nor, L, R)` is not invariant under `(nor → −nor, L↔R swap)` for its Rusanov dissipation term — the central-flux part is, but the dissipation has a sign that needs the simultaneous swap to cancel.
- **H_BFACE_NOR PARTIALLY CONFIRMED.** Sign-symmetrization `0.5*(Γ(n) + Γ(−n))` cleans the boundary **SXZ** drift (1.00 → 8.30e-19) but leaves the **SXY** drift at 1.00. The boundary branch has a second mechanism on SXY orthogonal to sign-symmetrization — the natural suspect is in-plane D4 rotational covariance of `FreeSurfaceTotal` within one side of the mesh (e.g. the x=0 vs x=L or z=0 vs z=L mirror pair, which differ by an in-plane rotation rather than an out-of-plane sign flip).

## 8. Next Step (next bisection)

For the interior non-fault branch the mechanism is pinned; the remaining design decision
is how to apply the fix in production. Two candidates:

- **Option A — Symmetrize in `flux_.Interior`**: modify `GodunovFlux::Interior(nor, L, R)` so it internally returns the sign-symmetrized mean. Lowest blast radius: a single `0.5*(...)` wrapper at the flux level. But CLAUDE.md flags `godunov_flux` / `wave_operator.inl` as Extreme Care; a production change needs a documented plan and explicit authorization.
- **Option B — Canonicalize `nor` per face pair before calling**: at assembly time, choose an orbit-invariant normal convention (e.g. vertex-lex ordering) and swap L/R to match. More invasive at the assembly loop level.

Option A is simpler and matches SeisSol's per-face Riemann convention (they build a
precomputed `fluxSolver` per face and then call it symmetrically from both sides); this
audit confirms that choice.

For the boundary branch, an additional diagnostic is required: dump per-face `F_h[SXY]`
across within-side orbits, grouped by boundary-attribute side (y=0 vs y=L vs x=0 vs
x=L vs z=0 vs z=L). If within a single side the drift is still 1.00, then the issue is
in-plane D4 rotational covariance of `FreeSurfaceTotal`. That's Test 3-bis below.

### Test 3-bis — `test_boundary_within_side_orbit_breakdown`  *(diagnostic)*

Extend the boundary-face per-face probe to split orbits by boundary side (six axis-
aligned sides on the outer cube) and print the within-side orbit drift. Expected:

- If drift-1.00 persists on single-side orbits (e.g., all z=0 faces alone) → bug is in-plane D4 rotational covariance of `FreeSurfaceTotal`.
- If single-side drift is ~0 but cross-side drift is 1.00 → bug is sign-anti-symmetry between normal orientations that differ by an in-plane rotation, a separable fix.

### Test 3-ter — swap `SetFreeSurfaceBCMode(Gamma)` for `Godunov` on the same fixture

Already partially exercised at Gate 6b. Gate 6b showed identical 1.876 drift in the
per-face diagnostic — but that's the uninformative probe. Repeat the lifted
Gate-7-style check with `Godunov` mode. If Godunov fixes bface SXY, the `Gamma` branch
is the root cause; if not, both BC modes share the bug.

## 9. Handoff

Interior non-fault branch: mechanism pinned (H_IFACE_NOR), Gate 7a green under
symmetrization. Ready to discuss fix direction (Option A vs Option B above) with the
user before touching production code.

Boundary branch: SXZ cleaned by sign-symmetrization; SXY residual stubborn at 1.00.
Propose to run Test 3-bis (within-side orbit breakdown) and Test 3-ter (Godunov vs
Gamma) next — both are pure audit extensions, no production touch — before committing
to a boundary-branch fix direction.
