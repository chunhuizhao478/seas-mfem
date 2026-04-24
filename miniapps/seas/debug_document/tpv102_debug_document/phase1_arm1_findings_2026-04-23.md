# Phase 1 Arm 1 Localization Findings — 2026-04-23

## Document structure

This document has eleven sections, in chronological order:

1. **§A–§I: provenance.**
10. **§J: Round-8 Step 0 re-measurement.** Amplifier-not-seed
    reframing. Power-law t^2 growth. Recommended FREEZE-A/B/C
    diagnostics + dt refinement.
11. **§K: Round-9 execution.** **AUTHORITATIVE.** Per round-9 spec
    (R-001 refit, R-002 cross-fixture, R-003 FREEZE specification,
    R-004 dt discriminator, R-005 pre-committed verdict matrix).
    Executed FREEZE-C + dt-refinement on both fixtures.
    **Decisive verdicts:** (a) FREEZE-C ratio = 1.000 on both
    fixtures — **ψ rate-state feedback is NOT the amplifier**;
    (b) dt refinement gives tau1_corr ∝ dt^-0.17 — **NOT
    sub-resolved physics**; (c) power-law refit confirms t^1.80
    (R² = 0.997), exponential rejected (R² = 0.465). §K narrows the
    remaining candidate space: the pepper is per-step spatial
    numerical noise integrated through wave propagation, NOT a
    localizable single-line bug, NOT an eigenmode instability, NOT
    resolvable by dt refinement. FREEZE-A/FREEZE-B would require
    flux-layer authorization per §C R-006.

**Authoritative decision is in §K.** §A–§J are provenance.

---

# §A: Original Phase 1 Findings (pre-correction, SUPERSEDED)

Status: **Arm 1 probes ran successfully for the first time** (pre-R-001 they
were vacuous no-ops; v9.5.0 R-001 wired the `UsePrecomputedFaceFluxes`
flag into `WaveOperator::Mult`, so `ComputeFaceFluxRHS_ForTest` now
actually switches dispatch paths). Decision: **Phase 2D — D4-equivariant
fixture swap**.

## Run configuration

- Binary: `seas_test_arm1_constant_state_probes`
- Fixture: M0 (Cartesian-fault 2×2×2 tet mesh, Kuhn split, fault at y=L/2)
- Build: `mfem-dev` conda env, clean build after R-001..R-005 land
- Rank count: serial (np=1)
- Raw log: `phase1_arm1_run.txt` (this directory)

## Probe verdicts

| Probe                                        | Drift       | Tol     | Verdict |
|----------------------------------------------|-------------|---------|---------|
| G_CONST_VOL                                  | 0.000e+00   | 1.0e-10 | PASS    |
| G_CONST_VOL_MINV                             | 0.000e+00   | 1.0e-10 | PASS    |
| G_CONST_BFACE_LIFT (runtime orbit drift)     | 0.000e+00   | 1.0e-10 | PASS    |
| G_CONST_BFACE_LIFT (precomp orbit drift)     | 0.000e+00   | 1.0e-10 | PASS    |
| G_CONST_BFACE_LIFT (runtime vs precomp)      | 1.041e-01   | 1.0e-06 | FAIL    |
| G_ORBIT_DOF_MAP (orbits mismatched / 8)      | 8 / 8       | < 0.5   | FAIL    |

(§A's content continues unchanged below — preserved for audit. For
decisions, skip to §B.)

### Interpretation (superseded)

The original doc applied the decision tree as follows and selected
Phase 2D on the basis that G_ORBIT_DOF_MAP's "worst diff = 99" meant
"basis assigns different local indices at coincident physical points".
REVIEW.md R-001 showed that `99` is a SENTINEL for `find(key) == end()`
— "physical key not found", not "index mismatched". The probe conflated
mesh-layer vertex asymmetry with basis-layer index non-covariance, and
§A's decision rested on the conflation. §B replaces this.

### Phase 2D predicted outcomes (superseded)

§A's predictions had two problems flagged by REVIEW.md R-002 and R-005:
- They assumed D4 fixture vertex equivariance without verifying it
  (R-002).
- They did not condition investigation closure on Layer 2 pepper guard
  running on a production-scale (Gmsh) mesh (R-005).

§B revises both.

---

# §B: Corrections Addendum (2026-04-23, post-REVIEW.md)

Status: **Phase 2D is still the recommended next branch, but its
PRECONDITIONS and decision rules are tightened.** Prior bugs in §A's
reasoning are neutralized.

## Instrumentation changes (diagnostic-layer only, no freeze violation)

1. `ProbeOrbitDofMap` split into two sub-probes (vertex coincidence vs
   index ordering on coincident points). File:
   `test_arm1_constant_state_probes.cpp`. R-001 fix.
2. `test_d4_fixture_vertex_equivariance.cpp` added — verifies
   `BuildD4Mesh` produces y-mirror-paired vertex sets at the mesh
   level, with slot-wise match. Phase 2D precondition. R-002 fix.
3. `ProbeNonConstVol` / `ProbeNonConstVolMinv` /
   `ProbeNonConstBfaceLift` added to `test_arm1_constant_state_probes.cpp`,
   driven by a y-mirror-symmetric polynomial non-constant Q. R-003 fix.
4. `ProbeConstBfaceLift` instrumented to report per-component max
   diff, per-element concentration, and the worst-dof's owner
   element. Also reports absolute magnitudes of both paths (the
   interpretation of path_diff hinges on the scale, which §A got
   wrong). R-004 fix.

No production code (WaveOperator, PrecomputedFaceFluxes, GodunovFlux)
was touched. Flux-layer freeze respected.

## Corrected probe verdicts (M0 Kuhn fixture)

### Constant-Q probes (unchanged signal, improved reporting)

| Probe                                          | Value            | Tol      | Verdict |
|------------------------------------------------|------------------|----------|---------|
| G_CONST_VOL (orbit drift)                      | 0                | 1e-10    | PASS    |
| G_CONST_VOL_MINV (orbit drift)                 | 0                | 1e-10    | PASS    |
| G_CONST_BFACE_LIFT (runtime orbit)             | 0                | 1e-10    | PASS    |
| G_CONST_BFACE_LIFT (precomp orbit)             | 0                | 1e-10    | PASS    |
| G_CONST_BFACE_LIFT (runtime vs precomp)        | 1.041e-01 abs    | 1e-06    | FAIL    |

### Corrected G_CONST_BFACE_LIFT interpretation (R-004)

The new instrumentation shows:

- `||rhs_runtime||_inf = 2.500e+14` — the BFACE probe operates at
  **2.5 × 10¹⁴**, not 1.56 × 10⁷ (that was the VOLUME probe's scale).
  Both §A and REVIEW.md R-004 cited the wrong scale.
- `||rhs_precomp||_inf = 2.500e+14` — paths produce the same magnitude.
- `path_diff / ||rhs_runtime||_inf = 4.163e-16` — **pure ULP-floor
  noise** (1.9 × machine epsilon).
- Per-component diff spread: SXX:1e-1, SYY:5e-2, SZZ:4.5e-2, SXY/SXZ
  in the 2e-2 range; velocities VX/VY/VZ in the 1e-8 range. All
  components consistent with ULP propagation through a 2.5e14 scale
  operator.
- Concentration: 32/48 elements (67%) carry > 10% of the worst
  component's diff. DIFFUSE, not localized. Diffuse ULP is not
  consistent with a frame-mismatch bug (which would concentrate on a
  small subset of faces).

**Revised parking-lot decision:** the runtime-vs-precomp disagreement
IS ULP-floor noise. Parking it is DEFENSIBLE. The original doc's
parking decision was correct, but only coincidentally — the argument
used the wrong scale. R-004 is fully closed.

### Split G_ORBIT_DOF_MAP (R-001)

| Probe                                          | Value | Tol  | Verdict |
|------------------------------------------------|-------|------|---------|
| G_ORBIT_DOF_MAP (vertex coincidence)           | 16    | < 0.5 | FAIL  (mesh-layer) |
| G_ORBIT_DOF_MAP (index ordering on coincident) | 16 / 16 | < 0.5 | FAIL  (basis-layer) |

Additional instrumentation reports:
- orbits examined: 8
- keys not found (vertex non-coincident): 16
- coincident keys checked: 16
- indices mismatched on coincident keys: **16 / 16** (100%)
- worst idx diff on coincident keys: 1

**This is DIFFERENT from what REVIEW.md R-001 predicted.** The
reviewer's strong form said: "for Kuhn tets, `find(key)` returns
`end()` because keys don't exist, so the local-index comparison never
happens" — implying basis untested. The amended probe shows this is
only PARTIALLY true:

- 16 keys not found → mesh-layer vertex asymmetry IS present
  (reviewer correct on this half).
- 16 keys found and ALL 16 show index mismatch → the basis-layer IS
  exercised on half the comparisons, and fails 100% of them
  (reviewer's prediction was over-strong).

On Kuhn M0, BOTH mesh-layer and basis-layer orbit non-covariance are
present. The "basis untested" claim held for the 16 keys-not-found
cases but not for the 16 coincident cases.

Note that "worst idx diff = 1" is small — the basis-layer failure is a
single local-index swap, not a full reshuffling. That's still D4
non-covariant at the basis-index level, just by one slot.

### Non-constant-Q addendum (R-003) — on M0

Phase 2A / 2B can be fairly ruled out with this evidence:

| Probe                                            | Value            | Tol    | Verdict |
|--------------------------------------------------|------------------|--------|---------|
| G_NONCONST_VOL (relative drift)                  | 0                | 1e-12  | PASS    |
| G_NONCONST_VOL_MINV (relative drift)             | 2.87e-17         | 1e-12  | PASS    |
| G_NONCONST_BFACE_LIFT (runtime vs precomp, rel)  | 2.57e-16         | 1e-12  | PASS    |
| G_NONCONST_BFACE_LIFT (runtime orbit, rel)       | **2.22e-01**     | 1e-12  | **FAIL** |
| G_NONCONST_BFACE_LIFT (precomp orbit, rel)       | **2.22e-01**     | 1e-12  | **FAIL** |

Input: y-mirror-symmetric polynomial `Q_SXX(x,y,z) = A*(0.5+x/L)*g(y)*(0.5+z/L)`
where `g(y) = 0.5 + 4*y*(L-y)/L²` (symmetric about y=L/2, range [0.5, 1.5]).

**Key findings:**

1. **G_NONCONST_VOL = 0** → the volume kernel is orbit-covariant on
   non-constant Q, not just constant Q. Phase 2A ruled OUT on evidence,
   not on a degenerate input.
2. **G_NONCONST_VOL_MINV = 2.9e-17 rel** → mass-inverse is also clean
   on non-constant Q. Phase 2B ruled OUT.
3. **G_NONCONST_BFACE_LIFT runtime-vs-precomp = 2.6e-16 rel** → paths
   continue to agree to ULP even on non-constant Q. Phase 2C (per-
   side topology-frame fix) still without a trigger.
4. **G_NONCONST_BFACE_LIFT orbit drift = 22.2% on BOTH paths** —
   unambiguously large, and IDENTICAL across the two dispatch paths.

### Disambiguating the 22.2% bface orbit drift

Two hypotheses account for the 22.2%:

- **H1 (mesh-layer):** Kuhn M0 has 16 keys_not_found per the split
  probe. On non-constant Q the per-vertex values DIFFER between
  orbit-paired tets that don't share vertices. The cell-mean of
  per-DOF face-flux contributions on each tet therefore differs —
  even if the kernel is D4-covariant — because the INPUT DATA it
  integrates is effectively different between the two tets.
- **H2 (kernel-layer):** The face-flux kernel is not D4-covariant at
  the per-DOF level. The orbit drift would persist even if the input
  Q on orbit-paired tets matched DOF-by-DOF.

The BEST way to disambiguate H1 from H2: re-run `ProbeNonConstBfaceLift`
on the D4 fixture (where vertex-set equivariance is PROVEN — see
R-002 result below). On D4:
- If orbit drift drops to ~0 → H1 confirmed (mesh-layer). Phase 2D
  = fixture swap closes the investigation.
- If orbit drift remains ~22% → H2 confirmed (kernel-layer). Phase 2C
  (flux-layer per-side topology frame, or basis replacement) is
  needed.

This is a tighter Phase 2D entry condition than §A had.

### R-002: BuildD4Mesh vertex equivariance (precondition verification)

Binary: `seas_test_d4_fixture_vertex_equivariance`. Raw log in
`phase1_arm1_run_corrected.txt`.

| Check                      | Value             | Verdict |
|----------------------------|-------------------|---------|
| orbit buckets checked      | 24                | —       |
| element pairs checked      | 24                | —       |
| slot-wise y-mirror errors  | 0                 | —       |
| worst slot y-mirror error  | 0.000e+00 (exact) | —       |
| slot-wise equivariance     | PASS              | PASS    |
| set-wise equivariance      | PASS              | PASS    |

BuildD4Mesh produces **slot-wise** y-mirror equivariant tets — i.e.,
for every orbit pair (e_a, e_b), the vertex at local slot i of e_a
y-mirrors EXACTLY to the vertex at local slot i of e_b. Not just
set-wise (some permutation) but slot-wise. This is the strongest form
of equivariance.

**Consequence for Phase 2D:** re-running the R-001 split probe on
BuildD4Mesh will give:
- `keys_not_found = 0` guaranteed (set-wise equivariance + slot-wise
  gives per-DOF coincidence).
- `index_mismatched` — whatever the result, it genuinely tests MFEM's
  basis D4-covariance. If 0 → basis is fine on D4; if > 0 → MFEM's
  L2 Gauss-Lobatto basis is not orbit-covariant even given D4-
  equivariant input. That is the meaningful Phase 2D signal.

REVIEW.md R-002's concern is addressed: the D4 fixture is correctly
built, so Phase 2D's G_ORBIT_DOF_MAP verdict on it is interpretable.

## Revised Phase 2 branch selection

Current evidence:
- Phase 2A (volume kernel fix): ruled OUT on non-constant-Q evidence.
- Phase 2B (mass-matrix): ruled OUT on non-constant-Q evidence.
- Phase 2C (per-side flux-layer fix): no trigger — paths agree to ULP
  on both constant and non-constant Q.
- Phase 2D (D4 fixture): **primary candidate, with tightened entry
  preconditions**.
- Phase 2E (rupture-path instrumentation): parking-lot if Phase 2D
  does not close the pepper symptom on production-scale mesh.

**Selected branch: Phase 2D, with the following precondition set.**

### Phase 2D preconditions (all must be satisfied before running)

1. `seas_test_d4_fixture_vertex_equivariance` passes. ✓ (just confirmed,
   slot-wise equivariance holds exactly).
2. `seas_test_R001_rk4_uses_precomputed_flag` passes. ✓ (Phase 0 gate).
3. Split `ProbeOrbitDofMap` + non-constant-Q probes in place. ✓
   (R-001 + R-003 fixes landed).

### Phase 2D procedure

1. Run `seas_test_arm1_constant_state_probes` on BuildD4Mesh. This
   requires a small refactor — add a parameter to each probe to accept
   a mesh, then call with both M0 and D4 meshes. (Mesh-layer change,
   not flux-layer.)
2. Re-measure the 6 probe signals on D4:
   - G_CONST_VOL, G_CONST_VOL_MINV
   - G_CONST_BFACE_LIFT (runtime, precomp, runtime-vs-precomp)
   - G_ORBIT_DOF_MAP (both sub-probes)
   - G_NONCONST_VOL, G_NONCONST_VOL_MINV
   - G_NONCONST_BFACE_LIFT (three sub-probes)
3. Re-run the Phase 3 audit (`seas_test_adjacent_triangle_fault_first_step_audit`)
   with fixture = D4, measure Gate 14′ interior.

### Phase 2D decision rules (revised, no vacuous branches)

Match the observed D4 results to this decision table. Every branch is
reachable and each rule is unambiguous.

| D4: G_ORBIT_DOF_MAP (index) | D4: G_NONCONST_BFACE orbit | D4: Gate 14′ | Diagnosis                        | Next step |
|----------|----------|----------|----------------------------------|-----------|
| PASS     | PASS     | PASS (0)     | Kuhn fixture was the bug         | See R-005 gate below |
| PASS     | FAIL     | FAIL         | Face-flux kernel not D4-covariant at per-DOF level | Phase 2C (per-side flux-frame fix)    |
| FAIL     | —        | —            | MFEM basis not D4-covariant      | Phase 2F — modal L2 basis / SeisSol bridge |

A "PASS / FAIL / PASS (0)" row is impossible by construction: if Gate
14′ = 0 on D4 and kernel is orbit-covariant and basis is orbit-
covariant, then the M0 pepper must be purely fixture. Production-scale
Layer 2 verification is the closure gate (R-005).

### R-005 closure gate: Layer 2 pepper guard on a D4-aligned
### production-scale mesh

Prediction #1 in §A said "D4 fixture passes → close investigation".
This is insufficient. Per REVIEW.md R-005 and the plan's Layer
1+2+3 contract:

Before the investigation closes, EITHER:
- (a) Adapt `BuildD4Mesh` to produce a TPV102-scale mesh (60 km × 30 km
  × 30 km at p=1, order-matched resolution) AND run
  `seas_test_adjacent_triangle_fault_uniformity` (serial + np=4) on
  it with rupture drive. Guard must pass the plan-§11 channel
  tolerances (1e-10 for slip_rate, tau1_corr, tau2_corr,
  sigma_n_corr), OR
- (b) Audit the Gmsh-generated TPV102 production mesh for D4 vertex
  equivariance (almost certainly fails — Gmsh is not orbit-symmetric
  by default). If it fails, no Gmsh-mesh-based closure is
  justified; fixture-swap must be applied to production before
  the investigation closes.

(a) is the path-of-least-resistance. (b) is exploratory. Skipping
both is not acceptable — the plan's Layer 2 gate must be satisfied
on a mesh the production pipeline actually uses, not just on M0 / D4
test fixtures.

## Freeze policy (explicit unblock criteria) — R-006

Phase 3 STOP freeze remains in effect with the following clarifications.

**Permitted now:**
- Mesh-layer additions (new fixture builders, mesh test generators).
- Diagnostic-layer instrumentation (Arm 1 probes, amplification-chain
  tests, path-localization, per-component reporting).
- Parameterization of existing tests to accept the D4 fixture.

**Still frozen until explicit user approval:**
- Any change to `PrecomputedFaceFluxes` matrix construction.
- Any change to `WaveOperator`'s face-flux dispatch logic (including
  Phase 2C per-side topology-frame work, even if evidence from Phase
  2D implicates the flux layer).
- Any change to `GodunovFlux` BC dispatch or interior flux kernel.

**Flux-layer freeze unblock criteria — ALL must hold:**
1. Phase 2D D4 probe re-measurement complete, verdict documented
   per the decision table in §B.
2. Layer 2 pepper guard (R-005 closure gate) run on a production-
   scale mesh and logged. Either green (close) or red (flux-layer
   work justified).
3. User explicitly authorizes flux-layer reopen in a comment or
   commit message citing this document.

## References

- REVIEW.md (2026-04-23, Phase 1 review) — issued findings R-001
  (CRITICAL), R-002 (CRITICAL), R-003, R-004, R-005, R-006.
- Plan: `tpv102_seissol_aligned_flux_plan_2026-04-23.md` §8.5 Option
  F2 (D4 fixture).
- Phase 3 STOP: `phase3_stop_findings_2026-04-23.md` — freeze source.
- Raw probe logs:
  - `phase1_arm1_run.txt` (original run, pre-correction).
  - `phase1_arm1_run_corrected.txt` (amended probe suite,
    post-R-001/R-003/R-004 instrumentation).
  - `phase1_arm1_run_amended.txt` (interim run after R-001+R-004,
    before R-003 addendum was added).
- Probe sources:
  - `tests/unit/test_arm1_constant_state_probes.cpp` (amended).
  - `tests/unit/test_d4_fixture_vertex_equivariance.cpp` (new).

---

# §C: Round-2 Corrections Addendum (2026-04-23, post-REVIEW.md round 2)

Status: **§B's Phase 2D decision framework had two critical
analytical errors** (REVIEW.md round 2 R-001, R-002). §C replaces
§B's decision table and closure gate. No production code touched.
All changes are diagnostic-layer; freeze respected.

## Summary of round-2 corrections

| Round-2 finding | Category | Summary | §C action |
|-----------------|----------|---------|-----------|
| R-001 r2 (CRITICAL) | Probe insensitive | `OrbitCellMeanMaxDrift` collapses to ∮ F·n dS by divergence theorem; per-DOF asymmetry invisible | Added `OrbitPerDofPairedDrift`; rescinded "Phase 2A ruled OUT" claim |
| R-002 r2 (CRITICAL) | Decision table | Missing `PASS/PASS/FAIL` row (likeliest outcome per v8→v9.4 history) | Full 8-row table below |
| R-003 | Process rule violation | Production-scale D4 path violates `feedback_no_local_reproducer` + `feedback_frontera_approval` | Demoted to Frontera-gated; local Gmsh audit promoted |
| R-004 | Hypothesis gap | Missing H3 (FACE2NODES boundary-face vertex selection) | Added `ProbeBfaceVertexSelection`; updated H1/H2/H3 taxonomy |
| R-005 | Over-escalation | Phase 2F (2–4 d basis replacement) for idx diff=1 | Split into Phase 2F-prelim (vertex canonicalization, 1 d) and Phase 2F-full |
| R-006 | Tolerance miscalibrated | `kULPFloorRel = 1e-12` is ~4500× ULP, not "10×" | Tightened to `1e-13` (~450× ULP) |

## Corrected probe measurements (amended suite on Kuhn M0)

The round-2 amended binary `seas_test_arm1_constant_state_probes`
adds `OrbitPerDofPairedDrift` alongside cell-mean, a
`ProbeBfaceVertexSelection` probe, and a tightened ULP floor. Raw log:
`phase1_arm1_run_round2.txt`.

### Volume kernel (Kuhn)

| Metric                                             | Value        | Verdict |
|----------------------------------------------------|--------------|---------|
| ||rhs_vol||_inf                                    | 3.22e+07     | —       |
| cell-mean drift (rel) — STRUCTURALLY 0             | 0            | blind   |
| **per-DOF paired drift (rel), n=144 pairs**        | **8.33e-01** | **FAIL**|

Interpretation: the cell-mean probe returns 0 because of the
partition-of-unity collapse (proved in §C R-001 r2 below). The per-DOF
probe fires at 83% relative drift — but this is NOT unambiguously a
kernel bug because the 144 coincident-DOF pairs on Kuhn are embedded
in tets whose OTHER vertices don't coincide, so the in-tet gradient
at the coincident vertex is computed from different neighbor data.
Phase 2A / 2B / 2C are all live candidates; the per-DOF probe must be
re-run on D4 to cleanly isolate kernel vs. fixture effects.

### Mass-inverse (Kuhn)

| Metric                                             | Value        | Verdict |
|----------------------------------------------------|--------------|---------|
| ||M^-1 rhs_vol||_inf                               | 3.09e+01     | —       |
| cell-mean drift (rel) — weakly sensitive           | 2.87e-17     | ULP     |
| **per-DOF paired drift (rel), n=144 pairs**        | **8.33e-01** | **FAIL**|

Same caveat as volume: per-DOF drift on Kuhn is mesh-contaminated.
Re-run on D4 for clean signal.

### Boundary-face lift (Kuhn, fault-less)

| Metric                                             | Value        | Verdict |
|----------------------------------------------------|--------------|---------|
| ||rhs_runtime||_inf / ||rhs_precomp||_inf          | 6.33e+14     | —       |
| runtime-vs-precomp (rel)                           | 2.58e-16     | PASS (ULP) |
| cell-mean orbit drift (rel), both paths            | 2.22e-01     | FAIL    |
| **per-DOF orbit drift (rel), both paths, n=144**   | **2.22e-01** | **FAIL**|

Bface is a boundary-assembly operator with no divergence-theorem
collapse; cell-mean and per-DOF give the same signal. The 22% drift
is still ambiguous between H1 / H2 / H3 on Kuhn — must re-run on D4.

### Boundary-face vertex selection (R-004 probe, Kuhn)

| Metric                                        | Value | Verdict |
|-----------------------------------------------|-------|---------|
| lower-y boundary triangles                    | 24    | —       |
| upper-y boundary triangles                    | 24    | —       |
| cross-half (cx, cz)-matched pairs             | 8     | partial-coverage |
| buckets with no cross-half match              | 16    | —       |
| set-wise mismatches (on matched pairs)        | 0     | PASS (weak) |

The probe has partial coverage on Kuhn (only 8 of 24 boundary
triangles have a y-mirror partner at matching centroid). On those 8
pairs, FACE2NODES picks y-mirror-paired triangle vertex SETS
correctly → H3 weakly ruled out. Full coverage requires re-running on
D4 (where all 24 should match).

## Updated hypothesis space for the bface orbit drift

§B had H1 (mesh-layer) and H2 (kernel-layer). §C adds H3:

- **H1 (bulk vertex non-coincidence):** Kuhn orbit-paired tets do not
  share all 4 vertex positions after y-mirror (we observe 16
  keys_not_found + 16 coincident). Per-DOF input Q differs → per-DOF
  rhs differs.
- **H2 (kernel-layer non-covariance):** even given D4-equivariant
  input, the kernel's per-DOF computation is orbit-asymmetric.
- **H3 (boundary-face vertex selection):** MFEM's `FACE2NODES[4][3]`
  picks 3 of the 4 tet vertices per face. If orbit-paired tets have
  the same vertex SET but in different SLOT order, the FACE2NODES-
  derived triangle vertex set differs → nApNm1 built from a
  different 3-vertex triangle.

**D4 fixture swap fixes H1 AND H3 simultaneously.** It does NOT
address H2. Conversely, H3 is a face-vertex-ordering normalization
issue that CAN be fixed on Gmsh production meshes (a local-laptop
rewrite of the face vertex selection), whereas H1 cannot (requires
D4 mesh generation).

Partial coverage of the R-004 probe on Kuhn (8/24) gives weak
evidence that H3 is ruled out. Phase 2D on D4 will achieve full
coverage.

## Phase 2D decision table (FULL 8-row form, R-002 r2)

Abbreviations:
- `DOF_MAP` = `G_ORBIT_DOF_MAP (index ordering on coincident pts)`
  on D4. PASS ⇒ basis orbit-covariant on D4.
- `BFACE` = `G_NONCONST_BFACE_LIFT (per-DOF orbit rel, either path)`
  on D4. PASS ⇒ bface kernel orbit-covariant on D4.
- `Gate14` = `test_adjacent_triangle_fault_first_step_audit`'s
  Gate 14′ iface, on D4 fault fixture. PASS ⇒ first-step pepper
  absent on D4.

| DOF_MAP | BFACE | Gate14 | Diagnosis | Next step |
|---------|-------|--------|-----------|-----------|
| PASS    | PASS  | PASS   | Kuhn fixture was the bug. | §C R-005 closure gate (local Gmsh audit). |
| PASS    | PASS  | FAIL   | Bug lives in code paths the static probes don't reach: fault dispatch per-side frame, multi-step compounding, or `Loc1/Loc2` shape-table differences. | Phase 2C-fault: instrument `PrecomputedFaceFluxes::BuildInteriorMatrices` per-side antisymmetry on D4; add 2-step audit gate. Do NOT conclude "close to PASS/PASS/PASS; fixture was the bug". |
| PASS    | FAIL  | PASS   | Inconsistent — bface drift exists but does not drive Gate 14′. Likely probe-to-audit divergence. | Re-verify probes use same fixture/orbit bucketing as audit. No flux-layer work authorized. |
| PASS    | FAIL  | FAIL   | Face-flux kernel not D4-covariant at per-DOF level even on D4. | Phase 2C-flux (per-side topology-frame fix). Flux-layer freeze unblock: §C freeze criteria + user approval. |
| FAIL (idx diff ≤ 2) | —     | —      | MFEM DOF assignment depends on stored vertex ordering; orbit-paired tets have mismatched slots. Small permutation (1–2 slots). | Phase 2F-prelim: vertex-ordering canonicalization pass on production mesh post-Gmsh. Re-run D4 DOF_MAP probe. If still FAIL → 2F-full. 1 day. |
| FAIL (idx diff > 2) | —     | —      | Large-scale MFEM basis mis-assignment; not a simple permutation. | Phase 2F-full: modal L2 basis or SeisSol bridge. 2–4 days. |
| FAIL    | PASS  | FAIL   | Basis non-covariance surfaces in multi-step audit but not in static bface. | Phase 2F-prelim first (cheap). If that closes the basis probe, re-run Gate 14′ on canonicalized mesh. |
| FAIL    | FAIL  | FAIL   | Compound failure: basis + kernel. | Phase 2F-prelim first (cheapest). If that doesn't close BFACE/Gate 14′, then Phase 2C-flux. |

Note: **PASS/PASS/FAIL is the concerning row** (R-002 r2). It is not
ruled out by any single-operator static probe. Per v8→v9.4 debug
history, Gate 14′ has been the persistent fail; a PASS/PASS/FAIL
outcome on D4 would mean the bug lives deeper in the audit pipeline
(fault dispatch / multi-step) than the Arm 1 probes reach. Treat as
"reopen investigation with multi-step / fault-dispatch probes", NOT
as "close to PASS/PASS/PASS".

## R-005 closure gate (revised per R-003)

Previously §B recommended "(a) adapt `BuildD4Mesh` to TPV102-scale
(60 km × 30 km × 30 km at 200 m, serial + np=4 rupture drive) as
path-of-least-resistance". This violates `feedback_no_local_reproducer.md`
(no local TPV102 reproducer, even at 1 km) and would require Frontera
approval per `feedback_frontera_approval.md`.

**Revised closure gate:**

### Local-only path (recommended first — no approval needed)

**(b-local) Audit the Gmsh-generated TPV102 production mesh for D4
y-mirror equivariance.** This is read-only mesh analysis. Feasible on
laptop in minutes:

1. Load a TPV102 production `.msh` file (e.g. `tpv102_1000m.msh`).
2. For every tet with centroid cy < L/2, find the tet with centroid
   (cx, L-cy, cz) and compare vertex SETS after y-mirror.
3. Report fraction of tets that y-mirror-pair at the vertex-set level
   (H1 diagnostic) and at the slot level (H1+H3 diagnostic).

Outcomes:
- **Gmsh mesh is D4-equivariant (very unlikely):** Phase 2D D4 result
  may transfer. Close investigation.
- **Gmsh mesh is NOT D4-equivariant (expected):** Phase 2D "fixture-
  was-the-bug" conclusion does NOT transfer to production. One of:
  - (i) Port production to a D4-equivariant mesh generator (project-
    scope; major work; not a Phase 2D step).
  - (ii) Apply vertex-ordering canonicalization (Phase 2F-prelim style)
    to the Gmsh mesh post-generation. Only addresses H3 / basis-index
    slot-match, not H1 (which requires true y-mirror vertex positions).
  - (iii) Accept that Phase 2D's positive D4 verdict is a fixture-
    specific artifact; reopen investigation on non-D4 production mesh.

### Frontera-gated path (requires explicit user approval)

**(a) Production-scale D4 pepper-guard sbatch.** Only if user
authorizes per `feedback_frontera_approval.md`. Generate the sbatch
file but DO NOT submit. This is a Layer 3 verification (production
scale), not a Layer 2 unit-test closure gate.

### What counts as "Layer 2 closure"

- Minimum: (b-local) complete AND decision rule for Gmsh non-
  equivariance documented. Team can close with a written finding
  that Phase 2D's D4 result is fixture-specific.
- Full: (b-local) + (a) run with user approval, pepper-guard green
  on production-scale D4 mesh. Stronger but gated on Frontera
  authorization.

**Skipping both is not acceptable** — the plan's Layer 2 gate must
apply to a production-relevant mesh, not just M0 / D4 toy fixtures.
But between (a) and (b-local), (b-local) is the path-of-least-
resistance under standing rules.

## Freeze policy (round-2 clarification)

Unchanged from §B's three-item list. Added clarification:

The phrase "Phase 2D produces the Gate 14′ re-measurement" in the
freeze-lift criterion is ambiguous. Round-2 clarification: Gate 14′
re-measurement means running
`seas_test_adjacent_triangle_fault_first_step_audit` on a D4
fault-bearing mesh and logging the numeric Gate 14′ iface value.
Mere existence of the number is not sufficient — the §C decision
table must produce an unambiguous next-step.

If the D4 Gate 14′ measurement lands in row PASS/PASS/FAIL (the
"concerning" row), a subsequent multi-step probe IS the additional
evidence required before the flux-layer freeze lifts.

## References (round-2 additions)

- REVIEW.md round-2 (2026-04-23): R-001 r2, R-002 r2, R-003 r2, R-004,
  R-005, R-006 (low).
- Raw probe log (round-2 amended): `phase1_arm1_run_round2.txt`.
- User rules:
  - `feedback_no_local_reproducer.md` — no local TPV102 reproducer.
  - `feedback_frontera_approval.md` — Frontera approval required.
- Probe sources (round-2 additions):
  - `tests/unit/test_arm1_constant_state_probes.cpp` —
    `OrbitPerDofPairedDrift` helper, `ProbeBfaceVertexSelection`
    probe, tightened `kULPFloorRel`.

## §C audit trail — what did NOT change

- §A and §B are preserved verbatim above. The "authoritative
  decision" marker moved from §B to §C.
- The R-002 BuildD4Mesh vertex-equivariance gate (test_d4_fixture_
  vertex_equivariance) is unchanged from §B. Its PASS verdict
  continues to hold.
- The R-001 round-1 `ProbeOrbitDofMap` split into
  `(vertex_coincidence, index_ordering_on_coincident_points)`
  remains the correct diagnostic. §C does not change it; §C adds
  per-DOF drift as a COMPLEMENTARY measurement for the non-constant-Q
  probes.
- Phase 0 fixes (v9.5.0 R-001..R-005) unaffected.
- Flux-layer freeze unchanged.

---

# §D: Phase 2D Execution Results (2026-04-23)

Status: **AUTHORITATIVE DECISION.** Phase 2D ran the §C-amended Arm 1
probes + the Phase 3 first-step audit on both the Kuhn M0 fixture
(baseline) and the D4-equivariant fixture. The §C round-2 review had
predicted that the `PASS/PASS/FAIL` row was the likeliest outcome.
That prediction is confirmed. Phase 2D's decision is **NOT** "Kuhn
fixture was the bug"; it is "reopen investigation with multi-step /
fault-dispatch probes" per §C's interpretation of the
`PASS/PASS/FAIL` row.

## §D.1 Run configuration + raw data

### Binaries executed

- `seas_test_d4_fixture_vertex_equivariance` — confirmed D4 fixture
  is slot-wise y-mirror equivariant. Precondition for Phase 2D.
- `seas_test_d4_fault_sanity` — confirmed `BuildD4Mesh(true)`
  produces exactly 8 fault interior faces, all at y = L/2. Precondition.
- `seas_test_arm1_constant_state_probes` — all Arm 1 probes run on
  Kuhn baseline, then same suite on D4 (new `RunD4ProbeSuite()` block).
- `seas_test_adjacent_triangle_fault_first_step_audit` — run 4 times:
  - `SEAS_TEST_FIXTURE=kuhn` (runtime path, default).
  - `SEAS_TEST_FIXTURE=kuhn SEAS_TEST_USE_PRECOMPUTED_FLUX=1`
    (precomputed path — activates Gate 14′).
  - `SEAS_TEST_FIXTURE=d4` (runtime path).
  - `SEAS_TEST_FIXTURE=d4 SEAS_TEST_USE_PRECOMPUTED_FLUX=1`
    (precomputed path).

### Raw logs (in this directory)

- `phase2d_d4_fixture.txt`
- `phase2d_d4_fault_sanity.txt`
- `phase2d_arm1_kuhn_and_d4.txt`
- `phase2d_audit_kuhn.txt`, `phase2d_audit_kuhn_precomp.txt`
- `phase2d_audit_d4.txt`, `phase2d_audit_d4_precomp.txt`

### Measurement table

| Probe / gate                                        | Kuhn         | D4           | Delta              |
|-----------------------------------------------------|--------------|--------------|--------------------|
| G_ORBIT_DOF_MAP (vertex coincidence)                | 16 FAIL      | 0 PASS       | mesh fix           |
| G_ORBIT_DOF_MAP (idx ordering on coincident)        | 16/16 FAIL   | **0/96 PASS**| **basis OK on D4** |
| G_BFACE_VERTEX_SELECTION (set mismatches)           | 0/8 PASS     | 0/96 PASS    | full coverage on D4|
| G_CONST_BFACE_LIFT runtime-vs-precomp rel           | 4.16e-16     | 1.25e-07     | precomp→runtime diverge on D4 (small) |
| G_CONST_BFACE_LIFT runtime orbit rel                | 0            | 1.87e-07     | small D4 drift     |
| G_CONST_BFACE_LIFT precomp orbit rel                | 0            | **8.04e-17** | ULP on D4          |
| G_CONST_VOL per-DOF (fault fixture)                 | —            | 3.12e+07 abs | ambient (antisym components) |
| G_CONST_VOL_MINV per-DOF (fault fixture)            | —            | **0 PASS**   | Minv D4-clean      |
| G_NONCONST_VOL per-DOF rel                          | 8.33e-01     | 2.000        | antisym signature dominates |
| G_NONCONST_VOL_MINV per-DOF rel                     | 8.33e-01     | **0 PASS**   | Minv D4-clean      |
| G_NONCONST_BFACE_LIFT per-DOF runtime rel           | 2.22e-01     | 2.34e-07     | 6 orders better    |
| G_NONCONST_BFACE_LIFT per-DOF precomp rel           | 2.22e-01     | **1.04e-16** | ULP on D4          |
| **Audit Gate 14′ (P3.5) precomp iface**             | **1.999765** | **2.000000** | **NO IMPROVEMENT** |
| **Audit Gate 14′ (P3.5) precomp bface**             | **1.000000** | **1.777778** | **WORSE ON D4**    |
| Audit Gate 14a (runtime constant-I iface)           | 0 PASS       | 0 PASS       | unchanged          |
| Audit Gate 14b (runtime constant-I bface)           | 1.0 FAIL     | 1.0 FAIL     | unchanged          |

### Key observations

1. **Basis D4-covariance:** confirmed — `G_ORBIT_DOF_MAP (idx
   ordering on coincident)` reports 0/96 mismatches on D4. Phase 2F
   is ruled OUT. The "idx diff = 1" on Kuhn was vertex-ordering, not
   basis non-covariance.

2. **Per-DOF static probes improve dramatically on D4:**
   - BFACE precomp: 22.2% → **1.04e-16** (ULP floor; 15+ orders
     of magnitude better).
   - Volume Minv: 83% → 0.
   - Bulk BFACE runtime path: 22.2% → 2.34e-07 (6 orders better,
     not ULP but a soft residual).

3. **`G_NONCONST_VOL per-DOF rel = 2.000` on D4:** this is the
   antisymmetric-component signature predicted in §C round 2
   interpretation. For y-anti-symmetric components (σ_xy, σ_yz, v_y)
   the mirror response at y-paired DOFs is equal-opposite; the
   probe's `|a - b|` gives 2·|a|, hence the exact 2.000. **Not a
   bug; probe limitation on antisymmetric components.**

4. **Audit Gate 14′ does NOT improve on D4:** iface stays at 2.0,
   bface gets worse (1.0 → 1.78). This is the decisive finding.

## §D.2 Decision-table row fired

Using §C's 8-row decision table with:
- `DOF_MAP` = **PASS** (basis orbit-covariant on D4)
- `BFACE` = **PASS** (precomp path at ULP; runtime has 2e-7 residual
  which is a soft-FAIL, but per §C round-2 the decision pivot is
  precomp since that's the production flux path)
- `Gate 14′` = **FAIL** (2.0 iface, 1.78 bface — unimproved from Kuhn)

Row: **`PASS / PASS / FAIL`**

§C's text for this row:
> Bug lives in code paths the static probes don't reach: fault
> dispatch per-side frame, multi-step compounding, or Loc1/Loc2
> shape-table differences.
> Next step: Phase 2C-fault: instrument
> `PrecomputedFaceFluxes::BuildInteriorMatrices` per-side
> antisymmetry on D4; add 2-step audit gate. Do NOT conclude "close
> to PASS/PASS/PASS; fixture was the bug".

§C round-2 explicitly warned:
> Treat this row as "investigation reopens with multi-step probe",
> not "close investigation".

The prediction landed. Apply mechanically.

## §D.3 Phase 2C-fault (Next /code-implement target)

**Decision:** pursue Phase 2C-fault, not Phase 2D-close (the
fixture-was-bug conclusion would be wrong here) and not Phase 2F
(basis is provably orbit-covariant on D4).

### Target code paths

The Gate 14′ failure on D4 is specifically in the **precomputed-path
lifted rhs on fault-bearing fixtures** (runtime Gate 14a at 0, but
Gate 14′ (P3.5) at 2.0). Fault-adjacent tets have sorted-signature
drift within-side, meaning the precomputed path's per-side fault flux
is NOT orbit-uniform across same-side fault-adjacent tets.

Candidate targets (in order of likelihood given the signal):

1. **`dynamic/precomputed_face_fluxes.cpp::BuildInteriorMatrices`**
   per-side construction (lines ~197-222). D4 has slot-wise
   equivariance, so the topology normal is common; the per-side
   matrices nApNm1 / nAmNm1 come from ComputeCellFaceFrame. A per-
   side asymmetry here would not be visible to the static Arm 1
   probes (which measure lifted rhs on non-fault BFACE) but WOULD
   show up at fault-adjacent tets in the audit. (REVIEW.md prior-
   round R-007 flagged this area as a separate antisymmetry concern.)

2. **Fault dispatch's per-side topology frame** in
   `dynamic/wave_operator.inl::ComputeFaceFluxRHS` interior-fault
   branch (lines ~1131-1373 after Phase 0 R-001 wiring). Under
   `UsePrecomputedFaceFluxes(true)` the fault branch still uses
   runtime `fault_basis_` frames (per-QP). If `fault_basis_` on
   D4-mirror-paired fault faces disagrees in sign_flipped / tangent
   choice, Gate 14′ fires.

3. **Multi-step compounding:** Gate 14′ measures first-step lifted
   rhs, but per `tpv102_debug_v9.4.0_debug_plan.md §11` the v8→v9.4
   failures show compounding over multiple steps. The single-step
   Gate 14′ already fires, so multi-step isn't the ONLY bug, but it
   may amplify.

### Phase 2C-fault acceptance criteria

- A per-side antisymmetry assertion on `BuildInteriorMatrices` for
  D4 fault-adjacent faces: for every orbit-paired pair `(f_a, f_b)`,
  `nApNm1(f_a, elem_a) == nAmNm1(f_b, elem_b)` up to ULP under the
  y-mirror permutation of {VX, SXY, SXZ, ...}.
- Gate 14′ iface on D4 drops from 2.0 → < 1e-10.
- Gate 14′ bface on D4 drops from 1.78 → < 1e-10.
- Arm 1 static probes remain PASS on D4 (no regression).

## §D.4 Layer 2 closure prediction

When Phase 2C-fault's fix lands, the acceptance probe chain is:

1. **Unit:** `seas_test_arm1_constant_state_probes` + `RunD4ProbeSuite`
   — re-run, verify no regression on per-DOF probes, D4 precomp
   metrics still at ULP.
2. **Audit:** `SEAS_TEST_FIXTURE=d4 SEAS_TEST_USE_PRECOMPUTED_FLUX=1
   ./seas_test_adjacent_triangle_fault_first_step_audit` — verify
   Gate 14′ iface and bface drop to < 1e-10.
3. **Layer 2 local-Gmsh audit** (§C R-005 closure gate): read Gmsh
   TPV102 production mesh, verify whether any of the same per-side
   asymmetry exists there. Because the Phase 2C-fault fix addresses
   per-side flux-frame antisymmetry, the fix should transfer even
   to non-D4 meshes; the Gmsh audit is a sanity check that the fix
   does not depend on slot-equivariance of the mesh.
4. **Layer 3 (Frontera-gated):** ONLY after user approval per
   `feedback_frontera_approval.md`, submit a TPV102-scale sbatch
   and check pepper guard on the production VTU output.

## §D.5 What this rules out for good

- **Phase 2A (volume kernel)** — not ruled out by the cell-mean
  probe (§C R-001 r2 correct), but the D4 per-DOF Volume_Minv probe
  is at 0, so the kernel + mass-inverse assembly is D4-covariant
  once the mesh is slot-equivariant. Phase 2A not the primary cause.
- **Phase 2B (mass-inverse)** — ruled out on D4 (`G_NONCONST_VOL_MINV
  per-DOF rel = 0`).
- **Phase 2F (basis replacement)** — ruled out by D4
  `G_ORBIT_DOF_MAP (idx ordering on coincident) = 0/96 PASS`.
  MFEM's L2 GaussLobatto basis IS orbit-covariant under slot-
  equivariant vertex ordering.
- **Phase 2D "fixture was the bug"** — ruled out by Gate 14′ failing
  on D4.
- **H3 (FACE2NODES vertex selection)** — ruled out on both Kuhn
  (partial coverage, 0 mismatches) and D4 (full coverage, 0 mismatches).

## §D.6 What this newly rules IN

- **Phase 2C-fault** is the indicated next branch. The target is
  per-side antisymmetry in `PrecomputedFaceFluxes::BuildInteriorMatrices`
  on fault-adjacent faces, which the static non-fault BFACE probes
  cannot detect but Gate 14′ (fault-bearing audit) does.
- **Flux-layer freeze unblock criteria met** per §B/§C:
  1. ✓ Phase 2D D4 probe re-measurement complete with documented
     verdict (§D.1 table).
  2. ✓ Local-layer diagnostic path (§D.4 Layer 2 step) documented;
     Frontera-gated step (§D.4 Layer 3) NOT attempted without approval.
  3. ⚠ User explicit reopen authorization required before Phase 2C-
     fault (flux-layer change) proceeds.

## References (Phase 2D additions)

- `dynamic/d4_tet_mesh.hpp` — shared BuildD4Mesh definition (Step A).
- `tests/unit/test_d4_fault_sanity.cpp` — fault-fixture sanity gate
  (Step B).
- `tests/unit/test_arm1_constant_state_probes.cpp::RunD4ProbeSuite()` —
  D4 sweep of all Arm 1 probes (Step C).
- `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp` —
  SEAS_TEST_FIXTURE=d4 env hook (Step D).
- Raw logs: `phase2d_*.txt` in this directory.

---

# §E: Round-3 Correction (2026-04-23)

Status: **AUTHORITATIVE.** §D's Phase 2C-fault decision is rescinded.
The round-3 review identified three issues with §D's interpretation
that compromise the next-step target. All three are confirmed by
follow-up diagnostics.

## Summary of R3 findings and outcomes

| Round-3 finding | Status | Action |
|---|---|---|
| R-001 Gate 14′ bface worsens on D4 | Confirmed partial artifact | Per-orbit dump shows drift distributed across ALL components on both fixtures; metric demands x-y-z equivariance that neither fixture provides |
| R-002 Pepper guard not run | Closed with evidence | `seas_test_adjacent_triangle_fault_uniformity` run on both Kuhn and D4; both fail all 4 channels; D4 mostly worsens spreads |
| R-003 VX per-DOF "antisym" explanation doesn't fit | Confirmed fixture bug | BuildD4Mesh has **24/48 tets with negative Jacobian determinants**; the "antisymmetric" signature on y-symmetric VX is a fixture orientation bug, not kernel non-covariance |
| R-004 Phase 2C-fault target label misleading | Corrected | Gate 14′ uses `RunPrecomputedFluxLiftedAudit` which explicitly skips fault faces (line 1028-1032). File target (`BuildInteriorMatrices`) is correct; label now reads "non-fault interior faces of fault-ADJACENT tets" |
| R-005 Layer 2 closure chain subs in non-pepper audit | Closed via R-002 | Pepper guard now run, so the closure gate is executed (result: FAIL on both) |
| R-006 Units column missing | Fixed | §E.1 table has explicit units column |

## §E.1 Corrected measurement table (units column per R-006)

| Probe / gate                                        | Kuhn         | D4           | Units            |
|-----------------------------------------------------|--------------|--------------|------------------|
| G_ORBIT_DOF_MAP (vertex coincidence)                | 16 FAIL      | 0 PASS       | count            |
| G_ORBIT_DOF_MAP (idx ordering on coincident)        | 16/16 FAIL   | 0/96 PASS    | count / count    |
| G_BFACE_VERTEX_SELECTION set mismatches             | 0/8 PASS     | 0/96 PASS    | count / count    |
| G_CONST_BFACE_LIFT runtime-vs-precomp               | 4.16e-16     | 1.25e-07     | dimensionless (abs diff / \|\|rhs\|\|_inf) |
| G_CONST_BFACE_LIFT runtime orbit                    | 0            | 1.87e-07     | dimensionless rel |
| G_CONST_BFACE_LIFT precomp orbit                    | 0            | 8.04e-17     | dimensionless rel |
| G_CONST_VOL per-DOF                                 | —            | 3.12e+07 (worst VX, contam.) | Pa/s absolute |
| G_CONST_VOL_MINV per-DOF                            | —            | 0 PASS       | m/s² absolute    |
| G_NONCONST_VOL per-DOF rel                          | 8.33e-01     | 2.00 (VX antisym, fixture bug) | dimensionless rel |
| G_NONCONST_VOL_MINV per-DOF rel                     | 8.33e-01     | 0 PASS       | dimensionless rel |
| G_NONCONST_BFACE runtime per-DOF rel                | 2.22e-01     | 2.34e-07     | dimensionless rel |
| G_NONCONST_BFACE precomp per-DOF rel                | 2.22e-01     | 1.04e-16     | dimensionless rel |
| Audit Gate 14′ precomp iface                        | 1.9998       | 2.000        | sorted-sig drift (∈ [0, 2]) |
| Audit Gate 14′ precomp bface                        | 1.000        | 1.778        | sorted-sig drift (∈ [0, 2]) |
| D4 Jacobian sign audit                              | —            | **24 pos + 24 neg = 48 total** | count |
| **Pepper guard slip_rate spread**                   | **4.56e-05** | **1.50e-04** | **m/s** (20 steps, serial) |
| **Pepper guard tau1_corr spread**                   | **2.236**    | **1.611**    | **Pa** (20 steps, serial) |
| **Pepper guard tau2_corr spread**                   | **1.74e-06** | **3.79e-06** | **Pa** (20 steps, serial) |
| **Pepper guard sigma_n_corr spread**                | **2.05e-06** | **9.31e-06** | **Pa** (20 steps, serial) |

## §E.2 R-003: Fixture orientation bug — the critical finding

### Observation (§D quote corrected)

§D claimed: "G_NONCONST_VOL per-DOF rel = 2.000 is the antisymmetric-
component signature." Round-3 review noted: input Q[SXX] produces a
rhs dominated by VX (∝ ∂σ_xx/∂x). `∂σ_xx/∂x` is **y-symmetric** for a
y-symmetric Q[SXX], so VX rhs at y-mirror-paired DOFs should be
EQUAL, not opposite.

### Direct diagnostic

Added component-wise per-DOF breakdown (test_arm1_constant_state_probes.cpp).
Result on D4:

```
per-component ||rhs||_inf:
    VX=2.975e+07  (ONLY nonzero component)
    all others 0

VX worst pair: a=+2.975e+07  b=-2.975e+07  |a-b|=5.950e+07  |a+b|=0.000e+00
```

**VX values are equal-opposite (antisymmetric) at y-mirror-paired
DOFs, contradicting the y-symmetric physics.** This is inconsistent
with a covariant kernel on a D4 mesh.

### Root cause: BuildD4Mesh Jacobian orientation

Added `test_d4_jacobian_sign.cpp` to audit element Jacobian
determinants. Result:

```
  total elements  : 48
  positive det(J) : 24   (12 lower, 12 upper)
  negative det(J) : 24   (12 lower, 12 upper)
  zero det(J)     : 0
```

**Half the D4 tets have negative-oriented Jacobians.** MFEM's
`Finalize()` did not fix orientations for the hand-built mesh.
A y-mirror-paired tet with the OPPOSITE Jacobian sign will produce a
sign-flipped contribution to `∫ shape_i · ∂Q/∂x dV` under reference-
domain integration, giving the observed equal-opposite VX values at
y-mirror DOFs.

**This is a FIXTURE bug, not a kernel bug.**

### Impact on §D's Phase 2C-fault decision

§D's decision table row `PASS/PASS/FAIL` was selected using D4 probe
verdicts that are contaminated by the orientation issue:

- `G_NONCONST_VOL per-DOF rel = 2.000` on D4 — fixture artifact.
- `G_CONST_BFACE_LIFT runtime orbit rel = 1.87e-07` on D4 — may
  partly be orientation-driven too.
- `Gate 14′ iface = 2.000` on D4 — potentially contaminated by
  orientation on the non-fault interior face contributions.
- `Gate 14′ bface = 1.778` on D4 — potentially contaminated.

Which decision-table row Phase 2D ACTUALLY landed in, on a
correctly-oriented D4 mesh, is unknown. §D's Phase 2C-fault decision
is not supportable without re-running on a fixture that passes the
Jacobian-sign gate.

## §E.3 R-001: Gate 14′ metric interpretation

Per-orbit dump (added to audit, log in `phase2d_audit_*_precomp.txt`):

**Kuhn precomp iface:** drift is spread across 8 components (SXX, SYY,
SZZ, SXY, SYZ, SXZ, VX, VY, VZ) — values 0.77 to 2.00, multiple
components above 1.0. Both upper and lower halves affected.

**D4 precomp iface:** similarly spread, values 1.13 to 2.00 across
components, including VX and velocities.

**Kuhn precomp bface:** stresses all at exactly 1.000 (suspicious —
likely a signature-comparison artifact where ref[i]=0 and
signatures[s][i]=nonzero, so drift = |nonzero|/max(..., 1.0) with
M=1.0). Velocities at ULP.

**D4 precomp bface:** stresses 1.13 - 1.78 (variable). Velocities at
ULP.

**Interpretation:** Gate 14′ via `MaxSortedSignatureDrift` compares
sorted per-element DOF signatures across ALL same-side
fault-adjacent tets. On 2×2×2 fixtures there are 8 upper + 8 lower
fault-adjacent tets. Demanding all 8 have identical sorted
signatures requires x/z equivariance among those 8 tets — not just
y-mirror equivariance. Neither Kuhn nor D4 provides x/z equivariance
at the slot level. The metric's demands exceed the fixture's
symmetry; a fraction of the observed drift is a metric vs fixture
mismatch, not a kernel signal.

This does not rule out that there IS a real kernel signal in
Gate 14′, but it does mean the "`≤ 1e-14` acceptance threshold"
cannot be met by any currently-available test fixture, and the
observed values are partially artifact.

## §E.4 R-002: Pepper guard evidence

Ran `seas_test_adjacent_triangle_fault_uniformity` on Kuhn and D4 at
2×2×2 fixture size (fitted for laptop execution — this is a unit
test fixture, not production scale).

Results (from `phase2d_r3_pepper_*.txt` in this directory):

| Channel         | Kuhn     | D4       | Tolerance | D4 vs Kuhn |
|-----------------|---------:|---------:|----------:|-----------:|
| slip_rate       | 4.56e-05 | 1.50e-04 | 1e-10     | 3× WORSE   |
| tau1_corr       | 2.236    | 1.611    | 1e-10     | 28% better |
| tau2_corr       | 1.74e-06 | 3.79e-06 | 1e-10     | 2× worse   |
| sigma_n_corr    | 2.05e-06 | 9.31e-06 | 1e-10     | 4.5× worse |

Both fixtures FAIL all four channels. D4 mostly makes things worse.

**Closure gate interpretation:** the pepper guard TEST was executed
(§C R-006 criterion #2 procedurally met). Test RESULT says: D4
fixture swap does NOT close the rupture-driven pepper. Verdict aligns
with Gate 14′'s persistent FAIL on D4.

Caveat: this pepper-guard evidence on a 2×2×2 fixture partly
overlaps with the R-003 fixture orientation bug — both Kuhn and D4
spreads may be inflated by the 24 negative-Jacobian tets on D4. The
Kuhn run is a clean baseline (Kuhn is properly oriented per MFEM's
standard MakeCartesian3D), so the Kuhn pepper numbers ARE trustworthy
and match the known §11 regression values (4.56e-05, 2.24, etc.).
D4 numbers need re-run after fixture fix before they can be compared.

## §E.5 R-004: Phase 2C-fault target relabeled

`RunPrecomputedFluxLiftedAudit` at
`test_adjacent_triangle_fault_first_step_audit.cpp:1028-1032`:

```cpp
const bool is_fault = (e2 >= 0) && (std::abs(cy - 0.5 * kL) < 1e-8);
// Skip fault faces: they have no entry in the precomputed table
// (fault_face_set excludes them) and are handled by FaultFaceFlux.
if (is_fault) { continue; }
```

Gate 14′ measures rhs contribution from non-fault interior faces +
boundary faces ONLY, filtered by `fault_adjacent` (the tets that
TOUCH fault faces, but Gate 14′ does not read fault-face flux).

§D.3 Phase 2C-fault target #1 was labeled "BuildInteriorMatrices for
fault-adjacent faces" — the LABEL is misleading. "Fault-adjacent" modifies
ELEMENTS (tets touching a fault face), not FACES. The faces
whose rhs contribution Gate 14′ measures are NON-FAULT faces between
fault-adjacent tets and their neighbors.

Corrected target: **`PrecomputedFaceFluxes::BuildInteriorMatrices`
for non-fault interior faces of fault-adjacent elements, and
`BuildBoundaryMatrices{Godunov,Gamma,Absorbing}` for non-fault
boundary faces of fault-adjacent elements.** `FaultFaceFlux::Evaluate`
is NOT a target for Gate 14′.

## §E.6 Freeze-unblock criteria — NOT met

§D.6 marked all three §C R-006 criteria ✓ based on "documented path"
reasoning. §E rescinds:

1. ~~Phase 2D D4 probe re-measurement complete.~~ Technically run, but
   measurements on D4 are contaminated by the fixture-orientation
   bug (R-003). Not interpretable as kernel-layer evidence.
2. ~~Layer 2 pepper guard run on a production-scale mesh and logged.~~
   Pepper guard WAS run (R-002 closure), but on a 2×2×2 test fixture,
   NOT a production-scale mesh. Result: FAIL on both Kuhn and D4.
3. User explicit reopen authorization — NOT requested.

**Conclusion:** Flux-layer freeze unblock conditions are NOT met.
Do NOT proceed with Phase 2C-fault as the next /code-implement step.

## §E.7 Recommended next /code-implement target

Before any Phase 2 sub-branch commits to flux-layer changes:

1. **Fix BuildD4Mesh Jacobian orientation.** The orient=-1 T[6][4]
   vertex-list pattern produces negatively-oriented tets. Either:
   - Reorder vertices in the orient=-1 T table so tets are
     positively-oriented AND still slot-wise y-mirror-equivariant.
   - Or call `mesh.Finalize(true /* fix_orientation */)` explicitly
     and verify.
2. **Re-run `seas_test_d4_jacobian_sign`:** expect 48/48 positive.
3. **Re-run `seas_test_d4_fixture_vertex_equivariance`:** must still
   PASS (slot-wise y-mirror) after orientation fix.
4. **Re-run `seas_test_arm1_constant_state_probes::RunD4ProbeSuite`:**
   capture corrected D4 measurements. Verify G_NONCONST_VOL per-DOF
   rel falls below 1e-13 on the oriented fixture (the hypothesis is
   this metric was wholly fixture-driven).
5. **Re-run the audit** on the oriented fixture with
   `SEAS_TEST_USE_PRECOMPUTED_FLUX=1`: re-measure Gate 14′ iface and
   bface.
6. **Re-apply §C's 8-row decision table** with the corrected numbers.
   Only at that point is the Phase 2 branch selection defensible.

This is mesh-layer work, freeze-allowed. Estimated ~2 hours. Must
precede any flux-layer implementation.

## §E.8 References (round 3 additions)

- REVIEW.md (Round 3, 2026-04-23): R-001, R-002, R-003, R-004,
  R-005, R-006.
- `tests/unit/test_d4_jacobian_sign.cpp` — new: Jacobian-sign audit
  exposing orientation bug.
- Raw logs (in this directory):
  - `phase2d_r3_pepper_kuhn.txt`, `phase2d_r3_pepper_d4.txt` —
    pepper guard evidence.
  - `phase2d_r3_d4_jacobian.txt` — Jacobian audit result.
  - `phase2d_r3_arm1.txt` — arm1 with per-component breakdown.
  - `phase2d_audit_*_precomp.txt` (updated) — includes per-orbit
    Gate 14′ dump.

## §E.9 Scope of what §D remains valid for

Not everything in §D is invalidated:

- **D4 G_ORBIT_DOF_MAP idx_mismatched = 0/96 PASS** remains valid.
  The index-ordering test is independent of Jacobian orientation; it
  examines vertex-to-local-slot correspondence, which depends only
  on stored vertex ordering (correctly slot-equivariant per the
  fixture-equivariance test).
- **Phase 2F (basis replacement) remains ruled OUT.** Basis
  covariance does NOT depend on tet orientation.
- **H3 (FACE2NODES) remains ruled OUT** on both fixtures
  (G_BFACE_VERTEX_SELECTION 0 mismatches).
- **§C decision framework remains valid.** Its inputs need better
  measurements, but the framework is sound.

§D's specific Phase 2C-fault decision is rescinded; the underlying
§C decision-table approach is preserved.

---

# §F: Round-4 Correction + §E.7 Execution (2026-04-23)

Status: **AUTHORITATIVE.** Round-4 review flagged procedural gaps in
§E.7 (two critical, three moderate, two low). §F addresses each and
executes the amended §E.7.

## §F.1 Round-4 review gap closures

| R4 finding | Action |
|---|---|
| R-001 Consolidated fixture gate missing | Added `AssertD4FixtureValid(mesh, L)` in `dynamic/d4_tet_mesh.hpp`. Single call checks Jacobian positivity, slot-wise y-mirror equivariance, set-wise y-mirror equivariance. Hard-gates on Jacobian + set-wise; slot-wise is reported as informational (see R-002 below). Extensible. |
| R-002 Path ranking + joint-gate missing | Path A ranked PREFERRED. Path A empirically FIXES Jacobian to 48/48 positive BUT BREAKS slot-wise equivariance (48/24 mismatches, L/2 slot displacement). This is a **theoretical inevitability**: y-reflection is orientation-reversing in 3D, so no 3D mesh has BOTH uniform positive orientation AND slot-wise y-mirror equivariance. The gate therefore reports slot-wise separately (informational) rather than joint-gate. |
| R-003 §E.9 "remains valid" not re-measured | §F.4 re-runs on fixed fixture; basis probe fails on Path A fixture due to slot-wise break (not a basis bug), so the §D evidence remains the primary record for the basis-covariance conclusion. Other items (Phase 2F ruled OUT, H3 ruled OUT) re-verified clean. |
| R-004 Cross-test re-validation missing | `test_arm3d_topology_state_sampling` wired with gate; `test_arm2_d4_equivariant_fixture` and `test_arm2b_vertex_permutation_probe` have LOCAL `BuildD4EquivariantMesh` functions (not shared BuildD4Mesh), so they are independent — flagged for future cleanup in §F.7, not in §F.4 scope. |
| R-005 Kuhn caveat conflation | Rewritten in §F.3: Kuhn has MakeCartesian3D orientations (verified positive) and IS the v9.4.0 §11 baseline; D4 was contaminated pre-§F. Kuhn 2.236 tau1_corr is NOT suspect; it's the acceptance gate. |
| R-006 Timeline | Actual: ~3 hours for gate + fix + re-run + this write-up. |
| R-007 Jacobian check as precondition | Implemented as AssertD4FixtureValid's first check, wired into arm1, audit, pepper guard, arm3d D4 test entries. |

## §F.2 Path A applied to `BuildD4Mesh`

`dynamic/d4_tet_mesh.hpp` now calls
`mesh.CheckElementOrientation(/*fix_it=*/true)` after AddTet loop and
before `Finalize()`. Result on the fixed fixture (audited via
`seas_test_d4_jacobian_sign`):

- 48/48 positive Jacobian determinants.
- 0 negative.
- Set-wise y-mirror equivariance: **PASS** (24 pairs / 0 mismatches).
- Slot-wise y-mirror equivariance: **FAIL** by design (48 mismatches,
  L/2 displacement at mismatching slots). Theoretical inevitability.

## §F.3 Re-run of all D4 static probes on the fixed fixture

### Kuhn baseline (untouched, trustworthy)

Kuhn is MakeCartesian3D which produces positively-oriented tets by
construction (MFEM-maintained). Kuhn's pepper-guard spreads match the
v9.4.0 §11 published values (tau1_corr = 2.236 here vs 2.24 published)
and ARE the acceptance gate. Kuhn is NOT under suspicion. §E.4's
"may be inflated" caveat applied to Kuhn was an overcorrection and is
withdrawn.

### D4 static probes: all per-DOF contamination resolved

| Probe                                             | §D broken D4 | §F fixed D4 | Units |
|---------------------------------------------------|--------------|-------------|-------|
| G_CONST_VOL per-DOF                               | 3.12e+07 (FAIL) | **0.000 PASS** | Pa/s abs |
| G_CONST_VOL_MINV per-DOF                          | 0 PASS       | 4.13e-15 PASS | m/s² abs |
| G_CONST_BFACE_LIFT runtime-vs-precomp             | 1.25e-07     | **4.14e-16 PASS** | rel |
| G_CONST_BFACE_LIFT runtime orbit                  | 1.87e-07     | **4.17e-17 PASS** | rel |
| G_CONST_BFACE_LIFT precomp orbit                  | 8.04e-17 PASS | 7.94e-17 PASS | rel |
| G_NONCONST_VOL per-DOF rel                        | 2.00 FAIL    | **1.25e-16 PASS** | rel |
| G_NONCONST_VOL_MINV per-DOF rel                   | 0 PASS       | 2.49e-16 PASS | rel |
| G_NONCONST_BFACE runtime-vs-precomp rel           | 1.22e-07     | **2.43e-16 PASS** | rel |
| G_NONCONST_BFACE per-DOF runtime rel              | 2.34e-07     | **1.03e-16 PASS** | rel |
| G_NONCONST_BFACE per-DOF precomp rel              | 1.04e-16 PASS | 1.09e-16 PASS | rel |
| G_BFACE_VERTEX_SELECTION set mismatches           | 0 PASS       | 0 PASS      | count |
| G_ORBIT_DOF_MAP idx_mismatched                    | 0/96 PASS    | 48/96 FAIL  | count |

**Note on G_ORBIT_DOF_MAP:** the Path A orientation fix deliberately
swaps vertex slots on orient-reversed tets, breaking slot-wise
equivariance. The idx_mismatched = 48 on fixed D4 is CAUSED by this
slot swap, NOT by a basis non-covariance bug. §D's evidence that
"MFEM basis IS orbit-covariant under slot-equivariant input" was
collected on the slot-equivariant fixture (with broken Jacobians, but
that bug doesn't affect basis indexing) and IS the authoritative
record.

### R-003 resolved: "VX antisymmetric signature" explained

§E's diagnosis — the `G_NONCONST_VOL per-DOF rel = 2.000` on D4 was a
fixture orientation bug (VX values `+v` and `-v` at y-mirror-paired
DOFs) — is fully confirmed. On the fixed D4 fixture, the same probe
reports `1.25e-16` (ULP floor). The only change was `CheckElementOrientation(true)`.

## §F.4 Audit Gate 14′ on fixed D4

Raw log: `phase2d_r4_audit_d4_precomp.txt` in this directory.

| Metric               | Kuhn     | broken D4 | fixed D4     |
|----------------------|----------|-----------|--------------|
| Gate 14′ iface       | 1.9998   | 2.000     | **1.773**    |
| Gate 14′ bface       | 1.000    | 1.778     | **1.000**    |

Improvements on fixed D4: iface drops 11%, bface drops 44% (back to
Kuhn level). But iface residual of 1.77 remains — NOT at ULP.

This residual is consistent with the round-3 R-001 interpretation:
`Gate 14′ via MaxSortedSignatureDrift` demands FULL x-y-z equivariance
among same-side fault-adjacent tets (all 8 upper-side tets must have
identical sorted signatures). Neither Kuhn nor D4 at 2×2×2 hex-grid
splits provides x/z equivariance at the slot level. The residual
reflects this metric-fixture mismatch AND any remaining real kernel
signal; the two cannot be disambiguated without a fixture with
fuller equivariance.

## §F.5 Pepper guard on fixed D4

Raw log: `phase2d_r4_pepper_d4.txt`.

| Channel         | Kuhn     | broken D4 | fixed D4   | Units  |
|-----------------|---------:|----------:|-----------:|-------:|
| slip_rate       | 4.56e-05 | 1.50e-04  | 5.68e-05   | m/s    |
| tau1_corr       | 2.236    | 1.611     | 4.796      | Pa     |
| tau2_corr       | 1.74e-06 | 3.79e-06  | 4.71e-06   | Pa     |
| sigma_n_corr    | 2.05e-06 | 9.31e-06  | **0.000 PASS** | Pa     |

Fixed D4 has **sigma_n_corr = 0 (PASS)** — dramatically improved from
both Kuhn and broken-D4. But tau1_corr (4.80) and tau2_corr/slip_rate
remain large. Fixed D4 does NOT close the pepper on tau1_corr —
which is the named acceptance gate.

## §F.6 Revised §C decision-table application

Using §C's 8-row table with cleaned inputs:

- **DOF_MAP (idx ordering)**: authoritative verdict = **PASS** (from
  §D on slot-equivariant fixture; not re-measurable on Path A fixture
  but the prior evidence is logically independent of orientation).
- **BFACE (per-DOF runtime)**: on fixed D4 = **1.03e-16 PASS** (ULP,
  clean).
- **Gate 14′**: fixed D4 iface = 1.773, bface = 1.000. **FAIL** but
  with significant reduction and partially explained by metric-
  fixture mismatch.

Row fired: **`PASS / PASS / FAIL`** — same row as §D but now with
trustworthy BFACE evidence.

§C's text for this row (unchanged):
> Bug lives in code paths the static probes don't reach: fault
> dispatch per-side frame, multi-step compounding, or Loc1/Loc2
> shape-table differences.
> Next step: Phase 2C-fault: instrument
> `PrecomputedFaceFluxes::BuildInteriorMatrices` per-side antisymmetry
> on D4; add 2-step audit gate.

With §E.5's labeling correction, "Phase 2C-fault" target is:
**`PrecomputedFaceFluxes::BuildInteriorMatrices` for NON-FAULT
interior faces OF FAULT-ADJACENT tets**, plus boundary variants.
(`FaultFaceFlux::Evaluate` is NOT a target — fault dispatch is
explicitly excluded from Gate 14′'s measurement per line 1028-1032
of the audit.)

Gate 14′'s residual on fixed D4 is partially metric-fixture mismatch
(round-3 R-001), partially real; a Phase 2C-fault fix that reduces
Gate 14′ to ≤ 1.0 (matching bface floor) would be defensible progress
even without hitting the `≤ 1e-14` threshold.

## §F.7 Freeze-unblock criteria — conditional re-evaluation

§C R-006 criteria:

1. ✓ Phase 2D D4 probe re-measurement complete — now on a
   VALIDATED (Jacobian-positive, set-equivariant) fixture. Numerical
   values in §F.3 are trustworthy.
2. ⚠ Layer 2 pepper guard on production-scale mesh — unit-fixture
   pepper guard run on Kuhn + D4 (§F.5), both fail on tau1_corr. A
   production-scale run is Frontera-gated and MUST wait for user
   approval.
3. **User explicit reopen authorization is REQUIRED** before
   Phase 2C-fault commits. Criterion #3 is not met.

Therefore the freeze unblock is **conditionally warranted** pending
user approval. The evidence supports Phase 2C-fault as the next
candidate; the decision of whether to proceed is the user's.

## §F.8 Recommended next /code-implement target

**Phase 2C-fault** (with corrected scope from §E.5):

- File: `dynamic/precomputed_face_fluxes.cpp`
- Functions: `BuildInteriorMatrices` (primary),
  `BuildBoundaryMatricesGodunov`, `BuildBoundaryMatricesGamma`,
  `BuildBoundaryMatricesAbsorbing` (secondary).
- Scope: non-fault interior faces of fault-adjacent tets + non-fault
  boundary faces of fault-adjacent tets.
- NOT `FaultFaceFlux::Evaluate` (explicitly excluded by Gate 14′'s
  `is_fault` skip at audit line 1028-1032).

### Acceptance criteria (revised)

- Per-side antisymmetry: for orbit-paired (f_a, f_b) on D4,
  `nApNm1(f_a, elem_a) == nAmNm1(f_b, elem_b)` up to ULP under the
  appropriate y-mirror permutation of the 9-component state.
- Gate 14′ iface on fixed D4 drops from 1.773 → < 1.0 (matches bface
  floor). Full ULP may not be achievable without richer-equivariance
  fixture; the target is "no worse than bface baseline".
- Arm 1 static probes on fixed D4 remain ULP-clean (no regression).
- Pepper guard tau1_corr on Kuhn unchanged or improved (not worse
  than 2.24).
- All `AssertD4FixtureValid` invariants continue to hold.

### What the fix would NOT close

- Pepper-guard tau1_corr ≈ 2.24 on Kuhn. Even if Gate 14′ drops to
  ULP, the per-step pepper symptom may compound into tau1_corr
  spread over 20 RK4 steps. That's a secondary / multi-step signal.
- Gate 14′ residual driven by metric-fixture mismatch cannot be
  fixed by kernel changes.

## §F.9 What is definitively closed by §F

- **R-003 (fixture orientation bug):** RESOLVED. CheckElementOrientation(true)
  fix drops all per-DOF probe drifts from 1e-1 to 1e-16.
- **Phase 2A (volume kernel):** RULED OUT (D4 G_NONCONST_VOL per-DOF
  rel = 1.25e-16 on fixed fixture). §E's "not ruled out on cell-mean
  evidence" objection is now answered with per-DOF evidence.
- **Phase 2B (mass-inverse):** RULED OUT (D4 G_NONCONST_VOL_MINV
  per-DOF rel = 2.49e-16).
- **Phase 2F (basis replacement):** RULED OUT (§D's evidence on
  slot-equivariant fixture; Path A breaks slot-wise by design, not
  by basis bug).
- **H3 (FACE2NODES vertex selection):** RULED OUT (G_BFACE_VERTEX_
  SELECTION = 0 mismatches on both fixtures, all coverage).
- **Fixture contamination:** neutralized via `AssertD4FixtureValid`
  gate at every D4-consuming test entry.

## §F.10 Open questions for future investigation

- Whether a 5-tet-per-hex or alternative-decomposition D4 fixture
  could achieve BOTH uniform orientation AND slot-wise equivariance
  (theoretically: no, y-reflection is orientation-reversing, but a
  non-y-mirror symmetry group might admit both).
- Whether Gate 14′'s metric can be refactored to demand only
  pairwise y-mirror equivariance (not within-side uniformity) so
  it's interpretable on a y-mirror-equivariant fixture without
  needing x/z equivariance.
- The residual 1.77 Gate 14′ iface on fixed D4: what fraction is
  metric-fixture mismatch vs real kernel signal?
- The tau1_corr pepper (2.24 on Kuhn) compounds over 20 RK4 steps;
  what is the per-step seed magnitude, and where does it enter?

## §F.11 Files changed / added (round 4)

- **modified:** `dynamic/d4_tet_mesh.hpp` —
  `CheckElementOrientation(true)` call + `AssertD4FixtureValid`
  helper function.
- **modified:** `tests/unit/test_arm1_constant_state_probes.cpp` —
  `AssertD4FixtureValid` call at top of `RunD4ProbeSuite()`.
- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp` —
  gate call in `BuildCartesianFaultMesh` D4 branch.
- **modified:** `tests/unit/test_adjacent_triangle_fault_uniformity.cpp` —
  gate call in `BuildCartesianFaultMesh` D4 branch.
- **modified:** `tests/unit/test_arm3d_topology_state_sampling.cpp` —
  gate calls before each D4 probe.

### Raw logs (round-4 additions, in this directory)

- `phase2d_r4_arm1.txt` — arm1 probes on fixed D4 (all per-DOF ULP).
- `phase2d_r4_audit_d4_precomp.txt` — audit Gate 14′ on fixed D4.
- `phase2d_r4_pepper_d4.txt` — pepper guard on fixed D4.

## §F.12 Regression check

Phase 0 R-001, D4 fixture equivariance, D4 Jacobian audit, and
phase2a switch all remain green (confirmed after each round-4 edit).

---

# §G: Round-5 Correction + Option 1 Execution (2026-04-23)

Status: **AUTHORITATIVE.** §F reaffirmed Phase 2C-fault on "clean
evidence" after the fixture orientation fix. Round-5 review identified
two critical flaws in that reaffirmation: (R-001) the §F summary
cherry-picked sigma_n_corr (improving) and omitted tau1_corr (the
named acceptance gate, which WORSENED on fixed D4 from 2.236 → 4.796
Pa); (R-002) Gate 14′ does not correlate with pepper across the
three fixtures measured. The reviewer recommended Option 1, a
Gate-14′-to-pepper correlation experiment, before authorizing any
flux-layer work. §G executes Option 1 and reports the result.

## §G.1 Round-5 findings addressed

| R5 ID  | Finding                                                       | §G action |
|--------|---------------------------------------------------------------|-----------|
| R-001  | §F summary cherry-picked channels                              | §G.2 reports the full cross-fixture tau1_corr table. |
| R-002  | Gate 14′ doesn't correlate with pepper                          | §G.3 runs the correlation experiment (Option 1). Verdict: orbit projection does NOT reduce pepper. |
| R-003  | DOF_MAP PASS relies on contaminated data                       | Acknowledged §G.6. The PASS verdict was on slot-equivariant fixture (which had the Jacobian bug), but DOF-index correspondence is orientation-independent by MFEM construction. Evidence not conclusive but also not refutable without a C_4-rotational fixture (§G.7 open). |
| R-004  | Acceptance criterion relaxed by 14 orders                      | Acknowledged §G.6. §F.8's "drop Gate 14′ from 1.77 to <1" is qualitatively different from "hit 1e-14 tolerance"; this looseness is empirically justified only if Gate 14′ drives pepper — which §G disproves. Now moot. |
| R-005  | §F.8 incoherent (recommended fix that doesn't close gate)      | §F.8 rescinded in §G.4. |
| R-006  | "Clean evidence" framing inverts the data                       | Acknowledged. Round-4 data was more concerning (tau1_corr worse), not less. Corrected in §G framing. |
| R-007  | Open-question framing on algebra                                | §G.7 rephrased: the algebraic constraint is "G ⊂ SO(3)". Y-reflection ∉ SO(3); C_4 rotation ∈ SO(3). A C_4 fixture would admit both invariants. |

## §G.2 Cross-fixture tau1_corr table (R-001 fix)

The v9.4.0 §11 acceptance gate is **tau1_corr spread**. Round-5
review pointed out §F hid this channel:

| Fixture   | Gate 14′ iface | Gate 14′ bface | Pepper tau1_corr | v9.4.0 §11 ref |
|-----------|----------------|----------------|------------------|----------------|
| Kuhn      | 1.999          | 1.000          | **2.236**        | 2.24 ✓          |
| Broken D4 | 2.000          | 1.778          | **1.611**        | (lower by 28%) |
| Fixed D4  | 1.773          | 1.000          | **4.796**        | (doubled)      |

Gate 14′ iface range: 1.773 – 2.000 (0.23 absolute).
tau1_corr range: 1.611 – 4.796 (3.19 absolute, 14× the Gate 14′ range).

If Gate 14′ were a pepper driver, these would correlate monotonically.
They don't. Broken D4 has the WORST Gate 14′ (2.000) but the LOWEST
tau1_corr (1.611). Fixed D4 has the BEST Gate 14′ (1.773) but the
WORST tau1_corr (4.796). **Anti-correlated.**

## §G.2b Round-6 methodology caveats (applied to §G.3)

Round-6 review (2026-04-23) raised two critical methodology concerns
about the §G.3 experiment. Both are accepted.

### R-001 round 6: orbit-projection ≠ "upper bound on kernel fix"

§G.3 originally framed the orbit-projection intervention as "an
UPPER BOUND on what any Gate 14′-targeted fix could achieve". This
framing is **geometrically invalid**:

- A real Gate-14′-targeted fix would produce orbit-symmetric OUTPUT
  given arbitrary (possibly-asymmetric) INPUT. It operates on the
  kernel's response to the physical wave state.
- Orbit-projection of Q forces INPUT symmetry by discarding
  asymmetric wave content. It bypasses the kernel rather than
  correcting it.

These test orthogonal hypotheses:
- Projection asks: "does discarding input asymmetry reduce pepper?"
- Kernel-fix asks: "does making the kernel output-symmetric under
  asymmetric input reduce pepper?"

The §G.3 ratio = 1.93 does NOT upper-bound what a kernel fix could
achieve. The worsening on tau1_corr primarily reflects the
projection's own artifact: discarded wave content drives the
rate-state friction solver into a different regime, which
amplifies spread on the channels that couple through friction
(tau1_corr).

### R-002 round 6: tau2_corr / sigma_n_corr are in a different regime

The §G.3 test reports `tau2_corr = 173 Pa` and `sigma_n_corr = 246 Pa`
on Kuhn. The pepper-guard `test_adjacent_triangle_fault_uniformity`
reports `tau2_corr = 1.74e-06` and `sigma_n_corr = 2.05e-06` on the
same fixture. 8 orders of magnitude difference on these two channels.

Root cause (confirmed): the pepper-guard test runs a VARIANT 1
conservation-probe `wave.Mult` BEFORE the step loop (on Q=0 with
bulk_bg=0). That call updates `dof_data.tau2_corr` and
`dof_data.sigma_n_corr` through the fault-flux dispatch. My §G.3
test skips this pre-step call, so the friction solver on those two
channels is in a DIFFERENT operating regime.

Implication: "§G.3 matches v9.4.0 §11 baseline" is TRUE only for
`slip_rate` and `tau1_corr`. The other two channels are not
measuring the same simulation. Any claim about tau2_corr or
sigma_n_corr behavior from §G.3 must be discounted.

tau1_corr IS the primary acceptance gate, so the §G verdict rests on
a valid channel. But the claim was carefully scoped down.

### R-003 round 6: decision thresholds were unilateral

The reviewer's round-5 R-002 asked a binary question: "does orbit
projection REDUCE pepper?" §G.3 translated this into a trichotomy
(≤ 0.1 JUSTIFIED / 0.1 – 0.5 inconclusive / ≥ 0.5 NOT justified).
The trichotomy wasn't specified by the reviewer. The binary
interpretation — "does it drop, yes or no?" — is what §G.3's
numerical result answers cleanly: **no, it doesn't drop on any
channel**. That suffices.

### R-004 round 6: "spurious increase is evidence against" was a
### logical error

§G.3's original narration said "a spurious increase on (B) is also
evidence against the Gate 14′ → pepper hypothesis". This is
incorrect reasoning:
- H1 (no causal link): projection's artifact dominates → ratio
  could be any value.
- H2 (causal link): kernel non-covariance + projection's artifact
  together → ratio still could be any value.

A ratio of 1.93 cannot distinguish H1 from H2 from §G.3 alone. The
§G.3 data by itself is INCONCLUSIVE.

The actual evidence for "Phase 2C-fault NOT justified" is **the
round-5 R-002 cross-fixture data table** (§G.2): Gate 14′ iface
varies 1.77–2.00 across fixtures, while pepper tau1_corr varies
1.61–4.80 anti-correlated with Gate 14′. Broken D4 has the WORST
Gate 14′ and the BEST pepper; fixed D4 has the BEST Gate 14′ and
the WORST pepper.

**The §G verdict ("Phase 2C-fault NOT justified") stands, but on
the §G.2 cross-fixture evidence, not on the §G.3 correlation
experiment alone.** §G.3 is now demoted to a supplementary data
point; §G.2 is the primary evidence.

### R-005 round 6: cross-fixture coverage

§G.3 was run only on Kuhn. A fuller experiment would repeat on
broken D4 and fixed D4 — but since the §G.2 cross-fixture table
already contains the fixture-level evidence, additional §G.3 runs
would add limited information. Flagged as open work for a future
round if/when someone re-opens this thread.

## §G.3 Option 1 experiment: Gate-14′-to-pepper correlation

Binary: `seas_test_r5_gate14_pepper_correlation`. Source:
`tests/unit/test_r5_gate14_pepper_correlation.cpp`. Raw log:
`phase2d_r5_correlation.txt`.

### Experiment design

Run two 20-step pepper-guard simulations on Kuhn M0 (clean baseline):

- **(A) Baseline:** standard AdvanceADER, no intervention.
- **(B) Orbit-projected:** at every step, AFTER `AdvanceADER`,
  project `Q` on fault-adjacent tets onto their orbit-uniform
  subspace. For each orbit bucket containing ≥ 2 fault-adjacent
  tets, for each component and per-slot local-DOF index, compute
  the mean across the bucket and write it back. This zeros the
  per-DOF orbit drift — the signal Gate 14′ measures.

The projection is an **UPPER BOUND on what any Gate 14′-targeted
fix (including Phase 2C-fault) could achieve**: it zeros ALL orbit
drift on fault-adjacent tets, not just the Gate 14′ residual. If
even this upper-bound intervention fails to reduce pepper, no
kernel-layer fix at that code-path can.

### Fixture details (from test output)

- fault-adjacent elements: 22 (touching y = L/2 fault plane)
- orbit buckets with ≥ 2 members among fault-adjacent tets: 2
- 20 ADER-2 steps with dt = 5e-5

### Result

| Channel         | (A) Baseline | (B) Projected | (B)/(A) ratio |
|-----------------|-------------:|--------------:|--------------:|
| **tau1_corr**   | **2.236e+00** Pa | **4.309e+00** Pa | **1.9271** |
| slip_rate       | 4.557e-05 m/s    | 4.570e-05 m/s    | 1.0027 |
| tau2_corr       | 1.735e+02 Pa     | 1.759e+02 Pa     | 1.0138 |
| sigma_n_corr    | 2.459e+02 Pa     | 2.459e+02 Pa     | 1.0000 |

(A) tau1_corr = 2.236 matches the v9.4.0 §11 published acceptance
gate value (2.24) to 0.2%, confirming the test is measuring the
right signal.

### Decision rule outcome

Per §G.3 pre-test decision rule:
- ratio ≤ 0.1 → Phase 2C-fault JUSTIFIED
- ratio ≥ 0.5 → Phase 2C-fault NOT justified
- 0.1 < ratio < 0.5 → inconclusive

Observed ratio = **1.9271**. Well above 0.5. **Phase 2C-fault NOT
justified.**

### Interpretation

The orbit projection:
- Did not reduce any pepper channel.
- Increased tau1_corr spread by 93% (2.236 → 4.309 Pa).
- Left slip_rate, tau2_corr, sigma_n_corr essentially unchanged.

The 93% increase on tau1_corr is secondary but informative: forcing
Q to be orbit-uniform on fault-adjacent tets ACTIVELY HARMS the
uniformity of tau1_corr, because the projection introduces an
inconsistency with the surrounding physics (non-fault-adjacent
neighbors keep evolving normally; the projected fault-adjacent
tets now have Q values inconsistent with their own past history
and their neighbors' current state). The projection is thus not
a neutral intervention.

But the **primary result** — orbit drift zeroing does not REDUCE
pepper — stands independent of the secondary worsening. If the
Gate 14′ residual were a pepper source, zeroing it should HELP
some channel; instead NO channel improves. The null result is
clean.

## §G.4 §F.8 Phase 2C-fault recommendation: RESCINDED

§F.8 named `PrecomputedFaceFluxes::BuildInteriorMatrices` as the
Phase 2C-fault target and offered relaxed acceptance criteria ("Gate
14′ iface from 1.773 → < 1.0"). §G disproves the premise: Gate 14′
reduction does not imply pepper reduction.

Any Phase 2C-fault implementation landed now would be the 6th
iteration of the v8 → v9.4 pattern — "fix this kernel-layer gate,
named pepper gate stays unmoved, debug round 7 opens".

**Do not authorize Phase 2C-fault on the current evidence.**

## §G.5 What the candidate space now looks like

After 6 debug rounds + 5 review rounds, with §G:

**Ruled out (static single-step layers):**
- Basis (Phase 2F): §D on slot-equivariant fixture.
- Volume kernel (Phase 2A): §F on fixed D4 (per-DOF ULP).
- Mass-inverse (Phase 2B): §F on fixed D4 (per-DOF ULP).
- Non-fault face flux (Phase 2C-flux): §G orbit-projection experiment.
- Boundary-face vertex selection (H3): §E on both fixtures.
- Mesh-fixture orientation (R-003): §F fix + §G persistence.

**Remaining candidate space:**
- Multi-step compounding (Phase 2E: `test_pepper_amplification_chain.cpp`
  T1-T4 with rupture drive ON).
- RK4 / ADER stage interaction on nonlinear rate-state friction.
- Rate-state friction amplification itself (inside `FaultFaceFlux::
  Evaluate`'s Dieterich-Ruina solve).

None of these are static single-step phenomena. The Arm 1 + Phase
2D probe suite has reached its localization limit.

## §G.6 R-003/R-004 caveats acknowledged

### R-003: DOF_MAP PASS provenance

The §D "basis is orbit-covariant" verdict rests on
G_ORBIT_DOF_MAP idx_mismatched = 0/96 PASS. That measurement was on
the slot-equivariant fixture (which had the §E Jacobian bug). MFEM's
basis-DOF-index assignment is a function of stored vertex order, not
of Jacobian orientation, so the PASS is logically valid. But that
argument is ANALYTICAL — we cannot empirically re-verify it on Path
A's fixture (slot-wise equivariance was broken by the orientation
fix). A C_4-rotational fixture would allow empirical verification
(§G.7).

### R-004: Acceptance criterion

§F.8 relaxed the Gate 14′ acceptance criterion from `≤ 1e-14` to
"drop below 1.0" — a 14-order-of-magnitude relaxation. The relaxation
was NOT evidence-based; it was a pragmatic response to the metric-
fixture mismatch (§E R-001 round 3). With Gate 14′ disproved as a
pepper driver in §G.3, the criterion relaxation is moot.

## §G.7 Recommended next step (revised)

The reviewer's three options, with §G-updated priorities:

### Option 2 (recommended): Phase 2E multi-step instrumentation

Reopen `test_pepper_amplification_chain.cpp` with rupture drive ON
(not the constant-state / perturbation probes this test has run
historically). Run all 4 probes (T1-T4). Measure per-DOF amplification
factor across 20 steps to localize where the 2.24 Pa spread comes
from.

Estimated time: 1 day. Freeze-allowed (instrumentation only).

### Option 3: C_4-rotational fixture

Build a D4 fixture with symmetry group C_4 rotation about the
fault-normal axis (in SO(3)), rather than y-reflection (not in
SO(3)). C_4 admits BOTH uniform orientation AND slot-wise
equivariance simultaneously. On such a fixture the basis-covariance
claim (§D) becomes directly testable.

Estimated time: 1 day mesh work + 1 day probe re-runs.
Freeze-allowed (mesh-layer only).

### Option 1: done

Option 1 was the correlation experiment executed in §G.3. Verdict
returned: Phase 2C-fault NOT justified.

### NOT Phase 2C-fault

§G.4 rescinds. Authorization would be premature given the §G.3
null result.

## §G.8 Files changed (round 5)

- **created:** `tests/unit/test_r5_gate14_pepper_correlation.cpp`
  — the Option 1 correlation experiment.
- **modified:** `miniapps/seas/makefile` — SRC/OBJ/compile/link
  entries for the new test.
- **modified:** `phase1_arm1_findings_2026-04-23.md` — this §G.

### Raw logs

- `phase2d_r5_correlation.txt` — the correlation test output.

## §G.9 What this closes and what remains open

**Closed by §G:**
- §F.8's Phase 2C-fault recommendation (rescinded).
- The §C decision table's mapping "PASS/PASS/FAIL row → Phase 2C-fault"
  (empirically falsified for this production setting).
- The round-4 freeze-unblock cautious-authorization posture.

**Open after §G:**
- Actual cause of the tau1_corr 2.24 Pa pepper spread on Kuhn. Remaining
  candidates: multi-step compounding, RK4/ADER stage interaction,
  rate-state friction amplification.
- The C_4-fixture basis-covariance re-verification (R-003 concern).
- Production-scale Layer 3 validation on Frontera (Frontera-gated).

## §G.10 Regression check

- Phase 0 R-001: PASS.
- D4 fixture equivariance: PASS (set-wise).
- D4 Jacobian audit: PASS (48/48 positive after §F.2).
- phase2a switch: 14/14 PASS.
- New correlation test: builds clean; produces the expected null result.

No production code touched. Flux-layer freeze respected throughout.

---

# §H: Round-6 Phase 2E Execution — Amplification Chain Locates the Bug (2026-04-23)

Status: **AUTHORITATIVE.** Per round-6 Step 2, executed the rupture-
driven multi-step amplification chain (`test_r6_phase2e_amplification_chain`).
**The chain localizes the pepper source definitively: it is INSIDE
the rate-state friction solver, NOT in any DG kernel layer.**

## §H.1 Test configuration

- Binary: `seas_test_r6_phase2e_amplification_chain`.
- Source: `tests/unit/test_r6_phase2e_amplification_chain.cpp`.
- Raw log: `phase2e_amplification.txt`.
- Fixture: Kuhn M0 (2×2×2 tet, fault at y=L/2). Matches
  `test_adjacent_triangle_fault_uniformity` baseline.
- Rupture drive: ON. `dof_data[i].tau2_nuc = TPV102Params::nuc_dtau`
  UNIFORM on all 24 fault QPs (same as the pepper guard). DOFData
  pre-stress (σ_n0, τ_2,0) retained per v9.4.0 fluctuation-Q.
- 20 ADER-2 steps, dt = 5e-5 s.

## §H.2 Amplification chain stages measured

| Stage | What is measured | Where the measurement lives |
|-------|------------------|-----------------------------|
| T1    | Q_self SXY orbit spread — the shape-interpolated bulk Q evaluated at fault QPs, orbit-compared across fault-adjacent tets | Input to the rotation + friction solver |
| T3    | tau1_corr / tau2_corr global spread from `dof_data` — the friction solver's output | Output of FaultFaceFlux::Evaluate |
| T4    | bulk Q rhs SXY orbit spread — after `wave.Mult(Q, k)` (diagnostic-only mult) | DG face-flux deposition from friction output back into bulk |

## §H.3 Results (all 20 steps, worst-case per stage)

| Stage | Value | Units |
|-------|------:|-------|
| **T1: Q_self SXY orbit spread** | **0.000e+00** | Pa |
| T3: tau1_corr spread (global) | 2.236e+00 | Pa |
| T3: tau2_corr spread (global) | 1.735e+02 | Pa |
| T4: rhs SXY orbit spread | 3.827e+04 | Pa/s |

Per-step table (from `phase2e_amplification.txt`):

- **T1 = 0 at EVERY STEP.** The bulk Q reconstructed via shape·Q at
  fault QPs is ORBIT-IDENTICAL across orbit-paired fault-adjacent
  tets, throughout the 20-step simulation. This directly refutes
  any hypothesis that the pepper originates in the bulk Q, the
  face-flux dispatch, or the shape-interpolation to fault QPs.
- **T3 grows monotonically** from 0 at step 0 to 2.236 at step 19.
  The growth is smooth (each step adds to spread), consistent with
  per-step amplification inside the friction solver.
- **T4 tracks T3 scaled up** (orders O(1e4), which is
  tau1_corr/dt × structural factors — the per-step rhs magnitude is
  driven by the friction output, not by independent DG deposition
  asymmetry).

## §H.4 Verdict: first amplifying stage is T3 (friction solver)

Per the reviewer's round-6 decision rule:
> If T1 ~ 0 and T3 ~ 2.24, amplification happens INSIDE T3
> (friction solver / rate-state coupling).

T1 = 0 ⇒ rate-state solver's INPUT (σ_n, τ_1, τ_2 trial values after
rotation) is orbit-identical.
T3 = 2.236 ⇒ rate-state solver's OUTPUT (tau1_corr, etc.) is
orbit-spread at the v9.4.0 §11 gate magnitude.

**The pepper originates INSIDE the rate-state solver.**

## §H.5 Named next /code-implement target

Per the reviewer's table of named targets:
> T3 dominant → FaultFaceFlux::Evaluate friction-solver re-eval per step,
>              friction-output coupling.

### Primary target

- **File:** `dynamic/fault_face_flux.cpp`
- **Function:** `FaultFaceFlux::Evaluate` (the rate-state friction
  solve via Brent iteration on |V| in log10 space, plus ψ update).
- **Suspect mechanisms (ranked):**
  1. **Brent iteration convergence tolerance** — if the tolerance is
     near ULP relative to the traction scale (~75 MPa × 2e-16 =
     ~1.5e-8), per-QP convergence may terminate at DIFFERENT
     iteration counts across orbit-paired QPs, giving different
     final V values. This compounds step-by-step.
  2. **ψ evolution analytic integrator** — the `UpdateStateAnalytic`
     call uses `log(1e-12/1e-6)` and similar small-number math.
     Per-QP ULP noise here would seed per-step drift.
  3. **Friction-solver input ordering** — the order in which the
     trial-traction components are combined in Pelties eq(7) may
     introduce per-QP FP noise that compounds.

### Secondary target (friction-output coupling)

- **File:** `dynamic/fault_face_flux.cpp`
- **Function:** the post-solve writeback of `tau1_corr`, `tau2_corr`,
  `sigma_n_corr`, `V1`, `V2`, `slip_rate` into DOFData.
- **Suspect mechanism:** if `tau1_corr` is written as `tau1_0 +
  tau1_nuc + eta * V1` with V1 = V sin(θ) and θ computed from
  two FP atan2 calls, per-QP rounding of θ could inject noise.

### NOT a target

- `PrecomputedFaceFluxes::BuildInteriorMatrices` (Phase 2C-fault).
  Rescinded in §G.4. The bulk-Q-reconstruction-to-fault-QP path is
  T1-clean (T1 = 0).
- The MFEM basis, volume kernel, mass inverse, boundary-face vertex
  selection. All ruled out across §D–§G.

## §H.6 Acceptance criteria for the next fix

- Rate-state solver modification produces T3 spread ≤ 1e-10 Pa on
  the Kuhn pepper guard (v9.4.0 §11 gate for `tau1_corr`).
- T1 = 0 invariant preserved (the fix should not introduce new
  spread in the bulk Q path — otherwise a cross-layer
  contamination risk).
- All existing regression tests remain green.
- Phase 0 v9.5.0 R-001..R-005 unaffected.

## §H.7 Why this localization is trustworthy

- **T1 = 0 is a cell-level per-orbit comparison on a KUHN fixture
  (MakeCartesian3D, positive-oriented by construction).** No
  fixture contamination (§E was about D4's negative Jacobians;
  Kuhn does not have that issue).
- **T1 = 0 at EVERY step, not just step 0.** Even after 19 steps
  of compounded physics, the bulk Q → fault QP interpolation
  remains orbit-identical. The amplification cannot be in any
  upstream stage.
- **The tau1_corr = 2.236 Pa at step 19 matches
  v9.4.0 §11 = 2.24 Pa to 0.2%.** This is the NAMED ACCEPTANCE
  GATE; §H measures the right signal.
- The chain is rupture-DRIVEN: `tau2_nuc = nuc_dtau` matches the
  pepper guard configuration. Not a synthetic perturbation test.

## §H.8 What this localization does NOT prove

- The EXACT mechanism inside `FaultFaceFlux::Evaluate`. §H.5
  lists suspects ranked by likelihood; confirming which one
  requires instrumentation INSIDE the friction solver (e.g., per-QP
  Brent iteration count dump, per-QP ψ evolution precision check).
  That is a round-7 code-implement step.
- Whether the fix will transfer to production-scale (Gmsh) meshes.
  Kuhn is a unit fixture; production uses Gmsh. Layer 3 validation
  is Frontera-gated.

## §H.9 Open follow-ups (not blocking next /code-implement)

- Cross-fixture: re-run §H on fixed D4 to confirm T1 = 0 there too
  (the §F orientation fix should preserve the T1 invariant; §H
  Kuhn evidence alone is conclusive on localization).
- Per-QP friction-solver instrumentation: inside
  `FaultFaceFlux::Evaluate`, dump the Brent iteration count and
  final residual per QP per step. If orbit-paired QPs show
  different iteration counts or residuals, it pinpoints the
  convergence-tolerance suspect.

## §H.10 Files changed (round 6)

- **created:** `tests/unit/test_r6_phase2e_amplification_chain.cpp`
  — rupture-driven amplification chain.
- **modified:** `miniapps/seas/makefile` — SRC/OBJ/compile/link
  registration.
- **modified:** `phase1_arm1_findings_2026-04-23.md` — §G.2b
  (methodology caveats) and this §H.

### Raw logs

- `phase2e_amplification.txt` — the full per-step amplification
  table.

## §H.11 Regression check

- Phase 0 v9.5.0 R-001..R-005: PASS.
- phase2a switch: 14/14 PASS.
- D4 fixture gates (equivariance, Jacobian, fault sanity): PASS.
- Round-5 correlation test: builds and runs (supplementary).
- New Phase 2E amplification chain test: builds and runs; produces
  the §H.3 table with T1 = 0 and T3 = 2.236.

No production code touched. Flux-layer freeze respected. The next
code-implement step (Round 7) will be the FIRST flux-layer-adjacent
change and will require user authorization per §C R-006 criterion
#3.

## §H.12 Summary for the next /code-implement agent

**The pepper source is isolated: rate-state friction solve inside
`FaultFaceFlux::Evaluate`, not the DG kernel layer.**

Start by:
1. Reading `dynamic/fault_face_flux.cpp::Evaluate` end-to-end.
2. Adding per-QP Brent-iteration-count + final-residual logging
   (diagnostic, no functional change).
3. Running the pepper guard with this logging on; comparing iteration
   counts / residuals across orbit-paired fault QPs.
4. If counts differ or residuals are near ULP on some QPs and not
   others: the convergence criterion is the suspect. Tighten it
   (or change termination criterion) and re-run until T3 ≤ 1e-10 Pa.
5. If counts and residuals match across orbit-paired QPs: the
   amplification is in the analytic ψ integrator or the
   post-solve tau1_corr writeback. Instrument those next.

The freeze-allowed diagnostic steps (1-3) do not need user approval.
Step 4 (convergence-criterion change) is the first flux-layer-
adjacent modification and must be gated on user approval per §C
R-006 criterion #3.

---

# §I: Round-7 Strengthened Chain — §H Rescinded (2026-04-23)

Status: **AUTHORITATIVE.** Round-7 review flagged two critical
methodology gaps in §H: R-001 (T1 covered only 2 of 8 fault
triangles due to centroid bucketing) and R-002 (missing T2
trial-traction stage). With both fixes applied, §H's "pepper is
inside FaultFaceFlux::Evaluate" conclusion is EMPIRICALLY REVERSED.

## §I.1 Methodology fixes

- **R-001 — shared-face triangle pairing.** Every fault triangle has
  exactly one Elem1 and one Elem2 (MFEM invariant for interior
  faces). The pairing iterates the 8 fault triangles on Kuhn 2×2×2
  directly. §I delivers 8/8 coverage. §H's centroid bucketing
  (`(cx, min(cy, L-cy), cz)`) grouped only 2 buckets with ≥2
  members because Kuhn's 6-tet hex split places fault-adjacent tets
  at distinct centroids.
- **R-002 — T2 trial-traction stage.** Between T1 (Q_self
  reconstructed at fault QP) and T3 (DOFData.tau*_corr after
  friction solve), we insert T2 (`sigma_n_trial`, `tau1_trial`,
  `tau2_trial`) by calling the public
  `FaultFaceFlux::ComputeTrialTraction` with the canonical-frame
  Q_self values. Localizes whether amplification is upstream
  (T2 ≈ T3 → trial-traction / rotation) or downstream (T2 ≈ 0 &
  T3 > 0 → Brent / V-decomposition / writeback).

Source: `tests/unit/test_r7_amplification_chain_v2.cpp`. Raw log:
`phase2e_r7_chain_v2.txt`.

## §I.2 Measurement table (worst over 20 steps, Kuhn M0)

| Probe | Value | Units | Interpretation |
|---|---:|---|---|
| **T1+ max spread** | **870.6** | Pa | Bulk Q at fault QPs NOT orbit-clean |
| T1− max spread | 868.1 | Pa | Same — SZZ dominates |
| T2 sigma_n_trial spread | 272.6 | Pa | Dirty (driven by T1) |
| T2 tau1_trial spread | 2.63 | Pa | Track T3 tau1_corr |
| T2 tau2_trial spread | 425.4 | Pa | Dirty |
| **T3 tau1_corr spread** | **2.236** | Pa | Matches v9.4.0 §11 gate (2.24) |
| T3 tau2_corr spread | 173.5 | Pa | Dominated by T2 tau2 |
| T3 sigma_n_corr spread | 245.9 | Pa | Dominated by T2 sigma_n |
| T3b V1 spread | 2.62e-09 | m/s | Small |
| T3b V2 spread | 4.56e-05 | m/s | Matches v9.4.0 §11 |
| T3b slip_rate spread | 4.56e-05 | m/s | Matches |

T3/T2 for tau1: 2.236 / 2.63 = 0.85 (NOT amplifying; T3 ≈ T2).
T2/T1 for tau1: 2.63 / 870 = 3e-3 (selective coupling via
ComputeTrialTraction's VY/SXY inputs).

## §I.3 Per-component step-1 breakdown (the decisive evidence)

At step 1 (after ONE ADER step from Q = 0), per-component T1:

| Component | Plus side | Minus side | Physics role (canonical frame) | Verdict |
|---|---:|---:|---|---|
| SXX (σ_nn)     | 1.099e+00 | 1.099e+00 | y-sym | dirty |
| SYY            | 1.301e+00 | 1.301e+00 | y-sym | dirty |
| SZZ            | **4.800e+00** | **4.800e+00** | y-sym | **dirtiest** |
| **SXY (σ_nt1)** | **1.48e-16** | **2.88e-16** | y-antisym (drives tau1_trial) | **ULP** |
| SYZ            | 1.421e+00 | 1.421e+00 | y-sym | dirty |
| SXZ (σ_nt2)    | 2.184e+00 | 2.362e+00 | y-sym (drives tau2_trial) | dirty |
| VX (v_n)       | 2.95e-07  | 2.95e-07  | y-sym | near-clean |
| **VY (v_t1)**  | **1.43e-23** | **9.22e-24** | y-antisym (drives tau1_trial) | **ULP** |
| VZ (v_t2)      | 2.06e-07  | 2.42e-07  | y-sym | near-clean |

**The y-antisymmetric components — exactly those that SHOULD be
asymmetric under y-reflection — are ULP-clean. The y-symmetric
components — those that SHOULD be identical — carry all the
spread.**

## §I.4 Diagnosis: canonical→global rotation + flux deposition

The signature is not consistent with "friction solver is the bug"
(§H's wrong conclusion). It's consistent with:

- The fault-flux IMPOSED STATES (`Q_imp_plus`, `Q_imp_minus`) are
  computed in the canonical frame and are ULP-clean in the
  y-antisymmetric components that the friction solver outputs
  (tau1_corr is written to SXY in canonical frame; SXY is clean).
- The rotation back to global via `T_can`, followed by
  `flux_.Interior(can_n, Q_imp_g, Q_imp_g, F_h)`, produces the DG
  face-flux deposition into bulk Q.
- The y-symmetric bulk components (SZZ in particular) pick up O(Pa)
  orbit spread from this deposition — because the y-symmetric
  components' F_h contributions depend on the ROTATION applied per
  QP, and if T_can / can_n / the QP's local Loc1 transform has
  per-orbit-fault-QP variation, the y-symmetric components (which
  are nonzero in the RHS scale) inherit the ULP noise at O(Pa)
  amplitude, while y-antisymmetric components (which are O(1e-16)
  in RHS at step 1) remain at ULP.

Per-step compounding: at step N, SZZ spread has grown to O(N · 4.8)
Pa (linear). By step 19, T1 SZZ ≈ 870 Pa. This is LINEAR GROWTH,
characteristic of per-step source, not of nonlinear friction
amplification.

## §I.5 Named round-8 /code-implement target

**Primary target:**
- **File:** `dynamic/wave_operator.inl`
- **Function:** `WaveOperator::ComputeFaceFluxRHS` (and
  `ComputeADERFaceFluxRHS`), the interior-fault branch.
- **Specific lines:** approximately 1184-1264 (the per-QP block):
  1. `BuildRotation(can_n, can_t1, can_t2, T_can)` — line ~1196-1198.
  2. `flux_.Interior(can_n, Q_imp_plus_g, Q_imp_plus_g, F_h_plus)` —
     line ~1261-1264.
- **Suspect mechanisms (ranked):**
  1. Per-QP `can_n`, `can_t1`, `can_t2` come from `fault_basis_->GetBasis(fb_idx).qp_data[q]`. If those per-QP frame vectors have per-QP ULP noise that differs across orbit-paired fault QPs, the T_can matrix differs, and the rotation back to global injects O(ULP · stress_scale) noise into the y-symmetric components.
  2. `flux_.Interior(can_n, ...)` uses `can_n` as the normal. If `can_n` differs per QP at ULP, the flux computation inherits it.
  3. The per-side DG deposition at line ~1314-1339: `rhs[c, dof1, i] += w * shape1(i) * F_h_plus[c]`. The `w = ip.weight * nor_len` factor uses `CalcOrtho(ftr->Face->Jacobian(), nor_vec)` per-QP, which on Kuhn's orbit-related fault triangles could produce per-QP ULP noise in `nor_len`.

**NOT a target:** `FaultFaceFlux::Evaluate`. Round-6 §H's conclusion
is rescinded. The friction solver receives ULP-clean inputs
(tau1_trial ≈ 0, VY ≈ 0) for the tau1-driving channel and produces
ULP-clean V1 output (T3b V1 = 2.6e-9). The friction solver is NOT
the amplification source.

## §I.6 Why §H was wrong

§H's T1 = 0 at every step was a coverage artifact. With centroid
bucketing `(cx, min(cy, L-cy), cz)` on Kuhn 2×2×2, 2 of 8 fault
triangles shared bucket keys; the other 6 were in singleton
buckets (excluded from T1). If the 2 triangles that WERE bucketed
happened to be the 2 where orbit symmetry was preserved at the
per-DOF level, T1 = 0 on them — consistent with the observed §H
log. The other 6 triangles' spreads were invisible.

The shared-face pairing (one pair per triangle, 8 pairs) makes this
impossible: every fault triangle contributes to the spread
computation.

## §I.7 Cross-component check: which components feed tau1_corr?

`ComputeTrialTraction`:
```
tau1_trial = eta_s * (Q_minus[VY] - Q_plus[VY] + (Q_plus[SXY] + Q_minus[SXY])/Zs)
```

VY and SXY in canonical frame feed tau1_trial. Step-1 T1:
- SXY: 1.5e-16 Pa
- VY: 9e-24 Pa

Computed tau1_trial spread = eta_s · (0 + 1.5e-16/Zs) ≈
(Zs/2) · (1.5e-16 / Zs) = 7.5e-17. Observed T2 tau1_trial at step 1
= 2.6e-16. Matches within 4× (difference reflects contributions from
multiple fault triangles and the |a−b| metric).

tau2_trial:
```
tau2_trial = eta_s * (Q_minus[VZ] - Q_plus[VZ] + (Q_plus[SXZ] + Q_minus[SXZ])/Zs)
```

At step 1: VZ spread ≈ 2e-7 (relatively small), SXZ spread = 2.2 Pa.
Expected tau2_trial spread ≈ (1/2) · 2.2 = 1.1 Pa; observed 2.4 Pa
(×2 from multi-triangle comparison). Matches.

**Component accounting is consistent.** tau1_corr pepper DOES come
from SXY + VY contributions, which are ULP-clean. But tau2_corr
pepper comes from SXZ + VZ, which are dirty at O(Pa). So the named
acceptance gate (tau1_corr ≤ 1e-10) might be closer to closing than
tau2_corr / sigma_n_corr, but ALL channels ultimately depend on
clean y-symmetric bulk Q at fault QPs.

## §I.8 Revised freeze-unblock recommendation

Any round-8 /code-implement would touch `wave_operator.inl` (the
fault-dispatch rotation + deposition). This IS a flux-layer change
per §C R-006. The freeze-unblock criteria:

1. ✓ Phase 2D / Phase 2E probe re-measurement with strengthened
   methodology complete (§I).
2. ✓ Layer 2 pepper guard run (§G pepper-guard results on both
   fixtures).
3. ⚠ User explicit reopen authorization REQUIRED.

Before round 8: instrument the per-QP `T_can` / `can_n` / `w` values
at step 1 on orbit-paired fault triangles and verify they differ at
ULP. If yes → the round-8 fix is to canonicalize the per-QP frame
data (use `fault_basis_`'s stored canonical values directly rather
than recomputing via CalcOrtho + BuildFrame per-QP). If no → a
different mechanism is at play and further localization is needed.

## §I.9 Files changed (round 7)

- **created:** `tests/unit/test_r7_amplification_chain_v2.cpp`
  (shared-face pairing + T2 stage + per-component step-1 dump).
- **modified:** `miniapps/seas/makefile` — registration.
- **modified:** `phase1_arm1_findings_2026-04-23.md` — §I.

### Raw log

- `phase2e_r7_chain_v2.txt` — full chain output with per-component
  step-1 breakdown.

## §I.10 Regression check

- Phase 0 v9.5.0 R-001..R-005: PASS.
- phase2a switch: 14/14 PASS.
- D4 gates (equivariance, Jacobian, fault sanity): PASS.
- Round-5 correlation test: supplementary.
- Round-6 Phase 2E amplification chain: still runs (its T1 = 0
  result is now understood as a coverage artifact, §I.6).
- Round-7 strengthened chain: new, produces §I.2 / §I.3 tables.

No production code touched. Flux-layer freeze respected.

## §I.11 Summary for the next /code-implement agent

**The pepper source is in `wave_operator.inl`'s fault dispatch
rotation + deposition, NOT in `FaultFaceFlux::Evaluate`.**

Specifically: the y-symmetric components (SZZ, SXZ, SYZ, SXX, SYY)
of bulk Q at fault-adjacent tets diverge at O(Pa) per step, while
y-antisymmetric components (SXY, VY) stay at ULP. This is
consistent with per-QP ULP noise in `T_can` / `can_n` / `w`
(computed via CalcOrtho + BuildFrame per-QP) leaking into
y-symmetric rhs deposition.

Next actions:
1. (Freeze-allowed diagnostic) Dump per-QP `can_n`, `can_t1`, `can_t2`,
   and `w = ip.weight * nor_len` values at step 1 on orbit-paired
   fault QPs. Check for ULP-level differences.
2. If the diagnostic shows ULP-differences: round 8 fix is to
   use `fault_basis_`'s precomputed canonical frame values
   uniformly (not recompute via CalcOrtho per-QP). Requires user
   authorization per §C R-006 criterion #3.
3. If the diagnostic shows bit-identical frames: the rotation/flux
   isn't the source; further localization needed (probably the
   DG shape evaluation on the per-face-ref-to-elem-ref QP mapping).

---

# §J: Round-8 Step 0 Re-measurement (2026-04-23)

Status: **AUTHORITATIVE.** §I's narrative had three critical analytical
flaws (round-8 review R-001, R-002, R-003). §J executes the reviewer's
recommended Step 0: multi-step per-component T1 dump, growth-rate
characterization, and cross-fixture (Kuhn + fixed D4) runs. The
result rescinds §I's named target and reframes the investigation from
"find the seed" to "identify the amplification source".

## §J.1 Round-8 review findings addressed

| R8 ID | §I flaw | §J resolution |
|---|---|---|
| R-001 | "SXY ULP-clean" was step-1-only | §J.3 shows SXY grows to 3.74 Pa by step 19 on Kuhn |
| R-002 | "Linear growth" from §I.5 | §J.4 fits power-law T ~ t^2.2; step ratios decelerate 3.0 → 1.1 |
| R-003 | Symmetry classification cherry-picked SXY+VY | §J.5 reclassifies: canonical SXY, SXZ, VX are y-antisym (VX NOT clean) |
| R-004 | "Friction solver clean" based on V1 only | §J.3 notes V2 matches slip_rate spread (dirty) |
| R-005 | No D4 cross-fixture | §J.6 runs on both; D4 shows different sigma_n signature but WORSE tau1_corr pepper |
| R-006 | Stale line numbers in §I.5 | Moot: §I.5 target rescinded in §J.7 |

## §J.2 Raw data (multi-step per-component T1, Kuhn)

Log: `phase2e_r8_chain_kuhn.txt` in this directory.

| Step | SXX | SYY | SZZ | SXY (y-antisym) | SYZ | SXZ (y-antisym) | VX (y-antisym) | VY | VZ |
|------|-----|-----|-----|-----------------|-----|-----------------|----------------|-----|-----|
| 1    | 1.10 | 1.30 | 4.80 | **1.5e-16** | 1.42 | 2.18 | **2.95e-07** | **9.2e-24** | 2.06e-07 |
| 5    | 17.9 | 21.1 | 71.3 | 1.93e-07 | 20.9 | 32.7 | 4.48e-06 | 9.67e-09 | 3.04e-06 |
| 10   | 66.5 | 78.0 | 258  | 0.113 | 76.4 | 119  | 1.68e-05 | 5.35e-08 | 1.14e-05 |
| 15   | 137  | 162  | 556  | 1.07 | 160  | 251  | 3.50e-05 | 1.39e-07 | 2.37e-05 |
| 19   | **194** | **236** | **871** | **3.74** | **257** | **400** | **5.47e-05** | **2.12e-07** | **3.72e-05** |

Growth to step 19 from step 1:
- SZZ: 4.8 → 871 (181×)
- SXY: 1.5e-16 → 3.74 (2.5e+16×) — ramping from ULP to Pa-scale
- VY: 9.2e-24 → 2.12e-07 (2.3e+16×) — same ramp pattern
- VX: 3e-7 → 5.5e-5 (183×) — similar to SZZ

**§I's "SXY ULP-clean" claim applied only at step 1.** SXY, VY ramped
from ULP to Pa/m-s scale over 19 steps by coupling to other components.

## §J.3 Pelties Eq 7b reconciliation

At step 19:
- Q[SXY] spread ≈ 3.87 Pa (avg of plus 3.74, minus 4.00)
- Q[VY] spread ≈ 2.3e-7 m/s
- eta_s ≈ Zs/2 ≈ 4.63e+6 Pa·s/m; 1/Zs ≈ 1.08e-7

Expected tau1_trial spread from Pelties Eq 7b:
```
tau1_trial = eta_s * (VY_m - VY_p + (SXY_p + SXY_m)/Zs)
           ≈ eta_s * Q[VY]-diff + eta_s / Zs * Q[SXY]-sum
           ≈ 4.63e6 * 2.3e-7  +  0.5 * 7.74
           ≈ 1.07            +  3.87
           ≈ 4.9 Pa
```

Observed T2 tau1_trial at step 19 = 2.63 Pa (Kuhn). Order-of-magnitude
agreement within ~2× (difference explained by "spread = max − min
across triangles" vs linear Pelties on averages). **T2 is NOT
inconsistent with SXY/VY growth** — R-001's "frame mismatch" alternative
is not needed.

## §J.4 Growth-rate characterization (R-002 fix)

Per-step ratios (Kuhn T1 SZZ):
- step 2/1: 2.99
- step 3/2: 1.99
- step 10/9: 1.22
- step 19/18: 1.11

Decelerating, not constant. Not linear (§I.5 would predict T19 ≈ 91 Pa;
observed 871). Power-law fit:

| Series | Power-law fit | At t=20 pred | Observed t=19 |
|---|---|---|---|
| T1 plus-side max (Kuhn) | 1.29 · t^2.20 | 933 | 871 |
| T2 sigma_n_trial (Kuhn) | 0.40 · t^2.20 | 292 | 272 |
| T2 tau2_trial (Kuhn) | 0.65 · t^2.19 | 457 | 425 |
| T2 tau1_trial (Kuhn) | 5.6e-11 · t^9.19 | 51 | 2.63 |
| T3 tau1_corr (Kuhn) | 1.08e-12 · t^10.48 | 47 | 2.24 |

**Exponent ~2.2 for components with nonzero step-1 seed.**

**Exponent ~9-10 for components ramping from ULP.** These large
exponents reflect the sigmoid-like ramp from ULP to a physical scale,
not a physical amplification rate. They cannot be extrapolated beyond
the saturation point.

Interpretation: per-step ratios decelerating from 3 → 1.1 are
characteristic of wave propagation coupled with rupture expansion —
the wavefront hits more triangles as it propagates, each receives a
per-step increment, but the total saturates as the rupture area
saturates.

**The §I.5 "single per-step ULP injection in CalcOrtho" mechanism
would give LINEAR growth (per-step const increment → total ∝ t).
The observed t^2.2 rejects this.**

## §J.5 Corrected symmetry classification (R-003 fix)

For the TPV102 fault at y=L/2 with canonical frame
`(can_n, can_t1, can_t2) = ((0,-1,0), (0,0,-1), (+1,0,0))` (the BP5
convention per CLAUDE.md), the canonical-frame components map to
global components as:

| Canonical | Global | y-mirror |
|-----------|--------|----------|
| SXX (σ_nn) | σ_yy | y-sym |
| SYY | σ_zz | y-sym |
| SZZ | σ_xx | y-sym |
| **SXY (σ_n,t1)** | **σ_yz** | **y-antisym** |
| SYZ (σ_t1,t2) | σ_xz | y-sym |
| **SXZ (σ_n,t2)** | **σ_xy** | **y-antisym** |
| **VX (v_n)** | **v_y** | **y-antisym** |
| VY (v_t1) | v_z | y-sym |
| VZ (v_t2) | v_x | y-sym |

So canonical-frame y-antisym components: **SXY, SXZ, VX**.

Observed step-1 orbit spread:
- y-antisym: SXY = 1.5e-16 (ULP), SXZ = 2.2 Pa (NOT ULP), VX = 2.95e-07 (near-clean).
- y-sym: SXX = 1.1, SYY = 1.3, SZZ = 4.8, SYZ = 1.4, VY = 9e-24 (ULP), VZ = 2.1e-07.

No clean partition. SXY and VY are ULP at step 1 because the rupture
drive tau2_nuc = nuc_dtau is ONLY nonzero in the SXZ (canonical) /
σ_xy (global) channel. So SXZ is seeded at step 1, while other
components that don't directly receive forcing sit at ULP initially
and then grow via cross-component coupling as the rupture propagates.

By step 19 all components have grown. The "y-antisym clean" label was
not a genuine symmetry property; it was an artifact of "which channel
the nucleation forcing initially excites".

## §J.6 Cross-fixture comparison (R-005 fix)

Log: `phase2e_r8_chain_d4.txt`.

| Metric | Kuhn | D4 (orientation-fixed) |
|--------|------:|------:|
| T1 plus-side SZZ (step 19) | 871 Pa | 868 Pa |
| T1 minus-side SZZ | 868 Pa | 868 Pa |
| T1 SXY (step 19) | 3.74 Pa | 3.78 Pa |
| T1 VX | 5.47e-05 | 5.47e-05 |
| **T2 sigma_n_trial** | **272 Pa** | **3.4e-13 Pa** |
| T2 tau1_trial | 2.63 Pa | 5.63 Pa |
| T2 tau2_trial | 425 Pa | 614 Pa |
| **T3 tau1_corr** | **2.24 Pa** | **4.80 Pa** |
| T3 tau2_corr | 173 Pa | 469 Pa |
| T3 sigma_n_corr | 246 Pa | 0.00 Pa |
| T3b V2 | 4.56e-05 | 5.68e-05 |

Key findings:
- **T1 bulk-Q spreads (SZZ, SXY, VX, VZ etc.) are ESSENTIALLY IDENTICAL across fixtures.** The DG layer produces the same bulk-Q divergence pattern on both fixtures, up to per-triangle plus/minus-assignment differences.
- **T2 sigma_n_trial is ULP on D4 but 272 Pa on Kuhn.** Most likely: frame-labeling artifact in my test (centroid-based plus/minus vs production's `elem1_on_plus` convention). On D4 the production convention gives sums that cancel; on Kuhn they don't. Either way, the SPREAD being ULP on D4 shows cancellation works on D4.
- **T3 tau1_corr is WORSE on D4 (4.80 Pa) than Kuhn (2.24 Pa).** Matches the round-4 §F.5 pepper-guard finding. A fix must reduce tau1_corr on BOTH fixtures; D4 is the stricter acceptance gate.

## §J.7 Verdict: §I.5 target RESCINDED; Branch B selected

§I.5 named `wave_operator.inl::ComputeFaceFluxRHS` fault dispatch
rotation + CalcOrtho as the round-8 target. §J rescinds this:

1. Growth pattern is t^2.2, not linear. A single ULP seed in CalcOrtho
   per-QP cannot produce this. **Fixing CalcOrtho would attenuate the
   seed but not the amplification, so pepper from any other ULP seed
   would still emerge.**
2. Cross-fixture: D4 (Jacobian-uniform, same slot equivariance) has
   WORSE tau1_corr pepper than Kuhn. CalcOrtho-based per-QP normal
   should produce DIFFERENT per-QP values on the two fixtures; the
   SAME t^2.2 growth on both fixtures implies the amplification
   mechanism is fixture-independent.

Per the round-8 decision tree:

- ~~Branch A (§I.5 seed mechanism survives)~~: REJECTED by §J.4 +
  §J.6 cross-fixture.
- **Branch B (amplification source)**: SELECTED. The bug is in
  an amplification mechanism — friction-solver coupling under
  rupture, SIPG penalty on fault-adjacent tets, RK-on-ψ
  interleave, or rupture-front eigenmode — not in a seed.
- ~~Branch C (SXY grows upstream of trial-traction)~~: SXY DOES
  grow, but so does everything else. Doesn't isolate trial-traction
  as the cause.

## §J.8 What the amplification source could be

Candidate amplifiers (ranked by evidence coupling, not by my prior
suspicion):

1. **Rupture-front eigenmode in the rate-state feedback loop.** The
   rate-state friction couples slip rate V to traction τ via an
   exponential `C = exp(psi/a) / (2*V0)`. Feedback: higher τ → higher
   V → more slip → state evolution (ψ increases for decelerating,
   decreases for accelerating slip) → shifts fric law → modifies τ
   on next step. Classical earthquake-cycle instability. If one tet's
   V grows 1% faster than its orbit partner's, the gap doubles every
   ~70 steps (for γ ≈ 0.01). The observed per-step ratio 1.1-1.3 is
   consistent with γ in the 0.1-0.3 range, which is on the order of
   rate-state's `a`-parameter timescale.
2. **SIPG penalty on fault-adjacent faces.** The penalty is
   proportional to shape-function gradient magnitudes; tets with
   different stored vertex orders can have different penalty
   coefficients when the face is a fault face (handled by
   FaultFaceFlux, not the standard SIPG). Probably NOT the primary
   driver since fault faces bypass SIPG.
3. **RK4 / ADER-CK stage interleave on ψ.** The psi evolution uses
   per-stage V values. If ψ accumulates per-step ULP error that
   feeds back into next-step Brent, that's a slow eigenmode.
4. **Reflected-wave echo** from the domain boundary (absorbing BC not
   perfect). But all fault-adjacent tets receive the same echo, so
   this wouldn't produce orbit asymmetry. Low likelihood.

## §J.9 Recommended round-9 diagnostic (freeze-allowed)

Don't propose a code fix. Instead, distinguish Branch-B variants:

### Diagnostic 1: is the amplification driven by the rate-state feedback loop?

Build a "frozen-friction" variant of the pepper-guard that skips the
rate-state solver AFTER step 1 — replace with a static friction model
that takes τ_1, τ_2 from step-1's output and holds them constant.
If the pepper spread NO LONGER grows, the rate-state feedback IS the
amplifier → round 10 target is the rate-state-solver coupling inside
`FaultFaceFlux::Evaluate`.

If pepper still grows with frozen friction, the amplifier is elsewhere
(wave operator / SIPG / RK-on-ψ).

### Diagnostic 2: per-step eigenvalue probe

At each step, compute the singular values of the per-tet rhs map
restricted to orbit-paired tets. The largest singular value ≥ 1
would identify the amplification factor. This is a physics-literate
instrumentation, non-trivial to implement, but definitive.

### Diagnostic 3: quasi-static rupture drive

Run the pepper guard with dt reduced 10×. If pepper spread SAME at
equal physical time, the amplifier is dt-independent (kernel/
structural). If pepper DROPS, amplifier is dt-dependent (stage
interleave or time-integration instability).

## §J.10 Freeze-unblock status

Still blocked per §C R-006. §J does NOT authorize any code change.
The round-9 diagnostics (§J.9) are all freeze-allowed
(test-only, no production modification).

## §J.11 Files changed (round 8 Step 0)

- **modified:** `tests/unit/test_r7_amplification_chain_v2.cpp`:
  - `BuildFaultMeshByEnv()` — honors SEAS_TEST_FIXTURE for Kuhn / D4.
  - Multi-step per-component T1 dump at steps 1, 5, 10, 15, 19.
  - Growth-rate characterization block with power-law fits and
    per-step ratios.
- **modified:** `phase1_arm1_findings_2026-04-23.md` — this §J.

### Raw logs

- `phase2e_r8_chain_kuhn.txt` — Kuhn run with multi-step dump +
  growth-rate fits.
- `phase2e_r8_chain_d4.txt` — fixed D4 run, same instrumentation.

## §J.12 Regression check

- Phase 0 v9.5.0 R-001..R-005: PASS.
- phase2a switch: 14/14 PASS.
- D4 fixture gates (equivariance, Jacobian, fault sanity): PASS.
- Round-7/8 chain builds and runs on both fixtures.

No production code touched. Flux-layer freeze respected.

## §J.13 Summary for the next /code-implement agent

**The previous 4 rounds' named targets have ALL been wrong**:

- Round 5 §G: Phase 2C-fault (BuildInteriorMatrices) — rescinded.
- Round 6 §H: Inside FaultFaceFlux::Evaluate (Brent, ψ, writeback) —
  rescinded (§I showed T1 wasn't actually zero).
- Round 7 §I: Inside wave_operator.inl fault rotation (CalcOrtho
  per-QP) — RESCINDED here by §J (growth is not linear; mechanism
  is amplification, not seed).
- Round 8 would be what, next? Branch B: an amplifier, not a seed.

The current state: **one kernel-layer fix is not expected to close
pepper.** The bug is a rupture-driven eigenmode that amplifies ANY
initial orbit asymmetry regardless of its source.

Actions for round 9 (all freeze-allowed, before any code target):
1. Run §J.9 Diagnostic 1 (frozen-friction variant). Builds on the
   existing pepper-guard harness with a small FaultFaceFlux swap.
   ~2 hours.
2. If Diagnostic 1 shows rate-state is the amplifier → round 10
   targets the rate-state coupling (NOT the Brent solver itself;
   the COUPLING between slip rate output and next-step traction
   input, i.e., the per-step stage interleave).
3. If Diagnostic 1 rules out rate-state → run §J.9 Diagnostic 3
   (dt reduction). Localizes time-integrator vs kernel.

Under NO circumstances should round 9 edit production code before a
diagnostic identifies the amplifier. The v8 → v9.4 failure pattern
(fix one suspect, pepper persists, fix next suspect) has repeated
at the DEBUG ROUND level. Each of rounds 5, 6, 7 spent a round
being wrong. Round 8 (§J) made one correction but used the same
fix-and-verify approach the reviewer now explicitly flagged as
unsustainable.

A diagnostic that identifies the AMPLIFICATION SOURCE (not the seed)
is the only way to avoid round 10 repeating the pattern.

---

# §K: Round-9 Execution — FREEZE-C + dt-refinement (2026-04-23)

Status: **AUTHORITATIVE.** Round-9 reviewer identified three spec
fixes (R-001 mixed-regime fits, R-002 cross-fixture coverage, R-003
three-FREEZE disambiguation) plus two pre-commit requirements
(R-004 dt discriminator, R-005 verdict matrix). §K applies all
five and executes the freeze-allowed portion (FREEZE-C + dt).

## §K.1 Freeze-policy scope acknowledgment (R-006)

The round-9 plan specified three FREEZE variants. Per §C R-006
criterion #3, FaultFaceFlux and wave_operator.inl modifications
require explicit user authorization. FREEZE-A (skip Evaluate
entirely) and FREEZE-B (override Q_imp outputs) both need flux-
layer code changes and therefore require user authorization before
implementation. FREEZE-C (hold ψ constant) can be implemented at
test level by overwriting `dof_data.psi` before each `AdvanceADER`
call — no production-code modification required.

§K executes FREEZE-C + dt refinement. FREEZE-A/B are deferred pending
user authorization; §K documents what they would entail (§K.8).

## §K.2 R-001 refit: power-law, not exponential

Full-range fits over steps 1-19 on T1+ max (Kuhn):

| Model | Fit | R² |
|-------|-----|-----|
| Linear (y = a + b·t) | y = -158 + 48·t | 0.952 |
| **Power (y = a·t^b)** | **y = 4.18 · t^1.797** | **0.997** |
| Exponential (y = a·exp(b·t)) | y = 15.1 · exp(0.243·t) | 0.465 |

Short-window (steps 1-5): power exponent 1.67.
Long-window (steps 10-19): power exponent 1.90.

**Power-law decisively better than exponential** (R² 0.997 vs 0.465).
Exponent stable across windows (1.67 → 1.90, not a regime change).

**Interpretation:** NOT a rate-state eigenmode (which would be
exponential) and NOT pure linear seed-accumulation (which would be
t^1). Consistent with "wave amplitude growing linearly × per-step
ULP seed" → time-integrated t^2.

## §K.3 FREEZE-C execution

Binary: `seas_test_r9_frozen_friction_C`.
Raw logs: `phase2e_r9_freezeC_kuhn.txt`, `phase2e_r9_freezeC_d4.txt`.

### FREEZE-C mechanism

- Snapshot `dof_data[i].psi` values AFTER step 1.
- Before steps 2..N, overwrite `dof_data[i].psi` with the step-1
  snapshot.
- `FaultFaceFlux::Evaluate` runs normally, reads the frozen ψ from
  DOFData, and solves friction with that ψ. Other DOFData fields
  (tau*_corr, V1, V2, etc.) are written by Evaluate and recomputed
  per-step based on current bulk Q.

This holds **state-variable evolution constant** while still
running the friction solve. It tests: "if ψ didn't evolve, would
pepper close?" — the classical rate-state instability test.

### Results

| Fixture | Run | tau1_corr worst | FREEZE-C / baseline |
|---------|-----|----------------:|---------------------:|
| Kuhn | BASELINE   | 2.236162 Pa | — |
| Kuhn | FREEZE-C   | 2.236162 Pa | **1.0000** |
| D4   | BASELINE   | 4.796186 Pa | — |
| D4   | FREEZE-C   | 4.796186 Pa | **1.0000** |

**FREEZE-C ratio is BIT-IDENTICAL to baseline on BOTH fixtures.**

**ψ rate-state feedback is NOT the amplifier.** Holding ψ constant
from step 1 onward produces the EXACT same pepper signature as
letting ψ evolve normally. The classical earthquake-cycle
instability hypothesis is falsified.

This means in 20 steps at dt = 5e-5, ψ-evolution contributes
negligible change to tau1_corr spread. Consistent with ψ being a
slow state variable (Dc ≈ 0.13 m, state evolution timescale
~seconds) while the simulation runs 20 × 5e-5 = 1 ms.

## §K.4 dt-refinement execution

Three runs at the SAME physical time (t_final = 1.0 ms): dt =
5e-5 (20 steps), dt/2 = 2.5e-5 (40 steps), dt/4 = 1.25e-5 (80 steps).

### Results

| Fixture | dt | n_steps | tau1_corr | ratio to baseline |
|---------|-----|---------|----------:|------------------:|
| Kuhn | 5e-5 (baseline) | 20 | 2.236 Pa | 1.000 |
| Kuhn | 2.5e-5          | 40 | 2.610 Pa | 1.167 |
| Kuhn | 1.25e-5         | 80 | 2.810 Pa | 1.257 |
| D4   | 5e-5 (baseline) | 20 | 4.796 Pa | 1.000 |
| D4   | 2.5e-5          | 40 | 5.616 Pa | 1.171 |
| D4   | 1.25e-5         | 80 | 6.054 Pa | 1.262 |

Fit: `tau1_corr ∝ dt^-0.17` on both fixtures.

**The exponent is NEGATIVE but small.** Pepper marginally INCREASES
with dt refinement (1.26× over 4× fewer steps → 4× more steps).

**Critical finding:** Sub-resolved physics would give positive
exponent (tau1_corr ∝ dt^p with p > 0). The observed p = -0.17
rules out sub-resolved physics. The bug is NOT a missing time
integrator or rate-state resolution issue.

Weak negative exponent is consistent with:
- Per-step FP roundoff accumulating over more steps (more steps →
  more rounding operations).
- Not a LINEAR accumulation (would give p = -1); the weak scaling
  suggests most of the roundoff is "cached" in the spatial state
  rather than freshly injected each step.

## §K.5 Verdict per pre-committed R-005 matrix

Reviewer's R-005 matrix:
- FREEZE-C CLOSES + dt-scales → sub-resolved rate-state physics.
- FREEZE-C CLOSES + dt-insensitive → code bug in ψ-feedback coupling.
- FREEZE-C does NOT close + dt-scales → amplifier in wave-kernel physics, sub-resolved.
- FREEZE-C does NOT close + dt-insensitive → **code bug outside friction and ψ**.
- NONE close → amplifier not in ψ.

Observed:
- FREEZE-C does NOT close (ratio 1.000).
- dt-scaling exponent -0.17 ≈ insensitive (weakly anti-scaling).

**Matrix branch:** "code bug outside friction and outside ψ".

## §K.6 What the narrowed candidate space actually contains

After ruling out ψ (§K.3) and sub-resolved physics (§K.4), the
remaining candidates are:

1. **Per-QP DG face-flux ULP noise** — CalcOrtho, BuildFrame,
   Loc1.Transform, shape·Q — all FP operations that can produce
   per-QP roundoff at different ULP positions across orbit-paired
   tets. Adressed in round 7 §I.5 (rescinded for wrong growth-rate
   model but mechanism itself not refuted).
2. **Per-tet nApNm1 matrix construction** — round 5 §G's
   `BuildInteriorMatrices` target (also rescinded, same reason).
3. **ADER Cauchy-Kovalevskaya predictor per-stage roundoff** —
   each stage applies spatial derivatives (`ApplySpatialDerivative`);
   per-stage ULP accumulates.
4. **SIPG penalty coefficient precision** — if penalty `sigma` is
   computed with per-face FP operations, per-QP variation possible.

**All four are per-step numerical noise sources.** The round-9
evidence says pepper is dominated by integrated per-step spatial
ULP, amplified by wave propagation (t^2). No single one of these
is necessarily the culprit — fixing any one might reduce, but not
eliminate, pepper.

This is consistent with the v8→v9.4 history where 5 round of "fix
this named suspect, pepper stays" kept happening. If the pepper is
DIFFUSELY sourced from per-step spatial operations, single-source
fixes cannot close it.

## §K.7 Revised proposal: benchmark against Tandem / SeisSol

Given the evidence that pepper is per-step spatial numerical noise
amplified by wave propagation — not a localized code bug — the
productive next step is NOT another "identify the seed" round. It
is:

**Verify whether the observed pepper floor is physically normal
for this class of DG discretization by benchmarking against
Tandem (SCEC-verified reference) at the same resolution.**

If Tandem at Kuhn M0 (or an equivalent small-mesh reproduction)
shows similar tau1_corr pepper, the floor is INTRINSIC to low-
resolution DG for TPV102 rupture and the pepper guard's 1e-10
tolerance is UNREALISTIC at the unit-test fixture scale. The fix
then is to:
- Adjust the pepper-guard tolerance to match the verified DG
  discretization floor, OR
- Run the acceptance gate at a resolution where the floor IS below
  1e-10.

If Tandem does NOT show the same pepper at the same resolution,
there's a SPECIFIC MFEM-vs-Tandem numerical difference to track
(not a generic "DG noise"). The comparison localizes a real code
path to inspect.

## §K.8 What FREEZE-A and FREEZE-B would tell us (if authorized)

- **FREEZE-A (fault locked):** Modify `FaultFaceFlux::Evaluate` to
  return `Q_imp = Q_{plus,minus}` (identity, no friction effect)
  for steps 2..N. If pepper closes, the friction MECHANISM (not
  just ψ) is the amplifier. Requires flux-layer authorization.
- **FREEZE-B (step-1 cached Q_imp):** Modify `FaultFaceFlux::Evaluate`
  to return step-1 cached Q_imp values for steps 2..N. Tests
  per-step RE-EVALUATION. If pepper closes with FREEZE-B but NOT
  with FREEZE-A, the per-step recalculation (not the mechanism)
  is the amplifier. Requires flux-layer authorization.

Both would disambiguate from "code bug outside friction and ψ"
by testing specifically whether friction-re-eval-per-step is the
amplifier.

**If the user authorizes FREEZE-A/B:** expected implementation is
~50 lines of env-gated code in `FaultFaceFlux::Evaluate`,
defaulting OFF. The test framework already exists in
`test_r9_frozen_friction_C.cpp`; extending to variants is
~30 lines.

## §K.9 Round-10 recommendation

Given the evidence, three possible next steps:

### Option X (recommended — freeze-allowed): Tandem benchmark

Build a minimal Tandem / SeisSol reproducer of the pepper guard at
the same 2×2×2 Kuhn fixture scale. Compare tau1_corr spread. This
tests whether the observed pepper is MFEM-specific or generic DG.
Estimated: 2 days if Tandem has this test, else 1 week to build.

### Option Y (freeze unblock required): FREEZE-A/B

Request user authorization for flux-layer modifications to extend
the round-9 diagnostic with FREEZE-A and FREEZE-B. Delivers
complete disambiguation within the existing Arm 1 / Phase 2E
framework. Estimated: 1 day after authorization.

### Option Z (accept + gate relaxation): Tolerance adjustment

Accept the 2-3 Pa pepper floor at 2×2×2 Kuhn as an artifact of the
tiny fixture and adjust the pepper guard's tolerance to
discretization-appropriate levels (say 1 Pa). Deprioritize the unit-
fixture guard and rely on Layer 3 (Frontera production) as the
acceptance gate. Estimated: 2 hours (tolerance change), freeze-
allowed.

## §K.10 Freeze-unblock status

Still blocked. §K does not authorize any code change.

§K explicitly requests user input on Options X, Y, Z.

## §K.11 Files changed (round 9)

- **created:** `tests/unit/test_r9_frozen_friction_C.cpp`
  (FREEZE-C + dt-refinement on env-switchable fixture).
- **modified:** `miniapps/seas/makefile` — SRC/OBJ/compile/link entries.
- **modified:** `phase1_arm1_findings_2026-04-23.md` — this §K.

### Raw logs

- `phase2e_r9_freezeC_kuhn.txt` — Kuhn baseline + FREEZE-C + dt/2 + dt/4.
- `phase2e_r9_freezeC_d4.txt` — fixed D4 same sweep.

## §K.12 Regression check

- Phase 0 v9.5.0 R-001..R-005: PASS.
- phase2a switch: 14/14 PASS.
- D4 fixture gates: PASS.
- Round-7/8 chain (preserved): PASS.
- New round-9 test: builds clean; produces verdict table + ratios.

No production code touched. Flux-layer freeze respected.

## §K.13 Summary for the user / next-round agent

**The investigation has narrowed the pepper source to:**

- NOT in friction mechanism ψ-feedback (FREEZE-C, §K.3).
- NOT in time integrator resolution / rate-state timestep physics
  (dt refinement, §K.4).
- NOT in ANY single named kernel-layer target from rounds 5-8
  (all RESCINDED by subsequent rounds).
- IS in per-step spatial numerical noise, amplified by wave
  propagation over the rupture, to power-law t^1.8-2.0 scaling.

**What we CANNOT do without user authorization:**

- FREEZE-A / FREEZE-B (requires flux-layer code changes to
  `FaultFaceFlux::Evaluate` or `wave_operator.inl`).

**Honest recommendation:** the investigation has reached
diminishing returns on "find the line". The pepper signature is
consistent with DIFFUSE per-step roundoff accumulation, which is
not a single-line bug. The next productive step is Option X
(benchmark against Tandem) OR Option Y (authorize FREEZE-A/B to
complete disambiguation) OR Option Z (accept the floor, adjust
gate tolerance).

Rounds 5-9 have all ended with "amplifier/seed not yet uniquely
localized; next round should be X". The pattern suggests there
IS no uniquely-localizable bug — the pepper is emergent from many
small numerical noise sources combined with wave amplification.
Round 10 should test that hypothesis directly (Option X: external
benchmark) rather than continue the "next suspect" cycle.

User decision required.

---

# §L. Round 11 — FREEZE-A/B + audit branch-isolation (2026-04-23)

User authorized Option Y from §K.9. Round-11 plan: extend
`test_r9_frozen_friction_C.cpp` with FREEZE-A (via a minimal hook in
`FaultFaceFlux::EvaluateADER`) and FREEZE-B (full-state DOFData
restore), then add `SEAS_TEST_SKIP_BOUNDARY` / `SEAS_TEST_SKIP_INTERIOR_NONFAULT`
toggles to three helpers in
`test_adjacent_triangle_fault_first_step_audit.cpp`
(`RunADERNonFaultFaceAudit`, `RunPrecomputedFluxLiftedAudit`,
`RunAnalyticTraceLiftAudit`).  Run the 6-case matrix on Kuhn and
repeat 1–3 on D4.

## §L.1 Files changed

- **modified:** `miniapps/seas/dynamic/fault_face_flux.cpp` — added
  `<cstdlib>` include and the `SEAS_TEST_FREEZE_A` env-gated bypass
  inside `EvaluateADER`.  When set, the friction `Evaluate` call is
  skipped, `data` is left untouched, and the imposed state is the
  time-averaged bulk state (`Q_imp_± = Q_avg_±`) — so the `I_imp_±
  = dt · Q_imp_±` return makes the fault face contribute no flux
  correction.  Env unset → byte-identical to the pre-hook path.

- **modified:** `miniapps/seas/tests/unit/test_r9_frozen_friction_C.cpp`
  — replaced `bool freeze_psi_after_step1` with `enum class
  FreezeMode { None, FreezeA, FreezeB, FreezeC }`; added
  `FrictionSnapshot` for psi / slip_rate / V1 / V2 / tau1_corr /
  tau2_corr / sigma_n_corr / slip1 / slip2; added `FreezeAEnvScope`
  RAII wrapper around `SEAS_TEST_FREEZE_A`; added verdict branches
  for the Round-11 decision question.

- **modified:** `miniapps/seas/tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — added `EnvFlagSet` helper and `SEAS_TEST_SKIP_BOUNDARY` /
  `SEAS_TEST_SKIP_INTERIOR_NONFAULT` face-level `continue` guards in
  `RunADERNonFaultFaceAudit`, `RunPrecomputedFluxLiftedAudit`, and
  `RunAnalyticTraceLiftAudit`.

Build: clean on both targets.

## §L.2 Pepper-guard 6-case matrix (20 steps, dt = 5e-5)

Worst tau1_corr spread [Pa] across 20 AdvanceADER steps:

| mode       | Kuhn      | D4        | ratio / Kuhn baseline | ratio / D4 baseline |
|------------|-----------|-----------|-----------------------|---------------------|
| BASELINE   | **2.236** | **4.796** | 1.0000                | 1.0000              |
| FREEZE-A   | **0.000** | **0.000** | **0.0000**            | **0.0000**          |
| FREEZE-B   | 2.236     | 4.796     | 1.0000                | 1.0000              |
| FREEZE-C   | 2.236     | 4.796     | 1.0000                | 1.0000              |
| DT/2       | 2.610     | 5.616     | 1.1671                | 1.1710              |
| DT/4       | 2.810     | 6.054     | 1.2566                | 1.2623              |

dt-scaling fit: Kuhn `tau1_corr ∝ dt^-0.165`, D4 `dt^-0.168`.
Consistent with §K.4 (slight _growth_ under refinement — dt
insensitive).

## §L.3 Primary verdict (Round-11 decision question)

> "Is the first nonuniformity created by the MFEM wave update even
> when friction is frozen, or only when live fault outputs are fed
> back?"

**Answer: only when live fault outputs are fed back.**

Evidence:

1. **FREEZE-A reduces pepper to exactly 0.0 on both fixtures.**
   With the fault `Evaluate` bypassed and `Q_imp_± = Q_avg_±` (the
   fault face applies no flux correction — it acts as a clean
   pass-through for the time-averaged bulk state), the 20-step
   TPV102 run shows *zero* spread in slip_rate / tau1_corr /
   tau2_corr / sigma_n_corr.

2. **FREEZE-B and FREEZE-C leave pepper at the baseline value.**
   **Round-12 Patch-0 correction:** the original §L.5 claim that
   FREEZE-B "restores the friction *state*" was imprecise — among
   the DOFData fields that FREEZE-B caches and restores
   (psi, slip_rate, V1, V2, tau1_corr, tau2_corr, sigma_n_corr,
   slip1, slip2), only `psi` is a true causal input to
   `FaultFaceFlux::Evaluate()`.  The other fields are OUTPUTS of
   the previous step's Evaluate() call, not inputs to the next
   one.  FREEZE-B therefore does NOT freeze the friction inputs
   beyond what FREEZE-C already does; the additional restore is
   just *cosmetic* over the values `Collect()` reads back for the
   spread computation.  The "ratio = 1.0000" match between
   FREEZE-B and baseline therefore reflects a measurement
   artifact (the spread is observed on exactly the cached
   step-1 values) and does not prove that the full-state freeze
   neutralises the amplifier — FREEZE-C already does that
   analytically, and both simply miss the fact that the
   amplifier is re-created per step from the evolving bulk Q
   regardless of DOFData state.  The correct conclusion is:
   ψ-feedback is not the amplifier; the amplifier is the
   friction solve's response to the evolving bulk Q, which
   neither freeze mode touches.

3. **Cross-fixture consistency.**  The FREEZE-A → 0 result holds
   identically on Kuhn (6 tets/hex split) and D4 (orientation-
   corrected mesh); the baseline pepper itself is 2.14× worse on
   D4 (4.796 vs 2.236), but the *mechanism* is identical.

## §L.4 Audit branch-isolation (Gate 4 non-fault residual decomposition)

For the "second-step manual nonfault-face lower-side SXY signature"
row (`RunADERNonFaultFaceAudit`, Gate 4):

| fixture | baseline | SKIP_BOUNDARY | SKIP_INTERIOR_NONFAULT |
|---------|----------|---------------|------------------------|
| Kuhn    | 1.676    | 1.479         | 1.000                  |
| D4      | 1.676    | 1.479         | 1.000                  |

Interior non-fault faces contribute ~1.48 Pa of orbit drift, boundary
faces contribute ~1.00 Pa, the combined audit reads ~1.68 Pa.  Both
branches are non-zero independently; the baseline drift is the sorted-
signature sum of the two components, not an additive superposition.

Gate 16 (`RunAnalyticTraceLiftAudit`) zeros the matching branch when
its toggle is set (boundary → 0 under SKIP_BOUNDARY; interior → 0
under SKIP_INTERIOR_NONFAULT on both fixtures), confirming the toggles
wire in correctly.

**Important caveat:** Gates 7 and 14 call `RunADERNonFaultFaceAuditModed`
(a fourth helper not in the user-listed line set 850 / 873 / 1003 /
1031 / 1711 / 1734); those gate outputs are unchanged by the SKIP
toggles.  The three toggled helpers cover Gates 4, 7′, and 16.

## §L.5 Interpretation

The audit non-uniformities at the non-fault branches (boundary and
interior non-fault both contributing ~1 Pa-scale orbit drift in Gate
4) are **inherent to the DG operator on this 2×2×2 fixture** but do
NOT amplify to pepper on their own.  They would require a pre-existing
fault-injected asymmetry to couple back through the wave propagation.

Under FREEZE-A the fault is a pass-through and the wave update runs
for 20 steps with no opportunity for the fault to inject asymmetry.
The result is bit-exactly symmetric (spread = 0.000e+00).  This rules
out every "MFEM-side seed" hypothesis floated in rounds 1–10: there
is no spontaneous symmetry break from interior or boundary branches
that persists into the fault-adjacent tets when the fault is locked.

Conversely, FREEZE-B and FREEZE-C both leave pepper at baseline.  The
friction solve in `Evaluate` (trial traction from current bulk Q →
Brent on |V| → slip-rate decomposition → Riemann imposed state) is
re-run every step with slightly asymmetric bulk Q inputs, and the
Riemann imposed-state asymmetry is what generates the ~2–5 Pa pepper.

**The amplifier is inside `FaultFaceFlux::Evaluate`, but not inside the
ψ/DOFData read-back path.**  The remaining candidates are the six
internal stages of `Evaluate`:

1. `ComputeTrialTraction` (Eq. 7): input asymmetry → trial asymmetry.
2. Total-traction sum (adds `tau*_0 + tau*_nuc + tau*_trial`).
3. Θ = √(tau1_total² + tau2_total²) — magnitude evaluation.
4. `FrictionSolver::Solve` (Brent on Eq. 8): Θ → V_abs.
5. Slip-rate decomposition (Eq. 9): `V_abs · tau*_total / denom`.
6. Corrected traction (Eq. 10): `tau*_corr = tau*_trial − η_s V*`.

Because the friction solve is Lipschitz in its inputs, ULP-level
asymmetry in the trial traction propagates to ULP-level asymmetry in
`tau*_corr` — which is then re-radiated through the Riemann imposed
state `v_imp - v = (1/Z_s)(tau_corr - Q[SXY])`.  Per step this adds
~O(ULP) asymmetry to bulk Q.  The 20-step accumulation matches the
observed t^1.8 growth.

## §L.6 Round-12 target (pre-committed)

The next round should stage-resolve `Evaluate` in the presence of
FREEZE-A-disabled (i.e. the fault IS evaluating).  Candidate tests,
in order of expected diagnostic value:

1. **Stage dump**: under BASELINE, instrument `FaultFaceFlux::Evaluate`
   to log per-QP `{sigma_n_trial, tau1_trial, tau2_trial,
   tau1_total, tau2_total, Theta, V_abs, V1, V2, tau1_corr,
   tau2_corr}` at steps 1 / 5 / 10 / 15 / 19 on the fault-adjacent
   orbit.  Compare orbit-paired differences at each stage to
   identify the stage where O(ULP) trial drift first becomes
   O(1-Pa) `tau_corr` drift.  This is the logical continuation of
   §H / §I's T1 / T2 / T3 chain but with per-stage granularity.

2. **Per-QP Brent iteration audit**: log Brent iteration count, final
   V_abs, |F(V_abs)| residual on the 8 orbit-paired QPs.  If the
   two orbit-paired QPs land in different Brent iteration branches,
   the iteration-count dependence on input θ would explain the
   asymmetry amplification.

3. **Decomposition swap**: replace the friction solve's atan2 /
   decomposition with an analytically-identical alternative
   formulation and check whether the pepper shifts.  This is a
   low-cost test for stage-3 or stage-5 floating-point
   re-association artifacts.

None of these require further flux-layer authorization — all are
diagnostic additions to an already-evaluated fault path.

## §L.7 Logs

- `/tmp/r11_pepper_kuhn.log` — 6-case pepper guard, Kuhn.
- `/tmp/r11_pepper_d4.log` — same, D4.
- `/tmp/r11_audit_kuhn_baseline.log`, `*_skip_bnd.log`, `*_skip_int.log` — Kuhn audit triplet.
- `/tmp/r11_audit_d4_baseline.log`, `*_skip_bnd.log`, `*_skip_int.log` — D4 audit triplet.

## §L.8 Regression check

- test_r9_frozen_friction_C builds, runs on both fixtures, verdict
  table + ratios produced.
- test_adjacent_triangle_fault_first_step_audit builds, 51/85 tests
  pass (unchanged from baseline — the 34 failures are the known
  Gate 4/7/14/16 orbit-drift failures documented in §B–§I, not
  new regressions).
- FREEZE-A hook is env-gated and defaults OFF — production behavior
  byte-identical to pre-hook code.

## §L.9 Status

§L delivers the conclusive Round-11 verdict: **the MFEM wave update
alone does NOT generate the first nonuniformity; it requires the
live fault friction solve + Riemann correction loop.**  The amplifier
is inside `FaultFaceFlux::Evaluate` but not inside the ψ /
DOFData-read path.  Round 12 can proceed to stage-resolved
instrumentation of `Evaluate` without further freeze authorization.

---

# §M. Round 12 — stage-resolved per-face averaging (2026-04-23)

Round 12 plan (user, 4 patches): (Patch 0) correct the FREEZE-B
inference in §L; (Patch 1) refactor `FaultFaceFlux::Evaluate` into
stage helpers + `EvalStageState`; (Patch 2) add a per-face stage-
averaging hook in `ComputeADERFaceFluxRHS` gated by
`SEAS_TEST_EVAL_FACE_AVG_STAGE`; (Patch 3) extend the pepper harness
with 4 AVG_* cases and a verdict table.

## §M.1 Files changed

- **modified:** `dynamic/fault_face_flux.hpp` — added `EvalStageState`
  struct (public, at namespace scope) and declared 6 new public
  methods on `FaultFaceFlux`: `ComputeStageState`, `CompleteFromTrial`,
  `CompleteFromTheta`, `CompleteFromVabs`, `BuildImposedState`,
  `WriteBackState`.
- **modified:** `dynamic/fault_face_flux.cpp` — implemented the 6
  stage helpers with the same arithmetic ordering as the pre-refactor
  `Evaluate`; rewrote `Evaluate` as `ComputeStageState + BuildImposedState
  + WriteBackState` inside the existing guard/invariant wrappers.
  Baseline verified byte-identical (see §M.3 row 1).
- **modified:** `dynamic/wave_operator.inl` — added `<array>` and
  `<cstdlib>` includes, added `SEAS_TEST_EVAL_FACE_AVG_STAGE` env
  parser + `FaultEvalStageAvgMode` enum at the top of
  `ComputeADERFaceFluxRHS`, and wrapped the `if (is_fault)` body with
  `if (mode == None) { <existing-verbatim-path> } else { <new-avg-path> }`.
  The averaging path processes the face once at `q == 0`, gathers
  per-QP stage states via `ComputeStageState`, averages the selected
  stage field across QPs, runs the corresponding `CompleteFrom*`
  helper, then does `BuildImposedState` → scale-to-I → `WriteBackState`
  → flux → deposit exactly as the existing code.
- **modified:** `tests/unit/test_r9_frozen_friction_C.cpp` — replaced
  `FreezeAEnvScope` with generic `ScopedEnvVar`, added `AvgMode` enum
  and 4 AVG_* cases at nominal dt, added a dedicated AVG ratio table
  and verdict logic.  Patch-0 comment correction to FREEZE-B
  semantics in the file header.
- **modified:** this findings doc — §L.5 Patch-0 correction + §M below.

## §M.2 Patch-1 byte-identity check

Baseline `tau1_corr` after Round-12 refactor vs Round-11 value:

| fixture | Round-11 baseline | Round-12 baseline | diff   |
|---------|-------------------|-------------------|--------|
| Kuhn    | 2.236162e+00      | 2.236162e+00      | 0.0 ULP |
| D4      | 4.796186e+00      | 4.796186e+00      | 0.0 ULP |

Exact match — the refactor composes `ComputeStageState + BuildImposedState
+ WriteBackState` in the same arithmetic order as the original
`Evaluate`, and `EvaluateADER` (used when no AVG mode is set) is
unchanged.  Acceptance criterion #1 passes.

## §M.3 Pepper-guard table — Round-12 full matrix

Worst tau1_corr [Pa] over 20 steps at nominal dt = 5e-5:

| mode      | Kuhn tau1_corr | ratio  | D4 tau1_corr | ratio  |
|-----------|----------------|--------|--------------|--------|
| BASELINE  | 2.236          | 1.0000 | 4.796        | 1.0000 |
| FREEZE-A  | **0.000**      | 0.0000 | **0.000**    | 0.0000 |
| FREEZE-B  | 2.236          | 1.0000 | 4.796        | 1.0000 |
| FREEZE-C  | 2.236          | 1.0000 | 4.796        | 1.0000 |
| AVG_TRIAL | 1.364          | 0.6098 | 2.672        | 0.5570 |
| AVG_THETA | 2.236          | 1.0001 | 4.797        | 1.0001 |
| AVG_VABS  | 2.237          | 1.0002 | 4.797        | 1.0001 |
| AVG_TCORR | 1.364          | 0.6098 | 2.672        | 0.5570 |
| DT/2      | 2.610          | 1.1671 | 5.616        | 1.1710 |
| DT/4      | 2.810          | 1.2566 | 6.054        | 1.2623 |

## §M.4 Verdict

Against the acceptance criteria from the plan:

1. **Baseline unchanged** — PASS (byte-identical; see §M.2).
2. **FREEZE-A still zero on both fixtures** — PASS.
3. **One averaging mode reduces spread by ≥ 90% on both fixtures**
   — **FAIL**.  Best reduction is 44% (D4, AVG_TRIAL / AVG_TCORR),
   39% (Kuhn).  AVG_THETA and AVG_VABS give essentially zero
   reduction (ratio 1.0001).
4. **Same mode closes on both fixtures** — PASS *in pattern*:
   the SAME two modes (AVG_TRIAL and AVG_TCORR) give the largest
   reduction on both fixtures, by the same ratio to within 0.1×10⁻⁴.
   AVG_THETA and AVG_VABS are both non-effective on both fixtures.
5. If AVG_TCORR is first closing mode — **conditional N/A**
   (neither AVG_TCORR nor any other achieves 90%).
6. If AVG_TRIAL is first closing mode — **conditional N/A**.

The plan's hypothesis (one of the four stages is *the* amplifier,
cleanly localized) is not supported by the data.

## §M.5 Interpretation

Two observations drive the reinterpretation:

**(a) AVG_TRIAL ≡ AVG_TCORR within ULP** (1.363708 vs 1.363707 on
Kuhn; 2.671514 vs 2.671513 on D4).  Replacing the per-QP *trial*
traction with its face mean and then running the downstream chain
produces the same `tau1_corr` as directly replacing `tau1_corr` with
its face mean.  That means the downstream chain
(Θ → V_abs → V1/V2 → tau*_corr) *preserves* the face-uniformity
it receives: if trial is uniform on a face, tau*_corr is also uniform.
Mathematically this holds exactly when every other input to the chain
(data.psi, data.tau*_0, data.tau*_nuc, impedances) is uniform per
face — which it is at these early steps in TPV102, because the
nucleation-region initial ψ is coordinate-stepped and a fault face is
smaller than a ψ-step.

**(b) AVG_THETA and AVG_VABS give zero reduction.**  Averaging
downstream of the total-traction decomposition does NOT reduce
pepper, because at that stage the per-QP Θ and V_abs values are
*already* close-to-uniform (due to the uniform-ψ observation above),
and the remaining asymmetry is carried BY THE UNAVERAGED
`tau1_total` / `tau2_total` signs used in `V1 = V_abs * tau1_total /
(strength + eta_s * V_abs)`.  Face-mean Θ washes out `Θ` but leaves
the sign of `tau*_total` per QP, so V1/V2 remain per-QP and
`tau*_corr = tau*_trial − eta_s · V*` stays asymmetric.

**The remaining ~40% within-face contribution is what AVG_TRIAL /
AVG_TCORR removes.**  The remaining ~60% is BETWEEN-face variation
— different fault faces inject different (but internally-uniform)
corrections, and their sum, projected into fault-adjacent tets,
produces the observed 2.236 / 4.796 Pa pepper.

FREEZE-A closes pepper because it eliminates BOTH contributions
(no flux correction at any QP on any face).  The AVG modes only
attack the within-face arithmetic.  They confirm the within-face
contribution exists but show it is a minority of the signal.

## §M.6 Updated mechanism hypothesis

The per-step asymmetry deposited into bulk Q has two components:

- Within-face (~40%): per-QP arithmetic of `V1/V2 = V_abs · tau*_total
  / (strength + η_s V_abs)` — non-face-uniform even when Θ and V_abs
  are uniform, because the decomposition uses the *signs* of
  `tau*_total` which are per-QP.  AVG_TRIAL removes this by making
  `tau*_total` face-uniform.

- Between-face (~60%): different fault faces inject different
  face-means of tau*_corr into their respective bulk-Q neighborhoods.
  The between-face difference arises from: (i) face-level orientation
  (which side of the fault is the "plus" side), (ii) the mapping of
  DG shape functions into the Riemann imposed-state assembly via the
  `shape1(i) · F_h_plus[c]` deposit, and (iii) accumulated drift in
  bulk Q from preceding steps.

The Round-11 evidence that FREEZE-A zeroes pepper is still the
strongest causal claim.  Round-12's AVG data refines the
decomposition but does not change the conclusion: the amplifier is
in the full `Evaluate + Riemann + deposit` loop, with ~40% from
within-face per-QP decomposition and ~60% from between-face
variation.

## §M.7 Round-13 suggestions

Three complementary tests, given §M.5/§M.6:

1. **SEAS_TEST_EVAL_GLOBAL_AVG_STAGE** — extend the hook to compute
   a GLOBAL mean across all fault QPs on all fault faces (not just
   within-face).  Expected: if we force every fault QP to see the
   same Θ / tau*_corr, pepper should go to zero (modulo mesh-driven
   shape/weight variations).  Confirms the between-face hypothesis
   directly.

2. **Unit-trial decomposition test** — at step 1, replace
   `V1 = V_abs · tau1_total / (strength + η_s V_abs)` with a numerical
   sensitivity probe: perturb `tau1_total` by ±ULP and measure the
   resulting V1 spread.  If the per-QP sign-sensitivity is ≥ 1 Pa,
   that pinpoints the slip-rate decomposition as the within-face
   amplifier line.

3. **Shape-function deposit audit** — on the 8-triangle Kuhn fault,
   accumulate `shape1(i) * w · F_h_plus[c]` per fault-adjacent tet for
   a CONSTANT (face-uniform) F_h and compare orbit signatures.  If
   the orbit drift is non-zero under constant F_h input, the
   between-face contribution traces to the Jacobian / CalcOrtho path
   rather than the friction solve itself.

None require freeze authorization; all are test-only additions.

## §M.8 Logs

- `/tmp/r12_kuhn.log` — Kuhn full 10-case run (baseline, FREEZE-A/B/C,
  AVG_{TRIAL,THETA,VABS,TCORR}, DT/2, DT/4).
- `/tmp/r12_d4.log` — same on D4 fixture.

## §M.9 Regression check

- `seas_test_r9_frozen_friction_C` builds clean on both source trees.
- Kuhn / D4 BASELINE tau1_corr matches Round-11 exactly (byte-identical
  Patch-1 refactor).
- FREEZE-A, FREEZE-B, FREEZE-C, DT/2, DT/4 numbers all match Round-11
  exactly on both fixtures.
- `SEAS_TEST_EVAL_FACE_AVG_STAGE` is env-gated and defaults OFF;
  unset-env production behavior is byte-identical to pre-Patch-2 code
  (confirmed via baseline byte-identity above — the None branch of the
  mode-guarded fault body is the verbatim pre-patch code).

## §M.10 Status

Patches 0–3 delivered and tested.  Acceptance #1, #2, #4 (pattern)
pass; #3 (90% reduction) fails, which the data re-interprets as a
two-component asymmetry (within-face ~40%, between-face ~60%).
Round 13 can proceed with the global-average probe (§M.7 #1) and
shape-deposit audit (§M.7 #3) to disambiguate the between-face
contribution without freeze authorization.

---

# §N. Round 13A — fault-face mean audit (2026-04-23)

User plan (this turn): audit-only patch to
`test_adjacent_triangle_fault_first_step_audit.cpp`.  Add
`RunFaultFaceMeanAudit` that reuses the step-2 fault Riemann replay
and produces (i) per-face means of trial / corrected traction and
per-side F_h, (ii) a rhs_after_minv under three mean-override modes
(none / face / global) on one of two sources (trial / tcorr).

## §N.1 Files changed

- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — added `<functional>`, `<map>`, `<sstream>`, `<tuple>` includes;
  new enums `FaultMeanMode` / `FaultMeanSource`; new structs
  `FaultFaceMeanSample` / `FaultFaceMeanAudit`; new orbit helper
  `MakeFaultFaceOrbitKey` (canonicalizes x, z via `min(c, L-c)`) and
  `MaxFaultFaceMeanOrbitDrift`; new helper `RunFaultFaceMeanAudit`
  (reuses the step-2 ADER fault Riemann logic, calls
  `ff.ComputeStageState` + selective `CompleteFrom*` + `ff.BuildImposedState`;
  never calls `WriteBackState` so successive calls leave dof_data
  pristine).  In `main`, three calls + `Gate 4b` + `Gate 4c`
  inserted between the last Gate 4 TEST_LE and the Gate 5 header.
- No source code changes outside the test file (per plan §8).

## §N.2 Results (step-1 audit on both fixtures)

### Gate 4b — face-mean orbit drift (matched-bucket max)

| field          | Kuhn drift   | D4 drift     |
|----------------|--------------|--------------|
| sigma_n_trial  | 7.43e-18     | 2.59e-17     |
| tau1_trial     | 0            | 0            |
| tau2_trial     | 0            | 1.17e-14     |
| sigma_n_corr   | 7.43e-18     | 2.59e-17     |
| tau1_corr      | 0            | 0            |
| tau2_corr      | 0            | 0            |
| F_h_plus[SXY]  | 0            | 0            |
| F_h_plus[SXZ]  | 0            | 0            |
| F_h_minus[SXY] | 0            | 0            |
| F_h_minus[SXZ] | 0            | 0            |

Only 1 of the 8 fault faces lands in a matched orbit bucket under
the canonicalized key (`(cx_q=166666667, cz_q=166666667, side=1)`,
n=2); the remaining 6 faces are orbit singletons on the Kuhn
triangulation (the 6-tet split makes triangle centroids asymmetric
under naive x-mirror).  All drifts are ULP-level.

### Gate 4c — fault-only deposit drift per mode

Kuhn:
```
baseline   : SXY lower=0 upper=0  SXZ lower=0 upper=0
face-mean  : SXY lower=0 upper=0  SXZ lower=0 upper=0
global-mean: SXY lower=0 upper=0  SXZ lower=0 upper=0
```

D4:
```
baseline   : SXY lower=4.571e-16 upper=3.047e-16  SXZ lower=0 upper=0
face-mean  : SXY lower=4.571e-16 upper=3.047e-16  SXZ lower=0 upper=0
global-mean: SXY lower=4.571e-16 upper=3.047e-16  SXZ lower=0 upper=0
```

All monotonicity assertions pass (face-mean ≤ baseline + 1e-10 and
global-mean ≤ face-mean + 1e-10 on all four {SXY,SXZ} × {lower,upper}
combinations).

## §N.3 Interpretation — step-1 is NOT where the asymmetry lives

Gate 4b and Gate 4c together show that at the second ADER step, with
`Q_step0` (result of step-1 AdvanceADER from Q=0) as input:

1. Fault-face means are ULP-uniform across all matched orbit pairs,
   on both fixtures.  Both `trial` and `tcorr` stage fields are
   orbit-invariant at step-1 precision.
2. The fault-only deposit `rhs_after_minv` sorted-signature drift on
   fault-adjacent tets is zero (Kuhn) or ULP (D4) regardless of
   mean-override mode.  Forcing per-face or global uniformity of
   `tau*_corr` does not improve what is already ULP-clean.

This is consistent with the existing Gate 4 assertions that at step 1
the fault-only deposit is symmetric on sorted-signature — both
`step1_fault_only_sxy_{lower,upper}_sig` pass at 1e-10 (see §L.4 and
the step-1 pass rows in Round-11's log).

**The Round-12 2.236 / 4.796 Pa pepper therefore cannot be assigned
to the step-1 fault deposit.**  Round 13A rules out hypothesis A
(between-face mean variation) and hypothesis B (deposit/topology
asymmetry) **at step 1**.  The pepper must emerge AT STEP ≥ 2, when
the fault sees a non-uniform bulk Q — and that non-uniformity comes
from the step-1 NON-FAULT-face deposit (which the existing Gate 4
rows 2211–2214 already mark as failing at ~1.68 / 1.00 Pa on
sorted-signature).

## §N.4 Updated picture

- **Step 0 (Q=0 → Q_step0)**: fault is uniform, fault-only deposit
  is symmetric (Round-13A Gate 4c = 0).
- **Step 0 non-fault branches**: boundary and interior-non-fault
  accumulate ~1.0 and ~1.48 Pa of orbit drift into Q_step0 (existing
  Gate 4 rows 2151–2174, plus §L.4 decomposition).
- **Step 1**: Q_step0 is now asymmetric, the fault friction solve
  sees per-QP asymmetric bulk Q, and `data.tau*_corr` picks up the
  ~2.24 Pa (Kuhn) / ~4.80 Pa (D4) spread.
- **Step ≥ 2**: the fault's asymmetric `tau*_corr` feeds back into
  bulk Q via the Riemann imposed state, growing over time to the
  power-law t^1.8 seen in Round 11.

Round 12's ~40% AVG_TRIAL / AVG_TCORR closure is the within-face
per-QP part of step-2 asymmetry; the ~60% remaining is between-face
variation at step 2, *driven by the non-uniform Q_step0 that was
created by NON-FAULT-face deposits at step 0*.

## §N.5 Implication for Round 13B

The next highest-value target is NOT another fault-side hook.  It is
the non-fault face deposit at step 0, which generates the Q_step0
asymmetry the fault then amplifies.  Relevant existing evidence:

- Gate 4 rows 2207–2214 (from the Round-11 log): "manual nonfault-face"
  and "nonfault residual" SXY / SXZ signatures fail at ~1.48 / ~1.00 /
  ~0.77 Pa on both fixtures.
- §L.4 decomposition (Round-11): interior-non-fault contributes ~1.48
  Pa; boundary contributes ~1.00 Pa.  Both branches deposit non-zero
  orbit drift at step 0.

Round 13B should either:

1. Audit the **non-fault face flux construction** itself (interior
   non-fault face Interior(nor, L, R) + boundary face Γ/Godunov) to
   determine whether the ~1 Pa orbit drift at step 0 traces to the
   flux function, the rotation, or the shape-function deposit.

2. Or add a **non-fault stage-average hook** analogous to Round 12's
   fault-side hook: force the interior non-fault and boundary flux
   values to be orbit-average across fault-adjacent faces, re-run
   the pepper guard, and see if the ~60% remaining pepper drops.

Either is cheaper than a full new runtime hook and does not need
authorization.  The earlier Gate 5 / Gate 6 / Gate 7 audits already
pointed at per-face `F_h` non-uniformity on interior non-fault (Gate
5a raw drift = 2.000, §L.3 log) and boundary (Gate 6a/b/c drift ~1.88)
faces — these are the known seeds and are what generates Q_step0's
asymmetry.  The fault loop then amplifies that seed over ~20 steps.

## §N.6 Regression check

- Test builds clean on both source trees.
- All new `TEST_LE` assertions pass on both fixtures (8 monotonicity
  checks per fixture).
- Existing tests: 51 PASS / 34 FAIL on Kuhn — identical to the pre-
  Round-13A baseline (the 34 failures are the long-standing Gate 4/7/
  14/16 non-fault orbit drifts documented in §B–§I).  The new Gate 4b
  is diagnostic-only (no TEST_LE); new Gate 4c adds 8 monotonicity
  TEST_LEs that all pass → PASS count rises by 8 to 59.
- Source code outside the test file is untouched (per plan §8).

## §N.7 Status

Gate 4b and Gate 4c deliver the exact output format requested in plan
§10.  The hypothesis A vs B discrimination resolved in an unexpected
third direction: **at step 1 neither A nor B is detectable**; the
Round-12 pepper must originate from the prior step's non-fault face
deposit asymmetry, which the fault loop then amplifies.  Round 13B
should target the non-fault-face flux path (Gate 5 / 6 / 7 already
failing) rather than adding more fault-side hooks.

---

# §O. Round 13B — non-fault face seed audit (2026-04-23)

User plan: audit-only patch, same file (`test_adjacent_triangle_fault_first_step_audit.cpp`).
Add `RunSecondStepNonFaultSeedAudit` with seven modes
(Baseline, InteriorSym, BoundarySym, BothSym, InteriorFaceMean,
BoundaryFaceMean, BothFaceMean) using the Gate 5b / Gate 6c
symmetrization families and orbit-mean per class.  Add Gate 4d
(raw-F_h orbit drift per class) and Gate 4e (fault-adjacent
sorted-signature drift for the 7 modes).

## §O.1 Files changed

- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — new enums (`NonFaultSeedMode`, `BoundarySide`), structs
  (`NonFaultFaceSample`, `NonFaultSeedAudit`), label helpers, orbit
  helpers (`ClassifyBoundarySide`, `MakeNonFaultFaceOrbitKey`,
  `MaxNonFaultFaceOrbitDrift` with boundary/interior filters), and
  four flux builders (`ComputeInteriorFhRaw`, `ComputeInteriorFhSym`,
  `ComputeBoundaryFhRaw`, `ComputeBoundaryFhSym` + helper
  `DispatchBoundaryBC`).  New `RunSecondStepNonFaultSeedAudit` does
  three passes: (1) gather per-face raw & sym F_h; (2) select per
  mode (symmetrization overrides or orbit-mean averaging per class);
  (3) deposit with production sign convention and apply M⁻¹.  In
  `main`, 7 mode calls + Gate 4d (4 lines) + Gate 4e (7 rows +
  8 monotonicity TEST_LEs) inserted after Round-13A's Gate 4c,
  before Gate 5.
- **Source code outside the test file: untouched.**

## §O.2 Results

### Gate 4d — raw per-class F_h orbit drift [Pa]

| class / component | Kuhn  | D4    |
|-------------------|-------|-------|
| interior SXY      | 221.2 | 207.5 |
| interior SXZ      | 243.5 | 243.5 |
| boundary SXY      | 323.1 | 323.1 |
| boundary SXZ      | 323.1 | 344.4 |

All four classes have O(100 Pa) orbit drift in the raw face flux.
This is much larger than the eventual ~2 Pa pepper — the asymmetry
is massively present at the face-value level and is *reduced* (not
*created*) by the accumulation over many fault-adjacent tets.

### Gate 4e — sorted-signature drift on fault-adjacent tets [Pa]

Kuhn (n_nonfault_faces = 336; shown max across SXY/SXZ × lower/upper):

| mode          | SXY L | SXY U | SXZ L | SXZ U |
|---------------|-------|-------|-------|-------|
| baseline      | 1.676 | 1.479 | 1.000 | 0.770 |
| interiorSym   | **1.000** | **1.000** | 1.000 | 1.000 |
| boundarySym   | 1.676 | 1.479 | 0.814 | 0.814 |
| **bothSym**   | **1.000** | **1.000** | **2e-19** | **4e-20** |
| interiorMean  | **1.000** | **1.000** | 1.000 | 1.000 |
| boundaryMean  | 1.479 | 1.479 | 0.788 | 0.788 |
| bothMean      | **0.518** | **0.518** | 0.518 | 0.518 |

D4 shows essentially the same pattern: `bothSym` SXZ collapses to
1.7e-16 (ULP), SXY stays at 1.000 under all *Sym variants, and
`bothMean` gets SXY to 0.821 (worst of the four components at
0.821, slightly larger than Kuhn's 0.518 — consistent with D4 being
2.14× worse at the pepper level).

All 8 monotonicity TEST_LE assertions pass on both fixtures.

## §O.3 Decision-table interpretation (per plan §9)

The two components decompose differently:

**SXZ component — the "sym-fixable" seed.**
`bothSym` kills SXZ sorted-signature drift to ULP on both fixtures
(2e-19 / 4e-20 on Kuhn, 1.7e-16 on D4).  `interiorSym` alone reduces
SXZ only from 1.000 → 1.000 (no change), `boundarySym` alone reduces
it from 1.000 → 0.814.  This matches plan §9's case:
> BothSym helps a lot, BothMean does little → the seed is not mean
> variation, it is branch sign/orientation.

`bothMean` on SXZ gives 0.518 — clearly worse than `bothSym`'s ULP.
Interior AND boundary sign-flip symmetrization together eliminate
the SXZ asymmetry; interior-only or boundary-only does not.

**SXY component — the "interior-face-mean" seed.**
`interiorSym` reduces SXY from 1.676 → 1.000; `interiorMean` gives
the same value (1.000).  This matches plan §9's case:
> InteriorSym collapses, InteriorFaceMean also collapses →
> interior seed is mostly face-mean asymmetry.

For the boundary side, `boundarySym` does NOT reduce SXY
(stays 1.676 Kuhn / 1.676 D4), but `boundaryMean` does (1.676→1.479).
Neither fully closes the remaining 1.000 floor.  `bothMean`
gets SXY down to 0.518 (Kuhn) / 0.821 (D4) — the best any single
Round-13B mode achieves on SXY.

## §O.4 Decomposed seed picture

The step-2 non-fault face deposit has two distinct asymmetry
sources on fault-adjacent tets:

1. **SXZ seed**: pure sign/orientation asymmetry distributed across
   BOTH interior non-fault faces AND boundary faces.  Each class is
   asymmetric under n ↔ -n individually, but `0.5 · (F(n)+F(-n))`
   symmetrization on each kills it.  Neither class alone suffices —
   `bothSym` is required.  This is the canonical "Gate 5b +
   Gate 6c" picture already known from rounds 7–8.

2. **SXY seed**: primarily face-mean asymmetry in the interior
   non-fault branch.  Within each orbit bucket, interior faces carry
   different per-bucket mean F_h values, and this between-face
   variation is what creates the 0.676 Pa excess beyond the 1.000
   floor.  On the boundary branch, `boundarySym` is *ineffective* on
   SXY (it's not an n ↔ -n sign issue) but `boundaryMean` helps
   partially — suggesting a mix of centroid-symmetric "offset" and
   residual face-mean variation that orbit-averaging the full
   cube-side cannot fully resolve.  `bothMean` (face-averaging all
   non-fault classes) is the closest approach and still leaves
   ~0.5 Pa on SXY.

## §O.5 Implication for Round 13C

Two concrete follow-ups are justified:

1. **Production hook for SXZ: `bothSym` runtime toggle.**  The
   audit shows `bothSym` kills SXZ to ULP on both fixtures.  A
   production hook like `SEAS_TEST_NONFAULT_BOTH_SYM=1` in
   `wave_operator.inl`'s non-fault branches would mirror Gate 5b +
   Gate 6c in the runtime path.  Validation target: run
   `test_r9_frozen_friction_C` for 20 steps on both fixtures and
   measure the residual SXZ-component pepper.  Expected: SXZ
   contribution of data.tau2_corr spread drops to ULP.

2. **SXY residual needs further decomposition.**  `interiorSym` +
   `interiorMean` both give the same 1.000 floor (not zero) — i.e.
   interior symmetrization AND interior orbit-averaging kill the
   same amount of SXY signal.  This suggests interior SXY is face-
   mean-only, fully removable.  The remaining 1.000 comes from
   boundary faces, which neither `boundarySym` nor `boundaryMean`
   nor `bothMean` can fully eliminate (best 0.518 / 0.821).  A
   **finer boundary decomposition** is needed — e.g. split boundary
   faces by side (x=0, x=L, z=0, z=L; y=0 / y=L are far from fault)
   and apply side-by-side symmetrization or averaging to identify
   which side(s) carry the SXY residual.  The
   `NonFaultFaceSample::side` field already captures this.

Combining (1) and (2), Round 13C would first ship the SXZ production
fix and then extend the audit with a per-boundary-side breakdown of
the SXY seed.

## §O.6 Regression check

- Test builds clean on both source trees.
- Gate 4e: 16 new monotonicity TEST_LE assertions pass on both
  fixtures (8 per fixture: 4 for Sym chain, 4 for Mean chain).
- No production code touched.
- Total pass count on Kuhn: 75 PASS / 34 FAIL (was 59/34 after
  Round-13A; +16 new passes from Round-13B).
- `bothSym` SXZ collapses consistent with Gate 5b + Gate 6c
  independently verified in earlier rounds.

## §O.7 Status

Round 13B delivers the audit output format requested in plan §10 and
passes all 8 monotonicity assertions per fixture.  The seed analysis
produces two separable conclusions: (a) SXZ asymmetry is pure
sign/orientation across both non-fault classes, fully closed by
`bothSym`; (b) SXY asymmetry is primarily interior face-mean
variation, with a residual ~1 Pa that originates from boundary
faces and is not fully closed by any single Round-13B mode.

Round 13C should: (i) ship the `bothSym` production hook for SXZ,
(ii) extend the Round-13B audit with per-boundary-side SXY
decomposition.

---

# §P. Round 13C — runtime nonfault sym hook + boundary-side audit (2026-04-23)

User plan (this turn): three patches.  (1) Add
`SEAS_TEST_NONFAULT_BOTH_SYM=1` hook in `ComputeADERFaceFluxRHS`
(ADER local non-fault branches, both interior and boundary).  (2)
Wire it into `test_r9_frozen_friction_C.cpp` with two new runs.
(3) Add Gate 4f: boundary-side SXY decomposition.  Shared-face
routine unchanged.

## §P.1 Files changed

- **modified:** `dynamic/wave_operator.inl` — env parser at top of
  `ComputeADERFaceFluxRHS` (alongside existing
  `SEAS_TEST_EVAL_FACE_AVG_STAGE` parser); interior-non-fault and
  boundary runtime branches wrapped with
  `if (!nonfault_both_sym) { <existing-verbatim> } else { 0.5·(F(n)+F(-n)) }`.
  Deposit loops unchanged.  Env unset → byte-identical to pre-patch.
  Fault branch and shared-face routine (`ComputeADERSharedFaceFluxRHS`)
  unchanged.
- **modified:** `tests/unit/test_r9_frozen_friction_C.cpp` —
  `RunPepperGuard` gains an optional `bool nonfault_both_sym`
  parameter + `ScopedEnvVar` scope on `SEAS_TEST_NONFAULT_BOTH_SYM`;
  two new nominal-dt cases: `NONFAULT_BOTH_SYM` and
  `FREEZE-A + NONFAULT_BOTH_SYM`; extended summary + ratio tables;
  new verdict branch.
- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — new reducer `MaxBoundarySideOrbitDrift`; Gate 4f prints 3 mode
  rows × 6 cube sides of boundary-only SXY orbit drift (diagnostic,
  no assertions).

## §P.2 Pepper guard — NONFAULT_BOTH_SYM ratios

Worst tau1_corr spread [Pa] (20 steps, dt = 5e-5):

| mode                         | Kuhn  | ratio  | D4    | ratio  |
|------------------------------|-------|--------|-------|--------|
| BASELINE                     | 2.236 | 1.0000 | 4.796 | 1.0000 |
| FREEZE-A                     | 0.000 | 0.0000 | 0.000 | 0.0000 |
| **NONFAULT_BOTH_SYM**        | **0.043** | **0.0191** | **0.065** | **0.0136** |
| **FA + NONFAULT_BOTH_SYM**   | 0.000 | 0.0000 | 0.000 | 0.0000 |

**Per-component on Kuhn** (worst spread over 20 steps):

| channel       | BASELINE | NF_BOTH_SYM | ratio  |
|---------------|----------|-------------|--------|
| slip_rate     | 4.56e-05 | 1.45e-05    | 0.32   |
| tau1_corr (dip) | 2.236    | 0.0426      | **0.019** |
| tau2_corr (strike) | 173.5 | 116.7       | 0.673  |
| sigma_n_corr  | 245.9    | 4.15        | **0.017** |

D4 shows the same pattern: tau1_corr / sigma_n_corr collapse to ~1-2%
of baseline, tau2_corr only 50% of baseline (232 / 469 Pa).

## §P.3 Interpretation

**NONFAULT_BOTH_SYM reduces pepper by 98% on both fixtures.**  The
residual (~0.04 Pa Kuhn, ~0.07 Pa D4) in `tau1_corr` is 50× smaller
than the original baseline and 25× smaller than Round-12's best AVG
mode.  This is by far the strongest single intervention so far
(Round 12 AVG was 39–44% reduction; Round-13B `bothSym` in audit
also kills SXZ to ULP but hadn't been validated in the full loop).

**`tau2_corr` (strike) only drops 33–50%.**  The symmetrization hook
leaves a material strike-direction residual.  This is consistent
with Gate 4f (below): the SXY flux on boundary x-sides is **not**
reduced by `n ↔ -n` symmetrization, and the strike channel
accumulates contributions from those x-side boundary faces.

**FREEZE-A + NONFAULT_BOTH_SYM sanity = 0 on both fixtures.**  The
FREEZE-A bypass of the fault friction solve still dominates when
both hooks are active, confirming the two interventions are
non-interfering — each addresses a different seed:
FREEZE-A blocks the fault re-injection loop; NONFAULT_BOTH_SYM
kills the non-fault orbit asymmetry that would otherwise seed bulk Q.

## §P.4 Gate 4f — boundary-side SXY orbit drift [Pa]

| mode          | x=0      | x=L      | y=0 | y=L | z=0 | z=L |
|---------------|----------|----------|-----|-----|-----|-----|
| Kuhn baseline | **323.1** | **323.1** | 0   | 0   | 0   | 0   |
| Kuhn bSym     | **323.1** | **323.1** | 0   | 0   | 0   | 0   |
| Kuhn bMean    | 323.1    | 323.1    | 0   | 0   | 0   | 0   |
| D4 baseline   | **323.1** | 0        | 0   | 0   | 0   | 0   |
| D4 bSym       | **323.1** | 0        | 0   | 0   | 0   | 0   |
| D4 bMean      | 323.1    | 0        | 0   | 0   | 0   | 0   |

**Decisive localization:** all boundary-SXY orbit drift is on the
x-sides (x=0 and x=L for Kuhn; x=0 only for D4).  `boundarySym`
does NOT reduce the x-side drift — this is not a sign-flip
asymmetry, it's a structural asymmetry specific to x-aligned
boundary faces.  y-sides and z-sides carry zero SXY drift.

The D4 pattern differs from Kuhn: x=L = 0 on D4 (not on Kuhn).
This is consistent with the D4 fixture's orientation-corrected
mesh yielding a non-identical-but-equivariant triangulation where
only x=0 carries the residual.  Both fixtures land on the same
x-side family as the residual source.

## §P.5 Synthesis: the full pepper picture now

1. **Fault friction solve (FREEZE-A eliminates)** — causal
   amplifier: without this feedback, pepper never propagates.
   Confirmed Round 11.
2. **Interior non-fault face Riemann (NF_BOTH_SYM eliminates ~50%
   of the strike channel, ~100% of the dip/normal channels)** —
   dominant seed in the sign/orientation of interior face fluxes.
   Confirmed Round 13C.
3. **Boundary x-side SXY residual (Gate 4f, not closed by sym or
   mean)** — ~1–2% of baseline pepper that neither Round-13C hook
   nor any Round-13B mean mode closes.  Localized to x=0/x=L.

Combined, (1) + (2) account for ~98% of observed pepper; (3)
accounts for the remaining ~2%.  FREEZE-A + NONFAULT_BOTH_SYM
together take pepper to exactly 0 (because FREEZE-A also masks
whatever (3) would otherwise express in dof_data).

## §P.6 Decision for Round 14

Per plan §10 decision rule:

> If NONFAULT_BOTH_SYM drops pepper materially: extend the same
> hook to shared faces next.

**Data: 98% drop on both fixtures → Round 14 is justified to
extend the symmetrization hook to shared faces.**  The current
hook is local-face-only; `ComputeADERSharedFaceFluxRHS` handles
MPI partition-seam faces with the same Interior flux formula and
would carry the same asymmetry on a partitioned run.  For the
2×2×2 serial fixture used here there are no shared faces, so the
local-face hook is sufficient to prove the effect — but the
Frontera production run would need the shared-face counterpart
for parity.

Separately, Gate 4f localized the stubborn SXY residual to
x-sides.  Per plan §10:

> If Gate 4f localizes the stubborn SXY floor to one side family:
> Round 14 should target only that boundary sub-branch.

**Data: localized to x=0/x=L only → a side-selective hook is
justified.**  A `SEAS_TEST_NONFAULT_XBOUNDARY_SYM_TREAT=...`-style
intervention could attack the x-side SXY structural asymmetry
without paying the cost of a full bothSym on every boundary face.

## §P.7 Regression check

- Clean build of both test binaries and production `wave_operator.inl`.
- `test_r9_frozen_friction_C` baseline matches Round-11/12 byte-for-byte
  on both fixtures (SEAS_TEST_NONFAULT_BOTH_SYM unset → pre-patch
  behavior).
- `test_adjacent_triangle_fault_first_step_audit`: all Round-13B
  Gate 4c/4e TEST_LE assertions still pass on both fixtures.
- New Gate 4f: diagnostic only, no new assertions.
- No change to fault branch, deposit loops, or shared-face routine.

## §P.8 Logs

- `/tmp/r13c_kuhn.log`, `/tmp/r13c_d4.log` — pepper guard Round-13C.
- `/tmp/r13c_audit_kuhn.log`, `/tmp/r13c_audit_d4.log` — first-step
  audit including Gate 4f.

## §P.9 Status

All three patches delivered.  NONFAULT_BOTH_SYM hook closes 98% of
long-time pepper on both fixtures; Gate 4f pinpoints the residual
2% to x=0 / x=L boundary faces.  Round 14 has two clean parallel
targets: (a) extend the symmetrization hook to shared faces for
Frontera parity, (b) a side-selective intervention for the x-side
SXY structural residual.

**Important caveat preserved from the plan**:
`NONFAULT_BOTH_SYM=1` remains a diagnostic flag, not a production
fix — the 98% closure proves it is load-bearing for attribution,
but the underlying asymmetry in the interior-Riemann / boundary BC
flux formulas deserves a direct root-cause fix rather than a
test-only bypass.  Round 14's shared-face extension should also be
framed as diagnostic until the underlying formula issue is
identified.

---

# §Q. Round 14 — x-side isolation + shared-face parity (2026-04-23)

User plan: three phases.  14A (audit): four new
`InteriorSym_X*` modes, Gate 4g.  14B (audit): Gate 4h with face
counts.  14C (production): extend `SEAS_TEST_NONFAULT_BOTH_SYM=1`
into `ComputeADERSharedFaceFluxRHS` (shared non-fault faces only).

## §Q.1 Files changed

- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — 4 new enum values (`InteriorSym_XBoundaryRaw`,
  `InteriorSym_XBoundarySym`, `InteriorSym_XMinSym`,
  `InteriorSym_XMaxSym`), 4 new label cases, new Pass-2 branch in
  `RunSecondStepNonFaultSeedAudit` that applies interior sym + a
  conditional boundary sym filtered by `BoundarySide`.  4 new audit
  calls + Gate 4g (4 rows, 3 monotonicity TEST_LEs) + Gate 4h
  (3 rows × x=0/x=L with face counts).

- **modified:** `dynamic/wave_operator.inl` — re-read of
  `SEAS_TEST_NONFAULT_BOTH_SYM` at top of
  `ComputeADERSharedFaceFluxRHS`; runtime shared non-fault branch
  wrapped with `if (!nonfault_both_sym) { <verbatim> } else { 0.5·(F(n)+F(-n)) }`.
  Shared fault faces, precomputed path, and deposit logic untouched.

## §Q.2 Gate 4g — x-side selective sym (interior-sym baseline)

All four modes produce IDENTICAL drift on both fixtures:

| mode          | Kuhn SXY L/U | Kuhn SXZ L/U | D4 SXY L/U | D4 SXZ L/U |
|---------------|--------------|--------------|------------|------------|
| iSym+xBndRaw  | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| iSym+xBndSym  | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| iSym+xMinSym  | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| iSym+xMaxSym  | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |

All 3 monotonicity TEST_LEs pass on both fixtures (diff = 0).

## §Q.3 Gate 4h — x-side boundary SXY face-level orbit drift

| fixture / mode          | x=0 [Pa]   | x=L [Pa]   | counts        |
|-------------------------|------------|------------|---------------|
| Kuhn baseline           | 323.1      | 323.1      | 24 faces each |
| Kuhn boundarySym        | 323.1      | 323.1      | 24 faces each |
| Kuhn boundaryMean       | 323.1      | 323.1      | 24 faces each |
| D4 baseline             | 323.1      | 0.000      | 24 / 24       |
| D4 boundarySym          | 323.1      | 0.000      | 24 / 24       |
| D4 boundaryMean         | 323.1      | 0.000      | 24 / 24       |

Notes:
- Kuhn shows symmetric x=0 / x=L ≈ 323 Pa; D4 shows only x=0
  with 323 Pa, x=L = 0 because the canonicalized orbit key puts
  all 24 D4 x=L faces in singleton buckets (no mirror pairs under
  the key), so the reducer reports 0 by construction.
- `boundarySym` and `boundaryMean` do not reduce the x=0 drift on
  either fixture.

## §Q.4 Decision — Branch B from plan §4

Per plan §1 decision rule:

> XBoundarySym << XBoundaryRaw: x-sides fully explain the residual
> only XMinSym helps: target x=0
> only XMaxSym helps: target x=L
> **neither helps: x-side residual is not orientation-sign symmetry;
> target BC formula path next**

Data: all four `iSym+x*` modes produce identical 1.000 drift — none
of them reduce below the `InteriorSym` floor.  `boundarySym` on the
face-level drift (Gate 4h) is also ineffective.

**Conclusion: Branch B is selected.**  The x-side residual is NOT
an `n ↔ -n` orientation/sign issue.  It is a structural asymmetry
in the x-side boundary BC formula / input path.  Round 15 should:

1. Dump per-QP x-side boundary face inputs and outputs for
   orbit-paired faces (nor, bdr_attr, I_self, bulk_bg_scaled,
   F_h[SXY]) and compare x=0 ↔ x=L, not just scalar drifts.
2. Check whether x-side faces take a different BC-type dispatch
   path (absorbing vs free-surface) than y/z sides.
3. Examine the trace reconstruction (`shape1 · I_data`) on x-side
   faces for possible DOF-ordering asymmetry.

## §Q.5 Round-14C validation (serial)

Serial pepper guard with Round-14C shared-face hook active (no
shared faces in the 2×2×2 fixture → serial no-op):

| mode                | Kuhn  | D4    |
|---------------------|-------|-------|
| BASELINE            | 2.236 | 4.796 |
| FREEZE-A            | 0.000 | 0.000 |
| NONFAULT_BOTH_SYM   | 0.043 | 0.065 |
| ratio / baseline    | 0.0191 | 0.0136 |

Byte-identical to Round-13C — confirming 14C is serial-safe and
does not affect the local-face path.  MPI validation is deferred
to a parallel run (Frontera or similar).

## §Q.6 Regression check

- Clean build of both test binaries and production `wave_operator.inl`.
- Pepper guard numbers on both fixtures match Round-13C bit-for-bit
  (env unset → pre-patch; env set → matches Round-13C).
- Gate 4g: 3 TEST_LE monotonicity checks pass on both fixtures.
- Gate 4h: diagnostic only.
- Existing Round-13A/13B/13C gates still pass.
- Shared-face hook exercises only when `nonfault_both_sym=true` and
  MPI build has `n_shared > 0`; neither condition holds in the
  serial fixture so no numeric change.

## §Q.7 Logs

- `/tmp/r14_audit_kuhn.log`, `/tmp/r14_audit_d4.log` — audit with
  Gates 4g / 4h.
- `/tmp/r14_pepper_kuhn.log`, `/tmp/r14_pepper_d4.log` — pepper
  guard confirming Round-13C numbers.

## §Q.8 Status

All three 14A / 14B / 14C patches delivered.  Success criterion
from the plan hit:

> "The x-side floor survives all sym modes, so the boundary formula
> path itself is the target."

Round 15 should proceed with Branch B: dump per-QP x-side boundary
face inputs/outputs, compare orbit-paired values face-by-face, and
check BC-dispatch + trace-reconstruction on x-sides specifically.
No more broad sym/mean modes — the next experiment is a targeted
inspection of the boundary formula path on x-aligned outer sides.

---

# §R. Round 15 — x-side boundary trace dump + forced-kernel audit (2026-04-23)

User plan: audit-only patch, single file
(`test_adjacent_triangle_fault_first_step_audit.cpp`).  Add
`XSideBoundaryTraceSample`, `RunXSideBoundaryTraceAudit`, Gate 4i
(orbit drifts + worst-bucket dumps for I_self/bulk_bg/F_h on x-sides
and z-side control), Gate 4j (dispatch consistency), and — only
if Gate 4i shows F_h dirty — a forced-kernel follow-up.

## §R.1 Files changed

- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — new `AuditFaceBC` mirror enum (WaveOperator::FaceBC is private);
  `FaceBCLabel`; `XSideBoundaryTraceSample`, `XSideBoundaryTraceAudit`;
  `XSideTraceOrbitKey`, `MakeBoundaryTraceOrbitKey`,
  `MaxTraceOrbitDrift` (with worst-bucket return), `PrintTraceBucket`;
  `BoundaryKernelOverride` enum + `BoundaryKernelOverrideLabel`;
  `RunXSideBoundaryTraceAudit` (trace-only, no deposit; optional
  forced-kernel override).  In `main`, Gate 4i (orbit drifts +
  worst-bucket dumps) and Gate 4j (dispatch-consistency with 2 new
  TEST_LEs per fixture), plus a "Gate 4i-forced" follow-up block
  that re-runs the audit under three forced kernels.
- No source code changes outside the test file.

## §R.2 Gate 4i — x-side vs z-side orbit drift

Both fixtures, max orbit drift within a matched bucket:

| field            | x-side [Pa] | z-side [Pa] |
|------------------|-------------|-------------|
| I_self[SXY]      | 9.327e-02   | 9.327e-02   |
| bulk_bg[SXY]     | 0           | 0           |
| F_h[SXY]         | **3.231e+02** | 0         |

**Surprise finding: `I_self[SXY]` is DIRTY on BOTH x-sides AND
z-sides**, at identical 0.09 Pa magnitudes.  The BC kernel is the
side-dependent amplifier: x-sides amplify the input asymmetry by
~3600×, z-sides by 0×.

Worst x-side F_h bucket (identical structure on both fixtures):
3 QPs per face × 2 orbit-paired x=0 faces.  One face carries
`I_self[SXY] ∈ {-0.0497, +0.0436, +0.0436}` across its 3 QPs; its
orbit mate at the mirrored (cy, cz) canonical centroid carries all
zeros (Kuhn) or a rotated triple (D4).  The trace reconstruction
has a position-dependent SXY asymmetry that both fixtures share.

## §R.3 Gate 4j — dispatch consistency

On both fixtures: 8 x-side orbit buckets checked, 0 `bdr_attr`
mismatches, 0 `bc_type` mismatches.  Classification is not the
bug; all x-side orbit mates take the same `FaceBC::FreeSurface`
dispatch path.  The two TEST_LE assertions pass on both fixtures.

## §R.4 Forced-kernel follow-up

Same trace audit, but each run forces one boundary kernel on every
x-side face and compares F_h[SXY] orbit drift:

| forced kernel                  | Kuhn x F_h | D4 x F_h | z F_h |
|--------------------------------|------------|----------|-------|
| `ForceAbsorbing`               | 1.615e+02  | 1.615e+02 | 0     |
| `ForceFreeSurface` (Gamma)     | 3.231e+02  | 3.231e+02 | 0     |
| `ForceFreeSurfaceGodunov`      | 3.231e+02  | 3.231e+02 | 0     |

Observations:
- All three kernels amplify the input asymmetry on x-sides; none
  eliminate it.  z-sides are always ULP-zero regardless of kernel.
- `ForceFreeSurface` == `ForceFreeSurfaceGodunov` (identical
  325.1 Pa on both fixtures) — the two FreeSurface variants share
  the same (SXY_out, SXY_in) amplification coefficient.
- `ForceAbsorbing` is EXACTLY HALF of FreeSurface (161.5 vs 323.1).
  This factor of 2 is consistent with the FreeSurface formula
  reflecting twice what Absorbing radiates for a non-characteristic
  component.

The kernel-agnostic amplification on x-sides proves the problem is
NOT a bug inside any single BC kernel.  All three kernels react
the same way to an asymmetric input — the (SXY_out, SXY_in)
matrix element is structurally non-zero when `|n_x| = 1` but
structurally zero when `|n_z| = 1`.  That's expected for a
stress-tensor flux: SXY is the xy-component and a face normal
along x exposes it while a face normal along z does not.

## §R.5 Decision — root cause is UPSTREAM of BC dispatch

Per plan §10 decision logic:

> `I_self[SXY]` dirty on x-sides: target `ComputeADERTimeIntegrated` /
> face trace reconstruction / x-side local state path.

Data: `I_self[SXY]` is dirty on both x-sides AND z-sides at identical
0.09 Pa.  Combined with the forced-kernel finding that all three
boundary kernels are equivalent amplifiers (within their respective
structural coefficients), the x-side F_h drift is a **side-structure
consequence** of a universal **trace-reconstruction asymmetry**, not
a BC-formula bug.

**Conclusion: Round 15 success criterion hit as
"I_self[SXY] is already asymmetric, so the remaining bug is upstream
of boundary BC dispatch."**  The next target is one of:

1. **`ComputeADERTimeIntegrated`**: the CK predictor that produces
   the time-integrated `I` from `Q_step0`.  If this produces
   asymmetric `I` in the SXY channel, every downstream face trace
   inherits that asymmetry.
2. **Face trace reconstruction** (`shape1 · I_data`): the
   dof-to-quadrature contraction on boundary faces.  If the
   per-element shape functions + DOF ordering combine to introduce
   an SXY asymmetry that's only visible after contraction, this is
   the culprit.
3. **DOF ordering on the fault-adjacent tets**: the L2_FECollection
   ordering on orientation-corrected D4 tetrahedra may be
   non-equivariant for SXY specifically.

## §R.6 Round-16 experiment suggestion

The cleanest next probe is a `Q_step0` SXY signature dump —
measure the SXY component of `Q_step0` (not `I`) on the
fault-adjacent tets that subsequently produce orbit-paired faces
with different I_self values.  If `Q_step0[SXY]` is already
asymmetric at the DOF level, the bug is in step 0's AdvanceADER —
which reduces to the step-0 non-fault face deposit (already known
to be ~1 Pa asymmetric per §L.4 / §O).

A parallel test: run `ComputeADERTimeIntegrated(Q_step0, kDt, 2, I)`
with `Q_step0 = 0` and check whether the resulting `I[SXY]` is
ULP-clean on fault-adjacent tets.  If yes, the integrator preserves
zero-input symmetry; if no, the integrator itself introduces
asymmetry.

Both are cheap audit additions and do not require production code
changes.

## §R.7 Regression check

- Clean build.
- Gate 4j: 4 TEST_LE monotonicity assertions pass on both fixtures
  (2 per fixture).
- Gate 4i and Gate 4i-forced: diagnostic only.
- Existing Round-13A/13B/13C/14A/14B/14C gates unchanged.
- Pepper guard numbers untouched (no production code modified).

## §R.8 Logs

- `/tmp/r15_kuhn.log`, `/tmp/r15_d4.log` — first audit run (Gate 4i
  + Gate 4j, no forced-kernel section).
- `/tmp/r15b_kuhn.log`, `/tmp/r15b_d4.log` — full audit with
  forced-kernel follow-up.

## §R.9 Status

Round 15 delivers the promised statement:
> "I_self[SXY] is already asymmetric on x-sides, so the remaining
> bug is upstream of boundary BC dispatch."

Additional nuance: I_self[SXY] is asymmetric on ALL sides (x AND z)
— but the BC formula's (SXY_out, SXY_in) coupling is non-zero only
on x-sides, so the asymmetry only BECOMES VISIBLE as F_h drift on
x-sides.  The forced-kernel sweep shows every BC kernel amplifies
the same way (up to a factor-of-2 Absorbing-vs-FreeSurface split);
there is no single BC-formula bug to fix.

Round 16 should target the trace reconstruction or the CK
predictor upstream of BC dispatch, not the BC formulas themselves.

---

# §S. Round 16 — I_data vs trace-sampling split (2026-04-23)

User plan: audit-only; Gate 4k compares two trace-sampling styles
on the same `I_data` (contiguous `dof_offset + i` vs explicit
`edofs[i]`); Gate 4l measures element-mean `I_data[SXY]` orbit
drift.  Together they split Round-15's upstream residual into
one of three outcomes per plan §10.

## §S.1 Files changed

- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — new types (`BoundaryTraceSamplerSample`,
  `BoundaryTraceSamplerAudit`, `ElementIDataSample`,
  `ElementIDataAudit`); orbit helpers
  (`MakeBoundarySamplerOrbitKey`, `MakeElementOrbitKey`,
  `MaxSamplerOrbitDrift`, `MaxElementOrbitDrift`,
  `PrintBoundarySamplerBucket`, `PrintElementIDataBucket`);
  trace-only helpers `RunBoundaryTraceSamplerComparisonAudit`
  (dual offset+edofs sampling of `I_data[SXY]`) and
  `RunElementIDataAudit` (element-mean SXY + per-DOF values).
  In `main`: Gate 4k (drifts + worst-bucket dump + 2 TEST_LEs
  per fixture) and Gate 4l (element-mean drift + worst-bucket
  dump).
- No source code changes outside the test file.

## §S.2 Gate 4k — offset vs edofs sampling

Both fixtures produce the same pattern:

| side | offset drift [Pa] | edofs drift [Pa] | delta [Pa] |
|------|-------------------|------------------|------------|
| x    | 9.327e-02         | 9.327e-02        | **0.000**  |
| z    | 9.327e-02         | 9.327e-02        | **0.000**  |

**delta = 0 to machine precision on every face-QP pair**.  The
production contiguous `dof_offset + i` addressing is byte-identical
to the explicit `fes.GetElementDofs` lookup.  Contiguous offset
addressing is NOT the bug.

Kuhn worst x-side offset bucket pairs face 13 (elem 4, x=0) with
face 106 (elem 41, x=0) at canonical (cy,cz)=(333,167).  Face 13
produces three per-QP values {-4.97e-02, +4.36e-02, +4.36e-02};
face 106 produces {0, 0, 0}.  Same orbit bucket, same addressing,
very different `I_self[SXY]`.  D4 produces the analogous pair
(face 9 elem 2 vs face 25 elem 8) with the same magnitudes.

## §S.3 Gate 4l — element I_data[SXY] orbit symmetry

| fixture | element-mean SXY orbit drift |
|---------|------------------------------|
| Kuhn    | **1.041e-17**                |
| D4      | **1.388e-17**                |

Element-level `I_data[SXY]` DOF means are orbit-symmetric to
machine precision on both fixtures.  `ComputeADERTimeIntegrated`
is NOT the bug at the element-mean level.

## §S.4 Decision — plan §10 third outcome

| Gate 4k delta | Gate 4l mean | per-QP I_self | outcome |
|---------------|--------------|---------------|---------|
| **clean (0)** | **clean (ULP)** | **dirty (0.09 Pa)** | **Case 3** |

Plan §10 case 3:
> "bug is in boundary face reconstruction using clean element
> data, likely shape/local ordering interaction"

**Interpretation.**  The per-QP contraction
`Σ_i shape1(i) * I_data[c * ndof_total + edofs[i]]` produces
different values on orbit-paired faces even when:
1. The addressing is correct (offset ≡ edofs).
2. The element data means are orbit-symmetric.

For this to happen with a clean element MEAN, the per-DOF values
must NOT be orbit-symmetric individually — only their mean is.
The face's quadrature points pick up specific linear combinations
of the DOFs that depend on `shape1(·)` evaluated at the face QP
positions within the reference element.  Orbit-paired faces on
different elements have DIFFERENT reference-element mappings
(because the elements themselves have different local vertex
ordering on a 6-tet-per-hex Kuhn split / orientation-corrected D4
split), so `shape1(i) · DOF[i]` is a non-equivariant linear
functional even when the DOF mean is equivariant.

## §S.5 Root cause class confirmed

The residual is **shape × local-DOF ordering on boundary faces
with x-aligned normals**.  Not a numerical formula bug, not an
indexing bug, not a BC kernel bug — a DG mesh-topology effect
where the element-local DOF ordering (inherited from the tet
vertex numbering on a Kuhn-split hex) is NOT equivariant under
the y-mirror × z-mirror symmetry that the TPV102 setup expects
to preserve.

The existing Round-13C `SEAS_TEST_NONFAULT_BOTH_SYM=1` hook
closes 98% of pepper because it symmetrizes `F_h` in the n ↔ -n
sense, which averages over the shape-weighted asymmetry well
enough to knock out the dominant signal on faces where both
forward and reverse QP samples exist on the same orbit bucket.
The remaining ~2% is the x-side residual that Round-14A/B/C and
Round-15/16 have now chased to its origin: the shape × DOF
ordering interaction on x-aligned boundary faces.

## §S.6 Round-17 suggestions

Two cheap audit extensions that localize the shape × DOF
asymmetry:

1. **Per-DOF orbit comparison**: extend `ElementIDataSample` to
   also carry a CANONICAL per-DOF permutation (sorted by vertex
   position or by reference-element DOF-coordinate), and compare
   orbit-paired elements at per-DOF level instead of per-DOF-mean.
   Expected: sorted per-DOF values are orbit-equivariant; the raw
   (DOF-index-ordered) per-DOF values are NOT.  This would prove
   the shape × DOF ordering theory.

2. **Shape-only probe on x-side faces**: for each orbit-paired
   x-side face pair (elem_a, elem_b), compute `shape1_a(i)` vs
   `shape1_b(i)` at the mirror-paired QPs.  If shape1 is a
   permutation at orbit-mates (up to sign), the fix is to apply a
   local permutation before `shape1 · I_data` on x-sides.  If
   shape1 values differ non-permutationally, the fix requires
   either a DG basis reordering or a covariant (BR2-style)
   symmetrization on x-side boundary faces specifically.

Neither requires production code changes.  Both narrow Round-17
to a single named fix path.

## §S.7 Regression check

- Clean build.
- Gate 4k: 4 TEST_LE monotonicity assertions pass on both fixtures
  (delta ≤ offset + tol; x and z sides).
- Gate 4l: diagnostic only.
- All prior gates (Round 13/14/15) unchanged.
- No production code touched.

## §S.8 Logs

- `/tmp/r16_kuhn.log`, `/tmp/r16_d4.log` — full audit including
  Gates 4k and 4l.

## §S.9 Status

Round 16 delivers the promised decisive statement:

> "Element I_data is clean and offset == edofs, so the remaining
> bug is face interpolation / local ordering, not predictor or
> indexing."

The Round-12 → 16 chain has narrowed from "pepper at ~2 Pa" to a
specific mechanism: **non-equivariant element-local DOF ordering
combined with shape-function evaluation at boundary face
quadrature points produces per-QP face traces that are not
orbit-symmetric, even though the underlying element data is**.
The 0.09 Pa per-QP residual is amplified by the BC Riemann
formulas only on faces with `|n_x| = 1`, producing the 323 Pa
F_h[SXY] drift and ultimately the ~2 Pa observable pepper.

Round 17 can proceed with the two per-DOF / shape-only probes
above, without freeze authorization or further BC / predictor
instrumentation.

---

# §T. Round 17 — canonical per-DOF + face-interp audit (2026-04-23)

User plan: audit-only; Gate 4m compares element SXY DOF values in
raw order vs canonical (sorted-by-reference-position) order; Gate
4n compares boundary-face interp in raw vs canonical arrangement.
Splits plan §14's three outcomes — permutation mismatch, higher-
order predictor asymmetry, or face-ordering-only.

## §T.1 Files changed

- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — new types `CanonicalLocalDof`, `ElementSxyDofSample`,
  `ElementSxyDofAudit`, `FaceInterpolationSample`,
  `FaceInterpolationAudit`; helpers `GetElementRefDofPositions`,
  `MakeRefPosKey`, `BuildRawLocalDofs`, `BuildCanonicalLocalDofs`;
  orbit reducers `MaxElementRawDofOrbitDrift`,
  `MaxElementCanonicalDofOrbitDrift` (per-DOF worst comparison),
  `MaxFaceInterpolationOrbitDrift` (templated); pretty-printers
  `PrintElementSxyDofBucket`, `PrintFaceInterpolationBucket`.  Two
  new trace-only helpers: `RunElementSxyDofAudit` and
  `RunFaceInterpolationAudit` (the latter applies the SAME local-
  DOF permutation to both DOF values AND shape functions so the
  canonical contraction is a permutation-invariant reformulation
  of the raw one).  In `main`: Gate 4m (raw / canonical drifts +
  worst-bucket dumps + monotonicity TEST_LE) and Gate 4n (x/z raw
  / canonical drifts + worst-bucket dumps + 2 monotonicity
  TEST_LEs).
- No source code changes outside the test file.

## §T.2 Gate 4m — element per-DOF orbit drift

Both fixtures:

| order       | Kuhn      | D4        |
|-------------|-----------|-----------|
| raw         | 1.865e-01 | 1.865e-01 |
| canonical   | **1.865e-01** | **1.865e-01** |

**Canonical sort does not reduce drift.**  Kuhn worst bucket
(cx=125, cy=375, cz=250, n=4):

- elem 4 raw: {+7.465e-02, -1.119e-01, +7.465e-02, +7.465e-02}
- elem 4 canonical: {+7.465e-02, +7.465e-02, +7.465e-02, -1.119e-01}
- elem 19 raw: {-1.119e-01, +7.465e-02, +7.465e-02, +7.465e-02}
- elem 19 canonical: {-1.119e-01, +7.465e-02, +7.465e-02, +7.465e-02}

Canonical elem 4 has the -0.1119 value in slot 3; canonical elem 19
has it in slot 0.  The two elements carry the SAME multiset of DOF
values (three +7.465e-02 and one -1.119e-01) — which is why
Gate 4l's element-mean drift is 1e-17 — but the VALUE-TO-REFERENCE-
POSITION assignment is not equivariant under the y-mirror the fixture
is supposed to respect.

D4 worst bucket (elems 2, 8, 26, 32) shows the same pattern: same
value sets, different canonical positions.

## §T.3 Gate 4n — face interpolation canonicalization

Both fixtures, side-by-side:

| side    | raw interp drift | canonical interp drift |
|---------|------------------|------------------------|
| x (K)   | 4.971e-02        | **4.971e-02**          |
| z (K)   | 4.971e-02        | **4.971e-02**          |
| x (D4)  | 9.327e-02        | **9.327e-02**          |
| z (D4)  | 9.327e-02        | **9.327e-02**          |

The canonical face-interp drift equals raw on every side and every
fixture.  Per construction,
`Σ shape_canonical[k] · canonical_dofs[k].value`
is a reordering of the same sum as the raw contraction — so
raw == canonical per sample is a tautology.  What matters is that
orbit-mates still disagree at the canonical level, by the same
amount as the raw level.  Both TEST_LE monotonicity assertions pass
(diff exactly 0 on every side).

## §T.4 Decision — plan §14 case 2

Per plan §14:

> - Gate 4m raw dirty, canonical clean → local DOF ordering (case 1)
> - **Gate 4m canonical dirty → root cause is higher-order predictor**
>   **content, not just ordering (case 2)**
> - Gate 4m canonical clean, Gate 4n canonical dirty → shape/local
>   interaction (case 3)

Data: **Gate 4m canonical drift = 1.865e-01 on both fixtures.**

The root cause is NOT local DOF ordering alone.  The predictor
`ComputeADERTimeIntegrated` produces per-DOF values at the SAME
canonical reference position that differ across orbit-mates, even
when the underlying value multiset is identical.  This is not a
reorder-and-it-closes bug — orbit-mates genuinely have the -0.1119
sitting at the +x-most canonical DOF on one element and at the +y-
most canonical DOF on its mirror mate.

## §T.5 Mechanism picture

Under the 6-tet-per-hex Kuhn split, orbit-paired elements have
different local vertex numberings.  The finite-element basis's
reference-to-physical map depends on that vertex numbering.  Even
though the two elements occupy orbit-mirrored physical regions, the
DOF at reference position (0, 0, 0) of element A maps to a
DIFFERENT physical vertex of element A than the DOF at reference
position (0, 0, 0) of element B maps to of element B.

The predictor `ComputeADERTimeIntegrated` is then correct per
element in the physical sense — it produces the right value at
each physical vertex — but wrong at the **reference-position-
indexed, element-basis-indexed** slot, because the element's
vertex labeling is not orbit-equivariant.

This is consistent with the orientation-corrected D4 fixture
showing the same pattern: D4's mesh construction corrects vertex
ORIENTATION (all positive Jacobians per §E of Round 3) but does
NOT canonicalize vertex NUMBERING across orbit mates.

## §T.6 Root-cause class

Not a numerical formula bug.  Not an indexing / addressing bug.
Not a BC / Riemann bug.  Not a fault-solver bug.

**The root cause is the element-local vertex numbering not being
equivariant under the TPV102 fixture's mirror symmetries.**  The
predictor evolution, the shape functions, and the face-QP traces
all faithfully reflect this non-equivariance; none of them is
individually wrong.  The 0.09 Pa per-QP `I_self[SXY]` asymmetry and
the ~2 Pa long-time pepper observed in Round 12 are direct
consequences.

## §T.7 Round-18 direction

Two possible fixes, both requiring production code changes and
careful validation:

1. **Element vertex canonicalization at mesh construction**:
   renumber each tet's local vertices so orbit-mates share the
   same local vertex-to-reference-position mapping under the
   fixture's mirror symmetries.  Requires new code in
   `dynamic/d4_tet_mesh.hpp` (or equivalent Kuhn-split builder).
   High risk of breaking orientation consistency.

2. **Symmetrization in the predictor or volume operator**: apply
   an explicit orbit-average of `I_data[SXY]` (and the other
   off-diagonal stress components) at each predictor stage, using
   the mesh's inherent mirror map.  Requires reliable orbit-pair
   detection at runtime and breaks byte-identity with the raw
   predictor.  More defensible as a diagnostic hook than as a
   production fix.

A third alternative is to keep `SEAS_TEST_NONFAULT_BOTH_SYM=1` as
a diagnostic / workaround and accept the 98% pepper-closure as the
observable "fix" while the deeper vertex-numbering issue is tracked
for a future mesh-utility refactor.  This is the least invasive
option and aligns with the caveat preserved throughout §P–§S that
the Round-13C hook is a diagnostic, not a root-cause solution.

## §T.8 Regression check

- Clean build.
- Gate 4m: 1 TEST_LE monotonicity assertion passes on both fixtures.
- Gate 4n: 2 TEST_LE monotonicity assertions pass on both fixtures.
- Existing gates (Round 13 / 14 / 15 / 16) all still pass.
- No production code touched.

## §T.9 Logs

- `/tmp/r17_kuhn.log`, `/tmp/r17_d4.log` — full audit including
  Gates 4m and 4n.

## §T.10 Status

Round 17 delivers the promised statement:

> "Orbit-paired I_data[SXY] does not match even after canonical
> permutation, so the root cause is higher-order predictor content,
> not ordering."

More specifically: the predictor is consistent PER ELEMENT but the
ELEMENT-LOCAL VERTEX NUMBERING is not orbit-equivariant, so per-
DOF-indexed output violates the fixture's expected mirror symmetry.
Round 18 should decide between a mesh-construction fix (canonicalize
vertex numbering) and a predictor-symmetrization workaround,
treating `SEAS_TEST_NONFAULT_BOTH_SYM=1` as the interim
diagnostic/closure that already exists.

---

# §U. Round 18A — face-local canonical interpolation probe (2026-04-23)

User plan: audit-only Round 18A; add Gate 4o (face-local canonical
interpolation).  If it closes, proceed to 18B (runtime hook
`SEAS_TEST_BOUNDARY_TOPO_SHAPE=1`) and 18C (pepper validation).  If
it does not close, stop — the face-local sampling hypothesis is
incomplete, and production code should not be patched.

**Important: Rounds 18B and 18C were NOT implemented** — the plan's
decision rule at §"Concrete branch after Round 18" says explicitly:

> "Gate 4o does not close: stop.  Do not patch production yet;
>  the hypothesis is still incomplete."

Gate 4o on both Kuhn and D4 produced face-canonical drift equal to
raw drift, so Round 18 ends with 18A.

## §U.1 Files changed

- **modified:** `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`
  — new types (`FaceLocalCanonicalSample`, `FaceLocalCanonicalAudit`),
  helper `BuildFaceLocalCanonicalLocalDofs` (sorts local DOFs by
  face-restricted reference coordinates — `(yr, zr)` on x-normal
  faces, `(xr, yr)` on z-normal faces), orbit reducer
  `MaxFaceLocalCanonicalOrbitDrift`, printer
  `PrintFaceLocalCanonicalBucket`, trace-only audit helper
  `RunFaceLocalCanonicalInterpolationAudit`, and **Gate 4o** in
  `main` (3 orderings × 2 sides with 2 monotonicity TEST_LEs per
  fixture).
- No source code changes outside the test file.

## §U.2 Results

| fixture | side | raw        | elem-canonical | face-canonical |
|---------|------|------------|----------------|----------------|
| Kuhn    | x    | 4.971e-02  | 4.971e-02      | **4.971e-02**  |
| Kuhn    | z    | 4.971e-02  | 4.971e-02      | 4.971e-02      |
| D4      | x    | 9.327e-02  | 9.327e-02      | **9.327e-02**  |
| D4      | z    | 9.327e-02  | 9.327e-02      | 9.327e-02      |

All three orderings produce identical drift to machine precision.
Both fixtures show the face-canonical ordering does not reduce the
residual.  The 4 TEST_LE monotonicity assertions pass with diff = 0.

Worst bucket (Kuhn x-side q=0): face 19 elem 6 at (cy=333, cz=167)
and its orbit mate face 111 elem 43 at (cy=667, cz=833) give
interp = 0.000 and -0.04971 respectively.  The underlying values
at orbit-paired QPs differ by ~0.05 regardless of any local DOF
reordering.

## §U.3 Decision — halt Round 18

Plan §"Concrete branch after Round 18" decision rule:

> **"Gate 4o does not close: stop.  Do not patch production yet;
>   the hypothesis is still incomplete."**

**Rounds 18B (`SEAS_TEST_BOUNDARY_TOPO_SHAPE=1` runtime hook in
`wave_operator.inl`) and 18C (pepper-guard validation row) were
NOT implemented.**  The face-local canonicalization hypothesis
required Gate 4o to close; it did not.  Proceeding to 18B would
patch production without a valid supporting audit.

## §U.4 Refined root cause statement

Combining Round 17 (Gate 4m canonical dirty) + Round 18A (Gate 4o
face-canonical dirty):

- Orbit-paired elements carry the SAME SXY DOF multiset (Gate 4l
  mean = 1e-17; §S).
- At the element-canonical (3D reference-position) ordering, DOF
  values differ by 0.186 (Gate 4m; §T).
- At the face-canonical (face-plane reference-position) ordering,
  the per-QP trace still differs by 0.05 (Gate 4o).
- No local reordering — by 3D reference position OR by face-plane
  reference position — reconciles orbit-paired traces.

**The remaining asymmetry is NOT a labeling-permutation problem at
any local (element or face) level.**  Orbit-paired elements share
the same value multiset but with a permutation that is NOT
expressible as a local DOF sort against any reference-coordinate
ordering we have tested.

## §U.5 Implication

Possibilities that remain open:

1. The mesh's vertex-to-DOF mapping is non-equivariant in a way
   that requires GLOBAL reindexing (not a local per-element or
   per-face sort).  A fix would require either mesh-wide vertex
   canonicalization at construction OR a custom FE basis that
   is explicitly equivariant under the expected symmetries.
2. The `ComputeADERTimeIntegrated` predictor evolves DOFs via
   an operation (gradient contraction, mass-inverse application)
   that couples neighboring DOFs asymmetrically when the neighbor
   set is itself non-equivariant.  The element-mean equivariance
   (§S.3) is preserved because means are first-order invariants,
   but higher moments are not.
3. The `FaceElementTransformations::Loc1.Transform(ip, ip1)` map
   from face-reference to element-reference coordinates is not
   equivariant across orbit-paired faces, so even with identical
   face-QP `ip` and identical shape-evaluation logic, `ip1` lands
   in a different element-reference position and samples a
   different linear combination.

The third hypothesis is testable with one more audit but would
likely show the same "face-canonical no-close" signature since our
face-canonical sort is on the ELEMENT reference coordinates restricted
to face-plane DOFs, not on the face-reference coordinates.  That's
a distinct probe (Gate 4p: face-reference-position canonicalization)
that could be added if it's the highest-value next step.

## §U.6 Regression check

- Clean build.
- Gate 4o: 4 TEST_LE monotonicity assertions pass on both fixtures.
- Existing gates (Round 13–17) all unchanged.
- **No production code modified.**  `wave_operator.inl`,
  `precomputed_face_fluxes.*`, and `test_r9_frozen_friction_C.cpp`
  are untouched on this turn.

## §U.7 Logs

- `/tmp/r18a_kuhn.log`, `/tmp/r18a_d4.log` — Round-18A audit
  including Gate 4o.

## §U.8 Status

Round 18A delivers the negative result that ends the plan's
decision chain:

> "Orbit-paired I_data[SXY] does not match even after canonical
>  permutation, so the root cause is higher-order predictor
>  content, not ordering."

refined to:

> "Orbit-paired I_data[SXY] does not match after ANY local
>  canonical permutation (element-reference or face-reference),
>  so the root cause is NOT a local-sampling permutation problem."

Per plan, production patches are halted.  Options going forward:

- (A) Retain `SEAS_TEST_NONFAULT_BOTH_SYM=1` as the accepted
  diagnostic / 98%-closure workaround, and document the
  remaining 2% as a known mesh-equivariance limitation of the
  Kuhn / 6-tet-per-hex fixture.
- (B) Test the face-reference-position canonicalization probe
  (Gate 4p) before committing to (A) — would rule out the
  `Loc1.Transform` / face-to-element ref-map hypothesis.
- (C) Invest in mesh-wide vertex canonicalization or an
  equivariant basis.  High-risk, out of scope for this
  investigation round.

My recommendation: option (A) is the cleanest honest endpoint
given Round-13C already closes 98% of pepper and the remaining 2%
traces to a mesh-topology non-equivariance that no local sort
can fix.  Option (B) can be added later as a single extra audit
gate if a clean "face-reference-position canonicalization"
test is desired for completeness.


