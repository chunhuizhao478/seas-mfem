# Implementation Plan: Bimaterial central-flux contrast guard (TPV31 sigma_n leak)

Worktree: `.claude/worktrees/fix-mixedflux-hetero` (branch `worktree-fix-mixedflux-hetero`)
Date: 2026-06-06
All math is ASCII (no LaTeX) per project convention.

Companion diagnosis (read first):
- `debug_document/tpv31_debug_document/drdg3d_central_flux_comparison_2026-06-06.md`
- this session's findings: sigma_n leak is EXCLUSIVE to the bimaterial step depths
  (2400/5000/10000 m); mid-layer flat; scalar TPV102 ~5 MPa bounded vs matrix TPV31
  >100% collapse.

## Overview

On fault-adjacent corridor faces that cross a STRONG material (impedance) contrast,
stop using the non-dissipative central flux and fall back to the EXISTING
impedance-weighted bi-material upwind (Pelties Godunov).  This removes the TPV31
on-fault sigma_n collapse (and the downstream tau_strike error) that is localized to
the depth-discontinuity stations, while leaving central flux on every
uniform-material corridor face (so the low-dispersion benefit and the across-fault
[[v_n]] protection are untouched).  Default OFF (byte-exact); opt-in + tunable.

### Why this is the right fix (one paragraph)

The fault is vertical (y=0); the material layers are horizontal (z=const).  The
faces where the material actually JUMPS are the horizontal layer-interface faces in
the VOLUME (z-normal), perpendicular to the fault.  Central flux there is both
non-dissipative AND not impedance-weighted, so it mishandles the impedance contrast,
seeds a spurious [[v_n]] at the fault near that depth, and LSW friction ratchets it
into the sigma_n collapse.  The bi-material upwind handles the contrast correctly
(verified: matrix+upwind+symmetric-mesh control is flat at all depths).  We touch
ONLY those volume contrast faces; the fault face flux and the fault-parallel
(y-normal) central corridor faces are unchanged, so the across-fault [[v_n]]
protection central provides is preserved.  This matches drdg3d's regime (central only
within uniform material; contrasts via impedance-weighted upwind).

## Constraints

- BYTE-EXACT DEFAULT (project regression contract; CLAUDE.md "every new flux mode
  MUST default OFF"): `contrast_tol < 0` => guard disabled => zero behavior change.
  TPV102 (scalar) and any homogeneous bi-material run are unaffected (contrast == 0
  everywhere => no face is ever reclassified, even with the guard enabled).
- DO NOT modify the fault-face imposed-state/friction flux, nor the scalar
  `WaveOperator` path.  The change is confined to the bi-material corridor
  classification.
- Central-set lifecycle (R-1408): `SetMixedFluxMode` (base, `wave_operator.inl:1639`)
  calls `BuildCentralFluxFaceSet_` (base) then, via the override
  `BimaterialWaveOperator::SetMixedFluxMode` (`bimaterial_wave_operator.hpp:138`),
  `BuildPerFaceCentralFluxMatrices_`.  The contrast filter MUST run inside the
  bi-material subclass (only it knows per-element material) and MUST erase reclassified
  faces from `central_flux_face_set_` so the dispatch in `InteriorFaceFlux_`
  (`bimaterial_wave_operator.inl:732`) / `SharedInteriorFaceFlux_` (`:761`) routes them
  to the existing upwind branch.
- MPI: shared (cross-rank) faces need the neighbour material from
  `shared_face_neighbour_material_` (`bimaterial_wave_operator.hpp:250`), populated by
  `ExchangeBiMaterialNeighbours_`.  Both ranks owning a shared face MUST classify it
  identically (deterministic contrast test on the same two material triples).
- NO local full-mesh runs (memory `feedback-no-local-mesh-runs`): local acceptance =
  compile + unit tests; production validation = Frontera A/B, user-submitted.
- Build in the worktree with the MFEM override (memory `worktree-build-mfem-dir-override`):
  `make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN <target>`,
  `MAIN=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver`.
- Do NOT revert a prior fix without citing the debug doc + approval (CLAUDE.md).

## Key math

Impedances per element: `Zp = rho*cp`, `Zs = rho*cs` (accessors `GodunovFlux::GetZp()`
/ `GetZs()`, `godunov_flux.hpp:172-173`).  Per-face contrast:

    cZp = |Zp_e1 - Zp_e2| / max(Zp_e1, Zp_e2)
    cZs = |Zs_e1 - Zs_e2| / max(Zs_e1, Zs_e2)
    contrast = max(cZp, cZs)
    strong  := (contrast_tol >= 0) && (contrast > contrast_tol)

TPV31 separation (this session, from the layer table; mesh h=50 m):
- STEP contrasts:    cZs = 12.4 / 14.8 / 17.7 %   (cZp 9.7-19.8 %)  -> must go UPWIND
- GRADIENT/element:  cZs ~ 0.36 %  (h=50 m; scales with h)          -> must stay CENTRAL
- ~35x gap => default `contrast_tol = 0.05` sits cleanly in the gap and is
  insensitive; safe down to h~200 m (gradient 1.4 %).  Steps are h-INDEPENDENT;
  gradient ~ |dZ/dz|*h/Z shrinks as h->0, so refining only widens the gap.

Flux dissipation identity (`godunov_flux.hpp:64-68`): `Interior - Central =
0.5*|A|*(Q_self - Q_nbr)`.  Central discards this SPD dissipation term; across an
impedance contrast it is also not impedance-weighted (the bi-material upwind uses the
Pelties star state `Q* = (2 Z_L/(Z_L+Z_R)) v`, `godunov_flux_bimaterial.cpp:24-33`).

---

## Phase 1: contrast utility + tolerance member/setter (foundation)

### Goal
A pure, tested impedance-contrast predicate exists, and `WaveOperator` carries a
`mixed_flux_contrast_tol_` (default -1 = disabled) with a public setter, with no
behavior change anywhere (no caller passes a non-negative tol yet).

### Files to Modify
- `dynamic/wave_operator.hpp` — add the member + setter/getter (base, so the driver
  can set it through a base `WaveOperator&`):
  - In the `protected:` block (near `cfl_rk_aware_`, `:943`):
    `real_t mixed_flux_contrast_tol_ = -1.0;`  (negative => guard disabled).
  - In the `public:` block (near `SetCflRkAware`, `:219`):
    `void SetMixedFluxContrastTol(real_t tol) { mixed_flux_contrast_tol_ = tol; }`
    `real_t GetMixedFluxContrastTol() const { return mixed_flux_contrast_tol_; }`
- `dynamic/godunov_flux_bimaterial.hpp` — add a static pure predicate (no state):
  `static bool IsStrongContrast(const GodunovFlux& a, const GodunovFlux& b, real_t tol);`
  returning `tol >= 0 && max(|Zp_a-Zp_b|/max(Zp_a,Zp_b), |Zs_a-Zs_b|/max(Zs_a,Zs_b)) > tol`.
- `dynamic/godunov_flux_bimaterial.cpp` — implement `IsStrongContrast` using
  `GetZp()/GetZs()`; guard `tol < 0 -> return false`; guard divide-by-zero
  (max impedance > 0 always for non-acoustic; if either is 0, return false).

### Detailed Requirements
1. `mixed_flux_contrast_tol_` default `-1.0` MUST mean "disabled" everywhere it is
   read (Phase 2).  `0.0` is a valid (very aggressive) tolerance, NOT "disabled".
2. `IsStrongContrast` is symmetric in its two arguments and side-effect free.
3. No existing call site sets a non-negative tol in Phase 1 (so the tree is
   byte-identical after Phase 1).

### Interfaces
- `WaveOperator<MeshType>::SetMixedFluxContrastTol(real_t)` / `GetMixedFluxContrastTol()`.
- `BimaterialFlux::IsStrongContrast(const GodunovFlux&, const GodunovFlux&, real_t) -> bool`.

### Edge Cases to Handle
- `tol < 0` => always false (disabled).
- Equal materials (homogeneous / within-layer gradient below tol) => false.
- Acoustic / zero impedance => false (SAFS has no acoustic regions; do not abort here —
  the existing acoustic guards live in the flux builders).

### Acceptance Criteria
- [ ] Compiles in the worktree (MFEM override) — both serial `Mesh` and `ParMesh` TUs.
- [ ] New unit test `seas_test_bimaterial_contrast_guard` Phase 1 block:
      `IsStrongContrast` returns false for tol<0; false for equal materials; true for
      the TPV31 5 km pair (cZs ~14.8%) at tol=0.05; false for the linear-gradient
      adjacent pair (cZs ~0.36%) at tol=0.05.
- [ ] `git diff` shows no change to any existing dispatch/flux numerics (member added,
      unused).

### Dependencies
- Depends on: nothing. Required by: Phase 2, 3.

---

## Phase 2: apply the guard in the bi-material central build + diagnostics

### Goal
With `mixed_flux_contrast_tol_ >= 0`, fault-adjacent corridor faces whose two elements
are a strong contrast are removed from `central_flux_face_set_` (so they dispatch to
the existing bi-material upwind), for BOTH local interior and shared (MPI) faces; a
one-line histogram of corridor-face contrasts is printed on rank 0.

### Files to Modify
- `dynamic/bimaterial_wave_operator.inl` — in `BuildPerFaceCentralFluxMatrices_`
  (`:515`):
  - LOCAL faces loop (`:550-571`, iterating `central_flux_face_set_`): after
    `ResolveFaceFluxOperands_(mesh_face, flux_e1, flux_e2, nor)` returns the two
    `const GodunovFlux*`, if `BimaterialFlux::IsStrongContrast(*flux_e1, *flux_e2,
    mixed_flux_contrast_tol_)` then record `mesh_face` in a local
    `std::vector<int> reclassified` and `continue` (do NOT build central matrices).
  - SHARED faces loop (`:660-690` region): get the local element's material via
    `FluxForElem_(elem1)` and the neighbour material from
    `shared_face_neighbour_material_.at(mesh_face_idx)` (a `{lambda,mu,rho}` or
    `{vp,vs,rho}` triple — confirm the stored convention at
    `bimaterial_wave_operator.hpp:250` and `ExchangeBiMaterialNeighbours_`), build a
    temporary `GodunovFlux` for the neighbour, and apply the SAME `IsStrongContrast`
    test; if strong, record + `continue`.
  - AFTER both loops: `for (int f : reclassified) central_flux_face_set_.erase(f);`
    Then the existing post-build invariant assert
    (`per_face_central_flux_.size() == central_flux_face_set_.size()`,
    `bimaterial_wave_operator.inl:877`) must still hold (the erased faces were never
    inserted into `per_face_central_flux_`).
- `dynamic/bimaterial_wave_operator.inl` — emit a rank-0 diagnostic: contrast
  histogram over ALL faces that WERE in `central_flux_face_set_` before filtering
  (bins e.g. [0,1%),[1,5%),[5,10%),[10,20%),[>=20%]), plus
  `[mixed-flux-contrast] tol=<t>  reclassified <N>/<M> corridor faces to upwind`.
  Use `MPI_Reduce` (sum bins to rank 0) so the histogram is global.

### Detailed Requirements
1. The filter MUST run only when `mixed_flux_contrast_tol_ >= 0`; otherwise the loops
   are byte-identical to today (no `IsStrongContrast` call, no erase, no histogram).
2. Reclassified faces MUST NOT have central matrices built (they fall through to the
   upwind branch in `InteriorFaceFlux_`/`SharedInteriorFaceFlux_` via the existing
   `central_flux_face_set_.count(mesh_face) > 0` check, which is now false for them).
3. Shared-face classification MUST be deterministic and rank-consistent: both ranks
   compute `IsStrongContrast` on the same ordered material pair (local vs neighbour);
   since the predicate is symmetric and both ranks see the same two triples, both
   reach the same verdict.  Add an assert/log if `shared_face_neighbour_material_`
   lacks an entry for a central shared face (mirror the existing IMPL-4 fail-loud at
   `:84`).
4. The histogram counts the PRE-filter corridor faces (so the bimodal gap is visible
   and one can confirm `tol` lands in the empty bin).

### Interfaces
- No new public interface; internal to `BuildPerFaceCentralFluxMatrices_`.
- Relies on Phase 1 `IsStrongContrast` and `mixed_flux_contrast_tol_`.

### Edge Cases to Handle
- `central_flux_face_set_` empty (mixed_flux=none) => loops no-op => histogram prints
  "0/0" only if tol>=0 (or skip the print when the set is empty).
- All corridor faces uniform (TPV102 scalar path does not reach this subclass; a
  homogeneous bi-material would) => 0 reclassified => byte-exact.
- A corridor face touching the fault-x-material-crossing line: still classified purely
  by its own two elements' contrast (no special-casing) — correct, that is exactly the
  z-normal layer face we want on upwind.
- Linear-gradient layer at coarse mesh: if a gradient face's contrast exceeds tol it is
  (harmlessly) put on upwind; the histogram makes this visible so tol can be raised.

### Acceptance Criteria
- [ ] `seas_test_bimaterial_contrast_guard` Phase 2 (extend
      `test_bimaterial_mixed_flux_dispatch.cpp` fixture, a small mesh with a material
      contrast crossing the corridor):
      - tol = -1 (disabled): `central_flux_face_set_` identical to a reference build
        (no faces removed) — BYTE-EXACT.
      - tol = 0.05 on a fixture with a strong-contrast corridor face: exactly the
        strong-contrast faces are removed from the set; uniform-material corridor faces
        retained.
      - homogeneous-material bi-material fixture, tol=0.05: zero faces removed
        (== disabled) — BYTE-EXACT.
- [ ] Parallel test (np=2, extend `tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp`
      or `test_bimaterial_mixed_flux_shared.cpp`): a strong-contrast SHARED face is
      removed on BOTH ranks (consistent classification); a uniform shared face retained.
- [ ] `make test` (serial unit suite) still green except the 3 known pre-existing
      failures (memory `preexisting-worktree-test-failures-2026-05`).
- [ ] Rank-0 histogram prints and shows the bimodal distribution on the test fixture.

### Dependencies
- Depends on: Phase 1. Required by: Phase 3 (CLI), Phase 4 (regression/Frontera).

---

## Phase 3: config + CLI wiring + banner

### Goal
The tolerance is reachable from the TOML `[numerics]` and a `--mixed-flux-contrast-tol`
CLI flag, printed in the startup banner / `--verify-dispatch`, default disabled.

### Files to Modify
- `spatial/code/spatial_friction.hpp` — add to the numerics config struct (near
  `mixed_flux`, `:133`): `real_t mixed_flux_contrast_tol = -1.0;` (disabled default).
- `spatial/code/spatial_friction.cpp` — parse it (near `:1096` where `mixed_flux` is
  read): `cfg.numerics.mixed_flux_contrast_tol = toml_real(n,
  "mixed_flux_contrast_tol", -1.0);`  Validation: allow any value; `< 0` => disabled;
  optionally warn if `>= 0` while `mixed_flux == "none"` (no effect).
- `drivers/spatial_dyn_driver.cpp`:
  - CLI parse (near `--mixed-flux` at `:576`):
    `const std::string cli_ctol = GetStringArg(argc, argv, "--mixed-flux-contrast-tol", "");`
    and after the config merge (`:692`): if non-empty, `cfg.numerics.mixed_flux_contrast_tol = std::stod(cli_ctol);`
  - Apply it BEFORE `wave.SetMixedFluxMode(...)` (`:1276`) so the build sees the tol:
    `wave.SetMixedFluxContrastTol(cfg.numerics.mixed_flux_contrast_tol);`
    (SetMixedFluxMode triggers the central build that reads the tol — order matters).
  - Banner (`:884` area): print
    `"mixed flux contrast tol: " << cfg.numerics.mixed_flux_contrast_tol
     << (tol < 0 ? " (disabled)" : "")`.

### Detailed Requirements
1. Precedence: CLI overrides TOML (same pattern as `--mixed-flux` at `:692`).
2. Ordering: `SetMixedFluxContrastTol` MUST be called before `SetMixedFluxMode`
   (the override that rebuilds the central set + matrices reads the tol).  If a later
   restart/static-cycle re-invokes `SetMixedFluxMode`, the tol persists (member), so
   the rebuild stays consistent (R-1408 lifecycle).
3. No effect on the scalar driver path / non-spatial drivers (they don't set it;
   member stays -1).

### Interfaces
- TOML key: `[numerics].mixed_flux_contrast_tol` (float, default -1).
- CLI: `--mixed-flux-contrast-tol <float>`.

### Edge Cases to Handle
- Flag given on a `mixed_flux=none` run: no central set => no effect; banner still
  prints the value; optional one-line warning.
- Malformed CLI value: `std::stod` throws -> catch and `MFEM_ABORT` with a clear
  message (match the `ParseMixedFlux` abort style at `:177`).

### Acceptance Criteria
- [ ] `seas_test_tpv_config_parse` extended: TOML `mixed_flux_contrast_tol = 0.05`
      parses to 0.05; absent => -1.
- [ ] `--mixed-flux-contrast-tol 0.05` overrides the TOML in a `--dry-run`
      `--verify-dispatch` (grep the banner line) on the TPV31 config (no full run).
- [ ] Default run (no flag, TOML absent) banner shows "(disabled)" and the dispatch is
      byte-identical to before this plan.

### Dependencies
- Depends on: Phase 1, 2. Required by: Phase 4.

---

## Phase 4: regression guards + Frontera A/B job

### Goal
A local hard-assert test locks the guard's behavior, and a Frontera sbatch performs the
decisive TPV31 A/B (guard off vs on, normal 50 m mesh) to confirm the sigma_n leak is
removed at the interface depths.

### Files to Create / Modify
- `tests/unit/test_bimaterial_contrast_guard.cpp` (NEW; Makefile target
  `seas_test_bimaterial_contrast_guard` mirroring `seas_test_bimaterial_central_flux`
  at `Makefile:5376`) — the Phase 1+2 assertions, promoted to hard `TEST_*`:
  - mechanism: across the TPV31 5 km contrast, for a normal-velocity jump,
    `dot(jump, F_upwind_self - F_central_self) > 0` strictly (central discards
    dissipation), and `-> 0.5*Zp*jump^2` in the homogeneous limit (within 1e-9 rel).
  - face-set: tol=0.05 removes exactly the strong-contrast corridor faces; tol=-1 and
    homogeneous-material are byte-exact.
- `jobs/tpv31_spatial/mixedflux_rk45/tpv31_p1_rk45_mixedflux_50m_normal_contrastguard.sbatch`
  (NEW) — copy of the gold job adding `--mixed-flux-contrast-tol 0.05`; OUT dir keyed
  by job id; same 12N/600r; `--verify-dispatch`; footer greps the histogram line and
  runs `visualize_results.py --tol-peak/--tol-rms`.
- `debug_document/tpv31_debug_document/tpv31_contrast_guard_findings_2026-06-06.md`
  (NEW) — record the before/after sigma_n-by-depth table once the Frontera A/B returns.

### Detailed Requirements
1. The unit test is the LOCAL gate; the Frontera A/B is the PRODUCTION gate (no local
   full-mesh).
2. Acceptance for the Frontera A/B (recorded in the findings doc, user-submitted):
   sigma_n excursion at dp024/dp050/dp100 drops from 16/9/64 MPa toward the mid-layer
   floor (target: < a few MPa, comparable to TPV102's bounded ~5 MPa generic leak);
   tau_strike at those stations tracks the SeisSol reference; V_strike/slip unchanged
   (the guard must not damp the rupture: peak V within ~5%).
3. Note explicitly in the findings doc that this guard does NOT address (a) the smaller
   generic asymmetric-mesh leak (~5 MPa, the upwind_unstructured plan target) nor
   (b) the static t=0 station-sampling mu-side offset at dp050 (a station-output issue).

### Acceptance Criteria
- [ ] `seas_test_bimaterial_contrast_guard` passes (mechanism + face-set + byte-exact).
- [ ] Frontera A/B (user-submitted): guard-on sigma_n flat-ish at 2.4/5/10 km vs the
      catastrophic guard-off baseline; benchmark tolerance gate improved at those
      stations; rupture (V_strike/slip) within 5% of guard-off.
- [ ] No existing TPV*/BP5 byte-exact regression broken (guard defaults disabled).

### Dependencies
- Depends on: Phase 1-3.

---

## Testing Strategy
- Per phase: build in the worktree (MFEM override), run the targeted unit test, then
  `make test` for the serial suite (accept the 3 known pre-existing failures).
- Mechanism oracle (local): the flux-primitive dissipation identity across the 5 km
  impedance contrast (Phase 1/4 test) — proves central discards dissipation there and
  upwind restores it, using the production flux builders directly.
- Classification oracle (local): the face-set membership tests (Phase 2) on a small
  bi-material fixture, serial + np=2 — proves exactly the right faces switch, byte-exact
  when disabled.
- Production oracle (Frontera, user-submitted): TPV31 sigma_n-by-depth A/B.
- Cross-check the scheme intent against drdg3d (central only within uniform material;
  contrasts via impedance-weighted upwind) — see the companion comparison doc.

## Risk Assessment
- R1 (medium): residual leak after the guard. The guard removes the catastrophic
  contrast-driven leak but NOT the generic asymmetric-mesh leak; expect TPV31 to land
  near TPV102's ~5 MPa level, not zero. Mitigation: documented as separate; combine
  with the upwind_unstructured theta-blend / symmetrization later if needed.
- R2 (low): coarse-mesh false positives (gradient faces flipped to upwind). Mitigation:
  the histogram makes it visible; tol is tunable; default 0.05 is safe to h~200 m.
- R3 (medium, MPI): shared-face classification disagreement across ranks ->
  non-conservative dispatch. Mitigation: symmetric predicate on the same two material
  triples + the np=2 test + fail-loud on a missing `shared_face_neighbour_material_`
  entry.
- R4 (low): ordering bug — tol applied AFTER SetMixedFluxMode -> guard silently inert.
  Mitigation: Phase 3 requirement #2 + the `--verify-dispatch` banner + the Frontera
  histogram grep (reclassified count must be > 0 for TPV31).
- R5 (low): transition reflections at the new central<->upwind boundaries on the layer
  faces. Mitigation: those faces are few (3 depths) and the boundary type already
  exists at the corridor edge; if visible, switch Phase-2 hard reclassification to a
  theta-blend (future extension, out of scope here).
- R6 (tricky code): `shared_face_neighbour_material_` storage convention
  ({lambda,mu,rho} vs {vp,vs,rho}) — confirm before building the temporary neighbour
  `GodunovFlux` (read `ExchangeBiMaterialNeighbours_` end-to-end first; Rule 5).
