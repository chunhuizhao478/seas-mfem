# Code Review: PLAN_spatial_seas_quasidynamic_driver_2026-05-31 (fresh review — round 2, 2026-06-02)

## Review Scope
- Plan: `miniapps/seas/document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.{md,pdf}`
  (PDF rebuilt from the `.md` via `build_pdf.sh`; 24 pp).
- Prior reviews: `REVIEW.md` (2026-06-01, R-001…R-006 — remediated) and `REVIEW_2026-06-02.md` (R-007…R-010).
- This round re-executes all three passes on the **changed sections** after the R-007…R-010 fixes, re-verifies
  every newly-cited source fact, and hunts for defects introduced by the edits.
- Status: still **plan-stage** (no implementation; the worktree holds only the plan + reviews). Findings are
  against the plan's normative contract.

## Pass 1 — Are R-007…R-010 resolved? (re-verified against source)

### [R-007] CRITICAL — BP5 prestress reproducibility → RESOLVED
- Phase 3 acceptance (plan §Phase 3, "Parity test") now **splits the comparison**: `a_values_/dc_values_/
  eta_values_/V_init_vec_` must match path(a) to 1e-10 (resolver reproduces them); `tau_pre_/sigma_n_per_dof_`
  are explicitly handled — it states `ComputeBP5Params` does **not** populate `sigma_n_per_dof_`
  (verified: `fault_geometry.hpp:1391` fills only `a/eta/dc/tau_pre/V_init`; `sigma_n_per_dof_` is written
  only at `:666/:685-691`) and that BP5 `tau_pre_` is spatially heterogeneous (verified: `bp5_params.hpp:333`
  `tau0_vec` uses `a_of_x2_x3` + a nucleation `delta_tau` taper), and forbids asserting a uniform
  `FaultLocalPrestress` `tau_pre_` against the analytic one.
- Phase 5 §B now carries a BP5-parity note ("prestress MUST be PER-DOF … NOT uniform FaultLocalPrestress").
- Phase 7 config now mandates a per-DOF BP5 prestress source (`bp5_analytic` or `tau0` SidecarHDF5) and
  documents that none of the four existing stress kinds (verified in `spatial_stress.hpp` +
  `spatial_dyn_driver.cpp:1337-1458`) computes the formula, and that `ResolveRateState` produces no `tau_pre`
  (verified: `RateStatePerDOFParams` has `a,b,Dc,V_init,f_0,V_0,eta,sigma_n_eff` — no prestress field,
  `spatial_friction.hpp:687`).
- **Verdict:** the parity acceptance criterion is now satisfiable, and the Phase-7 gate has a viable mechanism.
  Eta-source caveat (auto-from-material μ/c_s, `spatial_friction.cpp:2024`) is documented. Resolved.

### [R-008] MODERATE — `ParseQDSolverType` return type → RESOLVED
- Phase 1 now declares `mfem::seas::SolverType ParseQDSolverType(const std::string&)` returning enum values
  (`SolverType::CG_AMG`, …) with an explicit note that Phase 2 passes the result into the ctor's `SolverType`
  parameter. Verified the ctor parameter type (`elasticity_operator.hpp:111`) and that `SolverType` is a free
  enum in `mfem::seas` (`:46`). Phase 2 call site (plan §Phase 2, `bc, dg_method, ParseQDSolverType(...),
  domain_config`) now type-checks. Resolved.

### [R-009] MODERATE — checkpoint dynamic-code dependency → RESOLVED
- Phase 6 now specifies a QD-native checkpoint (`io/seas_qd_checkpoint.hpp`) and forbids reusing
  `io/tpv104_checkpoint.hpp`, citing the `dynamic/fault_face_flux.hpp`/`DOFData` include (verified
  `tpv104_checkpoint.hpp:40`) and the `Q`+`std::vector<DOFData>` API (verified `:74-77`, `:147-150`).
- **Cascading fix applied this round:** the "What is reused as-is (no changes)" list (plan §Architecture)
  still named `io/tpv104_checkpoint.hpp` — that contradiction is now removed and annotated. No other reuse
  reference remains (grep-confirmed). Resolved.

### [R-010] MODERATE — `SetRateStatePerDOF` V_init expansion → RESOLVED
- The method now takes `const mfem::Vector& init_vel_dir`; the size-`N` scalar `rs.V_init` (verified
  `spatial_friction.cpp:1941`) is decomposed into the `2N` `(dip,strike)` interleave along `init_vel_dir`
  (matching `bp5_params::V_init_vec`, `:279`), and the spec explicitly forbids keying off `tau_pre_` (NaN at
  that point — verified `fault_geometry.hpp:286-287`). The `2N` requirement matches the ctor assert
  (`rate_state_fault.hpp:161`). Signature is consistent across Architecture (§new code 3), Phase 3 step 2,
  the Phase 3 driver snippet (`Vector init_vel_dir(2); … (0,1) // BP5: pure strike` → `SetRateStatePerDOF(rs,
  init_vel_dir)`), the parity test, and the Appendix. Resolved.

## Pass 2 — New bugs introduced by the edits?
- **`init_vel_dir = (0,1)` for BP5** — correct: BP5 is pure strike-slip, layout is `(2i)=dip,(2i+1)=strike`,
  so magnitude → strike, dip → 0. Consistent with `bp5_params::V_init_vec` and `CLAUDE.md` ("dip slip ≈ 0").
- **Signature consistency** — all five `SetRateStatePerDOF` call/declaration sites carry `init_vel_dir`; no
  one-arg form remains (grep-confirmed). No stale `std::string ParseQDSolverType`, no stale
  `interleave consistent with tau_pre_`, no stale Phase-7 `FaultLocalPrestress (BP5 τ0/σ_n)` (grep-confirmed).
- **All newly-cited line numbers verified** against source: `fault_geometry.hpp:1391/:286-287`,
  `bp5_params.hpp:333/:279`, `spatial_friction.cpp:1941/:2024`, `rate_state_fault.hpp:161`,
  `tpv104_checkpoint.hpp:40/:74-77/:147-150`, `elasticity_operator.hpp:46`. No new mis-citation.
- **PDF** rebuilds (only the pre-existing 2× `U+2077` "⁷" missing-glyph warnings from the untouched
  `~10⁷-DOF` overview text — not introduced here).

No new defects found in the changed sections.

## Pass 3 — Residual / lower-priority observations (NOT blockers)
- **[OBS-1] LOW — `bp5_analytic` stress source is now on the BP5-parity critical path but under-specified.**
  Phase 7 introduces it as the fix for R-007 but leaves its exact form (new `StressSourceKind` enum value +
  parser + a `ComputeParams` functor wrapping `bp5_params::tau0_vec`) to the implementer. This is additive and
  low-risk, but since it now gates Phase 7 it deserves a one-line spec in Phase 3 or a dedicated sub-phase
  (e.g. "Phase 3c — `bp5_analytic` stress source"). Not a correctness bug; a planning-completeness note.
- **[OBS-2] LOW — `init_vel_dir` is uniform per call.** The plan passes a single `(dip,strike)` direction for
  the whole fault, which is right for BP5 (globally strike-slip). A future SAF problem with a curved fault may
  need a per-DOF loading direction; flagged only so it is not silently assumed global. Out of scope now.
- These are observations, not findings — no fix required to proceed.

## Summary
- Findings re-checked: R-007 (CRITICAL) ✔ resolved; R-008/R-009/R-010 (MODERATE) ✔ resolved.
- Cascading issues from the fixes: 1 found + fixed this round (reused-as-is list still naming `tpv104_checkpoint.hpp`).
- New CRITICAL/MODERATE defects: **0**.
- Low observations: 2 (OBS-1, OBS-2) — non-blocking.
- Plan compliance with the in-tree code it targets: **FULL** for every signature/line the plan now cites
  (re-verified this round); the remaining open items are genuine physics/design decisions the plan already
  flags as open questions (SAF loading function, heterogeneous material, per-DOF b/f0/V0).
- **Verdict: PASS** — R-001…R-010 are all resolved and source-consistent. The plan is implementable as a
  normative contract. Recommended (non-blocking) before coding Phase 7: add the one-line `bp5_analytic`
  stress-source spec (OBS-1).

## Unreviewed Areas (unchanged from round 1)
- Physics correctness of the BP5 golden itself (recurrence ~240 yr, dip slip ≈ 0) — gated by Phase 7.
- SAF non-planar far-field loading function (Phase 8 open question #1) — physics design, needs user/SCEC input.
- Makefile object-prerequisite list (Phase 0) — not audited line-by-line; verify no wave OBJs leak in.
- Whether any BP5-scalar read path (`sigma_n_bp5_`/`Vp_bp5_`) is reachable once `IsSAFSMode()` — recommend an implementer assert.
