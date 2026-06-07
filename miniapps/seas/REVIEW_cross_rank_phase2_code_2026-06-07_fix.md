# Fix Report — Cross-rank Phase 2 (against REVIEW_cross_rank_phase2_code_2026-06-07.md)

**Date:** 2026-06-07. Review verdict was "mergeable stepping stone", **0 critical**, 3 moderate,
2 minor, 2 suggestion. All 7 findings addressed. No critical ⇒ proceeding after fix per loop policy.

## Changes made
- **P2-001 (MODERATE) + P2-007 (SUGGESTION):** added `RunSharedUpwindMatchesSerial` (gate
  `shared_upwind_matches_serial`) — the BULK UPWIND shared Riemann matrix
  (`per_face_bimaterial_flux_`, the `interior_flux=matrix, mixed_flux=none` production path) on a
  fault-adjacent contrast seam now == the serial 2-sided upwind matrix to 1e-12, for a Coefficient
  mu(x) contrast (proves the re-baseline is the CORRECT bi-material flux) AND a Constant control
  (byte-exact). Corrected the plan's `Constraints` byte-exact claim: byte-exact for Constant always;
  for Coefficient only where the seam is material-continuous; a depth-profile seam cutting across
  depth RE-BASELINES the bulk shared upwind (O(grad·h), bulk only, fault unchanged) — to be
  re-baselined on Frontera.
- **P2-002 (MODERATE):** replaced the deleted stub WARNING with a one-line INFORMATIONAL rank-0
  banner in `ExchangeBiMaterialNeighbours_` (non-Constant + n_shared>0): "cross-rank seam material =
  TRUE peer … re-baselined vs the pre-2026-06 local-side stub."
- **P2-003 (MODERATE):** added a one-time rank-0 WARNING in `BuildPerFaceCentralFluxMatrices_` when
  the contrast guard is OFF (tol<0) + material non-Constant + central set non-empty — the strong
  contrast on a central corridor face builds a non-dissipative central flux silently. (The proper
  guard activation/test is Phase 3; the production config sets the tol in Phase 5.)
- **P2-004 (MINOR):** fixed the two now-false comments ("inert under the seam-continuous stub";
  histogram rationale) and changed the reclassify print from "(local+shared)" to "(local only;
  shared tallied separately in Phase 3)". The actual shared `n_reclass_shared` tally is Phase 3.
- **P2-005 (MINOR):** both accessors now `SetIntPoint` before `EvalAt` (`T.SetIntPoint(&ip)` in
  `MaterialAtLocal_`; `ftr->Elem2->SetIntPoint(&ip_peer)` in `MaterialAtNbr_`). Byte-exact for the
  explicit-ip FunctionCoefficients (parity 9/9), closes the GetIntPoint-reading-coefficient footgun.
  Added gate `coeff_peer_read_independent_of_stale_intpoint` (P2-T4): an `IPProbeCoefficient` reading
  `T.GetIntPoint()` returns the value at `ip_peer` even with the peer transform's int point poisoned.
- **P2-006 (SUGGESTION):** deferred to Phase 5 (drop/relabel the driver's `(seam_continuous=...)`
  banner) — flagged, no Phase-2 action per the reviewer.

## Decisions / deviations
- P2-002: did NOT re-add a "WRONG material" warning — the behavior is now CORRECT; the banner is
  informational only.
- P2-004 shared tally + P2-T2/P2-T3: the `n_reclass_shared` count and the guard-on reclassify gates
  are Phase 3 deliverables (the plan's Phase 3 = activate the contrast guard + tally); fixed only the
  stale wording in Phase 2.

## Verification (re-run after fixes)
- Gate `seas_test_bimaterial_seam_material_np2` (np=2): **119/119** (incl. P2-T1 upwind==serial 0.0
  for Coefficient + Constant; P2-T4 stale-intpoint 0.0; P2-003 warning fires).
- Byte-exact regressions: parity **9/9**, central-flux **6/6**, perside-material **11/11**,
  constant-parity **19/19**, dispatch **10/10**, mixed-flux-shared (np=2) **16/16**.

## Ready for next phase: YES (0 critical; all findings dispositioned; tests green).
