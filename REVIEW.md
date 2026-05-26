# Code Review: 2026-05-26 (Round 2 / final revision) — σ_n compressive normal-stress strength floor

> Round-2 fresh adversarial re-audit of the σ_n-strength-floor implementation
> AFTER the R-001/R-002/R-003 fix pass.  All three review passes re-run from
> scratch on the changed files (rule 10) — verifying the fixes AND hunting for
> new bugs the fixes (or round 1) may have introduced/missed.
>
> Round-1 findings: **R-001** (silent mis-nesting) FIXED; **R-002** (RS floor
> inert) RESOLVED-BY-DOC; **R-003** (no upper bound) FIXED.  One new residual
> finding this round: **R-004** (LOW).

## Review Scope
- Plan: `miniapps/seas/safs/project_7.0_alternative/debug_document/spatial_dynamic_rupture_sliver_blowup_2026-05-25.md` (§ Implementation Plan, Phases 1–4)
- Files reviewed: `spatial/code/spatial_friction.{hpp,cpp}`, `dynamic/fault_face_flux.{hpp,cpp}`, `dynamic/tpv205_friction.hpp`, `dynamic/tpv205_substep_iterator.cpp`, `drivers/spatial_dyn_driver.cpp`, the 7 SAFS configs, `spatial_friction_config_schema.md`, and tests `test_spatial_friction_config.cpp` / `test_tpv205_friction.cpp` / `test_fault_face_flux.cpp`.
- Domain context: `CLAUDE.md` (σ_n>0 compression; byte-exact mandate), `miniapps/seas/CLAUDE.md`, plan §2c/§6b/§Decisions.
- Build/run verification this round: `seas_test_spatial_friction_config` 61/61; `seas_test_tpv205_friction` 82/82; `seas_test_fault_face_flux` 30/30; `seas_test_tpv102_nuc_callback_parity` 18/18 (byte-exact RS, floor disabled). Driver builds; LSW Dc2 **and** RS depthprofile configs both dry-run with banner `sigma_n strength floor: 10 MPa` and no false guard trip.

## Verification of round-1 fixes

- **R-001 (MODERATE) — FIXED.** `spatial_friction.cpp:parse_root` now aborts if `[friction.slip_weakening]` or `[friction.rate_state]` contains `sigma_n_strength_floor_pa`. Verified safe: the loop only dereferences `fr.at(sub)` for sub-blocks that are *present*, and the upstream law-block `MFEM_VERIFY(!has_rs/!has_lsw)` guarantees exactly one sub-block is present and is a valid table by the time the guard runs — no out-of-range/`.at()` throw. Does not false-trip the shipped configs (LSW Dc2 + RS depthprofile both parse). Test `T_32_sigma_n_floor_misnested_aborts` passes (61/61).
- **R-002 (LOW) — RESOLVED BY DOCUMENTATION.** The fix agent correctly declined the v2 code change (flooring the solver argument) because it contradicts plan §Phase 3 decision 2 and widens the byte-exact surface. The schema now documents that the RS floor affects only the slip-rate decomposition and is largely inert (RS uses `abs(σ_n)`; tension strengthens RS friction). Accepted.
- **R-003 (LOW) — FIXED.** A `mfem::out` warning fires when the floor exceeds 1 GPa (exponent-typo guard). Verified silent at 10 MPa. The warning prints on all MPI ranks, consistent with the existing parser warnings (e.g. the `P_p_grad` warning at `:707`) — not a new defect.

## Findings (this round)

### [R-004] LOW [POSSIBLE] [spatial_friction.cpp:parse_root] — R-001 guard misses deeper mis-nesting (`[friction.rate_state.depth_profile]`, `[[...spatial]]`)

**Category:** EDGE_CASE

**Description:**
The R-001 guard checks only the two *direct* law sub-tables (`slip_weakening`, `rate_state`). The floor key nested one level deeper — under `[friction.rate_state.depth_profile]`, `[[friction.rate_state.spatial]]`, or `[[friction.slip_weakening.spatial]]` — is still silently ignored (those parsers also don't reject unknown keys), leaving the floor disabled with no error. This is the same silent-disable class as R-001 but for a less-likely nesting. Severity LOW because the overwhelmingly likely mistake (putting it with the other friction params, i.e. directly in the law sub-block) is now caught; deeper nesting is an obscure mistake.

**Trigger:** `sigma_n_strength_floor_pa = 10.0e6` placed under `[friction.rate_state.depth_profile]` (or a `[[...spatial]]` rule).

**Actual behavior:** parses silently; `cfg.sigma_n_strength_floor_pa == -1.0` (disabled).

**Expected behavior:** abort, as for the direct-sub-block case.

**Suggested fix (optional — extend the R-001 guard to the known deeper tables):**
```diff
       for (const char* sub : {"slip_weakening", "rate_state"})
       {
          if (fr.contains(sub)
              && fr.at(sub).contains("sigma_n_strength_floor_pa"))
          {
             MFEM_ABORT("[friction." << sub << "].sigma_n_strength_floor_pa is "
                        "mis-placed: ... Move it directly under [friction].");
          }
+         // Also catch the one-level-deeper tables that exist in the schema.
+         if (fr.contains(sub))
+         {
+            const auto& sb = fr.at(sub);
+            if (sb.contains("depth_profile")
+                && sb.at("depth_profile").contains("sigma_n_strength_floor_pa"))
+            {
+               MFEM_ABORT("[friction." << sub << ".depth_profile]."
+                          "sigma_n_strength_floor_pa is mis-placed: move it "
+                          "directly under [friction].");
+            }
+         }
       }
```
(`[[...spatial]]` arrays would need an array iteration; given the low likelihood, documenting "the floor is a top-level [friction] key" — already done in the schema — may be sufficient. The fix agent may close R-004 as "documented" rather than code it.)

**Test case:**
```cpp
static void T_33_sigma_n_floor_misnested_depthprofile_aborts()
{
   // RS depth-profile header + floor mis-nested under [friction.rate_state.depth_profile]
   std::string toml = MinimalRSHeaderWithDepthProfile();  // helper analogous to the RS configs
   toml += "[friction.rate_state.depth_profile]\n"
           "param_a_csv = \"a.csv\"\nparam_a_minus_b_csv = \"amb.csv\"\n"
           "depth_units = \"km\"\n"
           "sigma_n_strength_floor_pa = 10.0e6\n";   // WRONG nesting (2 levels deep)
   TEST_ASSERT(ParseAbortsInChild(toml),
               "floor under [friction.rate_state.depth_profile] must abort (R-004)");
}
```
*(If a depth-profile fixture is awkward in the test harness, downgrade/skip — this is LOW and the schema already documents the correct placement.)*

---

## Items checked and found correct (not bugs)

- **R-001 guard control flow is safe** (no `.at()` on an absent/non-table key — see Verification above).
- **Byte-exact linchpin intact:** `FaultFaceFlux::FaultFaceFlux` (`:36`) still never touches `sigma_n_strength_floor_` ⇒ `-1.0` default for TPV102/104/205/BP5; parity suites green this round.
- **All four strength sites still floored, no fifth:** unchanged from round 1 (LSW kernel `tpv205_friction.hpp:173–174`; RS `fault_face_flux.cpp:236`/`:580`); the round-2 change touched only the parser, not the kernels.
- **R-003 warning math:** `/1.0e3` is a heuristic hint only; double division, no truncation.
- **Both target sbatch configs parse post-fix** with the floor active (LSW Dc2 banner "10 MPa"; RS depthprofile banner "10 MPa"), confirming the new guard does not regress the shipped configs.
- **LSW Dc2 cap run is a genuine physics experiment, not a guaranteed fix:** per the plan §Risk-Assessment the floor bounds the *friction-mediated* free-slip feedback, but the central-flux/sliver σ_n-oscillation *source* (config uses `mixed_flux=adjacent`) remains. The job may still blow up if the oscillation itself (not the free-slip) drives it — this is the experiment the plan's Phase-4 acceptance asks for, not a code defect. (Documented in the sbatch script.)

## Summary
- Critical issues: 0
- Moderate issues: 0 (R-001 fixed)
- Low issues: 1 (R-004, POSSIBLE; R-002/R-003 resolved)
- Plan compliance: **FULL** for code Phases 1–3 + Phase 4 configs. Phase-4 cluster blow-up-arrest acceptance is the purpose of the sbatch jobs prepared alongside this review (needs the cluster).
- Verdict: **PASS** — the floor implementation and the round-1 fixes are correct and verified; R-004 is an optional LOW hardening of the R-001 guard (the schema already documents correct placement). Safe to submit the cluster jobs.

## Unreviewed Areas
- Prior uncommitted depth-profile (Phase 11b) feature + large pre-existing `spatial_dyn_driver.cpp` edits — out of scope (reviewed only for interaction with the floor; none found).
- Phase-4 cluster physics validation (V_max arrest) — the two sbatch jobs prepared with this review; cannot be run locally.
- Pre-existing `make test` failures (macOS MPI Bus errors; `ValidateFacetBCTables`/`fault_flux_=0`/sidecar serial aborts; `bp2_serial_smoke` ParaViewOutput API-drift build error) — all in files unmodified by this change.
