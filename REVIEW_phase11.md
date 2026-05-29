# Code Review: Phase 11 — depth-varying rate-and-state a(z)/b(z) — ROUND 2 (post-fix) — 2026-05-29

> Supersedes the round-1 Phase 11 review (recoverable via git history; durable copy
> `miniapps/seas/REVIEW_phase11_2026-05-29.md`).  Round 1 found 1 CRITICAL (R-001) + 4 LOW.
> This round re-executed all three passes from scratch on the changed files AND re-hunted for new
> bugs introduced by the fixes.

## Review Scope
- Plan: `miniapps/seas/document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md`,
  § "Phase 11" (~2047–2396).
- Files re-reviewed (the /code-fix diff):
  - `miniapps/seas/spatial/code/spatial_friction.cpp` — R-001 SRW guard in `resolve_rs_impl`; R-003
    `load_one_depth_csv` trailing-token guard.
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` — R-005 print-derived min/max(a,b) Allreduce + knots.
  - `…/config/spatial_friction_rate_state_safs_projected_stress_depthprofile.toml` — R-002 NOTE rewrite.
  - `miniapps/seas/tests/unit/test_resolve_rate_state_guards.cpp` — new `G_R001` test.
  - `miniapps/seas/tests/unit/test_friction_depth_profile.cpp` — new R-003 case.
- Domain context: `miniapps/seas/CLAUDE.md` (σ_n>0, z<0 below surface, no hardcoded constants, MPI
  collective discipline), round-1 `REVIEW.md`.
- Tests re-run locally (mfem-dev, main-repo libmfem): `resolve_rate_state_guards` 14/14,
  `friction_depth_profile` 24/24, `spatial_friction_resolver` 100/100, `seed_equilibrium_psi_rs` 21/21,
  `ader_tpv102_smoke` 4/4, `tpv102_setup` 24/24, `tpv102_nuc_callback_parity` 18/18,
  `friction_substep_iterator_parity` 36/36 — all PASS; `seas_spatial_dyn_driver` compiles; deliverable
  sbatch `bash -n` OK.

## Round-1 findings — re-verified as fixed (with evidence)
- **R-001 (was CRITICAL) — FIXED.** `resolve_rs_impl` now aborts when `state_evolution ==
  SlipLawStrongRateWeakening` and (`depth_profile.enabled` OR any spatial rule has non-NaN `b`)
  (`spatial_friction.cpp:1909-1929`). Placed after the BoxcarTaper guard, **before** the per-DOF loop
  and the V_w check, so it pre-empts the silent disequilibrium and any later abort. Verified it does
  **not** over-reject: the guard touches only `b` (and `depth_profile`), so scalar-b SRW and per-DOF
  `a`/`Dc`/`sigma_n`/`V_init`/`eta`/`V_w` rules still resolve (test `G_R001` control: scalar-b SRW +
  per-DOF `a` rule resolves with `p.b==b_default`, `p.a==0.012`). The per-rule check `isnan(r.b)` is
  correct (NaN = unset sentinel). Aging path untouched (still reads `d.b`). Always reached before the
  time loop (driver resolves `rs` unconditionally for the RS path), including `--dry-run`.
- **R-002 (LOW) — FIXED.** Config NOTE now states the real ~11 km transition and VW(0–11 km)/
  VS(11–16.6 km) split + super-critical patch, consistent with `[meta]` and the committed
  `*_vwvs11km.csv` (transition where `(a−b)=0` ≈ 11 km; patch `center_z_m=−4965 m` ⇒ VW). The only
  remaining "16.7 km" mention is the correct historical comparison to the *original* CSVs (line 150).
- **R-003 (LOW) — FIXED.** `load_one_depth_csv` reads the leftover token as `std::string`
  (`spatial_friction.cpp:199-208`), so a trailing non-numeric token aborts too. Verified trailing
  whitespace / trailing comma do **not** false-abort (the string extraction hits EOF → empty → pass).
  Test `D5` R-003 case passes.
- **R-004 (LOW) — no change, as recommended.** `DOFData.b` stays `quiet_NaN()` (guarded by R-028).
- **R-005 (LOW) — FIXED.** Driver print-derived block now prints the `a`/`(a−b)` knot tables and the
  resolved `min/max(a)`/`min/max(b)` over fault DOFs. The four `MPI_Allreduce` (`spatial_dyn_driver.cpp:
  1871-1874`) are **outside** the `if (rank==0)` block (line 1876) and inside the all-ranks-identical
  `if (print_derived){ if (!is_lsw && cfg.rate_state && depth_profile.enabled) }` guard, so every rank
  reaches them — no divergence/deadlock. Empty-fault ranks contribute `+inf`/`-inf` (absorbed by
  MIN/MAX). `rs` (declared `:1447`, "filled iff !is_lsw") is populated on this path; `rs.a.Size()==
  rs.b.Size()==N`, no OOB.

## Findings (this round)

### [R2-001] [LOW] [POSSIBLE] [spatial_friction.cpp:parse_rate_state] — SRW + depth_profile incompatibility is caught only at resolve time, not at parse time

**Category:** QUALITY (defense-in-depth / error-timing consistency)

**Description:**
`parse_rate_state` already reads `state_evolution` (sets `out.state_evolution`) and parses the
`depth_profile` block (sets `out.depth_profile.enabled`, loading the CSVs **at parse time** so a bad
file aborts at config load). The R-001 incompatibility (SRW + per-DOF `b`) is, however, enforced only
later in `resolve_rs_impl`. This is **correct for safety** — the driver always resolves `rs` before the
time loop (and before `--print-derived`), so the guard always fires — but it is inconsistent with the
depth-profile design intent ("read at parse time so a misconfig aborts at config load, early and
clear"). A parse-time guard would surface the error at the same point as a missing CSV. Marked POSSIBLE
because it is a hardening/consistency improvement, not a correctness gap (no path reaches the simulation
with SRW + per-DOF `b`).

**Trigger:** A config with `state_evolution = "slip_law_strong_rate_weakening"` + a
`[friction.rate_state.depth_profile]` block. Today: CSVs load successfully at parse, then
`ResolveRateState` aborts. Desired: abort at parse, alongside the other `[friction.rate_state]`
validators.

**Actual behavior:** Aborts at resolve time (still before any time stepping).

**Expected behavior:** Also abort at parse time for a faster, source-consistent error.

**Suggested fix:** add a guard at the end of `parse_rate_state` (after the `depth_profile` block is
parsed and `state_evolution` is known, near `spatial_friction.cpp:~768`):
```diff
       out.depth_profile.profile =
          LoadFrictionDepthProfileCSVs(out.depth_profile);
+      // R2-001: per-DOF b (depth profile) is aging-law only; the SRW iterator
+      // evolves psi with the scalar b_default (see the resolve-time R-001 guard
+      // in resolve_rs_impl).  Fail at config load, consistent with the other
+      // [friction.rate_state] validators.
+      MFEM_VERIFY(out.state_evolution != StateEvolutionKind::SlipLawStrongRateWeakening,
+                  "[friction.rate_state.depth_profile] is not supported with "
+                  "state_evolution=slip_law_strong_rate_weakening (per-DOF b is "
+                  "aging-law only; the SRW iterator uses the scalar b_default).");
    }
```
(Keep the resolve-time guard as the authoritative backstop — it also covers per-DOF `b` spatial
rules and the programmatic `ResolveRateState` API that bypasses TOML parsing.)

**Test case:**
```cpp
// extend tests/unit/test_spatial_friction_config.cpp (string-parse path)
void test_R2001_parse_rejects_srw_plus_depth_profile()
{
   const bool aborted = RunInChild([]() {
      // minimal RS toml string with state_evolution=SRW + a depth_profile block
      // pointing at two real temp CSVs; expect ParseSpatialFrictionConfigString
      // to MFEM_ABORT at parse (before resolve).
      (void) ParseSpatialFrictionConfigString(kSrwPlusDepthProfileToml);
   });
   TEST_ASSERT(aborted, "parse rejects SRW + depth_profile (R2-001)");
}
```

---

## New-bug hunt (fixes did NOT introduce regressions — specific checks)
- **Guard ordering:** R-001 runs after the BoxcarTaper guard, before the per-DOF loop / V_w check. For
  SRW + boxcar_taper + `b`, the BoxcarTaper guard fires first (both are correct rejections). No
  message-ordering hazard that masks a real error.
- **No over-rejection:** the guard does not block scalar-b SRW, nor per-DOF `a`/`Dc`/`sigma_n`/`V_init`/
  `eta`/`V_w` rules (only `b` + `depth_profile`). Confirmed by the `G_R001` control assertion.
- **R-003 false-positive check:** trailing whitespace and trailing comma do not abort (string
  extraction returns empty at EOF); only a genuine extra token aborts.
- **R-005 collective safety:** all four Allreduce are unconditional on every rank within the
  all-ranks-identical guard; not rank-0-gated; empty-fault ranks safe; no early `return`/`continue`
  between the per-rank reductions.
- **Byte-exactness preserved:** the spatial_friction.cpp edits are confined to the SRW guard and the
  CSV loader — neither is on the TPV102 aging path; `ader_tpv102_smoke` + parity tests still bit-identical.
- **Tests are not false-passing:** `G_R001` control resolves outside `RunInChild` (a spurious abort
  would crash the binary and fail the run); abort cases use the fork harness and the guard is the first
  reachable abort point.

## Summary
- Critical issues: 0
- Moderate issues: 0
- Low issues: 1 (R2-001 — optional parse-time guard for SRW + depth_profile; resolve-time guard already
  guarantees correctness)
- Plan compliance: **FULL** (11a/11b + 11c acceptance criteria; the round-1 CRITICAL is fixed and the
  SRW/per-DOF-`b` boundary the plan implied is now enforced with a test).
- Verdict: **PASS WITH FIXES** — the round-1 CRITICAL (R-001) is fixed and verified, no new
  CRITICAL/MODERATE issues were introduced, and the byte-exact TPV regression holds. The Phase-11
  deliverable (aging-law depth profile) is ready. R2-001 is an optional defense-in-depth hardening.

## Unreviewed Areas
- Full-mesh / Frontera runtime behavior (repo policy: no local full-mesh runs; verified via unit tests +
  driver compile + reading).
- `io/tpv104_checkpoint.hpp` comment-only update (Phase 11a step 6) — unchanged this round; `b` is a
  recomputed static field, not serialized.
- sbatch variants beyond the dev 8N/400r form (only `bash -n` parse-checked).
