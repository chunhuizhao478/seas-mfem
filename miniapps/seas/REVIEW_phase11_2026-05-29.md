# Code Review: Phase 11 — depth-varying rate-and-state a(z)/b(z) — 2026-05-29

> Supersedes the prior [DIAG-SIGN] speckle-instrumentation review (recoverable via git history).
> Scope is Phase 11 (11a per-DOF `b` channel, 11b two-CSV depth profile, 11c deliverable).
> Durable copy: `miniapps/seas/REVIEW_phase11_2026-05-29.md`.

## Review Scope
- Plan: `miniapps/seas/document/fullelasticity_dev/PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md`,
  § "Phase 11 — Depth-varying rate-and-state a(z)/b(z)" (lines ~2047–2396).
- Files reviewed:
  - `miniapps/seas/dynamic/fault_face_flux.hpp` (DOFData.b)
  - `miniapps/seas/dynamic/tpv102_substep_iterator.cpp` (3 ψ-update sites :186/:348/:530)
  - `miniapps/seas/dynamic/spatial_setup.hpp` (InitializeFaultDOFs_Spatial_RS, SeedEquilibriumPsi_RS)
  - `miniapps/seas/dynamic/tpv102_setup.hpp` (d.b = TPV102Params::b)
  - `miniapps/seas/friction/state_policies.hpp` (aging vs SRW UpdatePsi b-source)
  - `miniapps/seas/dynamic/friction_iterator_factory.cpp` (RS dispatch)
  - `miniapps/seas/spatial/code/spatial_friction.{hpp,cpp}` (interpolant, loader, parser, resolver, R-006 lift)
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` (print-derived depth summary :1817–1919)
  - `miniapps/seas/tests/unit/test_friction_depth_profile.cpp`
  - `miniapps/seas/safs/project_7.0_alternative/config/spatial_friction_rate_state_safs_projected_stress_depthprofile.toml`
- Domain context: `miniapps/seas/CLAUDE.md` (σ_n>0 compression, depth z<0, no hardcoded constants,
  Brent solver), repo memory (DOFData per-DOF channels, byte-exact TPV102 oracle).
- Tests run locally (mfem-dev, main-repo libmfem): `seas_test_friction_depth_profile` 23/23,
  `seas_test_spatial_friction_resolver` 100/100, `seas_test_seed_equilibrium_psi_rs` 21/21,
  `seas_test_resolve_rate_state_guards` 10/10, `seas_test_tpv102_setup` 24/24,
  `seas_test_ader_tpv102_smoke` 4/4, `seas_test_tpv102_nuc_callback_parity` 18/18,
  `seas_test_friction_substep_iterator_parity` 36/36 — all PASS.

## Verified correct (specific, not generic praise)
- **TPV102 byte-exactness holds.** `tpv102_setup.hpp:101` sets `d.b = TPV102Params::b`, and the
  aging policy reads `d.b`; since `d.b == AgingLawPsi(TPV102Params::b,…).GetB()` exactly, the
  iterator change is a true no-op for TPV102 (`ader_tpv102_smoke` + parity tests confirm bit-identity).
- **Per-DOF `b` is live end-to-end for the aging law.** `SpatialRule.b` exists
  (`spatial_friction.hpp:30`), `parse_spatial_rule` reads it (`:569`), `resolve_rs_impl` applies it
  (`:1952`), `InitializeFaultDOFs_Spatial_RS` sets `d.b=rs.b(i)` (`:342`), and BOTH the Phase-2
  adapter (`tpv102_substep_iterator.cpp` ×3) AND the Phase-5 unified `RateStateAgingPolicy`
  (`state_policies.hpp:58-62`) read `d.b`. The factory dispatches the unified `RateStateAgingIterator`.
- **R-006 lift is correct for the aging path:** per-DOF `b` accepted, per-DOF `f_0`/`V_0` still
  rejected (`spatial_friction.cpp:1946`). R-011 honored (no per-DOF `a<b` requirement, `:2011-2015`).
- **Interpolant flat-clamp (never linear extrapolation)** is correct (`:136-137`), binary search
  brackets correctly including at-knot queries, no div-by-zero (strictly increasing x).
- **`b = a − (a−b)` is computed at query time on independent grids** (`spatial_friction.hpp:464-466`);
  km→m via config `depth_to_m` (no hardcoded 1000, `:758-765`); `depth = max(0,−z)` (`:1907`).
- **Absent-profile path is byte-identical** to the scalar path (`rs.b(i)==b_default`; resolver test).
- **NaN-`b` default is safe:** `verify_dof_data_uninitialized` only checks `lsw_d_c`/`Zs_plus`; no
  LSW/TPV104/checkpoint path reads `d.b` (grep), and `SeedEquilibriumPsi_RS` guards it (R-028, `:433`).
- **Both `ResolveRateState` overloads** delegate to the single `resolve_rs_impl`, so both get profile
  seeding. The driver depth-summary block is gated under `if (print_derived)` and is MPI-safe
  (all ranks reach the Allreduce; rank-0-only printing).

## Findings

### [R-001] [CRITICAL] [spatial_friction.cpp:resolve_rs_impl + friction/state_policies.hpp] — per-DOF `b` is seeded but NOT evolved under the SRW state law → silent t=0 disequilibrium (Phase-11-INTRODUCED regression, config-reachable, no guard)

**Category:** BUG (silent wrong results; regression on a first-class friction law)

**Severity rationale (why CRITICAL, not MODERATE):**
- **Silent wrong physics.** Produces a t=0 disequilibrium with **no abort/warning** — the worst failure
  mode for a scientific code (CLAUDE.md regression criteria: V_max spikes / spurious slip).
- **It is a regression.** `git log -S "rs.b(ii)"` shows the `SeedEquilibriumPsi_RS` change to per-DOF
  `rs.b(ii)` was introduced by the Phase-11 commit `156ea34`. **Before Phase 11** the seed used scalar
  `blk.b_default`, which *matched* the SRW iterator's scalar `L.GetB()` — consistent. Phase 11 changed
  the seed but not the SRW evolution, breaking a previously-correct path.
- **Reachable with no code change.** One config key (`state_evolution =
  "slip_law_strong_rate_weakening"`) + a `depth_profile` (the Phase-11 feature) triggers it; the driver
  calls `SeedEquilibriumPsi_RS` **unconditionally** for the RS path (`spatial_dyn_driver.cpp:1592`),
  and the factory dispatches the SRW iterator for that law (`friction_iterator_factory.cpp:80`).
- **First-class law.** SRW = TPV104, a benchmark this very plan wires; depth-varying friction is the
  whole point of Phase 11, so an SRW user wanting depth variation hits exactly this.
- **Magnitude.** Seed offset `Δψ = (b_perDOF − b_default)·ln(V0/V_init) ≈ Δb·ln(1e-6/1e-12) ≈ 13.8·Δb`;
  for `Δb ~ 5e-3` that is `Δψ ~ 0.07` — large versus `f0 = 0.6`, i.e. a real disequilibrium, not noise.

**Description:**
Phase 11a lifted R-006 to accept per-DOF `b` for **all** rate-state configs, and changed
`SeedEquilibriumPsi_RS` to seed ψ from the per-DOF `rs.b(ii)` (`spatial_setup.hpp:456,473`). However,
only the **aging** policy reads `d.b` (`state_policies.hpp:58-62`); the **slip-law strong-rate-weakening
(SRW)** policy reads the **scalar** `L.GetB()` (= `blk.b_default`) at `state_policies.hpp:82-89`
(TPV104 left untouched, per the plan). When a config selects `state_evolution =
"slip_law_strong_rate_weakening"` **and** enables a `depth_profile` (or sets a per-DOF `b` spatial
rule), the resolver produces depth-varying `rs.b`, the equilibrium seed uses it, but the SRW
time-evolution uses scalar `b_default`. The seed and the dynamics disagree on `b`, so the fault is
**not in equilibrium at t=0** and the depth-varying `b` is silently ignored by the SRW dynamics. No
guard exists (the factory only checks SRW needs `V_w`; nothing checks `b`). Without a profile,
`rs.b(i) == b_default` everywhere so seed and evolution still agree — which is exactly why the bug is
latent until the Phase-11 depth profile (or a per-DOF `b` rule) makes `rs.b` non-scalar.

**Trigger:**
A `[friction.rate_state]` block with `state_evolution = "slip_law_strong_rate_weakening"` plus either
`[friction.rate_state.depth_profile]` enabled, or any spatial rule with a numeric `b`.

**Actual behavior:**
`resolve_rs_impl` accepts the per-DOF `b`; `SeedEquilibriumPsi_RS` seeds ψ with `rs.b(ii)` per-DOF;
the SRW iterator evolves ψ with scalar `blk.b_default`. No abort, no warning — silent disequilibrium.

**Expected behavior:**
Either reject the combination (faithful to the plan, which scopes SRW `b` as scalar and explicitly
leaves TPV104/SRW untouched), or thread `d.b` into the SRW path. Because the plan deliberately keeps
SRW `b` scalar (Phase 11a step 2, §Constraints) and routing `d.b` into SRW would also require setting
`d.b` in `tpv104_setup` and would risk the TPV104 byte-exact contract, the **reject** path is the
faithful fix.

**Suggested fix:** add a guard near the top of `resolve_rs_impl` (it has `cfg.state_evolution`,
`cfg.depth_profile`, and `cfg.spatial`), e.g. right after the existing BoxcarTaper guard
(`spatial_friction.cpp:~1891`):
```diff
+   // Phase 11a/b: per-DOF b is wired through the AGING iterator + equilibrium
+   // seed only.  The SRW policy (state_policies.hpp) evolves psi with the
+   // SCALAR b_default, while SeedEquilibriumPsi_RS seeds with rs.b(ii) per-DOF.
+   // A per-DOF b under SRW therefore seeds the fault out of the SRW iterator's
+   // equilibrium (silent t=0 disequilibrium).  Reject the combination loudly.
+   if (cfg.state_evolution == StateEvolutionKind::SlipLawStrongRateWeakening)
+   {
+      MFEM_VERIFY(!cfg.depth_profile.enabled,
+                  "ResolveRateState: [friction.rate_state.depth_profile] (per-DOF "
+                  "b) is not supported with state_evolution="
+                  "slip_law_strong_rate_weakening; the SRW iterator uses the "
+                  "scalar b_default (Phase 11a wired per-DOF b through the aging "
+                  "iterator only). Use the aging law or a scalar b.");
+      for (const auto& r : cfg.spatial)
+      {
+         MFEM_VERIFY(std::isnan(r.b),
+                     "ResolveRateState: a per-DOF 'b' spatial override is not "
+                     "supported with state_evolution="
+                     "slip_law_strong_rate_weakening (scalar b only).");
+      }
+   }
```

**Test case:**
```cpp
// tests/unit/test_resolve_rate_state_guards.cpp (extend)
void test_R001_srw_rejects_perdof_b()
{
   // SRW + depth_profile must abort (per-DOF b seeded but evolved scalar).
   const bool aborted_profile = RunInChild([]() {
      SpatialFrictionConfig cfg = /* RS defaults, all scalars valid */;
      cfg.state_evolution = StateEvolutionKind::SlipLawStrongRateWeakening;
      cfg.depth_profile.enabled = true;
      cfg.depth_profile.profile = /* any valid 2-knot a + a-b profile */;
      cfg.V_w_default = 0.1;                 // SRW needs V_w>0
      SpatialFrictionResolver R;
      R.ResolveRateState(cfg, coords, dof_to_elem, dof_to_attr, mat, mesh,
                         PorePressureSpec{}, sigma_n_total);   // expect MFEM_ABORT
   });
   TEST_ASSERT(aborted_profile, "SRW + depth_profile aborts");

   // SRW + a per-DOF b spatial rule must also abort.
   const bool aborted_rule = RunInChild([]() {
      SpatialFrictionConfig cfg = /* RS defaults */;
      cfg.state_evolution = StateEvolutionKind::SlipLawStrongRateWeakening;
      cfg.V_w_default = 0.1;
      SpatialRule r; r.kind = SpatialRule::Kind::Box; r.b = 0.02;  // per-DOF b
      cfg.spatial.push_back(r);
      SpatialFrictionResolver R;
      R.ResolveRateState(cfg, coords, dof_to_elem, dof_to_attr, mat, mesh,
                         PorePressureSpec{}, sigma_n_total);   // expect MFEM_ABORT
   });
   TEST_ASSERT(aborted_rule, "SRW + per-DOF b spatial rule aborts");
}
```

---

### [R-002] [LOW] [spatial_friction_rate_state_safs_projected_stress_depthprofile.toml] — stale "entirely velocity-weakening / ~16.7 km" NOTE contradicts the vwvs11km CSVs actually used (~11 km transition)

**Category:** DEVIATION / QUALITY (documentation that misstates the physics)

**Description:**
The config header NOTE (lines 24–28) states: *"the (a-b)=0 VW->VS transition is at ~16.7 km, just
below the fault bottom, so the meshed fault is (barely) entirely velocity-weakening."* That describes
the **original** `param_a.csv`/`param_a_minus_b.csv`. The config actually points at
`param_a_vwvs11km.csv` / `param_a_minus_b_vwvs11km.csv` (lines 149–150), whose `(a−b)=0` transition is
**~11 km** (correctly stated in `[meta]` line 14 and the `depth_profile` comment line 147). With the
vwvs11km data the meshed fault (≈16.6 km deep) has a velocity-**strengthening** band from ≈11–16.6 km,
so the NOTE's "entirely velocity-weakening" is wrong. `--print-derived` will print the true ≈11 km,
contradicting the in-file NOTE; a user reading only the comment would misjudge the VW/VS structure.

**Trigger:** Reading the config header to understand the model (no runtime effect — computation uses
the real CSVs and is correct).

**Actual behavior:** Header NOTE claims ~16.7 km / entirely VW.

**Expected behavior:** NOTE should describe ~11 km transition and the VW(0–11 km)/VS(11–16.6 km) split
the vwvs11km CSVs actually produce. (The nucleation patch at z≈−5 km is in the VW core → super-critical;
worth stating explicitly per plan §11c req 1.)

**Suggested fix:**
```diff
- # NOTE on coverage: the 500 m mesh fault reaches ~16.6 km depth, while the CSVs
- # sample to ~52-60 km; the profile is flat-clamped below the deeper sampled
- # depth.  With the committed CSVs the (a-b)=0 VW->VS transition is at ~16.7 km,
- # just below the fault bottom, so the meshed fault is (barely) entirely
- # velocity-weakening — the --print-derived depth-profile summary reports this.
+ # NOTE on coverage: the 500 m mesh fault reaches ~16.6 km depth, while the CSVs
+ # sample to ~52-60 km; the profile is flat-clamped below the deeper sampled
+ # depth.  With the *_vwvs11km.csv profiles the (a-b)=0 VW->VS transition is at
+ # ~11 km, so the meshed fault is velocity-WEAKENING from 0-~11 km and velocity-
+ # STRENGTHENING from ~11-16.6 km.  The nucleation patch (z ~= -5 km) sits in the
+ # VW core, so it is super-critical.  The --print-derived summary reports the
+ # resolved transition depth.
```

**Test case:** (documentation — optional integration assert) `--print-derived` output line
`VW->VS transition (a-b=0) at depth = <zc> m` should report `zc ≈ 11000 m`, matching `[meta]`.

---

### [R-003] [LOW] [spatial_friction.cpp:load_one_depth_csv] — "exactly 2 fields" guard misses a trailing NON-numeric token

**Category:** EDGE_CASE

**Description:**
The row-format guard `MFEM_VERIFY(!(ss >> extra), …)` (`:199-203`) declares `real_t extra;` and so
only fires when the 3rd whitespace token parses as a number. A row such as `0.1, 50, junk` extracts
`value=0.1`, `depth_km=50`, then `ss >> extra` fails on `"junk"` (failbit) → `!(ss>>extra)` is `true`
→ no abort. The junk token is silently ignored, so the stated intent ("expected exactly 2 fields") is
only partially enforced (numeric 3rd token aborts; non-numeric 3rd token is dropped).

**Trigger:** A malformed CSV row with a non-numeric trailing token after `value, depth_km`.

**Actual behavior:** Row silently parsed as the leading `(value, depth)`; trailing junk ignored.

**Expected behavior:** Abort on ANY leftover token, numeric or not.

**Suggested fix:** read the leftover as a string (catches both kinds):
```diff
-      real_t extra;
-      MFEM_VERIFY(!(ss >> extra),
-                  "[friction.rate_state.depth_profile] " << which << " '" << path
-                  << "' line " << lineno
-                  << ": expected exactly 2 fields 'value, depth_km'");
+      std::string extra;
+      MFEM_VERIFY(!(ss >> extra),
+                  "[friction.rate_state.depth_profile] " << which << " '" << path
+                  << "' line " << lineno
+                  << ": expected exactly 2 fields 'value, depth_km' (unexpected "
+                  "trailing token '" << extra << "')");
```

**Test case:**
```cpp
void test_R003_loader_rejects_trailing_nonnumeric()
{
   const bool aborted = RunInChild([]() {
      FrictionDepthProfileSpec s;
      s.param_a_csv = write_file(tmp_path("a_r003.csv"),
                                 "0.010 0 junk\n0.030 10\n");   // 3rd token non-numeric
      s.param_a_minus_b_csv = write_file(tmp_path("amb_r003.csv"),
                                         "-0.009 0\n0.130 60\n");
      s.depth_to_m = 1000.0;
      (void)LoadFrictionDepthProfileCSVs(s);
   });
   TEST_ASSERT(aborted, "loader aborts on a trailing non-numeric token");
}
```

---

### [R-004] [LOW] [fault_face_flux.hpp:DOFData] — `b` default is `quiet_NaN()`, a documented deviation from the plan's `0.0` (§11a step 1)

**Category:** DEVIATION (informational — not a bug)

**Description:**
Plan §11a step 1 specifies `real_t b = 0.0;` with the rationale that `0.0` is the unset/LSW sentinel
(an unset RS `b` then divides by zero in `UpdateStateAnalytic`, caught loudly). The implementation uses
`std::numeric_limits<real_t>::quiet_NaN()` and adds the R-028 guard in `SeedEquilibriumPsi_RS`
(`isfinite(d.b) && d.b > 0`, `spatial_setup.hpp:433`). This is internally consistent and arguably
**stronger** than the plan (NaN propagates loudly through any arithmetic; `0.0` only mis-evolves for
`psi<f0`). Verified that no LSW / TPV104 / checkpoint path reads `d.b`, so the NaN never reaches a
consumer unset. No correctness impact; flagged only so the deviation from plan text is on record.

**Trigger:** none (no failing input).

**Actual behavior:** `DOFData.b` defaults to NaN; guarded on the only RS read paths.

**Expected behavior:** Acceptable as-is. Recommend keeping NaN and noting the deviation in the plan
(or in a one-line code comment, already partly present at `fault_face_flux.hpp:50-60`).

**Suggested fix:** none required. If strict plan-conformance is desired, change to `0.0` — but this is
NOT recommended (it weakens the unset-detection). Prefer updating the plan to specify NaN.

**Test case:** n/a (not a bug; the existing R-028 guard test in `test_seed_equilibrium_psi_rs`
already exercises the unset-`b` path).

---

### [R-005] [LOW] [spatial_dyn_driver.cpp:print-derived + spatial_print_derived.cpp] — Phase 11c summary lives in the driver and omits the sample tables + explicit min/max(a,b) the plan's Detailed Requirements list (AC still met)

**Category:** DEVIATION

**Description:**
Plan §11c "Files to Modify" places the depth-profile summary in `PrintDerivedAndCheckRS`. The
implementation instead emits it from the driver (`spatial_dyn_driver.cpp:1839-1919`) — functionally
reasonable, since the driver has `cfg.rate_state->depth_profile` while `PrintDerivedAndCheckRS`
receives only the resolved per-DOF `rs`. However the summary omits two items in the plan's Detailed
Requirements: (a) the **sample tables** `a(depth_m)` and `(a−b)(depth_m)` (knot listings); (b) the
explicit **min/max of resolved `a` and `b`** over the fault DOFs. The plan's **Acceptance Criterion #1**
(CSV paths + VW↔VS transition depth) IS satisfied, and `L_nuc`/`f_ss` ranges from
`PrintDerivedAndCheckRS` indirectly reflect the `a`/`b` ranges. Low impact; mostly an auditability gap.

**Trigger:** `--print-derived` on the depth-profile config — the summary is less complete than specced.

**Actual behavior:** Prints CSV paths, depth ranges, fault max depth, R-029 mismatch warning, and the
VW→VS transition depth.

**Expected behavior:** Also print the knot tables and the resolved `min/max(a)`/`min/max(b)` over fault
DOFs (the latter needs an Allreduce over `rs.a`/`rs.b`).

**Suggested fix (optional):** extend the rank-0 block in `spatial_dyn_driver.cpp:1854-1918` to (1) list
`dp.profile.a_of_depth.{x,y}` and `dp.profile.amb_of_depth.{x,y}` knots, and (2) MPI_Allreduce
min/max over the local `rs.a`/`rs.b` and print them. No behavior change to the simulation.

**Test case:** (integration — optional) assert the `--print-derived` stdout contains an `a min/max`
line and a knot table; not required for unit suite.

---

## Summary
- Critical issues: 1 (R-001 — Phase-11-introduced SRW + per-DOF `b` silent t=0 disequilibrium; no guard)
- Moderate issues: 0
- Low issues: 4 (R-002 stale config NOTE; R-003 CSV trailing-token; R-004 NaN-vs-0.0 deviation [no fix]; R-005 summary location/completeness)
- Plan compliance: **FULL** for 11a/11b and the 11c acceptance criteria; the CRITICAL R-001 is a gap the
  plan itself did not anticipate (it scoped SRW `b` scalar but did not add a guard against per-DOF `b`
  reaching the SRW seed). Minor documented deviations otherwise (NaN sentinel, summary location,
  vwvs11km CSV swap) that do not break any AC.
- Verdict: **FAIL — must fix before proceeding.** The shipped Phase-11c deliverable config is aging-law
  and is correct/byte-exact, but R-001 is a CRITICAL, silent, config-reachable wrong-physics regression
  on the SRW (TPV104) law and must be fixed (add the guard in `resolve_rs_impl`) before any SRW run can
  safely use a depth profile or per-DOF `b`. R-002/R-003 are quick hardening/clarity fixes; R-004 needs
  no code change; R-005 is optional polish.

## Unreviewed Areas
- Full-mesh / Frontera runtime behavior of the depth-profile config (per repo policy, no local
  full-mesh runs; verified by unit tests + reading instead).
- `io/tpv104_checkpoint.hpp` comment-only update (Phase 11a step 6, doc-only) — not re-verified line by
  line; `b` is a recomputed static field and is not serialized, consistent with the plan.
- The committed sbatch variants beyond the dev 8N/400r form (multiple normalcap/pureupwind/triq
  variants exist; only the dev 8N/400r form was inspected against the Phase 11c spec).
