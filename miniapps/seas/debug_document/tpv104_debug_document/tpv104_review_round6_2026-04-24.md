# Code Review Round 6: R5-003 fix verification + pre-Frontera final check — 2026-04-24

Final adversarial audit after the R5-003 fix lands. Focus: confirm
`FrictionSolver::Method::NewtonRaphsonStable` wiring, CLI / banner /
smoke-test consistency, and Frontera scaffolding readiness.

## Review Scope

- Plan: `tpv104_debug_plan_2026-04-24.md` §4.10.X (R5-003 fix plan).
- Prior reviews (all earlier findings closed or carried forward):
  R-001..R-012, R2-001..R2-007, R3-001..R3-007, R4-001..R4-008,
  R5-001..R5-007.
- **New / changed since round 5:**
  - `dynamic/friction_solver.hpp` (enum + `SolveNRStable` decl)
  - `dynamic/friction_solver.cpp` (dispatch + `SolveNRStable` impl)
  - `dynamic/tpv104_substep_iterator.hpp` (docstring + 165 → 179 lines)
  - `dynamic/tpv104_substep_iterator.cpp` (457 → 496 lines)
  - `drivers/tpv104_driver.cpp` (R5-003 CLI + banner + `MapSolver`)
  - `tests/unit/test_friction_solver_stable.cpp` (NEW, 5 tests)
  - `tests/unit/test_tpv104_smoke.cpp` (banner expectation updated)
  - `jobs/tpv104/*` — 3 sbatch + 1 README (NEW, Frontera scaffolding)
- Pre-existing SeisSol-name scrub: one comment reference in
  `friction_solver.hpp:37` ("SeisSol SlowVelocityWeakeningLaw.h") —
  BP5/TPV102 shared code, untouched per [C2]. Acceptable.

## R5-003 fix plan compliance

| Plan item | Delivered | Status |
|---|---|---|
| `FrictionSolver::Method::NewtonRaphsonStable` enum value | `friction_solver.hpp:51-57` | ✅ |
| `FrictionSolver::SolveNRStable` dispatch method | `friction_solver.cpp:75-92` | ✅ |
| Brent warm-start → `SolveSlipRateNewtonStable` | `friction_solver.cpp:82-88` | ✅ |
| Iterator default = `NewtonRaphsonStable` | `tpv104_substep_iterator.hpp` (see R6-note-1) | ⚠ (see finding below) |
| Driver CLI `--friction-solver={brent,newton,newton-stable,hybrid}` | `tpv104_driver.cpp:119, 232-254` | ⚠ (see R6-001) |
| Default CLI value `newton-stable` | `tpv104_driver.cpp:119` | ✅ |
| Banner discloses exact code path | `tpv104_driver.cpp:175-207` | ✅ |
| `test_friction_solver_stable.cpp` T_FS_STABLE_1..5 | 5 tests, 16/16 passing | ✅ |
| Smoke test `TestBannerDefaults` updated to match banner | `test_tpv104_smoke.cpp:90-91` | ✅ |

**Net:** R5-003 is effectively closed on the enum + dispatch + test
layer. Two behavioral inconsistencies between the plan text and the
shipped CLI — flagged below.

## Findings

### [R6-001] MODERATE [drivers/tpv104_driver.cpp:242-247] — CLI value `--friction-solver=newton` aliases to LEGACY (MFEM-native μ), not the plan-canonical stable-asinh Newton

**Category:** DEVIATION (plan text vs implementation)

**Description:**
Plan §4.10.X's R5-003 fix plan lists the CLI options as
`{brent, newton, newton-stable, hybrid}` with `newton-stable` as the
default. The actual `MapSolver` implementation:
```cpp
if (s == "newton-legacy" || s == "newton")
{
   // Bare "newton" is an alias for newton-legacy — existing scripts
   // continue to work; banner disclosure keeps the distinction.
   return FrictionSolver::Method::NewtonRaphson;
}
```
Bare `--friction-solver=newton` resolves to **legacy** (MFEM-native μ
with the `psi/a > 700` branch), not the plan-canonical stable-asinh
Newton (`NewtonRaphsonStable`). A user who reads the plan text
"Newton solver as default" and types `--friction-solver newton`
expecting the canonical path gets the legacy path silently — the
banner does disclose "legacy, MFEM-native μ" but the CLI shortcut
points in the opposite direction from intent.

The round-5 R5-003 closure criteria asked for an unambiguous mapping
between CLI choice and dispatched solver. Currently:

| CLI arg | Expected (plan text) | Actual (MapSolver) |
|---|---|---|
| `newton-stable` | stable-asinh (canonical) | stable-asinh ✅ |
| `newton` | stable-asinh (canonical) | **legacy MFEM μ** ❌ |
| `newton-legacy` | legacy MFEM μ | legacy MFEM μ ✅ |
| `brent` | Brent | Brent ✅ |
| `hybrid` | Hybrid | Hybrid ✅ |

The intent for the alias was "existing scripts continue to work" —
but no shipped script today uses `--friction-solver=newton` (only
`newton-stable` does). The alias is gratuitous and surprises users.

**Trigger:** Any user typing `--friction-solver newton` expecting the
TPV104 canonical path.

**Actual behavior:** Legacy Newton runs. Banner prints "legacy,
MFEM-native μ" — the banner is correct but the CLI alias is
misleading.

**Expected behavior:** Either (a) `newton` → `NewtonRaphsonStable`
(the canonical shorthand), OR (b) remove the alias entirely and
require explicit `newton-stable` / `newton-legacy`.

**Suggested fix (option a — canonical shorthand):**
```diff
 static FrictionSolver::Method MapSolver(const std::string &s)
 {
-   if (s == "newton-stable")
+   if (s == "newton-stable" || s == "newton")
    {
+      // Bare "newton" is an alias for the plan-canonical
+      // `newton-stable` (§4.10.X) — the SeisSol-aligned
+      // stable-asinh Newton.  Typing `--friction-solver=newton`
+      // therefore routes to the TPV104-canonical path (not legacy).
       return FrictionSolver::Method::NewtonRaphsonStable;
    }
    if (s == "brent")
    {
       return FrictionSolver::Method::Brent;
    }
-   if (s == "newton-legacy" || s == "newton")
+   if (s == "newton-legacy")
    {
-      // Bare "newton" is an alias for newton-legacy — existing scripts
-      // continue to work; banner disclosure keeps the distinction.
       return FrictionSolver::Method::NewtonRaphson;
    }
```
Update the banner branch mirror (line 175-197) to match: `if
(friction_solver == "newton-stable" || friction_solver == "newton")`.

**Test case:**
```cpp
void test_R6_001_cli_newton_routes_to_stable() {
   const std::string binary = LocateDriver();
   if (binary.empty()) { return; }
   const std::string out = RunDriver(binary,
                                     "--dry-run --friction-solver newton");
   // After R6-001: `--friction-solver newton` routes to stable-asinh.
   EXPECT_NE(out.find("Friction solver: Newton-Raphson (stable-asinh"),
             std::string::npos);
   EXPECT_EQ(out.find("Friction solver: Newton-Raphson (legacy"),
             std::string::npos);
}
```

---

### [R6-002] MODERATE [drivers/tpv104_driver.cpp:203-207, 252-253] — Unknown `--friction-solver` string silently dispatches to stable-asinh while banner prints `"(unknown)"` — banner disagrees with what actually ran

**Category:** BUG (silent misconfiguration)

**Description:**
Banner branch (`tpv104_driver.cpp:203-207`):
```cpp
else
{
   os << "Friction solver: " << opts.friction_solver
      << " (unknown)\n";
}
```
`MapSolver` fallback (line 252-253):
```cpp
// Default (including unknown strings): plan-compliant stable-asinh.
return FrictionSolver::Method::NewtonRaphsonStable;
```
If a user makes a typo — `--friction-solver newton-stble` — the
banner prints `"Friction solver: newton-stble (unknown)"` but the
solver quietly runs `NewtonRaphsonStable`. The user's log says one
thing; the physics run did another. This is precisely the disclosure
failure R5-003 was meant to prevent.

Worse: the smoke test's `TestBannerDefaults` forbidden-on-default
list (`test_tpv104_smoke.cpp:109-115`) does NOT include
`"(unknown)"`, so a typo in any future sbatch script would:
1. Still run `NewtonRaphsonStable` (appears correct physically).
2. Print `"(unknown)"` in the banner.
3. Pass the banner test (required strings are present; `(unknown)` is
   not on the forbidden list).
4. Quietly mislead the downstream probe-diff analyst about what
   solver was actually used.

**Trigger:** Any typo or unrecognised value passed to
`--friction-solver`.

**Actual behavior:** Banner reports "(unknown)"; solver silently
runs stable-asinh; no error emitted.

**Expected behavior:** Either (a) abort with a loud usage message on
unknown strings, OR (b) banner reflects the actual dispatched method
(reading from the `Method` enum after `MapSolver`, not from the raw
CLI string).

**Suggested fix (option a — abort on unknown; recommended):**
```diff
 static FrictionSolver::Method MapSolver(const std::string &s)
 {
-   if (s == "newton-stable")
+   if (s == "newton-stable" || s == "newton")
    { return FrictionSolver::Method::NewtonRaphsonStable; }
    if (s == "brent")         { return FrictionSolver::Method::Brent; }
-   if (s == "newton-legacy" || s == "newton")
+   if (s == "newton-legacy")
    { return FrictionSolver::Method::NewtonRaphson; }
    if (s == "hybrid")        { return FrictionSolver::Method::HybridNRBisection; }
-   // Default (including unknown strings): plan-compliant stable-asinh.
-   return FrictionSolver::Method::NewtonRaphsonStable;
+   // R6-002: unknown strings abort loudly to keep the banner and the
+   // dispatched solver consistent.  No silent fallback.
+   throw std::runtime_error(
+      "Unknown --friction-solver value: '" + s +
+      "'.  Valid choices: brent | newton | newton-stable | "
+      "newton-legacy | hybrid.");
 }
```
Mirror the same early-out in the banner branch:
```diff
    else
    {
-      os << "Friction solver: " << opts.friction_solver
-         << " (unknown)\n";
+      // Should be unreachable after R6-002 — MapSolver throws on
+      // unknown strings before the banner is printed.  Kept as a
+      // defence-in-depth diagnostic.
+      os << "Friction solver: " << opts.friction_solver
+         << " (UNKNOWN — driver should have aborted earlier!)\n";
    }
```

**Test case:**
```cpp
void test_R6_002_unknown_solver_aborts() {
   const std::string binary = LocateDriver();
   if (binary.empty()) { return; }
   const std::string out = RunDriver(
      binary, "--dry-run --friction-solver newton-stble");   // typo
   // After R6-002: driver aborts with a loud error message.
   EXPECT_NE(out.find("Unknown --friction-solver value"),
             std::string::npos);
   // Banner must NOT advertise a silent default on a typo.
   EXPECT_EQ(out.find("Friction solver: Newton-Raphson (stable-asinh"),
             std::string::npos);
}
```

---

### [R6-003] LOW [dynamic/friction_solver.cpp:82-88] — `SolveNRStable` does Brent warm-start first, adding ~10× friction-solver cost vs native Newton with per-QP V_prev

**Category:** QUALITY (performance, documented tradeoff)

**Description:**
`SolveNRStable` warm-starts from a full Brent solve:
```cpp
const real_t V_warm = SolveBrent(tau, psi, sigma_n, eta, a);
int iterations = 0; bool converged = false;
const real_t V = SolveSlipRateNewtonStable(
   tau, psi, std::abs(sigma_n), eta, a, /*V0=*/V0,
   /*V_prev=*/std::max<real_t>(V_warm, static_cast<real_t>(0)),
   /*max_iter=*/60, /*tol=*/1e-8, &iterations, &converged);
return converged ? V : V_warm;
```
Brent on TPV104's log10-V space converges in ~10 iterations; the
stable-asinh Newton from V_warm converges in ~1 iteration. The
combined cost is ~11× a pure Newton from per-QP V_prev.

For the target Phase-2 run: `t = 3.0 s × 640 macro-steps × 5 sub-steps
× ~1M fault QPs` = ~3.2 billion friction-solver calls per rank-group,
across 400 ranks. At ~40 %  of stage time in friction solves (rough
estimate from TPV102 profiling), the 11× overhead translates to a
wall-clock penalty of ~4× on the friction-solver portion. Target
wall-clock 5 h → measured ~6 h likely.

Already flagged in comment `friction_solver.cpp:71-74`:
> "`FrictionSolver::Solve` has no V_prev parameter, so we warm-start
> from Brent's root… A future Step-9 optimisation may extend Solve's
> signature to accept V_prev from DOFData.slip_rate."

Not a bug — deliberate tradeoff. Flagged LOW so the Step-9 wiring
priority list includes "extend FrictionSolver::Solve with optional
V_prev override sourced from DOFData.slip_rate".

**Trigger:** Production t=3.0 s run at 400 ranks.

**Actual behavior:** ~11× friction-solver cost vs native-Newton-with-
V_prev.

**Expected behavior:** After the Step-9 extension, pure Newton with
per-QP V_prev; native TPV104 performance.

**Suggested fix (Step-9 follow-up):**
```diff
 real_t Solve(real_t tau, real_t psi, real_t sigma_n,
              real_t eta, real_t a, Method method = Method::Brent,
-             ) const;
+             real_t V_prev_override = -1.0) const;
+   // V_prev_override >= 0 seeds Newton-family methods directly and
+   // skips the Brent warm-start.  Pass
+   // DOFData.slip_rate (from the previous macro-step) for optimal
+   // convergence on TPV104.  -1.0 retains the current Brent-warm-
+   // start behaviour for backward compatibility.
```

**No test case required** (LOW performance concern, not a correctness
bug).

---

### [R6-004] LOW [test_tpv104_smoke.cpp:139-157] — `TestBannerNonDefault` covers `newton-legacy`, `brent`, `hybrid`, but NOT bare `newton` — if R6-001's alias is fixed, no test validates the new mapping

**Category:** QUALITY (test coverage gap)

**Description:**
`TestBannerNonDefault` at `test_tpv104_smoke.cpp:139-157` exercises:
- `--friction-solver newton-legacy` → banner `"Friction solver:
  Newton-Raphson (legacy"`.
- `--friction-solver brent` → banner `"Friction solver: Brent
  (log10-V"`.

Neither `newton-stable` (default-banner tested separately) nor bare
`newton` appears. If R6-001's fix lands (bare `newton` →
stable-asinh), no test would verify the re-routing. Similarly, no
test verifies `--friction-solver hybrid` produces `"Friction solver:
Hybrid"`.

**Trigger:** Any refactor of `MapSolver`'s string matching.

**Actual behavior:** 2 of 4 non-default CLI values are banner-tested.

**Expected behavior:** All 4 non-default CLI values banner-tested,
OR a parameterised test that runs the full `{brent, newton,
newton-stable, newton-legacy, hybrid}` matrix.

**Suggested fix:** add `TestBannerAllSolverOptions`:
```diff
+void TestBannerAllSolverOptions()
+{
+   const std::string binary = LocateDriver();
+   if (binary.empty()) { return; }
+
+   const std::vector<std::pair<std::string, std::string>> table = {
+      {"brent",         "Friction solver: Brent (log10-V"},
+      {"newton",        "Friction solver: Newton-Raphson (stable-asinh"},
+      {"newton-stable", "Friction solver: Newton-Raphson (stable-asinh"},
+      {"newton-legacy", "Friction solver: Newton-Raphson (legacy"},
+      {"hybrid",        "Friction solver: Hybrid NR+Bisection"},
+   };
+   for (const auto &[arg, expected] : table)
+   {
+      const std::string out = RunDriver(
+         binary, "--dry-run --friction-solver " + arg);
+      TEST_ASSERT(out.find(expected) != std::string::npos,
+                  ("--friction-solver=" + arg + " banner").c_str());
+   }
+}
```

---

## Frontera scaffolding — readiness audit

| Script | Ready today? | Notes |
|---|---|---|
| `tpv104_build_banner_sanity.sbatch` | **YES** | Runs `--dry-run`; does not need mesh. Safe to submit. |
| `tpv104_mesh_build.sbatch` | **YES** | Runs gmsh on the `.geo` templates. No driver dependency. |
| `tpv104_init_stations.sbatch` | **NO** | Driver's `--mesh` branch is not yet wired to the iterator; `--tfinal 0.0` hits `RunDryRun` regardless of `--mesh`. Header comment at `:30-38` correctly documents this as blocked. Submitting today will fail the post-run PASS gate. |
| `tpv104_t3_mesh_coupled.sbatch` | **NO** | Same driver-extension blocker + gated on Script 3 PASS (enforced at `:76-82`). Safe to have in the repo; submitting today would fail the grep-gate immediately (exit 1 before burning SUs). |

The 3 → 4 exit-1 guard is the most important defensive check — it
prevents a ~2000 SU burn if someone submits script 4 before script 3
passes. Correctly implemented at `tpv104_t3_mesh_coupled.sbatch:76-82`.

## Summary

- Critical issues: 0
- Moderate issues: 2 (R6-001 CLI alias direction, R6-002 unknown-
  string silent fallback)
- Low issues: 2 (R6-003 Brent-warm-start overhead, R6-004 test
  coverage gap)
- Plan compliance: FULL for R5-003 enum + dispatch + test layer.
  The CLI alias contract (plan §4.10.X listed `newton` as canonical)
  is PARTIAL — see R6-001.
- Verdict: PASS WITH FIXES — R6-001 and R6-002 should be closed
  before Step 9 ships. R6-003/R6-004 are polish; do not block the
  first two Frontera sbatch scripts (build-sanity + mesh-build), both
  of which are fit to submit today.

## Unreviewed Areas

- Step-9 driver mesh-coupled branch (the remaining blocker for
  Scripts 3 + 4). Correctly flagged as BLOCKED-PENDING-INTEGRATION in
  both sbatch headers; no code exists to review yet.
- `Tpv104SubStepIterator::Advance` behaviour under 400-rank MPI with
  shared fault faces — not reproducible locally per
  `feedback_no_local_reproducer`. First Frontera run will surface any
  issue.
- Wall-clock estimate for `tpv104_t3_mesh_coupled.sbatch` (header
  claims ~5 h; R6-003 suggests closer to 6 h due to Brent-warm-start
  overhead). Will be measured on the first real run.

---

## Validation checklist — Round 6

### Required closures before Step 9 driver lands

- [ ] **R6-001**: bare `--friction-solver=newton` routes to
  `NewtonRaphsonStable` (plan-canonical). Update `MapSolver` string
  match and the matching banner branch.
- [ ] **R6-002**: unknown `--friction-solver` string aborts with a
  loud error message (no silent default to stable-asinh).

### Recommended closures before Phase 3.B Frontera probe runs

- [ ] **R6-003**: Step 9's driver extension adds per-QP `V_prev`
  override threaded into `FrictionSolver::Solve`; skips the
  Brent warm-start inside `SolveNRStable` when `V_prev_override >=
  0`. Gives ~10× speedup on the friction-solver portion of the
  production run.
- [ ] **R6-004**: `TestBannerAllSolverOptions` covers all five CLI
  strings (including bare `newton` after R6-001 lands).

### Rolled forward from prior rounds

- [ ] R4-007 — station-writer MPI tie-break tolerance + no-owner
  abort. Still OPEN.
- [ ] R5-001 — `TestConvergenceUnderDtHalving` rewrite (not a no-op).
  Still OPEN.
- [ ] R5-002 — `TestSlipAccumulation` signbit directional guard. Still
  OPEN.
- [ ] Step-9 driver mesh-coupled branch (Scripts 3 + 4 prerequisite).
- [ ] Step 13 probe-diff tool normalises the reference → MFEM
  normal-stress sign flip.

### Frontera submission go/no-go

**GO today (no R6-001/R6-002 dependency):**
- Script 1 `tpv104_build_banner_sanity.sbatch` — banner + dry-run on
  Frontera HEAD. 5 min / ~0.04 SU.
- Script 2 `tpv104_mesh_build.sbatch` — gmsh mesh build (200 m /
  500 m / 1000 m). 30 min / ~0.08 SU.

**NO-GO today:**
- Script 3 `tpv104_init_stations.sbatch` — needs Step 9 mesh-coupled
  extension.
- Script 4 `tpv104_t3_mesh_coupled.sbatch` — needs Script 3 PASS
  (enforced by guard).

### Cross-round invariant anchors (updated)

- σ_n > 0 = compression — `T_TPV104_SIGN_1`.
- Production `SlipLawSRWPsi` aborts on base-virtual — R-001 + T_SRW_7.
- Raw V through state-evolution — R-002 + TestNoVsafeClamp.
- `(V/V_w)^8` unrolled integer power — R-004 + R2-002.
- Newton input validation — R2-004.
- Nucleation accumulator telescopes — T_TPV104_NUC_2/5.
- TPV104 strike-slip invariant (tau1_nuc = sigma_n_nuc = 0) — R3-001.
- TPV102 overwrite-pattern nucleation unaffected — T_TPV104_NUC_4.
- TPV104 mesh byte-identical to TPV102 — T_TPV104_MESH_0/0b.
- Slip accumulation live in iterator — R4-001 (sign-blind per R5-002).
- Probe-3 μ via `FrictionCoefficientStable` — R4-003.
- Probe output per-rank under MPI — R4-005.
- **Plan §4.10 Step 5 canonical solver reachable** — R5-003 enum +
  dispatch CLOSED; CLI alias contract PARTIAL — **NEEDS R6-001 CLOSURE**.
- **CLI / banner / solver tri-consistency** — banner matches solver
  on valid strings ✅; silent default on invalid strings — **NEEDS
  R6-002 CLOSURE**.
