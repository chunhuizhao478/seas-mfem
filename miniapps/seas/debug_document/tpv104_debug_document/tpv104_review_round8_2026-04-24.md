# Code Review Round 8: post-R7 honest-banner fix verification — 2026-04-24

Fresh adversarial audit after the round-7 R7-001 option-(b) fix
landed (banner downgrade + dispatch-tag disclosure). Primary focus:
(a) confirm R7-001..R7-007 closures are load-bearing and don't mask
new bugs, (b) hunt for regressions introduced by the banner rewrite
and the `--verify-dispatch` flag.

## Review Scope

- Plan: `tpv104_debug_plan_2026-04-24.md` §4.10 Step 9 + §4.10.X.
- Prior reviews: R-001..R-012, R2-001..R2-007, R3-001..R3-007,
  R4-001..R4-008, R5-001..R5-007, R6-001..R6-004, R7-001..R7-007.
- Files changed since round 7:
  - `drivers/tpv104_driver.cpp` (689 → 750, +61 lines).
  - `tests/unit/test_tpv104_smoke.cpp` (202 → 320, +118 lines).
  - `jobs/tpv104/*.sbatch` — PASS-criteria updates.
- Unchanged since round 7: iterator (hpp/cpp), friction_solver,
  fault_face_flux, wave_operator, all other test files.

## Prior-round closure verification

| ID | Prior severity | Round 8 status |
|---|---|---|
| R7-001 option (b) | CRITICAL | **CLOSED** — `drivers/tpv104_driver.cpp:266-287` banner now describes actual dispatch truthfully. `--verify-dispatch` flag (R7-004) prints machine-readable `[dispatch]` lines. |
| R7-002 | MODERATE | **CLOSED** — `drivers/tpv104_driver.cpp:523-531` removes the dead `SlipLawSRWPsi law` construction. Docstring points to where to reintroduce it under R7-001 option (a). |
| R7-003 | MODERATE | **CLOSED** — `jobs/tpv104/tpv104_init_stations.sbatch:17-29` and `tpv104_t3_mesh_coupled.sbatch:24-51` contain honest R7-001 disclosure blocks. `tpv104_build_banner_sanity.sbatch` PASS criteria include `[dispatch] friction_solver_actual=brent` check (`:87-93`). |
| R7-004 | MODERATE | **CLOSED** — `TestDispatchMatchesBanner` at `test_tpv104_smoke.cpp:196-243` verifies banner-vs-dispatch tri-consistency across 5 CLI variants. |
| R7-005 (= R6-001/R6-002) | MODERATE | **CLOSED** — `MapSolver` `(void)`-called at `:310` with `MFEM_ABORT` on unknown strings (`:123-125`). `newton` aliased to `newton-stable` (`:116-119`). `TestUnknownSolverAborts` at `:248-273` validates. |
| R7-006 | LOW | **CLOSED** — driver comment at `:226-237` accurately disclosures the CLI-vs-runtime mismatch. |
| R7-007 | LOW | **DOCUMENTED** — driver comment at `:641-650` and sbatch headers explicitly acknowledge the plan §3.12 deviation. Not code-fixable without R7-001 option (a). |

All round-7 findings CLOSED or properly documented. Frontera-sanity
path is now internally consistent.

## Findings

### [R8-001] MODERATE [drivers/tpv104_driver.cpp:146-152, tests/unit/test_tpv104_smoke.cpp:222-243] — `ActualDispatchSolverTag()` / `ActualDispatchIteratorTag()` / `ActualDispatchLawTag()` return hard-coded constants; tri-consistency test asserts those same constants — the test cannot detect a R7-001-option-(a) wiring that forgets to update the Tag functions

**Category:** BUG (test design — false-negative once the iterator wiring lands)

**Description:**
```cpp
// drivers/tpv104_driver.cpp:146-152
static std::string ActualDispatchSolverBanner()
{
   return "Brent (hard-coded via EvaluateADERTotal; --friction-solver flag IGNORED)";
}
static const char *ActualDispatchSolverTag() { return "brent"; }
static const char *ActualDispatchIteratorTag() { return "oneshot"; }
static const char *ActualDispatchLawTag() { return "slip-srw"; }
```
These functions return hard-coded strings that describe the CURRENT
state of the driver (one-shot + Brent). The R7-004 tri-consistency
test asserts the printed banner and `[dispatch]` lines match these
hard-coded values for every CLI variant.

**Problem:** when R7-001 option (a) lands (iterator wiring), the
dispatched method becomes CLI-dependent:
- `--friction-solver newton-stable` → `NewtonRaphsonStable` dispatched.
- `--friction-solver brent` → `Brent` dispatched.
- `--fault-iterator substep` → `Tpv104SubStepIterator::Advance` wired in.

At that point the Tag functions MUST change to return runtime values
derived from the actual dispatch, not hard-coded constants. If a
maintainer forgets to update them:
- `ActualDispatchSolverTag()` still returns `"brent"`.
- Banner line (line 275) still says `"Friction solver: Brent
  (hard-coded…)"`.
- Meanwhile the time loop dispatches `NewtonRaphsonStable`.

The tri-consistency test would pass (banner + dispatch tag agree on
"brent") even though the actual solver is NewtonRaphsonStable — the
exact disclosure failure R7-001 was filed to prevent, re-introduced
by the R7-004 fix's hard-coding.

**Trigger:** Future R7-001 option (a) wiring lands without updating
the `ActualDispatch*Tag` / `ActualDispatch*Banner` functions.

**Actual behavior:** Banner says "Brent + one-shot" is dispatched;
dispatch lines match; but actual runtime is something else.

**Expected behavior:** `ActualDispatch*Tag` values should be derived
from a runtime-computed enum of the dispatched path — not hard-coded.

**Suggested fix:** Introduce a `DispatchedPath` enum driven by a
single source of truth:
```diff
+enum class DispatchedPath {
+   OneShotBrent,          // Current R7-001 option (b).
+   SubStepIterator,       // Future R7-001 option (a).
+};
+
+// One source of truth for what the driver actually dispatches.
+// UPDATE THIS when you change the time loop, NOT the Tag functions.
+static constexpr DispatchedPath kDispatch = DispatchedPath::OneShotBrent;

-static std::string ActualDispatchSolverBanner()
-{
-   return "Brent (hard-coded via EvaluateADERTotal; --friction-solver flag IGNORED)";
-}
-static const char *ActualDispatchSolverTag() { return "brent"; }
-static const char *ActualDispatchIteratorTag() { return "oneshot"; }
-static const char *ActualDispatchLawTag() { return "slip-srw"; }
+static std::string ActualDispatchSolverBanner(const std::string &cli)
+{
+   switch (kDispatch)
+   {
+   case DispatchedPath::OneShotBrent:
+      return "Brent (hard-coded via EvaluateADERTotal; "
+             "--friction-solver flag IGNORED)";
+   case DispatchedPath::SubStepIterator:
+      return SolverBanner(cli) + " (via Tpv104SubStepIterator)";
+   }
+   return "unknown";
+}
+static const char *ActualDispatchSolverTag(const std::string &cli)
+{
+   switch (kDispatch)
+   {
+   case DispatchedPath::OneShotBrent:    return "brent";
+   case DispatchedPath::SubStepIterator:
+      if (cli == "newton-stable" || cli == "newton") return "newton-stable";
+      if (cli == "brent")                             return "brent";
+      if (cli == "newton-legacy")                     return "newton-legacy";
+      if (cli == "hybrid")                            return "hybrid";
+      return "unknown";
+   }
+   return "unknown";
+}
+static const char *ActualDispatchIteratorTag()
+{
+   return (kDispatch == DispatchedPath::OneShotBrent) ? "oneshot" : "substep";
+}
+static const char *ActualDispatchLawTag() { return "slip-srw"; }
```
Update `TestDispatchMatchesBanner` to derive the expected tag from the
CLI value + known `kDispatch`, rather than asserting a hard-coded
constant.

**Test case:**
```cpp
void test_R8_001_tri_consistency_survives_future_wiring() {
   // Compile with kDispatch = SubStepIterator and assert the banner
   // changes per --friction-solver value AND the dispatch tag matches.
   // With hard-coded tags, this fails.  With runtime-derived tags
   // (R8-001 fix), it passes.
   //
   // Today this is unreachable because kDispatch = OneShotBrent, but
   // the infrastructure must be ready before R7-001 option (a) lands.
}
```

---

### [R8-002] LOW [drivers/tpv104_driver.cpp:310] — `(void)MapSolver(friction_solver);` computes the Method then discards it; future iterator wiring needs to remove the `(void)` cast AND thread the Method down to the time loop, two edits that could fall out of sync

**Category:** QUALITY (fragile pattern)

**Description:**
```cpp
// R7-005: validate --friction-solver eagerly so typos abort before
// any simulation work (including --dry-run), even though the return
// value is not routed through the (one-shot) time loop — see R7-001.
(void)MapSolver(friction_solver);
```
The `(void)` cast is there for side-effect — `MapSolver` aborts on
unknown strings via `MFEM_ABORT`. The return `FrictionSolver::Method`
value is computed and discarded. When R7-001 option (a) wires the
iterator, someone must:
1. Replace `(void)MapSolver(...)` with `const auto method = MapSolver(...)`.
2. Pass `method` to the iterator's `Advance` call.

Both edits must land together. If a future maintainer only does step
(1) and forgets step (2), the wiring is inert — same class of bug as
R7-001.

**Trigger:** R7-001 option (a) landing with only partial wiring.

**Actual behavior:** `method` is unused silently today; will stay
unused if step (2) is missed.

**Expected behavior:** The `(void)` cast should be removed now and
replaced with a named local, so forgetting to thread it downstream is
a compiler warning (unused variable) rather than silent.

**Suggested fix:**
```diff
-   // R7-005: validate --friction-solver eagerly so typos abort before
-   // any simulation work (including --dry-run), even though the return
-   // value is not routed through the (one-shot) time loop — see R7-001.
-   (void)MapSolver(friction_solver);
+   // R7-005: validate --friction-solver eagerly so typos abort before
+   // any simulation work (including --dry-run).  Store the result as
+   // `method` even though the (one-shot) time loop does not consume
+   // it today — see R7-001.  When the iterator is wired (R7-001
+   // option a), `method` drives the per-sub-step friction solve.
+   const FrictionSolver::Method method = MapSolver(friction_solver);
+   (void)method;   // unused until iterator wiring lands
```
This is almost-equivalent (both work) but signals to the next
maintainer that `method` is the handoff point.

---

### [R8-003] LOW [dynamic/tpv104_setup.hpp:449-518] — `TPV104SurfaceStationWriter` has no `Close()` method; driver calls `Flush()` only, relying on destructor for file-close ordering that runs AFTER `MPI_Finalize()`

**Category:** QUALITY (asymmetric API + subtle MPI-ordering risk)

**Description:**
```cpp
// dynamic/tpv104_setup.hpp
class TPV104StationWriter {
   ...
   void Flush();
   void Close();   // explicit close
};

class TPV104SurfaceStationWriter {
   ...
   void Flush();   // only flush; no Close
};
```
`TPV104StationWriter::Close()` is an explicit API; the symmetrically-
named surface writer class has no Close(). In the driver:
```cpp
// Post-loop cleanup (line 730-733):
station_writer.Flush();
station_writer.Close();
surface_writer.Flush();                 // no Close call
Tpv104SubStepIterator::CloseAllProbeFiles();
...
MPI_Finalize();                         // line 747
return 0;                               // destructors fire here
```
`surface_writer`'s files close during its destructor when `main()`
returns. Since `MPI_Finalize()` was called at line 747 BEFORE the
destructor runs, the ofstream flush-and-close happens post-finalize.
For plain `std::ofstream` to a POSIX path this is fine — POSIX write
doesn't need MPI. But:
1. Asymmetry between the two writer classes is a code-smell.
2. If a future change routes surface output through MPI-IO (e.g.,
   per-rank files with rank-aware filename), the post-finalize
   destructor flush would be UB.

**Trigger:** Any future refactor that makes surface-writer output
MPI-aware; OR a filesystem under SLURM job-step teardown.

**Actual behavior:** Works today on plain `ofstream`; subtly risky
under MPI-IO.

**Expected behavior:** Add `Close()` to match `TPV104StationWriter`
and call it in the driver before `MPI_Finalize()`.

**Suggested fix:**
```diff
    void Flush()
    {
       for (auto &f : files_) { if (f.is_open()) { f.flush(); } }
    }
+
+   void Close()
+   {
+      for (auto &f : files_) { if (f.is_open()) { f.close(); } }
+   }

    private:
```
Driver update:
```diff
    surface_writer.Flush();
+   surface_writer.Close();
    Tpv104SubStepIterator::CloseAllProbeFiles();
```

---

### [R8-004] LOW [drivers/tpv104_driver.cpp:288-298] — `--verify-dispatch` `[dispatch]` lines print only on rank 0; multi-rank run inspections of non-rank-0 logs will miss the tri-consistency disclosure

**Category:** QUALITY (rank-0 disclosure limitation)

**Description:**
```cpp
if (rank == 0)
{
   ...
   if (verify_dispatch)
   {
      std::cout << "[dispatch] friction_solver_actual="
                << ActualDispatchSolverTag() << "\n";
      ...
   }
}
```
The `[dispatch]` lines are printed inside `if (rank == 0)`. Under a
multi-rank MPI run (e.g., the production Frontera `sbatch` with 400
ranks), only rank 0's stdout contains the dispatch disclosure. Other
ranks' logs don't show it. If a Frontera operator inspects a specific
rank's stderr/stdout file to debug a NaN on that rank, the
tri-consistency line is missing.

Not a bug — stdout is generally rank-0-only and this is the same
pattern TPV102 uses for banners. But when `SEAS_DIAG_TPV104_STATE` is
on AND the probe file for a non-rank-0 reveals an anomaly, the
operator has no way to confirm that rank's dispatched solver matches
the banner without cross-referencing rank 0's log.

**Trigger:** Multi-rank MPI run + per-rank stderr inspection.

**Actual behavior:** `[dispatch]` lines absent from non-rank-0 logs.

**Expected behavior:** Either (a) document that `[dispatch]` is
rank-0-only, OR (b) print per-rank `[dispatch]` lines gated by
`--verify-dispatch`.

**Suggested fix (option a — document):**
```diff
    if (verify_dispatch)
    {
-      // Machine-readable tri-consistency lines (R7-004).
+      // Machine-readable tri-consistency lines (R7-004).
+      // RANK-0-ONLY: per the `if (rank == 0)` guard above, these
+      // lines appear only in rank 0's stdout.  Non-rank-0 logs
+      // cannot self-verify dispatched solver under multi-rank runs;
+      // cross-reference rank 0's log instead.
       std::cout << "[dispatch] friction_solver_actual="
```
Or (option b — per-rank, richer output): move the `[dispatch]` block
outside the rank-0 guard and prepend rank index to each line.

---

## Summary

- Critical issues: 0
- Moderate issues: 1 (R8-001 — future-proofing the tri-consistency
  infrastructure against R7-001 option (a) wiring)
- Low issues: 3 (R8-002 `(void)` cast pattern, R8-003 surface-writer
  missing Close, R8-004 rank-0-only dispatch disclosure)
- Plan compliance: FULL (under R7-001 option (b)'s banner-downgrade
  semantics). Plan §4.10 Step 7 iterator remains UNWIRED by
  conscious choice; disclosure is honest in the banner + sbatch
  headers + smoke tests.
- Verdict: PASS WITH FIXES — none of the findings block Frontera
  sbatch submission today. R8-001 should close before R7-001 option
  (a) is pursued (prevents re-introducing the banner-vs-code
  divergence under the new wiring). R8-002..R8-004 are polish.

## Frontera readiness — round 8 call

| Script | Submit today? | Gate status |
|---|---|---|
| `tpv104_build_banner_sanity.sbatch` | **GO** | banner grep matches shipped strings; `[dispatch]` lines present |
| `tpv104_mesh_build.sbatch`           | **GO** | no driver dependency |
| `tpv104_init_stations.sbatch`        | **GO** (new this round) | driver's mesh-coupled branch works; tfinal=0 writes station t=0 row and exits; P2_D gate checks 9 fault stations |
| `tpv104_t3_mesh_coupled.sbatch`      | **GO** (new this round) | blocks on init-PASS (enforced at `:76-82`); production path runs Brent + one-shot honestly |

Round-7's blocked Scripts 3 + 4 are now unblocked by the R7-001
option (b) honest-banner fix. **All four scripts are submittable in
order**, with your explicit approval per
`feedback_frontera_approval.md`. Expected sequence with SU cost:

1. `tpv104_build_banner_sanity.sbatch` — ~0.04 SU, ~2 min.
2. `tpv104_mesh_build.sbatch` — ~0.08 SU, ~5 min.
3. `tpv104_init_stations.sbatch` — ~0.5 SU, ~5 min.
4. `tpv104_t3_mesh_coupled.sbatch` — ~2000 SU, ~5 h.

**Caveat**: the production run's "stable-asinh Newton" narrative is
an explicitly abandoned aspiration under R7-001 option (b). The
Phase-3 probe-diff analysis against the reference FVW runtime must
account for the fact that MFEM runs Brent + one-shot + macro-step
ψ, whereas the reference runs its own stable-asinh Newton + per-sub-
step accumulation. Probe thresholds (§5.10.2 Probe 2 tolerance, §5.10.4
Probe 4 tolerance) must be calibrated AT LEAST an order of magnitude
looser than originally planned. The plan §5.10 probe thresholds
assumed the iterator was wired — that assumption is no longer valid.

## Unreviewed Areas

- WaveOperator::AdvanceADER behavior on ParMesh — not reproducible
  locally. First Frontera run will surface any issue.
- Shared-fault MPI path in the production time loop — exercised by
  TPV102 but TPV104 has never run multi-rank.
- The fix report's claim of "263/263 tests passing" — not
  independently verified in this round; trust the implementer's
  report.

---

## Validation checklist — Round 8

### Required closures before R7-001 option (a) is pursued (future work)

- [ ] **R8-001**: `ActualDispatch*Tag` / `ActualDispatch*Banner`
  functions derive from a single `kDispatch` enum + CLI-aware switch;
  hard-coded constants replaced with runtime logic. Tri-consistency
  test updated to use the runtime logic.
- [ ] **R8-002**: `(void)MapSolver(...)` replaced with
  `const auto method = MapSolver(...); (void)method;` (or threaded
  into iterator call). Compiler warning catches forgotten wiring.

### Recommended polish

- [ ] **R8-003**: `TPV104SurfaceStationWriter::Close()` added;
  driver calls it before `MPI_Finalize()`.
- [ ] **R8-004**: `[dispatch]` lines documented as rank-0-only OR
  moved outside the rank-0 guard.

### Rolled forward from prior rounds

- [ ] R4-007 station-writer MPI tie-break tolerance — STILL OPEN.
- [ ] R5-007 Σ deltaT tolerance O ≤ 1000 guard — STILL OPEN.
- [ ] Step 13 probe-diff tool normalises reference → MFEM
  normal-stress sign flip (closes T_TPV104_SIGN_4).
- [ ] Phase 3 probe threshold recalibration — **now MORE important**
  per the round-8 caveat (iterator unwired → looser thresholds
  required).

### Cross-round invariant anchors (updated)

- σ_n > 0 compression convention — T_TPV104_SIGN_1.
- `SlipLawSRWPsi` production-mode guards — R-001 + T_SRW_7 (live in
  unit tests; NOT active in driver — driver uses free function).
- Raw V through state-evolution — R-002.
- Unrolled `(V/V_w)^8` — R-004.
- Newton input validation — R2-004.
- Nucleation accumulator telescopes — T_TPV104_NUC_2/5.
- TPV104 strike-slip invariant — R3-001.
- TPV104 mesh byte-identical to TPV102 — T_TPV104_MESH_0/0b.
- Iterator slip accumulation signbit + tight window — R5-002.
- Iterator dt-linearity test no longer self-diff — R5-001.
- Plan §4.10 Step 5 reachable via `FrictionSolver::Method::
  NewtonRaphsonStable` enum dispatch — R5-003 CLOSED at enum layer;
  **driver runs Brent regardless per R7-001 option (b)**.
- Banner-vs-code tri-consistency — R7-001 CLOSED; banner honest.
  **R8-001 open for future R7-001 option (a) work.**
- Iterator is reachable only from unit tests — not the driver
  production path. Plan §4.10 Step 7 deliverable is IN-REPO but
  UNWIRED, per R7-001 option (b).
