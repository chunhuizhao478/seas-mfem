# Code Review — Cross-rank material exchange, Phase 5 (deprecate seam_continuous + fault-locality optional)

**Date:** 2026-06-07
**Reviewer:** chunhui-code-reviewer agent (adversarial pass)
**Plan:** `document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md` — Phase 5
**Scope (Phase 5 diff — deprecation docs + one banner string + config-parse gate):**
- `dynamic/bimaterial_wave_operator.hpp` (setter + member doc comments)
- `dynamic/bimaterial_wave_operator.inl` (two comment blocks)
- `drivers/spatial_dyn_driver.cpp` (SetSeamContinuous call comment + banner string)
- `spatial/code/spatial_friction.hpp` / `.cpp` (`MaterialSpec::seam_continuous` field + TOML parse comments)
- `tests/unit/test_tpv_config_parse.cpp` (R-208 gate, pre-existing)

**Build/run verification done in this review (worktree override, conda mfem-dev):**
- `seas_test_tpv_config_parse` rebuilt (force-removed stale `spatial_friction.o` + test obj) and run:
  **184 / 184 passed, 0 failed, exit 0.** R-208 assertions present and green
  (see "Test verification" below).
- `tests/parallel/test_bimaterial_seam_material_np2.o` recompiled from scratch ⇒ the operator
  (`bimaterial_wave_operator.hpp` + `.inl`, the comment-only Phase-5 edits) **compiles cleanly**;
  no syntax error introduced by the rewritten comment blocks.

---

## THE HEADLINE CHECK — does `seam_continuous_` have ZERO live readers?

**VERDICT: YES. The deprecation's core claim ("stored but NEVER read / no behavioral effect")
is TRUE.** Confirmed by exhaustive grep of the member identifier `seam_continuous_` (trailing
underscore = the operator member, as distinct from the config field `seam_continuous`):

```
dynamic/bimaterial_wave_operator.hpp:190:   void SetSeamContinuous(bool v) { seam_continuous_ = v; }   <- SETTER (assignment)
dynamic/bimaterial_wave_operator.hpp:375:   bool seam_continuous_ = false;                              <- MEMBER DECL
dynamic/bimaterial_wave_operator.inl:876:  // (MFEM_VERIFY mode==Constant || seam_continuous_) is REMOVED  <- COMMENT
```

Those are the ONLY three occurrences in the entire codebase (excluding `extern/`). The references
are: (1) the setter's assignment, (2) the member declaration, (3) a comment documenting that the
former guard was removed. There is:
- **No** `if (seam_continuous_)` anywhere.
- **No** `MFEM_VERIFY(... seam_continuous_ ...)` / `MFEM_ASSERT` / `MFEM_ABORT` referencing it
  — confirmed by `grep "MFEM_VERIFY|MFEM_ASSERT|MFEM_ABORT" ... | grep -i seam` returning ONLY
  the line-876 comment.
- **No** getter (`GetSeamContinuous` / `IsSeamContinuous`) — grep is empty.
- **No** read of the member inside `BuildPerFaceCentralFluxMatrices_` (the former consumer): the
  only `seam_continuous` token inside that function body is the line-876 "is REMOVED" comment.

The former consumer (the central-build abort `MFEM_VERIFY(mode==Constant || seam_continuous_)` at
the shared arm of `BuildPerFaceCentralFluxMatrices_`) is genuinely gone (Phase 2). The flag is now
a pure write-only sink. **Doc claim corroborated by code.**

---

## Findings

| ID | Severity | Category | File / line | One-liner |
|----|----------|----------|-------------|-----------|
| P5-001 | MODERATE | Doc inconsistency | `tests/parallel/test_bimaterial_mixed_flux_shared.cpp:13-38` | Test header still describes `SetSeamContinuous(true)` as a LIVE central-flux gate / abort trigger — contradicts the Phase-5 "no effect / never read" deprecation. |
| P5-002 | LOW | Doc / scope gap | `drivers/spatial_dyn_driver.cpp` (whole file) | Phase-5 deliverable "`--partition-fault-locality` optional" is moot for the spatial driver: that flag is NOT wired into `spatial_dyn_driver.cpp` at all (it exists only in the native tpv*_driver.cpp). Not a regression, but the plan line is unsatisfiable as written for the driver under review; worth an explicit note. |
| P5-003 | LOW | Dead-ish test call | `tests/parallel/test_bimaterial_mixed_flux_shared.cpp:271` | `op.SetSeamContinuous(true)` is now a no-op store; the comment `// affirm seam-continuity` overstates it. Harmless (test still 16/16). Sub-case of P5-001. |

**Critical issues: 0.**

The four Phase-5 documentation sites that were the actual diff
(`bimaterial_wave_operator.hpp` setter + member, `spatial_dyn_driver.cpp` call + banner,
`spatial_friction.hpp` + `.cpp` field + parse) are all **accurate and mutually consistent** —
every one says "DEPRECATED … no effect / stored but never read / parsed for back-compat". No site
under the diff still claims the flag gates behavior. The one inconsistency (P5-001) is in a test
file that was NOT part of the Phase-5 diff but is the natural place a reader looks to learn what the
flag does, so the stale narrative directly undercuts the deprecation message.

---

### P5-001 (MODERATE) — `test_bimaterial_mixed_flux_shared.cpp` still narrates `SetSeamContinuous` as a live gate

**File:** `tests/parallel/test_bimaterial_mixed_flux_shared.cpp`, header comment lines 13-38
(and the inline comment at 271).

**What it says now (FALSE after Phase 2/5):**
- L14-15: "`SetMixedFluxMode(Adjacent)` with a het Mode::Coefficient material +
  `SetSeamContinuous(true)` **builds** the side-0 shared central matrices" — implies the build is
  CONDITIONED on the flag.
- L24-25: "The **declarative seam-continuity gate** does NOT conflate legitimate depth variation
  with unsupported lateral variation."
- L29-38 (Test 3.4b "…RequiresSeamContinuous"): "The guard is the single MFEM_VERIFY in the shared
  arm of `BuildPerFaceCentralFluxMatrices_` (`mode==Constant || seam_continuous_`) … the negative
  of 3.4a/3.4c: **WITHOUT `SetSeamContinuous(true)` the same Mode::Coefficient shared central face
  hits the verify**."

**Why it's wrong:** The `MFEM_VERIFY(mode==Constant || seam_continuous_)` it describes was REMOVED
in Phase 2 (confirmed: `bimaterial_wave_operator.inl:875-882` is now a comment block reading
"… is REMOVED"; no live verify references the flag anywhere). The central build now proceeds for
ANY material regardless of the flag; strong-contrast safety is the contrast guard
(`mixed_flux_contrast_tol`). So:
- Test 3.4a/3.4c would build identically WITHOUT `SetSeamContinuous(true)`.
- Test 3.4b's premise (an abort if you omit the flag) is no longer true — there is no abort to hit.

**Impact:** No runtime breakage — the test passes (16/16) because the setter is a harmless store
and the build no longer gates on it. The damage is purely documentary: this is the most detailed
in-tree description of what `seam_continuous` does, and it flatly contradicts the Phase-5
"no effect / never read" claim a future maintainer would be trying to trust. Phase 5's stated job
is "deprecation messaging accurate/consistent across all sites"; this site was missed.

**Trigger:** A maintainer reads this test to understand `seam_continuous`, concludes it still gates
the central build, and either (a) keeps the dead `SetSeamContinuous(true)` calls "to avoid the
abort", or (b) re-adds a guard, contradicting the Phase-2 design.

**Suggested fix (comment-only; pick one of: rewrite the header, or delete 3.4b's stale narrative):**

```diff
--- a/tests/parallel/test_bimaterial_mixed_flux_shared.cpp
+++ b/tests/parallel/test_bimaterial_mixed_flux_shared.cpp
@@
-//   Test 3.4a (Dispatch_SharedCentralFace_Build):
-//     On an np=2 ParMesh where a fault-adjacent non-fault interior face is split
-//     across the partition seam (a SHARED central face), SetMixedFluxMode(Adjacent)
-//     with a het Mode::Coefficient material + SetSeamContinuous(true) builds the
+//   Test 3.4a (Dispatch_SharedCentralFace_Build):
+//     On an np=2 ParMesh where a fault-adjacent non-fault interior face is split
+//     across the partition seam (a SHARED central face), SetMixedFluxMode(Adjacent)
+//     with a het Mode::Coefficient material builds the
//     side-0 shared central matrices and SharedInteriorFaceFlux_ dispatches the
//     central F* ...
//
//   Test 3.4c (Dispatch_DepthProfile_NoFalseAbort):
-//     A depth-only (mu(z)) Mode::Coefficient material + SetSeamContinuous(true) on a
+//     A depth-only (mu(z)) Mode::Coefficient material on a
//     2D (x,z) grid whose partition seam carries a fault-adjacent z-normal central
-//     SHARED face ... BUILDS with no false abort.  The declarative seam-continuity
-//     gate does NOT conflate legitimate depth variation with unsupported lateral
-//     variation.
+//     SHARED face ... BUILDS via the cross-rank TRUE-peer material read (Phase 2).
+//     (Historically gated by the now-DEPRECATED seam_continuous affirmation; that
+//     guard was removed in Cross-rank Phase 2 — the build no longer consults it.)
//
-//   Test 3.4b (Dispatch_SharedCentralFace_RequiresSeamContinuous) — OMITTED from
-//     the automated assertions.  The guard is the single MFEM_VERIFY in the shared
-//     arm of BuildPerFaceCentralFluxMatrices_ (mode==Constant || seam_continuous_).
-//     ... WITHOUT SetSeamContinuous(true) the same Mode::Coefficient shared central
-//     face hits the verify ... Revisit if death-test infra lands.
+//   Test 3.4b (formerly Dispatch_SharedCentralFace_RequiresSeamContinuous) — RETIRED.
+//     It exercised the seam_continuity MFEM_VERIFY guard, which Cross-rank Phase 2
+//     REMOVED: the central build now reads the TRUE peer material and proceeds for
+//     ANY material (strong contrast handled by the contrast guard, not this flag).
+//     There is no longer an abort to negate, so 3.4b is obsolete.
```

And at line 271:

```diff
-   op.SetSeamContinuous(true);                       // affirm seam-continuity
+   op.SetSeamContinuous(true);                       // DEPRECATED no-op (Phase 5); kept for back-compat exercise
```

**Test:** No new runtime test (comment-only). A grep guard would catch regressions:
`! grep -q "RequiresSeamContinuous\|declarative seam-continuity gate" tests/parallel/test_bimaterial_mixed_flux_shared.cpp`.

---

### P5-002 (LOW) — `--partition-fault-locality` is not a spatial-driver flag; the "make it optional" deliverable is moot there

**Plan Phase 5 line:** "`--partition-fault-locality` optional; note in jobs/configs."

**Finding:** In the worktree, `--partition-fault-locality` is parsed/used ONLY in the native
drivers — `drivers/tpv102_driver.cpp`, `drivers/tpv104_driver.cpp`, `drivers/tpv205_driver.cpp`
(each: `HasFlag(..., "--partition-fault-locality")` → optional `else if` branch → already an
opt-in, not required). It is **NOT present in `drivers/spatial_dyn_driver.cpp`** (grep for
`fault.locality|FaultLocality|partition` in that file returns only two unrelated comment lines).

The driver actually exercising the cross-rank bi-material path (the spatial driver) uses the
default ParMETIS partition with no fault-locality option at all. So:
- For the **native** drivers, the flag was ALREADY optional (an `else if`, not mandatory) before
  Phase 5 — nothing to change.
- For the **spatial** driver, "make it optional" is a no-op: there is no flag to make optional.

**Impact:** None functionally. This is a plan-vs-worktree scope note: the Phase-5 deliverable as
phrased ("make fault-locality optional") has no actionable target in the driver under review,
because (a) the spatial driver never required it, and (b) the only drivers that expose the flag
already treat it as optional. Worth recording so the deliverable isn't mistaken for incomplete.

**Cross-rank-does-NOT-require-fault-locality — POSITIVELY CONFIRMED (the whole point of the
feature):** The np=2 gate tests deliberately partition so the seam CUTS the fault, with NO
fault-locality, and pass:
- `tests/parallel/test_bimaterial_seam_fault_np2.cpp` — `PartitionByYSign(...)` partitions
  elements by centroid y-sign "so every y=0 fault face straddles the seam" (comment L156-157);
  reported 24/24.
- `tests/parallel/test_bimaterial_seam_material_np2.cpp` — EXPLICIT partition at `seam_x=2` so the
  chosen interior face straddles the rank seam (L764); reported 135/135.

These pass with the fault explicitly cut across ranks and no fault-locality applied, which is the
empirical proof that correctness does not depend on fault-locality. Borne out.

**Suggested fix:** Add one sentence to the Phase-5 section / jobs note:
> "`--partition-fault-locality` is an opt-in flag of the NATIVE tpv*_driver.cpp only (already an
> `else if`, never required); the spatial_dyn_driver.cpp does not expose it and never required it —
> the cross-rank exchange is correct under the default fault-cutting ParMETIS partition (gated by
> seam_fault_np2 / seam_material_np2)."

---

### P5-003 (LOW) — dead `SetSeamContinuous(true)` call + overstated comment in the np=2 mixed-flux test

`tests/parallel/test_bimaterial_mixed_flux_shared.cpp:271` `op.SetSeamContinuous(true); // affirm
seam-continuity`. Post-Phase-2 this is a no-op store; the comment overstates it. Sub-case of
P5-001; folded into that fix's diff. No separate action required.

---

## Back-compat (R-208) — VERIFIED at source AND by running the gate

**Source level:** `spatial/code/spatial_friction.cpp:1650`
`cfg.material.seam_continuous = toml_bool(m, "seam_continuous", false);`. The helper
`toml_bool` (spatial_friction.cpp:304-311):
```
bool toml_bool(... key, bool default_val) {
   if (!tbl.contains(key)) { return default_val; }   // key ABSENT -> default false, NO abort
   const auto& v = tbl.at(key);
   if (v.is_boolean()) { return v.as_boolean(); }     // =true / =false -> the value, NO abort
   MFEM_ABORT(... "must be a boolean");                // only aborts on a NON-boolean type
}
```
So a config OMITTING the key parses (default false), and a config with `seam_continuous = true`
parses (round-trips, no abort). Back-compat is structurally safe.

**Gate run (this review, fresh build):** `seas_test_tpv_config_parse` → **184 / 184, exit 0.**
The R-208 assertions are present and green:
```
PASSED: tpv31 material.seam_continuous defaults false (key omitted)
PASSED: 5.3a: matrix + mixed_flux=adjacent parses WITHOUT G1 abort (BUG-21)
PASSED: 5.3b: [material].seam_continuous=true round-trips into MaterialSpec (BUG-22)
PASSED: tpv31_rk_mixedflux [material].seam_continuous == true
PASSED: tpv31_p2_rk_mixedflux [material].seam_continuous == true
PASSED: tpv31 p-variant seam_continuous defaults false (key omitted)
```

---

## Banner string + compile check (no syntax error / unterminated literal)

`drivers/spatial_dyn_driver.cpp:1332-1337`:
```cpp
std::cout << "[mixed-flux] matrix (bi-material) + "
          << cfg.numerics.mixed_flux << " central flux enabled "
          << "(seam_continuous=" << (cfg.material.seam_continuous
                                     ? "true" : "false")
          << "; DEPRECATED, no effect — cross-rank seam material uses the "
             "TRUE peer)\n";
```
Well-formed: the ternary is fully parenthesized, every string literal is terminated, and the final
"…uses the " / "TRUE peer)\n" pair is valid adjacent-string-literal concatenation. No unterminated
literal, no stray operator. (Full driver TU not linked here — heavyweight — but the operator TU
`bimaterial_wave_operator.{hpp,inl}` that the comment-only Phase-5 edits touched WAS recompiled
clean via `test_bimaterial_seam_material_np2.o`, so the header/inl comment rewrites are
syntactically sound.)

---

## Deprecation-message consistency matrix (the 4 diff sites)

| Site | File:line | Says it has effect? | Wording |
|------|-----------|---------------------|---------|
| Setter | `bimaterial_wave_operator.hpp:181-190` | NO | "DEPRECATED … NO BEHAVIORAL EFFECT … stored but never read" |
| Member | `bimaterial_wave_operator.hpp:371-375` | NO | "DEPRECATED … stored but NEVER READ … config back-compat" |
| Driver call | `spatial_dyn_driver.cpp:1205-1211` | NO | "DEPRECATED … now has NO effect … flag is stored, never read" |
| Driver banner | `spatial_dyn_driver.cpp:1334-1337` | NO | "DEPRECATED, no effect — cross-rank seam material uses the TRUE peer" |
| Config field | `spatial_friction.hpp:603-611` | NO | "DEPRECATED … parsed for back-compat but has NO effect … no longer consulted" |
| Config parse | `spatial_friction.cpp:1647-1650` | NO | "DEPRECATED … parsed for back-compat but no longer has any effect" |
| .inl neighbour | `bimaterial_wave_operator.inl:288-304` | n/a (history) | accurately describes the stub→true-peer change |
| .inl central | `bimaterial_wave_operator.inl:875-882` | NO | "former … guard … is REMOVED" |

**All eight under-diff / adjacent sites are consistent and correct.** The ONE inconsistent site is
the test header in `test_bimaterial_mixed_flux_shared.cpp` (P5-001), which was outside the Phase-5
diff but still claims the flag gates the build.

---

## Regression risk from the comment-only edits

None expected and none found. The operator TU recompiled clean (no behavior change — comments
only). The config TU recompiled and its 184 tests pass. Reported regression suites
(seam_material np2 135/135, seam_fault np2 24/24, perside 11/11, parity 9/9, central 6/6,
dispatch 10/10, constant-parity 19/19, mixed-flux-shared np2 16/16) are consistent with
comment-only Phase-5 changes plus the already-landed Phases 1-4.

---

## Summary

**Counts:** 3 findings total — 0 CRITICAL, 1 MODERATE (P5-001), 2 LOW (P5-002, P5-003).

**Critical issues: 0.**

**The single most important check passed:** `seam_continuous_` (the operator member) has exactly
three references — setter assignment, member declaration, and one "is REMOVED" comment — and
**ZERO live readers**. No `if`, no `MFEM_VERIFY`, no getter consults it. The deprecation's core
claim ("stored but never read; no behavioral effect") is TRUE and corroborated by code.

- **Back-compat (R-208): SOLID.** `toml_bool(..., false)` never aborts on absent-or-boolean;
  `seam_continuous` omitted → false, `=true` → round-trips. Gate run: 184/184, exit 0.
- **Banner string: well-formed**, no unterminated literal; operator + config TUs recompile clean.
- **Cross-rank does NOT require fault-locality: confirmed** — np=2 gates cut the fault across ranks
  (PartitionByYSign / explicit seam_x partition) with no fault-locality and pass.
- **The 4 (8 incl. .inl) deprecation sites in the diff are accurate and mutually consistent.**

**The one substantive defect** is documentary: `test_bimaterial_mixed_flux_shared.cpp` still
narrates `SetSeamContinuous(true)` as a live central-flux gate / abort trigger (Tests 3.4a/3.4c
header + the retired-3.4b paragraph), directly contradicting the Phase-5 "no effect" message. The
test passes (the setter is a harmless store), so this is MODERATE not CRITICAL — but Phase 5's job
was message consistency, and this is the most detailed in-tree description of the flag, so it should
be fixed (comment-only diff provided in P5-001).

**Verdict: APPROVE (Phase 5 is sound). Merge-ready.** No code change required for correctness.
Recommend the comment-only P5-001 fix (and the P5-002 plan/jobs note) before closing the
deprecation, so the in-tree narrative stops contradicting the very deprecation Phase 5 lands.

### Recommended next steps
1. Apply the P5-001 comment-only diff (retire the 3.4b "RequiresSeamContinuous" narrative; drop the
   "declarative seam-continuity gate" language) so the test no longer claims the flag gates the build.
2. Add the P5-002 one-line note clarifying `--partition-fault-locality` is a native-driver opt-in
   only and the spatial driver never required fault-locality.
3. (Optional, LOW) Annotate the dead `SetSeamContinuous(true)` call at
   test_bimaterial_mixed_flux_shared.cpp:271 as a deprecated no-op (folded into step 1).
