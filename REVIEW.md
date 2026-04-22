# Code Review: TPV102 v9.3.0 — final audit + ADER+total-Q sbatch prep, 2026-04-22

## Review Scope

**User request:**
1. Final review on all code changes.
2. Prepare two sbatch job files similar to
   `jobs/tpv102/tpv102_200m_p1_0.01s_50rank_init.sbatch` and
   `jobs/tpv102/tpv102_200m_p1_2.0s_400rank_vmax_halt.sbatch`, but
   both on dev queue, 8 nodes / 400 cores, 2-hour budget, using ADER
   and total-Q formulation.

**Purpose-change check (carry-forward from round 6/7):**
- ADER is the default time integrator; RK4 is an alternative.
- The total-Q formulation is the only production path.

- Plans:
  - `miniapps/seas/debug_document/tpv102_debug_document/tpv102_ader_time_integration_plan.md`
  - `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v9.3.0_debug_plan.md`
- Files touched this round:
  - `drivers/tpv102_driver.cpp` (line 138: `--order` default flipped
    from 2 → 1; round-7 R-001 fix landed.)
  - `jobs/tpv102/tpv102_200m_p1_0.01s_400rank_ader_total_init.sbatch`
    (NEW — mirrors 50-rank init sanity, scaled to 8 nodes / 400 cores).
  - `jobs/tpv102/tpv102_200m_p1_2.0s_400rank_ader_total_vmax_halt.sbatch`
    (NEW — mirrors 400-rank V_max halt test, with ADER + total-Q
    banner assertions).
- Files reviewed:
  - `drivers/tpv102_driver.cpp` (CLI defaults, init path, nucleation,
    RK4 + ADER branches)
  - `dynamic/wave_operator.hpp` / `.inl`
  - `dynamic/fault_face_flux.hpp` / `.cpp`
  - `dynamic/godunov_flux.hpp` / `.cpp`
  - `dynamic/tpv102_setup_total.hpp`
  - Reference sbatch files (`tpv102_200m_p1_0.01s_50rank_init.sbatch`,
    `tpv102_200m_p1_2.0s_400rank_vmax_halt.sbatch`) for module / env
    matching.
- Domain context:
  - `miniapps/seas/CLAUDE.md`
  - `CLAUDE.md` (project)
  - Auto-memory: "Ask before Frontera runs — notify user and wait for
    explicit approval before any Frontera sbatch, ssh, or allocation
    check" — sbatch files CREATED but NOT submitted; Frontera approval
    gate remains on the user.

## Purpose-change compliance audit (final)

| Purpose | Status | Evidence |
|---|---|---|
| (1) ADER default, RK4 alternative | **REALIZED** | `tpv102_driver.cpp:155` — `--time-integrator` defaults to `"ader"`; unknown-flag fallback also targets `"ader"` and re-sets `use_ader = true` (round-6 R-005 fix). Banner at line ~1046 says "ADER-O(N) (default)". `tests/unit/test_R001_driver_defaults_to_ader.cpp` gates this. |
| (2) Total-Q only in the production driver | **REALIZED** | `tpv102_driver.cpp:543` calls `ZeroDOFDataPreStressTotal`; 507-512 `wave.SetAbsorbingBackground(Q_bg)`; 589 `InitializeStateTotal`; 602 `MFEM_VERIFY(wave.GetAbsorbingBackground() != nullptr)` fails loud if total-Q setup is skipped. All five in-loop nucleation sites use `ApplyNucleationTotal`. `tests/unit/test_R002_driver_init_total_q.cpp` gates this. |
| (3) Consistent ADER order vs DG order | **REALIZED (this round)** | `tpv102_driver.cpp:138` — `--order` default flipped from 2 → 1 so `--ader-order=2` default satisfies `O >= p + 1`. |

All three purpose-change items now land cleanly.

## Findings

### [R-001] LOW [dynamic/wave_operator.inl:928-1005, 1148-1159, 1641-1652, 1838-1911, 1988-2001, 2318-2331] — fluctuation-Q `else` branches still reachable at the wave-operator level; driver-level `MFEM_VERIFY` catches production misuse, but the branches remain live for BP5 / fluctuation-mode unit tests

**Category:** QUALITY / DEVIATION (carried forward from round 7 R-002)

**Description:**
Every total-Q-aware dispatch site in `wave_operator.inl` is
`if (has_bulk_bg_) { *Total(...) } else { *(...) }`. The driver's
`MFEM_VERIFY(wave.GetAbsorbingBackground() != nullptr)` at line 602
prevents the `else` branch from being reached under
`seas_tpv102_driver`. But:
- The `else` paths are still live in `wave_operator.inl` and used by
  BP5 + fluctuation-only unit tests (`test_R001_rk4_fluctuation_fault_dispatch`,
  `test_interior_fault_flux_path`, etc.).
- Any future driver that forgets `SetAbsorbingBackground` silently
  falls through to fluctuation semantics.

This was flagged in round 7 as MODERATE and carried forward to this
round at LOW because the user's purpose #2 is realized at the DRIVER
level; removing the wave-operator `else` branches would also require
migrating BP5 and several unit tests, which is out of scope for this
final review.

**Trigger:**
A future new driver (non-TPV102) that forgets `SetAbsorbingBackground`,
OR a refactor that drops the driver-level `MFEM_VERIFY` at line 602.

**Actual behavior:**
Wave operator silently takes fluctuation path.

**Expected behavior (deferred):**
Convert the `else` branches to `MFEM_VERIFY(has_bulk_bg_, "total-Q
only")` once BP5 and fluctuation-only tests are migrated.

**Suggested fix:** See round-7 REVIEW R-002 suggested fix (patch
listed 8 dispatch sites). Defer until BP5 migration.

**Test case:** Same as round-7 R-002.

---

### [R-002] LOW [drivers/tpv102_driver.cpp:1095] — `Vector Q_new(Q.Size())` allocated inside the per-step ADER loop; pre-allocate once before the loop for 10^5-step production runs

**Category:** QUALITY (micro-performance, carried from round 7 R-004)

**Description:**
ADER branch at line 1095 re-allocates a `9 · ndof_total`-size vector
every step. At 10^6 DOFs × 10^5 steps, that's 10^5 heap allocations of
~72 MB each. Declaring `Q_new` once before the loop (like `Q_tmp` at
line 995) amortizes the allocation.

**Suggested fix:**
```diff
@@ -994,6 +994,7 @@
    Vector k1(Q.Size()), k2(Q.Size()), k3(Q.Size()), k4(Q.Size());
    Vector Q_tmp(Q.Size());
+   Vector Q_new(Q.Size());   // reused by the ADER branch every step

    real_t t = 0.0;
    real_t V_max_global = 0.0;
@@ -1093,7 +1094,6 @@
-         Vector Q_new(Q.Size());
          wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
-         Q = Q_new;
+         Q.Swap(Q_new);   // O(1); Q_new is overwritten next iteration
```
Combines R-004 + R-005 from round 7. Defer — harmless for this
2-hour dev run but should land before longer production runs.

---

## sbatch files prepared

Two files created under `miniapps/seas/jobs/tpv102/`:

### 1. `tpv102_200m_p1_0.01s_400rank_ader_total_init.sbatch` (init sanity)

- **Queue:** `development`
- **Nodes / Ranks:** 8 / 400
- **Wall-clock budget:** 02:00:00 (actual usage ~5–10 min; the 2 hr
  ceiling is the user-requested allocation — `--tfinal=0.01` caps
  the run well below).
- **Integrator / formulation:** `--time-integrator ader
  --ader-order 2 --order 1` (explicit belt-and-braces even though
  these are now the driver defaults after round-6/7 fixes).
- **Acceptance criteria (automated in RESULT.txt):**
  - `[BUILD]` banner present.
  - `State representation:` banner contains `total-Q`.
  - `Time integrator:` banner contains `ADER-O`.
  - `[R-101 check]` line shows `pairs matched > 0` (or the documented
    "SKIPPED" message).
  - no `MFEM_ABORT` in the job log.

### 2. `tpv102_200m_p1_2.0s_400rank_ader_total_vmax_halt.sbatch` (V_max halt test)

- **Queue:** `development`
- **Nodes / Ranks:** 8 / 400
- **Wall-clock budget:** 02:00:00
- **Integrator / formulation:** `--time-integrator ader
  --ader-order 2 --order 1`.
- **Acceptance criteria (automated in RESULT.txt):**
  - Banner asserts `total-Q` + `ADER-O` modes (hard-failure if either
    string is missing — catches future driver-default regressions).
  - 9-station slip-propagation score (mirrors the RK4+fluctuation
    reference sbatch).
  - `V_max` first / peak / final values extracted for trend inspection.
  - Decision table:
    - 0/9 slipped + V_max peak ~ V_ini ⇒ "investigate (possible v9.3
      regression)".
    - ≥6/9 slipped + V_max peak > 0.1 m/s ⇒ "rupture evolved; v9.3
      ADER+total-Q validated".
    - Intermediate ⇒ "partial rupture; investigate per-station".

Both sbatches mirror the module load + LD_LIBRARY_PATH of the working
reference sbatches (per auto-memory rule
`feedback_sbatch_modules.md`). The rebuild block rm's
`drivers/tpv102_driver.o` and all `dynamic/*.o` to guarantee the
round-6/7 fixes compile in (important: the p=1 default flip applied
this round must be visible at run time).

**Frontera approval gate honoured (per auto-memory
`feedback_frontera_approval.md`):** files CREATED on disk; NOT
submitted. User submits manually with `sbatch <path>` after review.

## Summary
- Critical issues: 0
- Moderate issues: 0
- Low issues: 2 (R-001 fluctuation `else` branches deferred to BP5
  migration; R-002 ADER `Q_new` allocation pattern)
- Plan compliance: **FULL** for all three purpose-change items
  (ADER default, total-Q only, O >= p+1 consistent order).
- **Verdict: PASS.**

The code is ready for the two sbatch submissions. The remaining LOW
items do not affect the runs as-prepared: R-001 is a defensive-coding
hardening for non-TPV102 callers, and R-002 is a micro-performance
optimization.

## Next actions for the user

1. Commit the `--order` default flip (line 138 of
   `drivers/tpv102_driver.cpp`) and the two new sbatch files.
2. `git pull` on Frontera's `/scratch2/10024/zhaochun/seas-project/seas-mfem`.
3. Submit the init sanity sbatch first:
   ```
   sbatch miniapps/seas/jobs/tpv102/tpv102_200m_p1_0.01s_400rank_ader_total_init.sbatch
   ```
   Expect `PASS` in `RESULT.txt` within ~10 min.
4. On PASS, submit the 2 hr V_max halt test:
   ```
   sbatch miniapps/seas/jobs/tpv102/tpv102_200m_p1_2.0s_400rank_ader_total_vmax_halt.sbatch
   ```

## Unreviewed Areas

1. **End-to-end numerical comparison of ADER+total-Q vs the v8.0.0
   RK4+fluctuation baseline.** The station-slip score in the V_max
   halt sbatch is a pass/fail gate, not a correctness comparator.
   A cross-integrator diff of `flt_0_7.5.dat` between the v8.0.0
   reference and the v9.3 ADER run would identify any silent
   numerical drift.
2. **MFEM Vector::Swap availability** on the Frontera build's MFEM
   version. The R-002 fix sketch assumes `Vector::Swap(Vector &)` is
   available; verified in the local MFEM checkout at
   `/Users/chunhuizhao/projects/mfem/linalg/vector.hpp:379`. Should
   be present on Frontera's MFEM too (standard MFEM 4.x API), but
   deferred until the code-fix pass lands.
3. **BP5 regression under the `--order` default flip.** BP5's driver
   probably passes `--order` explicitly; verify before flipping the
   default affects any BP5 sbatch.
