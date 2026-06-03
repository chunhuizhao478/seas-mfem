# Code Review: Phase 0 — fault-dealiasing local reproduction harness (2026-06-02)

## Review Scope
- Plan: `miniapps/seas/document/fault_dealiasing_dev/seisol_overintegration_resample_speckle_2026-06-02.md`, **§6 Phase 0** ("Local reproduction harness (no production runs)").
- Files reviewed (the Phase-0 worktree changes on `feat/fault-overint-resample`):
  - `miniapps/seas/tests/unit/test_fault_planar_serial.cpp` (new, the harness)
  - `miniapps/seas/Makefile` (new build/test targets)
- Cross-referenced for API/convention fidelity:
  - `tests/unit/test_rupture_multistep_serial_vs_parallel.cpp` (the sibling test this is patterned on)
  - `dynamic/d4_tet_mesh.hpp` (`BuildD4Mesh` — the Kuhn 6-tet split + orient flag lifted here)
  - `dynamic/wave_operator.hpp` + `.inl` (ctor, `AdvanceADER`, fault dispatch at `:4117-4145`, `ComputeMaxDt` at `:5840-5911`, `SetAbsorbingBackground`)
  - `dynamic/fault_face_flux.hpp` (`DOFData::sigma_n_corr/slip_rate`), `dynamic/tpv102_setup.hpp` + `tpv102_setup_total.hpp` (`InitializeFaultDOFs`, `ZeroDOFDataPreStressTotal`, `InitializeStateTotal`), `config/tpv102_params.hpp`, `domain/boundary_config.hpp`
- Domain context consulted: `CLAUDE.md`, `miniapps/seas/CLAUDE.md` (sign conventions, total-Q dispatch, gmsh tet meshes), the plan's companion `REVIEW.md` (plan review, round 2).
- Build status: **the test compiles and links** — `seas_test_fault_planar_serial` (10 MB binary) and `tests/unit/test_fault_planar_serial.o` are present and current (built 15:59 vs source 15:59). No API/signature/link defects exist. The harness was **not executed** during this review (project policy: no local TPV102 reproducer runs — `[[feedback_no_local_reproducer]]`).

**Verification of the core metric (not a finding — confirms the design is sound):** the harness uses the total-Q setup (`InitializeStateTotal` puts the pre-stress in bulk `Q`, `ZeroDOFDataPreStressTotal` zeroes `DOFData.sigma_n0/tau*_0`). `AdvanceADER` with the default `FaultFrictionLaw::RateAndState` dispatches the interior-fault branch to `EvaluateADER` (`wave_operator.inl:4141`), which on total-`Q` writes `DOFData.sigma_n_corr` as the **TOTAL** normal traction (≈120 MPa). So `|sigma_n_corr − TPV102Params::sigma_n|` does measure the σ_n excursion as the header claims. `slip_rate` is likewise written. Both metrics are correctly sourced.

---

## Findings

### [R-001] MODERATE [POSSIBLE] [test_fault_planar_serial.cpp:RunRupture / Options::dt] — `dt` is hardcoded and decoupled from `--he` and the material wave speed: shrinking the element with `--he` (or `--nx`-with-finer-mesh intent) silently violates the ADER CFL limit and there is no guard

**Category:** ASSUMPTION / EDGE_CASE

**Description:**
`Options::dt = 1.0e-4` is a fixed constant (line 106) and is passed straight into `wave.AdvanceADER(Q, o.dt, 2, Q_new)` (line 523). `AdvanceADER` performs **no** internal CFL check — it integrates whatever `dt` it is handed. The CLI exposes `--he` (element edge) and `--dt` as **independent** knobs, and the header only documents that the default `dt` is "ADER-2-stable for he ~ 1 km" — an unenforced coupling.

The ADER-2 explicit stability bound scales as `dt ≲ cfl · h / cp`. With `cp = 6000 m/s` and the default `he = 1000 m`, `dt = 1e-4` sits in the stable band (the sibling `test_rupture_multistep_serial_vs_parallel` validates exactly this `dt`/`he`/`cp` combination on ~1 km tets). But the moment a user passes e.g. `--he 250`, the stable `dt` drops ~4×, the run becomes ADER-unstable, and the solution blows up. A blow-up makes `excursion_peak` huge (or `NaN`), so the acceptance assertions (`> 1e5`, `final > 3×early`) can **pass for the wrong reason** — the harness then "reproduces a growing σ_n excursion" that is pure numerical instability, not friction aliasing. That directly defeats the purpose of a Phase-0 A/B testbed.

The wave operator already exposes the right tool: `WaveOperator::ComputeMaxDt(cfl)` returns `cfl · h_min_ / cp` (`wave_operator.inl:5910`). The harness ignores it.

**Trigger:** `./seas_test_fault_planar_serial --he 250` (or any `--he` materially below ~700, or `--dt` raised above the stable band). Default invocation is safe.

**Actual behavior:** silent ADER blow-up; assertions may still report PASS off an instability-driven excursion, or print `NaN`/`inf` with no diagnostic.

**Expected behavior:** the harness should derive a CFL-safe `dt` from the operator (or hard-guard the supplied `dt`) so a mis-sized mesh fails loud instead of producing a misleading "reproduction."

**Suggested fix:** guard `dt` against the operator's CFL limit right after the operator is built in `RunRupture` (uses the already-constructed `wave`):
```diff
   WaveOperator<MeshT> wave(mesh, o.order, TPV102Params::lambda,
                            TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);
+
+  // CFL guard: AdvanceADER does no internal stability check, and --he/--dt
+  // are independent knobs.  ComputeMaxDt(cfl) = cfl * h_min / cp.  Use a
+  // conservative ADER-2 cfl=0.5; abort loud rather than blow up silently
+  // and report an instability-driven "excursion" as a reproduction.
+  const real_t dt_cfl = wave.ComputeMaxDt(/*cfl=*/0.5);
+  MFEM_VERIFY(o.dt <= dt_cfl,
+              "dt=" << o.dt << " s exceeds the ADER-2 CFL limit "
+              << dt_cfl << " s for he=" << o.he << " m (cp="
+              << TPV102Params::cp << " m/s).  Pass a smaller --dt or a "
+              "larger --he; a too-large dt makes the run blow up and the "
+              "sigma_n acceptance pass off a numerical instability, not "
+              "the friction aliasing this harness must isolate.");
```

**Test case:**
```python
def test_R001_cfl_guard_rejects_undersized_mesh():
    # Default mesh is stable; shrinking he without shrinking dt must NOT
    # silently "reproduce" a drift off an instability.
    ok = run_harness(args=[])                      # he=1000, dt=1e-4
    assert ok.returncode == 0
    bad = run_harness(args=["--he", "250"])        # dt now ~4x over CFL
    # With the guard: fail-loud (MFEM_VERIFY abort, nonzero exit, CFL message).
    assert bad.returncode != 0
    assert "exceeds the ADER-2 CFL limit" in bad.stderr
    # Without the guard (regression we are preventing): the run blows up and
    # may print NaN/inf excursions or PASS off the instability.
```

---

### [R-002] LOW [test_fault_planar_serial.cpp:PrintUsage / Options::nsteps] — documented default `--nsteps` (15000) contradicts the code default (8000); the no-arg acceptance path (`make test-fault-planar-serial`) runs 8000

**Category:** QUALITY (doc/code drift with an acceptance consequence)

**Description:**
`Options::nsteps = 8000` (line 108), but `PrintUsage` advertises `"--nsteps N number of ADER-2 steps (default 15000)"` (line 125). The Makefile target `test-fault-planar-serial` runs the binary with **no arguments** (Makefile, new target), so the acceptance path uses 8000 steps = 0.8 s of simulated time — not the 15000 (1.5 s) the docs imply. If the secular drift needs more than 0.8 s to clear the `final > 3×early` bar, the default acceptance run fails while a reader following the documented default (15000) would pass. At minimum the two numbers must agree; the choice should be the one that reliably reproduces the drift on the default mesh.

**Trigger:** `make test-fault-planar-serial` (uses defaults) vs reading `--help`.

**Actual behavior:** acceptance runs 8000 steps; help says 15000.

**Expected behavior:** one consistent default, chosen to satisfy the §6 Phase-0 acceptance on the default `8×8×2` mesh.

**Suggested fix (align the doc to the code, or bump the code if 8000 under-runs):**
```diff
-      "  --nsteps N       number of ADER-2 steps (default 15000)\n"
+      "  --nsteps N       number of ADER-2 steps (default 8000)\n"
```
The file-header usage example (line 51, `--nsteps 15000`) should be reconciled to the same number. If a smoke confirms 8000 is too few to reach `final > 3×early`, prefer raising `Options::nsteps` to the value that reproduces and update the help string to match.

**Test case:**
```python
def test_R002_documented_default_matches_code_default():
    help_txt = run_harness(args=["--help"]).stdout
    m = re.search(r"--nsteps N\s+number of ADER-2 steps \(default (\d+)\)", help_txt)
    assert m and int(m.group(1)) == NSTEPS_DEFAULT_IN_OPTIONS  # 8000
```

---

### [R-003] LOW [POSSIBLE] [test_fault_planar_serial.cpp:BuildFaultSlab] — the slab is only one element thick normal to the fault, with the absorbing boundary one element from the fault, which may suppress or contaminate the [[v_n]]-driven secular σ_n drift the harness must reproduce

**Category:** ASSUMPTION (geometry adequacy for the stated purpose)

**Description:**
`BuildFaultSlab` builds `nyv = 3` y-levels → exactly **one element layer on each side** of the fault (`Ly = 2*he`), and tags every exterior face — including the two faces only `he` away from the fault — as absorbing (attr 5). The drift mechanism the plan targets (§3.4) is: friction-aliased checkerboard `τ_corr` → radiated `v_n` into the near-fault bulk → `[[v_n]] ≠ 0` at the fault → `δσ_n = η_p·[[v_n]]`. With a single near-fault element and an absorber pressed right against it, the near-fault wavefield is heavily shaped by the absorbing BC (normal-radiated P-energy is absorbed within `he/cp ≈ 1.7e-4 s ≈ 2` steps). Plan §8 explicitly warns that near-fault bulk dissipation is a *distinct* speckle source from the friction aliasing — a 1-element-thick slab maximizes exactly that confound. It is uncertain (cannot be confirmed without running) whether this geometry reproduces the friction-aliasing drift faithfully, reproduces a boundary-driven artifact instead, or suppresses the drift below the `1e5 Pa` / `3×` bars.

This is flagged POSSIBLE: the in-plane extent (8×8 → 128 fault faces) is ample, the per-step speckle and the secular slip/state accumulation are largely fault-local, and the asymmetric main mesh is the correct choice to seed the leak — so it may well reproduce. But the normal-direction thickness is the weakest link in the harness's claim to isolate friction aliasing.

**Trigger:** the default acceptance run; sensitivity is structural, not input-dependent.

**Suggested fix:** parameterize the half-thickness so the absorber is several elements from the fault, e.g. add `int ny_half = 2;` to `Options` and generalize the y tiling (`nyv = 2*ny_half + 1`, fault plane at `y = ny_half*he`, hex layers `iy ∈ [0, 2*ny_half)` with the symmetric variant mirroring about the fault layer). Default `ny_half = 2–3`. Then re-confirm the asym run still clears the acceptance bars, and that the asym/sym contrast (Part B) widens with thickness (a positive control that the leak is friction-driven, not absorber-driven).

**Test case:** n/a for a hard pass/fail (POSSIBLE — would require running the harness, which policy forbids locally). Covered qualitatively by re-running the §7 A/B on Frontera with `ny_half ∈ {1,2,3}` and checking the drift is thickness-insensitive.

---

### [R-004] LOW [test_fault_planar_serial.cpp:VerifyMirrorSymmetry] — `std::array` is used without including `<array>` (compiles only by transitive include)

**Category:** QUALITY (latent include hygiene)

**Description:**
`VerifyMirrorSymmetry` declares `std::vector<std::array<real_t,3>> cen(ne);` (line 276) but the translation unit's explicit includes (lines 65–73) omit `<array>`. It compiles today only because `wave_operator.hpp` (included at line 58) pulls in `<array>` transitively. A future refactor that drops that transitive path breaks this file with a confusing error far from the cause.

**Suggested fix:**
```diff
 #include <algorithm>
+#include <array>
 #include <cmath>
```

**Test case:** n/a (compile-time hygiene).

---

### [R-005] LOW [Makefile:test-fault-planar-serial / test_fault_planar_serial.cpp:main] — the no-arg acceptance path runs ~16 000 ADER steps, far heavier than the rest of the unit-test suite

**Category:** QUALITY (test-suite fit)

**Description:**
The default run executes Part A (8000 steps on the 8×8×2 = 768-tet slab) **plus** the Part B symmetry probe (two 4000-step runs on 6×6×2 = 432-tet slabs) = ~16 000 ADER-2 macro-steps. `miniapps/seas/CLAUDE.md` budgets the unit suite at "~2 min total"; this single target plausibly dominates or exceeds that. It is closer to an integration test than a unit test. (Not a correctness defect — flagged only so the CI/`make test` wiring is a deliberate choice.)

**Suggested fix:** either give the no-arg / Makefile path smaller defaults (e.g. a CI mode `--nsteps 2000 --no-symmetry-probe` for `make test`, reserving the full 8000-step A/B for an explicit invocation), or keep this target out of the fast `make test` aggregate and document it as a heavier standalone reproduction. Example:
```diff
 test-fault-planar-serial: seas_test_fault_planar_serial
-	./seas_test_fault_planar_serial
+	./seas_test_fault_planar_serial            # full A/B reproduction (heavy)
+# CI-fast variant (kept out of the 2-min `make test` aggregate):
+# ./seas_test_fault_planar_serial --nsteps 2000 --no-symmetry-probe
```

**Test case:** n/a (suite-timing concern).

---

## Notes (reviewed, NOT findings)

- **Mesh generation is correct.** Vertex `vid` ↔ `AddVertex` ordering matches; the 6-tet Kuhn split + `orient` flag is lifted verbatim from `d4_tet_mesh.hpp`; fault detection (`|cy − he| < 1e-8` on interior faces) is exact (all three fault-face vertices lie at `y = he` exactly) and cannot misfire on non-fault faces (a single tet's vertices stay within one y-layer, so only genuine `y=he` faces average to `he`); exterior→attr 5, interior-non-fault left untagged. `CheckElementOrientation(true)` only swaps vertex slots, preserving the vertex/centroid positions that `VerifyMirrorSymmetry` checks.
- **Part B asym/sym discriminator is sound.** The vertex grid is mirror-symmetric for *both* variants, so the discriminator is the centroid set. I verified by hand that the orient(+1) Kuhn split's local tet centroids are **not** y-reflection-closed in full 3D (e.g. T1's reflection `(0.5,0.25,0.25)` has no partner in any column), so the asymmetric mesh fails the mirror check (`asym_not` holds) while the orient(−1)-upper symmetric mesh passes. The `rel_diff > 1e-3` assertion relies on the symmetric mesh producing measurably less drift; because vertex swaps are an exact change of basis, the symmetric mesh's *physical* solution is mirror-symmetric to round-off (the d4 "slot-wise equivariance" caveat is about a stricter DOF-level probe, not the physical solution), so the leak there is round-off-small and `rel_diff ≈ 1`. Reasonable.
- **The `--fault-overint` / `--fault-resample` knobs are parsed but unused.** This is **intended** per plan §6 (inert in Phase 0, wired in Phases 1–3) and is announced at runtime. Not a defect.
- **`MaxVnJumpInterior` / `EvalComp`** use the correct component-major layout `Q[c*ndof_total + dof]`, correct `±y` side assignment by element-centroid `y`, and `abs` makes the jump sign-agnostic. Correct.
- **MPI correctness:** per-step `GlobalMax` and the `do_print`-gated `MaxVnJumpInterior` reduction are collective on conditions identical across ranks (`verbose=true`, step-only predicates), so no rank-divergent deadlock. `TEST_TRUE` increments `num_tests` on all ranks but only rank 0's pass/fail counts feed the bcast exit code — consistent.
- **Makefile** new target mirrors the sibling `seas_test_rupture_multistep_serial_vs_parallel` object list and (correctly) appends `$(SEAS_POST_LINK_DEDUP_RPATH)` per the convention at Makefile:32 — the link succeeds (binary present).

---

## Summary
- Critical issues: 0
- Moderate issues: 1 (R-001 — no CFL guard; `--he`/`--dt` decoupled → silent ADER blow-up that can masquerade as a "reproduction")
- Low issues: 4 (R-002 doc/code default drift; R-003 one-element-thick slab vs the symptom it must isolate; R-004 missing `<array>`; R-005 heavy unit-test runtime)
- Plan compliance: **FULL** — all §6 Phase-0 requirements are present (serial + MPI-capable, homogeneous, planar fault, ≳20–40 fault faces, `max|σ_n−σ_n0|` and `max|[[v_n]]|` outputs, growing-excursion acceptance at p1, asymmetric-vs-mirror-symmetric variant). The two metrics are correctly sourced and the harness builds and links.
- Verdict: **PASS WITH FIXES.** No correctness bug in the harness logic and it compiles/links cleanly. R-001 should be fixed before the Phase-1/2/3 A/B sweeps begin, because an ungated `dt` lets a mis-sized mesh report instability as a drift and would corrupt the over-integration/resample comparison the harness exists to make. R-002–R-005 are low-risk cleanups.

## Unreviewed Areas
- **Runtime behavior / actual reproduction of the drift.** The harness was not executed (policy: no local TPV102 reproducer runs). Whether the default `8×8×2`, 8000-step, 1-element-thick configuration actually clears `excursion_peak > 1e5` and `final > 3×early` — and whether the absorbing-BC proximity (R-003) lets it — is unverified and must be confirmed on the §7 Frontera A/B (or an approved local run on a small fixture).
- **ADER-2 CFL constant.** R-001's `cfl=0.5` guard value is a conservative placeholder; the precise stable constant for the Kuhn-tet geometry was reasoned from the sibling test's validated `dt=1e-4`@`he=1000`, not measured. Calibrate against `ComputeMaxDt` output at run time.
- **FaultBasis sign/frame at `y = he`.** Assumed identical to the validated `y=0` TPV102 path (same `±y` normal, `|n·n_ref| = 1`, planar ⇒ frame exact per plan §1); not re-derived here.
