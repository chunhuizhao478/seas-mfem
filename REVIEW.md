# Code Review: TPV205 production-readiness audit + dev sbatch generation — 2026-04-27

## Review Scope
- Plan: `miniapps/seas/document/system_dev/tpv205_lsw_native_fields_plan_2026-04-27.md` (1087 lines).
- Prior review: REVIEW.md round 4 (R-001..R-005) closed by /code-fix in this session.
- Files audited (post-fix state):
  - `dynamic/fault_face_flux.{hpp,cpp}` — `EvaluateADER_LSW` (slip accumulation removed per R-001).
  - `dynamic/tpv205_friction.hpp`, `dynamic/tpv205_setup.hpp`, `dynamic/tpv205_substep_iterator.{hpp,cpp}`.
  - `drivers/tpv205_driver.cpp` — banner / disable_nucleation / paraview_write / SubStep init / time loop / diag-dump.
  - `dynamic/wave_operator.{hpp,inl}` — `FaultFrictionLaw` enum + dispatch sites at 3617–3631 (interior) and 4539–4553 (R-1601 shared-fault fallback) + `VerifySharedFaultDOFDataConsistency`.
  - `tests/unit/test_tpv205_evaluate_ader_lsw_parity.cpp`, `tests/verification/test_tpv205_mpi_rupture_crossing.cpp` (post-fix; `mpirun -np 2` reports R-001 ratio = 1.220 < 1.5 ✓).
  - Existing `jobs/tpv205/tpv205_mixed_flux_adjacent_200m_p2_O3_normal.sbatch` (P=2 + O=3, normal queue, 48h).
  - Reference template: `jobs/tpv102/tpv102_mixed_flux_adjacent_200m_p1_O2_dev.sbatch`.
- Domain context: `miniapps/seas/CLAUDE.md`, the tpv205 plan, prior REVIEW.md rounds, the SCEC TPV5 spec PDF (`tpv205/benchmark_document/TPV5_forwebsite.pdf`).

## Deliverable created
- **`miniapps/seas/jobs/tpv205/tpv205_mixed_flux_adjacent_200m_p1_O2_dev.sbatch`** — new dev-queue (8N × 400r × 2h) sbatch mirroring the TPV102 mixed-flux dev pattern, adapted for TPV205 (LSW friction, 4-patch pre-stress, 16-station SCEC layout, `tpv2053d_200m.msh`). Post-run summary parses the `x2_0_x3_7.5` hypocenter station and now also checks the `friction_solver_actual=lsw-closed-form` dispatch banner — a missing/`brent` value flags the R-001/R-016 fix as regressed.

## Findings

### [R-001] [LOW] [drivers/tpv205_driver.cpp:1713,2188–2199] — `psi` column in fault_qp_dump and ParaView comment is stale (TPV205 has no state variable)

**Category:** QUALITY / documentation

**Description:**
The end-of-run diagnostic dump (`SEAS_DIAG_DUMP_FAULT_QPS=1`) writes `d.psi` as the 12th column and labels the column header "psi" at line 2188:
```cpp
fs << "# columns: dof_idx x y z V1 V2 slip1 slip2 "
      "tau1_corr tau2_corr sigma_n_corr psi\n";
…
fs << … << " " << d.sigma_n_corr << " " << d.psi << "\n";
```
TPV205 uses linear slip-weakening, which has no state variable. `d.psi` is defensively zeroed at init (`tpv205_setup.hpp:135`) and never written. Every row of the dump therefore reports `psi = 0`, which (a) is misleading to anyone analysing the dump and (b) hides the actual TPV205-relevant per-QP scalar — μ_eff(δ).

Same channel mismatch appears in the ParaView setup comment at line 1713: "fault-surface PVD/VTU (slip, slip_rate, traction dip+strike, **psi**, sigma_n, …)". The actual ParaView state field carries `LSWFrictionCoefficient_TPV205(δ, μ_s, μ_d, d_c)` (driver line 1945–1948 — the R-005 fix from a prior round), not psi.

**Trigger:**
Run with `SEAS_DIAG_DUMP_FAULT_QPS=1`. Inspect `fault_qp_dump_rank<R>.txt`: every row's last column is 0.

**Actual behavior:**
`# columns: ... sigma_n_corr psi` and a column of zeros.

**Expected behavior:**
`# columns: ... sigma_n_corr mu_eff(delta)` where the dump computes `LSWFrictionCoefficient_TPV205(sqrt(slip1²+slip2²), d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c)` and writes that value.

**Suggested fix:**
```diff
--- a/miniapps/seas/drivers/tpv205_driver.cpp
+++ b/miniapps/seas/drivers/tpv205_driver.cpp
@@
             fs << std::scientific << std::setprecision(17);
             fs << "# columns: dof_idx x y z V1 V2 slip1 slip2 "
-                  "tau1_corr tau2_corr sigma_n_corr psi\n";
+                  "tau1_corr tau2_corr sigma_n_corr mu_eff_delta\n";
             const int n = static_cast<int>(dof_data.size());
             for (int i = 0; i < n; i++)
             {
                const Vector &c = fault_coords[i];
                const DOFData &d = dof_data[i];
+               const real_t delta = std::sqrt(d.slip1 * d.slip1
+                                              + d.slip2 * d.slip2);
+               const real_t mu_eff = mfem::seas::LSWFrictionCoefficient_TPV205(
+                                        delta, d.lsw_mu_s, d.lsw_mu_d,
+                                        d.lsw_d_c);
                fs << i << " " << c(0) << " " << c(1) << " " << c(2)
                   << " " << d.V1 << " " << d.V2
                   << " " << d.slip1 << " " << d.slip2
                   << " " << d.tau1_corr << " " << d.tau2_corr
-                  << " " << d.sigma_n_corr << " " << d.psi << "\n";
+                  << " " << d.sigma_n_corr << " " << mu_eff << "\n";
             }
```
Update the L1713 comment from "psi" to "mu_eff" likewise.

**Test case:**
None — diagnostic-only dump; the column header rename is verified by grep:
```bash
grep -n "psi\b" miniapps/seas/drivers/tpv205_driver.cpp
# expected: zero matches in the diag-dump and ParaView-comment blocks
# (only references in inherited TPV104 boilerplate / pre-existing comments
# that have nothing to do with TPV205 LSW).
```

---

### [R-002] [MODERATE] [POSSIBLE] [drivers/tpv205_driver.cpp:2138; dynamic/wave_operator.inl:5088–5466] — `VerifySharedFaultDOFDataConsistency()` after step 0 is uncovered for TPV205 at np > 1

**Category:** ASSUMPTION / test gap

**Description:**
The driver calls `wave.VerifySharedFaultDOFDataConsistency()` after the first macro-step (`tpv205_driver.cpp:2138`). The verify cross-rank-pairs every shared fault QP and aborts if any of the eight fields `{tau1_corr, tau2_corr, sigma_n_corr, V1, V2, psi, slip1, slip2}` differs by more than `tol = 1e-10` (relative).

The TPV205 dispatch on shared-fault QPs (R-1601 fallback, `wave_operator.inl:4539–4553`) calls `EvaluateADER_LSW(fdata, …)` independently on each rank that owns the shared QP. The function reads rank-local `I_plus_local / I_minus_local` derived from `elem1_on_plus = !qpd.sign_flipped`. The R-1601 comment (the ORIGIN of the SHARED FALLBACK) explicitly documents that "elem1_on_plus is rank-local while qpd.sign_flipped is rank-independent" — exactly the configuration that would make rank A and rank B compute opposite-sign V1/V2 for the same physical shared QP, which would trip the verify with relative diff ≈ 2.0.

The new `test_tpv205_mpi_rupture_crossing.cpp` does NOT call `VerifySharedFaultDOFDataConsistency()` — it bypasses this check. The TPV205 driver does call it (step==0), so a production np ≥ 2 run could abort on the first step with R-101 unpaired-entry messages and the user would not be able to predict it from the test suite.

**Why this is POSSIBLE rather than CRITICAL**: TPV104 production np ≥ 4 runs use the same R-1601 inline-`EvaluateADER` fallback on shared faces and have not reported verify aborts. Either (a) the rank-frame logic actually produces consistent output through `Evaluate` despite the comment's framing concern, OR (b) the verify tolerance happens to absorb whatever rank-disagreement does exist for rate-and-state outputs. Whether LSW outputs follow the same pattern is unverified — the only TPV205 MPI test bypasses the verify.

**Trigger:**
`mpirun -np N seas_tpv205_driver --mesh tpv205/mesh/tpv2053d_200m.msh ...` for N ≥ 2 with a partition that produces shared fault faces. The driver's first time-step completes, then `VerifySharedFaultDOFDataConsistency()` runs and may abort with `[R-101 rank-0 detail]` and `[UNPAIRED]` messages (or, on a less severe rank-mismatch, a value-tolerance abort with `max_rel_diff > 1e-10`).

**Actual behavior:**
Untested. The np=2 rupture-crossing test does not exercise this code path.

**Expected behavior:**
The verify must pass at np > 1 on a TPV205 mesh with shared-fault QPs (just as it does on TPV102/TPV104 production), OR the dispatch must be patched to write rank-consistent V1/V2/tau*_corr/slip* on shared QPs before the verify fires.

**Suggested fix:**
Extend `test_tpv205_mpi_rupture_crossing.cpp` to call the verify explicitly so any latent rank-mismatch surfaces in CI rather than at the first sbatch:

```diff
--- a/miniapps/seas/tests/verification/test_tpv205_mpi_rupture_crossing.cpp
+++ b/miniapps/seas/tests/verification/test_tpv205_mpi_rupture_crossing.cpp
@@
    AdvanceADERWithSubStep(wave, iterator, dof_data, fault_coords,
                           Q, dt, /*ader_order*/2, /*t_start*/0.0,
                           Q_new);
+
+   // R-002 (final review): the production driver calls
+   // VerifySharedFaultDOFDataConsistency() after step 0 (see
+   // tpv205_driver.cpp:2138).  Replicate that here so any latent
+   // rank-canonical-frame mismatch in EvaluateADER_LSW surfaces in
+   // CI rather than at the first sbatch on Frontera.  Default tol
+   // = 1e-10 (relative).
+   try
+   {
+      wave.VerifySharedFaultDOFDataConsistency();
+      num_tests++;
+      num_passed++;
+      if (rank == 0)
+      {
+         std::cout << "  PASSED: VerifySharedFaultDOFDataConsistency "
+                   << "after one ADER-O2 step (default tol 1e-10)\n";
+      }
+   }
+   catch (...)
+   {
+      // MFEM_ABORT cannot actually be caught — this branch is
+      // documentation only; if the verify trips it calls MPI_Abort
+      // and the whole test fails loudly.  Keeping the try/catch
+      // makes the intent explicit for future maintainers.
+      num_tests++;
+      num_failed++;
+   }
```
If the assertion trips on a real np > 1 run, the fix is in `EvaluateADER_LSW` — it must compute rank-consistent V1/V2 (e.g. by always orienting (V1, V2) along the canonical-+ tangent direction regardless of which rank owns Elem1). That fix is OUT OF SCOPE for /code-fix this round; the immediate action is to surface the issue in CI.

**Test case:**
The diff above is the test. To validate without applying the fix:
```bash
# Build current TPV205 binary, run np=2 mpi test, observe whether
# VerifySharedFaultDOFDataConsistency() trips:
make seas_test_tpv205_mpi_rupture_crossing -j8
mpirun -np 2 ./seas_test_tpv205_mpi_rupture_crossing
# Expected (post-extension): if the verify trips, the test exits
# with the [R-101 ...] abort message and a non-zero exit code;
# if it passes, the new test slot reports "PASSED:
# VerifySharedFaultDOFDataConsistency …".
```

---

### [R-003] [LOW] [tests/unit/test_tpv205_friction.cpp + Makefile] — Friction-helper unit test exists but has no Make rule, so it never runs in CI

**Category:** QUALITY / test wiring

**Description:**
`tests/unit/test_tpv205_friction.cpp` (203 lines) covers the strength-barrier short-circuits R-002/R-003 from `tpv205_friction.hpp` (`LSWFrictionCoefficient_TPV205`, `SolveLSW_TPV205`). The file is present in the repo but the Makefile has no `seas_test_tpv205_friction` target. `make test` and the `test-tpv205` aggregate (if any) do not exercise it, so a regression in `LSWFrictionCoefficient_TPV205` (e.g. removing the `mu_s >= 0.5 * mu_s_barrier` short-circuit) would compile and link, would not break the parity test or the MPI test, and would only surface as a wrong rupture-area boundary on Frontera.

**Trigger:**
`make seas_test_tpv205_friction` returns `make: *** No rule to make target 'seas_test_tpv205_friction'. Stop.` (verified during the build step of /code-fix in this session).

**Actual behavior:**
The test never runs.

**Expected behavior:**
Wire it into the Makefile alongside the existing TPV205 test targets at `Makefile:2851–2868`.

**Suggested fix:**
```diff
--- a/miniapps/seas/Makefile
+++ b/miniapps/seas/Makefile
@@
 # REVIEW R-016 Phase 4 — TPV205 LSW dispatch tests.
 TEST_TPV205_EVALUATE_ADER_LSW_PARITY_SRC = tests/unit/test_tpv205_evaluate_ader_lsw_parity.cpp
 TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ = $(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_SRC:.cpp=.o)
 TEST_TPV205_MPI_RUPTURE_CROSSING_SRC = tests/verification/test_tpv205_mpi_rupture_crossing.cpp
 TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ = $(TEST_TPV205_MPI_RUPTURE_CROSSING_SRC:.cpp=.o)
+TEST_TPV205_FRICTION_SRC = tests/unit/test_tpv205_friction.cpp
+TEST_TPV205_FRICTION_OBJ = $(TEST_TPV205_FRICTION_SRC:.cpp=.o)
@@
 $(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
 	@mkdir -p $(@D)
 	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@
+
+$(TEST_TPV205_FRICTION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
+	@mkdir -p $(@D)
+	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@
@@
 seas_test_tpv205_mpi_rupture_crossing: $(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ) $(TPV205_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ)
 	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ) $(TPV205_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(MFEM_LIBS)
+
+seas_test_tpv205_friction: $(TEST_TPV205_FRICTION_OBJ)
+	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV205_FRICTION_OBJ) $(MFEM_LIBS)
```

(The test depends only on `tpv205_friction.hpp` and `tpv205_params.hpp` — both header-only — plus `test_macros.hpp`, so no extra `.o` files are needed beyond MFEM itself.)

**Test case:**
```bash
make seas_test_tpv205_friction -j4 && ./seas_test_tpv205_friction
# Expected: 3 sub-tests pass — barrier_mu_no_collapse,
# barrier_locked_under_tensile, combined_barrier_stays_locked.
```

---

### [R-004] [LOW] [tpv205/mesh/] — `.msh` mesh files are not in the repo; production sbatch must build them on Frontera

**Category:** ASSUMPTION / production gap

**Description:**
`tpv205/mesh/` contains only `.geo` source (`tpv2053d_100m.geo`, `tpv2053d_200m.geo`); no pre-built `.msh` exists. Both the existing normal-queue sbatch and the new dev sbatch generated this round include an explicit `[ ! -s tpv205/mesh/tpv2053d_200m.msh ] && exit 1` guard, so a missing mesh fails the job loudly rather than silently producing empty output. But there is no helper sbatch / script in `jobs/tpv205/` that runs gmsh on the `.geo` to produce the `.msh` — the user must do this manually on Frontera, which is fragile.

TPV104 has `jobs/tpv104/tpv104_mesh_build.sbatch` for exactly this pattern. TPV205 lacks the equivalent.

**Trigger:**
Submit either TPV205 sbatch on a fresh Frontera workspace where the gmsh step has not been run.

**Actual behavior:**
The sbatch's `[ ! -s tpv205/mesh/tpv2053d_200m.msh ]` guard fires:
```
ERROR: tpv205/mesh/tpv2053d_200m.msh missing.  Run gmsh on tpv205/mesh/tpv2053d_200m.geo first.
```
Job aborts immediately. No data corruption, but a wasted submission.

**Expected behavior:**
A `jobs/tpv205/tpv205_mesh_build.sbatch` mirroring `jobs/tpv104/tpv104_mesh_build.sbatch` so the user can submit it once and get all standard `.msh` files built.

**Suggested fix:**
Out of scope for this review (the user asked for a code-running sbatch, not a mesh-build sbatch). Note as a follow-up: copy `jobs/tpv104/tpv104_mesh_build.sbatch` to `jobs/tpv205/tpv205_mesh_build.sbatch`, swap the gmsh inputs from `tpv104/mesh/*.geo` to `tpv205/mesh/*.geo`. The dev sbatch already documents the dependency.

**Test case:**
None — meta-tooling.

---

### [R-005] [LOW] [drivers/tpv205_driver.cpp:1664–1677] — Background banner reports `V_ini = 0 m/s` but does not echo the strength-barrier sentinel `μ_s = 10000`

**Category:** QUALITY

**Description:**
The startup banner at lines 1664–1677 prints `μ_s = 0.677 (rupture area), 10000 (strength barrier)` and `V_ini = 0 m/s`. That's correct content but the dev-queue post-run summary template parses `peak |σ_n − 120 MPa|` from the hypocenter station file as the dominant production sanity check. There is no **equivalent banner field** for `μ_s_barrier` so a future driver edit that overrode `μ_s_barrier` (e.g. via a plan that introduces a CLI override) could silently drop the barrier without showing up in the banner or the post-run summary.

This is a defensive-design comment — not a bug. Severity LOW because the existing `LSWFrictionCoefficient_TPV205` short-circuit + `SolveLSW_TPV205` short-circuit + `test_tpv205_friction` (once it's wired per R-003) cover the barrier semantics. The banner phrasing is fine for the dev sbatch.

**Suggested fix:**
None (note for future maintainers).

**Test case:**
None.

---

## Summary
- Critical issues: 0
- Moderate issues: 1  ([R-002] VerifySharedFaultDOFDataConsistency uncovered by tests for TPV205 — POSSIBLE)
- Low issues: 4  ([R-001] stale `psi` column, [R-003] unwired friction test, [R-004] mesh `.msh` not in repo, [R-005] banner echo gap)
- Plan compliance: FULL (all four phases of `tpv205_lsw_native_fields_plan_2026-04-27.md` landed; round-5 R-001..R-005 fixes verified by `seas_test_tpv205_evaluate_ader_lsw_parity` 14/14 and `mpirun -np 2 seas_test_tpv205_mpi_rupture_crossing` 6/6).
- Verdict: **PASS WITH FIXES** — production-ready for the dev-queue smoke run; address R-002 (test-coverage gap) before submitting np ≥ 64 production jobs.

## New artefact

`miniapps/seas/jobs/tpv205/tpv205_mixed_flux_adjacent_200m_p1_O2_dev.sbatch` (252 lines).

The dev sbatch mirrors `jobs/tpv102/tpv102_mixed_flux_adjacent_200m_p1_O2_dev.sbatch` line-for-line in structure (modules, `LD_LIBRARY_PATH`, build, `ibrun`, post-run summary) with the following TPV205-specific adaptations:

| TPV102 dev sbatch | TPV205 dev sbatch |
|---|---|
| driver `seas_tpv102_driver` | `seas_tpv205_driver` |
| mesh `tpv102/mesh/tpv102_200m.msh` | `tpv205/mesh/tpv2053d_200m.msh` |
| output prefix `tpv102_mfadj_p1_O2` | `tpv205_mfadj_p1_O2` |
| `--fric-law aging` | `--fric-law lsw` |
| hypocenter station `flt_0_7.5` | `x2_0_x3_7.5` |
| station column 5 = V2 strike | column 3 = V2 strike (different column order — TPV205 SCEC trace layout) |
| no LSW dispatch verdict | adds `friction_solver_actual=lsw-closed-form` post-run gate (R-016 regression sentinel) |
| no mesh existence guard | `[ ! -s tpv2053d_200m.msh ] && exit 1` (mesh `.msh` is not in the repo — see R-004) |

Same modules (`intel/19.1.1`, `impi/19.0.9`, `hypre/2.31.0`, `mumps/5.3`, `parmetis`, `petsc/3.15`, `fftw3/3.3.8`), same `LD_LIBRARY_PATH` order, same `8N x 400r x 2h` allocation, same `output-dt 0.05 / paraview-dt 0.5 / paraview-bulk-dt 0.5` cadence as the TPV102 dev pair.

## Unreviewed Areas
- `dynamic/wave_operator.inl` lines 1–3500 (~3500 lines of shared dispatch) — out of scope per `feedback_dynamic_folder_editable_for_tpv104`. Audited only the two TPV205 dispatch sites and the `VerifySharedFaultDOFDataConsistency` body for R-002.
- `tpv205/mesh/tpv2053d_*.geo` — assumed to emit Physical Surface 101/103/105 as the driver and `TPV205Params::bc_*_default` advertise; not parsed.
- The `jobs/tpv205/tpv205_mixed_flux_adjacent_200m_p2_O3_normal.sbatch` (existing) was inspected for pattern reuse but its production verdict is out of scope; the new dev sbatch is its low-cost counterpart for early-stage smoke testing.
- The `.msh` build flow (gmsh on `.geo` → `.msh`) — see R-004; a follow-up `jobs/tpv205/tpv205_mesh_build.sbatch` would close the gap but is out of scope for this review.
