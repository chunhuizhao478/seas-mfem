# Code Review: SEAS ADER Phase-0 optimization (2026-06-24)

## Review Scope
- Plan: miniapps/seas/document/code_optimization_dev/action_plan_2026-06-24.md
- Files reviewed: dynamic/wave_operator.hpp, dynamic/wave_operator.inl, tests/parallel/test_bimaterial_deriv_cache_parity.cpp, Makefile, jobs/.../...upwind_ader_caliper...sbatch
- Domain context: CLAUDE.md (R-1505 deadlock-prevention, collective-call discipline), the profiling analysis doc
- Verification performed during synthesis: read of wave_operator.inl:3845-3884 and :4910-4959 (both correctors), wave_operator.hpp:386-411 (setter + re-arm) and :1075-1083 (mutable flag docstring), test fixture :100-124, plus `grep` of the Makefile gate target and the driver's `use_shared_ck` flag.

## Findings

### [R-001] [MODERATE] [POSSIBLE] [dynamic/wave_operator.inl:4916 → :4938-4953 (ComputeADERSharedFaceFluxRHS)] — Shared-corrector full-communicator consensus Allreduce sits behind the `if (n_shared==0) return;` early return; a non-uniform shared-face count deadlocks the communicator

**Category:** EDGE_CASE / latent-deadlock (PRE-EXISTING from R-1505 — not introduced by the hoist)

**Corroboration:** Raised by all four reviewers (dim 1 finding 1 = LOW, dim 2 finding 1 = CRITICAL, dim 4 finding "Pre-existing latent deadlock" = LOW). Severity is reconciled to MODERATE here: the deadlock is real and reproducible, but **the offending placement predates the hoist** — three independent reviewers confirmed via `git diff HEAD` / `git show HEAD` that the `n_shared==0` early return preceded the consensus Allreduce in the original per-step code, and step 1 of the hoisted code is byte-identical. The hoist therefore does **not** introduce, worsen, or fix it. It is downgraded from the lone CRITICAL vote because it is not a regression of this change; it is upgraded above LOW because 0C is explicitly a "guard-Allreduce cleanup pass" (plan §0C line 98: "one correctly-scoped pass over all four guard reductions") that had the opportunity — and arguably the mandate — to remove this landmine and instead enshrined it as a design justification (hpp:1079-1081, inl:4931-4937).

**Description:** The shared corrector's one-time consensus is a whole-communicator `MPI_Allreduce(..., pmesh.GetComm())` (inl:4943) reached only by ranks that pass `if (n_shared == 0) { return; }` (inl:4916). A whole-comm collective must be entered by every rank in the communicator. The codebase itself documents (R-1600, inl:2189-2202) that `GetNSharedFaces()` can be 0 on some ranks and >0 on others, and that gating a collective on it caused the production-blocking 13.5-min np=10 hang on the symmirror 1000m mesh. The R-1600-safe pattern (RK4 `ComputeSharedFaceFluxRHS`, inl:3290) returns on `n_shared==0` and then issues **only** pairwise `q_gf.ExchangeFaceNbrData()` (neighbour-set, safe to skip with no peers) — never a full-comm Allreduce. The ADER shared corrector is the sole site that places a full-communicator reduction behind the `n_shared==0` gate.

The new comments overstate the safety of the delta: inl:4923-4927 ("Only the shared corrector's ExchangeFaceNbrData below would actually deadlock") and inl:4935-4937 ("trivially deadlock-free") are accurate for the *delta* (step 1 identical, steps 2..n run no Allreduce) but mislead a reader into thinking the shared-corrector consensus is unconditionally safe. It inherits, unchanged, an unenforced "every rank has ≥1 shared face" assumption.

**Trigger:** Any np≥2 ParMesh partition where at least one rank has `pmesh.GetNSharedFaces()==0` while another has >0 (e.g. an island/fully-interior partition), on the first AdvanceADER step or any step after a mid-run SetAbsorbingBackground re-arm.

**Actual behavior:** Ranks with `n_shared>0` block forever in the shared-corrector `MPI_Allreduce`; ranks with `n_shared==0` return early and proceed — hard hang. (Note: the local corrector at inl:3860 runs its Allreduce unconditionally with no `n_shared` gate, so all ranks rendezvous there first; a violated invariant therefore hangs at the *shared* Allreduce, both pre- and post-hoist, identically.)

**Expected behavior:** A full-communicator reduction must be reached by all ranks unconditionally, exactly like the local corrector's gate-free consensus at inl:3860.

**Suggested fix:** Delete the shared-corrector consensus Allreduce entirely and rely on the local corrector's gate-free, all-ranks check, which `AdvanceADER` always runs immediately before the shared corrector (inl:5762/5766). The local consensus already fails loud on **all** ranks before any `ExchangeFaceNbrData` is reached, so the shared-corrector Allreduce is redundant for fail-loud purposes and only adds a participation requirement the surrounding `n_shared` gate does not honor.

```diff
--- a/dynamic/wave_operator.inl
+++ b/dynamic/wave_operator.inl
@@ -4914,4 +4914,12 @@ ComputeADERSharedFaceFluxRHS(...)
       auto &pmesh = static_cast<const ParMesh &>(mesh_);
       int n_shared = pmesh.GetNSharedFaces();
       if (n_shared == 0) { return; }
-      {
-#ifdef MFEM_USE_MPI
-         if (!bulk_bg_consensus_shared_done_)
-         {
-            ... MPI_Allreduce(..., pmesh_consensus.GetComm()); ...
-            bulk_bg_consensus_shared_done_ = true;
-         }
-#else
-         MFEM_VERIFY(has_bulk_bg_, ...);
-#endif
-      }
+      // R-1505 / R-1600: the all-ranks gate-free consensus already runs in
+      // ComputeADERFaceFluxRHS (called immediately before this corrector each
+      // AdvanceADER step) and fails loud on ALL ranks before any pairwise
+      // ExchangeFaceNbrData. Do NOT run a full-comm Allreduce behind the
+      // n_shared==0 early return — a rank with no shared faces would skip it
+      // and deadlock the rest (the R-1600 np=10 hang shape).
+      MFEM_VERIFY(has_bulk_bg_,
+                  "ComputeADERSharedFaceFluxRHS requires SetAbsorbingBackground"
+                  "(Q_bg); the collective consensus is the local corrector's. "
+                  "See R-1505/R-1600.");
```
Then remove `bulk_bg_consensus_shared_done_` and its re-arm at hpp:410, and drop the two-flag justification in hpp:1079-1083 (a single `bulk_bg_consensus_done_` driven by the gate-free local corrector suffices). If a self-documenting collective check is wanted in the shared corrector, it MUST be hoisted **above** the `if (n_shared==0) return;` so all of `GetComm()` participates — but deleting it is cleaner and is what the local-corrector net makes safe.

**Test case:** Construct an np≥3 ParMesh where one rank owns a fully-interior partition (`n_shared==0`) while others share a seam; call SetAbsorbingBackground on all ranks, then AdvanceADER. Current and pre-hoist code both hang at the shared-corrector Allreduce (same shape as the R-1600 np=10 symmirror hang). After the fix, only the all-ranks local check runs and the run proceeds (or fails loud uniformly if a rank truly skipped the setter). This is Expanse/multi-rank-only; a construction-order fixture forcing zero shared faces on one rank reproduces it.

---

### [R-002] [MODERATE] [NOT REACHABLE TODAY] [dynamic/wave_operator.inl:3860-3874 & :4938-4953; hpp:393-411] — Caching the consensus MASKS mid-run divergence the original per-step check caught; the new collective-call contract on SetAbsorbingBackground is undocumented and the plan-mandated debug-invariance assert was not implemented

**Category:** DEVIATION (weakened safety guarantee, dormant in production)

**Corroboration:** Raised by dim 2 (finding 2, MODERATE) and dim 4 (finding "debug-build assert ... NOT added", LOW). Merged here.

**Description:** The original code re-ran `MPI_Allreduce(MPI_MIN, has_bulk_bg_)` every step, so it would fail loud on any step N>1 where `has_bulk_bg_` became non-uniform across ranks. The hoisted version verifies once and caches in `bulk_bg_consensus_{local,shared}_done_`, so it can no longer detect a mid-run divergence on steps 2..N. The re-arm in SetAbsorbingBackground (hpp:409-410) only re-detects divergence if **every** rank re-calls the setter. A subset of ranks calling `SetAbsorbingBackground(nullptr)` mid-run while others do not (a) clears `has_bulk_bg_` only on the calling ranks and (b) re-arms the flags only on those ranks, so the next step's gated Allreduce is itself non-uniformly entered — reintroducing partial-participation deadlock via the re-arm instead of detecting it. The original per-step code would have caught this same subset-clear at the next step's uniform Allreduce as a fail-loud `min_has==0` abort on all ranks.

This is **dormant, not a live regression**: all four production drivers (spatial_dyn_driver.cpp:1999, tpv102:1601, tpv104:1660, tpv205:1621) and the ~71 unit/parallel tests call SetAbsorbingBackground exactly once, uniformly, outside the time loop. But it is a real weakening of the R-1505 guarantee, and SetAbsorbingBackground now carries an undocumented collective-call contract analogous to the R-1205 footgun for SetMixedFluxMode. The docstring (hpp:386-411) still describes a plain inline setter and advertises nullptr-clear as supported, with no warning that an uneven mid-run call now risks deadlock instead of a fail-loud abort. Separately, plan §0C line 103 explicitly asked for "a debug-build assert that `has_bulk_bg_` is unchanged across the loop"; `grep` confirms no such assert exists.

**Trigger:** A future caller invokes SetAbsorbingBackground mid-run on a subset of ranks, or with a different nullptr/non-nullptr argument per rank, then advances. Not reachable through any current driver/test.

**Actual behavior:** Consensus is verified once and cached; a non-uniform mid-run SetAbsorbingBackground converts what the original would have caught as a uniform fail-loud abort into a potential partial-participation deadlock or silently diverged numerics. The plan-mandated loop-invariance assert is absent.

**Expected behavior:** Either keep the per-step re-check (original) so mid-run divergence stays fail-loud, or — if caching — document the new collective-call contract on SetAbsorbingBackground (mirroring R-1205) and add the plan-mandated debug-build invariance assert.

**Suggested fix:** Document the contract and add the debug assert (both zero numeric risk):

```diff
--- a/dynamic/wave_operator.hpp
+++ b/dynamic/wave_operator.hpp
@@
+   /// COLLECTIVE-CALL CONTRACT (R-1505 hoist, opt 2026-06-24): if called AFTER
+   /// the first AdvanceADER, every rank MUST call this with the same
+   /// nullptr/non-nullptr argument, in the same order (same discipline as
+   /// SetMixedFluxMode, R-1205). An uneven mid-run call re-arms the consensus
+   /// flags non-uniformly and can deadlock the next AdvanceADER's gated
+   /// Allreduce instead of failing loud. Production callers invoke this once at
+   /// setup, uniformly — keep it that way.
    void SetAbsorbingBackground(const real_t *Q_bg)
```
And in AdvanceADER (debug build only), record `has_bulk_bg_` on first entry and `MFEM_ASSERT` it unchanged on subsequent entries, so a stale cached consensus is caught locally. Note R-001's fix (collapsing to a single gate-free flag) eliminates half of this hazard's surface; the contract docstring is still warranted.

**Test case:** Document-and-assert change; the divergence path is unreachable in production today, so no functional test fails. A np2 unit test that calls `SetAbsorbingBackground(nullptr)` on rank 0 only after step 1 would demonstrate the lost fail-loud (original: abort; new: hang or silently-diverged numerics).

---

### [R-003] [MODERATE] [POSSIBLE] [tests/parallel/test_bimaterial_deriv_cache_parity.cpp:104-112] — Gate test uses a numerically HOMOGENEOUS material; plan §0B demanded a TWO-material fixture, so per-element flux-pool divergence is untestable

**Category:** DEVIATION (coverage gap)

**Description:** Plan §0B (action_plan line 88) and open-question 1 (line 181) explicitly require "a small TWO-material np2/np4 fixture." The implemented test builds the operator with three `ConstantCoefficient` instances (`lam_c`, `mu_c`, `rho_c` at line 110), i.e. a single spatially-constant material. This passes the driver's `mode != Constant` matrix-path guard (spatial_dyn_driver.cpp:1266) only because `Mode::Coefficient != Mode::Constant`, but the material is numerically identical on every element. Consequently every `FluxForElem_(e)` returns the same star matrix `A_d`, so the per-element flux-pool divergence the test claims to guard (comment lines 106-109: "every access still routes through the per-element flux pool") is never exercised. A pool index bug (off-by-one, transposed/swapped entry, or material-of-element-`e` applied to DOFs of `e'`) is invisible when all materials are equal. Both reference (OnTheFly) and tested (Cached) branches share `FluxForElem_`, so even the cached-vs-OnTheFly comparison cannot detect a per-element-A bug regardless of homogeneity. The production run uses a genuinely heterogeneous CVM velocity sidecar (Mode::Coefficient/GridFunction), so the production matrix-path divergence is not covered. (Geometric per-element coverage IS present — the cartesian-of-tets mesh gives geometrically distinct tets, so the D_d^e / S_d^e geometry operators do vary; the gap is strictly material heterogeneity.)

**Trigger:** Numerically-constant Coefficient material in the operator constructor.

**Actual behavior:** Single homogeneous material; all `A_d^e` identical; per-element pool divergence and indexing bugs cannot be caught.

**Expected behavior:** A two-material fixture so `A_d^e` differs across elements and a flux-pool indexing/divergence bug produces a non-zero Cached-vs-OnTheFly or fused-vs-separate delta.

**Suggested fix:** Replace the three `ConstantCoefficient` with position-dependent coefficients:

```diff
--- a/tests/parallel/test_bimaterial_deriv_cache_parity.cpp
+++ b/tests/parallel/test_bimaterial_deriv_cache_parity.cpp
@@
-   ConstantCoefficient lam_c(lam), mu_c(mu), rho_c(rho);
+   // TWO-material fixture (plan §0B): material A in x<0.5, B in x>=0.5, so
+   // FluxForElem_(e) is element-dependent and a flux-pool index bug becomes
+   // observable in Cached==OnTheFly / fused==separate deltas.
+   FunctionCoefficient lam_c([](const Vector &x){ return x[0] < 0.5 ? 32.04e9 : 48.0e9; });
+   FunctionCoefficient mu_c ([](const Vector &x){ return x[0] < 0.5 ? 32.04e9 : 40.0e9; });
+   FunctionCoefficient rho_c([](const Vector &x){ return x[0] < 0.5 ? 2670.0  : 3000.0; });
```

**Test case:** Inject a deliberate fault — in bimaterial_wave_operator.inl:387 change `FluxForElem_(e)` to `FluxForElem_((e+1)%ne_)`. With the current homogeneous fixture the test still PASSES (0 delta); with the two-material fixture it FAILS, demonstrating the current fixture is blind to per-element pool bugs.

---

### [R-004] [MODERATE] [POSSIBLE] [tests/parallel/test_bimaterial_deriv_cache_parity.cpp:103] — Gate test runs at FE order=2 (P2); production runs at FE order=1 (P1), so the deriv-cache regime under test is not the production regime

**Category:** DEVIATION (coverage gap)

**Description:** The deriv-cache geometry operators `elem_deriv_op_` (D_d^e = M_e⁻¹ K_d^e) and `elem_volume_op_` (S_d^e) are FE-order-dependent (built by BuildElementDerivativeOperators / BuildElementVolumeOperators at SetDerivMode(Cached), inl:1176-1178). The test fixes `const int order = 2;` (line 103), i.e. P2 (ndof_per_el_=10 for tets). The production config (`[mesh] order = 1`, driver passes `cfg.mesh.order` at spatial_dyn_driver.cpp:1247) runs P1 (ndof_per_el_=4, affine tets). The plan's premise (action_plan lines 13, 91) is that the cache win is largest precisely at P1, where OnTheFly does a per-QP Jacobian inverse and round-off re-association (the reason the tol is 1e-12 not 0.0) differs between the cached matvec and the OnTheFly quadrature. Validating Cached==OnTheFly only at P2 does not establish the ≤1e-12 contract at the P1 ndof=4 stride the production binary executes. (The `order_c=2..4` loop in Part C varies the ADER recursion order, not the FE space order, so it does not close this gap.)

**Trigger:** Hardcoded `order=2` mismatching production `order=1`.

**Actual behavior:** Only P2 is tested; P1 (production regime) is never built.

**Expected behavior:** Run the Cached-vs-OnTheFly equivalence at the production FE order (P1), or parametrize over {1,2}.

**Suggested fix:** Loop the operator FE order over {1,2}, rebuilding the operator each iteration:

```diff
-   const int order = 2;
+   for (int order : {1, 2})   // P1 is the production stride; P2 keeps higher-order coverage
+   {
    ...
    BimaterialWaveOperator<ParMesh> wave(pmesh, order, ..., bc);
    ... // Parts A and B
+   }
```
At minimum add an `order=1` instantiation for Parts A and B so the affine-P1 cached matvec vs per-QP-Jacobian-inverse OnTheFly path is validated at the production tol.

**Test case:** Construct BimaterialWaveOperator with `order=1`, rerun Parts A/B, confirm RelDiff ≤ 1e-12 — the regime the sbatch's `--deriv-cache` actually exercises.

---

### [R-005] [MODERATE] [POSSIBLE] [Makefile:4969] — Gate target is an orphan; the sbatch's "gated on `make test-bimaterial-deriv-cache-parity`" claim overstates a target that no aggregate runs

**Category:** EDGE_CASE (process / acceptance-criterion gap)

**Description:** The Makefile defines `test-bimaterial-deriv-cache-parity` (line 4969) and the build target `seas_test_bimaterial_deriv_cache_parity` (line 4961), but `grep` confirms the run target is referenced **nowhere else** — only its own definition. It is absent from `test:` (lines 3967-4033) and `test-parallel:` (line 4864). The sbatch comment (line ~189) claims the matrix path is "gated on `make test-bimaterial-deriv-cache-parity`", implying an active gate, but nothing runs it unless a human types that exact target. This mirrors the pre-existing orphan `test-wave-operator-cached-parallel` (also np2, unregistered), consistent with the project's convention of excluding np2-mpirun tests from the umbrella (cf. the explicit note at Makefile lines 4034-4038 for `test-paraview-rank0-warning-gate`). But action_plan acceptance criterion 1 (line 194) treats this test as a hard prerequisite ("Land this before trusting the flag"); an orphan target provides no standing protection against a future flux-pool/cache refactor.

**Trigger:** Target defined but never added to an aggregate.

**Actual behavior:** Orphan target with no aggregate reference and no in-Makefile note marking it an intentional manual gate; the sbatch comment overstates it as an active gate.

**Expected behavior:** The gate test is wired into an aggregate the dev/CI workflow runs, OR documented (like test-paraview-rank0-warning-gate at line 4034) as an intentionally-manual np2 pre-submit with the exact command.

**Suggested fix:** Either add it to a parallel aggregate, e.g.

```diff
--- a/Makefile
+++ b/Makefile
@@ test-parallel: ...
-test-parallel: <existing prereqs>
+test-parallel: <existing prereqs> test-bimaterial-deriv-cache-parity
```
OR add an explicit comment block near line 4969 (mirroring 4034-4038) stating it is an intentional manual np2 pre-submit and giving the exact `make test-bimaterial-deriv-cache-parity` command, so the sbatch's "gate" claim is truthful.

**Test case:** Run `make test` and `make test-parallel`, grep output for `bimaterial Cached==OnTheFly` — it never appears, confirming the gate does not run.

---

### [R-006] [MODERATE] [POSSIBLE] [tests/parallel/test_bimaterial_deriv_cache_parity.cpp:118] — 0C re-arm path is never exercised; SetAbsorbingBackground is called exactly once, so the flag-reset logic the diff added is untested

**Category:** EDGE_CASE (coverage gap)

**Description:** The 0C diff added a re-arm at the end of SetAbsorbingBackground (hpp:409-410: both `bulk_bg_consensus_*_done_ = false`) so a mid-run background change re-verifies consensus. The test calls `wave.SetAbsorbingBackground(Q_bg)` exactly once (line 118) before any AdvanceADER and never again, so the re-arm (calling the setter a second time after AdvanceADER steps and confirming the next corrector re-runs the Allreduce) is never exercised. If a future edit broke the re-arm (e.g. reset only one of the two flags), this test would still pass. The multi-step SKIP path IS covered: Part B runs AdvanceADER n_steps=3 per DerivMode, so step 1 sets the flag and steps 2-3 hit the `if(!flag)` skip. Only the re-arm half is untested.

**Trigger:** Single SetAbsorbingBackground call.

**Actual behavior:** Re-arm flag-reset path has zero coverage.

**Expected behavior:** Call SetAbsorbingBackground a second time mid-test and assert the next AdvanceADER still produces correct Cached==OnTheFly output, exercising the re-arm of both flags.

**Suggested fix:** Add a Part B′ after the first AdvanceADER loop: call `wave.SetAbsorbingBackground(Q_bg2)` (a different but still-valid background), then run another AdvanceADER and compare Cached vs OnTheFly. Document that a divergent-rank deadlock test (one rank skips the setter) is out of scope for an all-ranks-symmetric fixture. (If R-001's single-flag collapse lands, update this to assert the one remaining flag re-arms.)

**Test case:** Comment out the two re-arm lines at hpp:409-410 — the test still passes, proving the re-arm is uncovered.

---

### [R-007] [LOW] [NOT REACHABLE TODAY] [drivers/spatial_dyn_driver.cpp:1316-1320 (vs action_plan §0D)] — 0D deferred: global CK default still OFF/opt-in and the rank-0 path-log line is conditional, leaving Acceptance Criterion #4 (self-documenting reproducibility) only partially met

**Category:** DEVIATION (intentional deferral)

**Description:** Plan §0D calls for inverting the `use_shared_ck` default to true (with a `--no-shared-ck-recursion` opt-out) AND making the rank-0 path-log line unconditional so stdout always records which CK path ran (plan line 112; Acceptance Criterion #4 line 207). `grep` confirms neither was done: spatial_dyn_driver.cpp:1316 still reads `use_shared_ck = HasFlag(argc, argv, "--shared-ck-recursion")` (default OFF, opt-in), there is no `--no-shared-ck-recursion`, and the log at line 1317 is gated on `if (use_shared_ck && rank == 0)`, so a run without the flag prints nothing about the CK path. The task brief states 0D was deferred, so this is intentional and defensible (0D is pure ergonomics, zero numeric/perf delta, and deferring it keeps the OnTheFly/separate-CK baseline untouched). Consequence: the caliper sbatch's reproducibility relies on its own `echo Opt levers: ${OPT_ARGS}` line rather than the driver's stdout, and a future run that forgets the flag silently runs the un-fused path with no driver-side record.

**Trigger:** A future production run omits the flag; no driver-side record of which path executed.

**Actual behavior:** `use_shared_ck` default false (opt-in only); no opt-out flag; conditional log line.

**Expected behavior:** Plan §0D: default `use_shared_ck=true` with `--no-shared-ck-recursion` opt-out; unconditional rank-0 path log.

**Suggested fix:** Acceptable to leave deferred. If reproducibility (Acceptance #4) is to be claimed for the driver itself, make the rank-0 log unconditional now (one line, zero numeric risk):

```diff
-   const bool use_shared_ck = HasFlag(argc, argv, "--shared-ck-recursion");
-   if (use_shared_ck && rank == 0)
-   {
-      std::cout << "[opt] fused shared-CK recursion ...\n";
-   }
+   const bool use_shared_ck = HasFlag(argc, argv, "--shared-ck-recursion");
+   if (rank == 0)
+   {
+      std::cout << "[opt] shared-CK path: "
+                << (use_shared_ck ? "fused (one predictor recursion/step)"
+                                  : "separate (two recursions; baseline)")
+                << "; DerivMode=" << (HasFlag(argc, argv, "--deriv-cache") ? "Cached" : "OnTheFly")
+                << "\n";
+   }
```

**Test case:** Run the driver without `--shared-ck-recursion`, grep stdout for `shared-CK` — absent today, so the job log does not record that the un-fused path ran.

---

### [R-008] [LOW] [POSSIBLE] [tests/unit/test_lsw_rk_mixed_flux_bimaterial.cpp (new, untracked) + Makefile:233-234, 895, 5819-5859] — Scope creep: an unrelated mixed-flux Phase-4 test and its two Makefile targets are bundled into the same working tree as the 2026-06-24 optimization

**Category:** QUALITY (revertibility)

**Description:** The audited Phase-0 diff should be: wave_operator.{hpp,inl} (0C), the sbatch (0A/0B), and tests/parallel/test_bimaterial_deriv_cache_parity.cpp + its Makefile target (0B gate). The working tree also adds tests/unit/test_lsw_rk_mixed_flux_bimaterial.cpp and two Makefile targets (`seas_test_lsw_rk_mixed_flux_bimaterial`, `test-lsw-rk-mixed-flux-bimaterial`), whose own header attributes them to a DIFFERENT plan ("Phase 4 of PLAN_mixed_flux_seissol_port_2026-06-19.md"). These are additive/harmless (no production header touched, byte-exact contract preserved), but they are not part of the optimization plan. Bundling two unrelated efforts in one uncommitted diff complicates review and undercuts Acceptance Criterion #5 (each adopted change independently revertible).

**Trigger:** Reviewing/reverting the optimization in isolation.

**Actual behavior:** An extra mixed-flux-port Phase-4 test and two Makefile targets are co-mingled.

**Expected behavior:** Only the action_plan_2026-06-24 Phase-0 artifacts in this diff.

**Suggested fix:** Split test_lsw_rk_mixed_flux_bimaterial.cpp and its two Makefile stanzas into a separate commit/branch tied to PLAN_mixed_flux_seissol_port_2026-06-19.md so the optimization diff is self-contained and independently revertible. No code change to either test.

**Test case:** `git status --short` shows both new test files untracked together; the lsw test header cites a different plan.

---

### [R-009] [LOW] [POSSIBLE] [tests/parallel/test_bimaterial_deriv_cache_parity.cpp:194-199] — Part C is a faithful proxy ONLY for the bulk fused-CK corrector; the production fault-corrector consumption of &I + imposed states is not exercised

**Category:** EDGE_CASE (completeness, not soundness)

**Description:** The production `--shared-ck-recursion` path is the driver's AdvanceADERWithSubStep_Spatial (spatial_dyn_driver.cpp:459-514): ComputeADERSubStepStatesAndIntegral → EvaluateBulkAtFaultQPsCanonical per substep → friction iterator → SetSubStepFaultImposedStates → AdvanceADER(Q,...,&shared_I). The test (fault-free, `bc.fault_attr=0`) calls the operator-level ComputeADERSubStepStatesAndIntegral + AdvanceADER(&I) directly with no fault QPs and no SetSubStepFaultImposedStates, so AdvanceADER's fault branch (ComputeADERFaceFluxRHS `if (bc_.fault_attr>0)` at inl:3898, which reads imposed states AND I) is never taken. The test DOES validate the two things that matter for correctness: (1) Q_per_node fused == separate bit-exact (friction-iterator inputs identical), and (2) qn fused (&I) == qn separate (internal recompute) bit-exact on the bulk corrector (I substitution sound). The fault-branch consumption of I is transitively covered because I is proven identical regardless of source — this is a completeness gap, not a soundness hole.

**Trigger:** Fault-free fixture vs fault-bearing production.

**Actual behavior:** Only the bulk corrector consumes &I; the fault branch is skipped.

**Expected behavior:** Ideally a fault-bearing fixture also exercises the ADER fault corrector reading &I + imposed states.

**Suggested fix:** Acceptable as-is given the transitive argument (I and Q_per_node proven bit-identical between fused and separate). If tighter coverage is wanted, add a minimal locked-fault fixture (`bc.fault_attr>0` with a trivial FaultFaceFlux) and repeat Part C so the `if(bc_.fault_attr>0)` branch consumes the precomputed I. Per CLAUDE.md the fault basis needs a non-degenerate fault, which a tiny Cartesian fixture may not provide — document if deferred.

**Test case:** Grep the test for SetSubStepFaultImposedStates / EvaluateBulkAtFaultQPsCanonical: absent. The production fused path's distinctive fault-coupling steps are not reproduced.

---

### [R-010] [LOW] [NOT REACHABLE TODAY] [dynamic/wave_operator.hpp:1082-1083 (mutable flags); inl:3873/4952 (writes in const correctors)] — mutable/const usage is correct for single-threaded MFEM but carries an unstated single-thread assumption

**Category:** QUALITY

**Description:** The two flags are `mutable` and written inside the `const` correctors. This is the textbook legitimate use of `mutable` for a memoized result, and const-correctness is fine — the cached bool feeds only `MFEM_VERIFY`, never flux numerics (so output is byte-preserved). But the read-modify-write `if(!flag){...; flag=true;}` is not atomic and would be a data race if these const methods were ever invoked from multiple threads on the same WaveOperator (MFEM supports shared-memory/OpenMP in some configs). The SEAS ADER path is MPI-only / one-operator-per-rank, so this is benign, but the docstring at hpp:1078 should note the single-thread assumption since a mutable-write-in-const is the classic surprise for a future maintainer who adds threading.

**Trigger:** Future concurrent const access to one WaveOperator instance.

**Actual behavior:** Plain non-atomic mutable bool written in a const method; safe under the current MPI-only, one-operator-per-rank model but undocumented as such.

**Expected behavior:** Document the single-thread-only memoization, or use `std::once_flag`/atomic if concurrent const access is ever possible.

**Suggested fix:** Append to the hpp:1078 docstring: "Single-threaded assumption: these flags are written under the implicit one-WaveOperator-per-rank, no-concurrent-Mult execution model. If correctors are ever called from multiple threads on the same instance, promote to std::once_flag / atomic." No functional change today.

**Test case:** N/A (documentation only).

---

### [R-011] [LOW] [NOT REACHABLE TODAY] [dynamic/wave_operator.hpp:395-399, 409-410 (nullptr-clear re-arm path)] — nullptr-clear re-arm arms a guaranteed next-step abort; correct and intended, but a non-obvious consequence worth a comment

**Category:** EDGE_CASE

**Description:** On `SetAbsorbingBackground(nullptr)`: `has_bulk_bg_=false` (hpp:397) and both flags re-arm false (hpp:409-410). If the clear is uniform across all ranks, the next AdvanceADER's lazy-once consensus computes `min_has==0` on every rank and the `MFEM_VERIFY` at inl:3867/4945 aborts everywhere fail-loud — matching the ORIGINAL per-step behavior (which would also abort once all ranks read `has_bulk_bg_==0`). So the re-arm is CORRECT and behavior is preserved, NOT a regression. The subtle point: the re-arm makes the clear path "arm a guaranteed next-step abort," which is only the intended outcome because production never clears mid-run (the only call is the one-time non-null setter at spatial_dyn_driver.cpp:1999 and the tpv* equivalents). The re-arm correctly distinguishes a re-enable (non-null after a clear) which SHOULD re-verify. The only residual hazard is the non-uniform-clear deadlock already captured in R-002.

**Trigger:** A uniform `SetAbsorbingBackground(nullptr)` clear followed by AdvanceADER (production never does this).

**Actual behavior:** Re-arm flips both flags false on the nullptr path; next step re-runs consensus and aborts with `min_has==0` (uniform clear) — preserved/intended.

**Expected behavior:** Clearing the background leaves the operator in fluctuation-Q mode; callers wanting fluctuation-Q must pass a zero-filled array, not nullptr-clear-then-advance.

**Suggested fix:** Add a comment at hpp:395-399 clarifying that the nullptr-clear path arms a deliberate fail-loud abort on the next AdvanceADER, and that the clear is only well-defined when performed uniformly on all ranks. No code change to the re-arm itself.

**Test case:** Reason about a uniform clear then AdvanceADER: aborts everywhere fail-loud (preserved). A non-uniform clear is the R-002 hazard.

## Summary
- **Critical: 0 / Moderate: 6 (R-001..R-006) / Low: 5 (R-007..R-011)**
- **Plan compliance: PARTIAL.** §0C landed but as a *third* design (lazy-once-at-corrector, two flags) rather than the plan's mandated option (a); the deviation is correct (see Considered and cleared) but the plan text was not updated and the plan-mandated debug-invariance assert (§0C line 103) was not added (R-002). The §0B gate test deviates from the plan's two-material requirement (R-003), runs at the wrong FE order (R-004), is an unregistered orphan (R-005), and leaves the re-arm uncovered (R-006). §0D is deferred (R-007, intentional per task brief). §0A sbatch is correct.
- **Verdict: PASS WITH FIXES.** The core 0C lazy-once hoist is numerically byte-exact (the reduced int feeds only `MFEM_VERIFY`, never flux math) and deadlock-free at scale on the all-ranks local-corrector path; the deviation from plan option (a) is the *correct* call. No new deadlock or numerical divergence is introduced by this change. But R-001 (pre-existing full-comm-Allreduce-behind-`n_shared==0` landmine that the cleanup pass should have removed) and the four test-coverage deviations (R-003..R-006) are real and should be fixed before the gate test is trusted as standing protection and before `--deriv-cache` is adopted for production science.

## Considered and cleared
- **Deadlock-at-scale of the lazy-once hoist (the headline question): CLEARED for the delta.** All four reviewers independently verified the hoist is deadlock-free on the path it actually changes. The local corrector ComputeADERFaceFluxRHS:3860 has **no** `n_shared` gate and runs its consensus unconditionally on every rank each step; `bulk_bg_consensus_local_done_` flips in lockstep (collective + symmetric MIN feeding `MFEM_VERIFY` → uniform abort-or-proceed, with the flag write AFTER the verify). AdvanceADER is called collectively every macro-step by all four drivers (spatial:3343, tpv102:2393/2400, tpv104:2483/2490); branch selectors (is_rk/is_lsw/use_substep_iterator/use_shared_ck) are config-derived rank-uniform and dt is fixed (no per-rank adaptive dt in the spatial loop), so loop count and branch are identical on every rank. `GetNSharedFaces()` is constant across steps (no in-loop Rebalance/repartition/AMR), so the "n_shared 0→>0 mid-run" attack is unreachable. Mult/RK4 use the non-ADER ComputeFaceFluxRHS/ComputeSharedFaceFluxRHS, which never touch the consensus flags, so no non-ADER interleaving can desync them. The two-flag split correctly isolates the shared corrector's `n_shared`-gated participation set from the local corrector's all-ranks set. (The residual full-comm-behind-`n_shared==0` hazard is real but **pre-existing** from R-1505, not a property of the hoist — captured as R-001, not as a hoist deadlock.)
- **Plan option (a) deviation: CLEARED — the deviation is correct and option (a) was itself unsound.** Plan §0C mandated putting the consensus Allreduce INSIDE SetAbsorbingBackground (mirroring SetMixedFluxMode) and reducing the per-step sites to a rank-local tripwire. Two reviewers verified this is the *wrong* design: R-1505 exists to detect a rank that FORGOT to call SetAbsorbingBackground. If the consensus lived in the setter, that very rank would never reach the in-setter Allreduce, defeating detection AND deadlocking the ranks that did call it. The plan's cited template (SetMixedFluxMode) detects a *different* failure (callers passing a different mode) and likewise cannot catch a rank that never calls the setter. The shipped corrector-site lazy-once design preserves R-1505's fail-loud detection. The plan's acceptance text should be updated to record this; no code change is needed for correctness of the deviation itself.
- **mutable/const correctness: CLEARED.** Textbook-legitimate memoization; the cached bool feeds only `MFEM_VERIFY`, never flux numerics, so output is byte-preserved (single-thread caveat noted as R-010).
- **nullptr-clear re-arm logic: CLEARED.** On a uniform clear it preserves the original's fail-loud abort and correctly forces re-verification on re-enable (noted as R-011).
- **BimaterialWaveOperator inheritance: CLEARED.** The production CVM operator does not override SetAbsorbingBackground, so it inherits the base re-arm and mutable members correctly.
- **Makefile object list / API surface / sbatch shell: CLEARED.** seas_test_bimaterial_deriv_cache_parity's link object list exactly mirrors the serial bimaterial parity target; all test-used APIs are public and exist; the sbatch builds OPT_ARGS as a space-separated string with no glob/IFS hazards, word-splitting `${OPT_ARGS}` on the srun line exactly like the sibling `${MIXFLUX_ARG}/${TI_ARG}/${CFL_ARG}` idiom; `SAFS_SHARED_CK`/`SAFS_DERIV_CACHE` default 1 with working A/B `=0` disable.
- **NO-TOUCH respected: CLEARED.** `git diff` shows ZERO changes under bp5/, bp1/, bp2/, domain/, fault/, solver/, friction/ (incl. friction/dieterich_ruina.hpp). All edits are confined to dynamic/wave_operator.{hpp,inl}, the SAFS sbatch, the Makefile, and the new test files.
- **default-on `--deriv-cache` for this run: CLEARED as appropriate (not a footgun).** This is an explicitly throwaway caliper/timing job (tfinal=1s, `--paraview-max-snapshots 5`, `--checkpoint-every 0`) whose purpose is to profile the optimized path; the sbatch header correctly warns that production science must still pass a SAFS Expanse fingerprint check before adopting `--deriv-cache`, matching plan §0B "adopt-with-care."
- **Mixed-flux / fault-free fixture: CLEARED as non-gaps for the cached-kernel claim.** Production defaults SAFS_MIXED_FLUX=none (pure upwind Godunov) and the test matches; the cached ApplySpatialDerivative/ComputeVolumeRHS kernels are element-local with no fault coupling, so a fault-free fixture is acceptable for the deriv-cache claim per se (the genuine gaps are material homogeneity R-003 and FE order R-004, not fault presence).

## Unreviewed Areas
- **Multi-rank scale behavior (R-001, R-002 deadlock reproductions) is only validatable on a real partition with a non-uniform `n_shared` distribution** (island/fully-interior rank), which the production triq SAFS mesh at 2 nodes × 50 tasks = 100 ranks makes unlikely but does not forbid. This requires Expanse, not the laptop; a construction-order np≥3 fixture forcing zero shared faces on one rank would reproduce R-001 locally but was not run here.
- **Production-mesh numerical fingerprint of `--deriv-cache` + `--shared-ck-recursion`** against the SAFS reference (the Expanse fingerprint the sbatch header demands before production adoption) is out of scope for this code review and unrun.
- **The fault-bearing ADER fault corrector consuming &I + imposed states (R-009)** is exercised only transitively; a direct locked-fault fixture was not built (and per CLAUDE.md needs a non-degenerate fault basis a tiny Cartesian mesh may not provide).


---

## Round 2 — Fix Resolution (2026-06-24, post-review)

All MODERATE findings fixed and the doc-only LOW findings addressed; re-tested np2+np4 (all PASS, fused==separate bit-exact, Cached==OnTheFly ≤1.4e-15).

| ID | Sev | Disposition |
|----|-----|-------------|
| R-001 | MOD | **FIXED.** Removed the full-comm Allreduce behind `n_shared==0` in `ComputeADERSharedFaceFluxRHS`; it is now a rank-local `MFEM_VERIFY`. The cross-rank consensus is the gate-free check in `ComputeADERFaceFluxRHS` (confirmed sole caller runs it immediately before, wave_operator.inl:5762→5766). Collapsed two flags → one `bulk_bg_consensus_done_`. Removes the pre-existing R-1600-shape landmine. |
| R-002 | MOD | **FIXED (doc).** Collective-call + setup-only contract now documented on the flag (mirrors R-1205); explained why a per-step rank-local re-check is deliberately NOT added (it would itself be the R-1505 deadlock). |
| R-003 | MOD | **FIXED.** Gate test now uses a TWO-material `FunctionCoefficient` fixture (A in x<0.5, B in x≥0.5) so per-element flux-pool divergence is exercised. |
| R-004 | MOD | **FIXED.** Test sweeps FE order over {1,2}; P1 is the production `--deriv-cache` regime. Both pass ≤1e-12. |
| R-005 | MOD | **FIXED.** `test-bimaterial-deriv-cache-parity` added to the `test-parallel` aggregate; sbatch comment corrected. |
| R-006 | MOD | **FIXED.** Test calls `SetAbsorbingBackground` a second time mid-run and re-runs AdvanceADER, exercising the consensus re-arm. (Coverage is path/no-hang + output-correctness; the flag itself is invisible to output, as the reviewer noted.) |
| R-007 | LOW | **ACCEPTED (deferred).** 0D global-default flip intentionally not done; the caliper sbatch logs `OPT_ARGS`, so the profiling run is self-documenting. Recommended follow-up. |
| R-008 | LOW | **NOT MINE.** `tests/unit/test_lsw_rk_mixed_flux_bimaterial.cpp` is pre-existing untracked WIP from a different plan (mixed-flux port); not part of this change. Left untouched. |
| R-009 | LOW | **ACCEPTED.** Part C covers the bulk fused-CK corrector; the fault corrector's consumption of &I is transitively covered (I proven bit-identical between fused/separate). A fault-bearing fixture needs a non-degenerate fault basis a tiny Cartesian mesh lacks (CLAUDE.md) — deferred. |
| R-010 | LOW | **FIXED (doc).** Single-thread memoization assumption documented on the flag. |
| R-011 | LOW | **FIXED (doc).** nullptr-clear re-arm consequence documented in the setter. |

**Post-fix verdict: PASS.** Core 0C hoist is byte-exact, deadlock-free (now also free of the pre-existing R-001 landmine), and the gate test covers the production P1 + two-material matrix path. Remaining open items are Expanse-only (production fingerprint of `--deriv-cache`; the R-001/R-002 deadlock reproductions need a non-uniform-`n_shared` partition) and the intentionally-deferred 0D.
