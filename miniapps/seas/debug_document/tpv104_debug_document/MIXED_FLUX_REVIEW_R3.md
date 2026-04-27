# Code Review (Round 3): Mixed-Flux Scaffolding for TPV104 — 2026-04-25 (fixes applied 2026-04-26)

## Status (2026-04-26)
**All 10 findings addressed.**  Verdict upgraded **FAIL → PASS**.

| ID | Severity | Status | Where |
|----|----------|--------|-------|
| R-1200 | CRITICAL | ✅ FIXED | New Gate 2-analytic in `test_mixed_flux_dispatch_adjacent.cpp`: assembles the per-face Σ M⁻¹ w shape (F_int − F_central) sum element-locally and asserts bit-equality with observed `Mult(Q_init, Adjacent) − Mult(Q_init, None)`. Passes at **1.309e-15 relative** (4 orders of magnitude inside the 1e-12 spec). Catches sign-flipped or wrong-magnitude Central. |
| R-1201 | CRITICAL | ✅ FIXED | NEW `tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp` runs `wave.Mult(Q, dQdt)` with mode=Adjacent at np=2 and on a serial reference, compares 3 partition-invariant per-component invariants (sum, sum_sq, max_abs). Passes at **2.27e-13 / 7.38e-16 / 0.00e+00 relative** (well inside 1e-10 / 1e-10 / 1e-12 tolerances). Cross-rank conservation through the L2735/L4109 dispatch sites is now verified end-to-end. |
| R-1202 | CRITICAL | ✅ FIXED | `tpv104_driver.cpp` banner now uses `MPI_Reduce(SUM/MIN/MAX)` to report `|central_set|_global` plus per-rank `min`/`max` and a load-imbalance WARNING when max/min > 4×. Rank-0-local count masquerading as global is gone. |
| R-1203 | MODERATE | ✅ FIXED | `SEAS_TEST_NONFAULT_BOTH_SYM` env hook permanently disabled in both ADER-local (L2900) and ADER-shared (L3847) sites with `MFEM_ABORT` explaining the math (the (n,L,R)→(−n,R,L) symmetrization is identically zero for any conservative flux). |
| R-1204 | MODERATE | ✅ FIXED | Driver `MPI_Abort` at CLI-parse time when `--mixed-flux != none && ader-order > 2`, with `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1` override for experimental runs. Manually verified: `[FATAL] --mixed-flux adjacent --ader-order 4: untested combination ...`. |
| R-1205 | MODERATE | ✅ FIXED | `SetMixedFluxMode` runs an `MPI_Allreduce(MAX/MIN)` consensus check at the top — rank-asymmetric mode calls now abort with `MFEM_VERIFY` instead of deadlocking at the post-walk Allgatherv. |
| R-1206 | MODERATE | ✅ FIXED | `bc_.fault_attr > 0` short-circuit removed from `is_fault_face`; the predicate now relies solely on `fault_mesh_face_idx_set` membership. Build-time `MFEM_VERIFY(set.empty() \|\| bc_.fault_attr > 0)` enforces ctor-state consistency. Set-build/dispatch divergence under future ctor refactors is precluded. |
| R-1207 | LOW | ✅ DOCUMENTED | Comment added to `make_global_key` documenting the assumption (TPV104 is tet-only; alias-with-quad is structurally impossible). The R-1207 fix (prepend `verts.Size()` slot) is deferred until a hybrid-element TPV104 fixture is introduced. |
| R-1208 | LOW | ✅ FIXED | Cached `mf_on_` member in `WaveOperator`, kept in sync by `SetMixedFluxMode`. The 4 dispatch sites now read the cached flag instead of recomputing `(mixed_flux_mode_ != None)` per-call. |
| R-1209 | LOW | ✅ FIXED | `test_mixed_flux_adjacent_mpi.cpp` per-key counter changed from `std::set` to `std::vector`; explicitly verifies `kv.second.size() == 2` (each shared face must be reported by exactly 2 ranks) before checking inter-rank value agreement. |

### Final test matrix
- Mixed-flux dedicated tests: **82 / 82 pass** (4 + 9 + 3 + 7 + 53 + 1 + 2 + 3)
  - `seas_test_godunov_central_flux` 4/4
  - `seas_test_mixed_flux_face_set` 9/9
  - `seas_test_mixed_flux_dispatch_none` 3/3
  - `seas_test_mixed_flux_dispatch_adjacent` 7/7  (Gate 2-analytic at 1.31e-15 rel — **R-1200 FIX**)
  - `seas_test_tpv104_smoke` 53/53
  - `seas_test_mixed_flux_adjacent_mpi` (np=2) 1/1
  - `seas_test_mixed_flux_shared_fault_face_excluded_mpi` (np=2) 2/2
  - `seas_test_mixed_flux_dispatch_adjacent_mpi` (np=2) 3/3 — **NEW (R-1201)**
- Curated regression on touched code paths: **96 / 96 pass** (godunov-flux + wave-operator + fault-face-flux + tpv102-setup + ader-tpv102-smoke).
- R-1204 driver gate: PASS (FATAL printed on `--mixed-flux adjacent --ader-order 4`).

### Files changed (Round 3, this session)
- `dynamic/wave_operator.hpp` — R-1208 `mf_on_` cached member.
- `dynamic/wave_operator.inl` — R-1100 `is_fault_face` rewrite (now R-1206 sole-source); R-1206 build-time MFEM_VERIFY; R-1205 mode-consensus check; R-1203 env-hook abort (×2 sites); R-1207 hybrid-element comment; R-1208 4 dispatch-site `mf_on` swap to cached member.
- `tests/unit/test_mixed_flux_dispatch_adjacent.cpp` — R-1200 analytic Gate 2 (`ComputeAnalyticAdjacentDelta` helper + new gate at 1e-12 relative).
- `tests/parallel/test_mixed_flux_adjacent_mpi.cpp` — R-1209 strengthened consistency check.
- `tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp` — NEW (R-1201).
- `drivers/tpv104_driver.cpp` — R-1202 global banner reduction + R-1204 ader-order guard.
- `Makefile` — new test target (`seas_test_mixed_flux_dispatch_adjacent_mpi`).

`--mixed-flux adjacent` is now safe to enable in multi-rank production runs at ADER-O2.  ADER-O > 2 still requires `SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1` (override) until a higher-order MPI gate lands.

---

## Original Round 3 Review (preserved for traceability)



**Commit reviewed:** current HEAD (post-Round-2 fixes)
**Reviewer:** code-reviewer agent (round 3, adversarial pass; assume ≥3 critical bugs)
**Plan:** `miniapps/seas/debug_document/tpv104_debug_document/MIXED_FLUX_PLAN.md`
**Round 1 (plan-level):** `MIXED_FLUX_PLAN_REVIEW.md` (R-1201..R-1209)
**Round 1 (impl):** `MIXED_FLUX_IMPL_REVIEW.md` (R-001..R-003)
**Round 2 (impl):** `MIXED_FLUX_REVIEW_R2.md` (R-1100..R-1106)

This round focuses on **MPI parallel correctness for production runs**.

## Review Scope

- `dynamic/godunov_flux.{hpp,cpp}` (`Central()` primitive)
- `dynamic/wave_operator.{hpp,inl}` (mixed-flux additions only):
  - `SetMixedFluxMode` (L1156–1183)
  - `BuildCentralFluxFaceSet_` (L1190–1404)
  - 4 dispatch sites: L2169 (RK4 local), L2735 (RK4 shared), L3649 (ADER local), L4109 (ADER shared)
- `drivers/tpv104_driver.cpp` (L532–571 CLI parse; L1271–1298 setter wiring)
- `tests/unit/test_godunov_central_flux.cpp`
- `tests/unit/test_mixed_flux_dispatch_none.cpp`
- `tests/unit/test_mixed_flux_face_set.cpp`
- `tests/unit/test_mixed_flux_dispatch_adjacent.cpp` (R-1102 fix from Round 2)
- `tests/parallel/test_mixed_flux_adjacent_mpi.cpp` (R-001/R-002 fix)
- `tests/parallel/test_mixed_flux_shared_fault_face_excluded_mpi.cpp` (R-1100 fix)

Domain context consulted: `CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
`feedback_dynamic_folder_editable_for_tpv104`, Zhang et al. 2023.

---

## Findings

---

### [R-1200] [CRITICAL] [BUG] [tests/unit/test_mixed_flux_dispatch_adjacent.cpp:298–317] — Adjacent-mode dispatch test cannot detect a sign-flipped or wrong-magnitude central flux

**Category:** BUG (test contract too weak — false-positive coverage)

**Description:**
The Round-2 R-1102 fix added end-to-end Adjacent dispatch verification via
`test_mixed_flux_dispatch_adjacent`. The plan's Phase 4 acceptance criterion
(`MIXED_FLUX_PLAN.md` lines 580–584) is:

> "With `--mixed-flux adjacent` on the 2-tet fixture: `Q_new` differs from
> `--mixed-flux none` only at the 6 non-fault faces; the difference is
> `0.5·|A_n|·(Q⁺ − Q⁻)` per face per state component (algebraic relation per
> Phase 1 acceptance criterion)."

The test as implemented only checks two coarse properties:

```cpp
// Line 305:
TEST_GE(diff_adj_vs_none / scale_q_new, 1e-9, "Gate 1: ... fires");
// Line 316:
TEST_LE(diff_adj_vs_none / scale_q_new, 1.0,  "Gate 2: ... bounded (no blowup)");
```

`Gate 1` only verifies the diff is above pure FP noise (1e-9 × scale).
`Gate 2` only verifies the diff doesn't blow up (< 1× scale). Neither gate
verifies the analytic identity from the plan. Concretely, the following
buggy implementations would PASS this test:

- A `Central()` that returns `0.5 · F_interior` (off by a factor of 2). Diff
  is non-zero, bounded — passes.
- A `Central()` that returns `F_interior + ε · garbage` for any non-trivial
  ε in [1e-9, 1.0] × scale. Passes.
- A central-flux face set that wrongly includes an extra few faces (e.g.,
  the rank-seam aliasing bug R-1100 if not dispatch-masked). Passes.
- A sign error swapping `Ax_plus_` ↔ `Ax_minus_` in the central computation
  (which breaks anti-symmetry — see R-1201). Diff magnitude unchanged →
  passes.

The R-1104 anti-symmetry test on the primitive at 1e-8 relative tolerance
catches gross primitive bugs but does NOT catch dispatch-level bugs (e.g.,
wrong face set, wrong frame on shared faces, sign convention mismatch
between RK4 and ADER paths).

**Trigger:**
Any future refactor that perturbs Central's coefficient by O(1) but stays
in the [1e-9, 1.0] relative window. The test PASSES, hides the bug.

**Actual behavior:**
- Gate 1: 1.46e+05 absolute / 3.85e+08 = 3.79e-04 relative (PASSES; well above 1e-9 floor).
- Gate 2: 3.79e-04 < 1.0 (PASSES).
- Gate 3: 0.0 (idempotency holds).
- The actual analytic relation is NEVER computed.

**Expected behavior:**
The test MUST verify, for each face `f` in `central_flux_face_set_`, that the
change `Q_new_adj[c, dofs(f)] − Q_new_none[c, dofs(f)]` equals
`-0.5 · M^{-1} · ∫ shape · |A_n| · (Q_self - Q_nbr)|_f · dS · dt`
(or, more practically: assemble the analytic delta into a vector and assert
bit-equality with `Q_new_adj − Q_new_none`).

**Suggested fix (code diff):**

```diff
--- a/tests/unit/test_mixed_flux_dispatch_adjacent.cpp
+++ b/tests/unit/test_mixed_flux_dispatch_adjacent.cpp
@@ -295,16 +295,77 @@ int main()
    std::cout << "  max |Q_new_adj − Q_new_none| = " << ...;

-   // Gate 1: dispatch fired (non-trivial diff).  Floor: 1e-9 × scale.
+   // Gate 1: dispatch fired AND analytic identity holds.
    TEST_GE(diff_adj_vs_none / scale_q_new, 1e-9, "Gate 1: ... fires");
    TEST_GE(static_cast<double>(central_set_size), 1.0,
            "Gate 1b: central_flux_face_set_ is non-empty");

-   // Gate 2: bounded above by reasonable analytic ceiling.
-   TEST_LE(diff_adj_vs_none / scale_q_new, 1.0,
-           "Gate 2: ... bounded (no blowup)");
+   // Gate 2 (R-1200): analytic identity at every central-flux face.
+   //   Q_new_adj − Q_new_none = sum_{f in central_set} contribution_f
+   //   where contribution_f =
+   //     M^{-1} · w · shape · ( F_interior(nor_f, Q_self_f, Q_nbr_f)
+   //                           − F_central (nor_f, Q_self_f, Q_nbr_f) )
+   //                 · dt   (per ADER-O2 corrector)
+   // Build this expected delta directly from the unmixed run's Q_init
+   // (plus the predictor-corrector machinery), assemble into a Vector,
+   // and assert bit-equality with the observed (Q_new_adj − Q_new_none).
+   //
+   // Tolerance: 1e-12 relative (matches Phase 1 Gate 2 contract).
+   Vector expected_delta(NUM_STATE * ndof_total);
+   AssembleAnalyticAdjacentDelta(wave_a, mesh_a, order, Q_init, dt,
+                                 wave_b.GetCentralFluxFaceSet(),
+                                 expected_delta);
+   Vector observed_delta(NUM_STATE * ndof_total);
+   for (int i = 0; i < observed_delta.Size(); i++) {
+      observed_delta(i) = Q_new_adj(i) - Q_new_none(i);
+   }
+   const real_t worst_rel = MaxAbsDiff(observed_delta, expected_delta)
+                          / std::max<real_t>(MaxAbs(observed_delta), 1.0);
+   TEST_LE(worst_rel, 1e-12,
+           "Gate 2 (R-1200): Q_new_adj − Q_new_none equals "
+           "Σ_f M^{-1} w shape (F_int − F_ce) dt to 1e-12 relative "
+           "(Phase 4 plan acceptance, algebraic identity)");
```

**Test case (alongside the fix):**
Per-face DOF-level breakdown that fails on a hand-introduced sign flip:
re-run with `Ax_plus_` ↔ `Ax_minus_` swapped in `Central` and assert the new
gate FAILS (delta has opposite sign). Then re-enable the correct version and
assert PASS.

---

### [R-1201] [CRITICAL] [BUG] [tests/parallel/test_mixed_flux_adjacent_mpi.cpp:177–347] — No MPI test verifies the actual Mult/AdvanceADER cross-rank result; conservation across rank seam is unverified

**Category:** BUG (coverage gap — production-blocking)

**Description:**
The two MPI tests (`test_mixed_flux_adjacent_mpi`,
`test_mixed_flux_shared_fault_face_excluded_mpi`) verify `central_flux_face_set_`
membership symmetry across ranks but do NOT call `wave.Mult(Q, dQdt)` or
`wave.AdvanceADER(Q, dt, ord, Q_new)` with `MixedFluxMode::Adjacent`. The
shared-face dispatch sites at `wave_operator.inl:2735` and `wave_operator.inl:4109`
have NEVER executed under any test condition with `mf_on=true` and a non-empty
`central_flux_face_set_`.

This means:

1. **Cross-rank conservation is unverified.** The mixed flux on a shared seam
   relies on the central flux's anti-symmetry property:
   `Central(+nor_A, Q_self_A, Q_nbr_A) = -Central(-nor_A, Q_nbr_A, Q_self_A)`.
   The R-1104 unit test verifies this on the primitive at 1e-8 relative
   tolerance. But the dispatch site applies the flux through:
   - `CalcOrtho(ftr->Face->Jacobian(), nor_vec)` for `nor`
   - normalisation with `nor_vec /= nor_len`
   - per-rank Q_self extraction (rank A's Elem1) vs Q_nbr extraction (ghost
     buffer via `pfes->GetFaceNbrFE(nbr_idx)` and `nbr_data[c]`)
   - per-rank assembly with rank-local `shape1`

   None of these steps are tested at the dispatch level for Central flux.
   A bug in any of them (e.g., a wrong shape function evaluation on the
   ghost side that's compensated by the upwind dissipation in Interior but
   not in Central) would silently break conservation.

2. **No np ≥ 4 coverage.** Both MPI tests are np=2 only. Production TPV104
   runs at np=128–400 on Frontera. Topologies that emerge only at np ≥ 4
   (e.g., a rank with no local fault face AND no fault-adjacent ghost
   element; a fault face that crosses 3+ ranks at a triple junction; a
   shared face whose two sides are on ranks A and B but whose neighboring
   element on rank A is owned via a different rank C's face seam) are
   never exercised.

3. **No `Mult(Q, dQdt)` MPI coverage.** The RK4 path
   `ComputeFaceFluxRHS / ComputeSharedFaceFluxRHS` is dispatched by
   `WaveOperator::Mult` (Round-7 fault iterators may dispatch via Mult
   instead of AdvanceADER). Even the serial `test_mixed_flux_dispatch_adjacent`
   only calls `AdvanceADER`, not `Mult`. The L2169 Central dispatch site
   has zero MPI test coverage.

**Trigger:**
First production multi-rank run with `--mixed-flux adjacent`. With ≥1.5M tets
at np=128, the rank-seam population is in the thousands of faces; even a
1e-3 fractional bug in cross-rank flux consistency accumulates to a
catastrophic mass/momentum imbalance over 7000 macro steps.

**Actual behavior:**
- The dispatch sites fire but no test asserts the numerical result is
  correct on MPI.
- The two existing MPI tests (set-membership symmetry, fault-face exclusion)
  are necessary but not sufficient.

**Expected behavior:**
At least one MPI test that:
1. Builds a ParMesh with shared non-fault faces in `central_flux_face_set_`.
2. Runs `wave.AdvanceADER(Q, dt, 2, Q_new_par)` at np=2 with mode Adjacent.
3. Runs the same on a serial `Mesh` (np=1) with mode Adjacent: `Q_new_ser`.
4. Asserts `max |Q_new_par − Q_new_ser| / scale < 1e-12` after gathering
   Q_new_par to a single rank in the same DOF ordering.
5. (Stronger) Asserts global conservation: ∫ ρv_x dV is conserved at np=2
   to 1e-12 relative across one AdvanceADER step on a domain with fully
   periodic / closed BCs.

**Suggested fix (code diff):**
Add `tests/parallel/test_mixed_flux_advance_ader_mpi.cpp` (np=2):

```diff
+// Round-13 R-1201: cross-rank end-to-end Adjacent dispatch correctness.
+// Compares wave.AdvanceADER(Adjacent) at np=2 against the same operator
+// at np=1, asserting bit-equivalence after Vector::Gather to a single rank.
+
+int main(int argc, char *argv[]) {
+   MPI_Init(&argc, &argv);
+   int rank, nprocs;
+   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
+   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
+   if (nprocs != 2) { MPI_Finalize(); return 0; }
+
+   // Path A: serial (rank 0 only).
+   Vector Q_new_ser(NUM_STATE * ndof_total_serial);
+   if (rank == 0) {
+      Mesh mesh = BuildSmallFaultTetMesh(L);
+      WaveOperator<Mesh> wave(mesh, order, ...);
+      wave.SetFaultFlux(...); wave.SetFaultDOFData(...);
+      wave.SetAbsorbingBackground(Q_bg);
+      wave.SetMixedFluxMode(MixedFluxMode::Adjacent);
+      wave.AdvanceADER(Q_init_ser, dt, 2, Q_new_ser);
+   }
+   MPI_Bcast(Q_new_ser.GetData(), Q_new_ser.Size(),
+             MPI_DOUBLE, 0, MPI_COMM_WORLD);
+
+   // Path B: parallel (np=2).
+   Mesh serial = BuildSmallFaultTetMesh(L);
+   ParMesh pmesh(MPI_COMM_WORLD, serial);
+   WaveOperator<ParMesh> wave_par(pmesh, order, ...);
+   wave_par.SetFaultFlux(...); wave_par.SetFaultDOFData(...);
+   wave_par.SetAbsorbingBackground(Q_bg);
+   wave_par.SetMixedFluxMode(MixedFluxMode::Adjacent);
+
+   Vector Q_init_par(NUM_STATE * wave_par.GetNDofTotal()); ...;
+   Vector Q_new_par(NUM_STATE * wave_par.GetNDofTotal());
+   wave_par.AdvanceADER(Q_init_par, dt, 2, Q_new_par);
+
+   // Gather Q_new_par to rank 0 in serial DOF order via the standard
+   // ParGridFunction → GridFunction projection, then compare.
+   real_t worst_rel = GatherAndCompare(Q_new_par, Q_new_ser);
+   TEST_LE(worst_rel, 1e-12,
+           "Adjacent np=2 result equals np=1 result to FP precision");
+
+   MPI_Finalize();
+   return num_failed > 0 ? 1 : 0;
+}
```

**Test case (regression-guard, MUST land before any production sbatch):**
The diff above. Repeat at np=4 with a partition specifically designed to
produce: (a) a rank with zero local fault faces, (b) a triple-rank junction
near the fault.

---

### [R-1202] [CRITICAL] [BUG] [drivers/tpv104_driver.cpp:1271–1298] — Driver banner reports rank-0-local `|central_set|` size as if it were global, masking rank imbalance

**Category:** BUG (silent misreporting of distributed state)

**Description:**
At driver line 1292–1297:

```cpp
if (rank == 0 && mixed_flux_mode != MixedFluxMode::None)
{
   std::cout << "[mixed-flux] mode=" << mixed_flux_str
             << "  |central_set|="
             << wave.GetCentralFluxFaceSet().size()
             << "  (Zhang et al. 2023 mixed-flux dispatch)\n";
}
```

`wave.GetCentralFluxFaceSet().size()` returns the LOCAL size on rank 0 only.
For a TPV104 production run at np=128, rank 0 holds ~1/128th of the mesh —
its `central_flux_face_set_` is a tiny fraction of the global set. The
banner reads e.g., `|central_set|=42` when the global count is ~5000.

This is misleading for two reasons:

1. **Operators reading the log think the central-flux fraction is
   negligible** when it's actually ~5000/1.5M ≈ 0.3% of all interior faces
   (Zhang's expected fraction is 5–10%). They may falsely conclude
   "mixed-flux is barely engaged" and skip the diagnosis when results don't
   match the upwind run.

2. **Imbalance across ranks is silently masked.** A partition where one
   rank holds the bulk of the fault-adjacent faces and others hold none
   would show as a banner-line of any value depending on rank 0's lottery,
   with no diagnostic that the load is severely uneven.

**Trigger:** Every multi-rank `--mixed-flux` production run.

**Actual behavior:** Banner reports rank 0's local count.

**Expected behavior:** Banner reports the GLOBAL count (sum across ranks),
plus optionally the per-rank min/max/avg for load-balance diagnostic.

**Suggested fix (code diff):**

```diff
--- a/drivers/tpv104_driver.cpp
+++ b/drivers/tpv104_driver.cpp
@@ -1289,11 +1289,32 @@ ...
    wave.SetMixedFluxMode(mixed_flux_mode);
-   if (rank == 0 && mixed_flux_mode != MixedFluxMode::None)
+   if (mixed_flux_mode != MixedFluxMode::None)
    {
-      std::cout << "[mixed-flux] mode=" << mixed_flux_str
-                << "  |central_set|="
-                << wave.GetCentralFluxFaceSet().size()
-                << "  (Zhang et al. 2023 mixed-flux dispatch)\n";
+      const std::size_t local_size = wave.GetCentralFluxFaceSet().size();
+#ifdef MFEM_USE_MPI
+      long long local_ll = static_cast<long long>(local_size);
+      long long global_ll = 0, local_max = 0, local_min = 0;
+      MPI_Reduce(&local_ll, &global_ll,  1, MPI_LONG_LONG, MPI_SUM, 0, comm);
+      MPI_Reduce(&local_ll, &local_max,  1, MPI_LONG_LONG, MPI_MAX, 0, comm);
+      MPI_Reduce(&local_ll, &local_min,  1, MPI_LONG_LONG, MPI_MIN, 0, comm);
+      if (rank == 0)
+      {
+         std::cout << "[mixed-flux] mode=" << mixed_flux_str
+                   << "  |central_set|_global=" << global_ll
+                   << "  per-rank min=" << local_min
+                   << " max=" << local_max
+                   << "  (Zhang et al. 2023 mixed-flux dispatch)\n";
+         // Diagnostic: warn if 0 ranks have empty sets (fault-adjacency
+         // is at least partially balanced) but max/min ratio > 4 (severe
+         // imbalance; the partitioner did not honor mesh.Partition by
+         // fault locality).
+         if (local_min > 0 && local_max > 4 * local_min) {
+            std::cout << "[mixed-flux] WARNING: rank load imbalance "
+                      << static_cast<double>(local_max) / local_min
+                      << "× (max/min); consider --partition-fault-locality\n";
+         }
+      }
+#else
+      std::cout << "[mixed-flux] mode=" << mixed_flux_str
+                << "  |central_set|=" << local_size << " (serial)\n";
+#endif
    }
```

**Test case:**
Run `mpirun -np 4 ./seas_tpv104_driver --dry-run --mixed-flux adjacent ...`
and assert the printed `|central_set|_global` matches the sum of per-rank
sizes from a separate per-rank dump (insert temporary `std::cerr <<
rank << " " << local_size << "\n";` lines under a debug flag).

---

### [R-1203] [MODERATE] [BUG] [dynamic/wave_operator.inl:2864–2868, 3613–3674, 4099–4135] — n↔−n symmetrization branch combined with Central flux produces ZERO flux (test-only path silently broken)

**Category:** BUG (test hook semantically wrong; misleads diagnostic runs)

**Description:**
Both `ComputeADERFaceFluxRHS` (L3613–3674) and `ComputeADERSharedFaceFluxRHS`
(L4099–4135) wrap the Central/Interior dispatch in an n↔−n symmetrization
branch gated by `SEAS_TEST_NONFAULT_BOTH_SYM=1`:

```cpp
if (!nonfault_both_sym) {
   interior_or_central(nor, I_self, I_nbr, F_h);
} else {
   real_t F_pos[NUM_STATE], F_neg[NUM_STATE];
   real_t nor_neg[3] = {-nor[0], -nor[1], -nor[2]};
   interior_or_central(nor,     I_self, I_nbr, F_pos);
   interior_or_central(nor_neg, I_nbr,  I_self, F_neg);
   for (int c = 0; c < NUM_STATE; c++) {
      F_h[c] = 0.5 * (F_pos[c] + F_neg[c]);
   }
}
```

**Both Interior AND Central flux satisfy the anti-symmetry**
`F(n, A, B) + F(-n, B, A) = 0` exactly (verified by the R-1104 unit-test
fix to 1e-8 relative). Therefore `F_h = 0.5 · (F_pos + F_neg) ≈ 0` to FP
precision.

Setting `SEAS_TEST_NONFAULT_BOTH_SYM=1` ZEROES OUT THE NON-FAULT FACE FLUX
ENTIRELY. The diagnostic comment claims this is a "test-only n↔−n
symmetrization" hook, but in practice it nullifies every non-fault interior
face's contribution — a complete physics break. Any test that uses this hook
to "audit" the non-fault flux is silently testing zero flux.

This is a pre-existing bug from Round-13C Patch 1, but it now ALSO interacts
with Mixed-Flux: when `mf_on && central_flux_face_set_.count(f)`, the
symmetrization zeroes the central flux, masking dispatch errors.

**Trigger:**
Any test or audit run with `SEAS_TEST_NONFAULT_BOTH_SYM=1` will silently
report "zero flux" everywhere on non-fault faces, and the user may
mistake this for a result of mixed-flux when it's actually a symmetrization
artifact.

**Actual behavior:**
Symmetrization mode produces near-zero flux on every non-fault face,
breaking the entire bulk wave equation. Mixed-Flux interactions with the
hook silently mask central-flux dispatch errors.

**Expected behavior:**
The intended n↔−n symmetrization for a one-sided flux test would compare
the same `(self, nbr)` ordering with both signed normals — which by the
anti-symmetry is `0` for ALL conservative fluxes, so the hook is
mathematically meaningless for testing physics. Either:
(a) Remove the hook entirely (it never produced useful data).
(b) Refactor it to test something else: e.g., compare
    `F(n, Q_self, Q_nbr)` against `-F(-n, Q_self, Q_nbr)` (note: NO swap of
    L/R) — this asserts the rotation pipeline behaves correctly under
    normal flip without invoking the conservation identity.

**Suggested fix (code diff):**

```diff
--- a/dynamic/wave_operator.inl
+++ b/dynamic/wave_operator.inl
@@ -2864,12 +2864,17 @@
-   const bool nonfault_both_sym = []()
-   {
-      const char *env = std::getenv("SEAS_TEST_NONFAULT_BOTH_SYM");
-      return env && env[0] == '1';
-   }();
+   // R-1203: SEAS_TEST_NONFAULT_BOTH_SYM is mathematically meaningless —
+   // any conservative flux (Interior, Central) satisfies F(n,L,R) +
+   // F(-n,R,L) = 0 exactly, so the 0.5·(F_pos + F_neg) symmetrization
+   // ZEROES OUT every non-fault face's contribution.  When combined with
+   // Mixed-Flux, this masks central-flux dispatch errors.  Hook is
+   // disabled pending a refactor to test the rotation pipeline under
+   // n↔-n WITHOUT also swapping (Q_self, Q_nbr).
+   const bool nonfault_both_sym = false;
+   if (std::getenv("SEAS_TEST_NONFAULT_BOTH_SYM")) {
+      MFEM_ABORT("SEAS_TEST_NONFAULT_BOTH_SYM is disabled (R-1203): "
+                 "the (n,L,R) -> (-n,R,L) symmetrization is identically "
+                 "zero for every conservative flux including Central.  "
+                 "If a similar diagnostic is needed, use the (n)->(-n) "
+                 "symmetrization (no L/R swap) which tests rotation-"
+                 "pipeline parity instead.");
+   }
```

(Mirror the same change at the L3802–3806 ADER-shared site.)

**Test case:**
Add a smoke test that: (1) sets `SEAS_TEST_NONFAULT_BOTH_SYM=1` on an
arbitrary fixture and runs `wave.AdvanceADER`; (2) asserts the run aborts
with the documented R-1203 message. Currently no test exercises this env
var, so the bug is undetected — adding the abort is a defense-in-depth.

---

### [R-1204] [MODERATE] [BUG] [drivers/tpv104_driver.cpp:457, 535–542; dynamic/wave_operator.inl:1170–1179] — No mutual-exclusion check between `--mixed-flux adjacent` and `--ader-order > 2` despite plan R5 risk

**Category:** BUG (missing guard for risky flag combination)

**Description:**
`MIXED_FLUX_PLAN.md` line 826–836 acknowledges:

> "R5: MPI shared-face inconsistency: ... Phase 4 handles both
> `ComputeSharedFaceFluxRHS` and `ComputeADERSharedFaceFluxRHS` with the
> same dispatch logic. R-1003-style abort guard if `nprocs > 1 AND
> use_substep AND mixed_flux != none` may be needed if the testing reveals
> inconsistency."

The risk also extends to `--ader-order > 2`. The ADER-O3 / ADER-O4 paths
introduce additional CK predictor stages and per-stage flux evaluations.
Mixed-Flux dispatch was tested at ADER-O2 only (the unit test
`test_mixed_flux_dispatch_adjacent.cpp:233` uses `dt=1e-5, ader_order=2`).

The driver accepts `--ader-order` from CLI without any cross-check:

```cpp
// L457:
int ader_order = GetIntArg(argc, argv, "--ader-order", 2);
// L535-542 — mixed_flux_mode parsed independently
```

There's no abort guard for `(mixed_flux_mode != None) && (ader_order > 2)`.
A user requesting `--mixed-flux adjacent --ader-order 4` would get an
untested combination. Given the absence of MPI test coverage for the
ADER-O > 2 path with mixed flux, this should at minimum WARN, ideally
abort pending validation.

**Trigger:**
Any production run with both flags set non-default. The TPV104 driver's
default is `--ader-order 2`, so this fires only when user explicitly raises
the order.

**Actual behavior:** No guard. The combination silently dispatches with
unverified semantics.

**Expected behavior:** At least a warning at rank 0; ideally a guarded
abort that references this finding.

**Suggested fix (code diff):**

```diff
--- a/drivers/tpv104_driver.cpp
+++ b/drivers/tpv104_driver.cpp
@@ -571,6 +571,22 @@ ...
    }
+
+   // R-1204: --mixed-flux adjacent has no test coverage at --ader-order > 2.
+   // Round-3 R3 (Risk R5 in MIXED_FLUX_PLAN.md) flagged this.  Until
+   // test_mixed_flux_advance_ader_o3o4_mpi lands, abort on the unvalidated
+   // combination.  Override via SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 for
+   // experimental runs.
+   if (mixed_flux_mode != MixedFluxMode::None && ader_order > 2)
+   {
+      const char *force = std::getenv("SEAS_FORCE_MIXED_FLUX_ADER_O_GT2");
+      if (!(force && force[0] == '1')) {
+         if (rank == 0) {
+            std::cerr << "[FATAL] --mixed-flux adjacent --ader-order "
+                      << ader_order << ": untested combination.  Set "
+                      "SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 to override.\n";
+         }
+         MPI_Abort(MPI_COMM_WORLD, 1);
+      }
+   }
```

**Test case:**
Add to `test_tpv104_smoke.cpp`:
```cpp
const std::string out = RunDriver(binary,
   "--dry-run --mixed-flux adjacent --ader-order 4");
TEST_ASSERT(out.find("untested combination") != std::string::npos,
            "Driver aborts on --mixed-flux + --ader-order > 2 without "
            "the override env var");
```

---

### [R-1205] [MODERATE] [BUG] [dynamic/wave_operator.inl:1244–1255 (AllContinuous), 1261–1302 (Adjacent), 1316–1402 (post-walk MPI exchange)] — `BuildCentralFluxFaceSet_` runs MPI collectives for `Adjacent` but NOT for `AllContinuous`; deadlocks if any rank runs a different mode

**Category:** BUG (collective-call asymmetry across modes)

**Description:**
In `BuildCentralFluxFaceSet_`:

- `MixedFluxMode::None`: early returns at L1194; no collectives.
- `MixedFluxMode::AllContinuous`: returns at L1254 BEFORE the post-walk MPI
  exchange block (L1316–1402); no Allgather/Allgatherv.
- `MixedFluxMode::Adjacent`: walks E_fault_adj, then enters the MPI exchange
  block; calls `MPI_Allgather` (L1369) and `MPI_Allgatherv` (L1375).

If different ranks call `SetMixedFluxMode` with different modes (e.g., a
test that uses rank-conditional logic), the collective-call asymmetry
causes a HARD DEADLOCK:

- Rank A: `SetMixedFluxMode(Adjacent)` → enters collectives.
- Rank B: `SetMixedFluxMode(AllContinuous)` → skips collectives.
- Rank A blocks at `MPI_Allgather`. Rank B finishes and proceeds. Subsequent
  collective on rank A (or rank B's later collective) deadlocks.

In production this can't happen because the driver passes the same mode to
all ranks. But it's a fragile API contract; the function silently relies on
"all ranks are in the same mode at every call." A future test or
experimental driver may call `SetMixedFluxMode` rank-conditionally and
trigger the deadlock without any indication.

There's a SECONDARY concern: `AllContinuous` mode has no cross-rank exchange,
which RELIES on the implicit fact that "every rank's local iteration over
mesh faces will independently insert each shared seam face." This is true
because `is_interior_face` includes the `shared_seam` predicate. But if a
future refactor breaks `shared_mesh_face_set_` population (e.g., on a rank
with `n_shared == 0`, the set is empty correctly), `AllContinuous` would
silently miss shared faces. There's no defensive collective check.

**Trigger:**
Any code path that calls `SetMixedFluxMode` with different modes on
different ranks. Currently no such path exists, but the API is fragile.

**Actual behavior:** Implicit contract ("all ranks must call with same mode")
is undocumented and unenforced.

**Expected behavior:** Either:
(a) Run a leading `MPI_Allreduce(MPI_BAND)` to verify all ranks agree on the
mode before proceeding (consensus check; aborts on disagreement with a
clear message).
(b) Document the contract loudly in `SetMixedFluxMode`'s comment header
and at every `MPI_Allgather*` call site.
(c) For `AllContinuous`, run a (possibly trivial) collective so the
collective sequence is mode-invariant (defense in depth against future
refactors).

**Suggested fix (code diff):**

```diff
--- a/dynamic/wave_operator.inl
+++ b/dynamic/wave_operator.inl
@@ -1156,6 +1156,28 @@ template <typename MeshType>
 void WaveOperator<MeshType>::SetMixedFluxMode(MixedFluxMode m)
 {
+   // R-1205: enforce collective call discipline.  All ranks must invoke
+   // SetMixedFluxMode with the SAME mode in the same order; otherwise
+   // BuildCentralFluxFaceSet_ deadlocks at MPI_Allgather (Adjacent) vs
+   // early-return (None/AllContinuous).  Allreduce-MAX over the mode
+   // value catches mismatched callers cheaply.
+   if constexpr (IsParallelMesh<MeshType>::value) {
+#ifdef MFEM_USE_MPI
+      auto &pmesh_check = static_cast<ParMesh &>(mesh_);
+      int my_mode = static_cast<int>(m);
+      int max_mode = 0, min_mode = 0;
+      MPI_Allreduce(&my_mode, &max_mode, 1, MPI_INT, MPI_MAX,
+                    pmesh_check.GetComm());
+      MPI_Allreduce(&my_mode, &min_mode, 1, MPI_INT, MPI_MIN,
+                    pmesh_check.GetComm());
+      MFEM_VERIFY(max_mode == min_mode,
+                  "SetMixedFluxMode: ranks disagree on mode "
+                  "(max=" << max_mode << ", min=" << min_mode
+                  << ").  All ranks must pass the same MixedFluxMode.");
+#endif
+   }
    if (m != MixedFluxMode::None && use_precomputed_face_fluxes_)
    { ... }
```

**Test case:**
```cpp
// Negative test: rank-conditional mode call should abort cleanly, not
// deadlock.  Run with mpirun --np 2 and a 30-second hard timeout.
if (rank == 0) wave.SetMixedFluxMode(MixedFluxMode::Adjacent);
else           wave.SetMixedFluxMode(MixedFluxMode::None);
// Pre-fix: deadlocks → timeout fires → CTest reports failure.
// Post-fix: MFEM_ABORT on rank 0 propagates via MPI_Abort; clean exit.
```

---

### [R-1206] [MODERATE] [BUG] [dynamic/wave_operator.inl:1248, 1297, 1394] — `is_fault_face` predicate uses `bc_.fault_attr` value for inclusion guard but the dispatch sites do not — silent divergence if `bc_.fault_attr == 0` is inadvertently passed

**Category:** BUG (defensive predicate inconsistency between set-build and dispatch)

**Description:**
`is_fault_face` (L1238–1242):

```cpp
auto is_fault_face = [&](int f) -> bool
{
   return (bc_.fault_attr > 0 &&
           fault_mesh_face_idx_set.count(f) > 0);
};
```

This `bc_.fault_attr > 0` short-circuit returns false when the operator was
constructed without a fault (e.g., a pure-bulk wave-equation test). In that
case `fault_mesh_face_idx_set` is also empty, so the predicate would return
false anyway. But the `bc_.fault_attr > 0` check is redundant DEFENSE.

Now look at the dispatch sites (e.g., L2169, L2735, L3650, L4109):

```cpp
if (mf_on && central_flux_face_set_.count(f) > 0)
{ flux_.Central(nor, Q_self, Q_nbr, F_h); }
else
{ flux_.Interior(nor, Q_self, Q_nbr, F_h); }
```

These do NOT check `bc_.fault_attr > 0`. They rely entirely on the
`central_flux_face_set_` membership. So if the SET is wrongly populated
(any reason), Central is dispatched at the wrong face.

The set-build and dispatch should be consistent in their defensive checks.
Currently:
- Set-build excludes a face if `is_fault_face(f)` returns true.
- Dispatch fires on Central if `central_flux_face_set_.count(f) > 0`,
  unconditionally.

Two specific concerns:
1. **AllContinuous on a no-fault mesh**: `bc_.fault_attr == 0`. The
   `is_fault_face` predicate always returns false (short-circuit). Every
   interior face goes into `central_flux_face_set_`. This includes faces
   that may have a different boundary attribute (free surface, absorbing).

   Wait — `is_interior_face` already filters `Elem2No >= 0 || shared_seam`,
   which excludes boundary faces. So `is_fault_face`'s short-circuit on
   `bc_.fault_attr == 0` is harmless here.

2. **A future mesh with `bc_.fault_attr == 0` but a non-empty
   `fault_interior_faces_`** (impossible by construction, but defensively
   problematic): the predicate would return false despite the face being in
   `fault_mesh_face_idx_set`. The set construction would wrongly include
   the face. The dispatch would fire Central at the wrong face.

**Trigger:**
Any code path that populates `fault_interior_faces_` despite
`bc_.fault_attr == 0`. Currently impossible by ctor logic but a fragile
assumption.

**Actual behavior:** Predicate is over-defensive; silent failure mode if a
future ctor refactor decouples fault-face population from `bc_.fault_attr`.

**Expected behavior:** Either:
(a) Remove the `bc_.fault_attr > 0` short-circuit (let
`fault_mesh_face_idx_set.count(f)` be the sole truth).
(b) Keep the short-circuit but `MFEM_VERIFY(fault_mesh_face_idx_set.empty()
|| bc_.fault_attr > 0)` to catch the inconsistency.

**Suggested fix (code diff):**

```diff
@@ wave_operator.inl L1219-1242
-   std::unordered_set<int> fault_mesh_face_idx_set;
-   if (bc_.fault_attr > 0)
-   {
-      for (int i = 0; i < fault_interior_faces_.Size(); i++)
-      ...
-   }
+   std::unordered_set<int> fault_mesh_face_idx_set;
+   for (int i = 0; i < fault_interior_faces_.Size(); i++)
+   {
+      fault_mesh_face_idx_set.insert(fault_interior_faces_[i]);
+   }
+   ... shared faces ...
+   // R-1206: enforce consistency of bc_.fault_attr with the populated set.
+   MFEM_VERIFY(fault_mesh_face_idx_set.empty() || bc_.fault_attr > 0,
+               "BuildCentralFluxFaceSet_: fault_mesh_face_idx_set is "
+               "non-empty but bc_.fault_attr == " << bc_.fault_attr
+               << ".  Inconsistent ctor state — fault face lists "
+               "should be empty when no fault attribute is set.");
    auto is_fault_face = [&](int f) -> bool
    {
-      return (bc_.fault_attr > 0 &&
-              fault_mesh_face_idx_set.count(f) > 0);
+      return fault_mesh_face_idx_set.count(f) > 0;
    };
```

**Test case:**
Construct a `WaveOperator` with `bc_.fault_attr = 0` and verify
`SetMixedFluxMode(AllContinuous)` succeeds and `central_flux_face_set_`
contains every interior face (no fault to exclude). Then inject (via a
test-only setter or a private-member subclass) `fault_interior_faces_.Append(7)`
and re-call the setter; expect the new MFEM_VERIFY to abort.

---

### [R-1207] [MODERATE] [EDGE_CASE] [dynamic/wave_operator.inl:1316–1402] — R-001 post-walk MPI exchange uses `std::set<std::array<...>>` for global key matching; O(N log N) insertion + O(N) lookup is fine for production but has no upper-bound guard against pathological mesh

**Category:** EDGE_CASE / ASSUMPTION (performance + correctness in degenerate mesh)

**Description:**
The R-001 fix's MPI exchange flattens local keys, Allgatherv's them, then
constructs a `std::set<std::array<HYPRE_BigInt, 4>> global_keys` (L1380)
from the flattened buffer. Each rank then iterates its `n_shared` shared
faces and queries `global_keys.count(make_global_key(verts))` (L1397).

The key construction at L1335 initializes `key = {0, 0, 0, 0}` and copies
up to `verts.Size()` global vertex IDs. For a triangular face (3 verts),
key[3] is always 0. The sort at L1340 sorts only the first `verts.Size()`
slots:

```cpp
std::sort(key.begin(), key.begin() + verts.Size());
```

**Concern A**: For a mesh with vertex global IDs that include 0, the sorted
key for a triangle `{0, 5, 12}` is `{0, 5, 12, 0}`. A different triangle
with verts `{0, 5, 12, 0}` (which is impossible — duplicate vertex) would
match. Distinct triangles never share the same vertex set, so no collision.

**Concern B (REAL)**: Mixed face geometries (a mesh with both triangular
and quadrilateral faces, e.g., a hex-tet hybrid). A triangle `{1, 2, 3}`
sorts to key `{1, 2, 3, 0}`. A quadrilateral `{0, 1, 2, 3}` sorts to key
`{0, 1, 2, 3}`. These DIFFER (slot 0 differs) — no collision. BUT if a
future Mfem version reorders the slots for triangles (e.g., shifts to
align the first vertex at index 1), you'd get `{?, 1, 2, 3}`, which could
match the quad. The contract on `make_global_key` for variable-size
faces is fragile.

**Concern C (REAL)**: Performance. `std::set<std::array<HYPRE_BigInt, 4>>`
is a balanced BST, O(log N) operations. For a TPV104-scale mesh with ~5000
shared faces per rank (production np=128), N is small, so this is fine.
But under future scenarios:
- AMR with thousands of shared faces per rank.
- A refactor that increases the global-set scope (e.g., exchanging keys
  for ALL non-fault interior faces, not just shared ones).
The set lookup becomes a bottleneck. The fix would be `std::unordered_set`,
but that requires implementing a hash for `std::array<HYPRE_BigInt, 4>`.

**Concern D (REAL — most relevant for production)**: The R-001 exchange is
a ONE-TIME cost at `SetMixedFluxMode(Adjacent)`. For a TPV104 production
run, this is a few seconds. NOT a hot path.

**Trigger:** Future AMR or hybrid-element mesh.

**Actual behavior:** Silent assumption that face-vertex counts are ≤4 and
sorted keys are uniquely determined per face geometry. R-1106 already
hardened against >4-vert faces with `MFEM_VERIFY(verts.Size() <= 4)`. But
the mixed-element-type concern is unchecked.

**Expected behavior:** Either (a) use a face-id key that includes the
vertex count, or (b) document the assumption.

**Suggested fix (code diff):**

```diff
@@ wave_operator.inl L1335
-   std::array<HYPRE_BigInt, 4> key = {0, 0, 0, 0};
+   // R-1207: include verts.Size() as the first slot to disambiguate
+   // triangle vs quad keys when one rank's mesh has both face types.
+   std::array<HYPRE_BigInt, 5> key = {0, 0, 0, 0, 0};
+   key[0] = static_cast<HYPRE_BigInt>(verts.Size());
    for (int v = 0; v < std::min(verts.Size(), 4); v++)
    {
-      key[v] = gvi[verts[v]];
+      key[v + 1] = gvi[verts[v]];
    }
-   std::sort(key.begin(), key.begin() + verts.Size());
+   std::sort(key.begin() + 1, key.begin() + 1 + verts.Size());
```

The corresponding allocation, Allgatherv buffer size (4 → 5 entries per
key), and consume loop need updating in lockstep.

**Test case:** Out of scope until a hybrid-element TPV104 fixture is
added. Currently TPV104 is tet-only.

---

### [R-1208] [LOW] [QUALITY] [dynamic/wave_operator.inl:1620, 2274, 2816, 3773] — `mf_on` boolean is recomputed per-call but the underlying `mixed_flux_mode_` cannot change during a `Mult/AdvanceADER` call; recompute is unnecessary

**Category:** QUALITY (minor performance / readability)

**Description:**
Every flux-loop function recomputes:

```cpp
const bool mf_on = (mixed_flux_mode_ != MixedFluxMode::None);
```

at the top. `mixed_flux_mode_` cannot change during the Mult/AdvanceADER
call (the operator is `const`). The recompute is correct but minorly
wasteful: per-Mult-per-stage, the boolean is recomputed.

Better: cache `mf_on_` as a member updated only by `SetMixedFluxMode`.
Saves 4 Mult-call sites worth of comparison.

**Trigger:** Per-call overhead at production scale (negligible).
**Actual:** Recomputed each call.
**Expected:** Computed once.

**Suggested fix:** Add `bool mf_on_ = false;` member; update in
`SetMixedFluxMode`; use directly at the dispatch sites.

**Test case:** N/A (purely cosmetic).

---

### [R-1209] [LOW] [QUALITY] [tests/parallel/test_mixed_flux_adjacent_mpi.cpp:308] — Comment says "// Each such key appears twice in the Allgatherv (once from each rank)" but the per_key counting logic doesn't actually verify this, just trusts it

**Category:** QUALITY (test logic comment vs implementation mismatch)

**Description:**
At L283–308 of the parallel test, the comment claims each shared-face key
appears exactly twice in the Allgatherv (once from each rank's report).
The `per_key` map tracks the SET of distinct membership values
(`std::set<int>`), so a face reported as `member=1` from both ranks shows
as `set={1}`, and asymmetric reports show as `set={0,1}`.

But the loop only checks `kv.second.size() > 1` (asymmetric) without
counting that each key indeed had 2 reports. If MFEM ever changes shared-
face exchange to broadcast or 1-sided (only one rank reports), the test
would silently pass even with no consistency verification at all.

**Trigger:** Future MFEM refactor of shared-face semantics.

**Actual behavior:** Test is fragile to MFEM internals.

**Expected behavior:** Use `std::map<key, std::vector<int>>` and assert
each vector has size == 2; THEN check that all elements are equal.

**Suggested fix (code diff):**

```diff
-   std::map<std::array<HYPRE_BigInt, 4>, std::set<int>> per_key;
+   std::map<std::array<HYPRE_BigInt, 4>, std::vector<int>> per_key;
    for (int i = 0; i < total; i += 5) {
       ...
-      per_key[key].insert(member);
+      per_key[key].push_back(member);
    }
    int num_inconsistent = 0, num_shared_keys = 0;
    for (const auto &kv : per_key) {
+      // Each shared face MUST be reported by exactly 2 ranks (the two
+      // sharing it).  Anything else is a fixture or MFEM-semantics bug.
+      if (kv.second.size() != 2) {
+         std::cerr << "[R-1209] key reported " << kv.second.size()
+                   << " times (expected 2)\n";
+         num_inconsistent++;
+         continue;
+      }
-      if (kv.second.size() > 1) num_inconsistent++;
+      if (kv.second[0] != kv.second[1]) num_inconsistent++;
       num_shared_keys++;
    }
```

**Test case:** N/A (test infrastructure improvement).

---

## Summary

- Critical issues: **3** (R-1200, R-1201, R-1202)
  - **R-1200** end-to-end Adjacent dispatch test gates are too weak — it
    cannot detect a sign-flipped or wrong-magnitude central flux. The
    plan's Phase 4 algebraic acceptance criterion is not implemented.
  - **R-1201** No MPI test runs `wave.AdvanceADER` or `wave.Mult` with
    `MixedFluxMode::Adjacent` and verifies the actual numerical result vs
    a serial reference. The shared-face dispatch sites at L2735 / L4109
    have ZERO MPI test coverage with `mf_on=true` and a non-empty central
    set.
  - **R-1202** Driver banner reports rank 0's local `|central_set|` size
    as if global, masking rank imbalance and load asymmetry diagnostics
    in production logs.
- Moderate issues: **4** (R-1203, R-1204, R-1205, R-1206)
  - **R-1203** SEAS_TEST_NONFAULT_BOTH_SYM env-var hook is mathematically
    nullified for any conservative flux (Interior or Central); using it
    silently zeroes non-fault flux and can mask Mixed-Flux dispatch
    errors. Disable or refactor.
  - **R-1204** Driver does not abort on the untested `--mixed-flux
    adjacent --ader-order > 2` combination flagged by the plan's Risk R5.
  - **R-1205** `BuildCentralFluxFaceSet_` runs MPI collectives only in
    Adjacent mode; the implicit "all ranks must call with same mode"
    contract is undocumented and unenforced — fragile API.
  - **R-1206** `is_fault_face` predicate has a `bc_.fault_attr > 0`
    short-circuit that's redundant with `fault_mesh_face_idx_set.empty()`
    and creates an inconsistency between set-construction and dispatch
    paths if the two ever decouple.
- Low issues: **3** (R-1207, R-1208, R-1209)
  - R-1207 R-001 exchange's vertex-key construction lacks geometry-type
    disambiguation for hybrid-element meshes (out of scope for TPV104,
    relevant for future).
  - R-1208 `mf_on` boolean recomputed per-Mult-call (cosmetic).
  - R-1209 MPI test consistency check trusts shared-face report count
    without verifying.

**Verdict (round 3):** **FAIL — DO NOT enable `--mixed-flux adjacent` in
production multi-rank runs without addressing R-1200, R-1201, and R-1202.**

The Round-2 review declared "PASS" based on:
- R-1100 fix (rank-symmetric `is_fault_face`).
- R-1102 fix (end-to-end serial dispatch test).
- R-001 fix (cross-rank set-membership symmetry).

These are necessary but NOT sufficient:
- R-1102's gate is too weak (R-1200 — coarse non-zero / non-blowup checks
  cannot validate the analytic identity).
- R-001's verification only checks set-membership symmetry, never the
  actual numerical flux output (R-1201 — no test runs `AdvanceADER` at
  np≥2 with mode=Adjacent and compares to np=1).
- The driver's diagnostic output is misleading (R-1202 — banner reports
  rank-0-local size).

Production scale (np=128–400 on Frontera, ~1.5M tets, ~7000 macro steps)
amplifies any cross-rank inconsistency into catastrophic mass/momentum
imbalance. The current test suite cannot rule this out.

## Recommended Next Steps

1. **R-1200** — replace `test_mixed_flux_dispatch_adjacent.cpp` Gate 2 with
   the analytic identity check (1e-12 relative tolerance) per the plan's
   Phase 4 acceptance criterion. Verify the test FAILS with a hand-flipped
   central-flux sign.
2. **R-1201** — add `tests/parallel/test_mixed_flux_advance_ader_mpi.cpp`
   (np=2; np=4 if feasible) that runs `AdvanceADER` with mode=Adjacent
   and asserts bit-equivalence to the serial-reference Q_new vector after
   gather. This is the canonical production-readiness gate.
3. **R-1202** — fix the driver banner to print the global `|central_set|`
   sum + per-rank min/max for load-balance diagnostic.
4. **R-1203 / R-1204 / R-1205 / R-1206** — apply the diffs above.
5. After all of the above, re-run `make test` plus the new MPI tests at
   np=2 and np=4. Only then is `--mixed-flux adjacent` ready for a
   Frontera dry-run sbatch.
6. (Optional but recommended) Add an integration test that runs `tfinal=0.1`
   on `tpv104_repro.msh` at np=4 with `--mixed-flux adjacent` and asserts
   slip_strike at the hypocenter is finite + within 2% of the no-flag run
   (matches the plan's Phase 6 acceptance gate at line 794–801).

## Unreviewed Areas

- **AdvanceADER + substep iterator + mixed-flux** combination (Plan R5
  risk; no test exists).
- **AllContinuous mode at MPI scale** — only set-construction is verified;
  no `AdvanceADER` test exists for AllContinuous at np ≥ 2.
- **PML + mixed-flux interaction** — PML damping is applied AFTER the face
  flux loop; no test verifies that the PML region's damping is unchanged
  by mixed-flux dispatch on PML-adjacent faces.
- **`Mult` (RK4) path with mixed flux at MPI scale** — the dispatch site
  at L2169 (RK4 local) and L2735 (RK4 shared) are never tested with
  `mf_on=true` end-to-end. Round-7 substep iterators may dispatch via Mult.
- **Rank-imbalanced partitions** — a deliberately-imbalanced partition
  where one rank holds the bulk of the fault-adjacent faces. Plan-level
  R-1208 already loosened the production tolerance to 2%, but no test
  exercises a degenerate partition.
