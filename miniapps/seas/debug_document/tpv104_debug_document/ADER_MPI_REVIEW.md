# Code Review: ADER MPI/Parallel Support — TPV104

## Review Scope
- **Plan**: SCEC TPV104 spec + Pelties (2012) + `miniapps/seas/CLAUDE.md` (de-facto plan; no formal MPI plan doc)
- **Files reviewed (absolute paths)**:
  - `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.hpp` (lines 281–743)
  - `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/dynamic/wave_operator.inl` (predictor 977–1118; interior corrector 3070–4035; shared corrector 4040–4459; `AdvanceADER` 4460–4628)
  - `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/drivers/tpv104_driver.cpp` (time loop 1869–2018; substep wrapper 286–407)
  - `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/tests/unit/test_ader_*.cpp`, `test_fault_face_flux_ader_equivalence.cpp`, `test_rupture_multistep_serial_vs_parallel.cpp`
- **Domain context**: `CLAUDE.md` (sign conventions, FaultBasis t1=dip/t2=strike, shared-face Elem1-only assembly), R-101 / R-801 / R-305

ID range used: **R-1500 .. R-1510** (no collisions with R-1000..R-1009 driver, R-1100..R-1209 mixed-flux, R-1300..R-1306 substep iterator).

This review is distinct from `SUBSTEP_ITERATOR_MPI_REVIEW.md` (substep iterator) and `MIXED_FLUX_REVIEW_R*.md` (mixed flux). Focus is `wave.AdvanceADER` itself, the predictor, the corrector, and the ghost-data exchange between ADER stages.

---

## Pass 1 — De-facto plan compliance

| Plan requirement | Code site | Status |
|---|---|---|
| Predictor is element-local (CK recursion uses only `M_e^{-1} K_d^e`, no flux) | `ApplySpatialDerivative` `wave_operator.inl:822-924`; called from `ComputeADERTimeIntegrated:1013` | PASS — `MFEM_ASSERT(ndof == ndof_per_el_)` (line 860) plus the no-neighbour-read loop body confirm it. |
| Corrector exchanges I (time-integrated state) once across the rank seam | `ComputeADERSharedFaceFluxRHS:4099-4116` | PASS — single per-component `q_gf.ExchangeFaceNbrData()` loop with deep-copy guard. |
| Shared-fault uses canonical FaultBasis (BP5 convention) | `wave_operator.inl:4233-4253`; matches interior 3327-3344 | PASS — same `(qpd.normal, qpd.tangent1, qpd.tangent2, qpd.sign_flipped)` reconstruction. |
| Shared-face assembly only writes Elem1 RHS | `wave_operator.inl:4382-4391` | PASS — accumulates into `dof_offset1` only. |
| `bc_.fault_attr > 0` requires both `fault_flux_` and `fault_dof_data_` | shared `4078-4084`; interior `3122-3128` | PASS. |
| `MFEM_VERIFY(&Q != &Q_new)` aliasing guard | `AdvanceADER:4480` | PASS. |
| dt validation (>0) | `AdvanceADER:4474`; correctors echo it (3094, 4056) | PASS. |
| `pmesh.ExchangeFaceNbrData()` and `pfes->ExchangeFaceNbrData()` once at construction | `wave_operator.inl:137-140` | PASS — supported pattern; per-step `q_gf.ExchangeFaceNbrData()` exchanges values, not topology. |

**Plan-compliance gap**: parametric coverage at ADER-O ∈ {2, 3, 4} on the parallel path is missing. Test coverage is order-2 only (see R-1502).

---

## Findings

### [R-1500] (WITHDRAWN) Shared-fault corrector "one-sided flux" claim — re-analysis confirms correctness

**Status**: Withdrawn after re-analysis.

The side-selection at `wave_operator.inl:4376-4391` uses `I_imp_side = elem1_on_plus ? I_imp_plus_g : I_imp_minus_g` and `assemble_sign = elem1_on_plus ? -1.0 : +1.0`. Both ranks reconstruct bit-identical `can_n` from `FaultBasis`. Both call `EvaluateADER` deterministically on `(I_plus_local, I_minus_local)` and produce bit-identical `(I_imp_plus, I_imp_minus)`. Each rank then uses the half corresponding to its own Elem1. The redundant `EvaluateADER` call on each rank is **by design** to keep the friction-solver state writes (V1, slip_rate, theta, sigma_n_corr, tau_corr, Psi_diag) consistent across ranks per CLAUDE.md "shared-face semantics" / R-305. **Not a bug.**

---

### [R-1501] MODERATE [wave_operator.inl:AdvanceADER, ComputeADERTimeIntegrated] — Per-step heap allocation of NUM_STATE-sized buffers churns ~32 TB malloc traffic across a TPV104 production run

**Category:** RESOURCE / hot-loop allocation

**Description:**
`AdvanceADER` allocates `Vector I(NUM_STATE * ndof_total_)` and `Vector rhs(NUM_STATE * ndof_total_)` from scratch every macro-step. `ComputeADERTimeIntegrated` allocates THREE more (`D_curr` copy, `D_next`, `dQ_dxd`).

**Trigger:**
For TPV104 production at ~3 M tets, P=3 (20 DOFs/tet, 9 components):
- per-buffer ≈ 9 × 60 M = 540 MB
- 5 buffers per `AdvanceADER` call → 2.7 GB allocated+freed per macro-step
- 60 s run at dt=5e-3 s → 12 000 macro-steps → ~32 TB malloc/free traffic

This causes per-step page faults (RSS oscillation) and is borderline for jemalloc/glibc. The BP5/RK4 path already uses `mutable` member buffers (e.g. `ghost_gf_`) to avoid this.

**Actual behavior:** Re-allocate every step.
**Expected behavior:** Lazy-init `mutable` member buffers; reuse across macro-steps.

**Suggested fix:**
```diff
// wave_operator.hpp
+   mutable Vector ader_I_buf_;
+   mutable Vector ader_rhs_buf_;
+   mutable Vector ck_D_curr_buf_;
+   mutable Vector ck_D_next_buf_;
+   mutable Vector ck_dQ_dxd_buf_;
+   mutable std::vector<Vector> ader_nbr_data_;

// wave_operator.inl — AdvanceADER
-   Vector I(NUM_STATE * ndof_total_);
-   Vector rhs(NUM_STATE * ndof_total_);
+   const int N = NUM_STATE * ndof_total_;
+   if (ader_I_buf_.Size()   != N) { ader_I_buf_.SetSize(N); }
+   if (ader_rhs_buf_.Size() != N) { ader_rhs_buf_.SetSize(N); }
+   Vector &I = ader_I_buf_;
+   Vector &rhs = ader_rhs_buf_;
```
Mirror the same pattern in `ComputeADERTimeIntegrated` for `D_curr/D_next/dQ_dxd`. `mfem::Swap(D_curr, D_next)` swaps Vector handles in-place — safe with mutable members.

**Test case:**
```cpp
TEST_CASE("AdvanceADER scratch buffers reused across calls", "[ader][perf]") {
   ParMesh pmesh = BuildTwoTetFaultMesh();
   WaveOperator<ParMesh> wave(pmesh, /*p=*/2, lambda, mu, rho, bc);
   ...
   const auto rss0 = ResidentSetSizeKB();
   for (int step = 0; step < 100; ++step) {
      wave.AdvanceADER(Q, /*dt=*/1e-3, /*order=*/3, Q_new);
      Q.Swap(Q_new);
   }
   const auto rss1 = ResidentSetSizeKB();
   REQUIRE((rss1 - rss0) < 4096);  // 4 MB tolerance after 100 steps
}
```

---

### [R-1502] MODERATE [tests/unit/test_ader_interior_vs_shared_branch_live.cpp:191; test_fault_face_flux_ader_equivalence.cpp:318,401] — ADER parallel-equivalence tests pin order=2; production runs at O=3,4 are uncovered

**Category:** TEST COVERAGE

**Description:**
The only parallel-vs-serial ADER equivalence test (`test_ader_interior_vs_shared_branch_live.cpp`) hard-codes `order=2` at line 191:
```cpp
wave.AdvanceADER(Q, kDt, /*order=*/2, Q_new);
```
At ADER-O=2 the CK recursion exits after one pass (line 1007: `for (int k = 0; k < order - 1; ++k)`), masking any bug that fires at higher orders. Production sbatch scripts (`tpv104_hybrid_200m_O3.sbatch`, `_O4.sbatch`) run ADER-O=3 and 4, where:
- the recursion takes additional iterations
- the shared-face exchange of `I` carries higher-derivative information
- a freshness or coefficient bug would emerge but currently has zero parallel coverage

`test_fault_face_flux_ader_equivalence.cpp` is serial-only and also pins order=2.

**Trigger:** Production sbatch with `--ader-order 3` or `--ader-order 4` at np≥2.

**Actual behavior:** Tests pass at O=2; O=3/4 untested at np>1.
**Expected behavior:** Parametric sweep over O ∈ {2, 3, 4} on the parallel path.

**Suggested fix:** convert `main` into a templated helper and sweep:
```diff
+template <int ORDER>
+int RunSliceCheck(int rank, int nprocs) {
+   ... existing body, replacing 2 with ORDER ...
+   wave.AdvanceADER(Q, kDt, /*order=*/ORDER, Q_new);
+}
+
+int main(int argc, char *argv[]) {
+   MPI_Init(&argc, &argv); int rank=0, nprocs=1;
+   MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
+   if (nprocs != 2) { MPI_Finalize(); return 77; }
+   int rc = 0;
+   rc |= RunSliceCheck<2>(rank, nprocs);
+   rc |= RunSliceCheck<3>(rank, nprocs);
+   rc |= RunSliceCheck<4>(rank, nprocs);
+   MPI_Bcast(&rc, 1, MPI_INT, 0, MPI_COMM_WORLD);
+   MPI_Finalize();
+   return rc;
+}
```

**Test case:**
```cpp
TEST(ADER_MPI, InteriorVsSharedBranchEquivalenceAtAllOrders) {
   for (int O : {2, 3, 4}) {
      Vector slope_serial(NUM_STATE), slope_parallel(NUM_STATE);
      // serial reference run with wave.AdvanceADER(..., O, ...)
      // parallel run with wave.AdvanceADER(..., O, ...)
      for (int c = 0; c < NUM_STATE; ++c) {
         const auto rel = std::abs(slope_serial(c) - slope_parallel(c)) /
                          std::max({std::abs(slope_serial(c)),
                                    std::abs(slope_parallel(c)), 1.0});
         EXPECT_LT(rel, 1e-10) << "ADER-O=" << O << " comp=" << c;
      }
   }
}
```

---

### [R-1503] MODERATE [wave_operator.inl:1019 vs 1108] — CK recursion factorial denominator differs between TimeIntegrated `(k+2)` and SubStep `(k+1)`; both correct, but maintainability hazard

**Category:** QUALITY (maintainability)

**Description:**
Two CK recursions parameterise the factorial denominator differently:
- `ComputeADERTimeIntegrated:1019`: `fac *= dt / (k + 2)` (because `fac` was initialised to `dt` BEFORE the loop and the k=0 contribution was pre-added at line 1005, so at iteration k it represents `dt^{k+1}/(k+1)!`).
- `ComputeADERSubStepStates:1108`: `denom = k + 1` (because `fac[o]` was initialised to `1.0` and the k=0 contribution was pre-added outside the loop at line 1086-1089).

Both are mathematically correct, but a future maintainer will see `(k+2)` in one and `(k+1)` in the other and "fix" one to match.

**Trigger:** Future refactor or "correction" attempt.

**Suggested fix:** add cross-reference comments documenting why each recursion is correct.
```diff
// wave_operator.inl:1018 (TimeIntegrated)
+// fac was initialised to `dt` at line 1005 (k=0 pre-added), so at iter k
+// it represents dt^{k+1}/(k+1)! and updates to dt^{k+2}/(k+2)! for D(k+1).
+// ComputeADERSubStepStates uses `denom = k+1` because its k=0 pre-add
+// (line 1088) is at fac=1.0, so its `k` is one step behind ours.
 fac *= dt / static_cast<real_t>(k + 2);
```

**Test case:**
```cpp
TEST_CASE("CK: SubStep ∫ Q(τ)dτ matches TimeIntegrated I", "[ader][predictor]") {
   const real_t dt = 1e-3; const int O = 4;
   auto nodes = GaussLegendreNodes(O, 0.0, dt);
   auto weights = GaussLegendreWeights(O, 0.0, dt);

   Vector Q = MakeRandomState(...);
   Vector I_TI; wave.ComputeADERTimeIntegrated(Q, dt, O, I_TI);

   std::vector<Vector> Q_per_node;
   wave.ComputeADERSubStepStates(Q, dt, O, nodes, Q_per_node);

   Vector I_SS(I_TI.Size()); I_SS = 0.0;
   for (int o = 0; o < (int)nodes.size(); ++o)
      { I_SS.Add(weights[o], Q_per_node[o]); }

   for (int i = 0; i < I_TI.Size(); ++i) {
      const real_t rel = std::abs(I_TI(i) - I_SS(i)) /
                         std::max(std::abs(I_TI(i)), 1.0);
      REQUIRE(rel < 1e-12);
   }
}
```

---

### [R-1504] MODERATE [wave_operator.inl:4105-4116] — `q_gf.ExchangeFaceNbrData()` invoked NUM_STATE=9× per `AdvanceADER` — sequential point-to-point latency dominates parallel scaling

**Category:** PERFORMANCE / MPI

**Description:**
The shared-fault corrector (and siblings `ComputeSharedFaceFluxRHS:2624-2654`, `EvaluateBulkAtFaultQPsCanonical:1701-1712`) loops `c = 0..NUM_STATE-1` and calls `q_gf.ExchangeFaceNbrData()` per component. NUM_STATE=9. Each call is a separate MPI send/receive sequence over the rank topology graph (~6-26 neighbours on a 3D partition).

For a TPV104 production run: 12 000 macro-steps × 9 components × ~3 exchange sites per step ≈ 324 k MPI rounds. At 100 µs/round on Frontera that's 32 s of pure latency.

**Trigger:** Any np≥2 ADER run (always).

**Actual behavior:** 9 sequential per-component exchanges per shared corrector call.
**Expected behavior:** 1 vector-valued exchange per call using a `vdim=NUM_STATE` `ParGridFunction`.

**Suggested fix** (sketch — non-trivial refactor):
```diff
// wave_operator.hpp
+   mutable std::unique_ptr<ParFiniteElementSpace> pfes_full_;
+   mutable std::unique_ptr<ParGridFunction> ghost_gf_full_;

// wave_operator.inl — ctor (parallel branch)
+   pfes_full_ = std::make_unique<ParFiniteElementSpace>(
+      pfes->GetParMesh(), pfes->FEColl(), NUM_STATE,
+      mfem::Ordering::byNODES);
+   ghost_gf_full_ = std::make_unique<ParGridFunction>(pfes_full_.get());

// Replace the per-component loop:
-   for (int c = 0; c < NUM_STATE; ++c) {
-      ParGridFunction q_gf(pfes); q_gf.SetFromTrueDofs(I_c);
-      q_gf.ExchangeFaceNbrData();
-      nbr_data[c] = q_gf.FaceNbrData();
-   }
+   ghost_gf_full_->SetFromTrueDofs(I);   // I is NUM_STATE * ndof_total_
+   ghost_gf_full_->ExchangeFaceNbrData();
+   const Vector &nbr_full = ghost_gf_full_->FaceNbrData();
```
Caveat: shared with `ComputeSharedFaceFluxRHS` which has a deep-copy correctness guard at line 2636-2639 — the refactor must preserve that invariant. Treat as a separate ticket; provide microbenchmark first.

**Test case:**
```cpp
TEST_CASE("AdvanceADER shared-corrector exchanges scale O(1) in NUM_STATE",
          "[ader][mpi][.perf]") {
   ParMesh pmesh = ...;
   WaveOperator<ParMesh> wave(pmesh, ...);
   ResetMPISendCounter();
   for (int s = 0; s < 100; ++s) wave.AdvanceADER(Q, dt, /*O=*/3, Q_new);
   const int sends = MPISendCounter();
   const int neighbors = pmesh.GetNFaceNeighbors();
   // Pre-fix: ≤ 9 * 100 * neighbors. Post-fix: ≤ 1 * 100 * neighbors.
   REQUIRE(sends <= 100 * neighbors * 9);  // pre-fix bound
}
```

---

### [R-1505] LOW [wave_operator.inl:4073, 3105] — `MFEM_VERIFY(has_bulk_bg_, ...)` is rank-local; rank-asymmetric calls deadlock at the next collective

**Category:** DEFENSIVE PROGRAMMING

**Description:**
If a future driver path forgets to call `SetAbsorbingBackground(Q_bg)` on a subset of ranks, `MFEM_VERIFY` fires rank-locally → that rank aborts; other ranks proceed and DEADLOCK at `q_gf.ExchangeFaceNbrData()` (line 4111).

**Suggested fix:** mirror the `SetMixedFluxMode` collective-consensus pattern (line 1163-1183):
```diff
+#ifdef MFEM_USE_MPI
+   if constexpr (IsParallelMesh<MeshType>::value) {
+      auto &pmesh = static_cast<ParMesh &>(mesh_);
+      const int my_has = has_bulk_bg_ ? 1 : 0;
+      int min_has = 0;
+      MPI_Allreduce(&my_has, &min_has, 1, MPI_INT, MPI_MIN, pmesh.GetComm());
+      MFEM_VERIFY(min_has == 1, "AdvanceADER requires SetAbsorbingBackground "
+                  "on EVERY rank — rank-local call would deadlock at the "
+                  "next ExchangeFaceNbrData.");
+   } else
+#endif
    { MFEM_VERIFY(has_bulk_bg_, "..."); }
```

---

### [R-1506] LOW [wave_operator.inl:4471-4628; wave_operator.hpp:399-405] — `Q_new.SetSize(...)` hoisted to the bottom of `AdvanceADER` instead of the top

**Category:** API CONTRACT

**Description:**
`Q_new.SetSize(NUM_STATE * ndof_total_)` happens at line 4626 — AFTER all the corrector work. The PML branch at 4499-4563 modifies `rhs`, not `Q_new`, so this is benign today. But a future path that pre-stages into `Q_new` would silently misbehave.

**Suggested fix:**
```diff
    MFEM_VERIFY(&Q != &Q_new, "AdvanceADER: Q and Q_new must be distinct Vectors");
+   Q_new.SetSize(NUM_STATE * ndof_total_);  // hoisted from line 4626
```

---

### [R-1507] LOW [wave_operator.inl:822-924] (DECLINED OPTIMISATION) — `ApplySpatialDerivative` recomputes shape/dshape every CK iteration; precompute would cost 28 GB on production

Per-step recompute is correct and acceptable for production. Caching `K_d^e M_e^{-1}` per element per direction would cost `3 * ne * ndof_per_el² * sizeof(real_t)` ≈ 28 GB at 3 M tets P=3. Declined.

---

### [R-1508] LOW [wave_operator.inl:4200-4201, 860-863] — Mixed-element ghost-layer guard is `MFEM_ASSERT` (Debug-only); promote to `MFEM_VERIFY` for Release safety

**Category:** MESH-TYPE ASSUMPTION

**Suggested fix:**
```diff
- MFEM_ASSERT(ndof2 == ndof_per_el_, "Mixed element types in ghost");
+ MFEM_VERIFY(ndof2 == ndof_per_el_,
+             "ComputeADERSharedFaceFluxRHS: heterogeneous element types in "
+             "the ghost layer are not supported (ghost ndof=" << ndof2
+             << " vs ndof_per_el_=" << ndof_per_el_ << ").");
```
Mirror at `ApplySpatialDerivative:860`.

---

### [R-1509] LOW [wave_operator.inl:4474-4481] — Redundant size validation; declined

`AdvanceADER` validates `Q.Size()` and allocates `I/rhs` itself; downstream `ComputeADERVolumeUpdate / FaceFluxRHS` also `MFEM_VERIFY`. Redundant but harmless. No fix needed.

---

### [R-1510] LOW [tpv104_driver.cpp:1869-1928; wave_operator.inl:137-140, 4471] — Static-mesh assumption is undocumented; future AMR/moving-mesh extensions would silently consume stale topology

**Category:** DOCUMENTATION

**Suggested fix:**
```diff
+/// @warning Assumes a static ParMesh.  ParMesh::ExchangeFaceNbrData() and
+/// ParFiniteElementSpace::ExchangeFaceNbrData() are called ONCE in the
+/// ctor (line 137-138) and topology is not re-exchanged across macro-steps.
+/// A future AMR or moving-mesh extension MUST call them before every
+/// AdvanceADER if the partition or face-neighbor graph mutates.
 template <typename MeshType>
 void WaveOperator<MeshType>::AdvanceADER(...)
```

---

## Pass 3 — Quality strengths (notable absent issues)

1. **Consistent canonical-frame reconstruction** between interior and shared corrector branches (`3327-3344` vs `4239-4253`). The `qpd.sign_flipped` plus per-face `elem1_on_plus` correctly implements R-101 / R-801 across the partition seam.
2. **Element-locality of the predictor** enforced by `MFEM_ASSERT(ndof == ndof_per_el_)` at line 860 and absence of any `pfes->GetFaceNbrFE()` call inside `ApplySpatialDerivative` / `ComputeADERTimeIntegrated`.
3. **Deep-copy guard for ghost data** (line 4112-4115, mirror at 2636-2639): `MFEM_VERIFY(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(), "H2 regression")`.
4. **`SetMixedFluxMode` collective discipline** (line 1163-1183): MPI_Allreduce(MIN/MAX) catches mode-mismatched callers BEFORE the downstream `MPI_Allgatherv` deadlocks.
5. **R-1003 substep dispatch** correctly gates on absolute index range `dof_idx < substep_n_total_fault_qps_` on both interior and shared branches (line 3399-3402, 4286-4289).
6. **Shared-fault SBE assembly** correctly writes only Elem1 RHS (line 4382-4391); peer rank handles its own Elem1.
7. **No predictor stale-ghost bug.** Predictor is element-local; no neighbour data sampled.
8. **No double-application of absorbing BC at partition-cut boundary faces.** Dispatch at line 3193 correctly skips shared faces in the local boundary loop.
9. **No subcycling stage-exchange race.** `AdvanceADER` is single-stage; corrector exchanges `I` AFTER predictor CK recursion completes.
10. **CK recursion parameterised correctly for O ∈ {2, 3, 4}.** Test coverage gap (R-1502) hides regression risk but the code itself is correct at O ≤ 4.

---

## Findings summary

| ID | Severity | Category | File:Line |
|---|---|---|---|
| R-1500 | (withdrawn) | — | — |
| R-1501 | MODERATE | Resource management | `wave_operator.inl:4484, 4488, 999-1001` |
| R-1502 | MODERATE | Test coverage | `test_ader_interior_vs_shared_branch_live.cpp:191`; `test_fault_face_flux_ader_equivalence.cpp:318,401` |
| R-1503 | MODERATE | Maintainability | `wave_operator.inl:1019 vs 1108` |
| R-1504 | MODERATE | Performance/MPI | `wave_operator.inl:4105-4116` |
| R-1505 | LOW | Defensive programming | `wave_operator.inl:4073, 3105` |
| R-1506 | LOW | API contract | `wave_operator.inl:4471-4628`; header 399-405 |
| R-1507 | LOW | Performance (declined) | `wave_operator.inl:822-924` |
| R-1508 | LOW | Mesh-type assumption | `wave_operator.inl:4200-4201, 860-863` |
| R-1509 | LOW | API contract (declined) | `wave_operator.inl:4474-4481` |
| R-1510 | LOW | Documentation | `tpv104_driver.cpp:1869-1928`; `wave_operator.inl:137-140, 4471` |

- Critical issues: 0
- Moderate issues: 4
- Low issues: 6
- Plan compliance: FULL (test-coverage gap noted in R-1502)

**Verdict: PASS WITH FIXES** — ADER MPI/parallel support in `wave.AdvanceADER` is **production-ready for TPV104 at np>1**. No CRITICAL bugs found.

Two MODERATE physics-relevant items (R-1502, R-1503) and two MODERATE performance items (R-1501, R-1504) should be addressed before declaring full production maturity. The originally-suspected "shared-fault one-sided flux" finding (R-1500) was withdrawn after re-analysis: the apparent redundancy of running `EvaluateADER` on both ranks of a shared face is by design (deterministic friction-solver state writes mirrored on both sides per CLAUDE.md "shared-face semantics" / R-305).

## Unreviewed Areas

- **`AdvanceADERWithSubStep`** (substep iterator path) — covered by `SUBSTEP_ITERATOR_MPI_REVIEW.md`.
- **Mixed-flux ADER dispatch** — covered by `MIXED_FLUX_REVIEW_R*.md`.
- **PML branch** (`wave_operator.inl:4499-4563`) — out of scope for TPV104; pure-elastic half-space production path bypasses PML.
- **AMR / moving-mesh ADER** — not implemented; flagged via R-1510 as future risk.
