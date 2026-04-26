# Code Review: Mixed-Flux DG implementation (Round-11) — 2026-04-25

## Review Scope
- Plan: `miniapps/seas/debug_document/tpv104_debug_document/MIXED_FLUX_PLAN.md`
- Prior plan-level audit: `miniapps/seas/debug_document/tpv104_debug_document/MIXED_FLUX_PLAN_REVIEW.md` (R-1201–R-1209)
- Files reviewed:
  - `miniapps/seas/dynamic/godunov_flux.hpp` (Central decl)
  - `miniapps/seas/dynamic/godunov_flux.cpp` (Central impl, lines 382–423)
  - `miniapps/seas/dynamic/wave_operator.hpp` (MixedFluxMode enum + setter/accessors + private members)
  - `miniapps/seas/dynamic/wave_operator.inl` (SetMixedFluxMode, BuildCentralFluxFaceSet_, 4 dispatch sites at L2034 / L2604 / L3515 / L3973)
  - `miniapps/seas/drivers/tpv104_driver.cpp` (CLI + banner + setter call site)
  - `miniapps/seas/tests/unit/test_godunov_central_flux.cpp` (R-1101)
  - `miniapps/seas/tests/unit/test_mixed_flux_face_set.cpp` (R-1102)
  - `miniapps/seas/tests/unit/test_mixed_flux_dispatch_none.cpp` (R-1103)
  - `miniapps/seas/tests/unit/test_tpv104_smoke.cpp` (banner gate extension)
- Domain context consulted: `CLAUDE.md`, `miniapps/seas/CLAUDE.md` (Critical Numerical Details: DG method, MPI shared-face handling, fault-local frame conventions).

## Findings

### [R-001] [CRITICAL] [wave_operator.inl:BuildCentralFluxFaceSet_] — Adjacent mode is not MPI-consistent across rank seams

**Category:** BUG

**Description:**
`BuildCentralFluxFaceSet_` builds `central_flux_face_set_` from purely *local* information on each rank: `fault_interior_faces_` (local), `fault_shared_faces_` (local fault-side seams), and `mesh_.GetElementFaces` (local elements). For `MixedFluxMode::Adjacent` it constructs `E_fault_adj` from local fault-adjacent elements only (lines 1230–1255), then walks each such element's faces and inserts the *interior + non-fault* ones (lines 1261–1272). There is no MPI exchange.

**Trigger:**
Multi-rank ParMesh run with `--mixed-flux adjacent` whenever the partition cuts between a fault-adjacent element on one rank and a non-fault-adjacent element on another rank. Concretely: a *non-fault* shared face whose local element on rank A is fault-adjacent (touches a fault face on A), but whose remote element on rank B is **not** fault-adjacent on B (no fault face on B's local view of that element).

**Actual behavior:**
- Rank A: walks `E_fault_adj` ∋ X → walks X's faces → hits the shared face f → `is_interior_face(f) = true` (via `shared_mesh_face_set_.count(f) > 0`), `is_fault_face(f) = false` → `central_flux_face_set_.insert(f)`. At dispatch, rank A calls `flux_.Central(nor_A, Q_X, Q_W_ghost, F_h)`.
- Rank B: `E_fault_adj` does **not** include rank B's local element W (W has no fault face locally; only A's element is fault-adjacent, but A's fault-adjacency is not communicated to B). Consequently f ∉ rank-B's `central_flux_face_set_`. At dispatch, rank B calls `flux_.Interior(nor_B, Q_W, Q_X_ghost, F_h)`.

**Expected behavior:**
Both ranks must agree on Central vs Interior dispatch for every shared face, otherwise the inter-rank flux is non-conservative.

**Why it's a bug (algebra):**
Let `F_pos := flux_.Central(nor, Q_self, Q_nbr)` and recall the algebraic identity verified in `test_godunov_central_flux` Gate 2:
```
F_interior(nor, Q_self, Q_nbr) = F_central(nor, Q_self, Q_nbr) + 0.5·|A_n|·(Q_nbr − Q_self)
```
For a shared face with rank-A normal `nor_A` and rank-B normal `nor_B = −nor_A`:
- Both-Central: `F_A_central(nor_A, Q_X, Q_W) + F_B_central(nor_B, Q_W, Q_X) = 0` (anti-symmetric, conservative).
- Both-Interior: `F_A_int + F_B_int = 0` (anti-symmetric, conservative).
- Mixed (A=Central, B=Interior): `F_A_central + F_B_interior = 0.5·|A_n|·(Q_X − Q_W) ≠ 0`.

This produces a spurious volume source/sink of magnitude `0.5·|A_n|·||jump||` per ill-dispatched shared face per macro-step. For TPV104 production with `|A_n| ~ cp = 6×10³ m/s` and stress jumps of order `10⁷ Pa` between fault-adjacent and bulk-far elements, the per-face-per-step error is `~3×10¹⁰ Pa·m/s` of unbalanced flux. Over ~7000 macro-steps, the accumulated mass/momentum/stress imbalance is catastrophic.

**Why R-1207 doesn't fix this:**
R-1207's `is_interior_face` predicate only ensures shared seams are *recognized* as interior faces in the iteration. It does not ensure that the predicate's verdict on a face is the same on both ranks — and it can't, because the membership decision is made by walking E_fault_adj, which is a purely local set.

The seas_test_R002_empty_fault_rank_no_abort test exists for exactly this topology (one rank has no fault face) but only checks that the empty-fault rank doesn't abort — not that it dispatches consistently with its fault-adjacent neighbor.

**Suggested fix:**
After populating `E_fault_adj` and walking local elements (lines 1230–1272), exchange fault-adjacency information across shared face seams using the same MPI primitive `ParMesh` already exposes. Two-step plan:

1. Build a per-shared-face boolean `local_has_fault_adj_elem[sf]` indicating whether the local Elem1 on shared face `sf` is in `E_fault_adj`.
2. Exchange via `pmesh.GroupCommunicator` or pairwise MPI_Sendrecv along shared-face neighbor lists. After exchange, each rank also knows `remote_has_fault_adj_elem[sf]`.
3. For every shared face `sf` whose `local_has_fault_adj_elem[sf] || remote_has_fault_adj_elem[sf]` is true, insert `pmesh.GetSharedFace(sf)` into `central_flux_face_set_` (if it isn't already).

```diff
   // Step 2: for each element in E_fault_adj, walk its faces; insert any
   // interior non-fault face into central_flux_face_set_.
   Array<int> elem_faces, elem_orient;
   for (int e : E_fault_adj)
   {
      mesh_.GetElementFaces(e, elem_faces, elem_orient);
      for (int j = 0; j < elem_faces.Size(); j++)
      {
         const int f = elem_faces[j];
         if (is_interior_face(f) && !is_fault_face(f))
         {
            central_flux_face_set_.insert(f);
         }
      }
   }
+
+#ifdef MFEM_USE_MPI
+   // R-001 fix: a shared non-fault face whose REMOTE element is fault-
+   // adjacent (but local element isn't) must also be inserted, otherwise
+   // the two ranks dispatch different fluxes on the same physical face
+   // and conservation breaks across the rank seam.
+   if constexpr (IsParallelMesh<MeshType>::value)
+   {
+      auto &pmesh = static_cast<ParMesh &>(mesh_);
+      const int n_shared = pmesh.GetNSharedFaces();
+      std::vector<int> local_adj(n_shared, 0);
+      for (int sf = 0; sf < n_shared; sf++)
+      {
+         FaceElementTransformations *ftr =
+            pmesh.GetSharedFaceTransformations(sf);
+         if (ftr && ftr->Elem1No >= 0 &&
+             E_fault_adj.count(ftr->Elem1No) > 0)
+         {
+            local_adj[sf] = 1;
+         }
+      }
+      std::vector<int> remote_adj(n_shared, 0);
+      // Exchange local_adj → remote_adj across the shared-face graph.
+      // (Implementation: pmesh.GroupComm() or per-neighbor Sendrecv on
+      //  the shared-face neighbor list; details follow standard ParMesh
+      //  shared-face exchange used elsewhere in this file.)
+      ExchangeSharedFaceFlags_(pmesh, local_adj, remote_adj);
+      for (int sf = 0; sf < n_shared; sf++)
+      {
+         if (local_adj[sf] || remote_adj[sf])
+         {
+            const int f = pmesh.GetSharedFace(sf);
+            if (!is_fault_face(f))   // shared fault face still excluded
+            {
+               central_flux_face_set_.insert(f);
+            }
+         }
+      }
+   }
+#endif
}
```

The `ExchangeSharedFaceFlags_` helper is short (≤30 LOC) — pack `local_adj` per shared face, MPI_Sendrecv with each neighbor rank, unpack into `remote_adj`. Reuse the shared-face neighbor list already maintained by `ParMesh`.

**Test case** (multi-rank, currently absent):
```cpp
// tests/parallel/test_mixed_flux_adjacent_mpi_consistency.cpp
// Run with `mpirun -np 2 ./seas_test_mixed_flux_adjacent_mpi_consistency`.
//
// Build a 4-tet ParMesh with the fault on a face that is INTERIOR to
// rank 0 (so it lands in fault_interior_faces_ on rank 0 and rank 1
// has zero local fault faces).  Construct it so rank 0's fault-adjacent
// element shares a non-fault face with one of rank 1's elements.
//
// After SetMixedFluxMode(Adjacent):
//   - Gather central_flux_face_set_ membership for the shared face from
//     both ranks (Allreduce of a per-shared-face mask).
//   - Assert rank-0 membership == rank-1 membership for the shared face.
void test_R001_adjacent_consistency_across_seam() {
    // ... setup 4-tet ParMesh ...
    wave.SetMixedFluxMode(MixedFluxMode::Adjacent);
    int local_in = (rank == 0)
        ? wave.GetCentralFluxFaceSet().count(local_shared_face_idx)
        : wave.GetCentralFluxFaceSet().count(local_shared_face_idx);
    int peer_in = 0;
    MPI_Sendrecv(&local_in, 1, MPI_INT, peer_rank, 0,
                 &peer_in,  1, MPI_INT, peer_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    TEST_ASSERT(local_in == peer_in,
                "Adjacent: shared non-fault face has same membership on "
                "both sides of the rank seam");

    // Stronger test: run AdvanceADER one step in Adjacent mode and check
    // total stress integral is conserved up to FP precision (sum over
    // all rank-local SXX·dV minus the same sum at t=0 — should be 0
    // modulo boundary fluxes).  With the bug, the residual is
    // O(0.5·|A_n|·jump · |shared-face-set asymmetric|).
}
```

---

### [R-002] [MODERATE] [tests/unit/test_mixed_flux_*.cpp] — No multi-rank coverage of mixed-flux dispatch

**Category:** EDGE_CASE / coverage gap

**Description:**
The three mixed-flux unit tests (`test_godunov_central_flux`, `test_mixed_flux_face_set`, `test_mixed_flux_dispatch_none`) all use serial `Mesh` fixtures (`BuildSmallFaultTetMesh`, `BuildTwoTetFaultMesh`). None instantiates a ParMesh and none calls `mpirun -np N`. The dispatch flow at `ComputeSharedFaceFluxRHS:2604` and `ComputeADERSharedFaceFluxRHS:3973` is not exercised at all — neither for correctness on a single MPI rank nor for cross-rank consistency on multiple ranks.

**Trigger:**
Any production multi-rank run with `--mixed-flux adjacent`. The shared-face dispatch executes, the `central_flux_face_set_.count(mesh_face_idx)` lookup happens, but no test verifies its outcome.

**Actual behavior:**
The implementation passes all unit tests (3+9+2 = 14), but the test suite is silent on:
- Whether the shared-face dispatch site at L2604 / L3973 even *fires* with the right semantics (it may trivially fire because `mf_on=false` short-circuits in tests, masking dispatch logic bugs).
- Whether `BuildCentralFluxFaceSet_` correctly populates the set with shared mesh face indices (the predicate has a `shared_seam` branch, but the test fixture never has shared faces).
- Whether R-001's cross-rank inconsistency manifests.

**Expected behavior:**
At minimum, one parallel test exercising:
1. ParMesh construction with at least one shared *non-fault* face whose endpoints can both be fault-adjacent (or, for the R-001 bug, asymmetric in fault-adjacency).
2. `wave.SetMixedFluxMode(MixedFluxMode::Adjacent)`.
3. Cross-rank consistency assertion on the shared face's membership in `central_flux_face_set_`.
4. AdvanceADER bit-equivalence between `--mixed-flux none` (default) and an explicit `SetMixedFluxMode(MixedFluxMode::None)` on the same multi-rank fixture (R-1103 analog for ParMesh).

**Suggested fix:**
Create `tests/parallel/test_mixed_flux_mpi_dispatch.cpp` with two sub-tests:

```cpp
// Sub-test 1: dispatch-none bit-identity on ParMesh (R-1103 analog).
void test_R002_dispatch_none_parmesh()
{
    // Build a 12-tet ParMesh fault fixture; partition with METIS.
    // Path A: never call SetMixedFluxMode.
    // Path B: explicitly call SetMixedFluxMode(None).
    // Run AdvanceADER one step; assert max |Q_new_a - Q_new_b| == 0
    // bit-equal on every rank, then MPI_Allreduce to ensure global zero.
}

// Sub-test 2: shared-face set membership symmetry (R-001 verification).
void test_R002_shared_face_membership_symmetric()
{
    // Build 4-tet ParMesh with fault crossing rank seam such that one
    // rank has zero local fault faces but its element is the non-fault
    // partner of a fault-adjacent element on the other rank.
    // After SetMixedFluxMode(Adjacent), exchange membership flags for
    // the shared seam face; assert local == remote on every shared face.
}
```

These need a 4–12 tet hand-built fault fixture compatible with METIS partitioning across 2 MPI ranks. Add `seas_test_mixed_flux_mpi_dispatch` Makefile target following the existing `seas_test_R002_empty_fault_rank_no_abort` template, and a `mpirun -np 2` invocation in the test runner.

**Test case:**
See sub-tests 1 and 2 above — they ARE the test cases.

---

### [R-003] [LOW] [test_godunov_central_flux.cpp:166,235,282] — Test tolerances exceed plan-acceptance by 9 orders of magnitude

**Category:** QUALITY / DEVIATION

**Description:**
`MIXED_FLUX_PLAN.md` Phase 1 acceptance criterion specifies `1e-12 absolute` tolerance for the algebraic identity (Gate 2). The implemented test uses `1e-1 absolute` for Gate 2 (line 235) and Gate 3 (line 282), and `1e-3` for Gate 1 (line 166). The comments at lines 144–146, 228–234, 280–282 explain that the output scale is `~6e10` so 1 ULP in double precision is `~1e-5` absolute, making `1e-12` infeasible — this is a planning bug, not an implementation bug. The implementation is correct.

**Trigger:**
Re-reading the plan-vs-test tolerance comparison.

**Actual behavior:**
Gate 2 passes at `1.56e-2 absolute` (`~1e-12 relative`), Gate 1 passes at `1.36e-6 absolute`, Gate 3 passes at `5.25e-2 absolute`. Tolerances are 1e-1 / 1e-3 / 1e-1.

**Expected behavior:**
Either:
(a) Update `MIXED_FLUX_PLAN.md` Phase 1 acceptance to read "1e-12 *relative* (≈ 1e-1 absolute at typical TPV104 stress scales)" so plan and test agree;
(b) Change tests to express tolerance relatively: `worst / max(|F_up|, ...) < 1e-12`. This would auto-scale and document the contract more precisely.

**Suggested fix:**
Option (b) is more durable — it survives changes to the test material parameters or state magnitudes:

```diff
- TEST_LE(worst_g2, 1.0e-1,
-         "R-1204 algebraic identity: F_up − F_ce = 0.5·|A_n|·(Q_self−Q_nbr) "
-         "to FP precision (~1 ULP × scale)");
+ // Express tolerance relative to the per-pair output magnitude so the
+ // plan's "1e-12 relative" contract is checked directly.
+ real_t worst_rel_g2 = 0.0;
+ for (int p = 0; p < 3; p++)
+    for (int k = 0; k < 8; k++)
+    {
+       real_t F_up[NUM_STATE], F_ce[NUM_STATE];
+       flux.Interior(normals[k], pairs[p].Q_self, pairs[p].Q_nbr, F_up);
+       flux.Central (normals[k], pairs[p].Q_self, pairs[p].Q_nbr, F_ce);
+       real_t abs_an_jump_global[NUM_STATE];
+       ComputeAbsAnTimesJumpGlobal(flux, normals[k],
+                                   pairs[p].Q_self, pairs[p].Q_nbr,
+                                   abs_an_jump_global);
+       real_t scale = 0.0;
+       for (int c = 0; c < NUM_STATE; c++)
+       {
+          scale = std::max(scale, std::abs(F_up[c]));
+       }
+       if (scale == 0.0) { continue; }
+       for (int c = 0; c < NUM_STATE; c++)
+       {
+          real_t obs = F_up[c] - F_ce[c];
+          real_t exp_ = 0.5 * abs_an_jump_global[c];
+          worst_rel_g2 = std::max(worst_rel_g2,
+                                   std::abs(obs - exp_) / scale);
+       }
+    }
+ TEST_LE(worst_rel_g2, 1.0e-12,
+         "R-1204 algebraic identity to 1e-12 relative (per plan)");
```

**Test case:**
The test_godunov_central_flux test itself, with the relative tolerance change above. No new test needed.

---

## Summary
- Critical issues: 1 (R-001)
- Moderate issues: 1 (R-002)
- Low issues: 1 (R-003)
- Plan compliance: PARTIAL — Phase 1–6 are implemented and pass single-rank tests; the multi-rank correctness contract implicitly required by Phase 4 ("Phase 4 handles both `ComputeSharedFaceFluxRHS` and `ComputeADERSharedFaceFluxRHS` with the same dispatch logic" — Risk Assessment R5) is **not** met. Plan §Risk Assessment R5 anticipated this exact failure mode ("MPI shared-face inconsistency") and proposed an `R-1003-style abort guard if nprocs > 1 AND use_substep AND mixed_flux != none may be needed if the testing reveals inconsistency". No abort guard exists, and no MPI test was added that could surface the inconsistency.
- Verdict: **PASS WITH FIXES** — single-rank functionality is correct and well-tested. Before any production multi-rank `--mixed-flux adjacent` run, R-001 must be fixed (or an `nprocs > 1 && Adjacent` abort guard must be installed) and R-002 must be addressed with a parallel test that catches the bug. R-003 can be deferred.

## Unreviewed Areas

The following code paths are within scope but were inspected only briefly — bugs may remain:

- `wave_operator.inl:3340–3408` — fault-face-averaging avg-mode path inside `ComputeADERFaceFluxRHS`. This branch uses `flux_.Interior(can_n, Q_imp_*_g, Q_imp_*_g, ...)` (lines 3371, 3374) with `Q_self == Q_nbr` (Pelties per-side), so Interior == Central at zero jump and mixed-flux dispatch is moot. Verified by inspection but not by test.
- `wave_operator.inl:1480–2110` (RK4 `ComputeFaceFluxRHS` non-fault block) — the surrounding fault-classification and BC-classification branches around the dispatch site at L2034. Read and confirmed the dispatch is correctly nested inside `else { use_precomputed_face_fluxes_ }` so it only fires on the legacy runtime path; not a separate review pass on the BC branches themselves.
- `Makefile` Mixed-Flux test targets — three new build recipes (`seas_test_godunov_central_flux`, `seas_test_mixed_flux_face_set`, `seas_test_mixed_flux_dispatch_none`). Verified they compile and run on macOS + mfem-dev environment; not audited for missing object dependencies that would surface only on a clean Linux/Frontera build.
- Substep-iterator + Mixed-Flux interaction. Plan §Risk Assessment R5 raises this; no test exercises the cross-product. The fault dispatch (substep iterator path) and bulk non-fault dispatch (mixed-flux) are formally independent, but a regression suite that runs `--mixed-flux adjacent --fault-iterator substep` together would be prudent before TPV104 production submission.
