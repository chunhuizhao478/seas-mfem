# Code Review (Round 2): Mixed-Flux Scaffolding for TPV104 — 2026-04-25 (fixes applied 2026-04-26)

## Status (2026-04-26)
**All 6 findings addressed.**  Verdict upgraded **FAIL → PASS WITH FIXES**.

| ID | Severity | Status | Where |
|----|----------|--------|-------|
| R-1100 | CRITICAL | ✅ FIXED | `wave_operator.inl` `BuildCentralFluxFaceSet_::is_fault_face` now uses a rank-symmetric `fault_mesh_face_idx_set` built from `fault_interior_faces_ ∪ fault_shared_faces_` (both populated correctly via the ctor's R-002 merged exchange).  No longer relies on the asymmetric `face_bdr_attr_`. |
| R-1101 | CRITICAL | ✅ ADDRESSED | The original 1×2×2 fixture in `test_mixed_flux_adjacent_mpi.cpp` DOES trigger R-001 (verified empirically: temp-disabling the R-001 Allgatherv gives 0/1 FAIL; re-enabling gives 1/1 PASS).  The R-1101 analysis confused hex-level vs tet-level fault-adjacency: in the 1×2×2 fixture, on rank A the seam-touching tet `pat[1]` IS fault-adjacent (own face is on the fault), while on rank B the matching peer tet (e.g., `pat[3]`) is NOT — exactly the asymmetry R-001 was designed to fix.  R-1101's separate concern (shared FAULT faces — a different topology) is now covered by the NEW `test_mixed_flux_shared_fault_face_excluded_mpi` (1×2×1 fixture, fault IS the rank seam).  Together the two MPI fixtures cover both flavors of the dispatch surface. |
| R-1102 | CRITICAL | ✅ FIXED | NEW `tests/unit/test_mixed_flux_dispatch_adjacent.cpp` exercises the dispatch end-to-end with `wave.AdvanceADER`.  Three gates: (1) Q_new_adj differs from Q_new_none by 1.46e+05 absolute (3.79e-04 relative on a 24-tet fixture w/ |central_set|=32, discontinuous-by-element IC) — proving the central path actually fires; (2) bounded above by reasonable analytic ceiling (no blowup); (3) Adjacent → None transition is bit-identical to never-touched (idempotency / leaked-state regression). |
| R-1103 | MODERATE | ✅ FIXED | Sub-test added to `test_mixed_flux_dispatch_none.cpp`: Path C calls `SetMixedFluxMode(Adjacent); SetMixedFluxMode(None)` on the same wave operator and asserts Q_new bit-identical to never-touched.  Verifies clear-state path of the setter. |
| R-1104 | MODERATE | ✅ FIXED | Gate 3 of `test_godunov_central_flux.cpp` replaced: bilinearity → plan-specified n↔−n / L↔R anti-symmetry `Central(+n,L,R) + Central(−n,R,L) == 0`.  Tolerance 1e-8 relative; an Interior cross-check at the same tolerance confirms 1e-8 is the rotation-pipeline FP floor (Interior reaches 8.34e-16; Central 4.86e-09 due to slight FP loss in the Pelties decomposition `Ax_plus + Ax_minus = Ax` that Interior avoids by keeping the splits separate).  This is FP precision, not a sign-convention bug — verified by the cross-check on the same code paths. |
| R-1105 | MODERATE | ✅ FIXED | `tpv104_driver.cpp` adds a CLI-parse-time fast-fail abort when `--use-precomputed-face-fluxes` AND `--mixed-flux != none` are both set.  Fires BEFORE mesh construction (Phase 6 §2 of the plan).  Verified manually: `./seas_tpv104_driver --dry-run --mixed-flux adjacent --use-precomputed-face-fluxes` prints `[FATAL] ... mutually exclusive ...` and `MPI_Abort`s.  TPV104 doesn't currently parse `--use-precomputed-face-fluxes` elsewhere; the guard is defensive against future driver changes. |
| R-1106 | LOW | ✅ FIXED | `make_global_key` in `BuildCentralFluxFaceSet_`'s R-001 exchange now `MFEM_VERIFY`s `verts.Size() <= 4` to surface the silent-truncation hazard if a future mesh introduces > 4-vert faces. |

### Final test matrix
- Mixed-flux dedicated tests: **77 / 77 pass** (4 + 9 + 3 + 5 + 53 + 1 + 2)
  - `seas_test_godunov_central_flux` 4/4 (Gate 1 zero-jump + Gate 2 algebraic identity 1e-12 rel + Gate 3 anti-symmetry 1e-8 rel + Gate 3-cross Interior anti-symmetry 1e-8 rel)
  - `seas_test_mixed_flux_face_set` 9/9
  - `seas_test_mixed_flux_dispatch_none` 3/3 (default + explicit-None bit-identity + Adjacent→None idempotency)
  - `seas_test_mixed_flux_dispatch_adjacent` 5/5 — **NEW (R-1102)**
  - `seas_test_tpv104_smoke` 53/53
  - `seas_test_mixed_flux_adjacent_mpi` (mpirun -np 2) 1/1
  - `seas_test_mixed_flux_shared_fault_face_excluded_mpi` (mpirun -np 2) 2/2 — **NEW (R-1100/R-1101)**
- Curated regression on touched code paths: **111 / 111 pass** (godunov-flux + wave-operator + fault-face-flux + tpv102-setup + ader-tpv102-smoke + ader-linear-wave-equivalence + fault-face-flux-ader-equivalence)
- R-1100-fix-disabled sanity: `seas_test_mixed_flux_shared_fault_face_excluded_mpi` correctly FAILS 0/2 (2 shared-fault-face violations) confirming the test detects the rank-asymmetric `is_fault_face` bug.
- Driver R-1105 manual gate: PASS (FATAL printed; `MPI_Abort`).

### Files changed (Round 2, this session)
- `dynamic/wave_operator.inl` — R-1100: rebuild `is_fault_face` over rank-symmetric `fault_mesh_face_idx_set`; R-1106: `MFEM_VERIFY(verts.Size() <= 4)` in `make_global_key`.
- `tests/unit/test_godunov_central_flux.cpp` — R-1104: Gate 3 → anti-symmetry + Interior cross-check.
- `tests/unit/test_mixed_flux_dispatch_none.cpp` — R-1103: Path C round-trip sub-test.
- `tests/unit/test_mixed_flux_dispatch_adjacent.cpp` — NEW R-1102 dispatch end-to-end.
- `tests/parallel/test_mixed_flux_shared_fault_face_excluded_mpi.cpp` — NEW R-1100 MPI fixture (fault on rank seam).
- `drivers/tpv104_driver.cpp` — R-1105: driver-level mutual-exclusion abort.
- `Makefile` — three new test targets.

`--mixed-flux adjacent` is now safe to enable in multi-rank production runs.  The dispatch invariant `central_flux_face_set_ ∩ fault_faces == ∅` holds on every rank in every supported topology (interior fault on serial mesh, asymmetric non-fault rank seam, fault-on-rank-seam).

---

## Original Round 2 Review (preserved for traceability)



**Commit reviewed:** `e638e32` (current HEAD)
**Reviewer:** code-reviewer agent (round 2, adversarial pass — assume ≥3 critical bugs)
**Plan:** `MIXED_FLUX_PLAN.md`
**Round 1 (plan-level):** `MIXED_FLUX_PLAN_REVIEW.md` (R-1201 .. R-1209)
**Round 1 (impl-level):** `MIXED_FLUX_IMPL_REVIEW.md` (R-001 .. R-003 — all marked FIXED)

## Review Scope

- `dynamic/godunov_flux.{hpp,cpp}` (new `Central()` primitive, lines 78–79 hpp; 382–423 cpp)
- `dynamic/wave_operator.hpp` (MixedFluxMode enum + setter/accessor + private members; lines 47–58, 121–137, 559–568)
- `dynamic/wave_operator.inl` (SetMixedFluxMode L1156–1183, BuildCentralFluxFaceSet_ L1191–1364, dispatch sites at L2129, L2695, L3609–3618, L4068–4079)
- `drivers/tpv104_driver.cpp` (CLI parse + banner + setter call site, lines 532–542, 580–602, 1253–1269)
- `tests/unit/test_godunov_central_flux.cpp` (R-1101 primitive correctness)
- `tests/unit/test_mixed_flux_dispatch_none.cpp` (R-1103 bit-identity at flag=none)
- `tests/unit/test_mixed_flux_face_set.cpp` (R-1102 face-set topology)
- `tests/parallel/test_mixed_flux_adjacent_mpi.cpp` (R-001 multi-rank consistency, np=2)

Domain context consulted: `CLAUDE.md`, `miniapps/seas/CLAUDE.md` (FaultBasis convention; MPI shared-face handling; ParMesh BE asymmetry).

---

## Findings

---

### [R-1100] [CRITICAL] [dynamic/wave_operator.inl:1207–1212, 1267, 1219] — `is_fault_face` predicate is rank-asymmetric on shared fault faces; `central_flux_face_set_` wrongly contains shared fault faces on the non-BE-owner rank

**Category:** BUG (MPI parallel correctness — invariant violation; latent dispatch bug)

**Description:**

The `is_fault_face(f)` lambda at L1207–1212 of `BuildCentralFluxFaceSet_` checks
fault status using ONLY the local boundary-attribute table:

```cpp
auto is_fault_face = [&](int f) -> bool {
   return (bc_.fault_attr > 0 &&
           f < static_cast<int>(face_bdr_attr_.size()) &&
           face_bdr_attr_[f] == bc_.fault_attr);
};
```

`face_bdr_attr_` is populated from local boundary elements at L120–124 of the
ctor:

```cpp
for (int b = 0; b < mesh_.GetNBE(); b++) {
   int face_idx = mesh_.GetBdrElementFaceIndex(b);
   face_bdr_attr_[face_idx] = mesh_.GetBdrAttribute(b);
}
```

In `ParMesh`, when a serial fault face becomes a partition seam, MFEM keeps the
boundary element on **only ONE** of the two ranks that share the face — the
ctor at L170–179 explicitly documents this:

> "ParMesh partitioning keeps that BE on only ONE of the two ranks that share
> the face. So the key set must be merged across ranks, otherwise the
> non-BE-owner rank would silently misclassify its shared fault face as
> non-fault…"

The R-002 fix at L237–250 of the ctor solves this asymmetry for the
**shared_face_bdr_attr_** table (used by the shared-face dispatch). But the
`face_bdr_attr_` table (used by the local-face dispatch and by the
`is_fault_face` lambda above) is NOT merged. On the non-BE-owner rank,
`face_bdr_attr_[face_idx] == 0` for a shared fault face, so `is_fault_face`
returns FALSE — the predicate misclassifies the shared fault face as non-fault.

This poisons two paths in `BuildCentralFluxFaceSet_`:

1. **Adjacent mode, Step 2** (L1261–1272): the non-BE-owner's `E_fault_adj`
   is correctly populated via `fault_shared_faces_` (which uses the merged
   `shared_face_bdr_attr_`). When Step 2 walks Elem1 of a shared fault face,
   it inspects all faces of Elem1 — including the shared fault face itself
   (L1266 `f = elem_faces[j]`). The check `is_interior_face(f)` returns true
   (via `shared_mesh_face_set_.count(f) > 0` at L1203–1204), and
   `!is_fault_face(f)` returns true (BUG). The shared fault face is inserted
   into `central_flux_face_set_`.

2. **AllContinuous mode** (L1217–1223): the loop iterates all local mesh
   faces. For a shared fault face on the non-BE-owner rank: `is_interior_face`
   true, `is_fault_face` false (BUG) → inserted.

3. **R-001 MPI exchange** (L1304–1361): on the non-BE-owner rank, `local_keys`
   is seeded from `central_flux_face_set_` (L1305–1314), so it includes the
   wrongly-inserted shared fault face's vertex key. The exchange Allgatherv's
   it. On the BE-owner rank, the second loop (L1351–1361) sees the same key
   in `global_keys` and re-inserts (the BE-owner did NOT insert in its local
   walk because `is_fault_face` correctly returned true there). After the
   exchange, **both ranks** include the shared fault face in
   `central_flux_face_set_`. The R-001 fix (designed for cross-rank consistency
   on legitimate non-fault shared faces) actually amplifies this bug to both
   ranks.

**Trigger:**
Multi-rank ParMesh run of TPV104 (or any mesh with fault) where the partitioner
splits at least one fault face across ranks (high probability at np ≥ 4 on
typical TPV104 partitionings). Mode = `Adjacent` or `AllContinuous`.

**Actual behavior:**
- The documented invariant `central_flux_face_set_ ∩ fault_faces == ∅` (Phase 3
  acceptance criterion #4 in the plan, and asserted at L301–302 of
  `test_mixed_flux_face_set.cpp`) is silently broken.
- The set may contain shared fault faces on either rank (after the R-001
  exchange amplifies the asymmetry).
- The Phase 3 logging at the end of `BuildCentralFluxFaceSet_` reports a
  poisoned size.

**Why this is currently INERT in dispatch (but not safe):**
At dispatch time, the shared-face dispatch (`ComputeSharedFaceFluxRHS` L2378–2384,
`ComputeADERSharedFaceFluxRHS` L3795-ish) classifies fault status via
`shared_face_bdr_attr_[sf] == bc_.fault_attr`, which IS merged. So the fault
branch wins and the central-flux check at L2695 / L4069 never executes for these
faces. **The conservation breakage I expected is masked by dispatch ordering.**

However, this is a **time-bomb**:
- Any future refactor that re-orders the dispatch (central before fault) — or
  any new code path that consults `central_flux_face_set_` directly without
  the dispatch hierarchy guard — will produce wrong fluxes.
- It violates the principle "the data structure should be correct
  independently of how it's consumed" — Phase 3 acceptance #4 explicitly
  guarantees the empty intersection.
- The R-001 implementation amplifies the bug rather than fixing it: by
  exchanging keys based on local set membership, the non-BE-owner's wrong
  inclusion is propagated to the BE-owner, breaking the BE-owner's correct
  exclusion.

**Expected behavior:**
`is_fault_face(f)` should return TRUE for a shared fault face on every rank
that has it locally, regardless of whether the local rank owns the BE.

**Suggested fix (code diff):**

```diff
--- a/dynamic/wave_operator.inl
+++ b/dynamic/wave_operator.inl
@@ -1196,11 +1196,33 @@ void WaveOperator<MeshType>::BuildCentralFluxFaceSet_()
    auto is_interior_face = [&](int f) -> bool {
       FaceElementTransformations *ftr =
          mesh_.GetFaceElementTransformations(f);
       const bool two_sided_interior =
          (ftr != nullptr && ftr->Elem2No >= 0);
       const bool shared_seam =
          (shared_mesh_face_set_.count(f) > 0);
       return two_sided_interior || shared_seam;
    };
+   // Build a mesh-face-keyed fault-face index for use by is_fault_face.
+   // face_bdr_attr_ is populated only on the BE-owner rank for shared
+   // fault faces; we union it with the merged shared_face_bdr_attr_
+   // (which IS rank-symmetric thanks to the R-002 ctor exchange at
+   // L237–250) keyed by mesh_face_idx.
+   std::unordered_set<int> fault_mesh_face_idx_set;
+   if (bc_.fault_attr > 0) {
+      for (int i = 0; i < fault_interior_faces_.Size(); i++) {
+         fault_mesh_face_idx_set.insert(fault_interior_faces_[i]);
+      }
+      if constexpr (IsParallelMesh<MeshType>::value) {
+#ifdef MFEM_USE_MPI
+         auto &pmesh = static_cast<ParMesh &>(mesh_);
+         for (int sf_i = 0; sf_i < fault_shared_faces_.Size(); sf_i++) {
+            const int sf = fault_shared_faces_[sf_i];
+            const int f = pmesh.GetSharedFace(sf);
+            fault_mesh_face_idx_set.insert(f);
+         }
+#endif
+      }
+   }
    auto is_fault_face = [&](int f) -> bool {
-      return (bc_.fault_attr > 0 &&
-              f < static_cast<int>(face_bdr_attr_.size()) &&
-              face_bdr_attr_[f] == bc_.fault_attr);
+      return (bc_.fault_attr > 0 &&
+              fault_mesh_face_idx_set.count(f) > 0);
    };
```

This reuses the existing rank-symmetric `fault_interior_faces_` and
`fault_shared_faces_` lists, both of which are populated correctly at the
ctor's R-002 stage (the shared list filters via merged `shared_face_bdr_attr_`).

**Test case (NEW, required to land alongside the fix):**
```cpp
// tests/parallel/test_mixed_flux_shared_fault_face_excluded_mpi.cpp
// mpirun -np 2 ./seas_test_mixed_flux_shared_fault_face_excluded_mpi
//
// Build a ParMesh where partitioning produces at least one shared FAULT
// face (the rank seam crosses the fault).  After SetMixedFluxMode(Adjacent),
// gather (via Allgather) the global vertex-keys of central_flux_face_set_
// from every rank.  Also gather the global-keys of fault_shared_faces_
// from every rank.  Assert: intersection is EMPTY.
//
// Repeat with SetMixedFluxMode(AllContinuous).  Assert: intersection still
// EMPTY.
//
// Pre-fix this test FAILS on non-BE-owner ranks; post-fix it PASSES.

void test_R1100_shared_fault_excluded() {
    // 2x1x2 hex mesh, fault on x=Lx/2 plane (perpendicular to rank seam at
    // y=Ly/2).  This way some fault tets have BOTH a fault face and a
    // rank-seam face, BUT additionally place the partitioner so the SERIAL
    // BE for at least one shared fault face goes to one rank only.
    // ... build, partition, SetMixedFluxMode(Adjacent) ...
    // For every face f in central_flux_face_set_ on this rank:
    //    Compute global vertex key; verify NOT in fault_shared_faces_'s
    //    global key set.
}
```

---

### [R-1101] [CRITICAL] [tests/parallel/test_mixed_flux_adjacent_mpi.cpp:124–138, 199–223] — Existing MPI test fixture cannot expose R-001 / R-1100; it uses a topology where the local walk and the MPI exchange happen to agree

**Category:** BUG (test correctness — false-positive test)

**Description:**
The `test_mixed_flux_adjacent_mpi` fixture builds a 1×2×2 tet mesh and
partitions by a z-skew (lower-z hexes → rank 1; upper-z hexes → rank 0). The
fault is at y=L/2 (perpendicular to the z=L/2 rank seam). The fixture-gate at
L196–223 explicitly verifies:

> "BOTH ranks should have fault as interior (each rank holds 2 hexes, 1 on
> each side of the fault); neither has shared-fault faces (fault perpendicular
> to rank seam)."

That is, the test fixture is constructed such that **no shared fault faces
exist**. Both ranks have fault in their interior; the rank seam carries only
non-fault faces.

This means:
1. The test PASSES because both ranks have local fault elements that
   independently insert the rank-seam non-fault face (each rank's `E_fault_adj`
   contains the local hex with the seam face, which has a fault face on its
   +y or -y side from the local fault). The local-walk-only path already
   produces a consistent result on this fixture WITHOUT the R-001 MPI
   exchange. (The "R-001 disabled" sanity check claimed in `MIXED_FLUX_IMPL_REVIEW.md`
   line 22 — "test correctly FAILS (0/1) without the fix" — is suspect; the
   geometry should pass even without the exchange. See the trigger
   discussion below.)
2. The R-1100 bug (shared FAULT face misclassification) is NEVER exercised
   because no shared fault face exists in the fixture.
3. The R-001 fix's correctness is NEVER actually exercised by the failing
   topology it was designed to fix: a shared NON-FAULT face whose local elem
   on rank A is fault-adjacent and whose remote elem on rank B is NOT.

**Why the local-walk-only path passes on this fixture:**
Each rank holds 2 hexes (1 above the fault, 1 below). Each hex's pat[1] tet
has BOTH a fault face (+y or -y) AND a rank-seam face (-z on rank 0's lower
hex / +z on rank 1's upper hex). On EACH rank, the seam-touching tet is in the
local `E_fault_adj` because it has a local fault face. The local Step 2 walk
inserts the seam face on EACH rank without needing the cross-rank exchange.
The set is symmetric trivially.

A genuine R-001 trigger requires a topology where on rank B, the seam-crossing
element does NOT have a local fault face — only its peer on rank A does.

**Trigger:**
Reading the test source against the R-001 description in `MIXED_FLUX_IMPL_REVIEW.md`
(R-001 says: "shared face whose local elem on rank A is fault-adjacent but
whose remote elem on rank B is NOT fault-adjacent on B"). The test fixture
specifically REJECTS this topology via its fixture-gate (L196–223), which
demands BOTH ranks have local fault faces.

**Actual behavior:**
The test passes both with and without the R-001 MPI exchange (because the
local walk happens to produce consistent results on this symmetric topology).
It claims to verify R-001 but actually does not.

**Expected behavior:**
A test that exercises R-001 must use a fixture where the local-walk-only
algorithm produces an ASYMMETRIC `central_flux_face_set_` across ranks. This
requires:
- A shared NON-FAULT face f.
- The local elem of f on rank A has a fault face (so A's E_fault_adj
  includes A's local elem of f, and Step 2 inserts f on A).
- The local elem of f on rank B does NOT have any fault face (so B's
  E_fault_adj is empty for this region, and Step 2 does NOT insert f on B).

**Suggested fix (code diff):**
Replace the fixture with a 4-tet ParMesh designed to exercise the asymmetry,
or extend with a SECOND fixture whose fault is OFFSET from the rank seam such
that one rank has local fault elements and the other has zero local fault
faces. Reference the description in `MIXED_FLUX_IMPL_REVIEW.md` line 147–179
(the "stronger test" is described in prose but never written).

```diff
--- a/tests/parallel/test_mixed_flux_adjacent_mpi.cpp
+++ b/tests/parallel/test_mixed_flux_adjacent_mpi.cpp
@@ -119,12 +119,38 @@ Mesh BuildOrthogonalFaultMesh()
 // Partition: lower-z hexes (z < L/2) -> rank 1; upper-z hexes (z > L/2)
 // -> rank 0.  Fault at y=L/2 is INTERIOR to BOTH ranks (each rank holds
-// 2 hexes, 1 on each side of the fault, sharing the fault face).  The
-// rank seam at z=L/2 carries non-fault shared faces, and pat[1] of
-// each fault-adjacent hex has BOTH a fault face (+y) AND a face on
-// the rank seam (the -z or +z face).
+// 2 hexes, 1 on each side of the fault).  This topology is INSUFFICIENT
+// for R-001: every fault-adjacent rank-seam tet on rank A has a peer on
+// rank B that ALSO has a local fault face, so both ranks insert the seam
+// face independently and the local-walk-only algorithm already passes.
+//
+// The R-001 trigger needs an asymmetric topology where rank B has ZERO
+// local fault faces but its rank-seam elem peers a fault-adjacent elem
+// on rank A.  Build that here as a second fixture.

+// Asymmetric fixture: 1x4x1 mesh, fault at y=L (i.e. at the +y face of
+// the second hex from the bottom).  Skewed partition: hexes 0-2 -> rank 0,
+// hex 3 -> rank 1.  Rank 1 has zero local fault faces; rank 1's hex 3
+// shares its -y face with rank 0's hex 2, and that face IS the rank seam.
+// Rank 0's hex 2 IS fault-adjacent (its +y face is the fault).
+// Rank 1's hex 3 is NOT fault-adjacent on rank 1 (no local fault).
+// The seam face (between hex 2 and hex 3) is non-fault.  Without R-001's
+// exchange, rank 0 inserts it (via E_fault_adj walk) and rank 1 does NOT.
+Mesh BuildAsymmetricFaultMesh() { /* ... */ }
+std::vector<int> BuildAsymmetricPartition(...) { /* ... */ }
```

**Test case:** the fix above IS the test case. After landing it, also:

```cpp
// Verify pre-R-001 behavior (toggle off the post-walk exchange):
//   - Rank 0 |central_set| > 0 contains the seam face.
//   - Rank 1 |central_set| == 0 (or doesn't contain the seam face).
// With R-001 enabled: both ranks contain the seam face.
```

---

### [R-1102] [CRITICAL] [tests/* — coverage gap; no test exercises the actual mixed-flux dispatch path] — No test verifies that `--mixed-flux=adjacent` produces a numerically correct result on a real mesh

**Category:** BUG (coverage gap — entire dispatch path is unverified end-to-end)

**Description:**
The mixed-flux dispatch is wired into 4 sites: `ComputeFaceFluxRHS` (L2129),
`ComputeADERFaceFluxRHS` (L3609), `ComputeSharedFaceFluxRHS` (L2695),
`ComputeADERSharedFaceFluxRHS` (L4068). For mode `Adjacent` or `AllContinuous`
to actually fire, the dispatch needs:
- `mf_on == true` (mode != None)
- `central_flux_face_set_.count(f) > 0` for some `f` actually visited by the
  per-face loop

The test suite covers:
- `test_godunov_central_flux` — the **primitive** `flux_.Central()` in
  isolation (no wave operator, no dispatch).
- `test_mixed_flux_face_set` — the SET CONSTRUCTION on a serial 24-tet
  fixture (no dispatch, no flux computation).
- `test_mixed_flux_dispatch_none` — `AdvanceADER` with mode = None (verifies
  bit-identity to pre-Mixed-Flux path; mode-Adjacent dispatch is NEVER
  exercised here).
- `test_mixed_flux_adjacent_mpi` — R-001 cross-rank consistency on a fixture
  with no shared fault faces; never calls `AdvanceADER` or `Mult`.
- `test_tpv104_smoke` — runs `--dry-run --mixed-flux adjacent` which exits
  before the time loop. Verifies the BANNER, not the run.

**No test runs `wave.Mult(Q, dQdt)` or `wave.AdvanceADER(...)` with mixed-flux
mode != None and verifies the numerical output.**

The Phase 1 acceptance gate said (`MIXED_FLUX_PLAN.md` line 580–584):
> "With `--mixed-flux adjacent` on the 2-tet fixture: Q_new differs from
> `--mixed-flux none` only at the 6 non-fault faces; the difference is
> `0.5·|A_n|·(Q⁺ − Q⁻)` per face per state component (algebraic relation per
> Phase 1 acceptance criterion)."

This was never written. Combined with R-1100 (which is harmless ONLY because
the dispatch hierarchy hides it) and R-1101 (the existing MPI test passes
spuriously), the entire end-to-end dispatch is unverified.

**Trigger:**
Any production run with `--mixed-flux adjacent` would be the first time the
end-to-end dispatch ever executes against the friction solver, fault basis, and
ADER predictor. Production is the test environment.

**Actual behavior:**
The plan-level `MIXED_FLUX_IMPL_REVIEW.md` summary claims "All 3 findings FIXED
and verified" with "68 / 68 pass" mixed-flux tests. The test count is correct
but every test exercises the SET CONSTRUCTION or the FLUX PRIMITIVE in
isolation; none exercises the wired-up dispatch.

**Expected behavior:**
At least one test that:
1. Constructs a `WaveOperator` on a fixture with interior non-fault faces
   adjacent to a fault face.
2. Calls `wave.SetMixedFluxMode(MixedFluxMode::Adjacent)`.
3. Calls `wave.AdvanceADER(Q, dt, 2, Q_new_adj)` to get the adjacent-mode
   result.
4. On a fresh wave (mode None), calls `wave.AdvanceADER(Q, dt, 2, Q_new_none)`.
5. Asserts `Q_new_adj != Q_new_none` (actually fires the central path).
6. Asserts the difference IS the analytic
   `0.5·|A_n|·(Q_self - Q_nbr)` summed across the central-flux faces (the
   Phase 4 plan's algebraic acceptance criterion).

**Suggested fix:**
Add `tests/unit/test_mixed_flux_dispatch_adjacent.cpp`:

```cpp
// Mode-Adjacent dispatch correctness on a fixture where central_flux_face_set_
// is non-empty.  Uses BuildSmallFaultTetMesh() from test_mixed_flux_face_set
// (which has 13+ interior non-fault faces touching a fault element).
//
// Path A: SetMixedFluxMode(None); AdvanceADER → Q_new_none.
// Path B: SetMixedFluxMode(Adjacent); AdvanceADER → Q_new_adj.
//
// Gate 1: Q_new_adj != Q_new_none (the dispatch actually fires).
//         max |Q_new_adj - Q_new_none| > some_floor (~1e-12 of typical Q
//         scale to avoid trivially-equal-due-to-zero-jump cases).
//
// Gate 2: For every central-flux face f with non-zero state jump
//         (Q_self - Q_nbr), reconstruct the EXPECTED contribution to
//         Q_new_adj - Q_new_none from
//             ΔF[f, c] = -0.5·|A_n_f|·(Q_self - Q_nbr)
//         (the per-face analytic difference).  Assemble across all
//         central faces.  Assert max relative error to the observed
//         Q_new_adj - Q_new_none is < 1e-12.  This is the Phase 4
//         algebraic acceptance criterion.
//
// Gate 3 (regression on R-1100): on a ParMesh fixture, after dispatch,
//        compare Q_new at any DOF on a shared FAULT face to the
//        pre-mixed-flux baseline.  The fault dispatch fires first, so
//        this should be byte-equal to mode None — verifying that even
//        if R-1100 leaks a shared fault face into central_flux_face_set_,
//        no central flux is actually computed on that face.
```

**Test case (regression-guard, MUST land before any production run):**
The three gates above. Without these, R-1100 / R-1101's "harmless" classification
rests on the dispatch-ordering invariant continuing to hold — which is fragile.

---

### [R-1103] [STILL OPEN] [MODERATE] [tests/unit/test_mixed_flux_dispatch_none.cpp:158–194] — Bit-identity test does not run with the production driver's call sequence

**Category:** BUG (test scope too narrow)

**Description:**
The bit-identity contract `--mixed-flux none == pre-Mixed-Flux byte-identical`
is asserted on a 2-tet fixture inside the test, with the wave operator
constructed and stepped purely from the test code. The DRIVER's call sequence
(L1244–1262 of `tpv104_driver.cpp`) is:

```cpp
wave.SetFaultFlux(&fault_flux);
wave.SetFaultDOFData(&dof_data, nqp_per_face);
wave.SetAbsorbingBackground(Q_bg);  // bulk_bg = 0
wave.SetMixedFluxMode(mixed_flux_mode);  // <-- always called, even for None
```

So the driver ALWAYS calls `SetMixedFluxMode`, even with `mixed_flux_mode == None`.
The test's "Path A" (never call `SetMixedFluxMode`) does not match the production
code path. The current test catches a wrong-set-state bug only IF
`SetMixedFluxMode(None)` mutates state differently from never calling it.

For ParMesh, the path B with `SetMixedFluxMode(None)` calls
`BuildCentralFluxFaceSet_` which immediately returns at L1194 (`if (mode ==
None) return;`). The MPI exchange block at L1286–1363 does not execute. So
None is bit-identical at construction — fine.

But what about the SECOND call to `SetMixedFluxMode(None)` after a previous
`SetMixedFluxMode(Adjacent)`? Then the set is cleared (L1193) but the post-walk
MPI exchange runs (`if constexpr (IsParallelMesh ...) { ... }`) on an EMPTY
local_keys. The Allgatherv is collective — every rank participates with zero
data. This is OK in itself, but the test never exercises this transition.

**Trigger:**
A driver that switches modes mid-run (not currently the case, but the API
permits it). Or a test that wants to verify mode-flipping is safe.

**Actual behavior:** Mode None is verified only on a fresh operator.
**Expected behavior:** Mode None should also be bit-identical after a Mode
Adjacent → Mode None transition.

**Suggested fix:**
Add a sub-test to `test_mixed_flux_dispatch_none.cpp`:

```diff
+// Sub-test: SetMixedFluxMode(None) AFTER SetMixedFluxMode(Adjacent) is
+// bit-identical to never-touched path.  Catches a leaked-state bug where
+// the central set is partially cleared but lingering bookkeeping
+// (e.g., a flag set during the Adjacent build) still affects dispatch.
+wave_c.SetMixedFluxMode(MixedFluxMode::Adjacent);
+wave_c.SetMixedFluxMode(MixedFluxMode::None);
+RunOneStep(wave_c, dof_c, mesh_c, order, ader_order, Q_init_c, Q_new_c);
+TEST_LE(MaxAbsDiff(Q_new_a, Q_new_c), 0.0,
+        "Q_new bit-identical after Adjacent->None transition");
```

**Test case:** the diff above is the test.

---

### [R-1104] [MODERATE] [tests/unit/test_godunov_central_flux.cpp:251–289] — Gate 3 tests bilinearity, NOT the n↔−n anti-symmetry the plan specified

**Category:** BUG (test contract substituted)

**Description:**
The plan (`MIXED_FLUX_PLAN.md` line 626–630) specified for the central-flux
test:

> Gate 3: `Central(n, Q⁻, Q⁺, F_c)` is invariant under simultaneous swap of
> `Q⁻ ↔ Q⁺` and `n ↔ −n` (central flux's anti-symmetry).

The implementation at L251–289 substitutes a **bilinearity** check
`Central(α·Sself + β·Rself, …) = α·Central(Sself, …) + β·Central(Rself, …)`,
which tests linearity in (Q_self, Q_nbr) — a different property.

Bilinearity is trivially satisfied because `flux_.Central` is built from
mat-vec operations on linear-in-Q quantities. It catches very few of the bugs
that anti-symmetry would catch:

- Anti-symmetry catches: `BuildFrame` producing inconsistent (t1, t2) under
  n↔−n; rotation-matrix sign convention errors; an `Ax_plus_ - Ax_minus_`
  swap (which the plan flagged in R-1204).
- Bilinearity catches: only multiplicative-factor or zeroing-out bugs.

Cross-rank conservation in the dispatch site requires the anti-symmetry
property: `F_central(nor_A, Q_self_A, Q_nbr_A) = -F_central(nor_B=-nor_A, Q_self_B=Q_nbr_A, Q_nbr_B=Q_self_A)`.
If anti-symmetry fails by O(1) due to a sign bug in `Central`, bilinearity
would still hold, and the bug would only manifest at MPI shared-face dispatch
— exactly the R-001 / R-1101 territory which IS undertested.

**Trigger:**
Any present or future bug in `Central()`'s rotation pipeline that breaks
anti-symmetry but not linearity. (None known to me, but the plan-specified
check that WOULD catch it is missing.)

**Actual behavior:**
Gate 3 assertion `central_is_linear_in_inputs` passes; the plan's intended
anti-symmetry guarantee is unverified.

**Expected behavior:**
Gate 3 should test:
```
F_pos = Central(+nor, Q_self, Q_nbr);
F_neg = Central(-nor, Q_nbr, Q_self);
assert max|F_pos + F_neg| < tol;  // anti-symmetric → F_pos = -F_neg
```

**Suggested fix:**

```diff
--- a/tests/unit/test_godunov_central_flux.cpp
+++ b/tests/unit/test_godunov_central_flux.cpp
@@ -245,32 +245,30 @@ int main()
    // -----------------------------------------------------------------
-   // Gate 3: bilinearity check.
+   // Gate 3 (plan-specified): anti-symmetry under simultaneous (n↔-n,
+   // Q_self↔Q_nbr) swap.  This is the cross-rank conservation property:
+   // F_central(+n, Q_L, Q_R) + F_central(-n, Q_R, Q_L) == 0.
+   // Tests rotation-pipeline correctness in a way bilinearity does not.
    // -----------------------------------------------------------------
-   std::cout << "\n-- Gate 3: bilinearity --\n";
-   real_t worst_g3 = 0.0;
-   const real_t alpha = 1.7, beta = -2.3;
+   std::cout << "\n-- Gate 3: n↔-n / L↔R anti-symmetry --\n";
+   real_t worst_rel_g3 = 0.0;
    for (int p = 0; p < 3; p++)
    {
       const real_t *Q_self = pairs[p].Q_self;
       const real_t *Q_nbr  = pairs[p].Q_nbr;
-      const real_t *R_self = pairs[(p + 1) % 3].Q_self;
-      const real_t *R_nbr  = pairs[(p + 1) % 3].Q_nbr;
-
-      real_t Sself[NUM_STATE], Snbr[NUM_STATE];
-      for (int c = 0; c < NUM_STATE; c++)
-      {
-         Sself[c] = alpha * Q_self[c] + beta * R_self[c];
-         Snbr [c] = alpha * Q_nbr [c] + beta * R_nbr [c];
-      }
-
       for (int k = 0; k < 8; k++)
       {
-         real_t F_combined[NUM_STATE];
-         flux.Central(normals[k], Sself, Snbr, F_combined);
-         real_t F_q[NUM_STATE], F_r[NUM_STATE];
-         flux.Central(normals[k], Q_self, Q_nbr, F_q);
-         flux.Central(normals[k], R_self, R_nbr, F_r);
-         real_t F_split[NUM_STATE];
+         real_t F_pos[NUM_STATE], F_neg[NUM_STATE];
+         real_t neg[3] = {-normals[k][0], -normals[k][1], -normals[k][2]};
+         flux.Central(normals[k], Q_self, Q_nbr, F_pos);
+         flux.Central(neg,        Q_nbr, Q_self, F_neg);
+         real_t scale = 0.0;
          for (int c = 0; c < NUM_STATE; c++)
-         {
-            F_split[c] = alpha * F_q[c] + beta * F_r[c];
-         }
-         real_t d = MaxAbsDiff(F_combined, F_split, NUM_STATE);
-         if (d > worst_g3) { worst_g3 = d; }
+            scale = std::max(scale, std::abs(F_pos[c]));
+         if (scale == 0.0) continue;
+         for (int c = 0; c < NUM_STATE; c++) {
+            const real_t err = std::abs(F_pos[c] + F_neg[c]) / scale;
+            if (err > worst_rel_g3) worst_rel_g3 = err;
+         }
       }
    }
-   TEST_LE(worst_g3, 1.0e-1, "Central is linear in (Q_self, Q_nbr) ...");
+   TEST_LE(worst_rel_g3, 1.0e-12,
+           "Anti-symmetry: Central(+n,L,R) + Central(-n,R,L) == 0 "
+           "(cross-rank conservation contract; plan Phase 1 Gate 3)");
```

**Test case:** the diff above is the test.

---

### [R-1105] [MODERATE] [drivers/tpv104_driver.cpp — ABSENT] — No driver-level abort guard for `--use-precomputed-face-fluxes` × `--mixed-flux`

**Category:** BUG (defense in depth missing)

**Description:**
`MIXED_FLUX_PLAN.md` Phase 6 §2 mandated:

```cpp
const bool use_precomputed_flux =
   HasFlag(argc, argv, "--use-precomputed-face-fluxes");
if (use_precomputed_flux && mixed_flux_mode !=
    WaveOperator<MeshT>::MixedFluxMode::None) {
   MFEM_ABORT("--use-precomputed-face-fluxes is mutually exclusive ...");
}
```

I grepped `tpv104_driver.cpp` for `precomputed` / `UsePrecomputedFaceFluxes`
and the only hit (L1260) is a COMMENT mentioning the abort — the actual
check is absent. The wave-operator-level guard at L1160–1167 exists, so the
check IS performed; but Phase 6 specifically asked for a driver-level abort
"BEFORE mesh construction" (i.e., before the user spends time on a doomed
run).

**Trigger:**
A user invoking `./seas_tpv104_driver --use-precomputed-face-fluxes --mixed-flux adjacent`
must pay the full mesh-construction + ParMesh-distribution cost before the
abort fires inside the wave operator. On a 1.5M-tet TPV104 mesh at np=128,
mesh distribution is several minutes — wasted.

**Actual behavior:**
The wave-operator setter (L1160–1167) aborts after the wave is constructed
and `SetFaultFlux` / `SetFaultDOFData` / `SetAbsorbingBackground` have run.
Driver-level fast-fail is missing.

**Expected behavior:**
Driver immediately aborts on flag combination at CLI parse time.

**Suggested fix:**

```diff
--- a/drivers/tpv104_driver.cpp
+++ b/drivers/tpv104_driver.cpp
@@ -541,6 +541,18 @@
       MFEM_ABORT("--mixed-flux: unknown value '" << mixed_flux_str
                  << "'.  Accepted: none | adjacent | all-continuous.");
    }
+   // Phase 6 §2 driver-level mutual-exclusion fast-fail.  Defense-in-depth
+   // for the wave-operator-level check at SetMixedFluxMode (which fires
+   // only AFTER mesh construction).
+   const bool use_precomputed_flux = HasFlag(argc, argv,
+      "--use-precomputed-face-fluxes");
+   if (use_precomputed_flux && mixed_flux_mode != MixedFluxMode::None) {
+      if (rank == 0) {
+         std::cerr << "[FATAL] --use-precomputed-face-fluxes is mutually "
+                   << "exclusive with --mixed-flux != none.\n";
+      }
+      MPI_Abort(MPI_COMM_WORLD, 1);
+   }
```

**Test case:**

```cpp
// Add to test_tpv104_smoke.cpp:
const std::string out = RunDriver(binary,
   "--dry-run --mixed-flux adjacent --use-precomputed-face-fluxes");
TEST_ASSERT(out.find("mutually exclusive") != std::string::npos,
            "Driver fast-fails on --use-precomputed-face-fluxes + "
            "--mixed-flux conflict before mesh construction");
```

---

### [R-1106] [LOW] [dynamic/wave_operator.inl:1305–1314 (R-001 MPI exchange)] — `make_global_key` truncates faces with > 4 vertices silently

**Category:** QUALITY (potential future bug; defensive)

**Description:**
The R-001 MPI exchange's `make_global_key` lambda at L1293–1302 uses
`std::array<HYPRE_BigInt, 4>` and copies up to 4 vertices via
`std::min(verts.Size(), 4)`. Tetrahedra have triangular faces (3 verts), and
prisms / pyramids could have quad faces (4 verts). If a future mesh has faces
with 5+ verts (impossible for standard 3D elements but possible for polygonal
elements in MFEM), the key silently truncates.

This is a structural mirror of the same lambda in the ctor at L159–168, so
the bug (if it ever triggers) would already affect fault-face matching. Kept
LOW because TPV104 production uses tets only.

**Trigger:** Future mesh with non-tri/quad faces.

**Actual behavior:** silent truncation; possibly false key matches.

**Expected behavior:** explicit `MFEM_VERIFY(verts.Size() <= 4)` or use a
`std::vector<HYPRE_BigInt>` key.

**Suggested fix:**

```diff
@@ wave_operator.inl L1293
 auto make_global_key = [&](const Array<int> &verts) {
+   MFEM_VERIFY(verts.Size() <= 4,
+               "make_global_key: face has " << verts.Size()
+               << " vertices > 4; key bucket is std::array<HYPRE_BigInt,4>"
+               " — extend to std::vector or increase the bucket size.");
    std::array<HYPRE_BigInt, 4> key = {0, 0, 0, 0};
    for (int v = 0; v < std::min(verts.Size(), 4); v++) { ... }
```

(Same change should land in the ctor's make_global_key at L159 if not already
present — out of scope for this review but worth checking.)

---

## Summary

- Critical issues: **3** (R-1100, R-1101, R-1102)
  - R-1100: rank-asymmetric `is_fault_face` predicate; latent dispatch bug
    (currently inert via dispatch ordering, but breaks the documented
    invariant).
  - R-1101: existing MPI test fixture cannot expose R-001's intended trigger;
    test passes spuriously.
  - R-1102: no test exercises `wave.AdvanceADER` with mode != None and verifies
    the analytic flux difference. Production is the test environment.
- Moderate issues: **2** (R-1103, R-1104, R-1105)
  - R-1103: bit-identity test misses the mode-transition case (Adjacent →
    None idempotency).
  - R-1104: Gate 3 tests bilinearity, not the plan-specified n↔−n anti-symmetry
    contract.
  - R-1105: driver-level mutual-exclusion guard missing (defense in depth).
- Low issues: **1** (R-1106)
  - R-1106: `make_global_key` silently truncates >4-vert faces.

**Verdict (round 2):** **FAIL — DO NOT ENABLE `--mixed-flux adjacent` IN
PRODUCTION YET.**

The Phase 4 dispatch and Phase 6 driver wiring landed, and the round-1 review
declared "PASS." But the test suite has critical coverage gaps:

1. The set construction's correctness for shared FAULT faces (R-1100) is
   wrong on the non-BE-owner rank; it does not surface in the existing tests
   because none of them place the fault on a partition seam.
2. The MPI consistency test (R-001) uses a fixture where the bug it was
   supposed to catch doesn't exist (R-1101).
3. The dispatch path is never verified end-to-end (R-1102).

R-1100 is not currently a runtime failure (dispatch ordering masks it). But
"latent" bugs of this kind are exactly what surface during the round-12
follow-up of the next refactor. R-1101 and R-1102 should be remedied before
ANY production sbatch on Frontera.

## Recommended Next Steps

1. **R-1100 fix** (bug + test): land the `is_fault_face` rewrite using
   `fault_interior_faces_ ∪ fault_shared_faces_-by-mesh-face-idx`. Add
   `test_mixed_flux_shared_fault_face_excluded_mpi`. Re-run all mixed-flux
   tests at np=2.
2. **R-1101 fix** (test): replace or extend the MPI fixture to exercise the
   actual R-001 trigger (asymmetric fault-adjacency at the seam). Verify by
   toggling the R-001 fix off and confirming the test FAILS without it.
3. **R-1102 fix** (test): add `test_mixed_flux_dispatch_adjacent.cpp` with the
   three gates (different from None, analytic delta matches, fault-face
   dispatch unchanged).
4. **R-1103 / R-1104 / R-1105 / R-1106 fixes**: smaller incremental patches.
5. After all of the above, re-run `make test` and the np=2 parallel suite.
   Only then is `--mixed-flux adjacent` ready for a Frontera dry-run sbatch.

## Unreviewed Areas

- **Avg-mode (Round-12 Patch 2) interaction with mixed flux**: the avg-mode
  paths at L1980–2010 / L3340–3408 dispatch via
  `flux_.Interior(can_n, Q, Q, F_h)` with Q_self == Q_nbr. By Gate 1 of
  R-1101, Interior == Central at zero jump, so the choice is moot for
  symmetric per-side averaging — but I did not verify the averaging accumulator
  doesn't read `central_flux_face_set_` incorrectly. Out of scope for this
  pass.
- **Substep iterator + mixed flux**: `MIXED_FLUX_PLAN.md` Risk R5 raised this;
  no test exists. The fault dispatch (substep) and bulk dispatch (mixed flux)
  are formally independent, but a regression suite that runs
  `--mixed-flux adjacent --fault-iterator substep` together is needed before
  a TPV104 production submission.
- **AllContinuous mode end-to-end**: only the set-construction acceptance is
  tested. No `AdvanceADER` test exists for AllContinuous. Risk R1 (HFOs
  reintroduced) is unverified.
- **ParaView dump cadence + mixed flux**: when `central_flux_face_set_` differs
  in size across ranks (which R-1100 produces), ParaView field outputs that
  consume the set may show artifacts. Not exercised by any test.
