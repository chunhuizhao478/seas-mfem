# Code Review: PLAN_cross_rank_material_exchange_2026-06-06 (+ codebase generalization for CVM)

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md`
- Code explored: `dynamic/bimaterial_wave_operator.{hpp,inl}` (ExchangeBiMaterialNeighbours_,
  the two shared-face builds, AssignFaultSidePerMaterialImpedances), `dynamic/wave_operator.inl`
  (EvaluateBulkAtFaultQPsCanonical face-nbr template, ComputeADERSharedFaceFluxRHS,
  shared_fault_elem1_on_plus_/shared_fault_dof_offset_ indexing), `dynamic/fault_locality_partition.hpp`,
  the parallel test harness (`tests/parallel/test_bimaterial_mixed_flux_shared.cpp`,
  `test_r101_shared_fault.cpp`).
- Domain: CLAUDE.md, the B0 verdict, the Part-B/C memory.
- Goal restated by user: generalize EVERY part for arbitrary cross-element material
  contrast (CVM), and VERIFY against TPV31 (symmetric depth-profile) AND TPV6 (bimaterial).

## Findings

### [R-001] CRITICAL — Phase 4: per-side fault material on shared faces uses the peer CENTROID, which breaks TPV31 symmetry and is wrong for any depth-varying material

**Category:** BUG / DEVIATION

**Description:**
Phase 4 assigns the SHARED-fault peer side from the per-element CENTROID material delivered
by the Phase-1 exchange (`shared_face_neighbour_material_[sf]`), while the LOCAL side uses the
**eps-offset** material at the fault QP (the existing `offset_ip` helper).  Those are evaluated
at DIFFERENT depths:
- LOCAL eps-offset = material at the fault-QP point (depth z of the station), offset along the
  fault NORMAL (horizontal for a vertical fault) ⇒ stays at depth z.
- PEER centroid = material at the peer element's centroid (depth z' = z ± h/2).

For a depth-profile material (TPV31) the fault is SYMMETRIC (same material both sides), so the
correct result is `Zp_plus == Zp_minus`.  But LOCAL(z) vs PEER(z') gives `Zp_plus != Zp_minus`
⇒ a SPURIOUS bi-material contrast on every shared TPV31 fault face.  The plan's own Phase-4
"Edge Cases" claim ("peer material == local ⇒ Zp_plus == Zp_minus ⇒ byte-exact") is FALSE for
depth-profile: peer centroid != local fault-QP value.

This is a REGRESSION: today (stub) a shared TPV31 fault face keeps the homogeneous seed
(`Zp_plus == Zp_minus`, symmetric, correct); Phase 4 as planned makes it asymmetric (wrong) —
the exact σ_n-leak driver TPV31 must avoid.  It also gives the WRONG per-side material for a
genuine depth-varying bi-material fault (CVM), which is the whole point of this work.

**Trigger:** any run with a fault on a partition seam + a depth-varying material, i.e. TPV31 at
np>1 WITHOUT `--partition-fault-locality` (which the plan deliberately makes optional), or CVM.

**Actual behavior (planned):** shared-fault peer side = peer centroid material (wrong depth).

**Expected behavior:** shared-fault peer side = peer material evaluated at the SAME fault QP
(eps-offset into the peer element), so depth-profile stays symmetric and CVM is correct.

**Suggested fix:** add a dedicated FAULT-QP-resolved peer-material exchange (a NEW phase, or
fold into Phase 4) — do NOT reuse the per-element centroid exchange for the fault per-side:
- On each rank, for every shared fault face, evaluate the LOCAL material at the eps-offset of
  EACH fault QP (this is the same `offset_ip` the interior path uses), producing `nbf_per_face`
  triples per shared fault face.
- Exchange these per-shared-fault-QP triples with the peer (the peer's local side = our Elem2).
  Mirror the per-component shared-fault exchange already used for bulk Q
  (`EvaluateBulkAtFaultQPsCanonical`, wave_operator.inl:2242-2323) — key by `fault_shared_faces_`,
  use `shared_fault_dof_offset_`.
- Assign the peer side from the exchanged per-QP eps-offset material; the LOCAL side from the
  local eps-offset.  Then depth-profile ⇒ both sides at depth z ⇒ `Zp_plus == Zp_minus`
  (byte-exact symmetric); CVM ⇒ correct per-side contrast.
```diff
- // PEER side material: the peer's CENTROID material from shared_face_neighbour_material_[sf]
+ // PEER side material: the peer's EPS-OFFSET material AT THE FAULT QP, from a per-shared-
+ // fault-QP exchange (NOT the per-element centroid — that is a different depth and breaks
+ // depth-profile symmetry / CVM accuracy).
```

**Test case:**
```cpp
// np=2: a DEPTH-PROFILE (symmetric) fault on a partition seam must stay symmetric.
void test_R001_depthprofile_shared_fault_symmetric_np2() {
   // 2 elements, 1/rank, the inter-rank face tagged FAULT; Mode::Coefficient material
   // = function of z only (depth profile), SAME on both sides of the fault.
   // After AssignFaultSidePerMaterialImpedances, the shared fault DOF must have
   // Zp_plus == Zp_minus (to ~1e-12 rel) — NOT a spurious contrast.
   // (Fails under the centroid-for-peer plan; passes with the fault-QP exchange.)
}
```

---

### [R-002] MODERATE — Phase 4: `shared_fault_elem1_on_plus_` is POSITION-indexed but `shared_fault_dof_offset_` is RAW-sf-keyed; the plan conflates them

**Category:** ASSUMPTION / BUG-risk

**Description:**
The plan's Phase-4 loop uses `shared_fault_elem1_on_plus_[si]` and
`shared_fault_dof_offset_[sf]` without flagging that these use DIFFERENT index spaces:
- `shared_fault_elem1_on_plus_` is sized `fault_shared_faces_.Size()` and indexed by the
  POSITION `si` in `fault_shared_faces_` (wave_operator.inl:464, 521 use `sf_idx`).
- `shared_fault_dof_offset_` is a `std::map<int,int>` KEYED BY THE RAW shared-face index `sf`
  (= `fault_shared_faces_[si]`) (wave_operator.hpp:1063; looked up via `.find(sf)` at
  wave_operator.inl:2334).
An implementer who uses `sf` for the on_plus array (or `si` for the dof_offset map) gets a
silent off-by-mapping bug (wrong +/- side or wrong DOF slice).

**Trigger:** implementing Phase 4 from the plan as written.

**Suggested fix:** make the index spaces explicit in the plan and in the code:
```cpp
for (int si = 0; si < fault_shared_faces_.Size(); ++si) {
   const int sf = fault_shared_faces_[si];                 // RAW shared-face index
   const bool e1_plus = shared_fault_elem1_on_plus_[si];   // POSITION-indexed
   auto off_it = shared_fault_dof_offset_.find(sf);        // RAW-sf-keyed
   MFEM_VERIFY(off_it != shared_fault_dof_offset_.end(), "...");
   const int dof_off = off_it->second;
   ...
}
```

**Test case:** the R-001 np=2 fault test also asserts the shared fault DOF lands at the correct
`dof_off` (the DOFData index matches the station's expected location) — a wrong index space
writes to the wrong DOF and the assertion fails.

---

### [R-003] MODERATE — Plan: contradictory / under-specified TPV31 byte-exact contract; the "verify against TPV31" goal needs a precise statement

**Category:** DEVIATION (spec clarity)

**Description:**
The Constraints + Risk R1 say the exchange is NOT byte-exact for depth-profile (TPV31) and
recommend re-baseline, but Phase 4's acceptance claims TPV31 byte-exact ("peer == local").
With R-001 fixed, the precise contract should be:
- Phase 1 (bulk shared faces, centroid exchange): CORRECTS the neighbour material — TPV31
  bulk shared-face flux shifts by O(grad·h) (re-baseline; more correct). This is fine.
- Phase 4 (fault per-side, fault-QP exchange per R-001): the TPV31 FAULT stays SYMMETRIC
  (`Zp_plus == Zp_minus`) — this MUST be byte-exact / unchanged (the fault is where the σ_n
  leak lives; any spurious asymmetry is a regression, not a "re-baseline").
The plan must separate these: bulk = correction (re-baseline OK); fault per-side = strict
no-spurious-contrast (the R-001 fix guarantees it).

**Suggested fix:** in the plan, split the byte-exact constraint into (a) bulk shared faces
(correction, re-baseline) and (b) shared FAULT faces (MUST stay symmetric for depth-profile;
verified by the R-001 test).  Add an explicit acceptance: "np>1 TPV31 WITHOUT fault-locality:
every shared fault DOF has Zp_plus == Zp_minus to 1e-12."

**Test case:** the R-001 depth-profile np=2 fault test is the acceptance for (b).

---

### [R-004] MODERATE — Plan omits a TPV31 (symmetric) verification entirely; all new tests use a CONTRAST

**Category:** EDGE_CASE (test coverage)

**Description:**
The plan's tests (np=2 seam test, np=8 reproduction) all put a CONTRAST on the seam. None
verify the SYMMETRIC case (depth-profile / same-material-both-sides) stays symmetric on a
shared face — which is exactly where R-001 hides. The user explicitly asked to "verify against
TPV31 AND TPV6."

**Suggested fix:** add to the test plan a symmetric/depth-profile np=2 case for BOTH a shared
INTERIOR face (central + upwind build must reproduce the homogeneous flux to 1e-12) AND a shared
FAULT face (Zp_plus == Zp_minus). Add a Constant-material np=2 byte-exact case for the shared
fault per-side (Phase 4) too (the plan only has the Constant case for Phase 1).

**Test case:** see R-001 (fault) + an analogous shared-interior-face symmetric central-build
parity check.

---

### [R-005] LOW — Phase 3: shared-reclassify diagnostic count double-counts across the seam (plan acknowledges but leaves it fuzzy)

**Category:** QUALITY

**Description:**
A reclassified shared corridor face is seen by both ranks (both reclassify it — the predicate
is symmetric, so this is CONSISTENT and correct for the physics). The plan's printed
"reclassified" total adds `n_reclass_shared` per rank, double-counting. The plan says "document
the caveat" — acceptable, but a clean dedup is cheap.

**Suggested fix:** count a shared reclassified face only on the rank that owns the lower GLOBAL
element id of the pair (or simply report local + "shared (per-rank)" separately), so the printed
global total is exact. Diagnostic only — does not affect the (correct, consistent)
reclassification.

---

### [R-006] LOW [POSSIBLE] — Phase 1: the L2-order-0 vdim-3 byNODES FaceNbrData component order is assumed; assert it

**Category:** ASSUMPTION

**Description:**
Phase 1 reads `nbr[nbr_vdofs[0..2]]` as `{lambda, mu, rho}`. For L2-0 (1 dof/elem) vdim-3
byNODES this is correct, but it is the exact layout assumption that caused the R-004 bulk-Q
scramble. Also: `mat_gf.ExchangeFaceNbrData()` must trigger `pmesh.ExchangeFaceNbrData()` —
verify MFEM does this on first call for a fresh space (it does, but the plan should not assume
silently).

**Suggested fix:** keep the `MFEM_VERIFY(nbr_vdofs.Size() == 3, ...)`; have the np=2 test assert
the EXACT peer values (lambda, mu, rho separately, not just "different from local") so a
component transposition is caught.

---

### [R-007] LOW — Plan understates consumer B's reach: the shared UPWIND build feeds BOTH ADER and RK bulk shared-face flux

**Category:** QUALITY (clarity)

**Description:**
The plan frames consumer B (`BuildPerFaceBimaterialFluxMatrices_`, inl:442) as the "Arm 1
upwind" path. It actually populates `per_face_bimaterial_flux_`, which the shared-face bulk flux
uses on BOTH the ADER (`ComputeADERSharedFaceFluxRHS`) and RK (`SharedInteriorFaceFlux_`) paths,
and for the UPWIND portion of mixed flux too. So Phase 1 fixing the exchange corrects all bulk
shared-face flux regardless of integrator — good, but the plan should say so (and the
verification should include an ADER np=8 dry-run, not only RK).

**Suggested fix:** note in Phase 1 that the exchange corrects `per_face_bimaterial_flux_` for all
integrators; add an Arm-1 (ADER) np=8 dry-run to the acceptance alongside the RK one.

---

## Summary
- Critical issues: 1 (R-001 — fault per-side uses peer centroid; breaks TPV31 symmetry + CVM correctness)
- Moderate issues: 3 (R-002 index-space conflation; R-003 TPV31 contract split; R-004 missing symmetric tests)
- Low issues: 3 (R-005 diagnostic double-count; R-006 layout assertion; R-007 ADER coverage wording)
- Plan compliance: N/A (plan, not code) — the plan is ~80% sound; the per-element exchange
  (Phase 1) is correct for BULK shared faces (B/D), but Phase 4 (fault per-side) needs a
  SEPARATE fault-QP-resolved exchange — the centroid is the wrong quantity there.
- Verdict: **PASS WITH FIXES (plan must be revised before implementation)** — fold R-001 in as a
  distinct phase (fault-QP peer-material exchange) and add the symmetric/TPV31 tests (R-004)
  before coding.

## Codebase generalization audit (other parts checked — for the "generalize ALL parts" goal)
Verified these do NOT need changes (already general or not material-neighbour-dependent), so the
plan's 5-consumer inventory is complete EXCEPT the Phase-4 fault-QP gap:
- `ComputeMaxDt` (inl:946): per-element local material only — CFL is per-element wave speed, no
  neighbour needed. OK.
- Scalar `WaveOperator` path: homogeneous single material — no neighbour material. OK.
- Fault-flux conversion at shared faces (`FluxForElem_(e1)`, wave_operator.inl:3715): each rank
  applies its OWN local element's A to its own side — no peer A needed for the conversion (the
  peer applies its side on its rank). Only the per-side TRACTION (friction solve) needs the peer
  IMPEDANCE → that is exactly Phase 4 (covered, modulo R-001). OK.
- Free-surface / PML / absorbing fluxes: single-sided. OK.
- `GodunovFluxPool` dedup (round_sig 6): `flux_nbr` is built full-precision directly from the
  exchanged material (NOT via the deduped pool) — same as the stub did → byte-exact vs stub; a
  pre-existing local(deduped)-vs-nbr(full) asymmetry the exchange does not change. OK (note only).

## Unreviewed Areas
- The actual MFEM `ParGridFunction(L2-0, vdim=3).ExchangeFaceNbrData()` + `GetFaceNbrElementVDofs`
  round-trip — not run (plan, not code); the np=2 test (R-006) is the gate.
- Whether `shared_fault_elem1_on_plus_` / `shared_fault_dof_offset_` are populated before
  `AssignFaultSidePerMaterialImpedances` runs at np>1 (Phase 4 R6 in the plan) — confirmed the
  members exist; population ordering must be verified at implementation.
