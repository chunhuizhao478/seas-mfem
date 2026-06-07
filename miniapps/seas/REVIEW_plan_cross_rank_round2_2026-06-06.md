# Code Review (round 2): PLAN_cross_rank_material_exchange — TEST-GATING audit

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md`
- Round-1 review: `REVIEW_plan_cross_rank_material_exchange_2026-06-06.md` (R-001…R-007).
- Lens (this round): the user requirement — **EVERY change must be gated by a unit test**
  (a test that FAILS before the change and PASSES after). I audited each planned code change
  against the plan's Testing Strategy and flag every change whose only "gate" is an np=8
  `--dry-run` integration smoke (NOT a unit test, needs the 200 m mesh + 8 ranks, not
  CI-runnable) or that has no gate at all.
- Domain: CLAUDE.md ("verify via compile + unit tests"), the memory `feedback-no-local-mesh-runs`
  (no full-mesh runs — so the np=8 200 m dry-run is NOT a valid local/CI gate).

## Change-by-change gating matrix (the core audit)
| Plan change | Current gate in plan | Unit-gated? |
|-------------|----------------------|-------------|
| P1 exchange returns PEER material | np=2 neighbour-material test | YES |
| P1 Constant byte-exact | np=2 Constant test | YES |
| P1 L2-0 vdim-3 component order | (implied) | WEAK — see R-205 |
| P1 R-004 warning removed | none | NO — R-206 |
| P2 remove seam_continuous abort | "np=8 dry-run constructs" | NO (integration) — R-203 |
| P2 shared central build CORRECT values | np=2 item 4 "compare to serial" | PARTIAL — R-203 |
| P3 shared contrast guard reclassifies | np=2 strong-contrast test | YES |
| P3 IMPL-8 size invariant after SHARED reclass | "verify the erase loop" (prose) | NO — R-202 |
| P3 n_reclass_shared diagnostic | none | NO (LOW) — R-207 |
| P4 shared fault per-side CONTRAST | np=2 bimaterial-fault test | YES |
| **P4 shared fault per-side SYMMETRIC (TPV31)** | **none** | **NO — R-201 (CRITICAL)** |
| P4 sf-vs-position index correctness | (implied by the contrast test) | NO — R-204 |
| P5 seam_continuous now inert | "np=8 dry-run" | NO — R-208 |
| P5 jobs: fault-locality optional | np=8 dry-run | n/a (job files, not unit-testable) |

## Findings

### [R-201] CRITICAL — No unit test gates the depth-profile/symmetric SHARED fault face (the R-001 regression is undetectable)

**Category:** EDGE_CASE (missing gating test for a correctness-critical path)

**Description:**
Round-1 R-001 showed Phase 4's centroid-for-peer makes a symmetric (depth-profile / TPV31)
shared fault face spuriously bi-material (`Zp_plus != Zp_minus`). The plan has NO unit test that
would catch this — all Phase-4 tests use a CONTRAST. Per "every change gated by a unit test,"
the fault-QP-exchange fix (R-001) MUST ship with a test that fails under centroid-for-peer and
passes under the fault-QP exchange.

**Trigger:** implement Phase 4 (any variant) with no symmetric-fault gate → TPV31 regresses
silently.

**Expected behavior:** a depth-profile / same-material-both-sides fault on a partition seam
yields `Zp_plus == Zp_minus` (to 1e-12) at every shared fault DOF.

**Suggested fix (add to the plan's Testing Strategy as a REQUIRED gate for Phase 4):**
```cpp
// tests/parallel/test_bimaterial_seam_material_np2.cpp  (run np=2)
void test_R201_depthprofile_shared_fault_stays_symmetric() {
   // 2 hexes, 1/rank, inter-rank face tagged FAULT.  Mode::Coefficient material = f(z)
   // ONLY (depth profile), identical on both sides of the fault.
   // After AssignFaultSidePerMaterialImpedances:
   //   for the shared fault DOF: |Zp_plus - Zp_minus| <= 1e-12 * Zp_plus   (symmetric)
   // FAILS with peer=centroid (different depths); PASSES with the fault-QP eps-offset exchange.
}
```

**Test case:** the snippet above IS the gating test (and it is the acceptance for round-1 R-001).

---

### [R-202] MODERATE — No unit test gates the IMPL-8 size invariant after a SHARED reclassification

**Category:** EDGE_CASE

**Description:**
Phase 3 newly reclassifies SHARED corridor faces out of `central_flux_face_set_`. The IMPL-8
invariant (`per_face_central_flux_.size() == central_flux_face_set_.size()`, asserted in
`SetMixedFluxMode`) must still hold when a SHARED face is reclassified — the erase loop must
remove it from the set AND it must never have been inserted into `per_face_central_flux_`. The
plan only says "verify the erase loop handles shared mesh_face indices" in prose; no test gates it.

**Trigger:** a shared corridor face reclassified (strong contrast on a seam) — if the erase loop
misses it (or it was inserted into per_face_central_flux_), the SetMixedFluxMode MFEM_VERIFY
fires (or, worse, a silent size mismatch).

**Suggested fix:** in the np=2 strong-contrast test, after construction assert the invariant
directly:
```cpp
// (np=2, strong contrast on a shared corridor face)
TEST_ASSERT(op.GetPerFaceCentralFlux().size() == op.GetCentralFluxFaceSet().size(),
            "IMPL-8: per_face_central_flux_ size == central_flux_face_set_ size after "
            "SHARED reclassification");  // (add a const accessor for per_face_central_flux_ if absent)
```

**Test case:** the assertion above, added to the np=2 strong-contrast test.

---

### [R-203] MODERATE — Phase 2 "construction does not abort" is an np=8 integration smoke, not a unit gate; the shared CENTRAL build needs a numerical unit gate

**Category:** EDGE_CASE / DEVIATION (test rigor)

**Description:**
Phase 2's headline acceptance is "np=8 dry-run on the 200 m mesh constructs without abort." That
is (a) an integration smoke, not a unit test, (b) a full-mesh run the project forbids locally
(`feedback-no-local-mesh-runs`), and (c) only checks "no abort," not that the shared central
matrices are CORRECT. The plan's np=2 item 4 ("weak contrast → compare to serial") is the right
idea but is not called out as the REQUIRED Phase-2 gate.

**Trigger:** a wrong per-side A in the shared central build (e.g. local/peer swapped) would pass
"no abort" but be physically wrong.

**Suggested fix:** make the Phase-2 gate a NUMERICAL np=2 unit test: build the same 2-element
contrast face in SERIAL (both elements local) and at np=2 (the face shared); assert the central
flux matrices (or the resulting `InteriorFaceFlux_`/`SharedInteriorFaceFlux_` deposit for a fixed
(Q_self, Q_nbr)) are EQUAL to 1e-12. Demote the np=8 dry-run to a non-gating "production smoke."
```cpp
void test_R203_shared_central_matches_serial_np2() {
   // serial 2-elem central build  vs  np=2 shared central build, same materials + states
   // assert worst |F_serial - F_np2| <= 1e-12 * scale
}
```

**Test case:** the snippet above.

---

### [R-204] MODERATE — No unit test gates the Phase-4 sf-vs-position index mapping (round-1 R-002)

**Category:** ASSUMPTION

**Description:**
`shared_fault_elem1_on_plus_` is position-indexed (`si`), `shared_fault_dof_offset_` is raw-sf
keyed. The contrast test would catch a gross swap only if it happens to land on a different DOF;
a partition with one shared fault face would NOT distinguish `si` from `sf` (they coincide when
there is a single shared fault face). The gate must use a fixture with >=2 shared fault faces so
`si != sf`, and assert each shared fault DOF gets the value for ITS face (not a transposed one).

**Suggested fix:** in the fault test, use >=2 shared fault faces with DISTINCT per-side contrasts
and assert each shared fault DOF carries its own face's `Zp_plus/Zp_minus` (a transposition test).

**Test case:**
```cpp
void test_R204_shared_fault_index_mapping_np2() {
   // >=2 shared fault faces, each with a DIFFERENT contrast magnitude;
   // assert dof_data[shared_fault_dof_offset_[sf_a]+q].Zp_plus matches face a's value,
   // and likewise for sf_b — catches an si<->sf transposition.
}
```

---

### [R-205] MODERATE — Phase 1 component-order (lambda,mu,rho) not asserted per-component (R-204→R-006 round-1)

**Category:** ASSUMPTION

**Description:**
The np=2 neighbour-material test must assert the peer's lambda, mu, rho SEPARATELY (not just
"!= local" or a scalar Zp), or a component transposition in the L2-0 vdim-3 FaceNbrData unpack
passes silently — the exact R-004 failure mode.

**Suggested fix:** in the Phase-1 np=2 test, use a peer material with THREE distinct, non-
proportional values (e.g. lambda=1e10, mu=2e10, rho=2500) and assert each component of
`GetSharedFaceNeighbourMaterial()[sf]` equals the peer's, to 1e-12.

**Test case:** included in the Phase-1 np=2 test (per-component asserts).

---

### [R-206] LOW — Phase 1 removal of the R-004 warning is not gated

**Category:** QUALITY

**Description:**
The plan removes the one-time R-004 WARNING (the stub is gone). No test confirms it is no longer
printed for a Coefficient material in parallel (a stray warning would mislead users into thinking
the stub is still active).

**Suggested fix:** LOW — optionally capture stdout in the np=2 Coefficient test and assert the
"LOCAL-SIDE stub" warning string is absent. Acceptable to skip (the warning's removal is mechanical).

---

### [R-207] LOW — Phase 3 `n_reclass_shared` diagnostic count is not gated

**Category:** QUALITY

**Description:**
The new shared-reclassify diagnostic count has no test. It is print-only (no physics impact), but
the dedup/double-count logic (round-1 R-005) could mislead.

**Suggested fix:** LOW — optionally assert the per-rank `n_reclass_shared` in the np=2
strong-contrast test (each rank reclassifies the 1 shared face → 1). Acceptable to skip.

---

### [R-208] MODERATE — Phase 5 (seam_continuous now inert) has no unit gate

**Category:** EDGE_CASE

**Description:**
Phase 5 makes `seam_continuous` a no-op. Two behaviors need gating: (1) a bimaterial config with
`seam_continuous=false` + shared faces now CONSTRUCTS (no abort — the whole point); (2) a config
that still SETS `seam_continuous=true` still parses (no abort, backward-compat). The plan gates
neither with a unit test (only the np=8 dry-run).

**Suggested fix:** (1) is covered by the np=2 Phase-2 test if that test's config leaves
seam_continuous=false (assert no abort) — make that explicit. (2) add a config-parse assertion
(`seam_continuous=true` still parses to `MaterialSpec::seam_continuous==true`, the driver does not
abort) — extend `seas_test_tpv_config_parse`.
```cpp
// config-parse: a bimaterial-fault config with seam_continuous absent (=false) parses and
// (in the np=2 operator test) constructs with mixed flux + shared faces without abort.
```

---

### [R-209] MODERATE — The fault-QP peer-material exchange (R-001 fix) adds a COLLECTIVE inside AssignFaultSidePerMaterialImpedances; the plan must gate collective-safety

**Category:** BUG-risk (MPI)

**Description:**
The round-1 R-001 fix (exchange the peer's fault-QP eps-offset material) introduces an MPI
exchange INSIDE `AssignFaultSidePerMaterialImpedances` (called post-ctor by the driver). It must
be collective — EVERY rank must call it even with zero shared fault faces, or it deadlocks
(the same class as the contrast-guard MPI_Reduce). The plan currently puts the exchange at the
ctor (Phase 1) but the fault-QP variant must run after the fault setup (AssignFault) — a second
collective at a different point. No test gates the np>1 no-shared-fault-on-this-rank case.

**Trigger:** an MPI rank that owns fault DOFs but no SHARED fault face (common at np>1) — if the
exchange is guarded by `if (have shared fault faces)`, that rank skips the collective → deadlock.

**Suggested fix:** the fault-QP exchange must be UNCONDITIONAL on the parallel path (run on all
ranks). Gate with an np=2 test where ONE rank has a shared fault face and the partition is such
that the other rank's participation is required (the test simply completing without hang is the
gate; assert it returns).

**Test case:** the np=2 fault tests (R-201/R-204) already exercise both ranks; add an explicit
"all ranks reach AssignFault return" (no-hang) expectation, and an np=2 case where rank 1 has no
local shared fault face but must still participate in the collective.

---

## Summary
- Critical issues: 1 (R-201 — no gate for the symmetric/TPV31 shared fault face; the R-001 regression is invisible)
- Moderate issues: 5 (R-202 IMPL-8 invariant; R-203 numerical central gate vs np8 smoke; R-204 index-mapping; R-205 per-component; R-208 seam_continuous-inert; R-209 collective-safety)
- Low issues: 2 (R-206 warning removal; R-207 diagnostic count)
- Plan compliance: PARTIAL re: the "every change unit-gated" requirement — ~half the changes are
  gated only by the np=8 200 m `--dry-run`, which is an integration smoke and a FORBIDDEN local
  full-mesh run, not a unit test.
- Verdict: **PASS WITH FIXES** — before implementation the plan's Testing Strategy must (1) add the
  symmetric-fault gate (R-201), the IMPL-8 gate (R-202), the numerical central gate (R-203), the
  index-mapping gate (R-204), the per-component gate (R-205), the seam_continuous-inert gate
  (R-208), the collective-no-hang gate (R-209); and (2) reclassify EVERY np=8 200 m `--dry-run`
  acceptance from a GATE to a non-gating "production smoke," with a self-contained np<=2 unit test
  as the actual gate for each phase.

## Cross-cutting requirement (the user directive, made concrete)
Add a one-line "Gate:" field to EVERY phase's Acceptance Criteria naming the SPECIFIC unit test
(serial or np<=2) that fails-before/passes-after. No phase may list an np=8 full-mesh dry-run as
its gate (those are forbidden locally + are smokes, not gates). Consolidate the new parallel
tests into `tests/parallel/test_bimaterial_seam_material_np2.cpp` with the cases:
neighbour-material(per-component) · Constant-byte-exact · central-matches-serial ·
strong-contrast-reclassify(+IMPL-8) · fault-contrast · **fault-symmetric** · fault-index-mapping ·
collective-no-hang; plus a serial central-parity case and the config-parse seam_continuous case.

## Unreviewed Areas
- Whether MFEM `ParGridFunction(L2-0,vdim3).ExchangeFaceNbrData()` round-trips correctly — only
  the np=2 test (R-205) will confirm at implementation.
- Performance of a second collective in AssignFault (R-209) at production np — not gated by a unit
  test (acceptable; correctness first).
