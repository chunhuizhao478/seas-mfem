# Cross-rank Material Exchange — Phase 4 Adversarial Code Review

Plan:  document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md (Phase 4)
Scope: dynamic/bimaterial_wave_operator.inl (AssignFaultSidePerMaterialImpedances
       + MaterialAtLocal_/MaterialAtNbr_), dynamic/bimaterial_wave_operator.hpp,
       tests/parallel/test_bimaterial_seam_fault_np2.cpp, Makefile.
Reviewer: chunhui-code-reviewer (independent re-verification)
Date: 2026-06-07

Re-verified independently (built + ran in the worktree, np=2):
- seas_test_bimaterial_seam_fault_np2 .......... 24/24 (all 5 gates)
- seas_test_bimaterial_fault_perside_material .. 11/11 (interior reroute byte-exact)
- seas_test_bimaterial_seam_material_np2 ....... 135/135 (P1-P3 regressions)

Mandate: assume >= 3 bugs, find them; focus CRITICAL/correctness — wrong per-side
impedance, plus/minus transposition, R-204 index spaces, eps-offset peer Jacobian
(R-401), interior byte-exactness, collective hazards, gate vacuity.

================================================================================
VERDICT UP FRONT
================================================================================
The Phase-4 IMPLEMENTATION is, as far as I can verify, CORRECT. I found no live
critical bug in the production code path:
- The si/sf index spaces are correct and match the base canonical loop exactly.
- The plus/minus assignment is consistent with the interior loop AND with the
  serial interior-fault oracle (verified by an unsorted-value probe: +y/strong ->
  Zp_plus on BOTH serial and parallel).
- The eps-offset peer Jacobian (R-401) is mathematically sound; depth-profile
  symmetry holds to 0.0 on the skewed peer.
- The interior reroute is genuinely byte-exact (perside 11/11; the only new op is
  a benign T.SetIntPoint(&ip) before an explicit-ip EvalAt).
- No MPI collective inside AssignFault; no rank-conditional collective hazard.

The real defects are in the TEST SUITE: the contrast / gridfunction / index_mapping
gates are plus/minus-INSENSITIVE (sorted {min,max}), so a future plus/minus
transposition regression would pass silently. Two gates are also structurally
vacuity-prone. These are MODERATE (they let a whole bug-class through the gate),
not critical, because the current code is correct.

Critical issues: 0

================================================================================
FINDINGS
================================================================================

--------------------------------------------------------------------------------
P4-001  [MODERATE] [test-coverage / silent-regression]
        Plus/minus TRANSPOSITION is invisible to every value-checking gate
--------------------------------------------------------------------------------
File: tests/parallel/test_bimaterial_seam_fault_np2.cpp
      shared_fault_contrast (lines 353-379),
      shared_fault_gridfunction (lines 415-432),
      shared_fault_index_mapping (lines 470-499).

Description:
Every assertion that inspects per-side impedance does so through a SORTED pair:

      const real_t lo = std::min(dof[i].Zp_plus, dof[i].Zp_minus);
      const real_t hi = std::max(dof[i].Zp_plus, dof[i].Zp_minus);

and asserts {lo,hi} == {weakZp, strongZp}. The serial oracle in shared_fault_
contrast is likewise sorted (s_lo/s_hi, lines 344-345). Consequently if the
shared-loop e1_plus swap (inl:1427-1430) were inverted — i.e. the weak (peer)
material written to Zp_plus and the strong to Zp_minus — EVERY one of these gates
would still pass. The sorted pair has zero power to detect a plus/minus swap.

This is the single most likely place for a silent bug to land in a future edit
(e.g. a refactor of shared_fault_elem1_on_plus_, the canonical-frame sign logic,
or the swap branch itself). The plan explicitly calls this out as R-204/R4 risk;
the gate does not actually cover the +/- ASSIGNMENT, only the per-side SET.

Verification that the CURRENT code is correct (so this is a gap, not a live bug):
I patched a probe to print the UNSORTED values in shared_fault_contrast. With the
halfspace `x(1) >= 0 ? 120e9 : 30e9` (strong on +y):
      [PROBE2 serial]      Zp_plus=3.10032e+07  Zp_minus=1.55016e+07
      [PROBE2 rank 0 par]  Zp_plus=3.10032e+07  Zp_minus=1.55016e+07
      [PROBE2 rank 1 par]  Zp_plus=3.10032e+07  Zp_minus=1.55016e+07
i.e. +y/strong -> Zp_plus on BOTH the serial interior oracle and the parallel
shared face. The convention is consistent and correct today. But the gate would
not have noticed had it been wrong.

Trigger: a future change that inverts shared_fault_elem1_on_plus_ for shared
faces, or transposes the e1_plus branch in the shared loop, would mis-orient the
fault traction/Riemann (FaultFaceFlux consumes Zp_plus vs Zp_minus by physical
side) and the gate would stay green.

Suggested fix (add a plus/minus-SENSITIVE assertion; do NOT sort):
In shared_fault_contrast, after computing the parallel dof, assert the physical
side directly using the same +y-strong knowledge as the serial oracle, and tie it
to the serial oracle's UNSORTED pair:

```diff
@@ shared_fault_contrast, replace the sorted match block @@
-      real_t worst_contrast = 0.0, worst_match = 0.0;
-      int n_checked = 0;
-      for (int i = nint * nqp; i < static_cast<int>(dof.size()); ++i)
-      {
-         const real_t lo = std::min(dof[i].Zp_plus, dof[i].Zp_minus);
-         const real_t hi = std::max(dof[i].Zp_plus, dof[i].Zp_minus);
-         worst_contrast = std::max(worst_contrast, RelErr(lo, hi) > 0 ? 0.0 : 1.0);
-         worst_match = std::max(worst_match,
-                                std::max(RelErr(lo, weakZp), RelErr(hi, strongZp)));
-         ++n_checked;
-      }
+      // PLUS/MINUS-SENSITIVE (P4-001): the serial interior oracle and the
+      // parallel shared face MUST place the SAME physical material on the SAME
+      // signed side.  sdof[0] is the unsorted oracle pair; assert the parallel
+      // pair equals it component-for-component (Zp_plus==oracle Zp_plus, not just
+      // the sorted set).  A plus/minus transposition now FAILS.
+      const real_t oraclePlus  = sdof.empty() ? 0.0 : sdof[0].Zp_plus;
+      const real_t oracleMinus = sdof.empty() ? 0.0 : sdof[0].Zp_minus;
+      real_t worst_match = 0.0; int n_checked = 0; bool contrast_ok = true;
+      for (int i = nint * nqp; i < static_cast<int>(dof.size()); ++i)
+      {
+         worst_match = std::max(worst_match,
+            std::max(RelErr(dof[i].Zp_plus,  oraclePlus),
+                     RelErr(dof[i].Zp_minus, oracleMinus)));
+         if (RelErr(dof[i].Zp_plus, dof[i].Zp_minus) <= 1e-6) { contrast_ok = false; }
+         ++n_checked;
+      }
```
and assert `worst_match <= 1e-12 && contrast_ok` under `if (nshr > 0)`. Apply the
same unsorted pattern (compare Zp_plus to the +y-strong value, Zp_minus to the
-y-weak value) in shared_fault_gridfunction and shared_fault_index_mapping.

Priority: HIGH (closes the bug-class the whole phase is most exposed to).

Test impact: with the fix above the gate STILL passes on the current (correct)
code — I verified the unsorted values match the oracle — and would FAIL on a
transposition.

--------------------------------------------------------------------------------
P4-002  [MODERATE] [test-vacuity]
        Constant / depthprofile gates pass VACUOUSLY on a rank with nshr==0
--------------------------------------------------------------------------------
File: tests/parallel/test_bimaterial_seam_fault_np2.cpp
      shared_fault_constant_byteexact (lines 273-280),
      shared_fault_symmetric_depthprofile (lines 300-307).

Description:
These two gates fold their value check over `for (i = nint*nqp; i < dof.size(); ++i)`
with NO `if (nshr > 0)` guard and initialise `worst = 0.0`. On a rank where
nshr == 0 (so dof.size() == nint*nqp), the loop body never runs, worst stays 0.0,
and `TEST_ASSERT(worst <= ...)` passes WITHOUT having examined any shared-fault
DOF. The only cross-rank guard, `glob_shr >= 1`, proves SOME rank owns a shared
face — not that the rank executing the value assertion does.

For the current 2-tet PartitionByYSign meshes BOTH ranks happen to have nshr==1
(I probed: `constant nint=0 nshr=1 nqp=3 dofsize=3` on both ranks), so the gates
are NOT vacuous today. But the structure is fragile: any future partition / mesh
where the shared fault is owned asymmetrically (one rank nshr==0) would let that
rank's assertion pass while contributing 0 real checks — and the failing rank's
real check could be masked if the comparison is loose.

Trigger: a partition or mesh edit that leaves a participating rank with nshr==0,
or a 3+ rank generalisation.

Suggested fix: gate the value loop on the LOCAL nshr and require non-empty work,
mirroring the other three gates:
```diff
-      real_t worst = 0.0;
-      for (int i = nint * nqp; i < static_cast<int>(dof.size()); ++i)
-      { worst = std::max(worst, RelErr(dof[i].Zp_plus, Zp_const)); ... }
-      TEST_ASSERT(worst <= 1e-12, ...);
+      if (nshr > 0)
+      {
+         int n_checked = 0; real_t worst = 0.0;
+         for (int i = nint * nqp; i < static_cast<int>(dof.size()); ++i)
+         { worst = std::max(worst, RelErr(dof[i].Zp_plus, Zp_const));
+           worst = std::max(worst, RelErr(dof[i].Zp_minus, Zp_const)); ++n_checked; }
+         TEST_ASSERT(n_checked > 0 && worst <= 1e-12, ...);
+      }
```
(The `n_checked > 0` requirement is what makes a future vacuous pass FAIL.)

Priority: MEDIUM.

--------------------------------------------------------------------------------
P4-003  [LOW] [dead-code / readability]
        `v3y` skew ternary is a no-op (1.0 ? 1.0); the skew is real but mislabeled
--------------------------------------------------------------------------------
File: tests/parallel/test_bimaterial_seam_fault_np2.cpp, line 79:
      const real_t v3y = skew ? 1.0 : 1.0;

Description:
`v3y` is assigned 1.0 regardless of `skew`. The ACTUAL skew that makes the peer
asymmetric is applied via the +y apex x/z displacement on line 84
(`{skew ? 0.4 : 0.0, v3y, skew ? 0.3 : 0.0}`), which DOES skew the peer tet's
shape and exercise the peer Jacobian (R-401). So the test is still meaningful —
the depthprofile gate passed to 0.0 on a genuinely skewed peer. But the `v3y`
ternary is dead and misleading: it reads as if the y-extent changes with skew
when it does not. Harmless to correctness; it is a maintenance/readability trap.

Note (why y MUST stay at the same magnitude for this test): the fault face v0v1v2
lies in y=0 and its physical normal is purely +/-y. The eps-offset rides along
that normal, so z is preserved on both sides and f(z) is symmetric REGARDLESS of
the x/z apex skew. That is exactly what R-401 needs to demonstrate. Changing v3y
asymmetrically would NOT add coverage; the x/z skew is the right lever.

Suggested fix: drop the dead variable and inline the constant, with a comment that
the y-extent is intentionally symmetric (only the in-plane apex is skewed):
```diff
-   const real_t v3y = skew ? 1.0 : 1.0;
    real_t verts[5][3] = {
       {0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {0.0, 0.0, 1.0},
-      {skew ? 0.4 : 0.0, v3y, skew ? 0.3 : 0.0},        // v3 apex (+y)
+      // +y apex: skewed in-plane (x,z) so the peer tet has a DIFFERENT shape
+      // (peer Jacobian != local — R-401); y-extent kept = 1 on purpose so the
+      // fault face stays in y=0 (normal purely +/-y, eps-offset preserves z).
+      {skew ? 0.4 : 0.0, 1.0, skew ? 0.3 : 0.0},
```
Priority: LOW.

--------------------------------------------------------------------------------
P4-004  [LOW] [doc / plan-vs-code drift]
        Plan/comment say `.at(sf)` (throws); code correctly uses `.find()+continue`
--------------------------------------------------------------------------------
File: PLAN ...md Phase 4 ("base = shared_fault_dof_offset_.at(sf)") and the
      inl:1354-1356 comment ("shared_fault_dof_offset_ is keyed by the RAW
      shared-face index sf"); implementation inl:1370-1372.

Description:
The plan text and the prompt both reference `.at(sf)`, which would THROW
std::out_of_range if a shared fault face lacked an offset entry. The actual code
is SAFER than the spec — it uses `sh_dof_off.find(sf)` with a `continue` on miss
(inl:1370-1372), matching the interior loop's `.find()+continue` (inl:1278-1279).
This is the better choice (no hard throw on a benign missing entry; both ranks
behave locally). No code change needed; this is a NOTE that the implementation
diverged from the plan wording in the safe direction. Worth a one-line update to
the plan/comment so the audit trail reflects `.find()`, not `.at()`.

NOTE on the base canonical loop divergence: EvaluateBulkAtFaultQPsCanonical uses a
HARD MFEM_VERIFY on the missing shared_fault_dof_offset_ entry (wave_operator.inl
:2334-2338), whereas Phase-4 AssignFault silently `continue`s. This asymmetry is
defensible (AssignFault is a pre-step setup that can leave a DOF at its seeded
value; the canonical loop is the hot RHS where a missing entry IS a corruption),
but it means a genuinely-corrupt shared_fault_dof_offset_ would be caught by the
RHS loop, not by AssignFault. Acceptable; documented here for completeness.
Priority: LOW.

================================================================================
ITEMS ATTACKED AND CLEARED (no defect)
================================================================================

R-204 index spaces (si POSITION vs sf RAW) — CORRECT.
  Phase-4 shared loop: e1_plus = shared_fault_elem1_on_plus_[si] (POSITION si in
  fault_shared_faces_); dof_off via shared_fault_dof_offset_.find(sf) (RAW sf).
  Base EvaluateBulkAtFaultQPsCanonical: shared_fault_elem1_on_plus_[sf_idx]
  (POSITION) + shared_fault_dof_offset_.find(sf) (RAW). IDENTICAL mapping.
  SetFaultDOFData (wave_operator.hpp:744-749) builds
  shared_fault_dof_offset_[fault_shared_faces_[i]] = nint*nqp + i*nqp, so the
  test's (nint+si)*nqp block start == the operator's offset. Verified the
  index_mapping gate runs non-vacuously: n_blocks == nshr asserted, worst 0.0.

Plus/minus CONSISTENCY of the precomputed flag — CORRECT.
  shared_fault_elem1_on_plus_ (wave_operator.inl:464-521) uses the SAME canonical
  convention as interior_fault_elem1_on_plus_ ("smaller projection onto the
  canonical face normal => + side"), via NormalNeedsFlipToCanonical so both ranks
  agree on the canonical normal line. The shared loop applies the identical swap
  branch as the interior loop. Probe confirmed +y/strong -> Zp_plus on serial AND
  parallel (P4-001). (The GATE doesn't test this — that's P4-001 — but the CODE
  is right.)

eps-offset peer Jacobian (R-401) — SOUND.
  offset_ip normalises nr = J^{-1} n_phys to unit REFERENCE length, so the
  physical offset = eps_ref * n_phys / |J^{-1} n_phys| — i.e. PURELY along the
  exact face normal n_phys, with an element-dependent magnitude. For a fault in
  y=0 (normal +/-y) the z-coordinate is preserved on both sides => depth profile
  f(z) identical => exact symmetry on ANY mesh, including the skewed peer.
  Confirmed: shared_fault_symmetric_depthprofile worst rel = 0.0 on skew=true.
  The peer offset uses in2 = -in1 through *ftr->Elem2's own Jacobian and lands
  INSIDE the peer element (else the contrast gate would collapse to no-contrast;
  it shows real contrast). GetElement2IntPoint() after SetAllIntPoints is the
  correct peer reference base point (eip2 == Loc2.Transform(face_ip)).

Interior byte-exactness (EvalAt -> MaterialAtLocal_) — GENUINELY byte-exact.
  Old: material_->EvalAt(Elem1No, *ftr->Elem1, ip1) -> lambda_coef->Eval(T, ip1).
  New: MaterialAtLocal_ does T.SetIntPoint(&ip1) then the SAME EvalAt. The only
  added op is SetIntPoint; the production coefficients are explicit-ip
  FunctionCoefficients (heterogeneous_material: Eval reads the ip arg /
  T.Transform(ip), not T.GetIntPoint()), so the result is identical. Elem2 read
  is local (Elem2No < ne_ for an interior 2-sided face). Confirmed perside 11/11,
  seam_material 135/135.

Collective safety — CLEAN.
  AssignFaultSidePerMaterialImpedances performs NO MPI collective; it iterates
  fault_shared_faces_ purely locally. MaterialAtNbr_ reads already-exchanged
  FaceNbrData (GridFunction) or evaluates the coefficient locally (Coefficient/
  Constant) — no communication. The only collectives (ExchangeFaceNbrData + the
  3 GF exchanges) live in the ctor under if constexpr(IsParallelMesh) and run
  unconditionally on ALL ranks (R-5). A rank with nshr==0 just runs an empty loop
  — no rank-conditional collective. No deadlock hazard.

offset_ip lambda scope + ftr reuse ordering — CORRECT.
  The lambda (inl:1251) captures eps_ref by [&]; both are function-local and in
  scope for the shared loop. Per QP: ftr->SetAllIntPoints(&ip) (populates eip1/2)
  -> read xf/nvec from ftr->Face -> capture eip1/eip2 by value into offset_ip ->
  offset_ip internally SetIntPoint on the element transform (does NOT mutate
  ftr->eip1/2) -> MaterialAtLocal_/MaterialAtNbr_ re-SetIntPoint to the offset ip
  before EvalAt. No stale-ip aliasing; the const& to eip1 is stable across the
  internal element-transform SetIntPoint.

index_mapping has 2 DISTINCT contrasts — CONFIRMED non-vacuous.
  modf: -y always 30e9; +y col=1 -> 60e9 (x<1 column), col=2 -> 90e9 (x>1
  column). The two y=0 fault faces carry different +y moduli; each DOF matches
  ITS face's pair to 0.0 with n_blocks == nshr asserted.

================================================================================
PROPOSED UNIT TESTS (additions; not written into the suite)
================================================================================

Test NEW-1  shared_fault_plusminus_sensitive  (HIGH)
  Target: AssignFaultSidePerMaterialImpedances (shared loop e1_plus swap).
  Validates: the +y/strong material lands in Zp_plus (not just the sorted set),
  matching the serial interior oracle UNSORTED, per physical side.
  Sketch: reuse shared_fault_contrast; instead of min/max, assert
  RelErr(dof[i].Zp_plus, sdof[0].Zp_plus) <= 1e-12 AND
  RelErr(dof[i].Zp_minus, sdof[0].Zp_minus) <= 1e-12 (oracle pair UNSORTED).
  Negative control: temporarily TamperSharedFaultElem1OnPlus(0) (the existing
  SEAS_TEST_INTERNAL mutator, wave_operator.hpp:827) and assert the test now
  FAILS — proving the gate has plus/minus discrimination power.

Test NEW-2  shared_fault_nonvacuous_guard  (MEDIUM)
  Target: the constant/depthprofile gates.
  Validates: each rank with nshr>0 actually examines >=1 shared-fault DOF.
  Sketch: add `n_checked` counters under `if (nshr > 0)` and assert
  n_checked == nshr*nqp (see P4-002 diff).

Test NEW-3  shared_fault_transposition_negative_control  (MEDIUM)
  Target: regression guard against P4-001.
  Validates: with shared_fault_elem1_on_plus_[k] flipped via the test mutator,
  the plus/minus-sensitive assertion FAILS for the affected face and PASSES for
  the others — locking the gate's discrimination power into CI.

================================================================================
PRIORITY SUMMARY
================================================================================
| ID     | Severity | Status | Description                                       |
|--------|----------|--------|---------------------------------------------------|
| P4-001 | MODERATE | OPEN   | +/- transposition invisible to sorted gates       |
| P4-002 | MODERATE | OPEN   | constant/depthprofile gates vacuous if nshr==0    |
| P4-003 | LOW      | OPEN   | dead `v3y = skew?1.0:1.0` ternary (misleading)    |
| P4-004 | LOW      | OPEN   | plan says .at(sf); code uses safer .find()+continue|

Critical issues: 0

================================================================================
SUMMARY
================================================================================
Phase 4 is implemented correctly. The si/sf index spaces match the base canonical
loop exactly; the eps-offset peer Jacobian is mathematically sound (depth-profile
symmetric to 0.0 on a genuinely skewed peer); the plus/minus assignment agrees
with both the interior loop and the serial interior-fault oracle by physical side;
the interior reroute is byte-exact (perside 11/11); there is no MPI collective in
AssignFault and no rank-conditional collective hazard. All 24/24 Phase-4 gates and
the 11/11 + 135/135 regressions reproduce.

The defects are entirely in the TEST coverage. P4-001 (MODERATE) is the headline:
all per-side gates compare a SORTED {min,max} pair, so a plus/minus transposition
— the most physically dangerous regression and the one the plan flags as R-204/R4
— would pass silently. I confirmed the current assignment is correct via an
unsorted probe, but the gate has no power to keep it that way. P4-002 (MODERATE)
makes two gates structurally vacuity-prone. P4-003/P4-004 are LOW cleanups. None
block merge; P4-001's plus/minus-sensitive assertion (+ the tamper negative
control) should be added before this is relied on as a regression gate for the
fault sign convention.
