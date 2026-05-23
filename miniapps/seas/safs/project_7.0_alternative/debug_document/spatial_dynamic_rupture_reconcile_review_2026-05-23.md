# Code Review: 2026-05-23 — shared-fault cross-rank reconcile fix DESIGN

## Review Scope
- Plan reviewed: `PLAN_dip_channel_v1_cross_rank_stability_2026-05-22.md`
  → section "REVISED FIX DESIGN (owner-broadcast / reconcile)".
- Confirmed root cause: `spatial_dynamic_rupture_speckle_blowup_2026-05-22.md`
  (Frontera job 7747036 XRANK trace).
- Source context: `dynamic/wave_operator.inl` `ComputeADERSharedFaceFluxRHS`
  (LSW branch ~4758-4797; reconcile site), R-701 note (:3026-3043),
  `shared_fault_elem1_on_plus_` (:519-525), R-101 verify (:5364+),
  `dynamic/fault_face_flux.cpp` `EvaluateADER_LSW` / `BuildImposedState`.
- Directive under review: **the fix must be method-invariant** (not gated to LSW).
- No fix code exists yet; this reviews the DESIGN against the directive and for
  correctness, per the /code-review contract.

## Findings

### [R-001] CRITICAL [DEVIATION] design — the LSW gating violates method-invariance

**Category:** DEVIATION

**Description.** The drafted fix gates the reconcile to
`fault_friction_law_ ∈ {LSW, LSW_ForcedRupture}` "so rate-state stays
byte-exact." But the defect is **method-invariant**: two ranks compute the
shared-fault DOFData redundantly from inputs that differ at ~1e-14 (the same
physical QP interpolated through two element parametrizations, `shape1∘Loc1` vs
`shape2∘Loc2`). Rate-state has the *same* cross-rank inconsistency; it merely
does not *blow up* because the smooth `asinh` solve lacks the LSW `max(0,·)`
kink. Gating to LSW leaves rate-state cross-rank-inconsistent at ~1e-14 and one
non-smoothness (a `min`/`max` strength cap, a healing event, a future law) away
from the same desync. A method-specific fix to a method-invariant defect is the
exact "problem-specific" failure mode the directive forbids.

**Trigger.** Any friction law on a shared fault face at np>1 over a long run:
the two ranks' DOFData diverges at ~1e-14; for rate-state it stays bounded
(passes the 1e-10 R-101 tol today) but is never bit-consistent.

**Actual (design).** Reconcile only for LSW; rate-state keeps redundant,
inconsistent per-rank computation.

**Expected.** Reconcile the shared-fault DOFData for **all** friction laws in
`ComputeADERSharedFaceFluxRHS` (and the rate-state `else` branch), so the
shared-QP friction state is single-valued regardless of method.

**Consequence on the byte-exact contract (must be surfaced + accepted).**
Method-invariance and bit-exactness are mutually exclusive here: the reconcile
makes the non-owner adopt the owner's value, changing TPV102/104 by ~1e-14. This
is **physically negligible** and is the reconcile *removing* a latent cross-rank
inconsistency, i.e. TPV102/104 become **cross-rank bit-identical** where today
they differ at ~1e-14. So the regression contract changes:
- OLD: TPV102/104 bit-identical to the pre-fix binary.
- NEW: TPV102/104 (a) physically identical, `worst_rel ≤ 1e-13` vs pre-fix, AND
  (b) now cross-rank bit-identical (rank-A DOFData == rank-B DOFData to 0 ULP) —
  a strict improvement, verifiable by the R-101 verify reporting `max_rel_diff=0`.

**Suggested fix (design).**
```diff
- Gate the whole exchange on fault_friction_law_ ∈ {LSW, LSW_ForcedRupture};
- rate-state skips it (byte-exact TPV102/104).
+ Apply the reconcile for ALL friction laws (LSW, LSW_ForcedRupture, RateState)
+ in ComputeADERSharedFaceFluxRHS. The shared-fault DOFData is made single-valued
+ method-invariantly. Regression contract: TPV/BP5 physically-exact (≤1e-13) AND
+ now cross-rank bit-identical (R-101 max_rel_diff == 0), not bit-identical to the
+ pre-fix binary.
```

**Test case.**
```
test_shared_fault_reconcile_method_invariant (np=2):
  for law in {RateState, LSW}:
     drive a 2-tet shared-fault QP to a state where the two ranks' bulk-Q inputs
     differ by an INJECTED 1 ULP; run the reconcile;
     assert rankA.DOFData == rankB.DOFData to 0 ULP (bit-identical) for all 8
     R-101 fields. (RED before the fix for both laws; GREEN after — proves
     method-invariance, not just LSW.)
```

---

### [R-002] CRITICAL [BUG] sequencing — reconcile must precede flux-RHS assembly

**Category:** BUG

**Description.** The design reconciles the DOFData *after* `EvaluateADER_LSW`,
but the per-side flux RHS assembly (`rhs += assemble_sign · w · shape · F_h_side`,
`wave_operator.inl` ~3349-3357 / the ADER analog ~4788+) consumes the **imposed
state** `I_imp_*_g` derived from the *pre-reconcile* solve. If the reconcile
updates only the stored DOFData but the flux already used the non-owner's
divergent imposed state, the assembled fault RHS is still inconsistent across
ranks → the bulk re-diverges at ~1e-14 each step → the kink re-trips later. The
reconcile only fixes the R-101 *check* (which reads DOFData), not the physics.

**Trigger.** Every sub-step: the flux RHS uses `I_imp` computed before the
reconcile.

**Expected.** The reconcile (of both DOFData AND the imposed state) must complete
**before** the flux RHS assembly that consumes `I_imp_*_g`, OR the assembly must
read the reconciled (owner's) imposed state.

**Suggested fix.** Restructure the shared-fault QP loop into two passes, or defer
assembly: (1) solve + exchange + reconcile DOFData and `I_imp` for all shared
fault QPs; (2) assemble the RHS from the reconciled `I_imp`. Document the
ordering invariant in a comment.

**Test case.**
```
test_R002_reconcile_before_assembly (np=2, LSW at the kink, injected 1-ULP seed):
  after one sub-step, assert the assembled fault-face RHS contribution is
  bit-identical on both ranks (not just the DOFData). RED if assembly uses the
  pre-reconcile imposed state.
```

---

### [R-003] MODERATE [ASSUMPTION] owner selection must be globally unique + robust

**Category:** ASSUMPTION

**Description.** The design picks the owner as the rank with
`elem1_on_plus==true`. That relies on the side fix's "exactly one rank true per
shared face" invariant holding at *every* geometry — including the degenerate
band (fault normal ⟂ ref_normal, θ≈90°) where the projection margin is small. If
that invariant ever fails (two-true or zero-true), the reconcile has two owners
or none → undefined / deadlock.

**Expected.** A deterministic, globally-unique owner key independent of the
side-projection: e.g. owner = the rank holding the element with the smaller
**global element id** of the two sharing the face (or smaller sorted
global-vertex-key). Plus `MFEM_VERIFY` that exactly one rank claims ownership per
shared QP.

**Suggested fix.** Derive `is_owner` from the global element id comparison
(available via `ParMesh` face-neighbour global numbering), not from
`elem1_on_plus`. Keep `elem1_on_plus` only for the +/- flux role.

**Test case.**
```
test_R003_owner_unique (np=2, swept fault tilt incl. θ=90°):
  assert exactly one rank reports is_owner==true for every shared fault QP at
  every tilt (no two-owner / zero-owner).
```

---

### [R-004] MODERATE [BUG] imposed-state reconcile — rebuild leaves a re-seeding residual

**Category:** BUG

**Description.** Design choice (iii) offers "non-owner rebuilds `I_imp` from the
reconciled DOFData + its LOCAL `(Q_plus,Q_minus)`." Its local `(Q_plus,Q_minus)`
differ from the owner's at ~1e-14, so the rebuilt `I_imp` differs at ~1e-14 →
the assembled flux re-seeds the cross-rank bulk divergence every step. Over a
long run this re-accumulates and can re-trip the kink — i.e. the fix would only
delay the failure.

**Expected.** Broadcast the **owner's** canonical-frame imposed state
(`I_imp_plus`, `I_imp_minus`) alongside the DOFData; both ranks rotate the
owner's `I_imp` to global with `T_can` (bit-identical on both, since `can_*` is
bit-identical) for their own side's assembly. Then the assembled RHS is
cross-rank consistent to 0 ULP.

**Suggested fix.** Adopt design choice (ii)=broadcast `I_imp`; drop the "rebuild"
option.

**Test case.** Covered by R-002's test (assembled RHS bit-identical across ranks).

---

### [R-005] MODERATE [BUG] MPI collective safety — must not deadlock (R-1600 class)

**Category:** BUG

**Description.** The per-sub-step shared-fault DOFData exchange is a new MPI
communication on a path with a documented deadlock history (R-1600:
`EvaluateBulkAtFaultQPsCanonical` skipped its `ExchangeFaceNbrData` collective on
ranks without fault shared faces, deadlocking against ranks that have them). A
naive exchange repeats that bug: ranks with no shared fault faces must still
participate consistently, and the send/recv pairing must match the
face-neighbour topology exactly.

**Expected.** Use the same non-blocking point-to-point pattern over the
`ParMesh` face-neighbour group as the existing predictor-Q ghost exchange; every
rank posts sends/recvs for exactly its shared-fault-face neighbours; no
unguarded global collective whose participation depends on local fault presence.

**Suggested fix.** Mirror the existing `ExchangeFaceNbrData` topology; pair
shared QPs by the R-101 face-vertex key (the verify already builds this matcher).
Add an explicit "every rank reaches the same number of exchange rounds" guard.

**Test case.**
```
test_R005_no_deadlock (np=2, 4, 8; partition where some ranks have NO shared
  fault faces): the reconcile completes (no hang) and produces cross-rank
  bit-identical DOFData on the ranks that do share a fault face.
```

---

### [R-006] LOW [QUALITY] per-sub-step exchange cost

**Category:** QUALITY

**Description.** A DOFData+`I_imp` exchange every sub-step adds communication
proportional to the (small) shared-fault-QP count. Likely negligible vs the bulk
solve, but unmeasured.

**Suggested fix.** Measure on the Dc2 8N/400r run; if hot, pack all shared-fault
QPs into one buffer per neighbour (one message per neighbour per sub-step), not
one message per QP. Do NOT gate by a kink-proximity heuristic (that would
re-introduce method/problem specificity — R-001).

---

## Summary
- Critical: 2 (R-001 method-invariance; R-002 reconcile-before-assembly)
- Moderate: 3 (R-003 owner key; R-004 broadcast imposed state; R-005 deadlock)
- Low: 1 (R-006 perf)
- Plan compliance: PARTIAL — the design is sound in mechanism but (a) gated
  method-specifically (R-001) and (b) under-specified on sequencing (R-002) and
  imposed-state (R-004).
- Verdict: **PASS WITH FIXES** — adopt method-invariant reconcile (all laws);
  reconcile DOFData **and** broadcast the owner's imposed state **before** RHS
  assembly; owner by global element id; collective-safe exchange. Regression
  contract relaxes from bit-exact to physically-exact (≤1e-13) + cross-rank
  bit-identical (an improvement). **The byte-exact relaxation needs explicit
  user sign-off** (CLAUDE.md treats TPV/BP5 byte-exactness as non-negotiable;
  the directive overrides it here, but record the decision).

## Unreviewed
- The deleted v5 R-501 pack/unpack code (git history) — not re-read; the design
  recommends a new matcher-keyed pack instead, so it's out of scope.
