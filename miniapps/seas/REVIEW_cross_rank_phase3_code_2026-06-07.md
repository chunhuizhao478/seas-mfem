# Adversarial Code Review — Cross-rank Material Exchange, Phase 3

**Plan:**  `document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md` (Phase 3 section)
**Scope:** the Phase-3 diff only —
- `dynamic/bimaterial_wave_operator.inl` `BuildPerFaceCentralFluxMatrices_` (`n_reclass_shared`
  tally + the extra `MPI_Reduce` + the print-block changes).
- `tests/parallel/test_bimaterial_seam_material_np2.cpp` (`RunSharedContrastReclassify` +
  its two call sites).
**Reviewer:** chunhui-code-reviewer (adversarial pass; assume bugs)
**Date:** 2026-06-07

Phases 1+2 landed/reviewed/fixed (0 critical each). The contrast guard became ACTIVE in Phase 2
(the neighbour material is the TRUE peer). Phase 3 only adds the shared tally + the two gates.

---

## Verification re-run (independent, this review)

Built clean (forced rebuild of the test .o so it reflects the current `.inl`; the `.inl` is a
header-only template included via the `.hpp` into the test TU, so a stale-`.o` ABI hazard is the
only way a "green" run could be misleading — eliminated by `rm` + rebuild):

- **Gate `seas_test_bimaterial_seam_material_np2` (np=2): 131/131, 0 failed.**
  - STRONG (mod_hi=120e9, tol=0.1): log shows `reclassified ... local = 0, shared = 2`
    (= 1 seam face × 2 ranks); `op_off` baseline `shared central faces = 1` (face IS central);
    `op_on` `shared central faces = 0` (reclassified out). IMPL-8 PASSED on both ranks.
  - WEAK (mod_hi=33e9, tol=0.1): `reclassified ... shared = 0`; face STAYS central (`shared
    central faces = 1` in both `op_off` and `op_on`). IMPL-8 PASSED on both ranks.
- **Regressions:** parity **9/9**, central-flux **6/6**, dispatch **10/10**,
  mixed-flux-shared np=2 **16/16**. No regression.

Contrast sizing recomputed independently (Zp=sqrt(3·mod·rho), Zs=sqrt(mod·rho), rho=2670,
lam=mu=mod so Zp and Zs contrasts are *equal*):
- STRONG: ContrastValue(30e9,120e9) = **0.500000** > 0.1 → reclassify. Not borderline.
- WEAK:   ContrastValue(30e9,33e9) = **0.046537** < 0.1 → stays. Not borderline.
Both are an order of magnitude clear of the tol. The test is robustly non-borderline.

---

## Attack results on the named hazards

### 1. `MPI_Reduce(&n_reclass_shared, ...)` collective safety — CLEAN

`bimaterial_wave_operator.inl:975-976`. The reduce is inside `if (guard_on)` (line 953), which is
inside `if constexpr (IsParallelMesh<MeshType>::value)` (line 968). `guard_on = (contrast_tol >=
0.0)` (line 731) is derived from `mixed_flux_contrast_tol_`, a config/CLI scalar set identically on
all ranks before `SetMixedFluxMode`. So `guard_on` is rank-uniform and `IsParallelMesh` is a
compile-time constant identical on every rank ⇒ the reduce is reached by ALL ranks or NONE. No
rank-conditional reduce, no deadlock. A rank with `n_reclass_shared == 0` (e.g. a rank with no
shared reclassify) still calls the reduce with value 0 — correct collective participation, verified
by the strong-case log (rank 0 contributes 1, rank 1 contributes 1, SUM = 2). The four reduces
(hist, n_local_hist, n_reclass, n_reclass_shared) are issued in the same fixed order on all ranks.

### 2. Cross-seam symmetry of the reclassify decision — CLEAN (this is the load-bearing one)

The decision on each rank is `IsStrongContrast(flux_local, flux_nbr, tol)` (line 901-902) where
`flux_local = owned_flux_pool_->At(local_elem)` and `flux_nbr` is built from
`shared_face_neighbour_material_[sf]` (the TRUE peer, Phase 2). For a 2-rank seam between local
elements L_A (rank A) and L_B (rank B):

- Rank A evaluates `IsStrongContrast(flux(L_A), flux(peer=L_B))`.
- Rank B evaluates `IsStrongContrast(flux(L_B), flux(peer=L_A))`.

`ContrastValue(a,b)` (`godunov_flux_bimaterial.cpp:423-434`) is **symmetric in (a,b)**:
`max_zp = max(zp_a,zp_b)`, `c_zp = |zp_a - zp_b| / max_zp` (and the Zs analogue), `return
max(c_zp, c_zs)`. `IsStrongContrast` is `tol >= 0 && ContrastValue(a,b) > tol` (line 436-442) with
a rank-uniform `tol`. Therefore both ranks compute the IDENTICAL boolean ⇒

  (a) `n_reclass_shared` is symmetric (each adjacent rank counts the seam face once) ⇒ the
      MPI_SUM "counts twice" semantics is EXACTLY what the label claims (empirically `shared = 2`);
  (b) the **dispatch is consistent across the seam** — both ranks erase the face from their own
      `central_flux_face_set_`, so `SharedInteriorFaceFlux_` (inl:1080) takes the upwind branch on
      BOTH sides. There is NO "one side central, one side upwind" hazard.

Robustness of (a)/(b) further requires that rank A's `flux(L_A)` (from the LOCAL flux pool,
`MaterialAtLocal_` at L_A's reference centroid) equals rank B's `flux(peer=L_A)` (from
`MaterialAtNbr_` at L_A's reference centroid via the face-nbr transform). For Constant/Coefficient
this is the same coefficient evaluated at the same physical centroid (globally evaluable, no MPI);
for GridFunction it is the same nodal data interpolated at the same reference centroid through the
library `GetValue` face-nbr path. Phase 1's `accessor_local_peer_all_modes` gate already pins
local==peer to 1e-12 per component. So the two ranks' `ContrastValue` inputs match and the verdict
is bit-deterministic across the seam. **No asymmetry bug.**

### 3. IMPL-8 after a shared reclassify — CLEAN

A reclassified shared face takes `continue` at inl:913 BEFORE `per_face_central_flux_[mesh_face_idx]`
is ever touched (the insert is at line 919, after the guard) ⇒ it is **never** inserted into
`per_face_central_flux_`. The post-loop erase (line 936-948) removes it from
`central_flux_face_set_` AND first MFEM_VERIFYs the upwind fallback exists. The IMPL-8 assertion in
`SetMixedFluxMode` (inl:1196) is EXACT equality `per_face_central_flux_.size() ==
central_flux_face_set_.size()`, checked AFTER the build+erase, so a reclassified face is absent from
BOTH containers and the equality holds. Verified: IMPL-8 PASSED on both ranks in both strong and
weak cases. The gate's IMPL-8 assertion exercises the real reclassified state — in the strong case
`op_on`'s `central_flux_face_set_` genuinely shrank by the seam face (log: `shared central faces`
goes 1→0), so the equality is a live check, not vacuous.

### 4. Gate vacuity — CLEAN (with one MINOR note, P3-002)

- `op_off` (guard OFF) and `op_on` (guard ON, same tol) are constructed on the SAME `ppm` but each
  holds its OWN `central_flux_face_set_` / `per_face_central_flux_` (instance members inherited via
  `using Base::`). `op_off`'s construction + `SetMixedFluxMode` only READS the mesh (re-fetches its
  own shared-face transforms during its own build); it does not mutate `ppm` topology, so `op_on`
  sees an identical mesh. `mf = GetSharedFace(sf)` is computed ONCE before either operator
  (lines 680-681) from pure ParMesh topology — stable across both constructions.
- The reclassify assertion is non-vacuous: it is gated by `if (did)` where `did = baseline_central`
  (`op_off` has the face central on this rank). The `MPI_Allreduce(gdid>=1)` PROVES some rank had
  the baseline-central face, and the strong-case log confirms `op_off` `shared central faces = 1`
  (not 0) — so the baseline is real, and `op_on` then drops it to 0. The face cannot pass for the
  "wrong reason" of being non-central in the baseline.
- The `expect(false)` case asserting `gdid >= 1` would NOT mask a missing-face bug: if the face
  were silently absent on every rank, `gdid` would be 0 and that assert fails.
- IMPL-8 is asserted on EVERY rank unconditionally (line 695-698), independent of ownership, so even
  the rank that does not own the seam corridor face contributes a real IMPL-8 check.

### 5. guard-OFF (default tol<0) byte-exactness — CLEAN

Only ONE new statement lives outside `if (guard_on)`: the stack declaration `long long
n_reclass_shared = 0;` (line 772) — no allocation, no collective, no behavior. Both the
`++n_reclass_shared` (inside `if (guard_on && IsStrongContrast)`, line 901-912) and the
`MPI_Reduce(&n_reclass_shared, ...)` (inside `if (guard_on)`, line 975) are guarded. A default run
(tol<0 ⇒ guard_on=false) issues NO new collective and performs NO new work ⇒ byte-exact with
Phase 2. Confirmed by the unchanged parity 9/9 + central-flux 6/6 (those run guard-OFF).

### 6. Print correctness on rank 0 — CLEAN

`MPI_Reduce` writes the result only on root; the print is gated `hrank == 0` (line 983) so it reads
only the valid reduced `n_reclass_shared_g` / `n_reclass_g` / `n_local_g` / `hist_g`. On non-root
ranks those `_g` locals retain their pre-reduce local-init values but are never printed. The print
also fires only when `n_local_g > 0 || n_reclass_g > 0 || n_reclass_shared_g > 0`, so a guard-on +
mixed_flux=none + no-contrast run does not emit a noise "0/0" line — and, importantly, a run where
the ONLY reclassify is a shared seam face (local = 0) STILL prints, because the gate now includes
`n_reclass_shared_g > 0` (verified: the strong case prints with `local = 0, shared = 2`).

---

## Findings

### P3-001 — `n_reclass_shared_g` MPI_SUM double-counts asymmetrically when np≥3 at a tri-rank corner — INFORMATIONAL (label is honest for np=2; under-specified for np>2)

**File:** `dynamic/bimaterial_wave_operator.inl:959-993`
**Severity:** LOW (diagnostic-only; no physics, no dispatch impact)
**Category:** diagnostic accuracy / documentation

**Description.** The label printed is `"shared = per-rank sum; a 2-rank seam face counts twice"`.
For a standard interior partition seam a shared face is adjacent to exactly 2 ranks, so the
MPI_SUM is exactly `2 × (number of distinct reclassified seam faces)` and the label is correct.
But the value is a per-rank SUM with no de-dup; the plan's Phase-3 bullet explicitly allowed either
"dedup by lower global element id" OR "print local-exact + shared per-rank (R-005)". The
implementation chose the latter (acceptable per the plan). The residual concern is ONLY the label's
generality: a shared *mesh face* in MFEM is always exactly 2-element (so always exactly 2 adjacent
ranks for a cut face) — there is no tri-rank single-face case — so even at np=8 the "counts twice"
statement remains arithmetically correct per distinct reclassified face. **No correctness bug.** The
only thing a reader could misread is that `shared` is NOT a count of distinct faces; the label says
so. This is informational; recommend keeping as-is.

**Trigger:** none (cosmetic).
**Suggested fix (optional, clarity only):**
```diff
-                   << " (shared = per-rank sum; a 2-rank seam face counts twice)\n";
+                   << " (shared = per-rank sum over adjacent ranks; each distinct"
+                      " reclassified seam face is a 2-element cut => counts twice;"
+                      " divide by 2 for the distinct-face count)\n";
```

### P3-002 — `RunSharedContrastReclassify` never asserts the per-rank `n_reclass_shared` tally itself — MINOR (gate covers the EFFECT, not the COUNTER)

**File:** `tests/parallel/test_bimaterial_seam_material_np2.cpp:647-723`
**Severity:** LOW (test-coverage gap, not a code bug)
**Category:** test completeness

**Description.** Phase 3's headline new code path is the `n_reclass_shared` tally + its reduce. The
gate asserts the *consequence* of a reclassify (face leaves `central_flux_face_set_`, IMPL-8 holds)
but NOT the tally value. The "shared = 2" number is only ever eyeballed in the printed banner; no
assertion would fail if a future change made `n_reclass_shared` mis-count (e.g. double-increment, or
increment moved outside `if (guard_on)`) while still erasing the face. Because `n_reclass_shared` is
a private local with no accessor, the gate cannot read it today — so the counter is verified only by
manual log inspection, which is exactly the kind of thing that silently rots.

**Trigger:** a regression that corrupts the tally but preserves the erase ⇒ gate stays green,
banner silently wrong.
**Suggested fix (add a const accessor + assert it; sketch):**
```diff
 // bimaterial_wave_operator.hpp — add (parallel diagnostics; mutable member set in the build)
+   /// (Cross-rank Phase 3) Per-rank count of SHARED corridor faces the contrast
+   /// guard reclassified central->upwind in the last BuildPerFaceCentralFluxMatrices_.
+   long long GetNReclassShared() const { return n_reclass_shared_last_; }
```
```cpp
// in the test, strong case, on a rank owning the seam:
//   TEST_ASSERT(op_on.GetNReclassShared() == 1, "...: this rank reclassified its 1 seam face");
//   int g=0; MPI_Allreduce-SUM -> TEST_ASSERT(g == 2, "...: global per-rank sum == 2 (2-rank seam)");
// weak case: TEST_ASSERT(op_on.GetNReclassShared() == 0, ...).
```
(Promote the local `n_reclass_shared` to a mutable member `n_reclass_shared_last_` written at the
end of the build, so the value survives to the accessor.)

### P3-003 — histogram excludes shared faces, so a guard-on run with ONLY shared-contrast faces prints `histogrammed = 0` while still reporting a reclassify — INFORMATIONAL (intentional, R-005; document only)

**File:** `dynamic/bimaterial_wave_operator.inl:762-768, 909-910`
**Severity:** LOW (intentional design to avoid MPI double-count; worth a one-line note)
**Category:** diagnostic completeness

**Description.** Shared faces are deliberately NOT histogrammed (to avoid the across-seam
double-count, R-005), only locally. In the strong case the printed banner is
`local corridor faces histogrammed = 1 / bins = 1 0 0 0 0 / reclassified local = 0, shared = 2`.
A reader could be puzzled that the histogram shows the single *local* face in the [0,1)% bin (no
contrast there) while the reclassify line reports a strong (`shared`) contrast that never appears in
any bin. This is correct (the contrasting face is the shared seam face, excluded from the
histogram), but the connection is non-obvious. Not a bug; a comment in the print block would help an
operator reading the log on Frontera.

**Trigger:** none.
**Suggested fix:** add to the printed block (rank 0):
`"  (histogram is LOCAL corridor faces only; shared-seam contrasts are tallied as 'shared' above, "
` "not histogrammed, to avoid the across-seam MPI double-count)\n"`.

---

## What I explicitly tried to break and could NOT

- **Deadlock via rank-conditional reduce:** the reduce is under rank-uniform `guard_on` + a
  compile-time-uniform `if constexpr` ⇒ all-or-none. No deadlock.
- **Asymmetric reclassify across the seam:** `ContrastValue` is provably symmetric and `tol` is
  rank-uniform; both ranks read each other's TRUE material (Phase-1-gated local==peer to 1e-12) ⇒
  identical boolean ⇒ identical dispatch on both sides. No "central on one side, upwind on the
  other" inconsistency.
- **IMPL-8 violation:** reclassified shared face `continue`s before insert and is erased from the
  set; exact size-equality holds; verified on both ranks, both cases.
- **Vacuous gate pass:** baseline `op_off` proven to hold the face central (log `shared central
  faces = 1`); `gdid>=1` proves real ownership; `op_on` proven to drop it (1→0). Two independent
  operators on one mesh have independent state (instance members); no cross-contamination; `mf`
  stable.
- **guard-OFF regression:** only a stack init added outside `if (guard_on)`; parity/central-flux/
  dispatch all byte-exact green.
- **Borderline contrast sizing:** strong 0.500 / weak 0.0465 vs tol 0.1 — both ~5-10× clear.

---

## Proposed Unit Tests (additions, all LOW priority — current gate is sufficient for merge)

### Test NEW-1: `shared_reclassify_tally_value` (HIGH if P3-002 fix landed; else N/A)
- Target: `BimaterialWaveOperator::BuildPerFaceCentralFluxMatrices_` (the `n_reclass_shared` tally).
- Validates: the per-rank shared tally is exactly 1 on each adjacent rank (global SUM == 2) for the
  strong case and 0 for the weak case.
- Priority: MEDIUM. Requires the `GetNReclassShared()` accessor from P3-002. Closes the
  log-only-verification gap.
- Sketch: see P3-002 diff.

### Test NEW-2: `shared_reclassify_dispatch_consistent_both_sides`
- Target: `SharedInteriorFaceFlux_` after a strong-contrast reclassify.
- Validates: BOTH adjacent ranks route the seam face through the UPWIND branch (not central) — i.e.
  `central_flux_face_set_.count(mf) == 0` on the rank that owns Elem1 AND on the peer rank's view —
  directly exercising the cross-seam symmetry argument (currently only inferred).
- Priority: LOW (the symmetry is provable from `ContrastValue` symmetry; this would make it
  observable). Assert `op_on.GetCentralFluxFaceSet().count(mf) == 0` on every rank that has `mf` in
  its shared-face list, reduced.

### Test NEW-3: `guard_on_no_contrast_byteexact_vs_guard_off`
- Target: `BuildPerFaceCentralFluxMatrices_` with guard ON but a seam-continuous (no-contrast)
  material.
- Validates: turning the guard ON with no strong contrast leaves `central_flux_face_set_` and
  `per_face_central_flux_` bit-identical to the guard-OFF build (the guard adds the histogram/reduce
  but reclassifies nothing) — guards against a future change that erases or rebuilds central faces
  spuriously under guard-on.
- Priority: LOW.

---

## Summary

| ID | Severity | Category | Status |
|----|----------|----------|--------|
| P3-001 | LOW | diagnostic label generality | OPEN (cosmetic) |
| P3-002 | LOW | test coverage (tally value un-asserted) | OPEN |
| P3-003 | LOW | diagnostic/log clarity | OPEN (cosmetic) |

- **Critical issues: 0**
- Moderate: 0. Low/minor/informational: 3 (P3-001 label, P3-002 test gap, P3-003 log clarity).
- All named attack vectors (MPI collective safety, cross-seam reclassify symmetry, IMPL-8 after
  shared reclassify, gate vacuity, contrast sizing, guard-OFF byte-exactness) were exercised and
  found CLEAN.
- Independent re-run: gate **131/131**; regressions parity **9/9**, central-flux **6/6**, dispatch
  **10/10**, mixed-flux-shared np=2 **16/16**.

**Verdict: MERGEABLE.** Phase 3 is correct. The tally double-count is the honest, intended
per-rank-sum semantics (the label is accurate for the 2-element-cut shared-face reality at any np);
the contrast guard sees real contrast and reclassifies symmetrically on both sides of the seam; the
extra MPI_Reduce is collective-safe and fully gated under `guard_on`; default (guard-OFF) runs are
byte-exact. The 3 LOW findings are diagnostic-clarity and test-coverage hardening — none block
merge. Recommended (non-blocking): land P3-002 (assert the tally value, not just its effect) so the
"shared = 2" number stops being log-only.
