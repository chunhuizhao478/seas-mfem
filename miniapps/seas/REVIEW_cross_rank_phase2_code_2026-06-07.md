# Cross-rank Material Exchange — Phase 2 Code Review

Adversarial review of Phase 2 (route flux pool + shared neighbour through the
uniform accessor; remove the seam-continuity abort) against
`document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md`
(Phase 2 section).

Files reviewed (modified in Phase 2):
- `dynamic/bimaterial_wave_operator.inl` (ExchangeBiMaterialNeighbours_ rewrite;
  BuildPerFaceCentralFluxMatrices_ guard removal; ctor/accessor scaffolding)
- `dynamic/bimaterial_wave_operator.hpp` (accessor decls, exchange-count getter)
- `tests/parallel/test_bimaterial_seam_material_np2.cpp` (NEW Phase-2 gates)
- `tests/parallel/test_bimaterial_mixed_flux_shared.cpp` (check (c) flipped)

Cross-checked consumers:
- `dynamic/bimaterial_wave_operator.inl` `BuildPerFaceBimaterialFluxMatrices_`
  (the UPWIND shared build — the second consumer of the neighbour material)
- `drivers/spatial_dyn_driver.cpp` (seam_continuous propagation; banner)
- `dynamic/heterogeneous_material.{hpp,cpp}` (EvalAt / GetIntPoint contract)
- `dynamic/godunov_flux_bimaterial.cpp` (IsStrongContrast symmetry)
- `fem/pgridfunc.cpp` (ParGridFunction::GetValue face-nbr index convention)
- `tpv31/configs/tpv31.toml` (production flux/material mode)

**Review history:**
- **Rev 1** (2026-06-07): 7 findings — 0 CRITICAL, 3 MODERATE, 4 MINOR/SUGGESTION.

---

## Independent verification (re-run, not trusted from the brief)

Built with the worktree override + conda `mfem-dev`. All claims reproduced:

| Test | Result |
|------|--------|
| `seas_test_bimaterial_seam_material_np2` (np=2) | **109/109** |
| `seas_test_bimaterial_mixed_flux_shared` (np=2, flipped check c) | **16/16** |
| `seas_test_bimaterial_fault_perside_material` | **11/11** |
| `seas_test_phaseh_wave_operator_constant_parity` | **19/19** |
| `seas_test_bimaterial_wave_operator_parity` | **9/9** |

The `shared_central_matches_serial` gate matched to `worst rel 0.000000` (not a
false match — log shows "shared central faces = 1" and the serial reference built
its interior central face, so both entries are genuinely present and non-empty).

**Verdict on the build-level claims: confirmed.** The findings below are about
behavior the green tests do NOT cover, and about a re-baseline the test suite
silently endorses.

---

## Findings — Rev 1

### P2-001: Depth-profile parallel UPWIND seam matrices silently re-baseline (TPV31/TPV102) — the "byte-exact for Coefficient" guarantee does not hold for a cross-depth partition

**Severity:** MODERATE (correctness re-baseline of production gold; not a crash,
arguably an improvement — but undocumented-by-test and contradicting a plan claim)
**Category:** byte-exactness regression / correctness re-baseline
**File:** `dynamic/bimaterial_wave_operator.inl`
  - `ExchangeBiMaterialNeighbours_` (lines ~317-329, the rewrite)
  - consumed by `BuildPerFaceBimaterialFluxMatrices_` Pass 2 (lines 585-599)

**Description.** The task framing and the plan (`Constraints` line 66, `Regression`
lines 252-253) assert TPV31/TPV102 stay **byte-exact** and that the bulk re-baseline
"is the central path only / the fault is unchanged." That is **not** what the diff
does. `shared_face_neighbour_material_[sf]` is consumed by **two** builds:

1. `BuildPerFaceCentralFluxMatrices_` (line 825) — the central/mixed-flux path.
2. `BuildPerFaceBimaterialFluxMatrices_` Pass 2 (line 585) — the **UPWIND** shared
   Riemann matrix `per_face_bimaterial_flux_[mesh_face_idx][0]`.

TPV31 production (`tpv31/configs/tpv31.toml`) runs `mixed_flux = "none"` +
`interior_flux = "matrix"` with `kind = "depth_profile_1d"` (Mode::Coefficient,
`f(z)`). So the **central path is never built** (empty `central_flux_face_set_`),
but the **upwind** path IS — and it now consumes the rewritten neighbour material.

- OLD stub: `shared_face_neighbour_material_[sf] = per_elem_lmr_[local_elem]`
  ⇒ `flux_nbr == flux_local` ⇒ the shared upwind matrix was an effectively
  *homogeneous* Riemann at the seam (no impedance contrast).
- NEW: `flux_nbr` = the **peer element's centroid material** (peer's own `z`).
  For any ParMETIS seam that separates elements at **different depths** (the
  common 3D case), `flux_nbr != flux_local` ⇒ a genuine (depth) bi-material
  upwind Riemann ⇒ the shared-face numerical flux **changes**.

The new value is *more physically correct* (the old stub used the local element's
material as the neighbour's, which is wrong at a depth-varying seam). But:
- It **re-baselines** every TPV31/TPV102 **pure-upwind ADER** parallel gold run
  whose partition cuts across depth — i.e. the exact production configuration.
- The plan's "byte-exact for Coefficient (TPV31/TPV102)" is therefore only true
  for a *seam-vertical* (constant-z-across-seam) partition. The diff's own
  comment (inl line ~293, "CORRECT (was wrong) wherever it varies laterally")
  acknowledges the change, but the **headline plan constraint and the brief's
  byte-exact claim are overstated** — depth varies *across* seams in general,
  and that is the bulk path, not just the fault.

**Trigger.** TPV31/TPV102 at np>1 on any mesh where a partition seam separates
two elements at different `z`. With `interior_flux=matrix, mixed_flux=none`, the
shared upwind matrix at that seam changes vs. the pre-Phase-2 binary.

**Why not CRITICAL.** The change is in the *correct* direction and the fault per-
side impedances (the physically dominant quantity) are untouched in Phase 2.
But it IS a result-changing re-baseline of production gold that no test asserts,
and the constant-parity / wave-operator-parity suites do NOT catch it (they are
Constant-material and serial, respectively).

**Suggested fix (documentation + a gating test, not a behavior change):**
```diff
--- a/document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md
+++ b/document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md
@@ Constraints
-- **Byte-exact for Constant + Coefficient** (TPV205/BP5/TPV6/TPV31/TPV102): ...
+- **Byte-exact for Constant + Coefficient** ONLY where the partition seam is
+  material-continuous across it (Constant always; a depth profile only when the
+  seam does NOT separate elements at different depths).  A depth-profile seam
+  that cuts across depth RE-BASELINES the BULK shared UPWIND Riemann matrix
+  (per_face_bimaterial_flux_, used by interior_flux=matrix / mixed_flux=none —
+  i.e. the TPV31/TPV102 production path), from the old homogeneous stub to the
+  true (depth) bi-material flux.  The fault per-side impedances are unchanged.
```
Add the gating test in P2-T1 (below): build the depth-profile-on-a-cross-depth
z-seam at np=2 and assert the shared upwind matrix equals the serial 2-sided
upwind matrix (analogous to `shared_central_matches_serial`, but for the upwind
path that production actually uses). Re-baseline the TPV31/TPV102 parallel gold
on Frontera and document the shift.

---

### P2-002: Removed rank-0 WARNING leaves the seam re-baseline (P2-001) with no runtime signal

**Severity:** MODERATE (observability / silent result change)
**Category:** lost diagnostic
**File:** `dynamic/bimaterial_wave_operator.inl` `ExchangeBiMaterialNeighbours_`
(the deleted R-004 `if (Coefficient && n_shared>0) ... WARNING` block)

**Description.** Pre-Phase-2, a parallel Coefficient run printed a loud rank-0
WARNING that the seam neighbour material was a stub. Phase 2 deletes it (correctly
— the stub is gone). But because the new behavior **changes results** at cross-
depth seams (P2-001) AND silently builds central flux across strong contrasts
when the guard is off (P2-003), there is now **no runtime breadcrumb** that the
parallel material treatment differs from serial. A user comparing a new TPV31
parallel run against an old gold gets a silent numerical shift with nothing in the
log explaining it.

**Trigger.** Any parallel Coefficient run; the user only notices via a diff of
results, with no operator-log hint.

**Suggested fix.** Replace the deleted WARNING with a one-line *informational*
rank-0 banner (not a warning) that names the new, correct behavior, e.g.:
```diff
+      // (Cross-rank Phase 2) One-line informational banner: the seam neighbour
+      // material is now the TRUE peer (was a local-side stub).  Re-baselines the
+      // bulk shared Riemann at depth-/lateral-varying seams vs the old binary.
+      if (material_ != nullptr
+          && material_->mode != MaterialField::Mode::Constant
+          && n_shared > 0)
+      {
+         int rank = 0; MPI_Comm_rank(pmesh.GetComm(), &rank);
+         if (rank == 0)
+            mfem::out << "[wave_operator] cross-rank seam material = TRUE peer "
+                         "(MaterialAtNbr_); bulk seam Riemann re-baselined vs the "
+                         "pre-2026-06 local-side stub.\n";
+      }
```

---

### P2-003: Removing the seam_continuous abort opens a window where strong-contrast central flux on a SHARED face builds silently (guard OFF by default) — Phase-2/Phase-3 safety gap

**Severity:** MODERATE (latent stability; plan-deferred to Phase 3 but the safety
net is removed BEFORE its replacement lands)
**Category:** correctness/stability gap (central, non-dissipative flux across a
strong material jump)
**File:** `dynamic/bimaterial_wave_operator.inl` `BuildPerFaceCentralFluxMatrices_`
(removed `MFEM_VERIFY(mode==Constant || seam_continuous_)`, lines ~811-818;
contrast guard at lines 837-849)

**Description.** Pre-Phase-2, a central (mixed) flux on a non-Constant shared face
**aborted loudly** unless the user affirmed `seam_continuous=true`. Phase 2 removes
that abort. The replacement — the contrast guard (`IsStrongContrast` → reclassify
to upwind) — is **OFF by default** (`mixed_flux_contrast_tol_ = -1` ⇒
`guard_on=false`, line 702) and is only *intended to be activated in Phase 3*
(plan lines 171-179). So in the Phase-2 state of the tree, a parallel
`mixed_flux=adjacent` run with a genuine lateral/CVM contrast at a seam will now
**silently build a non-dissipative central flux across a strong material jump** —
no abort, no guard, no warning.

This is *consistent* with the interior 2-sided path (which also builds central flux
across strong contrasts when the guard is off), so it is not a new asymmetry, and
the plan explicitly defers the guard to Phase 3. But the loud abort that used to
force user acknowledgement is gone in the interim. A user who relied on the abort
(e.g. ran a lateral-CVM `mixed_flux=adjacent` parallel job expecting it to refuse)
now gets a silent, potentially unstable, central flux across the jump.

**Trigger.** `mixed_flux=adjacent` (or `all_continuous`) + a non-Constant material
that varies laterally across a partition seam + default `mixed_flux_contrast_tol`
(< 0). No abort; non-dissipative central flux at the strong-contrast seam.

**Suggested fix.** Either (a) land Phase 3's guard-on default for the matrix +
non-Constant + shared-central combination, or (b) keep a *soft* default: when
`mixed_flux != none` AND material is non-Constant AND `contrast_tol < 0`, emit a
rank-0 WARNING that the contrast guard is disabled and central flux may be built
across seam contrasts. Minimal interim diff:
```diff
+         // (Cross-rank Phase 2/3 gap) Guard is OFF by default; a strong seam
+         // contrast would build a non-dissipative central flux silently.
+         if (!guard_on && material_->mode != MaterialField::Mode::Constant)
+         {
+            static thread_local bool warned = false;
+            int r = 0; MPI_Comm_rank(pmesh.GetComm(), &r);
+            if (r == 0 && !warned) {
+               warned = true;
+               mfem::out << "[wave_operator] WARNING: central flux on shared "
+                            "faces with a non-Constant material and contrast "
+                            "guard DISABLED (mixed_flux_contrast_tol < 0): a "
+                            "strong seam contrast will build a non-dissipative "
+                            "central flux.  Set mixed_flux_contrast_tol >= 0.\n";
+            }
+         }
```

---

### P2-004: The reclassification COUNT undercounts shared reclassifications; the gating comment is stale and now factually wrong

**Severity:** MINOR (observability; only relevant when the guard is ON — Phase 3)
**Category:** stale comment + diagnostic undercount
**File:** `dynamic/bimaterial_wave_operator.inl` lines 845-846 (shared loop, no
`++n_reclass_local`), 893-896 + 920-921 (the reduction + the print)

**Description.** The shared-face contrast-guard branch (line 847)
`reclassified.push_back(mesh_face_idx)` but does **not** increment
`n_reclass_local`. The post-loop reduction (line 895) sets `n_reclass =
n_reclass_local` and the rank-0 print (line 920) reports
"reclassified central -> upwind = N (local+shared; shared excluded from the
histogram)". With Phase 2's TRUE peer material, a shared face CAN now be
reclassified — but it is **not counted**, so the printed total **undercounts** the
actual reclassifications, and the parenthetical "(local+shared)" is wrong (shared
is excluded from the *count*, not just the histogram).

Worse, the comment at lines 894-896 still reads
"shared reclassification is inert under the seam-continuous stub" — that premise
was deleted in Phase 2 (there is no stub anymore), so the justification for
omitting shared faces from the count is now **false**. The post-loop erase (line
882) DOES erase reclassified shared faces from `central_flux_face_set_`, so IMPL-8
holds; only the *reported* count is wrong.

**Trigger.** Guard ON (Phase 3) + a strong contrast on a shared central face: the
face is correctly reclassified and erased, but the operator log reports a smaller
"reclassified" count than reality and prints a contradictory "(local+shared)".

**Suggested fix.** Count shared reclassifications separately (disjoint-by-rank
double-count is the concern, so dedup by lower global element id OR report a
per-rank shared count) and fix the comment. Minimal honest version: add a
`long long n_reclass_shared = 0; ++n_reclass_shared;` in the shared branch, reduce
it with `MPI_SUM`, and print "local=… shared=… (shared possibly double-counted
across the seam)". This is squarely Phase 3 territory (the plan's `n_reclass_shared`
tally, lines 173-174) — flag so it is not lost.

---

### P2-005: `MaterialAtNbr_`/`MaterialAtLocal_` Coefficient peer/local reads do not `SetIntPoint` before `EvalAt` — fragile for any future coefficient that reads `T.GetIntPoint()`

**Severity:** MINOR (latent; NOT a bug for current TPV FunctionCoefficient-based
materials; would silently misbehave for a `GetIntPoint`-reading coefficient)
**Category:** robustness / latent UB-adjacent
**File:** `dynamic/bimaterial_wave_operator.inl` `MaterialAtNbr_` line 200,
`MaterialAtLocal_` line 150, `BuildGodunovFluxPool_` line 230;
contrasted with `heterogeneous_material.cpp` lines 93-117 (`MaxCpInElement`,
which DOES `T->SetIntPoint(&ip)` before `Eval` and documents the `GetIntPoint`
UB hazard).

**Description.** `MaterialField::EvalAt` (heterogeneous_material.hpp:201) calls
`lambda_coef->Eval(T, ip)` with no `SetIntPoint`. For a `FunctionCoefficient`,
`Eval(T, ip)` uses `T.Transform(ip, x)` (the explicit `ip`) — correct. But the
sibling helper `MaxCpInElement` explicitly warns (heterogeneous_material.cpp
lines 93-110) that "some coefficients rely on `T.GetIntPoint()` rather than the
explicit `ip`" and therefore calls `T->SetIntPoint(&ip)` first. The new accessors
(and the bulk pool loop) call `EvalAt` WITHOUT setting the int point, so:
- `MaterialAtNbr_` (line 200) passes `*ftr->Elem2` whose `IntPoint` is whatever
  the last `GetSharedFaceTransformations`/`SetAllIntPoints` left it at — NOT
  `ip_peer`.
- `MaterialAtLocal_` / `BuildGodunovFluxPool_` (lines 150, 230) likewise.

For the TPV materials (FunctionCoefficient / depth profile via `FunctionCoefficient`)
this is harmless and byte-exact with the legacy pool (which also did not
`SetIntPoint`). But it is a footgun: a future analytic CVM coefficient that reads
`T.GetIntPoint()` (e.g. a coefficient composed with a `GridFunctionCoefficient`)
would silently read the material at the WRONG point on the peer transform — a
silent wrong-material-at-seam bug, the exact class Phase 2 is supposed to close.

**Trigger.** Any Coefficient that internally reads `T.GetIntPoint()` instead of the
explicit `ip`. None today; a real hazard for the CVM/analytic future the plan
targets.

**Suggested fix.** Make the accessors set the int point before delegating, matching
`MaxCpInElement`:
```diff
   else
   {
+     // Match MaxCpInElement: some coefficients read T.GetIntPoint() rather than
+     // the explicit ip (heterogeneous_material.cpp:93-110).  Set it for safety.
+     ftr->Elem2->SetIntPoint(&ip_peer);
      material_->EvalAt(ftr->Elem2No, *ftr->Elem2, ip_peer, lam, mu, rho);
   }
```
and analogously `T.SetIntPoint(&ip)` before the `EvalAt` in `MaterialAtLocal_`.
(Byte-exact for FunctionCoefficient, since `Eval` uses the explicit `ip`; closes
the hazard for `GetIntPoint`-reading coefficients.)

---

### P2-006: `seam_continuous_` is now a dead behavioral flag (no consumer) but is still set + echoed by the driver

**Severity:** SUGGESTION (plan-sanctioned until Phase 5; flag so it is not lost)
**Category:** dead/inert state
**File:** `dynamic/bimaterial_wave_operator.{hpp,inl}` (`seam_continuous_`,
`SetSeamContinuous`); `drivers/spatial_dyn_driver.cpp` lines 1205-1210, 1333.

**Description.** The central-build `MFEM_VERIFY(mode==Constant || seam_continuous_)`
was the ONLY behavioral consumer of `seam_continuous_`. With it removed, the flag
is inert: the driver still parses `[material].seam_continuous`, calls
`SetSeamContinuous`, and prints `(seam_continuous=true/false)` in the mixed-flux
banner (line 1333) — but nothing reads it. The plan defers formal deprecation to
Phase 5 (lines 223-228), so this is expected; flagging only so the banner does not
mislead users into thinking the flag still does something. No action required in
Phase 2 beyond noting it; Phase 5 should drop the banner reference.

**Suggested fix.** Phase 5: drop the `(seam_continuous=...)` from the banner or
re-label it "(deprecated, no effect)".

---

### P2-007: New Phase-2 gates do not value-check the UPWIND shared matrix for a non-Constant material — the production path (P2-001) is only covered indirectly

**Severity:** SUGGESTION (coverage gap)
**Category:** test coverage
**File:** `tests/parallel/test_bimaterial_seam_material_np2.cpp`;
`tests/parallel/test_bimaterial_mixed_flux_shared.cpp` (Test 3.4c is
`value_check=false`).

**Description.** The Phase-2 gates assert (a) `shared_face_neighbour_material_`
matches the independent peer expectation, and (b) the **central** shared matrix
equals the serial interior central matrix (`shared_central_matches_serial`). But
the path TPV31/TPV102 production actually uses is the **upwind** shared matrix
`per_face_bimaterial_flux_` (interior_flux=matrix, mixed_flux=none) — and no gate
value-checks that for a depth/lateral material. Test 3.4c builds a depth-only z-seam
but with `value_check=false` (only "no false abort") and via the central path. The
re-baseline in P2-001 is therefore endorsed by green tests without being directly
verified. See proposed test P2-T1.

**Suggested fix.** Add P2-T1 (below).

---

## Proposed Unit Tests

### Test Suite: `tests/parallel/test_bimaterial_seam_material_np2.cpp`

**P2-T1: `shared_upwind_matches_serial_depthprofile` (HIGH)** — closes P2-001/P2-007.
- Target: `BimaterialWaveOperator::BuildPerFaceBimaterialFluxMatrices_` (Pass 2,
  shared upwind) via `GetPerFaceBimaterialFlux()`.
- Validates: a depth-profile `f(z)` on a partition seam that **cuts across depth**
  (a z-seam) builds the shared UPWIND Riemann matrix equal (1e-12) to the serial
  2-sided upwind matrix built from the two real per-element materials — proving the
  re-baseline is the *correct* bi-material flux, not garbage, AND pinning the new
  production gold so a future regression is caught.
- Setup/assert sketch (mirrors `RunSharedCentralMatchesSerial`, but upwind and
  with `mixed_flux=none`):
  ```cpp
  // depth profile lam=mu=30e9*(1+0.2*z); rho const
  auto fz = [](const Vector& x){ return 30.0e9*(1.0+0.2*x(2)); };
  // serial: 2x5 (x,z) grid, fault x=1; pick the z-interior face at z=3
  BimaterialWaveOperator<Mesh> op_s(serial, p, MakeCoefficient(...), bc);
  // op_s does NOT need SetMixedFluxMode (upwind built in ctor)
  const auto& us = op_s.GetPerFaceBimaterialFlux().at(fser); // [0][0],[0][1]
  // np=2: partition along z at z=3 (cross-depth seam)
  BimaterialWaveOperator<ParMesh> op_p(ppm, p, MakeCoefficient(...), bc);
  const auto& up = op_p.GetPerFaceBimaterialFlux().at(mf);   // shared side-0
  // on the rank owning the lower-z self side, assert MatRel(up[0], us[0][0]) and
  // MatRel(up[1], us[0][1]) <= 1e-12  (peer half == serial Elem2 half)
  ```
- Also assert (negative control) that under a CONSTANT material the same shared
  upwind matrix is byte-identical to the pre-Phase-2 value (== serial homogeneous),
  guarding the Constant byte-exact claim for the upwind path explicitly.

**P2-T2: `shared_strong_contrast_central_no_guard_builds_then_warns` (MEDIUM)** —
covers P2-003.
- Target: `BuildPerFaceCentralFluxMatrices_` with `mixed_flux=adjacent`, a strong
  lateral contrast at the seam, and `contrast_tol < 0`.
- Validates: with the guard OFF the central matrix is BUILT (no abort) on the
  strong-contrast shared face (documents the Phase-2/3 gap), and — once the
  suggested interim WARNING (P2-003) lands — that a rank-0 warning is emitted.
- Note: the BUILD assertion is the load-bearing one (it pins the documented
  Phase-2 behavior); the warning assertion is optional until the fix lands.

**P2-T3: `shared_reclass_count_includes_shared` (MEDIUM)** — covers P2-004.
- Target: the rank-0 reclassification tally in `BuildPerFaceCentralFluxMatrices_`.
- Validates: with the guard ON and a strong contrast on a SHARED central face, the
  reported "reclassified central -> upwind" count includes the shared face (or a
  separate shared count is printed), AND the face is absent from BOTH
  `central_flux_face_set_` and `per_face_central_flux_` (IMPL-8 still holds).
- This is mostly a Phase-3 gate; stub it now so the undercount does not survive.

**P2-T4: `coeff_peer_read_independent_of_stale_intpoint` (MEDIUM)** — covers P2-005.
- Target: `MaterialAtNbr_` / `MaterialAtLocal_` Coefficient path.
- Validates: a custom `Coefficient` subclass that reads `T.GetIntPoint()` (not the
  explicit `ip`) returns the value at `ip_peer`, NOT at a previously-set stale int
  point. Fails on the current code (no `SetIntPoint`), passes after the P2-005 fix.
- Sketch:
  ```cpp
  struct IPProbe : Coefficient {
    real_t Eval(ElementTransformation& T, const IntegrationPoint&) override {
      Vector x; T.Transform(T.GetIntPoint(), x); return x(2); } // uses GetIntPoint
  };
  // set ftr->Elem2->SetIntPoint to a WRONG ip, then MaterialAtNbr_(ftr, ip_peer)
  // must return the z at ip_peer, not at the wrong ip.
  ```

---

## Priority Summary

| ID | Severity | Status | Description |
|----|----------|--------|-------------|
| P2-001 | MODERATE | OPEN | Depth-profile parallel UPWIND seam re-baselines (TPV31/TPV102); byte-exact claim holds only for seam-vertical partitions |
| P2-002 | MODERATE | OPEN | Removed rank-0 WARNING ⇒ no runtime signal for the seam re-baseline |
| P2-003 | MODERATE | OPEN | seam_continuous abort removed before the Phase-3 guard is on by default ⇒ silent central flux across strong seam contrast |
| P2-004 | MINOR | OPEN | Reclassification count omits shared faces; "(local+shared)" + "inert stub" comments now false |
| P2-005 | MINOR | OPEN | Accessors don't SetIntPoint before EvalAt ⇒ fragile for GetIntPoint-reading coefficients (latent silent wrong-material) |
| P2-006 | SUGGESTION | OPEN | seam_continuous_ now dead but still echoed in driver banner (Phase-5 cleanup) |
| P2-007 | SUGGESTION | OPEN | No value-check gate for the shared UPWIND matrix (the production path) under a non-Constant material |

### Test Coverage Summary

| Suite | File Under Test | Proposed | Key Findings Covered |
|-------|-----------------|----------|----------------------|
| seam_material_np2 | `bimaterial_wave_operator.inl` (upwind shared) | P2-T1 | P2-001, P2-007 |
| seam_material_np2 | central build (guard off) | P2-T2 | P2-003 |
| seam_material_np2 | reclass tally | P2-T3 | P2-004 |
| seam_material_np2 | accessors (IntPoint) | P2-T4 | P2-005 |
| **Total** | | **4** | |

---

## Code Review Report (Rev 1)

**Date:** 2026-06-07
**Scope:** Phase 2 diff (ExchangeBiMaterialNeighbours_ rewrite; central guard
removal; two test files) + all consumers of `shared_face_neighbour_material_`.
**Reviewer:** chunhui-code-reviewer agent

### Summary

The Phase-2 mechanics are correct and the green-test claims reproduce exactly
(109/109 + all listed regressions). The accessor wiring is sound: `MaterialAtNbr_`
passes the full `Elem2No` to `ParGridFunction::GetValue` (verified against
`fem/pgridfunc.cpp:275`, which re-subtracts `GetNE()`), `ne_ == mesh_.GetNE()`
(wave_operator.inl:49), the `#ifdef MFEM_USE_MPI` / `if constexpr` / `else` nesting
in `MaterialAtNbr_` is well-formed for both serial and non-MPI builds, the ctor
order (`SetupMaterialFaceNbrExchange_` → `BuildGodunovFluxPool_` →
`ExchangeBiMaterialNeighbours_`) is right, and the GF + geometry exchanges are
unconditional/rank-uniform (collective-safe). IMPL-8 is preserved: a reclassified
shared face is pushed to `reclassified`, the central build is skipped (no map
entry), and the post-loop erase removes it from `central_flux_face_set_` — so it is
in NEITHER container (map.size() == set.size() holds). `IsStrongContrast` is
symmetric (`godunov_flux_bimaterial.cpp:423-442`), so both ranks reclassify a shared
contrast face identically ⇒ no dispatch asymmetry across the seam. The flipped
check (c) in `test_bimaterial_mixed_flux_shared.cpp` is a real assertion
(`dmat > 0` is robust: the lateral lam=mu=30e9*(1+0.2x) gives ~6e9 across the seam,
far above any tolerance; `c[0] != c[1]` follows). `shared_central_matches_serial`
matches to 0.0 for the *right* reason (both entries present, non-empty, self side
= lower-x on both serial and parallel).

**No CRITICAL defect found.** The substantive concerns are: (1) the diff
re-baselines the **bulk UPWIND** shared Riemann for depth profiles at cross-depth
seams — the TPV31/TPV102 production path with `mixed_flux=none, interior_flux=matrix`
— which the plan's "byte-exact for Coefficient" headline overstates (it is exact
only for seam-vertical partitions); (2) the deleted WARNING plus the
default-off contrast guard remove both safety nets (loud abort and dissipative
reclassify) for a non-Constant central shared face in the Phase-2/Phase-3 interim;
and (3) the reclassification *count* and two in-code comments are now stale/wrong.
None block merge of Phase 2 as a stepping stone, but P2-001 needs a doc correction
+ a value-checking gate (P2-T1) and a Frontera re-baseline before the new behavior
is treated as gold.

### Code Review Findings

| # | File | Line(s) | Severity | Description |
|---|------|---------|----------|-------------|
| 1 | bimaterial_wave_operator.inl | 317-329, 585-599 | MODERATE | Bulk UPWIND seam re-baseline for depth profiles (P2-001) |
| 2 | bimaterial_wave_operator.inl | ~330 (deleted block) | MODERATE | Lost rank-0 seam diagnostic (P2-002) |
| 3 | bimaterial_wave_operator.inl | 811-849, 702 | MODERATE | Abort removed; guard off by default (P2-003) |
| 4 | bimaterial_wave_operator.inl | 845-846, 894-896, 920-921 | MINOR | Shared reclass uncounted; stale comments (P2-004) |
| 5 | bimaterial_wave_operator.inl | 150, 200, 230 | MINOR | No SetIntPoint before EvalAt (P2-005) |
| 6 | spatial_dyn_driver.cpp | 1205-1210, 1333 | SUGGESTION | seam_continuous_ dead but echoed (P2-006) |
| 7 | tests/parallel/*np2.cpp | — | SUGGESTION | Upwind shared matrix not value-checked (P2-007) |

### Verified-OK (attacked, found sound)

- `Geometries.GetCenter(ftr->Elem2->GetGeometryType())` is the correct peer
  reference IP for the peer centroid (peer geometry type, peer transform). The
  Coefficient peer eval is byte-consistent with the peer rank's own pool for a
  globally-evaluable FunctionCoefficient (confirmed: `shared_central_matches_serial`
  = 0.0).
- GridFunction peer read via `GetValue(Elem2No, ip)` uses the library face-nbr path
  with correct DofTransformation/map-type handling; three GFs need not share a
  space; `material_gf_exchange_count_ == 3` (GF) / `0` (Constant/Coefficient) gated.
- Collective symmetry: `pmesh.ExchangeFaceNbrData()` + the 3 GF exchanges are
  unconditional and rank-uniform (mode is identical on all ranks). No
  rank-conditional collective.
- IMPL-8 (`per_face_central_flux_.size() == central_flux_face_set_.size()`) holds
  after shared reclassification (face in neither container).
- Cross-rank reclassify consistency (symmetric `IsStrongContrast`) ⇒ both ranks
  agree on central-set membership for a shared face.
- Constant byte-exactness (peer==local==const) ⇒ constant-parity 19/19 + the
  np2 `seam-continuous => neighbour == local (byte-exact)` gate.
- Test 3.4a check (c) flip is a genuine assertion (`dmat > 0`, `rel > 1e-9`).

### Documentation Status

- Plan `Constraints` (line 66) and `Regression` (lines 252-253) overstate
  byte-exactness for TPV31/TPV102 — should be scoped to seam-vertical partitions
  (P2-001 suggested diff).
- In-code comments at inl 894-896 ("shared reclassification is inert under the
  seam-continuous stub") and the "(local+shared)" print are now factually wrong
  (P2-004).
- The driver banner `(seam_continuous=...)` echoes a now-inert flag (P2-006).

### Test Status

- **Implemented:** Phase-2 gates present and green (109/109 np2; 16/16 mixed-flux-
  shared). Listed regressions reproduce (9/9, 11/11, 19/19, etc.).
- **Test framework:** Make + `mpirun -np 2`; builds clean under the worktree
  override + conda mfem-dev.
- **Priority tests for merge readiness:** P2-T1 (HIGH — value-check the upwind
  shared matrix for a depth profile, the production path) before relying on the
  re-baselined TPV31/TPV102 parallel gold.

### Recommended Next Steps

1. Correct the plan's byte-exactness claim (scope it to seam-vertical partitions)
   and add P2-T1 to value-check the shared UPWIND matrix for a depth profile.
2. Re-baseline TPV31/TPV102 parallel pure-upwind gold on Frontera; document the
   bulk seam shift as an intended correctness fix.
3. Restore a runtime signal (P2-002 informational banner) so the re-baseline is
   not silent; decide Phase-2/Phase-3 sequencing for the contrast guard default
   (P2-003) so the abort safety net is not absent in the interim.
4. Fix the reclassification count + stale comments (P2-004); add `SetIntPoint`
   before `EvalAt` in the accessors (P2-005) to harden the analytic-CVM future.

---

## Summary

- Findings: **7** (0 CRITICAL, 3 MODERATE, 2 MINOR, 2 SUGGESTION).
- Proposed tests: **4** (P2-T1 HIGH, P2-T2/T3/T4 MEDIUM).
- Independent re-verification: **PASS** (109/109 + all listed regressions).
- Verdict: **Phase 2 is mergeable as a stepping stone.** The mechanics are correct
  and IMPL-8 / collective symmetry / cross-rank reclassify consistency all hold.
  The substantive issues are a (correct, but undocumented-by-test and
  over-claimed-as-byte-exact) re-baseline of the TPV31/TPV102 **upwind** seam path,
  and the removal of two safety nets (abort + diagnostic) in the Phase-2/Phase-3
  interim. Address P2-001 (doc + P2-T1 + Frontera re-baseline) and P2-002/P2-003
  (restore a signal) before treating the new parallel results as gold.

**Critical issues: 0**
