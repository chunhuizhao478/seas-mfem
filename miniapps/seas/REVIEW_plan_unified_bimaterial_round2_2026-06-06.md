# Code Review (Round 2): PLAN_unified_bimaterial_volume_and_fault_2026-06-06 — revised plan

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_unified_bimaterial_volume_and_fault_2026-06-06.md`
  (after Round-1 fixes R-001..R-007 and the C1/C2/B2 run-arm additions).
- Cross-checked against: `dynamic/wave_operator.inl` (Mult + ADER fault conversion),
  `dynamic/spatial_setup.hpp`, `dynamic/fault_face_flux.{hpp,cpp}`,
  `dynamic/bimaterial_wave_operator.{hpp,inl}`, `drivers/spatial_dyn_driver.cpp`,
  `tpv6/benchmark_document/2007RuthRalphletter2.pdf`.
- Fresh adversarial pass (rule 10): re-verified Round-1 fixes AND hunted new issues in the
  added C1/C2/B2 content and the R-001 epsilon-offset rewrite.

## Round-1 fix verification (all present, not regressed)
- R-001 (eps-offset rule): present in B1.  R-002 (locked R/T, signs): present in B0 — the
  corrected coefficients `R_sigma=(Z2-Z1)/(Z1+Z2)`, `T_v=2Z1/(Z1+Z2)` are CORRECT.
  R-003 (square patch): present in C2.  R-004 (displacement accumulator): present in C2.
  R-005 (variant table + per-variant relaxation): present in B0/B2.  R-006 (halfspace
  mapping): present in B3.  R-007 (free-surface corner test): present in C2.

## Verified facts that STRENGTHEN the plan (not findings)
- The ADER fault-face conversion ALSO applies per-side A: `wave_operator.inl:4301-4303`
  `FluxForElem_(elem_plus).Interior(...)` / `FluxForElem_(elem_minus).Interior(...)`,
  parallel to the Mult/RK site at `:3048-3051`.  => TPV6 **Arm 1 (matrix + upwind + ADER
  + per-side fault) is viable** once B2 relaxes the `EvaluateADER_LSW` guard.  B0's table
  should record this row as per-side-A = YES.

## Findings

### [R-101] MODERATE — B1: shared (cross-rank) fault face per-side material must carry the eps-offset value, not the peer's centroid material

**Category:** BUG / ASSUMPTION

**Description:**
R-001 fixed the LOCAL fault per-side read with the eps-offset rule (evaluate each element
at `x_f +/- eps*n`).  But for a SHARED fault face the peer side lives on another rank: the
local conversion uses only the local element (`wave_operator.inl:3715`
`FluxForElem_(e1).Interior(...)`), and B1 says the peer material "comes from the cross-rank
exchange; reuse/extend `ExchangeBiMaterialNeighbours_`."  That exchange currently carries
each element's POOL (centroid-evaluated) material.  Using the peer's CENTROID material for
the fault per-side impedance reintroduces exactly the R-001 bug on shared fault faces: for
a depth-varying (or otherwise position-varying) material, the peer centroid depth != the
fault-DOF depth, so `Zp_minus` is wrong -> spurious across-fault contrast on the seam.
The spatial driver does NOT fault-locality-partition (verified earlier), so ParMETIS CAN
cut the fault -> shared fault faces DO occur.  Benign for TPV6 (uniform per side) but a
correctness bug for SAFS (CVM) and for any depth-varying material whose fault is cut.

Also: A2 (volume contrast guard) legitimately uses the CENTROID exchange (centroids ARE
the right per-element materials for an interior-face contrast), while B1 needs the
eps-offset-at-fault value.  These are DIFFERENT exchanged quantities — reusing one
exchange for both is wrong.

**Trigger:** np>1, bi-material/position-varying material, fault crossing a partition seam,
Part B enabled.

**Actual behavior (as planned):** peer side uses centroid material -> spurious contrast on
shared fault faces with a varying material.

**Expected behavior:** the peer rank evaluates ITS element at the fault-DOF eps-offset
(`x_f - eps*n` on its side) and sends that Zp/Zs; the local rank uses it verbatim.

**Suggested fix (B1):**
```diff
- 3. Shared fault faces: the neighbour element material comes from the cross-rank
-    exchange; reuse / extend `ExchangeBiMaterialNeighbours_` to cover fault shared faces
-    if it does not already.
+ 3. Shared fault faces: the peer rank MUST evaluate its own fault-adjacent element at the
+    fault-DOF eps-offset (`x_f - eps*n` on its side) and send the resulting (Zp,Zs) — NOT
+    its centroid/pool material.  Add a fault-specific exchange (or a new field in
+    `ExchangeBiMaterialNeighbours_`) carrying the per-fault-QP eps-offset impedances; do
+    NOT reuse A2's centroid neighbour-material exchange for the fault per-side read (they
+    are different quantities).  Assert (np=2, depth-varying material, fault on the seam)
+    that the shared-face `Zp_minus` equals the serial eps-offset value to 1e-12.
```

**Test case (intent):**
```cpp
// np=2, depth_profile material, fault cut by the partition near z=-5000:
//   shared-fault DOF Zp_minus(parallel) == Zp_minus(serial eps-offset)  (rel 1e-12)
//   => no spurious across-fault contrast on the seam.
```

---

### [R-102] MODERATE — B1: the eps-offset "map the probe point to the element reference IP" step has no concrete mechanism and no out-of-element guard

**Category:** ASSUMPTION / EDGE_CASE

**Description:**
B1 says "Map `pt_p` to its reference IP `ip_p` in `elem_plus`" but does not say HOW.
The naive route, `ElementTransformation::TransformBack` (physical->reference Newton
solve), can FAIL to converge or return a point OUTSIDE the reference element for a
distorted/sliver near-fault tet, or if `eps` is too large — yielding a garbage material
eval with no error. There is a cleaner, robust mechanism the plan should mandate: the
fault QP already HAS a known element-reference IP in each adjacent element (used to read
`Q_plus`/`Q_minus`); perturb THAT IP along the reference-space inward normal
`J^{-1} n` instead of round-tripping through physical space.

**Trigger:** sliver/distorted near-fault element, or `eps` large enough that
`TransformBack` leaves the reference element.

**Expected behavior:** a deterministic in-element reference IP for the probe, with a guard.

**Suggested fix (B1):**
```diff
- Map `pt_p` to its reference IP `ip_p` in `elem_plus` and `pt_m` to `ip_m` in
-   `elem_minus`, then `material.EvalAt(...)`
+ Obtain `ip_p` by perturbing the fault QP's KNOWN element-reference IP in `elem_plus`
+   along the inward reference normal: `ip_p = ip_face_p + eps_ref * normalize(J_plus^{-1}
+   n_inward)`, with `eps_ref` chosen so `ip_p` stays strictly inside the reference element
+   (clamp; assert `0 <= ip_p coords <= 1` for the simplex).  This uses only the Jacobian
+   at the face QP (no TransformBack Newton solve, no out-of-element risk).  If TransformBack
+   is used instead, MFEM_VERIFY it converged AND the result is in-element; abort loudly
+   otherwise.
```

**Test case (intent):**
```cpp
// distorted near-fault tet: the eps-offset probe IP is in-element (all bary coords in
// [0,1]); EvalAt returns the correct one-sided material (depth-symmetric => both sides
// equal; halfspace => the element's side).
```

---

### [R-103] LOW — C1: "independent of Part B" contradicts the C1 requirement that references B1

**Category:** QUALITY / DEVIATION

**Description:**
C1's Dependencies say "Independent of Part B," but C1 Detailed Requirement 1 requires
"the per-side read (B1) must yield `Zp_plus == Zp_minus` for every TPV31 fault DOF." If
Part B is genuinely not implemented for the TPV31 jobs, there is no per-side read to
constrain; if it IS implemented (unified branch), C1 is not independent of B1. The intent
("TPV31 needs only Part A; and IF B1 is present it must preserve TPV31") should be stated
without the contradiction.

**Suggested fix (C1):**
```diff
- 1. NO Part-B (per-side fault) wiring is needed or expected for TPV31; the per-side read
-    (B1) must yield `Zp_plus == Zp_minus` ...
+ 1. TPV31 needs only Part A.  On the Part-A-only tree the fault setup is unchanged
+    (single-material both sides).  IF this is built on the unified branch (Part B present),
+    B1's eps-offset read MUST yield `Zp_plus == Zp_minus` for every TPV31 fault DOF (the
+    B1 R-001 regression test enforces this), so TPV31 stays on the homogeneous-fault path
+    either way.
- ### Dependencies: depends on Part A (A1-A4).  Independent of Part B.
+ ### Dependencies: depends on Part A (A1-A4).  Does not REQUIRE Part B; if Part B is
+   present it must remain byte-exact for TPV31 (B1 R-001 test).
```

---

### [R-104] LOW — C2: under-specified coordinate frame for the nucleation patch, "which velocity" for stations, and restart of the displacement accumulator

**Category:** EDGE_CASE / QUALITY

**Description:**
Three small TPV6 gaps:
1. The patch is "centred at strike=15000/depth=7500" (SPEC frame, fault 0..30000). The
   canonical TPV config (TPV31) centers the fault at x=0. The plan must say which frame the
   `tpv6.toml` uses, or the implementer may hard-code 15000 in an x=0-centered mesh ->
   patch off the fault center.
2. R-004's per-side velocity/displacement: state that the per-side VELOCITY is the BULK
   particle velocity of the +/- element trace at the station (Q[VX..VZ]), not the
   imposed-state velocity — the spec wants ground motion on each side.
3. The displacement accumulator (R-004) is integrated state; if a TPV6 run checkpoints/
   restarts, the accumulator MUST be checkpointed and restored, else displacement resets
   to 0 on resume. (TPV6 is short — likely single-segment — but state it.)

**Suggested fix (C2):**
```diff
+ Coordinate frame: define the canonical origin in the tpv6 config header (match TPV31's
+ convention).  If the fault is centered at x=0, the patch center is x0=0 (NOT 15000);
+ depth center z0=-7500.  Patch test (R-003) must use the config's frame.
+ Stations: per-side VELOCITY = bulk Q[VX..VZ] of the +/- element at the station (ground
+ motion), not the imposed-state velocity.
+ Restart: include the per-station per-side displacement accumulator in the checkpoint
+ (and zero-init only on a fresh start), so resumed displacement is continuous.
```

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-101 shared-fault eps-offset exchange; R-102 eps-offset IP mechanism)
- Low issues: 2 (R-103 C1 wording; R-104 TPV6 frame/velocity/restart)
- Round-1 fixes: all present, none regressed; R-002 coefficients verified correct.
- Plan compliance: N/A (plan review).
- Verdict: PASS WITH FIXES — the run-arm additions are sound and Arm 1 (upwind+ADER+
  per-side fault) is now VERIFIED viable (ADER fault conversion is per-side at
  `wave_operator.inl:4301-4303`). Two MODERATE gaps remain in the per-side material read:
  the MPI shared-fault path must exchange eps-offset (not centroid) impedances (R-101), and
  the eps-offset IP construction needs a concrete, guarded mechanism (R-102). Fix both
  before B1 implementation; R-103/R-104 are quick clarifications.

## Unreviewed Areas
- Whether `ExchangeBiMaterialNeighbours_` can be extended to carry per-QP fault impedances
  vs needing a new exchange (R-101) — not traced; the implementer must confirm.
- The exact `tpv6` canonical origin convention (R-104.1) — to be fixed in the config.
- Frontera numerics for both arms and TPV31 — not evaluable locally (no full-mesh runs).
