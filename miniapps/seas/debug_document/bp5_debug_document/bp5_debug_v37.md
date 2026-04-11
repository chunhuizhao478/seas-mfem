# BP5 Debug v37: Interior Dirichlet face_int2 Sign Fix

**Date**: 2026-03-16
**Status**: Code applied, all 856 tests pass, TACC runs pending
**Previous**: v36 (per-quad-point traction + cross-element lifts + 0.5 factor)
**Branch**: `feature/elasticity`

---

## 1. Problem Statement

v36 results show p=1 recurrence at ~363 yr (v34: 277 yr, Tandem: 240 yr). The v36 cross-element lifting fix made things WORSE — the fault nearly locked up completely.

## 2. Root Cause: face_int2 Sign Error in Interior Dirichlet BR2 Loading

**File**: `domain/elasticity_operator.hpp`, `AssembleDirichletLoading()`, interior Dirichlet BR2 path

### 2.1 The Bug

The face integral for element 2 used the WRONG sign:

```cpp
// BUGGY (present since v34):
face_int1(u*dim+s, m) += shapes1(m,q) * factor;  // correct: positive
face_int2(u*dim+s, m) -= shapes2(m,q) * factor;  // BUG: should be +=
```

### 2.2 Why This Is Wrong

The bilinear form (`DGElasticityBR2Integrator`) and the slip assembly (`AssembleSlipContributionBR2`) both use the **SAME sign** for both elements' face integrals:

```cpp
// Bilinear form (DGElasticityBR2Integrator, lines 266-270):
IntFace_11(l*dim+s, m) += shapes1(m,q) * factor;  // lift into elem1
IntFace_21(l*dim+s, m) += shapes2(m,q) * factor;  // lift into elem2 (SAME sign)

// Slip assembly (AssembleSlipContributionBR2, lines 1190-1194):
face_int1(u*dim+s, m) += shapes1(m,q) * factor;
face_int2(u*dim+s, m) += shapes2(m,q) * factor;   // SAME sign
```

The physical reason: u_D_int is a prescribed JUMP (u1 - u2 = +/-Vp*t), identical in structure to the fault slip delta_u. The BR2 lifting of the jump into both elements uses the same reference normal. The sign difference for elem2's RHS contribution comes from `elvec2 -= penalty * ...` (opposite penalty sign), NOT from the face integral.

### 2.3 Impact in v34 (no cross-element terms)

In v34, each element only saw its own lift (no cross-element terms):
- elem1: `+penalty * fl_q1(positive)` → correct positive push
- elem2: `-penalty * fl_q2(negative)` → `-neg = positive` → WRONG! Should push elem2 negatively toward -Vp*t/2

The wrong sign on elem2's penalty meant elem2 was being pushed in the OPPOSITE direction from its prescribed value. This weakened the effective tectonic loading, explaining the 277 yr recurrence (vs Tandem's 240 yr).

### 2.4 Impact in v36 (cross-element terms + 0.5 factor)

The v36 fix combined both elements' lifts into a single value:
```
f_lifted_q = 0.5 * (eval1 + eval2)
```

With face_int2 negative:
- eval1 from f_lifted1 (via Minv1 * positive face_int1) → POSITIVE
- eval2 from f_lifted2 (via Minv2 * negative face_int2) → NEGATIVE
- eval1 + eval2 ≈ 0 for symmetric elements!
- f_lifted_q ≈ 0 → NO penalty enforcement!

This completely eliminated the Dirichlet penalty, causing the fault to lock up.

## 3. Fix Applied

Single-character change: `face_int2 -= ...` → `face_int2 += ...`

```cpp
// v37 FIXED:
face_int1(u*dim+s, m) += shapes1(m,q) * factor;
face_int2(u*dim+s, m) += shapes2(m,q) * factor;  // FIXED: same sign as bilinear form
```

Now:
- eval1 and eval2 are both POSITIVE (same sign)
- eval1 + eval2 ≈ 2*eval (constructive addition)
- f_lifted_q = 0.5 * 2*eval = eval
- elem1: +penalty * eval → POSITIVE (drives u1 toward +Vp*t/2) ✓
- elem2: -penalty * eval → NEGATIVE (drives u2 toward -Vp*t/2) ✓

## 4. Expected Impact

| Version | face_int2 sign | Cross-elem | 0.5 factor | Effective penalty | Recurrence |
|---------|---------------|------------|------------|-------------------|------------|
| v34     | WRONG (-)     | No         | No         | ~0.5x correct (elem2 wrong direction) | 277 yr |
| v36     | WRONG (-)     | Yes        | Yes        | ~0x (cancellation) | >363 yr (locked) |
| **v37** | **CORRECT (+)** | **Yes**  | **Yes**    | **1.0x correct** | **< 240 yr?** |

With correct penalty direction for both elements AND correct magnitude (matching bilinear form), the tectonic loading should be strongest in v37. This should produce the shortest recurrence interval, ideally matching or approaching Tandem's 240 yr.

## 5. Verification

### 5.1 Build
Succeeds with no errors (only benign duplicate -rpath linker warnings).

### 5.2 Unit Tests
All 856 tests pass, 0 failures.

## 6. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| v34 | Non-fault Y=0 interior Dirichlet BC (u = Vp*t) | Done (had sign bug) |
| v35 | Per-quad-point traction + BR2 skeleton fix | Reverted |
| v36 | Re-apply v35 + cross-element lifts + 0.5 factor | Done (sign bug amplified) |
| **v37** | **Fix face_int2 sign in interior Dirichlet loading** | **Applied** |
