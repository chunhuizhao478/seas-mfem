# Code Review v5: Diagnose Frontera 4-rank R-101 abort + propose fix

This is a diagnostic review prompted by the dispositive failure of the v4
4-rank sanity run on Frontera (job 7664983).  The run aborted during step 0
of the first RK4 time step with `MFEM abort: R-101 shared-fault DOFData
consistency FAILED` — exactly the guard R-305/R-101 was installed to raise.
This document examines whether the abort is a **real bug** or a **false
positive**, then proposes a fix.

## Summary of the Observed Failure

```
MFEM abort: R-101 shared-fault DOFData consistency FAILED.
Field 'tau1_corr' at centroid (-1.7117e+04, 9.71e-13, -1.5858e+04)
differs by 1.9371e-07 across the two ranks sharing the face
(tol=1.0000e-10).
R-001's (+,-) canonicalisation is not sufficient under the current MFEM
face-normal convention; the fix must be extended (e.g. by having the
owner rank broadcast its DOFData to the non-owner after Evaluate).
```

**Setup**: `tpv102_1000m_p1_1.5s_4rank_dev.sbatch`, 1000 m mesh, 4 MPI ranks,
order 1, tfinal = 1.5 s, `--debug-qnorm` on.  Abort occurs in the `step == 0`
post-Mult call of `VerifySharedFaultDOFDataConsistency` after all four RK4
stages have executed once.

## Review Scope
- **Input evidence**: Frontera stderr (quoted above), job 7664983.
- **Files reviewed**:
  - `miniapps/seas/dynamic/wave_operator.inl`            — R-001 swap + R-101 verifier
  - `miniapps/seas/dynamic/fault_face_flux.cpp`          — `Evaluate` / `ComputeTrialTraction`
  - `miniapps/seas/dynamic/godunov_flux.cpp`             — `BuildFrame` / `BuildRotation(Inverse)`
  - `miniapps/seas/dynamic/fault_face_flux.hpp`          — `DOFData` struct
  - `miniapps/seas/drivers/tpv102_driver.cpp`            — step==0 verify call
  - `miniapps/seas/tests/parallel/test_r101_shared_fault.cpp`  — R-302a inline test
  - `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v3_check.md` (R-001 origin)
  - `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v3_fix.md`   (R-001 fix record)
  - `miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v4_fix.md`   (v4 post-fix)
- **MFEM context**:
  - `mfem/mesh/pmesh.hpp:592-597` — documented shared-face `Elem1No/Elem2No`
    and normal convention (see R-501 root-cause analysis).

## Is this a false or real error?

**This is a REAL bug.**  The abort is not a tolerance artifact, not a
centroid-pairing artifact, and not a diagnostic overreach.  Three pieces
of evidence:

### 1. MFEM's documented face-normal convention forces opposite normals

From `mfem/mesh/pmesh.hpp` (lines 592-597):

> 1) Elem1No - the index of the first element that contains this face
>    this is the element that has the same outward unit normal vector as the face;
> 2) Elem2No - the index of the second element that contains this face
>    this element has outward unit normal vector as the face multiplied with -1;

For a shared face between rank A (owner) and rank B (non-owner):
- On rank A: `Elem1No = local_A`, `Elem2No = ghost_B` → face normal points **from A to B**.
- On rank B: `Elem1No = local_B`, `Elem2No = ghost_A` → face normal points **from B to A**.

These are **opposite**.  `CalcOrtho(ftr->Face->Jacobian(), nor)` in
`ComputeSharedFaceFluxRHS` therefore returns `nor_A = -nor_B` on the two
ranks sharing any shared fault face.

### 2. `BuildFrame(nor)` produces different frames under nor → -nor

`dynamic/godunov_flux.cpp:279-300`:

```cpp
void GodunovFlux::BuildFrame(const real_t *nor, real_t *t1, real_t *t2)
{
   real_t up[3] = {0.0, 0.0, 1.0};
   ...
   // t1 = up x nor
   t1[0] = up[1]*nor[2] - up[2]*nor[1];   // -up·cross vs nor sign
   t1[1] = up[2]*nor[0] - up[0]*nor[2];
   t1[2] = up[0]*nor[1] - up[1]*nor[0];
   ...
   // t2 = nor x t1
}
```

Under `nor → -nor`:
- `t1 → -t1` (cross product is anti-symmetric in nor)
- `t2 = nor × t1 → (-nor) × (-t1) = nor × t1 = t2` (unchanged)

So the fault-local frame `(n, t1, t2)_A = (+n, +t1, +t2)` and
`(n, t1, t2)_B = (-n, -t1, +t2)` differ.  Consequently, the rotation
matrices `Tinv` (global → fault-local) are different on the two ranks.

### 3. The R-001 (+,-) swap is not sufficient

R-001's fix (v3, preserved in v4) says: non-owner calls
`Evaluate(Q_minus_local, Q_plus_local, ...)` instead of
`Evaluate(Q_plus_local, Q_minus_local, ...)`.  This corrects the
"+/- side" labeling but **cannot** correct the rotation-matrix mismatch.

Let `Q(A-side)`, `Q(B-side)` = the same physical fields on each side of the
shared face.  After MFEM ghost exchange (R-005 deep copy, v3 guarded), both
ranks see both sides' raw Q data.  But:
- Rank A's `Q_plus_local = Tinv_A · Q_self_A = Tinv_A · Q(A-side)`.
- Rank A's `Q_minus_local = Tinv_A · Q_nbr_A   = Tinv_A · Q(B-side)`.
- Rank B's `Q_plus_local = Tinv_B · Q_self_B = Tinv_B · Q(B-side)`.
- Rank B's `Q_minus_local = Tinv_B · Q_nbr_B  = Tinv_B · Q(A-side)`.

Owner rank A calls `Evaluate(Tinv_A·Q(A-side), Tinv_A·Q(B-side))`.
Non-owner rank B (with swap) calls
`Evaluate(Tinv_B·Q(A-side), Tinv_B·Q(B-side))`.

Same Q data per slot, but `Tinv_A ≠ Tinv_B`.  Since Evaluate is
a non-linear function of its inputs (friction solve), the outputs differ.
Under `Tinv_B = diag(-1, -1, +1) · Tinv_A` for the velocity/stress components,
the `tau1_trial` computation inside `ComputeTrialTraction` picks up a sign
flip on the `VY, VZ` (and friction-driven `SXY`) components that the R-001
swap does not correct, because Eq. (7b) is:

```
tau1_trial = eta_s · (Q_minus[VY] - Q_plus[VY]
                      + Q_plus[SXY] / Zs_plus
                      + Q_minus[SXY] / Zs_minus)
```

The velocity-jump term `(VY^- - VY^+)` flips sign under the (+,-) swap,
but `Tinv_A · Q → Tinv_B · Q` flips the *other* direction on `VY, SXY`.
The two flips do not compose to the identity in the friction-solve
output.  Result: `tau1_corr` (and every downstream field) drifts between
ranks, by exactly the amount the abort reports.

### 4. Quantitative sanity check

At step 0, `Q = 0` initially, so the first `Mult` produces
`Evaluate(0, 0)` on both ranks — bit-identical output.  But stages 2-4
use `Q_tmp = Q + α · k_previous`, with `k_previous` carrying a small
fault-driven contribution that differs across ranks (see point 3).  The
observed 1.93e-7 Pa drift in `tau1_corr` after step 0 is the *accumulation
of 3 mismatched Evaluate calls*.  It is small at step 0 because initial
slip rate `V_ini ≈ 1e-12 m/s` makes the fault-driven `k_previous` tiny;
as the nucleation ramps V up to O(1 m/s), the drift will grow to
O(stress magnitude), i.e. MPa-scale.  **Ignoring the abort and relaxing
the tolerance is not a solution** — the underlying physics inconsistency
compounds.

### 5. Why did the v3 R-302a inline test pass with `max_diff = 0`?

The inline 2-tet test exercises exactly **one** Mult call with `Q = 0`:

```cpp
// test_r101_shared_fault.cpp:432
Vector Q(wave.Height());
Q = 0.0;
Vector k(Q.Size());
wave.Mult(Q, k);
wave.VerifySharedFaultDOFDataConsistency(1e-10);
```

`Q = 0` ⇒ `Tinv · Q = 0` regardless of frame.  Both ranks call
`Evaluate(fdata, 0, 0)` identically.  Output fields match to machine
precision.  **The test passed only because it happens to dodge the bug,
not because the bug is absent.**  The v3 check flagged R-302a as
"the only way to verify R-001" — but it is not; it verifies a trivial
case where the frame mismatch does not matter.  The driver's actual
multi-stage step 0 is the first run that exposes R-001's insufficiency.

## Findings

### [R-501] [CRITICAL] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — R-001 swap alone cannot produce identical `Evaluate` inputs on two ranks that share a fault face, because MFEM's `GetSharedFaceTransformations` returns opposite face normals and therefore different fault-local frames

**Category:** BUG (fundamental — incomplete R-001 design)

**Description:**
See the root-cause analysis above (§ "Is this a false or real error?").

The R-001 fix (v3) assumed that swapping `(Q_plus_local, Q_minus_local)`
between owner and non-owner produces bit-identical Evaluate inputs.
It does not.  The face-local frame `(nor, t1, t2)` is orientation-
dependent, and on two ranks sharing a face, `nor` is flipped (MFEM
convention) and `t1` therefore flipped too.  The swap corrects the
(+,-) labeling but not the rotation-matrix mismatch.

**Trigger:**
Any TPV102 run on ≥ 4 MPI ranks whose METIS partitioning places ≥ 1
fault face on a partition seam — which, per the Frontera 4-rank log,
happens on the 1000 m mesh.  Probably happens on every production
partitioning (200 m, 400-rank).

**Actual behavior:**
- Both ranks call `Evaluate(fdata, Q_plus_local, Q_minus_local)` with
  inputs that differ by the frame mismatch even after the (+,-) swap.
- DOFData drifts across ranks, step by step, accumulating with every
  RK4 stage.
- On the v4 code: `VerifySharedFaultDOFDataConsistency` at step 0
  aborts with `tau1_corr differs by ~1e-7` on a 1000 m / 4-rank mesh.
  If the verifier were removed (or tolerance loosened), the simulation
  would continue with two ranks' DOFData silently diverging — which
  is the v1 "rupture stops at partition seam" symptom transformed into
  a quieter "rupture fronts are asymmetric across ranks" symptom.

**Expected behavior:**
DOFData on both ranks of a shared fault face must be bit-identical
(or at least within round-off tolerance) after every Evaluate, so
the physical fault state is consistent across the partition seam.
Equivalently: both ranks must call `Evaluate` with identical inputs
AND accumulate identical flux contributions into their local RHS.

**Suggested fix:**
The cleanest fully-correct fix is to **have the owner rank broadcast
post-Evaluate state to the non-owner**, so the two ranks share one
authoritative set of updated DOFData and one authoritative set of
corrected fault-local stresses.  Each rank then locally constructs
its Q_imp in its own frame and accumulates its own flux.

Concretely, in `ComputeSharedFaceFluxRHS`, for each shared fault face QP:

1. **Both ranks** compute `Q_plus_local`, `Q_minus_local` in their own
   frame as today.  The owner (lower rank) calls `Evaluate`, obtaining
   `fdata` (updated) and `Q_imp_plus, Q_imp_minus` (in the owner's
   fault-local frame).  The non-owner does NOT call `Evaluate` yet.

2. Owner packs per-QP into a send buffer: the **frame-invariant**
   DOFData fields (`V1_local, V2_local, tau1_corr_local, tau2_corr_local,
   sigma_n_corr_local, slip_rate, psi, slip1, slip2`) plus the owner's
   `nor, t1, t2` unit vectors (9 doubles) so the non-owner can recover
   them in its own frame.

3. MPI point-to-point: owner sends to peer, non-owner receives.

4. Non-owner unpacks the received state.  To get DOFData in a form
   consistent across the frame difference, note that the *scalar*
   quantities `sigma_n_corr, psi, slip_rate, V_abs` are frame-invariant;
   the *tangent-plane* quantities `V1, V2, slip1, slip2, tau1_corr,
   tau2_corr` are 2-vectors in the face's tangent plane and must be
   rotated from the owner's (t1_A, t2_A) basis to the non-owner's
   (t1_B, t2_B) basis.  The rotation is 2x2 and derivable from the
   exchanged (t1, t2) unit vectors.  Non-owner writes the rotated
   values into its `fdata` and uses them to construct `Q_imp_plus,
   Q_imp_minus` in its own frame (via the Q_imp formulas in lines
   121-144 of `fault_face_flux.cpp`, but using `Q_self`, `Q_nbr` from
   non-owner's perspective).

5. Both ranks accumulate `F_h = flux_.Interior(nor_local, Q_imp_plus_global,
   Q_imp_minus_global)` with their own `nor`, into their own rhs.
   Because the underlying corrected stresses are the same across ranks
   (just rotated into each rank's frame), the flux magnitude is
   consistent and conservation is preserved.

**A simpler alternative that avoids the 2x2 rotation**: exchange the
owner's full `Q_imp_plus, Q_imp_minus` in the *global* frame (9 doubles
each), and have the non-owner use them directly without recomputing.
This costs 18+9 = 27 doubles of MPI exchange per shared fault QP per
Mult, but eliminates the tangent-plane rotation logic entirely.
Pseudocode:

```diff
 // ComputeSharedFaceFluxRHS, inside the shared-fault `if (dof_idx >= 0)` block
-const bool owner = (my_rank_ < peer_rank);
-if (owner)
-{
-   fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
-                         Q_imp_plus, Q_imp_minus);
-}
-else
-{
-   fault_flux_->Evaluate(fdata, Q_minus_local, Q_plus_local,
-                         Q_imp_minus, Q_imp_plus);
-}
+const bool owner = (my_rank_ < peer_rank);
+real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
+if (owner)
+{
+   // Owner computes authoritative state and global-frame imposed Qs.
+   fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
+                         Q_imp_plus, Q_imp_minus);
+   for (int c = 0; c < NUM_STATE; c++)
+   {
+      Q_imp_plus_g[c] = 0.0; Q_imp_minus_g[c] = 0.0;
+      for (int k = 0; k < NUM_STATE; k++)
+      {
+         Q_imp_plus_g[c]  += T(c, k) * Q_imp_plus[k];
+         Q_imp_minus_g[c] += T(c, k) * Q_imp_minus[k];
+      }
+   }
+}
+// Queue a point-to-point exchange (pseudo; real code must batch per-face
+// QPs into a contiguous send buffer and issue after the sf loop).
+ShareOwnerFaultState(peer_rank, owner,
+                     /*inout*/ fdata,
+                     /*inout*/ Q_imp_plus_g, Q_imp_minus_g);
+// Non-owner now has fdata and Q_imp_plus_g/minus_g populated with owner's
+// values.  Flux is computed in global frame directly.
+flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g, F_h);
+// Non-owner flip: its Elem1 is on the canonical-minus side, so the flux's
+// sign relative to its accumulation convention is reversed.
+if (!owner) { for (int c = 0; c < NUM_STATE; c++) { F_h[c] = -F_h[c]; } }
```

(Real implementation must batch all shared fault faces into a single
MPI exchange per Mult, not one per face — performance-critical.  Suggested
structure: before the `for (sf ...)` loop, pre-compute owner/non-owner
flags and send-buffer sizes; run the loop in two passes, where pass 1
computes owner's values and packs send buffers, then a single
MPI_Sendrecv per peer is issued, and pass 2 uses received data to
accumulate fluxes.  See v3 ctor's Allgatherv of fault keys for a nearby
batching template.)

**Test case:**
```cpp
// tests/parallel/test_r101_shared_fault.cpp — new test:
TEST_ASSERT_NO_ABORT(
   "R-501: DOFData identical after 4-stage RK4 on nonzero Q",
   {
      // Same 2-tet inline mesh as R-302a, but drive Mult 4 times on a
      // NON-ZERO Q (populated with a small random or sinusoidal
      // perturbation).  Current code aborts; fixed code reaches the
      // verify call with max_diff < 1e-10 relative.
      Vector Q(wave.Height());
      for (int i = 0; i < Q.Size(); i++) { Q(i) = 1e-3 * std::sin(i); }
      Vector k(Q.Size());
      for (int rk = 0; rk < 4; rk++) { wave.Mult(Q, k); Q.Add(0.25, k); }
      wave.VerifySharedFaultDOFDataConsistency(1e-8);
   });
```

---

### [R-502] [MODERATE] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — the absolute `tol = 1e-10` is unreasonable for fields whose physical magnitudes span 1e-20 to 1e+8; use per-field scale-relative tolerance

**Category:** QUALITY (diagnostic tuning), but blocks running the fix

**Description:**
The verifier uses a single absolute `tol = 1e-10` for all 8 mutable
DOFData fields.  Field magnitudes:

| Field              | Initial | Late-in-run |
|--------------------|---------|-------------|
| `sigma_n_corr`     | 1.2e+08 | 1.2e+08     |
| `tau1_corr`        | 7.5e+07 | 1e+08       |
| `tau2_corr`        | 0       | ~1e+07      |
| `V1`, `V2`         | 1e-12   | ~1 m/s      |
| `psi`              | 0.5     | 0.5 ± 0.1   |
| `slip1`, `slip2`   | 0       | ~m          |

An absolute `1e-10` tol on a field with magnitude `1.2e+08` demands
agreement to roughly `0.83 × machine_epsilon × magnitude` — i.e., less
than 1 ULP of double.  Even a *perfectly correct* fix with sub-ULP
compiler reassociation between ranks would trip this.  Once R-501 is
fixed, the diagnostic should still abort on meaningful physics drift
but tolerate bit-level reassociation noise.

**Trigger:**
Any field of magnitude > 1e-10 whose compiler-reassociated evaluation
differs between ranks by a few ULPs.

**Actual behavior:**
Absolute tol; false positives likely after R-501 fix on any field
with large magnitude.

**Expected behavior:**
Per-field relative tol (e.g., `max(1e-10, |field| × 1e-12)`) or explicit
per-field scales derived from `TPV102Params`.

**Suggested fix:**
```diff
 int fail_local = (max_diff > tol || n_unpaired > 0) ? 1 : 0;
+// R-502: absolute tol unreasonable for fields with magnitudes up to
+// ~1e8 Pa.  Use a scale-relative threshold that tolerates round-off
+// in field values and still catches meaningful drift.  `tol` is now
+// a *relative* tolerance.
+double scale = std::max({std::abs(all_data[max_diff_entry*REC + max_diff_field]),
+                         std::abs(max_diff_entry >= 0 ?
+                                  all_data[(max_diff_entry+1)*REC + max_diff_field]
+                                  : 0.0),
+                         1.0});
+const double effective_tol = tol * scale;
+int fail_local_scaled = (max_diff > effective_tol || n_unpaired > 0) ? 1 : 0;
```

(The exact form depends on how `max_diff_entry` is threaded through;
simpler is to apply the scale-relative compare *inside* the pair loop,
using `max(|all_data[a*REC+k]|, |all_data[b*REC+k]|, 1.0)` as the scale.)

**Test case:**
```cpp
// Synthetic: two ranks with fields (tau=1e8, tau + 1e-7) — should PASS
// with relative tol, FAIL with absolute 1e-10 tol.
TEST_ASSERT(
   "R-502: scale-relative tol tolerates ULP-level noise on large fields",
   ... );
```

---

### [R-503] [MODERATE] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — `same_centroid` uses `scale * DBL_EPSILON` where `scale = max(|va|, |vb|)`; at coordinates near zero (y = 9.7e-13 on the TPV102 fault plane) the tolerance collapses to ~1e-28, tighter than any realistic noise floor

**Category:** BUG (R-404 regression from v4)

**Description:**
R-404 (v4 fix) replaced exact-equality in `same_centroid` with a
scale-relative ULP tolerance:

```cpp
double scale = std::max(std::abs(va), std::abs(vb));
double ulp = scale * std::numeric_limits<double>::epsilon();
if (std::abs(va - vb) > ulp) { return false; }
```

This is correct for coordinates far from zero.  But the Frontera abort
reports `centroid y = 9.7e-13` — the TPV102 fault plane is y = 0, and
MFEM's `ftr->Face->Transform(ip, phys)` produces sub-machine-precision
round-off in the y-coordinate that differs between ranks by a few ULPs.
With `scale ≈ 9.7e-13`, the tolerance is `9.7e-13 × 2.2e-16 ≈ 2e-28`,
which is smaller than the round-off itself.  Two ranks with
`y_A = 9.70e-13` and `y_B = 9.71e-13` (1 ULP of the nonzero value)
would be classified as different centroids → R-305 unpaired abort.

**In the current Frontera abort this did not trigger** because the two
ranks happened to produce bit-identical y values (`ftr->Face->Transform`
is deterministic given same input IP and bit-identical vertex data).
But it is a latent bug: any future MFEM refactor or different partition
topology could produce sub-ULP y-coord drift and trigger a spurious
unpaired abort.

**Trigger:**
Any shared-fault QP whose physical centroid has at least one
coordinate near zero (fault plane), combined with non-bit-identical
inter-rank round-off in that coordinate.

**Actual behavior:**
Tight relative tol at near-zero coord → any noise fails `same_centroid` →
spurious R-305 unpaired abort.

**Expected behavior:**
Use a hybrid tolerance: `max(scale × DBL_EPSILON, absolute_floor)` where
absolute_floor is a mesh-scale-appropriate constant (e.g., `h_min * 1e-10`
or just a fixed `1e-9 m` — smaller than any mesh element but larger than
FP noise at mesh-scale coordinates).

**Suggested fix:**
```diff
 auto same_centroid = [&](int a, int b)
 {
    for (int k = 0; k < 3; k++)
    {
       double va = all_data[a*REC + k], vb = all_data[b*REC + k];
       double scale = std::max(std::abs(va), std::abs(vb));
-      double ulp = scale * std::numeric_limits<double>::epsilon();
-      if (std::abs(va - vb) > ulp) { return false; }
+      // R-503: mesh-scale floor so coordinates near zero (fault plane)
+      // are not compared with an ~1e-28 tolerance.
+      const double abs_floor = 1e-9;   // smaller than any fault-mesh element
+      double tol_k = std::max(scale * std::numeric_limits<double>::epsilon(),
+                              abs_floor * std::numeric_limits<double>::epsilon());
+      if (std::abs(va - vb) > tol_k) { return false; }
    }
    return true;
 };
```

An equally defensible fix: revert to exact `!=` compare (what v3 had
before R-404).  The R-404 motivation (compiler-reassociation round-off)
is hypothetical — the v3 `max_diff=0` test data argues the exact compare
is fine on real builds.  Either fix is acceptable; choose based on
whether the fix agent wants to keep the belt-and-suspenders.

**Test case:**
```cpp
TEST_ASSERT(
   "R-503: same_centroid tolerates ULP noise near zero coordinate",
   {
      // Synthetic centroids with y-coord noise at the sub-ULP-of-1e-13 level.
      // Pre-fix: UNPAIRED (false).  Post-fix: paired (true).
      double a[3] = {1.0, 9.7e-13, 2.0};
      double b[3] = {1.0, 9.7e-13 + 1e-15, 2.0};
      assert(same_centroid(a, b) == true);
   });
```

---

### [R-504] [LOW] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — abort-on-mismatch message blames R-001 specifically; once R-501 is fixed, the message wording will be misleading

**Category:** QUALITY (post-fix doc update)

**Description:**
The abort string in `VerifySharedFaultDOFDataConsistency` currently says:

```
R-001's (+,-) canonicalisation is not sufficient under the current MFEM
face-normal convention; the fix must be extended (e.g. by having the
owner rank broadcast its DOFData to the non-owner after Evaluate).
```

Once R-501 is fixed (owner broadcast implemented), this message should
still be printed on ANY future mismatch — but the user will be confused
because "R-001's (+,-) canonicalisation" is no longer the culprit.  Any
future drift after R-501 would be a new/different bug (e.g., an MPI
exchange miss, a field not broadcast, or a 2x2 rotation bug).

**Suggested fix:** update the message to describe the *symptom* in a
way that survives the R-501 fix:

```diff
-MFEM_ABORT("R-101 shared-fault DOFData consistency FAILED.  "
-           "Field '" << field_name << "' at centroid ("
-           << cx << ", " << cy << ", " << cz
-           << ") differs by " << max_diff
-           << " across the two ranks sharing the face (tol="
-           << tol << ").  R-001's (+,-) canonicalisation is not "
-           "sufficient under the current MFEM face-normal "
-           "convention; the fix must be extended (e.g. by having "
-           "the owner rank broadcast its DOFData to the non-owner "
-           "after Evaluate).");
+MFEM_ABORT("R-101 shared-fault DOFData consistency FAILED.  "
+           "Field '" << field_name << "' at centroid ("
+           << cx << ", " << cy << ", " << cz
+           << ") differs by " << max_diff
+           << " across the two ranks sharing the face (tol="
+           << tol << ").  The two ranks' DOFData diverged — check "
+           "that R-501 owner-broadcast covers every mutable field "
+           "written by FaultFaceFlux::Evaluate and the driver's RK4 "
+           "averaging step, and that the 2x2 tangent-plane rotation "
+           "is applied consistently.");
```

---

## Summary
- Critical issues: **1** (R-501 — R-001 swap is fundamentally insufficient)
- Moderate issues: **2** (R-502 tolerance scale; R-503 R-404 regression at near-zero coords)
- Low issues: **1** (R-504 post-fix message wording)
- Plan compliance: N/A (this is a re-review after a Frontera run surfaced a new bug).
- **Verdict: FAIL — must fix R-501 before proceeding.**
  The v3/v4 fix rounds made the diagnostic audit *better* (it caught the bug)
  but the underlying R-001 design is incomplete.  The Frontera abort is
  a **real, physically meaningful error**, not a tolerance artefact.

## Next Steps for the Fix Agent

1. **Do not** revert R-005 (deep copy), R-002 (global_fault_keys
   Allgatherv), R-101 (verifier), or R-305 (unpair abort).  All of these
   are correct and R-101 is what made the bug visible.
2. **Do** implement R-501 (owner broadcast of post-Evaluate state via
   MPI point-to-point, batched per-Mult).  The simpler global-frame
   Q_imp variant is preferred.
3. **Do** fix R-502 (scale-relative tol) so the fixed code can pass the
   verifier without a false positive on large-magnitude fields.
4. **Do** fix R-503 (hybrid tolerance in `same_centroid`) so coordinates
   near zero (fault plane) don't spuriously unpair.
5. **Do** update R-504 after R-501 lands.
6. **Add** the R-501 test case (nonzero-Q, multi-stage RK4) — the R-302a
   inline test must not be the only shared-fault consistency verification.

## Frontera Run Guidance

Do NOT re-submit any shared-fault-exercising job until R-501 is fixed.
Specifically:
- `tpv102_1000m_p1_1.5s_4rank_dev.sbatch` — WILL abort (just did).
- `tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch` — likely WILL abort.
- `tpv102_200m_p1_1.5s_50rank_dev.sbatch` — WILL abort (more shared seams).
- `tpv102_200m_p1_1.5s_400rank_dev.sbatch` — WILL abort (many more).

A cheap workaround for interim testing: disable the R-101 verifier call
in the driver temporarily (comment out `step == 0` branch in
`tpv102_driver.cpp:817-820`) to let the run proceed — but expect to see
the **v1 rupture-stops-at-partition-seam symptom** reappear, because
the DOFData drift R-101 is catching IS the physical mechanism of the
v1 bug.  That is: do not treat "disable the check" as a fix; treat it
as "prove the check is correctly identifying the v1 bug".

## Unreviewed Areas
- Whether the 2x2 tangent-plane rotation in R-501's fix preserves
  `slip1, slip2` conservation exactly under successive nor swaps.
  (Should — it's a unitary in-plane rotation — but the fix's unit test
  must verify.)
- Whether MFEM's `GetSharedFaceTransformations` gives bit-identical
  `ftr->Face->Jacobian()` entries on both ranks (up to the face-normal
  sign).  The v3 doc implicitly assumes yes; would need a direct test.
- Performance cost of owner-broadcast MPI per Mult at 400 ranks.  For
  TPV102 with ~hundreds of shared fault QPs globally, per-Mult MPI
  traffic is O(kB) which is negligible, but worth benchmarking on
  Frontera once the fix lands.
