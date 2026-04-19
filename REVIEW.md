# Code Review: 2026-04-19 v5 — Frontera 4-rank R-101 abort diagnosis

This review was triggered by the dispositive failure of the v4 4-rank
sanity run on Frontera (job 7664983).  The run aborted during step 0
of the first RK4 time step with `MFEM abort: R-101 shared-fault DOFData
consistency FAILED` — the exact guard R-305/R-101 was installed to raise.

Full analysis is in:
`miniapps/seas/debug_document/tpv102_debug_document/tpv102_debug_v5_check.md`

This file summarises the findings for the /code-fix agent.

## Verdict on the Frontera Abort

**Real bug, not a false positive.**  The R-001 (+,-) canonicalisation
design from v3 cannot produce identical `Evaluate` inputs on two ranks
that share a fault face, because MFEM's `GetSharedFaceTransformations`
returns opposite face normals on the two ranks (documented invariant in
`mfem/mesh/pmesh.hpp:592-597`), and `BuildFrame(nor)` therefore
produces different fault-local frames.  The (+,-) swap corrects the
labeling but not the rotation-matrix mismatch.  The abort correctly
reports real DOFData drift.

## Findings

### [R-501] [CRITICAL] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — R-001 swap alone cannot produce identical Evaluate inputs on two ranks sharing a fault face

**Category:** BUG (fundamental — incomplete R-001 design)

**Description:**
On shared face A↔B, MFEM gives `nor_A = -nor_B` (Elem1's outward normal
convention; pmesh.hpp:592-597).  `BuildFrame` then produces different
`(nor, t1, t2)` frames → different `Tinv` → `Tinv_A·Q ≠ Tinv_B·Q` for
any Q ≠ 0.  The R-001 swap of `(Q_plus_local, Q_minus_local)` on the
non-owner corrects the +/- labeling but cannot correct the frame
mismatch.  Result: `Evaluate(Q_plus_local, Q_minus_local, ...)` is
called with **different numerical inputs on the two ranks**, producing
different outputs and drifting DOFData.

The 4-rank Frontera log reports `tau1_corr` differs by 1.93e-7 Pa at
step 0.  At Q=0 (stage 1) both ranks compute identically; the drift
comes from stages 2-4 which drive Evaluate with Q_tmp ≠ 0.  Drift will
grow with simulation progress as V_abs ramps from 1e-12 to O(1 m/s),
producing O(MPa)-scale divergence by late-in-run.

**Trigger:**
Any TPV102 run whose METIS partitioning places ≥ 1 fault face on a
partition seam (observed on 4-rank 1000 m mesh, expected on all
higher-rank-count production configs).

**Actual behavior:**
Owner and non-owner call `Evaluate` with `Tinv_A·Q_self_A` vs
`Tinv_B·Q_self_B` as the first arg (values agree, transforms differ).
Swapping the arguments on non-owner does not fix this — what the swap
achieves is relabeling `(plus, minus)` within the Evaluate formula,
but the formula's velocity-jump term `(VY^- - VY^+)` is not symmetric
under the combined `(swap, Tinv_A → Tinv_B)` operation because
`Tinv_B = diag(-1,-1,+1) · Tinv_A` on the tangent components.

**Expected behavior:**
DOFData on both ranks of a shared fault face must be bit-identical
(or within true round-off) after every Evaluate.  Equivalently: each
shared fault QP has exactly one authoritative post-Evaluate state,
and both ranks' copies of that state agree.

**Suggested fix — owner-broadcasts, batched per Mult:**

Structure the fix in three phases inside `ComputeSharedFaceFluxRHS`:

1. **Owner Evaluate + pack phase (first loop over shared fault faces):**
   - For each shared fault QP the owner (lower rank ID) rotates its
     local Q into its frame, calls `Evaluate`, and rotates the resulting
     `Q_imp_plus, Q_imp_minus` back to the **global frame** (using T).
   - Owner packs per-QP into a flat send buffer: `fdata` (8 mutable
     scalars: V1_local, V2_local, tau1_corr_local, tau2_corr_local,
     sigma_n_corr, slip_rate, slip1, slip2, psi — per v4 R-301 field
     list) PLUS the full global-frame `Q_imp_plus_g, Q_imp_minus_g`
     (18 doubles per QP).  Record count = 8 + 18 = 26 doubles/QP.
   - Non-owner does not call Evaluate yet; leaves slot for received
     data.

2. **MPI point-to-point exchange (single collective per peer):**
   - Per peer rank, issue `MPI_Sendrecv` exchanging the per-peer
     buffers.  Owner sends, non-owner receives (inverse on the peer).
     Total traffic ~ 26 * n_shared_fault_qps_per_peer doubles.

3. **Accumulate flux phase (second loop):**
   - Both ranks now have authoritative `fdata`, `Q_imp_plus_g,
     Q_imp_minus_g` for every shared fault QP they own.
   - For each QP, compute `F_h = flux_.Interior(nor_local,
     Q_imp_plus_g, Q_imp_minus_g, F_h_out)` using the LOCAL `nor`.
   - **Sign correction:** non-owner's local `nor` is opposite owner's
     canonical normal.  Since `flux_.Interior`'s convention integrates
     `A_n^+ Q^+ + A_n^- Q^-` and both ranks end up with Q_imp buffers
     rotated via owner's T (global frame), the flux direction encoded
     in `Q_imp_plus_g / Q_imp_minus_g` matches owner's nor.  Non-owner
     must negate `F_h` before accumulating into its rhs to account for
     its Elem1 being on the `−canonical_nor` side.

Minimal code shape (pseudo — real implementation must batch and use
`MPI_Sendrecv` rather than the per-face call shown):

```diff
 // Inside ComputeSharedFaceFluxRHS, within the `if (is_fault && ...)` block:
-MFEM_VERIFY(sf < static_cast<int>(shared_face_peer_.size()) &&
-            shared_face_peer_[sf].resolved, "...");
-const int peer_rank = shared_face_peer_[sf].peer_rank;
-MFEM_VERIFY(peer_rank != my_rank_, "...");
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
+MFEM_VERIFY(sf < static_cast<int>(shared_face_peer_.size()) &&
+            shared_face_peer_[sf].resolved, "...");
+const int peer_rank = shared_face_peer_[sf].peer_rank;
+MFEM_VERIFY(peer_rank != my_rank_, "...");
+const bool owner = (my_rank_ < peer_rank);
+real_t Q_imp_plus_g[NUM_STATE] = {0}, Q_imp_minus_g[NUM_STATE] = {0};
+if (owner)
+{
+   // Owner computes authoritative state; rotates Q_imp back to GLOBAL
+   // frame so non-owner can use it verbatim.
+   fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
+                         Q_imp_plus, Q_imp_minus);
+   for (int c = 0; c < NUM_STATE; c++)
+   {
+      for (int k = 0; k < NUM_STATE; k++)
+      {
+         Q_imp_plus_g[c]  += T(c, k) * Q_imp_plus[k];
+         Q_imp_minus_g[c] += T(c, k) * Q_imp_minus[k];
+      }
+   }
+}
+// Owner packs (fdata + Q_imp_plus_g + Q_imp_minus_g) into a batched
+// send buffer keyed by peer_rank.  After the sf loop, one
+// MPI_Sendrecv per peer exchanges all QPs at once.  Non-owner
+// unpacks received data into fdata, Q_imp_plus_g, Q_imp_minus_g.
+//
+// (Implementation detail: add a new helper method
+//  `ExchangeSharedFaultState(peer_packs, ...)` that batches the
+//  per-face packs and runs the MPI collective; call it once between
+//  the two sf loops.)
+
+// Accumulate phase (second loop or inline after exchange):
+flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g, F_h);
+if (!owner)
+{
+   // Non-owner's Elem1 is on the opposite side of canonical nor.
+   // Flux computed with canonical (owner's) nor is inbound to its
+   // Elem1, not outbound; negate for the "rhs[Elem1] -= F_h" convention.
+   for (int c = 0; c < NUM_STATE; c++) { F_h[c] = -F_h[c]; }
+}
```

Because the MPI exchange is "between stages 1 and stage-N of the sf
loop", the simplest correct structure is:

```
// Pass 1: build owner/non-owner work-lists.  Owners pack per-peer send
// buffers (fdata + Q_imp_plus_g + Q_imp_minus_g).  Non-owners record
// their expected slot indices.
//
// Batched MPI: per peer, one MPI_Sendrecv (or Send/Recv pair).
//
// Pass 2: for each sf (owner and non-owner both), compute F_h in
// local frame using GLOBAL Q_imp buffers (authoritative after
// exchange), apply non-owner sign flip, accumulate into rhs.
```

**Test case:**
```cpp
// tests/parallel/test_r101_shared_fault.cpp — add R-501 regression test:
// Uses the same 2-tet inline mesh as R-302a but drives Mult on a
// NON-ZERO initial Q to surface the frame mismatch.  Current code
// aborts; fixed code completes with max_diff < scale × 1e-12.
TEST_R501_DOFData_identical_after_multistage_nonzero_Q:
   Vector Q(wave.Height());
   for (int i = 0; i < Q.Size(); i++) { Q(i) = 1e-3 * std::sin(i); }
   for (int rk = 0; rk < 4; rk++) { wave.Mult(Q, k); Q.Add(0.25, k); }
   wave.VerifySharedFaultDOFDataConsistency(1e-10);
   TEST_ASSERT("R-501: DOFData identical after 4 stages on nonzero Q",
               "reached (no abort)");
```

---

### [R-502] [MODERATE] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — absolute tol = 1e-10 is unreasonable for fields spanning 10^-20 … 10^+8

**Category:** QUALITY (diagnostic tuning — blocks validating the R-501 fix)

**Description:**
Field magnitudes: sigma_n_corr ≈ 1.2e+08 Pa, tau1_corr ≈ 7.5e+07 Pa,
V1 ≈ 1e-12 → 1 m/s, psi ≈ 0.5, slip ≈ 0 → 1 m.  A single absolute tol
of 1e-10 for all fields demands sub-ULP agreement on large-magnitude
stresses.  Even a perfectly correct R-501 fix with mild compiler
reassociation between ranks would trip this.

**Trigger:**
Any field with magnitude > 1e-10 whose inter-rank evaluation differs
by a few ULPs.

**Actual behavior:** absolute tol.

**Expected behavior:** scale-relative tol
(e.g., `max(1e-10, |field_max| × 1e-12)`).

**Suggested fix:**
```diff
 int fail_local = (max_diff > tol || n_unpaired > 0) ? 1 : 0;
+// R-502: scale-relative tol — absolute tol is unreasonable for fields
+// with magnitudes up to ~1e8 Pa.  Compute the scale from the maximum
+// absolute field value in the paired data.
+double max_abs_field = 0.0;
+for (int ii = 0; ii < n_entries; ii++)
+{
+   for (int kk = FIELD_BASE; kk < REC; kk++)
+   {
+      max_abs_field = std::max(max_abs_field, std::abs(all_data[ii*REC + kk]));
+   }
+}
+const double effective_tol = std::max(tol, max_abs_field * tol * 1e2);
+fail_local = (max_diff > effective_tol || n_unpaired > 0) ? 1 : 0;
```

(`tol * 1e2` keeps the tol interpretable: if the caller passes `1e-10`
they get `1e-8` relative on the largest field, which is what's needed
for Pa-scale stresses.)

**Test case:**
```cpp
TEST_R502_scale_relative_tol:
   // Inject a 1e-6 drift on tau1_corr ≈ 1e8 Pa → rel drift 1e-14,
   // should PASS.  Pre-fix: FAIL (abs 1e-6 > abs 1e-10).
```

---

### [R-503] [MODERATE] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — R-404 regression: `same_centroid` uses `scale × DBL_EPSILON` which collapses to ~1e-28 at near-zero coordinates (fault plane y = 9.7e-13)

**Category:** BUG (introduced by R-404 in v4)

**Description:**
R-404 added a sub-ULP tolerance `scale × DBL_EPSILON`.  For the fault
plane coordinate y = 9.7e-13, `scale = 9.7e-13` and `tol_k = 2.15e-28`,
far below any realistic inter-rank round-off floor.  The Frontera
abort shows the pairing worked (both ranks produced bit-identical y
values by pure luck of deterministic `Transform`), but any future
topology change or MFEM refactor could introduce sub-ULP drift in y
and produce a spurious R-305 unpaired abort.

**Trigger:**
Shared fault QPs with a near-zero coordinate component combined with
non-bit-identical inter-rank round-off in that coordinate.

**Actual behavior:** coordinate equality is asserted to sub-sub-ULP
precision near zero → fragile.

**Expected behavior:** hybrid tolerance with an absolute floor matched
to the physical mesh scale.

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
+      // R-503: hybrid scale-relative + absolute-floor tolerance.
+      // Fault-plane coordinates can be O(1e-13) after MFEM's affine
+      // face transform, where a pure scale-ULP tolerance collapses to
+      // ~1e-28 and is tighter than any realistic noise.  Floor at
+      // 1e-9 m (~1e-6 of the smallest plausible mesh element) scaled
+      // by epsilon so two ranks' centroid coordinates compare equal
+      // across near-zero axes.
+      const double abs_floor = 1e-9;
+      double tol_k = std::max(scale, abs_floor)
+                     * std::numeric_limits<double>::epsilon();
+      if (std::abs(va - vb) > tol_k) { return false; }
    }
    return true;
 };
```

**Test case:**
```cpp
TEST_R503_same_centroid_nearzero_coordinate:
   double a[3] = {1.0, 9.7e-13, 2.0};
   double b[3] = {1.0, 9.7e-13 + 1e-15, 2.0};
   assert(same_centroid(a, b) == true);   // pre-fix: false, post-fix: true
```

---

### [R-504] [LOW] `dynamic/wave_operator.inl:VerifySharedFaultDOFDataConsistency` — abort message blames R-001 specifically; post-R-501 wording will be misleading

**Category:** QUALITY (doc)

**Description:**
The abort text `"R-001's (+,-) canonicalisation is not sufficient ..."`
is correct today but will be misleading after R-501 is applied.  Any
future drift is no longer "R-001 insufficient" — it's a new bug (MPI
exchange miss, missing field, or rotation mismatch).

**Suggested fix:**
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
+           "averaging step.");
```

---

## Summary
- Critical issues: **1** (R-501 — R-001 design fundamentally incomplete)
- Moderate issues: **2** (R-502 abs tol scale; R-503 R-404 near-zero coord regression)
- Low issues: **1** (R-504 post-fix message wording)
- Plan compliance: N/A (post-deployment bug).
- **Verdict: FAIL — R-501 must be fixed before any further Frontera run.**
  The v4 code correctly catches the bug via the R-101 verifier; the
  underlying R-001 fix is not algorithmically sufficient for MFEM's
  shared-face normal convention.

## Do-Not-Do Guardrails
- Do NOT revert R-005 (deep copy), R-002 (`global_fault_keys` Allgatherv),
  R-101 (verifier), or R-305 (unpair abort).  These are all correct.
- Do NOT merely loosen the tol to hide the abort — R-502 is a legitimate
  scale-relative tol improvement, but R-501 is the actual physics bug.
  Loosening tol alone would make the simulation continue with silently
  drifting DOFData (= the v1 rupture-stops-at-seam symptom with more
  steps before it manifests).
- Do NOT repeatedly rewrite the R-001 swap logic to try to find a
  local-only algebraic fix.  The symmetry is fundamentally broken by
  MFEM's Elem1-outward-normal convention; no amount of sign-flipping
  inside a single rank's Evaluate call can fix it without knowing
  what the peer rank computed.

## Unreviewed Areas
- R-501 fix performance at 400 ranks (per-Mult MPI traffic for shared
  fault QPs).  Expected negligible (~kB/Mult) but must benchmark.
- Whether `ftr->Face->Jacobian()` is bit-identical on both ranks
  (modulo sign).  Documented invariant but not directly tested.
- Whether the 2x2 tangent-plane rotation alternative to R-501 (keep
  Q_imp in fault-local, rotate the 2-vector tangent fields per rank)
  is tractable.  The global-frame Q_imp approach in the suggested fix
  avoids the 2x2 entirely and is simpler.
