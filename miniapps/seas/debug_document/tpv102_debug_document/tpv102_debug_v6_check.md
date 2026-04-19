# Code Review v6: BP5 reuse audit — refactor TPV102 to reuse `fault/fault_basis.hpp` unchanged

User directive:
> *"if you can reuse the functions from BP5 workflow, do not create your own,
> BUT do NOT change any code of BP5 that would break it!  If the code needs
> to adjust somehow, create a new function but follows BP5 logic, later we
> need to wire quasi-dynamic and dynamic together, it is better to reduce
> code duplication."*

This review re-audits the v5 TPV102 code with a strict reuse lens:
**every line of TPV102 fault-frame code that has a BP5 equivalent must be
deleted and replaced with a call into BP5's existing, unmodified public API.**
Where BP5's API does not cover a TPV102-specific need, the new code must
live in the TPV102 namespace (`dynamic/`), not in BP5's (`fault/`), and must
follow BP5's algorithmic pattern rather than inventing a parallel one.

## TL;DR

1. `GodunovFlux::BuildFrame` is a duplicate of a subset of
   `FaultBasis::ComputeOrientedFrame` (BP5, private static helper in
   `fault/fault_basis.hpp:366-439`) **minus** the canonical-frame
   (`ref_normal`) orientation flip.  That missing flip is the root
   cause of R-001/R-101/R-501.
2. `FaultBasis` is already MPI-aware and templated on MeshType.  Its
   public methods `Compute`, `ComputeQPBasis`, `AppendSharedFaces`,
   `ComputeQPBasisShared` cover every frame need TPV102's flux code
   has.  **BP5 does not need to be modified.**
3. All three layers of TPV102-specific shared-fault machinery — R-001
   (+,-) swap, R-101 verifier, R-501 owner broadcast — can be
   **deleted** once `WaveOperator` owns a `FaultBasis` instance and
   reads `(normal, tangent1, tangent2)` from its per-QP table instead
   of calling `GodunovFlux::BuildFrame`.
4. The v5 check's R-602 (docstring) and v6's R-603 (frame-dependent
   broadcast fields) evaporate once canonical frame is adopted via
   `FaultBasis`.
5. What remains genuinely TPV102-specific (and must keep its own code,
   but can still reuse BP5 patterns): `GodunovFlux` (Riemann splits),
   `FaultFaceFlux::Evaluate` (trial-correct Riemann), bulk
   `ExchangeFaceNbrData` deep copy (R-005), and explicit RK4 time
   stepping.

## BP5 public API surface that TPV102 should call, unchanged

From `fault/fault_basis.hpp`:

| BP5 method | Signature | TPV102 call site (proposed) |
|---|---|---|
| `FaultBasis::Compute` | `void Compute(Mesh &mesh, const Array<int> &fault_faces, const Vector &ref_normal, const Vector &up)` (line 74) | `WaveOperator` ctor — build per-face basis for interior fault faces. |
| `FaultBasis::ComputeQPBasis` | `void ComputeQPBasis(Mesh &mesh, const Array<int> &fault_faces, const Vector &ref_normal, const Vector &up, const IntegrationRule &ir)` (line 116) | `WaveOperator` ctor — populate per-QP interior fault frame. |
| `FaultBasis::AppendSharedFaces` | `template <typename MeshType> void AppendSharedFaces(MeshType &mesh, const Array<int> &shared_faces, const Vector &ref_normal, const Vector &up)` (line 181) | `WaveOperator` ctor — append shared-fault per-face basis. |
| `FaultBasis::ComputeQPBasisShared` | `template <typename PMeshType> void ComputeQPBasisShared(PMeshType &mesh, const Array<int> &shared_faces, const Vector &ref_normal, const Vector &up, const IntegrationRule &ir, int interior_face_count)` (line 146) | `WaveOperator` ctor — per-QP shared-fault frame. |
| `FaultBasis::GetBasis(int fi)` | `const FaultBasisData &GetBasis(int fi) const` (line 228) | Flux inner loops — read `qp_data[q].normal / tangent1 / tangent2`. |

**No BP5 file is modified.**  `FaultBasis::ComputeOrientedFrame` stays
private (it already has public entry points that wrap it), and every
existing BP5 caller (`pseas.cpp`, `seas_bp5_full`, etc.) is unaffected.

## Findings

### [R-701] [CRITICAL] `dynamic/godunov_flux.cpp:BuildFrame` — duplicates a subset of `FaultBasis::ComputeOrientedFrame` without the canonical-orientation step; delete it (for fault faces) and call into BP5's `FaultBasis` instead

**Category:** DEVIATION (code duplication of BP5 machinery) + root cause of the v1–v5 shared-fault saga

**Description:**
`GodunovFlux::BuildFrame(nor, t1, t2)` (godunov_flux.cpp:279-300) builds a
rank-local tangent frame from the rank-local `nor` returned by
`CalcOrtho`.  It is a subset of BP5's `FaultBasis::ComputeOrientedFrame`
(`fault/fault_basis.hpp:366-439`) but drops the critical orientation
step:

```cpp
// BP5 lines 377-383:
real_t dot = 0.0;
for (int d = 0; d < dim; d++) { dot += n_raw(d) * ref_normal(d); }
sign_flipped = (dot < 0.0);
if (sign_flipped) { n_raw.Neg(); }
```

That `Neg()` is what makes `(n, t1, t2)` identical on both MPI ranks
sharing a face — because both ranks compare against the same
`ref_normal` and flip their own `n_raw` accordingly.  Without it,
rank A gets `nor = +y`, rank B gets `nor = -y` (MFEM Elem1-outward-
normal convention, `pmesh.hpp:592-597`), and every downstream computation
(Tinv, Q_plus_local, Evaluate output, DOFData) drifts.  This is
precisely the bug that triggered the Frontera abort (job 7664983) and
motivated all three workaround layers (R-001, R-101 verifier, R-501
owner-broadcast).

**Trigger:**
Any ≥ 2-rank TPV102 partitioning where METIS places ≥ 1 fault face on
a partition seam (confirmed on the 1000 m 4-rank run) **or** any
partitioning where two interior fault faces on the same rank have
inconsistent `CalcOrtho` vertex orderings (latent on any non-regular
mesh).

**Actual behavior:**
- Every shared fault QP has ≥ 2 independent tangent frames.
- Every workaround (R-001, R-101, R-501) consumes code and MPI
  bandwidth to patch the drift.
- `FaultBasis` exists and solves this exact problem, is already
  exercised by BP5 in production, but TPV102 doesn't call it.

**Expected behavior:**
TPV102's wave operator should own a `FaultBasis` instance.  In
`ComputeFaceFluxRHS` (fault branch) and `ComputeSharedFaceFluxRHS`
(fault branch), the fault-local tangent frame is read from
`fault_basis_->GetBasis(fi).qp_data[q]` instead of being rebuilt
from `CalcOrtho + BuildFrame` at each inner-loop iteration.  Bulk
(non-fault) Godunov fluxes keep using `GodunovFlux::BuildFrame`
because:

- Bulk interior flux is conservation-safe under rank-local frames
  (`F(L, R, +n) = −F(R, L, −n)` is a Godunov identity — verified in
  v5 fix).
- Only the **fault** Riemann solver (`FaultFaceFlux::Evaluate`) is
  nonlinear, which is where frame mismatch propagates into DOFData
  drift.

**Suggested fix (phased, minimal, no BP5 changes):**

Phase A — `WaveOperator` ctor owns and populates a `FaultBasis`:

```diff
 // dynamic/wave_operator.hpp
+#include "../fault/fault_basis.hpp"
 // ... existing includes ...
 
 template <typename MeshType = Mesh>
 class WaveOperator : public TimeDependentOperator
 {
 public:
    WaveOperator(MeshType &mesh, int order,
                 real_t lambda, real_t mu, real_t rho,
-                const BoundaryConfig &bc);
+                const BoundaryConfig &bc,
+                const Vector &fault_ref_normal = Vector(),
+                const Vector &fault_ref_up     = Vector());
 
 private:
    // ... existing members ...
+   // BP5 fault-frame canonicaliser.  Owned by the wave operator; built
+   // once in the ctor for all fault faces (interior-fault first, then
+   // shared-fault appended per BP5's layout contract).  Per-QP data is
+   // populated for both halves.  Flux routines read
+   // fault_basis_->GetBasis(fi).qp_data[q] instead of calling
+   // GodunovFlux::BuildFrame, so the (n, t1, t2) triple agrees across
+   // ranks sharing any fault face.  Not used on non-fault (bulk)
+   // shared faces, where the Godunov conservation identity makes
+   // rank-local frames safe.
+   std::unique_ptr<FaultBasis> fault_basis_;
+   Vector fault_ref_normal_;
+   Vector fault_ref_up_;
 };
```

```diff
 // dynamic/wave_operator.inl (ctor)
+   // BP5 fault frame (R-701).  Default: TPV102 y=0 fault with z-up.
+   fault_ref_normal_ = fault_ref_normal.Size() > 0 ? fault_ref_normal
+                                                   : Vector({0.0, 1.0, 0.0});
+   fault_ref_up_     = fault_ref_up.Size() > 0     ? fault_ref_up
+                                                   : Vector({0.0, 0.0, 1.0});
+   fault_basis_      = std::make_unique<FaultBasis>();
+   if (fault_interior_faces_.Size() > 0)
+   {
+      fault_basis_->Compute(static_cast<Mesh&>(mesh_),
+                            fault_interior_faces_,
+                            fault_ref_normal_, fault_ref_up_);
+   }
+   // Shared-fault portion (parallel only).  AppendSharedFaces uses the
+   // same ref_normal/up so shared QPs produce an identical frame on
+   // both ranks sharing the face (the v5 R-501 problem).
+   if constexpr (IsParallelMesh<MeshType>::value)
+   {
+#ifdef MFEM_USE_MPI
+      if (fault_shared_faces_.Size() > 0)
+      {
+         fault_basis_->AppendSharedFaces(static_cast<ParMesh&>(mesh_),
+                                         fault_shared_faces_,
+                                         fault_ref_normal_, fault_ref_up_);
+      }
+#endif
+   }
+   // Per-QP data.  ir is the same 2*order integration rule used by
+   // ComputeSharedFaceFluxRHS.
+   // Caveat: BP5's `ComputeQPBasis` expects ONE `IntegrationRule` for
+   // all fault faces.  For a mixed triangle/quad fault this must be
+   // extended — out of scope for TPV102 (all tet fault faces are tri).
+   if (fault_interior_faces_.Size() > 0 && nbf_per_face_ > 0)
+   {
+      const int fi0 = fault_interior_faces_[0];
+      auto *ftr = mesh_.GetInteriorFaceTransformations(fi0);
+      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
+                                               2 * order_);
+      fault_basis_->ComputeQPBasis(static_cast<Mesh&>(mesh_),
+                                   fault_interior_faces_,
+                                   fault_ref_normal_, fault_ref_up_, ir);
+   }
+   if constexpr (IsParallelMesh<MeshType>::value)
+   {
+#ifdef MFEM_USE_MPI
+      if (fault_shared_faces_.Size() > 0 && nbf_per_face_ > 0)
+      {
+         const int sf0 = fault_shared_faces_[0];
+         auto *ftr = static_cast<ParMesh&>(mesh_)
+                     .GetSharedFaceTransformations(sf0);
+         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
+                                                  2 * order_);
+         fault_basis_->ComputeQPBasisShared(static_cast<ParMesh&>(mesh_),
+                                            fault_shared_faces_,
+                                            fault_ref_normal_, fault_ref_up_,
+                                            ir,
+                                            fault_interior_faces_.Size());
+      }
+#endif
+   }
```

Note ordering caveat: `nbf_per_face_` is currently set by
`SetFaultDOFData` *after* the ctor, so the ctor's per-QP-basis
construction either (a) must be deferred to `SetFaultDOFData`, or
(b) nbf_per_face_ must be derivable in the ctor from the face geometry.
Option (b) is cleaner: compute from `IntRules.Get(geom, 2*order).GetNPoints()`
directly in the ctor, mirroring the driver's own computation at
`tpv102_driver.cpp:320-342`.

Phase B — fault branch of `ComputeFaceFluxRHS` reads from `fault_basis_`:

```diff
 // dynamic/wave_operator.inl:ComputeFaceFluxRHS, inside the fault branch
-real_t t1v[3], t2v[3];
-GodunovFlux::BuildFrame(nor, t1v, t2v);
-DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
-GodunovFlux::BuildRotation(nor, t1v, t2v, T);
-GodunovFlux::BuildRotationInverse(nor, t1v, t2v, Tinv);
+// Look up fault face's canonical QP frame (BP5 FaultBasis).
+int fi = ??? /* need a face_idx → fault_basis index map */;
+const FaultBasisQPData &qpd = fault_basis_->GetBasis(fi).qp_data[q];
+real_t t1v[3] = {qpd.tangent1[0], qpd.tangent1[1], qpd.tangent1[2]};
+real_t t2v[3] = {qpd.tangent2[0], qpd.tangent2[1], qpd.tangent2[2]};
+real_t nor_canonical[3] = {qpd.normal[0], qpd.normal[1], qpd.normal[2]};
+DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
+GodunovFlux::BuildRotation(nor_canonical, t1v, t2v, T);
+GodunovFlux::BuildRotationInverse(nor_canonical, t1v, t2v, Tinv);
```

(A side-effect: the final `flux_.Interior(nor_canonical, ...)` — not
the rank-local `nor` — must be used for conservation with the opposite
side.  Same change in the shared-face branch.  The Godunov identity
still holds under a globally-canonical `nor_canonical`; it is the
*local* accumulation sign that we must be careful about — see Phase D.)

Phase C — shared-face flux (`ComputeSharedFaceFluxRHS`) replaces the
v5 R-501 3-phase code with a local loop that uses `fault_basis_` for
the frame.  Owner-broadcast, send/recv buffers, auth_state, and
`sf_is_owner[]` all go away:

```diff
 // dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS
-const int nfs = fault_shared_faces_.Size();
-constexpr int AUTH_REC_SIZE = 9 + 2 * NUM_STATE;
-const bool fault_active = (fault_flux_ && fault_dof_data_ && nfs > 0 && nbf_per_face_ > 0);
-std::vector<bool> sf_is_owner(n_shared, false);
-std::map<int, std::vector<int>> sf_idx_by_peer;
-std::vector<real_t> auth_state;
-if (fault_active)
-{
-   // Phase 1 (owner): Evaluate + pack 27 doubles/QP ...
-   // Phase 2 (MPI): Isend/Irecv/Waitall ...
-   // Phase 3 builds F_h from auth_state ...
-}
-// ...single Phase-3-style face loop already here...
+const int nfs = fault_shared_faces_.Size();
+const bool fault_active = (fault_flux_ && fault_dof_data_ && nfs > 0 && nbf_per_face_ > 0);
+
+// Build sf → fault_basis index map (interior_faces first in the basis,
+// then shared).  One entry per shared face.
+std::vector<int> sf_to_basis_idx(n_shared, -1);
+if (fault_active)
+{
+   for (int sf_idx = 0; sf_idx < nfs; sf_idx++)
+   {
+      int sf = fault_shared_faces_[sf_idx];
+      if (sf >= 0 && sf < n_shared)
+      {
+         sf_to_basis_idx[sf] = fault_interior_faces_.Size() + sf_idx;
+      }
+   }
+}
+
+for (int sf = 0; sf < n_shared; sf++)
+{
+   // ... existing per-sf boilerplate (CalcShape, Q_self/Q_nbr) kept ...
+
+   bool sf_fault = fault_active && (sf < (int)shared_face_bdr_attr_.size())
+                   && (shared_face_bdr_attr_[sf] == bc_.fault_attr)
+                   && bc_.fault_attr > 0 && sf_to_basis_idx[sf] >= 0;
+
+   for (int q = 0; q < ir.GetNPoints(); q++)
+   {
+      // ... existing per-q boilerplate ...
+
+      if (sf_fault)
+      {
+         // Canonical per-QP frame from BP5 FaultBasis.
+         const FaultBasisQPData &qpd =
+            fault_basis_->GetBasis(sf_to_basis_idx[sf]).qp_data[q];
+         real_t can_nor[3] = {qpd.normal[0], qpd.normal[1], qpd.normal[2]};
+         real_t t1v[3]     = {qpd.tangent1[0], qpd.tangent1[1], qpd.tangent1[2]};
+         real_t t2v[3]     = {qpd.tangent2[0], qpd.tangent2[1], qpd.tangent2[2]};
+
+         DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
+         GodunovFlux::BuildRotation(can_nor, t1v, t2v, T);
+         GodunovFlux::BuildRotationInverse(can_nor, t1v, t2v, Tinv);
+
+         real_t Q_plus_local[NUM_STATE], Q_minus_local[NUM_STATE];
+         for (int c = 0; c < NUM_STATE; c++)
+         {
+            Q_plus_local[c] = 0.0; Q_minus_local[c] = 0.0;
+            for (int k = 0; k < NUM_STATE; k++)
+            {
+               Q_plus_local[c]  += Tinv(c, k) * Q_self[k];
+               Q_minus_local[c] += Tinv(c, k) * Q_nbr[k];
+            }
+         }
+
+         int dof_idx = shared_fault_dof_offset_[sf] + q;
+         DOFData &fdata = (*fault_dof_data_)[dof_idx];
+         real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
+         fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
+                               Q_imp_plus, Q_imp_minus);
+
+         real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
+         for (int c = 0; c < NUM_STATE; c++)
+         {
+            Q_imp_plus_g[c] = 0.0; Q_imp_minus_g[c] = 0.0;
+            for (int k = 0; k < NUM_STATE; k++)
+            {
+               Q_imp_plus_g[c]  += T(c, k) * Q_imp_plus[k];
+               Q_imp_minus_g[c] += T(c, k) * Q_imp_minus[k];
+            }
+         }
+         flux_.Interior(can_nor, Q_imp_plus_g, Q_imp_minus_g, F_h);
+      }
+      else
+      {
+         flux_.Interior(nor, Q_self, Q_nbr, F_h);
+      }
+
+      // Accumulate as before.  Sign is correct because can_nor is the
+      // SAME on both ranks and the Godunov conservation identity
+      // F(L, R, +n) = -F(R, L, -n) is no longer needed — both ranks
+      // now use the SAME n and SAME L, R (via swap based on who has
+      // Elem1 on canonical-+ side vs canonical-- side).
+      for (int c = 0; c < NUM_STATE; c++)
+         for (int i = 0; i < ndof; i++)
+            rhs[c*ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h[c];
+   }
+}
```

Phase D — per-rank self/nbr assignment under canonical frame.  When both
ranks use `can_nor` (identical) as the face normal but each rank's
`Elem1` is on a specific side of that normal, Evaluate needs to see
`(Q_plus, Q_minus) = (Q on +can_nor side, Q on -can_nor side)`
regardless of which rank's Elem1 is which side.  The v5 R-001 swap
logic **is still needed**, but now simplified: instead of a `my_rank <
peer_rank` heuristic (which was wrong in the first place — see
`tpv102_debug_v5_check.md`), swap based on `sign_flipped` from
`FaultBasis::ComputeOrientedFrame`:

```cpp
// Whether this rank's Elem1 is on the +can_nor or -can_nor side.
// `sign_flipped` is true iff this rank had to flip its local CalcOrtho
// to match can_nor — equivalently, this rank's Elem1 is on the -can_nor
// side.
const FaultBasisData &bdata = fault_basis_->GetBasis(sf_to_basis_idx[sf]);
const bool my_elem1_is_minus = bdata.qp_data[q].sign_flipped;

// Q_self is local Elem1's Q in canonical frame.  Q_nbr is ghost's Q.
// Evaluate(fdata, Q_plus, Q_minus).  If my Elem1 is on "+" side,
// Q_self is Q_plus; if on "-" side, Q_self is Q_minus.
const real_t *Q_plus_local_canonical  = my_elem1_is_minus ? Q_minus_local : Q_plus_local;
const real_t *Q_minus_local_canonical = my_elem1_is_minus ? Q_plus_local  : Q_minus_local;
// ... and Q_imp output buffers swap too, analogous to v5 R-501 ...
```

Both ranks now call `Evaluate(fdata, Q_plus_local_canonical, Q_minus_local_canonical, ...)`
with **bit-identical arguments** (because can_nor, t1, t2 are bit-identical
and the self/nbr → canonical-+/− mapping is a label change not a frame
change).  Evaluate produces bit-identical output.  DOFData is bit-identical.
**No MPI broadcast needed.**

Phase E — delete dead code:

- `GodunovFlux::BuildFrame` — keep (bulk still uses it).  No change.
- `shared_face_peer_` + `SharedFacePeer` struct — DELETE (no
  owner/non-owner concept anymore).
- R-501 3-phase MPI Isend/Irecv/Waitall block (~270 lines in
  `ComputeSharedFaceFluxRHS`) — DELETE.
- R-101 verifier — KEEP as insurance.  Expect `max_rel_diff=0` always
  under canonical frame.

**Test case:**
Reuse `TestR501_MultiStageRK4NonzeroQ` (nonzero-Q multi-stage).
Expected outcome after R-701: verifier still reports
`max_rel_diff=0, 3 pairs matched, 0 unpaired entries`.
Also: add `TestR601_TwoSharedFaultsPerPeer` (4-tet inline, 2 shared
fault faces between 2 ranks) to catch the R-601 ordering concern from
the previous v6 draft — which under R-701 trivially passes because
there is no packing at all.

---

### [R-702] [MODERATE] `dynamic/wave_operator.hpp:VerifySharedFaultDOFDataConsistency` — docstring says "Absolute tolerance", implementation (v5 R-502) is relative; carried over from v6 previous draft

**Category:** DEVIATION (docstring drift)

**Description:** Same as v6's previous R-602.  After R-701 is applied,
the verifier will pass trivially and the docstring inconsistency is
still a footgun for future callers.

**Suggested fix:**
```diff
-   /// @param[in] tol  Absolute tolerance for per-field equality; default
-   ///                 1e-10 matches the double-precision noise floor of
-   ///                 one Evaluate call.
+   /// @param[in] tol  RELATIVE tolerance for per-field equality:
+   ///                 `|a-b| / max(|a|, |b|, 1.0) > tol` triggers an
+   ///                 abort.  Default 1e-10 tolerates the double-
+   ///                 precision noise floor on Pa-scale fields.  Under
+   ///                 R-701's canonical fault frame, this check is
+   ///                 expected to report `max_rel_diff=0` always; it
+   ///                 remains as regression insurance.
```

**Test case:**
```cpp
TestR702_relative_tol_semantic_contract:
   // Two ranks with tau1_corr = 1e8 and 1e8 + 1e-2 (absolute diff 1e-2,
   // relative diff 1e-10).  Passing tol = 1e-9 must not abort under
   // the relative interpretation.
```

---

### [R-703] [LOW] `dynamic/godunov_flux.cpp:BuildFrame` — `up` direction hardcoded; after R-701 this function is called only on bulk (non-fault) faces where `up` is immaterial, but for future re-use a default-`up` parameter is cleaner

**Category:** QUALITY

**Description:**
Currently `GodunovFlux::BuildFrame` hardcodes `up = {0, 0, 1}` with a
`{1, 0, 0}` fallback.  After R-701, this function is only called for
non-fault bulk shared faces, where the specific tangent orientation
does not matter for conservation.  So this is not a correctness bug.
However, if a future benchmark needs a different `up`, the hardcoding
obstructs re-use.  Keep for now.  Flag only.

**No fix needed** in this round.  Left on the radar.

---

### [R-704] [LOW] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — `Q_self, Q_nbr` computed via shape-function evaluation for shared fault QPs but unused under R-701 (and under v5)

**Category:** QUALITY (perf)

Same as v5 R-604.  After R-701, the shared-fault branch uses
`Q_self, Q_nbr` (to compute `Q_plus_local, Q_minus_local`) — so in
the R-701 code, `Q_self, Q_nbr` are actually used.  **The v5 issue
disappears under R-701.**  No action needed.

---

### [R-705] [LOW] `dynamic/wave_operator.hpp` — `shared_face_peer_` and the `SharedFacePeer` struct become dead code after R-701 Phase E; delete them

**Category:** QUALITY (dead code)

**Description:**
Under R-701 Phase E, no rank is "owner" of a shared fault face any
more — both ranks do the full Evaluate locally, on identical inputs,
producing identical DOFData.  The peer-rank resolution in the ctor
(lines 223-253 of `wave_operator.inl`) and the `SharedFacePeer`
struct in `wave_operator.hpp` are now unused.

**Suggested fix:** delete the struct, the member, and the ctor loop
that populates it.  `shared_mesh_face_set_` (used by
`ComputeFaceFluxRHS` to skip partition-seam faces) stays.

---

## Summary
- Critical issues: **1** (R-701 — 250+ lines of TPV102-local frame
  machinery duplicates BP5 `FaultBasis` and produces the shared-fault
  saga)
- Moderate issues: **1** (R-702 docstring drift)
- Low issues: **3** (R-703 bulk-frame up-axis hardcoding;
  R-704 post-R-701 no-op note; R-705 dead peer-rank struct)
- Plan compliance (against user directive "reuse BP5, don't change
  BP5, add new code only where BP5 doesn't cover"): **INCOMPLETE** —
  current v5 code reuses no BP5 fault-frame machinery; R-701 closes
  this gap entirely without modifying BP5.
- **Verdict: PASS WITH FIXES** — R-701 is the only critical.  All
  MODERATE and LOW findings are small or subsumed by R-701.

## Reuse Matrix

To make the reuse decision explicit for the fix agent:

| TPV102 file | Function / code block | BP5 replacement | Action |
|---|---|---|---|
| `dynamic/godunov_flux.cpp:BuildFrame` | rank-local tangent frame | *(none — keep for bulk non-fault)* | keep |
| `dynamic/wave_operator.inl::ctor` | fault-face classification, nor resolution, peer-rank resolution | `FaultBasis::Compute/ComputeQPBasis/AppendSharedFaces/ComputeQPBasisShared` | **replace**: add `FaultBasis` instance; delete peer-rank loop + `shared_face_peer_` member. |
| `dynamic/wave_operator.inl::ComputeFaceFluxRHS` (fault branch) | `GodunovFlux::BuildFrame(nor, t1v, t2v)` | `fault_basis_->GetBasis(fi).qp_data[q]` | **replace** per-QP call with basis lookup. |
| `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS` (fault branch) | v5 3-phase R-501 | local loop using `fault_basis_->GetBasis(fi).qp_data[q]`; **DELETE** Phase 1 owner-Evaluate-pack, Phase 2 MPI Isend/Irecv/Waitall, Phase 3 auth_state unpack. | **delete + replace**. |
| `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency` | relative-tol diagnostic | *(keep)* | keep; expect trivial pass. |
| `dynamic/godunov_flux.cpp::Interior/Absorbing/FreeSurface` | Godunov Riemann | *(no BP5 analog; BP5 is elliptic)* | keep. |
| `dynamic/fault_face_flux.cpp::Evaluate` | trial-correct Riemann | *(no BP5 analog; BP5 post-processes friction)* | keep. |
| `dynamic/wave_operator.inl::ghost Q exchange (R-005)` | per-component `q_gf.ExchangeFaceNbrData` + deep copy | *(no BP5 analog; BP5 uses HYPRE global solve)* | keep. |

## Action for the /code-fix agent

1. Apply R-701 Phase A–E in one commit.  No BP5 file is to be modified.
   The only new #include inside `dynamic/` is
   `../fault/fault_basis.hpp`.
2. Apply R-702 docstring update.
3. R-703/R-704/R-705 flagged only; no implementation action required.
4. Implement the unit tests in the **Test Plan** below.  Every `MUST-PASS`
   test must green before committing R-701.  `NICE-TO-HAVE` tests can be
   merged in a follow-up if they expose scope creep.

---

## Test Plan (R-701 refactor acceptance criteria)

Because R-701 deletes ~270 lines of MPI exchange + R-001 swap logic and
re-routes fault-frame construction through BP5's `FaultBasis`, the risk
surface is broad: mis-indexing into `qp_data`, wrong interior/shared
offset in `ComputeQPBasisShared`, incorrect `sign_flipped` interpretation
for the (+,-) swap, accumulation-sign drift at the canonical-nor
convention, stale `shared_face_peer_` references after deletion, and
BP5 regressions from unintended header edits.

The tests below are organised into **8 categories** with explicit
severity.  Tests marked `MUST-PASS` block the R-701 merge; `NICE-TO-HAVE`
tests are valuable coverage but may be deferred.  All new tests live in
`miniapps/seas/tests/parallel/test_r701_canonical_frame.cpp` unless
stated otherwise — keeping them in one file eases triage and avoids
fragmenting the shared-fault-test story across binaries.

### Category A — FaultBasis ownership and construction invariants

#### [T-A1] [MUST-PASS] WaveOperator ctor builds FaultBasis on interior-only fault (serial)
**Purpose:** Verify that the refactored ctor populates `fault_basis_`
with correct `num_faces = fault_interior_faces_.Size()` and per-face
tangent vectors matching BP5's `FaultBasis::Compute` semantics.
**Preconditions:** serial build, mesh with ≥ 1 interior fault face.
**Setup:**
```cpp
Mesh mesh = BuildTwoTetSharedFaultMeshInline();  // reuse existing builder
WaveOperator<Mesh> wave(mesh, 1, lambda, mu, rho, bc);
```
**Assertions:**
- `wave.GetFaultBasis() != nullptr`.
- `wave.GetFaultBasis()->NumFaces() == wave.GetFaultInteriorFaces().Size()`.
- For each interior fault face, the normal is unit-length and consistent
  with the mesh geometry up to the canonical flip.
**Accessor needed:** expose `const FaultBasis *GetFaultBasis() const { return fault_basis_.get(); }`
in `wave_operator.hpp`.

#### [T-A2] [MUST-PASS] WaveOperator ctor extends FaultBasis for shared faces (2 ranks)
**Purpose:** Verify that shared faults are appended after interior in
`fault_basis_` following BP5's `(interior, then shared)` layout.
**Preconditions:** 2-rank inline 2-tet mesh, 1 shared fault face, 0
interior fault faces per rank (partition puts each tet on its own rank).
**Setup:**
```cpp
Mesh serial = BuildTwoTetSharedFaultMeshInline();
int part[2] = {0, 1};
ParMesh pmesh(MPI_COMM_WORLD, serial, part);
WaveOperator<ParMesh> wave(pmesh, 1, lambda, mu, rho, bc);
```
**Assertions:**
- `wave.GetFaultBasis()->NumFaces() == wave.GetFaultInteriorFaces().Size() + wave.GetFaultSharedFaces().Size()`
  on each rank.
- Per-QP `qp_data` populated for shared-face half of basis (size >= 1 QP
  per face).

#### [T-A3] [MUST-PASS] Frame agreement across ranks (the whole point of R-701)
**Purpose:** On a shared fault face, both ranks must produce bit-identical
`(normal, tangent1, tangent2)` for each QP.  This is the invariant that
makes R-501 MPI broadcast unnecessary.
**Preconditions:** same 2-rank inline mesh.
**Setup:**
```cpp
const FaultBasis *fb = wave.GetFaultBasis();
int n_int = wave.GetFaultInteriorFaces().Size();   // 0 in inline
// Shared basis lives in index [n_int, n_int + n_shared).
const FaultBasisData &bd = fb->GetBasis(n_int);    // first shared face
const FaultBasisQPData &qp = bd.qp_data[0];        // first QP
// Pack (normal, t1, t2) = 9 doubles into a buffer.
real_t my_frame[9] = {qp.normal[0], qp.normal[1], qp.normal[2],
                      qp.tangent1[0], qp.tangent1[1], qp.tangent1[2],
                      qp.tangent2[0], qp.tangent2[1], qp.tangent2[2]};
real_t peer_frame[9];
MPI_Status status;
MPI_Sendrecv(my_frame, 9, MPI_DOUBLE, 1-rank, 0,
             peer_frame, 9, MPI_DOUBLE, 1-rank, 0,
             MPI_COMM_WORLD, &status);
```
**Assertions:**
```cpp
for (int k = 0; k < 9; k++)
{
   TEST_ASSERT(my_frame[k] == peer_frame[k],
               "frame component " + std::to_string(k) + " differs: "
               + std::to_string(my_frame[k]) + " vs "
               + std::to_string(peer_frame[k]));
}
```
**Critical:** this test fails under the current v5 code (each rank
builds its own frame via GodunovFlux::BuildFrame).  After R-701 it must
pass with bitwise equality (not approximate).

#### [T-A4] [MUST-PASS] `sign_flipped` is opposite on the two ranks sharing a face
**Purpose:** The R-701 self/nbr swap logic (Phase D) relies on
`qp_data[q].sign_flipped == true` identifying the rank whose Elem1 is
on the canonical-minus side.  Exactly one of the two ranks must set
`sign_flipped=true` for any given shared QP.
**Setup:** reuse [T-A3] exchange; compare `qp.sign_flipped`.
**Assertions:**
```cpp
int my_flip = qp.sign_flipped ? 1 : 0;
int peer_flip;
MPI_Sendrecv(&my_flip, 1, MPI_INT, 1-rank, 0,
             &peer_flip, 1, MPI_INT, 1-rank, 0, MPI_COMM_WORLD, &status);
TEST_ASSERT(my_flip + peer_flip == 1,
            "exactly one rank must have sign_flipped=true per shared fault QP");
```

#### [T-A5] [NICE-TO-HAVE] No-fault mesh smoke test
**Purpose:** If `bc_.fault_attr == 0` (no fault), ctor must not attempt
to build `fault_basis_`, and flux routines must not dereference it.
**Setup:** simple box mesh with only absorbing/natural attrs.
**Assertions:** ctor completes; `wave.Mult(Q, k)` runs to completion;
no assertion failures.

---

### Category B — Shared-fault DOFData consistency under canonical frame

#### [T-B1] [MUST-PASS] TestR302_InlineTwoTetSharedFault still passes (Q=0 regression)
**Purpose:** pre-R-701 test must still pass with `max_rel_diff=0`.
**Preconditions:** existing test kept as-is.
**Assertions:** `VerifySharedFaultDOFDataConsistency(1e-10)` does not
abort.  `[R-101 check]` reports `max_rel_diff=0`.

#### [T-B2] [MUST-PASS] TestR501_MultiStageRK4NonzeroQ still passes (nonzero-Q regression)
**Purpose:** pre-R-701 test (which v5 was the first to make pass) must
continue to pass WITHOUT the owner-broadcast MPI code, proving R-701 is
algorithmically sufficient.
**Assertions:** `max_rel_diff=0` after 4 Mult calls on nonzero Q.

#### [T-B3] [MUST-PASS] New: DOFData bit-equality after 10 RK4 STEPS on nonzero Q
**Purpose:** Amortise any second-order drift from the canonical-frame
`Evaluate` path across repeated stages.
**Setup:**
```cpp
// Same 2-tet inline partition as B2; drive 10 full RK4 steps with a
// physically plausible dt (below CFL).
Vector Q(wave.Height());
for (int i = 0; i < Q.Size(); i++) { Q(i) = 1e3 * std::sin(0.37*i + 0.9*rank); }
const real_t dt = 0.5 * wave.ComputeMaxDt(0.5 / (3.0 * (2.0*order + 1.0)));
for (int step = 0; step < 10; step++) { RK4Step(wave, Q, dt); }
wave.VerifySharedFaultDOFDataConsistency(1e-10);   // RELATIVE tol
```
**Assertions:** no abort.

#### [T-B4] [NICE-TO-HAVE] 2-shared-faces-per-peer regression for R-601
**Purpose:** the v6-previous draft flagged R-601 (MPI ordering).  Under
R-701 there is no MPI packing, so the ordering concern is trivially
resolved.  This test proves it by construction.
**Setup:**
```cpp
// Build a 4-tet mesh with TWO shared fault faces between rank 0 and
// rank 1.  Partition 2 tets per rank.  Drive multi-stage Mult on
// nonzero Q.  Without R-701, v5 could silently mis-place data; with
// R-701, per-rank Evaluate on canonical frame produces identical
// DOFData on matching QPs.
Mesh serial = BuildFourTetTwoSharedFaultsMeshInline();  // NEW builder
// ... partition, configure, run ...
wave.VerifySharedFaultDOFDataConsistency(1e-10);
```
**Assertions:** no abort; `n_pairs == 6` (2 faces × 3 QPs).

---

### Category C — Bulk Godunov flux unchanged (regression guard)

#### [T-C1] [MUST-PASS] Serial wave-propagation smoke test
**Purpose:** `seas_test_wave_operator` (serial) must still pass after
R-701.  Proves that `ComputeFaceFluxRHS`'s non-fault branch (which
continues to use `GodunovFlux::BuildFrame`) is untouched.
**Assertions:** 17/17 tests pass.

#### [T-C2] [MUST-PASS] Parallel wave-propagation smoke test (2 and 4 ranks)
**Purpose:** `seas_test_parallel_wave_operator` must pass at 2 and 4
ranks after R-701.  Proves that `ComputeSharedFaceFluxRHS`'s non-fault
branch (unchanged from v5) and the R-005 deep-copy guard still work.
**Assertions:** 5/5 tests pass per rank count.

#### [T-C3] [MUST-PASS] Bulk energy conservation across partition seam
**Purpose:** Numerical energy `E = (1/2) sigma:eps + (1/2) rho |v|^2`
integrated over the full domain must be non-increasing step-to-step
(bulk Godunov is dissipative for Riemann upwind but conserves up to
boundary absorption).  Parallel build only.
**Setup:** reuse `P5_TestParallelEnergyConservation` in
`tests/parallel/test_parallel_wave_operator.cpp`.
**Assertions:** final energy ≤ initial (ratio < 1 + 1e-12).

---

### Category D — Conservation and correctness across shared fault face

#### [T-D1] [MUST-PASS] Flux symmetry: F_owner + F_nonowner = 0 when driven by Q=0
**Purpose:** At Q = 0 with nonzero initial DOFData (tau_ini, V_ini,
psi_ini), the Mult output at shared-fault DOFs must be consistent
across the interface: the sum of per-rank flux contributions at the
shared face's Elem1 on one rank plus Elem1 on the other rank must
equal twice the physical flux (standard DG consistency).
**Setup:** 2-tet inline mesh, 2 ranks, `Q = 0`.
```cpp
Vector k(Q.Size());
wave.Mult(Q, k);                  // k = dQ/dt after one Mult
// At each DOF of rank A's Elem1 (the shared-fault-adjacent element):
//   k_A[Elem1_A, local_dof_i] = accumulated contribution from A's view.
// Similarly on rank B.  The fault flux magnitude at the shared face
// should be equal and opposite across ranks.
// Pack A's Elem1 dof contributions for this QP, Sendrecv, compare.
```
**Assertions:** For the shared-fault flux contribution to dQ/dt at a
QP, `|k_A + k_B| < 1e-10 * max(|k_A|, |k_B|, 1)`.  (Relative tol; the
opposite-sign convention holds per the Godunov conservation identity.)

#### [T-D2] [MUST-PASS] Nucleation at the shared fault centre still triggers equally on both ranks
**Purpose:** At the nucleation zone on TPV102, both ranks' fdata must
see the same `tau1_0 += delta_tau` perturbation and produce identical
initial slip-rate ramp-up.  This exercises `ApplyNucleation` under the
new canonical-frame convention.
**Setup:** 2-tet inline mesh, configure fdata to look TPV102-like with
`tau_ini = 75e6 Pa, V_ini = 1e-12`, run `ApplyNucleation` then one
Mult, then verifier.
**Assertions:** verifier reports `max_rel_diff=0` after nucleation.

---

### Category E — Deletion of dead code

#### [T-E1] [MUST-PASS] `shared_face_peer_` member is deleted / unused
**Purpose:** confirm the R-705 dead code is actually removed.  Scanning
`wave_operator.hpp/inl` for `shared_face_peer_` must return zero hits
after R-701.
**Setup:** `grep -c shared_face_peer_ dynamic/wave_operator.*`.
**Assertions:** returns 0.  (Makefile-level or CI-level check; not a
runtime test.)

#### [T-E2] [MUST-PASS] No `MPI_Isend/Irecv/Waitall` inside `ComputeSharedFaceFluxRHS`
**Purpose:** confirm R-501's batched point-to-point exchange is gone.
**Setup:** `grep -E 'MPI_Isend|MPI_Irecv|MPI_Waitall' dynamic/wave_operator.inl`.
**Assertions:** returns 0 matches (aside from the allowed bulk
`q_gf.ExchangeFaceNbrData` internal calls, which are in MFEM not in
this file).

---

### Category F — BP5 non-regression

#### [T-F1] [MUST-PASS] No BP5 source file modified
**Purpose:** user directive: "do NOT change any code of BP5 that would
break it."  CI-level check.
**Setup:**
```bash
git diff --stat HEAD~1 -- 'miniapps/seas/fault/*' 'miniapps/seas/friction/*' \
   'miniapps/seas/solver/*' 'miniapps/seas/domain/*' \
   'miniapps/seas/config/bp5_params.hpp'
```
**Assertions:** output is empty (no BP5 files changed in the R-701
commit).  Exception: adding a public accessor or a public-friendly
method IS allowed if it's purely additive and doesn't rename/retype
anything that BP5 callers use.  Whoever reviews the R-701 commit must
manually approve any BP5 file change.

#### [T-F2] [MUST-PASS] BP5 test suite passes
**Purpose:** ensure no silent BP5 regression.
**Setup:**
```bash
source /Users/chunhuizhao/miniforge/etc/profile.d/conda.sh
conda activate mfem-dev
cd miniapps/seas
make seas_test_friction_law seas_test_elasticity_operator seas_test_bp5_params
for t in seas_test_friction_law seas_test_elasticity_operator \
         seas_test_bp5_params; do ./$t || exit 1; done
```
**Assertions:** all pass.

---

### Category G — Frontera-shape integration tests

#### [T-G1] [MUST-PASS] 4-rank TPV102 1000 m driver runs to t=0.1 s without R-101 abort
**Purpose:** the ORIGINAL Frontera failure signature was `tpv102_1000m.msh`
at 4 ranks aborting in `VerifySharedFaultDOFDataConsistency` at step 0.
After R-701 on the same config, the run must proceed without abort and
the `[R-101 check]` line must report `max_rel_diff=0`.
**Setup:**
```bash
mpirun -np 4 ./seas_tpv102_driver \
   --mesh tpv102/mesh/tpv102_1000m.msh --mesh-scale 1 \
   --order 1 --cfl 0.5 --tfinal 0.1 \
   --output-dir /tmp/tpv102_t_g1 --output-prefix tg1 --debug-qnorm
```
**Assertions:** exit 0, no `MFEM abort`, log contains
`[R-101 check] shared-fault DOFData consistency OK: ... max_rel_diff=0`.

#### [T-G2] [NICE-TO-HAVE] 2-rank driver 1-second run — station-output parity with serial
**Purpose:** compare TPV102 fault-surface station output at a central
station (`flt_0_3`) between 1-rank and 2-rank runs after R-701.  Should
be identical up to MPI round-off.
**Assertions:** slip, slip_rate, tau at t=1.0 s agree to 1e-10 relative.

---

### Category H — Edge cases that could silently regress

#### [T-H1] [MUST-PASS] `nbf_per_face_` is derivable in ctor before `SetFaultDOFData` is called
**Purpose:** R-701 Phase A moves per-QP FaultBasis construction into
the ctor, but the ctor runs BEFORE the driver calls `SetFaultDOFData`.
The ctor must be able to determine `nbf_per_face_` from face geometry
directly, independent of the driver's subsequent DOFData size choice.
**Setup:** construct `WaveOperator` on a mesh with triangle fault
faces at p=1 (3 QPs at 2*order integration), check
`wave.GetFaultBasis()->GetBasis(0).qp_data.size() == 3`.
**Assertions:** `qp_data.size() == IntRules.Get(Geometry::TRIANGLE, 2).GetNPoints()`.

#### [T-H2] [MUST-PASS] Driver's existing `SetFaultDOFData(data, nqp_per_face)` still validates `nqp_per_face == ctor-computed nbf`
**Purpose:** ensure the driver can't accidentally configure a different
QP count than what FaultBasis was built with.
**Assertions:** if driver passes `nqp_per_face != nbf_per_face_`,
SetFaultDOFData's MFEM_VERIFY aborts with a clear message.

#### [T-H3] [NICE-TO-HAVE] Ctor on a mesh where METIS places all fault interior to one rank
**Purpose:** the original TPV102 1000 m / 2-rank case (observed in
v2/v3 reviews: `rank 0: local=2307 shared=0, rank 1: local=2631 shared=0`).
Under R-701, ctor must successfully build an interior-only FaultBasis
on the rank(s) that have fault; other ranks have `num_faces=0`
FaultBasis.
**Assertions:** no abort; shared-fault code paths are dormant; serial
behaviour reproduced.

#### [T-H4] [MUST-PASS] Serial build path compiles and runs (`!IsParallelMesh<MeshType>`)
**Purpose:** `WaveOperator<Mesh>` (serial instantiation) must still
compile.  `ComputeQPBasisShared` / `AppendSharedFaces` calls are
gated by `if constexpr (IsParallelMesh<MeshType>::value)`.
**Setup:** `seas_test_wave_operator` (serial binary).
**Assertions:** builds and passes.

---

## Test matrix summary

| Category | # MUST-PASS | # NICE-TO-HAVE |
|---|---|---|
| A — FaultBasis ownership | 4 | 1 |
| B — DOFData consistency | 3 | 1 |
| C — Bulk flux regression | 3 | 0 |
| D — Conservation | 2 | 0 |
| E — Dead-code deletion | 2 | 0 |
| F — BP5 non-regression | 2 | 0 |
| G — Frontera-shape integration | 1 | 1 |
| H — Edge cases | 3 | 1 |
| **Total** | **20** | **4** |

The fix agent must make all 20 MUST-PASS tests green before the R-701
commit lands.  The 4 NICE-TO-HAVE tests are strongly recommended for
the same commit but can be deferred to a follow-up if scope is tight.

## Test infrastructure notes

- **`BuildFourTetTwoSharedFaultsMeshInline`** (for T-B4): new helper in
  `tests/parallel/test_r701_canonical_frame.cpp`.  Layout: 4 tets in a
  2×2 configuration where 2 fault-attr interior faces (between the
  "y<0" and "y>0" halves) become shared when partitioned 2 tets per
  rank.  Model after the existing 2-tet inline builder.
- **`TEST_ASSERT`**: reuse the existing macro pattern from
  `test_r101_shared_fault.cpp`.
- **MPI runner**: tests in Category A/B/D/G require `mpirun -np 2` or
  `-np 4`.  CI/local Makefile: `make test-r701-canonical`.
- **BP5 non-regression runner**: Category F tests should be included
  in the default `make test` target so any future BP5-touching refactor
  also exercises them.
- **Accessor additions** (additive to `wave_operator.hpp`, NOT in BP5):
  ```cpp
  const FaultBasis *GetFaultBasis() const { return fault_basis_.get(); }
  ```
  Used by T-A1, T-A2, T-A3, T-A4.  This is the only new TPV102-side
  public method needed by the test plan.

## Coverage gap analysis

What the test plan does NOT cover (flag as future work):

1. **Curved fault surfaces** — BP5's `ref_normal` is a global constant;
   a per-face `ref_normal` would be needed.  Out of scope.  [U-1]
2. **Non-conforming meshes (AMR)** — neither BP5 nor TPV102 supports
   NC fault faces.  Out of scope.  [U-2]
3. **Performance benchmarks** — T-G1/T-G2 are correctness tests, not
   performance regression tests.  Add a perf comparison between v5
   and post-R-701 at 400-rank on Frontera as a follow-up.  [U-3]
4. **Multi-peer (rank A has shared faults with both B and C)** — the
   inline tests only exercise a single peer pair.  At Frontera scale
   this is common; T-G1 exercises it incidentally.  Consider a targeted
   6-tet inline test as future work.  [U-4]
5. **Fault faces with degenerate geometry** (near-coplanar `up` and
   `n_raw`) — BP5's `ComputeOrientedFrame` already `MFEM_VERIFY`s
   `|up × n| > 1e-12`.  TPV102's `up = (0, 0, 1)` with y=0 fault plane
   gives `|up × n| = 1` always.  Safe.  [U-5]

## Unreviewed Areas
- Whether the eventual quasi-dynamic ↔ dynamic coupling will exercise
  non-planar faults, where `ref_normal` must be per-face.  Both BP5
  and TPV102 currently assume planar faults with a single
  `ref_normal`.  Flag for future work.
- Non-conforming mesh refinement on shared fault faces — neither BP5
  nor TPV102 handles NC faces.  Out of scope.
- Performance cost of BP5's `ComputeQPBasis` loop over shared fault
  QPs (done once in ctor, not hot path).  Negligible, not measured.
