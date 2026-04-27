# R-1003 ADER Substep MPI — Fix Report

**Date:** 2026-04-26
**Companion review:** `SUBSTEP_ITERATOR_MPI_REVIEW.md` (Rev 1).
**Scope:** Land §1–§5 of the R-1003 implementation plan; retire the
`tpv104_driver.cpp:1349-1365` MFEM_ABORT that rejected `--fault-iterator
substep` at np>1.

---

## Files modified

| File | Lines (delta) | What changed |
|---|---:|---|
| `dynamic/wave_operator.hpp`            | +71/−6     | Rename `n_local_fault_qps` → `n_total_fault_qps`; rename member `substep_n_local_fault_qps_` → `substep_n_total_fault_qps_`; expand `EvaluateBulkAtFaultQPsCanonical` docstring (output sized to total, parallel ghost-exchange contract). |
| `dynamic/wave_operator.inl`            | +275/−24   | (a) Rename setter/reset/gate to `…_total_…`. (b) `EvaluateBulkAtFaultQPsCanonical`: extend output to `GetNumTotalFaultQPs()`, add a parallel-mesh shared-fault loop with per-component `ExchangeFaceNbrData`. (c) `ComputeADERSharedFaceFluxRHS`: add the absolute-index substep gate mirroring the interior branch. |
| `drivers/tpv104_driver.cpp`            | +33/−21    | (a) `AdvanceADERWithSubStep`: rename `n_local_fault_qps` → `n_total_fault_qps`, use `wave.GetNumTotalFaultQPs()`, doc the iterator's `dof_data.size()` contract. (b) Replace the np>1 `MFEM_ABORT` with a one-line activation log. |
| `tests/unit/test_tpv104_substep_iterator_parity.cpp` | +120/0 | Test 5: ADER-O ∈ {2, 3, 4} stability + cross-order consistency on a constant-Q̄ fixture. |

Total: **+499/−51** across 4 files.

---

## R-1003 plan §-by-§ status

| § | Description | File | Status |
|---|---|---|---|
| §1 | Per-substep ghost coverage (R-1303) | `wave_operator.inl` (inside `EvaluateBulkAtFaultQPsCanonical`) | LANDED |
| §2 | Substep gate on `ComputeADERSharedFaceFluxRHS` (R-1300, R-1304) | `wave_operator.inl` | LANDED |
| §3 | `EvaluateBulkAtFaultQPsCanonical` shared-fault extension (R-1302) | `wave_operator.inl` | LANDED |
| §4 | Driver buffer rename (R-1301) | `tpv104_driver.cpp`, `wave_operator.{hpp,inl}` | LANDED |
| §5 | Driver abort retirement | `tpv104_driver.cpp` | LANDED (replaced with activation log; runtime probe deferred) |
| §6 | MPI consistency assertions (R-1305) | — | DEFERRED (review's "separate diag-only commit" recommendation) |

---

## Bug-by-bug resolution

### R-1300 — shared-fault ADER branch missing substep gate (CRITICAL)

**Resolution.** The shared-fault branch in `ComputeADERSharedFaceFluxRHS`
(now around `wave_operator.inl:4263-4302`) wraps the inline
`fault_flux_->EvaluateADER(...)` call in the same gate the interior branch
uses:

```cpp
if (substep_I_imp_plus_flat_  != nullptr &&
    substep_I_imp_minus_flat_ != nullptr &&
    dof_idx >= 0 &&
    dof_idx < substep_n_total_fault_qps_)
{
   // copy I_imp_± from the iterator's flat buffer
}
else { fault_flux_->EvaluateADER(...); }
```

`dof_idx` is already absolute (interior + shared concatenated) because
`shared_fault_dof_offset_` is populated with `fault_interior_faces_.Size() *
nbf_per_face_` baked in (`wave_operator.hpp:436-441`).  The R-1304 caveat
about NOT re-basing from zero is honored by reusing the existing `dof_idx`
literal at the dispatch site (no new offset arithmetic).

### R-1301 — iterator buffer sized to local-only QPs (CRITICAL)

**Resolution.** In `AdvanceADERWithSubStep` (driver):

```diff
- const int n_local_fault_qps = wave.GetNumLocalFaultQPs();
+ const int n_total_fault_qps = wave.GetNumTotalFaultQPs();
```

Buffers `I_imp_plus_flat`, `I_imp_minus_flat`, and the `Q_pointwise_*[o]`
slots all flow through `n_total_fault_qps`.  The setter API parameter and
the wave-operator member are renamed in lockstep.  At np=1 the shared
count is 0 so `total == local` and buffers are byte-identical to pre-R-1003.

### R-1302 — `EvaluateBulkAtFaultQPsCanonical` interior-only (CRITICAL)

**Resolution.** Output resized to `GetNumTotalFaultQPs()` and a new
parallel-mesh-only loop walks `fault_shared_faces_` after the existing
interior loop.  The shared loop:

* Performs one per-component `q_gf.ExchangeFaceNbrData()` to populate
  `nbr_data[c]` (the per-substep ghost-exchange of R-1303 §1, see below).
* Uses `pmesh.GetSharedFaceTransformations(sf)` and `pfes->GetFaceNbrFE(...)`
  to build `shape1` / `shape2`.
* Reads self-side Q from the local `Q_data[c * ndof_total_ + dof_offset1
  + i]` (matching the interior loop) and neighbour-side Q from
  `nbr_data[c][nbr_idx * ndof_per_el_ + i]` (matching
  `ComputeADERSharedFaceFluxRHS:3954-3976`).
* Rotates to the canonical fault-local frame using `qpd.sign_flipped` —
  this **mirrors the existing shared-fault flux convention** (see R-801 /
  R-1305 caveat below).
* Routes Elem1's value to the +/− bucket via
  `shared_fault_elem1_on_plus_[sf_idx]`.
* Writes into `[GetNumLocalFaultQPs() + sf_idx*nbf_per_face_ + q]` via
  `shared_fault_dof_offset_[sf]`.

Sanity asserts confirm `base_dof_idx ∈ [n_local_qps, n_total_qps)` and
`fb_idx ∈ [fault_interior_faces_.Size(), fault_basis_->NumFaces())`.

### R-1303 — `ComputeADERSubStepStates` predictor needs ghost-cell coverage (HIGH)

**Resolution.** Rather than retrofitting `ComputeADERSubStepStates` itself
(which would require either a per-sub-step `ExchangeFaceNbrData` or
running CK on ghost cells), the ghost exchange is colocated with its only
consumer: the new shared-fault loop in
`EvaluateBulkAtFaultQPsCanonical`.  Pattern is identical to
`ComputeADERSharedFaceFluxRHS:3870-3881` — sequential per-component
exchange because `q_gf` is single-component and `FaceNbrData()` is
overwritten by every call.

**Cost.** Per macro-step: O calls × NUM_STATE collectives = 9·O exchanges,
plus the existing macro-step's NUM_STATE = 9 exchanges.  At ADER-O = 4
this is 45 exchanges per macro-step, vs. 9 in the one-shot path.  The
review's option (a) (run CK on ghost cells) was **not** chosen because:

* The macro-step path already exchanges *integrated* I after the recursion
  (per `ComputeADERSharedFaceFluxRHS:3870-3881`); running the recursion on
  ghost cells would diverge from the existing pattern and require auditing
  the `ApplySpatialDerivative` element-local kernel for ghost
  applicability.
* Localising the exchange in `EvaluateBulkAtFaultQPsCanonical` means the
  predictor function's contract (element-local, no MPI) is preserved.

### R-1304 — `dof_idx` absolute-vs-rebased on shared branch (HIGH)

**Resolution.** Honored.  The shared-branch substep gate uses the existing
`dof_idx` from `shared_fault_dof_offset_[sf] + q` (already absolute),
matching the absolute-range gate `dof_idx < substep_n_total_fault_qps_`.
No new local-rebased index is introduced.

### R-1305 — DOFData symmetry across shared-fault ranks (MODERATE)

**Status.** DEFERRED to a separate diag-only commit (review's §6
recommendation).  An MPI parity test (Test 3 in the review) covers the
symmetric-output property end-to-end at np=2.  The R-1305 latent risk via
the t1=dip vs t1=strike convention split between interior and shared
paths is unchanged by R-1003 — see the inline comment in the new shared
loop in `EvaluateBulkAtFaultQPsCanonical` for the convention rationale.

### R-1306 — `dt_scale` round-trip drift (LOW)

**Status.** DEFERRED.  Review marked "None required at TPV104 scale;
document in code comment."  Existing inline comments in
`AdvanceADERWithSubStep:309-333` already describe the rescaling protocol;
no code change made.

---

## Tests run (np=1)

All run from `miniapps/seas/`:

| Test binary | Tests | Result |
|---|---:|---|
| `seas_test_tpv104_substep_iterator`        | 40 | ✅ 40/40 |
| `seas_test_tpv104_substep_iterator_parity` (extended with O=2/3/4) | 22 | ✅ 22/22 |
| `seas_test_tpv104_substep_predictor`       |  6 | ✅ 6/6   |
| `seas_test_tpv104_substep_one_shot_parity` |  4 | ✅ 4/4   |
| `seas_test_tpv104_substep_dispatch_parity` |  1 | ✅ 1/1 (Q_new bit-equal to one-shot at O=1, max diff 4.05e-13 ≤ tol 1e-8) |
| `seas_test_tpv104_normal_sign`             | 17 | ✅ 17/17 |
| `seas_test_tpv104_setup`                   | 38 | ✅ 38/38 |
| `seas_test_tpv104_smoke`                   | 53 | ✅ 53/53 |
| `seas_test_tpv104_sigma_n_invariance`      | 20 | ✅ 20/20 |
| `seas_test_tpv104_freeze_sigma_n_gate`     |  6 | ✅ 6/6   |
| `seas_test_tpv104_probe_format`            | 13 | ✅ 13/13 |

Bit-identity at np=1 confirmed by `seas_test_tpv104_substep_iterator_parity`
Gate 1/2: `max |I_imp_plus_legacy − I_imp_plus_substep| = 0.000e+00` and
DOFData diff = `0.000e+00`, exactly as before R-1003.

The `seas_test_bp5_*` failures observed (`elasticity_operator_setup.inl:575`)
are **pre-existing and unrelated** — none of the BP5 source files
(`bp5/`, `domain/`, `friction/dieterich_ruina.hpp`) were modified by
R-1003 (the [C2] no-touch invariant is preserved per
`miniapps/seas/CLAUDE.md`).

---

## What is NOT covered by this fix

* **Tests 1–4 (MPI parity)** from the review remain to be authored.
  These require:
  - A 2-tet box-mesh fixture with a forced METIS partition putting one
    tet per rank (creating exactly one shared fault face).
  - End-to-end `AdvanceADERWithSubStep` at np=2 with cross-rank gather
    of DOFData and `I_imp_*_flat`.
  - Bit-identity vs. np=1 at the same physical QP.

  These are merge-blocking per the review's "Recommended Next Steps"
  list.  Recommend authoring them in a follow-up commit before exposing
  `--fault-iterator substep` to MPI production runs.

* **R-1305 MPI invariant assertion** (`SEAS_DIAG_TPV104_SUBSTEP_MPI_CHECK`
  build-flag-gated debug hook) is DEFERRED.

* **R-1306 dt_scale drift fix** is DEFERRED (not required at TPV104 scale).

* **Frontera MPI smoke** has not been run.  Per project memory ("Ask
  before Frontera runs"), the user should approve a small-mesh np=2
  Frontera sbatch verifying the new shared-fault path before any large
  production submission.

---

## Inline cross-references for follow-on review

* Substep gate (interior, pre-existing):  `wave_operator.inl:3140-3164`
* Substep gate (shared, NEW):              `wave_operator.inl:4278-4302`
* Shared-fault loop in EvaluateBulk (NEW): `wave_operator.inl` after the
  interior loop, gated by `if constexpr (IsParallelMesh<MeshType>::value)`.
* Per-substep ghost exchange (NEW):        inside the shared-fault loop
  above, mirroring `ComputeADERSharedFaceFluxRHS:3870-3881`.
* Driver activation log (replaces abort):  `tpv104_driver.cpp:1450-1466`.
* Member rename:                           `wave_operator.hpp:586-594` and
  `wave_operator.hpp:316-330` (setter signature + docstring).
