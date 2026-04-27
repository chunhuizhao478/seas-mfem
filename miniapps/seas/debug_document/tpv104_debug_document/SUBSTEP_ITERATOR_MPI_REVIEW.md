# TPV104 Sub-step Iterator MPI Review (R-1300+)

## Summary

Adversarial review focused on the **MPI/parallel production limitation** (driver
aborts at np>1 with the R-1003 guard) of the ADER sub-step fault iterator.

Files reviewed:
- `dynamic/tpv104_substep_iterator.{hpp,cpp}` (272 + 688 lines)
- `drivers/tpv104_driver.cpp` — CLI parsing (~L240), R-1003 abort (L1349-1365),
  iterator wiring (L1318-1347, L1815-1825), R-1009 slip gate (L1889-1934),
  R-1008 nucleation gate (L1783-1798), `AdvanceADERWithSubStep` helper (L256-401)
- `dynamic/wave_operator.hpp` (L270-394, L540-555, L617-627)
- `dynamic/wave_operator.inl` — substep dispatch (L3075-3110), shared-fault
  ADER branch (L3924-3998), `EvaluateBulkAtFaultQPsCanonical` (L1422-1545),
  `ComputeADERSubStepStates` (L1043-1118), DOF layout assertions
- Existing tests under `tests/unit/test_tpv104_substep_*`

**Review history**
- **Rev 1** (2026-04-25): 7 new bugs (R-1300..R-1306); R-1003 implementation plan.

---

## Bugs Found — Rev 1

### R-1300: Shared-fault ADER branch never consults the substep side-channel

**Status:** OPEN
**Severity:** CRITICAL (root cause of the np>1 abort)
**File:** `dynamic/wave_operator.inl` lines 3980-3988
(shared-fault `ComputeADERSharedFaceFluxRHS`)

**Description.** The substep side-channel introduced in R-602/R-603
(`substep_I_imp_plus_flat_`, `substep_I_imp_minus_flat_`,
`substep_n_local_fault_qps_`) is consulted only on the *interior-fault* branch
of `ComputeADERFaceFluxRHS` (L3086-3110). The shared-fault branch
(`ComputeADERSharedFaceFluxRHS`, L3924-3988) **always** runs
`fault_flux_->EvaluateADER(...)` inline on the macro-step time-integrated I±,
with no `if (substep_I_imp_plus_flat_ != nullptr)` gate.

This is the literal divergence the driver guard (`tpv104_driver.cpp:1356-1365`)
warns about:

> "The substep guard fires only on the interior-fault branch; shared-fault
> faces would diverge from the iterator's per-sub-step semantic."

At np=1 there are zero shared-fault faces (`fault_shared_faces_.Size() == 0`),
so the bug is dormant. At np>1 every partition seam crossing the fault would
silently use macro-step EvaluateADER while interior-fault QPs use the iterator's
per-sub-step I_imp — producing a non-conservative Riemann solver with two
different time semantics on the same fault.

**Trigger.** Any `mpirun -np >1` with the fault crossing a partition boundary
under `--fault-iterator substep`. Fail-safe is currently the `MFEM_ABORT` at
`tpv104_driver.cpp:1358`.

**Fix.** Mirror the substep gate on the shared-fault branch (see R-1003
Implementation Plan §2). Apply the same `(substep_I_imp_plus_flat_ != nullptr
&& substep_I_imp_minus_flat_ != nullptr && dof_idx_global_in_substep_buffer >=
0 && dof_idx_global_in_substep_buffer < substep_n_local_fault_qps_)` guard
that wraps L3086-3110 of the interior branch, with the index translated to the
shared-fault sub-range (see R-1301).

**Test.** `Test 1` (R-1003-MPI-PARITY).

---

### R-1301: Iterator buffer is sized to local-only QPs but the substep semantic needs all (interior + shared) QPs

**Status:** OPEN
**Severity:** CRITICAL
**File:** `drivers/tpv104_driver.cpp` lines 362, 372-393
(`AdvanceADERWithSubStep` helper)

**Description.** The driver allocates the iterator's I_imp buffers at size
`NUM_STATE * GetNumLocalFaultQPs()` (L362, L372-375):

```cpp
const int n_local_fault_qps = wave.GetNumLocalFaultQPs();   // interior only
const size_t n_words = NUM_STATE * n_local_fault_qps;
std::vector<real_t> I_imp_plus_flat (n_words, 0.0);
std::vector<real_t> I_imp_minus_flat(n_words, 0.0);
```

`GetNumLocalFaultQPs()` is *interior-only* by the comment at
`wave_operator.hpp:388-389`:

```cpp
int GetNumLocalFaultQPs()  const
{ return fault_interior_faces_.Size() * nbf_per_face_; }
int GetNumSharedFaultQPs() const
{ return fault_shared_faces_.Size()   * nbf_per_face_; }
int GetNumTotalFaultQPs()  const { return GetNumLocalFaultQPs() + GetNumSharedFaultQPs(); }
```

The DOF layout (header L381-383) places shared-fault DOFs at offsets
`[GetNumLocalFaultQPs(), GetNumTotalFaultQPs())`. The iterator's `dof_data`
vector is sized `GetNumTotalFaultQPs()` by the driver init code, so when the
iterator runs `Advance(...)` it walks the entire range `[0, dof_data.size())`
and writes shared-fault entries into ψ/slip/V1/V2 — but the **per-sub-step
imposed-state accumulator `I_imp_*_flat` does not have storage for those
indices.** At np>1 this is an out-of-bounds write of NUM_STATE words for every
shared-fault QP per sub-step, into whatever lives past the buffer end.

Worse, even sizing the buffer correctly is not enough: `AdvanceWithSubStepStates`
reads `Q_pointwise_*[o]` of size `NUM_STATE * n` where `n =
dof_data.size()`, but `EvaluateBulkAtFaultQPsCanonical` (L1443) only iterates
`fault_interior_faces_` and outputs at size `NUM_STATE *
GetNumLocalFaultQPs()`. The two contracts are inconsistent the moment any
shared face exists.

**Trigger.** np>1 with the substep flag, after R-1003 lifts the guard. The
`MFEM_VERIFY` inside `AdvanceWithSubStepStates` at `tpv104_substep_iterator.cpp:577-589`
catches the size mismatch first — but only because the verify check is on
expected_words from `dof_data.size()`. The driver's `n_local_fault_qps` makes
the verify condition false → throws → unhelpful runtime error.

**Fix.** Replace `GetNumLocalFaultQPs()` with `GetNumTotalFaultQPs()` in
`tpv104_driver.cpp:362, 372-375, 383-393`, and extend
`EvaluateBulkAtFaultQPsCanonical` to also walk `fault_shared_faces_` and
populate the `[GetNumLocalFaultQPs(), GetNumTotalFaultQPs())` slice (see
R-1303). Also rename `wave_operator.hpp:304` parameter
`n_local_fault_qps` → `n_total_fault_qps` everywhere
(`SetSubStepFaultImposedStates`, `substep_n_local_fault_qps_`) — the *name* is
now actively misleading after the extension.

**Test.** `Test 1`, `Test 4`.

---

### R-1302: `EvaluateBulkAtFaultQPsCanonical` does not visit shared-fault faces

**Status:** OPEN
**Severity:** CRITICAL
**File:** `dynamic/wave_operator.inl` lines 1422-1545

**Description.** The helper iterates only `fault_interior_faces_` (L1443) and
sizes its output at `n_local_qps = GetNumLocalFaultQPs()` (L1433-1437). Under
`--fault-iterator substep` at np>1, every partition-seam fault QP would
appear in the iterator's input as a row of zeros (because the buffer is
zero-initialised), so the friction solve at every shared-fault QP would run
with `Q̃_+ = Q̃_- = 0` — i.e., the friction stage thinks both sides are at
ambient zero stress. This is a silent correctness failure (NaN-free, finite
psi update on the shared rows).

Reading the bulk Q on the neighbour element on a shared face requires
`pmesh->GetFaceNbrFE()` and the `face_nbr_data` ghost layer the driver
populates via `pmesh->ExchangeFaceNbrData()` before the predictor (per
`miniapps/seas/CLAUDE.md`). The helper currently has **no shared-face code
path at all**.

**Trigger.** Once R-1003 lifts the guard and shared faces exist.

**Fix.** Add a second loop over `fault_shared_faces_` after L1543 (see R-1003
Implementation Plan §3) that uses `GetSharedFaceTransformations()` and reads
Elem2 from the face-nbr space, packing into the
`[GetNumLocalFaultQPs(), GetNumTotalFaultQPs())` slice of the output.

**Test.** `Test 1`, `Test 4`.

---

### R-1303: `ComputeADERSubStepStates` predictor is only over locally-owned bulk DOFs — shared-face neighbour Taylor expansion is missing

**Status:** OPEN
**Severity:** HIGH
**File:** `dynamic/wave_operator.inl` line ~1043 (`ComputeADERSubStepStates`)

**Description.** The CK-recursion predictor `ComputeADERSubStepStates` produces
`Q_per_node[o]` of size `NUM_STATE * ndof_total_` — i.e., bulk Q on locally
owned elements. The shared-fault branch needs **both sides' Q at the
sub-step nodes**. On rank A the neighbour element of a shared fault face lives
in `pmesh->face_nbr_data` (ghost layer); the predictor must be run there as
well, OR the per-sub-step Q must be exchanged from the owning rank.

Two viable options:

(a) Each rank runs `ComputeADERSubStepStates` for its locally owned elements
    AND for its face-neighbour ghost layer (requires CK recursion to be
    well-defined on ghost cells — true since the recursion is element-local
    and uses element-local stiffness/mass matrices that are already built on
    ghost cells by MFEM's standard parallel assembly).

(b) Each rank runs the predictor only on owned cells, then exchanges
    `Q_per_node[o]` ghost-cell rows via a per-substep `ExchangeFaceNbrData()`
    call.

Option (a) is cheaper (no extra MPI per sub-step) and matches the existing
pattern that `ComputeADERTimeIntegrated`/`AdvanceADER` use the same I-vector
on ghost cells. Confirm via reading the existing time-integrated path's
ghost handling — if the existing macro-step path already runs the recursion
on ghost cells, (a) is the natural extension. If it relies on
`ExchangeFaceNbrData` of the time-integrated I, then (b) is needed.

**Trigger.** Pre-condition for shared-fault iterator dispatch.

**Fix.** See R-1003 Implementation Plan §1.

**Test.** `Test 2` (predictor parity np=1 vs np=2).

---

### R-1304: `dof_idx` indexing into the substep buffer assumes interior-only ordering

**Status:** OPEN
**Severity:** HIGH
**File:** `dynamic/wave_operator.inl` lines 3086-3094 (interior dispatch site)
                                       L3924-3940 (shared dispatch site, when
                                       R-1300 is fixed)

**Description.** On the interior branch, `dof_idx` is computed via
`fault_face_dof_offset_[f] + q` (~L2486-2487 for the average path, similar for
the ADER path). This offset map is filled for both interior and shared
faces (see L422 + L439 in `wave_operator.hpp`):

```cpp
shared_fault_dof_offset_[fault_shared_faces_[i]] =
   fault_interior_faces_.Size() * nbf_per_face_ + i * nbf_per_face_;
```

So shared-fault `dof_idx` already lives in the
`[GetNumLocalFaultQPs(), GetNumTotalFaultQPs())` range, and a single flat
index works across both branches **provided the iterator's I_imp buffer is
sized `GetNumTotalFaultQPs() * NUM_STATE` (R-1301).** Verify — when adding the
substep gate to the shared branch, do NOT compute a separate
`local_dof_idx = sf_idx_in_fault * nbf_per_face_ + q` and re-index from zero;
use the existing `dof_idx` directly so it lands in the upper slice. Easy
off-by-one trap: the code at L3933-3934 already computes
`sf_idx_in_fault = basis_idx - fault_interior_faces_.Size()` for the
shared-fault tables, but `dof_idx` itself is from `shared_fault_dof_offset_`
(L3926-3928) and is already global (interior + shared concatenated).

**Trigger.** R-1003 implementation; easy to write the wrong index.

**Fix.** When wiring the substep gate on the shared branch (R-1300), reuse the
existing `dof_idx` from L3927-3928 directly:

```cpp
const int dof_idx_global = it->second + q;   // already includes
                                             // GetNumLocalFaultQPs() offset
if (substep_I_imp_plus_flat_ != nullptr &&
    substep_I_imp_minus_flat_ != nullptr &&
    dof_idx_global >= 0 &&
    dof_idx_global < substep_n_local_fault_qps_)   // R-1301 rename
{
   const real_t *src_p = substep_I_imp_plus_flat_  + dof_idx_global * NUM_STATE;
   const real_t *src_m = substep_I_imp_minus_flat_ + dof_idx_global * NUM_STATE;
   for (int c = 0; c < NUM_STATE; c++)
   {
      I_imp_plus[c]  = src_p[c];
      I_imp_minus[c] = src_m[c];
   }
}
else { /* inline EvaluateADER as today */ }
```

**Test.** `Test 1`.

---

### R-1305: Iterator's per-sub-step state writes (slip1/slip2/psi) on shared-fault DOFs require MPI consistency, currently absent

**Status:** OPEN
**Severity:** MODERATE
**File:** `tpv104_substep_iterator.cpp` lines 410-411, 651-652, 413-419, 654-660,
          487, 677

**Description.** In the np=1 path the iterator's per-sub-step writes
`d.slip1 += s.V1 * dt_sub`, `d.psi = UpdateStateAnalyticSlipLawSRW(...)`, and
the final `WriteBackState` are all unambiguous because each fault QP is owned
by exactly one rank's DOFData entry. At np>1, **a shared-fault QP appears as a
DOFData row on BOTH ranks adjacent to the seam** (one as the
canonical-+ side, one as canonical-−). The current driver init populates
DOFData on both ranks identically.

If the iterator runs on both ranks independently with the same canonical
inputs Q̃_±, both ranks land the same final ψ, slip1, slip2 — the SRW analytic
update is deterministic in its inputs, the friction solve is deterministic in
its inputs. **But only if the inputs match across ranks bit-for-bit at every
sub-step.** The shared-fault Q̃_± values come from rank A's locally-owned
Elem1 plus rank A's ghost-layer Elem2; on rank B it's the mirror. Without an
MPI consistency check, a subtle mismatch (e.g., a missed ghost exchange after
nucleation, a sign-frame inconsistency on one rank only — see CLAUDE.md L13
on the FaultBasis t1=dip/t2=strike convention) silently produces ψ, slip1,
slip2 drift between the two ranks for the same physical QP. The resulting
DOFData asymmetry then feeds a wrong canonical-± dispatch on the next
macro-step.

**Trigger.** Any cross-rank inconsistency in the shared-fault frame, ghost
exchange order, or floating-point reduction. `tpv104_review_round*` notes
already discuss the BP5 vs Godunov frame mismatch (CLAUDE.md L13 R-801) — the
substep iterator inherits that risk.

**Fix.** Add a debug-mode invariant assertion: after each sub-step, on the
shared-fault QPs, MPI_Allreduce the per-rank `(d.psi, d.slip1, d.slip2,
d.slip_rate)` for the same physical QP and verify bit-identity. Fail-loud if
the two ranks disagree. Run under `MFEM_DEBUG` only — production cost is
prohibitive.

Suggested guard: `#ifdef SEAS_DIAG_TPV104_SUBSTEP_MPI_CHECK` block in the
iterator's per-sub-step loop (or a post-loop hook in the driver), iterating
over the `[GetNumLocalFaultQPs(), GetNumTotalFaultQPs())` DOFData range and
asserting cross-rank match via `MPI_Allreduce(MPI_MAX) == MPI_Allreduce(MPI_MIN)`
on each scalar.

**Test.** `Test 3` (cross-rank symmetry).

---

### R-1306: `dt_scale` re-`SetSubSteps` per macro-step has cumulative drift via repeated ratio→size→ratio round-trips

**Status:** OPEN
**Severity:** LOW (np=1 path; revisit at np>1)
**File:** `drivers/tpv104_driver.cpp` lines 314-333 (`AdvanceADERWithSubStep`)

**Description.** `AdvanceADERWithSubStep` reads the configured
`deltaT/weights`, computes `dt_scale = dt_step / configured_sum`, builds a
fresh `deltaT_scaled[o] = configured_deltaT[o] * dt_scale`, and re-calls
`iterator.SetSubSteps(deltaT_scaled, configured_weights)`. After the call,
`configured_deltaT` for the NEXT macro-step is the just-scaled vector (the
internal `deltaT_` is now the scaled one) — so the next iteration's
`configured_sum` = previous `dt_step`, not the original seed `dt`.

This works because the formula re-derives `deltaT_scaled[o] =
prev_deltaT[o] * (dt_new / prev_dt_sum)` which is mathematically
`(seed_deltaT[o] / seed_sum) * dt_new`, and the **ratios** are preserved.
But each macro-step accumulates one floating-point divide error of size
~eps. Over `N` macro-steps the relative drift is bounded by
`N * eps`, which is harmless at typical `N ~ 10^5`. The Σ deltaT == dt_macro
verify in `Advance` (L253-261) uses `sum_tol = max(1e-12, 10*O*eps)` so it
won't trip.

**Trigger.** Long simulations (N ≥ 10^11 macro-steps) — outside TPV104 scope.

**Fix.** Cache the seed `(deltaT, weights)` on the driver side once at L1325-1338,
and pass `seed_deltaT * (dt_step / seed_sum)` directly into the iterator
each macro-step rather than reading the iterator's mutated state. Or, add an
overload `iterator.AdvanceWithSubStepStatesScaled(scale=dt_step/seed_sum, ...)`
that scales internally per call without mutating `deltaT_`.

**Test.** None required at TPV104 scale; document in code comment.

---

## Pass-2 verification of np=1 path (no new bugs)

Re-checked the listed concerns:

- **Sub-step counter consistency.** `t_sub_cursor` (iterator L281) is updated
  exactly once per sub-step (L491, L681). The for-loop runs `o ∈ [0, O)` with
  `last_sub_step = (o == O - 1)` — only the last iteration calls
  `WriteBackState` (L487, L677). No stale counter.

- **Predictor / corrector residual.** The iterator commits the **last**
  sub-step's V1/V2/τ\*_corr/σ_n_corr to DOFData via `WriteBackState` (L487,
  L677), and the accumulated `I_imp_*_flat` is the time-weighted sum across
  all sub-steps. The corrector branch in `wave_operator.inl:3086-3110` then
  consumes `I_imp_*_flat` (not the per-sub-step traction) for the bulk RHS.
  The macro-step's bulk-side state matches the iterator's accumulator; the
  fault-side DOFData captures the last-sub-step terminal value, consistent
  with the WriteBack contract documented in the header (L91-94). No off-by-one.

- **R-1009 cross-check.** Driver L1894-1915 explicitly gates
  `dof_data[i].slip1/2 += V·dt_step` on `!use_substep_iterator`; iterator
  L410-411 / 651-652 owns the per-sub-step accumulation. One-shot path keeps
  the driver-level integration. R-1009 fix is correct.

- **Predictor cache invalidation.** The iterator does not cache predictor
  state across macro-steps; `Advance` zeros `I_imp_*_flat` at L264-266 / 611
  on every entry. No stale state if `dt` is rejected (TPV104 uses fixed dt;
  RK4 path doesn't reject).

- **ADER-O dependence.** `Advance` accepts arbitrary `O = deltaT_.size()`,
  and the driver wires equal-width sub-steps `1/O` for `O ∈ [1, ader_order]`.
  Only `ader_order = 2` is exercised in tests
  (`test_tpv104_substep_*_parity.cpp` line 192). Tests at O=3, O=4 are missing
  — see Test 5.

- **NumSubSteps equals ADER order?** Driver L1326-1337 sets
  `O = max(1, ader_order)`. Confirmed.

---

## Proposed Unit Tests

### Test Suite — MPI-aware substep iterator

**Test 1 (R-1003-MPI-PARITY): SubstepIterator np=2 vs np=1 bit-identity**
- Target: `AdvanceADERWithSubStep` end-to-end via `tpv104_driver`
- File: `tests/unit/test_tpv104_substep_iterator_mpi.cpp` (new)
- Validates: For a 2-element box mesh with one shared-fault face (METIS
  partition forced via `ParMesh::Partition`), `I_imp_*_flat[dof_idx]` and
  `dof_data[dof_idx].(psi, slip1, slip2, V1, V2, tau*_corr)` for the shared QP
  must match the np=1 result bit-for-bit at the same physical QP across all
  three ADER orders {2, 3, 4}.
- Priority: HIGH (blocks R-1003 merge)
- Detects: R-1300, R-1301, R-1302, R-1304
- Skeleton:
  ```cpp
  // Build 2-tet box, fault face = the interior face shared between the
  // two tets. Place METIS to put each tet on a different rank.
  ParMesh pmesh = MakeTwoTetFaultMesh(MPI_COMM_WORLD);  // |->np=2 splits
  Tpv104SubStepIterator iter(flux, srw); iter.SetSubSteps(...);
  Vector Q_init = ...;
  for (int O : {2,3,4}) {
     // np=2 reference run
     std::vector<real_t> I_imp_plus_p2, I_imp_minus_p2;
     std::vector<DOFData> dof_p2 = run_one_step(pmesh, /*np=2 active*/, O,
                                                I_imp_plus_p2, I_imp_minus_p2);
     // np=1 reference run with the same physical QP coordinates
     auto [I_p1, dof_p1] = run_one_step_np1(O);
     // Cross-rank gather of np=2 shared-fault DOFData onto rank 0
     // and bit-compare per physical QP.
     EXPECT_BITEQ(dof_p2_shared_qp.psi,   dof_p1_for_same_xyz.psi);
     EXPECT_BITEQ(I_imp_plus_p2[dof_idx], I_p1[dof_idx_p1_for_same_xyz]);
  }
  ```

**Test 2: ComputeADERSubStepStates ghost-cell parity (np=2 vs np=1)**
- Target: `WaveOperator::ComputeADERSubStepStates`
- File: `tests/unit/test_tpv104_compute_ader_substep_states_mpi.cpp` (new)
- Validates: After R-1003-§1, on a 2-element mesh the per-sub-step
  `Q_per_node[o]` evaluated on rank-A's owned cell for the cell on rank B
  matches what rank B computes locally for the same cell, bit-for-bit, in all
  NUM_STATE channels at all `tau_nodes`.
- Priority: HIGH
- Detects: R-1303
- Skeleton: pseudo:
  ```cpp
  pmesh = TwoTetMesh(MPI_COMM_WORLD);
  Q_init = FillSinusoidal();
  std::vector<Vector> Qpn; wave.ComputeADERSubStepStates(Q_init, dt, 4,
                                                        tau_nodes, Qpn);
  // Exchange Qpn[o] ghost-layer to compare cells across ranks.
  // Compare against np=1 reference for the same cell.
  ```

**Test 3 (R-1305): per-sub-step DOFData symmetry across shared-fault ranks**
- Target: `Tpv104SubStepIterator::Advance{,WithSubStepStates}` invoked under
  the driver's MPI dispatch
- File: `tests/unit/test_tpv104_substep_iterator_shared_face_symmetry.cpp` (new)
- Validates: For every shared-fault QP and every sub-step `o`, after the
  iterator returns, `(d.psi, d.slip1, d.slip2, d.V_abs)` on rank A equals
  rank B (the canonical-+/− pair for the same physical QP), bit-for-bit.
- Priority: HIGH
- Detects: R-1305 (mismatch in inputs would surface here as DOFData drift)
- Skeleton:
  ```cpp
  iter.AdvanceWithSubStepStates(...);
  // Gather per-rank DOFData[i] for shared-fault i to rank 0
  // and assert per-physical-QP equality.
  for (int i = wave.GetNumLocalFaultQPs(); i < wave.GetNumTotalFaultQPs(); i++) {
     real_t my_psi = dof_data[i].psi;
     real_t max_psi, min_psi;
     MPI_Allreduce(&my_psi, &max_psi, 1, MPI_DOUBLE, MPI_MAX, comm);
     MPI_Allreduce(&my_psi, &min_psi, 1, MPI_DOUBLE, MPI_MIN, comm);
     EXPECT_EQ(max_psi, min_psi)
       << "shared-fault DOFData drift at QP " << i;
  }
  ```

**Test 4: np=2 substep at O=1 == np=2 one-shot at O=1 (T_TPV104_SSI_3 lifted to MPI)**
- Target: substep dispatch with `--fault-iterator substep` + `--ader-order 1`
  (or O=1 quadrature) vs the default one-shot path
- File: `tests/unit/test_tpv104_substep_one_shot_parity_mpi.cpp` (new)
- Validates: Same contract as the existing
  `test_tpv104_substep_one_shot_parity.cpp` (which is np=1 only) but at np=2
  with shared-fault faces.
- Priority: HIGH
- Detects: R-1300, R-1301, R-1304 (all surface as a parity violation)

**Test 5: ADER orders 3 and 4 (np=1 first, then np=2)**
- Target: `iterator.SetSubSteps` + `AdvanceWithSubStepStates` + driver
- File: extend existing `test_tpv104_substep_iterator_parity.cpp`
- Validates: Currently only `ader_order = 2` is tested. Add parameterised
  cases for {2, 3, 4}, asserting the iterator path remains stable and
  `WriteBackState` writes consistent V1/V2 across orders.
- Priority: MEDIUM
- Detects: latent ADER-order bugs in the predictor's CK recursion.

**Test 6 (R-1306, optional): long-run drift of `dt_scale` round-trip**
- Target: `AdvanceADERWithSubStep` configured-deltaT mutation
- Priority: LOW
- Skeleton: run 10⁴ macro-steps; verify `Σ deltaT[o] / dt_step` drift is
  bounded by `N * eps`.

---

## Priority Summary

| ID | Severity | Status | Description |
|------|------|---|---|
| R-1300 | CRITICAL | OPEN | Shared-fault ADER branch never consults substep side-channel |
| R-1301 | CRITICAL | OPEN | Iterator buffer sized to local-only QPs; needs total |
| R-1302 | CRITICAL | OPEN | EvaluateBulkAtFaultQPsCanonical skips shared faces |
| R-1303 | HIGH     | OPEN | ComputeADERSubStepStates predictor needs ghost-cell coverage |
| R-1304 | HIGH     | OPEN | dof_idx into substep buffer assumes interior-only ordering |
| R-1305 | MODERATE | OPEN | Per-sub-step DOFData symmetry across shared-fault ranks |
| R-1306 | LOW      | OPEN | dt_scale re-SetSubSteps cumulative drift |

### Test Coverage Summary

| Suite | File | Tests | Detected bugs |
|---|---|---|---|
| MPI parity      | `test_tpv104_substep_iterator_mpi.cpp` (new)              | 1 | R-1300, R-1301, R-1302, R-1304 |
| Predictor MPI   | `test_tpv104_compute_ader_substep_states_mpi.cpp` (new)   | 1 | R-1303 |
| Symmetry        | `test_tpv104_substep_iterator_shared_face_symmetry.cpp` (new) | 1 | R-1305 |
| One-shot MPI    | `test_tpv104_substep_one_shot_parity_mpi.cpp` (new)       | 1 | R-1300/1/4 |
| ADER orders     | extend `test_tpv104_substep_iterator_parity.cpp`          | 2 | latent CK |
| Drift (opt)     | `test_tpv104_substep_iterator_dt_scale.cpp` (new)         | 1 | R-1306 |
| **Total**       |                                                           | **7** | |

---

## R-1003 Implementation Plan

The plan in five steps. Each step is independent; merge in order.

### §1 — Predictor extension to ghost cells (fixes R-1303)

**File:** `dynamic/wave_operator.inl` (`ComputeADERSubStepStates` body, ~L1043)

Verify the existing `ComputeADERTimeIntegrated` path: does the CK recursion
already run on `pmesh->face_nbr_data` ghost cells, or does the macro-step path
exchange the time-integrated I via `pmesh->ExchangeFaceNbrData()`?

- If the existing macro-step path runs the recursion on ghost cells:
  `ComputeADERSubStepStates` is **already correct** at the bulk level — the
  bug is downstream in `EvaluateBulkAtFaultQPsCanonical` (§3 below). No
  change here.
- If the existing macro-step path exchanges the integrated I after the
  recursion: add a ghost exchange at the END of `ComputeADERSubStepStates`,
  one per `tau_nodes[o]`:

  ```cpp
  // In ComputeADERSubStepStates, after the per-node Taylor evaluation:
  for (int o = 0; o < tau_nodes.size(); o++) {
     // … fill Q_per_node[o] with locally-owned bulk values …
     // R-1303: exchange ghost layer so EvaluateBulkAtFaultQPsCanonical
     // can read neighbour's Q_per_node[o] for shared-fault face Elem2.
     ParGridFunction gf(pfes_);
     gf.SetFromTrueDofs(Q_per_node[o]);   // or appropriate vector mapping
     gf.ExchangeFaceNbrData();
     // gf.FaceNbrData() now has neighbour Q_per_node[o].
     // Stash a pointer; consumer reads from it.
     face_nbr_substep_data_[o] = std::move(gf.FaceNbrData());
  }
  ```

  Add a `mutable std::vector<Vector> face_nbr_substep_data_` member to
  `WaveOperator` and a `GetFaceNbrSubStepData(int o)` accessor.

  **Cost:** one MPI exchange per sub-step per macro-step (`O` exchanges
  total). For O=4 with TPV104 typical dt this is 4 extra collectives per
  step — acceptable. The non-substep path is unchanged.

### §2 — Substep gate on the shared-fault ADER branch (fixes R-1300, R-1304)

**File:** `dynamic/wave_operator.inl` lines 3980-3988

**Current:**
```cpp
real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
fault_flux_->EvaluateADER(fdata,
                          I_plus_local, I_minus_local,
                          dt,
                          I_imp_plus, I_imp_minus);
```

**Replace with** (mirror the L3086-3110 interior-branch gate):
```cpp
real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
// R-1003: substep dispatch (mirrors interior-branch gate at L3086-3110).
// dof_idx is already absolute (interior + shared concatenated) per
// shared_fault_dof_offset_ population in the ctor.
if (substep_I_imp_plus_flat_ != nullptr &&
    substep_I_imp_minus_flat_ != nullptr &&
    dof_idx >= 0 &&
    dof_idx < substep_n_local_fault_qps_)   // see R-1301 rename
{
   const real_t *src_p =
      substep_I_imp_plus_flat_  + dof_idx * NUM_STATE;
   const real_t *src_m =
      substep_I_imp_minus_flat_ + dof_idx * NUM_STATE;
   for (int c = 0; c < NUM_STATE; c++)
   {
      I_imp_plus[c]  = src_p[c];
      I_imp_minus[c] = src_m[c];
   }
}
else
{
   fault_flux_->EvaluateADER(fdata,
                             I_plus_local, I_minus_local,
                             dt,
                             I_imp_plus, I_imp_minus);
}
```

The `dof_idx` lookup at L3926-3928 already produces an absolute index in the
total-fault range because `shared_fault_dof_offset_` is populated with the
interior offset baked in (`wave_operator.hpp:439`).

### §3 — Extend `EvaluateBulkAtFaultQPsCanonical` to shared faces (fixes R-1302)

**File:** `dynamic/wave_operator.inl` lines 1422-1545

After the existing interior loop (L1543, the closing `}` of the
`fi < fault_interior_faces_.Size()` loop), add:

```cpp
// R-1003 §3: shared-fault loop — populates the
// [GetNumLocalFaultQPs(), GetNumTotalFaultQPs()) slice of Q_*_flat.
for (int sf_idx = 0; sf_idx < fault_shared_faces_.Size(); sf_idx++)
{
   int sf = fault_shared_faces_[sf_idx];
   FaceElementTransformations *ftr =
      const_cast<ParMesh&>(*pmesh_).GetSharedFaceTransformations(sf);
   MFEM_VERIFY(ftr != nullptr,
               "EvaluateBulkAtFaultQPsCanonical: shared face "
               << sf << " has no transformations.");

   const int e1 = ftr->Elem1No;       // local
   const int e2 = ftr->Elem2No;       // ghost / face-nbr
   const FiniteElement *fe1 = fes_->GetFE(e1);
   const FiniteElement *fe2 = pfes_->GetFaceNbrFE(e2 - mesh_.GetNE());
   const int dof_offset1 = e1 * ndof_per_el_;
   // Face-nbr data is in pfes_->FaceNbrData(); offset within FaceNbrData is
   // pfes_->GetFaceNbrDofOffsets()[e2 - mesh_.GetNE()] (or equivalent).
   const int e2_local = e2 - mesh_.GetNE();
   const int dof_offset2_nbr = pfes_->GetFaceNbrDofOffset(e2_local);

   const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                            2 * order_);
   const int nqp = ir.GetNPoints();
   MFEM_VERIFY(nqp == nbf_per_face_, "shared face nqp mismatch");

   auto it_off = shared_fault_dof_offset_.find(sf);
   MFEM_VERIFY(it_off != shared_fault_dof_offset_.end(),
               "shared_fault_dof_offset_ missing entry for sf=" << sf);
   const int base_dof_idx = it_off->second;   // absolute (>= NumLocal)

   const bool elem1_on_plus = shared_fault_elem1_on_plus_[sf_idx];

   const int fb_idx = fault_interior_faces_.Size() + sf_idx;
   const FaultBasisData &bd = fault_basis_->GetBasis(fb_idx);

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      ftr->SetAllIntPoints(&ip);

      IntegrationPoint ip1, ip2;
      ftr->Loc1.Transform(ip, ip1);
      ftr->Loc2.Transform(ip, ip2);
      Vector shape1(fe1->GetDof()), shape2(fe2->GetDof());
      fe1->CalcShape(ip1, shape1);
      fe2->CalcShape(ip2, shape2);

      // Bulk Q at QP from each side.  Self side: locally-owned.
      // Neighbour side: pfes_->FaceNbrData() (ghost layer).
      const Vector &fnbr = pfes_->GetFaceNbrData();   // NUM_STATE * n_face_nbr_dofs
      real_t Q_self[NUM_STATE], Q_nbr[NUM_STATE];
      const real_t *Q_data = Q_bulk.GetData();
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t s_self = 0.0, s_nbr = 0.0;
         for (int i = 0; i < shape1.Size(); i++)
            s_self += shape1(i) * Q_data[c*ndof_total_ + dof_offset1 + i];
         for (int i = 0; i < shape2.Size(); i++)
            s_nbr  += shape2(i) * fnbr[c*pfes_->GetFaceNbrVSize()
                                       + dof_offset2_nbr + i];
         Q_self[c] = s_self;
         Q_nbr[c]  = s_nbr;
      }

      // Same canonical-frame rotation as L1525-1545 + L3964-3973
      // (interior + ADER shared branches).
      const FaultBasisQPData &qpd = bd.qp_data[q];
      real_t can_n[3], can_t1[3], can_t2[3];
      const bool should_negate_frame = !elem1_on_plus;
      for (int d = 0; d < 3; d++)
      {
         can_n[d]  = should_negate_frame ? -qpd.normal[d]   : qpd.normal[d];
         can_t1[d] = should_negate_frame ? -qpd.tangent1[d] : qpd.tangent1[d];
         can_t2[d] = should_negate_frame ? -qpd.tangent2[d] : qpd.tangent2[d];
      }
      DenseMatrix Tinv_can(NUM_STATE);
      GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);

      real_t Q_self_can[NUM_STATE], Q_nbr_can[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++)
      {
         Q_self_can[c] = 0.0; Q_nbr_can[c] = 0.0;
         for (int k = 0; k < NUM_STATE; k++)
         {
            Q_self_can[c] += Tinv_can(c, k) * Q_self[k];
            Q_nbr_can[c]  += Tinv_can(c, k) * Q_nbr[k];
         }
      }

      // Pack into the output flat arrays at absolute index.
      const int dof_idx_global = base_dof_idx + q;
      real_t *out_p = Q_plus_flat.data()  + dof_idx_global * NUM_STATE;
      real_t *out_m = Q_minus_flat.data() + dof_idx_global * NUM_STATE;
      const real_t *plus_src  = elem1_on_plus ? Q_self_can : Q_nbr_can;
      const real_t *minus_src = elem1_on_plus ? Q_nbr_can  : Q_self_can;
      for (int c = 0; c < NUM_STATE; c++)
      {
         out_p[c] = plus_src[c];
         out_m[c] = minus_src[c];
      }
   }
}
```

Also resize the output at L1433-1437 to use `GetNumTotalFaultQPs()`:

```cpp
const int n_total_qps = GetNumTotalFaultQPs();   // was: GetNumLocalFaultQPs()
const size_t expect_words = NUM_STATE * n_total_qps;
Q_plus_flat.assign(expect_words, 0.0);
Q_minus_flat.assign(expect_words, 0.0);
```

### §4 — Driver buffer sizing + parameter rename (fixes R-1301)

**File:** `drivers/tpv104_driver.cpp` line 362, 372-393

**Current:**
```cpp
const int n_local_fault_qps = wave.GetNumLocalFaultQPs();
const size_t n_words = NUM_STATE * n_local_fault_qps;
std::vector<real_t> I_imp_plus_flat(n_words, 0.0);
// ...
wave.SetSubStepFaultImposedStates(I_imp_plus_flat.data(),
                                  I_imp_minus_flat.data(),
                                  n_local_fault_qps);
```

**Replace with:**
```cpp
const int n_total_fault_qps = wave.GetNumTotalFaultQPs();   // R-1003
const size_t n_words = NUM_STATE * n_total_fault_qps;
std::vector<real_t> I_imp_plus_flat(n_words, 0.0);
std::vector<real_t> I_imp_minus_flat(n_words, 0.0);

if (n_total_fault_qps > 0)
{
   iterator.AdvanceWithSubStepStates(dof_data, fault_coords, V_w,
                                     Q_pointwise_plus, Q_pointwise_minus,
                                     dt_step, t_step_start,
                                     I_imp_plus_flat.data(),
                                     I_imp_minus_flat.data(),
                                     method);
}

wave.SetSubStepFaultImposedStates(
   n_total_fault_qps > 0 ? I_imp_plus_flat.data()  : nullptr,
   n_total_fault_qps > 0 ? I_imp_minus_flat.data() : nullptr,
   n_total_fault_qps);
```

In `wave_operator.hpp` rename:
- `int n_local_fault_qps` parameter → `int n_total_fault_qps` at L304
- `mutable int substep_n_local_fault_qps_ = 0;` → `substep_n_total_fault_qps_`
  at L552
- update all references in `wave_operator.inl` (L3089, L3937 once added).

### §5 — Driver guard removal + runtime probe

**File:** `drivers/tpv104_driver.cpp` lines 1349-1365

**Replace the abort with a runtime sanity check** that confirms the
shared-fault wiring has actually landed (defence in depth — catches a
half-merged R-1003):

```cpp
// R-1003 (round-12+, post-implementation): the shared-fault path now
// honours the substep side-channel.  Verify the wiring is in place by
// running a smoke test once at startup: with substep buffers set to a
// sentinel pattern, calling AdvanceADER on a Q_init must NOT produce the
// EvaluateADER-default I_imp signature on shared-fault QPs.  If it does,
// the build is half-merged and we abort loudly rather than silently
// regressing.
if (use_substep_iterator && nprocs > 1)
{
   bool shared_path_wired = ProbeSharedFaultSubStepWiring(wave, fault_flux);
   MFEM_VERIFY(shared_path_wired,
               "--fault-iterator substep is enabled at np=" << nprocs
               << " but the shared-fault ADER branch in "
               << "ComputeADERSharedFaceFluxRHS is NOT consuming the "
               << "substep side-channel (R-1003 half-merged).  Rebuild "
               << "with the wave_operator.inl shared-branch patch.");
}
```

`ProbeSharedFaultSubStepWiring` runs a single-step diagnostic at startup
(before the main time loop) that:
1. zeros Q,
2. installs a sentinel `I_imp_*_flat` pattern (e.g., I_imp_plus = 1.0, I_imp_minus = 2.0),
3. calls `wave.AdvanceADER(Q, dt_probe, ader_order, Q_new)`,
4. inspects the resulting `Q_new` at a known shared-fault QP via
   `wave.GetSharedFaultDof(...)` — if the substep gate fires, the bulk
   contribution at that QP must reflect the sentinel; if EvaluateADER fired
   it would reflect the friction-solver output (which is non-trivial).

Alternative: gate behind a build-time symbol `SEAS_TPV104_SUBSTEP_MPI` set by
the build system once §1-§4 are merged. Cheaper but less robust.

### §6 — MPI consistency assertions for DOFData (fixes R-1305)

**File:** new `dynamic/tpv104_substep_iterator_mpi_check.hpp` and `.cpp`

Add a debug-mode function:
```cpp
#ifdef SEAS_DIAG_TPV104_SUBSTEP_MPI_CHECK
void Tpv104SubStepIterator::AssertSharedFaultDOFDataSymmetry(
   const std::vector<DOFData> &dof_data,
   int num_local_fault_qps,
   int num_total_fault_qps,
   MPI_Comm comm) const;
#endif
```

Iterate `[num_local_fault_qps, num_total_fault_qps)`, for each scalar
(d.psi, d.slip1, d.slip2, d.V_abs, d.tau1_corr, d.tau2_corr,
d.sigma_n_corr) call MPI_Allreduce(MAX) and MPI_Allreduce(MIN) and assert
equal bit-for-bit. Hook from the driver after each
`AdvanceADERWithSubStep` call when the diag flag is set.

---

## Code Review Report (Rev 1)

**Date:** 2026-04-25
**Scope:** TPV104 sub-step iterator MPI extension (R-1003 readiness)
**Reviewer:** chunhui-code-reviewer agent

### Summary
The np=1 sub-step iterator is correct and well-tested. The driver-level
R-1003 abort is correctly placed and accurately describes the limitation.
The shared-fault MPI path requires four file edits (interior-branch mirror,
buffer sizing fix, EvaluateBulk extension, predictor ghost extension or
verification) plus one rename (`n_local` → `n_total`). All four are local
patches — no new MPI collectives required if the existing ADER predictor
already covers ghost cells. No design-level rework needed.

### Code Review Findings

| # | File | Lines | Severity | Description |
|---|---|---|---|---|
| 1 | `wave_operator.inl` | 3980-3988 | CRITICAL | Shared-fault ADER branch lacks substep gate (R-1300) |
| 2 | `tpv104_driver.cpp` | 362,372-393 | CRITICAL | Iterator buffer sized to local-only QPs (R-1301) |
| 3 | `wave_operator.inl` | 1422-1545 | CRITICAL | EvaluateBulkAtFaultQPsCanonical interior-only (R-1302) |
| 4 | `wave_operator.inl` | ~1043 | HIGH     | ComputeADERSubStepStates ghost coverage (R-1303) |
| 5 | `wave_operator.inl` | 3927-3940 | HIGH     | dof_idx routing into substep buffer on shared branch (R-1304) |
| 6 | `tpv104_substep_iterator.cpp` | 410-419, 651-660 | MODERATE | DOFData symmetry across shared-fault ranks (R-1305) |
| 7 | `tpv104_driver.cpp` | 314-333 | LOW      | dt_scale round-trip cumulative drift (R-1306) |

### Documentation Status
The driver guard at L1349-1365 is correctly worded. Once R-1003 lands,
update the comment block to a "wiring health probe" description (§5 above).
Header docstrings at `wave_operator.hpp:292-308, 310-333` need wording
updates: "interior fault QP" → "fault QP (interior + shared)" and the
parameter rename `n_local_fault_qps` → `n_total_fault_qps`.

### Test Status
- **Implemented:** 0 of 7 proposed MPI tests (all need to be authored).
- **Existing np=1 coverage:** strong — `test_tpv104_substep_iterator.cpp`,
  `*_parity.cpp`, `*_dispatch_parity.cpp`, `*_one_shot_parity.cpp`,
  `*_predictor.cpp` all run cleanly at np=1.
- **Priority for R-1003 merge:** Tests 1, 2, 4 (HIGH) are merge-blocking;
  Tests 3, 5 (MEDIUM) within one week; Test 6 (LOW) optional.

### Recommended Next Steps
1. Verify the existing macro-step path's ghost-cell coverage of the CK
   recursion (drives §1 vs §1-no-op). Read
   `wave_operator.inl::ComputeADERTimeIntegrated`.
2. Land §2 (shared-fault gate) and §4 (driver buffer rename) — these are
   trivial mechanical patches that can be verified by Test 4 at O=1
   (T_TPV104_SSI_3 lifted to MPI).
3. Land §3 (EvaluateBulk extension) — covered by Test 1 + Test 2.
4. Land §5 (guard replacement) only after Tests 1, 2, 4 pass at np=2 and
   np=4.
5. Defer §6 (MPI invariant checker) to a separate diag-only commit.
6. Add ADER-O ∈ {3, 4} to existing np=1 parity tests (Test 5) so R-1003
   merge does not regress higher-order behaviour.
