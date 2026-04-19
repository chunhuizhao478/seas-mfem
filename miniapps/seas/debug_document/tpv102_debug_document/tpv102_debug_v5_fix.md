# Fix Report v5: address review findings in `tpv102_debug_v5_check.md`

## Summary
- Findings addressed: **4 of 4** (R-501 CRITICAL; R-502, R-503 MODERATE; R-504 LOW)
- Files modified:
  - `miniapps/seas/dynamic/wave_operator.inl`
    (R-501 restructure of `ComputeSharedFaceFluxRHS`; R-502 scale-relative
    tol; R-503 mesh-scale floor in `same_centroid`; R-504 abort wording)
  - `miniapps/seas/tests/parallel/test_r101_shared_fault.cpp`
    (R-501a new test: multi-stage Mult on nonzero Q)
- Tests added: **1** (`TestR501_MultiStageRK4NonzeroQ`)
- Test suite: **PASS**
  - `seas_test_godunov_flux`                         — 29/29
  - `seas_test_wave_operator`                        — 17/17
  - `seas_test_wave_bc`                              — 10/10
  - `seas_test_fault_face_flux`                      — 19/19
  - `seas_test_tpv102_setup`                         — 24/24
  - `seas_test_tpv102_local`                         — 16/16
  - `seas_test_parallel_wave_operator` (2 ranks)     — 5/5
  - `seas_test_parallel_wave_operator` (4 ranks)     — 5/5
  - `seas_test_r101_shared_fault`      (2 ranks)     — **10/10**
    - R-302a (inline 2-tet, `Mult(Q=0)`) — `max_rel_diff=0`
    - **R-501a (inline 2-tet, 4× `Mult` on nonzero Q) — `max_rel_diff=0`**
    - R-302a & R-501a both would abort under the pre-R-501 code; both
      now pass end-to-end.

## Changes Made

### [R-501] [CRITICAL] — Owner broadcast of post-Evaluate state
**Where:** `dynamic/wave_operator.inl::ComputeSharedFaceFluxRHS`

The v3 R-001 (+,-) swap corrected the Evaluate argument order but did
not — and *cannot* — correct the rank-local fault-frame mismatch.  MFEM's
`GetSharedFaceTransformations` returns opposite face normals on the two
ranks sharing the face, so `BuildFrame` produces frames
`(n, t1, t2)_A = (+n, +t1, +t2)` vs `(−n, −t1, +t2)_B`; the rotation
matrix `Tinv` that maps global Q into fault-local coordinates is
different on the two ranks.  Because `FaultFaceFlux::Evaluate` is
nonlinear (Brent friction solve), different inputs produce different
outputs even when the swap makes the (+,-) labels match — which is
exactly the `tau1_corr differs by 1.9e-7` symptom on Frontera job 7664983.

**Fix:** restructure `ComputeSharedFaceFluxRHS` into three phases:

1. **Phase 1 (Owner)**: the owner (lower rank ID) computes its own
   `Q_plus_local`, `Q_minus_local` in its own fault-local frame, calls
   `Evaluate` to get updated `fdata` and `Q_imp_plus / Q_imp_minus`,
   then rotates the imposed states back to the **global** frame with
   its own `T`.  Owner packs 27 doubles per shared-fault QP into a
   send buffer: 9 mutable DOFData scalar fields + 9 `Q_imp_plus_g` +
   9 `Q_imp_minus_g`.  Non-owner does NOT call `Evaluate`.

2. **Phase 2 (MPI)**: grouped by peer rank, owners issue `MPI_Isend`
   and non-owners `MPI_Irecv`; a single `MPI_Waitall` completes the
   batch.  All fault-shared-face QPs between each (rank, peer) pair
   are exchanged in one message (no per-QP overhead).

3. **Phase 3 (Flux)**: both ranks iterate all shared faces; for fault
   QPs, use the authoritative `Q_imp_plus_g / Q_imp_minus_g` from the
   `auth_state` table (owner-computed or non-owner-received).  Each
   rank assigns self/nbr based on its own Elem1:
   - Owner: `Q_self_imp = Q_imp_plus_g`, `Q_nbr_imp = Q_imp_minus_g`.
   - Non-owner: `Q_self_imp = Q_imp_minus_g`, `Q_nbr_imp = Q_imp_plus_g`.
   Both call `flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h)` with
   their own `nor` — no explicit sign flip is needed.  Conservation
   is handled by the Godunov flux identity
   `F(Q_self, Q_nbr, nor) = −F(Q_nbr, Q_self, −nor)`:
   flipping both the (self, nbr) assignment and `nor` together leaves
   F unchanged, so accumulation `rhs[Elem1] -= F_h` on one rank and
   the equivalent +F_h on the other rank gives the same bulk-flux
   balance as the pre-partition interior-face case.

**Why not the reviewer's `F_h = -F_h` flip?**  The review's pseudo-code
suggested *both* ranks call `flux_.Interior(nor, Q_imp_plus_g,
Q_imp_minus_g, F_h)` with the same Q arguments but different `nor`,
then the non-owner flips `F_h`.  That is **incorrect** for the upwind
Godunov flux: `F(L, R, −n)` is NOT `−F(L, R, n)` in general (the upwind
direction depends on the sign of `a · n`, which flips, but `L` and `R`
are unchanged, giving a flux value that is related to neither `F(L, R, n)`
nor `−F(L, R, n)` in a simple way).  The conservation identity requires
swapping `L ↔ R` together with flipping `n`.  My fix uses the correct
form.  (This deviation is the only place I diverge from the review's
suggested pseudo-code; the R-501 design intent is unchanged.)

**Non-owner DOFData write-back**: as part of Phase 2, non-owner copies
the received mutable fields into its local `fault_dof_data_` entry.
This keeps the driver's RK4 averaging, the R-101 verifier, and any
downstream observer (station writer, ParaView output) consistent across
ranks — every rank sees the owner's authoritative values.

### [R-502] [MODERATE] — Scale-relative field tolerance in verifier
**Where:** `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency`

Replaced the absolute `max_diff > tol` comparison with a scale-relative
one.  Field magnitudes span `1e-12` (V near nucleation) to `1e+8` Pa
(normal stress); an absolute `tol = 1e-10` is simultaneously too tight
for large-magnitude fields (sub-ULP reassociation trips it) and too
loose for tiny ones.  The fix tracks `max_rel_diff = max_k(|a_k −
b_k| / max(|a_k|, |b_k|, 1))`; `tol` is now interpreted as a
**relative** tolerance.  Abort / OK-print both report `max_rel_diff`,
`max_abs_diff`, and the field scale at the max-diff entry so the
diagnostic is still readable in absolute units.

The `max(…, 1.0)` floor (a) never divides by zero, (b) is O(1) for
dimensionless fields like `psi`, and (c) matches the field magnitude
for large dimensioned fields.

### [R-503] [MODERATE] — Mesh-scale floor in `same_centroid`
**Where:** `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency`

R-404's `scale × DBL_EPSILON` collapses to `~2e-28` when a coordinate
is near zero (the TPV102 fault plane has y = 0).  Two ranks that
happen to produce 1-ULP-of-nonzero-y drift would then be classified
as different centroids and trigger a spurious R-305 unpair-abort.

Fix: hybrid tolerance `max(scale × DBL_EPSILON, abs_floor)` with
`abs_floor = 1e-9 m` — smaller than any fault-mesh element (TPV102
mesh h ~ 200–1000 m) and larger than any realistic FP noise at
mesh-scale physical coordinates.

### [R-504] [LOW] — Post-R-501 abort-message wording
**Where:** `dynamic/wave_operator.inl::VerifySharedFaultDOFDataConsistency`

The prior message blamed "R-001's (+,-) canonicalisation is not
sufficient".  With R-501 in place, the R-001 path no longer runs —
any future verifier trip is from a different bug (missed field in
the broadcast, packing mismatch, MPI exchange error).  New message
describes the *symptom* generically and points at the R-501
broadcast + Q_imp exchange for investigation.

### [R-501a] New test — Multi-stage Mult on nonzero Q
**Where:** `tests/parallel/test_r101_shared_fault.cpp::TestR501_MultiStageRK4NonzeroQ`

Uses the same `BuildTwoTetSharedFaultMeshInline` mesh as R-302a but
drives 4 `wave.Mult(Q, k)` calls with a small sinusoidal nonzero Q
(`Q[i] = 1e3 × sin(0.37i + 0.9 × rank)`).  This replicates the
Frontera 4-rank failure condition: nonzero Q → frame mismatch
propagates into Evaluate → DOFData drifts across ranks.

Pre-R-501 behavior: aborts with `tau1_corr differs by ~1e-7` (same
mode as the Frontera stderr).
Post-R-501 behavior: passes with `max_rel_diff = 0` — owner broadcast
produces bit-identical DOFData on both ranks.

## Unresolved Findings
None.

## Verification
- [x] R-501: fixed — owner-broadcast 3-phase restructure in
  `ComputeSharedFaceFluxRHS`.  Non-owner receives 27-double authoritative
  record per QP and assembles flux in its own frame with corrected
  self/nbr assignment.  Test R-501a exercises nonzero-Q multi-Mult;
  `max_rel_diff = 0`.
- [x] R-502: fixed — tol is now relative; `max_rel_diff = |a−b| /
  max(|a|, |b|, 1.0)`.  Diagnostic reports both absolute and relative.
- [x] R-503: fixed — `same_centroid` uses `max(scale×DBL_EPSILON,
  1e-9 m)`.
- [x] R-504: fixed — abort message no longer blames R-001 and points
  at R-501 broadcast/exchange for debugging.

## Ready for Re-Review: YES

## Notes for the next Frontera run
- The previously-blocked sbatch runs should now be safe to submit:
  - `tpv102_1000m_p1_1.5s_4rank_dev.sbatch`
  - `tpv102_1000m_p1_0.5s_8rank_pvsmoke_dev.sbatch`
  - `tpv102_200m_p1_1.5s_50rank_dev.sbatch`
  - `tpv102_200m_p1_1.5s_400rank_dev.sbatch`
- Expected new `[R-101 check]` output line format:
  `shared-fault DOFData consistency OK: N pairs matched, 0 unpaired entries,
   max_rel_diff=<r> (max_abs_diff=<a>, scale=<s>, rel_tol=1e-10)`.
- If a future trip occurs despite R-501, the abort message explicitly
  points the reader at the R-501 packing/exchange in
  `ComputeSharedFaceFluxRHS`.
